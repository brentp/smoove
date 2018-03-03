package lumpy

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	"sync"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/go-athenaeum/shpool"
	"github.com/brentp/goleft/covstats"
	"github.com/brentp/goleft/indexcov"
	"github.com/brentp/lumpy-smoother/shared"
	"github.com/brentp/xopen"
)

type cliargs struct {
	Fasta     string   `arg:"-f,required,help:fasta file."`
	Processes int      `arg:"-p,number of processors to parallelize."`
	OutDir    string   `arg:"-o,help:output directory."`
	Bams      []string `arg:"positional,required,help:path to bam to call."`
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

type filter struct {
	bam     string
	split   string
	disc    string
	command string
	sample  string
	stats   covstats.Stats
}

func (f filter) histpath(outdir string) string {
	return fmt.Sprintf("%s/%s.histo", outdir, f.sample)
}

func (fi filter) write_hist(outdir string) {
	f, err := os.Create(fi.histpath(outdir))
	check(err)
	for i, v := range fi.stats.H {
		fmt.Fprintf(f, "%d\t%.12f\n", i, v)
	}
	f.Close()
}

func lumpy_filter_cmd(bam string, outdir string, reference string) filter {
	// check if .split.bam and .disc.bam exist. if they do, then use.
	sm, err := indexcov.GetShortName(bam, strings.HasSuffix(bam, ".bam"))
	check(err)
	prefix := fmt.Sprintf("%s/%s", outdir, sm)

	if xopen.Exists(prefix+".split.bam") && xopen.Exists(prefix+".disc.bam") {
		return filter{bam: bam, split: prefix + ".split.bam", disc: prefix + ".disc.bam", sample: sm}
	}

	// symlink to out dir.
	olddir := filepath.Dir(bam)
	if xopen.Exists(fmt.Sprintf("%s/%s.split.bam", olddir, sm)) && xopen.Exists(fmt.Sprintf("%s/%s.disc.bam", olddir, sm)) {
		check(os.Symlink(fmt.Sprintf("%s/%s.split.bam", olddir, sm), fmt.Sprintf("%s/%s.split.bam", outdir, sm)))
		check(os.Symlink(fmt.Sprintf("%s/%s.disc.bam", olddir, sm), fmt.Sprintf("%s/%s.disc.bam", outdir, sm)))
		return filter{bam: bam, split: prefix + ".split.bam", disc: prefix + ".disc.bam", sample: sm}
	}

	f := filter{bam: bam, split: prefix + ".split.bam", disc: prefix + ".disc.bam", sample: sm}
	// use .tmp.bam in case of error while running lumpy filter.
	f.command = fmt.Sprintf("lumpy_filter -f %s %s %s.tmp.bam %s.tmp.bam %d && mv %s.tmp.bam %s && mv %s.tmp.bam %s", reference, bam, f.split, f.disc, 2, f.split, f.split, f.disc, f.disc)
	return f
}

func Lumpy(reference string, outdir string, bam_paths []string, out io.Writer, pool *shpool.Pool) {
	if pool == nil {
		pool = shpool.New(runtime.GOMAXPROCS(0), nil, &shpool.Options{LogPrefix: shared.Prefix})
	}
	filters := make([]filter, len(bam_paths))
	for i, p := range bam_paths {
		filter := lumpy_filter_cmd(p, outdir, reference)
		if filter.command != "" {
			pool.Add(shpool.Process{Command: filter.command, CPUs: 1, Prefix: "lumpy-filter"})
		}
		filters[i] = filter
	}
	bam_stats(filters, reference, outdir)
	if err := pool.Wait(); err != nil {
		panic(err)
	}

}

func bam_stats(bams []filter, fasta string, outdir string) {
	shared.Slogger.Printf("calculating bam stats for %d bams\n", len(bams))
	var wg sync.WaitGroup

	wg.Add(2)

	f := func(mod int) {
		for i, f := range bams {
			if i%2 == mod {
				continue
			}
			br, err := shared.NewReader(f.bam, 2, fasta)
			check(err)
			bams[i].stats = covstats.BamStats(br, 200000)
			bams[i].write_hist(outdir)
		}
		wg.Done()
	}

	go f(0)
	go f(1)
	wg.Wait()
	shared.Slogger.Print("done calculating bam stats\n")
}

func Main() {
	cli := cliargs{Processes: 3}
	arg.MustParse(&cli)
	runtime.GOMAXPROCS(cli.Processes)
	wtr := bufio.NewWriter(os.Stdout)
	Lumpy(cli.Fasta, cli.OutDir, cli.Bams, wtr, nil)
}
