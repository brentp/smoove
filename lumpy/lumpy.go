package lumpy

import (
	"bufio"
	"bytes"
	"fmt"
	"html/template"
	"io"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
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
	Name          string   `arg:"-n,required,help:project name used in output files."`
	Fasta         string   `arg:"-f,required,help:fasta file."`
	Exclude       string   `arg:"-e,help:BED of exclude regions."`
	ExcludeChroms string   `arg:"-C,help:ignore SVs with either end in this comma-delimited list of chroms. If this starts with ~ it is treated as a regular expression to exclude."`
	Processes     int      `arg:"-p,help:number of processors to parallelize."`
	OutDir        string   `arg:"-o,help:output directory."`
	Svtyper       bool     `arg:"-s,help:run svtyper directly on output"`
	Bams          []string `arg:"positional,required,help:path to bam(s) to call."`
}

func (c cliargs) Description() string {
	return "this runs lumpy ands sends output to STDOUT"
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

func Lumpy(project, reference string, outdir string, bam_paths []string, out io.Writer, pool *shpool.Pool, exclude_bed string, filter_chroms []string) *exec.Cmd {
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

	remove_sketchy_all(filters, 1000, reference, exclude_bed, filter_chroms)
	shared.Slogger.Print("starting lumpy")
	p := run_lumpy(filters, reference, outdir, false, project)
	return p
}

func run_lumpy(bams []filter, fa string, outdir string, has_cnvnator bool, name string) *exec.Cmd {
	if _, err := exec.LookPath("lumpy"); err != nil {
		log.Fatal(shared.Prefix + " lumpy not found on path")
	}

	lumpy_tmpl := "set -euo pipefail; lumpy -msw 3 -mw 4 -t $(mktemp) -tt 0 -P "
	pe_tmpl := "-pe id:{{.Sample}},bam_file:{{.DiscPath}},histo_file:{{.HistPath}},mean:{{.Mean}},stdev:{{.Std}},read_length:{{.ReadLength}},min_non_overlap:{{.ReadLength}},discordant_z:4,back_distance:30,weight:1,min_mapping_threshold:" + strconv.Itoa(int(MinMapQuality)) + " "
	sr_tmpl := "-sr id:{{.Sample}},bam_file:{{.SplitPath}},back_distance:10,weight:1,min_mapping_threshold:" + strconv.Itoa(int(MinMapQuality)) + " "

	del_tmpl := "-bedpe bedpe_file:%s/%s.del.bedpe,id:%s,weight:2 "
	dup_tmpl := "-bedpe bedpe_file:%s/%s.dup.bedpe,id:%s,weight:2 "

	var buf bytes.Buffer

	for _, sample := range bams {
		S := cs_from_filter(sample, outdir)
		t, err := template.New(sample.sample).Parse(pe_tmpl + sr_tmpl)
		check(err)
		check(t.Execute(&buf, S))
		if has_cnvnator {
			buf.Write([]byte(fmt.Sprintf(del_tmpl, outdir, sample.sample, sample.sample)))
			buf.Write([]byte(fmt.Sprintf(dup_tmpl, outdir, sample.sample, sample.sample)))
		}
	}

	cmdStr := lumpy_tmpl + buf.String()
	f, err := os.Create(outdir + "/" + name + "-lumpy-cmd.sh")
	shared.Slogger.Write([]byte("wrote lumpy command to " + f.Name()))
	check(err)
	f.WriteString(cmdStr + "\n")
	f.Close()

	return exec.Command("bash", "-c", cmdStr)
}

type cs struct {
	Sample     string
	DiscPath   string
	SplitPath  string
	HistPath   string
	Mean       string
	Std        string
	ReadLength string
}

func cs_from_filter(f filter, outdir string) cs {
	return cs{Sample: f.sample, DiscPath: f.disc, SplitPath: f.split, HistPath: f.histpath(outdir),
		Mean: fmt.Sprintf("%.3f", f.stats.TemplateMean), Std: fmt.Sprintf("%.3f", f.stats.TemplateSD),
		ReadLength: strconv.Itoa(f.stats.MaxReadLength),
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
	if _, err := exec.LookPath("lumpy"); err != nil {
		shared.Slogger.Fatal("lumpy executable not found in PATH")
	}
	cli := cliargs{Processes: 3, ExcludeChroms: "hs37d5,~:,~^GL,~decoy"}
	arg.MustParse(&cli)
	runtime.GOMAXPROCS(cli.Processes)
	wtr := bufio.NewWriter(os.Stdout)
	filter_chroms := strings.Split(strings.TrimSpace(cli.ExcludeChroms), ",")
	if cli.OutDir == "" {
		cli.OutDir = "./"
	}
	p := Lumpy(cli.Name, cli.Fasta, cli.OutDir, cli.Bams, wtr, nil, cli.Exclude, filter_chroms)
	p.Stderr = shared.Slogger
	p.Stdout = os.Stdout
	check(p.Run())
}
