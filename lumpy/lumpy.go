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
	"sort"
	"strconv"
	"strings"
	"sync"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/biogo/io/seqio/fai"
	"github.com/biogo/hts/bgzf"
	"github.com/brentp/faidx"
	"github.com/brentp/go-athenaeum/shpool"
	"github.com/brentp/goleft/covstats"
	"github.com/brentp/goleft/indexcov"
	"github.com/brentp/smoove/shared"
	"github.com/brentp/smoove/svtyper"
	"github.com/brentp/xopen"
	"github.com/pkg/errors"
)

type cliargs struct {
	Name           string   `arg:"-n,required,help:project name used in output files."`
	Fasta          string   `arg:"-f,required,help:fasta file."`
	Exclude        string   `arg:"-e,help:BED of exclude regions."`
	ExcludeChroms  string   `arg:"-C,help:ignore SVs with either end in this comma-delimited list of chroms. If this starts with ~ it is treated as a regular expression to exclude."`
	Processes      int      `arg:"-p,help:number of processors to parallelize."`
	OutDir         string   `arg:"-o,help:output directory."`
	NoExtraFilters bool     `arg:"-F,help:use lumpy_filter only without extra smoove filters."`
	Genotype       bool     `arg:"-s,help:stream output to svtyper for genotyping"`
	RemovePr       bool     `arg:"-x,help:remove PRPOS and PREND tags from INFO (only used with --gentoype)."`
	Bams           []string `arg:"positional,required,help:path to bam(s) to call."`
}

func (c cliargs) Description() string {
	return "this runs lumpy ands sends output to {outdir}/{name}-smoove.vcf.gz if --genotype is requested, the output goes to {outdir}/{name}-smoove.genotyped.vcf.gz"
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
	sm, err := indexcov.GetShortName(bam, strings.HasSuffix(bam, ".cram"))
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

func Lumpy(project, reference string, outdir string, bam_paths []string, pool *shpool.Pool, exclude_bed string, filter_chroms []string, extraFilters bool) *exec.Cmd {
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

	remove_sketchy_all(filters, 800, reference, exclude_bed, filter_chroms, extraFilters)
	shared.Slogger.Print("starting lumpy")
	p := run_lumpy(filters, reference, outdir, false, project)
	return p
}

func run_lumpy(bams []filter, fa string, outdir string, has_cnvnator bool, name string) *exec.Cmd {
	if _, err := exec.LookPath("lumpy"); err != nil {
		shared.Slogger.Fatal("lumpy not found on path")
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
		Mean: fmt.Sprintf("%.2f", f.stats.TemplateMean), Std: fmt.Sprintf("%.2f", f.stats.TemplateSD),
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

func writeContigs(b *bufio.Writer, fasta string) {
	fa, err := faidx.New(fasta)
	check(errors.Wrapf(err, "error opening fasta file: %s", fasta))
	ctgs := make([]fai.Record, 0, len(fa.Index))
	for _, idx := range fa.Index {
		ctgs = append(ctgs, idx)
	}
	sort.Slice(ctgs, func(i, j int) bool { return ctgs[i].Start < ctgs[j].Start })
	for _, ctg := range ctgs {
		fmt.Fprintf(b, "##contig=<ID=%s,length=%d>\n", ctg.Name, ctg.Length)
	}
}

// filter BND variants from in that have < bndSupport
// also sneak in contig header.
func bndFilter(in io.Reader, bndSupport int, fasta string) io.Reader {
	r, w := io.Pipe()
	b := bufio.NewReader(in)
	wb := bufio.NewWriter(w)
	contigsWritten := false
	go func() {
		defer w.Close()
		defer wb.Flush()
		for {
			line, err := b.ReadString('\n')
			if len(line) > 0 {
				if line[0] == '#' {
					wb.WriteString(line)
					if !contigsWritten {
						writeContigs(wb, fasta)
						contigsWritten = true
					}

					continue
				}
				// only filter out BNDs with < X pieces of evidence when we're streaming directly from lumpy.
				if strings.Contains(line, "SVTYPE=BND") {
					toks := strings.SplitN(line, "\t", 9)
					// BND elements have different filtering set by command-line.
					idx0 := strings.Index(toks[7], ";SU=") + 4
					if idx0 != 3 {
						idx1 := idx0 + strings.Index(toks[7][idx0:], ";")
						vstr := toks[7][idx0:idx1]
						v, err := strconv.Atoi(vstr)
						if err == nil && v < bndSupport {
							continue
						}
					}
				}
				wb.WriteString(line)
			}
			if err == io.EOF {
				break
			}
			check(err)
		}

	}()
	return r
}

const BndSupport = 6

func Main() {
	if _, err := exec.LookPath("lumpy"); err != nil {
		shared.Slogger.Fatal("lumpy executable not found in PATH")
	}
	cli := cliargs{Processes: 3, ExcludeChroms: "hs37d5,~:,~^GL,~decoy"}
	arg.MustParse(&cli)
	runtime.GOMAXPROCS(cli.Processes)
	filter_chroms := strings.Split(strings.TrimSpace(cli.ExcludeChroms), ",")
	if cli.OutDir == "" {
		cli.OutDir = "./"
	}
	p := Lumpy(cli.Name, cli.Fasta, cli.OutDir, cli.Bams, nil, cli.Exclude, filter_chroms, !cli.NoExtraFilters)
	p.Stderr = shared.Slogger
	var err error
	ivcf, err := p.StdoutPipe()
	if err != nil {
		log.Fatal(err)
	}
	p.Start()

	vcf := bndFilter(ivcf, BndSupport, cli.Fasta)

	if cli.Genotype {
		svtyper.Svtyper(vcf, cli.Fasta, cli.Bams, cli.OutDir, cli.Name, true, cli.RemovePr)
	} else {
		path := filepath.Join(cli.OutDir, cli.Name+"-smoove.vcf.gz")
		_, wtr := bgzopen(path)
		defer shared.Slogger.Printf("wrote to %s", path)
		check(err)
		defer wtr.Close()
		io.Copy(wtr, vcf)
	}
	p.Wait()
}

func bgzopen(path string) (*os.File, *bgzf.Writer) {

	f, err := os.Create(path)
	check(err)
	w := bgzf.NewWriter(f, 1)
	return f, w

}
