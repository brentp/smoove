package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"text/template"
	"time"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/biogo/io/seqio/fai"
	"github.com/biogo/hts/bam"
	"github.com/brentp/faidx"
	"github.com/brentp/gargs/process"
	"github.com/brentp/go-athenaeum/shpool"
	"github.com/brentp/go-athenaeum/tempclean"
	"github.com/brentp/goleft/covstats"
	"github.com/brentp/goleft/indexcov"
	"github.com/brentp/lumpy-smoother/evidencewindow"
	"github.com/brentp/lumpy-smoother/hipstr"
	"github.com/brentp/xopen"
	"github.com/pkg/errors"
	"github.com/valyala/fasttemplate"
)

type Logger struct {
	*log.Logger
}

func (l *Logger) Write(b []byte) (int, error) {
	l.Logger.Printf(string(b))
	return len(b), nil
}

var logger *Logger

func init() {
	l := log.New(os.Stderr, "[lumpy-smoother] ", log.Ldate|log.Ltime)
	logger = &Logger{Logger: l}
}

type cliargs struct {
	Processes     int      `arg:"-p,help:number of processes to use."`
	Name          string   `arg:"-n,required,help:project name used in output files."`
	Fasta         string   `arg:"-f,required,help:fasta file."`
	OutDir        string   `arg:"-o,help:output directory."`
	MaxDepth      int      `arg:"-d,help:maximum depth in splitters/discordant file."`
	Exclude       string   `arg:"-e,help:BED of exclude regions."`
	ExcludeChroms string   `arg:"-C,help:ignore SVs with either end in this comma-delimited list of chroms. If this starts with ~ it is treated as a regular expression to exclude."`
	Bams          []string `arg:"positional,required,help:path to bams to call."`
}

func has_prog(p string) string {
	if _, err := exec.LookPath(p); err == nil {
		return "Y"
	}
	return " "
}

func (cliargs) Description() string {
	tmpl := `
lumpy-smoother: call and genotype structural variants in parallel 

lumpy-smoother calls several programs. Those with 'Y' are found on
your $PATH. Only those with '*' are required.

  [{{lumpy}}] lumpy*
  [{{lumpy_filter}}] lumpy_filter*
  [{{cnvnator}}] cnvnator [per-sample CNV calls]
  [{{samtools}}] samtools [only required for CRAM input]
  [{{mosdepth}}] mosdepth [extra filtering of split and discordant files for better scaling]
  [{{svtyper}}] svtyper [genotype SV calls]
  [{{gsort}}] gsort [(sort)  ->  compress   ->  index ]
  [{{bgzip}}] bgzip [ sort   -> (compress) ->   index ]
  [{{tabix}}] tabix [ sort   ->  compress   -> (index)]
`
	t := fasttemplate.New(tmpl, "{{", "}}")

	vars := map[string]interface{}{
		"lumpy":        has_prog("lumpy"),
		"cnvnator":     has_prog("cnvnator"),
		"lumpy_filter": has_prog("lumpy_filter"),
		"samtools":     has_prog("samtools"),
		"mosdepth":     has_prog("mosdepth"),
		"svtyper":      has_prog("svtyper"),
		"gsort":        has_prog("gsort"),
		"bgzip":        has_prog("bgzip"),
		"tabix":        has_prog("tabix"),
	}

	return t.ExecuteString(vars)
}

type filtered struct {
	split   string
	disc    string
	project string
	sample  string
	command string
	xam     string
	stats   covstats.Stats
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

type cnvargs struct {
	Sample    string
	OutDir    string
	Achroms   string
	Bchroms   string
	Cchroms   string
	Bam       string
	Bin       int
	Reference string
}

const SVTYPER_LINES = 300

const cnvnator_cmd = `
set -euo pipefail
cd {{.OutDir}}
export REF_PATH={{.OutDir}}

# split to make 3 chunks to use less mem.
# use STDIN to get around https://github.com/abyzovlab/CNVnator/issues/101
# this requires a specific branch of cnvnator: https://github.com/brentp/CNVnator/tree/stdin
samtools view -T {{.Reference}} -u {{.Bam}} {{.Achroms}} | cnvnator -root {{.Sample}}.root -chrom {{.Achroms}} -unique -tree
samtools view -T {{.Reference}} -u {{.Bam}} {{.Bchroms}} | cnvnator -root {{.Sample}}.root -chrom {{.Bchroms}} -unique -tree
samtools view -T {{.Reference}} -u {{.Bam}} {{.Cchroms}} | cnvnator -root {{.Sample}}.root -chrom {{.Cchroms}} -unique -tree

cnvnator -root {{.Sample}}.root -his {{.Bin}}
cnvnator -root {{.Sample}}.root -stat {{.Bin}}
cnvnator -root {{.Sample}}.root -partition {{.Bin}}
cnvnator -root {{.Sample}}.root -call {{.Bin}} > {{.Sample}}.cnvnator
`

func cnvnator(bam filtered, chunks []string, ref string, outdir string, bin int) string {
	f := outdir + "/" + bam.sample + ".cnvnator"
	if xopen.Exists(f) {
		logger.Println("using existing cnvnator output:", f)
		return ""
	}

	bampath, err := filepath.Abs(bam.xam)
	check(err)

	var buf bytes.Buffer
	ca := cnvargs{Sample: bam.sample, OutDir: outdir, Reference: ref, Achroms: chunks[0], Bchroms: chunks[1],
		Cchroms: chunks[2], Bam: bampath, Bin: bin}
	t, err := template.New(bam.sample).Parse(cnvnator_cmd)
	check(err)
	check(t.Execute(&buf, ca))
	// use a bash script so we don't get an error about command too long.
	ft, err := tempclean.TempFile("", "cnvnator.sh")
	check(err)
	ft.Write(buf.Bytes())
	check(ft.Close())
	time.Sleep(100 * time.Millisecond)
	return "bash " + ft.Name()
}

// cnvnator requires 1.fa, 2.fa, ... in the working dir.
func writeFa(fa *faidx.Faidx, r fai.Record, outdir string) {
	if xopen.Exists(outdir + "/" + r.Name + ".fa") {
		return
	}
	f, err := ioutil.TempFile(outdir, "ls.fa")
	check(err)
	defer os.Remove(f.Name())

	b := bufio.NewWriter(f)
	b.WriteString(">" + r.Name + "\n")
	s, err := fa.GetRaw(r.Name, 0, r.Length)
	check(err)
	b.Write(s)
	b.Write([]byte{'\n'})
	b.Flush()

	f.Close()
	if xopen.Exists(outdir + "/" + r.Name + ".fa") {
		return
	}
	os.Rename(f.Name(), outdir+"/"+r.Name+".fa")
}

func split(fa *faidx.Faidx, n int, outdir string, excludeChroms []string) []string {
	sp := make([]string, n)
	// split the chroms into 3 chunks that are relatively even by total bases.
	seqs := make([]fai.Record, 0, len(fa.Index))
	var tot int
	for _, v := range fa.Index {
		// exclude chromosomes with ':' in the name because these screw with cnvnator.
		if strings.Contains(v.Name, ":") {
			continue
		}
		if contains(excludeChroms, v.Name) {
			continue
		}
		writeFa(fa, v, outdir)
		seqs = append(seqs, v)
		tot += v.Length
	}
	sort.Slice(seqs, func(i, j int) bool { return seqs[i].Start < seqs[j].Start })
	nper := tot / n

	var k, ktot int
	for i := 0; i < n; i++ {
		for k < len(seqs) && ktot < nper {
			sp[i] += "'" + seqs[k].Name + "' "
			ktot += seqs[k].Length
			k += 1
		}
		k += 1
		ktot = 0
	}
	return sp
}

func cnvnators(bams []filtered, ref string, outdir string, bin int, excludeChroms []string, pool *shpool.Pool) {

	fa, err := faidx.New(ref)
	check(err)
	sp := split(fa, 3, outdir, excludeChroms)
	fa.Close()

	for _, b := range bams {
		if cmd := cnvnator(b, sp, ref, outdir, bin); cmd != "" {
			//log.Println(cmd)
			pool.Add(shpool.Process{Command: cmd, Prefix: "cnvnator"})
		}
	}
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

type cnv struct {
	chrom string
	start int

	end   int
	i     int
	event string
}

func cnvFromLine(line string, i int) cnv {
	toks := strings.Split(strings.TrimSpace(line), "\t")
	chrom_interval := strings.Split(toks[1], ":")
	se := strings.Split(chrom_interval[1], "-")
	start, err := strconv.Atoi(se[0])
	check(err)
	stop, err := strconv.Atoi(se[1])
	return cnv{chrom: chrom_interval[0], start: start, end: stop, i: i, event: toks[0]}
}

func (c cnv) String(breakHalf int) string {
	return strings.Join([]string{c.chrom,
		strconv.Itoa(max(1, c.start-breakHalf)),
		strconv.Itoa(c.start + breakHalf),
		c.chrom,
		strconv.Itoa(max(1, c.end-breakHalf)),
		strconv.Itoa(c.end + breakHalf),
		strconv.Itoa(c.i),
		strconv.Itoa(c.end - c.start),
		"+", "+", "TYPE:" + strings.ToUpper(c.event)}, "\t")
}

func faOrder(fa *faidx.Faidx) map[string]int {
	recs := make([]fai.Record, 0, len(fa.Index))
	for _, r := range fa.Index {
		recs = append(recs, r)
	}
	sort.Slice(recs, func(i, j int) bool { return recs[i].Start < recs[j].Start })
	names := make(map[string]int, 2)
	for i, r := range recs {
		names[r.Name] = i
	}
	return names
}

// sort according to chrom, then start, then 2nd.
func sort2(a cnv, b cnv, chromOrder map[string]int) bool {
	if chromOrder[a.chrom] != chromOrder[b.chrom] {
		return chromOrder[a.chrom] < chromOrder[b.chrom]
	}
	if a.start != b.start {
		return a.start < b.start
	}
	return a.end < b.end
}

func cnvnatorToBedPe(b filtered, outdir string, fasta string) {
	fa, err := faidx.New(fasta)
	check(err)

	order := faOrder(fa)

	f, err := os.Open(outdir + "/" + b.sample + ".cnvnator")
	check(err)
	defer f.Close()

	delf, err := os.Create(outdir + "/" + b.sample + ".del.bedpe")
	check(err)
	defer delf.Close()

	dupf, err := os.Create(outdir + "/" + b.sample + ".dup.bedpe")
	check(err)
	defer dupf.Close()

	breakHalf := 150

	br := bufio.NewReader(f)
	i := 0

	dups := make([]cnv, 0, 1048)
	dels := make([]cnv, 0, 1048)

	for {
		line, err := br.ReadString('\n')
		i += 1
		if err == io.EOF {
			break
		}
		check(err)

		c := cnvFromLine(line, i)
		if c.event == "duplication" {
			dups = append(dups, c)
		} else if c.event == "deletion" {
			dels = append(dels, c)

		} else {
			panic("expecting 'duplication' or 'deletion' as 2nd column in:" + line)
		}
	}
	sort.Slice(dels, func(i, j int) bool { return sort2(dels[i], dels[j], order) })
	sort.Slice(dups, func(i, j int) bool { return sort2(dups[i], dups[j], order) })

	for _, d := range dels {
		delf.WriteString(d.String(breakHalf))
		delf.Write([]byte{'\n'})
	}
	for _, d := range dups {
		dupf.WriteString(d.String(breakHalf))
		dupf.Write([]byte{'\n'})
	}

}

func lumpy_filter_cmd(xam string, outdir string, threads int, project string, reference string) filtered {
	// check if .split.bam and .disc.bam exist. if they do, then use.
	sm, err := indexcov.GetShortName(xam, strings.HasSuffix(xam, ".cram"))
	check(err)
	prefix := fmt.Sprintf("%s/%s", outdir, sm)

	if xopen.Exists(prefix+".split.bam") && xopen.Exists(prefix+".disc.bam") {
		return filtered{split: prefix + ".split.bam", disc: prefix + ".disc.bam", sample: sm, xam: xam, project: project}
	}

	// symlink to out dir.
	olddir := filepath.Dir(xam)
	if xopen.Exists(fmt.Sprintf("%s/%s.split.bam", olddir, sm)) && xopen.Exists(fmt.Sprintf("%s/%s.disc.bam", olddir, sm)) {
		check(os.Symlink(fmt.Sprintf("%s/%s.split.bam", olddir, sm), fmt.Sprintf("%s/%s.split.bam", outdir, sm)))
		check(os.Symlink(fmt.Sprintf("%s/%s.disc.bam", olddir, sm), fmt.Sprintf("%s/%s.disc.bam", outdir, sm)))

		return filtered{split: prefix + ".split.bam", disc: prefix + ".disc.bam", sample: sm, xam: xam, project: project}
	}

	f := filtered{split: prefix + ".split.bam", disc: prefix + ".disc.bam", sample: sm, xam: xam, project: project}
	// use .tmp.bam in case of error while running lumpy filter.
	f.command = fmt.Sprintf("lumpy_filter -f %s %s %s.tmp.bam %s.tmp.bam %d && mv %s.tmp.bam %s && mv %s.tmp.bam %s", reference, xam, f.split, f.disc, threads, f.split, f.split, f.disc, f.disc)
	return f
}

func lumpy_filters(bams []string, outdir string, project string, reference string, pool *shpool.Pool) []filtered {
	filts := make([]filtered, len(bams))
	threads := 1
	if runtime.GOMAXPROCS(0) > len(bams) {
		threads = 2
	}
	for i, b := range bams {
		filts[i] = lumpy_filter_cmd(b, outdir, threads, project, reference)
		pool.Add(shpool.Process{Command: filts[i].command, CPUs: 1, Prefix: "lumpy-filter"})
	}

	return filts
}

func processor(cmds chan string) chan bool {
	done := make(chan bool)
	go func() {
		var anyError error
		for c := range process.Runner(cmds, make(chan bool), &process.Options{}) {
			if c.ExitCode() != 0 {
				anyError = c.Err
				logger.Println(c.Err, c.String())
				break
			}
			c.Cleanup()
			// TODO: add cnvnator time

		}
		if anyError != nil {
			panic(anyError)
		}
		done <- true
		close(done)

	}()
	return done
}

func bam_stats(bams []filtered, fasta string, outdir string) {
	logger.Printf("calculating bam stats for %d bams\n", len(bams))
	var wg sync.WaitGroup

	wg.Add(2)

	f := func(mod int) {
		for i, f := range bams {
			if i%2 == mod {
				continue
			}
			br, err := NewReader(f.xam, 2, fasta)
			check(err)
			bams[i].stats = covstats.BamStats(br, 200000)
			write_hist(bams[i], outdir)
		}
		wg.Done()
	}

	go f(0)
	go f(1)
	wg.Wait()
	logger.Print("done calculating bam stats\n")
}

func (f filtered) histpath(outdir string) string {
	return fmt.Sprintf("%s/%s.histo", outdir, f.sample)
}

func write_hist(fi filtered, outdir string) {
	f, err := os.Create(fi.histpath(outdir))
	check(err)
	for i, v := range fi.stats.H {
		fmt.Fprintf(f, "%d\t%.12f\n", i, v)
	}
	f.Close()
}

type reader struct {
	io.ReadCloser
	cmd *exec.Cmd
}

func (r *reader) Close() error {
	if err := r.cmd.Wait(); err != nil {
		return errors.Wrap(err, "error closing cram reader")
	}
	return r.ReadCloser.Close()
}

// NewReader returns a bam.Reader from any path that samtools can read.
func NewReader(path string, rd int, fasta string) (*bam.Reader, error) {
	var rdr io.Reader
	if strings.HasSuffix(path, ".bam") {
		var err error
		rdr, err = os.Open(path)
		if err != nil {
			return nil, err
		}

	} else {
		cmd := exec.Command("samtools", "view", "-T", fasta, "-u", path)
		cmd.Stderr = logger
		pipe, err := cmd.StdoutPipe()
		if err != nil {
			return nil, errors.Wrap(err, "error getting stdout for process")
		}
		if err = cmd.Start(); err != nil {
			pipe.Close()
			return nil, errors.Wrap(err, "error starting process")
		}
		rdr = &reader{ReadCloser: pipe, cmd: cmd}
	}
	return bam.NewReader(rdr, rd)
}

type CS struct {
	Sample     string
	DiscPath   string
	SplitPath  string
	HistPath   string
	Mean       string
	Std        string
	ReadLength string
}

func cs_from_filter(f filtered, outdir string) CS {
	return CS{Sample: f.sample, DiscPath: f.disc, SplitPath: f.split, HistPath: f.histpath(outdir),
		Mean: fmt.Sprintf("%.3f", f.stats.TemplateMean), Std: fmt.Sprintf("%.3f", f.stats.TemplateSD),
		ReadLength: strconv.Itoa(f.stats.MaxReadLength),
	}
}

func run_lumpy(bams []filtered, fa string, outdir string, has_cnvnator bool, exclude string, name string) *exec.Cmd {
	if exclude != "" {
		exclude = "-x " + exclude
	}
	if _, err := exec.LookPath("lumpy"); err != nil {
		panic("[lumpy-smoother] lumpy not found on path")
	}

	lumpy_tmpl := "set -euo pipefail; lumpy -msw 3 -mw 4 -t $(mktemp) -tt 0 "
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
	logger.Write([]byte("wrote lumpy command to " + f.Name()))
	check(err)
	f.WriteString(cmdStr + "\n")
	f.Close()

	return exec.Command("bash", "-c", cmdStr)
}

func svtyper(vcf, outdir, reference, exclude, lib string, bams []filtered) string {

	bam_paths := make([]string, len(bams))
	for i, b := range bams {
		bam_paths[i] = b.xam
	}

	if vcf != "" {
		defer os.Remove(vcf)
		vcf = "-i " + vcf
	}

	t, err := ioutil.TempFile("", "lumpy-smoother-svtyper-tmp-")
	check(err)
	t.Close()
	cmd := fmt.Sprintf("svtyper -B %s %s --max_reads 1000 -T %s -l %s -o %s", strings.Join(bam_paths, ","), vcf, reference, lib, t.Name())

	p := exec.Command("bash", "-c", cmd)
	p.Stderr = logger
	p.Stdout = logger
	check(p.Start())
	check(p.Wait())

	return t.Name()
}

func main() {

	if len(os.Args) > 1 && os.Args[1] == "hipstr" {
		os.Args = append(os.Args[:1], os.Args[2:]...)
		hipstr.Main()
		os.Exit(0)
	}
	if len(os.Args) > 1 && os.Args[1] == "evwin" {
		os.Args = append(os.Args[:1], os.Args[2:]...)
		evidencewindow.Main()
		os.Exit(0)
	}

	here, _ := filepath.Abs(".")
	cli := cliargs{Processes: runtime.GOMAXPROCS(0), OutDir: here, MaxDepth: 500, ExcludeChroms: "hs37d5,~:,~^GL"}

	arg.MustParse(&cli)
	var err error
	cli.Fasta, err = filepath.Abs(cli.Fasta)
	check(err)
	cli.OutDir, err = filepath.Abs(cli.OutDir)
	check(err)

	if err := os.MkdirAll(cli.OutDir, 0777); err != nil {
		panic(err)
	}
	runtime.GOMAXPROCS(cli.Processes)

	pool := shpool.New(runtime.GOMAXPROCS(0), nil, &shpool.Options{LogPrefix: "[lumpy-smoother]"})

	splits := lumpy_filters(cli.Bams, cli.OutDir, cli.Name, cli.Fasta, pool)

	fmt.Fprintln(logger, "sent lumpy_filter commands")
	bam_stats(splits, cli.Fasta, cli.OutDir)

	has_cnvnator := false
	if _, err := exec.LookPath("cnvnator"); err == nil {
		has_cnvnator = true
		fmt.Fprintln(logger, "sending cnvnator commands")
		cnvnators(splits, cli.Fasta, cli.OutDir, 500, strings.Split(cli.ExcludeChroms, ","), pool)
	}
	if err := pool.Wait(); err != nil {
		panic(err)
	}
	if has_cnvnator {
		for _, b := range splits {
			cnvnatorToBedPe(b, cli.OutDir, cli.Fasta)
		}
	}
	filter_chroms := strings.Split(strings.TrimSpace(cli.ExcludeChroms), ",")
	remove_high_depths(splits, cli.MaxDepth, cli.Fasta, cli.Exclude, filter_chroms)
	logger.Print("starting lumpy")

	p := run_lumpy(splits, cli.Fasta, cli.OutDir, has_cnvnator, cli.Exclude, cli.Name)
	vcf := run_svtypers(p, cli.OutDir, cli.Fasta, cli.Exclude, splits, cli.Name)

	vcf_sort(vcf, cli.Fasta)

}

func vcf_sort(vcf string, reference string) {
	if _, err := exec.LookPath("bgzip"); err != nil {
		logger.Printf("bgzip not found. can't sort and index")
		return
	}
	if _, err := exec.LookPath("gsort"); err != nil {
		logger.Printf("gsort not found. can't sort and index")
		return
	}
	if _, err := exec.LookPath("tabix"); err != nil {
		logger.Printf("tabix not found. can't index")
	}

	tmpl := fasttemplate.New(`
set -euo pipefail
gsort {{vcf}} {{reference}}.fai | bgzip -c > {{vcf}}.tmp.vcf.gz
mv {{vcf}}.tmp.vcf.gz {{vcf}}.gz
set +e
tabix {{vcf}}.gz
set -e
rm {{vcf}}`, "{{", "}}")
	cmd := tmpl.ExecuteString(map[string]interface{}{"vcf": vcf, "reference": reference})

	p := exec.Command("bash", "-c", cmd)
	p.Stdout = logger
	p.Stderr = logger
	check(p.Run())

	logger.Printf("sorted and indexed final vcf at: %s", vcf+".gz")

}

func writeTmp(header []string, lines []string) string {
	f, err := xopen.Wopen("tmp:lumpy-smoother-tmp")
	check(err)
	for _, h := range header {
		f.WriteString(h)
	}
	for _, l := range lines {
		_, err := f.WriteString(l)
		check(err)
	}
	f.Close()
	return f.Name()
}

// cant orphan BNDs or svtyper won't genotype.
func isFirstBnd(line string) bool {
	if !strings.Contains(line, "SVTYPE=BND") {
		return false
	}
	toks := strings.Split(line, "\t")
	return strings.HasSuffix(toks[2], "_1")
}

func fixReference(line string, fa *faidx.Faidx) string {
	toks := strings.SplitN(line, "\t", 5)
	pos, err := strconv.Atoi(toks[1])
	check(err)

	r, err := fa.At(toks[0], pos-1)
	check(err)
	toks[3] = string(r)

	return strings.Join(toks, "\t")
}

type lockedWriter struct {
	mu *sync.Mutex
	*xopen.Writer
	i int
}

func run_svtypers(p *exec.Cmd, outdir string, fasta string, exclude string, bams []filtered, name string) string {
	// make n svtyper workers
	vcfch := make(chan string, 16)

	fa, err := faidx.New(fasta)
	check(err)
	defer fa.Close()

	lumpyf, err := xopen.Wopen(outdir + "/" + name + "-lumpy.vcf")
	check(err)
	fmt.Fprintf(logger, "writing lumpy output to:"+lumpyf.Name()+"\n")
	defer lumpyf.Close()
	var wg sync.WaitGroup
	lib := outdir + "/" + name + "-svtyper.lib"

	hasSVTyper := false
	var svf lockedWriter
	if _, err := exec.LookPath("svtyper"); err == nil {
		hasSVTyper = true
		// processing is done, now write the final output.
		tmp, err := xopen.Wopen(outdir + "/" + name + "-svtyper.vcf")
		check(err)
		svf = lockedWriter{mu: &sync.Mutex{}, Writer: tmp}
		defer tmp.Close()

		fmt.Fprintln(logger, "writing svtyper output to:"+svf.Name()+"\n")
	} else {
		fmt.Fprint(logger, "svtyper not found on path, not genotyping\n")
	}

	// need to first run svtyper once to calculate library stats
	// we do this in the background while lumpy starts, but we use
	// libwg to make sure that no other svtyper processes start until
	// this one is finished.
	os.Remove(lib)
	var libwg sync.WaitGroup
	if hasSVTyper {
		libwg.Add(1)
		go func() {
			tmp := svtyper("", outdir, fasta, exclude, lib, bams)
			os.Remove(tmp)
			libwg.Done()
		}()
	}

	// do the actuall svtyping in parallel.
	procs := max(1, runtime.GOMAXPROCS(0)-1)
	if procs > 20 {
		procs -= 4
	}
	for i := 0; i < procs; i++ {
		wg.Add(1)
		go func() {
			libwg.Wait()
			for vcf := range vcfch {
				svtvcf := svtyper(vcf, outdir, fasta, exclude, lib, bams)
				svf.mu.Lock()
				writeSVT(svf, svtvcf, svf.i)
				svf.i++
				svf.mu.Unlock()
			}
			wg.Done()
		}()
	}

	stdout, err := p.StdoutPipe()
	check(err)
	go func() {
		bs := bufio.NewReaderSize(stdout, 2048)
		header := make([]string, 0, 200)
		chunk := make([]string, 0, SVTYPER_LINES+1)
		for {
			line, err := bs.ReadString('\n')
			if line != "" {
				// change to actual reference from 'N'.

				if line[0] == '#' {
					header = append(header, line)
					_, errx := lumpyf.WriteString(line)
					check(errx)
					continue
				}
				line = fixReference(line, fa)
				_, errx := lumpyf.WriteString(line)
				check(errx)
				chunk = append(chunk, line)
				if len(chunk) >= SVTYPER_LINES && !isFirstBnd(chunk[len(chunk)-1]) {
					if hasSVTyper {
						vcfch <- writeTmp(header, chunk)
					}
					chunk = chunk[:0]
				}

			}
			if err == io.EOF {
				break
			}
			check(err)

		}
		if len(chunk) > 0 {
			if hasSVTyper {
				vcfch <- writeTmp(header, chunk)
			}
			chunk = chunk[:0]
		}
		close(vcfch)

	}()

	check(p.Start())
	wg.Wait()
	if err := p.Wait(); err != nil {
		log.Fatal(err)
	}
	if hasSVTyper {
		return svf.Name()
	}
	return lumpyf.Name()
}

func writeSVT(w io.Writer, vcf string, i int) {

	defer os.Remove(vcf)
	br, err := xopen.Ropen(vcf)
	check(err)
	for {
		line, err := br.ReadString('\n')
		if line != "" {
			if (line[0] == '#' && i == 0) || line[0] != '#' {
				// leave in hom-ref for now so it doesn't look like missing data.
				//if line[0] == '#' || strings.Contains(line, "\t0/1:") || strings.Contains(line, "\t1/1:") {
				fmt.Fprint(w, line)
				//}
			}
		}
		if err == io.EOF {
			break
		}
		check(err)
	}

}
