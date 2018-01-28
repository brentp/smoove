package main

import (
	"bufio"
	"bytes"
	"fmt"
	"html/template"
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

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/biogo/io/seqio/fai"
	"github.com/biogo/hts/bam"
	"github.com/brentp/faidx"
	"github.com/brentp/gargs/process"
	"github.com/brentp/goleft/covstats"
	"github.com/brentp/goleft/indexcov"
	"github.com/brentp/xopen"
	"github.com/valyala/fasttemplate"
)

type cliargs struct {
	Processes int      `arg:"-p,help:number of processes to use."`
	Name      string   `arg:"-n,required,help:project name used in output files."`
	Fasta     string   `arg:"-f,required,help:fasta file."`
	OutDir    string   `arg:"-o,help:output directory."`
	MaxDepth  int      `arg:"-d,help:maximum depth in splitters/discordant file."`
	Exclude   string   `arg:"-e,help:BED of exclude regions."`
	Bams      []string `arg:"positional,required,help:path to bams to call."`
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
  [{{cnvnator}}] cnvnator
  [{{samtools}}] samtools
  [{{mosdepth}}] mosdepth
`
	t := fasttemplate.New(tmpl, "{{", "}}")

	vars := map[string]interface{}{
		"lumpy":        has_prog("lumpy"),
		"cnvnator":     has_prog("cnvnator"),
		"lumpy_filter": has_prog("lumpy_filter"),
		"samtools":     has_prog("samtools"),
		"mosdepth":     has_prog("mosdepth"),
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

const SVTYPER_LINES = 150

const cnvnator_cmd = `
set -euo pipefail
export REF_PATH={{.Reference}}
cd {{.OutDir}}

# split to make 3 chunks to use less mem.
cnvnator -root {{.Sample}}.root -chrom {{.Achroms}} -unique -tree {{.Bam}}
cnvnator -root {{.Sample}}.root -chrom {{.Bchroms}} -unique -tree {{.Bam}}
cnvnator -root {{.Sample}}.root -chrom {{.Cchroms}} -unique -tree {{.Bam}}

cnvnator -root {{.Sample}}.root -his {{.Bin}}
cnvnator -root {{.Sample}}.root -stat {{.Bin}}
cnvnator -root {{.Sample}}.root -partition {{.Bin}}
cnvnator -root {{.Sample}}.root -call {{.Bin}} > {{.Sample}}.cnvnator
`

func cnvnator(bam filtered, chunks []string, ref string, outdir string, bin int) string {
	f := outdir + "/" + bam.sample + ".cnvnator"
	if xopen.Exists(f) {
		log.Println("using existing cnvnator output:", f)
		return ""
	}

	var buf bytes.Buffer
	ca := cnvargs{Sample: bam.sample, OutDir: outdir, Reference: ref, Achroms: chunks[0], Bchroms: chunks[1],
		Cchroms: chunks[2], Bam: bam.xam, Bin: bin}
	t, err := template.New(bam.sample).Parse(cnvnator_cmd)
	check(err)
	check(t.Execute(&buf, ca))
	return buf.String()
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

func split(fa *faidx.Faidx, n int, outdir string) []string {
	sp := make([]string, n)
	// split the chroms into 3 chunks that are relatively even by total bases.
	seqs := make([]fai.Record, 0, len(fa.Index))
	var tot int
	for _, v := range fa.Index {
		writeFa(fa, v, outdir)
		seqs = append(seqs, v)
		tot += v.Length
	}
	sort.Slice(seqs, func(i, j int) bool { return seqs[i].Start < seqs[j].Start })
	nper := tot / n

	var k, ktot int
	for i := 0; i < n; i++ {
		for k < len(seqs) && ktot < nper {
			sp[i] += seqs[k].Name + " "
			ktot += seqs[k].Length
			k += 1
		}
		k += 1
		ktot = 0
	}
	return sp
}

func cnvnators(bams []filtered, ref string, outdir string, bin int, cmds chan string) chan bool {

	fa, err := faidx.New(ref)
	check(err)
	sp := split(fa, 3, outdir)
	fa.Close()

	done := make(chan bool)

	go func() {

		for _, b := range bams {
			if cmd := cnvnator(b, sp, ref, outdir, bin); cmd != "" {
				//log.Println(cmd)
				cmds <- cmd
			}
		}
		done <- true
		close(done)
	}()
	return done
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func cnvnatorToBedPe(b filtered, outdir string) {
	f, err := os.Open(outdir + "/" + b.sample + ".cnvnator")
	check(err)
	defer f.Close()

	delf, err := os.Create(outdir + "/" + b.sample + ".del.bedpe")
	check(err)
	defer delf.Close()

	dupf, err := os.Create(outdir + "/" + b.sample + ".dup.bedpe")
	check(err)
	defer dupf.Close()

	breakHalf := 125

	br := bufio.NewReader(f)
	i := 0

	for {
		line, err := br.ReadString('\n')
		i += 1
		if err == io.EOF {
			break
		}
		check(err)
		toks := strings.Split(strings.TrimSpace(line), "\t")
		chrom_interval := strings.Split(toks[1], ":")
		se := strings.Split(chrom_interval[1], "-")
		start, err := strconv.Atoi(se[0])
		check(err)
		stop, err := strconv.Atoi(se[1])

		var f *os.File
		if toks[0] == "duplication" {
			f = dupf
		} else {
			if toks[0] != "deletion" {
				panic("expecting 'DUPLICATION' or 'DELETION' as 2nd column in:" + line)
			}
			f = delf
		}
		f.WriteString(
			strings.Join([]string{chrom_interval[0],
				strconv.Itoa(max(1, start-breakHalf)),
				strconv.Itoa(start + breakHalf),
				chrom_interval[0],
				strconv.Itoa(max(1, stop-breakHalf)),
				strconv.Itoa(stop + breakHalf),
				strconv.Itoa(i),
				strconv.Itoa(stop - start),
				"+", "+", "TYPE:" + strings.ToUpper(toks[0])}, "\t"))
		f.Write([]byte{'\n'})
	}

}

func lumpy_filter_cmd(xam string, outdir string, threads int, project string) filtered {
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
	f.command = fmt.Sprintf("lumpy_filter %s %s %s %d", xam, f.split, f.disc, threads)
	return f
}

func lumpy_filters(bams []string, outdir string, project string, cmds chan string) (chan bool, []filtered) {
	filts := make([]filtered, len(bams))
	threads := 1
	if runtime.GOMAXPROCS(0) > len(bams) {
		threads = 2
	}
	done := make(chan bool)
	for i, b := range bams {
		filts[i] = lumpy_filter_cmd(b, outdir, threads, project)
	}

	go func() {
		for i := range bams {
			if filts[i].command != "" {
				cmds <- filts[i].command
			}
		}
		done <- true
		close(done)
	}()
	return done, filts
}

func processor(cmds chan string) chan bool {
	done := make(chan bool)
	go func() {
		var anyError error
		for c := range process.Runner(cmds, make(chan bool), &process.Options{}) {
			if c.ExitCode() != 0 {
				anyError = c.Err
				log.Println(c.Err)
				break
			}
			c.Cleanup()

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
	os.Stderr.Write([]byte(fmt.Sprintf("[lumpy-smoother] calculating bam stats for %d bams\n", len(bams))))
	for i, f := range bams {
		br, err := NewReader(f.xam, 2, fasta)
		check(err)
		bams[i].stats = covstats.BamStats(br, 200000)
		write_hist(bams[i], outdir)
	}
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
		return err
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
		cmd.Stderr = os.Stderr
		pipe, err := cmd.StdoutPipe()
		if err != nil {
			return nil, err
		}
		if err = cmd.Start(); err != nil {
			pipe.Close()
			return nil, err
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

	lumpy_tmpl := fmt.Sprintf("set -euo pipefail; lumpy -msw 3 -mw 4 -t $(mktemp) -tt 0 -P %s ", exclude)
	pe_tmpl := "-pe id:{{.Sample}},bam_file:{{.DiscPath}},histo_file:{{.HistPath}},mean:{{.Mean}},stdev:{{.Std}},read_length:{{.ReadLength}},min_non_overlap:{{.ReadLength}},discordant_z:4,back_distance:30,weight:1,min_mapping_threshold:20 "
	sr_tmpl := "-sr id:{{.Sample}},bam_file:{{.SplitPath}},back_distance:10,weight:1,min_mapping_threshold:20 "

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
	os.Stderr.WriteString("[lumpy-smoother] wrote lumpy command to " + f.Name() + "\n")
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
	p.Stderr = os.Stderr
	p.Stdout = os.Stdout
	check(p.Start())
	check(p.Wait())

	return t.Name()
}

func main() {
	here, _ := filepath.Abs(".")
	cli := cliargs{Processes: runtime.GOMAXPROCS(0), OutDir: here, MaxDepth: 500}

	arg.MustParse(&cli)
	if err := os.MkdirAll(cli.OutDir, 0777); err != nil {
		panic(err)
	}
	runtime.GOMAXPROCS(cli.Processes)

	cmds := make(chan string)
	pdone := processor(cmds) // processor executes the bash commands as it gets them.

	lumpy_filter_done, splits := lumpy_filters(cli.Bams, cli.OutDir, cli.Name, cmds)
	fmt.Fprintln(os.Stderr, "[lumpy-smoother] sent lumpy_filter commands")
	bam_stats(splits, cli.Fasta, cli.OutDir)

	has_cnvnator := false
	if _, err := exec.LookPath("cnvnator"); err == nil {
		has_cnvnator = true
		fmt.Fprintln(os.Stderr, "[lumpy-smoother] sending cnvnator commands")
		<-cnvnators(splits, cli.Fasta, cli.OutDir, 500, cmds)
	}
	<-lumpy_filter_done
	close(cmds)
	<-pdone
	if has_cnvnator {
		for _, b := range splits {
			cnvnatorToBedPe(b, cli.OutDir)
		}
	}
	remove_high_depths(splits, cli.MaxDepth)

	p := run_lumpy(splits, cli.Fasta, cli.OutDir, has_cnvnator, cli.Exclude, cli.Name)
	run_svtypers(p, cli.OutDir, cli.Fasta, cli.Exclude, splits, cli.Name)

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

func run_svtypers(p *exec.Cmd, outdir string, fasta string, exclude string, bams []filtered, name string) {
	// make n svtyper workers
	vcfch := make(chan string, 1)
	svtvcfs := make([]string, 0, 20)

	var mu sync.Mutex
	fa, err := faidx.New(fasta)
	check(err)
	defer fa.Close()

	lumpyf, err := xopen.Wopen(outdir + "/" + name + "-lumpy.vcf")
	check(err)
	os.Stderr.WriteString("[lumpy-smoother] writing lumpy output to:" + lumpyf.Name() + "\n")
	defer lumpyf.Close()
	var wg sync.WaitGroup
	lib := outdir + "/" + name + "-svtyper.lib"

	hasSVTyper := false
	var fsv *xopen.Writer
	if _, err := exec.LookPath("svtyper"); err == nil {
		hasSVTyper = true
		// processing is done, now write the final output.
		var err error
		fsv, err = xopen.Wopen(outdir + "/" + name + "-svtyper.vcf")
		check(err)
		defer fsv.Close()

		os.Stderr.WriteString("[lumpy-smoother] writing svtyper output to:" + fsv.Name() + "\n")
	} else {
		os.Stderr.WriteString("[lumpy-smoother] svtyper not found on path, not genotyping\n")

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
	for i := 0; i < max(1, runtime.GOMAXPROCS(0)-1); i++ {
		wg.Add(1)
		go func() {
			libwg.Wait()
			for vcf := range vcfch {
				svvcf := svtyper(vcf, outdir, fasta, exclude, lib, bams)
				mu.Lock()
				svtvcfs = append(svtvcfs, svvcf)
				mu.Unlock()
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

	for i, vcf := range svtvcfs {
		defer os.Remove(vcf)
		br, err := xopen.Ropen(vcf)
		check(err)
		for {
			line, err := br.ReadString('\n')
			if line != "" {
				if (line[0] == '#' && i == 0) || line[0] != '#' {
					fsv.WriteString(line)
				}
			}
			if err == io.EOF {
				break
			}
			check(err)
		}
	}
}
