package cnvnator

import (
	"bufio"
	"bytes"
	"html/template"
	"io"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"time"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/biogo/io/seqio/fai"
	"github.com/brentp/faidx"
	"github.com/brentp/go-athenaeum/tempclean"
	"github.com/brentp/lumpy-smoother/shared"
	"github.com/brentp/xopen"
	"github.com/pkg/errors"
)

const cnvnator_cmd = ` 
set -euo pipefail
cd {{.OutDir}}
export REF_PATH=

# split to make n chunks to use less mem.
# use STDIN to get around https://github.com/abyzovlab/CNVnator/issues/101
# this requires a specific branch of cnvnator: https://github.com/brentp/CNVnator/tree/stdin
{{range .Chroms}}
samtools view -T {{$.Reference}} -u {{$.Bam}} {{.}} | cnvnator -root {{$.Sample}}.root -chrom {{.}} -unique -tree
{{end}}

cnvnator -root {{.Sample}}.root -his {{.Bin}}
cnvnator -root {{.Sample}}.root -stat {{.Bin}}
cnvnator -root {{.Sample}}.root -partition {{.Bin}}
cnvnator -root {{.Sample}}.root -call {{.Bin}} > {{.Sample}}.cnvnator
`

func check(err error) {
	if err != nil {
		log.Fatal(err)
	}
}

type cnvargs struct {
	Sample    string
	OutDir    string
	Chroms    []string
	Bam       string
	Bin       int
	Reference string
}

func Cnvnator(bam string, name string, fasta string, outDir string, bin int, excludeChroms string) error {
	fa, err := faidx.New(fasta)
	if err != nil {
		return errors.Wrap(err, "error opening fasta")
	}
	sp := split(fa, 4, outDir, strings.Split(excludeChroms, ","))
	fa.Close()

	f := outDir + "/" + name + ".cnvnator"
	if xopen.Exists(f) {
		shared.Slogger.Println("using existing cnvnator output:", f)
		return nil
	}

	bampath, err := filepath.Abs(bam)
	if err != nil {
		return errors.Wrap(err, "error getting path")
	}

	var buf bytes.Buffer
	ca := cnvargs{Sample: name, OutDir: outDir, Reference: fasta, Chroms: sp, Bam: bampath, Bin: bin}
	t, err := template.New(name).Parse(cnvnator_cmd)
	check(err)
	check(t.Execute(&buf, ca))
	// use a bash script so we don't get an error about command too long.
	ft, err := tempclean.TempFile("", "cnvnator.sh")
	if err != nil {
		return errors.Wrap(err, "error getting temp file")
	}
	ft.Write(buf.Bytes())
	check(ft.Close())
	time.Sleep(100 * time.Millisecond)
	p := exec.Command("bash " + ft.Name())
	if err := p.Run(); err != nil {
		return err
	}
	cnvnatorToBedPe(name, outDir, fasta)
	return nil
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

func cnvnatorToBedPe(sample, outdir string, fasta string) {
	fa, err := faidx.New(fasta)
	check(err)

	order := faOrder(fa)

	f, err := os.Open(outdir + "/" + sample + ".cnvnator")
	check(err)
	defer f.Close()

	delf, err := os.Create(outdir + "/" + sample + ".del.bedpe")
	check(err)
	defer delf.Close()

	dupf, err := os.Create(outdir + "/" + sample + ".dup.bedpe")
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
	// split the chroms into n chunks that are relatively even by total bases.
	sp := make([]string, n)
	seqs := make([]fai.Record, 0, len(fa.Index))
	var tot int
	for _, v := range fa.Index {
		// exclude chromosomes with ':' in the name because these screw with cnvnator.
		if strings.Contains(v.Name, ":") {
			continue
		}
		if shared.Contains(excludeChroms, v.Name) {
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

type cliargs struct {
	Name          string `arg:"-n,required,help:project name used in output files."`
	Fasta         string `arg:"-f,required,help:fasta file."`
	OutDir        string `arg:"-o,help:output directory."`
	ExcludeChroms string `arg:"-C,help:ignore SVs with either end in this comma-delimited list of chroms. If this starts with ~ it is treated as a regular expression to exclude."`
	Bam           string `arg:"positional,required,help:path to bam to call."`
}

func Main() {
	cli := cliargs{OutDir: "", ExcludeChroms: "hs37d5,~:,~^GL,~decoy"}
	p := arg.MustParse(&cli)
	if _, err := exec.LookPath("cnvnator"); err != nil {
		p.Fail("cnvnator not found in path")
	}
	if err := Cnvnator(cli.Bam, cli.Name, cli.Fasta, cli.OutDir, 500, cli.ExcludeChroms); err != nil {
		log.Fatal(err)
	}
}
