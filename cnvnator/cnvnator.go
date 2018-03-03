package cnvnator

import (
	"bufio"
	"bytes"
	"html/template"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
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
	// TODO: convert to del/dup.bedpe
	return p.Run()
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
	arg.MustParse(&cli)
	if err := Cnvnator(cli.Bam, cli.Name, cli.Fasta, cli.OutDir, 500, cli.ExcludeChroms); err != nil {
		log.Fatal(err)
	}
}
