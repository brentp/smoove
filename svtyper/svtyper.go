package svtyper

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"os/exec"
	"runtime"
	"strings"
	"sync"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/go-athenaeum/tempclean"
	"github.com/brentp/lumpy-smoother/shared"
	"github.com/brentp/xopen"
)

type cliargs struct {
	Name      string   `arg:"-n,required,help:project name used in output files."`
	Fasta     string   `arg:"-f,required,help:fasta file."`
	Processes int      `arg:"-p,help:number of processors to use."`
	VCF       string   `arg:"-v,required,help:vcf to genotype (use - for stdin)."`
	Bams      []string `arg:"positional,required,help:path to bam to call."`
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

const chunkSize = 300

// cant orphan BNDs or svtyper won't genotype.
func isFirstBnd(line string) bool {
	if !strings.Contains(line, "SVTYPE=BND") {
		return false
	}
	toks := strings.Split(line, "\t")
	return strings.HasSuffix(toks[2], "_1")
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

func Svtyper(vcf io.Reader, outvcf io.Writer, reference string, bam_paths []string) {

	b := bufio.NewReader(vcf)
	header := make([]string, 0, 512)
	lines := make([]string, 0, chunkSize+1)
	ch := make(chan string, runtime.GOMAXPROCS(0))
	//
	go func() {
		defer close(ch)
		for {
			line, err := b.ReadString('\n')
			if len(line) > 0 {
				if line[0] == '#' {
					header = append(header, line)
					continue
				}
				lines = append(lines, line)
				if len(lines) < chunkSize || (len(lines) == chunkSize && isFirstBnd(line)) {
					continue
				}
				ch <- writeTmp(header, lines)
			}
			if err == io.EOF {
				break
			}
			check(err)
		}
	}()
	out := bufio.NewWriter(outvcf)
	_, err := exec.LookPath("svtyper")
	hasSvtyper := err == nil
	var lib string
	if hasSvtyper {
		flib, err := tempclean.TempFile("", "svtype-lib")
		check(err)
		lib = flib.Name()
		// run svtyper the first time to get the lib
		cmd := fmt.Sprintf("svtyper -B %s -T %s -l %s", strings.Join(bam_paths, ","), reference, lib)
		p := exec.Command("bash", "-c", cmd)
		check(p.Run())
	}

	var wg sync.WaitGroup
	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		wg.Add(1)
		go func() {
			for tmpf := range ch {
				f := tmpf
				if hasSvtyper { // genotype
					t, err := tempclean.TempFile("", "lumpy-smoother-svtyper-tmp-")
					check(err)
					t.Close()
					cmd := fmt.Sprintf("svtyper -i %s -B %s %s --max_reads 1000 -T %s -l %s -o %s", f, strings.Join(bam_paths, ","), vcf, reference, lib, t.Name())

					p := exec.Command("bash", "-c", cmd)
					p.Stderr = shared.Slogger
					p.Stdout = shared.Slogger
					check(p.Run())
					f = t.Name()
				}
				rdr, err := xopen.Ropen(f)
				check(err)
				io.Copy(out, rdr)
			}
			wg.Done()
		}()
	}
	wg.Wait()
}

func Main() {
	cli := cliargs{VCF: "-", Processes: 3}
	arg.MustParse(&cli)
	rdr, err := xopen.Ropen(cli.VCF)
	check(err)
	wtr := bufio.NewWriter(os.Stdout)
	//func Svtyper(vcf io.Reader, outvcf io.Writer, reference, exclude, bam_paths []string) {
	Svtyper(rdr, wtr, cli.Fasta, cli.Bams)
}
