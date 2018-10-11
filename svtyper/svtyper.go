package svtyper

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/go-athenaeum/tempclean"
	"github.com/brentp/smoove/shared"
	"github.com/brentp/xopen"
)

type cliargs struct {
	Name      string   `arg:"-n,required,help:project name used in output files."`
	OutDir    string   `arg:"-o,help:output directory."`
	Fasta     string   `arg:"-f,required,help:fasta file."`
	RemovePr  bool     `arg:"-x,help:remove PRPOS and PREND tags from INFO."`
	DupHold   bool     `arg:"-d,help:run duphold on output."`
	Processes int      `arg:"-p,help:number of processors to use."`
	VCF       string   `arg:"-v,required,help:vcf to genotype (use - for stdin)."`
	Bams      []string `arg:"positional,required,help:path to bam to call."`
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func checkWrite(i int, e error) {
	if e != nil {
		panic(e)
	}
}

const chunkSize = 600

// cant orphan BNDs or svtyper won't genotype.
func isFirstBnd(line string) bool {
	if !strings.Contains(line, "SVTYPE=BND") {
		return false
	}
	toks := strings.SplitN(line, "\t", 5)
	return strings.HasSuffix(toks[2], "_1")
}

func writeTmp(header []string, lines []string) string {
	f, err := xopen.Wopen("tmp:smoove-tmp")
	check(err)
	for _, h := range header {
		checkWrite(f.WriteString(h))
	}
	for _, l := range lines {
		checkWrite(f.WriteString(l))
	}
	f.Close()
	return f.Name()
}

// Svtyper parellelizes genotyping of the vcf and writes the the writer.
// p is optional. If set, then it is assumed that vcf is nil and the stdout of p will be used as the vcf.
func Svtyper(vcf io.Reader, reference string, bam_paths []string, outdir, name string, excludeNonRef bool, removePR bool, duphold bool) {
	b := bufio.NewReader(vcf)
	header := make([]string, 0, 512)
	lines := make([]string, 0, chunkSize+1)
	ch := make(chan string, runtime.GOMAXPROCS(0))

	if !xopen.Exists(outdir) {
		os.MkdirAll(outdir, 0755)
	}
	// read directly from the vcf (or process) and send off to a channel.
	// svtyper will receive from that channel to allow for parallelization.
	go func() {
		defer close(ch)
		for {
			line, err := b.ReadString('\n')
			if len(line) > 0 {
				if line[0] == '#' {
					if strings.HasPrefix(line, "#CHROM") {
						line = strings.TrimSpace(strings.Join(strings.Split(line, "\t")[:8], "\t")) + "\n"
					}
					header = append(header, line)
					continue
				}
				toks := strings.SplitN(line, "\t", 10)[:8]
				line = strings.TrimSpace(strings.Join(toks, "\t")) + "\n"
				lines = append(lines, line)
				if len(lines) < chunkSize || (len(lines) == chunkSize && isFirstBnd(line)) {
					continue
				}
				// send chunk off for genotyping
				ch <- writeTmp(header, lines)
				// reset lines array.
				lines = lines[:0]
			}
			if err == io.EOF {
				break
			}
			check(err)
		}
		if len(lines) > 0 {
			ch <- writeTmp(header, lines)
		}
	}()

	var psort *exec.Cmd
	var si io.WriteCloser

	// need a multiwriter to write to both stdout and potentially to the gsort+bcftools process

	if !(shared.HasProg("gsort") == "Y" && shared.HasProg("bcftools") == "Y") {
		shared.Slogger.Fatal("gsort and bcftools required for svtyper")
	}
	o := filepath.Join(outdir, name) + "-smoove.genotyped.vcf.gz"
	shared.Slogger.Printf("writing sorted, indexed file to %s", o)
	exRef := ""
	if excludeNonRef {
		shared.Slogger.Printf("excluding variants with all unknown or homozygous reference genotypes")
		exRef = " -c 1"
	}
	var cmd string
	if removePR {
		if excludeNonRef {
			cmd = fmt.Sprintf("set -euo pipefail; gsort /dev/stdin %s.fai | bcftools annotate -x INFO/PRPOS,INFO/PREND -Ou | bcftools view -c 1 -Oz %s -o %s", reference, exRef, o)
		} else {
			cmd = fmt.Sprintf("set -euo pipefail; gsort /dev/stdin %s.fai | bcftools annotate -x INFO/PRPOS,INFO/PREND -Oz -o %s", reference, o)
		}
	} else {
		cmd = fmt.Sprintf("set -euo pipefail; gsort /dev/stdin %s.fai | bcftools view -O z%s -o %s", reference, exRef, o)
	}
	cmd += fmt.Sprintf("; bcftools index --threads %d %s", 3, o)
	psort = exec.Command("bash", "-c", cmd)
	psort.Stderr = shared.Slogger
	var err error
	si, err = psort.StdinPipe()
	check(err)
	out := bufio.NewWriter(si)
	defer out.Flush()
	check(psort.Start())
	var mu sync.Mutex
	var headerPrinted = false
	var lib string
	flib, err := tempclean.TempFile("", "svtype-lib")
	check(err)
	check(flib.Close())
	check(os.Remove(flib.Name()))
	lib = flib.Name()

	// run svtyper the first time to get the lib
	p := exec.Command("svtyper", "-B", strings.Join(bam_paths, ","), "-T", reference, "-l", lib, "-o", "-")
	p.Stderr = shared.Slogger
	p.Stdout = shared.Slogger
	check(p.Run())

	var wg sync.WaitGroup
	// read from the channel to svtype in parallel.
	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		wg.Add(1)
		go func() {
			for tmpf := range ch {
				f := tmpf
				t, err := tempclean.TempFile("", "smoove-svtyper-tmp-")
				check(err)
				t.Close()
				p := exec.Command("svtyper", "-i", f, "-B", strings.Join(bam_paths, ","), "--max_reads", "25000", "-T", reference, "-l", lib, "-o", t.Name())
				p.Stderr = shared.Slogger
				p.Stdout = shared.Slogger
				check(p.Run())
				// TODO: add check here to make sure n output variants is same as n input variants
				f = t.Name()
				rdr, err := xopen.Ropen(f)
				check(err)
				mu.Lock()
				if !headerPrinted {
					_, err = io.Copy(out, rdr)
					headerPrinted = true
					out.Flush()
					mu.Unlock()
					continue
				}
				// already printed header...
				// so skip past it.
				for {
					line, err := rdr.ReadString('\n')
					if len(line) != 0 {
						if line[0] == '#' {
							continue
						}
						checkWrite(out.WriteString(line))
						_, err := io.Copy(out, rdr)
						check(err)
						break
					}
					if err == io.EOF {
						break
					}
					check(err)
				}
				check(out.Flush())
				mu.Unlock()
				os.Remove(t.Name())
				os.Remove(f)
			}
			check(out.Flush())
			wg.Done()
		}()
	}
	wg.Wait()
	out.Flush()
	if psort != nil {
		check(si.Close())
		check(psort.Wait())
	}
	if duphold {
		args := []string{"duphold", "-o", o + ".tmp.vcf.gz", "-v", o, "-f", reference}
		args = append(args, bam_paths...)
		cmd := exec.Command("smoove", args...)
		cmd.Stderr = shared.Slogger
		cmd.Stdout = shared.Slogger
		if err := cmd.Run(); err != nil {
			tempclean.Fatalf(err.Error())
		}

		cmd = exec.Command("mv", o+".tmp.vcf.gz", o)
		cmd.Stderr = shared.Slogger
		cmd.Stdout = shared.Slogger
		if err := cmd.Run(); err != nil {
			tempclean.Fatalf(err.Error())
		}

		_ = os.Remove(o + ".tmp.vcf.gz.csi")

		cmd = exec.Command("bcftools", "index", "--threads", "3", o)
		cmd.Stderr = shared.Slogger
		cmd.Stdout = shared.Slogger
		if err := cmd.Run(); err != nil {
			tempclean.Fatalf(err.Error())
		}
	}
	shared.Slogger.Printf("wrote sorted, indexed file to %s", o)
}

const BndSupport = 6

func Main() {
	cli := cliargs{VCF: "-", Processes: 3}
	p := arg.MustParse(&cli)
	if _, err := exec.LookPath("svtyper"); err != nil {
		p.Fail(shared.Prefix + " svtyper not found on PATH")
	}
	rdr, err := xopen.Ropen(cli.VCF)
	check(err)
	defer rdr.Close()
	runtime.GOMAXPROCS(cli.Processes)
	Svtyper(rdr, cli.Fasta, cli.Bams, cli.OutDir, cli.Name, false, cli.RemovePr, cli.DupHold)
}
