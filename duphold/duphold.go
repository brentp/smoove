package duphold

import (
	"os/exec"
	"runtime"
	"strings"
	"sync"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/go-athenaeum/tempclean"
	"github.com/brentp/smoove/shared"
)

type cliargs struct {
	Fasta     string   `arg:"-f,required,help:fasta file."`
	VCF       string   `arg:"-v,required,help:path to input SV VCF"`
	Processes int      `arg:"-p,help:number of threads ot use."`
	OutVCF    string   `arg:"-o,required,help:path to output SV VCF"`
	Bams      []string `arg:"required,help:path to sample bams"`
}

type pair struct {
	bamPath string
	i       int
}

func Main() {

	cli := cliargs{Processes: runtime.GOMAXPROCS(0)}
	arg.MustParse(&cli)
	shared.Slogger.Printf("running duphold on %d files in %d processes", len(cli.Bams), cli.Processes)

	ch := make(chan pair)
	var wg sync.WaitGroup
	wg.Add(len(cli.Bams))

	go func() {
		for i, b := range cli.Bams {
			ch <- pair{b, i}
		}
		close(ch)
	}()

	paths := make([]string, 0, len(cli.Bams))
	ftype := "v"
	if strings.HasSuffix(cli.OutVCF, "bcf") {
		ftype = "b"
	} else if strings.HasSuffix(cli.OutVCF, "vcf.gz") {
		ftype = "z"
	}
	if len(cli.Bams) == 1 {
		// if only a single bam, use the outfile directly.
		paths = append(paths, cli.OutVCF)
	} else {
		// for > 1 bams, use tmp bcf files.

		for _ = range cli.Bams {
			t, err := tempclean.TempFile("dh", "smoove-duphold.bcf")
			if err != nil {
				panic(err)
			}
			t.Close()
			paths = append(paths, t.Name())
		}
	}

	var t = "1"
	if cli.Processes > len(cli.Bams) {
		t = "2"
		if cli.Processes > 2*len(cli.Bams) {
			t = "3"
		}
		cli.Processes = len(cli.Bams)
	}
	for i := 0; i < cli.Processes; i++ {
		go func() {
			for b := range ch {
				sample_bcf := paths[b.i]
				cmd := exec.Command("duphold", "-d", "-t", t, "-o", sample_bcf, "-f", cli.Fasta, "-b", b.bamPath, "-v", cli.VCF)
				cmd.Stderr = shared.Slogger
				cmd.Stdout = shared.Slogger
				if err := cmd.Run(); err != nil {
					tempclean.Fatalf("%s", err)
				}
				cmd = exec.Command("bcftools", "index", "--csi", "--threads", t, sample_bcf)
				cmd.Stderr = shared.Slogger
				cmd.Stdout = shared.Slogger
				if err := cmd.Run(); err != nil {
					tempclean.Fatalf("%s", err)
				}
			}
			wg.Done()
		}()
	}
	wg.Wait()
	if len(cli.Bams) > 1 {
		shared.Slogger.Printf("starting bcftools merge")
		args := []string{"merge", "--threads", "3", "-o", cli.OutVCF, "-O", ftype}
		args = append(args, paths...)

		cmd := exec.Command("bcftools", args...)
		cmd.Stderr = shared.Slogger
		cmd.Stdout = shared.Slogger
		if err := cmd.Run(); err != nil {
			tempclean.Fatalf("%s", err)
		} else {
			tempclean.Cleanup()
		}
	}
	shared.Slogger.Printf("finished duphold")
}
