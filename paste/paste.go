package paste

import (
	"fmt"
	"io"
	"log"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/smoove/shared"
	"github.com/brentp/xopen"
)

type cliargs struct {
	Name   string   `arg:"-n,required,help:project name used in output files."`
	OutDir string   `arg:"-o,help:output directory."`
	VCFs   []string `arg:"positional,required,help:path to vcfs."`
}

func (c *cliargs) Description() string {
	return "square VCF files from different samples with the same number of records"
}

func countNonHeaderLines(path string) int {
	f, err := xopen.Ropen(path)
	defer f.Close()
	if err != nil {
		log.Fatal(err)
	}
	count := 0
	for {
		line, err := f.ReadBytes('\n')
		if len(line) > 0 {
			if line[0] != '#' {
				count++
			}
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
	}

	return count
}

func count(wg *sync.WaitGroup, procs int, vcfs []string) {

	m := make(map[int][]string, 5)
	L := sync.Mutex{}

	ch := make(chan string, 3)
	for i := 0; i < procs; i++ {
		go func() {
			for path := range ch {
				c := countNonHeaderLines(path)
				L.Lock()
				m[c] = append(m[c], path)
				L.Unlock()
				wg.Done()
			}

		}()
	}

	for _, v := range vcfs {
		ch <- v
	}
	close(ch)

	wg.Wait()
	if len(m) != 1 {
		for k, fs := range m {
			if k == 0 {
				shared.Slogger.Printf("files: %s had %d variants", strings.Join(fs, ","), k)
				continue
			}
			if len(fs) > 5 {
				shared.Slogger.Printf("%d files had %d variants", len(fs), k)
			} else {
				shared.Slogger.Printf("files: %s had %d variants", strings.Join(fs, ","), k)
			}

		}
		shared.Slogger.Fatal("please make sure that all files have the same number of variants")
	}
	for k := range m {
		shared.Slogger.Printf("all files had %d variants", k)
	}
}

func Main() {

	cli := cliargs{OutDir: "./"}
	arg.MustParse(&cli)
	outvcf := fmt.Sprintf(filepath.Join(cli.OutDir, cli.Name) + ".smoove.square.vcf.gz")
	var wg *sync.WaitGroup

	// TODO: check files in list
	if !strings.HasSuffix(cli.VCFs[0], ".list") {
		wg = &sync.WaitGroup{}
		wg.Add(len(cli.VCFs))
		procs := 5

		go count(wg, procs, cli.VCFs)
		shared.Slogger.Printf("squaring %d files to %s", len(cli.VCFs), outvcf)
	} else {
		shared.Slogger.Printf("squaring files from %s to %s", cli.VCFs[0], outvcf)
	}

	args := []string{"merge", "-o", outvcf, "-O", "z", "--threads", "3"}
	if strings.HasSuffix(cli.VCFs[0], ".list") && len(cli.VCFs) == 1 {
		args = append(args, []string{"-l", cli.VCFs[0]}...)
	} else {
		args = append(args, cli.VCFs...)
	}

	p := exec.Command("bcftools", args...)
	p.Stderr = shared.Slogger
	p.Stdout = shared.Slogger

	if err := p.Run(); err != nil {
		log.Fatal(err)
	}
	shared.Slogger.Printf("wrote squared file to %s", outvcf)
	if wg != nil {
		wg.Wait()
	}
}
