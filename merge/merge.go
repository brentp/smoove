package merge

import (
	"fmt"
	"log"
	"os"
	"os/exec"
	"path/filepath"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/smoove/shared"
	"github.com/brentp/xopen"
)

type cliargs struct {
	Name   string   `arg:"-n,required,help:project name used in output files."`
	OutDir string   `arg:"-o,help:output directory."`
	Fasta  string   `arg:"-f,required,help:fasta file."`
	VCFs   []string `arg:"positional,required,help:path to vcfs."`
}

func Main() {

	cli := cliargs{OutDir: "./"}
	arg.MustParse(&cli)
	shared.Slogger.Printf("merging %d files", len(cli.VCFs))

	f, err := xopen.Wopen(filepath.Join(cli.OutDir, cli.Name) + ".lsort.vcf")
	if err != nil {
		panic(err)
	}

	args := []string{"lsort", "-r", "-t", os.TempDir(), "-b", "1000"}
	args = append(args, cli.VCFs...)

	p := exec.Command("svtools", args...)
	p.Stderr = shared.Slogger
	p.Stdout = f

	if err := p.Run(); err != nil {
		log.Fatal(err)
	}
	f.Close()
	of := filepath.Join(cli.OutDir, cli.Name) + ".sites.vcf.gz"
	p = exec.Command("bash", "-c", fmt.Sprintf("svtools lmerge -f 4 -p 0.05 --sum -i %s | grep -v '^##bcftools_viewCommand' | bgzip -c > %s", f.Name(), of))
	p.Stderr = shared.Slogger
	p.Stdout = shared.Slogger
	if err := p.Run(); err != nil {
		log.Fatal(err)
	}
	os.Remove(f.Name())
	shared.Slogger.Printf("wrote sites file to %s", of)
}
