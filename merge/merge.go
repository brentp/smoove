package merge

import (
	"log"
	"os"
	"os/exec"
	"path/filepath"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/lumpy-smoother/shared"
	"github.com/brentp/xopen"
)

type cliargs struct {
	Name   string `arg:"-n,required,help:project name used in output files."`
	OutDir string `arg:"-o,help:output directory."`

	VCFs []string `arg:"positional,required,help:path to vcfs."`
}

func Main() {

	cli := cliargs{OutDir: "./"}
	arg.MustParse(&cli)

	f, err := xopen.Wopen(filepath.Join(cli.OutDir, cli.Name) + ".lmerge.vcf")
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
	args = []string{"lmerge", "-p", "0.05", "--sum", "-i", f.Name()}
	p = exec.Command("svtools", args...)
	p.Stderr = shared.Slogger
	f, err = xopen.Wopen(filepath.Join(cli.OutDir, cli.Name) + ".lsort.sites.vcf")
	if err != nil {
		panic(err)
	}
	p.Stdout = f
	if err := p.Run(); err != nil {
		log.Fatal(err)
	}
	f.Close()
	shared.Slogger.Printf("wrote sites file to %s", f.Name())
}
