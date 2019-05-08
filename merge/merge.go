package merge

import (
	"encoding/json"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

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

type cliplotargs struct {
	VCF  string `arg:"-v,required,help:path to input VCF from smoove 0.2.3 or greater."`
	HTML string `arg:"-h,required,help:path to output html file to be written."`
}

type count struct {
	sample       string
	split_before int
	disc_before  int
	split_after  int
	disc_after   int
}

func mustInt(s string) int {
	r, err := strconv.Atoi(s)
	if err != nil {
		panic(err)
	}
	return r

}

const tmpl = `
<!DOCTYPE html>
<html lang="en">
<head>
<script type="text/javascript" src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
</head>
<body>
<div id="before"></div>
<div id="after"></div>
<script>

var tracebefore = {
  x: %s,
  y: %s,
  text: %s,
  mode: 'markers',
  marker: {
    size: 10,
    color: %s,
  }
};
var traceafter = {
  x: %s,
  y: %s,
  text: %s,
  mode: 'markers',
  marker: {
    size: 10,
    color: %s,
  }
};

var layout = {title: 'number of reads before filtering', hovermode: 'closest'}
layout.xaxis = {title : "split reads"}
layout.yaxis = {title : "discordant reads"}

Plotly.newPlot('before', [tracebefore], layout);
layout.title = 'number of reads after filtering';
Plotly.newPlot('after', [traceafter], layout);
</script>
Created by <a href="https://github.com/brentp/smoove">smoove</a>
</body>
</html>
`

func mustMarshal(x interface{}) string {
	v, err := json.Marshal(x)
	if err != nil {
		panic(err)
	}
	return string(v)
}

func PlotCountsMain() {
	cli := cliplotargs{}
	arg.MustParse(&cli)
	plotCounts(cli.VCF, cli.HTML)
}

func plotCounts(path string, outpath string) {
	// see: https://plot.ly/javascript/line-and-scatter/
	f, err := xopen.Ropen(path)
	defer f.Close()
	if err != nil {
		panic(err)
	}

	split_before := make([]int, 0, 16)
	split_after := make([]int, 0, 16)
	disc_before := make([]int, 0, 16)
	disc_after := make([]int, 0, 16)
	samples := make([]string, 0, 16)
	patt := regexp.MustCompile("[,:]")

	for {
		line, err := f.ReadString('\n')
		if err == io.EOF {
			break
		}
		if line[0] != '#' {
			break
		}
		if !strings.HasPrefix(line, "##smoove_count_stats=") {
			continue
		}
		// ##smoove_count_stats=sample:4494396,7960884,747932,297854
		tmp := strings.Split(line, "=")
		info := patt.Split(strings.TrimRight(tmp[1], "\r\n"), -1)
		samples = append(samples, info[0])
		split_before = append(split_before, mustInt(info[1]))
		disc_before = append(disc_before, mustInt(info[2]))
		split_after = append(split_after, mustInt(info[3]))
		disc_after = append(disc_after, mustInt(info[4]))
	}
	rng := make([]int, len(disc_after))
	for i := 0; i < len(rng); i++ {
		rng[i] = i
	}
	wtr, err := xopen.Wopen(outpath)
	if err != nil {
		panic(err)
	}

	fmt.Fprintf(wtr, tmpl, mustMarshal(split_before), mustMarshal(disc_before), mustMarshal(samples), mustMarshal(rng), mustMarshal(split_after), mustMarshal(disc_after), mustMarshal(samples), mustMarshal(rng))
	wtr.Close()
	shared.Slogger.Printf("wrote html file of disc, split counts to %s", wtr.Name())

}

func Main() {

	cli := cliargs{OutDir: "./"}
	arg.MustParse(&cli)
	shared.Slogger.Printf("merging %d files", len(cli.VCFs))

	f, err := xopen.Wopen(filepath.Join(cli.OutDir, cli.Name) + ".lsort.vcf")
	if err != nil {
		panic(err)
	}

	args := []string{"lsort", "-r", "-t", os.TempDir(), "-b", "400"}
	args = append(args, cli.VCFs...)
	shared.Slogger.Printf("finished sorting %d files; merge starting.", len(cli.VCFs))

	p := exec.Command("svtools", args...)
	p.Stderr = shared.Slogger
	p.Stdout = f

	if err := p.Run(); err != nil {
		log.Fatal(err)
	}
	f.Close()
	of := filepath.Join(cli.OutDir, cli.Name) + ".sites.vcf.gz"
	p = exec.Command("bash", "-c", fmt.Sprintf("set -euo pipefail; svtools lmerge -f 20 -i %s | grep -v '^##bcftools_viewCommand' | bgzip -c > %s", f.Name(), of))
	p.Stderr = shared.Slogger
	p.Stdout = shared.Slogger
	if err := p.Run(); err != nil {
		if strings.Contains(err.Error(), "Required tag PREND") {
			log.Println("[smoove] use e.g.: `for f in *.genotyped.vcf.gz; do echo -n $f' '; bcftools view -H $f | grep -cv PREND; done | awk '$2 != 0'` to find the bad files.")
		}
		if strings.Contains(err.Error(), "Required tag PRPOS") {
			log.Println("[smoove] use e.g.: `for f in *.genotyped.vcf.gz; do echo -n $f' '; bcftools view -H $f | grep -cv PRPOS; done | awk '$2 != 0'` to find the bad files.")
		}
		log.Fatal(err)
	}
	os.Remove(f.Name())
	shared.Slogger.Printf("wrote sites file to %s", of)
	plotCounts(of, filepath.Join(cli.OutDir, cli.Name)+".smoove-counts.html")
}
