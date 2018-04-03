package main

import (
	"fmt"
	"io"
	"os"
	"strconv"

	"github.com/brentp/smoove"
	"github.com/brentp/smoove/annotate"
	"github.com/brentp/smoove/lumpy"
	"github.com/brentp/smoove/merge"
	"github.com/brentp/smoove/paste"
	"github.com/brentp/smoove/shared"
	"github.com/brentp/smoove/svtyper"
	"github.com/valyala/fasttemplate"
)

type progPair struct {
	name string
	help string
	main func()
}

var progs = []progPair{
	//"cnvnator": progPair{"call cnvnator and make bedpe files needed by lumpy", cnvnator.Main},
	progPair{"call", "call lumpy (and optionally svtyper)", lumpy.Main},
	progPair{"merge", "merge and sort (using svtools) calls from multiple samples", merge.Main},
	progPair{"genotype", "parallelize svtyper on an input VCF", svtyper.Main},
	progPair{"paste", "square final calls from multiple samples (each with same number of variants)", paste.Main},
	progPair{"annotate", "annotate a VCF with gene and quality of SV call", annotate.Main},
}

func Description() string {
	tmpl := `smoove version: {{version}}

smoove calls several programs. Those with 'Y' are found on your $PATH. Only those with '*' are required.

 *[{{bgzip}}] bgzip [ sort   -> (compress) ->   index ]
 *[{{gsort}}] gsort [(sort)  ->  compress   ->  index ]
 *[{{tabix}}] tabix [ sort   ->  compress   -> (index)]
 *[{{lumpy}}] lumpy
 *[{{lumpy_filter}}] lumpy_filter
 *[{{samtools}}] samtools [only required for CRAM input]

  [{{mosdepth}}] mosdepth [extra filtering of split and discordant files for better scaling]
  [{{svtyper}}] svtyper [required to genotype SV calls]
  [{{svtools}}] svtools [only needed for large cohorts].

Available sub-commands are below. Each can be run with -h for additional help.

`
	t := fasttemplate.New(tmpl, "{{", "}}")

	vars := map[string]interface{}{
		"version": smoove.Version,
		"lumpy":   shared.HasProg("lumpy"),
		//"cnvnator":     shared.HasProg("cnvnator"),
		"lumpy_filter": shared.HasProg("lumpy_filter"),
		"samtools":     shared.HasProg("samtools"),
		"mosdepth":     shared.HasProg("mosdepth"),
		"svtyper":      shared.HasProg("svtyper"),
		"svtools":      shared.HasProg("svtools"),
		"gsort":        shared.HasProg("gsort"),
		"bgzip":        shared.HasProg("bgzip"),
		"tabix":        shared.HasProg("tabix"),
	}
	return t.ExecuteString(vars)
}

func printProgs() {

	var wtr io.Writer = os.Stdout

	fmt.Fprintf(wtr, Description())
	var keys []string
	l := 5
	for _, p := range progs {
		keys = append(keys, p.name)
		if len(p.name) > l {
			l = len(p.name)
		}
	}
	fmtr := "%-" + strconv.Itoa(l) + "s : %s\n"

	for _, p := range progs {
		fmt.Fprintf(wtr, fmtr, p.name, p.help)
	}
	os.Exit(1)

}

func get(name string) (*progPair, bool) {
	for _, p := range progs {
		if p.name == name {
			return &p, true
		}
	}
	return nil, false
}

func main() {

	if len(os.Args) < 2 {
		printProgs()
	}
	var p *progPair
	var ok bool
	if p, ok = get(os.Args[1]); !ok {
		printProgs()
	}
	// remove the prog name from the call
	os.Args = append(os.Args[:1], os.Args[2:]...)
	shared.Slogger.Printf("starting with version %s", smoove.Version)
	(*p).main()
}
