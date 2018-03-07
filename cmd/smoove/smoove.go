package main

import (
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"

	"github.com/brentp/smoove"
	"github.com/brentp/smoove/cnvnator"
	"github.com/brentp/smoove/lumpy"
	"github.com/brentp/smoove/merge"
	"github.com/brentp/smoove/paste"
	"github.com/brentp/smoove/shared"
	"github.com/brentp/smoove/svtyper"
	"github.com/valyala/fasttemplate"
)

type progPair struct {
	help string
	main func()
}

var progs = map[string]progPair{
	"cnvnator": progPair{"call cnvnator and make bedpe files needed by lumpy", cnvnator.Main},
	"genotype": progPair{"parallelize svtyper on an input VCF", svtyper.Main},
	"call":     progPair{"call lumpy (and optionally svtyper) after filtering bams", lumpy.Main},
	"merge":    progPair{"merge and sort (using svtools) calls from multiple samples", merge.Main},
	"paste":    progPair{"square final calls from multiple samples (each with same number of variants)", paste.Main},
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

  [{{cnvnator}}] cnvnator [per-sample CNV calls]
  [{{mosdepth}}] mosdepth [extra filtering of split and discordant files for better scaling]
  [{{svtyper}}] svtyper [required to genotype SV calls]
  [{{svtools}}] only needed for large cohorts.

Available sub-commands are below. Each can be run with -h for additional help.

`
	t := fasttemplate.New(tmpl, "{{", "}}")

	vars := map[string]interface{}{
		"version":      smoove.Version,
		"lumpy":        shared.HasProg("lumpy"),
		"cnvnator":     shared.HasProg("cnvnator"),
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
	for k := range progs {
		keys = append(keys, k)
		if len(k) > l {
			l = len(k)
		}
	}
	fmtr := "%-" + strconv.Itoa(l) + "s : %s\n"
	sort.Strings(keys)
	for _, k := range keys {
		fmt.Fprintf(wtr, fmtr, k, progs[k].help)

	}
	os.Exit(1)

}

func main() {

	if len(os.Args) < 2 {
		printProgs()
	}
	var p progPair
	var ok bool
	if p, ok = progs[os.Args[1]]; !ok {
		printProgs()
	}
	// remove the prog name from the call
	os.Args = append(os.Args[:1], os.Args[2:]...)
	p.main()
}
