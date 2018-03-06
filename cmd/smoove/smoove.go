package main

import (
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"

	"github.com/brentp/goleft"
	"github.com/brentp/smoove/cnvnator"
	"github.com/brentp/smoove/lumpy"
	"github.com/brentp/smoove/merge"
	"github.com/brentp/smoove/svtyper"
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
}

func printProgs() {

	var wtr io.Writer = os.Stdout

	fmt.Fprintf(wtr, "smoove Version: %s\n\n", goleft.Version)
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
