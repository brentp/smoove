// Package evidence window is meant for filtering splitter and discordant files used by lumpy.
// There are a lot of spurious interchromosomal splitters and discordants across the genome that increase
// lumpy run-time and memory-use. This module removes any interchromosomal reads (including by a split read)
// that don't have other evidence to the same chromosome within a given window.
package evidencewindow

import (
	"io"
	"log"
	"os"
	"runtime"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/bgzf/index"
	"github.com/biogo/hts/sam"
	"github.com/brentp/bigly/bamat"
)

type Args struct {
	Window   int            `arg:"-w,required,help:size of sliding window"`
	Evidence int            `arg:"-e,required,help:required pieces of evidence in the window"`
	Bams     []string       `arg:"positional,required,help:bams in which to call STRs"`
	ibams    []*bamat.BamAt `arg:"-"`
	obams    []*bam.Writer  `arg:"-"`
	ofiles   []*os.File     `arg:"-"`
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func (cli *Args) setBams() {
	for _, path := range cli.Bams {

		ba, err := bamat.New(path)
		check(err)

		cli.ibams = append(cli.ibams, ba)

		fw, err := os.Create(path[:len(path)-4] + ".ev.bam")
		check(err)
		bw, err := bam.NewWriterLevel(fw, ba.Header(), 1, 1)
		check(err)

		cli.ofiles = append(cli.ofiles, fw)
		cli.obams = append(cli.obams, bw)
	}

}

func getChroms(bams []*bamat.BamAt) []string {
	chroms := make([]string, 0, 20)
	for i, b := range bams {
		if i == 0 {
			for _, chrom := range b.Header().Refs() {
				chroms = append(chroms, chrom.Name())
			}
		}
		// TODO: check all chroms have the same chrom set.
	}
	return chroms
}

type Aln struct {
	*sam.Record
	wtr  *bam.Writer
	refs map[string]*sam.Reference
}

func Main() {
	cli := &Args{Window: 1000, Evidence: 2}
	arg.MustParse(cli)
	runtime.GOMAXPROCS(1)
	EvidenceRemoval(cli)
}

func EvidenceRemoval(cli *Args) []string {

	cli.setBams()

	chroms := getChroms(cli.ibams)
	var skipped, written int

	for _, chrom := range chroms {
		countsByTid := make([]int, len(chroms))
		alns := make([]Aln, 0, 4096)
		for aln := range Merge(cli.ibams, cli.obams, chrom) {

			// remove items from the start of the queue.
			for len(alns) > 0 && alns[0].End() < aln.Start()-cli.Window {

				var o Aln
				// TODO: do something here to stop memory leak from growing underlying array.
				o, alns = alns[0], alns[1:]
				// NOTE: only tracking mate ids, might want to track posns as well.
				mates := o.MateIds()
				var used bool
				for _, m := range mates {
					if countsByTid[m] >= cli.Evidence || o.MateRef.ID() == o.Ref.ID() {
						written++
						used = true
						check(o.wtr.Write(o.Record))
						break
						// remove this one from the counts
					}
				}
				for _, m := range mates {
					countsByTid[m]--
				}
				if !used {
					skipped++
				}
			}

			mates := aln.MateIds()

			// if there are no mates from different chroms, we can write immediately.
			if len(mates) == 0 {
				written++
				check(aln.wtr.Write(aln.Record))
				continue
			}

			// increment all mates.
			alns = append(alns, aln)
			for _, mate := range mates {
				countsByTid[mate]++
			}

		}
		for _, aln := range alns {
			mates := aln.MateIds()
			var used bool
			for _, m := range mates {
				if countsByTid[m] >= cli.Evidence || aln.MateRef.ID() == aln.Ref.ID() {
					written++
					check(aln.wtr.Write(aln.Record))
					used = true
					break
				}
			}
			for _, m := range mates {
				countsByTid[m]--
			}
			if !used {
				skipped++
			}
		}
	}
	for _, bw := range cli.obams {
		check(bw.Close())
	}
	for _, br := range cli.ibams {
		check(br.Close())
	}
	log.Printf("skipped %d out of %d (%.2f%%)", skipped, skipped+written, 100.0*float64(skipped)/float64(skipped+written))
	names := make([]string, len(cli.ofiles))
	for i, f := range cli.ofiles {
		names[i] = f.Name()
	}
	return names
}

const maxUint = ^uint(0)
const maxInt = int(maxUint >> 1)

func (a Aln) MateIds() []int {
	var ids []int
	if a.MateRef.ID() != a.Ref.ID() {
		ids = append(ids, a.MateRef.ID())

	}
	tag, hasTag := a.Tag([]byte{'S', 'A'})
	// TODO: need lookup from string -> id.
	if hasTag {
		v := tag.Value().(string)
		if strings.Count(v, ";") > 1 {
			for _, vs := range strings.Split(v, ";") {
				if len(vs) == 0 {
					continue
				}
				mate := vs[:strings.Index(vs, ",")]
				ids = append(ids, a.refs[mate].ID())
			}
		} else {
			mate := v[:strings.Index(v, ",")]
			ids = append(ids, a.refs[mate].ID())
		}

	}
	return ids
}

func Merge(bams []*bamat.BamAt, obams []*bam.Writer, chrom string) chan Aln {
	ch := make(chan Aln, 2)
	q := make([]*sam.Record, len(bams))
	itrs := make([]*bam.Iterator, len(bams))
	var err error

	for i, b := range bams {

		itrs[i], err = b.Query(chrom, 0, -1)
		if err == io.EOF || err == index.ErrInvalid {
			q[i] = nil
			continue
		}
		if err != nil {
			panic(err)
		}
		if !itrs[i].Next() {
			if err := itrs[i].Error(); err != nil && err != io.EOF && err != index.ErrInvalid {
				panic(err)
			}
		} else {
			q[i] = itrs[i].Record()
		}
	}

	go func() {
		for {
			start := maxInt
			end := maxInt
			mini := maxInt
			for i, aln := range q {
				if aln == nil {
					q[i] = nil
					continue
				}
				if aln.Start() < start || (aln.Start() == start && aln.End() < end) {
					start = aln.Start()
					end = aln.End()
					mini = i
				}

			}
			if start == maxInt {
				break
			}
			ch <- Aln{Record: q[mini], wtr: obams[mini], refs: bams[mini].Refs}
			if !itrs[mini].Next() {
				if err := itrs[mini].Error(); err != nil && err != io.EOF {
					panic(err)
				}
				q[mini] = nil

			} else {
				q[mini] = itrs[mini].Record()
			}

		}
		close(ch)
	}()
	return ch
}
