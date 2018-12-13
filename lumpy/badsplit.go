package lumpy

import (
	"bytes"
	"io"
	"log"
	"os"
	"time"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
)

func isBad(counts []int8) bool {
	badEnds := 0
	nonOne := 0
	lo := 25
	if len(counts) < 50 {
		log.Println("short read")
		return false
	}

	// require first 25 and last 25 bases to be mapped quite well.
	// only allow 6 bad bases.
	hi := len(counts) - 25 - 1
	for i, c := range counts {
		if c != 1 {
			nonOne += 1
		}
		if i > lo && i < hi {
			continue
		}
		if c != 1 {
			badEnds += 1
		}
	}
	return badEnds > 5 || nonOne > 40
}

func reverse(s sam.Cigar) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}

// countBases increments bases in all cigars that are aligned.
func countBases(cigs []sam.Cigar, L int) []int8 {
	counts := make([]int8, L)
	for _, cig := range cigs {
		off := 0
		for _, op := range cig {
			if op.Type() == sam.CigarMatch || op.Type() == sam.CigarInsertion {
				for i := 0; i < op.Len(); i++ {
					counts[off+i]++
				}
			}

			off += op.Len() * op.Type().Consumes().Query
			if op.Type() == sam.CigarHardClipped {
				off += op.Len()
			}
		}
	}
	return counts
}

func getCigars(rec *sam.Record, maxL *int) []sam.Cigar {
	c := rec.Cigar
	tags, ok := rec.Tag([]byte{'S', 'A'})
	if !ok {
		panic("no sa tag")
	}
	strand := byte('+')
	if rec.Flags&sam.Reverse != 0 {
		strand = '-'
	}
	// remove SA and final ;
	b := []byte(tags[3 : len(tags)-1])
	cigs := make([]sam.Cigar, 0, 2)
	cigs = append(cigs, c)
	ref, read := c.Lengths()
	*maxL = ref
	if read > *maxL {
		*maxL = read
	}

	for _, part := range bytes.Split(b, []byte{';'}) {
		toks := bytes.Split(part, []byte{','})
		cig, err := sam.ParseCigar(toks[3])
		if toks[2][0] != strand {
			reverse(cig)
		}
		if err != nil {
			panic(err)
		}
		ref, read := cig.Lengths()
		if read > *maxL {
			*maxL = read
		}
		if ref > *maxL {
			*maxL = read
		}
		cigs = append(cigs, cig)
	}
	return cigs
}

// a bad splitter has soft clips that don't consume the entire read
// especially at the ends. e.g.
// good: 110M40S 110S40M
// bad: 110M40S 90M60S
// this function would remove `bad`. It requires
// that the outer 25 bases of each end of the read mostly
// corroborate the splitter and there can be no more than
// 40 conflicting bases. In the bad example above, there are
// 130 conflicting bases. This adjust for strand and allows
// multiple splitters (even though lumpy does not).
func badSplitter(rec *sam.Record) bool {
	var maxL int
	cigs := getCigars(rec, &maxL)
	counts := countBases(cigs, maxL)
	return isBad(counts)
}

func main() {
	// main function for testing only. called from lumpy_filter
	bampath := os.Args[1]
	removed := 0
	kept := 0
	f, err := os.Open(bampath)
	if err != nil {
		panic(err)
	}
	br, err := bam.NewReader(f, 2)
	t := time.Now()

	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		if badSplitter(rec) {
			removed++
			continue
		}

		kept++
		/*
			if err := bw.Write(rec); err != nil {
				panic(err)
			}
		*/
	}
	tot := removed + kept
	log.Printf("removed %d (%.1f%%) of %d sparse splitters in %.0f seconds", removed, 100*float64(removed)/float64(tot), tot, time.Since(t).Seconds())

	br.Close()
}
