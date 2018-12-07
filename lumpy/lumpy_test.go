package lumpy

import (
	"os"
	"testing"

	"github.com/biogo/hts/bam"

	. "gopkg.in/check.v1"
)

func TestStartEndFix(t *testing.T) {
	in := "1	1116265	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116258"
	out := fixStartEnd(in)
	if out != "1	1116258	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116265" {
		t.Errorf("didn't switch")
	}

	in = "1	1116265	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116258;"
	if fixStartEnd(in) != "1	1116258	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116265;" {
		t.Errorf("didn't switch")
	}

	in = "1	1116258	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116265;"
	if fixStartEnd(in) != in {
		t.Errorf("unneeded switch")
	}

	in = "1	1116265	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=111;"
	if fixStartEnd(in) != "1	111	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116265;" {
		t.Errorf("improper switch: %s", fixStartEnd(in))
	}

	in = "1	1116265	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=111"
	if fixStartEnd(in) != "1	111	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116265" {
		t.Errorf("improper switch: %s", fixStartEnd(in))
	}

}

func Test(t *testing.T) { TestingT(t) }

type LumpyTest struct{}

var _ = Suite(&LumpyTest{})

func (s *LumpyTest) TestNM(c *C) {
	f, err := os.Open("t.bam")
	c.Assert(err, IsNil)
	br, err := bam.NewReader(f, 1)
	c.Assert(err, IsNil)
	defer br.Close()

	r, err := br.Read()
	c.Assert(nm(r), Equals, 10)
	c.Assert(err, IsNil)

	// nm_above corrects for insertion, deletion events.
	c.Assert(nm_above(r, 5), Equals, false)
	c.Assert(nm_above(r, 4), Equals, false)
	c.Assert(nm_above(r, 3), Equals, true)

	r, err = br.Read()
	c.Assert(err, IsNil)
	c.Assert(nm(r), Equals, 0)

}

type InterChromosomalTest struct {
	data []inter
}

var _ = Suite(&InterChromosomalTest{})

func (s *InterChromosomalTest) SetUpSuite(c *C) {
	s.data = []inter{
		inter{read1_tid: 1, read2_tid: 6, read1_pos: 5156100, read2_pos: 154944085, qname: "a"},
		inter{read1_tid: 1, read2_tid: 5, read1_pos: 5160000, read2_pos: 106205130, qname: "b"},
		inter{read1_tid: 1, read2_tid: 18, read1_pos: 5192900, read2_pos: 10905578, qname: "c"},
		inter{read1_tid: 1, read2_tid: 23, read1_pos: 5208900, read2_pos: 50, qname: "d"},
		inter{read1_tid: 1, read2_tid: 8, read1_pos: 5467400, read2_pos: 149095348, qname: "e"},
		inter{read1_tid: 2, read2_tid: 8, read1_pos: 10000000, read2_pos: 149095348, qname: "f"},
		inter{read1_tid: 2, read2_tid: 8, read1_pos: 10000100, read2_pos: 149106348, qname: "g"},
	}
}

func (s *InterChromosomalTest) TestLocalRightOutOfWindow(c *C) {
	c.Assert(localRightOutOfWindow(s.data[1], s.data[5]), Equals, true)
	c.Assert(localRightOutOfWindow(s.data[5], s.data[1]), Equals, false)
	c.Assert(localRightOutOfWindow(s.data[0], s.data[1]), Equals, false)
	c.Assert(localRightOutOfWindow(s.data[0], s.data[4]), Equals, true)
	c.Assert(localRightOutOfWindow(s.data[0], s.data[0]), Equals, false)
}

func (s *InterChromosomalTest) TestLocalLeftOutOfWindow(c *C) {
	c.Assert(localLeftOutOfWindow(s.data[5], s.data[1]), Equals, true)
	c.Assert(localLeftOutOfWindow(s.data[1], s.data[5]), Equals, false)
	c.Assert(localLeftOutOfWindow(s.data[1], s.data[0]), Equals, false)
	c.Assert(localLeftOutOfWindow(s.data[4], s.data[0]), Equals, true)
	c.Assert(localLeftOutOfWindow(s.data[0], s.data[0]), Equals, false)
}

func (s *InterChromosomalTest) TestDistantRightOutOfWindow(c *C) {
	c.Assert(distantRightOutOfWindow(s.data[1], s.data[0]), Equals, true)
	c.Assert(distantRightOutOfWindow(s.data[0], s.data[1]), Equals, false)
	c.Assert(distantRightOutOfWindow(s.data[4], s.data[5]), Equals, false)
	c.Assert(distantRightOutOfWindow(s.data[5], s.data[6]), Equals, true)
	c.Assert(distantRightOutOfWindow(s.data[0], s.data[0]), Equals, false)
}

func (s *InterChromosomalTest) TestDistantLeftOutOfWindow(c *C) {
	c.Assert(distantLeftOutOfWindow(s.data[0], s.data[1]), Equals, true)
	c.Assert(distantLeftOutOfWindow(s.data[1], s.data[0]), Equals, false)
	c.Assert(distantLeftOutOfWindow(s.data[5], s.data[4]), Equals, false)
	c.Assert(distantLeftOutOfWindow(s.data[6], s.data[5]), Equals, true)
	c.Assert(distantLeftOutOfWindow(s.data[0], s.data[0]), Equals, false)
}

func (s *InterChromosomalTest) TestInterChromFilter(c *C) {
	counts := make(map[string]int)
	interchromfilter(counts, s.data[:2])

	var expected_counts = map[string]int{
		"a": -1,
		"b": -1,
	}
	c.Assert(counts["a"], Equals, expected_counts["a"])
	c.Assert(counts["b"], Equals, expected_counts["b"])

	retain := []inter{
		inter{read1_tid: 1, read2_tid: 6, read1_pos: 51561, read2_pos: 154944085, qname: "a"},
		inter{read1_tid: 1, read2_tid: 5, read1_pos: 51800, read2_pos: 106205130, qname: "b"},
		inter{read1_tid: 1, read2_tid: 18, read1_pos: 51829, read2_pos: 10905578, qname: "c"},
		inter{read1_tid: 1, read2_tid: 5, read1_pos: 52089, read2_pos: 106205110, qname: "d"},
	}
	counts = make(map[string]int)
	interchromfilter(counts, retain)
	c.Assert(counts["a"], Equals, -1)
	c.Assert(counts["c"], Equals, -1)
	c.Assert(counts["b"], Equals, 0)
	c.Assert(counts["d"], Equals, 0)
}

// Edge cases
func (s *InterChromosomalTest) TestInterChromFilterEdge1(c *C) {
	counts := make(map[string]int)
	edge1 := []inter{
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 5156100, read2_pos: 5156100, qname: "a"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 5160000, read2_pos: 5160000, qname: "b"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 5169000, read2_pos: 5169000, qname: "c"},
	}
	interchromfilter(counts, edge1)
	c.Assert(counts["a"], Equals, 0)
	c.Assert(counts["b"], Equals, 0)
	c.Assert(counts["c"], Equals, 0)
}

func (s *InterChromosomalTest) TestInterChromFilterEdge2(c *C) {
	counts := make(map[string]int)
	edge2 := []inter{
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 5156100, read2_pos: 5156100, qname: "a"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 5160000, read2_pos: 5169000, qname: "b"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 5169000, read2_pos: 5160000, qname: "c"},
	}
	interchromfilter(counts, edge2)
	c.Assert(counts["a"], Equals, -1)
	c.Assert(counts["b"], Equals, 0)
	c.Assert(counts["c"], Equals, 0)
}

func (s *InterChromosomalTest) TestInterChromFilterEdge3(c *C) {
	// test a run of in-window reads
	counts := make(map[string]int)
	edge3 := []inter{
		inter{read1_tid: 1, read2_tid: 12, read1_pos: 506000, read2_pos: 601020, qname: "0"},
		inter{read1_tid: 1, read2_tid: 22, read1_pos: 508000, read2_pos: 601020, qname: "1"},
		inter{read1_tid: 1, read2_tid: 24, read1_pos: 510000, read2_pos: 601020, qname: "2"},
		inter{read1_tid: 1, read2_tid: 24, read1_pos: 512000, read2_pos: 610020, qname: "3"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 514610, read2_pos: 415610, qname: "a0"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 515610, read2_pos: 515610, qname: "a"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 515610, read2_pos: 415620, qname: "b0"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 515620, read2_pos: 515610, qname: "b"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 515630, read2_pos: 515610, qname: "c"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 516000, read2_pos: 515610, qname: "d"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 517000, read2_pos: 517000, qname: "e"},
		inter{read1_tid: 1, read2_tid: 15, read1_pos: 519990, read2_pos: 600000, qname: "f"},
		inter{read1_tid: 1, read2_tid: 12, read1_pos: 519990, read2_pos: 600000, qname: "g"},
		inter{read1_tid: 1, read2_tid: 18, read1_pos: 519990, read2_pos: 600000, qname: "h"},
		inter{read1_tid: 1, read2_tid: 12, read1_pos: 522000, read2_pos: 600020, qname: "i"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 531990, read2_pos: 515610, qname: "j"},
		inter{read1_tid: 1, read2_tid: 12, read1_pos: 532000, read2_pos: 601020, qname: "k"},
		inter{read1_tid: 1, read2_tid: 22, read1_pos: 533000, read2_pos: 601020, qname: "l"},
		inter{read1_tid: 1, read2_tid: 24, read1_pos: 534000, read2_pos: 601020, qname: "m"},
		inter{read1_tid: 1, read2_tid: 24, read1_pos: 536000, read2_pos: 611020, qname: "n"},
	}
	interchromfilter(counts, edge3)
	c.Assert(counts["a0"], Equals, 0)
	c.Assert(counts["b0"], Equals, 0)
	c.Assert(counts["a"], Equals, 0)
	c.Assert(counts["b"], Equals, 0)
	c.Assert(counts["c"], Equals, 0)
	c.Assert(counts["d"], Equals, 0)
	c.Assert(counts["e"], Equals, 0)
	c.Assert(counts["f"], Equals, -1)
	c.Assert(counts["g"], Equals, 0)
	c.Assert(counts["h"], Equals, -1)
	c.Assert(counts["i"], Equals, 0)
	c.Assert(counts["j"], Equals, -1)
	c.Assert(counts["k"], Equals, 0)
}

func (s *InterChromosomalTest) TestInterChromFilterEdge4(c *C) {
	// test a run of in-window reads
	counts := make(map[string]int)
	edge4 := []inter{
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 51461, read2_pos: 41561, qname: "a0"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 51561, read2_pos: 51561, qname: "a"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 51461, read2_pos: 41562, qname: "b0"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 51562, read2_pos: 51561, qname: "b"},
	}
	interchromfilter(counts, edge4)
	c.Assert(counts["a0"], Equals, 0)
	c.Assert(counts["b0"], Equals, 0)
	c.Assert(counts["a"], Equals, 0)
	c.Assert(counts["b"], Equals, 0)
}

func (s *InterChromosomalTest) TestInterChromFilterEdge5(c *C) {
	// test a run of in-window reads
	counts := make(map[string]int)
	edge5 := []inter{
		inter{read1_tid: 1, read2_tid: 12, read1_pos: 51461, read2_pos: 41561, qname: "a0"},
		inter{read1_tid: 1, read2_tid: 2, read1_pos: 51561, read2_pos: 51561, qname: "a"},
		inter{read1_tid: 1, read2_tid: 13, read1_pos: 51461, read2_pos: 41562, qname: "b0"},
		inter{read1_tid: 1, read2_tid: 22, read1_pos: 51562, read2_pos: 51561, qname: "b"},
	}
	interchromfilter(counts, edge5)
	c.Assert(counts["a0"], Equals, -1)
	c.Assert(counts["b0"], Equals, -1)
	c.Assert(counts["a"], Equals, -1)
	c.Assert(counts["b"], Equals, -1)
}
