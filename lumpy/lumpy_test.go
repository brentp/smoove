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
