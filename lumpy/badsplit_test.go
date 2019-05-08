package lumpy

import (
	"github.com/biogo/hts/sam"
	. "gopkg.in/check.v1"
)

type BadSplitTest struct{}

var _ = Suite(&BadSplitTest{})

func (s *BadSplitTest) TestCounts(c *C) {

	a, err := sam.ParseCigar([]byte("100M50S"))
	c.Assert(err, IsNil)
	b, err := sam.ParseCigar([]byte("100S50M"))
	c.Assert(err, IsNil)

	cnts := countBases([]sam.Cigar{a, b}, 150, nil)
	c.Assert(len(cnts), Equals, 150)
	for i := 0; i < 150; i++ {
		c.Assert(cnts[i], Equals, int8(1))
	}

	c.Assert(isBad(cnts), Equals, false)
}

func (s *BadSplitTest) TestNotMatchingCounts(c *C) {

	a, err := sam.ParseCigar([]byte("100M10H"))
	c.Assert(err, IsNil)
	b, err := sam.ParseCigar([]byte("10M100S"))
	c.Assert(err, IsNil)

	cnts := countBases([]sam.Cigar{a, b}, 150, nil)
	for i := 0; i < 10; i++ {
		c.Assert(cnts[i], Equals, int8(2))
	}
	for i := 10; i < 100; i++ {
		c.Assert(cnts[i], Equals, int8(1))
	}
	for i := 100; i < 110; i++ {
		c.Assert(cnts[i], Equals, int8(0))
	}
	c.Assert(isBad(cnts), Equals, true)
}
