package lumpy

import (
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"
	"github.com/biogo/store/interval"
	"github.com/brentp/goleft/depth"
	"github.com/brentp/smoove/shared"
	"github.com/valyala/fasttemplate"
)

const MinMapQuality = byte(20)

type sm struct {
	soft  int
	hard  int
	match int
}

func (s *sm) total() float64 {
	return float64(s.match + s.hard + s.soft)
}

func (s *sm) pMatch() float64 {
	return float64(s.match) / s.total()
}

func (s *sm) pSkip() float64 {
	return float64(s.hard+s.soft) / s.total()
}

// count soft clips and matches in a read
func softMatchCount(r *sam.Record) sm {
	s := sm{}

	for _, cig := range r.Cigar {
		if cig.Type() == sam.CigarSoftClipped {
			s.soft += cig.Len()
		} else if cig.Type() == sam.CigarHardClipped {
			s.hard += cig.Len()
		} else if cig.Type() == sam.CigarMatch {
			s.match += cig.Len()
		}
	}
	return s
}

func abs(a int) int {
	if a > 0 {
		return a
	}
	return -a
}

func nm_above(r *sam.Record, max_mismatches int) bool {
	if nm, ok := r.Tag([]byte{'N', 'M'}); ok {
		if v, ok := nm.Value().(uint32); ok {
			if v > uint32(max_mismatches) {
				return true
			}
		} else if v, ok := nm.Value().(int32); ok {
			if v > int32(max_mismatches) {
				return true
			}
		}
	}
	return false

}

func interOrDistant(r *sam.Record) bool {
	return (r.Ref.ID() != r.MateRef.ID() && r.MateRef.ID() != -1) || (r.MateRef.ID() != -1 && r.Ref.ID() == r.MateRef.ID() && abs(r.Pos-r.MatePos) > 8000000)
}

// 1. interchromsomal discordant with an XA (maybe only with an XA that would be concordant)?
// 2. splitters where both end ops are soft-clips. e.g. discard 20S106M34S, but keep 87S123M
// 3. an interchromosomal where > 35% of the read is soft-clipped must have a splitter that goes to the same location as the other end.
// 4. an interchromosomal with NM tag and NM > 3 is skipped.
// 5. any read where both ends are soft clips of > 5 bases are skipped.
// 6. any read where there are more than 5 mismatches.
// NOTE: "interchromosomal" here includes same chrom with end > 8MB away.
func sketchyInterchromosomalOrSplit(r *sam.Record) bool {
	if interOrDistant(r) {
		// skip interchromosomal with XA
		if _, ok := r.Tag([]byte{'X', 'A'}); ok {
			return true
		}
		if nm_above(r, 4) {
			return true
		}

		// skip inter-chrom with >XX% soft if no SA
		if s := softMatchCount(r); s.pSkip() > 0.35 {
			tags, ok := r.Tag([]byte{'S', 'A'})
			if !ok {
				return true
			}
			hasSpl := false
			// look for splitter to matechrom
			atags := bytes.Split(tags[:len(tags)-1], []byte{';'})
			for _, t := range atags {
				pieces := bytes.Split(t, []byte{','})
				if string(pieces[0]) != r.MateRef.Name() {
					continue
				}
				pos, err := strconv.Atoi(string(pieces[1]))
				if err != nil {
					continue
				}
				if abs(pos-r.MatePos) < 1000 {
					hasSpl = true
					break
				}
			}
			if !hasSpl {
				return true
			}
		}

	}
	// if flanked by 'S'and right side is > 5 bases
	cig := r.Cigar
	if cig[0].Type() == sam.CigarSoftClipped && cig[len(cig)-1].Type() == sam.CigarSoftClipped && cig[len(cig)-1].Len() > 5 {
		return true
	}

	if nm_above(r, 5) {
		return true
	}

	// splitter flanked by 'S'
	// TODO: require 5 bases.
	if tags, ok := r.Tag([]byte{'S', 'A'}); ok {
		atags := bytes.Split(tags[:len(tags)-1], []byte{';'})
		if len(atags) > 1 {
			return false
		}
		for _, t := range atags {
			// skip things with a splitter that starts and ends with 'S'
			pieces := bytes.Split(t, []byte{','})
			cigar := string(pieces[3])
			if end := cigar[len(cigar)-1]; end != 'S' && end != 'H' {
				continue
			}
			var first rune
			for _, first = range cigar {
				if '0' <= first && first <= '9' {
					continue
				}
				break
			}
			if first != 'S' && first != 'H' {
				continue
			}
			return true
		}
	}
	return false
}

// run mosdepth to find high coverage regions
// read the bed file into an interval tree, iterate over the file,
// and only output reads that do not overlap high coverage intervals.
func remove_sketchy(fbam string, maxdepth int, fasta string, fexclude string, filter_chroms []string, extraFilters bool) {
	t0 := time.Now()

	var t map[string]*interval.IntTree

	if _, err := exec.LookPath("mosdepth"); err == nil && extraFilters {

		f, err := ioutil.TempFile("", "smoove-mosdepth-")
		check(err)
		defer f.Close()
		defer os.Remove(f.Name())
		cmd := `
export MOSDEPTH_Q0=OK
export MOSDEPTH Q1=HIGH
set -euo pipefail
samtools index {{bam}}
mosdepth -f {{fasta}} -n --quantize {{md1}}: {{prefix}} {{bam}}
rm -f {{prefix}}.mosdepth*.dist.txt
rm {{prefix}}.quantized.bed.gz.csi
`
		vars := map[string]interface{}{
			"md1":    strconv.Itoa(maxdepth + 1),
			"prefix": f.Name(),
			"bam":    fbam,
			"fasta":  fasta,
		}

		c := fasttemplate.New(cmd, "{{", "}}")
		s := c.ExecuteString(vars)
		p := exec.Command("bash", "-c", s)
		p.Stderr = os.Stderr
		p.Stdout = os.Stderr
		check(p.Run())
		defer os.Remove(f.Name() + ".quantized.bed.gz")
		defer os.Remove(fbam + ".bai")

		t = depth.ReadTree(f.Name()+".quantized.bed.gz", fexclude)
	} else {
		t = depth.ReadTree(fexclude)
	}

	fbr, err := os.Open(fbam)
	check(err)
	br, err := bam.NewReader(fbr, 1)
	check(err)

	fbw, err := ioutil.TempFile("", "smoove-mosdepth-bam")
	check(err)
	bw, err := bam.NewWriterLevel(fbw, br.Header(), 1, 1)

	// we know they are in order so avoid some lookups when filtering from remove chroms
	var last string
	var rmLast bool

	removed, tot := 0, 0
	badInter := 0
	for {
		rec, err := br.Read()
		if rec != nil {
			tot += 1
			rchrom := rec.Ref.Name()
			if rec.MapQ < MinMapQuality {
				removed++
				continue
			}
			if rec.Flags&(sam.QCFail|sam.Duplicate) != 0 {
				removed++
				continue
			}

			// remove it chrom is found and it overlaps a high-coverage region.
			if tt, ok := t[rec.Ref.Name()]; ok {
				if depth.Overlaps(tt, rec.Start(), rec.End()) {
					removed++
					continue
				}
			}
			if !extraFilters {
				check(bw.Write(rec))
				continue
			}
			if sketchyInterchromosomalOrSplit(rec) {
				badInter++
				continue
			}

			// block to check if it's in filter chroms
			// if same as last chrom ...
			if rchrom == last {
				// and we remove the last chrom
				if rmLast {
					// then skip.
					removed++
					continue
				}
			} else {
				// new chrom
				last = rchrom
				// if it's in the filtered...
				if shared.Contains(filter_chroms, last) {
					// then skip and set rmLast

					removed++
					rmLast = true
					continue
				} else {
					rmLast = false
				}
			}
			// END block to check if it's in filter chroms
			// if we made it here, we know the chrom is OK.
			// so check if mate is from a different chromosome and exclude if mate from filtered chroms
			if rec.MateRef.ID() != rec.Ref.ID() && shared.Contains(filter_chroms, rec.MateRef.Name()) {
				removed++
				continue
			}
			check(bw.Write(rec))

		}
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
	}
	check(bw.Close())
	check(br.Close())
	fbw.Close()
	fbr.Close()
	if err := os.Rename(fbw.Name(), fbam); err != nil {
		// rename won't work cross-device so we have to copy
		check(cp(fbam, fbw.Name()))
		os.Remove(fbw.Name())
	}

	pct := float64(removed) / float64(tot) * 100
	shared.Slogger.Printf("removed %d alignments out of %d (%.2f%%) with depth > %d or from excluded chroms from %s in %.0f seconds\n",
		removed, tot, pct, maxdepth, filepath.Base(fbam), time.Now().Sub(t0).Seconds())
	pct = float64(badInter) / float64(tot) * 100
	shared.Slogger.Printf("removed %d alignments out of %d (%.2f%%) that were bad interchromosomals or flanked-splitters from %s\n",
		badInter, tot, pct, filepath.Base(fbam))

	// TODO move first pass in singletonfilter that counts read names into the functon above as the reads are written.
	// this will avoid 1 pass on each bam.
	singletonfilter(fbam, strings.HasSuffix(fbam, ".split.bam"))

}

func remove_sketchy_all(bams []filter, maxdepth int, fasta string, fexclude string, filter_chroms []string, extraFilters bool) {

	if _, err := exec.LookPath("mosdepth"); err != nil {
		shared.Slogger.Print("mosdepth executable not found, proceeding without removing high-coverage regions.")
	}

	pch := make(chan []string, runtime.GOMAXPROCS(0))
	var wg sync.WaitGroup

	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		wg.Add(1)
		go func() {
			for pair := range pch {
				// do the split and disc serially to use less memory since mosdepth uses ~ 1GB per sample.
				remove_sketchy(pair[0], maxdepth, fasta, fexclude, filter_chroms, extraFilters)
				remove_sketchy(pair[1], maxdepth, fasta, fexclude, filter_chroms, extraFilters)
				proc := exec.Command("bash", "-c", fmt.Sprintf("samtools index %s && samtools index %s", pair[0], pair[1]))
				proc.Stderr = os.Stderr
				proc.Stdout = os.Stdout
				check(proc.Run())
			}
			wg.Done()
		}()
	}

	for _, b := range bams {
		pch <- []string{b.disc, b.split}
	}
	close(pch)

	wg.Wait()
}

// https://gist.github.com/elazarl/5507969#
func cp(dst, src string) error {
	s, err := os.Open(src)
	if err != nil {
		return err
	}
	// no need to check errors on read only file, we already got everything
	// we need from the filesystem, so nothing can go wrong now.
	defer s.Close()
	d, err := os.Create(dst)
	if err != nil {
		return err
	}
	if _, err := io.Copy(d, s); err != nil {
		d.Close()
		return err
	}
	return d.Close()
}

func singletonfilter(fbam string, split bool) {

	t0 := time.Now()

	f, err := os.Open(fbam)
	check(err)

	br, err := bam.NewReader(f, 1)
	check(err)

	counts := make(map[string]int, 100)

	for {
		rec, err := br.Read()
		if rec != nil {
			name := rec.Name
			if split {
				// split sets first letter to A or B. Here we normalize.
				name = "A" + name[1:]
			}
			counts[name]++
		}
		if err == io.EOF {
			break
		}
		check(err)
	}

	f.Seek(0, os.SEEK_SET)
	br, err = bam.NewReader(f, 1)
	check(err)

	fw, err := os.Create(fbam + ".tmp.bam")
	check(err)
	bw, err := bam.NewWriterLevel(fw, br.Header(), 1, 1)
	check(err)

	tot, removed := 0, 0
	for {
		rec, err := br.Read()
		// skip any singleton read as long as it's not a splitter.
		if rec != nil {
			tot += 1
			name := rec.Name
			if split {
				name = "A" + name[1:]
			}
			if counts[name] == 1 {
				removed++
				continue
			}
			check(bw.Write(rec))
		}

		if err == io.EOF {
			break
		}
		check(err)
	}
	check(bw.Close())
	check(br.Close())

	check(os.Rename(fw.Name(), f.Name()))
	pct := 100 * float64(removed) / float64(tot)
	shared.Slogger.Printf("removed %d singletons out of %d reads (%.2f%%) from %s in %.0f seconds", removed, tot, pct, filepath.Base(f.Name()), time.Now().Sub(t0).Seconds())
}
