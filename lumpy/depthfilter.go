package lumpy

import (
	"bytes"
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

func nm(r *sam.Record) int {
	if nm, ok := r.Tag([]byte{'N', 'M'}); ok {
		nmv := nm.Value()
		if v, ok := nmv.(uint32); ok {
			return int(uint32(v))
		} else if v, ok := nmv.(int32); ok {
			return int(int32(v))
		} else if v, ok := nmv.(uint8); ok {
			return int(uint8(v))
		} else if v, ok := nmv.(int8); ok {
			return int(int8(v))
		} else if v, ok := nmv.(int16); ok {
			return int(int16(v))
		} else if v, ok := nmv.(uint16); ok {
			return int(uint16(v))
		} else {
			shared.Slogger.Printf("bad NM type: %s was: %s", nmv, nm.Kind())
		}
	}
	return 0

}
func nm_above(r *sam.Record, max_mismatches int) bool {

	nmc := nm(r)
	if nmc <= max_mismatches {
		return false
	}
	nEvents := 0
	for _, cig := range r.Cigar {
		if cig.Type() == sam.CigarInsertion || cig.Type() == sam.CigarDeletion {
			nEvents += 1
			nmc -= (cig.Len() - 1)
		}
	}
	return nmc > max_mismatches || nEvents > 2
}

func interOrDistant(r *sam.Record) bool {
	return (r.Ref.ID() != r.MateRef.ID()) || (abs(r.Pos-r.MatePos) > 8000000)
}

// check if the alignment has a splitter nearby.
// useful for checking when small tandem dups of more than 3 copies
// are encompassed in a read.
func nearbySplitter(r *sam.Record) bool {
	tags, ok := r.Tag([]byte{'S', 'A'})
	if !ok {
		return false
	}
	atags := bytes.Split(tags[3:len(tags)-1], []byte{';'})
	for _, t := range atags {
		pieces := bytes.Split(t, []byte{','})
		if string(pieces[0]) != r.MateRef.Name() && string(pieces[0]) != r.Ref.Name() {
			continue
		}
		pos, err := strconv.Atoi(string(pieces[1]))
		if err != nil {
			continue
		}
		if abs(pos-r.MatePos) < 500 {
			return true
		}
		if abs(pos-r.Start()) < 500 {
			return true
		}
		if abs(pos-r.End()) < 500 {
			return true
		}
	}
	return false
}

// 1. any read with > 6 mismatches is removed.
// 2. any read with a splitter that goes to within 500 bases of itself or its mate is kept.
// 3. any interchromosomal with > 4 mismatches is discarded.
// 2. an interchromosomal where > 40% of the read is soft-clipped is removed
// NOTE: "interchromosomal" here includes same chrom with end > 8MB away.
func sketchyInterchromosomalOrSplit(r *sam.Record) bool {
	if nm_above(r, 6) {
		return true
	}
	if nearbySplitter(r) {
		return false
	}
	if interOrDistant(r) {
		if nm_above(r, 4) {
			return true
		}
		// skip inter-chrom with >XX% soft if no SA (checked above in nearby splitter)
		if s := softMatchCount(r); s.pSkip() > 0.40 {
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
	lowMQ := 0
	badInter := 0
	for {
		rec, err := br.Read()
		if rec != nil {
			tot += 1
			rchrom := rec.Ref.Name()
			if rec.MapQ < MinMapQuality {
				lowMQ++
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
			if !extraFilters {
				check(bw.Write(rec))
				continue
			}
			if sketchyInterchromosomalOrSplit(rec) {
				badInter++
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
	shared.Slogger.Printf("removed %d alignments out of %d (%.2f%%) with low mapq, depth > %d, or from excluded chroms from %s in %.0f seconds\n",
		removed, tot, pct, maxdepth, filepath.Base(fbam), time.Now().Sub(t0).Seconds())
	//shared.Slogger.Printf("of those, %d were removed due to low mapping quality\n", lowMQ)

	pct = float64(badInter) / float64(tot) * 100
	shared.Slogger.Printf("removed %d alignments out of %d (%.2f%%) that were bad interchromosomals or flanked-splitters from %s\n",
		badInter, tot, pct, filepath.Base(fbam))

	singletonfilter(fbam, strings.HasSuffix(fbam, ".split.bam"), tot)

}

func remove_sketchy_all(bams []filter, maxdepth int, fasta string, fexclude string, filter_chroms []string, extraFilters bool) {

	if _, err := exec.LookPath("mosdepth"); err != nil {
		shared.Slogger.Print("mosdepth executable not found, proceeding without removing high-coverage regions.")
	}

	pch := make(chan string, runtime.GOMAXPROCS(0))
	var wg sync.WaitGroup

	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		wg.Add(1)
		go func() {
			for bamp := range pch {
				remove_sketchy(bamp, maxdepth, fasta, fexclude, filter_chroms, extraFilters)
				proc := exec.Command("samtools", "index", bamp)
				proc.Stderr = os.Stderr
				proc.Stdout = os.Stdout
				check(proc.Run())
			}
			wg.Done()
		}()
	}

	for _, b := range bams {
		pch <- b.disc
	}
	for _, b := range bams {
		pch <- b.split
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

func singletonfilter(fbam string, split bool, originalCount int) {

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
	nwritten := 0
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
			nwritten++
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
	shared.Slogger.Printf("removed %d singletons of %d reads (%.2f%%) from %s in %.0f seconds", removed, tot, pct, filepath.Base(f.Name()), time.Now().Sub(t0).Seconds())
	pct = 100 * float64(nwritten) / float64(originalCount)
	shared.Slogger.Printf("%d reads (%.2f%%) of the original %d remain from %s", nwritten, pct, originalCount, filepath.Base(f.Name()))
}
