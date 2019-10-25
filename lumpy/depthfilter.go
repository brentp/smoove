package lumpy

import (
	"bytes"
	"io"
	"io/ioutil"
	"math"
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
	"github.com/kyroy/kdtree"
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
			shared.Slogger.Printf("bad NM type: %s was: %s", nmv, string(nm.Kind()))
		}
	}
	return 0

}
func nm_above(r *sam.Record, max_mismatches int) bool {

	nmc := nm(r)
	if nmc <= max_mismatches {
		return false
	}
	if nmc < 10 && len(r.Cigar) == 1 && r.MapQ >= 59 && abs(r.Pos-r.MatePos) < 1000000 && r.Ref.ID() == r.MateRef.ID() {
		// with a high mapq and a mate with a good cigar, we allow up to 10 mismatches.
		if mc, ok := r.Tag([]byte{'M', 'C'}); ok {
			if len(mc) > 3 {
				_, err := strconv.Atoi(string(mc[3 : len(mc)-1]))
				if err == nil && mc[len(mc)-1] == 'M' {
					return false
				}
			}
		}
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

func interChromosomalSplit(r *sam.Record) bool {
	tags, ok := r.Tag([]byte{'S', 'A'})
	if !ok {
		return false
	}
	atags := bytes.Split(tags[3:len(tags)-1], []byte{';'})
	for _, t := range atags {
		pieces := bytes.Split(t, []byte{','})
		if string(pieces[0]) == r.MateRef.Name() || string(pieces[0]) == r.Ref.Name() {
			return false
		}
	}
	return true
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

type readCount struct {
	before int
	after  int
}

// run mosdepth to find high coverage regions
// read the bed file into an interval tree, iterate over the file,
// and only output reads that do not overlap high coverage intervals.
func remove_sketchy(fbam string, maxdepth int, fasta string, fexclude string, filter_chroms []string, extraFilters bool) readCount {
	t0 := time.Now()

	var t map[string]*interval.IntTree

	if _, err := exec.LookPath("mosdepth"); err == nil && extraFilters {

		f, err := ioutil.TempFile("", "smoove-mosdepth-")
		check(err)
		defer f.Close()
		defer os.Remove(f.Name())
		cmd := `
export MOSDEPTH_Q0=OK
export MOSDEPTH_Q1=HIGH
set -euo pipefail
samtools index {{bam}}
mosdepth -f {{fasta}} --fast-mode -n --quantize {{md1}}: {{prefix}} {{bam}}
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
	isSplit := strings.HasSuffix(fbam, ".split.bam")

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
			if isSplit && badSplitter(rec) {
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

	result := readCount{before: tot}
	result.after = singletonfilter(fbam, isSplit, tot)
	return result
}

type sampleBam struct {
	sample      string
	bam         string
	splitOrDisc string
}

func mapToCounts(sm *sync.Map) map[string][4]int {
	result := make(map[string][4]int)
	sm.Range(func(key, value interface{}) bool {
		k := key.(sampleBam)
		v := value.(readCount)
		tmp := result[k.sample]
		if k.splitOrDisc == "split" {
			tmp[0] = v.before
			tmp[2] = v.after
		} else if k.splitOrDisc == "disc" {
			tmp[1] = v.before
			tmp[3] = v.after
		} else {
			panic("unknown type:" + k.splitOrDisc)
		}
		result[k.sample] = tmp

		return true
	})

	return result
}

func remove_sketchy_all(bams []filter, maxdepth int, fasta string, fexclude string, filter_chroms []string, extraFilters bool) map[string][4]int {

	if _, err := exec.LookPath("mosdepth"); err != nil {
		shared.Slogger.Print("mosdepth executable not found, proceeding without removing high-coverage regions.")
	}

	pch := make(chan sampleBam, runtime.GOMAXPROCS(0))
	var wg sync.WaitGroup
	var sm = &sync.Map{}

	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		wg.Add(1)
		go func() {
			for bamp := range pch {
				counts := remove_sketchy(bamp.bam, maxdepth, fasta, fexclude, filter_chroms, extraFilters)
				proc := exec.Command("samtools", "index", bamp.bam)
				proc.Stderr = os.Stderr
				proc.Stdout = os.Stdout
				check(proc.Run())
				sm.Store(bamp, counts)
			}
			wg.Done()
		}()
	}

	for _, b := range bams {
		pch <- sampleBam{bam: b.disc, sample: b.sample, splitOrDisc: "disc"}
	}
	for _, b := range bams {
		pch <- sampleBam{bam: b.split, sample: b.sample, splitOrDisc: "split"}
	}
	close(pch)

	wg.Wait()
	return mapToCounts(sm)
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

// Size around an interchromosomal to look for support during filtering
const inter_window_size = 10000
const empty_tid = math.MaxUint32

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

type point struct {
	posns [2]float64
}

func (p *point) Dimensions() int {
	return 2
}

func (p *point) Dimension(i int) float64 {
	return p.posns[i]
}

func (p *point) max_distance(posns [2]float64) int {
	//
	a := abs(int(posns[0]) - int(p.posns[0]))
	b := abs(int(posns[1]) - int(p.posns[1]))
	return max(a, b)
}

var _ kdtree.Point = (*point)(nil)

func tree_key_positions(rec *sam.Record) ([2]int, [2]float64) {
	// we really have a n-chroms x n-chroms space. here, we make
	// sure we use only the lower diagonal of that space.
	// then within a single chrom, we also use only the lower diagonal.
	tids := [2]int{rec.Ref.ID(), rec.MateRef.ID()}
	posns := [2]float64{float64(rec.Pos), float64(rec.MatePos)}
	if tids[0] > tids[1] {
		var tmp = tids[0]
		tids[0] = tids[1]
		tids[1] = tmp
		var tmpf = posns[0]
		posns[0] = posns[1]
		posns[1] = tmpf
	} else if tids[0] == tids[1] && posns[0] > posns[1] {
		var tmpf = posns[0]
		posns[0] = posns[1]
		posns[1] = tmpf
	}
	return tids, posns
}

func drop_orphans(br *bam.Reader, inters map[[2]int]*kdtree.KDTree, counts map[string]int, split bool) int {
	n_dropped := 0
	n_kept := 0
	for {
		rec, err := br.Read()
		if rec != nil {
			name := rec.Name
			if split {
				name = "A" + name[1:]
			}
			if counts[name] < 2 {
				continue
			}
			if !split && !(interOrDistant(rec) && rec.MateRef.ID() != -1 && (rec.Flags&sam.Read1 != 0)) {
				continue
			}

			if split && !interChromosomalSplit(rec) {
				continue
			}

			key, posns := tree_key_positions(rec)
			t, ok := inters[key]
			if !ok {
				continue
			}
			found := t.KNN(&point{posns: posns}, 2)
			distant := false
			self := false
			for _, f := range found {
				d := f.(*point).max_distance(posns)
				if d > 1000 {
					// should have only self and 1 other point so
					// fi we have 1 point that's far, we can skip
					distant = true
				} else if d == 0 {
					self = true
				}
			}
			// we test some reads that are not interchromosomal, so we require finding self
			if distant && self {
				counts[name]--
				n_dropped++
			} else {
				n_kept++
			}

		}
		if err == io.EOF {
			break
		}
		check(err)
	}

	shared.Slogger.Printf("kept %d putative orphans", n_kept)
	return n_dropped
}
func singletonfilter(fbam string, split bool, originalCount int) int {

	t0 := time.Now()

	// keyed by chroms
	inters := make(map[[2]int]*kdtree.KDTree, 20)

	f, err := os.Open(fbam)
	check(err)

	br, err := bam.NewReader(f, 1)
	check(err)

	counts := make(map[string]int, 100)

	// first pass finds pairs (so we can later exclude orphans and fills the kdtree
	// for interchromsomals
	for {
		rec, err := br.Read()
		if rec != nil {
			name := rec.Name
			if split {
				// split sets first letter to A or B. Here we normalize.
				name = "A" + name[1:]
				if interChromosomalSplit(rec) {
					key, posns := tree_key_positions(rec)
					t, ok := inters[key]
					if !ok {
						t = kdtree.New([]kdtree.Point{})
						inters[key] = t
					}
					// if the chroms were flipped, we have to flip the points as well
					t.Insert(&point{posns: posns})
				}
			} else {
				// Use tree to later eliminate interchromosomal reads are orphaned. i.e. there
				// aren't other nearby reads that could be signaling an SV.
				// Note that this only applies to discordants, not splitters.
				// There are sometimes reads with unmapped mates that get called as discordants
				// so we have to check tid of mate != -1
				if interOrDistant(rec) && rec.MateRef.ID() != -1 && (rec.Flags&sam.Read1 != 0) {
					key, posns := tree_key_positions(rec)
					t, ok := inters[key]
					if !ok {
						t = kdtree.New([]kdtree.Point{})
						inters[key] = t
					}
					// if the chroms were flipped, we have to flip the points as well
					t.Insert(&point{posns: posns})
				}
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
	td := time.Now()
	ndropped := drop_orphans(br, inters, counts, split)
	label := "discordant"
	if split {
		label = "split"
	}
	shared.Slogger.Printf("removed %d %s orphans in %.0f seconds", ndropped, label, time.Now().Sub(td).Seconds())

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
			if counts[name] < 2 {
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
	_ = f.Close()
	_ = fw.Close()

	check(os.Rename(fw.Name(), f.Name()))
	pct := 100 * float64(removed) / float64(tot)
	var additional string
	if !split {
		additional = "and isolated interchromosomals "
	}
	shared.Slogger.Printf("removed %d singletons %sof %d reads (%.2f%%) from %s in %.0f seconds", removed, additional, tot, pct, filepath.Base(f.Name()), time.Now().Sub(t0).Seconds())
	pct = 100 * float64(nwritten) / float64(originalCount)
	shared.Slogger.Printf("%d reads (%.2f%%) of the original %d remain from %s", nwritten, pct, originalCount, filepath.Base(f.Name()))
	return nwritten
}
