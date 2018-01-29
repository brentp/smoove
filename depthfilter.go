package main

import (
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"sync"
	"time"

	"github.com/biogo/hts/bam"
	"github.com/brentp/goleft/depth"
	"github.com/valyala/fasttemplate"
)

// run mosdepth to find high coverage regions
// read the bed file into an interval tree, iterate over the file,
// and only output reads that do not overlap high coverage intervals.
func remove_high_depth(fbam string, maxdepth int, filter_chroms []string) {
	t0 := time.Now()

	f, err := ioutil.TempFile("", "lumpy-smoother-mosdepth-")
	check(err)
	defer f.Close()
	defer os.Remove(f.Name())
	cmd := `
export MOSDEPTH_Q0=OK
export MOSDEPTH Q1=HIGH
set -euo pipefail
if [[ ! -e {{bam}}.bai ]]; then
  samtools index {{bam}}
fi
mosdepth -n --quantize {{md1}}: {{prefix}} {{bam}}
rm {{prefix}}.mosdepth.dist.txt
rm {{prefix}}.quantized.bed.gz.csi
`
	vars := map[string]interface{}{
		"md1":    strconv.Itoa(maxdepth + 1),
		"prefix": f.Name(),
		"bam":    fbam,
	}

	c := fasttemplate.New(cmd, "{{", "}}")
	s := c.ExecuteString(vars)
	p := exec.Command("bash", "-c", s)
	p.Stderr = os.Stderr
	p.Stdout = os.Stderr
	check(p.Run())
	defer os.Remove(f.Name() + ".quantized.bed.gz")
	defer os.Remove(fbam + ".bai")

	t := depth.ReadTree(f.Name() + ".quantized.bed.gz")

	fbr, err := os.Open(fbam)
	check(err)
	br, err := bam.NewReader(fbr, 1)
	check(err)

	fbw, err := ioutil.TempFile("", "lumpy-smoother-mosdepth-bam")
	check(err)
	bw, err := bam.NewWriterLevel(fbw, br.Header(), 1, 1)

	// we know they are in order so avoid some lookups when filtering from remove chroms
	var last string
	var rmLast bool

	removed, tot := 0, 0
	for {
		rec, err := br.Read()
		if rec != nil {
			tot += 1
			rchrom := rec.Ref.Name()

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
				if contains(filter_chroms, last) {
					// then skip and set rmLast

					removed++
					rmLast = true
					continue
				} else {
					rmLast = false
				}
			}
			// END block to check if it's in filter chroms

			// remove it chrom is found and it overlaps a high-coverage region.
			if tt, ok := t[rec.Ref.Name()]; ok {
				if depth.Overlaps(tt, rec.Start(), rec.End()) {
					removed++
					continue
				}
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
	fmt.Fprintf(os.Stderr, "[lumpy-smoother] removed %d alignments out of %d (%.2f%%) with depth > %d from %s in %.0f seconds\n",
		removed, tot, pct, maxdepth, filepath.Base(fbam), time.Now().Sub(t0).Seconds())
}

func contains(haystack []string, needle string) bool {
	for _, h := range haystack {
		if h == needle {
			return true
		}
	}
	return false
}

func remove_high_depths(bams []filtered, maxdepth int, filter_chroms []string) {

	if _, err := exec.LookPath("mosdepth"); err != nil {
		fmt.Fprintln(os.Stderr, "[lumpy-smoother] mosdepth executable not found, proceeding without removing high-coverage regions.")
		return
	}

	ch := make(chan string, runtime.GOMAXPROCS(0))
	var wg sync.WaitGroup

	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		wg.Add(1)
		go func() {
			for bam := range ch {
				remove_high_depth(bam, maxdepth, filter_chroms)
			}
			wg.Done()
		}()
	}

	for _, b := range bams {
		ch <- b.disc
		ch <- b.split
	}
	close(ch)

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
