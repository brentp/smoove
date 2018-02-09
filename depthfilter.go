package main

import (
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"time"

	"github.com/biogo/hts/bam"
	"github.com/brentp/goleft/depth"
	"github.com/valyala/fasttemplate"
)

const MinMapQuality = byte(20)

// run mosdepth to find high coverage regions
// read the bed file into an interval tree, iterate over the file,
// and only output reads that do not overlap high coverage intervals.
func remove_high_depth(fbam string, maxdepth int, fasta string, fexclude string, filter_chroms []string) {
	t0 := time.Now()

	f, err := ioutil.TempFile("", "lumpy-smoother-mosdepth-")
	check(err)
	defer f.Close()
	defer os.Remove(f.Name())
	cmd := `
export MOSDEPTH_Q0=OK
export MOSDEPTH Q1=HIGH
set -euo pipefail
samtools index {{bam}}
mosdepth -f {{fasta}} -n --quantize {{md1}}: {{prefix}} {{bam}}
rm {{prefix}}.mosdepth.dist.txt
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

	t := depth.ReadTree(f.Name()+".quantized.bed.gz", fexclude)

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
			if rec.MapQ < MinMapQuality {
				removed++
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
			// if we made it here, we know the chrom is OK.
			// so check if mate is from a different chromosome and exclude if mate from filtered chroms
			if rec.MateRef.ID() != rec.Ref.ID() && contains(filter_chroms, rec.MateRef.Name()) {
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
	logger.Printf("removed %d alignments out of %d (%.2f%%) with depth > %d or from excluded chroms from %s in %.0f seconds\n",
		removed, tot, pct, maxdepth, filepath.Base(fbam), time.Now().Sub(t0).Seconds())
	if !strings.HasSuffix(fbam, ".split.bam") {
		singletonfilter(fbam)
	}

}

func contains(haystack []string, needle string) bool {
	for _, h := range haystack {
		if h[0] != '~' && h == needle {
			return true
		}
		if h[0] == '~' {
			if match, _ := regexp.MatchString(h[1:], needle); match {
				return true
			}
		}
	}
	return false
}

func remove_high_depths(bams []filtered, maxdepth int, fasta string, fexclude string, filter_chroms []string) {

	if _, err := exec.LookPath("mosdepth"); err != nil {
		logger.Print("mosdepth executable not found, proceeding without removing high-coverage regions.")
		return
	}

	pch := make(chan []string, runtime.GOMAXPROCS(0))
	var wg sync.WaitGroup

	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		wg.Add(1)
		go func() {
			for pair := range pch {
				remove_high_depth(pair[0], maxdepth, fasta, fexclude, filter_chroms)
				remove_high_depth(pair[1], maxdepth, fasta, fexclude, filter_chroms)
				proc := exec.Command("bash", "-c", fmt.Sprintf("samtools index %s && samtools index %s", pair[0], pair[1]))
				proc.Stderr = os.Stderr
				proc.Stdout = os.Stdout
				check(proc.Run())
				/*
					newbams := evidencewindow.EvidenceRemoval(&evidencewindow.Args{Window: 1500, Evidence: 2, Bams: pair})
					check(os.Rename(newbams[0], pair[0]))
					check(os.Rename(newbams[1], pair[1]))
					os.Remove(pair[0] + ".bai")
					os.Remove(pair[1] + ".bai")
				*/
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

func singletonfilter(fbam string) {

	f, err := os.Open(fbam)
	check(err)

	br, err := bam.NewReader(f, 1)
	check(err)

	counts := make(map[string]int, 100)

	for {
		rec, err := br.Read()
		if rec != nil {
			counts[rec.Name]++
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
			if counts[rec.Name] == 1 {
				if _, ok := rec.Tag([]byte{'S', 'A'}); !ok {
					removed++
					continue
				}
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
	logger.Printf("removed %d singletons out of %d (%.2f%%), non-splitter reads from %s", removed, tot, pct, filepath.Base(f.Name()))
}
