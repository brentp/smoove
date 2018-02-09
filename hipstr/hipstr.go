// Package hipstr parallelizes calling of HipSTR and does best-practices filtering.
package hipstr

import (
	"bufio"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"os/exec"
	"runtime"
	"strconv"
	"strings"
	"sync"

	arg "github.com/alexflint/go-arg"
	"github.com/brentp/xopen"
)

type hargs struct {
	Regions string   `arg:"-r,required,help:BED file of regions containing STRs"`
	Fasta   string   `arg:"-f,required,help:path to reference fasta file"`
	Bams    []string `arg:"positional,required,help:bams in which to call STRs"`
}

func Main() {

	cli := &hargs{}
	arg.MustParse(cli)

	o := bufio.NewWriter(os.Stdout)
	defer o.Flush()

	HipStr(cli, o)
}

func check(e error) {
	if e != nil {
		panic(e)
	}
}

// TODO: use --bams-files <FILE> so that command-line doesnt get so large.
const cmd_tmpl = `HipSTR --silent --bams %s --fasta %s --regions %s --str-vcf %s`

func worker(args hargs, delete bool) (*xopen.Reader, string) {
	if delete {
		defer os.Remove(args.Regions)
	}
	t, err := ioutil.TempFile("", "hipstr-vcf-gz")
	check(err)
	check(t.Close())
	defer os.Remove(t.Name())
	tname := t.Name() + ".vcf.gz"
	cmd := fmt.Sprintf(cmd_tmpl, strings.Join(args.Bams, ","), args.Fasta, args.Regions, tname)
	p := exec.Command("bash", "-c", cmd)

	p.Stderr = os.Stderr
	p.Stdout = os.Stdout

	check(p.Run())

	x, err := xopen.Ropen(tname)
	check(err)
	return x, tname
}

const variantsPerRegion = 400

func makeRegions(regions string) chan string {
	f, err := xopen.Ropen(regions)
	check(err)
	ch := make(chan string, 10)
	go func() {
		defer f.Close()
		k, group := 0, 0
		w, err := xopen.Wopen("tmp:hipstr-" + strconv.Itoa(group) + "-")
		check(err)
		for {
			line, err := f.ReadString('\n')
			if len(line) != 0 {
				k++
				_, verr := w.WriteString(line)
				check(verr)
				if k == variantsPerRegion {
					k = 0
					w.Close()
					log.Println(w.Name())
					ch <- w.Name()
					group++
					w, err = xopen.Wopen("tmp:hipstr-" + strconv.Itoa(group) + "-")
					check(err)
				}
			}
			if err == io.EOF {
				break
			}
			check(err)
		}
		if k > 0 {
			w.Close()
			ch <- w.Name()
		}
		close(ch)
	}()
	return ch
}

type wl struct {
	mu *sync.Mutex
	io.Writer
	i int
}

func HipStr(args *hargs, wtr io.Writer) {
	if _, err := exec.LookPath("HipSTR"); err != nil {
		panic("error can't find hipstr executable")
	}
	var wg sync.WaitGroup

	out := wl{mu: &sync.Mutex{}, Writer: wtr, i: 0}

	ch := makeRegions(args.Regions)

	for i := 0; i < runtime.GOMAXPROCS(0); i++ {
		wg.Add(1)
		go func() {
			for r := range ch {
				gz, path := worker(hargs{Fasta: args.Fasta, Regions: r, Bams: args.Bams}, true)
				out.mu.Lock()
				if out.i == 0 {
					_, err := io.Copy(out, gz)
					check(err)
					out.i++
					out.mu.Unlock()
					gz.Close()
					os.Remove(path)
					continue
				}
				for {
					// only write lines that don't start with '#'
					line, err := gz.ReadString('\n')
					if len(line) != 0 {
						if line[0] != '#' {
							out.Write([]byte(line))
							_, err := io.Copy(out, gz)
							check(err)
							break
						}

					}
					if err == io.EOF {
						break
					}
					check(err)
				}

				out.i++
				log.Println("at:", out.i)
				out.mu.Unlock()
				gz.Close()
				os.Remove(path)
			}
			wg.Done()
		}()
	}
	wg.Wait()
}
