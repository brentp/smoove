package annotate

import (
	"bytes"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"

	arg "github.com/alexflint/go-arg"
	"github.com/biogo/store/interval"
	"github.com/brentp/go-athenaeum/unsplit"
	"github.com/brentp/vcfgo"
	"github.com/brentp/xopen"
)

type cliargs struct {
	GFF string `arg:"-g,help:path to GFF for gene annotation"`
	VCF string `arg:"positional,required,help:path to VCF(s) to annotate."`
}

func (c cliargs) Description() string {
	return `
GFF3 annotation files can be downloaded from Ensembl:
ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/
ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/ `
}

func id2name(toks [][]byte) (id, name string) {
	info := unsplit.New(toks[8], []byte{';'})
	for id == "" || name == "" {
		p := info.Next()
		if p == nil {
			break
		}
		if bytes.HasPrefix(p, []byte("ID=")) {
			tmp := bytes.SplitN(p, []byte{':'}, 2)
			id = string(tmp[len(tmp)-1])
			continue
		}
		if bytes.HasPrefix(p, []byte("Name=")) {
			name = string(p[5:])
		}

	}
	return id, name
}

// Integer-specific intervals
type irange struct {
	Start, End int
	UID        uintptr
	Ftype      string
	Name       string
}

func (i irange) Overlap(b interval.IntRange) bool {
	// Half-open interval indexing.
	return i.End > b.Start && i.Start < b.End
}
func (i irange) ID() uintptr              { return i.UID }
func (i irange) Range() interval.IntRange { return interval.IntRange{i.Start, i.End} }

// Overlaps checks for overlaps and fills result.
func Overlaps(trees map[string]*interval.IntTree, chrom string, start, end int, result *[]irange) {
	var tree = trees[chrom]
	if len(*result) != 0 {
		*result = (*result)[:0]
	}
	if tree == nil {
		if strings.HasPrefix(chrom, "chr") {
			tree = trees[chrom[3:len(chrom)]]
		}
		if tree == nil {
			return
		}
	}

	q := irange{Start: start, End: end, UID: uintptr(tree.Len() + 1)}

	tree.DoMatching(func(iv interval.IntInterface) bool {
		*result = append(*result, iv.(irange))
		return false
	}, q)

}

func chromStartEnd(toks [][]byte) (string, int, int) {
	s, err := strconv.Atoi(string(toks[3]))
	if err != nil {
		panic(err)
	}
	e, err := strconv.Atoi(string(toks[4]))
	if err != nil {
		panic(err)
	}
	return string(toks[0]), s, e
}

func getParent(info []byte) string {
	var s, e int
	if bytes.HasPrefix(info, []byte("Parent=")) {
		s = 7
		e = bytes.IndexRune(info, ';')
		if e == -1 {
			e = len(info)
		}
	} else {
		s = bytes.Index(info, []byte("Parent="))
		if s == -1 {
			panic(fmt.Sprintf("no parent in %s", info))
		}
		s += 7
		e = bytes.IndexRune(info[s:], ';')
		if e == -1 {
			e = len(info)
		} else {
			e += s
		}
	}
	if e > len(info) || s > len(info) || s < 0 || e < 0 {
		log.Println(s, e, string(info))
		log.Fatal("bad")
	}
	name := string(info[s:e])
	tmp := strings.SplitN(name, ":", 2)
	return tmp[len(tmp)-1]
}

func anyEqual(q []byte, db []string) bool {
	for _, d := range db {
		if string(q) == d {
			return true
		}
	}
	return false
}

func passOne(path string, ftypes []string) map[string]string {
	geneIdToNameMap := make(map[string]string, 64)
	f, err := xopen.Ropen(path)
	if err != nil {
		panic(err)
	}

	for {
		var id, name string
		line, err := f.ReadBytes('\n')
		if len(line) != 0 {
			if line[0] == '#' {
				continue
			}
			toks := bytes.SplitN(line, []byte{'\t'}, 11)
			if !anyEqual(toks[2], ftypes) {
				continue
			}

			id, name = id2name(toks)
			for _, ftype := range ftypes {
				if !strings.HasSuffix(ftype, "gene") {
					name = getParent(toks[8])
				}
				if id != "" && name != "" {
					geneIdToNameMap[id] = name
				}
			}
		}
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
	}
	return geneIdToNameMap
}

const upStreamDist = 5000
const downStreamDist = 5000

func readGff(path string) map[string]*interval.IntTree {
	t := make(map[string]*interval.IntTree, 20)
	geneToNameMap := passOne(path, []string{"gene"})
	transcriptToGeneMap := passOne(path, []string{"transcript", "mRNA"})
	if len(geneToNameMap) == 0 {
		log.Fatal("no records found with 'gene' type in gff")
	}
	if len(transcriptToGeneMap) == 0 {
		log.Fatal("no records found with 'transcript' or 'mRNA' type in gff")
	}

	f, err := xopen.Ropen(path)
	if err != nil {
		panic(err)
	}
	var k int
	for {
		var id, name string
		line, err := f.ReadBytes('\n')
		if len(line) != 0 {
			if line[0] == '#' {
				continue
			}
			line = bytes.TrimSpace(line)
			toks := bytes.SplitN(line, []byte{'\t'}, 11)
			if bytes.Equal(toks[2], []byte("gene")) {
				chrom, start, end := chromStartEnd(toks)
				id, name = id2name(toks)
				if _, ok := t[chrom]; !ok {
					t[chrom] = &interval.IntTree{}
				}
				if err := t[chrom].Insert(irange{Start: start, End: end, UID: uintptr(k), Ftype: "gene", Name: name}, false); err != nil {
					panic(err)
				}
				k++
				var dstart, dend int
				if toks[6][0] == '-' {
					start, end = end, end+upStreamDist
					dstart, dend = start-downStreamDist, start
					if dstart < 0 {
						dstart = 0
					}

				} else if toks[6][0] == '+' {
					if start > upStreamDist {
						start, end = start-upStreamDist, start
					} else {
						start, end = 0, start
					}
					dstart, dend = end, end+upStreamDist
				} else {
					log.Println("invalid strand: ", toks[6])
				}
				if err := t[chrom].Insert(irange{Start: start, End: end, UID: uintptr(k), Ftype: "upstream", Name: name}, false); err != nil {
					panic(err)
				}
				k++
				if err := t[chrom].Insert(irange{Start: dstart, End: dend, UID: uintptr(k), Ftype: "downstream", Name: name}, false); err != nil {
					panic(err)
				}
				k++
				continue
			}
			if bytes.Equal(toks[2], []byte("transcript")) {
				continue
			}

			if !(bytes.Equal(toks[2], []byte("exon")) || bytes.HasSuffix(toks[2], []byte("prime_UTR"))) {
				continue
			}
			// exon looks up transcript, but gene name is only in gene so we map transcript to gene and then store the gene name
			id = getParent(toks[8])
			//log.Println("id:", id)
			//log.Println("transcript map:", transcriptToGeneMap)
			//log.Println("gene map:", geneToNameMap)
			name = transcriptToGeneMap[id]
			if t, ok := geneToNameMap[name]; ok {
				name = t
			}
			if name == "" {
				continue // psuedo gene
			}
			chrom, start, end := chromStartEnd(toks)

			if _, ok := t[chrom]; !ok {
				t[chrom] = &interval.IntTree{}
			}
			if err := t[chrom].Insert(irange{Start: start, End: end, UID: uintptr(k), Ftype: string(toks[2]), Name: name}, false); err != nil {
				panic(err)
			}
			k++
		}

		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
	}
	//log.Fatal(geneToNameMap["ENST00000565691"])
	return t
}

func ostring(o irange, start, end int) string {
	ostr := o.Name + "|"
	ostr += o.Ftype
	return ostr
}

type counter struct {
	n     int
	bases int
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func getval(fields map[string]string, key string, vdefault float64) (float64, error) {
	sval, ok := fields[key]
	var val float64
	var err error
	if !ok {
		return vdefault, fmt.Errorf("field %s not found", key)
	}
	val, err = strconv.ParseFloat(sval, 64)
	if err != nil {
		return vdefault, fmt.Errorf("couldn't parse float from %s for field %s", sval, key)
	}
	return val, nil
}

func abs(a int) int {
	if a < 0 {
		return -a
	}
	return a
}

func getcisum(v *vcfgo.Variant) float64 {
	cipos, err := v.Info_.Get("CIPOS")
	if err != nil {
		log.Fatal(err)
	}
	ciend, err := v.Info_.Get("CIEND")
	_ = ciend
	if err != nil {
		log.Fatal(err)
	}
	cisum := abs(cipos.([]int)[0])
	cisum += abs(cipos.([]int)[1])
	cisum += abs(ciend.([]int)[0])
	cisum += abs(ciend.([]int)[1])
	return float64(cisum)
}

func smooveSampleQuality(fields map[string]string, cisum float64) int {
	ab, err := getval(fields, "AB", 0.5)
	if err != nil {
		return 0
	}
	as, err := getval(fields, "AS", 0)
	if err != nil {
		return 0
	}
	// high quality
	if ab > 0.167 {
		if as > 1.5 && cisum < 40 {
			return 4
		}
	}
	asc, err := getval(fields, "ASC", 0)
	if err != nil {
		return 0
	}
	// low quality
	if ab <= 0.167 && asc < 0.5 && as < 0.5 {
		return 1
	}
	ap, _ := getval(fields, "AP", 0)
	if ap > 7 && cisum < 200 {
		return 3
	}
	return 4
}

func setSmooveQuality(variant *vcfgo.Variant) float64 {
	var n, sum float64
	cisum := getcisum(variant)
	variant.Format = append(variant.Format, "SHQ")
	for _, s := range variant.Samples {
		if len(s.GT) == 2 && s.GT[0] == 0 && s.GT[1] == 1 {
			q := smooveSampleQuality(s.Fields, cisum)
			s.Fields["SHQ"] = strconv.Itoa(q)
			n++
			sum += float64(q)

		} else {
			s.Fields["SHQ"] = "-1"
		}
	}
	if n == 0 {
		return -1
	}
	return sum / n
}

func annotate(vcf *vcfgo.Reader, out *vcfgo.Writer, genes map[string]*interval.IntTree) {

	overlapping := make([]irange, 0, 24)

	for {
		variant := vcf.Read()
		if variant == nil {
			break
		}
		mq := setSmooveQuality(variant)
		variant.Info().Set("MSHQ", mq)
		Overlaps(genes, variant.Chromosome, int(variant.Start()), int(variant.End()), &overlapping)
		if len(overlapping) == 0 {
			out.WriteVariant(variant)
			continue
		}
		m := make(map[string]counter)
		for _, o := range overlapping {
			var c counter
			var ok bool
			var key = o.Name + "|" + o.Ftype
			if c, ok = m[key]; !ok {
				c = counter{}
			}
			c.n += 1
			c.bases += min(int(variant.End()), o.End) - max(int(variant.Start()), o.Start)
			//log.Printf("%s: %+v", key, c)
			m[key] = c
		}
		sg := ""
		for name, counts := range m {
			sg += fmt.Sprintf(",%s:%d:%d", name, counts.n, counts.bases)
		}

		variant.Info().Set("smoove_gene", sg[1:])
		out.WriteVariant(variant)
	}

	if err := vcf.Error(); err != nil {
		log.Println(err)
	}

}

func Main() {

	cli := &cliargs{}
	arg.MustParse(cli)

	genes := readGff(cli.GFF)

	f, err := xopen.Ropen(cli.VCF)
	if err != nil {
		panic(err)
	}
	vcf, err := vcfgo.NewReader(f, false)
	vcf.AddFormatToHeader("SHQ", "1", "Integer", "smoove het quality: -1==NOT HET 0==UNKNOWN, 1==VERYLOW, 3=MED, 4=HIGH")
	vcf.AddInfoToHeader("MSHQ", "1", "Float", "mean smoove het quality: -1==NOT HET 0==UNKNOWN, 1==VERYLOW, 3=MED, 4=HIGH")
	vcf.AddInfoToHeader("smoove_gene", ".", "String", "genes overlapping variants. format is gene|feature:nfeatures:nbases,...")
	if err != nil {
		panic(err)
	}

	out, err := vcfgo.NewWriter(os.Stdout, vcf.Header)
	if err != nil {
		panic(err)
	}

	annotate(vcf, out, genes)
}
