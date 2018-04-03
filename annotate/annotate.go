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

// Overlaps checks for overlaps without pulling intervals from the tree.
func Overlaps(tree *interval.IntTree, start, end int, result *[]irange) {
	if len(*result) != 0 {
		*result = (*result)[:0]
	}
	if tree == nil {
		return
	}

	q := irange{Start: start, End: end, UID: uintptr(tree.Len())}

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

func passOne(path string, ftype string) map[string]string {
	geneIdToNameMap := make(map[string]string, 64)
	f, err := xopen.Ropen(path)
	if err != nil {
		panic(err)
	}
	btype := []byte(ftype)

	for {
		var id, name string
		line, err := f.ReadBytes('\n')
		if len(line) != 0 {
			if line[0] == '#' {
				continue
			}
			toks := bytes.SplitN(line, []byte{'\t'}, 11)
			if !bytes.Equal(toks[2], btype) {
				continue
			}

			id, name = id2name(toks)
			if !strings.HasSuffix(ftype, "gene") {
				name = getParent(toks[8])
			}
			if id != "" && name != "" {
				geneIdToNameMap[id] = name
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

const upstreamDist = 1000

func readGff(path string) map[string]*interval.IntTree {
	t := make(map[string]*interval.IntTree, 20)
	geneToNameMap := passOne(path, "gene")
	transcriptToGeneMap := passOne(path, "transcript")
	/*
			x, err := os.Create("x.debug")
			for k, v := range transcriptToGeneMap {
				fmt.Fprintf(x, "%s\t%s\n", k, v)
			}
		x.Close()
	*/
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
				if toks[6][0] == '-' {
					start, end = end, end+upstreamDist

				} else if toks[6][0] == '+' {
					start, end = start-upstreamDist, start
				} else {
					log.Println("invalid strand: ", toks[6])
				}
				if err := t[chrom].Insert(irange{Start: start, End: end, UID: uintptr(k), Ftype: "upstream", Name: name}, false); err != nil {
					panic(err)
				}

				k += 1
			}
			if bytes.Equal(toks[2], []byte("transcript")) {
				continue
			}

			if !(bytes.Equal(toks[2], []byte("exon")) || bytes.HasSuffix(toks[2], []byte("prime_UTR"))) {
				continue
			}
			// exon looks up transcript, but gene name is only in gene so we map transcript to gene and then store the gene name
			id = getParent(toks[8])
			name = transcriptToGeneMap[id]
			if t, ok := geneToNameMap[name]; ok {
				name = t
			}
			//log.Println("name:", name)
			//log.Println("gene:", geneToNameMap[name])
			//log.Println(string(toks[2])+":", name, " from:", id, " info:", string(toks[8]))
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
			k += 1
		}
		/*
		   1       ensembl_havana  gene    65419   71585   .       +       .       ID=gene:ENSG00000186092;Name=OR4F5;biotype=protein_coding;description=olfactory receptor family 4 subfamily F member 5 [Source:HGNC Symbol%3BAcc:HGNC:14825];gene_id=ENSG00000186092;logic_name=ensembl_havana_gene;version=5
		   1       havana  mRNA    65419   71585   .       +       .       ID=transcript:ENST00000641515;Parent=gene:ENSG00000186092;Name=OR4F5-202;biotype=protein_coding;ccdsid=CCDS30547.1;tag=basic;transcript_id=ENST00000641515;version=1
		   1       havana  exon    65419   65433   .       +       .       Parent=transcript:ENST00000641515;Name=ENSE00003812156;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003812156;rank=1;version=1
		   1       havana  five_prime_UTR  65419   65433   .       +       .       Parent=transcript:ENST00000641515
		   1       havana  exon    65520   65573   .       +       .       Parent=transcript:ENST00000641515;Name=ENSE00003813641;constitutive=0;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=ENSE00003813641;rank=2;version=1
		   1       havana  five_prime_UTR  65520   65573   .       +       .       Parent=transcript:ENST00000641515
		   1       havana  five_prime_UTR  69037   69090   .       +       .       Parent=transcript:ENST00000641515

		*/

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

func annotate(vcf *vcfgo.Reader, out *vcfgo.Writer, genes map[string]*interval.IntTree) {

	overlapping := make([]irange, 0, 24)

	for {
		variant := vcf.Read()
		if variant == nil {
			break
		}
		Overlaps(genes[variant.Chromosome], int(variant.Start()), int(variant.End()), &overlapping)
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
	vcf, err := vcfgo.NewReader(f, true)
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
