package permuvals

import (
	"math"
	"github.com/montanaflynn/stats"
	"flag"
	"sort"
	"strconv"
	"math/big"
	"math/rand"
	"strings"
	"fmt"
	"github.com/jgbaldwinbrown/fasttsv"
	"io"
	"bufio"
	"os"
	"github.com/jgbaldwinbrown/go-intervals/intervalset"
)

type OvlStats struct {
	Ncomponents int
	MeanCounts float64
	MeanCoverage float64
}

type Comparison struct {
	Permutations OverlapSets
	Overlaps Overlaps
	IterCounts OverlapCounts
	Probs Probs
}

type Bspan struct {
	Chrom string
	intervalset.Span
}

func MakeBspan(chrom string, min, max int) Bspan {
	return Bspan{chrom, intervalset.Span{Min: min, Max: max}}
}

func (b Bspan) Width() int {
	return b.Max - b.Min
}

func Covered(spans []Bspan) (out int) {
	out=0
	for _, span := range spans {
		out += span.Width()
	}
	return
}

type Flags struct {
	BedPaths string
	GenomeBedPath string
	Iterations int
	Rseed int
	Verbose bool
	MaxComps int
}

type Bed struct {
	Intervals map[string]*intervalset.Set
	Chroms []string
	Name string
}

type Beds []Bed

func MakeBed(name string) (b Bed) {
	b.Intervals = make(map[string]*intervalset.Set)
	b.Name = name
	return
}

func (b Bed) Bspans() []Bspan {
	var bspans []Bspan
	for _, chrom := range b.Chroms {
		for _, ispan := range b.Intervals[chrom].AllIntervals() {
			span := ispan.(*intervalset.Span)
			bspan := Bspan{chrom, *span}
			bspans = append(bspans, bspan)
		}
	}
	return bspans
}

func WriteBspans(w io.Writer, bs ...Bspan) {
	for _, b := range bs {
		fmt.Fprintf(w, "%v\t%v\t%v\n", b.Chrom, b.Span.Min, b.Span.Max)
	}
}

func WriteBeds(w io.Writer, bs ...Bed) {
	for _, b := range bs {
		WriteBspans(w, b.Bspans()...)
	}
}

type Overlap struct {
	Bed
	Components []string
}

type Overlaps []Overlap
type OverlapSets []Overlaps

type OverlapCount struct {
	Count []int
	Covered []int
	Name string
	Components []string
}

type OverlapCounts []OverlapCount

type Prob struct {
	Name string
	CountProb float64
	CoveredProb float64
}

type Probs []Prob

func ParseBedEntry(line []string) (s Bspan, e error) {
	if len(line) < 3 { return s, fmt.Errorf("line too short") }
	s.Chrom = line[0]
	s.Min, e = strconv.Atoi(line[1])
	if e != nil { return }
	s.Max, e = strconv.Atoi(line[2])
	return
}

func (b *Bed) AddBspans(spans ...Bspan) {
	for _, s := range spans {
		if _, ok := b.Intervals[s.Chrom]; !ok {
			b.Intervals[s.Chrom] = intervalset.EmptyV1(intervalset.MakeZeroSpan)
			b.Chroms = append(b.Chroms, s.Chrom)
		}
		b.Intervals[s.Chrom].Add(intervalset.NewSetV1([]intervalset.Interval{&(s.Span)}, intervalset.MakeZeroSpan))
	}
}

func GetBed(r io.Reader, name string) (b Bed, err error) {
	b = MakeBed(name)
	s := fasttsv.NewScanner(r)
	for s.Scan() {
		var bspan Bspan
		bspan, err = ParseBedEntry(s.Line())
		if err != nil { return }
		b.AddBspans(bspan)
	}
	return
}

func GetBedpaths(bedpaths_path string) (paths []string, err error) {
	r, err := os.Open(bedpaths_path)
	if err != nil { return }
	s := bufio.NewScanner(r)
	for s.Scan() {
		paths = append(paths, s.Text())
	}
	return
}

func GetBeds(bedpaths_path string) (beds Beds, err error) {
	bedpaths, err := GetBedpaths(bedpaths_path)
	if err != nil { return }
	for _, path := range bedpaths {
		var r *os.File
		r, err = os.Open(path)
		if err != nil { return }
		defer r.Close()

		var bed Bed
		bed, err = GetBed(r, path)
		if err != nil { return }
		beds = append(beds, bed)
	}
	return
}

func AddBed(dest *Bed, src Bed) {
	// could hand-code this to be faster
	dest.AddBspans(AllBedSpans(src)...)
}

func CombineChroms(beds ...Bed) (out []string) {
	all_chroms := make(map[string]struct{})
	for _, bed := range beds {
		for _, chrom := range bed.Chroms {
			all_chroms[chrom] = struct{}{}
		}
	}
	for chrom, _ := range all_chroms {
		out = append(out, chrom)
	}
	return out
}

func (b *Bed) IntersectBed(src Bed) {
	// fmt.Println("start b:", *b)
	// fmt.Println("src:", src)
	all_chroms := CombineChroms(*b, src)
	for _, chrom := range all_chroms {
		_, b_has_chrom := b.Intervals[chrom]
		_, src_has_chrom := src.Intervals[chrom]
		if !b_has_chrom || !src_has_chrom {
			b.Intervals[chrom] = intervalset.EmptyV1(intervalset.MakeZeroSpan)
		} else {
			b.Intervals[chrom].Intersect(src.Intervals[chrom])
		}
	}
	// fmt.Println("end b:", *b)
}

func (b *Bed) SubtractBed(src Bed) {
	all_chroms := CombineChroms(*b, src)
	for _, chrom := range all_chroms {
		_, b_has_chrom := b.Intervals[chrom]
		_, src_has_chrom := src.Intervals[chrom]
		if !b_has_chrom {
			b.Intervals[chrom] = intervalset.EmptyV1(intervalset.MakeZeroSpan)
		} else if !src_has_chrom {
		}else {
			b.Intervals[chrom].Sub(src.Intervals[chrom])
		}
	}
}

func GetOverlap(beds Beds) Overlap {
	names := []string{}
	for _, bed := range beds {
		names = append(names, bed.Name)
	}
	name := strings.Join(names, ":")
	final := Overlap{MakeBed(name), names}
	if len(beds) > 1 {
		AddBed(&final.Bed, beds[0])
	}
	for i := 1; i < len(beds); i++ {
		final.IntersectBed(beds[i])
	}
	return final
}

func BigPow(y, x *big.Int) {
	y.SetInt64(1)
	two := big.NewInt(2)
	one := big.NewInt(1)
	for i:=big.NewInt(0); i.Cmp(x) < 0; i.Add(i, one) {
		y.Mul(y, two)
	}
}

func BedPerm(beds Beds, perm *big.Int) (bedperm Beds) {
	lperm := new(big.Int).Set(perm)
	two := big.NewInt(2)
	zero := big.NewInt(0)
	modret := big.NewInt(0)
	for _, bed := range beds {
		if modret.Mod(lperm, two).Cmp(zero) != 0 {
			bedperm = append(bedperm, bed)
		}
		lperm.Div(lperm, two)
	}
	return
}

func GetAllOverlaps(beds Beds) (overlaps Overlaps) {
	perm := big.NewInt(0)
	maxperm := big.NewInt(0)
	one := big.NewInt(1)
	BigPow(maxperm, big.NewInt(int64(len(beds))))

	for ; perm.Cmp(maxperm) < 0 ; perm.Add(perm, one) {
		to_overlap := BedPerm(beds, perm)
		overlaps = append(overlaps, GetOverlap(to_overlap))
	}
	return
}

func BinSum(sum *big.Int, src *big.Int) {
	zero := big.NewInt(0)
	two := big.NewInt(2)
	mod := big.NewInt(0)
	temp := big.NewInt(0).Set(src)
	for ; temp.Cmp(zero) > 0; temp.Div(temp, two) {
		mod.Mod(temp, two)
		sum.Add(sum, mod)
	}
}

func GetLimitedOverlaps(beds Beds, maxComps int) (overlaps Overlaps) {
	perm := big.NewInt(0)
	maxperm := big.NewInt(0)
	one := big.NewInt(1)
	BigPow(maxperm, big.NewInt(int64(len(beds))))
	binsum := big.NewInt(0)

	for ; perm.Cmp(maxperm) < 0 ; perm.Add(perm, one) {
		BinSum(binsum, perm)
		if int(binsum.Int64()) <= maxComps {
			to_overlap := BedPerm(beds, perm)
			overlaps = append(overlaps, GetOverlap(to_overlap))
		}
	}
	return
}

func GetOverlaps(beds Beds, maxComps int) (overlaps Overlaps) {
	if maxComps < 0 {
		return GetAllOverlaps(beds)
	}
	return GetLimitedOverlaps(beds, maxComps)
}

func GetFlags() (f Flags) {
	flag.StringVar(&f.BedPaths, "b", "", "File containing paths to all bed files to compare")
	flag.StringVar(&f.GenomeBedPath, "g", "", "Bed file containing the lengths of all chromosomes")
	flag.IntVar(&f.Iterations, "i", -1, "Number of permutation iterations to perform")
	flag.IntVar(&f.Rseed, "r", 0, "Random seed for permutations (default 0)")
	flag.Parse()
	if f.BedPaths == "" || f.GenomeBedPath == "" {
		panic(fmt.Errorf("Missing path"))
	}
	return
}

func AllBspans(chrom string, s intervalset.SetInput) []Bspan {
	result := []Bspan{}
	s.IntervalsBetween(s.Extent(), func(x intervalset.Interval) bool {
		sresult := x.(*intervalset.Span)
		result = append(result, Bspan{Chrom: chrom, Span: *sresult})
		return true
	})
	return result
}

func AllBedSpans(b Bed) (out []Bspan) {
	for _, chrom := range b.Chroms {
		out = append(out, AllBspans(chrom, b.Intervals[chrom])...)
	}
	return
}

func FprintBeds(w io.Writer, beds ...Bed) {
	for _, b := range beds {
		bspans := AllBedSpans(b)
		for _, span := range bspans {
			fmt.Fprintf(w, "%v\t%v\t%v\t%v\t%v\t%v\t%v\n", span.Chrom, span.Min, span.Max, ".", span.Width(), Covered(bspans), b.Name)
		}
	}
}

func FprintOvlsBed(w io.Writer, os Overlaps) {
	for _, o := range os {
		bspans := AllBedSpans(o.Bed)
		for _, span := range bspans {
			fmt.Fprintf(w, "%v\t%v\t%v\t%v\t%v\t%v\t%v\n", span.Chrom, span.Min, span.Max, len(o.Components), span.Width(), Covered(bspans), o.Name)
		}
	}
}

func OvlsStatsNcomp(os Overlaps, nComponents int) OvlStats {
	var ostats OvlStats
	var counts []float64
	var covs []float64
	for _, o := range os {
		if len(o.Components) == nComponents {
			bspans := AllBedSpans(o.Bed)
			counts = append(counts, float64(len(bspans)))
			covs = append(covs, float64(Covered(bspans)))
		}
	}
	var err error
	ostats.Ncomponents = nComponents
	ostats.MeanCounts, err = stats.Mean(counts)
	if err != nil {
		ostats.MeanCounts = math.NaN()
	}
	ostats.MeanCoverage, err = stats.Mean(covs)
	if err != nil {
		ostats.MeanCoverage = math.NaN()
	}
	return ostats
}

func AllOvlsStats(os Overlaps) []OvlStats {
	max_nc := -1
	for _, o := range os {
		if len(o.Components) > max_nc {
			max_nc = len(o.Components)
		}
	}
	var out []OvlStats
	for i:=0; i<=max_nc; i++ {
		out = append(out, OvlsStatsNcomp(os, i))
	}
	return out
}

func FprintOvlsStats(w io.Writer, os ...OvlStats) {
	for _, o := range os {
		fmt.Fprintf(w, "%v\t%v\t%v\n", o.Ncomponents, o.MeanCounts, o.MeanCoverage)
	}
}

func SpanNumPositions(span Bspan, genome Bed) (positions int) {
	chrbspans := AllBedSpans(genome)
	for _, chrspan := range chrbspans {
		cw, sw := chrspan.Width(), span.Width()
		if cw >= sw {
			positions += 1 + cw - sw
		}
	}
	return
}

func Raw2Bspan(rawpos int, span Bspan, genome Bed) (newspan Bspan) {
	chrbspans := AllBedSpans(genome)
	for _, chrspan := range chrbspans {
		width := chrspan.Width() - span.Width()
		if rawpos < width {
			newspan = Bspan{Chrom: chrspan.Chrom, Span: intervalset.Span{Min: rawpos, Max: rawpos + span.Width()}}
			return
		}
		rawpos -= width
	}
	return
}

func RandomizeSpan(span Bspan, genome Bed, randgen *rand.Rand) Bspan {
	npos := SpanNumPositions(span, genome)
	if npos < 1 {
		fmt.Fprintln(os.Stderr, "empty newspan")
		fmt.Fprintln(os.Stderr, span)
		fmt.Fprintln(os.Stderr, len(genome.Bspans()))
		fmt.Fprintln(os.Stderr, genome.Bspans())
	}
	rawpos := randgen.Intn(npos)
	return Raw2Bspan(rawpos, span, genome)
}

func RandomlyPlace(span Bspan, dest *Bed, genome Bed, randgen *rand.Rand) {
	newspan := RandomizeSpan(span, genome, randgen)
	// fmt.Println("newspan:")
	// fmt.Println(newspan)
	dest.AddBspans(newspan)
}

func Permute(beds Beds, genome Bed, randgen *rand.Rand, maxComps int) (ovls Overlaps) {
	var new_beds Beds
	for _, bed := range beds {
		new_bed := MakeBed(bed.Name)
		bspans := AllBedSpans(bed)
		for _, bspan := range bspans {
			RandomlyPlace(bspan, &new_bed, genome, randgen)
		}
		new_beds = append(new_beds, new_bed)
	}
	ovls = GetOverlaps(new_beds, maxComps)
	// fmt.Println("ovls:")
	// fmt.Println(ovls)
	return
}

func Permutations(beds Beds, genome Bed, iterations int, randgen *rand.Rand, maxComps int) (osets OverlapSets) {
	for i:=0; i<iterations; i++ {
		osets = append(osets, Permute(beds, genome, randgen, maxComps))
	}
	return
}

func CountPermutations(permutations OverlapSets) (counts OverlapCounts) {
	// fmt.Println("permutations:")
	// fmt.Println(permutations)
	for _, overlaps := range permutations {
		for i, overlap := range overlaps {
			if len(counts) <= i {
				counts = append(counts, OverlapCount{Name: overlap.Name, Components: overlap.Components})
			}
			spans := AllBedSpans(overlap.Bed)
			counts[i].Count = append(counts[i].Count, len(spans))
			counts[i].Covered = append(counts[i].Covered, Covered(spans))
		}
	}
	return
}

func FprintPermCounts(w io.Writer, counts OverlapCounts) {
	for _, count := range counts {
		fmt.Fprintf(w, "%v\t%v\t%v\n", count.Name, count.Count, count.Covered)
	}
}

func GetGenome(path string) (b Bed, e error) {
	var r *os.File
	r, e = os.Open(path)
	if e != nil { return }
	return GetBed(r, "genome")
}

func pcount(val int, dist []int) int {
	for i, dval := range dist {
		if val < dval { return i }
	}
	return len(dist)
}

func GetProb(ovl Overlap, count OverlapCount) (p Prob) {
	scount := append([]int{}, count.Count...)
	sort.Ints(scount)
	scovered := append([]int{}, count.Covered...)
	sort.Ints(scovered)
	p.Name = ovl.Name

	count_pcount := pcount(len(AllBedSpans(ovl.Bed)), scount)
	count_prob_num := len(scount) - count_pcount
	p.CountProb = float64(count_prob_num) / float64(len(scount))
	// p.CountProb = float64(len(scount) - pcount(len(AllBedSpans(ovl.Bed)), scount)) / float64(len(scount))

	covered_pcount := pcount(Covered(AllBedSpans(ovl.Bed)), scovered)
	covered_prob_num := len(scovered) - covered_pcount
	p.CoveredProb = float64(covered_prob_num) / float64(len(scovered))
	// p.CoveredProb = float64(len(scovered) - pcount(Covered(AllBedSpans(ovl.Bed)), scovered)) / float64(len(scovered))

	return
}

func GetProbs(ovls Overlaps, counts OverlapCounts) (out Probs) {
	countsmap := make(map[string]OverlapCount)
	for _, count := range counts {
		countsmap[count.Name] = count
	}
	for _, ovl := range ovls {
		out = append(out, GetProb(ovl, countsmap[ovl.Name]))
	}
	return
}

func FprintProbs(w io.Writer, probs Probs) {
	for _, prob := range probs {
		fmt.Fprintf(w, "%v\t%v\t%v\n", prob.CountProb, prob.CoveredProb, prob.Name)
	}
}

func FullCompare(flags Flags) (c Comparison, err error) {
	genome, err := GetGenome(flags.GenomeBedPath)
	if err != nil { return }
	beds, err := GetBeds(flags.BedPaths)
	if err != nil { return }

	if flags.Verbose {
		fmt.Println("inputs:")
		for _, bed := range beds {
			fmt.Println(bed.Name)
			bspans := AllBedSpans(bed)
			WriteBspans(os.Stdout, bspans...)
		}
	}

	c.Overlaps = GetOverlaps(beds, flags.MaxComps)
	if flags.Iterations > 0 {
		randgen := rand.New(rand.NewSource(int64(flags.Rseed)))
		c.Permutations = Permutations(beds, genome, flags.Iterations, randgen, flags.MaxComps)
		c.IterCounts = CountPermutations(c.Permutations)
		c.Probs = GetProbs(c.Overlaps, c.IterCounts)
	}
	return
}
