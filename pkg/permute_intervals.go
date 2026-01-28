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

// Unused
type OvlStats struct {
	Ncomponents int
	MeanCounts float64
	MeanCoverage float64
}

// All output from a round of permutations and comparison with real intervals
type Comparison struct {
	Permutations OverlapSets
	Overlaps Overlaps
	IterCounts OverlapCounts
	Probs Probs
}

// A span as used by a bed file, with a chromosome and a region
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

// Number of basepairs covered by all spans, not taking into account that some spans could overlap
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
	ToPermute []int
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

// All of the intervals that you get when you intersect a set of bed, plus the names of all the beds used in the intersection
type Overlap struct {
	Bed
	Components []string
}

type Overlaps []Overlap
type OverlapSets []Overlaps

// Count is the number of overlapped spans, Covered is the number of basepairs overlapped
type OverlapCount struct {
	Count []int
	Covered []int
	Name string
	Components []string
}

type OverlapCounts []OverlapCount

// The probability of getting the level of overlap by count or covered, when compared true data to the permuted distribution
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

// put all of the lines in a path into a []string
func GetBedpaths(bedpaths_path string) (paths []string, err error) {
	r, err := os.Open(bedpaths_path)
	if err != nil { return }
	s := bufio.NewScanner(r)
	for s.Scan() {
		paths = append(paths, s.Text())
	}
	return
}

// for each line in bedpaths_path, parse a Bed and add it to beds
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

// Just a wrapper for dest.AddBspans
func AddBed(dest *Bed, src Bed) {
	// could hand-code this to be faster
	dest.AddBspans(AllBedSpans(src)...)
}

// Get a unique list of all chromosomes in the bed file
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

// Intersect in place
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

// Subtract in place
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

// Simple power function using big.Int
func BigPow(y, x *big.Int) {
	y.SetInt64(1)
	two := big.NewInt(2)
	one := big.NewInt(1)
	for i:=big.NewInt(0); i.Cmp(x) < 0; i.Add(i, one) {
		y.Mul(y, two)
	}
}

// Generate one of the possible permutations of sets of beds to overlap based
// on its index being split into binary code, and only including the "1" digits
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

// Get all possible permutations of bed files as overlaps
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

// Add up all of the "1" digits in the binary representation of src, put the result in sum
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

// Get all permutations of beds, but only overlapping up to maxComps beds at a time
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

func parseIndices(s string) ([]int, error) {
	fields := strings.Split(s, ",")
	out := make([]int, 0, len(fields))
	for _, f := range fields {
		idx, e := strconv.Atoi(f)
		if e != nil {
			return out, e
		}
		out = append(out, idx)
	}
	return out, nil
}

func GetFlags() (f Flags) {
	flag.BoolVar(&f.Verbose, "v", false, "Print much more information while running")
	flag.StringVar(&f.BedPaths, "b", "", "File containing paths to all bed files to compare")
	flag.StringVar(&f.GenomeBedPath, "g", "", "Bed file containing the lengths of all chromosomes")
	flag.IntVar(&f.Iterations, "i", -1, "Number of permutation iterations to perform")
	flag.IntVar(&f.Rseed, "r", 0, "Random seed for permutations (default 0)")
	toPermuteStrp := flag.String("p", "", "comma-separated list of 0-indexed indices of beds to permute (default all)")
	flag.Parse()
	if f.BedPaths == "" || f.GenomeBedPath == "" {
		panic(fmt.Errorf("Missing path"))
	}

	if *toPermuteStrp != "" {
		var e error
		f.ToPermute, e = parseIndices(*toPermuteStrp)
		if e != nil {
			panic(e)
		}
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

// Unused
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

// Unused
func FprintOvlsStats(w io.Writer, os ...OvlStats) {
	for _, o := range os {
		fmt.Fprintf(w, "%v\t%v\t%v\n", o.Ncomponents, o.MeanCounts, o.MeanCoverage)
	}
}

// Count the number of possible locations you could place a span in a genome
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

// From a raw number indicating where to put the span, generate a new span at the indexed location in the genome
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

// Move a span to a random location somewhere in the genome
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

// RandomizeSpan, and put the result in dest
func RandomlyPlace(span Bspan, dest *Bed, genome Bed, randgen *rand.Rand) {
	newspan := RandomizeSpan(span, genome, randgen)
	// fmt.Println("newspan:")
	// fmt.Println(newspan)
	dest.AddBspans(newspan)
}

// Take all beds, then randomly permute all their span positions, then calculate overlaps for the permuted beds
func Permute(beds Beds, genome Bed, randgen *rand.Rand, maxComps int, toPermute []int) (ovls Overlaps) {
	toperm := make(map[int]struct{}, len(toPermute))
	for _, i := range toPermute {
		toperm[i] = struct{}{}
	}

	var new_beds Beds
	for i, bed := range beds {
		_, ok := toperm[i]
		if ok || len(toPermute) < 1 {
			new_bed := MakeBed(bed.Name)
			bspans := AllBedSpans(bed)
			for _, bspan := range bspans {
				RandomlyPlace(bspan, &new_bed, genome, randgen)
			}
			new_beds = append(new_beds, new_bed)
		} else {
			new_beds = append(new_beds, bed)
		}
	}
	ovls = GetOverlaps(new_beds, maxComps)
	// fmt.Println("ovls:")
	// fmt.Println(ovls)
	return
}

// run Permute as many times as specified in iterations
func Permutations(beds Beds, genome Bed, iterations int, randgen *rand.Rand, maxComps int, toPermute []int) (osets OverlapSets) {
	for i:=0; i<iterations; i++ {
		osets = append(osets, Permute(beds, genome, randgen, maxComps, toPermute))
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

// Find the probability that the true overlap was that much or more by chance
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
		c.Permutations = Permutations(beds, genome, flags.Iterations, randgen, flags.MaxComps, flags.ToPermute)
		c.IterCounts = CountPermutations(c.Permutations)
		c.Probs = GetProbs(c.Overlaps, c.IterCounts)
	}
	return
}

func Full() {
	flags := GetFlags()
	w := bufio.NewWriter(os.Stdout)
	defer w.Flush()

	comp, err := FullCompare(flags)
	if err != nil { panic(err) }
	FprintOvlsBed(w, comp.Overlaps)
	if flags.Iterations > 0 {
		FprintProbs(w, comp.Probs)
	}
}
