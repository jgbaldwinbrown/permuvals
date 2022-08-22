package permuvals

import (
	"github.com/jgbaldwinbrown/go-intervals/intervalset"
	"fmt"
)

func SubtractFull(bspans []Bspan, src *intervalset.Set) []Bspan {
	out := []Bspan{}
	for _, bspan := range bspans {
		fmt.Println("bspan:", bspan)
		fmt.Println("src:", src)
		intersected := src.Copy()
		intersected.Intersect(intervalset.NewSet([]intervalset.Interval{&bspan.Span}))
		intersections := intersected.AllIntervals()
		if len(intersections) == 0 {
			fmt.Println("no overlap")
			out = append(out, bspan)
		}
	}
	return out
}

func (b *Bed) SubtractFullsBed(src Bed) Bed {
	bspans := []Bspan{}
	all_chroms := CombineChroms(*b, src)
	for _, chrom := range all_chroms {
		_, b_has_chrom := b.Intervals[chrom]
		_, src_has_chrom := src.Intervals[chrom]
		if !b_has_chrom {
		} else if !src_has_chrom {
			bspans = append(bspans, AllBspans(chrom, b.Intervals[chrom])...)
		}else {
			goods := AllBspans(chrom, b.Intervals[chrom])
			bads := src.Intervals[chrom]
			bspans = append(bspans, SubtractFull(goods, bads)...)
			// b.Intervals[chrom].Sub(src.Intervals[chrom])
		}
	}
	out := MakeBed(b.Name)
	out.AddBspans(bspans...)
	return out
}
