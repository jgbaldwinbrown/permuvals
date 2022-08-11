package permuvals

import (
	"fmt"
	"testing"
	"github.com/jgbaldwinbrown/go-intervals/intervalset"
	"strings"
)

func TestGetOverlaps(t *testing.T) {
	expected := []Bspan {
		Bspan{"one", intervalset.Span{Min:5, Max:7}},
		Bspan{"one", intervalset.Span{Min:99, Max:105}},
		Bspan{"two", intervalset.Span{Min:3, Max:11}},
	}
	in1 := []Bspan {
		Bspan{"one", intervalset.Span{Min:2, Max:7}},
		Bspan{"one", intervalset.Span{Min:99, Max:110}},
		Bspan{"two", intervalset.Span{Min:0, Max:11}},
	}
	in2 := []Bspan {
		Bspan{"one", intervalset.Span{Min:5, Max:22}},
		Bspan{"one", intervalset.Span{Min:80, Max:105}},
		Bspan{"two", intervalset.Span{Min:3, Max:20}},
		Bspan{"three", intervalset.Span{Min:0, Max:11}},
	}
	in1b := MakeBed("first")
	in1b.AddBspans(in1...)
	in2b := MakeBed("second")
	in2b.AddBspans(in2...)
	ovls := GetOverlaps(Beds{in1b, in2b}, -1)
	ovlbspans := AllBedSpans(ovls[3].Bed)
	for i, b := range ovlbspans {
		if b != expected[i] {
			t.Errorf("actual and expected do not match. Actual: %v. Expected: %v.", ovlbspans, expected)
		}
	}
}


func TestGetProbs(t *testing.T) {
	bactual := MakeBed("actual")
	bactual.AddBspans(MakeBspan("one", 3, 5), MakeBspan("four", 22, 28))
	is := []string{"i1", "i2"}
	oactual := Overlaps{Overlap{bactual, is}}

	b1 := MakeBed("actual")
	b1.AddBspans(MakeBspan("one", 3, 5))

	b2 := MakeBed("actual")
	b2.AddBspans(MakeBspan("one", 3, 5), MakeBspan("four", 22, 28), MakeBspan("five", 0, 100))

	permutations := OverlapSets {
		Overlaps{Overlap{b1, is}},
		Overlaps{Overlap{b2, is}},
	}
	counts := CountPermutations(permutations)
	probs := GetProbs(oactual, counts)

	expected := []Prob {
		Prob{"actual", .5, .5},
	}
	for i, prob := range probs {
		if prob != expected[i] {
			t.Errorf("actual probs do not match expected. Actual: %v. Expected: %v.", probs, expected)
		}
	}
}

func TestFullCompare(t *testing.T) {
	flags := Flags{BedPaths: "bpaths.txt", GenomeBedPath: "g.bed", Iterations: 10, Rseed: 0}
	comp, err := FullCompare(flags)
	if err != nil {
		panic(err)
	}
	for _, prob := range comp.Probs {
		if prob.CountProb > .3 || prob.CoveredProb > .3 {
			t.Errorf("actual probs do not match expected. Actual: %v.", prob)
		}
	}
}

func TestFprintOvlsBed(t *testing.T) {
	in1 := []Bspan {
		Bspan{"one", intervalset.Span{Min:2, Max:7}},
		Bspan{"one", intervalset.Span{Min:99, Max:110}},
		Bspan{"two", intervalset.Span{Min:0, Max:11}},
	}
	in2 := []Bspan {
		Bspan{"one", intervalset.Span{Min:5, Max:22}},
		Bspan{"one", intervalset.Span{Min:80, Max:105}},
		Bspan{"two", intervalset.Span{Min:3, Max:20}},
		Bspan{"three", intervalset.Span{Min:0, Max:11}},
	}
	in1b := MakeBed("first")
	in1b.AddBspans(in1...)
	in2b := MakeBed("second")
	in2b.AddBspans(in2...)
	ovls := GetOverlaps(Beds{in1b, in2b}, -1)

	outbuild := strings.Builder{}
	FprintOvlsBed(&outbuild, ovls)
	fmt.Println(outbuild.String())
}


func expectedBspans() []Bspan {
	return []Bspan {
		Bspan{"one", intervalset.Span{Min:5, Max:7}},
		Bspan{"one", intervalset.Span{Min:99, Max:105}},
		Bspan{"two", intervalset.Span{Min:3, Max:11}},
	}
}

func in1Bspans() []Bspan {
	return []Bspan {
		Bspan{"one", intervalset.Span{Min:2, Max:7}},
		Bspan{"one", intervalset.Span{Min:99, Max:110}},
		Bspan{"two", intervalset.Span{Min:0, Max:11}},
	}
}

func in2Bspans() []Bspan {
	return []Bspan {
		Bspan{"one", intervalset.Span{Min:5, Max:22}},
		Bspan{"one", intervalset.Span{Min:80, Max:105}},
		Bspan{"two", intervalset.Span{Min:3, Max:20}},
		Bspan{"three", intervalset.Span{Min:0, Max:11}},
	}
}

func genomeBspans() []Bspan {
	return []Bspan {
		Bspan{"one", intervalset.Span{Min:0, Max:200}},
		Bspan{"two", intervalset.Span{Min:0, Max:300}},
	}
}

func toBed(name string, spans []Bspan) Bed {
	b := MakeBed(name)
	b.AddBspans(spans...)
	return b
}

func TestSpanNumPositions(t *testing.T) {
	span := MakeBspan("one", 0, 199)
	genome := toBed("genome",genomeBspans())
	npos := SpanNumPositions(span, genome)
	if npos != 104 {
		t.Errorf("npos (%v) not equal to 104", npos)
	}
}

func TestSpanNumPositions2(t *testing.T) {
	span := MakeBspan("two", 0, 298)
	genome := toBed("genome",genomeBspans())
	npos := SpanNumPositions(span, genome)
	if npos != 3 {
		t.Errorf("npos (%v) not equal to 3", npos)
	}
}
