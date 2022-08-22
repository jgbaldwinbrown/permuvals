package permuvals

import (
	"testing"
	"github.com/jgbaldwinbrown/go-intervals/intervalset"
)

func TestSubFull(t *testing.T) {
	expected := []Bspan {
		Bspan{"one", intervalset.Span{Min:2, Max:7}},
	}
	in1 := []Bspan {
		Bspan{"one", intervalset.Span{Min:2, Max:7}},
		Bspan{"one", intervalset.Span{Min:99, Max:110}},
		Bspan{"two", intervalset.Span{Min:0, Max:11}},
	}
	in2 := []Bspan {
		Bspan{"one", intervalset.Span{Min:80, Max:105}},
		Bspan{"two", intervalset.Span{Min:3, Max:20}},
		Bspan{"three", intervalset.Span{Min:0, Max:11}},
	}
	in1b := MakeBed("first")
	in1b.AddBspans(in1...)
	in2b := MakeBed("second")
	in2b.AddBspans(in2...)
	outb := in1b.SubtractFullsBed(in2b)
	subspans := AllBedSpans(outb)
	if len(subspans ) != len(expected) {
		t.Errorf("actual and expected do not match. Actual: %v. Expected: %v.", subspans, expected)
	}
	for i, b := range subspans {
		if b != expected[i] {
			t.Errorf("actual and expected do not match. Actual: %v. Expected: %v.", subspans, expected)
		}
	}
}
