package cpermuvals

//#include "bspan.h"
import "C"

import (
	"github.com/jgbaldwinbrown/permuvals/pkg"
)

func GoBspan(cs C.bspan) permuvals.Bspan {
	var gs permuvals.Bspan
	gs.Chrom = C.GoString(cs.chrom)
	gs.Min = int(cs.span.min)
	gs.Max = int(cs.span.max)
	return gs
}

func CBspan(gs permuvals.Bspan) C.bspan {
	var cs C.bspan
	cs.chrom = C.CString(gs.Chrom)
	cs.span.min = C.longlong(gs.Min)
	cs.span.max = C.longlong(gs.Max)
	return cs
}
