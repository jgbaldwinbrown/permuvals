package main

import (
	"bufio"
	"os"
	"github.com/jgbaldwinbrown/permuvals/pkg"
)

func main() {
	flags := permuvals.GetFlags()
	w := bufio.NewWriter(os.Stdout)
	defer w.Flush()

	comp, err := permuvals.FullCompare(flags)
	if err != nil { panic(err) }
	permuvals.FprintOvlsBed(w, comp.Overlaps)
	if flags.Iterations > 0 {
		permuvals.FprintProbs(w, comp.Probs)
	}
}
