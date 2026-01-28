# Permuvals

A Go library and executable for permuting the positions of intervals within a
genome and testing for deviations from a random permutation.

## Introduction

This library takes all of the spans specified in a set of bed files, permutes
the positions of those spans, and then asks if the amount of overlap between
spans from one bed to another is more than expected by chance. It tests the
probability of the degree of overlap for every possible set of input bed files.

## Installation

To install the command line executable, just run:

```sh
./install.sh /path/to/bin/dir/
```

## Usage

```
Usage of ./permute_intervals:
  -b string
    	File containing paths to all bed files to compare
  -c	Output raw overlap counts from each permutation
  -g string
    	Bed file containing the lengths of all chromosomes
  -i int
    	Number of permutation iterations to perform (default -1)
  -m int
    	Maximum number of beds to compare at once (default 4)
  -p string
    	comma-separated list of 0-indexed indices of beds to permute (default all)
  -r int
    	Random seed for permutations (default 0)
  -v	Print much more information while running
```

## Library

The library is documented internally and can be imported as follows:

```go
import (
	"github.com/jgbaldwinbrown/permuvals/pkg"
)
```
