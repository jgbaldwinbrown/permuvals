#!/bin/bash
set -e

go build

./permute_intervals -g g.bed -b <(ls b*.bed | grep 'b[1-3]') -i 5000 -r 0
./permute_intervals -g g.bed -b <(ls b*.bed | grep 'b[14]') -i 5000 -r 0
./permute_intervals -g g2.bed -b <(ls b*.bed | grep 'b[14]') -i 5000 -r 0
./permute_intervals -g g3.bed -b <(ls b*.bed | grep 'b[14]') -i 50000 -r 0
