#!/bin/bash
set -e

go build

time ( ./permute_intervals -g g.bed -b <(ls b*.bed | grep 'b[1-3]') -i 5000 -r 0 >/dev/null)
time ( ./permute_intervals -g g.bed -b <(ls b*.bed | grep 'b[14]') -i 5000 -r 0 >/dev/null)
time ( ./permute_intervals -g g2.bed -b <(ls b*.bed | grep 'b[14]') -i 5000 -r 0 >/dev/null)
time ( ./permute_intervals -g g3.bed -b <(ls b*.bed | grep 'b[14]') -i 50000 -r 0 >/dev/null)
time ( ./permute_intervals -g g3.bed -b <(ls b*.bed | grep 'b[14]') -i 500000 -r 0 >/dev/null)
time ( ./permute_intervals -g g3.bed -b <(ls b*.bed | grep 'b[14]') -i 5000000 -r 0 >/dev/null)
time ( ./permute_intervals -g g3.bed -b <(ls b*.bed | grep 'b[14]') -i 50000000 -r 0 >/dev/null)
