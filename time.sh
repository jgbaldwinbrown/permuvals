#!/bin/bash
set -e

go build intersect.go

time ( ./intersect -g g.bed -b <(ls b*.bed | grep 'b[1-3]') -i 5000 -r 0 >/dev/null)
time ( ./intersect -g g.bed -b <(ls b*.bed | grep 'b[14]') -i 5000 -r 0 >/dev/null)
time ( ./intersect -g g2.bed -b <(ls b*.bed | grep 'b[14]') -i 5000 -r 0 >/dev/null)
time ( ./intersect -g g3.bed -b <(ls b*.bed | grep 'b[14]') -i 50000 -r 0 >/dev/null)
time ( ./intersect -g g3.bed -b <(ls b*.bed | grep 'b[14]') -i 500000 -r 0 >/dev/null)
time ( ./intersect -g g3.bed -b <(ls b*.bed | grep 'b[14]') -i 5000000 -r 0 >/dev/null)
time ( ./intersect -g g3.bed -b <(ls b*.bed | grep 'b[14]') -i 50000000 -r 0 >/dev/null)
