#!/bin/bash
set -e

go build intersect.go

./intersect -g g.bed -b <(ls b*.bed | grep 'b[1-3]') -i 5000 -r 0
./intersect -g g.bed -b <(ls b*.bed | grep 'b[14]') -i 5000 -r 0
./intersect -g g2.bed -b <(ls b*.bed | grep 'b[14]') -i 5000 -r 0
./intersect -g g3.bed -b <(ls b*.bed | grep 'b[14]') -i 50000 -r 0
