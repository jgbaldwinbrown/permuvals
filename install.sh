#!/bin/sh

BINDIR="${1}"

go build cmd/permute_intervals.go
cp permute_intervals "${BINDIR}/"
