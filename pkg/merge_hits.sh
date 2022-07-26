#!/bin/bash
set -e

awk -F "\t" -v OFS="\t" '$'${1}' >= '${2} \
> ${3}_thresholded.bed

bedtools merge -i ${3}_thresholded.bed > ${3}_thresh_merge.bed
