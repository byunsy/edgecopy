#!/bin/bash

# skip this flag when analyzing duplicated genes
# --annotated-intervals exons.hg38.annot.tsv \

POP='EUR'
for CTR in BCM BGI BI WUGSC; do
    gatk FilterIntervals \
      -L exons.hg38.preproc.interval_list \
      --interval-merging-rule OVERLAPPING_ONLY \
      -I input.readcounts.${POP}.${CTR}.tsv.list \
      -O exons.hg38.${POP}.${CTR}.filt.interval_list
done

