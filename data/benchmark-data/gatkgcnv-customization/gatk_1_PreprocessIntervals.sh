#!/bin/bash

time gatk PreprocessIntervals \
  -R reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -L exons.hg38.interval_list \
  --interval-merging-rule OVERLAPPING_ONLY \
  --bin-length 0 \
  --padding 250 \
  -O exons.hg38.preproc.interval_list

