#!/bin/bash

time gatk AnnotateIntervals \
  -R reference/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -L exons.hg38.preproc.interval_list \
  --interval-merging-rule OVERLAPPING_ONLY \
  --mappability-track filter-tracks/hg38.k100.umap.single.merged.bed.gz \
  --segmental-duplication-track filter-tracks/hg38_seg_dup_track.bed \
  -O exons.hg38.annot.tsv

