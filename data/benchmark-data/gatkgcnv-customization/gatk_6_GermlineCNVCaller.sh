#!/bin/bash

POP='EUR'
for CTR in BCM BGI BI WUGSC; do
    time gatk --java-options "-Djava.util.concurrent.ForkJoinPool.common.parallelism=16" GermlineCNVCaller \
      --run-mode COHORT \
      -L exons.hg38.${POP}.${CTR}.filt.interval_list \
      --interval-merging-rule OVERLAPPING_ONLY \
      --contig-ploidy-calls results/ploidy-${POP}-${CTR}-calls \
      --input input.readcounts.${POP}.${CTR}.tsv.list \
      --output results \
      --output-prefix cohort-${POP}-${CTR} \
      --max-copy-number 10
done

