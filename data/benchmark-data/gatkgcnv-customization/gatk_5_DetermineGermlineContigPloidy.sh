#!/bin/bash

POP="EUR"
for CTR in BCM BGI BI WUGSC; do
    gatk DetermineGermlineContigPloidy \
      -L exons.hg38.${POP}.${CTR}.filt.interval_list \
      --interval-merging-rule OVERLAPPING_ONLY \
      --input input.readcounts.${POP}.${CTR}.tsv.list \
      --contig-ploidy-priors contig_ploidy_priors.${POP}.${CTR}.tsv \
      --output results \
      --output-prefix ploidy-${POP}-${CTR}
done

