#!/bin/bash

# Description
# - Calculate concordance between WES psv genotypes and WGS genotypes

POP="EUR"
CTR="BCM"

PSCNDIR="pscn-summary/${POP}-${CTR}"

for WES_PSCN in ${PSCNDIR}/*.psvs.psCN; do
    
    GENE=$(basename ${WES_PSCN})
    WGS_PSVS="data/WGS-psvs/${POP}.psvs.vcf.gz"

    python calculate_psv_concordance_new.py ${WES_PSCN} ${WGS_PSVS} ${GENE}

done

