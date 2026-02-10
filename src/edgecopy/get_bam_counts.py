import os
import sys
import glob
import pyreadr
import asyncio
import argparse
import subprocess
import pandas as pd
import pysam

from multiprocessing import Pool
from . import utilities as ut
from . import countreads_exonic as cre
from parascopy import examine_dupl

# -----------------------------------------------------------------------------
# Procedures
# -----------------------------------------------------------------------------

def run(inp):
    
    # read bam list (lines like /path/to.bam::SAMPLE)
    with open(inp.input_list, 'r') as fh:
        lines = [l.strip() for l in fh if l.strip()]
    file_list = [l.split("::")[0] for l in lines]

    os.makedirs(inp.all_cnts_dir, exist_ok=True)

    # build homolog bedfile for exonic read counting
    hom_table_dir = os.path.dirname(inp.hom_table)
    output_homolog_bed = os.path.join(hom_table_dir, 'homolog.bed')
    inp.homolog_bed = output_homolog_bed

    # bgzip and tabix-index the exon BED file
    exon_list_gz = f"{inp.exon_list}.gz"
    pysam.tabix_compress(inp.exon_list, exon_list_gz, force=True)
    pysam.tabix_index(exon_list_gz, preset="bed", force=True)

    print("Running Parascopy-examine to determine duplicated regions.")
    examine_dupl.main('parascopy examine', 
                     ['-t', inp.hom_table, '-R', exon_list_gz, '-o', output_homolog_bed])

    # build intervals and exon_list
    trees, exon_list, gene_table = cre.read_bedfile_pysam(bedfile=exon_list_gz, region=".", hombed=output_homolog_bed)

    # process BAMs
    bamstats = cre.process_files_in_parallel(file_list, int(inp.threads), trees, exon_list,
                                             region='.', minMQ=20, reference_file=inp.reference)

    # write count matrices (creates outdir/all/all.counts.tsv and meta)
    print(f"Saving counts matrices to: {inp.all_cnts_dir}")
    cre.print_count_matrix(file_list, exon_list, bamstats, inp.all_cnts_dir, reference_file=inp.reference)

