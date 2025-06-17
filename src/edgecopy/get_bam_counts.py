import os
import sys
import glob
import pyreadr
import asyncio
import argparse
import subprocess
import pandas as pd

from multiprocessing import Pool
from . import utilities as ut

# -----------------------------------------------------------------------------
# Procedures
# -----------------------------------------------------------------------------

def proc_count(bam_fp, sample_id, exons, outdir):
    """ Run the R script: run_ExomDepthCount.r """
    
    print(f'Running ExomeDepth function to count reads from BAM files [{sample_id}]')
    cwd = os.path.dirname(os.path.abspath(__file__))
    try:
        subprocess.check_call(
            [f'Rscript {cwd}/run_ExomeDepthCount.r -s {bam_fp} -o {outdir} -p {sample_id} -x {exons}'], 
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            shell=True
        )
    except subprocess.CalledProcessError:
        return sample_id
    

def proc_merge(outdir, gene_specific=False, ret=False):

    # Merge all individual RDS files
    rds_files = glob.glob(f"{outdir}/counts_df_*.rds")
    rds_files.sort()
    rds_list = [pyreadr.read_r(f)[None] for f in rds_files]

    filename = "all.counts.tsv"
    if gene_specific:
        filename = "gene.counts.tsv"
   
    outfile = os.path.join(outdir, filename)
    outfile = f"{outdir}/{filename}"
    counts = pd.concat(rds_list, axis=1)
    counts.astype(int).to_csv(outfile, sep="\t", index=None)
   
    # Clean up individual RDS files after merging
    for rds_f in rds_files:
        if os.path.isfile(rds_f):
            os.remove(rds_f)

    # Return filepath to merged counts file if True
    if ret:
        return outfile 
        

def run(inp):
    
    # Get a list of input BAM filepaths and corresponding sample_ids
    with open(inp.input_list, "r") as listfile:
        f_list = listfile.read().splitlines()

    # A list of tuples (bam_fp, sample_id)
    f_list = [(f.split("::")[0], f.split("::")[1]) for f in f_list]
    
    # Set up
    MAX_PROCESSES = int(inp.threads)
    RESULTS_DIR = inp.all_cnts_dir
    EXONS_FP = inp.exon_list
    
    if not os.path.isdir(RESULTS_DIR):
        os.makedirs(RESULTS_DIR)
    
    # Change into correct format, if necessary
    with open(EXONS_FP, 'r') as f:
        first_line = f.readline().strip()
    
    if len(first_line.split('\t'))!=4 or first_line!='#chr\tstart\tend\tname':
        new_exons = ut.read_bed(EXONS_FP)
        EXONS_FP  = ut.add_suffix(EXONS_FP, 'named')
        new_exons.to_csv(EXONS_FP, sep='\t', index=None)
        inp.exon_list = EXONS_FP

    # Run multiple processes (per bam files)
    with Pool(processes=MAX_PROCESSES) as pool:
        bc_pool_objs = [pool.apply_async(proc_count, args=(bam_fp, s_id, EXONS_FP, RESULTS_DIR)) for bam_fp,s_id in f_list]
        bc_ret = [obj.get() for obj in bc_pool_objs]

    # Merge the output files into one file
    proc_merge(RESULTS_DIR)
    
