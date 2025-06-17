import os
import sys
import glob
import pyreadr
import asyncio
import argparse
import subprocess
import pandas as pd

from multiprocessing import Pool

# -----------------------------------------------------------------------------
# Procedures
# -----------------------------------------------------------------------------

def proc_compute(sample_idx, num_samples, datadir):
    """
    Run ExomeDepthCompute.r scripts in parallel
    """
    print(f'Running ExomeDepth function to build reference sets [{sample_idx}/{num_samples}]')
    cwd = os.path.dirname(os.path.abspath(__file__))
    try:
        subprocess.check_call(
            [f'Rscript {cwd}/run_ExomeDepthCompute_parallel.r -i {sample_idx} -d {datadir}'], 
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            shell=True
        )
    except subprocess.CalledProcessError:
        return sample_idx
    

def proc_merge(datadir, ret=False):
    """
    If executed in parallel, merge individual files into one combined file.
    """
    
    outfile = os.path.join(datadir, "all.stats.tsv")
    all_tsvs = glob.glob(os.path.join(datadir, 'stat_params_*.tsv'))

    # Merge all tsv files together into one
    df = pd.concat([pd.read_csv(f, sep='\t') for f in all_tsvs], ignore_index=True)
    df.sort_values(by='sample', inplace=True)
    df.to_csv(outfile, sep='\t', index=False)

    # Clean up individual RDS files after merging
    for tsv_f in all_tsvs:
        if os.path.isfile(tsv_f):
            os.remove(tsv_f)

    # Return filepath to merged counts file if True
    if ret:
        return outfile 
        

def run(inp):
    
    # Get a list of input BAM filepaths
    with open(inp.input_list, "r") as listfile:
        f_list = listfile.read().splitlines()

    # Set up
    MAX_PROCESSES = int(inp.threads)
    ALL_CNT_DIR = inp.all_cnts_dir
    call_dir = os.path.join(ALL_CNT_DIR, 'ExomeDepth-Calls')
    os.makedirs(call_dir, exist_ok=True)

    print("This process might take some time.\n")

    # Run multiple processes (per bam files)
    with Pool(processes=MAX_PROCESSES) as pool:
        bc_pool_objs = [pool.apply_async(proc_compute, args=(s_idx+1, len(f_list), ALL_CNT_DIR)) for s_idx in range(len(f_list))]
        bc_ret = [obj.get() for obj in bc_pool_objs]

    # Merge the output files into one file
    proc_merge(ALL_CNT_DIR)
    
