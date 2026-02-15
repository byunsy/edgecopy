import os
import sys
import glob
import pyreadr
import asyncio
import argparse
import subprocess
import pandas as pd

from multiprocessing import Pool
from . import exomedepth_refsets as er

# -----------------------------------------------------------------------------
# Procedures
# -----------------------------------------------------------------------------       

def run(inp):
    
    count_matrix_file = os.path.join(inp.all_cnts_dir, 'all.counts.nondup.tsv')
    interval_file = os.path.join(inp.all_cnts_dir, 'all.counts.nondup.meta.tsv')

    out_fp = inp.all_stat_fp
    out_fp_log_file = open(out_fp + '.log', 'w') 
    
    er.referenceset_fit(count_matrix_file, interval_file, inp.input_samples, out_fp,
                        logfile=out_fp_log_file, exons_for_beta=10000, max_ref=30, threads=int(inp.threads))
                        
    out_fp_log_file.close()
