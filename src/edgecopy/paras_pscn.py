import os
import sys
import pandas as pd
import numpy as np
import argparse
import shutil
from glob import glob

from . import wes_psvs as wp
from . import estimate_pscn as ep

# Order of procedures
# 1. python run_parascopy_cn.py
# 2. python direct_pscn.py gene_psv gene_hmm > gene.psvs.new.pscn

def generate_wes_psvs(args):
    """
    Run parascopy cn to generate psvs.vcf.gz on exome bam files
    """
    print("\n(1) Running Parascopy-cn to generate WES PSVs.")

    # Artificial depth files
    with open(args.input,'r') as inp_f:
        sample_list = [line.strip().split('::')[1] for line in inp_f]
        outfile_path = (args.output + '/depth.csv')
        wp.generate_depth_csv(sample_list, outfile_path)

    # Get PSV counts on WES BAM files
    wp.get_psv_counts(args)
    
    # Clean up to save space
    print("Cleaning up pooled reads.")
    pooled_reads = glob(f'{args.output}/loci/*/pooled_reads/')
    for pr_dir in pooled_reads:
        shutil.rmtree(pr_dir)

    loci_dirs = os.path.join(args.output, os.path.join('loci', '*'))
    rm_dirs = glob(os.path.join(loci_dirs, 'bed')) + glob(os.path.join(loci_dirs, 'extra'))
    for rm_dir in rm_dirs:
        shutil.rmtree(rm_dir)

    rm_files = glob(os.path.join(loci_dirs, 'res.*'))
    for rm_f in rm_files:
        os.remove(rm_f)

    return os.path.join(args.output, 'psvs.vcf.gz')


def estimate_pscn(loci_list, agcn_dir, pscn_dir, gene, bias_fp=None):
    """
    Estimate psCN for each PSV (per gene)
    - Use bias parameters from different population but same center
    """
    gene_ = gene.split('_refcn')[0]
    with open(loci_list, 'r') as inpf:
        lines = inpf.read().splitlines()
    loci = [l.split('\t')[-1] for l in lines]

    if gene_ not in loci:
        return

    print(f"\nEstimating psCNs for each PSV of {gene_}.")
    
    if bias_fp:
        print(f" - Estimating with bias parameters: {bias_fp}.")
    else:
        print(f" - Estimating without bias parameters.")
    
    psv_vcf_fp = os.path.join(pscn_dir, f'loci/{gene_}/psvs.vcf.gz')
    hmm_out_fp = os.path.join(agcn_dir, f'{gene}/{gene}.hmm.out')
    pscn_out_fp = os.path.join(pscn_dir, f'loci/{gene_}/{gene}.psvs.psCN') 

    ep.main_function(psv_vcf_fp, hmm_out_fp, bias_file=bias_fp, outfile=pscn_out_fp)
    print(f"{gene} completed.")
    print(f'Saved to {pscn_out_fp}')


def run(args):
    
    # Make output directory if it does not already exist
    os.makedirs(args.output, exist_ok=True)    
    
    # (1) generate WES PSVs
    wes_psv = generate_wes_psvs(args)
    
    # (2) Estimate psCN per gene (with or without bias parameter)
    bias_fp = None
    if args.bias_fp:
        bias_fp = args.bias_fp

    for gene_dir in sorted(glob(f'{args.agcn_dir}/*/')):
        gene = os.path.basename(os.path.normpath(gene_dir)) 

        if gene == 'cc-objs':
            continue
    
        estimate_pscn(args.loci_list, args.agcn_dir, args.output, gene, bias_fp=bias_fp)

    # Clean up a number of Parascopy files that we don't need for Edgecopy
    rm_files1 = glob(os.path.join(args.output, f'res.*'))
    for f in rm_files1:
        os.remove(f)

    shutil.rmtree(os.path.join(args.output, 'model'))

