import re
import os
import sys
import numpy as np
import pandas as pd
import subprocess

from copy import deepcopy
from multiprocessing import Pool

from . import utilities as ut
from . import run_all as _run


def prep(inputwrapper, args):
    """ Prepares data for main analysis """

    # Check if gene has multiple unique refCN
    # - if so, create unique class for each copy
    # - and remove original
    to_append, to_remove = [], []
    for input_info in inputwrapper:
        if input_info.unique_cn:
            for info in input_info.unique_cn:
                _c, _s, _e = re.split(':|-', info[0]) 
                new_input_info = _run.InputInfo(args, ((_c,_s,_e), info[1]))
                new_input_info.refcn = info[2]
                new_input_info.gene_cnts_fp = info[3]
                new_exons_fp = ut.split_exon_file(input_info, new_input_info)
                new_input_info.gene_exon_fp = new_exons_fp
                to_append.append(new_input_info)
            to_remove.append(input_info)
  
    for old_input_info in to_remove:
        inputwrapper.remove(old_input_info)
    inputwrapper += to_append
   
    # Create output directory
    os.makedirs(inputwrapper[0].outdir, exist_ok=True)
    
    # update inputwrappers
    for inp in inputwrapper:
        inp.update(args)
    
    # Check quality, coverage, refcn, copy-number variance
    highrefcn_genes = []
    highrefcn = 8
    
    for info in inputwrapper:
        if int(info.refcn) >= highrefcn:
            highrefcn_genes.append(info.loci_name)

    lowqual_genes = []
    highcov_genes, lowcov_genes = ut.check_coverage(inputwrapper)
    lowcn_genes, highcn_genes = ut.check_cn_range(inputwrapper, args.high_refcn, args.maxcn)
    
    good_genes = list((set(highcov_genes) & set(lowcn_genes)) - set(lowqual_genes + highrefcn_genes))

    final_list_out = os.path.join(inputwrapper[0].outdir, "summary.input.info.tsv")
    stratified_out = os.path.join(inputwrapper[0].outdir, "summary.gene.category.tsv")
    
    final_list = ut.prep_final_list(good_genes, inputwrapper, final_list_out)
    gene_strat = ut.stratify_genes(inputwrapper, lowqual_genes, lowcov_genes, 
                                   highcn_genes, highrefcn_genes, good_genes, stratified_out)

    print("\nprep_data.py completed.")
    return final_list, gene_strat


if __name__ == '__main__':
    prep()
