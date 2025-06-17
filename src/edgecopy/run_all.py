import re
import os
import sys
import time
import json
import pickle
import numpy as np
import pandas as pd

import argparse
from argparse import RawTextHelpFormatter

from multiprocessing import Pool
from glob import glob

from . import get_bam_counts as _bamcounts
from . import get_bam_stats as _bamstats
from . import paras_pool as _pool
from . import paras_examine as _examine
from . import paras_pscn as _pscn
from . import prep_data as _prep
from . import analyze as _analyze
from . import headers as _h
from . import evaluate as _eval
from . import utilities as _util

# ----------------------------------------------------------------------------
# CLASS OBJECT FOR STORING INPUT INFO
# ----------------------------------------------------------------------------
class InputInfo:
    """
    Class for storing input information of a particular locus
    """    
    
    def __init__(self, args, loci_info):
       
        # User input variables
        self.input_list = args.input
        self.outdir = args.output
        self.loci_name = loci_info[1]
        self.loci_region = loci_info[0]
        self.reference = args.reference
        self.exon_list = args.exon_list
        self.threads = args.threads
        self.nondup = False

        self.all_cnts_dir = os.path.join(args.output, 'all')
        self.exon_dir = os.path.join(self.outdir, 'exons')
        
        if self.loci_name != 'depth':
            self.loci_region_tab = loci_info[0]
            self.loci_region_str = f'{loci_info[0][0]}:{loci_info[0][1]}-{loci_info[0][2]}' 
            self.cnts_dir = os.path.join(args.depth, f'counts-{loci_info[1]}')
            self.all_cnts_dir = os.path.join(args.depth, 'all')
            self.refcn_fp = os.path.join(args.depth, 'genes.167.refCN.list')
            self.hg_version = args.hgver
            self.exon_dir = os.path.join(args.depth, 'exons')
            self.exon_list = os.path.join(self.all_cnts_dir, 'all.counts.meta.tsv')

        # Fixed directories and filepaths
        self.cwd = os.path.dirname(os.path.abspath(__file__))
        self.data_dir = os.path.join(self.cwd, 'data')
        
        self.hom_table = args.hom_table

        self.allexons_fp = os.path.join(self.exon_dir, 'allexons.parascopy.examine.output')
        self.all_cnts_fp = os.path.join(self.all_cnts_dir, 'all.counts.tsv')
        self.all_meta_fp = os.path.join(self.all_cnts_dir, 'all.counts.meta.tsv')
        self.all_stat_fp = os.path.join(self.all_cnts_dir, 'all.stats.tsv')
        self.inpwrap_dir = os.path.join(self.outdir, "inputwrappers")
        
        self.gene_cnts_fp = ''
        self.gene_exon_fp = ''
        self.unique_cn = []

        # Get sample names
        with open(self.input_list, 'r') as inp:
            if '::' not in inp.readline():
                e = ("Sample names not found in input list. " + 
                     "Please specify sample names separated by '::'. " +
                     "For example, /path/to/file.bam::NA00001")
                raise IndexError(e)

            samples = [s.split('::')[1] for s in inp.read().splitlines()]
        self.input_samples = samples

    def update(self, args):
        """
        Update inputwrappers with new parameters for agcn
        """
        self.input_list = args.input
        self.ccobjdir = os.path.join(args.output, 'cc-objs')
        self.finaldir = os.path.join(args.output, self.loci_name)
        self.debug_fp = os.path.join(self.finaldir, f'{self.loci_name}.hmm.debug.log')
        self.debug_mode = args.debug
        
        if not args.priors:
            args.priors = os.path.join(self.data_dir, 'priors')
        self.priordir = args.priors
        self.priors_fp = os.path.join(self.priordir, f'{self.loci_name}.prior')
        #self.priors_fp = os.path.join('data/priors', f'default.dup.{self.refcn}.prior')
        
        self.all_cnts_fp = os.path.join(self.all_cnts_dir, 'all.counts.tsv')
        self.all_stat_fp = os.path.join(self.all_cnts_dir, 'all.stats.tsv')
        
        self.qual_thresh = args.qual_thresh
        self.t_prob = args.t_prob
        self.min_cc_size = args.min_cc_size
        self.max_iter = args.max_iter
        self.threads = args.threads
        self.skiphmm = args.skip_hmm

# ----------------------------------------------------------------------------
# ARGUMENT PARSERS
# ----------------------------------------------------------------------------

def parse_args_depth(in_argv):
    """ 
    Argument parser for depth pipeline
    """
    # Parse arguments from command-line
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, add_help=False)
    
    # Required arguments
    parser.add_argument("-i", "--input", required=True, 
                        help="Path to input BAM file list.\n"
                             "All entries should follow the format 'filename[::sample]'.\n\n")
    
    parser.add_argument("-o", "--output", required=True, 
                        help="Path to output directory.\n\n")
    
    parser.add_argument("-r", "--reference", required=True, 
                        help="Path to reference genome file.\n\n")
   
    parser.add_argument("-x", "--exon-list", required=True, 
                        help="Path to list of all exons.\n"
                             "All entries should be tab-separated,\n"
                             "following the format 'chr\\tstart\\tend\\tname'.\n\n")
   
    parser.add_argument("-t", "--hom-table", required=True, 
                        help="Path to homology table.\n\n")
    
    # Optional arguments
    parser.add_argument("-@", "--threads", default=8, 
                        help="Number of threads.\n\n")
    
    parser.add_argument("-h", "--help", action='help',
                        help="Show this help message.\n\n")
    
    args_parsed = parser.parse_args(in_argv)
    return args_parsed


def parse_args_agcn(in_argv):
    """
    Argument parser for agcn pipeline
    """
    # Parse arguments from command-line
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, add_help=False)

    # Required arguments
    parser.add_argument("-i", "--input", required=True, 
                        help="Path to input BAM file list.\n"
                             "All entries should follow the format 'filename[::sample]'.\n"
                             "Should match the argument used for exparascopy depth.\n\n")
    
    parser.add_argument("-d", "--depth", required=True, 
                        help="Path to exome-parascopy depth directory\n\n")
    
    parser.add_argument("-o", "--output", required=True, 
                        help="Path to output directory\n\n")
    
    loci_group = parser.add_mutually_exclusive_group(required=True)
    
    loci_group.add_argument("-l", "--locus",
                            help="A single locus following the format 'chr:start-end::name'.\n"
                                 "For example, 'chr5:70925087-70953015::SMN1'.\n"
                                 "Or one could pass in multiple regions, separating each with commas.\n"
                                 "For example, 'chr5:70925087-70953015::SMN1,chr7:74768962-74856551::NCF1'.\n\n")
    loci_group.add_argument("-L", "--loci-list",
                            help="Path to a file containing a list of loci and positions.\n"
                                 "Each region should be separated by a newline and,\n" 
                                 "all entries should match the format 'chr:start-end name'.\n"
                                 "For example, 'chr5:70925087-70953015 SMN1'.\n\n")
    
    parser.add_argument("-r", "--reference", required=True, 
                        help="Path to reference genome file.\n\n")
    
    parser.add_argument("-t", "--hom-table", required=True, 
                        help="Path to homology table.\n\n")
    
    # Optional arguments
    parser.add_argument("-x", "--exon-list", required=False, 
                        help="Path to a BED file containing a list of all exons.\n"
                             "All entries should be tab-separated,\n"
                             "following the format 'chr\\tstart\\tend\\tname'.\n"
                             "Should match the argument used for exparascopy depth.\n\n")
    
    parser.add_argument("-p", "--priors", required=False, 
                        help="Path to prior directory.\n\n")
    
    parser.add_argument("-@", "--threads", type=int, default=8, 
                        help="Number of threads.\n\n")
    
    parser.add_argument("--qual-thresh", type=int, default=20, 
                        help="Quality threshold for trusting CN estimates.\n\n")
    
    parser.add_argument("--high-refcn", type=int, default=8, 
                        help="Threshold for the highest reference copy number to examine.\n"
                             "Any gene with reference copy number higher will be skipped.\n\n")
    
    parser.add_argument("--maxcn", type=int, default=10, 
                        help="Threshold for the maximum copy number state.\n"
                             "Any gene with reference copy number higher will be skipped.\n\n")
    
    parser.add_argument("--t-prob", type=float, default=0.0001, 
                        help="Transition probability.\n"
                             "Default set at 0.0001.\n\n")
   
    parser.add_argument("--min-cc-size", type=int, default=5, 
                        help="Minimum size for a valid connected component.\n"
                             "Default set at 5.\n\n")
    
    parser.add_argument("--max-iter", type=int, default=5, 
                        help="Maximum number of iteration for estimating HMM path.\n"
                             "Default set at 5.\n\n")
    
    parser.add_argument("--hgver", type=str, default='hg38', 
                        help="hg version [hg38 | hg19]\n\n")
    
    parser.add_argument("--debug", action='store_true', 
                        help="Turn on debug mode for HMM.\n\n")
    
    parser.add_argument("--skip-hmm", action='store_true', 
                        help="Skip HMM estimation.\n\n")
    
    parser.add_argument("-h", "--help", action='help',
                        help="Show this help message.\n\n")
    
    args_parsed = parser.parse_args(in_argv)
    return args_parsed


def parse_args_pscn(in_argv):
    """ 
    Argument parser for pscn pipeline
    """
    # Parse arguments from command-line
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, add_help=False)
    
    # Required arguments
    parser.add_argument("-i", "--input", required=True, 
                        help="Path to input BAM file list.\n\n")
    parser.add_argument("-o", "--output", required=True, 
                        help="Path to output directory.\n\n")
    parser.add_argument("-l", "--loci-list", required=True,
                        help="Path to list of loci and position.\n\n")
    parser.add_argument("-f", "--reference", required=True, 
                        help="Path to reference genome file.\n\n")
    parser.add_argument("-t", "--hom-table", required=True, 
                        help="Path to homology table.\n\n")
    parser.add_argument("-a", "--agcn-dir", required=True, 
                        help="Path to agCN directory.\n\n")
    
    # Optional arguments
    parser.add_argument("-b", "--bias-fp", required=False, 
                        help="Path to output bias parameter.\n\n")
    parser.add_argument("-@", "--threads", default=8, 
                        help="Number of threads.\n\n")
    parser.add_argument("-h", "--help", action='help',
                        help="Show this help message.\n\n")
    
    args_parsed = parser.parse_args(in_argv)
    return args_parsed


def parse_args_eval(in_argv):
    """ 
    Argument parser for eval pipeline
    """
    # Parse arguments from command-line
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    
    # Required arguments
    parser.add_argument("-m", "--hmm-dir", required=True, 
                        help="path to HMM output directory")
    parser.add_argument("-t", "--truth-dir", required=True, 
                        help="path to ground-truth directory")
    parser.add_argument("-r", "--refcn-fp", required=True, 
                        help="path to list of genes with refCN (tab-sep)")
    
    args_parsed = parser.parse_args(in_argv)
    return args_parsed

# ----------------------------------------------------------------------------
# PIPELINES: depth --> agcn --> pscn --> eval
# ----------------------------------------------------------------------------

def run_depth_pipeline(inner_argv):
    """
    (1) Obtain read counts across all exons
    (2) Build reference set (use ExomeDepth's method)
    """
    
    # Obtain input information
    args = parse_args_depth(inner_argv)
    
    # Input information for each locus
    inputwrapper = InputInfo(args, ('init','depth'))
    
    # (1) Obtain all.counts.tsv
    _h.make_header("(1) Obtain all-counts matrix.")
    if not os.path.isfile(inputwrapper.all_cnts_fp):
        _bamcounts.run(inputwrapper)
    else:
        print(f"Found an existing all-counts matrix file: {inputwrapper.all_cnts_fp}.")

    # (2) Prepare allexons parascopy examine file
    #     - examines reference CNs of all exons based on homology table
    _h.make_header("(2) Examine reference CNs of all exons.")
    
    if not os.path.isfile(inputwrapper.allexons_fp):
        _examine.check_cnts_and_exons(inputwrapper)
        _examine.make_allexons(inputwrapper)
    else:
        print(f"Found an existing parascopy-examine file: {inputwrapper.allexons_fp}.")
    
    # (3) Build reference sets using ExomeDepth
    _h.make_header("(3) Build reference sets.")
    if not os.path.isfile(inputwrapper.all_stat_fp):
        _bamstats.run(inputwrapper)
    else:
        print(f"Found an existing all-stats file: {inputwrapper.all_stat_fp}.")
    
    # Save objects for later use in analysis pipeline
    wrapper_dir = inputwrapper.inpwrap_dir
    os.makedirs(wrapper_dir, exist_ok=True)

    print(f"Saving inputwrapper information to {wrapper_dir}.")
    inp_fp = os.path.join(wrapper_dir, f'{inputwrapper.loci_name}.obj')
    with open(inp_fp, 'wb') as outf:
        pickle.dump(inputwrapper, outf, protocol=pickle.HIGHEST_PROTOCOL)
    

def run_agcn_pipeline(inner_argv):
    """
    (1) Obtain read counts for duplicated exons of duplicated genes
    (2) Jointly analyze point CN estimates for all samples
    (3) Jointly analyze HMM-path CN estimates for all samples
    """
    
    # Obtain input information
    args = parse_args_agcn(inner_argv)
    
    # Passed in list of loci
    if args.loci_list:
        l_df = _util.read_bed(args.loci_list)
        loci_list = [[l[:3], l[3]] for l in l_df.values.tolist()]

    # Passed in a single locus
    else:
        loci_list = []
        for locus in args.locus.split(','):
            l_info = locus.split('::')
            l_info_tup = re.split(':|-', l_info[0])
            
            if len(l_info) >= 2:
                loci_list.append((l_info_tup, l_info[1]))
            else:
                loci_list.append((l_info_tup, '_'.join(l_info_tup)))
    
    # Load inputwrapper objects 
    inputwrapper = [InputInfo(args, loci) for loci in loci_list]

    # Obtain gene.counts.tsv
    _h.make_header("(1) Obtain gene-counts matrix.")
    _pool.run(inputwrapper)
    
    # Prepare data for analysis
    _h.make_header("(2) Prepare data for analysis.")
    final_list, gene_strat = _prep.prep(inputwrapper, args)
    inputwrapper = [inp for inp in inputwrapper if inp.loci_name in final_list]

    # Analyze point and vector (HMM) agCNs
    _h.make_header("(3) Analyze and estimate agCNs.")
    with Pool(processes=args.threads) as pool:
        pool_objs = [pool.apply_async(_analyze.run, args=(input_info,)) for input_info in inputwrapper]
        ret = [obj.get() for obj in pool_objs]
    
    #_util.concat_all_pnt(inputwrapper[0].outdir)
    _util.concat_all_hmm(inputwrapper[0].outdir)

def run_pscn_pipeline(inner_argv):
    """
    (1) Reformat Edgecopy's HMM results into Parascopy input BED
    (2) Pass the BED file into Parascopy to estimate psCNs
    """
    # Obtain input information
    args = parse_args_pscn(inner_argv)
    
    # Analyze psCNs using Parascopy's detect_cn module
    _h.make_header("Analyze and estimate psCNs.")
    _pscn.run(args)


def run_eval_pipeline(inner_argv):
    """
    """
    
    args = parse_args_eval(inner_argv)

    # Evaluate and compare estimated CNs against ground-truth
    # GENE.hmm.bed --> GENE.eval
    _h.make_header("Evaluating point and HMM estimates.")
    _eval.evaluate_exome_paras(args.hmm_dir, args.truth_dir, args.refcn_fp)
    

# ----------------------------------------------------------------------------
# MAIN FUNCTION
# ----------------------------------------------------------------------------
def main():
    
    program = os.path.basename(sys.argv[0])         # edgecopy
    
    # command is too short
    if len(sys.argv) < 2:
        print(_h.show_usage(program))
        return
    
    command = sys.argv[1]                           # depth, agcn, pscn, eval
    inner_argv = sys.argv[2:]                       # parameters
    inner_prog = '{} {}'.format(program, command)

    # Help command
    if command == '-h' or command == '--help' or command == 'help':
        _h.show_usage(program)
        return
    
    # Start Edgecopy
    _h.make_title()
    _h.make_header("Commands used")
    _h.print_command([program] + sys.argv[1:])

    # Take valid commands from user
    if command == 'depth':
        run_depth_pipeline(inner_argv)
        return

    elif command == 'agcn':
        run_agcn_pipeline(inner_argv)
        return

    elif command == 'pscn':
        run_pscn_pipeline(inner_argv)
        return

    elif command == 'eval':
        run_eval_pipeline(inner_argv)
        return
    
    else:
        _h.show_usage(program)
        return


if __name__ == '__main__':
    main()

