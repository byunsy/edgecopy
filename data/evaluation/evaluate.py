import os
import sys
import pandas as pd
import numpy as np
import argparse
from glob import glob
from collections import Counter
import cn_per_exon as cp


def get_diffs(estim, truth):
    
    try:
        estim = np.array([int(e) if e!='_' else 0 for e in str(estim).split(',')])
        truth = np.array([int(t) if t!='.' else 0 for t in str(truth).split(',')])
    except ValueError:
        print(estim, truth)
        sys.exit()
    
    diffs = []
    for i,d in enumerate(estim-truth):
        if estim[i]==0 or truth[i]==0:
            diffs.append('.')  
        elif d==0:
            diffs.append('.')
        else:
            diffs.append('X')
            
    return ','.join(diffs)


def is_ref(truth, refCN):
    if all(t=='.' for t in str(truth).split(',')):
        return False
    return Counter([t for t in str(truth).split(',') if t!='.']).most_common(1)[0][0]==str(refCN)


def evaluate_edgecopy(wes_dir, wgs_dir, out_dir, refcn_list_fp, QUAL_THRESH=20, MIN_CC_SIZE=5):
    """
    Evaluate and summarize the Edgecopy's WES HMM estimates using
    WGS estimates from Parascopy
    """
    fps = glob(os.path.join(wes_dir, f'*/*.hmm.out'))
    GENES = [fp.split('/')[-1].split('.hmm')[0] for fp in fps]
    
    REFCN = pd.read_csv(refcn_list_fp, sep='\t', header=None)
    REFCN = REFCN.set_index(0)[1].to_dict()
    
    for GENE in GENES:
        print(GENE)

        fp_wes = os.path.join(wes_dir, os.path.join(GENE, f'{GENE}.hmm.out'))
        fp_wgs = os.path.join(out_dir, os.path.join('wgs-estimates', f'{GENE}.exon.CN.bed'))
        
        refCN = REFCN.get(GENE)

        # Dataframe of WES agCN estimates
        df_wes = pd.read_csv(fp_wes, sep='\t', comment='#')
        df_wes.rename(columns={'sample':'name'}, inplace=True)

        # Dataframe of WGS agCN estimates
        try:
            df_wgs = pd.read_csv(fp_wgs, sep='\t')
            df_wgs = df_wgs[['name', 'agCN']]
        except FileNotFoundError:
            print(f"Could not find the file: {fp_wgs}")
            continue

        try:
            df_wgs['change'] = df_wgs.agCN.apply(lambda v: len(set([c for c in v.split(',') if c!='.']))>1)
        except AttributeError: # agCN truth is only one exon long
            df_wgs['change'] = False

        # Combined dataframe
        result = pd.merge(df_wes, df_wgs, how='left', on='name')
        result.rename(columns={'agCN':'truth', 'hmm_CN_filt':'hmm'}, inplace=True)

        # If agCN truth is only one exon long
        # if pd.api.types.is_float_dtype(result.truth):
        #     result.truth = np.floor(pd.to_numeric(result.truth, errors='coerce')).astype('Int64')

        EXONLEN = len(str(result.hmm[0]).split(','))
        result['truth']  = result.truth.fillna(','.join(['.']*EXONLEN))
        result['diffs']  = result.apply(lambda row: get_diffs(row.hmm, row.truth), axis=1)
        result['change'] = result['change'].astype('boolean').fillna(False)
        result['is_ref'] = result.truth.apply(lambda x: is_ref(x, refCN))

        # Clean up
        cols = ['locus', 'name', 'comp', 'is_ref', 'change', 'truth', 'grdy_CN', 
                'hmm_CN', 'hmm', 'diffs', 'hmm_qual', 'refset_size', 'alpha_beta']
        result = result[cols]

        # Export to output file
        fp_eval = os.path.join(out_dir, f'{GENE}.eval')
        result.to_csv(fp_eval, sep='\t', index=None)
          

# ============================================================================
# Run evaluations
# ============================================================================
def parse_args():
    """ Commandline argument parser """
    
    # Parse arguments from command-line
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--wes-dir", required=True, 
                        help="path to Edgecopy WES HMM estimate directory")
    parser.add_argument("-g", "--wgs-dir", required=True, 
                        help="path to Parascopy WGS estimate directory")
    parser.add_argument("-o", "--out-dir", required=True, 
                        help="path to directory to store all output")
    parser.add_argument("-r", "--refcn", required=True, 
                        help="path to list of genes with refCN (tab-sep)")
    parser.add_argument("-a", "--allexons", required=True, 
                        help="path to 'allexons.parascopy.examine.output'")
    parser.add_argument("--qual-thresh", required=False, type=int, default=20, 
                        help="Phred quality threshold")
    parser.add_argument("--min-cc-size", required=False, type=int, default=5,
                        help="Minimum size of connected component")
    parser.add_argument("--rerun", action='store_true', required=False,
                        help="Rerun regardless of existing files")
    args_parsed = parser.parse_args()
    return args_parsed


def main():
    args = parse_args()

    wes_dir = args.wes_dir
    wgs_dir = args.wgs_dir
    out_dir = args.out_dir
    refcn_list_fp = args.refcn
    allexons_fp = args.allexons
    qual_thresh = args.qual_thresh
    min_cc_size = args.min_cc_size
    rerun = args.rerun
    
    # (1) cn_per_exon.py
    # - Create WGS files that are easier to compare
    # - Divide genes based on reference CN if any has multi-reference CNs.
    print("\nExamining and preprocessing WGS files.\n")
    cp.run(wgs_dir, out_dir, allexons_fp, rerun=rerun)
    
    # (2) Evaluate and compare WES vs. WGS
    print("\nEvaluating WES HMM estimates based on WGS estimates.\n")
    evaluate_edgecopy(wes_dir, wgs_dir, out_dir, refcn_list_fp, 
                      QUAL_THRESH=qual_thresh, MIN_CC_SIZE=min_cc_size)
    
    print("\nEvaluation completed.")


if __name__ == '__main__':
    main()
