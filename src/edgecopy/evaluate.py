import os
import pandas as pd
import numpy as np
import argparse
from glob import glob
from collections import Counter


def get_diffs(estim, truth, trust):
    """
    Compare estimated CN vector against ground-truth and show
    at which positions the two differ
    - Mark '.' if estim==truth or if no ground-truth, or if no estim (no reads)
               or if estim[i] is not trustworthy
    - Mark 'X' if estim!=truth
    """
    if pd.isna(truth):
        return np.nan   
    try:
        estim = np.array([int(e) if e!='_' else 0 for e in estim.split(',')])
        truth = np.array([int(t) if t!='.' else 0 for t in truth.split(',')])
    except AttributeError:
        estim = np.array([int(estim)])
        if truth == '.':
            truth = np.array([0])
        else:
            truth = np.array([int(truth)])

    trust = [t=='True' for t in str(trust).split(',')]
    diffs = []

    for i,d in enumerate(estim-truth):
        if not trust[i]:
            diffs.append('.')
        elif d==0 or truth[i]==0:
            diffs.append('.')
        else:
            diffs.append('X')
    return ','.join(diffs)


def get_diffs2(estim, truth):
    """
    Compare estimated CN vector against ground-truth and show
    at which positions the two differ
    """
    if pd.isna(truth):
        return np.nan   
    try:
        estim = np.array([int(e) if e!='.' else 0 for e in estim.split(',')])
        truth = np.array([int(t) if t!='.' else 0 for t in truth.split(',')])
    except AttributeError:
        estim = np.array([int(estim)])
        if truth == '.':
            truth = np.array([0])
        else:
            truth = np.array([int(truth)])
    
    diffs = []
    for i,d in enumerate(estim-truth):
        if d==0 or truth[i]==0:
            diffs.append('.')
        else:
            diffs.append('X')
    return ','.join(diffs)


def is_ref(truth, refCN):
    if all(t=='.' for t in str(truth).split(',')):
        return False
    return Counter([t for t in str(truth).split(',') if t!='.']).most_common(1)[0][0]==str(refCN)


def evaluate_exome_paras(hmm_dir, truth_dir, refcn_list_fp, QUAL_THRESH=20, MIN_CC_SIZE=5):
    """
    Evaluate and summarize the point and HMM estimates using
    the WGS ground-truth from Parascopy
    """
    GENES = glob(f'{hmm_dir}/*.hmm.bed')
    GENES = [g.split('/')[-1].split('.hmm')[0] for g in GENES]
    
    REFCN = pd.read_csv(refcn_list_fp, sep='\t', header=None)
    REFCN = REFCN.set_index(0)[1].to_dict()
    
    EXCLD = os.path.join(hmm_dir, 'samples.exclude.list')
    samples_exclude = []
    if os.path.exists(EXCLD):
        with open(EXCLD, 'r') as inpf:
            samples_exclude = inpf.read().splitlines()
    
    for GENE in GENES:
        fp_hmm = f'{hmm_dir}/{GENE}.hmm.bed'
        fp_true = f'{truth_dir}/{GENE}.exon.CN.bed'
        refCN = REFCN.get(GENE)

        # Dataframe of output agCNs
        df_hmm = pd.read_csv(fp_hmm, sep='\t')
        if samples_exclude:
            df_hmm = df_hmm.loc[df_hmm['name'].apply(lambda x: x not in samples_exclude)]

        # Dataframe of ground-truth agCNs
        try:
            df_true = pd.read_csv(fp_true, sep='\t')
            df_true = df_true[['name', 'agCN']]
        except FileNotFoundError:
            print(f"Could not find the file: {fp_true}")
            continue

        try:
            df_true['change'] = df_true.agCN.apply(lambda v: len(set([c for c in v.split(',') if c!='.']))>1)
        except AttributeError: # agCN truth is only one exon long
            df_true['change'] = False

        # Combined dataframe
        result = pd.merge(df_hmm, df_true, how='left', on='name')
        result.rename(columns={'agCN':'truth', 'estim':'hmm'}, inplace=True)

        # If agCN truth is only one exon long
        if pd.api.types.is_float_dtype(result.truth):
            result.truth = np.floor(pd.to_numeric(result.truth, errors='coerce')).astype('Int64')

        #result['diffs'] = result.apply(lambda row: get_diffs(row.hmm, row.truth), axis=1)
        result['diffs'] = result.apply(lambda row: get_diffs(row.hmm, row.truth, row.trust), axis=1)
        
        result.change.fillna(False, inplace=True)
        result.fillna(pd.NA, inplace=True)
        try:
            EXONLEN = len(str(result.hmm[0]).split(','))
        except KeyError:
            print(GENE)
            print(result)
            continue

        #result['is_ref'] = result.truth.apply(lambda x: Counter(str(x).split(',')).most_common(1)[0][0]==str(refCN))
        result['is_ref'] = result.truth.apply(lambda x: is_ref(x, refCN))
        #result['trust'] = result.apply(lambda row: (not pd.isna(row.truth)) & (row.qual >= QUAL_THRESH) & (row.trust), axis=1)

        # Clean up
        cols = ['gene', 'name', 'comp', 'is_ref', 'change', 'truth', 'point', 'hmm', 'diffs', 'qual', 'refset_size','trust']
        result = result[cols]

        # Export to output file
        fp_eval = f'{hmm_dir}/{GENE}.eval'
        #result.to_csv(fp_eval, sep='\t', index=None, na_rep=' '*(2*EXONLEN-1))
        result.to_csv(fp_eval, sep='\t', index=None, na_rep=','.join(['.']*EXONLEN))
        #result.to_csv(fp_eval, sep='\t', index=None, na_rep='NA')
          

# ============================================================================
# Run evaluations
# ============================================================================
def parse_args():
    """ Commandline argument parser """
    
    # Parse arguments from command-line
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--hmm-dir", required=True, 
                        help="path to HMM output directory")
    parser.add_argument("-t", "--truth-dir", required=True, 
                        help="path to ground-truth directory")
    parser.add_argument("-r", "--refcn-fp", required=True, 
                        help="path to list of genes with refCN (tab-sep)")
    args_parsed = parser.parse_args()
    return args_parsed

def main():
    args = parse_args()

    hmm_dir = args.hmm_dir
    truth_dir = args.truth_dir
    refcn_list_fp = args.refcn_fp

    # Evaluate and compare estimated CNs against ground-truth
    # GENE.hmm.bed --> GENE.eval
    print("Evaluating point and HMM estimates.")
    evaluate_exome_paras(hmm_dir, truth_dir, refcn_list_fp)
    
    print("\nCompleted.")


if __name__ == '__main__':
    main()
