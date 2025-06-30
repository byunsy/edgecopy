import os
import sys
import pandas as pd
import numpy as np
from glob import glob
from scipy import stats

def concat_eval_files(eval_dir):
    """
    Concatenate all evaluation files in given evaluation directory
    """
    fps = sorted(glob(os.path.join(eval_dir, f'*.eval')))
    dfs = [pd.read_csv(fp, sep='\t') for fp in fps]
    concat_df = pd.concat(dfs, axis=0).reset_index(drop=True)
    return concat_df

def estim_vs_truth(e, t, refCN):
    """ 
    Input params: e (estim) and t (truth) are agCNs at a single position. 
    Returns TP, FP, FN, and TN (no event) based on e and t.
    """
    
    duptp, dupfp, dupfn = 0, 0, 0
    deltp, delfp, delfn = 0, 0, 0
    no_event = 0
    
    if t > refCN:         # Actual duplication
        if e == t:        # - predicted dup
            duptp += 1
        elif e > t:       # - predicted wrong dup
            dupfp += 1
        else:             # - predicted no dup
            dupfn += 1

    elif t < refCN:       # Actual deletion
        if e == t:        # - predicted del
            deltp += 1
        elif e < t:       # - predicted wrong del
            delfp += 1
        else:             # - predicted no del
            delfn += 1

    else:                 # Actual normalcy
        if e == t:        # - predicted normal
            no_event += 1
        elif e > t:       # - predicted dup
            dupfp += 1
        else:             # - predicted del
            delfp += 1

    return np.array([duptp, dupfp, dupfn, deltp, delfp, delfn, no_event], dtype='int32')

def v_tpfpfn(estim1, estim2, truth, refCN):
    """
    Input
    - estim1: HMM estimate without qual threshold
    - estim2: HMM estimate with q>20 threshold
    - truth : ground truth
    
    Compute the number of TPs, FPs, and FNs of estimated agCNs
    - separately for duplication and deletion events
    - also compute TP,TN,FP, and FN for NoCalls
    
    Returns array that sums to one:
    np.array([
        dup_tp, dup_fp, dup_fn, del_tp, del_fp, del_fn, no_events, 
        nc_duptp, nc_dupfp, nc_dupfn, nc_deltp, nc_delfp, nc_delfn, nc_no_event,
    ]) / total_num_exons
    """
    
    estim1 = str(estim1).split(',')
    estim2 = str(estim2).split(',')
    truth = str(truth).split(',')
    
    if len(estim1) != len(estim2) or len(estim1) != len(truth):
        raise ValueError("Input arrays must have the same length. "
                         f"e1={len(estim1)}, e2={len(estim2)}, t={len(truth)}.")
    
    tpfpfn = np.zeros(7, dtype='int32')
    nc_tpfpfn = np.zeros(7, dtype='int32')
    
    EXONLEN = len(estim1)
    num_notruth = 0
    
    for i in range(EXONLEN):
        
        if truth[i] == '.':
            num_notruth += 1
            continue
        
        if estim2[i] == '_':
            e, t = int(estim1[i]), int(truth[i])
            nc_tpfpfn += estim_vs_truth(e, t, refCN)
        else:
            e, t = int(estim2[i]), int(truth[i])
            tpfpfn += estim_vs_truth(e, t, refCN)

    EXONLEN_ANALYZED = EXONLEN - num_notruth
    if not EXONLEN_ANALYZED:
        return np.array([0]*6 + [1] + [0]*7)

    ret = np.concatenate((tpfpfn, nc_tpfpfn))

    # Result should sum to one
    ret = np.array(ret)/(EXONLEN_ANALYZED)
    assert np.isclose(np.sum(ret), 1, atol=1e-8)
    
    return ret.round(3)


def compute_score(rare_df, eval_dir):
    """
    Compute precision and recall values for genes with rare CNVs
    """
    save_df = []
    for unique_refcn in rare_df.refcn.unique():
        print('='*50)
        print('Reference copy number:', unique_refcn)
        print('='*50)
    
        df = rare_df.loc[rare_df.refcn == unique_refcn].copy()
        df['refcn'] = df.refcn.astype('Int32')

        # Compute TP, FP, and FNs, including NoCall values
        df['tpfpfn'] = df.apply(lambda row: v_tpfpfn(row.hmm_CN, row.hmm, row.truth, row.refcn), axis=1)
        df_tpfpfn = pd.DataFrame(df.groupby('locus').tpfpfn.sum()).reset_index()

        # Compute precision and recall
        df_tpfpfn['dup_prec'] = df_tpfpfn.tpfpfn.apply(lambda x: x[0] / (x[0] + x[1]))
        df_tpfpfn['dup_reca'] = df_tpfpfn.tpfpfn.apply(lambda x: x[0] / (x[0] + x[2] + x[7]))
        df_tpfpfn['del_prec'] = df_tpfpfn.tpfpfn.apply(lambda x: x[3] / (x[3] + x[4]))
        df_tpfpfn['del_reca'] = df_tpfpfn.tpfpfn.apply(lambda x: x[3] / (x[3] + x[5] + x[10]))
        df_tpfpfn['prec'] = df_tpfpfn.tpfpfn.apply(lambda x: (x[0]+x[3]) / (x[0]+x[3]+x[1]+x[4]))
        df_tpfpfn['reca'] = df_tpfpfn.tpfpfn.apply(lambda x: (x[0]+x[3]) / (x[0]+x[3]+x[2]+x[7]+x[5]+x[10]))

        # info = [str(s) for s in df_tpfpfn.tpfpfn.sum().round(2)]
        # print('Called Dup:', ', '.join(info[:3]))
        # print('Called Del:', ', '.join(info[3:6]))
        # print('Called No Event:', info[6], '\n')
        # print('NoCall Dup:', ', '.join(info[7:10]))
        # print('NoCall Del:', ', '.join(info[10:13]))
        # print('NoCall No Event:', info[13], '\n')

        # df_tpfpfn['nc_perc'] = df_tpfpfn.tpfpfn.apply(lambda x: sum(x[7:])/sum(x[:]))

        # total_nc = np.sum(np.sum(df_tpfpfn.tpfpfn)[7:])
        # total_ = np.sum(np.sum(df_tpfpfn.tpfpfn))
        # print('NoCall%:', round(total_nc/total_, 3), '\n')

        # print(*round(df_tpfpfn[['dup_prec','dup_reca','del_prec','del_reca','nc_perc']].mean(), 4).to_list(), sep='\n')
        # print()
        
        print(*round(df_tpfpfn[['prec','reca']].mean(), 4).to_list(), sep='\n')
        print()

        # sem_prec = round(stats.sem(df_tpfpfn['prec'], nan_policy='omit'),4)
        # sem_reca = round(stats.sem(df_tpfpfn['reca'], nan_policy='omit'),4)
        # print(sem_prec, sem_reca, sep='\n')

        df_tpfpfn['refcn'] = unique_refcn
        df_tpfpfn['method'] = 'Edgecopy'
        print(df_tpfpfn[['locus','prec','reca','refcn','method']])
        
        save_df.append(df_tpfpfn[['locus','prec','reca','refcn','method']])

    outfp = os.path.join(eval_dir, 'rare_cnvs.precision.recall.tsv')
    to_save = pd.concat(save_df, axis=0).reset_index(drop=True)
    to_save.to_csv(outfp, sep='\t', na_rep='NaN', index=None)
    
    return outfp

def run(eval_dir, refcn_fp):
    """
    eval_dir : Directory generated by running evaluate.py, containing all eval files
    refcn_fp : A TSV file with gene, refcn, and position (optional)
    """
    
    # Get reference copy numbers for genes to analyze
    refcns = pd.read_csv(refcn_fp, sep='\t', header=None, names=['gene', 'refcn', 'pos'])
    refcns = refcns[['gene', 'refcn']]
    refcns = dict(zip(refcns['gene'], refcns['refcn']))

    # Get evaluation files of genes to analyze
    concat_df = concat_eval_files(eval_dir)
    concat_df['refcn'] = concat_df.apply(lambda row: refcns.get(row.locus), axis=1)

    rare_cnvs = concat_df.locus.unique()
    rare_df = concat_df.loc[concat_df.alpha_beta > 1000].copy()

    print(f'Number of genes with rare CNVs in {eval_dir}:', len(rare_cnvs))
    print("DataFrame shape:", rare_df.shape)
    print("- Analyzing", rare_df['name'].unique().shape[0], "unique samples.")
    
    saved_fp = compute_score(rare_df, eval_dir)
    print(f"\nCompleted. Saved to {saved_fp}")
    
    
if __name__ == '__main__':
    eval_dir, refcn_fp = sys.argv[1], sys.argv[2]
    run(eval_dir, refcn_fp)