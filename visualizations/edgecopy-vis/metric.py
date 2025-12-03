import pandas as pd
import numpy as np


def v_tpfpfn(estim, truth, refCN):
    """
    Compute the number of TPs, FPs, and FNs of estimated agCNs
    - separately for duplication and deletion events
    """
    estim, truth = str(estim).split(','), str(truth).split(',')
    EXONLEN = len(estim)
    
    duptp, dupfp, dupfn = 0,0,0
    deltp, delfp, delfn = 0,0,0
    num_noevents = 0
    num_filtered = 0
    
    for i in range(EXONLEN):
        
        if truth[i] == '.' or estim[i]=='0':
            num_filtered += 1
            continue
        
        if estim[i] == '_':
            num_filtered += 1
            continue

        e, t = int(estim[i]), int(truth[i])

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
                num_noevents += 1
                continue
            elif e > t:       # - predicted dup
                dupfp += 1
            else:             # - predicted del
                delfp += 1

    tpfpfn = [duptp, dupfp, dupfn, deltp, delfp, delfn, num_noevents]
    
    if num_filtered == EXONLEN:
        return np.array([0]*6 + [1])

    return (np.array(tpfpfn)/(EXONLEN-num_filtered)).round(2)


def v_tpfpfn2(estim, truth, refCN, gt_miss, ec_miss, ab_miss):
    """
    Compute the number of TPs, FPs, and FNs of estimated agCNs
    - separately for duplication and deletion events
    """
    
    if gt_miss or ec_miss or ab_miss:
        return np.array([0]*6 + [1])
    
    estim, truth = str(estim).split(','), str(truth).split(',')
    EXONLEN = len(estim)
    
    duptp, dupfp, dupfn = 0,0,0
    deltp, delfp, delfn = 0,0,0
    num_noevents = 0
    num_filtered = 0
    
    for i in range(EXONLEN):
        
        if truth[i] == '.' or estim[i]=='0':
            num_filtered += 1
            continue
        
        if estim[i] == '_':
            num_filtered += 1
            continue

        e, t = int(estim[i]), int(truth[i])

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
                num_noevents += 1
                continue
            elif e > t:       # - predicted dup
                dupfp += 1
            else:             # - predicted del
                delfp += 1

    tpfpfn = [duptp, dupfp, dupfn, deltp, delfp, delfn, num_noevents]
    
    if num_filtered == EXONLEN:
        return np.array([0]*6 + [1])

    return (np.array(tpfpfn)/(EXONLEN-num_filtered)).round(2)


def estim_vs_truth(e, t, refCN):
    """ 
    INPUT: e (estim) and t (truth) are agCN at single position. 
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

def v_tpfpfn3(estim1, estim2, truth, refCN):
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


def v_correct(estim, truth, refCN):
    """
    Compute the number of correctly estimated agCNs
    - separately for reference vs non-reference CN
    """

    estim, truth = str(estim).split(','), str(truth).split(',')
    EXONLEN = len(estim)
    
    is_ref = all(t==str(refCN) for t in truth if t!='.')
        
    total, correct = 0, 0
    num_skipped = 0
    
    for i in range(EXONLEN):
        
        # If truth unavailable, skip
        if truth[i] == '.' or estim[i]=='_':
            num_skipped += 1
            continue
            
        if truth[i]:
            total += 1
            if estim[i] == truth[i]:
                correct += 1

    if num_skipped == EXONLEN:
        return np.array([np.nan, np.nan])
    
    if is_ref:            
        return np.array([correct/total, np.nan])
    
    return np.array([np.nan, correct/total])


def v_correct2(estim, truth, refCN, gt_miss, ec_miss, ab_miss):
    """
    Compute the number of correctly estimated agCNs
    - separately for reference vs non-reference CN
    """

    if gt_miss or ec_miss or ab_miss:
        return np.array([np.nan, np.nan])
    
    estim, truth = str(estim).split(','), str(truth).split(',')
    EXONLEN = len(estim)
    
    is_ref = all(t==str(refCN) for t in truth if t!='.')
        
    total, correct = 0, 0
    num_skipped = 0
    
    for i in range(EXONLEN):
        
        # If truth unavailable, skip
        if truth[i] == '.' or estim[i]=='_':
            num_skipped += 1
            continue
            
        if truth[i]:
            total += 1
            if estim[i] == truth[i]:
                correct += 1

    if num_skipped == EXONLEN:
        return np.array([np.nan, np.nan])
    
    if is_ref:            
        return np.array([correct/total, np.nan])
    
    return np.array([np.nan, correct/total])


def v_correct3(estim, truth, refCN):
    """
    Input
    - estim1: HMM estimate without qual threshold
    - estim2: HMM estimate with q>20 threshold
    - truth : ground truth
    
    Compute the number of correctly estimated agCNs
    - separately for reference vs non-reference CN
    """

    estim = str(estim).split(',')
    truth = str(truth).split(',')
    
    if len(estim) != len(truth):
        raise ValueError("Input arrays must have the same length. "
                         f"e={len(estim)}, t={len(truth)}.")
        
    EXONLEN = len(estim)
    is_ref = all(t==str(refCN) for t in truth if t!='.')
        
    total, corrc, error, filtr = 0, 0, 0, 0
    num_notruth = 0
    
    for i in range(EXONLEN):
        
        # If ground truth unavailable
        if truth[i] == '.':
            num_notruth += 1
            continue
        
        total += 1
        if estim[i]=='_':
            filtr += 1
        else:
            if estim[i] == truth[i]:
                corrc += 1
            else:
                error += 1

    EXONLEN_ANALYZED = EXONLEN - num_notruth
    if not EXONLEN_ANALYZED:
        return np.array([np.nan]*10)
    
    acc = corrc/total
    err = error/total
    fil = filtr/total
    assert np.isclose(acc+err+fil, 1, atol=1e-8)
    
    if not total-filtr:
        acc_f, err_f = np.nan, np.nan
    else:
        acc_f = corrc/(total-filtr)
        err_f = error/(total-filtr)
        assert np.isclose(acc_f+err_f, 1, atol=1e-8)
    
    if is_ref:
        return np.array([acc, err, fil, acc_f, err_f]+[np.nan]*5)
    return np.array([np.nan]*5+[acc, err, fil, acc_f, err_f])


def v_correct_combined(estim, truth, refCN):
    """
    Input
    - estim1: HMM estimate without qual threshold
    - estim2: HMM estimate with q>20 threshold
    - truth : ground truth
    
    Compute the number of correctly estimated agCNs
    - combined for reference vs non-reference CN
    """

    estim = str(estim).split(',')
    truth = str(truth).split(',')
    
    if len(estim) != len(truth):
        raise ValueError("Input arrays must have the same length. "
                         f"e={len(estim)}, t={len(truth)}.")
        
    EXONLEN = len(estim)
    is_ref = all(t==str(refCN) for t in truth if t!='.')
        
    total, corrc, error, filtr = 0, 0, 0, 0
    num_notruth = 0
    
    for i in range(EXONLEN):
        
        # If ground truth unavailable
        if truth[i] == '.':
            num_notruth += 1
            continue
        
        total += 1
        if estim[i]=='_':
            filtr += 1
        else:
            if estim[i] == truth[i]:
                corrc += 1
            else:
                error += 1

    EXONLEN_ANALYZED = EXONLEN - num_notruth
    if not EXONLEN_ANALYZED:
        return np.array([np.nan]*5)
    
    acc = corrc/total
    err = error/total
    fil = filtr/total
    assert np.isclose(acc+err+fil, 1, atol=1e-8)
    
    if not total-filtr:
        acc_f, err_f = np.nan, np.nan
    else:
        acc_f = corrc/(total-filtr)
        err_f = error/(total-filtr)
        assert np.isclose(acc_f+err_f, 1, atol=1e-8)
    
    return np.array([acc, err, fil, acc_f, err_f])


def examine_cnvs(df, cnv_type):

    # 1. Drop any NaN values
    df.dropna(inplace=True)

    # 2. Select genes with rare/common CNVs
    df = df.loc[df.gene.apply(lambda x: x in cnv_type)]
    
    # 3. Remove any ground-truth where more than half of them are dots
    # df = df.loc[~df.truth.apply(lambda x: sum(c=='.' for c in str(x).split(',')) >= len(str(x).split(','))/2)]
    
    # 4. Remove samples with all Falses for trust
    # df = df.loc[~df.trust.apply(lambda trst: all(t=='False' for t in str(trst).split(',')))]
    
    # 5. Select samples that has incorrect calls (X in diffs)
    df2 = df.loc[df.diffs.apply(lambda x: 'X' in x.split(','))]
    return df.sort_values(by='gene'), df2.sort_values(by='gene')


def examine_cnvs_(df, cnv_type):

    # 1. Drop any NaN values
    df.dropna(inplace=True)

    # 2. Select genes with rare/common CNVs
    df = df.loc[df.locus.apply(lambda x: x in cnv_type)]
    
    # 3. Remove samples with all '_'
    df = df.loc[~df.hmm.apply(lambda hmm: all(h=='_' for h in str(hmm).split(',')))]
    
    return df.sort_values(by='locus')

def examine_cnvs2(df, cnv_type):

    # 1. Drop any NaN values
    df.dropna(inplace=True)

    # 2. Select genes with rare/common CNVs
    df = df.loc[df.gene.apply(lambda x: x in cnv_type)].copy()
    
    # df['total'] = 1

    # 3. Remove any ground-truth where more than half of them are dots
    # df['gt_miss'] = df.truth.apply(lambda x: sum(c=='.' for c in str(x).split(',')) >= len(str(x).split(','))/2)
    df['gt_miss'] = df.apply(lambda r: sum(c=='.' for c in str(r.truth).split(',')) >= len(str(r.truth).split(','))/2 and not r.ab_filt, axis=1)
    
    # 4. Remove samples with all Falses for trust
    # df['ec_miss'] = df.trust.apply(lambda trst: all(t=='False' for t in str(trst).split(',')))
    df['ec_miss'] = df.apply(lambda r: all(t=='False' for t in str(r.trust).split(',')) and not r.gt_miss and not r.ab_filt, axis=1)

    return df.sort_values(by='gene')


def merge_alpha_beta_info(pop, center, df_fpfn):
    
    df_ab = pd.read_csv(f'alphabeta/all.alpha_beta.{pop}.{center}.tsv',sep='\t')
    df_ab['name'] = df_ab['sample'].apply(lambda x: x.split('.')[0])
    df_ab = df_ab[['name','alpha_beta']].copy()
    return pd.merge(df_fpfn, df_ab, on='name')


def compare_pscn(agcn1, agcn2, pscn1, pscn2, 
                 qual1, qual2, refCN, qual_thresh=20):
    """
    Compares two list of tuples containing psCNs
    Return reference and non-reference accuracy
    """
    agcn1 = agcn1.split(',')
    agcn2 = agcn2.split(',')
    pscn1 = [p.strip(')(') for p in pscn1.split(')(')]
    pscn2 = [p.strip(')(') for p in pscn2.split(')(')]
    qual1 = [q.strip(')(').split(',')[0] for q in qual1.split(')(')]
    qual2 = [q.strip(')(').split(',')[0] for q in qual2.split(')(')]
    
    e = 'lengths of psCN of ground truth and estimate differ.'
    assert len(agcn1) == len(agcn2), e
    assert len(pscn1) == len(pscn2), e

    EXONLEN = len(pscn1)
    
    # check ground truth to determine ref/non-ref
    refCN = int(refCN)
    is_ref = all(t==str(refCN) for t in agcn2 if t!='_')
    
    total, correct = 0, 0
    num_skipped = 0
    
    for i in range(EXONLEN):
        
        # if agcn is different
        if agcn1[i] != agcn2[i]:
            num_skipped += 1
            continue
        
        # if psCN is unreliable
        if '?' in pscn1[i] or '?' in pscn2[i]:
            num_skipped += 1
            continue

        # if quality is not high enough
        if int(qual1[i]) < qual_thresh or int(qual2[i]) < qual_thresh:
            num_skipped += 1
            continue

        total += 1
        if pscn1[i] == pscn2[i]:
            correct += 1

    if num_skipped == EXONLEN:
        return np.array([np.nan, np.nan])
            
    if is_ref:            
        return np.array([correct/total, np.nan])
    
    return np.array([np.nan, correct/total])


def compare_pscn2(agcn1, agcn2, pscn1, pscn2, 
                  qual1, qual2, qual_thresh=20):
    """
    Compares two list of tuples containing psCNs
    Checks if agCN and psCN perfectly matches
    Returns a list of three boolean values
    """
    agcn1 = agcn1.split(',')
    agcn2 = agcn2.split(',')
    pscn1 = [p.strip(')(') for p in pscn1.split(')(')]
    pscn2 = [p.strip(')(') for p in pscn2.split(')(')]
    qual1 = [q.strip(')(').split(',')[0] for q in qual1.split(')(')]
    qual2 = [q.strip(')(').split(',')[0] for q in qual2.split(')(')]
    
    e = 'lengths of psCN of ground truth and estimate differ.'
    assert len(agcn1) == len(agcn2), e
    assert len(pscn1) == len(pscn2), e

    EXONLEN = len(pscn1)
    
    valid, agcn_correct, pscn_correct = 1, 1, 1
    num_skipped = 0
    
    for i in range(EXONLEN):
        
        if '_' in agcn1[i] or '_' in agcn2[i]:
            num_skipped += 1
            continue
        
        # if agcn is different
        if agcn1[i] != agcn2[i]:
            agcn_correct = 0
            pscn_correct = 0
            continue
        
        # if psCN is unreliable
        if '?' in pscn1[i] or '?' in pscn2[i]:
            num_skipped += 1
            continue

        # if quality is not high enough
        if int(qual1[i]) < qual_thresh or int(qual2[i]) < qual_thresh:
            num_skipped += 1
            continue

        if pscn1[i] != pscn2[i]:
            pscn_correct = 0

    if num_skipped == EXONLEN:
        return np.array([0, 0, 0])
            
    return np.array([valid, agcn_correct, pscn_correct])