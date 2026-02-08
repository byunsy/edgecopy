"""Utilities for building ExomeDepth-style reference sets.

This module analyzes count matrices (exon x sample) and selects reference
sample sets for each test sample using a beta-binomial model. The code is
adapted from ExomeDepth-style logic and focuses on estimating beta-binomial
parameters and selecting reference samples that maximize an expected Bayes
Factor criterion.

The edits in this file are documentation-only (module docstring, function
docstrings, and inline comments). No computational logic has been changed.
"""

import sys
import os
import math
import random
import argparse
import numpy as np
import pandas as pd
import scipy.stats as stats

from scipy.special import betaln, logsumexp, gammaln
from math import lgamma
from scipy.optimize import minimize
from concurrent.futures import ProcessPoolExecutor
from scipy.optimize import nnls
from statsmodels.nonparametric.smoothers_lowess import lowess


# Replicates ExomeDepth function get.power.betabinom. Returns expected Bayes
# Factor E_{Y ~ BB(n, p1, phi)} [ L1(y) / L0(y) ]. Used to evaluate
# how well a candidate reference set fits under an alternative vs null.
def expected_BF(median_depth, alpha, beta):
    """Compute expected Bayes Factor under a beta-binomial model.

    Parameters
        median_depth (int): typical depth (used as n in the PMF).
        alpha (float): beta prior alpha.
        beta (float): beta prior beta.

    Returns
        float: expected Bayes Factor comparing a shifted beta (p1)
               against the null beta (p0) at the given depth.
    """

    def beta_binomial_pmf_efficient(n, alpha, beta):
        # compute log PMF iteratively for stability and speed
        log_pmf = np.zeros(n + 1)
        log_pmf[0] = betaln(alpha, beta + n) - betaln(alpha, beta)
        for y in range(n):
            log_pmf[y + 1] = (
                log_pmf[y]
                + np.log(n - y)
                - np.log(y + 1)
                + np.log(alpha + y)
                - np.log(beta + n - 1 - y)
            )
        return log_pmf

    # construct alternative prior (shifted towards p=0.5)
    alpha1 = (alpha + beta) * (alpha * 0.5) / (alpha * 0.5 + beta)
    beta1 = alpha + beta - alpha1
    L0 = beta_binomial_pmf_efficient(int(median_depth), alpha, beta)
    L1 = beta_binomial_pmf_efficient(int(median_depth), alpha1, beta1)
    # expected BF approximated from PMF logs
    return np.sum(np.exp(L1) * (L1 - L0))


## betaln works for large a and b...
def betabin_likelihood(params, data1, data2, exons):
    """Negative log-likelihood for the beta-binomial given parameters.

    This function is an internal objective used for numerical MLE. It
    expects `data1` and `data2` to be arrays containing counts per exon and
    `exons` to be the list/iterable of exon indices to include.
    """

    a, b = params
    if a < 0 or b < 0:
        # penalize invalid parameter regions (avoid negative space)
        return 1e12
    c = betaln(a, b)
    return -sum([betaln(data1[i] + a, data2[i] + b) - c for i in exons])

def MLE_betabb(vec1, vec2, exons, binomial=False):
    """Compute MLE for beta-binomial alpha/beta by numerical minimization.

    Returns estimated (alpha, beta, objective_value). Uses Nelder-Mead as a
    robust initial optimizer; more constrained optimizers are commented.
    """

    alpha, beta = 10, 10
    result = minimize(betabin_likelihood, [alpha, beta], args=(vec1, vec2, exons), method='Nelder-Mead')
    # alternative: L-BFGS-B with bounds (commented for now)
    # result = minimize(betabin_likelihood, [10,10], args=(vec1,vec2,exons,), method='L-BFGS-B', bounds=[(1,1000000),(1,1000000)])
    return result.x[0], result.x[1], result.fun
  
## use branch and bound?
def find_best_refset(vec1, counts, candidates, exons, binomial=False, DEBUG=False):
    """Greedy search to build a reference set maximizing expected BF.

    The routine tests each candidate by temporarily adding it to the
    reference set, estimating beta-binomial parameters via method-of-moments,
    computing an expected Bayes Factor, and selecting the candidate that
    improves the global criterion by a small margin. This is repeated until
    no candidate yields a meaningful improvement.

    Returns
        refset (list): indices selected as reference samples.
        g_alpha (float): alpha estimate for final refset.
        g_sum (float): alpha+beta (dispersion proxy) for final refset.
    """

    def obfn(alpha, beta, f=2):
        # optional objective variants (kept for historical reasons)
        if f == 1:
            return alpha * beta / (alpha + beta + 0.001)
        elif f == 2:
            return alpha + beta
        else:
            return (alpha + beta) * np.sqrt(alpha)

    selected = set()
    remaining = set(candidates)
    refset = []

    g_alpha, g_beta, g_BF = 0, 0, 0

    while True:
        best_alpha, best_beta, best_s, best_BF = 0, 0, -1, 0
        values = []
        for s in remaining:
            # try adding s
            refset.append(s)
            vec2 = np.sum(counts[:, refset], axis=1)
            alpha, beta, mean, flag = method_moments_estimate(vec1, vec2, exons)
            BF = expected_BF(np.median(vec1 + vec2), alpha, beta)
            if BF > best_BF:
                best_alpha, best_beta, best_s, best_BF = alpha, beta, s, BF
                if flag and DEBUG:
                    print('binomial fit found', s)
            refset.pop()
            values.append((s, round(alpha), round(beta)))

        # accept candidate only if it improves BF by a small margin
        if best_BF > g_BF + 0.02:
            remaining.remove(best_s)
            selected.add(best_s)
            refset.append(best_s)
            g_alpha, g_beta, g_BF = best_alpha, best_beta, best_BF
            if DEBUG:
                print('adding', best_s, best_alpha, best_beta, best_BF)
        else:
            if DEBUG:
                print('cannot improve fit,exiting')
            break

    if DEBUG:
        print('FINAL', refset, g_alpha, g_alpha + g_beta, g_BF)

    return refset, g_alpha, g_beta + g_alpha


def method_moments_estimate(vec1, vec2, exons, binomial=False):
    """Method-of-moments estimate for beta-binomial alpha and beta.

    This provides a quick (approximate) estimate used as an initial
    diagnostic and for guiding greedy refset selection. The function
    computes sample mean and a dispersion-derived sum_ab; if the estimate
    is negative it is clamped and a flag is returned.

    Returns
        alpha (float), beta (float), mean (float), flag (bool)
    """

    # totals across selected exons
    nsum = sum([vec1[i] for i in exons])
    dsum = sum([vec1[i] + vec2[i] for i in exons])
    mean = nsum / dsum

    # compute empirical sums used to estimate dispersion
    diffs1 = sum([(vec1[i] - mean * (vec1[i] + vec2[i])) ** 2 for i in exons])
    diffs2 = sum([(1.0 - mean) * mean * (vec1[i] + vec2[i]) for i in exons])
    diffs3 = sum([(1.0 - mean) * mean * (vec1[i] + vec2[i]) * (vec1[i] + vec2[i]) for i in exons])
    phi_num = diffs3 - diffs1
    phi_den = diffs1 - diffs2
    sum_ab = phi_num / phi_den
    flag = False
    if sum_ab < 0:
        # clamp to a large dispersion value and signal via flag
        sum_ab = 1000000
        flag = True

    alpha, beta = mean * sum_ab, sum_ab - mean * sum_ab
    return alpha, beta, mean, flag

    # print(mean,diffs1,diffs2,diffs3,phi_num,phi_den,phi)
    # nsum=0;dsum=0
    # for i in exons:
    #     n = vec1[i]+vec2[i]
    #     diff = vec1[i]-n*mean
    #     expected = n*mean*(1.0-mean) ## np(1.0-p)
    #     nsum += diff*diff - expected
    #     dsum += n*(n-1)
    # dsum *= mean*(1.0-mean)
    # dsum /= nsum 
    # alpha,beta = dsum*mean, dsum -dsum*mean
    # return alpha,beta,mean

def process_sample(n_exons, sample_names, i, corr_list, counts, maxpairs, DEBUG=True):
    """Process a single sample to build its reference set.

    Steps:
      - Build a non-negative least-squares (NNLS) fit using top correlated
        candidate samples to approximate the test sample.
      - Rank candidates by NNLS weights and run `find_best_refset` to greedily
        select a subset that improves the expected BF.

    Returns
        prev_alpha, prev_sum, refset, top_correlation, logs
    """

    exons = [e for e in range(n_exons)]
    vec1 = counts[:, i]
    logs = []

    x = np.array(counts[:, i])
    xmean = np.mean(x)

    # compute per-candidate means to normalize NNLS predictors
    Ameans = [np.mean(counts[:, corr_list[k][0]]) for k in range(1, maxpairs)]
    A = np.array([counts[:, corr_list[k][0]] / Ameans[k - 1] for k in range(1, maxpairs)])

    # solve non-negative least squares for weights
    w, residual = nnls(A.T, x / xmean)
    wlist = [(w[k - 1], corr_list[k][0], sample_names[corr_list[k][0]], corr_list[k][1]) for k in range(1, maxpairs)]
    wlist.sort(reverse=True)
    for w in wlist:
        logs.append(w)
    logs.append(['NNLS-residual', residual, residual * xmean])

    # candidates for greedy refset selection
    candidates = [w[1] for w in wlist]
    refset, prev_alpha, prev_sum = find_best_refset(vec1, counts, candidates, exons, binomial=False)

    logs.append(['build reference set for sample', i, sample_names[i], corr_list[1], len(refset), prev_sum, 1.0 / prev_sum, prev_alpha, '\n-----------'])
    if len(refset) == 0:
        logs.append(['ERROR-empty refset'])

    return prev_alpha, prev_sum, refset, corr_list[1][1], logs


## df is counts data table for all exons from ExomeDepth
def referenceset_fit(count_matrix_file, interval_file, output_file, 
                     logfile=sys.stdout, exons_for_beta=10000, max_ref=30, threads=8):
    """Main driver: build reference sets for all samples in a count matrix.

    Steps performed:
      - Read interval and count matrix files.
      - Sample ~10K exons (excluding extreme length tails) to estimate
        per-exon scaling and compute a correlation matrix across samples.
      - For each sample, compute candidate correlated samples and build a
        reference set using `process_sample` (parallelized).
      - Write a summary file with estimated dispersion (phi), mu and chosen
        reference samples.

    Parameters mirror the original script: paths to files, threading and
    tuning parameters for the number of exons used and maximum ref set size.
    """

    intervals_df = pd.read_csv(interval_file, sep='\t')
    allexons_df = pd.read_csv(count_matrix_file, sep='\t')
    n_exons = allexons_df.shape[0]
    n_samples = allexons_df.shape[1]
    print('read count matrix dimensions', allexons_df.shape, file=sys.stdout)

    # remove 5% tails by exon length and randomly select up to exons_for_beta
    length_exons = [(int(intervals_df.iloc[i, 2]) - int(intervals_df.iloc[i, 1]), i) for i in range(n_exons)]
    length_exons.sort()
    five_percent = int(n_exons * 5 / 100)
    random.seed(42)
    if 0.9 * n_exons > exons_for_beta:
        exons_keep = random.sample([length_exons[i] for i in range(five_percent, n_exons - five_percent)], k=exons_for_beta)
    else:
        exons_keep = [length_exons[i] for i in range(five_percent, n_exons - five_percent)]

    exons_dict = {e[1]: 1 for e in exons_keep}
    scales = np.array([e[0] for e in exons_keep])

    # subset count matrix to selected exons and scale by exon length
    df1 = allexons_df[allexons_df.index.isin(exons_dict.keys())]
    idf1 = allexons_df[allexons_df.index.isin(exons_dict.keys())]
    df1_scaled = df1.div(scales, axis=0)
    corrmat = df1_scaled.corr()

    maxpairs = min(max_ref, n_samples)
    counts = np.array([df1[df1.columns[i]].to_numpy() for i in range(n_samples)]).T
    means = [df1[df1.columns[i]].mean(axis=0) for i in range(n_samples)]

    sample_names = df1.columns
    n_exons = df1.shape[0]
    lists = []
    for i in range(n_samples):
        corr_list = list(enumerate(corrmat.iloc[i, :].tolist()))
        corr_list.sort(reverse=True, key=lambda x: x[1])
        lists.append(corr_list)

    print('building reference set using', threads, 'threads', file=sys.stdout)

    # parallel processing of samples
    if threads == 1:
        results = []
        for i in range(n_samples):
            print('potential refset samples for', sample_names[i])
            result = process_sample(n_exons, sample_names, i, lists[i], counts, maxpairs)
            prev_alpha, prev_sum, refset, corr, logs = result
            results.append(result)
    else:
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(process_sample, n_exons, sample_names, i, lists[i], counts, maxpairs) for i in range(n_samples)]
            results = [futures[i].result() for i in range(n_samples)]

    # write output summary with estimated parameters and chosen references
    out_fp = open(output_file, 'w')
    print('sample', 'correlation', 'phi', 'mu', 'logit_mu', 'mean_read_counts', 'num_reference_set', 'reference_set', sep='\t', file=out_fp)
    for i in range(n_samples):
        prev_alpha, prev_sum, refset, corr, logs = results[i]
        mu = prev_alpha / prev_sum
        logit_mu = math.log(mu) - math.log(1.0 - mu)
        mean = means[i]
        print(sample_names[i], round(corr, 4), 1.0 / prev_sum, mu, round(logit_mu, 2), mean, len(refset), ','.join([sample_names[s] for s in refset]), file=out_fp, sep='\t')
        for msg in logs:
            print(' '.join([str(a) for a in msg]), file=logfile)
    out_fp.close()

    print('finished building reference sets', file=sys.stdout)
            

#####################################################################################################
## only take directory as input | file1 = all.counts.nondup.tsv all.counts.nondup.meta.tsv | all.stats.tsv
## first= count matrix, second = exon intervals, third = output file

def parse_args():
    
    # Parse arguments from command-line
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", required=True, 
                        help="directory with count files ('all' directory)", default=None)
    parser.add_argument("-o","--out", required=False, 
                        help="output file with beta-binomial params and refsets", default=None)
    # parser.add_argument("--method", required=False, 
    #                     help="method for beta-binomial reference fit", default='HM')
    parser.add_argument("--exons", required=False, 
                        help="number of exons for beta-binomial reference fit", default=10000)
    parser.add_argument("--maxsize", required=False, 
                        help="max size of reference set", default=30)
    parser.add_argument("-@","--threads", required=False, 
                        help="threads for parallelization", default=4)
    args_parsed = parser.parse_args()
    return args_parsed

if __name__ == "__main__":
    args = parse_args()
    directory = args.dir
    count_matrix_file = directory + '/all.counts.nondup.tsv'
    interval_file = directory + '/all.counts.nondup.meta.tsv'
    
    if args.out == None: 
        stats_file= directory + '/all.stats.tsv'
        log_file= open(stats_file + '.log','w')
    else: 
        stats_file = args.out 
        log_file= open(args.out + '.log','w')

    referenceset_fit(count_matrix_file, interval_file, stats_file,
                    logfile=log_file, threads=int(args.threads), max_ref=args.maxsize, exons_for_beta=args.exons)
    log_file.close()
