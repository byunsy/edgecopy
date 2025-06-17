import sys
import math
import numpy as np
from scipy.special import betaln

#SIGMA=0.35

def fopt_prior_gradient(x,data,priors=None,epsilon=1e-8,sigma=0.35):
    """
    Same function will work for fopt and fopt_prior
    Return list or vector
    """
    n = len(priors)

    def prior_ll(cn):
        if priors==None: return 0.0
        penalty=0.0
        for j in range(n):
            delta = (cn-priors[j][0])*(cn-priors[j][0])/(2*sigma*sigma)
            penalty += math.exp(-delta)*priors[j][1]    
        return -math.log(penalty)
        
    derivates = [0.0 for i in range(data.n)]
    for i in range(data.n):
        f1 = fopt(x,data,sample_list=data.samples_include[i]) + prior_ll(x[i])
        x[i] += epsilon
        f2 = fopt(x,data,sample_list=data.samples_include[i]) + prior_ll(x[i])
        derivates[i] = (f2-f1)/(epsilon)
        x[i] -= epsilon
    
    return derivates


def fopt_prior(x, data, priors, sigma=0.35):
    """
    Fractional CN prior, mixture of normal distributions centred at integer 
    values with prior weights
    - priors is a tuplelist
    """
    LL = fopt(x, data)
    n = len(priors)
    
    for i in range(data.n):
        penalty = 0.0
        for j in range(n):
            delta = (x[i]-priors[j][0])*(x[i]-priors[j][0])/(2*sigma*sigma)
            penalty += math.exp(-delta)*priors[j][1]
        LL -= math.log(penalty)
    
    return LL


def fopt(x, data, sample_list=None, scale_factor=1):
    """
    Beta-binomial log-likelihood function
    x : copy-number vector
    """
    n_exons = len(data.ecounts[0])
    if sample_list ==None: sample_list = [i for i in range(data.n)]
    #sample_list = [i for i in range(data.n)]
    
    # For each sample, compute beta log-likelihood
    total_ll = 0.0
    for i in sample_list:
        
        # Obtain alpha and beta from ExomeDepth
        alpha, beta = data.alpha[i], data.beta[i]
        
        # Given alpha and beta, compute new alpha and beta
        a1 = alpha*x[i]
        b1_denom = sum(data.means[j] for j in data.reference_sets[i])
        b1 = beta*sum([x[j]*data.means[j]/b1_denom for j in data.reference_sets[i]])
        
        # Ensure a1 and b1 is non-neg (e.g. when x[i]==0)
        # - if a1==0, then alpha_new==0, which causes error
        a1, b1 = max(a1, 1e-6), max(b1, 1e-6)
        
        try:
            alpha_new = a1*(alpha+beta)/(a1+b1)
            beta_new = b1*(alpha+beta)/(a1+b1)
        except ZeroDivisionError:
            print("ERROR", data.genename)
            alpha_new = 1e-6
            beta_new = 1e-6
        
        # Compute beta log-likelihoods (sum beta_ll across exons) 
        beta_ll = 0.0
        for e in range(n_exons):
            c = data.ecounts[i][e]
            Rc = sum([data.ecounts[j][e] for j in data.reference_sets[i]])
            beta_ll += betaln(c+alpha_new, Rc+beta_new) - betaln(alpha_new, beta_new)
        total_ll += beta_ll+data.comb_log[i]
        
    return -1*total_ll


def fopt_exon(hmm_data, cn_mat, scale_factor=1):
    """
    Beta-binomial log-likelihood function
    cn_mat : copy-number matrix (sample x exons)
    """
    total_ll = 0.0 
    for s in range(hmm_data.num_samples):
        beta_ll = 0.0
        for obs in range(hmm_data.num_obsrvs):
    
            # Obtain alpha and beta
            alpha, beta = hmm_data.alpha[s], hmm_data.beta[s]

            # Given alpha and beta, compute new alpha and beta
            a1 = alpha*cn_mat[s,obs]
            b1_denom = sum(hmm_data.means[j] for j in hmm_data.refset[s])
            b1 = beta*sum([cn_mat[j,obs]*hmm_data.means[j]/b1_denom for j in hmm_data.refset[s]])
            
            # Ensure a1 and b1 is non-neg (e.g. when x[i]==0)
            # - if a1==0, then alpha_new==0, which causes error
            a1, b1 = max(a1, 1e-6), max(b1, 1e-6)
            
            try:
                alpha_new = a1*(alpha+beta)/(a1+b1)
                beta_new = b1*(alpha+beta)/(a1+b1)
            except ZeroDivisionError:
                print("Error", data.genename)
                alpha_new = 1e-6
                beta_new = 1e-6

            # Compute beta log-likelihood 
            c = hmm_data.obsrvs[s,obs]
            Rc = hmm_data.refsum[s,obs]
            beta_ll += betaln(c+alpha_new, Rc+beta_new) - betaln(alpha_new, beta_new)# + comb_log(c,Rc)
        
        total_ll += beta_ll + hmm_data.comb_log[s]
    #print('hmm_data-scale',hmm_data.scale_factor,file=sys.stderr)
    
    return -1*total_ll

def pairwise_graph(data, thresh=5): ## calculate likelihoods for each edge
    return 1
