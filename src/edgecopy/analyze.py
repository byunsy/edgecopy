import os
import sys
import math
import pickle
import random
import argparse
import numpy as np
import pandas as pd
import scipy.optimize as opt
import networkx as nx

from . import cnhmm_run as _hmm
from . import exomecounts
from .exomecounts import *
from .optimize_functions import *

from time import perf_counter
from datetime import timedelta
from multiprocessing import Pool
from scipy.special import logsumexp
from bisect import bisect_right
from collections import Counter


def find_closest_vector(x, possible_vals=[2,3,4,5], prior={}): 
    """
    prior is a dictionary
    """
    n = len(x)
    M = len(possible_vals)
    dp = np.full((n + 1, M), np.inf)
    backtrack = np.zeros((n, M), dtype=int)

    dp[0, :] = 0  # base case

    for i in range(1, n + 1):
        for j, v in enumerate(possible_vals):
            best = np.inf
            best_k = -1
            for k in range(j + 1):  # u <= v to maintain non-decreasing
                u = possible_vals[k]
                if prior == {}: 
                    cost = dp[i-1][k] + abs(x[i-1]-v)
                else: 
                    cost = dp[i-1][k] + abs(x[i-1]-v) - 0.1*np.log(prior[v]/prior[x[i-1]])
                if cost < best:
                    best = cost
                    best_k = k
            dp[i][j] = best
            backtrack[i-1][j] = best_k

    # Reconstruct x'
    idx = np.argmin(dp[n])
    x_prime = [0]*n

    for i in range(n-1, -1, -1):
        x_prime[i] = possible_vals[idx]
        idx = backtrack[i][idx]

    #diffs = sum([1 for i in range(n) if x[i] != x_prime[i]])
    return x_prime

def closest_vector_full(frac_vec, integer_vec, distance=False, prior={}):
    """
    """
    n=len(frac_vec)
    seq = [[integer_vec[i],float(frac_vec[i]),i] for i in range(n)]
    seq.sort(key=lambda x:x[1])

    Xorig = [s[0] for s in seq]
    possible_vals = list(set(Xorig))
    possible_vals.sort()

    Xnew = find_closest_vector(Xorig, possible_vals=possible_vals, prior=prior)
    diffs = [i for i in range(n) if Xorig[i] != Xnew[i]]
    
    if distance: 
        return len(diffs)

    newvec = [(Xnew[i],seq[i][2]) for i in range(n)]
    newvec.sort(key=lambda x:x[1])
    newvec_final = [a[0] for a in newvec]
    
    ## Xnew needs to be in order of original integer_vec
    return newvec_final, diffs
    print(integer_vec, file=sys.stderr)
    print(newvec_final, file=sys.stderr)

def vector_prior(data1,cn_vector, scale=False):
    """
    """
    return -np.sum([math.log(data1.prior.get(round(cn_vector[i]), 1e-9)) for i in range(data1.n)])
    #if scale:  return -np.sum([math.log(data1.prior.get(round(cn_vector[i]), 1e-9))*len(data1.samples_include[i]) for i in range(data1.n)])

def subset_shift(frac_vec, integer_vec, data1, frac_ll, integer_ll):
    """
    """
    new_vec,diffs = closest_vector_full(frac_vec, integer_vec)
    new_ll = fopt(new_vec, data1)
    
    integer_prior = vector_prior(data1, integer_vec)
    new_prior = vector_prior(data1, new_vec)
    
    if new_ll + new_prior < integer_ll + integer_prior: 
        print('better-sol', file=sys.stdout, end=' ')
    else: 
        print('worse', file=sys.stdout, end=' ')

    print('LL',frac_ll,integer_ll-frac_ll,new_ll-frac_ll,'|priorLLs',integer_prior,new_prior,file=sys.stdout)

    LLP, newgreedy_vec, initial_ll, newgreedy_ll, stats = greedy_search(new_vec,data1, order='optimized',mindelta=10)
    distance  = closest_vector_full(frac_vec, newgreedy_vec, distance=True)
    newgreedy_prior = LLP-newgreedy_ll
    #counts = Counter(updated[1]);          print(counts,file=sys.stderr)
    print('new-greedy',newgreedy_ll,newgreedy_prior,'distance-frac',distance,file=sys.stdout)
    
    if LLP < integer_ll + integer_prior: 
        print('final solution improved', LLP, file=sys.stdout)
        return LLP, newgreedy_vec, initial_ll, newgreedy_ll
    return None


def greedy_search(init_vec, exome_data, order='optimized', mindelta=0):
    """ 
    Iteratively update CN estimates in the order of fractional copy numbers
    - local greedy search to update until likelihood cannot be improved further
    """
    MAX_ITERATION = 100
    MINCN, MAXCN = exome_data.minCN, exome_data.maxCN
    REFCN = exome_data.refCN
    
    # Order of update
    L = list(range(exome_data.n))
    if order == 'optimized' and exome_data.order != None: 
        L = exome_data.order
   
    # Compute initial likelihood
    prior_ll = vector_prior(exome_data, init_vec) # -np.sum([math.log(exome_data.prior.get(int(init_vec[i]), 1e-9)) for i in range(exome_data.n)])
    current_ll = fopt(init_vec, exome_data) + prior_ll
    initial_ll = current_ll
    iterations, total_updates = 0, 0
    beam = [init_vec[i] for i in range(exome_data.n)]
    
    while iterations < MAX_ITERATION:
        updated = 0
        
        for i in L: # for each sample in order L
            ll_list = []
            current_c = round(beam[i])
            current_ll = fopt(beam, exome_data, sample_list=exome_data.samples_include[i]) - math.log(exome_data.prior.get(int(current_c), 1e-9))
            
            for c in range(max(MINCN,current_c-2), min(current_c+3, MAXCN+1)):  # for each CN of CN range
                beam[i] = max(c, 0.05)
                if c == current_c: 
                    ll = current_ll
                else: 
                    ll = fopt(beam, exome_data, sample_list=exome_data.samples_include[i]) - math.log(exome_data.prior.get(c, 1e-9))
                ll_list.append((-ll,c))

            ll_list.sort(reverse=True)
            llsum =  logsumexp([l[0] for l in ll_list])
            weighted = sum([math.exp(l[0]-llsum)*l[1] for l  in ll_list])
            delta = ll_list[0][0]-ll_list[1][0]
            #print('weigt',llsum,weighted,ll_list,i,file=sys.stderr)
            
            if ll_list[0][1] != current_c and delta > mindelta: 
                #print('updating sample',i,exome_data.samplenames[i],current_c,delta,file=sys.stderr)
                updated += 1
                beam[i] = ll_list[0][1]
            else: beam[i] = current_c
            #beam[i] = weighted # ll_list[0][1]

        iterations += 1
        if updated == 0:      
            break

        total_updates += updated

    #print('fractional',[round(b,2) for b in beam],updated,file=sys.stderr)
    #beam = [round(b) for b in beam]

    ll = fopt(beam, exome_data)
    prior_ll = vector_prior(exome_data, beam) #-np.sum([math.log(exome_data.prior.get(int(beam[i]), 1e-9)) for i in range(exome_data.n)])
    delta = initial_ll - (ll + prior_ll)
    print(f"\ngreedy-Search took {iterations} iterations to complete {updated}.")
    
    return (ll+prior_ll, beam, initial_ll, ll,(iterations,total_updates))


def rounding_rand(x, strength=1.):
    """
    """
    floor = int(np.floor(x))
    ceil = int(np.ceil(x))
    frac = x - floor
    prob_up = frac**strength / (frac**strength + (1 - frac)**strength)
    return floor + (np.random.rand() < prob_up)


def run_greedy_search(exome_data, fractional=None, random_starts=10):
    """
    Runs greedy/beam search for multiple initial CN vectors
    """
    R_cn = range(exome_data.minCN, exome_data.maxCN+1)
    prob_cn = [exome_data.prior.get(i, 1e-9) for i in R_cn]
    s = sum(prob_cn)
    prob_cn = [p/s for p in prob_cn]
    
    def generate_random_start():
        rand_vec = np.random.choice(list(R_cn), size=exome_data.n, p=prob_cn)
        rand_vec.sort()
        init_vec = [1 for  i in range(exome_data.n)]
        k=0
        for index in exome_data.order: 
            init_vec[index] = int(rand_vec[k])
            k +=1 
        #init_vec = [rounding_rand(x) for x in fractional.x]
        return init_vec
        print(Counter(init_vec),'random init', file=sys.stderr)

    prior_dict = {cn:exome_data.prior.get(cn, 1e-9) for cn in R_cn}

    priors = [(exome_data.prior.get(cn, 1e-9), cn) for cn in R_cn]
    priors.sort(reverse=True)
    #print('priors',priors,file=sys.stderr)

    top_candidates= []
    n_random=0
    previous_list = []; previous_lltable = {}
    
    for it in range(random_starts + len(priors)):
        
        print("\n---------------------------------------------------------------------------------", file=sys.stdout)
        if it < random_starts:
            if it==0: 
                init_vec = [round(x) for x in fractional.x] ## simple rounded fractional vector
                init_type = 'frac2int_HMM'
            else:
                init_vec = generate_random_start()
                init_type='random'
            n_random +=1
            print(f"+++ Starting greedy search with Init vector: {init_type}", file=sys.stdout)
        
        else: 
            (prob,copy_num) = priors[it-random_starts]
            if prob > 0.2:  
                init_vec = [copy_num for  i in range(exome_data.n)]
                init_type = 'constant'
                print(f"+++ Starting greedy search with Init vector: {init_type} {copy_num} {prob}", file=sys.stdout)
            else: break

        top_ret = greedy_search(init_vec, exome_data, order='optimized', mindelta=0)
        
        if round(top_ret[0],2) not in previous_lltable: 
            previous_list.append((top_ret[0],top_ret[1],top_ret[0]-top_ret[3]))
            previous_lltable[round(top_ret[0],2)] = 1
            old_solution = False
            (iterations,total_updates) = top_ret[4]

            print(f"greedy-Search took {iterations} iterations, updates= {total_updates}",file=sys.stdout)
            print("Best result:", top_ret[0],top_ret[3],top_ret[2]-top_ret[0],file=sys.stdout)
        
        else: 
            previous_lltable[round(top_ret[0],2)] +=1
            old_solution = True

        if old_solution: 
            continue  ## don't need to run the subset-shift and greedy

        if fractional != None and exome_data.n >=5:
            updated_best = None
            updated_best = subset_shift(fractional.x, top_ret[1], exome_data, fractional.fun, top_ret[3])
            
            if updated_best != None: 
                top_candidates.append(updated_best)
                
                if round(updated_best[0],2) not in previous_lltable: 
                    previous_list.append((updated_best[0], updated_best[1], updated_best[0]-updated_best[3]))
                    previous_lltable[round(updated_best[0], 2)] = 1
                else: 
                    previous_lltable[round(updated_best[0], 2)] +=1
            else: 
                top_candidates.append(top_ret)

        else: 
            top_candidates.append(top_ret)

    NBDlist = []
    print("\n final comparison-----------------best",top_candidates[0][0], file=sys.stdout)
    previous_list.sort()
    
    for i in range(len(previous_list)):
        
        sol = previous_list[i]
        ll = round(sol[0], 2)
        print('cnstring', ll, previous_lltable[ll], sol[2], Counter(sol[1]), file=sys.stdout) 
        
        if i == 0: 
            continue
        
        diff_samples = [k for k in range(exome_data.n) if previous_list[0][1][k] != previous_list[i][1][k]]
        delta1 = previous_list[0][0]-previous_list[i][0]
        delta2 = previous_list[0][2]-previous_list[i][2]
        
        if -delta1/math.log(10) < 20: 
            NBDlist.append((diff_samples, ll, 'loptima'))
        
        if  delta1-delta2 > 0 or delta1 > -30:  
            print('diffs', len(diff_samples), 'LLdelta(prior)', delta1, delta2, 'pair', i, file=sys.stdout)
            #print(diff_samples,[exome_data.samplenames[k] for k in diff_samples],file=sys.stderr)

    top_candidates.sort()
    print("\nGreedy search completed.")
    print("FINAL BEST:", top_candidates[0]) 
    
    delta = sum([ abs(fractional.x[i]-top_candidates[0][1][i]) for i in range(exome_data.n)])
    delta1 = sum([ abs(round(fractional.x[i])-top_candidates[0][1][i]) for i in range(exome_data.n)])
    
    print('difference between fractional and best integer solution', delta, delta1, file=sys.stdout)
    #init_vec_HMM,'\n',top_candidates[0][1],file=sys.stderr)
   
    # Set the CN vector to the best CN vector among the three beams
    exome_data.cn = [c for c in top_candidates[0][1]]
    initial_ll, current_ll = top_candidates[0][2], top_candidates[0][0]
    
    return initial_ll-current_ll, current_ll, NBDlist

def best_fractional(self,random_start=False,max_ab=100000):
    """
    Most likely fractional estimate by optimizing the joint log-likelihood function
    - minimize negative log-likelihood using L-BFGS-B algorithm
    """

    options_bfgs = {'maxfun': 50*self.n, 'maxiter':25*self.n}
    for i in range(self.n):
        s = self.alpha[i] + self.beta[i]
        if s > max_ab: 
            self.alpha[i] *= max_ab/s
            self.beta[i] *= max_ab/s
        #print('alpha,beta',self.alpha[i],self.beta[i],i,self.reference_sets[i],file=sys.stderr)
    
    minCN = max(self.minCN, 1.01)
    maxCN = self.maxCN
    self.cn_bounds = [(minCN,maxCN) for i in range(self.n)]
    
    """
    x0 = [self.refCN] * self.n
    result0 = opt.minimize(fopt,x0,bounds=self.cn_bounds,tol=1e-8,args=(self),method='L-BFGS-B',options=options_bfgs)
    print("Success:", result0.success,result0.message,result0.nfev,result0.njev,result0.nit,file=sys.stderr)
    ll_nopenalty = result0.fun
    print('\nFractional results: Likelihood:', round(result0.fun, 3))
    print('CN-vector:', ' '.join([str(a) for a in sorted(result0.x)]), sep='\n')
    self.compare = sorted([[result0.x[i], i] for i in range(self.n)])
    self.order = [self.compare[i][1] for i in range(self.n)]
    print('frac-results',round(result0.fun,2),[round(float(x),2) for x in result0.x],file=sys.stderr)
    """
    
    #x0 = np.random.choice([i for i in range(1,11)], size=self.n,p=[self.prior.get(cn, 1e-9) for cn in range(1,11)])
    x0 = [self.refCN] * self.n
    priors = [(cn,self.prior.get(cn, 1e-9)) for cn in range(max(1,self.minCN),self.maxCN+1)]
    
    #result1 = opt.minimize(fopt_prior,x0,bounds=self.cn_bounds,tol=1e-7,args=(self,priors),method='L-BFGS-B')
    result1 = opt.minimize(
        fopt_prior,
        x0,
        bounds=self.cn_bounds,
        tol=1e-8,
        args=(self,priors),
        jac=fopt_prior_gradient,
        method='L-BFGS-B',
        options=options_bfgs,
    )
    
    print('frac-results-prior', round(result1.fun,2), [round(float(x),2) for x in result1.x], file=sys.stdout)
    print("Success:", result1.success, result1.message, result1.nfev, result1.njev, result1.nit, file=sys.stdout)
    
    self.compare = sorted([[result1.x[i], i] for i in range(self.n)])
    self.order = [self.compare[i][1] for i in range(self.n)]
    #print('frac-results-prior',round(result1.fun,2),[round(float(x)) for x in result1.x],file=sys.stderr)

    distinct_vals = list(Counter([round(result1.x[i],2) for i in range(self.n)]))
    result1.fun = fopt(result1.x, self) ## calculate likelihood without prior
    
    print('distinct-values', len(distinct_vals), self.n, 'LL-noprior', result1.fun, file=sys.stdout)
    #for s in range(self.n): print('sample',s,self.reference_sets[s],self.samplenames[s],self.alpha[s]+self.beta[s],file=sys.stderr)
    
    return result1


def graph_cut(data1):
    """
    For adding graph-cut based subsets
    """
    G = nx.Graph()
    for i in range(data1.n):
        G.add_edges_from([(i,s) for s in data1.reference_sets[i]])
    #recursive_min_cut(G)
    
    cut_value, partition = nx.stoer_wagner(G)
    cut_edges = [(u, v) for u in partition[0] for v in G[u] if v in partition[1]]
    
    if cut_value < 3: 
        print('edges',cut_edges,file=sys.stderr)
    print('cut',cut_value,partition,file=sys.stderr)
    spartition = partition[0]
    
    if len(partition[0]) > len(partition[1]): 
        spartition = partition[1]
    
    if len(spartition) > 1: 
        
        for s in spartition: 
            best[s] = best[s]-1
        ll = round(float(fopt(best, data1) +vector_prior(data1,best) ),2)
        print('cutLL',spartition,ll,file=sys.stderr)
        
        for s in spartition: 
            best[s] = best[s]+2
        ll = round(float(fopt(best, data1) +vector_prior(data1,best) ),2)
        print('cutLL',spartition,ll,file=sys.stderr)
        
        for s in spartition: 
            best[s] = best[s]-1

def compute_qual_final(data1, best, NBDlist=None, UPDATE=True):
    """
    Consider neighborhoods (single-sample) + local optima from greedy heuristic + reference-set-based
    """
    sf = -10/np.log(10)
    global_ll = fopt(best, data1) + vector_prior(data1, best)
    LLtable = [[] for i in range(data1.n)]

    NBDlist_add = []
    ## use referencesets to find more neighborhood solutions
    updates=0
    newvec = [best[j] for j in range(data1.n)]
    
    for i in range(data1.n):
        
        slist = [i] + data1.reference_sets[i]
        
        for s in slist: 
            newvec[s] +=1 
        new_ll1 = fopt(newvec, data1) + vector_prior(data1, newvec)
        
        for s in slist: 
            newvec[s] = max(newvec[s]-2,0.05)
        new_ll2 = fopt(newvec, data1) + vector_prior(data1, newvec)
        
        for s in slist: 
            newvec[s] = best[s]
        delta = min(new_ll1, new_ll2)-global_ll
        
        if delta < -1 and UPDATE: 
            #print('ERROR new likelihood is better than best',slist,i,delta,global_ll,file=sys.stdout)
            #print('names',[(data1.samplenames[s],best[s]) for s in slist],file=sys.stderr)
            if new_ll1 < new_ll2: 
                for s in slist: 
                    best[s] +=1 
                i = 0
                global_ll = new_ll1
                updates +=1
                continue 
            else: 
                for s in slist: 
                    best[s] -=1 
                i =  0
                global_ll = new_ll2
                updates +=1
                continue 

        if delta > 0 and delta/np.log(10) < 20 and not UPDATE: ## qualvalue = 200
            #print('adding sample+refset as NBD',slist,i,round(delta,2),file=sys.stderr)
            NBDlist_add.append((slist, delta+global_ll, 'refset'))
        if updates >= 5: break
  
    if not UPDATE: 
        NBDlist += NBDlist_add 
        NBDlist.sort(key=lambda x: x[1])
        prev = (0,0)
        
        for nbd in NBDlist: 
            
            ## we don't want to add the same subset twice... or use singleton sets
            if len(nbd[0]) ==1 or set(nbd[0]) == set(prev): 
                continue 
            
            print('nbdsol log10-diff',round((nbd[1]-global_ll)/math.log(10),1),nbd[0],nbd[2],file=sys.stdout)
            for s in nbd[0]: LLtable[s].append(float(global_ll-nbd[1]))
            prev = nbd[0]

    for it in range(10):

        updates=0
        for i in range(data1.n):
            cn = best[i]
            ll_list = []
            ll_list1 = []
            
            for j in range(max(data1.minCN,cn-2), min(data1.maxCN+1,cn+3)):
                best[i] = max(j, 0.05)
                ll = fopt(best, data1, sample_list=data1.samples_include[i]) - math.log(data1.prior.get(int(cn), 1e-9))
                ll_list.append(ll)
                ll_list1.append((round(float(ll),1),j))
                if j == cn: 
                    refll = ll 
            ll_list1.sort()
            
            if int(cn) != ll_list1[0][1] and UPDATE: ## local update
                cn = ll_list1[0][1]
                refll = ll_list1[0][0]
                updates += 1
                #print('not local optima',i,data1.samplenames[i],best[i],cn,ll_list1,file=sys.stdout)
            #print('LL for sample',i,best[i],ll_list1,file=sys.stderr)

            for ll in ll_list: 
                LLtable[i].append(float(refll-ll))
            best[i] = int(cn)
            data1.cn[i] = cn
        
        if updates == 0: 
            break

    if UPDATE: 
        return best
    
    ## quality is always calculated using sorted list of likleihoods...
    lowquality = 0
    
    for i in range(data1.n):
        LLtable[i].sort(reverse=True)
        quality = round(sf* (logsumexp(LLtable[i][1:])-logsumexp(LLtable[i])) )
        data1.cn_quality[i] = quality
        if quality < 20: 
            print('sample ',i,'CN',best[i],'QUAL',data1.cn_quality[i],data1.samplenames[i],data1.alpha[i]+data1.beta[i],LLtable[i][0:2],file=sys.stdout)
            lowquality += 1
    
    print('lowquality',lowquality,data1.n,global_ll,file=sys.stdout)


def analyze_component(inp, data1, refCN, c, comp, outname, logfile=sys.stdout):
    """
    Perform all analyses on a given connected component
    """
    print()
    print('='*50)
    print('Analyzing connected component', c, comp, data1.n)
    print('='*50)
    timer = perf_counter()

    data1.refCN = refCN
    data1.cn = [data1.refCN] * data1.n

    # Display reference set for each sample
    for j,refset in enumerate(data1.reference_sets):
        print(j,data1.samplenames[j], end='\t')
        for i in refset:
            print(i, data1.samplenames[i], end='; ')
        print()
    data1.samples_include = [] # for efficient ()
    
    for i in range(data1.n): 
        data1.samples_include.append([i])
    
    for i in range(data1.n):
        for j in data1.reference_sets[i]: 
            data1.samples_include[j].append(i)

    #scale_factor = sum([math.log(len(data1.reference_sets[i])+2,2) for i in range(data1.n)])/data1.n
    #data1.scale_factor = 2.0/scale_factor
    data1.scale_factor = 1
    #print('scalefac',data1.scale_factor,file=sys.stderr)
    
    n_exons = len(data1.ecounts[0])
    data1.comb_log = [0.0 for i in range(data1.n)]
    for i in range(data1.n):
        for e in range(n_exons):
            c = data1.ecounts[i][e]
            Rc = sum([data1.ecounts[j][e] for j in data1.reference_sets[i]])
            if c > 0 and Rc > 0:        
                data1.comb_log[i] += sum([ math.log((Rc+c-k+1)/k) for k in range(1,c+1)])
            #data1.comb_log += gammaln(c+Rc + 1) - gammaln(c + 1) - gammaln(Rc + 1)
    
    print(sum([c for c in data1.comb_log]), 'LLconst for C(n,r)', file=sys.stdout)


    timer0 = perf_counter()
    print(f"step-0: [{str(timedelta(seconds = timer0 - timer))[:-3]}]")

    # ========================================================================
    # (1) Compute fractional estimates
    
    data1.cn = [data1.refCN] * data1.n
    result = best_fractional(data1)
    #print("fractional-order:", data1.order)
    #print("fractional-cn:", [x for x in result.x])

    #timer1 = perf_counter()
    #print(f"step-1: [{str(timedelta(seconds = timer1 - timer0))[:-3]}]",file=sys.stderr) 
    # ========================================================================
    # (2) Compute integer estimates (local updates using greedy/beam search)
    data1.cn = [refCN] * data1.n
    delta, current,NBDlist = run_greedy_search(data1,fractional=result)  
    best, best_score = [c for c in data1.cn], current
    print(best, best_score)

    greedy_order = [a[1] for a in sorted([(data1.cn[i], i) for i in data1.order], key=lambda x: x[0])]
    #print("greedy-order:", greedy_order)
    print("greedy-cn:", best)

    #timer2 = perf_counter()
    #print(f"step-2: [{str(timedelta(seconds = timer2 - timer1))[:-3]}]",file=sys.stderr) 

    #print('debug',best,best_score,file=sys.stderr)
    #best1 = compute_qual_final(data1,best,NBDlist=None,UPDATE=True)
    compute_qual_final(data1,best,NBDlist,UPDATE=False)
    
    fr = np.array([x for x in result.x]) 
    gr = np.array([b for b in data1.cn]) 
    fr_ll, gr_ll = fopt(fr, data1), fopt(gr, data1)


    #print('compare',ll1,ll1-prior1,'|',ll2,ll2-prior2,'|',ll3,ll3-prior3,file=sys.stderr)
    #for i in range(data1.n): print('finalQ',i,data1.cn[i],data1.cn_quality[i],data1.alpha[i],data1.beta[i])
    #print("Quality:", data1.cn_quality)

    # ========================================================================
    # (6) Reorganize data for export
    
    QUAL_THRESH, MIN_CC_SIZE = inp.qual_thresh, inp.min_cc_size
    data1.cn = [0 if cn==0.05 else cn for cn in data1.cn]
    data1.fracCN = fr
    
    #df_data = pd.DataFrame([[data1.genename]*len(fr), data1.samplenames, [c]*len(fr), 
    #                        [round(f,6) for f in fr], data1.cn, data1.cn_quality, 
    #                        [len(r) for r in data1.reference_sets]]).T
    
    #df_data.columns = ['gene', 'name', 'comp', 'frac', 'grdy', 'qual', 'refset_size']
    #df_data.frac = df_data.frac.map('{:.6f}'.format)
    #df_data['trust'] = np.where((df_data['qual'] >= QUAL_THRESH) & 
    #                            (comp > MIN_CC_SIZE), True, False)
    
    #if c == len(data1.comp_list)-1:
    #    for s in data1.samples_noref:
    #        noref_row = [data1.genename, s, '.', '.', '.', 0, 0, False]
    #        df_data = pd.concat([df_data, pd.DataFrame([noref_row], 
    #                            columns=df_data.columns)], 
    #                            ignore_index=True)

    #df_data2 = pd.DataFrame([[data1.genename], [c], [fr_ll], [gr_ll]]).T
    #df_data2.columns = ['gene', 'comp', 'frac_lkhd', 'grdy_lkhd']
    #df_data2['trust'] = np.where(comp > 5, True, False)
    
    #outdir, outfname = os.path.dirname(outname), os.path.basename(outname)
    #outfp = f"{outdir}/{outfname}"
    #outfp2 = f"{outdir}/{outfname}.lkhd"
    #
    #if c == 0:
    #    df_data.to_csv(outfp, sep='\t', index=False)
    #    df_data2.to_csv(outfp2, sep='\t', index=False)
    #else:
    #    df_data.to_csv(outfp, sep='\t', mode='a', index=False, header=False)
    #    df_data2.to_csv(outfp2, sep='\t', mode='a', index=False, header=False)
    
    # ========================================================================
    # (7) Summarize results
    
    print('\nSTATS', delta, current, fopt(data1.cn, data1), data1.n)
    print('- bestsol:', ''.join([str(a) for a in data1.cn]))
    #print('- subs:', ''.join(['^' if i in diff_idx else '.' for i in range(data1.n)])) 
    

# ============================================================================
# MAIN FUNCTION
# ============================================================================

def run(inp):

    # Set up exome data class
    data = ExomeData()
    data.get_parameters2(inp.all_stat_fp, inp.priors_fp, inp.loci_name)
    data.gene_counts(inp.gene_cnts_fp)
    refCN = inp.refcn

    # Make final output directory
    os.makedirs(inp.finaldir, exist_ok=True)
    final_out_fp = os.path.join(inp.finaldir, f'{inp.loci_name}.out')
    final_log_fp = os.path.join(inp.finaldir, f'{inp.loci_name}.log')
    
    # Make directory to store connected-component objs
    os.makedirs(inp.ccobjdir, exist_ok=True)
    
    timer_start = perf_counter()
    
    # print all stdout outputs to logfile
    original_stdout = sys.stdout
    with open(final_log_fp, 'w') as logf:
        
        # For each connected component
        c = 0
        cc_list = []
        
        print(f"{data.genename} has {len(data.comp_list)} connected components.")
        for comp in data.comp_list:
            
            # Check if already estimated (stored as file)
            cc_fp = os.path.join(inp.ccobjdir, f'{inp.loci_name}.{c}.est') 
            if os.path.exists(cc_fp):
                print(f"- Found existing object for {inp.loci_name}'s component-{c}.")
                with open(cc_fp, 'rb') as cc_in:
                    cc_obj = pickle.load(cc_in)
                cc_list.append(cc_obj)
                c += 1
                continue
            
            data_cc = data.construct_subset(c)
            print(f"- Estimating agCNs for {data_cc.genename} component {c} [{comp} samples].")
            
            # Estimate agCN
            sys.stdout = logf
            analyze_component(inp, data_cc, refCN, c, comp, final_out_fp)

            # Save connected component object
            with open(cc_fp, 'wb') as cc_out:
                pickle.dump(data_cc, cc_out, protocol=pickle.HIGHEST_PROTOCOL)
            
            c += 1
            cc_list.append(data_cc)
            sys.stdout = original_stdout

        # Estiamte HMM path
        if not inp.skiphmm:
            try:
                _hmm.run_hmm(inp, cc_list)
            except IndexError as e:
                print(f"ERROR [INDEX_ERROR] {inp.loci_name}")
                print(e)

        timer_end = perf_counter()
        print(f"ELAPSED: [{str(timedelta(seconds = timer_end - timer_start))[:-3]}] ({data.genename})")

