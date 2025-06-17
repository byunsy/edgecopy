import os
import re
import math
import pandas as pd
import numpy as np
from scipy.special import logsumexp

from .exomecounts import *
from .optimize_functions import *

class cnHMM:
    
    # ========================================================================
    # INITIALIZATION
    # ========================================================================
    def __init__(self, hidden, obsrvs, num_samples, reference_sets, refCN, 
                 means, alpha, beta, genename):

        self.genename = genename
        self.hidden = np.array(hidden) # hidden states (copy number states)
        self.obsrvs = np.array(obsrvs) # observations (read counts for each exon)
        self.num_hidden = self.hidden.shape[0]
        self.num_obsrvs = self.obsrvs.shape[1]
        self.num_samples = num_samples
        self.init_p = None
        self.tran_p = None
        self.emis_p = None
        self.cn_mat = None
        self.refset = reference_sets
        self.refsum = np.array([self.obsrvs[refset].sum(axis=0) for refset in reference_sets])
        self.refCN  = refCN
        self.refmade = {}
        self.means = np.array(means)
        self.data_order = None
        self.debug_fp = None
        self.point_cn = None
        self.ground_truth = None
        self.scores = [[]]*self.num_samples
        self.vpaths = [[]]*self.num_samples
        self.backtr = [[]]*self.num_samples
        self.bestLL = [[]]*self.num_samples
        self.forward_p = [[]]*self.num_samples
        self.backward_p = [[]]*self.num_samples
        self.posterior = [[]]*self.num_samples
        self.cn_mat_fb = [[]]*self.num_samples
        self.alpha = alpha
        self.beta = beta
        self.scale_factor = 1.0
        self.comb_log = []

    # ========================================================================
    # SET AND SAVE FUNCTIONS
    # ========================================================================
    def set_initial_p(self, init_prob):
        self.init_p = init_prob    
        
    def set_transition_p(self, tran_matrix):
        self.tran_p = tran_matrix
        
    def set_emission_p(self, emis_matrix):
        self.emis_p = emis_matrix

    def init_emission_p(self):
        self.emis_p = np.full((self.num_samples, self.num_hidden, self.num_obsrvs), -np.inf)
        
    def set_cn_mat(self, cn_mat):
        self.cn_mat = cn_mat
        
    def set_data_order(self, cc):
        sorted_order = [s[0] for s in sorted(enumerate(cc.cn), key=lambda x: x[1])]
        self.data_order = sorted_order
        self.point_cn = cc.cn

    def set_debug_fp(self, inp):
        self.debug_fp = inp.debug_fp

    def set_references_made(self):
        for i in range(len(self.refset)):
            self.refmade[i] = []
            for j,j_refset in enumerate(self.refset):
                if i in j_refset:
                    self.refmade[i].append(j) 
        
    def save_scores(self, scores, s_idx):
        if self.scores[s_idx]:
            self.scores[s_idx].append([s for s in scores])
        else:
            self.scores[s_idx] = [[s for s in scores]]
    
    def save_hidden_path(self, hidden_path, s_idx):
        if self.vpaths[s_idx]:
            self.vpaths[s_idx].append(hidden_path)
        else:
            self.vpaths[s_idx] = [hidden_path]
            
    def save_backtrack(self, backtrack, s_idx):
        if self.backtr[s_idx]:
            self.backtr[s_idx].append(backtrack)
        else:
            self.backtr[s_idx] = [backtrack]
    
    def save_bestLL(self, LL, s_idx):
        if self.bestLL[s_idx]:
            self.bestLL[s_idx].append(LL)
        else:
            self.bestLL[s_idx] = [LL]
            
    # ========================================================================
    # COMPUTATION FUNCTIONS
    # ========================================================================        
    def compute_initial_p(self, prior_fp, cc):
        """
        prior_fp: path to CN distribution of a particular gene
        """
        prior_d = np.array(pd.read_csv(prior_fp, sep='\t'))[:,1]
        i_probs = np.full((self.num_samples, self.num_hidden), prior_d)
        #for s_idx in range(self.num_samples):
        #    
        #    # Also take into account our point estimate of a sample
        #    # If point_estim != refcn, give more weight to point_estim
        #    # and less to refcn
        #    if cc.cn[s_idx] != self.refCN:
        #        DELTA = i_probs[s_idx, self.refCN-1] - 0.1
        #        i_probs[s_idx, cc.cn[s_idx]-1] += DELTA
        #        i_probs[s_idx, self.refCN-1] -= DELTA
        
        self.set_initial_p(np.log(i_probs))
    
    def compute_transition_p(self, t):
        """
        Compute transition probabilities for each sample
        -  t: probability of transitioning from one state to another state
        """
        all_trans_p = []
        for s_idx in range(self.num_samples):
            all_trans_p.append(self.compute_sample_transition_p(s_idx, t))
        self.set_transition_p(np.array(all_trans_p))
    
    def compute_sample_transition_p(self, s_idx, t):
        """
        Compute matrix of probability of moving from one state (ROW) to another (COL). 
        Returns a 10x10 matrix of negative log-likelihoods.
        -  t: probability of transitioning from one state to another state
        -  STARTCN: index of point CN estimate
        -  CNRANGE: indices of CN range given MAXJUMP
        """
        MAXJUMP = 2
        STARTCN = int(self.point_cn[s_idx])
        CNRANGE = range(max(STARTCN-MAXJUMP,0), min(STARTCN+MAXJUMP+1,10))
        t_probs = np.full((self.num_obsrvs-1, self.num_hidden, self.num_hidden), np.log(1e-20))

        for from_cn in CNRANGE:
            to_cn_range = range(max(from_cn-MAXJUMP,0), min(from_cn+MAXJUMP+1,10))
            for to_cn in to_cn_range:
                if from_cn == to_cn:  # Remain in current CN state
                    t_probs[:,from_cn,to_cn] = np.log(1-t)
                else:  # Jump to nearby CN states
                    t_probs[:,from_cn,to_cn] = np.log(t/(len(to_cn_range)-1))
        return t_probs
    def update_transition_p(self, t):
        
        quals = []
        for i,row in enumerate(self.cn_mat.astype(np.int32)):
            quals.append(self.compute_path_qual(i,row))
        hmm_qual = np.array(quals)

        all_trans_p = []
        for s_idx in range(self.num_samples):
            #all_trans_p.append(self.udpate_sample_transition_p(s_idx, t, hmm_qual))
            all_trans_p.append(self.udpate_sample_transition_p2(s_idx, t, hmm_qual))
        self.set_transition_p(np.array(all_trans_p))
    
    def udpate_sample_transition_p(self, s_idx, t, hmm_qual, qual_thresh=20):

        MAXJUMP = 2
        STARTCN = int(self.point_cn[s_idx])
        CNRANGE = range(max(STARTCN-MAXJUMP,0), min(STARTCN+MAXJUMP+1,10))
        t_probs = self.tran_p[s_idx]
        #t /= 4

        # For each exon a -> exon a+1 transition, find all unique CN transitions
        for i in range(self.cn_mat.shape[1]-1):

            hmm_cn_mat = self.cn_mat[:, i:i+2].copy()
            
            # Both from and to CN states must have qual > 20
            conditions = (hmm_qual[:, i] > qual_thresh) & (hmm_qual[:, i+1] > qual_thresh)
            hmm_cn_mat = hmm_cn_mat[conditions]

            unique_trans, counts = np.unique(hmm_cn_mat, axis=0, return_counts=True)
            #unique_trans = np.array([unique_trans[x] for x,c in enumerate(counts) if c>=2])
            uniq_from_ = np.unique(unique_trans[:,0])

            # possible 'from_cn's based on CN profiles
            for from_cn in uniq_from_:

                # possible 'to_cn's and their frequencies based on CN profiles
                uniq_to_cn = {}
                for j,trans in enumerate(unique_trans):
                    if trans[0]==from_cn:
                        uniq_to_cn[int(trans[1])] = int(counts[j])

                from_cn = int(from_cn)
                to_cn_range = range(max(from_cn-MAXJUMP,0), min(from_cn+MAXJUMP+1,10))
                for to_cn in to_cn_range:

                    # Number of valid to_cn states not observed in CN profiles
                    NUM_UNOBS = (len(to_cn_range)-len(uniq_to_cn))

                    # Observed in CN profiles                
                    if to_cn in uniq_to_cn:
                        to_cn = int(to_cn)
                        p = uniq_to_cn[to_cn] / sum(uniq_to_cn.values())
                        t_probs[i,from_cn,to_cn] = np.log(p - (t*NUM_UNOBS/len(uniq_to_cn)))

                    # Not observed in CN profiles
                    else:
                        t_probs[i,from_cn,to_cn] = np.log(t)
        return t_probs  
        
    def udpate_sample_transition_p2(self, s_idx, t, hmm_qual, qual_thresh=20):
        """
        """
        MAXJUMP = 2
        STARTCN = int(self.point_cn[s_idx])
        CNRANGE = range(max(STARTCN-MAXJUMP,0), min(STARTCN+MAXJUMP+1,10))
        t_probs = self.tran_p[s_idx]
        MIN_HQ_CNT = 1

        # For each exon a -> exon a+1 transition, find all unique CN transitions
        for i in range(self.cn_mat.shape[1]-1):
            hmm_cn_mat = self.cn_mat[:, i:i+2].copy()

            # Both from and to CN states must have qual > 20
            conditions = (hmm_qual[:, i] > qual_thresh) & (hmm_qual[:, i+1] > qual_thresh)
            hmm_cn_mat = hmm_cn_mat[conditions]

            # Count all instances CN changes
            diffs = hmm_cn_mat[:,1] - hmm_cn_mat[:,0]
            diffs = diffs[diffs!=0]
            unique_diffs, cnts = np.unique(diffs, return_counts=True)
            uniq_d_dict = {int(d):c for d,c in zip(unique_diffs, cnts) if c >= MIN_HQ_CNT}

            for from_cn in CNRANGE:

                to_cn_range = range(max(from_cn-MAXJUMP,0), min(from_cn+MAXJUMP+1,10))
                for to_cn in to_cn_range:

                    if from_cn == to_cn:
                        continue
                    
                    # Update existing tprob based on observed CN dups/dels
                    if to_cn - from_cn in uniq_d_dict:
                        
                        p = t + (uniq_d_dict[to_cn - from_cn] / sum(uniq_d_dict.values()))
                        t_probs[i,from_cn,to_cn] = np.log(p)
                        
                # Normalize the row (in log-space)
                t_probs[i,from_cn,:] -= logsumexp(t_probs[i,from_cn,:])
                
        return t_probs

    def compute_emission_p(self, means, cn_mat, soi=None, MIN_CN=0.05):
        """
        Compute emission probability using beta-binomial model
        means  : mean read depth of each sample
                 - list of FLOAT values (mean read depths)
        cn_mat : exon-level CN estimate of each sample
                 - matrix (sample x exons) of INT values (CN estimates)
        """ 
        def compute_ll(self, s, obs):
                
            # Obtain alpha and beta
            alpha, beta = self.alpha[s], self.beta[s]

            # Given alpha and beta, compute new alpha and beta
            a1 = alpha*cn_mat[s,obs]
            b1_denom = sum(means[j] for j in self.refset[s])
            b1 = beta*sum([cn_mat[j,obs]*means[j]/b1_denom for j in self.refset[s]])
            alpha_new = a1*(alpha+beta)/(a1+b1)
            beta_new = b1*(alpha+beta)/(a1+b1)

            # Ensure a1 and b1 is non-neg (e.g. when cn_mat[i]==0)
            #a1, b1 = max(a1, 1e-6), max(b1, 1e-6)
            
            # Compute beta log-likelihood 
            c = self.obsrvs[s,obs]
            Rc = self.refsum[s,obs]
            if not alpha_new or not beta_new:
                print("!!!!", c, alpha_new, Rc, beta_new)
            beta_ll = betaln(c+alpha_new, Rc+beta_new) - betaln(alpha_new, beta_new)
            return beta_ll

        emis_matrix = self.emis_p
        samples_to_update = range(self.num_samples)
        if soi:
            samples_to_update = [soi]
        
        for s_idx in samples_to_update: 
            for obs_idx in range(self.num_obsrvs):
                c_tmp = cn_mat[s_idx, obs_idx]
                for cn_idx,cn in enumerate(self.hidden):
                    cn_mat[s_idx,obs_idx] = max(cn, MIN_CN)
                    sample_ll = compute_ll(self, s_idx, obs_idx)
                    refmade_ll = sum([compute_ll(self, r, obs_idx) for r in self.refmade[s_idx]])
                    emis_matrix[s_idx, cn_idx, obs_idx] = sample_ll + refmade_ll
                cn_mat[s_idx, obs_idx] = max(c_tmp, MIN_CN)
                
        self.set_emission_p(emis_matrix)
        
    def compute_emission_p_binom(self, means, cn_mat, soi=None):
        """
        Compute emission probabilities using binomial model
        means  : mean read depth of each sample
                 - list of FLOAT values (mean read depths)
        cn_mat : exon-level CN estimate of each sample
                 - matrix (sample x exons) of INT values (CN estimates)
        """ 
        def compute_ll(self, s, obs):
            numer = means[s]*cn_mat[s,obs]
            denom = numer + sum([means[j]*cn_mat[j,obs] for j in self.refset[s]])
            p = numer/denom
            ll = np.log(p)*self.obsrvs[s,obs] + np.log(1.0-p)*self.refsum[s,obs]
            return ll

        emis_matrix = self.emis_p
        samples_to_update = range(self.num_samples)
        if soi:
            samples_to_update = [soi]
        
        K = 10
        for s_idx in samples_to_update: 
            for obs_idx in range(self.num_obsrvs):
                c_tmp = cn_mat[s_idx,obs_idx]
                for cn_idx,cn in enumerate(self.hidden):
                    cn_mat[s_idx,obs_idx] = cn
                    sample_ll = compute_ll(self, s_idx, obs_idx)
                    refmade_ll = sum([compute_ll(self, r, obs_idx) for r in self.refmade[s_idx]])
                    emis_matrix[s_idx, cn_idx, obs_idx] = (sample_ll + refmade_ll)/scale_f(K)
                cn_mat[s_idx,obs_idx] = c_tmp
                
        self.set_emission_p(emis_matrix)
        
    def compute_first_LL(self, s_idx, hs_idx):
        """
        Compute the initial log-likelihood of path before running viterbi
        """
        score = self.init_p[s_idx, hs_idx] + self.emis_p[s_idx, hs_idx, 0]
        for obs_idx in range(1, self.num_obsrvs):
            trans_p = self.tran_p[s_idx, obs_idx-1, hs_idx, hs_idx]
            emiss_p = self.emis_p[s_idx, hs_idx, obs_idx]
            score += (trans_p + emiss_p)
        return score
    
    def reconstruct_path(self, score_of, backtrack):
        """
        Starting from the last column, find the max state
        and work backwards to reconstruct the path with the highest likelihood
        """
        # Find the state with total max sum weight
        last_col = [score_of[hs_idx][self.num_obsrvs-1] for hs_idx in range(self.num_hidden)]
        max_idx = last_col.index(max(last_col))
        path = self.hidden[max_idx]
        
        # Reconstruct path starting from max_state        
        for j in range(self.num_obsrvs)[:0:-1]:
            path = str(self.hidden[backtrack[j,max_idx]]) + ',' + str(path)
            max_idx = backtrack[j,max_idx]       

        return path 
    
    def viterbi(self, s_idx):
        """
        Runs Viterbi algorithm on HMM for a single sample
        """
        backtrack = np.full((self.num_obsrvs, self.num_hidden), -1)
        score_of = np.full((self.num_hidden, self.num_obsrvs), -np.inf)

        # Set initial states/probs
        for hs_idx in range(self.num_hidden):
            score_of[hs_idx, 0] = self.init_p[s_idx, hs_idx] + self.emis_p[s_idx, hs_idx, 0]
            
        for obs_idx in range(1, self.num_obsrvs):
            for curr_hs in range(self.num_hidden):
                emiss_p = self.emis_p[s_idx, curr_hs, obs_idx]
                for prev_hs in range(self.num_hidden):

                    # Compute score
                    trans_p = self.tran_p[s_idx, obs_idx-1, prev_hs, curr_hs]
                    score = score_of[prev_hs, obs_idx-1] + trans_p + emiss_p

                    # Get the max sum weight (since positive log-likelihood)
                    if score > score_of[curr_hs][obs_idx]:
                        score_of[curr_hs][obs_idx] = score
                        backtrack[obs_idx][curr_hs] = prev_hs
                        
        return score_of, backtrack
    
    def print_debug(self, cc, cn_matrix, cn_matrix_tmp):
        """
        Print out key information for debugging purposes
        """
        with open(self.debug_fp, 'a+') as debug_f:
            for i in range(cn_matrix.shape[0]):
                if not (cn_matrix[i] == cn_matrix_tmp[i]).all():
                    sample = cc.samplenames[i]
                    gene = cc.genename
                    refmade_len = len(self.refmade[i])
                    
                    print("="*60, file=debug_f)
                    print(sample, file=debug_f)
                    print("="*60, file=debug_f)
                    print('', file=debug_f)
                    
                    print(f'Sample index: {i}', file=debug_f)
                    print(f'Point agCN estimate : {cc.cn[i]}', file=debug_f)
                    print(f'Vector agCN estimate: {cn_matrix[i]} --> {cn_matrix_tmp[i]}', file=debug_f)
                    print('', file=debug_f)
                    
                    print(f"Number of samples that {sample} is in the reference of:", file=debug_f)
                    print(refmade_len, file=debug_f)
                    print('', file=debug_f)
                    
                    print(f"Emission matrix (CNs x exons) of {sample}:", file=debug_f)
                    print(np.round(self.emis_p[i]), file=debug_f)
                    print('', file=debug_f)

                    print(f"Read counts for {sample}", file=debug_f)
                    print(self.obsrvs[i], file=debug_f)
                    print('', file=debug_f)
                    
                    print(f"Read counts for {sample}'s reference set", file=debug_f)
                    for r in self.refset[i]:
                        print(self.obsrvs[r], file=debug_f)
                    print('', file=debug_f)

        print(f"Debug log file saved to: {self.debug_fp}")

    def compute_exon_CN(self, cc, MAX_ITER):
        """
        Check global LL after running viterbi on all samples 
        """
        # Initialize CN matrix based on point CN estimates
        cn_matrix = np.array([[s_cn]*self.num_obsrvs for s_cn in cc.cn], dtype=np.float64)
        
        # For each sample, find for which samples it is the reference of
        self.set_references_made()
        
        # Initialize emission probabilities
        self.init_emission_p()
        self.compute_emission_p(cc.means, cn_matrix)
        
        # Initialize hidden path and best log-likelihood
        # - updated paths should ideally have improved LL
        for s_idx,s_cn in enumerate(cc.cn):
            init_hidden_path = ','.join([f'{s_cn}']*self.num_obsrvs)
            self.save_hidden_path(init_hidden_path, s_idx)
            init_bestLL = self.compute_first_LL(s_idx, s_cn-1)
            self.save_bestLL(init_bestLL, s_idx)

        fopt_sofar = fopt_exon(self, cn_matrix)
        print("Initial:", fopt_sofar)

        ITER = 0
        cn_matrix_tmp = cn_matrix.copy()

        while ITER < MAX_ITER:
            for s_idx in self.data_order:
                
                # Run Viterbi algorithm for sample s_idx
                score_of, backtrack = self.viterbi(s_idx)
                self.save_scores(score_of, s_idx)
                self.save_backtrack(backtrack, s_idx)
                
                # Iff score improved from latest bestLL, then update path and score
                score = max([s[-1] for s in score_of])
                
                if score > self.bestLL[s_idx][-1]:
                    hidden_path = self.reconstruct_path(score_of, backtrack)
                    self.save_hidden_path(hidden_path, s_idx)
                    self.save_bestLL(score, s_idx)
                else:
                    hidden_path = self.vpaths[s_idx][-1]
                    
                #try:
                #    cn_matrix_tmp[s_idx] = np.array(list(map(int, hidden_path.split(','))))
                #except AttributeError:
                #    cn_matrix_tmp[s_idx] = np.array([hidden_path])

                cn_matrix_tmp[s_idx] = np.array(list(map(int, str(hidden_path).split(','))))
                #self.compute_emission_p(cc.means, cn_matrix, s_idx)
            
            # Check which vectors have changed
            for i in range(cn_matrix.shape[0]):
                if not (cn_matrix[i] == cn_matrix_tmp[i]).all():
                    print(f'{i}\t{cc.cn[i]}\t{cn_matrix[i]} --> {cn_matrix_tmp[i]}')
                    
            if self.debug_fp:
                self.print_debug(cc, cn_matrix, cn_matrix_tmp)
            
            # If global LL of new matrix is better, update
            fopt_new = fopt_exon(self, cn_matrix_tmp)
            if fopt_new < fopt_sofar:
                cn_matrix = cn_matrix_tmp.copy()
                self.compute_emission_p(cc.means, cn_matrix, s_idx)
                print(f"Iteration {ITER+1}:", fopt_new)
                ITER += 1
                fopt_sofar = fopt_new
            else:
                print(f'old: {fopt_sofar}\tnew: {fopt_new}')
                print(f"LL did not improve after iteration {ITER}. Stop iterating.")
                break
        print()
        cn_matrix = self.check_empty_rc(cn_matrix)
        self.set_cn_mat(cn_matrix)
    
    def compute_exon_CN2(self, cc, MAX_ITER, MIN_CN=0.05):
        """
        Check global LL after running viterbi on all samples 
        """
        # Initialize CN matrix based on point CN estimates
        cn_matrix = np.array([[s_cn]*self.num_obsrvs for s_cn in cc.cn], dtype=np.float64)
        cn_matrix[cn_matrix == 0] = MIN_CN  # change CN=0 to MIN_CN for better computation

        # For each sample, find for which samples it is the reference of
        self.set_references_made()
        
        # Initialize emission probabilities
        self.init_emission_p()
        self.compute_emission_p(cc.means, cn_matrix)
        
        # Initialize hidden path and best log-likelihood
        # - updated paths should ideally have improved LL
        for s_idx,s_cn in enumerate(cc.cn):
            init_hidden_path = ','.join([f'{s_cn}']*self.num_obsrvs)
            self.save_hidden_path(init_hidden_path, s_idx)
            init_bestLL = self.compute_first_LL(s_idx, s_cn) #####
            self.save_bestLL(init_bestLL, s_idx)

        fopt_sofar = fopt_exon(self, cn_matrix)
        print(f"[{self.genename}] Initial LL:", fopt_sofar)

        ITER = 0
        while ITER < MAX_ITER:
            for s_idx in self.data_order:
                
                # Run Viterbi algorithm for sample s_idx
                score_of, backtrack = self.viterbi(s_idx)
                self.save_scores(score_of, s_idx)
                self.save_backtrack(backtrack, s_idx)
                
                # Iff score improved from latest bestLL, then update path and score
                score = max([s[-1] for s in score_of])
                
                if score > self.bestLL[s_idx][-1]:
                    hidden_path = self.reconstruct_path(score_of, backtrack)
                    self.save_hidden_path(hidden_path, s_idx)  # update path
                    self.save_bestLL(score, s_idx)             # update score
                else:
                    hidden_path = self.vpaths[s_idx][-1]
                
                # Update CN path of sample s_idx right away
                cn_matrix[s_idx] = np.array(list(map(int, str(hidden_path).split(','))))
                cn_matrix[cn_matrix == 0] = MIN_CN
                self.compute_emission_p(cc.means, cn_matrix, s_idx)
             
            fopt_new = fopt_exon(self, cn_matrix)
            if fopt_new < fopt_sofar: 
                ITER += 1
                print(f"[{self.genename}] Iteration {ITER}: {fopt_new}")
                fopt_sofar = fopt_new
            else:
                print(f"[{self.genename}] Likelihood did not improve after iteration {ITER}.")
                break

        cn_matrix[cn_matrix == MIN_CN] = 0
        self.set_cn_mat(cn_matrix)


    def forward(self, s_idx):
        """
        Runs Forward algorithm on HMM
        """
        # Initialize
        alpha = np.full((self.num_hidden, self.num_obsrvs), -np.inf, dtype=np.float64)

        # Set initial states/probs
        for hs_idx in range(self.num_hidden):
            alpha[hs_idx, 0] = self.init_p[s_idx, hs_idx] + self.emis_p[s_idx, hs_idx, 0]

        # Forward
        for obs_idx in range(1, self.num_obsrvs):
            for curr_hs in range(self.num_hidden):
                a_score = []
                emiss_p = self.emis_p[s_idx, curr_hs, obs_idx]
                for prev_hs in range(self.num_hidden):

                    # Compute score
                    trans_p = self.tran_p[s_idx, obs_idx-1, prev_hs, curr_hs]
                    a_score.append(alpha[prev_hs, obs_idx-1] + trans_p)
                
                alpha[curr_hs, obs_idx] = logsumexp(a_score) + emiss_p

        return alpha
    
    def backward(self, s_idx):
        """
        Runs Backward algorithm on HMM
        """
        # Initialize
        beta = np.full((self.num_hidden, self.num_obsrvs), -np.inf, dtype=np.float64)

        # Set initial states/probs
        beta[:, self.num_obsrvs-1] = 0.0

        # Backward
        for obs_idx in reversed(range(self.num_obsrvs-1)):
            for curr_hs in range(self.num_hidden):
                b_score = []
            
                for prev_hs in range(self.num_hidden):

                    # Compute score
                    emiss_p = self.emis_p[s_idx, prev_hs, obs_idx+1]
                    trans_p = self.tran_p[s_idx, obs_idx, curr_hs, prev_hs]
                    b_score.append(beta[prev_hs, obs_idx+1] + trans_p + emiss_p)
                
                beta[curr_hs, obs_idx] = logsumexp(b_score)
                
        return beta
    
    def forward_backward(self, s_idx):
        """
        alpha : forward probs
        beta  : backward probs
        sigma : forward_sink (sum of last forward probs)
        gamma : (alpha*beta)/sigma or in logspace, alpha+beta-sigma
        """
        alpha = self.forward(s_idx)
        beta  = self.backward(s_idx)
        self.forward_p[s_idx] = alpha
        self.backward_p[s_idx] = beta
        
        mid_obs = self.num_obsrvs//2
        sigma1 = logsumexp(alpha[:, self.num_obsrvs-1])
        sigma2 = logsumexp(beta[:,0] + self.emis_p[s_idx,:,0] + self.init_p[s_idx,:])
        sigma3 = logsumexp(alpha[:,mid_obs] + beta[:,mid_obs])
        
        try:
            assert round(sigma1,4) == round(sigma2,4)
            assert round(sigma2,4) == round(sigma3,4)
        except AssertionError:
            print("[ERROR]", self.genename, s_idx, sigma1, sigma2, sigma3, sep='\t')
            # raise

        gamma = (alpha + beta) - sigma1
        return gamma
    
    def run_forward_backward(self):
        """
        Runs forward and backward algorithm for each sample
        and outputs posterior probability (matrix)
        """
        for s_idx in self.data_order:
            post_prob = self.forward_backward(s_idx)
            self.posterior[s_idx] = post_prob
    
    def compute_path_qual(self, s_idx, cn_path):
        """
        Computes Phred quality value for each CN of a given path
        """
        def compute_qual(ll_arr, best_idx):
            lse_all = logsumexp(ll_arr)
            ll_arr2 = np.delete(ll_arr, best_idx)
            lse_others = logsumexp(ll_arr2)
            quality = -10 * ((lse_others - lse_all)/np.log(10))
            return quality
        
        path_qual = []
        for o_idx in range(self.num_obsrvs):
            ll_arr = np.array(self.posterior)[s_idx,:,o_idx]
            best_cn = cn_path[o_idx]
            exon_qual = abs(round(compute_qual(ll_arr, best_cn), 2))
            #if exon_qual == -0.0:
            #    print(self.genename, s_idx, best_cn, ll_arr.round(), sep='\t')
            path_qual.append(exon_qual)     
        
        return path_qual

    def compute_exon_CN_fwbw(self):
        """
        Compute the most likely path of CNs based on Forward-Backward algorithm
        """
        for s_idx in range(self.num_samples):
            most_likely_path = []
            for o_idx in range(self.num_obsrvs):
                most_likely_cn = np.argmax(self.posterior[s_idx][:, o_idx])
                most_likely_path.append(most_likely_cn)
            self.cn_mat_fb[s_idx] = most_likely_path
            
        self.cn_mat_fb = np.array(self.cn_mat_fb, dtype=np.int32)

