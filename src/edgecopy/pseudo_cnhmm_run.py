import os
import re
import sys
import math
import pandas as pd
import numpy as np
import subprocess

from scipy.special import logsumexp
from collections import namedtuple

from . import pseudo_cnhmm as ch

def run_hmm(inp, data_cc_list):
    """
    Run HMM estimation pipeline on given input samples 
    """

    gene = inp.loci_name
    hmmout_f = os.path.join(inp.finaldir, f'{gene}.agcn.bed')
    hmmmat_f = os.path.join(inp.finaldir, f'{gene}.agcn.tsv')
    hmmout2_f = os.path.join(inp.finaldir, f'{gene}.hmm.out')

    with open(hmmout_f, "w") as hmmout, open(hmmout2_f, "w") as hmmout2:
        
        # Columns for agcn.bed
        hmm_cols = ['#chrom', 'start', 'end', 'locus', 'sample', 'frac', 'point', 'point_qual', 
                    'hmm', 'hmm_filt', 'hmm_qual', 'refset_size', 'alpha_beta']
        print('\t'.join(hmm_cols), file=hmmout)
        
        # Columns for hmm.out
        hmm2_cols = ['locus', 'sample', 'comp', 'comp_size', 'frac_CN', 'grdy_CN', 'grdy_qual',
                     'hmm_CN', 'hmm_CN_filt', 'hmm_qual', 'refset_size', 'alpha_beta'] 
        
        cnt_df = pd.read_csv(inp.gene_cnts_fp, sep='\t')
        cnt_exons = cnt_df.apply(lambda row: f'#{row["#chrom"]}:{row["start"]}-{row["end"]}', axis=1)
        cnt_exons = '\n'.join(cnt_exons.to_list())
        
        print(cnt_exons, file=hmmout2) # gene.hmm.out (exon information)
        print('\t'.join(hmm2_cols), file=hmmout2)
        
        # For each connected component analyzed
        for cc_i, data_cc in enumerate(data_cc_list):
            print(f"- Estimating agCN-HMM paths for {inp.loci_name} component {cc_i} [{data_cc.n} samples].")
            
            # Set up
            hidden_states = range(0,11)
            hmm = ch.cnHMM(hidden_states, 
                           data_cc.ecounts, 
                           data_cc.n, 
                           data_cc.reference_sets, 
                           data_cc.refCN, 
                           data_cc.means,
                           data_cc.alpha,
                           data_cc.beta,
                           data_cc.genename,)
            
            hmm.comb_log = [data_cc.comb_log[i] for i in range(data_cc.n)] ## set scale factor 
            hmm.set_data_order(data_cc)
            
            if inp.debug_mode:
                hmm.set_debug_fp(inp)
                if os.path.exists(inp.debug_fp):
                    os.remove(inp.debug_fp)
            
            hmm.compute_initial_p(inp.priors_fp, data_cc)

            # Round one of HMM
            hmm.compute_transition_p_pseudo(inp.t_prob)
            hmm.compute_exon_CN_pseudo(data_cc, inp.max_iter)
            hmm.run_forward_backward()

            # Round two of HMM
            #hmm.update_transition_p(inp.t_prob)
            #hmm.compute_exon_CN2(data_cc, inp.max_iter)
            #hmm.run_forward_backward()

            # Estimate CNs using Forward-Backward
            hmm.compute_exon_CN_fwbw()
            
            # Reorganize and report summary
            refset_sizes = [len(r) for r in data_cc.reference_sets]
            
            #cnt_df = pd.read_csv(inp.gene_cnts_fp, sep='\t')
            #cnt_exons = cnt_df.apply(lambda row: f'#{row["#chrom"]}:{row["start"]}-{row["end"]}', axis=1)
            #cnt_exons = '\n'.join(cnt_exons.to_list())
            #print(cnt_exons, file=hmmout2) # gene.hmm.out (exon information)

            # (1) gene.agcn.tsv
            # Matrix format
            mat = pd.DataFrame(hmm.cn_mat.T, dtype='Int32')
            mat.columns = data_cc.samplenames
            final_mat = pd.concat([cnt_df.iloc[:,:4], mat], axis=1)
            final_mat.to_csv(hmmmat_f, sep='\t', index=None)
            
            # (2) gene.agcn.bed
            # BED format
            for i in range(hmm.num_obsrvs):
                chrom = cnt_df.iloc[i]['#chrom'] 
                start = cnt_df.iloc[i]['start']
                end   = cnt_df.iloc[i]['end']
                alpha_beta = np.add(hmm.alpha, hmm.beta)

                for j,row in enumerate(hmm.cn_mat_fb.astype(np.int32)):
                    
                    gname = data_cc.genename
                    sname = data_cc.samplenames[j]
                    fract = round(data_cc.fracCN[j], 2)
                    point = hmm.point_cn[j]
                    refsize = refset_sizes[j]
                    
                    pnt_qual = round(data_cc.cn_quality[j]) 
                    hmm_qual = hmm.compute_path_qual(j, row)

                    trust = []
                    for q in hmm_qual:
                        if q < inp.qual_thresh or pnt_qual < inp.qual_thresh:
                            trust.append(False)
                        else:
                            trust.append(True)

                    estim  = row.copy()
                    estim1 = estim[i] 
                    estim2 = estim[i] if trust[i] else '_'

                    hmm_qual_i = round(hmm_qual[i])
                    ab = round(alpha_beta[j], 2)

                    hmm_result = [chrom, start, end, gname, sname, fract, point, pnt_qual, estim1, estim2, hmm_qual_i, refsize, ab]
                    hmm_result = list(map(str, hmm_result))
                    print('\t'.join(hmm_result), file=hmmout)

                    # (3) gene.hmm.out
                    # Only need to run once (when interating i)
                    if i == 0:
                        qual   = ','.join([str(round(q)) for q in hmm_qual])
                        estim3 = ','.join([str(e) for e in estim])
                        estim4 = ','.join([str(e) if trust[k] else '_' for k,e in enumerate(estim)])
                    
                        hmm_result2 = [gname, sname, cc_i, data_cc.n, fract, point, pnt_qual, estim3, estim4, qual, refsize, ab] 
                        hmm_result2 = list(map(str, hmm_result2))
                        print('\t'.join(hmm_result2), file=hmmout2)

    hmmtmp_f = os.path.join(inp.finaldir, f'{gene}.agcn.sorted.bed')
    subprocess.call([f'(grep "^#" {hmmout_f}; grep -v "^#" {hmmout_f} | bedtools sort -i -) > {hmmtmp_f}'], shell=True)
    subprocess.call([f'mv {hmmtmp_f} {hmmout_f}'], shell=True)

