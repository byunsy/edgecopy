import os
import sys
import pandas as pd
import numpy as np
from glob import glob
from scipy import stats

POP = sys.argv[1]
CENTERS = ['BCM','BGI','BI','WUGSC']

col = ['psv','gene','chrom','start','nr_correct','nr_total','nr_ccrd','r_correct','r_total','r_ccrd']
MIN_SAMPLES = 1

low_ccrds, high_ccrds = [], []
r_low_ccrds, nr_low_ccrds = [], []
r_high_ccrds, nr_high_ccrds = [], []
print()

all_df_low, all_df_high = [], []

for CTR in CENTERS:
    
    # Create dictionary {gene:refcn,...}
    refcns = {}
    refcn_fp = f'../{POP}-{CTR}-agcn/summary.input.info.tsv'
    with open(refcn_fp, 'r') as file:
        for line in file:
            info = line.strip().split('\t')
            refcns[info[0]] = int(info[1])
    
    fp = f'./psv.concord.{POP}.{CTR}.out'
    df = pd.read_csv(fp, sep=' ', header=None)
    df.columns = col
    df['refcn'] = df.gene.apply(lambda x: refcns.get(x.split('.')[0], 4))
    df = df.loc[df.psv == 'PSV'].copy()

    # Low RefCN ==============================================================
    df_low = df.loc[(df.nr_total>MIN_SAMPLES) & (df.r_total>MIN_SAMPLES) & (df.refcn<=4)].copy()
    r_ccrd, nr_ccrd = df_low[['r_ccrd','nr_ccrd']].mean()
    
    df_low['correct'] = df_low['nr_correct'] + df_low['r_correct']
    df_low['total']   = df_low['nr_total'] + df_low['r_total']
    df_low['ccrd']    = df_low['correct'] / df_low['total']

    comb_ccrd = df_low['ccrd'].mean()
    print(fp, 'RefCN=4', round(r_ccrd,4), round(nr_ccrd,4), round(comb_ccrd,4), sep='\t')
    
    low_ccrds.append(comb_ccrd)
    r_low_ccrds.append(r_ccrd)
    nr_low_ccrds.append(nr_ccrd)

    all_df_low.append(df_low)

    # High RefCN =============================================================
    df_high = df.loc[(df.nr_total>MIN_SAMPLES) & (df.r_total>MIN_SAMPLES) & (df.refcn>4) & (df.refcn<8)].copy()
    r_ccrd, nr_ccrd = df_high[['r_ccrd','nr_ccrd']].mean()
    
    df_high['correct'] = df_high['nr_correct'] + df_high['r_correct']
    df_high['total']   = df_high['nr_total'] + df_high['r_total']
    df_high['ccrd']    = df_high['correct'] / df_high['total']

    comb_ccrd = df_high['ccrd'].mean()
    print(fp, 'RefCN=6', round(r_ccrd,4), round(nr_ccrd,4), round(comb_ccrd,4), sep='\t')
    
    high_ccrds.append(comb_ccrd)
    r_high_ccrds.append(r_ccrd)
    nr_high_ccrds.append(nr_ccrd)
    print()
    
    all_df_high.append(df_high)


print('='*60)
print('TOTAL RefCN=4', np.array(r_low_ccrds).mean().round(4), 
                       np.array(nr_low_ccrds).mean().round(4),
                       np.array(low_ccrds).mean().round(4), sep='\t')
print('TOTAL RefCN=6', np.array(r_high_ccrds).mean().round(4), 
                       np.array(nr_high_ccrds).mean().round(4),
                       np.array(high_ccrds).mean().round(4), sep='\t')

final_low  = pd.concat(all_df_low, axis=0).reset_index(drop=True)
final_high = pd.concat(all_df_high, axis=0).reset_index(drop=True)

print(final_low)
print(final_high)

print('='*60)
print('COMBINED CONCORDANCE')
print('TOTAL RefCN=4', #round(final_low.r_ccrd.mean(),4), 
                       #round(final_low.nr_ccrd.mean(),4),
                       round(final_low.ccrd.mean(),4),
                       round(stats.sem(final_low.ccrd, nan_policy='omit'),4), sep='\t')
print('TOTAL RefCN=6', #round(final_high.r_ccrd.mean(),4), 
                       #round(final_high.nr_ccrd.mean(),4),
                       round(final_high.ccrd.mean(),4),
                       round(stats.sem(final_high.ccrd, nan_policy='omit'),4), sep='\t')


