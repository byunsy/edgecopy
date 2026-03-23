import os
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

# usage: python draw_fractionalCN_histograms.py

df = pd.read_csv('fracCN_vs_wgsCN.tsv', sep='\t')

# ============================================================================
# FCGR3A
# ============================================================================

gene = 'FCGR3A'
df_g = df.loc[df['locus']==gene]

plt.figure(figsize=(6,5), dpi=300)
my_bins = np.arange(1.00, 8.00, 0.05)
ax = sns.histplot(df_g, x='frac_CN', bins=my_bins, hue='True CN', palette='tab10',
                  alpha=1, linewidth=0.1, zorder=10)

ax.grid(which='major', color='#DDDDDD', linewidth=0.8,)
ax.grid(which='minor', color='#DDDDDD', linestyle=':', linewidth=0.5,)
ax.minorticks_on()
ax.legend_.set_title("WGS CN")
ax.set_xlabel(None)
ax.set_ylabel(None)
ax.set_xticks(np.arange(1.5, 6.6, 0.5))

plt.title(gene, fontsize=16)
plt.xlim(1.5, 6.5)
plt.ylim(0, 35)
plt.savefig(f'frac_CN.{gene}.hist.png')
plt.close()

# ============================================================================
# OTOA
# ============================================================================

gene = 'OTOA'
df_g = df.loc[df['locus']==gene]

plt.figure(figsize=(6,5), dpi=300)
my_bins = np.arange(3.00, 8.00, 0.05)
ax = sns.histplot(df_g, x='frac_CN', bins=my_bins, hue='True CN', palette='tab10',
                  alpha=1, linewidth=0.1, zorder=10)

ax.grid(which='major', color='#DDDDDD', linewidth=0.8,)
ax.grid(which='minor', color='#DDDDDD', linestyle=':', linewidth=0.5,)
ax.minorticks_on()
ax.legend_.set_title("WGS CN")
ax.set_xlabel(None)
ax.set_ylabel(None)

plt.title(gene, fontsize=16)
plt.xlim(3.5, 7.5)
plt.ylim(0, 35)
plt.savefig(f'frac_CN.{gene}.hist.png')
plt.close()

# ============================================================================
# RHCE
# ============================================================================

gene = 'RHCE'
df_g = df.loc[df['locus']==gene]

plt.figure(figsize=(6,5), dpi=300)
my_bins = np.arange(1.00, 7.00, 0.05)
ax = sns.histplot(df_g, x='frac_CN', bins=my_bins, hue='True CN', palette='tab10',
                  alpha=1, linewidth=0.1, zorder=10)

ax.grid(which='major', color='#DDDDDD', linewidth=0.8,)
ax.grid(which='minor', color='#DDDDDD', linestyle=':', linewidth=0.5,)
ax.minorticks_on()
ax.legend_.set_title("WGS CN")
ax.set_xlabel(None)
ax.set_ylabel(None)
ax.set_xticks(np.arange(1.5, 5.6, 0.5))

plt.title(gene, fontsize=16)
plt.xlim(1.5, 5.5)
plt.ylim(0, 35)
plt.savefig(f'frac_CN.{gene}.hist.png')
plt.close()

# ============================================================================
# SMN1
# ============================================================================

gene = 'SMN1'
df_g = df.loc[df['locus']==gene]

plt.figure(figsize=(6,5), dpi=300)
my_bins = np.arange(1.00, 7.00, 0.05)
ax = sns.histplot(df_g, x='frac_CN', bins=my_bins, hue='True CN', palette='tab10',
                  alpha=1, linewidth=0.1, zorder=10)

ax.grid(which='major', color='#DDDDDD', linewidth=0.8,)
ax.grid(which='minor', color='#DDDDDD', linestyle=':', linewidth=0.5,)
ax.minorticks_on()
ax.legend_.set_title("WGS CN")
ax.set_xlabel(None)
ax.set_ylabel(None)
ax.set_xticks(np.arange(1.5, 6.5, 0.5))

plt.title(gene, fontsize=16)
plt.xlim(1.5, 6.0)
plt.ylim(0, 35)
plt.savefig(f'frac_CN.{gene}.hist.png')
plt.close()