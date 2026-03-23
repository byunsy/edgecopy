import os
import math
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

# usage: python draw_concordance_gene_subset.py

df = pd.read_csv('edgecopy.concordance.gene_subset.tsv', sep='\t')

plt.figure(figsize=(9, 6), dpi=300)
my_colors = ['#D8315B','#0A2463', '#3E92CC']

ax = sns.barplot(df, x="Gene Type", y="Concordance", hue="Method", 
                 zorder=100, palette=my_colors, saturation=1, width=0.8)

for p in ax.patches:
    height = p.get_height()
    ax.annotate(f'{height:.2f}', (p.get_x() + p.get_width() / 2., height-0.05),
                ha='center', va='bottom', fontsize=11, xytext=(0, 3), zorder=200,
                textcoords='offset points', color='white', fontweight='bold',)

# Formatting
ax.grid(which='major', color='#BBBBBB', linewidth=0.8, zorder=-100)
ax.set_xticks(range(4))
ax.set_xticklabels(["All genes\n(n=108)", 
                    "All genes\n(non-ref)\n(n=108)",
                    "High seq. sim.\ngenes (>0.99)\n(n=40)", 
                    "Common CNV\ngenes\n(n=39)"])

handles, labels = ax.get_legend_handles_labels()

new_labels = [
    "Edgecopy",
    "Edgecopy (No optim.)",
    "ExomeDepth (Custom)"
]

ax.legend(handles, new_labels,
          fontsize=13.7,
          loc='upper center', 
          bbox_to_anchor=(0.5, 1.13), 
          ncol=3,
          handletextpad=0.4,
          columnspacing=0.8,)

plt.grid(axis='x')
plt.ylim(0.35, 1.0)
plt.ylabel("Mean Concordance", fontsize=16)
plt.xlabel(None)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.savefig('concordance.gene_subset.png')