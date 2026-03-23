import os
import math
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

# usage: python draw_concordance_experimental.py

df = pd.read_csv('edgecopy.concordance.experimental.tsv', sep='\t')

plt.figure(figsize=(5, 5.5), dpi=300)

ax = sns.barplot(df, x="Accuracy", y="Gene", color='#0A2463',
                 zorder=100, saturation=1, width=0.6)
    
for p in ax.patches:
    width = p.get_width()
    y = p.get_y() + p.get_height() / 2

    ax.annotate(
        f'{width:.3f}',
        (width-0.12, y),
        ha='left',
        va='center',
        fontsize=12,
        fontweight='bold',
        xytext=(3, 0),
        textcoords='offset points',
        color='white', zorder=200
    )

# Formatting
ax.spines["right"].set_visible(False)

plt.xlim(0.5, 1.02)
plt.xlabel('Concordance',fontsize=15)
plt.ylabel('')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.tight_layout()
plt.savefig('concordance.experimental.png')