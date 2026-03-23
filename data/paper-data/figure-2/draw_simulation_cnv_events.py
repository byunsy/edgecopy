import os
import sys
import math
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
from glob import glob

# usage: python3 draw_simulation_cnv_events.py

df1 = pd.read_csv('edgecopy.simulation.cnv_events.tsv', sep='\t')
df1['Prop'] = df1['Prop'].map(lambda x: f"{x:.2f}")

events = ['Deletion', 'Duplication']
metrics = ['Precision','Recall']

for event in events:
    for metric in metrics:
        df_event = df1.loc[df1['CNV Event'].str.contains(event, na=False)].copy()
        
        plt.figure(figsize=(7, 4), dpi=300)
        ax = sns.lineplot(df_event, x='Prop', y=metric, hue='CNV Event', linewidth=3,
                  marker='o', markersize=10, palette='Set1')

        # Formatting
        ax.grid(which='major', color='#BBBBBB', linewidth=0.8)
        ax.minorticks_on()
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))

        plt.legend(loc='lower left', fontsize=14, title_fontsize=14)
        plt.xlabel('')
        plt.ylabel('')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.ylim(0.994, 1.0003)
        plt.savefig(f'simulation.{event.lower()}.{metric.lower()}.png')
        plt.close()