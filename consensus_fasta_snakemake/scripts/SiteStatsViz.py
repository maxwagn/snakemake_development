import sys
import os
jn = os.path.join
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# USAGE python SiteStatsViz.py <input_path_file> 


filename = sys.argv[1]

site_stats_df = pd.read_csv(filename, sep='\t',na_values=['.','NA','NaN'], index_col=[0,1])
site_stats_df = site_stats_df.drop_duplicates(keep='first')
site_stats_df = site_stat_df.sample(frac = 0.01) # ramdomly samples 10% of the dataframe

stats = ['DP','MQ', 'MQ0F','MQSB','NS']
cutoff_quantile = 0.01
n_cols = 3
n_rows = int(np.ceil(len(stats)*1./n_cols))

fig = plt.figure(figsize= (16,3*n_rows) )

for i,stat in enumerate(stats):
    ax = fig.add_subplot(n_rows,n_cols,i+1)
    ss = site_stats_df[stat]
    mx = ss.quantile(1-cutoff_quantile)
    ss.loc[ss>mx] = mx
    ss.plot(kind='hist', ax=ax,bins=100)
    ax.set_xlabel(stat)

file_title = filename.split("/")[-1]
fig.suptitle(file_title, fontsize=16)
fig.savefig(filename[:-7] + '.pdf')

