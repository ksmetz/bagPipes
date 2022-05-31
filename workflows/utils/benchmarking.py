import os
import glob
import pandas as pd
from snakemake.io import *

## Find + combine benchmark files
dfs = list()
tables = list()
files = glob.glob(os.path.join('output/benchmarks/*.tsv'))
for f in files:
	table = pd.read_csv(f, sep="\t")
	id = os.path.splitext(os.path.basename(f))[0]
	table['rule'] = id.split("_")[0]
	table['sample'] = "_".join(id.split("_")[1:])
	tables.append(table[['sample', 'rule', 's', 'max_rss', 'max_vms', 'max_uss', 'max_pss', 'io_in', 'io_out', 'mean_load']])
dfs.append(pd.concat(tables))

df = pd.concat(dfs)

## Function to convert seconds to hh:mm:ss
def convert(seconds): 
	min, sec = divmod(seconds, 60) 
	hour, min = divmod(min, 60) 
	return "%d:%02d:%02d" % (hour, min, sec)

## Aggregate statistics (various ways)
agg = df.groupby('rule').agg('max')
agg = df.groupby(['rule', 'sample']).agg('max') ## For some reason this looks nicer

## Convert seconds to hh:mm:ss
agg['s'] = agg['s'].apply(convert)

## Write out as tsv
agg.to_csv("benchmarking.tsv", sep="\t")