#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import glob
import shutil
import os
from utils.namer import namer

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep='\t')

## Convert all columns to strings
samples = samples.astype(str)

## Decide between paired-end or single-end sequencing workflow
if config['paired'] == 'TRUE':

	## Concatenate the sequencing directory to Read1 and Read2 for full paths
	samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
	samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

	## Set sample names
	samples['sn'] = samples[config['fileNamesFrom']].apply('_'.join, axis=1)

	## Set run summary name using helper script
	runName = namer(samples, config['fileNamesFrom'])

	## Write samplesheet with sampleNames to new file
	newSamplesheet = ('output/{name}_ChIPpipeSamplesheet.txt').format(name = runName)
	samples.to_csv(newSamplesheet, sep="\t", index=False)

	## Group by id and extract Read1 & Read2
	read1 = samples.groupby('sn')['Read1'].apply(list).to_dict()
	read2 = samples.groupby('sn')['Read2'].apply(list).to_dict()

	## Decide between workflow with merging and without, based on config file
	if config['mergeBy'] == config['fileNamesFrom'] or config['mergeBy'] == '':
		include: "ChIPpipe.snakefile"
	else:
		## Merge according to mergeBy parameter, define merge name (mn)
		samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)

		## Build dictionary of merged BAM files
		mergeSamples = samples.groupby('mn')['sn'].apply(list).to_dict()
	
		include: "ChIPpipeMerged.snakefile"

else:

	## Concatenate the sequencing directory to Read1
	samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)

	## Set sample names
	samples['sn'] = samples[config['fileNamesFrom']].apply('_'.join, axis=1)

	## Set run summary name using helper script
	runName = namer(samples, config['fileNamesFrom'])

	## Write samplesheet with sampleNames to new file
	newSamplesheet = ('output/{name}_ChIPpipeSamplesheet.txt').format(name = runName)
	samples.to_csv(newSamplesheet, sep="\t", index=False)

	## Group by id and extract Read1
	read1 = samples.groupby('sn')['Read1'].apply(list).to_dict()

	## Decide between workflow with merging and without, based on config file
	if config['mergeBy'] == config['fileNamesFrom'] or config['mergeBy'] == '':
		include: "ChIPpipeSE.snakefile"

	else:
		## Merge according to mergeBy parameter, define merge name (mn)
		samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)

		## Build dictionary of merged BAM files
		mergeSamples = samples.groupby('mn')['sn'].apply(list).to_dict()
	
		include: "ChIPpipeMergedSE.snakefile"


## Define actions on success
onsuccess:
	## Success message
	print("ChIPpipe completed successfully! Wahoo!")

##### Define rules #####
rule all:
	input:
		("output/QC/{name}_multiqc_report.html").format(name=runName),
		("output/peaks/{name}_peakCounts.tsv").format(name=runName)