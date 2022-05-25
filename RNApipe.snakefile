#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil
import os

##### Load config and sample sheets #####
configfile: "config/config.yaml"

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"], sep='\t')

## Convert all columns to strings
samples = samples.astype(str)

## Concatenate the sequencing directory to Read1 and Read2 for full paths
samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

## Concatenate columns to identify which groups to run (i.e. Seq_Rep will be run together)
samples['id'] = samples[config['mergeBy']].apply('_'.join, axis=1)

## Group by id and extract Read1 & Read2
read1 = samples.groupby('id')['Read1'].apply(list).to_dict()
read2 = samples.groupby('id')['Read2'].apply(list).to_dict()

## Define actions on success
onsuccess:
    ## Success message
    print("RNApipe completed successfully! Wahoo!")

##### Define rules #####
rule all:
	input:
		[expand("output/{group}.txt", group=key) for key in samples['id']]

rule test:
	input:
		read1 = lambda wildcards: read1.get(wildcards.group),
		read2 = lambda wildcards: read2.get(wildcards.group)
	output:
		test = "output/{group}.txt" # temp
	params:
		dir = "output/"
	shell:
		'echo {input.read1} {input.read2} >> {params.dir}/{wildcards.group}.txt'
