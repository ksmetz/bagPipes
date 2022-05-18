#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import shutil

##### Load config and sample sheets #####
configfile: "config/config.yaml"

## Read in samplesheet
samples = pd.read_table(config["samplesheet"])

## Convert all columns to strings
samples = samples.astype(str)

## Concatenate the sequencing directory to Read1 and Read2 for full paths
samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

## Concatenate columns to identify which groups to run (i.e. Seq_Rep will be run together)
samples['id'] = samples[config['mergeBy']].agg('_'.join, axis=1)

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
        "output/test.txt"

rule test:
	input:
		read1 = lambda wildcards: samples.loc[wildcards.sample]["Read1"],
		read2 = lambda wildcards: samples.loc[wildcards.sample]["Read2"]
	output:
		test = "output/{sample}.txt" # temp
	params:
		dir = "output/",
		name = '{sample}'
	shell:
        'echo {input.read1} {input.read2} >> {params.dir}/{params.name}.txt'
