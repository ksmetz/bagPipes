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

## Set sample names
samples['sn'] = samples[config['fileNamesFrom']].apply('_'.join, axis=1)

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
		[expand("output/QC/{sampleName}_{read}_fastqc.{ext}", sampleName=key, read=['R1', 'R2'], ext=['zip', 'html']) for key in samples['sn']],
		# [expand("output/trim/{sampleName}_{read}_{ext}", sampleName=key, read=['R1', 'R2'], ext=['trimming_report.txt', 'trimmed.fastq.gz']) for key in samples['sn']]
		[expand("output/quant/{sampleName}/quant.sf", sampleName=key) for key in samples['sn']],
		# [expand("output/align/{sampleName}_{ext}", sampleName=key, ext=['sorted.bam', 'sorted.bam.bai', 'stats.txt']) for key in samples['sn']]
		[expand("output/signal/unstranded/{sampleName}.bw", sampleName=key) for key in samples['sn']],
		[expand("output/signal/stranded/{sampleName}_{dir}.bw", sampleName=key, dir=['fwd', 'rev']) for key in samples['sn']]

rule fastqc:
	input:
		R1 = lambda wildcards: read1.get(wildcards.sampleName),
		R2 = lambda wildcards: read2.get(wildcards.sampleName)
	output:
		zip1 = "output/QC/{sampleName}_R1_fastqc.zip", # temp
		zip2 = "output/QC/{sampleName}_R2_fastqc.zip", # temp
		html1 = "output/QC/{sampleName}_R1_fastqc.html", # temp
		html2 = "output/QC/{sampleName}_R2_fastqc.html" # temp
	log:
		err = 'output/logs/fastqc_{sampleName}.err',
		out = 'output/logs/fastqc_{sampleName}.out'
	params:
		dir = "output/QC"
	shell:
		"""
		module load fastqc/0.11.5;
		fastqc -o {params.dir} {input.R1} {input.R2} 1> {log.out} 2> {log.err};
		mv {params.dir}/$(basename {input.R1} .fastq.gz)_fastqc.zip  {output.zip1};
		mv {params.dir}/$(basename {input.R2} .fastq.gz)_fastqc.zip  {output.zip2};
		mv {params.dir}/$(basename {input.R1} .fastq.gz)_fastqc.html  {output.html1};
		mv {params.dir}/$(basename {input.R2} .fastq.gz)_fastqc.html  {output.html2}
		"""

rule trim:
	input:
		R1 = lambda wildcards: read1.get(wildcards.sampleName),
		R2 = lambda wildcards: read2.get(wildcards.sampleName)
	output:
		report1 = "output/trim/{sampleName}_R1_trimming_report.txt", # temp
		report2 = "output/trim/{sampleName}_R2_trimming_report.txt", # temp
		trim1 = "output/trim/{sampleName}_R1_trimmed.fastq.gz", # temp
		trim2 = "output/trim/{sampleName}_R2_trimmed.fastq.gz" # temp
	log:
		err = 'output/logs/trim_{sampleName}.err',
		out = 'output/logs/trim_{sampleName}.out'
	params:
		dir = "output/trim"
	shell:
		"""
		module load trim_galore/0.4.3;
		trim_galore -o {params.dir} --paired {input.R1} {input.R2} 1> {log.out} 2> {log.err};
		mv {params.dir}/$(basename {input.R1} .fastq.gz)_val_1.fq.gz  {output.trim1};
		mv {params.dir}/$(basename {input.R2} .fastq.gz)_val_2.fq.gz  {output.trim2};
		mv {params.dir}/$(basename {input.R1})_trimming_report.txt  {output.report1};
		mv {params.dir}/$(basename {input.R2})_trimming_report.txt  {output.report2}
		"""

rule quant:
	input:
		trim1 = rules.trim.output.trim1,
		trim2 = rules.trim.output.trim2
	output:
		quant = "output/quant/{sampleName}/quant.sf"
	log:
		err = 'output/logs/quant_{sampleName}.err',
		out = 'output/logs/quant_{sampleName}.out'
	params:
		dir = "output/quant",
		index = config['salmon']
	shell:
		"""
		module load salmon/1.4.0;
		salmon quant --writeUnmappedNames --threads 1 -i {params.index} -l A -1 {input.trim1} -2 {input.trim2} -o {params.dir}/{wildcards.sampleName} 1> {log.out} 2> {log.err};
		"""

rule align:
	input:
		trim1 = rules.trim.output.trim1,
		trim2 = rules.trim.output.trim2
	output:
		index = "output/align/{sampleName}_sorted.bam.bai",
		stats = "output/align/{sampleName}_stats.txt",
		bam = "output/align/{sampleName}_sorted.bam"
	log:
		err = 'output/logs/align_{sampleName}.err',
		out = 'output/logs/align_{sampleName}.out'
	params:
		index = config['hisat2']
	shell:
		"""
		module load hisat2/2.1.0;
		module load samtools/1.9;
		hisat2 -q -x {params.index} -1 {input.trim1} -2 {input.trim2} | samtools view -u | samtools sort -o {output.bam} 1> {log.out} 2> {log.err};
		samtools flagstat {output.bam} > {output.stats} 2>> {log.err};
		samtools index {output.bam} 1>> {log.out} 2>> {log.err}
		"""

rule signal:
	input:
		bam = rules.align.output.bam
	output:
		"output/signal/unstranded/{sampleName}.bw"
	log:
		err = 'output/logs/signal_{sampleName}.err',
		out = 'output/logs/signal_{sampleName}.out'
	shell:
		"""
		module load deeptools/3.0.1;
		bamCoverage -b {input.bam} -o {output} 1> {log.out} 2> {log.err}
		"""

rule strandedSignal:
	input:
		bam = rules.align.output.bam
	output:
		fwd = "output/signal/stranded/{sampleName}_fwd.bw",
		rev = "output/signal/stranded/{sampleName}_rev.bw"
	log:
		err = 'output/logs/strandedSignal_{sampleName}.err',
		out = 'output/logs/strandedSignal_{sampleName}.out'
	shell:
		"""
		module load deeptools/3.0.1;
		bamCoverage --filterRNAstrand forward -b {input.bam} -o {output.fwd} 1> {log.out} 2> {log.err};
		bamCoverage --filterRNAstrand reverse -b {input.bam} -o {output.rev} 1>> {log.out} 2>> {log.err}
		"""

