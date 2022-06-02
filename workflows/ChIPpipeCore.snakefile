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

## Optional + conditional outputs
## if mergeBy supplied, call peaks from merged samples
mergeOutputs = dict()
peaksOutput = dict()
peaksInput = dict()

if config['mergeBy'] == config['fileNamesFrom'] or config['mergeBy'] == '':
	mergeOutputs["align"] = []
	mergeOutputs["signal"] = []
	peaksInput['bam'] = [expand("output/align/{sampleName}_nodups_sorted.bam", sampleName=key) for key in samples['sn']]
	peaksInput['index'] = [expand("output/align/{sampleName}_nodups_sorted.bam.bai", sampleName=key) for key in samples['sn']]
	peaksOutput["peaks"] = [expand("output/peaks/{sampleName}_peaks.narrowPeak", sampleName=key) for key in samples['sn']]
	peakName = [expand(key) for key in samples['sn']]
else:
	## Merge according to mergeBy parameter, define merge name (mn)
	samples['mn'] = samples[config['mergeBy']].agg('_'.join, axis=1)

	## Build dictionary of merged BAM files
	mergeSamples = samples.groupby('mn')['sn'].apply(list).to_dict()

	## Define rule all inputs
	mergeOutputs["align"] = [expand("output/mergeAlign/{mergeName}_{ext}", mergeName=key, ext=['nodups_sorted.bam', 'nodups_sorted.bam.bai', 'stats.txt']) for key in mergeSamples]
	mergeOutputs["signal"] = [expand("output/mergeSignal/{mergeName}.bw", mergeName=key) for key in mergeSamples]
	peaksInput['bam'] = [expand("output/mergeAlign/{mergeName}_nodups_sorted.bam", mergeName=key) for key in mergeSamples]
	peaksInput['index'] = [expand("output/mergeAlign/{mergeName}_nodups_sorted.bam.bai", mergeName=key) for key in mergeSamples]
	peaksOutput["peaks"] = [expand("output/peaks/{mergeName}_peaks.narrowPeak", mergeName=key) for key in mergeSamples]
	peakName = [expand(key) for key in mergeSamples]

## Define actions on success
onsuccess:
	## Success message
	print("ChIPpipe completed successfully! Wahoo!")

##### Define rules #####
rule all:
	input:
		("output/QC/{name}_multiqc_report.html").format(name=runName),
		("output/peaks/{name}_peakCounts.tsv").format(name=runName)

rule fastqc:
	input:
		R1 = lambda wildcards: read1.get(wildcards.sampleName),
		R2 = lambda wildcards: read2.get(wildcards.sampleName)
	output:
		zip1 = temp("output/QC/{sampleName}_R1_fastqc.zip"), # temp
		zip2 = temp("output/QC/{sampleName}_R2_fastqc.zip"), # temp
		html1 = temp("output/QC/{sampleName}_R1_fastqc.html"), # temp
		html2 = temp("output/QC/{sampleName}_R2_fastqc.html") # temp
	log:
		err = 'output/logs/fastqc_{sampleName}.err',
		out = 'output/logs/fastqc_{sampleName}.out'
	params:
		dir = "output/QC",
		version = config['fastqcVers']
	benchmark: 
		'output/benchmarks/fastqc_{sampleName}.tsv'
	shell:
		"""
		module load fastqc/{params.version};
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
		report1 = temp("output/trim/{sampleName}_R1_trimming_report.txt"), # temp
		report2 = temp("output/trim/{sampleName}_R2_trimming_report.txt"), # temp
		trim1 = temp("output/trim/{sampleName}_R1_trimmed.fastq.gz"), # temp
		trim2 = temp("output/trim/{sampleName}_R2_trimmed.fastq.gz") # temp
	log:
		err = 'output/logs/trim_{sampleName}.err',
		out = 'output/logs/trim_{sampleName}.out'
	benchmark: 
		'output/benchmarks/trim_{sampleName}.tsv'
	params:
		dir = "output/trim",
		version = config['trimVers']
	shell:
		"""
		module load trim_galore/{params.version};
		trim_galore -o {params.dir} --paired {input.R1} {input.R2} 1> {log.out} 2> {log.err};
		mv {params.dir}/$(basename {input.R1} .fastq.gz)_val_1.fq.gz  {output.trim1};
		mv {params.dir}/$(basename {input.R2} .fastq.gz)_val_2.fq.gz  {output.trim2};
		mv {params.dir}/$(basename {input.R1})_trimming_report.txt  {output.report1};
		mv {params.dir}/$(basename {input.R2})_trimming_report.txt  {output.report2}
		"""

rule align:
	input:
		trim1 = rules.trim.output.trim1,
		trim2 = rules.trim.output.trim2
	output:
		bam = temp("output/align/{sampleName}_sorted.bam"),
		stats = "output/align/{sampleName}_stats.txt",
		filteredBam = "output/align/{sampleName}_nodups_sorted.bam",
		dupStats = "output/align/{sampleName}_dup_metrics.txt",
		index = "output/align/{sampleName}_nodups_sorted.bam.bai"
	log:
		err = 'output/logs/align_{sampleName}.err',
		out = 'output/logs/align_{sampleName}.out'
	benchmark: 
		'output/benchmarks/align_{sampleName}.tsv'
	params:
		index = config['bwamem'],
		bwaVersion = config['bwaVers'],
		samtoolsVersion = config['samtoolsVers'],
		javaVersion = config['javaVers']
	shell:
		"""
		module load bwa/{params.bwaVersion};
		module load samtools/{params.samtoolsVersion};
		module load java/{params.javaVersion};
		bwa mem -t 8 {params.index} {input.trim1} {input.trim2} | samtools view -u | samtools sort -o {output.bam} 1> {log.out} 2> {log.err};
		samtools flagstat {output.bam} > {output.stats} 2>> {log.err};
		java -Xmx16g -jar /nas/longleaf/apps/picard/2.10.3/picard-2.10.3/picard.jar MarkDuplicates I={output.bam} O={output.filteredBam} M={output.dupStats} REMOVE_SEQUENCING_DUPLICATES=true;
		samtools index {output.filteredBam} 1>> {log.out} 2>> {log.err}
		"""

rule signal:
	input:
		bam = rules.align.output.filteredBam
	output:
		"output/signal/{sampleName}.bw"
	log:
		err = 'output/logs/signal_{sampleName}.err',
		out = 'output/logs/signal_{sampleName}.out'
	benchmark: 
		'output/benchmarks/signal_{sampleName}.tsv'
	params:
		version = config['deeptoolsVers']
	shell:
		"""
		module load deeptools/{params.version};
		bamCoverage -b {input.bam} -o {output} 1> {log.out} 2> {log.err}
		"""

rule mergeAlign:
	input:
		bams = lambda wildcards: ["output/align/{sampleName}_nodups_sorted.bam".format(sampleName=value) for value in mergeSamples[wildcards.mergeName]],
		bais = lambda wildcards: ["output/align/{sampleName}_nodups_sorted.bam.bai".format(sampleName=value) for value in mergeSamples[wildcards.mergeName]]
	output:
		bam = "output/mergeAlign/{mergeName}_nodups_sorted.bam",
		bai = "output/mergeAlign/{mergeName}_nodups_sorted.bam.bai",
		stats = "output/mergeAlign/{mergeName}_stats.txt"
	log:
		err = 'output/logs/mergeAlign_{mergeName}.err',
		out = 'output/logs/mergeAlign_{mergeName}.out'
	benchmark: 
		'output/benchmarks/mergeAlign_{mergeName}.tsv'
	params:
		version = config['samtoolsVers']
	shell:
		"""
		module load samtools/{params.version};
		samtools merge {output.bam} {input.bams} 1>> {log.out} 2>> {log.err};
		samtools flagstat {output.bam} > {output.stats} 2>> {log.err};
		samtools index {output.bam} 1>> {log.out} 2>> {log.err}
		"""

rule mergeSignal:
	input:
		bam = rules.mergeAlign.output.bam
	output:
		"output/mergeSignal/{mergeName}.bw"
	log:
		err = 'output/logs/mergeSignal_{mergeName}.err',
		out = 'output/logs/mergeSignal_{mergeName}.out'
	benchmark: 
		'output/benchmarks/mergeSignal_{mergeName}.tsv'
	params:
		version = config['deeptoolsVers']
	shell:
		"""
		module load deeptools/{params.version};
		bamCoverage -b {input.bam} -o {output} 1> {log.out} 2> {log.err}
		"""

## NEEDS FIX: uses full list of peak names as output for each macs2 call, causing error
rule peaks:
	input:
		bam = peaksInput["bam"],
		index = peaksInput["index"]
	output:
		peaksOutput["peaks"]
	params:
		dir = "output/peaks",
		version = config['macsVers'],
		name = peakName
	log:
		err = 'output/logs/peaks_{params.name}.err',
		out = 'output/logs/peaks_{params.name}.out'
	benchmark: 
		'output/benchmarks/peaks_{params.name}.tsv'
	shell:
		"""
		module load macs/{params.version};
		macs2 callpeak -t {input.bam} -f BAM -q 0.01 -g hs --nomodel --shift 0 --extsize 200 --keep-dup all -B --SPMR --outdir {params.dir} -n {params.name}
		"""

rule multiqc:
	input:
		[expand("output/QC/{sampleName}_{read}_fastqc.{ext}", sampleName=key, read=['R1', 'R2'], ext=['zip', 'html']) for key in samples['sn']],
		[expand("output/trim/{sampleName}_{read}_{ext}", sampleName=key, read=['R1', 'R2'], ext=['trimming_report.txt', 'trimmed.fastq.gz']) for key in samples['sn']],
		[expand("output/align/{sampleName}_{ext}", sampleName=key, ext=['nodups_sorted.bam', 'nodups_sorted.bam.bai', 'stats.txt', 'dup_metrics.txt']) for key in samples['sn']],
		[expand("output/signal/{sampleName}.bw", sampleName=key) for key in samples['sn']],
		peaksOutput["peaks"],
		mergeOutputs["align"],
		mergeOutputs["signal"]
	output:
		("output/QC/{name}_multiqc_report.html").format(name=runName)
	params:
		name = runName,
		version = config['multiqcVers']
	log:
		err = 'output/logs/multiqc.err',
		out = 'output/logs/multiqc.out'
	benchmark: 
		'output/benchmarks/multiqc.tsv'
	shell:
		"""
		module load multiqc/{params.version};
		multiqc -f output/* -o output/QC 1> {log.out} 2> {log.err}
		mv output/QC/multiqc_report.html output/QC/{params.name}_multiqc_report.html
		mv output/QC/multiqc_data output/QC/{params.name}_multiqc_data
		"""

# Try adjusting header line to for NAME in {wildcards.sampleName}
rule countMatrix:
	input:
		peaks = peaksOutput["peaks"],
		other = [expand("output/align/{sampleName}_{ext}", sampleName=key, ext=['nodups_sorted.bam.bai', 'stats.txt', 'dup_metrics.txt']) for key in samples['sn']],
		bams = [expand("output/align/{sampleName}_nodups_sorted.bam", sampleName=key) for key in samples['sn']]
	output:
		tmp = temp("output/peaks/TEMP.bed"),
		peakMerge = temp(("output/peaks/{name}_mergedPeaks.bed").format(name=runName)),
		peakCounts = ("output/peaks/{name}_peakCounts.tsv").format(name=runName)
	params:
		samplesheet = newSamplesheet,
		name = runName,
		version = config['bedtoolsVers']
	log:
		err = 'output/logs/countMatrix.err',
		out = 'output/logs/countMatrix.out'
	benchmark: 
		'output/benchmarks/countMatrix.tsv'
	shell:
		"""
		module load bedtools/{params.version}
		cat {input.peaks} | awk '{{ OFS="\t"}};{{ print $1, $2, $3, $4 }}' | sort -k1,1 -k2,2n | bedtools merge > {output.tmp};
		grep -ve "-1" {output.tmp} > {output.peakMerge}; 
		printf "chr\tstart\tstop\t" > {output.peakCounts};
		for f in {input.bams}; do NAME=$(basename $f _filter_sorted.bam); printf '%s\t' "$NAME" >> {output.peakCounts}; done;
		printf "\n" >> {output.peakCounts};
		bedtools multicov -bams {input.bams} -bed {output.peakMerge} >> {output.peakCounts}
		"""