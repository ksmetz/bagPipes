#!/usr/bin/env python3

import pandas as pd
import os

##### Load config and sample sheets #####

configfile: "config.yaml"

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)

##### Define rules #####

rule all:
	input:
		'peaks/peakReadCounts.tsv',
		'peaks/FRiP.txt'

rule trim_trimGalore_PE:
	input:
		read1 = lambda wildcards: samples.loc[wildcards.sample]["Read1"],
		read2 = lambda wildcards: samples.loc[wildcards.sample]["Read2"]
	output:
		read1 = "trim/{sample}_R1_val_1.fq.gz", # temp
		read2 = "trim/{sample}_R2_val_2.fq.gz"  # temp
	params:
		dir = "trim/",
		name = '{sample}'
	shell:
		'module load trim_galore/0.6.2; '
		'trim_galore --gzip -o {params.dir} --basename {params.name} --paired {input.read1} {input.read2}; '

rule align_bwa_PE:
	input:
		trimread1 = "trim/{sample}_R1_val_1.fq.gz",
		trimread2 = "trim/{sample}_R2_val_2.fq.gz"
	output:
		sortedbam = "align/{sample}_sorted.bam" # temp
	params:
		index = config['bwa_index']
	threads: 8
	shell:
		'module load bwa/0.7.17; '
		'module load samtools/1.8; '
		'bwa mem -t {threads} {params.index} {input.trimread1} {input.trimread2} | samtools view -u | samtools sort -o {output}'

rule stats_samtools_PE:
	input:
		"align/{sample}_sorted.bam"
	output:
		"align/{sample}_stats.txt"
	threads: 8
	shell:
		'module load samtools/1.8; '
		'samtools flagstat -@{threads} {input} > {output}'

rule filter_picardtools_PE:
	input:
		"align/{sample}_sorted.bam"
	output:
		filteredbam = "align/{sample}_filtersorted.bam",
		metrics = "align/{sample}_dup_metrics.txt"
	threads: 8
	shell:
		'module load java/10.0.2; '
		'java -Xmx16g -jar /nas/longleaf/apps/picard/2.10.3/picard-2.10.3/picard.jar '
		'MarkDuplicates I={input} O={output.filteredbam} M={output.metrics} '
		'REMOVE_SEQUENCING_DUPLICATES=true'

rule index_samtools_PE:
	input:
		"align/{sample}_filtersorted.bam"
	output:
		"align/{sample}_filtersorted.bam.bai"
	threads: 8
	shell:
		'module load samtools/1.8; '
		'samtools index {input}'

rule signal_bedgraphs_PEunmerged:
	input:
		'align/{sample}_filtersorted.bam'
	output:
		'signal/{sample}.bedgraph'
	threads: 8
	shell:
		'module load bedtools/2.26; '
		'bedtools genomecov -bga -ibam {input} > {output}'

rule peakcall_macs2_PEunmerged:
	input:
		'align/{sample}_filtersorted.bam'
	output:
		file = 'peaks/{sample}_peaks.narrowPeak'
	params:
		name = '{sample}'
	threads: 8
	shell:
		'module load macs/2016-02-15; '
		'macs2 callpeak -t {input} -f BAM -g hs -B --outdir peaks -n {params.name}'

rule FRiPpart1_samtools_PE:
	input:
		peaks = "peaks/{sample}_peaks.narrowPeak",
		alignments = "align/{sample}_filtersorted.bam"		
	output:
		inPeaks = temp('peaks/{sample}_peaks.txt'),
		inTotal = temp('peaks/{sample}_total.txt')
	threads: 8
	shell:
		'module load samtools/1.8; '
		'samtools view -F 260 -c -L {input.peaks} {input.alignments} >> {output.inPeaks}; '
		'samtools view -F 260 -c {input.alignments} >> {output.inTotal}'

rule FRiPpart2_samtools_PE:
	input:
		inPeaks = 'peaks/{sample}_peaks.txt',
		inTotal = 'peaks/{sample}_total.txt'
	output:
		temp('peaks/{sample}_FRiP.txt')
	params:
		name = '{sample}'
	shell: 
		'cat {input.inPeaks} {input.inTotal} | tr "\n" "\t" | awk -v sample="{params.name}" \'{{OFS="\\t"}}; {{print sample, $1, $2, $1/$2}}\' >> {output}'

rule FRiPpart3_samtools_PE:
	input:
		expand('peaks/{sample}_FRiP.txt', sample=samples.index)
	output:
		'peaks/FRiP.txt'
	shell:
		'cat {input} | awk -v OFS="\\t" \'BEGIN{{print "sample", "inPeaks", "inTotal", "fraction"}}; {{print $0}}\' > {output}'

rule peakmerge_bedtools_PE:
	input:
		expand("peaks/{sample}_peaks.narrowPeak", sample=samples.index)
	output:
		'peaks/peakMerge.narrowPeak'
	threads: 8
	shell:
		'module load bedtools/2.26; '
		'cat {input} | awk \'{{ OFS="\\t" }};{{ print $1, $2, $3, $4 }}\' | sort -k1,1 -k2,2n | bedtools merge > {output}'

rule count_bedtools_PE:
	input:
		bam = expand("align/{sample}_filtersorted.bam", sample=samples.index),
		bai = expand("align/{sample}_filtersorted.bam.bai", sample=samples.index),
		bed = 'peaks/peakMerge.narrowPeak',
	output:
		initial = temp('peaks/temp_counts.tsv'),
		final = 'peaks/peakReadCounts.tsv'
	params:
		header = expand('{sample}', sample=samples.index)
	threads: 8
	shell:
		'module load bedtools/2.26; '
		'bedtools multicov -bams {input.bam} -bed {input.bed} > {output.initial}; '
		'awk -v OFS="\\t" \'BEGIN{{print "chr", "start", "end", "{params.header}" }}; {{print}}\' {output.initial} > {output.final}'
