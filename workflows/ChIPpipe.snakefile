#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
	threads: 8
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

rule peaks:
	input:
		bam = rules.align.output.filteredBam,
		index = rules.align.output.index
	output:
		"output/peaks/{sampleName}_peaks.narrowPeak"
	params:
		dir = "output/peaks",
		version = config['macsVers']
	log:
		err = 'output/logs/peaks_{sampleName}.err',
		out = 'output/logs/peaks_{sampleName}.out'
	benchmark: 
		'output/benchmarks/peaks_{sampleName}.tsv'
	shell:
		"""
		module load macs/{params.version};
		macs2 callpeak -t {input.bam} -f BAM -q 0.01 -g hs --nomodel --shift 0 --extsize 200 --keep-dup all -B --SPMR --outdir {params.dir} -n {wildcards.sampleName}
		"""

rule multiqc:
	input:
		[expand("output/QC/{sampleName}_{read}_fastqc.{ext}", sampleName=key, read=['R1', 'R2'], ext=['zip', 'html']) for key in samples['sn']],
		[expand("output/trim/{sampleName}_{read}_{ext}", sampleName=key, read=['R1', 'R2'], ext=['trimming_report.txt', 'trimmed.fastq.gz']) for key in samples['sn']],
		[expand("output/align/{sampleName}_{ext}", sampleName=key, ext=['nodups_sorted.bam', 'nodups_sorted.bam.bai', 'stats.txt', 'dup_metrics.txt']) for key in samples['sn']],
		[expand("output/signal/{sampleName}.bw", sampleName=key) for key in samples['sn']],
		[expand("output/peaks/{sampleName}_peaks.narrowPeak", sampleName=key) for key in samples['sn']]
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

rule countMatrix:
	input:
		peaks = [expand("output/peaks/{sampleName}_peaks.narrowPeak", sampleName=key) for key in samples['sn']],
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