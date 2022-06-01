# bagPipes
***********************

### Overview
***********************
bagPipes is a collection of next-gen sequencing processing pipelines for the Phanstiel Lab.

bagPipes includes pipelines for paired-end RNA-seq data, ChIP-seq/CUT&RUN data, and ATAC-seq data using the Snakemake framework.

This pipeline is inteded to be run on the UNC HPC longleaf cluster with SLURM.

## Running
-----------------------
The `RNApipe` SLURM wrapper will sequentially launch the `RNApipeCore` and `mergeSignal` workflows using the `RNAconfig.yaml` file.

`RNApipeCore` creates quantification files, alignments and signal tracks. Output files are named according to `fileNamesFrom` list in the config file.

`mergeSignal` generates merged alignments and signal tracks based on the `mergeBy` list in the config file. 

For both workflows, stranded signal tracks will be created by default, but can be toggled off using the `stranded` parameter in the config file.

Launch the pipeline using the following command:
```bash
sbatch RNApipe.sh
```

## Workflow
-----------------------
![](dags/RNApipeCoreDAG.png)

## Merging
-----------------------
Follow these steps to create merged signal tracks and alignment files. Requires existing `.bam` and `.bai` files for individual samples from a core run. Determines merging using the `mergeBy` parameter in `config/config.yaml`.  

## Unlocking
-----------------------
Use this command when you want to delete a previously generated output and re-run.
```bash
./unlock.sh RNApipeCore
```

## To-do
-----------------------
**RNApipe**
- Make `fastqc` rule more robust to different file extensions (currently hardcoded to accept `.fastq.gz` only)
- Make creation of the samplesheet a rule
- (Optional) Make output for quant rule (Salmon) accept the directories themselves as outputs
- Consider trying to incorporate the mergeSignal into RNApipeCore 

**ChIPpipe + ATACpipe**
- Create single pipeline for each, compatible with signalMerge script
- Include logic for using either merged- or non-merged aligns for peak calling, depending on `fileNamesFrom`/`mergeBy` params

**This pipeline is still in development. More information to come upon completion!**