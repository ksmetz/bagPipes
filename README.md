# bagPipes
***********************

### Overview
***********************
bagPipes is a collection of next-gen sequencing processing pipelines for the Phanstiel Lab.

bagPipes includes pipelines for paired-end RNA-seq data, ChIP-seq/CUT&RUN data, and ATAC-seq data using the Snakemake framework.

This pipeline is inteded to be run on the UNC HPC longleaf cluster with SLURM.

## Running RNApipe
-----------------------
The `RNApipe` SLURM wrapper will sequentially launch the `RNApipeCore` and `mergeSignal` workflows using the `RNAconfig.yaml` file.

`RNApipeCore` creates quantification files, alignments and signal tracks. Output files are named according to `fileNamesFrom` list in the config file.

`mergeSignal` generates merged alignments and signal tracks based on the `mergeBy` list in the config file. 

For both workflows, stranded signal tracks will be created by default, but can be toggled off using the `stranded` parameter in the config file.

Launch the pipeline using the following command:
```bash
sbatch RNApipe.sh
```

## Running RNApipe
-----------------------
The `ChIPpipe` SLURM wrapper will launch either the `ChIPpipe` or `ChIPpipeMerged` workflows via the `ChIPpipeLauncher` decision workflow, using the `ChIPconfig.yaml` file.

`ChIPpipe` creates peak calls, alignments and signal tracks. Output files are named according to `fileNamesFrom` list in the config file.

`ChIPpipeMerged` creates merged alignments and signal tracks, and calls peaks from the merged alignments. This will be run if the `mergeBy` entry of the config is not blank, and different from the `fileNamesFrom` entry.

Both workflows generate a final summary `peakCounts.tsv` file containing the counts from each individual alignment at a merged set of peaks.

Launch the pipeline using the following command:
```bash
sbatch ChIPpipe.sh
```

## Workflows
-----------------------
# RNApipe
![](dags/RNApipeCoreDAG.png)

# ChIPpipe
![](dags/ChIPpipeDAG.png)

## Merging
-----------------------
Follow these steps to create merged signal tracks and alignment files. Requires existing `.bam` and `.bai` files for individual samples from a core run. Determines merging using the `mergeBy` parameter in `config/config.yaml`.  

## Benchmarking
-----------------------
Run the `benchmarking.py` script within `workflows/utils` in order to summarize all existing benchmark files.
```
module load python/3.6.6
python3 ./workflows/utils/benchmarking.py
```

## Unlocking
-----------------------
Use this command when you want to delete a previously generated output and re-run.
```bash
./unlock.sh RNApipe
./unlock.sh ChIPpipe
```

## To-do
-----------------------
**RNApipe**
- Make `fastqc` rule more robust to different file extensions (currently hardcoded to accept `.fastq.gz` only; new param in config?)
- Make creation of the samplesheet a rule
- (Optional) Make output for quant rule (Salmon) accept the directories themselves as outputs
- Consider trying to incorporate the mergeSignal into RNApipeCore 

**ChIPpipe + ATACpipe**
- Create ATACpipe based on ChIPpipe
- Consider editing ChIPpipe rule countMatrix to create header line from `{wildcards.sampleName}`
- Move MACS2 peak calling settings into config

**This pipeline is still in development. More information to come upon completion!**