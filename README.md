# bagPipes
***********************

### Overview
***********************
bagPipes is a collection of next-gen sequencing processing pipelines for the Phanstiel Lab.

bagPipes includes pipelines for paired-end RNA-seq data, ChIP-seq/CUT&RUN data, and ATAC-seq data using the Snakemake framework.

This pipeline is inteded to be run on the UNC HPC longleaf cluster with SLURM.

## Running
-----------------------
```bash
sbatch RNApipeCore.sh
```

## Workflow
-----------------------
![](dags/RNApipeCoreDAG.png)

## Unlocking
-----------------------
```bash
./unlock.sh RNApipeCore
```

## To-do
-----------------------
RNApipe
- Add adjustable version #s
- Add merge workflow

ChIPpipe + ATACpipe
- Create

**This pipeline is still in development. More information to come upon completion!**