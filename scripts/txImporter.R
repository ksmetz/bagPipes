#!/usr/bin/R

## This script summarizes Salmon quantifications for RNA-Seq count matrix input to DESeq2
## (Note: Ultimately updated to tximeta!)

### INITIALIZE ===================================
## Load required libraries
library(tximport)
library(yaml)
library(GenomicFeatures)

## Parse inputs (1) project name and (2) the path to the relevant transcriptome annotation file
args <- commandArgs(trailingOnly=T)
samplesheet <- (args[1]) # samplesheet.txt
gtfFile <- (args[2]) # /proj/phanstiel_lab/References/genomes/GENCODE.v19/gencode.v19.annotation.gtf_withproteinids
quantDir <- (args[3]) # output/quant/
runName <- (args[4]) # runName from snakemake

## Input errors
if (length(args)!=4) 
{
  stop("Please provide 
        (1) a path to a tsv config file containing all samples used, 
        (2) a transcriptome annotation (GTF) file path, 
        (3) a path to a salmon output directory where the output will be saved, and 
        (4) a base name for the final file")
}


### READ IN ===================================
## Read in samplesheet file, name settings
samplesheet <- read.delim(samplesheet, sep = "\t", header=T, stringsAsFactors=F)
samplesheet <- apply(samplesheet, 2, as.character)
config <- read_yaml("./config/config.yaml")
fileNamesFrom <- config$fileNamesFrom

## Find salmon files 
samples <- apply(samplesheet[, fileNamesFrom], 1, paste, collapse="_")
files <- file.path(quantDir, samples, "quant.sf")
names(files) <- samples


### RUN ===================================
## Make TxDb from GENCODE Annotation file 
txdb <- makeTxDbFromGFF(gtfFile, format="gtf")
k <- keys(txdb, keytype="GENEID")
df <- select(txdb, keys=k, keytype="GENEID", columns="TXNAME")

## Assign proper TxDb columns to tx2gene object
tx2gene <- df[, 2:1]

## Run tximport
txi <- tximport(files, type='salmon', tx2gene=tx2gene)

## Write resulting txi object to RDS file
rdsFile <- paste0(quantDir, "/", runName, "_tximport.rds")
saveRDS(txi, file=rdsFile)