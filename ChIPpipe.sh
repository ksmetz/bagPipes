#! /bin/bash -login
#SBATCH -J ChIPpipe
#SBATCH -t 4320
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1  
#SBATCH -p general
#SBATCH --mem=2gb
#SBATCH -o "%x-%j.out"

## Exit if any command fails
set -e

## Load required modules
module load python/3.6.6

## Create and activate virtual environment with requirements
python3 -m venv env && source env/bin/activate && pip3 install -r config/requirements.txt

## Make directory for slurm logs
mkdir -p output/logs_slurm

## Execute ChIPpipeCore snakemake workflow
snakemake -s workflows/ChIPpipeLauncher.snakefile --configfile "config/ChIPconfig.yaml" --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status ./workflows/utils/status.py -j 100 --max-jobs-per-second 5 --max-status-checks-per-second 5 --rerun-incomplete -p --latency-wait 500 

## Success message
echo "Workflow completed successfully!"

