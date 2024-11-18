#!/bin/bash
#SBATCH --job-name=kallisto_index
#SBATCH --partition=pe2
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=kallisto_index_%j.log

source ~/.bashrc
mamba activate ALS_analysis_env

base_dir="/gpfs/commons/groups/knowles_lab/sadamson/JAX_Atp1a3_project/references"
kallisto index -i $base_dir/kallisto.index $base_dir/gencode.vM35.transcripts.fa 

