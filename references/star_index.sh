#!/bin/bash
#SBATCH --job-name=star_index
#SBATCH --partition=pe2
#SBATCH --mail-type=END
#SBATCH --mem=100G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=star_index_%j.log


ref_dir="/gpfs/commons/groups/knowles_lab/sadamson/JAX_Atp1a3_project/references"


module load star/2.7.10b
star --runThreadN 8 --runMode genomeGenerate --genomeDir $ref_dir/star_index --genomeFastaFiles $ref_dir/GRCm39.primary_assembly.genome.fa --sjdbGTFfile $ref_dir/gencode.vM35.basic.annotation.gff3

