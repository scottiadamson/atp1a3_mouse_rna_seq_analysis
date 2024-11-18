#!/bin/bash
#SBATCH --job-name=JAX_Atp1a3_leafcutter
#SBATCH --partition=pe2
#SBATCH --mail-type=END
#SBATCH --mem=15G
#SBATCH --time=12:00:00
#SBATCH --output=JAX_Atp1a3_leafcutter%j.log
#SBATCH --cpus-per-task=4

module load leafcutter/knowles

base_dir="/gpfs/commons/groups/knowles_lab/sadamson/JAX_Atp1a3_project/differential_splicing"
gtf="/gpfs/commons/groups/knowles_lab/sadamson/JAX_Atp1a3_project/references/gencode.vM35.basic.annotation_collapsed.gtf"
leafcutter_dir="/nfs/sw/leafcutter/leafcutter-0.2.9/leafcutter"

python $leafcutter_dir/clustering/leafcutter_cluster_regtools.py -j $base_dir/junction_samples.txt -m 50 -o ./ -l 500000

mkdir $base_dir/BS_D801N_WT
$leafcutter_dir/scripts/leafcutter_ds.R -s 20 --min_samples_per_intron 3 -e $gtf -o $base_dir/BS_D801N_WT/leafcutter_ds --num_threads 4 $base_dir/_perind_numers.counts.gz $base_dir/BS_D801N_WT_metadata.txt
mkdir $base_dir/BS_E815K_WT
$leafcutter_dir/scripts/leafcutter_ds.R -s 20 --min_samples_per_intron 3 -e $gtf -o $base_dir/BS_E815K_WT/leafcutter_ds --num_threads 4 $base_dir/_perind_numers.counts.gz $base_dir/BS_E815K_WT_metadata.txt

mkdir $base_dir/CM_D801N_WT
$leafcutter_dir/scripts/leafcutter_ds.R -s 20 --min_samples_per_intron 3 -e $gtf -o $base_dir/CM_D801N_WT/leafcutter_ds --num_threads 4 $base_dir/_perind_numers.counts.gz $base_dir/CM_D801N_WT_metadata.txt
mkdir $base_dir/CM_E815K_WT
$leafcutter_dir/scripts/leafcutter_ds.R -s 20 --min_samples_per_intron 3 -e $gtf -o $base_dir/CM_E815K_WT/leafcutter_ds --num_threads 4 $base_dir/_perind_numers.counts.gz $base_dir/CM_E815K_WT_metadata.txt

mkdir $base_dir/CX_D801N_WT
$leafcutter_dir/scripts/leafcutter_ds.R -s 20 --min_samples_per_intron 3 -e $gtf -o $base_dir/CX_D801N_WT/leafcutter_ds --num_threads 4 $base_dir/_perind_numers.counts.gz $base_dir/CX_D801N_WT_metadata.txt
mkdir $base_dir/CX_E815K_WT
$leafcutter_dir/scripts/leafcutter_ds.R -s 20 --min_samples_per_intron 3 -e $gtf -o $base_dir/CX_E815K_WT/leafcutter_ds --num_threads 4 $base_dir/_perind_numers.counts.gz $base_dir/CX_E815K_WT_metadata.txt

mkdir $base_dir/HP_D801N_WT
$leafcutter_dir/scripts/leafcutter_ds.R -s 20 --min_samples_per_intron 3 -e $gtf -o $base_dir/HP_D801N_WT/leafcutter_ds --num_threads 4 $base_dir/_perind_numers.counts.gz $base_dir/HP_D801N_WT_metadata.txt
mkdir $base_dir/HP_E815K_WT
$leafcutter_dir/scripts/leafcutter_ds.R -s 20 --min_samples_per_intron 3 -e $gtf -o $base_dir/HP_E815K_WT/leafcutter_ds --num_threads 4 $base_dir/_perind_numers.counts.gz $base_dir/HP_E815K_WT_metadata.txt

