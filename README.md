This github repo has the code to analyze RNA-seq from the following paper:  
M. Terrey et al., (2024) Alternating hemiplegia of childhood associated mutations in *Atp1a3* reveal diverse neurological alterations in mice. Manuscript in preparation.   
Note that there is some differential splicing analysis in this repo which is not part of the publication currently.
 
In order to reproduce this analysis first clone this repo:  
`git clone https://github.com/scottiadamson/atp1a3_mouse_rna_seq_analysis.git
cd atp1a3_mouse_rna_seq_analysis`  
Then install the required dependecies using conda:  
`mamba env create -f environment.yml`
Additional QC info requires [https://github.com/getzlab/rnaseqc](RNA-SeQC2) which is best installed with pre-compiled binary files [https://github.com/getzlab/rnaseqc/releases](here).  
  
  
For RNA-seq analysis, install R version 4.3.3 and the following dependencies with CRAN:
`install.packages('tidyverse', 'cowplot', 'UpSetR', 'BiocManager')`
Then install other R dependencies with BiocManager:
`BiocManager::install('DESeq2')
BiocManager::install('tximport')
BiocManager::install('EnhancedVolcano')`
  
Install SRA toolkit to download the fastq files using the instructions [https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit](here).  
Then download the fastq files:
`fastq-dump -X 2 ./fastq/SRAXXXXXXXX`

Then download the references and build references:
`cd references
bash download_annotations.sh
bash kallisto_index.sh
bash star_index.sh # only if doing splicing workflow
cd ../`

Then adjust the variables in the config/config.yaml file to reflect the paths in your setup.  

Then run the snakemake workflow:
`conda activate atp1a3_mouse_rna_seq_env

slurm_config_folder="/config/slurm_config"
snakemake --snakefile Snakefile --profile $slurm_config_folder --nolock -k --rerun-triggers mtime
`

After this has run successfully perform the differential expression analysis:
`cd differential_expression
Rscript DESeq.R` 

If performing differential splicing analysis, first set "splicing_workflow" in the config file to True. 
Then install [http://davidaknowles.github.io/leafcutter/](leafcutter). 
Then change the paths in the leafcutter.sh file and run differential splicing analysis:
`cd differential_splicing
bash leafcutter.sh`
