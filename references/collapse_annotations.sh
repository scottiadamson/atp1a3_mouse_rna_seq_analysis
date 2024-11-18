#for rnaseqc I need to collapse annotations...
#they recommend this script https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py

python collapse_annotation.py gencode.vM35.basic.annotation.gtf gencode.vM35.basic.annotation_rnaseqc_collapsed.gtf

