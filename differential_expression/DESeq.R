#module load R/4.3.3

library('DESeq2')
library('tximport')
library('tidyverse')
library('cowplot')
library('EnhancedVolcano')
library('UpSetR')

add_gene_names <- function(res_df, transcript_info, gene_biotypes){
  res_df$gene_id <- rownames(res_df)
  res_df <- left_join(data.frame(res_df), transcript_info %>% dplyr::select(gene_id, gene_name) %>% distinct, by = c('gene_id' = 'gene_id'))
  res_df <- res_df[order(res_df$pvalue),]
  return(res_df)
}

#to load the workspace
base_dir <- '/gpfs/commons/groups/knowles_lab/sadamson/JAX_Atp1a3_project/'
#load(paste0(base_dir, 'differential_expression/DEseq2.RData'))

theme_set(theme_bw() + theme(text = element_text(family = 'Helvetica')))

#read in metadata
metadata_df <- read_csv(paste0(base_dir, 'metadata.csv'))

#find kallisto files
files <- file.path(base_dir, "kallisto_out/", metadata_df$Sample, "abundance.h5")

#verify that all of the files exist
all(file.exists(files))

# parse gene and transcript info from sample tsv
sample_df <- read_tsv(files[1] %>% str_replace('.h5', '.tsv'))
transcript_info <- separate_wider_delim(data.frame(sample_df$target_id), delim = "|",
                                       cols = colnames(data.frame(sample_df$target_id)),
                                       names = c('transcript_id','gene_id', 'affy_gene_id', 'affy_transcript_id',
                                                 'transcript_name', 'gene_name', 'transcript_length', 'transcript_type', 'filler'))
transcript_info <- transcript_info %>% dplyr::select(transcript_id, gene_id, transcript_name, gene_name, transcript_length, transcript_type)
transcript_info$complete_transcript_name <- sample_df$target_id
tx2gene <- transcript_info %>% dplyr::select(complete_transcript_name, gene_id, transcript_length)
gene_id2gene <- as.data.frame(transcript_info %>% dplyr::select(gene_id, gene_name) %>% unique())
gene_id2gene_list <- gene_id2gene$gene_name
names(gene_id2gene_list) <- gene_id2gene$gene_id

#format experimental covariates
colData <- data.frame(metadata_df)
rownames(colData) <- metadata_df$Sample
keeper_cols <- c('Tissue', 'Genotype', 'Sex')
colData <- colData[,keeper_cols]
colData$Genotype <- relevel(as.factor(colData$Genotype), ref = "WT")

#import data into DESeq2
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, dropInfReps = TRUE, lengthCol = 'transcript_length')
count_df <- as.data.frame(txi$counts)
colnames(count_df) <- rownames(colData)
count_df$gene_name <- as.vector(gene_id2gene_list[rownames(count_df)])
count_df <- count_df[,c('gene_name', rownames(colData))]
write.table(count_df, paste0(base_dir, 'differential_expression/gene_count_matrix.tsv'), quote = FALSE, sep = '\t', row.names = FALSE)

dds <- DESeqDataSetFromTximport(txi,
                               colData = colData,
                               design = ~ Genotype + Tissue + Sex)

#output for GEO here
sample_nickname_tib <- read_tsv(paste0(base_dir, 'sample_nicknames.tsv'))
nicknames <- sample_nickname_tib$nickname
names(nicknames) <- sample_nickname_tib$original_name
gene_counts <- as.data.frame(counts(dds))
colnames(gene_counts) <- nicknames[colnames(gene_counts)]
original_cols <- colnames(gene_counts)
gene_counts$gene_id <- rownames(gene_counts)
gene_counts$gene_name <- gene_id2gene_list[gene_counts$gene_id]
gene_counts <- gene_counts[,c('gene_id', 'gene_name', original_cols)]
as_tibble(gene_counts) %>% write_tsv(paste0(base_dir, 'for_GEO/Gene_counts.tsv'))

#PCA of all RNA-seq samples
PCA_df_prep <- function(deseq_object, metatdata_df){
  vsd <- vst(deseq_object, blind=TRUE)
  #perform PCA
  PCs <- prcomp(t(assay(vsd)),
                center = TRUE,
                scale. = FALSE)
  n_PCs <- 9
  important_PCs <- as_tibble(PCs$x[,paste0('PC', 1:n_PCs)])
  important_PCs$sample <- rownames(PCs$x)
  percent_var <- PCs$sdev^2/sum( PCs$sdev^2 )
  var_explained_df <- data.frame(percent_var*100, paste0('PC', 1:length(percent_var)))
  colnames(var_explained_df) <- c('percent_variance_explained', 'PC')
  
  PCA_plot_df <- left_join(important_PCs, metadata_df, by = c('sample' = 'Sample'))
  output_list <- list()
  output_list$var_explained_df <- var_explained_df
  output_list$PCA_plot_df <- PCA_plot_df
  return(output_list)
}

PCA_stuff <- PCA_df_prep(dds, metatdata_df)
var_explained_df <- PCA_stuff$var_explained_df
PCA_plot_df <- PCA_stuff$PCA_plot_df

#plot and save PCA
ggplot(PCA_plot_df, aes(x = PC1, y = PC2, color = Genotype, shape = factor(Tissue))) +
  geom_point(size = 2) +
  scale_color_manual(values=c("blue", "red", 'black')) +
  labs(x = paste0("PC1:   ", round(var_explained_df[1,'percent_variance_explained'],2), '% variance explained'),
       y = paste0("PC2:   ", round(var_explained_df[2,'percent_variance_explained'],2), '% variance explained'),
       color = "Genotype", shape = 'Tissue')

ggsave(paste0(base_dir, 'differential_expression/plots/', 'all_sample_PCA.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'all_sample_PCA.pdf'), width = 8, height = 6)

##Analyze Brainstem
BS_samples <- pull(metadata_df %>% filter(Tissue == 'Brainstem') %>% dplyr::select(Sample))
BS_files <- file.path(base_dir, "kallisto_out/", BS_samples, "abundance.h5")
BS_txi <- tximport(BS_files, type = "kallisto", tx2gene = tx2gene, dropInfReps = TRUE, lengthCol = 'transcript_length')
BS_colData <- colData[colData$Tissue == 'Brainstem', c('Genotype', 'Sex')]

BS_dds <- DESeqDataSetFromTximport(BS_txi,
                                    colData = BS_colData,
                                    design = ~ Genotype + Sex)
BS_dds_expressed <- rowSums(counts(BS_dds)>=1) >= 0.5 * nrow(BS_colData)
BS_dds <- BS_dds[BS_dds_expressed,]
BS_dds <- DESeq(BS_dds, fitType='local')

BS_PCA_stuff <- PCA_df_prep(BS_dds, metatdata_df)
BS_var_explained_df <- BS_PCA_stuff$var_explained_df
BS_PCA_plot_df <- BS_PCA_stuff$PC

ggplot(BS_PCA_plot_df, aes(x = PC1, y = PC2, label = sample, color = Genotype, shape = Sex)) +
 geom_point(size = 2) +
 scale_color_manual(values=c("blue", "red", 'black')) +
 labs(x = paste0("PC1:   ", round(BS_var_explained_df[1,'percent_variance_explained'],2), '% variance explained'),
      y = paste0("PC2:   ", round(BS_var_explained_df[2,'percent_variance_explained'],2), '% variance explained'),
      color = "Molecular model", shape = 'Sex')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_PCA.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_PCA.pdf'), width = 8, height = 6)

BS_D801N_WT_res <- lfcShrink(BS_dds, contrast = c('Genotype', 'D801N', 'WT'), type = 'ashr')
BS_D801N_WT_res <- add_gene_names(BS_D801N_WT_res, transcript_info)
write_tsv(BS_D801N_WT_res, paste0(base_dir, 'differential_expression/', 'BS_D801N_vs_WT.tsv'))

EnhancedVolcano(BS_D801N_WT_res,
                lab = BS_D801N_WT_res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = 'Brainstem D801N vs WT')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_D801N_vs_WT_volcano_plot.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_D801N_vs_WT_volcano_plot.pdf'), width = 8, height = 6)

ggplot(BS_D801N_WT_res, aes(pvalue)) + geom_histogram(bins = 25) + theme_bw()
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_D801N_vs_WT_pvalue_histogram.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_D801N_vs_WT_pvalue_histogram.pdf'), width = 8, height = 6)

BS_E815K_WT_res <- lfcShrink(BS_dds, contrast = c('Genotype', 'E815K', 'WT'), type = 'ashr')
BS_E815K_WT_res <- add_gene_names(BS_E815K_WT_res, transcript_info)
write_tsv(BS_E815K_WT_res, paste0(base_dir, 'differential_expression/', 'BS_E815K_vs_WT.tsv'))

EnhancedVolcano(BS_E815K_WT_res,
                lab = BS_E815K_WT_res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Brainstem E815K vs WT')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_E815K_vs_WT_volcano_plot.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_E815K_vs_WT_volcano_plot.pdf'), width = 8, height = 6)

ggplot(BS_E815K_WT_res, aes(pvalue)) + geom_histogram(bins = 25) + theme_bw()
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_E815K_vs_WT_pvalue_histogram.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_E815K_vs_WT_pvalue_histogram.pdf'), width = 8, height = 6)

BS_D801N_E815K_res <- lfcShrink(BS_dds, contrast = c('Genotype', 'D801N', 'E815K'), type = 'ashr')
BS_D801N_E815K_res <- add_gene_names(BS_D801N_E815K_res, transcript_info)
write_tsv(BS_D801N_E815K_res, paste0(base_dir, 'differential_expression/', 'BS_D801N_vs_E815K.tsv'))

EnhancedVolcano(BS_D801N_E815K_res,
                lab = BS_D801N_E815K_res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Brainstem D801N vs E815K')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_D801N_vs_E815K_volcano_plot.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_D801N_vs_E815K_volcano_plot.pdf'), width = 8, height = 6)

ggplot(BS_D801N_E815K_res, aes(pvalue)) + geom_histogram(bins = 25) + theme_bw()
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_D801N_vs_E815K_pvalue_histogram.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'BS_D801N_vs_E815K_pvalue_histogram.pdf'), width = 8, height = 6)


##Analyze Cerebellum
CM_samples <- pull(metadata_df %>% filter(Tissue == 'Cerebellum') %>% dplyr::select(Sample))
CM_files <- file.path(base_dir, "kallisto_out/", CM_samples, "abundance.h5")
CM_txi <- tximport(CM_files, type = "kallisto", tx2gene = tx2gene, dropInfReps = TRUE, lengthCol = 'transcript_length')
CM_colData <- colData[colData$Tissue == 'Cerebellum', c('Genotype', 'Sex')]

CM_dds <- DESeqDataSetFromTximport(CM_txi,
                                   colData = CM_colData,
                                   design = ~ Genotype + Sex)
CM_dds_expressed <- rowSums(counts(CM_dds)>=1) >= 0.5 * nrow(CM_colData)
CM_dds <- CM_dds[CM_dds_expressed,]
CM_dds <- DESeq(CM_dds, fitType='local')
plotDispEsts(CM_dds)

CM_PCA_stuff <- PCA_df_prep(CM_dds, metatdata_df)
CM_var_explained_df <- CM_PCA_stuff$var_explained_df
CM_PCA_plot_df <- CM_PCA_stuff$PC

ggplot(CM_PCA_plot_df, aes(x = PC1, y = PC2, label = sample, color = Genotype, shape = Sex)) +
  geom_point(size = 2) +
  scale_color_manual(values=c("blue", "red", 'black')) +
  labs(x = paste0("PC1:   ", round(CM_var_explained_df[1,'percent_variance_explained'],2), '% variance explained'),
       y = paste0("PC2:   ", round(CM_var_explained_df[2,'percent_variance_explained'],2), '% variance explained'),
       color = "Molecular model", shape = 'Sex')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_PCA.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_PCA.pdf'), width = 8, height = 6)

CM_D801N_WT_res <- lfcShrink(CM_dds, contrast = c('Genotype', 'D801N', 'WT'), type = 'ashr')
CM_D801N_WT_res <- add_gene_names(CM_D801N_WT_res, transcript_info)
write_tsv(CM_D801N_WT_res, paste0(base_dir, 'differential_expression/', 'CM_D801N_vs_WT.tsv'))

EnhancedVolcano(CM_D801N_WT_res,
                lab = CM_D801N_WT_res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Cerebellum D801N vs WT')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_D801N_vs_WT_volcano_plot.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_D801N_vs_WT_volcano_plot.pdf'), width = 8, height = 6)

ggplot(CM_D801N_WT_res, aes(pvalue)) + geom_histogram(bins = 25) + theme_bw()
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_D801N_vs_WT_pvalue_histogram.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_D801N_vs_WT_pvalue_histogram.pdf'), width = 8, height = 6)

CM_E815K_WT_res <- lfcShrink(CM_dds, contrast = c('Genotype', 'E815K', 'WT'), type = 'ashr')
CM_E815K_WT_res <- add_gene_names(CM_E815K_WT_res, transcript_info)
write_tsv(CM_E815K_WT_res, paste0(base_dir, 'differential_expression/', 'CM_E815K_vs_WT.tsv'))

EnhancedVolcano(CM_E815K_WT_res,
                lab = CM_E815K_WT_res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Cerebellum E815K vs WT')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_E815K_vs_WT_volcano_plot.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_E815K_vs_WT_volcano_plot.pdf'), width = 8, height = 6)

ggplot(CM_E815K_WT_res, aes(pvalue)) + geom_histogram(bins = 25) + theme_bw()
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_E815K_vs_WT_pvalue_histogram.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_E815K_vs_WT_pvalue_histogram.pdf'), width = 8, height = 6)

CM_D801N_E815K_res <- lfcShrink(CM_dds, contrast = c('Genotype', 'D801N', 'E815K'), type = 'ashr')
CM_D801N_E815K_res <- add_gene_names(CM_D801N_E815K_res, transcript_info)
write_tsv(CM_D801N_E815K_res, paste0(base_dir, 'differential_expression/', 'CM_D801N_vs_E815K.tsv'))

EnhancedVolcano(CM_D801N_E815K_res,
                lab = CM_D801N_E815K_res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Cerebellum D801N vs E815K')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_D801N_vs_E815K_volcano_plot.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_D801N_vs_E815K_volcano_plot.pdf'), width = 8, height = 6)

ggplot(CM_D801N_E815K_res, aes(pvalue)) + geom_histogram(bins = 25) + theme_bw()
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_D801N_vs_E815K_pvalue_histogram.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CM_D801N_vs_E815K_pvalue_histogram.pdf'), width = 8, height = 6)


##Analyze Cortex
CX_samples <- pull(metadata_df %>% filter(Tissue == 'Cortex') %>% dplyr::select(Sample))
CX_files <- file.path(base_dir, "kallisto_out/", CX_samples, "abundance.h5")
CX_txi <- tximport(CX_files, type = "kallisto", tx2gene = tx2gene, dropInfReps = TRUE, lengthCol = 'transcript_length')
CX_colData <- colData[colData$Tissue == 'Cortex', c('Genotype', 'Sex')]

CX_dds <- DESeqDataSetFromTximport(CX_txi,
                                   colData = CX_colData,
                                   design = ~ Genotype + Sex)
CX_dds_expressed <- rowSums(counts(CX_dds)>=1) >= 0.5 * nrow(CX_colData)
CX_dds <- CX_dds[CX_dds_expressed,]
CX_dds <- DESeq(CX_dds, fitType='local')
plotDispEsts(CX_dds)

CX_PCA_stuff <- PCA_df_prep(CX_dds, metatdata_df)
CX_var_explained_df <- CX_PCA_stuff$var_explained_df
CX_PCA_plot_df <- CX_PCA_stuff$PC

ggplot(CX_PCA_plot_df, aes(x = PC1, y = PC2, label = sample, color = Genotype, shape = Sex)) +
  geom_point(size = 2) +
  scale_color_manual(values=c("blue", "red", 'black')) +
  labs(x = paste0("PC1:   ", round(CX_var_explained_df[1,'percent_variance_explained'],2), '% variance explained'),
       y = paste0("PC2:   ", round(CX_var_explained_df[2,'percent_variance_explained'],2), '% variance explained'),
       color = "Molecular model", shape = 'Sex')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_PCA.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_PCA.pdf'), width = 8, height = 6)

CX_D801N_WT_res <- lfcShrink(CX_dds, contrast = c('Genotype', 'D801N', 'WT'), type = 'ashr')
CX_D801N_WT_res <- add_gene_names(CX_D801N_WT_res, transcript_info)
write_tsv(CX_D801N_WT_res, paste0(base_dir, 'differential_expression/', 'CX_D801N_vs_WT.tsv'))

EnhancedVolcano(CX_D801N_WT_res,
                lab = CX_D801N_WT_res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = 'Cortex D801N vs WT')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_D801N_vs_WT_volcano_plot.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_D801N_vs_WT_volcano_plot.pdf'), width = 8, height = 6)

ggplot(CX_D801N_WT_res, aes(pvalue)) + geom_histogram(bins = 25) + theme_bw()
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_D801N_vs_WT_pvalue_histogram.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_D801N_vs_WT_pvalue_histogram.pdf'), width = 8, height = 6)

CX_E815K_WT_res <- lfcShrink(CX_dds, contrast = c('Genotype', 'E815K', 'WT'), type = 'ashr')
CX_E815K_WT_res <- add_gene_names(CX_E815K_WT_res, transcript_info)
write_tsv(CX_E815K_WT_res, paste0(base_dir, 'differential_expression/', 'CX_E815K_vs_WT.tsv'))

EnhancedVolcano(CX_E815K_WT_res,
                lab = CX_E815K_WT_res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = 'Cortex E815K vs WT')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_E815K_vs_WT_volcano_plot.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_E815K_vs_WT_volcano_plot.pdf'), width = 8, height = 6)

ggplot(CX_E815K_WT_res, aes(pvalue)) + geom_histogram(bins = 25) + theme_bw()
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_E815K_vs_WT_pvalue_histogram.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_E815K_vs_WT_pvalue_histogram.pdf'), width = 8, height = 6)

CX_D801N_E815K_res <- lfcShrink(CX_dds, contrast = c('Genotype', 'D801N', 'E815K'), type = 'ashr')
CX_D801N_E815K_res <- add_gene_names(CX_D801N_E815K_res, transcript_info)
write_tsv(CX_D801N_E815K_res, paste0(base_dir, 'differential_expression/', 'CX_D801N_vs_E815K.tsv'))

EnhancedVolcano(CX_D801N_E815K_res,
                lab = CX_D801N_E815K_res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = 'Cortex D801N vs E815K')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_D801N_vs_E815K_volcano_plot.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_D801N_vs_E815K_volcano_plot.pdf'), width = 8, height = 6)

ggplot(CX_D801N_E815K_res, aes(pvalue)) + geom_histogram(bins = 25) + theme_bw()
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_D801N_vs_E815K_pvalue_histogram.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'CX_D801N_vs_E815K_pvalue_histogram.pdf'), width = 8, height = 6)


##Analyze Hippocampus
HP_samples <- pull(metadata_df %>% filter(Tissue == 'Hippocampus') %>% dplyr::select(Sample))
HP_files <- file.path(base_dir, "kallisto_out/", HP_samples, "abundance.h5")
HP_txi <- tximport(HP_files, type = "kallisto", tx2gene = tx2gene, dropInfReps = TRUE, lengthCol = 'transcript_length')
HP_colData <- colData[colData$Tissue == 'Hippocampus', c('Genotype', 'Sex')]

HP_dds <- DESeqDataSetFromTximport(HP_txi,
                                   colData = HP_colData,
                                   design = ~ Genotype + Sex)
HP_dds_expressed <- rowSums(counts(HP_dds)>=1) >= 0.5 * nrow(HP_colData)
HP_dds <- HP_dds[HP_dds_expressed,]
HP_dds <- DESeq(HP_dds, fitType='local')
plotDispEsts(HP_dds)

HP_PCA_stuff <- PCA_df_prep(HP_dds, metatdata_df)
HP_var_explained_df <- HP_PCA_stuff$var_explained_df
HP_PCA_plot_df <- HP_PCA_stuff$PC

ggplot(HP_PCA_plot_df, aes(x = PC1, y = PC2, label = sample, color = Genotype, shape = Sex)) +
  geom_point(size = 2) +
  scale_color_manual(values=c("blue", "red", 'black')) +
  labs(x = paste0("PC1:   ", round(HP_var_explained_df[1,'percent_variance_explained'],2), '% variance explained'),
       y = paste0("PC2:   ", round(HP_var_explained_df[2,'percent_variance_explained'],2), '% variance explained'),
       color = "Molecular model", shape = 'Sex')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_PCA.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_PCA.pdf'), width = 8, height = 6)

HP_D801N_WT_res <- lfcShrink(HP_dds, contrast = c('Genotype', 'D801N', 'WT'), type = 'ashr')
HP_D801N_WT_res <- add_gene_names(HP_D801N_WT_res, transcript_info)
write_tsv(HP_D801N_WT_res, paste0(base_dir, 'differential_expression/', 'HP_D801N_vs_WT.tsv'))

EnhancedVolcano(HP_D801N_WT_res,
                lab = HP_D801N_WT_res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = 'Hippocampus D801N vs WT')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_D801N_vs_WT_volcano_plot.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_D801N_vs_WT_volcano_plot.pdf'), width = 8, height = 6)

ggplot(HP_D801N_WT_res, aes(pvalue)) + geom_histogram(bins = 25) + theme_bw()
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_D801N_vs_WT_pvalue_histogram.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_D801N_vs_WT_pvalue_histogram.pdf'), width = 8, height = 6)

HP_E815K_WT_res <- lfcShrink(HP_dds, contrast = c('Genotype', 'E815K', 'WT'), type = 'ashr')
HP_E815K_WT_res <- add_gene_names(HP_E815K_WT_res, transcript_info)
write_tsv(HP_E815K_WT_res, paste0(base_dir, 'differential_expression/', 'HP_E815K_vs_WT.tsv'))

EnhancedVolcano(HP_E815K_WT_res,
                lab = HP_E815K_WT_res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = 'Hippocampus E815K vs WT')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_E815K_vs_WT_volcano_plot.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_E815K_vs_WT_volcano_plot.pdf'), width = 8, height = 6)

ggplot(HP_E815K_WT_res, aes(pvalue)) + geom_histogram(bins = 25) + theme_bw()
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_E815K_vs_WT_pvalue_histogram.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_E815K_vs_WT_pvalue_histogram.pdf'), width = 8, height = 6)

HP_D801N_E815K_res <- lfcShrink(HP_dds, contrast = c('Genotype', 'D801N', 'E815K'), type = 'ashr')
HP_D801N_E815K_res <- add_gene_names(HP_D801N_E815K_res, transcript_info)
write_tsv(HP_D801N_E815K_res, paste0(base_dir, 'differential_expression/', 'HP_D801N_vs_E815K.tsv'))

EnhancedVolcano(HP_D801N_E815K_res,
                lab = HP_D801N_E815K_res$gene_name,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = 'Hippocampus D801N vs E815K')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_D801N_vs_E815K_volcano_plot.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_D801N_vs_E815K_volcano_plot.pdf'), width = 8, height = 6)

ggplot(HP_D801N_E815K_res, aes(pvalue)) + geom_histogram(bins = 25) + theme_bw()
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_D801N_vs_E815K_pvalue_histogram.svg'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'HP_D801N_vs_E815K_pvalue_histogram.pdf'), width = 8, height = 6)


#make a matrix of DE genes in each comparison
all_comps <- list(BS_D801N_WT_res, CM_D801N_WT_res, CX_D801N_WT_res, HP_D801N_WT_res,
               BS_E815K_WT_res, CM_E815K_WT_res, CX_E815K_WT_res, HP_E815K_WT_res,
               BS_D801N_E815K_res, CM_D801N_E815K_res, CX_D801N_E815K_res, HP_D801N_E815K_res)
extract_DE_genes <- function(df){
  DE_genes <- df[df[,'padj'] <= 0.05, 'gene_id']
  DE_genes <- DE_genes[is.na(DE_genes) == FALSE]
  return(DE_genes)
}

DE_lookup <- function(df, gene_list){
  rownames(df) <- df[,'gene_id']
  df <- df[gene_list,]
  is_sig <- df[,'padj'] <= 0.05
  return(is_sig)
}

DE_genes <- unique(unlist(lapply(all_comps, extract_DE_genes)))
DE_sig <- lapply(all_comps, DE_lookup, gene_list = DE_genes)
comp_vector <- c('Brainstem_D801N_vs_WT', 'Cerebellum_D801N_vs_WT', 'Cortex_D801N_vs_WT', 'Hippocampus_D801N_vs_WT',
                 'Brainstem_E815K_vs_WT', 'Cerebellum_E815K_vs_WT', 'Cortex_E815K_vs_WT', 'Hippocampus_E815K_vs_WT',
                 'Brainstem_D801N_vs_E815K', 'Cerebellum_D801N_vs_E815K', 'Cortex_D801N_vs_E815K', 'Hippocampus_D801N_vs_E815K')
sig_df <- as.data.frame(DE_sig, row.names = DE_genes, col.names = comp_vector)
sig_df[is.na(sig_df)] <- FALSE
sig_df <- 1*sig_df

pdf(file = paste0(base_dir, 'differential_expression/plots/', 'DE_gene_upset_plot.pdf'), width = 8, height = 6, onefile = FALSE)
upset(sig_df, nsets = 12, order.by = "freq")
dev.off()

pdf(file = paste0(base_dir, 'differential_expression/plots/', 'DE_gene_upset_plot_WT_comps.pdf'), width = 8, height = 6, onefile = FALSE)
upset(sig_df[,colnames(sig_df)[grepl("_WT", colnames(sig_df))]], nsets = 9, order.by = "freq")
dev.off()

DE_gene_count <- as.data.frame(colSums(sig_df))
DE_gene_count$comparison <- rownames(DE_gene_count)
colnames(DE_gene_count) <- c('n_DE_genes', 'comparison')
ggplot(DE_gene_count %>%
         mutate(comparison = fct_reorder(comparison, desc(n_DE_genes)))) +
  geom_bar(aes(x=comparison, y = n_DE_genes), stat = 'identity') +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.margin=margin(0, 0, 0, 1, "cm")) +
  ylab('Number of differentially expressed genes')
ggsave(paste0(base_dir, 'differential_expression/plots/', 'DE_gene_histogram.pdf'), width = 8, height = 6)
ggsave(paste0(base_dir, 'differential_expression/plots/', 'DE_gene_histogram.svg'), width = 8, height = 6)

save.image(paste0(base_dir, 'differential_expression/DEseq2.RData'))
