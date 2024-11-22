#module load R/4.3.3

#please run DESeq.R first to get the Rdata image to load

library('tidyverse')
library('UpSetR')

#to load the workspace
#change this to where the RData object exists in your path after running DESeq.R
out_dir <- '/gpfs/commons/groups/knowles_lab/sadamson/JAX_Atp1a3_project/github_repo/atp1a3_mouse_rna_seq_analysis/differential_expression/'
load(paste0(out_dir, 'DEseq2.RData'))

#make upset plots of DE genes for each comparison
all_comps <- list(BS_D801N_WT_res, CM_D801N_WT_res, CX_D801N_WT_res, HP_D801N_WT_res,
               BS_E815K_WT_res, CM_E815K_WT_res, CX_E815K_WT_res, HP_E815K_WT_res)
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
                 'Brainstem_E815K_vs_WT', 'Cerebellum_E815K_vs_WT', 'Cortex_E815K_vs_WT', 'Hippocampus_E815K_vs_WT')
sig_df <- as.data.frame(DE_sig, row.names = DE_genes, col.names = comp_vector)
sig_df[is.na(sig_df)] <- FALSE
sig_df <- 1*sig_df

pdf(file = paste0('plots/', 'Figure_4A_DE_gene_upset_plots.pdf'), width = 8, height = 6, onefile = FALSE)
upset(sig_df, nsets = 9, order.by = "freq")
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
ggsave(paste0('plots/QC_plots/', 'DE_gene_histogram.pdf'), width = 8, height = 6)
ggsave(paste0('plots/QC_plots/', 'DE_gene_histogram.svg'), width = 8, height = 6)

#make volcano plots looking at neuroinflammation genes
neuroinflammation_gene_tib <- read_csv('Supplementary_file_7_neuroinflammation_genes.csv')
#gene_name,gene_type
#0610009B22RIK,GOAD_microglial
neuroinflammation_genes <- neuroinflammation_gene_tib$gene_name
CM_D801N_WT_res <- CM_D801N_WT_res %>% 
    mutate(gene_name = toupper(gene_name)) %>%
    mutate(Flam_gene = gene_name %in% neuroinflammation_genes)

y_bound <- 10
p <- ggplot()
p <- p + geom_point(data = CM_D801N_WT_res %>%  filter(padj <= 0.05, Flam_gene == FALSE),
                    aes(x = log2FoldChange, y = pmin(-log10(padj), y_bound)), pch = 19, color = 'aquamarine3',size = 3, alpha= 0.5, stroke = 0)
p <- p + geom_point(data = CM_D801N_WT_res %>%  filter(padj <= 0.05, Flam_gene  == TRUE),
                    aes(x = log2FoldChange, y = pmin(-log10(padj), y_bound)), colour = 'red1')

p <- p + theme(panel.background = element_blank(), legend.title = element_blank(), legend.position = "none")
p <- p + geom_vline(xintercept = 1, colour = 'grey50', linetype = 'longdash') + 
    geom_vline(xintercept = -1, colour = 'grey50', linetype = 'longdash')
p <- p + xlim(-2.5,7) #+ ylim(0, y_bound)
p
ggsave(paste0('plots/publication_plots/', "Supplemental_Figure_S6A_neuroinflammation_volcano_CM_D801N_WT.pdf"), width = 6.3, height = 6.3) 
ggsave(paste0('plots/publication_plots/', "Supplemental_Figure_S6A_neuroinflammation_volcano_CM_D801N_WT.svg"), width = 6.3, height = 6.3) 

CM_E815K_WT_res <- CM_E815K_WT_res %>% 
    mutate(gene_name = toupper(gene_name)) %>%
    mutate(Flam_gene = gene_name %in% neuroinflammation_genes)

y_bound <- 10
p <- ggplot()
p <- p + geom_point(data = CM_E815K_WT_res %>%  filter(padj <= 0.05, Flam_gene == FALSE),
                    aes(x = log2FoldChange, y = pmin(-log10(padj), y_bound)), pch = 19, color = 'orange1',size = 3, alpha= 0.5, stroke = 0)
p <- p + geom_point(data = CM_E815K_WT_res %>%  filter(padj <= 0.05, Flam_gene  == TRUE),
                    aes(x = log2FoldChange, y = pmin(-log10(padj), y_bound)), colour = 'red1')

p <- p + theme(panel.background = element_blank(), legend.title = element_blank(), legend.position = "none")
p <- p + geom_vline(xintercept = 1, colour = 'grey50', linetype = 'longdash') + 
    geom_vline(xintercept = -1, colour = 'grey50', linetype = 'longdash')
p <- p  + xlim(-2.5,7) # + ylim(0, y_bound)
p
ggsave(paste0('plots/publication_plots/', "Supplemental_Figure_S6B_neuroinflammation_volcano_CM_E815K_WT.pdf"), width = 6.3, height = 6.3) 
ggsave(paste0('plots/publication_plots/', "Supplemental_Figure_S6B_neuroinflammation_volcano_CM_E815K_WT.svg"), width = 6.3, height = 6.3) 

#boxplots of neuroinflammation genes
p <- ggplot(data = CM_D801N_WT_res %>% filter(padj <=0.05), aes(x = Flam_gene, y = log2FoldChange)) 
p <- p + geom_boxplot()
p <- p + ylim(-3.5, 6.5) + theme_bw()
p
ggsave(paste0('plots/publication_plots/', "Supplemental_Figure_S6C_neuroinflammation_volcano_CM_D801N_WT.pdf"), width = 6.3, height = 6.3) 
ggsave(paste0('plots/publication_plots/', "Supplemental_Figure_S6C_neuroinflammation_volcano_CM_D801N_WT.svg"), width = 6.3, height = 6.3) 

wilcox.test(log2FoldChange ~ Flam_gene, data = CM_D801N_WT_res, alternative = "less") # p-value < 2.2e-16

p <- ggplot(data = CM_E815K_WT_res %>% filter(padj <=0.05), aes(x = Flam_gene, y = log2FoldChange)) 
p <- p + geom_boxplot()
p <- p + ylim(-3.5, 6.5) +theme_bw()
p
ggsave(paste0('plots/publication_plots/', "Supplemental_Figure_S6C_neuroinflammation_volcano_CM_E815K_WT.pdf"), width = 6.3, height = 6.3) 
ggsave(paste0('plots/publication_plots/', "Supplemental_Figure_S6C_neuroinflammation_volcano_CM_E815K_WT.svg"), width = 6.3, height = 6.3) 

wilcox.test(log2FoldChange ~ Flam_gene, data = CM_E815K_WT_res, alternative = "less") #p-value < 2.2e-16

CM_E815K_neuroinflammation_lfc <- pull(CM_E815K_WT_res %>% filter(padj <=0.05, Flam_gene == TRUE) %>% dplyr::select(log2FoldChange))
CM_D801N_neuroinflammation_lfc <- pull(CM_D801N_WT_res %>% filter(padj <=0.05, Flam_gene == TRUE) %>% dplyr::select(log2FoldChange))
wilcox.test(CM_D801N_neuroinflammation_lfc, CM_E815K_neuroinflammation_lfc, alternative = "less") #p-value = 1.249e-06

#brainstem boxplots
BS_D801N_WT_res <- BS_D801N_WT_res %>% 
    mutate(gene_name = toupper(gene_name)) %>%
    mutate(Flam_gene = gene_name %in% neuroinflammation_genes)

p <- ggplot(data = BS_D801N_WT_res %>% filter(padj <=0.05), aes(x = Flam_gene, y = log2FoldChange)) 
p <- p + geom_boxplot()
p <- p + ylim(-3.5, 6.5) + theme_bw()
p
ggsave(paste0('plots/publication_plots/', "Supplemental_Figure_S6D_neuroinflammation_volcano_BS_D801N_WT.pdf"), width = 6.3, height = 6.3) 
ggsave(paste0('plots/publication_plots/', "Supplemental_Figure_S6D_neuroinflammation_volcano_BS_D801N_WT.svg"), width = 6.3, height = 6.3) 

wilcox.test(log2FoldChange ~ Flam_gene, data = BS_D801N_WT_res %>% filter(padj <=0.05), alternative = "less") #p-value = 0.5249

BS_E815K_WT_res <- BS_E815K_WT_res %>% 
    mutate(gene_name = toupper(gene_name)) %>%
    mutate(Flam_gene = gene_name %in% neuroinflammation_genes)

p <- ggplot(data = BS_E815K_WT_res %>% filter(padj <= 0.05), aes(x = Flam_gene, y = log2FoldChange)) 
p <- p + geom_boxplot()
p <- p + ylim(-3.5, 6.5) + theme_bw()
p
ggsave(paste0('plots/publication_plots/', "Supplemental_Figure_S6D_neuroinflammation_volcano_BS_E815K_WT.pdf"), width = 6.3, height = 6.3) 
ggsave(paste0('plots/publication_plots/', "Supplemental_Figure_S6D_neuroinflammation_volcano_BS_E815K_WT.svg"), width = 6.3, height = 6.3) 

wilcox.test(log2FoldChange ~ Flam_gene, data = BS_E815K_WT_res %>% filter(padj <=0.05), alternative = "less") # p-value = 4.872e-07

BS_E815K_neuroinflammation_lfc <- pull(BS_E815K_WT_res %>% filter(padj <= 0.05, Flam_gene == TRUE) %>% dplyr::select(log2FoldChange))
BS_D801N_neuroinflammation_lfc <- pull(BS_D801N_WT_res %>% filter(padj <= 0.05, Flam_gene == TRUE) %>% dplyr::select(log2FoldChange))
wilcox.test(BS_D801N_neuroinflammation_lfc, BS_E815K_neuroinflammation_lfc, alternative = "less") # p-value = 6.891e-07


