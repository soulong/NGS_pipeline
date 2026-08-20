
rstudioapi::getActiveDocumentContext()$path |>
  dirname() |> setwd()
getwd()

library(tidyverse)
library(magrittr)
# library(ChIPseeker)
# library(ChIPQC)
# library(DiffBind)
# library(plyranges)
# library(EnrichedHeatmap)
# library(profileplyr)
# library(eulerr)
library(ggrastr)
library(readxl)
library(fs)
library(ggpubr)



# RNAseq on MSTO, 1uM ET0825 treated at 24h
rna_et825 <- read_xlsx(path('rnaseq_deseq2', '2022-11-23_wald_test_ET0825.0.3uM_Veh.Ctr.xlsx'), sheet=1)
rna_et823 <- read_xlsx(path('rnaseq_deseq2', '2022-11-23_wald_test_ET0823.1uM_Veh.Ctr.xlsx'), sheet=1)

# read yap peaks
peak_et825 <- 'db_analysis/2022-11-22_db_result_YAP.xlsx' %>% read_xlsx(sheet='ET0825_vs_DMSO') %>%
  filter(str_detect(annotation, 'Promoter'))
peak_et823 <- 'db_analysis/2022-11-22_db_result_YAP.xlsx' %>% read_xlsx(sheet='ET0823_vs_DMSO') %>%
  filter(str_detect(annotation, 'Promoter'))


# common genes
all_genes <- intersect(pull(rna_et825, gene_name), pull(peak_et825, SYMBOL))

# filter by common genes
rna_et825 <- filter(rna_et825, gene_name %in% all_genes) %>%
  dplyr::select(gene=gene_name, lfc=log2FoldChange, pvalue)
rna_et823 <- filter(rna_et823, gene_name %in% all_genes) %>%
  dplyr::select(gene=gene_name, lfc=log2FoldChange, pvalue)

peak_et825 <- filter(peak_et825, SYMBOL %in% all_genes) %>%
  group_by(SYMBOL) %>%
  summarise(lfc=mean(Fold), pvalue=mean(p.value)) %>%
  dplyr::rename(gene=SYMBOL)
peak_et823 <- filter(peak_et823, SYMBOL %in% all_genes) %>%
  group_by(SYMBOL) %>%
  summarise(lfc=mean(Fold), pvalue=mean(p.value)) %>%
  dplyr::rename(gene=SYMBOL)

# merge
et825 <- left_join(peak_et825, rna_et825, by='gene', suffix=c('_peak', '_rna'))
et823 <- left_join(peak_et823, rna_et823, by='gene', suffix=c('_peak', '_rna'))


yap_targeted_genes <- c("CCN1","ANKRD1","CCN2","TOP2A","KIF14","CCNA2","CDCA8","CENPF","KIF23",
                        "KIF20B","KNTC1","RRM2","MCM3","TUBB","MYBL1","RAD18",
                        "ZWILCH","SGO1","TIMELESS","GINS1","SMC3","TK1","MRE11",
                        "MCM7","SUV39H2","GADD45B","FOSL1","CENPV","RUVBL2","MYC","GLI2","AXL",
                        "ABCB1","CAT","GPATCH4","LMNB2","TXN","WSB2","AREG","FOXF2","IGFBP3",
                        "RASSF2","AMOTL2","NPPB","CCND1")


# plot
et825 %>%
  filter(gene %in% yap_targeted_genes) %>%
  # filter(lfc_peak < 0, lfc_rna > 0) %>%
  # arrange(lfc_peak)
  # filter(gene %in% c('ANKRD1'))
  mutate(label=ifelse(gene %in% c('CCN1','CCN2','ANKRD1'), gene, NA)) %>%
  # filter(pvalue_peak < 0.01, pvalue_rna < 0.01) %>%
  ggplot(aes(lfc_rna, lfc_peak, label=label)) +
  geom_point(size=3) +
  # geom_label() +
  stat_cor() +
  geom_smooth(method='lm', size=2, color='#e41a1c', fill='#e41a1c', alpha=0.1) +
  theme_bw(16) +
  labs(x='log2(FoldChange) - RNA level', 
       y='log2(FoldChange) - Peak level',
       title='Correlation of typical YAP targeted gene subset')
cowplot::ggsave2('rnaseq_cross_correlation_yap_tageted_genes.pdf')


