

rstudioapi::getActiveDocumentContext()$path |>
  dirname() |> setwd()
getwd()

library(tidyverse)
library(magrittr)
library(ChIPseeker)
# library(ChIPQC)
library(DiffBind)
library(plyranges)
library(EnrichedHeatmap)
library(profileplyr)
library(eulerr)
library(ggrastr)
# library(cowplot)

## multiprocessing -------------------------
ncore <- max(1, trunc(parallel::detectCores() / 2))

# library(furrr)
# plan(multisession, workers=ncore)
# plan()

# library(BiocParallel)
# if(Sys.getenv()['OS'] == 'Windows_NT') {
#   bpparam <- SnowParam(ncore) 
# } else {
#   bpparam <-  MulticoreParam(ncore)
# }
# print(bpparam)


## read config from yml ----------------------
config <- yaml::read_yaml("bin/config.yml")
bam_dir <- config$bam_dir # "bam"
peak_dir <-  config$peak_dir # "peak"
species <- config$species
cat(bam_dir, peak_dir, species)

# load txdb
if(species == 'hs') {
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  annodb <- 'org.Hs.eg.db'
} else {
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  annodb <- 'org.Mm.eg.db'
}


## edit metadata info ----------------------
c1 <- str_c('Y', 1:28)
c2 <- c(rep('DMSO', 12), rep('ET0825', 8), rep('ET0823', 8))
c3 <- c(rep('YAP', 3), rep('TAZ', 3), rep('panTEAD', 3), rep('H3K27AC', 3),
        rep('YAP', 2), rep('TAZ', 2), rep('panTEAD', 2), rep('H3K27AC', 2),
        rep('YAP', 2), rep('TAZ', 2), rep('panTEAD', 2), rep('H3K27AC', 2))
c4 <- c(rep(19, 6), rep(15, 3), rep(12, 3),
        rep(19, 4), rep(15, 2), rep(12, 2),
        rep(19, 4), rep(15, 2), rep(12, 2))

info <- tibble(Sample = c1, Factor = c3, Treatment = c2, PCR_Cycle = c4) %>%
  mutate(Condition = str_c(Factor, Treatment, sep = '_'),
         Label = str_c(Condition, Sample, sep = '_')) %>%
  mutate(Bam_path = str_c(Sample, '.clean.bam', sep = '') %>% file.path(bam_dir, .),
         Peak_path = str_c(Sample, '_peaks.narrowPeak', sep = '') %>% file.path(peak_dir, .))
# sort samples if wanted
order <- c('YAP_DMSO','YAP_ET0825','YAP_ET0823',
           'TAZ_DMSO','TAZ_ET0825','TAZ_ET0823',
           'panTEAD_DMSO','panTEAD_ET0825','panTEAD_ET0823',
           'H3K27AC_DMSO','H3K27AC_ET0825','H3K27AC_ET0823')
info <- map(order, ~ str_detect(.x, info$Condition) %>% which()) %>%
  purrr::reduce(c) %>%
  info[., ]
print(info)

# construct bam, peak file paths
bam_files <- info$Bam_path %>% set_names(., nm=info$Label)
peak_files <- info$Peak_path %>% set_names(., nm=info$Label)





## ********************************-------------------------------------
## 1.construct dba object  ---------------------
obj <- info %>% 
  filter(Factor == 'YAP') %>%
  dplyr::transmute(SampleID=Label, 
                   Condition=Condition, 
                   Factor=Factor, 
                   Treatment=Treatment, 
                   bamReads=Bam_path, 
                   Peaks=Peak_path, 
                   PeakCaller="narrow", 
                   # Replicate=PCR_Cycle,
                   ControlID=NA, bamControl=NA) %>%
  # filter(Factor %in% c("YAP")) %>%
  dba(sampleSheet=.) %>%
  dba.blacklist(blacklist=ifelse(species=='hs', DBA_BLACKLIST_GRCH38, DBA_BLACKLIST_MM10))

# set config
obj$config$fragmentSize <- 0
obj$config$minQCth <- obj$config$mapQCth <- 10
obj$config$cores <- ncore
obj$config$RunParallel <- ifelse(Sys.getenv()['OS'] == 'Windows_NT', FALSE, TRUE)

## 2.consensus peak within each TF ---------------------

#### method 1, not recommend ####
# consensus <- dba.peakset(obj, consensus=DBA_CONDITION, minOverlap=2)
# consensus <- dba(consensus, mask=consensus$masks$Consensus)
# consensus_peaks <- dba.peakset(consensus, bRetrieve=TRUE)
# mcols(consensus_peaks) <- consensus$called
# colnames(mcols(consensus_peaks)) <- colnames(consensus$class)

#### method 2, recommend #### 
peak_list <- obj$samples$Peaks %>%
  set_names(obj$samples$SampleID) %>%
  map(ChIPseeker::readPeakFile) %>%
  # filter out low quality peak 
  # foldchange > 2, -log10(qvalue) > 2
  map(~ plyranges::filter(.x, V7 > 2, V9 > 2)) %>%
  GRangesList()
# peak_list$YAP_DMSO_Y1$V7 %>% summary()
# peak_list$YAP_DMSO_Y1$V9 %>% summary()

consensus_peaks <- peak_list %>%
  endoapply(mutate, center=start + V10) %>% 
  extraChIPs::makeConsensus(p=0, var=c('V5','V7','center')) %>%
  mutate(score=vapply(V5, mean, numeric(1)),
         fold=vapply(V7, mean, numeric(1)),
         center=vapply(center, mean, numeric(1))) %>%
  plyranges::select(-V5, -V7) %>%
  plyranges::filter(n > 1) # overlap = 2
# relocate
mcols(consensus_peaks) <- mcols(consensus_peaks) %>%
  as.data.frame() %>%
  dplyr::relocate(center, .after=last_col())

print('n:')
print(consensus_peaks$n %>% table())
print('fold:')
print(consensus_peaks$fold %>% summary())
print('score:')
print(consensus_peaks$score %>% summary())
print('width:')
print(consensus_peaks %>% width() %>% summary())

# # filter consensus peak by width
# width(consensus_peaks) %>% density() %T>% plot()
# consensus_peaks <- consensus_peaks %>% 
#   plyranges::filter(., width(.) < 3000, width(.) > 50)

## keep standard chromsome
seqlevels(consensus_peaks)
consensus_peaks <- consensus_peaks %>%
  keepStandardChromosomes(pruning.mode='coarse') %>%
  dropSeqlevels('MT', pruning.mode='coarse') %>%
  dropSeqlevels('chrMT', pruning.mode='coarse')


## 3.count binding matrix  ---------------------
# calculate consensus peaks with fixed width
obj <- dba.count(obj, summits=F,
                 score=DBA_SCORE_READS,
                 peaks=consensus_peaks, 
                 bRemoveDuplicates=T)

obj <- dba.normalize(obj, method=DBA_DESEQ2, 
                     normalize=DBA_NORM_LIB, library=DBA_LIBSIZE_FULL)

binding_gr <- dba.peakset(obj, bRetrieve=T)

fc <- rowSums(as.data.frame(mcols(binding_gr)[,4:5])) / 
  rowSums(as.data.frame(mcols(binding_mat_raw[,1:3])))
row_mean <- rowMeans(as.data.frame(mcols(binding_gr)[,1:5]))

# normalize
norm_factor <- dba.normalize(obj, bRetrieve=T)[['norm.factors']]

binding_norm <- mcols(binding_gr) %>% as.data.frame()
for(i in seq_along(norm_factor)) {
  binding_norm[, i] <- binding_norm[, i] / norm_factor[i]
}
binding_norm[1:5,]

fc_norm <- rowSums(binding_norm[,4:5]) / 
  rowSums(binding_norm[,1:3])
row_mean_norm <- rowMeans(binding_norm[,1:5])


# merge
mcols(binding_gr)$fc <- fc
mcols(binding_gr)$row_mean <- row_mean
mcols(binding_gr)$fc_norm <- fc_norm
mcols(binding_gr)$row_mean_norm <- row_mean_norm


# annotate
seqlevelsStyle(binding_gr) <- 'UCSC'
seqlevels(binding_gr)

peak_et825 <- binding_gr %>% 
  annotatePeak(level='gene', 
               TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
               annoDb='org.Hs.eg.db')
peak_et825 <- peak_et825@anno %>%
  as.data.frame()

# save
as.data.frame(peak_et825) %>%
  writexl::write_xlsx(str_glue('{Sys.Date()}_yap_foldchange.xlsx'))





# RNAseq on MSTO, 1uM ET0825 treated at 24h
rna_et825 <- read_xlsx(path('rnaseq_deseq2', '2022-11-23_wald_test_ET0825.0.3uM_Veh.Ctr.xlsx'), sheet=1)
# rna_et823 <- read_xlsx(path('rnaseq_deseq2', '2022-11-23_wald_test_ET0823.1uM_Veh.Ctr.xlsx'), sheet=1)
# 
# # read yap peaks
# peak_et825 <- '../db_analysis/2022-11-22_db_result_YAP.xlsx' %>% read_xlsx(sheet='ET0825_vs_DMSO') %>%
#   filter(str_detect(annotation, 'Promoter'))
# peak_et823 <- '../db_analysis/2022-11-22_db_result_YAP.xlsx' %>% read_xlsx(sheet='ET0823_vs_DMSO') %>%
#   filter(str_detect(annotation, 'Promoter'))

# common genes
all_genes <- intersect(pull(rna_et825, gene_name), pull(peak_et825, SYMBOL))

# filter by common genes
rna_et825 <- filter(rna_et825, gene_name %in% all_genes) %>%
  dplyr::select(gene=gene_name, lfc=log2FoldChange, pvalue)

peak_et825 <- filter(peak_et825, SYMBOL %in% all_genes) %>%
  group_by(SYMBOL) %>%
  summarise(lfc=mean(log2(fc))) %>%
  dplyr::rename(gene=SYMBOL)

# merge
et825 <- left_join(peak_et825, rna_et825, by='gene', suffix=c('_peak', '_rna'))


# plot
et825 %>%
  # filter(lfc_peak < 0, lfc_rna > 0) %>%
  # arrange(lfc_peak)
  # filter(gene %in% c('ANKRD1'))
  mutate(label=ifelse(gene %in% c('CCN1','CCN2','ANKRD1'), gene, NA)) %>%
  # filter(pvalue_peak < 0.01, pvalue_rna < 0.01) %>%
  ggplot(aes(lfc_rna, lfc_peak, label=label)) +
  geom_point() +
  geom_label() +
  stat_cor(method='kendall')


















