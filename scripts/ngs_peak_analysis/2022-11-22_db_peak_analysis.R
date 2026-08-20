

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








dba_analysis <- function(info, factor='YAP', color=c("white","#e41a1c")) {
  
  ## ********************************-------------------------------------
  ## 1.construct dba object  ---------------------
  obj <- info %>% 
    filter(Factor == factor) %>%
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
    makeConsensus(p=0, var=c('V5','V7','center')) %>%
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
  
  ## save consensus peak
  writexl::write_xlsx(as.data.frame(consensus_peaks), 
                      str_glue('{Sys.Date()}_consensus_peaks_{factor}.xlsx'))
  # rtracklayer::export.bed(consensus_peaks, 
  #                         str_glue('{Sys.Date()}_consensus_peaks_{factor}.bed'))
  
  
  ## 3.consensus peak overlap analysis  ---------------------
  pdf(str_glue("{Sys.Date()}_consensus_peaks_overlap_{factor}.pdf"), width=7, height=7)
  dat <- mcols(consensus_peaks) %>%
    as.data.frame() %>%
    dplyr::select(-n:-center)
  
  dat0 <- data.frame(DMSO=rowSums(dat[, 1:3]) > 0, 
                    ET0825=rowSums(dat[, 4:5]) > 0, 
                    ET0823=rowSums(dat[, 6:7]) > 0)
  print(plot(euler(dat0),
       legend=T, quantities=T, 
       type=c('counts','percent'),
       fills=list(fill=RColorBrewer::brewer.pal(n=ncol(dat0), 'Dark2'), alpha=0.8),
       labels=list(col='black', font=2),
       main='Overlap of at least 1 samples per condtion' ) )
  
  dat1 <- data.frame(DMSO=rowSums(dat[, 1:3]) > 1, 
                     ET0825=rowSums(dat[, 4:5]) > 1, 
                     ET0823=rowSums(dat[, 6:7]) > 1)
  print(plot(euler(dat1),
       legend=T, quantities=T, 
       type=c('counts','percent'),
       fills=list(fill=RColorBrewer::brewer.pal(n=ncol(dat1), 'Dark2'), alpha=0.8),
       labels=list(col='black', font=2),
       main='Overlap of at least 2 samples per condtion'))
  while(dev.cur()!=1) dev.off()
  
  
  ## 4.count binding matrix  ---------------------
  # calculate consensus peaks with fixed width
  obj <- dba.count(obj, summits=F,
                   peaks=consensus_peaks, 
                   bRemoveDuplicates=T)
  # normalize
  obj <- dba.normalize(obj, method=DBA_DESEQ2, 
                       normalize=DBA_NORM_LIB, library=DBA_LIBSIZE_FULL)
  dba.normalize(obj, bRetrieve=T)
  
  dba.show(obj) %>% as_tibble() %>%
    mutate(peakReads = round(Reads * FRiP), outReads=Reads - peakReads) %>%
    pivot_longer(c(outReads, peakReads), names_to="type", values_to="value") %>%
    ggplot(aes(ID, value, fill=type)) +
    geom_bar(stat='identity', position='fill', width=.7) +
    scale_fill_manual(values=c('outReads'='grey80', 'peakReads'='steelblue')) +
    labs(x=NULL, y="Read Counts") +
    ggpubr::theme_pubr(x.text.angle = 45)
  cowplot::ggsave2(str_glue('{Sys.Date()}_consensus_peaks_frip_{factor}.pdf'), height=4, width=4)

  
  ## 5.dba qc   ---------------------
  pdf(str_glue('{Sys.Date()}_db_qc_{factor}.pdf'))
  dba.plotPCA(obj, attributes=DBA_CONDITION, score=DBA_SCORE_NORMALIZED, label=DBA_CONDITION)
  # heatmap correlation
  dba.plotHeatmap(obj, DBA_ID)
  # heatmap peak counts
  # dba.plotHeatmap(obj, DBA_ID, maxSites=obj$totalMerged, correlations=F, scale="row")
  while(dev.cur()!=1) dev.off()
  
  ## 6.consensus peak profile  ---------------------
  # all samples share same color scale by merge samples within condition
  consensus_profiles <- dba.plotProfile(obj, 
                                      samples=obj$masks$All, 
                                      sites=consensus_peaks,
                                      normalize=TRUE, 
                                      maxSites=obj$totalMerged,
                                      labels=DBA_TREATMENT, 
                                      merge=DBA_ID,
                                      style='point', distanceAround=3000, nOfWindows=50)
  # modify parameter
  sampleData(consensus_profiles)$Condition <- rownames(sampleData(consensus_profiles))
  sampleData(consensus_profiles)$Factor <- sampleData(consensus_profiles)$Condition %>%
    map_chr(~ str_split(.x, '_', simplify=T)[1])
  sampleData(consensus_profiles)$Treatment <- sampleData(consensus_profiles)$Condition %>%
    map_chr(~ str_split(.x, '_', simplify=T)[2])
  consensus_profiles@rowRanges@elementMetadata$Peak <- 'Consensus'
  consensus_profiles@params$rowGroupsInUse <- 'Peak'
  # plot 
  pdf(str_glue('{Sys.Date()}_consensus_peaks_profile_{factor}.pdf'), width=5, height=7)
  generateEnrichedHeatmap(consensus_profiles,
                          matrices_axis_name=c('-3k', '', '3k'),
                          all_color_scales_equal = T,
                          # color_by_sample_group = "Treatment",
                          matrices_color = color,
                          use_raster=T, raster_quality=3)
  while(dev.cur()!=1) dev.off()
  


  
  
  
  ## ********************************-------------------------------------
  
  ## 1.db analysis setup    ---------------------
  # obj$contrasts <- NULL
  obj <- obj %>%
    dba.contrast(contrast=c('Treatment', 'ET0825', 'DMSO')) %>%
    dba.contrast(contrast=c('Treatment', 'ET0823', 'DMSO')) %>%
    dba.analyze()
  
  # get result
  compares <- dba.show(obj, bContrasts=T, bDesign=T)

  compares_results <- map(seq_len(nrow(compares)), 
                          ~ dba.report(obj, contrast=.x, th=1, fold=0, 
                                       bCounts=T, bNormalized=T, 
                                       bDB=T, bNotDB=T,
                                       bAll=T, bGain=F, bLoss=F))
  compares_results <- map(compares_results, ~ .x %>% .[['peaks']] %>% .[[1]] %>% GRanges)
  names(compares_results) <- str_c(compares$Group, '_vs_', compares$Group2)
  
  # 2.db peaks annotate ---------------------
  compares_results <- map(
    compares_results,
    function(res) {
      seqlevelsStyle(res) <- 'UCSC'
      res <- ChIPseeker::annotatePeak(res, TxDb=txdb, annoDb=annodb, level='gene')
      # rowRanges(profile) %>% mcols() %>% .[['annotation_short']] %>% unique()
      return(res) })
  # get Ganges
  compares_results <- map(compares_results, ~ .x@anno)
  
  # 3.db results save ---------------------
  map(compares_results, ~ as.data.frame(.x)) %>%
    writexl::write_xlsx(str_glue('{Sys.Date()}_db_result_{factor}.xlsx'))
  
  
  # 4.set threshold parameter   ---------------------
  th <- 0.05 # FDR threshold
  fold <- round(log2(1.5), 2) # FoldChange=1.5, fold=log2(FoldChange)
  
  # 5.db volcano plot  ---------------------
  pdf(str_glue('{Sys.Date()}_db_volcano_{factor}.pdf'), width=4, height=5)
  
  res <- compares_results$ET0825_vs_DMSO
  volcano_data <- as.data.frame(res) %>%
    mutate(Direction=ifelse(Fold > fold & FDR < th, 'Gain', 
                            ifelse(Fold < -fold & FDR < th, 'Loss', 'NS')))
  N_Gain <- filter(volcano_data, Direction=="Gain") %>% nrow()
  N_Loss <- filter(volcano_data, Direction=="Loss") %>% nrow()
  x_max <- range(volcano_data$Fold) %>% abs %>% max()
  p <- volcano_data %>%
    ggplot(aes(Fold, -log10(FDR), color=Direction)) +
    geom_point(show.legend=F, alpha=0.5) +
    scale_color_manual(values=c('Gain'='#ca0020', 'Loss'='#0571b0', 'NS'='grey50')) +
    coord_cartesian(xlim=c(-x_max, x_max)) +
    labs(x='log2(FoldChange)', y='-log10(FDR)', 
         title=names(compares_results)[1], 
         subtitle=str_glue('FDR < {th} & abs(log2(FC)) > {fold}, N_Gain={N_Gain}, N_Loss={N_Loss}')) +
    theme_bw(12) # +  theme(panel.grid=element_blank())
  ggrastr::rasterise(p, dpi=600) %>% print()
  
  res <- compares_results$ET0823_vs_DMSO
  volcano_data <- as.data.frame(res) %>%
    mutate(Direction=ifelse(Fold > fold & FDR < th, 'Gain', 
                            ifelse(Fold < -fold & FDR < th, 'Loss', 'NS')))
  N_Gain <- filter(volcano_data, Direction=="Gain") %>% nrow()
  N_Loss <- filter(volcano_data, Direction=="Loss") %>% nrow()
  x_max <- range(volcano_data$Fold) %>% abs %>% max()
  p <- volcano_data %>%
    ggplot(aes(Fold, -log10(FDR), color=Direction)) +
    geom_point(show.legend=F, alpha=0.5) +
    scale_color_manual(values=c('Gain'='#ca0020', 'Loss'='#0571b0', 'NS'='grey50')) +
    coord_cartesian(xlim=c(-x_max, x_max)) +
    labs(x='log2(FoldChange)', y='-log10(FDR)', 
         title=names(compares_results)[2], 
         subtitle=str_glue('FDR < {th} & abs(log2(FC)) > {fold}, N_Gain={N_Gain}, N_Loss={N_Loss}')) +
    theme_bw(12) # +  theme(panel.grid=element_blank())
  ggrastr::rasterise(p, dpi=600) %>% print()
  
  while(dev.cur()!=1) dev.off()
  
  
  # 6.db significant peak overlap ---------------------
  compares_results_sig <- compares_results %>%
    map(~ plyranges::filter(.x, FDR < th, abs(Fold) > fold)) # %>% lengths() %>% sum()
  
  pdf(str_glue('{Sys.Date()}_db_overlap_{factor}.pdf'), width=8, height=7)
  set.seed(123)
  list(
    ET0825_Gain=as.data.frame(compares_results_sig$ET0825_vs_DMSO) %>% 
      plyranges::filter(Fold > 0) %>% 
      unite('id', seqnames, start, end) %>% pull(id),
    ET0825_Loss=as.data.frame(compares_results_sig$ET0825_vs_DMSO) %>% 
      plyranges::filter(Fold < 0) %>% 
      unite('id', seqnames, start, end) %>% pull(id),
    ET0823_Gain=as.data.frame(compares_results_sig$ET0823_vs_DMSO) %>% 
      plyranges::filter(Fold > 0) %>% 
      unite('id', seqnames, start, end) %>% pull(id),
    ET0823_Loss=as.data.frame(compares_results_sig$ET0823_vs_DMSO) %>% 
      plyranges::filter(Fold < 0) %>% 
      unite('id', seqnames, start, end) %>% pull(id)
  ) %>%
    euler(shape='ellipse') %>% 
    plot(legend=T, quantities=T,
         labels=T, #labels=list(col='black', font=2,
         fills=list(fill=RColorBrewer::brewer.pal(4, 'Paired'), alpha=0.8)) %>%
    print()
  while(dev.cur()!=1) dev.off()
  
  
  
  # 7.db significant peak profile ---------------------
  # peak profile
  compares_results_sig <- map(
    compares_results_sig,
    function(res) { 
      seqlevelsStyle(res) <- 'NCBI'
      return(res) })

  compares_results_profiles <- compares_results_sig %>%
    map(~ dba.plotProfile(obj, 
                          samples=obj$masks$Consensus,
                          sites=.x, 
                          maxSites=max(lengths(compares_results_sig)),
                          labels=DBA_TREATMENT, 
                          merge=DBA_ID,
                          style='point', distanceAround=3000, nOfWindows=50))
  # add attributes
  compares_results_profiles <- map(
    compares_results_profiles,
    function(profile) {
      # col annotation
      sampleData(profile)$Condition <- rownames(sampleData(profile))
      sampleData(profile)$Factor <- sampleData(profile)$Condition %>%
        map_chr(~ str_split(.x, '_', simplify=T)[1])
      sampleData(profile)$Treatment <- sampleData(profile)$Condition %>%
        map_chr(~ str_split(.x, '_', simplify=T)[2])
      # row annotation
      rowRanges(profile)$`Binding` <- rowRanges(profile)$Fold %>%
        map_chr(~ ifelse(.x > 0, 'Gain', 'Loss'))
      # rowRanges(profile)$`annotation_short` <- rowRanges(profile)$annotation %>%
      #   map_chr(~ str_split(.x, ' ', simplify=T)[1]) %>% 
      #   map_chr(~ ifelse(.x %in% c('Exon','Intron',"3'","5'"), 'Gene Body',
      #                    ifelse(.x %in% c('Distal','Downstream'), 'Distal', 'Promoter'))) # %>% unique()
      # profile@params$rowGroupsInUse <- 'Binding'
      return(profile) })
  
  # 8.db significant peak plot  ---------------------
  pdf(str_glue('{Sys.Date()}_db_profile_{factor}.pdf'), width=5, height=7)
  # by db gain and loss
  profile <- compares_results_profiles[[1]]
  profile@params$rowGroupsInUse <- 'Binding'
  if(length(profile)>0) {
    generateEnrichedHeatmap(profile,
                            matrices_axis_name=c('-3k', '', '3k'),
                            decreasing=T,
                            # color_by_sample_group = "Treatment",
                            matrices_color = color,
                            genes_to_label=c('CCN1','CCN2','ANKRD1'),
                            use_raster=T, raster_quality=3)
  }

  # # by annotation_short
  # profile@params$rowGroupsInUse <- 'annotation_short'
  # generateEnrichedHeatmap(profile,
  #                         matrices_axis_name=c('-3k', '', '3k'),
  #                         decreasing=T,
  #                         # color_by_sample_group = "Treatment",
  #                         matrices_color = color,
  #                         genes_to_label=c('CCN1','CCN2','ANKRD1'),
  #                         use_raster=T, raster_quality=3)
  
  # by db gain and loss
  profile <- compares_results_profiles[[2]]
  if(length(profile)>0) {
    profile@params$rowGroupsInUse <- 'Binding'
    generateEnrichedHeatmap(profile,
                            matrices_axis_name=c('-3k', '', '3k'),
                            decreasing=T,
                            # color_by_sample_group = "Treatment",
                            matrices_color = color,
                            genes_to_label=c('CCN1','CCN2','ANKRD1'),
                            use_raster=T, raster_quality=3)
  }

  while(dev.cur()!=1) dev.off()
  
}



dba_analysis(info, 'YAP', c("white","#e41a1c"))
dba_analysis(info, 'TAZ', c("white","#377eb8"))
dba_analysis(info, 'panTEAD', c("white","#4daf4a"))
dba_analysis(info, 'H3K27AC', c("white","#984ea3"))





