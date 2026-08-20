
rstudioapi::getActiveDocumentContext()$path |>
  dirname() |> setwd()
getwd()

library(tidyverse)
library(magrittr)
# library(ChIPseeker)
# library(ChIPQC)
# library(DiffBind)
library(plyranges)
library(EnrichedHeatmap)
library(profileplyr)
# library(eulerr)
library(ggrastr)
library(readxl)
library(fs)
library(ggpubr)
library(BiocParallel)


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

ncore <- max(8, trunc(parallel::detectCores() / 8))


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

# # construct bam, peak file paths
# bam_files <- info$Bam_path
# peak_files <- info$Peak_path



## db peaks
get_profile <- function(factor='YAP') {
  db_peaks <- readxl::read_xlsx(str_glue('db_analysis/2022-11-22_db_result_{factor}.xlsx')) %>%
    as_granges()
  # seqlevelsStyle(db_peaks_yap) <- 'NCBI'
  bed_file <- file_temp(ext='bed')
  rtracklayer::export.bed(db_peaks, bed_file)
  
  info_subset <- filter(info, Factor == factor)
  bam_files <- info_subset$Bam_path
  chipProfile <- BamBigwig_to_chipProfile(bam_files, # bw_files 
                                          bed_file,
                                          format = "bam",
                                          style = "point", 
                                          nOfWindows = 100, bin_size = 50, distanceAround = 3000,
                                          removeDup = TRUE, paired = TRUE, normalize = "RPM",
                                          quant_params = SnowParam(ncore) # for win
                                          # quant_params = MulticoreParam(round(parallel::detectCores()/2)) 
  )
  # naming assays
  names(assays(chipProfile)) <- info_subset$Label
  
  # convert to profileplyr
  profile <- as_profileplyr(chipProfile)
  
  # convert style to UCSC
  seqlevelsStyle(rowRanges(profile)) <- "UCSC" # rowRanges(profile)
  rowRanges(profile)
  
  # rename sample
  rownames(sampleData(profile)) <- rownames(sampleData(profile)) %>%
    map_chr(~ str_split(.x, fixed(".clean.bam"), simplify = T)[1]) %>%
    match(info_subset$Sample) %>%
    info_subset$Label[.]
  
  # annotate regions
  profile <- annotateRanges(profile, 
                            TxDb = "hg38",
                            # annoDb = org.Hs.eg.db,
                            tssRegion = c(-3000, 3000))
  # merge exon and intron
  profile@rowRanges$Type <- 
    profile@rowRanges$annotation_short %>%
    str_replace_all("(Exon)|(Intron)|(5p UTR)|(3p UTR)|(Downstream)", "Gene Body") %>%
    factor(levels = c("Promoter", "Gene Body", "Distal Intergenic"))
  table(profile@rowRanges$Type)
  # group row by annotation
  profile@params$rowGroupsInUse <- 'Type'
  
  # group sample by target
  sampleData(profile)$Factor <- rownames(sampleData(profile)) %>%
    match(info_subset$Label) %>% info_subset$Factor[.]
  sampleData(profile)$Treatment <- rownames(sampleData(profile)) %>%
    match(info_subset$Label) %>% info_subset$Treatment[.]
  
  # raw counts
  if(T) {
    profile_raw <- profile

    # bam reads
    bam_reads <- chipProfile@metadata[["AlignedReadsInBam"]] %>%
      set_names(nm=info_subset$Label)
    norm_factor <- bam_reads/1e6
    
    # back to raw reads
    for(i in names(bam_reads)) {
      assays(profile_raw)[[i]] <- assays(chipProfile)[[i]] * norm_factor[i]
    }
  }
  
  return(list(rpm=profile, raw=profile_raw))
}



merge_profile <- function(profile) {
  # merge sample
  sample_data <-  sampleData(profile) %>%
    as_tibble(rownames='id') %>%
    dplyr::select(id, Factor, Treatment) %>%
    unite('Condition', Factor, Treatment, remove=F)
  sample_0 <- map_chr(split(sample_data, sample_data$Condition), 
                      ~ pull(.x, id)[1])
  # sort according to raw order
  sample_0 <- names(sample_0) %>%
    match(sample_data$Condition, .) %>%
    unique() %>%
    sample_0[.]
  
  profile_merge <- profile[,,sample_0]
  names(assays(profile_merge)) <- names(sample_0)
  rownames(sampleData(profile_merge)) <- names(assays(profile_merge))
  
  for(i in seq_along(sample_0)) {
    # print(sample_0[i])
    # print(names(sample_0[i]))
    sample_name <- sample_0[i]
    group_name <- names(sample_0[i])
    id_rep <- filter(sample_data, Condition == group_name) %>%
      pull(id)
    assay_rep <- profile[,,id_rep] %>%
      assays() %>%
      as.list()
    assays(profile_merge)[[group_name]] <- Reduce("+", assay_rep)/length(assay_rep)
  }
  
  return(profile_merge)
}



plot_profile <- function(profile, name='YAP_all', 
                         color=c("white","#e41a1c"),
                         genes_label=c('CCN1','CCN2','ANKRD1'),
                         w=8
                         ) {
  # show all sample
  pdf(str_glue('{Sys.Date()}_consensus_peaks_profile_{name}.pdf'), 
      width=w, height=7)
  generateEnrichedHeatmap(
    profile,
    include_group_annotation=T,
    matrices_axis_name=c("-3k", "0", "3k"),
    group_anno_color=RColorBrewer::brewer.pal(length(unique(mcols(profile)[['Type']])), "Dark2"),
    all_color_scales_equal = T,
    # color_by_sample_group='Factor',
    # matrices_color = as.list(rep(list(c("white","red2")), nrow(profile@sampleData))),
    # matrices_color=list(YAP=c("white","#e41a1c"),
    #                     TAZ=c("white","#377eb8"),
    #                     panTEAD=c("white","#4daf4a"),
    #                     H3K27AC=c("white","#984ea3")),
    matrices_color = color,
    genes_to_label=genes_label,
    use_raster=T,
    raster_quality=2,
    return_ht_list=F)
  while(dev.cur() != 1) dev.off()
}



yap <- get_profile('YAP')
plot_profile(yap$rpm, 'YAP_rpm_all', color=c("white","#e41a1c"))
plot_profile(yap$raw, 'YAP_raw_all', color=c("white","#e41a1c"))
# keep part of samples
yap_profile <- yap$raw
sampleData(yap_profile) %>% rownames()
keep_sample <- setdiff(seq_len(nrow(sampleData(yap_profile))), 3)
yap_profile <- yap_profile[,,keep_sample]
# merge samples
yap_profile <- merge_profile(yap_profile)
plot_profile(yap_profile, 'YAP_raw_subset_merged', color=c("white","#e41a1c"), w=5)


taz <- get_profile('TAZ')
plot_profile(taz$rpm, 'TAZ_rpm_all', color=c("white","#377eb8"))
plot_profile(taz$raw, 'TAZ_raw_all', color=c("white","#377eb8"))
# keep part of samples
taz_profile <- taz$raw
sampleData(taz_profile) %>% rownames()
keep_sample <- setdiff(seq_len(nrow(sampleData(taz_profile))), 3)
taz_profile <- taz_profile[,,keep_sample]
# merge samples
taz_profile <- merge_profile(taz_profile)
plot_profile(taz_profile, 'TAZ_raw_subset_merged', color=c("white","#377eb8"), w=5)


tead <- get_profile('panTEAD')
plot_profile(tead$rpm, 'panTEAD_rpm_all', color=c("white","#4daf4a"))
plot_profile(tead$raw, 'panTEAD_raw_all', color=c("white","#4daf4a"))
# keep part of samples
tead_profile <- tead$raw
sampleData(tead_profile) %>% rownames()
keep_sample <- setdiff(seq_len(nrow(sampleData(tead_profile))), 1)
tead_profile <- tead_profile[,,keep_sample]
# merge samples
tead_profile <- merge_profile(tead_profile)
plot_profile(tead_profile, 'panTEAD_raw_subset_merged', color=c("white","#4daf4a"), w=5)


h3k27ac <- get_profile('H3K27AC')
plot_profile(h3k27ac$rpm, 'H3K27AC_rpm_all', color=c("white","#984ea3"))
plot_profile(h3k27ac$raw, 'H3K27AC_raw_all', color=c("white","#984ea3"))
# keep part of samples
h3k27ac_profile <- h3k27ac$raw
sampleData(h3k27ac_profile) %>% rownames()
keep_sample <- setdiff(seq_len(nrow(sampleData(h3k27ac_profile))), 3)
h3k27ac_profile <- h3k27ac_profile[,,keep_sample]
# merge samples
h3k27ac_profile <- merge_profile(h3k27ac_profile)
plot_profile(h3k27ac_profile, 'H3K27AC_raw_subset_merged', color=c("white","#984ea3"), w=5)






x <- mcols(db_peaks_yap) %>% as_tibble()
sig_index <- which(abs(x$Fold) > 0.58 & x$FDR < 0.05)





