
library(tidyverse)
library(GenomicAlignments)
library(GenomicRanges)
library(GenomeInfoDb)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)   # hs
library(TxDb.Mmusculus.UCSC.mm39.knownGene)  # mm
library(TxDb.Celegans.UCSC.ce11.refGene)     # celegans
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene) # fly
library(rtracklayer)
library(furrr)
options(future.globals.maxSize=2 * 1024^3)  # 2GB if needed


# config ------------
setwd("D:\\EED_CutTag\\cel_08-25\\result")
# setwd("F:/workspace/J009_IMR_0_0.3_0.6/result")

txdb <- TxDb.Celegans.UCSC.ce11.refGene
# txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

bam_dir <- "03_bam"
peak_dir <- "04_peaks"
peak_pattern <- "broadPeak$"
consensus_peakfile <- "04_peaks/consensus_H3K27me3.bed"


peak_files <- list.files(peak_dir, pattern=peak_pattern, full.names=TRUE, recursive=FALSE)

df <- tibble(
  peak_file=peak_files,
  sample=basename(peak_files) %>% str_replace("_peaks.*Peak$", ""),
  bam_file=str_c(bam_dir, "/", sample, ".filtered.bam")
) |> print()



# genomic data ------------
gene <- genes(txdb) %>%
  keepStandardChromosomes(pruning.mode="coarse") %>%
  GenomicRanges::reduce() |> 
  sort()
seqlevelsStyle(gene) <- "UCSC"
seqlevels(gene)

gene_gr <- gene %>%
  # resize(width=width(.) + 2000, fix="center") %>%
  GenomicRanges::reduce()

tss_gr <- promoters(gene, 1000, 1000) %>%
  GenomicRanges::reduce()

consensus_gr <- rtracklayer::import(consensus_peakfile) %>%
  keepStandardChromosomes(pruning.mode="coarse") |> 
  GenomicRanges::reduce()
seqlevelsStyle(consensus_gr) <- "UCSC"
seqlevels(consensus_gr)


# functions ------------
calculate_ratio_separated <- function(
    target_gr, query_gr_list, sample_ratio=0.25
) {
  
  res <- map(seq_along(query_gr_list), \(idx) {
    # idx <- 5
    roi <- query_gr_list[[idx]]
    roi_name <- names(query_gr_list)[idx]
    
    if(sample_ratio==0)  roi_sampled <- roi else 
      roi_sampled <- sample(roi, floor(length(roi) * sample_ratio), replace=FALSE) |> sort()
    
    roi_sampled_L <- flank(roi_sampled, width=floor(1.0*width(roi_sampled)), 
                           start=TRUE, both=FALSE, ignore.strand=T) |> 
      sort()
    roi_sampled_R <- flank(roi_sampled, width=floor(1.0*width(roi_sampled)), 
                           start=FALSE, both=FALSE, ignore.strand=T) |> 
      sort()
    
    cnt <- list(roi_sampled, roi_sampled_L, roi_sampled_R) |> 
      set_names(c('_IN', '_L','_R')) |> 
      map(\(x) sum(countOverlaps(target_gr, x, ignore.strand=TRUE))) %>% 
      list_c()
    
    enframe(cnt, name="region", value="count") %>%
      pivot_wider(names_from=region, values_from=count) |> 
      dplyr::rename_with(\(x) str_c(roi_name, x))
  }  )
  
  return(purrr::reduce(res, cbind))
}


get_sample_stat <- function(
    bam_path, peak_path,
    sample_times=100, sample_ratio=0.25) {
  
  # Read BAM ONCE per sample
  bam_gr <- GenomicAlignments::readGAlignmentPairs(bam_path) %>% 
    granges() |> 
    keepStandardChromosomes(pruning.mode="coarse") |> 
    sort()
  seqlevelsStyle(bam_gr) <- "UCSC"
  
  # peak
  peak_gr <- rtracklayer::import(peak_path) %>% 
    keepStandardChromosomes(pruning.mode="coarse")
  seqlevelsStyle(peak_gr) <- "UCSC"
  
  query_gr_list = lst(bam=bam_gr, peak=peak_gr, gene=gene_gr, tss=tss_gr, consensus=consensus_gr)
  
  res <- future_map(
    seq_len(sample_times), \(x) calculate_ratio_separated(
      bam_gr, query_gr_list, sample_ratio)
  )
  
  return(list_rbind(res, names_to="sample_times"))
}



# run ------------
plan(multisession, workers=min(6, nrow(df))) 

stats <- list()
for(idx in seq_len(nrow(df))) {
  print(str_c('processing: ', df$sample[idx]))
  stats[[df$sample[idx]]] <- get_sample_stat(
    df$bam_file[idx], df$peak_file[idx], sample_times=100, sample_ratio=0.25)
}

plan(sequential)




# plot ------------
stats_tidy <- stats %>%
  list_rbind(names_to="sample") %>% #glimpse()
  mutate(
    random = bam1 * 2 / (bam2 + bam3),
    consensus = consensus1 * 2 / (consensus2 + consensus3),
    peak = peak1 * 2 / (peak2 + peak3),
    gene = gene1 * 2 / (gene2 + gene3),
    tss = tss2 * 2 / (tss2 + tss3)
  ) |> 
  as_tibble()
summary_stats <- stats_tidy %>%
  summarise(across(everything(), mean), .by=sample) %>% 
  glimpse()

writexl::write_xlsx(stats_tidy, glue::glue("{Sys.Date()}_stats_tidy.xlsx"))



# Plot
metrics <- colnames(stats_tidy)[18:21]
plots <- map(metrics, ~ {
  stats_tidy %>%
    separate_wider_delim(sample, delim="_", names=c("group","rep"), cols_remove=F) %>% 
    # mutate(sample_split=str_split(sample, "_", simplify=TRUE)) %>%
    # mutate(group=ifelse(ncol(sample_split) >= 2, sample_split[,1], sample)) %>%
    ggplot(aes(sample, .data[[.x]], fill=group)) +
    geom_violin() +
    geom_boxplot(width=0.3, outlier.shape=NA, alpha=0.7) +
    labs(title=.x) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=90, hjust=1))
})

patchwork::wrap_plots(plots, ncol=1) %>%
  ggsave(glue::glue("{Sys.Date()}_stats_tidy.pdf"), ., 
         width=8, height=3 * length(plots))





