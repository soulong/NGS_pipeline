library(GenomicAlignments)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(tidyverse)
library(DESeq2)

setwd("F:/workspace/2026-01-16_EED_NMN_LH/IMR/result")

source("F:/GitHub/myScript/NGS/peak_process.R")


samples <- crossing(c("0_","05_"), c("1","2","3")) %>% 
  unite("uni", sep="") %>% pull() %>% print()
bam_files <- str_c("03_bam/IMR", samples, ".filtered.bam") %>% 
  print()
bw_files <- str_c("03_bam/IMR", samples, ".CPM.bw") %>% 
  print()

gr  <- import("04_peaks/consensus_H3K27me3.bed")



# get coverage over bam
result <- list()
for(idx in seq_along(bam_files)) {
  # bam <- bam_files[idx]
  res <- get_region_and_flanks_coverage(
    bam_path   = bam_files[idx],
    gr         = gr,
    flank_size = 3000,
    mode       = "coverage"          # or "coverage"
  )
  result[[samples[idx]]] <- res
}

# get count matrx
df <- list_rbind(result, names_to="group") %>% 
  as_tibble() %>% 
  separate_wider_delim(group, "_", names=c("group","rep")) %>% 
  print()

# merge replicates
df_merge <- df %>% 
  reframe(across(starts_with("count_"), sum),
          across(starts_with("cov_"), mean),
          .by=c(group, seqnames:width)) %>% 
  pivot_wider(names_from=group, values_from=6:last_col()) %>% 
  print()
df_merge <- df_merge %>% 
  mutate(cov_ratio_0=(cov_region_0/(cov_left_0 + cov_right_0)),
         cov_ratio_05=(cov_region_05/(cov_left_05 + cov_right_05))) %>% 
  mutate(relative_diff=(cov_ratio_05) / cov_ratio_0) %>% 
  arrange(desc(relative_diff)) %>% 
  print()
# export(df_merge %>% GRanges(), str_glue("{Sys.Date()}_df_merge.gtf"))


# diff analysis
df_peak_count_wide <- df %>% 
  pivot_wider(id_cols=seqnames:width, names_from=c(group, rep), values_from=count_region) %>% 
  unite(uid, seqnames:width, sep='_') %>% 
  print()

col_data <- dplyr::select(df, group:rep) %>% 
  unite('uid', everything(), sep='_', remove=F) %>% 
  distinct() %>% 
  column_to_rownames('uid') %>% 
  print()

dds <- DESeqDataSetFromMatrix(column_to_rownames(df_peak_count_wide, 'uid'),
                              col_data, design= ~ group)
dds <- DESeq(dds)
dds_result <- results(dds, contrast=c("group", "05", "0")) %>% 
  as_tibble(rownames='uid') %>% 
  arrange(pvalue) %>% 
  print()
dds_result_sig <- dds_result %>% 
  dplyr::filter(padj < 0.05) %>% 
  print()
dds_result_sig_grange <- dds_result_sig %>% 
  dplyr::select(uid) %>% 
  separate(uid, c("seqnames","start","end","width"), sep="_") %>% 
  GRanges() %>% 
  print()




# peak heatmap
gr_use <- list(sig=dds_result_sig_grange,
               all=df_merge %>% GRanges
)
for (idx in seq_along(gr_use)) {
  # hp <- heatmap_profile(
  #   bw_files[c(1,4)],
  #   gr,
  #   mode='reference_point',
  #   upstream=3000, downstream=3000, bin_size=100,
  #   colcor_scales=c(-1, 0, 1),
  #   kmeans=NA,
  #   sort_by='mean',
  #   scale='row',
  #   return_granges=T
  # )
  # # hp$heatmap
  # pdf(str_glue("{Sys.Date()}_refpoint_scaleby_row.pdf"), width=7, height=9)
  # print(hp$heatmap)
  # dev.off()
  
  hp2 <- heatmap_profile(
    bw_files[c(1,4)],
    gr_use[[idx]],
    mode='reference_point',
    reference_point='center',
    upstream=3000, downstream=3000, bin_size=100,
    color_scales=c(0, NA),
    colors=c("white","#295072"),
    kmeans=NA,
    sort_by='mean',
    scale='none',
    return_granges=T
  )
  # hp$heatmap
  pdf(str_glue("{Sys.Date()}_refpoint_scaleby_none_{names(gr_use)[idx]}.pdf"), width=7, height=9)
  print(hp2$heatmap)
  dev.off()
}


# gr <- df_merge %>% as_tibble() %>% 
#   mutate(type=as.factor(sign(log2(as.numeric(relative_diff) + 1e-8)))) %>% 
#   dplyr::filter(abs(log2(as.numeric(relative_diff) + 1e-8)) > 1) %>% 
#   # dplyr::filter((as.numeric(cov_region_DMSO) + as.numeric(cov_region_IMR)) < 50) %>%
#   # dplyr::slice_sample(n=2000, replace=F) %>% 
#   # dplyr::slice_max(relative_diff, n=x, with_ties=F) %>%
#   # dplyr::slice_min(relative_diff, n=10000, with_ties=F) %>%
#   GRanges() %>% list()

for(x in c(2000, 10000, 50000)) {
  
  gr <- df_merge %>% as_tibble() %>% 
    # dplyr::filter((as.numeric(cov_region_DMSO) + as.numeric(cov_region_IMR)) < 50) %>%
    # dplyr::slice_sample(n=2000, replace=F) %>% 
    dplyr::slice_max(relative_diff, n=x, with_ties=F) %>%
    # dplyr::slice_min(relative_diff, n=10000, with_ties=F) %>%
    GRanges() %>% list()
  
  
  # hp <- heatmap_profile(
  #   bw_files[c(1,4)],
  #   gr,
  #   mode='reference_point',
  #   upstream=3000, downstream=3000, bin_size=100,
  #   colcor_scales=c(-1, 0, 1),
  #   kmeans=NA,
  #   sort_by='mean',
  #   scale='row',
  #   return_granges=T
  # )
  # # hp$heatmap
  # pdf(str_glue("{Sys.Date()}_refpoint_scaleby_row_{x}.pdf"), width=7, height=9)
  # print(hp$heatmap)
  # dev.off()
  
  hp2 <- heatmap_profile(
    bw_files[c(1,4)],
    gr,
    mode='reference_point',
    upstream=3000, downstream=3000, bin_size=100,
    color_scales=c(0, NA),
    colors=c("white","#295072"),
    kmeans=NA,
    sort_by='mean',
    scale='none',
    return_granges=T
  )
  # hp$heatmap
  pdf(str_glue("{Sys.Date()}_refpoint_scaleby_none_{x}.pdf"), width=7, height=9)
  print(hp2$heatmap)
  dev.off()
}



