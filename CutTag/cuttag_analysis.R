
library(tidyverse)
library(GenomicRanges)
library(ChIPseeker)
library(GenomeInfoDb)
library(EnrichedHeatmap)
ht_opt$message = FALSE
library(circlize)
library(DESeq2)
require(org.Hs.eg.db)

# load helper
parent_path <- rstudioapi::getActiveDocumentContext()$path |> 
  dirname()|> dirname()
source(file.path(parent_path, 'ngs_helper.R'))


data_dir <- norm_path(r"(E:\NGS\2026-01-08_CutTag_293T_ZTX\result)") %>% 
  setwd()

annodb <- "org.Hs.eg.db"

set_style <- function(x) { 
  seqlevelsStyle(x) <- "UCSC"
  x <- keepStandardChromosomes(x, pruning.mode="tidy")
  return(x)
  }

# txdb <- TxDb.Hsapiens.UCSC.T2T-CHM13v2.knownGene::TxDb.Hsapiens.UCSC.T2T-CHM13v2.knownGene
txdb_file <- 'F:/index/hs/chm13/t2t_chm13v2_txdb.rds'
if(file.exists(txdb_file)) {
  txdb <- read_rds(txdb_file)
} else {
  # t2t_gtf <- norm_path(r"(F:\index\hs\chm13\chm13.draft_v2.0.gene_annotation.gff3.gz)")
  t2t_gtf <- norm_path(r"(F:\index\hs\chm13\hs1.ncbiRefSeq.gtf.gz)")
  
  # gtf <- rtracklayer::import(t2t_gtf)
  # gtf_df <- as_tibble(gtf)
  # gtf_df_clean <- gtf_df %>% 
  #   dplyr::select(seqnames:type, phase, 
  #                 gene_biotype,
  #                 gene_name, 
  #                 gene_id=source_gene,
  #                 transcript_biotype,
  #                 transcript_name=source_transcript_name,
  #                 transcript_id=source_transcript
  #                 )
  # gtf_df_clean_gr <- gtf_df_clean %>% 
  #   makeGRangesFromDataFrame(keep.extra.columns=T) %>% 
  #   print()
  # clean_gtf_file <- 'F:/index/hs/chm13/chm13v2_clean.gtf'
  # rtracklayer::export(gtf_df_clean_gr, clean_gtf_file)
  
  txdb <- txdbmaker::makeTxDbFromGFF(t2t_gtf)
  seqlevelsStyle(txdb) <- 'UCSC'
  seqlevels(txdb)
  write_rds(txdb, txdb_file)
}



# prepare peak ------------------------------------------------

## read peak ------------------------------------------------
peak_suffix <- '.narrowPeak'
peak_file <- '04_peaks' %>% 
  list.files(pattern=peak_suffix, full.names=T) %>% 
  set_names(., nm=map_chr(., \(x) basename(x) %>% 
                            str_replace(fixed(peak_suffix),'') %>% 
                            str_replace(fixed('_peaks'),''))) %>% 
  print()

# import 
peak <- peak_file %>% 
  map(rtracklayer::import) %>% 
  GRangesList() %>% 
  set_style() %>% 
  print()


### consensus ------------------------------------------------
consensus <- peak %>% 
  # subset samples
  .[which(str_detect(names(.), 'TEAD'))] %>% 
  # , group=names(.) %>% str_replace_all("_\\d","")
  compute_consensus(., min_occurrence=2, min_gapwidth=1)

consensus <- consensus$all
# remove small fragment
consensus <- consensus[width(consensus) > 20]

# # annotate consensus
# anno <- consensus %>%
#   ChIPseeker::annotatePeak(tssRegion=c(-3000, 3000),
#                            TxDb=txdb,
#                            level='gene',
#                            annoDb=annodb)
# consensus <- anno@anno
# # peak within TSS
# consensus_within_tss <- 
#   (abs(consensus$distanceToTSS) <= 3000) %>% 
#   consensus_peak[.] %>% 
#   print()
# rtracklayer::export(consensus_peak_within_tss, 'consensus_peak_within_tss.bed')


### tss overlap ------------------------------------------------
# tss <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
# 
# hit <- findOverlaps(tss, consensus, minoverlap=1000)
# tss_overlap_with_consensus <- tss[unique(queryHits(hit))]
# 
# hit <- findOverlaps(consensus, tss, minoverlap=1000)
# consensus_overlap_with_tss <- consensus[unique(queryHits(hit))]



## read bam ------------------------------------------------
bam_suffix <- '.filtered.bam'
bam_file <- '03_bam' %>% 
  list.files(pattern=bam_suffix, full.names=T) %>%
  str_subset(".bai", negate=T) %>% 
  set_names(., nm=map_chr(., \(x) basename(x) %>% 
                            str_replace(fixed(bam_suffix),''))) %>% 
  # subset samples
  .[which(str_detect(names(.), 'BRD4'))] %>% 
  print()

# # get count over bam
# bam_count <- list()
# for(f in bam_file) {
#   # bam <- bam_file[idx]
#   res <- compute_count(f, consensus, flank_size=0, mode="count")
#   bam_count[[f]] <- res$count
# }
# count <- list_rbind(bam_count, names_to="group") %>% 
#   as_tibble() %>% 
#   unite('uid', seqnames:end) %>% 
#   dplyr::select(uid, group, count_region) %>% 
#   pivot_wider(id_cols=uid, 
#               names_from=group, 
#               names_prefix='03_bam/', 
#               values_from=count_region) %>% 
#   rename_with(\(x) str_replace_all(x, '03_bam/', '') %>% 
#                 str_replace_all('.filtered.bam', '')) %>% 
#   column_to_rownames("uid")
# get count over bam
bam_count <- chromVAR::getCounts(
  bam_file, format = "bam", consensus,
  paired = TRUE, by_rg = FALSE)
count <- assay(bam_count, 'counts')
rownames(count) <- as_tibble(consensus) %>% 
  unite('uid', 1:3) %>% pull(uid)
colnames(count) <- names(bam_file)

count[1:2,1:2]


### deseq2 ------------------------------------------------
# coldata
col_data <- colnames(count) %>% 
  tibble(sample=.) %>% 
  mutate(group=sapply(strsplit(sample, "_"), `[`, 1)) %>% 
  column_to_rownames("sample") %>% 
  print()
# deseq2
dds <- DESeqDataSetFromMatrix(
  count, col_data, design = ~ group) %>% 
  DESeq()
dds_result <- results(dds, contrast=c("group", "P62", "DMSO")) %>% 
  as_tibble(rownames='uid') %>% 
  separate(uid, c("seqnames","start","end","width"), sep="_") %>% 
  arrange(pvalue) %>% 
  print()

# anno deg
dds_result_gr_anno <- dds_result %>% 
  makeGRangesFromDataFrame(keep.extra.columns=T) %>%
  ChIPseeker::annotatePeak(tssRegion=c(-3000, 3000),
                           TxDb=txdb,
                           level='gene',
                           annoDb=annodb) %>% 
  { .@anno } %>% 
  # trim() %>% # fix out of genome range problem
  print()

# tidy deg
yap_targets <- c(
  "CTGF","CCN2","ANKRD1","CYR61","CCN1","TOP2A","KIF14","CCNA2",
  "CDCA8","CENPF","KIF23","KIF20B","KNTC1","RRM2",
  "MCM3","SGOL1.AS1","TUBB","MYBL1","RAD18",
  "ZWILCH","SGOL1","TIMELESS","GINS1","SMC3","TK1",
  "MRE11A","MCM7","SUV39H2","GADD45B","FOSL1","CENPV",
  "RUVBL2","MYC","GLI2","AXL","ABCB1","CAT","GPATCH4",
  "LMNB2","TXN","WSB2","AREG","FOXF2","IGFBP3","RASSF2",
  "AMOTL2","NPPB","CCND1")
deg_tidy <- as_tibble(dds_result_gr_anno) %>% 
  janitor::clean_names() %>% 
  mutate(yap_targets=ifelse(gene_id %in% yap_targets, gene_id, NA)) %>% 
  mutate(type=ifelse(
    log2fold_change > log2(1.5) & padj < 0.05, 'pos',
    ifelse(log2fold_change < log2(1/1.5) & padj < 0.05, 'neg', NA))) %>% 
  mutate(label=ifelse(!is.na(yap_targets) & type %in% 'pos', yap_targets, NA)) %>% 
  glimpse()
# save
writexl::write_xlsx(
  deg_tidy, str_glue('{Sys.Date()}_valcano_consensus_tead_deseq2_deg.xlsx'))

# valcono
p <- deg_tidy %>% 
  ggplot(aes(log2fold_change, -log10(pvalue))) +
  geom_point(aes(color=type), alpha=0.7, show.legend=F) +
  scale_color_discrete(direction = -1) +
  ggrepel::geom_text_repel(aes(label=label), size=5, show.legend=F) +
  theme_bw()
p <- ggrastr::rasterize(p, dpi=600)
ggsave(str_glue('{Sys.Date()}_valcano_consensus_tead_deseq2_deg.pdf'),
       p, width=5, height=5)



## read bw ------------------------------------------------
bw_suffix <- '.CPM.bw'
bw_file <- '03_bam' %>% 
  list.files(pattern=bw_suffix, full.names=T) %>% # bw
  set_names(., nm=map_chr(., \(x) basename(x) %>% 
                            str_replace(fixed(bw_suffix),''))) %>% 
  # subset samples
  .[which(str_detect(names(.), 'TEAD|BRD4'))] %>% 
  print()

# reorder as heatmap needs
new_order <- names(bw_file)[order(factor(sapply(strsplit(names(bw_file), "_"), `[`, 2)))]
bw_file <- bw_file[new_order]

# read bw
bw <- bw_file %>% 
  map(rtracklayer::import) %>% 
  GRangesList() %>% 
  set_style() %>% 
  print()




# heatmap ------------------------------------------------

# signal of interest
soi <- bw %>% .[which(str_detect(names(.), 'TEAD|BRD4'))]
rois <- dplyr::filter(deg_tidy, !is.na(type)) %>% 
  split(., .[['type']]) %>% 
  map(\(x) makeGRangesFromDataFrame(x, keep.extra.columns=T)) %>% 
  print()
  
for(idx in seq_along(rois)) {
  # region of interest
  roi <- rois[[idx]]
  # calculate mat over signal
  mat_list <- compute_signal_matrix(
    signal=soi, 
    region=roi, 
    mode='reference_point', reference_point="center",
    bin_size=50, scale='none'
  )
  # plot
  result <- heatmap_profile(mat_list,
                            color_scales = c(0, NA), 
                            colors = c("white", "#b12923"))
  # save
  f_name <- str_glue('consensus_tead_and_deg_{names(rois)[idx]}')
  pdf(str_glue("{Sys.Date()}_heatmap_{f_name}.pdf"), 
      width = length(result$heatmap)*1.6, height = 8)
  print(result$heatmap)
  while(dev.cur() != 1) dev.off()
}




# 
# 
# 
# 
# # calculate mat over signal
# mat_list <- soi %>% 
#   map(\(x) 
#       normalizeToMatrix(
#         signal=x, 
#         target=resize(roi, width=1, fix=pos), 
#         extend=3000, 
#         value_column="score", 
#         mean_mode="w0", w=50, keep=c(0, 0.99))
#       )
# 
# p <- heatmap_profile(
#   mat_list, roi, kmeans=2, return_granges=F)
# 
# # set global scale bar
# all_values <- unlist(lapply(mat_list, function(m) as.vector(m)))
# global_min <- quantile(all_values, 0.01)
# global_max <- quantile(all_values, 0.99)
# if(global_min < 0 & global_max > 0) {
#   col_fun <- colorRamp2(
#     breaks = c(global_min, 0, global_max),
#     colors = c("#295072", "white", "#b12923"))
# } else {
#   col_fun <- colorRamp2(
#     breaks = c(global_min, global_max),
#     colors = c("white", "#a50026"))
# }
# 
# # set same ylim for top annotation
# all_ylim <- lapply(mat_list, function(mat) {
#   # Extract the mean enrichment profile for this sample
#   avg <- colMeans(mat, na.rm = TRUE)  # or colMeans if you used target_ratio
#   c(min(avg), max(avg))
# })
# global_ymin <- min(sapply(all_ylim, `[`, 1)) * 0
# global_ymax <- max(sapply(all_ylim, `[`, 2)) * 1.2
# # rebuild the top_annotation with fixed ylim
# fixed_anno <- HeatmapAnnotation(
#   lines = anno_enriched(
#     ylim = c(global_ymin, global_ymax),   # <-- THIS IS THE KEY
#     gp = gpar(col = "black", lwd = 2),
#     axis_param = list(
#       at = c(global_ymin, 0, global_ymax),
#       labels = c(round(global_ymin, 2), "0", round(global_ymax, 2)),
#       side = "left"
#     )),
#   height = unit(2, "cm")  # make taller for visibility
# )
# # create heatmap over signal
# ht_list <- NULL
# for(i in seq_along(mat_list)) {
#   mat <- mat_list[[i]]
#   sample_name <- names(mat_list)[i] %>% 
#     str_split_1('[.]') %>% 
#     .[1]
#   # plot
#   ht <- EnrichedHeatmap(
#     mat,
#     name = sample_name,
#     col = col_fun,
#     axis_name = c("-3kb", "TSS", "3kb"),
#     column_title = sample_name,
#     top_annotation = fixed_anno,
#     show_heatmap_legend = ifelse(i==1, T, F),
#     heatmap_legend_param = list(
#       title = "normalized",
#       legend_direction = "horizontal",
#       at = c(global_min, global_max),
#       labels = c(round(global_min, 2), round(global_max, 2))),
#     use_raster=TRUE
#   )
#   if(is.null(ht_list)) {
#     ht_list <- ht
#   } else ht_list <- ht_list + ht
# }
# 
# # save plot
# pdf(str_glue("{Sys.Date()}_heatmap_{vis_name}.pdf"), 
#     width = length(ht_list)*1.6, height = 8)
# draw(ht_list, 
#      row_order = NULL,
#      # column_split = sample_groups,
#      # gap = unit(5, "mm"),
#      heatmap_legend_side = "bottom",
#      annotation_legend_side = "bottom",
#      padding = unit(c(10, 10, 10, 10), "mm"),
#      use_raster = TRUE)
# while(dev.cur() != 1) dev.off()



