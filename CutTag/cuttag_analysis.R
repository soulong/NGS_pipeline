
library(tidyverse)
library(GenomicRanges)
library(ChIPseeker)
library(chromVAR)
library(EnrichedHeatmap)
library(circlize)
library(DESeq2)



setwd('C:/Users/haohe/Desktop/CUTTAG/results')
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
annodb <- "org.Hs.eg.db"



## load peak -----
# bam files
bam_files <- '02_alignment/bowtie2/target/markdup' %>% 
  list.files(pattern='.bam$', full.names=T) %>% 
  set_names(., nm=map_chr(., basename)) %>% 
  print()

# bw files
bw_files <- '03_peak_calling/03_bed_to_bigwig' %>% 
  list.files(pattern='.bigWig$', full.names=T) %>% 
  set_names(., nm=map_chr(., basename)) %>% 
  print()
  
# peak files
peak_files <- '03_peak_calling/04_called_peaks/macs2' %>% 
  list.files(pattern='.broadPeak$', full.names=T) %>% 
  set_names(., nm=map_chr(., basename)) %>% 
  print()


# import peaks
peaks <- peak_files %>% 
  map(rtracklayer::import) %>% 
  GRangesList() %>% 
  print()
seqlevelsStyle(peaks) <- "UCSC"
genome(peaks) <- "hg38"
# print(standardChromosomes(peaks))

# keep autosomes only
peaks <- peaks %>% 
  keepStandardChromosomes(pruning.mode="tidy")
print(peaks[[1]]@seqnames)



## consensus peak -----
# get consensus peak
consensus_peak <- peaks %>% 
  GRangesList() %>% 
  unlist() %>% 
  GenomicRanges::reduce()

# annotate consensus peak
consensus_peak_anno <- consensus_peak %>% 
  ChIPseeker::annotatePeak(tssRegion=c(-3000, 3000),
                           TxDb=txdb, 
                           level='gene',
                           annoDb=annodb)
consensus_peak <- consensus_peak_anno@anno

# peak within TSS
consensus_peak_within_tss <- 
  (abs(consensus_peak$distanceToTSS) <= 3000) %>% 
  consensus_peak[.]
# rtracklayer::export(consensus_peak_within_tss, 'consensus_peak_within_tss.bed')




## visualize TSS -----
# check TSS
tss <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tss_hit <- findOverlaps(tss, consensus_peak, minoverlap=1000)
tss_overlap_with_peak <- tss[unique(queryHits(tss_hit))]

# mat <- normalizeToMatrix(peaks[[1]], 
#                          tss_overlap_with_peak[1:100], 
#                          extend = 3000,
#                          value_column="score", 
#                          mean_mode = "w0", 
#                          w = 50,
#                          keep = c(0, 0.99))
# EnrichedHeatmap(mat)

# loop over samples
mat_list <- peaks %>% 
  map(\(x) 
      normalizeToMatrix(
        signal=x, 
        target=resize(tss_overlap_with_peak, width=1, fix="center"), 
        extend=3000, 
        value_column="score", 
        mean_mode="w0", w=50, keep=c(0, 0.99)) )

# set global scale bar
all_values <- unlist(lapply(mat_list, function(m) as.vector(m)))
global_min <- quantile(all_values, 0.01)
global_max <- quantile(all_values, 0.99)
col_fun <- colorRamp2(
  breaks = c(global_min, global_max),
  colors = c("white", "#a50026"))

# set same ylim for top annotation
all_ylim <- lapply(mat_list, function(mat) {
  # Extract the mean enrichment profile for this sample
  avg <- colMeans(mat, na.rm = TRUE)  # or colMeans if you used target_ratio
  c(min(avg), max(avg))
})
global_ymin <- min(sapply(all_ylim, `[`, 1)) * 0
global_ymax <- max(sapply(all_ylim, `[`, 2)) * 1.2
# rebuild the top_annotation with fixed ylim
fixed_anno <- HeatmapAnnotation(
  lines = anno_enriched(
    ylim = c(global_ymin, global_ymax),   # <-- THIS IS THE KEY
    gp = gpar(col = "black", lwd = 2),
    axis_param = list(
      at = c(global_ymin, 0, global_ymax),
      labels = c(round(global_ymin, 2), "0", round(global_ymax, 2)),
      side = "left"
    )),
  height = unit(2, "cm")  # make taller for visibility
)

# cluster only once on the selected matrix
row_ord <- mat_list[[13]] %>% 
  EnrichedHeatmap() %>% 
  row_order()

# create heatmap over samples
ht_list <- NULL
for(i in seq_along(mat_list)) {
  mat <- mat_list[[i]]
  sample_name <- names(mat_list)[i] %>% 
    str_split_1('[.]') %>% 
    .[1]
  
  ht <- EnrichedHeatmap(
    mat,
    name = sample_name,
    col = col_fun,
    axis_name = c("-3kb", "TSS", "3kb"),
    column_title = sample_name,
    top_annotation = fixed_anno,
    show_heatmap_legend = ifelse(i==1, T, F),
    heatmap_legend_param = list(
      title = "normalized",
      legend_direction = "horizontal",
      at = c(global_min, global_max),
      labels = c(round(global_min, 2), round(global_max, 2))),
    use_raster=TRUE
  )
  if(is.null(ht_list)) {
    ht_list <- ht
  } else ht_list <- ht_list + ht
}


# draw and save
ht_opt$message = FALSE
pdf("multi_sample_peak_heatmap_roword.pdf", width = 20, height = 8)
draw(ht_list, 
     row_order = row_ord,
     # column_split = sample_groups,
     # gap = unit(5, "mm"),
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     padding = unit(c(10, 10, 10, 10), "mm"),
     use_raster = TRUE)
while(dev.cur() != 1) dev.off()

pdf("multi_sample_peak_heatmap.pdf", width = 20, height = 8)
draw(ht_list, 
     row_order = NULL,
     # column_split = sample_groups,
     # gap = unit(5, "mm"),
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     padding = unit(c(10, 10, 10, 10), "mm"),
     use_raster = TRUE)
while(dev.cur() != 1) dev.off()





## DEG peak ----
# count over bam files
fragment_counts <- getCounts(bam_files, 
                             format = "bam",
                             consensus_peak_within_tss, 
                             paired = TRUE, by_rg = FALSE)
colnames(fragment_counts) <- colnames(fragment_counts) %>% 
  str_replace_all('(.target.markdup).*','')

count_mat <- assay(fragment_counts)
dim(count_mat)
# remove low count genes
selected_row <- which(rowSums(count_mat) > 10)
count_mat_selected = count_mat[selected_row, ]
dim(count_mat_selected)

col_data <- tibble(sample=colnames(fragment_counts),
                   condition=c(rep('HU',3), rep('IMR', 6), rep('IMRJ009',6))) %>% 
  column_to_rownames('sample') %>% 
  print()

dds <- DESeqDataSetFromMatrix(
  countData = count_mat_selected,
  colData = col_data, design = ~ condition)
dds <- DESeq(dds)
# count_norm <- counts(dds, normalized = TRUE)
count_vst <- vst(dds, blind=FALSE)

pdf("DESeq2_pca.pdf", width=7, height=7)
DESeq2::plotPCA(count_vst, "condition", ntop=2000) + 
  ggrepel::geom_label_repel(aes(label=name), size=2, max.overlaps=20) +
  theme_bw()
while(dev.cur() != 1) dev.off()





