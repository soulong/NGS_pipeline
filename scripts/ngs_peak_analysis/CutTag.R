
library(tidyverse)
library(GenomicAlignments)
library(GenomicRanges)
library(plyranges)
library(rtracklayer)
library(ChIPseeker)
library(clusterProfiler)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(org.Hs.eg.db)
library(edgeR)
# library(ggpointdensity)
# library(lemon)
library(EnrichedHeatmap)
library(circlize)
library(rlist)
library(profileplyr)
# library(soGGi)
library(NbClust)
library(RColorBrewer)
library(BiocParallel)


rstudioapi::getActiveDocumentContext()$path %>%
  dirname() %>% setwd()

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb_for_annotateRanges <- "hg38"
species <- "Homo_sapiens"
PE_reads <- TRUE



# get peak files
peak_files <- list.files(pattern = ".narrowPeak", recursive = T, full.names = T)
names(peak_files) <- map_chr(peak_files, ~ str_split(.x, "/", simplify = T)[4])

# order
peak_files <- peak_files[c(4,1,5,2,6,3)]
  
# read peaks
peak_list <- map(peak_files, rtracklayer::import)

# filter peaks with autosomal
peak_list <- map(peak_list, keepStandardChromosomes, 
                 species = species, pruning.mode = "coarse")

# change seqlevel style
seqstyle <- seqlevelsStyle(peak_list[[1]])
for(i in seq_along(peak_list)) seqlevelsStyle(peak_list[[i]]) <- "UCSC" # NCBI
seqlevels(peak_list[[1]])

# get merged peaks
peak_merge <- GRangesList(peak_list) %>% 
  unlist() %>% 
  GenomicRanges::reduce()

# sample peaks overlap
for(i in names(peak_list)) {
  mcols(peak_merge)[[i]] <- overlapsAny(peak_merge, peak_list[[i]])
}

# convert to original style 
seqlevelsStyle(peak_merge) <- "NCBI" # seqlevels(peak_merge)

# subset smaples
if(F) {
  select <- c(5, 6)
  peak_merge <- peak_merge[rowSums(as.data.frame(mcols(peak_merge)[, select])) > 0, select]
  # export bed file
  export(peak_merge, "peak_merge.bed")
}

# save to bed for chiprofile use
export.bed(peak_merge, "peak_merge.bed")

# convert back to UCSC
seqlevelsStyle(peak_merge) <- "UCSC" # seqlevels(peak_merge)


# bam files
bam_files <- list.files(pattern = ".bam$", recursive = T, full.names = T)
bam_files <- bam_files[c(4,1,5,2,6,3)]
# # bigwig files
# bw_files <- list.files(pattern = ".bw$", recursive = T, full.names = T)
# bw_files <- bw_files[c(4,1,5,2,6,3)]

# coverage to chiprofile
# here, BamBigwig_to_chipProfile can't take named signalFiles as input, 
# otherwise, error will appear
names(bam_files) <- NULL
chipProfile <- BamBigwig_to_chipProfile(bam_files, # bw_files 
                                        "peak_merge.bed",
                                        format = "bam",
                                        paired = PE_reads, normalize = "RPM",
                                        quant_params = SnowParam(parallel::detectCores()-4) # for win
                                        # quant_params = MulticoreParam(round(parallel::detectCores()/2)
                                        )

# save in case of corruption
file_name <- "chipProfile_bam.rds"
if(!file.exists(file_name)) {
  write_rds(chipProfile, file_name)
  # chipProfile <- read_rds(file_name)
}

# convert to profileplyr
profile <- as_profileplyr(chipProfile)

# convert seqstyle back ti UCSC
seqlevelsStyle(rowRanges(profile)) <- "UCSC" # rowRanges(profile)

# rename sample
rownames(sampleData(profile)) <- rownames(sampleData(profile)) %>%
  map_chr(~ str_split(.x, "[.]", simplify = T)[1])

# annotate regions
profile <- annotateRanges(profile, 
                          TxDb = txdb_for_annotateRanges,
                          tssRegion = c(-3000, 3000))
# merge exon and intron
profile@rowRanges$Type <- 
  profile@rowRanges$annotation_short %>%
  str_replace_all("(Exon)|(Intron)|(5p UTR)|(3p UTR)|(Downstream)", "Gene Body") %>%
  factor(levels = c("Promoter", "Gene Body", "Distal Intergenic"))
table(profile@rowRanges$Type)

# summarize profile for clustering
profile_rowmeans <- assays(profile) %>%
  as.list() %>%
  map_dfc(~ rowMeans(.x))

# remove low valuerows
low_rows <- (rowSums(profile_rowmeans) < 1)
if(sum(low_rows) > 0) {
  print(str_glue("{sum(low_rows)} rows have low values [rowMeans < 1]"))
  profile_rowmeans <- profile_rowmeans[!low_rows, ]
  profile <- profile[!low_rows, ]
}


# heatmap by annotation
profile@params$rowGroupsInUse <- "Type"
pdf(file="profile_by_annotation.pdf", width = 9, heigh = 7, useDingbats = T)
generateEnrichedHeatmap(
  profile,
  group_anno_color = brewer.pal(3, "Dark2"),
  matrices_color = as.list(rep(list(c("white","red2")), nrow(profile@sampleData))),
  all_color_scales_equal = T,
  # matrices_axis_name = c("-3000", "TSS", "3000"),
  # extra_annotation_columns = c("cluster"),
  # extra_anno_width = c(6),
  # extra_anno_color = as.list(list(brewer.pal(9, "Paired")[seq_len(length(unique(profile1@rowRanges$cluster)))])),
  # genes_to_label = c("SMAD2","PSMG2","MAPK4"),
  gene_label_font_size = 6,
  use_raster = T,
  raster_quality = 3,
  return_ht_list = F)
dev.off()


# subset H3
profile_h3 <- profile[,,1:2]
tmp <- assays(profile_h3) %>%
  as.list() %>%
  map_dfc(~ rowMeans(.x))
profile_h3 <- profile_h3[rowSums(tmp) > 1, ]
pdf(file="profile_H3_by_annotation.pdf", width = 9, heigh = 7, useDingbats = T)
generateEnrichedHeatmap(
  profile_h3,
  group_anno_color = brewer.pal(3, "Dark2"),
  matrices_color = as.list(rep(list(c("white","#e41a1c")), nrow(profile_h3@sampleData))),
  all_color_scales_equal = T,
  # matrices_axis_name = c("-3000", "TSS", "3000"),
  # extra_annotation_columns = c("cluster"),
  # extra_anno_width = c(6),
  # extra_anno_color = as.list(list(brewer.pal(9, "Paired")[seq_len(length(unique(profile1@rowRanges$cluster)))])),
  # genes_to_label = c("SMAD2","PSMG2","MAPK4"),
  gene_label_font_size = 6,
  use_raster = T,
  raster_quality = 3,
  return_ht_list = F)
dev.off()

# subset TEAD
profile_tead <- profile[,,3:4]
tmp <- assays(profile_tead) %>%
  as.list() %>%
  map_dfc(~ rowMeans(.x))
profile_tead <- profile_tead[rowSums(tmp) > 1, ]
pdf(file="profile_TEAD_by_annotation.pdf", width = 9, heigh = 7, useDingbats = T)
generateEnrichedHeatmap(
  profile_tead,
  group_anno_color = brewer.pal(3, "Dark2"),
  matrices_color = as.list(rep(list(c("white","#377eb8")), nrow(profile_tead@sampleData))),
  all_color_scales_equal = T,
  # matrices_axis_name = c("-3000", "TSS", "3000"),
  # extra_annotation_columns = c("cluster"),
  # extra_anno_width = c(6),
  # extra_anno_color = as.list(list(brewer.pal(9, "Paired")[seq_len(length(unique(profile1@rowRanges$cluster)))])),
  # genes_to_label = c("SMAD2","PSMG2","MAPK4"),
  gene_label_font_size = 6,
  use_raster = T,
  raster_quality = 3,
  return_ht_list = F)
dev.off()

# subset YAP
profile_yap <- profile[,,5:6]
tmp <- assays(profile_yap) %>%
  as.list() %>%
  map_dfc(~ rowMeans(.x))
profile_yap <- profile_yap[rowSums(tmp) > 1, ]
pdf(file="profile_YAP_by_annotation.pdf", width = 9, heigh = 7, useDingbats = T)
generateEnrichedHeatmap(
  profile_yap,
  group_anno_color = brewer.pal(3, "Dark2"),
  matrices_color = as.list(rep(list(c("white","#4daf4a")), nrow(profile_yap@sampleData))),
  all_color_scales_equal = T,
  # matrices_axis_name = c("-3000", "TSS", "3000"),
  # extra_annotation_columns = c("cluster"),
  # extra_anno_width = c(6),
  # extra_anno_color = as.list(list(brewer.pal(9, "Paired")[seq_len(length(unique(profile1@rowRanges$cluster)))])),
  # genes_to_label = c("SMAD2","PSMG2","MAPK4"),
  gene_label_font_size = 6,
  use_raster = T,
  raster_quality = 3,
  return_ht_list = F)
dev.off()


# close all dev
while(dev.cur() != 1) dev.off()

`

# subset promoters
gc(reset = T)
profile_promoter <- profile[profile@rowRanges$Type == "Promoter",]

optimal_cluster <- 5
profile_promoter_rowmeans <- assays(profile_promoter) %>%
  as.list() %>%
  map_dfc(~ rowMeans(.x))
# get optimal cluster number
if(is_null(optimal_cluster)) {
  # kmeans
  nbcluster <- NbClust(profile_promoter_rowmeans, 
                       distance = "euclidean", method = "kmeans", 
                       min.nc = 4, max.nc = 12, index = "ch")
  optimal_cluster <- nbcluster$Best.nc[1]
  print(str_glue("NBcluster optimal best number: {optimal_cluster}"))
  # add cluster id to profiler
  mcols(profile_promoter)$cluster <- nbcluster$Best.partition # rowRanges(profile_promoter)
} else {
  # using indicated clusters
  cluster <- kmeans(profile_promoter_rowmeans, centers = optimal_cluster)
  # add cluster id to profiler
  mcols(profile_promoter)$cluster <- as.factor(cluster$cluster)
}

profile_promoter@params$rowGroupsInUse <- "sgGroup" # "annotation_short"
pdf(file="profile_promoter_by_cluster.pdf", width = 9, heigh = 7, useDingbats = T)
generateEnrichedHeatmap(
  profile_promoter,
  include_group_annotation = F,
  # group_anno_color = brewer.pal(optimal_cluster, "Dark2"),
  # matrices_color = as.list(rep(list(c("white","red2")), nrow(profile_promoter@sampleData))), # list(c("white","blur2"), c("white","red2")),
  matrices_color = list(c("white","#e41a1c"), c("white","#e41a1c"), 
                        c("white","#377eb8"), c("white","#377eb8"),
                        c("white","#4daf4a"), c("white","#4daf4a")),
  all_color_scales_equal = T,
  use_raster = T,
  raster_quality = 3,
  raster_device = "png",
  return_ht_list = F)
dev.off()

pdf(file="profile_promoter_by_cluster_no_equal.pdf", width = 9, heigh = 7, useDingbats = T)
generateEnrichedHeatmap(
  profile_promoter,
  include_group_annotation = F,
  # group_anno_color = brewer.pal(optimal_cluster, "Dark2"),
  matrices_color = list(c("white","#e41a1c"), c("white","#e41a1c"), 
                        c("white","#377eb8"), c("white","#377eb8"),
                        c("white","#4daf4a"), c("white","#4daf4a")),
  all_color_scales_equal = F,
  use_raster = T,
  raster_quality = 3,
  raster_device = "png",
  return_ht_list = F)
dev.off()




# subset promoters H3
gc(reset = T)
profile_promoter_h3 <- profile_promoter[,,1:2]
profile_promoter_h3 <- profile_promoter_h3[profile_promoter_h3 %over% peak_merge[rowSums(as.data.frame(mcols(peak_merge[,1:2]))) > 0, ],]
tmp <- assays(profile_promoter_h3) %>%
  as.list() %>%
  map_dfc(~ rowMeans(.x))
profile_promoter_h3 <- profile_promoter_h3[rowSums(tmp) > 1, ]
cluster <- kmeans(tmp[rowSums(tmp) > 1, ], centers = 1)
mcols(profile_promoter_h3)$cluster <- as.factor(cluster$cluster)
profile_promoter_h3@params$rowGroupsInUse <- "sgGroup" # "annotation_short"
pdf(file="profile_promoter_H3_by_cluster.pdf", width = 9, heigh = 7, useDingbats = T)
generateEnrichedHeatmap(
  profile_promoter_h3,
  include_group_annotation = F,
  # group_anno_color = brewer.pal(3, "Dark2")[1:2],
  matrices_color = as.list(rep(list(c("white","#e41a1c")), nrow(profile_promoter_h3@sampleData))), # list(c("white","blur2"), c("white","red2")),
  all_color_scales_equal = T,
  use_raster = T,
  raster_quality = 3,
  raster_device = "png",
  return_ht_list = F)
dev.off()

# subset promoters TEAD
gc(reset = T)
profile_promoter_tead <- profile_promoter[,,3:4]
profile_promoter_tead <- profile_promoter_tead[profile_promoter_tead %over% peak_merge[rowSums(as.data.frame(mcols(peak_merge[,3:4]))) > 0, ],]
tmp <- assays(profile_promoter_tead) %>%
  as.list() %>%
  map_dfc(~ rowMeans(.x))
profile_promoter_tead <- profile_promoter_tead[rowSums(tmp) > 1, ]
cluster <- kmeans(tmp[rowSums(tmp) > 1, ], centers = 1)
mcols(profile_promoter_tead)$cluster <- as.factor(cluster$cluster)
profile_promoter_tead@params$rowGroupsInUse <- "sgGroup" # "annotation_short"
pdf(file="profile_promoter_TEAD_by_cluster.pdf", width = 9, heigh = 7, useDingbats = T)
generateEnrichedHeatmap(
  profile_promoter_tead,
  include_group_annotation = F,
  # group_anno_color = brewer.pal(3, "Dark2")[1:2],
  matrices_color = as.list(rep(list(c("white","#377eb8")), nrow(profile_promoter_tead@sampleData))), # list(c("white","blur2"), c("white","red2")),
  all_color_scales_equal = T,
  use_raster = T,
  raster_quality = 3,
  raster_device = "png",
  return_ht_list = F)
dev.off()

# subset promoters YAP
gc(reset = T)
profile_promoter_yap <- profile_promoter[,,5:6]
profile_promoter_yap <- profile_promoter_yap[profile_promoter_yap %over% peak_merge[rowSums(as.data.frame(mcols(peak_merge[,5:6]))) > 0, ],]
tmp <- assays(profile_promoter_yap) %>%
  as.list() %>%
  map_dfc(~ rowMeans(.x))
profile_promoter_yap <- profile_promoter_yap[rowSums(tmp) > 1, ]
cluster <- kmeans(tmp[rowSums(tmp) > 1, ], centers = 2)
mcols(profile_promoter_yap)$cluster <- as.factor(cluster$cluster)
profile_promoter_yap@params$rowGroupsInUse <- "sgGroup" # cluster
pdf(file="profile_promoter_YAP_by_cluster.pdf", width = 9, heigh = 7, useDingbats = T)
generateEnrichedHeatmap(
  profile_promoter_yap,
  include_group_annotation = F,
  # group_anno_color = brewer.pal(3, "Dark2")[1:2],
  matrices_color = as.list(rep(list(c("white","#4daf4a")), nrow(profile_promoter_yap@sampleData))), # list(c("white","blur2"), c("white","red2")),
  all_color_scales_equal = T,
  use_raster = T,
  raster_quality = 3,
  raster_device = "png",
  return_ht_list = F)
dev.off()



# close all dev
while(dev.cur() != 1) dev.off()`






# subset promoter peaks
# close all dev


















# get counts from bam
# convert seqstyle back
seqlevelsStyle(peak_merge) <- seqstyle[1]
peak_count <- summarizeOverlaps(peak_merge, # GenomicAlignments
                                BamFileList(bam_files, asMates = PE_reads), 
                                singleEnd = !PE_reads, 
                                fragments = F)
as_tibble(assay(peak_count))
as_tibble(assay(peak_count)) %>% colSums()

# bam counts (seq depth)
bam_count <- map_dfr(bam_files,
                     ~ Rsamtools::idxstatsBam(.x) %>% 
                       as_tibble(), 
                     .id = "sample") %>%
  group_by(sample) %>%
  summarise(counts = sum(mapped)) %>%
  deframe()


# # edgeR norm
# peak_count_norm <- edgeR::calcNormFactors(peak_count, method = "none")
# peak_count_norm_cpm <- cpm(peak_count_norm)
# 
# as_tibble(assay(peak_count))
# colSums(assay(peak_count))

((assay(peak_count) / colSums(assay(peak_count)) ) * 1e6) %>% as_tibble() %>% colSums()

peak_count_norm_cpm <- as_tibble(assay(peak_count))
# use bam count to normalize peak count
for(i in names(bam_count)) {
  peak_count_norm_cpm[, i] <- peak_count_norm_cpm[, i] / (bam_count[i] / 1e6)
}
peak_count_norm_cpm %>% colSums()

# add nrom cpm to peak_merge
peak_merge_norm <- granges(peak_count)
seqlevelsStyle(peak_merge_norm) <- "UCSC"
mcols(peak_merge_norm) <- peak_count_norm_cpm # assay(peak_count)

# norm_factor <- peak_count_norm$samples %>%
#   rownames_to_column("sample") %>%
#   transmute(sample,
#             norm_factor = lib.size * norm.factors) %>%
#   deframe()

gene <- genes(txdb)
tss <- promoters(gene, 0, 1)
tss <- keepStandardChromosomes(tss, "Homo_sapiens", "coarse")

find_overlaps(tss, peak_merge_norm) %>% mcols() %>% .[, -1] %>% as_tibble() %>% colSums()

tss <- tss[100:1000,]

rowsums <- mcols(peak_merge_norm) %>% as.data.frame() %>% rowSums()
tss <- peak_merge_norm[rowsums > 10]

mat_list <- list()
for(i in colnames(mcols(peak_merge_norm))) {
  mat_list <- list.append(mat_list, 
    normalizeToMatrix(peak_merge_norm, tss, value_column = i, 
                      extend = 3000, mean_mode = "w0", w = 10, 
                      smooth = T, background = 0)
    )
}
names(mat_list) <- colnames(mcols(peak_merge_norm))

plist <- list()
for(i in colnames(mcols(peak_merge_norm))) {
  p <- EnrichedHeatmap(mat_list[[i]], name = i, column_title = i,
                       # quantile(peak_merge_norm$ET825_H3, c(0.02, 0.98)),
                       top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = seq_len(3)))), 
                       colorRamp2(c(0, 100), c("white", "red")),
                       axis_name_rot = 90)
  plist <- list.append(plist, p)
}

plots <- plist[[4]] + plist[[1]] + plist[[5]] + plist[[2]] + plist[[6]] + plist[[3]]



partition = paste0("cluster", kmeans(mat_list$ET825_Y, centers = 3)$cluster)
lgd = Legend(at = c("cluster1", "cluster2", "cluster3"), title = "Clusters", 
             type = "lines", legend_gp = gpar(col = 2:4))
cluster_anno <- Heatmap(partition, name = "partition",
                        col = structure(2:4, names = paste0("cluster", 1:3)), 
                        show_row_names = FALSE, width = unit(3, "mm"))
plots_anno <- cluster_anno + plots
draw(plots_anno, split = partition, annotation_legend_list = list(lgd), 
     # ht_gap = unit(c(2, 8, 8), "mm")
     )



draw(cluster_anno + plist[[4]] + plist[[1]], 
     split = partition, annotation_legend_list = list(lgd), 
     ht_gap = unit(c(2, 8, 8), "mm"))




mat1 = normalizeToMatrix(peak_merge_norm, tss, value_column = "V_Y", 
                         extend = 3000, mean_mode = "w0", w = 50, 
                         smooth = T, background = 0)
mat2 = normalizeToMatrix(peak_merge_norm, tss, value_column = "ET825_Y", 
                         extend = 3000, mean_mode = "w0", w = 50, 
                         smooth = T, background = 0)
mat3 = normalizeToMatrix(peak_merge_norm, tss, value_column = "V_Y", 
                         extend = 3000, mean_mode = "w0", w = 50, 
                         smooth = T, background = 0)
mat4 = normalizeToMatrix(peak_merge_norm, tss, value_column = "ET825_Y", 
                         extend = 3000, mean_mode = "w0", w = 50, 
                         smooth = T, background = 0)
mat5 = normalizeToMatrix(peak_merge_norm, tss, value_column = "V_Y", 
                         extend = 3000, mean_mode = "w0", w = 50, 
                         smooth = T, background = 0)
mat6 = normalizeToMatrix(peak_merge_norm, tss, value_column = "ET825_Y", 
                         extend = 3000, mean_mode = "w0", w = 50, 
                         smooth = T, background = 0)

EnrichedHeatmap(mat1, name = "V_Y",
                colorRamp2(quantile(mat1, c(0, 0.95), na.rm = T), c("white", "red")),
                axis_name_rot = 90) +
EnrichedHeatmap(mat2, name = "ET825_Y",
                colorRamp2(quantile(mat1, c(0, 0.95), na.rm = T), c("white", "red")),
                axis_name_rot = 90)



## for single rep only
# get de peaks
peak_de <- peak_count_norm_cpm %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(H3_cpm = mean(c(ET825_H3, V_H3)), 
         H3_log2FC = log2((ET825_H3 + 1)/(V_H3+1)),
         TEAD_cpm = mean(c(ET825_T, V_T)), 
         TEAD_log2FC = log2((ET825_T + 1)/(V_T+1)),
         YAP_cpm = mean(c(ET825_Y, V_Y)), 
         YAP_log2FC = log2((ET825_Y + 1)/(V_Y+1))
         ) %>%
  ungroup()

# annotate and tidy
cpm <- select(peak_de, c(contains("_cpm"))) %>%
  mutate(id = row_number()) %>%
  pivot_longer(!id, names_to = c("marker", ".value"), names_sep = "_")
log2FC <- select(peak_de, c(contains("_log2FC"))) %>%
  mutate(id = row_number()) %>%
  pivot_longer(!id, names_to = c("marker", ".value"), names_sep = "_")
peak_de_tidy <- left_join(cpm, log2FC) %>%
  mutate(direction = case_when(
    log2FC >= 2 & cpm >= 5 ~ "up",
    log2FC <= -2 & cpm >= 5 ~ "down",
    log2FC >= 2 & cpm > 1 & cpm < 5 ~ "up_low",
    log2FC <= -2 & cpm > 1 & cpm < 5 ~ "down_low",
    TRUE ~ "others"))
# count
peak_de_tidy %>% dplyr::count(marker, direction)

# density plot
peak_de_tidy %>%
  ggplot(aes(log10(cpm), color = marker)) +
  geom_density(size = 1) +
  scale_x_log10() +
  theme_minimal()
# MA plot
peak_de_tidy %>%
  ggplot(aes(log10(cpm), log2FC)) +
  geom_pointdensity() +
  lemon::facet_rep_wrap(~ marker, repeat.tick.labels = T) +
  theme_minimal()



# use normalized score
norm_factor <- peak_count_norm$samples %>%
  rownames_to_column("sample") %>%
  transmute(sample,
            norm_factor = lib.size * norm.factors) %>%
  deframe()

# add cpm
for(i in names(peak_list)) {
  peak_list[[i]]$cpm <- as.integer((peak_list[[i]]$score / 
                                      norm_factor[which(names(norm_factor)==i)]) *1e6) # for cpm, multuply 1e6
}

# for(i in names(peak_list)) {
#   peak_list[[i]] <- find_overlaps(peak_list[[i]], peak_merge_norm_split[[i]])
# }


# overall coverage
tagMatrix <- map(peak_list[c("V_Y","ET825_Y")],
                 ~ getTagMatrix(.x, # windows=promoter
                                weightCol = "cpm",
                                type = "start_site", by = "gene", #nbin = 10,
                                upstream = 3000, downstream = 3000,
                                TxDb = txdb))
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red2")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), facet = "none", free_y = F) # for TSS
# plotPeakProf(tagMatrix, conf = F, facet = "none", free_y = F) # for any region

# de peaks only
peak_de_index <- which(peak_de$YAP_log2FC < -2 & # abs(peak_de$YAP_log2FC)
                         peak_de$YAP_cpm > 5)
peak_merge_de <- peak_merge_norm[peak_de_index, ]
peak_de_select <- list(
  V_Y = find_overlaps(peak_list[["V_Y"]], peak_merge_de),
  ET825_Y = find_overlaps(peak_list[["ET825_Y"]], peak_merge_de))

tagMatrix <- map(peak_de_select,
                 ~ getTagMatrix(.x, # windows=promoter
                                weightCol = "cpm",
                                type = "start_site", by = "gene", #nbin = 10,
                                upstream = 3000, downstream = 3000,
                                TxDb = txdb))
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red2")
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), facet = "row", free_y = F)


tagMatrix <- map(peak_list[c("V_Y","ET825_Y")],
                 ~ getTagMatrix(.x, # windows=promoter
                                weightCol = "cpm",
                                type = "start_site", by = "gene", #nbin = 10,
                                upstream = 3000, downstream = 3000,
                                TxDb = txdb))
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = F, facet = "row", free_y = F)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red2")
# plotPeakProf(tagMatrix, conf = F, facet = "row", free_y = F)





# # prepare data
# peak_de_select <- list(
#   V_Y = peak_merge_norm[,"V_Y"],
#   ET825_Y = peak_merge_norm[,"ET825_Y"])
# colnames(mcols(peak_de_select$V_Y)) <- "cpm"
# colnames(mcols(peak_de_select$ET825_Y)) <- "cpm"

# tagMatrix <- getTagMatrix(peak_de_select, # windows=promoter
#                           weightCol = "cpm",
#                           type = "start_site", by = "gene", #nbin = 10,
#                           upstream = 3000, downstream = 3000,
#                           TxDb = txdb)
# tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red2")


# prof plot
tagMatrix <- map(peak_de_select,
                 ~ getTagMatrix(.x, # windows=promoter
                                weightCol = "cpm",
                                type = "start_site", by = "gene", #nbin = 10,
                                upstream = 3000, downstream = 3000,
                                TxDb = txdb))
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), facet = "none", free_y = TRUE)
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red2")



# de peaks coverage
comparisons <- list(c("ET825_H3","V_H3"), c("ET825_T","V_T"), c("ET825_Y","V_Y"))







## using MAnorm
library(MAnorm2)
browseVignettes("MAnorm2")

# prepare occupancy
for(p in names(narrowpeaks)) {
  # print(p)
  hits <- findOverlaps(merged_peaks, narrowpeaks[[p]]) %>% 
    as_tibble()
  mcols(merged_peaks) <- mutate(as_tibble(mcols(merged_peaks)), 
                                "{p}.occupancy" := ifelse(row_number() %in% hits$queryHits, 1, 0))
}

# get bam files
bam_files <- list.files(pattern = ".bam$", recursive = T, full.names = T)
names(bam_files) <- map_chr(bam_files, 
                             ~ str_split(.x, "BAM/", simplify = T) %>%
                               .[2] %>%
                               str_split(pattern = "[.]", simplify = T) %>%
                               .[1])
PE_reads <- TRUE
# to bamlist
bam_files <- BamFileList(bam_files, asMates = PE_reads)

# get counts from bam
counted_peaks <- GenomicAlignments::summarizeOverlaps(merged_peaks, 
                                                      bam_files, 
                                                      singleEnd = !PE_reads, fragments = F)
# counts <- counted_peaks
counted_peaks <- bind_cols(as_tibble(granges(counted_peaks)),
                           as_tibble(assay(counted_peaks))) %>%
  as_granges()
mcols(counted_peaks) %>% as_tibble()


comparisons <- list(c("ET825_H3","V_H3"), c("ET825_T","V_T"), c("ET825_Y","V_Y"))

# Perform within-group normalization
counted_peaks_norm <- as_tibble(counted_peaks)
contrast_list <- list()
for(c in comparisons) {
  counted_peaks_norm <- counted_peaks_norm %>% 
    normalize(count = c[1], occupancy = str_glue("{c[1]}.occupancy")) %>%
    normalize(count = c[2], occupancy = str_glue("{c[2]}.occupancy"))
  
  contrast_list[[str_glue("{c[1]}_vs_{c[2]}")]] <- 
    list(bioCond(counted_peaks_norm[c[1]], counted_peaks_norm[str_glue("{c[1]}.occupancy")], name = c[1]),
         bioCond(counted_peaks_norm[c[2]], counted_peaks_norm[str_glue("{c[2]}.occupancy")], name = c[2])) %>%
    set_names(c(c[1], c[2]))
}
contrast_list



contrast_list_fit <- list()
contrast_list_res <- list()
for(c in names(contrast_list)) {
  print(c)
  
  # Perform between-group normalization
  # contrast_list_norm[[c]] <- normBioCond(contrast_list[[c]]) # for multiple samples # multiple reps
  cns <- names(contrast_list[[c]])
  contrast_list[[c]]$blind <- 
    bioCond(counted_peaks_norm[c(cns[1], cns[2])], 
            counted_peaks_norm[c(str_glue("{cns[1]}.occupancy"), str_glue("{cns[2]}.occupancy"))], 
            occupy.num = 2, name = "blind") # for single rep
  # summary
  summary(contrast_list_fit[[c]]$blind)
  
  # Fit a mean-variance curve
  contrast_list_fit[[c]] <- 
    fitMeanVarCurve(contrast_list[[c]], method = "parametric", init.coef = c(0.1, 10))
  
  # Perform differential tests
  contrast_list_res[[c]] <- diffTest(contrast_list_fit[[c]][[1]], contrast_list_fit[[c]][[2]])
}

MAplot(contrast_list_res[[c]], pval = 0.01)
abline(h = 0, lwd = 2, lty = 5, col = "green3")

contrast_list_res[[1]] %>%
  as_tibble() %>%
  filter(padj < 0.05) %>%
  group_by(sign(Mval.t)) %>%
  summarise(n = n())








## annotate peaks
annoated_peaks <- granges(counted_peaks)
seqlevelsStyle(merged_peaks)
seqlevelsStyle(annoated_peaks) <- "UCSC"
seqlevels(annoated_peaks)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annoated_peaks <- ChIPseeker::annotatePeak(annoated_peaks, 
                                           TxDb = txdb, 
                                           level = "gene")
annoated_peaks <- as_granges(annoated_peaks)
seqlevelsStyle(annoated_peaks) <- "NCBI"

# merge with counts
data <- left_join(as_tibble(annoated_peaks), counted_peaks)

# columns(org.Hs.eg.db)
anno_symbol <- bitr(unique(data$geneId), 
                    "ENTREZID", "SYMBOL", org.Hs.eg.db) %>%
  as_tibble() %>%
  rename(geneId = ENTREZID, symbol = SYMBOL) %>%
  distinct(geneId, .keep_all = T)

# add gene symbol
data <- left_join(data, anno_symbol) %>%
  relocate(symbol, .after = distanceToTSS)

# write counts matrix
writexl::write_xlsx(data, str_glue("{Sys.Date()}_counts.xlsx"))






## diff test
library(edgeR)

# to see distribution
data %>%
  mutate(tyoe = ifelse(str_detect(annotation, "Promoter"), "Promoter",
                       ifelse(str_detect(annotation, "Intergenic"), "Intergenic", "GeneBody"))) %>%
  group_by(tyoe) %>%
  summarise(across(ET825_T:last_col(), sum))
  
min_counts <- 10
# summarize to gene promoters
data_promoter <- data %>%
  filter(str_detect(annotation, "Promoter")) %>%
  # filter(!str_detect(annotation, "Intergenic")) %>%
  group_by(symbol) %>%
  summarise(across(ET825_T:last_col(), sum)) %>%
  column_to_rownames("symbol") %>%
  .[rowSums(.) > min_counts, ]
dim(data_promoter)

# summarize to gene promoters + gene bodies
data_gene <- data %>%
  filter(!str_detect(annotation, "Intergenic")) %>%
  # filter(!str_detect(annotation, "Intergenic")) %>%
  group_by(symbol) %>%
  summarise(across(ET825_T:last_col(), sum)) %>%
  column_to_rownames("symbol") %>%
  .[rowSums(.) > min_counts, ]
dim(data_gene)


# construct deglist
deglist <- data_promoter %>%
  DGEList(., group = colnames(.))






library(DiffBind)

# get peak files
peak_files <- list.files(pattern = ".narrowPeak", recursive = T, full.names = T)
names(peak_files) <- map_chr(peak_files, ~ str_split(.x, "/", simplify = T)[4])
# get bam files
bam_files <- list.files(pattern = ".bam$", recursive = T, full.names = T)
names(bam_files) <- map_chr(bam_files, 
                            ~ str_split(.x, "BAM/", simplify = T) %>%
                              .[2] %>%
                              str_split(pattern = "[.]", simplify = T) %>%
                              .[1])


samples <- enframe(peak_files, "SampleID", "Peaks") %>%
  mutate(PeakCaller = "narrow") %>%
  left_join(enframe(bam_files, "SampleID", "bamReads")) %>%
  mutate(Condition = rep(c("H3", "TEAD", "YAP"), times = 2),
         Treatment = rep(c("ET825", "Veh"), each = 3)) %>%
  group_by(Condition, Treatment) %>%
  mutate(Replicate = row_number())

# construct diffbind obj
df <- dba(sampleSheet = samples)

# correlation
plot(df)

# bam counts
df <- dba.count(df, bParallel = T)


# normalize
df <- dba.normalize(df)
# df <- dba.normalize(df, method = DBA_EDGER, 
#                     normalize = DBA_NORM_LIB, library=DBA_LIBSIZE_FULL)

# contrast
df <- dba.contrast(df)

# run analysis
df <- dba.analyze(df)

# correlation using de peaks
plot(df, contrast=1)

# get de peaks
df_de <- dba.report(df)



## QC plot
# venn
dba.plotVenn(df, contrast=1, bDB=TRUE,
             bGain=TRUE, bLoss=TRUE, bAll=FALSE)

# pca
dba.plotPCA(df, attributes=c(DBA_CONDITION,DBA_TREATMENT), label = DBA_TREATMENT)

# pca using de peaks
dba.plotPCA(df, DBA_CONDITION, label = DBA_TREATMENT, contrast=1)

# ma
dba.plotMA(df, bXY = F)

# volcano
dba.plotVolcano(df)

# heatmap
dba.plotHeatmap(df, correlations=F,
                # contrast=1, 
                scale="row", 
                score=DBA_SCORE_NORMALIZED, 
                colScheme = colorRampPalette(c("#034e7b", "white", "#cb181d"))(n = 100))



## peak plot
profiles <- dba.plotProfile(df)
dba.plotProfile(profiles)
dba.plotProfile(profiles, merge=c(DBA_TISSUE, DBA_REPLICATE))





