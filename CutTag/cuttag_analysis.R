
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
txdb_file <- 'F:/index/hs/chm13/t2t_chm13v2_txdb.sqlite'
if(file.exists(txdb_file)) {
  txdb <- AnnotationDbi::loadDb(txdb_file)
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
  AnnotationDbi::saveDb(txdb, txdb_file)
}

# geneset
geneset_interest <- c(
  "ADM","AKAP2","ALCAM","AMOTL2","ANKRD1","ANKRD13A","ANLN","ARHGAP29",
  "ARHGEF28","AXL","C15orf52","CCDC80","COL12A1","CPA4","CRIM1","CTGF",
  "CTNNA1","CYR61","CCN1","CCN2","DKK1","DUSP1","EFEMP1","EPB41L2","EVA1A","F3",
  "FOSL1","FSP1","GADD45A","GADD45B","GAS6","IGFBP3","IRS1","LEPREL1",
  "LIMA1","MGLL","NABP1","NEDD4","OPHN1","PALM2-AKAP2","PDLIM2","PLK2",
  "ROR1","SEMA3C","SMURF2","STARD13","STOM","TAGLN","TBC1D23","TEAD1",
  "TGM2","THBS1","UBASH3B","WWC2"
)


# bigwig ------------------------------------------------
# run used in heatmap profiling
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



# peak ------------------------------------------------
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


## consensus ------------------------------------------------

# subset sample peak
consensus <- peak %>% 
  # .[which(str_detect(names(.), 'TEAD'))] %>% 
  # , group=names(.) %>% str_replace_all("_\\d","")
  compute_consensus(., min_occurrence=2, min_gapwidth=1) %>% 
  {.$all} %>% 
  # remove small fragment
  {.[width(.) > 20]}

# annotate peak
anno <- consensus %>% 
  ChIPseeker::annotatePeak(tssRegion=c(-3000, 3000),
                           TxDb=txdb,
                           level='gene',
                           annoDb=annodb) %>% 
  {.@anno}

# filter
consensus <- anno %>% 
  {.[.$geneId %in% geneset_interest]} %>%
  {.[abs(.$distanceToTSS) <= 3000]} %>%
  print()



## subset overlap ------------------------------------------------
# tss <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
# 
# hit <- findOverlaps(tss, consensus, minoverlap=1000)
# tss_overlap_with_consensus <- tss[unique(queryHits(hit))]
# 
# hit <- findOverlaps(consensus, tss, minoverlap=1000)
# consensus_overlap_with_tss <- consensus[unique(queryHits(hit))]



# bam ------------------------------------------------
# run only wnat to get DESeq2 degs using bam count
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


## deseq2 ------------------------------------------------
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
deg_tidy <- as_tibble(dds_result_gr_anno) %>% 
  mutate(geneset_interest=ifelse(
    geneId %in% geneset_interest, geneId, NA)) %>% 
  mutate(type=ifelse(
    log2FoldChange > log2(1.5) & padj < 0.05, 'pos',
    ifelse(log2FoldChange < log2(1/1.5) & padj < 0.05, 'neg', NA))) %>% 
  mutate(label=ifelse(
    !is.na(geneset_interest) & type %in% 'pos', geneset_interest, NA)) %>% 
  glimpse()
dplyr::count(deg_tidy, type)
# save
writexl::write_xlsx(
  deg_tidy, str_glue('{Sys.Date()}_valcano_consensus_tead_deseq2_deg.xlsx'))

# valcono
p <- deg_tidy %>% 
  ggplot(aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(color=type), alpha=0.7, show.legend=F) +
  scale_color_discrete(direction = -1) +
  ggrepel::geom_text_repel(aes(label=label), size=5, show.legend=F) +
  theme_bw()
p <- ggrastr::rasterize(p, dpi=600)
ggsave(str_glue('{Sys.Date()}_valcano_consensus_tead_deseq2_deg.pdf'),
       p, width=5, height=5)

# p <- deg_tidy %>% 
#   ggplot(aes(base_mean, log2fold_change)) +
#   geom_point(aes(color=type), alpha=0.7, show.legend=F) +
#   scale_color_discrete(direction = -1) +
#   ggrepel::geom_text_repel(aes(label=label), size=5, show.legend=F) +
#   theme_bw()
# p <- ggrastr::rasterize(p, dpi=600)
# ggsave(str_glue('{Sys.Date()}_ma_consensus_tead_deseq2_deg.pdf'),
#        p, width=5, height=5)




# heatmap ------------------------------------------------

# signal of interest
soi <- bw %>% .[which(str_detect(names(.), 'TEAD|BRD4'))]

# region of interest
# rois <- dplyr::filter(deg_tidy, !is.na(type)) %>%
#   split(., .[['type']]) %>%
#   map(\(x) makeGRangesFromDataFrame(x, keep.extra.columns=T)) %>%
#   print()
rois <- GRangesList(yap_target=consensus) %>% print()
  
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
  # set rowname
  consensus_rowname <- roi@elementMetadata %>% as_tibble() %>% 
    mutate(uid=str_glue('{geneId}_{row_number()}'),.by=geneId) %>% 
    pull(uid)
  labels <- str_subset(
    consensus_rowname, paste0(geneset_interest, collapse="_|"))
  # plot
  result <- heatmap_profile(
    mat_list,
    color_scales = c(0, NA), 
    colors = c("white", "#b12923"),
    rownames=consensus_rowname,
    labels=labels
  )
  # save
  # f_name <- str_glue('consensus_tead_and_deg_{names(rois)[idx]}')
  f_name <- str_glue('{names(rois)[idx]}')
  pdf(str_glue("{Sys.Date()}_heatmap_{f_name}.pdf"), 
      width = length(result$heatmap)*1.6, height = 8)
  print(result$heatmap)
  while(dev.cur() != 1) dev.off()
}






