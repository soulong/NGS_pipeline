library(tidyverse)
library(GenomicRanges)
library(ChIPseeker)
library(GenomeInfoDb)
library(EnrichedHeatmap)
ht_opt$message = FALSE
library(circlize)
library(DESeq2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(org.Hs.eg.db)

# load helper from parent directory
parent_path <- rstudioapi::getActiveDocumentContext()$path |> 
  dirname() |> dirname()
source(file.path(parent_path, 'ngs_helper.R'))


# ============================= Configuration =============================
data_dir <- norm_path(r"(E:\NGS\ATACseq_dataset\result)") %>% 
  setwd()

annodb <- "org.Hs.eg.db"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

set_style <- function(x) { 
  seqlevelsStyle(x) <- "UCSC"
  x <- keepStandardChromosomes(x, pruning.mode="tidy")
  return(x)
}

# Geneset of interest (example: chromatin regulators)
geneset_interest <- c(
  "EZH2","EED","SUZ12","JARID2","AEBP2","PCL1","PCL2","PCL3",
  "RING1A","RING1B","CBX2","CBX4","CBX6","CBX7","CBX8",
  "KDM1A","KDM2B","KDM4A","KDM4B","KDM4C","KDM5A","KDM5B",
  "HDAC1","HDAC2","HDAC3","SIN3A","SIN3B","NCOR1","NCOR2",
  "BRD2","BRD3","BRD4","BRDT",
  "CTCF","RAD21","SMC1A","SMC3","STAG1","STAG2",
  "YY1","TEAD1","TEAD2","TEAD3","TEAD4"
)


# ============================= Import BigWig Files =============================
# run used in heatmap profiling
bw_suffix <- '.CPM.bw'
bw_file <- '03_bam' %>% 
  list.files(pattern=bw_suffix, full.names=T) %>%
  set_names(., nm=map_chr(., \(x) basename(x) %>% 
                            str_replace(fixed(bw_suffix),''))) %>% 
  # subset samples (example: select specific groups)
  # .[which(str_detect(names(.), 'WT|KO'))] %>% 
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


# ============================= Import Peak Files =============================
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


# ============================= Consensus Peaks =============================
## Compute consensus peaks from replicates
consensus <- peak %>% 
  # subset to specific target if needed
  # .[which(str_detect(names(.), 'WT'))] %>% 
  compute_consensus(., min_occurrence=2, min_gapwidth=1) %>% 
  {.$all} %>% 
  # remove small fragments
  {.[width(.) > 20]} %>% 
  print()

# Annotate consensus peaks
anno <- consensus %>% 
  ChIPseeker::annotatePeak(tssRegion=c(-3000, 3000),
                           TxDb=txdb,
                           level='gene',
                           annoDb=annodb) %>% 
  {.@anno}

# Filter peaks near TSS of genes of interest
consensus_filtered <- anno %>% 
  {.[.$geneId %in% geneset_interest]} %>% 
  {.[abs(.$distanceToTSS) <= 3000]} %>% 
  print()


# ============================= Count Reads in Peaks =============================
## Get count matrix over consensus peaks
bam_suffix <- '.filtered.bam'
bam_file <- '03_bam' %>% 
  list.files(pattern=bam_suffix, full.names=T) %>%
  str_subset(".bai", negate=T) %>% 
  set_names(., nm=map_chr(., \(x) basename(x) %>% 
                            str_replace(fixed(bam_suffix),''))) %>% 
  # subset samples
  # .[which(str_detect(names(.), 'WT|KO'))] %>% 
  print()

# Get counts using chromVAR (recommended for ATACseq)
count <- chromVAR::getCounts(
  bam_file, format = "bam", consensus,
  paired = TRUE, by_rg = FALSE)
count <- assay(count, 'counts')

# Set rownames
rownames(count) <- as_tibble(consensus) %>% 
  unite('uid', 1:3) %>% pull(uid)
colnames(count) <- names(bam_file)
count[1:2,1:2]


# ============================= Differential Accessibility Analysis =============================
## Create sample metadata
col_data <- colnames(count) %>% 
  tibble(sample=.) %>% 
  mutate(
    group = sapply(strsplit(sample, "_"), `[`, 1),
    condition = case_when(
      str_detect(group, "WT") ~ "Control",
      str_detect(group, "KO") ~ "Knockout",
      str_detect(group, "Treat") ~ "Treatment",
      TRUE ~ group
    )
  ) %>% 
  column_to_rownames("sample") %>% 
  print()

## Run DESeq2
dds <- DESeqDataSetFromMatrix(
  count, col_data, design = ~ condition) %>% 
  DESeq()

## Extract results
dds_result <- results(dds, contrast=c("condition", "Knockout", "Control")) %>% 
  as_tibble(rownames='uid') %>% 
  separate(uid, c("seqnames","start","end","width"), sep="_") %>% 
  arrange(pvalue) %>% 
  print()

## Annotate DARs (Differential Accessible Regions)
dds_result_gr_anno <- dds_result %>% 
  makeGRangesFromDataFrame(keep.extra.columns=T) %>%
  ChIPseeker::annotatePeak(tssRegion=c(-3000, 3000),
                           TxDb=txdb,
                           level='gene',
                           annoDb=annodb) %>% 
  { .@anno } %>% 
  print()

## Tidy results for visualization
deg_tidy <- as_tibble(dds_result_gr_anno) %>% 
  mutate(
    geneset_interest = ifelse(geneId %in% geneset_interest, geneId, NA),
    type = case_when(
      log2FoldChange > log2(1.5) & padj < 0.05 ~ 'gain',
      log2FoldChange < log2(1/1.5) & padj < 0.05 ~ 'loss',
      TRUE ~ NA_character_
    ),
    label = ifelse(!is.na(geneset_interest) & !is.na(type), geneset_interest, NA)
  ) %>% 
  glimpse()

dplyr::count(deg_tidy, type)

## Save results
writexl::write_xlsx(
  deg_tidy, str_glue('{Sys.Date()}_ATACseq_DARs.xlsx'))


# ============================= Volcano Plot =============================
p <- deg_tidy %>% 
  ggplot(aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(color=type), alpha=0.7, show.legend=F) +
  scale_color_manual(values=c(gain="#e41a1c", loss="#377eb8"), na.value="gray80") +
  ggrepel::geom_text_repel(aes(label=label), size=5, show.legend=F) +
  theme_bw() +
  labs(x="Log2 Fold Change", y="-Log10 P-value", title="Differential Accessibility")
p <- ggrastr::rasterize(p, dpi=600)
ggsave(str_glue('{Sys.Date()}_ATACseq_volcano.pdf'),
       p, width=6, height=6)


# ============================= MA Plot =============================
p_ma <- deg_tidy %>% 
  ggplot(aes(baseMean, log2FoldChange)) +
  geom_point(aes(color=type), alpha=0.7, show.legend=F) +
  scale_color_manual(values=c(gain="#e41a1c", loss="#377eb8"), na.value="gray80") +
  ggrepel::geom_text_repel(aes(label=label), size=5, show.legend=F) +
  theme_bw() +
  labs(x="Mean Normalized Counts", y="Log2 Fold Change", title="MA Plot")
p_ma <- ggrastr::rasterize(p_ma, dpi=600)
ggsave(str_glue('{Sys.Date()}_ATACseq_MA.pdf'),
       p_ma, width=6, height=6)


# ============================= Heatmap Visualization =============================
# Signal of interest
soi <- bw %>% .[which(str_detect(names(.), 'WT|KO|Treat'))]

# Regions of interest
rois <- deg_tidy %>%
  dplyr::filter(!is.na(type)) %>%
  split(., .[['type']]) %>%
  map(\(x) makeGRangesFromDataFrame(x, keep.extra.columns=T)) %>%
  {c(consensus=consensus_filtered, .)} %>%
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
  
  # set rowname
  consensus_rowname <- roi@elementMetadata %>% as_tibble() %>% 
    mutate(uid=str_glue('{geneId}_{row_number()}'),.by=geneId) %>% 
    pull(uid)
  
  # label genes of interest
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
  f_name <- names(rois)[idx]
  pdf(str_glue("{Sys.Date()}_heatmap_{f_name}.pdf"), 
      width = length(result$heatmap)*1.6, height = 8)
  print(result$heatmap)
  while(dev.cur() != 1) dev.off()
}


# ============================= Nucleosome Positioning Analysis =============================
# Fragment length distribution from BAM QC
frag_file <- 'multiqc/BAM_fragment_length.txt'

if(file.exists(frag_file)) {
  frag_data <- read_tsv(frag_file, skip=1) %>%
    pivot_longer(-1, names_to="sample", values_to="count") %>%
    mutate(
      fragment_length = as.numeric(`# fragment length`),
      category = case_when(
        fragment_length < 100 ~ "Nucleosome-free",
        fragment_length < 250 ~ "Mononucleosome",
        fragment_length < 450 ~ "Dinucleosome",
        TRUE ~ "Trinucleosome+"
      )
    )
  
  p_frag <- frag_data %>%
    ggplot(aes(fragment_length, count, fill=category)) +
    geom_area(alpha=0.7, position="identity") +
    facet_wrap(~sample, scales="free_y", ncol=2) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(x="Fragment Length (bp)", y="Read Count", title="Nucleosome Positioning")
  
  ggsave(str_glue('{Sys.Date()}_nucleosome_positioning.pdf'),
         p_frag, width=10, height=8)
}


# ============================= TSS Enrichment Score =============================
# Calculate TSS enrichment for each sample
tss <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)

tss_enrichment <- list()
for(i in seq_along(bw)) {
  sample_name <- names(bw)[i]
  signal <- bw[[i]]
  
  # Calculate signal around TSS
  tss_signal <- signalOverRegions(signal, tss, bin_size=100)
  
  # TSS enrichment = signal at TSS / signal at flanking regions
  tss_center <- rowMeans(tss_signal[, 45:55])
  tss_flank <- rowMeans(cbind(rowMeans(tss_signal[, 1:10]), 
                               rowMeans(tss_signal[, 90:100])))
  
  enrichment <- mean(tss_center / tss_flank, na.rm=TRUE)
  tss_enrichment[[sample_name]] <- enrichment
}

tss_enrich_df <- tibble(
  sample = names(tss_enrichment),
  TSS_enrichment = unlist(tss_enrichment)
) %>% 
  arrange(desc(TSS_enrichment))

print(tss_enrich_df)

# Save TSS enrichment
writexl::write_xlsx(tss_enrich_df, str_glue('{Sys.Date()}_TSS_enrichment.xlsx'))

# Plot TSS enrichment
p_tss <- tss_enrich_df %>%
  ggplot(aes(reorder(sample, TSS_enrichment), TSS_enrichment)) +
  geom_col(fill="#377eb8") +
  coord_flip() +
  theme_bw() +
  labs(x="Sample", y="TSS Enrichment Score", title="ATACseq Quality Metric") +
  geom_hline(yintercept=5, linetype="dashed", color="red") +
  annotate("text", x=Inf, y=5, label="Threshold (5x)", hjust=1.1, vjust=-0.5, color="red")

ggsave(str_glue('{Sys.Date()}_TSS_enrichment.pdf'),
       p_tss, width=8, height=6)
