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

# locate this script (RStudio active document, Rscript --file=, or cwd);
# ngs_helper.R lives in the parent directory
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  doc <- rstudioapi::getActiveDocumentContext()$path
  script_dir <- if (!is.null(doc) && nzchar(doc)) dirname(doc) else getwd()
} else {
  ca <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", ca, value = TRUE)
  script_dir <- if (length(fa) > 0) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
}
source(file.path(dirname(script_dir), 'ngs_helper.R'))


# ============================= Configuration =============================
# Dataset result directory (produced by run_atacseq.sh).
# Override via the ATACSEQ_RESULT_DIR environment variable if needed.
data_dir <- Sys.getenv("ATACSEQ_RESULT_DIR",
                       unset = norm_path(r"(E:\NGS\ATACseq_dataset\result)"))
setwd(data_dir)

# Annotation resources MUST match the genome the reads were aligned against.
# NOTE: peaks aligned to chm13 (T2T) coordinates must NOT be annotated with
# this hg38 TxDb (Mb-level coordinate differences) — align to hs/GRCh38 or
# wire up a chm13-compatible TxDb here.
species <- "hg38"
anno_pkg <- switch(species,
  hg38 = list(annodb = "org.Hs.eg.db",
              txdb   = TxDb.Hsapiens.UCSC.hg38.knownGene),
  stop("No annotation database wired up for species: ", species)
)
require(anno_pkg$annodb, character.only = TRUE)
annodb <- anno_pkg$annodb
txdb   <- anno_pkg$txdb

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
# run used in heatmap profiling — pick up bigwigs from any normalization
# the pipeline produced (config norm_method: CPM/Spike/SpikeFree)
bw_file <- '03_bam' %>%
  list.files(pattern='\\.(CPM|Spike|SpikeFree)\\.bw$', full.names=T) %>%
  set_names(., nm=map_chr(., \(x) basename(x) %>%
                            str_remove('\\.(CPM|Spike|SpikeFree)\\.bw$'))) %>%
  # subset samples (example: select specific groups)
  # .[which(str_detect(names(.), 'WT|KO'))] %>%
  print()

# reorder as heatmap needs — natural sort on the 2nd "_"-token so numeric
# series sort 2 < 10 (plain alphabetical order would give "10" < "2")
group_token <- sapply(strsplit(names(bw_file), "_"), `[`, 2)
new_order <- order(match(group_token, str_sort(unique(group_token), numeric = TRUE)))
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

# Filter peaks near TSS of genes of interest.
# NOTE: ChIPseeker's `geneId` column holds ENTREZ ids; gene SYMBOLs live in
# the `SYMBOL` column. Matching symbols against geneId always yields empty.
consensus_filtered <- anno %>%
  {.[.$SYMBOL %in% geneset_interest]} %>%
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

## Extract results.
## uid was built from 3 fields (seqnames_start_end) — separating it into 4
## names produced an all-NA `width` column.
contrast_pair <- c("Knockout", "Control")   # adjust to your design
dds_result <- results(dds, contrast=c("condition", contrast_pair[1], contrast_pair[2])) %>%
  as_tibble(rownames='uid') %>%
  separate(uid, c("seqnames","start","end"), sep="_") %>%
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

## Tidy results for visualization (SYMBOL, not the Entrez geneId, is matched
## against the geneset symbols)
deg_tidy <- as_tibble(dds_result_gr_anno) %>%
  mutate(
    geneset_hit = ifelse(SYMBOL %in% geneset_interest, SYMBOL, NA),
    type = case_when(
      log2FoldChange > log2(1.5) & padj < 0.05 ~ 'gain',
      log2FoldChange < log2(1/1.5) & padj < 0.05 ~ 'loss',
      TRUE ~ NA_character_
    ),
    label = ifelse(!is.na(geneset_hit) & !is.na(type), geneset_hit, NA)
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

  # set rowname (SYMBOL + occurrence index so duplicate genes stay unique)
  consensus_rowname <- roi@elementMetadata %>% as_tibble() %>%
    mutate(uid=str_glue('{SYMBOL}_{row_number()}'),.by=SYMBOL) %>%
    pull(uid)

  # label genes of interest — exact prefix match on the uid (SYMBOL_n);
  # the old regex alternation (paste0(..., collapse="_|")) could never match
  genes_in_rois <- str_replace(consensus_rowname, "_[0-9]+$", "")
  labels <- consensus_rowname[genes_in_rois %in% geneset_interest]

  # plot
  result <- heatmap_profile(
    mat_list,
    color_scales = c(0, NA),
    colors = c("white", "#b12923"),
    rownames=consensus_rowname,
    labels=labels
  )

  # save (close only the device opened above, not every open device)
  f_name <- names(rois)[idx]
  pdf(str_glue("{Sys.Date()}_heatmap_{f_name}.pdf"),
      width = length(result$heatmap)*1.6, height = 8)
  dev_id <- dev.cur()
  print(result$heatmap)
  dev.off(dev_id)
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
# TSS enrichment = signal at the TSS vs. the flanks of the +/-2kb promoter
# window (standard ATACseq QC metric; ENCODE threshold ~ >5).
# The promoters come from getPromoters(upstream=2000, downstream=2000), i.e.
# the window ALREADY spans TSS-2000..TSS+2000, so normalizeToMatrix must use
# extend=c(0,0) (adding another 2kb on each side would double the window).
# Bin count depends on how normalizeToMatrix settles the bin width, so the
# center/flank indices are derived from ncol() as fractions of the window:
#   center = the middle bins straddling the TSS, flanks = outermost 10%/side.
tss <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)

bin_size <- 100

tss_enrichment <- list()
for(i in seq_along(bw)) {
  sample_name <- names(bw)[i]
  signal <- bw[[i]]

  # signal matrix over the promoter windows (one row per TSS)
  tss_signal <- EnrichedHeatmap::normalizeToMatrix(
    signal, tss,
    value_column = "score",
    extend = c(0, 0),
    mean_mode = "w0",
    w = bin_size
  )

  n_bin <- ncol(tss_signal)
  half  <- n_bin %/% 2
  center_ix   <- unique(round(half + c(-1.5, -0.5, 0.5, 1.5)))
  n_flank     <- max(1L, round(n_bin * 0.1))
  flank_l     <- 1:n_flank
  flank_r     <- (n_bin - n_flank + 1):n_bin

  # TSS enrichment = signal at TSS / signal at flanking regions
  tss_center <- rowMeans(tss_signal[, center_ix])
  tss_flank  <- rowMeans(cbind(rowMeans(tss_signal[, flank_l]),
                               rowMeans(tss_signal[, flank_r])))

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
