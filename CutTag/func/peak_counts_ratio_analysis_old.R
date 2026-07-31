
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

setwd("D:\\EED_CutTag\\fly\\result")

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

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


# ==== Genomic resources (minimal, standard chromosomes only) ====
genes <- genes(txdb) %>%
  keepStandardChromosomes(pruning.mode="coarse") %>%
  sort()
seqlevelsStyle(genes) <- "UCSC"
seqlevels(genes)

genes_reduced <- genes %>%
  resize(width=width(.) + 2000, fix="center") %>%
  GenomicRanges::reduce()

tss_reduced <- promoters(genes, 1000, 1000) %>%
  GenomicRanges::reduce()

consensus <- rtracklayer::import(consensus_peakfile) %>%
  keepStandardChromosomes(pruning.mode="coarse")
seqlevelsStyle(consensus) <- "UCSC"
seqlevels(consensus)

# ==== Functions (optimized for parallelization) ====
compute_region <- function(gr, g, t) {
  in_gene <- GenomicRanges::reduce(subsetByOverlaps(gr, g, ignore.strand=FALSE))
  flank_gene <- GenomicRanges::reduce(c(
    flank(in_gene, width=10000, start=TRUE,  both=FALSE),
    flank(in_gene, width=10000, start=FALSE, both=FALSE)
  ))
  
  in_tss <- GenomicRanges::reduce(subsetByOverlaps(gr, t, ignore.strand=FALSE))
  flank_tss <- GenomicRanges::reduce(c(
    flank(in_tss, width=3000, start=TRUE,  both=FALSE),
    flank(in_tss, width=3000, start=FALSE, both=FALSE)
  ))
  
  list(
    in_gene=in_gene, flank_gene=flank_gene,
    in_tss=in_tss, flank_tss=flank_tss
  )
}

do_sampling <- function(peak, peak_name, g, t, sample_size) {
  if (length(peak) < sample_size) sample_size <- length(peak)
  sampled <- sample(peak, size=sample_size, replace=FALSE)
  regions <- compute_region(sampled, g, t)
  names(regions) <- paste0(peak_name, "_", names(regions))
  regions[[peak_name]] <- sampled
  regions
}

# Main worker function: minimal input, no side effects
get_stats_one_rep <- function(
    bam_gr, peak,
    g, t, consens,
    sample_size
) {
  if (sample_size==0) {
    gene_sampled <- g
    tss_sampled  <- t
  } else {
    gene_sampled <- sample(g, sample_size, replace=FALSE)
    tss_sampled  <- sample(t, sample_size, replace=FALSE)
  }

  
  gene_flank <- GenomicRanges::reduce(c(
    flank(gene_sampled, 10000, start=TRUE,  both=FALSE),
    flank(gene_sampled, 10000, start=FALSE, both=FALSE)
  ))
  tss_flank <- GenomicRanges::reduce(c(
    flank(tss_sampled, 3000, start=TRUE,  both=FALSE),
    flank(tss_sampled, 3000, start=FALSE, both=FALSE)
  ))
  
  regions <- list(
    gene=gene_sampled,
    tss=tss_sampled,
    gene_flank=gene_flank,
    tss_flank=tss_flank
  )
  
  if (!is.null(peak))      regions <- c(regions, do_sampling(peak, "peak", g, t, sample_size))
  if (!is.null(consens))   regions <- c(regions, do_sampling(consens, "consensus", g, t, sample_size))
  
  counts <- map_int(regions, ~ sum(countOverlaps(bam_gr, .x, ignore.strand=TRUE)))
  enframe(counts, name="region", value="count") %>%
    pivot_wider(names_from=region, values_from=count)
}

# Wrapper: read BAM once
get_stats <- function(bam_path, peak_path,
                      g, t, consens,
                      sample_rep=200, sample_size=3000) {
  
  # Read BAM ONCE per sample (outside future_map)
  bam_gr <- GenomicAlignments::readGAlignmentPairs(bam_path) %>% granges()
  seqlevelsStyle(bam_gr) <- "UCSC"

  
  # Clip sample_size if too large
  sample_size <- min(sample_size, length(bam_gr), length(g), length(t))
  
  # peak
  peak_gr <- rtracklayer::import(peak_path) %>% 
    keepStandardChromosomes(pruning.mode="coarse") %>% 
    GenomicRanges::reduce()
  seqlevelsStyle(peak_gr) <- "UCSC"
  
  res <- map(
    seq_len(sample_rep),
    \(x) get_stats_one_rep(bam_gr, peak_gr, g, t, consens, sample_size))
  
  data <- list_rbind(res, names_to="sampling_rep") %>%
    mutate(total=length(bam_gr), .before=1)
  
  return(data)
}


# ==== Run in parallel across samples ====
# Flatten: one future per (sample, rep) is too fine; better: one future per sample
plan(multisession, workers=min(10, nrow(df)))  # one worker per sample

stats <- future_map2(
  df$bam_file, df$peak_file,
  \(x, y) get_stats(x, y, genes_reduced, tss_reduced, consensus,
                    sample_rep=100, sample_size=3000)) %>%
  set_names(df$sample)

plan(sequential)


# ==== Post-processing ====
stats_tidy <- stats %>%
  list_rbind(names_to="sample") %>%
  mutate(
    peak_of_total=peak / total,
    consensus_of_total=consensus / total,
    gene_of_total=gene / total,
    tss_of_total=tss / total,
    gene_in_of_flank=gene / gene_flank,
    tss_in_of_flank=tss / tss_flank,
    consensus_gene_in_of_flank=consensus_in_gene / consensus_flank_gene,
    consensus_tss_in_of_flank=consensus_in_tss / consensus_flank_tss,
    peak_in_of_flank=peak_in_gene / peak_flank_gene,
    peak_tss_in_of_flank=peak_in_tss / peak_flank_tss
  )

# Save
writexl::write_xlsx(stats_tidy, glue::glue("{Sys.Date()}_stats_tidy_old.xlsx"))

# Summary
summary_stats <- stats_tidy %>%
  summarise(across(everything(), mean), .by=sample) %>% 
  glimpse()


# Plot
metrics <- colnames(stats_tidy) %>% str_subset("_of_")
plots <- map(metrics, ~ {
  stats_tidy %>%
    separate_wider_delim(sample, delim="-", names=c("group","rep"), cols_remove=F) %>% 
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
  ggsave(glue::glue("{Sys.Date()}_stats_tidy_old.pdf"), ., 
         width=8, height=3 * length(plots))




# stats <- future_map2(
#   df$bam_file, df$peak_file,
#   \(x, y) get_stats(x, y, genes_reduced, tss_reduced, consensus,
#                     sample_rep=1, sample_size=0)) %>%
#   set_names(df$sample)



