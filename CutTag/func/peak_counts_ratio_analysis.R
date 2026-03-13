
# ==== Setup ====
rstudioapi::getActiveDocumentContext()$path |>
  dirname() |>
  setwd()

library(tidyverse)
library(GenomicAlignments)
library(GenomicRanges)
library(GenomeInfoDb)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)  # Only mouse needed
library(rtracklayer)
library(furrr)

options(future.globals.maxSize=2 * 1024^3)  # 2GB if needed

setwd("F:/workspace/IMR/result")
# setwd("F:/workspace/J009_IMR_0_0.3_0.6/result")


# ==== Input ====
peak_files <- list.files("04_peaks", pattern="broadPeak$", full.names=TRUE, recursive=FALSE)
df <- tibble(
  peak_file=peak_files,
  sample=basename(peak_files) %>% str_replace("_peaks.*Peak$", ""),
  bam_file=str_c("03_alignment/", sample, ".filtered.bam")
)

# ==== Genomic resources (minimal, standard chromosomes only) ====
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene) %>%
# genes <- genes(TxDb.Mmusculus.UCSC.mm39.knownGene) %>%
  keepStandardChromosomes(pruning.mode="coarse") %>%
  sort()

genes_reduced <- genes %>%
  resize(width=width(.) + 2000, fix="center") %>%
  reduce()

tss_reduced <- promoters(genes, 1000, 1000) %>%
  reduce()

consensus <- rtracklayer::import("04_peaks/consensus.bed") %>%
  keepStandardChromosomes(pruning.mode="coarse")


# ==== Functions (optimized for parallelization) ====
compute_region <- function(gr, g, t) {
  in_gene <- reduce(subsetByOverlaps(gr, g, ignore.strand=FALSE))
  flank_gene <- reduce(c(
    flank(in_gene, width=10000, start=TRUE,  both=FALSE),
    flank(in_gene, width=10000, start=FALSE, both=FALSE)
  ))
  
  in_tss <- reduce(subsetByOverlaps(gr, t, ignore.strand=FALSE))
  flank_tss <- reduce(c(
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

  
  gene_flank <- reduce(c(
    flank(gene_sampled, 10000, start=TRUE,  both=FALSE),
    flank(gene_sampled, 10000, start=FALSE, both=FALSE)
  ))
  tss_flank <- reduce(c(
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
  
  # Clip sample_size if too large
  sample_size <- min(sample_size, length(bam_gr), length(g), length(t))
  
  # peak
  peak_gr <- rtracklayer::import(peak_path) %>% 
    keepStandardChromosomes(pruning.mode="coarse") %>% 
    reduce()
  
  res <- map(
    seq_len(sample_rep),
    \(x) get_stats_one_rep(bam_gr, peak_gr, g, t, consens, sample_size))
  
  data <- list_rbind(res, names_to="sampling_rep") %>%
    mutate(total=length(bam_gr), .before=1)
  
  return(data)
}


# ==== Run in parallel across samples ====
# Flatten: one future per (sample, rep) is too fine; better: one future per sample
plan(multisession, workers=min(6, nrow(df)))  # one worker per sample

stats <- future_map2(
  df$bam_file, df$peak_file,
  \(x, y) get_stats(x, y, genes_reduced, tss_reduced, consensus,
                    sample_rep=200, sample_size=3000)) %>%
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
writexl::write_xlsx(stats_tidy, glue::glue("{Sys.Date()}_stats_tidy.xlsx"))

# Summary
summary_stats <- stats_tidy %>%
  summarise(across(everything(), mean), .by=sample) %>% 
  glimpse()


# Plot
metrics <- colnames(stats_tidy) %>% str_subset("_of_")
plots <- map(metrics, ~ {
  stats_tidy %>%
    separate_wider_delim(sample, delim="_", names=c("group","replicagte"), cols_remove=F) %>% 
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





stats <- future_map2(
  df$bam_file, df$peak_file,
  \(x, y) get_stats(x, y, genes_reduced, tss_reduced, consensus,
                    sample_rep=1, sample_size=0)) %>%
  set_names(df$sample)

# plan(sequential)


