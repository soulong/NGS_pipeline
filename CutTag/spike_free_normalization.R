
# Spike-free normalization (ChIPseqSpikeInFree)
#
# Reads CutTag/config.yml and the dataset samplesheet, then estimates
# spike-free scaling factors from the aligned BAMs in <root_dir>/result/03_bam.
# Output: <prefix>_SF.txt (and QC plots) in result/03_bam, which run_cuttag.sh
# consumes via the `spike_free_sf_file` config key when norm_method: SpikeFree.
#
# Usage: run from RStudio (open this file and source it) or
#        Rscript spike_free_normalization.R
# Requires: the CutTag pipeline has completed alignment (03_bam/*.filtered.bam).

# devtools::install_github("stjude/ChIPseqSpikeInFree")
library(ChIPseqSpikeInFree)
library(tidyverse)
library(yaml)

# Locate this script (RStudio active document, Rscript --file=, or cwd)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  doc <- rstudioapi::getActiveDocument()$path
  script_dir <- if (!is.null(doc) && nzchar(doc)) dirname(doc) else getwd()
} else {
  ca <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", ca, value = TRUE)
  script_dir <- if (length(fa) > 0) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
}
setwd(script_dir)

config <- read_yaml("config.yml")

# Sample sheet and BAM directory derive from the pipeline's own root_dir,
# NOT from a hardcoded relative path (the two used to drift apart)
root_dir    <- config$root_dir
samplesheet <- file.path(root_dir, config$samplesheet)
bam_dir     <- file.path(root_dir, config$outdir, "03_bam")

# Chrom sizes for the species that was actually aligned against
# (was hardcoded to human chm13, which breaks every non-human dataset)
fai_file <- switch(config$species,
  chm13 = file.path(config$index_rootdir, "hs/chm13/chm13v2.0.fa.gz.fai"),
  hs    = file.path(config$index_rootdir, "hs/v49/GRCh38.primary_assembly.genome.fa.gz.fai"),
  mm    = file.path(config$index_rootdir, "mm/vM38/GRCm39.primary_assembly.genome.fa.gz.fai"),
  fly   = file.path(config$index_rootdir, "fly/Drosophila_melanogaster.BDGP6.54.dna.toplevel.fa.gz.fai"),
  cel   = file.path(config$index_rootdir, "cel/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz.fai"),
  stop("Unknown species: ", config$species)
)
if (!file.exists(fai_file)) stop("Chrom sizes file not found: ", fai_file)

sample_info <- read_csv(samplesheet) %>%
  mutate(ID = str_c(sample, ".filtered.bam"))

# ChIPseqSpikeInFree groups samples by ANTIBODY: use the samplesheet 'target'
# column (a hardcoded antibody would pool all targets into one model)
if (any(is.na(sample_info$target) | !nzchar(sample_info$target))) {
  stop("Empty 'target' in samplesheet; fill it before running spike-free normalization")
}

# Drop samples whose filtered BAM is missing (pipeline not finished / rerun later)
missing_bam <- !file.exists(file.path(bam_dir, sample_info$ID))
if (any(missing_bam)) {
  warning("Missing filtered BAMs (excluded): ",
          paste(sample_info$ID[missing_bam], collapse = ", "))
  sample_info <- sample_info[!missing_bam, ]
}
if (nrow(sample_info) == 0) stop("No usable filtered BAMs in ", bam_dir)

save_prefix <- "SpikeFree"
# Absolute path: the meta file is written here but read AFTER switching to bam_dir
meta_file <- normalizePath(file.path(script_dir, str_c(save_prefix, "_metaFile.txt")),
                           mustWork = FALSE)

tibble(ID = sample_info$ID, ANTIBODY = sample_info$target, GROUP = sample_info$group) %>%
  write_tsv(meta_file)

# ChIPseqSpikeInFree resolves bamFiles (and writes output) relative to the
# working directory, so run it from the BAM directory
setwd(bam_dir)

ChIPseqSpikeInFree(bamFiles = sample_info$ID,
                   chromFile = fai_file,
                   metaFile = meta_file,
                   prefix = save_prefix,
                   cutoff_QC = 1.2,
                   maxLastTurn = 0.99,
                   ncores = 6)


# # transform sf for bigwig/bedgraph normalization
# statistics <- "../statistics_CPM.csv" %>%
#   read_csv() %>% #glimpse()
#   # transmute(sample, filtered_reads) %>%
#   left_join(sample_info) %>%
#   print()
#
# scale_factor <- str_c(save_prefix, "_SF.txt") %>%
#   read_table() %>%
#   left_join(statistics) %>%
#   mutate(scale_factor=1000000 / (filtered_reads * SF)) %>%
#   glimpse()
#
# # save
# scale_factor %>%
#   write_csv(str_c(save_prefix, ".csv"))
#
# scale_factor %>%
#   dplyr::select(sample, scale_factor=scale_factor, SF) %>%
#   write_csv(str_c(save_prefix, "_scale_factors.csv"))
