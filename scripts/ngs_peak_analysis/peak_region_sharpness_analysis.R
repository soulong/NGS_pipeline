# set path
# rstudioapi::getActiveDocumentContext()$path |>
#   dirname() |>
#   setwd()
setwd("/mnt/f/workspace/MEF/result/")

# load packages
library(tidyverse)
library(GenomicAlignments)
library(GenomicRanges)
library(GenomeInfoDb)
library(furrr)
# library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
options(future.globals.maxSize = 2 * 1024 * 1024^2) # 2 GB


## define input -------------------------------
peak_files <- list.files("04_peaks", pattern = "broadPeak$", full.names = T, recursive = F)
print(peak_files)
df <- tibble(
  peak_file = peak_files,
  sample = map_chr(peak_files, basename) %>%
    str_replace_all("_peaks.*Peak$", ""),
  bam_file = str_c("03_alignment/", sample, ".filtered.bam"),
) %>% print()


## consensus peak -------------------------------
consensus <- rtracklayer::import("04_peaks/consensus.bed")
# consensus_annno <- ChIPseeker::annotatePeak(
#   consensus,
#   tssRegion = c(-3000, 3000),
#   TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#   level = "gene")
# consensus <- consensus_annno@anno
consensus <- keepStandardChromosomes(
  consensus,
  pruning.mode = "coarse"
)


## gene/tss -------------------------------
genes <- 
  genes(TxDb.Mmusculus.UCSC.mm39.knownGene) %>% 
# genes(TxDb.Hsapiens.UCSC.hg38.knownGene) %>%
  keepStandardChromosomes(pruning.mode = "coarse") %>%
  sort() %>%
  print()

# used in following functions
genes_reduced <- genes %>%
  resize(width = width(genes) + 2000, fix = "center") %>%
  reduce() # %>% print()
tss_reduced <- promoters(genes, 1000, 1000) %>%
  reduce() # %>% print()


compute_region <- function(gr, g = genes_reduced, t = tss_reduced) {
  in_gene <- subsetByOverlaps(gr, g, ignore.strand = FALSE) %>%
    reduce() # %>% print()
  flank_gene <- c(
    flank(in_gene, width = 10000, start = T, both = F),
    flank(in_gene, width = 10000, start = F, both = F)
  ) %>% reduce() # %>% print()

  in_tss <- subsetByOverlaps(gr, t, ignore.strand = FALSE) %>%
    reduce() # %>% print()
  flank_tss <- c(
    flank(in_tss, width = 3000, start = T, both = F),
    flank(in_tss, width = 3000, start = F, both = F)
  ) %>% reduce() # %>% print()

  return(list(in_gene, flank_gene, in_tss, flank_tss) %>%
    set_names(nm = c("in_gene", "flank_gene", "in_tss", "flank_tss")))
}

do_sampling <- function(peak, peak_name = "peak", sample_size = 1000) {
  peak_sampled <- sample(consensus, size = sample_size, replace = F)
  peak_regions <- compute_region(peak_sampled)
  names(peak_regions) <- str_c(peak_name, "_", names(peak_regions))
  peak_regions[[peak_name]] <- peak_sampled
  return(peak_regions)
}
# do_sampling(consensus, "consensus")

get_stats_summary <- function(bam_path,
                              peak_path,
                              g = genes_reduced,
                              t = tss_reduced,
                              consens = consensus,
                              sample_rep = NULL,
                              sample_size = 1000 # if <= 1, then run all regions, only used when sample_rep != NULL
) {
  # bam_path <- df$bam_file[1]
  bam_gr <- GenomicAlignments::readGAlignmentPairs(bam_path) %>%
    granges()

  if (sample_size <= 1) sample_size <- floor(length(bam_gr) * sample_size)
  if (sample_size > length(bam_gr)) sample_size <- length(sample_size)

  # peak_path <- df$peak_file[1]
  peak <- rtracklayer::import(peak_path) %>% reduce()

  # set rep
  if (is.null(sample_rep)) {
    sample_rep <- 1
    sample_size <- length(bam_gr)
  }

  # run here
  res <- list()
  for (i in seq_len(sample_rep)) {
    gene_sampled <- sample(g, sample_size, replace = F)
    tss_sampled <- sample(t, sample_size, replace = F)
    regions <- c(
      gene = list(gene_sampled),
      tss = list(tss_sampled),
      gene_flank = c(
        flank(gene_sampled, width = 10000, start = T, both = F),
        flank(gene_sampled, width = 10000, start = F, both = F)
      ) %>% reduce() %>% list(),
      tss_flank = c(
        flank(tss_sampled, width = 3000, start = T, both = F),
        flank(tss_sampled, width = 3000, start = F, both = F)
      ) %>% reduce() %>% list()
    )
    if (!is.null(peak_path)) regions <- c(regions, do_sampling(peak, "peak", sample_size))
    if (!is.null(consens)) regions <- c(regions, do_sampling(consens, "consensus", sample_size))
    # print(regions)

    gr_counts <- map_int(
      regions, \(x) countOverlaps(bam_gr, x, ignore.strand = T) %>% sum()
    ) %>%
      enframe() %>%
      pivot_wider()

    gr_counts_list <- list(gr_counts) %>% set_names(i)

    res <- c(res, gr_counts_list)
  }

  data <- list_rbind(res, names_to = "rep") %>%
    mutate(total = length(bam_gr), .before = 1)

  return(data)
}

stats <- map2(
  df2$bam_file, df2$peak_file,
  \(x, y) get_stats_summary(x, y, genes_reduced, tss_reduced,
    consensus,
    sample_rep = NULL, sample_size = 1000, return_sum = FALSE
  )
) %>%
  set_names(nm = df2$sample)


stats_tidy <- stats %>%
  list_rbind(names_to = "sample") %>% # glimpse()
  # dplyr::rename(gene_flank=genes_flank) %>%
  mutate(
    peak_of_total = peak / total,
    consensus_of_total = consensus / total,
    gene_of_total = gene / total,
    tss_of_total = tss / total,
    gene_in_of_flank = gene / gene_flank,
    tss_in_of_flank = tss / tss_flank,
    consensus_gene_in_of_flank = consensus_in_gene / consensus_flank_gene,
    consensus_tss_in_of_flank = consensus_in_tss / consensus_flank_tss,
    peak_in_of_flank = peak_in_gene / peak_flank_gene,
    peak_tss_in_of_flank = peak_in_tss / peak_flank_tss
  ) %>%
  # view()
  glimpse()

stats_tidy %>%
  reframe(across(everything(), mean), .by = sample) %>%
  view()
writexl::write_xlsx(stats_tidy, str_glue("stats_tidy.xlsx"))




## count_bams_over_grange --------------------

count_bams_over_grange <- function(bam_files,
                                   gr,
                                   paired = TRUE,
                                   ignore.strand = TRUE,
                                   sample_names = NULL) {
  # validate inputs
  if (!is.character(bam_files) || length(bam_files) == 0) stop("bam_files must be a non-empty character vector")
  if (!is(gr, "GRanges")) stop("gr must be a GRanges object")
  # derive sample names
  if (is.null(sample_names)) sample_names <- tools::file_path_sans_ext(basename(bam_files))
  if (length(sample_names) != length(bam_files)) stop("sample_names must match length of bam_files")
  # read bam alignments (each once) and convert to GRanges
  bam_grs <- set_names(lapply(bam_files, function(bf) {
    if (paired) {
      GenomicAlignments::readGAlignmentPairs(bf) %>% GenomicAlignments::granges()
    } else {
      GenomicAlignments::readGAlignments(bf) %>% GenomicAlignments::granges()
    }
  }), sample_names)
  # compute counts for the region and attach as metadata columns
  counts_list <- lapply(bam_grs, function(bgr) {
    as.integer(GenomicRanges::countOverlaps(gr, bgr, ignore.strand = ignore.strand))
  })
  # ensure names on list -> column names in DataFrame
  counts_df <- S4Vectors::DataFrame(counts_list)
  # preserve existing metadata and bind new columns
  mcols(gr) <- S4Vectors::DataFrame(cbind(as.data.frame(mcols(gr)), as.data.frame(counts_df)))
  return(gr)
}

region_list <- list(
  consensus = consensus,
  consensus_l = flank(consensus, width = 3000, start = T, both = F),
  consensus_r = flank(consensus, width = 3000, start = F, both = F),
  tss = tss_reduced,
  tss_l = flank(tss_reduced, width = 3000, start = T, both = F),
  tss_r = flank(tss_reduced, width = 3000, start = F, both = F)
)

region_count_list <- map(
  region_list, \(x) count_bams_over_grange(df$bam_file, x, paired = T, ignore.strand = T)
)

write_rds(region_count_list, "region_count_list.rds")
# region_count_list <- read_rds("region_count_list.rds")


plan(cluster, workers = 6)
region_count <- region_count_list %>%
  future_map(
    \(x) as_tibble(x) %>%
      rename_with(\(x) str_replace(x, ".filtered", ""), everything()) %>%
      dplyr::select(!c(seqnames:strand)) %>%
      mutate(region_id = row_number(), .before = 1) %>%
      pivot_longer(!c(region_id), names_to = "sample", values_to = "count")
  ) %>%
  glimpse()

region_count_tidy <- list(
  consensus =
    region_count[1:3] %>%
      list_rbind(names_to = "region") %>%
      pivot_wider(names_from = region, values_from = count) %>%
      rename_with(\(x) str_replace(x, "(consensus)|(tss)", "peak")) %>%
      separate_wider_regex("sample", patterns = c(group="[A-Z]+\\d+", "_\\d+"), cols_remove = F) %>% 
      reframe(across(peak:peak_r, sum), .by = c(region_id, group)) %>%
          mutate(peak_ratio = round(peak / (peak_l + peak_r), 4)) %>%
      left_join(as_tibble(consensus) %>% mutate(region_id = row_number(), .before = 1), .) %>% 
      mutate(group = str_replace_all(group, "(IMR)|(MEF)", "Conc")) %>% 
      pivot_wider(names_from = group, values_from = c(peak, peak_l, peak_r, peak_ratio)) %>% 
      mutate(ratio_diff_05=peak_ratio_Conc05 - peak_ratio_Conc0,
      ratio_diff_1=peak_ratio_Conc1 - peak_ratio_Conc0) %>% 
      dplyr::arrange(dplyr::desc(ratio_diff_05)),
  tss =
    region_count[4:6] %>%
      list_rbind(names_to = "region") %>%
      pivot_wider(names_from = region, values_from = count) %>%
      rename_with(\(x) str_replace(x, "(consensus)|(tss)", "peak")) %>%
      separate_wider_regex("sample", patterns = c(group="[A-Z]+\\d+", "_\\d+"), cols_remove = F) %>% 
      reframe(across(peak:peak_r, sum), .by = c(region_id, group)) %>%
          mutate(peak_ratio = round(peak / (peak_l + peak_r), 4)) %>%
      left_join(as_tibble(tss_reduced) %>% mutate(region_id = row_number(), .before = 1), .) %>% 
      mutate(group = str_replace_all(group, "(IMR)|(MEF)", "Conc")) %>% 
      pivot_wider(names_from = group, values_from = c(peak, peak_l, peak_r, peak_ratio)) %>% 
      mutate(ratio_diff_05=peak_ratio_Conc05 - peak_ratio_Conc0,
      ratio_diff_1=peak_ratio_Conc1 - peak_ratio_Conc0) %>% 
      dplyr::arrange(dplyr::desc(ratio_diff_05))
)

writexl::write_xlsx(region_count_tidy, str_glue("{Sys.Date()}_region_count_list.xlsx"))


