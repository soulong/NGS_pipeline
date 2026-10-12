# ==============================================================================
# genome_fingerprint_cdf.R
# Redraw deepTools plotFingerprint results (BAM_fingerprint.txt) as:
#   1. the fingerprint curve per sample (replicates averaged per bin; the CDF
#      is scale-invariant, so mean and sum pooling give the same curve)
#   2. the cumulative distribution of per-bin coverage, plotted against the
#      actual coverage value (log scale), with a low-coverage zoom
#   3. a numeric summary: zero-coverage bin fraction, median coverage, and the
#      share of reads captured by the top 1 / 5 / 10 percent of bins
#
# Outputs (into out_dir, the parent of result_dir):
#   <date>_cdf.pdf                           fingerprint curve
#   <date>_cdf_cov_ecdf.pdf / _lowcov.pdf    coverage cumulative distributions,
#                                            bin-rank axis + coverage axis side by side
#   <date>_fingerprint_summary.csv           numeric summary per sample
#   <date>_cdf*.xlsx                         thinned curve / distribution tables
# ==============================================================================

library(tidyverse)

# ---- config ------------------------------------------------------------------
result_dir       <- "E:\\EED\\2026-07-21_CutTag_cel_fly_mm\\fly\\result\\multiqc"
out_dir          <- dirname(result_dir)   # plots / tables go next to the multiqc folder
fingerprint_file <- "BAM_fingerprint.txt"

# sample name pattern, adjust per project (column names after clean_names):
sample_pattern <- c(sample = "[a-zA-Z0-9]+", "[-_]", rep = "\\d+")   # e.g. d1_1
# sample_pattern <- c(sample = ".*", "[-_]", rep = "\\d+")           # e.g. MN0_1

# sample_levels <- c("d1", "d5", "l", "m", "ctrl", "nmn")  # plot/legend order; NULL = alphabetical
sample_levels <- NULL

bin_width  <- 10      # coverage bin width for the cumulative coverage plot
low_cov_max <- 10000  # coverage cutoff for the low-coverage zoom
thin_by    <- 200     # keep every n-th point when drawing curves
top_fracs  <- c(0.01, 0.05, 0.1)  # "signal in top x% of bins" summary metrics

# ---- functions ---------------------------------------------------------------

# Read deepTools plotFingerprint raw counts into long format: one row per
# replicate x genome bin
read_fingerprint <- function(path, sample_pattern) {
  read_delim(path, delim = "\t", skip = 1) |>
    janitor::clean_names() |>
    mutate(bin_id = row_number()) |>
    pivot_longer(-bin_id, names_to = "name", values_to = "coverage") |>
    separate_wider_regex(name, patterns = sample_pattern)
}

# Order samples by the configured factor levels (no-op if levels are NULL)
apply_sample_levels <- function(tbl) {
  if (is.null(sample_levels)) tbl else mutate(tbl, sample = factor(sample, sample_levels))
}

# deepTools-style fingerprint curve: bins ordered by descending coverage,
# cumulative fraction of signal vs fraction of genome bins.
# arrange() sorts globally by group column first, so the per-group cumsum
# below always runs over descending coverage within each group
fingerprint_curve <- function(coverage_tbl, .by) {
  coverage_tbl |>
    arrange(pick(all_of(.by)), desc(coverage)) |>
    mutate(
      cumfrac_reads = cumsum(coverage) / sum(coverage),
      frac_genome   = row_number() / n(),
      .by           = all_of(.by)
    ) |>
    ungroup()
}

# Cumulative distribution of per-bin coverage on fixed-width coverage bins.
# x is the actual coverage value (bin lower edge), not a rank
coverage_hist <- function(coverage_tbl, bin_width, max_cov = Inf) {
  coverage_tbl |>
    filter(coverage < max_cov) |>
    mutate(coverage_bin = floor(coverage / bin_width) * bin_width, .by = sample) |>
    summarise(n_bins = n(), .by = c(sample, coverage_bin)) |>
    arrange(sample, coverage_bin) |>
    mutate(cumfrac_bins = cumsum(n_bins) / sum(n_bins),
           bin_rank = row_number(), .by = sample) |>
    ungroup()
}

# Numeric fingerprint summary: signal concentration per sample
fingerprint_summary <- function(coverage_tbl, top_fracs = c(0.01, 0.05, 0.1)) {
  coverage_tbl |>
    arrange(sample, desc(coverage)) |>
    mutate(rank = row_number(), .by = sample) |>
    group_by(sample) |>
    group_modify(\(d, .grp) {
      total <- sum(d$coverage)
      out   <- tibble(
        n_bins        = nrow(d),
        total_reads   = total,
        zero_cov_frac = mean(d$coverage == 0),
        median_cov    = median(d$coverage)
      )
      for (f in top_fracs) {
        k <- max(1L, floor(nrow(d) * f))
        out[[sprintf("top_%g_pct_reads_frac", f * 100)]] <-
          sum(d$coverage[d$rank <= k]) / total
      }
      out
    }) |>
    ungroup()
}

# Keep every n-th row per sample (thinning dense curves for plotting / export).
# Must be per sample: a global slice drops each sample's last point (the CDF
# would never reach 1.0) and shifts the sampling phase between samples whose
# row counts are not multiples of `by`.
thin_rows <- function(tbl, by = thin_by) {
  sampled <- slice(tbl, seq(1, n(), by = by), .by = sample)
  last    <- slice(tbl, n(), .by = sample)          # keep the final CDF point
  bind_rows(sampled, last) |>
    distinct(sample, frac_genome, .keep_all = TRUE) |>
    arrange(sample, frac_genome)
}

# Shared look: line + palette + labels + theme
common_style <- function(x_lab, y_lab) {
  list(
    geom_line(linewidth = 1),
    scale_color_brewer(palette = "Paired"),
    labs(x = x_lab, y = y_lab, color = "Sample"),
    theme_bw(base_size = 8)
  )
}

# Side-by-side cumulative coverage distribution: bin-rank x axis (old-style
# fingerprint-like ordering, 1,2,3,...) next to the actual coverage value axis
cov_ecdf_figure <- function(hist_tbl, subtitle = NULL) {
  p_rank <- ggplot(hist_tbl, aes(bin_rank, cumfrac_bins, color = sample)) +
    common_style("Coverage bin rank (ascending)", "Cumulative Fraction of Bins")
  p_cov  <- ggplot(hist_tbl, aes(coverage_bin + 1, cumfrac_bins, color = sample)) +
    common_style("Coverage per bin (log scale)", "Cumulative Fraction of Bins") +
    scale_x_log10()
  patchwork::wrap_plots(p_rank, p_cov, nrow = 1) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = subtitle)
}

# ---- run ---------------------------------------------------------------------
raw <- read_fingerprint(file.path(result_dir, fingerprint_file), sample_pattern)
raw <- apply_sample_levels(raw)

# Sample-level table: replicates averaged per bin (the CDF is scale-invariant,
# so this equals pooling the reads)
sample_cov <- raw |>
  summarise(coverage = mean(coverage), .by = c(sample, bin_id)) |>
  select(-bin_id)

curve_sample <- sample_cov |> fingerprint_curve(.by = "sample")
hist_all <- sample_cov |> coverage_hist(bin_width)
hist_low <- sample_cov |> coverage_hist(bin_width, max_cov = low_cov_max)
# NOTE: summary statistics describe the per-bin replicate-mean table (e.g.
# median_cov = median of fractional means), which matches the pooled-reads
# curve but not absolute pooled values — keep that in mind when comparing runs
summ     <- sample_cov |>
  fingerprint_summary(top_fracs) |>
  print()

# warn when a sample vanishes from the low-coverage zoom (all its bins are
# >= low_cov_max), which would otherwise read as a missing curve
lowcov_dropped <- setdiff(unique(sample_cov$sample), unique(hist_low$sample))
if (length(lowcov_dropped) > 0) {
  warning("Samples with all bins >= low_cov_max (", low_cov_max,
          ") are absent from the low-cov plot: ",
          paste(lowcov_dropped, collapse = ", "))
}

# ---- plots -------------------------------------------------------------------
stamp <- format(Sys.Date(), "%Y-%m-%d")

p1 <- curve_sample |>
  thin_rows() |>
  ggplot(aes(frac_genome, cumfrac_reads, color = sample)) +
  common_style("Fraction of genome bins (ordered by descending coverage)",
               "Cumulative Fraction of Reads")
ggsave(file.path(out_dir, paste0(stamp, "_cdf.pdf")), p1, width = 5, height = 3)

p_ecdf     <- cov_ecdf_figure(hist_all)
p_ecdf_low <- cov_ecdf_figure(hist_low, str_c("bins with coverage <", low_cov_max))
ggsave(file.path(out_dir, paste0(stamp, "_cdf_cov_ecdf.pdf")), p_ecdf,
       width = 8, height = 3)
ggsave(file.path(out_dir, paste0(stamp, "_cdf_cov_ecdf_lowcov.pdf")), p_ecdf_low,
       width = 8, height = 3)

# ---- export ------------------------------------------------------------------
write_csv(summ, file.path(out_dir, paste0(stamp, "_fingerprint_summary.csv")))
curve_sample |>
  thin_rows() |>
  rio::export(file.path(out_dir, paste0(stamp, "_cdf.xlsx")))
hist_all |>
  rio::export(file.path(out_dir, paste0(stamp, "_cdf_cov_ecdf.xlsx")))
hist_low |>
  rio::export(file.path(out_dir, paste0(stamp, "_cdf_cov_ecdf_lowcov.xlsx")))

cat("Done. Outputs in:", out_dir, "\n")


