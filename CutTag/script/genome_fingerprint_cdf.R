# ==============================================================================
# genome_fingerprint_cdf.R
# Redraw deepTools plotFingerprint results (BAM_fingerprint.txt) as:
#   1. the fingerprint curve per sample (replicates averaged per bin; the CDF
#      is scale-invariant, so mean and sum pooling give the same curve)
#   2. the cumulative distribution of per-bin coverage on the ABSOLUTE axis
#      (bin-rank axis + log coverage axis; depth-DEPENDENT — deepening a
#      library shifts it right — reads as library depth / complexity QC),
#      full range and low-coverage zoom
#   3. the same cumulative distribution on a depth-INDEPENDENT axis: per-bin
#      coverage as fold of the uniform expectation mu = mean bin coverage
#      (bin counts scale with depth and so does mu, so (c+1)/mu is invariant;
#      fold = 1 marks "exactly uniform"), full range and below-uniform zoom
#   4. a numeric summary: zero-coverage bin fraction, median coverage, and the
#      share of reads captured by the top 1 / 5 / 10 percent of bins
#
# Outputs (into out_dir, the parent of result_dir):
#   <date>_cdf_all.pdf    multi-page: fingerprint curve, absolute-coverage
#                         ECDF (full + low-cov zoom), depth-normalised fold
#                         ECDF (full + below-uniform zoom)
#   <date>_cdf_all.xlsx   all tables as sheets: summary, fingerprint curve,
#                         and the four cumulative distribution tables
# ==============================================================================

library(tidyverse)

# ---- config ------------------------------------------------------------------
result_dir       <- "D:\\EED_CutTag\\X101SC26064008-Z01-J011\\mice\\result\\multiqc"
out_dir          <- dirname(result_dir)   # plots / tables go next to the multiqc folder
fingerprint_file <- "BAM_fingerprint.txt"

# sample name pattern, adjust per project (column names after clean_names):
# sample_pattern <- c(sample = ".*", "[_]", gender=".*", "[_]", rep = "\\d+")   # e.g. d1_1
sample_pattern <- c(sample = ".*", "[-_]", rep = "\\d+") 
grouping <- setdiff(str_subset(names(sample_pattern), ".+"), "rep")
# sample_pattern <- c(sample = ".*", "[-_]", rep = "\\d+")           # e.g. MN0_1

# sample_levels <- c("d1", "d5", "l", "m", "ctrl", "nmn")  # plot/legend order; NULL = alphabetical
sample_levels <- NULL

bin_width  <- 10      # coverage bin width for the cumulative coverage plot
low_cov_max <- 10000  # coverage cutoff for the low-coverage zoom
fold_log_bin <- 0.05  # log10 bin width for the depth-normalised fold axis
bg_fold_max  <- 1     # fold cutoff for the below-uniform (background) zoom
thin_by    <- 200     # keep every n-th point when drawing curves
top_fracs  <- c(0.01, 0.05, 0.1)  # "signal in top x% of bins" summary metrics

# ---- functions ---------------------------------------------------------------

# Read deepTools plotFingerprint raw counts into long format: one row per
# replicate x genome bin
read_fingerprint <- function(path, sample_pattern) {
  # path <- file.path(result_dir, fingerprint_file)
  read_delim(path, delim = "\t", skip = 1) |>
    janitor::clean_names() |>
    mutate(bin_id = row_number()) |>
    pivot_longer(-bin_id, names_to = "name", values_to = "coverage") |>
    separate_wider_regex(name, patterns = sample_pattern, too_few = "align_start")
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
    mutate(coverage_bin = floor(coverage / bin_width) * bin_width, .by = grouping) |>
    summarise(n_bins = n(), .by = c(grouping, coverage_bin)) |>
    arrange(grouping, coverage_bin) |>
    mutate(cumfrac_bins = cumsum(n_bins) / sum(n_bins),
           bin_rank = row_number(), .by = grouping) |>
    ungroup()
}

# Numeric fingerprint summary: signal concentration per sample
fingerprint_summary <- function(coverage_tbl, top_fracs = c(0.01, 0.05, 0.1)) {
  coverage_tbl |>
    arrange(grouping, desc(coverage)) |>
    mutate(rank = row_number(), .by = grouping) |>
    group_by(grouping) |>
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
  sampled <- slice(tbl, seq(1, n(), by = by), .by = grouping)
  last    <- slice(tbl, n(), .by = grouping)          # keep the final CDF point
  bind_rows(sampled, last) |>
    distinct(grouping, frac_genome, .keep_all = TRUE) |>
    arrange(grouping, frac_genome)
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
  p_rank <- ggplot(hist_tbl, aes(bin_rank, cumfrac_bins, color = grouping)) +
    common_style("Coverage bin rank (ascending)", "Cumulative Fraction of Bins")
  p_cov  <- ggplot(hist_tbl, aes(coverage_bin + 1, cumfrac_bins, color = grouping)) +
    common_style("Coverage per bin (log scale)", "Cumulative Fraction of Bins") +
    scale_x_log10()
  patchwork::wrap_plots(p_rank, p_cov, nrow = 1) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = subtitle)
}

# Depth-normalised coverage distribution: per-bin coverage expressed as fold
# of the uniform expectation mu = mean bin coverage. Bin counts scale with
# sequencing depth (c_i ~ N * q_i) and so does mu, hence the ECDF of
# (c_i + 1) / mu is depth-INDEPENDENT — deepening a library rescales
# numerator and denominator together. The +1 keeps zero-coverage bins on the
# log axis, where they pile up at 1/mu = n_bins / N_reads. x is binned on
# log10 fold, so the axis is multiplicative.
coverage_fold_hist <- function(coverage_tbl, log_bin = fold_log_bin,
                               max_fold = Inf) {
  coverage_tbl |>
    mutate(mu   = sum(coverage) / n(),
           fold = (coverage + 1) / mu, .by = grouping) |>
    filter(fold < max_fold) |>
    mutate(fold_bin = 10^(floor(log10(fold) / log_bin) * log_bin),
           .by = grouping) |>
    summarise(n_bins = n(), .by = c(grouping, fold_bin)) |>
    arrange(grouping, fold_bin) |>
    mutate(cumfrac_bins = cumsum(n_bins) / sum(n_bins),
           bin_rank = row_number(), .by = grouping) |>
    ungroup()
}

# Single-panel cumulative distribution on the depth-normalised fold axis;
# the dashed line marks the uniform expectation (fold = 1)
fold_figure <- function(fold_tbl, subtitle = NULL) {
  ggplot(fold_tbl, aes(fold_bin, cumfrac_bins, color = grouping)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    common_style("Coverage per bin as fold of uniform expectation (log scale)",
                 "Cumulative Fraction of Bins") +
    scale_x_log10() +
    patchwork::plot_annotation(title = subtitle)
}

# ---- run ---------------------------------------------------------------------
raw <- read_fingerprint(file.path(result_dir, fingerprint_file), sample_pattern)
raw <- apply_sample_levels(raw) %>% 
  unite("grouping", any_of(grouping))

# Sample-level table: replicates averaged per bin (the CDF is scale-invariant,
# so this equals pooling the reads)
sample_cov <- raw |>
  summarise(coverage = mean(coverage), 
            .by = c(grouping, bin_id)) |>
  select(-bin_id)

curve_sample <- sample_cov |> fingerprint_curve(.by = "grouping")
hist_all <- sample_cov |> coverage_hist(bin_width)
hist_low <- sample_cov |> coverage_hist(bin_width, max_cov = low_cov_max)
hist_norm     <- sample_cov |> coverage_fold_hist()
hist_norm_bg  <- sample_cov |> coverage_fold_hist(max_fold = bg_fold_max)
# NOTE: summary statistics describe the per-bin replicate-mean table (e.g.
# median_cov = median of fractional means), which matches the pooled-reads
# curve but not absolute pooled values — keep that in mind when comparing runs
summ     <- sample_cov |>
  fingerprint_summary(top_fracs) |>
  print()

# warn when a sample vanishes from the low-coverage zoom (all its bins are
# >= low_cov_max), which would otherwise read as a missing curve
lowcov_dropped <- setdiff(unique(sample_cov$grouping), unique(hist_low$grouping))
if (length(lowcov_dropped) > 0) {
  warning("Samples with all bins >= low_cov_max (", low_cov_max,
          ") are absent from the low-cov plot: ",
          paste(lowcov_dropped, collapse = ", "))
}

# ---- plots -------------------------------------------------------------------
stamp <- format(Sys.Date(), "%Y-%m-%d")

p1 <- curve_sample |>
  thin_rows() |>
  ggplot(aes(frac_genome, cumfrac_reads, color = grouping)) +
  common_style("Fraction of genome bins (ordered by descending coverage)",
               "Cumulative Fraction of Reads") +
  ggtitle("Genome fingerprint: signal concentration (depth-independent)")

p_ecdf     <- cov_ecdf_figure(hist_all,
                              "all bins (absolute coverage; depth-dependent)")
p_ecdf_low <- cov_ecdf_figure(hist_low,
                              str_c("low-coverage zoom: bins with coverage <",
                                    low_cov_max))
p_norm     <- fold_figure(hist_norm,
                          "all bins (fold of uniform expectation; depth-independent)")
p_norm_low <- fold_figure(hist_norm_bg,
                          str_c("below-uniform zoom: fold <", bg_fold_max))

# all five pages merged into one multi-page PDF
pdf(file.path(out_dir, paste0(stamp, "_cdf_all.pdf")), width = 8, height = 3.5)
print(p1)
print(p_ecdf)
print(p_ecdf_low)
print(p_norm)
print(p_norm_low)
dev.off()

# ---- export ------------------------------------------------------------------
export_tbl <- list(
  summary           = summ,
  fingerprint_curve = curve_sample |> thin_rows(),
  cov_ecdf          = hist_all,
  cov_ecdf_lowcov   = hist_low,
  cov_ecdf_norm     = hist_norm,
  cov_ecdf_norm_bg  = hist_norm_bg
)
rio::export(export_tbl, file.path(out_dir, paste0(stamp, "_cdf_all.xlsx")))

cat("Done. Outputs in:", out_dir, "\n")


