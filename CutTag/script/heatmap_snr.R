# ==============================================================================
# heatmap_snr.R
# Compute signal-to-noise ratio (center / flank) from deepTools computeMatrix
# .mat.gz files (reference-point mode: TSS or center).
#
#   SNR(region) = sum(center bins) + pseudocount
#               / sum(flank bins)  + pseudocount
#
# where "center" = bins within +-center_width bp of the reference point and
# "flank" = the outermost flank_width bp on each side of the window.
#
# Run from RStudio: set the config section below, then Ctrl+A / Run (Source).
# Outputs (into out_dir):
#   heatmap_snr_summary.csv      one row per sample (median/mean/IQR of SNR)
#   heatmap_snr_per_region.csv   long table: region x sample center/flank/SNR
#   heatmap_snr_boxplot.pdf      log2(SNR) per sample, grouped by prefix
# ==============================================================================

library(tidyverse)
library(data.table)
library(jsonlite)

# ============================== config ======================================
mat_file <- "D:\\EED_CutTag\\c_elegans\\result\\multiqc\\consensus_H3K27me3_H3K27me3_center_CPM_mat.gz"
out_dir  <- "D:\\EED_CutTag\\c_elegans\\result\\multiqc"

center_width <- 1000    # bp around reference point = signal
flank_width  <- 4000   # bp at outermost window edges = background
pseudocount  <- 1      # added to center and flank sums (avoids Inf / NaN)
write_per_region <- TRUE  # also write per-region table (large: n_regions x n_samples rows)
# =============================================================================

# ============================ parse .mat.gz =================================
# New deepTools format: one JSON line starting with '@'
#   @{"upstream":[...],"downstream":[...],"bin size":[...],"ref point":[...],
#     "sample_labels":[...],"sample_boundaries":[0,600,...], ...}
# 'sample_boundaries' = 0-based column offsets of the signal values per sample.
# Legacy format: one @key value line per header entry
#   @group_boundaries [0,600,...]  (column offsets per sample)
#   @sample_labels [s1, s2, ...]   @upstream 3000  @downstream 3000  @bin_size 10

read_mat_meta <- function(mat_file) {
  con <- gzfile(mat_file, open = "rt")
  on.exit(close(con))
  header_lines <- character()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) stop("No data rows found in matrix file")
    if (!startsWith(line, "@")) break
    header_lines <- c(header_lines, line)
  }
  n_header <- length(header_lines)

  if (grepl("\\{", header_lines[1])) {
    h <- fromJSON(sub("^@", "", header_lines[1]))
    list(
      n_header      = n_header,
      sample_labels = h$sample_labels,
      sample_bounds = as.integer(h$sample_boundaries),   # 0-based, per sample
      upstream      = as.integer(h$upstream),
      downstream    = as.integer(h$downstream),
      bin_size      = as.integer(h$`bin size`),
      ref_point     = h$`ref point`
    )
  } else {
    get_val <- function(key, as_int = TRUE) {
      ln <- header_lines[startsWith(header_lines, paste0("@", key))]
      if (length(ln) == 0) return(NA)
      v <- sub(paste0("^@", key, "[[:space:]]*"), "", ln)
      if (grepl("^\\[", v)) {
        v <- str_remove_all(v, "^\\[|\\]$|'|\"")
        out <- as.integer(strsplit(v, ",\\s*")[[1]])
      } else if (as_int) out <- as.integer(v) else out <- v
      out
    }
    list(
      n_header      = n_header,
      sample_labels = get_val("sample_labels", as_int = FALSE),
      sample_bounds = get_val("group_boundaries"),        # legacy: cols per sample
      upstream      = get_val("upstream"),
      downstream    = get_val("downstream"),
      bin_size      = get_val("bin_size"),
      ref_point     = get_val("ref_point", as_int = FALSE)
    )
  }
}

meta <- read_mat_meta(mat_file)
n_samples <- length(meta$sample_labels)
if (length(meta$sample_bounds) != n_samples + 1) {
  stop("sample_boundaries length does not match sample_labels")
}

# ---- read only the signal columns we need (center + flank of every sample) --
# first 6 data columns are coordinates: chrom, start, end, name, score, strand
coord_cols <- 1:6
needed <- coord_cols
per_sample_cols <- list()
for (i in seq_len(n_samples)) {
  c0 <- meta$sample_bounds[i] + 1        # 1-based, within signal columns
  c1 <- meta$sample_bounds[i + 1]
  bsize  <- meta$bin_size[i]
  up     <- meta$upstream[i]
  nbin   <- c1 - c0 + 1

  ref_col <- up / bsize + 1              # reference point column (1-based)
  n_center <- round(center_width / bsize)
  n_flank  <- round(flank_width / bsize)
  center_cols <- (ref_col - n_center):(ref_col + n_center)
  flank_cols  <- c(seq_len(n_flank), (nbin - n_flank + 1):nbin)
  # guard against windows narrower than the requested widths
  center_cols <- center_cols[center_cols %in% seq_len(nbin)]
  flank_cols  <- flank_cols[flank_cols %in% seq_len(nbin)]

  abs_center <- center_cols + c0 - 1 + 6 # absolute column in the data table
  abs_flank  <- flank_cols + c0 - 1 + 6
  needed <- c(needed, abs_center, abs_flank)
  per_sample_cols[[i]] <- list(center = abs_center, flank = abs_flank)
}
needed <- sort(unique(needed))
# fread(select=) renumbers columns: remap file positions -> dt column indices
per_sample_cols <- lapply(per_sample_cols, function(x)
  list(center = match(x$center, needed), flank = match(x$flank, needed)))

cat(sprintf("Reading %s (header=%d, samples=%d, signal cols=%d)...\n",
            basename(mat_file), meta$n_header, n_samples, length(needed) - 6))
dt <- fread(mat_file, skip = meta$n_header, header = FALSE,
            select = needed, na.strings = c("nan", "NA", "NaN"))
dt[is.na(dt)] <- 0
cat(sprintf("Rows read: %d\n", nrow(dt)))

# ============================ compute SNR ====================================
res <- vector("list", n_samples)
for (i in seq_len(n_samples)) {
  c_cols <- per_sample_cols[[i]]$center
  f_cols <- per_sample_cols[[i]]$flank
  sub_center <- dt[, ..c_cols]
  sub_flank  <- dt[, ..f_cols]
  center_sum <- rowSums(as.matrix(sub_center))  # as.matrix: handles 1-column case
  flank_sum  <- rowSums(as.matrix(sub_flank))
  n_c <- length(c_cols)
  n_f <- length(f_cols)
  res[[i]] <- tibble(
    sample      = meta$sample_labels[i],
    chrom       = dt[[1]],
    start       = dt[[2]],
    end         = dt[[3]],
    center_mean = center_sum / n_c,
    flank_mean  = flank_sum / n_f,
    snr         = (center_sum + pseudocount) / (flank_sum + pseudocount)
  )
}
df <- bind_rows(res) |> print()

# ============================ summary + output ==============================
summary_tbl <- df %>%
  summarise(
    n_regions     = n(),
    n_zero_flank  = sum(flank_mean == 0),
    median_snr    = median(snr),
    mean_snr      = mean(snr),
    sd_snr        = sd(snr),
    q25_snr       = quantile(snr, 0.25),
    q75_snr       = quantile(snr, 0.75),
    .by = sample
  ) |>
  arrange(desc(median_snr)) |>
  print()

write_csv(summary_tbl, file.path(out_dir, "heatmap_snr_summary.csv"))
if (write_per_region) {
  fwrite(df, file.path(out_dir, "heatmap_snr_per_region.csv"))
}

# boxplot of log2(SNR) per sample, colored by sample name prefix (group)
p <- df %>%
  mutate(
    snr2  = log2(snr),
    group = replace_na(str_extract(sample, "^[^0-9-]+"), "NA")
  ) %>%
  ggplot(aes(x = sample, y = snr, fill = group)) +
  geom_boxplot(outlier.size = 0.3, outliers = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = str_glue("log2 SNR ({basename(mat_file)})"),
       x = NULL, y = "log2(SNR)")
p
ggsave(file.path(out_dir, "heatmap_snr_boxplot.pdf"), p, width = 14, height = 7)
cat("Done. Outputs written to:", out_dir, "\n")
