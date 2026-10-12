# ==============================================================================
# heatmap_snr.R
# Compute signal-to-noise ratio (center / flank) from deepTools-style count
# matrices (.mat.gz, reference-point mode: TSS or center).
#
#   SNR(region) = mean(center bins) + pseudocount
#               / mean(flank bins)  + pseudocount
#
# Regions are fractions of the whole matrix window (works for any bin size /
# upstream / downstream): "center" = central `center_frac` of the window,
# centered on the reference point; "flank" = the remaining fraction, split
# evenly between the outermost left and right window edges
# (center_frac 0.2 -> flank 0.4 + 0.4, always disjoint from the center).
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
mat_file <- "E:\\EED\\2026-08-25_CutTag_cel\\result\\multiqc\\consensus_H3K27me3_H3K27me3_center_CPM_mat.gz"
out_dir  <- "E:\\EED\\2026-08-25_CutTag_cel\\result"

center_frac <- 0.2     # fraction of the window used as signal, centered on ref point
                       # flank = the rest, split evenly: (1 - center_frac)/2 per side
pseudocount  <- 0.0001   # added to center and flank bin means (avoids Inf / NaN)
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
        parts <- strsplit(v, ",\\s*")[[1]]
        # only coerce bracket lists to integer when the caller wants numbers —
        # as.integer on string labels (sample_labels) silently produced NA
        out <- if (as_int) as.integer(parts) else parts
      } else if (as_int) out <- as.integer(v) else out <- v
      out
    }
    bounds <- get_val("sample_boundaries")   # column offsets per sample
    if (length(bounds) == 1 && is.na(bounds)) {
      # group_boundaries delimits REGION groups (matrix rows), NOT the
      # per-sample signal columns — using it would put center/flank windows
      # on entirely wrong columns. Derive the column layout from the matrix
      # geometry instead, or refuse the file.
      ups <- get_val("upstream"); dws <- get_val("downstream")
      bs  <- get_val("bin size")
      labels <- get_val("sample_labels", as_int = FALSE)
      if (anyNA(ups) || anyNA(dws) || anyNA(bs) || length(labels) == 0 ||
          !(length(ups) == 1 && length(dws) == 1 && length(bs) == 1)) {
        stop("legacy matrix without sample_boundaries and without a usable ",
             "upstream/downstream/bin size triple; regenerate the matrix ",
             "with a current deepTools computeMatrix")
      }
      nbin <- as.integer((ups + dws) %/% bs)
      bounds <- as.integer(cumsum(c(0L, rep(nbin, length(labels)))))
      warning("legacy matrix: sample_boundaries missing; derived per-sample ",
              "column offsets from upstream/downstream/bin size (",
              nbin, " bins/sample)")
    }
    list(
      n_header      = n_header,
      sample_labels = get_val("sample_labels", as_int = FALSE),
      sample_bounds = bounds,
      upstream      = get_val("upstream"),
      downstream    = get_val("downstream"),
      bin_size      = get_val("bin size"),
      ref_point     = get_val("ref point", as_int = FALSE)
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
# upstream / bin size are scalars in deepTools computeMatrix output but
# per-sample vectors in nf-core consensus matrices; normalize to per-sample
ups    <- rep_len(meta$upstream, n_samples)
bsizes <- rep_len(meta$bin_size, n_samples)
for (i in seq_len(n_samples)) {
  c0 <- meta$sample_bounds[i] + 1        # 1-based, within signal columns
  c1 <- meta$sample_bounds[i + 1]
  bsize  <- bsizes[i]
  up     <- ups[i]
  nbin   <- c1 - c0 + 1

  ref_col <- up / bsize + 1              # reference point column (1-based)
  h <- floor(center_frac / 2 * nbin)     # center = +-h bins around ref point
  center_cols <- (ref_col - h):(ref_col + h)
  center_cols <- center_cols[center_cols %in% seq_len(nbin)]
  if (length(center_cols) == 0) {
    stop(sprintf("sample %s: center fraction does not fit a %d-bin window",
                 meta$sample_labels[i], nbin))
  }

  n_flank <- floor((1 - center_frac) / 2 * nbin)  # flank bins per side
  c_first <- min(center_cols)
  c_last  <- max(center_cols)
  left_flank  <- seq_len(min(n_flank, c_first - 1))    # outermost left bins
  right_start <- max(c_last + 1, nbin - n_flank + 1)   # outermost right bins
  right_flank <- if (right_start <= nbin) right_start:nbin else integer(0)
  flank_cols  <- c(left_flank, right_flank)
  if (length(left_flank) == 0 || length(right_flank) == 0) {
    stop(sprintf("sample %s: flank fraction does not fit beside the center in a %d-bin window",
                 meta$sample_labels[i], nbin))
  }

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
    center_nbins = n_c,
    flank_nbins  = n_f,
    center_mean  = center_sum / n_c,
    flank_mean   = flank_sum / n_f,
    snr          = (center_mean + pseudocount) / (flank_mean + pseudocount)
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

# flag zero-flank regions: with the tiny pseudocount their SNR explodes and
# can dominate mean_snr — report median/IQR for such samples
if (any(summary_tbl$n_zero_flank > 0)) {
  warning("zero-flank regions present; mean_snr is inflated for those samples — prefer median/q25/q75")
}

write_csv(summary_tbl, file.path(out_dir, "heatmap_snr_summary.csv"))
if (write_per_region) {
  fwrite(df, file.path(out_dir, "heatmap_snr_per_region.csv"))
}

# boxplot of log2(SNR) per sample, colored by sample name prefix (group)
p <- df %>%
  mutate(group = replace_na(str_extract(sample, "^[^0-9-]+"), "NA")) %>%
  ggplot(aes(x = sample, y = log2(snr), fill = group)) +
  # outlier.shape = NA is the ggplot2 way to hide outliers (the argument
  # `outliers` does not exist and its extremes were silently kept plotted)
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = str_glue("log2 SNR ({basename(mat_file)})"),
       x = NULL, y = "log2(SNR)")
p
ggsave(file.path(out_dir, "heatmap_snr_boxplot.pdf"), p, width = 14, height = 7)
cat("Done. Outputs written to:", out_dir, "\n")



