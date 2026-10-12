# ==============================================================================
# peak_counts_ratio_v3.R
# Peak-vs-background enrichment of CUT&Tag fragment counts.
#
# For every sample the filtered BAM is read once and each region set is
# counted once; every repetition then bootstrap-resamples sample_ratio of
# the regions (WITH replacement) and just sums the precomputed per-region
# counts, so repetitions are cheap and the spread across reps is a real
# bootstrap distribution.
#
# Region sets: per-sample MACS peaks, consensus peaks, genes, TSS windows.
# Background: equal-width adjacent flanks (flank_scale * region width per
# side), with the region set itself punched out — a neighbour of the same
# set falling inside a flank is contamination, not background. A fragment
# counts as background iff its midpoint falls on the leftover flank, so it
# is counted exactly once even when it straddles a punched-out region.
# Bases lost to chromosome edges are compensated by rescaling the cleaned
# count to the requested flank width.
#
#   ratio_flank = (in + pc) / (flank count normalised to the regions' width + pc)
#       Density ratio: flank counts are normalised to the sampled regions'
#       total width using the actual kept flank width, so the value is
#       INDEPENDENT of flank_scale (it only sets how local the background
#       is) and a width-matched null region scores ~1. Regions whose flanks
#       are entirely covered by same-set neighbours give NA (no estimable
#       background).
#
# Outputs (written into result_dir):
#   <date>_peak_ratio_per_rep.csv   one row per sample x rep x region set
#   <date>_peak_ratio_summary.csv   per sample x set: n regions, mean/median
#                                   flank ratio
#   <date>_peak_ratio_samples.csv   per sample: n_fragments, FRiP over own
#                                   MACS peaks, FRiP over the consensus set
#   <date>_peak_ratio.pdf           violin/boxplot per set
# ==============================================================================

library(tidyverse)
library(GenomicAlignments)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(furrr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)   # hs
library(TxDb.Mmusculus.UCSC.mm39.knownGene)  # mm
library(TxDb.Celegans.UCSC.ce11.refGene)     # celegans
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene) # fly

# ---- config ------------------------------------------------------------------
result_dir    <- "E:/EED/2026-07-21_CutTag_cel_fly_mm/fly/result"
bam_dir       <- "03_bam"
# peak_dir      <- "04_peaks/macs/peak"        # or "04_peaks/seacr/peak"
peak_dir <- "04_peaks"
peak_pattern  <- "_peaks\\.broadPeak$"       # or "_stringent\\.bed$" for seacr
# consensus_bed <- "04_peaks/macs/consensus_H3K27me3.bed"
consensus_bed <- "04_peaks/consensus_H3K27me3.bed"
# txdb          <- TxDb.Celegans.UCSC.ce11.refGene
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene

sample_subset <- NULL  # NULL = all
# sample_subset <- c("MN0_1", "MN10_1")  # two samples for testing;
flank_scale   <- 2        # flank width = flank_scale * region width, per side
sample_ratio  <- 0.30     # fraction of regions bootstrap-resampled WITH
                          # replacement per repetition
n_reps        <- 100
pseudocount   <- 1        # per-side pseudocount in ratio_flank:
                          # ratio_flank = (in + pc) / (flank_per_side + pc);
                          # zero-flank regions give NA (unestimable).
set_seed      <- 20260923
n_workers     <- 4        # 1 = sequential (low memory, clear errors); 4-6 for full runs
write_per_rep <- TRUE

options(future.globals.maxSize = 3 * 1024^3)

# ---- functions ---------------------------------------------------------------

# Import peak/BED regions, keep standard chromosomes, use UCSC seqname style
import_regions <- function(path) {
  gr <- rtracklayer::import(path)
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  seqlevelsStyle(gr) <- "UCSC"
  sort(gr)
}

# Paired-end fragments of a filtered BAM as a GRanges
read_fragments <- function(bam_path) {
  gr <- readGAlignmentPairs(bam_path) |> granges()
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  seqlevelsStyle(gr) <- "UCSC"
  sort(gr)
}

# Reduced gene bodies and TSS windows from a TxDb. TSS windows are built from
# the ALREADY reduced genes (purrr::lst was sequential, but it's gone from
# newer purrr — spell the dependency out with an intermediate variable).
build_annotation <- function(txdb) {
  gene <- genes(txdb) |>
    keepStandardChromosomes(pruning.mode = "coarse")
  seqlevelsStyle(gene) <- "UCSC"
  sl   <- seqlengths(txdb)
  gene <- set_genome(gene, sl[!is.na(sl) & sl > 0])
  gene <- reduce(sort(gene))
  list(gene = gene,
       tss  = reduce(sort(suppressWarnings(trim(promoters(gene, 1000, 1000))))))
}

# Restrict to chromosomes with known length and attach seqlengths (for trim())
set_genome <- function(gr, seqlens) {
  gr <- gr[as.character(seqnames(gr)) %in% names(seqlens)]
  seqlengths(gr) <- seqlens[as.character(seqlevels(gr))]
  gr
}

# Flank windows of gr, clipped to chromosome bounds, with gr itself punched
# out: a same-set neighbour inside a flank is contamination, not background.
# scale_* compensates ONLY chromosome-edge loss (requested / clipped width);
# the punch-out is a genuine removal and is not rescaled back.
make_flanks <- function(gr, flank_scale = 1) {
  w    <- pmax(1L, floor(flank_scale * width(gr)))
  roi  <- reduce(gr, ignore.strand = TRUE)
  strand(roi) <- "*"
  fl_L <- suppressWarnings(trim(flank(gr, w, start = TRUE)))
  fl_R <- suppressWarnings(trim(flank(gr, w, start = FALSE)))
  list(
    flank_L = fl_L,
    flank_R = fl_R,
    scale_L = w / pmax(1L, width(fl_L)),
    scale_R = w / pmax(1L, width(fl_R)),
    w       = w,
    roi     = roi
  )
}

# Per-region fragment count over flank \ roi (1:1 with the flanks' regions).
# A fragment counts as background iff its MIDPOINT falls on a leftover piece:
# overlap counting cannot be subtracted piecewise (a fragment straddling a
# hole would be lost by inclusion-exclusion, or double-counted per piece),
# while the midpoint attributes every fragment exactly once.
# psetdiff keeps the 1:1 region mapping that a plain setdiff() would destroy.
# Returns list(cnt, width): cnt is edge-rescaled; width is the EFFECTIVE
# flank width = requested width - punched-out width (chromosome-edge loss is
# already restored in cnt via the scale factor, so it is not deducted here).
count_flank <- function(frag_idx, fl, side = c("L", "R")) {
  side     <- match.arg(side)
  flank_gr <- if (side == "L") fl$flank_L else fl$flank_R
  scale    <- if (side == "L") fl$scale_L   else fl$scale_R
  n        <- length(flank_gr)
  hits <- findOverlaps(flank_gr, fl$roi, ignore.strand = TRUE)
  q <- flank_gr[queryHits(hits)]
  s <- fl$roi[subjectHits(hits)]
  strand(q) <- "*"
  strand(s) <- "*"
  holes <- pintersect(q, s)
  keep  <- width(holes) > 0L
  holes <- holes[keep]
  qh    <- queryHits(hits)[keep]
  punched <- integer(n)
  if (length(holes) > 0L) {
    pw <- rowsum(width(holes), qh)
    punched[as.integer(rownames(pw))] <- as.integer(pw[, 1])
  }
  strand(flank_gr) <- "*"
  kept  <- psetdiff(flank_gr, split(holes, factor(qh, levels = seq_len(n))))
  piece <- unlist(kept, use.names = FALSE)
  reg   <- rep(seq_len(n), lengths(kept))
  ok    <- width(piece) > 0L
  piece <- piece[ok]
  reg   <- reg[ok]
  cnt <- integer(n)
  rs  <- rowsum(count_fragment_mids(frag_idx, piece), reg)
  cnt[as.integer(rownames(rs))] <- as.integer(rs[, 1])
  list(cnt   = as.integer(cnt * scale),
       width = pmax(0L, as.integer(fl$w) - punched))
}

# Number of fragments overlapping each region of `gr`. Equivalent to
# countOverlaps(gr, bam) but computed with findInterval on the sorted fragment
# coordinates, which is much faster for millions of fragments:
# a fragment [fs, fe] overlaps [s, e] iff fs <= e and fe >= s, hence
#   count = #(fs <= e) - #(fe < s)
#
# The sorted per-chromosome coordinates are prepared ONCE per sample
# (prepare_fragment_index): count_fragments runs several times per region set,
# and re-splitting + re-sorting millions of coordinates on every call was the
# dominant avoidable cost of this script.
prepare_fragment_index <- function(bam) {
  chr <- as.character(seqnames(bam))
  # fragment midpoints: the flank background attributes each fragment by its
  # midpoint so a fragment straddling a punched-out region is neither lost
  # nor double-counted (see count_flank)
  mid <- start(bam) + (width(bam) - 1L) %/% 2L
  list(
    st  = lapply(split(start(bam), chr), sort),
    en  = lapply(split(end(bam),   chr), sort),
    mid = lapply(split(mid, chr), sort)
  )
}

count_fragments <- function(frag_idx, gr) {
  st_by <- frag_idx$st
  en_by <- frag_idx$en
  chr   <- as.character(seqnames(gr))
  vapply(seq_along(gr), \(i) {
    st <- st_by[[chr[i]]]
    if (is.null(st)) return(0L)
    as.integer(findInterval(end(gr)[i], st) -
                 findInterval(start(gr)[i] - 1L, en_by[[chr[i]]]))
  }, integer(1))
}

# number of fragment MIDPOINTS inside each region (same findInterval trick)
count_fragment_mids <- function(frag_idx, gr) {
  mid_by <- frag_idx$mid
  chr    <- as.character(seqnames(gr))
  vapply(seq_along(gr), \(i) {
    mid <- mid_by[[chr[i]]]
    if (is.null(mid)) return(0L)
    as.integer(findInterval(end(gr)[i], mid) -
                 findInterval(start(gr)[i] - 1L, mid))
  }, integer(1))
}

# One bootstrap draw: sample_ratio of each region set, WITH replacement.
# ratio_flank is a DENSITY ratio: flank counts are normalised to the sampled
# regions' total width using the actual kept flank width, so it is
# independent of flank_scale and a width-matched null region scores ~1.
resample_ratios <- function(counted) {
  counted |>
    group_by(type) |>
    group_modify(\(d, ...) {
      k   <- max(1L, floor(nrow(d) * sample_ratio))
      itr <- sample.int(nrow(d), k, replace = TRUE)
      in_cnt    <- sum(d$in_cnt[itr])
      flank_cnt <- sum(d$l_cnt[itr]) + sum(d$r_cnt[itr])
      w_reg     <- sum(d$w[itr])
      w_flank   <- sum(d$l_w[itr]) + sum(d$r_w[itr])
      # float division first: flank_cnt * w_reg overflows 32-bit integers
      flank_norm <- if (w_flank == 0) NA_real_
                    else (as.numeric(flank_cnt) / w_flank) * w_reg
      tibble(
        n_sampled   = k,
        in_cnt      = in_cnt,
        flank_cnt   = flank_cnt,
        ratio_flank = if (is.na(flank_norm) || flank_norm == 0) NA_real_
                      else (in_cnt + pseudocount) / (flank_norm + pseudocount)
      )
    }) |>
    ungroup()
}

# Per-sample pipeline: read the BAM once, count each region set once, then
# bootstrap from the precomputed per-region counts
analyse_sample <- function(bam_path, peak_path, consensus_gr, annotation, seqlens) {
  message("processing: ", basename(bam_path))
  bam   <- read_fragments(bam_path)
  peaks <- import_regions(peak_path)
  if (length(peaks) == 0) stop("No peaks found in ", basename(peak_path))
  if (length(consensus_gr) == 0) stop("Consensus peak set is empty")

  sets <- list(peak = peaks, consensus = consensus_gr,
               gene = annotation$gene, tss = annotation$tss)

  frag_idx <- prepare_fragment_index(bam)
  counted <- imap(sets, \(gr, type) {
    gr <- set_genome(gr, seqlens)
    fl <- make_flanks(gr, flank_scale)
    message(sprintf("  %s: %d regions", type, length(gr)))
    fl_l <- count_flank(frag_idx, fl, "L")
    fl_r <- count_flank(frag_idx, fl, "R")
    tibble(type   = type,
           in_cnt = count_fragments(frag_idx, gr),
           w      = width(gr),
           l_cnt  = fl_l$cnt, l_w = fl_l$width,
           r_cnt  = fl_r$cnt, r_w = fl_r$width)
  }) |>
    list_rbind()

  per_rep <- seq_len(n_reps) |>
    map(\(r) resample_ratios(counted) |> mutate(rep = r)) |>
    list_rbind()

  list(
    per_rep     = per_rep,
    # dplyr::count is masked by matrixStats::count via the SummarizedExperiment
    # load chain, which errors on a tibble argument — use the namespace form
    region_n    = dplyr::count(counted, type, name = "n_regions"),
    n_fragments = length(bam),
    # FRiP over the sample's own MACS peaks and over the consensus set
    frip_own_peak  = sum(count_fragments(frag_idx, peaks)) / length(bam),
    frip_consensus = sum(count_fragments(frag_idx, consensus_gr)) / length(bam)
  )
}

# ---- run ---------------------------------------------------------------------
set.seed(set_seed)

peak_files <- list.files(file.path(result_dir, peak_dir), pattern = peak_pattern)
samples <- tibble(
  sample    = str_replace(peak_files, peak_pattern, ""),
  peak_file = file.path(result_dir, peak_dir, peak_files),
  bam_file  = file.path(result_dir, bam_dir, str_c(sample, ".filtered.bam"))
)
if (!is.null(sample_subset)) samples <- filter(samples, sample %in% sample_subset)
stopifnot(nrow(samples) > 0,
          all(file.exists(samples$peak_file)),
          all(file.exists(samples$bam_file)))
print(samples, width = Inf)

seqlens <- {
  sl <- seqlengths(txdb)
  sl[!is.na(sl) & sl > 0]
}
consensus  <- import_regions(file.path(result_dir, consensus_bed))
annotation <- build_annotation(txdb)

if (n_workers == 1) {
  plan(sequential)
} else {
  plan(multisession, workers = min(n_workers, nrow(samples)))
}
res <- future_map(seq_len(nrow(samples)), \(i) {
  analyse_sample(samples$bam_file[i], samples$peak_file[i],
                 consensus, annotation, seqlens)
}, .options = furrr_options(seed = TRUE)) |>
  set_names(samples$sample)
plan(sequential)

# ---- output ------------------------------------------------------------------
stamp <- format(Sys.Date(), "%Y-%m-%d")

per_rep  <- imap(res, \(r, s) mutate(r$per_rep, sample = s)) |> list_rbind()
region_n <- imap(res, \(r, s) mutate(r$region_n, sample = s)) |> list_rbind()
summary_ratio <- per_rep |>
  group_by(sample, type) |>
  summarise(
    mean_flank   = mean(ratio_flank, na.rm = TRUE),
    median_flank = median(ratio_flank, na.rm = TRUE),
    n_na_flank   = sum(is.na(ratio_flank)),
    .groups      = "drop"
  ) |>
  left_join(region_n, by = c("sample", "type"))
sample_stats <- imap(res, \(r, s) {
  tibble(sample         = s,
         n_fragments    = r$n_fragments,
         frip_own_peak  = r$frip_own_peak,
         frip_consensus = r$frip_consensus)
}) |> list_rbind()

if (write_per_rep) {
  write_csv(per_rep, file.path(result_dir, paste0(stamp, "_peak_ratio_per_rep.csv")))
}
write_csv(summary_ratio, file.path(result_dir, paste0(stamp, "_peak_ratio_summary.csv")))
write_csv(sample_stats,  file.path(result_dir, paste0(stamp, "_peak_ratio_samples.csv")))

plot_tbl <- per_rep |>
  # strip the replicate suffix whether samples are named with "_" (MN0_1)
  # or "-" (CTRL-1)
  mutate(group = str_remove(sample, "[_-][^-_]*$"))

p <- ggplot(plot_tbl, aes(sample, ratio_flank, fill = group)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_violin() +
  geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~type, scales = "free_y") +
  labs(title = "Peak enrichment over local flank background",
       subtitle = paste("density ratio (width-normalised); flank =", flank_scale,
                        "x width per side, same-set regions punched out;",
                        "bootstrap", sample_ratio, "x", n_reps, "reps"),
       y = "density ratio (region / flank)", x = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(file.path(result_dir, paste0(stamp, "_peak_ratio.pdf")), p,
       width = 10, height = 6)

cat("Done. Outputs in:", result_dir, "\n")
