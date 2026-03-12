
library(GenomicRanges)
# library(consensusSeekeR)
library(EnrichedHeatmap)
library(rtracklayer)
library(circlize)
library(ComplexHeatmap)
library(grid)


#' Compute consensus peak regions from a list of peak calls
#'
#' Creates consensus peak sets from multiple peak calling results
#' (typically from ChIP-seq, ATAC-seq or similar experiments).
#' Supports grouping of samples and two different computation strategies.
#'
#' @param peak_list Named list of GRanges objects (each element = peak calls from one sample/condition)
#' @param groups Optional character vector of the same length as peak_list,
#'   specifying group membership for each sample. If NULL, all samples are
#'   treated as one group.
#' @param method Character. Peak consensus strategy to use:
#'   \itemize{
#'     \item{"granges"} {Simple coverage-based approach (faster, less precise)}
#'     \item{"consensusseeker"} {More sophisticated algorithm from consensusSeekeR package (slower, more accurate)}
#'   }
#' @param genome_build String or BSgenome object (required when method = "consensusseeker")
#' @param min_occurrence Numeric (only for "granges" method). Minimum number of peak files
#'   in which a region must be present to be included in consensus. Default: 2
#' @param max_coverage Numeric (only for "granges" method). Upper coverage threshold
#'   for slice (usually Inf). Rarely needs changing.
#' @param min_gapwidth Integer (only for "granges" method). Minimum gap width to merge
#'   nearby regions during the reduce step. Default: 1L (effectively merges touching regions)
#' @param ... Additional arguments passed to consensusSeekeR::findConsensusPeakRegions()
#'   when method = "consensusseeker"
#'
#' @return Named list of GRanges objects — one consensus peak set per group
#'
#' @note When method = "consensusseeker", genome_build is required.
#'       Common values: "hg38", "hg19", "mm10", "mm39", BSgenome.Hsapiens.UCSC.hg38, etc.
#'
#' @examples
#' \dontrun{
#' consensus <- compute_consensus_peaks(
#'   peak_list = list(CnR = CnR_peaks, CnT = CnT_peaks, ENCODE = encode_peaks),
#'   groups  = c("labA", "labA", "public"),
#'   method  = "granges",
#'   min_occurrence = 2,
#'   min_gapwidth   = 50
#' )
#' }
#'
#' @importFrom GenomicRanges coverage GRangesList GRanges reduce seqinfo seqlevelsInUse
#' @importFrom IRanges slice
#' @export
compute_consensus_peaks <- function(
    peak_list,
    groups = NULL,
    method = c("granges", "consensusseeker"),
    genome_build = NULL,
    min_occurrence = 2L,
    max_coverage   = Inf,
    min_gapwidth   = 1L,
    ...) {
  
  # ── Input validation ────────────────────────────────────────────────────────
  if (!inherits(peak_list, "list") || length(peak_list) == 0) {
    stop("peak_list must be a non-empty named list of GRanges objects")
  }
  
  if (is.null(names(peak_list)) || any(names(peak_list) == "")) {
    stop("peak_list must be named")
  }
  
  if (!all(vapply(peak_list, \(x) inherits(x, "GRanges"), logical(1)))) {
    stop("All elements of peak_list must be GRanges objects")
  }
  
  method <- match.arg(tolower(method[1]), c("granges", "consensusseeker"))
  
  if (method == "consensusseeker" && is.null(genome_build)) {
    stop("genome_build is required when method = 'consensusseeker'")
  }
  
  if (!is.null(groups)) {
    if (length(groups) != length(peak_list)) {
      stop("groups must have same length as peak_list or be NULL")
    }
    if (anyNA(groups)) {
      stop("groups vector contains NA values")
    }
  } else {
    groups <- rep("all", length(peak_list))
  }
  
  if (!is.numeric(min_occurrence) || length(min_occurrence) != 1 ||
      min_occurrence < 1 || min_occurrence > length(peak_list)) {
    stop("min_occurrence must be a number between 1 and length(peak_list)")
  }
  min_occurrence <- as.integer(min_occurrence)
  
  # ── Processing ──────────────────────────────────────────────────────────────
  message(sprintf("Computing consensus peaks using method = '%s' ...", method))
  
  consensus_by_group <- lapply(unique(groups), function(current_group) {
    
    idx <- which(groups == current_group)
    current_peaks <- peak_list[idx]
    
    n_samples <- length(current_peaks)
    if (n_samples < 2) {
      message(sprintf(
        "Group '%s' has only %d file → returning original peaks (no consensus computed)",
        current_group, n_samples
      ))
      return(current_peaks[[1]])
    }
    
    grl <- GenomicRanges::GRangesList(current_peaks)
    
    if (method == "granges") {
      
      cov <- GenomicRanges::coverage(grl)
      sliced <- IRanges::slice(
        x          = cov,
        lower      = min_occurrence,
        upper      = max_coverage,
        rangesOnly = TRUE
      )
      merged <- GenomicRanges::GRanges(sliced)
      reduced <- GenomicRanges::reduce(merged, min.gapwidth = min_gapwidth)
      
      return(reduced)
      
    } else {  # consensusseeker
      
      if (!requireNamespace("consensusSeekeR", quietly = TRUE)) {
        stop("Package 'consensusSeekeR' is required for method='consensusseeker'", call. = FALSE)
      }
      
      # Get seqinfo only for chromosomes actually used
      used_chr <- GenomeInfoDb::seqlevelsInUse(grl)
      chr_info <- GenomicRanges::seqinfo(genome_build)[used_chr]
      
      res <- consensusSeekeR::findConsensusPeakRegions(
        peaks    = unlist(grl, use.names = FALSE),
        chrInfo  = chr_info,
        ...
      )
      
      return(res$consensusRanges)
    }
  })
  
  names(consensus_by_group) <- unique(groups)
  
  message(sprintf("Consensus peak computation finished"))
  
  invisible(consensus_by_group)
}





get_region_and_flanks_coverage <- function(
    bam_path,
    gr,
    flank_size = 0,
    mode = c("count", "coverage"),
    ...
) {
  
  mode <- match.arg(mode)
  
  if (!file.exists(paste0(bam_path, ".bai"))) {
    stop("BAM index (.bai) not found.")
  }
  
  df <- data.frame(
    seqnames = as.character(seqnames(gr)),
    start    = start(gr),
    end      = end(gr),
    width    = width(gr),
    stringsAsFactors = FALSE
  )
  if (!is.null(names(gr)) && all(nzchar(names(gr)))) {
    df$name <- names(gr)
  }
  
  # ── 关键改动：不丢弃任何 flank，即使宽度<=0 ───────────────────────
  gr_region <- gr
  
  # left flank：允许 start<1
  if(flank_size > 0) {
    gr_left   <- flank(gr, width = flank_size, start = TRUE, both = FALSE)
    
    # right flank：允许 end > seqlengths
    gr_right  <- flank(gr, width = flank_size, start = FALSE, both = FALSE)
    
    # 合并所有要查询的区间（但仍分开统计）
    all_regions <- c(gr_region, gr_left, gr_right)
  } else {
    all_regions <- gr_region
  }
  
  param <- ScanBamParam(
    what = c("rname", "strand", "pos", "qwidth"),
    which = IRanges::reduce(all_regions, ignore.strand = TRUE),  # 合并重叠区间加速读取
    ...
  )
  
  message("Reading alignments...")
  gal <- readGAlignmentPairs(bam_path, param = param)
  
  message("Computing overlaps...")
  
  # 分别计数（即使宽度<=0，findOverlaps也会返回0）
  count_region <- countOverlaps(gr_region, gal, ignore.strand = TRUE)
  df$count_region <- count_region
  df$cov_region <- count_region / (width(gr_region) / 1000)
  
  if(flank_size > 0) {
    count_left   <- countOverlaps(gr_left,   gal, ignore.strand = TRUE)
    df$count_left   <- count_left
    df$cov_left   <- count_left   / (flank_size / 1000)
    
    count_right  <- countOverlaps(gr_right,  gal, ignore.strand = TRUE)
    df$count_right  <- count_right
    df$cov_right  <- count_right  / (flank_size / 1000)
  } 
  
  message("Done")
  return(df)
}




#' Create DeepTools-style Heatmap with Profile Plot
#'
#' This function mimics the behavior of deepTools' computeMatrix and plotHeatmap,
#' creating publication-quality heatmaps with profile plots for ChIP-seq, ATAC-seq,
#' or other genomic signal data. It supports two main modes:
#' 
#' - **scale_regions**: For visualizing signal across genomic features of varying lengths
#'   (e.g., gene bodies, enhancers). Input GRanges should have width > 1. The function
#'   will scale all regions to the same length defined by the 'body' parameter.
#'   
#' - **reference_point**: For visualizing signal around specific genomic positions
#'   (e.g., TSS, TES, peak summits). Input GRanges should have width == 1. The function
#'   anchors visualization at these points and extends upstream/downstream.
#'
#' @param bw_files Character vector of BigWig file paths. Multiple samples will be
#'   displayed as separate side-by-side heatmap panels.
#'   
#' @param granges_list List of GRanges objects, or a single GRanges object. All regions
#'   should follow the same mode convention (width > 1 for scale_regions, width == 1
#'   for reference_point). Multiple region sets will be concatenated and can be split
#'   using the 'region_names' parameter for visual grouping.
#'   
#' @param mode Character: Either "scale_regions" or "reference_point".
#'   - "scale_regions": Scales genomic regions (width > 1) to uniform length. Use for
#'     features like gene bodies where you want to compare signal patterns regardless
#'     of actual region length.
#'   - "reference_point": Anchors at a specific point (width == 1). Use for point
#'     features like TSS or peak summits where absolute distance matters.
#'     
#' @param upstream Integer: Number of base pairs to extend upstream of the region start
#'   (scale_regions mode) or reference point (reference_point mode). Default: 3000.
#'   
#' @param downstream Integer: Number of base pairs to extend downstream of the region end
#'   (scale_regions mode) or reference point (reference_point mode). Default: 3000.
#'   
#' @param body Integer: For scale_regions mode only. The length in bins to which all
#'   regions will be scaled. Higher values give more detail but increase computation.
#'   Default: 3000. Ignored in reference_point mode.
#'   
#' @param bin_size Integer: Size of bins in base pairs for signal calculation. Smaller
#'   bins give higher resolution but increase computation time and memory usage.
#'   Default: 50.
#'   
#' @param reference_point Character: For reference_point mode only. Specifies which
#'   point to anchor on if input regions have width > 1 (though width == 1 is expected).
#'   Options: "center" (default), "TSS" (start), or "TES" (end).
#'   
#' @param sample_names Character vector: Custom names for each sample. If NULL, will
#'   use BigWig filenames with extensions removed. Length must match 'bw_files'.
#'   
#' @param region_names Character vector: Names for each region set in 'granges_list'.
#'   Used for row splitting in the heatmap. If NULL, will use "Region1", "Region2", etc.
#'   
#' @param colors Color mapping for heatmap. Can be:
#'   - NULL (default): Uses blue-white-red gradient from -2 to 2
#'   - A colorRamp2 object for custom color mapping
#'   - A character vector of colors (will be mapped linearly from -2 to 2)
#'   
#' @param kmeans Integer or NA: Number of k-means clusters to split rows. If provided,
#'   performs k-means clustering on the signal matrix and groups similar patterns together.
#'   Useful for identifying different regulatory patterns. Default: NULL (no clustering).
#'   
#' @param sort_by Character: How to sort rows within each group. Options:
#'   - "mean": Sort by mean signal across the region (default)
#'   - "median": Sort by median signal
#'   - "max": Sort by maximum signal
#'   - "none": No sorting (maintains input order)
#'   
#' @param scale Character: How to scale the signal values. Options:
#'   - "none" (default): No scaling, use raw signal values
#'   - "row": Z-score normalization per region (useful for comparing patterns)
#'   - "column": Z-score normalization per position (useful for comparing positions)
#'   
#' @param profile_height Unit object specifying the height of profile plot above each
#'   heatmap. Default: unit(2, "cm"). Increase for more prominent profile plots.
#'   
#' @param return_granges Logical: If TRUE, returns a list with both the heatmap object
#'   and the re-ordered GRanges that match the heatmap row order. If FALSE, returns
#'   only the heatmap. Default: FALSE. Useful for extracting regions in specific
#'   clusters or exporting sorted region lists.
#'   
#' @param ... Additional arguments passed to normalizeToMatrix() function from
#'   EnrichedHeatmap package. Advanced users can control smoothing, trimming, etc.
#'
#' @return 
#'   - If return_granges = FALSE: A HeatmapList object that can be drawn with draw()
#'   - If return_granges = TRUE: A list with elements:
#'     - $heatmap: The HeatmapList object
#'     - $granges: GRanges object with regions in the same order as heatmap rows
#'     - $clusters: If kmeans was used, cluster assignments for each region
#'
#' @examples
#' \dontrun{
#' # Example 1: scale_regions mode with gene bodies
#' # Visualize H3K36me3 across gene bodies (width > 1)
#' result <- heatmap_profile(
#'   bw_files = c("H3K36me3_rep1.bw", "H3K36me3_rep2.bw"),
#'   granges_list = genes,  # GRanges with actual gene coordinates
#'   mode = "scale_regions",
#'   upstream = 2000,
#'   downstream = 2000,
#'   body = 5000,
#'   sample_names = c("Replicate 1", "Replicate 2"),
#'   sort_by = "mean",
#'   kmeans = 3
#' )
#' draw(result)
#'
#' # Example 2: reference_point mode at TSS
#' # Visualize H3K4me3 around transcription start sites
#' tss <- resize(genes, width = 1, fix = "start")  # Create TSS points
#' result <- heatmap_profile(
#'   bw_files = c("H3K4me3.bw", "H3K27ac.bw"),
#'   granges_list = tss,
#'   mode = "reference_point",
#'   reference_point = "TSS",
#'   upstream = 3000,
#'   downstream = 3000,
#'   sample_names = c("H3K4me3", "H3K27ac"),
#'   return_granges = TRUE  # Get ordered gene list
#' )
#' draw(result$heatmap)
#' 
#' # Example 3: Multiple region sets with clustering
#' enhancers_center <- resize(enhancers, width = 1, fix = "center")
#' promoters_tss <- resize(promoters, width = 1, fix = "start")
#' result <- heatmap_profile(
#'   bw_files = "ATAC.bw",
#'   granges_list = list(enhancers_center, promoters_tss),
#'   mode = "reference_point",
#'   region_names = c("Enhancers", "Promoters"),
#'   kmeans = 4,
#'   return_granges = TRUE
#' )
#' # Export top cluster regions
#' top_cluster <- result$granges[result$clusters == 1]
#' }
#'
#' @export
# heatmap_profile <- function(
#     bw_files,
#     granges_list,
#     mode = c("scale_regions", "reference_point"),
#     upstream = 3000,
#     downstream = 3000,
#     body = 3000,
#     bin_size = 50,
#     reference_point = c("center", "tss", "tes"),
#     sample_names = NULL,
#     region_names = NULL,
#     color_scales = c(-2, 0, 2),
#     colors = c("#295072", "white", "#b12923"),
#     kmeans = NULL,
#     sort_by = c("mean", "median", "max", "none"),
#     scale = c("none", "row", "column"),
#     profile_height = unit(2, "cm"),
#     return_granges = FALSE,
#     ...
# ) {
#   
# 
#   # Input validation and setup
# 
#   
#   mode <- match.arg(mode)
#   reference_point <- match.arg(reference_point)
#   sort_by <- match.arg(sort_by)
#   scale <- match.arg(scale)
#   
#   # Convert single GRanges to list
#   if (!is.list(granges_list)) {
#     granges_list <- list(granges_list)
#   }
#   
#   # Set default sample names (remove file extensions)
#   if (is.null(sample_names)) {
#     sample_names <- basename(bw_files)
#     sample_names <- gsub("\\.bw$|\\.bigwig$|\\.bigWig$", "", sample_names)
#   }
#   
#   # Set default region names
#   if (is.null(region_names)) {
#     if (length(granges_list) == 1) {
#       region_names <- NULL  # Don't split if only one region set
#     } else {
#       region_names <- paste0("Region", seq_along(granges_list))
#     }
#   }
#   
#   # Setup color scheme
#   colors <- colorRamp2(color_scales, colors)
#   
# 
#   # Prepare genomic regions based on mode
# 
#   
#   # Combine all region sets into one GRanges object
#   all_regions <- do.call(c, granges_list)
#   
#   # Adjust regions based on mode
#   if (mode == "scale_regions") {
#     # For scale_regions: use regions as-is (assumes width > 1)
#     target <- all_regions
#   } else {
#     # For reference_point: ensure width == 1 at specified anchor
#     if (reference_point == "tss") {
#       target <- resize(all_regions, width = 1, fix = "start")
#     } else if (reference_point == "tes") {
#       target <- resize(all_regions, width = 1, fix = "end")
#     } else {
#       target <- resize(all_regions, width = 1, fix = "center")
#     }
#   }
#   
#   # Create region split factor for row annotation
#   if (!is.null(region_names)) {
#     region_lengths <- sapply(granges_list, length)
#     region_split <- rep(region_names, region_lengths)
#   } else {
#     region_split <- NULL
#   }
#   
# 
#   # Compute signal matrices for all samples
# 
#   
#   mat_list <- list()
#   
#   for (j in seq_along(bw_files)) {
#     message("Processing sample ", j, "/", length(bw_files), ": ", sample_names[j])
#     
#     # Import BigWig file
#     bw <- import(bw_files[j], format = "BigWig")
#     
#     # Compute normalized matrix based on mode
#     if (mode == "scale_regions") {
#       mat <- normalizeToMatrix(
#         bw, target,
#         value_column = "score",
#         extend = c(upstream, downstream),
#         mean_mode = "w0",
#         w = bin_size,
#         target_ratio = body / (upstream + downstream + body),
#         ...
#       )
#     } else {
#       mat <- normalizeToMatrix(
#         bw, target,
#         value_column = "score",
#         extend = c(upstream, downstream),
#         mean_mode = "w0",
#         w = bin_size,
#         ...
#       )
#     }
#     
#     # Apply scaling if requested
#     if (scale == "row") {
#       mat <- t(scale(t(mat)))
#     } else if (scale == "column") {
#       mat <- scale(mat)
#     }
#     
#     mat_list[[j]] <- mat
#   }
#   
#   names(mat_list) <- sample_names
#   
# 
#   # Row ordering and clustering
# 
#   
#   # Use first sample for ordering/clustering
#   reference_mat <- mat_list[[1]]
#   
#   # Calculate row order based on sorting method
#   if (sort_by != "none") {
#     row_order <- switch(
#       sort_by,
#       mean = order(rowMeans(reference_mat, na.rm = TRUE), decreasing = TRUE),
#       median = order(apply(reference_mat, 1, median, na.rm = TRUE), decreasing = TRUE),
#       max = order(apply(reference_mat, 1, max, na.rm = TRUE), decreasing = TRUE)
#     )
#   } else {
#     row_order <- seq_len(nrow(reference_mat))
#   }
#   
#   # Apply k-means clustering if requested
#   cluster_assignment <- NULL
#   if (!is.null(kmeans)) {
#     message("Performing k-means clustering with k=", kmeans)
#     set.seed(123)
#     km <- kmeans(reference_mat, centers = kmeans, nstart = 25)
#     row_order <- order(km$cluster)
#     cluster_assignment <- km$cluster[row_order]
#     cluster_split <- paste0("Cluster", cluster_assignment)
#   } else {
#     cluster_split <- region_split
#   }
#   
#   # Apply row ordering to all matrices
#   mat_list <- lapply(mat_list, function(m) m[row_order, , drop = FALSE])
#   
#   # Store ordered GRanges for output
#   ordered_granges <- all_regions[row_order]
#   
# 
#   # Create axis labels
#   if (mode == "scale_regions") {
#     axis_name <- c(
#       paste0("-", upstream/1000, "kb"),
#       "Start",
#       "End",
#       paste0("+", downstream/1000, "kb")
#     )
#   } else {
#     ref_label <- switch(reference_point,
#                         TSS = "TSS",
#                         TES = "TES",
#                         center = "Center")
#     axis_name <- c(
#       paste0("-", upstream/1000, "kb"),
#       ref_label,
#       paste0("+", downstream/1000, "kb")
#     )
#   }
#   
# 
#   # Calculate shared y-axis limits for profile plots
#   all_means <- lapply(mat_list, function(m) colMeans(m, na.rm = TRUE))
#   ylim_shared <- range(unlist(all_means), na.rm = TRUE)
#   print(ylim_shared)
#   ylim_shared <- c(0, max(ylim_shared))
#   print(ylim_shared)
#   
# 
#   # Prepare profile plot colors
#   if (!is.null(kmeans)) {
#     # For k-means: create colors for each cluster
#     cluster_colors <- rainbow(kmeans, v = 0.8, s = 0.7)
#     
#     # Create profile annotation function that colors by cluster
#     create_profile_anno <- function(mat, clusters) {
#       function(index) {
#         # Calculate mean for each cluster
#         cluster_profiles <- lapply(1:kmeans, function(k) {
#           cluster_idx <- which(clusters == k)
#           if (length(cluster_idx) > 0) {
#             colMeans(mat[cluster_idx, , drop = FALSE], na.rm = TRUE)
#           } else {
#             rep(NA, ncol(mat))
#           }
#         })
#         
#         # Create plot
#         pushViewport(viewport())
#         grid.rect(gp = gpar(fill = "white", col = NA))
#         
#         # Plot each cluster profile
#         for (k in 1:kmeans) {
#           profile <- cluster_profiles[[k]]
#           if (!all(is.na(profile))) {
#             x_coords <- seq(0, 1, length.out = length(profile))
#             y_coords <- (profile - ylim_shared[1]) / (ylim_shared[2] - ylim_shared[1])
#             grid.lines(x_coords, y_coords, gp = gpar(col = cluster_colors[k], lwd = 2))
#           }
#         }
#         
#         # Add axis
#         grid.yaxis(at = c(0, 0.5, 1), 
#                    label = round(c(ylim_shared[1], mean(ylim_shared), ylim_shared[2]), 2),
#                    gp = gpar(fontsize = 8))
#         
#         popViewport()
#       }
#     }
#   } else {
#     # For no k-means: use single color for all samples
#     profile_color <- "#E41A1C"  # Consistent red color
#   }
#   
# 
#   # Create heatmaps for all samples
# 
#   
#   ht_list <- NULL
#   
#   for (j in seq_along(mat_list)) {
#     
#     # Create top annotation based on whether k-means is used
#     if (!is.null(kmeans)) {
#       top_anno <- HeatmapAnnotation(
#         profile = create_profile_anno(mat_list[[j]], cluster_assignment),
#         height = profile_height,
#         annotation_name_side = "left"
#       )
#     } else {
#       top_anno <- HeatmapAnnotation(
#         profile = anno_enriched(
#           gp = gpar(col = profile_color, lwd = 2),
#           ylim = ylim_shared,
#           height = profile_height
#         )
#       )
#     }
#     
#     ht <- EnrichedHeatmap(
#       mat_list[[j]],
#       col = colors,
#       name = if (j == 1) "Signal" else paste0("Signal_", j),
#       column_title = sample_names[j],
#       axis_name = axis_name,
#       axis_name_rot = 0,
#       border = TRUE,
#       cluster_rows = FALSE,
#       show_row_names = FALSE,
#       show_heatmap_legend = (j == 1),
#       use_raster = TRUE,
#       raster_quality = 2,
#       split = cluster_split,
#       top_annotation = top_anno,
#       width = unit(4, "cm")
#     )
#     
#     # Concatenate heatmaps
#     if (is.null(ht_list)) {
#       ht_list <- ht
#     } else {
#       ht_list <- ht_list + ht
#     }
#   }
#   
# 
#   # Return results
# 
#   
#   if (return_granges) {
#     result <- list(
#       heatmap = ht_list,
#       granges = ordered_granges
#     )
#     if (!is.null(cluster_assignment)) {
#       result$clusters <- cluster_assignment
#     }
#     return(result)
#   } else {
#     return(ht_list)
#   }
# }
heatmap_profile <- function(
    bw_files,
    granges_list,
    mode = c("scale_regions", "reference_point"),
    upstream = 3000,
    downstream = 3000,
    body = 3000,
    bin_size = 50,
    reference_point = c("center", "tss", "tes"),
    sample_names = NULL,
    region_names = NULL,
    color_scales = c(-2, 0, 2), 
    colors = c("#295072", "white", "#b12923"), 
    kmeans = NA,
    sort_by = c("mean", "median", "max", "none"),
    scale = c("none", "row", "column"),
    profile_height = unit(2, "cm"),
    return_granges = FALSE,
    ...
) {
  # 1. 基础验证 ----------------------------------------------------------------
  mode <- match.arg(mode)
  reference_point <- match.arg(reference_point)
  sort_by <- match.arg(sort_by)
  scale <- match.arg(scale)
  
  if (length(color_scales) != length(colors)) stop("color_scales 与 colors 长度必须一致")
  if (!is.list(granges_list)) granges_list <- list(granges_list)
  if (is.null(sample_names)) {
    sample_names <- gsub("\\.bw$|\\.bigwig$", "", basename(bw_files), ignore.case = TRUE)
  }
  
  # 2. 准备区域 ----------------------------------------------------------------
  all_regions <- do.call(c, granges_list)
  fix_map <- c("tss" = "start", "tes" = "end", "center" = "center")
  target <- if (mode == "scale_regions") all_regions else IRanges::resize(all_regions, width = 1, fix = fix_map[reference_point])
  
  # 3. 计算矩阵 ----------------------------------------------------------------
  mat_list <- lapply(seq_along(bw_files), function(j) {
    message("Processing: ", sample_names[j])
    bw <- rtracklayer::import(bw_files[j], format = "BigWig")
    mat <- normalizeToMatrix(
      bw, target, value_column = "score", extend = c(upstream, downstream),
      mean_mode = "w0", w = bin_size
      # target_ratio = if(mode == "scale_regions") body / (upstream + downstream + body) else 0.4#,
      # ...
    )
    if (scale == "row") mat <- t(scale(t(mat)))
    if (scale == "column") mat <- scale(mat)
    return(mat)
  })
  
  # 4. 排序与分群 --------------------------------------------------------------
  ref_mat <- mat_list[[1]]
  row_order <- seq_len(nrow(ref_mat))
  cluster_split <- if(!is.null(region_names)) rep(region_names, sapply(granges_list, length)) else NULL
  cluster_assignment <- NULL
  
  if (sort_by != "none") {
    val <- switch(sort_by, mean = rowMeans(ref_mat, na.rm = TRUE),
                  median = apply(ref_mat, 1, median, na.rm = TRUE),
                  max = apply(ref_mat, 1, max, na.rm = TRUE))
    row_order <- order(val, decreasing = TRUE)
  }
  
  if (!is.na(kmeans)) {
    set.seed(42)
    km <- kmeans(ref_mat, centers = kmeans, nstart = 25)
    cluster_assignment <- km$cluster
    row_order <- order(cluster_assignment, rowMeans(ref_mat, na.rm = TRUE), decreasing = c(FALSE, TRUE))
    cluster_split <- paste0("Cluster ", cluster_assignment[row_order])
    cluster_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                        "#FFFF33", "#A65628", "#F781BF", "#999999")[1:kmeans]
  } else {
    cluster_colors <- "firebrick" # Default color if no kmeans
  }
  
  mat_list <- lapply(mat_list, function(m) m[row_order, , drop = FALSE])
  
  # 5. 修复 Y 轴范围 (解决 unit 长度为 0 的错误) -------------------------------
  all_vals <- unlist(lapply(mat_list, function(x) colMeans(x, na.rm=T)))
  min_limit <- quantile(all_vals, 0.01, na.rm = TRUE)
  max_limit <- quantile(all_vals, 0.99, na.rm = TRUE)
  # all_means <- unlist(lapply(mat_list, function(m) colMeans(m, na.rm = TRUE)))
  # ylim_shared <- c(min_limit, max(max_limit, max(all_means, na.rm = TRUE)))
  ylim_shared <- c(min_limit - abs(min_limit) * 0.1, max_limit * 1.2)
  print(ylim_shared)
  
  if (is.null(color_scales[1]) | is.na(is.null(color_scales[1]))) 
    color_scales[1] <- quantile(all_vals, 0.01, na.rm = TRUE)
  if (is.null(color_scales[length(color_scales)]) | is.na(color_scales[length(color_scales)]))
    color_scales[length(color_scales)] <- quantile(all_vals, 0.99, na.rm = TRUE)
  
  col_fun <- circlize::colorRamp2(color_scales, colors)
  axis_labels <- if(mode == "scale_regions") {
    c(paste0("-", upstream/1000, "kb"), "Start", "End", paste0("+", downstream/1000, "kb"))
  } else {
    c(paste0("-", upstream/1000, "kb"), toupper(reference_point), paste0("+", downstream/1000, "kb"))
  }
  
  # 6. 绘图 --------------------------------------------------------------------
  ht_list <- NULL
  for (j in seq_along(mat_list)) {
    
    top_anno <- HeatmapAnnotation(
      enriched = anno_enriched(
        gp = gpar(col = cluster_colors, lwd = 2),
        # gp = if(!is.na(kmeans)) gpar(col = 1:kmeans) else gpar(col = "firebrick"),
        ylim = ylim_shared,
        height = profile_height,
        # 修复点 1: 仅第一个样本显示 Y 轴刻度和标签
        axis_param = list(
          side = "left",
          gp = gpar(fontsize = 6)
        )
      ),
      show_annotation_name = (j == 1),
      annotation_name_side = "left"
    )

    # 修复点 2: 共享 Legend 名，仅第一个显示
    ht <- EnrichedHeatmap(
      mat_list[[j]],
      col = col_fun,
      name = if (j > 1) NULL else "Signal", 
      column_title = sample_names[j],
      axis_name = axis_labels,
      axis_name_rot = 0,
      split = cluster_split,
      cluster_rows = FALSE,
      top_annotation = top_anno,
      show_row_names = FALSE,
      use_raster = TRUE,
      show_heatmap_legend = (j == 1), 
      width = unit(4, "cm"),
      height = unit(8, "cm")
    )
    ht_list <- if (is.null(ht_list)) ht else ht_list + ht
  }
  
  # 6. Create Cluster Legend
  annotation_legends <- list()
  if (!is.na(kmeans)) {
    annotation_legends[[1]] <- Legend(
      labels = paste0("Cluster ", 1:kmeans),
      title = "Clusters",
      type = "lines",
      legend_gp = gpar(col = cluster_colors, lwd = 4)
    )
  }
  
  p <- draw(ht_list, 
            heatmap_legend_side = "right", 
            annotation_legend_list = annotation_legends)
  if (return_granges) {
    return(list(heatmap = p, granges = all_regions[row_order], clusters = cluster_assignment[row_order]))
  } else {
    draw(ht_list, 
         heatmap_legend_side = "right", 
         annotation_legend_list = annotation_legends)
  }
}


