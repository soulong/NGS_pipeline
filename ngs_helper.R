
library(GenomicRanges)
# library(consensusSeekeR)
library(EnrichedHeatmap)
library(rtracklayer)
library(circlize)
library(ComplexHeatmap)
library(grid)




#' @Title Detect Operating System
#' @description Identify the current OS (windows, wsl, or linux)
#' @return Character string: "windows", "wsl", "linux", or system name
detect_os <- function() {
  sysinfo <- Sys.info()
  
  if (grepl("microsoft|Microsoft", sysinfo["release"], ignore.case = TRUE)) {
    return("wsl")
  }
  
  if (.Platform$OS.type == "windows") {
    return("windows")
  }
  
  if (sysinfo["sysname"] == "Linux") {
    return("linux")
  }
  
  return(sysinfo["sysname"])
}


#' @Title Safe Path Normalization
#' @description Normalize file paths and convert Windows paths to WSL format when needed
#' @param path Character; the file path to process
#' @param standardize Logical; whether to convert Windows paths to WSL format (default: TRUE)
#' @param os Character; the operating system (default: auto-detect)
#' @return Character string with normalized path
norm_path <- function(path=r"(C:\Users\haohe\GitHub\shinyLab)", 
                      standardize=TRUE, 
                      os=detect_os()
) {
  path <- path |> 
    gsub("\\\\", "/", x = _) |> 
    gsub("/+", "/", x = _) |> 
    trimws()
  
  if(standardize & grepl("^[A-Za-z]:[/\\\\]", path) & os == "wsl") {
    path <- file.path(
      "/mnt", 
      tolower(sub("^([A-Za-z]):.*", "\\1", path)), 
      gsub("\\\\", "/", sub("^[A-Za-z]:[/\\\\]", "", path)))
  }
  
  return(path)
}




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
compute_consensus <- function(
    peak_list,
    groups = NULL,
    method = c("granges", "consensusseeker"),
    genome_build = NULL,
    min_occurrence = 2L,
    max_coverage   = Inf,
    min_gapwidth   = 1L,
    ...) {
  
  # ── Input validation ────────────────────────────────────────────────────────
  if (!inherits(peak_list, c("list", "CompressedGRangesList")) || length(peak_list) == 0) {
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
  message(sprintf("using method = '%s' ...", method))
  
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
  
  return(consensus_by_group)
} 


#' Compute read counts for genomic regions and their flanking regions
#'
#' This function calculates read counts and coverage for genomic regions
#' and optionally their left/right flanking regions from BAM files.
#'
#' @param bam_path Path to BAM file (with .bai index)
#' @param gr GRanges object with genomic regions
#' @param flank_size Size of flanking regions to include (default: 0)
#' @param mode Character vector specifying what to return: "count" for raw counts,
#'        "coverage" for normalized coverage, or c("count", "coverage") for both
#' @param ... Additional arguments passed to ScanBamParam
#' @return List with two data frames: $count and $coverage. Each data frame contains
#'         columns: seqnames, start, end, width, name (optional),
#'         count_region, count_left, count_right (for count) or
#'         cov_region, cov_left, cov_right (for coverage)
#' @export
compute_count <- function(
    bam_path,
    gr,
    flank_size = 0,
    mode = c("count", "coverage"),
    ...
) {
  
  mode <- match.arg(mode, choices = c("count", "coverage"), several.ok = TRUE)
  
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
  
  # Key change: do not drop any flank regions, even if width <= 0
  gr_region <- gr
  
  # Left flank: allows start < 1
  if(flank_size > 0) {
    gr_left   <- flank(gr, width = flank_size, start = TRUE, both = FALSE)
    
    # Right flank: allows end > seqlengths
    gr_right  <- flank(gr, width = flank_size, start = FALSE, both = FALSE)
    
    # Combine all regions to query (but still count separately)
    all_regions <- c(gr_region, gr_left, gr_right)
  } else {
    all_regions <- gr_region
  }
  
  param <- ScanBamParam(
    what = c("rname", "strand", "pos", "qwidth"),
    which = IRanges::reduce(all_regions, ignore.strand = TRUE),  # Merge overlapping regions for faster reading
    ...
  )
  
  message("Reading alignments...")
  gal <- readGAlignmentPairs(bam_path, param = param)
  
  message("Computing overlaps...")
  
  # Count separately (even if width <= 0, findOverlaps will return 0)
  count_region <- countOverlaps(gr_region, gal, ignore.strand = TRUE)
  
  # Initialize count df with base columns
  df_count <- df
  df_count$count_region <- count_region
  
  # Initialize coverage df with base columns
  df_coverage <- df
  df_coverage$cov_region <- count_region / (width(gr_region) / 1000)
  
  if(flank_size > 0) {
    count_left   <- countOverlaps(gr_left,   gal, ignore.strand = TRUE)
    df_count$count_left   <- count_left
    
    count_right  <- countOverlaps(gr_right,  gal, ignore.strand = TRUE)
    df_count$count_right  <- count_right
    
    # Coverage calculations
    df_coverage$cov_left   <- count_left   / (flank_size / 1000)
    df_coverage$cov_right  <- count_right  / (flank_size / 1000)
  } else {
    # Add NA columns for consistency when no flanks
    df_count$count_left <- NA_integer_
    df_count$count_right <- NA_integer_
    df_coverage$cov_left <- NA_real_
    df_coverage$cov_right <- NA_real_
  }
  
  # Build result based on mode
  result <- list()
  
  if ("count" %in% mode) {
    result$count <- df_count
  } else {
    result$count <- data.frame()
  }
  
  if ("coverage" %in% mode) {
    result$coverage <- df_coverage
  } else {
    result$coverage <- data.frame()
  }
  
  message("Done")
  return(result)
}




#' Compute signal matrices from bigWig files or GRangesList
#'
#' This function imports signal data from bigWig/BigWig files or GRangesList objects
#' and computes normalized matrices over genomic regions using EnrichedHeatmap's
#' normalizeToMatrix function.
#'
#' @param signal Named character vector of bigWig/bigwig file paths, or
#'        named GRangesList/CompressedGRangesList object. Names are used as sample names.
#' @param region Genomic regions (GRanges object). Required when signal is file paths or GRangesList.
#' @param mode "scale_regions" or "reference_point"
#' @param reference_point "center", "tss", or "tes"
#' @param upstream Upstream extension in bp (default: 3000)
#' @param downstream Downstream extension in bp (default: 3000)
#' @param bin_size Bin size for matrix computation (default: 50)
#' @param scale Scaling method: "none", "row", or "column"
#' @return Named list of matrices (names derived from signal names)
#' @export
compute_signal_matrix <- function(
    signal,
    region,
    mode = c("scale_regions", "reference_point"),
    reference_point = c("center", "tss", "tes"),
    upstream = 3000,
    downstream = 3000,
    bin_size = 50,
    scale = c("none", "row", "column"),
    ...
) {
  # Input validation
  mode <- match.arg(mode)
  reference_point <- match.arg(reference_point)
  scale <- match.arg(scale)
  
  # Validate signal: must be named file paths or named GRangesList
  if (is.character(signal)) {
    # Named file paths
    if (is.null(names(signal)) || any(names(signal) == "")) {
      stop("signal must be named character vector of file paths")
    }
    if (!all(file.exists(signal))) {
      stop("Some signal files do not exist: ", paste(signal[!file.exists(signal)], collapse = ", "))
    }
    signal_input_type <- "files"
  } else if (inherits(signal, "CompressedGRangesList") || inherits(signal, "GRangesList")) {
    # Named GRangesList
    if (is.null(names(signal)) || any(names(signal) == "")) {
      stop("signal must be named GRangesList")
    }
    signal_input_type <- "granges_list"
  } else {
    stop("signal must be: (1) named character vector of bw/bigwig file paths, ",
         "or (2) named GRangesList/CompressedGRangesList")
  }
  
  # Sample names come from names(signal)
  sample_names <- names(signal)
  
  # Validate region: must be GRanges
  if (is.null(region)) {
    stop("region must be provided as GRanges")
  }
  
  if (!is(region, "GRanges")) {
    stop("region must be a GRanges object")
  }
  
  fix_map <- c("tss" = "start", "tes" = "end", "center" = "center")
  
  # Prepare target based on mode
  target <- if (mode == "scale_regions") {
    region
  } else {
    IRanges::resize(region, width = 1, fix = fix_map[reference_point])
  }
  
  
  
  # Compute matrices
  
  mat_list <- NULL
  
  if (signal_input_type == "files") {
    # Import from bigWig files
    mat_list <- lapply(seq_along(signal), function(j) {
      message("Importing: ", sample_names[j], " from ", basename(signal[j]))
      bw <- rtracklayer::import(signal[j])
      mat <- normalizeToMatrix(
        bw, target,
        value_column = "score",
        extend = c(upstream, downstream),
        mean_mode = "w0",
        w = bin_size
      )
      return(mat)
    })
    names(mat_list) <- sample_names
    
  } else if (signal_input_type == "granges_list") {
    # Import from GRangesList
    mat_list <- lapply(seq_along(signal), function(j) {
      message("Processing: ", sample_names[j])
      bw <- signal[[j]]
      mat <- normalizeToMatrix(
        bw, target,
        value_column = "score",
        extend = c(upstream, downstream),
        mean_mode = "w0",
        w = bin_size,
        ...
      )
      return(mat)
    })
    names(mat_list) <- sample_names
  }
  
  # Apply scaling across ALL matrices using global mean/sd
  if (scale == "row") {
    # Combine all matrices to compute global row statistics
    all_values <- do.call(rbind, lapply(mat_list, function(m) as.vector(m)))
    global_mean <- mean(all_values, na.rm = TRUE)
    global_sd <- sd(all_values, na.rm = TRUE)
    mat_list <- lapply(mat_list, function(mat) (mat - global_mean) / global_sd)
  } else if (scale == "column") {
    # Combine all matrices to compute global column statistics
    all_mat <- do.call(cbind, mat_list)
    col_means <- apply(all_mat, 2, mean, na.rm = TRUE)
    col_sds <- apply(all_mat, 2, sd, na.rm = TRUE)
    mat_list <- lapply(mat_list, function(mat) {
      for (j in seq_len(ncol(mat))) {
        mat[, j] <- (mat[, j] - col_means[j]) / col_sds[j]
      }
      return(mat)
    })
  }
  
  # Return named matrix list
  return(mat_list)
}


#' Generate enriched heatmap profiles for genomic signal data
#'
#' This function creates publication-ready heatmaps showing signal enrichment
#' over genomic regions using EnrichedHeatmap. Supports multiple signal inputs,
#' region grouping, k-means clustering, and various visualization options.
#'
#' @param signal Named character vector of bigWig/bigwig file paths, or
#'        named list of pre-computed matrices. Names are used as sample names.
#' @param region Genomic regions (GRanges object). Required when signal is file paths.
#'        Not required when signal is named matrix list.
#' @param mode "scale_regions" to scale all regions to same length, or "reference_point"
#'        to anchor at a reference point (default: "reference_point")
#' @param reference_point Reference point for alignment: "center", "tss", or "tes"
#'        (default: "center")
#' @param upstream Length of upstream region in bp (default: 3000)
#' @param downstream Length of downstream region in bp (default: 3000)
#' @param body Length of body region for scale_regions mode in bp (default: 3000)
#' @param bin_size Bin size for matrix computation (default: 50)
#' @param color_scales Numeric vector for color scale breakpoints (default: c(NA, 0, NA))
#' @param colors Character vector of colors for color scale (default: c("#295072", "white", "#b12923"))
#' @param kmeans Number of clusters for k-means clustering (default: NA, no clustering)
#' @param sort_by_sample Sample index to use for sorting/clustering (default: 1)
#' @param sort_by Method for sorting rows: "mean", "median", "max", or "none" (default: "mean")
#' @param scale Scaling method: "none", "row", or "column" (default: "none")
#' @param profile_height Height of the top enrichment profile annotation (default: unit(2, "cm"))
#' @param return_granges Whether to return granges along with the plot (default: FALSE)
#' @param ... Additional arguments passed to EnrichedHeatmap
#' @return Returns list with heatmap, granges, and clusters.
#' @export
heatmap_profile <- function(
    signal,
    region = NULL,
    mode = c("reference_point", "scale_regions"),
    reference_point = c("center", "tss", "tes"),
    upstream = 3000,
    downstream = 3000,
    body = 3000,
    bin_size = 50,
    color_scales = c(NA, 0, NA), 
    colors = c("#295072", "white", "#b12923"), 
    kmeans = NA,
    sort_by_sample = 1,
    sort_by = c("mean", "median", "max", "none"),
    scale = c("none", "row", "column"),
    profile_height = unit(2, "cm"),
    ...
) {
  
  # Input Format Checking
  
  mode <- match.arg(mode)
  reference_point <- match.arg(reference_point)
  sort_by <- match.arg(sort_by)
  scale <- match.arg(scale)
  
  # Check signal input type
  signal_input_type <- NULL
  if (is.character(signal)) {
    # Named file paths to bw/bigwig files
    if (is.null(names(signal)) || any(names(signal) == "")) {
      stop("signal must be named character vector of file paths")
    }
    if (!all(file.exists(signal))) {
      stop("Some signal files do not exist: ", paste(signal[!file.exists(signal)], collapse = ", "))
    }
    signal_input_type <- "files"
  } else if (is.list(signal) && !is.null(names(signal))) {
    # Named list of pre-computed matrices
    first_elem <- signal[[1]]
    if (is.matrix(first_elem) || inherits(first_elem, "matrix")) {
      signal_input_type <- "mat_list"
    } else {
      stop("Unrecognized list type in signal. Must be named list of matrices")
    }
  } else {
    stop("signal must be: (1) named character vector of bw/bigwig file paths, ",
         "or (2) named list of matrices")
  }
  
  # Check region input type
  if (is.null(region)) {
    # If region is NULL, require signal to be pre-computed matrices
    if (signal_input_type == "files") {
      stop("region must be provided when signal is file paths")
    }
  } else {
    # Validate region input - must be GRanges
    if (!is(region, "GRanges")) {
      stop("region must be a GRanges object")
    }
  }
  
  # Color scales validation
  if (length(color_scales) != length(colors)) {
    stop("color_scales and colors must have the same length")
  }
  
  
  # Set up sample names from signal
  
  sample_names <- names(signal)
  
  
  # Step 1-3: Compute matrices if needed (skip if pre-computed)
  
  if (signal_input_type == "mat_list") {
    # Use pre-computed matrices directly - start from step 4
    mat_list <- signal
    # Ensure names
    if (is.null(names(mat_list))) {
      names(mat_list) <- sample_names
    }
    # Skip to step 4 - get info from pre-computed
    mode <- "reference_point"  # Assume reference_point for pre-computed
    all_regions <- NULL
  } else {
    # Need to compute matrices from files
    # Use separate function for computing signal matrices
    mat_list <- compute_signal_matrix(
      signal = signal,
      region = region,
      mode = mode,
      reference_point = reference_point,
      upstream = upstream,
      downstream = downstream,
      bin_size = bin_size,
      scale = scale
    )
    
    # Update all_regions from the computed result
    all_regions <- region
  }
  
  
  # Step 4: Sorting and Clustering
  
  # Validate sort_by_sample index
  if (sort_by_sample < 1 || sort_by_sample > length(mat_list)) {
    warning("sort_by_sample (", sort_by_sample, ") is out of range, using 1")
    sort_by_sample <- 1
  }
  
  ref_mat <- mat_list[[sort_by_sample]]
  ref_sample_name <- names(mat_list)[sort_by_sample]
  message("Using sample '", ref_sample_name, "' (index ", sort_by_sample, ") for sorting/clustering")
  
  row_order <- seq_len(nrow(ref_mat))
  
  # Create cluster_split - no region grouping in new design
  if (!is.null(all_regions)) {
    cluster_split <- rep("Region", length(all_regions))
  } else {
    cluster_split <- NULL
  }
  cluster_assignment <- NULL
  
  if (sort_by != "none") {
    val <- switch(sort_by, 
                  mean = rowMeans(ref_mat, na.rm = TRUE),
                  median = apply(ref_mat, 1, median, na.rm = TRUE),
                  max = apply(ref_mat, 1, max, na.rm = TRUE))
    row_order <- order(val, decreasing = TRUE)
  }
  
  if (!is.na(kmeans)) {
    # Remove rows with NA/NaN/Inf values before kmeans
    valid_rows <- apply(ref_mat, 1, function(x) all(is.finite(x)))
    if (!all(valid_rows)) {
      message("Removing ", sum(!valid_rows), " rows with NA/NaN/Inf values before kmeans clustering")
      ref_mat_clean <- ref_mat[valid_rows, , drop = FALSE]
      valid_indices <- which(valid_rows)
    } else {
      ref_mat_clean <- ref_mat
      valid_indices <- seq_len(nrow(ref_mat))
    }
    set.seed(42)
    km <- kmeans(ref_mat_clean, centers = kmeans, nstart = 25)
    cluster_assignment <- km$cluster
    row_order <- valid_indices[order(cluster_assignment, rowMeans(ref_mat_clean, na.rm = TRUE), decreasing = c(FALSE, TRUE))]
    cluster_split <- paste0("Cluster ", cluster_assignment[order(cluster_assignment, rowMeans(ref_mat_clean, na.rm = TRUE), decreasing = c(FALSE, TRUE))])
    cluster_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                        "#FFFF33", "#A65628", "#F781BF", "#999999")[1:kmeans]
  } else {
    cluster_colors <- "firebrick"
    cluster_assignment <- NULL
  }
  
  mat_list <- lapply(mat_list, function(m) m[row_order, , drop = FALSE])
  
  
  # Step 5: Y-axis range and color scale
  # Bug fix 1: Calculate ylim per cluster when kmeans is enabled
  if (!is.na(kmeans) && !is.null(cluster_assignment)) {
    # Calculate ylim for each cluster separately
    cluster_ylims <- lapply(unique(cluster_assignment), function(k) {
      cluster_rows <- which(cluster_assignment == k)
      cluster_mats <- lapply(mat_list, function(m) m[cluster_rows, , drop = FALSE])
      cluster_vals <- unlist(lapply(cluster_mats, function(x) colMeans(x, na.rm = TRUE)))
      c(quantile(cluster_vals, 0.01, na.rm = TRUE), quantile(cluster_vals, 0.99, na.rm = TRUE))
    })
    # Take global min/max across all clusters
    global_min <- min(sapply(cluster_ylims, `[`, 1))
    global_max <- max(sapply(cluster_ylims, `[`, 2))
    ylim_shared <- c(global_min - abs(global_min) * 0.1, global_max * 1.2)
  } else {
    all_vals <- unlist(lapply(mat_list, function(x) colMeans(x, na.rm = TRUE)))
    min_limit <- quantile(all_vals, 0.01, na.rm = TRUE)
    max_limit <- quantile(all_vals, 0.99, na.rm = TRUE)
    ylim_shared <- c(min_limit - abs(min_limit) * 0.1, max_limit * 1.2)
  }
  print(ylim_shared)
  
  if (is.null(color_scales[1]) || is.na(color_scales[1])) {
    color_scales[1] <- quantile(all_vals, 0.01, na.rm = TRUE)
  }
  if (is.null(color_scales[length(color_scales)]) || is.na(color_scales[length(color_scales)])) {
    color_scales[length(color_scales)] <- quantile(all_vals, 0.99, na.rm = TRUE)
  }
  
  col_fun <- circlize::colorRamp2(color_scales, colors)
  axis_labels <- if (mode == "scale_regions") {
    c(paste0("-", upstream/1000, "kb"), "Start", "End", paste0("+", downstream/1000, "kb"))
  } else {
    c(paste0("-", upstream/1000, "kb"), str_to_sentence(reference_point), paste0("+", downstream/1000, "kb"))
  }
  
  
  # Step 6: Plotting
  ht_list <- NULL
  for (j in seq_along(mat_list)) {
    
    top_anno <- HeatmapAnnotation(
      enriched = anno_enriched(
        gp = gpar(col = cluster_colors, lwd = 2),
        ylim = ylim_shared,
        axis_param = list(
          at = c(ylim_shared[1], 0, ylim_shared[2]),
          labels = c(round(ylim_shared[1], 2), "0", round(ylim_shared[2], 2)),
          side = "left",
          gp = gpar(fontsize = 6))),
      height = profile_height,
      show_annotation_name = (j == 1),
      annotation_name_side = "left"
    )
    
    ht <- EnrichedHeatmap(
      mat_list[[j]],
      col = col_fun,
      name = if (j > 1) NULL else "Signal",
      column_title = names(mat_list)[j],
      axis_name = axis_labels,
      axis_name_rot = 0,
      split = cluster_split,
      cluster_rows = FALSE,
      top_annotation = top_anno,
      show_row_names = FALSE,
      use_raster = TRUE,
      show_heatmap_legend = (j == 1),
      heatmap_legend_param = list(
        legend_direction = "horizontal",
        at = c(ylim_shared[1], ylim_shared[2]),
        labels = c(round(ylim_shared[1], 2), round(ylim_shared[2], 2)))
    )
    ht_list <- if (is.null(ht_list)) ht else ht_list + ht
  }
  
  # Create Cluster Legend
  # Bug fix 2: Add direction = "horizontal" to legend
  annotation_legends <- list()
  if (!is.na(kmeans)) {
    annotation_legends[[1]] <- Legend(
      labels = paste0("Cluster ", 1:kmeans),
      title = "Clusters",
      type = "lines",
      direction = "horizontal",
      legend_gp = gpar(col = cluster_colors, lwd = 4)
    )
  }
  
  p <- draw(ht_list,
            annotation_legend_list = annotation_legends,
            heatmap_legend_side = "bottom",
            annotation_legend_side = "bottom",
            padding = unit(c(10, 10, 10, 10), "mm"),
            use_raster = TRUE)
  
  return(list(heatmap = p,
              granges = if (!is.null(all_regions)) all_regions[row_order] else NULL,
              clusters = cluster_assignment[row_order]
  )
  )
}

