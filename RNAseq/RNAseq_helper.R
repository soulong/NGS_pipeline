
library(tidyverse)
library(furrr)
library(clusterProfiler)
library(ComplexHeatmap)

#' @param df a data.frame
#' @param group named vector, names are columns in df, values are groups
#' @param threshold min threshold to keep
filter_df_by_group <- function(df, group, threshold=5) {
  filtered_index <- as.character(unique(group)) %>% 
    map(\(x) group[group==x] %>%
          names() %>% 
          df[, ., drop=FALSE] %>% 
          rowMeans(na.rm=T) %>% 
          # at least one group has mean counts > X
          magrittr::is_greater_than(threshold) %>% 
          which()) %>% #lengths()
    unlist() %>%
    unique()
  filtered <- df[filtered_index, ]
  print(str_glue("before filter: {nrow(df)}\n after filter: {nrow(filtered)}"))
  return(filtered)
}



#' @param x a numeric vector
se <- function(x) {sd(x)/sqrt(length(x))}



#' @param mat a matrix, a data.frame will be convert to matrix
#' @param by scale by which order
scale_mat <- function(mat, by=c("row","column","none")) {
  
  if(is.data.frame(mat)) {
    print("The input is a data frame, convert it to the matrix")
    mat = as.matrix(mat) }
  
  if(by[1]=="row") {
    if(any(is.na(mat))) {
      mat = (mat - rowMeans(mat, na.rm=T))/rowSds(mat, na.rm=T)
    } else mat = t(scale(t(mat)))
  }
  
  if(by[1]=="column") {
    if(any(is.na(mat))) {
      mat = t((t(mat) - colMeans(mat, na.rm=T))/colSds(mat, na.rm=T))
    } else mat = scale(mat)
  }
  
  if(by[1]=="none") {}
  
  return(mat)
}


#' @description
#' ComplexHeatmap column annotation from df
#' 
#' @param df data.frame, first column must be sample names
#' @param split_by df column, character, split heatmap by sample group
#' @param annotate_col df columns, character vector, only works if split_by is NULL
annoCol_from_sample_info <- function(
    sample_df, split_by=NULL, annotate_col=NULL) {
  
  sample_df <- select(sample_df, 1, any_of(c(split_by, annotate_col)))
  
  if(!is.null(split_by)) {
    # split heatmap column by sample group
    if(!is.factor(sample_df[[split_by]][1])) {
      print("convert contrast to factor")
      sample_df[[split_by]] <- sample_df[[split_by]] %>% 
        factor(levels=unique(.[[split_by]])) }
    #
    sample_split <- sample_df[["contrast"]]
    sample_split_color <- 
      colorRampPalette(brewer.pal(8,"Dark2"))(length(unique(sample_split))) %>% 
      set_names(unique(sample_split))
    #
    anno_sample <- HeatmapAnnotation(contrast=anno_block(
      labels=unique(sample_split), labels_gp=gpar(fontsize=5),
      gp=gpar(fill=sample_split_color)), which="column")
    
    return(list(anno_sample=anno_sample, 
                sample_split=sample_split, 
                sample_split_color=sample_split_color))
  } else {
    # anno for each sample
    sample_anno <- sample_df[,-1,drop=F]
    sample_anno_color <- map(
      colnames(sample_anno), \(x) sample_anno[[x]] %>% 
        {colorRampPalette(brewer.pal(8,"Dark2"))(length(unique(.)))[match(., unique(.))]} %>% 
        set_names(sample_anno[[x]]) ) %>% set_names(colnames(sample_anno))
    #
    anno_sample <- HeatmapAnnotation(
      df=sample_anno, which="column",
      col=sample_anno_color, annotation_name_side="left")
    
    return(list(anno_sample=anno_sample, sample_anno_color=sample_anno_color))
  }
}

#' @description
#' ComplexHeatmap row annotation for cluster
#' 
#' @param mat expr mat used for cluster
#' @param k kmean cluster k
annoRow_from_cluster <- function(mat, k=NULL) {
  set.seed(42)
  cl <- kmeans(mat, k)
  # cl$cluster <- str_c("C", cl$cluster)
  cl_uni <- unique(cl$cluster) %>% sort()
  cl_color <- colorRampPalette(brewer.pal(12,"Paired"))(length(cl_uni))
  # anno cluster
  align_to = split(1:nrow(mat), cl$cluster)
  anno_cl <- anno_block(
    align_to = align_to, 
    panel_fun = function(index, nm) {
      npos = as.numeric(nm)
      # npos = as.numeric(unlist(strsplit(nm, split="C"))[2])
      grid.rect(gp=gpar(fill=cl_color[npos], col=NA))
      grid.text(label=paste("C",nm,":",length(index), sep=""), 
                rot=90, gp=gpar(col="grey20", fontsize=6))
    }, which = "row")
  
  return(list(anno_cl=anno_cl, align_to=align_to, cl=cl, cl_color=cl_color))
}



#' @param align_to named vector, come from annoRow_from_cluster[["align_to"]]
#' @param gglist list of ggplot object, names should be same with heatmap cluter name
#' @param ggplot.panel.arg ggplot args
anno_gglist <- function(
    align_to, gglist, 
    # panel size, gap, width, fill, color
    ggplot.panel.arg = c(2, 0.25, 5, "grey90","black")) {
  anno_enrich <- anno_zoom(
    align_to = align_to, which = "row", 
    panel_fun = function(index, nm) {
      g <- gglist[[nm]]
      g <- grid.grabExpr(print(g))
      grid::pushViewport(viewport())
      grid::grid.rect()
      grid::grid.draw(g)
      grid::popViewport() },
    size = unit(as.numeric(ggplot.panel.arg[1]), "cm"), 
    gap = unit(as.numeric(ggplot.panel.arg[2]), "cm"), 
    width = unit(as.numeric(ggplot.panel.arg[3]),  "cm"),
    side = "right", 
    link_gp = gpar(fill=ggplot.panel.arg[4], col=ggplot.panel.arg[5]))
  
  return(anno_enrich)
}



#' @param gene_list list of character genes
#' @param geneset_list list of dataframe geneset
enrich_ora_list <- function(gene_list, geneset_list) {
  enrich_result <- list()
  for(i in seq_along(gene_list)) {
    print(str_glue(
      ">>>>>>>> ORA enrichment on: {names(gene_list)[i]}"))
    genes <- gene_list[[i]]
    
    if(length(genes) < 20) next
    ora_list <- future_map(
      geneset_list, 
      \(x) enricher(gene = genes, 
                    pvalueCutoff = 0.2, qvalueCutoff = 0.7, 
                    minGSSize = 10, maxGSSize = 500,
                    TERM2GENE = x$term2gene,
                    TERM2NAME = x$term2name),
      .options = furrr_options(seed = T)
    )
    # merge different geneset_list 
    ora_merge <- ora_list[[1]]
    ora_merge@organism <- species
    ora_merge@result <- ora_list %>% 
      map(\(x) slot(x, "result")) %>%
      purrr::reduce(.f = bind_rows)
    ora_merge@geneSets <- ora_list %>% 
      map(\(x) slot(x, "geneSets")) %>%
      purrr::reduce(.f = c)
    # split each database
    ora_merge@result <- ora_merge@result %>% 
      mutate(Category = map_chr(
        Description, \(x) str_split(x, "[_]", n=2, simplify=T) %>% .[1]),
        Description = map_chr(
          Description, \(x) str_split(x, "[_]", n=2, simplify=T) %>% .[2]) %>%
          map_chr(\(x) str_replace_all(x, "_", " ")) %>%
          map_chr(str_to_sentence) ) #%>%
    # mutate(Category = factor(Category, levels=unique(.$Category)))
    
    enrich_result[[names(gene_list)[i]]] <- ora_merge
  }
  
  return(enrich_result)
}



#' @param gene_list_ranked list of ranked named gene character 
#' @param geneset_list list of dataframe geneset
enrich_gsea_list <- function(gene_list_ranked, geneset_list) {
  enrich_result <- list()
  for(i in seq_along(gene_list_ranked)) {
    print(str_glue(">>>>>>>> GSEA enrichment on: {names(gene_list_ranked)[i]}"))
    
    genelist <- gene_list_ranked[[i]] %>% na.omit() %>% .[is.finite(.)]
    if(length(genelist) < 50) next
    
    gsea_list <- future_map(geneset_list, \(x) GSEA(
      genelist, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, 
      TERM2GENE = x$term2gene, TERM2NAME = x$term2name, 
      seed=T, verbose=F),
      .options = furrr_options(seed = T) )
    # merge geneset_list 
    gsea <- gsea_list[[1]]
    gsea@organism <- species
    gsea@result <- gsea_list %>% 
      map(\(x) slot(x, "result")) %>%
      purrr::reduce(.f = bind_rows)
    gsea@geneSets <- gsea_list %>% 
      map(\(x) slot(x, "geneSets")) %>%
      purrr::reduce(.f = c)
    # split each database
    gsea@result <- gsea@result %>%
      mutate(Category = map_chr(
        Description, \(x) str_split(x, "_", n=2, simplify=T) %>% .[1]),
        Description = map_chr(
          Description, \(x) str_split(x, "_", n=2, simplify=T) %>% .[2]) %>%
          map_chr(\(x) str_replace_all(x, "_", " ")) %>%
          map_chr(str_to_sentence) ) #%>%
    # mutate(Category = factor(Category, levels=unique(.$Category)))
    
    # to list
    enrich_result[[names(gene_list_ranked)[i]]] <- gsea
  }
  
  return(enrich_result)
}




#' @name plot_volcano
#' @title plot volcano
#' @description plot volcano plot for DESeq2 pariwise compare result
#'
#' @param res a data.drame from DESeq2 result
#' @param fc_thred foldchange cutoff, corelated to log2FoldChange column, high value mean less filter
#' @param padj_thred p.adjust cutoff, corelated to padj column, high value mean less filter
#' @param size point size
#' @param alpha point transparency
#' @param title plot title
#' @param xlim two element numeric vector, restrict x axis, default: NULL
#' @param ylim two element numeric vector, restrict y axis, default: NULL
#' @param ns_resampling numeric, downsampling NS points
#' @param color three element character vector, map point color, ordered by down-regulated, ns, up-regulated
#'
#' @importFrom dplyr filter mutate if_else sample_n bind_rows %>%
#' @importFrom rlang sym !!
#' @import ggplot2
#'
#' @return ggplot2 object
#'
#' @export
#'
plot_volcano <- function(res,
                         fc_thred=1.5,
                         padj_thred=0.05,
                         size=1.5,
                         alpha=0.7,
                         title="Volcano of DEGs",
                         xlim=NULL, # c(-5, 5)
                         ylim=NULL, # c(0, 20)
                         ns_resampling=2000,
                         color=c("#0571b0", "#bababa", "#ca0020")
) {
  df <- filter(res, !is.na(symbol), !is.na(padj)) %>%
    mutate(`-log10(padj)`=-log10(padj))
  # take sig genes
  df.sig <- filter(df, padj < padj_thred, abs(log2FoldChange) > log2(fc_thred)) %>%
    mutate(direction=if_else(log2FoldChange > 0, "up", "down"))
  # sample un-sig genes,  fo reducing rendering points
  df.unsig <- filter(df, padj >= padj_thred) %>%
    sample_n(size=ns_resampling) %>%
    mutate(direction="ns")
  # merge sig and un-sig
  df.merge <- bind_rows(df.sig, df.unsig)
  
  # set lims
  df.merge$shape <- "inrange"
  if(!is.null(xlim)) {
    df.merge <- df.merge %>%
      mutate(shape=if_else(log2FoldChange > xlim[2] | log2FoldChange < xlim[1], "outrange", "inrange")) %>%
      mutate(log2FoldChange=if_else(log2FoldChange > xlim[2], xlim[2],
                                    if_else(log2FoldChange < xlim[1], xlim[1], log2FoldChange)))
  }
  if(!is.null(ylim)) {
    df.merge <- df.merge %>%
      mutate(shape=if_else(`-log10(padj)` > ylim[2] | `-log10(padj)` < ylim[1], "outrange", shape)) %>%
      mutate(`-log10(padj)`=if_else(`-log10(padj)` > ylim[2], ylim[2],
                                    if_else(`-log10(padj)` < ylim[1], ylim[1], `-log10(padj)`)))
  }
  
  # plot
  plot <- ggplot(df.merge, aes(log2FoldChange, `-log10(padj)`)) +
    geom_point(aes(shape=shape, fill=direction), size=size, alpha=alpha, stroke=0) +
    scale_shape_manual(values=c(inrange=21, outrange=24)) +
    scale_fill_manual(values=c(down=color[1], ns=color[2], up=color[3])) +
    geom_vline(xintercept=c(-log2(fc_thred), log2(fc_thred)), linetype="dashed") +
    geom_hline(yintercept=-log10(padj_thred), linetype="dashed") +
    labs(x="log2 (FoldChange)", y="-log10 (p.adj)", 
         title=title, 
         subtitle=str_glue(
           "padj_threshold: {padj_thred}, FoldChange_threshold: {fc_thred}")) +
    theme_bw()
  
  if("ggrastr" %in% installed.packages()) plot <- ggrastr::rasterize(plot)
  
  return(plot)
}



#' @name plot_enrich_ora
#' @title plot_enrich_ora
#' @description dotplot for enrich result from enrich_ORA
#'
#' @param enrichment data.frame, from ORA enrichment,
#'        result can be merged from multi enrichment, grouping info can be viewed by facet
#' @param show_category integer, show enriched category number
#' @param plot_type character, one of dot, bar
#' @param axis_x character, plot axis x by  which column,
#'        one of GeneRatio, -log10(p.adjust), enrichFactor, Count
#' @param order_by character, filter show_category by which column,
#'        one of GeneRatio, -log10(p.adjust), enrichFactor, Count
#' @param color_by character, map point color, one of GeneRatio, -log10(p.adjust), enrichFactor, Count
#' @param color_range two element character vector, mapping color range, from low to high
#' @param size_by character, map point size, one of GeneRatio, -log10(p.adjust), enrichFactor, Count
#'        only works if plot_type is dot
#' @param size_range two element numeric vector, mapping point size, from low to high
#'        only works if plot_type is dot
#' @param text_len_limit integer, wrap y-axis text length
#' @param facet_by character, facet plot by which column in enrichment
#' @param facet_scales character, facet scales, one of free_y, free_x, free, fixed,
#'        only works if facet_by is not NULL
#'
#' @import ggplot2
#' @importFrom dplyr arrange mutate slice vars group_by
#' @importFrom stringr str_wrap
#' @importFrom rlang sym !!
#'
#' @return ggplot2 object
#'
#' @export
#'
plot_enrich_ora <- function(enrichment,
                            show_category=10,
                            plot_type="dot", # bar
                            axis_x="GeneRatio",
                            order_by="-log10(p.adjust)",
                            color_by="-log10(p.adjust)",
                            color_range=c("#377eb8", "#e41a1c"),
                            size_by="Count",
                            size_range=c(2, 6),
                            text_len_limit=40,
                            facet_by=NULL,
                            facet_scales="free"
) {
  
  if(!(class(enrichment) == "data.frame"))
    stop("enrichment must be a data.frame (result slot from clusterProfiler)")
  
  if(!(plot_type %in% c("dot", "bar")))
    stop("plot_type must be one of dor, bar")
  
  if(!(axis_x %in% c("GeneRatio", "-log10(p.adjust)", "enrichFactor", "Count")))
    stop("axis_x must be one of GeneRatio, -log10(p.adjust), enrichFactor, Count")
  
  if(!(order_by %in% c("GeneRatio", "-log10(p.adjust)", "enrichFactor", "Count")))
    stop("order_by must be one of GeneRatio, -log10(p.adjust), enrichFactor, Count")
  
  if(!(size_by %in% c("GeneRatio", "-log10(p.adjust)", "enrichFactor", "Count")))
    stop("size_by must be one of GeneRatio, -log10(p.adjust), enrichFactor, Count")
  
  if(!(color_by %in% c("GeneRatio", "-log10(p.adjust)", "enrichFactor", "Count")))
    stop("color_by must be one of GeneRatio, -log10(p.adjust), enrichFactor, Count")
  
  if(!is.null(facet_by)) {
    if(!(facet_scales %in% c("free", "free_x", "free_y", "fixed")))
      stop("facet_scales must be one of free, free_x, free_y, fixed")
  }
  
  # pharse fraction to decimal
  enrichment$GeneRatio <- sapply(enrichment$GeneRatio, function(x) eval(parse(text=x)))
  
  # group facet
  if(!is.null(facet_by)) {
    data <- group_by(enrichment, !!sym(facet_by)) %>%
      mutate(`-log10(p.adjust)`=-log10(p.adjust)) %>%
      arrange(desc(!!sym(order_by))) %>%
      dplyr::slice(seq_len(show_category)) %>%
      arrange(!!sym(axis_x)) %>%
      mutate(Description=factor(Description, levels=unique(.$Description)))
  } else {
    data <- mutate(enrichment, `-log10(p.adjust)`=-log10(p.adjust)) %>%
      arrange(desc(!!sym(order_by))) %>%
      dplyr::slice(seq_len(show_category)) %>%
      arrange(!!sym(axis_x)) %>%
      mutate(Description=factor(Description, levels=unique(.$Description)))
  }
  
  if(plot_type=="dot") {
    plot <- ggplot(data, aes(!!sym(axis_x), Description)) +
      geom_segment(aes(yend=Description), xend=0, linewidth=0.7, color="grey50", alpha=0.7) +
      geom_point(aes(size=!!sym(size_by), fill=!!sym(color_by)), shape=21) +
      scale_size_continuous(range=size_range) +
      scale_fill_gradient(low=color_range[1], high=color_range[2]) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=text_len_limit)) +
      theme_classic() +
      ylab("")
  } else {
    plot <- ggplot(data, aes(!!sym(axis_x), Description)) +
      geom_bar(aes(fill=!!sym(color_by)), stat="identity", color=NA, width=0.75) +
      scale_size_continuous(range=size_range) +
      scale_fill_gradient(low=color_range[1], high=color_range[2]) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=text_len_limit)) +
      theme_classic() +
      ylab("")
  }
  
  # facet plot
  if(!is.null(facet_by)) {
    plot <- plot + facet_wrap(vars(!!sym(facet_by)), 
                              scales=facet_scales, ncol=1)
  }
  
  # return
  return(plot)
}




#' @name plot_enrich_gsea
#' @title plot_enrich_gsea
#' @description dotplot for enrich result from enrich_GSEA
#'
#' @param enrichment data.frame, from GSEA enrichment,
#'        result can be merged from multi enrichment, grouping info can be viewed by facet
#' @param show_category integer, show enriched category number
#' @param plot_type character, one of dot, bar
#' @param filter_by character, filter show_category by which column,
#'        one of -log10(p.adjust), NES, Count
#' @param color_by character, map point color, one of -log10(p.adjust), NES, Count
#' @param color_range two element character vector, mapping color range, from low to high,
#'        for gsea plot, the value of NSE = 0 will be white color
#' @param size_by character, map point size, currently, only Count avaible for gseaResult
#' @param size_range two element numeric vector, mapping point size, from low to high
#'        only works if plot_type is dot
#' @param text_len_limit integer, wrap y-axis text length
#' @param facet_by character, facet plot by which column in enrichment,
#'        default avaible value is direction, user use add aditional column in enrichment
#' @param facet_scales character, facet scales, one of free_y, free_x, "free, fixed,
#'        only works if facet_by is not NULL
#'
#' @import ggplot2
#' @importFrom dplyr arrange mutate slice vars group_by
#' @importFrom stringr str_wrap str_split
#' @importFrom rlang sym !!
#'
#' @return ggplot2 object
#'
#' @export
#
plot_enrich_gsea <- function(enrichment,
                             show_category=10,
                             plot_type="dot",
                             filter_by="-log10(p.adjust)",
                             color_by="-log10(p.adjust)",
                             color_range=c("#377eb8", "#e41a1c"),
                             size_by="Core_Count", # not changable by now
                             size_range=c(2, 6),
                             text_len_limit=40,
                             facet_by=NULL, # direction
                             facet_scales="free"
) {
  
  if(!(class(enrichment) == "data.frame"))
    stop("enrichment must be a data.frame (result slot from clusterProfiler)")
  
  if(!(plot_type %in% c("dot", "bar")))
    stop("plot_type must be one of dor, bar")
  
  if(!(filter_by %in% c("-log10(p.adjust)", "NES", "Count")))
    stop("filter_by must be one of -log10(p.adjust), NES, Count")
  
  if(!(color_by %in% c("-log10(p.adjust)", "NES", "Count")))
    stop("color_by must be one of -log10(p.adjust), NES, Count")
  
  if(!is.null(facet_by)) {
    if(!(facet_scales %in% c("free", "free_x", "free_y", "fixed")))
      stop("facet_scales must be one of free, free_x, free_y, fixed")
  }
  
  # count core enrichment genes
  enrichment$Core_Count <- map_int(enrichment$core_enrichment, ~ length(str_split(.x, fixed("/"), simplify = T)))
  
  # group facet
  if(!is.null(facet_by)) {
    data <- mutate(enrichment,
                   `-log10(p.adjust)`=-log10(p.adjust),
                   direction=if_else(NES>=0, "Up-regulated", "Down-regulated")) %>%
      group_by(direction, !!sym(facet_by)) %>%
      arrange(desc(abs(!!sym(filter_by)))) %>%
      dplyr::slice(seq_len(show_category)) %>%
      arrange(NES) %>%
      mutate(Description=factor(Description, levels=unique(.$Description)))
  } else {
    data <- mutate(enrichment, `-log10(p.adjust)`=-log10(p.adjust),
                   direction=if_else(NES>=0, "Up-regulated", "Down-regulated")) %>%
      group_by(direction) %>%
      arrange(desc(abs(!!sym(filter_by)))) %>%
      dplyr::slice(seq_len(show_category)) %>%
      arrange(NES) %>%
      mutate(Description=factor(Description, levels=unique(.$Description)))
  }
  
  # glimpse(data)
  
  if(plot_type=="dot") {
    plot <- ggplot(data, aes(NES, Description)) +
      geom_segment(aes(yend=Description), xend=0, linewidth=0.7, color="grey50", alpha=0.7) +
      geom_point(aes(size=Core_Count, fill=!!sym(color_by)), shape=21) +
      scale_size_continuous(range=size_range) +
      scale_fill_gradient2(low=color_range[1], mid="white", high=color_range[2], midpoint=0) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=text_len_limit)) +
      theme_classic() +
      ylab("")
  } else {
    plot <- ggplot(data, aes(NES, Description)) +
      geom_bar(aes(fill=!!sym(color_by)), stat="identity", color=NA, width=0.75) +
      scale_fill_gradient2(low=color_range[1], mid="white", high=color_range[2], midpoint=0) +
      scale_y_discrete(labels=function(x) str_wrap(x, width=text_len_limit)) +
      theme_classic() +
      ylab("")
  }
  
  # facet plot
  if(!is.null(facet_by)) {
    plot <- plot + facet_wrap(vars(!!sym(facet_by)),
                              scales=facet_scales, ncol = 1)
  }
  
  # return
  return(plot)
}
