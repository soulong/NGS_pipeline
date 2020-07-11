
# speed up if installing packages
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

#################### parameters here
# sig genes thredhold
deg_padj <- 0.05
deg_fc <- 1.5
mincount <- 100
lrt_padj <- 0.05
# write to file
savefile <- T
#################### parameters here

# if need to install pacman
if (!require("pacman")) install.packages("pacman")
# load packages
package_need <- c("yaml", "here", "rlist", "rjson", "tximport", 
                  "readxl", "writexl", "DESeq2", "DEGreport", 
                  "pheatmap", "tidyverse", "corrplot", "RColorBrewer", 
                  "ggplot2", "ggrepel", "cowplot", "ggfortify", "cluster")
pacman::p_load(char=package_need)

# setup working directory
work_dir <- ifelse(interactive(), 
                   dirname(dirname(rstudioapi::getActiveDocumentContext()$path)), 
                   here::here())
setwd(work_dir)
cat(paste0("\n  working directory is set on: ", getwd()))

# read config file from bin/config.yml
config <- suppressWarnings(read_yaml("bin/config.yml"))

# load annotation data
if(config$species=="hs") {
  annotation <- read_xlsx("bin/annotation.ensembl.95.xlsx", sheet="GRCh38_annotation") %>% .[, c(-3, -7, -8)]
  orgdb <- "oeg.Hs.eg.db"
} else {
  if(config$species=="mm") {
    annotation <- read_xlsx("bin/annotation.ensembl.95.xlsx", sheet="GRCm38_annotation") %>% .[, c(-3, -7, -8)]
    orgdb <- "oeg.Mm.eg.db"
  } else {
    stop("unsupported species setting, please check species option in `bin/config.yml` !")
  }
}

# read contrast info
contrast_info <- read_xlsx("bin/user_parameters.xlsx", sheet="contrast", range="A1:C100", col_types="text")
comparison <- contrast_info$comparison[!is.na(contrast_info$comparison)]
contrasts <- list()
for (i in 1:length(na.omit(contrast_info$contrast))) {
  contrasts <- list.append(contrasts, c(contrast_info$contrast[i], contrast_info$control[i]))
}


# read sample info
sample <- read_xlsx("bin/user_parameters.xlsx", sheet="sample", col_types="text")
sample <- sample[ , c("project_id", "sample_name", "include", "batch", comparison)] %>%  filter(include=="yes")

# construct sample contrasts accroding to sample info colData
sample$contrasts <- apply(sample[, comparison], 1, paste, collapse=".")
# define factor levels acorrding to sample info order
sample$contrasts <-  factor(sample$contrasts, levels=unique(sample$contrasts))

# check if contrast has its relative items in sample_info$condition
if(!all(as.factor(unlist(contrasts)) %in% sample$contrasts)) {
  stop("not all 'contrast'_vs_'control' were found in sample definition, please check it!")
}

# prepare colData
colData <- data.frame(condition=sample$contrasts, batch=sample$batch)
rownames(colData) <- sample$sample_name
print(colData)


# get gene counts file path
expr_file <- try(grep("^[0-9].*_gene_expression.xlsx", list.files(config$output_dir, ".xlsx$"), value=T))
# check file exist
if(length(expr_file)==0) {
  stop("can't find `_gene_expression.xlsx` file, make sure it's exist in output directory!")
} else {
  # check unique
  if(length(expr_file)>1) {
    stop("there's more than one `_gene_expression.xlsx` file, please keep only one wanted file!")
  }
}

# read counts matrix
counts <- read_xlsx(file.path(config$output_dir, expr_file), sheet="counts")

# calculatye mitocondria counts ratio
mt.genes.prefix <- if_else(config$species=="hs", "^MT-", "^mt-")
mt.genes <- str_replace_na(counts$external_gene_name) %>% 
  str_detect(mt.genes.prefix) %>% 
  str_replace_all("TRUE", "mt.counts") %>%
  str_replace_all("FALSE", "nuclei.counts") %>%
  enframe(name=NULL, value="mt.genes")
counts.mt <- counts[, c(-2:-5)] %>%
  bind_cols(mt.genes) %>%
  gather(key="sample", value="counts", -ensembl_gene_id, -mt.genes) %>% 
  group_by(mt.genes, sample) %>% summarise(sum=sum(counts))
counts.mt$sample <- factor(counts.mt$sample, levels=rev(colnames(counts)[c(-1:-5)]))

# remove samples not include (inculde=no) and convert to matrix
counts <- counts[, c(colnames(counts)[1], sample$sample_name)] %>% column_to_rownames(var="ensembl_gene_id")


## import to deseq2
deseq.matrix <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design= ~ batch + condition)

## pair-wise wald test
cat("  perform wald test\n")
deseq.wald <- DESeq(deseq.matrix, parallel=T)

# vst transformation
vst <- vst(deseq.wald, blind=F)
vst_mat <- assay(vst)
pca <- plotPCA(vst, intgroup="condition", ntop=2000, returnData=F) + 
  theme_cowplot() +
  stat_ellipse(aes(fill=group), geom="polygon", alpha=1/4, level=0.95) +
  theme(legend.position="top") +
  geom_text_repel(aes(label=name))

vst_mat.batch <- vst
assay(vst_mat.batch) <- limma::removeBatchEffect(assay(vst_mat.batch), vst_mat.batch$batch)
pca.batch <- plotPCA(vst_mat.batch, intgroup="condition", ntop=2000, returnData=F) +
  theme_cowplot() +
  stat_ellipse(aes(fill=group), geom="polygon", alpha=1/4, level=0.95) +
  theme(legend.position="top") +
  geom_text_repel(aes(label=name))


# deseq2 wald test result extract
deg_result_extract <- function(dds, treatment, control, anno) {
  con <- c("condition", treatment, control)
  deg_unshrink <- results(dds, contrast=con)
  deg_shrink <- lfcShrink(dds, contrast=con, res=deg_unshrink, parallel=T, quiet=T) %>% 
    data.frame() %>% rownames_to_column(var="ensembl_gene_id") %>% as_tibble()
  deg_shrink <- right_join(anno, deg_shrink, by=c("ensembl_gene_id"="ensembl_gene_id")) %>% arrange(padj)
  deg_shrink_sig <- filter(deg_shrink, baseMean > mincount & abs(log2FoldChange) > log2(deg_fc) & padj < deg_padj)
  return(list(deg=deg_shrink, deg_sig=deg_shrink_sig)) 
  }

# statistics comparation
cat("  perform deg analysis\n")
deg_results <- list()
for (i in 1:length(contrasts)) {
  re <- deg_result_extract(deseq.wald, treatment=contrasts[[i]][1], control=contrasts[[i]][2], anno=annotation)
  deg_result <- list(deg=re$deg, deg_sig=re$deg_sig, treatment=contrasts[[i]][1], control=contrasts[[i]][2])
  deg_results <- list.append(deg_results, deg_result) 
  }


## lrt test
cat("  perform likehood ratio test\n")
deseq.lrt <- DESeq(deseq.matrix, test="LRT", reduced= ~ 1, parallel=T)

# lrt result
lrt <- results(deseq.lrt) %>% data.frame() %>% 
  rownames_to_column(var="ensembl_gene_id") %>% as_tibble()
lrt <- right_join(annotation, lrt, by=c("ensembl_gene_id"="ensembl_gene_id"))
lrt_sig <- filter(lrt, !is.na(padj) & baseMean > mincount & padj < lrt_padj) %>% arrange(padj)

# lrt rlog transformation for degPatterns
rld <- rlog(deseq.lrt, blind=F)
rld_mat <- assay(rld)

# deg patterns
cat("  perform degPattern cluster\n")
lrt_pattern_clusters <- degPatterns(rld_mat[lrt_sig$ensembl_gene_id, ], 
                                    metadata=colData, time="condition", reduce=T, plot=F) 

# save cslusters
clusters.anno <- right_join(annotation, lrt_pattern_clusters$df, by=c("ensembl_gene_id"="genes"))
tpm <- read_xlsx(file.path(config$output_dir, expr_file), sheet="tpm")
clusters.anno.tpm <- left_join(clusters.anno, tpm[, c(-2, -3, -4)], by=c("ensembl_gene_id"="ensembl_gene_id"))
  

## save to file
if(isTRUE(savefile)) {
  
  cat("  save data\n\n")

  # save R result
  save(deseq.wald, deg_results, deseq.lrt, lrt_sig, lrt_pattern_clusters,
       file=file.path(config$output_dir, str_c(Sys.Date(),"_deseq2.RData")))
  
  # create dir
  statistics_dir <- file.path(config$output_dir, "statistics")
  if(!dir.exists(statistics_dir)) dir.create(statistics_dir)
  
  # write deg result
  for (i in 1:length(deg_results)) {
      write_xlsx(list(all=deg_results[[i]]$deg, sig=deg_results[[i]]$deg_sig), 
                 path=file.path(statistics_dir, str_c(Sys.Date(), "_wald_TEST_", 
                      deg_results[[i]]$treatment, "_", deg_results[[i]]$control, ".xlsx")))
  }
  
  # write lrt result
  write_xlsx(list(all=lrt, sig=lrt_sig), 
             path=file.path(statistics_dir, str_c(Sys.Date(),"_likelihood_ratio_TEST.xlsx")))
  
  # write cluster result
  write_xlsx(clusters.anno.tpm, 
             path=file.path(statistics_dir, str_c(Sys.Date(),"_degPattern_clusters.xlsx")))
  
  # close device
  while (!is.null(dev.list())) {
    dev.off()
  }

  # open pdf
  pdf(file.path(config$output_dir, str_c(Sys.Date(),"_deseq2_summary.pdf")), useDingbats=F)
  
  p1 <- ggplot(counts.mt, aes(x=sample, y=sum, fill=mt.genes)) + 
    geom_bar(stat="identity", position="stack") + 
    coord_flip() + 
    theme(legend.position="top") +
    labs(x="", y="Counts")
  p2 <- ggplot(counts.mt, aes(x=sample, y=sum, fill=mt.genes)) + 
    geom_bar(stat="identity", position="fill") + 
    geom_hline(yintercept=0.9, size=0.4, color="grey50") +
    geom_hline(yintercept=0.95, size=0.4, color="grey50") +
    coord_flip() + 
    theme(legend.position="top") + 
    labs(x="", y="Percentage")
  mt_plot <- plot_grid(p1, p2, nrow=1)
  
  print(mt_plot)
  
  color_palette <- colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
  
  # PCA plot
  print(pca)
  
  # PVA plot for correct batch
  print(pca.batch)

  # correlation plot
  sample_cor <- cor(vst_mat)
  corrplot(sample_cor, cl.lim=c(min(sample_cor), max(sample_cor)), is.corr=F)
  
  # top variable genes heatmap
  topVarGenes <- head(order(rowVars(vst_mat), decreasing=T), 1000)
  topVarGenes_mat <- vst_mat[topVarGenes, ]
  # Annotate our heatmap
  annotation_col <- data.frame(row.names=row.names(colData), sampe_group=as.character(colData$condition))
  pheatmap(topVarGenes_mat, color=color_palette, show_rownames=F,
           annotation_col=annotation_col, border_color=NA, scale="row",
           fontsize=10,  main="top 1000 variable genes")
  
  # deg pattern clusters
  print(lrt_pattern_clusters$plot + theme(legend.position="none"))
  
  # disperation plot
  plotDispEsts(deseq.wald)
  
  dev.off()
}

