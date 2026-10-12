
# processing gtf to bed — builds the BED region files that helper.sh
# setup_genome expects in each index directory:
#   genes.bed, genes_protein_coding.bed, genes_lncRNA.bed
# BED column 4 (name) carries the gene_name.
library(tidyverse)
library(GenomicRanges)
library(plyranges)

# genome directory to process (edit, or override via the GTF_TO_BED_DIR
# environment variable)
work_dir <- Sys.getenv("GTF_TO_BED_DIR", unset = "Y:/index/hs/chm13")
setwd(work_dir)

# pick THE gtf.gz of the directory — warn instead of silently taking the
# first match when several are present
gtf_files <- list.files(pattern = "gtf.gz")
if (length(gtf_files) == 0) stop("no gtf.gz found in ", getwd())
if (length(gtf_files) > 1) {
  warning("multiple gtf.gz in ", getwd(), "; using the first: ", gtf_files[1],
          " (all: ", paste(gtf_files, collapse = ", "), ")")
}
gtf <- gtf_files[1] %>%
  rtracklayer::import() %>%
  print()

gene_name <- 'gene_name'


gtf$source %>% table()
gtf$level %>% table()
gtf$type %>% table()
gtf$gene_type %>% table()
gtf$tag %>% table()

# drop the mitochondrial contig under any common naming (UCSC chrM, Ensembl
# MtDNA / mitochondrion_genome) — the old filter only knew "chrM", so fly/cel
# assemblies kept their mitochondrial gene silently
mito_pattern <- "chrM$|MtDNA|mitochondrion"

genes <- gtf %>%
  filter(type=='gene' & (! str_detect(as.character(seqnames), mito_pattern))) %>%
  filter(! is.na(!!as.name(gene_name))) %>%
  print()

genes@elementMetadata[[gene_name]] %>% unique() %>% length()

# tss <- promoters(genes, 0, 0)
table(genes$gene_biotype)
table(genes$source)


# save to UCSC style bed.
# rtracklayer only maps an mcol literally named `name` to BED column 4 —
# exporting with the gene_name mcol left the name column empty.
export_genes_bed <- function(genes, file) {
  out <- genes[, gene_name]
  names(out) <- mcols(out)[[gene_name]]
  rtracklayer::export(out, file)
  message("wrote ", file, ": ", length(out), " regions")
}

# helper.sh expects genes.bed (all genes) as the first heatmap BED_REGIONS
# entry; it was never produced before, so that QC heatmap silently never ran
export_genes_bed(genes, "genes.bed")

genes %>%
  filter(gene_biotype=='protein_coding') %>%
  export_genes_bed("genes_protein_coding.bed")

genes %>%
  filter(gene_biotype=='lncRNA') %>%
  export_genes_bed("genes_lncRNA.bed")
