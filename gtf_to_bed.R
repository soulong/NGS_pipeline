
# processing gtf to bed
library(tidyverse)
library(GenomicRanges)
library(plyranges)

setwd("Y:/index/hs/chm13")
# setwd("D:/Index/hs/v49")
# setwd("D:/Index/mm/vM38")

gtf <- list.files(pattern="gtf.gz") %>% .[1] %>% 
  rtracklayer::import() %>% 
  print()

gene_name <- 'gene_name'


gtf$source %>% table()
gtf$level %>% table()
gtf$type %>% table()
gtf$gene_type %>% table()
gtf$tag %>% table()

genes <- gtf %>% 
  filter(type=='gene' & (! seqnames %in% c("chrM"))) %>% 
  filter(! is.na(!!as.name(gene_name))) %>% 
  print()

genes@elementMetadata[[gene_name]] %>% unique() %>% length()

# tss <- promoters(genes, 0, 0)
table(genes$gene_biotype)
table(genes$source)


# save to UCSC style bed
genes %>%
  filter(gene_biotype=='protein_coding') %>% 
  .[, gene_name] %>% 
  rtracklayer::export("genes_protein_coding.bed")

genes %>%
  filter(gene_biotype=='lncRNA') %>% 
  .[, gene_name] %>% 
  rtracklayer::export("genes_lncRNA.bed")

