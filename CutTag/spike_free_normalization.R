
# devtools::install_github("stjude/ChIPseqSpikeInFree")
library(ChIPseqSpikeInFree)
library(tidyverse)
library(yaml)

rstudioapi::getActiveDocument()$path |> 
  dirname() |> setwd()

config <- "config.yaml" %>% read_yaml()
gsize_file <- config$index_rootdir %>% 
  file.path('hs/chm13/chm13v2.0.fa.gz.fai')

sample_info <- "../../samplesheet.csv" %>% 
  read_csv() %>% #glimpse()
  mutate(ID=str_c(sample, ".filtered.bam")) %>% 
  glimpse()

save_prefix="SpikeFree"
meta_file <- str_c(save_prefix, "_metaFile.txt")
tibble(ID=sample_info$ID, ANTIBODY="H3K27me3", GROUP=sample_info$group) %>% 
  write_tsv(meta_file)


# data will saved to bam directory
setwd("result/03_bam")

# run
ChIPseqSpikeInFree(bamFiles=sample_info$ID, 
                   chromFile=gsize_file,
                   metaFile=meta_file,
                   prefix=save_prefix,
                   cutoff_QC=1.2,
                   maxLastTurn=0.99,
                   ncores=6)


# # transform sf for bigwig/bedgraph normalization
# statistics <- "../statistics_CPM.csv" %>% 
#   read_csv() %>% #glimpse()
#   # transmute(sample, filtered_reads) %>% 
#   left_join(sample_info) %>%
#   print()
# 
# scale_factor <- str_c(save_prefix, "_SF.txt") %>% 
#   read_table() %>%
#   left_join(statistics) %>% 
#   mutate(scale_factor=1000000 / (filtered_reads * SF)) %>% 
#   glimpse()
# 
# # save
# scale_factor %>% 
#   write_csv(str_c(save_prefix, ".csv"))
# 
# scale_factor %>% 
#   dplyr::select(sample, scale_factor=scale_factor, SF) %>% 
#   write_csv(str_c(save_prefix, "_scale_factors.csv"))


