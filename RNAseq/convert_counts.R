
rstudioapi::getActiveDocumentContext()$path |>
  dirname() |>
  setwd()

library(tidyverse)
library(org.Hs.eg.db)
library(GenomicFeatures)

# get counts
raw <- read_csv("raw_counts_matrix.csv")
# strip dot version suffix
raw$Geneid <- str_replace(raw$Geneid,
                          pattern=".[0-9]+$",
                          replacement="")
# # merge
# raw <- left_join(raw1, raw2, by=join_by(ensembl==Geneid)) %>% 
#   filter(!if_any(everything(), is.na))
glimpse(raw)



# ensembl to gene symbol
gene_symbol <- AnnotationDbi::select(
  org.Hs.eg.db, keys=raw$Geneid, keytype="ENSEMBL", 
  columns=c("SYMBOL"),
  multiVals = "first") %>% 
  as_tibble() %>% 
  magrittr::set_colnames(c("ensembl", "symbol"))
gene_symbol


# get gene length (skip intron)
txdb <- makeTxDbFromGFF("D:/Index/hs/GRCh38.gtf.gz")
exons_list_per_gene <- exonsBy(txdb, by="gene")
gene_length <- exons_list_per_gene %>% 
  reduce() %>% 
  width() %>% 
  sum() %>% 
  as.data.frame(col.names="length") %>% 
  as_tibble(rownames="ensembl") %>% 
  magrittr::set_colnames(c("ensembl", "length")) %>% 
  mutate(ensembl=str_replace_all(ensembl, ".[0-9]+$", ""))
gene_length




# merge
counts <- right_join(gene_symbol, raw, join_by(ensembl==Geneid)) %>%
  right_join(gene_length, ., join_by(ensembl==ensembl)) %>% 
  mutate(symbol=ifelse(is.na(symbol), ensembl, symbol)) %>% 
  dplyr::select(!ensembl) %>% 
  filter(!is.na(length)) %>% 
  distinct(symbol, .keep_all = T)
counts


# rename sample
sample_info <- read_excel("sample_info.xlsx", "sample")
names(counts)[-1:-2] <- sample_info$sample_name[
  match(names(counts)[-1:-2], sample_info$project_id)]
counts



# calculate TPM
get_tpm <- function(counts, lengths) {
  rate <- counts / lengths
  return(round(rate / sum(rate) * 1e6, 2))
  }

tpm <- counts %>%
  mutate(across(3:last_col(), \(x) get_tpm(x, length)))
tpm



# save
writexl::write_xlsx(
  list(counts=dplyr::select(counts, !length),
       tpm=dplyr::select(tpm, !length)),
  str_glue("quantification/{Sys.Date()}_gene_expression.xlsx"))




  
