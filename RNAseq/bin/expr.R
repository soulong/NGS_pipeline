
# speed up if installing packages
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")


# load packages
package_need <- c("yaml", "here", "tximport", "tidyverse", "writexl", "readxl")
sapply(package_need, function(x) require(x, character.only=T))

# setup working directory
work_dir <- ifelse(interactive(), 
                   dirname(dirname(rstudioapi::getActiveDocumentContext()$path)), 
                   here::here())
setwd(work_dir)
cat(paste0("\n  working directory is set on: \n  ", getwd(), "\n\n"))

# read config file from bin/config.yml
config <- suppressWarnings(read_yaml("bin/config.yml"))

# create output directory
suppressWarnings(if(!dir.exists(config$output_dir))  dir.create(config$output_dir))


# load annotation data
if(config$species=="hs") {
  t2g <- read_xlsx("bin/annotation.ensembl.95.xlsx", sheet="GRCh38_tx2gene")
  annotation <- read_xlsx("bin/annotation.ensembl.95.xlsx", sheet="GRCh38_annotation") %>% .[, c(-3, -7, -8)]
} else {
  if(config$species=="mm") {
    t2g <- read_xlsx("bin/annotation.ensembl.95.xlsx", sheet="GRCm38_tx2gene")
    annotation <- read_xlsx("bin/annotation.ensembl.95.xlsx", sheet="GRCm38_annotation") %>% .[, c(-3, -7, -8)]
  } else {
    stop("unsupported species setting, please check species option in `bin/config.yml` !")
  }
}


# import salmon
# read parameters
# read sample information
sample <- read_xlsx("bin/user_parameters.xlsx", sheet="sample", col_types="text")

# check user_parameter if project_id and sample_name were unique
if(length(sample$project_id)!=unique(length(sample$project_id))) {
  cat("\n check user_parameter file, `project_id` were not unique!")
  exit()
} else {
  if(length(sample$sample_name)!=unique(length(sample$sample_name))) {
    cat("\n check user_parameter file, `sample_name` were not unique!")
    exit()
  }
}

if(dir.exists(config$salmon_dir)) {
  
  # read salmon files
  cat("\n read samlon files ...\n")
  salmon_dirs <- list.dirs(config$salmon_dir, recursive=F, full.names=F)
  salmon_quant_files <- file.path(config$salmon_dir, salmon_dirs, "quant.sf")
  names(salmon_quant_files) <- salmon_dirs
  
  # set project_id to sample_name in expression matrix
  for(i in 1:length(salmon_quant_files)) {
    names(salmon_quant_files)[i] <- sample$sample_name[match(names(salmon_quant_files)[i], sample$project_id)]
  }
  
  # re-order sample according to sample info in excel
  salmon_quant_files <- salmon_quant_files[match(sample$sample_name, names(salmon_quant_files))]
  
  # get transcript level
  tx_salmon.transcript <- tximport(salmon_quant_files, type="salmon", txOut=T, countsFromAbundance="no")
  
  cat("\n extract trancript level expression\n")
  transcript.counts.anno <- tx_salmon.transcript$counts %>% as_tibble(rownames="ensembl_transcript_id") %>%
    left_join(., t2g, by=c("ensembl_transcript_id"="TXNAME")) %>% 
    left_join(., annotation[, c("ensembl_gene_id", "external_gene_name")], by=c("GENEID"="ensembl_gene_id")) %>%
    dplyr::select(ensembl_transcript_id, ensembl_gene_id=GENEID, external_gene_name, everything()) %>%
    mutate_if(is.numeric, round, 0)
  transcript.tpm.anno <- tx_salmon.transcript$abundance %>% as_tibble(rownames="ensembl_transcript_id") %>%
    left_join(., t2g, by=c("ensembl_transcript_id"="TXNAME")) %>% 
    left_join(., annotation[, c("ensembl_gene_id", "external_gene_name")], by=c("GENEID"="ensembl_gene_id")) %>%
    dplyr::select(ensembl_transcript_id, ensembl_gene_id=GENEID, external_gene_name, everything()) %>%
    mutate_if(is.numeric, round, 2)
  
  # get gene level
  cat("\n summarize to gene level expression\n")
  tx_salmon.gene <- summarizeToGene(tx_salmon.transcript, t2g, ignoreTxVersion=T, countsFromAbundance="no")
  
  # annotate gene
  cat("\n annotate genes")
  tx_salmon.gene.counts.anno <- tx_salmon.gene$counts  %>% 
    as_tibble(rownames="ensembl_gene_id") %>%
    left_join(., annotation, by=c("ensembl_gene_id"="ensembl_gene_id")) %>%
    dplyr::select(ensembl_gene_id, external_gene_name, entrezgene_id, gene_biotype, description, everything()) %>%
    mutate_if(is.numeric, round, 0)
  tx_salmon.gene.tpm.anno <- tx_salmon.gene$abundance %>% 
    as_tibble(rownames="ensembl_gene_id") %>%
    left_join(., annotation, by=c("ensembl_gene_id"="ensembl_gene_id")) %>%
    dplyr::select(ensembl_gene_id, external_gene_name, entrezgene_id, gene_biotype, description, everything()) %>%
    mutate_if(is.numeric, round, 2)
  
  # save to file
  cat("\n write to disk")
  write_xlsx(list(tpm=transcript.tpm.anno, counts=transcript.counts.anno), 
             path=file.path(config$output_dir, str_c(Sys.Date(),"_transcript_expression.xlsx")))
  write_xlsx(list(tpm=tx_salmon.gene.tpm.anno, counts=tx_salmon.gene.counts.anno),
             path=file.path(config$output_dir, str_c(Sys.Date(),"_gene_expression.xlsx")))
  
  cat("\n salmon quantification were done!\n\n")
  
} else { 
  cat("\n no salmon output directory found, please check related files!")
  exit()
 }



##### import featurecounts #####
if(file.exists(str_c(config$featurecounts_dir,"/gene_featurecounts.txt"))) {
  print("pharse featurecounts ...")
  
  # read featurecounts result
  gene_count_featurecounts <- suppressMessages(read_table2(paste0(config$featurecounts_dir,"/gene_featurecounts.txt"), comment = "#"))
  
  colnames(gene_count_featurecounts)[-1:-6] <- sapply(colnames(gene_count_featurecounts)[-1:-6], 
                                                          function(x) str_split(x, "/", simplify = T)[3])
  gene_count_featurecounts$Geneid <- sapply(gene_count_featurecounts$Geneid, function(x) str_split(x,"[.]",simplify=T)[1])
  
  # repolace sanmple name
  name_index <- match(colnames(gene_count_featurecounts)[-1:-6], sample$project_id)
  if(all(is.na(name_index))) stop("incompatible project_id and sample_name, please check it!")
  colnames(gene_count_featurecounts)[-1:-6] <- sample$sample_name[name_index]


  # calculate TPM values
  tpm <- function(counts, lengths) {
    rate <- counts / lengths
    rate / sum(rate) * 1e6
  }
  tpm_featurecounts <- gene_count_featurecounts %>%
    gather(sample, counts, 7:ncol(gene_count_featurecounts)) %>%
    group_by(sample) %>%
    mutate(tpm=tpm(counts, Length)) %>%
    dplyr::select(-counts) %>%
    spread(sample, tpm) %>% 
    dplyr::select(-2:-6) %>%
    left_join(annotation, by=c("Geneid"="ensembl_gene_id")) %>%
    dplyr::select(ensembl_gene_id=Geneid, external_gene_name, entrezgene_id,  
                  gene_biotype, description, everything()) %>%
    mutate_if(is.numeric, round, 2)
  tpm_featurecounts <- dplyr::select(tpm_featurecounts, 1:5, 
                                     match(sample$sample_name, colnames(tpm_featurecounts)))
  
  
  counts_featurecounts <- dplyr::select(gene_count_featurecounts, -2:-6) %>%
    left_join(annotation, by=c("Geneid"="ensembl_gene_id")) %>%
    dplyr::select(ensembl_gene_id=Geneid, external_gene_name, entrezgene_id, 
                  gene_biotype, description, everything()) %>%
    mutate_if(is.numeric, round, 0)

  counts_featurecounts <- dplyr::select(counts_featurecounts, 1:5, 
                                     match(sample$sample_name, colnames(counts_featurecounts)))
  
  
  cat("\n save file\n")
  write_xlsx(list(tpm=tpm_featurecounts, counts=counts_featurecounts),
             path=file.path(config$output_dir, str_c(Sys.Date(),"_gene_expression_bam.xlsx")))
  
  cat("\n featurecounts quantification were done!\n\n")

} else {
  print("no featurecounts output file found, skip it!")
 }


cat("\n all jobs were done!\n\n")


