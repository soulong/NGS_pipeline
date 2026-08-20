

rstudioapi::getActiveDocumentContext()$path |>
  dirname() |> setwd()
getwd()

library(tidyverse)
library(magrittr)
library(ChIPseeker)
# library(ChIPQC)
# library(DiffBind)
library(plyranges)
# library(EnrichedHeatmap)
# library(profileplyr)
library(eulerr)
library(ggrastr)
# library(cowplot)

## multiprocessing -------------------------
ncore <- max(1, trunc(parallel::detectCores() / 2))

# library(furrr)
# plan(multisession, workers=ncore)
# plan()

library(BiocParallel)
if(Sys.getenv()['OS'] == 'Windows_NT') {
  bpparam <- SnowParam(ncore) 
} else {
  bpparam <-  MulticoreParam(ncore)
}
print(bpparam)


## read config from yml ----------------------
config <- yaml::read_yaml("bin/config.yml")
bam_dir <- config$bam_dir # "bam"
peak_dir <-  config$peak_dir # "peak"
species <- config$species
cat(bam_dir, peak_dir, species)

# load txdb
if(species == 'hs') {
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  annodb <- 'org.Hs.eg.db'
} else {
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
  annodb <- 'org.Mm.eg.db'
}


library(fs)
peak_files <- dir_ls('merged_peak', type='file', glob='*.narrowPeak')
names(peak_files) <- peak_files %>%
  path_file() %>%
  str_replace_all('_peaks.narrowPeak','')
names(peak_files)


## ********************************-------------------------------------
## ChipSeeker            ----------------------------------------------
## ********************************-------------------------------------
list(YAP=c("white","#e41a1c"),
     TAZ=c("white","#377eb8"),
     panTEAD=c("white","#4daf4a"),
     H3K27AC=c("white","#984ea3"))

peaks <- map(peak_files, readPeakFile)

peaks <- map(peaks, function(x) {seqlevelsStyle(x) <- 'UCSC'; return(x)})

tss <- getBioRegion(txdb, 
                    upstream=3000, downstream=3000,
                    by='gene', type='start_site')

tag_matrix_list <- map(peaks, 
                       ~ getTagMatrix(.x, windows=tss, weightCol='V5'))

# plot
pdf(str_glue('{Sys.Date()}_merged_peaks_tss_prof.pdf'), width=6, height=2)
plotPeakProf(tag_matrix_list[1:3], nbin=200, facet='column', free_y=F)
plotPeakProf(tag_matrix_list[4:6], nbin=200, facet='column', free_y=F)
plotPeakProf(tag_matrix_list[7:9], nbin=200, facet='column', free_y=F)
plotPeakProf(tag_matrix_list[10:12], nbin=200, facet='column', free_y=F)
while(dev.cur() != 1) dev.off()

pdf(str_glue('{Sys.Date()}_merged_peaks_tss_heatmap.pdf'), width=24, height=6)
tagHeatmap(tag_matrix_list, xlim=c(-3000, 3000), 
           color=c(rep('#984ea3',3), rep('#4daf4a',3), rep('#377eb8',3), rep('#e41a1c',3)))
while(dev.cur() != 1) dev.off()


