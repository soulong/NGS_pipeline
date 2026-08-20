
rstudioapi::getActiveDocumentContext()$path |>
  dirname() |> setwd()
getwd()

library(tidyverse)
library(EnrichedHeatmap)
library(profileplyr)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(RColorBrewer)


config <- yaml::read_yaml("bin/config.yml")
bam_dir <- config$bam_dir # "bam"
peak_dir <-  config$peak_dir # "peak"
species <- config$species

peak_dir <- "callpeak" # this is call masc2 using clean.bam
print(str_glue("{bam_dir}, {peak_dir}, {species}"))



###################### get peak ###########################
peak_files <- list.files(peak_dir, ".narrowPeak$", full.names = T) %>%
  set_names(., map_chr(., ~ basename(.x) %>% str_replace_all("_peaks.narrowPeak", "")))

peak_files <- peak_files[c(4,9,5,11)]

print(peak_files)

# extract info
info <- enframe(peak_files, "sample", "path") %>%
  separate(sample, c("treatment", "target","factor"), sep="_", remove = F) %>%
  # group is used to identify replicates
  unite("group", target, treatment, sep="_", remove = F) %>%
  arrange(target, treatment, factor)
# reorder files according to info
peak_files <- peak_files[info$sample]

# get peak grange list
peak_list <- map(peak_files, readPeakFile) # using chipseeker
# peaks <- GRangesList(peak_list)

# keep autosome
peak_list <- map(peak_list, function(x) { 
  keepStandardChromosomes(x, species=ifelse(species=="hs", "Homo_sapiens", "Mus_musculus"), pruning.mode="coarse") %>%
    dropSeqlevels("MT", pruning.mode="coarse") })
seqlevels(peak_list[[1]])

# get merged peaks
peak_yap_taz <- GRangesList(peak_list) %>% 
  unlist() %>% 
  GenomicRanges::reduce()

# sample peaks overlap
for(i in names(peak_list)) {
  mcols(peak_yap_taz)[[i]] <- overlapsAny(peak_yap_taz, peak_list[[i]]) }

# subset peaks based on sample consensus
if(T) {
  select <- 1:4 #c(5, 6)
  peak_yap_taz <- peak_yap_taz[rowSums(as.data.frame(mcols(peak_yap_taz)[, select])) > 0, select]

}

# # change seq style
# seqlevelsStyle(peak_yap_taz)
# peak_yap_taz <- map(peak_yap_taz, function(x) { 
#   seqlevelsStyle(x) <- "UCSC"; return(x) })
# seqlevelsStyle(peak_yap_taz)
# seqlevels(peak_yap_taz)

# save to bed for chiprofile use
rtracklayer::export(peak_yap_taz, "peak_yap_taz.bed")





###################### profile bam ###########################
bam_files <- list.files(bam_dir, ".bam$", full.names = T) %>%
  # str_subset(".clean", negate = F) %>%
  set_names(., map_chr(., ~ basename(.x) %>% str_replace_all(".bam", "")))

bam_files <- bam_files[c(4,9,5,11)]

print(bam_files)

# coverage to chiprofile
# here, BamBigwig_to_chipProfile can't take named signalFiles as input, 
# otherwise, error will appear
names(bam_files) <- NULL
chipProfile <- BamBigwig_to_chipProfile(bam_files, # bw_files 
                                        "peak_yap_taz.bed",
                                        format = "bam",
                                        style = "point", 
                                        nOfWindows = 100, bin_size = 20, distanceAround = 3000,
                                        removeDup = TRUE, paired = TRUE, normalize = "RPM",
                                        quant_params = SnowParam(parallel::detectCores()-4) # for win
                                        # quant_params = MulticoreParam(round(parallel::detectCores()/2))
                                        )
# save in case of corruption
file_name <- "chipProfile_bam.rds"
if(!file.exists(file_name)) {
  write_rds(chipProfile, file_name)
  # chipProfile <- read_rds(file_name) \
  }

# convert to profileplyr
profile <- as_profileplyr(chipProfile)

# convert style to UCSC
seqlevelsStyle(rowRanges(profile)) <- "UCSC" # rowRanges(profile)
rowRanges(profile)

# rename sample
rownames(sampleData(profile)) <- rownames(sampleData(profile)) %>%
  map_chr(~ str_split(.x, fixed("_Fix.bam"), simplify = T)[1])

# annotate regions
profile <- annotateRanges(profile, 
                          TxDb = "hg38",
                          # annoDb = org.Hs.eg.db,
                          tssRegion = c(-3000, 3000))

# merge exon and intron
profile@rowRanges$Type <- 
  profile@rowRanges$annotation_short %>%
  str_replace_all("(Exon)|(Intron)|(5p UTR)|(3p UTR)|(Downstream)", "Gene Body") %>%
  factor(levels = c("Promoter", "Gene Body", "Distal Intergenic"))
table(profile@rowRanges$Type)

# summarize profile for clustering
profile_rowmeans <- assays(profile) %>%
  as.list() %>%
  map_dfc(~ rowMeans(.x))

# remove low value rows
low_rows <- (rowSums(profile_rowmeans) < 1)
if(sum(low_rows) > 0) {
  print(str_glue("{sum(low_rows)} rows have low values [rowMeans < 1]"))
  profile_rowmeans <- profile_rowmeans[!low_rows, ]
  profile <- profile[!low_rows, ]
}

# # group sample by target
# sampleData(profile)$target <- factor(c("TAZ","TAZ","YAP","YAP"))

# group row by annotation
profile@params$rowGroupsInUse <- "Type"


# all color equal
pdf(file="profile_by_annotation.pdf", width = 7, heigh = 7, useDingbats = T)
generateEnrichedHeatmap(
  profile,
  include_group_annotation = T,
  matrices_axis_name = c("-3000", "Center", "3000"),
  group_anno_color = brewer.pal(3, "Dark2"),
  all_color_scales_equal = T,
  matrices_color = as.list(rep(list(c("white","red2")), nrow(profile@sampleData))),
  use_raster = T,
  raster_quality = 2,
  return_ht_list = F)
dev.off()

# color equal by each target
pdf(file="profile_by_annotation_and_target.pdf", width = 7, heigh = 7, useDingbats = T)
generateEnrichedHeatmap(
  profile,
  include_group_annotation = T,
  matrices_axis_name = c("-3000", "Center", "3000"),
  group_anno_color = brewer.pal(3, "Dark2"),
  # color_by_sample_group = "target",
  all_color_scales_equal = T,
  matrices_color = list(c("white","#e41a1c"),c("white","#e41a1c"), c("white","#377eb8"),  c("white","#377eb8")),
  genes_to_label = c("CCN1","CCN2","ANKRD1","CTCF"),
  use_raster = T,
  raster_quality = 2,
  return_ht_list = F)
dev.off()



# subset promoters & YAP
# profile_yap_promoter <- profile[profile@rowRanges$Type == "Promoter",][,,3:4]
profile_yap <- profile[,,3:4]
profiler_yap <- profile_yap[profile_yap %over% peak_merge[rowSums(as.data.frame(mcols(peak_merge[,3:4]))) > 0, ],]
tmp <- assays(profiler_yap) %>%
  as.list() %>%
  map_dfc(~ rowMeans(.x))
profiler_yap <- profiler_yap[rowSums(tmp) > 1, ]
cluster <- kmeans(tmp[rowSums(tmp) > 1, ], centers = 1)
mcols(profiler_yap)$cluster <- as.factor(cluster$cluster)
profiler_yap@params$rowGroupsInUse <- "sgGroup" # "annotation_short"
pdf(file="profile_yap_by_cluster.pdf", width = 9, heigh = 7, useDingbats = T)
generateEnrichedHeatmap(
  profiler_yap,
  include_group_annotation = F,
  # group_anno_color = brewer.pal(3, "Dark2")[1:2],
  matrices_color = as.list(rep(list(c("white","#e41a1c")), nrow(profiler_yap@sampleData))), # list(c("white","blur2"), c("white","red2")),
  all_color_scales_equal = T,
  use_raster = T,
  raster_quality = 3,
  raster_device = "png",
  return_ht_list = F)
dev.off()


