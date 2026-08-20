

rstudioapi::getActiveDocumentContext()$path |>
  dirname() |> setwd()
getwd()

library(tidyverse)
library(ChIPseeker)
library(ChIPQC)
library(DiffBind)
library(yaml)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggrastr)
library(cowplot)
library(EnrichedHeatmap)
library(profileplyr)

library(furrr)
plan(multisession, workers = 4)


config <- yaml::read_yaml("bin/config.yml")
bam_dir <- config$bam_dir # "bam"
peak_dir <-  config$peak_dir # "peak"
species <- config$species

peak_dir <- "peak" # this is call masc2 using clean.bam
print(str_glue("{bam_dir}, {peak_dir}, {species}"))



###################### bam qc ###########################
bam_files <- list.files(bam_dir, ".bam$", full.names = T) %>%
  str_subset(".clean", negate = F) %>%
  set_names(., map_chr(., ~ basename(.x) %>% str_replace_all(".bam", "")))
print(bam_files)




###################### peak vis [ChIPseeker] ###########################
peak_files <- list.files(peak_dir, ".narrowPeak$", full.names = T) %>%
  set_names(., map_chr(., ~ basename(.x) %>% str_replace_all("_peaks.narrowPeak", "")))
print(peak_files)

# remove fix samples [mssing related ctrl sample]
# peak_files <- str_subset(peak_files, "_Fix_", negate = T)

# metedata info
c1 <- str_c('Y', 1:28)
c2 <- c(rep('DMSO', 12), rep('ET0825', 8), rep('ET0823', 8)) %>%
  factor(levels = c('DMSO', 'ET0825', 'ET0823'))
c3 <- c(rep('YAP', 3), rep('TAZ', 3), rep('panTEAD', 3), rep('H3K27AC', 3),
        rep('YAP', 2), rep('TAZ', 2), rep('panTEAD', 2), rep('H3K27AC', 2),
        rep('YAP', 2), rep('TAZ', 2), rep('panTEAD', 2), rep('H3K27AC', 2)) %>%
  factor(levels = c('YAP', 'TAZ', 'panTEAD', 'H3K27AC'))
c4 <- c(rep(19, 6), rep(15, 3), rep(12, 3),
        rep(19, 4), rep(15, 2), rep(12, 2),
        rep(19, 4), rep(15, 2), rep(12, 2)) %>%
  factor()

info <- tibble(Sample = c1, Target = c3, Treatment = c2, PCR_Cycle = c4) %>%
  mutate(Condition = str_c(Target, Treatment, sep = '_'),
         Label = str_c(Condition, Sample, sep = '_')) %>%
  mutate(peak_path = str_c(Sample, '.clean.bam', sep = ''),
         bam_path = str_c(Sample, '_peaks.narrowPeak', sep = ''))

  left_join(enframe(peak_files, 'Sample', 'peak_path'), by='Sample') %>%
  left_join(enframe(bam_files, 'Sample', 'bam_path'), by='Sample')

# reorder files according to info
peak_files <- info
names(peak_files) <- info$Label

# get peak grange list
peak_list <- future_map(peak_files, readPeakFile) # using chipseeker
# peaks <- GRangesList(peak_list)

# change seq style
seqlevelsStyle(peak_list[[1]])
peak_list <- map(peak_list, function(x) {seqlevelsStyle(x) <- "UCSC"; return(x)} )
seqlevelsStyle(peak_list[[1]])
seqlevels(peak_list[[1]])

# keep autosome
peak_list <- future_map(peak_list, function(x) 
  { keepStandardChromosomes(x, pruning.mode="coarse") %>% 
      dropSeqlevels("chrM", pruning.mode="coarse")} )
seqlevels(peak_list[[1]])


# set txdb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene # seqlevelsStyle(txdb) <- "UCSC"


# get tss region
tss_tagmatrix_list <- list()
tss <- getBioRegion(TxDb=txdb, 
                    upstream=3000, downstream=3000, 
                    type="start_site", by="gene")
for(i in unique(info$Group)) {
  cat("\n##### get tagMatrix on group: ", i, "\n")
  subset <- filter(info, Group==i) %>% 
    pull(Label) %>% peak_list[.]
  tag_matrix <- future_map(subset, function(x) {getTagMatrix(x, windows=tss, weightCol="V5")})
  tss_tagmatrix_list[[i]] <- tag_matrix
}

# plot tss region
for(i in names(tss_tagmatrix_list)) {
  cat("\n##### plot profile on group: ", i, "\n")
  # peak average profile
  p <- plotPeakProf(tss_tagmatrix_list[[i]], conf=0.95, facet="column", free_y=T)
  ggsave2(str_glue("peak_profile_tss_{fs::path_sanitize(i)}.pdf"), p,
          width = length(tss_tagmatrix_list[[i]])*4, height = 4, limitsize = F)
  # peak heatmap
  pdf(str_glue("peak_heatmap_tss_{fs::path_sanitize(i)}.pdf"), width = length(tss_tagmatrix_list[[i]])*2.5, height = 5)
  tagHeatmap(tss_tagmatrix_list[[i]], xlim=c(-3000, 3000))
  dev.off()
}
while(dev.cur()!=1) dev.off()


if(T) {
  # body region
  for(i in unique(info$Group)) {
    cat("\n##### plot profile body on group: ", i, "\n")
    subset <- filter(info, Group==i) %>% 
      pull(Label) %>% peak_list[.]
    p <- plotPeakProf2(subset, conf=0.95, weightCol="V5", 
                       facet="column", free_y=T,
                       upstream=rel(0.2), downstream=rel(0.2),
                       TxDb=txdb, type="body", by="gene", nbin=1000)
    ggsave2(str_glue("peak_profile_body_{fs::path_sanitize(i)}.pdf"), p,
            width = length(subset)*4, height = 4, limitsize = F) }
}




###################### peak annotation [ChIPseeker] ###########################
# reduce peaks among group
# info_filter <- filter(info, !factor %in% c("NF"))
peak_list_reduced <- list()
for(i in unique(info$Group)) {
  cat('\n##### make consensus peaks of:', i)
  peak_list_reduced[[i]] <- filter(info, Group==i) %>%
    pull(Label) %>% 
    peak_list[.] %>%
    GRangesList() %>%
    GenomicRanges::coverage() %>%
    # keep regions of two sample overlap
    IRanges::slice(lower=2, rangesOnly=F) %>%
    GRanges() %>%
    # merge peaks of gap < 30
    GenomicRanges::reduce(min.gapwidth=30)
  }

names(peak_list_reduced)
# plot overlap
pdf(str_glue("peak_annotation_peak_overlap_DMSO.pdf"), width = 6, height = 6)
ChIPseeker::vennplot(peak_list_reduced[1:3], by="Vennerable")
dev.off()

# plot overlap
pdf(str_glue("peak_annotation_peak_overlap_ET0825.pdf"), width = 6, height = 6)
ChIPseeker::vennplot(peak_list_reduced[5:7], by="Vennerable")
dev.off()

# plot overlap
pdf(str_glue("peak_annotation_peak_overlap_ET0823.pdf"), width = 6, height = 6)
ChIPseeker::vennplot(peak_list_reduced[9:11], by="Vennerable")
dev.off()

# plot overlap
pdf(str_glue("peak_annotation_peak_overlap_YAP.pdf"), width = 6, height = 6)
ChIPseeker::vennplot(peak_list_reduced[c(1,5,9)], by="Vennerable")
dev.off()

# plot overlap
pdf(str_glue("peak_annotation_peak_overlap_TAZ.pdf"), width = 6, height = 6)
ChIPseeker::vennplot(peak_list_reduced[c(2,6,10)], by="Vennerable")
dev.off()

# plot overlap
pdf(str_glue("peak_annotation_peak_overlap_panTEAD.pdf"), width = 6, height = 6)
ChIPseeker::vennplot(peak_list_reduced[c(3,7,11)], by="Vennerable")
dev.off()

# plot overlap
pdf(str_glue("peak_annotation_peak_overlap_H3K27AC.pdf"), width = 6, height = 6)
ChIPseeker::vennplot(peak_list_reduced[c(4,8,12)], by="Vennerable")
dev.off()

while(dev.cur()!=1) dev.off()


# annotate peaks
peak_list_reduced_anno <- map(peak_list_reduced,
                              ~ annotatePeak(.x,
                                             tssRegion = c(-3000, 3000),
                                             TxDb = txdb,
                                             annoDb = "org.Hs.eg.db",
                                             level = "gene"))
# peak annoated genes overlap
pdf(str_glue("peak_annotation_gene_overlap.pdf"), width = 6, height = 6)
gene_overlap <- map(peak_list_reduced_anno, function(i) as.data.frame(i)$SYMBOL)
# vennplot(gene_overlap, by="Vennerable")
vennplot(gene_overlap[3:6], by="Vennerable")
dev.off()
while(dev.cur()!=1) dev.off()

# plot annotation
p1 <- plotAnnoBar(peak_list_reduced_anno) %>% 
  ggplotify::as.ggplot()
p2 <- plotDistToTSS(peak_list_reduced_anno, title = "Distance to TSS") %>% 
  ggplotify::as.ggplot()
ggsave2("peak_annotation_profile.pdf", plot_grid(p1, p2),
        width=12, height = length(peak_list_reduced_anno)*0.5)

# write mapped genes
peak_list_reduced_anno_df <- map(peak_list_reduced_anno, as.data.frame)
writexl::write_xlsx(peak_list_reduced_anno_df, "peak_annotation_genes.xlsx")




###################### peak vis ###########################






###################### peak diff [DiffBind] ###########################
use_bam_suffix <- ".clean.bam"
blacklist <- DBA_BLACKLIST_GRCH38

# construct sample_info from peak info files
sample_info <- info %>% 
  dplyr::transmute(SampleID=Label, 
                   bamReads=file.path(bam_dir, str_c(Sample, {{use_bam_suffix}})), 
                   Peaks=Path, PeakCaller="narrow", 
                   Condition=Group, Treatment=Treatment, Factor=Target, Replicate=PCR_Cycle,
                   ControlID=NA, bamControl=NA)

# sample_file <- "sample_info.csv"
# if(file.exists(sample_file)) {
#   print("file exist, please modify it")
# } else {
#   write_csv(sample_info, sample_file) 
#   print("please modify [sample_info.csv]") }
# 
# ## edit the sample_info.csv to fullfill requirements
# sample_info <- read_csv("sample_info.csv")
# sample_info


## construct dba object
cut_tag <- sample_info %>%
  # filter(Condition %in% c("H3K27ac_ctrl","H3K27ac_ET825")) %>%
  dba(sampleSheet=., minOverlap=2) %>%
  dba.blacklist(blacklist=blacklist)
## count matrix
cut_tag <- dba.count(cut_tag, score=DBA_SCORE_NORMALIZED, 
                     bRemoveDuplicates=T)

if(T) {
  write_rds(cut_tag, "cut_tag.rds", compress = "bz2")
  cut_tag <- read_rds("cut_tag.rds")
  }


pdf("cut_tag_qc.pdf")
dba.plotPCA(cut_tag, attributes=DBA_CONDITION, score=DBA_SCORE_NORMALIZED)
dba.plotPCA(cut_tag, attributes=DBA_TREATMENT, score=DBA_SCORE_NORMALIZED)
dba.plotPCA(cut_tag, attributes=DBA_FACTOR, score=DBA_SCORE_NORMALIZED)
# heatmap correlation
dba.plotHeatmap(cut_tag, DBA_ID)
# heatmap peak counts
dba.plotHeatmap(cut_tag, DBA_ID, correlations=F, scale="row")
dev.off()


# profile
sample_profiles <- dba.plotProfile(cut_tag, normalize=TRUE,
                                   labels=DBA_CONDITION, merge=DBA_ID)
write_rds(sample_profiles, "cut_tag_profiles.rds", compress = "bz2")

pdf("cut_tag_profile.pdf", width = 12, height = 7)
dba.plotProfile(sample_profiles)
dev.off()

dba.save(cut_tag, "cut_tag")

dba.show(cut_tag) %>% as_tibble() %>%
  mutate(peakReads = round(Reads * FRiP), outReads = Reads - peakReads) %>%
  pivot_longer(c(outReads, peakReads), names_to = "type", values_to = "value") %>%
  ggplot(aes(ID, value, fill = type)) +
  geom_col() +
  labs(x=NULL, y="Read Counts") +
  ggpubr::theme_pubr(x.text.angle = 45)
ggsave("cut_tag_frip.pdf")






###########################################################
# ## subset yap analysis
cut_tag_yap <- cut_tag %>%
  dba.contrast(contrast=c('Condition', 'YAP_ET0825', 'YAP_DMSO')) %>%
  dba.contrast(contrast=c('Condition', 'YAP_ET0823', 'YAP_DMSO')) %>%
  # dba.contrast(contrast=c('Treatment', 'ET0825', 'ET0823')) %>%
  dba.analyze()
cut_tag_yap_ET0825 <- dba.report(cut_tag_yap, contrast=1, th=0.05, fold = 1, bCalled=T, bCounts=T, bCalledDetail=F)

cut_tag_yap_ET0825_profiles <- dba.plotProfile(cut_tag, sites=cut_tag_yap_ET0825, labels=DBA_CONDITION, merge=DBA_ID)
dba.plotProfile(cut_tag_yap_ET0825_profiles)

cut_tag_yap_ET0823 <- dba.report(cut_tag_yap, contrast=2, th=0.05, fold = 1, bCalled=T, bCounts=T, bCalledDetail=F)
writexl::write_xlsx(list(ET0825_DMSO=cut_tag_yap_ET0825 %>% as.data.frame(),
                         ET0823_DMSO=cut_tag_yap_ET0823 %>% as.data.frame()),
                    'cut_tag_yap_db_analysis.xlsx')

library(profileplyr)
load(system.file("extdata", "chr21_test_data.RData", package = "EnrichedHeatmap"))
ls()
genes


# x <- bam_files[1] %>% set_names(nm=NULL)
# cov <- coverage(x) %>% GRanges()
# seqlevelsStyle(cov) <- 'UCSC'

# genes <- genes(txdb)
# tss <- promoters(genes, 0, 1)
# tss <- getBioRegion(TxDb=txdb,
#                     upstream=1000, downstream=1000,
#                     type="start_site", by="gene")
# tss <- resize(tss, 1, fix='center')

# tss <- tss[1:5000]
# mat <- EnrichedHeatmap::normalizeToMatrix(cov, tss, extend=3000, value_column='score', 
#                                           mean_mode='w0', smooth=T, keep=c(0, 0.99), verbose=T)
# EnrichedHeatmap::EnrichedHeatmap(mat, col=c('white', 'red2'))

soGGi::regionPlot()
rtracklayer::export.bed(cut_tag_yap_ET0825, 'cut_tag_yap_ET0825.bed')
profile <- profileplyr::BamBigwig_to_chipProfile(bam_files, 'cut_tag_yap_ET0825.bed', 
                                                 samplename=
                                                 format='bam', paired=T,
                                                 style='percentOfRegion', nOfWindows=100,
                                                 # style='point', bin_size=20,
                                                 quant_params=MulticoreParam(8)
                                                 )


profileplyr::generateEnrichedHeatmap()
profileplyr::generateProfilePlot()
pdf('cut_tag_yap_db_ET0825_profile.pdf')

dba.plotProfile(cut_tag, 
                samples=list(cut_tag$masks$YAP, cut_tag$masks$TAZ, cut_tag$masks$panTEAD, cut_tag$masks$H3K27AC),
                sites=cut_tag_yap_ET0825,
                merge=DBA_ID, doPlot=T)
while(dev.cur()!=1) dev.off()

pdf("cut_tag_yap_db_analysis.pdf")

dba.plotMA(cut_tag_yap, contrast=1)
dba.plotMA(cut_tag_yap, contrast=2)

dba.plotVolcano(cut_tag_yap, contrast=1, fold = 1, dotSize=2)
dba.plotVolcano(cut_tag_yap, contrast=2, fold = 1, dotSize=2)

dba.plotHeatmap(cut_tag_yap, attributes=DBA_CONDITION, maxSites=1000, 
                report = cut_tag_yap_ET0825,
                contrast=1,
                correlations=F, scale="row")
# dba.plotHeatmap(cut_tag_yap, attributes=DBA_CONDITION, maxSites=1000, 
#                 report = cut_tag_yap_ET0823,
#                 contrast=2,
#                 correlations=F, scale="row")
while(dev.cur()!=1) dev.off()



###########################################################
# ## compare only yap
yap <- sample_info %>%
  filter(Factor %in% c("YAP")) %>%
  dba(sampleSheet=.,
      config=data.frame(th=0.05, RunParallel=TRUE, minQCth=10)) %>%
  dba.blacklist(blacklist=blacklist)

# calculate consensus peak within each group
yap_consensus <- dba.peakset(yap, consensus = DBA_TREATMENT, minOverlap=2)
yap_consensus <- dba(yap_consensus,
                     mask=yap_consensus$masks$Consensus,
                     minOverlap=1)
# overall consensus peak
consensus_peaks <- dba.peakset(yap_consensus, bRetrieve=TRUE)
# binding matrix
yap <- dba.count(yap, peaks=consensus_peaks)


yap_consensus_overlap <- dba.overlap(yap_consensus, yap_consensus$masks$Consensus)
writexl::write_xlsx(list(only_DMSO=yap_consensus_overlap$onlyA %>% as.data.frame(),
                         only_ET0825=yap_consensus_overlap$onlyB %>% as.data.frame(),
                         only_ET0823=yap_consensus_overlap$onlyC %>% as.data.frame(),
                         not_DMSO=yap_consensus_overlap$notA %>% as.data.frame(),
                         not_ET0825=yap_consensus_overlap$notB %>% as.data.frame(),
                         not_ET0823=yap_consensus_overlap$notC %>% as.data.frame(),
                         in_ALL=yap_consensus_overlap$inAll %>% as.data.frame()),
                    'yap_peak_overlap.xlsx')
pdf("yap_qc.pdf")
dba.plotVenn(yap_consensus, yap_consensus$masks$Consensus)
dba.plotPCA(yap, attributes = DBA_TREATMENT)
dba.plotProfile(yap, labels=DBA_TREATMENT, merge=DBA_ID, doPlot=T)

dba.plotHeatmap(yap, attributes=DBA_TREATMENT, maxSites=1000, 
                correlations=T, scale="row")
dba.plotHeatmap(yap, attributes=DBA_TREATMENT, maxSites=1000, 
                correlations=F, scale="row")
while(dev.cur()!=1) dev.off()


if(T) {
  write_rds(yap, "yap.rds", compress = "bz2")
  yap <- read_rds("yap.rds")
}


## yap design
yap <- yap %>%
  dba.contrast(contrast=c('Treatment', 'ET0825', 'DMSO')) %>%
  dba.contrast(contrast=c('Treatment', 'ET0823', 'DMSO')) %>%
  # dba.contrast(contrast=c('Treatment', 'ET0825', 'ET0823')) %>%
  dba.analyze()

dba.show(yap, bContrast=TRUE)
yap_ET0825 <- dba.report(yap, contrast=1, th=1, bCalled=T, bCounts=T, bCalledDetail=F)
yap_ET0823 <- dba.report(yap, contrast=2, th=1, bCalled=T, bCounts=T, bCalledDetail=F)
writexl::write_xlsx(list(ET0825_DMSO=yap_ET0825 %>% as.data.frame(), 
                         ET0823_DMSO=yap_ET0823 %>% as.data.frame()),
                    'yap_db_analysis.xlsx')

pdf("yap_db_analysis.pdf")
dba.plotProfile(yap, labels=DBA_TREATMENT, merge=DBA_ID, doPlot=T)

dba.plotMA(yap, contrast=1)
dba.plotMA(yap, contrast=2)

dba.plotVolcano(yap, contrast=1, fold = 1, dotSize=2)
dba.plotVolcano(yap, contrast=2, fold = 1, dotSize=2)

dba.plotHeatmap(yap, attributes=DBA_CONDITION, maxSites=1000, 
                correlations=F, scale="row")
dba.plotHeatmap(yap, attributes=DBA_CONDITION, maxSites=1000, 
                report = yap_ET0825,
                contrast=1,
                correlations=F, scale="row")
dba.plotHeatmap(yap, attributes=DBA_CONDITION, maxSites=1000, 
                contrast=2,
                report = yap_ET0823,
                correlations=F, scale="row")
while(dev.cur()!=1) dev.off()













###################### track vis [trackViewer] ###########################
library(trackViewer)
wg_files <- list.files(bam_dir, ".bigwig$", full.names = T) %>%
  set_names(., map_chr(., ~ basename(.x) %>% str_replace_all(".rpgc.bigwig", "")))
print(wg_files)

# gene region
roi <- peak_list$Ctrl_YAP_Fix[1:2]

# change style
seqlevelsStyle(txdb) <- "NCBI"
seqlevelsStyle(txdb) <- "UCSC"

# build gene model
gene_model <- geneModelFromTxdb(txdb,
                                org.Hs.eg.db,
                                gr=GRanges(Rle("chr1"), IRanges(start=1, end=10000000), Rle("*")))
entrez_id <- AnnotationDbi::get("CCN1", org.Hs.egSYMBOL2EG)
gene_track <- geneTrack(entrez_id, txdb)[[1]]


grW <- parse2GRanges("chr1:85,578,761-85,585,950")
ids <- getGeneIDsFromTxDb(grW, txdb)
symbols <- AnnotationDbi::mget(ids, org.Hs.egSYMBOL)
gene_track <- geneTrack(ids, txdb, symbols, asList=F)



# which bams to show
f <- group_by(info, group) %>%
  slice_head()
fs <- str_c(bam_dir, "/", f$sample, ".bam") %>%
  set_names(nm=f$group)
print(fs)

genes <- c("CCN1", "CCN2") # ANKRD1
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

res <- list()
grs <- list()
for(g in genes) {
  cat("get gene:", g, "\n")
  entrez_id <- AnnotationDbi::get(g, org.Hs.egSYMBOL2EG)
  symbols <- AnnotationDbi::mget(entrez_id, org.Hs.egSYMBOL)
  gene_track <- geneTrack(entrez_id, txdb, symbols, asList=F)
  
  trackList <- list()
  trackList[["GRCh38"]] <- gene_track
  
  grW <- gene_track@dat
  seqlevelsStyle(grW) <- "NCBI"
  # import data
  for(f in names(fs)) {
    target <- importBam(fs[f], ranges=grW, pairs =T)
    target$dat <- coverageGR(target$dat)
    seqlevelsStyle(target) <- "UCSC"
    trackList[[f]] <- target
  }
  
  seqlevelsStyle(grW) <- "UCSC"
  grs[[g]] <- grW
  res[[g]] <- trackList
}


viewTracks(res$CCN2, gr=grs$CCN2, autoOptimizeStyle=T, smooth=F)

optSty <- optimizeStyle(res$CCN1)
viewTracks(trackList, gr=grW, viewerStyle = optSty$style)


getGeneIDsFromTxDb(grW, txdb)


seqlevelsStyle(grW) <- "NCBI"
target <- importBam(f, ranges=grW, pairs =T)
target$dat <- coverageGR(target$dat)
seqlevelsStyle(target) <- "UCSC"
seqlevelsStyle(grW) <- "UCSC"

trackList <- trackList(gene_track, target)
optSty <- optimizeStyle(trackList, theme="safe")

viewTracks(trackList, gr=grW)



viewerStyle <- trackViewerStyle()
# setTrackViewerStyleParam(viewerStyle, "margin", c(.1, .05, .02, .02))
vp <- viewTracks(trackList, 
                 gr=gr, viewerStyle=viewerStyle, 
                 autoOptimizeStyle=TRUE)






browseTracks(trackList)

vp <- viewTracks(trackList, 
                 gr=gr, viewerStyle=viewerStyle, 
                 autoOptimizeStyle=TRUE)



viewTracks(trackList(fox2), gr=gr, autoOptimizeStyle=TRUE, newpage=FALSE)





###########################################################
# chipseq qc
# chipqc <- ChIPQC(cut_tag, annotation="hg38")






