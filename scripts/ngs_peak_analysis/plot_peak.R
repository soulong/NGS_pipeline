
# visualize peak based on peak center

rstudioapi::getActiveDocumentContext()$path |>
  dirname() |>
  setwd()

library(tidyverse)
library(rtracklayer)
library(EnrichedHeatmap)


if(F) {
  c("HUVEC_5mM","IMR_0.5mM","MEF_5mM")
  setwd("MEF_5mM")
  setwd("..")
  getwd()
}


## 1.use consensus peak --------------------------------
# load peak
peak <- import("consensus_peaks.bed")

# set peak to 6kb window centered on the middle 
peak_resized <- resize(peak, width=6000, fix = "center")

# set the center for resized peak
peak_center <- resize(peak_resized, width=1, fix = "center")

# read bigwig
bw_files <- list.files(".", ".(bigwig)|(bw)") %>%
  set_names(., nm=basename(.))
signals <- map(bw_files, \(x) import(x, format = "bw", which=peak_resized))
  
# normalization for resized reads
signals_mat <- map(signals, 
                   \(x) normalizeToMatrix(x, peak_center, 
                                          value_column="score",
                                          # empty_value = 0,
                                          # keep = c(0, 1),
                                          mean_mode="w0", w=40, 
                                          extend = 1000))
# check quantile range outliers for the matrix
walk(signals_mat, \(x) print(quantile(x, probs = c(0.005, 0.5,0.995))))


# mapping haetmap colors and set the range
# col_fun<- circlize::colorRamp2(c(0, 1.3), c("white", "red"))
col_fun <- circlize::colorRamp2(quantile(signals_mat[[1]], c(0.01, 0.99)), c("#3288bd", "#d53e4f"))
# set haetmap parameters
ht_global_opt(heatmap_column_names_gp = gpar(fontsize = 14),
              heatmap_legend_title_gp = gpar(fontsize = 14),
              heatmap_legend_labels_gp = gpar(fontsize = 14),
              ADD = TRUE)
ht_list <- imap(signals_mat, 
                \(x,y) EnrichedHeatmap(x, 
                                       name = y, #axis_name_rot = 45,
                                       column_title = y, 
                                       use_raster = T, pos_line = F,
                                       col = col_fun, 
                                       top_annotation = HeatmapAnnotation(
                                         enriched = anno_enriched(
                                           #ylim = c(0.1, 0.6),
                                           axis_param = list(
                                             # at = c(0.1, 0.3, 0.5),
                                             # labels = c("0", "0.25", "0.5"),
                                             side = "left",
                                             facing = "outside"
                                           ))),
                                       axis_name = c("-1kb", "center", "1kb")))
# draw the heatmap
pdf(str_glue("{Sys.Date()}_heatmap_peak_center_1kb.pdf"), 
    width = 5, height = 8)
draw(ht_list[[1]] + ht_list[[2]])
# walk(ht_list, \(x) print(x))
dev.off()
while(dev.cur() != 1) dev.off()




## 2.use separated peaks -------------------------------
# read peak
# peak_files <- list.files(".", "narrowPeak", recursive = F, full.names = F)
peak_files <- list.files(".", ".narrowPeak") %>%
  set_names(., nm=basename(.))
peaks <- map(peak_files, import)

# set peaks to 6kb window centered on the middle 
peaks_resized <- map(peaks, \(x) resize(x, width=6000, fix = "center"))

# read bigwig
bw_files <- list.files(".", ".(bigwig)|(bw)") %>%
  set_names(., nm=basename(.))
signals <- map2(bw_files, peaks_resized, \(x,y) import(x, format = "bw", which=y))

# set the center for resized peaks
peaks_center <- map(peaks_resized, \(x) resize(x, width=1, fix = "center"))

# normalization for resized reads
signals_mat <- map2(signals, peaks_center, 
                    \(x,y) normalizeToMatrix(x, y, 
                                             value_column="score",
                                             # empty_value = 0,
                                             # keep = c(0, 1),
                                             mean_mode="w0", w=40, 
                                             extend = 1000))

# check quantile range outliers for the matrix
walk(signals_mat, \(x) print(quantile(x, probs = c(0.005, 0.5,0.995))))



# mapping haetmap colors and set the range
# col_fun<- circlize::colorRamp2(c(0, 1.2), c("white", "red"))
cols_fun <- map(signals_mat, 
               \(x) circlize::colorRamp2(quantile(x, c(0.005, 0.995)), c("#e41a1c", "#386cb0")))
# set haetmap parameters
ht_global_opt(heatmap_column_names_gp = gpar(fontsize = 14),
              heatmap_legend_title_gp = gpar(fontsize = 14),
              heatmap_legend_labels_gp = gpar(fontsize = 14),
              ADD = TRUE)
ht_list <- imap(signals_mat, 
                \(x,y) EnrichedHeatmap(x, 
                                       name = y, #axis_name_rot = 45,
                                       column_title = y, 
                                       use_raster = T, pos_line = F,
                                       col = cols_fun[[y]], 
                                       axis_name = c("-3kb", "center", "3kb")))
# draw the heatmap
pdf(str_glue("{Sys.Date()}_enrichedHeatmap.pdf"), width = 3, height = 7)
walk(ht_list, \(x) print(x))
dev.off()
while(dev.cur() != 1) dev.off()






