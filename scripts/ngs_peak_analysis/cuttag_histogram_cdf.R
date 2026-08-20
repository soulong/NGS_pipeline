
library(tidyverse)
library(RSQLite)
library(ggpubr)

# setwd("C:/Users/haohe/Documents/GitHub")
# setwd("D:/LGLab/Project_NMN/DAPI_analysis_2025-12-10")
rstudioapi::getActiveDocumentContext()$path %>% 
  dirname() %>% setwd()
getwd()

# use_db <- "hist_100_60000_quantile_50_logspace" 
# use_db <-  "hist_100_8000_quantile_50_linear"
# use_db <- "hist_100_60000_quantile_100_linear"
# use_db <- "hist_150_65535_mean_50_linear"
# use_db <- "hist_150_65535_mean_50_sqrt"
# use_db <- "hist_150_65535_quantile_50_sqrt"
# use_db <- "hist_150_65535_mean_50_geomspace"
# use_db <- "hist_150_20000_mean_200_linear"
# use_db <- "hist_150_15000_mean_100_linear"
# use_db <- "hist_150_20000_mean_50_linear"


for(use_db in c(
  # "hist_100_60000_quantile_50_logspace",
  # "hist_100_8000_quantile_50_linear",
  # "hist_100_60000_quantile_100_linear",
  # "hist_150_65535_mean_50_linear",
  # "hist_150_65535_mean_50_sqrt",
  # "hist_150_65535_quantile_50_sqrt",
  "hist_150_65535_mean_50_geomspace",
  # "hist_150_12000_mean_50_geomspace",
  # "hist_150_20000_mean_200_linear",
  # "hist_150_15000_mean_100_linear",
  # "hist_150_20000_mean_50_linear"
)) {
  
  db_files <- list.files(".", pattern=".db", recursive=T, full.names=T) %>% 
    str_subset(use_db) %>% print()
  
  read_db <- function(f, table_name="cell") {
    conn <- DBI::dbConnect(RSQLite::SQLite(), f)
    df <- tbl(conn, table_name) %>% collect()
    DBI::dbDisconnect(conn)
    return(df)
  }
  
  
  raw_list <- map(db_files, read_db, table_name="cell")
  
  raw <- list_rbind(raw_list) %>% glimpse()
  # count(raw, directory)
  
  hist <- raw %>% 
    mutate(
      # Extract clean directory name
      directory=basename(dirname(directory)) %>% str_remove_all("__20.*"),
      # Compute bin edge mean
      bin_edge_mean=(bin_edge_low + bin_edge_high) / 2,
      # Cell type
      celltype=ifelse(str_detect(directory, fixed("MEF")), "MEF", "IMR90"),
      # Assign group — NO reference to 'group' inside!
      group=case_when(
        directory == "20251114 MEF KO H3K27me3 IF" ~
          case_when(
            well == "B6" ~ "DMSO",
            well == "C6" ~ "NMN_5mM",
            well == "D6" ~ "NMN_10mM",
            TRUE ~ NA_character_
          ),
        directory == "LH 20231107 IMR90 NMN treatment H3K27me3" ~
          case_when(
            well == "B1" ~ "DMSO",
            well %in% c("C1") ~ "NMN",
            TRUE ~ NA_character_
          ),
        directory == "LH 20231107 IMR90 young H3K27me3" ~
          ifelse(well == "D4", "DMSO", NA_character_),
        directory == "20251201 MEF WT H3K27me3 IF NMN 226" ~
          case_when(
            well == "C2" ~ "DMSO",
            well == "C3" ~ "NMN_5mM",
            well == "C4" ~ "NMN_10mM",
            TRUE ~ "EED26"
          ),
        directory == "20251211 MEF WT NMN EED226 H3K27me3 IF" ~
          case_when(
            well == "D2" ~ "DMSO",
            well == "D3" ~ "NMN_3mM",
            well == "D4" ~ "NMN_10mM",
            TRUE ~ "EED26"
          ),
        directory == "20251214 IMR90 H3K27me3 IF NMN 226" ~
          case_when(
            well == "B2"  ~ "DMSO",
            well == "B3"  ~ "EED226_1uM",
            well == "B4"  ~ "EED226_3uM",
            well == "B5"  ~ "EED226_10uM",
            well == "C1"  ~ "NMN_1mM",
            well == "D2" ~ "NMN_5mM",
            well == "C3"  ~ "NMN_10mM",
            TRUE ~ NA_character_
          ),
        directory == "20251222 old plate IMR90 NMN" ~
          case_when(
            well == "B4" ~ "DMSO",
            well == "A5" ~ "NMN_0.3mM",
            well == "B6" ~ "NMN_0.6mM",
            well == "C5" ~ "NMN_1mM",
            TRUE ~ NA_character_
          ),
        directory == "20251222 old plate2 IMR90 NMN" ~
          case_when(
            well == "C5" ~ "DMSO",
            well == "C3" ~ "NMN_0.5mM",
            well == "C4" ~ "NMN_1mM",
            TRUE ~ "NAM_1mM"
          ),
        directory == "20251223 MEF H3K27me3 IF NMN" ~
          case_when(
            well == "C2" ~ "DMSO",
            well == "C3" ~ "NMN_0.5mM",
            TRUE ~ NA_character_
          ),
        TRUE ~ NA_character_
      )
      
    ) %>%
    select(-timepoint, -stack)
  
  
  glimpse(hist)
  hist_filter <- hist %>% 
    # filter(directory==unique(hist_selected$directory)) %>%
    # filter(between(cell_mean_intensity, 300, 6000)) %>%
    # filter(between(expm1(cell_mean_intensity), 300, 6000)) %>%
    # filter(directory=="20251114 MEF KO H3K27me3 IF") %>%
    # filter(well=="B6") %>% 
    # filter(field %in% c("1","2")) %>% 
    filter(between(total_pixels, 4500, 9000))
  
  # distinct(hist_filter, bin_index, bin_edge_mean)
  hist_filter %>% 
    distinct(directory, group, well, field, cell_id) %>% 
    dplyr::count(directory, group, field, well) %>% 
    dplyr::count(directory, group, well) %>% 
    rio::export(str_glue("{Sys.Date()}_field_count_per_group.xlsx"))
  

  hist_cdf <- hist_filter %>%
    reframe(bin_mean=mean(bin_mean),
            bin_edge_mean=mean(bin_edge_mean),
            bin_count=sum(bin_count),
            .by=c(directory, group, celltype, bin_index)) %>%
    group_by(directory, group, celltype) %>%
    # Ensure bins are in order
    arrange(bin_index, .by_group=T) %>%
    # weighted CDF for binned data
    mutate(cdf=cumsum(bin_count) / sum(bin_count) ) %>%
    ungroup() %>%
    arrange(celltype) %>%
    split(., .[['directory']])
  p2 <- hist_cdf %>%
    imap(\(x, y) x %>%
           ggplot(aes(bin_index, cdf, color=group)) +
           # geom_step(direction="hv") +
           geom_line() +
           labs(x="Intensity Bins", y="Cumulative Fraction of Pixels", title=y) +
           theme_bw(8) +
           theme(panel.grid=element_blank())
    ) %>%
    patchwork::wrap_plots(ncol=3)
  # print(p2)
  ggsave(str_glue("{Sys.Date()}_cdf_{use_db}.pdf"),
         p2, width=12, heigh=6)

  # statistics
  stats <- hist_filter %>%
    split(., .[['directory']]) %>%
    map(function(x) {
      if("DMSO" %in% x[["group"]] & length(unique(x[["group"]])) > 1) {
        dmso <- filter(x, group == "DMSO") %>% pull(bin_count)
        setdiff(unique(x[["group"]]), "DMSO") %>%
          set_names(., nm=.) %>%
          map(\(y) filter(x, group == y) %>% pull(bin_count) %>%
                ks.test(., dmso) %>%
                broom::tidy()) %>%
          list_rbind(names_to="group")
      } else return(NULL)
    }) %>%
    list_rbind(names_to="directory") %>%
    dplyr::rename(ks_distance=statistic)
  writexl::write_xlsx(stats, str_glue("{Sys.Date()}_cdf_{use_db}.xlsx"))
}

