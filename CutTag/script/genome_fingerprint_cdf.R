
library(tidyverse)

setwd("D:\\EED_CutTag\\cel_08-25\\result\\multiqc")

raw <- read_delim('multiqc/BAM_fingerprint.txt', delim="\t", skip=1) %>% 
  janitor::clean_names() %>% 
  glimpse()

data <- raw %>% 
  mutate(bin_id=row_number()) |> 
  pivot_longer(!bin_id, values_to='coverage') %>% 
  # separate_wider_regex(name, patterns=c(sample="[a-zA-Z]+", "\\d", "[-_]", rep="\\d+")) |>    # mice
  # separate_wider_regex(name, patterns=c(sample=".*", "[-_]", rep="\\d+")) |>     #  cel
  separate_wider_regex(name, patterns=c(sample="[a-zA-Z]+", "[-_]", rep="\\d+")) |>    # fly
  reframe(coverage=mean(coverage), .by=c(sample, bin_id)) |>    # merge by sample
  dplyr::select(!bin_id) |>     # discard old bin id
  print()
# > data
# # A tibble: 2,937,774 × 2
# sample coverage
# <fct>     <dbl>
#   1 d1       36164.
# 2 d5       14174 
# 3 l         4081 
# 4 m        20500.
# 5 ctrl     25026 
# 6 nmn      33852 
# 7 d1       47125 
# 8 d5       16618.
# 9 l         5040.
# 10 m        27281 

#  |> unique()
df <- data %>% 
  # mutate(sample=factor(sample, c("c","n","on" ))) %>%  # ms
  # mutate(sample=factor(sample, c("d1","d5","l","m","ctrl","nmn" ))) %>%   # cel
  print()

common_pal <- function() {
  list(
    geom_line(linewidth = 1),
    scale_color_brewer(palette = 'Paired'),
    labs(x = "Fraction of genome bins (ordered by descending coverage)",
         y = "Cumulative Fraction of Reads",
         color = "Sample"),
    theme_bw(base_size = 8)
  )
}


# df %>% 
#   ggplot(aes(coverage, color = sample)) +
#   stat_ecdf() +
#   common_pal() +  
#   coord_cartesian(xlim=c(0, 0.25e6))


df1 <- df %>% 
  group_by(sample) %>%
  # arrange(coverage) %>%
  arrange(desc(coverage)) %>%
  mutate(cumsum_reads = cumsum(coverage),
         total_reads = sum(coverage),
         cumfrac_reads = cumsum_reads / total_reads,
         frac_genome = row_number() / n()) %>%   # Fraction of bins
  ungroup() %>% 
  print()
# df$sample %>% unique()

# ks_result_1 <- ks.test(
#   filter(df1, sample=='d1') %>% pull(coverage),
#   filter(df1, sample=='d5') %>% pull(coverage)
#   ) %>% print()

p1 <- df1 %>% 
  dplyr::slice(seq(1, n(), by = 200)) %>% # reduce data point
  ggplot(aes(frac_genome, cumfrac_reads, color = sample)) +
  common_pal()
p1
ggsave(str_glue("{Sys.Date()}_cdf.pdf"), p1, width=5, height=3)

df1 %>% 
  dplyr::slice(seq(1, n(), by = 200)) %>% # reduce data point
  rio::export(str_glue("{Sys.Date()}_cdf.xlsx"))



df2 <- df |> 
  # set coverage bin width = 200
  # mutate(coverage_bin=cut_width(coverage, width=0.001, boundary=0, closed="left"), .by=sample) %>%
  mutate(coverage_bin=cut_width(coverage, width=10, boundary=, closed="left"), .by=sample) %>%
  mutate(coverage_bin=as.numeric(coverage_bin)) %>%
  reframe(coverage_bin_count=n(), .by=c(sample, coverage_bin)) %>% 
  arrange(coverage_bin) %>%          # Sort by coverage (high to low)
  mutate(cumsum_bin = cumsum(coverage_bin_count),
         total_bin = sum(coverage_bin_count),
         cumfrac_bin = cumsum_bin / total_bin,
         frac_bin = row_number() / n(), 
         .by=c(sample)) %>%   # Fraction of bins
  glimpse()
# df$sample %>% unique()

# ks_result_2 <- ks.test(
#   filter(df2, sample=='d1') %>% pull(cumfrac_bin),
#   filter(df2, sample=='d5') %>% pull(cumfrac_bin)
# ) %>% print()

p2 <- df2 %>% 
  # slice(seq(1, n(), by = 200)) %>% # reduce data point
  ggplot(aes(frac_bin, cumfrac_bin, color = sample)) +
  common_pal()
p2
ggsave(str_glue("{Sys.Date()}_cdf_bins.pdf"), p2, width=5, height=3)

df2 %>% 
  # slice(seq(1, n(), by = 200)) %>% # reduce data point
  rio::export(str_glue("{Sys.Date()}_cdf_bins.xlsx"))



df3 <- df %>% 
  filter(coverage < 10000) %>%
  # set coverage bin width = 200
  # mutate(coverage_bin=cut_width(coverage, width=0.001, boundary=0, closed="left"), .by=sample) %>%
  mutate(coverage_bin=cut_width(coverage, width=10, boundary=, closed="left"), .by=sample) %>%
  mutate(coverage_bin=as.numeric(coverage_bin)) %>%
  reframe(coverage_bin_count=n(), .by=c(sample, coverage_bin)) %>% 
  arrange(coverage_bin) %>%          # Sort by coverage (high to low)
  mutate(cumsum_bin = cumsum(coverage_bin_count),
         total_bin = sum(coverage_bin_count),
         cumfrac_bin = cumsum_bin / total_bin,
         frac_bin = row_number() / n(), 
         .by=c(sample)) %>%   # Fraction of bins
  glimpse()
# df$sample %>% unique()

# ks_result_3 <- ks.test(
#   filter(df3, sample=='dmso') %>% pull(cumfrac_bin),
#   filter(df3, sample=='imr05') %>% pull(cumfrac_bin)
# ) %>% print()

p3 <- df3 %>% 
  # slice(seq(1, n(), by = 200)) %>% # reduce data point
  ggplot(aes(frac_bin, cumfrac_bin, color = sample)) +
  common_pal() +
  labs(y='Cumulative Fraction of Bins')
p3
ggsave(str_glue("{Sys.Date()}_cdf_bins_cutoff.pdf"), p, width=5, height=3)

