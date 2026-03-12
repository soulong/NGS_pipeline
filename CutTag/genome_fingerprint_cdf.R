
library(tidyverse)
setwd("E:/Project_NMN/IMR/result/05_multiqc")


raw <- read_delim('BAM_fingerprint.txt', delim="\t", skip=1) %>% 
  janitor::clean_names() %>% 
  glimpse()

df <- raw %>% 
  select(1:6) %>% 
  transmute(dmso=(imr0_1 + imr0_2 + imr0_3),
            imr05=(imr05_1 + imr05_2 + imr05_3)) %>%
  mutate(., mean=rowMeans(select(., 1:2))) %>%
  # mutate(dmso=dmso/sum(dmso), 
  #        IMR05=IMR05/sum(IMR05)) %>% 
  # filter(mean > 100) %>%
  # arrange(desc(mean)) %>%
  # filter(row_number() > 2000, row_number() < 10000) %>%
  # slice_max(mean, prop=0.9) %>%
  select(!mean) %>%
  
  pivot_longer(everything(), names_to='sample', values_to='coverage') %>% 
  group_by(sample) %>%
  # arrange(coverage) %>%
  arrange(desc(coverage)) %>%          # Sort by coverage (high to low)
  mutate(cumsum_reads = cumsum(coverage),
         total_reads = sum(coverage),
         cumfrac_reads = cumsum_reads / total_reads,
         frac_genome = row_number() / n()) %>%   # Fraction of bins
  glimpse()
# df$sample %>% unique()

ks_result <- ks.test(
  filter(df, sample=='dmso') %>% pull(coverage),
  filter(df, sample=='imr05') %>% pull(coverage)
  ) %>% print()


p <- df %>% 
  slice(seq(1, n(), by = 200)) %>% # reduce data point
  ggplot(aes(frac_genome, cumfrac_reads, color = sample)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Fraction of genome bins (ordered by descending coverage)",
    y = "Cumulative Fraction of Reads",
    color = "Sample"
  ) +
  theme_bw(base_size = 8)
p
ggsave(str_glue("{Sys.Date()}_cdf.pdf"), p, width=5, height=3)

# p +
#   coord_cartesian(xlim=c(0.001, 0.7), y=c(0.25, 1))
# ggsave(str_glue("{Sys.Date()}_cdf2.pdf"), width=5, height=3)







df <- raw %>% 
  select(1:6) %>% 
  
  transmute(dmso=(imr0_1 + imr0_2 + imr0_3),
            imr05=(imr05_1 + imr05_2 + imr05_3)) %>%

  mutate(., mean=rowMeans(select(., 1:2))) %>%
  # filter(mean < 5000) %>%
  select(!mean) %>%
  # mutate(dmso=dmso/sum(dmso),
  #        IMR05=IMR05/sum(IMR05)) %>%
  pivot_longer(everything(), names_to='sample', values_to='coverage') %>% 
  # set coverage bin width = 200
  # mutate(coverage_bin=cut_width(coverage, width=0.001, boundary=0, closed="left"), .by=sample) %>%
  mutate(coverage_bin=cut_width(coverage, width=10, boundary=, closed="left"), .by=sample) %>%
  mutate(coverage_bin=as.numeric(coverage_bin)) %>%
  reframe(coverage_bin_count=n(), .by=c(sample, coverage_bin)) %>% 
  # group_by(sample) %>%
  arrange(coverage_bin) %>%          # Sort by coverage (high to low)
  mutate(cumsum_bin = cumsum(coverage_bin_count),
         total_bin = sum(coverage_bin_count),
         cumfrac_bin = cumsum_bin / total_bin,
         frac_bin = row_number() / n(), 
         .by=sample) %>%   # Fraction of bins
  glimpse()
# df$sample %>% unique()

ks_result <- ks.test(
  filter(df, sample=='dmso') %>% pull(cumfrac_bin),
  filter(df, sample=='imr05') %>% pull(cumfrac_bin)
) %>% print()


p <- df %>% 
  # slice(seq(1, n(), by = 200)) %>% # reduce data point
  ggplot(aes(frac_bin, cumfrac_bin, color = sample)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Genome Coverage",
    y = "Cumulative Fraction of Coverage Bins",
    color = "Sample"
  ) +
  # coord_cartesian(xlim=c(0, 0.05)) +
  theme_bw(base_size = 8)
p
ggsave(str_glue("{Sys.Date()}_cdf_bins.pdf"), p, width=5, height=3)


df <- raw %>% 
  select(1:6) %>% 
  
  transmute(dmso=(imr0_1 + imr0_2 + imr0_3),
            imr05=(imr05_1 + imr05_2 + imr05_3)) %>%
  
  mutate(., mean=rowMeans(select(., 1:2))) %>%
  filter(mean < 5000) %>%
  select(!mean) %>%
  # mutate(dmso=dmso/sum(dmso),
  #        IMR05=IMR05/sum(IMR05)) %>%
  pivot_longer(everything(), names_to='sample', values_to='coverage') %>% 
  # set coverage bin width = 200
  # mutate(coverage_bin=cut_width(coverage, width=0.001, boundary=0, closed="left"), .by=sample) %>%
  mutate(coverage_bin=cut_width(coverage, width=10, boundary=, closed="left"), .by=sample) %>%
  mutate(coverage_bin=as.numeric(coverage_bin)) %>%
  reframe(coverage_bin_count=n(), .by=c(sample, coverage_bin)) %>% 
  # group_by(sample) %>%
  arrange(coverage_bin) %>%          # Sort by coverage (high to low)
  mutate(cumsum_bin = cumsum(coverage_bin_count),
         total_bin = sum(coverage_bin_count),
         cumfrac_bin = cumsum_bin / total_bin,
         frac_bin = row_number() / n(), 
         .by=sample) %>%   # Fraction of bins
  glimpse()
# df$sample %>% unique()

ks_result <- ks.test(
  filter(df, sample=='dmso') %>% pull(cumfrac_bin),
  filter(df, sample=='imr05') %>% pull(cumfrac_bin)
) %>% print()


p <- df %>% 
  # slice(seq(1, n(), by = 200)) %>% # reduce data point
  ggplot(aes(frac_bin, cumfrac_bin, color = sample)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Genome Coverage",
    y = "Cumulative Fraction of Coverage Bins",
    color = "Sample"
  ) +
  # coord_cartesian(xlim=c(0, 0.05)) +
  theme_bw(base_size = 8)
ggsave(str_glue("{Sys.Date()}_cdf_bins_5000.pdf"), p, width=5, height=3)

