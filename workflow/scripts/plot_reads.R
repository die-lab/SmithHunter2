#!/usr/bin/env Rscript
# Reads through the pipeline: funnel (all samples and per sample), genome of origin,
# read length by origin. Adapted from plot_funnel.R of the smithhunter24 port of v0.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript plot_reads.R <read_funnel.tsv> <origin_summary.tsv> ",
       "<length_distribution.tsv> <reads.pdf>")
}
read_table <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = character(0))
}
funnel <- read_table(args[1])
origin <- read_table(args[2])
lengths <- read_table(args[3])
out_pdf <- args[4]

samples <- setdiff(colnames(funnel), c("stage", "label", "total"))
group_label <- function(class, genomes) ifelse(genomes == "", class, paste0(class, " (", genomes, ")"))

# 1. funnel
tot <- funnel %>%
  mutate(idx = row_number(), pct = 100 * total / total[1])
p_funnel <- ggplot(tot) +
  geom_rect(aes(xmin = idx - 0.4, xmax = idx + 0.4, ymin = -total / 2, ymax = total / 2),
            fill = "#9DC3E6") +
  geom_text(aes(x = idx, y = 0, label = paste0(comma(round(total)), "  (", round(pct, 1), "%)")),
            size = 3.2, fontface = "bold") +
  scale_x_continuous(breaks = tot$idx, labels = tot$label, trans = "reverse") +
  coord_flip() +
  labs(title = "Read funnel, all samples", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"))

per <- funnel %>%
  mutate(idx = row_number()) %>%
  select(idx, label, all_of(samples)) %>%
  pivot_longer(all_of(samples), names_to = "sample", values_to = "reads") %>%
  group_by(sample) %>%
  mutate(pct = 100 * reads / reads[idx == 1]) %>%
  ungroup()
p_pct <- ggplot(per, aes(idx, pct, colour = sample)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = tot$idx, labels = tot$label) +
  labs(x = NULL, y = "% of raw reads", colour = NULL, title = "Per sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")

# 2. origin
orig <- origin %>%
  mutate(group = group_label(class, genomes)) %>%
  pivot_longer(all_of(samples), names_to = "sample", values_to = "reads")
p_origin <- ggplot(orig, aes(sample, reads, fill = group)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = percent) +
  labs(x = NULL, y = "share of reads in the length range", fill = NULL,
       title = "Genome of origin") +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 2))

# 3. lengths
len <- lengths %>%
  mutate(group = group_label(class, genomes)) %>%
  pivot_longer(all_of(samples), names_to = "sample", values_to = "reads") %>%
  group_by(sample, group, length) %>%
  summarise(reads = sum(reads), .groups = "drop")
p_len <- ggplot(len, aes(length, reads, colour = group)) +
  geom_line() +
  geom_point(size = 0.8) +
  facet_wrap(~sample, ncol = 1, scales = "free_y") +
  labs(x = "read length (nt)", y = "reads", colour = NULL, title = "Read length by origin") +
  theme_bw() +
  theme(legend.position = "bottom")

pdf(out_pdf, onefile = TRUE, width = 21 / 2.54, height = 29.7 / 2.54)
print(p_funnel / p_pct)
print(p_origin)
print(p_len)
invisible(dev.off())
