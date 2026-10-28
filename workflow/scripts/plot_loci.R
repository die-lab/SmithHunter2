#!/usr/bin/env Rscript
# Profile of each passing locus: depth (dashed line) and reads starting (5') or ending (3')
# at each position (bars), six loci per page. Adapted from plot_clusters.R of the
# smithhunter24 port of v0 (itself after Makeplots.R).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(ggforce)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript plot_loci.R <locus_reads.tsv.gz> <loci.tsv> <loci.pdf>")
}
reads <- read.delim(gzfile(args[1]), stringsAsFactors = FALSE)
loci <- read.delim(args[2], stringsAsFactors = FALSE, check.names = FALSE,
                   na.strings = character(0))
out_pdf <- args[3]

RED <- "#DC3220"
BLUE <- "#005AB5"
PER_PAGE <- 6

pass <- loci %>%
  filter(pass == 1) %>%
  mutate(facet = paste0(locus_id, " (", strand, ")  ", comma(round(total_count)), " reads",
                        "  5' ", five_score, ", 3' ", three_score, "\n", rep_sequence))

pdf(out_pdf, onefile = TRUE, paper = "a4", width = 21 / 2.54, height = 29.7 / 2.54)
if (nrow(pass) == 0) {
  plot.new()
  text(0.5, 0.5, "No passing loci")
} else {
  r <- inner_join(reads, select(pass, locus_id, strand, facet), by = "locus_id")
  ends <- bind_rows(
    r %>% mutate(pos = ifelse(strand == "+", start, end), end_type = "5' ends"),
    r %>% mutate(pos = ifelse(strand == "+", end, start), end_type = "3' ends")) %>%
    group_by(facet, end_type, pos) %>%
    summarise(n = sum(count), .groups = "drop")
  span <- r$end - r$start + 1
  idx <- rep(seq_len(nrow(r)), span)
  depth <- data.frame(facet = r$facet[idx],
                      pos = r$start[idx] + sequence(span) - 1,
                      n = r$count[idx]) %>%
    group_by(facet, pos) %>%
    summarise(n = sum(n), .groups = "drop")
  levels <- pass$facet[order(pass$locus_id)]
  ends$facet <- factor(ends$facet, levels = levels)
  depth$facet <- factor(depth$facet, levels = levels)

  base <- ggplot() +
    geom_col(data = ends, aes(pos, n, fill = end_type, colour = end_type), width = 0.9,
             alpha = 0.5, position = "identity") +
    geom_line(data = depth, aes(pos, n), linetype = "dashed", colour = "grey20") +
    scale_fill_manual(name = NULL, values = c("5' ends" = RED, "3' ends" = BLUE)) +
    scale_colour_manual(name = NULL, values = c("5' ends" = RED, "3' ends" = BLUE)) +
    expand_limits(y = 0) +
    labs(x = "position", y = "reads") +
    theme_bw() +
    theme(legend.position = "bottom", strip.text = element_text(size = 6, family = "mono"))
  pages <- ceiling(length(levels) / PER_PAGE)
  for (pg in seq_len(pages)) {
    print(base + facet_wrap_paginate(~facet, ncol = 2, nrow = 3, scales = "free", page = pg))
  }
}
invisible(dev.off())
