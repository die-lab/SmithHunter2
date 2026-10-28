#!/usr/bin/env Rscript
# Coverage of each focus genome contig (+ strand up, - strand down, pseudo-log scale) with
# the passing loci labelled by number. Adapted from the coverage overview and candidate map
# of plot_coverage.R in the smithhunter24 port of v0 (four genome segments, ggrepel labels).

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript plot_genomes.R <coverage.tsv.gz> <loci.tsv> <genomes.pdf>")
}
cov <- read.delim(gzfile(args[1]), stringsAsFactors = FALSE, check.names = FALSE)
loci <- read.delim(args[2], stringsAsFactors = FALSE, check.names = FALSE,
                   na.strings = character(0))
out_pdf <- args[3]

RED <- "#DC3220"
BLUE <- "#005AB5"
MAX_CONTIGS <- 30
BINS <- 3000

pass <- loci %>%
  filter(pass == 1) %>%
  mutate(sign = ifelse(strand == "+", 1, -1),
         label = sub("^.*_0*", "", locus_id))

contigs <- cov %>%
  group_by(genome, contig) %>%
  summarise(reads = sum(depth), .groups = "drop") %>%
  left_join(count(pass, genome, contig, name = "n_pass"), by = c("genome", "contig")) %>%
  mutate(n_pass = replace_na(n_pass, 0)) %>%
  arrange(desc(n_pass), desc(reads))
if (nrow(contigs) > MAX_CONTIGS) {
  contigs <- head(filter(contigs, n_pass > 0 | row_number() <= MAX_CONTIGS), MAX_CONTIGS)
}

signed_breaks <- function(m) {
  p <- 10^(1:max(1, ceiling(log10(m))))
  c(-rev(p), 0, p)
}

plot_contig <- function(g, ctg) {
  d <- filter(cov, genome == g, contig == ctg)
  l <- filter(pass, genome == g, contig == ctg)
  len <- max(d$position, l$start, l$end)
  w <- max(1, ceiling(len / BINS))
  b <- d %>%
    mutate(bin = (position - 1) %/% w) %>%
    group_by(strand, bin) %>%
    summarise(depth = sum(depth) / w, .groups = "drop") %>%
    complete(strand = c("+", "-"), bin = 0:((len - 1) %/% w), fill = list(depth = 0)) %>%
    mutate(x = bin * w + 1, y = ifelse(strand == "+", depth, -depth))
  l <- mutate(l, x = start, y = sign * total_count)
  m <- max(abs(c(b$y, l$y)), 1)
  n_seg <- if (len > 4000) 4 else 1
  edges <- round(seq(0, len, length.out = n_seg + 1))
  segment <- function(i) {
    lo <- edges[i]
    hi <- edges[i + 1]
    bs <- filter(b, x > lo - w, x <= hi)
    ls <- filter(l, x > lo, x <= hi)
    p <- ggplot() +
      geom_col(data = bs, aes(x, y, fill = strand), width = w, just = 0) +
      geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.2) +
      scale_fill_manual(values = c("+" = RED, "-" = BLUE)) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 1, base = 10),
                         breaks = signed_breaks(m), labels = function(v) comma(abs(v))) +
      coord_cartesian(xlim = c(lo, hi), ylim = c(-m, m) * 1.15) +
      labs(x = paste("position on", ctg), y = "depth", fill = "strand") +
      theme_bw() +
      theme(legend.position = if (i == 1) "top" else "none")
    if (nrow(ls) > 0) {
      p <- p +
        geom_point(data = ls, aes(x, y), shape = 21, fill = "white", size = 1.6) +
        geom_text_repel(data = ls, aes(x, y, label = label), size = 2.5,
                        min.segment.length = 0, max.overlaps = Inf, seed = 1)
    }
    p
  }
  wrap_plots(lapply(seq_len(n_seg), segment), ncol = 1) +
    plot_annotation(
      title = paste0(g, " / ", ctg, ": coverage and passing loci (", nrow(l), ")"),
      subtitle = "points: passing loci at their start, height = reads; labels = locus number")
}

pdf(out_pdf, onefile = TRUE, width = 21 / 2.54, height = 29.7 / 2.54)
if (nrow(contigs) == 0) {
  plot.new()
  text(0.5, 0.5, "No reads on focus genomes")
} else {
  for (k in seq_len(nrow(contigs))) print(plot_contig(contigs$genome[k], contigs$contig[k]))
}
invisible(dev.off())
