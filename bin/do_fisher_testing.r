#!/usr/bin/env Rscript
library(argparse)
library(tidyr)
library(dplyr)

parser <- argparse::ArgumentParser(description="")
parser$add_argument('-i1', type="character")
parser$add_argument('-i2', type="character")

args <- parser$parse_args()


kmers1 <- read.table(args$i1, header = TRUE, sep = ';')
kmers2 <- read.table(args$i2, header = TRUE, sep = ';')
kmers <- kmers1 %>%
	left_join(kmers2, by="kmer") %>%
	mutate(sum.x = sum(count.x),
	       sum.y = sum(count.y))

library(stats)
fisher <- function(a, b, c, d) {
	data <- matrix(c(a,b,c,d), ncol=2)

	fisher.test(data, alternative="less")
}

options(digits=22)
fisher_tab <- kmers %>%
	rowwise() %>%
	mutate(p = fisher(count.x, sum.x, count.y, sum.y)$p.value) %>%
	mutate(Z = qnorm(p), Z_ = Z*(-1))
write.table(fisher_tab, "fisher.csv", sep=";", row.names=FALSE)
	
