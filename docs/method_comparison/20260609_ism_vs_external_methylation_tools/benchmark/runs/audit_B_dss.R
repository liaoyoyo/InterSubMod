#!/usr/bin/env Rscript
# Audit B: DSS beta-binomial (dispersion-shrunk Wald) on the SAME per-CpG counts as ISM Fisher.
# No-replicate (1 sample HP1 vs 1 sample HP1-1) -> DSS smoothing=TRUE borrows strength across CpGs
# (= the dispersion-aware behaviour ISM Fisher lacks). Output per-CpG DSS p for comparison.
suppressMessages({library(DSS); library(bsseq)})
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]; outfile <- args[2]
d <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
d <- d[d$hp1_N >= 5 & d$hp2_N >= 5, ]
M <- matrix(c(d$hp1_X, d$hp2_X), ncol = 2); Cov <- matrix(c(d$hp1_N, d$hp2_N), ncol = 2)
colnames(M) <- colnames(Cov) <- c("HP1", "HP1m1")
bs <- BSseq(chr = as.character(d$chr), pos = d$pos, M = M, Cov = Cov, sampleNames = c("HP1", "HP1m1"))
dml <- tryCatch(DMLtest(bs, group1 = "HP1", group2 = "HP1m1", smoothing = TRUE),
                error = function(e) { cat("smoothing=TRUE failed:", conditionMessage(e), "; retry FALSE\n");
                                      DMLtest(bs, group1 = "HP1", group2 = "HP1m1", smoothing = FALSE) })
write.table(dml[, c("chr", "pos", "mu1", "mu2", "diff", "stat", "pval", "fdr")],
            outfile, sep = "\t", row.names = FALSE, quote = FALSE)
cat("DSS_DONE rows:", nrow(dml), " median_pval:", median(dml$pval, na.rm = TRUE), "\n")
