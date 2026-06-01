#!/usr/bin/env Rscript
# =============================================================================
# make_example.R  --  R version of the synthetic worked-example generator.
# Produces harmonised-format sumstats with a PLANTED answer (see README.md).
# The committed *.h.tsv.gz are the source of truth; regenerate only if needed,
# and if you do, re-check tests/expected_truth.tsv.
# Run from data/example/:   Rscript make_example.R
# =============================================================================
set.seed(20260601)
SE <- 0.02; alleles <- c("A", "C", "G", "T")

region <- function(start, n, spacing = 1000) start + (0:(n - 1)) * spacing
bump <- function(pos, causal, zpeak, width = 40000)
  zpeak * exp(-((pos - causal)^2) / (2 * width^2)) + rnorm(length(pos), 0, 0.25)

posA <- region(50e6, 400); posB <- region(60e6, 400)
pos <- c(posA, posB); n <- length(pos)
ea <- sample(alleles, n, TRUE); oa <- sample(alleles, n, TRUE)
same <- ea == oa; oa[same] <- alleles[(match(ea[same], alleles) %% 4) + 1]
eaf <- runif(n, 0.10, 0.45)
rsid <- paste0("rs", 9000000 + 0:(n - 1)); chrom <- rep(1L, n)

make_trait <- function(causalA, causalB) {
  z <- numeric(n); inA <- pos < 55e6
  z[inA]  <- bump(pos[inA],  causalA, 9.5)
  z[!inA] <- bump(pos[!inA], causalB, 9.0)
  list(beta = z * SE, se = rep(SE, n), p = pmax(2 * pnorm(-abs(z)), 1e-300))
}

write_h <- function(path, tr) {
  df <- data.frame(hm_rsid = rsid, hm_chrom = chrom, hm_pos = pos,
                   hm_other_allele = oa, hm_effect_allele = ea,
                   hm_effect_allele_frequency = round(eaf, 4),
                   hm_beta = round(tr$beta, 6), standard_error = round(tr$se, 6),
                   p_value = formatC(tr$p, format = "e", digits = 3))
  gz <- gzfile(path, "w"); write.table(df, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
}

sharedA <- posA[201]; b1 <- posB[151]; b2 <- posB[351]
write_h("trait1.h.tsv.gz", make_trait(sharedA, b1))
write_h("trait2.h.tsv.gz", make_trait(sharedA, b2))
cat("wrote trait1.h.tsv.gz, trait2.h.tsv.gz\n")
cat(sprintf("Locus A shared causal chr1:%d -> PP.H4 high; Locus B distinct %d/%d -> PP.H3 high\n",
            sharedA, b1, b2))
