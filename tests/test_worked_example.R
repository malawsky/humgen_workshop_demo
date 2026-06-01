#!/usr/bin/env Rscript
# =============================================================================
# test_worked_example.R  --  check the trusted coloc wrapper against a planted
#                            answer. This passes BEFORE any pipeline is built,
#                            because it calls coloc_wrapper.R directly on the
#                            synthetic data in data/example/.
#
# Run from the repo root:   Rscript tests/test_worked_example.R
#
# The synthetic data (see data/example/README.md) contains:
#   Locus A  -> both traits share one causal variant  -> expect PP.H4 high
#   Locus B  -> traits have distinct causal variants  -> expect H3 to dominate
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})
source(file.path("R", "coloc_wrapper.R"))

read_std <- function(path) {
  d <- fread(path)
  data.frame(
    snp  = d$hm_rsid,
    pos  = d$hm_pos,
    beta = d$hm_beta,
    se   = d$standard_error,
    maf  = d$hm_effect_allele_frequency,
    p    = d$p_value,
    stringsAsFactors = FALSE
  )
}

ss1 <- read_std("data/example/trait1.h.tsv.gz")
ss2 <- read_std("data/example/trait2.h.tsv.gz")
truth <- fread("tests/expected_truth.tsv")

region <- function(ss, lo, hi) ss[ss$pos >= lo & ss$pos <= hi, , drop = FALSE]

fail <- 0L
check <- function(cond, msg) {
  cat(ifelse(cond, "PASS  ", "FAIL  "), msg, "\n", sep = "")
  if (!cond) fail <<- fail + 1L
}

# --- Locus A: expect a shared signal (PP.H4 high) ---------------------------
A <- truth[truth$locus == "A", ]
resA <- run_coloc_abf(region(ss1, A$region_start, A$region_end),
                      region(ss2, A$region_start, A$region_end),
                      type1 = "quant", type2 = "quant",
                      n1 = 50000, n2 = 50000)
sumA <- summarise_coloc(resA, "A")
cat(sprintf("Locus A: PP.H3=%.3f  PP.H4=%.3f\n", sumA$PP.H3, sumA$PP.H4))
check(sumA$PP.H4 > 0.8, "Locus A colocalises (PP.H4 > 0.8)")

# --- Locus B: expect distinct signals (H3 dominates) ------------------------
B <- truth[truth$locus == "B", ]
resB <- run_coloc_abf(region(ss1, B$region_start, B$region_end),
                      region(ss2, B$region_start, B$region_end),
                      type1 = "quant", type2 = "quant",
                      n1 = 50000, n2 = 50000)
sumB <- summarise_coloc(resB, "B")
cat(sprintf("Locus B: PP.H3=%.3f  PP.H4=%.3f\n", sumB$PP.H3, sumB$PP.H4))
check(sumB$PP.H3 > sumB$PP.H4, "Locus B does NOT colocalise (PP.H3 > PP.H4)")

cat(if (fail == 0L) "\nAll worked-example checks passed.\n" else sprintf("\n%d check(s) failed.\n", fail))
quit(status = if (fail == 0L) 0L else 1L)
