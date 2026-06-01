#!/usr/bin/env Rscript
# =============================================================================
# test_fst_cache.R  --  offline round-trip test for the parsed-frame .fst cache.
#
# load_sumstats() caches the parsed, standardised frame as <study>.std.v<N>.fst
# so re-runs skip the multi-minute parse of full-genome files. That shortcut is
# only safe if fst::read_fst returns a frame IDENTICAL to a fresh parse. This
# test pins that: parse the synthetic example, round-trip it through fst, and
# assert the data and per-column types match exactly.
#
# Run from the repo root:   Rscript tests/test_fst_cache.R
# =============================================================================

suppressPackageStartupMessages({
  library(fst)
})
source("R/download_sumstats.R")

fail <- 0L
check <- function(cond, msg) {
  cat(ifelse(isTRUE(cond), "PASS  ", "FAIL  "), msg, "\n", sep = "")
  if (!isTRUE(cond)) fail <<- fail + 1L
}

fresh <- read_harmonised("data/example/trait1.h.tsv.gz")

tmp <- tempfile(fileext = ".fst")
fst::write_fst(fresh, tmp)
cached <- fst::read_fst(tmp)

check(identical(fresh, cached), "fst round-trip is identical() to a fresh parse")
check(
  identical(sapply(fresh, class), sapply(cached, class)),
  "per-column classes survive the round-trip"
)
check(identical(names(cached), STD_COLS), "cached frame has STD_COLS in order")
check(is.character(cached$chr), "chr stays character (X/Y survive)")
check(is.numeric(cached$p), "p stays numeric")

if (fail > 0L) {
  cat(sprintf("\n%d check(s) FAILED.\n", fail))
  quit(status = 1L)
}
cat("\nAll fst-cache checks passed.\n")
