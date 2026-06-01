#!/usr/bin/env Rscript
# =============================================================================
# test_extract_region.R  --  the keyed (data.table) extract_region must return
# EXACTLY the same rows as the plain data.frame scan, for every locus.
#
# This guards the (chr, pos)-keyed binary-search lookup that replaced the old
# full-genome scan. The bug it pins: an unquoted `chr` inside ss[.(chr)] would
# resolve to the table's chr COLUMN and pull in every chromosome. We assert the
# keyed result equals the boolean scan on a MULTI-chromosome fixture.
#
# Run from the repo root:   Rscript tests/test_extract_region.R
# =============================================================================

suppressPackageStartupMessages(library(data.table))
source("R/define_loci.R")

fail <- 0L
check <- function(cond, msg) {
  cat(ifelse(isTRUE(cond), "PASS  ", "FAIL  "), msg, "\n", sep = "")
  if (!isTRUE(cond)) fail <<- fail + 1L
}

set.seed(1)
n <- 5000
df <- data.frame(
  snp = paste0("rs", seq_len(n)),
  chr = sample(c(as.character(1:22), "X"), n, replace = TRUE),
  pos = sample(1:1000000L, n, replace = TRUE),
  ea = "A", oa = "G", maf = 0.2, beta = 0.1, se = 0.05, p = 0.01,
  stringsAsFactors = FALSE
)

# Reference: the plain data.frame scan (the trusted, pre-optimisation logic).
scan_region <- function(ss, chr, start, end) {
  ss[ss$chr == chr & ss$pos >= start & ss$pos <= end, , drop = FALSE]
}

# Keyed input, as the pipeline builds it.
dt <- data.table::as.data.table(df)
data.table::setkey(dt, chr, pos)

norm <- function(x) {
  x <- as.data.frame(x)
  x <- x[order(x$snp), , drop = FALSE]
  rownames(x) <- NULL
  x
}

queries <- list(
  c("1", 1, 1000000), c("X", 200000, 800000),
  c("7", 400000, 450000), c("22", 1, 50000),
  c("3", 999999, 1000000), c("9", 1, 1) # likely-empty windows too
)
all_ok <- TRUE
for (q in queries) {
  ch <- q[1]
  st <- as.numeric(q[2])
  en <- as.numeric(q[3])
  got <- norm(extract_region(dt, ch, st, en))
  exp <- norm(scan_region(df, ch, st, en))
  ok <- identical(got, exp)
  if (!ok) all_ok <- FALSE
  # Every returned row must be on the requested chromosome (the bug's signature).
  if (nrow(got) > 0 && !all(got$chr == ch)) all_ok <- FALSE
}
check(all_ok, "keyed extract_region == boolean scan across chromosomes")

# Explicit duplicate / wrong-chr guard on a known window.
g <- extract_region(dt, "1", 1, 1000000)
check(all(g$chr == "1"), "no wrong-chromosome rows leak in")
check(!any(duplicated(g$snp)), "no duplicated rows")

if (fail > 0L) {
  cat(sprintf("\n%d check(s) FAILED.\n", fail))
  quit(status = 1L)
}
cat("\nAll extract_region checks passed.\n")
