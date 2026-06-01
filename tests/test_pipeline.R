#!/usr/bin/env Rscript
# =============================================================================
# test_pipeline.R  --  known-answer INTEGRATION test for the coloc-pair CLI.
#
# Runs R/coloc_pair.R end-to-end on the synthetic worked example
# (data/example/) into a TEMPORARY outdir, then checks the written
# coloc_results.tsv against the planted answer:
#   Locus A (~chr1:50Mb) -> shared causal variant  -> PP.H4 > 0.8 (top row)
#   Locus B (~chr1:60Mb) -> distinct causals        -> PP.H3 > PP.H4
#
# Run from the repo root:   Rscript tests/test_pipeline.R
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

fail <- 0L
check <- function(cond, msg) {
  cat(ifelse(isTRUE(cond), "PASS  ", "FAIL  "), msg, "\n", sep = "")
  if (!isTRUE(cond)) fail <<- fail + 1L
}

# --- 1. Run the actual CLI into a temporary output dir ----------------------
outdir <- tempfile("coloc_pipeline_test_")
status <- system2(
  "Rscript",
  c(
    "R/coloc_pair.R",
    "--study1", "data/example/trait1.h.tsv.gz", "--type1", "quant",
    "--study2", "data/example/trait2.h.tsv.gz", "--type2", "quant",
    "--n1", "50000", "--n2", "50000",
    "--outdir", outdir
  )
)
cat(sprintf("CLI exit status: %d\n", status))
check(status == 0L, "CLI exited 0")

# --- 2. Results table was written -------------------------------------------
res_path <- file.path(outdir, "coloc_results.tsv")
check(file.exists(res_path), "coloc_results.tsv written")

if (!file.exists(res_path)) {
  cat(sprintf("\n%d check(s) failed.\n", fail))
  quit(status = 1L)
}

# --- 3. Known-answer checks on the table ------------------------------------
res <- as.data.frame(fread(res_path))
cat("Rows in results table: ", nrow(res), "\n", sep = "")
check(nrow(res) == 2L, "exactly 2 locus rows")

pp_cols <- c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")
have_cols <- all(c("locus", pp_cols) %in% names(res))
check(have_cols, "results table has locus + PP.H0..H4 columns")

if (nrow(res) >= 1L && have_cols) {
  # Table is written sorted by PP.H4 desc; enforce here for robustness.
  res <- res[order(-res$PP.H4), , drop = FALSE]

  # Sanity: each PP in [0,1], and the five PPs sum to ~1 per row.
  pp_in_range <- all(vapply(pp_cols, function(c) all(res[[c]] >= 0 & res[[c]] <= 1), logical(1)))
  check(pp_in_range, "all PP.H0..H4 within [0,1]")
  row_sums <- rowSums(res[, pp_cols])
  cat("Per-row PP sums: ", paste(sprintf("%.6f", row_sums), collapse = ", "), "\n", sep = "")
  check(all(abs(row_sums - 1) < 1e-6), "each row's PP.H0..H4 sums to ~1 (within 1e-6)")

  # No duplicate locus ids.
  check(!any(duplicated(res$locus)), "no duplicate locus ids")

  # Identify loci by position: A ~50Mb, B ~60Mb (label contains the Mb value).
  topfn <- res[1, ]
  cat(sprintf("Top row (max PP.H4): locus=%s  PP.H3=%.3f  PP.H4=%.3f\n",
              topfn$locus, topfn$PP.H3, topfn$PP.H4))
  is_A <- grepl("50", topfn$locus)
  check(is_A && topfn$PP.H4 > 0.8,
        sprintf("top locus is the ~chr1:50Mb locus and colocalises (PP.H4=%.3f > 0.8)", topfn$PP.H4))

  if (nrow(res) >= 2L) {
    other <- res[2, ]
    cat(sprintf("Other row: locus=%s  PP.H3=%.3f  PP.H4=%.3f\n",
                other$locus, other$PP.H3, other$PP.H4))
    is_B <- grepl("60", other$locus)
    check(is_B && other$PP.H3 > other$PP.H4,
          sprintf("other locus is the ~chr1:60Mb locus and does NOT colocalise (PP.H3=%.3f > PP.H4=%.3f)",
                  other$PP.H3, other$PP.H4))
  }
}

cat(if (fail == 0L) "\nAll pipeline checks passed.\n" else sprintf("\n%d check(s) failed.\n", fail))
quit(status = if (fail == 0L) 0L else 1L)
