# =============================================================================
# coloc_wrapper.R  --  TRUSTED TEMPLATE. Do not change the core of this file.
# -----------------------------------------------------------------------------
# This is the statistical heart of the tool and it has been checked by hand
# against the worked example in data/example/. Treat it as a known-good
# pattern: call these functions, build the pipeline around them, but do not
# rewrite the colocalisation logic. If you think something here is wrong,
# STOP and raise it with the human rather than editing it.
#
# It wraps coloc::coloc.abf() (Giambartolomei et al. 2014, the single-causal-
# variant enumeration / approximate Bayes factor method).
#
# Expected input: each trait as a data.frame with columns
#     snp   (character)  - variant id, used to align the two traits
#     beta  (numeric)    - effect size
#     se    (numeric)    - standard error of beta
#     maf   (numeric)    - minor/effect allele frequency, in (0, 1)
# =============================================================================

suppressPackageStartupMessages(library(coloc))

# Build one coloc dataset list from a standardised data.frame.
.coloc_dataset <- function(df, type, N, s = NULL, sdY = NULL) {
  d <- list(
    snp     = as.character(df$snp),
    beta    = as.numeric(df$beta),
    varbeta = as.numeric(df$se)^2,      # coloc wants varbeta = se^2
    MAF     = as.numeric(df$maf),
    N       = N,
    type    = type
  )
  if (type == "cc") {
    if (is.null(s)) stop("type='cc' needs s = proportion of samples that are cases")
    d$s <- s
  }
  if (type == "quant" && !is.null(sdY)) d$sdY <- sdY
  d
}

# Light, scientist-friendly cleaning: drop unusable rows, keep shared SNPs.
.prepare <- function(df) {
  ok <- is.finite(df$beta) & is.finite(df$se) & df$se > 0 &
        is.finite(df$maf) & df$maf > 0 & df$maf < 1
  df <- df[ok, , drop = FALSE]
  df[!duplicated(df$snp), , drop = FALSE]
}

#' Run coloc.abf on one locus for two traits.
#'
#' @return the full coloc.abf result (a list with $summary and $results),
#'         or NULL if the locus has too little overlapping data.
run_coloc_abf <- function(df1, df2,
                          type1, type2,
                          n1, n2,
                          s1 = NULL, s2 = NULL,
                          sdY1 = NULL, sdY2 = NULL,
                          p1 = 1e-4, p2 = 1e-4, p12 = 1e-5) {
  df1 <- .prepare(df1)
  df2 <- .prepare(df2)
  shared <- intersect(df1$snp, df2$snp)
  if (length(shared) < 2L) {
    warning("locus skipped: fewer than 2 shared variants")
    return(NULL)
  }
  df1 <- df1[match(shared, df1$snp), , drop = FALSE]
  df2 <- df2[match(shared, df2$snp), , drop = FALSE]

  d1 <- .coloc_dataset(df1, type1, n1, s1, sdY1)
  d2 <- .coloc_dataset(df2, type2, n2, s2, sdY2)

  # coloc.abf is chatty; keep the console clean for the demo.
  suppressWarnings(suppressMessages(
    coloc::coloc.abf(dataset1 = d1, dataset2 = d2, p1 = p1, p2 = p2, p12 = p12)
  ))
}

#' One-row tidy summary of a coloc result.
#'
#' @param res    output of run_coloc_abf()
#' @param locus  a label for the locus (e.g. "chr1:50.0-50.4Mb")
#' @return a one-row data.frame, or NULL if res is NULL.
summarise_coloc <- function(res, locus) {
  if (is.null(res)) return(NULL)
  s  <- res$summary
  rr <- res$results
  # Lead variant for the shared signal = max posterior of being THE shared SNP.
  top <- rr$snp[which.max(rr$SNP.PP.H4)]
  data.frame(
    locus      = locus,
    nsnps      = as.integer(s[["nsnps"]]),
    PP.H0      = unname(s[["PP.H0.abf"]]),
    PP.H1      = unname(s[["PP.H1.abf"]]),
    PP.H2      = unname(s[["PP.H2.abf"]]),
    PP.H3      = unname(s[["PP.H3.abf"]]),
    PP.H4      = unname(s[["PP.H4.abf"]]),
    lead_snp   = top,
    stringsAsFactors = FALSE
  )
}
