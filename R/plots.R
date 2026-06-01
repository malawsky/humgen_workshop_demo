# =============================================================================
# plots.R  --  the Miami plot and the per-locus locus-zoom plots.
#
# STATUS: STUB. Implement the TODOs during the build.
# =============================================================================
#
# Miami plot: two Manhattan plots back to back -- trait 1 pointing up, trait 2
# pointing down, sharing the x-axis (genomic position). Highlight the loci that
# colocalise (PP.H4 above threshold).
#   * The `topr` package (CRAN) draws Manhattan/Miami plots directly from a
#     data.frame with CHROM, POS, P columns -- a good starting point.
#   * A ggplot2 version is fine too if you prefer full control of the highlight.
#
# Locus zoom: a regional association plot for one locus with a gene track.
#   * The `locuszoomr` package (CRAN): build a locus object with locus(),
#     then plot with locus_plot(). Gene tracks need an Ensembl annotation
#     package, e.g. EnsDb.Hsapiens.v86 (GRCh38) from Bioconductor.
#   * LD overlay is optional (LDlink API) and NOT needed for the demo.
# -----------------------------------------------------------------------------

#' Miami plot of the two traits with colocalising loci highlighted.
#' @param ss1,ss2       standardised sumstats (cols incl. chr,pos,p)
#' @param coloc_table   summarised coloc results (one row per locus)
#' @param outfile       path to write the PNG
#' @param pp4_threshold loci with PP.H4 >= this are highlighted
make_miami <- function(ss1, ss2, coloc_table, outfile, pp4_threshold = 0.8) {
  stop("TODO: draw a Miami plot (trait1 up / trait2 down) and highlight coloc hits")
}

#' Locus-zoom plot for a single locus, for one or both traits.
#' @param ss        standardised sumstats for the trait
#' @param chr,start,end  region
#' @param outfile   path to write the PNG
#' @param ens_db    an EnsDb object for the gene track, or NULL to skip it
make_locuszoom <- function(ss, chr, start, end, outfile, ens_db = NULL) {
  stop("TODO: regional association plot via locuszoomr::locus() + locus_plot()")
}
