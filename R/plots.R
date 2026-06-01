# =============================================================================
# plots.R  --  the Miami plot and the per-locus locus-zoom plots.
#
# STATUS: make_miami() implemented; make_locuszoom() is the remaining stub.
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
#'
#' Trait 1 is drawn as a Manhattan pointing UP (y = -log10(p)) and trait 2
#' pointing DOWN (y = +log10(p) = the negated -log10(p)). Both share the same
#' genomic x-axis (position), faceted by chromosome. The genome-wide threshold
#' is shown as a dashed line on each half, and the loci that colocalise
#' (PP.H4 >= pp4_threshold) are boxed with a semi-transparent rectangle and
#' labelled with their locus id and PP.H4.
#'
#' @param ss1,ss2       standardised sumstats (cols incl. chr,pos,p)
#' @param loci          data.frame(locus_id, chr, start, end, PP.H4) of loci
#' @param outfile       path to write the PNG
#' @param pp4_threshold loci with PP.H4 >= this are highlighted
#' @return invisibly, the path written
make_miami <- function(ss1, ss2, loci, outfile, pp4_threshold = 0.8) {
  suppressPackageStartupMessages(library(ggplot2))

  # Guard empty inputs: nothing to draw.
  if (is.null(ss1) || is.null(ss2) || nrow(ss1) == 0 || nrow(ss2) == 0) {
    warning("make_miami: empty sumstats; skipping plot")
    return(invisible(outfile))
  }

  gw_line <- -log10(5e-8) # genome-wide significance line height

  # Build the point data: trait 1 up (+), trait 2 down (-). Coerce p away from
  # 0 first so -log10(p) stays finite.
  mk <- function(ss, trait, sign) {
    p <- pmax(as.numeric(ss$p), 1e-300)
    data.frame(
      chr = as.character(ss$chr),
      pos = as.numeric(ss$pos),
      y = sign * -log10(p),
      trait = trait,
      stringsAsFactors = FALSE
    )
  }
  pts <- rbind(
    mk(ss1, "trait 1", 1),
    mk(ss2, "trait 2", -1)
  )
  pts <- pts[is.finite(pts$pos) & is.finite(pts$y), , drop = FALSE]

  # The loci to highlight: those that colocalise above threshold.
  hit <- loci[!is.na(loci$PP.H4) & loci$PP.H4 >= pp4_threshold, , drop = FALSE]

  p <- ggplot(pts, aes(x = .data$pos, y = .data$y, colour = .data$trait))

  # Shade colocalising loci across both halves first, so points sit on top.
  if (nrow(hit) > 0) {
    p <- p + geom_rect(
      data = hit, inherit.aes = FALSE,
      aes(xmin = .data$start, xmax = .data$end, ymin = -Inf, ymax = Inf),
      fill = "gold", alpha = 0.25
    )
  }

  p <- p +
    geom_point(size = 0.6, alpha = 0.6) +
    # Threshold lines on both halves, and a separator at y = 0.
    geom_hline(
      yintercept = c(-gw_line, gw_line),
      linetype = "dashed", colour = "grey40"
    ) +
    geom_hline(yintercept = 0, colour = "black") +
    facet_grid(~chr, scales = "free_x", space = "free_x") +
    scale_colour_manual(
      values = c("trait 1" = "#1f77b4", "trait 2" = "#d62728")
    ) +
    labs(
      x = "Position (bp)",
      y = expression(-log[10](p) ~ "  (trait 1 up / trait 2 down)"),
      colour = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")

  # Label each highlighted locus near the top of the upper half.
  if (nrow(hit) > 0) {
    lab <- data.frame(
      chr = as.character(hit$chr),
      pos = (as.numeric(hit$start) + as.numeric(hit$end)) / 2,
      y = max(pts$y, na.rm = TRUE),
      label = sprintf("%s\nPP.H4=%.2f", hit$locus_id, hit$PP.H4),
      stringsAsFactors = FALSE
    )
    p <- p + geom_text(
      data = lab, inherit.aes = FALSE,
      aes(x = .data$pos, y = .data$y, label = .data$label),
      size = 3, vjust = 1, colour = "grey20"
    )
  }

  # Render headlessly to PNG.
  png(outfile, width = 1400, height = 800, res = 120)
  on.exit(dev.off(), add = TRUE)
  print(p)
  invisible(outfile)
}

#' Locus-zoom plot for a single locus, for one or both traits.
#' @param ss        standardised sumstats for the trait
#' @param chr,start,end  region
#' @param outfile   path to write the PNG
#' @param ens_db    an EnsDb object for the gene track, or NULL to skip it
make_locuszoom <- function(ss, chr, start, end, outfile, ens_db = NULL) {
  stop("TODO: regional association plot via locuszoomr::locus() + locus_plot()")
}
