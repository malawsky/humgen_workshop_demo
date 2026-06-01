# =============================================================================
# plots.R  --  the Miami plot and the per-locus locus-zoom plots.
#
# STATUS: implemented. Both plots are self-contained ggplot2 (no Bioconductor).
# =============================================================================
#
# Miami plot: two Manhattan plots back to back -- trait 1 pointing up, trait 2
# pointing down, sharing the x-axis (genomic position), with colocalising loci
# (PP.H4 above threshold) shaded and labelled.
#
# Locus zoom: a regional -log10(p) plot for one locus, overlaying both traits so
# coinciding peaks are visible, with the lead SNP marked. Gene tracks are
# deliberately omitted (they would require the EnsDb / Bioconductor annotation
# stack); the demo data is synthetic so a gene track adds nothing here.
# -----------------------------------------------------------------------------

#' Thin a Manhattan/Miami point cloud so it renders in seconds, not minutes.
#'
#' Genome-wide sumstats hold tens of millions of variants; almost all have a
#' tiny -log10(p) and pile into a dense band at the axis that is visually
#' saturated long before every point is drawn, yet still costs the renderer.
#' Keep every point with abs(y) >= keep_y (all peaks and near-peaks, where y is
#' the signed -log10(p)) and randomly thin the sub-threshold band to at most
#' max_low points. Deterministic via `seed`; restores the global RNG state.
#'
#' @param pts     data.frame with a numeric `y` column (signed -log10 p)
#' @param keep_y  points with abs(y) >= keep_y are always kept
#' @param max_low cap on the number of sub-threshold (abs(y) < keep_y) points
#' @param seed    RNG seed for reproducible thinning
#' @return pts with the low band downsampled
downsample_miami_points <- function(pts, keep_y = 4, max_low = 5e5L, seed = 1L) {
  if (is.null(pts) || nrow(pts) == 0L) {
    return(pts)
  }
  low <- which(abs(pts$y) < keep_y)
  if (length(low) <= max_low) {
    return(pts)
  }
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)
  drop <- sample(low, length(low) - max_low)
  pts[-drop, , drop = FALSE]
}

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

  # Thin the dense sub-threshold cloud so a genome-wide plot renders in seconds.
  pts <- downsample_miami_points(pts)

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

#' Regional (locus-zoom) plot: the two traits stacked for visual comparison.
#'
#' Subsets both traits to the locus, then ZOOMS to a window of +/- `flank` bp
#' around the lead SNP (clamped to the locus) so the signal is readable instead
#' of lost in a wide merged window. The two traits are drawn in SEPARATE
#' stacked panels (trait 1 on top, trait 2 below) sharing the x-axis, so
#' coinciding peaks line up vertically. The lead SNP is a vertical dashed line
#' in both panels and the colocalisation posterior is annotated as a subtitle.
#' No gene track (deliberately self-contained ggplot2; no Bioconductor stack).
#'
#' @param ss1,ss2   standardised sumstats (cols incl. snp,chr,pos,p)
#' @param chr,start,end  locus region
#' @param outfile   path to write the PNG
#' @param lead_snp  optional lead SNP id; centres the window and is marked
#' @param locus_id  optional locus label used as the plot title
#' @param pp4,pp3   optional coloc posteriors, annotated at the top
#' @param flank     half-width in bp of the focused window around the lead SNP
#' @return invisibly, the path written
make_locuszoom <- function(ss1, ss2, chr, start, end, outfile,
                           lead_snp = NULL, locus_id = NULL,
                           pp4 = NA, pp3 = NA, flank = 2.5e5) {
  suppressPackageStartupMessages(library(ggplot2))

  gw_line <- -log10(5e-8) # genome-wide significance line height

  # Subset each trait to the region (same logic as extract_region).
  region <- function(ss) {
    if (is.null(ss) || nrow(ss) == 0) {
      return(ss[0, , drop = FALSE])
    }
    ss[ss$chr == chr & ss$pos >= start & ss$pos <= end, , drop = FALSE]
  }
  r1 <- region(ss1)
  r2 <- region(ss2)

  # Build the point data, coercing p away from 0 so -log10(p) stays finite.
  mk <- function(ss, trait) {
    if (nrow(ss) == 0) {
      return(NULL)
    }
    p <- pmax(as.numeric(ss$p), 1e-300)
    data.frame(
      snp = as.character(ss$snp),
      pos = as.numeric(ss$pos),
      y = -log10(p),
      trait = trait,
      stringsAsFactors = FALSE
    )
  }
  pts <- rbind(mk(r1, "trait 1"), mk(r2, "trait 2"))

  # Guard empty region: nothing to draw.
  if (is.null(pts) || nrow(pts) == 0) {
    warning("make_locuszoom: no variants in region; skipping plot")
    return(invisible(outfile))
  }
  pts <- pts[is.finite(pts$pos) & is.finite(pts$y), , drop = FALSE]
  if (nrow(pts) == 0) {
    warning("make_locuszoom: no finite points in region; skipping plot")
    return(invisible(outfile))
  }

  # Centre on the lead SNP (fall back to the strongest point), then zoom to a
  # focused window so the locus is readable rather than a wide merged span.
  lead_pos <- NA_real_
  if (!is.null(lead_snp) && !is.na(lead_snp)) {
    hit <- pts$pos[pts$snp == lead_snp]
    if (length(hit) > 0) lead_pos <- hit[1]
  }
  if (is.na(lead_pos)) lead_pos <- pts$pos[which.max(pts$y)]

  win_lo <- max(start, lead_pos - flank)
  win_hi <- min(end, lead_pos + flank)
  pts <- pts[pts$pos >= win_lo & pts$pos <= win_hi, , drop = FALSE]
  if (nrow(pts) == 0) {
    warning("make_locuszoom: no points in focused window; skipping plot")
    return(invisible(outfile))
  }

  # Stack the two traits (trait 1 on top) so coinciding peaks line up.
  pts$trait <- factor(pts$trait, levels = c("trait 1", "trait 2"))

  title <- if (!is.null(locus_id)) {
    locus_id
  } else {
    sprintf("chr%s:%g-%g", chr, start, end)
  }

  # Annotate the colocalisation posterior (and lead SNP) at the top.
  sub_bits <- character(0)
  if (!is.na(pp4)) sub_bits <- c(sub_bits, sprintf("PP.H4 = %.3f", pp4))
  if (!is.na(pp3)) sub_bits <- c(sub_bits, sprintf("PP.H3 = %.3f", pp3))
  if (!is.null(lead_snp) && !is.na(lead_snp)) {
    sub_bits <- c(sub_bits, sprintf("lead %s", lead_snp))
  }
  subtitle <- if (length(sub_bits)) paste(sub_bits, collapse = "   |   ") else NULL

  p <- ggplot(pts, aes(x = .data$pos, y = .data$y, colour = .data$trait)) +
    geom_point(size = 1.1, alpha = 0.7) +
    geom_hline(yintercept = gw_line, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = lead_pos, linetype = "dashed", colour = "grey20") +
    facet_grid(trait ~ ., scales = "free_y") +
    scale_colour_manual(
      values = c("trait 1" = "#1f77b4", "trait 2" = "#d62728"),
      guide = "none"
    ) +
    labs(
      x = "Position (bp)", y = expression(-log[10](p)),
      title = title, subtitle = subtitle
    ) +
    theme_bw(base_size = 12) +
    theme(strip.text = element_text(face = "bold"))

  # Render headlessly to PNG.
  png(outfile, width = 1100, height = 800, res = 120)
  on.exit(dev.off(), add = TRUE)
  print(p)
  invisible(outfile)
}
