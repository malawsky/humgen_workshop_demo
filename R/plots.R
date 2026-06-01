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

#' Regional association (locus-zoom) plot overlaying both traits.
#'
#' Subsets each trait's standardised sumstats to the region [start,end] on
#' `chr` and draws a single panel of -log10(p) vs position, with the two traits
#' overlaid in distinct colours so coinciding peaks are visible. The
#' genome-wide threshold is a dashed line; if `lead_snp` is supplied and found
#' in either trait, its position is marked with a vertical dashed line. No gene
#' track (deliberately self-contained ggplot2; no Bioconductor stack).
#'
#' @param ss1,ss2   standardised sumstats (cols incl. snp,chr,pos,p)
#' @param chr,start,end  region
#' @param outfile   path to write the PNG
#' @param lead_snp  optional lead SNP id to mark with a vertical line
#' @param locus_id  optional locus label used as the plot title
#' @return invisibly, the path written
make_locuszoom <- function(ss1, ss2, chr, start, end, outfile,
                           lead_snp = NULL, locus_id = NULL) {
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

  title <- if (!is.null(locus_id)) {
    locus_id
  } else {
    sprintf("chr%s:%g-%g", chr, start, end)
  }

  p <- ggplot(pts, aes(x = .data$pos, y = .data$y, colour = .data$trait)) +
    geom_point(size = 1.1, alpha = 0.7) +
    geom_hline(yintercept = gw_line, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(
      values = c("trait 1" = "#1f77b4", "trait 2" = "#d62728")
    ) +
    labs(
      x = "Position (bp)", y = expression(-log[10](p)),
      colour = NULL, title = title
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")

  # Mark the lead SNP's position if we can locate it in either trait.
  if (!is.null(lead_snp) && !is.na(lead_snp)) {
    hit <- pts$pos[pts$snp == lead_snp]
    if (length(hit) > 0) {
      p <- p + geom_vline(
        xintercept = hit[1], linetype = "dashed", colour = "grey20"
      )
    }
  }

  # Render headlessly to PNG.
  png(outfile, width = 1100, height = 750, res = 120)
  on.exit(dev.off(), add = TRUE)
  print(p)
  invisible(outfile)
}
