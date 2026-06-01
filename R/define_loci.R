# =============================================================================
# define_loci.R  --  turn genome-wide significant variants into testable loci.
#
# STATUS: implemented.
# =============================================================================
#
# Approach (keep it simple, no LD reference needed for the demo):
#   1. Find variants with p < p_threshold (default 5e-8) in trait 1 and/or
#      trait 2 depending on `sig_mode`.
#   2. Window each hit by +/- `window` bp and MERGE overlapping windows into
#      loci (per chromosome).
#   3. Return a data.frame of loci: locus_id, chr, start, end.
#
# `sig_mode` (set by the --sig-mode CLI flag):
#   "either" - locus is significant in trait 1 OR trait 2   (default)
#   "both"   - significant in BOTH traits
#   "1"      - significant in trait 1
#   "2"      - significant in trait 2
# -----------------------------------------------------------------------------

#' @param ss1,ss2    standardised sumstats data.frames (cols incl. chr,pos,p)
#' @param sig_mode   one of "either","both","1","2"
#' @param p_threshold genome-wide significance threshold
#' @param window     half-width in bp around each lead hit
#' @return data.frame(locus_id, chr, start, end)
define_loci <- function(ss1, ss2,
                        sig_mode = "either",
                        p_threshold = 5e-8,
                        window = 5e5) {
  sig_mode <- match.arg(as.character(sig_mode), c("either", "both", "1", "2"))

  # An empty result has exactly these four columns; reuse it everywhere.
  empty <- data.frame(
    locus_id = character(0), chr = character(0),
    start = numeric(0), end = numeric(0),
    stringsAsFactors = FALSE
  )

  # Significant hits in one trait: data.frame(chr, pos), finite p < threshold.
  sig_hits <- function(ss) {
    keep <- is.finite(ss$p) & ss$p < p_threshold
    data.frame(
      chr = as.character(ss$chr[keep]),
      pos = as.numeric(ss$pos[keep]),
      stringsAsFactors = FALSE
    )
  }

  h1 <- sig_hits(ss1)
  h2 <- sig_hits(ss2)

  # Build the candidate hit set, tagging each hit with its trait of origin.
  if (sig_mode == "1") {
    hits <- cbind(h1, trait = rep("1", nrow(h1)), stringsAsFactors = FALSE)
  } else if (sig_mode == "2") {
    hits <- cbind(h2, trait = rep("2", nrow(h2)), stringsAsFactors = FALSE)
  } else {
    # "either" / "both": union of hits from both traits.
    hits <- rbind(
      cbind(h1, trait = rep("1", nrow(h1)), stringsAsFactors = FALSE),
      cbind(h2, trait = rep("2", nrow(h2)), stringsAsFactors = FALSE)
    )
  }

  if (nrow(hits) == 0) {
    return(empty)
  }

  # Window each hit. Numeric arithmetic avoids integer overflow at large pos.
  hits$start <- pmax(1, hits$pos - window)
  hits$end <- hits$pos + window

  # Merge overlapping windows per chromosome (chr is character).
  merged_list <- lapply(split(hits, hits$chr), function(chr_hits) {
    chr_hits <- chr_hits[order(chr_hits$start), , drop = FALSE]

    chr <- chr_hits$chr[1]
    starts <- chr_hits$start
    ends <- chr_hits$end
    traits <- chr_hits$trait

    # Sweep-merge: extend the current locus while the next window overlaps it.
    cur_start <- starts[1]
    cur_end <- ends[1]
    cur_traits <- traits[1]
    out_start <- out_end <- numeric(0)
    out_traits <- character(0)

    for (i in seq_along(starts)[-1]) {
      if (starts[i] <= cur_end) {
        cur_end <- max(cur_end, ends[i])
        cur_traits <- paste(unique(c(cur_traits, traits[i])), collapse = "")
      } else {
        out_start <- c(out_start, cur_start)
        out_end <- c(out_end, cur_end)
        out_traits <- c(out_traits, cur_traits)
        cur_start <- starts[i]
        cur_end <- ends[i]
        cur_traits <- traits[i]
      }
    }
    out_start <- c(out_start, cur_start)
    out_end <- c(out_end, cur_end)
    out_traits <- c(out_traits, cur_traits)

    data.frame(
      chr = chr, start = out_start, end = out_end,
      traits = out_traits, stringsAsFactors = FALSE
    )
  })

  loci <- do.call(rbind, merged_list)
  rownames(loci) <- NULL

  # "both": keep loci with a significant variant from BOTH traits (regional).
  if (sig_mode == "both") {
    has_both <- grepl("1", loci$traits) & grepl("2", loci$traits)
    loci <- loci[has_both, , drop = FALSE]
  }

  if (nrow(loci) == 0) {
    return(empty)
  }

  # Order by chromosome (numeric chrs numerically; X/Y/MT after), then start.
  chr_num <- suppressWarnings(as.numeric(loci$chr))
  # Non-numeric chrs (X, Y, MT) sort after numeric ones, alphabetically.
  chr_rank <- ifelse(is.na(chr_num), Inf, chr_num)
  ord <- order(chr_rank, loci$chr, loci$start)
  loci <- loci[ord, , drop = FALSE]

  data.frame(
    locus_id = sprintf(
      "chr%s:%.1f-%.1fMb", loci$chr, loci$start / 1e6, loci$end / 1e6
    ),
    chr = loci$chr,
    start = loci$start,
    end = loci$end,
    stringsAsFactors = FALSE
  )
}

#' Subset a standardised sumstats data.frame to one locus region.
extract_region <- function(ss, chr, start, end) {
  ss[ss$chr == chr & ss$pos >= start & ss$pos <= end, , drop = FALSE]
}
