# =============================================================================
# define_loci.R  --  turn genome-wide significant variants into testable loci.
#
# STATUS: STUB. Implement the TODOs during the build.
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
  stop("TODO: select significant variants per sig_mode, window and merge into loci")
}

#' Subset a standardised sumstats data.frame to one locus region.
extract_region <- function(ss, chr, start, end) {
  ss[ss$chr == chr & ss$pos >= start & ss$pos <= end, , drop = FALSE]
}
