# =============================================================================
# download_sumstats.R  --  resolve a GWAS Catalog study to harmonised sumstats
#                          and return a standardised data.frame.
#
# STATUS: local-file path implemented; accession download is deferred (step 7).
# =============================================================================
#
# Harmonised summary statistics live on the EBI FTP/HTTPS server, e.g.:
#   https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/
#       GCST004001-GCST005000/GCST004988/harmonised/
#       29059683-GCST004988-EFO_0000305.h.tsv.gz        (GRCh38, the ".h" file)
#
# The 1000-wide range bucket and the PMID/EFO in the filename make the path
# awkward to build by hand. Prefer one of:
#   * gwascatftp::get_file_paths() / get_harmonised_list()  to resolve the URL
#   * or list the .../harmonised/ directory and pick the "*.h.tsv.gz" file.
# Study metadata (sample size, trait, case/control counts) comes from
#   gwasrapidd::get_studies(study_id = ...).
#
# Harmonised -> standardised column mapping (this is fixed, put it in CLAUDE.md):
#   hm_rsid                       -> snp
#   hm_chrom                      -> chr
#   hm_pos                        -> pos
#   hm_effect_allele              -> ea
#   hm_other_allele               -> oa
#   hm_effect_allele_frequency    -> maf
#   hm_beta                       -> beta
#   standard_error                -> se
#   p_value                       -> p
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

STD_COLS <- c("snp", "chr", "pos", "ea", "oa", "maf", "beta", "se", "p")

#' Resolve a study accession to a harmonised sumstats URL.
#' @param accession e.g. "GCST004988"
#' @return a single URL string to the "*.h.tsv.gz" file.
resolve_harmonised_url <- function(accession) {
  stop("TODO: resolve the harmonised .h.tsv.gz URL for the accession")
}

# Harmonised column name -> standardised column name (fixed mapping, see CLAUDE.md).
HM_MAP <- c(
  snp  = "hm_rsid",
  chr  = "hm_chrom",
  pos  = "hm_pos",
  ea   = "hm_effect_allele",
  oa   = "hm_other_allele",
  maf  = "hm_effect_allele_frequency",
  beta = "hm_beta",
  se   = "standard_error",
  p    = "p_value"
)

#' Load sumstats for one trait as a standardised data.frame.
#'
#' Accepts EITHER a GWAS Catalog accession (downloads + caches to `cache_dir`)
#' OR a path to a local harmonised-format file (used directly, no download).
#' This dual mode is what lets the worked example run with no network.
#'
#' @param study      accession string OR local file path
#' @param cache_dir  where downloads are cached (skip re-download if present)
#' @return data.frame with columns STD_COLS (one row per variant)
load_sumstats <- function(study, cache_dir = "cache") {
  if (!file.exists(study)) {
    # TODO (step 7): resolve_harmonised_url(study) -> download to cache_dir.
    stop("accession download not implemented yet (see step 7)")
  }

  raw <- data.table::fread(study)

  # Pull the harmonised columns into the standardised STD_COLS order.
  df <- data.frame(
    snp  = as.character(raw[[HM_MAP[["snp"]]]]),
    chr  = as.character(raw[[HM_MAP[["chr"]]]]),   # character so X/Y survive
    pos  = as.integer(raw[[HM_MAP[["pos"]]]]),
    ea   = as.character(raw[[HM_MAP[["ea"]]]]),
    oa   = as.character(raw[[HM_MAP[["oa"]]]]),
    maf  = as.numeric(raw[[HM_MAP[["maf"]]]]),
    beta = as.numeric(raw[[HM_MAP[["beta"]]]]),
    se   = as.numeric(raw[[HM_MAP[["se"]]]]),
    p    = as.numeric(raw[[HM_MAP[["p"]]]]),        # p_value is a sci-notation string
    stringsAsFactors = FALSE
  )

  # Drop variants with an unusable effect size or standard error.
  keep <- is.finite(df$beta) & df$beta != 0 & is.finite(df$se) & df$se != 0
  df <- df[keep, , drop = FALSE]

  as.data.frame(df[, STD_COLS], stringsAsFactors = FALSE)
}
