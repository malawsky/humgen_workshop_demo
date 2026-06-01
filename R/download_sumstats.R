# =============================================================================
# download_sumstats.R  --  resolve a GWAS Catalog study to harmonised sumstats
#                          and return a standardised data.frame.
#
# STATUS: local-file path implemented; accession download IMPLEMENTED over
#         HTTPS (we list the harmonised dir and pick the "*.h.tsv.gz" file).
#         We do NOT use gwascatftp: it needs the lftp binary, which is absent.
# =============================================================================
#
# Harmonised summary statistics live on the EBI FTP/HTTPS server, e.g.:
#   https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/
#       GCST004001-GCST005000/GCST004988/harmonised/
#       29059683-GCST004988-EFO_0000305.h.tsv.gz        (GRCh38, the ".h" file)
#
# The 1000-wide range bucket and the PMID/EFO in the filename make the path
# awkward to build by hand, so we compute the bucket (accession_bucket) and
# list the .../harmonised/ directory to pick the single "*.h.tsv.gz" file.
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
suppressPackageStartupMessages(library(curl))
suppressPackageStartupMessages(library(gwasrapidd))

STD_COLS <- c("snp", "chr", "pos", "ea", "oa", "maf", "beta", "se", "p")

# Compute the 1000-wide range bucket that prefixes the accession on the FTP
# tree, preserving the accession's zero-padding width. For example accession
# GCST004988 falls in bucket "GCST004001-GCST005000", and GCST90012345 falls
# in bucket "GCST90012001-GCST90013000".
accession_bucket <- function(accession) {
  n <- sub("^GCST", "", accession)
  width <- nchar(n)
  num <- as.numeric(n) # numeric, not integer: large accessions overflow int
  lo <- (num - 1) %/% 1000 * 1000 + 1
  hi <- lo + 999
  pad <- function(x) formatC(x, width = width, flag = "0", format = "d")
  sprintf("GCST%s-GCST%s", pad(lo), pad(hi))
}

#' Resolve a study accession to a harmonised sumstats URL.
#'
#' Builds the bucketed harmonised directory URL, lists it over HTTPS, and
#' returns the single "*.h.tsv.gz" file found there.
#' @param accession e.g. "GCST004988"
#' @return a single URL string to the "*.h.tsv.gz" file.
resolve_harmonised_url <- function(accession) {
  base <- "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
  dir_url <- sprintf(
    "%s/%s/%s/harmonised/", base, accession_bucket(accession), accession
  )
  resp <- curl::curl_fetch_memory(dir_url)
  if (resp$status_code >= 400) {
    stop(sprintf(
      "no harmonised directory for %s (HTTP %d): %s",
      accession, resp$status_code, dir_url
    ))
  }
  html <- rawToChar(resp$content)
  # The regex captures only up to ".gz", so a "*.h.tsv.gz.tbi" index yields the
  # same string as its data file; unique() then collapses the two to one entry.
  pattern <- "[A-Za-z0-9._-]+[.]h[.]tsv[.]gz"
  files <- unique(regmatches(html, gregexpr(pattern, html))[[1]])
  if (length(files) == 0) {
    stop(sprintf(
      "no .h.tsv.gz harmonised file found for %s at %s", accession, dir_url
    ))
  }
  if (length(files) > 1) {
    stop(sprintf(
      "multiple .h.tsv.gz files for %s: %s",
      accession, paste(files, collapse = ", ")
    ))
  }
  paste0(dir_url, files[1])
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

#' Read a harmonised-format file into a standardised data.frame.
#'
#' Reads the harmonised columns into the standardised STD_COLS order and drops
#' variants with an unusable effect size or standard error.
#' @param path  path to a harmonised "*.h.tsv.gz" file
#' @return data.frame with columns STD_COLS (one row per variant)
read_harmonised <- function(path) {
  raw <- data.table::fread(path)

  # Pull the harmonised columns into the standardised STD_COLS order.
  df <- data.frame(
    snp = as.character(raw[[HM_MAP[["snp"]]]]),
    chr = as.character(raw[[HM_MAP[["chr"]]]]), # character so X/Y survive
    pos = as.integer(raw[[HM_MAP[["pos"]]]]),
    ea = as.character(raw[[HM_MAP[["ea"]]]]),
    oa = as.character(raw[[HM_MAP[["oa"]]]]),
    maf = as.numeric(raw[[HM_MAP[["maf"]]]]),
    beta = as.numeric(raw[[HM_MAP[["beta"]]]]),
    se = as.numeric(raw[[HM_MAP[["se"]]]]),
    p = as.numeric(raw[[HM_MAP[["p"]]]]), # p_value is a sci-notation string
    stringsAsFactors = FALSE
  )

  # Drop variants with an unusable effect size or standard error.
  keep <- is.finite(df$beta) & df$beta != 0 & is.finite(df$se) & df$se != 0
  df <- df[keep, , drop = FALSE]

  as.data.frame(df[, STD_COLS], stringsAsFactors = FALSE)
}

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
  if (file.exists(study)) {
    return(read_harmonised(study))
  }

  # Accession: download (once) to the cache, then read from there.
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  dest <- file.path(cache_dir, paste0(study, ".h.tsv.gz"))
  if (!(file.exists(dest) && file.info(dest)$size > 0)) {
    url <- resolve_harmonised_url(study)
    # Download to a .part file then rename, so a truncated download never
    # leaves a non-empty file the cache check would later trust.
    part <- paste0(dest, ".part")
    curl::curl_download(url, part)
    file.rename(part, dest)
  }
  read_harmonised(dest)
}

#' Resolve a study accession to its sample size N.
#'
#' N is the total initial-ancestry sample size, summed across all "initial"
#' ancestry rows (handles multi-ancestry studies).
#' @param accession e.g. "GCST004988"
#' @return numeric N (> 0)
resolve_study_n <- function(accession) {
  st <- gwasrapidd::get_studies(study_id = accession)
  anc <- st@ancestries
  if (is.null(anc) || nrow(anc) == 0) {
    stop(sprintf(
      "no ancestry metadata for %s; pass --n1/--n2 explicitly", accession
    ))
  }
  init <- anc[anc$type == "initial", , drop = FALSE]
  if (nrow(init) == 0) {
    stop(sprintf(
      "no 'initial' sample size for %s; pass --n1/--n2 explicitly", accession
    ))
  }
  n <- sum(init$number_of_individuals, na.rm = TRUE)
  if (!is.finite(n) || n <= 0) {
    stop(sprintf(
      "resolved a non-positive N for %s; pass --n1/--n2 explicitly", accession
    ))
  }
  n
}
