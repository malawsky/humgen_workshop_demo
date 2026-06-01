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
# Harmonised -> standardised column mapping. Two layouts exist in the wild and
# both are supported (see HM_MAP / pick_col): the older hm_* names and the
# modern GWAS-SSF names. We try the hm_* name first, then the plain name:
#   hm_rsid                    | rsid                     -> snp
#   hm_chrom                   | chromosome               -> chr
#   hm_pos                     | base_pair_location       -> pos
#   hm_effect_allele           | effect_allele            -> ea
#   hm_other_allele            | other_allele             -> oa
#   hm_effect_allele_frequency | effect_allele_frequency  -> maf
#   hm_beta                    | beta                     -> beta
#   standard_error                                        -> se
#   p_value                                               -> p
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

# Standardised column -> candidate harmonised column names, in priority order.
# Older harmonised files use the hm_* names; modern GWAS-SSF files use the plain
# names. We try each candidate in turn so both layouts load (see CLAUDE.md).
HM_MAP <- list(
  snp  = c("hm_rsid", "rsid"),
  chr  = c("hm_chrom", "chromosome"),
  pos  = c("hm_pos", "base_pair_location"),
  ea   = c("hm_effect_allele", "effect_allele"),
  oa   = c("hm_other_allele", "other_allele"),
  maf  = c("hm_effect_allele_frequency", "effect_allele_frequency"),
  beta = c("hm_beta", "beta"),
  se   = c("standard_error", "hm_standard_error"),
  p    = c("p_value", "hm_p_value")
)

# Return the first candidate column present in `raw`, or stop with guidance.
pick_col <- function(raw, std_name) {
  hit <- HM_MAP[[std_name]][HM_MAP[[std_name]] %in% names(raw)]
  if (length(hit) == 0) {
    stop(sprintf(
      "harmonised file has no column for '%s' (looked for: %s); present: %s",
      std_name, paste(HM_MAP[[std_name]], collapse = ", "),
      paste(names(raw), collapse = ", ")
    ))
  }
  raw[[hit[1]]]
}

#' Read a harmonised-format file into a standardised data.frame.
#'
#' Reads the harmonised columns into the standardised STD_COLS order and drops
#' variants with an unusable effect size or standard error. Handles both the
#' older hm_* column layout and the modern GWAS-SSF layout (see HM_MAP).
#' @param path  path to a harmonised "*.h.tsv.gz" file
#' @return data.frame with columns STD_COLS (one row per variant)
read_harmonised <- function(path) {
  raw <- data.table::fread(path)

  # Pull the harmonised columns into the standardised STD_COLS order.
  df <- data.frame(
    snp = as.character(pick_col(raw, "snp")),
    chr = as.character(pick_col(raw, "chr")), # character so X/Y survive
    pos = as.integer(pick_col(raw, "pos")),
    ea = as.character(pick_col(raw, "ea")),
    oa = as.character(pick_col(raw, "oa")),
    maf = as.numeric(pick_col(raw, "maf")),
    beta = as.numeric(pick_col(raw, "beta")),
    se = as.numeric(pick_col(raw, "se")),
    p = as.numeric(pick_col(raw, "p")), # p_value is a sci-notation string
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
