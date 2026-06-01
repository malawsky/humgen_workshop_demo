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
# both are supported (see HM_MAP): the older hm_* names and the
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
suppressPackageStartupMessages(library(fst))

STD_COLS <- c("snp", "chr", "pos", "ea", "oa", "maf", "beta", "se", "p")

# Bump this whenever read_harmonised's columns, types, or row filter change, so
# old .std.fst caches (keyed by accession) are not silently reused as stale.
STD_PARSE_VERSION <- 1L

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

#' Read a harmonised-format file into a standardised data.frame.
#'
#' Reads the harmonised columns into the standardised STD_COLS order and drops
#' variants with an unusable effect size or standard error. Handles both the
#' older hm_* column layout and the modern GWAS-SSF layout (see HM_MAP).
#' Reads with explicit colClasses + threading and renames/reorders in place, so
#' a second full copy of the 9 columns is never made (lower peak RAM, faster).
#' @param path  path to a harmonised "*.h.tsv.gz" file
#' @return data.frame with columns STD_COLS (one row per variant)
read_harmonised <- function(path) {
  # Read the header only (nrows = 0) to learn which candidate column names this
  # file uses, then re-read selecting just those. Harmonised files carry many
  # columns we never use; on a >1 GB file, parsing only the 9 we need is far
  # faster and uses far less memory than reading everything.
  header <- data.table::fread(path, nrows = 0)
  want <- vapply(STD_COLS, function(s) {
    hit <- HM_MAP[[s]][HM_MAP[[s]] %in% names(header)]
    if (length(hit) == 0) "" else hit[1]
  }, character(1))
  missing <- STD_COLS[want == ""]
  if (length(missing) > 0) {
    stop(sprintf(
      "harmonised file %s lacks column(s) for: %s; present: %s",
      path, paste(missing, collapse = ", "),
      paste(names(header), collapse = ", ")
    ))
  }

  # Force the right type per column at parse time so fread never guesses and we
  # avoid an as.*()-driven second copy of every column. chr stays character so
  # X/Y survive; pos as integer is safe (positions < 2^31). Keyed by the ACTUAL
  # present names (unname(want)), which is what fread sees in the file.
  classes <- c(
    snp = "character", chr = "character", pos = "integer",
    ea = "character", oa = "character", maf = "double",
    beta = "double", se = "double", p = "double"
  )
  col_classes <- setNames(classes[STD_COLS], unname(want))

  # Read once, only the 9 columns, using all available threads.
  raw <- data.table::fread(
    path,
    select = unname(want),
    colClasses = col_classes,
    nThread = data.table::getDTthreads()
  )

  # `want` is named by STD_COLS in priority order, so unname(want) aligns
  # positionally with STD_COLS; rename + reorder in place (no copy).
  data.table::setnames(raw, unname(want), STD_COLS)
  data.table::setcolorder(raw, STD_COLS)

  # Drop variants with an unusable effect size or standard error.
  keep <- is.finite(raw$beta) & raw$beta != 0 &
    is.finite(raw$se) & raw$se != 0
  raw <- raw[keep]

  # Downstream code expects a plain data.frame; columns already in STD_COLS.
  as.data.frame(raw)
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

  # Parsed-frame cache, keyed by accession (same policy as the .h.tsv.gz cache).
  # Parsing a full-genome harmonised file takes minutes; caching the already-
  # standardised frame as .fst lets re-runs skip both download AND parse.
  std_cache <- file.path(
    cache_dir, sprintf("%s.std.v%d.fst", study, STD_PARSE_VERSION)
  )
  if (file.exists(std_cache) && file.info(std_cache)$size > 0) {
    return(fst::read_fst(std_cache))
  }

  dest <- file.path(cache_dir, paste0(study, ".h.tsv.gz"))
  if (!(file.exists(dest) && file.info(dest)$size > 0)) {
    url <- resolve_harmonised_url(study)
    # Download to a .part file then rename, so a truncated download never
    # leaves a non-empty file the cache check would later trust.
    part <- paste0(dest, ".part")
    curl::curl_download(url, part)
    file.rename(part, dest)
  }
  df <- read_harmonised(dest)

  # Write parsed frame via temp-then-rename, so a crash mid-write never leaves a
  # truncated .std.fst that the cache check above would later trust.
  tmp <- paste0(std_cache, ".part")
  fst::write_fst(df, tmp)
  file.rename(tmp, std_cache)
  df
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
