#!/usr/bin/env Rscript
# =============================================================================
# test_accession_online.R  --  ONLINE test for the accession download path.
#
# Exercises ONLY metadata / directory-listing requests (resolve_harmonised_url
# and resolve_study_n). It NEVER downloads the full multi-hundred-MB sumstats
# file. Pure offline checks (accession_bucket) run first and always; when the
# network is unavailable the test prints SKIP and exits without the online
# assertions, so offline runs still pass.
#
# Run from the repo root:   Rscript tests/test_accession_online.R
# =============================================================================

source("R/download_sumstats.R")

fail <- 0L
check <- function(cond, msg) {
  cat(ifelse(isTRUE(cond), "PASS  ", "FAIL  "), msg, "\n", sep = "")
  if (!isTRUE(cond)) fail <<- fail + 1L
}

# --- 1. Pure offline checks (no network) ------------------------------------
check(
  accession_bucket("GCST004988") == "GCST004001-GCST005000",
  "accession_bucket(GCST004988) == GCST004001-GCST005000"
)
check(
  accession_bucket("GCST90012345") == "GCST90012001-GCST90013000",
  "accession_bucket(GCST90012345) == GCST90012001-GCST90013000"
)

# --- 2. Connectivity probe --------------------------------------------------
# Any error or a >=400 status means we are effectively offline: skip the
# online assertions (the offline checks above have already run).
options(timeout = 30)
online <- tryCatch(
  {
    probe <- curl::curl_fetch_memory(
      "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
    )
    probe$status_code < 400
  },
  error = function(e) FALSE
)
if (!isTRUE(online)) {
  cat("SKIP: offline (network unavailable)\n")
  quit(status = if (fail == 0L) 0L else 1L)
}

# --- 3. Online assertions (metadata / listing only) -------------------------
check(
  endsWith(
    resolve_harmonised_url("GCST004988"),
    "29059683-GCST004988-EFO_0000305.h.tsv.gz"
  ),
  "resolve_harmonised_url(GCST004988) ends with the expected .h.tsv.gz file"
)
check(
  resolve_study_n("GCST004988") == 139274,
  "resolve_study_n(GCST004988) == 139274"
)

cat(if (fail == 0L) "\nAll accession-online checks passed.\n" else sprintf("\n%d check(s) failed.\n", fail))
quit(status = if (fail == 0L) 0L else 1L)
