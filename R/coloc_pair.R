#!/usr/bin/env Rscript
# =============================================================================
# coloc_pair.R  --  CLI: colocalise two GWAS Catalog studies.
#
# Usage:
#   Rscript R/coloc_pair.R \
#     --study1 GCST90012345 --type1 cc   --s1 0.2 \
#     --study2 GCST90067890 --type2 quant \
#     --sig-mode either --outdir results/run1
#
# The worked example (no network needed):
#   Rscript R/coloc_pair.R \
#     --study1 data/example/trait1.h.tsv.gz --type1 quant \
#     --study2 data/example/trait2.h.tsv.gz --type2 quant \
#     --outdir results/example
#
# Argument parsing is complete; the pipeline below loads, defines loci, runs
# coloc, and writes the results table (plots are wired in later steps).
# =============================================================================

suppressPackageStartupMessages(library(optparse))

opts <- list(
  make_option("--study1", type = "character", help = "Trait 1: GCST accession OR local harmonised file"),
  make_option("--study2", type = "character", help = "Trait 2: GCST accession OR local harmonised file"),
  make_option("--outdir", type = "character", default = "results/run", help = "Output directory [%default]"),
  make_option("--sig-mode",
    type = "character", default = "either",
    help = "Define loci significant in: either|both|1|2 [%default]", dest = "sig_mode"
  ),
  make_option("--type1", type = "character", default = "quant", help = "Trait 1 type: quant|cc [%default]"),
  make_option("--type2", type = "character", default = "quant", help = "Trait 2 type: quant|cc [%default]"),
  make_option("--s1", type = "double", default = NA, help = "Trait 1 case proportion (if cc)"),
  make_option("--s2", type = "double", default = NA, help = "Trait 2 case proportion (if cc)"),
  make_option("--n1", type = "double", default = NA, help = "Trait 1 sample size (else from study metadata)"),
  make_option("--n2", type = "double", default = NA, help = "Trait 2 sample size (else from study metadata)"),
  make_option("--window", type = "double", default = 5e5, help = "Locus half-width in bp [%default]"),
  make_option("--p-threshold", type = "double", default = 5e-8, dest = "p_threshold", help = "GWAS sig threshold [%default]"),
  make_option("--pp4-threshold", type = "double", default = 0.8, dest = "pp4_threshold", help = "PP.H4 to call a coloc [%default]"),
  make_option("--top-n", type = "integer", default = 5, dest = "top_n", help = "Locus zooms for top-N loci by PP.H4 [%default]"),
  # NOTE: cached sumstats are currently always kept under <outdir>/cache;
  # this flag is reserved and performs no cleanup (no deletion is done).
  make_option("--keep-sumstats", action = "store_true", default = FALSE, dest = "keep_sumstats", help = "Keep downloaded sumstats")
)
args <- parse_args(OptionParser(option_list = opts))
if (is.null(args$study1) || is.null(args$study2)) stop("--study1 and --study2 are required")

# Load the project's R modules (paths relative to repo root).
here <- function(...) file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "..", ...)
source(here("R", "coloc_wrapper.R")) # TRUSTED - do not modify
source(here("R", "download_sumstats.R"))
source(here("R", "define_loci.R"))
source(here("R", "plots.R"))
source(here("R", "report.R"))

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# Column names of an (empty) results table, taken from summarise_coloc().
RESULT_COLS <- c(
  "locus", "nsnps", "PP.H0", "PP.H1", "PP.H2",
  "PP.H3", "PP.H4", "lead_snp"
)
results_path <- file.path(args$outdir, "coloc_results.tsv")

# Write a header-only results table (used when there is nothing to report).
write_empty_results <- function() {
  empty <- as.data.frame(setNames(rep(list(character(0)), length(RESULT_COLS)), RESULT_COLS))
  write.table(empty, results_path, row.names = FALSE, quote = FALSE, sep = "\t")
  message("Results table: ", results_path)
}

# Resolve a trait's sample size: use --nX if given, else fail with guidance.
# Local files carry no sample-size metadata, so N must be supplied explicitly.
resolve_n <- function(n, study, flag) {
  if (!is.na(n)) {
    return(as.numeric(n))
  }
  if (file.exists(study)) {
    stop(sprintf("local file '%s' has no sample-size metadata; pass %s", study, flag))
  }
  return(resolve_study_n(study))
}

# --- PIPELINE ----------------------------------------------------------------

# 1. Load both traits' sumstats (local file read, or cached download).
ss1 <- load_sumstats(args$study1, cache_dir = file.path(args$outdir, "cache"))
ss2 <- load_sumstats(args$study2, cache_dir = file.path(args$outdir, "cache"))

# 2. Resolve sample sizes and case proportions for each trait.
n1 <- resolve_n(args$n1, args$study1, "--n1")
n2 <- resolve_n(args$n2, args$study2, "--n2")
s1 <- if (args$type1 == "cc" && !is.na(args$s1)) args$s1 else NULL
s2 <- if (args$type2 == "cc" && !is.na(args$s2)) args$s2 else NULL

# 3. Define significant loci shared/tested across the two traits.
loci <- define_loci(ss1, ss2,
  sig_mode = args$sig_mode,
  p_threshold = args$p_threshold, window = args$window
)

if (nrow(loci) == 0) {
  message("No significant loci found; writing an empty results table.")
  write_empty_results()
  quit(save = "no", status = 0)
}

# Index both frames once by (chr, pos) so each per-locus extract_region call is
# a keyed binary-search lookup instead of a full ~tens-of-millions-row scan.
# setkey sorts in place (no full copy); define_loci already ran and the only
# later consumers (extract_region, make_miami) do not depend on row order.
data.table::setDT(ss1)
data.table::setkey(ss1, chr, pos)
data.table::setDT(ss2)
data.table::setkey(ss2, chr, pos)

# 4. Run coloc on each locus, collecting one summary row each.
rows <- list()
for (i in seq_len(nrow(loci))) {
  locus <- loci[i, ]
  r1 <- extract_region(ss1, locus$chr, locus$start, locus$end)
  r2 <- extract_region(ss2, locus$chr, locus$start, locus$end)
  res <- run_coloc_abf(r1, r2,
    type1 = args$type1, type2 = args$type2,
    n1 = n1, n2 = n2, s1 = s1, s2 = s2
  )
  row <- summarise_coloc(res, locus$locus_id)
  if (!is.null(row)) rows[[length(rows) + 1]] <- row
}
if (length(rows) == 0) {
  message("No locus had enough overlapping variants to test; writing an empty table.")
  write_empty_results()
  quit(save = "no", status = 0)
}
results <- do.call(rbind, rows)

# 5. Sort by PP.H4 (strongest colocalisation first) and write the table.
results <- results[order(-results$PP.H4), , drop = FALSE]
write.table(results, results_path, row.names = FALSE, quote = FALSE, sep = "\t")

# Miami plot: attach each locus's PP.H4 (matched by id) and draw it. Plot
# failures warn but never abort the run.
loci$PP.H4 <- results$PP.H4[match(loci$locus_id, results$locus)]
miami_path <- file.path(args$outdir, "miami.png")
tryCatch(
  make_miami(ss1, ss2, loci, miami_path, args$pp4_threshold),
  error = function(e) warning("Miami plot failed: ", conditionMessage(e))
)

# Locus-zoom plots: top-N loci by PP.H4 (results is already sorted desc). Map
# each back to its coordinates via `loci` (locus_id == locus). Failures warn but
# never abort the run.
zoom_paths <- character(0)
tryCatch(
  {
    top <- head(results, args$top_n)
    for (i in seq_len(nrow(top))) {
      lid <- top$locus[i]
      lc <- loci[loci$locus_id == lid, , drop = FALSE]
      if (nrow(lc) == 0) next
      safe <- gsub("[^A-Za-z0-9]+", "_", lid)
      zoom_path <- file.path(args$outdir, sprintf("locuszoom_%s.png", safe))
      make_locuszoom(ss1, ss2, lc$chr[1], lc$start[1], lc$end[1], zoom_path,
        lead_snp = top$lead_snp[i], locus_id = lid
      )
      if (file.exists(zoom_path)) zoom_paths <- c(zoom_paths, zoom_path)
    }
  },
  error = function(e) warning("Locus-zoom plots failed: ", conditionMessage(e))
)

# Self-contained HTML report (final output). A failing report warns but never
# aborts the run, matching the Miami / locus-zoom blocks above.
report_path <- file.path(args$outdir, "report.html")
tryCatch(
  make_report(args$outdir, results, loci, args, n1, n2, miami_path, zoom_paths, report_path),
  error = function(e) warning("HTML report failed: ", conditionMessage(e))
)

# 6. Console summary.
n_coloc <- sum(results$PP.H4 >= args$pp4_threshold)
message(sprintf("Loci tested: %d", nrow(results)))
message(sprintf("Colocalising (PP.H4 >= %.2f): %d", args$pp4_threshold, n_coloc))
message("Results table: ", results_path)
if (file.exists(miami_path)) message("Miami plot: ", miami_path)
message(sprintf("Locus-zoom plots written: %d", length(zoom_paths)))
if (file.exists(report_path)) message("HTML report: ", report_path)
