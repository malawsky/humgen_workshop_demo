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
# Argument parsing is complete. The PIPELINE below is the build target.
# =============================================================================

suppressPackageStartupMessages(library(optparse))

opts <- list(
  make_option("--study1", type = "character", help = "Trait 1: GCST accession OR local harmonised file"),
  make_option("--study2", type = "character", help = "Trait 2: GCST accession OR local harmonised file"),
  make_option("--outdir", type = "character", default = "results/run", help = "Output directory [%default]"),
  make_option("--sig-mode", type = "character", default = "either",
              help = "Define loci significant in: either|both|1|2 [%default]", dest = "sig_mode"),
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
  make_option("--keep-sumstats", action = "store_true", default = FALSE, dest = "keep_sumstats", help = "Keep downloaded sumstats")
)
args <- parse_args(OptionParser(option_list = opts))
if (is.null(args$study1) || is.null(args$study2)) stop("--study1 and --study2 are required")

# Load the project's R modules (paths relative to repo root).
here <- function(...) file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), "..", ...)
source(here("R", "coloc_wrapper.R"))     # TRUSTED - do not modify
source(here("R", "download_sumstats.R"))
source(here("R", "define_loci.R"))
source(here("R", "plots.R"))

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# --- PIPELINE (build target) -------------------------------------------------
# 1. Load both traits' sumstats (download + cache, or read local file).
#       ss1 <- load_sumstats(args$study1, cache_dir = file.path(args$outdir, "cache"))
#       ss2 <- load_sumstats(args$study2, ...)
#       Resolve n1/n2 from study metadata if not supplied.
#
# 2. Define significant loci.
#       loci <- define_loci(ss1, ss2, sig_mode = args$sig_mode,
#                           p_threshold = args$p_threshold, window = args$window)
#
# 3. For each locus: extract the region from both traits, run coloc via the
#    TRUSTED wrapper, and collect one summary row per locus.
#       res  <- run_coloc_abf(extract_region(ss1, ...), extract_region(ss2, ...),
#                             type1 = args$type1, type2 = args$type2,
#                             n1 = n1, n2 = n2, s1 = args$s1, s2 = args$s2)
#       row  <- summarise_coloc(res, locus_label)
#    Write the combined table to file.path(args$outdir, "coloc_results.tsv").
#
# 4. Miami plot of both traits, colocalising loci highlighted ->
#       file.path(args$outdir, "miami.png")
#
# 5. Locus zooms for the top-N loci by PP.H4 ->
#       file.path(args$outdir, "locuszoom_<locus>.png")
#
# 6. Print a short summary to the console: how many loci tested, how many
#    colocalise at PP.H4 >= args$pp4_threshold, and where the outputs are.
# -----------------------------------------------------------------------------

stop("PIPELINE NOT YET IMPLEMENTED -- this is the build target. See SPEC.md.")
