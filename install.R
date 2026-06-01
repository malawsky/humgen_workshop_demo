#!/usr/bin/env Rscript
# Install the R packages coloc-pair needs. Run once: Rscript install.R
# Everything here is on CRAN -- no Bioconductor, no system tools required.
repos <- "https://cloud.r-project.org"

cran <- c(
  "optparse",     # CLI argument parsing
  "data.table",   # fast sumstats I/O
  "R.utils",      # lets data.table::fread read gzipped (.h.tsv.gz) files
  "coloc",        # colocalisation (coloc.abf)
  "gwasrapidd",   # GWAS Catalog REST API (study metadata / sample size)
  "curl",         # resolve + download harmonised sumstats over HTTPS
  "ggplot2",      # Miami and locus-zoom plots
  "knitr"         # image_uri: embed plots inline in the HTML report
)
to_get <- setdiff(cran, rownames(installed.packages()))
if (length(to_get)) install.packages(to_get, repos = repos)

cat("Done. Sanity check:\n")
for (p in cran) {
  cat(sprintf("  %-12s %s\n", p,
              tryCatch(as.character(packageVersion(p)),
                       error = function(e) "MISSING")))
}
