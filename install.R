#!/usr/bin/env Rscript
# Install the R packages coloc-pair needs. Run once: Rscript install.R
repos <- "https://cloud.r-project.org"

cran <- c(
  "optparse",     # CLI argument parsing
  "data.table",   # fast sumstats I/O
  "coloc",        # colocalisation (coloc.abf)
  "topr",         # Manhattan / Miami plots
  "locuszoomr",   # regional locus-zoom plots
  "gwasrapidd",   # GWAS Catalog REST API (study metadata, top associations)
  "remotes"       # to install gwascatftp from GitHub
)
to_get <- setdiff(cran, rownames(installed.packages()))
if (length(to_get)) install.packages(to_get, repos = repos)

# gwascatftp (resolve harmonised file paths on the EBI FTP server) is on GitHub.
if (!requireNamespace("gwascatftp", quietly = TRUE)) {
  remotes::install_github("cfbeuchel/gwascatftp")
}

# Gene tracks for locus-zoom plots (GRCh38) come from Bioconductor.
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = repos)
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
  BiocManager::install("EnsDb.Hsapiens.v86", update = FALSE, ask = FALSE)
}

cat("Done. Sanity check:\n")
for (p in c("coloc", "topr", "locuszoomr", "data.table", "optparse")) {
  cat(sprintf("  %-12s %s\n", p,
              tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")))
}
