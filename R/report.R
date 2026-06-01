# =============================================================================
# report.R  --  the self-contained HTML report (the final output).
#
# STATUS: implemented. No pandoc / rmarkdown::render needed -- the HTML is built
# directly as a string and written out, with plot PNGs embedded INLINE as
# base64 data URIs (knitr::image_uri) so report.html is fully self-contained
# (one clickable file, no external assets).
# =============================================================================
#
# The report gathers, for one coloc-pair run:
#   1. a header block of the run inputs (traits, types, Ns, thresholds),
#   2. a "Colocalising loci" table (the subset with PP.H4 >= threshold),
#   3. an "All tested loci" table (the full results, coloc rows highlighted),
#   4. the Miami plot and each per-locus locus-zoom plot, embedded inline.
#
# It is deliberately dependency-light (base R + knitr for image_uri) and never
# errors on empty/degenerate inputs: at worst it writes a minimal valid page.
# -----------------------------------------------------------------------------

#' Build a self-contained HTML report for one coloc-pair run.
#'
#' Writes `report.html` at `outfile`, embedding the Miami and locus-zoom PNGs
#' inline as base64 data URIs so the file is fully portable. Tolerates 0-row
#' `results`, empty `zoom_paths`, and a missing `miami_path`.
#'
#' @param outdir      output directory of the run (recorded for reference)
#' @param results     data.frame of coloc summaries (sorted by PP.H4 desc),
#'                    cols incl. locus,nsnps,PP.H0..PP.H4,lead_snp; can be empty
#' @param loci        data.frame of loci (unused for layout, kept for context)
#' @param args        parsed CLI args (study1/2, type1/2, sig_mode, window,
#'                    p_threshold, pp4_threshold, ...)
#' @param n1,n2       resolved sample sizes (NOT args$n1/args$n2)
#' @param miami_path  path to the Miami PNG (may be NULL or not exist)
#' @param zoom_paths  character vector of locus-zoom PNG paths (may be empty)
#' @param outfile     path to write report.html
#' @return invisibly, the outfile path
make_report <- function(outdir, results, loci, args, n1, n2,
                        miami_path, zoom_paths, outfile) {
  # --- helpers ---------------------------------------------------------------

  # Escape text interpolated into HTML (study ids, locus ids, etc.).
  esc <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }

  # Format a possibly-missing scalar for display.
  fmt <- function(x) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return("NA")
    }
    esc(x)
  }

  # Render a data.frame as an HTML <table>. `row_class` is an optional
  # character vector (one per row, "" for none) applied to each <tr>.
  html_table <- function(df, row_class = NULL) {
    if (is.null(df) || nrow(df) == 0) {
      return("<p>Nothing to report.</p>")
    }
    if (is.null(row_class)) row_class <- rep("", nrow(df))

    head_cells <- paste0("<th>", esc(names(df)), "</th>", collapse = "")
    header <- paste0("<thead><tr>", head_cells, "</tr></thead>")

    body_rows <- vapply(seq_len(nrow(df)), function(i) {
      cells <- paste0(
        "<td>",
        vapply(df[i, ], function(v) esc(format(v)), character(1)),
        "</td>",
        collapse = ""
      )
      cls <- row_class[i]
      tr_open <- if (nzchar(cls)) {
        sprintf("<tr class=\"%s\">", esc(cls))
      } else {
        "<tr>"
      }
      paste0(tr_open, cells, "</tr>")
    }, character(1))

    paste0(
      "<table>", header, "<tbody>",
      paste(body_rows, collapse = ""), "</tbody></table>"
    )
  }

  # Embed a PNG inline as a base64 data URI; return "" if it can't be read.
  embed_png <- function(path, caption) {
    if (is.null(path) || !is.character(path) || length(path) != 1) {
      return("")
    }
    if (!file.exists(path)) {
      return("")
    }
    uri <- tryCatch(knitr::image_uri(path), error = function(e) NULL)
    if (is.null(uri)) {
      return("")
    }
    paste0(
      "<h2>", esc(caption), "</h2>\n",
      "<img src=\"", uri, "\" style=\"max-width:100%\" alt=\"",
      esc(caption), "\">\n"
    )
  }

  if (is.null(results)) {
    results <- data.frame()
  }
  pp4_thr <- args$pp4_threshold

  # --- inline style ----------------------------------------------------------

  style <- paste0(
    "<style>",
    "body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;",
    "margin:2rem auto;max-width:1100px;color:#222;line-height:1.4;}",
    "h1{font-size:1.6rem;}h2{font-size:1.2rem;margin-top:2rem;}",
    "table{border-collapse:collapse;margin:0.75rem 0;width:100%;}",
    "th,td{border:1px solid #ccc;padding:6px 10px;text-align:left;",
    "font-size:0.9rem;}",
    "th{background:#f2f2f2;}",
    "tr.coloc{background:#fff4cc;}",
    "dl.meta{display:grid;grid-template-columns:max-content 1fr;",
    "gap:2px 1rem;}dl.meta dt{font-weight:600;}dl.meta dd{margin:0;}",
    "img{border:1px solid #eee;margin:0.5rem 0;}",
    "</style>"
  )

  # --- metadata / header block ----------------------------------------------

  meta_items <- list(
    c("Trait 1 (study1)", fmt(args$study1)),
    c("Trait 2 (study2)", fmt(args$study2)),
    c("Type 1", fmt(args$type1)),
    c("Type 2", fmt(args$type2)),
    c("Sample size N1", fmt(n1)),
    c("Sample size N2", fmt(n2)),
    c("sig-mode", fmt(args$sig_mode)),
    c("Window (bp)", fmt(args$window)),
    c("p-threshold", fmt(args$p_threshold)),
    c("PP.H4 threshold", fmt(pp4_thr)),
    c("Date generated", fmt(as.character(Sys.Date())))
  )
  meta_html <- paste0(
    "<dl class=\"meta\">",
    paste0(
      vapply(meta_items, function(it) {
        paste0("<dt>", esc(it[1]), "</dt><dd>", it[2], "</dd>")
      }, character(1)),
      collapse = ""
    ),
    "</dl>"
  )

  # --- colocalising-loci table ----------------------------------------------

  if (nrow(results) > 0) {
    is_coloc <- !is.na(results$PP.H4) & results$PP.H4 >= pp4_thr
  } else {
    is_coloc <- logical(0)
  }

  if (nrow(results) > 0 && any(is_coloc)) {
    cl <- results[is_coloc, , drop = FALSE]
    ratio <- cl$PP.H4 / (cl$PP.H3 + cl$PP.H4)
    coloc_tbl_df <- data.frame(
      locus = cl$locus,
      nsnps = cl$nsnps,
      PP.H3 = round(cl$PP.H3, 3),
      PP.H4 = round(cl$PP.H4, 3),
      `PP.H4/(PP.H3+PP.H4)` = round(ratio, 3),
      lead_snp = cl$lead_snp,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    coloc_section <- html_table(coloc_tbl_df)
  } else {
    coloc_section <- sprintf(
      "<p>No loci colocalised at PP.H4 &gt;= %.2f.</p>", pp4_thr
    )
  }

  # --- all-tested-loci table -------------------------------------------------

  if (nrow(results) > 0) {
    all_df <- results
    # Round PP columns for display where present.
    pp_cols <- intersect(
      c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"), names(all_df)
    )
    for (cc in pp_cols) all_df[[cc]] <- round(all_df[[cc]], 3)
    row_class <- ifelse(is_coloc, "coloc", "")
    all_section <- html_table(all_df, row_class = row_class)
  } else {
    all_section <- "<p>No loci were tested.</p>"
  }

  # --- plots section ---------------------------------------------------------

  plots_html <- embed_png(miami_path, "Miami plot")
  if (!is.null(zoom_paths) && length(zoom_paths) > 0) {
    znames <- names(zoom_paths)
    for (k in seq_along(zoom_paths)) {
      zp <- zoom_paths[[k]]
      lid <- if (!is.null(znames) && nzchar(znames[k])) znames[k] else basename(zp)
      plots_html <- paste0(plots_html, embed_png(zp, paste0("Locus zoom — ", lid)))
    }
  }
  if (!nzchar(plots_html)) {
    plots_html <- "<p>No plots available.</p>"
  }

  # --- assemble --------------------------------------------------------------

  html <- paste0(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n",
    "<meta charset=\"utf-8\">\n",
    "<title>coloc-pair report</title>\n",
    style, "\n</head>\n<body>\n",
    "<h1>coloc-pair report</h1>\n",
    "<h2>Run inputs</h2>\n", meta_html, "\n",
    "<h2>Colocalising loci</h2>\n", coloc_section, "\n",
    "<h2>All tested loci</h2>\n", all_section, "\n",
    "<h2>Plots</h2>\n", plots_html, "\n",
    "</body>\n</html>\n"
  )

  writeLines(html, outfile)
  invisible(outfile)
}
