# SPEC — coloc-pair

A command-line tool that colocalises two GWAS Catalog studies.

## What it does

Given two studies (each a GWAS Catalog accession **or** a local harmonised
sumstats file) and an output directory, it:

1. **Loads** harmonised summary statistics for both traits and standardises the
   columns (`snp, chr, pos, ea, oa, maf, beta, se, p`). Accessions are
   downloaded and cached under `<outdir>/cache/`; local files are read directly.
2. **Defines significant loci** — variants with `p < --p-threshold` (default
   5e-8), windowed by `--window` (default 500 kb) and merged where they overlap.
   `--sig-mode` chooses whether a locus must be significant in trait 1, trait 2,
   either (default), or both.
3. **Runs colocalisation** per locus with `coloc::coloc.abf` via the trusted
   wrapper in `R/coloc_wrapper.R`, producing PP.H0–H4 and the lead shared SNP.
4. **Writes a results table** `<outdir>/coloc_results.tsv`, one row per locus,
   sorted by PP.H4 descending.
5. **Draws a Miami plot** `<outdir>/miami.png` — trait 1 up, trait 2 down,
   shared x-axis, with colocalising loci (PP.H4 ≥ `--pp4-threshold`) highlighted.
6. **Draws locus-zoom plots** `<outdir>/locuszoom_<locus>.png` for the top
   `--top-n` loci by PP.H4.
7. Prints a short summary: loci tested, loci that colocalise, output paths.

## Interface

Arguments are already defined in `R/coloc_pair.R` (`--study1/2`, `--type1/2`,
`--s1/2`, `--n1/2`, `--sig-mode`, `--window`, `--p-threshold`,
`--pp4-threshold`, `--top-n`, `--outdir`, `--keep-sumstats`). Implement the
pipeline that those arguments drive; do not change the flag names.

## Build order (suggested, one commit each)

1. `load_sumstats()` in `download_sumstats.R` — make the **local-file** path
   work first so the worked example runs, then add accession download/caching.
2. `define_loci()` in `define_loci.R`.
3. Wire steps 1–4 of the pipeline in `coloc_pair.R`; produce the results table.
4. Add an **integration test** (`tests/test_pipeline.R`) that runs the CLI on
   `data/example/` and checks PP.H4 high at locus A, plus that the table was
   written.
5. `make_miami()` in `plots.R`, then wire step 5.
6. `make_locuszoom()` in `plots.R`, then wire step 6.

## Done when

- `Rscript tests/test_worked_example.R` is green (it already is).
- `Rscript R/coloc_pair.R --study1 data/example/trait1.h.tsv.gz --type1 quant
  --study2 data/example/trait2.h.tsv.gz --type2 quant --outdir results/example`
  writes a results table, a Miami plot, and at least one locus-zoom, and reports
  locus A as colocalising.
- The reviewer subagent has looked over the diff and its points are addressed.

## Deliberately out of scope

LD-based clumping/conditioning, multi-causal coloc (SuSiE), cross-build
liftover, and heavy input validation. Single-causal `coloc.abf` is enough here.
