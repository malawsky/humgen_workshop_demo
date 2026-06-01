# CLAUDE.md

This project is a small command-line tool, **coloc-pair**, that takes two GWAS
Catalog studies, runs colocalisation on their significant loci, and produces a
Miami plot, a results table, and locus-zoom plots. The build target is
described in `SPEC.md`. The audience is research scientists, so favour clear,
readable code over defensive engineering.

## Workflow rules
- Work in small steps. After each working step, `git add -A && git commit`.
- Read `SPEC.md` and make a plan before writing code. Use the **planner**
  subagent to plan, the **coder** to implement, the **reviewer** to check, and
  the **tester** to run checks. Roles are in `.claude/agents/`.
- The analysis is in **R**. The CLI entry point is `R/coloc_pair.R`.

## Trusted template — do not modify
- `R/coloc_wrapper.R` is a **trusted template**: its colocalisation logic has
  been checked by hand against the worked example. Build around it; call its
  functions; **do not rewrite it**. If you believe it is wrong, stop and ask.

## Facts and conventions
- "Significant" / genome-wide significant means **p < 5e-8**.
- "Colocalises" means **PP.H4 >= 0.8** (the `--pp4-threshold` flag).
- A locus is significant variants windowed by `--window` (default 500 kb) and
  merged where they overlap.
- Harmonised sumstats live on EBI:
  `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/<RANGE>/<GCST>/harmonised/<PMID>-<GCST>-<EFO>.h.tsv.gz` (GRCh38).
  Resolve the path with `gwascatftp`; get study metadata with `gwasrapidd`.
- Harmonised → standardised columns (fixed mapping):
  `hm_rsid→snp, hm_chrom→chr, hm_pos→pos, hm_effect_allele→ea,
  hm_other_allele→oa, hm_effect_allele_frequency→maf, hm_beta→beta,
  standard_error→se, p_value→p`.
- `coloc::coloc.abf` wants `varbeta = se^2`; case/control traits also need `s`.

## Data and outputs
- `data/example/` holds a synthetic worked example with a **planted answer**;
  the tests check against it. Treat it as read-only.
- Real downloaded sumstats are **temporary**: cache under the run's `--outdir`,
  never commit them (`.gitignore` already excludes `cache/` and `results/`).
- The real datasets are large; cache them once and reuse.

## Checks
- `Rscript tests/test_worked_example.R` must stay green. It runs the trusted
  wrapper on the worked example and asserts the planted PP.H4 / PP.H3 outcome.
- When you add the end-to-end pipeline, add an integration test that runs the
  CLI on `data/example/` and checks the same outcome plus that the expected
  output files were written.
