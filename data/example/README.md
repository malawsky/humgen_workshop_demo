# Worked example: a planted answer

These two files are **synthetic** GWAS summary statistics in GWAS Catalog
*harmonised* column format. They exist so you can check the tool against an
answer you already know, with no network access.

| file              | trait                         |
|-------------------|-------------------------------|
| `trait1.h.tsv.gz` | simulated trait 1 (quant)     |
| `trait2.h.tsv.gz` | simulated trait 2 (quant)     |

Both files share an identical variant panel on a fake `chr1`, with two
well-separated regions:

- **Locus A (`chr1:~50.0–50.4 Mb`)** — both traits have their peak at the
  **same** causal variant. Colocalisation should return a **high PP.H4**
  (one shared causal variant).
- **Locus B (`chr1:~60.0–60.4 Mb`)** — the two traits peak at **different**
  variants ~200 kb apart. Colocalisation should return a **high PP.H3**
  (two distinct causal variants).

The expected outcome is encoded in `../../tests/expected_truth.tsv` and checked
by `tests/test_worked_example.R`. With the trusted wrapper the planted answer
reproduces as roughly **PP.H4 ≈ 0.93** at locus A and **PP.H3 ≈ 1.0** at locus B.

Regenerate with `python3 make_example.py` or `Rscript make_example.R`
(equivalent). If you regenerate, re-check `expected_truth.tsv`.
