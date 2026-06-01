# coloc-pair

A command-line tool that runs **colocalisation** on two GWAS studies and tells
you, for each shared locus, whether they are driven by the **same causal
variant**. Point it at two GWAS Catalog accessions (or two local harmonised
files) and it produces:

- `coloc_results.tsv` — one row per locus with the coloc posteriors `PP.H0`–`PP.H4`
- `miami.png` — both traits back-to-back, colocalising loci highlighted
- `locuszoom_<locus>.png` — a regional plot for each of the top loci
- `report.html` — a single self-contained page with the tables and every plot

Colocalisation uses `coloc::coloc.abf` (one causal variant per locus). A locus
**colocalises** when `PP.H4 >= 0.8`.

---

## QuickStart: T2D vs LDL from scratch

This runs end to end on two real GWAS Catalog studies and needs nothing but R
and a network connection. It was tested from a clean R library — install,
run, open the report.

**1. Get the code and install dependencies** (one time; all from CRAN):

```bash
git clone git@github.com:malawsky/humgen_workshop_demo.git
cd humgen_workshop_demo
Rscript install.R
```

**2. Run the colocalisation** — Type 2 diabetes vs LDL cholesterol, both from
the Genes & Health British-Bangladeshi/Pakistani cohort:

```bash
Rscript R/coloc_pair.R \
  --study1 GCST90727286 --type1 cc --s1 0.2578 \
  --study2 GCST90727331 --type2 quant \
  --sig-mode either \
  --outdir results/t2d_ldl
```

- `--study1 GCST90727286` — T2D, a **case/control** trait (`--type1 cc`).
  `--s1 0.2578` is the case fraction (11,348 cases / 44,026 total).
- `--study2 GCST90727331` — LDL cholesterol, a **quantitative** trait.
- **Sample sizes are pulled automatically** from GWAS Catalog metadata
  (N = 44,026 and 25,080), so no `--n` is needed for accessions.
- The first run downloads the two harmonised files (~48 MB + ~36 MB) into
  `results/t2d_ldl/cache/` and reuses them on later runs.

It takes roughly **a minute** after the download. You should see:

```
Loci tested: 20
Colocalising (PP.H4 >= 0.80): 0
Results table: results/t2d_ldl/coloc_results.tsv
Miami plot: results/t2d_ldl/miami.png
Locus-zoom plots written: 5
HTML report: results/t2d_ldl/report.html
```

**3. Open the report:**

```bash
open results/t2d_ldl/report.html        # macOS  (Linux: xdg-open)
```

**Reading the result.** Across the 20 genome-wide-significant loci, none reach
the `PP.H4 >= 0.8` colocalisation threshold — T2D and LDL are largely driven by
*different* causal variants here, which is the expected biology. The strongest
shared-signal evidence is at **chr19:18.8–19.8 Mb** (`PP.H4 ≈ 0.40`, near
*TM6SF2*); the well-known *APOE* locus (`rs7412`, chr19:44–45 Mb) shows a clear
LDL-only signal (`PP.H2 ≈ 0.85`). To see what a *positive* colocalisation looks
like, run the worked example below.

---

## 30-second smoke test (no network)

A synthetic pair ships in `data/example/` with a planted answer: one shared
locus, one distinct. Use it to confirm your install works.

```bash
Rscript R/coloc_pair.R \
  --study1 data/example/trait1.h.tsv.gz --type1 quant \
  --study2 data/example/trait2.h.tsv.gz --type2 quant \
  --n1 50000 --n2 50000 \
  --outdir results/example
```

Locus A colocalises (`PP.H4 ≈ 0.93`); locus B does not (`PP.H3 ≈ 1.0`). Local
files carry no metadata, so `--n1`/`--n2` are required for them.

To check the trusted statistical core on its own: `Rscript tests/test_worked_example.R`.

---

## Running your own studies

Either argument can be a **GWAS Catalog accession** (downloaded and cached) or a
**local harmonised `*.h.tsv.gz` file** (read directly). Both the older `hm_*`
and the modern GWAS-SSF harmonised column layouts are supported.

Key options (`Rscript R/coloc_pair.R --help` for all):

| Flag | Meaning | Default |
|---|---|---|
| `--study1`, `--study2` | accession or local file (required) | — |
| `--type1`, `--type2` | `quant` or `cc` (case/control) | `quant` |
| `--s1`, `--s2` | case fraction; **required for `cc`** | — |
| `--n1`, `--n2` | sample size; **required for local files** (auto for accessions) | metadata |
| `--sig-mode` | define loci significant in `either`, `both`, `1`, or `2` | `either` |
| `--window` | locus half-width in bp (merged where they overlap) | `5e5` |
| `--p-threshold` | genome-wide significance | `5e-8` |
| `--pp4-threshold` | `PP.H4` to call a colocalisation | `0.8` |
| `--top-n` | how many top loci get a locus-zoom plot | `5` |
| `--outdir` | where results, plots, cache are written | `results/run` |

Notes:

- Pick two studies that **have harmonised summary statistics** on the GWAS
  Catalog. Full files are large; the first run caches them under `<outdir>/cache/`.
  The first run also caches the *parsed* sumstats as `<study>.std.v1.fst`, so a
  re-run skips both the download and the multi-minute parse of full-genome files.
- `cache/` and `results/` are git-ignored.
- Requirements: R (tested on 4.5) and the CRAN packages in `install.R`
  (`coloc`, `data.table`, `R.utils`, `gwasrapidd`, `curl`, `ggplot2`, `knitr`,
  `fst`, `optparse`). No Bioconductor or system tools needed.

See `SPEC.md` for the full specification and `CLAUDE.md` for conventions.
