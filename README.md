# coloc-pair

A small command-line tool that takes **two GWAS Catalog studies**, runs
**colocalisation** on their significant loci, and produces a **Miami plot**, a
**results table**, and **locus-zoom plots** for the top hits.

This repository is a **template for a live build**: the structure, the trusted
statistical core, the worked example, and the agent roles are in place; the
download/loci/plot/pipeline code is left as the build target (`SPEC.md`).

## Why it is laid out this way

It follows the workshop principles:

- **Context as layout** — files sit where the agent expects (`R/`, `data/`,
  `tests/`, `results/`); `CLAUDE.md` is the standing brief.
- **A trusted template** — `R/coloc_wrapper.R` is known-good and must not be
  rewritten; the agent builds around it.
- **An answer you already know** — `data/example/` has a synthetic pair with a
  *planted* result (one shared locus, one distinct), checked by the test.
- **Specialised agents** — planner, coder, reviewer, tester in `.claude/agents/`.
- **Version control as memory** — commit after each working step.

## Setup

```bash
Rscript install.R        # coloc, topr, locuszoomr, gwasrapidd, gwascatftp, ...
```

## The worked example (no network needed)

```bash
Rscript tests/test_worked_example.R
```

This runs the trusted wrapper on `data/example/` and checks the planted answer:
locus A colocalises (PP.H4 ≈ 0.93), locus B does not (PP.H3 ≈ 1.0). It passes
*before* any of the build target is written — it is checking the trusted core.

## The build target

See `SPEC.md`. In short, make this work end to end:

```bash
Rscript R/coloc_pair.R \
  --study1 data/example/trait1.h.tsv.gz --type1 quant \
  --study2 data/example/trait2.h.tsv.gz --type2 quant \
  --outdir results/example
```

and then on real studies by passing GWAS Catalog accessions instead of paths:

```bash
Rscript R/coloc_pair.R \
  --study1 GCST90012345 --type1 cc --s1 0.2 \
  --study2 GCST90067890 --type2 quant \
  --sig-mode either --outdir results/run1
```

## The loop to run during the build

1. Write one clear, scoped task (a step from `SPEC.md`).
2. Have the **planner** plan it; read the plan.
3. Approve it, or correct it and plan again.
4. Let the **coder** implement that one step.
5. Run your checks (`tests/`, and the **tester** for new features).
6. **Review** the diff (`/review`, the **reviewer** subagent), then commit.
7. Clear context, start the next step.

## Notes for real data

- Harmonised sumstats are large; the first run downloads and caches them under
  `<outdir>/cache/`. **Pre-warm the cache before a live demo** so the run is
  quick. Cache and outputs are git-ignored.
- Pick two studies that have *harmonised* summary statistics on the GWAS Catalog
  and that plausibly share a locus (e.g. two related lipid or glycaemic traits)
  so there is a real PP.H4 to show. Verify availability beforehand.
- `--type` is `cc` (case/control) or `quant`; case/control traits also need
  `--s` (case proportion).
