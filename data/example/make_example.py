#!/usr/bin/env python3
"""
Generate a synthetic pair of GWAS summary-statistics files in GWAS Catalog
*harmonised* format, with a deliberately PLANTED answer.

The point of this file is verification, not realism. Two small traits are
simulated over two well-separated regions on a fake chromosome:

  Locus A (~chr1:50.0-50.4 Mb):  BOTH traits peak at the SAME variant.
                                 -> coloc should return a high PP.H4 (shared).

  Locus B (~chr1:60.0-60.4 Mb):  the two traits peak at DIFFERENT variants,
                                 ~200 kb apart.
                                 -> coloc should return a high PP.H3 (distinct).

Both files share an identical variant panel (same SNP ids, positions, alleles
and allele frequencies), so coloc aligns them cleanly. The signal is encoded
as a smooth Gaussian "LD bump" of z-scores around each causal position.

The committed *.h.tsv.gz files are the source of truth for the worked-example
test. Re-run this script (or make_example.R) only if you want to regenerate
them; if you do, update tests/expected_truth.tsv to match.
"""
import gzip
import numpy as np

RNG = np.random.default_rng(20260601)
SE = 0.020                      # constant standard error per SNP
ALLELES = np.array(["A", "C", "G", "T"])


def region(start, n, spacing=1000):
    """Evenly spaced positions for a region, returned as an int array."""
    return start + np.arange(n) * spacing


def bump(pos, causal_pos, z_peak, width=40_000):
    """A smooth Gaussian peak of z-scores centred on causal_pos."""
    z = z_peak * np.exp(-((pos - causal_pos) ** 2) / (2.0 * width ** 2))
    z += RNG.normal(0, 0.25, size=pos.shape)        # mild noise
    return z


def build_panel():
    """Shared variant panel across both traits (two regions)."""
    posA = region(50_000_000, 400)
    posB = region(60_000_000, 400)
    pos = np.concatenate([posA, posB])
    n = pos.size

    # Random but fixed alleles and allele frequencies, shared by both traits.
    ea = RNG.choice(ALLELES, size=n)
    oa = RNG.choice(ALLELES, size=n)
    same = ea == oa                                 # avoid ea == oa
    oa[same] = ALLELES[(np.searchsorted(ALLELES, ea[same]) + 1) % 4]
    eaf = RNG.uniform(0.10, 0.45, size=n)
    rsid = np.array([f"rs{9_000_000 + i}" for i in range(n)])
    chrom = np.ones(n, dtype=int)
    return dict(rsid=rsid, chrom=chrom, pos=pos, ea=ea, oa=oa, eaf=eaf,
                posA=posA, posB=posB)


def make_trait(panel, causalA, causalB):
    """z-scores -> beta/se/p for one trait given its causal positions."""
    pos = panel["pos"]
    z = np.zeros(pos.size)
    inA = pos < 55_000_000
    inB = ~inA
    z[inA] = bump(pos[inA], causalA, z_peak=9.5)    # lead p ~ 1e-21
    z[inB] = bump(pos[inB], causalB, z_peak=9.0)
    beta = z * SE
    se = np.full(pos.size, SE)
    # two-sided p from the normal; clip to avoid exact zero
    from math import erfc, sqrt
    p = np.array([erfc(abs(zi) / sqrt(2.0)) for zi in z])
    p = np.clip(p, 1e-300, 1.0)
    return beta, se, p


def write_harmonised(path, panel, beta, se, p):
    cols = ["hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele",
            "hm_effect_allele", "hm_effect_allele_frequency",
            "hm_beta", "standard_error", "p_value"]
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for i in range(panel["pos"].size):
            fh.write("\t".join([
                panel["rsid"][i],
                str(panel["chrom"][i]),
                str(int(panel["pos"][i])),
                panel["oa"][i],
                panel["ea"][i],
                f"{panel['eaf'][i]:.4f}",
                f"{beta[i]:.6f}",
                f"{se[i]:.6f}",
                f"{p[i]:.3e}",
            ]) + "\n")


def main():
    panel = build_panel()

    # Locus A: SHARED causal variant (same position in both traits).
    sharedA = panel["posA"][200]            # centre of region A
    # Locus B: DISTINCT causal variants, ~200 kb apart.
    b1 = panel["posB"][150]
    b2 = panel["posB"][350]

    beta1, se1, p1 = make_trait(panel, causalA=sharedA, causalB=b1)
    beta2, se2, p2 = make_trait(panel, causalA=sharedA, causalB=b2)

    write_harmonised("trait1.h.tsv.gz", panel, beta1, se1, p1)
    write_harmonised("trait2.h.tsv.gz", panel, beta2, se2, p2)

    # Emit the planted truth for the test to check against.
    with open("../../tests/expected_truth.tsv", "w") as fh:
        fh.write("locus\tregion_start\tregion_end\texpected\tlead_pos\n")
        fh.write(f"A\t49990000\t50410000\tH4\t{int(sharedA)}\n")
        fh.write(f"B\t59990000\t60410000\tH3\t{int(b1)},{int(b2)}\n")

    print("wrote trait1.h.tsv.gz, trait2.h.tsv.gz and tests/expected_truth.tsv")
    print(f"  Locus A shared causal at chr1:{int(sharedA)}  -> expect PP.H4 high")
    print(f"  Locus B distinct causals at chr1:{int(b1)} / {int(b2)} -> expect PP.H3 high")


if __name__ == "__main__":
    main()
