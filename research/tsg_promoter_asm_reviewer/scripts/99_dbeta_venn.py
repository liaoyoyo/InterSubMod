#!/usr/bin/env python3
"""
99 - ISM-clustering vs Δβ-magnitude discovery-channel VENN (T-CODE-3, ISM-3).

ISM gates on discrete read-clustering (CramersV/PERMANOVA); Δβ is a per-position allele
mean difference. They are orthogonal discovery channels. This quantifies the union:
  ISM channel  = clustering detected (reliable OR latent)
  Δβ channel   = dbeta_max >= LARGE (0.20)
3 cells: ISM-only / Δβ-only / both, each with TP/FP composition + TP:FP ratio.

GUARDRAIL: this is NOT a naive Δβ add (Δβ-only is FP-enriched = D5 regression-to-extreme).
The honest union discovery = ISM-clustered OR (large-Δβ AND has-clustering). The Venn proves
why: Δβ-only (no clustering) carries the FP enrichment; the intersection is the clean signal.

Reads landed catalog (script 98). Output: catalog/discovery_venn.tsv + discovery_venn.json
"""
import csv, json
from collections import Counter

GS = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer/genome_survey_v2"
CAT = f"{GS}/catalog/catalog_skeleton.tsv"
LARGE = 0.20


def ratio(tp, fp):
    return round(tp / fp, 1) if fp else float("inf")


def main():
    rows = list(csv.DictReader(open(CAT), delimiter="\t"))
    cells = {"ISM_only": Counter(), "dbeta_only": Counter(), "both": Counter(), "neither": Counter()}
    # restrict to TP/FP (FN has no FP comparator); genome-wide
    for r in rows:
        cl = r["class"]
        if cl not in ("TP", "FP"):
            continue
        ism = r["clustering_reliability"] in ("reliable", "latent")
        dbig = float(r["dbeta_max"] or 0) >= LARGE
        cell = ("both" if ism and dbig else "ISM_only" if ism else "dbeta_only" if dbig else "neither")
        cells[cell][cl] += 1

    base = Counter()
    for r in rows:
        if r["class"] in ("TP", "FP"):
            base[r["class"]] += 1
    base_ratio = ratio(base["TP"], base["FP"])

    out = {"base_TP_FP": dict(base), "base_ratio": base_ratio, "large_dbeta_cut": LARGE, "cells": {}}
    print(f"base TP:FP = {base['TP']}:{base['FP']} = {base_ratio}:1\n")
    print(f"{'cell':>12} {'TP':>8} {'FP':>7} {'TP:FP':>9} {'vs base':>9}")
    for k in ("both", "ISM_only", "dbeta_only", "neither"):
        c = cells[k]
        rt = ratio(c["TP"], c["FP"])
        enr = round(rt / base_ratio, 2) if base_ratio and rt != float("inf") else None
        out["cells"][k] = {"TP": c["TP"], "FP": c["FP"], "TP_FP_ratio": rt, "enrich_vs_base": enr}
        print(f"{k:>12} {c['TP']:>8} {c['FP']:>7} {str(rt):>9} {str(enr):>9}")

    out["interpretation"] = (
        "dbeta_only (large Δβ, NO clustering) is the FP-enriched cell (enrich<1 vs base = "
        "MORE FP per TP) = D5 regression-to-extreme. ISM_only catches clustering-only loci "
        "(e.g. BRCA2: low per-position Δβ, real clustering). The honest union discovery channel "
        "= ISM-clustered OR (large-Δβ AND clustering), i.e. excludes the dbeta_only FP cell. "
        "Δβ does NOT replace CramersV; it is a complementary channel that must be clustering-gated."
    )

    with open(f"{GS}/catalog/discovery_venn.json", "w") as f:
        json.dump(out, f, indent=1)
    with open(f"{GS}/catalog/discovery_venn.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["cell", "TP", "FP", "TP_FP_ratio", "enrich_vs_base"])
        for k in ("both", "ISM_only", "dbeta_only", "neither"):
            c = out["cells"][k]
            w.writerow([k, c["TP"], c["FP"], c["TP_FP_ratio"], c["enrich_vs_base"]])
    print(f"\n[99] wrote {GS}/catalog/discovery_venn.{{tsv,json}}")
    print("DONE")


if __name__ == "__main__":
    main()
