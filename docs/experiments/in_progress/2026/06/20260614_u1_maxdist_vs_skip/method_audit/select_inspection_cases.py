#!/usr/bin/env python3
"""Select 4-layer + edge-case regions from chr1 significance_summary for manual visual inspection.

Layers (task-gate 見樹也見林): canonical / subclone-candidate / somatic-residual / edge / outlier.
Output: JSON list of selected regions with Pos + key metrics + region dir name. No fabrication —
every value read from the summary CSV this run.
"""
import csv
import json
import math
import sys

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dbeta_chr1_s2/significance_summary.csv"


def fnum(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (ValueError, TypeError):
        return None


def boolc(x):
    return x == "true"


rows = []
with open(path) as f:
    for r in csv.DictReader(f):
        rows.append(r)


def g(r, k):
    return fnum(r.get(k))


def pick(r):
    """Compact metric dict for a selected region."""
    return {
        "pos": r.get("Pos"),
        "chr": r.get("Chr"),
        "num_reads": int(g(r, "NumReads") or 0),
        "nreads_valid": int(g(r, "NReadsValid") or 0),
        "hp_auc_normal": g(r, "HP_AUC_Normal"),
        "hp_auc_tumor": g(r, "HP_AUC_Tumor"),
        "global_p": g(r, "GlobalP"),
        "cramers_v": g(r, "CramersV"),
        "permanova_f": g(r, "ClusterPermanovaF"),
        "permanova_p": g(r, "ClusterPermanovaP"),
        "permanova_valid": r.get("ClusterPermanovaValid"),
        "tumor_hp1": int(g(r, "Tumor_HP1") or 0),
        "tumor_hp2": int(g(r, "Tumor_HP2") or 0),
        "normal_hp1": int(g(r, "Normal_HP1") or 0),
        "normal_hp2": int(g(r, "Normal_HP2") or 0),
        "germline_asm_dbeta": g(r, "GermlineAsmDbeta"),
        "germline_asm_sig": r.get("GermlineAsmDbeta_Sig"),
        "subclone_hp1": g(r, "SubcloneDbeta_HP1"),
        "subclone_hp1_sig": r.get("SubcloneDbeta_HP1_Sig"),
        "subclone_hp2": g(r, "SubcloneDbeta_HP2"),
        "subclone_hp2_sig": r.get("SubcloneDbeta_HP2_Sig"),
        "somatic_residual_dbeta": g(r, "SomaticResidualDbeta"),
        "somatic_residual_sig": r.get("SomaticResidualDbeta_Sig"),
        "verification_class": r.get("VerificationClass"),
        "significant": r.get("Significant"),
    }


sel = {}

# 1. CANONICAL: strong germline ASM, everything agrees, decent reads
cand = [
    r for r in rows
    if boolc(r.get("GermlineAsmDbeta_Sig", "")) and (g(r, "NumReads") or 0) >= 40
    and (g(r, "CramersV") or 0) >= 0.4 and abs(g(r, "GermlineAsmDbeta") or 0) >= 0.15
    and (g(r, "HP_AUC_Normal") or 0) >= 0.75
]
cand.sort(key=lambda r: -abs(g(r, "GermlineAsmDbeta") or 0))
sel["canonical_germline_asm"] = [pick(r) for r in cand[:2]]

# 2. SUBCLONE candidate: same-hap subclone Δβ sig + strong, both germline+carrier present, decent reads
cand = [
    r for r in rows
    if (boolc(r.get("SubcloneDbeta_HP1_Sig", "")) and abs(g(r, "SubcloneDbeta_HP1") or 0) >= 0.15)
    or (boolc(r.get("SubcloneDbeta_HP2_Sig", "")) and abs(g(r, "SubcloneDbeta_HP2") or 0) >= 0.15)
]
cand = [r for r in cand if (g(r, "NumReads") or 0) >= 40]
cand.sort(key=lambda r: -max(abs(g(r, "SubcloneDbeta_HP1") or 0), abs(g(r, "SubcloneDbeta_HP2") or 0)))
sel["subclone_candidate"] = [pick(r) for r in cand[:2]]

# 3. SOMATIC residual sig: tumor differs from normal baseline, strong
cand = [
    r for r in rows
    if boolc(r.get("SomaticResidualDbeta_Sig", "")) and abs(g(r, "SomaticResidualDbeta") or 0) >= 0.15
    and (g(r, "NumReads") or 0) >= 40
]
cand.sort(key=lambda r: -abs(g(r, "SomaticResidualDbeta") or 0))
sel["somatic_residual_sig"] = [pick(r) for r in cand[:2]]

# 4a. EDGE tiny-N: smallest NumReads that still produced clustering output (valid permanova)
cand = [r for r in rows if (g(r, "NumReads") or 0) >= 6 and r.get("ClusterPermanovaValid") == "true"]
cand.sort(key=lambda r: (g(r, "NumReads") or 0))
sel["edge_tiny_n"] = [pick(r) for r in cand[:2]]

# 4b. EDGE single-HP tumor (somatic haplotag dominance: one tumor HP empty but other populated)
cand = [
    r for r in rows
    if ((g(r, "Tumor_HP1") or 0) == 0) != ((g(r, "Tumor_HP2") or 0) == 0)  # exactly one is zero
    and ((g(r, "Tumor_HP1") or 0) + (g(r, "Tumor_HP2") or 0)) >= 10
]
cand.sort(key=lambda r: -((g(r, "Tumor_HP1") or 0) + (g(r, "Tumor_HP2") or 0)))
sel["edge_single_hp_tumor"] = [pick(r) for r in cand[:2]]

# 4c. EDGE NaN-heavy: largest drop NumReads -> NReadsValid (SKIP removed many)
cand = [r for r in rows if (g(r, "NumReads") or 0) >= 30 and (g(r, "NReadsValid") or 0) >= 6]
cand.sort(key=lambda r: -((g(r, "NumReads") or 0) - (g(r, "NReadsValid") or 0)))
sel["edge_nan_heavy"] = [pick(r) for r in cand[:2]]

# 5. OUTLIER: CramersV == 1.0 (perfect cluster-label match) with FEW reads (overfit suspicion)
cand = [
    r for r in rows
    if (g(r, "CramersV") or 0) >= 0.999 and (g(r, "NumReads") or 0) <= 25 and (g(r, "NumReads") or 0) >= 6
]
cand.sort(key=lambda r: (g(r, "NumReads") or 0))
sel["outlier_cramersv1_lowN"] = [pick(r) for r in cand[:2]]

# summary counts
counts = {k: len(v) for k, v in sel.items()}
out = {"source": path, "n_total_regions": len(rows), "selected_counts": counts, "cases": sel}
print(json.dumps(out, indent=2, ensure_ascii=False))
