#!/usr/bin/env python3
"""Pilot: augmented VerificationClass (fix A+B+C+D) vs original — rescue & false-positive quantification.

Original VC = 2x2(label_sig × cluster_sig). Derive orig flags from the original VC column, then augment:
  A: cluster_sig OR HP-AUC>=0.7 (continuous distance structure, not just discrete cluster×label match)
  C: cluster_sig OR raw GlobalP<=0.05 (ignore CramersV-reliability gate)
  B: label_sig OR Δβ sig (calibrated mean-shift enters the verdict)
  D: structure (via A/C) without label match & no label signal -> new class "StructureNoLabel" (not Noise)

Reports: how many original Weak/Noise are rescued to a valid-structure class, and how many original
Strong change (should be ~0 — fixes only ADD recognition). No source edit; characterization only.
"""
import csv
import math
import sys
from collections import Counter

CSV = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dbeta_wg_abc/significance_summary.csv"
AUC_T = 0.7

rows = list(csv.DictReader(open(CSV)))
n = len(rows)


def f(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (ValueError, TypeError):
        return None


def tg(r, c):
    return r.get(c) == "true"


def new_vc(r):
    ovc = r["VerificationClass"]
    orig_cluster = ovc in ("Strong", "Subclone")
    orig_label = ovc in ("Strong", "Weak")
    # A: continuous HP structure
    aucn, auca = f(r.get("HP_AUC_Normal")), f(r.get("HP_AUC_All"))
    a_struct = (aucn is not None and aucn >= AUC_T) or (auca is not None and auca >= AUC_T)
    # C: raw association ignoring reliability gate
    gp = f(r.get("GlobalP"))
    c_raw = gp is not None and gp <= 0.05
    # B: calibrated Δβ signal
    b_dbeta = any(tg(r, c) for c in ["GermlineAsmDbeta_Sig", "SubcloneDbeta_HP1_Sig", "SubcloneDbeta_HP2_Sig",
                                     "SomaticResidualDbeta_Sig", "AltSubcloneDbeta_HP1_Sig", "AltSubcloneDbeta_HP2_Sig"])
    new_cluster = orig_cluster or a_struct or c_raw
    new_label = orig_label or b_dbeta
    if new_cluster and new_label:
        return "Strong"
    if new_cluster and not new_label:
        return "Subclone" if orig_cluster else "StructureNoLabel"  # D
    if not new_cluster and new_label:
        return "Weak"
    return "Noise"


orig = [r["VerificationClass"] for r in rows]
new = [new_vc(r) for r in rows]

valid = ("Strong", "Subclone", "StructureNoLabel")
rescued = sum(1 for i in range(n) if orig[i] in ("Weak", "Noise") and new[i] in valid)
weaknoise = sum(1 for i in range(n) if orig[i] in ("Weak", "Noise"))
strong_changed = sum(1 for i in range(n) if orig[i] == "Strong" and new[i] != "Strong")
strong_total = sum(1 for x in orig if x == "Strong")
# transitions
trans = Counter((orig[i], new[i]) for i in range(n) if orig[i] != new[i])

out = {
    "n": n,
    "orig_dist": dict(Counter(orig)),
    "new_dist": dict(Counter(new)),
    "rescued_weaknoise_to_valid": rescued,
    "weaknoise_total": weaknoise,
    "rescue_rate": round(rescued / weaknoise, 4) if weaknoise else None,
    "orig_Strong_changed (false-positive check)": strong_changed,
    "orig_Strong_total": strong_total,
    "top_transitions": {f"{k[0]}->{k[1]}": v for k, v in sorted(trans.items(), key=lambda kv: -kv[1])[:8]},
}
import json  # noqa: E402

print(json.dumps(out, indent=2, ensure_ascii=False))
