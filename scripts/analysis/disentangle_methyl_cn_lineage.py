#!/usr/bin/env python3
"""Disentangle CN-gain vs sSNV-lineage contributions to within-genotype-group methyl bimodality.

methyl_auxiliary flags every bimodal genotype-group as hp-explained / cn-flagged / residual,
but cn came from topology (=unknown) so cn-flagged was near-empty. We re-classify ALL bimodal
groups using the SAVANA SM_CNBED cn (usable-purity samples):

  hp-explained (unchanged: bimodality along germline HP)
  -> of the REST (non-HP bimodal):
       my_cn == gain      -> CN-GAIN effect (amplification multiplicity)  [Q3 proof]
       my_cn == neutral   -> beyond-CN residual, PUREST sSNV-lineage/epi signal [Q2 proof]
       my_cn == loh/loss  -> beyond-gain but LOH has its own confound

This isolates each axis (Q1/Q4). Δβ (|means diff|) reported per class.
NOTE: this is the BIMODAL-level split (numerator). Proper RATE (bimodality rate | CN state,
needing unimodal denominators) requires the CN-aware methyl_auxiliary re-run (separate).
Read-only.
"""
import json
import os
from collections import defaultdict, Counter

MSR = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
CNV = "/big7_disk/liaoyoyo2001/cnv_sv_work"
SAMPLES = ["H1437", "H2009", "HCC1954"]


def load_cn(bed):
    iv = defaultdict(list)
    for ln in open(bed):
        p = ln.split()
        if len(p) < 4 or not p[1].isdigit():
            continue
        iv[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for c in iv:
        iv[c].sort()
    return iv


def cn_at(iv, chrom, pos):
    for s, e, lab in iv.get(chrom, []):
        if s <= pos <= e:
            return lab
    return "neutral"


def parse_region(rg):
    chrom, rng = rg.split(":")
    s, e = rng.split("-")
    return chrom, (int(s) + int(e)) // 2


def run(sample):
    aux = json.load(open(f"{MSR}/{sample}/methyl_auxiliary_annotation.json"))
    iv = load_cn(f"{CNV}/{sample}/savana_wgs/{sample}_smcnbed.bed")
    bim = aux["hidden_candidates"]  # ALL bimodal groups
    n_hp = 0
    cls = Counter()   # my-CN-based class of non-HP bimodal
    db_by = defaultdict(list)
    has_alt = Counter()  # non-HP bimodal with ALT in genotype, by class
    for h in bim:
        db = abs(h["means"][0] - h["means"][1]) if len(h.get("means", [])) == 2 else None
        if h["flag"] == "hp-explained":
            n_hp += 1
            continue
        chrom, mid = parse_region(h["region"])
        cn = cn_at(iv, chrom, mid)
        cls_key = "cn_gain" if cn == "gain" else ("residual_neutral" if cn == "neutral"
                  else ("residual_loh" if cn == "loh" else "residual_loss"))
        cls[cls_key] += 1
        if db is not None:
            db_by[cls_key].append(db)
        if "A" in h["genotype"]:
            has_alt[cls_key] += 1
    med = lambda x: round(sorted(x)[len(x) // 2], 3) if x else None
    return {
        "sample": sample,
        "n_bimodal_total": len(bim),
        "n_hp_explained": n_hp,
        "nonHP_by_myCN": dict(cls),
        "nonHP_median_delta_beta": {k: med(v) for k, v in db_by.items()},
        "nonHP_with_ALT_genotype": dict(has_alt),
    }


def main():
    res = [run(s) for s in SAMPLES]
    json.dump(res, open(f"{CNV}/methyl_lineage/disentangle_methyl_cn_lineage.json", "w"),
              ensure_ascii=False, indent=2)
    print(f"{'sample':9} {'bimodal':>7} {'hp_exp':>6} | non-HP by my SAVANA CN:")
    for r in res:
        c = r["nonHP_by_myCN"]
        print(f"{r['sample']:9} {r['n_bimodal_total']:>7} {r['n_hp_explained']:>6} | "
              f"CN-gain={c.get('cn_gain',0)} neutral={c.get('residual_neutral',0)} "
              f"loh={c.get('residual_loh',0)} loss={c.get('residual_loss',0)}")
        print(f"          Δβ median: gain={r['nonHP_median_delta_beta'].get('cn_gain')} "
              f"neutral={r['nonHP_median_delta_beta'].get('residual_neutral')} "
              f"loh={r['nonHP_median_delta_beta'].get('residual_loh')}")
        print(f"          帶 ALT genotype 的: gain={r['nonHP_with_ALT_genotype'].get('cn_gain',0)} "
              f"neutral={r['nonHP_with_ALT_genotype'].get('residual_neutral',0)} "
              f"loh={r['nonHP_with_ALT_genotype'].get('residual_loh',0)}")
    print(f"\nwrote {CNV}/methyl_lineage/disentangle_methyl_cn_lineage.json")


if __name__ == "__main__":
    main()
