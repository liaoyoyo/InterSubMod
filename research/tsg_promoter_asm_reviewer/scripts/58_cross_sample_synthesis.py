#!/usr/bin/env python3
"""
58 - Deterministic cross-sample synthesis of HCC1395 key ASM positions (v2).

Reads the 6 per-sample JSONs from 57_cross_sample_keypos.py and classifies each of
the 38 key positions. NO LLM in the number path (anti-fabrication by construction,
CLAUDE.md §13 layer A).

v2 CHANGES (adversarial-eval NEEDS_WORK fixes, 2026-06-03):
  1. TWO orthogonal axes per position (the v1 single label conflated them):
       somatic_status  : recurrent_somatic / hcc1395_private_somatic / no_somatic_asm
       germline_status : recurrent_germline_asm / single / none
     A position can be BOTH hcc1395_private_somatic AND recurrent_germline_asm in
     others (e.g. SOX2) -- v1 hid this; v2 reports both.
  2. "imprinting" is RESERVED for SIGN-CONCORDANT recurrence only. Opposite-sign
     germline deltas across samples argue AGAINST shared parent-of-origin imprinting
     (= sporadic allelic methylation / noise). Each recurrent_germline_asm hit gets
     sign_concordant + imprinting_consistent flags. (eval found 3/4 v1 "imprinting"
     hits were sign-DISCORDANT.)
  3. SUBGROUP N surfaced: every delta carries the min N of the two groups it rests on
     (HP-axis: main & sub; germline: HP1 & HP2) + low_subgroup_n flag (<COV_MIN).

THRESHOLDS (a priori): ASM_MIN=0.05 (HP-axis), GERM_MIN=0.10 (germline), COV_MIN=10.

Output: genome_survey_v2/cn_confound/cross_sample/cross_sample_synthesis.json
"""
import os, json
import numpy as np

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
CS = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample"

REF_SAMPLE = "HCC1395"
OTHER = ["COLO829", "H1437", "H2009", "HCC1937", "HCC1954"]
ALL = [REF_SAMPLE] + OTHER
CANCER = {
    "HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
    "H1437": "lung", "H2009": "lung", "COLO829": "melanoma",
}

ASM_MIN = 0.05
GERM_MIN = 0.10
COV_MIN = 10


def load(sample):
    p = f"{CS}/{sample}_keypos.json"
    return json.load(open(p)) if os.path.exists(p) else None


def pk(r):
    return f"{r['chrom']}:{r['pos']}"


def grp_n(r, tag):
    return r.get("hp_groups", {}).get(tag, {}).get("n", 0)


def hp_axis_subgroup_minN(r):
    """min N of (main, sub) for the credible axis 'HPk_vs_HPk-1'."""
    axis = r.get("axis", "")
    if "_vs_" not in axis:
        return None
    main, sub = (t.replace("HP", "") for t in axis.split("_vs_"))
    nm, ns = grp_n(r, main), grp_n(r, sub)
    if nm == 0 or ns == 0:
        return None
    return min(nm, ns)


def germ_subgroup_minN(r):
    n1, n2 = grp_n(r, "1"), grp_n(r, "2")
    if n1 == 0 or n2 == 0:
        return None
    return min(n1, n2)


def somatic_asm_state(r):
    """True if this sample shows somatic-controlled ASM at the locus."""
    hp_d = r.get("hp_axis_delta")
    return bool(r.get("somatic_in_sample") and r.get("has_somatic_subhap")
                and hp_d is not None and abs(hp_d) >= ASM_MIN)


def germline_asm_state(r):
    """True if HP1-vs-HP2 |Δβ| >= GERM_MIN (allelic / candidate-imprinting)."""
    g = r.get("germline_delta")
    return g is not None and abs(g) >= GERM_MIN


def main():
    data = {s: load(s) for s in ALL}
    missing = [s for s in ALL if data[s] is None]
    idx = {s: {pk(r): r for r in data[s]["per_position"]} for s in ALL if data[s]}

    ref_records = data[REF_SAMPLE]["per_position"]
    positions = []
    somatic_counts = {}
    germline_counts = {}
    imprinting_consistent_hits = []
    discordant_recurrent_hits = []
    recurrent_somatic_hits = []

    for ref in ref_records:
        key = pk(ref)
        gene = ref["gene"]
        ref_hp_d = ref.get("hp_axis_delta")
        ref_hp_minN = hp_axis_subgroup_minN(ref)

        per_sample = {}
        for s in ALL:
            r = idx.get(s, {}).get(key)
            if r is None:
                per_sample[s] = dict(state="missing")
                continue
            som = somatic_asm_state(r)
            germ = germline_asm_state(r)
            per_sample[s] = dict(
                cancer=CANCER[s],
                somatic_asm=som, germline_asm=germ,
                somatic_in_sample=bool(r.get("somatic_in_sample")),
                has_subhap=bool(r.get("has_somatic_subhap")),
                hp_axis_delta=r.get("hp_axis_delta"),
                hp_axis_minN=hp_axis_subgroup_minN(r),
                germline_delta=r.get("germline_delta"),
                germline_minN=germ_subgroup_minN(r),
                n_reads=r.get("n_reads_window"),
            )

        # ---- somatic axis ----
        others_som = [s for s in OTHER if per_sample[s].get("somatic_asm")]
        ref_som = per_sample[REF_SAMPLE].get("somatic_asm", False)
        if others_som:
            somatic_status = "recurrent_somatic"
            recurrent_somatic_hits.append(dict(
                pos=key, gene=gene, samples=others_som,
                deltas={s: per_sample[s]["hp_axis_delta"] for s in others_som}))
        elif ref_som:
            somatic_status = "hcc1395_private_somatic"
        else:
            somatic_status = "no_somatic_asm"
        somatic_counts[somatic_status] = somatic_counts.get(somatic_status, 0) + 1

        # ---- germline axis (direction-aware) ----
        germ_samples = [s for s in ALL if per_sample[s].get("germline_asm")]
        germ_deltas = {s: per_sample[s]["germline_delta"] for s in germ_samples}
        germ_signs = set(np.sign(v) for v in germ_deltas.values() if v is not None)
        sign_concordant = len(germ_signs) == 1
        n_germ = len(germ_samples)
        imprinting_consistent = (n_germ >= 2 and sign_concordant)
        low_germ_n = any((per_sample[s]["germline_minN"] or 0) < COV_MIN
                         for s in germ_samples)
        if n_germ >= 2:
            germline_status = "recurrent_germline_asm"
            hit = dict(pos=key, gene=gene, n_samples=n_germ, samples=germ_samples,
                       deltas=germ_deltas,
                       cancer_types=sorted(set(CANCER[s] for s in germ_samples)),
                       sign_concordant=sign_concordant,
                       imprinting_consistent=imprinting_consistent,
                       low_subgroup_n=low_germ_n)
            (imprinting_consistent_hits if imprinting_consistent
             else discordant_recurrent_hits).append(hit)
        elif n_germ == 1:
            germline_status = "single_germline_asm"
        else:
            germline_status = "none"
        germline_counts[germline_status] = germline_counts.get(germline_status, 0) + 1

        positions.append(dict(
            pos=key, gene=gene, set=ref["set"], axis=ref["axis"],
            cn_class=ref.get("cn_class"),
            hcc1395_hp_axis_delta=ref_hp_d, hcc1395_hp_axis_minN=ref_hp_minN,
            hcc1395_low_subgroup_n=(ref_hp_minN is not None and ref_hp_minN < COV_MIN),
            somatic_status=somatic_status, germline_status=germline_status,
            n_other_somatic_asm=len(others_som),
            other_somatic_asm_samples=others_som,
            n_germline_asm_samples=n_germ, germline_asm_samples=germ_samples,
            germline_sign_concordant=sign_concordant if n_germ >= 2 else None,
            imprinting_consistent=imprinting_consistent,
            per_sample=per_sample,
        ))

    # per-sample somatic overlap
    somatic_overlap = {}
    for s in ALL:
        if not data[s]:
            continue
        recs = data[s]["per_position"]
        somatic_overlap[s] = dict(
            cancer=CANCER[s],
            n_somatic_at_keypos=sum(1 for r in recs if r.get("somatic_in_sample")),
            n_with_subhap=sum(1 for r in recs if r.get("has_somatic_subhap")),
            n_somatic_asm=sum(1 for r in recs if somatic_asm_state(r)),
            of_total=len(recs))

    out = dict(
        meta=dict(
            script="58_cross_sample_synthesis.py (v2)",
            ref_sample=REF_SAMPLE, other_samples=OTHER, cancer_types=CANCER,
            samples_available=[s for s in ALL if data[s]], samples_missing=missing,
            n_key_positions=len(positions),
            thresholds=dict(ASM_MIN=ASM_MIN, GERM_MIN=GERM_MIN, COV_MIN=COV_MIN),
            note=("Two orthogonal axes (somatic / germline). 'imprinting_consistent' "
                  "reserved for SIGN-CONCORDANT germline recurrence (opposite-sign "
                  "deltas => sporadic allelic methylation, NOT shared imprinting). "
                  "Numbers deterministic from 57 JSONs (no LLM). Cross-sample method = "
                  "window-mean 5mC beta (consistent across samples); reproduces "
                  "DIRECTION of HCC1395 credible-loci deltas, not exact magnitude "
                  "(table used paired-CpG MAX-collapse Wilcoxon). low_subgroup_n flags "
                  "deltas resting on a group with <COV_MIN reads."),
        ),
        headline=dict(
            somatic_status_counts=somatic_counts,
            germline_status_counts=germline_counts,
            n_recurrent_somatic=len(recurrent_somatic_hits),
            n_recurrent_germline_asm=len(imprinting_consistent_hits) + len(discordant_recurrent_hits),
            n_imprinting_consistent=len(imprinting_consistent_hits),
            n_recurrent_discordant=len(discordant_recurrent_hits),
            recurrent_somatic_hits=recurrent_somatic_hits,
            imprinting_consistent_hits=imprinting_consistent_hits,
            discordant_recurrent_hits=discordant_recurrent_hits,
        ),
        somatic_overlap=somatic_overlap,
        per_position=positions,
    )
    outp = f"{CS}/cross_sample_synthesis.json"
    with open(outp, "w") as f:
        json.dump(out, f, indent=2,
                  default=lambda o: None if isinstance(o, float) and np.isnan(o) else o)

    print(f"[58v2] wrote {outp}")
    print(f"[58v2] somatic_status: {somatic_counts}")
    print(f"[58v2] germline_status: {germline_counts}")
    print(f"[58v2] recurrent_somatic={len(recurrent_somatic_hits)} "
          f"recurrent_germline_asm={len(imprinting_consistent_hits)+len(discordant_recurrent_hits)} "
          f"(imprinting_consistent={len(imprinting_consistent_hits)} "
          f"discordant={len(discordant_recurrent_hits)})")
    print("[58v2] per-sample somatic overlap @ key positions:")
    for s, v in somatic_overlap.items():
        print(f"     {s:9s} ({v['cancer']:8s}): somatic={v['n_somatic_at_keypos']:2d}/"
              f"{v['of_total']}  subhap={v['n_with_subhap']:2d}  somatic_ASM={v['n_somatic_asm']:2d}")
    print("[58v2] IMPRINTING-CONSISTENT (sign-concordant) hits:")
    for h in imprinting_consistent_hits:
        print(f"     {h['pos']:18s} {h['gene'][:14]:14s} n={h['n_samples']} "
              f"deltas={h['deltas']} lowN={h['low_subgroup_n']}")
    print("[58v2] DISCORDANT recurrent (NOT imprinting) hits:")
    for h in discordant_recurrent_hits:
        print(f"     {h['pos']:18s} {h['gene'][:14]:14s} n={h['n_samples']} "
              f"deltas={h['deltas']} cancers={h['cancer_types']}")


if __name__ == "__main__":
    main()
