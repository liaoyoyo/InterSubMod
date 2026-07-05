#!/usr/bin/env python3
"""Fill the missing n_residual_candidate_cnclean using the new SAVANA CN.

The parallel session's methyl_auxiliary_annotation.py finds, per sSNV-genotype group,
methylation bimodality that survives HP exclusion = 'residual-candidate' (L3 hidden
epigenetic substructure candidate). But it could not compute cn-clean (cn=unknown) so
n_residual_candidate_cnclean=0 for every sample.

Here we apply the SAVANA SM_CNBED (usable-purity samples) to each residual-candidate
region -> real CN state -> count how many are cn-clean (cn != gain, matching methyl_aux's
cn_flag=(cn=='gain')). cn-clean residual = methylation substructure NOT explained by
HP and NOT in a CN-amplified (multiplicity) region = the strongest 'beyond genetic' signal.
Read-only.
"""
import json
import os
from collections import Counter, defaultdict

MSR = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
CNV = "/big7_disk/liaoyoyo2001/cnv_sv_work"
SAMPLES = ["H1437", "H2009", "HCC1954"]  # usable SAVANA purity


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
    return "neutral"  # SM_CN_SPARSE=0 convention


def parse_region(rg):
    chrom, rng = rg.split(":")
    s, e = rng.split("-")
    return chrom, int(s), int(e)


def run(sample):
    aux = json.load(open(f"{MSR}/{sample}/methyl_auxiliary_annotation.json"))
    iv = load_cn(f"{CNV}/{sample}/savana_wgs/{sample}_smcnbed.bed")
    resid = [h for h in aux["hidden_candidates"] if h["flag"] == "residual-candidate"]
    cn_dist = Counter()
    clean = []
    for h in resid:
        chrom, s, e = parse_region(h["region"])
        cn = cn_at(iv, chrom, (s + e) // 2)
        cn_dist[cn] += 1
        h["_savana_cn"] = cn
        if cn != "gain":  # methyl_aux cn_flag = (cn=='gain'); cn-clean = not gain
            clean.append(h)
    # delta_beta of the two methylation modes (means)
    def db(h):
        m = h.get("means") or []
        return round(abs(m[0] - m[1]), 3) if len(m) == 2 else None
    clean_db = [db(h) for h in clean if db(h) is not None]
    # strict clean = neutral or loh only (exclude loss too, the safest)
    strict = [h for h in clean if h["_savana_cn"] in ("neutral", "loh")]
    return {
        "sample": sample,
        "n_tested": aux["summary"]["b_hidden_substructure"]["n_groups_tested"],
        "n_residual_candidate": len(resid),
        "cn_dist_of_residual": dict(cn_dist),
        "n_residual_cnclean_notgain": len(clean),
        "n_residual_strict_neutral_or_loh": len(strict),
        "clean_median_delta_beta": (sorted(clean_db)[len(clean_db) // 2] if clean_db else None),
        "clean_examples": [{"region": h["region"], "geno": h["genotype"], "cn": h["_savana_cn"],
                            "means": h["means"], "sizes": h["sizes"], "n_reads": h["n_reads"]}
                           for h in sorted(clean, key=lambda x: -(abs(x["means"][0] - x["means"][1]) if len(x.get("means", [])) == 2 else 0))[:5]],
    }


def main():
    res = [run(s) for s in SAMPLES]
    outp = "/big7_disk/liaoyoyo2001/cnv_sv_work/methyl_lineage/residual_cnclean_savana.json"
    os.makedirs(os.path.dirname(outp), exist_ok=True)
    json.dump(res, open(outp, "w"), ensure_ascii=False, indent=2)
    print(f"{'sample':9} {'tested':>7} {'resid':>6} {'cn_dist(of residual)':40} {'cnclean(≠gain)':>14} {'strict(neut/loh)':>16}")
    for r in res:
        print(f"{r['sample']:9} {r['n_tested']:>7} {r['n_residual_candidate']:>6} "
              f"{str(r['cn_dist_of_residual']):40} {r['n_residual_cnclean_notgain']:>14} "
              f"{r['n_residual_strict_neutral_or_loh']:>16}")
    print(f"\nwrote {outp}")
    for r in res:
        print(f"\n{r['sample']} clean(≠gain) median Δβ={r['clean_median_delta_beta']}; top examples:")
        for e in r["clean_examples"]:
            print(f"  {e['region']} geno={e['geno']} cn={e['cn']} means={e['means']} sizes={e['sizes']}")


if __name__ == "__main__":
    main()
