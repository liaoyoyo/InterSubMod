#!/usr/bin/env python3
"""
87 - Rebuild manifest.json correctly (no re-render needed; the PNGs are fine).

Bug fixed: script 85 took CramersV from clustering/significance.json, which holds the
RAW (pre-reliability) cramers_v. The summary's CramersV (what the Significant gate uses)
is `cramers_v_reliable ? cramers_v : 0` where reliable = min expected cell count >= 5
(Cochran). They DIFFER for sparse tables. We now keep BOTH, clearly separated:

  cv_* (gated)  : from existence-scan summary via curated_loci.json -> drives PASS/MISSED
  raw_* (raw)   : from significance.json -> latent structure the reliability gate hid
  perm_*        : PERMANOVA (distance-based, robust to sparsity) -> is the structure real?

Derived `latent` flag marks MISSED loci whose CramersV was gated to ~0 but PERMANOVA
still says the reads cluster (valid & p<=0.05 & raw effect present) = prime "did we miss
a real one?" candidates for the user to judge.

Output: display_v2/manifest.json  (overwrites 85's, same figs/)
"""
import json, os, glob

ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
DV = f"{ROOT}/genome_survey_v2/cn_confound/cross_sample/display_v2"
SCAN = "/big7_disk/liaoyoyo2001/ism_display_scan"
FIGDIR = f"{DV}/figs"


def region_dir(cls, chrom, pos):
    for d in glob.glob(f"{SCAN}/HCC1395_{cls}/curated_{cls}/{chrom}/{chrom}_{pos}/*"):
        if os.path.exists(f"{d}/clustering/significance.json"):
            return d
    return None


def jstats(rd):
    s = json.load(open(f"{rd}/clustering/significance.json"))
    g = lambda k: s.get(k, {})
    raw = {a: round(g("global_" + a).get("cramers_v", 0.0), 3)
           for a in ("alt", "hp", "hp_family", "hp_fine")}
    cs = g("cluster_structure"); ls = g("label_structure")
    return dict(
        raw_alt=raw["alt"], raw_hp=raw["hp"], raw_hpfam=raw["hp_family"], raw_hpfine=raw["hp_fine"],
        raw_max=round(max(raw.values()), 3), optimal_k=s.get("optimal_k"),
        perm_f=round(cs.get("permanova_f", 0), 1), perm_p=cs.get("permanova_p"),
        perm_valid=bool(cs.get("permanova_valid", False)),
        disp_warn=bool(cs.get("dispersion_warning", False)),
        hp_perm_f=round(ls.get("hp_permanova_f", 0), 1), hp_perm_p=ls.get("hp_permanova_p"),
        hp_perm_valid=bool(ls.get("hp_permanova_valid", False)),
    )


def main():
    curated = json.load(open(f"{DV}/curated_loci.json"))
    out = []
    n_latent = 0
    for d in curated:
        chrom, pos, cls = d["chr"], d["pos"], d["cls"]
        key = f"{chrom}_{pos}"
        rec = dict(d)  # gated cv_* + db + gp + reads + ncpg + loh + passed_gating + reject(base) + category
        rec["has_fig"] = os.path.exists(f"{FIGDIR}/{key}_meth.png")
        rd = region_dir(cls, chrom, pos)
        if rd:
            st = jstats(rd)
            rec.update(st)
            # latent real structure that the CramersV reliability gate hid
            latent = (d["category"] == "MISSED" and d["cv_max"] < 0.1 and st["raw_max"] >= 0.1
                      and st["perm_valid"] and (st["perm_p"] is not None and st["perm_p"] <= 0.05))
            rec["latent"] = bool(latent)
            if latent:
                n_latent += 1
                rec["reject"] = (f"CramersV 卡方不可靠(表稀疏 expected<5)→歸零；"
                                 f"但 PERMANOVA F={st['perm_f']} p={st['perm_p']:.2g} 有結構（原始 CramersV_max={st['raw_max']:.2f}）— 可能漏掉")
        else:
            rec["latent"] = False
            rec.setdefault("raw_max", 0.0)
        out.append(rec)

    with open(f"{DV}/manifest.json", "w") as f:
        json.dump(out, f)
    fig = sum(1 for x in out if x["has_fig"])
    print(f"[87] manifest rebuilt: {len(out)} loci, has_fig={fig}, latent(被可靠性閘控的真結構)={n_latent}")
    # show a few latent examples
    lat = [x for x in out if x.get("latent")]
    lat.sort(key=lambda x: -x.get("perm_f", 0))
    for x in lat[:6]:
        print(f"   LATENT {x['chr']}:{x['pos']} gatedCV={x['cv_max']:.2f} rawCV={x['raw_max']:.2f} "
              f"permF={x['perm_f']} Δβ={x['db']:+.2f} reads={x['reads']}")


if __name__ == "__main__":
    main()
