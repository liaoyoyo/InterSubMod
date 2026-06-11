#!/usr/bin/env python3
"""
98 - Build the per-locus methylation-clustering CATALOG skeleton (T-CODE-1, ISM-4).

One locus = one row. Tags each of the 332,705 ISM-scanned loci (6 samples x TP/FP/FN)
with the 16-column schema (9 classification + 7 numeric) defined in
docs/paper_focus/02_paper_framework/位點甲基分群catalog_schema提案.md (定稿 2026-06-09),
then collapses each locus to ONE of 7 TAGs (A clean-cis / B copy-artifact / C germline-
allelic / D high-dbeta FP-prone / E latent / F untestable / G no-signal).

GUARDRAIL (characterization != filter): the catalog is a DESCRIPTIVE directory.
TAG-D (FP-prone) is annotation only, NEVER a filter (methylation-filter is concluded
NEGATIVE, 4-ways dead). minN<15 latent (TAG-E) is FP-leaning -> describe, do not filter.

Provenance of every column value:
  - clustering_reliability/cramersV/permanova/dbeta/axis/LOH/coverage/n_CpG : genome-wide,
    direct from ism_existence_scan/<sample>_<cls>/significance_summary.csv (each loci a row).
  - reliable gate = `Significant` column (= PassedGating & GlobalP<=0.05 & CramersV>=0.1 &
    NumReads>=20), exactly the gate in script 88; latent = script-88 PERMANOVA recovery (minN=5).
  - cis_status/within_d/copy_d/cis_perm_p : HCC1395 only, from cis_scan_full.json (816 HP-axis
    loci) + copy_partition_confirm.json (strongest loci d_within vs d_copy). Reconciled with the
    2026-06-07 capstone (chr17=clean cis, BRCA2=copy/subclone, 6 strongest=untestable pure-ALT).
  - gene : curated known loci only (no genome-wide gene model bundled); else NA.

Outputs (genome_survey_v2/catalog/):
  catalog_skeleton.tsv      - one row per locus, all schema columns
  catalog_tag_counts.json   - per-TAG counts overall + per sample x class
  catalog_audit.json        - T-CODE-2 feed: genome-wide reliable/latent/none + latent TP/FP minN sweep
"""
import os, csv, json, glob
import numpy as np
from collections import Counter, defaultdict

EX = "/big7_disk/liaoyoyo2001/ism_existence_scan"
ROOT = "/big7_disk/liaoyoyo2001/InterSubMod/research/tsg_promoter_asm_reviewer"
GS = f"{ROOT}/genome_survey_v2"
OUTDIR = f"{GS}/catalog"
os.makedirs(OUTDIR, exist_ok=True)

SAMPLES = ["HCC1395", "HCC1937", "HCC1954", "H1437", "H2009", "COLO829"]
CANCER = {"HCC1395": "breast", "HCC1937": "breast", "HCC1954": "breast",
          "H1437": "lung", "H2009": "lung", "COLO829": "melanoma"}
CLS = ["tp", "fp", "fn"]

LATENT_MINN = 5            # script-88 default floor (avoids 2-3 read degeneracy; minN<15 FP-leaning)
DBETA_LARGE, DBETA_MED = 0.20, 0.10   # schema TBD cutoffs (large>=.2 / med .1-.2 / small<.1)
LOWCOV = 20               # low coverage for TAG-D regression-to-extreme regime

# curated gene annotation (no genome-wide gene model bundled)
GENE = {
    "chr17:79991120": ("TBC1D16", False),
    "chr13:32315128": ("ZAR1L/BRCA2", True),
    "chr13:32339132": ("BRCA2", True),
    "chr13:32337160": ("BRCA2", True),
}
# 2026-06-07 capstone within-control verdicts (verified): loci flagged for a documented
# mechanical/CpG-creation caveat -> cannot be claimed as clean cis even if within-dominant.
MECH_CAVEAT = {"chr18:11741161"}      # CREATES a CpG -> methyl diff may be sequence artifact
CIS_MARGIN = 1.5                       # cis requires within-allele d clearly > copy d (not a tie)


def num(r, k):
    try:
        v = float(r.get(k, ""))
        return None if v != v else v
    except (ValueError, TypeError):
        return None


def tb(r, k):
    return str(r.get(k, "")).lower() == "true"


def is_latent(r, minn):
    """script-88 canonical latent = real cluster structure the CramersV Cochran gate zeroed."""
    if tb(r, "Significant"):
        return False
    if not tb(r, "PassedGating"):
        return False
    if (num(r, "NumReads") or 0) < 20:
        return False
    if not (tb(r, "ClusterPermanovaValid") and (num(r, "ClusterPermanovaP") or 1) <= 0.05):
        return False
    if not (tb(r, "LabelHPPermanovaValid") and (num(r, "LabelHPPermanovaP") or 1) <= 0.05):
        return False
    mn = min(num(r, "HP1FamilyN") or 0, num(r, "HP2FamilyN") or 0)
    return mn >= minn


def reliability(r):
    if tb(r, "Significant"):
        return "reliable"
    if is_latent(r, LATENT_MINN):
        return "latent"
    return "none"


def dbeta_max(r):
    cand = [abs(num(r, k) or 0) for k in ("AlleleDelta", "HPMergedDelta")]
    return round(max(cand), 4) if cand else 0.0


def dbeta_mag(v):
    return "large" if v >= DBETA_LARGE else ("medium" if v >= DBETA_MED else "small")


def axis_of(r):
    d = str(r.get("DominantLabel", "")).lower()
    return {"hp": "HP-axis", "allele": "ALLELE-axis"}.get(d, "none")


# ---- HCC1395 cis / copy join sources ----
cis = {d["locus"]: d for d in json.load(open(f"{GS}/cis_scan_full.json"))}
cp = {d["locus"]: d for d in json.load(open(f"{GS}/copy_partition_confirm.json"))["rows"]}
BONF = 0.05 / max(len(cis), 1)

# per-locus blind-ARI (high-ARI subset only) for clustering_origin
ari = {}
amap = f"{GS}/q2_high_ari_perlocus_map.tsv"
if os.path.exists(amap):
    for row in csv.DictReader(open(amap), delimiter="\t"):
        try:
            ari[row["locus"]] = float(row["ari_blind_allele"])
        except (ValueError, KeyError):
            pass


def cis_status(locus):
    """capstone-reconciled cis judgement. Only HCC1395 HP-axis loci are cis-tested.
    NOTE: p_cis can be 0.0 (strongest) -> MUST use `is not None`, never truthiness (0.0 is falsy)."""
    d = cis.get(locus)
    if d is None:
        return "NA", None, None, None         # not cis-tested
    pcis = d.get("p_cis")
    psig = (pcis is not None and pcis <= 0.05)
    c = cp.get(locus)
    if c is not None:                          # has within-subclone control (leaky tag)
        dw, dc = abs(c["d_within"]), abs(c["d_copy"])
        if dw > CIS_MARGIN * dc and psig:      # within-allele CLEARLY drives > copy (not a tie)
            if locus in MECH_CAVEAT:
                return "cis-mechanical", c["d_within"], c["d_copy"], pcis  # CpG-creation caveat
            return "cis", c["d_within"], c["d_copy"], pcis
        return "copy", c["d_within"], c["d_copy"], pcis      # copy/subclone dominant or marginal
    # no within-subclone control (pure-ALT clonal) -> cannot separate cis vs copy
    if d["cis_tier"] == "T3":
        return "untestable", None, None, pcis  # T3 cis-candidate but needs within-control / COLO829
    return "not-cis", None, None, pcis


def clustering_origin(locus, rel, ax, cstat):
    if cstat == "cis":
        return "somatic-linked"
    if locus in ari:
        # blind-ARI ruler: high ARI vs germline-het background -> germline-allelic (B2 default)
        return "germline-allelic"
    if rel in ("reliable", "latent"):
        return "germline-allelic"            # B2: read-clustering is germline-allelic background
    return "drift"


def assign_tag(rel, cstat, dmag, loh, cov):
    """collapse to ONE of 7 TAGs. cis-based tags (A/B/F) are HCC1395-only & most specific."""
    if cstat == "cis":
        return "A"                            # clean cis somatic-ASM (chr17/TBC1D16)
    if cstat == "copy":
        return "B"                            # copy/subclone-confounded (BRCA2)
    if cstat in ("untestable", "cis-mechanical"):
        return "F"                            # pure-ALT clonal OR CpG-creation caveat -> needs new method / COLO829
    if rel == "latent":
        return "E"                            # CramersV-gated but PERMANOVA-real clustering
    if rel == "reliable":
        return "C"                            # reliable clustering = germline-allelic background
    # rel == none
    if dmag == "large" and loh and cov < LOWCOV:
        return "D"                            # high-dbeta FP-prone regression-to-extreme (ANNOTATE ONLY)
    return "G"                                # no signal


def main():
    cols = ["locus_id", "chrom", "pos", "ref", "alt", "sample", "cancer", "class",
            "somatic_status", "clustering_reliability", "cramersV_value", "permanova_p",
            "passed_gating", "dbeta_max", "dbeta_magnitude", "axis", "clustering_origin",
            "ari_value", "cis_status", "within_d", "copy_d", "cis_perm_p",
            "mod_type", "recurrence", "LOH_status", "loh_subtype", "coverage", "n_CpG",
            "gene", "is_TSG_promoter", "provenance_tier", "TAG"]
    fout = open(f"{OUTDIR}/catalog_skeleton.tsv", "w", newline="")
    w = csv.writer(fout, delimiter="\t", lineterminator="\n")  # LF (avoid CRLF tripping awk $NF)
    w.writerow(cols)

    tag_overall = Counter()
    tag_by_sc = defaultdict(Counter)          # (sample,class) -> TAG counts
    rel_by_sc = defaultdict(Counter)          # (sample,class) -> reliability counts
    # T-CODE-2 audit: latent minN sweep, genome-wide TP/FP
    minn_sweep = {m: {"tp": 0, "fp": 0} for m in (0, 3, 5, 8, 10, 15)}
    reliable_tpfp = {"tp": 0, "fp": 0}
    n_rows = 0
    examples = defaultdict(list)              # TAG -> up to 5 canonical loci (HCC1395 preferred)

    for s in SAMPLES:
        for cl in CLS:
            f = f"{EX}/{s}_{cl}/significance_summary.csv"
            if not os.path.exists(f):
                continue
            for r in csv.DictReader(open(f)):
                chrom, pos = r.get("Chr"), r.get("Pos")
                if not chrom or not pos:
                    continue
                locus = f"{chrom}:{pos}"
                rel = reliability(r)
                cv = num(r, "CramersV")
                pp = num(r, "ClusterPermanovaP")
                dmaxv = dbeta_max(r)
                dmag = dbeta_mag(dmaxv)
                ax = axis_of(r)
                cov = int(num(r, "NumReads") or 0)
                ncpg = int(num(r, "NumCpGs") or 0)
                loh = tb(r, "Potential_LOH")
                loh_sub = r.get("LOH_Subtype", "None")
                # cis only for HCC1395 (cis_scan is HCC1395 HP-axis)
                if s == "HCC1395":
                    cstat, wd, cd, pc = cis_status(locus)
                else:
                    cstat, wd, cd, pc = "NA", None, None, None
                origin = clustering_origin(locus, rel, ax, cstat)
                tag = assign_tag(rel, cstat, dmag, loh, cov)
                som = {"tp": "somatic", "fn": "somatic-missed", "fp": "FP"}[cl]
                gene, tsg = GENE.get(locus, ("NA", False))
                ptier = "🟢" if s == "HCC1395" and locus in cis else "🟢"  # all skeleton cols grep-able
                w.writerow([locus, chrom, pos, r.get("Ref"), r.get("Alt"), s, CANCER[s], cl.upper(),
                            som, rel, cv, pp, tb(r, "PassedGating"), dmaxv, dmag, ax, origin,
                            ari.get(locus, ""), cstat, wd, cd, pc,
                            "NA", "single", ("LOH" if loh else "nonLOH"), loh_sub, cov, ncpg,
                            gene, tsg, ptier, tag])
                tag_overall[tag] += 1
                tag_by_sc[(s, cl)][tag] += 1
                rel_by_sc[(s, cl)][rel] += 1
                if rel == "reliable" and cl in ("tp", "fp"):
                    reliable_tpfp[cl] += 1
                if cl in ("tp", "fp"):
                    for m in minn_sweep:
                        if (not tb(r, "Significant")) and is_latent(r, m):
                            minn_sweep[m][cl] += 1
                if len(examples[tag]) < 5 and (s == "HCC1395" or len(examples[tag]) < 2):
                    examples[tag].append({"locus": locus, "sample": s, "class": cl.upper(),
                                          "gene": gene, "cramersV": cv, "dbeta_max": dmaxv,
                                          "cis_status": cstat, "ari": ari.get(locus)})
                n_rows += 1
    fout.close()

    # TAG counts
    tag_counts = {
        "n_loci": n_rows,
        "overall": dict(tag_overall),
        "tag_legend": {"A": "clean cis somatic-ASM", "B": "copy/subclone-confounded",
                       "C": "germline-allelic background", "D": "high-dbeta FP-prone (ANNOTATE ONLY)",
                       "E": "latent (PERMANOVA-recovered)", "F": "untestable (pure-ALT clonal)",
                       "G": "no signal"},
        "by_sample_class": {f"{s}_{cl}": dict(tag_by_sc[(s, cl)]) for s in SAMPLES for cl in CLS
                            if (s, cl) in tag_by_sc},
        "examples": {k: v for k, v in examples.items()},
    }
    json.dump(tag_counts, open(f"{OUTDIR}/catalog_tag_counts.json", "w"), indent=1, ensure_ascii=False)

    # T-CODE-2 audit feed
    audit = {
        "reliability_definition": {
            "reliable": "Significant = PassedGating & GlobalP<=0.05 & CramersV>=0.1 & NumReads>=20 (CramersV Cochran-gated)",
            "latent": f"NOT Significant & PassedGating & NumReads>=20 & ClusterPermanovaValid&P<=0.05 & LabelHPPermanovaValid&P<=0.05 & min(HP1FamilyN,HP2FamilyN)>={LATENT_MINN}",
            "none": "neither",
        },
        "reliability_by_sample_class": {f"{s}_{cl}": dict(rel_by_sc[(s, cl)])
                                        for s in SAMPLES for cl in CLS if (s, cl) in rel_by_sc},
        "reliable_TP_FP_genomewide": reliable_tpfp,
        "latent_minN_sweep_TP_FP_genomewide": {str(m): minn_sweep[m] for m in minn_sweep},
        "interpretation": ("CramersV Cochran gate zeroes sparse-table CramersV; latent recovers "
                           "real HP clustering via PERMANOVA. minN sweep: rec TP:FP stays ~base 47:1 "
                           "(enrich ~1) => characterization gain, NOT a TP/FP filter. Recommend 3-state "
                           "reliable/latent/unreliable, used for description only."),
    }
    json.dump(audit, open(f"{OUTDIR}/catalog_audit.json", "w"), indent=1, ensure_ascii=False)

    # console summary (for Read-back verification)
    print(f"[98] catalog rows = {n_rows}")
    print(f"[98] TAG overall = {dict(tag_overall)}")
    print(f"[98] reliable TP/FP genome-wide = {reliable_tpfp}")
    print(f"[98] latent minN sweep (TP/FP) = {{m: minn_sweep[m] for m in minn_sweep}}".replace("minn_sweep[m]", "") or "")
    for m in (0, 5, 15):
        print(f"      minN={m}: TP={minn_sweep[m]['tp']} FP={minn_sweep[m]['fp']}")
    print(f"[98] chr17 TAG = {[e for e in examples.get('A', []) if e['locus']=='chr17:79991120']}")
    print(f"[98] wrote {OUTDIR}/catalog_skeleton.tsv + catalog_tag_counts.json + catalog_audit.json")
    print("DONE")


if __name__ == "__main__":
    main()
