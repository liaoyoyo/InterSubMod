"""
scan_cis_candidates.py — P3b systematic cis-candidate scan.

Runs the full pipeline per locus (extract -> collapse -> score -> genomic-context ->
cis-ladder -> mechanical-cis), ranks, and writes cis-candidate records. Answers
"how many BRCA2-like (T3) cis-candidates exist". Parallel via multiprocessing (多資源).

characterization only — NO TP/FP discriminator, NO F1 filter.

    python3 stages/scan_cis_candidates.py --candidates <json> --out <json> [--limit N --workers W]
"""
import os
import sys
import json
import argparse
from collections import Counter
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from stages.stage_extract import extract_or_cache  # noqa: E402
from lib.cis_asm_core import (load_level1_plus, collapse_modtype, cis_test, cohesion_per_tag,  # noqa: E402
                              dbeta_axis, dbeta_mod, classify_cis_tier, power_class)
from lib.genomic_context import GenomicContext  # noqa: E402
from lib.causation import RefContext  # noqa: E402

TUMOR = ("/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/"
         "20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam")
NORMAL = "/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam"
CACHE = os.path.join(HERE, "..", "cache", "level1")
LOH_BED = os.path.join(HERE, "..", "..", "output", "seqc2_loh_only.bed")

_gc = None
_rc = None


def _init():
    global _gc, _rc
    _gc = GenomicContext(loh_bed=LOH_BED)
    _rc = RefContext()


def score_one(locus):
    chrom, pos, ref, alt = locus["chrom"], locus["pos"], locus["ref"], locus["alt"]
    try:
        path, _ = extract_or_cache(chrom, pos, ref, alt, TUMOR, NORMAL, cache_dir=CACHE)
        Dp = load_level1_plus(path)
        D = collapse_modtype(Dp)
        s = str(pos)
        ctx = _gc.annotate(chrom, pos)
        ct = cis_test(D, s) if ctx["hp_axis_valid"] else {"testable": False, "n_shared_cpg": 0}
        coh = cohesion_per_tag(D, s)
        d_hp = dbeta_axis(D, s, {"1"}, {"1-1"})
        d_som = ct.get("d_somatic") if ct.get("testable") else d_hp
        tier = classify_cis_tier(d_som, ct.get("p_cis"), ct, ctx["hp_axis_valid"])
        return {
            "locus": f"{chrom}:{pos}", "ref": ref, "alt": alt,
            "total_cn": ctx["total_cn"], "loh": ctx["loh_status"], "primary_axis": ctx["primary_axis"],
            "dbeta_HP": d_hp, "dbeta_5mC": dbeta_mod(Dp, s, "m"), "dbeta_5hmC": dbeta_mod(Dp, s, "h"),
            "sil_HP11": coh.get("1-1"),
            "d_cis": ct.get("d_cis"), "d_drift": ct.get("d_drift"), "p_cis": ct.get("p_cis"),
            "n_shared_cpg": ct.get("n_shared_cpg", 0),
            "power_class": power_class(ct.get("n_shared_cpg", 0), ct.get("p_cis")),
            "cis_tier": tier["cis_tier"], "tier_label": tier["label"],
            "mechanical_cis": _rc.mechanical_cis(chrom, pos, ref, alt),
        }
    except Exception as e:
        return {"locus": f"{chrom}:{pos}", "error": str(e)[:120]}


def _rank_key(r):
    torder = {"T3": 0, "T2": 1, "T1": 2, "T0": 3}.get(r.get("cis_tier"), 4)
    return (torder, -(abs(r.get("d_cis") or 0)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    loci = json.load(open(a.candidates))
    if a.limit:
        loci = loci[:a.limit]
    with Pool(a.workers, initializer=_init) as pool:
        recs = pool.map(score_one, loci)
    errs = [r for r in recs if "error" in r]
    recs = [r for r in recs if "error" not in r]
    recs.sort(key=_rank_key)
    json.dump(recs, open(a.out, "w"), indent=1)
    tc = Counter(r.get("cis_tier") for r in recs)
    print(f"scanned {len(recs)} loci ({len(errs)} errors) -> {a.out}")
    print(f"cis_tier distribution: {dict(tc)}")
    t3 = [r for r in recs if r.get("cis_tier") == "T3"]
    print(f"\nT3 cis-candidates: {len(t3)}")
    for r in t3[:20]:
        print(f"  {r['locus']:18s} d_cis={r.get('d_cis')} d_drift={r.get('d_drift')} "
              f"cn={r['total_cn']} mech={r['mechanical_cis']} 5mC={r.get('dbeta_5mC')} 5hmC={r.get('dbeta_5hmC')}")


if __name__ == "__main__":
    main()
