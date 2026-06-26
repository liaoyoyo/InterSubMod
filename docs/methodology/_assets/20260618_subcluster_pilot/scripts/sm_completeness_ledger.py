#!/usr/bin/env python3
"""[S4 完整性帳本 — no-miss 保證]

讀 sm_linkage_genomewide.json (census + universe)。每個 union sSNV 歸入恰好一個 bucket：
- linked            : n_realized_links>=1 (有 >=1 powered 同讀 partner, coread>=COREAD_MIN)
- underpowered      : 有 <=50kb partner 但無 powered link (覆蓋/共讀不足)
- isolated_singleton: <=50kb 內無任何 union sSNV (single-molecule 不可連)
sum-check: 三 bucket 加總 == universe 總 sSNV 數。
另量「union 修補的漏失」: FP-source 且 somatic-confirmed 且 linked = 若只用 TP 會漏掉的 clonally-informative
(γ 類)；以及 phase_gap 候選 (isolated 但可能有 same-PS 遠 partner, 待 Tier-PS 補)。
輸出 sm_completeness_ledger.json。數字落檔。
"""
import json
from collections import Counter

A = "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot"
IN = f"{A}/sm_linkage_genomewide.json"
OUT = f"{A}/sm_completeness_ledger.json"


def main():
    d = json.load(open(IN))
    census = d["census"]
    uni = d["universe"]

    bucket = {}
    by_src = Counter()
    by_bucket = Counter()
    fp_somatic_linked = []     # union 修補的漏失 (TP-only 會漏)
    linked_somatic = Counter()
    for key, c in census.items():
        np_ = c["n_partners_le50k"]
        nr = c["n_realized_links"]
        if np_ == 0:
            b = "isolated_singleton"
        elif nr >= 1:
            b = "linked"
        else:
            b = "underpowered"
        bucket[key] = b
        by_bucket[b] += 1
        by_src[(c["src"], b)] += 1
        if b == "linked" and c.get("somatic"):
            linked_somatic[c["src"]] += 1
            if c["src"] == "FP":
                fp_somatic_linked.append({"locus": key, "vaf": c.get("vaf"),
                                          "normal": [c.get("normal_ref"), c.get("normal_alt")],
                                          "n_links": nr, "max_coread": c.get("max_coread")})

    total = len(census)
    s = sum(by_bucket.values())
    assert s == total, f"sum-check FAIL {s} != {total}"

    out = {
        "universe_total": total,
        "universe_tp": uni["n_tp"], "universe_fp": uni["n_fp"],
        "buckets": dict(by_bucket),
        "bucket_sum_check": {"sum": s, "total": total, "ok": s == total},
        "by_source_bucket": {f"{k[0]}|{k[1]}": v for k, v in by_src.items()},
        "linked_somatic_by_source": dict(linked_somatic),
        "union_recovered_fp_somatic_linked": {
            "count": len(fp_somatic_linked),
            "note": "FP-source somatic-confirmed & single-molecule-linked = TP-only 會漏的 clonally-informative (γ 類)",
            "examples": sorted(fp_somatic_linked, key=lambda x: -(x["max_coread"] or 0))[:30],
        },
        "phase_gap_candidates_note": "isolated_singleton 中可能有 same-PS 遠 partner → 待 sm_linkage_phaseset.py(Tier-PS) 細分",
    }
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print(f"LEDGER total={total} buckets={dict(by_bucket)} sum_ok={s==total} "
          f"fp_somatic_linked={len(fp_somatic_linked)} linked_somatic={dict(linked_somatic)} -> {OUT}")


if __name__ == "__main__":
    main()
