#!/usr/bin/env python3
"""LCA 增益的 A/B 對照（擴充能力）。

LongLineage 的 tag_bam 有 pre/post-LCA 兩套 receipt。本模組不預設它們是
乾淨 A/B：會逐染色體核對 input identity、threads 與非 LCA stats，通過才能
回答「LCA 到底買到什麼」。HCC1395 現有歷史 receipt 雖顯示 lv 覆蓋
86,037 → 427,519 = 4.97×，但有 in_bam / threads 不一致，因此只能當描述數，
不能當純 LCA causal effect。

先前一份查證誤判「無來源檔可 grep、不能搬」，因為它用
`lv_written − lca_resolved` 去湊 86,037（85,700，差 337）。那不是那個數的定義；
86,037 是 pre-LCA 那套 receipt 的 lv_written 加總，可逐檔 grep。
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from capability import (MALFORMED, Capability, Registry, file_ref, finalize, make,
                        probe, sample_in_path)

LCA_FIELDS = {"lv_written", "lca_resolved", "lca_candidates_sum"}
# threads 雖非生物輸入，仍是重現性身分的一部分。下列任一欄位不同，
# 就不把 pre/post 解讀為只差 LCA 開關的 causal A/B。
A_B_IDENTITY_FIELDS = {
    "assignments", "in_bam", "lineage_paths", "region", "schema_name",
    "schema_version", "sidecar", "threads", "topology",
}


def _agg(files):
    tot = {}
    for f in files:
        try:
            st = json.loads(Path(f).read_text(encoding="utf-8")).get("stats") or {}
        except (OSError, ValueError):
            continue
        for k, v in st.items():
            try:
                tot[k] = tot.get(k, 0) + int(v)
            except (TypeError, ValueError):
                pass
    return tot


def _read_doc(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _identity_error(path: Path, doc: dict, expected_sample: str) -> str:
    """檢查 receipt 檔名與上游 input identity。

    只核對 receipt 實際提供的身分欄位；至少 topology / in_bam /
    assignments 之一必須存在，否則不能單靠目錄名猜樣本。
    """
    if not sample_in_path(path.name, expected_sample):
        return f"檔名不屬於 {expected_sample}：{path.name}"
    keys = ("topology", "in_bam", "assignments", "lineage_paths", "sidecar", "out_bam")
    values = [(k, doc.get(k)) for k in keys if doc.get(k)]
    if not values:
        return f"{path.name} 缺 input identity 欄位"
    wrong = [f"{k}={v}" for k, v in values if not sample_in_path(v, expected_sample)]
    if wrong:
        return f"{path.name} 含異樣本路徑：" + "；".join(wrong)
    return ""


def load(reg: Registry, pre_dir, post_dir, expected_sample=None) -> Capability:
    cap = make(reg, "lca_ab", "LCA 增益 A/B（pre / post receipt）",
               enables=["kpi:lca_gain", "selfcheck:C11"])
    pre = Path(pre_dir) if pre_dir else None
    post = Path(post_dir) if post_dir else None
    if not pre or not pre.is_dir():
        cap.state = "ABSENT"
        probe(cap, "S0", False, f"pre-LCA receipt 目錄不存在：{pre}")
        return cap
    all_pre = sorted(pre.glob("*.tag_bam.receipt.json"))
    all_post = sorted(post.glob("*.tag_bam.receipt.json")) if post and post.is_dir() else []
    if expected_sample:
        pre_f = [f for f in all_pre if sample_in_path(f.name, str(expected_sample))]
        post_f = [f for f in all_post if sample_in_path(f.name, str(expected_sample))]
    else:
        pre_f, post_f = all_pre, all_post
    if not pre_f or not post_f:
        cap.state = MALFORMED if expected_sample and (all_pre or all_post) else "ABSENT"
        probe(cap, "S0", False,
              f"sample={expected_sample or '*'}：pre {len(pre_f)} 檔 / post {len(post_f)} 檔，"
              "至少一邊為空；不會借用其他樣本 receipt")
        return cap
    probe(cap, "S0", True, f"pre {len(pre_f)} 檔、post {len(post_f)} 檔")

    docs = {}
    for f in pre_f + post_f:
        cap.paths.append(file_ref(f))
        try:
            doc = _read_doc(f)
        except (OSError, ValueError) as exc:
            cap.state = MALFORMED
            probe(cap, "S1", False, f"{f.name} 無法解析：{exc}")
            return cap
        if expected_sample:
            err = _identity_error(f, doc, str(expected_sample))
            if err:
                cap.state = MALFORMED
                probe(cap, "S1", False, err + "；已 fail-closed")
                return cap
        docs[str(f)] = doc
    if expected_sample:
        probe(cap, "S1", True,
              f"pre/post 共 {len(pre_f) + len(post_f)} 份 receipt 的檔名與 input identity "
              f"均屬於 {expected_sample}")

    # 兩側檔數不必相等 —— post 端可能多跑了幾條染色體（例如 2026-08-08 補上
    # chrX/chrY 讓輸出 BAM 成為來源的完整複本）。若直接把兩側全部加總相比，
    # 多出來的染色體會讓 reads_total / hp_written 等非 LCA 欄位也不同，
    # 看起來像「兩套不是同輸入」，其實只是分母不同。
    # 正確做法是**按染色體取交集**後再比，這才是真正的同輸入 A/B。
    def _by_chrom(files):
        out = {}
        for f in files:
            m = re.search(r"\.(chr[0-9XYM]+)\.", Path(f).name)
            if m:
                out[m.group(1)] = f
        return out

    pre_c, post_c = _by_chrom(pre_f), _by_chrom(post_f)
    shared = sorted(set(pre_c) & set(post_c))
    only_post = sorted(set(post_c) - set(pre_c))
    only_pre = sorted(set(pre_c) - set(post_c))
    if not shared:
        cap.state = "PARTIAL"
        probe(cap, "S1", False, "兩側沒有共同的染色體，無法對照")
        return cap
    extra = ""
    if only_post:
        extra += f"；只有 post 有：{'、'.join(only_post)}（已排除）"
    if only_pre:
        extra += f"；只有 pre 有：{'、'.join(only_pre)}（已排除）"
    probe(cap, "S1", True, f"按染色體取交集 {len(shared)} 條，可做同輸入對照{extra}")

    a, b = _agg([pre_c[c] for c in shared]), _agg([post_c[c] for c in shared])
    keys = set(a) | set(b)
    diff = sorted(k for k in keys if a.get(k) != b.get(k))
    same = len(keys) - len(diff)
    unexpected = [k for k in diff if k not in LCA_FIELDS]
    # 聚合值相同不代表逐染色體相同（不同 chromosome 可相互抵銷）。
    # 因此另對每個 shared chromosome 比較非 LCA 欄位。
    per_chrom_unexpected = {}
    per_chrom_input_mismatch = {}
    for chrom in shared:
        da, db = docs[str(pre_c[chrom])], docs[str(post_c[chrom])]
        sa = da.get("stats") or {}
        sb = db.get("stats") or {}
        bad = sorted(k for k in set(sa) | set(sb)
                     if k not in LCA_FIELDS and sa.get(k) != sb.get(k))
        if bad:
            per_chrom_unexpected[chrom] = bad
        bad_inputs = sorted(k for k in A_B_IDENTITY_FIELDS if da.get(k) != db.get(k))
        if bad_inputs:
            per_chrom_input_mismatch[chrom] = bad_inputs
    check_ok = not unexpected and not per_chrom_unexpected and not per_chrom_input_mismatch
    probe(cap, "S2", check_ok,
          f"交集 {len(shared)} 條染色體；聚合後 {same} 個欄位相同，"
          f"不同：{'、'.join(diff) or '無'}；"
          f"逐染色體非 LCA stats mismatch {len(per_chrom_unexpected)} 條；"
          f"input identity mismatch {len(per_chrom_input_mismatch)} 條"
          + (f"；stats mismatch chromosome：{'、'.join(per_chrom_unexpected)}"
             if per_chrom_unexpected else "")
          + (f"；input mismatch chromosome：{'、'.join(per_chrom_input_mismatch)}"
             if per_chrom_input_mismatch else ""))
    if not check_ok:
        cap.state = "PARTIAL"
        fields = sorted(set(unexpected) |
                        {x for xs in per_chrom_unexpected.values() for x in xs} |
                        {x for xs in per_chrom_input_mismatch.values() for x in xs})
        cap.reason = ("除了 LCA 相關欄位，" + "、".join(fields) +
                      " 也不同 —— 兩套不是同輸入，增益數字不可解讀為純 LCA 效果。")

    pre_lv, post_lv = a.get("lv_written", 0), b.get("lv_written", 0)
    gain = (post_lv / pre_lv) if pre_lv else 0
    probe(cap, "S3", True,
          f"lv 覆蓋 {pre_lv:,} → {post_lv:,}（{gain:.3f}×）；"
          f"LCA 解出 {b.get('lca_resolved', 0):,} 條")

    cap.counts.update({"pre_files": len(pre_f), "post_files": len(post_f),
                       "shared_files": len(shared), "only_pre_files": len(only_pre),
                       "only_post_files": len(only_post),
                       "pre_lv": pre_lv, "post_lv": post_lv,
                       "lca_resolved": b.get("lca_resolved", 0)})
    cap.payload = {"n_files": len(shared), "pre_files": len(pre_f),
                   "post_files": len(post_f), "shared_files": len(shared),
                   "only_pre": only_pre, "only_post": only_post,
                   "identical_fields": same,
                   "differing_fields": diff, "unexpected": unexpected,
                   "per_chrom_unexpected": per_chrom_unexpected,
                   "per_chrom_input_mismatch": per_chrom_input_mismatch,
                   "pre_lv": pre_lv, "post_lv": post_lv, "gain": gain,
                   "pre": a, "post": b}
    return finalize(cap)
