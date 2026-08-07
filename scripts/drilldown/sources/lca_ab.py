#!/usr/bin/env python3
"""LCA 增益的 A/B 對照（擴充能力）。

LongLineage 的 tag_bam 跑了兩輪：一輪不做 LCA、一輪做。兩套 receipt 除了
LCA 相關欄位外完全相同，所以是乾淨的同輸入 A/B —— 這是唯一能回答
「LCA 到底買到什麼」的統計（實測 lv 覆蓋 86,037 → 427,519 = 4.97×）。

先前一份查證誤判「無來源檔可 grep、不能搬」，因為它用
`lv_written − lca_resolved` 去湊 86,037（85,700，差 337）。那不是那個數的定義；
86,037 是 pre-LCA 那套 receipt 的 lv_written 加總，可逐檔 grep。
"""
from __future__ import annotations

import json
from pathlib import Path

from capability import Capability, Registry, finalize, make, probe

LCA_FIELDS = {"lv_written", "lca_resolved", "lca_candidates_sum"}


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


def load(reg: Registry, pre_dir, post_dir) -> Capability:
    cap = make(reg, "lca_ab", "LCA 增益 A/B（pre / post receipt）",
               enables=["kpi:lca_gain", "selfcheck:C11"])
    pre = Path(pre_dir) if pre_dir else None
    post = Path(post_dir) if post_dir else None
    if not pre or not pre.is_dir():
        cap.state = "ABSENT"
        probe(cap, "S0", False, f"pre-LCA receipt 目錄不存在：{pre}")
        return cap
    pre_f = sorted(pre.glob("*.json"))
    post_f = sorted(post.glob("*.tag_bam.receipt.json")) if post and post.is_dir() else []
    if not pre_f or not post_f:
        cap.state = "ABSENT"
        probe(cap, "S0", False, f"pre {len(pre_f)} 檔 / post {len(post_f)} 檔，至少一邊為空")
        return cap
    probe(cap, "S0", True, f"pre {len(pre_f)} 檔、post {len(post_f)} 檔")

    if len(pre_f) != len(post_f):
        cap.state = "PARTIAL"
        probe(cap, "S1", False, f"檔數不等（{len(pre_f)} vs {len(post_f)}），A/B 不對稱")
    else:
        probe(cap, "S1", True, "檔數相等，可做同輸入對照")

    a, b = _agg(pre_f), _agg(post_f)
    keys = set(a) | set(b)
    diff = sorted(k for k in keys if a.get(k) != b.get(k))
    same = len(keys) - len(diff)
    unexpected = [k for k in diff if k not in LCA_FIELDS]
    probe(cap, "S2", True,
          f"{same} 個欄位完全相同；不同的：{'、'.join(diff) or '無'}"
          + (f"　⚠ 非 LCA 欄位也不同：{'、'.join(unexpected)}" if unexpected else ""))
    if unexpected:
        cap.state = "PARTIAL"
        cap.reason = ("除了 LCA 相關欄位，" + "、".join(unexpected) +
                      " 也不同 —— 兩套不是同輸入，增益數字不可解讀為純 LCA 效果。")

    pre_lv, post_lv = a.get("lv_written", 0), b.get("lv_written", 0)
    gain = (post_lv / pre_lv) if pre_lv else 0
    probe(cap, "S3", True,
          f"lv 覆蓋 {pre_lv:,} → {post_lv:,}（{gain:.3f}×）；"
          f"LCA 解出 {b.get('lca_resolved', 0):,} 條")

    cap.counts.update({"pre_lv": pre_lv, "post_lv": post_lv,
                       "lca_resolved": b.get("lca_resolved", 0)})
    cap.payload = {"n_files": len(pre_f), "identical_fields": same,
                   "differing_fields": diff, "unexpected": unexpected,
                   "pre_lv": pre_lv, "post_lv": post_lv, "gain": gain,
                   "pre": a, "post": b}
    return finalize(cap)
