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
import re
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
