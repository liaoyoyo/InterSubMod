#!/usr/bin/env python3
"""
unphase 救援原型 (a+c)：對 longphase-S tagged BAM 中『真正無 HP tag』的 read，
用甲基把它救回 HP1/HP2，產 reassigned_hp + confidence（不覆寫原 tag，可回溯）。

與外推驗證的差異：
  外推 = 「假裝 unphase」(held-out，read 自己有 germline 證據被遮，有真 HP 可對照)
  本原型 = 對『真正無 HP tag』的 read 救援（沒有真 HP，這才是實際要救的對象）

方法（per ±2kb 局部窗，沿用已驗證的甲基質心法）：
  1. 窗內 read 分兩群：anchor = 有 HP:Z:1/2 的 read（已知 HP）；target = 無 HP tag（真 unphase）。
  2. anchor 建 HP1/HP2 甲基質心。
  3. 每個 target read：算到 HP1/HP2 質心的甲基距離 → 預測 HP + confidence。
     confidence = |d1-d2|/(d1+d2)（affinity margin, [0,1]）× coverage_factor(min(1, n_cpg/8))。
  4. 輸出：reassigned_hp / original_hp(="."無tag) / reassign_conf / n_cpg / n_anchor。
  5. 守恆：救回數 + 無法救(無甲基/anchor不足)數 = 總 target 數。

設計對齊 design_02：不覆寫原 tag、保留 original_hp + reassign_source + conf 三欄、雙口徑(全 target vs 有甲基覆蓋)。
唯讀 BAM。輸出 unphase_rescue_{chrom}.json + per-read TSV（供熱圖/IGV）。
"""
import sys, json, argparse, bisect
import numpy as np
import pysam

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"
GVCF = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz"
AD = "/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets"

def get_hp(read):
    try: return str(read.get_tag("HP"))
    except KeyError: return "."  # 真 unphase = 無 HP tag

def read_meth(read, start, end):
    mods = read.modified_bases or {}
    mc = None
    for k, calls in mods.items():
        if k[0] in ("C", b"C") and k[2] in ("m", 27551):
            mc = calls; break
    if not mc: return None
    q2r = {q: r for q, r in read.get_aligned_pairs(matches_only=True)}
    meth = {}
    for qpos, qual in mc:
        rpos = q2r.get(qpos)
        if rpos is None or rpos < start or rpos > end: continue
        meth[rpos] = qual / 255.0
    return meth if meth else None

def rescue_window(bam, chrom, ws, we, min_anchor=8, min_cpg=5):
    """回傳 (anchor_reads, target_results, positions)。"""
    anchors = []  # (hp, read_id, meth)
    targets = []  # (read_id, meth, mapq)
    seen = set()
    allpos = set()
    for read in bam.fetch(chrom, ws, we):
        if read.is_secondary or read.is_supplementary or read.is_unmapped: continue
        if read.query_name in seen: continue
        hp = get_hp(read)
        meth = read_meth(read, ws, we)
        if meth is None or len(meth) < min_cpg:
            if hp == "." :  # 真 unphase 但無甲基 → 記為無法救
                pass
            continue
        seen.add(read.query_name)
        allpos.update(meth.keys())
        if hp in ("1", "2"):
            anchors.append((hp, read.query_name, meth))
        elif hp == ".":
            targets.append((read.query_name, meth, read.mapping_quality))
    n1 = sum(1 for h, _, _ in anchors if h == "1")
    n2 = sum(1 for h, _, _ in anchors if h == "2")
    if n1 < min_anchor or n2 < min_anchor or not targets:
        return None
    positions = sorted(allpos)
    cov = {p: sum(1 for _, _, m in anchors if p in m) for p in positions}
    positions = [p for p in positions if cov[p] >= 0.3 * len(anchors)]
    if len(positions) < min_cpg:
        return None
    def vec(meth):
        return np.array([meth.get(p, np.nan) for p in positions])
    c1 = np.nanmean(np.array([vec(m) for h, _, m in anchors if h == "1"]), axis=0)
    c2 = np.nanmean(np.array([vec(m) for h, _, m in anchors if h == "2"]), axis=0)
    results = []
    for rid, meth, mapq in targets:
        v = vec(meth)
        mask = ~np.isnan(v) & ~np.isnan(c1) & ~np.isnan(c2)
        if mask.sum() < 3:
            results.append({"read_id": rid, "reassigned_hp": ".", "original_hp": ".",
                            "reassign_conf": 0.0, "n_cpg_used": int(mask.sum()),
                            "reason": "insufficient_shared_cpg", "mapq": mapq})
            continue
        d1 = np.sqrt(np.mean((v[mask]-c1[mask])**2))
        d2 = np.sqrt(np.mean((v[mask]-c2[mask])**2))
        pred = "1" if d1 < d2 else "2"
        margin = abs(d1-d2)/(d1+d2) if (d1+d2) > 0 else 0
        cov_factor = min(1.0, mask.sum()/8.0)
        conf = round(margin * cov_factor, 4)
        results.append({"read_id": rid, "reassigned_hp": pred, "original_hp": ".",
                        "reassign_source": "methyl_centroid", "reassign_conf": conf,
                        "n_cpg_used": int(mask.sum()), "mapq": mapq,
                        "dist_hp1": round(float(d1), 4), "dist_hp2": round(float(d2), 4)})
    return {"n_anchor_hp1": n1, "n_anchor_hp2": n2, "n_cpg": len(positions),
            "n_target": len(targets), "results": results}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--n-windows", type=int, default=30)
    ap.add_argument("--win", type=int, default=2000)
    ap.add_argument("--conf-threshold", type=float, default=0.1, help="conf>=此值才算『救回』")
    args = ap.parse_args()
    chrom = args.chrom
    rng = np.random.RandomState(20260603)

    # 取窗：用 germline het 密集處（確保有 anchor）
    vcf = pysam.VariantFile(GVCF)
    hets = []
    for rec in vcf.fetch(chrom):
        s0 = rec.samples[0]; gt = s0.get("GT")
        if gt and len(gt) == 2 and gt[0] != gt[1] and None not in gt:
            hets.append(rec.pos)
    if len(hets) < 10:
        json.dump({"chrom": chrom, "n_windows": 0, "note": "too few hets"},
                  open(f"{AD}/unphase_rescue_{chrom}.json", "w"))
        print(f"[rescue] {chrom}: too few hets"); return
    centers = [hets[i] for i in rng.choice(len(hets), min(args.n_windows*3, len(hets)), replace=False)]

    bam = pysam.AlignmentFile(BAM, "rb")
    windows = []
    all_rescued = []
    for c in centers:
        if len(windows) >= args.n_windows: break
        r = rescue_window(bam, chrom, c - args.win, c + args.win)
        if r:
            r["center"] = int(c)
            windows.append(r)
            all_rescued.extend(r["results"])
    bam.close()

    # 統計
    rescuable = [x for x in all_rescued if x["reassigned_hp"] in ("1", "2")]
    high_conf = [x for x in rescuable if x["reassign_conf"] >= args.conf_threshold]
    out = {
        "chrom": chrom, "n_windows": len(windows),
        "n_target_total": len(all_rescued),
        "n_assigned": len(rescuable),
        "n_high_conf": len(high_conf),
        "frac_assigned": round(len(rescuable)/len(all_rescued), 3) if all_rescued else None,
        "frac_high_conf": round(len(high_conf)/len(all_rescued), 3) if all_rescued else None,
        "conf_median": round(float(np.median([x["reassign_conf"] for x in rescuable])), 4) if rescuable else None,
        "hp_balance": {"to_HP1": sum(1 for x in rescuable if x["reassigned_hp"] == "1"),
                       "to_HP2": sum(1 for x in rescuable if x["reassigned_hp"] == "2")},
        "windows": windows,
    }
    json.dump(out, open(f"{AD}/unphase_rescue_{chrom}.json", "w"), ensure_ascii=False, indent=2)
    print(f"[rescue] {chrom}: {len(windows)}窗 target={len(all_rescued)} "
          f"assigned={len(rescuable)}({out['frac_assigned']}) high_conf={len(high_conf)} "
          f"conf_med={out['conf_median']} -> unphase_rescue_{chrom}.json")

if __name__ == "__main__":
    main()
