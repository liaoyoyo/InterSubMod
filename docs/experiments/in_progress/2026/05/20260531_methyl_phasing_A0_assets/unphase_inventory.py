#!/usr/bin/env python3
"""
real unphase / HP3 read 可救援量清點（漏斗分類，read-only）。

回答用戶：實際 unphase 量、可修正歸類量、無法區分量 + HP3 統計。
單 pass 掃 chr，2kb bin 記每 read：tag(1/2/unphase/3) + 是否有 >=5 CpG 甲基。
每 bin: n_hp1 / n_hp2 anchor + n_unphase / n_unphase_cpg / n_h3 / n_h3_cpg。
漏斗（per chr → 後合計推全基因組）：
  U_total → U_有甲基(>=5CpG) → U_有本地參考(bin 內 hp1>=8 & hp2>=8) → U_兩者皆備(可嘗試救援)
  無法區分 = U_total − U_兩者皆備。
數字落 JSON；報告 Read 後才引用。
"""
import sys, json, argparse
import pysam

BAM = "/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam"

def get_hp(read):
    try: return str(read.get_tag("HP"))
    except KeyError: return "unphase"

def has_min_cpg(read, k=5):
    mods = read.modified_bases or {}
    for key, calls in mods.items():
        if key[0] in ("C", b"C") and key[2] in ("m", 27551):
            return len(calls) >= k
    return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--bin", type=int, default=2000)
    ap.add_argument("--min-anchor", type=int, default=8)
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()
    bam = pysam.AlignmentFile(BAM, "rb")
    bins = {}  # binidx -> dict counts
    seen = set()
    for read in bam.fetch(args.chrom):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        if read.query_name in seen:
            continue
        seen.add(read.query_name)
        hp = get_hp(read)
        b = (read.reference_start) // args.bin
        d = bins.setdefault(b, {"hp1": 0, "hp2": 0, "u": 0, "u_cpg": 0, "h3": 0, "h3_cpg": 0})
        if hp == "1": d["hp1"] += 1
        elif hp == "2": d["hp2"] += 1
        elif hp == "unphase":
            d["u"] += 1
            if has_min_cpg(read): d["u_cpg"] += 1
        elif hp == "3":
            d["h3"] += 1
            if has_min_cpg(read): d["h3_cpg"] += 1
    bam.close()
    U = U_cpg = U_anchored = U_cpg_anchored = 0
    H3 = H3_cpg = H3_anchored = H3_cpg_anchored = 0
    for b, d in bins.items():
        anchored = (d["hp1"] >= args.min_anchor and d["hp2"] >= args.min_anchor)
        U += d["u"]; U_cpg += d["u_cpg"]
        H3 += d["h3"]; H3_cpg += d["h3_cpg"]
        if anchored:
            U_anchored += d["u"]; U_cpg_anchored += d["u_cpg"]
            H3_anchored += d["h3"]; H3_cpg_anchored += d["h3_cpg"]
    out = {"chrom": args.chrom, "bin_size": args.bin, "min_anchor": args.min_anchor,
           "unphase": {"total": U, "with_cpg>=5": U_cpg, "in_anchored_bin": U_anchored,
                       "attemptable(cpg+anchored)": U_cpg_anchored,
                       "cannot_distinguish": U - U_cpg_anchored},
           "h3": {"total": H3, "with_cpg>=5": H3_cpg, "in_anchored_bin": H3_anchored,
                  "attemptable(cpg+anchored)": H3_cpg_anchored}}
    outp = f"{args.out_dir}/unphase_inventory_{args.chrom}.json"
    json.dump(out, open(outp, "w"), indent=1, ensure_ascii=False)
    u = out["unphase"]
    print(f"[inv] {args.chrom}: unphase total={U} cpg={U_cpg}({100*U_cpg/max(1,U):.0f}%) anchored={U_anchored}({100*U_anchored/max(1,U):.0f}%) attemptable={U_cpg_anchored}({100*U_cpg_anchored/max(1,U):.0f}%) | H3 total={H3} attemptable={H3_cpg_anchored}({100*H3_cpg_anchored/max(1,H3):.0f}%) -> {outp}")

if __name__ == "__main__":
    main()
