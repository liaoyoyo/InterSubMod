#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
從 COSMIC Cell Lines Project CompleteCNA 建 per-sample CN BED(2026-06-29 多樣本 CN 補洞)。
⚠ CLP CN = gene-level 稀疏 aberration 清單(非全基因組 segment)→ BED 只含 gain/loss/loh 區間;
   未列基因 = 未評估(下游 SM_CN_SPARSE=1 → 無重疊回 'unknown' 非 'neutral',誠實)。
label: MUT_TYPE=gain → gain;loss & MINOR_ALLELE=0 → loh(hemizygous);loss & minor>0 → loss。
用法: python3 clp_cn_bed.py  → 每樣本寫 {workdir}/{sample}_clp_cn.bed + 印覆蓋統計。
"""
import os, gzip, csv

ANN = os.environ.get("SM_ANNOT_DIR", "/big7_disk/liaoyoyo2001/gene_annotation")
CNA = f"{ANN}/cosmic_v104/CellLinesProject_CompleteCNA_v104_GRCh38.tsv.gz"
MSROOT = "/big7_disk/liaoyoyo2001/big7_disk_output/multisample_subclone"
# 我們的樣本名 → COSMIC CLP SAMPLE_NAME
NAME_MAP = {"COLO829": "COLO-829", "H1437": "NCI-H1437", "H2009": "NCI-H2009",
            "HCC1937": "HCC1937", "HCC1954": "HCC1954", "HCC1395_DORADO": "HCC1395"}


def label(mut_type, minor):
    if mut_type == "gain":
        return "gain"
    if mut_type == "loss":
        return "loh" if minor == "0" else "loss"
    return None


def main():
    # 收集每個 COSMIC 名的 CN 區間
    by_cos = {}
    with gzip.open(CNA, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cn = r.get("SAMPLE_NAME", "").strip()
            lab = label(r.get("MUT_TYPE", "").strip(), r.get("MINOR_ALLELE", "").strip())
            chrom = r.get("CHROMOSOME", "").strip()
            try:
                s = int(r.get("GENOME_START", "")); e = int(r.get("GENOME_STOP", ""))
            except ValueError:
                continue
            if not lab or not chrom:
                continue
            by_cos.setdefault(cn, []).append((f"chr{chrom}", s, e, lab))
    # 每樣本寫 BED
    for ours, cos in NAME_MAP.items():
        d = os.path.join(MSROOT, ours)
        if not os.path.isdir(d):
            print(f"  {ours}: workdir 不存在,跳過"); continue
        ivs = sorted(by_cos.get(cos, []))
        bed = os.path.join(d, f"{ours}_clp_cn.bed")
        with open(bed, "w") as fh:
            for chrom, s, e, lab in ivs:
                fh.write(f"{chrom}\t{s}\t{e}\t{lab}\n")
        from collections import Counter
        c = Counter(x[3] for x in ivs)
        print(f"  {ours} (COSMIC={cos}): {len(ivs)} CN 區間 {dict(c)} → {bed}")


if __name__ == "__main__":
    main()
