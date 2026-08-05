"""共用 region 載入器 (子任務 #1).
讀 ISM per-region 輸出 (methylation.csv read×CpG + reads.tsv 標籤 + cpg_sites.tsv),
以 read_id join。供想法1/想法3 複用。read-only。
"""
import glob
import numpy as np
import pandas as pd

RUN_TP = "output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2/intersubmod_tp"

# BRCA2 pilot: RegionID -> (chr, snv_pos)  (5 TP region, canonical=22305)
BRCA2_REGIONS = {
    22305: ("chr13", 32315128),
    22306: ("chr13", 32317522),
    22307: ("chr13", 32324831),
    22308: ("chr13", 32339132),
    22309: ("chr13", 32345630),
}


def find_region_dir(chrom, pos, run_tp=RUN_TP):
    hits = glob.glob(f"{run_tp}/filtered_snv_tp/{chrom}/{chrom}_{pos}/{chrom}_*")
    return hits[0] if hits else None


def load_region(region_dir):
    """回傳 dict: meth(DataFrame read_id×CpGpos), reads(DataFrame), snv_pos(int)."""
    meth = pd.read_csv(f"{region_dir}/methylation/methylation.csv")
    meth = meth.set_index("read_id")
    meth.columns = [int(c) for c in meth.columns]  # CpG genomic positions
    reads = pd.read_csv(f"{region_dir}/reads/reads.tsv", sep="\t")
    reads["hp"] = reads["hp"].astype(str)
    snv_pos = None
    for line in open(f"{region_dir}/metadata.txt"):
        if line.startswith("SNV Position"):
            snv_pos = int(line.split(":")[-1].strip())
    return {"meth": meth, "reads": reads, "snv_pos": snv_pos, "dir": region_dir}


def tumor_axis_groups(reads):
    """tumor reads (is_tumor==1) 的三軸分群 -> read_id list。
    HP-family: HP1={1,1-1} vs HP2={2}; HP-fine: 1/1-1/2; subclone: HP1germ(1) vs HP1carrier(1-1);
    allele: ALT vs REF (排除 UNKNOWN)。"""
    t = reads[reads.is_tumor == 1]
    g1 = t[t.hp == "1"].read_id.tolist()
    g11 = t[t.hp == "1-1"].read_id.tolist()
    g2 = t[t.hp == "2"].read_id.tolist()
    alt = t[t.alt_support == "ALT"].read_id.tolist()
    ref = t[t.alt_support == "REF"].read_id.tolist()
    return {
        "hpfam": [("HP1", g1 + g11), ("HP2", g2)],
        "hpfine": [("HP1", g1), ("HP1-1", g11), ("HP2", g2)],
        "subclone": [("HP1germ", g1), ("HP1carrier", g11)],
        "allele": [("ALT", alt), ("REF", ref)],
    }


def normal_hp_groups(reads):
    """normal reads (is_tumor==0) 的 HP1 vs HP2 分群 (normal-anchored ASM-control 用)。"""
    n = reads[reads.is_tumor == 0]
    hp1 = n[n.hp.isin(["1", "1-1"])].read_id.tolist()
    hp2 = n[n.hp == "2"].read_id.tolist()
    return hp1, hp2
