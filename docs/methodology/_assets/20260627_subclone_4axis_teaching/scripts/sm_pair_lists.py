"""[pairwise 列表檔] 把 powered+somatic 的 2-locus 對, 按 config × HP 寫成 TSV 列表檔 (lists/)。
另出 CN-clean 子集 (LOH+neutral, 排除 CN-gain 混淆) = 可寫論文的乾淨集。+ per-sSNV census TSV。
"""
import json
import os
from collections import defaultdict, Counter

A = os.environ.get("SM_WORKDIR", "/big7_disk/liaoyoyo2001/InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot")
CNBED = os.environ.get("SM_CNBED", "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed")
LD = f"{A}/lists"
os.makedirs(LD, exist_ok=True)


def load_cn():
    iv = defaultdict(list)
    if not CNBED or not os.path.exists(CNBED):
        return iv
    for ln in open(CNBED):
        p = ln.split()
        if len(p) < 4 or not p[1].isdigit():
            continue
        iv[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for c in iv:
        iv[c].sort()
    return iv


def cn_state(cn, chrom, pos):
    for s, e, lab in cn.get(chrom, []):
        if s <= pos <= e:
            return lab
    return "neutral"


def config_of(c):
    k = (int(c["RA"] > 0), int(c["AR"] > 0), int(c["AA"] > 0))
    return {(1, 1, 0): "mutual_excl", (1, 0, 1): "nested_a_in_b", (0, 1, 1): "nested_b_in_a",
            (0, 0, 1): "co_linked", (1, 1, 1): "independent"}.get(k, "degenerate")


HDR = ["chrom", "pos_a", "pos_b", "dist", "RR", "RA", "AR", "AA", "rel", "config",
       "hp_a", "hp_b", "same_hp", "coread", "vaf_a", "vaf_b", "cn_a", "src_a", "src_b"]


def main():
    d = json.load(open(f"{A}/sm_linkage_genomewide.json"))
    cen = d["census"]
    cn = load_cn()
    files = defaultdict(list)
    files_clean = defaultdict(list)
    cnt = Counter()
    cnt_clean = Counter()
    for p in d["pairs"]:
        if not (p["powered"] and p["somatic_a"] and p["somatic_b"]):
            continue
        cfg = config_of(p["cells"])
        if cfg == "degenerate":
            continue
        hp = "sameHP" if p["same_hp"] else "diffHP"
        cst = cn_state(cn, p["chrom"], p["a"])
        va = cen[f"{p['chrom']}:{p['a']}"].get("vaf")
        vb = cen[f"{p['chrom']}:{p['b']}"].get("vaf")
        row = [p["chrom"], p["a"], p["b"], p["dist"], p["cells"]["RR"], p["cells"]["RA"],
               p["cells"]["AR"], p["cells"]["AA"], p["rel"], cfg, p["hp_a"], p["hp_b"],
               p["same_hp"], p["coread"], va, vb, cst, p["src_a"], p["src_b"]]
        key = f"{cfg}_{hp}"
        files[key].append(row)
        cnt[key] += 1
        if cst in ("loh", "neutral"):
            files_clean[key].append(row)
            cnt_clean[key] += 1

    def write(fn, rows):
        with open(f"{LD}/{fn}", "w") as fh:
            fh.write("\t".join(HDR) + "\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")

    for key, rows in files.items():
        write(f"{key}.tsv", rows)
    for key, rows in files_clean.items():
        write(f"{key}.CLEAN_loh_neutral.tsv", rows)

    # per-sSNV census TSV
    with open(f"{LD}/per_sSNV_census.tsv", "w") as fh:
        fh.write("locus\tsrc\tnormal_ref\tnormal_alt\tsomatic\tvaf\tn_partners_le50k\tn_realized_links\tmax_coread\tcn_state\n")
        for k, c in cen.items():
            cst = cn_state(cn, c["chrom"], c["pos"])
            fh.write(f'{k}\t{c["src"]}\t{c["normal_ref"]}\t{c["normal_alt"]}\t{c["somatic"]}\t{c["vaf"]}\t'
                     f'{c["n_partners_le50k"]}\t{c["n_realized_links"]}\t{c["max_coread"]}\t{cst}\n')

    print("=== pairwise list files (powered+somatic) → lists/ ===")
    print(f"{'config_HP':24s} {'all':>7} {'CLEAN(loh+neutral)':>18}")
    for key in sorted(cnt):
        print(f"{key:24s} {cnt[key]:>7} {cnt_clean.get(key,0):>18}")
    print(f"\ntotal files: {len(files)+len(files_clean)+1}  + per_sSNV_census.tsv ({len(cen)} sSNV)")
    json.dump({"counts_all": dict(cnt), "counts_clean_loh_neutral": dict(cnt_clean)},
              open(f"{A}/sm_pair_lists_counts.json", "w"), ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
