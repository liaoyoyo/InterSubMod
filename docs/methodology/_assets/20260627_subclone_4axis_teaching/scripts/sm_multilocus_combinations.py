"""[多位點組合 census] 從 2-locus pairwise 升級到 multi-locus (≥3 sSNV) 的 per-read 基因型組合。
每個 group 取 somatic sSNV (census somatic==True)，tumor per-read 抽各 sSNV 等位 → 對「覆蓋全部 sSNV」的 read
建基因型向量字串 (R/A) → 不同向量 = 局部細胞 population (如 chr17: RRRR=ancestral / γ / α-only / α+β)。
每 group 記 n_sSNV / n_full_cov_reads / n_populations / 各 combo:count / CN state / sameHP。
genome-wide 聚合: 每 group population 數分布 + 比例 (≥2/≥3/≥4) + CN 分層 (clean LOH/neutral vs gain)。
argv: chroms out_path。輸出含 per-group list。
"""
import sys
import json
from collections import Counter, defaultdict
import os
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # import 同夾的(env-var 化) sm_linkage
import sm_linkage_genomewide as M

# CNBED env-var(2026-06-29);空字串或缺檔 → 無 CN(全 unknown,graceful)
CNBED = os.environ.get("SM_CNBED", "/big8_disk/data/HCC1395/SEQC2/CNV/ngs_benchmark_cnvs_gain_loss_loh.bed")
MINREAD = 3          # 一個 population 需 >=3 支持 read
MAX_SNV = 8          # group 內 somatic sSNV 上限 (>8 spread 太大 full-cov~0; 取最密 8 個)


def load_cn():
    iv = defaultdict(list)
    if not CNBED or not os.path.exists(CNBED):
        return iv  # 無 CN → cn_state 回 'neutral'(此處無 LOH/gain 標,下游可標 unknown)
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


def main(chroms, out_path):
    cen = json.load(open(f"{M.A}/sm_linkage_genomewide.json"))["census"]
    cn = load_cn()
    tb = pysam.AlignmentFile(M.TBAM, "rb")
    groups_out = []
    agg = Counter()
    agg_clean = Counter()
    for chrom in chroms:
        snvs = M.load_union(chrom)
        for g in M.make_groups(snvs):
            som = [s for s in g if cen.get(f"{chrom}:{s[0]}", {}).get("somatic") is True]
            if len(som) < 2:
                continue
            # 若 >MAX_SNV, 取最密集連續 MAX_SNV 個 (span 最小) 以利 full-coverage
            if len(som) > MAX_SNV:
                best = min(range(len(som) - MAX_SNV + 1),
                           key=lambda i: som[i + MAX_SNV - 1][0] - som[i][0])
                som = som[best:best + MAX_SNV]
            calls, hp, ps = M.per_read_calls(tb, chrom, som)
            positions = [s[0] for s in som]
            combo = Counter()
            hpset = Counter()
            for rn, c in calls.items():
                if all(c.get(p) in ("REF", "ALT") for p in positions):
                    sstr = "".join("A" if c[p] == "ALT" else "R" for p in positions)
                    combo[sstr] += 1
                    if rn in hp and "A" in sstr:
                        hpset[hp[rn]] += 1
            nfull = sum(combo.values())
            if nfull < MINREAD:
                continue
            pops = {k: v for k, v in combo.items() if v >= MINREAD}
            npop = len(pops)
            # populations with >=1 ALT (non-ancestral) = mutation-bearing 群
            npop_alt = sum(1 for k in pops if "A" in k)
            cst = cn_state(cn, chrom, positions[len(positions) // 2])
            same_hp = (len(hpset) == 1)
            rec = {"chrom": chrom, "start": positions[0], "end": positions[-1],
                   "span": positions[-1] - positions[0], "n_sSNV": len(som),
                   "positions": positions, "n_full_cov_reads": nfull,
                   "n_populations": npop, "n_populations_with_ALT": npop_alt,
                   "populations": dict(sorted(pops.items(), key=lambda x: -x[1])),
                   "cn": cst, "dominant_hp_count": len(hpset), "same_hp": same_hp}
            groups_out.append(rec)
            agg[f"npop_{min(npop, 5)}"] += 1
            agg[f"nsnv_{min(len(som), 6)}"] += 1
            if npop >= 2:
                agg["groups_ge2_pop"] += 1
            if npop >= 3:
                agg["groups_ge3_pop"] += 1
            if cst in ("loh", "neutral"):
                agg_clean[f"npop_{min(npop, 5)}"] += 1
                if npop >= 2:
                    agg_clean["groups_ge2_pop"] += 1
                if npop >= 3:
                    agg_clean["groups_ge3_pop"] += 1
    tb.close()
    out = {"params": {"MINREAD": MINREAD, "MAX_SNV": MAX_SNV},
           "n_groups_analyzed": len(groups_out),
           "aggregate_all": dict(agg), "aggregate_clean_loh_neutral": dict(agg_clean),
           "groups": groups_out}
    json.dump(out, open(out_path, "w"), ensure_ascii=False)
    print(f"DONE chroms={chroms} groups={len(groups_out)} ge2pop={agg['groups_ge2_pop']} "
          f"ge3pop={agg['groups_ge3_pop']} -> {out_path}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1].split(","), sys.argv[2])
