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


def germ_family(hptag):
    """L0 germline HP 家族(2026-07-06 使用者定案:家族優先於算法)。
    HP tag 語義(reference_hp_tag_definition):1/2=germline hap;1-1/2-1/3=integrated somatic。
    germline 家族 = 第一碼:1/1-1→'1'、2/2-1→'2';'3'=somatic-integrated 無明確 germline 家族;None/其他→'none'(unphased)。
    🔴 只有同 germline 家族的 read 才依 sSNV 連接建樹(不同家族=不同親代染色體=分開樹)。"""
    if hptag is None:
        return "none"
    s = str(hptag)
    if s.startswith("1"):
        return "1"
    if s.startswith("2"):
        return "2"
    if s == "3":
        return "3"
    return "none"


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
            # gap#1(2026-07-04):同時記錄 partial-read subcube-groups('X'=未覆蓋位點=cube 自由維度)+
            #   per-column coverage(含 partial),供 group-Steiner 覆蓋條件。full-cov 路徑(combo)不變 → populations 向後相容。
            subread = Counter()
            colcov = {p: [0, 0] for p in positions}  # pos -> [nREF, nALT]
            # L0 per-HP-family split(2026-07-06 使用者定案):樹 per germline 家族分開建。
            #   pooled(combo/subread)保留供對照;下游 layered driver 用 *_by_hp 建 per-family 樹。
            combo_fam = defaultdict(Counter)
            subread_fam = defaultdict(Counter)
            colcov_fam = defaultdict(lambda: {p: [0, 0] for p in positions})
            reads_fam = Counter()
            for rn, c in calls.items():
                vec = "".join("A" if c.get(p) == "ALT" else ("R" if c.get(p) == "REF" else "X") for p in positions)
                ncov = len(positions) - vec.count("X")
                if ncov == 0:
                    continue
                fam = germ_family(hp.get(rn))  # L0: read 的 germline 家族
                reads_fam[fam] += 1
                subread[vec] += 1
                subread_fam[fam][vec] += 1
                for i, p in enumerate(positions):
                    if vec[i] == "R":
                        colcov[p][0] += 1; colcov_fam[fam][p][0] += 1
                    elif vec[i] == "A":
                        colcov[p][1] += 1; colcov_fam[fam][p][1] += 1
                if ncov == len(positions):  # full-cov 路徑不變
                    combo[vec] += 1
                    combo_fam[fam][vec] += 1
                    if rn in hp and "A" in vec:
                        hpset[hp[rn]] += 1
            nfull = sum(combo.values())
            subread_groups = {v: n for v, n in subread.items() if n >= MINREAD}
            # gap#1: 只要有 full-cov population 或 partial subcube-group(≥MINREAD)就保留(原:無 full-cov 即整組丟 → 40.4% 空)
            if nfull < MINREAD and not subread_groups:
                continue
            pops = {k: v for k, v in combo.items() if v >= MINREAD}
            # L0 per-family(≥MINREAD 過濾同 pooled;去空家族)
            pops_by_hp = {f: {k: v for k, v in cc.items() if v >= MINREAD} for f, cc in combo_fam.items()}
            pops_by_hp = {f: d for f, d in pops_by_hp.items() if d}
            subread_by_hp = {f: {k: v for k, v in cc.items() if v >= MINREAD} for f, cc in subread_fam.items()}
            subread_by_hp = {f: d for f, d in subread_by_hp.items() if d}
            colcov_by_hp = {f: {str(p): colcov_fam[f][p] for p in positions} for f in set(pops_by_hp) | set(subread_by_hp)}
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
                   "subread_groups": dict(sorted(subread_groups.items(), key=lambda x: -x[1])),  # gap#1: partial subcube-groups
                   "col_coverage": {str(p): colcov[p] for p in positions},                        # gap#1: per-locus REF/ALT 覆蓋(含 partial)
                   # L0 per-HP-family(2026-07-06):樹 per germline 家族分開建的底料
                   "populations_by_hp": {f: dict(sorted(d.items(), key=lambda x: -x[1])) for f, d in pops_by_hp.items()},
                   "subread_groups_by_hp": {f: dict(sorted(d.items(), key=lambda x: -x[1])) for f, d in subread_by_hp.items()},
                   "col_coverage_by_hp": colcov_by_hp,
                   "reads_by_hp": dict(reads_fam),
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
