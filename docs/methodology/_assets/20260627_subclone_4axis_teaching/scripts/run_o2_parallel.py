#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""sm_linkage O2 分片平行 driver(2026-07-01):
切片(gap>50kb 邊界,不切群=byte-identical)→ 平行跑 worker(吃滿核)→ 合併 → 記錄 time/memory 觀察數據。
用法: run_o2_parallel.py <chroms|all> <n_shards> <out_json> [maxpar]
env: SM_TBAM/SM_NBAM/SM_VD/SM_VCF_MODE。輸出: out_json(完整 sm_linkage 格式) + out_json.obs(觀察數據)。"""
import os, sys, glob, time, json, subprocess
from collections import defaultdict, Counter
SCD = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCD)
import sm_linkage_genomewide as M

CHROMS = sys.argv[1].split(",") if sys.argv[1] != "all" else [f"chr{c}" for c in range(1, 23)]
N_SHARDS = int(sys.argv[2]); OUT = sys.argv[3]
MAXPAR = int(sys.argv[4]) if len(sys.argv) > 4 else 44
WORKER = f"{SCD}/sm_linkage_o2_shard.py"
PY = sys.executable
TMP = os.path.dirname(OUT) or "."

def shard_ranges(pos, k):
    n = len(pos)
    if n == 0: return []
    if k <= 1: return [(pos[0], pos[-1])]
    gapcuts = [i for i in range(n - 1) if pos[i + 1] - pos[i] > 50000]
    if not gapcuts: return [(pos[0], pos[-1])]
    k = min(k, len(gapcuts) + 1)
    chosen = []
    for j in range(1, k):
        ideal = n * j / k
        c = min((c for c in gapcuts if c not in chosen), key=lambda c: abs(c - ideal))
        chosen.append(c)
    chosen = sorted(set(chosen))
    ranges = []; prev = 0
    for c in chosen:
        ranges.append((pos[prev], pos[c])); prev = c + 1
    ranges.append((pos[prev], pos[-1]))
    return ranges

# 1) 計算切片:每 chrom sSNV 數比例分配 shard,總計 ~N_SHARDS
counts = {}; poslist = {}
for ch in CHROMS:
    snvs = M.load_union(ch); poslist[ch] = [s[0] for s in snvs]; counts[ch] = len(snvs)
tot = sum(counts.values()) or 1
shards = []  # (chrom, lo, hi)
for ch in CHROMS:
    if counts[ch] == 0: continue
    k = max(1, round(N_SHARDS * counts[ch] / tot))
    for lo, hi in shard_ranges(poslist[ch], k):
        shards.append((ch, lo, hi))
print(f"切片: {len(shards)} 片 (target {N_SHARDS}), chroms={len(CHROMS)}, maxpar={MAXPAR}", flush=True)

# 2) 平行跑 worker(cap MAXPAR)
t0 = time.time()
running = []; done = []; i = 0
shard_out = lambda j: f"{TMP}/_o2shard_{j}.json"
def launch(j):
    ch, lo, hi = shards[j]
    return subprocess.Popen([PY, WORKER, ch, str(lo), str(hi), shard_out(j)],
                            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
while i < len(shards) or running:
    while i < len(shards) and len(running) < MAXPAR:
        running.append((i, launch(i))); i += 1
    time.sleep(0.5)
    still = []
    for j, p in running:
        if p.poll() is None: still.append((j, p))
        else:
            if p.returncode != 0: print(f"⚠ shard {j} rc={p.returncode}: {p.stderr.read().decode()[:200]}", flush=True)
            done.append(j)
    running = still
total_wall = time.time() - t0

# 3) 合併
merged = {"params": {"TIER_R": M.TIER_R, "COREAD_MIN": M.COREAD_MIN, "MAPQ_MIN": M.MAPQ_MIN, "NORM_REF_MIN": M.NORM_REF_MIN},
          "universe": {"n_tp": 0, "n_fp": 0, "n_union": 0, "per_chrom": {}},
          "n_groups_multi": 0, "n_singletons": 0, "pairs": [], "census": {}, "tally_powered_somatic": Counter()}
benches = []
for j in range(len(shards)):
    p = shard_out(j)
    if not os.path.exists(p): continue
    d = json.load(open(p))
    merged["pairs"].extend(d["pairs"]); merged["census"].update(d["census"])
    merged["n_groups_multi"] += d["n_groups_multi"]; merged["n_singletons"] += d["n_singletons"]
    u = d["universe"]; merged["universe"]["n_tp"] += u["n_tp"]; merged["universe"]["n_fp"] += u["n_fp"]; merged["universe"]["n_union"] += u["n_union"]
    for ch, v in u["per_chrom"].items():
        pc = merged["universe"]["per_chrom"].setdefault(ch, {"tp": 0, "fp": 0, "union": 0})
        for kk in ("tp", "fp", "union"): pc[kk] += v[kk]
    for k, v in d["tally_powered_somatic"].items(): merged["tally_powered_somatic"][k] += v
    if os.path.exists(p + ".bench"): benches.append(json.load(open(p + ".bench")))
merged["pairs"].sort(key=lambda x: (x["chrom"], x["a"], x["b"]))   # 決定性排序
merged["tally_powered_somatic"] = dict(merged["tally_powered_somatic"])
merged["n_pairs_recorded"] = len(merged["pairs"])
json.dump(merged, open(OUT, "w"), ensure_ascii=False)

# 4) 觀察數據(time/memory)
sum_wall = sum(b["wall_sec"] for b in benches); max_wall = max((b["wall_sec"] for b in benches), default=0)
peaks = [b["peak_rss_mb"] for b in benches]
obs = {"n_shards": len(shards), "maxpar": MAXPAR, "total_wall_sec": round(total_wall, 1),
       "sum_shard_wall_sec": round(sum_wall, 1), "max_shard_wall_sec": round(max_wall, 1),
       "parallel_efficiency": round(sum_wall / max(total_wall, 1e-9), 1),
       "peak_rss_mb_max_single": round(max(peaks), 1) if peaks else None,
       "peak_rss_mb_sum": round(sum(peaks), 1) if peaks else None,
       "peak_rss_mb_median": round(sorted(peaks)[len(peaks)//2], 1) if peaks else None,
       "total_pairs": merged["n_pairs_recorded"], "universe": merged["universe"]["n_union"],
       "slowest_shards": sorted(benches, key=lambda b: -b["wall_sec"])[:5]}
json.dump(obs, open(OUT + ".obs", "w"), ensure_ascii=False, indent=1)
print(f"DONE total_wall={total_wall:.0f}s sum_shard_wall={sum_wall:.0f}s parallel_eff={obs['parallel_efficiency']}x "
      f"peak_max={obs['peak_rss_mb_max_single']}MB pairs={merged['n_pairs_recorded']} -> {OUT}", flush=True)
# 清理 shard 暫存
for j in range(len(shards)):
    for suf in ("", ".bench"):
        try: os.remove(shard_out(j) + suf)
        except OSError: pass
