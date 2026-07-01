<!--
建立時間: 2026-07-01
類型: 程式觀察紀錄 — sSNV 骨幹 O2(region-pileup)優化 + 分片平行 全基因組 time/memory benchmark
狀態: O2 已套 production(byte-identical 三層驗證);全基因組平行 benchmark 完成
build_branch: research/subclonal-reconstruction-202606
provenance: HCC1395。所有 time/memory 皆本輪實測(o2par_full.json.obs + /usr/bin/time -v + 各測試 log);O2 正確性以 original vs O2 同輸入 diff 證明。
-->

# sSNV 骨幹 O2 優化 + 分片平行 — 程式觀察紀錄（time/memory）

> 目的：記錄 O2(per-group region pileup)優化的**正確性驗證** + **時間/記憶體觀察數據**,作為程式效能基準。
> 對象：`sm_linkage_genomewide.py`(per_read_calls)+ 新增 `sm_linkage_o2_shard.py` / `run_o2_parallel.py`。

---

## §1 O2 優化內容

**瓶頸**（效能審計定論）：`per_read_calls` 對 group 內**每個 sSNV 各發一次 pileup**（per-locus）。ONT read span p95~34kb、group 由 gap≤50kb 串成 → 一條長 read 橫跨 group 內 k 個 sSNV 時被 pileup 引擎**重複解壓 k 次**。大 BAM 隨機讀是機器主瓶頸。

**O2**：改對 group span 發**單次 region pileup**，走訪 columns 時只在 sSNV 位點 tabulate（`col.reference_pos+1 in posset`）→ 每條 read 只解壓一次。語意與逐位點版等價。

---

## §2 正確性驗證（byte-identical,三層,全在 original vs O2 同輸入）

| 層 | 測試 | 結果 |
|---|---|---|
| Unit | chr20-22 **397 個真實 group**,逐 group 比對 calls/hp/ps | **397/397 byte-identical, 0 mismatch** |
| 端到端 | chr22 完整 sm_linkage(原版 main vs O2 main),含 normal is_somatic | **7/7 檢查全過**（pairs/census/tally/universe/n_groups_multi/n_singletons/n_pairs）|
| 分片合併 | chr22 切 2 片平行 + 合併 vs reference | **7/7 byte-identical** |

→ O2 **保正確性**（input-agnostic:等價性是 pileup 函式性質,與 BAM/VCF 內容無關）。

---

## §3 加速觀察（time）

| 情境 | 原版 | O2 | 加速 |
|---|--:|--:|--:|
| per_read_calls unit(chr20-22, 397 group) | 348.9s | 205.3s | **1.70x** |
| 端到端 chr22(含 normal pileup 稀釋) | 153s | 109s | 1.41x |
| 超寬 group chr8 503kb/33-sSNV | 8.7s | 2.9s | **2.98x** |
| 超寬 group chr8 258kb/43-sSNV | 12.0s | 1.9s | **6.19x** |

**關鍵**：即使最寬 group（503kb）O2 也更快（2.98x）——region pileup 對長 read 只解壓一次的省法,勝過「走訪較多 column」的 overhead。無需 span guard。

---

## §4 全基因組分片平行 benchmark（time + memory 觀察數據）

**設定**：45 片（gap>50kb 邊界切,不切群=byte-identical）/ maxpar 44 核 / O2 / canonical longphase_s 輸入。
**觀察數據**（`o2par_full.json.obs` + `/usr/bin/time -v`）：

| 指標 | 值 |
|---|---|
| 總 wall clock | **349.3s（5:54）** |
| sum_shard_wall（所有片 CPU 時間和）| 9968.9s |
| **parallel_efficiency**（=sum/total）| **28.5x**（44 核中 ~65%）|
| CPU 使用率（/usr/bin/time）| **2424%**（~24 核有效平均）|
| peak RSS 單片最大 | 408.7 MB |
| peak RSS 中位 | 366.9 MB |
| peak RSS 總和（並發上界）| 16,647 MB（~16.6 GB）|
| driver 層 peak RSS | 418 MB |
| 最慢片 | chr7 349s / chr15 329s / chr7 318s（大 chrom = 尾端 straggler）|

**解讀**：
- **記憶體充裕**：單片 ~367MB、44 片並發 ~16.6GB,遠低於 503G → 可再加平行度或加大 READ_CAP。
- **平行效率 28.5/44 ≈ 65%**：非 I/O 完全飽和崩潰（若 I/O-bound 到底,效率會遠低）,但受**大 chrom straggler**（chr7/chr15 單片 5+ 分鐘,拖住尾端）+ 部分 I/O 稀釋。
- **改進方向**：把最慢的 chr7/chr15/chr8 再切更細（目前它們仍是單片整條臂）→ 進一步逼近 44x。

---

## §5 輸入說明（誠實邊界）

本 benchmark 用 **canonical longphase_s 輸入**（`filtered_snv_tp/fp.vcf.gz` single mode:TP=29754/FP=627）。
先前 topology reference（worktree）用的是 **`/big8/.../HCC1395/pileup` per-chrom VCF**（不同 TP/FP 集,FP 較多）+ 產它的 20260618 腳本已不存在。故**未**對舊 reference 做 genome-scale diff（雙重版本 confound）；O2 正確性以**同輸入 original vs O2** 證明,與輸入選擇正交。
→ 若日後要 canonical 生產重跑（如 C3/C2 修正 baseline）,須先確認正確的輸入 VCF 集 + BAM。

---

## §6 產物

- `sm_linkage_genomewide.py`：`per_read_calls` 已套 O2（byte-identical）。
- `sm_linkage_o2_shard.py`：分片 worker（O2 + load_union 過濾到 [lo,hi] + 記 wall/peak RSS）。
- `run_o2_parallel.py`：分片平行 driver（gap>50kb 切片 + cap 平行 + 合併 + `.obs` 觀察數據）。
- 用法：`run_o2_parallel.py all 44 <out.json> 44`（env: SM_TBAM/SM_NBAM/SM_VD/SM_VCF_MODE）。

**增量正確性修正（已完成,每步小規模驗證,2026-07-01）**：D4(甲基 overclaim 字串,commit 979c54d)→ C3(去噪 eps floor + gate,979c54d)→ C2(4-gamete 計數旗標,44ed19c)。canonical determinacy A_determined 1812→1741 / incompatible 12→118。詳見 `20260701_ssnv_backbone_method_spec_and_correctness_audit_01.md` 修正落地。
