<!--
建立時間: 2026-07-18 06:10 +08:00
目標: 實作並驗證 HCC1395 k>8 positional components 的 read-support bounded segmentation
處理範圍: HCC1395 chr1-22；LongPhase-S recalibrated FILTER=PASS；既有 50 kb positional components
topic: 20260718_k_gt8_read_supported_segmentation
status: verdict_PROBE
audit_version: 0.1
關聯檔案:
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/pre-decision-audit.md
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/extract_lossless_read_linkage.py
  - InterSubMod/research/20260716_read_linked_hypercube_exact_likelihood_validation/scripts/lossless_read_contract.py
  - InterSubMod/docs/CURRENT_FOCUS.md
-->

# Pre-Decision Audit：HCC1395 k>8 read-support 分塊

> **Verdict：PROBE（70/100；red-team 由 GO 降一級）**。ordered hypergraph dynamic programming（DP）可精確求得固定 `k≤8` 非重疊分塊下的最大完整 read-support；但 full HCC1395 的 read extraction runtime、coverage bias、單 read 固定判讀超過 8 點的不可保存比例仍未實測。先過 toy、chr22、密集區與極端鏈 probe，才啟動 chr1-22。

> **服務目標**：G4（同一樣本完整可重現）與 G5（可外部審核的方法、數量與 runtime receipt）。

## §0 Task Type、scope 與研究啟動五問

- **Task Type**：B — Comprehensive validation。
- **固定 final scope**：HCC1395 chr1-22；probe 不取代 final。
- **Domain**：Complex；以 safe-to-fail probe 校準再 full。
- **比較基準**：canonical densest-8 舊流程；408 個 `k>8` components、47,570 sites、44,306 cap-excluded sites。

1. **Thread D？** 否；本輪使用 read-level genotype/linkage，不使用 methylation。
2. **Thread B 撤回範圍？** 否；不建立或宣稱 production FP filter。
3. **KDE-corrected？** 不適用；切割目標不使用 KDE。
4. **VCF caller AF？** 不使用；VAF 只允許在後續候選樹排序，禁止作切割權重或與相同 reads 雙計數。
5. **長計算／C++／搬移／NO-GO？** 是長 BAM 計算與新 spec；使用獨立 Python research cycle，不改 canonical、不刪搬檔，full run 受 resource gate 控制。

## §1 問題與主張上限

要解決的工程問題是：對一個按座標排序、因相鄰距離 `≤50 kb` 串聯成的 `k>8` component，將所有位點分配到連續、非重疊且每塊 `k≤8` 的 blocks，同時最大化「完整落入某一 block」的 read/pattern constraint 權重。

對 read/pattern group \(g\)，令 \(I_g\) 為固定可判讀 `R/A` 位點集合，權重 \(w_g\) 為支持 molecule 數：

\[
L(P)=\sum_g w_g\cdot I\left(\nexists B\in P:I_g\subseteq B\right)
\]

這是 ordered hypergraph partition，不是 pairwise clique-expanded min-cut。相同 constraint 即使跨多個 boundaries，也只損失一次。

**主張上限**：

- 可以稱「在固定 `k≤8` 工程限制下，read-support 最佳的局部非重疊分塊」。
- 不可稱「唯一真實演化樹」。
- `PARTIAL_CUT`、`UNAVOIDABLE_READ_GT8`、`INSUFFICIENT_BRIDGE_OBSERVATION` 只能進局部子結構分母，不能進 complete-tree 分母。
- VAF 可在分塊之後比較候選樹相對可能性，但結果仍是推測，不能回頭證明切割為 biological truth。

## §2 Credibility Score

| 維度 | 分數 | 證據 |
|---|---:|---|
| 理論基礎 | 20 | ordered hypergraph block-reward DP 可形式化且避免多 boundary double count |
| 觀察支撐 | 20 | canonical HCC1395 component/site census、frozen M2 extraction receipt 均已存在 |
| 機制清晰度 | 20 | BAM → unique molecule sparse calls → constraints → DP → blocks → downstream tree 分層清楚 |
| 反例風險 | 10 | 已預註冊跨多切點、HP/PS 混合、single-read >8、coverage weighting 反例 |
| 所需資源 | 0 | full HCC1395 尚未實測；當前 load/CPU 超過 production gate |
| **TOTAL** | **70 / 100** | 分數達 GO，但 red-team 因 runtime 與不可保存比例未知降為 **PROBE** |

## §3 假設與 falsifiers

| 假設 | 重要性 | 現況 | 動作 |
|---|---:|---:|---|
| primary unique molecule sparse calls 足以表示切割 constraint | High | 已有 frozen extractor receipt | probe 重算 conservation |
| HP family × exact known PS 不可混合 | High | 已知成立 | schema 與 verifier hard gate |
| 一條 fixed R/A read >8 可安全截短 | High | False | 標 `UNAVOIDABLE_READ_GT8`，不可稱 lossless |
| absence of crossing read 等於 biological independence | High | False | 另算 span/callable opportunity；不足則 abstain |
| pairwise boundary cost 可直接相加 | High | False | 主算法用 block reward，每 constraint loss 一次 |
| VAF 可決定切點 | High | False | 禁止；VAF 只作 downstream candidate-tree ranking |

若出現任何一項，方法維持 PROBE 或 NO-GO：

1. toy oracle 與 brute force optimum 不一致；
2. 打亂 constraint 輸入後 digest 改變；
3. site、constraint 或 HP/PS conservation 失敗；
4. 任一輸出 block `k>8`；
5. 同一 constraint 的 loss 被計算超過一次；
6. `k≤8` 原 components 發生非預期變化；
7. full runtime 超過舊 baseline 4 倍且無可行 sparse/cache 修正。

## §4 Step → Verify

1. **實作純 DP core**  
   → 驗證：toy、brute-force、shuffle、HP/PS isolation、single-read >8 測試全部 exit 0。
2. **實作 chromosome-streaming I/O 與 receipt**  
   → 驗證：輸入 hash、stage time、output hash、site/constraint conservation 均寫入 JSON。
3. **chr21 negative + chr22 E2E probe**  
   → 驗證：chr21 `SKIP_NO_TARGET`；chr22 6 個 `k>8` components、98 sites 全部有 fate。
4. **chr10 dense + chr6 extreme probe**  
   → 驗證：`k=12` 與最大 `k=3574` 均完成；peak RSS、wall time 與不可保存 constraint 比例可量測。
5. **HCC1395 chr1-22 full**  
   → 驗證：408/408 components、47,570/47,570 sites、每 block `k≤8`、`MAX_SNV` exclusion=0。
6. **與 densest-8 比較**  
   → 驗證：site retention、constraint retention、status、blocks、wall/CPU/RSS/I/O 均有相同 denominator 與可追溯 receipt。
7. **HTML 技術報告**  
   → 驗證：artifact validator、portable HTML packaging 與 browser/structural receipt 通過。

## §5 預註冊 gates

- **Correctness hard gates**：site conservation=100%；constraint mass conservation=100%；每 block `k≤8`；HP/PS 零混合；rerun digest 一致。
- **實用性 GO**：至少 80% 舊排除位點進入 `NO_LOSS` 或 read-support retention≥95% 且 weight-stable 的 partial blocks。
- **PROBE**：50–80%。
- **不可取代 densest-8**：<50%。
- **Runtime**：新 end-to-end `≤2×` 舊約 84.78 分為可接受；`2–4×` 需優化；`>4×` 工程 NO-GO。

## §6 Resource gate

2026-07-18 06:10 +08:00 實測：

- 48 logical CPUs；MemAvailable 約 359 GiB。
- `/big7_disk` 約 716 GiB free，但使用率 99%。
- load1 約 79.7，`load/CPU≈1.66`，高於 full-run gate `1.25`（load 60）。

因此允許低負載程式開發與小 probe；full HCC1395 必須等 load≤60 且確認無同類 BAM/NFS 高 I/O 工作後才啟動。

## §7 Decision

- **允許**：獨立 research scripts/tests、既有 frozen read extractor 的唯讀重用、probe 與 receipt。
- **暫停**：在當前 contention 下啟動 full；直接改 canonical；把 partial block 稱完整真實樹；用 VAF 或 methylation 決定切點。
- **最簡替代方案**：維持 densest-8 可節省 runtime，但會繼續排除 44,306 sites，不能回答使用者要求的全位點 fate，因此只保留作 baseline。

