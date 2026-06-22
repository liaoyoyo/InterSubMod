---
title: 切群重設計 — 三閘(balance+null-excess+alignment) + coarse/fine 雙輸出 + 邊緣群 + per-read
date: 2026-06-22
status: in_progress
tier: 3
sample: HCC1395 tumor-only（單樣本 ⭐2-3 characterization）
scope: kprofile_examples 29 代表位點（小規模驗證；全基因組 sweep 未跑）
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/cluster_redesign.json
build_commit: 578763c
observation_standard: true
---

# 切群重設計 — 三閘 + coarse/fine（修 FM1/FM2/FM3 + 精細度標準）

> **HTML**：`InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/20260622_cluster_redesign_observation_01.standalone.html`
> **動機**：用戶觀察舊切群 3 缺陷 — (FM1) 單離群吃群位→整位點當「無法切」；(FM2) maxclust 強切碎/距離門檻才對；(FM3) 多解析度被折成標籤數。

## 0. 三個失敗模式（實測確認，29 例）
| FM | 現象 | 實測 |
|---|---|---|
| FM1 單離群吃群位 | 1-2 read 離群占 maxclust 群名額 → min<3 整 k 被丟 | 8/29（chr4 4 大群只接受到 k=2） |
| FM2 maxclust 強切 | 強制 k 群會切碎緊密群 | 2/29（且純 gap 自身過切→需配 FM1 隔離） |
| FM3 多解析度折疊 | k=2/k=3 都清楚卻只輸出一個 | 26/29 有多有效 k |

## 1. 解法 = 四閘「切更細」標準
切更細到 k 必須四閘齊備：
1. **balance**：群 ≥ **3 read**（業界 caller 最低支持標準，「多數足以確認非雜訊」）。**不用大小比例**。
2. **null-excess（「明顯/真實」）**：clusterboot Jaccard **扣 within-1-group null（逐 CpG 欄內非 NaN 重排→重算 BERNOULLI）≥ 0.10**。防丟 read 湊對齊。⚠ excess 在高覆蓋位點隨 k **自動上升**、不能單獨定天花板。
3. **gap-significance（精細度天花板，scale-invariant）**：分支跳躍 ≥ **max(8×中位 gap, 0.4×該樹最大 gap)**。絕對地板確保樹有結構；相對 max 防大樹（n 大→中位 gap 極小）過切。
4. **alignment（歸因）**：可靠對齊 germline 軸（CramérV≥0.3 & p<0.05 & Cochran e≥5）。
- **supported = real & (big_gap 或 reliable_align)**；精細度只納 supported。
- **4 類**：`CONFIRMED`（真實+對齊）/ `REAL_NOVEL`（真實+大跳不對齊=**subclone 候選**）/ `REAL_DIFFUSE`（真實但無大跳+無可靠對齊=低信心候選，**不丟棄**）/ `NO_CLEAR`（全無真實=**真無法切**）/ edge（≥3 離群）/ outlier。
- **雙輸出**：`coarse`=最粗 CONFIRMED 骨幹 · `fine`=最細 supported + 各自 per-read 群歸屬。

### 🔧 小規模驗證抓到的迭代 + 2 個 bug（先小規模的價值）
1. 初版純 gap 切過度碎裂(chr4 27群) → 加 silhouette 仍列[2..6](選不出 k) → 校準發現「可靠對齊」乾淨選群數。
2. balance 8% 修 chr1[77,3]太鬆但過度促進 → **整合 null-excess** 修過切 → 用戶定 balance=3、「明顯」交給 null。
3. 🐛 **bug1 stab_excess 用 maxclust(單離群)→ excess=None 誤判 NOT_REAL**（用戶 chr4 ALT 子群觀察揪出）→ 改評估實際核心群 → chr4 ALT 群(k=5 gap 138×)浮出，但 excess 隨 k 上升使 fine 過衝 k=6。
4. 🐛 **bug2 NO_CLEAR 把「真實但 diffuse」誤判「無結構」**（用戶質疑「無法切是真無法還是方法問題」揪出）→ 4 個原 NO_CLEAR 全有真實結構(excess 0.3-1.0)→ 加 **REAL_DIFFUSE** 類。
5. **gap-significance scale-invariant**（max(8×中位, 0.4×最大)）定天花板 → chr4 拉回 fine k=3(coarse2+ALT novel)。

## 2. 結果（29 例，全 L1 from json；gap scale-invariant 收斂後）
| 指標 | 值 | 來源 key |
|---|---|---|
| fine confidence | CONFIRMED 14 / REAL_NOVEL 11 / REAL_DIFFUSE 4 / NO_CLEAR 0 | loci[].fine_confidence |
| 邊緣群位點 | 25 | has_edge_group |
| 多解析度 confirmed（≥2） | 10 | FM3_multiresolution_confirmed(≥2) |

### A/B/C/D 四類（對 subclone 目標；HTML cls_of）
| 類 | coarse/fine | n | 意義 |
|---|---|---:|---|
| **A** germline-cis | CONFIRMED=CONFIRMED（各 k 對齊） | 13 | cis-ASM，**非 subclone** |
| **B** 骨幹+候選 | CONFIRMED(粗)+novel/diffuse(細) | 4 | germline 骨幹 + **subclone 候選子結構** |
| **C** 純 novel/diffuse | NONE(粗)+REAL_NOVEL/DIFFUSE(細) | 12 | 不對齊 germline 的真實結構 = **subclone 候選 / 需 cis-control** |
| **D** 真無法切 | fine NO_CLEAR(全無真實) | 0 | 29 例皆預選 splittable → 真 NO_CLEAR 要全基因組才出現 |

🔑 **B+C = 16 有 subclone 候選**：「真實但不跟 germline 走」正是 subclone 藏身處（真 subclone 甲基不被 germline haplotype 決定）。

## 3. 精細度確認標準（用戶問「完整/快速」）
**有** = gap-significance 階梯：每位點掃 k=2..6，**天花板 = 最細「big_gap(scale-invariant) 且 real(excess≥0.10)」或可靠對齊那一階**。完整（掃所有 k）+ 快速（僅 null 成本秒級）+ 客觀（門檻明寫）+ **不會隨 k 無限上升**（excess 會，gap-ratio 不會）。「小群沒切」= 看卡哪閘：excess<0.10→真無結構；real 但 gap 弱+不對齊→REAL_DIFFUSE(記錄非丟棄)；群<3→邊緣。

### 🔑 「無法切是真無法還是方法問題」驗證（用戶問）
4 個原判 NO_CLEAR 的位點（chr1_119446237/chr3_1026519/chr2_44798355/chr3_131113525）**沒有一個是真無法切** — 全有真實結構（excess 0.3-1.0、熱圖可見多塊、部分 allele 對齊），只是 gap 弱 + 單樣本 e<5 無法歸因 → **是方法太嚴的誤判，已加 REAL_DIFFUSE 修正**。UPGMA 樹圖 `figs_dendro/dendro_*.png` 肉眼可證。流程圖 `figs_dendro/method_flowchart.png`。

## 4. 限制
- 代表位點皆已 splittable → **FM1 的 CANT_SPLIT→可切 rescue 要全基因組才看得到**。
- 單樣本 ⭐2-3：**REAL_NOVEL = subclone 候選非確認**，需 normal cis-control 排 cis-ASM。
- granularity_diagnostic.py 用 maxclust（單離群干擾），正式版應先 quarantine 再評估小群。
- 純分析、未改 C++（若有效再走 methodology-audit→cpp-change）。

## 5. 重生
```bash
cd InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot
python3 scripts/cluster_redesign.py            # 三閘 + coarse/fine + per-read → cluster_redesign.json
python3 scripts/plot_cluster_redesign.py       # old vs new 對照圖
python3 scripts/granularity_diagnostic.py      # 精細度階梯曲線
python3 scripts/build_redesign_observation.py  # 觀察 HTML（§13-A 注入）
```

## 6. 下一步
- **全基因組**：改指 **big7 本機資料**（已有完整副本，big8=NFS 慢）重跑 binary → 套三閘切群 → 量化全域 A/B/C 比例 + old vs new delta + REAL_NOVEL(subclone 候選)位點數。
- 關聯：[[project_subcluster_cluster_count_determination]]、[[project_distance_matrix_cluster_validation_methods]]、[[project_tumor_only_axis_negative_subclone_classification]]、[[feedback_baseline_dependence_not_result]]。
