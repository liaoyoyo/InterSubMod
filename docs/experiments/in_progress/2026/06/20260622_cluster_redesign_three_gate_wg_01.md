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

## 1. 解法 = 三閘「切更細」標準
切更細到 k 必須三閘齊備：
1. **balance**：群 ≥ **3 read**（業界 caller 最低支持標準，「多數足以確認非雜訊」）。**不用大小比例**。
2. **null-excess（「明顯」）**：該切法 clusterboot Jaccard **扣 within-1-group null（逐 CpG 欄內非 NaN 重排→重算 BERNOULLI）≥ 0.10**。防「過度碎裂 / 丟 read 湊對齊」。
3. **alignment（歸因）**：可靠對齊 germline 軸（CramérV≥0.3 & p<0.05 & Cochran e≥5）。
- **每群標**：`CONFIRMED`（對齊 germline）/ `REAL_NOVEL`（真實但不對齊 = **subclone 候選**）/ `edge`（≥3 離群成一組）/ `outlier`。
- **雙輸出**：`coarse`=最粗 CONFIRMED 骨幹 · `fine`=最細全真實結構（confirmed∪novel）+ 各自 per-read 群歸屬。

### 🔧 小規模驗證抓到的迭代（先小規模的價值）
1. 初版純 gap 切 → **過度碎裂**（chr4 切 27 群）。
2. 加 silhouette 閘 → 仍每位點列 [2..6]（silhouette 跨 k 平、選不出 k）。
3. 資料校準 → **「可靠對齊」乾淨選出有意義群數**（收斂）。
4. 加 balance 8% → 修 chr1「[77,3] 太鬆」但**過度促進**（chr4 k2→k5，丟 20% read 湊對齊）。
5. **整合 null-excess** → 修過切（chr4 退回 k=2）。
6. 用戶定 balance=3（業界）→「明顯」交給 null；雙輸出。

## 2. 結果（29 例，全 L1 from json）
| 指標 | 值 | 來源 key |
|---|---|---|
| fine confidence | CONFIRMED 12 / REAL_NOVEL 17 | loci[].fine_confidence |
| 有 REAL_NOVEL（subclone 候選） | 17 | has_REAL_NOVEL |
| 邊緣群位點 | 25 | has_edge_group |
| 多解析度 confirmed（≥2） | 10 | FM3_multiresolution_confirmed(≥2) |

### A/B/C 三類（對 subclone 目標）
| 類 | coarse/fine | n | 意義 |
|---|---|---:|---|
| **A** germline-cis | CONFIRMED=CONFIRMED（各 k 對齊） | 12 | cis-ASM，**非 subclone** |
| **B** 骨幹+候選 | CONFIRMED(粗)+REAL_NOVEL(細) | 7 | germline 骨幹 + **subclone 候選子結構** |
| **C** 純 novel | NONE(粗)+REAL_NOVEL(細) | 10 | 不對齊 germline 的真實結構 = **subclone 候選 / 需 cis-control** |

🔑 **B+C = 17 有 subclone 候選**：「真實但不跟 germline 走」正是 subclone 藏身處（真 subclone 甲基不被 germline haplotype 決定）。

## 3. 精細度確認標準（用戶問「完整/快速」）
**有** = 三閘精細度階梯：每位點掃 k=2..6 的 excess(扣 null)曲線 → **天花板 = excess 掉到 0.10 以下、或出現 <3 單離群那一階**。完整（掃所有 k）+ 快速（僅 null 成本秒級）+ 客觀（門檻明寫）。「小群沒切」= 看卡哪閘：excess<0.10→不夠明顯;群<3→歸邊緣;不對齊→REAL_NOVEL 記錄。圖：`figs_redesign/granularity_curve.png`。

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
