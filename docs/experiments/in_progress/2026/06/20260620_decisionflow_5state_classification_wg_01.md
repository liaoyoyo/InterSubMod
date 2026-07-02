---
title: HCC1395 全基因組「位點甲基能否切成不同群」判別流程 — 5 態分類結果
date: 2026-06-20
status: in_progress
tier: 3
sample: HCC1395 (單樣本; ⭐2-3 偵測非驗證)
scope: 全基因組 (22 autosome, TP SNV, w=±5000, BERNOULLI, paired_full)
read_sets: tumor-only (subclone 向) + merged tumor+normal (cis-control 對照) 分開計算
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/decisionflow_summary.json
build_commit: e57a7ab
binary: build/bin/inter_sub_mod (e57a7ab, --distance-metric BERNOULLI --nan-distance-strategy SKIP)
observation_standard: true
---

# HCC1395 判別流程 5-態分類（全基因組）— 觀察標準

> **這份回答三件事**：① valid 群 size≥多少合適（門檻依據）· ② 依判別流程，HCC1395 每位點落哪一態、比例多少（tumor/merged 分開）· ③ 流程是否有效。
> **數字 SoT**：`decisionflow_summary.json`（由 `decisionflow_analyze.py` 從 `decisionflow_records_{tumor,merged}.json` 算出；records 由 `decisionflow_wg.py` 逐 chr 重跑 binary 產生）。
> **HTML 視覺版**：`InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/20260620_decisionflow_5state_classification_01.standalone.html`

## 0. 判別流程（5 態定義）

每個位點（read×CpG 甲基矩陣 + read×read BERNOULLI 距離）依序判：

1. **precondition**：去 NaN（peel）後 `n_complete≥6` 可聚類？否 → **① 不可驗證**。
2. **outlier-tolerant 無監督切群**（UPGMA average，T=valid 群門檻）：是否 ≥2 個 valid 群（size≥T）、outlier(<T)≤20%？
   - 切得出 → 是否對齊 HP/carrier/allele（CramérV≥0.3 & χ² p<0.05 & Cochran e≥5）？
     - 對齊 → **⑤ 確認真結構**
     - 不對齊 → **④ 可分未對齊**（epiallele? 雜訊?）
   - 切不出 → a-priori PERMANOVA（distance-based pseudo-F，permutation，BH-FDR per 軸）是否顯著？
     - 顯著 → **③ 監督可分 mean-shift**（內以 PERMDISP `dispP≥0.05` 分 location-clean / dispersion-confounded）
     - 不顯著 → **② 1群/無訊號**

headline 門檻 **valid≥4（≤3 當 outlier）**；T=3/4/5 敏感度一併報。NPERM=99。

## 1. 5-態組成與比例（valid≥4；N=30490）

| 狀態 | tumor 位點 | tumor % | merged 位點 | merged % |
|---|---:|---:|---:|---:|
| ① 不可驗證 | 413 | 1.35% | 377 | 1.24% |
| ② 1群/無訊號 | 5,051 | 16.57% | 715 | 2.35% |
| ③ 監督可分 mean-shift | 5,525 | 18.12% | 2,061 | 6.76% |
| &nbsp;&nbsp;└ location-clean / dispersion | 4,907 / 618 | | 1,836 / 225 | |
| ④ 可分未對齊 | 13,390 | 43.92% | 16,748 | 54.93% |
| ⑤ 確認真結構 | 6,111 | 20.04% | 10,589 | 34.73% |
| precond-pass | 30,077 | 98.65% | 30,113 | 98.76% |

## 2. valid 群門檻 size≥多少合適（數據依據）

**核心證據：對齊切群的最小群（minority）遠大於未對齊切群。**

| read-set | 對齊(→⑤) minority 中位 | 未對齊(→④) minority 中位 |
|---|---:|---:|
| tumor | **18** | **7** |
| merged | **22** | **9** |

→ 小 minority（3-5）的切群幾乎都「未對齊」（雜訊/單一離群鏈）；大 minority 才對齊真生物軸。

**門檻敏感度（T=3/4/5；⑤對齊 / ④未對齊 / 對齊率）**

| 門檻 | tumor ⑤ | tumor ④ | tumor 對齊率 | merged ⑤ | merged ④ | merged 對齊率 |
|---|---:|---:|---:|---:|---:|---:|
| valid≥3 | 5,590 | 17,528 | 24.18% | 10,852 | 17,325 | 38.51% |
| **valid≥4 ★** | 6,111 | 13,390 | **31.34%** | 10,589 | 16,748 | **38.74%** |
| valid≥5 | 6,324 | 11,244 | 36.00% | 10,400 | 16,465 | 38.71% |

→ tumor：門檻 3→5 把 ④（未對齊）從 17,528 砍到 11,244，對齊率 24%→36%；證明小群多為雜訊。**valid≥4** 是平衡點（保住 ⑤ 數量、砍掉大半小群雜訊）。先前個案測試另證：**群大小（minority≥6 vs =3）比 silhouette 區分度更乾淨地分開「該救/不該救」**，故以 size 為主判準。merged 因 germline 結構強，門檻不敏感（28177→26865）。

## 3. 「切不出 ≠ 沒訊號」

切不出離散群的位點，再測 a-priori mean-shift：

| read-set | ②真無訊號（=100−③） | ③ mean-shift（總） | ③ location-clean |
|---|---:|---:|---:|
| tumor | 47.76%（100−52.24） | **52.24%** | 46.4% |
| merged | 25.76%（100−74.24） | **74.24%** | 66.14% |

→ tumor 切不出的位點仍有近半（52%，location-clean 46%）有 a-priori 甲基均值差（PERMDISP 已濾散開度混淆）。**確認「切不出 ≠ 沒訊號」在 tumor 也成立**；merged 更高（74%）因 germline cis 結構普遍。

## 4. 判別流程是否有效

| 指標 | tumor | merged | 解讀 |
|---|---:|---:|---|
| 切群對齊率 ⑤/(④+⑤) | 31.34% | 38.74% | 無監督切得出的群，~1/3 對應已知生物軸 |
| 切不出有訊號 ③/(②+③) | 52.24% | 74.24% | 「切不出≠沒訊號」量化 |

**有效，但須讀對軸**：
- 流程**有效地把 30,490 位點分成可操作的 5 桶**，且 partial(chr1-10)↔全量高度一致（tumor ⑤ 19.5%→20.04%），跨基因組穩定。
- 🔴 **merged=cis-control**：merged 加 normal read → ② 17%→2%、⑤ 20%→35%。這是預期：normal 補強 germline-haplotype（cis-ASM）結構 → 更多位點切得出且對齊 **HP/allele 軸 = germline cis 效應，非 somatic subclone**。**tumor-only 才是 subclone 向**。
- 🔴 **④「可分未對齊」是最大桶**（tumor 43.9% / merged 54.9%）：含 (a) 小 minority + 無 perm 訊號 = 偏雜訊；(b) perm 顯著但離散切法不對 a-priori = 可能真 epiallele 子結構或更細層。門檻敏感度顯示 ④ 隨 T 上升大幅收縮（tumor −35%），故 ④ 相當部分是門檻邊緣的弱切群。
- 🔴 **單樣本 ⑤≠subclone**：⑤ 對齊 HP/carrier/allele 多為 germline cis-ASM；somatic subclone 需 normal cis-control 排除 + 多樣本，本表 ⭐2-3 偵測層。

## 5. 驗證表（每數字 → 來源 + 重算 + L 級）

| 數字 | 值 | 來源 key | 重算/加總 | L |
|---|---|---|---|---|
| tumor N | 30,490 | summary.tumor.n_total | =records_tumor.json 筆數；與 memory 切群總帳 30490 一致 ✓ | L1 |
| tumor 5 態加總 | 413+5051+5525+13390+6111 | states_T4 | =30,490 = N ✓ | L1 |
| tumor ③ loc+disp | 4907+618 | S3_loc+S3_disp | =5,525 = S3 ✓ | L1 |
| tumor 切群對齊率 | 31.34% | split_align_rate | 6111/(13390+6111)=6111/19501=0.3134 ✓ | L1 |
| tumor 切不出有訊號 | 52.24% | nonsplit_meanshift_rate | 5525/(5051+5525)=5525/10576=0.5224 ✓ | L1 |
| tumor minority 對齊/未對齊中位 | 18 / 7 | minority_dist.aligned/unaligned.median | records byT[4].minority 分組中位 | L1 |
| tumor 門檻 T3/T5 ④ | 17528 / 11244 | threshold_sensitivity | classify(T) 重算 | L1 |
| merged N | 30,490 | summary.merged.n_total | =records_merged.json 筆數 ✓ | L1 |
| merged 5 態加總 | 377+715+2061+16748+10589 | states_T4 | =30,490 = N ✓ | L1 |
| merged 切群對齊率 | 38.74% | split_align_rate | 10589/(16748+10589)=10589/27337=0.3874 ✓ | L1 |
| merged 切不出有訊號 | 74.24% | nonsplit_meanshift_rate | 2061/(715+2061)=2061/2776=0.7424 ✓ | L1 |

全數字 L1（從 records 機械重算）。PERMANOVA/PERMDISP F 值於 smoke test 向量化前後完全一致（22816187 allele F=38.455 不變）。

### 5b. 對抗驗證（2-agent workflow `wf_e3d0bc78`，2026-06-20）

| 軸 | verdict | 重點 |
|---|---|---|
| 獨立重算 | **VERIFIED_MATCH** | fresh agent 自寫 BH-FDR、未複用 analyze，從原始 records 重算 tumor 5 態全數一致（413/5051/5525/13390/6111）；partition 互斥窮盡 413+19501(split)+10576(nonsplit)=30490 ✓ |
| 方法學 | **SOUND_WITH_CAVEATS**（0 blocking） | state③ **非** tumor-only NEGATIVE 陷阱；PERMANOVA/PERMDISP 實作數值正確（向量化 2 群路徑與暴力 SSW bit-identical；betadisper 對非歐 BERNOULLI 的 real/imag 處理正確）；valid≥4 / FDR pool / merged=cis-control 皆成立 |

> **🔑 state③ 為何不是「tumor-only 非監督+PERMANOVA NEGATIVE」陷阱**（[[project_tumor_only_axis_negative_subclone_classification]]）：
> 那個 NEGATIVE 是 **cluster-first**（silhouette 先挑距離矩陣最分離的 partition → 再 PERMANOVA 測它 = **循環**，被 read-內甲基相關打敗，Noise→TI 100%）。
> 本 state③ 是結構上**相反**的 **label-first**：只對「**切不出**」的位點、用**外部 BAM tag**（HP/carrier/allele，**非**甲基 clustering 產生）做 **label-permutation PERMANOVA（距離矩陣固定）**。label 置換在 H0 下**保留 read-內相關結構** → exchangeability 成立 → p 值**非** anti-conservative。故 **52% 是真的 omnibus 標籤↔距離關聯，非置換假陽性**。
> ⚠ **但「真關聯」≠ subclone**：單樣本下主導軸是 allele/germline-haplotype = **germline cis-ASM / read-identity**，非 somatic subclone。本表 ⑤/③ 維持 ⭐2-3 偵測層、不宣稱 subclone（既有非循環判別率僅個位數 %，見 [[project_tumor_only_axis_negative_subclone_classification]]）。

## 6. 限制（含對抗驗證提出的 non-blocking concerns）

- **單樣本 HCC1395**：⭐2-3 偵測層；⑤/③ 對齊 a-priori 多為 germline cis-ASM，**非 somatic subclone**（需 normal cis-control + 多樣本）。tumor-only 的「subclone-向」僅是**方向**，非正向 subclone 宣稱。
- 🔴 **allele 軸是 3 軸中最 read-identity-confounded**（ALT/REF 標籤鄰近 somatic SNV，甲基 co-segregate 可能來自 cis-haplotype / mapping-bias 非 subclone）→ **next step = normal-anchored cis-control 專測 allele 軸**（對齊 FN-audit §2.4「allele 軸顯著需 cis-control」）。
- **FDR family 門檻依賴**：BH 套在「該 T 下非 split + precond-pass」位點池；位點會隨 T 進出 pool（非固定 family）。內部自洽（每 T 重 classify+重 pool）且已報，但 reviewer 須知非固定家庭。
- **NPERM=99 → p 下限 0.01**（51% hp p 落在 floor）：對 FDR-pooled discovery 足夠、52% 作為比例 robust，但**單位點 q 粗粒**，引用單點 p 須併報 F；增至 999 可銳化尾部（成本權衡非正確性）。
- **valid≥4 / dispP≥0.05 是約定門檻**；報告已標。dispP-clean 僅表「非散開度混淆」，**不蘊含 location→subclone**。
- **④ 桶異質**：含雜訊 + 可能真 epiallele，本表未進一步細分（門檻敏感度提供上界）。
- **merged HP 軸 = germline haplotype 層**（cis-ASM），其高 ⑤/③ 為 cis-control 預期，非 tumor subclone 證據。

## 7. 重生指令

```bash
cd InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot
python3 scripts/decisionflow_wg.py        # ~31min: binary→records (dual-mode)
python3 scripts/decisionflow_analyze.py   # records→decisionflow_summary.json
python3 scripts/decisionflow_plot.py      # →figs_decisionflow/fig1-4.png
python3 scripts/build_decisionflow_report.py  # →standalone HTML (§13-A 注入)
```
