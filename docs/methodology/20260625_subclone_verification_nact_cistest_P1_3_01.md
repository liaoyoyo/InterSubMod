<!--
建立時間: 2026-06-25
類型: methodology — P1.3 Normal-Anchored Concordant cis-Test (NACT) + 對抗驗證 + genotype census
狀態: in_progress（單樣本 HCC1395 ⭐2-3 觀察層；mostly_artifact 裁決經 4-agent 對抗驗證）
build_branch: docs/method-comparison-ism-external-202606
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/nact_summary.json, docs/methodology/_assets/20260618_subcluster_pilot/nact_stats.json, docs/methodology/_assets/20260618_subcluster_pilot/survivor_census.json, docs/methodology/_assets/20260618_subcluster_pilot/anchored_retest.json, docs/methodology/_assets/20260618_subcluster_pilot/candidate_characterization.json
-->

# P1.3 — 用 normal-anchored cis-test 洗 subclone 候選（HCC1395，雙軌；mostly_artifact 裁決）

> **TL;DR**：對 1,139 個 methylation subclone 候選跑 normal-anchored cis-test（NACT）。**DEMOTE 側可信**（345 個正確判 cis-ASM）。但 **102 個「存活候選」經 4-agent 對抗驗證 = mostly_artifact，根因是 double-dip 循環（非 CN/LOH）**。決定性 genotype census + genotype-軸 re-test 把可救的收斂到 **9 個非循環、genotype-錨定、tumor-specific somatic-ASM 訊號**（仍 = characterization 非 confirmed subclone）。**深層結論：單樣本甲基無法非循環地證明 subclone**；真確認只能靠 Phase 2 sSNV read-level 連鎖。⭐2-3 單樣本。**數字全來自上列 JSON（本輪 Read 回真值）。**

## §0 NACT 方法（Normal-Anchored Concordant cis-Test）

設計經獨立 workflow（3 透鏡 + 對抗合成，spec=`_assets/.../scripts/` + 工作流紀錄）。對每候選位點：

- **Track A（純 tumor reads，去 T/N 組成 confound）**：結構軸優先序 cluster_split → 平衡 a-priori somatic 軸 fallback → 否則 STRUCTURE_TOO_THIN。per-CpG signature `S_struct`（MWU+BH-FDR<0.05、|Δβ|≥0.2、coherent run≥3），否則 NO_TUMOR_SIGNATURE。
- **Track B（純 normal reads，獨立 germline baseline，不碰 clustering）**：在 S_struct CpG 上三讀數 — **R1** same-signature 投影 germline HP 軸（+ containment）/ **R2** per-CpG residual 相減 / **R3** phase-free normal-bimodality。
- **AND-gate verdict**：任一讀數喊 cis → **cis_asm**（demote-on-any）；全不喊 cis 且 R2 residual 存活 且 R3 uniform → **candidate_subclone**；任一覆蓋失敗/不一致 → **undetermined**（寧 undetermined 不誤判）。
- CN/LOH guard（候選 96.7% 落 CN/LOH）+ STEP 9 permutation null（設計含、**首跑未實作**，見 §6 caveat）。

## §1 主結果（全 1,139 候選，nact_summary.json）

| verdict | n | % | 意義 |
|---|---|---|---|
| **cis_asm**（任一讀數 demote） | 345 | 30.3% | normal 沿同 signature 共分離 → germline-cis |
| **undetermined** | 561 | 49.3% | 405 no-consensus + 92 no-baseline + 52 thin + 12 sparse |
| **no_tumor_signature** | 131 | 11.5% | tumor-only 無 coherent signature = T/N 假結構 |
| **candidate_subclone** | 102 | 9.0% | AND-gate 存活 → **下節證實 mostly artifact** |

by set：TP candidate 91 / FP 11。

## §2 對抗驗證裁決（4-agent workflow，verdict = mostly_artifact）

獨立 review（3 透鏡各 spot-check 真 region + 合成）裁決 **102 mostly_artifact，根因 double-dip 循環非 CN**：

- **DEMOTE 側 work**：345 cis_asm 的 germline_baseline 0.139 / containment 0.196 / frac_normal_bimodal 0.316 — 全顯著高於 survivors → 正確移除真 germline-cis。
- **PROMOTE 側 = double-dip**：100/102 用 cluster_split（甲基衍生軸，專案已知 tumor-only 無監督 clustering = NEGATIVE double-dip）。spot-check survivors chr5:88881646 / chr2:153468381 = **tumor reads 100% genotype 同質**（全 hp=2-1+ALT）卻被甲基切 ~2:1 → 無 a-priori 軸可錨 → 純 collider；R1/R2/R3 三讀數**結構性必然一致**（零判別力）。

## §3 決定性 census（survivor_census.json — §13.7 fresh 重算，獨立確認 verification）

**candidate vs cis_asm 對比**（證實「signature 閘不判別 + R1 結構性必然」）：

| 指標 | candidate(102) | cis_asm(345) |
|---|---|---|
| struct_dbeta_median（tumor signature 強度） | **0.718** | 0.712 ← **相同，不判別** |
| germline_baseline_strength | **0.036** | 0.139 ← survivors 近零 germline ASM |
| containment（與 normal germline 軸重疊） | **0.0** | 0.196 |
| direction_concordance | **0.0** | 0.333 |
| frac_normal_bimodal | 0.083 | 0.316 |

**genotype-homogeneity census**：102 = **homogeneous 72 / anchored 30**。dominant_geno_frac median **1.0**（63/102 完全同質，71/102 ≥0.9）→ 多數 survivors 是純甲基切群、無 genotype 基礎。

## §4 genotype-軸 re-test（anchored_retest.json — 把 30 anchored 改用 a-priori 軸非循環重測）

| 結果 | n | 意義 |
|---|---|---|
| NO_GENO_SIGNATURE | 15 | cluster split 與 genotype 軸無關 = double-dip（即使有軸） |
| genotype-軸 OK | 15 | → candidate 9 / undetermined 5 / cis_asm 1 |

**淨 9 個 candidate**（8 TP / 1 FP；7 REF_vs_ALT + 2 HP1_vs_HP1-1 軸）= genotype-錨定（非循環）、normal 不解釋（tumor-specific）的 somatic-allele 甲基訊號。🔴 但這 = **somatic-ASM（somatic 變異等位帶不同甲基）= characterization，非 confirmed subclone**（subclone 需 ALT 群內再有結構 → 又回 double-dip）。

## §5 完整誠實 funnel

```
34,736 全位點
└─ 1,139 subclone 候選 (cat8 A 723 ∪ subclone_novel 523)
   └─ NACT cis-test:
      ├─ cis_asm 345 ······· DEMOTE 正確（germline-cis 移除）
      ├─ undetermined 561 ··· 覆蓋/一致性不足
      ├─ no_tumor_signature 131 ··· T/N 假結構
      └─ candidate 102 ····· 經驗證 mostly artifact
         └─ genotype census: 72 double-dip(同質) + 30 anchored
            └─ genotype-軸 re-test: 15 無 sig + 9 tumor-specific somatic-ASM + 5 undet + 1 cis
               └─ 淨 9 個非循環 genotype-錨定 tumor-specific 訊號（characterization 級）
```

## §6 統計誠實（nact_stats.json）+ 機制更正

- **🔴 neutral=0/11 underpowered，不可用**：candidate rate neutral 0/11=0%（95% CI [0, 28.5]，涵蓋 LOH 率）；Fisher neutral vs LOH **p=0.615 不顯著**。先前「neutral=0 → CN-driven」推論**撤回**（powerless）。
- **LOH 特異升高**：LOH 10.3% vs gain/loss 5.8%，Fisher **p=0.016 顯著**（但 survivors germline_baseline 0.036 → 機制**非** LOH-unmask germline ASM；是 double-dip）。
- **candidate 不分 TP/FP**：TP 9.3% vs FP 6.7%，Fisher **p=0.373 不顯著** → 非 somatic-specific = artifact-consistent。
- **count 脆弱**：R3 門檻 102@<0.15 但 70@0.12 / 58@0.10 / 35@0.05 → 「102」threshold+axis-dependent，非穩健。

## §7 能宣稱 vs 不能宣稱（對抗驗證定）

**可宣稱**：① NACT DEMOTE 側 work（345 cis-ASM 正確移除）；② 102 survivors mostly double-dip artifact（72/102 genotype 同質、signature 與 cis_asm 相同 0.718/0.712）；③ 淨 9 個非循環 genotype-錨定 tumor-specific somatic-ASM 訊號 → Phase 2 候選；④ candidate = 「not refuted as germline-cis」⭐3，非偵測。

**🔴 不可宣稱**：① 「102 subclone」/任何偵測語；② 「整池 CN-driven / LOH-unmask」（germline_baseline 0.036 反駁）；③ neutral=0 證 CN-driven（p=0.615 powerless）；④ R1/R2/R3 為獨立佐證（double-dip loci 結構性必然一致）；⑤ 單樣本任何 survivor 升 subclone（⭐3 characterization-only）。

## §8 深層結論（回答 Q-B / Q-C）

**單樣本甲基無法非循環地證明 somatic subclone**：genotype-同質群內的甲基子結構**永遠 double-dip**（無 a-priori 錨）；唯一非循環的甲基訊號是 **genotype-aligned = cis/somatic-ASM（characterization）**，定義上非 subclone。→ 「結構×標籤異常 = subclone 候選」（Q-B）**必要非充分且多為循環**；subclone 真確認（Q-C）**只能靠 Phase 2 sSNV read-level 連鎖**（正交、非循環）。甲基的角色 = **characterize/corroborate 已被 sSNV 連鎖定義的 lineage**，非偵測 subclone。

## §9 下一步

- **Phase 2**：把 **9 個 tumor-specific（+ 30 anchored）** 交給 sSNV read-level 連鎖（feasibility 掃描 ≥2 sSNV 共現 → 互斥驗證）做真確認。
- **可選補強**（不改裁決）：STEP 9 改 held-out reproducibility null（非 relabel）；修 R1 direction-concordance bug（line 228，使更保守、不救 102）；negative-control census（甲基切群機器在同質 normal 上的 manufactured-signature FDR）。

## §10 Provenance + 自驗

| 數字群 | 來源檔（grep 可得） | 腳本 |
|---|---|---|
| NACT verdict 分佈 | `nact_summary.json` | `p1_nact_cistest.py` |
| CI + Fisher | `nact_stats.json` | `p1_nact_stats.py` |
| genotype census + 對比 | `survivor_census.json` | `p1_survivor_census.py` |
| genotype-軸 re-test | `anchored_retest.json` | `p1_anchored_retest.py` |
| 候選 profile | `candidate_characterization.json` | `p1_candidate_characterize.py` |

> 紅線：⭐3 單樣本 HCC1395；candidate=「not refuted as germline-cis」非偵測；DEMOTE 側可信、PROMOTE 側 mostly double-dip artifact；機制 = double-dip 循環非 CN（germline_baseline 0.036 反駁 LOH-unmask）；單樣本甲基無法非循環證 subclone → Phase 2 sSNV 連鎖。§13 數字不手打，全 JSON 溯源。
