<!--
建立時間: 2026-06-24
類型: methodology — 逐 CpG 差異甲基歸因 + subclone 驗證 + site-first vs read-level 分工
狀態: in_progress（單樣本 HCC1395 探索 ⭐2-3）
build_branch: feat/summary-nreadsvalid
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/percpg_summary.json, docs/methodology/_assets/20260618_subcluster_pilot/phylo_cpp_wg_full_records_v6.json, docs/methodology/_assets/20260618_subcluster_pilot/phylo_cpp_wg_full_percpg.json
-->

# 逐 CpG 差異甲基歸因 + subclone 驗證 — site-first vs read-level 分工（HCC1395 single-sample）

> **回答用戶 5 問**：① label_structure PERMANOVA 進對齊欄可行；② 每位點輸出哪些 CpG 對齊哪軸；③ 無法對齊但有 coherent CpG 的有效驗證；④ cluster-first 是否等同 site-first / site-first 是否足夠；⑤ 多切法是否都要驗。**⭐2-3 觀察層，非 subclone caller。**

## §0 核心方法學答覆（Q5/Q6）

| 問題 | 答覆（grounded `docs/method_comparison/20260609_ism_vs_external_methylation_tools/`）|
|---|---|
| **site-first（逐 CpG×標籤）是否足夠？** | **對「標籤已定義」問題（cis-ASM: HP/allele/sub-hap/T-N）足夠且更標準**（= DSS/modkit 範式）。read-level 在此**不加值**（間接）。本次實證：cis-ASM 例的 cluster-split 差異 CpG 與標籤 CpG **containment=1.0**（cluster 只是複述標籤）。|
| **read-level（cluster-first）唯一價值？** | **無監督發現 a-priori 標籤沒有的子群**。site-first 測不了沒定義的群。本次找到 **523 個 subclone_novel**（coherent CpG 不被任何標籤解釋）= site-first 抓不到。|
| **多切法（coarse/fine/geometry）都要驗？** | **不需各自獨立驗** — 它們是 k/stability 問題。驗主切（coarse 主分叉）+ 切法相依即 S3 unstable（已標）。verdict 對主分叉計算。|
| 🔴 **誠實邊界** | coherent CpG 只驗到「**甲基子群**」非「**somatic subclone**」。本次 523 subclone_novel **FP 富集 1.78×** → 非 somatic-specific（需 normal-control/single-cell 才可宣稱 subclone）。|

## §1 方法（Q2/Q3 — 逐 CpG 差異 + 歸因）

**site-first 逐 CpG**：每位點對各二元分組{HP1/HP2、HP1/HP1-1、HP2/HP2-1、REF/ALT、tumor/normal、**cluster-split 主分叉**}，逐 CpG 算 Δβ + **兩統計並列**：
- **Fisher 2×2**（meth>0.5/unmeth × 組）+ BH-FDR — 對齊既有 PerCpgAsm 口徑。
- **Mann-Whitney U**（連續 β 秩和）+ BH-FDR — **rank-based 對 read 非獨立穩健**。
- sig = FDR<0.05 且 |Δβ|≥0.2（判別用 MWU）。

> 🔴 **beta-binomial（DSS）在單樣本 read-level 不適用**：DSS 靠跨 replicate/全基因組借 dispersion，單樣本單位點**無 replicate**。pilot 證實「單位點共享 φ」會把組間真訊號當過離散→循環→殺光 power。改用 **MWU**（method_comparison 亦指 rank-based 繞過非獨立）。Fisher 標 sound-for-purpose（量級偏 anti-conservative）。

**🔑 判別器（cis-ASM vs subclone）= CpG-set overlap**：cluster-split 的 sig CpG 集合 vs 各 label sig CpG 集合的 **containment**：
- containment ≥0.5 → **cis-ASM(cluster≈該軸)**（cluster 只複述標籤）。
- 全 label containment <0.5 且 cluster coherent run ≥3 → **subclone_novel**（distinct 新訊號）。
- coherent run = 依基因座標連續且同向 |Δβ|≥0.2 顯著 CpG 串（DMR-like, DSS callDMR ≥2-5 標準）。

## §2 結果（全 34,736 位點 verdict 分佈，絕對比例）

```
no_cluster (coarse=1, 無切群可測) ········ 24,180 (69.61%)
有結構 coarse≥2 ·························· 10,556 (30.39%)
├─ cis-ASM(cluster≈標籤) ················ 7,098 (20.43%)  ← site-first 足夠
│    REF/ALT 3,804 · tumor/normal 1,406 · HP1/HP2 1,239 · HP1/HP1-1 339 · HP2/HP2-1 310
├─ ⭐ subclone_novel ···················· 523 ( 1.51%)  ← read-level 獨有, coherent CpG 不被標籤解釋
├─ cluster_no_sigCpG(切群但無 coherent) ·· 2,284 ( 6.57%)  ← 多為不平衡(小群灌 Δβ群), 逐 CpG 打臉
└─ weak_scattered ······················· 651 ( 1.87%)
驗算 7098+523+2284+651 = 10,556 ✓
```

## §3 🔴 523 subclone_novel 的誠實裁決

| 指標 | 值 | 意義 |
|---|---|---|
| 數量 | 523 (1.51%) | read-level + 逐 CpG 驗證雙過的「新甲基子群」|
| TP/FP | 410 / 113，**FP率/TP率 = 1.78** | **反判別 — FP 富集 → 非 somatic-specific** |
| CN | LOH 307 / gain 207 | 多在 LOH/CN 變動區（cis-residual/CN 鏡子嫌疑）|
| cat8 | C-S1 219 / C-S2 144 / A 107 / C-S3 53 | 散在各態 |
| 範例 | chr17:70314452 (run=35!) / chr8:71956115 (run=28) / chr19:4967388 (run=20) | 大 coherent run、皆 LOH |

**結論**：read-level 確實能**發現並驗證** site-first 抓不到的 novel 甲基子群（523 個，coherent run 達 35 CpG）；**但它們 FP 富集 1.78× → 在單樣本 tumor-only 不對應 somatic subclone**（與整體反判別一致）。**「甲基子群存在」≠「somatic subclone」** — 要宣稱後者需 normal-control 扣 cis / multi-region / single-cell。

## §4 label_structure 對齊欄（Q1，已落地）
significance.json 的 `label_structure`（hp/allele PERMANOVA p/F + dispersion 旗）全 34,736 有值，已併入 records_v6（`ls_hp_p/ls_hp_disp/ls_al_p/ls_al_disp`）。**HP PERMANOVA 顯著且無 dispersion warning = 8,360 位點**（乾淨 HP-aligned；其餘顯著者帶 dispersion 警告 = 散度驅動須存疑）。取代純 CramérV 對齊欄。

## §5 Pilot 驗證（7 範例，方法可信度）
- cis-ASM 例（①⑦）：cluster-split CpG 與正確軸 **containment=1.0**、coherent run 3-6 → 方法正確識別「cluster 複述標籤」。
- no-structure（⑤⑥ B1/D）：verdict=no_cluster ✓。
- 「subclone 候選」②④：split 不平衡(53/3, 225/3)→ MWU 無 power → cluster_no_sigCpG → **逐 CpG 比 Δβ群嚴謹**（Δβ群=0.396 被 3-read 小群灌出，per-CpG 打臉）。

## §6 Provenance + 自驗
| 產物 | 腳本 |
|---|---|
| 逐 CpG 核心(Fisher+MWU+overlap) | `scripts/per_cpg_diff.py` |
| 全量 verdict + label_structure | `scripts/extract_per_cpg_full.py` → `phylo_cpp_wg_full_percpg.json` / `percpg_summary.json` |
| 併入 records | `scripts/merge_percpg_v6.py` → `phylo_cpp_wg_full_records_v6.json` |

```bash
cd InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot
python3 -c "import json;print(json.load(open('percpg_summary.json'))['verdict_dist'])"   # verdict 分佈
python3 -c "import json;r=json.load(open('phylo_cpp_wg_full_records_v6.json'));sn=[x for x in r if x['pc_verdict']=='subclone_novel'];print('subclone_novel',len(sn),'TP/FP',sum(x['set']=='TP' for x in sn),sum(x['set']=='FP' for x in sn))"
```

## §7 差異甲基「定位點」量化 + FP/TP 歸因 + 研究可用性（2026-06-24）

全 34,736 位點 site-first 逐 CpG（MWU FDR<0.05 且 |Δβ|≥0.2）→ 共 **2,819,893 差異 CpG instance**（跨軸）。

### 各軸差異 CpG（數量 + 方向）
| 軸 | 有差異位點 | 差異 CpG 數 | hyper%（增甲基）|
|---|---|---|---|
| tumor/normal | 34,253 | 1,240,666 | **90.1%** |
| HP1/HP2 | 29,935 | 582,728 | 50.5% |
| HP1/HP1-1 | 14,940 | 312,168 | 85.3% |
| HP2/HP2-1 | 14,647 | 305,808 | 85.4% |
| REF/ALT | 19,727 | 214,056 | 46.7% |
| cluster-split | 8,272 | 164,467 | 42.9% |

🔑 **somatic 軸（T/N・HP-1）hyper 偏高（85-90%）= somatic 甲基多為「增甲基(gain)」**；germline/cluster 平衡(~50%)。

### 「定位點」覆蓋 + FP/TP 歸因 + 研究可用性（≥3 差異 CpG）
| 軸/狀況 | %位點 | FP率/TP率 | 機制 | 研究可用 |
|---|---|---|---|---|
| 任一軸 union | **96.8%** | ~1 | 普遍存在 | ⚠ 不判別 |
| tumor/normal | 95.8% | 0.98 | 普遍 somatic 偏移 | ⚠ 描述可、判別不可 |
| germline HP1/HP2 | 77.2% | 0.98 | cis-ASM 非 somatic | ⚠ |
| **somatic HP-1**（HP1/HP1-1 ∪ HP2/HP2-1）| **70.8%** | **0.33** | somatic 變異定義子單倍型 | ✅ **TP 專一, 可做研究** |
| REF/ALT | 43.3% | 0.50 | TP 有真 ALT；FP=germline/artifact | ✅ TP 專一 |
| cluster/subclone | 23.4% | 1.03 | 子群 cis/CN 驅動 | ✗ FP 富集, 單樣本不可宣稱 subclone |

### 結論（研究可用性裁決）
- **可做研究的差異定位點 = somatic 甲基 marker（HP1/HP1-1, 70.8%, TP 專一 0.33×）+ REF/ALT（43.3%, 0.50）** → 標記「哪些 read 帶 somatic 變異 + 其甲基」。
- 🔴 **但「TP 專一」未扣 normal-control** → 可能含 germline cis-residual（真 somatic 須 normal-anchored，COLO829 缺）= **gap#2**。
- ✗ **subclone 定位（cluster, 23.4%）FP 富集 → 單樣本不可當 somatic subclone 研究基礎**。
- ⚠ 「有差異甲基」幾乎全有（96.8%）但不判別 → **量本身無研究價值，TP 專一 + 方向(hyper) + 強度才有**。
- 已知 gap：① 逐 CpG **座標未存**（只數量/方向/強度）；② TP 專一未扣 normal-control；③ strand 未分（5mC/5hmC 混）。

### 圖表 + 自驗
9 張圖（4 主直方 X=差異CpG數/位點 Y=SNV數 TP/FP + per-軸 + 方向 + 強度 + FP/TP 歸因 + CDF）：`obs_ws/cpp_wg/20260624_differential_methylation_marker_charts.standalone.html`。資料 `percpg_cpg_classification.json` / `percpg_per_locus.json` / `fptp_attribution.json`。腳本 `extract_percpg_quant.py` / `build_quant_charts.py`。

> 紅線：site-first 足夠 cis-ASM ≠ read-level 無用（它做發現）；523 subclone_novel = 甲基子群非 somatic subclone（FP 富集）；Fisher per-CpG sound-for-purpose 量級偏 anti-conservative；多切法=sensitivity。**「有差異甲基」普遍(96.8%)不判別，TP 專一軸(somatic HP-1 70.8%)才是真定位點**；單樣本 HCC1395 ⭐2-3。
