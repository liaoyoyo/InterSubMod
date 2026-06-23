<!--
建立時間: 2026-06-24
類型: methodology — 全位點 subclone 分類決策樹（修正版 canonical）
狀態: in_progress（單樣本 HCC1395 探索 ⭐2-3，全程 partial）
build_branch: feat/summary-nreadsvalid
build_commit: d6bd243
data_sources: docs/methodology/_assets/20260618_subcluster_pilot/corrected_tree.json, docs/methodology/_assets/20260618_subcluster_pilot/contingency_summary.json, docs/methodology/_assets/20260618_subcluster_pilot/no_structure_breakdown.json, docs/methodology/_assets/20260618_subcluster_pilot/phylo_cpp_wg_full_records_v4.json
取代: 本檔取代先前對話中「無結構=24.7% unevaluable」的錯誤拆解（誤把 {HP1,HP1-1} somatic 軸當單標籤）
-->

# 全位點 subclone 分類決策樹（修正版 canonical）— HCC1395 single-sample

> **一句話**：把 34,736 個 somatic-SNV 窗依「**標籤組成 × 甲基結構 × cluster×label 列聯型**」三維分類，逐節點給絕對比例（/34,736）+ 排除理由 + 自驗指令。**⭐2-3 單樣本探索，觀察/刻畫層，非 subclone caller。**

## 0. Scope / 紅線

- 單樣本 HCC1395（tumor-only + SEQC2 TP/FP 標籤）；分群來自 C++ phylo-v4.1（null95 + CramérV），距離 = BERNOULLI 加權算術平均。
- **對齊（CramérV≥0.3）= cis-ASM 跡象，非 subclone**；**確認 subclone 需單細胞/multi-region/第二 somatic 共現**（單樣本無）。
- 比例一律**絕對於 34,736**（非條件比例），方便自驗。

## 1. 生物學前提驗證（你的斷言：HP1-1 ⟹ ALT）

`-1` 子單倍型由 somatic 變異定義 → `{HP1,HP1-1}` = `{REF,ALT}` somatic 軸。實測：

| 項 | 數 |
|---|---|
| 有 somatic 子單倍型（HP1-1/HP2-1）的位點 | 34,595 |
| 其中同時有 ALT | 34,372 |
| **一致率** | **99.36%** ✓ |

→ 「somatic 軸（REF vs ALT）可測」的定義成立；99.6% 位點都有 somatic 軸，**真正單標籤罕見**。

## 2. 第一層：8 類分類（標籤組成 × 甲基結構），絕對比例

```
ROOT 34,736 (100.00%)  [TP 30,077 / FP 4,659]
│
├─ 甲基結構 coarse≥2 ················· 10,556 (30.39%)
│   ├─ 有軸（germline/somatic）
│   │   ├─ C·S1 對齊=cis-ASM ········· 6,857 (19.74%)
│   │   ├─ C·S2 未對齊 ··············· 1,569 ( 4.52%)
│   │   └─ C·S3 不穩定(modal<0.7) ···· 1,407 ( 4.05%)
│   └─ 單標籤（單克隆內再分群）
│       └─ A·⭐SUBCLONE 候選 ········· 723 ( 2.08%)
│
└─ 甲基結構 coarse=1 ················· 24,180 (69.61%)
    ├─ 有軸 → D 可測無結構 ········· 21,963 (63.23%)
    └─ 單標籤 → B 真單群/無法區分 ··· 2,217 ( 6.38%)
        ├─ B1 LOH（生物同質）········ 1,978 ( 5.69%)
        ├─ B3 非LOH足覆蓋（生物同質）· 193 ( 0.56%)
        └─ B2 低覆蓋（技術無法處理）·· 46 ( 0.13%)

sum-check: 6857+1569+1407+723+21963+1978+193+46 = 34,736 ✓
```

### 各類排除理由（以「確認 somatic subclone」為目標）

| 類 | n (絕對%) | 主因 → 細項 |
|---|---|---|
| **A subclone 候選** | 723 (2.08%) | **保留**：單一克隆標籤內甲基再分群 = tumor-only 唯一真 subclone 訊號軸（within-HP）|
| C·S1 對齊 | 6,857 (19.74%) | 排除為 subclone：結構可歸因 a-priori 軸 = **cis-ASM**（germline/somatic 等位特異甲基）|
| C·S2 未對齊 | 1,569 (4.52%) | 混合桶（見 §3 列聯型細拆）|
| C·S3 不穩定 | 1,407 (4.05%) | modal_frac<0.7 = seed-borderline |
| D 可測無結構 | 21,963 (63.23%) | **真單群（甲基層）**：有 REF/ALT 或 HP1/HP2 可測、但單一甲基群體 |
| B1 LOH | 1,978 (5.69%) | **無法區分（生物性）**：LOH 丟一單倍型 → 只剩單標籤 → 形不成比較 |
| B3 非LOH足覆蓋 | 193 (0.56%) | 無法區分（生物性）：天生單標籤 |
| B2 低覆蓋 | 46 (0.13%) | **無法區分（技術性）**：reads 太少才只取到單標籤 |

**🔴 修正紀錄**：先前誤把「無結構」拆成 24.7% unevaluable；真相 = 其中 6,374 個是 `{HP1,HP1-1}`（有 somatic 軸）= 可測（歸 D），**真正單標籤僅 B=2,217 (6.38%)**，且**幾乎全生物性**（B1+B3=2,171）、技術性僅 B2=46 (0.13%)。

## 3. 第二層：有結構位點的 cluster×label 列聯型（10,556 個再細拆）

對每個 coarse≥2 位點建 coarse_label×hp 列聯表（MIN_CELL=3）：

| 列聯型 | n (絕對%) | /有結構 | 意義 |
|---|---|---|---|
| mixed（①②並存）| 3,980 (11.46%) | 37.7% | 最複雜：部分標籤裂群 + 部分群混標籤 |
| **many:1 結構>標籤 ①** | 3,187 (9.17%) | 30.2% | **subclone-like**（結構超出標籤）|
| one2one ✅ | 1,739 (5.01%) | 16.5% | 乾淨對齊 = cis-ASM（驗證：99.6% 是 S1）|
| **1:many 跨標籤 ②** | 1,650 (4.75%) | 15.6% | 甲基跨標籤 = 無-ASM/trans |

sum-check 3980+3187+1739+1650 = 10,556 ✓

### 列聯型 × 標籤主軸（定位最乾淨 subclone）

| | many:1 ①（subclone-like）| 1:many ②（跨標籤）|
|---|---|---|
| single_label somatic（純克隆內）| **1,043 ⭐ 最乾淨 subclone 候選** | — |
| somatic 單家族 | 642 | 98 |
| germline 兩家族 | 1,453（單倍型內異質, 模糊）| 1,552（跨單倍型合併=無 ASM）|

**兩個關鍵不對稱**：① many:1×single-somatic=1,043 = 一個 somatic 克隆內再分群、無 germline 污染。② 1:many 幾乎全 germline（1,552 vs somatic 98）→ 甲基跨標籤合併只發生在跨 germline 單倍型（無 ASM），極少跨 somatic 邊界。

### S2「未對齊」其實過半是 subclone-like

S2 (2,190) = many:1 **1,210** + mixed 594 + 1:many 380 + one2one 6 → **過半（55%）是「結構超出標籤」= subclone 訊號**，不該與「跨標籤無-ASM」混桶。

## 4. 精煉 subclone 候選層級（取代單一「A=723」）

| 級 | 定義 | n |
|---|---|---|
| 🥇 最乾淨 | many:1 × single-label-somatic | 1,043 |
| 🥈 次之 | many:1 × somatic 單家族 | 642（合 somatic-軸 subclone-like **1,685**）|
| ⚠ 模糊 | many:1 × germline（需 normal-control 排表觀噪音）| 1,453 |
| 正交佳例 | LOH + somatic 軸 + 結構（癌突變對齊好說明）| 3,553 (10.23%) |

> **仍須**：① normal-control 扣 germline；② 第二 somatic 共現/CN-segment/single-cell 才可宣稱 subclone。本層只「定位候選」非「確認」。

## 5. Provenance（資料/腳本溯源）

| 產物 | 來源腳本 | 落檔 |
|---|---|---|
| per-locus 真值（34,736）| C++ `inter_sub_mod`（phylo-v4.1, commit c186092）| `_assets/20260618_subcluster_pilot/phylo_cpp_wg_full_records_v4.json` |
| 標籤組成 + 8 類 | `scripts/extract_label_composition.py` | `corrected_tree.json` / `label_composition.json` |
| 列聯型 | `scripts/extract_contingency_type.py` | `contingency_summary.json` / `contingency_type.json` |
| 無結構拆解 | `scripts/extract_hp_allele_diversity.py` | `no_structure_breakdown.json` |

路徑前綴 `InterSubMod/docs/methodology/_assets/20260618_subcluster_pilot/`。

## 6. 自驗指令（任何時候可重算驗證）

```bash
cd InterSubMod/.claude/worktrees/ism-review-infra/docs/methodology/_assets/20260618_subcluster_pilot
# 8 類計數 + sum-check
python3 -c "import json;d=json.load(open('corrected_tree.json'));print({k:v['n'] for k,v in d['categories_abs'].items()});print('sum',d['sum_check'])"
# 列聯型計數
python3 -c "import json;d=json.load(open('contingency_summary.json'));print({k:v['n'] for k,v in d['type_abs'].items()})"
# 生物學一致率
python3 -c "import json;print(json.load(open('corrected_tree.json'))['biology_check_HP1_1_implies_ALT'])"
# 從 records 重算任一類（例: A subclone 候選 = 單標籤 + coarse≥2）
python3 -c "import json;r=json.load(open('phylo_cpp_wg_full_records_v4.json'));print('coarse>=2:',sum(1 for x in r if x['coarse_ng']>=2))"
```

> 每個比例 = n / 34,736。重跑腳本需 `output/_phylo_wg_full/` 保留的 region dirs（matrices/reads/phylo_groups.tsv）。
