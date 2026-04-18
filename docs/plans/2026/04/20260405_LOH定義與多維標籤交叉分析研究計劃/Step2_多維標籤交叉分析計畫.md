<!--
建立時間: 2026-04-05 23:30
目標: Step 2 — HP × ALT/REF 多維標籤交叉分析
處理範圍: 從 reads.tsv 提取交叉計數、self-phasing 量化、甲基化距離分析
關聯檔案:
  - docs/plans/2026/04/20260405_LOH定義與多維標籤交叉分析研究計劃/00_總覽與執行順序.md
-->

# Step 2: 多維標籤交叉分析計畫

## 啟動條件

1. **Step 1 完成並經使用者確認方向**才啟動 — Step 1 的四象限定義（Q1-Q4）是本 Step 的分層基礎
2. **必須有 Step 1.2 的四象限欄位** `loh_quadrant`（Q1_both / Q2_ism_only / Q3_bed_only / Q4_neither）已合併至 master dataset
3. **LongPhase-TO baseline** 為唯一基準（與 Step 1 一致），PON-only 僅作觀察對照

---

## 子任務 2.1: HP x ALT/REF 交叉計數提取

### 目標

從每個 variant 對應的 per-region `reads.tsv` 提取 HP x ALT/REF 交叉計數，擴充 master dataset。

### 資料來源

- **Master dataset**: `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz`
  - 748,391 rows x 116 columns
  - 關鍵欄位:
    - `source_region_root` — 指向 region 根目錄
    - `Chr`, `Pos` — 用於定位 reads.tsv 路徑
    - `variant_key` — 用於 merge key（格式: `chr1:1288697:C:T`）
    - `HP1FamilyN`, `HP2FamilyN`, `NHP3`, `NHP0` — 用於驗證合計一致性
- **reads.tsv 路徑模式**: `{source_region_root}/{Chr}/{Chr}_{Pos}/{Chr}_{Pos-5000}_{Pos+5000}/reads/reads.tsv`
- **reads.tsv 欄位**:

  | 欄位 | 說明 | 範例值 |
  |------|------|--------|
  | `hp` | HP 標籤（LongPhase 輸出） | `"1"`, `"2"`, `"1-1"`, `"2-1"`, `"3"`, `"0"` |
  | `alt_support` | 該 read 在 variant 位置的 allele 支持 | `"ALT"`, `"REF"`, `"UNKNOWN"` |

### 提取方法

#### HP 分組定義

```python
HP_GROUPS = {
    "HP1":   lambda hp: hp == "1",
    "HP2":   lambda hp: hp == "2",
    "HP1_1": lambda hp: hp == "1-1",
    "HP2_1": lambda hp: hp == "2-1",
    "HP3":   lambda hp: hp == "3",
    "HP0":   lambda hp: hp == "0",
}

ALT_GROUPS = ["ALT", "REF", "UNKNOWN"]
```

#### 新增 18 欄交叉計數

| HP group | ALT 欄 | REF 欄 | UNK 欄 |
|----------|--------|--------|--------|
| HP1 | `HP1_ALT` | `HP1_REF` | `HP1_UNK` |
| HP2 | `HP2_ALT` | `HP2_REF` | `HP2_UNK` |
| HP1_1 | `HP1_1_ALT` | `HP1_1_REF` | `HP1_1_UNK` |
| HP2_1 | `HP2_1_ALT` | `HP2_1_REF` | `HP2_1_UNK` |
| HP3 | `HP3_ALT` | `HP3_REF` | `HP3_UNK` |
| HP0 | `HP0_ALT` | `HP0_REF` | `HP0_UNK` |

#### 並行處理架構

```
multiprocessing.Pool(cpu_count())
    ├── Worker 1: 處理 chunk 1 (N/cpu_count rows)
    ├── Worker 2: 處理 chunk 2
    ├── ...
    └── Worker K: 處理 chunk K
每個 Worker:
    for row in chunk:
        reads_tsv_path = build_path(row)
        if not exists(reads_tsv_path): record missing, skip
        df_reads = pd.read_csv(reads_tsv_path, sep='\t', usecols=['hp','alt_support'])
        cross_tab = pd.crosstab(df_reads['hp'], df_reads['alt_support'])
        extract 18 counts → append to result
    return DataFrame(results)
```

### Python 實作模板

```python
#!/usr/bin/env python3
"""Step 2.1: Extract HP x ALT/REF cross-tabulation from per-region reads.tsv.

Reads the master dataset, locates each region's reads.tsv,
builds an 18-column cross-count table, and merges back.

Output: augmented_master.tsv.gz
"""

import sys
import warnings
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import (
    MASTER_DATASET_PATH, OUTPUT_ROOT, ensure_dir, load_master_dataset,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

OUT_DIR = ensure_dir(OUTPUT_ROOT / "20260407_multidim_tag_cross")

HP_LABELS = ["1", "2", "1-1", "2-1", "3", "0"]
HP_COL_PREFIXES = ["HP1", "HP2", "HP1_1", "HP2_1", "HP3", "HP0"]
ALT_LABELS = ["ALT", "REF", "UNKNOWN"]
ALT_COL_SUFFIXES = ["ALT", "REF", "UNK"]

CROSS_COLUMNS = [
    f"{hp}_{alt}"
    for hp in HP_COL_PREFIXES
    for alt in ALT_COL_SUFFIXES
]  # 18 columns total

WINDOW_HALF = 5000  # ISM default: position +/- 5000 bp


# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

def build_reads_tsv_path(source_region_root: str, chrom: str, pos: int) -> Path:
    """Construct the reads.tsv path from region root, chromosome, and position.

    Pattern: {root}/{chr}/{chr}_{pos}/{chr}_{pos-5000}_{pos+5000}/reads/reads.tsv
    """
    start = pos - WINDOW_HALF
    end = pos + WINDOW_HALF
    return (
        Path(source_region_root)
        / chrom
        / f"{chrom}_{pos}"
        / f"{chrom}_{start}_{end}"
        / "reads"
        / "reads.tsv"
    )


# ---------------------------------------------------------------------------
# Single-region extraction
# ---------------------------------------------------------------------------

def extract_cross_counts(reads_tsv_path: Path) -> Dict[str, int]:
    """Read a single reads.tsv and return the 18-element cross-count dict."""
    result = {col: 0 for col in CROSS_COLUMNS}

    if not reads_tsv_path.exists():
        return result  # all zeros; caller tracks missing count

    try:
        df = pd.read_csv(
            reads_tsv_path, sep="\t",
            usecols=["hp", "alt_support"],
            dtype={"hp": str, "alt_support": str},
        )
    except Exception:
        return result

    for hp_label, hp_prefix in zip(HP_LABELS, HP_COL_PREFIXES):
        mask_hp = df["hp"] == hp_label
        for alt_label, alt_suffix in zip(ALT_LABELS, ALT_COL_SUFFIXES):
            col_name = f"{hp_prefix}_{alt_suffix}"
            result[col_name] = int((mask_hp & (df["alt_support"] == alt_label)).sum())

    return result


# ---------------------------------------------------------------------------
# Worker function (processes a chunk of rows)
# ---------------------------------------------------------------------------

def _worker(args: Tuple[int, pd.DataFrame]) -> pd.DataFrame:
    """Process a chunk of master rows; return DataFrame with cross columns."""
    chunk_id, chunk_df = args
    records = []
    missing_count = 0

    for idx, row in chunk_df.iterrows():
        path = build_reads_tsv_path(
            row["source_region_root"], row["Chr"], int(row["Pos"])
        )
        counts = extract_cross_counts(path)

        if not path.exists():
            missing_count += 1

        counts["_master_idx"] = idx
        records.append(counts)

    result_df = pd.DataFrame(records)
    if missing_count > 0:
        warnings.warn(
            f"[Chunk {chunk_id}] {missing_count}/{len(chunk_df)} "
            f"reads.tsv files not found"
        )
    return result_df


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 72)
    print("Step 2.1: HP x ALT/REF Cross-Count Extraction")
    print("=" * 72)

    # 1. Load master dataset
    df = load_master_dataset()

    required_cols = ["source_region_root", "Chr", "Pos", "variant_key",
                     "HP1FamilyN", "HP2FamilyN", "NHP3", "NHP0"]
    for col in required_cols:
        assert col in df.columns, f"Missing required column: {col}"

    # 2. Split into chunks for multiprocessing
    n_workers = cpu_count()
    chunks = np.array_split(df, n_workers)
    args_list = [(i, chunk) for i, chunk in enumerate(chunks)]

    print(f"[Main] Dispatching {len(df):,} rows to {n_workers} workers ...")

    # 3. Parallel extraction
    with Pool(n_workers) as pool:
        results = list(tqdm(
            pool.imap_unordered(_worker, args_list),
            total=n_workers,
            desc="Extracting cross-counts",
        ))

    # 4. Combine results and merge back
    cross_df = pd.concat(results, ignore_index=False).sort_values("_master_idx")
    cross_df = cross_df.set_index("_master_idx")
    cross_df.index.name = None

    for col in CROSS_COLUMNS:
        df[col] = cross_df[col].values

    # 5. Validation
    print("\n[Validation] Checking HP family count consistency ...")

    df["_HP1_sum"] = df["HP1_ALT"] + df["HP1_REF"] + df["HP1_UNK"]
    df["_HP2_sum"] = df["HP2_ALT"] + df["HP2_REF"] + df["HP2_UNK"]
    df["_HP1_1_sum"] = df["HP1_1_ALT"] + df["HP1_1_REF"] + df["HP1_1_UNK"]
    df["_HP2_1_sum"] = df["HP2_1_ALT"] + df["HP2_1_REF"] + df["HP2_1_UNK"]
    df["_HP3_sum"] = df["HP3_ALT"] + df["HP3_REF"] + df["HP3_UNK"]
    df["_HP0_sum"] = df["HP0_ALT"] + df["HP0_REF"] + df["HP0_UNK"]

    # HP1FamilyN should equal HP1 + HP1_1
    df["_hp1_family_check"] = df["_HP1_sum"] + df["_HP1_1_sum"]
    hp1_mismatch = (df["_hp1_family_check"] != df["HP1FamilyN"]).sum()

    # HP2FamilyN should equal HP2 + HP2_1
    df["_hp2_family_check"] = df["_HP2_sum"] + df["_HP2_1_sum"]
    hp2_mismatch = (df["_hp2_family_check"] != df["HP2FamilyN"]).sum()

    print(f"  HP1FamilyN mismatch: {hp1_mismatch:,} / {len(df):,}")
    print(f"  HP2FamilyN mismatch: {hp2_mismatch:,} / {len(df):,}")

    if hp1_mismatch > 0 or hp2_mismatch > 0:
        print("  [WARNING] Mismatches detected — investigate before proceeding")

    # Drop validation columns
    df.drop(columns=[c for c in df.columns if c.startswith("_")], inplace=True)

    # 6. Save
    out_path = OUT_DIR / "augmented_master.tsv.gz"
    print(f"\n[Save] Writing augmented master to {out_path} ...")
    df.to_csv(out_path, sep="\t", index=False, compression="gzip")
    print(f"[Save] Done — {len(df):,} rows x {len(df.columns)} cols")

    # 7. Summary statistics
    print("\n[Summary] Cross-count column statistics (TO mode):")
    df_to = df[df["mode"] == "to"]
    for col in CROSS_COLUMNS:
        vals = df_to[col]
        print(f"  {col:>12s}: mean={vals.mean():.2f}, median={vals.median():.0f}, "
              f"max={vals.max():.0f}, zero_frac={( vals == 0).mean():.3f}")


if __name__ == "__main__":
    main()
```

### 輸出

| 檔案 | 說明 |
|------|------|
| `{OUTPUT_ROOT}/20260407_multidim_tag_cross/augmented_master.tsv.gz` | 原始 116 欄 + 18 欄交叉計數 |

### 效能預估

| 項目 | 估計 |
|------|------|
| 需讀取的 reads.tsv 檔案數 | ~748K（每行一個） |
| 單檔讀取時間（SSD） | ~1-3 ms（usecols 2 欄） |
| 總 I/O 時間（單核） | ~12-37 分鐘 |
| 並行加速（32 核） | ~0.5-1.5 分鐘 |
| 記憶體峰值 | ~2-4 GB（master dataset + cross counts） |

### 驗證步驟

1. **合計一致性**: `HP1_ALT + HP1_REF + HP1_UNK + HP1_1_ALT + HP1_1_REF + HP1_1_UNK == HP1FamilyN`
2. **合計一致性**: `HP2_ALT + HP2_REF + HP2_UNK + HP2_1_ALT + HP2_1_REF + HP2_1_UNK == HP2FamilyN`
3. **合計一致性**: `HP3_ALT + HP3_REF + HP3_UNK == NHP3`
4. **合計一致性**: `HP0_ALT + HP0_REF + HP0_UNK == NHP0`
5. **Missing file rate**: 記錄 reads.tsv 找不到的比例（預期 < 1%）
6. **Random sample 手動驗證**: 隨機抽 10 個 region，手動比對 reads.tsv 內容

---

## 子任務 2.2: Self-Phasing 偏差量化

### 目標

用 HP x ALT 交叉表量化 self-phasing 偏差程度，驗證 TO 模式中 ALT reads 偏向單一 HP 的現象。

### 前置依賴

- 子任務 2.1 完成（需要 18 欄交叉計數）

### 核心指標定義

```python
EPSILON = 1e-6  # avoid division by zero

# ---- Primary HP（HP1, HP2）的 skew ----

# ALT reads 在 HP1 vs HP2 之間的偏斜程度
# 理想 paired: ALT 均勻分佈在兩個 HP → skew ≈ 0
# Self-phasing TO: ALT 幾乎全部在一個 HP → skew → 1.0
ALT_HP_skew = |HP1_ALT - HP2_ALT| / (HP1_ALT + HP2_ALT + EPSILON)

# REF reads 的 HP 偏斜（對照組）
REF_HP_skew = |HP1_REF - HP2_REF| / (HP1_REF + HP2_REF + EPSILON)

# ---- 衍生指標 ----

# ALT 與 REF 的 skew 差異
# Self-phasing: ALT skew >> REF skew
skew_gap = ALT_HP_skew - REF_HP_skew

# 含 sub-HP 的完整 ALT skew（HP1+HP1_1 vs HP2+HP2_1）
ALT_family_skew = |(HP1_ALT + HP1_1_ALT) - (HP2_ALT + HP2_1_ALT)|
                  / (HP1_ALT + HP1_1_ALT + HP2_ALT + HP2_1_ALT + EPSILON)

# ALT reads 佔主 HP 的比例（self-phasing 強度）
dominant_HP_ALT_frac = max(HP1_ALT, HP2_ALT) / (HP1_ALT + HP2_ALT + EPSILON)
```

### 預期行為

| 模式 | ALT_HP_skew | REF_HP_skew | skew_gap | 解釋 |
|------|-------------|-------------|----------|------|
| **Paired** | 低（~0.0-0.3） | 低（~0.0-0.3） | ~0 | 真正 phasing：ALT/REF 都分配到正確 HP |
| **TO（真 somatic）** | 高（→1.0） | 低-中 | 高 | Self-phasing：ALT 被偏向一個 HP，REF 相對均勻 |
| **TO（germline FP）** | 高（→1.0） | 高（→1.0） | ~0 | 真 heterozygous：兩邊都高 skew 但方向一致 |
| **LOH 區域** | 高 | 高 | 低-高 | 生物學 LOH：只有一個 haplotype 有 reads |

### 分析步驟

1. **計算指標**：對 augmented master 所有行計算上述 5 個 skew 指標
2. **TO vs Paired 分佈比較**：確認 self-phasing 效應的模式差異
3. **TP vs FP 差異**：在 TO 模式內比較 TP/FP 的 skew 分佈
4. **LOH 四象限交叉**：在 Step 1 的 Q1-Q4 四象限內分別觀察 skew 分佈
5. **與 HPMergedDelta 的相關性**：驗證 skew 指標是否提供 HPMergedDelta 之外的額外資訊
6. **AUC 評估**：計算各 skew 指標區分 TP/FP 的 AUC（per sample, per mode）

### 圖表規格

| 圖編號 | 內容 | 格式 | 維度 |
|--------|------|------|------|
| F1 | ALT_HP_skew 分佈：Paired vs TO | Split violin plot | 7 samples, 2x1 panel (mode) |
| F2 | ALT_HP_skew 分佈：TP vs FP (TO only) | Split violin plot | 7 samples panel |
| F3 | skew_gap (ALT_skew - REF_skew) 分佈：Paired vs TO, TP vs FP | 2x2 facet violin | 7 samples |
| F4 | ALT_HP_skew vs HPMergedDelta 散佈圖 | Scatter + density | Paired/TO 各一，colored by truth_label |
| F5 | ALT_HP_skew 在四象限內的 boxplot (TO only) | Grouped boxplot | 7 samples x 4 quadrants |
| F6 | AUC heatmap: sample x skew_metric (TO only) | Heatmap | 7 samples x 5 metrics |

**圖表共通規格**:
- 解析度: 180 dpi
- 7 samples 使用 `SAMPLE_ORDER` 固定順序
- 色盤: `TRUTH_PALETTE`（TP/FP）、`MODE_PALETTE`（Paired/TO）
- 四象限色盤: Q1=紫、Q2=橙、Q3=綠、Q4=灰

### 輸出

| 檔案 | 說明 |
|------|------|
| `{OUT_DIR}/self_phasing_skew_statistics.tsv` | Per-sample, per-mode skew 統計（mean, median, std, AUC） |
| `{OUT_DIR}/F1_alt_hp_skew_mode_comparison.png` | 圖 F1 |
| `{OUT_DIR}/F2_alt_hp_skew_tp_fp_to.png` | 圖 F2 |
| `{OUT_DIR}/F3_skew_gap_facet.png` | 圖 F3 |
| `{OUT_DIR}/F4_skew_vs_hpmerged_scatter.png` | 圖 F4 |
| `{OUT_DIR}/F5_skew_by_quadrant_boxplot.png` | 圖 F5 |
| `{OUT_DIR}/F6_skew_auc_heatmap.png` | 圖 F6 |

---

## 子任務 2.3: LOH x HP x ALT 聯合分析

### 目標

在 Step 1 的 LOH 四象限（Q1-Q4）內深入分析 HP x ALT/REF 交叉模式，理解不同 LOH 定義下的 read 分配行為差異。

### 前置依賴

- 子任務 2.1 完成（18 欄交叉計數）
- Step 1.2 完成（`loh_quadrant` 欄位可用）

### 分析步驟

#### 2.3.1 四象限交叉計數 Profile

針對每個象限計算標準化的 HP x ALT 交叉 Profile：

```python
# 每象限內的平均交叉計數（標準化為佔 NumReads 的比例）
for quadrant in ['Q1_both', 'Q2_ism_only', 'Q3_bed_only', 'Q4_neither']:
    subset = df_to[df_to['loh_quadrant'] == quadrant]
    for col in CROSS_COLUMNS:
        profile[quadrant][col] = (subset[col] / subset['NumReads']).median()
```

**預期結果**:

| 象限 | 預期 ALT 分佈模式 | 預期 REF 分佈模式 |
|------|-------------------|-------------------|
| Q1_both | ALT 極度偏向一個 HP | REF 也偏向（生物學 LOH） |
| Q2_ism_only | ALT 極度偏向一個 HP | REF 相對均勻（self-phasing） |
| Q3_bed_only | ALT 較均勻 | REF 較均勻（在 LOH.bed 但 HP 平衡） |
| Q4_neither | ALT 較均勻 | REF 較均勻（無 LOH） |

#### 2.3.2 Self-Phasing Signature 驗證

Q2_ism_only 象限是 self-phasing 的高嫌疑區。驗證特徵：

1. `ALT_HP_skew >> REF_HP_skew`（ALT 偏但 REF 不偏）
2. `dominant_HP_ALT_frac > 0.9`（ALT 幾乎全在一個 HP）
3. `HP1_REF ≈ HP2_REF`（REF 大致平衡）
4. 此模式在 **TO** 顯著但在 **Paired** 不顯著

計算 Q2 內滿足以上 4 個條件的比例 → **self-phasing purity ratio**

#### 2.3.3 TP/FP 分離度 per 象限

在每個象限內，用 `ALT_HP_skew` 和 `skew_gap` 計算 TP/FP 的 AUC：

```python
for quadrant in quadrants:
    for metric in ['ALT_HP_skew', 'skew_gap', 'ALT_family_skew']:
        auc = compute_auc(subset_tp[metric], subset_fp[metric])
        # Record per (sample, quadrant, metric)
```

**關鍵問題**: skew 指標在哪個象限能最好地區分 TP/FP？

#### 2.3.4 Sub-HP 分析

探索 HP1-1 和 HP2-1（sub-HP）在不同象限中的行為：

- Q1（LOH.bed 內）: sub-HP 比例是否更高？
- 是否有 ALT reads 集中在 sub-HP 的模式（例如 `HP2_1_ALT >> HP2_ALT`）？

### 圖表規格

| 圖編號 | 內容 | 格式 |
|--------|------|------|
| F7 | 四象限 HP x ALT normalized profile heatmap | 4 panel heatmap（Q1-Q4），x=ALT/REF/UNK, y=HP groups |
| F8 | Q2 self-phasing signature 驗證：ALT_skew vs REF_skew scatter | Scatter (TO only)，colored by truth_label |
| F9 | AUC per quadrant x metric heatmap | Heatmap: 4 quadrants x 5 metrics，7 samples 小多圖 |
| F10 | Sub-HP ALT fraction by quadrant | Grouped bar chart |

### 輸出

| 檔案 | 說明 |
|------|------|
| `{OUT_DIR}/quadrant_cross_profile.tsv` | 四象限標準化交叉 Profile 統計 |
| `{OUT_DIR}/self_phasing_purity_ratio.tsv` | Q2 self-phasing 純度比例 per sample |
| `{OUT_DIR}/quadrant_auc_matrix.tsv` | 每象限 x 每 metric 的 AUC |
| `{OUT_DIR}/F7_quadrant_profile_heatmap.png` | 圖 F7 |
| `{OUT_DIR}/F8_q2_self_phasing_scatter.png` | 圖 F8 |
| `{OUT_DIR}/F9_quadrant_auc_heatmap.png` | 圖 F9 |
| `{OUT_DIR}/F10_sub_hp_alt_fraction.png` | 圖 F10 |

---

## 子任務 2.4: 甲基化距離 x 標籤組合

### 目標

分析甲基化距離指標（PairwiseMedianDist, CramersV, HPMergedDelta, AlleleDelta）在不同 HP x ALT 組合下的行為，理解 self-phasing 如何影響甲基化分析結果。

### 前置依賴

- 子任務 2.1 完成（18 欄交叉計數）
- 子任務 2.2 完成（skew 指標）

### 分析步驟

#### 2.4.1 甲基化指標 vs Skew 相關性

```python
methyl_metrics = [
    'PairwiseMedianDist', 'CramersV',
    'HPMergedDelta', 'AlleleDelta',
]
skew_metrics = [
    'ALT_HP_skew', 'REF_HP_skew', 'skew_gap',
    'ALT_family_skew', 'dominant_HP_ALT_frac',
]

# Spearman correlation matrix: methyl_metrics x skew_metrics
# Per mode (Paired/TO), per truth_label (TP/FP)
```

**預期**:
- TO 模式: `HPMergedDelta` 與 `ALT_HP_skew` 高度正相關（self-phasing 驅動 HP 差異）
- Paired 模式: 相關性低（真正的 phasing，HP 差異反映生物學）
- `AlleleDelta` 與 `skew_gap` 正相關（ALT 比 REF 更偏斜 → allele 差異更大）

#### 2.4.2 Skew 分層下的甲基化指標分佈

將 TO 模式按 `ALT_HP_skew` 分為三層：

| 層級 | ALT_HP_skew 範圍 | 預期含義 |
|------|-------------------|----------|
| Low skew | [0.0, 0.3) | HP 平衡（類似 paired） |
| Medium skew | [0.3, 0.7) | 中等偏斜 |
| High skew | [0.7, 1.0] | 極度偏斜（strong self-phasing） |

在每層內比較：
1. HPMergedDelta 分佈（TP vs FP）
2. CramersV 分佈（TP vs FP）
3. TP/FP 的 AUC（用 HPMergedDelta 和 CramersV）

**關鍵假說**: 在 low skew 層，TO 的甲基化指標行為應接近 Paired，TP/FP 區分度可能更好。

#### 2.4.3 LOH 四象限 x Skew 層 x 甲基化指標

三維交叉分析（僅在資料量足夠時執行）：

```
for quadrant in [Q1, Q2, Q3, Q4]:
    for skew_tier in [Low, Medium, High]:
        for metric in methyl_metrics:
            compare TP vs FP distribution
            compute AUC
```

此分析可揭示：某些象限+skew 組合是否有足夠的 TP/FP 區分度來建立局部篩選規則。

### 圖表規格

| 圖編號 | 內容 | 格式 |
|--------|------|------|
| F11 | Spearman correlation heatmap: methyl x skew | 2x2 facet (mode x truth), heatmap |
| F12 | HPMergedDelta by skew tier: TP vs FP violin (TO only) | 3-panel violin (Low/Med/High skew) |
| F13 | CramersV by skew tier: TP vs FP violin (TO only) | 3-panel violin |
| F14 | AUC trend: skew tier x methyl metric (TO only) | Line plot, per sample |

### 輸出

| 檔案 | 說明 |
|------|------|
| `{OUT_DIR}/methyl_skew_correlation.tsv` | Spearman 相關矩陣 per mode x truth |
| `{OUT_DIR}/skew_tier_methyl_statistics.tsv` | 每 skew 層的甲基化指標統計與 AUC |
| `{OUT_DIR}/F11_methyl_skew_correlation_heatmap.png` | 圖 F11 |
| `{OUT_DIR}/F12_hpmerged_by_skew_tier.png` | 圖 F12 |
| `{OUT_DIR}/F13_cramersv_by_skew_tier.png` | 圖 F13 |
| `{OUT_DIR}/F14_auc_by_skew_tier.png` | 圖 F14 |

---

## 腳本對應

| 腳本檔案 | 對應子任務 | 輸出目錄 |
|----------|-----------|----------|
| `scripts/analysis/build_hp_alt_cross_extraction.py` | 2.1 | `20260407_multidim_tag_cross/` |
| `scripts/analysis/build_self_phasing_quantification.py` | 2.2, 2.3, 2.4 | `20260407_multidim_tag_cross/` |

**分拆理由**: 2.1 是 I/O 密集型（讀取 ~748K 個 reads.tsv），執行一次後輸出 augmented master，後續 2.2-2.4 都基於此 augmented master 進行純計算分析，無需重複 I/O。

---

## 驗證清單

### 資料完整性

- [ ] `HP1_ALT + HP1_REF + HP1_UNK + HP1_1_ALT + HP1_1_REF + HP1_1_UNK == HP1FamilyN`（per row，容許 0 mismatch）
- [ ] `HP2_ALT + HP2_REF + HP2_UNK + HP2_1_ALT + HP2_1_REF + HP2_1_UNK == HP2FamilyN`（per row，容許 0 mismatch）
- [ ] `HP3_ALT + HP3_REF + HP3_UNK == NHP3`（per row）
- [ ] `HP0_ALT + HP0_REF + HP0_UNK == NHP0`（per row）
- [ ] reads.tsv missing rate < 1%
- [ ] Random sample 10 regions 手動對照 reads.tsv 原檔

### Self-Phasing 量化

- [ ] ALT_HP_skew 在 Paired 模式的 median < 0.3
- [ ] ALT_HP_skew 在 TO 模式的 median > 0.5（預期因 self-phasing 偏高）
- [ ] skew_gap（ALT_skew - REF_skew）在 TO > Paired（方向一致性）
- [ ] 以上三點跨 **7 samples 方向一致**

### 交叉分析

- [ ] Q2_ism_only 的 self-phasing purity ratio > 50%（預期 self-phasing 是主因）
- [ ] Q1_both 的 ALT_HP_skew 與 REF_HP_skew 同時高（生物學 LOH）
- [ ] Q4_neither 的所有 skew 指標最低

### 甲基化關聯

- [ ] TO 模式 HPMergedDelta 與 ALT_HP_skew Spearman rho > 0.3
- [ ] Paired 模式 HPMergedDelta 與 ALT_HP_skew Spearman rho < 0.2
- [ ] Low skew tier 的 TO 甲基化指標行為接近 Paired

---

## 風險與緩解

| 風險 | 影響 | 緩解措施 |
|------|------|----------|
| reads.tsv 路徑模式不一致 | 2.1 missing rate 過高 | 先跑小樣本（1 sample）驗證路徑；建立 fallback 搜尋機制 |
| HP 標籤格式未知值 | 交叉計數遺漏 | 先統計所有 unique HP 值，確認無遺漏 |
| NumReads 極低（< 5） | Skew 指標噪音大 | 設定 `min_reads=5` 篩選門檻，分別報告篩選前後結果 |
| Q2/Q3 象限資料量不足 | 統計不穩定 | 報告每象限 N，若 N < 50 則標記為「insufficient」 |
| 三維交叉分析（2.4.3）cell size 過小 | 過度擬合風險 | 僅在 cell N > 30 時報告 AUC |

---

## 時間預估

| 子任務 | 預估時間 | 瓶頸 |
|--------|----------|------|
| 2.1 I/O 提取 | 2-5 分鐘（32 核並行） | 磁碟 I/O |
| 2.1 驗證 | 5 分鐘 | 手動抽查 |
| 2.2 分析 + 圖表 | 15-20 分鐘 | 圖表繪製 |
| 2.3 分析 + 圖表 | 15-20 分鐘 | 圖表繪製 |
| 2.4 分析 + 圖表 | 20-30 分鐘 | 三維交叉計算 |
| **合計** | **~60-80 分鐘** | |

---

## 後續銜接

### 向 Step 3.4 輸出

- Self-phasing 量化結果 → 方法修正提案的基礎（如何校正 self-phasing 對 PERMANOVA 的影響）
- Skew tier 分層策略 → 可能成為 QS v3 的輸入特徵

### 向 Step 4 輸出

- Q1_both（生物學 LOH）的 ALT+REF 同時偏斜模式 → 驗證癌症 LOH 的表觀遺傳效應
- Q2 self-phasing purity ratio → 量化 TO 模式的可靠性上限

### 關鍵決策點

Step 2 完成後需回答的問題：

1. **Self-phasing 是否是 TO 模式 HP Imbalance 的主因？** → 如果 Q2 self-phasing purity ratio > 70%，則確認
2. **Skew 指標能否作為 TP/FP filter？** → 如果任何 skew 指標 AUC > 0.65，則值得進一步探索
3. **甲基化指標在 low skew 層是否恢復區分度？** → 如果是，則 "skew-aware" 分析策略有前景
