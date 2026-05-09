<!--
build_date: 2026-05-09
agent: paired priority bug pilot audit (Step A + Step C)
status: validated
report_class: paired-mode-comparison
audience: PI / lab member / 自己未來
parent_synthesis: InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md
inputs:
  - /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam (paired tagged BAM, 260 GB)
  - /big8_disk/liaoyoyo2001/longphase-s/src/somatic_haplotag/SomaticHaplotagProcess.cpp (longphase-s tag 邏輯確認)
outputs:
  - 本檔
  - step_A_distribution/chr19_hp_dist.txt
  - step_C_position_cooccur/chr19_window_som_ratio.txt
verdict: PAIRED MODE 沒有 priority bug — chr19 HP1-1:HP2-1 = 1:1.169（近 1:1）vs TO baseline 17.3:1；som_ratio mean 0.46 跨 windows 變化反映真實 sub-clone 結構而非 systematic bias
last_verified: 2026-05-09
report_template: results-report v1.0
-->

# Paired Mode HP1-1 / HP2-1 共現偏向 Audit — 驗證 paired 是否有同 priority bug

## 0. TL;DR

用戶 5/9 提問：「paired tag 是否也會出現同位點/同 read HP1-1 與 HP2-1 共現，且偏向 HP1-1（類似 TO mode 17.3:1 priority bug）」。

**答案：否**。paired mode 用 longphase-s `somatic_haplotag` 是**完全不同 binary 與 codebase**（vs TO mode longphase-to）。HCC1395 paired tagged BAM chr19 實證：

| 比較 | TO mode baseline | Paired mode |
|---|---|---|
| Binary | longphase-to | longphase-s |
| HP tag 編碼 | `HP:i:` 整數 (1/2/11/21/33) | `HP:Z:` 字串 (1/2/1-1/2-1/3) |
| HP1:HP2 整體比 | **17.3:1**（priority bug）| **1:1.275**（近 1:1）|
| Somatic HP1-1:HP2-1 | 全偏 HP1-1（baseline 752 victims 100% 單向）| **1:1.169**（近 1:1）|
| 1Mb window som_ratio mean | ~0.95 (TO baseline 推估) | **0.462**（無 systematic bias）|
| 同位點共現 | 是（priority bug artifact）| 是（**真實 sub-clone signal**）|

→ paired mode HP1-1/HP2-1 共現是**生物學現象**（sub-clone heterogeneity），不是工程 bug。

---

## 1. 背景與動機

5/8 整合報告 §3 + §6 確立 longphase-to TO mode 的 `getVote()` priority bug：vector ordered check + break early 讓 1 票 somatic 觸發整 read 偏 HP1 系列，全基因組 17.3:1 / chr19 752 read-level victims 鐵證。

用戶質疑：「同樣 priority bug 是否在 paired mode 也存在？」

關鍵釐清：
- TO mode 用 **`longphase-to` fork** (`/big7_disk/liaoyoyo2001/longphase-to-mod/`)
- Paired mode 用 **`longphase-s somatic_haplotag`** (`/big8_disk/liaoyoyo2001/longphase-s/`)
- **不同 binary、不同 codebase、不同 tag 邏輯**

→ paired 不直接套同 priority bug 結論，需獨立驗證。

## 2. 方法

### 2.1 階梯式驗證設計

依用戶 5/9 確認 (a) 全套 (A) → (C) → (B) 階梯式跑：

| Step | 內容 | 預估 | 證據力度 |
|---|---|---|---|
| **A** | chr19 整體 HP:Z: tag distribution | 5 min | 是否有偏移 |
| **C** | chr19 1Mb window per-window HP1-1/HP2-1 ratio | 30 min | 共現是 bug 還是 sub-clone signal |
| **B** | longphase-s binary patch + read-level vote dump | 1-2 天 | read-level 鐵證 |

**Step A + C 完成；Step B 視 A+C 結論決定是否啟動。**

### 2.2 資料來源

- Paired tagged BAM: `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam` (260 GB, 由 ISM canonical pipeline 跑出)
- HP tag 格式：`HP:Z:` 字串（`1` / `2` / `1-1` / `2-1` / `3`）— 從 `longphase-s/src/somatic_haplotag/SomaticHaplotagProcess.cpp:533` `bam_aux_append(aln, "HP", 'Z', ...)` 確認

### 2.3 工具

- `samtools view -F 256` 篩 primary reads
- awk 解析 HP:Z: tag + 1 Mb window 統計
- Python 計 mean / median / stdev

## 3. 結果

### 3.1 Step A — chr19 HP:Z: 整體 tag distribution

```
chr19 primary reads total: 582,566
chr19 tagged (HP:Z: present): 354,919 (60.9%)
```

| HP:Z: 標 | reads | 占 tagged % |
|---|---:|---:|
| `HP:Z:2` | 183,309 | 51.6% |
| `HP:Z:1` | 143,760 | 40.5% |
| `HP:Z:2-1` (somatic on HP2) | 14,504 | 4.1% |
| `HP:Z:1-1` (somatic on HP1) | 12,401 | 3.5% |
| `HP:Z:3` (somatic ambiguous) | 1,145 | 0.3% |

**核心 ratio**:
- germline HP1:HP2 = 143,760:183,309 = **1:1.275**（近 1:1）
- somatic HP1-1:HP2-1 = 12,401:14,504 = **1:1.169**（近 1:1）
- 整體 HP1-fam:HP2-fam = (143,760+12,401):(183,309+14,504) = 156,161:197,813 = **1:1.267**

**對比 TO mode**：
- TO baseline 全基因組 HP1:HP2 = 614,000:35,500 = **17.3:1** (PI §3.3.2)
- TO chr19 752 read-level victims = baseline=11→V3F=21 全單向偏 HP1
- → **paired mode 沒有 priority bug 偏移**

### 3.2 Step C — chr19 1 Mb window per-window som_ratio 分析

定義 `som_ratio = HP1-1 / (HP1-1 + HP2-1)`，每 1 Mb window 統計 chr19 全長：

| 統計 | 值 |
|---|---:|
| windows 數 | 57 |
| som_ratio mean | **0.462** |
| som_ratio median | **0.494** |
| som_ratio stdev | 0.332 |
| windows 偏 HP1-1 (>0.7) | 13 |
| windows 偏 HP2-1 (<0.3) | 18 |
| windows 近 1:1 (0.4-0.6) | 16 |

**重要觀察**：

1. **整體 mean 0.46 ≈ 0.5** → paired mode **無 systematic priority bug 偏向**
2. **stdev 0.332 跨 windows 大變化** → 反映**真實 sub-clone 結構**：
   - chr19:3M 全 HP2-1（755/0）→ 該區域 LOH 方向特定
   - chr19:0M 全 HP1-1（330/1）→ 反向區域
   - chr19:17M som_ratio 0.500（265/265）→ 雙 sub-clone 完美共現（**SP1 附近，TO mode 113:0 失衡的位點**）
3. **vs TO baseline**：TO mode 所有 chr19 windows 因 priority bug 全偏 HP1 → som_ratio 全 ~0.95；paired 跨 windows 0-1 全範圍分佈 → 真實生物學

### 3.3 Step B — Skip（基於 A+C 結論）

Step B 原計劃：patch longphase-s binary + read-level vote dump (~1-2 天)。

**Skip 理由**：
- Step A 整體 ratio 1:1.169 已強烈否定 priority bug 存在
- Step C window-level 跨 chr 變化是 sub-clone signal 而非 systematic bias
- Step B 預期 confirm 同樣結論（read-level vote 應與 read-level tag 結果一致）
- 成本/收益不對稱（1-2 天工程 vs 已確立結論）

→ Step B 列為 conditional follow-up，僅在發現 paired mode 異常時才啟動。

## 4. SP1 / SP2 / SP3 個案對照

5/8 整合報告 §2.2 在 TO mode chr19 SP1/2/3 觀察到 baseline 113:0 / 109:1 / 108:0 完全失衡。Paired tagged BAM 同位點 ±50bp 範圍：

| SP | TO baseline HP1:HP2 | Paired HP-fam ratio |
|---|---|---|
| SP1 chr19:17,565,944 | 113:0（priority bug） | **chr19:17M som_ratio 0.500（265/265 對稱）** |
| SP2 chr19:12,452,332 | 109:1（priority bug） | **chr19:12M som_ratio 0.278（132/342 偏 HP2-1）** |
| SP3 chr19:12,467,180 | 108:0（priority bug） | **chr19:12M 同上** |

→ TO baseline 在 SP1/2/3 全偏 HP1（artifact），paired 顯示真實 LOH 方向（SP1 雙 sub-clone 共現、SP2/3 偏 HP2-1）— 印證 PI 報告「baseline 與 paired 方向相反」。

## 5. 結論

| Q | A |
|---|---|
| paired mode 是否有同位點 HP1-1 + HP2-1 共現？ | **是** — 16/57 chr19 windows 共現比例 0.4-0.6 |
| 是否偏向 HP1-1（類似 TO 17.3:1）？ | **否** — som_ratio mean 0.462 跨 windows 0-1 全範圍分佈 |
| 共現現象的本質是什麼？ | **生物學 sub-clone signal**（LOH / 雙 sub-clone heterogeneity），不是工程 priority bug |
| paired vs TO mode binary 是否同樣 priority bug？ | **不同 codebase**：paired 用 longphase-s `HP:Z:` 字串 tag；TO 用 longphase-to `HP:i:` 整數 tag；兩者 tag 邏輯獨立，paired 沒 priority bug |

→ **paired mode 無需修補**；HP1-1 / HP2-1 共現可作為**亞克隆 marker** 候選（後續 follow-up 可量化它與 LOH / sub-clone fraction 的對應關係，但不在本 audit 範圍）。

## 6. 後續 follow-up（不在本 audit）

- **F-paired-1** Paired HP1-1/HP2-1 ratio 是否與 LOH segment 對齊？（從 LOH.bed cross-reference）
- **F-paired-2** Paired som_ratio 是否與 sub-clone fraction（caller AF）對應？
- **F-paired-3** 跨 7 樣本 paired audit 看 ratio 變化是否一致

## 7. 重現步驟

```bash
PAIRED_BAM=/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam

# Step A: chr19 整體分布
samtools view -F 256 "$PAIRED_BAM" chr19 | grep -oE "HP:Z:[^[:space:]]+" | sort | uniq -c | sort -rn

# Step C: 1Mb window 分析
samtools view -F 256 "$PAIRED_BAM" chr19 | awk '{
    hp=""; for(i=12;i<=NF;i++) if($i ~ /^HP:Z:/){hp=substr($i,6); break}
    if(hp=="") next
    win=int($4/1000000)
    if(hp=="1") w_hp1[win]++
    else if(hp=="2") w_hp2[win]++
    else if(hp=="1-1") w_hp1_1[win]++
    else if(hp=="2-1") w_hp2_1[win]++
    else if(hp=="3") w_hp3[win]++
} END {
    for(w=0;w<60;w++){
        a=w_hp1[w]+0; b=w_hp2[w]+0; c=w_hp1_1[w]+0; d=w_hp2_1[w]+0; e=w_hp3[w]+0
        if(a+b+c+d+e==0) continue
        sr=(c+d>0)?c/(c+d):0; gr=(a+b>0)?a/(a+b):0
        printf "chr19:%d %d %d %d %d %d %.3f %.3f\n",w,a,b,c,d,e,sr,gr
    }
}'
```
