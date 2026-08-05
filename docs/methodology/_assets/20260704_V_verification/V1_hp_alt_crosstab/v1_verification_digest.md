# V1 — HP tag × ALT-support 交叉表(unphase / HP3 operational 定義驗證)

**軌道**: V1_hp_alt_crosstab · **狀態**: ✅ completed(genome-wide 已驗證)
**日期**: 2026-07-04 · **樣本**: HCC1395(單樣本 ⭐3 上限)

## 1. 資料與方法(每個數字皆可追檔)

| 項目 | 值 / 路徑 |
|---|---|
| Tagged tumor BAM | `/big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam` (283 GB, → `HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam`) |
| Somatic sSNV(TP) | `data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz` — **30,490** sites |
| Somatic sSNV(FP) | `data/vcf/HCC1395/pileup/filtered_snv_fp.vcf.gz` — **4,842** sites |
| 合併分析位點 | **35,332** SNV positions(chr1–22;VCF 只含 autosome) |
| Read filter | MAPQ≥20,排除 unmapped/secondary/supplementary(同 ReadParser + ref script) |
| ALT 判定 | **忠實複刻 `ReadParser::determine_alt_support_with_reason`**:走 CIGAR 取 SNV 位鹼基;`==ALT→ALT`、`==REF→REF`、其他→UNKNOWN;BQ<20 或落在 deletion/refskip → UNKNOWN(`min_base_quality=20`) |
| HP 映射 | `1→1, 2→2, 11→1-1, 21→2-1, 33→3`;無 HP tag → `0`(unphase)。同 `ReadParser.cpp:133–152` |
| 覆蓋讀取數(union TP+FP) | **2,961,158** reads |
| Read 分類 | 每條 read 對其覆蓋的所有 somatic 位點:**ALT-touching**(≥1 個 ALT)/ **REF-only**(有 confident REF、0 ALT)/ **UNKNOWN-only**(有覆蓋但無 confident REF/ALT) |

**位點覆蓋完整性**:22 條染色體每一條 `positions_covered == positions_total`(35,332/35,332),無漏掃。

## 2. Genome-wide 交叉表(union TP+FP,已驗證)

| HP tag | family | reads (%) | ALT-touching | REF-only | UNKNOWN-only | **%ALT/covered** | **%ALT/confident** |
|---|---|---:|---:|---:|---:|---:|---:|
| `0` | unphase | 296,674 (10.0%) | 2,331 | 232,968 | 61,375 | **0.79%** | 0.99% |
| `1` | germline HP1 | 672,491 (22.7%) | **0** | 533,020 | 139,471 | **0.0%** | 0.0% |
| `2` | germline HP2 | 672,217 (22.7%) | **0** | 534,302 | 137,915 | **0.0%** | 0.0% |
| `1-1` | HP1-1 | 629,591 (21.3%) | 525,400 | 15,762 | 88,429 | **83.5%** | **97.1%** |
| `2-1` | HP2-1 | 615,402 (20.8%) | 513,104 | 15,738 | 86,560 | **83.4%** | **97.0%** |
| `3` | HP3 | 74,783 (2.5%) | 62,379 | 1,591 | 10,813 | **83.4%** | **97.5%** |

*%ALT/covered = ALT-touching ÷ 全部覆蓋讀取;%ALT/confident = ALT ÷ (ALT+REF),排除 UNKNOWN-only。*
*完整輸出:`v1_hp_alt_crosstab_ALL.json`、`v1_hp_alt_crosstab_ALL_crosstab.tsv`。*

## 3. 對 PI 三項斷言的裁決

### ✅ 斷言 1:unphase(tag 0/空)≈100% REF-only,不 touch 任何 ALT — **驗證通過**
- unphase 讀取中 **0.79%** touch ALT(其餘覆蓋讀取中 REF-only=**78.5%**、UNKNOWN-only=**20.7%**),confident 判定中僅 **0.99%** 為 ALT。
- 即 **99.2%** 的 unphase 讀取不碰任何 somatic ALT,與斷言一致。
- ⚠️ **非字面 100%**:仍有 2,331 條 unphase 讀取 touch ALT(其中 2,269 條 touch 的是 **TP** 真 somatic 位點)。這批是「過了 ALT 卻仍被判 unphase」的少數例外 → 與 PI「unphase=沒碰 germline 或 germline 衝突」的定義相容(unphase 的判別軸是 germline,不是 ALT),但說明 unphase **不是**純粹「只過 REF」。

### ✅ 斷言 2:HP1-1/HP2-1 都經過 somatic;HP3 都經過 ALT — **方向驗證通過,字面「都」需修正為「壓倒性多數」**
- confident 判定中 ALT 佔比:HP1-1 **97.1%**、HP2-1 **97.0%**、HP3 **97.5%** — 三者一致 ~97%,遠高於 germline 的 0%。
- ⚠️ **非字面 100%**:三家族仍各有 confident **REF-only** 讀取(HP1-1=15,762、HP2-1=15,738、HP3=1,591),以及大量 **UNKNOWN-only**(88,429/86,560/10,813;佔各家族 ~14%)。UNKNOWN-only = 讀取有覆蓋 somatic 位點但該位點 BQ<20 / 落在 indel / 非 REF 非 ALT。
- 判別軸「是否經過 ALT」在 **confident** 層級成立(97% vs 0%);但「HP3=**都** touch ALT」若當硬性 operational 過濾條件會誤殺 ~2.5% confident-REF-only + ~14% UNKNOWN-only 讀取。

### ✅ 斷言 3:germline HP1/HP2 的 REF/ALT 分布 — **驗證通過(強結果)**
- germline HP1/HP2 **ALT-touching = 0 條**(genome-wide,1,344,708 條讀取,exactly 0.0%)。
- 完全符合 LongPhase-S 設計:germline haplotype (1/2) 是「純 germline phase、未碰任何 somatic ALT」的讀取;碰了 somatic ALT 的同 phase 讀取會被升格為 1-1/2-1。
- **cross-HP germline 衝突比例**:本交叉表以 raw HP tag 分類,germline 1 與 2 互斥(一條 read 只有一個 HP tag),故無「跨 HP」讀取;germline 衝突發生在 *tagging 上游*(LongPhase-S 決定不給 1/2 或給 0),本 BAM 中呈現為 unphase(0) 類別。

## 4. 判別軸總結(一句話)

**「是否經過(confident support)任一 somatic ALT」是 HP 家族的乾淨判別軸**:germline(1/2)= 0% ALT、somatic(1-1/2-1/3)= 97% ALT(confident)。unphase(0)靠 germline 軸定義,附帶效果是 99.2% 不碰 ALT。PI 修正方向正確。

## 5. P0 文件修正建議

**缺陷(原始碼佐證)**:`include/core/DataStructs.hpp:32` 註解寫
`hp_tag: "1", "2", "1-1", "2-1", "unphase", etc.` — 但 `ReadParser.cpp:125,159` 對「無 HP tag / 降級」一律設 **`"0"`**,**從不產生字面字串 `"unphase"`**。
- 下游比對確實接受 `"unphased"`(`LabelTest.cpp:459,772`:`hp=="0" || hp=="HP0" || hp=="unphased"`),但 ReadParser 只會輸出 `"0"`,`"unphase"`/`"unphased"` 是 dead branch。

**建議修正**:
1. `DataStructs.hpp:32` 註解把 `"unphase"` 改為 **`"0"`(= unphase/unphased,無 HP tag 或 germline-only 降級)**,與實際輸出一致。
2. 文件(HP tag 語意表)明確:**unphase 的判別軸是 germline,不是 ALT**;實測 unphase 有 0.79%(genome-wide 2,331 條)讀取仍 touch ALT,故不可寫「unphase = 只過 REF」的硬定義,應寫「unphase = germline phase 無法確定(未碰 germline 或 germline 衝突);絕大多數(99.2%)不碰 somatic ALT」。
3. HP3 / 1-1 / 2-1 定義文件把「**都**經過 ALT」改為「**經過 confident ALT support(~97% of confident calls;family 內另有 ~14% 讀取在 somatic 位點無 confident 判定)**」,避免當硬過濾條件時誤殺。

## 6. 管線驗證軌跡(chr22 pilot → genome-wide)

- 先以 **chr22**(484 sites, 902,335 reads)驗管線:unphase 0.66% ALT、germline 0.0%、HP1-1/2-1/HP3 ~84–89% ALT — 與 genome-wide 同向(chr22 為 pilot,非最終數字)。
- 再全掃 **chr1–22**(35,332 sites, 2.96M reads, 16 workers, ~35 min);本 digest 所有數字取自 genome-wide 輸出。

## 7. 檔案清單(track dir: `docs/methodology/_assets/20260704_V_verification/V1_hp_alt_crosstab/`)

- `v1_hp_alt_crosstab.py` — 分析腳本(CIGAR-walk, multiprocessing)
- `v1_hp_alt_crosstab_ALL.json` / `_crosstab.tsv` — genome-wide 交叉表(union + TP-only view)
- `val_chr22.json` / `_crosstab.tsv` — chr22 pilot
- `v1_hp_alt_crosstab.png` — 圖(本 digest)
