<!--
建立時間: 2026-07-17
目標: 在論文投影片引用前，從 HCC1395 ONT BAM 可重現地驗證 primary-read N50
處理範圍: 投影片 headline 所指 HCC1395 5 kHz tumor BAM 全檔；paired normal／Dorado 僅作可選敏感度，不混入主 claim
狀態: verdict_GO_VALIDATED
audit_version: 0.1
關聯檔案:
  - InterSubMod/tools/bam_primary_read_n50.cpp
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260717_bam_read_n50_validation/
-->

# Pre-Decision Audit：HCC1395 ONT BAM primary-read N50

> **任務類型**：B — Comprehensive validation。  
> **服務目標**：G3（read-level evidence）＋G5（可外部重算的論文數據）。

## §0 Cynefin domain

**Complicated**。N50 定義與計算是確定性的；主要風險不是未知生物機制，而是 record 去重、輸入版本與 denominator 語意。

## §1 Observation completeness

| Observation | 狀態 | Tier | 來源 |
|---|---:|---:|---|
| 5 kHz tumor／normal 與 Dorado tumor／normal BAM 路徑已定位 | ✓ | L1 | 實際檔案＋symlink readback |
| 四個候選 BAM 通過 `samtools quickcheck -v` | ✓ | L1 | 2026-07-17 shell execution |
| Repo／KB 沒有可引用的 sample-specific N50 | ✓ | L1 | `rg`＋KB 查核 |
| Primary record 定義為排除 `0x100/0x800` | ✓ | L1 | SAM flag contract |
| 1 Mb pilot 的 custom tool 與獨立 `samtools view | awk` N50 完全一致 | ✓ | L1 | records=12,314；bases=119,466,726；N50=21,039 bp |
| 投影片 headline 的 HCC1395 5 kHz tumor BAM 全檔結果 | ✓ | L1 | exit 0；N50=21,047 bp；30,629,994 retained records |

## §2 Credibility score

| 維度 | 分數 | 理由 |
|---|---:|---|
| 理論基礎 | 20 | 標準 weighted read-length N50 |
| 觀察支撐 | 20 | Pilot 三路徑一致；全 BAM exit 0 並保存結果 checksum |
| 機制清晰度 | 20 | BAM record → sequence length → histogram → cumulative bases |
| 反例風險 | 20 | secondary／supplementary 重複與 paired-end 語意皆已預先處理 |
| 所需資源 | 10 | 約 1–2 小時只讀掃描；不產生大型中間檔 |
| **TOTAL** | **90 / 100** | **GO / VALIDATED** |

**Falsifier observable**：若 custom tool 與獨立 SAM text extraction 對相同 BAM 得到不同 records、total bases 或 N50，或全檔讀取非 exit 0，則結果不得進投影片。

**Reality-test 三反例**：

1. 若未排 supplementary，record count 與長度會被重複計算。
2. 若用 mean 代替 N50，長尾 ultra-long reads 會改變統計語意。
3. 若把 `N50/300` 說成連續 2×150 read 倍數，會錯誤忽略 paired-end insert gap。

## §3 Assumption map

| Importance | Known | Unknown／處理 |
|---|---|---|
| High | `l_qseq` 是 BAM 中保存的 query sequence length | BAM-based N50 不等於 raw POD5／FASTQ library N50；報告中明示 |
| High | Primary record 每 read 至多保留一個主 alignment | 若來源 pipeline 有非標準 hard-clipped primary，需另作 FASTQ sensitivity；本輪不外推 raw-library claim |
| Low | ONT duplication flag 對 N50 影響預期小 | 保留 duplicate primary records，因目標是 BAM read distribution；另報 filter 定義 |

## §4 Quick pilot

1. 擷取 HCC1395 5 kHz tumor `chr20:1-1,000,000`。  
   → 驗證：custom tool `records_used=12,314`、`total_bases=119,466,726`。
2. 以 `samtools stats -F 0x900` 比對 count、total length、mean 與 max。  
   → 驗證：四欄一致。
3. 以 `samtools view -F 0x900 | awk` 獨立重算 N50。  
   → 驗證：兩路徑皆為 `21,039 bp`。

Checkpoint：pilot 全 PASS → 允許全檔只讀掃描。

## §5 Gap diagnosis

| Missing | Impact | Effort | Priority |
|---|---:|---:|---:|
| Raw FASTQ／POD5 N50 | 只影響「原始 library N50」claim | 資料目前未定位 | P2 |
| Headline BAM 全檔結果 | 決定投影片精確數字 | 約 40 分鐘 | P0 |
| Paired normal／Dorado 全檔敏感度 | 只影響跨資料集泛化；不影響 sample-specific headline | 額外 1–2 小時、共享儲存負載高 | P2 |
| Illumina matched-library empirical N50 | 不影響 150 bp protocol 比較 | 需另取短讀資料 | P2 |

## §6 Conflict scan

| Prior statement | Relation | 處理 |
|---|---|---|
| 一般 ONT read N50 常大於 10 kb | 非本資料證據 | 不引用為本研究數值 |
| 原投影片寫 10–50 kb、100–500 倍 | 未驗證範圍 | 由實測 sample-specific N50 取代 |
| ClairS Extended Data Fig. 7 有 read-length distributions | 支持資料背景，但無文字 N50 | 只作背景來源，不代替本次計算 |

## §7 Red team 與 verdict

最強反方：

1. Aligned BAM 可能排除部分 raw reads，因此不能把結果稱為原始測序 library N50。
2. 不同 sample、chemistry、basecaller 與 DNA extraction 的 N50 可不同，單一數字不能代表所有 ONT。

處理：投影片主句固定寫 `HCC1395 5 kHz tumor BAM primary-read N50`。本輪完整掃描主句所指的整個 BAM；paired normal／Dorado 不混入 headline，未完成全檔敏感度時也不作跨資料集泛化。

**Verdict：GO / VALIDATED（90/100）**。全檔 exit 0；輸入、filter、runtime、結果與 checksum 均已保存。投影片 headline 可填 `21.0 kb；約 140×單端150 bp`。
