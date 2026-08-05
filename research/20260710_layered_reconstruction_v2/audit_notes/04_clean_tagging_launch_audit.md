<!--
建立時間：2026-07-11（Asia/Taipei）
目標：審計 LongPhase-S 7-dataset clean production tagging 的程式、輸入、命令、命名、resume/skip 與驗收契約
處理範圍：scripts/pipeline/steps/01_longphase_s.sh、scripts/pipeline/run_benchmark.sh、7 個 paired_full canonical roots、LongPhase-S v1.0.0 KB/runtime/source、既有 production sidecar collaborator run
關聯檔案：InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/02_layered_chain_audit.md；InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/layered_v2_input_manifest.json
限制：本文件只做唯讀審計與 launch design；未修改 pipeline/source，未啟動、重啟或終止任何 LongPhase-S 長計算
-->

# 7-dataset clean LongPhase-S production tagging 啟動審計

> **用 SCQA + Gate review：可定義一條不帶 truth、角色不混淆的 tag-only 命令；但 7 份 persistent BAM 現在因容量不足必須 HOLD，既有 2-sample FIFO sidecar run 僅 attach，不得重啟（影響：高；信心：高）。**

## 0. 結論先行

### 0.1 Verdict

| Gate | Verdict | 理由 |
|---|---|---|
| 7-dataset 輸入可用性 | **PASS** | 7/7 tumor/normal BAM `samtools quickcheck` PASS、BAI 存在；raw VCF、germline phased VCF 及 index 存在；24 primary contigs 字典一致。 |
| caller backbone 身分 | **PASS，但需重新命名鎖定** | 7 組 canonical `somatic_pass.vcf.gz` record body 與各 raw ClairS VCF 的 `FILTER=PASS` stream SHA-256 完全相同；但現行 step 的同名檔案有角色覆寫風險，正式 run 應建立 `*.caller_pass.snv.vcf.gz` snapshot。 |
| 無 truth 的 genome-wide tagging | **PASS（source/runtime 證實）** | `--truth-vcf`、`--truth-bed` 都是 optional；不傳時不進 benchmark/HC BED removal。 |
| caller-PASS-preserving tag-only | **PASS，命令必須帶 `--disableFilter`** | LongPhase-S 預設會做內部 variant filtering；僅移除 `--output-somatic-vcf` 不會停用它。`--disableFilter` 才會接受所有供應的 caller-PASS variants。 |
| smoke probe | **PROBE** | 建議 HCC1395 `chr1:1-5000000`；實測輸入有 37 個 caller-PASS SNV，足以驗證 CLI、`@PG`、HP vocabulary、BAI 與 record-count contract。此 probe 尚未由本審計執行。 |
| 7 份 persistent BAM full launch | **HOLD / capacity NO-GO** | 歷史輸出合計 1.841 TB；含 20% headroom 至少需約 2.21 TB，建議 hard gate 2.3 TB free。2026-07-11 實測 `/big7_disk` 僅 932 G free、98% used。 |
| 既有 collaborator sidecar run | **ATTACH-ONLY** | `20260711_longphase_s_production_sidecars_v1` 已有 HCC1395/HCC1395_DORADO START；其 BAM 是 FIFO 且 payload 不保留，不能當作 persistent BAM/BAI 完成證據，也不得由本流程重啟。 |

**本輪沒有啟動任何長計算。** 可執行順序只能是：容量解除 → bounded smoke PASS → full sample-by-sample launch → per-sample acceptance → 7/7 completeness → 才更新 layered manifest 或 promotion。

### 0.2 任務分類與研究目標

- Task type：**(C) Production deployment launch audit（design-only）**。
- Scope：7 datasets、24 primary contigs genome-wide tagging；下游 layered reconstruction 仍只消費 chr1–22。
- 服務目標：**G3**（read-level epigenetic/phasing artifact）、**G4**（7-dataset reproducibility）、**G5**（可外部驗證的 production evidence chain）。
- 不是 Thread B 撤回方向；KDE、VCF caller AF 不適用；涉及長計算與大型 BAM，因此必須有 launch/capacity gate。

## 1. 角色先分清楚：禁止再使用一個 `somatic_pass.vcf.gz` 表示兩件事

| Canonical role | 建議名稱 | 定義 | 可否餵給本次 tag-only |
|---|---|---|---|
| `caller_raw_vcf` | 原 ClairS `output.vcf.gz` | caller 原始輸出，含 PASS/LowQual | 不能直接；先抽 `FILTER=PASS` |
| `caller_pass_input` | `<sample>.clairs_paired.caller_pass.snv.vcf.gz` | raw ClairS 中 `FILTER=PASS && TYPE=SNV` 的 immutable snapshot | **是；唯一 somatic input** |
| `lps_tag_only_bam` | `<sample>.lps_tag_only.bam` | 無 truth、`--disableFilter`、無 `--output-somatic-vcf` 的 persistent tagged BAM | 本次目標 |
| `lps_dna_fraction` | `<sample>.lps_tag_only_purity.out` | CLI 名稱仍為 purity；概念應寫 tumor DNA fraction | 診斷 sidecar，不是 cellular purity |
| `lps_recalibrated_vcf` | `<sample>.lps_recalibrated.all.vcf.gz` | 預設 LongPhase-S filter + `--output-somatic-vcf` 的 PASS/LowQual recalibration | **另一個實驗角色，不得冒充 caller PASS** |
| `benchmark_truth_conditioned_bam` | `<sample>.legacy_truth_conditioned.bam` | 帶 `--truth-vcf`，多數還帶 `--truth-bed` 的歷史 BAM | 不可作 clean production backbone |
| `streamed_read_tag_sidecar` | `<sample>.read_tags.tsv.gz` | 從 FIFO BAM 擷取 HP/PS；沒有 persistent BAM/BAI | 可供下游 sidecar 路徑，但不是本次 BAM deliverable |

本次「clean」有四個同時成立的條件：

1. `--truth-vcf` **不存在**。
2. `--truth-bed` **不存在**。
3. `--output-somatic-vcf` **不存在**，避免把 tagging 與 VCF recalibration 混成同一 role。
4. `--disableFilter` **存在**，因為輸入已是 caller PASS；否則 LongPhase-S 仍會在記憶體中二次決定哪些 variant 可參與 somatic tagging。

## 2. 程式碼 control flow 審計

### 2.1 `01_longphase_s.sh` 不能直接拿來跑本次 production tag-only

權威檔案：`InterSubMod/scripts/pipeline/steps/01_longphase_s.sh`。

| Code range | 實際行為 | Production 風險 |
|---|---|---|
| 81–99 | 輸出 truth path，並強制 `validate_file "${TRUTH_VCF}"` | 無 truth 的 production mode 無法從此入口表達。 |
| 124–152 | 從 caller VCF 抽 `FILTER=PASS` 到 `longphase_s/somatic_pass.vcf.gz` | 此刻檔案 role 是 caller PASS。 |
| 155–172 | 固定加入 `--output-somatic-vcf --truth-vcf`，有 BED 再加 `--truth-bed` | tagging、recalibration、benchmark 三個角色耦合。 |
| 192–209 | 產生/索引 `<sample>_tagged.bam`，跑 haplotag QC | 這一段本身合理，但前述命令已污染語意。 |
| 220–246 | 再從 `<prefix>_sc.vcf` 抽 PASS，覆寫同一路徑 `somatic_pass.vcf.gz` | 同名檔案由 caller PASS 變成 LongPhase recalibrated PASS。 |
| 248–280 | 用 truth VCF/BED 做 `bcftools isec` 產 TP/FP | benchmark-only，不屬 production tagging。 |
| 289–291 | 直接 `rm` isec temp 與 `somatic_pass.vcf.gz` | 與 repo 不刪檔 governance 衝突，也破壞 exact upstream 可重現性。 |

附帶問題：line 183 用 `/usr/bin/time ... | tee`；因全檔有 `set -euo pipefail`，LongPhase nonzero 通常可傳遞，但沒有把 LongPhase exit、index exit、acceptance exit 分別寫入結構化 status。

### 2.2 `run_benchmark.sh` 沒有 clean-tagging 或真正 resume mode

權威檔案：`InterSubMod/scripts/pipeline/run_benchmark.sh`。

- line 110–123 每次建立新 canonical run id；已存在就加 `_1`, `_2`。這是 **avoid overwrite**，不是 resume。
- line 205–235 的正常路徑固定呼叫上述 `01_longphase_s.sh`，所以仍會帶 truth/recalibration。
- `--skip-longphase`（165–204）同時要求 tagged BAM、TP VCF、FP VCF；它是「以既有 benchmark artifacts 繼續 InterSubMod」，不是驗收 clean production BAM。
- `--skip-longphase` 找不到明示輸入時用「lexically latest canonical directory」，若存在不同 role 的 run，可能拿錯 artifact。
- pipeline 沒有 `ACCEPTED` sentinel、輸入 hash match、partial BAM quarantine 或 per-sample checkpoint；LongPhase-S CLI help 也沒有 checkpoint/resume option。

**判定**：本次不可用 `run_benchmark.sh` 啟動。應先用 LongPhase-S 原生 CLI 產生一個獨立、role-specific 的 production bundle；通過完整驗收後，後續流程才可用明示 `--tagged-bam` 消費，且 TP/FP benchmark 另行處理。

### 2.3 LongPhase-S source 對參數的真正語意

本機 source：`/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s`；commit `420e0ff73188126da8bb2122725fbffb2d076ccc`，工作樹 clean。binary：

```text
/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/longphase-s
version: 1.0.0
sha256: 07cbd0aa0c9f33ed59c5baff45fbe3554ef96d55457635de7348c4501b283f54
```

| Source | 判斷 |
|---|---|
| `src/somatic_haplotag/SomaticHaplotag.cpp:58–68` | filter 預設 enabled；automatic DNA-fraction estimation 預設 enabled；`--output-somatic-vcf` 預設 false。 |
| 同檔 100–108 | truth VCF/BED 用 optional validator；不是 required input。 |
| `SomaticHaplotagProcess.cpp:70–78` | **不論是否輸出 `_sc.vcf`，variantCalling/filtering 都會執行**。 |
| 同檔 80–88 | `--output-somatic-vcf` 只控制是否把 recalibrated VCF 寫出。 |
| 同檔 90–101 | 有 truth VCF 且 truth BED 時，先移除 BED 外 tumor/truth variants，再呼叫 `tagRead()`；故 truth BED 會改變 tagging universe，不只是報表範圍。 |
| `SomaticVarCaller.cpp:1206–1225` | variants 可被多種規則標 `isFilterOut`；只有 `enableFilter=false` 時不會 continue，並把所有供應 variants 標為 high-confidence somatic。 |
| 同檔 1366–1405 | read HP calibration 依 `isHighConSomaticSNP`；所以 `--disableFilter` 確實讓 caller-PASS variants 保持可參與。 |
| `HaplotagParsingBam.cpp:40–49` | 讀入 tumor BAM header，新增 `@PG ID:longphase-s VN:<version> CL:<完整命令>`，並要求 BAM index。 |
| 同檔 56–70 | `-o PREFIX` 寫出 `PREFIX.bam`；LongPhase-S 本身不建立 BAI。 |

### 2.4 基本判定與設定，不保留隱含預設

| 設定 | 本次固定值 | 最基本語意 |
|---|---:|---|
| `-q` | 20 | MAPQ <20 alignment 保留在 BAM，但不給新的 haplotag；不是 variant QUAL。 |
| `-p` | 0.6 | read 所見 allele 中，主 haplotype 比例未達 0.6 就不做該 haplotype assignment。明寫避免版本預設漂移。 |
| `--tagSupplementary` | on | supplementary alignment 也進 tagging。 |
| `-t` | 12 | 與現有 7-sample sidecar run 一致；歷史 CPU utilization 顯示主要瓶頸不是純 CPU。 |
| `--disableFilter` | on | 接受所有已供應 caller-PASS SNV；不代表停用 MAPQ 或 `-p` read assignment rules。 |
| `--output-somatic-vcf` | off | 不產生 `_sc.vcf`；避免 tag-only 與 recalibration role 混淆。 |
| `--truth-vcf/--truth-bed` | absent | 不進 benchmark、不用 HC BED 裁掉 tagging variants。 |
| `--region` | absent（full） | production BAM 全 genome；smoke 才有 region。 |
| `--tumor-purity` | absent | 保留 automatic tumor DNA fraction estimation作診斷；因 filter disabled，不用它刪 caller-PASS variants。不可把結果稱為 cellular purity。 |
| output format | BAM | 本次 deliverable 要 BAM + BAI；CRAM 雖省空間，但屬另一個 output contract。 |

## 3. 七個 dataset 的實際輸入鎖定

共同 reference：

```text
/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta
size_bytes=3144230986
FAI=/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta.fai
FAI_sha256=502a1b8fb73ccd53285c28a0f12df90c818b4fe3de1e862ef47c593ef1a0a4b4
```

### 3.1 BAM 與 caller raw VCF

| Dataset | Tumor BAM（config path；target size） | Normal BAM（config path；target size） | ClairS raw VCF |
|---|---|---|---|
| HCC1395 | `/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam`；292.1 GB | `/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam`；149.4 GB | `/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/output.vcf.gz` |
| HCC1395_DORADO | `/big8_disk/data/HCC1395/ONT_Dorado/HCC1395.bam` → `/big8_disk/Google_somatic_data/bams/HCC1395/HCC1395_Tumor_ONT.GRCh38.sorted.bam`；250.2 GB | `/big8_disk/data/HCC1395/ONT_Dorado/HCC1395BL.bam` → `/big8_disk/Google_somatic_data/bams/HCC1395/HCC1395_Normal_ONT.GRCh38.sorted.bam`；91.4 GB | `/big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/output.vcf.gz` |
| COLO829 | `/big8_disk/data/COLO829/ONT_PAO/PAO29420.bam`；100.4 GB | `/big8_disk/data/COLO829/ONT_PAO/PAO33946.bam`；124.7 GB | `/big8_disk/data/COLO829/ONT_PAO/ClairS_v0_4_1/output.vcf.gz` |
| H1437 | `/big8_disk/data/H1437/ONT/H1437.bam` → `/big8_disk/Google_somatic_data/bams/H1437/H1437_Tumor_ONT.GRCh38.sorted.bam`；243.7 GB | `/big8_disk/data/H1437/ONT/H1437BL.bam` → `/big8_disk/Google_somatic_data/bams/H1437/H1437_Normal_ONT.GRCh38.sorted.bam`；149.6 GB | `/big8_disk/data/H1437/ONT/ClairS_v0_4_1/output.vcf.gz` |
| H2009 | `/big8_disk/data/H2009/ONT/H2009.bam` → `/big8_disk/Google_somatic_data/bams/H2009/H2009_Tumor_ONT.GRCh38.sorted.bam`；327.7 GB | `/big8_disk/data/H2009/ONT/H2009BL.bam` → `/big8_disk/Google_somatic_data/bams/H2009/H2009_Normal_ONT.GRCh38.sorted.bam`；92.6 GB | `/big8_disk/data/H2009/ONT/ClairS_v0_4_1/output.vcf.gz` |
| HCC1937 | `/big8_disk/data/HCC1937/ONT/HCC1937.bam` → `/big8_disk/Google_somatic_data/bams/HCC1937/HCC1937_Tumor_ONT.GRCh38.sorted.bam`；472.1 GB | `/big8_disk/data/HCC1937/ONT/HCC1937BL.bam` → `/big8_disk/Google_somatic_data/bams/HCC1937/HCC1937_Normal_ONT.GRCh38.sorted.bam`；80.2 GB | `/big8_disk/data/HCC1937/ONT/ClairS_v0_4_1/output.vcf.gz` |
| HCC1954 | `/big8_disk/data/HCC1954/ONT/HCC1954.bam` → `/big8_disk/Google_somatic_data/bams/HCC1954/HCC1954_Tumor_ONT.GRCh38.sorted.bam`；253.2 GB | `/big8_disk/data/HCC1954/ONT/HCC1954BL.bam` → `/big8_disk/Google_somatic_data/bams/HCC1954/HCC1954_Normal_ONT.GRCh38.sorted.bam`；105.1 GB | `/big8_disk/data/HCC1954/ONT/ClairS_v0_4_1/output.vcf.gz` |

所有 BAM config path 與 resolved target 都必須寫進 `input_files.tsv`；使用 `stat -Lc`，不可用不 follow symlink 的 `stat -c`。

### 3.2 Phased germline VCF 與 caller-PASS 數

| Dataset | Germline phased VCF | records / phased+PS | caller PASS SNV | smoke chr1:1–5 Mb |
|---|---|---:|---:|---:|
| HCC1395 | `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz` | 1,698,116 / 1,687,835 | 113,997 | 37 |
| HCC1395_DORADO | `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395_DORADO/paired_full/20260315_HCC1395_DORADO_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz` | 1,681,556 / 1,664,980 | 112,387 | 38 |
| COLO829 | `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/COLO829/paired_full/20260315_COLO829_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz` | 1,769,629 / 1,761,450 | 38,196 | 41 |
| H1437 | `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/H1437/paired_full/20260315_H1437_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz` | 1,701,898 / 1,694,629 | 75,578 | 34 |
| H2009 | `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/H2009/paired_full/20260315_H2009_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz` | 1,783,991 / 1,774,501 | 157,405 | 103 |
| HCC1937 | `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1937/paired_full/20260315_HCC1937_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz` | 1,780,458 / 1,769,002 | 49,548 | 108 |
| HCC1954 | `/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz` | 1,850,528 / 1,841,976 | 20,969 | 37 |

QA 實測：

- 7/7 tumor/normal BAM：`samtools quickcheck -v` exit 0，`.bai` 存在。
- 7/7 raw VCF：`.tbi` 存在；7/7 germline VCF：`.csi` 存在。
- reference、tumor BAM、normal BAM、raw VCF、germline VCF 的 chr1–22/X/Y 名稱與長度：7/7 全部一致。
- 7 組 caller PASS 全是 biallelic SNV；PASS count 如表。
- 對每一組做 `bcftools view -f PASS -H raw | sha256sum` 與 canonical `somatic_pass.vcf.gz` record stream hash：7/7 相同。因此現有 canonical 檔案內容可被證實是 caller PASS；仍建議由 raw 重建 role-specific snapshot，避免未來被同名覆寫。
- 大型 BAM 本輪未做 full-content SHA-256；目前 identity strength 是 resolved path + inode/size/mtime + BAI checksum + header dictionary + quickcheck。若要對外 handoff，應另排 I/O window 產生 full checksum，不能假裝已完成。

## 4. 輸出 layout 與不可混淆命名

計畫中的 persistent run（**尚未建立**）：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/
└── 20260711_longphase_s_persistent_tag_only_no_truth_v1/
    ├── README.md
    ├── input_manifest.json
    ├── params.json
    ├── code.sha256
    ├── run_status.tsv
    ├── completeness.tsv
    ├── smoke/
    │   └── HCC1395/
    │       ├── HCC1395.lps_tag_only.chr1_1_5Mb.bam
    │       ├── HCC1395.lps_tag_only.chr1_1_5Mb.bam.bai
    │       ├── command.sh.txt
    │       ├── longphase_s.log
    │       └── acceptance.tsv
    └── samples/<SAMPLE>/
        ├── inputs/
        │   ├── <SAMPLE>.clairs_paired.caller_pass.snv.vcf.gz
        │   ├── <SAMPLE>.clairs_paired.caller_pass.snv.vcf.gz.csi
        │   ├── input_files.tsv
        │   └── input.sha256
        ├── tag_only/
        │   ├── <SAMPLE>.lps_tag_only.bam
        │   ├── <SAMPLE>.lps_tag_only.bam.bai
        │   └── <SAMPLE>.lps_tag_only_purity.out
        ├── command.sh.txt
        ├── longphase_s.log
        ├── resource.time.txt
        ├── acceptance.json
        ├── acceptance.tsv
        ├── output.sha256
        └── ACCEPTED
```

此 root 在 7/7 acceptance 前是 research round，不得直接宣稱 canonical。若要 canonical promotion，應建立新的 immutable run bundle及 completeness/provenance，而不是覆蓋 20260314/15 historical runs或在舊目錄換 BAM。

## 5. Caller-PASS snapshot 命令

每個 dataset 先由 §3.1 的 raw VCF 建立新 snapshot。以 HCC1395 為例：

```bash
RUN_ROOT=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_persistent_tag_only_no_truth_v1
SAMPLE=HCC1395
RAW=/big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/output.vcf.gz
INPUT_DIR="$RUN_ROOT/samples/$SAMPLE/inputs"
CALLER_PASS="$INPUT_DIR/$SAMPLE.clairs_paired.caller_pass.snv.vcf.gz"

mkdir -p "$INPUT_DIR"
/usr/bin/bcftools view -f PASS -v snps -Oz -o "$CALLER_PASS" "$RAW"
/usr/bin/bcftools index -c "$CALLER_PASS"
test "$(/usr/bin/bcftools view -H "$CALLER_PASS" | wc -l)" -eq 113997
test "$(/usr/bin/bcftools view -f PASS -v snps -H "$RAW" | sha256sum | cut -d' ' -f1)" = \
     "$(/usr/bin/bcftools view -H "$CALLER_PASS" | sha256sum | cut -d' ' -f1)"
```

其餘 expected count：HCC1395_DORADO 112387、COLO829 38196、H1437 75578、H2009 157405、HCC1937 49548、HCC1954 20969。任一 count/hash 不等即 exit nonzero；不能「繼續看看」。

## 6. Bounded smoke probe（PARTIAL / 非 validation evidence）

### 6.1 輸入、命令、輸出

- 輸入：HCC1395 的五個實際 input，somatic snapshot 應為 113,997 records。
- 限制區域：`chr1:1-5000000`，已驗證有 37 個 caller-PASS SNV。
- 輸出：`$RUN_ROOT/smoke/HCC1395/HCC1395.lps_tag_only.chr1_1_5Mb.bam[.bai]`。
- 此輸出必須明標 `PARTIAL_SMOKE_ONLY`，不可寫入 layered production manifest。

```bash
RUN_ROOT=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_persistent_tag_only_no_truth_v1
SMOKE="$RUN_ROOT/smoke/HCC1395"
mkdir -p "$SMOKE"

set -o pipefail
/usr/bin/time -v -o "$SMOKE/resource.time.txt" \
  /big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/longphase-s somatic_haplotag \
  -s /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz \
  -b /big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam \
  --tumor-snv-file "$RUN_ROOT/samples/HCC1395/inputs/HCC1395.clairs_paired.caller_pass.snv.vcf.gz" \
  --tumor-bam-file /big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam \
  -r /big8_disk/ref/GRCh38_no_alt_analysis_set.fasta \
  --region chr1:1-5000000 \
  -t 12 --tagSupplementary -q 20 -p 0.6 --disableFilter \
  -o "$SMOKE/HCC1395.lps_tag_only.chr1_1_5Mb" \
  2>&1 | tee "$SMOKE/longphase_s.log"
LPS_RC=${PIPESTATUS[0]}
test "$LPS_RC" -eq 0
/usr/local/bin/samtools index -@ 8 "$SMOKE/HCC1395.lps_tag_only.chr1_1_5Mb.bam"
```

### 6.2 Smoke acceptance

1. LongPhase exit 0 → 驗證：`LPS_RC=0`。
2. Log role 正確 → 驗證：truth VCF/BED 為空、`tag region=chr1:1-5000000`、`variant filtering=disabled`、`write calling VCF=disabled`。
3. BAM 是 regular file → 驗證：`test -f BAM && test ! -p BAM`。
4. BAM 完整 → 驗證：`samtools quickcheck -v BAM` 無輸出且 exit 0。
5. Index → 驗證：`samtools index` exit 0，`.bam.bai` 非空，`idxstats` exit 0。
6. `@PG` → 驗證：最新 `PN:longphase-s` CL 含 `--region ... --disableFilter -q 20 -p 0.6 --tagSupplementary`，且不含 `--truth-*`/`--output-somatic-vcf`。
7. HP vocabulary → 驗證：所有實際 `HP:Z` 值只在 `{1,2,3,4,1-1,2-1,1-2,2-2}`；沒有 HP tag代表 unassigned，不要求寫 `HP:Z:.`。
8. MM/ML preservation → 驗證：同一批 deterministic qname 的 input/output MM/ML presence 與值 hash 相同。

預期 log 片段：

```text
[Benchmark Files]
truth VCF file               :
truth BED file               :
tag region                   : chr1:1-5000000
variant filtering            : disabled
write calling VCF            : disabled
```

## 7. Full production 命令契約

### 7.1 唯一允許的 command shape

```bash
/usr/bin/time -v -o "$RESOURCE" \
  "$LPS" somatic_haplotag \
  -s "$GERMLINE" \
  -b "$NORMAL" \
  --tumor-snv-file "$CALLER_PASS" \
  --tumor-bam-file "$TUMOR" \
  -r /big8_disk/ref/GRCh38_no_alt_analysis_set.fasta \
  -t 12 --tagSupplementary -q 20 -p 0.6 --disableFilter \
  -o "$PREFIX" \
  2>&1 | tee "$LOG"
```

Full command 中禁止出現：`--truth-vcf`、`--truth-bed`、`--benchmark-log`、`--output-somatic-vcf`、`--somatic-calling-log`、`--region`。`PREFIX` 必須是 `<sample>/tag_only/<sample>.lps_tag_only`。

### 7.2 七組完整參數渲染

下面呼叫可直接餵給一個實作上述 command shape 的 `launch_one SAMPLE GERMLINE NORMAL RAW TUMOR EXPECTED` function；所有 input/output role 都已展開。這是 launch specification，**容量 gate 未解除前禁止執行**。

```bash
launch_one HCC1395 \
  /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz \
  /big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam \
  /big8_disk/data/HCC1395/ONT/ClairS_v0_4_0/output.vcf.gz \
  /big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam 113997

launch_one HCC1395_DORADO \
  /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395_DORADO/paired_full/20260315_HCC1395_DORADO_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz \
  /big8_disk/data/HCC1395/ONT_Dorado/HCC1395BL.bam \
  /big8_disk/data/HCC1395/ONT_Dorado/ClairS_v0_4_0/output.vcf.gz \
  /big8_disk/data/HCC1395/ONT_Dorado/HCC1395.bam 112387

launch_one COLO829 \
  /big7_disk/liaoyoyo2001/big7_disk_output/canonical/COLO829/paired_full/20260315_COLO829_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz \
  /big8_disk/data/COLO829/ONT_PAO/PAO33946.bam \
  /big8_disk/data/COLO829/ONT_PAO/ClairS_v0_4_1/output.vcf.gz \
  /big8_disk/data/COLO829/ONT_PAO/PAO29420.bam 38196

launch_one H1437 \
  /big7_disk/liaoyoyo2001/big7_disk_output/canonical/H1437/paired_full/20260315_H1437_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz \
  /big8_disk/data/H1437/ONT/H1437BL.bam \
  /big8_disk/data/H1437/ONT/ClairS_v0_4_1/output.vcf.gz \
  /big8_disk/data/H1437/ONT/H1437.bam 75578

launch_one H2009 \
  /big7_disk/liaoyoyo2001/big7_disk_output/canonical/H2009/paired_full/20260315_H2009_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz \
  /big8_disk/data/H2009/ONT/H2009BL.bam \
  /big8_disk/data/H2009/ONT/ClairS_v0_4_1/output.vcf.gz \
  /big8_disk/data/H2009/ONT/H2009.bam 157405

launch_one HCC1937 \
  /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1937/paired_full/20260315_HCC1937_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz \
  /big8_disk/data/HCC1937/ONT/HCC1937BL.bam \
  /big8_disk/data/HCC1937/ONT/ClairS_v0_4_1/output.vcf.gz \
  /big8_disk/data/HCC1937/ONT/HCC1937.bam 49548

launch_one HCC1954 \
  /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1954/paired_full/20260315_HCC1954_paired_full_full_complete_matrix/longphase_s/germline_phased_merged.vcf.gz \
  /big8_disk/data/HCC1954/ONT/HCC1954BL.bam \
  /big8_disk/data/HCC1954/ONT/ClairS_v0_4_1/output.vcf.gz \
  /big8_disk/data/HCC1954/ONT/HCC1954.bam 20969
```

每次呼叫的實際 output：

```text
$RUN_ROOT/samples/<SAMPLE>/tag_only/<SAMPLE>.lps_tag_only.bam
$RUN_ROOT/samples/<SAMPLE>/tag_only/<SAMPLE>.lps_tag_only.bam.bai
$RUN_ROOT/samples/<SAMPLE>/tag_only/<SAMPLE>.lps_tag_only_purity.out
$RUN_ROOT/samples/<SAMPLE>/longphase_s.log
$RUN_ROOT/samples/<SAMPLE>/resource.time.txt
```

## 8. Resume / skip / collision 規則

LongPhase-S 無 block-level checkpoint，不能對 partial BAM append/resume。只能在 **sample 粒度**重新開始，而且不能覆寫現有檔案。

| 狀態 | 判定條件 | 行為 |
|---|---|---|
| `UNSTARTED` | sample dir 無 output、status 無 START | 可取得 sample lock 後啟動。 |
| `RUNNING` | status START 且 owner/session 明確仍 active；或 output 是 live FIFO | **ATTACH-ONLY**；不得啟動同 sample。 |
| `COMPLETE` | `ACCEPTED` 存在、acceptance 13 項全 PASS、BAM/BAI checksum 重驗 PASS | SKIP；只讀取，不重跑。 |
| `PARTIAL` | BAM 存在但無 ACCEPTED；LongPhase nonzero；BAI/quickcheck 失敗 | STOP。不得把「檔案存在」當完成；由 operator 搬到 archive/incomplete 後，用新 `_retryNN` output prefix 重跑。 |
| `STALE START` | status START，但 owner 無法確認 | 不可自行 kill/重跑；先人工確認外部 session。 |
| `ROLE MISMATCH` | `@PG` 含 truth/output-somatic，或缺 `--disableFilter` | 不可接受為 tag-only；另列 legacy/recalibrated role。 |
| `FIFO` | `test -p output.bam` | 只可能是 streaming sidecar transport；永遠不能滿足 persistent BAM/BAI contract。 |

必要防重：

- 每 sample 使用 `flock -n $sample/.launch.lock`。
- `run_status.tsv` 必須記 `START`、`LONGPHASE_EXIT`、`INDEX_EXIT`、`ACCEPTANCE` 四個獨立事件。
- launch 前若任何預期 output path 已存在，立即 exit 20；不覆寫、不自動刪除。
- resume 時先驗 input snapshot hash 與 original manifest 相同；input 漂移就建立新 run id，不能沿用舊 root。

### 8.1 可直接落地的 fail-closed、frozen provenance、atomic launch gate

本節是**新 wrapper 的規格與命令**，不是對現有 `run_longphase_production_sidecars.sh` 的修改。本輪不得把它貼進或 source 進 live runner/consumer。

#### A. Run root 與 manifest 一次性建立

```bash
set -euo pipefail
umask 002

RUN_ROOT=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_persistent_tag_only_no_truth_v1
test ! -e "$RUN_ROOT"
mkdir "$RUN_ROOT"                         # 無 -p：名稱衝突立即 fail
mkdir "$RUN_ROOT/samples" "$RUN_ROOT/smoke"

{
  printf 'run_id\t20260711_longphase_s_persistent_tag_only_no_truth_v1\n'
  printf 'created_at\t%s\n' "$(date --iso-8601=seconds)"
  printf 'host\t%s\n' "$(hostname)"
  printf 'owner_pid\t%s\n' "$$"
  printf 'contract\tpersistent_tag_only_no_truth\n'
} > "$RUN_ROOT/run_owner.tsv.pending"
mv "$RUN_ROOT/run_owner.tsv.pending" "$RUN_ROOT/run_owner.tsv"
```

Gate 意義：`test ! -e` + non-`-p` `mkdir` 讓同名 root 不能被兩個 session 同時「沿用」。若任一命令失敗，沒有 LongPhase process 被啟動。

#### B. Frozen input/code provenance

先生成 `input_manifest.json.pending`，內容至少對每 sample 固定：logical path、`readlink -f` target、`stat -Lc` size/mtime/inode、BAM BAI SHA-256、raw/caller-pass VCF及 index SHA-256、germline及 index SHA-256、expected PASS count、reference FAI SHA-256。binary 同時鎖 version/commit/hash。

```bash
LPS=/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/longphase-s
test "$(sha256sum "$LPS" | cut -d' ' -f1)" = \
  07cbd0aa0c9f33ed59c5baff45fbe3554ef96d55457635de7348c4501b283f54
test "$(sha256sum /big8_disk/ref/GRCh38_no_alt_analysis_set.fasta.fai | cut -d' ' -f1)" = \
  502a1b8fb73ccd53285c28a0f12df90c818b4fe3de1e862ef47c593ef1a0a4b4

MANIFEST_SOURCE="$RUN_ROOT/input_manifest.rendered.json"
mv "$RUN_ROOT/input_manifest.json.pending" "$MANIFEST_SOURCE"
jq -S . "$MANIFEST_SOURCE" > "$RUN_ROOT/input_manifest.json.canonical.pending"
mv "$RUN_ROOT/input_manifest.json.canonical.pending" "$RUN_ROOT/input_manifest.json"
sha256sum "$RUN_ROOT/input_manifest.json" > "$RUN_ROOT/input_manifest.sha256.pending"
mv "$RUN_ROOT/input_manifest.sha256.pending" "$RUN_ROOT/input_manifest.sha256"
chmod 0444 "$MANIFEST_SOURCE" "$RUN_ROOT/input_manifest.json" "$RUN_ROOT/input_manifest.sha256"

{
  sha256sum "$LPS"
  printf '%s  %s\n' 420e0ff73188126da8bb2122725fbffb2d076ccc longphase_s_git_commit
  /usr/local/bin/samtools --version | sed -n '1,2p'
  /usr/bin/bcftools --version | sed -n '1,2p'
} > "$RUN_ROOT/code_provenance.txt.pending"
mv "$RUN_ROOT/code_provenance.txt.pending" "$RUN_ROOT/code_provenance.txt"
chmod 0444 "$RUN_ROOT/code_provenance.txt"
```

`*.pending → mv` 只在同一 filesystem 內使用；這讓 reader 不會看見半寫 manifest。若 manifest 已存在，wrapper 必須用 `sha256sum -c` 驗證，不能重新 render 後默默替換。

#### C. Capacity 與 7-sample completeness preflight

```bash
AVAILABLE_BYTES=$(df -PB1 --output=avail "$RUN_ROOT" | tail -n 1 | tr -d ' ')
REQUIRED_BYTES=2300000000000
test "$AVAILABLE_BYTES" -ge "$REQUIRED_BYTES" || {
  printf 'CAPACITY_FAIL\tavailable=%s\trequired=%s\n' "$AVAILABLE_BYTES" "$REQUIRED_BYTES" >&2
  exit 30
}

jq -e '.dataset_count == 7 and (.samples | length) == 7' \
  "$RUN_ROOT/input_manifest.json" >/dev/null
sha256sum -c "$RUN_ROOT/input_manifest.sha256"
```

現在的 932 G free 會在 `exit 30` 停止；這正是預期 fail-closed 行為。

#### D. Command argv 先凍結、再驗 forbidden/required flags

不要只保存容易被 quoting 改寫的 command string；每 sample 先寫 JSON argv array。例如：

```json
[
  "/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/longphase-s",
  "somatic_haplotag",
  "-s", "<germline>",
  "-b", "<normal_bam>",
  "--tumor-snv-file", "<caller_pass_snapshot>",
  "--tumor-bam-file", "<tumor_bam>",
  "-r", "/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta",
  "-t", "12",
  "--tagSupplementary",
  "-q", "20",
  "-p", "0.6",
  "--disableFilter",
  "-o", "<sample_root>/tag_only/<sample>.lps_tag_only"
]
```

Fail-closed argv validator：

```bash
ARGV_JSON="$SAMPLE_ROOT/command.argv.json"
jq -e '
  (.[0] | endswith("/longphase-s")) and
  .[1] == "somatic_haplotag" and
  ([.[] | select(. == "--disableFilter")] | length == 1) and
  ([.[] | select(. == "--tagSupplementary")] | length == 1) and
  ([.[] | select(startswith("--truth-vcf") or startswith("--truth-bed") or
                  . == "--output-somatic-vcf" or . == "--somatic-calling-log" or
                  startswith("--region"))] | length == 0)
' "$ARGV_JSON" >/dev/null || exit 31

jq -r '.[] | @sh' "$ARGV_JSON" > "$SAMPLE_ROOT/command.argv.quoted.txt.pending"
mv "$SAMPLE_ROOT/command.argv.quoted.txt.pending" "$SAMPLE_ROOT/command.argv.quoted.txt"
sha256sum "$ARGV_JSON" > "$SAMPLE_ROOT/command.argv.sha256"
chmod 0444 "$ARGV_JSON" "$SAMPLE_ROOT/command.argv.quoted.txt" "$SAMPLE_ROOT/command.argv.sha256"
```

注意：實際執行端應以 argv array 還原執行，不可 `eval` quoted text。`quoted.txt` 只給人審閱。

#### E. Sample lock、狀態原子轉換與 exit capture

```bash
exec {LOCK_FD}>"$SAMPLE_ROOT/.launch.lock"
flock -n "$LOCK_FD" || exit 32

OUT_BAM="$SAMPLE_ROOT/tag_only/$SAMPLE.lps_tag_only.bam"
OUT_BAI="$OUT_BAM.bai"
test ! -e "$OUT_BAM" && test ! -e "$OUT_BAI" || exit 33
sha256sum -c "$RUN_ROOT/input_manifest.sha256" || exit 34
sha256sum -c "$SAMPLE_ROOT/command.argv.sha256" || exit 35

jq -n --arg state RUNNING --arg at "$(date --iso-8601=seconds)" \
  --arg manifest_sha "$(cut -d' ' -f1 "$RUN_ROOT/input_manifest.sha256")" \
  '{state:$state,at:$at,manifest_sha256:$manifest_sha}' \
  > "$SAMPLE_ROOT/state.json.pending"
mv "$SAMPLE_ROOT/state.json.pending" "$SAMPLE_ROOT/state.json"

set +e
set -o pipefail
/usr/bin/time -v -o "$SAMPLE_ROOT/resource.time.txt" \
  "${ARGV[@]}" 2>&1 | tee "$SAMPLE_ROOT/longphase_s.log"
LPS_RC=${PIPESTATUS[0]}
set -e

if [[ "$LPS_RC" -ne 0 ]]; then
  jq -n --arg state FAILED --argjson rc "$LPS_RC" --arg at "$(date --iso-8601=seconds)" \
    '{state:$state,longphase_exit:$rc,at:$at}' > "$SAMPLE_ROOT/state.json.pending"
  mv "$SAMPLE_ROOT/state.json.pending" "$SAMPLE_ROOT/state.json"
  exit "$LPS_RC"
fi
```

`flock` descriptor必須持有到 acceptance 完成；失敗時只寫 FAILED，不刪 partial BAM。wrapper 的 trap 也只能寫狀態，不能 `rm`。

#### F. 驗收後才原子 publish `ACCEPTED`

```bash
/usr/local/bin/samtools quickcheck -v "$OUT_BAM"
/usr/local/bin/samtools index -@ 8 "$OUT_BAM"
test -s "$OUT_BAI"
# 接著跑 §9 的 idxstats/@PG/HP/MM-ML/output-hash checks。
jq -e '.all_pass == true and .check_count >= 13' "$SAMPLE_ROOT/acceptance.json" >/dev/null

{
  printf 'sample\t%s\n' "$SAMPLE"
  printf 'accepted_at\t%s\n' "$(date --iso-8601=seconds)"
  printf 'manifest_sha256\t%s\n' "$(cut -d' ' -f1 "$RUN_ROOT/input_manifest.sha256")"
  printf 'command_sha256\t%s\n' "$(cut -d' ' -f1 "$SAMPLE_ROOT/command.argv.sha256")"
  sha256sum "$OUT_BAM" "$OUT_BAI"
} > "$SAMPLE_ROOT/ACCEPTED.pending"
mv "$SAMPLE_ROOT/ACCEPTED.pending" "$SAMPLE_ROOT/ACCEPTED"
chmod 0444 "$SAMPLE_ROOT/ACCEPTED"

jq -n --arg state ACCEPTED --arg at "$(date --iso-8601=seconds)" \
  '{state:$state,at:$at}' > "$SAMPLE_ROOT/state.json.pending"
mv "$SAMPLE_ROOT/state.json.pending" "$SAMPLE_ROOT/state.json"
```

Fail-closed 關係是：**沒有完整 acceptance → 沒有 `ACCEPTED`；沒有 7 個可重驗 `ACCEPTED` → 沒有 matrix completeness；沒有 completeness → 不得改 layered manifest。**

## 9. Persistent BAM/BAI acceptance contract

每個 sample 同時產 machine-readable `acceptance.json` 與 human-readable `acceptance.tsv`；兩者來自同一次檢查且摘要一致。下列項目全部 PASS 才能建立 `ACCEPTED`：

| # | Check | 強驗證 |
|---:|---|---|
| 1 | input existence/index | 5 inputs 存在；BAM BAI、VCF CSI/FAI 存在。 |
| 2 | input identity | resolved path、`stat -Lc` size/mtime、BAI/VCF/index SHA-256 與 manifest 相同。 |
| 3 | caller universe | snapshot records = expected，record-stream hash = raw `FILTER=PASS && SNV`。 |
| 4 | LongPhase exit | `longphase_exit=0`。 |
| 5 | role log | truth paths empty；filter disabled；write calling VCF disabled；tag region all。 |
| 6 | output type | BAM 是非空 regular file、不是 FIFO/symlink。 |
| 7 | BAM integrity | `samtools quickcheck -v` exit 0 且無訊息。 |
| 8 | coordinate/index | `@HD SO:coordinate`；`samtools index -@ 8` exit 0；BAI nonempty；`idxstats` exit 0。 |
| 9 | alignment conservation | input/output `samtools idxstats` 每 contig mapped/unmapped counts 完全相同。 |
| 10 | current `@PG` | 恰有本次新增 LongPhase-S PG；VN 1.0.0；CL 有 `--disableFilter/-q20/-p0.6/--tagSupplementary`，無 truth/output-somatic/region。 |
| 11 | HP vocabulary | 所有被掃描 alignment 的 HP 只屬 8 個實際字串；`.` 表示 no tag。full acceptance 建議掃完整 BAM。 |
| 12 | MM/ML preservation | 對 input/output exact alignment identity 比對，非 HP/PS auxiliary tags不得被意外移除；至少完整比較 MM/ML。 |
| 13 | output identity | BAM、BAI、purity sidecar、command、log、acceptance 產 SHA-256；寫入 `output.sha256`。 |

### 9.1 可直接重驗的核心 acceptance commands

下列變數取自 frozen manifest，不允許由「latest directory」推測：

```bash
INPUT_BAM=<frozen_tumor_bam>
OUT_BAM=<sample_root>/tag_only/<sample>.lps_tag_only.bam
OUT_BAI="$OUT_BAM.bai"
PREFIX="${OUT_BAM%.bam}"
LOG=<sample_root>/longphase_s.log
```

檔案型別、完整性、BAI 與 alignment conservation：

```bash
test -f "$OUT_BAM" && test ! -L "$OUT_BAM" && test ! -p "$OUT_BAM" && test -s "$OUT_BAM"
/usr/local/bin/samtools quickcheck -v "$OUT_BAM"
test -s "$OUT_BAI"
/usr/local/bin/samtools idxstats "$OUT_BAM" >/dev/null
diff -u \
  <(/usr/local/bin/samtools idxstats "$INPUT_BAM") \
  <(/usr/local/bin/samtools idxstats "$OUT_BAM")
```

Log role 與禁止副產物：

```bash
grep -Eq '^truth VCF file[[:space:]]*:[[:space:]]*$' "$LOG"
grep -Eq '^truth BED file[[:space:]]*:[[:space:]]*$' "$LOG"
grep -Eq '^tag region[[:space:]]*:[[:space:]]*all$' "$LOG"
grep -Eq '^variant filtering[[:space:]]*:[[:space:]]*disabled$' "$LOG"
grep -Eq '^write calling VCF[[:space:]]*:[[:space:]]*disabled$' "$LOG"
! grep -q 'removing tumor & truth somatic variants outside bed regions' "$LOG"
test ! -e "${PREFIX}_sc.vcf"
test ! -e "${PREFIX}_somatic_haplotag.metrics"
```

`@PG`（input preflight 已驗證 7 個 raw tumor BAM 都沒有既有 LongPhase-S PG，因此 output 應恰有 1 行）：

```bash
PG_LINES=$(/usr/local/bin/samtools view -H "$OUT_BAM" | awk '$1=="@PG" && $0~/PN:longphase-s/')
test "$(printf '%s\n' "$PG_LINES" | grep -c .)" -eq 1
printf '%s\n' "$PG_LINES" | grep -q -- '--disableFilter'
printf '%s\n' "$PG_LINES" | grep -q -- '--tagSupplementary'
printf '%s\n' "$PG_LINES" | grep -q -- '-q 20'
printf '%s\n' "$PG_LINES" | grep -q -- '-p 0.6'
! printf '%s\n' "$PG_LINES" | grep -Eq -- '--truth-vcf|--truth-bed|--output-somatic-vcf|--region'
```

完整 HP vocabulary scan：

```bash
/usr/local/bin/samtools view "$OUT_BAM" | awk -F '\t' '
BEGIN {
  ok["1"]; ok["2"]; ok["3"]; ok["4"];
  ok["1-1"]; ok["2-1"]; ok["1-2"]; ok["2-2"]
}
{
  for (i=12; i<=NF; i++) {
    if ($i ~ /^HP:Z:/) {
      hp=substr($i,6)
      if (!(hp in ok)) { print "UNKNOWN_HP\t" hp "\t" $1 > "/dev/stderr"; bad=1 }
    }
  }
}
END { exit bad }
'
```

MM/ML preservation stream digest。這會完整掃 input/output，成本必須列入 acceptance I/O window；digest 同時含 alignment identity，避免只比 tag multiset：

```bash
aux_digest() {
  /usr/local/bin/samtools view "$1" | awk -F '\t' 'BEGIN{OFS="\t"}
  {
    mm="."; ml="."
    for (i=12; i<=NF; i++) {
      if ($i ~ /^MM:Z:/) mm=$i
      else if ($i ~ /^ML:B:C,/) ml=$i
    }
    print $1,$2,$3,$4,$6,mm,ml
  }' | sha256sum | cut -d" " -f1
}
test "$(aux_digest "$INPUT_BAM")" = "$(aux_digest "$OUT_BAM")"
```

若 input 使用新版 `Mm/Ml` 大小寫 alias，validator 必須先從實際 header/records擴充匹配；不可因 regex 沒找到就把兩個空 digest 誤判為 preservation PASS。前置條件是 input 的 MM/ML count >0；HCC1395 5 kHz 等 methylation BAM 應強制此 gate。

Matrix 完成條件：

```text
dataset_count=7
n_accepted=7
n_failed=0
n_partial=0
all_truth_flags_absent=true
all_persistent_bam=true
all_bai=true
all_idxstats_conserved=true
```

只有達到上述條件，`completeness.tsv` 才能標 `ready_for_layered_manifest=true`。

## 10. 資源估算與 capacity hard gate

歷史 run 帶 truth/recalibration，故只作 production tag-only 的近似估算；其 BAM read payload 與本次相同量級。

| Dataset | 歷史 wall | Max RSS | 歷史 tagged BAM |
|---|---:|---:|---:|
| HCC1395（24 threads） | 1:38:56 | 12.0 GiB | 278.3 GB |
| HCC1395_DORADO | 2:14:45 | 10.2 GiB | 238.4 GB |
| COLO829 | 1:44:48 | 6.8 GiB | 95.4 GB |
| H1437 | 1:42:40 | 8.4 GiB | 230.1 GB |
| H2009 | 4:45:35 | 12.7 GiB | 310.8 GB |
| HCC1937 | 3:45:06 | 8.5 GiB | 447.2 GB |
| HCC1954 | 3:04:42 | 7.2 GiB | 240.8 GB |
| **合計** | **18:56:32 serial** | — | **1.841 TB** |

- 12 threads、最多 2 samples parallel：預估 10–14 h；CPU utilization 歷史僅 209%–550%，I/O contention 是主要不確定性。
- 2 parallel RAM 建議 ≥32 GiB available；若與其他研究 job 共存，建議 ≥48 GiB headroom。
- persistent output 預估 1.84 TB，±5% 約 1.75–1.93 TB；加 20% headroom 的 hard gate 取 **2.3 TB free**。
- 2026-07-11 實測：`/dev/sdc1 42T, 39T used, 932G available, 98%`。故缺口約 1.2–1.3 TB；且磁碟已 98%，不應用「先跑到滿」試探。
- 現行 `MIN_FREE_GB_DEFAULT=800` 只適合單一 run 的粗 guard，不足以保護 7-output matrix；本次要用 projected-output gate。

**容量未解除時的合法替代只有兩種**：保留既有 streaming sidecar（無 BAM/BAI），或另立 CRAM contract。兩者都不能假稱完成本次 persistent BAM deliverable。

## 11. 已存在 collaborator run：只讀狀態與角色差異

既有 root：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_production_sidecars_v1
```

2026-07-11 本次讀取到：

```text
HCC1395_DORADO  production_tagging  START  expected_records=112387;truth_flags=absent
HCC1395         production_tagging  START  expected_records=113997;truth_flags=absent
```

該 run 的正確定位：

- `--truth-vcf/--truth-bed` 確實 absent，這部分是 clean genome-wide。
- command **有** `--output-somatic-vcf`，且 log 明示 `variant filtering: enabled`、`write calling VCF: enabled`；這是 `tagging + recalibration` role，不是本文件定義的 caller-PASS-preserving tag-only。
- `<sample>_production.bam` 是 named pipe；consumer 擷取 read-tag TSV，BAM payload 不保留，完成後 FIFO 改名。它不能產 BAI，也不能通過 §9。
- launcher line 18 若 run root 已存在便整體 abort，因此沒有原地 resume；partial failure需新的 run id 或人工設計 resume。
- `input_files.tsv` 用 `stat -c`，對 symlink 記到 79/80-byte link 本身而非 91–250 GB target；正式 provenance 要改用 `stat -Lc` + `readlink -f`。本審計只紀錄，不修改該 script。
- sidecar capture沒有保存 BAM header `@PG`；command file/log可證命令，但若要完整外部 audit，應另存 streamed header。

**Ownership rule**：不啟動、重啟、kill、搬移、覆寫此 root；本文件的 planned persistent root 使用不同名稱。其未來 PASS 只能算 sidecar/recalibration evidence，不能替代 persistent tag-only completeness。

## 12. 風險登錄

| Priority | 風險 | 已知證據 | Control |
|---|---|---|---|
| P0 | 7 BAM 容量不足 | 需約 2.3 TB free，現有 932 G | HOLD full launch；先擴容/配置其他合法 volume。 |
| P0 | truth BED 改變 tag universe | source 在 `tagRead` 前移除 BED 外 tumor variants | production CL 禁止 truth；`@PG`/log雙驗。 |
| P0 | 預設 filter 改變 caller PASS | source `variantCalling` 永遠跑；filter 預設 true | tag-only 必帶 `--disableFilter`。 |
| P1 | 同名 VCF role 覆寫 | current step 先寫 caller PASS、後寫 recalibrated PASS | role-specific filename + immutable snapshot + stream hash。 |
| P1 | partial BAM 被當成可 resume | LongPhase 無 checkpoint；pipeline 只看 path/新 run | ACCEPTED sentinel、quickcheck/BAI、partial stop/quarantine。 |
| P1 | external run collision | 已有 2 sample START/FIFO | attach-only + distinct run root + `flock`。 |
| P1 | symlink provenance錯記 | sidecar `stat -c` 得 79/80 bytes | `readlink -f` + `stat -Lc`。 |
| P1 | 大 BAM 無 content hash | 本輪為避免額外 1.9 TB I/O 未算 full SHA | 明示 limitation；外部 handoff 前另排 checksum。 |
| P2 | `purity`術語誤讀 | KB 已正名 tumor DNA fraction | 報告只用 DNA fraction；不升級 cellular purity。 |
| P2 | output role 被 promotion 過早 | synthesis round 尚無 7/7 completeness | acceptance 13 項 + matrix completeness 後才 promotion。 |

## 13. Step → Verify launch checklist

1. 擴容或指定合法 persistent volume → 驗證：同一 filesystem free ≥2.3 TB；輸出完整路徑記錄於 manifest。
2. 不干擾既有 sidecar collaborator run → 驗證：planned root 名稱不同，HCC1395/HCC1395_DORADO 無 duplicate process/lock。
3. 建立 7 組 role-specific caller-PASS snapshot → 驗證：counts 等於 §3.2，record hashes 7/7 match raw PASS。
4. 跑 HCC1395 bounded smoke → 驗證：§6 八項全 PASS，且標 PARTIAL_SMOKE_ONLY。
5. 人工檢閱 smoke `@PG` 與 log → 驗證：無 truth/output-somatic，filter disabled，HP/MM/ML contract PASS。
6. 以 sample-level lock 執行 full command → 驗證：每個 LongPhase exit 0、實際 output/log/resource path 可見。
7. 建 BAI 與完整性驗收 → 驗證：§9 13 項 per sample 全 PASS；任何一項 FAIL 即無 ACCEPTED。
8. 建 matrix completeness → 驗證：7/7 accepted、0 partial、0 role mismatch。
9. 僅在 7/7 後更新 layered manifest/canonical pointer → 驗證：manifest 指向新 `*.lps_tag_only.bam`，不再指向 truth-conditioned historical BAM。

## 14. 本次實際執行的唯讀驗證與輸出片段

### 14.1 LongPhase-S help/source/binary

輸入：已知 binary/source exact paths。命令：

```bash
/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/longphase-s somatic_haplotag --help
sha256sum /big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/longphase-s
git -c safe.directory=/big8_disk/liaoyoyo2001/knowledge/codebase/longphase-s \
  -C /big8_disk/liaoyoyo2001/knowledge/codebase/longphase-s rev-parse HEAD
```

實際片段：

```text
--disableFilter        disable somatic variant filtering and accept all tumor VCF variants
--output-somatic-vcf   default: false
--truth-vcf            benchmark argument
--truth-bed            only using variants within these regions for tagging and evaluation
version 1.0.0; commit 420e0ff...; sha256 07cbd0aa...
```

### 14.2 7-input checks

命令類型：`samtools quickcheck -v`、BAI/TBI/CSI existence、`bcftools query` GT/PS count、BAM/VCF/FAI 24-primary-contig dictionary diff、caller PASS record hash comparison。

實際片段：

```text
7/7 tumor quickcheck PASS
7/7 normal quickcheck PASS
HCC1395 germline 1698116 total / 1687835 phased+PS
...
HCC1954 germline 1850528 total / 1841976 phased+PS
7/7 tumor=PASS normal=PASS rawvcf=PASS germline=PASS (chr1-22,X,Y dictionary)
7/7 raw PASS record hash == canonical caller-PASS record hash
```

### 14.3 歷史資源與磁碟

輸入：7 個 exact canonical `longphase_s.log`、7 個 exact BAM。命令：bounded `awk/stat/df`。

實際片段：

```text
historical wall range: 1:38:56–4:45:35
historical max RSS: 7,088,124–13,316,560 kB
historical tagged BAM sum: 1,840,983,466,353 bytes
/big7_disk available: 932G; use: 98%
```

本次文件輸出：

```text
InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/04_clean_tagging_launch_audit.md
```

## 15. 權威來源

1. `InterSubMod/scripts/pipeline/steps/01_longphase_s.sh`（實際 benchmark/tagging wrapper）。
2. `InterSubMod/scripts/pipeline/run_benchmark.sh`（orchestration/skip 行為）。
3. `InterSubMod/scripts/pipeline/config.sh`（7 sample paths與 reference）。
4. `/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-s.md`（active runtime-fact，last verified 2026-06-11）。
5. `/big8_disk/liaoyoyo2001/knowledge/06_workflows/phasing-workflow.md`（active runtime-fact，last verified 2026-06-11）。
6. `/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/src/somatic_haplotag/{SomaticHaplotag.cpp,SomaticHaplotagProcess.cpp,SomaticVarCaller.cpp}`。
7. `/big8_disk/liaoyoyo2001/Knowledge/codebase/longphase-s/src/haplotag/HaplotagParsingBam.cpp`。
8. 7 個 §3 所列 canonical run logs/BAM/VCF 與 `layered_v2_input_manifest.json`。
9. 既有 collaborator root `20260711_longphase_s_production_sidecars_v1` 的 manifest、command、log、status與 bounded launcher source。

---

**最終 handoff**：命令與驗收契約已可審核；目前不能誠實宣稱「可立即 full launch」，因 persistent BAM capacity gate 尚未通過。既有 sidecar run要繼續由其 owner 管理；本文件不改變其狀態。
