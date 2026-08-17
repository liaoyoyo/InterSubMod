<!--
建立時間: 2026-08-13
目標: 記錄一頁式多 BAM dashboard、bounded ingestion contract、甲基化觀察圖與瀏覽器終驗
處理範圍: 7 technical datasets（6 biological IDs + 1 technical replicate）；HCC1395-only downstream 保持隔離
關聯檔案:
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_input_manifest.schema.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_input_manifest.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_artifact.json
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_multi_bam_analysis_overview.standalone.html
  - InterSubMod/research/20260813_hcc1395_drilldown_validation/results/multi_bam_dashboard_browser_qa.json
-->

用 SCQA + verification chain：**多 BAM 一頁式 dashboard 與 v1.1 ingestion contract 已通過 bounded data-product QA；全 7 datasets 可選取比較，甲基化圖已補 exact N/D 與 invalid-axis gate，但 science 仍 `PARTIAL`、production 尚未完成（影響：高，信心：高）。**

# 多 BAM 一頁式觀察面板驗證紀錄

> [!WARNING]
> **Task Type B — Comprehensive validation。** 本輪納入 7/7 technical datasets＝6 biological IDs＋1 technical replicate；不是 subset。工程與 snapshot data-product 可驗收，但沒有 truth、KDE provenance 或等價全樣本 ISM／lineage，因此不可提升為 caller accuracy、phasing accuracy、甲基化因果或 cellular lineage 結論。

## 1. 最終 verdict

| 層級 | Verdict | 可安全宣稱 | 不可宣稱 |
|---|---|---|---|
| Ingestion contract | **PASS_BOUNDED_WITH_METADATA_DRIFT** | 來源 v3/schema、同樣本 producer evidence chain、固定三段 1 MiB、BAI/FAI、header、quickcheck 可稽核 | full BAM content identity、RG biological identity |
| Dashboard engineering | **PASS** | standalone、離線、selector isolation、responsive、keyboard、ARIA、無外連 | production live ingestion、URL state、async retry/race |
| Data product | **PASS_BOUNDED** | 7-dataset topology、bounded BAM readiness、existing producer tag receipt、HCC-only diagnostic | 完整 BAM QC、全樣本 ISM/lineage、truth-aware cohort analysis |
| Science | **PARTIAL / BLOCKED FOR UPGRADE** | 描述性 observation 與 denominator diagnostics | F1、phasing accuracy、KDE-corrected effect、platform effect、LCA causality |

服務目標：**G3** read-level epigenetic 可稽核呈現、**G4** 多樣本一致性與 reproducibility、**G5** 可被外部重跑的資料產品契約。

## 2. 最終交付與內容雜湊

| 交付 | InterSubMod 路徑 | SHA256 | 大小 |
|---|---|---|---:|
| Output JSON Schema v1.1 | `InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_input_manifest.schema.json` | `786126a683cbdad09ff0f890b460f0eedf4d9ab50a77ee99c0fc6ebc02dc86f0` | 18,455 B |
| Canonical BAM manifest | `InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_input_manifest.json` | `2e290c06be9e98be480ea913b62b3de30ad827eb9e6894a5a68b4e5dfa53f190` | 44,329 B |
| Manifest validation receipt | `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/multi_bam_input_manifest_validation.json` | `0a18c63c35920d9c089ccf651850b1a4f3242caa681bcd8e4b87c447415d130b` | 20,287 B |
| Dashboard artifact | `InterSubMod/research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_artifact.json` | `e8e9c2bf91d56de269a3d6868ac1e001e27c35dd93fd43a920e1d0d5f8a7ea0b` | 323,107 B |
| Standalone HTML | `InterSubMod/research/20260813_hcc1395_drilldown_validation/20260813_multi_bam_analysis_overview.standalone.html` | `224122cfc41ba24acc886f7e04b7f6d3553887be1118125b1f541a9646d0476c` | 1,134,761 B |
| Browser QA receipt | `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/multi_bam_dashboard_browser_qa.json` | `0b72cb659f918913e2f25f6a6cb997354b031387cab9ca8a40c7fdfcf319b5ef` | 18,502 B |

上述 contracts 與 6 個直接 dashboard source files 已由 `.gitignore` 精確 negation 放行，canonical source schema 亦已放行。**目前工作樹仍未 stage／commit，因此只可宣稱「可被版本控制納入」；不可宣稱另一個 clean checkout 已具備這些新檔。** 外部 frozen source manifest、BAM/BAI/FAI 與 producer receipts仍是環境依賴。

## 3. 輸入與權威邊界

### 3.1 Ingestion inputs

- Frozen v3 source manifest：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json`
- Canonical source schema：`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/schemas/layered_input_manifest_v3.schema.json`
- Cohort topology：`InterSubMod/research/20260813_hcc1395_drilldown_validation/results/cohort_topology_metrics.csv`
- 七份 producer receipts：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/<dataset>/producer_capture_receipt_v2.json`

### 3.2 Dashboard analytical inputs

- `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/bundle_overview.csv`
- `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/methylation_coverage_by_k.csv`
- `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/methylation_axis_metrics.csv`
- `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/metrics_audit.json`
- `InterSubMod/research/20260813_hcc1395_drilldown_validation/results/legacy_browser_visual_audit.generated.json`

資料分成三層：7-dataset topology；7-dataset bounded BAM/tag readiness；HCC1395-only v1/v3 downstream。第三層對其他 6 datasets 一律 `ABSENT_NO_EQUIVALENT_BUNDLE`／null，不借用 HCC1395。

## 4. Ingestion contract 改進與驗證

### Step 1 → Verify：schema 與 semantic validation

Builder 現在會在任何 BAM I/O 前，實際執行 bounded Draft 2020-12 validator：支援兩份 schema 使用的 `$ref`、`allOf`、`oneOf`、`prefixItems`、`items:false`、`const`、`required`、`additionalProperties`、type/range/pattern/date-time 等；遇到未支援 keyword 直接拒絕，不會靜默略過。comprehensive intake 另固定 canonical source schema `$id`＋SHA256，不能用 permissive `--source-schema` 繞過；之後再跑跨欄、跨檔 semantic validation。

負向測試確認下列都會拒絕：claim ceiling、task type、authority／assurance、denominator 被改寫；required 欄缺失；unexpected property；2 MiB chunk；任意 chunk label；跨樣本 receipt；output/receipt/input path collision。

### Step 2 → Verify：同樣本 producer evidence chain

每一列用 `sample` 綁定，不用 `biological_id`（避免 HCC1395 與 HCC1395_DORADO 交叉）：

1. source sample＝subject-binding sample＝producer receipt sample；receipt paths 7/7 唯一。
2. receipt file SHA＝source receipt SHA＝subject binding SHA。
3. sidecar、index、native validation 與 LongPhase-S all/pass VCF 必須完整等於 `capture_outputs`。
4. tumor storage identity、caller/raw-all/baseline VCF 必須等於 `producer_inputs`。
5. command argv、producer inputs、effective options 以 canonical JSON 重算 SHA。
6. mapped＝unique＋duplicate，conflict＝0，且 subject/receipt/validation counts 逐欄一致。
7. Dashboard 會再次核對 receipt status、完整 `verification_summary`、claim ceiling、source/output schema、builder、source manifest 與 topology SHA；不接受只保留 output SHA 的偽造 receipt。

### Step 3 → Verify：bounded BAM policy

- 固定 `chunk_size_bytes=1,048,576`。
- 固定三段且依序為 first／middle／last；offset 由 frozen size 重算；每段必須恰為 1 MiB。
- BAM header、`samtools quickcheck`、BAI full SHA 與 reference FAI/dictionary 均核對。
- quickcheck 前後重讀 bounded identity，若檔案在驗證中改變即拒絕。
- subprocess timeout＝120 秒；沒有 flagstat、mosdepth、read-length、MM/ML 或 full BAM scan。
- output 與 receipt 原子寫入，且不能互撞或覆寫任一 input；dashboard artifact 也改為原子寫入，並保護 implicit source/output schema 與兩個 builder scripts。
- strict identity 聚合由每個 manifest row 動態重算；all-MATCH 與 mixed MATCH/drift 都不再被目前 0/7、7/7 observation 寫死。

### 完整命令與實際輸出

```bash
python3 research/20260813_hcc1395_drilldown_validation/scripts/build_multi_bam_input_manifest.py \
  --source-manifest /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json \
  --source-schema docs/methodology/_assets/20260627_subclone_4axis_teaching/schemas/layered_input_manifest_v3.schema.json \
  --topology-csv research/20260813_hcc1395_drilldown_validation/results/cohort_topology_metrics.csv \
  --schema research/20260813_hcc1395_drilldown_validation/multi_bam_input_manifest.schema.json \
  --output research/20260813_hcc1395_drilldown_validation/multi_bam_input_manifest.json \
  --receipt research/20260813_hcc1395_drilldown_validation/results/multi_bam_input_manifest_validation.json
```

實際：exit 0；7 datasets／6 biological；tumor/normal quickcheck 7/7；tumor/normal/reference bounded payload 7/7；dictionary 7/7；sidecar receipt/conservation 7/7；strict all-role identity 0/7；mount-device-only drift 7/7；full BAM SHA 0/7；無 `@RG` 7/7；HCC1395 cross-directory pairing 1/7。

```bash
python3 -m unittest -v \
  research/20260813_hcc1395_drilldown_validation/scripts/test_multi_bam_input_manifest.py
```

實際：**32/32 tests PASS，exit 0**。新增 regression 覆蓋 canonical/permissive source schema、receipt summary coherence、all-MATCH、mixed MATCH/drift dashboard、動態文案與 implicit-input collision。

## 5. 已接入的 alignment-tag 指標

所有分母都是 producer capture 的 mapped **alignment records**，包含 primary／secondary／supplementary；不是 raw reads、unique reads 或 primary-read rate。

| Dataset | HP / all | HP+PS / HP-tagged | Duplicate identities / all | 解讀 |
|---|---:|---:|---:|---|
| COLO829 | 53.33% | 97.63% | 0.000036% | receipt denominator diagnostic |
| H1437 | 52.32% | 95.73% | 0.000083% | receipt denominator diagnostic |
| H2009 | 60.82% | 93.43% | 0.000041% | receipt denominator diagnostic |
| HCC1395 | 43.32% | 94.29% | **16.14%** | 高 duplicate composition；不可直接作樣本品質排名 |
| HCC1395_DORADO | 49.16% | 93.88% | 0.000027% | technical replicate；不可宣稱 platform effect |
| HCC1937 | 61.55% | 97.54% | **25.94%** | 高 duplicate composition；不可直接作 phasing accuracy |
| HCC1954 | 60.23% | 99.42% | 0.000045% | receipt denominator diagnostic |

`duplicate_identity_conflicts=0` 只表示重複 alignment identity rows 沒有互相衝突，不代表 duplicate rows 可忽略。因此 dashboard 把 duplicate/all 畫成第三條 denominator guard。

## 6. Dashboard 圖表與版型改進

### 主層

- Authority bar 常駐 `PARTIAL`；兩個 blocker 名稱常駐、原因按需展開。
- All 預設顯示 7 datasets、6 biological、1 technical；macro 排除 HCC1395_DORADO。
- Topology 主圖保留每 dataset exact numerator／denominator；不 pooled loci。
- 另建 no-script horizontal opportunity panel，7 個 dataset 名稱、bar 與 exact count 同列，`HCC1395_DORADO · technical` 不只靠顏色辨識。
- 四層 evidence anchor：aggregate、HCC1395 canonical、H2009 observed extreme、HCC1395_DORADO technical control。

### 診斷層

- Active-k 圖旁新增 8-row exact N/D rail。
- 甲基化有效 axis 圖排除 circular cluster；cluster exact value 保留在 validity rail，直接標 `INVALID · DOUBLE-DIPPING`。
- Axis rail 壓縮成 `N / D · %`，v1 raw/BH、v3 raw/BH 四欄可直接比較。
- HCC1395 v1/v3 仍為 `BLOCKED` observation：舊版方法只參考 sticky selection、progressive disclosure、證據與案例共置；不套用 legacy A/B class、V threshold、472 candidates 或 14 displayed cases。

### 細節層

- 10 張精確表格；6 張大型 audit tables 預設收合，button 有 `aria-controls` 且指向 labeled region。
- 公式、缺失欄位、BAM identity caveat 用四個原生 `<details>` 按需展開。
- 桌機表格保留欄位；窄螢幕只在 table container 局部橫捲，不造成 page overflow。

![All scope：PARTIAL authority、selector、priority cards 與主圖](figures/19_multi_bam_dashboard_all_desktop.png)

![390px：blocker details 收合、短摘要與優先卡](figures/21_multi_bam_dashboard_mobile.png)

![7 dataset horizontal opportunity panel](figures/28_multi_bam_dashboard_opportunity_fixed_labels.png)

![甲基化 axis exact N/D 與 validity rail](figures/29_multi_bam_dashboard_denominator_rails.png)

## 7. Artifact、standalone 與 browser QA

### 產生 artifact

```bash
python3 research/20260813_hcc1395_drilldown_validation/scripts/build_multi_bam_dashboard_artifact.py \
  --results-dir research/20260813_hcc1395_drilldown_validation/results \
  --bam-manifest research/20260813_hcc1395_drilldown_validation/multi_bam_input_manifest.json \
  --bam-manifest-receipt research/20260813_hcc1395_drilldown_validation/results/multi_bam_input_manifest_validation.json \
  --output research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_artifact.json
```

實際：exit 0；18 snapshot datasets、19 sources、27 manifest blocks、10 cards、6 charts、10 tables、2 required access issues。Builder 會重新跑 manifest v1.1 schema＋semantic contract，並核對 builder/schema/output SHA，不只相信 receipt 名稱。

### 封裝 standalone

```bash
node research/20260813_hcc1395_drilldown_validation/scripts/deliver_multi_bam_dashboard.mjs \
  --input research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_artifact.json \
  --output research/20260813_hcc1395_drilldown_validation/20260813_multi_bam_analysis_overview.standalone.html \
  --screenshot research/20260813_hcc1395_drilldown_validation/figures/multi_bam_delivery_failure.png \
  --timeout-ms 120000 --ready-timeout-ms 45000 --action-timeout-ms 15000
```

實際：exit 0；`validation/package/verification=passed`；35 rendered blocks、6 charts、10 tables、3 no-script HTML blocks；1440/390 與 source dialog PASS。

### Playwright 終驗

```bash
python3 research/20260813_hcc1395_drilldown_validation/scripts/run_multi_bam_dashboard_qa.py \
  --html research/20260813_hcc1395_drilldown_validation/20260813_multi_bam_analysis_overview.standalone.html \
  --artifact research/20260813_hcc1395_drilldown_validation/multi_bam_dashboard_artifact.json \
  --receipt research/20260813_hcc1395_drilldown_validation/results/multi_bam_dashboard_browser_qa.json \
  --figures-dir research/20260813_hcc1395_drilldown_validation/figures
```

實際：**40/40 assertions PASS**；All＋7 selector states；1440、1024、512、390、320；console error 0、page error 0、external HTTP(S) request 0；11 screenshots。各窄螢幕 `document.scrollWidth == clientWidth`，寬表格維持 local `overflow-x:auto`。

## 8. 獨立審查與 Claude Code

三個獨立 reviewer 的終版重查涵蓋 ingestion P1 regression、來源鏈/數值重算與 UI/metric semantics：

- 來源/數值終版：**PASS，682/682** 輕量 checks＋artifact contract PASS；8/8 input pins、42-row tag chart、14-row tag rail、14-row identity rail、8-row BAM rail皆與 manifest/receipt 一致，未重讀大型 BAM payload。
- UI/metrics 終版：**ACCEPT**；最新 artifact/HTML/QA hash chain一致，官方 40/40 PASS；selector、cluster invalid gate、active-k rail、duplicate guard、ARIA 與 1024/512/390/320 responsive均無 P0/P1。
- Ingestion 對抗終版：**ACCEPT**；32/32 PASS；all-MATCH=`strict 7/drift 0`、mixed All=`1/6`、HCC1395=`1/0`、H1437=`0/1`，舊的寫死 drift 文案 0 hits；permissive schema、receipt tampering 與 implicit collision 皆 fail-closed。唯一 handoff 注意為新檔仍待 stage/commit。

Claude Code 2.1.229 已以終版檔案、high effort、唯讀 allowed tools 再嘗試一次。實際輸出：缺少 `socat` 的 sandbox warning，接著在任何 Read/Bash 前回傳 `You've hit your session limit · resets 2:20pm (Asia/Taipei)`，exit 1。因此狀態是 **NOT EXECUTED**，不是 ACCEPT／REJECT；`InterSubMod/research/20260813_hcc1395_drilldown_validation/claude_code_round5_review.md` 保存收據。

## 9. 尚未完成、下一輪優先順序

1. **P0 science**：truth VCF＋HighConf BED＋som.py/hap.py receipt；沒有它不能談 F1/accuracy。
2. **P0 methylation**：MM/ML、CpG call/site/read grain、threshold、KDE input/output hashes與方法參數。
3. **P0 multi-sample**：生成 7 個 sample-specific、同 contract 的 ISM/lineage bundles；目前只有 HCC1395。
4. **P1 BAM QC**：明確授權後收集 primary mapping、depth/IQR/breadth、read N50、phase blocks；不可由 BAM size 或 alignment receipt 代理。
5. **P1 identity**：補可靠 RG/SM 或 canonical dataset manifest；目前 14 BAM 無 `@RG`，HCC1395 BAM/caller 跨目錄須人工確認。
6. **P1 production**：biological→dataset cascade、URL state、async/retry/race、refresh ownership、locus-detail view 與 release smoke。
7. **P2 storage attestation**：若 release 需要 full BAM SHA，另開有 I/O 預算的受控工作；本輪明確沒有掃描 14 個完整 BAM。

---

**PARTIAL footer — claim ceiling**：本頁可標 `bounded snapshot data-product QA PASS / science PARTIAL / production NOT COMPLETE`。不得把 bounded payload、tag receipt、視覺 QA 或 cluster 100% 數值改寫成生物真值、phasing accuracy、KDE-corrected effect、platform effect或 lineage causality。

Evidence ledger：`InterSubMod/research/autoresearch/evidence_ledger.jsonl`；本次更新以 superseding entry 保存，不覆寫舊快照。
