<!--
建立時間: 2026-07-22
目標: 保存 HCC1395 exact-PS technical report 的來源、圖表契約與重現命令
處理範圍: PARTIAL / exploratory pilot / HCC1395 chr1-22 only
關聯檔案: InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/report/artifact.json
-->

# HCC1395 exact-PS report source notes

> **PARTIAL / EXPLORATORY PILOT**：此檔為 report supporting artifact，不是第二份 reader-facing report。

## Reporting job

- 問題：每個 exact PS 是否應分開保留 read constraint，以及 k≤12 分割在 HCC1395 的量化影響。
- 主要讀者：technical / 教授與研究方法 reviewer。
- 決策：是否接受 exact-PS fail-closed policy，並允許進入 HCC1395 downstream adapter 階段。
- Claim ceiling：extraction、exact-PS segmentation、C++ parity；不包含 T/Topo/VAF/clone/subclone。

## Input paths

- Run root：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2`
- Run receipt：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2/run_receipt.json`
- SEQC2 truth VCF：`/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`
- SEQC2 HC BED：`/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed`
- Input manifest：`/big7_disk/liaoyoyo2001/InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/input_contract/big7_hcc1395_pilot_manifest.json`

## Recompute command

```bash
PYTHONDONTWRITEBYTECODE=1 python research/20260722_exact_ps_k12_hcc1395_pilot/scripts/build_hcc1395_validation_report.py \
  --run-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2 \
  --truth-vcf /big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz \
  --hc-bed /big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed \
  --input-manifest /big7_disk/liaoyoyo2001/InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/input_contract/big7_hcc1395_pilot_manifest.json \
  --output-dir /big7_disk/liaoyoyo2001/InterSubMod/research/20260722_exact_ps_k12_hcc1395_pilot/report
```

Expected output excerpt:

```text
aggregate_exact_match=true
S=79687
primary_unique=36384
units=39544
cross_ps=0
python_cpp_mismatch=0
```

## Chart map

| Section | Question | Family/type | Fields | Supported claim | Palette |
|---|---|---|---|---|---|
| Support disposition | k≤12 loss concentrates where? | Composition / 100% stacked bar | scope, molecule_weight, disposition | Overall and k>12 denominators differ materially | categorical, direct legend |
| Chromosome coverage | Which chromosomes lose primary eligibility? | Comparison / bar | chrom, primary_rate | chr6/chr16 are upstream-scope outliers | single-root blue, no redundant legend |

## Portable delivery exception and QA

Canonical `deliver_portable_artifact.mjs` 先完成 artifact validation、native reader packaging 與 static-chart extraction，但共用 reader 的 `.analytics-top-bar { width: 100vw; margin: calc(50% - 50vw) }` 在 Linux non-overlay scrollbar 下固定造成 8 px document overflow。這是 renderer chrome 問題，不是 report block、table 或 chart data overflow。

最終 delivery 保留同一 canonical artifact payload、native reader、chart runtime、light/dark static SVG extraction 與 canonical verifier；只加入具名 CSS correction：

```css
.analytics-top-bar {
  width: 100% !important;
  max-width: 100% !important;
  margin-right: 0 !important;
  margin-left: 0 !important;
}
```

執行命令：

```bash
node research/20260722_exact_ps_k12_hcc1395_pilot/scripts/deliver_report_with_topbar_fix.mjs \
  --artifact research/20260722_exact_ps_k12_hcc1395_pilot/report/artifact.json \
  --output research/20260722_exact_ps_k12_hcc1395_pilot/report/report.html \
  --receipt research/20260722_exact_ps_k12_hcc1395_pilot/report/report_delivery_receipt.json
```

實際 QA：`ok=true`；static charts=2；blocks=27、charts=2、tables=7、metrics=6；viewports=1440/390；source dialog=`passed`、keyboard interaction=`keyboard_menu_semantic_click`。

## Evidence and caveats

- Raw gzip TSV is recomputed before comparison with `run_receipt.json`.
- `molecule_weight` is molecule–unit incidence, not globally unique reads.
- SEQC2 HC is a benchmark scope. Loci outside it cannot be labeled false positive solely from non-overlap.
- chr6p/chr16q are documented complete-loss/benchmark-excluded regions. The current pilot manifest's unlisted-CN=`neutral` policy is unsafe for downstream CN/clone interpretation, although segmentation does not consume CN.
- The direct Big7 BAM supplies read evidence at a manifest-bound candidate catalog; it did not de novo call the 79,687 loci.

## Audience-spec mapping

1. Title → report title block.
2. Technical summary → `technical_summary`.
3. Key findings with visual evidence → exact-PS risk, support disposition, chromosome outlier sections.
4. Scope/data/metric definitions → definitions table and visible PARTIAL notice.
5. Methodology → exact-PS pipeline section.
6. Limitations/robustness → validation table and limitations section.
7. Recommended next steps → evidence-contract-first next steps.
8. Further questions → signed bridge, k>12 decomposition, and outlier classification questions.
