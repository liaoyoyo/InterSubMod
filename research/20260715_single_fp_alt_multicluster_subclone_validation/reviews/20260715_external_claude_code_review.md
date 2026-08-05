<!--
建立時間: 2026-07-15T05:30:00+08:00
目標: 保存外部 Claude Code 對單一 FP focal-ALT 甲基多群報告的完整 read-only adversarial audit
處理範圍: 主報告、report dataset、latest FP/TP/strict/topology/legacy/KB
關聯檔案: InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/reviews/20260715_external_claude_code_prompt.md
-->

# External Claude Code audit: full retained review

## 執行收據

- Claude Code executable：`/bip7_disk/liaoyoyo2001/.local/bin/claude`
- Claude Code version：`2.1.202`
- Model / effort：`opus` / `high`
- 第 1 輪：project plan mode；exit code `0`，只保留 summary，故不作正式完整稽核。
- 第 2 輪：`--safe-mode --permission-mode dontAsk --tools Read,Grep,Glob,Bash --disallowedTools Edit,Write,NotebookEdit`；exit code `0`。
- 第 2 輪沒有寫檔工具；內建 sandbox 因本機缺 `socat` 未啟用，完成後另以 git status / mtime 稽核外部修改。
- Prompt：`InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/reviews/20260715_external_claude_code_prompt.md`

## 1. Overall verdict

**VALID WITH CORRECTIONS**

主結論「僅支持 strict epigenetic follow-up candidate，不支持 high-probability subclone，亦不支持 linear evolution」可信且保守。14 個強制審查項目全部通過，未發現會翻轉 verdict 的數據錯誤。修正集中在呈現對稱性與 provenance，不影響結論。

## 2. Findings

### M1：both-evaluable stable endpoint 必須與 high-threshold 對稱揭露

- Source：`report_dataset.json -> matched_tp_specificity.both_evaluable`。
- Stable：FP `572/4,738 = 12.07%`、TP `469/4,738 = 9.90%`、paired `p=0.000599`。
- High threshold：FP `101/4,738 = 2.132%`、TP `96/4,738 = 2.026%`、`p=0.770`。
- 修正：主報告已補 stable 結果，並用同一 both-evaluable estimand 報告 6-sample 等權 `+0.419 pp`、leave-HCC1395 `+0.120 pp, p=0.902`。

### M2：legacy 21.55% 必須附可重算分子 / 分母

- Source：`InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/canonical_reference_v1/canonical_summary.json`。
- 重算：`687/3,188 = 21.5496%`。
- common-site legacy 另為 `670/3,302 = 20.29%`，兩者分母不同。
- 修正：主報告已明列 `687/3,188 evaluable` 與 source；仍只作 legacy sensitivity。

### M3：ALT readset 成員 100% 相同，不等於 reads.tsv 檔案 100% 相同

- `n_same_alt_readset=3,302/3,302`。
- `reads.tsv` byte-level identical `3,016/3,302 = 91.34%`；`286` 個不同。
- methylation matrix identical `3,302/3,302`；Bernoulli distance identical `2,081/3,302 = 63.02%`。
- 修正：主報告已區分 read member set、reads metadata、methylation 與 distance 四層。

### 高風險項目均通過

- 輸入確為同一次 `20260711_longphase_s...recalibrated.pass.vcf.gz`，truth 只在 producer 後 post hoc split。
- focal ALT 明確來自 `alt_support=ALT`，未以 HP tag 代替。
- 新 BAM 全在 observation workspace，沒有 canonical overwrite；sidecar missing `0`，MM+ML 完整。
- KB 語義一致：`1-1/2-1` 是第一 somatic haplotype，`1-2/2-2` 是第二 somatic haplotype；後者全體為 `0`。
- strict 10 明確標成 sensitivity intersection，不是 FDR/PPV、independent sites 或 confirmed clones。

## 3. Independent recalculation table

| # | 指標 | 報告值 | 外部重算 | 判定 |
|---|---|---:|---:|---|
| 1 | evaluable / all FP | 64.13% | `4,967/7,745 = 64.1317%` | PASS |
| 2 | stable / evaluable | 11.74% | `583/4,967 = 11.7375%` | PASS |
| 3 | residual / evaluable | 9.12% | `453/4,967 = 9.1202%` | PASS |
| 4 | high threshold / evaluable | 2.07% | `103/4,967 = 2.0737%` | PASS |
| 5 | strict / evaluable | 0.201% | `10/4,967 = 0.20133%` | PASS |
| 6 | high threshold FP-TP | +0.065 pp | `1.3299%-1.2653%=+0.0646 pp`, `p=0.77234` | PASS |
| 7 | leave-HCC1395 stable | -0.153 pp | `-0.15273 pp`, `p=0.77968` | PASS |
| 8 | explicit HP 1-2/2-2 | 0 | all/stable/high/strict 均為 `0` | PASS |
| 9 | topology linear Fisher | OR 1.217, p 0.441 | OR `1.21735`, p `0.44056` | PASS |
| 10 | legacy distance identity | 63.0% | `2,081/3,302 = 63.022%` | PASS |
| 11 | forced∩stable Jaccard | 0.328 | `364/(891+583-364)=0.32793` | PASS |
| 12 | 6-sample equal stable RD | +0.113 pp | `+0.11312 pp`；4 positive / 2 negative | PASS |
| 13 | coverage power | 6.02% -> 19.23% | `88/1,463`; `100/520` | PASS |
| 14 | per-sample denominator sums | FP 7,745 / eval 4,967 / stable 583 / high 103 / strict 10 | 全部逐列相加一致 | PASS |
| 15 | both-evaluable stable | 初稿未列 | FP 12.07% vs TP 9.90%, `p=0.000599` | 修正後 PASS |

Strict candidate 的 cluster sizes、AF、topology 在主報告、strict TSV 與 strict summary 三處一致；case catalog 的 HP-axis 與 technical-axis 反例亦一致。

## 4. Claim boundary audit

### Read-level methylation heterogeneity：可支持到 L2

同 pipeline 直接量測，並有 column / row-circular null、10-seed modal-K、assignment ARI。可稱 `focal-ALT read-level methylation heterogeneity candidate`。但 coverage power、shared reads 與相鄰位點使 prevalence 只能描述，不能把 site 當獨立 biological replicate。

### Subclone：不支持，只能 hypothesis-generating

- orthogonal cellular confirmation `0/7,745`
- explicit second somatic haplotype `0`
- high-threshold matched-TP specificity `p=0.772`
- tumor REF 也有 `112/450 = 24.89%` stable
- matched-normal methylation `NOT_TESTED`

最多稱 `strict epigenetic follow-up candidate`，不可稱 subclone。

### Linear / branching topology：結構上不可識別

單一 focal sSNV 只有 `ROOT -> ALT` trivial edge；focal linear identifiable `0/7,745`。高門檻沒有 linear enrichment，且同一 endpoint 同時出現在 linear 與 branching regional context。Same-pipeline tree 只能作 context。

## 5. Missing evidence

1. 同 region matched-normal methylation matrix，使用相同 row-circular null。
2. 至少 2 個獨立 sSNV / CN markers 的 group-specific genetic co-membership。
3. purity / CN-adjusted CCF 與 methyl group fraction 的一致性。
4. single-cell、colony、spatial 或 lineage barcode 的 cellular identity。
5. 要談 topology，另需多 mutation states 的共現 / 互斥與 cellular prevalence constraints。

HCC1395 與 Dorado 是同一 biological sample；COLO829 `truth_bed=null`，均不能作未加 caveat 的 biological replication / truth confidence。

## 6. External reviewer final wording

> 在最新一次 LongPhase-S recalibrated `FILTER=PASS`、post hoc truth-split 的 7,745 個 truth-FP sSNV 中，約 11.74%（583/4,967 evaluable；Wilson 10.87-12.66%）的 focal-ALT readset 呈現 null-gated、seed-stable 的 read-level 甲基異質性。此異質性不具 truth-FP specificity：高門檻 endpoint FP 1.330% vs matched TP 1.265%（p=0.772）。經 row-circular null、tumor-REF 背景與 assignment 穩定性交集，僅存 10 sites / 9 unique ALT readsets，屬保守 sensitivity 交集，非 FDR/PPV、非 confirmed clone。所有 focal single-site linear-topology identifiable 與 orthogonal subclone confirmation 都是 0，同一 endpoint 可落在 linear 或 branching regional context。因此只能稱為需正交驗證的 read-linked epigenetic heterogeneity candidates；缺乏 matched-normal methylation、genetic co-membership、CCF 與 cellular identity前，不能視為 high-probability subclone，更不能由單一 sSNV 推定 linear evolution。

**外部 reviewer 結論：未發現會翻轉 verdict 的數據錯誤。**
