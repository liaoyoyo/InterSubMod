# cis-ASM characterization pipeline — phase plan

> **定位**：系統性「找出 + 嚴格分析 + 分層輸出」cis-candidate ASM 位點（像 BRCA2）。
> **characterization only** — 不產 TP/FP discriminator、不餵 F1 filter（與 concluded-dead ASM-as-filter 區隔）。
> **Python-first**；C++ 改動只在重複需要時走 `/cpp-change` PDD。單樣本 HCC1395（cross-sample = future）。
> 完整計劃書：`InterSubMod/docs/plans/20260601_ISM_cis_ASM_characterization_pipeline_plan_01.md`

## 結構

```
pipeline/
├── loci.yaml                      # 位點 SoT（10 loci；role canonical|control|candidate）
├── lib/cis_asm_core.py            # 重構自 scripts 18/34/35（行為不變）：load_level1 / dbeta_axis / cis_test / cohesion_per_tag / score_locus
├── fixtures/brca2_level1.tsv.gz   # 凍結 BRCA2 Level-1（regression fixture，26,684 rows）
├── tests/test_brca2_regression.py # golden guard（7 checks；自含 + pytest 相容）
├── Makefile                       # make test / make features
└── 00_PIPELINE_PLAN.md            # 本檔
```

## 已驗證的 schema（trial + fix + PoC 三輪 value-first 收斂，2026-06-02/03）

- **保留**：per-CpG 5mC / 5hmC（各自 Δβ）、always-on 特徵（vaf/cn/dist_to_tss/mechanical_cis/gene/power_class）、HP-cis 特徵（dbeta_HP/d_cis/d_drift/cohesion，testability-gated 4/10）。
- **drop**：read-level 5hmC clustering（PoC 36 證實 5hmC ~6% 太稀疏 → cluster-ARI=0.04 雜訊）。
- **歸 P3 pooled**：去甲基化轉換 signature（Δ5mC vs Δ5hmC 負相關，單位點 underpowered）。
- **修正**：dbeta_5mC paired + window-restricted（ONT 長 read 不限窗會稀釋）；kstar/crosstag gate on n_som。

## Golden（凍結 2026-06-03，cross-extractor validated；改動需 [golden-update] commit）

| metric | golden | 說明 |
|--------|-------:|------|
| dbeta_HP | −0.122 | HP1 vs HP1-1 paired；MSA Level-1 == BAM-windowed（diff 0.000）|
| dbeta_5mC / 5hmC | −0.121 / −0.001 | ASM 純 5mC |
| d_cis / d_drift / p_cis | −0.1418 / −0.0217 / 0.0 | somatic 偏離 baseline，germline 不漂移 |
| sil_HP11 / sil_HP1 | 0.313 / 0.119 | somatic 子克隆最內聚 |

## Phase 路線（exit criteria）

| Phase | 內容 | 狀態 | exit |
|-------|------|------|------|
| **P0 地基** | pipeline/ + lib + loci.yaml + **BRCA2 regression test** | ✅ **完成** (99f052d) | `make test` PASS；故意 break(min-collapse)→7/7 FAIL ✓ |
| **P1 抽取 substrate** | `stages/stage_extract.py`（BAM-direct，mod_code 欄，window-restricted）+ Level-1-plus 持久化 cache + `lib` 加 load_level1_plus/dbeta_mod | ✅ **完成** | `test_brca2_bam_extract` PASS：5mC=−0.121 / 5hmC=−0.001 / **any=−0.122==Level-1(cross-extractor invariant)** / cache MISS→HIT ✓ |
| **P2a LOH/CNV + axis-gate** | `lib/genomic_context.py`：整數 CN（SEQC2 gain/loss_cn bed）+ loh_status + **axis_validity gate**（LOH→ALLELE primary，修 43% HP-in-LOH 違規）| ✅ **完成** | `test_genomic_context` PASS：BRCA2 CN=5/nonLOH/HP-valid；LOH→ALLELE ✓ |
| P2b cis-core | cis 欄 genome-wide + cis-ladder T0-T3 verdict + power-class | TODO | het-null 落 T0/T1；BRCA2 T3 |
| P3 scan + causation | two-stage scan + ranking + 證據卡 + mechanical-cis + (pooled 轉換 signature) | TODO | scan grep 到 BRCA2；het-null 落 T0/T1 |
| P4 持續驗證閉環 | git-diff hook（advisory exit-0）+ provenance manifest | TODO | 編輯 stage 觸發 advisory |

## git-diff 持續驗證迴圈

1. 改 `lib/` 或 `stages/` 前後跑 `make test`。
2. 數值若變 → diff 舊輸出 + commit message 標 `[golden-update]` + 理由；否則 golden 不動。
3. **cross-extractor 一致性**（Level-1 dbeta_HP == BAM dbeta_any windowed）= 一條 invariant；背離即抽取層出問題。

## Guardrails

characterization only（不重開 ASM-as-filter / TO-germline-FP）· Python-first · 磁碟寫 `/` 非 big7/big8 · golden 來自實測非捏造。
