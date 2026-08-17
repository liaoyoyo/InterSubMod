<!--
建立時間: 2026-08-13
目標: 保存 GitHub Pages P0 敘述校正的逐 claim before/after/evidence/disposition、驗證命令與外部發布邊界
處理範圍: InterSubMod/docs/explain/*.standalone.html；2026-08-12 inventory 中具有 Pages occurrence 的指定 P0 claims，另含分母與 methyl-phasing-assist 防誤讀補強
關聯檔案: InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv; InterSubMod/research/20260813_public_docs_p0_correction/scripts/p0_claim_registry.json; InterSubMod/research/20260813_public_docs_p0_correction/claim_guard_receipt.md
-->

# GitHub Pages P0 敘述校正收據

用 Claim–Evidence–Disposition：**Pages P0 source 校正完成 — 指定 22 個 claim 中 21 個有 Pages occurrence 並已逐處修正；C150 在 Pages 無 occurrence，明列 N/A 而非略過**（影響：高；本地 source 語意信心：高；遠端發布信心：未發布）。

## 1. 任務邊界與結論

- 任務類型：**(B) Comprehensive validation 的 GitHub Pages P0 子範圍**；指定 Pages occurrence 無抽樣。
- 服務目標：G2（工具邊界）、G3（read-level epigenetic claim ceiling）、G4（可重現與版本界線）、G5（可被外部稽核的公開文件）。
- 本地結果：15 份 `InterSubMod/docs/explain/*.standalone.html` 已修正；`P009`、`P013` 無本批指定 P0 mutation，因此未改。
- 科學邊界：本批修正文案與圖例，不把文件修正冒充新的生物學實驗。權威 ceiling 仍是「local recurrence-allowed minimum mutation-state candidate arborescences」；禁止 confirmed cellular clone/lineage、CN/LOH-corrected CCF、methylation-rescored topology。
- 發布邊界：**沒有 commit、push、GitHub Pages deploy 或 live URL 修改**；本地 source PASS 不等於遠端頁面已更新。

## 2. 權威輸入

1. Claim inventory：
   `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv`
2. Source scope：
   `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/source_scope.tsv`
3. 模型與 claim ceiling：
   `InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json`，尤其 `model_contract` 與 `claim_boundary`。
4. 分母 registry：
   `InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv`
5. CLI effective default：
   `InterSubMod/include/utils/ArgParser.hpp:181-194`；`Config.hpp` initializer 位於 `InterSubMod/include/core/Config.hpp:40`，但會被 parser 清空後於無參數時補入 NHD。

## 3. 逐 claim 校正矩陣

Disposition `SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED` 表示本地 Pages source 已符合 bounded wording，但遠端 GitHub Pages 仍需另行發布及 live-byte re-fetch。

| Claim | 校正前問題 | 校正後 Pages 敘述 | 主要證據 | Disposition |
|---|---|---|---|---|
| C043 | 把 InterSubMod 與 LongLineage 一概寫成不能輸出 tagged BAM | InterSubMod 不寫 BAM；LongLineage public main `5daf50f` 不寫，但 feature `b9aaa12` 有 `longlineage-tag-bam` | inventory C043；版本邊界實查 | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C059 | `significance_summary.csv` 寫成 193 欄且無 schema version | commit `73afaea` 口徑為 199 欄；col187 `VerificationSchemaVersion=2`、col196 `RegionStratificationSchemaVersion=1`；仍無單一 whole-file layout version | `InterSubMod/src/core/RegionProcessor.cpp:1880-1968`；inventory C059 | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C109 | 把「同 cell lineage」寫成直接觀測 | 直接觀測止於同一 physical molecule；cellular co-membership/lineage 是 model-dependent inference | `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md:39-49` | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C112 | 直接稱 tree/reconstruction，缺 local、candidate、recurrence-allowed、model-conditional | 統一為 local recurrence-allowed minimum mutation-state candidate arborescence/topology；不是 biological lineage tree | `authority_manifest.json:63-97` | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C124 | 同一樣本三個 allele-PERMANOVA `p=0.001` 稱為 methylation corroboration 的「硬證據」 | 改為 strong locus-specific allele-associated methylation evidence；同一生物樣本，非 independent clone validation | inventory C124；`authority_manifest.json:76,85,88-95` | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C125 | normal-ASM-screen-negative 被寫成 tumor-acquired subclone methylation | 改為 tumor-associated、normal-ASM-screen-negative methylation；陰性 screen 不證明 acquisition 或 clone identity | inventory C125；matched-normal audit | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C126 | 「clean subset 證明 subclone methylation 真正存在」 | 改為在目前 controls 下相容的 lineage/mutation-state-associated methylation；保留未排除 confound | inventory C126 | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C127 | 多頁把區域結果稱 regional subclonal states | 改為 regional molecular-state candidates；僅在歷史主張被引用時保留 subclone 字樣並明示否定/校正 | `authority_manifest.json:78-97` | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C128 | 以 SEQC2-TVAF 近似 CCF 梯度佐證 subclone | 標為 historical case-specific approximate TVAF transform；noncanonical、未校正 CN/LOH/purity、不是 clone corroboration | `authority_manifest.json:83,91` | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C129 | 把一條 long read 等同一個 tumor cell chromosome | 一條 read 取樣一個 chromosome copy 的一個 DNA molecule；originating cell 沒有 barcode、不可知 | claim contract；inventory C129 | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C130 | 甲基被稱為對同一 subclone split 的獨立/正交印證 | 改為 mutation-defined genetic grouping 後的 bounded concordant association；同一樣本/reads，非獨立 | `authority_manifest.json:76,85` | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C131 | 宣稱得到 correct tag | 只輸出帶不確定性的 candidate molecular-state label；遺傳標記與尚未驗證的 methyl rescue 分開 | inventory C131；claim contract | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C143 | methyl/filter 頁把 production default 寫成 BERNOULLI | effective no-argument default = NHD；BERNOULLI 只在明示要求時使用 | `InterSubMod/include/utils/ArgParser.hpp:181-194` | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C144 | glossary/core 也把 BERNOULLI 當實際 production default | 01–03 三頁同步：`Config.hpp` 只是 initializer，CLI effective behavior 是 NHD | `ArgParser.hpp:181-194`; `Config.hpp:40` | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C148 | `verify_pipeline_numbers.py` 被說成驗證方法文件每個數字 | 限定為 historical 35,332-site pipeline verifier；不涵蓋 exact-PS、LongLineage、storage、code/test counts | inventory C148；verifier 實際輸出範圍 | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C149 | 公開 clone 後可用 tracked paths 在約 10 分鐘得到結果 | 顯著標 `PUBLIC REPRODUCTION BLOCKED` / `INTERNAL_ONLY`：公開 Git objects 缺 BAM/BAI、FASTA/FAI、VCF；列出 fixture、license、checksum、CI preflight 補件 | inventory C149 | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C150 | blanket assurance 出現在 A003/A004，inventory 無 Pages occurrence | Pages 不製造替代 occurrence；此 claim 由 entrypoint correction stream 處理，Pages disposition 明列 N/A | inventory C150 occurrences | **N/A_PAGES**；非遺漏 |
| C151 | 說 `inter_sub_mod` 產出 paper core exact-PS funnel | `inter_sub_mod` 產 per-region methylation/statistics；exact-PS funnel 來自獨立 research solver + Python runners | `authority_manifest.json` artifacts；inventory C151 | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C152 | 說 LOH 使 read 差異必然 somatic、非 error | LOH 只降低 simple two-parent-haplotype 解釋；sequencing error、CN/purity、stochastic methylation 等仍需守門 | inventory C152；normal error/homopolymer observations | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C153 | 0 observed alpha/beta co-alt 被說成只能是 sister branches | 改為與 branching candidate 相容，但 recurrence、coverage、error 下非唯一可識別；遠距 blocks 無直接 bridge | `authority_manifest.json:65,72-76`；inventory C153 | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C154 | 甲基被說成可知 subclone 關閉哪些 genes | 只描述 local methylation state；gene regulation/function 需 expression 或 functional assay | inventory C154 | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |
| C155 | 顯示 methylation 已 rescue untagged reads 並得到 correct final tag | 移到 unvalidated future hypothesis；目前不自動補標，需獨立 per-read truth、預先固定 classifier 與 error-rate benchmark | inventory C155；claim contract | SOURCE_FIXED / EXTERNAL_PUBLISH_REQUIRED |

### 額外防誤讀補強

- `P011` funnel 同時顯示兩個不可混用的分母：`170,131 / 255,752 strict components = 66.52%`；相對 `469,849 dataset records = 36.21%`。依據 `denominator_registry.tsv:2-5`。
- `index` 的 `methyl-phasing-assist` 卡改為未驗證研究假說；`88.5%` 僅屬模擬條件回推，真實 unphase reads 的困難區分布不可外推成 production rescue 成功率。
- `P008` 與 `P010` 的舊 CCF 顯示補上「歷史 TVAF 近似轉換」界線；不再暗示 corrected CCF 或 clone corroboration。

## 4. 修改檔案

1. Effective default：`InterSubMod/docs/explain/01_background-glossary.standalone.html`、`02_ism-core.standalone.html`、`03_methylation-read-filter.standalone.html`
2. chr2 claim ceiling：`04_subclone-reconstruction-chr2-18M.standalone.html`、`05_subclone-correction-audit-chr2-18M.standalone.html`、`06_ism-subclone-pipeline-concept.standalone.html`、`07_subclone-judgment-workstation-chr2-18M.standalone.html`、`08_subclone-logic-chain-chr2-18M.standalone.html`、`10_ism-cpp-vs-chr2-subclone-capability.standalone.html`
3. 系統、schema、版本：`11_system-map-overview.standalone.html`、`12_intersubmod-io.standalone.html`、`14_upstream-data.standalone.html`
4. 驗證與重現性：`15_python-html-layer.standalone.html`、`16_how-to-run.standalone.html`
5. 導覽與研究狀態：`InterSubMod/docs/explain/index.standalone.html`

頁 04 與 index 的最後硬化由 parent stream 在共享 worktree 協調完成；本收據驗證的是**最終 bytes**，未覆蓋該 stream 的修改。

## 5. Step → Verify 與實際輸出

### R01 — 34-P0 fail-closed source guard

- 輸入：claim inventory、machine-readable registry、27 份本地 README/Wiki/Pages/generator source。
- 命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py
```

- 實際輸出片段：

```text
inventory_p0=34
registry_p0=34
checked_target_rules=77
required_anchors=140
forbidden_anchors=79
unique_documents=27
errors=0
verdict=PASS
```

- 判定：本批 Pages checks 與其他公開 source 一起通過；guard 仍把外部發布面鎖為 `external_action_required`。

### R02 — Standalone HTML/SVG 靜態 QA

- 輸入：`InterSubMod/docs/explain/*.standalone.html`（17 份）。
- 命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 tools/explain_page_qa.py docs/explain/*.standalone.html
```

- 實際輸出：`總計 0 個 FAIL`。
- 注意：01、10、index 無 SVG；02–09 的既有 SVG 多數缺 `<title>/<desc>`，工具列為 WARN 而非 FAIL。這些 accessibility debt 沒有冒充為已修完。

### R03 — Strict parser、ID、details、anchors、JSON

- 輸入：相同 17 份 HTML。
- 命令：以 `lxml.html.HTMLParser(recover=False)` 解析，檢查 duplicate IDs、每個 `details` 恰有一個 `summary`、local targets/fragments、embedded JSON。
- 實際輸出：

```text
files 17 parsed 17 ids=167 details=92 svg=37 json=1
PASS lxml HTML parse; unique IDs; details/summary; local anchors; embedded JSON
```

### R04 — Standalone resource dependency

- 輸入：17 份 HTML 的 `script/style/img/source/iframe/video` resource 與 CSS `url(...)`。
- 命令：lxml + regex 掃描外部 resource URL；一般 `<a href="https://...">` 引用連結不算 runtime dependency。
- 實際輸出：

```text
files=17 external_resource_dependencies=0
PASS: standalone pages have no external script/style/image/font resource dependencies
```

### R05 — Legacy phrase residual `rg`

- 輸入：指定 P0 Pages。
- 命令：

```bash
rg -n -i 'production default.{0,30}BERNOULLI|BERNOULLI.{0,30}production default|regional operational subclonal states|CLEAN（tumor-acquired）|tumor-acquired subclone methylation|甲基.{0,24}(獨立印證|正交印證)|每一個數字都|方法文件裡的每個數字|約 ?10 ?分鐘|10 ?分鐘.*第一個' docs/explain/<P0-targets>
```

- 實際輸出：無命中；`rg exit=1`（預期的 no-match）。

### R06 — index methyl-phasing-assist semantic assertion

- 輸入：`InterSubMod/docs/explain/index.standalone.html`。
- 驗證：必須同時出現「沒有已驗證 per-read rescue/correct-tag truth」、「88.5% 僅模擬」與「不可外推 production rescue」，並拒絕正向 production claim。
- 實際輸出：

```text
index_required=3/3 forbidden_hits=0
PASS: methyl-phasing-assist is an unvalidated hypothesis; 88.5% is simulation-only and not production-generalizable
```

## 6. 已知 residual 與後續必要動作

1. **P007 generator drift 已關閉**：parent stream 已修正
   `InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/scripts/build_workstation_html.py`，
   實際重建 page-07 與 `fig4_corrected.svg`；page-07 為 3 SVG、0 FAIL。C125/C127/C128 已新增
   3 個 generator target rules（6 required、5 forbidden），終局 guard 仍 PASS。
2. **Browser QA**：本收據的早期 snapshot 只含 static/parser checks；總報告另以 Chromium 執行 desktop/mobile/print/no-external-runtime QA，結果以 final HTML QA receipt 為準。
3. **SVG accessibility debt**：page-07 的 3 SVG 已補 `<title>/<desc>`；全 17 頁共 37 SVG，其中其他 7 頁仍有 26 SVG 同時缺 `<title>/<desc>`，未把整站 accessibility debt 宣稱為歸零。
4. **Public fixture 缺失**：P016 已誠實標 `PUBLIC REPRODUCTION BLOCKED`；若要關閉 C149 的外部限制，仍須發布可授權 tiny fixture、index/checksum 與 CI smoke test。
5. **外部發布**：需 commit/push/deploy 後，重新抓取 17 個 live URL，逐檔比較 HTTP 200、content hash 與新 bounded anchors；完成前不能把 disposition 改為 `PUBLISHED_VERIFIED`。

## 7. Deployment action

`NONE` — 本子任務未 commit、未 push、未 deploy、未改 GitHub settings。遠端 GitHub Pages 仍是外部待辦。
