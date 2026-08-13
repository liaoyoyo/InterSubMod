<!--
建立時間: 2026-08-13 11:30 +08:00
目標: 讓無對話背景的新AI在cold start不誤讀InterSubMod／LongLineage研究狀態
處理範圍: research handoff snapshot facts、authority與禁止推論
關聯檔案:
  - InterSubMod/docs/handoff/20260813_完整研究資料與軟體交接_01/00_INDEX.md
  - InterSubMod/docs/handoff/20260813_完整研究資料與軟體交接_01/registries/artifact_registry.json
驗證方式: reader acceptance test；所有數值回指frozen denominator registry
證據等級: context routing，非獨立science evidence
狀態: ACTIVE_HANDOFF_CONTEXT
-->

# AI Cold-Start Context

## 身分與任務類型

- Project：InterSubMod／LongLineage，癌症ONT long-read read-level研究。
- Snapshot：2026-08-13 research handoff，**不是production release**。
- Scope：Task B comprehensive inventory + D external handoff；7 technical datasets、6 biological IDs、chr1-22。
- Goals：G1/G3/G4/G5。

## 先讀順序

1. `InterSubMod/AGENTS.md`
2. `InterSubMod/docs/CURRENT_FOCUS.md`
3. 本目錄上一層`00_INDEX.md`
4. 本包`evidence/authority_manifest.json`與`evidence/denominator_registry.tsv`（2026-08-01 SoT exact copies）
5. 本包`registries/artifact_registry.json`、`registries/claim_registry.json`與`registries/authority_superseded_crosswalk.json`

## 不可改寫的事實

- `confirmed cellular subclone = 0`
- `confirmed linear ancestry = 0`
- methylation = association only；不能選edge、建立clone或當獨立clone validation。
- CN/LOH沒有整合進frozen reconstruction。
- `88.2579% = 63,506 / 71,955`，只是一種rooted-unlabeled graph-shape signature的model-conditional比例；不是accuracy、prevalence、clone fraction或biological correctness。
- HCC1395與HCC1395_DORADO是一個biological ID的兩個technical datasets。
- master manifest 19 file lines = 1 header +18 logical run rows；「19 canonical runs」已INVALIDATED。Physical inventory另為35 current＋16 pending archive＝51。
- `tree.nwk` leaves是reads；它是read dendrogram，不是cellular lineage。
- InterSubMod現行summary contract為199欄；按header/schema讀，不能用舊的59/114/117欄位置。
- effective no-argument distance metric是NHD。

## 軟體責任

- LongPhase-S：paired tumor-normal phasing/haplotagging/tumor-DNA-fraction context/recalibrated VCF。
- LongPhase-TO：tumor-only candidate phasing、PON/LOH/HI context、tagged reads。
- exact-PS／LongLineage：從exact PS×HP read assignments建local candidate mutation-state families／graph shapes／abstains。
- InterSubMod：BAM+FASTA+VCF → per-region methylation、read distance、read clustering與statistics。
- Python／HTML：呈現validated data；不重新計算或創造science。

## Git freeze

- InterSubMod include baseline：`ddd8909a838318d8a77969313e9561c8ff9d01c2`。
- InterSubMod exclude：`73afaeac`、2026-08-13 drilldown/CNV dirty work；標`IN_PROGRESS/PARTIAL`。
- LongLineage preview include：`b9aaa12`；該 immutable candidate 的 frozen baseline為47/47 tests。後續private safety stack `f60b5f3`為49/49 foundation PASS，但不是preview baseline，不能把49歸給`b9aaa12`。
- LongLineage exclude：`9ad976b`、`6ce62b2`、dirty work。
- LongLineage production `run` intentionally returns`KernelBlocked` exit 6；P3/P4/P5/P7/P8 BLOCKED。

## Finality演算法

只有artifact registry列為`AUTHORITY + FULL + FINAL_FOR_SCOPE`、且producer/hash完整者，才是其明示scope內的science/source-byte final；目前是19筆。另1筆`VALIDATED_DERIVED + FINAL_FOR_SCOPE`只對append-only provenance adjudication final，不是science authority。任何`PARTIAL/HISTORICAL/IN_PROGRESS/INVALIDATED`或`NON_FINAL/SUPERSEDED`不得升格。目錄名、mtime、HTML標題、最新commit都不是finality證據。

## 執行前檢查

1. 你是在bip7還是bip8？只信`hostname`與該host receipt。
2. site profile是否以local file提供，且doctor通過BAM/VCF/FASTA/index/tool hash？
3. 要跑的是SUPPORTED、REPRODUCIBLE_LEGACY或ARCHIVED workflow？
4. scope是FULL、PARTIAL或DEMO？DEMO不能寫入science ledger。
5. input artifact是否current，是否被superseded？
6. command、commit、schema、denominator、output與hash會寫到哪個receipt？

## 六問標準答案

1. **專案是什麼？** Read-level sSNV linkage、local mutation-state candidate reconstruction與methylation association研究。
2. **目前結論？** 有frozen technical/model evidence；沒有confirmed cellular clone或ancestry；methyl association-only；CN/LOH未整合。
3. **哪些資料final？** 查artifact registry，不看資料夾名；只有AUTHORITY+FULL+FINAL_FOR_SCOPE+hash。
4. **三套工具如何分工？** LongPhase上游tag/phasing；exact-PS/LongLineage candidate reconstruction；InterSubMod read-methylation/statistics。
5. **bip7/bip8如何跑？** bootstrap profile→doctor→plan→clean build/test→synthetic smoke→real-data read-only preflight；每台host各自receipt。
6. **如何驗證/繼續？** replay authority、跑registry/claim/link/HTML gates；新研究先pre-decision audit、固定scope/inputs/stop rules，另開PR。

## 回答時禁止

- 禁止將molecule/read-level結果改寫為cell-level truth。
- 禁止把88.2579%稱為準確率。
- 禁止把7 technical datasets稱為7 biological replicates。
- 禁止把feature branch能力稱為main能力。
- 禁止把local source corrected稱為live published。
- 禁止把bip7看到的NFS mount當作在bip8成功執行。
- 禁止從舊InterSubMod報告推論LongPhase外部工具的一般F1/CLI行為；先查KB/paper。

## 目前阻塞

- 158 public claims中33/33 local P0 guards已PASS；第34筆P0 `C108` GitHub About已由live re-fetch取得bounded live `CONFIRMED_WITH_LIMITS`。Registry現為69 `CONFIRMED`、69 `CONFIRMED_WITH_LIMITS`、20 `UNVERIFIED`；P1/P2、default branch、Wiki與Pages仍未閉合，publication/release維持`BLOCKED`。
- bip8 host-local receipt缺失。
- full-history secret scan工具/gate未閉合。
- LongLineage後續private safety stack `f60b5f3` foundation 49/49 tests已PASS（`b9aaa12` frozen baseline是47/47），但4 source rows、21 source-license rows、11 dependency license determinations與7 history findings仍blocked；repository目前private，曾公開暴露不能視為未發生。
- tag／GitHub Release不得建立。

本context只負責避免誤解；真正數字與狀態以registries及receipts為準。
