<!--
建立時間: 2026-08-13
目標: 記錄 README、Quickstart 與 Project Summary 的 P0 公開敘述修正
處理範圍: InterSubMod/README.md、README.zh-TW.md、QUICKSTART.md、README_PROJECT_SUMMARY.md
關聯檔案:
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md
  - InterSubMod/research/20260813_public_docs_p0_correction/00_INDEX.md
-->

# Repository entrypoint P0 correction receipt

## 結論

Task type **B — Comprehensive validation continuation**；服務 G2／G3／G4／G5。入口文件涉及的
17 個 P0 claim families 已在本機 source 給出逐項 disposition：16 個完成 local correction，
C108 GitHub About 因不屬於 tracked files，保留可直接套用文案並標
`EXTERNAL_ACTION_REQUIRED`。本輪不 commit、push、merge、修改 About 或部署。

這裡的 `LOCAL_SOURCE_CORRECTED` 只代表 working-tree source 已修，**不代表** default branch 或
GitHub live surface 已更新；C135、C136、C141 仍須在同一個 remote change 原子合併英文與中文
README，才能轉為 live resolved。

## Step → Verify

1. 盤點入口 P0 → 驗證：16 個 ID 均有 before／after／evidence／disposition。
2. 修正科學 claim ceiling → 驗證：read 只被稱為 physical molecule；cell／lineage 明列為模型推論。
3. 修正數字、schema、版本與 executable provenance → 驗證：分母重算吻合，版本與 solver 路徑可追。
4. 修正 quick start → 驗證：公開命令不再引用缺失 fixture；內部 wrapper 明標 `INTERNAL ONLY`。
5. 執行 target-only QA → 驗證：Markdown links、residual scan、`git diff --check` 皆須 exit 0。
6. 保留發布邊界 → 驗證：remote About／default branch／Wiki／Pages 不以本地變更冒充 resolved。

## 輸入、命令、輸出

### 輸入路徑

- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md`
- `InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json`
- `InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv`
- `InterSubMod/README.md`
- `InterSubMod/README.zh-TW.md`
- `InterSubMod/QUICKSTART.md`
- `InterSubMod/README_PROJECT_SUMMARY.md`

### 執行命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod

python3 scripts/analysis/check_md_links.py \
  README.md README.zh-TW.md QUICKSTART.md README_PROJECT_SUMMARY.md

awk 'BEGIN {
  printf "170131/255752=%.4f%%\n", 100*170131/255752;
  printf "170131/469849=%.4f%%\n", 100*170131/469849;
  printf "63506/71955=%.4f%%\n", 100*63506/71955
}'

rg -n -i \
  '重建腫瘤亞群譜系|same cell lineage.*direct observation|兩支引擎都不能寫 BAM|265 tests / 38|193 columns|193 欄|每個數字、指令與檔案格式|every number, command and file format|data/bam/HCC1395/tumor.bam|data/ref/hg38.fa|one_snv.vcf' \
  README.md README.zh-TW.md QUICKSTART.md README_PROJECT_SUMMARY.md

git diff --check -- \
  README.md README.zh-TW.md QUICKSTART.md README_PROJECT_SUMMARY.md
```

### 輸出路徑

- `InterSubMod/README.md`
- `InterSubMod/README.zh-TW.md`
- `InterSubMod/QUICKSTART.md`
- `InterSubMod/README_PROJECT_SUMMARY.md`
- 本收據：`InterSubMod/research/20260813_public_docs_p0_correction/entrypoint_correction_receipt.md`

## Claim-by-claim correction matrix

| ID | Before（問題） | After（本機 source） | 依據 | Disposition |
|---|---|---|---|---|
| C043 | 兩支引擎一概不能寫 tagged BAM。 | 分列 `inter_sub_mod` 不寫 BAM、LongLineage main `5daf50f` 無 writer、feature `b9aaa12` 有 `longlineage-tag-bam`。 | Clean revision comparison／audit R16–R18。 | LOCAL_SOURCE_CORRECTED；live publish pending |
| C059 | current summary 被寫成 193 欄且無 schema。 | 釘 `73afaeac` 為 199 欄，列 `VerificationSchemaVersion=2` 與 `RegionStratificationSchemaVersion=1`，並說明尚非全檔 layout version。 | Fresh source/runtime header。 | LOCAL_SOURCE_CORRECTED |
| C107 | 首頁標題宣稱 subclonal reconstruction。 | 改為 single-molecule sSNV co-occurrence 的 local mutation-state candidate reconstruction。 | Authority claim boundary。 | LOCAL_SOURCE_CORRECTED；default-branch publish pending |
| C108 | GitHub About 宣稱 subclone resolution。 | 建議 About：`Read-level somatic-mutation co-occurrence and bounded methylation association for local mutation-state candidate analysis in ONT data.` | Authority claim boundary。 | EXTERNAL_ACTION_REQUIRED；未修改 live About |
| C109 | same molecule 被說成 same cell lineage 的直接觀察。 | 直接觀測只限 physical molecule；cellular co-membership／lineage 明列 model-dependent。 | Claim contract v5。 | LOCAL_SOURCE_CORRECTED |
| C112 | tree／reconstruction 未加 candidate、local、recurrence-allowed 限定。 | 首次定義與方法表統一用 local recurrence-allowed candidate arborescence／mathematical topology signature。 | Authority model contract。 | LOCAL_SOURCE_CORRECTED |
| C113 | 170,131 被標成一般 single sites，grain 隱藏。 | 明列 `170,131 / 255,752` strict read-linkage components，單位為 `k=1` component。 | Denominator registry。 | LOCAL_SOURCE_CORRECTED |
| C134 | README blanket 宣稱每個數字、命令、格式都已驗證。 | 改為 artifact-level table：identity/date、command/result、scope、known failure。 | Audit C145–C151 與 command receipts。 | LOCAL_SOURCE_CORRECTED |
| C135 | Default README 把整合輸出說成已解析 cellular subclone。 | EN/ZH source 都採 bounded candidate framing。 | Confirmed cellular subclone count=0。 | LOCAL_SOURCE_CORRECTED；須原子 merge 至 default branch |
| C136 | 甲基化 read clusters 被說成不同 cell populations。 | 改為 genetic grouping 後的 pattern-conditioned association；不作 cell identity。 | Authority methylation role。 | LOCAL_SOURCE_CORRECTED；default-branch publish pending |
| C137 | 宣稱 C++ 自動產 heatmap。 | 分開 C++ per-region data 與 downstream Python rendering。 | Fresh minimal-run artifact census。 | LOCAL_SOURCE_CORRECTED |
| C138 | 把含 `/big7`／`/big8` 的 wrapper 當 public standard command。 | Quickstart 明標 `INTERNAL ONLY`；public analysis command 使用顯式自備 inputs。 | Public Git objects 缺必要 inputs。 | LOCAL_SOURCE_CORRECTED；portable fixture 仍缺 |
| C139 | Project Summary 稱一鍵從 BAM 到圖。 | 改稱內部 orchestration，列 hard-coded path 與 public fixture 缺口。 | Script source audit。 | LOCAL_SOURCE_CORRECTED |
| C140 | Cluster heatmap 被稱為 subclone groups。 | 改為 read-level methylation clusters，明示不是 cellular subclones。 | Current claim boundary。 | LOCAL_SOURCE_CORRECTED |
| C141 | 宣稱雙語，但 default branch 中文 endpoint 404。 | 本機 EN/ZH 同步修正並互連；不可在尚未 merge 時宣稱 endpoint 已修。 | Locked main raw endpoint 404。 | LOCAL_SOURCE_CORRECTED；atomic remote merge required |
| C142 | 把 66.52% 套到全部 mutations。 | 明列 `170,131/255,752=66.5219%` 是 strict components；對 469,849 records 為 36.2097%。 | Independent denominator recomputation。 | LOCAL_SOURCE_CORRECTED |
| C150 | 另一處 blanket verification assurance。 | 改成 bounded verification scope 與 per-artifact receipt links。 | Public audit contradiction set。 | LOCAL_SOURCE_CORRECTED |

## 驗證結果

- `check_md_links.py`：exit 0、stdout 空，代表 4 個入口文件的 Markdown links 無錯。
- 分母重算：

  ```text
  170131/255752=66.5219%
  170131/469849=36.2097%
  63506/71955=88.2579%
  ```

- 舊錯誤句型 residual scan：0 命中；`rg` exit 1 是預期的 no-match 狀態。
- `git diff --check`：exit 0、stdout 空。
- C++／CMake build inputs 自 audited core `73afaeac` 起未變；2026-08-13 重跑既有 build artifact：
  GoogleTest 270 tests／39 suites 全過（exit 0），CTest 270/270 全過（exit 0，7.55 s）。
- 共享 branch 執行中由 `83741469` 前進至 `95d420f6`；唯讀 conflict scan 顯示兩個新增
  commit 對本收據四個 target files 與 C++ core build inputs皆無差異，因此未 reset、未覆寫。

## 尚未完成的外部動作

1. GitHub About 套用 C108 建議文字。
2. 英文／繁中 README 與 Quickstart／Summary 原子合併至 default branch，確認 raw endpoints HTTP 200。
3. Wiki source 發布到 `InterSubMod.wiki.git`；Pages source 經正式 deploy 後核對 live bytes／commit。
4. 若要把 public quick start 升為 end-to-end acceptance test，需發布有授權的 BAM／FASTA／VCF tiny fixture、checksum、portable wrapper 與 output hash。
