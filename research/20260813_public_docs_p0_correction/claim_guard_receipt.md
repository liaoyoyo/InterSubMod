<!--
建立時間: 2026-08-13
目標: 保存公開文件 P0 claim guard 的輸入、命令、實際輸出與 fail-closed 驗證證據
處理範圍: 2026-08-12 claim inventory 的 34 個 P0；本地 README/Wiki/Pages source；外部發布狀態只登記不變更
關聯檔案: InterSubMod/research/20260813_public_docs_p0_correction/scripts/p0_claim_registry.json; InterSubMod/research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py
-->

# P0 claim guard 執行收據

## 結論

本地 source guard 最終 **PASS（exit 0）**：原始 inventory 的 34/34 個 P0 都有唯一 disposition；其中
33 個有本地 source 規則，C108 GitHub About 是純外部動作。Guard 實際掃描 27 份文件、77 個
claim-target 規則、140 個正向 anchor 與 79 個舊錯誤 residual；錯誤數為 0。新增的 3 個
target 規則直接守住 page-07 generator，避免重新產生時回寫舊語意。

這個 PASS **不等於已發布**。GitHub About、default branch、GitHub Wiki 與 GitHub Pages 四個公開面
仍全部標為 `external_action_required`；必須另行 push/deploy/edit 並重新抓 live bytes，才能關閉公開狀態。

## 輸入與輸出

- 權威 P0 輸入：
  `/big7_disk/liaoyoyo2001/InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv`
- Machine-readable registry：
  `/big7_disk/liaoyoyo2001/InterSubMod/research/20260813_public_docs_p0_correction/scripts/p0_claim_registry.json`
- Validator：
  `/big7_disk/liaoyoyo2001/InterSubMod/research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py`
- Version-control boundary：`InterSubMod/.gitignore` 對此單一 registry 加精確 negation；
  `git status --short --untracked-files=all` 可見該檔，無須靠未記錄的 `git add -f`。
- 掃描範圍：repo 內 27 份 README／Quickstart／Summary／Wiki／Pages／page-07 generator source；逐檔路徑與規則見 registry。
- 輸出：stdout JSON summary、stderr 錯誤清單、process exit code；不寫入或改變被掃描文件。

## Step → Verify

1. Registry 反查 inventory P0 集合
   → 驗證：`inventory_p0=34`、`registry_p0=34`，缺少或多出 claim ID 直接 exit 1。
2. 驗證 disposition 與外部狀態
   → 驗證：33 個 `local_source_fixed_external_publish_required`、1 個
   `external_action_required`；四個公開面全部保留外部動作狀態。
3. 掃描本地 source
   → 驗證：每個 local rule target 必須存在；每個正向 regex 必須命中；任何舊錯誤 residual 命中即失敗。
4. 做負向 mutation test
   → 驗證：只在記憶體移除 C155 後，guard 明確回報 `registry missing inventory P0 IDs: C155` 並 exit 1。

## R01 — 語法、JSON 與正向 guard

### 執行命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 -m py_compile \
  research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py
python3 -m json.tool \
  research/20260813_public_docs_p0_correction/scripts/p0_claim_registry.json >/dev/null
python3 research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py
```

### 實際輸出片段

```json
{
  "checked_target_rules": 77,
  "errors": 0,
  "external_actions": 35,
  "external_only_claims": 1,
  "external_public_surfaces": 4,
  "forbidden_anchors": 79,
  "inventory_p0": 34,
  "local_claims": 33,
  "registry_p0": 34,
  "required_anchors": 140,
  "unique_documents": 27,
  "verdict": "PASS"
}
FINAL_GUARD_EXIT=0
```

判定：語法與 JSON 均有效，正向 source scan exit 0。

## R02 — Fail-closed mutation test

### 執行命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py \
  --simulate-drop-claim C155
```

### 實際輸出片段

```text
ERRORS:
- registry missing inventory P0 IDs: C155
...
"inventory_p0": 34,
"registry_p0": 33,
"verdict": "FAIL"
MUTATION_GUARD_EXIT=1
```

判定：漏一筆 inventory P0 時確實 fail closed，不是裝飾性永真檢查。

終局負向摘要：`registry_p0=33`、`checked_target_rules=75`、`required_anchors=137`、
`forbidden_anchors=77`、`unique_documents=27`、`errors=1`、exit 1。

## 收斂紀錄

- 首次 source scan：`FAIL`，59 個錯誤。它同時抓到尚未完成的 Pages 修正與數個過度寬鬆的初版 regex。
- 將 residual 改成具語境範圍的舊錯誤片語，避免把「不支持 X」誤判成仍宣稱 X。
- Pages source 完成後，guard 尚保留 3 個真實 gap：頁 04 的「subclone 重建樹」殘留，以及頁 11
  未完整並列 66.52% strict-component 與 36.21% dataset-record 分母。兩處 source 修正後才轉為 PASS。
- 曾嘗試用 process substitution 做 mutation probe，但 `Path.resolve()` 取得已關閉的 `/proc/.../fd`，該次只驗到
  registry load failure，**不列為有效 fail-closed 證據**；R02 是取代它的有效純記憶體測試。

## 能與不能證明的事

Guard 能證明：指定 source 目前含必要的 bounded wording、沒有登記的舊錯誤片語、34 個 P0 無 registry 漏項，
且外部發布狀態沒有被 source edit 冒充為已完成。

Guard 不能單獨證明：底層生物結論為真、精確數字已重新運算、所有未登記同義 overclaim 都不存在，或遠端 GitHub
已更新。這些仍須由 frozen artifact／runtime receipt、人工科學審查與發布後 live re-fetch 分別驗證。
