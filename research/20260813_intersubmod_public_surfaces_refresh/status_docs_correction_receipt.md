<!--
建立時間: 2026-08-13
目標: 記錄公開資訊全面校正後 README、CURRENT_FOCUS 與文件索引的狀態同步及驗證結果
處理範圍: InterSubMod/README.md、InterSubMod/README.zh-TW.md、InterSubMod/docs/CURRENT_FOCUS.md、InterSubMod/docs/README.md、InterSubMod/docs/experiments/INDEX.md
關聯檔案:
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/00_INDEX.md
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/remote_state_receipt.md
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/github_correction_receipt.md
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/pages_correction_receipt.md
-->

# Status docs correction receipt — 2026-08-13

## 結論

五個狀態入口已同步到同一套 2026-08-13 權威敘述：158 個 claim families 的 verdict、58 個
問題的 P0/P1/P2 分級、scientific claim ceiling、精確資料量、Pages 本地 QA 與 CCU no-change
邊界均一致。GitHub About 明列為 `RESOLVED_LIVE`；default `main`、Wiki、Pages 明列為待發布，
沒有把本地 source correction 說成 live deployment。

本次服務 **G4（跨入口一致性與可重現性）／G5（可由外部收據驗證）**。未修改 canonical 科學
結果、演算法、HTML 報告、CCU source／live site，也未 commit、push 或 deploy。

## 修改清單

| 檔案 | 修改 | 驗收重點 |
|---|---|---|
| `InterSubMod/README.md` | 更新 GitHub surfaces 狀態並連結 refresh cycle／remote receipt | About live resolved；default main／Wiki／Pages pending |
| `InterSubMod/README.zh-TW.md` | 與英文 README 同步 publication boundary | EN／ZH 語意一致 |
| `InterSubMod/docs/CURRENT_FOCUS.md` | 新增 2026-08-13 全面校正主區段；將 2026-08-10 標為歷史快照；限縮 2026-08-05 容量敘述；連結最終 HTML/MD | 母體、claim ceiling、精確 bytes、live 狀態、CCU no-change 全列出 |
| `InterSubMod/docs/README.md` | 新增 2026-08-13 稽核、修正循環、P0 與最終 HTML＋SVG 報告入口 | 新入口目標皆存在 |
| `InterSubMod/docs/experiments/INDEX.md` | 新增 2026-08-13 條目與最終報告；移除不可重現的 `<300ms` performance claim | OpenMP 只保留 implementation status，不宣稱 runtime／speedup 數字 |

## 固定數值與判讀邊界

- Claim-family verdict：`69 + 31 + 28 + 26 + 4 = 158`。
- 問題集合：`28 + 26 + 4 = 58`，優先級為 `P0 34 + P1 20 + P2 4`。
- Pages 本地 QA：17 頁、37 個 inline SVG、68 個 page-profile checks、0 failure。
- 7 個 current compressed sidecars：`6,256,168,164 bytes = 5.826510641724 GiB`。
- HCC1395 tumor BAM：`283,071,595,503 bytes = 263.63096712436527 GiB`。
- 歷史 tagged-BAM `1.67 TiB` 缺 committed per-file bytes/hash registry，只能作 then-recorded
  snapshot，不是目前容量證據，也不得推導節省倍率。
- Scientific ceiling：direct evidence 到 same physical molecule；結構到 local、recurrence-allowed、
  model-conditional mutation-state candidate；甲基到 pattern-conditioned association。Cellular
  clone／lineage／causality 仍未被確認。

## Step → Verify

1. 同步五個狀態入口 → 驗證：指定 verdict、priority、Pages QA、exact bytes、About、CCU 與
   benchmark boundary 八類 required anchors 全存在，`required_missing=0`。
2. 排除過時狀態文字 → 驗證：舊的「About／Wiki／Pages 全待外部處理」與
   `single Region <300ms` pattern 搜尋結果為 `stale_status_patterns=0`。
3. 驗證 claim correction 未回歸 → 驗證：P0 registry `34/34`、P1/P2 registry `24/24`，兩者
   `errors=0`、`verdict=PASS`。
4. 驗證連結與格式 → 驗證：9 個本輪新增內部 target 全存在；`git diff --check` exit 0。

## 執行命令與實際輸出

### P0 claim guard

輸入：`InterSubMod/research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py`

```bash
python3 research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py
```

輸出：stdout。

```text
inventory_p0=34 registry_p0=34 checked_target_rules=77
required_anchors=140 forbidden_anchors=79 errors=0 verdict=PASS
```

### P1/P2 claim guard

輸入：`InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/scripts/validate_public_p1_p2_claims.py`

```bash
python3 research/20260813_intersubmod_public_surfaces_refresh/scripts/validate_public_p1_p2_claims.py
```

輸出：stdout。

```text
inventory_problem_p1_p2=24 registry_problem_p1_p2=24 checked_target_rules=26
required_anchors=80 forbidden_anchors=50 errors=0 verdict=PASS
```

### 新增內部連結

輸入：五個 status docs 內本輪新增的九個 local targets。

```bash
python3 - <<'PY'
from pathlib import Path
targets = [
    'research/20260813_intersubmod_public_surfaces_refresh/00_INDEX.md',
    'research/20260813_intersubmod_public_surfaces_refresh/remote_state_receipt.md',
    'research/20260813_intersubmod_public_surfaces_refresh/20260813_CCU教學站重點改進清單_01.md',
    'docs/reports/validated/2026/08/20260813_InterSubMod公開資訊更新與CCU改進清單_01.md',
    'docs/reports/validated/2026/08/20260813_InterSubMod公開資訊更新與CCU改進清單_01.standalone.html',
    'docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.md',
    'docs/reports/validated/2026/08/20260813_公開文件P0修正與驗證_01.md',
    'research/20260813_public_docs_p0_correction/00_INDEX.md',
    'research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md',
]
missing = [p for p in targets if not Path(p).is_file()]
print(f'new_internal_link_targets={len(targets)} missing={len(missing)}')
raise SystemExit(bool(missing))
PY
```

輸出：stdout。

```text
new_internal_link_targets=9 missing=0
```

### 格式與 stale-string guard

輸入：五個 status docs。

```bash
git diff --check -- README.md README.zh-TW.md docs/CURRENT_FOCUS.md \
  docs/README.md docs/experiments/INDEX.md
rg -n 'Main-repo source corrected locally|Main-repo source 於 2026-08-13 完成本地修正|About、Wiki、Pages 與 default branch 發布|About, Wiki, Pages and default-branch publication|單 Region\s*<\s*300ms' \
  README.md README.zh-TW.md docs/CURRENT_FOCUS.md docs/README.md docs/experiments/INDEX.md
```

輸出：`git diff --check` exit 0；`rg` 無命中，wrapper 回報如下。

```text
stale_status_patterns=0
```

## 最終判定

**PASS — status documentation correction completed locally.** 五個指定入口已採用一致且可稽核的
2026-08-13 狀態；About 的 live resolution 與仍待發布的 default main／Wiki／Pages 已分開。
下一個 publication gate 仍需明確授權後 merge／push／publish，再以 remote refs 與 live bytes
重驗；本收據不把該外部動作標為已完成。
