<!--
建立時間: 2026-08-13
目標: 保存 InterSubMod 公開文件 P1/P2 24-claim correction registry 的輸入、規則、正向 guard、負控制與 P0 回歸證據
處理範圍: 2026-08-12 claim inventory 中 verdict=NEEDS_CORRECTION|CONTRADICTED|UNVERIFIABLE 且 priority=P1|P2 的完整 24 個 claim families
task_type: B_comprehensive_validation
服務目標: G2/G3/G4/G5
狀態: LOCAL_GUARD_PASS_PUBLICATION_PENDING
關聯檔案:
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/scripts/p1_p2_claim_registry.json
  - InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/scripts/validate_public_p1_p2_claims.py
  - InterSubMod/research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py
-->

# InterSubMod P1/P2 public-claim guard 收據

## 結論

P1/P2 本地 source guard 最終 **PASS（exit 0）**。權威 inventory 的 **24/24** 個 problem
claim families 都有唯一 registry entry，其中 P1=20、P2=4；validator 實際掃描 12 份文件、
26 個 claim-target rules、80 個 required anchors 與 50 個 forbidden anchors，errors=0。

Disposition 為：

- `local_source_corrected_publication_pending`：21 個。
- `unsupported_public_claim_removed_publication_pending`：3 個（C048、C090、C157）。

24/24 都有至少一個 target check、非空 `required_all`、非空 `forbidden_any`、非空 evidence、
非空 external actions，且 `live_status` 全部維持 `publication_pending`。這個 PASS 只證明目前
working tree 的指定文字 contract，**不等於 GitHub default branch、Wiki 或 Pages 已發布**。

負控制在記憶體移除 C157 後如預期 **FAIL（exit 1）**，明確回報
`missing=['C157']`。既有 P0 guard 同時回歸 **PASS（34/34、errors=0）**，故這個 P1/P2
registry 沒有以犧牲 P0 contract 換取通過。

## 輸入與輸出路徑

### 輸入

1. `/big7_disk/liaoyoyo2001/InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv`
2. `/big7_disk/liaoyoyo2001/InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/github_correction_receipt.md`
3. `/big7_disk/liaoyoyo2001/InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/pages_correction_receipt.md`
4. `/big7_disk/liaoyoyo2001/InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/scripts/validate_public_p1_p2_claims.py`
5. Registry 所列的 12 份 README／Summary／Wiki／Pages source。

### 輸出

1. `/big7_disk/liaoyoyo2001/InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/scripts/p1_p2_claim_registry.json`
2. `/big7_disk/liaoyoyo2001/InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/claim_guard_receipt.md`
3. Guard stdout JSON、stderr errors 與 process exit code；validator 本身不改被掃描文件。

### Frozen checksums

```text
d5d8794c2467ac6ac466711f91890666f3dbc571d9b3e27b07ba5f10d642ddbc  claim_inventory.tsv
4097a0afcb8e44754832a24c804c349ceadce4188484636b37caf60e011bba23  validate_public_p1_p2_claims.py
a5c4022865124747dea073b1d54850c7f79a3d6e3e4572ebe9925513e3a1069c  p1_p2_claim_registry.json
```

Registry 大小：653 lines、26,323 bytes。

## Step → Verify

1. 從 inventory 重建 problem P1/P2 集合
   → 驗證：精確為 24 IDs；registry 必須集合相等，duplicate／missing／extra 皆 exit 1。
2. 驗證每 claim 的 metadata contract
   → 驗證：source verdict 與 priority 必須逐列等於 inventory；disposition 僅能為兩個允許值；
   `live_status` 必須是 `publication_pending`；evidence／external actions／checks 不得為空。
3. 掃描 corrected source
   → 驗證：26 個 target 都存在，80 個 bounded anchors 全命中，50 個舊 overclaim anchors
   全不命中。
4. 做 fail-closed 負控制
   → 驗證：`--simulate-drop-claim C157` 產生 registry=23、errors=1、verdict=FAIL、exit 1。
5. 做 P0 regression
   → 驗證：34/34 P0、77 rules、140 required、79 forbidden、errors=0、exit 0。

## 24-claim disposition

| Claim | Inventory | Local disposition | Guard target 與核心 contract | Live |
|---|---|---|---|---|
| C047 | UNVERIFIABLE / P1 | local source corrected | page14：exact sidecar bytes；tagged-BAM total `UNVERIFIED`，缺 per-file receipt | publication pending |
| C048 | CONTRADICTED / P1 | unsupported claim removed | page14：沒有同對象 exact bytes 就不報 287×／293.34× | publication pending |
| C089 | CONTRADICTED / P1 | local source corrected | Upstream Wiki：只確認 sidecar exact sum，不把 storage snapshot 稱 fully verified | publication pending |
| C090 | UNVERIFIABLE / P1 | unsupported claim removed | Upstream Wiki：九欄不含 SEQ/QUAL；無 byte census／utility audit | publication pending |
| C110 | NEEDS_CORRECTION / P1 | local source corrected | README：per-locus marginal VAF、無 linkage／extra assumptions 的條件式不可識別 | publication pending |
| C111 | NEEDS_CORRECTION / P1 | local source corrected | README：current single-bulk measurement set；缺 orthogonal data／extra assumptions | publication pending |
| C114 | NEEDS_CORRECTION / P1 | local source corrected | page11：flags 在 `cohort_receipt.json`、`summary/all7_summary.json`，非 authority top-level | publication pending |
| C115 | NEEDS_CORRECTION / P1 | local source corrected | README：`docs/explain` editorial upstream；Wiki manually synchronized derivative | publication pending |
| C116 | NEEDS_CORRECTION / P2 | local source corrected | README：deploy `fbdf7c7` 為 17 pages／37 inline SVG elements，非 37 semantic figures | publication pending |
| C117 | NEEDS_CORRECTION / P1 | local source corrected | page02：tested formulation negative；materially new hypothesis + pre-decision audit 才重啟 | publication pending |
| C118 | NEEDS_CORRECTION / P2 | local source corrected | page03：當時範圍與日期的 convergent negative tests，非 bulletproof | publication pending |
| C119 | NEEDS_CORRECTION / P1 | local source corrected | page01：paired calling 通常改善但仍誤判，須 truth-set benchmark | publication pending |
| C120 | NEEDS_CORRECTION / P1 | local source corrected | page02：1.8×drift+0.02 heuristic；cis-compatible candidate；HCC1395-only completeness | publication pending |
| C121 | NEEDS_CORRECTION / P1 | local source corrected | page02：PERMANOVA-positive candidate missed by V gate；非 latent truth | publication pending |
| C122 | NEEDS_CORRECTION / P1 | local source corrected | page09：p 是與 specified null 不相容度；V 是 association magnitude | publication pending |
| C123 | NEEDS_CORRECTION / P1 | local source corrected | page09：label-associated centroid separation；須 PERMDISP／design limits | publication pending |
| C132 | NEEDS_CORRECTION / P2 | local source corrected | page02：project-specific combination；novelty／superiority 待 systematic review + benchmark | publication pending |
| C133 | NEEDS_CORRECTION / P1 | local source corrected | README + page15：exit 3 只 guard declared required metrics，不驗證 truth／science | publication pending |
| C145 | CONTRADICTED / P1 | local source corrected | README + page16：version/date-scoped 270 tests／39 suites，future 以 command output 為準 | publication pending |
| C146 | CONTRADICTED / P1 | local source corrected | page15：deploy `fbdf7c7` 全 repo 2,147 tracked `.py` + exact command | publication pending |
| C147 | CONTRADICTED / P1 | local source corrected | page15：`scripts/` 291 tracked `.py`，未排 generated/test/helper | publication pending |
| C156 | CONTRADICTED / P1 | local source corrected | page10：normal 無六個 expected ALT，但有其他 occasional non-reference errors | publication pending |
| C157 | UNVERIFIABLE / P1 | unsupported claim removed | Project Summary：無完整 benchmark 就不報 32-core speedup／<300 ms latency | publication pending |
| C158 | UNVERIFIABLE / P2 | local source corrected | page02：無 dated systematic search／matched benchmark 就不稱 first／unique | publication pending |

## R01 — Inventory 與 registry exact-set 對帳

### 執行命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 - <<'PY'
import csv,json
inv='research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv'
reg='research/20260813_intersubmod_public_surfaces_refresh/scripts/p1_p2_claim_registry.json'
rows=list(csv.DictReader(open(inv,encoding='utf-8'),delimiter='\t'))
sel=[r for r in rows if r['verdict'] in
     {'NEEDS_CORRECTION','CONTRADICTED','UNVERIFIABLE'}
     and r['priority'] in {'P1','P2'}]
data=json.load(open(reg,encoding='utf-8'))
claims=data['claims']
print('inventory='+str(len(sel)))
print('registry='+str(len(claims)))
print('ids='+','.join(c['claim_id'] for c in claims))
print('p1='+str(sum(c['priority']=='P1' for c in claims))+
      ' p2='+str(sum(c['priority']=='P2' for c in claims)))
print('live_pending='+str(sum(c['live_status']=='publication_pending' for c in claims)))
print('claims_with_checks='+str(sum(bool(c['checks']) for c in claims)))
print('claims_with_external_actions='+str(sum(bool(c['external_actions']) for c in claims)))
print('nonempty_required='+str(sum(bool(ch['required_all'])
      for c in claims for ch in c['checks'])))
print('nonempty_forbidden='+str(sum(bool(ch['forbidden_any'])
      for c in claims for ch in c['checks'])))
PY
```

### 實際輸出

```text
inventory=24
registry=24
ids=C047,C048,C089,C090,C110,C111,C114,C115,C116,C117,C118,C119,C120,C121,C122,C123,C132,C133,C145,C146,C147,C156,C157,C158
p1=20 p2=4
live_pending=24
claims_with_checks=24
claims_with_external_actions=24
nonempty_required=26
nonempty_forbidden=26
```

## R02 — JSON／validator 語法與正向 guard

### 執行命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 -c "from pathlib import Path; compile(Path(
  'research/20260813_intersubmod_public_surfaces_refresh/scripts/validate_public_p1_p2_claims.py'
).read_text(), 'validate_public_p1_p2_claims.py', 'exec'); print('validator_syntax=PASS')"
python3 -m json.tool \
  research/20260813_intersubmod_public_surfaces_refresh/scripts/p1_p2_claim_registry.json
python3 research/20260813_intersubmod_public_surfaces_refresh/scripts/validate_public_p1_p2_claims.py
```

### 實際輸出片段

```text
validator_syntax=PASS
JSON_TOOL_EXIT=0
```

```json
{
  "checked_target_rules": 26,
  "dispositions": {
    "local_source_corrected_publication_pending": 21,
    "unsupported_public_claim_removed_publication_pending": 3
  },
  "errors": 0,
  "forbidden_anchors": 50,
  "inventory_problem_p1_p2": 24,
  "live_status": "publication_pending",
  "registry_problem_p1_p2": 24,
  "required_anchors": 80,
  "unique_documents": 12,
  "verdict": "PASS"
}
FINAL_GUARD_EXIT=0
```

## R03 — Fail-closed 負控制

### 執行命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/20260813_intersubmod_public_surfaces_refresh/scripts/validate_public_p1_p2_claims.py \
  --simulate-drop-claim C157
```

### 實際輸出片段

```text
ERRORS:
- claim-set mismatch: missing=['C157'], extra=[]
```

```json
{
  "checked_target_rules": 25,
  "errors": 1,
  "forbidden_anchors": 46,
  "inventory_problem_p1_p2": 24,
  "registry_problem_p1_p2": 23,
  "required_anchors": 76,
  "unique_documents": 11,
  "verdict": "FAIL"
}
NEGATIVE_CONTROL_EXIT=1
```

判定：漏一個 claim 時確實 fail closed；guard 不是永真檢查。

## R04 — P0 regression guard

### 執行命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py
```

### 實際輸出

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
P0_REGRESSION_EXIT=0
```

## R05 — Whitespace 與可追蹤性檢查

### 執行命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
git diff --check -- \
  research/20260813_intersubmod_public_surfaces_refresh/scripts/p1_p2_claim_registry.json \
  docs/explain/15_python-html-layer.standalone.html
git check-ignore -v \
  research/20260813_intersubmod_public_surfaces_refresh/scripts/p1_p2_claim_registry.json \
  research/20260813_intersubmod_public_surfaces_refresh/claim_guard_receipt.md
```

### 實際輸出

```text
DIFF_CHECK_EXIT=0
.gitignore:366:research/**/*.json research/20260813_intersubmod_public_surfaces_refresh/scripts/p1_p2_claim_registry.json
```

注意：registry 目前被既有 broad JSON ignore 規則遮蔽。檔案已存在且 guard 直接由該路徑讀取，
但若要納入 Git，主流程必須另加**精確 negation**或明確 force-add；本子任務依授權沒有修改
`.gitignore`。

## Guard 實際抓到的問題與收斂紀錄

1. 建 registry 前的人工 residual scan 發現 C133 在 page15 圖題仍寫「由構造防止捏造數字」，
   與同頁 desc／caption 的有界說明衝突。Pages stream 將圖題改為「必填指標缺失時
   fail-closed」後，C133 的 forbidden anchor 才通過；guard 沒有為求綠燈而弱化規則。
2. 第一輪 machine guard 有 5 個 errors：C090 的 forbidden regex 誤把「不宣稱 >99%」當
   正向 overclaim；C111、C116、C157 的 required regex 沒容許 Markdown 換行或中英文標點。
   這四類均是 validator contract 的 whitespace／negation calibration，沒有隱藏 source 缺口；
   收窄語境後重跑為 0 errors。
3. C047／C048／C089 刻意不以 current raw BAM bytes 代替 hypothetical tagged-BAM bytes；
   guard 要求「精確 sidecar + tagged total UNVERIFIED + 不報倍率」三者同時成立。

## 能與不能證明

這個 guard 能證明：24 個 inventory claim 無漏項；每個 claim 有明確 local disposition 與外部
動作；指定 source 含必要的 bounded wording且沒有登記的舊錯誤片語；P0 contract 未回歸。

它不能單獨證明：底層生物學結論為真、所有同義 overclaim 都已被自然語言枚舉、遠端公開面
已部署，或未來 commit 的測試／檔案數仍相同。因此 share status 是
**local source ready; public publication pending**，不是 live resolved。

## Provenance footer

- 工作路徑：`/big7_disk/liaoyoyo2001/InterSubMod`
- 日期：2026-08-13（Asia/Taipei）
- 寫入範圍：只新增／更新本 cycle 的 registry 與本收據；未修改 README、Wiki、Pages 或 CCU
- 外部動作：NONE（未 commit、push、merge、PR 或 deploy）
