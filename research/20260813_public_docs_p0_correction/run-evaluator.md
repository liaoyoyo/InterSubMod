<!--
建立時間: 2026-08-13
目標: 判定公開文件 P0 修正循環是否適用 /run-evaluator P5 gate，並保存不執行 mechanical score 的證據
處理範圍: Task B comprehensive validation continuation；服務 G2/G3/G4/G5
cycle_id: 20260813_public_docs_p0_correction
status: NOT_APPLICABLE
run_evaluator_status: NOT_RUN
關聯檔案:
  - InterSubMod/.claude/skills/run-evaluator/SKILL.md
  - InterSubMod/research/20260813_public_docs_p0_correction/00_INDEX.md
  - InterSubMod/research/20260813_public_docs_p0_correction/claim_guard_receipt.md
  - InterSubMod/research/20260813_public_docs_p0_correction/remote_state_receipt.md
  - InterSubMod/research/20260813_public_docs_p0_correction/ccu_correction_receipt.md
-->

# `/run-evaluator` gate 判定：NOT_APPLICABLE／NOT_RUN

> **結論：本輪不可套用 `/run-evaluator` 的 ⭐1–5 與 6-component mechanical score。**
> 這是公開文件 P0 correction cycle，不是已完成 P4 的科學跨樣本 generalization cycle；必要的
> `state.json`、`plan.json`、`pilot.json`、`generalize.json` 全部不存在。依 skill 的
> **DO NOT USE WHEN**，若以 mock／缺值代算會產生虛假低風險。因此本文件只保存 eligibility
> gate 與實際 correction evidence，不建立 `evaluation.json`、不更新 state mirror、不評星。

## 1. Gate 判定

| 項目 | 必要條件 | 實際觀察 | 判定 |
|---|---|---|---|
| Cycle 類型 | 已完成 P4 的研究／實驗 cycle | Task B 文件與教學敘述修正；本輪不重跑 7 樣本 BAM pipeline | 不適用 |
| `state.json` | 必須存在且 phase 可進 P5 | 缺少 | 不可執行 |
| `plan.json` | 必須含 hypothesis、tier、preconditions | 缺少 | 不可執行 |
| `pilot.json` | 必須含 subgroup metric results | 缺少 | 不可執行 |
| `generalize.json` | 必須含跨樣本 consistency 與 sample metrics | 缺少 | 不可執行 |
| `precheck.json` | 建議；缺少時單一 component 可給 0.5 | 缺少；但其他四個必要 artifact 也缺，不可只靠中性值代算 | 不可執行 |
| Tier upgrade | evaluator 是 ⭐4／⭐5 升級 gate | 本輪沒有科學 tier upgrade request | 不評級 |

**Final gate verdict：`NOT_APPLICABLE / NOT_RUN`。** 這不是 `approve_tier`、
`downgrade_tier` 或 `pending_review` 三者之一；因為尚未進入可計分母體。任何既有文件中的
`L2 ⭐⭐⭐⭐` 標示均**未經本 evaluator 核准**，不可把本文件解讀成 tier endorsement。

## 2. 六個 risk components：全部 N/A，不填 0 或 0.5

| Component | Skill 所需資料 | 本輪狀況 | 值 |
|---|---|---|---|
| `multi_sample_consistency` | `n_samples_passed / n_samples_total` | 文件 surface 不是 biological samples；沒有 `generalize.consistency` | N/A |
| `effect_size_stability` | 各樣本 `metric_value` 的 min/max | 本輪沒有 effect-size metric，只有 claim/QA assertions | N/A |
| `precondition_freshness` | `precheck.verdict` | `precheck.json` 缺；不能用單一 0.5 掩蓋整個 P4 artifact 缺失 | N/A |
| `subgroup_homogeneity` | pilot subgroup metric 的 stddev/mean | 沒有 pilot subgroup | N/A |
| `pitfall_coverage_score` | `plan.hypothesis` keyword sweep + 13 條 rule 結果 | `plan.json` 缺；不能聲稱 13/13 pitfalls 已檢查 | N/A |
| `tier_support_alignment` | ledger stability 與 `plan.tier_used` | 沒有 `plan.tier_used`，且沒有本輪 tier 升級 | N/A |

- `retraction_risk`：**NOT_COMPUTED**。
- `base_verdict`／`final_verdict`：**NOT_COMPUTED**。
- `tier_recommendation`：**NOT_RATED**（不是 1–5 中任一級）。
- `precedent_similarity`：Path A 原本即為 `null`；本輪不建立 evaluation payload。

### Manual pitfall observations（不是自動 sweep 分數）

- P-13 working-tree dirty 是實際風險，但本輪沒有產出新的 binary／生物數據結論；文件 correction
  以 target-scoped diff、明示 build identity 與 patch checksum 管理。它不能被誤寫成
  `pitfall_coverage_score=1.0`。
- P-07 single-track tier gate 沒有正式觸發，因為本輪沒有 `plan.tier_used` 或 ⭐4／⭐5 升級申請；
  既有證據星等若要保留，仍應由 scientific-rigor／人工 reviewer 另行核定。
- 其餘 P-01～P-12 多數針對統計、caller、KDE、AF、跨樣本與 feature semantics；沒有
  `plan.hypothesis` 時不能用「無 keyword 命中」冒充已檢查或不相關。

## 3. 實際輸入、命令、輸出與觀察

### R01 — P5 artifact eligibility

輸入路徑：
`InterSubMod/state/cycles/20260813_public_docs_p0_correction/`。

執行命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
for n in state plan precheck pilot generalize evaluation; do
  p="state/cycles/20260813_public_docs_p0_correction/$n.json"
  if test -f "$p"; then echo "$n=PRESENT"; else echo "$n=MISSING"; fi
done
```

實際輸出：

```text
state=MISSING
plan=MISSING
precheck=MISSING
pilot=MISSING
generalize=MISSING
evaluation=MISSING
```

觀察：skill 的最小輸入不成立；停止 mechanical evaluation 是正確 fail-closed 行為。

### R02 — 本地 InterSubMod P0 claim guard

輸入路徑：

- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv`
- `InterSubMod/research/20260813_public_docs_p0_correction/scripts/p0_claim_registry.json`
- registry 登記的 27 份 README／Quickstart／Wiki／Pages／page-07 generator source

執行命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py
python3 research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py \
  --simulate-drop-claim C155
```

實際輸出片段：

```text
inventory_p0=34; registry_p0=34; checked_target_rules=77
required_anchors=140; forbidden_anchors=79; unique_documents=27
errors=0; verdict=PASS
GUARD_EXIT=0

ERRORS:
- registry missing inventory P0 IDs: C155
verdict=FAIL
NEGATIVE_GUARD_EXIT=1
```

觀察：34/34 P0 有 disposition；33 個本地 claim 有 source rules，C108 為純外部 About action。
正向 guard 通過，且漏掉 C155 的 mutation probe 確實 fail closed。這證明本地登記範圍的 source
contract，不證明 live publication 或底層生物結論。

### R03 — CCU patch identity 與 apply-check

輸入路徑：

- `/tmp/lab_tutorial_p0.GIEknc/repo`，index baseline 為 CCU `9eb1618d...`
- `InterSubMod/research/20260813_public_docs_p0_correction/ccu_source_correction_9eb1618.patch`

執行命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
sha256sum research/20260813_public_docs_p0_correction/ccu_source_correction_9eb1618.patch
git -C /tmp/lab_tutorial_p0.GIEknc/repo rev-parse HEAD
git -C /tmp/lab_tutorial_p0.GIEknc/repo apply --check --cached \
  /big7_disk/liaoyoyo2001/InterSubMod/research/20260813_public_docs_p0_correction/ccu_source_correction_9eb1618.patch
```

實際輸出：

```text
aee944790ad8b60634818fca8888b0f980b3c1e5522ef64c00015a0a21d6c3e6  .../ccu_source_correction_9eb1618.patch
9eb1618d359e602d9c528675952b20d051fb2346
CCU_PATCH_CACHED_APPLY_CHECK_EXIT=0
```

觀察：source-only patch 可對 baseline index 乾淨套用；CCU receipt 另記錄本地 build/verify 與
13/13 finding assertions PASS。精確狀態是 `PATCH_VALIDATED_ON_PINNED_CLONE / NOT_APPLIED / NOT_DEPLOYED`，不是 live fixed。

### R04 — Build／tests 證據邊界

輸入路徑：
`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/command_receipts.md`
的 R07／R08。

既有執行命令：

```bash
cmake -S /big7_disk/liaoyoyo2001/InterSubMod \
  -B /tmp/intersubmod-public-audit.XNde4y -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/intersubmod-public-audit.XNde4y -j$(nproc)
/tmp/intersubmod-public-audit.XNde4y/bin/run_tests
ctest --test-dir /tmp/intersubmod-public-audit.XNde4y --output-on-failure
```

收據中的實際輸出：

```text
build identity: 73afaeac-dirty; GCC 11.4.0; htslib 1.18; build exit 0
Running 270 tests from 39 test suites.
PASSED: 270 tests.
CTest: 100% tests passed, 0 tests failed out of 270.
```

觀察：這是 2026-08-12 全面稽核對 current public-doc test-count 的 runtime 證據，不是本文件
correction cycle 的 `pilot.json`／`generalize.json`，也沒有在本 evaluator 重跑。當前 HEAD 為
`95d420f6`，且共享 dirty tree 中另有非本輪的 `tests/` 修改；因此不可把舊 build receipt
升格成「目前整個 dirty worktree 已重新驗證」。

### R05 — Live remote state 重查

輸入：InterSubMod main/develop、Wiki master、CCU main 與四個 live content endpoints。

執行命令：

```bash
git ls-remote https://github.com/liaoyoyo/InterSubMod.git \
  refs/heads/main refs/heads/develop
git ls-remote https://github.com/liaoyoyo/InterSubMod.wiki.git refs/heads/master
git ls-remote https://github.com/ccu-bioinformatics-lab/lab-tutorial.git refs/heads/main
curl -L -sS -o /dev/null -w '%{http_code}' <URL>
curl -L -fsS <URL> | sha256sum
```

實際輸出摘要：

```text
InterSubMod main    635437a65a33f8ba698acf85b22ebb069455c6cc
InterSubMod develop ddd8909a838318d8a77969313e9561c8ff9d01c2
Wiki master         6cfc990f1fbe0b0aa76b8adbe56bbe645ef7ee7b
CCU main            9eb1618d359e602d9c528675952b20d051fb2346
main README         HTTP 200, sha256 8187a83223d538d0ef25b12722535397f67c0dade394298785ee626e52e54efa
main README.zh-TW    HTTP 404
Wiki Home           HTTP 200, sha256 a17d9ee25119062c05ae0e101aaf7d5b201ce1436b5b4f9378ece1134a64175b
Pages index         HTTP 200, sha256 deecd51ed7446185551de47c84903f8cfb8477de2a5317b36c1202375b9941b6
CCU index           HTTP 200, sha256 bc0821bcd5ea6df2326d1b7bb6560d518f8e31599b4ddf34e0d1f529a4a0466f
```

觀察：與 `InterSubMod/research/20260813_public_docs_p0_correction/remote_state_receipt.md`
相同；本機修正尚未發布。`LOCAL_SOURCE_CORRECTED`／`PATCH_VALIDATED_ON_PINNED_CLONE` 不得改寫為
`LIVE_RESOLVED`。

### R06 — Final report／browser artifact readiness

輸入與命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
test -f docs/reports/validated/2026/08/20260813_公開文件P0修正與驗證_01.md
test -f docs/reports/validated/2026/08/20260813_公開文件P0修正與驗證_01.standalone.html
test -f research/20260813_public_docs_p0_correction/report_qa/qa_receipt.json
```

評估快照的實際結果：三者皆 `MISSING`。

觀察：在本 gate 執行時間點，總報告、offline HTML 與 desktop/mobile/print browser QA 尚未形成；
即使後續補齊，這只能關閉 document-delivery gate，不會把本 cycle 轉成可套用 6-component
scientific evaluator 的 P4 cycle。

**Post-evaluator follow-up（2026-08-13）**：上述三個 artifact 後續均已形成；終局 browser QA 為
desktop/mobile 0 overflow、0 console/page/external-request error、4/4 SVG accessibility PASS、
print details 0 hidden、exit 0。Fresh-reader 修訂後另加入 33-entry correction manifest 與
final-report SHA-256 sidecar，將未 commit target／guard／patch／MD／HTML bytes 固定；完整結果見
`InterSubMod/research/20260813_public_docs_p0_correction/html_qa_receipt.md`。此補件關閉文件交付 gate，
但不改變本 evaluator 的 `NOT_APPLICABLE / NOT_RUN` 判定。

## 4. Correction-cycle completion risks（非 retraction-risk score）

| Risk | 當前證據 | 影響 | 關閉條件 |
|---|---|---|---|
| Live publication 未做 | refs/hash 仍是 audited old state；繁中 README 仍 404 | 高 | owner 授權 merge/publish/deploy，逐面 live re-fetch |
| GitHub About 未改 | live About 仍含 `subclone resolution` | 高 | repo owner 改文案並保存 live receipt |
| P1/P2 未處理 | InterSubMod 58 個問題中本輪只處理 34 P0；另 24 個仍在 queue | 中 | 開獨立 P1/P2 cycle，不混入 P0 numerator |
| CCU 非本輪 finding | OLD-P1 11 項與 non-REGRESSED delta 5 項仍待辦 | 中 | finding-by-finding evidence/correction cycle |
| P007 generator drift | **已關閉**：generator 已修、重建 P007，page-07 3 SVG／0 FAIL；C125/C127/C128 已加 generator guard | 低 | 後續 CI 固定重跑 generator + claim guard，防止復發 |
| SVG accessibility debt | 02–09 多數既有 SVG 缺 `<title>/<desc>`，目前是 WARN | 中 | 補 title/desc/ARIA 後重跑全 17 頁 QA |
| Public fixture 缺失 | 無可授權 tiny BAM/FASTA/VCF，Quickstart 只能標 internal/self-provided | 高 | fixture + checksum + license + portable wrapper + CI smoke |
| Final report QA | **已關閉**：MD／HTML／browser receipt 已形成；desktop/mobile/print/SVG gate PASS | 低 | 後續內容異動時重跑相同 QA |

## 5. Decision 與下一個合法 gate

1. 本輪 `/run-evaluator` 保持 **NOT_APPLICABLE／NOT_RUN**；不建立或 mock
   `InterSubMod/state/cycles/20260813_public_docs_p0_correction/evaluation.json`。
2. 文件 correction 的本地 acceptance 應由 claim guard、source receipts、patch apply-check、HTML/browser
   QA 與獨立 reader test 判定；這些是 document QA，不是跨樣本 retraction-risk。
3. 對外完成必須另過 publication gate：commit/review → push/merge → Wiki/Pages/CCU deploy → live
   ref/hash/anchor re-fetch。未完成前，唯一允許狀態是 `LOCAL_SOURCE_CORRECTED` 或
   `PATCH_VALIDATED_ON_PINNED_CLONE / NOT_APPLIED / NOT_DEPLOYED`。
4. 若未來真的要使用 `/run-evaluator`，必須是具有實際 hypothesis、pilot metrics、跨樣本
   generalization 與 tier proposal 的科學 cycle；不得為了取得星級替本輪文件修正製造 mock artifacts。

## Provenance footer

- 工作路徑：`/big7_disk/liaoyoyo2001/InterSubMod`
- 評估時 HEAD：`95d420f6`
- 評估日期：2026-08-13（Asia/Taipei）
- 寫入輸出：`InterSubMod/research/20260813_public_docs_p0_correction/run-evaluator.md`
- 外部動作：NONE（未 commit、push、merge、PR、deploy 或修改 GitHub About）
