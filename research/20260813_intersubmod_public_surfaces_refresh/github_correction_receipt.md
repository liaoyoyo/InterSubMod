<!--
建立時間: 2026-08-13
目標: 修正 InterSubMod GitHub README、Quickstart、Project Summary 與 Wiki 的 P1/P2 公開敘述，並保存逐 claim disposition 與驗證證據
處理範圍: claim_inventory.tsv 中 verdict=NEEDS_CORRECTION|CONTRADICTED|UNVERIFIABLE 且 priority=P1|P2 的 24 個 claim families；禁止修改 docs/explain Pages 與 CCU
task_type: B_comprehensive_validation
服務目標: G2/G3/G4/G5
狀態: LOCAL_SOURCE_CORRECTED_EXTERNAL_PUBLISH_REQUIRED
關聯檔案:
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/source_scope.tsv
  - InterSubMod/docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.md
  - InterSubMod/research/20260813_public_docs_p0_correction/claim_guard_receipt.md
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json
-->

# InterSubMod GitHub public surfaces P1/P2 correction receipt

## 結論

鎖定母體為 **24 個 P1/P2 problem claim families**。依 2026-08-12 occurrence inventory：

- 14 個有 README／Summary／Wiki occurrence，已逐 occurrence 做本機 source correction。
- 10 個 inventory occurrence 只在 `docs/explain` Pages；本子任務未修改 Pages，列為
  `DELEGATED_TO_PAGES`。
- 另在 current derived Wiki 發現 C117 的新／未登記 recurrence，已將該 Wiki wording 收斂；
  C117 的 Pages occurrences 仍交由 Pages stream。
- 因此 24/24 都有 disposition：**15 個至少有一個 GitHub source correction**（含 C117
  Wiki drift），**9 個純 Pages delegation**；這不代表 15 個 claim 的 live publication 已完成。
- 既有 P0 guard 最終 PASS：34/34 P0 registry 完整，77 個 target rules、140 個 required
  anchors、79 個 forbidden anchors，errors=0。

`LOCAL_SOURCE_CORRECTED` 只表示共享 working tree 的指定來源已修正。GitHub default branch、
Wiki remote 與 Pages deploy 都是獨立外部狀態；本輪沒有 commit、push、merge、PR 或 deploy。

## Scope 與假設

- 只修改 README、Project Summary 與 `docs/wiki/*.md`；`QUICKSTART.md` 已在 P0 cycle 收斂且
  本輪沒有新增差異。
- 不修改 `InterSubMod/docs/explain/*.standalone.html`，即使 claim 仍存在於 Pages。
- 不修改 CCU lab-tutorial。
- 不把 current raw BAM bytes 代替 hypothetical tagged-BAM bytes。
- 不把 shared dirty tree 的整體 `git diff` 冒充本子任務獨有 diff；下列 changed-files 清單來自
  本子任務實際 `apply_patch` targets。

## 輸入路徑

1. `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv`
2. `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/source_scope.tsv`
3. `InterSubMod/docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.md`
4. `InterSubMod/research/20260813_public_docs_p0_correction/claim_guard_receipt.md`
5. `InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json`
6. 本機 source：`InterSubMod/{README.md,README.zh-TW.md,QUICKSTART.md,README_PROJECT_SUMMARY.md}`
   與 `InterSubMod/docs/wiki/*.md`

## Changed files

本子任務含 supplemental drift correction，共以 `apply_patch` 修改 11 個檔案：

1. `InterSubMod/README.md`
2. `InterSubMod/README.zh-TW.md`
3. `InterSubMod/README_PROJECT_SUMMARY.md`
4. `InterSubMod/docs/wiki/Home.md`
5. `InterSubMod/docs/wiki/System-Overview.md`
6. `InterSubMod/docs/wiki/Upstream-and-Data.md`
7. `InterSubMod/docs/wiki/Analysis-and-Presentation.md`
8. `InterSubMod/docs/wiki/How-to-Run.md`
9. `InterSubMod/docs/wiki/README.md`
10. `InterSubMod/docs/wiki/InterSubMod-Engine.md`
11. `InterSubMod/docs/wiki/LongLineage-Engine.md`

新增本收據：
`InterSubMod/research/20260813_intersubmod_public_surfaces_refresh/github_correction_receipt.md`。

未修改：`InterSubMod/QUICKSTART.md`、`InterSubMod/docs/explain/`、CCU repo。

## Claim-by-claim disposition

| Claim | GitHub source 修正 | Pages 狀態 | Disposition |
|---|---|---|---|
| C047 | README EN/ZH、System、Upstream：1.67 TiB 標 UNVERIFIED；列缺少的 7-path/exact-byte/hash/compression receipt | P014 未動 | LOCAL_SOURCE_CORRECTED / DELEGATED_TO_PAGES / PUBLISH_REQUIRED |
| C048 | README EN/ZH、Upstream：移除 287×；沒有兩端 exact-byte receipt 就不報倍率 | P014 未動 | LOCAL_SOURCE_CORRECTED / DELEGATED_TO_PAGES / PUBLISH_REQUIRED |
| C089 | README EN/ZH、System、Upstream：只確認 sidecar exact bytes，不稱整個 storage snapshot fully verified | P014 未動 | LOCAL_SOURCE_CORRECTED / DELEGATED_TO_PAGES / PUBLISH_REQUIRED |
| C090 | README EN/ZH、Upstream：只說九欄 contract 不保留 SEQ/QUAL；不稱 >99% payload 或 0% utility | P014 未動 | LOCAL_SOURCE_CORRECTED / DELEGATED_TO_PAGES / PUBLISH_REQUIRED |
| C110 | README EN/ZH、Wiki Home：限於 per-locus marginal VAF only、無 linkage／extra assumptions 的不可識別性 | P011 未動 | LOCAL_SOURCE_CORRECTED / DELEGATED_TO_PAGES / PUBLISH_REQUIRED |
| C111 | README EN/ZH、Wiki Home：限於 current single-bulk measurement set，無 orthogonal data／extra assumptions | P011 未動 | LOCAL_SOURCE_CORRECTED / DELEGATED_TO_PAGES / PUBLISH_REQUIRED |
| C114 | README EN/ZH、System：明指 `cohort_receipt.json` 與 `summary/all7_summary.json`，不是 authority manifest 頂層 | P000/P011 未動 | LOCAL_SOURCE_CORRECTED / DELEGATED_TO_PAGES / PUBLISH_REQUIRED |
| C115 | README EN/ZH、Wiki README：`docs/explain` 是 editorial upstream；Wiki 是人工同步 derivative，可漂移 | inventory 無 Pages claim occurrence | LOCAL_SOURCE_CORRECTED / PUBLISH_REQUIRED |
| C116 | README EN/ZH：deploy `fbdf7c7` 的 17 頁有 37 inline SVG elements；不等同 37 張語意獨立圖 | inventory 無 Pages correction target | LOCAL_SOURCE_CORRECTED / PUBLISH_REQUIRED |
| C117 | Current System Wiki 的 `DEAD／勿再開` 漂移改成「tested formulation negative；materially new hypothesis + pre-decision audit 才重開」 | P000–P003 未動 | WIKI_DRIFT_CORRECTED / DELEGATED_TO_PAGES / PUBLISH_REQUIRED |
| C118 | inventory 無 A/W occurrence | P000/P002 未動 | DELEGATED_TO_PAGES |
| C119 | inventory 無 A/W occurrence | P001 未動 | DELEGATED_TO_PAGES |
| C120 | inventory 無 A/W occurrence | P001/P002 未動 | DELEGATED_TO_PAGES |
| C121 | inventory 無 A/W occurrence | P001/P002 未動 | DELEGATED_TO_PAGES |
| C122 | inventory 無 A/W occurrence | P009 未動 | DELEGATED_TO_PAGES |
| C123 | inventory 無 A/W occurrence | P002/P009 未動 | DELEGATED_TO_PAGES |
| C132 | inventory 無 A/W occurrence | P002/P005 未動 | DELEGATED_TO_PAGES |
| C133 | README EN/ZH、Analysis：exit 3 只防已宣告 required metrics 靜默缺失；不保證 truth、denominator、source 或 optional fields | P015 未動 | LOCAL_SOURCE_CORRECTED / DELEGATED_TO_PAGES / PUBLISH_REQUIRED |
| C145 | README 已由 P0 cycle 釘 `73afaeac-dirty`／date／command；Home/System/How-to 移除未限定 265/38，How-to 保留 version-scoped 270/39 snapshot 並要求 current output/exit | P000/P016 未動 | LOCAL_SOURCE_CORRECTED / DELEGATED_TO_PAGES / PUBLISH_REQUIRED |
| C146 | Analysis：改為 deploy `fbdf7c7` 的 2,147 tracked `.py`，附 exact command | P015 未動 | LOCAL_SOURCE_CORRECTED / DELEGATED_TO_PAGES / PUBLISH_REQUIRED |
| C147 | Analysis：改為 deploy `fbdf7c7` 的 `scripts/` 291 tracked `.py`；明示不排 tests/generated-named sources，附 exact command | P015 未動 | LOCAL_SOURCE_CORRECTED / DELEGATED_TO_PAGES / PUBLISH_REQUIRED |
| C156 | inventory 無 A/W occurrence | P010 未動 | DELEGATED_TO_PAGES |
| C157 | Project Summary 移除數值 speedup／latency；只保留 OpenMP 設計與 benchmark 缺口 | inventory 無 Pages occurrence | LOCAL_SOURCE_CORRECTED / PUBLISH_REQUIRED |
| C158 | inventory 無 A/W occurrence | P002/P005 未動 | DELEGATED_TO_PAGES |

## 關鍵數值與算法邊界

- Current 7 sidecars：`6,256,168,164 bytes = 5.826510641724 GiB`，可由既有 command receipt 重播。
- `1.67 TiB`：沒有 committed tagged-BAM per-file registry，維持 `UNVERIFIED`。
- `287×`：已從 GitHub source 的正向 claim 移除；顯示值本身算術為 293.34×，但沒有 tagged
  exact receipt，故本輪不以 293.34× 取代。
- Audited Pages deploy `fbdf7c7`：全 repo 2,147 個 tracked `.py`；`scripts/` 291 個；17 個
  standalone HTML 合計 37 個 inline `<svg>` elements。這些都是 version-scoped corpus counts。
- Test count：270/39 只屬 tracked core `73afaeac`、2026-08-12 的 snapshot；future commit 的
  authority 是實際 command output 與 exit code，不是 hard-coded count。

## 執行命令與實際輸出

### R01 — 母體與 surface partition

完整命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 - <<'PY'
import csv
rows=list(csv.DictReader(open(
    'research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv',
    encoding='utf-8'), delimiter='\t'))
sel=[r for r in rows if r['verdict'] in
     {'NEEDS_CORRECTION','CONTRADICTED','UNVERIFIABLE'}
     and r['priority'] in {'P1','P2'}]
gh=[r for r in sel if any(x.startswith(('A','W'))
    for x in r['occurrences'].split('|'))]
print(f'problem_p1_p2={len(sel)} locked_github_occurrence_claims={len(gh)} '
      f'locked_pages_only_claims={len(sel)-len(gh)}')
print('github_ids='+','.join(r['claim_id'] for r in gh))
print('pages_only_ids='+','.join(r['claim_id'] for r in sel if r not in gh))
PY
```

實際輸出：

```text
problem_p1_p2=24 locked_github_occurrence_claims=14 locked_pages_only_claims=10
github_ids=C047,C048,C089,C090,C110,C111,C114,C115,C116,C133,C145,C146,C147,C157
pages_only_ids=C117,C118,C119,C120,C121,C122,C123,C132,C156,C158
```

### R02 — P0 regression guard

完整命令：

```bash
python3 research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py
```

實際輸出：

```text
inventory_p0=34; registry_p0=34; checked_target_rules=77
required_anchors=140; forbidden_anchors=79; unique_documents=27
errors=0; verdict=PASS
```

### R03 — Markdown relative links

完整命令：

```bash
python3 scripts/analysis/check_md_links.py \
  README.md README.zh-TW.md QUICKSTART.md README_PROJECT_SUMMARY.md docs/wiki/*.md
```

實際輸出：無；exit 0。

### R04 — Bash fenced-block syntax（只 parse、不執行範例）

完整命令：

```bash
python3 - <<'PY'
from pathlib import Path
import re, subprocess, sys
files=[Path('README.md'),Path('README.zh-TW.md'),Path('QUICKSTART.md'),
       Path('README_PROJECT_SUMMARY.md'),*sorted(Path('docs/wiki').glob('*.md'))]
errors=[]; blocks=0
for p in files:
    for i,m in enumerate(re.finditer(
            r'```(?:bash|sh)\n(.*?)\n```',p.read_text(encoding='utf-8'),re.S),1):
        blocks+=1
        proc=subprocess.run(['bash','-n','-c',m.group(1)],capture_output=True,text=True)
        if proc.returncode: errors.append(f'{p}:block-{i}:{proc.stderr.strip()}')
print(f'markdown_files={len(files)} bash_or_sh_blocks={blocks} syntax_errors={len(errors)}')
for e in errors: print('ERROR',e)
sys.exit(bool(errors))
PY
```

實際輸出：

```text
markdown_files=13 bash_or_sh_blocks=19 syntax_errors=0
```

### R05 — EN/ZH parity contracts

完整命令：

```bash
python3 - <<'PY'
from pathlib import Path
import sys
checks={
 'README.md':['not identifiable from those marginals alone',
   'current single-bulk measurement set','cohort_receipt.json','editorial upstream',
   '37 inline `<svg>` elements','6,256,168,164 bytes','UNVERIFIED',
   'omitted optional fields'],
 'README.zh-TW.md':['無法只從這些邊際值識別唯一聯合結構',
   'single-bulk measurement set','cohort_receipt.json','editorial upstream',
   '37 個 inline `<svg>` elements','6,256,168,164 bytes','UNVERIFIED',
   '未宣告的 optional fields'],
}
errors=[]
for p,needles in checks.items():
    s=Path(p).read_text(encoding='utf-8')
    for n in needles:
        if n not in s: errors.append(f'{p}: missing {n}')
print(f'en_zh_contracts={sum(map(len,checks.values()))} errors={len(errors)}')
for e in errors: print('ERROR',e)
sys.exit(bool(errors))
PY
```

實際輸出：`en_zh_contracts=16 errors=0`。

### R06 — Version-scoped inventory counts

完整命令：

```bash
git ls-tree -r --name-only fbdf7c7 | rg -c '\.py$'
git ls-tree -r --name-only fbdf7c7 -- scripts | rg -c '\.py$'
git grep -h -o '<svg[ >]' fbdf7c7 -- 'docs/explain/*.standalone.html' | wc -l
```

實際輸出：

```text
2147
291
37
```

### R07 — 舊 claim residual 與 whitespace

完整命令：

```bash
rg -n -i \
 'provably non-identifiable|已被證明無解|29 hand-drawn|29 張|287×|287x|265 tests|2,144|scripts/.*268|structurally impossible|構造上不可能|DEAD|勿再開|32 核.*線性加速|< ?300 ?ms' \
 README.md README.zh-TW.md QUICKSTART.md README_PROJECT_SUMMARY.md docs/wiki --glob '*.md'
git diff --check -- \
  README.md README.zh-TW.md QUICKSTART.md README_PROJECT_SUMMARY.md docs/wiki
```

實際輸出：residual `rg` 無命中、exit 1（預期）；`git diff --check` 無輸出、exit 0。

## Supplemental drift correction：schema、LongLineage scope 與 runtime

2026-08-13 的追加全文掃描發現三組相鄰敘述 drift；以下均已在本機 Wiki source 修正，
未修改 `docs/explain` 或 CCU。

| Drift | 修正 | 證據／限制 | Disposition |
|---|---|---|---|
| `InterSubMod-Engine.md` 同時寫 source 193 欄與 audited 199 欄 | 統一為 historical files 59/114/117/157/180；tracked core `73afaeac`、2026-08-12 runtime/source header 199 欄；col187 `VerificationSchemaVersion=2`、col196 `RegionStratificationSchemaVersion=1` | 199 是版本快照；兩欄是 component-level schema，仍無 whole-file layout version | LOCAL_SOURCE_CORRECTED / PUBLISH_REQUIRED |
| `LongLineage-Engine.md` 把 0 topology 寫成一般 current real-data 結論 | headline、intro、section、table 與 artifact wording 全部限定 frozen HCC1395-only dataset-gate、79,687 loci、0 topology units | 符合 C031 minimum rewrite；不可外推每個 real-data run、其他樣本或所有版本 | LOCAL_SOURCE_CORRECTED / PUBLISH_REQUIRED |
| System／How-to／ISM Engine 報 2.9 秒但缺完整 provenance | 移除秒數，改為 historical internal single-locus smoke exit 0；若要報 runtime，須補 hardware、input locus、commit、date、repetitions | 符合 C088/C157 的 version/provenance ceiling | LOCAL_SOURCE_CORRECTED / PUBLISH_REQUIRED |
| ISM binary `STALE` 以 mtime 作永久敘述 | 改為每次記錄 source commit、dirty diff、binary SHA-256、compiler、build command/date；不以 mtime 當公開 stale/current authority | 工作目錄瞬時狀態不能泛化到 public clone | LOCAL_SOURCE_CORRECTED / PUBLISH_REQUIRED |

Supplemental 驗證命令：

```bash
rg -n -i '193|199|59|114|117|157|180|VerificationSchemaVersion|RegionStratificationSchemaVersion' \
  docs/wiki/InterSubMod-Engine.md
rg -n -i '2\.9.{0,12}秒|約 3 秒|binary 是 \*\*STALE|執行檔是 \*\*STALE|5 個原始碼' \
  docs/wiki/InterSubMod-Engine.md docs/wiki/System-Overview.md docs/wiki/How-to-Run.md
rg -n -i 'HCC1395-only|79,687|0 topology units|每個 real-data run|所有版本都為 0' \
  docs/wiki/LongLineage-Engine.md docs/wiki/System-Overview.md
python3 research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py
python3 scripts/analysis/check_md_links.py docs/wiki/*.md
git diff --check -- docs/wiki \
  research/20260813_intersubmod_public_surfaces_refresh/github_correction_receipt.md
```

實際結果：

- ISM 欄數只剩兩組一致、version-scoped 說明；沒有 193-column current claim。
- `2.9 秒`、`約 3 秒`、mtime-based `STALE`／`5 個原始碼` 無正向殘留；
  `InterSubMod-Engine.md` 的資料範例仍合法包含 GlobalP `2.990000e-01`，它不是 runtime。
- LongLineage 的 0-unit wording 都與 frozen HCC1395-only scope 或明確禁止外推同段出現。
- P0 guard、relative links 與 `git diff --check` 全部 exit 0；完整 guard counters 見最終重跑。

## 尚未完成

1. 10 個含 Pages occurrence 的 claim 都仍需 Pages stream 按原 occurrence 修正；C117 雖已修
   current Wiki drift，Pages 仍未完成。
2. GitHub default branch、Wiki repo 與 Pages 未發布；local source 不等於 live resolved。
3. Tagged-BAM exact receipt、field-level BAM census、performance benchmark 與 dated systematic
   novelty review 仍缺；在補件前應維持 UNVERIFIED／不報數字，而不是推估替代。
4. `docs/images/howto-six-steps.png` 可能仍含 historical test-count label；本輪未獲圖片修改範圍，
   Wiki caption 已明示圖中固定數字只屬 historical illustration。

## Provenance footer

- 工作路徑：`/big7_disk/liaoyoyo2001/InterSubMod`
- 日期：2026-08-13（Asia/Taipei）
- 寫入方式：只用 `apply_patch`
- 外部動作：NONE（未 commit、push、merge、PR 或 deploy）
