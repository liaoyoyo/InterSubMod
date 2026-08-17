<!--
建立時間: 2026-08-13
目標: 修正 docs/wiki 全部 P0 公開敘述並留下可重驗收據
處理範圍: docs/wiki/*.md（8 個公開 Wiki source pages）
關聯檔案: research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
-->

# Wiki P0 correction receipt

## 結論

Task type **B — Comprehensive validation continuation**；服務 G2／G3／G4／G5。`claim_inventory.tsv` 中所有 `priority=P0` 且 occurrence 含 `W001`–`W008` 的 10 個 claim 已在 `InterSubMod/docs/wiki/*.md` source 修正：**10/10 LOCAL_SOURCE_CORRECTED**。本輪沒有 commit、push 或發布 GitHub Wiki；live Wiki 仍需另行同步，因此不能標成 live resolved。

## Step → Verify

1. 盤點 Wiki P0 claim → 驗證：機械查得 10 個 ID（C043、C059、C107、C109、C112、C113、C142、C148、C149、C151）。
2. 修正 8 個 Wiki source pages → 驗證：`git diff --numstat -- docs/wiki` 為 121 insertions／81 deletions，且未改 `docs/wiki/README.md`。
3. 做 residual、結構與 shell snippet 檢查 → 驗證：舊錯誤句型 0 命中；9 檔／115 headings／86 links／0 errors；13 個 bash fences／0 syntax errors；`git diff --check` exit 0。
4. 保留外部發布邊界 → 驗證：本輪未執行 `git commit`、`git push` 或 GitHub Wiki clone/push。

## 輸入、命令、輸出

### 輸入路徑

- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/source_scope.tsv`
- `InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json`
- `InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv`
- `InterSubMod/docs/wiki/*.md`

註：任務文字所列 audit research root 內沒有 `authority_manifest.json`／`denominator_registry.tsv` 副本；`claim_inventory.tsv` 的 evidence 欄明確指向上述 handoff canonical files，本輪依該實際路徑讀取。

### 初始潔淨度

```bash
git status --short -- docs/wiki \
  research/20260813_public_docs_p0_correction/wiki_correction_receipt.md
git diff -- docs/wiki
git diff --cached -- docs/wiki
```

實際輸出：空（target 在編輯前無 tracked／staged／untracked 變更）。

### P0 inventory

```bash
awk -F '\t' '
  NR==1 {for(i=1;i<=NF;i++) h[$i]=i; next}
  $h["priority"]=="P0" && $h["occurrences"] ~ /W00[1-8]/ {
    print $h["claim_id"], $h["category"], $h["occurrences"]
  }' \
  research/20260812_intersubmod_github_public_docs_full_validation/claim_inventory.tsv
```

實際輸出 ID：`C043 C059 C107 C109 C112 C113 C142 C148 C149 C151`（10 個）。

### 輸出路徑

- `InterSubMod/docs/wiki/Analysis-and-Presentation.md`
- `InterSubMod/docs/wiki/Home.md`
- `InterSubMod/docs/wiki/How-to-Run.md`
- `InterSubMod/docs/wiki/InterSubMod-Engine.md`
- `InterSubMod/docs/wiki/LongLineage-Engine.md`
- `InterSubMod/docs/wiki/System-Overview.md`
- `InterSubMod/docs/wiki/Upstream-and-Data.md`
- `InterSubMod/docs/wiki/_Sidebar.md`
- 本收據：`InterSubMod/research/20260813_public_docs_p0_correction/wiki_correction_receipt.md`

## Claim-by-claim before／after／evidence／disposition

| ID | Before（問題） | After（source 修正） | Evidence | Disposition |
|---|---|---|---|---|
| C043 | 「InterSubMod 與 LongLineage 都不能寫 tagged BAM」忽略 revision 差異。 | 明列 `inter_sub_mod` 不寫 BAM；LongLineage public main `5daf50f` 無 writer，feature `b9aaa12` 有 `longlineage-tag-bam`；資料來源仍須看 provenance。同步修正 System／Upstream／LongLineage。 | claim inventory：feature `b9aaa12` clean build 有 binary、main `5daf50f` 無。 | LOCAL_SOURCE_CORRECTED；live publish pending |
| C059 | current `significance_summary.csv` 被寫成 193 欄且完全無版本欄。 | 改成 commit `73afaea` current 199 欄；第 187 欄 `VerificationSchemaVersion=2`、第 196 欄 `RegionStratificationSchemaVersion=1`；仍缺 whole-file layout version。同步修正 Engine／How-to。 | fresh runtime/source header；C059 evidence。 | LOCAL_SOURCE_CORRECTED |
| C107 | 公開 framing 指向 subclonal／lineage reconstruction。 | Home 與 Sidebar 採用「Local mutation-state candidate reconstruction from single-molecule sSNV co-occurrence」。 | `authority_manifest.json` `claim_boundary`。 | LOCAL_SOURCE_CORRECTED |
| C109 | 將「同一細胞譜系」寫成直接觀測。 | 改成只直接觀測同一 physical molecule 共現；cellular co-membership、順序與 lineage 仍是 model-dependent inference。同步修正 Home／System／LongLineage。 | `claim-contract-v5.md` 與 C109 evidence。 | LOCAL_SOURCE_CORRECTED |
| C112 | exact-PS 結果常省略 local／recurrence-allowed／minimum／candidate 限定。 | Home、System、Engine、LongLineage 的 scope banner 與 funnel 改用 canonical candidate terminology；88.26% 明標數學 topology、非 cellular lineage prevalence。 | `authority_manifest.json` model contract／claim boundary。 | LOCAL_SOURCE_CORRECTED |
| C113 | `170,131` 被標為一般「single sites」，grain 不明。 | 明列 `k=1 strict read-linkage components: 170,131 / 255,752 (66.52%)`。 | `denominator_registry.tsv` `k1_linkage_abstain`。 | LOCAL_SOURCE_CORRECTED |
| C142 | 66.52% 被套到「全部突變」。 | 明列 66.52% 只適用 strict-component denominator；不可套到 469,849 dataset-records。 | `170131/255752=66.5219%`；`170131/469849=36.2097%`。 | LOCAL_SOURCE_CORRECTED |
| C148 | `verify_pipeline_numbers.py` 被說成重算所有方法／Wiki 數字。 | 限定為 historical 35,332-site pipeline；列出 exact-PS、LongLineage、storage/code/test 的分開驗證來源與未覆蓋狀態。同步修正 Analysis／How-to。 | 腳本實際讀取 2026-06-27 teaching data，shown output 為 35,332-site family。 | LOCAL_SOURCE_CORRECTED |
| C149 | public clone 被承諾約 10 分鐘即可用 repo 路徑跑分析。 | 移除 public timing promise；標示 HCC1395 為 internal receipt；改成使用者自備絕對路徑與輸入 preflight；Engine 的同語意 internal paths 一併修正。 | public Git objects 缺 tumor BAM、reference FASTA、example VCF。 | LOCAL_SOURCE_CORRECTED；portable licensed fixture 仍缺 |
| C151 | `inter_sub_mod` 被說成直接產生 exact-PS funnel 核心數字。 | 分離 `inter_sub_mod`（per-region methylation/statistics）與 research `exact_ps_topology_af` C++ solver + Python runners（2026-07-24 funnel）。同步修正 Home／System／Engine。 | `authority_manifest.json` 與 solver artifact provenance。 | LOCAL_SOURCE_CORRECTED |

## Historical validator fresh run

輸入路徑：

- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/verify_pipeline_numbers.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/*.json`

執行命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts
python3 verify_pipeline_numbers.py
```

輸出：stdout-only，無新檔案；exit 0。實際輸出片段：

```text
sSNV total = 35332
TP / FP = 30490/4842  ✓
linked = 21554
underpowered = 5458
isolated_singleton = 8320
三桶加總 == total = 35332  ✓
```

這次 fresh run 確認腳本驗證的是 HCC1395 historical teaching-data family；其讀取契約與輸出中沒有 2026-07-24 exact-PS 的 469,849／255,752／98,955／71,955 funnel，也沒有 LongLineage revision、storage、source-count 或 current test-suite 的全域驗證。因此 C148 的縮限改寫有執行證據支持。

## 數值驗算

命令：

```bash
awk 'BEGIN {
  printf "170131/255752=%.4f%%\n", 100*170131/255752;
  printf "170131/469849=%.4f%%\n", 100*170131/469849;
  printf "63506/71955=%.4f%%\n", 100*63506/71955
}'
```

實際輸出：

```text
170131/255752=66.5219%
170131/469849=36.2097%
63506/71955=88.2579%
```

## Residual scan

命令：

```bash
rg -n -i \
  '重建腫瘤亞群譜系|同一個細胞譜系.*直接觀測|全部突變中有三分之二|單點無共現.*66|兩支.*都沒有寫 BAM|真正產出論文核心數字.*那支程式|原始碼 193|且無版本欄位|方法.*每一個數字重新算|約 \*\*10 分鐘' \
  docs/wiki
```

實際輸出：空；exit 1 代表沒有舊錯誤句型。正向 anchor scan 可找到：canonical headline、physical-molecule ceiling、170,131 / 255,752、`5daf50f`／`b9aaa12`、199 欄、historical 35,332-site 與 `exact_ps_topology_af`。

## Markdown／shell／diff QA

執行命令：

~~~~bash
cd /big7_disk/liaoyoyo2001/InterSubMod

# 檢查 Wiki fenced blocks、heading jumps、duplicate headings、empty links
# 與 missing relative targets；只讀，不寫檔。
python3 - <<'PY'
from pathlib import Path
import re, sys
files = sorted(Path('docs/wiki').glob('*.md'))
errors = []
links = headings = 0
for p in files:
    s = p.read_text(encoding='utf-8')
    if s.count('```') % 2:
        errors.append(f'{p}: unbalanced fenced code blocks')
    prev = 0
    seen = {}
    in_fence = False
    for n, line in enumerate(s.splitlines(), 1):
        if line.startswith('```'):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        m = re.match(r'^(#{1,6})\s+(.+?)\s*$', line)
        if m:
            level = len(m.group(1))
            headings += 1
            if prev and level > prev + 1:
                errors.append(f'{p}:{n}: heading jump H{prev}->H{level}')
            prev = level
            key = re.sub(r'\s+', ' ', m.group(2).strip().lower())
            seen[key] = seen.get(key, 0) + 1
        for label, target in re.findall(r'!?\[([^\]]*)\]\(([^)]+)\)', line):
            links += 1
            if not label.strip() or not target.strip():
                errors.append(f'{p}:{n}: empty markdown link')
            if not re.match(r'^(?:https?://|#|mailto:)', target):
                base = (p.parent / target.split('#', 1)[0]).resolve()
                if not base.exists():
                    errors.append(f'{p}:{n}: missing relative target {target}')
    for key, count in seen.items():
        if count > 1:
            errors.append(f'{p}: duplicate heading x{count}: {key}')
print(f'wiki_files={len(files)} headings={headings} markdown_links={links} errors={len(errors)}')
for e in errors:
    print('ERROR', e)
sys.exit(bool(errors))
PY

# 對所有 bash-language fenced blocks 逐一做語法剖析，不執行範例命令。
python3 - <<'PY'
from pathlib import Path
import re, subprocess, sys
errors = []
blocks = 0
for p in sorted(Path('docs/wiki').glob('*.md')):
    s = p.read_text(encoding='utf-8')
    for i, m in enumerate(re.finditer(r'```bash\n(.*?)\n```', s, re.S), 1):
        blocks += 1
        proc = subprocess.run(['bash', '-n', '-c', m.group(1)], text=True, capture_output=True)
        if proc.returncode:
            errors.append(f'{p}:bash-block-{i}: {proc.stderr.strip()}')
print(f'bash_blocks={blocks} syntax_errors={len(errors)}')
for e in errors:
    print('ERROR', e)
sys.exit(bool(errors))
PY

git diff --check -- docs/wiki
git diff --numstat -- docs/wiki
~~~~

實際輸出片段：

```text
wiki_files=9 headings=115 markdown_links=86 errors=0
bash_blocks=13 syntax_errors=0

12  5   docs/wiki/Analysis-and-Presentation.md
10  8   docs/wiki/Home.md
38  23  docs/wiki/How-to-Run.md
38  28  docs/wiki/InterSubMod-Engine.md
3   1   docs/wiki/LongLineage-Engine.md
15  13  docs/wiki/System-Overview.md
2   2   docs/wiki/Upstream-and-Data.md
3   1   docs/wiki/_Sidebar.md
```

`git diff --check -- docs/wiki`：exit 0、無輸出。

全 repo 的既有 `scripts/analysis/validate_docs_structure.py` 也曾試跑，但它無 target-only CLI，掃描 1,277 個歷史 Markdown 並回報大量既有 naming／metadata／link baseline，無法用來判定本次 Wiki delta；因此本輪另採上述 target-only checker，且不把全 repo baseline 誤記成本次 regression。

## 尚未完成的外部動作

1. 將 `InterSubMod/docs/wiki/*.md` 同步至 `InterSubMod.wiki.git` 並發布；本輪未獲 push 授權，未執行。
2. 發布有授權、含 SHA-256 與下載說明的小型 BAM／FASTA／VCF fixture；未完成前，公開 How-to 的 analysis steps 必須維持 internal-data 標示。
3. 若要公開宣稱 LongLineage tagged-BAM，需把 feature `b9aaa12` 的能力正式整合／發布到 public main，並補 schema/provenance/run receipt；目前只能寫版本差異。
