# 風格 Skills 遷移清單（17 個 baseline）

> 來源：`/big7_disk/liaoyoyo2001/InterSubMod/.claude/skills/`
> 目標：`<new-project>/.claude/skills/`

---

## 拷貝指令

```bash
SRC=/big7_disk/liaoyoyo2001/InterSubMod/.claude/skills
DST=<new-project>/.claude/skills

PORTABLE_SKILLS=(
  # 元方法論（7）
  confirmation-protocol
  scientific-rigor
  fast-learning-coach
  problem-framing-ideation
  pre-decision-audit
  known-pitfalls
  infra-ops

  # 敘述框架（1）
  narrative-frame

  # 報告生成（7）
  weekly-report
  structured-tech-report
  results-report
  report
  myPPT
  pptx-build
  html-report-build

  # 視覺化（2）
  image-gen
  image-vision-check

  # 文件管理（4）
  doc-standards
  memory-consolidation
  data-audit
  citation-verification
)

mkdir -p "$DST"
for s in "${PORTABLE_SKILLS[@]}"; do
  if [ -d "$SRC/$s" ]; then
    cp -r "$SRC/$s" "$DST/"
    echo "  cp -r $s"
  else
    echo "  MISSING: $s"
  fi
done
```

---

## 清洗檢查（拷貝後執行）

```bash
cd <new-project>/.claude/skills

# 1. 找出殘留專案 reference
grep -rn "InterSubMod\|HCC1395\|COLO829\|H1437\|H2009\|HCC1937\|HCC1954\|HCC1395_5kHz\|big7_disk\|big8_disk\|bip7_disk\|longphase-to-mod\|ClairS-TO\|paired_pileup\|evidence_ledger\|cycle_id" *.md */SKILL.md 2>/dev/null

# 2. 人工檢視並改寫為通用敘述
#    - 多在 known-pitfalls 與 narrative-frame 的 reference 案例
#    - 在 scientific-rigor §8/§9 postmortem 範例

# 3. 改寫 hook 引用路徑（若 skill 內提到 scripts/hooks/...）
grep -rn "scripts/hooks/" *.md */SKILL.md 2>/dev/null
```

---

## Skill 用途速查

| Skill | 一行用途 |
|-------|---------|
| `confirmation-protocol` | Hard Gate / Gate / Review / FYI 4-tier 確認協議；CLAUDE.md §1 暫停判斷矩陣的細節版 |
| `scientific-rigor` | 證據分級 L1-L5 / Effect Size / Pre-registration / Postmortem / Spaced recall 元方法論 |
| `fast-learning-coach` | 費曼 + 主動回想 + 80/20 + 間隔重複，dialog-driven 學習新主題 |
| `problem-framing-ideation` | 5W1H + gap analysis 收斂 1-3 可註冊假說；限 <2hr |
| `pre-decision-audit` | 決定前審查（≤30min 7 outputs + Cynefin gate + 5-dim credibility + GO/PROBE/NO-GO） |
| `known-pitfalls` | 已知陷阱清單；6 大類速查；外部 claim 對照必查 |
| `infra-ops` | 磁碟/基礎設施 preflight（防 /tmp 寫滿、長 pipeline OOM） |
| `narrative-frame` | **取代固定範本** — 50+ framework catalog + N1-N6 動態挑選 |
| `weekly-report` | 7 Phase W1-W7 週進度報告（PI 受眾） |
| `structured-tech-report` | 13 段技術報告（Toyota A3 + ADR + SRE Postmortem + Diátaxis） |
| `results-report` | 實驗結果報告（決策導向；前置需 results-analysis） |
| `report` | AI 對話執行報告（會話結束前撰寫） |
| `myPPT` | PPT 場景路由器（識別週報 / 母稿 / 新製作 → 委派 sub-skill） |
| `pptx-build` | PPTX 製作（outline → section → slide + Tier 1/2/3 分流 + Vision 10-checkpoint review） |
| `html-report-build` | LLM-direct HTML 3 模式（report / slide / standalone PI 終版） |
| `image-gen` | 雙軌生圖（codex AI 概念圖 + pycairo flow/icon） |
| `image-vision-check` | 6 維 vision quality check（4/6 通過；fail 觸發重生） |
| `doc-standards` | 文件命名規範 / partial flag ribbon / metadata template |
| `memory-consolidation` | 記憶生命週期管理（status / merge / 降級 / MEMORY.md <200 行） |
| `data-audit` | 研究輸出組織完整性檢核（圖片連結 / 索引覆蓋 / 命名合規） |
| `citation-verification` | 學術引用驗證（WebSearch + Google Scholar 驗證後才寫入 .bib） |

---

## 新專案可選 — 域 specific skill 補位建議

| 新任務領域 | 建議補位 skill |
|-----------|---------------|
| Web app | `feature-dev` (plugin) + `frontend-design` (plugin) + 自製 `e2e-test` skill |
| 數據分析 | 自製 `data-pipeline` + `eda-checklist` skill |
| 學術寫作 | 自製 `paper-outline` + 既有 `citation-verification` + `structured-tech-report` |
| 工具 / CLI 開發 | `agent-sdk-dev` (plugin) + 自製 `cli-ux-review` skill |

---

## 注意 — 不遷移的 skills（共 23 個 ISM 研究專用）

`auc-confound-guard` `feature-layered-observation` `multi-sample-consistency` `pivot-direction` `inject-hypothesis` `init-research` `review-evidence` `observation-analysis` `cycle-init` `research-loop` `check-staleness` `run-evaluator` `conclude-research` `validation-protocol` `cycle-state` `research-dashboard` `research-context-loader` `provenance-tier-audit` `methodology-audit` `verification-loop` `cpp-change` `results-analysis` `implementation-notes` `html-preview`(deprecated)

→ 純 ISM 7-Phase Waterfall + C++ 流程；新任務若是 C++ 專案可參考 `cpp-change` 改寫，但**不直接拷貝**。
