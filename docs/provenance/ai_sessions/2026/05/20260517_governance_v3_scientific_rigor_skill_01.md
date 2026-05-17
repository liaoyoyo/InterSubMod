<!--
建立時間: 2026-05-17
目標: AI session 執行報告 — governance v3 分流 + /scientific-rigor 元 skill 建立
處理範圍: 2026-05-15 至 2026-05-17 跨 session 工作
關聯檔案:
  - InterSubMod/.claude/CLAUDE.md (v3 D2 分流後)
  - InterSubMod/AGENTS.md (v3 D2 分流後)
  - InterSubMod/.claude/skills/scientific-rigor/SKILL.md (新建)
  - InterSubMod/templates/postmortem.md + research_index.md (新建)
  - /bip7_disk/liaoyoyo2001/.claude/plans/scientific-rigor-skill-draft-lazy-muffin.md (plan)
session: 跨 2 天 (2026-05-15 →ish→ 2026-05-17)
classification: ai_session_report
type: meta_governance_refactor
-->

# AI Session Report — Governance v3 分流 + /scientific-rigor 元 skill

## TL;DR（1 行）

跨 2 天會話完成 InterSubMod governance 大幅重構（CLAUDE.md/AGENTS.md D2 分流）+ 建立元方法論 skill `/scientific-rigor`（整合啟發式學習工作流映射），共 **6 commits** + **驗證 7/7 通過 4.5/5**。

---

## 摘要（3 行）

**為何**: 用戶要求重新審視 Claude Code AI harness 全架構，發現 always-loaded ~20K 違反業界 ≤5K acceptable + 6 既有嚴謹度 skill 缺級聯順序 + 啟發式學習未滲透到工作流。

**怎麼做**: 業界 4 個 researcher（context engineering / agent harness / @-import 規格 / 整合性框架）+ 3 輪 reviewer + 7 query 驗證 + plan mode 2 次 + AskUserQuestion 3 問。

**結果**: always-loaded 從 ~20K 降到 ~5.5K（-72%）；新 skill 8 cross-ref 全部雙向；驗證 7/7 通過；plan 完整記錄後續輪次。

---

## §1 完成事項清單（6 commits）

| Commit | 內容 | 規模 |
|--------|------|------|
| `696c7c1` | refactor(governance): v3 D2 分流 CLAUDE.md + AGENTS.md | 430 ins / 511 del / 2 files |
| `42217cf` | feat(skill): /scientific-rigor 元方法論 | 656 ins / 2 files（SKILL.md 513 + evals.json 143）|
| `a743b23` | chore(claude): §3 移除「（規劃中）」標 | 1 ins / 1 del |
| `c1dde00` | docs(skills): 6 既有 skill 加 /scientific-rigor cross-reference | 76 ins / 6 files |
| `a031d21` | feat(templates): postmortem + research_index | 247 ins / 2 files |
| `dce837f` | fix(skills): 驗證後補強 + 2 skill cross-ref 補齊 | 55 ins / 1 del / 3 files |

**共 +1465 行 / -513 行 / 16 個檔案**

---

## §2 關鍵決策與業界對照

### 2.1 D2 分流（CLAUDE.md vs AGENTS.md）

**業界依據**：
- AGENTS.md 為 [Linux Foundation 標準](https://agents.md/)（2025-11 開源，60k+ repos 採用）
- Anthropic Claude Code 官方支援 `@-import` 但**不省 token**（imports 計入 always-loaded）
- Cline Memory Bank 6-file dependency hierarchy 為最具體範本

**本專案決策**：
- **CLAUDE.md** 161 行：僅 Claude Code 特定（確認矩陣 / Skills / Hooks / Rules / Opus 4.7）
- **AGENTS.md** 286 行：跨 agent governance（G1-G5 / Phase 主軸 / KB / Step→Verify / IO 顯示 / 任務切割 / 回應分級 / 封存 Hard Gate）

### 2.2 /scientific-rigor 元 skill（業界創新）

**業界對照確認**：Anthropic Skills / Voyager / AutoGen 均無「meta-skill 元方法論層」概念 — 本專案因「研究領域複雜度高 + Opus 4.7 literal 行為」需要顯式 meta-layer。

**14 章結構**:
- §0.5 最小可用子集（Cognitive Load 對策）
- §1 Foundation
- §2 證據分級 5 levels + §2.1 Checklist
- §3 Effect Size 量化（Cohen's d / AUC / F1 / RRR / NNT）
- §4 因果推論 + DAG（Pearl）
- §5 對照組（引用 validation-protocol + auc-confound-guard）
- §6 消融實驗（引用 methodology-audit）
- §7 Pre-registration + 可重現性 7 項
- §8 任務壓縮 + Reflexion buffer
- §9 PDCA + Blameless Postmortem
- §10 啟發式學習工作流映射 ⭐⭐（用戶強調）
- §11 與既有 skill 12 步協作圖
- §11.5 Governance Tier（JMIR 2026 類比）
- §11.6 雙環學習（Argyris 1977）

### 2.3 啟發式學習工作流映射（§10）

**用戶要求**：「啟發式學習應廣泛套用所有工作流，不只學新主題」

**§10.1 對應映射表**:
| 啟發式概念 | 工作流套用 |
|----------|---------|
| Active Recall | 改動前先口述方案 + 預期 |
| Spaced Repetition | 結論 7d / 30d / 90d 後 spaced check |
| Feynman | 重大決策後用簡單話重述 |
| Interleaving | 跨樣本穿插驗證 |
| Deliberate Practice | NEGATIVE 走 postmortem |
| Retrieval Practice | 任務開始先回憶 |

---

## §3 業界框架引用清單

### 一手資料
- [Anthropic Claude Code Memory docs](https://code.claude.com/docs/en/memory) — @-import + ≤200 行建議
- [Anthropic Agent Skills](https://platform.claude.com/docs/en/agents-and-tools/agent-skills/overview)
- [Linux Foundation AGENTS.md](https://agents.md/)
- [Cline Memory Bank](https://docs.cline.bot/features/memory-bank)
- [Reflexion (arxiv 2303.11366)](https://arxiv.org/abs/2303.11366)
- [Voyager (arxiv 2305.16291)](https://arxiv.org/abs/2305.16291)
- [Generative Agents (arxiv 2304.03442)](https://arxiv.org/abs/2304.03442)
- [Goal Drift in LLM Agents (arxiv 2505.02709)](https://arxiv.org/abs/2505.02709)
- [Replication Crisis (OSC 2015)](https://en.wikipedia.org/wiki/Replication_crisis)
- [Argyris Double-loop Learning 1977](https://hbr.org/1977/09/double-loop-learning-in-organizations)
- [Roediger & Karpicke 2006 Active Recall](https://pubmed.ncbi.nlm.nih.gov/16507066/)
- [Paas & van Merriënboer 2020 CLT](https://journals.sagepub.com/doi/10.1177/0963721420922183)
- [Google SRE Postmortem Culture](https://sre.google/sre-book/postmortem-culture/)
- [JMIR 2026 Three-tier AI competency](https://www.jmir.org/2026/1/e86550)

### 二手評論
- Karpathy 「context engineering」(2025-06 X 推文)
- Simon Willison context engineering 部落格
- Cognition「Don't Build Multi-Agents」
- Cem Karaca 「CLAUDE.md 42K token disaster」

---

## §4 Token Budget 變化

| 階段 | always-loaded | 對照業界 ≤2K ideal / ≤5K acceptable |
|------|--------------|-----------------------------------|
| **改前** | ~10K（用戶估）→ 真實 ~20K（含 rules/ 永遠載入）| 違反 |
| **改後** | ~5.5K | acceptable ✓ |
| **降幅** | **-72%** | 仍可進步空間：rules paths frontmatter 加入後可降至 ~3K |

---

## §5 驗證結果（Step D）

**7 fresh-session queries 通過率 7/7** + reviewer 評 **4.5/5**

| Query | 章節 | 評分 |
|------|------|------|
| 1 證據分級 5 levels | §2 | ⭐⭐⭐⭐⭐ |
| 2 F1 +0.0112 結論驗證 | §3+§4+§5+§0.5 | ⭐⭐⭐⭐ |
| 3 30d Spaced Recall | §10.1+§10.2+§2 | ⭐⭐⭐⭐ |
| 4 新研究方向準備 | §7.1+templates | ⭐⭐⭐⭐ |
| 5 NEGATIVE postmortem | templates+§9.2+§11.6 | ⭐⭐⭐⭐ |
| 6 程式修改前 active recall | §10.2+AGENTS §6 | ⭐⭐⭐⭐⭐ |
| 7 DAG collider audit | §4+/known-pitfalls P-01 | ⭐⭐⭐⭐⭐ |

---

## §6 雙向 cross-reference 完整性（8 skills）

| Skill | scientific-rigor cross-ref |
|-------|---------------------------|
| /known-pitfalls | ✅ 6 refs |
| /methodology-audit | ✅ 6 refs |
| /verification-loop | ✅ 6 refs |
| /validation-protocol | ✅ 6 refs |
| /fast-learning-coach | ✅ 4 refs |
| /memory-consolidation | ✅ 7 refs |
| /check-staleness | ✅ 4 refs |
| /auc-confound-guard | ✅ 7 refs |

---

## §7 plan 完整記錄（後續輪次）

`/bip7_disk/liaoyoyo2001/.claude/plans/scientific-rigor-skill-draft-lazy-muffin.md`

### Step F.3（評估後不建議立刻做）
- `evidence_level_lint.sh` / `causal_claim_check.sh` / `skill_change_audit.sh` 三 hook
- 拒絕原因: 22 hooks 已逾管理閾值 + Cognitive Load 風險
- skill_change_audit 可之後考慮

### Step F.4（P2/P3 業界空缺）
- Cynefin 域顯式化（融入 /confirmation-protocol）
- Productive Failure 顯式化（已部分在 templates/postmortem.md Reopen Threshold）
- Living-update 機制（git commit prefix + monthly audit）

### M2-M5 後續矛盾
- M2: 41 skills 視覺分層
- M3: 4 觸發機制可預測性
- M4: Memory expiry
- M5: KB ↔ Memory 邊界

### Q4-Q10 技術債
- Q4: 152 allow list 重構
- Q5: biorxiv/ensembl MCP 僵屍驗證
- Q9: archive/ 過時 agent 清理
- Q10: MEMORY.md 精簡（目前 113 行，未達截斷但可降 Concluded 區）

---

## §8 教訓（feedback memory 候選）

### 8.1 對 always-loaded 估計不準 — Rules 永遠載入

**事實**: `.claude/rules/*.md` 預設**永遠載入**（不是條件式），加 `paths:` frontmatter 才條件化。

**影響**: 原以為 CLAUDE.md 15KB → 真實 ~20K（含 4 rules）。

**教訓**: 任何 always-loaded budget 估計必含 rules/。

### 8.2 Spawn subagent 不能用來測 CLAUDE.md 即時載入

**事實**: subagent 繼承 main agent 啟動快取的 CLAUDE.md，不會 fresh read。

**影響**: 我曾嘗試 Edit CLAUDE.md + spawn subagent 驗證 → subagent 仍看舊版。

**教訓**: 載入機制驗證必須開新 session（fresh boot），不能靠 subagent。

### 8.3 「定論」措辭陷阱

**用戶教訓**: 「Phase 1A F1 已鎖定 + TO 模式甲基化增益為負」我寫成定論 — 用戶糾正「歷史觀察非定論，仍應深入研究」。

**修正**: SKILL.md §2 + AGENTS.md §3 全文禁用「鎖定」「定論」「已證實」假裝 L1 等級。

### 8.4 元 skill 「meta-methodology」是業界創新

**事實**: Anthropic Skills / Voyager / AutoGen 都是 flat skill library，無 meta-skill 層。

**本專案需求**: 研究領域複雜度高 + 6 嚴謹度 skill 缺級聯順序 → 元 skill 提供協作圖。

**警告**: Cognitive Load 風險 — §0.5 最小可用子集為對策。

---

## §9 給下個 session 的交接

### 立即可推進（按優先級）

| 優先 | 動作 | 預估 |
|------|------|------|
| **P0** | F.4 P2 Cynefin 顯式化（融入 /confirmation-protocol）| 30 min |
| **P0** | Q5 biorxiv/ensembl MCP 僵屍驗證 | 15 min |
| **P1** | M2 41 skills 視覺分層 | 60 min |
| **P1** | Q10 MEMORY.md Concluded 區降級 grep-only | 30 min |
| **P2** | Q4 152 allow list 重構 | 60 min |
| **P2** | M4 Memory expiry 機制 | 60 min |
| **P3** | M3 4 觸發機制可預測性 | 60 min |

### 暫不做（評估後拒絕或留待）
- F.3 三 hook（Cognitive Load 風險）
- M5 KB ↔ Memory 邊界（短期內無實質衝突）

### 用戶驗證待跑（建議下次 fresh session）

3 個 query 跑新 skill：
1. 「請列出 /scientific-rigor §2 證據分級 5 levels」
2. 「我要驗證 F1 提升 +0.0112 的結論」
3. 「30 天前的 NEGATIVE 結論該怎麼複查？」

---

## §10 相關檔案路徑速查

- **Plan**: `/bip7_disk/liaoyoyo2001/.claude/plans/scientific-rigor-skill-draft-lazy-muffin.md`
- **新 skill**: `InterSubMod/.claude/skills/scientific-rigor/SKILL.md` + `evals.json`
- **Templates**: `InterSubMod/templates/postmortem.md` + `research_index.md`
- **Governance**: `InterSubMod/.claude/CLAUDE.md` (161 行) + `InterSubMod/AGENTS.md` (286 行)
- **Backup**: `InterSubMod/docs/drafts/backup/{CLAUDE,AGENTS}_md_prod_20260517_pre_v3deploy.md`
- **Drafts (中間設計版本)**: `InterSubMod/docs/drafts/20260516_*` + `20260517_*`
- **舊 HTML 整體審查**: `InterSubMod/docs/reports/validated/2026/05/20260515_claude_code_AI_harness_full_architecture_01.standalone.html`

---

> **Session 完成時間**: 2026-05-17（橫跨 2 天）
> **Branch**: `refactor/phase1-safety`
> **HEAD**: `dce837f`
> **未 push**: 等用戶決定（依規則不主動 push）
