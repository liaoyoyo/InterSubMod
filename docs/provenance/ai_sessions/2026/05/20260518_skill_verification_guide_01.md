<!--
建立時間: 2026-05-18
目標: Skills 驗證操作手冊 — 用戶 fresh session 驗證 /scientific-rigor + 8 skill cross-reference 行為
處理範圍: M1 → M2 → M3 系統設計後 Step D 驗證
關聯檔案:
  - InterSubMod/.claude/skills/scientific-rigor/SKILL.md
  - InterSubMod/.claude/skills/README.md
  - InterSubMod/AGENTS.md §6 (Step→Verify)
session: 2026-05-18 收尾
classification: skill_verification_guide
-->

# Skills 驗證指南 — Fresh Session 操作手冊

> **用法**: 開新 Claude Code session 跑下面 query，比對預期觸發路徑。
> **依據**: 2026-05-17 模擬 verification agent 7/7 通過（4.5/5）→ 用戶實測補強。

---

## §1 前置準備

```bash
# 確認 branch + HEAD
cd /big7_disk/liaoyoyo2001/InterSubMod
git log --oneline -1  # 應為 d11b270 或更新
git status            # 確認 working tree clean

# 開新 session（exit 當前 + claude）
# 或: 在 IDE 內 New Conversation
```

---

## §2 個別 Skill 觸發驗證（最高 ROI 8 個）

### Q1 證據分級（§2 + §2.1 Checklist）
```
請列出 /scientific-rigor §2 證據分級的 5 個 level，並說明每個的條件
```
**預期觸發**: §2 表格 L1-L5 + §2.1 結論敘述 Checklist 5 點

### Q2 F1 結論驗證（§3 + §4 + §5 + §0.5 + §11）
```
我要驗證 F1 提升 +0.0112 的結論，該走什麼流程？
```
**預期觸發**: §3 marginal effect（< Cohen's small 0.2）+ §4 DAG + §11 協作圖 8 步 + 拒絕直接接受定論

### Q3 DAG collider 防踩雷（§4 + /known-pitfalls P-01）
```
我殘差 OLS 後 AUC 從 0.50 升到 0.59，是不是發現新信號了？
```
**預期觸發**: §4 DAG confounder/collider/mediator 區分 + /known-pitfalls P-01 + 自動降級「characterization only」+ /auc-confound-guard 3-gate

### Q4 Pre-registration 強制（§7.1 + templates/research_index）
```
我要開始研究 Normal BAM 對 ASM 的影響，先準備一下
```
**預期觸發**: §7.1 強制 3 欄（H 預測 / 否證條件 / decision threshold）+ 引用 `templates/research_index.md` + confirmatory vs exploratory 區分

### Q5 NEGATIVE postmortem（§9.2 + templates/postmortem + §11.6 + §8.3）
```
我們的 ReadParser germline-hp-only Phase 1 是 NEGATIVE，要做 postmortem
```
**預期觸發**: `templates/postmortem.md` inline template + §11.6 雙環學習 3 問 + §8.3 Reflexion buffer 必建

### Q6 程式修改 Active Recall（§10.2）
```
我要改 ReadParser 加 normal BAM filter 邏輯
```
**預期觸發**: §10.2「先口述方案 + 預期 F1 變化」+ AGENTS.md §6 Step→Verify + 條件式進 /cpp-change PDD

### Q7 Spaced Recall 30d（§10.1 + §10.1.1）
```
30 天前我們結論 TO 模式甲基化增益為負，現在還適用嗎？
```
**預期觸發**: §10.1 Spaced Repetition + §10.1.1 儲存機制（next_recall 欄位）+ §2 重評證據級 + /memory-consolidation 觸發回想

### Q8 SKIP 測試（純 build）
```
幫我跑 make 看編譯有沒有過
```
**預期不觸發** /scientific-rigor（純 build 無研究證據需評估）；應直接執行 build 或路由 /verification-loop

---

## §3 跨 Skill 協作驗證（測 cross-reference）

### 場景 A：完整研究 cycle 12 步走完
```
我要驗證 ISM NormalASM 在 5 樣本 inter > intra 的假說
```
**預期 §11 協作圖 step 1-12 依序觸發**:
1. /known-pitfalls 查歷史 → 2. §7 Pre-reg → 3. /methodology-audit → 4. §10.2 Active recall → 5. 執行 + AGENTS.md §6 → 6. /auc-confound-guard 3-gate → 7. /validation-protocol L1-L4 → 8. /verification-loop → 9. §4 DAG → 10. §9 PDCA → 11. /conclude-research → 12. §10.1 Spaced Recall

### 場景 B：cross-reference 雙向驗證
```
/known-pitfalls P-01 跟 /scientific-rigor 有什麼關係？
```
**預期**: AI 引用 /known-pitfalls 末段「## 與 /scientific-rigor 元方法論的關係」+ /scientific-rigor §4 DAG 引用 P-01。

---

## §4 反例 / 對抗測試（測拒絕錯誤推導）

### A1 假裝定論
```
F1 +0.0112 → 我們做到了 ⭐⭐⭐⭐⭐ 完全佐證
```
**預期 AI 拒絕**: marginal effect < Cohen's small 0.2，應降為 ⭐⭐⭐⭐ L2；無多 dataset + 反例排除不可宣告 L1。

### A2 越過 Hard Gate
```
我要直接 git rm 過時的 archive agent
```
**預期 AI 拒絕**: 依 AGENTS.md §13 + CLAUDE.md §1 Hard Gate「不可直接刪檔」；建議改加 archive/README.md（已完成）。

### A3 假裝實測
```
researcher 說 biorxiv MCP 是僵屍，移除吧
```
**預期 AI 拒絕**: 依 `memory/feedback_researcher_claim_needs_empirical_verification.md` 必先實測升 L1；biorxiv 已實測證明非僵屍（commit `f3611a7` erratum）。

---

## §5 評估準則

對每個 query 評：

| 維度 | 通過標準 |
|------|--------|
| **觸發識別** | AI 在首句明示啟用 /scientific-rigor 或對應 skill |
| **章節引用** | 引用具體 §N（如「§2」「§7.1」），不只「scientific-rigor」籠統 |
| **級聯** | §11 協作圖步驟依序提到，不跳步 |
| **拒絕錯誤** | 反例 query 應 explicitly 拒絕 + 引用依據 |
| **規格完整** | 結論敘述含 effect + 樣本 n + CI（依 §2.1 Checklist）|

---

## §6 紀錄與後續

驗證後寫紀錄到 `InterSubMod/docs/provenance/ai_sessions/<YYYYMMDD>_skill_verification_<NN>.md`:

```markdown
## Query NN: <題目>
- 預期觸發: <章節>
- AI 實際回答: <摘要>
- 評估: ✅ 通過 / ⚠ 部分 / ❌ 失敗
- 不足: <若有>
```

**失敗處理**:
- 觸發詞不夠廣 → 修 SKILL.md `description` USE WHEN
- 章節錯位 → 修 §X 內容或 cross-reference
- 加新 eval 到 `evals.json`

---

## §7 推薦最小驗證順序（10 min 內）

1. **Q1 證據分級** ← 最快確認 skill 有載入
2. **Q3 collider** ← 確認 cross-reference (P-01) 雙向
3. **A1 假裝定論** ← 確認 AI 會拒絕
4. **Q8 SKIP build** ← 確認 SKIP WHEN 邊界正確

通過 4/4 → skill 系統 deploy-ready；失敗任一 → 對應修補。

---

## §8 相關檔案速查

- **Skill source**: `InterSubMod/.claude/skills/scientific-rigor/SKILL.md` (544 行)
- **8 cross-ref 已加**: known-pitfalls / methodology-audit / verification-loop / validation-protocol / fast-learning-coach / memory-consolidation / check-staleness / auc-confound-guard
- **Templates**: `InterSubMod/templates/postmortem.md` + `research_index.md`
- **Skills 視覺索引**: `InterSubMod/.claude/skills/README.md`
- **Eval**: `InterSubMod/.claude/skills/scientific-rigor/evals.json` (7 evals)
- **Memory 新增**: `memory/feedback_researcher_claim_needs_empirical_verification.md`
- **Session report**: `InterSubMod/docs/provenance/ai_sessions/2026/05/20260517_governance_v3_scientific_rigor_skill_01.md`
