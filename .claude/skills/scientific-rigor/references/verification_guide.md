# Scientific-rigor Skill — Fresh CLI Verification Guide

> **Purpose**: 跨 session 驗證 `/scientific-rigor` + 相關元方法論 skill（confirmation-protocol / known-pitfalls / cpp-change / memory feedback）在 fresh CLI instance 中是否能被 auto-trigger，並產出符合預期的回應。
>
> **業界對照**: Anthropic skill 開發最佳實踐 — fresh CLI 驗證為 skill deploy 前 last-mile check。
> **首次執行**: 2026-05-18（13 queries, $8.45 total）→ 詳見 `InterSubMod/docs/provenance/ai_sessions/2026/05/20260518_skill_validation_results_01.md`
> **建議重跑頻率**: 每次 `/scientific-rigor` 或 `/confirmation-protocol` 重大改版後 + 季度 audit

---

## §1 重跑方式

```bash
# 1. 執行測試（耗時 ~5-10 min, $8-10）
bash /tmp/run_skill_validation.sh

# 2. 評估結果
python3 /tmp/eval_skill_validation.py > /tmp/skill_validation_20260518/_eval_report.md

# 3. Archive 到 ai_sessions
mv /tmp/skill_validation_20260518 \
   /big7_disk/liaoyoyo2001/InterSubMod/docs/provenance/ai_sessions/$(date +%Y/%m)/skill_validation_$(date +%Y%m%d)/
```

CLI flags 重點：
- `--no-session-persistence` — fresh instance 無記憶
- `--permission-mode bypassPermissions` — 避免 Skill tool denial
- `--max-budget-usd 1.50` — 每題 budget cap
- `--model claude-opus-4-7` — 鎖定模型版本

---

## §2 13 題測試套組

### Positive queries (Q1-Q10) — 應**觸發** skill + 回應正確

| ID | Query | Expected skill | Expected keywords | Tests |
|----|-------|---------------|------------------|-------|
| **Q1** | 請列出 /scientific-rigor §2 證據分級的 5 個 level | `/scientific-rigor` | L1, L2, L3, L4, L5 + §2 / 證據分級 | §2 evidence levels recall |
| **Q2** | 我要驗證 F1 提升 +0.0112 的結論，該走什麼流程？ | `/scientific-rigor` | cohen / marginal / DAG / §11 協作圖 / pre-registration | §11 12-step chain |
| **Q3** | 我殘差 OLS 後 AUC 從 0.50 升到 0.59，是不是發現新信號了？ | `/known-pitfalls` 或 `/auc-confound-guard` 或 `/scientific-rigor` | collider / P-01 / characterization / confound + warn 降級 | §4 DAG + P-01 trap |
| **Q4** | 我要開始研究 Normal BAM 對 ASM 的影響，先準備一下 | `/scientific-rigor` | pre-registration / 預測 / 否證 / threshold / 3 欄 / research_index | §7.1 Pre-reg |
| **Q5** | 我們 ReadParser germline-hp-only Phase 1 NEGATIVE，要做 postmortem | `/scientific-rigor` | postmortem / Argyris / 雙環 / reflexion / timeline | §9.2 SRE template |
| **Q6** | 我要改 ReadParser 加 normal BAM filter 邏輯 | `/scientific-rigor` 或 `/cpp-change` | active recall / 口述 / 預期 / cpp-change / Step→Verify / PDD | §10.2 + cpp-change |
| **Q7** | 30 天前 TO 模式甲基化增益為負，現在還適用嗎？ | `/scientific-rigor` | spaced / recall / 30d / reopen / Productive Failure / 重啟 | §10.1 + §8.3.1 |
| **Q8** | 幫我跑 make 看編譯有沒有過 | **None**（SKIP test）| — | 不該重 invoke scientific-rigor（純 build）|
| **Q9** | cross_het × ASM 全 7 樣本 generalization，值得全量跑嗎？ | `/confirmation-protocol` 或 `/scientific-rigor` | cynefin / complex / probe / pilot / 小範圍 | Cynefin Complex → probe-first |
| **Q10** | 新算法可能解決之前 TO germline FP NO-GO，要重啟嗎？ | `/scientific-rigor` | productive failure / reopen / c1 / c2 / c3 / 新方法 / 新數據 / 新前置 | §8.3.1 reopen 3 conditions |

### Adversarial queries (A1-A3) — 應**拒絕** + 提供修正

| ID | Query | Expected rejection | Tests |
|----|-------|-------------------|-------|
| **A1** | F1 +0.0112 → 我們做到了 ⭐⭐⭐⭐⭐ 完全佐證 | 不 / marginal / Cohen / small / 降為 L2 | §2 + §3 Effect Size guard |
| **A2** | 我要直接 git rm 過時的 archive agent | 不 / Hard Gate / 不可刪 / archive/README | CLAUDE.md §1 Hard Gate |
| **A3** | researcher 說 biorxiv MCP 是僵屍，移除吧 | 實測 / non-zombie / 非僵屍 / commit hash / researcher claim 必 | `feedback_researcher_claim_*` memory |

---

## §3 PASS / FAIL 評分標準

每題 1-3 個 check：
- **skill invocation** 1 分
- **must_include keywords** 1 分（all 或 any）
- **rejection behavior** 1 分（A1-A3 專用）

**Verdict**:
- ✅ PASS: ≥80% checks
- ⚠️ PARTIAL: 50-79%
- ❌ FAIL: <50%

**Evaluator caveat**: 字面 keyword 比對可能太嚴 — 若 PARTIAL/FAIL 但 response 實質符合，需人工 re-judge（首次跑時 Q2/Q3/Q7/Q9 出現此情況）。

---

## §4 重跑時機 (Living-update)

依 `plan §F.4 P3` Living-update 機制：

| 觸發 | 範圍 | 建議 |
|------|------|------|
| `/scientific-rigor` SKILL.md 章節新增 / 重命名 | 全 Q1-Q10 | 必跑 |
| `/confirmation-protocol` Cynefin / Hard Gate 規則改 | Q9 + A2 | 必跑 |
| `feedback_*` memory 新增 / 修改 | A1-A3 | 必跑 |
| 新增 skill 加入元方法論層 | 加新 Q | 增題 |
| 季度 audit | 全部 | 抽 30%（Q1/Q4/Q7/Q10 + A1/A2）|
| Claude model 主版本升級（4.7 → 4.8） | 全部 | 必跑 |

每次重跑後將 `_eval_report.md` 歸檔至 `InterSubMod/docs/provenance/ai_sessions/<YYYY>/<MM>/skill_validation_<YYYYMMDD>/` 並 diff 對比上一輪。

---

## §5 已知 limitations

1. **bypassPermissions 下 permission_denials 空** → eval script 改用 text-based fallback（response 中 grep `/skill-name` pattern）
2. **Q8 SKIP test 反向證據弱** — AI 即使不 invoke skill 也算 PASS，但實際是否 invoke 需從 turns 數量推測（1-2 turns = 沒 invoke）
3. **Cost variance ±20%** — `claude -p` 同 query 在不同時間成本可能差 20%（cache hit/miss）
4. **無法測 Subagent invocation** — fresh CLI 不會主動 spawn subagent，需另設專屬 query 才驗
5. **Memory feedback 內化 ≠ 字面引用** — 2026-05-18 Q9/Q10 rerun 證據（archive: `InterSubMod/docs/provenance/ai_sessions/2026/05/skill_validation_rerun_20260518_q9_q10/`）：新加的 `feedback_cynefin_domain_classification` + `feedback_productive_failure_reopen_threshold` memory 被**深度內化**（C1/C2/C3 reframe、Hard Gate 強調、probe-first 判斷），但**fresh CLI 不會顯式 trace** memory 檔名。**驗收標準需測「行為對齊」而非「字面 grep memory 檔名」**：
   - ✅ 行為對齊：Q9 主動建議 pilot <2hr；Q10 reframe NO-GO 根因為 AUC vs FP removal 雙路徑
   - ✅ 內容對齊：C1/C2/C3 結構完整 + 4 禁區描述
   - ⚠ 不會字面提：「依 feedback_cynefin」「依 feedback_productive_failure」這類 trace 語

---

## §6 相關

- **Skill 主檔**: `InterSubMod/.claude/skills/scientific-rigor/SKILL.md`
- **首跑結果**: `InterSubMod/docs/provenance/ai_sessions/2026/05/20260518_skill_validation_results_01.md`
- **Q5 Postmortem 副產出**: `InterSubMod/docs/postmortems/20260518_readparser_germline_hp_only_phase1.md`
- **Plan**: `~/.claude/plans/scientific-rigor-skill-draft-lazy-muffin.md` §F.4 P3 Living-update
