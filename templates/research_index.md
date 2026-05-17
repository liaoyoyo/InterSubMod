<!--
建立時間: 2026-05-17
目標: 新研究專案 00_INDEX.md 標準範本 — 含 /scientific-rigor §7 Pre-registration 3 欄
處理範圍: research/<topic>/00_INDEX.md 起手樣板
關聯檔案:
  - InterSubMod/.claude/skills/scientific-rigor/SKILL.md §7 Pre-registration + §10.2 retrieval practice
  - InterSubMod/.claude/skills/init-research/SKILL.md (多週級長期專案 scaffolding)
  - InterSubMod/.claude/skills/cycle-init/SKILL.md (短週期 hypothesis cycle)
-->

# Research Project: <topic> — 00_INDEX.md Template

> **用法**: 複製本檔到 `InterSubMod/research/<topic>/00_INDEX.md` 並填入。
> **觸發**: 新研究方向開跑前強制走（`/scientific-rigor §7.1 Pre-registration`）。

---

# <Topic>: <one-line description>

**Created**: YYYY-MM-DD
**Owner**: <user / agent>
**Status**: planning / in-progress / paused / concluded
**Tier**: ⭐⭐ (planning) → ⭐⭐⭐ (preliminary) → ⭐⭐⭐⭐ (validated) → ⭐⭐⭐⭐⭐ (publishable)
**Cycle IDs**: <list>

---

## §1 Pre-registration（Confirmatory — `/scientific-rigor §7.1` 強制 3 欄）

> **規則**: 達 NO-GO 條件 → verdict **不可事後改寫**（呼應 `CLAUDE.md §1 Hard Gate`）
> **commit hash 為錨點**: 註冊後第一個 commit hash 是「事先」的證明

| H 預測（confirmatory） | 否證條件（falsification） | Decision threshold |
|------------------------|-------------------------|---------------------|
| <如: ISM 在 NormalASM 場景 inter > intra> | <如: 跨 5 樣本 ≥3 樣本反向> | <如: inter-intra delta < 0 in ≥3 samples → NO-GO> |
| <H2 ...> | <... > | <... > |

**Exploratory（事後分析）**: 任何不在上表的觀察 → 必須標 `exploratory` tag，**不能宣告 confirmatory 結論**。

---

## §2 服務的研究目標（G1-G5 對應 — AGENTS.md §3）

- [ ] **G1**: longphase 家族 + ISM 整合突破
- [ ] **G2**: longphase-s + longphase-to 協同 > 單獨任一
- [ ] **G3**: ISM read-level epigenetic 給領域突破
- [ ] **G4**: 多樣本一致性 + reproducibility
- [ ] **G5**: 業界級貢獻（可被外部驗證）

（勾選對應的 G — 每小任務必須回答「服務哪個 Gx」）

---

## §3 可重現性檢核（`/scientific-rigor §7.2`）

- [ ] **seed 固定**（若涉及隨機）: `<seed_value>`
- [ ] **環境記錄**: commit hash `<hash>` + conda env `<env_name>` / Docker `<image>`
- [ ] **數據版本**: `<dataset>` + KDE-corrected `<yes/no>` + subsample ratio `<%>`
- [ ] **參數檔案 snapshot**: `params.json` / `config.yaml` 路徑
- [ ] **中間產出存檔位置**: `<path>`
- [ ] **命令 + 退出碼**: 記錄於 `scripts/run.sh` + log
- [ ] **期望輸出片段**: 對比實際（依 `AGENTS.md §11 IO 顯示規則`）

---

## §4 啟動研究 5 問檢核（`AGENTS.md §4`）

- [ ] Thread D（read-level epigenetic）相關？
- [ ] Thread B 撤回範圍內？
- [ ] 資料 KDE-corrected？
- [ ] 需要 VCF caller AF（非 merged AF）？
- [ ] 觸及 長計算 / C++ / 檔案搬移 / NO-GO gate？

---

## §5 預期高影響場景檢核（`AGENTS.md §5`）

- [ ] 研究重點排序 → 列候選方向 + 收益/風險
- [ ] 假說選擇 → 寫出前提 + 可能 confound
- [ ] 統計方法選擇 → 說明為何選 + 替代方案
- [ ] 改進/優化模糊指令 → 要求用戶定義成功標準
- [ ] 多檔案/多步驟重構 → 列受影響檔案 + 預期改動

---

## §6 子目錄結構

```
research/<topic>/
├── 00_INDEX.md              ← 本檔
├── README.md                ← 給人類讀者的入口
├── figures/                 ← 圖表
├── data/                    ← 中間數據
├── scripts/                 ← 分析腳本
├── reports/                 ← 報告
├── REFLECTION.md            ← `/scientific-rigor §8.3` Reflexion buffer（NEGATIVE 必建）
└── DAG/                     ← `/scientific-rigor §4` 因果圖 mermaid
```

---

## §7 進度與證據鏈

| Cycle ID | Phase | Verdict | Tier | Evidence | Notes |
|---------|-------|---------|------|----------|-------|
| cycle_YYYYMMDD-HHMM-slug | P0-P6 | POSITIVE / NEGATIVE / CONDITIONAL | ⭐ | `state/cycles/.../evaluation.json` | ... |

---

## §8 後續

**若 verdict = POSITIVE**: → `/conclude-research` P6 收尾 + memory 更新 + KB 紀錄
**若 verdict = NEGATIVE / NO-GO**:
1. 走 `InterSubMod/templates/postmortem.md` 範本
2. 寫 `REFLECTION.md` (`/scientific-rigor §8.3`)
3. 若要 reopen → 滿足 postmortem 中的 Reopen Threshold 才可

---

## §9 相關 skill 觸發索引

| 階段 | 觸發 skill |
|------|----------|
| 啟動 | `/known-pitfalls` 查歷史陷阱 + `/scientific-rigor §7` Pre-reg |
| 假說設計 | `/inject-hypothesis` + `/research-loop` P1 |
| 執行 | `/cycle-init` P0 → `/check-staleness` P2 → `/feature-layered-observation` P3 |
| 驗證 | `/auc-confound-guard` + `/validation-protocol` L1-L4 |
| 跨樣本 | `/multi-sample-consistency` P4 |
| 評估 | `/run-evaluator` P5 |
| 收尾 | `/conclude-research` P6 |
| 失敗 | `templates/postmortem.md` + `/scientific-rigor §8.3 §9.2 §11.6` |
| 30d 後 | `/scientific-rigor §10.1 Spaced Repetition` |
