<!--
建立時間: 2026-05-18
目標: V6 production tag Tier 1.2 — 5-day workflow 壓縮到 4 個工作日（W3 deadline 5/22 週五）
parent: InterSubMod/research/selfphasing_v6_production/00_PLAN.md
classification: execution_workflow
-->

# V6 Production Tag — 4-Day Compressed Workflow

> **Deadline**: 2026-05-22 週五（W3 結束）
> **今天**: 2026-05-18 週一
> **可用工作日**: Mon/Tue/Wed/Thu/Fri = 5 天但 5-day workflow 壓縮為 **4 active days + Fri buffer**
> **Critical path**: COLO829 V6 ISM (~2hr) + git tag (🔴) + PI email (🔴)

---

## §0 壓縮策略 — Day 1+2 並行

原 5-day workflow 順序執行；壓縮版讓 **COLO829 ISM 計算（2hr） 與 errata content review 並行**，省一個工作日。

```
原計劃 (5 day):  COLO829 → marker → F1 → tag → email
                ─────────────────────────────────────►

壓縮版 (4 day):  COLO829 ──┐  → tag ─┐
                          ├─ marker/F1 ─┐
                  errata ──┘            ├─ email
                                        ─────►
```

---

## §1 Day 1 (週一 5/18) — 並行啟動

### §1.1 Track A: COLO829 V6 ISM 啟動（背景跑）
**Deliverable**: `output/canonical/COLO829/V6/` 含 `step1_master` row

**AI 動作**：
1. 確認 V6 binary commit hash（last green commit on V6 branch）
2. 確認 Archive TO source: `output/canonical/COLO829/TO/` exists
3. 啟動 ISM run（依 `InterSubMod/scripts/run_vcf_all_snv.sh --sample COLO829 --binary v6 --kde-corrected`）
4. `run_in_background: true` 讓 ISM 跑（~2hr）

**🟡 列假設**: 假設 COLO829 ONT R10 無 MM/ML tag → V6 marker coverage 可能偏低 → sample-level gating 標註

### §1.2 Track B: Manifest 草稿 + Errata content review
**Deliverable**:
- `manifest.yaml` v0 草稿（commit hash 待填）
- `InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01/` 5 條 finalize

**AI 動作**：
1. Read 現有 `manifest.yaml` + 寫 `v6_binary_commit: <pending>` placeholder + 其餘欄位 fill
2. Read 5 條 errata candidate（含 5/9 paired audit E5）逐條 review
3. 列出每條 errata 是否需 V6 evidence 補強

**🟢 一行告知**: errata content review 是文字工作，可在 ISM 跑同時並行（影響: 中, 信心: 高）

### §1.3 Day 1 收尾 — 用戶 checkpoint
**狀態檢查**:
- [ ] COLO829 ISM 啟動成功（檢查 stderr 無 error）
- [ ] manifest.yaml v0 placeholder ready
- [ ] 5 errata content reviewed + decision: 需補 V6 evidence 的條目標記

---

## §2 Day 2 (週二 5/19) — 證據收斂

### §2.1 COLO829 ISM 完成 + 整合到 step1_master
**Deliverable**: `step1_master_three_way.tsv` 加 COLO829 V6 row（**6/7 done**）

**AI 動作**：
1. 等 Day 1 background COLO829 ISM 完成
2. 加 COLO829 V6 row 到 step1_master TSV
3. Sanity check: ColumnCount, NaN rate, region_count vs 5 樣本一致

### §2.2 6-sample marker coverage 比較（暫不含 H1395 paired pure）
**Deliverable**: `research/v6_bam_tpfp_hp_loh_cn/step6_marker_coverage_6sample.tsv`

**AI 動作**：
1. 跑 marker coverage 比較 script: V3F vs V5 vs V6 across 6 samples
2. 統計：
   - Per-sample Δcoverage (V6-V3F)
   - Sign test direction: n_positive / n_negative
3. 視覺化：boxplot Δcoverage by sample

**🟠 節點暫停 — 階段性產出**:
- 若 V6 marker coverage **退步** 任一樣本 → 🔴 STOP，重新評估 V6 binary（最壞情況 abort tag）
- 若 V6 ≥ V3F 全部 6 樣本 → 推進 Day 3

### §2.3 Day 2 收尾 — 用戶 checkpoint
**狀態檢查**:
- [ ] COLO829 ISM 完成
- [ ] step1_master 6/7 done
- [ ] marker coverage 結果：n_positive / 6
- [ ] **若任一樣本退步 → 用戶 review + abort/continue 決定**

---

## §3 Day 3 (週三 5/20) — Caller F1 + Manifest Finalize

### §3.1 7-sample caller F1 比較
**Deliverable**: `research/v6_bam_tpfp_hp_loh_cn/step7_caller_f1_7sample.tsv`

**AI 動作**：
1. 跑 caller F1 comparison: V3F / V5 / V6 vs SEQC2 truth set
2. F1 metric 口徑明示（依 `feedback_outside_claim_must_query_kb`）:
   - caller-level FILTER=PASS F1（內部標準）
   - hap.py F1（如需對外比較）
3. **Equivalence check**: |F1_V6 - F1_V3F| < 0.005 across all 7 samples
4. 如有 regression sample → STOP

### §3.2 Manifest finalize（填 V6 commit hash）
**Deliverable**: `manifest.yaml` v1 ready for tag

**AI 動作**：
1. 確認 V6 binary 在 HCC1395 + 6 樣本 ISM 都用同一 commit hash
2. `git log -1 --format=%H` on V6 binary branch
3. 寫入 `v6_binary_commit: <full_sha>` + `v6_baseline_v5_commit: 938f0df`
4. 加 sample list + ISM run timestamps

### §3.3 Day 3 收尾 — 用戶 checkpoint
**狀態檢查**:
- [ ] 7-sample caller F1 結果：no regression
- [ ] manifest.yaml v1 ready
- [ ] **🔴 Hard Gate prep**: Day 4 git tag 前需用戶確認 manifest 內容

---

## §4 Day 4 (週四 5/21) — 🔴 Hard Gate 1 + Email Draft

### §4.1 🔴 Hard Gate 1: `git tag v6-prod-{YYYYMMDD}`
**Deliverable**: Git tag created + pushed to remote

**用戶必須親自確認**:
1. AI 顯示 manifest.yaml 完整內容
2. AI 顯示 git tag 命令 + 預期影響（不可逆）
3. **用戶明示「執行 git tag」**才動作

**AI 動作**（用戶 ack 後）:
```bash
# 提示用戶執行（AI 不直接執行 git tag — Hard Gate 不可逆）
git tag -a v6-prod-20260521 -m "V6 production binary..."
git push origin v6-prod-20260521
```

### §4.2 PI Errata Email Draft
**Deliverable**: Email draft 在 `InterSubMod/docs/reports/validated/2026/05/20260521_PI_V6_signoff_email_draft_01.md`

**AI 動作**：
1. 寫 email draft 結構：
   - Subject: V6 Production Binary Sign-Off + 4-29 PI Report Errata (5 items)
   - §1 Summary（V6 binary 動機 + 通過 7 樣本驗證）
   - §2 Errata 5 條（逐條：原文 → 更正 → 影響）
   - §3 V6 sign-off ask（lab meeting 或 written ack）
2. 引用 v6-prod tag + manifest commit hash
3. 估計 reading time 標明（PI 友善）

### §4.3 Day 4 收尾
**狀態檢查**:
- [ ] git tag v6-prod-20260521 created
- [ ] PI email draft ready for user review

---

## §5 Day 5 (週五 5/22) — 🔴 Hard Gate 2 + Buffer

### §5.1 🔴 Hard Gate 2: PI Email Send
**Deliverable**: Email sent to PI

**用戶必須親自確認**:
1. AI 顯示 email draft 完整內容
2. 用戶 review + edit（若需要）
3. **用戶明示「Send」**才執行
4. **AI 不能 send email**（無 SMTP access）— 用戶 copy-paste 到 mail client send

### §5.2 W3 Deadline Buffer
若 Day 1-4 有 1 天 slip → 用 Day 5 補。

### §5.3 完成後（T1.2 ✅）
**Unlock**:
- thread_d_paper Tier 2 Archive TO 7-sample rerun（V6 binary）可啟動
- T4.3 PI errata package 完成（5 條 + V6 sign-off）

**AI 動作**（T1.2 ✅ 後）：
1. Update `docs/CURRENT_FOCUS.md §2026-05-17`: T1.2 ✅ DONE 2026-05-22
2. Update `InterSubMod/knowledge/10_research_status/01_current_focus_snapshot.md`: T1.2 mark done
3. 寫 T1.2 完成 session report

---

## §6 風險 + Fallback

| 風險 | 機率 | 緩解 |
|------|------|------|
| **R1**: COLO829 V6 ISM 跑 > 4hr | 中 | Day 1 啟動後若 4hr 未完成 → 評估是否跳過 COLO829（accept limitation）|
| **R2**: V6 任一樣本 marker coverage 退步 | 低（5/7 已驗 +9% over V3F）| Day 2 stop + 評估 V6 binary 是否 revert |
| **R3**: Caller F1 regression > 0.005 | 低 | Day 3 stop + 找 root cause（可能是 SEQC2 truth set 對 V6 region 標記差異）|
| **R4**: PI 在 sign-off 前要求新 audit | 中 | T4.3 workload 擴大；Day 5 暫不 send，等 PI 提問澄清後 Day 5+ 再 send |
| **R5**: Day 1+2 並行任務記憶交叉污染 | 低 | 兩 track 不共寫同檔案；track A 寫 output/ , track B 寫 docs/ |

---

## §7 AI vs 用戶責任分工

| 任務 | AI 能做 | 用戶必須 |
|------|---------|---------|
| COLO829 V6 ISM 跑 | ✅ 啟動 + 監看 | — |
| Marker coverage / F1 分析 | ✅ script + 表 + 圖 | — |
| Manifest.yaml 寫入 | ✅ | — |
| Errata content review + finalize | ✅ draft + suggest | ack content |
| **🔴 git tag** | ❌ 提示命令 | **執行 `git tag` + `git push`** |
| **🔴 PI email send** | ❌ draft email | **review + send via mail client** |
| CURRENT_FOCUS update（T1.2 ✅）| ✅ | — |

---

## §8 每日標準收尾格式

每天收尾 AI 必須輸出（依 `/scientific-rigor §2.1` Checklist + AGENTS.md §6 Step→Verify）：

```markdown
## Day N (YYYY-MM-DD) — 完成狀況

### Deliverables
- [x] D1: ...（artifact path）
- [x] D2: ...（artifact path）
- [ ] D3: ...（pending, blocked by X）

### Step→Verify
- Step 1: <動作>
  - Verify: <command + 期望輸出>
  - 實際: <實際輸出>
  - 結果: ✅/❌

### Tomorrow's Critical Path
- 最重要的 1 件事
- Blocked? Yes/No

### Hard Gate 距離
- Day 4 git tag: <剩 N 天 / on track / at risk>
- Day 5 email: <剩 N 天 / on track / at risk>
```

---

## §9 相關

- Parent plan: `~/.claude/plans/tender-pondering-blossom.md`
- 00_PLAN: `InterSubMod/research/selfphasing_v6_production/00_PLAN.md`
- Self-phasing 整合: `InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md`
- PI Report 4-29 Errata companion: `InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01/`
- 2026-05-15 V6 evidence: `InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md`
- V6 framework dir: `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/`
