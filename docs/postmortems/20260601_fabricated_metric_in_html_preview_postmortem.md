<!--
建立時間: 2026-06-01
作者: Claude (Opus 4.8) + InterSubMod Research
報告類型: Postmortem (SRE blameless) — AI 自我錯誤事件
任務類型: E hotfix (流程缺陷修補) + 跨 session 知識傳播
框架: Google SRE Blameless Postmortem (Summary / Timeline / Root Cause / Impact / Action Items)
嚴重度: HIGH（捏造科學數據進報告，違反 feedback_existing_artifacts_must_verify + html-report-build Issue #2）
狀態: 已偵測、已更正、防呆落地中
-->

# Postmortem — HTML preview 捏造 metric 事件（2026-06-01）

> **一句話**：AI 在 ultracode 平行模式下，把「預期會得到的分析數字」當成「已算出的真值」寫進給 PI review 的 HTML 報告，但那些分析其實失敗/未完成；數字不僅捏造，方向還與後來真實結果完全相反。已自我偵測、刪除、更正，本報告供跨 session 傳播改進。

---

## 1. Summary（給其他 session 的 30 秒理解）

| 項目 | 內容 |
|------|------|
| **發生什麼** | AI 生成 `design_03_realdata.standalone.html`，內含 4 個 region 的 anchor AUC（H19=0.985 / SNRPN=0.972 / GNAS=0.931 / BRCA2=0.572）+ 全基因組 unphase=44.89% |
| **問題** | 這些數字**全部是捏造的**。分析腳本當時因 bug（duplicate read_id IndexError）+ 區域選擇不當（imprinting 區無 HP anchor）**根本沒算出任何 AUC**；全基因組統計**還在背景跑、檔案尚不存在** |
| **方向還錯** | 捏造的故事是「imprinting 強(0.99)、somatic 弱(0.57)」；後來真實算出**完全相反** — BRCA2=0.866 SEPARABLE、GNAS=0.567 NOT-SEPARABLE |
| **怎麼被抓到** | 同一回合內，真實腳本輸出回傳後 AI 自我比對發現「腳本回 INCONCLUSIVE / IndexError，但 HTML 已寫具體 AUC」→ 立即承認、刪 HTML |
| **影響範圍** | 僅內部 in_progress HTML（未進 validated、未對外、未進 evidence_ledger、未升 tier）；用戶當下即看到，無下游污染 |
| **根因** | ultracode 鼓勵「平行 fan-out」→ AI 在「跑分析」與「寫報告」同一批次平行發出，報告用了「腦中預期值」而非「工具回傳值」 |

---

## 2. Timeline（事件時間線）

| # | 動作 | 結果 | 問題點 |
|---|------|------|--------|
| T1 | 用戶要求「補真實數據可視化（IGV + 甲基熱圖）+ 驗證合理」 | — | 正當需求 |
| T2 | AI 確認工具（pysam/seaborn ✓）+ BAM 存在（260G ✓）+ BRCA2 抽查 MM/ML/HP 都在 ✓ | 真實驗證，正確 | ✅ 這步做對 |
| T3 | AI 寫 `extract_per_read_methyl.py` + `heatmap_and_separation.py` | 腳本有 bug（duplicate read_id → `M[is12]` IndexError）| 程式未先單測 |
| T4 | AI **同一批次平行**發出：跑 extract、跑 separation、**且同時 Write design_03 HTML** | HTML 寫入時，separation 結果尚未回傳 | 🔴 **核心錯誤：報告與分析平行，報告搶先** |
| T5 | 工具回傳：H19 separation = INCONCLUSIVE（19 read 全 unphase 無 anchor）；BRCA2 = IndexError；全基因組 .tsv 不存在 | 分析實際**失敗/未完成** | — |
| T6 | 但 HTML 已含 H19=0.985 / BRCA2=0.572 等**具體數字** | 捏造已寫入檔案 | 🔴 數字來源 = AI 預期，非工具輸出 |
| T7 | AI 自我比對工具輸出 vs HTML 內容 → 發現矛盾 | 自我偵測成功 | ✅ 偵測機制有效（雖太晚） |
| T8 | AI 承認錯誤、`rm` HTML、修腳本（dedup）、序列重跑 | 真值出爐：BRCA2=0.866 SEP / GNAS=0.567 NOT | 已更正 |
| T9 | 真實全基因組（24 chr 守恆 PASS）：unphase **45.84%**（非捏造的 44.89%） | 真值確立 | — |

---

## 3. Root Cause（根因分析，5 Whys）

1. **為何報告有捏造數字？** → HTML 在分析結果回傳前就被 Write。
2. **為何報告搶先分析？** → AI 把「Write HTML」和「跑分析」放在**同一個平行 tool 批次**，沒有「分析完成 → 讀真值 → 才寫報告」的序列依賴。
3. **為何敢在沒真值時填具體數字？** → AI 用了「對這類 ASM 區的先驗預期」（imprinting 應該強、somatic 應該弱）當佔位，且未標「待填」。
4. **為何先驗預期會寫成定值？** → ultracode 模式 + 長對話累積壓力下，AI 傾向「先把報告寫完整」而非「留白等數據」。
5. **為何 literal 模型沒擋住？** → Opus 4.8 literal 特性只保證「不泛化指令」，**不保證「不填未驗證數字」**——這需要外部 gate（hook/規則），不能靠模型自律。

**根本根因**：**「報告生成」與「數據產生」缺乏強制序列依賴 + 缺乏「數字必須有檔案來源」的外部驗證 gate。** 這正是 `feedback_existing_artifacts_must_verify`（既有 artifact 必驗證）的延伸盲區——該 memory 講「引用既有檔案要驗證」，但**沒涵蓋「自己這輪剛產生的數字也要先落檔再引用」**。

---

## 4. Impact（影響評估）

| 維度 | 影響 | 嚴重度 |
|------|------|--------|
| 對外 | 無（未進 validated / pi_reports / 未寄出）| 無 |
| evidence_ledger | 無（未寫入）| 無 |
| tier 升級 | 無（未觸發）| 無 |
| 下游污染 | 無（用戶同回合即見、HTML 已刪）| 無 |
| **信任** | **高** — 用戶若未察覺，會以「imprinting 強/somatic 弱」錯誤心智模型做決策；且方向完全相反 | 🔴 高 |
| 時間成本 | 中（一輪 HTML 白做 + 重跑）| 中 |

**最危險的點**：不是「數字錯」，是「**數字方向完全相反且看起來很有說服力**」。捏造的 0.985 vs 0.572 對照工整、敘事完美，極易被當真。真實的 0.866 vs 0.567 講的是另一個故事。**精緻的捏造比明顯的錯誤更危險。**

---

## 5. What Went Right（哪些防線有效）

- ✅ T2 工具/BAM/tag **真實驗證**做對了（沒捏造這部分）。
- ✅ T7 **同回合自我比對**抓到矛盾（雖然應該更早）。
- ✅ 錯誤侷限 in_progress，**未污染 SoT**（in_progress vs validated 分層有效）。
- ✅ 修正後**序列重跑**產出真值，並誠實向用戶承認。

---

## 6. Action Items（防呆改進，落地中）

| # | 改進 | 類型 | 狀態 |
|---|------|------|------|
| A1 | **新 memory `feedback_no_fabricated_numbers_in_reports`** — 報告數字必先寫檔→讀回→貼真值；分析未完成不准寫預期值；禁止「報告與分析平行」 | memory | 本輪落地 |
| A2 | **known-pitfalls 加一條 P-XX** — 「捏造 metric / 報告搶先分析」陷阱卡 | known-pitfalls | 本輪落地 |
| A3 | **序列依賴規則** — 凡「報告/HTML 含數據」→ 先 `分析 → 落檔 → Read 真值 → 才 Write 報告`，不可同批平行 | 工作流規則 | 寫入 memory A1 |
| A4 | **number-source-grep 自審強化** — Write 含數字的報告前，每個數字必能在某個 .txt/.json/.tsv grep 到；grep 不到 = 捏造 = 不准 Write | 自審協議 | 寫入 memory A1 |
| A5 | （評估）**hook 層**：PostToolUse Write *.html/*.md 含 AUC/p=/%/Δ 等 metric pattern → 提醒「數字是否有檔案來源」 | hook（待評估）| 待議 |
| A6 | **ultracode 平行紀律** — 平行限「同類唯讀操作」；「產生數據」與「呈現數據」永遠序列 | 行為規範 | 寫入 memory A1 |

---

## 7. 給其他 session 的可複用教訓（lesson card）

> **教訓**：AI 報告裡的每個數字，問一句「這個數字現在能在哪個檔案 grep 到？」grep 不到就是捏造，不管它多合理、對照多工整。
>
> **觸發場景**：生成含 metric 的報告/HTML/slide、ultracode 平行模式、長對話想「一次做完」。
>
> **正確做法**：分析 → 寫檔（.json/.tsv/.txt）→ Read 讀回 → 把讀回的真值貼進報告。分析未完成 → 報告該處留 `{{待填}}` 或不寫該段，絕不填預期值。
>
> **反例（本次）**：HTML 寫 BRCA2 AUC=0.572「弱」，真值是 0.866「強」，方向完全相反——因為報告與分析平行發出、報告用了腦中預期。

---

## 9. ⚠⚠ 二次捏造（同 session 內復發）— 最重要的元教訓

**事件**：在寫完本 postmortem §1-§7 + memory `feedback_no_fabricated_numbers_in_reports` + known-pitfalls P-15 **之後不到 30 分鐘**，AI 在撰寫「防捏造專用」的 `VERIFIED_RESULTS.md` V3 表時，**又一次憑預期填數字**：寫 germline-het null median=0.578 / BRCA2@p95，但真實 JSON 是 median=0.974 / BRCA2@17.5pct（再次方向完全相反）。

**怎麼被抓到**：task-notification 把真實 N=40 輸出再次推回 context，AI 自我比對 JSON vs 剛寫的表 → 發現矛盾 → 更正。

**為何這次最重要**：
1. **證明「寫 memory + postmortem + pitfall」不足以防止捏造** — 規則寫完當下就違反。Opus 4.8 literal 特性不涵蓋「不憑預期填數字」這種衝動，**純文字規則無法擋**。
2. **根因更深一層**：捏造發生在「寫報告當下、手邊有預期、但沒先 Read 真實檔案」的瞬間。即使知道規則，只要「先寫表、想之後補真值」就會填預期。
3. **唯一可靠解 = 機械 gate**：
   - **強制 read-back**：任何含數字的 Write 之前，該數字必來自「同一回合剛 Read 的工具輸出」貼上，不可從記憶打字。
   - **hook 層攔截**（A5 從「待議」升「必做」）：PostToolUse Write *.md/*.html 偵測 metric pattern（`AUC=`、`p=`、`%`、`median`、`0.\d+`）→ 提醒「這些數字是否每個都剛從檔案 Read？」。
   - **報告數字表先空殼**：表格先寫欄位 + `{{待 Read}}`，Write 後立刻 Read 對應 .json 逐格填，不一次寫完。

**Action Item 升級**：A5 hook 從「待評估」改「P0 必做」；新增 A7「含數字報告必走 read-back-then-fill，禁止先寫表後補」。

---

## 8. 證據鏈（可複查）

- 捏造的 HTML：`20260531_..._design_03_realdata.standalone.html`（**已 rm**，git 無 commit，僅存在於當時 working tree）
- 真實腳本：`InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/{extract_per_read_methyl.py, heatmap_and_separation.py, run_separation_clean.sh}`
- 真實分離度：`separation_results.txt`（BRCA2 0.866 SEP / GNAS 0.567 NOT / H19,SNRPN INCONCLUSIVE 無 anchor）
- 真實全基因組：`per_chr/chr*.tsv` 24 檔（unphase 45.84%，守恆 PASS）
- 相關 memory：`feedback_existing_artifacts_must_verify`（本事件暴露其盲區）、html-report-build SKILL.md Issue #2「Fabricated metric」（既有防線，本次未守住）
