# Example: Debug postmortem — A3 + 5 Whys + Inverted Pyramid + Blameless 5-段

> 場景: bug fix postmortem, Who=工程同儕, Why=紀錄根因, What=失敗事件, When=紙本

```markdown
---
framework: A3 Report + 5 Whys + Inverted Pyramid (lead) + SRE Blameless Postmortem
audience: 工程同儕 / future me
when: 紙本詳述
applied_by: narrative-frame v0.1
cycle_id: 20260522-XXXX
---

# A3 — chr19 HP tag priority bug 17.3× amplification 根因解析

## TL;DR (Inverted Pyramid lead, BLUF)

V3F V1.7-I 修正後仍見 17.3× HP1 偏置 — 經追查為 **assignment 1.77× × priority bug 9.8×** 疊加，而非單一 priority bug 17.5×。修正路徑: V6 同時處理 assignment + priority 兩層。

## 1. Background

ONT BAM phasing 由 longphase 內部 HP tag 決定 caller direction；priority bug 會讓 chr19 全 read 標錯 HP。

## 2. Current State (before fix)

- chr19 752/752 reads hp=33（V3F 保守處理）
- V5 Layer 1.5 在 germline-absent 區仍 4.19:1 偏 HP1
- HCC1395 paired-pileup F1=0.7153（卡住）

## 3. Goal

- chr19 victim 100% rescue（hp ∈ {1, 2}）
- HP1:HP2 ratio 在 germline-absent 區 ≤ 1.5:1（合理 phasing baseline）

## 4. Analysis（根因）

### 4.1 5 Whys

- Why 1: V1.7-I 修正後仍見 17.3× HP1 → because 既有 model 假設 priority bug 是唯一 amplifier
- Why 2: 為何 priority bug 單獨不夠 → because 上游 assignment 已 1.77× 偏 HP1
- Why 3: 為何 assignment 偏 HP1 → because longphase getVote() countMap 在 per-read 內 reset 但 cumulative across reads 偏 HP1
- Why 4: 為何 cumulative 偏 HP1 → because germline het positions 在 chr19 q-arm 密度遞減 (LOH 區)
- Why 5: 為何沒早發現 → because **既有 metric 只看 chr-level ratio 不看 per-region 分佈**

**Root cause**: Step 5 — metric 設計盲點導致 1.77× 上游偏置被忽略 12+ 月

### 4.2 Fishbone（6M）

- **Method**: priority bug 修正只動 priority 不動 assignment
- **Measurement**: chr-level ratio 太粗（per-region 才能抓 1.77×）
- **Material**: chr19 q-arm LOH 區 germline het 密度低（外因）
- **Man**: 早期 PI report 沒要求 per-region 拆解
- Machine / Milieu: 非主因

## 5. Proposal

- Option A: 只修 priority bug（V5 路徑）→ 17.3× 無法解（不完整）
- Option B: 修 priority + assignment（V6 路徑）✅ chosen
- Option C: 改用 PON-only phasing → 副作用大（LOH detection 損失）

## 6. Plan

| Step | Who | When | Verify |
|------|-----|------|--------|
| 修 V6 binary | 工程 | 2026-05-09 | unit test pass |
| Run V6 paired pileup HCC1395 | 跑 pipeline | 2026-05-12 | chr19 victim rescue ✓ |
| Cross-sample LOSO | 跑 analysis | 2026-05-20 | F1 uplift 量化 |

## 7. Follow-up

- 加 per-region HP ratio 進 monitoring dashboard（避免再盲）
- audit 過去 12 月所有 V3F-based 結論

---

## §SRE Blameless 5-段（補根因敘事）

### Timeline

- 2026-Q1: chr19 priority bug 首次發現（V3F）
- 2026-Q1-Q2: V1.7-A → V1.7-I 多次 patch
- 2026-05-09: V6 落地（修 assignment + priority 兩層）
- 2026-05-12: V6 paired pileup HCC1395 verified
- 2026-05-20: LOSO H_NEW_4 confirmed +0.00699 HCC1395

### What Went Well

- Priority bug 機制本身已正確診斷（V3F 已部分修）
- 17.3× 量化指標讓問題量化 — 不會被掩蓋

### What Went Poorly

- 12+ 月才發現 assignment 1.77× 上游偏置 — metric 盲點
- 5 個 V1.7 patch 都假設 priority bug 是唯一 amplifier — confirmation bias

### What We Got Lucky With

- chr19 hotspot 夠 visible（752 reads 集中）— 若分散在 N chr 可能更難發現
- 用戶（PI）早期要求 per-region 視覺化才促成此 audit

### Action Items

- [x] V6 落地修兩層
- [ ] 加 per-region HP ratio dashboard（防止再盲）
- [ ] Audit 12 月 V3F 結論（是否其他樣本也受影響）

---

Framework: A3 Report (Toyota TPS; Shook 2008) + 5 Whys (Ohno) + Inverted Pyramid (journalism) + SRE Blameless Postmortem (Google)
```

---

## 為什麼選此 hybrid（N2 推薦理由）

**5W 識別**:
- Who: 工程同儕 + future me
- Why: 紀錄 + 教訓
- What: 失敗事件 + 根因
- When: 紙本詳述（無時長壓力）
- How: .md 報告

**主框架 A3** — Toyota 工程改善敘事完整覆蓋背景 / 根因 / 提案 / 計畫
**Sub: 5 Whys** — §4.1 用，追到 root cause
**Sub: Fishbone 6M** — §4.2 用，補多因素分類
**Lead: Inverted Pyramid + BLUF** — TL;DR 給沒時間讀全文的人
**Wrap: SRE Blameless** — Postmortem 不歸咎個人，著重系統性失敗

**不選擇**:
- SCQA: 缺 timeline + action items 結構
- AIDA: 行銷框架不適合根因敘事
- Hero's Journey: 個人化敘事 vs blameless 衝突
