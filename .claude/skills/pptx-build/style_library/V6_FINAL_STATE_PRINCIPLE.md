# V6 為終態對照原則（self-phasing 系列 PPT）

> 2026-05-12 用戶確認。Self-phasing 系列 slide 預設用 **baseline vs V6** 對照，V3F/V5 只在推演必要時使用。

## 核心規則

| 場景 | 對照標的 | 何時用 V3F / V5 |
|---|---|---|
| **預設** PPT slide（slide 06 / 13 / 14 / 15 等指標、機制、結果展示）| **baseline vs V6** | 不用 V3F / V5 |
| **推演必要**（slide 07 / 10 / 16 commits 鏈、Layer 1.5 caveat）| baseline / V3F / V5 / V6 全列 | 顯示 commit 演進過程 |
| **歷史個案**（slide 04a/04b SP1/2/3）| baseline / V5 / V6 全列 | V5 翻 HP2 對齊 paired 是事實鐵證，V6 不能取代 |

## 為什麼

1. **V6 是 production candidate**（current HEAD）— 對 PI 報告就是「我們最終實作什麼」
2. **V3F/V5 是中間階段** — 對 PI 而言不需要每張都看 commit 演進
3. **V5 vs V6 對照（SP1/2/3）是例外** — Layer 1.5 trade-off audit 必要保留 V5 個案

## 應用 checklist

寫 slide 時自問：
- ✅ 這張 slide 目的是「展示最終實作效果」？→ 只列 baseline vs V6
- ✅ 這張 slide 目的是「解釋為何修法分多階段」？→ 列 baseline / V3F / V5 / V6 完整鏈
- ✅ 這張 slide 目的是「揭露 V5 vs V6 trade-off」？→ 列 baseline / V5 / V6
- ❌ 不要因「先寫 V5 後加 V6」就保留 V5 並列 — 預設清理為 V6 only

## 已套用的 slide

| Slide | 對照方式 | 原因 |
|---|---|---|
| 04a/04b SP1/2/3 | baseline / V5 / V6 IGV | V5 個案鐵證必要 + V6 對照 |
| 05 phasing 球員兼裁判 | baseline 概念 + PON-only 修法 | V6 = PON-only flag + tag 修，slide 05 只講 phasing 修 = 8b8c1fd 即可 |
| **06 priority bug** | **baseline vs V6** | 預設原則；V3F/V5 在 speaker note 末尾一句帶過 |
| 07 兩層修補 | baseline / V3F / V5 / V6 commits 完整 | 推演必要（commit 鏈時序展示） |
| 08-15 metrics | baseline vs V6 | 預設原則 |
| 16 V5 caveat | baseline / V5 / V6 三向 | trade-off audit 必要 |
| 17 errata | 視 errata 條目 | 依各條 errata 講解需要 |

## 注意

V5 / V6 在 caller F1 對 SEQC2 truth set **三版本完全相同**（HCC1395 0.93 = 0.7166, 0.6 = 0.6273）— 不是 caller F1 變化驅動 V6，是 marker engineering + 4 樣本 ratio + germline-absent 區行為驅動。

V3F → V5 → V6 是「演算法精化過程」而非「F1 改善過程」。
