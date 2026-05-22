# Freytag's Pyramid Skeleton

> Gustav Freytag《Die Technik des Dramas》(1863)

```markdown
---
framework: Freytag's Pyramid
when: <戲劇 / 小說 / 長篇 case / 紀錄影片>
---

# <Title>

## 1. Exposition（鋪陳）

<人物 + 設定 + 起始狀態 + 隱含張力>

## 2. Rising Action（上升）

<衝突逐步升級 — 累積>

- Event 1
- Event 2
- Event 3 (張力最強)

## 3. Climax（高潮）

<關鍵轉折點 — 主角面對最大挑戰 / 重大決定>

## 4. Falling Action（下降）

<衝突開始解決 — 後果展開>

## 5. Resolution / Dénouement（結局）

<新平衡 + loose ends 收尾>

---

範例（chr19 priority bug 戰役）:

1. Exposition: 2024 ClairS-TO 用 longphase HP tag — 假設一切正常
2. Rising: chr19 752 reads hp=33 異常 → V3F 修一部分 → V5 仍 4.19:1
3. Climax: 12 月 audit — 發現 assignment 1.77× 才是隱藏主因
4. Falling: 設計 V6 修兩層 → HCC1395 LOSO 驗證
5. Resolution: 100% chr19 修復 + 學到 metric 盲點教訓 + 未來 monitoring 加強

---

Framework: Freytag's Pyramid (Freytag《Die Technik des Dramas》1863)
```
