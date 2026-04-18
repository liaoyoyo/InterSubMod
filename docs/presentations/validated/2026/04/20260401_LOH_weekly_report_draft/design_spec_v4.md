<!--
建立時間: 2026-04-02
目標: PPTX v4 設計規範 — 業界最佳實踐 + 偵探故事線
-->

# PPTX v4 設計規範

## 1. 核心原則

- **One Slide, One Message**：每頁只傳達一個核心訊息
- **視覺優先**：圖片/圖表佔 60-70% 面積，文字 30-40%
- **最少文字**：每頁 ≤ 6 行，每行 ≤ 15 個中文字
- **放大展示**：關鍵數字 40-60pt，讓後排也能看清

## 2. 字體規範

| 元素 | 最小字體 | 建議字體 | 備註 |
|------|:-------:|:-------:|------|
| 大標題 | 28pt | 32-36pt | 中文用 Microsoft JhengHei Bold |
| 副標題 | 18pt | 20-22pt | |
| 內文 bullets | 16pt | 18-20pt | 最多 5-6 行 |
| 強調數字 | 36pt | 40-60pt | accent 色 |
| 表格文字 | 12pt | 13-14pt | |
| 圖片標題/caption | 10pt | 11-12pt | |
| 頁碼/footer | 9pt | 10pt | |
| **絕對最小值** | **10pt** | — | 任何元素不得低於 10pt |

## 3. 配色（深色主題）

- Background: #1a1a2e
- Title: #FFFFFF
- Body text: #E0E0E0（對比度 > 7:1）
- Accent cyan: #00d4ff（數據、正面）
- Accent red: #ff6b6b（警告、否決）
- Accent teal: #4ecdc4（結論、教訓）
- Accent gold: #ffd93d（轉折、發現）
- Card fill: #22224a
- Subtle text: #8888AA

## 4. 佈局常數（16:9, 13.333" x 7.5"）

- 安全邊距：左右 0.6"、上 0.4"
- 標題區：y=0.4, h=0.85"
- 內容區：y=1.4 ~ 6.2"（高度 4.8"）
- Footer 區：y=7.05"
- 留白目標：≥ 20% 面積

## 5. 新版 Slide Types（v4 新增）

### 5.1 section_divider — 段落分隔頁
- 全版深色背景 + 居中大標題 (36pt)
- 副標題 (18pt) + 階段進度條
- 用途：三幕之間的呼吸空間

### 5.2 big_number — 大數字強調頁
- 一個超大數字 (60pt) + 一行說明 (20pt)
- 可選：左側數字 + 右側圖片
- 用途：關鍵統計數字的視覺衝擊

### 5.3 full_image — 全版圖片頁
- 圖片佔 80%+ 面積
- 標題在上方半透明條上
- 用途：最重要的數據圖表

### 5.4 process_flow — 流程圖頁
- 水平或垂直的步驟流程
- 圓角矩形 + 箭頭連接
- 用途：ISM pipeline、因果鏈

### 5.5 before_after — 前後對比頁
- 左右分割，中間分隔線 + 箭頭
- 用途：HP bug fix、QS 修正

### 5.6 funnel — 漏斗/篩選頁
- 遞減寬度的橫條 + strike-through
- 用途：假說否決、驗證架構

## 6. 投影片總數目標

v3: 16 頁 → v4: **25 頁**（含 3 張 section divider）

## 7. 故事結構（三幕）

### Act I: 發現異常（Slides 1-8）
- Cover + Agenda
- 動機（為何研究 TO）
- 技術背景（ISM pipeline + Paired vs TO）
- HP Bug 發現與修正
- LOH enrichment 方向翻轉

### Act II: 深入追查（Slides 9-17）
- 288K 同位點配對（大數字 + 全版圖）
- 系統性觀察總覽
- 排除 depth 假說
- HP Ratio 直接證據（全版圖）
- Copy Number 排除
- 根因確認：Self-Phasing（機制圖 + 量化）

### Act III: 修正與展望（Slides 18-25）
- 甲基化假說全部否決（funnel）
- QS mode-aware 修正（before/after）
- 修正方向路線圖
- 核心結論（對比總結）
- 討論與下週方向
