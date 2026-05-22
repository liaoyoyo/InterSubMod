# Fishbone / Ishikawa Skeleton

> 石川馨 (Kaoru Ishikawa) 1968 — 6M 分類

```markdown
---
framework: Fishbone / Ishikawa
when: <複雜根因 / 多因素並存 / 工程 RCA / product defect>
---

# <Effect / Problem>

## 6M 分類因素

### Man（人因）

- <人因素 1>
- <人因素 2>

### Machine（機器 / tool）

- <機器因素 1>
- <機器因素 2>

### Method（方法）

- <方法因素 1>
- <方法因素 2>

### Material（資料 / input）

- <資料因素 1>
- <資料因素 2>

### Measurement（測量）

- <測量因素 1>
- <測量因素 2>

### Milieu / Environment（環境）

- <環境因素 1>
- <環境因素 2>

## 影響權重排序

| Category | 主因子 | 影響度 | 可修正性 |
|----------|--------|-------|--------|
| Method | <X> | 高 | 高 → 優先 |
| Measurement | <Y> | 高 | 高 → 優先 |
| Material | <Z> | 中 | 低 (外因) |

## Action（針對 高×高 因子）

<行動列表>

---

範例（chr19 HP tag bias）:

- **Method**: priority bug 修正只動 priority 不動 assignment
- **Measurement**: chr-level ratio 太粗（per-region 才能抓 1.77×）
- **Material**: chr19 q-arm LOH 區 germline het 密度低（外因）
- **Man**: 早期 PI report 沒要求 per-region 拆解
- Machine / Milieu: 非主因

→ 高×高 = Method + Measurement → 兩個必同時修

---

Framework: Fishbone / Ishikawa (石川馨《Guide to Quality Control》1976)
```
