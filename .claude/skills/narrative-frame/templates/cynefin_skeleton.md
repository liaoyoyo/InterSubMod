# Cynefin Skeleton — 5 域分類

> Dave Snowden, IBM Research (1999); Snowden & Boone HBR 2007

```markdown
---
framework: Cynefin Framework
when: <因果性質判斷 / 危機分類 / 組織決策>
---

# Domain Gate

## 5 Domains

### 1. Clear（清晰）
- **因果**: 明顯 + 重複 + universal
- **Practice**: Best practice 直接套
- **動作**: Sense → Categorize → Respond
- **範例**: 標準 SOP；已知 ATM 操作

### 2. Complicated（複雜）
- **因果**: 需要 expertise 但可確認
- **Practice**: Good practice + expert judgment
- **動作**: Sense → Analyze → Respond
- **範例**: bug debugging；汽車修理

### 3. Complex（複雜系統） ⚠
- **因果**: emergent，事前不可預測，事後可解釋
- **Practice**: **probe-first**（小規模試 + 觀察）
- **動作**: Probe → Sense → Respond
- **範例**: 新研究方向；組織文化變革；新 framework cross-sample 驗證
- **⚠ 禁套 best-practice — deterministic framework 會誤導**

### 4. Chaotic（混亂）
- **因果**: 無；需要先穩定
- **Practice**: Act first → 之後再分析
- **動作**: Act → Sense → Respond
- **範例**: 火警；data center down

### 5. Disorder（無法分類）
- 不知道在哪個 domain → 先做 domain classification

## Domain 判斷問句

**「相同行動是否曾重複產生可預測結果？」**
- Yes（已知 cause-effect）→ Clear / Complicated
- Maybe（partial pattern）→ Complex
- No（混亂）→ Chaotic

---

範例（research direction 決策）:

「Cycle 5 chr8 zone gate 屬於哪個 domain？」

- 相同 chr8 zone rule 是否曾 cross-sample 一致？→ NO（只有 HCC1395 single-sample）
- **Domain = Complex** → 強制 probe-first
- 不可套 deterministic LR threshold framework
- 必須跑小規模 pilot (HCC1937 + HCC1954) 觀察 → 再決定全推

---

業界引用: 「Match practice to causality nature — not the other way around」(Snowden)

---

Framework: Cynefin (Snowden, IBM 1999; Snowden & Boone HBR 2007)
```
