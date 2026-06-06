# 方法解釋圖 verify checklist（12-criteria + 圖例細節迭代 + SciFig + 5 秒）

> 三層疊用：**(A) 12-criteria 領域 rubric（主）** + **(B) SciFig R1-R6 通用 rubric** + **(C) 5 秒測試**。
> 每輪 FIX 後重跑 renderer（diff SVG）→ 對「圖例細節」段逐項複查 → 就地修 data/spec 再重生（**持續檢視→驗證→修正**）。

## A. 12-criteria（5 維 / 逐項判 PASS/PARTIAL/FAIL）

### 生物忠實 BIO
- **BIO-1** read→CpG→haplotype 三層正確；HP1-1 = HP1 germline 下 somatic-allele 子標（**非第三 haplotype**）；**label flip：HP1-1（分析側）≡ HP2-1（圖側）須一致或註明**。FAIL：CpG 當 read 屬性 / HP1-1 當第三 haplotype / 甲基方向標反 / label-flip 漏標。
- **BIO-2** β＝某群在某 CpG 的甲基比例（0-1，per-CpG 群率）；Δβ 有方向（負=tumor 低甲基）。數字須與 verified 真值一致。

### 工程嚴謹 ENG
- **ENG-1** Δβ 五步齊全且順序對：(1)per-CpG 兩群算 β (2)取 shared CpG (3)per-CpG Δ (4)**mean=Δβ 主鍵 + max|Δ| 副鍵齊列** (5)paired Wilcoxon over CpG。FAIL：漏 shared-CpG / Wilcoxon 講成對 read / 漏副鍵。
- **ENG-2** 舊 ISM 描述正確：10kb + k=4 分群 + Cramér's V(hp/alt) + gating；**須含 V alt 與 n_cpg（非只 V hp + n_reads）**；指出 germline_hp_only 丟 somatic 子標後果。
- **ENG-3** 每數字可追溯（L1，附 src），無新編造；衍生/合成數字明標。

### 清晰 CLR
- **CLR-1** 一圖一概念、先直覺後公式、過 5 秒測試。
- **CLR-2** 新舊用對照結構呈**互補**（非取代）。

### 誠實 HON
- **HON-1** 示意 vs 真值明標且視覺可區分（示意=synthetic watermark / 文字註記；真值=附 src）。
- **HON-2** 不過度宣稱（如 BRCA2 真 cis 小 d_within marginal、單樣本 ★3）；不可隱 d_within 只報 d_cis。
- **HON-3** 誠實標方法盲點（如 Δβ-mean 交叉型盲點 + 副鍵救援）。

### 完整 CMP
- **CMP-1** 各段皆有實質且自洽。
- **CMP-2** 同一 anchor 位點貫穿；合成例與真實例明確區分。

## B. SciFig R1-R6（通用 method-figure rubric，L3 借用當 pre-publish checklist）
- **R1 技術正確性**（含 label-flip 與數字核對）· **R2 視覺清晰可讀** · **R3 結構連貫** · **R4 設計一致** · **R5 可詮釋 + accessibility**（色盲友善 + 形狀冗餘 + caption 可獨立讀）· **R6 實作品質**（viewBox/role/title/desc、無 chartjunk）。

## C. 5 秒測試
第一眼能否抓到主訊息？標題給 thesis、對照並列、色彩+形狀冗餘。扣分：單圖塞太多子圖（>5 秒才跟完）。

---

## 🔁 圖例細節迭代 verify loop（持續檢視→驗證→修正）

每輪改完重跑 `render_figure_spec.py`，對下列**圖例/legend 細節**逐項複查（這些最容易在迭代中 drift）：

- [ ] **label-flip 註明**（HP1-1≡HP2-1）— 漏 = BIO-1 FAIL（領域捏造高風險）
- [ ] **HP1-1 白話定義**（「HP1 germline 下 somatic-allele 子標，非第三 haplotype」）出現一次
- [ ] **色標 + 形狀冗餘一致**（實心紅=甲基 / 空心藍框=未甲基；HP1 綠 / HP1-1 橙紅 / normal 灰）跨所有子圖一致
- [ ] **單位/座標標清楚**（β 0-1、CpG 軸、SNV 標）
- [ ] **示意標註**（schematic reads = synthetic watermark + 文字註記）
- [ ] **真值附 src**（provenance 稽核表每 metric 有來源；HTML 嵌入時 frontmatter `data_sources:` 齊）
- [ ] **Δβ 符號方向不畫反**（負值 → tumor 低甲基；bar 長度/方向一致）
- [ ] **anchor 位點跨圖一致**（chr13:32,315,128 在每張子圖同一位點）
- [ ] **口徑不混**（cis-test β 0.25/0.228/0.108 ≠ deepdive 5mC-only 0.627/0.592；同圖只用同口徑）
- [ ] **caption 可獨立讀懂**（不看正文也懂這張圖在說什麼）

**修法分級**：圖例層問題（標籤/色/註記）→ 改 `data.json`(label_flip_note / schematic) 或 `figure_spec.json`(shared.color_scale / annotations) 再重生，**不必重畫**。結構層問題（缺物件/佈局）→ 加/調 primitive。數字問題 → 回 LOCK-AND-GATHER 補 verified 真值（**勿手打**）。

## 反捏造硬底線（與 §13 對齊）
- 數字 grep 不到來源 = 捏造 → renderer 已 refuse（verified 缺值 exit 3）。
- 示意當量化（read 甲基點當真資料）= 非數字型捏造 → 必 synthetic watermark + 文字註記。
- 「不在白名單」≠「刪已驗證真值」（§13 要 traceable 不是 delete）。
- 轉 validated/PI 前：frontmatter `data_sources:` 齊（否則 `number_provenance_check.sh` exit 2）。
