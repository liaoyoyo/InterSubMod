<!--
建立時間: 2026-08-17
狀態: validated
目標: 讓任何人 clone 後能實際跑完一次 inter_sub_mod，看到輸入與輸出的形狀
處理範圍: 完全合成資料，無生物學意義
build_branch: chore/handoff-20260813
驗證方式: 實跑 RUN.sh，確認 exit 0 且 CpG 數 > 0
-->

# 合成 fixture

**這組資料完全是合成的，沒有任何生物學意義。** 不要拿它的數字做科學解讀。

它存在的理由只有一個：讓人**實際跑完一次**，確認自己的建置正確，並看到輸入與輸出
長什麼樣。公開 repo 沒有、也不應該有真實的 tumor BAM／reference／somatic VCF，
所以在這之前任何人 clone 下來都只能編譯和跑單元測試，永遠不知道真正跑起來是什麼樣子。

## 怎麼用

```bash
# 1. 產生（需要 pysam、samtools、bgzip、tabix）
python3 scripts/make_synthetic_fixture.py

# 2. 跑
bash tests/fixtures/synthetic/RUN.sh
```

檔案不進版本控制 —— 產生器 `scripts/make_synthetic_fixture.py` 才是唯一真值來源，
用固定亂數種子（`SEED = 20260817`），每次產生的結果完全相同。

## 2026-08-17 的實測結果

| 項目 | 值 |
|---|---|
| exit code | 0 |
| Total regions | 2 |
| Total CpG sites found | 244 |
| Average reads per region | 24 |
| 全部輸入檔合計 | 約 14 KB |

輸出目錄會包含 `methylation.csv`（read × CpG 的 beta 值矩陣）、`significance_summary.csv`、
`subclone_structure.txt`、`region_stratification*.tsv`，以及 `run_params.json` 與
`run_summary.json`（完整參數與執行摘要）。欄位意義見
`InterSubMod/.claude/rules/output-structure.md`。

## 🔴 這組 fixture 順便抓出來的一個坑

第一版產生器把 MM 標籤寫成 `C+m,...`，結果是：

```
exit code 0
Total CpG sites found: 0
```

**程式正常結束、不報錯，但甲基化資料全部消失。**

根因：`src/core/MethylationParser.cpp` 比對的目標字串是字面的 `"C+m?"`（含問號）。
SAM 規格中 `?` 表示「未列出的鹼基修飾狀態未知」，`.` 表示「確定未修飾」——
InterSubMod 目前只認 `C+m?` 這個 flavor。

**用自己資料的人若看到 `Total CpG sites found: 0`，第一件要查的就是這個：**

```bash
samtools view your.bam | head -1 | tr '\t' '\n' | grep '^MM'
# 需要看到 MM:Z:C+m?,...   若是 C+m. 或裸 C+m，甲基化會靜默地全部讀不到
```

這類「跑完不報錯但結果是空的」失敗模式，比直接報錯難除錯得多，因為它看起來像是
資料本身沒有訊號。這也正是這組 fixture 值得存在的理由 —— 有一個已知會成功的
基準，才能區分「我的資料有問題」與「我的環境有問題」。
