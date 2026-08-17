<!--
建立時間：2026-08-13
目標：提供可公開、可重生的 InterSubMod tiny synthetic E2E fixture
處理範圍：1 synthetic contig、1 SNV、12 synthetic reads；DEMO only
關聯檔案：InterSubMod/scripts/handoff/build_tiny_public_fixture.sh；InterSubMod/scripts/handoff/run_tiny_public_e2e.sh
-->

# Tiny public synthetic fixture — DEMO ONLY

> **PARTIAL / DEMO：此 fixture 只驗證軟體執行、I/O contract 與 199-column schema；不可作科學 validation evidence。**

- `reference.fa`：200 bp synthetic contig `chrTiny`。
- `variants.vcf`：單一 `chrTiny:101 A>T` synthetic SNV。
- `tumor.sam`：12 條 synthetic reads；6 條 HP1/REF、6 條 HP2/ALT，全部具有 Dorado-style `MM:Z:C+m?` 與 `ML:B:C` tags。
- BAM、BAI 與 FAI 由 `InterSubMod/scripts/handoff/build_tiny_public_fixture.sh` 在隔離工作目錄重生，Regular Git 不保存 binary。

`clustering/tree.nwk` 是由 read–read methylation distance 產生的 **read dendrogram**。它不是 cellular lineage、clonal ancestry、phylogeny truth、subclone prevalence 或生物正確率的證據。

本 fixture 服務 G4（reproducibility）與 G5（外部可驗證交接），不驗證 G1/G3 的生物學 claim。

---

## 🔴 `MM:Z:C+m?` 的問號是必要的 —— 換自己的資料時最容易踩的坑

上面第三點寫的 `MM:Z:C+m?` 那個 `?` 不是可有可無的裝飾。
`src/core/MethylationParser.cpp` 比對的目標字串是**字面的 `"C+m?"`**。

SAM 規格中 `?` 表示「未列出的鹼基修飾狀態未知」，`.` 表示「確定未修飾」。
**InterSubMod 目前只認 `C+m?` 這個 flavor。** 若你的 BAM 寫成 `C+m.` 或裸 `C+m`：

```
exit code 0
Total CpG sites found: 0
```

程式正常結束、不報錯、甲基化資料全部消失。這種失敗比直接報錯難除錯得多，
因為它看起來像是資料本身沒有訊號。

**看到 `Total CpG sites found: 0` 時的第一個檢查：**

```bash
samtools view your.bam | head -1 | tr '\t' '\n' | grep '^MM'
# 需要看到 MM:Z:C+m?,...
```

這也是這組 fixture 存在的價值之一 —— 有一個已知會成功的基準，才能區分
「我的資料沒有訊號」與「我的環境或輸入格式有問題」。

（本註記由 2026-08-17 的獨立重現實驗補上：另做一組合成 fixture 時第一版誤用
`C+m,`，得到 exit 0 但 CpG 數為 0，因而定位到這個 parser 限制。那組重複的
fixture 已移除，本 fixture 為唯一入口。）
