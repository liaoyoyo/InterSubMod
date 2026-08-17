<!--
建立時間: 2026-08-13
目標: 記錄 CCU lab-tutorial OLD-P0 與 REGRESSED findings 的本機來源修正、generated-site 同步與驗證
處理範圍: /tmp/lab_tutorial_p0.GIEknc/repo（baseline 9eb1618d359e602d9c528675952b20d051fb2346）
關聯檔案:
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/lab_tutorial_delta_reaudit.tsv
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json
  - InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md
  - InterSubMod/research/20260813_public_docs_p0_correction/00_INDEX.md
-->

# CCU lab-tutorial P0／REGRESSED correction receipt

## 結論

Task type **B — Comprehensive validation continuation**；服務 G2／G3／G4／G5。
本子任務指定的 13 個 findings 已全部完成本機 source 與 generated site correction：

- OLD-P0：INDEX、PRINTALL、M09、M12、M13、SR1、SR2B，共 7 項。
- status=REGRESSED：DELTA-NEW-001、002、004、005、006、008，共 6 項。
- Finding-by-finding positive／negative assertion：**13/13 PASS**。
- 全站 build、file-safe verify、SVG layout、DOM／內部連結、SVG XML 與 git diff --check：皆 exit 0。

LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED 只表示 clean clone working tree 已修正，
**不代表 GitHub default branch 或 live Pages 已更新**。本輪沒有 commit、push、merge、PR 或 deploy。
OLD-P1 與非 REGRESSED delta findings 不在本子任務的完成宣告內，仍須由總稽核追蹤。

## Step → Verify

1. 鎖定 baseline → 驗證：HEAD=9eb1618d359e602d9c528675952b20d051fb2346，開始修改前 worktree clean。
2. 依 authority 修正文與公式 → 驗證：read／HP family 不再升格成 cell、親源或實體 copy。
3. 依 denominator registry 修正數字 → 驗證：M1、formal methylation、nested funnel 與 0.909 語義可逐項重算。
4. 修正 candidate model 語氣 → 驗證：likelihood、error floor、runtime、triplet 與 local topology 均列明假設及校準需求。
5. 同步 figure／caption／quiz／glossary → 驗證：不是只加警語；正文、式子、SVG 與測驗均同步。
6. 重建 generated site → 驗證：23 modules、site/print-all.html 與 data bundles 建置成功。
7. 執行獨立 QA → 驗證：13/13 findings、26 HTML、2,119 refs、163 standalone SVG 與版面檢查皆通過。

## 輸入、命令與輸出

### 輸入路徑

- /tmp/lab_tutorial_p0.GIEknc/repo
- InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/lab_tutorial_delta_reaudit.tsv
- InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json
- InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv
- InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md

### 核心執行命令

~~~bash
cd /tmp/lab_tutorial_p0.GIEknc/repo

git rev-parse HEAD
git status --short

node tools/build.mjs
node tools/verify.mjs
node tools/audit_order.mjs
python3 tools/check_svg_layout.py
git diff --check
~~~

另以 Python standard library 執行：

- 26 個 generated HTML 的 document skeleton、head title、duplicate ID、local file 與 fragment 檢查。
- 163 個 standalone SVG 的 XML 與 viewBox 檢查；src/svg/_defs.svg 是 fragment library，未誤當 standalone SVG。
- 13 個 finding 的 required-presence 與 forbidden-residual assertions。

### 實際輸出摘要

~~~text
HEAD 9eb1618d359e602d9c528675952b20d051fb2346

建置完成：23 個模組，387 ms
詞彙 97 條 · guardrail 21/21 就位
draft 1 篇：sr6.html
✓ 無錯誤

verify：66 個檔案（26 HTML · 20 JS · 7 CSS）· 6.38 MB
✓ file:// 安全檢查全數通過

SVG layout：✓ 無文字重疊、無超出畫布
HTML DOM/link：26 files、2,119 local refs，PASS
SVG XML：163 standalone files，PASS
finding assertions：13/13 PASS
git diff --check：exit 0
~~~

非阻塞 observation：

- verify.mjs 警告 glossary 的 epimutation 尚未被模組以 glossary token 引用。
- audit_order.mjs 指出 tumour purity、tumour DNA fraction、copy number 首次使用早於正式講解章節；命令仍 exit 0。
- 第一輪 SVG layout 曾找出 5 個文字溢出／重疊；逐一修正後重跑為 0。

## Finding-by-finding correction matrix

| Finding | Before（問題） | After（本機 source／generated site） | 推論與證據 | Disposition |
|---|---|---|---|---|
| OLD-P0-INDEX | 首頁把 read 還原成正常／癌細胞，並把 DNA fraction 當 cellular purity。 | 改成局部 haplotype 指派、somatic-supporting molecules 與 tumour DNA fraction；明示 read≠cell、DNA fraction≠cellular purity。 | Authority 只允許局部分子候選；bulk read 沒有來源細胞標籤。 | LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED；live pending |
| OLD-P0-PRINTALL | Aggregate 繼承來源頁 P0／regression。 | 先修 source，再由 build 重建；required／forbidden scan 全過。 | Build 成功不等同科學正確，因此另做 8 個 print-all assertions。 | LOCAL_GENERATED_SITE_CORRECTED；live pending |
| OLD-P0-M09 | methyl/read pattern 被升格為 cell group，且 M1 與 formal denominator 混用。 | 改為 fixed exact-PS×HP×genetic-pattern regional association；M1=102,842/469,849=21.8883%，formal=1,045 frozen → 811 evaluable → 3 robust=0.3699%；6/6 有 M1 screen、5/6 有 formal units、COLO829 formal=0。 | Denominator registry；formal association 不識別 clone、cell 或 lineage。 | LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED；live pending |
| OLD-P0-M12 | Somatic haplotag 被說成同一癌細胞群，benchmark scope 過廣。 | 改為同一局部 germline／somatic haplotype-carrier read pattern；HP label 不識別細胞。數字限於 bioRxiv v1、8 technical datasets、6 cell lines 與 BAM synthetic mixtures，非病人 cohort／臨床或外部驗證。 | Paper scope 與 authority cellular ceiling。 | LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED；live pending |
| OLD-P0-M13 | 真變異被說成必然形成一格路徑；GT3 當 cell groups；synthetic f=p。 | 改為指定 triplet heuristic 下的 model-compatible candidate；GT3 是 technical read-pattern branch；加入 f=pκ/[pκ+2(1-p)]，僅 κ=2 等條件成立時 f=p，並同步 quiz／SVG。 | Claim contract recurrence allowed；DNA fraction 與 cellular purity 的分母不同。 | LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED；live pending |
| OLD-P0-SR1 | 長期 clone-tree 願景被寫成 current output；methylation 被當 universal clock。 | 首屏分開 long-term vision 與 current exact-PS×HP read-compatible candidates；confirmed cellular clone／linear ancestry 皆為 0。甲基化只稱 bounded pattern-conditioned regional association，需 site／tissue／assay／model calibration。 | Authority manifest 與 claim contract；跨單位速率分母不可直接相除成通用時鐘。 | LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED；live pending |
| OLD-P0-SR2B | 舊 estimator 未 redirect；71,955／85,941 被當衝突；0.909 被當 truth。 | 加 SUPERSEDED banner 並導向 SR2C；同步 98,955 → 85,941 → 75,224 → 71,955 nested funnel、63,506/71,955=88.2579%；0.909 只稱 technical reproducibility。 | Denominator registry；重現性不等於 truth、accuracy 或 biological validity。 | LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED；live pending |
| DELTA-NEW-001 | HP families 被等同實體 chromosome copies／homologs。 | 統一為 phase-tag read groups；僅在 diploid、CN-neutral、正確 phasing、無 artifact 等條件下近似 homolog backgrounds；不識別實體 copy、親源或細胞。 | Exact-PS×HP 是技術容器；CN／LOH 未校正且 HP 方向任意。 | LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED；live pending |
| DELTA-NEW-002 | 斷言 family likelihood product 必錯，misassignment 會 duplicate molecule 並固定高估。 | 明示 read 只觀測一次；共享參數本身不使乘積失效。joint-window factor 是 proposal；exact 或 composite 取決於完整 observation model、條件獨立與 latent marginalization，偏誤方向須 simulation／calibration。 | 機率分解條件；misclassification 不等於 duplicate observation。 | LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED；live pending |
| DELTA-NEW-004 | 宣稱 ξ 高一數量級且主導 k≥2，ε 大致無害。 | 改列 allele error、family transition、mapping、chimera、partial coverage、phasing、correlated error 等 candidate channels；量級、region dependence、identifiability 與主導性均待 truth／simulation calibration。 | 無 ξ decomposition receipt；多個 error mechanisms 可產生同形狀假狀態。 | LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED；live pending |
| DELTA-NEW-005 | 把 δ·r background probability 直接當 detection limit。 | 只稱 modelled／expected background probability；detection limit 另依賴 n、α、target power、multiple testing、state competition、overdispersion、detector 與 calibration，須逐 unit 做 simulation／power analysis。 | 背景機率不是決策規則的 power boundary。 | LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED；live pending |
| DELTA-NEW-006 | 宣稱五參數 optimizer 幾毫秒收斂。 | 移除 milliseconds；明示低維／convex domain 不保證 objective concave，須報 optimizer、multi-start、硬體、unit size、wall time、iterations、convergence／failure rate 與 memory。 | 無 runtime benchmark receipt，且原頁自身要求 multi-start。 | LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED；live pending |
| DELTA-NEW-008 | 未驗證 LiLT 式子被稱為 production likelihood。 | 改為 proposed implementation likelihood；須先通過 simulation truth、calibration、independent truth set 與 versioned deployment gate。 | Authority 禁止 calibrated likelihood/posterior claim。 | LOCAL_SOURCE_AND_GENERATED_SITE_CORRECTED；live pending |

## 變更檔案與輸出位置

Clean clone 共 54 個 tracked files 有本機修改：

- src/modules/：6 個 authoritative module sources。
- src/data/：3 個 glossary／module metadata／quiz sources。
- src/svg/：28 個 scientific diagrams。
- tools/build.mjs：1 個首頁 source builder。
- site/：16 個 build-generated HTML／data bundles，含 site/print-all.html。

為避免 `/tmp` clone 消失後無法交付，另保存**authoritative source-only patch**；generated
`site/` 不納入 patch，maintainer 套用後應以 `node tools/build.mjs` 重生，避免雙 source-of-truth：

- `InterSubMod/research/20260813_public_docs_p0_correction/ccu_source_correction_9eb1618.patch`
- baseline：`9eb1618d359e602d9c528675952b20d051fb2346`
- patch：38 source files，799 insertions／1,117 deletions，260,072 bytes
- SHA-256：`aee944790ad8b60634818fca8888b0f980b3c1e5522ef64c00015a0a21d6c3e6`
- fresh detached baseline worktree：`git apply --check` exit 0

套用命令：

```bash
cd /path/to/clean/lab-tutorial
test "$(git rev-parse HEAD)" = 9eb1618d359e602d9c528675952b20d051fb2346
test -z "$(git status --short)"
git apply --check \
  /big7_disk/liaoyoyo2001/InterSubMod/research/20260813_public_docs_p0_correction/ccu_source_correction_9eb1618.patch
git apply \
  /big7_disk/liaoyoyo2001/InterSubMod/research/20260813_public_docs_p0_correction/ccu_source_correction_9eb1618.patch
node tools/build.mjs
node tools/verify.mjs
```

本收據：

- InterSubMod/research/20260813_public_docs_p0_correction/ccu_correction_receipt.md

## 未納入本子任務的剩餘項目

以下仍不可被本收據標成 resolved：

- OLD-P1：GLOSSARY、M02、M03、M05、M06、M07、M11、SR2、SR2C、SR4、SR6。
- DELTA-NEW-003、007、009、010、011。
- 其中 live source 仍可搜尋到待查數字／敘述，例如 SR2C 的 64.89%、n=60 → 0.07、
  H2009 唯一到 k=12／19× cost、LOH「安靜地錯」等；需依各 finding 另做 evidence receipt，
  不能因本輪 P0 build 通過而視為已驗證。
- SR5、DELTA-NUM-012、DELTA-NUM-013 在 delta re-audit 中原已標 RESOLVED，本輪未重開。

## Maintainer 後續動作

1. Review clean clone 的 54-file diff，特別檢查 source 與 generated site/print-all.html 是否同一批納入。
2. 透過正常 PR／review 流程 commit 與 push；本輪沒有代做。
3. 正式部署 GitHub Pages，記錄 deploy commit SHA、workflow run 與 live URL。
4. 對 live Pages 重跑 HTTP、byte/commit provenance、13 finding residual scan 與 browser smoke test。
5. 另開 P1／non-REGRESSED correction cycle；不得把本收據的 13/13 誤解成整站所有 finding 已清零。
