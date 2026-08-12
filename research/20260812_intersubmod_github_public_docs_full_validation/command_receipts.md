<!--
建立時間: 2026-08-12 Asia/Taipei
目標: 保存 InterSubMod 公開 README／Wiki／Explanation Center 的 fresh 機械驗證命令與輸出片段
處理範圍: Task Type B；current local、remote main/develop/feature、lab-tutorial deployed Pages、LongLineage remote feature
關聯檔案: InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/verification_results.tsv
-->

# 公開文件完整驗證 command receipts

這份紀錄服務 G2／G3／G4／G5。判定只綁定下列 commit 與輸入；未把 technical PASS 升格成 biological validation。所有 repo 內既有變更均保留，本輪只新增本檔與 `verification_results.tsv`；build／smoke outputs 放在 `/tmp`。

## 0. Authority／版本矩陣

| 對象 | 本輪 authority | commit／狀態 |
|---|---|---|
| 使用者指定的 lab-tutorial URL | GitHub Pages Actions artifact；以 HTTP blob hash 對 remote main `site/` | `9eb1618d359e602d9c528675952b20d051fb2346` |
| lab-tutorial 舊 `gh-pages` ref | 非本輪 deployed authority | `a4b45a8e5238df5013ffdc14a13da55da7752aa8` |
| InterSubMod current local | 實際 build/run source | `73afaeac8e61c767241fa59c1ca6043a1c95290c-dirty` |
| InterSubMod remote feature／develop | 公開 README/Wiki 文件軌 | `ddd8909a838318d8a77969313e9561c8ff9d01c2` |
| InterSubMod remote main | GitHub default main | `635437a65a33f8ba698acf85b22ebb069455c6cc` |
| LongLineage clean-tested remote feature | clean clone build/run | `b9aaa12a11fa00606bd174dabd0f172a5d112359` |
| LongLineage remote develop／main | 未在本輪 build | `cec375da5cdff87c2b2d66563cf55a6f0eb957bb`／`5daf50f04cbe233abfade816ce9e0903f6b38954` |

### R00 — remote ref pin

輸入：三個 GitHub repo 的 heads。

完整命令：

```bash
git ls-remote --heads https://github.com/liaoyoyo/InterSubMod.git main develop feat/lineage-tag-methylation-axes
git ls-remote --heads https://github.com/liaoyoyo/LongLineage.git main develop feat/lineage-tagging-toolchain
git ls-remote --symref https://github.com/CCU-Bioinformatics-Lab/lab-tutorial.git HEAD
git ls-remote --heads https://github.com/CCU-Bioinformatics-Lab/lab-tutorial.git
```

輸出路徑：無（read-only remote query）。退出碼：`0`。

實際輸出片段：

```text
ddd8909a838318d8a77969313e9561c8ff9d01c2 refs/heads/develop
ddd8909a838318d8a77969313e9561c8ff9d01c2 refs/heads/feat/lineage-tag-methylation-axes
635437a65a33f8ba698acf85b22ebb069455c6cc refs/heads/main
cec375da5cdff87c2b2d66563cf55a6f0eb957bb refs/heads/develop
b9aaa12a11fa00606bd174dabd0f172a5d112359 refs/heads/feat/lineage-tagging-toolchain
5daf50f04cbe233abfade816ce9e0903f6b38954 refs/heads/main
ref: refs/heads/main HEAD
9eb1618d359e602d9c528675952b20d051fb2346 refs/heads/main
a4b45a8e5238df5013ffdc14a13da55da7752aa8 refs/heads/gh-pages
```

## 1. lab-tutorial deployed scope

### R01 — deployed blob 是否等於 remote main

輸入：`https://ccu-bioinformatics-lab.github.io/lab-tutorial/index.html`、temp clone `/tmp/intersubmod-tutorial-audit.bZP2Ov/lab-tutorial`。

完整命令：

```bash
URL=https://ccu-bioinformatics-lab.github.io/lab-tutorial/index.html
REPO=/tmp/intersubmod-tutorial-audit.bZP2Ov/lab-tutorial
curl -L -sS -o /dev/null -w '%{http_code}\n' "$URL"
curl -L -sS "$URL" | sha256sum
git -C "$REPO" show origin/main:site/index.html | sha256sum
git -C "$REPO" show origin/gh-pages:index.html | sha256sum
```

輸出路徑：無。退出碼：`0`。

實際輸出：

```text
HTTP_STATUS=200
DEPLOYED_SHA256=bc0821bcd5ea6df2326d1b7bb6560d518f8e31599b4ddf34e0d1f529a4a0466f
REMOTE_MAIN_SITE_SHA256=bc0821bcd5ea6df2326d1b7bb6560d518f8e31599b4ddf34e0d1f529a4a0466f
GH_PAGES_SHA256=b47c6abbc13f60df4708ba74358dfb26987927ba042ccc094d6e0cb8802a69a9
```

判定：目前 URL 是 remote main 的 Actions deployment，不是舊 `gh-pages` branch blob。

### R02 — 公開頁面數與 draft

輸入：remote main 的 `site/*.html` 共 26 個。

完整命令：

```bash
REPO=/tmp/intersubmod-tutorial-audit.bZP2Ov/lab-tutorial
base=https://ccu-bioinformatics-lab.github.io/lab-tutorial
for p in $(git -C "$REPO" ls-tree -r --name-only origin/main site | awk -F/ '/\.html$/{print $2}'); do
  http=$(curl -L -sS -o "/tmp/intersubmod-public-audit.XNde4y/published-$p" -w '%{http_code}' "$base/$p")
  # 對 HTTP 200 頁逐一比較 deployed 與 origin/main:site blob SHA-256
done
```

輸出路徑：`/tmp/intersubmod-public-audit.XNde4y/published-*.html`。退出碼：`0`。

實際輸出摘要：

```text
index.html HTTP=200
m00.html ... m13.html HTTP=200
sr1.html ... sr5.html HTTP=200
sr6.html HTTP=404
HTTP_200=25 NON_200=1 HASH_MISMATCH_AMONG_200=0
```

`src/data/modules.json` 將唯一一頁標成 `"draft": true`；deploy workflow 的 `strip_drafts.mjs` 會在 Pages artifact 移除它。因此公開正式頁為 25，不是部署遺漏。

## 2. Canonical authority 與分母

### R03 — 13 個 authority artefact hash

輸入：`InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json`。

完整命令：

```bash
AUTH=docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/authority_manifest.json
jq -r '.artifacts[] | [.artifact_id,.path,.sha256] | @tsv' "$AUTH" |
while IFS=$'\t' read -r id path expected; do
  actual=$(sha256sum "$path" | cut -d' ' -f1)
  [ "$actual" = "$expected" ] && status=MATCH || status=MISMATCH
  printf '%s\t%s\t%s\n' "$id" "$status" "$actual"
done
```

輸出路徑：stdout。退出碼：`0`。

實際輸出片段：

```text
strict_linkage_ready MATCH 1c95ed3c6509e9477b2ab38dc380d16e8072e0ba0ce8e17f82842dee7a0dd2ee
all7_summary MATCH d2edca043564f9da4a5372dc70542b89024040a9b41a6290ccf938d262b46bcd
topology_summary MATCH 3f77559d57b34eebda05d482176738a31ad983e47b7ee26499d69aa2d2984194
methyl_report_data MATCH 98a980c9fa28ed9859fe9656f3ab089ccfddd7eb907bd3580947ae375d56e3fc
methyl_reader_report MATCH b78ad76a76799b0d116fa2d4f36a3fda9a0602e715e6ccda56bc7e2afb6cc8aa
```

總計：`13/13 MATCH`，0 missing，0 mismatch。

### R04 — 7 datasets／6 biological IDs／469,849

輸入：`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/latest_tree_input_contract_audit.json`（SHA-256 `526f89fc1d5b87c052c6a48fb70820aa87349ee8508ba4f938bd864328678d42`）。

完整命令：

```bash
jq '{dataset_count:(.samples|length), biological_ids:([.samples[].biological_id]|unique), biological_count:([.samples[].biological_id]|unique|length), sum_all_ssnv:([.samples[].all_ssnv_count]|add), all_pass:([.samples[].pass]|all)}' \
  research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/latest_tree_input_contract_audit.json
```

輸出路徑：stdout。退出碼：`0`。

實際輸出：

```json
{"dataset_count":7,"biological_ids":["COLO829","H1437","H2009","HCC1395","HCC1937","HCC1954"],"biological_count":6,"sum_all_ssnv":469849,"all_pass":true}
```

### R05 — funnel 分母重算

輸入：`InterSubMod/docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv`（SHA-256 `a41f726cbf66f22b7e95ddb44216a08bf480d12961758ebd5b1ab6e3e61db6a4`）。

完整命令：

```bash
awk -F '\t' 'NR>1 && ($1=="k1_linkage_abstain" || $1=="one_rooted_unlabeled_topology" || $1=="methyl_evaluable_units" || $1=="methyl_robust_associations") {printf "%s\t%s/%s=%.10f%%\trecorded=%s%%\n",$1,$2,$3,100*$2/$3,$4}' \
  docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv
awk 'BEGIN{printf "50432/469849=%.10f%%\n",100*50432/469849; printf "170131/469849=%.10f%%\n",100*170131/469849}'
```

輸出路徑：stdout。退出碼：`0`。

實際輸出：

```text
k1_linkage_abstain 170131/255752=66.5218649317% recorded=66.5219%
one_rooted_unlabeled_topology 63506/71955=88.2579389896% recorded=88.2579%
methyl_evaluable_units 811/1045=77.6076555024% recorded=77.6077%
methyl_robust_associations 3/811=0.3699136868% recorded=0.3699%
50432/469849=10.7336612401%
170131/469849=36.2097184415%
```

關鍵更正：README 把 `170,131/255,752` 的 component-level 66.52% 寫成「全部 mutations 約三分之二 isolated」，分母與 grain 都錯。Current positional-singleton authority 是 `50,432/469,849=10.7337%`。

## 3. InterSubMod fresh build／tests／smoke

### R06 — source 與 remote feature build-input 邊界

輸入：current local worktree 與 `origin/feat/lineage-tag-methylation-axes=ddd8909a`。

完整命令：

```bash
git diff --quiet origin/feat/lineage-tag-methylation-axes -- CMakeLists.txt cmake src include tests external
printf 'TRACKED_BUILD_INPUT_DIFF_EXIT=%s\n' "$?"
git status --porcelain --untracked-files=all -- CMakeLists.txt cmake src include tests external
git diff --name-status origin/feat/lineage-tag-methylation-axes..HEAD
```

輸出路徑：stdout。退出碼：`0`。

實際輸出：

```text
TRACKED_BUILD_INPUT_DIFF_EXIT=0
UNTRACKED_BUILD_INPUT_ROWS=0
A scripts/build_cn_hp_observation_pilot.py
```

這只證明 C++／CMake build inputs 對 remote feature byte-equivalent；下一節仍明標實際 build identity 是 `73afaeac-dirty`，不冒充 clean remote checkout build。

### R07 — configure／build

輸入路徑：`/big7_disk/liaoyoyo2001/InterSubMod`。

完整命令：

```bash
cmake -S /big7_disk/liaoyoyo2001/InterSubMod \
  -B /tmp/intersubmod-public-audit.XNde4y \
  -DCMAKE_BUILD_TYPE=Release
/usr/bin/time -f 'AUDIT_WALL_SECONDS=%e\nAUDIT_EXIT=%x' \
  cmake --build /tmp/intersubmod-public-audit.XNde4y -j"$(nproc)"
```

輸出路徑：`/tmp/intersubmod-public-audit.XNde4y/bin/`。Configure exit：`0`；build exit：`0`。

實際輸出片段：

```text
GNU 11.4.0
htslib 1.18
InterSubMod build identity: 73afaeac-dirty
[100%] Built target run_tests
AUDIT_WALL_SECONDS=242.81
AUDIT_EXIT=0
```

實際 warnings：`LabelTest.cpp:899 unused n_total_family`、`DistanceMatrix.cpp:562 unused parameter reads`、`DistanceMatrix.cpp:354 unused function calculate_distance`、`test_phase1_2.cpp` 的 `argc/argv` unused。Fresh first build也下載／編譯 jemalloc、Eigen、GoogleTest，故約 4 分鐘符合文件的 2–5 分鐘量級。

### R08 — direct GoogleTest 與 CTest

輸入：`/tmp/intersubmod-public-audit.XNde4y/bin/run_tests`。

完整命令：

```bash
/usr/bin/time -f 'AUDIT_WALL_SECONDS=%e\nAUDIT_EXIT=%x' \
  /tmp/intersubmod-public-audit.XNde4y/bin/run_tests
ctest --test-dir /tmp/intersubmod-public-audit.XNde4y -N
/usr/bin/time -f 'AUDIT_WALL_SECONDS=%e\nAUDIT_EXIT=%x' \
  ctest --test-dir /tmp/intersubmod-public-audit.XNde4y --output-on-failure
```

輸出路徑：stdout。兩命令 exit：`0`。

實際輸出：

```text
[==========] Running 270 tests from 39 test suites.
[==========] 270 tests from 39 test suites ran. (1972 ms total)
[  PASSED  ] 270 tests.
Total Tests: 270
100% tests passed, 0 tests failed out of 270
Total Test time (real) = 6.20 sec
```

判定：公開 `265 tests / 38 suites` 已少報 `5 tests / 1 suite`；「全過」仍成立。

### R09 — CLI defaults

輸入：`/tmp/intersubmod-public-audit.XNde4y/bin/inter_sub_mod`。

完整命令：

```bash
/tmp/intersubmod-public-audit.XNde4y/bin/inter_sub_mod --help
printf 'HELP_EXIT=%s\n' "$?"
/tmp/intersubmod-public-audit.XNde4y/bin/inter_sub_mod
printf 'NO_ARGS_EXIT=%s\n' "$?"
```

輸出路徑：stdout/stderr。實際 exit：help `1`；no args `1`。

實際輸出片段：

```text
--threads ... Default: 1
--distance-metric ... Default: NHD
HELP_EXIT=1
Error: --tumor-bam is required
NO_ARGS_EXIT=1
```

Source contract：`include/core/Config.hpp` 初始 `distance_metrics={BERNOULLI}`、`threads=16`；`ArgParser.hpp` local default `NHD` 並清除／重設 metrics。故文件的 thread mismatch 正確；若寫成「help 顯示 BERNOULLI」則錯，現行 help 顯示 NHD。

### R10 — single-locus smoke

輸入路徑：

- tumor BAM：`/big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam`
- reference：`/big7_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa`
- source VCF：`/big7_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz`
- temp VCF：`/tmp/intersubmod-public-audit.XNde4y/smoke/one_snv.vcf`

完整命令：

```bash
bcftools view -r chr19:29283968-29283968 -Ov \
  -o /tmp/intersubmod-public-audit.XNde4y/smoke/one_snv.vcf \
  /big7_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz
/usr/bin/time -f 'AUDIT_WALL_SECONDS=%e\nAUDIT_EXIT=%x' \
  /tmp/intersubmod-public-audit.XNde4y/bin/inter_sub_mod \
  --tumor-bam /big7_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/tumor.bam \
  --reference /big7_disk/liaoyoyo2001/InterSubMod/data/ref/hg38.fa \
  --vcf /tmp/intersubmod-public-audit.XNde4y/smoke/one_snv.vcf \
  --output-dir /tmp/intersubmod-public-audit.XNde4y/smoke/out_min
```

輸出路徑：`/tmp/intersubmod-public-audit.XNde4y/smoke/out_min/`。兩命令 exit：`0`。

實際輸出：

```text
Threads: 16
Distance Metrics: NHD
Total regions: 1 / Successful: 1 / Failed: 0
Total reads processed: 85
Forward strand (+): 40 / Reverse strand (-): 45
Total CpG sites found: 11
Metric: NHD / Total valid read pairs: 3443 / Total invalid pairs: 127
Warning: Region has fewer than 100 reads; KDE correction unavailable; using 75x fallback
AUDIT_WALL_SECONDS=1.90
AUDIT_EXIT=0
```

### R11 — fresh remote-develop clone inputs

輸入：`/tmp/intersubmod-github-audit.nvGQP6/develop-tree`，detached `ddd8909a838318d8a77969313e9561c8ff9d01c2`。

完整命令：

```bash
for p in data/bam/HCC1395/tumor.bam data/ref/hg38.fa \
  data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz one_snv.vcf \
  build/bin/inter_sub_mod scripts/run_vcf_all_snv.sh requirements.txt; do
  test -e "/tmp/intersubmod-github-audit.nvGQP6/develop-tree/$p" && echo "$p PRESENT" || echo "$p MISSING"
done
```

輸出路徑：stdout。退出碼：`0`。

實際輸出：

```text
data/bam/HCC1395/tumor.bam MISSING
data/ref/hg38.fa MISSING
data/vcf/HCC1395/pileup/filtered_snv_tp.vcf.gz MISSING
one_snv.vcf MISSING
build/bin/inter_sub_mod MISSING
scripts/run_vcf_all_snv.sh PRESENT
requirements.txt PRESENT
```

判定：README 的 build 可做，但範例分析不是 fresh-clone self-contained smoke test。

## 4. 文件數量與 evidence-chain

### R12 — standalone pages／inline SVG

輸入：InterSubMod current worktree、HEAD、remote refs。

完整命令：

```bash
for ref in HEAD origin/feat/lineage-tag-methylation-axes origin/main origin/develop; do
  pages=$(git ls-tree -r --name-only "$ref" docs/explain | rg '\.standalone\.html$' | wc -l)
  svg=$(git grep -h -o -E '<svg([[:space:]]|>)' "$ref" -- 'docs/explain/*.standalone.html' | wc -l)
  printf '%s\t%s\tpages=%s\tsvg=%s\n' "$ref" "$(git rev-parse "$ref")" "$pages" "$svg"
done
```

輸出路徑：stdout。退出碼：`0`。

實際輸出：

```text
HEAD 73afaeac... pages=17 svg=37
origin/feat/lineage-tag-methylation-axes ddd8909a... pages=17 svg=37
origin/develop ddd8909a... pages=17 svg=37
origin/main 635437a6... pages=11 svg=29
WORKTREE 73afaeac-dirty pages=17 svg=37
```

判定：17 pages 對 feature/develop 成立；README 的 29 SVG 對目前文件軌已 stale，37 是可重現 raw `<svg>` count。

### R13 — 334 source refs +111 path claims receipt 缺失

輸入：HEAD、remote develop/main、worktree。

完整命令：

```bash
for ref in HEAD origin/develop origin/main; do
  for p in scratchpad/collect_references.json scratchpad/collect_paths.json; do
    git cat-file -e "$ref:$p" 2>/dev/null
    printf '%s\t%s\texists=%s\n' "$ref" "$p" "$?"
  done
done
for p in scratchpad/collect_references.json scratchpad/collect_paths.json; do
  test -e "$p"
  printf 'WORKTREE\t%s\texists=%s\n' "$p" "$?"
done
```

輸出路徑：stdout。命令最終 exit：`0`（個別 `git cat-file` 為 `128`，`test` 為 `1`）。

實際輸出：

```text
HEAD scratchpad/collect_references.json exists=128
HEAD scratchpad/collect_paths.json exists=128
origin/develop scratchpad/collect_references.json exists=128
origin/develop scratchpad/collect_paths.json exists=128
origin/main scratchpad/collect_references.json exists=128
origin/main scratchpad/collect_paths.json exists=128
WORKTREE scratchpad/collect_references.json exists=1
WORKTREE scratchpad/collect_paths.json exists=1
```

判定：目前只能找到 commit message 的總結，不能重建逐項 population 或驗證 0 fabricated／0 out-of-range。

## 5. Storage claims

### R14 — 7 current raw BAM／sidecar byte inventory

輸入：`latest_tree_input_contract_audit.json` 的 7 組 `raw_bam` 與 `latest_hp_ps_sidecar`。

完整命令：

```bash
MAN=research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/latest_tree_input_contract_audit.json
jq -r '.samples[] | [.sample,.raw_bam,.latest_hp_ps_sidecar] | @tsv' "$MAN" |
while IFS=$'\t' read -r sample bam sidecar; do
  printf '%s\t%s\t%s\n' "$sample" "$(stat -Lc '%s' "$bam")" "$(stat -Lc '%s' "$sidecar")"
done
# 對 7+7 筆 bytes 加總，再換算 TiB/GiB 與 ratio
```

輸出路徑：stdout；中間 inventory `/tmp/intersubmod-public-audit.XNde4y/size_inventory.tsv`。退出碼：`0`。

實際輸出：

```text
HCC1395        292055926761  1433423193
HCC1395_DORADO 250150882482  1646911913
COLO829        100423632767   339540321
H1437          243719059370   585465558
H2009          327679299934   886667909
HCC1937        472119759082   728870441
HCC1954        253211723892   635288829
raw_bytes=1939360284288 raw_TiB=1.763837903389
side_bytes=6256168164 side_GiB=5.826510641724
ratio=309.991712730438
```

判定：`5.83 GiB` fresh verified。`1.67 TiB` 宣稱的是沒有落地的 hypothetical tagged BAM，repo 未提供逐檔 registry；current raw total 是 1.764 TiB，不能冒充 tagged total。僅按展示數字計算，`1.67 TiB / 5.83 GiB = 293.34x`，不是 287x。`>99%` SEQ+QUAL 也沒有 component byte census。

## 6. LongLineage clean feature build／real-data frozen report

### R15 — clean source pin、configure、build

輸入：`/tmp/longlineage-public-audit.NhzTdr/src`；HEAD `b9aaa12a11fa00606bd174dabd0f172a5d112359`；tracked status rows `0`。

完整命令：

```bash
/usr/bin/cmake -S /tmp/longlineage-public-audit.NhzTdr/src \
  -B /tmp/longlineage-public-audit.NhzTdr/build \
  -DCMAKE_BUILD_TYPE=Debug -DLONGLINEAGE_BUILD_TESTS=ON
/usr/bin/time -f 'AUDIT_WALL_SECONDS=%e\nAUDIT_EXIT=%x' \
  /usr/bin/cmake --build /tmp/longlineage-public-audit.NhzTdr/build -j4
```

輸出路徑：`/tmp/longlineage-public-audit.NhzTdr/build/bin/`。兩命令 exit：`0`。

實際輸出片段：

```text
GNU 11.4.0; Boost 1.74; HTSlib 1.18; Jansson 2.13.1; OpenSSL 3.0.2; zlib 1.2.11
[100%] Built target test_validator_fault
AUDIT_WALL_SECONDS=85.10
AUDIT_EXIT=0
```

### R16 — LongLineage tests／check_all

輸入：同一 clean build。

完整命令：

```bash
ctest --test-dir /tmp/longlineage-public-audit.NhzTdr/build -N
/usr/bin/time -f 'AUDIT_WALL_SECONDS=%e\nAUDIT_EXIT=%x' \
  ctest --test-dir /tmp/longlineage-public-audit.NhzTdr/build --output-on-failure
/usr/bin/time -f 'AUDIT_WALL_SECONDS=%e\nAUDIT_EXIT=%x' \
  bash scripts/ci/check_all.sh /tmp/longlineage-public-audit.NhzTdr/build
```

輸出路徑：stdout。兩個執行 exit：`0`。

實際輸出：

```text
Total Tests: 47
100% tests passed, 0 tests failed out of 47
Total Test time (real) = 148.87 sec
AUDIT_EXIT=0
NO-NETWORK RESULT: PASS authority=linux-network-namespace
CHECK ALL RESULT: FOUNDATION_PASS build=/tmp/longlineage-public-audit.NhzTdr/build (release-phase blockers remain governed)
AUDIT_WALL_SECONDS=195.77
AUDIT_EXIT=0
```

### R17 — LongLineage CLI boundary

完整命令：

```bash
/tmp/longlineage-public-audit.NhzTdr/build/bin/longlineage --version
/tmp/longlineage-public-audit.NhzTdr/build/bin/longlineage --help
/tmp/longlineage-public-audit.NhzTdr/build/bin/longlineage
```

輸出路徑：stdout/stderr。實際 exit：`0`、`0`、`2`。

實際輸出片段：

```text
longlineage 0.1.0
Usage:
  longlineage preflight --manifest FILE [--repo DIR]
  longlineage run --manifest FILE [--repo DIR]
  longlineage dataset-gate --manifest FILE [--repo DIR]
  longlineage probe --manifest FILE --partial-output DIR [--repo DIR]
Production run is fail-closed until P3, P4 and P5 are VERIFIED.
HCC1395 dataset-gate writes READY_FOR_VALIDATION staging only;
longlineage-validate must independently validate and freeze it.
```

### R18 — HCC1395 frozen real-data report

輸入：clean clone 的 `docs/reports/20260720_HCC1395完整科學運算與parity報告_01.json`；SHA-256 `7e0b650b14c3a4fe6ed4dace4c9ef2cb5e56b8257a559cbf6d7d8ecc5fa562eb`。

完整命令：

```bash
jq '{scope,production_claim_allowed,
  counts:{site_keys:.counts.site_keys,m1_stable_assignments:.counts.m1_stable_assignments,m2_eligible:.counts.m2_eligible,topology_primary_hp_units:.counts.topology_primary_hp_units},
  artifact_rows:[.artifacts[]|{role:(.role // .artifact_id // .path),logical_rows}]}' \
  docs/reports/20260720_HCC1395完整科學運算與parity報告_01.json
awk 'BEGIN{printf "NOT_STABLE=%d\n",79687-12851; printf "NOT_STABLE_RATE=%.10f%%\n",100*(79687-12851)/79687}'
```

輸出路徑：stdout。退出碼：`0`。

實際輸出片段：

```text
scope.dataset_count=1
scope.dataset_ids=["HCC1395"]
scope.completeness="FULL"
production_claim_allowed=false
site_keys=79687
m1_stable_assignments=12851
m2_eligible=5
topology_primary_hp_units=0
cooccurrence_pairs logical_rows=134278
cooccurrence_sites logical_rows=79687
topology_units logical_rows=0
NOT_STABLE=66836
NOT_STABLE_RATE=83.8731537139%
```

判定：0 topology 對這個 frozen one-dataset gate 成立；但 README 的「66,836 sites never produce a co-occurrence row at all」被 artifact census 直接反駁。正確是 66,836 未進 M1 stable assignments；co-occurrence site/pair artefacts仍有 79,687／134,278 rows。

## 7. 最終機械結論

1. Canonical funnel 的大多數整數與正確分母都可重算；最嚴重的數據敘述錯誤是把 `170,131/255,752` component statistic 寫成「全部 469,849 mutations 的三分之二 isolated」。
2. `88.26%` 本身正確，但只能讀成 `63,506/71,955` ranked mathematical shape statistic。
3. Fresh C++ tests 全綠，但公開 `265/38` 已是舊快照；content-equivalent build 現為 `270/39`。
4. Fresh-clone README 缺 BAM／FASTA／VCF／one-site fixture，範例不是可攜 smoke test。
5. `17 pages` 對 feature/develop 正確；`29 SVG` 對該文件軌已增至 `37`。
6. `334+111, 0 errors` 沒有 committed item registry，不能 fresh 重驗。
7. `5.83 GiB` 可由 current 7 sidecars 驗證；`1.67 TiB`、`287x`、`SEQ+QUAL >99%` 的證據鏈不完整，其中 `287x` 與顯示單位的算術直接不相容。
8. LongLineage remote feature clean foundation build/tests pass；同時明示 release blockers。HCC1395 formal report為 one-dataset、non-production、0 topology，不可泛化。

## 8. CCU current delta、LongPhase contract 與 HTML QA 補充收據

### R19 — CCU lab-tutorial current-delta 完整 finding census

輸入路徑：

- baseline website source：`/tmp/lab_tutorial_delta.XI3q4E/repo` commit `46b6f5b3016c187ad742fecbfa813f835b09e605`
- current website source：同 repo commit `9eb1618d359e602d9c528675952b20d051fb2346`
- 舊稽核：`InterSubMod/research/20260811_ccu_lab_tutorial_claim_audit/`

完整命令：

```bash
git -C /tmp/lab_tutorial_delta.XI3q4E/repo diff --stat \
  46b6f5b3016c187ad742fecbfa813f835b09e605..9eb1618d359e602d9c528675952b20d051fb2346
python3 - <<'PY'
import csv, collections
p='research/20260812_intersubmod_github_public_docs_full_validation/lab_tutorial_delta_reaudit.tsv'
with open(p, encoding='utf-8') as f:
    rows=list(csv.DictReader(f, delimiter='\t'))
print(len(rows), collections.Counter(r['status'] for r in rows))
PY
sha256sum research/20260812_intersubmod_github_public_docs_full_validation/lab_tutorial_delta_reaudit.tsv
```

輸出路徑：`InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/lab_tutorial_delta_reaudit.tsv`。

實際輸出片段：

```text
57 files changed, 8105 insertions(+), 2774 deletions(-)
32 Counter({'OPEN': 17, 'REGRESSED': 6, 'PARTIAL': 6, 'RESOLVED': 3})
04f88a4ceb46644ee287940699575f001c0f89e387616ea6f3f1b5c6e259df06
```

判定：29/32 problem-focused findings 尚未完全解決；此分母不是網站全部句子。

### R20 — LongPhase-S/TO KB → pinned source → runtime 三層核對

輸入路徑：

- `/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-s.md`
- `/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-to.md`
- `/big8_disk/liaoyoyo2001/longphase-s/src/haplotag/HaplotagType.h`，local commit `420e0ff73188126da8bb2122725fbffb2d076ccc`、dirty
- `/big8_disk/liaoyoyo2001/longphase-to/HaplotagProcess.cpp`，local commit `9d74f8b4176ac1d64160f3ea5346a04a800815da`、dirty
- binaries SHA-256：LongPhase-S `754f996a...`；LongPhase-TO `5186d3d3...`

完整命令：

```bash
rg -n 'HP:|--pb|--ont|PacBio|ONT|haplotag' \
  /big8_disk/liaoyoyo2001/knowledge/05_tools/longphase-{s,to}.md
/big8_disk/liaoyoyo2001/longphase-to/longphase-to phase --help
sed -n '97,108p;327,340p' \
  /big8_disk/liaoyoyo2001/longphase-s/src/haplotag/HaplotagType.h
sed -n '430,450p' \
  /big8_disk/liaoyoyo2001/longphase-to/HaplotagProcess.cpp
sha256sum /big8_disk/liaoyoyo2001/longphase-to/longphase-to \
  /big8_disk/liaoyoyo2001/longphase-s/longphase-s
```

輸出路徑：stdout；本輪未改外部工具。退出碼：`0`。

實際輸出片段：

```text
--ont, --pb    ont: Oxford Nanopore genomic reads.
               pb: PacBio HiFi/CCS genomic reads.
LongPhase-S complete ReadHP strings: ./1/2/3/4/1-1/2-1/1-2/2-2
bam_aux_append(aln, "HP", 'i', sizeof(haplotype), ...)
```

判定：M6 的 `HP1-1-1` 與 TO-only-ONT 均被直接反駁；M7/glossary 的「整個 LongPhase-TO 為 `HP:Z`」也被 source contract 反駁。因 local external repos dirty，結論只綁上述 commit/worktree/binary hash；公開教材仍應釘正式 release。

### R21 — Standalone HTML desktop/mobile/print browser QA

輸入路徑：`InterSubMod/docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.standalone.html`。

完整命令：

```bash
python3 research/20260812_intersubmod_github_public_docs_full_validation/qa_html_report.py \
  docs/reports/validated/2026/08/20260812_InterSubMod_GitHub公開說明與教學完整驗證_01.standalone.html \
  --out-dir research/20260812_intersubmod_github_public_docs_full_validation/html_qa
```

輸出路徑：

- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/html_qa/qa_receipt.json`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/html_qa/desktop.png`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/html_qa/mobile.png`
- `InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/html_qa/print.pdf`

實際輸出片段（最終重播應一致）：

```text
desktop viewport=1440 documentScroll=1440
mobile viewport=390 documentScroll=390
console_errors=[]; page_errors=[]; broken_internal_links=[]
SVG count=4; all viewBox/role/title/desc=true
prohibited external assets/gradient/text-shadow/backdrop-filter=false
pass=true; exit=0
```
