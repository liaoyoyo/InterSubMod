<!--
建立時間: 2026-08-13 17:49 +08:00
目標: 凍結本輪修改前 InterSubMod GitHub 與 Pages 的 live 狀態，分開 local correction 與 remote publication
處理範圍: GitHub repository metadata、main/develop/feature heads、17 個 GitHub Pages routes；全程只讀
關聯檔案:
  - InterSubMod/research/20260813_public_docs_p0_correction/remote_state_receipt.md
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/source_scope.tsv
-->

# Remote state receipt — 2026-08-13

## 結論

- GitHub About description 已更新為：`ONT read-level somatic-mutation integration for local mutation-state candidate analysis and methylation association; research software.`；C108 的 remote wording 本輪重查為 **RESOLVED_LIVE**。
- GitHub homepage 欄位仍為空字串；不是內容錯誤，但 Pages 導覽可見性仍可改善。
- Default `main` 仍為 `635437a65a33f8ba698acf85b22ebb069455c6cc`；live README 仍是舊版高估能力文字。
- `develop` 與 `feat/lineage-tag-methylation-axes` 均仍為 `ddd8909a838318d8a77969313e9561c8ff9d01c2`。
- Wiki `master` 仍為 `6cfc990f1fbe0b0aa76b8adbe56bbe645ef7ee7b`。
- 17/17 Pages routes 為 HTTP 200，但 17/17 SHA-256 均與 2026-08-12 frozen deployed bytes 相同；本輪 local correction 尚未發布。

## 輸入

- `https://api.github.com/repos/liaoyoyo/InterSubMod`
- `https://github.com/liaoyoyo/InterSubMod.git`
- `https://github.com/liaoyoyo/InterSubMod.wiki.git`
- `https://liaoyoyo.github.io/InterSubMod/explain/{index,01..16}.standalone.html`

## 執行命令

```bash
curl -fsSL https://api.github.com/repos/liaoyoyo/InterSubMod \
  | jq '{description,homepage,default_branch,pushed_at,updated_at}'
git ls-remote https://github.com/liaoyoyo/InterSubMod.git \
  refs/heads/main refs/heads/develop refs/heads/feat/lineage-tag-methylation-axes
git ls-remote https://github.com/liaoyoyo/InterSubMod.wiki.git refs/heads/master refs/heads/main
# 對 index 與 01..16 的固定 route 逐頁執行：
curl -L -sS -o /dev/null -w '%{http_code}' "$url"
curl -L -fsS "$url" | sha256sum
```

Exit code：`0`。

## 實際輸出片段

```text
description: ONT read-level somatic-mutation integration for local mutation-state candidate analysis and methylation association; research software.
homepage: ""
default_branch: main
635437a65a33f8ba698acf85b22ebb069455c6cc  refs/heads/main
ddd8909a838318d8a77969313e9561c8ff9d01c2  refs/heads/develop
ddd8909a838318d8a77969313e9561c8ff9d01c2  refs/heads/feat/lineage-tag-methylation-axes
6cfc990f1fbe0b0aa76b8adbe56bbe645ef7ee7b  refs/heads/master (Wiki)
index  200  deecd51ed7446185551de47c84903f8cfb8477de2a5317b36c1202375b9941b6
01_background-glossary  200  5e3206cbaa4d67a60c93904ec8c198006cea722df6b29198be401548d48af4a2
...
16_how-to-run  200  41d41079455f976d4b18cec97a51f0bf5a5de9b895349abe35ee944edf64bd0e
```

## 狀態語意

`LOCAL_SOURCE_CORRECTED` 不等於 `LIVE_RESOLVED`。只有 About 已由 live API 與 GitHub repository page 共同確認；default README、Wiki 與 Pages 必須在未來取得明確發布授權後 push／merge／publish，再重新抓取 live bytes 才能閉環。

## 封版前唯讀重查（2026-08-13 19:30:26 +08:00）

- GitHub API 描述、homepage 與 default branch 狀態未變。
- `main=635437a65a33f8ba698acf85b22ebb069455c6cc`；`develop` 與 feature 均為 `ddd8909a838318d8a77969313e9561c8ff9d01c2`；Wiki `master=6cfc990f1fbe0b0aa76b8adbe56bbe645ef7ee7b`。
- 17/17 Pages 仍為 HTTP 200，且每頁 SHA-256 均與 2026-08-12 frozen deployed bytes 相同。
- 結論不變：About 已上線；default `main`、Wiki 與 Pages 仍待發布。
