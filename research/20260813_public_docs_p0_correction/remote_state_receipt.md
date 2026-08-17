<!--
建立時間: 2026-08-13 10:25 +08:00
目標: 將本地 P0 修正與 GitHub／CCU live publication state 分開驗證
處理範圍: InterSubMod main/develop/Wiki/Pages、GitHub About、CCU lab-tutorial main/live index
關聯檔案:
  - InterSubMod/research/20260812_intersubmod_github_public_docs_full_validation/source_scope.tsv
  - InterSubMod/research/20260813_public_docs_p0_correction/00_INDEX.md
-->

# Remote publication state receipt

## 結論

2026-08-13 只讀重查顯示，公開版本平面仍與 2026-08-12 locked audit 相同；本輪 working-tree
修正尚未發布。GitHub About 仍含 `subclone resolution`，default main 仍是舊 README，
`README.zh-TW.md` 仍 HTTP 404，Wiki／Pages live bytes 仍等於舊 hash，CCU remote main 仍是
`9eb1618`。因此本輪任何 local source／patch 成功都不能標成 `LIVE_RESOLVED`。

## 輸入與命令

輸入是下列公開端點；命令只讀，沒有 push／deploy：

```bash
git ls-remote https://github.com/liaoyoyo/InterSubMod.git \
  refs/heads/main refs/heads/develop
git ls-remote https://github.com/liaoyoyo/InterSubMod.wiki.git refs/heads/master
git ls-remote https://github.com/ccu-bioinformatics-lab/lab-tutorial.git \
  refs/heads/main refs/heads/master

curl -L -sS -o /dev/null -w '%{http_code}' <URL>
curl -L -fsS <URL> | sha256sum
```

## 實際輸出

### Remote refs

```text
635437a65a33f8ba698acf85b22ebb069455c6cc  InterSubMod main
ddd8909a838318d8a77969313e9561c8ff9d01c2  InterSubMod develop
6cfc990f1fbe0b0aa76b8adbe56bbe645ef7ee7b  InterSubMod Wiki master
9eb1618d359e602d9c528675952b20d051fb2346  CCU lab-tutorial main
```

### HTTP／hash

| Endpoint | HTTP | SHA-256 | 判定 |
|---|---:|---|---|
| `raw.githubusercontent.com/liaoyoyo/InterSubMod/main/README.md` | 200 | `8187a83223d538d0ef25b12722535397f67c0dade394298785ee626e52e54efa` | 舊 main，未發布本輪修正 |
| `raw.githubusercontent.com/liaoyoyo/InterSubMod/main/README.zh-TW.md` | 404 | `NA` | C141 仍為 live gap |
| `raw.githubusercontent.com/liaoyoyo/InterSubMod/develop/README.md` | 200 | `24929e94c663ec28a90cf7c337f11702c42af47eb68068c68d4372e96b3990eb` | 舊 feature README |
| `raw.githubusercontent.com/wiki/liaoyoyo/InterSubMod/Home.md` | 200 | `a17d9ee25119062c05ae0e101aaf7d5b201ce1436b5b4f9378ece1134a64175b` | 舊 Wiki bytes |
| `liaoyoyo.github.io/InterSubMod/explain/index.standalone.html` | 200 | `deecd51ed7446185551de47c84903f8cfb8477de2a5317b36c1202375b9941b6` | 舊 Pages bytes |
| `.../01_background-glossary.standalone.html` | 200 | `5e3206cbaa4d67a60c93904ec8c198006cea722df6b29198be401548d48af4a2` | 舊 Pages bytes |
| `.../08_subclone-logic-chain-chr2-18M.standalone.html` | 200 | `3f18b87edda37184e78ddd7386b4a25587871ea44ad38ea3d626413463030c6b` | 舊 Pages bytes |
| `ccu-bioinformatics-lab.github.io/lab-tutorial/index.html` | 200 | `bc0821bcd5ea6df2326d1b7bb6560d518f8e31599b4ddf34e0d1f529a4a0466f` | 仍對應 CCU `9eb1618` source |

GitHub repository 頁的 live About 文字仍為：

> InterSubMod - performs read-level integration of methylation profiles with haplotypes,
> somatic alleles, and tumor/normal labels, enabling somatic variant validation, subclone
> resolution, and single-molecule epigenomic clustering.

## 輸出與限制

- 本收據輸出：`InterSubMod/research/20260813_public_docs_p0_correction/remote_state_receipt.md`。
- 沒有下載或保存 remote body；hash 由 stdout pipeline 計算。
- 這是 2026-08-13 的 point-in-time check；發布後必須重跑，不可永久沿用。

