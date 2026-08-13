# `docs/wiki/` — GitHub Wiki 的來源檔

這個目錄是 **GitHub Wiki 頁面的版控來源**。
Wiki 本身是一個獨立的 git repo（`InterSubMod.wiki.git`），內容不會隨主 repo 一起走，
所以在這裡保留一份可審查、可重新產生的來源。

> **不要直接在 GitHub 網頁上編輯 Wiki。**
> 那樣改動不會回流到這裡，下次發布就會被覆蓋。請改這裡的檔案再重新發布。

---

## 檔案對應

| 本地檔案 | Wiki 頁面 | 來源 |
|---|---|---|
| `Home.md` | Wiki 首頁 | 手寫（整體導覽） |
| `_Sidebar.md` | 所有頁面的側邊欄 | 手寫 |
| `System-Overview.md` | System Overview | `docs/explain/11_system-map-overview.standalone.html` |
| `InterSubMod-Engine.md` | InterSubMod Engine | `docs/explain/12_intersubmod-io.standalone.html` |
| `LongLineage-Engine.md` | LongLineage Engine | `docs/explain/13_longlineage-io.standalone.html` |
| `Upstream-and-Data.md` | Upstream & Data | `docs/explain/14_upstream-data.standalone.html` |
| `Analysis-and-Presentation.md` | Analysis & Presentation | `docs/explain/15_python-html-layer.standalone.html` |
| `How-to-Run.md` | How to Run | `docs/explain/16_how-to-run.standalone.html` |

---

## 🔴 發布前必須先滿足的前提

Wiki 的圖片是用 **raw URL** 引用主 repo 的檔案：

```
https://raw.githubusercontent.com/liaoyoyo/InterSubMod/<IMMUTABLE_40HEX_SOURCE_COMMIT>/docs/images/<name>.png
```

**所以 `docs/images/` 必須先存在於同一個 immutable source commit，圖才會顯示。**
禁止以會移動的 `main`／`develop` branch 名稱發布；Wiki source receipt 必須記錄同一個
40-hex source commit 與 claim-registry hash。

`.gitignore` 原本有 `*.png` 全域忽略，已加入例外規則
（`!docs/images/*.png`、`!docs/images/*.svg`）讓這些展示圖能入版控。

---

## 發布步驟

### 第一次發布

GitHub 的 Wiki repo **要先在網頁上建立第一個頁面才會存在**，否則 clone 會失敗。

1. 到 `https://github.com/liaoyoyo/InterSubMod/wiki`
2. 點 **Create the first page**，隨便存一個內容（等下會被覆蓋）
3. 接著執行下面的指令

### 每次發布（建議用腳本）

```bash
bash scripts/publish_wiki.sh \
  --source-commit <IMMUTABLE_40HEX_SOURCE_COMMIT> \
  --dry-run \
  --receipt /tmp/intersubmod-wiki-receipt.json

# 只有 release candidate 的 publication gate 全部通過後，才可另加 --push。
```

腳本必須從指定 commit materialize `docs/wiki/`，並擋掉最常見的失敗：source commit
不是 40-hex／圖片不在同一 commit／claim receipt 不一致／Wiki repo 尚未初始化。

不提供從 dirty worktree 手動複製的發布備援；若腳本失敗，應修復腳本或維持
`publication_state=BLOCKED`，不可繞過 source-commit 與 receipt gate。

推完後到 `https://github.com/liaoyoyo/InterSubMod/wiki` 確認：
**每一頁的圖都有顯示**、側邊欄有出現、頁面之間的連結都能點。

---

## 修改流程

Wiki 內容的真實來源是 `docs/explain/` 的 HTML 頁面。正確的修改順序是：

```
改 docs/explain/*.standalone.html
  → python3 tools/explain_page_qa.py docs/explain/*.standalone.html   # 必須 0 FAIL
  → python3 tools/extract_svg_for_github.py                            # 圖有變才需要
  → 更新 docs/wiki/ 對應的 .md
  → 依上面的步驟發布
```

**不要只改 Wiki 而不改上游 HTML** —— 兩邊會漂移，而 HTML 版才是給教授看的完整版。

---

## 內容紀律

這些頁面對外公開，且描述的是研究方法與能力邊界，因此：

- **一個數字都不准憑記憶寫。** 所有數字必須能在來源 HTML 或原始資料檔中 grep 到。
- **警語與限制不可為了精簡而刪。**「已知缺口」「陷阱」「誠實標註」是這份文件最有價值的部分。
- **不可宣稱超出 `claim_boundary` 的能力**（例如「confirmed subclone」）。
  可主張與禁止主張的完整清單見 System Overview 頁。

本目錄的頁面在產生時已經過一次獨立的數字比對驗證（逐一 grep 回來源）。
日後手動編輯時請維持同樣標準。
