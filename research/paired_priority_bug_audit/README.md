<!--
建立時間: 2026-05-20
類型: directory deprecation header
觸發: 4-agent multi-agent audit (5/19) + user verify (5/20) 揭露 path 命名誤導
-->

# Directory Deprecation Header — paired_priority_bug_audit

> ⚠️ **Path 命名 alias 通知（2026-05-20 加入）**
>
> 本目錄 **`research/paired_priority_bug_audit/`** path 名含 "paired" 為**歷史命名沿用**，但目錄內容實質為：
>
> **ClairS-TO + LongPhase-TO TO mode pipeline 的 phasing priority bug audit**（**全程 tumor-only mode**，無 paired data 參與）。

## Canonical naming（docs / reports / 新 cross-reference 用）

**`clairs_to_phasing_audit`**

對應內容：
- ClairS-TO ssrs model output (`/big7_disk/.../canonical/HCC1395/to_pileup/` = canonical `clairs_to_ssrs`)
- LongPhase-TO V3F / V5 / V6 binary chain (`longphase-to-mod` fork)
- HP priority bug `17.3:1 → 1.85:1` 修補機制 audit
- phaseC genome three-way comparison（12 ISM runs）
- phaseD 5-sample expansion（V6 cross-sample）
- V6 caller F1 verification（檔案 09_V6_caller_F1_verification.md）

## 為何 path 不 rename

- git history 保留（既有 commit history 全引用此 path）
- 跨檔 reference 數十處引用 `research/paired_priority_bug_audit/`，rename 風險 broken link
- 對應 `~/.claude/plans/tender-pondering-blossom.md` plan v2 R-MENTAL-DRIFT 紀律：重大 path rename 需 48hr cooling-off + plan v3 規劃
- Alias 策略已滿足 docs canonical naming + 0 broken link risk

## Naming 對應表

| Filesystem path (legacy) | Canonical name (docs) | 含義 |
|--------------------------|---------------------|------|
| `research/paired_priority_bug_audit/` | `clairs_to_phasing_audit` | ClairS-TO + LongPhase-TO phasing audit |
| `/big7/.../canonical/HCC1395/to_pileup/` | `clairs_to_ssrs` | ClairS-TO v0.3.0 ssrs model output |

## 5/19-5/20 framing 修正紀錄

| Item | Status |
|------|--------|
| 4-agent multi-agent audit 揭露 framing 錯誤 | 2026-05-19 |
| User explicit verify「baseline 與 V6 都是 TO tag.bam」 | 2026-05-20 |
| 5 file framing patch `paired-pileup` → `clairs_to_ssrs` | 2026-05-20 |
| 真實 F1=0.7166 source verify（ClairS-TO ssrs @ 0.93 purity vs SEQC2）| 2026-05-20 |
| Alias 策略採用（path 不 rename）| 2026-05-20 |

## Cross-reference

- `InterSubMod/docs/data_specs/20260411_工作區命名與目錄結構_01.md` §1.5 Canonical Naming Alias Table
- `InterSubMod/research/methyl_augmented_filter_phase2/cycle2/cycle2_NEGATIVE_postmortem.md` §5 Lessons & Action Items
- `InterSubMod/research/autoresearch/evidence_ledger.jsonl` cycle_id `20260520_caller_mode_framing_correction`
- `InterSubMod/research/paired_priority_bug_audit/09_V6_caller_F1_verification.md` F1=0.7166 真實 source
- KB `/big8_disk/liaoyoyo2001/Knowledge/05_tools/variant-callers.md` ClairS-TO 官方 model name (ss/ssrs)

---

**本檔目的**：未來任何讀此目錄的人（新 collaborator / paper reviewer / 未來 self / AI agent）立即明白 path 名 vs 實質內容的 mismatch，並對應到 canonical naming。
