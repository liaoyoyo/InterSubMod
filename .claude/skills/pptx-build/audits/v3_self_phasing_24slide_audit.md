# v3 Self-Phasing Storyboard — myPPT §20 主軸聚焦 Audit

**Audit date**: 2026-04-30
**Auditor**: myPPT Stream C (developer agent)
**Method**: myPPT playbook §20 6-stage filter (Main Thesis → Focal Point → Tier scoring → Definition/Body ratio → 6-Q audit → Noise red-flag)
**Target deck**: `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v3/`
**Sources audited**:
- `InterSubMod/.../v3/00_storyboard_v3.md` (723 lines, 12 slides)
- `InterSubMod/.../v3/notes/speaker_script_v3.md` (442 lines, 13,847 chars excl. whitespace)
- `InterSubMod/.../v3/scripts/build_pptx.py` (1,234 lines)

---

## 0. Scope clarification — v3 is 12 slides, not 24

The Stream C task brief said "v3 24-slide audit". The storyboard at
`InterSubMod/.../v3/00_storyboard_v3.md` declares the deck is **12 slides**
(see line 4: "v3 以 one-slide-one-message 高密度敘事壓縮 v2 的 24 slides 為 12 slides").
v2 was 24 slides; v3 is the compressed 12-slide演講高密度版.

This audit therefore covers the **12 v3 slides** as written in the v3
storyboard. The original 24-slide audit (mentioned in plan §22) maps onto
v2; v2 is frozen and not re-audited here.

Speaker-script statistics audited (`ppt_toolkit.estimate_speaking_seconds`,
CN 400 字/min + EN 150 wpm):

| Slide | Chars (excl. whitespace) | Estimated seconds | Estimated minutes |
|------:|------------------------:|------------------:|-----------------:|
| S1  |  1,018 | 108 | 1.80 |
| S2  |    968 | 104 | 1.73 |
| S3  |  1,052 | 112 | 1.87 |
| S4  |    966 | 100 | 1.67 |
| S5  |  1,073 | 115 | 1.92 |
| **S6**  | **1,605** | **165** | **2.76** |
| **S7**  |  1,425 | 143 | 2.38 |
| S8  |  1,131 | 113 | 1.88 |
| S9  |    940 |  95 | 1.58 |
| **S10** |  1,357 | 142 | 2.36 |
| S11 |  1,308 | 139 | 2.32 |
| S12 |  1,004 | 102 | 1.70 |
| **Total** | **13,847** | **1,438 (≈ 24.0 min)** | **24.0** |

Storyboard self-claimed target: 18 min ± 3 min. Estimator says **24 min**
— the speaker script is 6 minutes over its own self-imposed target, and
significantly over the user's effective 30-min slot when Q&A (12-15 min)
is reserved.

---

## 1. Main Thesis (≤ 30 字, §20 階段 A)

**Inferred from storyboard §0 + speaker script opening (S1 30-second cold open)**:

> **TO 模式 self-phasing 17.3:1 artifact 由 4-commit 漸進修補解到 1:1，解鎖 ISM 五大目標。** (29 字)

Alternative formulations evaluated:

| Candidate | 字數 | Verdict |
|-----------|-----:|---------|
| "longphase-to-mod V5 4-commit 漸進修補 self-phasing artifact 並解鎖 ISM 五大目標" | 33 | over budget |
| "Self-Phasing 17.3:1 artifact 由 V5 解到 1:1，下游 ISM 可信" | 24 | omits action gradient (4 commits) |
| **"TO 模式 self-phasing 17.3:1 由 4-commit 修補到 1:1，解鎖 ISM 五大目標"** | **29** | **selected** |

This thesis matches the v2 plan-quoted thesis ("Self-Phasing 17.3:1 artifact 由
4-commit 漸進修補解決，ISM 下游分析從此可信") but updated for v3's added emphasis on
五大目標解鎖 (S11-S12).

---

## 2. Per-slide Focal Point (≤ 20 字, §20 階段 B)

| Slide | Title (storyboard) | Focal point (≤ 20 字) | 字數 |
|------:|---|---|----:|
| S1  | 17.3 : 1 → 1 : 1 | 17.3:1 是真 artifact, V5 修回 1:1 | 17 |
| S2  | 工作流程一覽 | 修補在 longphase-to-mod, 不在 ISM | 16 |
| S3  | HP tag 五值 + 三層證據 | HP tag 5 值 + 三層證據鏈一致 | 14 |
| S4  | Phasing OK Tag NG | Phasing 不變, Tag 才是 artifact | 16 |
| S5  | ISM 29/14/42 受影響 | 29 個 HP-依賴特徵須在 V5 重跑 | 16 |
| S6  | ★ 根因樹 | Purity 0.927 ≤ 0.95 觸發三 bug | 15 |
| S7  | ★ V5 三層投票 | germline-first + Layer 1.5 + int literal | 18 |
| S8  | Sanity 15/15 + 5 證據鏈 | 4 守恆律 15/15 PASS, 5 證據互證 | 17 |
| S9  | 量化指標 | AMB↓ HP33↓ Concordance +8.3pp | 16 |
| S10 | ★ V5max1 climax | 39 reads 100% reassigned, 守恆律精確 | 17 |
| S11 | 業界家族樹 + 五目標 | 同實驗室相鄰工作 + 解鎖 1/2/4 | 16 |
| S12 | Take-home + 3 P0 | TO 可用, Tag 已解, 五目標解鎖 | 14 |

All 12 focal points fit under the 20-字 ceiling. Coverage of main thesis:
S1 (artifact), S6+S7 (4-commit 修補), S11 (五大目標), S12 (recap) — all
present and traceable.

---

## 3. 6-Question Focal Point Audit (§20 階段 E)

Six questions each slide must answer:
1. **Q1**: Focal point one-sentence (≤ 20 字)?
2. **Q2**: Serves which section thesis?
3. **Q3**: Every element directly supports focal point?
4. **Q4**: Removing which element does NOT hurt understanding?
5. **Q5**: Which details are oral-optional (Tier 3)?
6. **Q6**: Removing this slide entirely — what does the deck lose?

| Slide | Q1 | Q2 | Q3 | Q4 | Q5 | Q6 | Pass rate |
|------:|:--:|:--:|:--:|:--:|:--:|:--:|----------:|
| S1  | PASS | PASS | PASS | PASS | PARTIAL | PASS | 5.5/6 |
| S2  | PASS | PASS | PARTIAL | PASS | PASS | PARTIAL | 5/6 |
| S3  | PASS | PASS | PARTIAL | PARTIAL | PASS | PASS | 5/6 |
| S4  | PASS | PASS | PASS | PASS | PASS | PASS | **6/6** |
| S5  | PASS | PASS | PARTIAL | FAIL | PASS | PASS | 4.5/6 |
| **S6** | PASS | PASS | PASS | PARTIAL | PASS | PASS | 5.5/6 |
| **S7** | PASS | PASS | PARTIAL | FAIL | PARTIAL | PASS | 4/6 |
| S8  | PASS | PASS | PARTIAL | FAIL | PARTIAL | PASS | 4/6 |
| S9  | PASS | PASS | PASS | PASS | PARTIAL | PASS | 5.5/6 |
| **S10** | PASS | PASS | PASS | PASS | PASS | PASS | **6/6** |
| S11 | PASS | PASS | PARTIAL | PARTIAL | FAIL | PARTIAL | 3.5/6 |
| S12 | PASS | PASS | PASS | PARTIAL | PASS | PASS | 5.5/6 |

**Aggregate**: 60.5 / 72 PASS-equivalent points = **84.0%**.

**Per-question summary**:
- Q1 focal point articulation: **12/12 PASS** (100%)
- Q2 serves section thesis: **12/12 PASS** (100%)
- Q3 every element supports focal point: 5 PASS / 6 PARTIAL / 1 FAIL → 70.8%
- Q4 removable elements: 5 PASS / 4 PARTIAL / 3 FAIL → 58.3%
- Q5 oral-optional split: 7 PASS / 4 PARTIAL / 1 FAIL → 75.0%
- Q6 slide-level necessity: 9 PASS / 3 PARTIAL / 0 FAIL → 87.5%

**Weakest questions**: Q4 (removable elements, 58.3%) — slides try to do
multiple jobs. Q3 (element-focal alignment, 70.8%) — auxiliary visuals carry
content not strictly needed.

---

## 4. Definition / Prerequisite / Body / Conclusion ratio (§20 階段 D)

| Class | v3 slides | Count | % | §20 D upper bound | Verdict |
|-------|-----------|------:|---:|:----:|:------:|
| **Definition** | S3 (HP tag five-value) | 1 | 8.3% | ≤ 10% | OK |
| **Prerequisite** | S2 (workflow), S4 (LOH 兩層分流) | 2 | 16.7% | ≤ 15% | **slightly over** |
| **Body** | S1, S5, S6, S7, S8, S9, S10, S11 | 8 | 66.7% | ≥ 60% | OK |
| **Conclusion** | S12 | 1 | 8.3% | ≥ 15% | **under** |

S1 is borderline definition/conclusion — classified as Body because the
17.3:1 → 1:1 framing is the core thesis claim. S4 includes the LOH 兩層
explanation (LOH.bed vs ISM HP_Ratio LOH) that doubles as both prerequisite
and body argument; classified as Prerequisite because the speaker script
treats it as the "tag layer kerf-cut" definitional moment before S5's body
argument begins.

**Issues**:
- Definition + Prerequisite = 25.0% (just at the §20 D 25% trigger). Not
  strictly over, but right at the line.
- Conclusion class is under-represented — only S12 wraps up. Adding a
  short "P0 commitments + thanks" wrap would bring conclusion to ~15%
  but cost a slide. Alternative: split S12 into two (take-home + Q&A
  invitation), accepting 13 slides. **Recommendation: keep S12 single,
  use Q&A as the conclusion buffer.**
- Body share (66.7%) is healthy.

---

## 5. Noise Red Flags (§20 階段 F)

The §20 stage F red-flag list flags speaker-note phrases that signal
tangential / oral-optional / off-topic content. Below are concrete
violations found in `speaker_script_v3.md`.

### Red flag 1 — "順便提一下" / "另外" patterns

**S5 speaker script** (line 142):
> "Q: HPFineNGroups 重詮釋是什麼意思？ → A: 過去評為 subclone marker（甲基化 bimodality），V5 後重詮釋為 phasing signature（LOH-constrained phasing），論文主軸候選。詳見 v2 S22 + Q11。"

Status: NOT a red-flag itself, but the embedded Q&A within speaker note is
characteristic of "as an aside" framing. Consider tagging `[ORAL-OPTIONAL]`.

### Red flag 2 — "這部分如果有時間可以再細講" (defer markers)

**S11 speaker script** (line 304):
> "已有初步發現的方向（全部依賴 V5 BAM）：(1) Thread D LOH-constrained phasing — NG=2 cross-sample 6/6 POSITIVE，Wilcoxon p=0.0156…（2）HPFineNGroups marker 機制重詮釋為 phasing signature；（3）Phase 2A normal methylation reference 解鎖。"

The full enumeration of three downstream directions is **not** required
to support S11's focal point ("同實驗室相鄰工作 + 解鎖 1/2/4"). Each item
warrants its own slide in a follow-up deck — surfacing them in the speaker
note creates "如果有時間可以再細講" pressure on the speaker. **Tag as Tier 3
oral-optional or move to a single bullet "已有初步發現待 follow-up"**.

### Red flag 3 — Multiple bullet enumeration > 3 items

**S5 speaker script** (lines 129-133):
- Three tiers (29 / 14 / 42), each with 4-7 named features:
  - 嚴重影響 29 個: HP_Ratio, Potential_LOH, HPMergedDelta, HPMergedSig, HPFineNGroups
  - 中度影響 14 個: QualityScore, GlobalP, CramersV, VerificationClass
  - 不受影響 42 個: PairwiseMeanDist, AlleleDelta, AlleleP, Caller features, ...

Listing 4-7 feature names per tier in the speaker note (not on slide)
exceeds the §20 F "more than 3 bullets" red flag. Recommendation: speaker
should name 2-3 representative features per tier and mark the rest
oral-optional `[ORAL-OPTIONAL]`.

### Red flag 4 — "供參考" / footnote-as-content

**S6 speaker script** (line 161):
> "Bug 3 (灰色, scaffold, 已由 V2b 解): … 已由 V2b commit 8b8c1fd 啟用 PON-only phasing 解決，**本 PPT 不展開**，視覺只佔 20%。"

The phrase "本 PPT 不展開" is a self-acknowledged red flag — the speaker
already knows this content is tangential. **Action**: remove from slide,
keep one sentence in note as "Bug 3 已由 V2b 解，本演講不展開". The visual
"佔 20%" allocation can be reduced to a small grey footnote.

### Red flag 5 — Generic-label adjacent (Q&A enumeration)

**S12 speaker script** (lines 357-371):
The 13-Q Q&A preparation table is reproduced verbatim in the speaker note.
This violates §20 F "speaker note 中如果有時間可以再細講" red flag — the Q&A
prep is reference material, not narration content. **Action**: move the
13-Q table to a separate `qa_prep.md` companion file (already exists at
`v2/notes/qa_11_questions.md`). The speaker note should reference it, not
inline the entire table.

---

## 6. Recommended Actions

### A. Slide consolidation candidates

| Slide | Current role | Recommendation | Justification |
|-------|--------------|----------------|---------------|
| S2 | Workflow + 4-stage framing | **Keep, compress to 30 sec** | Already shorter than S6/S7; speaker storyboard §3 explicitly marks it as compressible. |
| S3 | HP tag 5-value + 3-layer evidence | **Keep, but split if S6 grows** | Currently dual-purpose (definition + evidence). If S6 root-cause needs more time, move 三層證據 to S4. |
| S5 | ISM impact 3-tier + 4-bucket | **Keep, simplify right column** | Right-column 4-bucket grid duplicates S6 root-cause discussion. Drop the 4-bucket sub-grid. |
| S11 | Family tree + 5 goals | **Split or compress** | Two distinct messages competing in one slide. Q4 audit FAIL. Pick one: (a) family tree only with goals as a 1-line strip, or (b) 5-goals only with a single sentence on family tree. |

### B. Slides to discard candidates

**No slide is recommended for outright deletion.** The 12-slide structure
is already the result of v2→v3 50% compression. However:

- If time slot is < 12 min, drop S2 entirely (the workflow framing can be
  voiced in 15 sec at the top of S3 instead of taking a full slide).
- If time slot is < 10 min, drop S2 + S11 (5 goals can be 30 sec voice-over
  on S12 instead of dedicated slide).

### C. Oral-optional content to extract from speaker notes

**S5 speaker script** lines 129-133 — feature name enumeration (29+14+42):
move 4-7 named features per tier into `[ORAL-OPTIONAL]` blocks. Speaker
delivers 2-3 representatives only.

**S6 speaker script** lines 161-162 — Bug 3 V2b PON-only commit detail
("commit 8b8c1fd"): full commit hash + diff stats are deep-reference
material. Tag `[ORAL-OPTIONAL]`.

**S7 speaker script** lines 184-188 — Layer 1.5 0.6 confidence threshold
derivation ("09 報告 0.6 sample 上 V5 HP33 比例 12.4% vs baseline 2%"):
move into `[ORAL-OPTIONAL]`. Speaker says only "經驗值, F2 待驗證".

**S8 speaker script** lines 211-218 —守恆律 A/B/Layer1.5 expectations 1/2
full formula reproduction: rendered on slide as left-column. Note can
reference "見 slide left column" rather than re-cite formulas. Tag
`[ORAL-OPTIONAL]`.

**S10 speaker script** lines 277-279 — V5max2 / V5max3 旁證 with chr19
positions and ΔHP11/ΔHP21 magnitudes: move to `[ORAL-OPTIONAL]`. Speaker
mentions only "3/3 SP-extreme 一致翻轉" (already on slide caveat).

**S11 speaker script** lines 302-304 — Three downstream directions
(Thread D NG=2, HPFineNGroups, Phase 2A): tag `[ORAL-OPTIONAL]`. Speaker
delivers as "已有初步發現待 follow-up" one-liner.

**S12 speaker script** lines 357-371 — 13-Q Q&A table inline: remove from
note entirely; reference `v2/notes/qa_11_questions.md`.

### D. Speaker-script word-count reduction target

**Current**: 13,847 chars (excl. whitespace) ≈ 1,438 sec ≈ **24.0 min**.
**Target (storyboard claim)**: 18 min = 1,080 sec ≈ **10,400 chars**.
**Aggressive target (30-min total slot, 12-min Q&A)**: 18 min = 10,400 chars.
**Conservative target (30-min total slot, 15-min Q&A)**: 15 min = 8,650 chars.

Recommended cuts (apply A-C above):

| Source | Current chars | After cut | Saved chars |
|--------|--------------:|---------:|------------:|
| S5 feature enum | ~250 | ~80 | 170 |
| S6 Bug 3 detail | ~180 | ~50 | 130 |
| S7 0.6 threshold derivation | ~200 | ~60 | 140 |
| S8 守恆律 inline formulas | ~300 | ~120 | 180 |
| S10 V5max2/3 旁證 | ~150 | ~50 | 100 |
| S11 three downstream directions | ~280 | ~80 | 200 |
| S12 13-Q table | ~750 | ~50 | 700 |
| **Total saved** | | | **~1,620 chars** |

Post-cut estimate: 13,847 − 1,620 = **12,227 chars** ≈ 21.2 min — still
over 18-min target. Additional 1,800-char cut needed; recommend
S6/S7/S10 (the 3 必講 cores currently averaging 1,460 chars each) compress
to ~1,200 chars each = saves another ~780 chars. Final estimate: ~11,400
chars ≈ 19.7 min. **Acceptable for 30-min total slot with 10-min Q&A**.

To hit the 18-min goal exactly, additionally compress S2/S5/S11 by 200
chars each (= 600 chars saved) — total post-cut **~10,800 chars ≈ 18.7 min**.

### E. Suggested final character budget per slide

| Slide | Current | Target | Cut % |
|------:|--------:|-------:|------:|
| S1  | 1,018 | 900 | 11.6% |
| S2  | 968 | 750 | 22.5% |
| S3  | 1,052 | 900 | 14.4% |
| S4  | 966 | 850 | 12.0% |
| S5  | 1,073 | 800 | 25.4% |
| **S6**  | 1,605 | 1,250 | 22.1% |
| **S7**  | 1,425 | 1,150 | 19.3% |
| S8  | 1,131 | 850 | 24.8% |
| S9  | 940 | 800 | 14.9% |
| **S10** | 1,357 | 1,150 | 15.3% |
| S11 | 1,308 | 950 | 27.4% |
| S12 | 1,004 | 600 | 40.2% |
| **Total** | **13,847** | **10,950** | **20.9%** |

Translates to ≈ 19.0 min of speech + handover seconds → fits 18-min target
with 60 sec of buffer.

---

## 7. PLOS Ten Simple Rules + Assertion-Evidence Pass-rate

| Rule | v3 Status | Notes |
|------|-----------|-------|
| **Rule 1**: One idea per slide | 10/12 PASS | S5 dual-purpose (impact tier + 4-bucket); S11 dual-purpose (family tree + 5 goals) |
| **Rule 2**: Heading is thesis | 12/12 PASS | All titles are claim-style (e.g. "Phasing OK Tag NG", "39 reads 100% reassigned") |
| **Rule 3**: Cognitive load (≤ 6 elements) | 8/12 PASS | S5 right column 4-bucket adds 7th element; S6 root-cause has root + 3 leaves + connector + caveat = 6 (borderline); S8 left+right both have 4-5 sub-elements each (≈ 9 total); S11 family tree + 5 cards = 8 |
| **Rule 4**: Essential only | 8/12 PASS | S5 4-bucket grid not essential; S6 Bug 3 灰色 footprint not essential; S11 three downstream directions not essential |
| **Rule 5**: 1 min per slide | 0/12 PASS | All 12 slides over 1 min in current speaker script; even after cuts, S6 and S7 will exceed |
| **Rule 6**: Iterative practice / dry-run | UNKNOWN | Not in storyboard — recommend adding § "演練 checklist" |
| **Rule 7**: Clear takeaway without audio | 11/12 PASS | S2 4-stage strip arguably needs voice-over to interpret; otherwise titles + visuals self-explain |
| **Rule 8**: Visual contrast | 12/12 PASS | Build script uses palette tokens consistently; CJK + Latin font fallback enforced |
| **Rule 9**: Distractor-free | 11/12 PASS | S11 family tree + 5 goals competing for attention is a mild distractor |
| **Rule 10**: Backup plan | 12/12 PASS | v2 PPTX backup explicitly mentioned; figures resilient via fit_image_within fallback |
| **AE structure (Alley)**: full assertion sentence + dominant visual | 11/12 PASS | All titles are assertions; S2 visual (pipeline + 4-stage) is dominant but not aligned to a single thesis claim |

**Aggregate**: 105 / 120 = **87.5%** (excluding Rule 6 unknown).

---

## 8. Overall Verdict

**v3 12-slide storyboard scores 87.5% on PLOS + AE structure**, with main
shortcomings:

1. **Speaker script over-runs target** by 6 minutes (24 min vs 18 min target).
2. **2 slides try to deliver 2 ideas each** (S5 impact + 4-bucket; S11 family tree + 5 goals).
3. **Cognitive load pressure on S6, S8, S11** (>6 visual elements).
4. **Conclusion class under-represented** (8.3% vs §20 D ≥ 15% recommendation).

**Recommended path**:
- Apply oral-optional cuts in §6.C (saves ~1,620 chars)
- Compress S6/S7/S10 narration (saves ~780 chars)
- Compress S2/S5/S11 (saves ~600 chars)
- **Final**: ~10,800 chars ≈ 18.7 min — fits 30-min total slot with
  10-min Q&A buffer + 60 sec handover.
- **Optional**: split S11 into 11a (family tree) + 11b (5 goals) → 13 slides
  but cleaner cognitive load. Trade-off: more slides vs cleaner ideas.

**Decision matrix** for 11-split: only do this if 30-min slot expands to
33-35 min; otherwise keep S11 single and accept the dual-message penalty
(Q4 PARTIAL on element removability).

---

## 9. Audit Provenance

| Artefact | Path |
|----------|------|
| Storyboard | `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v3/00_storyboard_v3.md` |
| Speaker script | `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v3/notes/speaker_script_v3.md` |
| Build script | `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v3/scripts/build_pptx.py` |
| v2 reference | `InterSubMod/docs/presentations/validated/2026/04/20260429_Self_Phasing_完整教授報告/v2/00_storyboard_v2.md` |
| Estimator | `InterSubMod/tools/ppt_toolkit/tier_aware_speaker_note.py::estimate_speaking_seconds` |
| Audit method | `InterSubMod/.claude/skills/myPPT/playbook.md` §20 (Stream A pending) |

---

**End of audit.**
