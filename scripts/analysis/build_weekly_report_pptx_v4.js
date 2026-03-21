#!/usr/bin/env node
/**
 * build_weekly_report_pptx_v4.js
 * Fresh design from scratch using PptxGenJS (Midnight Executive palette)
 * 根據:
 *   - docs/references/manual/20260311_個人PPT設計風格規範_01.md
 *   - docs/references/manual/20260311_研究週報PPTX客製化設定與製作手冊_01.md
 * 來源: docs/reports/validated/2026/03/20260311_研究主線週報_..._01.md
 */

const PPTXGENJS_PATH = "/home/liaoyoyo2001/.nvm/versions/node/v18.20.8/lib/node_modules/pptxgenjs";
const pptxgen = require(PPTXGENJS_PATH);

const OUTPUT_PATH = "/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/draft/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_v04.pptx";

// ── COLOR PALETTE (Midnight Executive adapted for research) ──────────────
const C = {
  dark:     "1A237E",   // deep indigo-navy for dark slides
  dark2:    "141B4D",   // darker navy accent
  navy:     "1E2761",   // navy text / headers
  mid:      "1565C0",   // mid blue
  teal:     "00897B",   // teal accent (P1)
  amber:    "F9A825",   // amber highlight (P2, warnings)
  indigo:   "5E35B1",   // indigo (P3)
  bg:       "F0F5FF",   // very light blue-white for content slides
  bg_alt:   "E3ECFF",   // slightly darker alt row
  text:     "0D1B2A",   // near black
  white:    "FFFFFF",
  muted:    "546E7A",
  border:   "90A4AE",
  green:    "2E7D32",
  red:      "C62828",
  light:    "ECEFF1",   // light for card backgrounds on dark slides
};

// ── FONTS ─────────────────────────────────────────────────────────────────
const F = { zh: "Noto Sans TC", en: "Arial" };

// ── SLIDE DIMENSIONS (LAYOUT_16x9) ────────────────────────────────────────
const W = 10;       // width inches
const H = 5.625;    // height inches

// ─────────────────────────────────────────────────────────────────────────
// HELPERS
// ─────────────────────────────────────────────────────────────────────────

/** Add title + english subtitle + separator to a content (light) slide */
function addHeader(pres, slide, zhTitle, enTitle) {
  // Left accent bar
  slide.addShape(pres.ShapeType.rect, {
    x: 0, y: 0, w: 0.1, h: H,
    fill: { color: C.navy },
    line: { color: C.navy },
  });
  // Chinese title
  slide.addText(zhTitle, {
    x: 0.25, y: 0.16, w: 9.6, h: 0.62,
    fontSize: 24, fontFace: F.zh, bold: true,
    color: C.navy, margin: 0,
  });
  // English subtitle
  slide.addText(enTitle, {
    x: 0.25, y: 0.76, w: 9.6, h: 0.22,
    fontSize: 10.5, fontFace: F.en,
    color: C.muted, margin: 0,
  });
  // Separator line
  slide.addShape(pres.ShapeType.line, {
    x: 0.25, y: 1.01, w: 9.5, h: 0,
    line: { color: C.border, width: 0.75 },
  });
}

/** Header cell */
function th(text) {
  return {
    text,
    options: {
      fill: { color: C.navy },
      color: C.white,
      bold: true,
      fontSize: 10,
      fontFace: F.zh,
      align: "center",
      valign: "middle",
    },
  };
}

/** Data cell */
function td(text, { bold = false, color, fill, align = "left" } = {}) {
  const opts = {
    fontSize: 10,
    fontFace: F.zh,
    color: color || C.text,
    bold,
    align,
    valign: "middle",
  };
  if (fill) opts.fill = fill;
  return { text: String(text), options: opts };
}

/** Numeric cell (Arial, center) */
function tn(text, { bold = false, color, fill } = {}) {
  const opts = {
    fontSize: 10,
    fontFace: F.en,
    color: color || C.text,
    bold,
    align: "center",
    valign: "middle",
  };
  if (fill) opts.fill = fill;
  return { text: String(text), options: opts };
}

/** Alternating row fill */
function rowFill(i) {
  return i % 2 === 1 ? { color: C.bg_alt } : undefined;
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 1 – COVER
// ─────────────────────────────────────────────────────────────────────────
function buildCover(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.dark };

  // Left teal accent bar
  slide.addShape(pres.ShapeType.rect, {
    x: 0, y: 0, w: 0.55, h: H,
    fill: { color: C.teal }, line: { color: C.teal },
  });

  // Vertical label on bar
  slide.addText("週報", {
    x: 0, y: 2.0, w: 0.55, h: 1.5,
    fontSize: 12, fontFace: F.zh, bold: true,
    color: C.white, rotate: 270, align: "center", valign: "middle", margin: 0,
  });

  // Main title
  slide.addText("InterSubMod 研究主線整合週報", {
    x: 0.8, y: 1.0, w: 8.8, h: 1.1,
    fontSize: 34, fontFace: F.zh, bold: true,
    color: C.white, margin: 0, wrap: true,
  });

  // English subtitle
  slide.addText("InterSubMod Weekly Research Report", {
    x: 0.8, y: 2.15, w: 8, h: 0.4,
    fontSize: 16, fontFace: F.en,
    color: "CADCFC", margin: 0,
  });

  // Date badge
  slide.addShape(pres.ShapeType.rect, {
    x: 0.8, y: 2.7, w: 4.8, h: 0.48,
    fill: { color: C.teal }, line: { color: C.teal },
  });
  slide.addText("2026-03-05 至 2026-03-11　Phase 2 Annotation 整合", {
    x: 0.85, y: 2.72, w: 4.7, h: 0.44,
    fontSize: 13, fontFace: F.zh, bold: true,
    color: C.white, margin: 0, valign: "middle",
  });

  // Researcher + date
  slide.addText("廖子游　2026-03-11", {
    x: 0.8, y: 5.1, w: 4, h: 0.3,
    fontSize: 12, fontFace: F.zh,
    color: "90A4AE", margin: 0,
  });

  slide.addNotes(
    "【封面】\n" +
    "這頁目的：提供報告的基本定位與識別資訊。\n" +
    "口頭重點：本週報告覆蓋 3/5 至 3/11，核心新增是 phase 2 benchmark 與 annotation layer 回接，\n" +
    "最終確立三層框架（caller-first → methylation-support → annotation/QC）。\n" +
    "來源：docs/reports/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 2 – AGENDA
// ─────────────────────────────────────────────────────────────────────────
function buildAgenda(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.bg };
  addHeader(pres, slide, "本週報告大綱", "Weekly Report Agenda — 14 Slides");

  const items = [
    ["03", "核心結論：三層框架已收斂"],
    ["04", "研究背景與場景定義"],
    ["05", "本週主要目標"],
    ["06", "三層框架：特徵定位表"],
    ["07", "Paired Pure Benchmark"],
    ["08", "Candidate-Specific Rescue"],
    ["09", "Coverage Ceiling 解讀"],
    ["10", "Phase 2 Raw Benchmark"],
    ["11", "Annotation Layer 最佳策略"],
    ["12", "待補齊事項與限制"],
    ["13", "下一步建議"],
    ["14", "附錄路徑"],
  ];

  const col1 = items.slice(0, 6);
  const col2 = items.slice(6);
  const startY = 1.15;
  const itemH = 0.62;

  [col1, col2].forEach((col, ci) => {
    const xBase = ci === 0 ? 0.35 : 5.15;
    col.forEach((item, i) => {
      const y = startY + i * itemH;
      // Number badge
      slide.addShape(pres.ShapeType.rect, {
        x: xBase, y: y + 0.04, w: 0.44, h: 0.4,
        fill: { color: C.navy }, line: { color: C.navy },
      });
      slide.addText(item[0], {
        x: xBase, y: y + 0.04, w: 0.44, h: 0.4,
        fontSize: 12, fontFace: F.en, bold: true,
        color: C.white, align: "center", valign: "middle", margin: 0,
      });
      // Text
      slide.addText(item[1], {
        x: xBase + 0.52, y: y + 0.07, w: 4.2, h: 0.36,
        fontSize: 13, fontFace: F.zh,
        color: C.text, margin: 0,
      });
    });
  });

  slide.addNotes(
    "【大綱】\n" +
    "共 14 頁。從核心結論開始，依序鋪陳背景、目標、方法框架、數據證據（4 張），" +
    "最後是待補齊事項、下一步與附錄。"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 3 – OPENING THESIS (dark)
// ─────────────────────────────────────────────────────────────────────────
function buildOpeningThesis(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.dark };

  // Title
  slide.addText("核心結論：三層框架已收斂", {
    x: 0.5, y: 0.18, w: 9, h: 0.68,
    fontSize: 26, fontFace: F.zh, bold: true,
    color: C.white, margin: 0,
  });
  slide.addText("Key Findings: Three-Tier Framework Established", {
    x: 0.5, y: 0.84, w: 9, h: 0.24,
    fontSize: 11, fontFace: F.en,
    color: "90A4AE", margin: 0,
  });

  // 3 tier cards
  const tiers = [
    { tier: "第一層", en: "Caller-First", color: C.teal, body: "GQ / QUAL gate\n最穩定、跨 dataset 一致\n→ 最佳 rescue 的守門條件" },
    { tier: "第二層", en: "Methylation-Support", color: "1E88E5", body: "Quality_Score / Pairwise\nagreement_positive\n→ TO 下有真實正向 support" },
    { tier: "第三層", en: "Annotation / QC / Triage", color: C.indigo, body: "hp_assign / artifact-veto\ndataset-aware 標記\n→ 制度化分析層，不做 hard rule" },
  ];

  tiers.forEach((t, i) => {
    const x = 0.4 + i * 3.12;
    const y = 1.2;
    // Card background
    slide.addShape(pres.ShapeType.rect, {
      x, y, w: 2.95, h: 2.85,
      fill: { color: "1A2A5E" },
      line: { color: t.color, width: 2 },
    });
    // Top color bar
    slide.addShape(pres.ShapeType.rect, {
      x, y, w: 2.95, h: 0.12,
      fill: { color: t.color }, line: { color: t.color },
    });
    // Tier label
    slide.addText(t.tier, {
      x: x + 0.14, y: y + 0.2, w: 2.7, h: 0.38,
      fontSize: 15, fontFace: F.zh, bold: true,
      color: t.color, margin: 0,
    });
    slide.addText(t.en, {
      x: x + 0.14, y: y + 0.57, w: 2.7, h: 0.24,
      fontSize: 9.5, fontFace: F.en,
      color: "78909C", margin: 0,
    });
    // Body
    slide.addText(t.body, {
      x: x + 0.14, y: y + 0.86, w: 2.7, h: 1.8,
      fontSize: 12, fontFace: F.zh,
      color: C.white, margin: 0, wrap: true,
    });
  });

  // Bottom takeaway bar
  slide.addShape(pres.ShapeType.rect, {
    x: 0, y: 4.35, w: W, h: 1.28,
    fill: { color: "0D3349" }, line: { color: "0D3349" },
  });
  // Teal left accent
  slide.addShape(pres.ShapeType.rect, {
    x: 0, y: 4.35, w: 0.12, h: 1.28,
    fill: { color: C.teal }, line: { color: C.teal },
  });
  slide.addText(
    "甲基訊號在 TO 場景下有真實 support 價值，但最穩定的 rescue 仍在 caller 端（GQ/QUAL gate）。\n" +
    "annotation layer 已定型為第三層制度化框架。下一步是 annotation export，不是改 hard filter。",
    {
      x: 0.25, y: 4.37, w: 9.6, h: 1.22,
      fontSize: 12.5, fontFace: F.zh, bold: true,
      color: C.white, margin: 0, wrap: true, valign: "middle",
    }
  );

  slide.addNotes(
    "【核心結論】\n" +
    "這頁目的：讓聽眾在 5 秒內知道本週最終定論。\n" +
    "三層框架：\n" +
    "  第一層 GQ/QUAL — 最穩定、跨 dataset 最一致（paired / TO 均成立）\n" +
    "  第二層甲基 support — TO 下成立但不超過 caller-only，paired 因 coverage ceiling 結論保守\n" +
    "  第三層 annotation — 制度化標記，讓分析層訊號可被後續工作使用\n" +
    "4 個核心證據：paired pure F1 增益、candidate-specific rescue、coverage ceiling、annotation best policy\n" +
    "關鍵限制：GQ 仍是 caller 自身的信心分數，用它 rescue 本質上是放寬 caller threshold，非 InterSubMod 獨立甲基貢獻。"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 4 – BACKGROUND
// ─────────────────────────────────────────────────────────────────────────
function buildBackground(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.bg };
  addHeader(pres, slide, "研究背景：InterSubMod 的核心問題", "Research Context: What InterSubMod Addresses");

  // ─ Left column: problem ─
  slide.addText("InterSubMod 整合三種資訊來源", {
    x: 0.3, y: 1.15, w: 4.5, h: 0.3,
    fontSize: 13, fontFace: F.zh, bold: true, color: C.navy, margin: 0,
  });

  const inputs = [
    "somatic caller 的輸出（PASS 與 non-PASS 邊緣）",
    "LongPhase 的 phase / haplotag 結果",
    "read-level 甲基距離、label / cluster 訊號",
  ];
  inputs.forEach((s, i) => {
    slide.addShape(pres.ShapeType.rect, {
      x: 0.3, y: 1.5 + i * 0.58, w: 0.32, h: 0.32,
      fill: { color: C.teal }, line: { color: C.teal },
    });
    slide.addText(String(i + 1), {
      x: 0.3, y: 1.5 + i * 0.58, w: 0.32, h: 0.32,
      fontSize: 11, fontFace: F.en, bold: true,
      color: C.white, align: "center", valign: "middle", margin: 0,
    });
    slide.addText(s, {
      x: 0.68, y: 1.52 + i * 0.58, w: 4.05, h: 0.28,
      fontSize: 12, fontFace: F.zh, color: C.text, margin: 0, wrap: true,
    });
  });

  slide.addText("→ 回答兩件事：", {
    x: 0.3, y: 3.25, w: 4.5, h: 0.28,
    fontSize: 12, fontFace: F.zh, color: C.muted, margin: 0,
  });
  const answers = [
    "哪些 caller 邊緣位點其實是 TP，值得救回",
    "哪些位點雖有結構訊號，但更像 artifact",
  ];
  answers.forEach((a, i) => {
    slide.addText([
      { text: "→ ", options: { color: C.amber, bold: true } },
      { text: a },
    ], {
      x: 0.35, y: 3.55 + i * 0.42, w: 4.4, h: 0.36,
      fontSize: 12, fontFace: F.zh, color: C.text, margin: 0, wrap: true,
    });
  });

  // Vertical divider
  slide.addShape(pres.ShapeType.line, {
    x: 5.05, y: 1.12, w: 0, h: 4.2,
    line: { color: C.border, width: 0.75 },
  });

  // ─ Right column: 3 scenarios ─
  slide.addText("本週三個核心場景", {
    x: 5.2, y: 1.15, w: 4.5, h: 0.3,
    fontSize: 13, fontFace: F.zh, bold: true, color: C.navy, margin: 0,
  });

  const scenarios = [
    { name: "paired tumor-normal", zh: "配對腫瘤/正常樣本", role: "主驗證場景（最穩定）", color: C.teal },
    { name: "tumor-only (TO)", zh: "僅腫瘤樣本", role: "驗證 borderline rescue（最有甲基空間）", color: "1E88E5" },
    { name: "mixed purity", zh: "混和純度樣本", role: "已延後，不是本週主證據", color: C.muted },
  ];

  scenarios.forEach((s, i) => {
    const y = 1.52 + i * 1.04;
    slide.addShape(pres.ShapeType.rect, {
      x: 5.2, y, w: 4.5, h: 0.92,
      fill: { color: i === 2 ? "ECEFF1" : "E8F5FE" },
      line: { color: s.color, width: 1.5 },
    });
    slide.addText(s.name, {
      x: 5.3, y: y + 0.06, w: 4.3, h: 0.28,
      fontSize: 12.5, fontFace: F.en, bold: true,
      color: i === 2 ? C.muted : C.navy, margin: 0,
    });
    slide.addText(s.zh + "　" + s.role, {
      x: 5.3, y: y + 0.38, w: 4.3, h: 0.42,
      fontSize: 11, fontFace: F.zh,
      color: i === 2 ? C.muted : C.text, margin: 0, wrap: true,
    });
  });

  slide.addNotes(
    "【研究背景】\n" +
    "InterSubMod 的核心作用是把 somatic caller / phasing / 甲基資訊整合，回答兩個問題：救什麼（TP rescue）、排除什麼（artifact triage）。\n" +
    "三個場景中，paired 是最穩定的主驗證環境，TO 是甲基 rescue 空間最大的場景（coverage 完整），mixed purity 本週延後。"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 5 – GOALS
// ─────────────────────────────────────────────────────────────────────────
function buildGoals(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.bg };
  addHeader(pres, slide, "本週主要目標與完成度", "Weekly Goals & Completion Status");

  const goals = [
    ["paired pure 正式重跑", "5kHz & DORADO 新版欄位 rerun", "100%", "5kHz 有正增益，DORADO 幾乎持平略負"],
    ["TO 主流程打通", "ClairS-TO → LongPhase-TO → InterSubMod", "100%", "TO pipeline 已可重跑"],
    ["candidate-specific rescue", "5kHz TO、DORADO paired、DORADO TO", "100%", "TO 甲基 support 成立，不超過 caller-only"],
    ["phase 2 benchmark", "raw pileup/full 對照、finer interval、orthogonality", "90%", "paired 已有定論，TO 雙模型尚未補齊"],
    ["annotation layer 回接", "phase 2 結果回接 analysis-layer", "90%", "policy 已完成，dashboard 整合尚未完成"],
    ["文件制度化", "Agent / skills / 週報流程", "100%", "可持續用同一套格式追蹤"],
  ];

  const rows = goals.map((g, i) => {
    const is100 = g[2] === "100%";
    const fill = rowFill(i);
    const pctFill = is100 ? { color: "C8E6C9" } : { color: "FFF9C4" };
    const pctColor = is100 ? C.green : "8B6914";
    return [
      td(g[0], { bold: true, fill }),
      td(g[1], { fill }),
      { text: g[2], options: { fill: pctFill, color: pctColor, bold: true, fontSize: 10.5, fontFace: F.en, align: "center", valign: "middle" } },
      td(g[3], { color: C.muted, fill }),
    ];
  });

  slide.addTable(
    [[th("主目標"), th("子目標"), th("完成度"), th("本週結論")], ...rows],
    { x: 0.3, y: 1.12, w: 9.4, colW: [1.8, 2.6, 0.9, 4.1], rowH: 0.52, border: { pt: 0.5, color: C.border } }
  );

  slide.addNotes(
    "【本週目標】\n" +
    "6 個主目標中 4 個 100%，2 個 90%。\n" +
    "90% 的兩項：\n" +
    "  phase 2 benchmark — TO 的 pileup vs full 雙模型對照尚未補齊\n" +
    "  annotation layer 回接 — policy 已完成，但尚未接到 dashboard / round summary\n" +
    "文件制度化 100%：之後可用 Agent/skills 流程持續追蹤。"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 6 – FEATURE ROLES
// ─────────────────────────────────────────────────────────────────────────
function buildFeatureRoles(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.bg };
  addHeader(pres, slide, "三層框架：各特徵定位與使用限制", "Feature Role Definitions in the Three-Tier Framework");

  const features = [
    ["GQ", "caller-first 主規則", "最穩定、跨 dataset 最一致", "不可解讀成 caller 應全面放寬；仍是 caller-side ranking 訊號"],
    ["Quality_Score", "support annotation", "兩個 TO dataset 都有正向 support", "不可直接 hard keep"],
    ["PairwiseMedianDist", "dataset-aware support", "5kHz TO 偏高、DORADO TO 偏低（方向相反）", "不可升級成跨樣本全域規則"],
    ["agreement_positive", "5kHz TO label-support", "比單純 Strong/Subclone 更乾淨的 label 訊號", "不可當所有 dataset 通用 support"],
    ["hp_assign_rate", "phase/QC annotation", "反映 HP 可分派性與相位完整度", "不宜當 truth-level keep 規則"],
    ["low VAF+high AD+low CV", "artifact triage", "在 5kHz 對刪 FP 有實際價值", "不應拿來做 TP rescue 主規則"],
  ];

  const rows = features.map((f, i) => {
    const fill = rowFill(i);
    return [
      td(f[0], { bold: true, fill }),
      td(f[1], { fill }),
      td(f[2], { fill }),
      td(f[3], { color: C.muted, fill }),
    ];
  });

  slide.addTable(
    [[th("特徵"), th("定位"), th("已驗證的意義"), th("不該怎麼用")], ...rows],
    { x: 0.3, y: 1.12, w: 9.4, colW: [1.55, 1.8, 2.85, 3.2], rowH: 0.52, border: { pt: 0.5, color: C.border } }
  );

  slide.addNotes(
    "【特徵定位表】\n" +
    "GQ：最穩定，但本質是 caller 信心分數，用它 rescue = 放寬 caller threshold。\n" +
    "Quality_Score：兩個 TO dataset 都有正向，適合作為 support annotation。\n" +
    "PairwiseMedianDist：5kHz TO >=0.20 有效，DORADO TO <=0.20 有效，方向相反。\n" +
    "→ 這意味著 Pairwise 不是穩定的生物學訊號，而是 dataset-specific 結構/noise 代理。\n" +
    "hp_assign_rate：TP/FP 在 TO 場景幾乎都達 1.0，無鑑別力，降階為 QC 標記。"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 7 – PAIRED PURE BENCHMARK
// ─────────────────────────────────────────────────────────────────────────
function buildPairedPure(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.bg };
  addHeader(pres, slide, "Paired Pure Rerun：兩平台正式 Benchmark", "HCC1395 5kHz / DORADO Paired Pure Benchmark Results");

  const data = [
    ["HCC1395 5kHz paired",   "ClairS",      "29865", "1430", "9582", "0.8443", false],
    ["HCC1395 5kHz paired",   "LongPhase-S", "29754", "627",  "9693", "0.8522", false],
    ["HCC1395 5kHz paired",   "InterSubMod", "29752", "544",  "9695", "0.8532", true],
    ["HCC1395_DORADO paired", "ClairS",      "29986", "588",  "9461", "0.8565", false],
    ["HCC1395_DORADO paired", "LongPhase-S", "29889", "240",  "9558", "0.8592", false],
    ["HCC1395_DORADO paired", "InterSubMod", "29877", "238",  "9570", "0.8590", true],
  ];

  const rows = data.map((r) => {
    const ism = r[6];
    const fill = ism ? { color: "DDEEFF" } : undefined;
    return [
      td(r[0], { fill }),
      td(r[1], { bold: ism, fill }),
      tn(r[2], { fill }),
      tn(r[3], { fill }),
      tn(r[4], { fill }),
      tn(r[5], { bold: ism, color: ism ? C.navy : C.text, fill }),
    ];
  });

  slide.addTable(
    [[th("Dataset"), th("Method"), th("TP"), th("FP"), th("FN"), th("F1")], ...rows],
    { x: 0.3, y: 1.12, w: 9.4, colW: [2.55, 1.75, 1.05, 1.05, 1.05, 1.95], rowH: 0.44, border: { pt: 0.5, color: C.border } }
  );

  // Takeaways
  const ta = [
    "5kHz：InterSubMod F1=0.8532 > LongPhase-S 0.8522（+0.0010）。增益主要來自削 FP：627→544（-83），TP 僅少 2。在最 noisy 的場景有小幅但真實的 precision 改善。",
    "DORADO：InterSubMod 0.8590 ≈ LongPhase-S 0.8592，幾乎持平略負。DORADO 本身精確度已高，甲基整合邊際效益小。",
    "結論：5kHz 是主要 study dataset；同一套規則在 DORADO 可攜性受限，後續結論均需 cross-platform 驗證。",
  ];

  ta.forEach((t, i) => {
    slide.addText([
      { text: `${i + 1}. `, options: { bold: true, color: C.teal } },
      { text: t },
    ], {
      x: 0.3, y: 4.0 + i * 0.42, w: 9.4, h: 0.38,
      fontSize: 11, fontFace: F.zh, color: C.text, margin: 0, wrap: true,
    });
  });

  slide.addNotes(
    "【Paired Pure Benchmark】\n" +
    "5kHz 增益幅度不大（+0.001），但來源清楚：削 83 個 FP，而非大量救 TP。\n" +
    "DORADO 幾乎持平，可能因為 DORADO base quality 更好，甲基整合空間小。\n" +
    "這也說明為什麼後續主力轉向 TO 場景做 rescue 驗證。\n" +
    "數據來源：20260307_純樣本round執行與進度報告_01.md"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 8 – CANDIDATE-SPECIFIC RESCUE
// ─────────────────────────────────────────────────────────────────────────
function buildCandidateRescue(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.bg };
  addHeader(pres, slide, "Candidate-Specific Rescue：甲基 Support 四象限驗證", "Candidate-Specific Methylation Rescue Across 4 Quadrants");

  const data = [
    ["5kHz TO",       "qual>=10 or gq>=20",          "491", "118", "+0.006824", "caller-only 最佳",             C.green],
    ["5kHz TO",       "gq>=10 + Pairwise>=0.20",      "300", "68",  "+0.004219", "甲基 support 有效",             C.teal],
    ["5kHz TO",       "gq>=10 + agreement_positive",  "148", "25",  "+0.002163", "label/cluster 交叉訊號有價值",  C.teal],
    ["DORADO TO",     "gq>=10",                       "40",  "11",  "+0.000540", "caller-only 最佳",             C.green],
    ["DORADO TO",     "gq>=10 + Quality>=60",          "36",  "11",  "+0.000476", "有正向 support，未超過 caller", C.teal],
    ["DORADO TO",     "gq>=10 + Pairwise<=0.20",       "30",  "9",   "+0.000398", "方向與 5kHz TO 相反（⚠）",     C.amber],
    ["DORADO paired", "gq>=15",                       "97",  "88",  "+0.000502", "caller-first 成立",            C.green],
    ["DORADO paired", "gq>=10 + Pairwise>=0.20",      "10",  "36",  "-0.000280", "甲基 support 不成立（⚠）",     C.red],
    ["DORADO paired", "gq>=10 + agreement_positive",  "36",  "72",  "-0.000298", "甲基 support 不成立（⚠）",     C.red],
  ];

  const rows = data.map((r, i) => {
    const fill = rowFill(i);
    const f1IsNeg = r[4].startsWith("-");
    return [
      td(r[0], { bold: true, fill }),
      td(r[1], { color: C.muted, fill }),
      tn(r[2], { fill }),
      tn(r[3], { fill }),
      tn(r[4], { bold: true, color: f1IsNeg ? C.red : C.green, fill }),
      td(r[5], { color: r[6], fill }),
    ];
  });

  slide.addTable(
    [[th("Dataset"), th("Rule"), th("TP↑"), th("FP↑"), th("ΔF1"), th("判讀")], ...rows],
    { x: 0.3, y: 1.12, w: 9.4, colW: [1.5, 2.35, 0.68, 0.68, 0.88, 3.31], rowH: 0.35, border: { pt: 0.5, color: C.border } }
  );

  // Bottom note box
  slide.addShape(pres.ShapeType.rect, {
    x: 0.3, y: 4.64, w: 9.4, h: 0.68,
    fill: { color: C.bg_alt }, line: { color: C.border, width: 0.5 },
  });
  slide.addText(
    "注意：各行 delta F1 來自不同評估框架（複合規則 / GQ sweep / annotation policy），數值不可直接跨行比較。" +
    "  paired 結論保守（Coverage Ceiling 見下頁）。",
    {
      x: 0.4, y: 4.66, w: 9.2, h: 0.64,
      fontSize: 10, fontFace: F.zh, color: C.muted, margin: 0, wrap: true, valign: "middle",
    }
  );

  slide.addNotes(
    "【Candidate-Specific Rescue】\n" +
    "5kHz TO：甲基 support 成立。最佳 caller gate（qual>=10 or gq>=20）+491 TP/118 FP，F1+0.006824。\n" +
    "  gq>=10 + Pairwise>=0.20 有效（+0.004219）；agreement_positive 也有正向（+0.002163）。\n" +
    "DORADO TO：也有正向 support，但 Pairwise 最佳方向相反（<=0.20），代表 Pairwise 是 dataset-specific 訊號。\n" +
    "DORADO paired：兩條甲基規則都負效果。但需參照 coverage ceiling 解讀（下一頁）。\n" +
    "數據來源：20260308_HCC1395_5kHz_TO_candidate_specific甲基rescue驗證_01.md 等"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 9 – COVERAGE CEILING
// ─────────────────────────────────────────────────────────────────────────
function buildCoverage(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.bg };
  addHeader(pres, slide, "Coverage Ceiling：為何 Paired 結論要保守？", "Why Paired Methylation Results Should Be Interpreted Conservatively");

  // Big callout numbers — left
  slide.addShape(pres.ShapeType.rect, {
    x: 0.3, y: 1.12, w: 4.3, h: 1.8,
    fill: { color: "FFF8E1" }, line: { color: C.amber, width: 2 },
  });
  slide.addText("8 ~ 20%", {
    x: 0.35, y: 1.2, w: 4.2, h: 1.0,
    fontSize: 52, fontFace: F.en, bold: true,
    color: C.amber, align: "center", margin: 0,
  });
  slide.addText("Paired 的 candidate-specific 甲基 coverage 上限", {
    x: 0.35, y: 2.18, w: 4.2, h: 0.52,
    fontSize: 12, fontFace: F.zh, color: C.muted, align: "center", margin: 0, wrap: true,
  });

  // Big callout numbers — right
  slide.addShape(pres.ShapeType.rect, {
    x: 5.1, y: 1.12, w: 4.3, h: 1.8,
    fill: { color: "E0F2F1" }, line: { color: C.teal, width: 2 },
  });
  slide.addText("87 ~ 100%", {
    x: 5.15, y: 1.2, w: 4.2, h: 1.0,
    fontSize: 52, fontFace: F.en, bold: true,
    color: C.teal, align: "center", margin: 0,
  });
  slide.addText("TO 的 candidate-specific 甲基 coverage 上限", {
    x: 5.15, y: 2.18, w: 4.2, h: 0.52,
    fontSize: 12, fontFace: F.zh, color: C.muted, align: "center", margin: 0, wrap: true,
  });

  // Table
  const data = [
    ["HCC1395 5kHz paired",   "11.41%", "10.20%", "coverage ceiling 偏低，甲基無效結論需保守"],
    ["HCC1395_DORADO paired", "8.29%",  "19.66%", "coverage ceiling 更低，不能解讀成 paired 無甲基價值"],
    ["HCC1395 5kHz TO",       "87.32%", "71.48%", "TO coverage 足夠，support 結論較可信"],
    ["HCC1395_DORADO TO",     "100%",   "100%",   "TO coverage 完整，可直接判讀 support 關係"],
  ];

  const rows = data.map((r, i) => {
    const isTO = r[0].includes("TO");
    const fill = rowFill(i);
    const valColor = isTO ? C.green : C.amber;
    return [
      td(r[0], { bold: true, fill }),
      tn(r[1], { bold: true, color: valColor, fill }),
      tn(r[2], { bold: true, color: valColor, fill }),
      td(r[3], { color: C.muted, fill }),
    ];
  });

  slide.addTable(
    [[th("Dataset"), th("Lost TP coverage"), th("Removed FP coverage"), th("解讀")], ...rows],
    { x: 0.3, y: 3.08, w: 9.4, colW: [2.1, 1.7, 1.85, 3.75], rowH: 0.36, border: { pt: 0.5, color: C.border } }
  );

  slide.addNotes(
    "【Coverage Ceiling】\n" +
    "關鍵：paired 只有 8~20% 的 candidate-specific methylation coverage，TO 有 87~100%。\n" +
    "因此 paired 的「甲基無效」可能部分源自覆蓋率不足，而非甲基訊號本身無效。\n" +
    "TO 的 coverage 幾乎完整，因此 TO 的 support/non-support 結論更可信。\n" +
    "DORADO TO 用了 candidate-window subset BAM + subset haplotag 方法，避免生成 200G+ full BAM。"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 10 – PHASE 2 RAW BENCHMARK
// ─────────────────────────────────────────────────────────────────────────
function buildPhase2Raw(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.bg };
  addHeader(pres, slide, "Phase 2：Paired Raw Pileup/Full 直接 Benchmark", "Phase 2: Raw Caller Models vs Merged Output Benchmark");

  const data = [
    ["5kHz paired",   "raw pileup_filter",       "30595", "7450",  "8852", "0.789630", false, true],
    ["5kHz paired",   "raw full_alignment_filter","30759", "12684", "8688", "0.742164", false, true],
    ["5kHz paired",   "merged output",            "29865", "1430",  "9582", "0.844336", true,  false],
    ["DORADO paired", "raw pileup_filter",        "30872", "2232",  "8575", "0.851043", false, true],
    ["DORADO paired", "raw full_alignment_filter","31332", "2559",  "8115", "0.854455", false, true],
    ["DORADO paired", "merged output",            "29986", "588",   "9461", "0.856486", true,  false],
  ];

  const rows = data.map((r) => {
    const isMerged = r[6];
    const isRaw = r[7];
    const fill = isMerged ? { color: "DDEEFF" } : undefined;
    return [
      td(r[0], { fill }),
      td(r[1], { bold: isMerged, fill }),
      tn(r[2], { fill }),
      tn(r[3], { bold: isRaw, color: isRaw ? C.red : C.text, fill }),
      tn(r[4], { fill }),
      tn(r[5], { bold: isMerged, color: isMerged ? C.navy : C.text, fill }),
    ];
  });

  slide.addTable(
    [[th("Dataset"), th("Source"), th("TP"), th("FP"), th("FN"), th("F1")], ...rows],
    { x: 0.3, y: 1.12, w: 9.4, colW: [1.9, 2.5, 1.0, 1.0, 1.0, 2.0], rowH: 0.44, border: { pt: 0.5, color: C.border } }
  );

  const ta = [
    "Raw caller 確有額外 recall：5kHz pileup TP=30595 vs merged 29865（+730 TP），DORADO pileup +886 TP",
    "但 precision 代價極大：5kHz pileup FP=7450 vs merged 1430（5× 以上），F1 從 0.844 降至 0.790",
    "結論：paired 場景下，raw pileup/full 不適合當最終輸出；但可考慮用 raw 的額外 TP 作為 rescue candidate pool",
  ];

  ta.forEach((t, i) => {
    slide.addText([
      { text: `${i + 1}. `, options: { bold: true, color: C.teal } },
      { text: t },
    ], {
      x: 0.3, y: 4.0 + i * 0.42, w: 9.4, h: 0.38,
      fontSize: 11, fontFace: F.zh, color: C.text, margin: 0, wrap: true,
    });
  });

  slide.addNotes(
    "【Phase 2 Raw Benchmark】\n" +
    "核心問題：raw pileup/full 有多少額外 recall？precision 代價是多少？\n" +
    "5kHz：raw pileup +730 TP，但 FP 增加 6020（5× 倍增），F1 從 0.844 掉到 0.790。\n" +
    "DORADO：raw model FP 增加量相對較小（pileup +1644 FP），F1 從 0.856 降到 0.851，幅度小很多。\n" +
    "這反映 DORADO caller 本身更精確，raw 的噪音較少。\n" +
    "未來方向：用 raw 的額外 TP 作為 candidate pool 做 candidate-specific rescue，可能比直接拿 raw 輸出更合適。\n" +
    "數據來源：phase2_paired_model_feature_analysis/ 目錄"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 11 – ANNOTATION LAYER
// ─────────────────────────────────────────────────────────────────────────
function buildAnnotation(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.bg };
  addHeader(pres, slide, "Annotation Layer：最佳策略已定型", "Annotation Layer: Best Policy by Dataset & Its Value");

  const data = [
    ["HCC1395 5kHz TO",       "caller_primary_only (gq>=10)",  "+0.006754", "annotation 有價值，但沒超過 caller"],
    ["HCC1395 5kHz paired",   "caller_primary_only (gq>=15)",  "+0.000873", "paired 仍由 caller 主導"],
    ["HCC1395_DORADO paired", "caller_strict_only (gq>=20)",   "+0.000635", "更嚴格 caller gate 最佳"],
    ["HCC1395_DORADO TO",     "caller_primary_only (gq>=10)",  "+0.000540", "annotation 幾乎追平，但仍略低"],
  ];

  const rows = data.map((r, i) => {
    const fill = rowFill(i);
    return [
      td(r[0], { bold: true, fill }),
      td(r[1], { color: C.navy, fill }),
      tn(r[2], { bold: true, color: C.green, fill }),
      td(r[3], { color: C.muted, fill }),
    ];
  });

  slide.addTable(
    [[th("Dataset"), th("Best Policy"), th("ΔF1 vs baseline"), th("判讀")], ...rows],
    { x: 0.3, y: 1.12, w: 9.4, colW: [2.0, 2.7, 1.5, 3.2], rowH: 0.42, border: { pt: 0.5, color: C.border } }
  );

  // Value explanation
  slide.addText("Annotation Layer 的真正制度化價值", {
    x: 0.3, y: 3.02, w: 9, h: 0.32,
    fontSize: 13.5, fontFace: F.zh, bold: true, color: C.navy, margin: 0,
  });

  const pts = [
    "讓 Quality_Score [55,60)、Pairwise dataset-specific top-bin 等區間訊號可被穩定輸出，不需每次重算",
    "讓後續報告與人工審閱能直接判讀「此位點是否落在驗證過的 support 帶」",
    "讓 cross-dataset 解讀時能明確區分：support（甲基正向）/ QC（HP 完整度）/ artifact（低 VAF+高 AD+低 CV）三種標記",
  ];

  pts.forEach((p, i) => {
    slide.addText([
      { text: `${i + 1}. `, options: { bold: true, color: C.teal } },
      { text: p },
    ], {
      x: 0.35, y: 3.42 + i * 0.58, w: 9.2, h: 0.5,
      fontSize: 11.5, fontFace: F.zh, color: C.text, margin: 0, wrap: true,
    });
  });

  slide.addNotes(
    "【Annotation Layer】\n" +
    "四個象限的最佳 policy 都是 caller-primary 或 caller-strict，說明 annotation 的 best case 仍以 caller gate 為主。\n" +
    "但 annotation 的價值不在於取代 caller，而是把有效但不穩定的訊號制度化：\n" +
    "  1. 區間訊號穩定輸出\n" +
    "  2. 人工審閱效率提高\n" +
    "  3. 三種標記（support/QC/artifact）讓跨 dataset 解讀更清晰\n" +
    "數據來源：annotation_best_policy.tsv（研究輸出目錄）"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 12 – GAPS / LIMITATIONS
// ─────────────────────────────────────────────────────────────────────────
function buildGaps(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.bg };
  addHeader(pres, slide, "目前尚缺與待補齊事項", "Current Research Gaps & Limitations");

  const gaps = [
    {
      num: "1", severity: "高",
      title: "Paired Coverage Ceiling 尚未量化呈現",
      detail: "paired 的 candidate-specific coverage 僅 8~20%，「甲基無效」結論需保守解讀；需補量化數字到主報告正文而非只放在 caveat 段。",
    },
    {
      num: "2", severity: "中",
      title: "TO 雙模型（Pileup/Full）尚未對比",
      detail: "目前 TO 結論完全以 pileup 路線為主，full model 在 TO 場景的表現尚未測試，存在未知的 model 差異。",
    },
    {
      num: "3", severity: "中",
      title: "Annotation Layer 尚未接入整合輸出",
      detail: "policy evaluation 已完成，但尚未正式接到：dashboard / round summary / candidate diagnostics 統一輸出欄位。",
    },
    {
      num: "4", severity: "低",
      title: "Pool B FN 探勘結論尚未納入主報告",
      detail: "ClairS non-PASS FN 確有可研究空間，但甲基不適合單獨作為 Pool B rescue 主規則，結論尚未整合進本週主線敘事。",
    },
    {
      num: "5", severity: "低",
      title: "F1 delta 統計顯著性未驗證",
      detail: "多數增益在第 4 位小數（+0.001~+0.007），未做 bootstrap 95% CI 或多重比較校正，難以區分真實增益與隨機波動。",
    },
  ];

  const severityColor = { "高": C.red, "中": C.amber, "低": C.muted };

  gaps.forEach((g, i) => {
    const y = 1.15 + i * 0.73;
    // Number badge
    slide.addShape(pres.ShapeType.rect, {
      x: 0.3, y: y + 0.04, w: 0.44, h: 0.44,
      fill: { color: severityColor[g.severity] }, line: { color: severityColor[g.severity] },
    });
    slide.addText(g.num, {
      x: 0.3, y: y + 0.04, w: 0.44, h: 0.44,
      fontSize: 14, fontFace: F.en, bold: true,
      color: C.white, align: "center", valign: "middle", margin: 0,
    });
    // Severity label
    slide.addText(`[${g.severity}]`, {
      x: 0.82, y: y + 0.04, w: 0.5, h: 0.26,
      fontSize: 9, fontFace: F.zh, bold: true,
      color: severityColor[g.severity], margin: 0,
    });
    // Title
    slide.addText(g.title, {
      x: 0.82, y: y + 0.1, w: 8.85, h: 0.3,
      fontSize: 12.5, fontFace: F.zh, bold: true,
      color: C.navy, margin: 0,
    });
    // Detail
    slide.addText(g.detail, {
      x: 0.82, y: y + 0.4, w: 8.85, h: 0.28,
      fontSize: 10.5, fontFace: F.zh,
      color: C.muted, margin: 0, wrap: true,
    });
  });

  slide.addNotes(
    "【待補齊事項】\n" +
    "優先順序：\n" +
    "  高（coverage ceiling）— 影響結論的可信度邊界\n" +
    "  中（TO 雙模型 + annotation 接入）— 影響研究速度和可重複性\n" +
    "  低（Pool B FN + 統計顯著性）— 完整性和可重複性改進\n" +
    "最高優先動作：annotation layer 接入整合矩陣（對下一步直接可行）"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 13 – NEXT STEPS (dark)
// ─────────────────────────────────────────────────────────────────────────
function buildNextSteps(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.dark };

  slide.addText("下一步建議", {
    x: 0.5, y: 0.18, w: 9, h: 0.65,
    fontSize: 26, fontFace: F.zh, bold: true,
    color: C.white, margin: 0,
  });
  slide.addText("Recommended Next Steps — Priority Actions", {
    x: 0.5, y: 0.8, w: 9, h: 0.24,
    fontSize: 11, fontFace: F.en, color: "90A4AE", margin: 0,
  });

  const steps = [
    {
      p: "P1", color: C.teal,
      title: "Annotation Layer 接入整合矩陣",
      detail: "把本輪 annotation layer 結果正式接到 round summary / candidate diagnostics 統一輸出，不急著改 core filter。",
    },
    {
      p: "P1", color: C.teal,
      title: "TO Annotation Score/Ranking 實驗",
      detail: "在 TO 路線優先做 annotation ranking 實驗，驗證是否能改善人工審閱效率與候選排序品質。",
    },
    {
      p: "P2", color: C.amber,
      title: "Paired Coverage Ceiling 量化補充",
      detail: "若要繼續深挖 paired，先補 coverage ceiling 的量化呈現，再決定是否試新甲基規則。",
    },
    {
      p: "P3", color: C.indigo,
      title: "Artifact Triage 維持支線",
      detail: "artifact_lowvaf_highad_lowcv 維持 triage 支線，不與 rescue 主線混用。",
    },
  ];

  steps.forEach((s, i) => {
    const x = i < 2 ? 0.4 : 5.2;
    const y = i < 2 ? (1.15 + i * 1.7) : (1.15 + (i - 2) * 1.7);

    // Card
    slide.addShape(pres.ShapeType.rect, {
      x, y, w: 4.4, h: 1.55,
      fill: { color: "1A2A5E" },
      line: { color: s.color, width: 2 },
    });
    // Top bar
    slide.addShape(pres.ShapeType.rect, {
      x, y, w: 4.4, h: 0.1,
      fill: { color: s.color }, line: { color: s.color },
    });
    // Priority badge
    slide.addShape(pres.ShapeType.rect, {
      x: x + 0.14, y: y + 0.16, w: 0.52, h: 0.3,
      fill: { color: s.color }, line: { color: s.color },
    });
    slide.addText(s.p, {
      x: x + 0.14, y: y + 0.16, w: 0.52, h: 0.3,
      fontSize: 10, fontFace: F.en, bold: true,
      color: C.white, align: "center", valign: "middle", margin: 0,
    });
    // Title
    slide.addText(s.title, {
      x: x + 0.74, y: y + 0.16, w: 3.52, h: 0.34,
      fontSize: 13, fontFace: F.zh, bold: true,
      color: C.white, margin: 0,
    });
    // Detail
    slide.addText(s.detail, {
      x: x + 0.14, y: y + 0.54, w: 4.12, h: 0.88,
      fontSize: 11, fontFace: F.zh,
      color: "B0BEC5", margin: 0, wrap: true,
    });
  });

  slide.addNotes(
    "【下一步建議】\n" +
    "P1（最優先，兩項）：\n" +
    "  1. Annotation layer 接入 — 最直接可行的落地行動。\n" +
    "  2. TO annotation ranking — 驗證 annotation score 能否提升排序品質。\n" +
    "P2：先量化 paired coverage ceiling 問題規模，再決定是否繼續深挖。\n" +
    "P3：artifact triage 維持支線，不與主線混用。\n" +
    "核心原則：先做 annotation export，再做更大範圍的 cross-dataset ranking 驗證，不急著改 hard rule。"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// SLIDE 14 – APPENDIX
// ─────────────────────────────────────────────────────────────────────────
function buildAppendix(pres) {
  const slide = pres.addSlide();
  slide.background = { color: C.bg };
  addHeader(pres, slide, "附錄：來源文件與數據路徑", "Appendix: Source Documents & Data Output Paths");

  const BASE = "/big8_disk/liaoyoyo2001/InterSubMod/";
  const RUNS = "/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/";

  const sources = [
    { label: "主報告", path: BASE + "docs/reports/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md" },
    { label: "Phase 2 分析", path: BASE + "docs/experiments/in_progress/2026/03/20260311_phase2_paired_raw與feature_interval_orthogonality分析_01.md" },
    { label: "Annotation 驗證", path: BASE + "docs/experiments/in_progress/2026/03/20260311_phase2_finer_interval回接annotation_layer驗證_01.md" },
    { label: "四象限整合報告", path: BASE + "docs/reports/validated/2026/03/20260311_HCC1395_四象限甲基rescue整合報告_01.md" },
    { label: "Phase 2 輸出目錄", path: RUNS + "20260311_phase2_paired_model_feature_analysis/" },
    { label: "Annotation 輸出目錄", path: RUNS + "20260311_phase2_annotation_layer/" },
    { label: "Policy Evaluation TSV", path: RUNS + "20260311_phase2_annotation_layer/annotation_policy_evaluation.tsv" },
    { label: "Best Policy TSV", path: RUNS + "20260311_phase2_annotation_layer/annotation_best_policy.tsv" },
    { label: "Presence Summary TSV", path: RUNS + "20260311_phase2_annotation_layer/annotation_presence_summary.tsv" },
  ];

  sources.forEach((s, i) => {
    const y = 1.12 + i * 0.46;
    const bgFill = i % 2 === 0 ? C.bg : C.bg_alt;
    slide.addShape(pres.ShapeType.rect, {
      x: 0.28, y: y - 0.01, w: 9.44, h: 0.42,
      fill: { color: bgFill }, line: { color: C.border, width: 0.3 },
    });
    slide.addText(s.label + "：", {
      x: 0.35, y: y + 0.02, w: 1.55, h: 0.3,
      fontSize: 10, fontFace: F.zh, bold: true,
      color: C.navy, margin: 0,
    });
    // Path - truncate if too long for display
    const displayPath = s.path.length > 90 ? "..." + s.path.slice(-87) : s.path;
    slide.addText(displayPath, {
      x: 1.9, y: y + 0.03, w: 7.7, h: 0.28,
      fontSize: 8.5, fontFace: F.en,
      color: C.muted, margin: 0,
    });
  });

  slide.addNotes(
    "【附錄路徑】\n" +
    "所有路徑均為絕對路徑，可直接在 /big8_disk/liaoyoyo2001/ 環境下查閱。\n" +
    "主報告是本週 PPTX 的直接來源文件。\n" +
    "Phase 2 分析 + Annotation 驗證是本週兩個最新實驗文件。\n" +
    "四象限整合報告整合了 5kHz TO / DORADO TO / DORADO paired 三個象限的 rescue 驗證。"
  );
}

// ─────────────────────────────────────────────────────────────────────────
// MAIN
// ─────────────────────────────────────────────────────────────────────────
async function main() {
  const pres = new pptxgen();
  pres.layout = "LAYOUT_16x9";
  pres.author = "廖子游";
  pres.title = "InterSubMod 研究主線整合週報 2026-03-05 至 2026-03-11";

  buildCover(pres);
  buildAgenda(pres);
  buildOpeningThesis(pres);
  buildBackground(pres);
  buildGoals(pres);
  buildFeatureRoles(pres);
  buildPairedPure(pres);
  buildCandidateRescue(pres);
  buildCoverage(pres);
  buildPhase2Raw(pres);
  buildAnnotation(pres);
  buildGaps(pres);
  buildNextSteps(pres);
  buildAppendix(pres);

  await pres.writeFile({ fileName: OUTPUT_PATH });
  console.log("✓ PPTX written to:", OUTPUT_PATH);
}

main().catch((err) => {
  console.error("Error:", err);
  process.exit(1);
});
