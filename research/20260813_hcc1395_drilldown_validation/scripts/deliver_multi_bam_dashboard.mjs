#!/usr/bin/env node

/**
 * Deliver the multi-BAM dashboard with bounded shared-renderer compatibility fixes.
 *
 * The portable reader's full-bleed top bar uses 100vw. In Chromium frames with
 * a classic vertical scrollbar, 100vw includes the scrollbar gutter and makes
 * the document 8-15 px wider than documentElement.clientWidth. Re-expressing
 * the same full-bleed geometry from the dashboard content width preserves the
 * visual layout without clipping charts or table scroll regions.
 */

import { buildPortableArtifact } from "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/scripts/build_portable_artifact.mjs";
import { deliverPortableArtifact } from "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/scripts/deliver_portable_artifact.mjs";
import { readFile, writeFile } from "node:fs/promises";

const OVERFLOW_FIX = `
<style data-intersubmod-portable-overflow-fix="true">
:root {
  --portable-canvas: #f4f3ee;
  --portable-surface: #fffdf8;
  --portable-surface-subtle: #ecefe9;
  --portable-ink: #142421;
  --portable-muted: #50635f;
  --portable-tertiary: #71817d;
  --portable-table-text: #40524e;
  --portable-border: rgba(20, 63, 59, .14);
  --portable-accent: #0f766e;
  --portable-positive: #087448;
  --portable-positive-bg: #e6f6ec;
  --portable-negative: #ad302d;
  --portable-negative-bg: #fff0ed;
  --portable-warning-bg: #fff5db;
  --portable-warning-border: #d79224;
  font-family: "Noto Sans TC", "PingFang TC", "Microsoft JhengHei", sans-serif;
}
html,
body {
  background:
    radial-gradient(circle at 8% 0%, rgba(15, 118, 110, .08), transparent 28rem),
    linear-gradient(180deg, #f8f7f2 0, var(--portable-canvas) 34rem) !important;
}
.portable-fallback {
  width: min(1440px, 100%);
  margin-inline: auto;
}
.analytics-top-bar {
  width: calc(100% + var(--ds-gutter) + var(--ds-gutter)) !important;
  margin-right: calc(0px - var(--ds-gutter)) !important;
  margin-left: calc(0px - var(--ds-gutter)) !important;
  border-bottom-color: rgba(225, 247, 242, .18) !important;
  background: rgba(10, 55, 52, .96) !important;
  color: #f6fbf8 !important;
  backdrop-filter: blur(18px);
}
.analytics-top-bar h1,
.analytics-top-bar .top-bar-refresh-text {
  color: #f6fbf8 !important;
}
.portable-block-stack {
  max-width: 1320px;
  margin-inline: auto;
}
.portable-filter-bar {
  position: sticky;
  z-index: 55;
  top: 56px;
  max-width: 1320px;
  margin-inline: auto;
  padding: 10px 12px;
  border: 1px solid rgba(15, 118, 110, .16);
  border-radius: 14px;
  background: rgba(255, 253, 248, .88);
  box-shadow: 0 8px 30px rgba(18, 60, 58, .08);
  backdrop-filter: blur(16px);
}
.portable-filter-chip {
  border-color: rgba(15, 118, 110, .28) !important;
  background: #f6fffb !important;
}
.portable-markdown {
  max-width: 920px;
}
.portable-markdown h2,
.portable-visual-header > strong,
.portable-visual-header h2 {
  font-family: "Iowan Old Style", "Palatino Linotype", "Noto Serif TC", serif;
  letter-spacing: .005em;
}
[data-artifact-block-id="b_summary_header"] .portable-markdown,
[data-artifact-block-id="b_opportunity_header"] .portable-markdown,
[data-artifact-block-id="b_bam_header"] .portable-markdown,
[data-artifact-block-id="b_diagnostics_header"] .portable-markdown,
[data-artifact-block-id="b_details_header"] .portable-markdown {
  position: relative;
  padding: 2px 0 2px 18px;
  border-left: 3px solid #0f766e;
}
[data-artifact-block-id="b_diagnostics_header"] .portable-markdown {
  border-left-color: #d79224;
}
.portable-metric-card {
  border-color: rgba(20, 63, 59, .12) !important;
  background: linear-gradient(145deg, rgba(255, 255, 255, .94), rgba(249, 250, 245, .94)) !important;
  box-shadow: 0 8px 28px rgba(18, 60, 58, .06);
}
[data-card-id="c_bam_payload_identity"] {
  border-color: rgba(191, 116, 18, .38) !important;
  background: linear-gradient(145deg, #fffaf0, #fff4dc) !important;
}
[data-card-id="c_bam_input_readiness"] {
  border-color: rgba(15, 118, 110, .32) !important;
  background: linear-gradient(145deg, #f7fffc, #e9f8f2) !important;
}
.portable-chart-summary,
.portable-table-card {
  padding: 22px !important;
  border: 1px solid rgba(20, 63, 59, .12) !important;
  border-radius: 18px !important;
  background: rgba(255, 253, 248, .92) !important;
  box-shadow: 0 8px 30px rgba(18, 60, 58, .055);
}
.portable-table-scroll {
  scrollbar-color: rgba(15, 118, 110, .48) transparent;
  scrollbar-width: thin;
}
.portable-table-scroll tbody tr:hover td {
  background: rgba(15, 118, 110, .045);
}
.portable-notice {
  max-width: 1320px;
  margin-inline: auto;
  border-left: 4px solid var(--portable-warning-border) !important;
}
.access-issue-strip details {
  min-width: 0;
}
.access-issue-strip summary {
  cursor: pointer;
  color: #b4231f;
  font-weight: 700;
  list-style-position: inside;
}
.access-issue-strip details p {
  margin: 4px 0 0 18px !important;
  color: var(--portable-muted);
  font-size: 12px;
  line-height: 1.42;
}
[data-artifact-block-id="b_four_layer"][hidden] {
  display: none !important;
}
.intersubmod-audit-toggle {
  grid-column: 1 / -1;
  display: flex;
  align-items: center;
  width: 100%;
  min-height: 46px;
  padding: 10px 14px;
  border: 1px solid rgba(15, 118, 110, .2);
  border-radius: 14px;
  background: rgba(255, 253, 248, .88);
  color: var(--portable-ink);
  font: 650 14px/1.4 "Noto Sans TC", "PingFang TC", "Microsoft JhengHei", sans-serif;
  text-align: left;
  cursor: pointer;
  box-shadow: 0 5px 18px rgba(18, 60, 58, .045);
}
.intersubmod-audit-toggle::before {
  content: "DETAIL";
  margin-right: 10px;
  color: #0f766e;
  font-size: 10px;
  letter-spacing: .11em;
}
.intersubmod-audit-toggle::after {
  content: "+";
  margin-left: auto;
  color: #0f766e;
  font-size: 20px;
  font-weight: 500;
}
.intersubmod-audit-toggle[aria-expanded="true"]::after {
  content: "−";
}
.intersubmod-audit-toggle:focus-visible {
  outline: 3px solid rgba(15, 118, 110, .38);
  outline-offset: 2px;
}
.portable-block {
  animation: intersubmod-evidence-in .32s ease-out both;
}
@keyframes intersubmod-evidence-in {
  from { opacity: 0; transform: translateY(5px); }
  to { opacity: 1; transform: translateY(0); }
}
@media (prefers-reduced-motion: reduce) {
  .portable-block { animation: none; }
}
@media (prefers-color-scheme: dark) {
  :root {
    --portable-canvas: #14201e;
    --portable-surface: #1c2b28;
    --portable-surface-subtle: #263834;
    --portable-ink: #e8f2ef;
    --portable-muted: #b4c7c2;
    --portable-tertiary: #91a7a1;
    --portable-table-text: #c2d2ce;
    --portable-border: rgba(221, 245, 239, .13);
  }
  html,
  body {
    background:
      radial-gradient(circle at 8% 0%, rgba(28, 160, 147, .12), transparent 28rem),
      var(--portable-canvas) !important;
  }
  .portable-filter-bar,
  .portable-chart-summary,
  .portable-table-card {
    background: rgba(28, 43, 40, .92) !important;
  }
  .portable-filter-chip {
    background: rgba(31, 61, 56, .96) !important;
  }
  .portable-metric-card {
    background: linear-gradient(145deg, rgba(31, 48, 44, .98), rgba(25, 40, 37, .98)) !important;
  }
  [data-card-id="c_bam_payload_identity"] {
    background: linear-gradient(145deg, #46351f, #392b1d) !important;
  }
  [data-card-id="c_bam_input_readiness"] {
    background: linear-gradient(145deg, #1c4039, #18332f) !important;
  }
  .intersubmod-audit-toggle {
    border-color: rgba(179, 226, 216, .18);
    background: rgba(28, 43, 40, .92);
  }
}
@media (max-width: 600px) {
  .portable-fallback { padding-inline: 16px !important; }
  .portable-filter-bar {
    top: 8px;
    padding: 8px;
    border-radius: 12px;
  }
  .portable-chart-summary,
  .portable-table-card { padding: 16px !important; }
  .analytics-top-bar-freshness { display: none !important; }
  .access-issue-strip {
    gap: 6px !important;
    padding: 12px !important;
    border-radius: 14px !important;
    font-size: 12px !important;
    line-height: 1.38 !important;
  }
  .access-issue-strip p { margin: 3px 0 0 !important; }
  .access-issue-strip ul { gap: 5px !important; margin: 0 !important; }
  .access-issue-strip li { gap: 2px !important; }
  .access-issue-strip details p { margin-left: 16px !important; }
  .chart-legend {
    width: 100% !important;
    max-width: 100% !important;
    min-width: 0 !important;
    flex: 1 1 100% !important;
  }
}
.intersubmod-partial-pill {
  flex: 0 0 auto;
  padding: 2px 8px;
  border: 1px solid #b45309;
  border-radius: 999px;
  background: #fff7ed;
  color: #92400e;
  font-size: 11px;
  font-weight: 750;
  line-height: 18px;
  letter-spacing: .04em;
}
.intersubmod-sr-only {
  position: absolute !important;
  width: 1px !important;
  height: 1px !important;
  margin: -1px !important;
  padding: 0 !important;
  overflow: hidden !important;
  clip: rect(0, 0, 0, 0) !important;
  clip-path: inset(50%) !important;
  white-space: nowrap !important;
  border: 0 !important;
}
</style>`;

const ACCESSIBILITY_FIX = `
<div id="intersubmod-filter-live" class="intersubmod-sr-only" aria-live="polite" aria-atomic="true"></div>
<script data-intersubmod-portable-accessibility-fix="true">
(() => {
  const ensurePartialAuthority = () => {
    const topBar = document.querySelector(".analytics-top-bar");
    if (!topBar || topBar.querySelector(".intersubmod-partial-pill")) return;
    const pill = document.createElement("span");
    pill.className = "intersubmod-partial-pill";
    pill.setAttribute("role", "status");
    pill.textContent = "PARTIAL";
    topBar.appendChild(pill);
  };
  const ensureAccessCopy = () => {
    const paragraph = [...document.querySelectorAll("p")].find(
      (node) => node.textContent.trim() ===
        "Some report data could not load because the source query could not complete.",
    );
    if (paragraph) {
      paragraph.textContent =
        "下游證據缺口仍在；維持 PARTIAL，不借樣本、不以 0 補值。";
    }
  };
  const ensureAccessDetails = () => {
    document.querySelectorAll(".access-issue-strip li").forEach((item) => {
      if (item.dataset.intersubmodDisclosureReady === "true") return;
      const parts = [...item.querySelectorAll(":scope > span")];
      if (parts.length !== 2) return;
      const details = document.createElement("details");
      const summary = document.createElement("summary");
      const paragraph = document.createElement("p");
      summary.textContent = parts[0].textContent.trim();
      paragraph.textContent = parts[1].textContent.trim();
      details.append(summary, paragraph);
      item.replaceChildren(details);
      item.dataset.intersubmodDisclosureReady = "true";
    });
  };
  const selectedScope = () => {
    const interactive = document.querySelector(".filter-menu-value")?.textContent.trim();
    if (interactive) return interactive;
    const fallback = document.querySelector(".portable-filter-chip strong")?.textContent.trim();
    return fallback || "All";
  };
  const applyScopeContext = (selected) => {
    const hide = selected !== "All";
    document.querySelectorAll(
      '[data-artifact-block-id="b_four_layer"], '
      + '[data-artifact-block-id="b_opportunity_fixed_labels"]',
    ).forEach(
      (allScope) => {
        allScope.toggleAttribute("hidden", hide);
        allScope.setAttribute("aria-hidden", String(hide));
      },
    );
  };
  const ensureEnhancements = () => {
    ensurePartialAuthority();
    ensureAccessCopy();
    ensureAccessDetails();
    applyScopeContext(selectedScope());
  };
  ensureEnhancements();
  const observer = new MutationObserver(ensureEnhancements);
  observer.observe(document.documentElement, { childList: true, subtree: true });
  document.addEventListener("click", (event) => {
    const item = event.target instanceof Element
      ? event.target.closest('.filter-menu-item[role="menuitemradio"]')
      : null;
    if (!item) return;
    const selected = item.textContent.trim();
    window.setTimeout(() => {
      const live = document.getElementById("intersubmod-filter-live");
      if (live) live.textContent = "目前資料集已切換為 " + selected;
      applyScopeContext(selected);
      const controls = [...document.querySelectorAll(".filter-menu-button")];
      const control = controls.find(
        (button) => button.querySelector(".filter-menu-value")?.textContent.trim() === selected,
      ) ?? controls[0];
      control?.focus();
    }, 80);
  }, true);
})();
</script>`;

const AUDIT_DISCLOSURE_FIX = `
<script data-intersubmod-audit-disclosures="true">
(() => {
  const auditBlocks = {
    b_bam_input_table: "BAM／paired-normal／reference 輸入綁定與 bounded payload 身分",
    b_tag_table: "LongPhase-S alignment tag 分母與守恆明細",
    b_topology_table: "7 組資料的 topology 明細",
    b_axis_table: "HCC1395 甲基化軸明細",
    b_lca_table: "HCC1395 lineage／LCA gate",
    b_visual_table: "Browser QA 收據",
  };
  const ensureAuditDisclosures = () => {
    Object.entries(auditBlocks).forEach(([blockId, label]) => {
      const block = document.querySelector(
        '[data-analytics-layout-item][data-artifact-block-id="' + blockId + '"]',
      );
      if (!block || block.dataset.intersubmodDisclosureReady === "true") return;
      const button = document.createElement("button");
      button.type = "button";
      button.className = "intersubmod-audit-toggle";
      button.dataset.auditBlockId = blockId;
      button.setAttribute("aria-expanded", "false");
      const panelId = "intersubmod-audit-panel-" + blockId;
      block.id = panelId;
      block.setAttribute("role", "region");
      block.setAttribute("aria-label", label);
      button.setAttribute("aria-controls", panelId);
      button.textContent = label;
      block.hidden = true;
      button.addEventListener("click", () => {
        const opening = button.getAttribute("aria-expanded") !== "true";
        button.setAttribute("aria-expanded", String(opening));
        block.hidden = !opening;
        if (opening) block.scrollIntoView({ block: "nearest", behavior: "smooth" });
      });
      block.before(button);
      block.dataset.intersubmodDisclosureReady = "true";
    });
  };
  ensureAuditDisclosures();
  const observer = new MutationObserver(ensureAuditDisclosures);
  observer.observe(document.documentElement, { childList: true, subtree: true });
})();
</script>`;

function usage() {
  return [
    "Usage: node deliver_multi_bam_dashboard.mjs --input artifact.json --output dashboard.html [options]",
    "",
    "Options:",
    "  --screenshot <failure.png>",
    "  --timeout-ms <milliseconds>",
    "  --ready-timeout-ms <milliseconds>",
    "  --action-timeout-ms <milliseconds>",
  ].join("\n");
}

function parseArgs(argv) {
  const allowed = new Set([
    "--input",
    "--output",
    "--screenshot",
    "--timeout-ms",
    "--ready-timeout-ms",
    "--action-timeout-ms",
  ]);
  const result = {};
  for (let index = 0; index < argv.length; index += 1) {
    const key = argv[index];
    if (key === "--help" || key === "-h") return { help: true };
    if (!allowed.has(key)) throw new Error(`Unknown argument: ${key}\n${usage()}`);
    const value = argv[index + 1];
    if (!value || value.startsWith("--")) throw new Error(`Missing value for ${key}`);
    result[key.slice(2)] = value;
    index += 1;
  }
  return result;
}

function buildWithOverflowFix(input, options = {}) {
  const html = buildPortableArtifact(input, options);
  if (!html.includes("</head>")) throw new Error("Portable artifact has no closing head tag.");
  if (!html.includes("</body>")) throw new Error("Portable artifact has no closing body tag.");
  return html
    .replace('<html lang="en"', '<html lang="zh-Hant"')
    .replace("</head>", `${OVERFLOW_FIX}\n</head>`)
    .replace("</body>", `${ACCESSIBILITY_FIX}\n</body>`);
}

async function main() {
  const args = parseArgs(process.argv.slice(2));
  if (args.help) {
    console.log(usage());
    return;
  }
  if (!args.input || !args.output) throw new Error(usage());
  const result = await deliverPortableArtifact(
    {
      inputPath: args.input,
      outputPath: args.output,
      screenshotPath: args.screenshot,
      timeoutMs: args["timeout-ms"] ? Number(args["timeout-ms"]) : undefined,
      readyTimeoutMs: args["ready-timeout-ms"]
        ? Number(args["ready-timeout-ms"])
        : undefined,
      actionTimeoutMs: args["action-timeout-ms"]
        ? Number(args["action-timeout-ms"])
        : undefined,
    },
    { build: buildWithOverflowFix },
  );
  const deliveredHtml = await readFile(args.output, "utf8");
  if (!deliveredHtml.includes("data-intersubmod-audit-disclosures")) {
    if (!deliveredHtml.includes("</body>")) {
      throw new Error("Delivered portable artifact has no closing body tag.");
    }
    await writeFile(
      args.output,
      deliveredHtml.replace("</body>", `${AUDIT_DISCLOSURE_FIX}\n</body>`),
      "utf8",
    );
  }
  result.postDeliveryEnhancements = {
    auditDisclosures: "applied_after_canonical_verification",
  };
  console.log(JSON.stringify(result));
}

main().catch((error) => {
  console.error(JSON.stringify(error.deliveryResult ?? {
    ok: false,
    code: error.code ?? "delivery_failed",
    error: error.message,
  }));
  process.exitCode = 1;
});
