#!/usr/bin/env node

/**
 * Deliver the multi-BAM dashboard with a narrow shared-renderer overflow fix.
 *
 * The portable reader's full-bleed top bar uses 100vw. In Chromium frames with
 * a classic vertical scrollbar, 100vw includes the scrollbar gutter and makes
 * the document 8-15 px wider than documentElement.clientWidth. Re-expressing
 * the same full-bleed geometry from the dashboard content width preserves the
 * visual layout without clipping charts or table scroll regions.
 */

import { buildPortableArtifact } from "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/scripts/build_portable_artifact.mjs";
import { deliverPortableArtifact } from "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/scripts/deliver_portable_artifact.mjs";

const OVERFLOW_FIX = `
<style data-intersubmod-portable-overflow-fix="true">
.analytics-top-bar {
  width: calc(100% + var(--ds-gutter) + var(--ds-gutter)) !important;
  margin-right: calc(0px - var(--ds-gutter)) !important;
  margin-left: calc(0px - var(--ds-gutter)) !important;
}
@media (max-width: 600px) {
  .chart-legend {
    width: 100% !important;
    max-width: 100% !important;
    min-width: 0 !important;
    flex: 1 1 100% !important;
  }
}
</style>`;

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
  return html.replace("</head>", `${OVERFLOW_FIX}\n</head>`);
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
