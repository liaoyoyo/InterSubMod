#!/usr/bin/env node

/**
 * Deliver a canonical Data Analytics portable report with one compatibility fix.
 *
 * The upstream reader makes its sticky header full-bleed with `width:100vw`.
 * Chromium with a classic vertical scrollbar includes the scrollbar gutter in
 * `vw`, which creates an 8 px document-level horizontal overflow at the
 * verifier's desktop viewport. The report content itself remains in bounds.
 *
 * This wrapper keeps the canonical renderer, chart extraction, structural
 * checks, source-dialog check, desktop check, and mobile check. It changes only
 * the sticky header from viewport width/negative viewport margins to the width
 * of its containing report surface.
 */

import { resolve } from "node:path";
import { pathToFileURL } from "node:url";

import { buildPortableArtifact } from
  "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/scripts/build_portable_artifact.mjs";
import { deliverPortableArtifact } from
  "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/scripts/deliver_portable_artifact.mjs";

const UPSTREAM_HEADER_CSS =
  "width:100vw;height:48px;min-height:48px;" +
  "margin-right:calc(50% - 50vw);margin-left:calc(50% - 50vw);";
const SCROLLBAR_SAFE_HEADER_CSS =
  "width:100%;height:48px;min-height:48px;margin-right:0;margin-left:0;";
const RUNTIME_SCROLLBAR_SAFE_STYLE = [
  '<style data-hcc1395-scrollbar-compat="true">',
  ".analytics-top-bar{width:100%!important;margin-right:0!important;margin-left:0!important}",
  "</style>",
].join("");

function buildScrollbarSafeArtifact(input, options = {}) {
  const html = buildPortableArtifact(input, options);
  const matches = html.split(UPSTREAM_HEADER_CSS).length - 1;
  if (matches !== 1) {
    throw new Error(
      `Expected one upstream 100vw header rule, found ${matches}; ` +
      "the portable reader CSS may have changed.",
    );
  }
  const patched = html.replace(UPSTREAM_HEADER_CSS, SCROLLBAR_SAFE_HEADER_CSS);
  if (!patched.includes("</head>")) {
    throw new Error("Portable artifact has no closing head tag.");
  }
  return patched.replace("</head>", `${RUNTIME_SCROLLBAR_SAFE_STYLE}</head>`);
}

function positiveNumber(value, fallback) {
  const parsed = Number(value);
  return Number.isFinite(parsed) && parsed > 0 ? parsed : fallback;
}

async function main(argv = process.argv.slice(2)) {
  const [inputPath, outputPath] = argv;
  if (!inputPath || !outputPath) {
    throw new Error(
      "Usage: deliver_portable_artifact_scrollbar_safe.mjs <artifact.json> <report.html>",
    );
  }
  const result = await deliverPortableArtifact(
    {
      actionTimeoutMs: positiveNumber(process.env.ACTION_TIMEOUT_MS, 10_000),
      inputPath,
      outputPath,
      readyTimeoutMs: positiveNumber(process.env.READY_TIMEOUT_MS, 30_000),
      screenshotPath: `${outputPath}.verification-failure.png`,
      timeoutMs: positiveNumber(process.env.VERIFICATION_TIMEOUT_MS, 60_000),
    },
    { build: buildScrollbarSafeArtifact },
  );
  process.stdout.write(`${JSON.stringify(result)}\n`);
}

const isMain =
  process.argv[1] && pathToFileURL(resolve(process.argv[1])).href === import.meta.url;
if (isMain) {
  try {
    await main();
  } catch (error) {
    const result = error?.verificationResult ?? {
      ok: false,
      code: error?.code ?? "delivery_failed",
      error: error?.message ?? String(error),
    };
    process.stderr.write(`${JSON.stringify(result)}\n`);
    process.exitCode = 1;
  }
}
