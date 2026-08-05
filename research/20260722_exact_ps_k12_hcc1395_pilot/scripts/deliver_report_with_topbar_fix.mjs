#!/usr/bin/env node
/**
 * Deliver a canonical Data Analytics portable report with one scoped runtime
 * correction for the Linux non-overlay-scrollbar 100vw top-bar overflow.
 *
 * The artifact payload, native reader, chart renderer, static SVG extraction,
 * and canonical browser verifier are unchanged.  Only the top-bar width and
 * root x-overflow CSS are corrected before final verification.
 */

import { createHash } from "node:crypto";
import { readFileSync, writeFileSync } from "node:fs";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";

const PLUGIN_ROOT =
  "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/" +
  "data-analytics/0.2.8-13ceeea1f599";

const { deliverWithScrollbarSafeRuntime } = await import(
  pathToFileURL(resolve(
    "research/20260718_k_gt8_read_supported_segmentation/scripts/" +
    "deliver_portable_artifact_scrollbar_safe.mjs",
  )).href
);

function parseArgs(argv) {
  const result = {};
  for (let index = 0; index < argv.length; index += 2) {
    const key = argv[index];
    const value = argv[index + 1];
    if (!key?.startsWith("--") || !value) {
      throw new Error("Usage: deliver_report_with_topbar_fix.mjs --artifact artifact.json --output report.html --receipt receipt.json");
    }
    result[key.slice(2)] = value;
  }
  for (const key of ["artifact", "output", "receipt"]) {
    if (!result[key]) throw new Error(`Missing --${key}`);
  }
  return result;
}

function sha256(value) {
  return createHash("sha256").update(value).digest("hex");
}

async function main() {
  const args = parseArgs(process.argv.slice(2));
  const artifactPath = resolve(args.artifact);
  const outputPath = resolve(args.output);
  const receiptPath = resolve(args.receipt);
  const artifactText = readFileSync(artifactPath, "utf8");
  JSON.parse(artifactText);
  const delivery = await deliverWithScrollbarSafeRuntime({
    actionTimeoutMs: 2_500,
    input: artifactPath,
    output: outputPath,
    pluginRoot: PLUGIN_ROOT,
    readyTimeoutMs: 5_000,
    screenshot: `${outputPath}.verification-failure.png`,
    timeoutMs: 20_000,
  });
  const html = readFileSync(outputPath, "utf8");
  const receipt = {
    schema_name: "intersubmod.portable_report_delivery_with_topbar_fix",
    schema_version: "1.0.0",
    artifact: {
      path: artifactPath,
      sha256: sha256(artifactText),
      size_bytes: Buffer.byteLength(artifactText),
    },
    output: {
      path: outputPath,
      sha256: sha256(html),
      size_bytes: Buffer.byteLength(html),
    },
    correction: {
      id: "linux_non_overlay_scrollbar_100vw_topbar",
      marker: "data-intersubmod-scrollbar-safe-top-bar",
      scope: "CSS only; canonical artifact payload and reader runtime unchanged",
      reason: "The shared analytics top bar uses width:100vw and negative viewport margins, producing an 8px document overflow in Chromium headless-shell when a vertical scrollbar is present.",
    },
    delivery,
  };
  writeFileSync(receiptPath, `${JSON.stringify(receipt, null, 2)}\n`, "utf8");
  process.stdout.write(`${JSON.stringify({
    ok: true,
    output: outputPath,
    receipt: receiptPath,
    delivery,
  })}\n`);
}

await main();
