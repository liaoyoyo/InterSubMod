#!/usr/bin/env node

/**
 * Package a report with the official Data Analytics runtime and apply one
 * audited non-overlay-scrollbar width correction before official verification.
 */

import { createHash, randomUUID } from "node:crypto";
import {
  chmodSync,
  existsSync,
  mkdtempSync,
  readFileSync,
  renameSync,
  statSync,
  writeFileSync,
} from "node:fs";
import { tmpdir } from "node:os";
import { basename, join, resolve } from "node:path";
import { pathToFileURL } from "node:url";

const PATCH_ID = "intersubmod-portable-scrollbar-width-fix-v1";
const PATCH_CSS = `<style id="${PATCH_ID}">
@media screen {
  .analytics-top-bar,
  .portable-page-header {
    width: 100% !important;
    margin-left: 0 !important;
    margin-right: 0 !important;
  }
}
</style>`;

function parseArguments(argv) {
  const options = {};
  const valueOptions = new Set([
    "action-timeout-ms",
    "input",
    "output",
    "plugin-root",
    "ready-timeout-ms",
    "receipt",
    "screenshot",
    "timeout-ms",
  ]);
  for (let index = 0; index < argv.length; index += 1) {
    const argument = argv[index];
    if (argument === "--help" || argument === "-h") return { help: true };
    if (!argument.startsWith("--")) throw new Error(`Unexpected argument: ${argument}`);
    const key = argument.slice(2);
    if (!valueOptions.has(key)) throw new Error(`Unknown option: ${argument}`);
    const value = argv[index + 1];
    if (!value || value.startsWith("--")) throw new Error(`Missing value for ${argument}`);
    if (options[key] !== undefined) throw new Error(`${argument} may only be supplied once.`);
    options[key] = value;
    index += 1;
  }
  for (const key of ["input", "output", "plugin-root", "receipt"]) {
    if (!options[key]) throw new Error(`--${key} is required.`);
  }
  return options;
}

function usage() {
  return [
    "Usage: deliver_portable_report_scrollbar_safe.mjs \\",
    "  --input artifact.json --output report.html --receipt receipt.json \\",
    "  --plugin-root /path/to/data-analytics-plugin [timeout options]",
  ].join("\n");
}

function positiveNumber(value, fallback, label) {
  if (value === undefined) return fallback;
  const parsed = Number(value);
  if (!Number.isFinite(parsed) || parsed <= 0) throw new Error(`${label} must be positive.`);
  return parsed;
}

function sha256Buffer(value) {
  return createHash("sha256").update(value).digest("hex");
}

function artifactIdentity(path) {
  const payload = readFileSync(path);
  return {
    path: resolve(path),
    size_bytes: payload.length,
    sha256: sha256Buffer(payload),
  };
}

function injectPatch(html) {
  if (html.includes(PATCH_ID)) throw new Error(`${PATCH_ID} is already present.`);
  const closingHead = html.indexOf("</head>");
  if (closingHead < 0) throw new Error("Official portable HTML has no closing head element.");
  return `${html.slice(0, closingHead)}${PATCH_CSS}\n${html.slice(closingHead)}`;
}

async function main() {
  const options = parseArguments(process.argv.slice(2));
  if (options.help) {
    process.stdout.write(`${usage()}\n`);
    return;
  }

  const inputPath = resolve(options.input);
  const outputPath = resolve(options.output);
  const receiptPath = resolve(options.receipt);
  const pluginRoot = resolve(options["plugin-root"]);
  for (const path of [outputPath, receiptPath]) {
    if (existsSync(path)) throw new Error(`Refusing to overwrite existing output: ${path}`);
  }

  const scriptRoot = join(pluginRoot, "skills", "build-report", "scripts");
  const builderPath = join(scriptRoot, "build_portable_artifact.mjs");
  const extractorPath = join(scriptRoot, "extract_portable_chart_svgs.mjs");
  const verifierPath = join(scriptRoot, "verify_portable_artifact.mjs");
  const [{ buildPortableArtifact }, { extractPortableChartSvgs }, { verifyPortableArtifact }] =
    await Promise.all([
      import(pathToFileURL(builderPath).href),
      import(pathToFileURL(extractorPath).href),
      import(pathToFileURL(verifierPath).href),
    ]);

  const artifact = JSON.parse(readFileSync(inputPath, "utf8"));
  const temporaryDirectory = mkdtempSync(join(tmpdir(), "intersubmod-portable-delivery-"));
  const initialPath = join(temporaryDirectory, `initial-${randomUUID()}.html`);
  const candidatePath = join(temporaryDirectory, `candidate-${randomUUID()}.html`);
  const initialHtml = buildPortableArtifact(artifact);
  writeFileSync(initialPath, initialHtml, "utf8");

  const readyTimeoutMs = positiveNumber(
    options["ready-timeout-ms"],
    45_000,
    "--ready-timeout-ms",
  );
  const actionTimeoutMs = positiveNumber(
    options["action-timeout-ms"],
    20_000,
    "--action-timeout-ms",
  );
  const timeoutMs = positiveNumber(options["timeout-ms"], 120_000, "--timeout-ms");
  const staticCharts = await extractPortableChartSvgs({
    actionTimeoutMs,
    htmlPath: initialPath,
    readyTimeoutMs,
  });
  const officialHtml = buildPortableArtifact(artifact, { staticCharts });
  const patchedHtml = injectPatch(officialHtml);
  writeFileSync(candidatePath, patchedHtml, "utf8");

  const verification = await verifyPortableArtifact({
    actionTimeoutMs,
    artifactPath: inputPath,
    htmlPath: candidatePath,
    readyTimeoutMs,
    screenshotPath: options.screenshot ? resolve(options.screenshot) : undefined,
    timeoutMs,
  });
  renameSync(candidatePath, outputPath);

  const receipt = {
    schema_name: "intersubmod.portable_report_delivery_receipt",
    schema_version: "1.0.0",
    status: "PASS",
    pass: true,
    artifact: artifactIdentity(inputPath),
    output: artifactIdentity(outputPath),
    official_runtime: {
      plugin_root: pluginRoot,
      builder: artifactIdentity(builderPath),
      chart_extractor: artifactIdentity(extractorPath),
      verifier: artifactIdentity(verifierPath),
    },
    static_chart_blocks: Object.keys(staticCharts).sort(),
    patch: {
      id: PATCH_ID,
      css_sha256: sha256Buffer(PATCH_CSS),
      scope: [".analytics-top-bar", ".portable-page-header"],
      change: "100vw full-bleed width to parent 100% with zero horizontal margins",
      reason: (
        "Official iframe verification uses non-overlay vertical scrollbars; 100vw includes " +
        "the scrollbar and causes an 8 px page-level overflow at the 1440 px desktop probe."
      ),
    },
    verification,
    retained_initial_html: {
      path: initialPath,
      basename: basename(initialPath),
      size_bytes: statSync(initialPath).size,
      sha256: artifactIdentity(initialPath).sha256,
      lifecycle: "disposable /tmp diagnostic; not a report deliverable",
    },
  };
  writeFileSync(receiptPath, `${JSON.stringify(receipt, null, 2)}\n`, "utf8");
  for (const path of [outputPath, receiptPath, options.screenshot ? resolve(options.screenshot) : undefined]) {
    if (path && existsSync(path)) chmodSync(path, 0o444);
  }
  process.stdout.write(`${JSON.stringify({ ok: true, output: outputPath, receipt: receiptPath, verification })}\n`);
}

main().catch((error) => {
  process.stderr.write(`${JSON.stringify({ ok: false, error: error?.message ?? String(error) })}\n`);
  process.exitCode = 1;
});
