#!/usr/bin/env node

import { createHash } from "node:crypto";
import {
  chmodSync,
  existsSync,
  mkdirSync,
  readFileSync,
  renameSync,
  statSync,
  writeFileSync,
} from "node:fs";
import { dirname, join, resolve, sep } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";
import { performance } from "node:perf_hooks";

const TOPIC_ROOT = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const DEFAULT_PLUGIN_ROOT =
  "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/" +
  "data-analytics/0.2.8-13ceeea1f599";
const UPSTREAM_TOP_BAR_RULE =
  "width:100vw;height:48px;min-height:48px;" +
  "margin-right:calc(50% - 50vw);margin-left:calc(50% - 50vw);";
const PATCHED_TOP_BAR_RULE =
  "width:100%;height:48px;min-height:48px;margin-right:0;margin-left:0;";
const RUNTIME_PATCH =
  '<style data-intersubmod-classic-scrollbar-compat="true">' +
  ".analytics-top-bar{width:100%!important;margin-right:0!important;margin-left:0!important}" +
  "</style>";

function parseArguments(argv) {
  const values = new Map();
  for (let index = 0; index < argv.length; index += 2) {
    const key = argv[index];
    const value = argv[index + 1];
    if (!key?.startsWith("--") || value == null || value.startsWith("--")) {
      throw new Error(`Expected --key value pairs, observed ${JSON.stringify(argv.slice(index, index + 2))}`);
    }
    if (values.has(key)) throw new Error(`${key} may only be provided once`);
    values.set(key, value);
  }
  const required = ["--input", "--output", "--receipt"];
  for (const key of required) {
    if (!values.has(key)) throw new Error(`${key} is required`);
  }
  const positiveInteger = (key, fallback) => {
    const value = Number(values.get(key) ?? fallback);
    if (!Number.isInteger(value) || value < 1_000 || value > 60_000) {
      throw new Error(`${key} must be an integer between 1000 and 60000`);
    }
    return value;
  };
  return {
    inputPath: resolve(values.get("--input")),
    outputPath: resolve(values.get("--output")),
    receiptPath: resolve(values.get("--receipt")),
    screenshotPath: values.has("--failure-screenshot")
      ? resolve(values.get("--failure-screenshot"))
      : null,
    pluginRoot: resolve(values.get("--plugin-root") ?? DEFAULT_PLUGIN_ROOT),
    readyTimeoutMs: positiveInteger("--ready-timeout-ms", 15_000),
    actionTimeoutMs: positiveInteger("--action-timeout-ms", 5_000),
    timeoutMs: positiveInteger("--timeout-ms", 30_000),
  };
}

function assertTopicOutput(path, label) {
  if (!(path === TOPIC_ROOT || path.startsWith(`${TOPIC_ROOT}${sep}`))) {
    throw new Error(`${label} must be within ${TOPIC_ROOT}`);
  }
}

function assertAbsent(path, label) {
  if (existsSync(path)) throw new Error(`${label} already exists: ${path}`);
}

function sha256Path(path) {
  return createHash("sha256").update(readFileSync(path)).digest("hex");
}

function writeJsonExclusive(path, payload) {
  writeFileSync(path, `${JSON.stringify(payload, null, 2)}\n`, {
    encoding: "utf8",
    flag: "wx",
  });
}

function countOccurrences(text, token) {
  return text.split(token).length - 1;
}

function applyClassicScrollbarCompatibility(html) {
  const upstreamRuleCount = countOccurrences(html, UPSTREAM_TOP_BAR_RULE);
  if (upstreamRuleCount !== 1) {
    throw new Error(`Expected one canonical top-bar rule, observed ${upstreamRuleCount}`);
  }
  if (countOccurrences(html, "</head>") !== 1) {
    throw new Error("Expected exactly one closing head tag");
  }
  return {
    html: html
      .replace(UPSTREAM_TOP_BAR_RULE, PATCHED_TOP_BAR_RULE)
      .replace("</head>", `${RUNTIME_PATCH}</head>`),
    upstreamRuleCount,
    runtimePatchCount: 1,
  };
}

function compactError(error) {
  return {
    code: error?.code ?? error?.verificationResult?.code ?? "delivery_failed",
    message: String(error?.message ?? error),
    verificationResult: error?.verificationResult ?? null,
  };
}

async function main() {
  const options = parseArguments(process.argv.slice(2));
  const candidatePath = `${options.outputPath}.compat-candidate.html`;
  for (const [path, label] of [
    [options.outputPath, "output"],
    [options.receiptPath, "receipt"],
    [candidatePath, "candidate"],
  ]) {
    assertTopicOutput(path, label);
    assertAbsent(path, label);
  }
  if (options.screenshotPath) {
    assertTopicOutput(options.screenshotPath, "failure screenshot");
    assertAbsent(options.screenshotPath, "failure screenshot");
  }
  if (!existsSync(options.inputPath)) throw new Error(`artifact does not exist: ${options.inputPath}`);

  const modulePaths = {
    builder: join(options.pluginRoot, "skills/build-report/scripts/build_portable_artifact.mjs"),
    extractor: join(options.pluginRoot, "skills/build-report/scripts/extract_portable_chart_svgs.mjs"),
    verifier: join(options.pluginRoot, "skills/build-report/scripts/verify_portable_artifact.mjs"),
  };
  for (const [label, path] of Object.entries(modulePaths)) {
    if (!existsSync(path)) throw new Error(`${label} module does not exist: ${path}`);
  }

  const [{ buildPortableArtifact }, { extractPortableChartSvgs }, { verifyPortableArtifact }] =
    await Promise.all([
      import(pathToFileURL(modulePaths.builder)),
      import(pathToFileURL(modulePaths.extractor)),
      import(pathToFileURL(modulePaths.verifier)),
    ]);

  const startedAt = performance.now();
  let stage = "input";
  let patchReceipt = null;
  let staticChartCount = 0;
  let verification = null;
  const baseReceipt = {
    schema_name: "intersubmod.singleton_sidecar_portable_compat_delivery_receipt",
    schema_version: "1.0.0",
    created_at_utc: new Date().toISOString(),
    command: process.argv,
    artifact: {
      path: options.inputPath,
      sha256: sha256Path(options.inputPath),
    },
    plugin_modules: Object.fromEntries(
      Object.entries(modulePaths).map(([label, path]) => [
        label,
        { path, sha256: sha256Path(path) },
      ]),
    ),
    timeouts_ms: {
      ready: options.readyTimeoutMs,
      action: options.actionTimeoutMs,
      total: options.timeoutMs,
    },
    compatibility_scope:
      "Canonical plugin runtime only; analytics-top-bar classic-scrollbar width compatibility.",
  };

  try {
    const artifact = JSON.parse(readFileSync(options.inputPath, "utf8"));
    mkdirSync(dirname(options.outputPath), { recursive: true });

    stage = "canonical_build";
    let html = buildPortableArtifact(artifact);
    writeFileSync(candidatePath, html, { encoding: "utf8", flag: "wx" });

    stage = "canonical_static_chart_extraction";
    const staticCharts = await extractPortableChartSvgs({
      actionTimeoutMs: options.actionTimeoutMs,
      htmlPath: candidatePath,
      readyTimeoutMs: options.readyTimeoutMs,
    });
    staticChartCount = Object.keys(staticCharts).length;
    html = buildPortableArtifact(artifact, { staticCharts });

    stage = "classic_scrollbar_compatibility";
    const patched = applyClassicScrollbarCompatibility(html);
    html = patched.html;
    patchReceipt = {
      upstream_rule_occurrences: patched.upstreamRuleCount,
      runtime_patch_occurrences: patched.runtimePatchCount,
      forbidden_scientific_content_changes: 0,
    };
    writeFileSync(candidatePath, html, { encoding: "utf8", flag: "w" });

    stage = "canonical_verification";
    verification = await verifyPortableArtifact({
      actionTimeoutMs: options.actionTimeoutMs,
      artifactPath: options.inputPath,
      htmlPath: candidatePath,
      readyTimeoutMs: options.readyTimeoutMs,
      screenshotPath: options.screenshotPath,
      timeoutMs: options.timeoutMs,
    });

    stage = "publish";
    renameSync(candidatePath, options.outputPath);
    const receipt = {
      ...baseReceipt,
      ok: true,
      stages: {
        canonical_build: "passed",
        canonical_static_chart_extraction: "passed",
        classic_scrollbar_compatibility: "passed",
        canonical_verification: "passed",
        publish: "passed",
      },
      static_chart_count: staticChartCount,
      compatibility_patch: patchReceipt,
      verification,
      html: {
        path: options.outputPath,
        sha256: sha256Path(options.outputPath),
        size_bytes: statSync(options.outputPath).size,
      },
      elapsed_ms: Math.round((performance.now() - startedAt) * 10) / 10,
    };
    writeJsonExclusive(options.receiptPath, receipt);
    chmodSync(options.outputPath, 0o444);
    chmodSync(options.receiptPath, 0o444);
    console.log(JSON.stringify(receipt));
    return;
  } catch (error) {
    const receipt = {
      ...baseReceipt,
      ok: false,
      failed_stage: stage,
      static_chart_count: staticChartCount,
      compatibility_patch: patchReceipt,
      verification,
      error: compactError(error),
      candidate: existsSync(candidatePath)
        ? {
            path: candidatePath,
            sha256: sha256Path(candidatePath),
            size_bytes: statSync(candidatePath).size,
          }
        : null,
      elapsed_ms: Math.round((performance.now() - startedAt) * 10) / 10,
    };
    writeJsonExclusive(options.receiptPath, receipt);
    chmodSync(options.receiptPath, 0o444);
    if (existsSync(candidatePath)) chmodSync(candidatePath, 0o444);
    if (options.screenshotPath && existsSync(options.screenshotPath)) {
      chmodSync(options.screenshotPath, 0o444);
    }
    console.error(JSON.stringify(receipt));
    process.exitCode = 1;
  }
}

await main();
