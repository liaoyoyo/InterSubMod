#!/usr/bin/env node

import { createHash, randomUUID } from "node:crypto";
import { readFileSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { basename, join, resolve } from "node:path";
import { pathToFileURL } from "node:url";


function parseArguments(argv) {
  const values = {};
  for (let index = 0; index < argv.length; index += 2) {
    const key = argv[index];
    const value = argv[index + 1];
    if (!key?.startsWith("--") || !value) throw new Error(`Invalid argument pair: ${key ?? ""}`);
    values[key.slice(2)] = value;
  }
  for (const key of ["input", "output", "plugin-root", "receipt"]) {
    if (!values[key]) throw new Error(`--${key} is required`);
  }
  return values;
}


function sha256Text(text) {
  return createHash("sha256").update(text).digest("hex");
}


function injectOverflowCompatibilityStyle(html) {
  const style = [
    '<style id="intersubmod-portable-overflow-fix">',
    ".analytics-top-bar{width:100%!important;margin-left:0!important;margin-right:0!important}",
    ".chart-legend{width:100%!important;max-width:100%!important}",
    "</style>",
  ].join("");
  if (!html.includes("</head>")) throw new Error("Portable HTML has no closing head tag");
  return html.replace("</head>", `${style}</head>`);
}


async function main() {
  const args = parseArguments(process.argv.slice(2));
  const inputPath = resolve(args.input);
  const outputPath = resolve(args.output);
  const receiptPath = resolve(args.receipt);
  const scriptsRoot = resolve(args["plugin-root"], "skills", "build-report", "scripts");
  const [{ buildPortableArtifact }, { extractPortableChartSvgs }, { verifyPortableArtifact }] =
    await Promise.all([
      import(pathToFileURL(join(scriptsRoot, "build_portable_artifact.mjs"))),
      import(pathToFileURL(join(scriptsRoot, "extract_portable_chart_svgs.mjs"))),
      import(pathToFileURL(join(scriptsRoot, "verify_portable_artifact.mjs"))),
    ]);

  const artifactText = readFileSync(inputPath, "utf8");
  const artifact = JSON.parse(artifactText);
  const baseHtml = buildPortableArtifact(artifact);
  const temporaryBase = join(
    tmpdir(),
    `${basename(outputPath)}.${process.pid}.${randomUUID()}.base.html`,
  );
  writeFileSync(temporaryBase, baseHtml, "utf8");
  const staticCharts = await extractPortableChartSvgs({
    actionTimeoutMs: 5000,
    htmlPath: temporaryBase,
    readyTimeoutMs: 20000,
  });
  const finalHtml = injectOverflowCompatibilityStyle(
    buildPortableArtifact(artifact, { staticCharts }),
  );
  writeFileSync(outputPath, finalHtml, "utf8");
  const verification = await verifyPortableArtifact({
    actionTimeoutMs: 5000,
    artifactPath: inputPath,
    htmlPath: outputPath,
    readyTimeoutMs: 20000,
    screenshotPath: `${outputPath}.verification-failure.png`,
    timeoutMs: 60000,
  });
  const receipt = {
    schema_name: "intersubmod.portable_report_package_receipt",
    schema_version: "1.0.0",
    status: "PASS",
    input: inputPath,
    input_sha256: sha256Text(artifactText),
    output: outputPath,
    output_sha256: sha256Text(finalHtml),
    static_chart_block_ids: Object.keys(staticCharts).sort(),
    compatibility_override: {
      reason: "Portable reader 100vw topbar and unconstrained Recharts legend overflow with non-overlay scrollbars",
      scope: [".analytics-top-bar", ".chart-legend"],
      data_or_chart_spec_changed: false,
    },
    verification,
  };
  writeFileSync(receiptPath, `${JSON.stringify(receipt, null, 2)}\n`, "utf8");
  process.stdout.write(`${JSON.stringify(receipt)}\n`);
}


await main();
