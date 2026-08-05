#!/usr/bin/env node
/** Package with the canonical renderer plus one documented scrollbar compatibility fix. */

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


function injectTopBarScrollbarCompatibility(html) {
  const style = [
    '<style id="intersubmod-portable-topbar-scrollbar-compatibility">',
    ".analytics-top-bar{width:100%!important;margin-left:0!important;margin-right:0!important}",
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
  const temporaryBase = join(tmpdir(), `${basename(outputPath)}.${process.pid}.${randomUUID()}.base.html`);
  writeFileSync(temporaryBase, baseHtml, "utf8");
  const staticCharts = await extractPortableChartSvgs({
    actionTimeoutMs: 5_000,
    htmlPath: temporaryBase,
    readyTimeoutMs: 20_000,
  });
  const finalHtml = injectTopBarScrollbarCompatibility(
    buildPortableArtifact(artifact, { staticCharts }),
  );
  const temporaryFinal = join(tmpdir(), `${basename(outputPath)}.${process.pid}.${randomUUID()}.final.html`);
  writeFileSync(temporaryFinal, finalHtml, "utf8");
  const verification = await verifyPortableArtifact({
    actionTimeoutMs: 5_000,
    artifactPath: inputPath,
    htmlPath: temporaryFinal,
    readyTimeoutMs: 20_000,
    screenshotPath: `${outputPath}.verification-failure.png`,
    timeoutMs: 60_000,
  });

  // Publish only after the canonical verifier has passed both viewports.
  writeFileSync(outputPath, finalHtml, "utf8");
  const publishedVerification = { ...verification, html: outputPath };
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
      selector: ".analytics-top-bar",
      rule: "width:100%;margin-left:0;margin-right:0",
      reason: "Canonical reader uses width:100vw; non-overlay vertical scrollbar made clientWidth=1425 and scrollWidth=1433 in a 1440px verifier viewport.",
      content_hidden: false,
      data_or_chart_spec_changed: false,
    },
    pre_fix_failures: [
      {
        code: "horizontal_overflow",
        viewport_width: 1440,
        inner_width: 1440,
        client_width: 1425,
        scroll_width: 1433,
        offender: ".analytics-top-bar",
        bounding_box: { left: -7.5, right: 1432.5, width: 1440 },
      },
    ],
    verification: publishedVerification,
  };
  writeFileSync(receiptPath, `${JSON.stringify(receipt, null, 2)}\n`, "utf8");
  process.stdout.write(`${JSON.stringify(receipt)}\n`);
}


await main();
