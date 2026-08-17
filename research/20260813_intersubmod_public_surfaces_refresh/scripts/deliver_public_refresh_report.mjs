#!/usr/bin/env node
/**
 * Deliver a canonical Data Analytics report while correcting the packaged
 * reader's 100vw/classic-scrollbar interaction. This does not replace the
 * official renderer: validation, HTML generation, chart-SVG extraction, and
 * final browser verification all come from the selected plugin package.
 */

import { copyFileSync, mkdtempSync, readFileSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { dirname, join, resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { mkdirSync } from "node:fs";


function usage() {
  return [
    "Usage: node deliver_public_refresh_report.mjs --plugin-root <dir> --input <artifact.json> --output <report.html>",
    "",
    "The adapter preserves the canonical artifact and official renderer, then",
    "adds a scoped width rule before running the official browser verifier.",
  ].join("\n");
}


function parseArguments(argv) {
  const values = {};
  for (let index = 0; index < argv.length; index += 1) {
    const key = argv[index];
    if (key === "--help" || key === "-h") return { help: true };
    if (!["--plugin-root", "--input", "--output"].includes(key)) {
      throw new Error(`Unknown argument: ${key}\n${usage()}`);
    }
    const value = argv[index + 1];
    if (!value || value.startsWith("--")) throw new Error(`Missing value for ${key}.`);
    values[key.slice(2)] = value;
    index += 1;
  }
  for (const key of ["plugin-root", "input", "output"]) {
    if (!values[key]) throw new Error(`--${key} is required.\n${usage()}`);
  }
  return values;
}


function addScrollbarWidthFix(html) {
  const marker = '<style id="intersubmod-portable-scrollbar-width-fix">';
  if (html.includes(marker)) return html;
  const style = `${marker}
html{scrollbar-gutter:stable}
@media screen and (min-width:601px){
  .analytics-top-bar{
    width:calc(100% + var(--ds-gutter) + var(--ds-gutter))!important;
    margin-right:calc(0px - var(--ds-gutter))!important;
    margin-left:calc(0px - var(--ds-gutter))!important;
  }
}
@media screen and (min-width:761px){
  .portable-page-header{
    width:calc(100% + 64px)!important;
    margin-right:-32px!important;
    margin-left:-32px!important;
  }
}
</style>`;
  if (!html.includes("</head>")) throw new Error("Packaged report has no </head> marker.");
  return html.replace("</head>", `${style}</head>`);
}


async function main() {
  const options = parseArguments(process.argv.slice(2));
  if (options.help) {
    process.stdout.write(`${usage()}\n`);
    return;
  }
  const pluginRoot = resolve(options["plugin-root"]);
  const inputPath = resolve(options.input);
  const outputPath = resolve(options.output);
  const scriptRoot = join(pluginRoot, "skills/build-report/scripts");
  const [{ buildPortableArtifact }, { extractPortableChartSvgs }, { verifyPortableArtifact }] =
    await Promise.all([
      import(pathToFileURL(join(scriptRoot, "build_portable_artifact.mjs"))),
      import(pathToFileURL(join(scriptRoot, "extract_portable_chart_svgs.mjs"))),
      import(pathToFileURL(join(scriptRoot, "verify_portable_artifact.mjs"))),
    ]);

  const artifact = JSON.parse(readFileSync(inputPath, "utf8"));
  const temporaryDirectory = mkdtempSync(join(tmpdir(), "intersubmod-public-report-"));
  const candidatePath = join(temporaryDirectory, "candidate.html");

  let html = buildPortableArtifact(artifact);
  writeFileSync(candidatePath, html, "utf8");
  const staticCharts = await extractPortableChartSvgs({ htmlPath: candidatePath });
  html = addScrollbarWidthFix(buildPortableArtifact(artifact, { staticCharts }));
  writeFileSync(candidatePath, html, "utf8");

  const verification = await verifyPortableArtifact({
    actionTimeoutMs: 3_000,
    artifactPath: inputPath,
    htmlPath: candidatePath,
    readyTimeoutMs: 7_500,
    timeoutMs: 20_000,
  });
  mkdirSync(dirname(outputPath), { recursive: true });
  copyFileSync(candidatePath, outputPath);
  process.stdout.write(`${JSON.stringify({
    ok: true,
    html: outputPath,
    compatibilityFix: "100vw-to-content-width for classic-scrollbar iframes",
    officialStages: {
      validationAndPackage: "passed",
      staticChartExtraction: "passed",
      verification: "passed",
    },
    verification,
  })}\n`);
}


main().catch((error) => {
  const verification = error?.verificationResult;
  process.stderr.write(`${JSON.stringify({
    ok: false,
    code: verification?.code ?? error?.code ?? "delivery_failed",
    error: verification?.error ?? error?.message ?? String(error),
    screenshot: verification?.screenshot,
  })}\n`);
  process.exitCode = 1;
});
