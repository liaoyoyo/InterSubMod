#!/usr/bin/env node

import { createHash, randomUUID } from "node:crypto";
import {
  existsSync,
  mkdirSync,
  readFileSync,
  readdirSync,
  renameSync,
  statSync,
  writeFileSync,
} from "node:fs";
import { homedir } from "node:os";
import { basename, dirname, join, resolve } from "node:path";
import { performance } from "node:perf_hooks";
import { fileURLToPath, pathToFileURL } from "node:url";

const DEFAULT_PLUGIN_ROOT = join(
  homedir(),
  ".codex",
  "plugins",
  "cache",
  "openai-curated-remote",
  "data-analytics",
  "0.2.8-13ceeea1f599",
);
const DEFAULT_BROWSER_EXECUTABLE = join(
  homedir(),
  ".cache",
  "ms-playwright",
  "chromium_headless_shell-1223",
  "chrome-headless-shell-linux64",
  "chrome-headless-shell",
);

const PATCH_MARKER = "data-inter-sub-mod-portable-topbar-fix";
const PATCH_CSS = [
  ".analytics-top-bar{width:100%!important;margin-right:0!important;margin-left:0!important}",
  ".chart-legend{width:100%!important;max-width:100%!important;min-width:0!important}",
  ".chart-legend-item{min-width:0!important}",
  ".chart-legend-button{max-width:100%!important}",
].join("");
const EXPECTED_VIEWPORTS = Object.freeze([1440, 390]);
const EXPECTED_PLUGIN_NAME = "datascience-mcp-widgets";
const EXPECTED_PLUGIN_VERSION = "0.2.8";
const EXPECTED_MODULE_SHA256 = Object.freeze({
  builder: "0b86883810cf23c81f563a48a7708283a6a4cadf11228ad261efb64426ef728e",
  extractor: "77e78593d8971563b6394334128dd3ee243f2182d3d087591d3a5499ece9283b",
  verifier: "b495c4cc34113fb2918118eac302f4e7e2152c1b1b9b63e3646cea97ecbf9b3f",
  readerAsset: "6c5ed0d30e002c9e0cbbb11933296e4f831fb4de5933e21f34bb44f3999be7b5",
  packageJson: "856f7cc8a4737576167795ca11431a5b64721caae2137f5bc614d4e7ea6855fa",
  runtimePart001: "154f1d561bab28174f88d71ae710599709e5a5eda64dee08b025a699449dbbfc",
  runtimePart002: "9459ed4b76bd825daba9637564dd3122ca617a88dd16531d3e9bca5c7ced080e",
  runtimePart003: "cdfe2e6787faa61d37f043df6996d5f988c07360bb2b3d8af2aa3c18e8db3ac0",
  validatorServer: "eff59c6085d2ab6b6153c80a03749e764e160f8c6711da8433f7bd6762e1db66",
  browserHelpers: "84aa4f8a2a11376ebee6942d2d7e10a083d16c65dfdb73114d7b34e51c27f69d",
  browserCli: "aac3b12fc12c7ad2e044533f881791dd3b23bd0ee31ddf6682dd7f6de99e6596",
});
const EXPECTED_BROWSER_SHA256 =
  "7b8e92dca0acf9c24b5974507b3031d6bd18cc009cd431c2595e521430ea747a";
const WRAPPER_PATH = fileURLToPath(import.meta.url);

function usage() {
  return [
    "Usage:",
    "  node deliver_portable_artifact_strict.mjs \\",
    "    --input <artifact.json> \\",
    "    --output <report.html> \\",
    "    --receipt <delivery_receipt.json>",
    "",
    "The renderer is pinned to the repository-audited data-analytics plugin build.",
  ].join("\n");
}

function parseArguments(argv) {
  const values = {};
  const supported = new Set(["input", "output", "receipt"]);
  for (let index = 0; index < argv.length; index += 1) {
    const argument = argv[index];
    if (argument === "--help" || argument === "-h") return { help: true };
    if (!argument.startsWith("--")) {
      throw new Error(`Unexpected argument: ${argument}\n${usage()}`);
    }
    const key = argument.slice(2);
    if (!supported.has(key)) throw new Error(`Unknown argument: ${argument}\n${usage()}`);
    const value = argv[index + 1];
    if (!value || value.startsWith("--")) {
      throw new Error(`Missing value for ${argument}.\n${usage()}`);
    }
    if (values[key] !== undefined) {
      throw new Error(`${argument} may only be specified once.`);
    }
    values[key] = value;
    index += 1;
  }
  for (const key of ["input", "output", "receipt"]) {
    if (!values[key]) throw new Error(`--${key} is required.\n${usage()}`);
  }
  return values;
}

function sha256(value) {
  return createHash("sha256").update(value).digest("hex");
}

function fileFingerprint(path) {
  const bytes = readFileSync(path);
  return {
    path,
    bytes: bytes.byteLength,
    sha256: sha256(bytes),
  };
}

function countOccurrences(text, needle) {
  if (!needle) return 0;
  let count = 0;
  let offset = 0;
  while ((offset = text.indexOf(needle, offset)) !== -1) {
    count += 1;
    offset += needle.length;
  }
  return count;
}

function applyTopbarOverflowFix(html) {
  if (countOccurrences(html, PATCH_MARKER) !== 0) {
    throw new Error("Portable HTML already contains the strict top-bar patch marker.");
  }
  if (countOccurrences(html, "</head>") !== 1) {
    throw new Error("Portable HTML must contain exactly one closing head tag.");
  }
  const style = `<style ${PATCH_MARKER}="true">${PATCH_CSS}</style>`;
  const patched = html.replace("</head>", `${style}</head>`);
  if (countOccurrences(patched, PATCH_MARKER) !== 1) {
    throw new Error("Strict top-bar patch injection did not produce exactly one marker.");
  }
  return patched;
}

function requireRegularFile(path, label) {
  if (!existsSync(path) || !statSync(path).isFile()) {
    throw new Error(`${label} is not a regular file: ${path}`);
  }
}

function requireUnusedPath(path, label) {
  if (existsSync(path)) {
    throw new Error(`${label} already exists; refusing to overwrite: ${path}`);
  }
}

function staticChartContract(staticCharts, expectedChartIds) {
  const entries = Object.entries(staticCharts);
  if (entries.length !== expectedChartIds.length) {
    throw new Error(
      `Static SVG extraction returned ${entries.length} charts; expected ${expectedChartIds.length}.`,
    );
  }
  const extractedChartIds = entries.map(([, chart]) => chart?.chartId).sort();
  const expectedSorted = [...expectedChartIds].sort();
  if (JSON.stringify(extractedChartIds) !== JSON.stringify(expectedSorted)) {
    throw new Error(
      `Static SVG chart ids do not match visible chart blocks: ` +
        `actual=${JSON.stringify(extractedChartIds)} expected=${JSON.stringify(expectedSorted)}`,
    );
  }
  for (const [chartKey, chart] of entries) {
    for (const theme of ["light", "dark"]) {
      if (typeof chart?.[theme]?.svg !== "string" || !chart[theme].svg.startsWith("<svg")) {
        throw new Error(`Static chart ${chartKey} is missing its ${theme} SVG variant.`);
      }
    }
  }
  return {
    chartKeys: entries.map(([chartKey]) => chartKey).sort(),
    chartIds: extractedChartIds,
    count: entries.length,
    themeVariants: entries.length * 2,
  };
}

function visibleChartIds(artifact) {
  const ids = (artifact?.manifest?.blocks ?? [])
    .filter((block) => block?.type === "chart")
    .map((block) => block.chartId);
  if (ids.length === 0) {
    throw new Error("Report must contain at least one visible chart block.");
  }
  if (ids.some((id) => typeof id !== "string" || !id)) {
    throw new Error("Every visible chart block must declare chartId.");
  }
  if (new Set(ids).size !== ids.length) {
    throw new Error("Visible chart blocks must reference unique chart ids for strict SVG accounting.");
  }
  return ids;
}

function hasSourceBackedContent(artifact) {
  const manifest = artifact?.manifest ?? {};
  const visible = {
    card: new Set(
      (manifest.blocks ?? [])
        .filter((block) => block?.type === "metric-strip")
        .flatMap((block) => block.cardIds ?? []),
    ),
    chart: new Set(
      (manifest.blocks ?? [])
        .filter((block) => block?.type === "chart")
        .map((block) => block.chartId),
    ),
    table: new Set(
      (manifest.blocks ?? [])
        .filter((block) => block?.type === "table")
        .map((block) => block.tableId),
    ),
  };
  return [
    ["card", manifest.cards ?? []],
    ["chart", manifest.charts ?? []],
    ["table", manifest.tables ?? []],
  ].some(([kind, items]) =>
    items.some(
      (item) => visible[kind].has(item?.id) && Boolean(item?.sourceId || item?.source),
    ),
  );
}

function assertOfficialVerification(verification, artifact, expectedChartCount) {
  if (!verification?.ok) throw new Error("Official portable verification did not return ok=true.");
  if (verification.counts?.charts !== expectedChartCount) {
    throw new Error(
      `Official verification counted ${verification.counts?.charts} charts; ` +
        `expected ${expectedChartCount}.`,
    );
  }
  if (JSON.stringify(verification.viewports) !== JSON.stringify(EXPECTED_VIEWPORTS)) {
    throw new Error(
      `Official verification viewports differ: ` +
        `actual=${JSON.stringify(verification.viewports)} ` +
        `expected=${JSON.stringify(EXPECTED_VIEWPORTS)}`,
    );
  }
  if (
    hasSourceBackedContent(artifact) &&
    (verification.sourceDialog !== "passed" ||
      !["keyboard_menu_semantic_click", "semantic_click"].includes(
        verification.sourceInteraction,
      ))
  ) {
    throw new Error(
      `Source interaction verification failed: dialog=${verification.sourceDialog} ` +
        `interaction=${verification.sourceInteraction}`,
    );
  }
}

function writeJson(path, value) {
  writeFileSync(path, `${JSON.stringify(value, null, 2)}\n`, "utf8");
}

function compactError(error) {
  const text = error?.message ?? String(error);
  return text.length > 1_000 ? `${text.slice(0, 997)}...` : text;
}

async function run() {
  const startedAt = performance.now();
  const startedAtIso = new Date().toISOString();
  let stage = "arguments";
  let inputPath;
  let outputPath;
  let receiptPath;
  let candidateOutputPath;
  let candidateReceiptPath;
  let publishedOutput = false;
  let inputFingerprint;
  let pluginRoot;
  let packageMetadata;
  let browserExecutable;
  let browserFingerprint;
  const failureArchivePaths = [];

  try {
    const parsed = parseArguments(process.argv.slice(2));
    if (parsed.help) {
      process.stdout.write(`${usage()}\n`);
      return;
    }

    inputPath = resolve(parsed.input);
    outputPath = resolve(parsed.output);
    receiptPath = resolve(parsed.receipt);
    pluginRoot = resolve(DEFAULT_PLUGIN_ROOT);

    if (new Set([inputPath, outputPath, receiptPath]).size !== 3) {
      throw new Error("--input, --output, and --receipt must resolve to three distinct paths.");
    }
    stage = "path_validation";
    requireRegularFile(inputPath, "Artifact input");
    requireUnusedPath(outputPath, "Report output");
    requireUnusedPath(receiptPath, "Delivery receipt");

    const modulePaths = {
      builder: join(
        pluginRoot,
        "skills",
        "build-report",
        "scripts",
        "build_portable_artifact.mjs",
      ),
      extractor: join(
        pluginRoot,
        "skills",
        "build-report",
        "scripts",
        "extract_portable_chart_svgs.mjs",
      ),
      verifier: join(
        pluginRoot,
        "skills",
        "build-report",
        "scripts",
        "verify_portable_artifact.mjs",
      ),
      readerAsset: join(pluginRoot, "assets", "portable-artifact-reader.html"),
      packageJson: join(pluginRoot, "package.json"),
      runtimePart001: join(
        pluginRoot,
        "assets",
        "portable-artifact-reader.html.gz.b64.part001",
      ),
      runtimePart002: join(
        pluginRoot,
        "assets",
        "portable-artifact-reader.html.gz.b64.part002",
      ),
      runtimePart003: join(
        pluginRoot,
        "assets",
        "portable-artifact-reader.html.gz.b64.part003",
      ),
      validatorServer: join(pluginRoot, "mcp", "server.cjs"),
      browserHelpers: join(
        pluginRoot,
        "skills",
        "build-report",
        "scripts",
        "portable_browser_helpers.mjs",
      ),
      browserCli: join(
        pluginRoot,
        "skills",
        "build-report",
        "scripts",
        "portable_browser_cli.mjs",
      ),
    };
    const expectedRuntimePartNames = [
      "portable-artifact-reader.html.gz.b64.part001",
      "portable-artifact-reader.html.gz.b64.part002",
      "portable-artifact-reader.html.gz.b64.part003",
    ];
    const actualRuntimePartNames = readdirSync(join(pluginRoot, "assets"))
      .filter((name) => name.startsWith("portable-artifact-reader.html.gz.b64.part"))
      .sort();
    if (
      actualRuntimePartNames.length !== expectedRuntimePartNames.length ||
      actualRuntimePartNames.some(
        (name, index) => name !== expectedRuntimePartNames[index],
      )
    ) {
      throw new Error(
        `Embedded reader runtime part set mismatch: ` +
          `${JSON.stringify(actualRuntimePartNames)}; expected ` +
          `${JSON.stringify(expectedRuntimePartNames)}.`,
      );
    }
    for (const [label, path] of Object.entries(modulePaths)) {
      requireRegularFile(path, `Official plugin ${label}`);
      const actualSha256 = fileFingerprint(path).sha256;
      if (actualSha256 !== EXPECTED_MODULE_SHA256[label]) {
        throw new Error(
          `Official plugin ${label} checksum mismatch: ${actualSha256}; expected ` +
            `${EXPECTED_MODULE_SHA256[label]}.`,
        );
      }
    }
    packageMetadata = JSON.parse(readFileSync(modulePaths.packageJson, "utf8"));
    if (
      packageMetadata.name !== EXPECTED_PLUGIN_NAME ||
      packageMetadata.version !== EXPECTED_PLUGIN_VERSION
    ) {
      throw new Error(
        `Unsupported data-analytics renderer package: ` +
          `${packageMetadata.name}@${packageMetadata.version}; expected ` +
          `${EXPECTED_PLUGIN_NAME}@${EXPECTED_PLUGIN_VERSION}.`,
      );
    }
    browserExecutable = resolve(DEFAULT_BROWSER_EXECUTABLE);
    requireRegularFile(browserExecutable, "Pinned Chromium headless-shell");
    browserFingerprint = fileFingerprint(browserExecutable);
    if (browserFingerprint.sha256 !== EXPECTED_BROWSER_SHA256) {
      throw new Error(
        `Pinned Chromium checksum mismatch: ${browserFingerprint.sha256}; expected ` +
          `${EXPECTED_BROWSER_SHA256}.`,
      );
    }

    mkdirSync(dirname(outputPath), { recursive: true });
    mkdirSync(dirname(receiptPath), { recursive: true });
    candidateOutputPath = join(
      dirname(outputPath),
      `.${basename(outputPath)}.tmp-${process.pid}-${randomUUID()}`,
    );
    candidateReceiptPath = join(
      dirname(receiptPath),
      `.${basename(receiptPath)}.tmp-${process.pid}-${randomUUID()}`,
    );

    stage = "official_module_import";
    const [
      { buildPortableArtifact },
      { extractPortableChartSvgs },
      { verifyPortableArtifact },
    ] = await Promise.all([
      import(pathToFileURL(modulePaths.builder).href),
      import(pathToFileURL(modulePaths.extractor).href),
      import(pathToFileURL(modulePaths.verifier).href),
    ]);
    if (
      typeof buildPortableArtifact !== "function" ||
      typeof extractPortableChartSvgs !== "function" ||
      typeof verifyPortableArtifact !== "function"
    ) {
      throw new Error("Official plugin modules do not expose the required functions.");
    }

    stage = "artifact_read";
    const artifactBytes = readFileSync(inputPath);
    inputFingerprint = {
      path: inputPath,
      bytes: artifactBytes.byteLength,
      sha256: sha256(artifactBytes),
    };
    const artifact = JSON.parse(artifactBytes.toString("utf8"));
    const expectedChartIds = visibleChartIds(artifact);

    stage = "canonical_validation_and_initial_package";
    let html = applyTopbarOverflowFix(buildPortableArtifact(artifact));
    writeFileSync(candidateOutputPath, html, "utf8");

    stage = "official_static_svg_extraction";
    const staticCharts = await extractPortableChartSvgs({
      browserExecutable,
      htmlPath: candidateOutputPath,
    });
    const staticSvg = staticChartContract(staticCharts, expectedChartIds);

    stage = "canonical_repackage_with_static_svg";
    html = applyTopbarOverflowFix(buildPortableArtifact(artifact, { staticCharts }));
    const patchMarkerCount = countOccurrences(html, PATCH_MARKER);
    const embeddedSvgVariants = countOccurrences(html, "portable-static-chart-svg");
    if (patchMarkerCount !== 1) {
      throw new Error(`Final HTML contains ${patchMarkerCount} top-bar patch markers; expected 1.`);
    }
    if (embeddedSvgVariants !== staticSvg.themeVariants) {
      throw new Error(
        `Final HTML contains ${embeddedSvgVariants} static SVG variants; ` +
          `expected ${staticSvg.themeVariants}.`,
      );
    }
    writeFileSync(candidateOutputPath, html, "utf8");

    stage = "official_desktop_and_narrow_verification";
    const verification = await verifyPortableArtifact({
      artifactPath: inputPath,
      browserExecutable,
      htmlPath: candidateOutputPath,
      timeoutMs: 30_000,
    });
    assertOfficialVerification(verification, artifact, expectedChartIds.length);

    stage = "receipt_prepare";
    const outputFingerprint = fileFingerprint(candidateOutputPath);
    outputFingerprint.path = outputPath;
    const receipt = {
      version: 1,
      ok: true,
      generatedAt: new Date().toISOString(),
      command: {
        argv: process.argv,
        node: process.version,
      },
      wrapper: fileFingerprint(WRAPPER_PATH),
      input: inputFingerprint,
      output: outputFingerprint,
      plugin: {
        root: pluginRoot,
        name: packageMetadata.name ?? null,
        version: packageMetadata.version ?? null,
        files: Object.fromEntries(
          Object.entries(modulePaths).map(([label, path]) => [label, fileFingerprint(path)]),
        ),
      },
      browser: browserFingerprint,
      patch: {
        marker: PATCH_MARKER,
        markerCount: patchMarkerCount,
        css: PATCH_CSS,
        sha256: sha256(PATCH_CSS),
        scope: "single deterministic override outside the canonical embedded artifact payload",
      },
      stages: {
        canonicalValidation: "passed",
        package: "passed",
        staticSvgExtraction: "passed",
        officialBrowserVerification: "passed",
        atomicOutputRename: "pending",
      },
      staticSvg: {
        ...staticSvg,
        embeddedThemeVariants: embeddedSvgVariants,
      },
      verification: {
        counts: verification.counts,
        sourceDialog: verification.sourceDialog,
        sourceInteraction: verification.sourceInteraction,
        viewports: verification.viewports,
        timings: verification.timings,
      },
      timing: {
        startedAt: startedAtIso,
        elapsedMsBeforePublish: Math.round((performance.now() - startedAt) * 10) / 10,
      },
    };
    writeJson(candidateReceiptPath, receipt);

    stage = "atomic_publish";
    renameSync(candidateOutputPath, outputPath);
    publishedOutput = true;
    receipt.stages.atomicOutputRename = "passed";
    receipt.timing.completedAt = new Date().toISOString();
    receipt.timing.totalMs = Math.round((performance.now() - startedAt) * 10) / 10;
    writeJson(candidateReceiptPath, receipt);
    try {
      renameSync(candidateReceiptPath, receiptPath);
    } catch (error) {
      renameSync(outputPath, candidateOutputPath);
      publishedOutput = false;
      throw error;
    }

    process.stdout.write(
      `${JSON.stringify({
        ok: true,
        input: inputPath,
        output: outputPath,
        receipt: receiptPath,
        stages: receipt.stages,
        staticSvg: {
          charts: staticSvg.count,
          themeVariants: staticSvg.themeVariants,
        },
        verification: receipt.verification,
      })}\n`,
    );
  } catch (error) {
    if (publishedOutput && outputPath && candidateOutputPath) {
      try {
        renameSync(outputPath, candidateOutputPath);
        publishedOutput = false;
      } catch {
        // Preserve the original failure; stdout/stderr still reports whether final output remains.
      }
    }
    const existingCandidates = [candidateOutputPath, candidateReceiptPath].filter(
      (candidate) => candidate && existsSync(candidate),
    );
    if (existingCandidates.length > 0) {
      try {
        const archiveParent = join(
          dirname(receiptPath ?? outputPath),
          "archive",
          "failed_delivery_staging",
          `${new Date().toISOString().replaceAll(":", "")}-${randomUUID()}`,
        );
        mkdirSync(archiveParent, { recursive: true });
        for (const candidate of existingCandidates) {
          const archived = join(archiveParent, basename(candidate));
          renameSync(candidate, archived);
          failureArchivePaths.push(archived);
        }
      } catch {
        // Preserve candidates in place when archival itself fails.
      }
    }

    const failure = {
      version: 1,
      ok: false,
      generatedAt: new Date().toISOString(),
      stage,
      error: compactError(error),
      input: inputFingerprint ?? (inputPath ? { path: inputPath } : null),
      output: outputPath ?? null,
      outputExists: Boolean(outputPath && existsSync(outputPath)),
      pluginRoot: pluginRoot ?? null,
      plugin:
        packageMetadata == null
          ? null
          : { name: packageMetadata.name ?? null, version: packageMetadata.version ?? null },
      failureArchivePaths,
      patch: {
        marker: PATCH_MARKER,
        css: PATCH_CSS,
        sha256: sha256(PATCH_CSS),
      },
      timing: {
        startedAt: startedAtIso,
        totalMs: Math.round((performance.now() - startedAt) * 10) / 10,
      },
    };
    if (receiptPath && !existsSync(receiptPath)) {
      try {
        mkdirSync(dirname(receiptPath), { recursive: true });
        writeJson(receiptPath, failure);
      } catch {
        // stderr remains the fail-closed audit surface when receipt publication also fails.
      }
    }
    process.stderr.write(`${JSON.stringify(failure)}\n`);
    process.exitCode = 1;
  }
}

await run();
