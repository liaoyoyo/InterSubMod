#!/usr/bin/env node

import { createHash } from "node:crypto";
import { copyFileSync, readFileSync, statSync, writeFileSync } from "node:fs";
import { resolve } from "node:path";

import { buildPortableArtifact } from "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/scripts/build_portable_artifact.mjs";
import { extractPortableChartSvgs } from "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/scripts/extract_portable_chart_svgs.mjs";
import { verifyPortableArtifact } from "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599/skills/build-report/scripts/verify_portable_artifact.mjs";

const root = "/big7_disk/liaoyoyo2001/InterSubMod";
const artifactPath = resolve(
  root,
  "research/20260725_methyl_alt_ref_topology_overlay_validation/artifact.json",
);
const outputPath = resolve(
  root,
  "research/20260725_methyl_alt_ref_topology_overlay_validation/" +
    "20260725_ALTREF甲基差異與latest候選拓撲對應驗證_01.html",
);
const receiptPath = resolve(
  root,
  "research/20260725_methyl_alt_ref_topology_overlay_validation/" +
    "results/html_delivery_receipt.json",
);
const dynamicPath = `/tmp/alt-ref-topology-${process.pid}-dynamic.html`;
const candidatePath = `/tmp/alt-ref-topology-${process.pid}-verified.html`;
const screenshotPath = `/tmp/alt-ref-topology-${process.pid}-verification-failure.png`;
const topbarFix = `
<style id="intersubmod-portable-topbar-scrollbar-fix">
/* Packaged reader 0.2.8 uses 100vw for the sticky bar. Classic scrollbars make
   100vw wider than documentElement.clientWidth, so use the shell width. */
.analytics-top-bar {
  width: 100% !important;
  margin-right: 0 !important;
  margin-left: 0 !important;
}
</style>`;

process.stderr.write("[1/5] validate and build dynamic portable HTML\n");
const artifact = JSON.parse(readFileSync(artifactPath, "utf8"));
writeFileSync(dynamicPath, buildPortableArtifact(artifact), "utf8");
const skipStaticCharts = process.env.SKIP_STATIC_CHARTS === "1";
process.stderr.write(
  skipStaticCharts
    ? "[2/5] retain the offline dynamic chart runtime (static extraction explicitly skipped)\n"
    : "[2/5] extract light/dark static chart SVGs\n",
);
const staticCharts = skipStaticCharts
  ? {}
  : await extractPortableChartSvgs({
      htmlPath: dynamicPath,
      readyTimeoutMs: 60_000,
      actionTimeoutMs: 30_000,
    });
process.stderr.write("[3/5] rebuild portable HTML and inject scoped top-bar fix\n");
const packaged = buildPortableArtifact(artifact, { staticCharts });
if (!packaged.includes("</head>")) throw new Error("portable HTML has no closing head");
const patched = packaged.replace("</head>", `${topbarFix}\n</head>`);
if ((patched.match(/intersubmod-portable-topbar-scrollbar-fix/g) ?? []).length !== 1) {
  throw new Error("topbar fix was not injected exactly once");
}
writeFileSync(candidatePath, patched, "utf8");
process.stderr.write("[4/5] verify desktop/mobile layout and source interaction\n");
const verification = await verifyPortableArtifact({
  artifactPath,
  htmlPath: candidatePath,
  readyTimeoutMs: 5000,
  actionTimeoutMs: 2500,
  timeoutMs: 20000,
  screenshotPath,
});
process.stderr.write("[5/5] replace output only after the candidate passes verification\n");
copyFileSync(candidatePath, outputPath);
const sha256 = (path) =>
  createHash("sha256").update(readFileSync(path)).digest("hex");
const staticChartCount = Object.keys(staticCharts).length;
const allPass =
  verification.ok === true &&
  verification.counts?.blocks === 32 &&
  verification.counts?.charts === 6 &&
  verification.counts?.metrics === 5 &&
  verification.counts?.tables === 6 &&
  verification.sourceDialog === "passed" &&
  JSON.stringify(verification.viewports) === JSON.stringify([1440, 390]) &&
  (skipStaticCharts || staticChartCount === 6);
if (!allPass) {
  throw new Error(
    `portable delivery contract failed: ${JSON.stringify({
      staticChartCount,
      verification,
    })}`,
  );
}
const receipt = {
  schema_name: "intersubmod.portable_html_delivery_receipt",
  schema_version: "1.1.0",
  created_at: new Date().toISOString(),
  task_type: "B_comprehensive_validation",
  goals: ["G3", "G4"],
  input: {
    path: artifactPath,
    sha256: sha256(artifactPath),
  },
  command:
    "CHROMIUM_EXECUTABLE_PATH=" +
    String(process.env.CHROMIUM_EXECUTABLE_PATH || "<playwright-default>") +
    " node research/20260725_methyl_alt_ref_topology_overlay_validation/" +
    "scripts/deliver_portable_with_topbar_fix.mjs",
  output: {
    path: outputPath,
    sha256: sha256(outputPath),
    size_bytes: statSync(outputPath).size,
  },
  packaging: {
    portable_builder: "data-analytics build-report 0.2.8",
    static_chart_count: staticChartCount,
    chart_delivery: skipStaticCharts
      ? "offline_dynamic_runtime"
      : "static_svg_with_dynamic_fallback",
    css_workaround: {
      id: "intersubmod-portable-topbar-scrollbar-fix",
      reason:
        "Packaged reader 0.2.8 uses width:100vw for the sticky top bar; " +
        "classic scrollbars made scrollWidth exceed clientWidth by 8px.",
      scope: ".analytics-top-bar width and horizontal margins only",
      analytical_data_changed: false,
    },
  },
  verification: {
    ok: verification.ok,
    counts: verification.counts,
    viewports: verification.viewports,
    source_dialog: verification.sourceDialog,
    source_interaction: verification.sourceInteraction,
  },
  all_pass: allPass,
};
writeFileSync(receiptPath, `${JSON.stringify(receipt, null, 2)}\n`, "utf8");
process.stdout.write(`${JSON.stringify({
  ok: true,
  artifact: artifactPath,
  output: outputPath,
  receipt: receiptPath,
  staticChartCount,
  chartDelivery: skipStaticCharts ? "offline_dynamic_runtime" : "static_svg_with_dynamic_fallback",
  cssWorkaround: "portable_reader_0.2.8_100vw_classic_scrollbar",
  verification,
})}\n`);
