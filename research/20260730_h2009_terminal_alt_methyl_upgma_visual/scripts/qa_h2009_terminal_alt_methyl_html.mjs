#!/usr/bin/env node
/** Independent Playwright QA for the H2009 terminal-ALT portable HTML report. */

import { readFileSync, writeFileSync } from "node:fs";
import { pathToFileURL } from "node:url";

import { chromium } from "/bip7_disk/liaoyoyo2001/miniconda3/lib/python3.9/site-packages/playwright/driver/package/index.mjs";

const WORK = "/big7_disk/liaoyoyo2001/InterSubMod/research/20260730_h2009_terminal_alt_methyl_upgma_visual";
const ARTIFACT = `${WORK}/artifact.json`;
const HTML = `${WORK}/20260730_H2009末端ALT甲基距離與UPGMA圖解_01.html`;
const RESULTS = `${WORK}/results`;
const RECEIPT = `${RESULTS}/20260730_H2009末端ALT甲基UPGMA_HTML_Playwright_QA_01.json`;
const TOP_SCREENSHOT = `${RESULTS}/02_H2009_terminal_ALT_HTML_top_preview.png`;
const FIGURE_SCREENSHOT = `${RESULTS}/03_H2009_terminal_ALT_UPGMA_preview.png`;
const NATIVE_HEATMAP_SCREENSHOT = `${RESULTS}/04_H2009_terminal_ALT_native_heatmap_preview.png`;
const CHROMIUM = "/bip7_disk/liaoyoyo2001/.cache/ms-playwright/chromium-1223/chrome-linux64/chrome";

function expectedCounts(artifact) {
  const blocks = artifact.manifest.blocks;
  const metricCards = blocks
    .filter((block) => block.type === "metric-strip")
    .reduce((sum, block) => sum + (block.cardIds?.length ?? 0), 0);
  return {
    blocks: blocks.filter((block) => block.type !== "metric-strip").length + metricCards,
    charts: blocks.filter((block) => block.type === "chart").length,
    html: blocks.filter((block) => block.type === "html").length,
    metrics: metricCards,
    tables: blocks.filter((block) => block.type === "table").length,
  };
}

async function qaViewport(browser, artifact, viewport, capture) {
  const consoleErrors = [];
  const pageErrors = [];
  const externalRequests = [];
  const context = await browser.newContext({
    viewport: { width: viewport.width, height: viewport.height },
    colorScheme: "light",
    reducedMotion: "reduce",
  });
  const page = await context.newPage();
  page.setDefaultTimeout(60_000);
  page.on("console", (message) => {
    if (message.type() === "error") consoleErrors.push(message.text());
  });
  page.on("pageerror", (error) => pageErrors.push(String(error)));
  page.on("request", (request) => {
    if (/^(?:https?|wss?):\/\//.test(request.url())) externalRequests.push(request.url());
  });

  await page.goto(pathToFileURL(HTML).href, { waitUntil: "load", timeout: 120_000 });
  await page.waitForFunction(
    () => {
      const state = document.documentElement.dataset.dataAnalyticsPortableReader ?? "";
      return ["ready", "failed", "missing-runtime", "unsupported"].includes(state);
    },
    null,
    { timeout: 120_000 },
  );
  const readerState = await page.evaluate(
    () => document.documentElement.dataset.dataAnalyticsPortableReader ?? "",
  );
  if (readerState !== "ready") {
    throw new Error(`Portable reader entered terminal state ${readerState} at ${viewport.name}`);
  }
  await page.waitForTimeout(750);

  const title = (await page.locator("#data-analytics-portable-reader h1").first().textContent()).trim();
  if (title !== artifact.manifest.title) {
    throw new Error(`Title mismatch at ${viewport.name}: ${title} != ${artifact.manifest.title}`);
  }

  const bodyText = await page.locator("#data-analytics-portable-reader").innerText();
  const requiredText = [
    "DEMO / PARTIAL SCOPE",
    "81＋5",
    "2.894439",
    "6/6",
    "UPGMA",
  ];
  const missingText = requiredText.filter((text) => !bodyText.includes(text));
  if (missingText.length) {
    throw new Error(`Missing explanatory text at ${viewport.name}: ${missingText.join(", ")}`);
  }

  const actualCounts = {
    blocks: await page.locator("#data-analytics-portable-reader [data-analytics-layout-item]").count(),
    charts: await page.locator("#data-analytics-portable-reader .chart-frame").count(),
    html: await page.locator("#data-analytics-portable-reader iframe.report-html-frame").count(),
    metrics: await page.locator("#data-analytics-portable-reader .report-metric-card").count(),
    tables: await page.locator("#data-analytics-portable-reader .table-card").count(),
  };
  const expected = expectedCounts(artifact);
  if (JSON.stringify(actualCounts) !== JSON.stringify(expected)) {
    throw new Error(
      `Rendered counts differ at ${viewport.name}: ${JSON.stringify(actualCounts)} != ${JSON.stringify(expected)}`,
    );
  }

  const geometry = await page.evaluate(() => {
    const selectors = [
      "[data-analytics-layout-item]",
      ".chart-frame",
      ".table-card",
      ".report-metric-card",
    ];
    const bad = [];
    for (const selector of selectors) {
      for (const node of document.querySelectorAll(`#data-analytics-portable-reader ${selector}`)) {
        const rect = node.getBoundingClientRect();
        if (!(rect.width > 0 && rect.height > 0)) {
          bad.push({ selector, width: rect.width, height: rect.height });
        }
      }
    }
    return {
      bad,
      clientWidth: document.documentElement.clientWidth,
      scrollWidth: document.documentElement.scrollWidth,
    };
  });
  if (geometry.bad.length) {
    throw new Error(`Zero-size content at ${viewport.name}: ${JSON.stringify(geometry.bad.slice(0, 5))}`);
  }
  if (geometry.scrollWidth > geometry.clientWidth + 1) {
    throw new Error(`Horizontal page overflow at ${viewport.name}: ${JSON.stringify(geometry)}`);
  }

  const chartSvgs = await page
    .locator("#data-analytics-portable-reader .chart-frame svg.recharts-surface")
    .count();
  const chartRenderDetails = await page
    .locator("#data-analytics-portable-reader .chart-frame")
    .evaluateAll((frames) =>
      frames.map((frame) => {
        const plot = frame.querySelector(".chart-plot");
        const rect = plot?.getBoundingClientRect();
        return {
          childElements: plot?.childElementCount ?? 0,
          canvases: plot?.querySelectorAll("canvas").length ?? 0,
          height: rect?.height ?? 0,
          rects: plot?.querySelectorAll("rect").length ?? 0,
          svgs: plot?.querySelectorAll("svg").length ?? 0,
          textLength: plot?.textContent?.trim().length ?? 0,
          width: rect?.width ?? 0,
        };
      }),
    );
  if (
    chartRenderDetails.some(
      (detail) =>
        detail.childElements < 1 || detail.width <= 0 || detail.height <= 0,
    )
  ) {
    throw new Error(`Native heatmap did not render at ${viewport.name}: ${JSON.stringify(chartRenderDetails)}`);
  }

  const iframeLocator = page.locator("iframe.report-html-frame");
  if ((await iframeLocator.count()) !== 1) {
    throw new Error(`Expected one UPGMA iframe at ${viewport.name}`);
  }
  const iframeHandle = await iframeLocator.first().elementHandle();
  const frame = await iframeHandle.contentFrame();
  const image = frame.locator("img");
  await image.waitFor({ state: "visible", timeout: 60_000 });
  const imageStatus = await image.evaluate((img) => ({
    complete: img.complete,
    naturalWidth: img.naturalWidth,
    naturalHeight: img.naturalHeight,
    alt: img.alt,
    dataUri: img.src.startsWith("data:image/png;base64,"),
    clientWidth: img.clientWidth,
  }));
  if (
    !imageStatus.complete ||
    imageStatus.naturalWidth !== 2336 ||
    imageStatus.naturalHeight !== 1804 ||
    !imageStatus.dataUri ||
    imageStatus.clientWidth < 1
  ) {
    throw new Error(`Embedded UPGMA figure failed at ${viewport.name}: ${JSON.stringify(imageStatus)}`);
  }
  const iframeGeometry = await frame.locator("body").evaluate((body) => ({
    clientWidth: body.clientWidth,
    scrollWidth: body.scrollWidth,
  }));
  const iframeText = await frame.locator("body").innerText();
  if (!iframeText.includes("不是 mutation ancestry")) {
    throw new Error(`UPGMA claim-ceiling caption is missing at ${viewport.name}`);
  }
  if (iframeGeometry.scrollWidth > iframeGeometry.clientWidth + 1) {
    throw new Error(`UPGMA iframe overflow at ${viewport.name}: ${JSON.stringify(iframeGeometry)}`);
  }

  let sourceInteraction = "not_run";
  if (viewport.name === "desktop") {
    const sourceButton = page
      .locator(
        '#data-analytics-portable-reader button[data-artifact-action="open-options"]' +
          '[data-artifact-has-source="true"]',
      )
      .first();
    await sourceButton.scrollIntoViewIfNeeded();
    await sourceButton.evaluate((button) => button.click());
    const sourceAction = page
      .locator('[role="menu"] [data-artifact-action="view-source"]')
      .first();
    await sourceAction.waitFor({ state: "visible" });
    await sourceAction.evaluate((button) => button.click());
    const dialog = page.locator('[data-artifact-dialog="source"]').first();
    await dialog.waitFor({ state: "visible" });
    if ((await dialog.getByRole("tab", { name: "Overview", exact: true }).count()) !== 1) {
      throw new Error("Source dialog lacks its Overview tab");
    }
    await page.keyboard.press("Escape");
    sourceInteraction = "passed";
  }

  if (capture) {
    await page.evaluate(() => window.scrollTo(0, 0));
    await page.waitForTimeout(250);
    await page.screenshot({ path: TOP_SCREENSHOT, fullPage: false });
    await iframeLocator.first().evaluate((node) => node.scrollIntoView({ block: "start" }));
    await page.waitForTimeout(500);
    await page.screenshot({ path: FIGURE_SCREENSHOT, fullPage: false });
    const chartFrame = page.locator("#data-analytics-portable-reader .chart-frame").first();
    await chartFrame.scrollIntoViewIfNeeded();
    await page.waitForTimeout(300);
    await chartFrame.screenshot({ path: NATIVE_HEATMAP_SCREENSHOT });
  }

  if (externalRequests.length) {
    throw new Error(`External requests at ${viewport.name}: ${externalRequests.join(", ")}`);
  }
  if (consoleErrors.length || pageErrors.length) {
    throw new Error(
      `Browser errors at ${viewport.name}: console=${JSON.stringify(consoleErrors)}, ` +
        `page=${JSON.stringify(pageErrors)}`,
    );
  }

  const result = {
    name: viewport.name,
    viewport: { width: viewport.width, height: viewport.height },
    reader_state: readerState,
    title,
    required_text_present: true,
    counts: actualCounts,
    chart_svgs: chartSvgs,
    chart_render_details: chartRenderDetails,
    embedded_upgma: imageStatus,
    upgma_claim_ceiling_present: true,
    iframe_overflow: iframeGeometry,
    page_overflow: geometry,
    source_interaction: sourceInteraction,
    external_requests: externalRequests,
    console_errors: consoleErrors,
    page_errors: pageErrors,
  };
  await context.close();
  return result;
}

async function main() {
  const artifact = JSON.parse(readFileSync(ARTIFACT, "utf8"));
  const browser = await chromium.launch({
    headless: true,
    executablePath: CHROMIUM,
    args: ["--no-sandbox", "--disable-dev-shm-usage", "--disable-gpu"],
  });
  let viewports;
  try {
    viewports = [
      await qaViewport(browser, artifact, { name: "desktop", width: 1440, height: 1000 }, true),
      await qaViewport(browser, artifact, { name: "mobile", width: 390, height: 844 }, false),
    ];
  } finally {
    await browser.close();
  }

  const receipt = {
    ok: true,
    created_at: new Date().toISOString(),
    html: HTML,
    artifact: ARTIFACT,
    browser: CHROMIUM,
    viewports,
    screenshots: [TOP_SCREENSHOT, FIGURE_SCREENSHOT, NATIVE_HEATMAP_SCREENSHOT],
    note:
      "Canonical builder passed validation, packaging, and structural verification. " +
      "The host dump-dom extractor exited early; this independent Node Playwright run verifies the actual dynamic " +
      "heatmap, tables, metric cards, sandboxed UPGMA image, source dialog, desktop/mobile geometry, " +
      "console/page errors, and external requests.",
  };
  writeFileSync(RECEIPT, `${JSON.stringify(receipt, null, 2)}\n`, "utf8");
  process.stdout.write(`${JSON.stringify(receipt, null, 2)}\n`);
}

await main();
