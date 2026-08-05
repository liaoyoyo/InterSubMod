#!/usr/bin/env node

/**
 * Deliver a canonical Data Analytics portable report while fixing the reader's
 * 100vw top-bar overflow in non-overlay-scrollbar Chromium.
 */

import { resolve } from "node:path";
import { pathToFileURL } from "node:url";

const VALUE_OPTIONS = new Set([
  "action-timeout-ms",
  "input",
  "output",
  "plugin-root",
  "ready-timeout-ms",
  "screenshot",
  "timeout-ms",
]);

function usage() {
  return [
    "Usage: node deliver_portable_artifact_scrollbar_safe.mjs",
    "  --plugin-root <data-analytics-plugin>",
    "  --input <artifact.json>",
    "  --output <report.html>",
    "  [--ready-timeout-ms <ms>]",
    "  [--action-timeout-ms <ms>]",
    "  [--timeout-ms <ms>]",
    "  [--screenshot <failure.png>]",
  ].join("\n");
}

export function parseArguments(argv) {
  const options = {};
  for (let index = 0; index < argv.length; index += 1) {
    const argument = argv[index];
    if (argument === "--help" || argument === "-h") return { help: true };
    if (!argument.startsWith("--")) {
      throw new Error(`Unexpected argument: ${argument}\n${usage()}`);
    }
    const key = argument.slice(2);
    if (!VALUE_OPTIONS.has(key)) {
      throw new Error(`Unknown argument: ${argument}\n${usage()}`);
    }
    const value = argv[index + 1];
    if (!value || value.startsWith("--")) {
      throw new Error(`Missing value for ${argument}.\n${usage()}`);
    }
    if (options[key] !== undefined) {
      throw new Error(`${argument} may only be specified once.`);
    }
    options[key] = value;
    index += 1;
  }
  for (const key of ["plugin-root", "input", "output"]) {
    if (!options[key]) throw new Error(`--${key} is required.\n${usage()}`);
  }
  return options;
}

function positiveNumber(value, label) {
  if (value === undefined) return undefined;
  const parsed = Number(value);
  if (!Number.isFinite(parsed) || parsed <= 0) {
    throw new Error(`${label} must be a positive number.`);
  }
  return parsed;
}

export function patchPortableRuntime(runtimeHtml) {
  const marker = "</head>";
  const index = runtimeHtml.toLowerCase().lastIndexOf(marker);
  if (index < 0) throw new Error("Portable reader runtime has no closing head element.");
  const style = [
    "<style data-intersubmod-scrollbar-safe-top-bar>",
    ".analytics-top-bar{",
    "width:100%!important;",
    "max-width:100%!important;",
    "margin-right:0!important;",
    "margin-left:0!important;",
    "}",
    "</style>",
  ].join("");
  return `${runtimeHtml.slice(0, index)}${style}${runtimeHtml.slice(index)}`;
}

export async function deliverWithScrollbarSafeRuntime(options) {
  const pluginRoot = resolve(options.pluginRoot);
  const buildModule = await import(pathToFileURL(resolve(
    pluginRoot,
    "skills/build-report/scripts/build_portable_artifact.mjs",
  )).href);
  const deliveryModule = await import(pathToFileURL(resolve(
    pluginRoot,
    "skills/build-report/scripts/deliver_portable_artifact.mjs",
  )).href);

  const runtimeHtml = patchPortableRuntime(
    buildModule.readPackagedReaderRuntime().html,
  );
  const build = (input, buildOptions = {}) => buildModule.buildPortableArtifact(
    input,
    { ...buildOptions, runtimeHtml },
  );

  return deliveryModule.deliverPortableArtifact(
    {
      actionTimeoutMs: positiveNumber(
        options.actionTimeoutMs,
        "--action-timeout-ms",
      ),
      inputPath: resolve(options.input),
      outputPath: resolve(options.output),
      readyTimeoutMs: positiveNumber(
        options.readyTimeoutMs,
        "--ready-timeout-ms",
      ),
      screenshotPath: options.screenshot
        ? resolve(options.screenshot)
        : undefined,
      timeoutMs: positiveNumber(options.timeoutMs, "--timeout-ms"),
    },
    { build },
  );
}

export async function runCli(argv = process.argv.slice(2)) {
  try {
    const parsed = parseArguments(argv);
    if (parsed.help) {
      process.stdout.write(`${usage()}\n`);
      return;
    }
    const result = await deliverWithScrollbarSafeRuntime({
      actionTimeoutMs: parsed["action-timeout-ms"],
      input: parsed.input,
      output: parsed.output,
      pluginRoot: parsed["plugin-root"],
      readyTimeoutMs: parsed["ready-timeout-ms"],
      screenshot: parsed.screenshot,
      timeoutMs: parsed["timeout-ms"],
    });
    process.stdout.write(`${JSON.stringify({
      ...result,
      runtimePatch: {
        id: "intersubmod-scrollbar-safe-top-bar",
        reason: "100vw includes classic vertical scrollbar width",
      },
    })}\n`);
  } catch (error) {
    const result = error?.deliveryResult ?? {
      ok: false,
      stage: "invocation",
      code: error?.code ?? "invalid_invocation",
      error: error?.message ?? String(error),
    };
    process.stderr.write(`${JSON.stringify(result)}\n`);
    process.exitCode = 1;
  }
}

const isMain = process.argv[1]
  && pathToFileURL(resolve(process.argv[1])).href === import.meta.url;
if (isMain) await runCli();
