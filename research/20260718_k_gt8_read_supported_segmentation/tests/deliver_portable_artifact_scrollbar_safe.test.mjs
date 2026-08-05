import assert from "node:assert/strict";
import test from "node:test";

import {
  parseArguments,
  patchPortableRuntime,
} from "../scripts/deliver_portable_artifact_scrollbar_safe.mjs";

test("patchPortableRuntime appends a scoped final head style", () => {
  const runtime = "<!doctype html><html><head><style>.x{color:red}</style></head><body></body></html>";
  const patched = patchPortableRuntime(runtime);
  assert.match(patched, /data-intersubmod-scrollbar-safe-top-bar/);
  assert.match(patched, /\.analytics-top-bar\{/);
  assert.match(patched, /width:100%!important/);
  assert.match(patched, /margin-left:0!important/);
  assert.ok(
    patched.indexOf("data-intersubmod-scrollbar-safe-top-bar")
      > patched.indexOf(".x{color:red}"),
  );
});

test("patchPortableRuntime rejects a runtime without head", () => {
  assert.throws(
    () => patchPortableRuntime("<html><body></body></html>"),
    /no closing head element/,
  );
});

test("parseArguments requires canonical delivery paths", () => {
  assert.deepEqual(
    parseArguments([
      "--plugin-root", "/plugin",
      "--input", "/input.json",
      "--output", "/report.html",
      "--timeout-ms", "60000",
    ]),
    {
      "plugin-root": "/plugin",
      input: "/input.json",
      output: "/report.html",
      "timeout-ms": "60000",
    },
  );
  assert.throws(() => parseArguments(["--input", "a", "--output", "b"]), /plugin-root/);
});
