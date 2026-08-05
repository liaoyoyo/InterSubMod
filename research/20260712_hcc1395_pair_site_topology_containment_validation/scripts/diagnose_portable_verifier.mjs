#!/usr/bin/env node
// Print the raw portable-verifier result, including overflow dimensions.

import { randomUUID } from "node:crypto";
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";
import { pathToFileURL } from "node:url";

const [pluginRootArg, htmlArg, artifactArg] = process.argv.slice(2);
if (!pluginRootArg || !htmlArg || !artifactArg) {
  throw new Error("usage: diagnose_portable_verifier.mjs <plugin-root> <html> <artifact>");
}

const pluginRoot = resolve(pluginRootArg);
const browserCli = await import(
  pathToFileURL(join(pluginRoot, "skills/build-report/scripts/portable_browser_cli.mjs"))
);
const browserHelpers = await import(
  pathToFileURL(join(pluginRoot, "skills/build-report/scripts/portable_browser_helpers.mjs"))
);
const verifier = await import(
  pathToFileURL(join(pluginRoot, "skills/build-report/scripts/verify_portable_artifact.mjs"))
);

const html = readFileSync(resolve(htmlArg), "utf8");
const artifact = JSON.parse(readFileSync(resolve(artifactArg), "utf8"));
const viewport = { name: "desktop", width: 1440, height: 1000 };
const channel = randomUUID();
const instrumented = browserCli.injectPortableVerifierProbe(html, {
  actionTimeoutMs: 10_000,
  channel,
  checkSource: false,
  expectedCounts: verifier.expectedPortableCounts(artifact),
  readerRoot: "#data-analytics-portable-reader-root",
  readyTimeoutMs: 20_000,
  title: artifact.manifest.title,
  viewport,
});

const temporaryDirectory = mkdtempSync(join(tmpdir(), "hcc1395-portable-diagnostic-"));
try {
  const harnessPath = join(temporaryDirectory, "harness.html");
  writeFileSync(
    harnessPath,
    browserCli.buildPortableVerifierHarness({
      channel,
      frames: [{ html: instrumented, viewport }],
      timeoutMs: 31_000,
    }),
    "utf8",
  );
  const result = await browserCli.spawnChromiumDump({
    arguments: browserCli.chromiumDumpArguments({
      height: viewport.height,
      profilePath: join(temporaryDirectory, "profile"),
      url: pathToFileURL(harnessPath).href,
      virtualTimeBudgetMs: 31_250,
      width: viewport.width,
    }),
    executablePath: browserHelpers.resolveChromiumExecutable(),
    timeoutMs: 60_000,
  });
  const verifierResult = browserCli.parsePortableVerifierDump(result.stdout);

  const diagnosticScript = `<script>(()=>{
    const finish=()=>{
      if(document.documentElement.dataset.dataAnalyticsPortableReader!=="ready"){
        setTimeout(finish,25); return;
      }
      setTimeout(()=>{
        const vw=document.documentElement.clientWidth;
        const rows=[...document.querySelectorAll("body *")].map(el=>{
          const r=el.getBoundingClientRect();
          const s=getComputedStyle(el);
          return {tag:el.tagName.toLowerCase(),id:el.id||"",cls:el.className?.baseVal||el.className||"",text:(el.textContent||"").replace(/\\s+/g," ").trim().slice(0,100),left:r.left,right:r.right,width:r.width,clientWidth:el.clientWidth,scrollWidth:el.scrollWidth,overflowX:s.overflowX};
        }).filter(x=>x.width>0&&((x.right>vw+1&&x.right<vw+120)||x.left<-1)).sort((a,b)=>(b.right-vw)-(a.right-vw)).slice(0,120);
        parent.postMessage({kind:"layout-diagnostic",payload:{clientWidth:vw,scrollWidth:document.documentElement.scrollWidth,bodyScrollWidth:document.body.scrollWidth,rows}},"*");
      },100);
    }; finish();
  })();</script>`;
  const headMatch = /<head(?:\\s[^>]*)?>/i.exec(html);
  const diagnosticHtml = `${html.slice(0, headMatch.index + headMatch[0].length)}${diagnosticScript}${html.slice(headMatch.index + headMatch[0].length)}`;
  const diagnosticHtmlBase64 = Buffer.from(diagnosticHtml, "utf8").toString("base64");
  const diagnosticHarness = `<!doctype html><html><head><meta charset="utf-8"></head><body style="margin:0"><iframe id="f" sandbox="allow-scripts" style="border:0;display:block;width:1440px;height:1000px" width="1440" height="1000"></iframe><script>
    const encode=v=>btoa(unescape(encodeURIComponent(JSON.stringify(v))));
    addEventListener("message",e=>{if(e.data?.kind!=="layout-diagnostic")return;const m=document.createElement("meta");m.id="layout-diagnostic";m.dataset.result=encode(e.data.payload);document.head.append(m);document.querySelector("#f").removeAttribute("srcdoc");});
    document.querySelector("#f").srcdoc=new TextDecoder().decode(Uint8Array.from(atob(${JSON.stringify(diagnosticHtmlBase64)}),c=>c.charCodeAt(0)));
  </script></body></html>`;
  const diagnosticHarnessPath = join(temporaryDirectory, "layout-harness.html");
  writeFileSync(diagnosticHarnessPath, diagnosticHarness, "utf8");
  const diagnosticDump = await browserCli.spawnChromiumDump({
    arguments: browserCli.chromiumDumpArguments({
      height: viewport.height,
      profilePath: join(temporaryDirectory, "diagnostic-profile"),
      url: pathToFileURL(diagnosticHarnessPath).href,
      virtualTimeBudgetMs: 15_000,
      width: viewport.width,
    }),
    executablePath: browserHelpers.resolveChromiumExecutable(),
    timeoutMs: 45_000,
  });
  const match = diagnosticDump.stdout.match(/<meta[^>]+id="layout-diagnostic"[^>]+data-result="([^"]+)"/i);
  const layoutDiagnostic = match
    ? JSON.parse(decodeURIComponent(escape(Buffer.from(match[1], "base64").toString("binary"))))
    : {
        error: "layout diagnostic marker missing",
        dumpHasMarkerText: diagnosticDump.stdout.includes("layout-diagnostic"),
        dumpTail: diagnosticDump.stdout.slice(-2_000),
        stderrTail: diagnosticDump.stderr.slice(-2_000),
      };
  console.log(JSON.stringify({ verifierResult, layoutDiagnostic }, null, 2));
} finally {
  rmSync(temporaryDirectory, { force: true, recursive: true });
}
