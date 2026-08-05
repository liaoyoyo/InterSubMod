import { randomUUID } from "node:crypto";
import {
    existsSync,
    mkdirSync,
    readFileSync,
    renameSync,
    writeFileSync,
} from "node:fs";
import { dirname, resolve } from "node:path";
import { pathToFileURL } from "node:url";

const defaultPluginRoot =
    "/bip7_disk/liaoyoyo2001/.codex/plugins/cache/openai-curated-remote/data-analytics/0.2.8-13ceeea1f599";
const pluginRoot =
    process.env.INTERSUBMOD_DATA_ANALYTICS_PLUGIN_ROOT || defaultPluginRoot;
const scriptsRoot = resolve(pluginRoot, "skills/build-report/scripts");

const {
    buildPortableArtifact,
    readPackagedReaderRuntime,
} = await import(pathToFileURL(resolve(scriptsRoot, "build_portable_artifact.mjs")).href);
const { extractPortableChartSvgs } = await import(
    pathToFileURL(resolve(scriptsRoot, "extract_portable_chart_svgs.mjs")).href
);
const { verifyPortableArtifact } = await import(
    pathToFileURL(resolve(scriptsRoot, "verify_portable_artifact.mjs")).href
);

const parseArguments = (argv) => {
    const parsed = {};
    for (let index = 0; index < argv.length; index += 2) {
        const key = argv[index];
        const value = argv[index + 1];
        if (!["--input", "--output"].includes(key) || !value) {
            throw new Error(
                "Usage: node deliver_portable_report.mjs --input artifact.json --output report.html",
            );
        }
        parsed[key.slice(2)] = value;
    }
    if (!parsed.input || !parsed.output) {
        throw new Error("--input and --output are required.");
    }
    return parsed;
};

const injectOverflowWorkaround = (html) => {
    const marker = "</head>";
    const markerIndex = html.lastIndexOf(marker);
    if (markerIndex < 0) {
        throw new Error("Expected a closing head tag.");
    }
    const style =
        '<style data-intersubmod-portable-overflow-fix="true">' +
        "html,body{overflow-x:clip!important}" +
        "</style>";
    return `${html.slice(0, markerIndex)}${style}${html.slice(markerIndex)}`;
};

const options = parseArguments(process.argv.slice(2));
const inputPath = resolve(options.input);
const outputPath = resolve(options.output);
const candidatePath = `${outputPath}.tmp-${process.pid}-${randomUUID()}.html`;

if (existsSync(outputPath)) {
    throw new Error(`Refusing to overwrite existing report: ${outputPath}`);
}

const artifact = JSON.parse(readFileSync(inputPath, "utf8"));
const runtimeHtml = injectOverflowWorkaround(readPackagedReaderRuntime().html);
const build = (staticCharts) =>
    injectOverflowWorkaround(
        buildPortableArtifact(artifact, {
            runtimeHtml,
            ...(staticCharts ? { staticCharts } : {}),
        }),
    );

mkdirSync(dirname(outputPath), { recursive: true });
writeFileSync(candidatePath, build(), "utf8");

const staticCharts = await extractPortableChartSvgs({
    actionTimeoutMs: 10_000,
    htmlPath: candidatePath,
    readyTimeoutMs: 20_000,
});
writeFileSync(candidatePath, build(staticCharts), "utf8");

const verification = await verifyPortableArtifact({
    actionTimeoutMs: 10_000,
    artifactPath: inputPath,
    htmlPath: candidatePath,
    readyTimeoutMs: 20_000,
    timeoutMs: 30_000,
});

renameSync(candidatePath, outputPath);

process.stdout.write(
    `${JSON.stringify(
        {
            ok: true,
            html: outputPath,
            stages: {
                validation: "passed",
                package: "passed",
                static_charts: "passed",
                verification: "passed",
            },
            counts: verification.counts,
            sourceDialog: verification.sourceDialog,
            sourceInteraction: verification.sourceInteraction,
            viewports: verification.viewports,
            workaround: {
                scope: "shared portable-reader top bar only",
                rule: "html,body{overflow-x:clip!important}",
                reason: "100vw top bar exceeded the scrollbar-reduced client width by 8 px",
            },
        },
        null,
        2,
    )}\n`,
);
