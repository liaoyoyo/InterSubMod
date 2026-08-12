#!/usr/bin/env python3
"""drilldown 的契約測試。只用 stdlib，不需要真實資料就能跑大部分。

跑法：python3 -m pytest tests/test_drilldown_contract.py -v
      （或直接 python3 tests/test_drilldown_contract.py）
"""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
for p in ("scripts/drilldown", "scripts/drilldown/sources",
          "scripts/drilldown/panels", "scripts/drilldown/emit", "scripts"):
    sys.path.insert(0, str(ROOT / p))

CLI = ROOT / "scripts" / "build_drilldown_dashboard.py"


class TestNoHardcodedSitePaths(unittest.TestCase):
    """站點絕對路徑寫死 = 別人 clone 下來跑不動。這條曾經違反過。"""

    def test_no_site_paths_in_source(self):
        bad = []
        for f in (ROOT / "scripts" / "drilldown").rglob("*.py"):
            txt = f.read_text(encoding="utf-8")
            for line in txt.splitlines():
                s = line.strip()
                if s.startswith("#") or s.startswith("*"):
                    continue        # 註解裡舉例可以
                for pat in ("/big7_disk/", "/bip7_disk/", "/big8_disk/"):
                    if pat in s and "example" not in s.lower():
                        bad.append(f"{f.relative_to(ROOT)}: {s[:90]}")
        self.assertEqual(bad, [], "原始碼含寫死的站點路徑：\n" + "\n".join(bad))

    def test_cli_has_no_site_defaults(self):
        txt = CLI.read_text(encoding="utf-8")
        for line in txt.splitlines():
            s = line.strip()
            if s.startswith("#"):
                continue
            if "default=" in s:
                for pat in ("/big7_disk/", "/bip7_disk/", "/big8_disk/"):
                    self.assertNotIn(pat, s, f"CLI 預設值寫死站點路徑：{s[:90]}")


class TestRefuseOnMissingCore(unittest.TestCase):
    """硬核心缺席必須 exit 3 且訊息帶絕對路徑 —— 不可靜默產半成品。"""

    def test_exit_code_and_message(self):
        r = subprocess.run(
            [sys.executable, str(CLI), "--sample", "X",
             "--topology", "/nonexistent/path.jsonl", "--probe-only"],
            capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 3, f"期望 exit 3，實得 {r.returncode}")
        self.assertIn("/nonexistent/path.jsonl", r.stderr + r.stdout)

    def test_missing_config_gives_actionable_error(self):
        """沒設定站點路徑時，錯誤訊息要講清楚三種提供方式。"""
        env = {"PATH": "/usr/bin:/bin", "HOME": "/tmp"}
        with tempfile.TemporaryDirectory() as d:
            r = subprocess.run(
                [sys.executable, str(CLI), "--sample", "X", "--probe-only"],
                capture_output=True, text=True, timeout=120, cwd=d, env=env)
            out = r.stderr + r.stdout
            if "找不到" in out:
                for kw in ("CLI", "環境變數", "設定檔"):
                    self.assertIn(kw, out, "錯誤訊息缺少提供方式：" + kw)


class TestSelfcheckSemantics(unittest.TestCase):
    """SKIP 必須是獨立第三態 —— 算成 PASS 會讓人以為驗過了。"""

    def test_three_states(self):
        import selfcheck
        self.assertEqual({selfcheck.PASS, selfcheck.FAIL, selfcheck.SKIP},
                         set(selfcheck.LABEL))
        self.assertNotEqual(selfcheck.SKIP, selfcheck.PASS)

    def test_skip_not_counted_as_pass(self):
        import selfcheck
        c = selfcheck._chk("Cx", "t", selfcheck.SKIP, "d", need="nope")
        self.assertEqual(c["status"], selfcheck.SKIP)
        self.assertEqual(c["statusLabel"], "無法檢查")

    def test_c12_names_original_ssnv_denominator(self):
        """C12 denominator is MLHP original-sSNV blocks, not topology active k."""
        import selfcheck

        class Cap:
            usable = True
            linkage = None

            def __init__(self, payload):
                self.payload = payload

        topology = Cap({
            "regions": [],
            "l1": {"n": 0, "chrom": [], "csrVal": []},
        })
        strict_edges = Cap({
            "broken": 0, "checked": 2, "broken_ids": [],
            "n_pass": 3, "n_edge": 4, "broken_active": 1,
        })
        result = selfcheck.run({"topology": topology, "strict_edges": strict_edges})
        c12 = next(c for c in result["checks"] if c["id"] == "C12")
        self.assertIn("MLHP original-sSNV block", c12["title"])
        self.assertIn("不是 topology active k", c12["detail"])


class TestUiHonestyAndMobileContract(unittest.TestCase):
    """UI may not overclaim validation or expose controls that do nothing."""

    @staticmethod
    def _render(summary):
        import shell
        return shell.render(
            {"sample": "S", "l1": {"n": 1}}, [],
            {"regions": 1, "cap_summary": "test", "built_at": "now"},
            {"selfcheck": {"summary": summary}},
        )

    def test_fail_is_blocking_in_sticky_badge_and_heading(self):
        html = self._render({"pass": 10, "fail": 1, "skip": 1, "total": 12})
        self.assertIn("自檢 BLOCKED · 1 FAIL · 1 SKIP", html)
        self.assertIn("自檢未通過：1 項不成立", html)
        self.assertNotIn("確認整體沒有問題", html)

    def test_skip_is_incomplete_not_pass(self):
        html = self._render({"pass": 11, "fail": 0, "skip": 1, "total": 12})
        self.assertIn("自檢 INCOMPLETE · 1 SKIP", html)
        self.assertIn("SKIP 不是 PASS", html)

    def test_noscript_fallback_is_visible_in_source(self):
        html = self._render({"pass": 12, "fail": 0, "skip": 0, "total": 12})
        self.assertIn("<noscript>", html)
        self.assertIn("SELFCHECK.md", html)

    def test_descriptive_matrix_has_no_universal_click_claim(self):
        src = (ROOT / "scripts/drilldown/assets/analysis.js").read_text(encoding="utf-8")
        self.assertIn("var interactive = !!spec.filterMap", src)
        self.assertIn("目前矩陣是描述性摘要", src)
        self.assertNotIn("點任一格", src)
        self.assertIn("role='button' tabindex='0'", src,
                      "real filterMap interactions must remain keyboard accessible")

    def test_no_unproven_cn_altered_percentage(self):
        src = (ROOT / "scripts/drilldown/assets/locus.js").read_text(encoding="utf-8")
        self.assertNotIn("94%", src)
        self.assertIn("sample-matched CN provenance", src)

    def test_sticky_offset_and_panel_overflow_contract(self):
        css = (ROOT / "scripts/drilldown/assets/dashboard.css").read_text(encoding="utf-8")
        shell_src = (ROOT / "scripts/drilldown/emit/shell.py").read_text(encoding="utf-8")
        self.assertIn("top: var(--authority-height)", css)
        self.assertIn("ResizeObserver", shell_src)
        self.assertIn('data-mode="zoom"', css)
        self.assertIn("contain: inline-size", css)
        self.assertIn('data-mode="fit"] .panel-fig { width: 100% !important', css)
        self.assertIn(".authority .authority-secondary { display: none; }", css)
        self.assertIn(".levelbar .crumb { display: none; }", css)
        self.assertIn("@media (max-width: 420px)", css)
        self.assertIn('class="authority-primary"', shell_src)

    def test_analysis_mounts_use_real_payload_or_honest_unavailable(self):
        shell_src = (ROOT / "scripts/drilldown/emit/shell.py").read_text(encoding="utf-8")
        analysis_src = (ROOT / "scripts/drilldown/assets/analysis.js").read_text(encoding="utf-8")
        self.assertNotIn("待 P3 實作", shell_src)
        self.assertIn('id=\'cooccurContent\'', shell_src)
        self.assertIn('id=\'selfcheckContent\'', shell_src)
        self.assertIn("空白不代表交集為 0", analysis_src)
        self.assertIn("SELFCHECK.md", analysis_src)

    def test_observation_scope_precedes_filters_and_exposes_claim_ceiling(self):
        import shell
        html = shell.render(
            {
                "sample": "S", "l1": {"n": 10},
                "observationScope": {
                    "available": False,
                    "claim_ceiling": "OBSERVATION_ONLY_NO_TRUTH_SET",
                    "universe": {"available": True, "num": 10, "den": 10,
                                 "label": "topology 觀察母體"},
                    "summary_join": {"available": False, "num": None, "den": 10},
                    "axis_testable": {"available": False, "num": None, "den": 10},
                    "axis_significant_raw": {"available": False, "num": None, "den": None},
                    "detail_available": {"available": False, "num": None, "den": 10},
                    "reason": "missing ISM",
                    "multiplicity": "raw p only",
                    "missingness": "not random",
                    "legacy_crosswalk": "legacy class is not current axis",
                    "axis_definitions": [],
                },
            }, [], {"regions": 1, "cap_summary": "test", "built_at": "now"},
            {"cooccur": [], "selfcheck": {"summary": {"pass": 0, "fail": 0,
                                                         "skip": 0, "total": 0}}},
        )
        self.assertLess(html.index("interpretation-rail"), html.index('id="scope"'))
        self.assertIn("Claim ceiling：觀察性下鑽", html)
        self.assertIn("UNAVAILABLE", html)
        self.assertIn("不可解讀為 0", html)
        self.assertIn("data-view-target='selfcheck'", html)

    def test_axis_missingness_and_legacy_classes_are_not_collapsed(self):
        import payload
        l1 = {"n": 2, "chroms": ["chr1"], "chrom": [0, 0], "dpos": [100, 100]}

        class Cap:
            usable = True
            reason = ""
            payload = {
                "rows": {
                    ("chr1", 100): {
                        "hp_p": None, "hp_valid": None, "hp_n": None,
                        "allele_p": None, "allele_valid": None, "allele_n": None,
                        "lineage_p": None, "lineage_valid": False, "lineage_n": None,
                    }
                },
                "axes": [{"id": "hp"}, {"id": "allele"}, {"id": "lineage"}],
                "hit_dir": {0},
            }

        cap = Cap()
        payload.fill_axis_codes(l1, cap)
        self.assertEqual(l1["axis"], [payload.AXIS_UNTESTED, 0],
                         "summary present-but-untested must differ from no summary row")
        scope = payload.observation_scope(l1, cap)
        self.assertEqual(scope["summary_join"]["num"], 1)
        self.assertEqual(scope["axis_testable"]["num"], 0)
        self.assertIn("不把舊 standalone", scope["legacy_crosswalk"])
        self.assertIn("raw p", scope["multiplicity"])

    def test_no_igv_dead_branch(self):
        cli_src = CLI.read_text(encoding="utf-8")
        self.assertNotIn("if False else", cli_src)

    def test_lightweight_igv_skip_is_explicit_in_receipt_and_ui(self):
        cli_src = CLI.read_text(encoding="utf-8")
        shell_src = (ROOT / "scripts/drilldown/emit/shell.py").read_text(encoding="utf-8")
        self.assertIn('"--igv-mode"', cli_src)
        self.assertIn('"asset_policy": boot["assetPolicy"]', cli_src)
        self.assertIn("此 build 為輕量觀察版", shell_src)
        self.assertIn("缺圖不可解讀為沒有 read-level 訊號", shell_src)

    def test_analysis_and_panel_runtime_guards(self):
        analysis_src = (ROOT / "scripts/drilldown/assets/analysis.js").read_text(encoding="utf-8")
        panel_src = (ROOT / "scripts/drilldown/assets/panel.js").read_text(encoding="utf-8")
        self.assertIn('inner.querySelectorAll("td[data-x]")', analysis_src)
        self.assertNotIn('host.querySelectorAll("td[data-x]")', analysis_src)
        self.assertIn("if (!boot.panels)", panel_src)
        self.assertIn("輸出資產政策", panel_src)


class TestNCHSSuppression(unittest.TestCase):
    """n < 30 不報百分比（NCHS 標準）。"""

    def test_small_n_suppressed(self):
        import cooccurrence as co
        m = co._matrix([("a", "x")] * 5 + [("b", "y")] * 5, ["a", "b"], ["x", "y"])
        for row in m["cells"]:
            for cell in row:
                if cell["n"] < co.MIN_N_FOR_PCT:
                    self.assertTrue(cell["small"])
                    self.assertIsNone(cell["rowPct"])


class TestPaletteIsSingleSource(unittest.TestCase):
    """色盤只能有一份 SoT —— 本專案已有三套配色，不可再加第四套。"""

    def test_composite_imports_std(self):
        src = (ROOT / "scripts/drilldown/panels/composite.py").read_text(encoding="utf-8")
        self.assertIn("import ism_heatmap_std", src)

    def test_palette_export_shape(self):
        import composite
        p = composite.palette_export()
        self.assertIn("tracks", p)
        self.assertIn("scales", p)
        ids = {t["id"] for t in p["tracks"]}
        self.assertTrue({"HP", "ALT", "cluster"} <= ids)
        for sc in p["scales"]:
            self.assertTrue(sc["stops"], "色階不可為空")


class TestPngEncoder(unittest.TestCase):
    """零依賴 PNG 編碼器要產出合法檔頭。"""

    def test_signature_and_size(self):
        import pngenc
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "t.png"
            cv = pngenc.Canvas(4, 3, "#000000")
            cv.rect(1, 1, 2, 1, (255, 0, 0))
            n = cv.save(f)
            self.assertGreater(n, 0)
            head = f.read_bytes()[:8]
            self.assertEqual(head, b"\x89PNG\r\n\x1a\n")


class TestAnnotationDropIn(unittest.TestCase):
    """SAVANA 原生格式要能直接吃，門檻與 savana_to_smcnbed.py 一致。"""

    def test_savana_thresholds_match_converter(self):
        import annotation as ann
        conv = (ROOT / "scripts/analysis/savana_to_smcnbed.py").read_text(encoding="utf-8")
        self.assertIn("gain_cn", conv)
        self.assertEqual(ann.GAIN_CN, 3.0)
        self.assertEqual(ann.LOSS_CN, 1.0)
        self.assertEqual(ann.LOH_MINOR, 0.5)

    def test_segment_classification(self):
        import annotation as ann
        with tempfile.TemporaryDirectory() as d:
            f = Path(d) / "seg.tsv"
            f.write_text("chromosome\tstart\tend\tcopyNumber\tminorAlleleCopyNumber\n"
                         "chr1\t100\t200\t4.0\t1.5\n"      # gain
                         "chr1\t200\t300\t1.0\t0.9\n"      # loss
                         "chr1\t300\t400\t2.0\t0.1\n",     # loh
                         encoding="utf-8")
            seg, _cols = ann._load_segments(f)
            labs = seg["chr1"][2]
            self.assertEqual(labs, ["gain", "loss", "loh"])


class TestOutputContract(unittest.TestCase):
    """已產出的儀表板（若存在）必須符合輸出契約。"""

    OUT = Path("/bip7_disk/liaoyoyo2001/drilldown_out/HCC1395_v1")

    def setUp(self):
        if not (self.OUT / "receipt.json").is_file():
            self.skipTest("本機沒有已產出的儀表板，跳過")

    def test_receipt_has_input_hashes(self):
        r = json.loads((self.OUT / "receipt.json").read_text(encoding="utf-8"))
        self.assertTrue(r.get("inputs"), "receipt 缺 inputs")
        for i in r["inputs"]:
            if i.get("exists") and i.get("size_bytes"):
                self.assertTrue(i.get("sha256"), f"{i['path']} 缺 sha256")

    def test_selfcheck_present(self):
        self.assertTrue((self.OUT / "SELFCHECK.md").is_file())

    def test_shards_exist_for_each_chrom(self):
        for i in range(1, 23):
            self.assertTrue((self.OUT / "data" / f"L2.chr{i}.js").is_file(),
                            f"缺 chr{i} 分片")


if __name__ == "__main__":
    unittest.main(verbosity=2)
