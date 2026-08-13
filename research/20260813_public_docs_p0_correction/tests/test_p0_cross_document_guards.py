#!/usr/bin/env python3
"""Regression tests for registry-owned public-document semantic guards."""

from __future__ import annotations

import copy
import importlib.util
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "research/20260813_public_docs_p0_correction/scripts/validate_public_p0_claims.py"
REGISTRY = ROOT / "research/20260813_public_docs_p0_correction/scripts/p0_claim_registry.json"
SPEC = importlib.util.spec_from_file_location("validate_public_p0_claims_for_test", SCRIPT)
if SPEC is None or SPEC.loader is None:  # pragma: no cover - import infrastructure failure
    raise RuntimeError(f"cannot load validator module: {SCRIPT}")
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)
REGISTRY_DATA = json.loads(REGISTRY.read_text(encoding="utf-8"))
AUTHORITY = REGISTRY_DATA["cross_document_guards"]


class P0CrossDocumentGuardTest(unittest.TestCase):
    @staticmethod
    def authority_rule(group: str, guard_id: str) -> dict[str, object]:
        return copy.deepcopy(next(item for item in AUTHORITY[group] if item["guard_id"] == guard_id))

    def validate_text(self, group: str, guard_id: str, content: str) -> list[str]:
        selected = self.authority_rule(group, guard_id)
        selected.pop("target_globs", None)
        selected["targets"] = ["sample.md"]
        dummy_document = {
            "guard_id": "DUMMY_DOCUMENT",
            "targets": ["sample.md"],
            "required_all": ["SAFE_SENTINEL"],
            "forbidden_any": ["IMPOSSIBLE_DUMMY_FORBIDDEN"],
        }
        dummy_corpus = {
            "guard_id": "DUMMY_CORPUS",
            "targets": ["sample.md"],
            "forbidden_any": ["IMPOSSIBLE_DUMMY_FORBIDDEN"],
        }
        dummy_context = {
            "guard_id": "DUMMY_CONTEXT",
            "targets": ["sample.md"],
            "anchor_pattern": "SAFE_SENTINEL",
            "window_chars": 50,
            "required_all": ["SAFE_SENTINEL"],
        }
        contract = {
            "public_source_globs": ["sample.md"],
            "required_tracked_paths": ["sample.md"],
            "required_executable_paths": [],
            "document_rules": [selected] if group == "document_rules" else [dummy_document],
            "corpus_rules": [selected] if group == "corpus_rules" else [dummy_corpus],
            "context_rules": [selected] if group == "context_rules" else [dummy_context],
        }
        with tempfile.TemporaryDirectory(prefix="ism-p0-guard-") as temporary:
            root = Path(temporary)
            (root / "sample.md").write_text(content + "\nSAFE_SENTINEL\n", encoding="utf-8")
            subprocess.run(["git", "init", "-q"], cwd=root, check=True)
            subprocess.run(["git", "add", "--", "sample.md"], cwd=root, check=True)
            errors, _ = VALIDATOR.validate_cross_document_guards(
                {"cross_document_guards": contract}, root
            )
        return errors

    def assert_guard_rejects(self, group: str, guard_id: str, content: str) -> None:
        errors = self.validate_text(group, guard_id, content)
        self.assertTrue(any(guard_id in item for item in errors), "\n".join(errors))

    def assert_guard_accepts(self, group: str, guard_id: str, content: str) -> None:
        errors = self.validate_text(group, guard_id, content)
        self.assertEqual(errors, [], "\n".join(errors))

    def test_current_schema_rejects_59_or_193_and_accepts_199(self) -> None:
        for old_columns in (59, 193):
            self.assert_guard_rejects(
                "document_rules",
                "G_SCHEMA_CURRENT_199",
                f"Current significance_summary.csv has {old_columns} columns.",
            )
        self.assert_guard_accepts(
            "document_rules", "G_SCHEMA_CURRENT_199", "Current significance_summary.csv has 199 columns."
        )

    def test_markdown_emphasis_cannot_hide_stale_current_schema(self) -> None:
        for stale in (
            "Current `significance_summary.csv` has **193** columns.",
            "目前 `significance_summary.csv` 固定 **59** 欄。",
        ):
            self.assert_guard_rejects("document_rules", "G_SCHEMA_CURRENT_199", stale)

    def test_legacy_schema_number_requires_explicit_history_and_current_199(self) -> None:
        self.assert_guard_rejects(
            "context_rules",
            "G_LEGACY_SCHEMA_MENTION_CONTEXT",
            "The current significance_summary.csv has 59 columns.",
        )
        self.assert_guard_accepts(
            "context_rules",
            "G_LEGACY_SCHEMA_MENTION_CONTEXT",
            "Historical significance_summary.csv had 59 columns; the current schema has 199 columns.",
        )

    def test_significance_csv_alias_cannot_present_59_as_current(self) -> None:
        self.assert_guard_rejects(
            "context_rules",
            "G_LEGACY_SCHEMA_MENTION_CONTEXT",
            "Current significance.csv has 59 columns.",
        )
        self.assert_guard_accepts(
            "context_rules",
            "G_LEGACY_SCHEMA_MENTION_CONTEXT",
            "Historical significance.csv had 59 columns; current significance_summary.csv has 199 columns.",
        )

    def test_confirmed_cellular_and_linear_ancestry_ceiling_is_zero(self) -> None:
        self.assert_guard_rejects(
            "document_rules",
            "G_CONFIRMED_CLAIM_CEILING_ZERO",
            "Confirmed cellular subclones: 2. Confirmed linear ancestry relationships: 1.",
        )
        self.assert_guard_accepts(
            "document_rules",
            "G_CONFIRMED_CLAIM_CEILING_ZERO",
            "Confirmed cellular subclones: 0. Confirmed linear ancestry relationships: 0.",
        )

    def test_single_molecule_evidence_does_not_prove_cellular_ancestry(self) -> None:
        for overclaim in (
            "Same-read co-occurrence directly proves a cellular branch.",
            "Same-read co-occurrence validates an evolutionary lineage.",
            "The data confirmed child branch A.",
            "The difference must be somatic.",
            "It cannot be two parental haplotypes.",
        ):
            self.assert_guard_rejects(
                "corpus_rules", "G_NO_CELLULAR_OR_ANCESTRY_OVERCLAIM", overclaim
            )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_CELLULAR_OR_ANCESTRY_OVERCLAIM",
            "Same-read co-occurrence supports a molecular-state candidate, not a cellular lineage.",
        )

    def test_cellular_ancestry_mention_requires_negation_or_zero(self) -> None:
        self.assert_guard_rejects(
            "context_rules",
            "G_CELLULAR_ANCESTRY_MENTION_CONTEXT",
            "A confirmed cellular lineage was recovered.",
        )
        self.assert_guard_accepts(
            "context_rules",
            "G_CELLULAR_ANCESTRY_MENTION_CONTEXT",
            "Confirmed cellular lineage = 0; this is only a candidate.",
        )

    def test_tested_negative_filter_is_not_a_universal_dead_end(self) -> None:
        self.assert_guard_rejects(
            "corpus_rules",
            "G_NO_UNIVERSAL_FILTER_ABSOLUTISM",
            "Every methylation variant filter is DEAD and impossible.",
        )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_UNIVERSAL_FILTER_ABSOLUTISM",
            "The tested formulation was negative; a materially new hypothesis can be audited.",
        )

    def test_longlineage_revision_status_rejects_old_public_feature_labels(self) -> None:
        for stale_status in (
            "LongLineage public main 5daf50f has no writer.",
            "The feature branch b9aaa12 adds longlineage-tag-bam.",
        ):
            self.assert_guard_rejects(
                "corpus_rules", "G_LONGLINEAGE_REVISION_STATUS", stale_status
            )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_LONGLINEAGE_REVISION_STATUS",
            "5daf50f is a private baseline; b9aaa12 is a private public-preview candidate.",
        )

    def test_c043_uses_private_revision_scoped_status(self) -> None:
        claim = next(item for item in REGISTRY_DATA["claims"] if item["claim_id"] == "C043")
        serialized = json.dumps(claim, ensure_ascii=False)
        self.assertIn("private baseline/main snapshot", serialized)
        self.assertIn("private public-preview candidate", serialized)
        self.assertIn("NOT_READY", serialized)
        self.assertNotIn('"feature `?b9aaa12`?.*longlineage-tag-bam"', serialized)

    def test_strengthscore_current_contract_is_five_component_equal_weight(self) -> None:
        self.assert_guard_rejects(
            "document_rules",
            "G_STRENGTHSCORE_CURRENT_FIVE_COMPONENTS",
            "Current StrengthScore is a four-term weighted score.",
        )
        self.assert_guard_accepts(
            "document_rules",
            "G_STRENGTHSCORE_CURRENT_FIVE_COMPONENTS",
            "StrengthScore is the equal average of five component sub-scores; the four-term 0.30/0.25/0.25/0.20 formula is historical.",
        )

    def test_production_permanova_is_999_with_point_001_floor(self) -> None:
        self.assert_guard_rejects(
            "document_rules",
            "G_PERMANOVA_PRODUCTION_999",
            "Production uses 99 permutations and p-floor 0.01; see RegionProcessor.cpp.",
        )
        self.assert_guard_accepts(
            "document_rules",
            "G_PERMANOVA_PRODUCTION_999",
            "Production uses 999 permutations and p-floor 0.001; see RegionProcessor.cpp and LabelTest.hpp.",
        )

    def test_build_outputs_are_not_tracked_and_executable_count_is_not_fixed(self) -> None:
        for stale_statement in (
            "The checked-in build/bin binary is ready to run.",
            "The repository ships exactly 6 executables.",
        ):
            self.assert_guard_rejects(
                "corpus_rules", "G_NO_TRACKED_BINARY_OR_FIXED_EXECUTABLE_COUNT", stale_statement
            )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_TRACKED_BINARY_OR_FIXED_EXECUTABLE_COUNT",
            "Build outputs are not version-controlled; enumerate targets from the checked commit.",
        )

    def test_run_log_is_necessary_but_not_sufficient(self) -> None:
        self.assert_guard_rejects(
            "corpus_rules",
            "G_RUN_LOG_NOT_SUFFICIENT_ALONE",
            "run.log alone is sufficient to reproduce the result.",
        )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_RUN_LOG_NOT_SUFFICIENT_ALONE",
            "run.log is necessary but not sufficient; hashes and a command receipt are also required.",
        )

    def test_public_docs_reject_hardcoded_big7_copy_commands(self) -> None:
        for copy_command in (
            "cp /big7_disk/liaoyoyo2001/input.bam ./data/input.bam",
            "$ rsync /big7_disk/liaoyoyo2001/input.bam ./data/input.bam",
            "<code>cp /big7_disk/liaoyoyo2001/input.bam ./data/input.bam</code>",
        ):
            self.assert_guard_rejects(
                "corpus_rules", "G_NO_HARDCODED_BIG7_COPY_COMMAND", copy_command
            )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_HARDCODED_BIG7_COPY_COMMAND",
            "cp <DATA_ROOT>/input.bam ./data/input.bam",
        )

    def test_claim_guard_cannot_guarantee_truth_or_prevent_fabrication(self) -> None:
        for absolute_claim in (
            "The claim guard prevents fabrication.",
            "The validator guarantees scientific correctness.",
            "拒絕渲染可以防捏造。",
        ):
            self.assert_guard_rejects(
                "corpus_rules", "G_NO_ANTI_FABRICATION_GUARANTEE", absolute_claim
            )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_ANTI_FABRICATION_GUARANTEE",
            "The validator blocks omission of declared required fields; it does not guarantee scientific correctness.",
        )

    def test_inter_sub_mod_and_exact_ps_solver_must_not_be_conflated(self) -> None:
        self.assert_guard_rejects(
            "corpus_rules",
            "G_NO_ISM_EXACT_PS_CONFLATION",
            "InterSubMod produces the exact-PS topology funnel.",
        )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_ISM_EXACT_PS_CONFLATION",
            "inter_sub_mod does not produce exact-PS results; a separate exact-PS solver produces the funnel.",
        )

    def test_665_and_63506_require_their_exact_grains(self) -> None:
        self.assert_guard_rejects(
            "context_rules", "G_665_DENOMINATOR_CONTEXT", "66.52% of all mutations are isolated."
        )
        self.assert_guard_accepts(
            "context_rules",
            "G_665_DENOMINATOR_CONTEXT",
            "170,131 / 255,752 strict components are k=1 (66.52%).",
        )
        self.assert_guard_rejects(
            "context_rules", "G_63506_DENOMINATOR_CONTEXT", "The result contains 63,506 unique trees."
        )
        self.assert_guard_accepts(
            "context_rules",
            "G_63506_DENOMINATOR_CONTEXT",
            "63,506 / 71,955 = 88.26% is a rooted-unlabeled model-conditional graph-shape statistic.",
        )

    def test_unique_tree_requires_negation_or_bounded_model_context(self) -> None:
        self.assert_guard_rejects(
            "context_rules", "G_UNIQUE_TREE_SEMANTIC_CONTEXT", "A unique tree was biologically confirmed."
        )
        self.assert_guard_accepts(
            "context_rules", "G_UNIQUE_TREE_SEMANTIC_CONTEXT", "This is not a unique tree; it is a candidate."
        )

    def test_readme_test_counts_must_be_dynamic_not_copied(self) -> None:
        self.assert_guard_rejects(
            "document_rules",
            "G_README_DYNAMIC_TEST_COUNTS",
            "Current test/suite counts are dynamically generated. Audit: 270 tests / 39 suites; CTest 270/270.",
        )
        self.assert_guard_accepts(
            "document_rules",
            "G_README_DYNAMIC_TEST_COUNTS",
            "Current test/suite counts are dynamically generated from the checked commit's CTest receipt.",
        )

    def test_python_39_is_rejected_and_310_is_required(self) -> None:
        self.assert_guard_rejects(
            "document_rules", "G_PYTHON_MINIMUM_310", "Dependencies: Python 3.9."
        )
        self.assert_guard_accepts(
            "document_rules", "G_PYTHON_MINIMUM_310", "Dependencies: Python >=3.10."
        )

    def test_develop_raw_url_is_rejected(self) -> None:
        self.assert_guard_rejects(
            "corpus_rules",
            "G_NO_DEVELOP_RAW_URL",
            "https://raw.githubusercontent.com/liaoyoyo/InterSubMod/develop/docs/images/figure.png",
        )
        self.assert_guard_accepts(
            "corpus_rules", "G_NO_DEVELOP_RAW_URL", "Use the repository-relative path docs/images/figure.png."
        )

    def test_permanova_does_not_establish_cluster_truth(self) -> None:
        for overclaim in (
            "PERMANOVA proves that the clusters are real.",
            "Significant PERMANOVA establishes that the clusters are real.",
            "PERMANOVA 顯著，所以真有 HP 分群。",
            "有沒有結構要看 PERMANOVA。",
        ):
            self.assert_guard_rejects(
                "corpus_rules", "G_NO_PERMANOVA_TRUTH_OR_CLUSTER_EXISTENCE", overclaim
            )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_PERMANOVA_TRUTH_OR_CLUSTER_EXISTENCE",
            "PERMANOVA tests label-associated distance variation under a permutation null; it does not prove cluster truth.",
        )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_PERMANOVA_TRUTH_OR_CLUSTER_EXISTENCE",
            "PERMANOVA does not prove that clusters are real.",
        )

    def test_methylation_cannot_bridge_topology(self) -> None:
        for overclaim in (
            "The result uses a methylation-bridged topology.",
            "Methylation establishes the candidate tree.",
            "甲基化建立了候選樹。",
            "甲基橋接建立候選樹。",
        ):
            self.assert_guard_rejects(
                "corpus_rules", "G_NO_METHYLATION_TOPOLOGY_BRIDGE", overclaim
            )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_METHYLATION_TOPOLOGY_BRIDGE",
            "Methylation is association-only and cannot establish, merge, or direct an edge.",
        )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_METHYLATION_TOPOLOGY_BRIDGE",
            "This is not a methylation-bridged topology.",
        )

    def test_loh_annotation_does_not_confirm_mechanism(self) -> None:
        for overclaim in (
            "This is a SEQC2-confirmed LOH region.",
            "LOH was confirmed by the observed imbalance.",
            "LOH 獲證，只剩單一親本。",
        ):
            self.assert_guard_rejects(
                "corpus_rules", "G_NO_LOH_MECHANISM_OVERCLAIM", overclaim
            )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_LOH_MECHANISM_OVERCLAIM",
            "The SEQC2 LOH annotation and HP2 imbalance are compatible with an LOH context; allele-specific CN is not integrated.",
        )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_LOH_MECHANISM_OVERCLAIM",
            "LOH is not a confirmed mechanism.",
        )

    def test_seven_dataset_cohort_is_not_seven_biological_samples(self) -> None:
        for overclaim in (
            "7 samples",
            "7-sample cohort",
            "A cohort of seven biological samples was validated.",
            "全 7 樣本 funnel",
        ):
            self.assert_guard_rejects(
                "corpus_rules", "G_NO_SEVEN_SAMPLE_COHORT_SHORTHAND", overclaim
            )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_SEVEN_SAMPLE_COHORT_SHORTHAND",
            "7 technical datasets / 6 biological IDs",
        )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_SEVEN_SAMPLE_COHORT_SHORTHAND",
            "This is not a 7-sample cohort.",
        )

    def test_longlineage_blockers_are_not_reduced_to_parity(self) -> None:
        for overclaim in (
            "BLOCKED means the parity evidence does not exist.",
            "M1/M2/topology 核心都已實作。",
        ):
            self.assert_guard_rejects(
                "corpus_rules", "G_LONGLINEAGE_BLOCKER_SET_NOT_SIMPLIFIED", overclaim
            )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_LONGLINEAGE_BLOCKER_SET_NOT_SIMPLIFIED",
            "P3/P4/P5/P7/P8 remain blocked across parity, source-origin, license, dependency, and release-safety gates.",
        )

    def test_cis_heuristic_cannot_claim_confounds_removed(self) -> None:
        for overclaim in (
            "The normal-anchored cis-test removes copy, germline, and drift confounds.",
            "Copy confounding is eliminated by the normal-anchored control.",
        ):
            self.assert_guard_rejects(
                "corpus_rules", "G_NO_CIS_HEURISTIC_GUARANTEE", overclaim
            )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_CIS_HEURISTIC_GUARANTEE",
            "The normal-anchored cis-compatible heuristic can flag some observed confounds but cannot prove they were removed.",
        )
        self.assert_guard_accepts(
            "corpus_rules",
            "G_NO_CIS_HEURISTIC_GUARANTEE",
            "The normal-anchored control does not eliminate copy confounding.",
        )

    def test_untracked_fixture_file_is_rejected(self) -> None:
        dummy_document = {
            "guard_id": "DUMMY_DOCUMENT",
            "targets": ["sample.md"],
            "required_all": ["SAFE_SENTINEL"],
            "forbidden_any": ["IMPOSSIBLE_DUMMY_FORBIDDEN"],
        }
        dummy_corpus = {
            "guard_id": "DUMMY_CORPUS",
            "targets": ["sample.md"],
            "forbidden_any": ["IMPOSSIBLE_DUMMY_FORBIDDEN"],
        }
        dummy_context = {
            "guard_id": "DUMMY_CONTEXT",
            "targets": ["sample.md"],
            "anchor_pattern": "SAFE_SENTINEL",
            "window_chars": 50,
            "required_all": ["SAFE_SENTINEL"],
        }
        with tempfile.TemporaryDirectory(prefix="ism-p0-untracked-") as temporary:
            root = Path(temporary)
            (root / "sample.md").write_text("SAFE_SENTINEL\n", encoding="utf-8")
            (root / "fixture.dat").write_text("fixture\n", encoding="utf-8")
            subprocess.run(["git", "init", "-q"], cwd=root, check=True)
            subprocess.run(["git", "add", "--", "sample.md"], cwd=root, check=True)
            contract = {
                "public_source_globs": ["sample.md"],
                "required_tracked_paths": ["fixture.dat"],
                "required_executable_paths": [],
                "document_rules": [dummy_document],
                "corpus_rules": [dummy_corpus],
                "context_rules": [dummy_context],
            }
            errors, _ = VALIDATOR.validate_cross_document_guards(
                {"cross_document_guards": contract}, root
            )
        self.assertTrue(any("not exactly one stage-0 Git index entry" in item for item in errors), errors)

    def test_individual_claim_target_cannot_escape_repository(self) -> None:
        for escaped in ("/etc/hosts", "../outside.md"):
            registry = copy.deepcopy(REGISTRY_DATA)
            claim = next(item for item in registry["claims"] if item.get("checks"))
            claim["checks"][0] = {
                "target": escaped,
                "required_all": ["localhost"],
                "forbidden_any": ["IMPOSSIBLE_ESCAPE_SENTINEL"],
            }
            errors, summary = VALIDATOR.validate_registry(registry, ROOT)
            self.assertEqual(summary["verdict"], "FAIL")
            self.assertTrue(
                any("repo-relative without '..'" in item for item in errors), errors
            )


if __name__ == "__main__":
    unittest.main()
