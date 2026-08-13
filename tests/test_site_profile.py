import importlib.util
import hashlib
import json
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "scripts" / "site" / "site_profile.py"
SPEC = importlib.util.spec_from_file_location("site_profile", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)
TOOL_FIXTURE_BYTES = b"#!/bin/sh\nexit 0\n"


def valid_profile(root: Path) -> dict:
    file_spec = lambda name: {
        "path": str(root / name),
        "checksum_policy": "EXISTENCE_ONLY",
    }
    inputs = {
        role: file_spec(role)
        for role in (
            "tumor_bam", "tumor_bam_index", "normal_bam", "normal_bam_index",
            "somatic_vcf", "somatic_vcf_index", "somatic_vcf_pileup",
            "somatic_vcf_pileup_index", "somatic_vcf_indel", "somatic_vcf_indel_index",
            "to_somatic_vcf", "to_somatic_vcf_index", "to_somatic_vcf_pileup",
            "to_somatic_vcf_pileup_index", "to_somatic_vcf_indel",
            "to_somatic_vcf_indel_index", "truth_vcf", "truth_vcf_index", "truth_bed",
        )
    }
    inputs["germline_phased_dir"] = {"path": str(root / "germline")}
    return {
        "schema_version": "1.1.0",
        "site_id": "fixture",
        "expected_hostname": MODULE.platform.node(),
        "project_root": str(ROOT),
        "data_roots": {"primary": str(root / "data"), "output": str(root / "output")},
        "reference": {
            "genome_build": "GRCh38",
            "contig_contract": {
                "naming": "CHR_PREFIXED",
                "scope": "AUTOSOMES_1_22",
                "required_contigs": [f"chr{index}" for index in range(1, 23)],
            },
            "fasta": file_spec("reference.fa"),
            "indexes": [file_spec("reference.fa.fai")],
        },
        "tools": {
            name: {
                "path": str(root / name),
                "checksum_policy": "EXPECTED_SHA256",
                "sha256": hashlib.sha256(TOOL_FIXTURE_BYTES).hexdigest(),
                "required": True,
            }
            for name in ("samtools", "bcftools", "longphase_s", "intersubmod")
        },
        "datasets": {
            "operator-alias": {
                "biological_id": "S1-bio",
                "technical_dataset_id": "S1",
                "genome_build": "GRCh38",
                "platform_label": "fixture-platform",
                "truth_set_label": "fixture-truth",
                "truth_total": 7,
                "inputs": inputs,
            }
        },
    }


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def materialize_profile(profile: dict) -> None:
    for path in profile["data_roots"].values():
        Path(path).mkdir(parents=True, exist_ok=True)
    for spec in profile["tools"].values():
        path = Path(spec["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(TOOL_FIXTURE_BYTES)
        path.chmod(0o755)
    for spec in (profile["reference"]["fasta"], *profile["reference"]["indexes"]):
        path = Path(spec["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".fai":
            path.write_text(
                "".join(
                    f"{contig}\t100\t0\t100\t101\n"
                    for contig in profile["reference"]["contig_contract"]["required_contigs"]
                ),
                encoding="utf-8",
            )
        else:
            path.write_text("fixture\n", encoding="utf-8")
    for dataset in profile["datasets"].values():
        for role, spec in dataset["inputs"].items():
            path = Path(spec["path"])
            if role == "germline_phased_dir":
                path.mkdir(parents=True, exist_ok=True)
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("fixture\n", encoding="utf-8")


class SiteProfileTest(unittest.TestCase):
    def test_example_has_supported_schema(self):
        profile = json.loads((ROOT / "config" / "site-profile.example.json").read_text())
        self.assertEqual(MODULE.validate_profile(profile), [])

    def test_shell_assignments_join_record_field_and_are_quoted(self):
        with tempfile.TemporaryDirectory(prefix="profile with spaces ") as directory:
            profile = valid_profile(Path(directory))
            output = MODULE.shell_assignments(profile, "S1")
            self.assertIn("TUMOR_BAM=", output)
            self.assertIn("GERMLINE_PHASED_DIR=", output)
            self.assertIn("TRUTH_TOTAL=7", output)
            self.assertIn("'", output)

            del profile["datasets"]["operator-alias"]["inputs"]["somatic_vcf_pileup_index"]
            output = MODULE.shell_assignments(profile, "S1")
            assignments = dict(line.split("=", 1) for line in output.splitlines())
            self.assertEqual(assignments["SOMATIC_VCF_PILEUP_INDEX"], "''")
            isolated = subprocess.run(
                [
                    "bash", "-c",
                    'SOMATIC_VCF_PILEUP_INDEX=INHERITED_POISON; eval "$1"; '
                    'test -z "$SOMATIC_VCF_PILEUP_INDEX"',
                    "profile-shell", output,
                ],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertEqual(isolated.returncode, 0, isolated.stderr)

    def test_reference_contract_and_dataset_build_must_match(self):
        with tempfile.TemporaryDirectory() as directory:
            profile = valid_profile(Path(directory))
            profile["datasets"]["operator-alias"]["genome_build"] = "GRCh37"
            self.assertTrue(any("must equal reference.genome_build" in error for error in MODULE.validate_profile(profile)))

            profile = valid_profile(Path(directory))
            profile["reference"]["contig_contract"]["required_contigs"][-1] = "chrX"
            self.assertTrue(any("exactly list autosomes" in error for error in MODULE.validate_profile(profile)))

    def test_profile_completeness_and_unknown_role_fail_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            profile = valid_profile(Path(directory))
            del profile["datasets"]["operator-alias"]["inputs"]["germline_phased_dir"]
            errors = MODULE.validate_profile(profile)
            self.assertTrue(any("germline_phased_dir" in error for error in errors))

            profile = valid_profile(Path(directory))
            profile["datasets"]["operator-alias"]["inputs"]["mystery_role"] = {
                "path": "/tmp/mystery", "checksum_policy": "EXISTENCE_ONLY"
            }
            errors = MODULE.validate_profile(profile)
            self.assertTrue(any("unknown roles: mystery_role" in error for error in errors))

            profile = valid_profile(Path(directory))
            profile["tools"]["bcftool_typo"] = profile["tools"]["bcftools"].copy()
            errors = MODULE.validate_profile(profile)
            self.assertTrue(any("tools has unknown roles: bcftool_typo" in error for error in errors))

            profile = valid_profile(Path(directory))
            del profile["tools"]["samtools"]["required"]
            errors = MODULE.validate_profile(profile)
            self.assertTrue(any("tools.samtools.required must be explicitly true" in error for error in errors))

            profile = valid_profile(Path(directory))
            profile["tools"]["samtools"] = {
                "path": str(Path(directory) / "samtools"),
                "checksum_policy": "EXISTENCE_ONLY",
                "required": True,
            }
            errors = MODULE.validate_profile(profile)
            self.assertTrue(any("cannot be EXISTENCE_ONLY" in error for error in errors))

    def test_bam_and_random_access_vcf_indexes_are_required(self):
        with tempfile.TemporaryDirectory() as directory:
            profile = valid_profile(Path(directory))
            del profile["datasets"]["operator-alias"]["inputs"]["tumor_bam_index"]
            self.assertTrue(any("tumor_bam_index" in error for error in MODULE.validate_profile(profile)))

            profile = valid_profile(Path(directory))
            inputs = profile["datasets"]["operator-alias"]["inputs"]
            inputs["somatic_vcf"]["path"] = str(Path(directory) / "somatic.vcf.gz")
            del inputs["somatic_vcf_index"]
            errors = MODULE.validate_profile(profile)
            self.assertTrue(any("required for random-access somatic_vcf" in error for error in errors))

    def test_size_only_requires_positive_expected_size_and_matches_exactly(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.bam"
            path.write_bytes(b"fixture")
            profile = valid_profile(Path(directory))
            spec = profile["datasets"]["operator-alias"]["inputs"]["tumor_bam"]
            spec["path"] = str(path)
            spec["checksum_policy"] = "SIZE_ONLY"
            self.assertTrue(any("size_bytes" in error for error in MODULE.validate_profile(profile)))
            spec["size_bytes"] = 0
            self.assertTrue(any("positive integer" in error for error in MODULE.validate_profile(profile)))
            spec["size_bytes"] = 7
            self.assertEqual(MODULE.validate_profile(profile), [])
            self.assertEqual(MODULE.inspect_file("input", spec)["status"], "SIZE_MATCH")
            spec["size_bytes"] = 8
            self.assertEqual(MODULE.inspect_file("input", spec)["status"], "SIZE_MISMATCH")

    def test_zero_byte_file_fails_even_existence_only(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "empty"
            path.touch()
            result = MODULE.inspect_file(
                "fixture", {"path": str(path), "checksum_policy": "EXISTENCE_ONLY"}
            )
            self.assertEqual(result["status"], "EMPTY")

    def test_expected_sha_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tool"
            path.write_text("fixture")
            path.chmod(0o755)
            result = MODULE.inspect_file(
                "tool:fixture",
                {"path": str(path), "checksum_policy": "EXPECTED_SHA256", "sha256": "0" * 64},
                executable=True,
            )
            self.assertEqual(result["status"], "HASH_MISMATCH")

    def test_tool_hash_is_measured_under_non_hash_policy(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tool"
            path.write_bytes(b"x")
            path.chmod(0o755)
            result = MODULE.inspect_file(
                "tool:fixture",
                {"path": str(path), "checksum_policy": "SIZE_ONLY", "size_bytes": 1},
                executable=True,
                always_hash=True,
            )
            self.assertEqual(result["status"], "SIZE_MATCH")
            self.assertEqual(result["actual_sha256"], MODULE.sha256_file(path))

    def test_machine_preflight_receipt_binds_profile_mounts_and_tool_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            profile = valid_profile(root)
            (root / "data").mkdir()
            (root / "output").mkdir()
            (root / "germline").mkdir()
            for spec in profile["reference"]["indexes"] + [profile["reference"]["fasta"]]:
                path = Path(spec["path"])
                if path.suffix == ".fai":
                    path.write_text(
                        "".join(
                            f"{contig}\t100\t0\t100\t101\n"
                            for contig in profile["reference"]["contig_contract"]["required_contigs"]
                        ),
                        encoding="utf-8",
                    )
                else:
                    path.write_bytes(b"x")
            for spec in profile["tools"].values():
                Path(spec["path"]).write_bytes(TOOL_FIXTURE_BYTES)
                Path(spec["path"]).chmod(0o755)
            for role, spec in profile["datasets"]["operator-alias"]["inputs"].items():
                if role != "germline_phased_dir":
                    Path(spec["path"]).write_bytes(b"x")
            profile_path = root / "profile.json"
            profile_path.write_text(json.dumps(profile), encoding="utf-8")

            completed = subprocess.run(
                [
                    str(SCRIPT), "preflight", "--profile", str(profile_path),
                    "--sample", "S1", "--mode", "real-preflight",
                ],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            receipt = json.loads(completed.stdout)
            self.assertTrue(receipt["pass"])
            self.assertEqual(receipt["receipt_type"], "site_preflight")
            self.assertEqual(receipt["profile"]["path"], str(profile_path.resolve()))
            self.assertEqual(receipt["profile"]["sha256"], MODULE.sha256_file(profile_path))
            self.assertEqual(receipt["selected_datasets"], ["operator-alias"])
            self.assertEqual(set(receipt["tool_hashes"]), {"samtools", "bcftools", "longphase_s", "intersubmod"})
            self.assertTrue(all(receipt["tool_hashes"].values()))
            self.assertIn("primary", receipt["mounts"])
            contigs = next(
                check for check in receipt["checks"]
                if check["label"] == "reference:contig_contract"
            )
            self.assertEqual(contigs["status"], "CONTIG_CONTRACT_MATCH")

    def test_machine_preflight_rejects_reference_fai_contig_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            profile = valid_profile(root)
            materialize_profile(profile)
            Path(profile["reference"]["indexes"][0]["path"]).write_text(
                "chr1\t100\t0\t100\t101\n", encoding="utf-8"
            )
            result = MODULE.inspect_profile(profile, include_real_data=True)
            contract = next(
                check for check in result["checks"]
                if check["label"] == "reference:contig_contract"
            )
            self.assertEqual(contract["status"], "CONTIG_CONTRACT_MISMATCH")
            self.assertIn("chr22", contract["missing_contigs"])
            self.assertFalse(result["pass"])

    def test_missing_local_checksum_locator_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.bam"
            path.write_bytes(b"fixture")
            result = MODULE.inspect_file(
                "dataset:fixture:tumor_bam",
                {
                    "path": str(path),
                    "checksum_policy": "LOCATOR_ONLY",
                    "checksum_uri": str(Path(directory) / "missing.sha256"),
                },
            )
            self.assertEqual(result["status"], "CHECKSUM_LOCATOR_MISSING")

    def test_plan_uses_only_profile_roles_and_has_zero_side_effects(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            profile = valid_profile(root)
            profile_path = root / "profile.json"
            profile_path.write_text(json.dumps(profile), encoding="utf-8")
            tagged_bam = root / "existing.tagged.bam"
            tp_vcf = root / "existing.tp.vcf"
            fp_vcf = root / "existing.fp.vcf"
            for path in (tagged_bam, tp_vcf, fp_vcf):
                path.write_text("fixture\n", encoding="utf-8")
            before = sorted(str(path.relative_to(root)) for path in root.rglob("*"))
            completed = subprocess.run(
                [
                    str(ROOT / "scripts" / "pipeline" / "run_benchmark.sh"),
                    "--site-profile", str(profile_path), "--sample", "S1",
                    "--mode", "to-pure", "--vcf-source", "pileup",
                    "--output-root", str(root / "output" / "canonical"),
                    "--skip-longphase", "--skip-cleanup", "--plan",
                    "--tagged-bam", str(tagged_bam), "--tagged-bam-sha256", sha256(tagged_bam),
                    "--tp-vcf", str(tp_vcf), "--tp-vcf-sha256", sha256(tp_vcf),
                    "--fp-vcf", str(fp_vcf), "--fp-vcf-sha256", sha256(fp_vcf),
                ],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            after = sorted(str(path.relative_to(root)) for path in root.rglob("*"))
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn(f"selected_vcf={root / 'to_somatic_vcf_pileup'}", completed.stdout)
            self.assertIn("stages=InterSubMod,filter-analysis", completed.stdout)
            self.assertIn(f"existing_tagged_bam={tagged_bam}", completed.stdout)
            self.assertIn(f"existing_tp_vcf={tp_vcf}", completed.stdout)
            self.assertIn(f"existing_fp_vcf={fp_vcf}", completed.stdout)
            self.assertIn("single_parent_parse", completed.stdout)
            self.assertNotIn("/big8_disk", completed.stdout + completed.stderr)
            self.assertEqual(after, before)

    def test_profile_skip_longphase_requires_three_sha_bound_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            profile = valid_profile(root)
            profile_path = root / "profile.json"
            profile_path.write_text(json.dumps(profile), encoding="utf-8")
            completed = subprocess.run(
                [
                    str(ROOT / "scripts" / "pipeline" / "run_benchmark.sh"),
                    "--site-profile", str(profile_path), "--sample", "S1",
                    "--skip-longphase", "--plan",
                ],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertEqual(completed.returncode, 3)
            self.assertIn("SHA-256", completed.stderr)
            self.assertNotIn("latest", completed.stdout + completed.stderr)

    def test_profile_run_writes_locked_profile_and_preflight_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            profile = valid_profile(root)
            materialize_profile(profile)
            profile_path = root / "profile.json"
            profile_path.write_text(json.dumps(profile), encoding="utf-8")
            tagged_bam = root / "existing.tagged.bam"
            tp_vcf = root / "existing.tp.vcf"
            fp_vcf = root / "existing.fp.vcf"
            for path in (tagged_bam, tp_vcf, fp_vcf):
                path.write_text("##fixture\n", encoding="utf-8")
            completed = subprocess.run(
                [
                    str(ROOT / "scripts" / "pipeline" / "run_benchmark.sh"),
                    "--site-profile", str(profile_path), "--sample", "S1",
                    "--skip-longphase", "--skip-intersubmod", "--skip-cleanup",
                    "--dry-run", "--min-free-gb", "0", "--run-tag", "lock-test",
                    "--tagged-bam", str(tagged_bam), "--tagged-bam-sha256", sha256(tagged_bam),
                    "--tp-vcf", str(tp_vcf), "--tp-vcf-sha256", sha256(tp_vcf),
                    "--fp-vcf", str(fp_vcf), "--fp-vcf-sha256", sha256(fp_vcf),
                ],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            runs = list((root / "output" / "canonical" / "S1").glob("*/*_lock-test"))
            self.assertEqual(len(runs), 1)
            provenance = runs[0] / "provenance"
            locked = provenance / "site-profile.locked.json"
            receipt = provenance / "real_data_preflight_receipt.json"
            self.assertEqual(locked.read_bytes(), profile_path.read_bytes())
            self.assertEqual(locked.stat().st_mode & 0o777, 0o444)
            self.assertTrue(json.loads(receipt.read_text())["pass"])
            command_receipt = (provenance / "command_receipt.txt").read_text()
            self.assertIn("min_free_gb=CLI_OVERRIDE:0", command_receipt)
            self.assertIn("run_tag=CLI_OVERRIDE:lock-test", command_receipt)

    def test_profile_run_rejects_existing_artifact_sha_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            profile = valid_profile(root)
            materialize_profile(profile)
            profile_path = root / "profile.json"
            profile_path.write_text(json.dumps(profile), encoding="utf-8")
            paths = [root / name for name in ("existing.tagged.bam", "existing.tp.vcf", "existing.fp.vcf")]
            for path in paths:
                path.write_text("fixture\n", encoding="utf-8")
            completed = subprocess.run(
                [
                    str(ROOT / "scripts" / "pipeline" / "run_benchmark.sh"),
                    "--site-profile", str(profile_path), "--sample", "S1",
                    "--skip-longphase", "--skip-intersubmod", "--dry-run", "--min-free-gb", "0",
                    "--tagged-bam", str(paths[0]), "--tagged-bam-sha256", "0" * 64,
                    "--tp-vcf", str(paths[1]), "--tp-vcf-sha256", sha256(paths[1]),
                    "--fp-vcf", str(paths[2]), "--fp-vcf-sha256", sha256(paths[2]),
                ],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertEqual(completed.returncode, 3)
            self.assertIn("Existing tagged BAM SHA-256 mismatch", completed.stdout + completed.stderr)

    def test_profile_mode_forces_legacy_cleanup_off(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            profile = valid_profile(root)
            profile_path = root / "profile.json"
            profile_path.write_text(json.dumps(profile), encoding="utf-8")
            completed = subprocess.run(
                [
                    str(ROOT / "scripts" / "pipeline" / "run_benchmark.sh"),
                    "--site-profile", str(profile_path), "--sample", "S1", "--plan",
                ],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            stage_line = next(line for line in completed.stdout.splitlines() if line.startswith("[PLAN] stages="))
            self.assertNotIn("cleanup", stage_line)
            self.assertIn("forces --skip-cleanup", completed.stderr)

    def test_incomplete_profile_never_falls_back_to_builtin_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            profile = valid_profile(root)
            del profile["datasets"]["operator-alias"]["inputs"]["to_somatic_vcf_pileup"]
            profile_path = root / "profile.json"
            profile_path.write_text(json.dumps(profile), encoding="utf-8")
            completed = subprocess.run(
                [
                    str(ROOT / "scripts" / "pipeline" / "run_benchmark.sh"),
                    "--site-profile", str(profile_path), "--sample", "S1",
                    "--mode", "to-pure", "--vcf-source", "pileup", "--plan",
                ],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("to_somatic_vcf_pileup", completed.stderr)
            self.assertNotIn("/big8_disk", completed.stdout + completed.stderr)

    def test_profile_output_root_traversal_is_rejected_by_realpath(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            profile = valid_profile(root)
            profile_path = root / "profile.json"
            profile_path.write_text(json.dumps(profile), encoding="utf-8")
            completed = subprocess.run(
                [
                    str(ROOT / "scripts" / "pipeline" / "run_benchmark.sh"),
                    "--site-profile", str(profile_path), "--sample", "S1", "--plan",
                    "--output-root", str(root / "output" / ".." / "escape"),
                ],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertEqual(completed.returncode, 3)
            self.assertIn("escapes data_roots.output", completed.stderr)

    def test_profile_project_root_must_match_executing_clone(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            profile = valid_profile(root)
            profile["project_root"] = str(root)
            profile_path = root / "profile.json"
            profile_path.write_text(json.dumps(profile), encoding="utf-8")
            completed = subprocess.run(
                [
                    str(ROOT / "scripts" / "pipeline" / "run_benchmark.sh"),
                    "--site-profile", str(profile_path), "--sample", "S1", "--plan",
                ],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertEqual(completed.returncode, 3)
            self.assertIn("does not match the executing clone", completed.stderr)

    def test_run_tag_traversal_is_rejected(self):
        completed = subprocess.run(
            [str(ROOT / "scripts" / "pipeline" / "run_benchmark.sh"), "--run-tag", "../../escape", "--plan"],
            text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        self.assertEqual(completed.returncode, 2)
        self.assertIn("Unsafe --run-tag", completed.stderr)

    def test_mode_and_vcf_source_are_closed_enums(self):
        runner = str(ROOT / "scripts" / "pipeline" / "run_benchmark.sh")
        for args, message in (
            (["--mode", "mystery", "--plan"], "Unsupported --mode"),
            (["--vcf-source", "mystery", "--plan"], "Unsupported --vcf-source"),
        ):
            completed = subprocess.run(
                [runner, *args], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False
            )
            self.assertEqual(completed.returncode, 2)
            self.assertIn(message, completed.stderr)

    def test_cleanup_requires_bound_sentinel_and_rejects_symlink_target(self):
        cleanup = ROOT / "scripts" / "pipeline" / "steps" / "04_cleanup.sh"
        with tempfile.TemporaryDirectory() as directory:
            allowed = Path(directory) / "allowed"
            run = allowed / "sample" / "run"
            lp = run / "longphase_s"
            lp.mkdir(parents=True)
            victim = Path(directory) / "victim.bam"
            victim.write_bytes(b"keep")
            (lp / "S1_tagged.bam").symlink_to(victim)

            missing = subprocess.run(
                [str(cleanup), "--output-dir", str(run), "--allowed-root", str(allowed)],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertEqual(missing.returncode, 3)
            self.assertEqual(victim.read_bytes(), b"keep")

            (run / ".intersubmod-run-root").write_text(
                f"schema_version=1\nrun_root={run.resolve()}\nallowed_root={allowed.resolve()}\n"
            )
            linked = subprocess.run(
                [str(cleanup), "--output-dir", str(run), "--allowed-root", str(allowed)],
                text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertEqual(linked.returncode, 3)
            self.assertIn("symlink", linked.stderr.lower())
            self.assertEqual(victim.read_bytes(), b"keep")

    def test_bootstrap_defaults_do_not_depend_on_cwd(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "created.json"
            completed = subprocess.run(
                [str(ROOT / "scripts" / "site" / "bootstrap"), "--output", str(output)],
                cwd="/", text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(json.loads(output.read_text())["project_root"], str(ROOT))


if __name__ == "__main__":
    unittest.main()
