#!/usr/bin/env python3
import hashlib
import json
import os
import pathlib
import subprocess
import sys
import tempfile
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_hypercube_reductions_exhaustive.py"


class ExhaustiveReductionAuditTest(unittest.TestCase):
    def test_full_deterministic_audit_and_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            output = pathlib.Path(directory) / "receipt.json"
            environment = dict(os.environ)
            environment.update({
                "OPENBLAS_NUM_THREADS": "1",
                "OMP_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
                "NUMEXPR_NUM_THREADS": "1",
            })
            completed = subprocess.run(
                [sys.executable, str(SCRIPT), "--json-output", str(output)],
                cwd=ROOT.parents[1],
                env=environment,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
            receipt = json.loads(output.read_text(encoding="utf-8"))
            self.assertTrue(receipt["all_pass"])
            self.assertEqual(receipt["presolve"]["presolve_cases"], 61_340)
            self.assertEqual(
                receipt["presolve"]["selected_set_predicate_checks"], 1_979_356
            )
            self.assertEqual(
                receipt["sparse_fixed_cardinality_no_good"]["pairs"], 23_909
            )
            self.assertEqual(receipt["general_dense_no_good"]["pairs"], 21_844)
            self.assertEqual(receipt["presolve"]["mismatches"], 0)
            sidecar = output.with_name(output.name + ".sha256")
            expected = hashlib.sha256(output.read_bytes()).hexdigest()
            self.assertEqual(
                sidecar.read_text(encoding="utf-8").strip(),
                f"{expected}  {output.name}",
            )


if __name__ == "__main__":
    unittest.main()
