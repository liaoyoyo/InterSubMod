from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_matched_normal_methyl_tags.py"
SPEC = importlib.util.spec_from_file_location("audit_matched_normal_methyl_tags", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class FakeRead:
    def __init__(self, tags: dict[str, object]) -> None:
        self.tags = tags

    def has_tag(self, name: str) -> bool:
        return name in self.tags

    def get_tag(self, name: str) -> object:
        return self.tags[name]


def test_mm_aliases_are_supported() -> None:
    assert MODULE.has_any_tag(FakeRead({"MM": "C+m,0;"}), ("MM", "Mm"))
    assert MODULE.has_any_tag(FakeRead({"Mm": "C+m,0;"}), ("MM", "Mm"))
    assert not MODULE.has_any_tag(FakeRead({}), ("MM", "Mm"))


def test_mm_codes_are_parsed_without_deltas() -> None:
    assert MODULE.summarize_mm_tag("C+m?,1,2;C+h?,3;") == ["C+m?", "C+h?"]


def test_index_candidates_include_bam_dot_bai(tmp_path: Path) -> None:
    bam = tmp_path / "normal.bam"
    bam.write_bytes(b"bam")
    index = tmp_path / "normal.bam.bai"
    index.write_bytes(b"bai")
    assert MODULE.find_index(bam) == index


def test_default_output_is_a_versioned_receipt() -> None:
    source = SCRIPT.read_text(encoding="utf-8")
    assert "Refusing to overwrite existing audit receipt" in source
