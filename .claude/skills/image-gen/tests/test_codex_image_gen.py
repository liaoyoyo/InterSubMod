"""Test codex_image_gen.build_command builds correct codex CLI invocation."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "tools"))
from codex_image_gen import build_command, extract_prompt_text  # noqa

FIXTURES = Path(__file__).parent / "fixtures"

def test_extract_concept_prompt():
    yaml_path = FIXTURES / "valid_concept.yaml"
    text = extract_prompt_text(yaml_path)
    assert "LOH region" in text, f"Expected LOH in prompt, got: {text!r}"
    assert "$imagegen" in text, f"Expected $imagegen suffix, got: {text!r}"
    assert "no gradient" in text.lower(), "Expected constraints embedded"

def test_build_command():
    yaml_path = FIXTURES / "valid_concept.yaml"
    out_dir = Path("/tmp/test_codex_out")
    cmd = build_command(yaml_path, out_dir)
    assert cmd[0] == "codex", f"First token must be 'codex', got {cmd[0]}"
    assert "exec" in cmd, "Must use 'exec' subcommand"
    assert "--image-dir" in cmd, "Must specify --image-dir"
    idx = cmd.index("--image-dir")
    assert cmd[idx+1] == str(out_dir), f"--image-dir value wrong: {cmd[idx+1]}"
    last = cmd[-1]
    assert "$imagegen" in last, f"Prompt argument must end with $imagegen: {last!r}"

def main():
    tests = [test_extract_concept_prompt, test_build_command]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS: {t.__name__}")
        except AssertionError as e:
            print(f"FAIL: {t.__name__}: {e}")
            failed += 1
    return 1 if failed else 0

if __name__ == "__main__":
    sys.exit(main())
