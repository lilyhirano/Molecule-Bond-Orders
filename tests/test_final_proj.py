from pathlib import Path

import numpy as np
import pytest


def _read_csv_matrix(csv_path: Path) -> np.ndarray:
    rows = []
    for line in csv_path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append([float(value) for value in line.split(",")])
    if not rows:
        raise AssertionError(f"Reference CSV is empty: {csv_path}")
    return np.array(rows, dtype=float)


def _parse_stdout_matrix(stdout_text: str, label: str) -> np.ndarray:
    lines = stdout_text.splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if line.strip() == label) + 1
    except StopIteration as exc:
        raise AssertionError(f"Could not find matrix label '{label}' in stdout") from exc

    matrix_rows = []
    for line in lines[start:]:
        if not line.strip():
            break
        matrix_rows.append([float(value) for value in line.split()])

    if not matrix_rows:
        raise AssertionError(f"No matrix rows found under label '{label}'")

    return np.array(matrix_rows, dtype=float)


@pytest.mark.parametrize("molecule_id", ["H2", "HF", "HO", "H2O"])
def test_wiberg_bond_order_near_reference(student_output_dir, molecule_id):
    base_dir = Path(__file__).resolve().parent.parent
    csv_path = base_dir / "tests" / "test_checks" / f"{molecule_id}.csv"
    stdout_path = Path(student_output_dir) / f"{molecule_id}.stdout"

    assert csv_path.exists(), f"Missing reference file: {csv_path}"
    assert stdout_path.exists(), f"Missing stdout file: {stdout_path}"

    expected = _read_csv_matrix(csv_path)
    actual = _parse_stdout_matrix(
        stdout_path.read_text(encoding="utf-8"),
        label="Wiberg bond order matrix:",
    )

    assert actual.shape == expected.shape, (
        f"Shape mismatch for {molecule_id}: got {actual.shape}, expected {expected.shape}"
    )
    np.testing.assert_allclose(actual, expected, atol=0.35, rtol=0.25)

