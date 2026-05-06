import pathlib
import subprocess
import json

import pytest

BASE_DIR = pathlib.Path(__file__).parent.parent
INPUT_DIR = BASE_DIR / "sample_input" / "final_proj"
BINARY_DIR = BASE_DIR / "build"
EXECUTABLE_PATH = BINARY_DIR / "molecule_bond_order"
PREFERRED_OUTPUT_SUBDIR = BASE_DIR / "project_output" / "final_proj"
PREFERRED_HDF5_SUBDIR = BASE_DIR / "project_output" / "final_proj_hdf5"
FALLBACK_OUTPUT_SUBDIR = pathlib.Path("/tmp/qc279_pytest_output/final_proj")
FALLBACK_HDF5_SUBDIR = pathlib.Path("/tmp/qc279_pytest_output/final_proj_hdf5")
RUN_CONFIG_SUBDIR = pathlib.Path("/tmp/qc279_pytest_output/configs")
RUN_OUTPUT_SUBDIR = PREFERRED_OUTPUT_SUBDIR
RUN_HDF5_SUBDIR = PREFERRED_HDF5_SUBDIR

@pytest.fixture(scope="session")
def student_output_dir():
    return RUN_OUTPUT_SUBDIR


@pytest.fixture(scope="session")
def student_hdf5_dir():
    return RUN_HDF5_SUBDIR

def pytest_report_header(config):
    return "Running QC279 pytest suite"

def pytest_configure():
    global RUN_OUTPUT_SUBDIR, RUN_HDF5_SUBDIR
    try:
        PREFERRED_OUTPUT_SUBDIR.mkdir(parents=True, exist_ok=True)
        PREFERRED_HDF5_SUBDIR.mkdir(parents=True, exist_ok=True)
        RUN_OUTPUT_SUBDIR = PREFERRED_OUTPUT_SUBDIR
        RUN_HDF5_SUBDIR = PREFERRED_HDF5_SUBDIR
    except PermissionError:
        RUN_OUTPUT_SUBDIR = FALLBACK_OUTPUT_SUBDIR
        RUN_HDF5_SUBDIR = FALLBACK_HDF5_SUBDIR
    run_all_test_cases()

def run_all_test_cases():
    RUN_OUTPUT_SUBDIR.mkdir(parents=True, exist_ok=True)
    RUN_HDF5_SUBDIR.mkdir(parents=True, exist_ok=True)
    RUN_CONFIG_SUBDIR.mkdir(parents=True, exist_ok=True)

    if not EXECUTABLE_PATH.exists():
        raise FileNotFoundError(f"Expected executable not found: {EXECUTABLE_PATH}")
    if not INPUT_DIR.exists():
        raise FileNotFoundError(f"Expected input directory not found: {INPUT_DIR}")

    for test_file in sorted(INPUT_DIR.glob("*.json")):
        with test_file.open("r", encoding="utf-8") as f:
            config = json.load(f)

        # Force writable output path inside /tmp during test execution.
        config["output_file_path"] = str((RUN_HDF5_SUBDIR / f"{test_file.stem}.hdf5").resolve())
        run_config_path = RUN_CONFIG_SUBDIR / test_file.name
        run_config_path.write_text(json.dumps(config, indent=4), encoding="utf-8")

        result = subprocess.run(
            args=[str(EXECUTABLE_PATH), str(run_config_path)],
            capture_output=True,
            text=True,
            check=True,
        )
        student_output_path = RUN_OUTPUT_SUBDIR / test_file.stem
        student_output_path.with_suffix(".stdout").write_text(result.stdout)