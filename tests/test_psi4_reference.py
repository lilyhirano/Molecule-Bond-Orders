import json
from pathlib import Path

import numpy as np
import pytest


Z_TO_SYMBOL = {
    1: "H",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
}


def _parse_matrix_block(stdout_text: str, label: str) -> np.ndarray:
    lines = stdout_text.splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if line.strip() == label) + 1
    except StopIteration as exc:
        raise AssertionError(f"Could not find matrix label: {label}") from exc

    matrix_lines = []
    for line in lines[start:]:
        if not line.strip():
            break
        matrix_lines.append(line.strip())

    if not matrix_lines:
        raise AssertionError(f"No matrix rows found for label: {label}")

    rows = [[float(x) for x in row.split()] for row in matrix_lines]
    return np.array(rows, dtype=float)


def _read_hdf5_matrix(hdf5_path: Path, dataset_name: str) -> np.ndarray:
    h5py = pytest.importorskip("h5py")
    if not hdf5_path.exists():
        raise AssertionError(f"Missing hdf5 file: {hdf5_path}")
    with h5py.File(hdf5_path, "r") as f:
        if dataset_name not in f:
            raise AssertionError(f"Missing dataset '{dataset_name}' in {hdf5_path}")
        return np.asarray(f[dataset_name], dtype=float)


def _read_xyz_atomic_numbers_and_coords(xyz_path: Path):
    lines = xyz_path.read_text(encoding="utf-8").splitlines()
    atom_lines = lines[2:]
    atoms = []
    for line in atom_lines:
        if not line.strip():
            continue
        z_str, x_str, y_str, zc_str = line.split()
        z = int(z_str)
        atoms.append((z, float(x_str), float(y_str), float(zc_str)))
    return atoms


def _build_psi4_molecule_string(atoms, charge: int, multiplicity: int) -> str:
    atom_lines = []
    for z, x, y, zc in atoms:
        symbol = Z_TO_SYMBOL[z]
        atom_lines.append(f"{symbol} {x} {y} {zc}")
    return "\n".join(
        [
            f"{charge} {multiplicity}",
            *atom_lines,
            "units bohr",
            "symmetry c1",
            "no_reorient",
            "no_com",
        ]
    )


def _compute_reference_bond_orders(psi4, atoms, nalpha: int, nbeta: int):
    total_z = sum(z for z, *_ in atoms)
    charge = total_z - (nalpha + nbeta)
    multiplicity = nalpha - nbeta + 1
    if multiplicity < 1:
        raise AssertionError(f"Invalid multiplicity computed: {multiplicity}")

    psi4.core.be_quiet()
    psi4.set_options(
        {
            "basis": "sto-3g",
            "scf_type": "pk",
            "e_convergence": 1e-10,
            "d_convergence": 1e-10,
            "reference": "rhf" if nalpha == nbeta else "uhf",
        }
    )

    mol = psi4.geometry(_build_psi4_molecule_string(atoms, charge, multiplicity))
    _, wfn = psi4.energy("scf", molecule=mol, return_wfn=True)

    S = np.asarray(wfn.S())
    Ca = np.asarray(wfn.Ca())
    Cb = np.asarray(wfn.Cb())
    Pa = Ca[:, :nalpha] @ Ca[:, :nalpha].T
    Pb = Cb[:, :nbeta] @ Cb[:, :nbeta].T
    Ptot = Pa + Pb

    n_atoms = len(atoms)
    nbf = S.shape[0]
    atom_to_aos = [[] for _ in range(n_atoms)]
    basis = wfn.basisset()
    for mu in range(nbf):
        atom_to_aos[basis.function_to_center(mu)].append(mu)

    mulliken = np.zeros((n_atoms, n_atoms))
    mayer = np.zeros((n_atoms, n_atoms))
    wiberg = np.zeros((n_atoms, n_atoms))

    PSa = Pa @ S
    PSb = Pb @ S
    PStot = Ptot @ S

    for i in range(n_atoms):
        for j in range(n_atoms):
            if i == j:
                continue
            for mu in atom_to_aos[i]:
                for nu in atom_to_aos[j]:
                    mulliken[i, j] += Ptot[mu, nu] * S[mu, nu]
                    mayer[i, j] += PSa[mu, nu] * PSa[nu, mu] + PSb[mu, nu] * PSb[nu, mu]
                    wiberg[i, j] += PStot[mu, nu] * PStot[nu, mu]

    return mulliken, mayer, wiberg


def _load_cpp_and_reference_bond_orders(student_hdf5_dir, molecule_id):
    psi4 = pytest.importorskip("psi4")
    repo_root = Path(__file__).resolve().parent.parent

    hdf5_path = Path(student_hdf5_dir) / f"{molecule_id}.hdf5"
    mulliken_cpp = _read_hdf5_matrix(hdf5_path, "mulliken_bond_order_matrix")
    mayer_cpp = _read_hdf5_matrix(hdf5_path, "mayer_bond_order_matrix")
    wiberg_cpp = _read_hdf5_matrix(hdf5_path, "wiberg_bond_order_matrix")

    config_path = repo_root / "sample_input" / "final_proj" / f"{molecule_id}.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    atoms_path = repo_root / config["atoms_file_path"]
    atoms = _read_xyz_atomic_numbers_and_coords(atoms_path)

    mulliken_ref, mayer_ref, wiberg_ref = _compute_reference_bond_orders(
        psi4,
        atoms,
        int(config["num_alpha_electrons"]),
        int(config["num_beta_electrons"]),
    )

    return (
        mulliken_cpp,
        mayer_cpp,
        wiberg_cpp,
        mulliken_ref,
        mayer_ref,
        wiberg_ref,
    )


@pytest.mark.parametrize("molecule_id", ["H2", "HF", "HO", "H2O"])
def test_mulliken_bond_orders_against_psi4(student_hdf5_dir, molecule_id):
    mulliken_cpp, _, _, mulliken_ref, _, _ = _load_cpp_and_reference_bond_orders(
        student_hdf5_dir, molecule_id
    )
    np.testing.assert_allclose(mulliken_cpp, mulliken_ref, atol=0.35, rtol=0.25)


@pytest.mark.parametrize("molecule_id", ["H2", "HF", "HO", "H2O"])
def test_mayer_bond_orders_against_psi4(student_hdf5_dir, molecule_id):
    _, mayer_cpp, _, _, mayer_ref, _ = _load_cpp_and_reference_bond_orders(
        student_hdf5_dir, molecule_id
    )
    np.testing.assert_allclose(mayer_cpp, mayer_ref, atol=0.35, rtol=0.25)


@pytest.mark.parametrize("molecule_id", ["H2", "HF", "HO", "H2O"])
def test_wiberg_bond_orders_against_psi4(student_hdf5_dir, molecule_id):
    _, _, wiberg_cpp, _, _, wiberg_ref = _load_cpp_and_reference_bond_orders(
        student_hdf5_dir, molecule_id
    )
    np.testing.assert_allclose(wiberg_cpp, wiberg_ref, atol=0.35, rtol=0.25)
