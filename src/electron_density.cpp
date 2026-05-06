# include "molecule_bond_order.hpp"

void ElectronDensity::calculate_electron_density(){
    const int num_atoms = molecule.num_atoms;
    const int n_ao = static_cast<int>(molecule.aos.size());
    const arma::mat& S = molecule.overlap_matrix;

    electron_density_matrix.zeros(num_atoms, 1);

    // Mulliken atomic population:
    // rho_A = sum_{mu in A} sum_nu P_{mu,nu} S_{mu,nu}
    for (int atom = 0; atom < num_atoms; ++atom) {
        double atom_density = 0.0;
        for (int mu = 0; mu < n_ao; ++mu) {
            if (molecule.aos[mu].atom_idx != atom) {
                continue;
            }
            for (int nu = 0; nu < n_ao; ++nu) {
                atom_density += density_matrix_total(mu, nu) * S(mu, nu);
            }
        }
        electron_density_matrix(atom, 0) = atom_density;
    }
}