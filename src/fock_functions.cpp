#include <armadillo>
#include <cmath>
#include <vector>

#include "molecule_visualization.hpp"


constexpr double EV_PER_HARTREE = 27.211324570273;

/// Fock Matrix Constructor ///

FockMatrix::FockMatrix(const Molecule& mol)
    : molecule(mol),
      p(mol.p),
      q(mol.q),
      gamma_matrix(static_cast<arma::uword>(mol.num_atoms),
                     static_cast<arma::uword>(mol.num_atoms),
                     arma::fill::zeros) {
    calculate_gamma_matrix();

    density_matrix_alpha.zeros(molecule.num_basis_functions, molecule.num_basis_functions);
    density_matrix_beta.zeros(molecule.num_basis_functions, molecule.num_basis_functions);
    fock_matrix_alpha.zeros(molecule.num_basis_functions, molecule.num_basis_functions);
    fock_matrix_beta.zeros(molecule.num_basis_functions, molecule.num_basis_functions);
    hamiltonian_matrix.zeros(molecule.num_basis_functions, molecule.num_basis_functions);

    assemble_hamiltonian_matrix();
}

/// Fock Matrix Creation ///

//boys function, where if T is 0, then the function is 1.0
double boys_f0(double T) {
    if (T < 1e-18) {
        return 1.0;
    }
    const double rt = std::sqrt(T);
    return std::sqrt(M_PI / (4.0 * T)) * std::erf(rt);
}

//uses boys function to calculate the coulomb integral
double coulomb_merged_ss(double P, double Q, double R) {
    const double rho = P * Q / (P + Q);
    const double T = rho * R * R;
    return (2.0 * M_PI * M_PI) / (P * Q) * std::sqrt(M_PI / (P + Q)) * boys_f0(T);
}

static double primitive_s_norm(double alpha) {
    return primitive_gaussian_normalization(alpha, 0, 0, 0);
}

//calculates the gamma matrix element for the valence s orbital
static double gamma_valence_ss_ev(const AOInfo& a, const AOInfo& b, double R_bohr) {
    double sum_au = 0.0;
    for (arma::uword k = 0; k < a.exponents.n_elem; ++k) {
        for (arma::uword kp = 0; kp < a.exponents.n_elem; ++kp) {
            const double P = a.exponents(k) + a.exponents(kp);
            const double wA = a.coefficients(k) * primitive_s_norm(a.exponents(k)) * a.coefficients(kp) *
                              primitive_s_norm(a.exponents(kp));
            for (arma::uword l = 0; l < b.exponents.n_elem; ++l) {
                for (arma::uword lp = 0; lp < b.exponents.n_elem; ++lp) {
                    const double Q = b.exponents(l) + b.exponents(lp);
                    const double wB = b.coefficients(l) * primitive_s_norm(b.exponents(l)) *
                                      b.coefficients(lp) * primitive_s_norm(b.exponents(lp));
                    sum_au += wA * wB * coulomb_merged_ss(P, Q, R_bohr);
                }
            }
        }
    }
    return sum_au * EV_PER_HARTREE;
}

static int valence_s_ao_index(const Molecule& m, int atom_idx) {
    int first_on_atom = -1;
    for (int i = 0; i < static_cast<int>(m.aos.size()); ++i) {
        if (m.aos[i].atom_idx != atom_idx) {
            continue;
        }
        if (first_on_atom < 0) {
            first_on_atom = i;
        }
        if (m.aos[i].l(0) == 0 && m.aos[i].l(1) == 0 && m.aos[i].l(2) == 0) {
            return i;
        }
    }
    return first_on_atom;
}

//finds valence s orbital and calculates the gamma matrix element for those 2 orbitals
static double gamma_matrix_element(const Molecule& m, int I, int J) {
    const int i_mu = valence_s_ao_index(m, I);
    const int i_nu = valence_s_ao_index(m, J);
    if (i_mu < 0 || i_nu < 0) {
            return 0.0;
        }
        const AOInfo& oi = m.aos[i_mu];
        const AOInfo& oj = m.aos[i_nu];
    const double R = arma::norm(m.atoms[I].center - m.atoms[J].center, 2);
    return gamma_valence_ss_ev(oi, oj, R);
}

//creates matrix
void FockMatrix::calculate_gamma_matrix() {
    const int n_atom = molecule.num_atoms;
    if (n_atom <= 0) {
        gamma_matrix.set_size(0, 0);
        return;
    }
    gamma_matrix.zeros(n_atom, n_atom);
    for (int I = 0; I < n_atom; ++I) {
        for (int J = 0; J < n_atom; ++J) {
            gamma_matrix(I, J) = gamma_matrix_element(molecule, I, J);
        }
    }
}



//s orbital is l = 0,0,0
static bool ao_is_valence_s(const AOInfo& ao) {
    return ao.l(0) == 0 && ao.l(1) == 0 && ao.l(2) == 0;
}

//determine which half core to use (s/p) from paramaters in atom
static double cndo2_core_diagonal_ev(const Atom& atom, bool is_s_shell) {
    const double half_sum_ia = is_s_shell ? atom.parameters.at(0) : atom.parameters.at(1);
    return -half_sum_ia;
}

void fock_diagonal_elements(const Molecule& mol, const arma::mat& gamma_atom, const arma::mat& P_tot,
                            const arma::mat& P_spin, arma::mat& F_spin) {
    const int n_ao = static_cast<int>(mol.aos.size());
    const int n_atom = mol.num_atoms;

    std::vector<double> p_on_atom(static_cast<size_t>(n_atom), 0.0);
    for (int mu = 0; mu < n_ao; ++mu) {
        const int A = mol.aos[mu].atom_idx;
        p_on_atom[static_cast<size_t>(A)] += P_tot(mu, mu);
    }

    for (int mu = 0; mu < n_ao; ++mu) {
        const int A = mol.aos[mu].atom_idx;
        const Atom& atom_a = mol.atoms[A];
        const bool is_s = ao_is_valence_s(mol.aos[mu]);

        const double core = cndo2_core_diagonal_ev(atom_a, is_s);
        const double z_a = atom_a.parameters.at(3);

        const double same_atom_factor =
            (p_on_atom[static_cast<size_t>(A)] - z_a) - (P_spin(mu, mu) - 0.5);
        double f_mu = core + same_atom_factor * gamma_atom(A, A);

        for (int C = 0; C < n_atom; ++C) {
            if (C == A) {
                continue;
            }
            const double z_c = mol.atoms[C].parameters.at(3);
            f_mu += (p_on_atom[static_cast<size_t>(C)] - z_c) * gamma_atom(A, C);
        }
        F_spin(mu, mu) = f_mu;
    }
}

void fock_off_diagonal_elements(const Molecule& mol, const arma::mat& gamma_atom, const arma::mat& overlap,
                                const arma::mat& P_spin, arma::mat& F_spin) {
    const int n_ao = static_cast<int>(mol.aos.size());

    for (int mu = 0; mu < n_ao; ++mu) {
        for (int nu = mu + 1; nu < n_ao; ++nu) {
            const int A = mol.aos[mu].atom_idx;
            const int B = mol.aos[nu].atom_idx;
            const double beta_a = mol.atoms[A].parameters.at(2);
            const double beta_b = mol.atoms[B].parameters.at(2);
            const double s_mu_nu = overlap(mu, nu);
            
            const double val =
                -0.5 * (beta_a + beta_b) * s_mu_nu - P_spin(mu, nu) * gamma_atom(A, B);
            F_spin(mu, nu) = val;
            F_spin(nu, mu) = val;
        }
    }
}

//p spin is the density matrix and f out is the fock matrix :D
void FockMatrix::assemble_fock_spin(const arma::mat& P_spin, arma::mat& F_out) {
    const arma::mat P_tot = density_matrix_alpha + density_matrix_beta;
    fock_diagonal_elements(molecule, gamma_matrix, P_tot, P_spin, F_out);
    fock_off_diagonal_elements(molecule, gamma_matrix, molecule.overlap_matrix, P_spin, F_out);
}




/// Eigenvalue solver ///

//get the updated molecular coefficients (𝑪𝜶, 𝑪𝜷) and corresponding eigenvalues (𝝐𝜶, 𝝐𝜷) from the fock matrix
void FockMatrix::solve_eigenvalue_problem() {

    prev_density_matrix_alpha = density_matrix_alpha;
    prev_density_matrix_beta = density_matrix_beta;

    const arma::uword n = fock_matrix_alpha.n_rows;

    arma::mat eigenvectors_alpha;
    arma::mat eigenvectors_beta;

    arma::eig_sym(eigenvalues_alpha, eigenvectors_alpha, fock_matrix_alpha);
    arma::eig_sym(eigenvalues_beta, eigenvectors_beta, fock_matrix_beta);

    density_matrix_alpha.zeros(n, n);
    density_matrix_beta.zeros(n, n);

    if (p > 0) {
        const arma::uword pu = static_cast<arma::uword>(p);
        const arma::mat C_occ_a = eigenvectors_alpha.cols(0, pu - 1);
        density_matrix_alpha = C_occ_a * C_occ_a.t();
    }
    if (q > 0) {
        const arma::uword qu = static_cast<arma::uword>(q);
        const arma::mat C_occ_b = eigenvectors_beta.cols(0, qu - 1);
        density_matrix_beta = C_occ_b * C_occ_b.t();
    }

    eigen++;
}


bool FockMatrix::check_convergence() {
    const double convergence_threshold = 1e-6;
    //need both to be near equal
    return arma::approx_equal(density_matrix_alpha, prev_density_matrix_alpha, "absdiff", convergence_threshold) &&
           arma::approx_equal(density_matrix_beta, prev_density_matrix_beta, "absdiff", convergence_threshold);
}


void FockMatrix::find_convergence(){
    while (true) {
        calculate_fock_spin();
        solve_eigenvalue_problem();
        if (check_convergence()) {
            break;
        }
    }
}


/// Energy Calculation ///

void hamiltonian_diagonal_elements(const Molecule& mol, const arma::mat& gamma_matrix, arma::mat& hamiltonian_matrix) {
    const int n_ao = static_cast<int>(mol.aos.size());
    const int n_atom = mol.num_atoms;

    for (int mu = 0; mu < n_ao; ++mu) {
        const int A = mol.aos[mu].atom_idx;
        const Atom& atom_a = mol.atoms[A];
        
        const bool is_s = ao_is_valence_s(mol.aos[mu]);
        const double core = cndo2_core_diagonal_ev(atom_a, is_s);

        const double z_a = atom_a.parameters.at(3);

        const double same_atom_factor = (z_a - 0.5) * gamma_matrix(A, A);
        double h_mu = core - same_atom_factor;

        for (int B = 0; B < n_atom; ++B) {
            if (B == A) {
            continue;
            }
            const double z_b = mol.atoms[B].parameters.at(3);
            h_mu -= z_b * gamma_matrix(A, B);
        }
        hamiltonian_matrix(mu, mu) = h_mu;
    }
}

void hamiltonian_off_diagonal_elements(const Molecule& mol, const arma::mat& gamma_matrix, const arma::mat& overlap, arma::mat& hamiltonian_matrix) {
    
    const int n_ao = static_cast<int>(mol.aos.size());

    for (int mu = 0; mu < n_ao; ++mu) {
        for (int nu = mu + 1; nu < n_ao; ++nu) {
            const int A = mol.aos[mu].atom_idx;
            const int B = mol.aos[nu].atom_idx;
            const double beta_a = mol.atoms[A].parameters.at(2);
            const double beta_b = mol.atoms[B].parameters.at(2);

            const double s_mu_nu = overlap(mu, nu);

            const double val = - 0.5 * (beta_a + beta_b) * s_mu_nu;
            hamiltonian_matrix(mu, nu) = val;
            hamiltonian_matrix(nu, mu) = val;
        }
    }
}


void FockMatrix::assemble_hamiltonian_matrix(){
    hamiltonian_diagonal_elements(molecule, gamma_matrix, hamiltonian_matrix);
    hamiltonian_off_diagonal_elements(molecule, gamma_matrix, molecule.overlap_matrix, hamiltonian_matrix);
}


void FockMatrix::calculate_repulsion_energy(){
    for(int A = 0; A < molecule.num_atoms; ++A){
        for(int B = 0; B < A; ++B){

            const double z_a = molecule.atoms[A].parameters.at(3);
            const double z_b = molecule.atoms[B].parameters.at(3);
            const double R = arma::norm(molecule.atoms[A].center - molecule.atoms[B].center, 2);
            repulsion_energy += z_a * z_b / R * EV_PER_HARTREE;
        }
    }
}

void FockMatrix::calculate_electon_energy(){
    electron_energy = 0.0;

    for(int mu = 0; mu < molecule.num_basis_functions; ++mu){
        for(int nu = 0; nu < molecule.num_basis_functions; ++nu){
            electron_energy += (hamiltonian_matrix(mu, nu) + fock_matrix_alpha(mu, nu)) * density_matrix_alpha(mu, nu);
            electron_energy += (hamiltonian_matrix(mu, nu) + fock_matrix_beta(mu, nu)) * density_matrix_beta(mu, nu);
        }
    }
    electron_energy *= 0.5;
}

void FockMatrix::calculate_total_energy(){

    electron_energy = 0.0;
    repulsion_energy = 0.0;
    total_energy = 0.0;

    calculate_repulsion_energy();
    calculate_electon_energy();
    total_energy = electron_energy + repulsion_energy;
}


/////////////////////////////////////////////////

double cndo_gamma_ij_ev(const Molecule& mol, int I, int J) { 
    return gamma_matrix_element(mol, I, J); 
}

double cndo_gamma_ij_dgamma_dR_ev(const Molecule& mol, int I, int J) {
    if (I == J) {
        return 0.0;
    }
    const int i_mu = valence_s_ao_index(mol, I);
    const int i_nu = valence_s_ao_index(mol, J);
    if (i_mu < 0 || i_nu < 0) {
        return 0.0;
    }
    const AOInfo& oi = mol.aos[i_mu];
    const AOInfo& oj = mol.aos[i_nu];
    const double R = arma::norm(mol.atoms[I].center - mol.atoms[J].center, 2);
    constexpr double h = 1e-8;
    return (gamma_valence_ss_ev(oi, oj, R + h) - gamma_valence_ss_ev(oi, oj, R - h)) / (2.0 * h);
}



////// 
arma::mat FockMatrix::find_gammaAB_RA()
{
  const int num_atoms = molecule.num_atoms;
  arma::mat gammaAB_RA(3, num_atoms * num_atoms, arma::fill::zeros);

for (int I = 0; I < molecule.num_atoms; ++I) { //iterate through each atom
  for (int J = 0; J < molecule.num_atoms; ++J) {
    
    const int col = I * molecule.num_atoms + J;
    if (I == J) {
      continue;
    }
    arma::vec r_ab = molecule.atoms[I].center - molecule.atoms[J].center;
    const double r = arma::norm(r_ab, 2);
    if (r < 1e-14) { //skip if the distance is too small
      continue;
    }
    r_ab /= r;
    const double dg = cndo_gamma_ij_dgamma_dR_ev(molecule, I, J);

    gammaAB_RA(0, col) = dg * r_ab(0);  //
    gammaAB_RA(1, col) = dg * r_ab(1);
    gammaAB_RA(2, col) = dg * r_ab(2);
    }
    }

    return gammaAB_RA;
}