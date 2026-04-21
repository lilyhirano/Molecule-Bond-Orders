#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <armadillo>
#include <nlohmann/json.hpp>
#include "molecule_visualization.hpp"

using namespace std;
using json = nlohmann::json; 

void Molecule::read_atoms_file(){
    ifstream atoms_file(atoms_file_path_);
    string line;
    
    getline(atoms_file, line);
    num_atoms = stoi(line);

    getline(atoms_file, line); //ignore  comment line

    atoms.resize(num_atoms);
    num_element_symbols.resize(5, 0);  // [H, C, N, O, F]
    int i = 0;
    while(getline(atoms_file, line)){
        istringstream iss(line);
        int element_symbol;
        double x, y, z;
        iss >> element_symbol >> x >> y >> z;
        atoms[i].element_symbol = element_symbol;
        atoms[i].center = arma::vec({x, y, z});
        i++;
    }
    //determine num of each element
    for (int i = 0; i < num_atoms; ++i){
        if (atoms[i].element_symbol == 1){ //H
            num_element_symbols[0]++;
            atoms[i].parameters = {7.176, 0, 9.0, 1.0};
        } else if (atoms[i].element_symbol == 6){ //C
            num_element_symbols[1]++;
            atoms[i].parameters = {14.051, 5.572, 21.0, 4.0};
        } else if (atoms[i].element_symbol == 7){ //N
            num_element_symbols[2]++;
            atoms[i].parameters = {19.316, 7.275, 25.0, 5.0};
        } else if (atoms[i].element_symbol == 8){ //O
            num_element_symbols[3]++;
            atoms[i].parameters = {25.390, 9.111, 31.0, 6.0};
        } else if (atoms[i].element_symbol == 9){ //F
            atoms[i].parameters = {32.272, 11.080, 39.0, 7.0};
            num_element_symbols[4]++;
        }
        else {
            throw runtime_error("Invalid element symbol: " + to_string(atoms[i].element_symbol));
        }
    }
}

void Molecule::determine_num_basis_functions(){ 
    const int nheavy =
        num_element_symbols[1] + num_element_symbols[2] + num_element_symbols[3] +
        num_element_symbols[4];
    num_basis_functions = num_element_symbols[0] + 4 * nheavy;
}

//creates the AOInfo struct for each atom
void Molecule::build_aos(const STO3G_maps& sto3g){
    vector<AOInfo> new_aos;
    new_aos.reserve(static_cast<size_t>(num_basis_functions));

    static const vector<arma::ivec> l_p = {
        arma::ivec({1, 0, 0}), arma::ivec({0, 1, 0}), arma::ivec({0, 0, 1})};

    for (int i = 0; i < num_atoms; ++i) {
        const Atom& atom = atoms[i];

        if (atom.element_symbol == 1) { // H: 1xs
            AOInfo ao;
            ao.atom_idx = i;
            ao.type = "H_s";
            ao.center = atom.center;
            ao.exponents = sto3g.H_s.at("exponents");
            ao.coefficients = sto3g.H_s.at("contraction_coefficients");
            ao.l = arma::ivec({0, 0, 0});
            new_aos.push_back(std::move(ao));
        } else if (atom.element_symbol == 6) { // C: 1×s + 3×p
            AOInfo ao_s;
            ao_s.atom_idx = i;
            ao_s.type = "C_s";
            ao_s.center = atom.center;
            ao_s.exponents = sto3g.C_s.at("exponents");
            ao_s.coefficients = sto3g.C_s.at("contraction_coefficients");
            ao_s.l = arma::ivec({0, 0, 0});
            new_aos.push_back(std::move(ao_s));

            static const vector<string> p_types = {"C_p_x", "C_p_y", "C_p_z"};
            for (int pi = 0; pi < 3; ++pi) {
                AOInfo ao_p;
                ao_p.atom_idx = i;
                ao_p.type = p_types[pi];
                ao_p.center = atom.center;
                ao_p.exponents = sto3g.C_p.at("exponents");
                ao_p.coefficients = sto3g.C_p.at("contraction_coefficients");
                ao_p.l = l_p[pi];
                new_aos.push_back(std::move(ao_p));
            }
        } else if (atom.element_symbol == 7) { // N: 1×s + 3×p (minimal STO-3G)
            AOInfo ao_s;
            ao_s.atom_idx = i;
            ao_s.type = "N_s";
            ao_s.center = atom.center;
            ao_s.exponents = sto3g.N_s.at("exponents");
            ao_s.coefficients = sto3g.N_s.at("contraction_coefficients");
            ao_s.l = arma::ivec({0, 0, 0});
            new_aos.push_back(std::move(ao_s));

            static const vector<string> p_types = {"N_p_x", "N_p_y", "N_p_z"};
            for (int pi = 0; pi < 3; ++pi) {
                AOInfo ao_p;
                ao_p.atom_idx = i;
                ao_p.type = p_types[pi];
                ao_p.center = atom.center;
                ao_p.exponents = sto3g.N_p.at("exponents");
                ao_p.coefficients = sto3g.N_p.at("contraction_coefficients");
                ao_p.l = l_p[pi];
                new_aos.push_back(std::move(ao_p));
            }
        } else if (atom.element_symbol == 8) { // O: 1×s + 3×p (minimal STO-3G)
            AOInfo ao_s;
            ao_s.atom_idx = i;
            ao_s.type = "O_s";
            ao_s.center = atom.center;
            ao_s.exponents = sto3g.O_s.at("exponents");
            ao_s.coefficients = sto3g.O_s.at("contraction_coefficients");
            ao_s.l = arma::ivec({0, 0, 0});
            new_aos.push_back(std::move(ao_s));

            static const vector<string> p_types = {"O_p_x", "O_p_y", "O_p_z"};
            for (int pi = 0; pi < 3; ++pi) {
                AOInfo ao_p;
                ao_p.atom_idx = i;
                ao_p.type = p_types[pi];
                ao_p.center = atom.center;
                ao_p.exponents = sto3g.O_p.at("exponents");
                ao_p.coefficients = sto3g.O_p.at("contraction_coefficients");
                ao_p.l = l_p[pi];
                new_aos.push_back(std::move(ao_p));
            }
        } else if (atom.element_symbol == 9) { // F: 1×s + 3×p (minimal STO-3G)
            AOInfo ao_s;
            ao_s.atom_idx = i;
            ao_s.type = "F_s";
            ao_s.center = atom.center;
            ao_s.exponents = sto3g.F_s.at("exponents");
            ao_s.coefficients = sto3g.F_s.at("contraction_coefficients");
            ao_s.l = arma::ivec({0, 0, 0});
            new_aos.push_back(std::move(ao_s));

            static const vector<string> p_types = {"F_p_x", "F_p_y", "F_p_z"};
            for (int pi = 0; pi < 3; ++pi) {
                AOInfo ao_p;
                ao_p.atom_idx = i;
                ao_p.type = p_types[pi];
                ao_p.center = atom.center;
                ao_p.exponents = sto3g.F_p.at("exponents");
                ao_p.coefficients = sto3g.F_p.at("contraction_coefficients");
                ao_p.l = l_p[pi];
                new_aos.push_back(std::move(ao_p));
            }
        } else {
            throw runtime_error("Unknown atomic number: " + to_string(atom.element_symbol));
        }
    }

    aos = std::move(new_aos);
}

double primitive_gaussian_normalization(double alpha, int l, int m, int n) {
    int L = l + m + n;
    double fac = std::pow(2.0, 2 * L + 1.5)
               * std::pow(alpha, L + 1.5)
               / (std::pow(M_PI, 1.5)
                  * std::tgamma(l + 1)
                  * std::tgamma(m + 1)
                  * std::tgamma(n + 1));
    return std::sqrt(fac);
}

double cartesian_overlap_1d(int la, int lb, double PAx, double PBx, double alpha, double beta, double ABx) {
    if (la < 0 || lb < 0) return 0.0;
    if (la == 0 && lb == 0) return 1.0;
    if (la == 0) return PBx * cartesian_overlap_1d(la, lb-1, PAx, PBx, alpha, beta, ABx) + (lb-1)/(2*(alpha+beta)) * cartesian_overlap_1d(la, lb-2, PAx, PBx, alpha, beta, ABx);
    if (lb == 0) return PAx * cartesian_overlap_1d(la-1, lb, PAx, PBx, alpha, beta, ABx) + (la-1)/(2*(alpha+beta)) * cartesian_overlap_1d(la-2, lb, PAx, PBx, alpha, beta, ABx);
    return PAx * cartesian_overlap_1d(la-1, lb, PAx, PBx, alpha, beta, ABx) + (la-1)/(2*(alpha+beta)) * cartesian_overlap_1d(la-2, lb, PAx, PBx, alpha, beta, ABx)
        + lb/(2*(alpha+beta)) * cartesian_overlap_1d(la-1, lb-1, PAx, PBx, alpha, beta, ABx);
}

double primitive_gaussian_overlap(
    double alpha, const arma::vec &A, int l1, int m1, int n1,
    double beta, const arma::vec &B, int l2, int m2, int n2
) {
    double gamma = alpha + beta;
    arma::vec P = (alpha * A + beta * B) / gamma;
    double AB2 = arma::accu(arma::square(A - B));
    double pre = std::pow(M_PI / gamma, 1.5) * std::exp(-alpha * beta / gamma * AB2);

    double Sx = cartesian_overlap_1d(l1, l2, P(0)-A(0), P(0)-B(0), alpha, beta, A(0)-B(0));
    double Sy = cartesian_overlap_1d(m1, m2, P(1)-A(1), P(1)-B(1), alpha, beta, A(1)-B(1));
    double Sz = cartesian_overlap_1d(n1, n2, P(2)-A(2), P(2)-B(2), alpha, beta, A(2)-B(2));

    double N1 = primitive_gaussian_normalization(alpha, l1, m1, n1);
    double N2 = primitive_gaussian_normalization(beta,  l2, m2, n2);

    return N1 * N2 * pre * Sx * Sy * Sz;
}

double contracted_gaussian_overlap(
    const arma::vec &A, const arma::vec &alphas, const arma::vec &coeffsA, int lA, int mA, int nA,
    const arma::vec &B, const arma::vec &betas,  const arma::vec &coeffsB, int lB, int mB, int nB
) {
    double S = 0.0;
    for (size_t i = 0; i < alphas.n_elem; ++i){
        for (size_t j = 0; j < betas.n_elem; ++j){
            S += coeffsA(i) * coeffsB(j) * primitive_gaussian_overlap(
                alphas(i), A, lA, mA, nA,
                betas(j), B, lB, mB, nB
            );
        }
    }
    return S;
}


void Molecule::build_overlap_matrix(){
    int N = aos.size();
    arma::mat S(N, N, arma::fill::zeros);

    for (int mu = 0; mu < N; ++mu) {
        for (int nu = 0; nu < N; ++nu) {
            int l1=aos[mu].l(0), m1=aos[mu].l(1), n1=aos[mu].l(2);
            int l2=aos[nu].l(0), m2=aos[nu].l(1), n2=aos[nu].l(2);
            S(mu,nu) = contracted_gaussian_overlap( // overlap of contracted gto
                aos[mu].center, aos[mu].exponents, aos[mu].coefficients, l1,m1,n1,
                aos[nu].center, aos[nu].exponents, aos[nu].coefficients, l2,m2,n2
            );
        }
    }
    overlap_matrix = std::move(S);
}


// ∂S_uv/∂(center_u) along one Cartesian axis (u and v are AOs).
// χ = N(l,m,n) g_unnorm; ∂g/∂A = 2α g(l+1) − l g(l−1) on unnormalized g, N fixed in α.
// ∫ (∂χ_u/∂A) χ_v = 2α (N_l/N_{l+1}) S(l+1) − l (N_l/N_{l−1}) S(l−1) with normalized overlaps S.
static double contracted_overlap_d_u_v(const AOInfo& u, const AOInfo& v, int dir)
{
    int l1 = u.l(0), m1 = u.l(1), n1 = u.l(2);
    int l2 = v.l(0), m2 = v.l(1), n2 = v.l(2);
    double sum = 0.0;

    for (size_t i = 0; i < u.exponents.n_elem; ++i) {  //iterate through each primitive gaussian function for u anf v
        double alpha = u.exponents(i);
        for (size_t j = 0; j < v.exponents.n_elem; ++j) {
            double beta = v.exponents(j);
            double cij = u.coefficients(i) * v.coefficients(j);

            if (dir == 0) { //x
                double Nl = primitive_gaussian_normalization(alpha, l1, m1, n1);
                double tminus = 0.0;
                if (l1 > 0) {
                    double Nlm = primitive_gaussian_normalization(alpha, l1 - 1, m1, n1);
                    tminus = -static_cast<double>(l1) * (Nl / Nlm) * primitive_gaussian_overlap(
                        alpha, u.center, l1 - 1, m1, n1,
                        beta, v.center, l2, m2, n2);
                }
                double Nlp = primitive_gaussian_normalization(alpha, l1 + 1, m1, n1);
                double tplus = 2.0 * alpha * (Nl / Nlp) * primitive_gaussian_overlap(
                    alpha, u.center, l1 + 1, m1, n1,
                    beta, v.center, l2, m2, n2);
                sum += cij * (tminus + tplus);
            }    
            else if (dir == 1) { //y
                double Nl = primitive_gaussian_normalization(alpha, l1, m1, n1);
                double tminus = 0.0;
                if (m1 > 0) {
                    double Nmm = primitive_gaussian_normalization(alpha, l1, m1 - 1, n1);
                    tminus = -static_cast<double>(m1) * (Nl / Nmm) * primitive_gaussian_overlap(
                        alpha, u.center, l1, m1 - 1, n1,
                        beta, v.center, l2, m2, n2);
                }
                double Nmp = primitive_gaussian_normalization(alpha, l1, m1 + 1, n1);
                double tplus = 2.0 * alpha * (Nl / Nmp) * primitive_gaussian_overlap(
                    alpha, u.center, l1, m1 + 1, n1,
                    beta, v.center, l2, m2, n2);
                sum += cij * (tminus + tplus);
            } 
            else { //z
                double Nl = primitive_gaussian_normalization(alpha, l1, m1, n1);
                double tminus = 0.0;
                if (n1 > 0) {
                    double Nnm = primitive_gaussian_normalization(alpha, l1, m1, n1 - 1);
                    tminus = -static_cast<double>(n1) * (Nl / Nnm) * primitive_gaussian_overlap(
                        alpha, u.center, l1, m1, n1 - 1,
                        beta, v.center, l2, m2, n2);
                }
                double Nnp = primitive_gaussian_normalization(alpha, l1, m1, n1 + 1);
                double tplus = 2.0 * alpha * (Nl / Nnp) * primitive_gaussian_overlap(
                    alpha, u.center, l1, m1, n1 + 1,
                    beta, v.center, l2, m2, n2);
                sum += cij * (tminus + tplus);
            }
        }
    }
    return sum;
}



arma::mat Molecule::find_SR_A()
{
    const int N = static_cast<int>(aos.size());
    arma::mat SR_A(3, N * N, arma::fill::zeros);

    for (int u = 0; u < N; ++u) { //iterate through each AO
        for (int v = 0; v < N; ++v) {

            const int col = u * N + v;
            if (aos[u].atom_idx == aos[v].atom_idx) { //skip if the same atom (0)
                continue;
            }
            for (int direction = 0; direction < 3; ++direction) { //iterate through each direction
                SR_A(direction, col) =
                    contracted_overlap_d_u_v(aos[u], aos[v], direction);
            }
        }
    }

    return SR_A;
}


