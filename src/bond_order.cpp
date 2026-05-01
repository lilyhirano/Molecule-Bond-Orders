#include "molecule_visualization.hpp"
#include <math.h>
#include <armadillo>

//orthogonalization transformation
void BondOrder::calculate_X_matrix() {
    arma::vec s;
    arma::mat U;
    arma::eig_sym(s, U, molecule.overlap_matrix);

    X_matrix = U * arma::diagmat(arma::pow(s, -0.5)) * U.t();
}


//need to build a map of atoms to their basis functions
std::vector<std::vector<int>> build_atom_to_ao_map(const Molecule& molecule) {
    std::vector<std::vector<int>> atom_to_aos(static_cast<size_t>(molecule.num_atoms));

    for (int ao = 0; ao < static_cast<int>(molecule.aos.size()); ++ao) {
        const int atom_idx = molecule.aos[ao].atom_idx;
        atom_to_aos[static_cast<size_t>(atom_idx)].push_back(ao);
    }
    return atom_to_aos;
}


void BondOrder::calculate_wiberg_bond_order(){
    int num_atoms = molecule.num_atoms;
    const auto atom_to_aos = build_atom_to_ao_map(molecule);

    for(int i = 0; i < num_atoms; ++i){
        for(int j = i + 1; j < num_atoms; ++j){
            if(i == j){
                continue;
            }

            double bond_order = 0.0;
            for (int mu : atom_to_aos[static_cast<size_t>(i)]) {
                for (int nu : atom_to_aos[static_cast<size_t>(j)]) {
                    bond_order += density_matrix_total(mu, nu) * density_matrix_total(nu, mu);
                }
            }
            wiberg_bond_order_matrix(i, j) = bond_order;
            wiberg_bond_order_matrix(j, i) = bond_order;
        }

    }
}


void BondOrder::calculate_mayer_bond_order(){
    int num_atoms = molecule.num_atoms;
    const arma::mat& S = molecule.overlap_matrix;
    const arma::mat PS_alpha = density_matrix_alpha * S;
    const arma::mat PS_beta = density_matrix_beta * S;

    const auto atom_to_aos = build_atom_to_ao_map(molecule);

    for(int i = 0; i < num_atoms; ++i){
        for(int j = 0; j < num_atoms; ++j){
            if(i == j){
                continue;
            }

            double bond_order = 0.0;
            for (int mu : atom_to_aos[static_cast<size_t>(i)]) {
                for (int nu : atom_to_aos[static_cast<size_t>(j)]) {
                    bond_order += PS_alpha(mu, nu) * PS_alpha(nu, mu);
                    bond_order += PS_beta(mu, nu) * PS_beta(nu, mu);
                }
            }
            mayer_bond_order_matrix(i, j) =  bond_order;
        }
    }
}


void BondOrder::calculate_mulliken_bond_order(){
    int num_atoms = molecule.num_atoms;
    const auto atom_to_aos = build_atom_to_ao_map(molecule);
    arma::mat overlap_matrix = molecule.overlap_matrix;

    for(int i = 0; i < num_atoms; ++i){
        for(int j = i + 1; j < num_atoms; ++j){
            if(i == j){
                continue;
            }

            double bond_order = 0.0;
            for (int mu : atom_to_aos[static_cast<size_t>(i)]) {
                for (int nu : atom_to_aos[static_cast<size_t>(j)]) {
                    bond_order += density_matrix_total(mu, nu) * overlap_matrix(mu, nu);
                }
            }
            mulliken_bond_order_matrix(i, j) = 2 * bond_order;
            mulliken_bond_order_matrix(j, i) = 2 * bond_order;
        }
    }
}