#include <fstream>
#include <armadillo>
#include <vector>
#include <unordered_map>
#include <nlohmann/json.hpp>
#include "molecule_visualization.hpp"

using namespace std;
using json = nlohmann::json; 

/*
    Reads the STO3G basis functions from a JSON file and stores them in a map.
    The map's form is :
        atomic_number: [atomic_number],
        shell_momentum: [shell_momentum],
        exponents: [exponents],
        contraction_coefficients: [contraction_coefficients]
*/
void STO3G_maps::read_STO3G_basis_functions(
    const string &file_path, unordered_map<string, arma::vec> &basis_function_map)
{
    ifstream basis_function_file(file_path);
    json basis_function_json = json::parse(basis_function_file);

    basis_function_map["atomic_number"] =
        arma::vec({static_cast<double>(basis_function_json["atomic_number"])});
    basis_function_map["shell_momentum"] =
        arma::vec({static_cast<double>(basis_function_json["shell_momentum"])});

    const auto& cgs = basis_function_json["contracted_gaussians"];
    arma::vec exps(cgs.size(), arma::fill::zeros);
    arma::vec coeffs(cgs.size(), arma::fill::zeros);
    for (size_t i = 0; i < cgs.size(); i++) {
        exps(i) = static_cast<double>(cgs[i]["exponent"]);
        coeffs(i) = static_cast<double>(cgs[i]["contraction_coefficient"]);
    }
    basis_function_map["exponents"] = std::move(exps);
    basis_function_map["contraction_coefficients"] = std::move(coeffs);
}