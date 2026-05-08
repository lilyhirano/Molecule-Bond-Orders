#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <armadillo>
#include <nlohmann/json.hpp>
#include "molecule_bond_order.hpp"

using namespace std;
using json = nlohmann::json; 

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " path/to/config" << std::endl;
        return EXIT_FAILURE;
    }

    // parse the config file
    fs::path config_file_path(argv[1]);
    if (!fs::exists(config_file_path)) {
        std::cerr << "Path: " << config_file_path << " does not exist" << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream config_file(config_file_path);
    json config = json::parse(config_file);

    fs::path atoms_file_path = config["atoms_file_path"];
    fs::path output_file_path = config["output_file_path"];
    int p = config["num_alpha_electrons"];
    int q = config["num_beta_electrons"];
        
    STO3G_maps sto3g;
    Molecule molecule(atoms_file_path, p, q);
    molecule.build_aos(sto3g);
    molecule.build_overlap_matrix();

    BondOrder bond_order(molecule);
    
    std::cout << std::defaultfloat << std::setprecision(6); 

    cout << "Overlap matrix:" << endl;
    cout << molecule.overlap_matrix << endl;
    cout << "Orthogonal overlap matrix:" << endl;
    cout<< bond_order.get_X_matrix() << endl;
    cout << "Total density matrix:" << endl;
    cout<< bond_order.get_density_matrix_total() << endl;
    
    cout<<"Wiberg bond order matrix:" << endl;
    cout<<bond_order.get_wiberg_bond_order_matrix() << endl;
    cout<<"Mayer bond order matrix:" << endl;
    cout<<bond_order.get_mayer_bond_order_matrix() << endl;
    cout<<"Mulliken bond order matrix:" << endl;
    cout<<bond_order.get_mulliken_bond_order_matrix() << endl;



      // check that output dir exists
    if (!fs::exists(output_file_path.parent_path())){
        fs::create_directories(output_file_path.parent_path()); 
    }

    // delete the file if it does exist (so that no old answers stay there by accident)
        if (fs::exists(output_file_path)){
        fs::remove(output_file_path); 
    }

    HighFive::File file(output_file_path.string(), HighFive::File::Overwrite);

    bond_order.get_wiberg_bond_order_matrix().save(arma::hdf5_name(output_file_path, "wiberg_bond_order_matrix",
        arma::hdf5_opts::append + arma::hdf5_opts::trans));
    bond_order.get_mayer_bond_order_matrix().save(arma::hdf5_name(output_file_path, "mayer_bond_order_matrix",
        arma::hdf5_opts::append + arma::hdf5_opts::trans));
    bond_order.get_mulliken_bond_order_matrix().save(arma::hdf5_name(output_file_path, "mulliken_bond_order_matrix",
        arma::hdf5_opts::append + arma::hdf5_opts::trans));

    return EXIT_SUCCESS;
}
