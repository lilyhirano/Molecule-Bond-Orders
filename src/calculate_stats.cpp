#define ARMA_USE_HDF5
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <armadillo>
#include <nlohmann/json.hpp>
#include <filesystem>

using json = nlohmann::json;
namespace fs = std::filesystem;

struct BondStats
{
    double mae;
    double stdev;
};

BondStats calculate_mae_stdev(const arma::mat& mat1, const arma::mat& mat2)
{
    std::vector<double> diffs;
    double total_diff = 0.0;
    int n = mat1.n_rows;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i == j) continue;
            double d = std::abs(mat1(i, j) - mat2(i, j));
            diffs.push_back(d);
            total_diff += d;
        }
    }

    double mae = total_diff / diffs.size();
    double var_sum = 0.0;
    for (double d : diffs)
    {
        var_sum += std::pow(d - mae, 2);
    }
    
    return {mae, std::sqrt(var_sum / diffs.size())};
}

void printStatRow(std::string m1, std::string m2, BondStats s)
{
    std::cout.precision(6);
    std::cout.setf(std::ios::fixed); 

    std::cout << m1 << " vs " << m2 << "\t | MAE: " << s.mae << " | STDEV: " << s.stdev << std::endl;
}



int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " path/to/config" << std::endl;
        return EXIT_FAILURE;
    }

    // parse config
    fs::path config_file_path(argv[1]);
    if (!fs::exists(config_file_path)) {
        std::cerr << "Path: " << config_file_path << " does not exist" << std::endl;
        return EXIT_FAILURE;
    }

    std::ifstream config_file(config_file_path);
    json config = json::parse(config_file);

    // get output path for HDF5 file
    fs::path hdf5_path = config["output_file_path"];

    if (!fs::exists(hdf5_path)) {
        std::cerr << "HDF5 file not found at: " << hdf5_path << std::endl;
        return EXIT_FAILURE;
    }

    std::string mol_name = fs::path(config["atoms_file_path"]).stem().string();
    std::string csv_path = "tests/test_checks/" + mol_name + ".csv";

    arma::mat wiberg, mayer, mulliken, reference;

    // load matrices from HDF5 file
    try
    {
        wiberg.load(arma::hdf5_name(hdf5_path.string(), "wiberg_bond_order_matrix"));
        mayer.load(arma::hdf5_name(hdf5_path.string(), "mayer_bond_order_matrix"));
        mulliken.load(arma::hdf5_name(hdf5_path.string(), "mulliken_bond_order_matrix"));

        if (!reference.load(csv_path, arma::csv_ascii))
        {
            csv_path = "../tests/test_checks/" + mol_name + ".csv";
            reference.load(csv_path, arma::csv_ascii);
        }
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error loading HDF5 file: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "\n=== STATS: REFERENCE vs. METHODS (" << mol_name << ") ===" << std::endl;
    
    // compare each method to reference
    printStatRow("Reference", "Wiberg",   calculate_mae_stdev(reference, wiberg));
    printStatRow("Reference", "Mayer",    calculate_mae_stdev(reference, mayer));
    printStatRow("Reference", "Mulliken", calculate_mae_stdev(reference, mulliken));
    
    std::cout << "=================================================\n" << std::endl;
    return EXIT_SUCCESS;
}