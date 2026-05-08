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
    std::string molecule;
    double mae;
};

double calculate_mae(const arma::mat& mat1, const arma::mat& mat2)
{
    double total_diff = 0.0;
    int n = mat1.n_rows;
    int count = 0;

    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        { // Only upper triangle to avoid double-counting
            total_diff += std::abs(mat1(i, j) - mat2(i, j));
            count++;
        }
    }
    return (count > 0) ? (total_diff / count) : 0.0;
}

void printMethodTable(std::string methodName, const std::vector<BondStats>& stats)
{
    std::cout << "\n--- TABLE: " << methodName << " ---" << std::endl;
    std::cout << "Molecule\t| Mean Difference" << std::endl;
    std::cout << "---------------------------------" << std::endl;
    for (const auto& s : stats)
    {
        std::cout << s.molecule << "\t\t| " << s.mae << std::endl;
    }
}


int main() {
    std::vector<std::string> molecules = {"H2", "H2O", "HF", "HO", "HCN", "C2H4"};
    std::vector<BondStats> wiberg_table, mayer_table, mulliken_table;

    for (const std::string& mol : molecules) {
        // Construct paths
        std::string config_path = "sample_input/final_proj/" + mol + ".json";
        
        // Skip if config doesn't exist
        if (!fs::exists(config_path)) continue;

        std::ifstream f(config_path);
        json config = json::parse(f);
        fs::path hdf5_path = config["output_file_path"];
        std::string csv_path = "tests/test_checks/" + mol + ".csv";

        arma::mat wiberg, mayer, mulliken, reference;

        try {
            wiberg.load(arma::hdf5_name(hdf5_path.string(), "wiberg_bond_order_matrix"));
            mayer.load(arma::hdf5_name(hdf5_path.string(), "mayer_bond_order_matrix"));
            mulliken.load(arma::hdf5_name(hdf5_path.string(), "mulliken_bond_order_matrix"));
            
            if (!reference.load(csv_path, arma::csv_ascii)) {
                reference.load("../" + csv_path, arma::csv_ascii);
            }

            wiberg_table.push_back({mol, calculate_mae(wiberg, reference)});
            mayer_table.push_back({mol, calculate_mae(mayer, reference)});
            mulliken_table.push_back({mol, calculate_mae(mulliken, reference)});

        } catch (...) {
            continue; 
        }
    }

    // Output the three separate tables
    printMethodTable("WIBERG", wiberg_table);
    printMethodTable("MAYER", mayer_table);
    printMethodTable("MULLIKEN", mulliken_table);

    return 0;
}

// int main(int argc, char* argv[]) {
//     if (argc != 2) {
//         std::cerr << "Usage: " << argv[0] << " path/to/config" << std::endl;
//         return EXIT_FAILURE;
//     }

//     // parse config
//     fs::path config_file_path(argv[1]);
//     if (!fs::exists(config_file_path)) {
//         std::cerr << "Path: " << config_file_path << " does not exist" << std::endl;
//         return EXIT_FAILURE;
//     }

//     std::ifstream config_file(config_file_path);
//     json config = json::parse(config_file);

//     // get output path for HDF5 file
//     fs::path hdf5_path = config["output_file_path"];

//     if (!fs::exists(hdf5_path)) {
//         std::cerr << "HDF5 file not found at: " << hdf5_path << std::endl;
//         return EXIT_FAILURE;
//     }

//     std::string mol_name = fs::path(config["atoms_file_path"]).stem().string();
//     std::string csv_path = "tests/test_checks/" + mol_name + ".csv";

//     arma::mat wiberg, mayer, mulliken, reference;

//     // load matrices from HDF5 file
//     try
//     {
//         wiberg.load(arma::hdf5_name(hdf5_path.string(), "wiberg_bond_order_matrix"));
//         mayer.load(arma::hdf5_name(hdf5_path.string(), "mayer_bond_order_matrix"));
//         mulliken.load(arma::hdf5_name(hdf5_path.string(), "mulliken_bond_order_matrix"));

//         if (!reference.load(csv_path, arma::csv_ascii))
//         {
//             csv_path = "../tests/test_checks/" + mol_name + ".csv";
//             reference.load(csv_path, arma::csv_ascii);
//         }
//     }
//     catch (const std::exception& e)
//     {
//         std::cerr << "Error loading HDF5 file: " << e.what() << std::endl;
//         return EXIT_FAILURE;
//     }

//     std::cout << "\n=== STATS: REFERENCE vs. METHODS (" << mol_name << ") ===" << std::endl;
    
//     // compare each method to reference
//     printStatRow("Reference", "Wiberg",   calculate_stats(reference, wiberg));
//     printStatRow("Reference", "Mayer",    calculate_stats(reference, mayer));
//     printStatRow("Reference", "Mulliken", calculate_stats(reference, mulliken));
    
//     std::cout << "=================================================\n" << std::endl;
//     return EXIT_SUCCESS;
// }