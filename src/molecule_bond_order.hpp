#include <fstream>
#include <iostream>
#include <stdexcept>
#include <armadillo>
#include <vector>
#include <unordered_map>
#include <highfive/H5File.hpp>
#include <nlohmann/json.hpp>



using json = nlohmann::json; 
namespace fs = std::filesystem;
using namespace std;

struct Atom
{
    int element_symbol;
    arma::vec center;
    // [0]=(I_s+A_s)/2, [1]=(I_p+A_p)/2, [2]=|−β| from Table 1.1 (eV), [3]=Z*
    vector<double> parameters;
};

struct AOInfo {
    int atom_idx;
    string type; 
    arma::vec center;
    arma::vec exponents;
    arma::vec coefficients;
    arma::ivec l; // (lx,ly,lz)

};

// used by overlap  and gamma integrals
double primitive_gaussian_normalization(double alpha, int l, int m, int n);

class STO3G_maps;

class Molecule{
    public:
        int num_atoms;
        int num_basis_functions;
        vector<int> num_element_symbols; // [H, C, N, O, F]
        string atoms_file_path_; 
        vector<Atom> atoms;
        int p;
        int q;
        vector<AOInfo> aos;
        arma::mat overlap_matrix;
    
        Molecule(const string& atoms_file_path, int p, int q): atoms_file_path_(atoms_file_path), p(p), q(q){
            read_atoms_file();
            determine_num_basis_functions();
        }

        //copy constructor
        Molecule(const Molecule& other): atoms_file_path_(other.atoms_file_path_), p(other.p), q(other.q), 
        num_atoms(other.num_atoms), num_basis_functions(other.num_basis_functions), 
        num_element_symbols(other.num_element_symbols), atoms(other.atoms),
         aos(other.aos), overlap_matrix(other.overlap_matrix){}

        ~Molecule(){}

        void build_aos(const STO3G_maps& sto3g);
        void build_overlap_matrix();
    
    private:
        //implemented in molecule_functions.cpp
        void read_atoms_file(); 
        void determine_num_basis_functions(); 
    };

class STO3G_maps{
    public:
        unordered_map<string, arma::vec> H_s;
        unordered_map<string, arma::vec> C_p;
        unordered_map<string, arma::vec> C_s;
        unordered_map<string, arma::vec> N_p;
        unordered_map<string, arma::vec> N_s;
        unordered_map<string, arma::vec> F_p;
        unordered_map<string, arma::vec> F_s;
        unordered_map<string, arma::vec> O_p;
        unordered_map<string, arma::vec> O_s;
    
        STO3G_maps(){
            read_STO3G_basis_functions("./basis/H_s_STO3G.json", H_s);
            read_STO3G_basis_functions("./basis/C_s_STO3G.json", C_s);
            read_STO3G_basis_functions("./basis/C_p_STO3G.json", C_p);
            read_STO3G_basis_functions("./basis/N_s_STO3G.json", N_s);
            read_STO3G_basis_functions("./basis/N_p_STO3G.json", N_p);
            read_STO3G_basis_functions("./basis/F_s_STO3G.json", F_s);
            read_STO3G_basis_functions("./basis/F_p_STO3G.json", F_p);
            read_STO3G_basis_functions("./basis/O_s_STO3G.json", O_s);
            read_STO3G_basis_functions("./basis/O_p_STO3G.json", O_p);
        }
        ~STO3G_maps(){}
    
    private:
        //implemented in sto3g_functions.cpp
        void read_STO3G_basis_functions(const string &file_path, unordered_map<string, arma::vec> &basis_function_map);
    };



class FockMatrix {
    public:
        explicit FockMatrix(const Molecule& mol);

        ~FockMatrix() = default;

        arma::mat get_density_matrix_alpha() const { return density_matrix_alpha; }
        arma::mat get_density_matrix_beta() const { return density_matrix_beta; }
        
        // calcualte for both alpha and beta
        void calculate_fock_spin(){
            assemble_fock_spin(density_matrix_alpha, fock_matrix_alpha);
            assemble_fock_spin(density_matrix_beta, fock_matrix_beta);
            iterations++;
        }

        void solve_eigenvalue_problem();
        void find_convergence();
    private:
        const Molecule& molecule;
        int p;
        int q;
        int iterations = -1;

        arma::mat density_matrix_alpha;
        arma::mat density_matrix_beta;
        arma::mat fock_matrix_alpha;
        arma::mat fock_matrix_beta;
        arma::mat prev_density_matrix_alpha;
        arma::mat prev_density_matrix_beta;


        arma::vec eigenvalues_alpha;
        arma::vec eigenvalues_beta;
        int eigen = 0;

        void calculate_gamma_matrix();

        void assemble_fock_spin(const arma::mat& P_spin, arma::mat& F_out);
        bool check_convergence();
};


void fock_diagonal_elements(const Molecule& mol, const arma::mat& gamma_atom, const arma::mat& P_tot,
                            const arma::mat& P_spin, arma::mat& F_spin);
void fock_off_diagonal_elements(const Molecule& mol, const arma::mat& gamma_atom, const arma::mat& overlap,
                                const arma::mat& P_spin, arma::mat& F_spin);


class BondOrder{
    public:
         BondOrder(const Molecule& mol) : molecule(mol)
         {
            FockMatrix fock_matrix(molecule);
            fock_matrix.find_convergence(); //read converged density matrix
            density_matrix_alpha = fock_matrix.get_density_matrix_alpha();
            density_matrix_beta = fock_matrix.get_density_matrix_beta();
            density_matrix_total = density_matrix_alpha + density_matrix_beta;

            int num_atoms = molecule.num_atoms;
            mulliken_bond_order_matrix = arma::zeros(num_atoms, num_atoms);
            mayer_bond_order_matrix = arma::zeros(num_atoms, num_atoms);
            wiberg_bond_order_matrix = arma::zeros(num_atoms, num_atoms);
            X_matrix = arma::zeros(molecule.num_basis_functions, molecule.num_basis_functions);

            calculate_X_matrix();
            calculate_wiberg_bond_order(); 
            calculate_mayer_bond_order();
            calculate_mulliken_bond_order();
         };
        ~BondOrder() = default;

        arma::mat get_mulliken_bond_order_matrix() const { return mulliken_bond_order_matrix; }
        arma::mat get_mayer_bond_order_matrix() const { return mayer_bond_order_matrix; }
        arma::mat get_wiberg_bond_order_matrix() const { return wiberg_bond_order_matrix; }
        arma::mat get_density_matrix_total() const { return density_matrix_total; }
        arma::mat get_X_matrix() const { return X_matrix; }

    private:
        const Molecule& molecule;
        arma::mat density_matrix_alpha;
        arma::mat density_matrix_beta;
        arma::mat density_matrix_total;
        arma::mat mulliken_bond_order_matrix;
        arma::mat mayer_bond_order_matrix;
        arma::mat wiberg_bond_order_matrix;
        arma::mat X_matrix;

        void calculate_wiberg_bond_order();
        void calculate_mayer_bond_order();
        void calculate_mulliken_bond_order();
        void calculate_X_matrix();

};
