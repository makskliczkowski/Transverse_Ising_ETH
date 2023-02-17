#pragma once

#ifndef _XXZ
#define _XXZ

#include "_base.hpp"

/// @brief Model for XXZ chain with next-nearest neighbour terms and disorder
/// @tparam U1_sector U(1) symmetry sector
template <int U1_sector = 0>
class XXZ : 
    public hamiltonian_base<double, U1_hilbert_space<U1::spin>>
{
    //<! ----------------------------------------------------- INHERIT TYPEDEFs FROM BASE
    typedef typename hamiltonian_base<double, U1_hilbert_space<U1::spin>>::matrix matrix;
    typedef typename hamiltonian_base<double, U1_hilbert_space<U1::spin>>::sparse_matrix sparse_matrix;

    //<! ----------------------------------------------------- MODEL PARAMETERS
private:
    disorder<double> disorder_generator;    // generator for random disorder and couplings
    arma::vec _disorder;                    // disorder array on Z field
    double _w = 1.0;                        // disorder strength on top of uniform field
    u64 seed;                               // seed for random generator

protected:
    double _J1 = 1.0;                       // nearest neighbour coupling amplitude
    double _J2 = 0.0;                       // next-nearest neighbour coupling amplitude
    double _delta1 = 0.55;                  // nearest neighbour interaction amplitude
    double _delta2 = 0.0;                   // next-nearest neighbour interaction amplitude

    float total_Sz;                         // total magnetization symmetry sector
    
    //<! ----------------------------------------------------- INITIALIZE MODEL
    virtual void init() override
    {   
        // initialize hilbert space
        this->total_Sz = U1_sector + (this->system_size % 2) / 2.;
        this->_hilbert_space = U1_hilbert_space<U1::spin>(this->system_size, this->total_Sz);
        this->dim = this->_hilbert_space.get_hilbert_space_size();

        // initialize disorder
        disorder_generator = disorder<double>(this->seed);

        // create hamiltonian
        this->create_hamiltonian();
    }
public:
    //<! ----------------------------------------------------- CONSTRUCTORS
    XXZ() = default;
    XXZ(std::istream& os);
    XXZ(int _BC, int L, double J1, double delta1, 
            double w, const u64 _seed = std::random_device{}(), 
            double J2 = 0, double delta2 = 0);

    //<! ----------------------------------------------------- HAMILTONIAN BUILDERS
    virtual void create_hamiltonian() override;
    virtual sparse_matrix create_local_hamiltonian(int site) override;
    virtual void set_hamiltonian_elements(u64 k, double value, u64 new_idx) override;

    //<! ----------------------------------------------------- OVERRIDEN OPERATORS
    virtual std::ostream& write(std::ostream&) const override;
    virtual std::istream& read(std::istream&) override;
};

//<! ---------------------------------------------------------------------------------------------------------------------------------------
//<! ------------------------------------------------------------------------------------------------------------------------ IMPLEMENTATION

//<! ------------------------------------------------------------------------------ CONSTRUCTORS
/// @brief Constructor of XXZ hamiltonian class
/// @tparam U1_sector U(1) symmetry sector as teamplate input
/// @param _BC boundary conditions (=0 for PBC)
/// @param L system size
/// @param J1 nearest neighbour coupling
/// @param delta1 nearest neighbour interaction
/// @param w on-site disorder
/// @param _seed seed for random device
/// @param J2 next-nearest neighbour coupling
/// @param delta2 next-nearest neighbour interaction
template <int U1_sector>
XXZ<U1_sector>::XXZ(int _BC, int L, double J1, double delta1, double w, 
                        const u64 _seed, double J2, double delta2)
{ 
    this->_boundary_condition = _BC;
    this->system_size = L; 
    this->_J1 = J1;
    this->_delta1 = delta1;
    
    //<! disorder terms
    this->_w = w;
    this->seed = _seed;

    //<! NNN terms
    this->_J2 = J2;
    this->_delta2 = delta2;
    init(); 
}


/// @brief Constructor of XXZ hamiltonian class
/// @tparam U1_sector U(1) symmetry sector as teamplate input
/// @param os Input stream to read model parameters
template <int U1_sector>
XXZ<U1_sector>::XXZ(std::istream& os)
    { os >> *this; }

//<! ------------------------------------------------------------------------------ HAMILTONIAN BUILDERS
/// @brief Set hamiltonian matrix element given with value and new index
/// @tparam U1_sector U(1) symmetry sector as teamplate input 
/// @param k current basis state
/// @param value value of matrix element
/// @param new_idx new index to be found in hilbert space
template <int U1_sector>
void XXZ<U1_sector>::set_hamiltonian_elements(u64 k, double value, u64 new_idx)
{
    u64 idx = this->_hilbert_space.find(new_idx);
    try {
        H(idx, k) += value;
        H(k, idx) += value;
    } 
    catch (const std::exception& err) {
        std::cout << "Exception:\t" << err.what() << "\n";
        std::cout << "SHit ehhh..." << std::endl;
        printSeparated(std::cout, "\t", 14, true, new_idx, idx, this->_hilbert_space(k), value);
    }
}


/// @brief Method to create hamiltonian within the class
/// @tparam U1_sector U(1) symmetry sector as teamplate input
template <int U1_sector>
void XXZ<U1_sector>::create_hamiltonian()
{
    this->H = sparse_matrix(this->dim, this->dim);
    this->_disorder = disorder_generator.uniform(system_size, this->_w);    
	std::cout << "disorder: \t\t" << this->_disorder.t() << std::endl;

    std::vector<double> coupling = {this->_J1, this->_J2};
    std::vector<double> interaction = {this->_delta1, this->_delta2};
    
    std::vector<int> neighbor_distance = {1, 2};
    for (u64 k = 0; k < this->dim; k++) {
		double s_i, s_j;
		u64 base_state = this->_hilbert_space(k);
		for (int j = 0; j < this->system_size; j++) {
			s_i = checkBit(base_state, this->system_size - 1 - j) ? 0.5 : -0.5;				// true - spin up, false - spin down
			
            //<! longitudinal field with disorder
			this->H(k, k) += this->_disorder(j) * s_i;                            // diagonal elements setting

			for(int a = 0; a < neighbor_distance.size(); a++){
                int r = neighbor_distance[a];
                int nei = j + r;
                if(nei >= this->system_size)
                    nei = (this->_boundary_condition)? -1 : nei % this->system_size;

                
                if (nei >= 0) //<! boundary conditions
                {
                    s_j = checkBit(base_state, this->system_size - 1 - nei) ? 0.5 : -0.5;
                    if(s_i < 0 && s_j > 0){
                        u64 new_idx =  flip(base_state, BinaryPowers[this->system_size - 1 - nei], this->system_size - 1 - nei);
                        new_idx =  flip(new_idx, BinaryPowers[this->system_size - 1 - j], this->system_size - 1 - j);
                        
                        // 0.5 cause flip 0.5*(S+S- + S-S+)
                        this->set_hamiltonian_elements(k, 0.5 * coupling[a], new_idx);
                    }
                    
                    //<! Interaction (spin correlations) with neighbour at distance r
                    this->H(k, k) += interaction[a] * s_i * s_j;
                }
            }
		}
		//std::cout << std::bitset<4>(base_state) << "\t";
	}
}


/// @brief Method to create hamiltonian within the class
/// @tparam U1_sector U(1) symmetry sector as teamplate input
/// @param site site index where the local hamiltonian acts
/// @return the local hamiltonian at site site
template <int U1_sector>
typename XXZ<U1_sector>::sparse_matrix XXZ<U1_sector>::create_local_hamiltonian(int site)
{
    sparse_matrix H_local(dim, dim);
    v_2d<double> parameters = { {this->_J1, this->_J1, this->_delta1},
                                {this->_J2, this->_J2, this->_delta2}
                                };
    
    std::vector<op_type> _operators = {operators::sigma_x, operators::sigma_y, operators::sigma_z };
    std::vector<int> neighbor_distance = {1, 2};

    for (size_t k = 0; k < this->dim; k++) 
    {
        int base_state = this->_hilbert_space(k);
        //<! go over all neighbours
        for(int a = 0; a < neighbor_distance.size(); a++)
        {
            int r = neighbor_distance[a];
            int nei = site + r;
            if(nei >= this->system_size)
                nei = (this->_boundary_condition)? -1 : nei % this->system_size;

            if (nei >= 0) {
                //<! coupling and interaction
                for(int b = 0; b < _operators.size(); b++){
                    op_type op = _operators[b];
                    auto [val1, op_k] = op(base_state, this->system_size, { site });
                    auto [val2, opop_k] = op(op_k, this->system_size, { nei });
                    
                    u64 idx = this->_hilbert_space.find(opop_k);
                    H_local(idx, k) += parameters[a][b] * real(val1 * val2);
                }

                if(a < 1)
                {
                    //<! longitudinal field with disorder
                    auto [val, op_k] = operators::sigma_z(base_state, this->system_size, { site });
                    u64 idx = this->_hilbert_space.find(op_k);
                    H_local(idx, k) += this->_disorder(site) * real(val);

                    auto [val2, op_k2] = operators::sigma_z(base_state, this->system_size, { nei });
                    idx = this->_hilbert_space.find(op_k2);
                    H_local(idx, k) += this->_disorder(nei) * real(val2);
                }
            }
        }
    }

    return H_local;
}


//<! ------------------------------------------------------------------------------ OVVERRIDEN OPERATORS AND OPERATOR KERNELS
/// @brief Read model parameters from input stream
/// @tparam U1_sector U(1) symmetry sector as teamplate input 
/// @param os input stream to read parameters
template <int U1_sector>
std::istream& XXZ<U1_sector>::read(std::istream& os)
{
    
    return os;
}

/// @brief Write hamiltonian to stream as human readable
/// @tparam U1_sector U(1) symmetry sector as teamplate input 
/// @param os input stream to read parameters
template <int U1_sector>
std::ostream& XXZ<U1_sector>::write(std::ostream& os) const
{
    printSeparated(os, "\t", 16, true, "Model:", "XXZ - 1D Anisotropic Heisenberg model");
    os << std::endl;
    printSeparated(os, "\t", 16, true, "Hamiltonian:", "H = \u03A3_r J_r\u03A3_i [(S^x_i S^x_i+1 + S^y_i S^y_i+1) + \u0394_r S^z_iS^z_i+1] + \u03A3_i h_i S^z_i");
    printSeparated(os, "\t", 16, true, "----------------------------------------------------------------------------------------------------");
    printSeparated(os, "\t", 16, true, "Parameters:");
    printSeparated(os, "\t", 16, true, "    L", this->system_size);
    printSeparated(os, "\t", 16, true, "    J_1", this->_J1);
    printSeparated(os, "\t", 16, true, "    \u0394_1", this->_delta1);
    printSeparated(os, "\t", 16, true, "    w", this->_w);
    //printSeparated(os, "\t", 16, true, "disorder", this->_disorder.t());

    printSeparated(os, "\t", 16, true, "    J_2", this->_J2);
    printSeparated(os, "\t", 16, true, "    \u0394_2", this->_delta2);
    printSeparated(os, "\t", 16, true, "    seed", this->seed);
    printSeparated(os, "\t", 16, true, "----------------------------------------------------------------------------------------------------");
    printSeparated(os, "\t", 16, true, "----------------------------------------------------------------------------------------------------");

    return os;
}






#endif