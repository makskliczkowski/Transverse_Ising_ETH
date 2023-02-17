#pragma once

#ifndef _QUANTUM_SUN
#define _QUANTUM_SUN

#include "_base.hpp"

/// @brief Model for EBT, QuantumSun model
class QuantumSun : 
    public hamiltonian_base<double, full_hilbert_space>
{
    //<! ----------------------------------------------------- INHERIT TYPEDEFs FROM BASE
    typedef typename hamiltonian_base<double, full_hilbert_space>::matrix        matrix;
    typedef typename hamiltonian_base<double, full_hilbert_space>::sparse_matrix sparse_matrix;

    //<! ----------------------------------------------------- MODEL PARAMETERS
private:
    disorder<double> disorder_generator;    // generator for random disorder and couplings
    
    arma::vec _long_range_couplings;        // random coupling, i.e. distance of spins to grain
    arma::vec _disorder;                    // disorder array on Z field
    
    double _w = 0.5;                        // disorder value on top of uniform field
    double _J = 1.;                         // coupling amplitude
    double _hz = 1.0;                       // longitudinal uniform field
    double _alfa = 0.75;                    // coupling base -- controls long-range interaction
    double _zeta = 0.2;                     // randomness in long range coupling (random distance between localised spins abd grain)
    u64 _seed = std::random_device{}();     // seed for random generator
    
    int grain_size = 3;                     // grain size
    bool start_counting = 1;                // start first coupling =1 ? (i.e. exponent u_0 = 0)
    
    //<! ----------------------------------------------------- INITIALIZE MODEL
    virtual void init() override
    {   
        // initialize hilbert space
        this->_hilbert_space = full_hilbert_space(this->system_size);
        this->dim = this->_hilbert_space.get_hilbert_space_size();

        // initialize disorder
        disorder_generator = disorder<double>(this->_seed);

        // create hamiltonian
        this->create_hamiltonian();
    }

public:
    //<! ----------------------------------------------------- CONSTRUCTORS
    QuantumSun() = default;
    QuantumSun(std::istream& os);
    QuantumSun(int L, double J, double alfa, 
            double w, double hz, const u64 seed = std::random_device{}(), int M = 3, double zeta = 0.2);

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

QuantumSun::QuantumSun(int L, double J, double alfa, 
            double w, double hz, const u64 seed, int M, double zeta)
{ 
    this->system_size = L; 
    this->grain_size = M;
    this->_J = J;
    this->_alfa = alfa;
    this->_zeta = zeta;
    
    this->_hz = hz;
    //<! disorder terms
    this->_w = w;
    this->_seed = seed;
    init(); 
}

QuantumSun::QuantumSun(std::istream& os)
    { os >> *this; }

//<! ------------------------------------------------------------------------------ HAMILTONIAN BUILDERS
/// @brief Set hamiltonian matrix element given with value and new index
/// @tparam U1_sector U(1) symmetry sector as teamplate input 
/// @param k current basis state
/// @param value value of matrix element
/// @param new_idx new index to be found in hilbert space
void QuantumSun::set_hamiltonian_elements(u64 k, double value, u64 new_idx)
{
    u64 idx = this->_hilbert_space.find(new_idx);
    try {
        H(idx, k) += value;
    } 
    catch (const std::exception& err) {
        std::cout << "Exception:\t" << err.what() << "\n";
        std::cout << "SHit ehhh..." << std::endl;
        printSeparated(std::cout, "\t", 14, true, new_idx, idx, this->_hilbert_space(k), value);
    }
}


/// @brief Method to create hamiltonian within the class
void QuantumSun::create_hamiltonian()
{
    this->H = sparse_matrix(this->dim, this->dim);
    this->_disorder = disorder_generator.uniform(system_size, this->_hz - this->_w, this->_hz + this->_w); 
    
    const size_t dim_loc = ULLPOW( (this->system_size - this->grain_size) );
	const size_t dim_erg = ULLPOW( (this->grain_size) );
	
    auto neighbor_generator = disorder<int>(this->_seed);

    /* Create random neighbours for coupling hamiltonian */
    auto random_neigh = neighbor_generator.uniform(this->system_size, 0, this->grain_size - 1);
    
	#ifndef ENSEMBLE
		#define ENSEMBLE GOE
		#pragma message ("--> Using implicit random matrix ensemble: i.e., Gaussian Orthogonal Ensemble")
	#endif

	/* Create GOE Matrix */
	auto grain = ENSEMBLE(my_disorder.get_seed());
	arma::mat H_grain = grain.generate_matrix(dim_erg);

    
    /* Create random couplings */
    this->_long_range_couplings = arma::vec(this->system_size, arma::fill::zeros);
    if(this->_alfa > 0){
        if(this->_alfa < 1.0){
            for (int j = this->grain_size; j < this->system_size; j++){ 
                int pos = j - this->grain_size + 1;
                double u_j = pos + disorder_generator.random_uni<double>(-this->_zeta, this->_zeta);
                this->_long_range_couplings(j) = std::pow(this->_alfa, u_j);
            }
        } else {
            this->_long_range_couplings = arma::vec(this->system_size, arma::fill::ones);
        }
    }
    #ifdef EXTRA_DEBUG
	    std::cout << "disorder: \t\t" << this->_disorder.t() << std::endl;   
	    std::cout << "couplings: \t\t" << this->_long_range_couplings.t() << std::endl;
	    std::cout << "random_neigh: \t\t" << random_neigh.t() << std::endl;
        std::cout << "Grain matrix: \t\t" << H_grain << std::endl;
    #endif

    /* Generate coupling and spin hamiltonian */
    for (u64 k = 0; k < this->dim; k++) {
		u64 base_state = this->_hilbert_space(k);
		for (int j = this->grain_size; j < this->system_size; j++) {

			/* disorder on localised spins */
            auto [val, Sz_k] = operators::sigma_z(base_state, this->system_size, { j });
			this->set_hamiltonian_elements(k, this->_disorder(j) * real(val), Sz_k);

			/* coupling of localised spins to GOE grain */
			int nei = random_neigh(j);
		    auto [val1, Sx_k] = operators::sigma_x(base_state, this->system_size, { j });
		    auto [val2, SxSx_k] = operators::sigma_x(Sx_k, this->system_size, { nei });
			this->set_hamiltonian_elements(k, this->_J * this->_long_range_couplings(j) * real(val1 * val2), SxSx_k);
		}
	}
	this->H = this->H + arma::kron(H_grain, arma::eye(dim_loc, dim_loc));
}


/// @brief Method to create hamiltonian within the class
/// @param site site index where the local hamiltonian acts
/// @return the local hamiltonian at site site
typename QuantumSun::sparse_matrix QuantumSun::create_local_hamiltonian(int site)
{
    sparse_matrix H_local(dim, dim);
    
    return H_local;
}


//<! ------------------------------------------------------------------------------ OVVERRIDEN OPERATORS AND OPERATOR KERNELS
/// @brief Read model parameters from input stream
/// @tparam U1_sector U(1) symmetry sector as teamplate input 
/// @param os input stream to read parameters
std::istream& QuantumSun::read(std::istream& os)
{
    
    return os;
}

/// @brief Write hamiltonian to stream as human readable
/// @tparam U1_sector U(1) symmetry sector as teamplate input 
/// @param os input stream to read parameters
std::ostream& QuantumSun::write(std::ostream& os) const
{
    printSeparated(os, "\t", 16, true, "Model:", "Quantum Sun model - O-dimensional EBT toy model");
    os << std::endl;
    printSeparated(os, "\t", 16, true, "Hamiltonian:", "H = R + J\u03A3_i \u03B1^{u_j} S^x_i S^x_i+1 + \u03A3_i h_i S^z_i\t\t u_j in [j - \u03B6, j + \u03B6]");
    printSeparated(os, "\t", 16, true, "----------------------------------------------------------------------------------------------------");
    printSeparated(os, "\t", 16, true, "Parameters:");
    printSeparated(os, "\t", 16, true, "    L", this->system_size);
    printSeparated(os, "\t", 16, true, "    grain size", this->grain_size);

    printSeparated(os, "\t", 16, true, "    J", this->_J);
    printSeparated(os, "\t", 16, true, "    \u03B1", this->_alfa);
    printSeparated(os, "\t", 16, true, "    w", this->_w);
    printSeparated(os, "\t", 16, true, "    hz", this->_hz);
    printSeparated(os, "\t", 16, true, "    \u03B6", this->_zeta);
    //printSeparated(os, "\t", 16, true, "disorder", this->_disorder.t());

    printSeparated(os, "\t", 16, true, "    seed", this->_seed);
    
    return os;
}





#endif