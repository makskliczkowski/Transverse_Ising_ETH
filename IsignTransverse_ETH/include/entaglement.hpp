
template <typename _ty>
inline
_ty my_conjungate(_ty x) { return std::conj(x); }

template <>
inline
int my_conjungate(int x) { return x; }

template <>
inline
float my_conjungate(float x) { return x; }

template <>
inline
double my_conjungate(double x) { return x; }

template <typename _type>
inline
std::complex<_type> my_conjungate(std::complex<_type> x) { return std::conj(x); }

namespace entaglement{

    //<! calculates the reduced density matrix of a subsystem with set size
    template <typename _ty>
    inline
    auto reduced_density_matrix(
        const arma::Col<_ty>& state,//<! input state in full Hilbert space
        int A_size,                 //<! subsystem size
        int L,                      //<! system size
        int config = 2              //<! on-site configuration number
        ) 
        -> arma::Mat<_ty> 
        {
        int num_of_bits = log2(config);
    	// set subsytsems size
    	const long long dimA = (ULLPOW( (num_of_bits *      A_size ) ));
    	const long long dimB = (ULLPOW( (num_of_bits * (L - A_size)) ));
        const long long N = dimA * dimB;
    	arma::Mat<_ty> rho(dimA, dimA, arma::fill::zeros);
    	for (long long n = 0; n < N; n++) {						// loop over configurational basis
    		long long counter = 0;
    		for (long long m = n % dimB; m < N; m += dimB) {		// pick out state with same B side (last L-A_size bits)
    			long idx = n / dimB;							        // find index of state with same B-side (by dividing the last bits are discarded)
    			rho(idx, counter) += my_conjungate(state(n)) * state(m);
    			counter++;										        // increase counter to move along reduced basis
    		}
    	}
    	return rho;	
    }


    //<! calculates the reduced density matrix of a subsystem with set size for symmetic case (particle number conservation for now)
    template <typename _ty>
    inline
    auto reduced_density_matrix_sym(
        const arma::Col<_ty>& state, //<! input state in full Hilbert space
        int A_size,                  //<! subsystem size
        int L,                       //<! system size
        const v_1d<u64>& full_map,   //<! mapping to specific symmetry sector
        int config = 2               //<! on-site configuration number
        ) 
        -> arma::Mat<_ty>
        {
    	// set subsytsems size
        int num_of_bits = log2(config);
    	const long long dimA = (ULLPOW( (num_of_bits *      A_size ) ));
    	const long long dimB = (ULLPOW( (num_of_bits * (L - A_size)) ));
        const long long full_dim = dimA * dimB;
        const long long N = full_map.size();

        auto find_index = [&](u64 index){   return binary_search(full_map, 0, N - 1, index);  };

    	arma::Mat<_ty> rho(dimA, dimA, arma::fill::zeros);
    	for (long long n = 0; n < N; n++) {						// loop over configurational basis
    		long long counter = 0;
            const u64 true_n = full_map[n];
    		for (long long j = true_n % dimB; j < full_dim; j += dimB) {	// pick out state with same B side (last L-A_size bits)
    			long idx = true_n / dimB;
                long long m = find_index(j);
                if(m >= 0)
                    rho(idx, counter) += my_conjungate(state(n)) * state(m);
                counter++;  // increase counter to move along reduced basis
    		}
    	}
    	return rho;	
    }

};


//<! entropy calculations from the reduced density matrix
namespace entropy{

    /// @brief von Neumann entropy
    /// @tparam _ty input state type
    /// @param state input state in full Hilbert space
    /// @param A_size subsystem size
    /// @param L system size
    /// @param full_map mapping to specific global symmetry sector (total spin symmetry, particle number symmetry, etc.)
    /// @return entanglement entropy
    template <typename _ty>
    inline 
    double vonNeumann(
        const arma::Col<_ty>& state,                //<! input state in full Hilbert space
        int A_size,                                 //<! subsystem size
        int L,                                      //<! system size
        const v_1d<u64>& full_map = v_1d<u64>()     //<! mapping to specific symmetry sector
        ){
    	arma::Mat<_ty> rho = full_map.empty()? entaglement::reduced_density_matrix(state, A_size, L)
                                                 : entaglement::reduced_density_matrix_sym(state, A_size, L, full_map);
    	arma::vec probabilities;
    	arma::eig_sym(probabilities, rho); //diagonalize to find probabilities and calculate trace in rho's eigenbasis
    	double entropy = 0;
    #pragma omp parallel for reduction(+: entropy)
    	for (int i = 0; i < probabilities.size(); i++) {
    		auto value = probabilities(i);
    		entropy += (abs(value) > 0) ? -value * log(abs(value)) : 0;
    	}
    	//double entropy = -real(trace(rho * real(logmat(rho))));
    	return entropy;
    }
    
    //<! von Neumann entropy for each subsystem size

    /// @brief von Neumann entropy for all subsystem sizes
    /// @tparam _ty input state type
    /// @param state input state in full Hilbert space
    /// @param L system size
    /// @param full_map mapping to specific global symmetry sector (total spin symmetry, particle number symmetry, etc.)
    /// @return entanglement entropy
    template <typename _ty>
    inline
    arma::vec vonNeumann(
        const arma::Col<_ty>& state,                //<! input state in full Hilbert space
        int L,                                      //<! system size
        const v_1d<u64>& full_map = v_1d<u64>()     //<! mapping to specific symmetry sector
    ){
    	arma::vec _entropy(L - 1, arma::fill::zeros);
    //#pragma omp parallel for
    	for (int i = 0; i < L - 1; i++)
    		_entropy(i) = vonNeumann(state, i + 1, L);
    	return _entropy;
    }
    
    
    //<! reyni entropy related to multifractality of the Hilbert space
    //<!    S = 1 / (1 - a) * log2( Tr( q^a ) )

    /// @brief reyni entropy related to multifractality of the Hilbert space: S = 1 / (1 - a) * log2( Tr( q^a ) )
    /// @tparam _ty input state type
    /// @param state input state in full Hilbert space
    /// @param A_size subsystem size
    /// @param L system size
    /// @param alfa reyni entropy order
    /// @param full_map mapping to specific global symmetry sector (total spin symmetry, particle number symmetry, etc.)
    /// @return entanglement entropy
    template <typename _ty>
    inline
    double reyni_entropy(
        const arma::Col<_ty>& state,
        int A_size,
        int L,
        int alfa,
        const v_1d<u64>& full_map = v_1d<u64>()
        ) {
        assert(alfa > 1 && "Only alfa>=2 powers are possible");
    	arma::Mat<_ty> rho = full_map.empty()? entaglement::reduced_density_matrix(state, A_size, L)
                                                : entaglement::reduced_density_matrix_sym(state, A_size, L, full_map);
        rho = arma::powmat(rho, alfa);
        return log2(real(arma::trace(rho))) / (1.0 - alfa);
    }
    
    
    //<! Shannon entropy used in information theory
    //<!    q = \sum_n w_n |n><n|
    //<!    S = -\sum_n Tr( |w_n|^2 log2(|w_n|^2) )
    
    /// @brief Shannon entropy used in information theory
    /// @tparam _ty input state type
    /// @param state input state in full Hilbert space
    /// @param A_size subsystem size
    /// @param L system size
    /// @param full_map mapping to specific global symmetry sector (total spin symmetry, particle number symmetry, etc.)
    /// @return entanglement entropy
    template <typename _ty>
    inline
    double shannon_entropy(
        const arma::Col<_ty>& state,                //<! input state in full Hilbert space
        int A_size,                                 //<! subsystem size
        int L,                                      //<! system size
        const v_1d<u64>& full_map = v_1d<u64>()     //<! mapping to specific symmetry sector
        ) {
    	arma::Mat<_ty> rho = full_map.empty()? entaglement::reduced_density_matrix(state, A_size, L)
                                              : entaglement::reduced_density_matrix_sym(state, A_size, L, full_map);
        arma::vec probabilities;
        arma::eig_sym(probabilities, rho); //diagonalize to find probabilities and calculate trace in rho's eigenbasis
        double _entropy = 0;
    #pragma omp parallel for reduction(+: _entropy)
        for (int i = 0; i < probabilities.size(); i++) {
        	auto value = abs(probabilities(i) * probabilities(i));
        	_entropy += ((value) < 1e-10) ? 0 : -value * log2(value);
        }
        return _entropy;
    }
    
    /// @brief Calculates the entropy using the Schmidt decomposition of a wavefunction
    /// @tparam _ty input state type
    /// @param state input state in full Hilbert space
    /// @param A_size subsystem size
    /// @param L system size
    /// @param config on-site configuration number (local hilbert space)
    /// @return entanglement entropy
    template <typename _ty>
    inline
    auto schmidt_decomposition(
        const arma::Col<_ty>& state,
        int A_size,
        int L,
        int config = 2
        )
    {
        int num_of_bits = std::log2(config);
    	const long long dimA = (ULLPOW( (num_of_bits *      A_size ) ));
    	const long long dimB = (ULLPOW( (num_of_bits * (L - A_size)) ));

        // reshape array to matrix
        arma::Mat<_ty> rho = arma::reshape(state, dimA, dimB);

        // get schmidt coefficients from singular-value-decomposition
        arma::vec schmidt_coeff = arma::svd(rho);

        //calculate entropy
        double entropy = 0;
    #pragma omp parallel for reduction(+: entropy)
    	for (int i = 0; i < schmidt_coeff.size(); i++) {
    		auto value = schmidt_coeff(i) * schmidt_coeff(i);
    		entropy += (abs(value) > 0) ? -value * std::log(value) : 0;
    	}
        return entropy;
    }
};