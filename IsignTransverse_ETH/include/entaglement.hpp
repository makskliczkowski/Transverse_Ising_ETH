
namespace entaglement{

    //<! calculates the reduced density matrix of a subsystem with set size
    inline
    auto reduced_density_matrix(
        const arma::cx_vec& state,  //<! input state in full Hilbert space
        int A_size,                 //<! subsystem size
        int L,                      //<! system size
        int config = 2              //<! on-site configuration number
        ) 
        -> arma::cx_mat 
        {
        int num_of_bits = log2(config);
    	// set subsytsems size
    	const long long dimA = (ULLPOW( (num_of_bits *      A_size ) ));
    	const long long dimB = (ULLPOW( (num_of_bits * (L - A_size)) ));
        const long long N = dimA * dimB;
    	arma::cx_mat rho(dimA, dimA, arma::fill::zeros);
    	for (long long n = 0; n < N; n++) {						// loop over configurational basis
    		long long counter = 0;
    		for (long long m = n % dimB; m < N; m += dimB) {		// pick out state with same B side (last L-A_size bits)
    			long idx = n / dimB;							        // find index of state with same B-side (by dividing the last bits are discarded)
    			rho(idx, counter) += conj(state(n)) * state(m);
    			counter++;										        // increase counter to move along reduced basis
    		}
    	}
    	return rho;	
    }

};


//<! entropy calculations from the reduced density matrix
namespace entropy{

    //<! typical von Neumann entropy used in statistical physics
    //<!    q = \sum_n w_n |n><n|
    //<!    S = -\sum_n Tr( w_n ln(w_n) )
    inline 
    double vonNeumann(
        const arma::cx_vec& state,  //<! input state in full Hilbert space
        int A_size,                 //<! subsystem size
        int L                       //<! system size
        ){
    	arma::cx_mat rho = entaglement::reduced_density_matrix(state, A_size, L);
    	arma::vec probabilities;
    	arma::eig_sym(probabilities, rho); //diagonalize to find probabilities and calculate trace in rho's eigenbasis
    	double entropy = 0;
    #pragma omp parallel for reduction(+: entropy)
    	for (int i = 0; i < probabilities.size(); i++) {
    		auto value = probabilities(i);
    		entropy += (abs(value) < 1e-10) ? 0 : -value * log(abs(value));
    	}
    	//double entropy = -real(trace(rho * real(logmat(rho))));
    	return entropy;
    }
    
    //<! von Neumann entropy for each subsystem size
    inline
    arma::vec vonNeumann(
        const arma::cx_vec& state,  //<! input state in full Hilbert space
        int L                       //<! system size
    ){
    	arma::vec _entropy(L - 1, arma::fill::zeros);
    //#pragma omp parallel for
    	for (int i = 0; i < L - 1; i++)
    		_entropy(i) = vonNeumann(state, i + 1, L);
    	return _entropy;
    }
    
    
    //<! reyni entropy related to multifractality of the Hilbert space
    //<!    S = 1 / (1 - a) * log2( Tr( q^a ) )
    inline
    double reyni_entropy(
        const arma::cx_vec& state,  //<! input state in full Hilbert space
        int A_size,                 //<! subsystem size
        int L,                      //<! system size
        int alfa               //<! 
        ) {
        assert(alfa > 1 && "Only alfa>=2 powers are possible");
        arma::cx_mat rho = arma::powmat(entaglement::reduced_density_matrix(state, A_size, L), alfa);
        return log2(real(arma::trace(rho))) / (1.0 - alfa);
    }
    
    
    //<! Shannon entropy used in information theory
    //<!    q = \sum_n w_n |n><n|
    //<!    S = -\sum_n Tr( |w_n|^2 log2(|w_n|^2) )
    inline
    double shannon_entropy(
        const arma::cx_vec& state,  //<! input state in full Hilbert space
        int A_size,                 //<! subsystem size
        int L                       //<! system size
        ) {
        arma::cx_mat rho = entaglement::reduced_density_matrix(state, A_size, L);
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
    

};