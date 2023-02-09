#include "include/IsingModel.h"

// ----------------------------------------------------------------------------- CONSTRUCTORS -----------------------------------------------------------------------------
IsingModel_disorder::IsingModel_disorder(int L, double J, double J0, double g, double g0, double h, double w, int _BC) {
	this->L = L; this->J = J; this->g = g; this->h = h;
	this->J0 = J0; this->g0 = g0;  this->w = w;
	this->_BC = _BC;
	
	this->use_real_matrix = true;

	//change info
	#ifdef ANDERSON
		this->info = "_L=" + std::to_string(this->L) + \
			",J=" + to_string_prec(this->J, 2) + \
			",J0=" + to_string_prec(this->J0, 2) + \
			",w=" + to_string_prec(this->w, 2);
		this->N = this->L * this->L * this->L;
	#else
		this->info = "_L=" + std::to_string(this->L) + \
			",J=" + to_string_prec(this->J) + \
			",J0=" + to_string_prec(this->J0) + \
			",g=" + to_string_prec(this->g) + \
			",g0=" + to_string_prec(this->g0) + \
			",h=" + to_string_prec(this->h) + \
			",w=" + to_string_prec(this->w);
		this->N = ULLPOW(this->L);
	#endif
	this->set_neighbors();
	this->use_Sz_sym = false;
	#ifdef HEISENBERG		// HEISENBERG
		this->use_Sz_sym = true;
		std::cout << "Using Heisenberg with Sz sym" << std::endl;
	#elif defined(XYZ)		// XYZ
		if(this->J0 == 0 && this->g0 == 0){
			this->use_Sz_sym = true;
			std::cout << "Using XYZ with Sz sym" << std::endl;
		} else
			std::cout << "Using XYZ" << std::endl;
	#elif defined(ANDERSON)
			std::cout << "Using 3D Anderson model" << std::endl;
	#elif defined(QUANTUM_SUN)
		this->use_Sz_sym = false;
		std::cout << "Using Quantum Sun" << std::endl;
	#else					//ISING
		if(this->g == 0 && this->g0 == 0){
			this->use_Sz_sym = true;
			std::cout << "Using Ising with Sz sym" << std::endl;
		} else
			std::cout << "Using Ising" << std::endl;
	#endif	
	if(use_Sz_sym)
		generate_mapping();
	std::cout << "Spectrum size:\t" << this->N << std::endl;
	this->hamiltonian();
}

// ----------------------------------------------------------------------------- BASE GENERATION AND RAPPING -----------------------------------------------------------------------------

/// <summary>
/// Return the index in the case of no mapping in disorder
/// </summary>
/// <param name="index"> index to take</param>
/// <returns>index</returns>
u64 IsingModel_disorder::map(u64 index) const {
	if (index < 0 || index >= std::pow(2, L)) throw "Element out of range\n No such index in map\n";
	return this->use_Sz_sym? mapping[index] : index;
}

u64 IsingModel_disorder::find_in_map(u64 index) const {
	return this->use_Sz_sym? binary_search(this->mapping, 0, this->N - 1, index) : index;
}

/// <summary>
/// W razie gdyby�my robili Sz symetri�, dla picu
/// </summary>
void IsingModel_disorder::generate_mapping() {
	this->mapping = std::vector<u64>();
	for (u64 j = 0; j < (ULLPOW(this->L)); j++){
		if (__builtin_popcountll(j) == int(this->L / 2.))
			this->mapping.push_back(j);
	}
	this->N = this->mapping.size();
}

// ----------------------------------------------------------------------------- BUILDING HAMILTONIAN -----------------------------------------------------------------------------

/// <summary>
/// Sets the non-diagonal elements of the Hamimltonian matrix, by acting with the operator on the k-th state
/// </summary>
/// <param name="k"> index of the basis state acted upon with the Hamiltonian </param>
/// <param name="value"> value of the given matrix element to be set </param>
/// <param name="temp"> resulting vector form acting with the Hamiltonian operator on the k-th basis state </param>
void IsingModel_disorder::setHamiltonianElem(u64 k, double value, u64 new_idx) {
	u64 idx = find_in_map(new_idx);
	#ifdef XYZ
		if(this->use_Sz_sym && idx >= this->N)
			return;
	#endif
	try{
		H(idx, k) += value;
		#ifdef HEISENBERG
			H(k, idx) += value;
		#endif
	} 
	catch (const std::exception& err) {
		stout << "Exception:\t" << err.what() << "\n";
		stout << "SHit ehhh..." << std::endl;
		printSeparated(std::cout, "\t", 14, true, new_idx, idx, map(k), value);
	}

}

void IsingModel_disorder::hamiltonian() {
	#ifdef ANDERSON
		auto lattice = std::make_unique<lattice3D>(this->L);
		this->N = lattice->volume;
	#endif
	try {
		this->H = arma::sp_mat(N, N);                                //  hamiltonian memory reservation
	}
	catch (const bad_alloc& e) {
		std::cout << "Memory exceeded" << e.what() << "\n";
		assert(false);
	}
	#pragma omp critical
	{
		this->dh = disorder::uniform<double>(this->L, 0);
		this->dJ = disorder::uniform<double>(this->L, 0);
		this->dg = disorder::uniform<double>(this->L, 0);
		#if defined(LOCAL_PERT)
			this->dh(this->L / 2.) = this->w;
			this->dh(0) = 0.1;
		#else
			this->dh.reset(this->L, this->w);                               // creates random disorder vector
			this->dJ.reset(this->L, this->J0);                              // creates random exchange vector
			this->dg.reset(this->L, this->g0);                              // creates random transverse field vector
		#endif
	}
	#if MODEL == 0
		this->hamiltonian_Ising();
	#elif MODEL == 1
		this->hamiltonian_heisenberg();
	#elif MODEL == 2
		this->H = anderson::hamiltonian(*lattice, this->J, this->w);
	#elif MODEL == 3
		this->hamiltonian_xyz();
	#elif MODEL == 4
		this->hamiltonian_qsun();
	#else
		this->hamiltonian_Ising();
	#endif
}

/// <summary>
/// Generates the total Hamiltonian of the system. The diagonal part is straightforward,
/// while the non-diagonal terms need the specialized setHamiltonainElem(...) function
/// </summary>
void IsingModel_disorder::hamiltonian_qsun()
{
	int M = 3; // bubble size
	const size_t dim_loc = ULLPOW((L - M));
	const size_t dim_erg = ULLPOW((M));
	
	/* Create GOE Matrix */
	arma::mat H_grain = my_gen.create_goe_matrix<double>(dim_erg);
	

	/* Create Rest of Hamiltonian */
	arma::Col<int> random_neigh(this->L, arma::fill::zeros);
	for(int i = 0; i < this->L; i++)
		random_neigh(i) = i < M? -1 : my_gen.randomInt_uni(0, M - 1);

	std::cout << random_neigh.t() << std::endl;

	for (int j = 0; j < L; j++) 
	 	this->dJ(j) = j < M? 0.0 : std::pow(this->g, j - M + 1 + this->dJ(j));

	for (u64 k = 0; k < N; k++) {
		double s_i, s_j;
		u64 base_state = map(k);
		for (int j = M; j < L; j++) {
			/* disorder on localised spins */
            auto [val, Sz_k] = operators::sigma_z(base_state, this->L, { j });
			this->setHamiltonianElem(k, (this->h + this->dh(j)) * real(val), Sz_k);

			/* coupling of localised spins to GOE grain */
			int nei = random_neigh(j);
		    auto [val1, Sx_k] = operators::sigma_x(base_state, this->L, { j });
		    auto [val2, SxSx_k] = operators::sigma_x(Sx_k, this->L, { nei });
			this->setHamiltonianElem(k, this->J * this->dJ(j) * real(val1 * val2), SxSx_k);
		}
	}
	this->H = this->H + arma::kron(H_grain, arma::eye(dim_loc, dim_loc));
}

void IsingModel_disorder::hamiltonian_xyz(){
	std::cout << this->N << std::endl;
	std::vector<std::vector<double>> parameters = { { 1.0 * (1 - this->J0), 1.0 * (1 + this->J0), this->g},
                                                    {this->J * (1 - this->J0), this->J * (1 + this->J0), this->J * this->g }
                                                };
    for(auto& x : parameters)
        std::cout << x << std::endl;
    std::vector<op_type> XYZoperators = {operators::sigma_x, operators::sigma_y, operators::sigma_z };
    std::vector<int> neighbor_distance = {1, 2};

    for (size_t k = 0; k < this->N; k++) {
		int base_state = map(k);
	    for (int j = 0; j < this->L; j++) {
            cpx val = 0.0;
            u64 op_k;
            std::tie(val, op_k) = operators::sigma_z(base_state, this->L, { j });
			double fieldZ = (j == this->L - 1)? this->w : this->h;
            this->setHamiltonianElem(k, fieldZ * real(val), op_k);
	    	
            std::tie(val, op_k) = operators::sigma_x(base_state, this->L, { j });			
            this->setHamiltonianElem(k, this->g0 * real(val), op_k);

            for(int a = 0; a < neighbor_distance.size(); a++){
                int r = neighbor_distance[a];
				int nei = j + r;
				if(nei >= this->L)
					nei = (this->_BC)? -1 : nei % this->L;

	    	    if (nei >= 0) {
                    for(int b = 0; b < XYZoperators.size(); b++){
                        op_type op = XYZoperators[b];
		                auto [val1, op_k] = op(base_state, this->L, { j });
		                auto [val2, opop_k] = op(op_k, this->L, { nei });
						this->setHamiltonianElem(k, parameters[a][b] * real(val1 * val2), opop_k);
                    }
	    	    }
            }
	    }
	}
}
void IsingModel_disorder::hamiltonian_Ising() {
	for (u64 k = 0; k < N; k++) {
		double s_i, s_j;
		u64 base_state = map(k);
		for (int j = 0; j <= L - 1; j++) {
			s_i = checkBit(base_state, L - 1 - j) ? 1.0 : -1.0;;							 // true - spin up, false - spin down

			if( this->g != 0 || this->dg(j) != 0){
				u64 new_idx = flip(base_state, BinaryPowers[this->L - 1 - j], this->L - 1 - j);
				setHamiltonianElem(k, this->g + this->dg(j), new_idx);
			}
			/* disorder */
			H(k, k) += (this->h + dh(j)) * s_i;                             // diagonal elements setting

			if (nearest_neighbors[j] >= 0) {
				/* Ising-like spin correlation */
				s_j = checkBit(base_state, this->L - 1 - nearest_neighbors[j]) ? 1.0 : -1.0;
				this->H(k, k) += (this->J + this->dJ(j)) * s_i * s_j;		// setting the neighbors elements
			}
		}
	}
}
void IsingModel_disorder::hamiltonian_heisenberg(){
	for (u64 k = 0; k < N; k++) {
		double s_i, s_j;
		u64 base_state = map(k);
		for (int j = 0; j <= L - 1; j++) {
			s_i = checkBit(base_state, L - 1 - j) ? 0.5 : -0.5;							 // true - spin up, false - spin down
			/* disorder */
			H(k, k) += (this->h + this->dh(j)) * s_i;                             // diagonal elements setting

			int nei = this->nearest_neighbors[j];
			if (nei >= 0) {
				s_j = checkBit(base_state, L - 1 - nei) ? 0.5 : -0.5;
				if(s_i < 0 && s_j > 0){
					u64 new_idx =  flip(base_state, BinaryPowers[this->L - 1 - nei], this->L - 1 - nei);
					new_idx =  flip(new_idx, BinaryPowers[this->L - 1 - j], this->L - 1 - j);
					// 0.5 cause flip 0.5*(S+S- + S-S+)
					setHamiltonianElem(k, 0.5 * (this->J + this->dJ(j)), new_idx);
				}
				
				/* Ising-like spin correlation */
				H(k, k) += (this->g + this->dg(j)) * s_i * s_j;
			}
		}
		//std::cout << std::bitset<4>(base_state) << "\t";
	}
	//std::cout << std::endl << arma::mat(this->H) << std::endl;
}


void IsingModel_disorder::diag_sparse(bool get_eigenvectors, int maxiter, double tol, double sigma){
	#ifdef ARMA_USE_SUPERLU
		arma::eigs_opts opts;
		opts.maxiter = maxiter;
		opts.tol = tol;
		int num = 500;
		arma::eigs_sym(this->eigenvalues, this->eigenvectors, this->H, num, sigma, opts);
	#else
		#pragma message ("No SuperLu defined. Ignoring shift-invert possibility.")
	#endif
}
// ----------------------------------------------------------------------------- PHYSICAL QUANTITES -----------------------------------------------------------------------------



// ----------------------------------------------------------------------------- CREATE OPERATOR TO CALCULATE MATRIX ELEMENTS -----------------------------------------------------------------------------
arma::sp_cx_mat IsingModel_disorder::create_operator(std::initializer_list<op_type> operators, std::vector<int> sites) const {
	arma::sp_cx_mat opMatrix(this->N, this->N);
#pragma omp parallel for
	for (long int k = 0; k < N; k++) {
		u64 base_state = map(k);
		for (auto& op : operators) {
			cpx value; u64 new_idx;
			std::tie(value, new_idx) = op(base_state, this->L, sites);
			u64 idx = find_in_map(new_idx);
			if(idx > this->N) continue;
		#pragma omp critical
			opMatrix(idx, k) += value;
		}
	}
	return opMatrix;
}
arma::sp_cx_mat IsingModel_disorder::create_operator(std::initializer_list<op_type> operators, arma::cx_vec prefactors) const {
	arma::sp_cx_mat opMatrix(this->N, this->N);
	arma::cx_vec pre = prefactors.is_empty()? arma::cx_vec(this->L, arma::fill::ones) : prefactors;
	assert(pre.size() == this->L && "Input array of different size than system size!");
#pragma omp parallel for
	for (long int k = 0; k < N; k++) {
		u64 base_state = map(k);
		for (int j = 0; j < this->L; j++) {
			for (auto& op : operators) {
				cpx value; u64 new_idx;
				std::tie(value, new_idx) = op(base_state, this->L, { j });
				u64 idx = find_in_map(new_idx);
				if(idx > this->N) continue;
			#pragma omp critical
				opMatrix(idx, k) += value * pre(j);
			}
		}
	}
	return opMatrix / sqrt(this->L);
}
arma::sp_cx_mat IsingModel_disorder::create_operator(std::initializer_list<op_type> operators, int corr_len, arma::cx_vec prefactors) const {
	arma::sp_cx_mat opMatrix(this->N, this->N);
	auto neis = get_neigh_vector(this->_BC, this->L, corr_len);
	arma::cx_vec pre = prefactors.is_empty()? arma::cx_vec(this->L, arma::fill::ones) : prefactors;
	assert(pre.size() == this->L && "Input array of different size than system size!");
#pragma omp parallel for
	for (long int k = 0; k < N; k++) {
		u64 base_state = map(k);
		for (int j = 0; j < this->L; j++) {
			const int nei = neis[j];
			if (nei < 0) continue;
			for (auto& op : operators) {
				cpx value; u64 new_idx;
				std::tie(value, new_idx) = op(base_state, this->L, { j, nei });
				u64 idx = find_in_map(new_idx);
				if(idx > this->N) continue;
			#pragma omp critical
				opMatrix(idx, k) += value * pre(j);
			}
		}
	}
	return opMatrix / sqrt(this->L);
}

arma::sp_cx_mat IsingModel_disorder::spin_current() const {
	arma::sp_cx_mat opMatrix(this->N, this->N);
	auto neis = get_neigh_vector(this->_BC, this->L, 1);
#pragma omp parallel for
	for (long int k = 0; k < N; k++) {
		u64 base_state = map(k);
		for (int j = 0; j < this->L; j++) {
			const int nei = neis[j];
			if(nei < 0) continue;
			cpx value_x, value_y; 
			u64 new_idx_x, new_idx;
			// sigma^x_j sigma^y_l+1
			std::tie(value_x, new_idx_x) = IsingModel_disorder::sigma_y(base_state, this->L, { nei });
			std::tie(value_y, new_idx) = IsingModel_disorder::sigma_x(new_idx_x, this->L, { j });

			u64 idx = find_in_map(new_idx);
			if(idx > this->N) continue;
		#pragma omp critical
			opMatrix(idx, k) += value_x * value_y;

			// sigma^y_j sigma^x_l+1
			std::tie(value_x, new_idx_x) = IsingModel_disorder::sigma_x(base_state, this->L, { nei });
			std::tie(value_y, new_idx) = IsingModel_disorder::sigma_y(new_idx_x, this->L, { j });

			idx = find_in_map(new_idx);
			if(idx > this->N) continue;
		#pragma omp critical
			opMatrix(idx, k) -= value_x * value_y;
		}
	}
	return 2.0 * S * this->J * opMatrix / sqrt(this->L); // factor of 2 cause of commutation of pauli operators
}
arma::sp_cx_mat IsingModel_disorder::fourierTransform(op_type op, int k) const {
	const double q = k * two_pi / double(this->L);
	arma::cx_vec k_exp(this->L);
	for (int j = 0; j < this->L; j++) {
		k_exp(j) = std::exp(im * q * double(j + 1));
	}
	return create_operator({op}, k_exp);
}
arma::sp_cx_mat IsingModel_disorder::createHq(int k) const {
	const double q = k * two_pi / double(this->L);
	arma::sp_cx_mat opMatrix(this->N, this->N);
	std::vector<cpx> k_exp(this->L);
	for (int j = 0; j < this->L; j++)
		k_exp[j] = std::cos(q * double(j + 1));
#pragma omp parallel for
	for (long int n = 0; n < this->N; n++) {
		u64 base_state = map(N);
		cpx __2 = 0.0, __3 = 0.0;
		u64 new_idx1 = 0, new_idx2 = 0;
		for (int j = 0; j < this->L; j++) {
			auto nei = this->nearest_neighbors[j];			
			auto [s_i, __1] = IsingModel_disorder::sigma_z(base_state, this->L, { j });
			opMatrix(n, n) += this->h / 2. * s_i * k_exp[j];

			std::tie(__2, new_idx1) = IsingModel_disorder::sigma_x(base_state, this->L, { j });
			u64 idx1 = find_in_map(new_idx1);
		#pragma omp critical
			{
				opMatrix(idx1, n) += this->g / 2. * k_exp[j];
			}
			if (nei >= 0) {
				std::tie(__3, new_idx2) = IsingModel_disorder::sigma_x(base_state, this->L, { nei });
				u64 idx2 = find_in_map(new_idx2);
			#pragma omp critical
				{
					opMatrix(idx2, n) += this->g / 2. * k_exp[j];
				}
				auto [s_j, __4] = IsingModel_disorder::sigma_z(base_state, this->L, { nei });

				opMatrix(n, n) += (this->J * s_i + this->h / 2.) * s_j * k_exp[j];
			}
		}
	}
	return opMatrix / std::sqrt(this->L);
}
arma::sp_cx_mat IsingModel_disorder::createHlocal(int k) const {
	arma::sp_cx_mat opMatrix(this->N, this->N);
	#ifdef HEISENBERG
#pragma omp parallel for
	for (long int n = 0; n < this->N; n++) {
		u64 base_state = map(n);
		cpx __2 = 0.0, __3 = 0.0;
		u64 new_idx1 = 0, new_idx2 = 0;
		auto nei = this->nearest_neighbors[k];
		auto [s_iii, __1] = IsingModel_disorder::sigma_z(base_state, this->L, { k });
		double s_i = real(s_iii);
		opMatrix(n, n) += (this->h + this->dh(k)) / 2. * s_i;

		
		if (nei >= 0) {
			auto [s_jjj, __4] = IsingModel_disorder::sigma_z(base_state, this->L, { nei });
			double s_j = real(s_jjj);
			opMatrix(n, n) += ((this->g + this->g0) * s_i + (this->h + this->dh(nei)) / 2.) * s_j;
			if(s_i < 0 && s_j > 0){
				u64 new_idx =  flip(base_state, BinaryPowers[this->L - 1 - nei], this->L - 1 - nei);
				new_idx =  flip(new_idx, BinaryPowers[this->L - 1 - k], this->L - 1 - k);
				// 0.5 cause flip 0.5*(S+S- + S-S+)
				u64 idx = find_in_map(new_idx);
				opMatrix(idx, n) += 0.5 * (this->J + this->dJ(k));
				opMatrix(n, idx) += 0.5 * (this->J + this->dJ(k));
			}
		}
	}

	#else
#pragma omp parallel for
	for (long int n = 0; n < this->N; n++) {
		u64 base_state = map(n);
		cpx __2 = 0.0, __3 = 0.0;
		u64 new_idx1 = 0, new_idx2 = 0;
		auto nei = this->nearest_neighbors[k];
		auto [s_i, __1] = IsingModel_disorder::sigma_z(base_state, this->L, { k });
		opMatrix(n, n) += this->h / 2. * s_i;

		std::tie(__2, new_idx1) = IsingModel_disorder::sigma_x(base_state, this->L, { k });
		u64 idx1 = find_in_map(new_idx1);
	#pragma omp critical
		{
			opMatrix(idx1, n) += this->g / 2.;
		}
		if (nei >= 0) {
			std::tie(__3, new_idx2) = IsingModel_disorder::sigma_x(base_state, this->L, { nei });
			u64 idx2 = find_in_map(new_idx2);
		#pragma omp critical
			{
				opMatrix(idx2, n) += this->g / 2.;
			}
			auto [s_j, __4] = IsingModel_disorder::sigma_z(base_state, this->L, { nei });

			opMatrix(n, n) += (this->J * s_i + this->h / 2.) * s_j;
		}
	}
	#endif
	return opMatrix;
}


// ----------------------------------------------------------------------------- INTEGRABLE SOLUTIONS -----------------------------------------------------------------------------
arma::vec IsingModel_disorder::get_non_interacting_energies(){
	const u64 dim = ULLPOW(this->L);
	arma::vec energies(dim);
	auto epsilon = [this](int site){
		double hi = this->h + this->dh(site);
		return sqrt(this->g * this->g + hi * hi);
	};
	int counter = 0;
	for(int k = 0; k <= L; k++) // number of flipped bits
	{
		std::vector<int> bitmask(k, 1); 	// string with k-leading 1's
		bitmask.resize(this->L, 0);			// L-k trailing 0's

		// -------- permute all binary representations to get all combinations of subset of size k
    	do {
			double E = 0.0;
    	    for (int i = 0; i < this->L; i++)  {
    	        if (bitmask[i]) E += epsilon(i);
				else 			E -= epsilon(i);
    	    }
			energies(counter++) = E;
    	} while (std::prev_permutation(bitmask.begin(), bitmask.end()));
	}
	sort(energies.begin(), energies.end());
	return energies;
}

