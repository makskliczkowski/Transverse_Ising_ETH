#include "include/IsingModel.h"

// ----------------------------------------------------------------------------- CONSTRUCTORS -----------------------------------------------------------------------------
IsingModel_disorder::IsingModel_disorder(int L, double J, double J0, double g, double g0, double h, double w, int _BC) {
	this->L = L; this->J = J; this->g = g; this->h = h;
	this->J0 = J0; this->g0 = g0;  this->w = w;
	this->N = ULLPOW(this->L);
	this->_BC = _BC;

	//change info
	this->info = "_L=" + std::to_string(this->L) + \
		",J=" + to_string_prec(this->J, 2) + \
		",J0=" + to_string_prec(this->J0, 2) + \
		",g=" + to_string_prec(this->g, 2) + \
		",g0=" + to_string_prec(this->g0, 2) + \
		",h=" + to_string_prec(this->h, 2) + \
		",w=" + to_string_prec(this->w, 2);
	this->set_neighbors();
	#ifndef HEISENBERG
		if(this->g == 0 && this->g0 == 0)
			generate_mapping();
	#else
		generate_mapping();
	#endif	
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
	#ifndef HEISENBERG
		return index;
		return (this->g == 0 && this->g0 == 0? mapping[index] : index);
	#else
		return mapping[index];
	#endif
}

u64 IsingModel_disorder::find_in_map(u64 index) const {
	#ifdef HEISENBERG
		return binary_search(this->mapping, 0, this->N - 1, index);
	#else
		return index;
		return (this->g == 0 && this->g0 == 0)? binary_search(this->mapping, 0, this->N - 1, index) : index;
	#endif
}

/// <summary>
/// W razie gdyby�my robili Sz symetri�, dla picu
/// </summary>
void IsingModel_disorder::generate_mapping() {
	this->mapping = std::vector<u64>();
	for (u64 j = 0; j < (ULLPOW(this->L)); j++){
		if (__builtin_popcountll(j) == this->L / 2.)
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
	try{
		H(idx, k) += value;
		#ifdef HEISENBERG
			H(k, idx) += value;
		#endif
	} 
	catch (const std::exception& err) {
		stout << "Exception:\t" << err.what() << "\n";
		stout << "SHit ehhh..." << std::endl;
		printSeparated(std::cout, "\t", 14, true, new_idx, idx, k, value);
	}

}

/// <summary>
/// Generates the total Hamiltonian of the system. The diagonal part is straightforward,
/// while the non-diagonal terms need the specialized setHamiltonainElem(...) function
/// </summary>
void IsingModel_disorder::hamiltonian() {
	try {
		this->H = arma::sp_mat(N, N);                                //  hamiltonian memory reservation
	}
	catch (const bad_alloc& e) {
		std::cout << "Memory exceeded" << e.what() << "\n";
		assert(false);
	}
	#pragma omp critical
	{
		this->dh = arma::vec(this->L, arma::fill::zeros);
		this->dJ = arma::vec(this->L, arma::fill::zeros);
		this->dg = arma::vec(this->L, arma::fill::zeros);
		#if defined(LOCAL_PERT)
			this->dh(this->L / 2.) = this->w;
			this->dh(0) = 0.1;
		#else
			this->dh = create_random_vec(L, this->w);                               // creates random disorder vector
			this->dJ = create_random_vec(L, this->J0);                              // creates random exchange vector
			this->dg = create_random_vec(L, this->g0);                              // creates random transverse field vector
		#endif
	}
	#ifdef ANDERSON
		auto lattice = std::make_unique<lattice3D>(this->L);
		this->H = (arma::mat)anderson::hamiltonian(*lattice, this->J, this->w);
	#else
		#ifdef HEISENBERG
			this->hamiltonian_heisenberg();
		#else
			this->hamiltonian_Ising();
		#endif
	#endif
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

