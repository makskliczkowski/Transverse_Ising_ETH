#include "include/IsingModel.h"

// ------------------------------------------------------------------------------------------------ CONSTRUCTORS ------------------------------------------------------------------------------------------------
IsingModel_sym::IsingModel_sym(int L, double J, double g, double h, int k_sym, bool p_sym, bool x_sym, int _BC) {
	this->L = L; this->J = J; this->g = g; this->h = h; this->_BC = _BC;
	
	symmetries.p_sym = (p_sym) ? 1 : -1;
	symmetries.x_sym = (x_sym) ? 1 : -1;
	this->info = "_L=" + std::to_string(this->L) + \
		",J=" + to_string_prec(this->J) + \
		",g=" + to_string_prec(this->g) + \
		",h=" + to_string_prec(this->h) + \
		",k=" + std::to_string(k_sym) + \
		",p=" + std::to_string(symmetries.p_sym) + \
		",x=" + std::to_string(symmetries.x_sym);
	if(k_sym == -1)
		k_sym = this->L / 2;
	symmetries.k_sym = k_sym * two_pi / double(this->L);

	k_sector = abs(this->symmetries.k_sym) < 1e-4 || abs(this->symmetries.k_sym - pi) < 1e-4;
	// precalculate the exponents
	this->k_exponents = v_1d<cpx>(this->L, 0.0);
#pragma omp parallel for
	for (int l = 0; l < this->L; l++)
		this->k_exponents[l] = std::exp(-im * this->symmetries.k_sym * double(l));
	this->createSymmetryGroup();


	this->mapping = std::vector<u64>();
	this->normalisation = std::vector<cpx>();

	this->use_Sz_sym = false;
	#ifdef HEISENBERG		// HEISENBERG
		this->use_Sz_sym = true;
	#elif defined(XYZ)		// XYZ
		this->use_Sz_sym = false;
	#else					//ISING
		if(this->g == 0)
			this->use_Sz_sym = true;
	#endif	
	
	this->generate_mapping();
	if (this->N <= 0) {
		std::cout << "No states in Hilbert space" << std::endl;
		return;
	}
	this->set_neighbors(); // generate neighbors
	this->hamiltonian();
}

// ------------------------------------------------------------------------------------------------ BASE GENERATION, SEC CLASSES AND RAPPING ------------------------------------------------------------------------------------------------
void IsingModel_sym::createSymmetryGroup() {
	this->symmetry_group = std::vector<std::function<u64(u64, int)>>();
	this->symmetry_eigVal = std::vector<cpx>();
	std::function<u64(u64, int)> e = [](u64 n, int L)->u64 {return n; };		// neutral element
	std::function<u64(u64, int)> T = e;											// in 1st iteration is neutral element
	std::function<u64(u64, int)> Z = static_cast<u64(*)(u64, int)>(&flip);		// spin flip operator (all spins)
	std::function<u64(u64, int)> P = reverseBits;								// parity operator
	if(this->_BC){
		this->symmetry_group.push_back(e);
		this->symmetry_eigVal.push_back(1.0);
		this->symmetry_group.push_back(P);
		this->symmetry_eigVal.push_back(double(this->symmetries.p_sym));
		if (this->h == 0) {
				this->symmetry_group.push_back(Z);
				this->symmetry_eigVal.push_back(double(this->symmetries.x_sym));

				this->symmetry_group.push_back(multiply_operators(P, Z));
				this->symmetry_eigVal.push_back(double(this->symmetries.p_sym * (long)this->symmetries.x_sym));
		}

	}
	else{
		for (int k = 0; k < this->L; k++) {
					this->symmetry_group.push_back(T);
					this->symmetry_eigVal.push_back(this->k_exponents[k]);
			if (this->h == 0) {
					this->symmetry_group.push_back(multiply_operators(Z, T));
					this->symmetry_eigVal.push_back(this->k_exponents[k] * double(this->symmetries.x_sym));
			}
			if (this->k_sector) {
					this->symmetry_group.push_back(multiply_operators(P, T));
					this->symmetry_eigVal.push_back(this->k_exponents[k] * double(this->symmetries.p_sym));
				if (this->h == 0) {
					this->symmetry_group.push_back(multiply_operators(multiply_operators(P, Z), T));
					NO_OVERFLOW(this->symmetry_eigVal.push_back(this->k_exponents[k] * double(this->symmetries.p_sym * (long)this->symmetries.x_sym));)
				}
			}
			T = multiply_operators(rotate_left, T);
		}
	}
}


/// <summary>
/// Takes the state from the index in mapping
/// </summary>
u64 IsingModel_sym::map(u64 index) const {
	if (index >= this->N) throw "Element out of range\n No such index in map\n";
	return this->mapping[index];
}

/// <summary>
/// Find representatives of other EC generated by reflection, spin-flip and (reflection x spin-flip) symmetry
/// </summary>
/// <param name="base_vector"> current base vector to act with symmetries </param>
/// <param name="min"> index of EC class representative by translation symmetry </param>
/// <returns></returns>
std::pair<u64, cpx> IsingModel_sym::find_SEC_representative(u64 base_idx) const {
	u64 SEC = INT64_MAX;
	int _min = INT_MAX;
	for (int l = 0; l < this->symmetry_group.size(); l++) {
		u64 new_idx = this->symmetry_group[l](base_idx, this->L);
		if (new_idx < SEC) {
			SEC = new_idx;
			_min = l;
		}
	}
	return std::make_pair(SEC, this->symmetry_eigVal[_min]);
}

/// <summary>
/// From applying symmetry operators the function finds the normalisation for a given state
/// </summary>
cpx IsingModel_sym::get_symmetry_normalization(u64 base_idx) const {
	cpx normalisation = cpx(0.0, 0.0);
	for (int l = 0; l < this->symmetry_group.size(); l++) {
		if (this->symmetry_group[l](base_idx, this->L) == base_idx)
			normalisation += this->symmetry_eigVal[l];
	}
	return std::sqrt(normalisation);
}

/// <summary>
/// Generates the mapping to the reduced Hilbert space (reduced by symmetries: translation, spin-flip and reflection symmetry
/// adding Sz=0 (total spin) symmetry is straightforward, however, the transverse field breaks the SU(2) symmetry;
///
/// The procedure hase been successfully optimized using multithreading:
/// - each thread functions in the range [start, stop)
/// </summary>
/// <param name="start"> first index for a given thread from the original Hilbert space </param>
/// <param name="stop"> last index for a given thread from the original Hilbert space </param>
/// <param name="map_threaded"> vector containing the mapping from the reduced basis to the original Hilbert space
///                             for a given thread, the whole mapping will be merged in the generate_mapping() procedure </param>
/// <param name="_id"> identificator for a given thread </param>
void IsingModel_sym::mapping_kernel(u64 start, u64 stop, std::vector<u64>& map_threaded, std::vector<cpx>& norm_threaded, int _id){
	for (u64 j = start; j < stop; j++) {
		if (this->use_Sz_sym)
			if (__builtin_popcountll(j) != this->L / 2.) 
				continue;
		
		auto [SEC, some_value] = find_SEC_representative(j);
		if (SEC == j) {
			cpx N = get_symmetry_normalization(j);					// normalisation condition -- check if state in basis
			if (std::abs(N) > 1e-6) {
				map_threaded.push_back(j);
				norm_threaded.push_back(N);
			}
		}
	}
}

/// <summary>
/// Splits the mapping onto threads, where each finds basis states in the reduced Hilbert space within a given range.
/// The mapping is retrieved by concatenating the resulting maps from each thread
/// </summary>
void IsingModel_sym::generate_mapping() {
	u64 start = 0, stop = static_cast<u64>(std::pow(2, this->L));
	u64 two_powL = BinaryPowers[L];
	if (num_of_threads == 1)
		mapping_kernel(start, stop, this->mapping, this->normalisation, 0);
	else {
		//Threaded
		v_2d<u64> map_threaded(num_of_threads);
		v_2d<cpx> norm_threaded(num_of_threads);
		std::vector<std::thread> threads;
		threads.reserve(num_of_threads);
		for (int t = 0; t < num_of_threads; t++) {
			start = (u64)(two_powL / (double)num_of_threads * t);
			NO_OVERFLOW(stop = ((t + 1) == num_of_threads ? two_powL : u64(two_powL / (double)num_of_threads * (double)(t + 1)));)
			map_threaded[t] = v_1d<u64>();
			norm_threaded[t] = v_1d<cpx>();
			threads.emplace_back(&IsingModel_sym::mapping_kernel, this, start, stop, ref(map_threaded[t]), ref(norm_threaded[t]), t);
		}
		for (auto& t : threads) t.join();

		for (auto& t : map_threaded)
			this->mapping.insert(this->mapping.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));

		for (auto& t : norm_threaded)
			this->normalisation.insert(this->normalisation.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));
	}
	this->N = this->mapping.size();
	//assert(mapping.size() > 0 && "Not possible number of electrons - no. of states < 1");
}

// ------------------------------------------------------------------------------------------------ BUILDING HAMILTONIAN ------------------------------------------------------------------------------------------------

/// <summary>
/// Finds the representative for a given base_idx in sector_alfa normalisation potentailly from other symmetry sector beta (the same is used creating the Hamiltonian with beta=alfa)
/// </summary>
/// <returns>Representative binary number and eigenvalue from symmetries to return to that state from base_idx</returns>
std::pair<u64, cpx> find_rep_and_sym_eigval(u64 base_idx, const IsingModel_sym& sector_alfa, cpx normalisation_beta) {
	u64 idx = binary_search(sector_alfa.mapping, 0, sector_alfa.N - 1, base_idx);
	if (idx < sector_alfa.N)	return std::make_pair(idx, sector_alfa.normalisation[idx] / normalisation_beta);
	
	auto [min, sym_eig] = sector_alfa.find_SEC_representative(base_idx);
	idx = binary_search(sector_alfa.mapping, 0, sector_alfa.N - 1, min);
	if (idx < sector_alfa.N)	return std::make_pair(idx, sector_alfa.normalisation[idx] / normalisation_beta * conj(sym_eig));
	else						return std::make_pair(0, 0);
}

/// <summary>
/// Sets the non-diagonal elements of the Hamimltonian matrix with symmetry sectors: therefore the matrix elements are summed over the SEC
/// </summary>
/// <param name="k"> index of the basis state acted upon with the Hamiltonian </param>
/// <param name="value"> value of the given matrix element to be set </param>
/// <param name="temp"> resulting vector form acting with the Hamiltonian operator on the k-th basis state </param>
void IsingModel_sym::setHamiltonianElem(u64 k, double value, u64 new_idx) {
	auto [idx, sym_eig] = find_rep_and_sym_eigval(new_idx, *this, this->normalisation[k]);
	H(idx, k) += value * sym_eig;
}

/// <summary>
/// Generates the total Hamiltonian of the system. The diagonal part is straightforward,
/// while the non-diagonal terms need the specialized setHamiltonainElem(...) function
/// </summary>
void IsingModel_sym::hamiltonian() {
	try {
		this->H = arma::sp_cx_mat(this->N, this->N); //hamiltonian
		//this->H = arma::conv_to<arma::Mat<double>>::from(this->H);
	}
	catch (const bad_alloc& e) {
		std::cout << "Memory exceeded" << e.what() << "\n";
		assert(false);
	}
	#if MODEL == 0
		this->hamiltonian_Ising();
	#elif MODEL == 1
		this->hamiltonian_heisenberg();
	#elif MODEL == 3
		this->hamiltonian_xyz();
		std::cout << "Using XYZ " << std::endl;
	#else
		this->hamiltonian_Ising();
	#endif
}

void IsingModel_sym::hamiltonian_xyz(){
	std::cout << this->N << std::endl;
	double eta = this->use_Sz_sym? 0.0 : 0.5;
	std::vector<std::vector<double>> parameters = {{1.0 * (1 - eta), 1.0 * (1 + eta), this->g},
                                                    {this->J * (1 - eta), this->J * (1 + eta), this->J * this->g}
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
            std::tie(val, op_k) = operators::sigma_z(base_state, L, { j });
            this->setHamiltonianElem(k, this->h * real(val), op_k);

            std::tie(val, op_k) = operators::sigma_x(base_state, L, { j });
            this->setHamiltonianElem(k, 0.4 * real(val), op_k);
	    	
            //std::tie(val, op_k) = operators::sigma_x(base_state, L, { j });
			//this->setHamiltonianElem(k, this->g0 * real(val), op_k);
            for(int a = 0; a < neighbor_distance.size(); a++){
                int r = neighbor_distance[a];
				int nei = j + r;
				if(nei >= this->L)
					nei = (this->_BC)? -1 : nei % this->L;

	    	    if (nei >= 0) {
                    for(int b = 0; b < XYZoperators.size(); b++){
                        op_type op = XYZoperators[b];
		                auto [val1, op_k] = op(base_state, L, { j });
		                auto [val2, opop_k] = op(op_k, L, { nei });
						this->setHamiltonianElem(k, parameters[a][b] * real(val1 * val2), opop_k);
                    }
	    	    }
            }
	    }
	}
}

void IsingModel_sym::hamiltonian_Ising(){
	for (long int k = 0; k < this->N; k++) {
		double s_i;
		double s_j;
		for (int j = 0; j <= this->L - 1; j++) {
			s_i = checkBit(this->mapping[k], L - 1 - j) ? 1.0 : -1.0;              // true - spin up, false - spin down
			/* transverse field */
			if (this->g != 0) {
				NO_OVERFLOW(u64 new_idx = flip(this->mapping[k], BinaryPowers[L - 1 - j], (L - 1 - j));)
				this->setHamiltonianElem(k, this->g, new_idx);
			}
			/* disorder */
			this->H(k, k) += this->h * s_i;                                    // diagonal

			if (this->nearest_neighbors[j] >= 0) {
				/* Ising-like spin correlation */
				s_j = checkBit(this->mapping[k], this->L - 1 - this->nearest_neighbors[j]) ? 1.0 : -1.0;
				this->H(k, k) += this->J * s_i * s_j;
			}
		}
	}
}
void IsingModel_sym::hamiltonian_heisenberg() {
	for (long int k = 0; k < N; k++) {
		double s_i, s_j;
		int base_state = map(k);
		for (int j = 0; j <= L - 1; j++) {
			s_i = checkBit(base_state, L - 1 - j) ? 0.5 : -0.5;							 // true - spin up, false - spin down
			/* disorder */
			H(k, k) += this->h * s_i;                             // diagonal elements setting

			int nei = this->nearest_neighbors[j];
			if (nei >= 0) {
				s_j = checkBit(base_state, L - 1 - nei) ? 0.5 : -0.5;
				if(s_i * s_j < 0){
					u64 new_idx =  flip(base_state, BinaryPowers[this->L - 1 - nei], this->L - 1 - nei);
					new_idx =  flip(new_idx, BinaryPowers[this->L - 1 - j], this->L - 1 - j);
					setHamiltonianElem(k, 0.5 * this->J, new_idx);
				}
				/* Ising-like spin correlation */
				H(k, k) += this->g * s_i * s_j;
			}
		}
		//std::cout << std::bitset<4>(base_state) << "\t";
	}
	//std::cout << std::endl << arma::mat(this->H) << std::endl;
}

// ------------------------------------------------------------------------------------------------ PHYSICAL QUANTITTIES ------------------------------------------------------------------------------------------------


// ----------------------------------------------------------------------------- WRAPPERS FOR SIGMA OPERATORS - creating matrix -----------------------------------------------------------------------------
arma::sp_cx_mat IsingModel_sym::create_operator(std::initializer_list<op_type> operators, std::vector<int> sites) const {
	const double G = double(this->symmetry_group.size());
	//// throwables
	//for (auto& site : sites)
	//	if ((site < 0 || site >= this->L)) throw "Site index exceeds chain\n";
	arma::sp_cx_mat operator_matrix(this->N, this->N);
	for (int i = 0; i < this->N; i++) {
		const u64 base_vec = this->mapping[i];
		set_OperatorElem(operators, 1.0, sites, operator_matrix, base_vec, i);
	}
	return operator_matrix / G;
}
arma::sp_cx_mat IsingModel_sym::create_operator(std::initializer_list<op_type> operators, arma::cx_vec prefactors) const {// calculating normalisation for both sector symmetry groups
	const double G = double(this->symmetry_group.size());

	arma::cx_vec pre = prefactors.is_empty()? arma::cx_vec(this->L, arma::fill::ones) : prefactors;
	assert(pre.size() == this->L && "Input array of different size than system size!");
	
	arma::sp_cx_mat operator_matrix(this->N, this->N);
	for (int i = 0; i < this->N; i++) {
		const u64 base_vec = this->mapping[i];
		for (int j = 0; j < this->L; j++)
			set_OperatorElem(operators, pre(j), { j }, operator_matrix, base_vec, i);
	}
	return operator_matrix / (G * sqrt(this->L));
};
arma::sp_cx_mat IsingModel_sym::create_operator(std::initializer_list<op_type> operators, int corr_len, arma::cx_vec prefactors) const {
	const double G = double(this->symmetry_group.size());
	auto neis = get_neigh_vector(this->_BC, this->L, corr_len);
	arma::sp_cx_mat operator_matrix(this->N, this->N);
	
	arma::cx_vec pre = prefactors.is_empty()? arma::cx_vec(this->L, arma::fill::ones) : prefactors;
	assert(pre.size() == this->L && "Input array of different size than system size!");

	for (int i = 0; i < this->N; i++) {
		const u64 base_vec = this->mapping[i];
		for (int j = 0; j < this->L; j++) {
			const int nei = neis[j];
			if (nei >= 0) {
				set_OperatorElem(operators, pre(j), { j, nei }, operator_matrix, base_vec, i);
			}
		}
	}
	return operator_matrix / (G * sqrt(this->L));
}

void IsingModel_sym::set_OperatorElem(std::vector<op_type> operators, cpx prefactor, std::vector<int> sites, arma::sp_cx_mat& operator_matrix, u64 base_vec, u64 cur_idx) const{
	auto set_MatrixElem = [&](u64 vec_sym, cpx sym_eig, op_type op) {
		auto [op_value, vec_sym_tmp] = op(vec_sym, L, sites);
		if (abs(op_value) > 1e-14) {
			auto [idx_sym, sym_eigVal] = (vec_sym == vec_sym_tmp) ? \
				std::make_pair(cur_idx, conj(sym_eig)) :\
				find_rep_and_sym_eigval(vec_sym_tmp, *this, this->normalisation[cur_idx]);
			if (idx_sym < this->N){
				operator_matrix(idx_sym, cur_idx) += prefactor * conj(sym_eigVal) * op_value * sym_eig;
			}
		}
	};
	for (auto& op : operators) {
		for (int k = 0; k < this->symmetry_group.size(); k++) {
			auto sym_operation = this->symmetry_group[k];
			auto symRep = this->symmetry_eigVal[k];
			set_MatrixElem(sym_operation(base_vec, this->L), symRep, op);
		}
	}
}

arma::sp_cx_mat IsingModel_sym::symmetryRotation() const {
	std::vector<u64> full_map;
	if (this->use_Sz_sym)
		full_map = generate_full_map(this->L, true);
	u64 dim_tot = this->use_Sz_sym? full_map.size() : ULLPOW(this->L);
	arma::cx_vec output(dim_tot, arma::fill::zeros);

	auto find_index = [&](u64 index){
		return this->use_Sz_sym? binary_search(full_map, 0, dim_tot - 1, index) : index;
	};

	arma::sp_cx_mat U(dim_tot, this->N);

#pragma omp parallel for
	for (long int k = 0; k < this->N; k++) {
		for (int i = 0; i < this->symmetry_group.size(); i++) {
			auto new_idx = this->symmetry_group[i](this->mapping[k], this->L);
			u64 idx = find_index(new_idx);
			if(idx < dim_tot) // only if exists in sector
				U(idx, k) += conj(this->symmetry_eigVal[i] / (this->normalisation[k] * sqrt(this->symmetry_group.size())));
			// CONJUNGATE YOU MORON CAUSE YOU RETURN TO FULL STATE, I.E. INVERSE MAPPING!!!!!! 
		}
	}
	return U;
}
arma::cx_vec IsingModel_sym::symmetryRotation(const arma::cx_vec& state, std::vector<u64> full_map) const {
	if(full_map.empty() && this->use_Sz_sym)	
		full_map = generate_full_map(this->L, true);
	u64 dim_tot = this->use_Sz_sym? full_map.size() : ULLPOW(this->L);
	arma::cx_vec output(dim_tot, arma::fill::zeros);

	auto find_index = [&](u64 index){
		return this->use_Sz_sym? binary_search(full_map, 0, dim_tot - 1, index) : index;
	};

#pragma omp parallel for
	for (long int k = 0; k < this->N; k++) {
		for (int i = 0; i < this->symmetry_group.size(); i++) {
			auto new_idx = this->symmetry_group[i](this->mapping[k], this->L);
			u64 idx = find_index(new_idx);
			if(idx < dim_tot) // only if exists in sector
				output(idx) += conj(this->symmetry_eigVal[i] / (this->normalisation[k] * sqrt(this->symmetry_group.size())) * state(k));
		}
	}
	return output;
}
arma::cx_vec IsingModel_sym::symmetryRotation(u64 state_idx, std::vector<u64> full_map) const {
	if(this->use_Sz_sym || full_map.empty())	
		full_map = generate_full_map(this->L, true);
	u64 dim_tot = this->use_Sz_sym? full_map.size() : ULLPOW(this->L);
	arma::cx_vec output(dim_tot, arma::fill::zeros);

	auto find_index = [&](u64 index){
		return this->use_Sz_sym? binary_search(full_map, 0, dim_tot - 1, index) : index;
	};
#pragma omp parallel for
	for (long int k = 0; k < this->N; k++) {
		for (int i = 0; i < this->symmetry_group.size(); i++) {
			auto new_idx = this->symmetry_group[i](this->mapping[k], this->L);
			u64 idx = find_index(new_idx);
			if(idx < dim_tot) // only if exists in sector
				output(idx) += conj(this->symmetry_eigVal[i] / (this->normalisation[k] * sqrt(this->symmetry_group.size())) * this->eigenvectors(k, state_idx));
		}
	}
	return output;
}
arma::sp_cx_mat IsingModel_sym::fourierTransform(op_type op, int k) const {
	const double q = k * two_pi / double(this->L);
	arma::cx_vec k_exp(this->L);
	for (int j = 0; j < this->L; j++) {
		k_exp(j) = std::exp(im * q * double(j + 1));
	}
	return create_operator({op}, k_exp);
}
arma::sp_cx_mat IsingModel_sym::spin_current() const {
	auto beta = std::make_unique<IsingModel_disorder>(this->L, this->J, 0.0, this->g, 0.0, this->h, 0.0, this->_BC);
	auto op_full = beta->spin_current();
	auto U = this->symmetryRotation();
	return U.t() * op_full * U;
}

// ----------------------------------------------------------------------------- INTEGRABLE SOLUTIONS -----------------------------------------------------------------------------
arma::vec IsingModel_sym::get_non_interacting_energies(){
	const u64 dim = ULLPOW(this->L);
	arma::vec energies(dim);
	int counter = 0;
	double epsilon = sqrt(this->g * this->g + this->h + this->h);
	for(int k = 0; k <= L; k++)
	{
		int degeneracy = std::tgamma(this->L + 1) / ( std::tgamma(k + 1) * std::tgamma(this->L - k + 1) ); // tgamma(n) = (n-1)!
		for(int d = 0; d < degeneracy; d++)
			energies(counter++) = -(this->L - 2 * k) * epsilon;
	}
	return energies;
}

