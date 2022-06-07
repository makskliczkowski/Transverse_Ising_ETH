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
	#ifndef USE_HEISENBERG
		if(this->g == 0 && this->g0 == 0)
			generate_mapping();
	#else
		generate_mapping();
	#endif
		
	#ifdef USE_HEISENBERG
		this->hamiltonian_heisenberg();
	#else
		this->hamiltonian();
	#endif
}

// ----------------------------------------------------------------------------- BASE GENERATION AND RAPPING -----------------------------------------------------------------------------

/// <summary>
/// Return the index in the case of no mapping in disorder
/// </summary>
/// <param name="index"> index to take</param>
/// <returns>index</returns>
u64 IsingModel_disorder::map(u64 index) const {
	if (index < 0 || index >= std::pow(2, L)) throw "Element out of range\n No such index in map\n";
	#ifndef USE_HEISENBERG
		return (this->g == 0? mapping[index] : index);
	#else
		return mapping[index];
	#endif
}

/// <summary>
/// W razie gdyby�my robili Sz symetri�, dla picu
/// </summary>
void IsingModel_disorder::generate_mapping() {
	this->mapping = std::vector<u64>();
	for (u64 j = 0; j < (ULLPOW(this->L)); j++)
		if (__builtin_popcountll(j) != this->L / 2.)
			this->mapping.push_back(j);
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
	H(new_idx, k) += value;
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
		this->dh = create_random_vec(L, this->w);                               // creates random disorder vector
		this->dJ = create_random_vec(L, this->J0);                              // creates random exchange vector
		this->dg = create_random_vec(L, this->g0);                              // creates random transverse field vector
	}
	//this->dh.zeros();
	//dh(1) = 0.165; dh(4) = -0.24;
	for (long int k = 0; k < N; k++) {
		double s_i, s_j;
		for (int j = 0; j <= L - 1; j++) {
			s_i = checkBit(k, L - 1 - j) ? 1.0 : -1.0;;							 // true - spin up, false - spin down

			u64 new_idx = flip(k, BinaryPowers[this->L - 1 - j], this->L - 1 - j);
			setHamiltonianElem(k, this->g + this->dg(j), new_idx);

			/* disorder */
			H(k, k) += (this->h + dh(j)) * s_i;                             // diagonal elements setting

			if (nearest_neighbors[j] >= 0) {
				/* Ising-like spin correlation */
				s_j = checkBit(k, this->L - 1 - nearest_neighbors[j]) ? 1.0 : -1.0;
				this->H(k, k) += (this->J + this->dJ(j)) * s_i * s_j;		// setting the neighbors elements
			}
		}
	}
}

void IsingModel_disorder::hamiltonian_heisenberg(){
	this->H = arma::sp_mat(N, N); 

	this->dh = create_random_vec(L, this->w);                               // creates random disorder vector
	this->dJ = create_random_vec(L, this->J0);                              // creates random exchange vector
	this->dg = create_random_vec(L, this->g0);                              // creates random transverse field vector
	for (long int kk = 0; kk < N; kk++) {
		double s_i, s_j;
		int k = map(k);
		for (int j = 0; j <= L - 1; j++) {
			s_i = checkBit(k, L - 1 - j) ? 0.5 : -0.5;							 // true - spin up, false - spin down
			/* disorder */
			H(kk, kk) += (this->h + this->dh(j)) * s_i;                             // diagonal elements setting

			if (nearest_neighbors[j] >= 0) {
				auto [value_x, new_idx_x] = IsingModel_disorder::sigma_x(k, this->L, {j, nearest_neighbors[j]});
				setHamiltonianElem(kk, 0.25 * (this->J + this->dJ(j)), new_idx_x);
				
				auto [value_y, new_idx_y] = IsingModel_disorder::sigma_y(k, this->L, {j, nearest_neighbors[j]});
				setHamiltonianElem(kk, 0.25 * real(value_y) * (this->J + this->dJ(j)), new_idx_y);

				/* Ising-like spin correlation */
				auto [value_z, new_idx_z] = IsingModel_disorder::sigma_z(k, this->L, {j, nearest_neighbors[j]});
				setHamiltonianElem(kk, 0.25 * real(value_z) * (this->g + this->dg(j)), new_idx_z);
			}
		}
	}
}
// ----------------------------------------------------------------------------- PHYSICAL QUANTITES -----------------------------------------------------------------------------

// ----------------------------------------------------------------------------- CREATE OPERATOR TO CALCULATE MATRIX ELEMENTS -----------------------------------------------------------------------------
arma::sp_cx_mat IsingModel_disorder::create_operator(std::initializer_list<op_type> operators, std::vector<int> sites) const {
	arma::sp_cx_mat opMatrix(this->N, this->N);
#pragma omp parallel for
	for (long int k = 0; k < N; k++) {
		for (auto& op : operators) {
			cpx value; u64 idx;
			std::tie(value, idx) = op(k, this->L, sites);
#pragma omp critical
			opMatrix(idx, k) += value;
		}
	}
	return opMatrix;
}
arma::sp_cx_mat IsingModel_disorder::create_operator(std::initializer_list<op_type> operators) const {
	arma::sp_cx_mat opMatrix(this->N, this->N);
#pragma omp parallel for
	for (long int k = 0; k < N; k++) {
		for (int j = 0; j < this->L; j++) {
			for (auto& op : operators) {
				cpx value; u64 idx;
				std::tie(value, idx) = op(k, this->L, { j });
#pragma omp critical
				opMatrix(idx, k) += value;
			}
		}
	}
	return opMatrix / sqrt(this->L);
}
arma::sp_cx_mat IsingModel_disorder::create_operator(std::initializer_list<op_type> operators, int corr_len) const {
	arma::sp_cx_mat opMatrix(this->N, this->N);
	auto neis = get_neigh_vector(this->_BC, this->L, corr_len);
#pragma omp parallel for
	for (long int k = 0; k < N; k++) {
		for (int j = 0; j < this->L; j++) {
			const int nei = neis[j];
			if (nei < 0) continue;
			for (auto& op : operators) {
				cpx value; u64 idx;
				std::tie(value, idx) = op(k, this->L, { j, nei });
#pragma omp critical
				opMatrix(idx, k) += value;
			}
		}
	}
	return opMatrix / sqrt(this->L);
}

arma::sp_cx_mat IsingModel_disorder::fourierTransform(op_type op, int k) const {
	const double q = k * two_pi / double(this->L);
	arma::sp_cx_mat opMatrix(this->N, this->N);

	std::vector<cpx> k_exp(this->L);
	for (int j = 0; j < this->L; j++) {
		NO_OVERFLOW(k_exp[j] = std::exp(im * q * double(j + 1));)
	}
#pragma omp parallel for
	for (long int n = 0; n < N; n++) {
		cpx value = 0.0;
		u64 idx = 0;
		for (int j = 0; j < this->L; j++) {
			std::tie(value, idx) = op(n, this->L, { j });
#pragma omp critical
			{
				opMatrix(idx, n) += value * k_exp[j];
			}
		}
	}
	return opMatrix / sqrt(this->L);
}
arma::sp_cx_mat IsingModel_disorder::createHq(int k) const {
	const double q = k * two_pi / double(this->L);
	arma::sp_cx_mat opMatrix(this->N, this->N);
	std::vector<cpx> k_exp(this->L);
	for (int j = 0; j < this->L; j++)
		k_exp[j] = std::cos(q * double(j + 1));
#pragma omp parallel for
	for (long int n = 0; n < this->N; n++) {
		cpx __2 = 0.0, __3 = 0.0;
		u64 idx1 = 0, idx2 = 0;
		for (int j = 0; j < this->L; j++) {
			auto nei = this->nearest_neighbors[j];			
			auto [s_i, __1] = IsingModel_disorder::sigma_z(n, this->L, { j });
			opMatrix(n, n) += this->h / 2. * s_i * k_exp[j];

			std::tie(__2, idx1) = IsingModel_disorder::sigma_x(n, this->L, { j });
#pragma omp critical
			{
				opMatrix(idx1, n) += this->g / 2. * k_exp[j];
			}
			if (nei >= 0) {
				std::tie(__3, idx2) = IsingModel_disorder::sigma_x(n, this->L, { nei });
#pragma omp critical
				{
					opMatrix(idx2, n) += this->g / 2. * k_exp[j];
				}
				auto [s_j, __4] = IsingModel_disorder::sigma_z(n, this->L, { nei });

				opMatrix(n, n) += (this->J * s_i + this->h / 2.) * s_j * k_exp[j];
			}
		}
	}
	return opMatrix / std::sqrt(this->L);
}
arma::sp_cx_mat IsingModel_disorder::createHlocal(int k) const {
	arma::sp_cx_mat opMatrix(this->N, this->N);
#pragma omp parallel for
	for (long int n = 0; n < this->N; n++) {
		cpx __2 = 0.0, __3 = 0.0;
		u64 idx1 = 0, idx2 = 0;
		auto nei = this->nearest_neighbors[k];
		auto [s_i, __1] = IsingModel_disorder::sigma_z(n, this->L, { k });
		opMatrix(n, n) += this->h / 2. * s_i;

		std::tie(__2, idx1) = IsingModel_disorder::sigma_x(n, this->L, { k });
#pragma omp critical
		{
			opMatrix(idx1, n) += this->g / 2.;
		}
		if (nei >= 0) {
			std::tie(__3, idx2) = IsingModel_disorder::sigma_x(n, this->L, { nei });
#pragma omp critical
			{
				opMatrix(idx2, n) += this->g / 2.;
			}
			auto [s_j, __4] = IsingModel_disorder::sigma_z(n, this->L, { nei });

			opMatrix(n, n) += (this->J * s_i + this->h / 2.) * s_j;
		}
	}
	return opMatrix / std::sqrt(this->L);
}

// ----------------------------------------------------------------------------- TO REFACTOR AND CREATE DESCRIPTION -----------------------------------------------------------------------------

/// <summary>
/// Calculates the spin correlation matrix within a given state (non-equilibrium average)
/// </summary>
/// <param name="state_id"> index of given state </param>
/// <returns> correlation matrix </returns>
arma::mat IsingModel_disorder::correlation_matrix(u64 state_id) const {
	arma::mat corr_mat(L, L, arma::fill::zeros);
	const arma::subview_col state = this->eigenvectors.col(state_id);
#pragma omp parallel shared(corr_mat)
	{
		u64 idx;
		vector<bool> vect(L), temp(L);
#pragma omp parallel
		for (int p = 0; p < N; p++) {
			int_to_binary(p, vect);
			for (int m = 0; m < L; m++) {
				double Szm = 0;
				double S2_tmp = 0;

				if (vect[m] == 1) Szm = 0.5;
				else Szm = -0.5;

				S2_tmp = (Szm * Szm + 0.5) * state(p) * state(p); //  <Sz(m) Sz(m)> + 2*( <S(m)^+ S(m)^-> + <S(m)^- S(m)^+> ) on-site
				corr_mat(m, m) += S2_tmp;
				for (int k = m + 1; k < L; k++) {
					double Szk = 0;
					if (vect[k] == 1) Szk = 0.5;
					else Szk = -0.5;
					S2_tmp = Szm * Szk * state(p) * state(p); //  <Sz(m) Sz(k)> for k > m

					// <S^+ S^->
					temp = vect;
					if (vect[m] == 1 && vect[k] == 0) {
						temp[m] = 0;
						temp[k] = 1;
						idx = binary_to_int(temp);
						S2_tmp += 0.5 * state(idx) * state(p);
					}
					//<S^- S^+>
					else if (vect[m] == 0 && vect[k] == 1) {
						temp[m] = 1;
						temp[k] = 0;
						idx = binary_to_int(temp);
						S2_tmp += 0.5 * state(idx) * state(p);
					}
					corr_mat(m, k) += S2_tmp;
					corr_mat(k, m) += S2_tmp;
				}
			}
		}
	}
	return corr_mat;
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

// ----------------------------------------------------------------------------- ENTAGLEMENT -----------------------------------------------------------------------------
auto IsingModel_disorder::reduced_density_matrix(const arma::cx_vec& state, int A_size) const -> arma::cx_mat {
	// set subsytsems size
	const long long dimA = ULLPOW(A_size);
	const long long dimB = ULLPOW((L - A_size));
	arma::cx_mat rho(dimA, dimA, arma::fill::zeros);
	for (long long n = 0; n < this->N; n++) {						// loop over configurational basis
		long long counter = 0;
		for (long long m = n % dimB; m < this->N; m += dimB) {			// pick out state with same B side (last L-A_size bits)
			long idx = n / dimB;							// find index of state with same B-side (by dividing the last bits are discarded)
			rho(idx, counter) += conj(state(n)) * state(m);
			counter++;										// increase counter to move along reduced basis
		}
	}
	return rho;	
}
