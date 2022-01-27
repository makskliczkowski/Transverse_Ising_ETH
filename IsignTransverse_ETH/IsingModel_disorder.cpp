#include "include/IsingModel.h"

// ----------------------------------------------------------------------------- CONSTRUCTORS -----------------------------------------------------------------------------
IsingModel_disorder::IsingModel_disorder(int L, double J, double J0, double g, double g0, double h, double w, int _BC) {
	this->L = L; this->J = J; this->g = g; this->h = h;
	this->J0 = J0; this->g0 = g0;  this->w = w;
	this->N = ULLPOW(this->L);
	this->_BC = _BC;

	//this->ran.reset();
	this->reset_random();
	//change info
	this->info = "_L=" + std::to_string(this->L) + \
		",J0=" + to_string_prec(this->J0, 2) + \
		",g=" + to_string_prec(this->g, 2) + \
		",g0=" + to_string_prec(this->g0, 2) + \
		",h=" + to_string_prec(this->h, 2) + \
		",w=" + to_string_prec(this->w, 2);
	this->set_neighbors();
	this->hamiltonian();
}

// ----------------------------------------------------------------------------- BASE GENERATION AND RAPPING -----------------------------------------------------------------------------

/// <summary>
/// Return the index in the case of no mapping in disorder
/// </summary>
/// <param name="index"> index to take</param>
/// <returns>index</returns>
u64 IsingModel_disorder::map(u64 index) {
	if (index < 0 || index >= std::pow(2, L)) throw "Element out of range\n No such index in map\n";
	return index;
}

/// <summary>
/// W razie gdybyœmy robili Sz symetriê, dla picu
/// </summary>
void IsingModel_disorder::generate_mapping() {
	this->mapping = std::vector<u64>();
	std::vector<bool> base_vector(L);
	int_to_binary(u64(std::pow(2, L / 2) - 1), base_vector);
	while (next_permutation(base_vector.begin(), base_vector.end()))
		this->mapping.push_back(binary_to_int(base_vector));
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
		this->H = sp_mat(N, N);                                //  hamiltonian memory reservation
	}
	catch (const bad_alloc& e) {
		std::cout << "Memory exceeded" << e.what() << "\n";
		assert(false);
	}

	this->dh = create_random_vec(L, this->w);                               // creates random disorder vector
	this->dJ = create_random_vec(L, this->J0);                              // creates random exchange vector
	this->dg = create_random_vec(L, this->g0);                              // creates random transverse field vector
	//this->dh.zeros();
	//dh(1) = 0.165; dh(4) = -0.24;
	for (long int k = 0; k < N; k++) {
		double s_i, s_j;
		for (int j = 0; j <= L - 1; j++) {
			s_i = checkBit(k, L - 1 - j) ? 1.0 : -1.0;;							 // true - spin up, false - spin down

			NO_OVERFLOW(u64 new_idx = flip(k, BinaryPowers[this->L - 1 - j], this->L - 1 - j);)
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

// ----------------------------------------------------------------------------- PHYSICAL QUANTITES -----------------------------------------------------------------------------

/// <summary>
/// Calculates the matrix element for sigma_z Pauli matrix
/// </summary>
/// <param name="sites">Sites the matrix works on</param>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <returns>The matrix element</returns>
double IsingModel_disorder::av_sigma_z(u64 alfa, u64 beta, std::vector<int> sites) {
	for (auto& site : sites)
		if (site < 0 || site >= L) throw "Site index exceeds chain";

	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
	double value = 0;
#pragma omp parallel for reduction (+: value)
	for (int k = 0; k < N; k++) {
		double S_z = 1;
		for (auto& site : sites)
			S_z *= checkBit(k, L - 1 - site) ? 1.0 : -1.0;
		value += S_z * state_alfa(k) * state_beta(k);
	}
	return value;
}

/// <summary>
/// Calculates the matrix element for sigma_z extensive
/// </summary>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <returns>The matrix element</returns>
double IsingModel_disorder::av_sigma_z(u64 alfa, u64 beta)
{
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
	double value = 0;
#pragma omp parallel for reduction (+: value)
	for (int k = 0; k < N; k++) {
		for (int l = 0; l < this->L; l++) {
			double Sz = checkBit(k, L - 1 - l) ? 1.0 : -1.0;
			value += Sz * state_alfa(k) * state_beta(k);
		}
	}
	return value / sqrt(this->L);
}

/// <summary>
/// Calculates the matrix element for sigma_z extensive correlations
/// </summary>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <param name="corr_length">correlation length</param>
/// <returns>The matrix element</returns>
double IsingModel_disorder::av_sigma_z(u64 alfa, u64 beta, int corr_length)
{
	if (corr_length >= L) throw "exceeding correlation length\n";

	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
	double value = 0;
	auto neis = get_neigh_vector(this->_BC, this->L, corr_length);
#pragma omp parallel for reduction (+: value) collapse(2)
	for (int k = 0; k < N; k++) {
		for (int l = 0; l < this->L; l++) {
			int nei = neis[l];
			if (nei < 0) continue;
			double Sz = checkBit(k, L - 1 - l) ? 1.0 : -1.0;
			double Sz_corr = checkBit(k, L - 1 - nei) ? 1.0 : -1.0;
			value += Sz * Sz_corr * state_alfa(k) * state_beta(k);
		}
	}
	return value / sqrt(this->L);
}

/// <summary>
/// Calculates the matrix element for sigma_x Pauli matrix
/// </summary>
/// <param name="sites">Sites the matrix works on</param>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <returns>The matrix element</returns>
double IsingModel_disorder::av_sigma_x(u64 alfa, u64 beta, std::vector<int> sites) {
	for (auto& site : sites) {
		if (site < 0 || site >= L) throw "Site index exceeds chain";
	}
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
	double value = 0;
#pragma omp parallel for reduction (+: value)
	for (int k = 0; k < N; k++) {
		for (auto& site : sites) {
			NO_OVERFLOW(u64 idx = flip(k, BinaryPowers[this->L - 1 - site], this->L - 1 - site);)
			value += state_alfa(idx) * state_beta(k);
		}
	}
	return value;
}

/// <summary>
/// Calculates the matrix element for sigma_x extensive (sum over the system)
/// </summary>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <returns>The matrix element</returns>
double IsingModel_disorder::av_sigma_x(u64 alfa, u64 beta) {
	double overlap = 0;
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
#pragma omp parallel for reduction(+: overlap)
	for (long int k = 0; k < N; k++) {
		for (int j = 0; j < this->L; j++) {
			NO_OVERFLOW(u64 new_idx = flip(k, BinaryPowers[this->L - 1 - j], this->L - 1 - j);)
			overlap += state_alfa(new_idx) * state_beta(k);
		}
	}
	return overlap / sqrt(this->L);
}

/// <summary>
/// Calculates the matrix element for sigma_x extensive correlations : s^x_i s^x_i+1
/// </summary>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <param name="corr_length">correlation length</param>
/// <returns>The matrix element</returns>
double IsingModel_disorder::av_sigma_x(u64 alfa, u64 beta, int corr_length)
{
	if (corr_length >= L) throw "exceeding correlation length\n";

	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
	double value = 0;
	auto neis = get_neigh_vector(this->_BC, this->L, corr_length);
#pragma omp parallel for reduction (+: value)
	for (int k = 0; k < N; k++) {
		for (int l = 0; l < this->L; l++) {
			int nei = neis[l];
			if (nei < 0) continue;
			NO_OVERFLOW(
				u64 idx = flip(k, BinaryPowers[this->L - 1 - nei], this->L - 1 - nei);
				u64 new_idx = flip(idx, BinaryPowers[this->L - 1 - l], this->L - 1 - l);
			);
			value += state_alfa(new_idx) * state_beta(k);
		}
	}
	return value / sqrt(this->L);
}

/// <summary>
///
/// </summary>
/// <param name="alfa"></param>
/// <param name="beta"></param>
/// <param name="sites"></param>
/// <returns></returns>
double IsingModel_disorder::av_spin_flip(u64 alfa, u64 beta, std::vector<int> sites) {
	if (sites.size() != 2) throw "Not implemented such exotic operators, choose 1 or 2 sites\n";
	for (auto& site : sites) {
		if (site < 0 || site >= L) throw "Site index exceeds chain";
	}
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
	double value = 0;
#pragma omp parallel
	{
		std::vector<bool> base_vector(L), temp(L);
#pragma omp for reduction (+: value)
		for (int k = 0; k < N; k++) {
			int_to_binary(map(k), base_vector);
			auto it = sites.begin();
			auto it2 = it + 1;
			temp = base_vector;
			if ((base_vector[*it] == 0 && base_vector[*it2] == 1) || (base_vector[*it] == 1 && base_vector[*it2] == 0)) {
				temp[*it] = !base_vector[*it];
				temp[*it2] = !base_vector[*it2];
				const u64 idx = binary_to_int(temp);
				value += state_alfa(idx) * state_beta(k);
			}
		}
	}
	return 2.0 * value;
}

/// <summary>
///
/// </summary>
/// <param name="alfa"></param>
/// <param name="beta"></param>
/// <returns></returns>
double IsingModel_disorder::av_spin_flip(u64 alfa, u64 beta) {
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
	double value = 0;
#pragma omp parallel
	{
		std::vector<bool> base_vector(L), temp(L);
#pragma omp for reduction (+: value)
		for (int k = 0; k < N; k++) {
			int_to_binary(map(k), base_vector);
			for (int l = 0; l < this->L; l++) {
				temp = base_vector;
				const int nei = this->nearest_neighbors[l];
				if (nei < 0) continue;
				if ((base_vector[l] == 0 && base_vector[nei] == 1) || (base_vector[nei] == 0 && base_vector[l])) {
					temp[l] = !base_vector[l];
					temp[nei] = !base_vector[nei];
					const u64 idx = binary_to_int(temp);
					value += state_alfa(idx) * state_beta(k);
				}
			}
		}
	}
	return 2.0 * value / sqrt(L);
}

/// <summary>
///
/// </summary>
/// <param name="alfa"></param>
/// <param name="beta"></param>
/// <param name="sites"></param>
/// <returns></returns>
cpx IsingModel_disorder::av_spin_current(u64 alfa, u64 beta, std::vector<int> sites) {
	if (sites.size() != 2) throw "Not implemented such exotic operators, choose 1 or 2 sites\n";
	for (auto& site : sites) {
		if (site < 0 || site >= L) throw "Site index exceeds chain";
	}
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
	double value_real = 0, value_imag = 0;
#pragma omp parallel
	{
		std::vector<bool> base_vector(L), temp(L);
#pragma omp for reduction (+: value_real, value_imag)
		for (int k = 0; k < N; k++) {
			int_to_binary(map(k), base_vector);
			auto l = *(sites.begin());
			auto nei = *(sites.begin() + 1);
			temp = base_vector;
			if (nei < 0) continue;
			cpx value = 0.0;
			if (base_vector[l] && !base_vector[nei]) {
				temp[l] = 0;
				temp[nei] = 1;
				const u64 idx = binary_to_int(temp);
				value = state_alfa(idx) * state_beta(k) * im;
			}
			else if (!base_vector[l] && base_vector[nei]) {
				temp[l] = 1;
				temp[nei] = 0;
				const u64 idx = binary_to_int(temp);
				value = -state_alfa(idx) * state_beta(k) * im;
			}
			value_real += value.real();
			value_imag += value.imag();
		}
	}
	return 2i * cpx(value_real, value_imag);
}
/// <summary>
///
/// </summary>
/// <param name="alfa"></param>
/// <param name="beta"></param>
/// <returns></returns>
cpx IsingModel_disorder::av_spin_current(u64 alfa, u64 beta) {
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
	double value_real = 0, value_imag = 0;
#pragma omp parallel
	{
		std::vector<bool> base_vector(L), temp(L);
#pragma omp for reduction (+: value_real, value_imag)
		for (int k = 0; k < N; k++) {
			int_to_binary(map(k), base_vector);
			for (int l = 0; l < this->L; l++) {
				temp = base_vector;
				const int nei = this->nearest_neighbors[l];
				if (nei < 0) continue;
				cpx value = 0.0;
				if (base_vector[l] && !base_vector[nei]) {
					temp[l] = 0;
					temp[nei] = 1;
					const u64 idx = binary_to_int(temp);
					value = state_alfa(idx) * state_beta(k) * im;
				}
				else if (!base_vector[l] && base_vector[nei]) {
					temp[l] = 1;
					temp[nei] = 0;
					const u64 idx = binary_to_int(temp);
					value = -state_alfa(idx) * state_beta(k) * im;
				}
				value_real += value.real();
				value_imag += value.imag();
			}
		}
	}
	return 2i * cpx(value_real, value_imag) / sqrt(this->L);
}
// ----------------------------------------------------------------------------- CALCULATE MATRIX ELEMENTS VIA INPUT OPERATOR -----------------------------------------------------------------------------
cpx IsingModel_disorder::av_operator(u64 alfa, u64 beta, op_type op, std::vector<int> sites) {
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
	double overlap_real = 0, overlap_imag = 0;
#pragma omp parallel for reduction(+: overlap_real, overlap_imag)
	for (int k = 0; k < N; k++) {
		auto [ret_val, new_idx] = op(k, this->L, sites);
		cpx overlap = ret_val * state_alfa(new_idx) * state_beta(k);
		overlap_real += real(overlap);
		overlap_imag += imag(overlap);
	}
	return cpx(overlap_real, overlap_imag);
}
cpx IsingModel_disorder::av_operator(u64 alfa, u64 beta, op_type op) {
	double overlap_real = 0, overlap_imag = 0;
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
#pragma omp parallel for reduction(+: overlap_real, overlap_imag)
	for (long int k = 0; k < N; k++) {
		for (int j = 0; j < this->L; j++) {
			auto [ret_val, new_idx] = op(k, this->L, { j });
			cpx overlap = ret_val * state_alfa(new_idx) * state_beta(k);
			overlap_real += real(overlap);
			overlap_imag += imag(overlap);
		}
	}
	return cpx(overlap_real, overlap_imag) / sqrt(this->L);
}
cpx IsingModel_disorder::av_operator(u64 alfa, u64 beta, op_type op, int corr_len) {
	if (corr_len >= L) throw "exceeding correlation length\n";

	double overlap_real = 0, overlap_imag = 0;
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);
	auto neis = get_neigh_vector(this->_BC, this->L, corr_len);
#pragma omp parallel for reduction(+: overlap_real, overlap_imag)
	for (int k = 0; k < N; k++) {
		for (int l = 0; l < this->L; l++) {
			int nei = neis[l];
			if (nei < 0) continue;
			auto [ret_val, new_idx] = op(k, this->L, { l, nei });
			cpx overlap = state_alfa(new_idx) * state_beta(k);
			overlap_real += real(overlap);
			overlap_imag += imag(overlap);
		}
	}
	return cpx(overlap_real, overlap_imag) / sqrt(this->L);
}

// ----------------------------------------------------------------------------- CREATE OPERATOR TO CALCULATE MATRIX ELEMENTS -----------------------------------------------------------------------------
sp_cx_mat IsingModel_disorder::create_operator(std::initializer_list<op_type> operators, std::vector<int> sites) const {
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
sp_cx_mat IsingModel_disorder::create_operator(std::initializer_list<op_type> operators) const {
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
sp_cx_mat IsingModel_disorder::create_operator(std::initializer_list<op_type> operators, int corr_len) const {
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

sp_cx_mat IsingModel_disorder::createSq(int k) const {
	const double q = k * two_pi / double(this->L);
	arma::sp_cx_mat opMatrix(this->N, this->N);
	
	std::vector<cpx> k_exp(this->L);
	for (int j = 0; j < this->L; j++) {
		NO_OVERFLOW(k_exp[j] = std::exp(im * q * double(j + 1));)
	}
#pragma omp parallel for
	for (long int n = 0; n < N; n++) {
		for (int j = 0; j < this->L; j++) {
			auto [value, idx] = IsingModel::sigma_z(n, this->L, { j });
			opMatrix(n, n) += value * k_exp[j];
		}
	}
	return opMatrix / sqrt(this->L);
}

sp_cx_mat IsingModel_disorder::createHq(int k) const {
	const double q = k * two_pi / double(this->L);
	arma::sp_cx_mat opMatrix(this->N, this->N);
	std::vector<cpx> k_exp(this->L);
	for (int j = 0; j < this->L; j++)
		k_exp[j] = std::cos(q * double(j + 1));
#pragma omp parallel for
	for (long int n = 0; n < this->N; n++) {
		for (int j = 0; j < this->L; j++) {
			auto nei = this->nearest_neighbors[j];			
			auto [s_i, __1] = IsingModel_disorder::sigma_z(n, this->L, { j });
			opMatrix(n, n) += this->h / 2. * s_i * k_exp[j];

			auto [__2, idx1] = IsingModel_disorder::sigma_x(n, this->L, { j });
			opMatrix(idx1, n) += this->g / 2. * k_exp[j];

			if (nei >= 0) {
				auto [__3, idx2] = IsingModel_disorder::sigma_x(n, this->L, { nei });
				opMatrix(idx2, n) += this->g / 2. * k_exp[j];
				auto [s_j, __4] = IsingModel_disorder::sigma_z(n, this->L, { nei });

				opMatrix(n, n) += (this->J * s_i + this->h / 2.) * s_j * k_exp[j];
			}
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
mat IsingModel_disorder::correlation_matrix(u64 state_id) const {
	mat corr_mat(L, L, fill::zeros);
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





// ----------------------------------------------------------------------------------------- entaglement
auto IsingModel_disorder::reduced_density_matrix(const arma::cx_vec& state, int A_size) const -> arma::cx_mat {
	// set subsytsems size
	const u64 dimA = ULLPOW(A_size);
	const u64 dimB = ULLPOW(L - A_size);
	cx_mat rho(dimA, dimA, fill::zeros);
#pragma omp parallel for shared(rho)
	for (long int n = 0; n < N; n++) {						// loop over configurational basis
		if (abs(state(n)) < 1e-10) continue;				// discard non-essential terms
		long counter = 0;
		for (long m = n % dimB; m < N; m += dimB) {			// pick out state with same B side (last L-A_size bits)
			long idx = (long)std::floor(1.0 * n / dimB);	// find index of state with same B-side (by dividing the last bits are discarded)
			rho(idx, counter) += conj(state(n)) * state(m);
			counter++;										// increase counter to move along reduced basis
		}
	}
	return rho;
}
/// <summary>
/// Calculates the entropy of the system via the mixed density matrix
/// </summary>
/// <param name="state_id"> state index to produce the density matrix </param>
/// <param name="A_size"> size of subsystem </param>
/// <returns> entropy of considered systsem </returns>
double IsingModel_disorder::entaglement_entropy(const arma::cx_vec& state, int A_size) const {
	auto rho = reduced_density_matrix(state, A_size);
	vec probabilities;
	eig_sym(probabilities, rho); //diagonalize to find probabilities and calculate trace in rho's eigenbasis
	double entropy = 0;
#pragma omp parallel for reduction(+: entropy)
	for (int i = 0; i < probabilities.size(); i++) {
		auto value = probabilities(i);
		entropy += (abs(value) < 1e-10) ? 0 : -value * log2(abs(value));
	}
	//double entropy = -real(trace(rho * real(logmat(rho))));
	return entropy;
}

template <typename _ty>
inline arma::Mat<_ty> matrix_pow(const arma::Mat<_ty>& matrix, int exponent) {
	if (exponent == 1) 
		return matrix;
	else 
		return matrix * matrix_pow(matrix, exponent - 1);
}
double IsingModel_disorder::reyni_entropy(const arma::cx_vec& state, int A_size, unsigned alfa) const {
	assert(alfa > 1 && "Only alfa>=2 powers are possible");
	auto rho = reduced_density_matrix(state, A_size);
	return real(log(trace(matrix_pow(rho, alfa))) / (1.0 - alfa));
}

double IsingModel_disorder::shannon_entropy(const arma::cx_vec& state, int A_size) const {
	auto rho = reduced_density_matrix(state, A_size);
	vec probabilities;
	eig_sym(probabilities, rho); //diagonalize to find probabilities and calculate trace in rho's eigenbasis
	double entropy = 0;
#pragma omp parallel for reduction(+: entropy)
	for (int i = 0; i < probabilities.size(); i++) {
		auto value = probabilities(i) * probabilities(i);
		entropy += (abs(value) < 1e-10) ? 0 : -value * log2(abs(value));
	}
	return entropy;
}