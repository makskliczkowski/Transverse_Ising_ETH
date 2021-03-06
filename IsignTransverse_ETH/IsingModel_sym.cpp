#include "include/IsingModel.h"

// ------------------------------------------------------------------------------------------------ CONSTRUCTORS ------------------------------------------------------------------------------------------------
IsingModel_sym::IsingModel_sym(int L, double J, double g, double h, int k_sym, bool p_sym, bool x_sym, int _BC) {
	this->L = L; this->J = J; this->g = g; this->h = h; this->_BC = _BC;
	symmetries.k_sym = k_sym * two_pi / double(this->L);
	symmetries.p_sym = (p_sym) ? 1 : -1;
	symmetries.x_sym = (x_sym) ? 1 : -1;
	k_sector = abs(this->symmetries.k_sym) < 1e-4 || abs(this->symmetries.k_sym - pi) < 1e-4;
	// precalculate the exponents
	this->k_exponents = v_1d<cpx>(this->L, 0.0);
#pragma omp parallel for
	for (int l = 0; l < this->L; l++)
		this->k_exponents[l] = std::exp(-im * this->symmetries.k_sym * double(l));
	this->createSymmetryGroup();

	this->info = "_L=" + std::to_string(this->L) + \
		",J=" + to_string_prec(this->J, 2) + \
		",g=" + to_string_prec(this->g, 2) + \
		",h=" + to_string_prec(this->h, 2) + \
		",k=" + std::to_string(k_sym) + \
		",p=" + std::to_string(symmetries.p_sym) + \
		",x=" + std::to_string(symmetries.x_sym);

	this->mapping = std::vector<u64>();
	this->normalisation = std::vector<cpx>();
	this->generate_mapping();
	if (this->N <= 0) return;
	this->set_neighbors(); // generate neighbors
	#ifdef USE_HEISENBERG
		this->hamiltonian_heisenberg();
	#else
		this->hamiltonian();
	#endif
}

// ------------------------------------------------------------------------------------------------ BASE GENERATION, SEC CLASSES AND RAPPING ------------------------------------------------------------------------------------------------
void IsingModel_sym::createSymmetryGroup() {
	this->symmetry_group = std::vector<std::function<u64(u64, int)>>();
	this->symmetry_eigVal = std::vector<cpx>();
	std::function<u64(u64, int)> e = [](u64 n, int L)->u64 {return n; };		// neutral element
	std::function<u64(u64, int)> T = e;											// in 1st iteration is neutral element
	std::function<u64(u64, int)> Z = static_cast<u64(*)(u64, int)>(&flip);		// spin flip operator (all spins)
	std::function<u64(u64, int)> P = reverseBits;								// parity operator
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
		if (this->g == 0 && __builtin_popcountll(j) != this->L / 2.) continue;
		auto [SEC, some_value] = find_SEC_representative(j);
		if (SEC == j) {
			cpx N = get_symmetry_normalization(j);					// normalisation condition -- check wether state in basis
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
	#ifdef USE_HEISENBERG
		double g_temp = this->g;
		this->g = 0;
	#endif
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
	#ifdef USE_HEISENBERG
		this->g = g_temp;
	#endif
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
	this->H = arma::sp_cx_mat(this->N, this->N); //hamiltonian
	for (long int kk = 0; kk < this->N; kk++) {
		double s_i;
		double s_j;
		int k = this->mapping[kk];
		for (int j = 0; j <= this->L - 1; j++) { 
			s_i = checkBit(k, L - 1 - j) ? 1.0 : -1.0;	// true - spin up, false - spin down
			
			/* longitudal field */
			H(kk, kk) += (this->h) * s_i;                             // diagonal elements setting

			if (nearest_neighbors[j] >= 0) {
				auto [value_x, new_idx_x] = IsingModel_disorder::sigma_x(k, this->L, {j, nearest_neighbors[j]});
				setHamiltonianElem(kk, this->J, new_idx_x);
				
				auto [value_y, new_idx_y] = IsingModel_disorder::sigma_y(k, this->L, {j, nearest_neighbors[j]});
				setHamiltonianElem(kk, real(value_y) * (this->J), new_idx_y);

				/* Ising-like spin correlation */
				auto [value_z, new_idx_z] = IsingModel_disorder::sigma_z(k, this->L, {j, nearest_neighbors[j]});
				setHamiltonianElem(kk, real(value_z) * (this->g), new_idx_z);
			}
		}
	}
}

// ------------------------------------------------------------------------------------------------ PHYSICAL QUANTITTIES ------------------------------------------------------------------------------------------------

/// <summary>
/// <summary>
/// Calculates the matrix element for sigma_z Pauli matrix
/// </summary>
/// <param name="sites">Sites the matrix works on</param>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <returns>The matrix element</returns>
double IsingModel_sym::av_sigma_z(u64 alfa, u64 beta, std::vector<int> sites) {
	auto sig_z = IsingModel::sigma_z;
	return real(av_operator(alfa, beta, *this, *this, sig_z, sites));
}

/// <summary>
/// Calculates the matrix element for sigma_z extensive
/// </summary>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <returns>The matrix element</returns>
double IsingModel_sym::av_sigma_z(u64 alfa, u64 beta) {
	auto sig_z = IsingModel::sigma_z;
	return real(av_operator(alfa, beta, *this, *this, sig_z));
}

/// <summary>
/// Calculates the matrix element for sigma_z extensive correlations
/// </summary>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <param name="corr_length">correlation length</param>
/// <returns>The matrix element</returns>
double IsingModel_sym::av_sigma_z(u64 alfa, u64 beta, int corr_length) {
	auto sig_z = IsingModel::sigma_z;
	return real(av_operator(alfa, beta, *this, *this, sig_z, corr_length));
}

/// <summary>
/// Calculates the matrix element for sigma_x Pauli matrix
/// </summary>
/// <param name="sites">Sites the matrix works on</param>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <returns>The matrix element</returns>
double IsingModel_sym::av_sigma_x(u64 alfa, u64 beta, std::vector<int> sites) {
	auto sig_x = IsingModel::sigma_x;
	return real(av_operator(alfa, beta, *this, *this, sig_x, sites));
}

/// <summary>
/// Calculates the matrix element for sigma_x extensive
/// </summary>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <returns>The matrix element</returns>
double IsingModel_sym::av_sigma_x(const u64 alfa, const u64 beta) {
	auto sig_x = IsingModel::sigma_x;
	return real(av_operator(alfa, beta, *this, *this, sig_x));
}

/// <summary>
/// Calculates the matrix element for sigma_x extensive correlations
/// </summary>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <param name="corr_length">correlation length</param>
/// <returns>The matrix element</returns>
double IsingModel_sym::av_sigma_x(u64 alfa, u64 beta, int corr_length) {
	auto sig_x = IsingModel::sigma_x;
	return real(av_operator(alfa, beta, *this, *this, sig_x, corr_length));
}

/// <summary>
///
/// </summary>
/// <param name="alfa"></param>
/// <param name="beta"></param>
/// <returns></returns>
double IsingModel_sym::av_spin_flip(u64 alfa, u64 beta) {
	auto spin_flip = IsingModel::spin_flip;
	auto value = av_operator(alfa, beta, *this, *this, spin_flip);
	value += conj(av_operator(beta, alfa, *this, *this, spin_flip));
	return 0.5 * real(value);
}

/// <summary>
///
/// </summary>
/// <param name="alfa"></param>
/// <param name="beta"></param>
/// <param name="sites"></param>
/// <returns></returns>
double IsingModel_sym::av_spin_flip(u64 alfa, u64 beta, std::vector<int> sites) {
	auto spin_flip = IsingModel::spin_flip;
	auto value = av_operator(alfa, beta, *this, *this, spin_flip, sites);
	value += conj(av_operator(beta, alfa, *this, *this, spin_flip, sites));
	return 0.5 * real(value);
}

/// <summary>
///
/// </summary>
/// <param name="alfa"></param>
/// <param name="beta"></param>
/// <returns></returns>
cpx IsingModel_sym::av_spin_current(u64 alfa, u64 beta) {
	auto spin_flip = IsingModel::spin_flip;
	auto value = im * av_operator(alfa, beta, *this, *this, spin_flip, 1);
	value -= conj(im * av_operator(beta, alfa, *this, *this, spin_flip, 1));
	return 0.5i * value;
}

/// <summary>
///
/// </summary>
/// <param name="alfa"></param>
/// <param name="beta"></param>
/// <param name="sites"></param>
/// <returns></returns>
cpx IsingModel_sym::av_spin_current(u64 alfa, u64 beta, std::vector<int> sites) {
	auto spin_flip = IsingModel::spin_flip;
	auto value = im * av_operator(alfa, beta, *this, *this, spin_flip, sites);
	value += conj(im * av_operator(beta, alfa, *this, *this, spin_flip, sites));
	return 0.5i * value;
}

/// <summary>
/// calculates spin-spin correlation matrix within a given state
/// </summary>
/// <param name="state_idx"> index of eigenstate to calculate correlation matrix </param>
/// <returns> correlation matrix </returns>
arma::mat IsingModel_sym::correlation_matrix(u64 state_idx) const {
	arma::mat corr_mat(this->L, this->L, arma::fill::zeros);
	arma::sp_cx_mat SigmaZOp = 0.5 * this->create_operator({ IsingModel_sym::sigma_z });
	for (int i = 0; i < L; i++) {
		for (int j = i; j < L; j++) {
			auto kinetic = 0.5 * arma::cdot(this->eigenvectors.col(state_idx), this->create_operator({ IsingModel_sym::spin_flip }, { i, j }) * this->eigenvectors.col(state_idx));
			kinetic += conj(kinetic);
			auto corr_zz = arma::cdot(this->eigenvectors.col(state_idx), this->create_operator({ IsingModel_sym::sigma_z }, { i, j }) * this->eigenvectors.col(state_idx));
			corr_mat(i, j) = real(kinetic + corr_zz);
		}
		//for (int j = 0; j < i; j++) {
		//	corr_mat(i, j) = corr_mat(j, i);
		//}
	}
	return corr_mat;
}

// ----------------------------------------------------------------------------- WRAPPERS FOR SIGMA OPERATORS - getting overlap -----------------------------------------------------------------------------

/// <summary>
/// This form of average operators takes potentailly two symmetry sectors alfa and beta with
/// according states and tries to find a matrix element <alfa|op|beta>. Additionally it takes the
/// list of sites that the operator works on as a multiplication (f.e. op_0 * op_1 |beta> with {0, 1}
/// </summary>
cpx av_operator(u64 alfa, u64 beta, const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta, op_type op, std::vector<int> sites) {
	// throwables
	for (auto& site : sites)
		if ((site < 0 || site >= sec_alfa.L) && sec_alfa.L != sec_beta.L) throw "Site index exceeds chain or incompatible chain lengths. Your choice\n";
	if (sec_alfa.J != sec_beta.J || sec_alfa.h != sec_beta.h || sec_alfa.g != sec_beta.g) throw "incompatible model parameters, sun \\('.')// \n";

	arma::subview_col state_alfa = sec_alfa.eigenvectors.col(alfa);
	arma::subview_col state_beta = sec_beta.eigenvectors.col(beta);

	// calculating normalisation for both sector symmetry groups
	double G = 0;
	double G_alfa = (double)sec_alfa.get_sym_group().size();
	double G_beta = (double)sec_beta.get_sym_group().size();
	G = std::sqrt(G_alfa * G_beta);

	// going through all sector beta states
	double overlap_real = 0;
	double overlap_imag = 0;
#pragma omp parallel for reduction(+:overlap_real, overlap_imag)
	for (int i = 0; i < sec_beta.N; i++) {
		const u64 base_vec = sec_beta.mapping[i];
		cpx overlap = apply_sym_overlap(state_alfa, state_beta, base_vec, i, sec_alfa, sec_beta, op, sites);
		overlap_real += overlap.real();
		overlap_imag += overlap.imag();
	}
	return cpx(overlap_real, overlap_imag) / G;
}

/// <summary>
/// This form of average operators takes potentailly two symmetry sectors alfa and beta with
/// according states and tries to find a matrix element <alfa|op|beta>. Without taking states it calculates
/// <alfa| \sum _i = 0 ^N op_i |beta>
/// </summary>
cpx av_operator(u64 alfa, u64 beta, const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta, op_type op) {
	if (sec_alfa.L != sec_beta.L) throw "Incompatible chain lengths. Your choice\n";
	if (sec_alfa.J != sec_beta.J || sec_alfa.h != sec_beta.h || sec_alfa.g != sec_beta.g) throw "incompatible model parameters, sun \\('.')// \n";

	arma::subview_col state_alfa = sec_alfa.eigenvectors.col(alfa);
	arma::subview_col state_beta = sec_beta.eigenvectors.col(beta);

	// calculating normalisation for both sector symmetry groups
	double G = 0;
	double G_alfa = (double)sec_alfa.get_sym_group().size();
	double G_beta = (double)sec_beta.get_sym_group().size();
	G = std::sqrt(G_alfa * G_beta);

	// going through all sector beta states
	double overlap_real = 0;
	double overlap_imag = 0;
#pragma omp parallel for reduction(+:overlap_real, overlap_imag)
	for (int i = 0; i < sec_beta.N; i++) {
		for (int l = 0; l < sec_beta.L; l++) {
			const u64 base_vec = sec_beta.mapping[i];
			cpx overlap = apply_sym_overlap(state_alfa, state_beta, base_vec, i, sec_alfa, sec_beta, op, { l });
			overlap_real += overlap.real();
			overlap_imag += overlap.imag();
		}
	}
	return cpx(overlap_real, overlap_imag) / (G * sqrt(sec_alfa.L));
}

/// <summary>
/// This form of average operators takes potentailly two symmetry sectors alfa and beta with
/// according states and tries to find a matrix element <alfa|op|beta>. Without taking states it calculates
/// <alfa| \sum _i = 0 ^N op_i op_{i + corr_len} |beta>
/// </summary>
cpx av_operator(u64 alfa, u64 beta, const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta, op_type op, int corr_len) {
	if (sec_alfa.L != sec_beta.L) throw "Incompatible chain lengths. Your choice\n";
	if (corr_len >= sec_alfa.L) throw "Exceeding correlation length\n";
	if (sec_alfa.J != sec_beta.J || sec_alfa.h != sec_beta.h || sec_alfa.g != sec_beta.g) throw "incompatible model parameters, sun \\('.')// \n";

	arma::subview_col state_alfa = sec_alfa.eigenvectors.col(alfa);
	arma::subview_col state_beta = sec_beta.eigenvectors.col(beta);

	// calculating normalisation for both sector symmetry groups
	double G = 0;
	double G_alfa = (double)sec_alfa.get_sym_group().size();
	double G_beta = (double)sec_beta.get_sym_group().size();
	G = std::sqrt(G_alfa * G_beta);

	auto neis = get_neigh_vector(sec_alfa._BC, sec_alfa.L, corr_len);
	// going through all sector beta states

	double overlap_real = 0;
	double overlap_imag = 0;
#pragma omp parallel for reduction(+:overlap_real, overlap_imag)
	for (int i = 0; i < sec_beta.N; i++) {
		for (int l = 0; l < sec_beta.L; l++) {
			const u64 base_vec = sec_beta.mapping[i];
			const int nei = neis[l];
			if (nei >= 0) {
				cpx overlap = apply_sym_overlap(state_alfa, state_beta, base_vec, i, sec_alfa, sec_beta, op, { l, nei });
				overlap_real += overlap.real();
				overlap_imag += overlap.imag();
			}
		}
	}
	return cpx(overlap_real, overlap_imag) / (G * sqrt(sec_alfa.L));
}

/// <summary>
/// It applies the symmetry overlap for given base vector from sector beta
/// </summary>
/// <param name="alfa">sector alfa eigenstate (left) </param>
/// <param name="beta">sector beta eigenstate (right) </param>
/// <param name="base_vec">base vector in numerical representation</param>
/// <param name="k">index in mapping of beta sector</param>
/// <param name="op">Operator that shall be applied on beta state</param>
/// <param name="sites">sites that stand for multiplication of operators</param>
/// <returns> overlap of <alfa| coeff(k) |base_vec> </returns>
cpx apply_sym_overlap(const arma::subview_col<cpx>& alfa, const arma::subview_col<cpx>& beta, u64 base_vec, u64 k, const IsingModel_sym& sec_alfa, \
	const IsingModel_sym& sec_beta, op_type op, std::vector<int> sites)
{
	const int L = sec_alfa.L;
	auto get_overlap_sym = [&](u64 vec_sym, cpx sym_eig) {
		cpx overlap = 0.0;
		auto [val_sym_beta, vec_sym_tmp] = op(vec_sym, L, sites);
		if (abs(val_sym_beta) > 1e-14) {
			auto [idx_sym, val_sym_alfa] = (vec_sym == vec_sym_tmp) ? \
				std::make_pair(k, sym_eig) : \
				find_rep_and_sym_eigval(vec_sym_tmp, sec_alfa, sec_beta.normalisation[k]);
			if (idx_sym < sec_alfa.N)
				overlap = sym_eig * conj(val_sym_alfa * alfa(idx_sym)) * beta(k) * val_sym_beta;
		}
		return overlap;
	};
	cpx overlap = 0.0;
	auto symmetry_group = sec_beta.get_sym_group();
	auto symmetry_eigVal = sec_beta.get_sym_eigVal();
	for (int k = 0; k < symmetry_group.size(); k++) {
		auto sym_operation	= symmetry_group[k];
		auto symRepr		= symmetry_eigVal[k];
		overlap += get_overlap_sym(sym_operation(base_vec, L), symRepr);
	}
	return overlap;
}


//--------------------------------------------------------- dummy functions
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
arma::vec IsingModel_sym::first_interacting_correction(){
	const u64 dim = ULLPOW(this->L);
	arma::vec energies(dim);
	int counter = 0;
	const double epsilon = sqrt(this->g * this->g + this->h + this->h);
	const double lambda = this->h * this->h / (epsilon * epsilon);
	for(int k = 0; k <= L; k++)
	{
		std::vector<int> bitmask(k, 1); 	// string with k-leading 1's
		bitmask.resize(this->L, 0);			// L-k trailing 0's

		// -------- permute all binary representations to get all combinations of subset of size k
    	do {
			double E = 0.0;
    	    for (int i = 0; i < this->L; i++)  {
				const int nei = this->nearest_neighbors[i];
				if(nei >= 0){
    	        	if (bitmask[i] == bitmask[nei]) E += this->J * lambda;
					else 							E -= this->J * lambda;
				}
    	    }
			energies(counter++) = E;
    	} while (std::prev_permutation(bitmask.begin(), bitmask.end()));
	}
	auto E = this->get_non_interacting_energies();
	return E + energies;
}

// ----------------------------------------------------------------------------------------- entaglement
auto IsingModel_sym::reduced_density_matrix(const arma::cx_vec& state, int A_size) const -> arma::cx_mat {
	// set subsytsems size
	const long long dimA    = ULLPOW(A_size);
	const long long dimB    = ULLPOW((this->L - A_size));
	const long long dim_tot = ULLPOW(this->L);
	arma::cx_mat rho(dimA, dimA, arma::fill::zeros);
	const arma::cx_vec state_full_hilbert = this->symmetryRotation(state);
	for (long long n = 0; n < dim_tot; n++) {							// loop over whole configurational basis
		long long counter = 0;
		for (long long m = n % dimB; m < dim_tot; m += dimB) {			// pick out state with same B side (last L-A_size bits)
			long idx = n / dimB;										// find index of state with same B-side (by dividing the last bits are discarded)
			rho(idx, counter) += conj(state_full_hilbert(n)) * state_full_hilbert(m);
			counter++;													// increase counter to move along reduced basis
		}
	}
	return rho;
}



// ----------------------------------------------------------------------------- WRAPPERS FOR SIGMA OPERATORS - creating matrix -----------------------------------------------------------------------------
arma::sp_cx_mat IsingModel_sym::create_operator(std::initializer_list<op_type> operators, std::vector<int> sites) const {
	const double G = double(this->symmetry_group.size());
	//// throwables
	//for (auto& site : sites)
	//	if ((site < 0 || site >= this->L)) throw "Site index exceeds chain\n";
	arma::sp_cx_mat operator_matrix(this->N, this->N);
	for (int i = 0; i < this->N; i++) {
		const u64 base_vec = this->mapping[i];
		set_OperatorElem(operators, sites, operator_matrix, base_vec, i);
	}
	return operator_matrix / G;
}
arma::sp_cx_mat IsingModel_sym::create_operator(std::initializer_list<op_type> operators) const {// calculating normalisation for both sector symmetry groups
	const double G = double(this->symmetry_group.size());
	arma::sp_cx_mat operator_matrix(this->N, this->N);
	for (int i = 0; i < this->N; i++) {
		const u64 base_vec = this->mapping[i];
		for (int j = 0; j < this->L; j++)
			set_OperatorElem(operators, { j }, operator_matrix, base_vec, i);
	}
	return operator_matrix / (G * sqrt(this->L));
};
arma::sp_cx_mat IsingModel_sym::create_operator(std::initializer_list<op_type> operators, int corr_len) const {
	const double G = double(this->symmetry_group.size());
	auto neis = get_neigh_vector(this->_BC, this->L, corr_len);
	arma::sp_cx_mat operator_matrix(this->N, this->N);
	for (int i = 0; i < this->N; i++) {
		const u64 base_vec = this->mapping[i];
		for (int j = 0; j < this->L; j++) {
			const int nei = neis[j];
			if (nei >= 0) {
				set_OperatorElem(operators, { j, nei }, operator_matrix, base_vec, i);
			}
		}
	}
	return operator_matrix / (G * sqrt(this->L));
}

void IsingModel_sym::set_OperatorElem(std::vector<op_type> operators, std::vector<int> sites, arma::sp_cx_mat& operator_matrix, u64 base_vec, u64 cur_idx) const{
	auto set_MatrixElem = [&](u64 vec_sym, cpx sym_eig, op_type op) {
		auto [op_value, vec_sym_tmp] = op(vec_sym, L, sites);
		if (abs(op_value) > 1e-14) {
			auto [idx_sym, sym_eigVal] = (vec_sym == vec_sym_tmp) ? \
				std::make_pair(cur_idx, conj(sym_eig)) :\
				find_rep_and_sym_eigval(vec_sym_tmp, *this, this->normalisation[cur_idx]);
			if (idx_sym < this->N)
				operator_matrix(idx_sym, cur_idx) += conj(sym_eigVal) * op_value * sym_eig;
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
	u64 dim_tot = ULLPOW(this->L);
	arma::sp_cx_mat U(dim_tot, this->N);
#pragma omp parallel for
	for (long int k = 0; k < this->N; k++) {
		for (int i = 0; i < this->symmetry_group.size(); i++) {
			auto idx = this->symmetry_group[i](this->mapping[k], this->L);
			U(idx, k) += this->symmetry_eigVal[i] / (this->normalisation[k] * sqrt(this->symmetry_group.size()));
		}
	}
	return U;
}
arma::cx_vec IsingModel_sym::symmetryRotation(const arma::cx_vec& state) const {
	arma::cx_vec output(ULLPOW(this->L), arma::fill::zeros);
#pragma omp parallel for
	for (long int k = 0; k < this->N; k++) {
		for (int i = 0; i < this->symmetry_group.size(); i++) {
			auto idx = this->symmetry_group[i](this->mapping[k], this->L);
			output(idx) += this->symmetry_eigVal[i] / (this->normalisation[k] * sqrt(this->symmetry_group.size())) * state(k);
		}
	}
	return output;
}
arma::sp_cx_mat IsingModel_sym::fourierTransform(op_type op, int q) const {
	auto beta = std::make_unique<IsingModel_disorder>(this->L, this->J, 0, this->g, 0, this->h, 0, this->_BC);
	auto fullMatrix = beta->fourierTransform(op, q);
	auto U = this->symmetryRotation();
	return U.t() * fullMatrix * U;
}