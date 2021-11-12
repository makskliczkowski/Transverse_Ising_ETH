#include "include/IsingModel.h"

// ------------------------------------------------------------------------------------------------ CONSTRUCTORS ------------------------------------------------------------------------------------------------
IsingModel_sym::IsingModel_sym(int L, double J, double g, double h, int k_sym, bool p_sym, bool x_sym, int _BC) {
	this->L = L; this->J = J; this->g = g; this->h = h; this->_BC = _BC;
	symmetries.k_sym = k_sym * two_pi / double(this->L);
	symmetries.p_sym = (p_sym) ? 1.0 : -1.0;
	symmetries.x_sym = (x_sym) ? 1.0 : -1.0;
	k_sector = abs(this->symmetries.k_sym) < 1e-4 || abs(this->symmetries.k_sym - pi) < 1e-4;
	// precalculate the exponents
	this->k_exponents = v_1d<cpx>(this->L, 0.0);
#pragma omp parallel for
	for (int l = 0; l < this->L; l++)
		this->k_exponents[l] = std::exp(-1i * this->symmetries.k_sym * double(l));
	this->createSymmetryGroup();

	this->info = "_L=" + std::to_string(this->L) + \
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
				this->symmetry_eigVal.push_back(this->k_exponents[k] * double(this->symmetries.p_sym * this->symmetries.x_sym));
			}
		}
		T = multiply_operators(rotate_left, T);
	}
}

/// <summary>
/// Takes the state from the index in mapping
/// </summary>
u64 IsingModel_sym::map(u64 index) {
	if (index >= this->N) throw "Element out of range\n No such index in map\n";
	return this->mapping[index];
}

/// <summary>
/// Finds the representative of the equivalent class after vector rotation (mainly used after applied another symmetry)
/// </summary>
/// <param name="base_vector"> vector from EC to find representative </param>
/// <returns> index of the representative state in the EC </returns>
std::tuple<u64, int> IsingModel_sym::find_translation_representative(u64 base_idx) const {
	u64 EC_sym = base_idx;
	u64 idx = INT_MAX;
	u64 base_rotate = EC_sym;
	int counter = 1, count_to_rep = 1;
	while (idx != base_idx) {
		idx = rotate_left(base_rotate, this->L); // new rotate by decimal representation
		base_rotate = idx;
		counter++;
		if (idx < EC_sym) {
			EC_sym = idx;
			count_to_rep = counter;
		}
	}
	return std::make_tuple(EC_sym, count_to_rep);
}

/// <summary>
/// Find representatives of other EC generated by reflection, spin-flip and (reflection x spin-flip) symmetry
/// </summary>
/// <param name="base_vector"> current base vector to act with symmetries </param>
/// <param name="min"> index of EC class representative by translation symmetry </param>
/// <returns></returns>
std::tuple<u64, int> IsingModel_sym::find_SEC_representative(u64 base_idx, std::vector<u64>& minima) const {
	u64 R = INT_MAX;
	int eig_Rsym = 0;
	u64 idx_R = 0;
	if (this->k_sector) {
		idx_R = reverseBits(base_idx, this->L);
		auto tupleR = find_translation_representative(idx_R);
		R = std::get<0>(tupleR);
		eig_Rsym = std::get<1>(tupleR);
	}
	minima[0] = R;

	if (this->h == 0) {
		// check spin-flip
		auto tupleX = find_translation_representative(flip(base_idx, this->L));
		minima[1] = std::get<0>(tupleX);

		// check spin-flip and reflection
		int eig_RXsym = INT_MAX;
		if (this->k_sector) {
			auto tupleRX = find_translation_representative(flip(idx_R, this->L));
			minima[2] = std::get<0>(tupleRX);
			eig_RXsym = std::get<1>(tupleRX);
		}

		switch (std::min_element(minima.begin(), minima.end()) - minima.begin()) {
		case 1:
			return std::make_tuple(minima[1], symmetries.x_sym * std::get<1>(tupleX)); break;
		case 0:
			return std::make_tuple(minima[0], symmetries.p_sym * eig_Rsym); break;
		case 2:
			return std::make_tuple(minima[2], symmetries.p_sym * symmetries.x_sym * eig_RXsym); break;
		default: throw "Index out of range\n";
		}
	}
	else
		return std::make_tuple(R, symmetries.p_sym * eig_Rsym);
}

/// <summary>
/// From applying symmetry operators the function finds the normalisation for a given state
/// </summary>
cpx IsingModel_sym::get_symmetry_normalization(u64 base_idx) {
	cpx normalisation = cpx(0.0, 0.0);

	u64 Translation = base_idx;
	u64 PT = (k_sector) ? reverseBits(base_idx, this->L) : INT_MAX;
	u64 ZT = flip(base_idx, this->L);
	u64 PZT = (k_sector) ? flip(PT, this->L) : INT_MAX;

	for (int l = 0; l < this->L; l++) {
		if (Translation == base_idx)
			normalisation += this->k_exponents[l];
		Translation = rotate_left(Translation, this->L);
		if (this->k_sector) {
			if (PT == base_idx)
				normalisation += ((double)this->symmetries.p_sym) * this->k_exponents[l];
			PT = rotate_left(PT, this->L);
		}
		if ((this->h == 0)) {
			if (ZT == base_idx)
				normalisation += ((double)this->symmetries.x_sym) * this->k_exponents[l];
			ZT = rotate_left(ZT, this->L);
			if (this->k_sector) {
				if (PZT == base_idx)
					normalisation += ((double)this->symmetries.p_sym * (double)this->symmetries.x_sym) * this->k_exponents[l];
				PZT = rotate_left(PZT, this->L);
			}
		}
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
	std::vector<u64> minima(3);
	for (u64 j = start; j < stop; j++) {
		if (this->g == 0) {
			std::vector<bool> base_vector(this->L); // temporary dirac-notation base vector
			int_to_binary(j, base_vector);
			if (std::accumulate(base_vector.begin(), base_vector.end(), 0) != this->L / 2.) continue;
		}
		//check translation
		//auto [min, trans_eig] = find_translation_representative(j);
		//
		//u64 min_R_RX = INT_MAX;
		//if (min == j) {
		//	minima.assign(3, INT_MAX);
		//	auto tuple = find_SEC_representative(j, minima);
		//	min_R_RX = std::get<0>(tuple);
		//}
		//if (min_R_RX < j) continue;
		//std::min(min, min_R_RX)
		u64 SEC = INT64_MAX;
		for (auto& sym_op : this->symmetry_group) {
			u64 new_idx = sym_op(j, this->L);
			if (new_idx < SEC) SEC = new_idx;
		}

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
			stop = ((t + 1) == num_of_threads ? two_powL : u64(two_powL / (double)num_of_threads * (double)(t + 1)));
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
	int sym_eig = 1;
	if (idx > sector_alfa.N) {
		std::vector<u64> minima(3);
		auto tup_T = sector_alfa.find_translation_representative(base_idx);
		auto tup_S = sector_alfa.find_SEC_representative(base_idx, minima);
		auto [min, trans_eig] = (std::get<0>(tup_T) > std::get<0>(tup_S)) ? tup_S : tup_T;
		sym_eig = trans_eig;
		//finding index in reduced Hilbert space
		idx = binary_search(sector_alfa.mapping, 0, sector_alfa.N - 1, min);
	}
	if (idx < sector_alfa.N) {
		cpx translation_eig = conj(sector_alfa.k_exponents[abs(sym_eig) - 1]);
		cpx val = translation_eig * (sector_alfa.normalisation[idx] / normalisation_beta) * double(sgn(sym_eig));
		return std::make_pair(idx, val);
	}
	else
		return std::make_pair(INT_MAX, 0);
}

/// <summary>
/// Sets the non-diagonal elements of the Hamimltonian matrix with symmetry sectors: therefore the matrix elements are summed over the SEC
/// </summary>
/// <param name="k"> index of the basis state acted upon with the Hamiltonian </param>
/// <param name="value"> value of the given matrix element to be set </param>
/// <param name="temp"> resulting vector form acting with the Hamiltonian operator on the k-th basis state </param>
void IsingModel_sym::setHamiltonianElem(u64 k, double value, u64 new_idx) {
	auto [idx, sym_eig] = find_rep_and_sym_eigval(new_idx, *this, this->normalisation[k]);
	if (idx < this->N)
		H(idx, k) += value * sym_eig;
}

/// <summary>
/// Generates the total Hamiltonian of the system. The diagonal part is straightforward,
/// while the non-diagonal terms need the specialized setHamiltonainElem(...) function
/// </summary>
void IsingModel_sym::hamiltonian() {
	try {
		this->H = cx_mat(this->N, this->N, fill::zeros); //hamiltonian
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
				u64 new_idx = flip(this->mapping[k], BinaryPowers[L - 1 - j], (L - 1 - j));
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

// ------------------------------------------------------------------------------------------------ PHYSICAL QUANTITTIES ------------------------------------------------------------------------------------------------

/// <summary>
/// <summary>
/// Calculates the matrix element for sigma_z Pauli matrix
/// </summary>
/// <param name="sites">Sites the matrix works on</param>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <returns>The matrix element</returns>
double IsingModel_sym::av_sigma_z(u64 alfa, u64 beta, std::initializer_list<int> sites) {
	auto sig_z = IsingModel_sym::sigma_z;
	return real(av_operator(alfa, beta, *this, *this, sig_z, sites));
}

/// <summary>
/// Calculates the matrix element for sigma_z extensive
/// </summary>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <returns>The matrix element</returns>
double IsingModel_sym::av_sigma_z(u64 alfa, u64 beta) {
	auto sig_z = IsingModel_sym::sigma_z;
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
	auto sig_z = IsingModel_sym::sigma_z;
	return real(av_operator(alfa, beta, *this, *this, sig_z, corr_length));
}

/// <summary>
/// Calculates the matrix element for sigma_x Pauli matrix
/// </summary>
/// <param name="sites">Sites the matrix works on</param>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <returns>The matrix element</returns>
double IsingModel_sym::av_sigma_x(u64 alfa, u64 beta, std::initializer_list<int> sites) {
	auto sig_x = IsingModel_sym::sigma_x;
	return real(av_operator(alfa, beta, *this, *this, sig_x, sites));
}

/// <summary>
/// Calculates the matrix element for sigma_x extensive
/// </summary>
/// <param name="alfa">Left state</param>
/// <param name="beta">Right state</param>
/// <returns>The matrix element</returns>
double IsingModel_sym::av_sigma_x(const u64 alfa, const u64 beta) {
	auto sig_x = IsingModel_sym::sigma_x;
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
	auto sig_x = IsingModel_sym::sigma_x;
	return real(av_operator(alfa, beta, *this, *this, sig_x, corr_length));
}

/// <summary>
///
/// </summary>
/// <param name="alfa"></param>
/// <param name="beta"></param>
/// <returns></returns>
double IsingModel_sym::av_spin_flip(u64 alfa, u64 beta) {
	auto spin_flip = IsingModel_sym::spin_flip;
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
double IsingModel_sym::av_spin_flip(u64 alfa, u64 beta, std::initializer_list<int> sites) {
	auto spin_flip = IsingModel_sym::spin_flip;
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
	auto spin_flip = IsingModel_sym::spin_flip;
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
cpx IsingModel_sym::av_spin_current(u64 alfa, u64 beta, std::initializer_list<int> sites) {
	auto spin_flip = IsingModel_sym::spin_flip;
	auto value = im * av_operator(alfa, beta, *this, *this, spin_flip, sites);
	value += conj(im * av_operator(beta, alfa, *this, *this, spin_flip, sites));
	return 0.5i * value;
}

/// <summary>
/// calculates spin-spin correlation matrix within a given state
/// </summary>
/// <param name="state_idx"> index of eigenstate to calculate correlation matrix </param>
/// <returns> correlation matrix </returns>
mat IsingModel_sym::correlation_matrix(u64 state_idx) {
	mat corr_mat(this->L, this->L, arma::fill::zeros);
	auto spin_flip = IsingModel_sym::spin_flip;
	auto sig_z_nn = IsingModel_sym::sigma_z;
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			auto kinetic = 0.5 * av_operator(state_idx, state_idx, *this, *this, spin_flip, { i, j });
			kinetic += conj(kinetic);
			auto corr_zz = av_operator(state_idx, state_idx, *this, *this, sig_z_nn, { i, j });
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
cpx av_operator(u64 alfa, u64 beta, const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta, op_type op, std::initializer_list<int> sites) {
	// throwables
	for (auto& site : sites)
		if ((site < 0 || site >= sec_alfa.L) && sec_alfa.L != sec_beta.L) throw "Site index exceeds chain or incompatible chain lengths. Your choice\n";
	if (sec_alfa.J != sec_beta.J || sec_alfa.h != sec_beta.h || sec_alfa.g != sec_beta.g) throw "incompatible model parameters, sun \\('.')// \n";

	arma::subview_col state_alfa = sec_alfa.eigenvectors.col(alfa);
	arma::subview_col state_beta = sec_beta.eigenvectors.col(beta);

	// calculating normalisation for both sector symmetry groups
	double G = 0;
	double G_alfa = sec_alfa.L;
	double G_beta = G_alfa;
	if (sec_alfa.h == 0) {
		G_alfa += sec_alfa.L;
		G_beta += sec_beta.L;
	}
	if (sec_alfa.k_sector) G_alfa += (sec_alfa.h == 0) ? 2 * sec_alfa.L : sec_alfa.L;
	if (sec_beta.k_sector) G_beta += (sec_beta.h == 0) ? 2 * sec_beta.L : sec_beta.L;
	G = std::sqrt(G_alfa * G_beta);

	// going through all sector beta states
	double overlap_real = 0;
	double overlap_imag = 0;
#pragma omp parallel
	{
#pragma omp for reduction(+:overlap_real, overlap_imag)
		for (int i = 0; i < sec_beta.N; i++) {
			const u64 base_vec = sec_beta.mapping[i];
			cpx overlap = apply_sym_overlap(state_alfa, state_beta, base_vec, i, sec_alfa, sec_beta, op, sites);
			overlap_real += overlap.real();
			overlap_imag += overlap.imag();
		}
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
	double G_alfa = sec_alfa.L;
	double G_beta = G_alfa;
	if (sec_alfa.h == 0) {
		G_alfa += sec_alfa.L;
		G_beta += sec_beta.L;
	}
	if (sec_alfa.k_sector) G_alfa += (sec_alfa.h == 0) ? 2 * sec_alfa.L : sec_alfa.L;
	if (sec_beta.k_sector) G_beta += (sec_beta.h == 0) ? 2 * sec_beta.L : sec_beta.L;
	G = std::sqrt(G_alfa * G_beta);

	// going through all sector beta states
	double overlap_real = 0;
	double overlap_imag = 0;
#pragma omp parallel
	{
		std::vector<bool> base_vector(sec_beta.L, 0);
#pragma omp for reduction(+:overlap_real, overlap_imag)
		for (int i = 0; i < sec_beta.N; i++) {
			for (int l = 0; l < sec_beta.L; l++) {
				const u64 base_vec = sec_beta.mapping[i];
				cpx overlap = apply_sym_overlap(state_alfa, state_beta, base_vec, i, sec_alfa, sec_beta, op, { l });
				overlap_real += overlap.real();
				overlap_imag += overlap.imag();
			}
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
	double G_alfa = sec_alfa.L;
	double G_beta = G_alfa;
	if (sec_alfa.h == 0) {
		G_alfa += sec_alfa.L;
		G_beta += sec_beta.L;
	}
	if (sec_alfa.k_sector) G_alfa += (sec_alfa.h == 0) ? 2 * sec_alfa.L : sec_alfa.L;
	if (sec_beta.k_sector) G_beta += (sec_beta.h == 0) ? 2 * sec_beta.L : sec_beta.L;
	G = std::sqrt(G_alfa * G_beta);

	auto neis = get_neigh_vector(sec_alfa._BC, sec_alfa.L, corr_len);
	// going through all sector beta states

	double overlap_real = 0;
	double overlap_imag = 0;
#pragma omp parallel
	{
#pragma omp for reduction(+:overlap_real, overlap_imag)
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
	const IsingModel_sym& sec_beta, op_type op, std::initializer_list<int> sites)
{
	const int L = sec_alfa.L;
	auto get_overlap_sym = [&](u64 vec_sym, cpx sym_eig) {
		cpx overlap = 0.0;
		auto [val_sym_beta, vec_sym_tmp] = op(vec_sym, L, sites);
		if (abs(val_sym_beta) > 1e-14) {
			auto [idx_sym, val_sym_alfa] = !(vec_sym == vec_sym_tmp) ? \
				find_rep_and_sym_eigval(vec_sym_tmp, sec_alfa, sec_beta.normalisation[k]) : \
				std::make_pair(vec_sym, conj(sym_eig));
			if (idx_sym < sec_alfa.N)
				overlap = sym_eig * conj(val_sym_alfa * alfa(idx_sym)) * beta(k) * val_sym_beta;
		}
		return overlap;
	};
	cpx overlap = 0.0;
	auto symmetry_group = sec_beta.get_sym_group();
	auto symmetry_eigVal = sec_beta.get_sym_eigVal();
	for (int k = 0; k < symmetry_group.size(); k++) {
		auto sym_operation = symmetry_group[k];
		auto symRepr = symmetry_eigVal[k];
		overlap += get_overlap_sym(sym_operation(base_vec, L), symRepr);
	}
	//u64 Translation = base_vec;					// Translation
	//u64 PT;										// Parity translation
	//u64 ZT;										// Flip translation
	//u64 PZT;									// Parity Flip translation
	//if (sec_beta.k_sector) {
	//	PT = reverseBits(base_vec, L);
	//}
	//if (sec_beta.h == 0) {
	//	ZT = flip(base_vec, L);
	//	if (sec_beta.k_sector) {
	//		PZT = reverseBits(ZT, L);
	//	}
	//}
	//
	//for (int l = 0; l < sec_alfa.L; l++) {
	//	auto T_eig = sec_beta.k_exponents[l];
	//	overlap += get_overlap_sym(Translation, T_eig);
	//	Translation = rotate_left(Translation, L);
	//
	//	if (sec_beta.k_sector) {
	//		auto PT_eig = T_eig * double(sec_beta.symmetries.p_sym);
	//		overlap += get_overlap_sym(PT, PT_eig);
	//		PT = rotate_left(PT, L);
	//	}
	//
	//	if (sec_beta.h == 0) {
	//		auto ZT_eig = T_eig * double(sec_beta.symmetries.x_sym);
	//		overlap += get_overlap_sym(ZT, ZT_eig);
	//		ZT = rotate_left(ZT, L);
	//		if (sec_beta.k_sector) {
	//			auto PZT_eig = ZT_eig * double(sec_beta.symmetries.p_sym);
	//			overlap += get_overlap_sym(PZT, PZT_eig);
	//			PZT = rotate_left(PZT, L);
	//		}
	//	}
	//}
	return overlap;
}


// ----------------------------------------------------------------------------- WRAPPERS FOR SIGMA OPERATORS - creating matrix -----------------------------------------------------------------------------
sp_cx_mat IsingModel_sym::create_operator(op_type op, std::initializer_list<int> sites) {
	// throwables
	for (auto& site : sites)
		if ((site < 0 || site >= this->L)) throw "Site index exceeds chain or incompatible chain lengths. Your choice\n";
	double G = this->L;
	if (this->h == 0)
		G += this->L;
	if (this->k_sector) G += (this->h == 0) ? 2 * this->L : this->L;
	// going through all sector beta states
	sp_cx_mat operator_matrix(this->N, this->N);
	for (int i = 0; i < this->N; i++) {
		const u64 base_vec = this->mapping[i];
		set_OperatorElem(op, sites, operator_matrix, base_vec, i);
	}
	return operator_matrix / (G * sqrt(this->L));
}
sp_cx_mat IsingModel_sym::create_operator(op_type op) {// calculating normalisation for both sector symmetry groups
	double G = this->L;
	if (this->h == 0)
		G += this->L;
	if (this->k_sector) G += (this->h == 0) ? 2 * this->L : this->L;
	// going through all sector beta states
	sp_cx_mat operator_matrix(this->N, this->N);
	for (int i = 0; i < this->N; i++) {
		const u64 base_vec = this->mapping[i];
		for (int j = 0; j < this->L; j++)
			set_OperatorElem(op, { j }, operator_matrix, base_vec, i);
	}
	return operator_matrix / (G * sqrt(this->L));
}
sp_cx_mat IsingModel_sym::create_operator(op_type op, int corr_len) {
	double G = this->L;
	if (this->h == 0)
		G += this->L;
	if (this->k_sector) G += (this->h == 0) ? 2 * this->L : this->L;
	auto neis = get_neigh_vector(this->_BC, this->L, corr_len);
	sp_cx_mat operator_matrix(this->N, this->N);
	for (int i = 0; i < this->N; i++) {
		const u64 base_vec = this->mapping[i];
		for (int j = 0; j < this->L; j++) {
			const int nei = neis[j];
			if (nei >= 0) {
				set_OperatorElem(op, { j, nei }, operator_matrix, base_vec, i);
			}
		}
	}
	return operator_matrix / (G * sqrt(this->L));
}

void IsingModel_sym::set_OperatorElem(op_type op, std::initializer_list<int> sites, sp_cx_mat& operator_matrix, u64 base_vec, u64 cur_idx){
	auto set_MatrixElem = [&](u64 vec_sym, cpx sym_eig) {
		auto [op_value, vec_sym_tmp] = op(vec_sym, L, sites);
		if (abs(op_value) > 1e-14) {
			auto [idx_sym, sym_eigVal] = !(vec_sym == vec_sym_tmp) ? \
				find_rep_and_sym_eigval(vec_sym_tmp, *this, this->normalisation[cur_idx]) : std::make_pair(vec_sym, conj(sym_eig));
			if (idx_sym < N)
				operator_matrix(idx_sym, cur_idx) += sym_eigVal * op_value;
		}
	};
	for (int k = 0; k < this->symmetry_group.size(); k++) {
		auto sym_operation	= this->symmetry_group[k];
		auto symRepr		= this->symmetry_eigVal[k];
		set_MatrixElem(sym_operation(base_vec, this->L), symRepr);
	}

	//u64 Translation = base_vec;					// Translation
	//u64 PT;										// Parity translation
	//u64 ZT;										// Flip translation
	//u64 PZT;									// Parity Flip translation
	//if (this->k_sector)
	//	PT = reverseBits(base_vec, L);
	//if (this->h == 0) {
	//	ZT = flip(base_vec, this->L);
	//	if (this->k_sector)
	//		PZT = reverseBits(ZT, L);
	//}
	//cpx overlap = 0.0;
	//for (int l = 0; l < L; l++) {
	//	auto T_eig = this->k_exponents[l];
	//	set_MatrixElem(Translation, T_eig);
	//	Translation = rotate_left(Translation, this->L);
	//
	//	if (this->k_sector) {
	//		auto PT_eig = T_eig * double(this->symmetries.p_sym);
	//		set_MatrixElem(PT, PT_eig);
	//		PT = rotate_left(PT, this->L);
	//	}
	//	if (this->h == 0) {
	//		auto ZT_eig = T_eig * double(this->symmetries.x_sym);
	//		set_MatrixElem(ZT, ZT_eig);
	//		ZT = rotate_left(ZT, this->L);
	//		if (this->k_sector) {
	//			auto PZT_eig = ZT_eig * double(this->symmetries.p_sym);
	//			set_MatrixElem(PZT, PZT_eig);
	//			PZT = rotate_left(PZT, this->L);
	//		}
	//	}
	//}
}