#include "include/IsingModel.h"
template<typename T> IsingModel<T>::~IsingModel() {}

// ------------------------------------------------------------------------------------------------ PRINTERS ------------------------------------------------------------------------------------------------

/// <summary>
/// prints the basis vector in the given spin-symmetry block
/// </summary>
/// <param name="Sz"> custom spin symmetry block </param>
template <typename T> void IsingModel<T>::print_base_spin_sector(int Sz) {
	std::vector<bool> temp(L);
	int summm = 0;
	int p = +1;
	int z = +1;
	Col<int> check_sectors(std::pow(2, L), arma::fill::zeros);
	for (int k = 0; k < N; k++) {
		int_to_binary(map(k), temp);
		if (std::accumulate(temp.begin(), temp.end(), 0) == Sz)
			stout << k << "\t\t" << temp << endl;
	}
}

/// <summary>
/// prints the state in the regular basis (only the highest coefficients are present)
/// </summary>
/// <param name="_id"> index of the printed state </param>
template <typename T> void IsingModel<T>::print_state(u64 _id) {
	vec state = abs(eigenvectors.col(_id));
	double max = arma::max(state);
	std::vector<bool> base_vector(L);
	for (int k = 0; k < N; k++) {
		int_to_binary(map(k), base_vector);
		if (abs(state(k)) >= 0.01 * max) {
			stout << state(k) << " * |" << base_vector << "> + ";
		}
	}
	stout << endl;
}

/// <summary>
/// Sets the neigbours depending on the Boundary condition (BC) defined as a makro in the 'headers.h' file
/// </summary>
template <typename T> void IsingModel<T>::set_neighbors() {
	this->nearest_neighbors = std::vector<int>(L, 0);
	this->next_nearest_neighbors = std::vector<int>(L, 0);
	switch (this->_BC) {
	case 0:
		// periodic boundary conditions
		for (int i = 0; i < this->L; i++) {
			this->nearest_neighbors[i] = (i + 1) % this->L;
			this->next_nearest_neighbors[i] = (i + 2) % this->L;
		}
		break;
	case 1:
		// open boundary conditions
		for (int i = 0; i < this->L; i++) {
			this->nearest_neighbors[i] = (i + 1) % this->L;
			this->next_nearest_neighbors[i] = (i + 2) % this->L;
		}
		NO_OVERFLOW(
			this->nearest_neighbors[L - 1] = -1;
			this->next_nearest_neighbors[L - 2] = -1;
			this->next_nearest_neighbors[L - 1] = -1;
		);
		break;
	default:
		for (int i = 0; i < this->L; i++) {
			this->nearest_neighbors[i] = (i + 1) % this->L;
			this->next_nearest_neighbors[i] = (i + 2) % this->L;
		}
		break;
	}
}

// ------------------------------------------------------------------------------------------------ DIAGONALIZATIONS ------------------------------------------------------------------------------------------------

/// <summary>
/// General procedure to diagonalize the Hamiltonian using eig_sym from the Armadillo library
/// </summary>
template <typename T> void IsingModel<T>::diagonalization(bool withoutEigenVec) {
	//out << real(H) << endl;
	try {
		if (withoutEigenVec) arma::eig_sym(this->eigenvalues, arma::Mat<T>(this->H));
		else				 arma::eig_sym(this->eigenvalues, this->eigenvectors, arma::Mat<T>(this->H));
	}
	catch (const bad_alloc& e) {
		stout << "Memory exceeded" << e.what() << "\n";
		stout << "dim(H) = " << H.size() * sizeof(H(0, 0)) << "\n";
		assert(false);
	}
	//for (long int i = 0; i < N; i++)
	//	this->eigenvectors.col(i) = arma::normalise(this->eigenvectors.col(i));

	double E_av = trace(eigenvalues) / double(N);
	auto i = min_element(begin(eigenvalues), end(eigenvalues), [=](int x, int y) {
		return abs(x - E_av) < abs(y - E_av);
		});
	this->E_av_idx = i - begin(eigenvalues);
}

/// <summary>
/// calculates the total spin from the correlation matrix
/// </summary>
/// <param name="corr_mat"> spin correlation matrix </param>
/// <returns></returns>
template <typename T> double IsingModel<T>::total_spin(const mat& corr_mat) {
	double S2 = arma::accu(corr_mat);
	if (S2 < -0.25) return 0;
	return (sqrt(1 + 4 * S2) - 1.0) / 2.0;
}
// ------------------------------------------------------------------------------------------------ THERMODYNAMIC QUANTITIES ------------------------------------------------------------------------------------------------

/// <summary>
/// 
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="temperature"></param>
/// <returns></returns>
template <typename Type> std::tuple<arma::vec, arma::vec, arma::vec> IsingModel<Type>::thermal_quantities(const arma::vec& temperature) {
	arma::vec Cv(temperature.size(), arma::fill::zeros);
	arma::vec S(temperature.size(), arma::fill::zeros);
	arma::vec E(temperature.size(), arma::fill::zeros);
#pragma omp parallel for shared(temperature)
	for (int k = 0; k < temperature.size(); k++) {
		const double T = temperature(k);
		double Z = 0;
		double E_av = 0, E_av2 = 0;
		for (long int i = 0; i < this->N; i++) {
			double gibbs = std::exp(-(this->eigenvalues(i) - this->eigenvalues(0)) / T);
			Z += gibbs;
			E_av += this->eigenvalues(i) * gibbs;
			E_av2 += this->eigenvalues(i) * this->eigenvalues(i) * gibbs;
		}
		E_av /= Z;
		E_av2 /= Z;
		E(k) = E_av;
		S(k) = (std::log(Z) + (E_av - this->eigenvalues(0)) / T) / (double)L;
		Cv(k) = (E_av2 - E_av * E_av) / double(L * T * T);
	}
	return std::make_tuple(Cv, S, E);
}

// ------------------------------------------------------------------------------------------------ ERGODIC QUANTITIES ------------------------------------------------------------------------------------------------

/// <summary>
/// The IPR, also called the participation ratio, quantifies how delocalized a state is in a certain basis.
/// The state is completely delocalized, when:
/// IPR=dim(hilbert space)
/// </summary>
/// <param name="state_idx"> index of the eigenvector used to calculate this quantity </param>
/// <returns> returns the IPR value </returns>
template <typename T> double IsingModel<T>::ipr(int state_idx) const {
	double ipr = 0;
	arma::subview_col state = eigenvectors.col(state_idx);
#pragma omp parallel for reduction(+: ipr)
	for (int n = 0; n < N; n++) {
		double value = abs(conj(state(n)) * state(n));
		ipr += value * value;
	}
	return 1.0 / ipr;
}

/// <summary>
/// The information entropy is basicaly calculated similarly to the inverse participation ratio but changes the way we add the values multiplying each
/// by the logarithm in information entropy way
/// </summary>
/// <param name="_id">index of the eigenstate</param>
/// <returns>The information entropy</returns>
template <typename T> double IsingModel<T>::information_entropy(u64 _id) const {
	arma::subview_col state = this->eigenvectors.col(_id);
	double ent = 0;
#pragma omp parallel for reduction(+: ent)
	for (int k = 0; k < this->N; k++) {
		double val = abs(conj(state(k)) * state(k));
		ent += val * log(val);
	}
	return -ent / log(0.48 * this->N);
}

/// <summary>
/// Information entropy caluclated in the basis of the other model, quantity needed to check how te perturbation
/// effects the information entropy at all and what is the difference between them
/// </summary>
/// <param name="_id">Index of the state in alfa sector</param>
/// <param name="beta">second model</param>
/// <param name="min">minimum state for beta basis</param>
/// <param name="max">maximum state for beta basis</param>
/// <returns>information entropy in beta model basis</returns>
template <typename T> double IsingModel<T>::information_entropy(u64 _id, const IsingModel<T>& beta, u64 _min, u64 _max) const {
	arma::subview_col state_alfa = this->eigenvectors.col(_id);
	double ent = 0;
#pragma omp parallel for reduction(+: ent)
	for (long k = (long)_min; k < (long)_max; k++) {
		cpx c_k = cdot(beta.get_eigenState(k), state_alfa);
		double val = abs(conj(c_k) * c_k);
		ent += val * log(val);
	}
	return -ent / log(0.48 * this->N);
}

/// <summary>
/// Calculates the energy-level statistics within the energy window denoted by the indices _min and _max
/// computed as in: PHYSICAL REVIEW B 91, 081103(R) (2015)
/// </summary>
/// <param name="_min"> index of eigenenergy, being the lower bound of energy window </param>
/// <param name="_max"> index of eigenenergy, being the upper bound of energy window </param>
/// <returns></returns>
template <typename T> double IsingModel<T>::eigenlevel_statistics(u64 _min, u64 _max) const {
	double r = 0;
	if (_min <= 0) assert(false && "too low index");
	if (_max >= N) assert(false && "index exceeding Hilbert space");
#pragma omp parallel for reduction(+: r)
	for (long k = (long)_min; k < (long)_max; k++) {
		NO_OVERFLOW(
			const double delta_n = eigenvalues(k) - eigenvalues(k - 1);
			const double delta_n_next = eigenvalues(k + 1) - eigenvalues(k);
		);
		const double min = std::min(delta_n, delta_n_next);
		const double max = std::max(delta_n, delta_n_next);
		if (abs(delta_n) <= 1e-15) assert(false && "Degeneracy!!!\n");
		r += min / max;
	}
	return r / double(_max - _min);
}


/// <summary>
/// Calculates the energy-level statistics within the energy window denoted by the indices _min and _max
/// computed as in: PHYSICAL REVIEW B 91, 081103(R) (2015)
/// </summary>
/// <returns>Vector for whole spectrum eigenlevel statistics</returns>
template <typename T> vec IsingModel<T>::eigenlevel_statistics_with_return() const {
	vec r(N - 2);
#pragma omp parallel for shared(r)
	for (int k = 1; k < N - 1; k++) {
		NO_OVERFLOW(r(k - 1) = eigenlevel_statistics(k, k + 1);)
	}
	return r;
}

/// <summary>
/// Calculates the average spectrum repulsion in the system as the average of
/// the x-component spin matrix in the eigenstates
/// </summary>
/// <typeparam name="T"> typename as class: model with disorder or symmetries </typeparam>
/// <param name="op"> operator as function acting on a specific eigenstate </param>
/// <param name="A"> class instance </param>
/// <param name="site"> position of the spin, where the operator is acted upon </param>
/// <returns> returns the average spectrum repulsion </returns>
template <typename T>
double IsingModel<T>::spectrum_repulsion(double (IsingModel::* op)(int, int), IsingModel& A, int site) {
	//double rn_next = 0, rn = (A.*op)(0, 1);
	double average = 0;
#pragma omp parallel for reduction(+:average)
	for (int k = 1; k < A.N; k++) {
		const double rn = (A.*op)(k - 1, 1);
		const double rn_next = (A.*op)(k, 1);
		average += abs(rn_next - rn);
		//rn = rn_next;
	}
	return average / (A.N - 1.0);
}

/// <summary>
/// averages the distance between consecutive energy levels in the system
/// </summary>
/// <typeparam name="T"> does not matter </typeparam>
/// <returns> mean level spacing </returns>
template <typename T> double IsingModel<T>::mean_level_spacing_av(u64 _min, u64 _max) const {
	if (_min <= 0) throw "too low index";
	if (_max > N) throw "index exceeding Hilbert space";
	double omega_H = 0;
#pragma omp parallel for reduction(+: omega_H)
	for (long int k = (long)_min; k < (long)_max; k++) {
		NO_OVERFLOW(omega_H += this->eigenvalues(k) - this->eigenvalues(k - 1););
	}
	return omega_H / double(_max - _min);
}

/// <summary>
/// 
/// </summary>
/// <typeparam name="T"></typeparam>
/// <returns></returns>
template <typename T> double IsingModel<T>::mean_level_spacing_trace() const {
	const double chi = 0.341345;
	double trace_H2 = 0;
	double trace_H = 0;
#pragma omp parallel for reduction(+: trace_H, trace_H2)
	for (int k = 0; k < N; k++) {
		trace_H += this->eigenvalues(k);
		trace_H2 += this->eigenvalues(k) * this->eigenvalues(k);
	}
	return sqrt(trace_H2 / double(N) - trace_H * trace_H / double(N * N)) / (chi * N);
}

// ----------------------------------------------------------- OPERATORS AND AVERAGES -------------------------------------------------------

/// <summary>
/// Prints to file the average of a given operator in each eigenstate as a function of eigenenergies
/// </summary>
/// <param name="op"> operator as function acting on a specific eigenstate </param>
/// <param name="A"> class instance </param>
/// <param name="site"> position of the spin, where the operator is acted upon </param>
/// <param name="name"> name of the file to print data </param>
/// <param name="separator"> separator between columns in file </param>
template <typename T>
void IsingModel<T>::operator_av_in_eigenstates(double (IsingModel<T>::* op)(int, int), IsingModel<T>& A, int site, std::string name, std::string separator) {
	std::ofstream file(name);                                                                  // file to write the average to
	if (!file.is_open()) throw "Can't open file " + name + "\n Choose another file\n";
	// vec res(A.get_hilbert_size(), fill::zeros);

#pragma omp parallel for
	for (int k = 0; k < A.get_hilbert_size(); k++) {
		const double res = (A.*op)(k, site);
#pragma omp critical
		file << A.eigenvalues(k) / (double)A.L << separator << res << endl;
	}
	//for (int k = 0; k < A.get_hilbert_size(); k++)
	//    file << A.eigenvalues(k) / (double)A.L << separator << res(k) << endl;
	file.close();
}

/// <summary>
///
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="op"></param>
/// <param name="A"></param>
/// <param name="site"></param>
/// <returns></returns>
template <typename T>
vec IsingModel<T>::operator_av_in_eigenstates_return(double (IsingModel<T>::* op)(int, int), IsingModel<T>& A, int site) {
	vec temp(A.get_hilbert_size(), fill::zeros);
#pragma omp parallel for shared(temp)
	for (int k = 0; k < A.get_hilbert_size(); k++)
		temp(k) = (A.*op)(k, site);
	return temp;
}

//-------------------------------------------------------------------------------------------------LIOMs

template <typename _type>
sp_cx_mat IsingModel<_type>::create_StringOperator(coordinate alfa, coordinate beta, int j, int ell) const {
	if (ell < 0) assert(false && "last argument is positive, duh");
	op_type op_alfa, op_beta;
	switch (alfa) {
	case coordinate::x: op_alfa = IsingModel::sigma_x; break;
	case coordinate::y: op_alfa = IsingModel::sigma_y; break;
	case coordinate::z: op_alfa = IsingModel::sigma_z; break;
	}
	switch (beta) {
	case coordinate::x: op_beta = IsingModel::sigma_x; break;
	case coordinate::y: op_beta = IsingModel::sigma_y; break;
	case coordinate::z: op_beta = IsingModel::sigma_z; break;
	}
	sp_cx_mat SigmaAlfa = create_operator({ op_alfa }, { properSite(j) });
	sp_cx_mat SigmaBeta = create_operator({ op_beta }, { properSite(j + ell) });
	sp_cx_mat SigmaZstring;
	if (ell == 0) {
		if (alfa == beta && alfa == coordinate::y) return -create_operator({ IsingModel_sym::sigma_z });
		else assert(false && "Don't know what this could possibly be, check again if you need this");
	}
	else if (ell == 1) {
		SigmaZstring = arma::eye<sp_cx_mat>(this->N, this->N);
	}
	else {
		std::vector<int> sites;
		for (int k = 1; k <= ell - 1; k++)
			sites.push_back(properSite(j + k));
		SigmaZstring = create_operator({ IsingModel_sym::sigma_z }, sites);
	}
	return SigmaAlfa * SigmaZstring * SigmaBeta;
}

template <typename _type>
sp_cx_mat IsingModel<_type>::create_LIOMoperator_densities(int n, int j) const {
	if (n < 0) assert(false && "Only positive integers for LIOMs");
	auto S = [this](coordinate alfa, coordinate beta, int j, int ell) {
		return create_StringOperator(alfa, beta, j, ell);
	};
	if (n % 2 == 0) {
		//if (n == 0)
		//	return SpMat<cpx>(arma::real(this->H), arma::imag(this->H));
		//else
			return this->J * (S(coordinate::x, coordinate::x, j, j + n) + S(coordinate::y, coordinate::y, j, j + n - 2))
			+ this->g * (S(coordinate::x, coordinate::x, j, j + n - 1) + S(coordinate::y, coordinate::y, j, j + n - 1));
	}
	else {
		return S(coordinate::x, coordinate::y, j, j + n) - S(coordinate::y, coordinate::x, j, j + n);
	}
}

// ------------------------------------------------------------------------------------------------ STATISTICS AND PROBABILITIES ------------------------------------------------------------------------------------------------

/// <summary>
/// Creates a probabilty distribution of data and saves it in the directory
/// </summary>
void probability_distribution(std::string dir, std::string name, const arma::vec& data, int n_bins) {
	std::ofstream file(dir + name + ".dat");
	if (n_bins <= 0)
		n_bins = 1 + long(3.322 * log(data.size()));
	double _min = arma::min(data);
	double _max = arma::max(data);
	NO_OVERFLOW(arma::vec prob_dist(n_bins + 1, arma::fill::zeros);)
	prob_dist = normalise_dist(arma::conv_to<arma::vec>::from(arma::hist(data, n_bins)), _min, _max);
	const double std_dev = arma::stddev(data);
	const double mean = 0.0;// arma::mean(data);
	const double step = abs(_max - _min) / (double)n_bins;
	for (int p = 0; p < n_bins; p++) {
		double x = p * step + _min;
		file << x << "\t\t" << prob_dist(p) << "\n";// << "\t\t" << gaussian(x, mean, std_dev) << std::endl;
	}
	file.close();
}

/// <summary>
/// Creates a probabilty distribution of data and returns vector from it
/// </summary>
arma::vec probability_distribution_with_return(const arma::vec& data, int n_bins) {
	if (n_bins <= 0)
		n_bins = 1 + long(3.322 * log(data.size()));
	return normalise_dist(arma::conv_to<arma::vec>::from(arma::hist(data, n_bins)),
		arma::min(data), arma::max(data));
}

/// <summary>
/// Takes the vector and averages over small buckets of size mu symmetricaly in order to maintain only fluctuations
/// </summary>
/// <param name="data">vector to be characterized</param>
/// <param name="mu">size of the bucket</param>
/// <returns></returns>
arma::vec data_fluctuations(const arma::vec& data, int mu) {
	arma::vec fluct(data.size() - mu, arma::fill::zeros);
	assert(mu < data.size() && "Bucket exceeds data container\nTry again\n");
	int end = (int)data.size() - mu / 2;
#pragma omp parallel for shared(fluct, end, mu, data)
	for (int k = mu / 2; k < end; k++) {
		double average = 0;
		for (int n = k - mu / 2; n < k + mu / 2; n++)
			average += data(n);
		NO_OVERFLOW(fluct(k - mu / 2) = data(k) - average / double(mu);)
	}
	return fluct;
}

/// <summary>
/// Calculates a quantity similar to spectrum repulsion and finds the average of it and the biggest outliers
/// </summary>
/// <param name="data">vector of data for the repulsion to be performed on</param>
/// <returns></returns>
arma::vec statistics_average(const arma::vec& data, int num_of_outliers) {
	NO_OVERFLOW(
		std::vector<double> spec_rep(num_of_outliers + 1, INT_MIN);
		double average = 0;
	for (int k = 1; k < data.size() - 1; k++) {
		double repulsion = abs(data(k) - data(k - 1));
		average += repulsion;
		int i = num_of_outliers + 1;
		while (i > 1) {
			if (repulsion < spec_rep[i - 1]) break;
			i--;
		}
		if (i > 0 && i <= num_of_outliers) {
			spec_rep.insert(spec_rep.begin() + i, repulsion);
			spec_rep.pop_back();
		}
	}
	);
	spec_rep[0] = average / double(data.size() - 2);
	return (vec)spec_rep;
}

// ------------------------------------------------------------------------------------------------ TOOLS ------------------------------------------------------------------------------------------------

/// <summary>
/// Overlapping of two eigenstates of possibly different matrices A and B
/// </summary>
/// <param name="A">matrix A</param>
/// <param name="B">matrix B</param>
/// <param name="n_a">number of A eigenstate</param>
/// <param name="n_b">number of B eigenstate</param>
/// <returns>A_n_a dot n_b_B</returns>
template <typename T>
T overlap(const IsingModel<T>& A, const IsingModel<T>& B, int n_a, int n_b) {
	if (A.get_hilbert_size() != B.get_hilbert_size()) throw "Incompatible Hilbert dimensions\n";
	if (n_a >= A.get_hilbert_size() || n_b >= B.get_hilbert_size() || n_a < 0 || n_b < 0) throw "Cannot create an overlap between non-existing states\n";
	return arma::cdot(A.get_eigenState(n_a), B.get_eigenState(n_b));
}

/// <summary>
/// Calculates the overlap between two given states
/// </summary>
/// <param name="A">Model for the left eigenvector</param>
/// <param name="B">Model for the right eigenvector</param>
/// <param name="n_a">number of the left eigenvector</param>
/// <param name="n_b">number of the right eigenvector</param>
/// <returns></returns>
template<> cpx overlap<cpx>(const IsingModel<cpx>& A, const IsingModel<cpx>& B, int n_a, int n_b) {
	if (A.get_hilbert_size() != B.get_hilbert_size()) throw "Incompatible Hilbert dimensions\n";
	if (n_a >= A.get_hilbert_size() || n_b >= B.get_hilbert_size() || n_a < 0 || n_b < 0) throw "Cannot create an overlap between non-existing states\n";
	double overlap_real = 0, overlap_imag = 0;
	auto state_A = A.get_eigenState(n_a);
	auto state_B = B.get_eigenState(n_b);
	state_A = arma::normalise(state_A);
	state_B = arma::normalise(state_B);
#pragma omp parallel for reduction(+: overlap_real, overlap_imag)
	for (int k = 0; k < A.get_hilbert_size(); k++) {
		cpx over = conj(state_A(k)) * state_B(k);
		overlap_real += real(over);
		overlap_imag += imag(over);
	}
	return cpx(overlap_real, overlap_imag);
}


template <typename _type>
void IsingModel<_type>::time_evolve_state(arma::cx_vec& state, double time) {
	//assert(state.size() == this->N, ("state is of wrong dimensions! set dimension to:" + std::to_string(this->N)));
	if (this->eigenvectors.is_empty()) {
		stout << "Model not diagonalized! Did it for you, but dude...";
		this->diagonalization();
	}
	for (long k = 0; k < this->N; k++)
		state += std::exp(-im * eigenvalues(k) * time) * arma::dot(eigenvectors.col(k), state) * eigenvectors.col(k);
	state = arma::normalise(state);
}

// UNRESOLVED TEMPLATE EXTERNALS <- COMPILER DOESN"T KNOW ABOUT THEM SADLY
template IsingModel<double>::~IsingModel();
template IsingModel<cpx>::~IsingModel();
template void IsingModel<double>::set_neighbors();
template void IsingModel<cpx>::set_neighbors();
template void IsingModel<cpx>::diagonalization(bool);
template void IsingModel<double>::diagonalization(bool);
template double IsingModel<cpx>::eigenlevel_statistics(u64, u64) const;
template double IsingModel<double>::eigenlevel_statistics(u64, u64) const;
template vec IsingModel<cpx>::eigenlevel_statistics_with_return() const;
template vec IsingModel<double>::eigenlevel_statistics_with_return() const;
template double IsingModel<double>::ipr(int) const;
template double IsingModel<cpx>::ipr(int) const;
template double IsingModel<double>::information_entropy(u64) const;
template double IsingModel<cpx>::information_entropy(u64) const;
template double IsingModel<double>::information_entropy(u64, const IsingModel<double>&, u64, u64) const;
template double IsingModel<cpx>::information_entropy(u64, const IsingModel<cpx>&, u64, u64) const;
template cpx overlap(const IsingModel<cpx>&, const IsingModel<cpx>&, int, int);
template double overlap(const IsingModel<double>&, const IsingModel<double>&, int, int);
template double IsingModel<cpx>::total_spin(const mat&);
template double IsingModel<double>::total_spin(const mat&);
template std::tuple<arma::vec, arma::vec, arma::vec> IsingModel<double>::thermal_quantities(const arma::vec&);
template std::tuple<arma::vec, arma::vec, arma::vec> IsingModel<cpx>::thermal_quantities(const arma::vec&);
template double IsingModel<double>::mean_level_spacing_av(u64, u64) const;
template double IsingModel<cpx>::mean_level_spacing_av(u64, u64) const;
template double IsingModel<double>::mean_level_spacing_trace() const;
template double IsingModel<cpx>::mean_level_spacing_trace() const;
template sp_cx_mat IsingModel<cpx>::create_StringOperator(coordinate, coordinate, int,int) const;
template sp_cx_mat IsingModel<double>::create_StringOperator(coordinate, coordinate, int, int) const;
template sp_cx_mat IsingModel<cpx>::create_LIOMoperator_densities(int, int) const;
template sp_cx_mat IsingModel<double>::create_LIOMoperator_densities(int, int) const;
template void IsingModel<double>::time_evolve_state(arma::cx_vec&, double);
template void IsingModel<cpx>::time_evolve_state(arma::cx_vec&, double);