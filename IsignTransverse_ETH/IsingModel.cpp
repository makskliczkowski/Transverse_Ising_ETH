#include "include/IsingModel.h"
template<typename T> IsingModel<T>::~IsingModel() {}

// ------------------------------------------------------------------------------------------------ INITIALIZE HELPERS ------------------------------------------------------------------------------------------------
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
template <typename T> void IsingModel<T>::diagonalization(bool get_eigenvectors, const char* method) {
	//out << real(H) << endl;
	arma::Mat<T> H_temp;
	try {
		H_temp = arma::Mat<T>(this->H);
		if (get_eigenvectors) arma::eig_sym(this->eigenvalues, this->eigenvectors, H_temp, method);
		else arma::eig_sym(this->eigenvalues, H_temp);
	}
	catch (...) {
		handle_exception(std::current_exception(), 
			"sparse - dim(H) = " + std::to_string(H.n_nonzero * sizeof(H(0, 0)))
			+ " bytes\ndense - dim(H) = " + std::to_string(H_temp.n_alloc * sizeof(H_temp(0, 0))) + " bytes"
		);
	}
	//for (long int i = 0; i < N; i++)
	//	this->eigenvectors.col(i) = arma::normalise(this->eigenvectors.col(i));

	double E_av = arma::trace(eigenvalues) / double(N);
	auto i = min_element(begin(eigenvalues), end(eigenvalues), [=](int x, int y) {
		return abs(x - E_av) < abs(y - E_av);
		});
	this->E_av_idx = i - begin(eigenvalues);
}

template <typename T> void IsingModel<T>::diagonalization_sparse(int num_of_eigvals, bool get_eigenvectors, const char* form){
	arma::eigs_opts opts;
	opts.maxiter =  10000;
	//if(get_eigenvectors) arma::eigs_sym(this->eigenvalues, this->eigenvectors, this->H, num_of_eigvals, form, opts);
	//else arma::eigs_sym(this->eigenvalues, this->H, num_of_eigvals, form, opts);
}
template <typename T> void IsingModel<T>::diagonalization_sparse(int num_of_eigvals, bool get_eigenvectors, double sigma){
	arma::eigs_opts opts;
	opts.maxiter =  10000;
	//if(get_eigenvectors) arma::eigs_sym(this->eigenvalues, this->eigenvectors, this->H, num_of_eigvals, sigma, opts);
	//else arma::eigs_sym(this->eigenvalues, this->H, num_of_eigvals, sigma, opts);
}


/// <summary>
/// calculates the total spin from the correlation matrix
/// </summary>
/// <param name="corr_mat"> spin correlation matrix </param>
/// <returns></returns>
template <typename T> double IsingModel<T>::total_spin(const arma::mat& corr_mat) {
	double S2 = arma::accu(corr_mat);
	if (S2 < -0.25) return 0;
	return (sqrt(1 + 4 * S2) - 1.0) / 2.0;
}

// ------------------------------------------------------------------------------------------------ ERGODIC QUANTITIES ------------------------------------------------------------------------------------------------

/// <summary>
/// The IPR, also called the participation ratio, quantifies how delocalized a state is in a certain basis.
/// The state is completely delocalized, when:
/// IPR=dim(hilbert space)
/// </summary>
/// <param name="state_idx"> index of the eigenvector used to calculate this quantity </param>
/// <returns> returns the IPR value </returns>
template <typename T> 
double IsingModel<T>::ipr(int state_idx) const {
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
template <typename T> 
double IsingModel<T>::information_entropy(u64 _id) const {
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
template <typename T> 
double IsingModel<T>::information_entropy(u64 _id, const IsingModel<T>& beta, u64 _min, u64 _max) const {
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
template <typename T> 
double IsingModel<T>::eigenlevel_statistics(u64 _min, u64 _max) const {
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
template <typename T> 
arma::vec IsingModel<T>::eigenlevel_statistics_with_return() const {
	arma::vec r(N - 2);
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


// ----------------------------------------------------------- ENTAGLEMENT -------------------------------------------------------
/// <summary>
/// Calculates the entropy of the system via the mixed density matrix
/// </summary>
/// <param name="state_id"> state index to produce the density matrix </param>
/// <param name="A_size"> size of subsystem </param>
/// <returns> entropy of considered systsem </returns>
template <typename T> 
double IsingModel<T>::entaglement_entropy(const arma::cx_vec& state, int A_size) const {
	auto rho = reduced_density_matrix(state, A_size);
	arma::vec probabilities;
	arma::eig_sym(probabilities, rho); //diagonalize to find probabilities and calculate trace in rho's eigenbasis
	double entropy = 0;
	//#pragma omp parallel for reduction(+: entropy)
	for (int i = 0; i < probabilities.size(); i++) {
		auto value = probabilities(i);
		entropy += (abs(value) < 1e-10) ? 0 : -value * log(abs(value));
	}
	//double entropy = -real(trace(rho * real(logmat(rho))));
	return entropy;
}
template <typename T> 
arma::vec IsingModel<T>::entaglement_entropy(const arma::cx_vec& state) const {
	arma::vec entropy(this->L - 1, arma::fill::zeros);
#pragma omp parallel for
	for (int i = 0; i < this->L - 1; i++)
		entropy(i) = entaglement_entropy(state, i + 1);
	return entropy;
}

template <typename _ty>
arma::Mat<_ty> matrix_pow(const arma::Mat<_ty>& matrix, int exponent) {
	if (exponent < 0)
		assert(false && "Support only posotive exponents");
	else if (exponent == 0) {
		auto X = arma::eye(matrix.n_rows, matrix.n_cols);
		return arma::cx_mat(X, X);
	}
	else if (exponent == 1)
		return matrix;
	else
		return matrix * matrix_pow(matrix, exponent - 1);
}
template <typename T> 
double IsingModel<T>::reyni_entropy(const arma::cx_vec& state, int A_size, unsigned alfa) const {
	assert(alfa > 1 && "Only alfa>=2 powers are possible");
	auto rho = reduced_density_matrix(state, A_size);
	return log2(real(trace(matrix_pow(rho, alfa)))) / (1.0 - alfa);
}

template <typename T> 
double IsingModel<T>::shannon_entropy(const arma::cx_vec& state, int A_size) const {
	auto rho = reduced_density_matrix(state, A_size);
	arma::vec probabilities;
	arma::eig_sym(probabilities, rho); //diagonalize to find probabilities and calculate trace in rho's eigenbasis
	double entropy = 0;
#pragma omp parallel for reduction(+: entropy)
	for (int i = 0; i < probabilities.size(); i++) {
		auto value = probabilities(i) * probabilities(i);
		entropy += (abs(value) < 1e-10) ? 0 : -value * log2(abs(value));
	}
	return entropy;
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
arma::vec IsingModel<T>::operator_av_in_eigenstates_return(double (IsingModel<T>::* op)(int, int), IsingModel<T>& A, int site) {
	arma::vec temp(A.get_hilbert_size(), arma::fill::zeros);
#pragma omp parallel for shared(temp)
	for (int k = 0; k < A.get_hilbert_size(); k++)
		temp(k) = (A.*op)(k, site);
	return temp;
}

//-------------------------------------------------------------------------------------------------LIOMs

template <typename _type>
arma::sp_cx_mat IsingModel<_type>::create_StringOperator(coordinate alfa, coordinate beta, int j, int ell) const {
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
	arma::sp_cx_mat SigmaAlfa = create_operator({ op_alfa }, { properSite(j) });
	arma::sp_cx_mat SigmaBeta = create_operator({ op_beta }, { properSite(j + ell) });
	arma::sp_cx_mat SigmaZstring;
	if (ell == 0) {
		if (alfa == beta && alfa == coordinate::y) return -create_operator({ IsingModel_sym::sigma_z });
		else assert(false && "Don't know what this could possibly be, check again if you need this");
	}
	else if (ell == 1) {
		SigmaZstring = arma::eye<arma::sp_cx_mat>(this->N, this->N);
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
arma::sp_cx_mat IsingModel<_type>::create_LIOMoperator_densities(int n, int j) const {
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
	return (arma::vec)spec_rep;
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
void IsingModel<_type>::set_coefficients(const arma::cx_vec& initial_state){
	this->coeff = arma::cx_vec(N, arma::fill::zeros);
	for (long k = 0; k < this->N; k++)
		this->coeff(k) = arma::dot(eigenvectors.col(k), initial_state);
}

template <typename _type>
auto IsingModel<_type>::time_evolve_state(const arma::cx_vec& state, double time)
	-> arma::cx_vec
{
	arma::cx_vec state_evolved(this->N, arma::fill::zeros);
	for (long k = 0; k < this->N; k++)
		state_evolved += std::exp(-im * eigenvalues(k) * time) * this->coeff(k) * eigenvectors.col(k);
	return arma::normalise(state_evolved);
}

template <typename _type>
void IsingModel<_type>::time_evolve_state_ns(
	arma::cx_vec& state,	//<! state at time dt
	double dt, 				//<! time step, dt << 1
	int order				//<! maximal order of expansion
	) {
	arma::cx_vec temp_state = state;
	for(int i = 1; i <= order; i++){
		const cpx prefactor = -im * dt /  double(i);
		temp_state = prefactor * this->H * temp_state;
		state += temp_state;
	}
	//state = arma::normalise(state);
}

// UNRESOLVED TEMPLATE EXTERNALS <- COMPILER DOESN"T KNOW ABOUT THEM SADLY
template IsingModel<double>::~IsingModel();
template IsingModel<cpx>::~IsingModel();
template void IsingModel<double>::set_neighbors();
template void IsingModel<cpx>::set_neighbors();
template void IsingModel<cpx>::diagonalization(bool, const char*);
template void IsingModel<double>::diagonalization(bool, const char*);
template void IsingModel<double>::diagonalization_sparse(int, bool, double);
template void IsingModel<cpx>::diagonalization_sparse(int, bool, double);
template void IsingModel<double>::diagonalization_sparse(int, bool, const char*);
template void IsingModel<cpx>::diagonalization_sparse(int, bool, const char*);
template double IsingModel<cpx>::eigenlevel_statistics(u64, u64) const;
template double IsingModel<double>::eigenlevel_statistics(u64, u64) const;
template arma::vec IsingModel<cpx>::eigenlevel_statistics_with_return() const;
template arma::vec IsingModel<double>::eigenlevel_statistics_with_return() const;
template double IsingModel<double>::ipr(int) const;
template double IsingModel<cpx>::ipr(int) const;
template double IsingModel<double>::information_entropy(u64) const;
template double IsingModel<cpx>::information_entropy(u64) const;
template double IsingModel<double>::information_entropy(u64, const IsingModel<double>&, u64, u64) const;
template double IsingModel<cpx>::information_entropy(u64, const IsingModel<cpx>&, u64, u64) const;
template double overlap(const IsingModel<double>&, const IsingModel<double>&, int, int);
template double IsingModel<cpx>::total_spin(const arma::mat&);
template double IsingModel<double>::total_spin(const arma::mat&);

template arma::sp_cx_mat IsingModel<cpx>::create_StringOperator(coordinate, coordinate, int,int) const;
template arma::sp_cx_mat IsingModel<double>::create_StringOperator(coordinate, coordinate, int, int) const;
template arma::sp_cx_mat IsingModel<cpx>::create_LIOMoperator_densities(int, int) const;
template arma::sp_cx_mat IsingModel<double>::create_LIOMoperator_densities(int, int) const;
template arma::cx_vec IsingModel<double>::time_evolve_state(const arma::cx_vec&, double);
template arma::cx_vec IsingModel<cpx>::time_evolve_state(const arma::cx_vec&, double);
template void IsingModel<double>::time_evolve_state_ns(arma::cx_vec&, double, int);
template void IsingModel<cpx>::time_evolve_state_ns(arma::cx_vec&, double, int);
template void IsingModel<double>::set_coefficients(const arma::cx_vec&);
template void IsingModel<cpx>::set_coefficients(const arma::cx_vec&);

template double IsingModel<double>::shannon_entropy(const arma::cx_vec&, int) const;
template double IsingModel<cpx>::shannon_entropy(const arma::cx_vec&, int) const;
template double IsingModel<double>::reyni_entropy(const arma::cx_vec&, int, unsigned) const;
template double IsingModel<cpx>::reyni_entropy(const arma::cx_vec&, int, unsigned) const;
template double IsingModel<double>::entaglement_entropy(const arma::cx_vec&, int) const;
template double IsingModel<cpx>::entaglement_entropy(const arma::cx_vec&, int) const;
template arma::vec IsingModel<double>::entaglement_entropy(const arma::cx_vec&) const;
template arma::vec IsingModel<cpx>::entaglement_entropy(const arma::cx_vec&) const;