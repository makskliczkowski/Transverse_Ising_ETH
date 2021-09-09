#include "include/IsingModel.h"
template<typename T> IsingModel<T>::~IsingModel() {
}

// ---- PRINTERS
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
	this->next_nearest_neighbors = std::vector<int>(L,0);
	switch (this->_BC) {
	case 0:
		// periodic boundary conditions
		for(int i = 0; i < this->L; i++){
			this->nearest_neighbors[i] = (i+1) % this->L;
			this->next_nearest_neighbors[i] = (i+2) % this->L;
		}
		break;

	case 1:
		// open boundary conditions
		for(int i = 0; i < this->L; i++){
			this->nearest_neighbors[i] = (i+1) % this->L;
			this->next_nearest_neighbors[i] = (i+2) % this->L;
		}
		this->nearest_neighbors[L-1] = -1;
		this->next_nearest_neighbors[L-2] = -1;
		this->next_nearest_neighbors[L-1] = -1;
		break;
	default:
		for(int i = 0; i < this->L; i++){
			this->nearest_neighbors[i] = (i+1) % this->L;
			this->next_nearest_neighbors[i] = (i+2) % this->L;
		}
		break;
	}
}

/// <summary>
/// General procedure to diagonalize the Hamiltonian using eig_sym from the Armadillo library
/// </summary>
template <typename T> void IsingModel<T>::diagonalization() {
	//out << real(H) << endl;
	try {
		arma::eig_sym(eigenvalues, eigenvectors, H);
	}
	catch (const bad_alloc& e) {
		stout << "Memory exceeded" << e.what() << "\n";
		stout << "dim(H) = " << H.size() * sizeof(H(0, 0)) << "\n";
		assert(false);
	}
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
	return (sqrt(1 + 4 * arma::accu(corr_mat)) - 1.0) / 2.0;
}

/// <summary>
/// The IPR, also called the participation ratio, quantifies how delocalized a state is in a certain basis.
/// The state is completely delocalized, when:
/// IPR=dim(hilbert space)
/// </summary>
/// <param name="state_idx"> index of the eigenvector used to calculate this quantity </param>
/// <returns> returns the IPR value </returns>
template <typename T> double IsingModel<T>::ipr(int state_idx) {
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
///
/// </summary>
/// <param name="_id"></param>
/// <returns></returns>
template <typename T> double IsingModel<T>::information_entropy(u64 _id) {
	arma::subview_col state = this->eigenvectors.col(_id);
	double ent = 0;
#pragma omp parallel for reduction(+: ent)
	for (int k = 0; k < N; k++) {
		double val = abs(conj(state(k)) * state(k));
		ent += val * log(val);
	}
	return -ent / log(0.48 * this->N);
}

/// <summary>
/// 
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="_id"></param>
/// <param name="beta"></param>
/// <returns></returns>
template <typename T> double IsingModel<T>::information_entropy(u64 _id, const IsingModel<T>& beta, u64 _min, u64 _max) {
	arma::subview_col state_alfa = this->eigenvectors.col(_id);
	double ent = 0;
//#pragma omp parallel for reduction(+: ent)
	for (int k = _min; k < _max; k++) {
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
template <typename T> double IsingModel<T>::eigenlevel_statistics(u64 _min, u64 _max) {
	double r = 0;
	if (_min <= 0) throw "too low index";
	if (_max >= N) throw "index exceeding Hilbert space";
	//double delta_n = eigenvalues(_min) - eigenvalues(_min - 1);
	//double delta_n_next = 0;
#pragma omp parallel for reduction(+: r)
	for (int k = _min; k < _max; k++) {
		const double delta_n = eigenvalues(k) - eigenvalues(k - 1);
		const double delta_n_next = eigenvalues(k + 1) - eigenvalues(k);
		const double min = std::min(delta_n, delta_n_next);
		const double max = std::max(delta_n, delta_n_next);
		if (abs(delta_n) <= 1e-14) continue;// throw "Degeneracy!!!\n";
		r += min / max;
		//delta_n = delta_n_next;
	}
	return r / double(_max - _min);
}
template <typename T> vec IsingModel<T>::eigenlevel_statistics_with_return() {
	vec r(N - 2);
#pragma omp parallel for shared(r)
	for (int k = 1; k < N - 1; k++)
		r(k - 1) = eigenlevel_statistics(k, k + 1);
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

// - - - - - OPERATOR AVERAGES - - - - -

// -----------------------------------------------------------> TO REFACTOR AND CREATE DESCRIPTION -------------------------------------------------------

/// <summary>
/// Prints to file the average of a given operator in each eigenstate as a function of eigenenergies
/// </summary>
/// <param name="op"> operator as function acting on a specific eigenstate </param>
/// <param name="A"> class instance </param>
/// <param name="site"> position of the spin, where the operator is acted upon </param>
/// <param name="name"> name of the file to print data </param>
/// <param name="separator"> separator between columns in file </param>
template <typename T>
void IsingModel<T>::operator_av_in_eigenstates(double (IsingModel<T>::* op)(int, int), \
	IsingModel<T>& A, int site, std::string name, std::string separator) {
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
template <typename T>
vec IsingModel<T>::operator_av_in_eigenstates_return(double (IsingModel<T>::* op)(int, int), IsingModel<T>& A, int site) {
	vec temp(A.get_hilbert_size(), fill::zeros);
#pragma omp parallel for shared(temp)
	for (int k = 0; k < A.get_hilbert_size(); k++)
		temp(k) = (A.*op)(k, site);
	return temp;
}

// - - - - - STATISTICS AND PROBABILITIES - - - - - 
/// <summary>
///
/// </summary>
/// <param name="dir"></param>
/// <param name="name"></param>
/// <param name="data"></param>
/// <param name="_min"></param>
/// <param name="_max"></param>
/// <param name="step"></param>
void probability_distribution(std::string dir, std::string name, const arma::vec& data, double _min, double _max, double step) {
	std::ofstream file(dir + name + ".dat");
	int size = static_cast<int>((_max - _min) / step + 1);
	arma::vec prob_dist(size, arma::fill::zeros);
	for (int k = 1; k < data.size(); k++) {
		if (data(k) > _min && data(k) < _max) {
			const int bucket = static_cast<int>((data(k) + abs(_min)) / step);
			// out << "data(k) + _min: " << data(k) + _min<< "bucket: " << bucket << std::endl;
			prob_dist(bucket) += 1;
		}
	}
	prob_dist = normalise_dist(prob_dist, _min, _max);
	for (int p = 0; p < size; p++)
		file << p * step + _min << "\t" << prob_dist(p) << std::endl;
	file.close();
}
arma::vec probability_distribution_with_return(const arma::vec& data, double _min, double _max, double step) {
	int size = static_cast<int>((_max - _min) / step + 1);
	arma::vec prob_dist(size, arma::fill::zeros);
	for (int k = 1; k < data.size(); k++) {
		if (data(k) > _min && data(k) < _max) {
			const int bucket = static_cast<int>((data(k) + abs(_min)) / step);
			prob_dist(bucket) += 1;
		}
	}
	return normalise_dist(prob_dist, _min, _max);
}
/// <summary>
///
/// </summary>
/// <param name="data"></param>
/// <param name="mu"></param>
/// <returns></returns>
arma::vec data_fluctuations(const arma::vec& data, int mu) {
	arma::vec fluct(data.size() - mu, arma::fill::zeros);
	assert(mu < data.size() && "Bucket exceeds data container\nTry again\n");
	int end = data.size() - mu / 2.;
#pragma omp parallel for shared(fluct, end, mu, data)
	for (int k = mu / 2.; k < end; k++) {
		double average = 0;
		for (int n = k - mu / 2; n < k + mu / 2; n++)
			average += data(n);
		fluct(k - mu / 2.) = data(k) - average / double(mu);
	}
	return fluct;
}

/// <summary>
///
/// </summary>
/// <param name="data"></param>
/// <returns></returns>
arma::vec statistics_average(const arma::vec& data, int num_of_outliers) {
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
	spec_rep[0] = average / double(data.size() - 2);
	return (vec)spec_rep;
}

// - - - - - TOOLS - - - - -

/// <summary>
/// Overlapping of two eigenstates of possibly different matrices A and B
/// </summary>
/// <param name="A">matrix A</param>
/// <param name="B">matrix B</param>
/// <param name="n_a">number of A eigenstate</param>
/// <param name="n_b">number of B eigenstate</param>
/// <returns>A_n_a dot n_b_B</returns>
template <typename T>
T overlap(const IsingModel<T>& A, const IsingModel<T>& B, int n_a, int n_b){
	if (A.get_hilbert_size() != B.get_hilbert_size()) throw "Incompatible Hilbert dimensions\n";
	if(n_a >= A.get_hilbert_size() || n_b >= B.get_hilbert_size() || n_a < 0 || n_b < 0) throw "Cannot create an overlap between non-existing states\n";
	return arma::cdot(A.get_eigenState(n_a), B.get_eigenState(n_b));
}
template<> cpx overlap<cpx>(const IsingModel<cpx>& A, const IsingModel<cpx>& B, int n_a, int n_b) {
	if (A.get_hilbert_size() != B.get_hilbert_size()) throw "Incompatible Hilbert dimensions\n";
	if (n_a >= A.get_hilbert_size() || n_b >= B.get_hilbert_size() || n_a < 0 || n_b < 0) throw "Cannot create an overlap between non-existing states\n";
	double overlap_real = 0, overlap_imag = 0;
	auto state_A = A.get_eigenState(n_a);
	auto state_B = B.get_eigenState(n_b);
	arma::normalise(state_A);
	arma::normalise(state_B);
#pragma omp parallel for reduction(+: overlap_real, overlap_imag)
	for (int k = 0; k < A.get_hilbert_size(); k++) {
		cpx over = conj(state_A(k)) * state_B(k);
		overlap_real += real(over);
		overlap_imag += imag(over);
	}
	return cpx(overlap_real, overlap_imag);
}




// SOME SHIT WE HAVE TO DO
template IsingModel<double>::~IsingModel();
template IsingModel<cpx>::~IsingModel();
template void IsingModel<double>::set_neighbors();
template void IsingModel<cpx>::set_neighbors();
template void IsingModel<cpx>::diagonalization();
template void IsingModel<double>::diagonalization();
template double IsingModel<cpx>::eigenlevel_statistics(u64, u64);
template double IsingModel<double>::eigenlevel_statistics(u64, u64);
template vec IsingModel<cpx>::eigenlevel_statistics_with_return();
template vec IsingModel<double>::eigenlevel_statistics_with_return();
template double IsingModel<double>::ipr(int);
template double IsingModel<cpx>::ipr(int);
template double IsingModel<double>::information_entropy(u64);
template double IsingModel<cpx>::information_entropy(u64);
template double IsingModel<double>::information_entropy(u64, const IsingModel<double>&, u64, u64);
template double IsingModel<cpx>::information_entropy(u64, const IsingModel<cpx>&, u64, u64);
template cpx overlap(const IsingModel<cpx>&, const IsingModel<cpx>&, int, int);
template double overlap(const IsingModel<double>&, const IsingModel<double>&, int, int);
