#include "include/IsingModel.h"

IsingModel::~IsingModel() {
}
// - - - - - GETTERS and SETTERS - - - - -
std::string IsingModel::get_info() const {
	return this->info;
}
u64 IsingModel::get_hilbert_size() const {
	return this->N;
}
const cx_mat& IsingModel::get_hamiltonian() const {
	return this->H;
}
const vec& IsingModel::get_eigenvalues() const {
	return this->eigenvalues;
}
const cx_mat& IsingModel::get_eigenvectors() const {
	return this->eigenvectors;
}
const std::vector<u64>& IsingModel::get_mapping() const {
	return this->mapping;
}
/// <summary>
/// Return given eigenenergy at index idx
/// </summary>
/// <param name="idx"> index of eigenvalue </param>
/// <returns> eigenvalue at idx </returns>
double IsingModel::get_eigenEnergy(u64 idx) const
{
	return this->eigenvalues(idx);
}
/// <summary>
/// Returns given eigenstate at index idx
/// </summary>
/// <param name="idx"> index of eigenstate </param>
/// <returns> eigenstate at idx </returns>
const cx_vec& IsingModel::get_eigenState(u64 idx) const {
	return this->eigenvectors.col(idx);
}

// ---- PRINTERS
/// <summary>
/// prints the basis vector in the given spin-symmetry block
/// </summary>
/// <param name="Sz"> custom spin symmetry block </param>
void IsingModel::print_base_spin_sector(int Sz) {
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
void IsingModel::print_state(u64 _id) {
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
void IsingModel::set_neighbors() {
	this->nearest_neighbors = std::vector<int>(L);
	switch (this->_BC) {
	case 0:
		// periodic boundary conditions
		std::iota(nearest_neighbors.begin(), nearest_neighbors.end(), 1);
		nearest_neighbors[L - 1] = 0;
		break;

	case 1:
		// open boundary conditions
		std::iota(nearest_neighbors.begin(), nearest_neighbors.end(), 1);
		nearest_neighbors[L - 1] = -1;
		break;
	default:
		std::iota(nearest_neighbors.begin(), nearest_neighbors.end(), 1);
		nearest_neighbors[L - 1] = 0;
		break;
	}
}

/// <summary>
/// General procedure to diagonalize the Hamiltonian using eig_sym from the Armadillo library
/// </summary>
void IsingModel::diagonalization() {
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
double IsingModel::total_spin(const mat& corr_mat) {
	return (sqrt(1 + 4 * arma::accu(corr_mat)) - 1.0) / 2.0;
}

/// <summary>
/// The IPR, also called the participation ratio, quantifies how delocalized a state is in a certain basis.
/// The state is completely delocalized, when:
/// IPR=dim(hilbert space)
/// </summary>
/// <param name="state_idx"> index of the eigenvector used to calculate this quantity </param>
/// <returns> returns the IPR value</returns>
double IsingModel::ipr(int state_idx) {
	double ipr = 0;
	cx_vec state = eigenvectors.col(state_idx);
#pragma omp parallel for reduction(+: ipr)
	for (int n = 0; n < N; n++) {
		double value = abs(state(n) * state(n));
		ipr += value * value;
	}
	return 1.0 / ipr;
}

/// <summary>
///
/// </summary>
/// <param name="_id"></param>
/// <returns></returns>
double IsingModel::information_entropy(const u64 _id) {
	cx_vec state = this->eigenvectors.col(_id);
	double ent = 0;
#pragma omp parallel for reduction(+: ent)
	for (int k = 0; k < N; k++) {
		double val = abs(conj(state(k)) * state(k));
		ent += val * log(val);
	}
	return -ent / log(0.48 * N);
}

/// <summary>
/// Calculates the energy-level statistics within the energy window denoted by the indices _min and _max
/// computed as in: PHYSICAL REVIEW B 91, 081103(R) (2015)
/// </summary>
/// <param name="_min"> index of eigenenergy, being the lower bound of energy window </param>
/// <param name="_max"> index of eigenenergy, being the upper bound of energy window </param>
/// <returns></returns>
double IsingModel::eigenlevel_statistics(u64 _min, u64 _max) {
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
		if (max == 0) throw "Degeneracy!!!\n";
		r += min / max;
		//delta_n = delta_n_next;
	}
	return r / double(_max - _min);
}
vec IsingModel::eigenlevel_statistics_with_return() {
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
double IsingModel::spectrum_repulsion(double (IsingModel::* op)(int, int), IsingModel& A, int site) {
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
///
/// </summary>
/// <param name="_min"></param>
/// <param name="_max"></param>
/// <param name="symmetry"></param>
/// <param name="original"></param>
/// <returns></returns>
std::unordered_map<u64, u64> mapping_sym_to_original(u64 _min, u64 _max, const IsingModel& symmetry, const IsingModel& original) {
	std::unordered_map<u64, u64> map;
	std::vector<double> E_dis = arma::conv_to<std::vector<double>>::from(original.eigenvalues);
	for (int k = 0; k < symmetry.N; k++) {
		double E = symmetry.eigenvalues(k);
		if (E < original.eigenvalues(_min) && E >= original.eigenvalues(_max)) continue;
		auto idx = binary_search(E_dis, _min, _max, E);
		if (idx < 0 || idx >= original.N) continue;
		double E_prev = (idx == 0) ? (original.eigenvalues(0) - 1.0) : original.eigenvalues(idx - 1);
		double E_next = (idx == original.N - 1) ? (original.eigenvalues(original.N - 1) + 1.0) : original.eigenvalues(idx + 1);
		if (abs(E - E_prev) > 1e-8 && abs(E - E_next) > 1e-8)
			map[k] = idx;
	}
	return map;
}

// - - - - - OPERATOR AVERAGES - - - - -

// -----------------------------------------------------------> TO REFACTOR AND CREATE DESCRIPTION -------------------------------------------------------


/// <summary>
///
/// </summary>
/// <param name="site"></param>
/// <param name="beta"></param>
/// <param name="alfa"></param>
/// <param name="sector_alfa"></param>
/// <param name="sector_alfa"></param>
/// <returns></returns>
double av_sigma_x_sym_sectors(int site, const u64 beta, const u64 alfa, const IsingModel_sym& sector_alfa, const IsingModel_sym& sector_beta) {
	assert(sector_beta.L == sector_alfa.L && "incompatible chain lengths");
	const int L = sector_beta.L;
	const u64 N_alfa = sector_alfa.N;
	const u64 N_beta = sector_beta.N;
	const arma::subview_col<cpx>& state_beta = sector_beta.eigenvectors.col(beta);
	const arma::subview_col<cpx>& state_alfa = sector_alfa.eigenvectors.col(alfa);
	std::vector<bool> base(L);
	cpx overlap = 0;
	for (long int k = 0; k < N_beta; k++) {
		if (abs(state_beta(k)) < 1e-12) continue;
		int_to_binary(sector_beta.mapping[k], base);
		base[site] = !base[site];
		u64 idx = binary_search(sector_alfa.mapping, 0, N_alfa - 1, binary_to_int(base));
		int sym_eig = 1;
		if (idx == -1) {
			auto tup_T = sector_alfa.find_translation_representative(base);
			auto tup_S = sector_alfa.find_SEC_representative(base);
			auto [min, trans_eig] = (std::get<0>(tup_T) > std::get<0>(tup_S)) ? tup_S : tup_T;
			sym_eig = trans_eig;
			idx = binary_search(sector_alfa.mapping, 0, N_alfa - 1, min);
			if (idx < 0 || idx >= N_alfa) continue; // min is not present in symmetry sector
		}
		cpx translation_eig = (abs(sym_eig) == 1 || sector_alfa.symmetries.k_sym == 0) ? \
			cpx(1.0) : std::exp(-1i * sector_alfa.symmetries.k_sym * double(abs(sym_eig) - 1));
		cpx value_new = double(sgn(sym_eig)) * translation_eig * (sector_alfa.normalisation[idx] / sector_beta.normalisation[k]);
		overlap += conj(state_alfa(idx)) * value_new * state_beta(k);
	}
	return real(overlap);
}

/// <summary>
///
/// </summary>
/// <param name="n"></param>
/// <param name="m"></param>
/// <returns></returns>
double IsingModel::av_sigma_z_extensive(u64 alfa, u64 beta) {
	std::vector<bool> base_vector(L);
	cpx overlap = 0.0;
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);

	//stout << "n="<< real(state_n.st()) << "m=" << real(state_m.st());
	for (long int k = 0; k < N; k++) {
		int_to_binary(map(k), base_vector);
		double val = 0;
		for (int j = 0; j < this->L; j++) {
			val += base_vector[j] ? 1.0 : -1.0;
		}
		//stout << val << endl;
		overlap += conj(state_alfa(k)) * val * state_beta(k);
	}
	return real(overlap) / double(this->L * this->L);
}

/// <summary>
///
/// </summary>
/// <param name="n"></param>
/// <param name="m"></param>
/// <returns></returns>
double IsingModel::av_sigma_z_extensive_corr(u64 alfa, u64 beta, int corr_len) {
	std::vector<bool> base_vector(L);
	cpx overlap = 0.0;
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);

	std::vector<int> corr_nei(this->L, 0);
	for (int k = 0; k < this->L; k++)
	{
		if (this->_BC == 1) {
			if (k + corr_len > this->L)
				corr_nei[k] = -1;																				// OBC
			else
				corr_nei[k] = k + corr_len;
		}
		else
			corr_nei[k] = (k + corr_len) % this->L;																// PBC
	}
	//stout << "n="<< real(state_n.st()) << "m=" << real(state_m.st());
	for (long int k = 0; k < N; k++) {
		int_to_binary(map(k), base_vector);
		double val = 0;
		for (int j = 0; j < this->L; j++) {
			double S_i = base_vector[j] ? 1 : -1;
			double S_i_corr = 0;
			if (corr_nei[j] >= 0) S_i_corr = double(base_vector[corr_nei[j]] ? 1 : -1);
			val += S_i * S_i_corr;
		}
		//stout << val << endl;
		overlap += conj(state_alfa(k)) * val * state_beta(k);
	}
	return real(overlap) / double(this->L * this->L);
}

/// <summary>
/// Prints to file the average of a given operator in each eigenstate as a function of eigenenergies
/// </summary>
/// <param name="op"> operator as function acting on a specific eigenstate </param>
/// <param name="A"> class instance </param>
/// <param name="site"> position of the spin, where the operator is acted upon </param>
/// <param name="name"> name of the file to print data </param>
/// <param name="separator"> separator between columns in file </param>
void IsingModel::operator_av_in_eigenstates(double (IsingModel::* op)(int, int), \
	IsingModel& A, int site, std::string name, std::string separator) {
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
vec IsingModel::operator_av_in_eigenstates_return(double (IsingModel::* op)(int, int), IsingModel& A, int site) {
	vec temp(A.get_hilbert_size(), fill::zeros);
#pragma omp parallel for shared(temp)
	for (int k = 0; k < A.get_hilbert_size(); k++)
		temp(k) = (A.*op)(k, site);
	return temp;
}

/// <summary>
/// Calculates the quantum fidelity of two eigenstates with slightly different parameters
/// </summary>
/// <param name="_min"> lower bound of averaging bucket (over energy) </param>
/// <param name="_max"> upper bound of averaging bucket (over energy) </param>
/// <param name="Hamil"> Hamiltonian with given set of parameters </param>
/// <param name="J"> new exchange integral </param>
/// <param name="g"> new trasnverse field </param>
/// <param name="h"> new uniform perpendicular field </param>
/// <param name="w"> new disorder stregth </param>
/// <returns></returns>
double quantum_fidelity(u64 _min, u64 _max, const IsingModel& Hamil, double J, double J0, double g, double g0, double h, double w) {
	if (_min < 0) throw "too low index";
	if (_max >= Hamil.get_hilbert_size()) throw "index exceeding Hilbert space";

	std::unique_ptr<IsingModel> Hamil2;
	if (typeid(Hamil) == typeid(IsingModel_disorder)) Hamil2 = std::make_unique<IsingModel_disorder>(Hamil.L, J, J0, g, g0, h, w);
	else Hamil2 = std::make_unique<IsingModel_sym>(Hamil.L, J, g, h);

	double fidelity = 0;
#pragma omp parallel for reduction(+: fidelity)
	for (long int k = _min; k < _max; k++)
		fidelity += overlap(Hamil, *Hamil2, k, k);
	return fidelity / double(_max - _min);
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
	arma::normalise(prob_dist);
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
			// out << "data(k) + _min: " << data(k) + _min<< "bucket: " << bucket << std::endl;
			prob_dist(bucket) += 1;
		}
	}
	return arma::normalise(prob_dist);;
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
double overlap(const IsingModel & A, const IsingModel & B, int n_a, int n_b)
{
	if(n_a >= A.N || n_b >= B.N || n_a < 0 || n_b < 0) throw "Cannot create an overlap between non-existing states\n";
	return abs(arma::cdot(A.eigenvectors.col(n_a), B.eigenvectors.col(n_b)));
}
