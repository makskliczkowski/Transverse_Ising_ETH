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
	auto i = min_element(begin(eigenvalues), end(eigenvalues), [=](double x, double y) {
		return abs(x - E_av) < abs(y - E_av);
		});
	this->E_av_idx = i - begin(eigenvalues);
	printSeparated(std::cout, "\t", 16, true, "mean energy", "energies close to this value (-1,0,+1 around found index");
	printSeparated(std::cout, "\t", 16, true, E_av, eigenvalues(this->E_av_idx - 1), eigenvalues(this->E_av_idx),  eigenvalues(this->E_av_idx + 1));
}

// ----------------------------------------------------------- OPERATORS AND AVERAGES -------------------------------------------------------


//------------------------------------------------------------------------------------------------- TFIM LIOMs OPERATORS
template <typename _type>
arma::sp_cx_mat IsingModel<_type>::create_tfim_liom_plus(int n) const {
	if(n == 0) return cast_cx_sparse(this->H);
	const u64 size = ULLPOW(this->L);
	arma::sp_cx_mat tfim_liom(size, size);
	std::vector<int> sites;
	for(u64 kk = 0; kk < size; kk++){
		u64 k = map(kk);
		for(int j = 0; j < this->L; j++){

			// S^zz_{j, j+n+1}
			auto [right_value, idx0] = IsingModel::sigma_z(
				k, this->L, 
						{ this->properSite(j + n + 1) }); 
			sites = std::vector<int>(n); iota(sites.begin(), sites.end(), j + 1);
			auto [middle_value, idx1] = IsingModel::sigma_x(idx0, this->L, this->properSite(sites));
			auto [left_value, idx] = IsingModel::sigma_z(idx1, this->L, { this->properSite(j) } );
			tfim_liom(idx, k) += this->J * left_value * middle_value * right_value;

			// S^yy_{j,j+n-1}
			if(n > 1){
				std::tie(right_value, idx0) = IsingModel::sigma_y(k, this->L, { this->properSite(j + n - 1) });
				middle_value = 1.0;
				idx1 = idx0;
				if(n > 2){
					sites = std::vector<int>(n - 2); iota(sites.begin(), sites.end(), j + 1);
					std::tie(middle_value, idx1) = IsingModel::sigma_x(idx0, this->L, this->properSite(sites));
				}
				std::tie(left_value, idx) = IsingModel::sigma_y(idx1, this->L, { this->properSite(j) });
				tfim_liom(idx, k) += this->J * left_value * middle_value * right_value;
			} else {
				std::tie(middle_value, idx) = IsingModel::sigma_x(k, this->L, { this->properSite(j) });
				tfim_liom(idx, k) -= this->J * middle_value;
			}

			// S^zz_{j, j+n}
			std::tie(right_value, idx0) = IsingModel::sigma_z(
				k, this->L, 
						{ this->properSite(j + n) });
			idx1 = idx0;
			if(n > 1){
				sites = std::vector<int>(n - 1); iota(sites.begin(), sites.end(), j + 1);
				std::tie(middle_value, idx1) = IsingModel::sigma_x(idx0, this->L, this->properSite(sites));
			}
			std::tie(left_value, idx) = IsingModel::sigma_z(
						idx1, this->L, 
						{ this->properSite(j) });
			tfim_liom(idx, k) -= this->g * left_value * middle_value * right_value;

			// S^yy_{j, j+n}
			std::tie(right_value, idx0) = IsingModel::sigma_y(
				k, this->L, 
						{ this->properSite(j + n) });
			idx1 = idx0;
			if(n > 1){
				sites = std::vector<int>(n - 1); iota(sites.begin(), sites.end(), j + 1);
				std::tie(middle_value, idx1) = IsingModel::sigma_x(idx0, this->L, this->properSite(sites));
			}
			std::tie(left_value, idx) = IsingModel::sigma_y(
						idx1, this->L, 
						{ this->properSite(j) });
			tfim_liom(idx, k) -= this->g * left_value * middle_value * right_value;
		}
	}
	return tfim_liom;
}

template <typename _type>
arma::sp_cx_mat IsingModel<_type>::create_tfim_liom_minus(int n) const {
	const u64 size = ULLPOW(this->L);
	arma::sp_cx_mat tfim_liom(size, size);
	std::vector<int> sites;
	for(u64 k = 0; k < size; k++){
		for(int j = 0; j < this->L; j++){
			
			// S^yz_{j, j+n+1}
			auto [right_value, idx0] = IsingModel::sigma_z(
				k, this->L, 
						{ this->properSite(j + n + 1) }
						);
			cpx middle_value = 1.0;
			if(n > 0){
				sites = std::vector<int>(n); iota(sites.begin(), sites.end(), j + 1);
			 	std::tie(middle_value, idx0) = IsingModel::sigma_x(k, this->L, this->properSite(sites));
			}
			auto [left_value, idx] = IsingModel::sigma_y(idx0, this->L, { this->properSite(j) } );
			tfim_liom(idx, k) += this->J * left_value * middle_value * right_value;

			// S^zy_{j, j+n+1}
			std::tie(right_value, idx0) = IsingModel::sigma_y(
				k, this->L, 
						{ this->properSite(j + n + 1) }
						);
			if(n > 0){
				sites = std::vector<int>(n); iota(sites.begin(), sites.end(), j + 1);
			 	std::tie(middle_value, idx0) = IsingModel::sigma_x(idx0, this->L, this->properSite(sites));
			}
			std::tie(left_value, idx) = IsingModel::sigma_z(idx0, this->L, { this->properSite(j) } );
			tfim_liom(idx, k) -= this->J * left_value * middle_value * right_value;
		}
	}
	return tfim_liom;
}


template <typename _type>
arma::sp_cx_mat IsingModel<_type>::spin_imbalance() const {
	arma::cx_vec imbal(this->L, arma::fill::zeros);
	for(int i = 0; i < this->L; i++){
		imbal(i) = (i % 2 == 0)? 1. : -1.;
	}
	return this->create_operator({IsingModel::sigma_z}, imbal);			
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
template double overlap(const IsingModel<double>&, const IsingModel<double>&, int, int);

template arma::sp_cx_mat IsingModel<cpx>::create_tfim_liom_plus(int) const;
template arma::sp_cx_mat IsingModel<double>::create_tfim_liom_plus(int) const;
template arma::sp_cx_mat IsingModel<cpx>::create_tfim_liom_minus(int) const;
template arma::sp_cx_mat IsingModel<double>::create_tfim_liom_minus(int) const;
template arma::sp_cx_mat IsingModel<double>::spin_imbalance() const;
template arma::sp_cx_mat IsingModel<cpx>::spin_imbalance() const;
template arma::cx_vec IsingModel<double>::time_evolve_state(const arma::cx_vec&, double);
template arma::cx_vec IsingModel<cpx>::time_evolve_state(const arma::cx_vec&, double);
template void IsingModel<double>::time_evolve_state_ns(arma::cx_vec&, double, int);
template void IsingModel<cpx>::time_evolve_state_ns(arma::cx_vec&, double, int);
template void IsingModel<double>::set_coefficients(const arma::cx_vec&);
template void IsingModel<cpx>::set_coefficients(const arma::cx_vec&);