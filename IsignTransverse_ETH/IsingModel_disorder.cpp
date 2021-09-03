#include "include/IsingModel.h"

/* CONSTRUCTORS */
IsingModel_disorder::IsingModel_disorder(int L, double J, double J0, double g, double g0, double h, double w, int _BC) {
	this->L = L; this->J = J; this->g = g; this->h = h;
	this->J0 = J0; this->g0 = g0;  this->w = w;
	this->N = static_cast<u64>(std::pow(2, L));
	this->_BC = _BC;

	//change info
	this->info = "_L=" + std::to_string(this->L) + \
		",J0=" + to_string_prec(this->J0, 2) + \
		",g=" + to_string_prec(this->g, 2) + \
		",g0=" + to_string_prec(this->g0, 2) + \
		",h=" + to_string_prec(this->h, 2) + \
		",w=" + to_string_prec(this->w, 2);
	set_neighbors();
	hamiltonian();
}
/*
IsingModel_disorder::IsingModel_disorder(const IsingModel_disorder& A) {
	this->L = A.L; this->J = A.J; this->g = A.g; this->h = A.h; this->disorder_strength = A.disorder_strength;
	this->N = A.N; this->mapping = A.mapping;
	this->H = A.H; this->eigenvectors = A.eigenvectors; this->eigenvalues = A.eigenvalues;
}
IsingModel_disorder::IsingModel_disorder(IsingModel_disorder&& A) noexcept {
	this->L = A.L; this->J = A.J; this->g = A.g; this->h = A.h; this->disorder_strength = A.disorder_strength;
	this->N = A.N; this->mapping = A.mapping;
	this->H = A.H; this->eigenvectors = A.eigenvectors; this->eigenvalues = A.eigenvalues;
}
IsingModel_disorder::~IsingModel_disorder()
{
	//out << "Destroying the disordered Ising model\n";
}
*/

/* BASE GENERATION AND RAPPING*/
/// <summary>
/// Return the index in the case of no mapping in disorder
/// </summary>
/// <param name="index"> index to take</param>
/// <returns>index</returns>
u64 IsingModel_disorder::map(u64 index) {
	if (index < 0 || index >= std::pow(2, L)) throw "Element out of range\n No such index in map\n";
	return index;
}
void IsingModel_disorder::generate_mapping() {
	this->mapping = std::vector<u64>();
	std::vector<bool> base_vector(L);
	int_to_binary(u64(std::pow(2, L / 2) - 1), base_vector);
	while (next_permutation(base_vector.begin(), base_vector.end()))
		this->mapping.push_back(binary_to_int(base_vector));
	this->N = this->mapping.size();
}

//-------------------------------------------------------------------------------------------------------------------------------
// BUILDING HAMILTONIAN
/// <summary>
/// Sets the non-diagonal elements of the Hamimltonian matrix, by acting with the operator on the k-th state
/// </summary>
/// <param name="k"> index of the basis state acted upon with the Hamiltonian </param>
/// <param name="value"> value of the given matrix element to be set </param>
/// <param name="temp"> resulting vector form acting with the Hamiltonian operator on the k-th basis state </param>
void IsingModel_disorder::setHamiltonianElem(u64 k, double value, std::vector<bool>& temp) {
	u64 idx = binary_to_int(temp);
	H(idx, k) += value;
}
/// <summary>
/// Generates the total Hamiltonian of the system. The diagonal part is straightforward,
/// while the non-diagonal terms need the specialized setHamiltonainElem(...) function
/// </summary>
void IsingModel_disorder::hamiltonian() {
	try {
		this->H = cx_mat(N, N, fill::zeros);                                //  hamiltonian
	}
	catch (const bad_alloc& e) {
		std::cout << "Memory exceeded" << e.what() << "\n";
		assert(false);
	}

	this->dh = create_random_vec(L, this->w);                               // creates random disorder vector
	this->dJ = create_random_vec(L, this->J0);                              // creates random exchange vector
	this->dg = create_random_vec(L, this->g0);                              // creates random transverse field vector

	std::vector<bool> base_vector(L);
	std::vector<bool> temp(base_vector);                                    // changes under H action
	for (long int k = 0; k < N; k++) {
		int_to_binary(k, base_vector);                                      // check state number
		double s_i, s_j;
		for (int j = 0; j <= L - 1; j++) {
			s_i = base_vector[j] ? 1.0 : -1.0;                              // true - spin up, false - spin down
			/* transverse field */
			temp = base_vector;
			temp[j] = !base_vector[j];                                      // negates on that site
			setHamiltonianElem(k, this->g + this->dg(j), temp);
			/* disorder */
			H(k, k) += (this->h + dh(j)) * s_i;                             // diagonal

			if (nearest_neighbors[j] >= 0) {
				/* Ising-like spin correlation */
				s_j = base_vector[nearest_neighbors[j]] ? 1.0 : -1.0;
				this->H(k, k) += (this->J + this->dJ(j)) * s_i * s_j;
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------------------------------------
/* PHYSICAL QUANTITES */

/// <summary>
///
/// </summary>
/// <param name="state_id"></param>
/// <param name="site"></param>
/// <returns></returns>
double IsingModel_disorder::av_sigma_z(int site, u64 alfa, u64 beta) {
	if (site < 0 || site >= L) throw "Site index exceeds chain";
	arma::subview_col state_alfa = this->eigenvectors.col(alfa);
	arma::subview_col state_beta = this->eigenvectors.col(beta);

	cpx value = 0;
	std::vector<bool> base_vector(L);
	for (int k = 0; k < N; k++) {
		int_to_binary(map(k), base_vector);
		double Sz = base_vector[site] ? 1.0 : -1.0;
		value += Sz * conj(state_alfa(k)) * state_beta(k);
	}
	return real(value);
}

/// <summary>
/// Calculates the matrix element of the x-component spin matrix within the eigenstate state_id
/// </summary>
/// <param name="state_id"> index of eigenstate </param>
/// <param name="site"> position of the spin, where the operator is acted upon </param>
/// <returns> return the average in the eigenstate </returns>
double IsingModel_disorder::av_sigma_x(int site, u64 alfa, u64 beta) {
	const vec state_alfa = real(eigenvectors.col(alfa));
	const vec state_beta = real(eigenvectors.col(beta));
	double value = 0;
#pragma omp parallel
	{
		std::vector<bool> base_vector(L);
#pragma omp for reduction(+:value)
		for (int k = 0; k < N; k++) {
			int_to_binary(k, base_vector);
			base_vector[site] = !base_vector[site];
			const u64 idx = binary_to_int(base_vector);
			value += state_alfa(idx) * state_beta(k);
		}
	}
	return value;
}

/// <summary>
///
/// </summary>
/// <param name="n"></param>
/// <param name="m"></param>
/// <returns></returns>
double IsingModel_disorder::av_sigma_x_extensive(const u64 n, const u64 m) {
	std::vector<bool> base_vector(L), temp(L);
	cpx overlap = 0;
	cx_vec state_n = this->eigenvectors.col(n);
	cx_vec state_m = this->eigenvectors.col(m);
	for (long int k = 0; k < N; k++) {
		int_to_binary(k, base_vector);
		for (int j = 0; j < this->L; j++) {
			temp = base_vector;
			temp[j] = !base_vector[j];
			overlap += conj(state_n(binary_to_int(temp))) * state_m(k);
		}
	}
	return real(overlap) / double(this->L * this->L);
}

/// <summary>
/// Calculates the spin correlation matrix within a given state (non-equilibrium average)
/// </summary>
/// <param name="state_id"> index of given state </param>
/// <returns> correlation matrix </returns>
mat IsingModel_disorder::correlation_matrix(u64 state_id) {
	mat corr_mat(L, L, fill::zeros);
	const cx_vec state = this->eigenvectors.col(state_id);
#pragma omp parallel shared(state, corr_mat)
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

				S2_tmp = real((Szm * Szm + 0.5) * conj(state(p)) * state(p)); //  <Sz(m) Sz(m)> + 2*( <S(m)^+ S(m)^-> + <S(m)^- S(m)^+> ) on-site
				corr_mat(m, m) += S2_tmp;
				for (int k = m + 1; k < L; k++) {
					double Szk = 0;
					if (vect[k] == 1) Szk = 0.5;
					else Szk = -0.5;
					S2_tmp = Szm * Szk * real(conj(state(p)) * state(p)); //  <Sz(m) Sz(k)> for k > m

					// <S^+ S^->
					temp = vect;
					if (vect[m] == 1 && vect[k] == 0) {
						temp[m] = 0;
						temp[k] = 1;
						idx = binary_to_int(temp);
						S2_tmp += 0.5 * real(conj(state(idx)) * state(p));
					}
					//<S^- S^+>
					else if (vect[m] == 0 && vect[k] == 1) {
						temp[m] = 1;
						temp[k] = 0;
						idx = binary_to_int(temp);
						S2_tmp += 0.5 * real(conj(state(idx)) * state(p));
					}
					corr_mat(m, k) += S2_tmp;
					corr_mat(k, m) += S2_tmp;
				}
			}
		}
	}

	return corr_mat;
}

/// <summary>
/// Calculates the entropy of the system via the mixed density matrix
/// </summary>
/// <param name="state_id"> state index to produce the density matrix </param>
/// <param name="A_size"> size of subsystem </param>
/// <returns> entropy of considered systsem </returns>
double IsingModel_disorder::entaglement_entropy(u64 state_id, int A_size) {
	std::vector<bool> base_vector(L);
	const cx_vec state = this->eigenvectors.col(state_id);
	u64 dimA = std::pow(2, A_size);
	u64 dimB = std::pow(2, L - A_size);
	cx_mat rho(dimA, dimA, fill::zeros);
#pragma omp parallel for shared(rho,state, dimA, dimB)
	for (long int n = 0; n < N; n++) {
		if (abs(state(n)) < 1e-10) continue;
		u64 counter = 0;
		for (u64 m = n % dimB; m < N; m += dimB) {
			u64 idx = std::floor(1.0 * n / dimB);
			rho(idx, counter) += conj(state(n)) * state(m);
			counter++;
		}
	}
	vec probabilities;
	eig_sym(probabilities, rho);
	double entropy = 0;
#pragma omp parallel for reduction(+: entropy)
	for (int some_index = 0; some_index < dimA; some_index++) {
		auto value = probabilities(some_index);
		entropy += (abs(value) < 1e-12) ? 0 : -value * log(abs(value));
	}
	//entropy = -trace(rho * real(logmat(rho)));
	return entropy;
}