#include "include/IsingModel.h"

/* CONSTRUCTORS */
IsingModel_disorder::IsingModel_disorder(int L, vector<double>& J, double g, double h, double disorder_strength) {
    this->L = L; this->J = J; this->g = g; this->h = h;
    this->disorder_strength = disorder_strength;
    this->N = static_cast<u64>(std::pow(2, L));
    set_neighbors();
    hamiltonian();

    /* SYMMETRY MATRICES */
    this->X = sp_mat(N, N);
    this->P = sp_mat(N, N);
    this->T = sp_mat(N, N);
}
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

/* BASE GENERATION AND RAPPING*/
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
/* SYMMETRY FUNCTIONS */

void IsingModel_disorder::create_X_matrix() {
    for (u64 k = 0; k < N; k++) {
        std::vector<bool> base_vector(L);
        int_to_binary(k, base_vector);
        base_vector.flip();
        u64 idx = binary_to_int(base_vector);
        X(idx, k) = 1;
    }
}



//-------------------------------------------------------------------------------------------------------------------------------
/* BUILDING HAMILTONIAN */
/// <summary>
/// Sets the non-diagonal elements of the Hamimltonian matrix, by acting with the operator on the k-th state
/// </summary>
/// <param name="k"> index of the basis state acted upon with the Hamiltonian </param>
/// <param name="value"> value of the given matrix element to be set </param>
/// <param name="temp"> resulting vector form acting with the Hamiltonian operator on the k-th basis state </param>
void IsingModel_disorder::setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) {
    u64 idx = binary_to_int(temp);
    H(idx, k) += value;
    //H(k, idx) += value;
}
/// <summary>
/// Generates the total Hamiltonian of the system. The diagonal part is straightforward, 
/// while the non-diagonal terms need the specialized setHamiltonainElem(...) function
/// </summary>
void IsingModel_disorder::hamiltonian() {
    try {
        this->H = mat(N, N, fill::zeros); //hamiltonian
    }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }
    this->dh = create_random_vec(L, this->disorder_strength);
    //out << dh.t();
    std::vector<bool> base_vector(L);
    std::vector<bool> temp(base_vector); // changes under H action
    for (u64 k = 0; k < N; k++) {
        int_to_binary(k, base_vector);
        double s_i, s_j;
        for (int j = 0; j <= L - 1; j++) {
            s_i = (base_vector[j] == 1) ? 1.0 : -1.0;
            
            /* transverse field */
            temp = base_vector;
            temp[j] = !base_vector[j];
            setHamiltonianElem(k, g, std::move(temp));

            /* disorder */
            H(k, k) += (this->h + dh(j)) * s_i;

            if (nearest_neighbors[j] >= 0) {
                /* Ising-like spin correlation */
                s_j = (base_vector[nearest_neighbors[j]] == 1) ? 1.0 : -1.0;
                H(k, k) += this->J[j] * s_i * s_j;
            }
        }
    }
}


//-------------------------------------------------------------------------------------------------------------------------------
/* PHYSICAL QUANTITES */
/// <summary>
/// Calculates the matrix element of the x-component spin matrix within the eigenstate state_id
/// </summary>
/// <param name="state_id"> index of eigenstate </param>
/// <param name="site"> position of the spin, where the operator is acted upon </param>
/// <returns> return the average in the eigenstate </returns>
double IsingModel_disorder::av_sigma_x(int state_id, int site) {
    vec state = eigenvectors.col(state_id); 
    std::vector<bool> base_vector(L), temp;
    double value = 0;
    for (int k = 0; k < N; k++) {
        int_to_binary(k, base_vector);
        temp = base_vector;
        temp[site] = !base_vector[site];
        u64 idx = binary_to_int(temp);
        value += state(idx) * state(k);
    }
    return value;
}

/// <summary>
/// Calculates the spin correlation matrix within a given state (non-equilibrium average)
/// </summary>
/// <param name="state_id"> index of given state </param>
/// <returns> correlation matrix </returns>
mat IsingModel_disorder::correlation_matrix(u64 state_id) {
    mat corr_mat(L, L, fill::zeros); 
    u64 idx;
    vector<bool> vect(L), temp(L);
    vec state = this->eigenvectors.col(state_id);
    for (int p = 0; p < N; p++) { //! 1 -localized, 2 - electrons
        int_to_binary(p, vect);
        for (int m = 0; m < L; m++) {
            double Szm = 0;
            double S2_tmp = 0;

            if (vect[m] == 1) Szm = 0.5;
            else Szm = -0.5;

            S2_tmp = Szm * Szm * state(p) * state(p) + 0.5 * state(p) * state(p); //  <Sz(m) Sz(m)> + 2*( <S(m)^+ S(m)^-> + <S(m)^- S(m)^+> ) on-site

            corr_mat(m, m) += S2_tmp;
            for (int k = m + 1; k < L; k++) {
                double Szk = 0;
                if (vect[k] == 1) Szk = 0.5;
                else Szk = -0.5;
                S2_tmp = Szm * Szk * (state(p) * state(p)); //  <Sz(m) Sz(k)> for k > m

                // <S^+ S^->
                temp = vect;
                if (vect[m] == 1 && vect[k] == 0) {
                    temp[m] = 0;
                    temp[k] = 1;
                    idx = binary_to_int(temp);
                    S2_tmp += 0.5 * state(idx) * state(p);
                }
                //<S^- S^+>
                temp = vect;
                if (vect[m] == 0 && vect[k] == 1) {
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
    vec state = this->eigenvectors.col(state_id);
    u64 dimA = std::pow(2, A_size);
    u64 dimB = std::pow(2, L - A_size);
    mat rho(dimA, dimA, fill::zeros);
    for (u64 n = 0; n < N; n++) {
        u64 counter = 0;
        for (u64 m = n % dimB; m < N; m += dimB) {
            u64 idx = std::floor(n / dimB);
            rho(idx, counter++) += state(n) * state(m);
        }
    }
    vec probabilities;
    eig_sym(probabilities, rho);
    double entropy = 0;
#pragma omp parallel for reduction(+: entropy)
    for (int some_index = 0; some_index < dimA; some_index++) {
        auto value = probabilities(some_index);
        entropy += (abs(value) < 1e-12)? 0 : -value * log(abs(value));
    }
    //entropy = -trace(rho * real(logmat(rho)));
    return entropy;
}
