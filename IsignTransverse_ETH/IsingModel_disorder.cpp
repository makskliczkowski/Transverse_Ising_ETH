#include "include/IsingModel.h"

/* CONSTRUCTORS */
IsingModel_disorder::IsingModel_disorder(int L, vector<double>& J, double g, double h) {
    use_mapping = false;
    this->L = L; this->J = J; this->g = g; this->h = h;
    this->N = static_cast<u64>(std::pow(2, L));
    try {
        this->H = mat(N, N, fill::zeros); //hamiltonian
    }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }
    set_neighbors();
    //print_base_spin_sector(0);
    hamiltonian();
    //out << H;
}
IsingModel_disorder::IsingModel_disorder(const IsingModel_disorder& A) {
    this->L = A.L; this->J = A.J; this->g = A.g; this->h = A.h; this->N = A.N;
    this->H = A.H; this->eigenvectors = A.eigenvectors; this->eigenvalues = A.eigenvalues;
}
IsingModel_disorder::IsingModel_disorder(IsingModel_disorder&& A) noexcept {
    this->L = A.L; this->J = A.J; this->g = A.g; this->h = A.h; this->N = A.N;
    this->H = A.H; this->eigenvectors = A.eigenvectors; this->eigenvalues = A.eigenvalues;
}
IsingModel_disorder::~IsingModel_disorder()
{
	//out << "Destroying the disordered Ising model\n";
}
/* HELPER FUNCTIONS */

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
    H.zeros();
    this->dh = disorder_strength * create_random_vec(L);
    //out << dh.t();
    std::vector<bool> base_vector(L);
    std::vector<bool> temp(base_vector); // changes under H action
    for (u64 k = 0; k < N; k++) {
        int_to_binary(k, base_vector);
        int s_i, s_j;
        for (int j = 0; j <= L - 1; j++) {
            s_i = base_vector[j];
            
            /* transverse field */
            temp = base_vector;
            temp[j] = !base_vector[j];
            setHamiltonianElem(k, g, std::move(temp));

            /* disorder */
            H(k, k) += (this->h + dh(j)) * (s_i - 0.5);

            if (nearest_neighbors[j] >= 0) {
                /* Ising-like spin correlation */
                s_j = base_vector[nearest_neighbors[j]];
                H(k, k) += this->J[j] * (s_i - 0.5) * (s_j - 0.5);
            }
        }
    }
}

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





