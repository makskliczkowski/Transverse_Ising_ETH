#include "include/IsingModel.h"

#include "include/IsingModel.h"

IsingModel_disorder::IsingModel_disorder(int L, vector<double>& J, double g, double h) {
    this->L = L; this->J = J; this->g = g; this->h = h;
    this->N = static_cast<u64>(std::pow(2, L));
    try {
        this->H = mat(N, N, fill::zeros); //hamiltonian
    }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }
    set_neighbours();
    hamiltonian();
}

IsingModel_disorder::IsingModel_disorder(const IsingModel_disorder& A) {
    this->L = A.L; this->J = A.J; this->g = A.g; this->h = A.h; this->N = A.N;
    this->H = A.H; this->eigenvectors = A.eigenvectors; this->eigenvalues = A.eigenvalues;
}
IsingModel_disorder::IsingModel_disorder(IsingModel_disorder&& A) noexcept {
    this->L = A.L; this->J = A.J; this->g = A.g; this->h = A.h; this->N = A.N;
    this->H = A.H; this->eigenvectors = A.eigenvectors; this->eigenvalues = A.eigenvalues;
}

void IsingModel_disorder::setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) {
    u64 idx = binary_to_int(temp);
    H(idx, k) += value;
    H(k, idx) += value;
}
void IsingModel_disorder::hamiltonian() {
    H.zeros();
    this->dh = w * create_random_vec(L);
    std::vector<bool> base_vector(L);
    std::vector<bool> temp(base_vector); // changes under H action
    for (u64 k = 0; k < N; k++) {
        int_to_binary(k, base_vector);
        int s_i, s_j;
        for (int j = 0; j <= L - 1; j++) {
            s_i = static_cast<double>(base_vector[j]);
            
            /* transverse field */
            temp = base_vector;
            temp[j] = ~base_vector[j];
            setHamiltonianElem(k, g, std::move(temp));

            /* disorder */
            H(k, k) += (this->h + dh(j)) * (s_i - 0.5);

            if (nearest_neighbours[j] >= 0) {
                /* Ising-like spin correlation */
                s_j = static_cast<double>(base_vector[nearest_neighbours[j]]);
                H(k, k) += this->J[j] * (s_i - 0.5) * (s_j - 0.5);
            }
        }
    }
}