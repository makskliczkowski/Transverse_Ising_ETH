#include "include/IsingModel.h"

/* GETTERS and SETTERS*/
u64 IsingModel::get_hilbert_size() {
    return N;
}
mat IsingModel::get_hamiltonian() {
    return H;
}
vec IsingModel::get_eigenvalues() {
    return eigenvalues;
}
mat IsingModel::get_eigenvectors() {
    return eigenvectors;
}

void IsingModel::set_neighbours() {
    this->nearest_neighbours = std::vector<int>(L);
    switch (_BC) {
    case 0:
        std::iota(nearest_neighbours.begin(), nearest_neighbours.end(), 1);
        nearest_neighbours[L - 1] = 0;
        break;

    case 1:
        std::iota(nearest_neighbours.begin(), nearest_neighbours.end(), 1);
        nearest_neighbours[L - 1] = -100;
        break;
    }
}

void IsingModel::diagonalization() {
    try {
        arma::eig_sym(eigenvalues, eigenvectors, H);
    }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        out << "dim(H) = " << H.size() * sizeof(H(0, 0)) << "\n";
        assert(false);
    }
}

double IsingModel::ipr(int state_idx) {
    double ipr = 0;
#pragma omp parallel for shared (ipr) reduction(+: ipr)
    for (int n = 0; n < N; n++) {
        ipr += std::pow(eigenvectors.col(state_idx)(n), 4);
    }
    return 1.0 / ipr;
}

