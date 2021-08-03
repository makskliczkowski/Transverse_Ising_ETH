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

/* HELPER FUNCTIONS */

/// <summary>
/// 
/// </summary>
/// <param name="Sz"></param>
void IsingModel::print_base_spin_sector(int Sz) {
    std::vector<bool> temp(L);
    for (int k = 0; k < N; k++) {
        int_to_binary((use_mapping)? mapping[k] : k, temp);
        int Sz = 0;
        for (int l = 0; l < L; l++)
            Sz += temp[l] == 0 ? -1 : 1;
        if (Sz == 0)
            out << temp << "\t\t" << ((use_mapping)? periodicity[k] : 1) << std::endl;
    }
}

/// <summary>
/// Sets the neigbours depending on the Boundary condition (BC) defined as a makro in the 'headers.h' file
/// </summary>
void IsingModel::set_neighbours() {
    this->nearest_neighbours = std::vector<int>(L);
    switch (_BC) {
    case 0:
        std::iota(nearest_neighbours.begin(), nearest_neighbours.end(), 1);
        nearest_neighbours[L - 1] = 0;
        break;

    case 1:
        std::iota(nearest_neighbours.begin(), nearest_neighbours.end(), 1);
        nearest_neighbours[L - 1] = -1;
        break;
    }
}

/// <summary>
/// General procedure to diagonalize the Hamiltonian using eig_sym from the Armadillo library
/// </summary>
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

/// <summary>
/// The IPR, also called the participation ratio, quantifies how delocalized a state is in a certain basis. 
/// The state is completely delocalized, when:    IPR=dim(hilbert space)
/// </summary>
/// <param name="state_idx"> index of the eigenvector used to calculate this quantity </param>
/// <returns> returns the IPR value</returns>
double IsingModel::ipr(int state_idx) {
    double ipr = 0;
#pragma omp parallel for shared (ipr) reduction(+: ipr)
    for (int n = 0; n < N; n++) {
        ipr += std::pow(eigenvectors.col(state_idx)(n), 4);
    }
    return 1.0 / ipr;
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
    double rn_next = 0, rn = (A.*op)(0, 1);
    double average = 0;
    for (int k = 1; k < A.N; k++) {
        rn_next = (A.*op)(k, 1);
        average += abs(rn_next - rn);
        rn = rn_next;
    }
    return average / (A.N - 1.0);
}

/// <summary>
/// Prints to file the average of a given operator in each eigenstate as a function of eigenenergies
/// </summary>
/// <typeparam name="T"> typename as class: model with disorder or symmetries </typeparam>
/// <param name="op"> operator as function acting on a specific eigenstate </param>
/// <param name="A"> class instance </param>
/// <param name="site"> position of the spin, where the operator is acted upon </param>
/// <param name="name"> name of the file to print data </param>
/// <param name="separator"> separator between columns in file </param>
void IsingModel::operator_av_in_eigenstates(double (IsingModel::* op)(int, int),\
    IsingModel& A, int site, std::string name, std::string separator) {
    std::ofstream file(name);
    if (!file.is_open()) {
        throw "Can't open file " + name + "\n Choose another file\n";
    }
    for (int k = 0; k < A.N; k++)
        file << A.eigenvalues(k) / (double)A.L << separator << (A.*op)(k, site) << endl;
    file.close();
}
