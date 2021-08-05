#include "include/IsingModel.h"

IsingModel::~IsingModel(){
}
/* GETTERS and SETTERS*/
u64 IsingModel::get_hilbert_size() const {
    return N;
}
mat IsingModel::get_hamiltonian() const {
    return H;
}
vec IsingModel::get_eigenvalues() const {
    return eigenvalues;
}
mat IsingModel::get_eigenvectors() const {
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
void IsingModel::set_neighbors() {
    this->nearest_neighbors = std::vector<int>(L);
    switch (_BC) {
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
#pragma omp parallel for reduction(+: ipr)
    for (int n = 0; n < N; n++) {
        ipr += std::pow(eigenvectors.col(state_idx)(n), 4);
    }
    return 1.0 / ipr;
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
    double delta_n = eigenvalues(_min) - eigenvalues(_min - 1);
    double delta_n_next = 0;
    for (int k = _min + 1; k < _max; k++) {
        delta_n_next = eigenvalues(k) - eigenvalues(k - 1);
        double min = std::min(delta_n, delta_n_next), max = std::max(delta_n, delta_n_next);
        if (max == 0) throw "Degeneracy!!!\n";
        r += min / max;
        delta_n = delta_n_next;
    }
    return r / double(_max - _min);
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
