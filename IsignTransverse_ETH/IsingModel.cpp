#include "include/IsingModel.h"

IsingModel::~IsingModel(){
}
/* GETTERS and SETTERS*/
std::string IsingModel::get_info() const{
    return this->info;
}
u64 IsingModel::get_hilbert_size() const {
    return this->N;
}
const mat& IsingModel::get_hamiltonian() const {
    return this->H;
}
const vec& IsingModel::get_eigenvalues() const {
    return this->eigenvalues;
}
const mat& IsingModel::get_eigenvectors() const {
    return this->eigenvectors;
}
const std::vector<u64>& IsingModel::get_mapping() const{
    return this->mapping;
}
/// <summary>
/// Return given eigenenergy at index idx
/// </summary>
/// <param name="idx"> index of eigenvalue </param>
/// <returns> eigenvalue at idx </returns>
double IsingModel::get_eigenEnergy(int idx) const
{
    return this->eigenvalues(idx);
}
/// <summary>
/// Returns given eigenstate at index idx
/// </summary>
/// <param name="idx"> index of eigenstate </param>
/// <returns> eigenstate at idx </returns>
const vec& IsingModel::get_eigenState(int idx) const {
    return this->eigenvectors.col(idx);
}

// ---- PRINTERS
/// <summary>
/// prints the basis vector in the given spin-symmetry block
/// </summary>
/// <param name="Sz"> custom spin symmetry block </param>
void IsingModel::print_base_spin_sector(int Sz) {
    std::vector<bool> temp(L);
    for (int k = 0; k < N; k++) {
        int_to_binary(map(k), temp);
        const int sum = std::accumulate(temp.begin(), temp.end(), 0);
            //,[](int a, bool b)
            //{return a + ((b) ? 1:-1);});
        //for (int l = 0; l < L; l++)
        //    Sz += temp[l] == 0 ? -1 : 1;
        if(sum == Sz)
            out << temp << std::endl;
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
/// calculates the commutator of two operatos (one of the symmetry operators or hamiltonian) defined in enum: oeprators
/// </summary>
/// <param name="A"> first operator </param>
/// <param name="B"> second operator</param>
/// <returns> false if commute and true otherwise </returns>
bool IsingModel::commutator(IsingModel ::operators A, IsingModel::operators B) {
    sp_mat A1 = choose_operator(A);
    sp_mat A2 = choose_operator(B);
    if (A1.is_zero(1e-12)) throw "First operator matrix is empty:\n either not created or wrong parameters\n";
    if (A2.is_zero(1e-12)) throw "Second operator matrix is empty:\n either not created or wrong parameters\n";
    return (A1 * A2 - A2 * A1).is_zero(1e-10);
}

/// <summary>
/// Choosing the right operator
/// </summary>
/// <param name="A">the option</param>
/// <returns>sparse operator</returns>
sp_mat IsingModel::choose_operator(IsingModel::operators A) {
    switch (A) {
    case IsingModel::operators::H:
         return sp_mat(this->H); break;
    case IsingModel::operators::X:
        return this->X; break;
    case IsingModel::operators::P:
        return this->P; break;
    case IsingModel::operators::T:
        return this->T; break;
    default:
        out << "Wrong input value:\nSetting by default unit matrix!!!" << endl;
        return (sp_mat)eye(N, N);
    }
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
/// prints the state in the regular basis (only the highest coefficients are present)
/// </summary>
/// <param name="_id"> index of the printed state </param>
void IsingModel::print_state(u64 _id) {
    vec state = eigenvectors.col(_id);
    double max = arma::max(state);
    std::vector<bool> base_vector(L);
    for (int k = 0; k < N; k++) {
        int_to_binary(map(k), base_vector);
        if (abs(state(k)) >= 0.01 * max) {
            out << state(k) << " * |" << base_vector << "> + ";
        }
    }
    out << endl;
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
    //double delta_n = eigenvalues(_min) - eigenvalues(_min - 1);
    //double delta_n_next = 0;
#pragma omp parallel for reduction(+: r)
    for (int k = _min; k < _max - 1; k++) {
        const double delta_n = eigenvalues(k) - eigenvalues(k-1);
        const double delta_n_next = eigenvalues(k+1) - eigenvalues(k);
        const double min = std::min(delta_n, delta_n_next);
        const double max = std::max(delta_n, delta_n_next);
        if (max == 0) throw "Degeneracy!!!\n";
        r += min / max;
        //delta_n = delta_n_next;
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
    //double rn_next = 0, rn = (A.*op)(0, 1);
    double average = 0;
#pragma omp parallel for reduction(+:average)
    for (int k = 1; k < A.N; k++) {
        const double rn = (A.*op)(k-1,1);
        const double rn_next = (A.*op)(k, 1);
        average += abs(rn_next - rn);
        //rn = rn_next;
    }
    return average / (A.N - 1.0);
}

/// <summary>
/// Prints to file the average of a given operator in each eigenstate as a function of eigenenergies
/// </summary>
/// <param name="op"> operator as function acting on a specific eigenstate </param>
/// <param name="A"> class instance </param>
/// <param name="site"> position of the spin, where the operator is acted upon </param>
/// <param name="name"> name of the file to print data </param>
/// <param name="separator"> separator between columns in file </param>
void IsingModel::operator_av_in_eigenstates(double (IsingModel::* op)(int, int),\
    IsingModel& A, int site, std::string name, std::string separator) {
    std::ofstream file(name);                                                                  // file to write the average to
    if (!file.is_open()) throw "Can't open file " + name + "\n Choose another file\n";
    // vec res(A.get_hilbert_size(), fill::zeros);

#pragma omp parallel for
    for (int k = 0; k < A.get_hilbert_size(); k++){
        const double res = (A.*op)(k, site);
#pragma omp critical
        file << A.eigenvalues(k) / (double)A.L << separator << res << endl;
    }
    //for (int k = 0; k < A.get_hilbert_size(); k++)
    //    file << A.eigenvalues(k) / (double)A.L << separator << res(k) << endl;
    file.close();
}

/// <summary>
/// Prints to vector the average of a given operator in each eigenstate as a function of eigenenergies
/// </summary>
/// <param name="op"> operator as function acting on a specific eigenstate </param>
/// <param name="A"> class instance </param>
/// <param name="site"> position of the spin, where the operator is acted upon </param>
/// <param name="separator"> separator between columns in file </param>
vec IsingModel::operator_av_in_eigenstates_return(double (IsingModel::* op)(int, int), IsingModel& A, int site) {
    vec temp(A.get_hilbert_size(),fill::zeros);
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
double quantum_fidelity(u64 _min, u64 _max, const std::unique_ptr<IsingModel>& Hamil, vector<double>& J, double g, double h, double w){
    if (_min < 0) throw "too low index";
    if (_max >= Hamil->get_hilbert_size()) throw "index exceeding Hilbert space";

    std::unique_ptr<IsingModel> Hamil2;
    if (typeid(Hamil) == typeid(IsingModel_disorder)) Hamil2.reset(new IsingModel_disorder(Hamil->L, J, g, h, w));
    else Hamil2.reset(new IsingModel_sym(Hamil->L, J, g, h));


    double fidelity = 0;
#pragma omp parallel for reduction(+: fidelity)
    for (long int k = _min; k < _max; k++)
        fidelity += overlap(Hamil, Hamil2, k, k);
    return fidelity / double(_max - _min);
}
