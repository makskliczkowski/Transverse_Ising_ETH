#include "include/IsingModel.h"

/* CONSTRUCTORS */
IsingModel_sym::IsingModel_sym(int L, vector<double>& J, double g, double h){
    this->L = L; this->J = J; this->g = g; this->h = h;

    this->mapping = std::vector<u64>();
    generate_mapping();
    this->periodicity = std::vector<int>(N);
    check_periodicity(); // generate the periodicity of the basis states
    //print_base_spin_sector(0);
    set_neighbors(); // generate neighbors
    hamiltonian();
}
IsingModel_sym::IsingModel_sym(const IsingModel_sym& A) {
    this->L = A.L; this->J = A.J; this->g = A.g; this->h = A.h;
    this->N = A.N; this->mapping = A.mapping;
    this->H = A.H; this->eigenvectors = A.eigenvectors; this->eigenvalues = A.eigenvalues;
}
IsingModel_sym::IsingModel_sym(IsingModel_sym&& A) noexcept {
    this->L = A.L; this->J = A.J; this->g = A.g; this->h = A.h;
    this->N = A.N; this->mapping = A.mapping;
    this->H = A.H; this->eigenvectors = A.eigenvectors; this->eigenvalues = A.eigenvalues;
}
IsingModel_sym::~IsingModel_sym()
{
	//out << "Destroying the Ising model with symmetries\n";
}

/* BASE GENERATION, SEC CLASSES AND RAPPING*/
u64 IsingModel_sym::map(u64 index) {
    if (index < 0 || index >= this->N) throw "Element out of range\n No such index in map\n";
    return mapping[index];
}
/// <summary>
/// Finds the representative of the equivalent class after vector rotation (mainly used after applied another symmetry)
/// </summary>
/// <param name="base_vector"> vector from EC to find representative </param>
/// <returns> index of the representative state in the EC </returns>
u64 IsingModel_sym::find_translation_representative(std::vector<bool>& base_vector) {
    u64 EC_symmetry = binary_to_int(base_vector), current_idx = EC_symmetry;
    u64 idx = INT_MAX;
    do {
        std::rotate(base_vector.begin(), base_vector.begin() + 1, base_vector.end());
        idx = binary_to_int(base_vector);
        if (idx < EC_symmetry)
            EC_symmetry = idx;
    } while (idx != current_idx);
    return EC_symmetry;
}
/// <summary>
/// Find representatives of other EC generated by reflection, spin-flip and (reflection x spin-flip) symmetry
/// </summary>
/// <param name="base_vector"> current base vector to act with symmetries </param>
/// <param name="min"> index of EC class representative by translation symmetry </param>
/// <returns></returns>
std::vector<u64> IsingModel_sym::find_SEC_representative(std::vector<bool>& base_vector) {
    std::vector<u64> minima;
    std::vector<bool> temp = base_vector;

    //check reflection symmetry
    std::reverse(temp.begin(), temp.end());
    minima.push_back(find_translation_representative(temp));
    //minima.push_back(INT_MAX);
    if (h == 0) {
        temp = base_vector;
        
        // check spin-flip and reflection
        temp.flip();
        minima.push_back(find_translation_representative(temp));

        // check spin-flip
        std::reverse(temp.begin(), temp.end());
        minima.push_back(find_translation_representative(temp));

    }
    return minima;
}

/// <summary>
/// Generates the mapping to the reduced Hilbert space (reduced by symmetries: translation, spin-flip and reflection symmetry
/// adding Sz=0 (total spin) symmetry is straightforward, however, the transverse field breaks the SU(2) symmetry;
/// 
/// The procedure hase been successfully optimized using multithreading:
/// - each thread functions in the range [start, stop)
/// </summary>
/// <param name="start"> first index for a given thread from the original Hilbert space </param>
/// <param name="stop"> last index for a given thread from the original Hilbert space </param>
/// <param name="map_threaded"> vector containing the mapping from the reduced basis to the original Hilbert space
///                             for a given thread, the whole mapping will be merged in the generate_mapping() procedure </param>
/// <param name="_id"> identificator for a given thread </param>
void IsingModel_sym::mapping_kernel(u64 start, u64 stop, std::vector<u64>& map_threaded, int _id) {
    std::vector<bool> base_vector(L); // temporary dirac-notation base vector
    for (u64 j = start; j < stop; j++) {
        int_to_binary(j, base_vector);

        //check translation
        u64 min = find_translation_representative(base_vector);

        u64 min_R_RX = INT_MAX;
        if (min == j) {
            auto minima = find_SEC_representative(base_vector );
            min_R_RX = *std::min_element(minima.begin(), minima.end());
            if (min_R_RX < j) continue;
        }

        if (std::min(min, min_R_RX) == j)
            map_threaded.push_back(j);
    }
}
/// <summary>
/// Splits the mapping onto threads, where each finds basis states in the reduced Hilbert space within a given range.
/// The mapping is retrieved by concatenating the resulting maps from each thread
/// </summary>
void IsingModel_sym::generate_mapping() {
    u64 start = 0, stop = static_cast<u64>(std::pow(2, L)); // because parity reflects the other half
    if (num_of_threads == 1)
        mapping_kernel(start, stop, mapping, 0);
    else {
        //Threaded
        //std::vector<my_uniq_ptr> map_threaded(num_of_threads);
        std::vector<std::vector<u64>> map_threaded(num_of_threads);
        std::vector<std::thread> threads;
        threads.reserve(num_of_threads);
        for (int t = 0; t < num_of_threads; t++) {
            start = (u64)(std::pow(2, L) / (double)num_of_threads * t);
            stop = ((t + 1) == num_of_threads ? (u64)std::pow(2, L) : u64(std::pow(2, L) / (double)num_of_threads * (double)(t + 1)));
            map_threaded[t] = std::vector<u64>();
            threads.emplace_back(&IsingModel_sym::mapping_kernel, this, start, stop, ref(map_threaded[t]), t);
        }
        for (auto& t : threads) t.join();

        for (auto& t : map_threaded)
            mapping.insert(mapping.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));
        mapping.shrink_to_fit();
        sort(mapping.begin(), mapping.end());
    }
    this->N = this->mapping.size();
    assert(mapping.size() > 0 && "Not possible number of electrons - no. of states < 1");
}

/// <summary>
/// Iterates over the reduced Hilbert space and checks the periodicity under translation 
/// of each basis state and the multipliocity of the SEC states
/// </summary>
void IsingModel_sym::check_periodicity() {
#pragma omp parallel for
    for (int k = 0; k < N; k++) {
        std::vector<bool> base_vector(L), temp(L);
        int_to_binary(mapping[k], base_vector);
        temp = base_vector;
        u64 idx = 0;
        int period_EC = 0; // because it always eneter do-while loop
        do {
            std::rotate(temp.begin(), temp.begin() + 1, temp.end());
            idx = binary_to_int(temp);
            period_EC++;
        } while (idx != mapping[k]); 

        auto SEC = find_SEC_representative(base_vector);
        SEC.push_back(mapping[k]);
        std::sort(SEC.begin(), SEC.end());
        SEC.erase(unique(SEC.begin(), SEC.end()), SEC.end());
        int multiplicity_in_SEC = SEC.size();
        //if (h == 0) multiplicity_in_SEC = SEC.size();
        //else multiplicity_in_SEC = 1;
        this->periodicity[k] = period_EC * multiplicity_in_SEC;
    }
}


/* BUILDING HAMILTONIAN */
/// <summary>
/// Sets the non-diagonal elements of the Hamimltonian matrix with symmetry sectors: therefore the matrix elements are summed over the SEC
/// </summary>
/// <param name="k"> index of the basis state acted upon with the Hamiltonian </param>
/// <param name="value"> value of the given matrix element to be set </param>
/// <param name="temp"> resulting vector form acting with the Hamiltonian operator on the k-th basis state </param>
void IsingModel_sym::setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) {

    u64 min = find_translation_representative(temp);

    auto minima = find_SEC_representative(temp);

    //finding index in reduced Hilbert space
    u64 idx = binary_search(mapping, 0, N - 1,\
        std::min(min, *std::min_element(minima.begin(), minima.end())));

    value = value * sqrt((double)periodicity[k] / (double)periodicity[idx]);
    H(idx, k) += value;
    //H(k, idx) += value;
}
/// <summary>
/// Generates the total Hamiltonian of the system. The diagonal part is straightforward, 
/// while the non-diagonal terms need the specialized setHamiltonainElem(...) function
/// </summary>
void IsingModel_sym::hamiltonian() {
    try {
        this->H = mat(N, N, fill::zeros); //hamiltonian
    }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }    
    std::vector<bool> base_vector(L);
    std::vector<bool> temp(base_vector); // changes under H action
    for (u64 k = 0; k < N; k++) {
        int_to_binary(mapping[k], base_vector);
        double s_i, s_j;
        for (int j = 0; j <= L - 1; j++) {
            s_i = base_vector[j];

            /* transverse field */
            temp = base_vector;
            temp[j] = !base_vector[j];
            setHamiltonianElem(k, 0.5 * g, std::move(temp));

            /* disorder */
            H(k, k) += this->h * (s_i - 0.5);

            if (nearest_neighbors[j] >= 0) {
                /* Ising-like spin correlation */
                s_j = base_vector[nearest_neighbors[j]];
                H(k, k) += this->J[j] * (s_i - 0.5) * (s_j - 0.5);
            }
        }
    }
}

/* PHYSICAL QUANTITTIES */
/// <summary>
/// Calculates the matrix element of the x-component spin matrix within the eigenstate state_id
/// </summary>
/// <param name="state_id"> index of eigenstate </param>
/// <param name="site"> position of the spin, where the operator is acted upon </param>
/// <returns> return the average in the eigenstate </returns>
double IsingModel_sym::av_sigma_x(int state_id, int site) {
    vec state = eigenvectors.col(state_id);
    double value = 0;
#pragma omp parallel for reduction(+: value)
    for (int k = 0; k < N; k++) {
        std::vector<bool> base_vector(L), temp;
        int_to_binary(k, base_vector);
        temp = base_vector;
        temp[site] = !base_vector[site];

        // FIND SEC represetnative!!!
        u64 min = find_translation_representative(temp);

        auto minima = find_SEC_representative(temp);

        //finding index in reduced Hilbert space
        u64 idx = binary_search(mapping, 0, N - 1, \
            std::min(min, *std::min_element(minima.begin(), minima.end())));
        //u64 idx = binary_to_int(temp);

        value += 0.5 * state(idx) * state(k) * sqrt((double)periodicity[k] / (double)periodicity[idx]);
    }
    return value;
}