#include "include/IsingModel.h"

/* CONSTRUCTORS */
IsingModel_sym::IsingModel_sym(int L, double J, double g, double h, int k_sym, bool p_sym, bool x_sym, int _BC){
    this->L = L; this->J = J; this->g = g; this->h = h;
    this->info = "_L="+std::to_string(this->L) + \
        ",g=" + to_string_prec(this->g) + \
        ",h=" + to_string_prec(this->h);
    this->_BC = _BC;
    
    symmetries.k_sym = k_sym * two_pi / double(this->L);
    symmetries.p_sym = (p_sym) ? 1.0 : -1.0;
    symmetries.x_sym = (x_sym) ? 1.0 : -1.0;

    this->mapping = std::vector<u64>();
    this->normalisation = std::vector<cpx>();
    generate_mapping();
    /*for (int k = 0; k < N; k++) {
        std::vector<bool> temp(L);
        int_to_binary(mapping[k], temp);
        out << mapping[k] << "\t\t" << temp << "\t\t" << normalisation[k] << endl;
    }*/

    set_neighbors(); // generate neighbors
    hamiltonian();
}

/* BASE GENERATION, SEC CLASSES AND RAPPING*/
u64 IsingModel_sym::map(u64 index) {
    if (index >= this->N) throw "Element out of range\n No such index in map\n";
    return mapping[index];
}
/// <summary>
/// Finds the representative of the equivalent class after vector rotation (mainly used after applied another symmetry)
/// </summary>
/// <param name="base_vector"> vector from EC to find representative </param>
/// <returns> index of the representative state in the EC </returns>
std::tuple<u64, int> IsingModel_sym::find_translation_representative(std::vector<bool>& base_vector) const {
    u64 EC_symmetry = binary_to_int(base_vector);
    u64 current_idx = EC_symmetry;
    u64 idx = INT_MAX;
    int counter = 1, count_to_rep = 1;
    while(idx != current_idx) {
        std::rotate(base_vector.begin(), base_vector.begin() + 1, base_vector.end());
        idx = binary_to_int(base_vector);
        counter++;
        if (idx < EC_symmetry) {
            EC_symmetry = idx;
            count_to_rep = counter;
        }
    }
    return std::make_tuple(EC_symmetry, count_to_rep);
}
/// <summary>
/// Find representatives of other EC generated by reflection, spin-flip and (reflection x spin-flip) symmetry
/// </summary>
/// <param name="base_vector"> current base vector to act with symmetries </param>
/// <param name="min"> index of EC class representative by translation symmetry </param>
/// <returns></returns>
std::tuple<u64, int> IsingModel_sym::find_SEC_representative(const std::vector<bool>& base_vector) const {
    std::vector<u64> minima;
    std::vector<bool> temp = base_vector;
    bool k_sector = abs(symmetries.k_sym) < 1e-4 || abs(symmetries.k_sym - pi) < 1e-4;

    //check reflection symmetry
    std::reverse(temp.begin(), temp.end());
    auto tupleR = find_translation_representative(temp);
    minima.push_back((k_sector ? std::get<0>(tupleR) : INT_MAX));

    if (h == 0) {
        temp = base_vector;
        
        // check spin-flip
        temp.flip();
        auto tupleX = find_translation_representative(temp);
        minima.push_back(std::get<0>(tupleX));

        // check spin-flip and reflection
        std::reverse(temp.begin(), temp.end());
        auto tupleRX = find_translation_representative(temp);
        minima.push_back((k_sector ? std::get<0>(tupleRX) : INT_MAX));

        switch (std::min_element(minima.begin(), minima.end()) - minima.begin()) {
        case 1:
            return std::make_tuple(minima[1], symmetries.x_sym * std::get<1>(tupleX)); break;
        case 0:
            return std::make_tuple(minima[0], symmetries.p_sym * std::get<1>(tupleR)); break;
        case 2:
            return std::make_tuple(minima[2], symmetries.p_sym * symmetries.x_sym * std::get<1>(tupleRX)); break;
        default: throw "Index out of range\n";
        }
    }
    else
        return std::make_tuple(std::get<0>(tupleR), symmetries.p_sym * std::get<1>(tupleR));
}

/// <summary>
/// 
/// </summary>
/// <param name="base_vector"></param>
/// <param name="k"></param>
/// <returns></returns>
cpx IsingModel_sym::get_symmetry_normalization(std::vector<bool>& base_vector, u64 k) {
    cpx normalisation = cpx(0.0, 0.0);
    std::vector<bool> PT = base_vector;
    std::reverse(PT.begin(), PT.end());
    std::vector<bool> ZT = base_vector;
    ZT.flip();
    std::vector<bool> PZT = ZT;
    std::reverse(PZT.begin(), PZT.end());
    bool k_sector = abs(symmetries.k_sym) < 1e-4 || abs(symmetries.k_sym - pi) < 1e-4;
    for (int l = 0; l < L; l++) {
        if (binary_to_int(base_vector) == k) 
            normalisation += std::exp(-1i * symmetries.k_sym * double(l));
        if (k_sector && binary_to_int(PT) == k)
            normalisation += ((double)symmetries.p_sym) * std::exp(-1i * symmetries.k_sym * double(l));
        if ((h == 0) && binary_to_int(ZT) == k)
            normalisation += ((double)symmetries.x_sym) * std::exp(-1i * symmetries.k_sym * double(l));
        if (k_sector && (h == 0) && binary_to_int(PZT) == k)
            normalisation += ((double)symmetries.p_sym * (double)symmetries.x_sym) * std::exp(-1i * symmetries.k_sym * double(l));
        std::rotate(PT.begin(), PT.begin() + 1, PT.end());
        std::rotate(ZT.begin(), ZT.begin() + 1, ZT.end());
        std::rotate(PZT.begin(), PZT.begin() + 1, PZT.end());
        std::rotate(base_vector.begin(), base_vector.begin() + 1, base_vector.end());
    }
    return std::sqrt(normalisation);
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
void IsingModel_sym::mapping_kernel(u64 start, u64 stop, std::vector<u64>& map_threaded, std::vector<cpx>& norm_threaded, int _id) {
    std::vector<bool> base_vector(L); // temporary dirac-notation base vector
    for (u64 j = start; j < stop; j++) {
        int_to_binary(j, base_vector);
        if (g == 0 && std::accumulate(base_vector.begin(), base_vector.end(), 0) != L / 2.) continue;
        //check translation
        auto [min, trans_eig] = find_translation_representative(base_vector);

        u64 min_R_RX = INT_MAX;
        if (min == j) {
            auto tuple = find_SEC_representative(base_vector);
            min_R_RX = std::get<0>(tuple);
        }
        if (min_R_RX < j) continue;
        
        if (std::min(min, min_R_RX) == j) {
            cpx N = get_symmetry_normalization(base_vector, j);             // normalisation condition -- check wether state in basis
            if (std::abs(N) > 1e-6) {
                //out << base_vector << "\t\t" << N << endl;
                map_threaded.push_back(j);
                norm_threaded.push_back(N);
            }
        }
    }
}
/// <summary>
/// Splits the mapping onto threads, where each finds basis states in the reduced Hilbert space within a given range.
/// The mapping is retrieved by concatenating the resulting maps from each thread
/// </summary>
void IsingModel_sym::generate_mapping() {
    u64 start = 0, stop = static_cast<u64>(std::pow(2, L)); // because parity reflects the other half
    if (num_of_threads == 1)
        mapping_kernel(start, stop, mapping, normalisation, 0);
    else {
        //Threaded
        //std::vector<my_uniq_ptr> map_threaded(num_of_threads);
        v_2d<u64> map_threaded(num_of_threads);
        v_2d<cpx> norm_threaded(num_of_threads);
        std::vector<std::thread> threads;
        threads.reserve(num_of_threads);
        for (int t = 0; t < num_of_threads; t++) {
            start = (u64)(std::pow(2, L) / (double)num_of_threads * t);
            stop = ((t + 1) == num_of_threads ? (u64)std::pow(2, L) : u64(std::pow(2, L) / (double)num_of_threads * (double)(t + 1)));
            map_threaded[t] = v_1d<u64>();
            norm_threaded[t] = v_1d<cpx>();
            threads.emplace_back(&IsingModel_sym::mapping_kernel, this, start, stop, ref(map_threaded[t]), ref(norm_threaded[t]), t);
        }
        for (auto& t : threads) t.join();

        for (auto& t : map_threaded)
            mapping.insert(mapping.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));

        for (auto& t : norm_threaded)
            normalisation.insert(normalisation.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));
    }
    this->N = this->mapping.size();
    //assert(mapping.size() > 0 && "Not possible number of electrons - no. of states < 1");
}


/* BUILDING HAMILTONIAN */
/// <summary>
/// Sets the non-diagonal elements of the Hamimltonian matrix with symmetry sectors: therefore the matrix elements are summed over the SEC
/// </summary>
/// <param name="k"> index of the basis state acted upon with the Hamiltonian </param>
/// <param name="value"> value of the given matrix element to be set </param>
/// <param name="temp"> resulting vector form acting with the Hamiltonian operator on the k-th basis state </param>
void IsingModel_sym::setHamiltonianElem(u64 k, double value, std::vector<bool>& temp) {
    auto idx_temp = binary_to_int(temp);
    u64 idx = binary_search(mapping, 0, N - 1, idx_temp);
    int sym_eig = 1;
    if (idx == -1) {
        auto tup_T = find_translation_representative(temp);
        auto tup_S = find_SEC_representative(temp);
        auto [min, trans_eig] = (std::get<0>(tup_T) > std::get<0>(tup_S)) ? tup_S : tup_T;
        sym_eig = trans_eig;
        //finding index in reduced Hilbert space
        idx = binary_search(mapping, 0, N - 1, min);
        if (idx < 0 || idx >= N) return;// out << "Element do not exist\n";
    }
    cpx translation_eig = (abs(sym_eig) == 1 || symmetries.k_sym == 0) ? \
        cpx(1.0) : std::exp(-1i * symmetries.k_sym * double(abs(sym_eig) - 1));
    cpx value_new = value * translation_eig * (normalisation[idx] / normalisation[k]) * double(sgn(sym_eig));
    H(idx, k) += value_new;
}
/// <summary>
/// Generates the total Hamiltonian of the system. The diagonal part is straightforward, 
/// while the non-diagonal terms need the specialized setHamiltonainElem(...) function
/// </summary>
void IsingModel_sym::hamiltonian() {
    try {
        this->H = cx_mat(N, N, fill::zeros); //hamiltonian
    }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }
    std::vector<bool> base_vector(L);
    std::vector<bool> temp(base_vector); // changes under H action
    for (long int k = 0; k < N; k++) {
        int_to_binary(mapping[k], base_vector);
        double s_i, s_j;
        for (int j = 0; j <= L - 1; j++) {
            s_i = base_vector[j] ? 1.0 : -1.0;                              // true - spin up, false - spin down
            /* transverse field */
            if (g != 0) {
                temp = base_vector;
                temp[j] = !base_vector[j];                                  // negates on that site
                setHamiltonianElem(k, this->g, temp);
            }
            /* disorder */
            H(k, k) += this->h * s_i;                                       // diagonal

            if (nearest_neighbors[j] >= 0) {
                /* Ising-like spin correlation */
                s_j = base_vector[nearest_neighbors[j]] ? 1.0 : -1.0;
                this->H(k, k) += this->J * s_i * s_j;
            }
        }
    }
}

/* PHYSICAL QUANTITTIES */
/// <summary>
/// 
/// </summary>
/// <param name="n"></param>
/// <param name="m"></param>
/// <returns></returns>
double IsingModel_sym::spin_flip_mat_element(const u64 n, const u64 m) {
    std::vector<bool> base_vector(L), temp(L);
    cpx overlap = 0;
    cx_vec state_n = this->eigenvectors.col(n);
    cx_vec state_m = this->eigenvectors.col(m);
    for (long int k = 0; k < N; k++) {
        int_to_binary(mapping[k], base_vector);
        for (int j = 0; j < this->L; j++) {
            temp = base_vector;
            temp[j] = !base_vector[j];
            auto idx_temp = binary_to_int(temp);
            u64 idx = binary_search(mapping, 0, N - 1, idx_temp);
            int sym_eig = 1;
            if (idx == -1) {
                auto tup_T = find_translation_representative(temp);
                auto tup_S = find_SEC_representative(temp);
                auto [min, trans_eig] = (std::get<0>(tup_T) > std::get<0>(tup_S)) ? tup_S : tup_T;
                sym_eig = trans_eig;
                idx = binary_search(mapping, 0, N - 1, min);
                if (idx < 0 || idx >= N) continue; // out << "Element do not exist\n";
            }
            cpx translation_eig = (abs(sym_eig) == 1 || symmetries.k_sym == 0) ? \
                cpx(1.0) : std::exp(-1i * symmetries.k_sym * double(abs(sym_eig) - 1));
            cpx value_new = translation_eig * (normalisation[idx] / normalisation[k]) * double(sgn(sym_eig));
            overlap += conj(state_n(idx)) * value_new * state_m(k);
        }
    }
    return real(overlap);// / double(this->L * this->L);
}