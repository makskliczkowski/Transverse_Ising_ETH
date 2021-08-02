#include "include/IsingModel.h"

IsingModel_sym::IsingModel_sym(int L, vector<double>& J, double g, double h){
	this->L = L; this->J = J; this->g = g; this->h = h;

    this->mapping = std::vector<u64>();
    generate_mapping();
    this->N = this->mapping.size();
    this->periodicity = std::vector<int>(N);
    check_periodicity();

    try {
        this->H = mat(N, N, fill::zeros); //hamiltonian
    }
    catch (const bad_alloc& e) {
        std::cout << "Memory exceeded" << e.what() << "\n";
        assert(false);
    }
    std::vector<bool> temp(L);
    for (int k = 0; k < N; k++) {
        int_to_binary(mapping[k], temp);
        int Sz = 0;
        for (int l = 0; l < L; l++) 
            Sz += temp[l] == 0 ? -1 : 1;
        if(Sz == 0)
            out << temp << endl;
    }
    set_neighbours();
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

/// <summary>
/// 
/// </summary>
/// <param name="start"></param>
/// <param name="stop"></param>
/// <param name="map_threaded"></param>
/// <param name="_id"></param>
void IsingModel_sym::mapping_kernel(u64 start, u64 stop, std::vector<u64>& map_threaded, int _id) {
    std::vector<bool> base_vector(L), temp(L); // temporary dirac-notation base vector
    u64 idx = 0, min = INT_MAX;

    for (u64 j = start; j < stop; j++) {
        int_to_binary(j, base_vector);
        
        //check transaltion
        min = j;
        do {
            std::rotate(base_vector.begin(), base_vector.begin() + 1, base_vector.end());
            idx = binary_to_int(base_vector);
            if (idx <= min)
                min = idx;
        } while (idx != j);
        int_to_binary(j, base_vector);

        u64 min_R_RX = INT_MAX;
        if (min == j) {
            temp = base_vector;
            std::reverse(temp.begin(), temp.end());
            u64 minR = binary_to_int(temp), idx_R = minR;
            if (minR < min) continue;
            do {
                std::rotate(temp.begin(), temp.begin() + 1, temp.end());
                idx = binary_to_int(temp);
                if (idx <= minR)
                    minR = idx;
            } while (idx != idx_R);
            if (minR < min) continue;

            if (h == 0) {
                temp = base_vector;
                // check flip and reverse
                temp.flip();
                u64 minRX = binary_to_int(temp), idx_RX = minRX;
                if (minRX < min) continue;
                do {
                    std::rotate(temp.begin(), temp.begin() + 1, temp.end());
                    idx = binary_to_int(temp);
                    if (idx <= minRX)
                        minRX = idx;
                } while (idx != idx_RX);
                if (minRX < min) continue;

                // check flip
                std::reverse(temp.begin(), temp.end());
                u64 minX = binary_to_int(temp), idx_X = minX;
                if (minX < min) continue;
                do {
                    std::rotate(temp.begin(), temp.begin() + 1, temp.end());
                    idx = binary_to_int(temp);
                    if (idx <= minX)
                        minX = idx;
                } while (idx != idx_X);

                min_R_RX = std::min(minX, minRX);
            }
            else
                min_R_RX = j;
            min_R_RX = std::min(minR, min_R_RX);
        }
            
        if (std::min(min, min_R_RX) == j)
            map_threaded.push_back(j);
    }
}
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
    assert(mapping.size() > 0 && "Not possible number of electrons - no. of states < 1");
}

/// <summary>
/// 
/// </summary>
void IsingModel_sym::check_periodicity() {
#pragma omp parallel for
    for (int k = 0; k < N; k++) {
        std::vector<bool> base_vector(L), temp(L);
        int_to_binary(mapping[k], base_vector);
        u64 idx = 0;
        int period_EC = 0; // because it always eneter do-while loop
        int multiplicity_in_SEC = 1;
        do {
            std::rotate(base_vector.begin(), base_vector.begin() + 1, base_vector.end());
            idx = binary_to_int(base_vector);
            period_EC++;
        } while (idx != mapping[k]); 

        int_to_binary(mapping[k], base_vector);
        temp = base_vector;
        u64 min_R_RX = INT_MAX;
        std::reverse(temp.begin(), temp.end());
        u64 minR = binary_to_int(temp), idx_R = minR;
        do {
            std::rotate(temp.begin(), temp.begin() + 1, temp.end());
            idx = binary_to_int(temp);
            if (idx <= minR)
                minR = idx;
        } while (idx != idx_R);
        u64 minX = INT_MAX, minRX = INT_MAX;
        if (h == 0) {
            temp = base_vector;
            // check flip and reverse
            temp.flip();
            minRX = binary_to_int(temp);
            u64 idx_RX = minRX;
            do {
                std::rotate(temp.begin(), temp.begin() + 1, temp.end());
                idx = binary_to_int(temp);
                if (idx <= minRX)
                    minRX = idx;
            } while (idx != idx_RX);

            // check flip
            std::reverse(temp.begin(), temp.end());
            minX = binary_to_int(temp);
            u64 idx_X = minX;
            do {
                std::rotate(temp.begin(), temp.begin() + 1, temp.end());
                idx = binary_to_int(temp);
                if (idx <= minX)
                    minX = idx;
            } while (idx != idx_X);

        }
        else
            multiplicity_in_SEC = 1;
        if (mapping[k] == minRX) multiplicity_in_SEC++;
        this->periodicity[k] = period_EC * multiplicity_in_SEC;
    }
}

/// <summary>
/// 
/// </summary>
/// <param name="k"></param>
/// <param name="value"></param>
/// <param name="temp"></param>
void IsingModel_sym::setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) {
    u64 min = INT_MAX, idx = 0;
    u64 tmp_idx = binary_to_int(temp); // saving the index of current temporary after H action
    // checking periodicity of left vector
    do {
        std::rotate(temp.begin(), temp.begin() + 1, temp.end());
        idx = binary_to_int(temp);
        if (idx <= min)
            min = idx;
    } while (idx != tmp_idx);

    u64 min_R_RX = INT_MAX;
    std::vector<bool> base_vector(L);
    int_to_binary(tmp_idx, temp);
    base_vector = temp;
    std::reverse(temp.begin(), temp.end());
    u64 minR = binary_to_int(temp), idx_R = minR;
    do {
        std::rotate(temp.begin(), temp.begin() + 1, temp.end());
        idx = binary_to_int(temp);
        if (idx <= minR)
            minR = idx;
    } while (idx != idx_R);

    if (h == 0) {
        temp = base_vector;
        // check flip and reverse
        temp.flip();
        u64 minRX = binary_to_int(temp), idx_RX = minRX;
        do {
            std::rotate(temp.begin(), temp.begin() + 1, temp.end());
            idx = binary_to_int(temp);
            if (idx <= minRX)
                minRX = idx;
        } while (idx != idx_RX);

        // check flip
        std::reverse(temp.begin(), temp.end());
        u64 minX = binary_to_int(temp), idx_X = minX;
        do {
            std::rotate(temp.begin(), temp.begin() + 1, temp.end());
            idx = binary_to_int(temp);
            if (idx <= minX)
                minX = idx;
        } while (idx != idx_X);

        min_R_RX = std::min(minX, minRX);
    }
    else
        min_R_RX = tmp_idx;
    min_R_RX = std::min(minR, min_R_RX);
    idx = std::min(min, tmp_idx);
    idx = std::min(idx, min_R_RX);
    idx = binary_search(mapping, 0, N - 1, idx);

    value = value * sqrt((double)periodicity[k] / (double)periodicity[idx]);
    H(idx, k) += value;
    H(k, idx) += value;
}
/// <summary>
/// 
/// </summary>
void IsingModel_sym::hamiltonian() {
    H.zeros();
    std::vector<bool> base_vector(L);
    std::vector<bool> temp(base_vector); // changes under H action
    for (u64 k = 0; k < N; k++) {
        int_to_binary(mapping[k], base_vector);
        int s_i, s_j;
        for (int j = 0; j <= L - 1; j++) {
            s_i = static_cast<double>(base_vector[j]);

            /* transverse field */
            temp = base_vector;
            temp[j] = ~base_vector[j];
            setHamiltonianElem(k, g, std::move(temp));

            /* disorder */
            H(k, k) += this->h * (s_i - 0.5);

            if (nearest_neighbours[j] >= 0) {
                /* Ising-like spin correlation */
                s_j = static_cast<double>(base_vector[nearest_neighbours[j]]);
                H(k, k) += this->J[j] * (s_i - 0.5) * (s_j - 0.5);
            }
        }
    }
}