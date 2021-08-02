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
    //std::vector<bool> temp(L);
    //for (int k = 0; k < N; k++) {
    //    int_to_binary(mapping[k], temp);
    //    out << temp << endl;
    //}
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
    int n = 1;
    //out << "A new thread joined tha party! from " << start << " to " << stop << endl;
    std::vector<bool> temp(L), temp_flipped(L), temp_reversed(L); // temporary dirac-notation base vector
    u64 minX = INT_MAX, idx = 0, minR = INT_MAX, minRX = INT_MAX, min = INT_MAX;
    for (u64 j = start; j < stop; j++) {
        int_to_binary(j, temp);

        u64 min_R_RX;
        /*if (h == 0) {
            temp_flipped = temp;
            temp_reversed = temp;
            std::reverse(temp_reversed.begin(), temp_reversed.end());
            //minR = std::min(j, binary_to_int(temp_reversed));

            // check both reverse (R) and flip (X)
            temp_reversed.flip();
            minRX = std::min(j, binary_to_int(temp_reversed));

            // check only flip (X)
            temp_flipped.flip();
            minX = std::min(j, binary_to_int(temp_flipped));
            min_R_RX = std::min(minX, minRX);
        }
        else*/
            min_R_RX = j;

        //check transaltion
         min = j;
        do {
            std::rotate(temp.begin(), temp.begin() + 1, temp.end());
            idx = binary_to_int(temp);
            if (idx <= min)
                min = idx;
        } while (idx != j);
           
            
        if (min == j)//(std::min(min, min_R_RX) == j)
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
        std::vector<bool> base_vector(L);
        int_to_binary(mapping[k], base_vector);
        u64 idx = 0;
        int period = 0; // because it always eneter do-while loop
        do {
            std::rotate(base_vector.begin(), base_vector.begin() + 1, base_vector.end());
            idx = binary_to_int(base_vector);
            period++;
        } while (idx != mapping[k]); 
        this->periodicity[k] = period;
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
    } while (idx != tmp_idx); // if it is tmp_idx again
    /*temp.flip();
    min = std::min(min, binary_to_int(temp));

    std::reverse(temp.begin(), temp.end());
    min = std::min(min, binary_to_int(temp));

    temp.flip(); */
    idx = std::min(min, tmp_idx);
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