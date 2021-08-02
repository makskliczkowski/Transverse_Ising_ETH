#include "include/IsingModel.h"

IsingModel_sym::IsingModel_sym(int L, vector<double>& J, double g, double h){
	this->L = L; this->J = J; this->g = g; this->h = h;
    this->mapping = std::vector<u64>();
    generate_mapping();
    this->N = this->mapping.size();
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
    std::vector<bool> temp(L), temp_flipped(L); // temporary dirac-notation base vector
    u64 min = N, idx = 0;
    for (u64 j = start; j < stop; j++) {
        int_to_binary(j, temp);
        temp_flipped = temp;
        temp_flipped.flip();
        min = binary_to_int(temp_flipped);
        if (j < min) {
            while (next_permutation(temp.begin(), temp.end())) { // translation
                idx = binary_to_int(temp);
                if (idx < min)
                    min = idx;
            }
            if (min == j)
                map_threaded.push_back(j);
        }
        else continue;
    }
}
/// <summary>
/// 
/// </summary>
void IsingModel_sym::generate_mapping() {
    u64 start = 0, stop = static_cast<u64>(std::pow(2, L - 1)); // because parity reflects the other half
    if (num_of_threads == 1)
        mapping_kernel(start, stop, mapping, 0);
    else {
        //Threaded
        //std::vector<my_uniq_ptr> map_threaded(num_of_threads);
        std::vector<std::vector<u64>> map_threaded(num_of_threads);
        std::vector<std::thread> threads;
        threads.reserve(num_of_threads);
        for (int t = 0; t < num_of_threads; t++) {
            start = (u64)(std::pow(2, L - 1) / (double)num_of_threads * t);
            stop = ((t + 1) == num_of_threads ? (u64)std::pow(2, L - 1) : u64(std::pow(2, L - 1) / (double)num_of_threads * (double)(t + 1)));
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
/// <param name="k"></param>
/// <param name="value"></param>
/// <param name="temp"></param>
void IsingModel_sym::setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) {

}
/// <summary>
/// 
/// </summary>
void IsingModel_sym::hamiltonian() {
    H.zeros();
    for (u64 k = 0; k < N; k++) {
        std::vector<bool> base_vector(L);
        int_to_binary(mapping[k], base_vector);

        std::vector<bool> temp(base_vector); // changes under H action
        int s_i, s_j;
        for (int j = 0; j <= L - 1; j++) {
            s_i = static_cast<double>(base_vector[j]);

            /* transverse field */
            temp[j] = ~base_vector[j];
            setHamiltonianElem(k, g, std::move(temp));
            //dsdsa

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