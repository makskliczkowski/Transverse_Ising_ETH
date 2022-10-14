
#include "SpinOperators.hpp"
//<! namespace with hamiltonian generators in full Hilbert space
namespace hamiltonian{

    inline
    auto
    XYZ_nnn_OBC(int L,                          //<! system size
                double J1, double J2,           //<! kinetic energy terms for nearest and next-nearest neighbours
                double eta1, double eta2,       //<! XY anisotropy strength for nearest and next-nearest neighbours
                double delta1, double delta2,   //<! interaction for nearest and next-nearest neighbours
                double hz, double hx            //<! magnetic fields
        )
    {
        size_t N = ULLPOW(L);
        arma::sp_mat H(N, N);

        std::vector<std::vector<double>> parameters = {{J1 * (1 - eta1), J1 * (1 + eta1), J1 * delta1},
                                                        {J2 * (1 - eta2), J2 * (1 + eta2), J2 * delta2}
                                                    };
        for(auto& x : parameters)
            std::cout << x << std::endl;
        std::vector<op_type> XYZoperators = {operators::sigma_x, operators::sigma_y, operators::sigma_z };
        std::vector<int> neighbor_distance = {1, 2};
        for (size_t k = 0; k < N; k++) {
		    double s_i, s_j;
		    for (int j = 0; j <= L - 1; j++) {
                cpx val = 0.0;
                u64 op_k;
                std::tie(val, op_k) = operators::sigma_z(k, L, { j });
                double fieldZ = (j != L - 1)? hz : 0.5 * hz;
		    	H(op_k, k) += fieldZ * real(val);

                std::tie(val, op_k) = operators::sigma_x(k, L, { j });
                double fieldX = (j != 0)? hx : 0.5 * hx;
		    	H(op_k, k) += fieldX * real(val);

                for(int a = 0; a < neighbor_distance.size(); a++){
                    int r = neighbor_distance[a];
		    	    if (j < L - r && j >= r - 1) {
                        for(int b = 0; b < XYZoperators.size(); b++){
                            op_type op = XYZoperators[b];
			                auto [val1, op_k] = op(k, L, { j });
			                auto [val2, opop_k] = op(op_k, L, { j + r });
		    	    	    H(opop_k, k) += parameters[a][b] * real(val1 * val2);
                        }
		    	    }
                }
		    }
	    }

        return H;
    }
}