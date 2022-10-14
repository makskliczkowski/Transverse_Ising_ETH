
#if defined(HEISENBERG) || defined(XYZ)
	const double S = 0.5;
#else
	const double S = 1.0;
#endif

namespace operators{
    
    inline
    std::pair<cpx, u64> 
    sigma_x(u64 base_vec, int L, std::vector<int> sites) {
		for (auto& site : sites) {
			//site = properSite(site);
			base_vec = flip(base_vec, BinaryPowers[L - 1 - site], L - 1 - site);
		}
		return std::make_pair(S, base_vec);
	};

    inline
	std::pair<cpx, u64> 
    sigma_y(u64 base_vec, int L, std::vector<int> sites) {
		auto tmp = base_vec;
		cpx val = 1.0;
		for (auto& site : sites) {
			//site = properSite(site);
			val *= S * (checkBit(tmp, L - 1 - site) ? im : -im);
			tmp = flip(tmp, BinaryPowers[L - 1 - site], L - 1 - site);
		}
		return std::make_pair(val, tmp);
	};
	
    inline
    std::pair<cpx, u64> 
    sigma_z(u64 base_vec, int L, std::vector<int> sites) {
		auto tmp = base_vec;
		double val = 1.0;
		for (auto& site : sites) {
			//site = properSite(site);
			val *= checkBit(tmp, L - 1 - site) ? S : -S;
		}
		return std::make_pair(val, base_vec);
	};

}