
namespace adiabatics{

	template <typename _ty>
	inline
	auto 
	gauge_potential(
    	const arma::Mat<_ty>& mat_elem,
    	const arma::vec& eigenvalues,
    	int L
    ) -> std::tuple<double, double, double, arma::vec> 
	{
        const size_t N = eigenvalues.size();
		const double lambda = double(L) / double(N);
		const size_t mu = long(0.5 * N);

		double E_av = arma::trace(eigenvalues) / double(N);
		auto i = min_element(begin(eigenvalues), end(eigenvalues), [=](double x, double y) {
			return abs(x - E_av) < abs(y - E_av);
		});
		const long E_av_idx = i - begin(eigenvalues);
		long int E_min = E_av_idx - long(mu / 2);
		long int E_max = E_av_idx + long(mu / 2);

        double AGP = 0.0;
		double typ_susc = 0.0;
		double susc = 0.0;
		arma::vec susc_vec(N, arma::fill::zeros);
    #pragma omp parallel for reduction(+ : AGP, susc, typ_susc)
		for (long int i = 0; i < N; i++)
		{
			double susc_tmp = 0;
			for (long int j = 0; j < N && j != i; j++)
			{
				const double nominator = std::abs(mat_elem(i, j) * conj(mat_elem(i, j)));
				const double omega_ij = eigenvalues(j) - eigenvalues(i);
				const double denominator = omega_ij * omega_ij + lambda * lambda;
				
				AGP += omega_ij * omega_ij * nominator / (denominator * denominator);
				susc_tmp += nominator / (omega_ij * omega_ij);
			}
			susc_vec(i) = susc_tmp;
			if (susc_tmp > 0 && (i > E_min && i < E_max))
			{
				typ_susc += std::log(susc_tmp);
				susc += susc_tmp;
			}
		}
        return std::make_tuple(AGP / double(N), exp(typ_susc / double(mu)), susc / double(mu), susc_vec);
    }



};