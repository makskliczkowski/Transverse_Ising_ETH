#pragma once
namespace spectrals{
	
	// ---------------------------------------------------------------------------------- RESPONSE FUNCTION
	//<! set the omega bins for calculating the spectral function
	//<! class with energy differences for given eigenvalues,
	class preset_omega {
	private:
		arma::vec eigenvalues;
		std::vector<double> energy_diferences;
		std::vector<size_t> idx_alfa;
		std::vector<size_t> idx_beta;
		double tol = 2.0;			//<! width of antidiagonal to get matrix elements
		double E_av = 0.0;			//<! mean energy (E_n + E_m)/2 ~ E_av 
	public:
		preset_omega() = default;
		explicit preset_omega(
			const arma::vec& _eigenvalues,	//<! input eigenvalues
			double tolerance,				//<! width of antidiagonal
			double mean_energy				//<! mean energy
			) 
				: eigenvalues(_eigenvalues), tol(tolerance), E_av(mean_energy) 
				{ set_elements(); }
		
		//<! method to fill energy_differences and indices according to input settings
		void set_elements()
		{
			const u64 N = this->eigenvalues.size();
			this->idx_beta = std::vector<size_t>(); 
			this->idx_alfa = std::vector<size_t>();
			this->energy_diferences = std::vector<double>();
			for (long int i = 0; i < N; i++)
			{
				for (long int j = 0; j < N && j != i; j++)
				{
					if (abs((this->eigenvalues(j) + this->eigenvalues(i)) / 2. - this->E_av) < this->tol / 2.)
					{
						this->idx_alfa.push_back(i);
						this->idx_beta.push_back(j);
						this->energy_diferences.push_back(abs(this->eigenvalues(j) - this->eigenvalues(i)));
					}
				}
			}
			auto permut = sort_permutation(this->energy_diferences, [](const double a, const double b)
										   { return a < b; });
			apply_permutation(this->energy_diferences, permut);
			apply_permutation(this->idx_beta, permut);
			apply_permutation(this->idx_alfa, permut);
		}

		//<! get all elements
		auto get_elements() const
			{ return std::make_tuple(this->energy_diferences, this->idx_alfa, this->idx_beta); }

		//<! get number of energy differences
		auto get_size() const
			{ return this->energy_diferences.size(); }
	};

	//<! calculate response function for input model on self-built 'log' scale 
	//<! (fixed number of matrix elementes in omega bucket)
	inline 
	auto spectralFunction(
	    const arma::cx_mat &mat_elem,   	//<! input calculated matrix elements for any operator
		const preset_omega& input_omegas,	//<! input omegas and indices as tuple
		int M								//<! number of matrix elements in omega bucket
	    ) -> std::pair<arma::vec, arma::vec>
		{
			
			auto [energy_diff, idx_alfa, idx_beta] = input_omegas.get_elements();
			long int size = (int)energy_diff.size();
			long int bucket_num = int(size / (double)M);
			arma::vec response_fun(bucket_num + 1, arma::fill::zeros);
			arma::vec omegas(bucket_num + 1, arma::fill::zeros);
			for (int k = 0; k < bucket_num; k++)
			{
				double element = 0;
				double omega = 0;
				for (long int p = k * M; p < (k + 1) * M; p++)
				{
					cpx overlap = mat_elem(idx_alfa[p], idx_beta[p]);
					element += abs(overlap * overlap);
					omega += energy_diff[p];
				}
				response_fun(k) = element / double(M);
				omegas(k) = omega / double(M);
			}
			double element = 0;
			double omega = 0;
			int counter = 0;
			for (long int p = bucket_num * M; p < size; p++)
			{
				cpx overlap = mat_elem(idx_alfa[p], idx_beta[p]);
				element += abs(overlap * overlap);
				omega += energy_diff[p];
				counter++;
			}
			response_fun(bucket_num) = element / double(counter);
			omegas(bucket_num) = omega / (double)counter;
			return std::make_pair(omegas, response_fun);
		}

	//<! calculate response function for input model on self-built 'log' scale 
	//<! (fixed number of matrix elementes in omega bucket)
	template <typename _type>
	inline 
	void spectralFunction(
	    const IsingModel<_type> &alfa,  //<! input model (diagonalized)
	    const arma::cx_mat &mat_elem,   //<! input calculated matrix elements for any operator
	    std::string name                //<! filename to write data (contains directory)
	    ){
		const u64 N = alfa.get_hilbert_size();
		v_1d<long int> idx_beta, idx_alfa; // indices satysfying first condition in sum
		v_1d<double> energy_diff;		   // energy differnece(omega) of the above indices
		const double tol = 0.25 * alfa.L;
		for (long int i = 0; i < N; i++)
		{
			for (long int j = 0; j < N && j != i; j++)
			{
				if (abs((alfa.get_eigenEnergy(j) + alfa.get_eigenEnergy(i)) / 2. - alfa.get_eigenEnergy(alfa.E_av_idx)) < tol / 2.)
				{
					idx_alfa.push_back(i);
					idx_beta.push_back(j);
					energy_diff.push_back(abs(alfa.get_eigenEnergy(j) - alfa.get_eigenEnergy(i)));
				}
			}
		}
		auto permut = sort_permutation(energy_diff, [](const double a, const double b)
									   { return a < b; });
		apply_permutation(energy_diff, permut);
		apply_permutation(idx_beta, permut);
		apply_permutation(idx_alfa, permut);
		long int size = (int)energy_diff.size();
		// from L=12
		v_1d<int> Mx = v_1d<int>({100, 400, 700, 1200, 2000, 5000, 8000, 12000, 15000, 6000, 8000, 10000, 12000, 14000, 16000});
		NO_OVERFLOW(int M = Mx[alfa.L - 8];);
		std::ofstream reponse_fun;
		openFile(reponse_fun, name + alfa.get_info({}) + ".dat", ios::out);
		long int bucket_num = int(size / (double)M);
		for (int k = 0; k < bucket_num; k++)
		{
			double element = 0;
			double omega = 0;
			for (long int p = k * M; p < (k + 1) * M; p++)
			{
				cpx overlap = mat_elem(idx_alfa[p], idx_beta[p]);
				element += abs(overlap * overlap);
				omega += energy_diff[p];
			}
			reponse_fun << omega / (double)M << "\t\t" << element / double(M) << endl;
			reponse_fun.flush();
		}
		double element = 0;
		double omega = 0;
		int counter = 0;
		for (long int p = bucket_num * M; p < size; p++)
		{
			cpx overlap = mat_elem(idx_alfa[p], idx_beta[p]);
			element += abs(overlap * overlap);
			omega += energy_diff[p];
			counter++;
		}
		reponse_fun << omega / (double)counter << "\t\t" << element / double(counter) << endl;
		reponse_fun.flush();
		reponse_fun.close();
	}

	// ---------------------------------------------------------------------------------- INTEGRATED RESPONSE FUNCTION
	//<! calculate integrated respone function and writes to file given by
	//<!  input 'name' on omega log scale
	template <typename _type>
	inline 
	void 
	integratedSpectralFunction(
	    const IsingModel<_type> &alfa,   //<! input model (diagonalized) 
	    const arma::cx_mat &mat_elem,    //<! input calculated matrix elements for any operator
	    std::string name                 //<! filename to write data (contains directory)
	    ){
		const u64 N = alfa.get_hilbert_size();
		const double wH = alfa.mean_level_spacing_analytical();
		auto omegas = arma::logspace(std::floor(log10(wH)) - 1, 2, 5000);
		std::ofstream reponse_fun;
		openFile(reponse_fun, name + alfa.get_info({}) + ".dat", ios::out);
		for (auto &w : omegas)
		{
			double overlap = 0.;
	#pragma omp parallel for reduction(+: overlap)
			for (long int n = 0; n < N; n++)
			{
				for (long int m = 0; m < N; m++)
				{
					const double w_nm = abs(alfa.get_eigenEnergy(n) - alfa.get_eigenEnergy(m));
					if (w >= w_nm)
						overlap += abs(mat_elem(n, m) * conj(mat_elem(n, m)));
				}
			}
			overlap *= 1. / double(N);
			if (w == omegas(0))
				printSeparated(reponse_fun, "\t", 12, true, w, overlap, wH);
			else
				printSeparated(reponse_fun, "\t", 12, true, w, overlap);
		}
		reponse_fun.close();
	}

	//<! same as above but with output data and not written to file
	template <typename _type>
	[[nodiscard]]
	inline
	auto integratedSpectralFunction(
	    const IsingModel<_type> &alfa,   //<! input model (diagonalized)  
	    const arma::cx_mat &mat_elem,    //<! input calculated matrix elements for any operator 
	    const arma::vec &omegas          //<! omega range to caluclate data on
	    ) -> arma::vec 
	    {
		const u64 N = alfa.get_hilbert_size();
		arma::vec intSpec(omegas.size());
		for (int i = 0; i < omegas.size(); i++)
		{
			const double w = omegas(i);
			double overlap = 0.;
	#pragma omp parallel for reduction(+: overlap)
			for (long int n = 0; n < N; n++)
			{
				for (long int m = 0; m < N; m++)
				{
					const double w_nm = abs(alfa.get_eigenEnergy(n) - alfa.get_eigenEnergy(m));
					if (w >= w_nm)
						overlap += abs(mat_elem(n, m) * conj(mat_elem(n, m)));
				}
			}
			intSpec(i) = overlap / double(N);
		}
		return intSpec;
	}

	// ---------------------------------------------------------------------------------- TIME EVOLUTION
	// ----------------------------------------- OUT-OF-TIME CORRELATION FUNCTION (OTOC)
#define OTOC out_o_time_correlator
	//<! out of time correlator: < A(t) * B >
	template <typename _type>
	inline 
	void out_of_time_correlator(
	    const IsingModel<_type> &alfa,  //<! input model (diagonalized) 
	    const arma::cx_mat &A,    		//<! evolving operator
	    const arma::cx_mat &B,    		//<! static operator
	    std::string name                //<! filename to write data (contains directory)
	    );
	// ----------------------------------------- AUTOCORRELATION FUNCTION
	//<! time evolution of autocorrelation function
	template <typename _type>
	inline 
	void autocorrelation_function(
	    const IsingModel<_type> &alfa,   //<! input model (diagonalized) 
	    const arma::cx_mat &mat_elem,    //<! input calculated matrix elements for any operator
	    std::string name                 //<! filename to write data (contains directory)
	    ){
		const u64 N = alfa.get_hilbert_size();
		const double tH = 1. / alfa.mean_level_spacing_analytical();

		std::ofstream tEvolution;
		openFile(tEvolution, name + alfa.get_info({}) + ".dat", ios::out);
		double norm_diag = 0;
	#pragma omp parallel for reduction(+: norm_diag)
		for (long int k = 0; k < N; k++)
		{
			cpx temp = mat_elem(k, k);
			norm_diag += abs(temp * temp);
		}
		norm_diag /= double(N);
		const int t_max = (int)std::ceil(std::log(tH));
		auto times = arma::logspace(-2, t_max, 5000);
		for (auto &t : times)
		{
			double overlap = 0.;
			if (t > 5 * tH)
				break;
	#pragma omp parallel for reduction(+: overlap)
			for (long int n = 0; n < N; n++)
			{
				overlap += abs(mat_elem(n, n) * conj(mat_elem(n, n)));
				for (long int m = n + 1; m < N; m++)
				{
					const double w_nm = alfa.get_eigenEnergy(n) - alfa.get_eigenEnergy(m);
					double value = std::cos(w_nm * t);
					overlap += 2. * abs(mat_elem(n, m) * conj(mat_elem(n, m))) * value;
				}
			}
			overlap *= 1. / double(N);
			// tEvolution << t << "\t\t" << overlap << "\t\t" << tH << std::endl;
			if (t == times(0))
				printSeparated(tEvolution, "\t", 12, true, t, overlap, tH, norm_diag);
			else
				printSeparated(tEvolution, "\t", 12, true, t, overlap);
		}
		tEvolution.close();
	}

	//<! time evolution of autocorrelation function with output data
	template <typename _type>
	[[nodiscard]]
	inline
	auto autocorrelation_function(
	    const IsingModel<_type> &alfa,   //<! input model (diagonalized)  
	    const arma::cx_mat &mat_elem,    //<! input calculated matrix elements for any operator 
	    const arma::vec &times           //<! time range to caluclate data on
	    ) -> std::pair<arma::vec, double>
	{
		const u64 N = alfa.get_hilbert_size();
		const double tH = 1. / alfa.mean_level_spacing_analytical();
		double LTA = 0;
	#pragma omp parallel for reduction(+: LTA)
		for (long int k = 0; k < N; k++)
		{
			cpx temp = mat_elem(k, k);
			LTA += abs(temp * temp);
		}
		LTA /= double(N);
		arma::vec timeEv(times.size(), arma::fill::zeros);
		for (long int k = 0; k < times.size(); k++)
		{
			auto t = times(k);
			double overlap = 0.;
	#pragma omp parallel for reduction(+: overlap)
			for (long int n = 0; n < N; n++)
			{
				overlap += abs(mat_elem(n, n) * conj(mat_elem(n, n)));
				for (long int m = n + 1; m < N; m++)
				{
					const double w_nm = alfa.get_eigenEnergy(n) - alfa.get_eigenEnergy(m);
					overlap += 2. * abs(mat_elem(n, m) * conj(mat_elem(n, m))) * std::cos(w_nm * t);
				}
			}
			overlap *= 1. / double(N);
			timeEv(k) = overlap;
		}
		return std::make_pair(timeEv, LTA);
	}

	
	// -----------------------------------------  QUANTUM QUENCH
	[[nodiscard]]
	inline
	auto time_evolution(
		const arma::cx_vec& init_state,	//<! initial state
		const arma::sp_cx_mat& op_mat,	//<! operator matrix to fine time evolution
		double time
	) -> double
	{
		double output = 0;


		return output;
	}


	// -----------------------------------------  GREEN'S FUNCTION AND DOS
	inline
	auto density_of_states(
		const arma::vec& energies,
		double dw,
		double eta = -1
	) -> arma::vec
	{
		if(eta <= 0) eta = dw / 5.;
		auto omegas = arma::regspace(1.25 * arma::min(energies), dw, 1.25 * arma::min(energies));
		arma::vec DOS(omegas.size(), arma::fill::zeros);
		for(int k = 0; k < omegas.size(); k++){
			double w = omegas(k);
			for(int n = 0; n < energies.size(); n++)
				DOS(k) += -1. / pi * imag(1. / (w + eta * 1i - (energies(n) - energies(0))));
		}
		return DOS;
	};
	
}