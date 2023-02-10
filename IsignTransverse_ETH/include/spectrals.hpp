#pragma once

#ifndef __COMMONS
	#include "commons.hpp"
#endif

#ifndef ARMA
	#include "armadillo_wrapper.hpp"
#endif

namespace spectrals{

	//<! calculate index of mean energy <E> = Tr(H) / dim
	inline 
	auto get_mean_energy_index(
		const arma::vec& eigenvalues
		) -> u64
		{
			const u64 N = eigenvalues.size();
			double E_av = arma::trace(eigenvalues) / double(N);
			
			auto i = min_element(begin(eigenvalues), end(eigenvalues), [=](double x, double y) {
				return abs(x - E_av) < abs(y - E_av);
				});
			return i - eigenvalues.begin();
		};

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


		//<! set omega bins
		auto set_omega_bins(int num_of_points = 1000) const -> arma::vec 
		{
			long int size = (int)this->energy_diferences.size();
			long int M = int(size / double(num_of_points));
			long int bucket_num = int(size / (double)M);
			arma::vec omegas(bucket_num + 1, arma::fill::zeros);
			for (int k = 0; k <= bucket_num; k++)
			{
				double omega = 0;
				long p_end = (k == bucket_num)? size : (k + 1) * M;
				int counter = 0;
				for (long int p = k * M; p < p_end; p++){
					omega += this->energy_diferences[p];
					counter++;
				}
				omegas(k) = omega / (double)counter;
			}
			return omegas;
		}

		//--------------------------- GETTERS
		//<! get all elements
		auto get_elements() const
			{ return std::make_tuple(this->energy_diferences, this->idx_alfa, this->idx_beta); }

		//<! get number of energy differences
		auto get_size() const
			{ return this->energy_diferences.size(); }

		//<! get eigenvalues
		auto get_eigenvalues() const
			{ return this->eigenvalues; }


		//--------------------------- SAVERS
		template <typename _ty>
		inline
		auto save_matrix_elements(std::string filename, const arma::Mat<_ty>& mat_elem) const {
			arma::vec omegas = this->energy_diferences;
			arma::vec mat_elem_to_save = arma::vec(this->energy_diferences.size(), arma::fill::zeros);
			for(int i = 0; i < this->energy_diferences.size(); i++){
				double element = abs(mat_elem(this->idx_alfa[i], this->idx_beta[i]));
				mat_elem_to_save(i) = element * element;
			}
			omegas.save(arma::hdf5_name(filename + ".hdf5", "omegas"));
			mat_elem_to_save.save(arma::hdf5_name(filename + ".hdf5", "mat_elem",arma::hdf5_opts::append));
		}
	};

	//<! calculate response function for input model on self-built 'log' scale 
	//<! (fixed number of matrix elementes in omega bucket)
	template <typename matrix>
	inline 
	auto spectralFunction(
	    const matrix &mat_elem,				//<! input calculated matrix elements for any operator
		const preset_omega& input_omegas,	//<! input energy differences and indices as tuple
	    const arma::vec& omegas	        	//<! omega range to caluclate data on
	    ) -> arma::vec
		{
			std::vector<double> energy_diff;
			std::vector<size_t> idx_alfa, idx_beta;
			std::tie(energy_diff, idx_alfa, idx_beta) = input_omegas.get_elements();
			long int size = (int)energy_diff.size();

			arma::vec response_fun(omegas.size(), arma::fill::zeros);
			arma::vec counter(omegas.size(), arma::fill::zeros);
		#pragma omp parallel for
			for (long int p = 0; p < size; p++)
			{
				double element = abs(mat_elem(idx_alfa[p], idx_beta[p]) * conj(mat_elem(idx_alfa[p], idx_beta[p])));

				u64 idx = min_element(begin(omegas), end(omegas), [=](double x, double y) {
				return abs(x - energy_diff[p]) < abs(y - energy_diff[p]);
				}) - omegas.begin();
				#pragma omp critical
				{
					response_fun(idx) += element;
					counter(idx)++;
				}
			}
			counter.elem( arma::find(counter == 0) ).ones();
			return response_fun / counter;
		}

	//<! calculate response function for input model on self-built 'log' scale 
	//<! (fixed number of matrix elementes in omega bucket)
	template <typename matrix>
	inline 
	void spectralFunction(
	    const matrix& mat_elem,			//<! input calculated matrix elements for any operator
		const arma::vec& eigenvalues,	//<! eigenenergies
		double tolerance,				//<! width of antidiagonal
		double mean_energy,				//<! mean energy
	    std::string name                //<! filename to write data (contains directory)
	    ){
		spectrals::preset_omega set_omega(eigenvalues, tolerance, mean_energy);
		auto [energy_diff, idx_alfa, idx_beta] = set_omega.get_elements();
		int num_of_points = 3000;
		long int size = (int)energy_diff.size();
		long int M = int(size / double(num_of_points));
		long int bucket_num = int(size / (double)M);

		std::ofstream reponse_fun;
		openFile(reponse_fun, name + ".dat", ios::out);
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
			if(element < 1e-36) continue;
			printSeparated(reponse_fun, "\t", 16, true, omega / (double)M, element / double(M));
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
		printSeparated(reponse_fun, "\t", 16, true, omega / (double)counter, element / double(counter));
		reponse_fun.flush();
		reponse_fun.close();
	}

	//<! calculate response function for input model on self-built 'log' scale 
	//<! (fixed number of matrix elementes in omega bucket)
	inline 
	void spectralFunction(
		arma::vec& omegas,			//<! omega values
	    arma::vec& mat_elem,		//<! input calculated matrix elements for any operator
	    std::string name,           //<! filename to write data (contains directory)
		int num_of_points = 3000	//<! total number of points
	    ){

		std::cout << arma::min(omegas) << "\t" << arma::max(omegas) << std::endl;
		auto p = sort_permutation(omegas, [](const double a, const double b)
								   { return a < b; });
		std::cout << arma::min(omegas) << "\t" << arma::max(omegas) << std::endl;
		apply_permutation(omegas, p);
		apply_permutation(mat_elem, p);	
		
		long int size = (int)omegas.size();
		long int M = int(size / double(num_of_points));
		long int bucket_num = int(size / (double)M);
		
		std::ofstream reponse_fun;
		openFile(reponse_fun, name + ".dat", ios::out);
		for (int k = 0; k < bucket_num; k++)
		{
			double element = 0;
			double omega = 0;
			for (long int p = k * M; p < (k + 1) * M; p++)
			{
				element += mat_elem[p];
				omega += omegas[p];
			}
			if(element < 1e-32) continue;
			printSeparated(reponse_fun, "\t", 16, true, omega / (double)M, element / double(M));
			reponse_fun.flush();
		}
		double element = 0;
		double omega = 0;
		int counter = 0;
		for (long int p = bucket_num * M; p < size; p++)
		{
			element += abs(mat_elem[p]);
			omega += omegas[p];
			counter++;
		}
		printSeparated(reponse_fun, "\t", 16, true, omega / (double)counter, element / double(counter));
		reponse_fun.flush();
		reponse_fun.close();
	}
	// ---------------------------------------------------------------------------------- INTEGRATED RESPONSE FUNCTION
	//<! calculate integrated respone function and writes to file given by
	//<!  input 'name' on omega log scale
	template <typename matrix>
	inline 
	void 
	integratedSpectralFunction(
	    const matrix& mat_elem,			//<! input calculated matrix elements for any operator
		const arma::vec& eigenvalues,	//<! eigenenergies
	    std::string name                //<! filename to write data (contains directory)
	    ){
		const u64 N = eigenvalues.size();
		const u64 E_av_idx = spectrals::get_mean_energy_index(eigenvalues);
		const double wH = statistics::mean_level_spacing(eigenvalues.begin() + u64(E_av_idx - 0.25 * N), eigenvalues.begin() + u64(E_av_idx + 0.25 * N));
		auto omegas = arma::logspace(std::floor(log10(wH)) - 1, 2, 3000);
		std::ofstream reponse_fun;
		openFile(reponse_fun, name + ".dat", ios::out);
		for (auto &w : omegas)
		{
			double overlap = 0.;
	#pragma omp parallel for reduction(+: overlap)
			for (long int n = 0; n < N; n++)
				for (long int m = 0; m < N; m++)
					if (w >= abs(eigenvalues(n) - eigenvalues(m)))
						overlap += abs(mat_elem(n, m) * conj(mat_elem(n, m)));

			overlap *= 1. / double(N);
			if (w == omegas(0))
				printSeparated(reponse_fun, "\t", 12, true, w, overlap, wH);
			else
				printSeparated(reponse_fun, "\t", 12, true, w, overlap);
		}
		reponse_fun.close();
	}

	//<! same as above but with output data and not written to file
	template <typename matrix>
	[[nodiscard]]
	inline
	auto integratedSpectralFunction(
	    const matrix& mat_elem,			//<! input calculated matrix elements for any operator
		const arma::vec& eigenvalues, 	//<! eigenenergies 
	    const arma::vec &omegas         //<! omega range to caluclate data on
	    ) -> arma::vec 
	    {
		const u64 N = eigenvalues.size();
		arma::vec intSpec(omegas.size());
	#pragma omp parallel for
		for (int i = 0; i < omegas.size(); i++)
		{
			const double w = omegas(i);
			double overlap = 0.;
			for (long int n = 0; n < N; n++)
				for (long int m = 0; m < N; m++)
					if (w >= abs(eigenvalues(n) - eigenvalues(m)))
						overlap += abs(mat_elem(n, m) * conj(mat_elem(n, m)));
			intSpec(i) = overlap / double(N);
		}
		return intSpec;
	}

	// ---------------------------------------------------------------------------------- TIME EVOLUTION
	// ----------------------------------------- OUT-OF-TIME CORRELATION FUNCTION (OTOC)
#define OTOC out_of_time_correlator
	//<! out of time correlator: < A(t) * B >
	template <typename matrix>
	inline 
	void out_of_time_correlator(
	    const matrix& A,    			//<! evolving operator
	    const matrix& B,    			//<! static operator
		const arma::vec& eigenvalues, 	//<! eigenenergies 
	    std::string name                //<! filename to write data (contains directory)
	    );

	// ----------------------------------------- AUTOCORRELATION FUNCTION
	//<! time evolution of autocorrelation function
	template <typename matrix>
	inline 
	void autocorrelation_function(
	    const matrix& mat_elem,			//<! input calculated matrix elements for any operator
		const arma::vec& eigenvalues, 	//<! eigenenergies 
	    std::string name                //<! filename to write data (contains directory)
	    ){
		const u64 N = eigenvalues.size();
		auto E_av_idx = spectrals::get_mean_energy_index(eigenvalues);
		const double tH = 1. / statistics::mean_level_spacing(eigenvalues.begin() + u64(E_av_idx - 0.25 * N), eigenvalues.begin() + u64(E_av_idx + 0.25 * N));

		std::ofstream tEvolution;
		openFile(tEvolution, name + ".dat", ios::out);
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
					const double w_nm = eigenvalues(n) - eigenvalues(m);
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
	template <typename matrix>
	[[nodiscard]]
	inline
	auto autocorrelation_function(
	    const matrix& mat_elem,			//<! input calculated matrix elements for any operator
		const arma::vec& eigenvalues, 	//<! eigenenergies 
	    const arma::vec& times          //<! time range to caluclate data on
	    ) -> std::pair<arma::vec, double>
	{
		const u64 N = eigenvalues.size();
		double LTA = 0;
	#pragma omp parallel for reduction(+: LTA)
		for (long int k = 0; k < N; k++)
		{
			cpx temp = mat_elem(k, k);
			LTA += abs(temp * temp);
		}
		LTA /= double(N);
		arma::vec timeEv(times.size(), arma::fill::zeros);
	#pragma omp parallel for
		for (long int k = 0; k < times.size(); k++)
		{
			auto t = times(k);
			double overlap = 0.;
			for (long int n = 0; n < N; n++)
			{
				overlap += abs(mat_elem(n, n) * conj(mat_elem(n, n)));
				for (long int m = n + 1; m < N; m++)
				{
					const double w_nm = eigenvalues(n) - eigenvalues(m);
					overlap += 2. * abs(mat_elem(n, m) * conj(mat_elem(n, m))) * std::cos(w_nm * t);
				}
			}
			overlap *= 1. / double(N);
			timeEv(k) = overlap;
		}
		return std::make_pair(timeEv, LTA);
	}

	
	// -----------------------------------------  QUANTUM QUENCH
	template <typename matrix>[[nodiscard]]
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