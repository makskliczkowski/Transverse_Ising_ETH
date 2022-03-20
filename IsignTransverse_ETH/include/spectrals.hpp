// ---------------------------------------------------------------------------------- RESPONSE FUNCTION
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
	openFile(reponse_fun, name + alfa.get_info({}) + "_M=" + std::to_string(M) + ".dat", ios::out);
	// long int M = std::pow(N, 0.75);
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
	auto omegas = arma::logspace(std::floor(log10(wH)) - 1, 2, 300);
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
// ----------------------------------------- AUTOCORRELATION FUNCTION
//<! time evolution of autocorrelation function
template <typename _type>
inline 
void timeEvolution(
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
	auto times = arma::logspace(-2, t_max, 500);
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
auto timeEvolution(
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


// ---------------------------------------------------------------------------------- SPECTRAL FORM FACTOR
// ----------------------------------------- MEAN LEVEL SPACING
//<! mean level spacing between iterators
template <typename iterator_type>
[[nodiscard]]
inline
double
mean_level_spacing(
    iterator_type begin,  //<! first iterator to consider
    iterator_type end     //<! last iterator
    ){
	double omega_H = 0;
    u64 size = std::distance(begin, end);
#pragma omp parallel for reduction(+: omega_H)
	for (auto it = begin; it != end; ++it) {
		omega_H += *std::next(it) - *it;
	}
	return omega_H / double(size);
}

//<! mean level spacing for whole eigenvalue array
[[nodiscard]]
inline
double
mean_level_spacing(const arma::vec& eigenvalues)
{
    const size_t N = eigenvalues.size();
	const double chi = 0.341345;
	double trace_H2 = 0;
	double trace_H = 0;
#pragma omp parallel for reduction(+: trace_H, trace_H2)
	for (int k = 0; k < N; k++) {
		trace_H += eigenvalues(k);
		trace_H2 += eigenvalues(k) * eigenvalues(k);
	}
	return sqrt(trace_H2 / double(N) - trace_H * trace_H / double(N * N)) / (chi * N);
}

// ----------------------------------------- SPECTRAL FORM FACTOR (SFF)
[[nodiscard]]
inline 
double 
spectral_structure_factor_folded(
    const arma::vec& eigenvalues,   //<! eigenvalues to generate SFF
    double t                        //<! time point at which SFF calculated
    ){
    const size_t N = eigenvalues.size();
	double ssf_re = 0, ssf_im = 0;
//#pragma omp parallel for reduction(+: ssf_re, ssf_im)
	for (long n = 0; n < N; n++) {
		cpx ssf = std::exp(-im * eigenvalues(n) * t);
		ssf_re += real(ssf);
		ssf_im += imag(ssf);
	}
	double ssf = abs(cpx(ssf_re, ssf_im));
	ssf *= ssf;
	return ssf / double(N);
}

[[nodiscard]]
inline
arma::vec
spectral_structure_factor_folded(
    const arma::vec& eigenvalues,   //<! eigenvalues to generate SFF
    const arma::vec& times          //<! time range to calculate within
    ){
	arma::vec ssf(times.size(), arma::fill::zeros);
#pragma omp parallel for
	for (long i = 0; i < ssf.size(); i++)
		ssf(i) = spectral_structure_factor_folded(eigenvalues, times(i));
	return ssf;
}
