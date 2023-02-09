#pragma once

namespace statistics{

//! ---------------------------------------------------------------- IPR
//<! calculate participation ratio of input state for any q (q=1 is the typically used one)
template <typename _type> 
[[nodiscard]]
inline
double participation_ratio_plane_wave(
    const arma::Col<_type>& _state,   //<! input state
    double q = 1.0
    ) {
    double pr = 0;
    const size_t N = _state.size();
#pragma omp parallel for reduction(+: pr)
	for (int n = 0; n < N; n++) {
        cpx overlap = 0;
        double k = two_pi * n / double(N);
        for(int l = 0; l < N; l++)
            overlap += std::exp(1i * k * double(l)) * _state(l);
		double value = abs(conj(overlap) * overlap);
		pr += std::pow(value, q);
	}
	return pr  / std::pow(double(N), q);
}

//<! calculate participation ratio of input state for any q (q=1 is the typically used one)
template <typename _type> 
[[nodiscard]]
inline
double participation_ratio(
    const arma::Col<_type>& _state,   //<! input state
    double q = 1.0
    ) {
    double pr = 0;
    const size_t N = _state.size();
#pragma omp parallel for reduction(+: pr)
	for (int n = 0; n < N; n++) {
		double value = abs(conj(_state(n)) * _state(n));
		pr += std::pow(value, q);
	}
	return pr;
}

//<! calculate inverse participation ratio of input state
template <typename _type> 
[[nodiscard]]
inline
double inverse_participation_ratio(
    const arma::Col<_type>& _state   //<! input state
    ) {
    double ipr = 0;
    const size_t N = _state.size();
#pragma omp parallel for reduction(+: ipr)
	for (int n = 0; n < N; n++) {
		double value = abs(conj(_state(n)) * _state(n));
		ipr += value * value;
	}
	return 1.0 / ipr;
}


//! ---------------------------------------------------------------- INFORMATION ENTROPY
//<! calculate information entropy of input state in computational basis (full)
template <typename _type> 
[[nodiscard]]
inline
double information_entropy(
    const arma::Col<_type>& _state   //<! inout state
    ) {
	double ent = 0;
    const size_t N = _state.size();
#pragma omp parallel for reduction(+: ent)
	for (int k = 0; k < N; k++) {
		double val = abs(conj(_state(k)) * _state(k));
		ent += val * log(val);
	}
	return -ent / log(0.48 * N);
}

//<! calculate information entropy of input state in another eigenbasis set within range
template <typename _type>
[[nodiscard]]
inline 
double information_entropy(
    const arma::Col<_type>& _state,     //<! inout state
    const arma::Mat<_type>& new_basis,  //<! new eigenbasis to find overlap with _state
    u64 _min,                           //<! first state in new basis
    u64 _max                            //<! last eigenstate in new basis
    ) {
    const size_t N = _state.size();
	double ent = 0;
#pragma omp parallel for reduction(+: ent)
	for (long k = (long)_min; k < (long)_max; k++) 
    {
		cpx c_k = cdot(new_basis.col(k), _state);
		double val = abs(conj(c_k) * c_k);
		ent += val * log(val);
	}
	return -ent / log(0.48 * N);
}
//<! same but without range (taken full space)
template <typename _type>
[[nodiscard]]
inline 
double information_entropy(
    const arma::Col<_type>& _state,     //<! inout state
    const arma::Mat<_type>& new_basis  //<! new eigenbasis to find overlap with _state
    ) 
    { return information_entropy(_state, new_basis, 0, new_basis.n_cols); }


//! ---------------------------------------------------------------- LEVEL SPACING
//<! level spacing (gap ratio) between iterators (no bound checks!)
template <typename iterator_type>
[[nodiscard]]
inline
double eigenlevel_statistics(
    iterator_type begin,
    iterator_type end
    ) {
	const size_t size           = std::distance(begin, end);
    const iterator_type first   = std::next(begin);
    const iterator_type last    = std::prev(end);
    double E_prev = *begin;
	
    double r = 0;
//#pragma omp parallel for reduction(+: r)
    int counter = 0;
	for (auto it = first; it != last; ++it) 
    {
        double E_next = *std::next(it);
		
        const double delta_n        = (*it) - E_prev;
		const double delta_n_next   = E_next - (*it);

		const double min = std::min(delta_n, delta_n_next);
		const double max = std::max(delta_n, delta_n_next);
		
        if (abs(delta_n) <= 1e-15){ 
            std::cout << "Index: " << it << std::endl;
            for (auto it = first; it != last; ++it)
                std::cout << *it << std::endl; 
            assert(false && "Degeneracy!!!\n");
        }
		r += min / max;
        
        E_prev = (*it);
        counter++;
	}
	return r / double(counter);
}
//<! level statistics for full spectrum
inline
double eigenlevel_statistics(const arma::vec& energies)
    { return eigenlevel_statistics(energies.begin(), energies.end()); };

//<! level statistics with return for distribution calculation
template <typename iterator_type>
[[nodiscard]]
inline
arma::vec eigenlevel_statistics_return(
    iterator_type begin,
    iterator_type end
    ){
	const size_t size           = std::distance(begin, end);
    const iterator_type last    = end - 2;
    arma::vec gap_ratio(size - 2);
    int counter = 0;
	for (auto it = begin; it != last; ++it){
        gap_ratio(counter) = eigenlevel_statistics(it, it + 3); // *2 because
        counter++;
    }
    return gap_ratio;
};

//<! level statistics with return for distribution calculation (all eigenvalues) 
[[nodiscard]]
inline
arma::vec eigenlevel_statistics_return(const arma::vec& energies)
    { return eigenlevel_statistics_return(energies.begin(), energies.end()); }


//! ---------------------------------------------------------------- UNFOLDING
//<! spectral unfolding as return
[[nodiscard]]
inline
arma::vec unfolding(const arma::vec& eigenvalues, int n = 6){
    const size_t N = eigenvalues.size();

    // calculate cummulative distribution function (cdf)
    arma::vec cdf(eigenvalues.size(), arma::fill::zeros);
    std::iota(cdf.begin(), cdf.end(), 0);
    
    // fit polynomial order 10 to cdf
    auto p = arma::polyfit(eigenvalues, cdf, n);

    // evaluate fit at each energy: result is the unfolded energy
    arma::vec res = arma::polyval(p, eigenvalues);
    
    return res;
}

// ---------------------------------------------------------------------------------- SPECTRAL STATISTICS
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
// ----------------------------------------- TYPICAL LEVEL SPACING
//<! typical level spacing between iterators
template <typename iterator_type>
[[nodiscard]]
inline
double
typical_level_spacing(
    iterator_type begin,  //<! first iterator to consider
    iterator_type end     //<! last iterator
    ){
	double omega_H_typ = 0;
    u64 size = std::distance(begin, end);
#pragma omp parallel for reduction(+: omega_H_typ)
	for (auto it = begin; it != end; ++it) {
		omega_H_typ += log(*std::next(it) - *it);
	}
	return exp(omega_H_typ / double(size));
}
// ----------------------------------------- SPECTRAL FORM FACTOR (SFF)
//<! Finite-Temperature sff at time-point
[[nodiscard]]
inline 
std::pair<double, double> 
spectral_form_factor(
    const arma::vec& eigenvalues,   //<! eigenvalues to generate SFF
    double t,                       //<! time point at which SFF calculated
    double beta = 0.0               //<! inverse temperature 
    ){
    const size_t N = eigenvalues.size();
	double sff_re = 0, sff_im = 0;
    double Z = 0;
	for (long n = 0; n < N; n++) {
        Z += std::exp(- beta * eigenvalues(n));
		cpx sff = std::exp(- (beta + im * t) * eigenvalues(n));
		sff_re += real(sff);
		sff_im += imag(sff);
	}
	double sff = abs(cpx(sff_re, sff_im));
	sff = sff * sff / (Z * Z);
	return std::make_pair(sff, 1.0 / double(N));
}

//<! sff at time-point with gaussian filter
[[nodiscard]]
inline 
std::pair<double, double>
spectral_form_factor_filter(
    const arma::vec& eigenvalues,   //<! eigenvalues to generate SFF
    double t,                       //<! time point at which SFF calculated
    double eta = 0.5                //<! filter controling fractioon of eigenstates
    ){
    const size_t N = eigenvalues.size();
	double sff_re = 0, sff_im = 0;
    const double mean = arma::mean(eigenvalues);
    const double stddev = arma::stddev(eigenvalues);
    const double denom = 2.0 * eta * eta * stddev * stddev;
    double Z = 0;
	for (long n = 0; n < N; n++) {
        const double filter = exp( -(eigenvalues(n) - mean) * (eigenvalues(n) - mean) / denom );
        Z += abs(filter * filter);
		cpx sff = filter * std::exp(-im * double(two_pi) * eigenvalues(n) * t);
		sff_re += real(sff);
		sff_im += imag(sff);
	}
	double sff = abs(cpx(sff_re, sff_im));
	sff *= sff;
	return std::make_pair(sff, Z);
}

//<! sff for time range
[[nodiscard]]
inline
std::pair<arma::vec, double>
spectral_form_factor(
    const arma::vec& eigenvalues,   //<! eigenvalues to generate SFF
    const arma::vec& times,         //<! time range to calculate within
    double beta = 0.0,              //<! inverse temperature
    double eta = 0.5                //<! filter parameter
    ){
    double Z = 0;
	arma::vec sff(times.size(), arma::fill::zeros);
#pragma omp parallel for
	for (long i = 0; i < sff.size(); i++){
        double sff_temp = 0;
		if((eta < 0.0 || eta > 1.0) || beta > 0)
            std::tie(sff_temp, Z) = spectral_form_factor(eigenvalues, times(i), beta);
        else
            std::tie(sff_temp, Z) = spectral_form_factor_filter(eigenvalues, times(i), eta);
        sff(i) = sff_temp;
    }
	return std::make_pair(sff, Z);
}


};