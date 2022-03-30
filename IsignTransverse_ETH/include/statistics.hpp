#pragma once

namespace statistics{

//! ---------------------------------------------------------------- IPR
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
		cpx c_k = cdot(new_basis(k), _state);
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
		
        if (abs(delta_n) <= 1e-15) 
            assert(false && "Degeneracy!!!\n");
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
//<! spectral unfolding in-place
inline
void unfolding(arma::vec& eigenvalues){

}

//<! spectral unfolding as return
inline
arma::vec unfolding(const arma::vec& eigenvalues){
    const size_t N = eigenvalues.size();
    arma::vec unfolded_energies(N, arma::fill::zeros);


    return unfolded_energies;
}

//! ---------------------------------------------------------------- DISTRIBUTION FUNCTIONS AND PROBABILITIES
//<! Creates a probabilty distribution of data and saves it in the directory
inline
void probability_distribution(
    std::string dir,        //<! directory to save distribution
    std::string name,       //<! name of distribution
    const arma::vec& data,  //<! random data
    int n_bins = -1         //<! number of bins in the histogram
    ) {
	if (n_bins <= 0)
		n_bins = 1 + long(3.322 * log(data.size()));
	const double _min = arma::min(data);
	const double _max = arma::max(data);
	auto prob_dist = normalise_dist(
        arma::conv_to<arma::vec>::from(
            arma::hist(data, n_bins)
            ), _min, _max);
	auto x = arma::linspace(_min, _max, prob_dist.size());

    createDirs(dir);
    save_to_file(dir + name + ".dat", x, prob_dist);
}

//<! Creates a probabilty distribution of data and saves it in the directory
inline
auto probability_distribution(
    const arma::vec& data,  //<! random data
    int n_bins = -1         //<! number of bins in the histogram
    ) -> arma::vec
    {
    if (n_bins <= 0)
		n_bins = 1 + long(3.322 * log(data.size()));
	const double _min = arma::min(data);
	const double _max = arma::max(data);
	return normalise_dist(
        arma::conv_to<arma::vec>::from(
            arma::hist(data, n_bins)
            ), _min, _max);
    }

//<! data fluctuations around average structure
inline
auto data_fluctuations(
    const arma::vec& data,  //<! data to find fluctuations
    int mu                  //<! number of datapoints in bucket to find structure
) -> arma::vec 
{
	arma::vec fluct(data.size(), arma::fill::zeros);
	assert(mu < data.size() && "Bucket exceeds data container\nTry again\n");

#pragma omp parallel for
	for (int k = mu / 2; k < data.size() - mu / 2; k++) {
		double average = 0;
        int num = mu / 2;
		for (int n = k - num; n < k + num; n++)
			average += data(n);
		fluct(k - mu / 2) = data(k) - average / double(mu);
	}
	return fluct;
}

//<! Calculates a quantity similar to spectrum repulsion and finds
//<! the average of it and the biggest outliers
inline
auto statistics_average(
    const arma::vec& diag_mat_elem, //<! diagonal matrix element of observable (in eigenbasis)
    int num_of_outliers             //<! number of outliers of data to print
) -> arma::vec
{
	std::vector<double> spec_rep(num_of_outliers + 1, INT_MIN);
	double average = 0;
	for (int k = 1; k < diag_mat_elem.size(); k++) {
		double repulsion = abs(diag_mat_elem(k) - diag_mat_elem(k - 1));
		average += repulsion;
		int i = num_of_outliers + 1;
		while (i > 1) {
			if (repulsion < spec_rep[i - 1]) break;
			i--;
		}
		if (i > 0 && i <= num_of_outliers) {
			spec_rep.insert(spec_rep.begin() + i, repulsion);
			spec_rep.pop_back();
		}
	}
	spec_rep[0] = average / double(diag_mat_elem.size() - 2);
	return (arma::vec)spec_rep;
}


};