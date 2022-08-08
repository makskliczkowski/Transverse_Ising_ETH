#pragma once

namespace statistics{
//! ---------------------------------------------------------------- DISTRIBUTION FUNCTIONS AND PROBABILITIES
//! ---------------------------------- HISTOGRAM OF COUNTS
//<! Creates a probabilty distribution of data and saves it in the directory
template <typename ... _ty>
inline
void probability_distribution(
    std::string dir,        //<! directory to save distribution
    std::string name,       //<! name of distribution
    const arma::vec& data,  //<! random data
    int n_bins,		        //<! number of bins in the histogram
	_ty... args
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
    save_to_file(dir + name + ".dat", x, prob_dist, args...);
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
//! ---------------------------------- HISTOGRAM OF COUNTS IN SPECIFIED BINS


//! ---------------------------------- FLUCTUATIONS
//<! data fluctuations around average structure
inline
auto remove_fluctuations(
    const arma::vec& data,  //<! data to find fluctuations
    int mu                  //<! number of datapoints in bucket to find structure
) -> arma::vec 
{
	arma::vec new_data = data;
	assert(mu < data.size() && "Bucket exceeds data container\nTry again\n");

#pragma omp parallel for
	for (int k = mu / 2; k < data.size() - mu / 2; k++) {
		double average = 0;
		for (int n = k - mu / 2; n < k + mu / 2; n++)
			average += data(n);
		new_data(k - mu / 2) =  average / double(mu);
	}
	return new_data;

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