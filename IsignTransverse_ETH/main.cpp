#include "include/IsingModel.h"


std::random_device rd;
std::mt19937::result_type seed = rd() ^ (
	(std::mt19937::result_type)
	std::chrono::duration_cast<std::chrono::seconds>(
		std::chrono::system_clock::now().time_since_epoch()
		).count() +
	(std::mt19937::result_type)
	std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::high_resolution_clock::now().time_since_epoch()
		).count());
std::mt19937_64 gen(seed);

double disorder_strength = 5.0;


int main(const int argc, char* argv[]) {
	auto start = std::chrono::high_resolution_clock::now();

	int L = 10;
	double g = 0.5;
	double h = 0.02;
	std::vector<double> J(L);
	std::fill(J.begin(), J.end(), 1.0);

	for (h = 0.0; h <= 5.0; h += 0.1) {
		std::unique_ptr<IsingModel> B(new IsingModel_sym(L, J, g, h));
		u64 N = B->get_hilbert_size();
		double r = 0;
		for (int av = 0; av < 10; av++) {
			B->hamiltonian();
			B->diagonalization();
			r += B->eigenlevel_statistics(N / 2, 6 * N / 10);
		}
		out << h << "\t\t" << r / 10.0 << endl;
	}

	//A->operator_av_in_eigenstates(&IsingModel::av_sigma_x, *A, 1, "results/sigma_x_average.txt", "\t\t");

	/*for (int L = 2; L <= 20; L += 2) {
		std::vector<double> J(L);
		std::fill(J.begin(), J.end(), 1.0);
		std::unique_ptr<IsingModel> A(new IsingModel_disorder(L, J, g, h));
		A->diagonalization();
		out << L << "\t\t" << A->spectrum_repulsion(&IsingModel::av_sigma_x, *A, 1) << endl;
	}*/
	auto stop1 = std::chrono::high_resolution_clock::now();
	out << "Time to finish simulation with symmetry: "<<\
		double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(stop1 - start)).count()) / 1000.0 << " seconds" << endl;
	
	auto stop2 = std::chrono::high_resolution_clock::now();
	out << double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(stop2 - stop1)).count()) / 1000.0 << " seconds" << endl;
	return 0;
}