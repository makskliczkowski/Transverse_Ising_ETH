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


int main(const int argc, char* argv[]) {
	auto start = std::chrono::high_resolution_clock::now();

	int L = 6;
	double g = 0.5;
	double h = 0.0;
	double w = 0.0;


	std::vector<double> J(L, 1);
		std::unique_ptr<IsingModel> Hamil(new IsingModel_disorder(L, J, g, h, w));
		Hamil->diagonalization();
		u64 N = Hamil->get_hilbert_size();
		std::vector<double> vectttt = quantum_fidelity(0 * N / 20.0, 0 * N / 20.0+1, L, J, g, h, 0.02);
	mat corr_mat = Hamil->correlation_matrix(0);
	out << "<Sz(i)Sz(j)> = \n" << corr_mat << endl << endl;
	out << "S = " << Hamil->total_spin(corr_mat) << endl;
	out << "S^2 = " << arma::accu(corr_mat) << endl;
	vec Sq(L + 1, fill::zeros);
#pragma omp parallel for shared (Sq)
	for (int k = 0; k <= L; k++) {
		double q;
		if (!_BC) q = 2 * k * pi / (double)L;
		else q = k * pi / ((double)L + 1.0);
		cpx Sq_temp;
		for (int i = 0; i < L; i++) {
			for (int j = 0; j < L; j++) {
				Sq_temp += std::exp(1i * q * double(j - i)) * corr_mat(i, j);
			}
		}
		Sq(k) = real(2.0 * Sq_temp / pi / (L + 1.0));
	}
	out << Sq << endl;


	auto stop1 = std::chrono::high_resolution_clock::now();
	out << "Time to finish simulation with symmetry: "<<\
		double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(stop1 - start)).count()) / 1000.0 << " seconds" << endl;
	
	auto stop2 = std::chrono::high_resolution_clock::now();
	out << double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(stop2 - stop1)).count()) / 1000.0 << " seconds" << endl;
	return 0;
}

