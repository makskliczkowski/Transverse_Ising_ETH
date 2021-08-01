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

double w = 0.2;


int main(const int argc, char* argv[]) {
	cout << "am ready" << endl;

	int L = 8;
	double g = 0.5;
	double h = 0.2;
	int averages = 1;

	for (L = 2; L <= 12; L+=2) {
		std::vector<double> J(L);
		std::fill(J.begin(), J.end(), 1.0);
		std::unique_ptr<IsingModel_disorder> A(new IsingModel_disorder(L, J, g, h));
		double ipr1 = 0, ipr2 = 0;
		A->diagonalization();
		for (int k = 0; k < A->get_hilbert_size() - 1; k++) {
			if (abs(A->get_eigenvalues()(k) - A->get_eigenvalues()(k + 1)) <= 1e-10)
				out << " fock" << endl;
		}
		for (int m = 0; m < averages; m++) {
			A->hamiltonian();
			A->diagonalization();
			ipr1 += A->ipr(0) / (double)averages;
			ipr2 += A->ipr(static_cast<int>((double)A->get_hilbert_size() / 2.0)) / (double)averages;
		}
		out << ipr1 << "\t\t" << ipr2 << endl;
	}
	return 0;
}