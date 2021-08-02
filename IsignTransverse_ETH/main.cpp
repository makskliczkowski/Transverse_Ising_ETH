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

double w = 0.0;


int main(const int argc, char* argv[]) {
	cout << "am ready" << endl;

	int L = 8;
	double g = 0.5;
	double h = 0.2;
	std::vector<double> J(L);
	std::fill(J.begin(), J.end(), 1.0);
	int averages = 1;

	std::unique_ptr<IsingModel_sym> A(new IsingModel_sym(L, J, g, h));
	A->diagonalization();
	out << A->get_eigenvalues()(0) << endl;

		std::unique_ptr<IsingModel_disorder> B(new IsingModel_disorder(L, J, g, h));
		
		B->diagonalization();
		out << B->get_eigenvalues()(0) << endl;

	return 0;
}