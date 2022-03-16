#include <armadillo>
#include <cassert> 
class foo {
public:
	arma::mat matrix;
	arma::mat eigenvectors;
	arma::vec eigenvalues;

	int _L;
	size_t N;
	foo(int L) : _L(L) {
		this->N = (1ULL << this->_L);
		this->matrix = arma::randu<arma::mat>(this->N, this->N);
		this->matrix += this->matrix.t();
	};
	~foo() = default;
	void diagonalization() {
		std::string message = "dim = " + std::to_string(matrix.n_cols * sizeof(matrix(0, 0)));
		try {
			arma::eig_sym(this->eigenvalues, this->eigenvectors, this->matrix);
		}
		catch (const std::runtime_error& err) {
			std::cout << "Runtime error:\t" << err.what() << "\n";
			std::cout << message << std::endl;
			assert(false);
		}
		catch (const std::bad_alloc& err) {
			std::cout << "Bad alloc error:\t" << err.what() << "\n";
			std::cout << message << std::endl;
			assert(false);
		}
		catch (const std::exception& err) {
			std::cout << "Exception:\t" << err.what() << "\n";
			std::cout << message << std::endl;
			assert(false);
		}
		catch (...) {
			std::cout << "Unknown error...!" << "\n";
			std::cout << message << std::endl;
			assert(false);
		}
	}
};


int main(const int argc, char* argv[]) {
	int L = 10;
	foo obj(L);
	obj.diagonalization();
	for (int k = 0; k < 10; k++) {
		std::cout << obj.eigenvalues(k) << std::endl;
	}
	return 0;
}