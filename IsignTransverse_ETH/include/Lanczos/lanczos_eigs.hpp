#pragma once

namespace lanczos {

	template <typename _type>
	void Lanczos<_type>::diagonalization(
		const arma::Col<_type>& rand	//<! random input
	) 
	{
		try {
			this->build(rand);
		}
		catch (const std::bad_alloc& e) {
			std::cout << "Memory exceeded" << e.what() << "\n";
			assert(false);
		}
		// no need to catch here as lanczos_steps << 10 000
		arma::eig_sym(
			this->eigenvalues,
			this->eigenvectors,
			this->H_lanczos
		);
	}

	template <typename _type>
	void Lanczos<_type>::diagonalization() {
		arma::Col<_type> rand = this->ran.template create_random_vec<_type>(this->N);
		this->diagonalization(rand);
	}

}