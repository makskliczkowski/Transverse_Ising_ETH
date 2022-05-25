#pragma once

#define try_alloc(code) try{ code;}\
						catch(const std::bad_alloc& e) {						\
							std::cout << "Memory exceeded" << e.what() << "\n";	\
							assert(false);										\
						}

namespace lanczos {
	inline
		void Lanczos::diagonalization()
	{
		try_alloc(this->build(););
		arma::eig_sym(
			this->eigenvalues,
			this->eigenvectors,
			this->H_lanczos
		);

	}
	inline
	void Lanczos::diagonalization(
		const arma::cx_vec& random	//<! random input
	) 
	{
		try_alloc(this->build(random););
		arma::eig_sym(
				this->eigenvalues,
				this->eigenvectors,
				this->H_lanczos
			);
		//if (!random.is_empty())
		//	this->initial_random_vec = random;
		//auto eigs = [this](arma::cx_mat& V) {
		//	arma::eig_sym(
		//		this->eigenvalues,
		//		V,
		//		this->H_lanczos
		//	);
		//};
		//if (this->use_krylov) {
		//	try_alloc(this->build_krylov(););
		//	arma::cx_mat V; 
		//	arma::eig_sym(
		//		this->eigenvalues,
		//		V,
		//		this->H_lanczos
		//	);
		//	this->eigenvectors = this->krylov_space * V;
		//} else {
		//	try_alloc(this->build_lanczos(););
		//	//eigs(this->eigenvectors);
		//	arma::eig_sym(
		//		this->eigenvalues,
		//		this->eigenvectors,
		//		this->H_lanczos
		//	);
		//}
	}

}