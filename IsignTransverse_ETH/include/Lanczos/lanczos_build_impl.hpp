#pragma once
#ifndef _LANCZOSBUILD
#define _LANCZOSBUILD

 
namespace lanczos 
{

	//<! builds lanczos tridiagonal matrix with or without
	//<! orthogonalization and no krylov space in memory
	inline
	void Lanczos::build_lanczos()
	{
		this->randVec_inKrylovSpace = arma::cx_vec(
			params.lanczos_steps,
			arma::fill::zeros
			);
		this->H_lanczos = arma::cx_mat(
			params.lanczos_steps,
			params.lanczos_steps,
			arma::fill::zeros
			);

		//<! set intial steps
		const u64 N = H.n_cols; //<! dimension of Hilbert space
		randVec_inKrylovSpace(0) = arma::cdot(this->initial_random_vec, this->initial_random_vec); // =1

		arma::cx_vec fi_next(N, arma::fill::zeros);
		//if (this->myParams.memory_over_performance)
		//	this->model->hamil_vec_onthefly(random_vec, fi_next);
		//else
		fi_next = H * this->initial_random_vec;

		arma::cx_vec fi_prev = this->initial_random_vec;
		cpx alfa = arma::cdot(this->initial_random_vec, fi_next);
		fi_next = fi_next - alfa * this->initial_random_vec;
		H_lanczos(0, 0) = alfa;

		//<! lanczos procedure
		for (int j = 1; j < params.lanczos_steps; j++) {
			cpx beta = arma::norm(fi_next);
			arma::cx_vec fi = fi_next / beta;
			randVec_inKrylovSpace(j) = arma::cdot(fi, this->initial_random_vec);

			//if (this->myParams.memory_over_performance)
			//	this->model->hamil_vec_onthefly(fi, fi_next);
			//else
			fi_next = H * fi;


			alfa = arma::cdot(fi, fi_next);
			fi_next = fi_next - alfa * fi - beta * fi_prev;

			H_lanczos(j, j) = alfa;
			H_lanczos(j, j - 1) = beta;
			H_lanczos(j - 1, j) = conj(beta);

			fi_prev = fi;
		}
		randVec_inKrylovSpace = arma::normalise(randVec_inKrylovSpace);
	}

	//<! builds lanczos tridiagonal matrix
	//<! with orthogonalization and krylov space
	inline 
	void Lanczos::build_krylov()
	{
		this->krylov_space = arma::cx_mat(
			this->N,
			this->params.lanczos_steps,
			arma::fill::zeros
			);
		this->H_lanczos = arma::cx_mat(
			this->params.lanczos_steps,
			this->params.lanczos_steps,
			arma::fill::zeros
			);

		this->krylov_space.col(0) = this->initial_random_vec;
		arma::cx_vec fi_next = this->H * krylov_space.col(0);

		cpx alfa = arma::cdot(this->krylov_space.col(0), fi_next);
		fi_next = fi_next - alfa * this->krylov_space.col(0);
		H_lanczos(0, 0) = alfa;

		for (int j = 1; j < this->params.lanczos_steps; j++) {
			cpx beta = arma::norm(fi_next);
			this->krylov_space.col(j) = fi_next / beta;

			fi_next = this->H * this->krylov_space.col(j);

			alfa = arma::cdot(this->krylov_space.col(j), fi_next);
			this->orthogonalize(fi_next, j);

			this->H_lanczos(j, j) = alfa;
			this->H_lanczos(j, j - 1) = beta;
			this->H_lanczos(j - 1, j) = conj(beta);
		}
		this->randVec_inKrylovSpace = this->krylov_space.t() * this->initial_random_vec;
	}

	//<! general lanczos build selecting either memory efficient or with krylov space
	inline
	void Lanczos::build(const arma::cx_vec& random) {
		this->initial_random_vec = random;
		this->build();
	}
	inline
	void Lanczos::build() {
		if (this->use_krylov)
			this->build_krylov();
		else
			this->build_lanczos();
	}

	//<! orthogonalizes input vector to the krylov space,
	//<! spanned by the first j vectors
	inline
	void Lanczos::orthogonalize(
			arma::cx_vec& vec_to_ortho,			//<! vector to orthogonalize
			int j									//<! current dimension of Krylov space
		) {
		arma::cx_vec temporary(this->N, arma::fill::zeros);
		for (int k = 0; k <= j; k++)
			temporary += arma::cdot(this->krylov_space.col(k), vec_to_ortho) * this->krylov_space.col(k);

		vec_to_ortho = vec_to_ortho - temporary;

		temporary.zeros();
		for (int k = 0; k <= j; k++)
			temporary += arma::cdot(this->krylov_space.col(k), vec_to_ortho) * this->krylov_space.col(k);

		vec_to_ortho = vec_to_ortho - temporary;
	};

	

}

#endif