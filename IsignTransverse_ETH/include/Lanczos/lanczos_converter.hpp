#pragma once

namespace lanczos {

	//! ------------------------------------------------------ from: KRYLOV -> to: HILBERT
	//<! conversion of input state to original Hilbert space
	template <typename _type>
	inline
	auto Lanczos<_type>::conv_to_hilbert_space(
			const arma::Col<_type>& rand,			//<! input random vector
			const arma::Col<_type>& state_lanczos	//<! index of state to transform
		) -> arma::Col<_type>
	{
		assert(state_lanczos.size() == this->params.lanczos_steps
			&& "Wrong state dimensions! Required dim = "
			+ std::to_string(this->params.lanczos_steps)
		);
		arma::Col<_type> state(this->N, arma::fill::zeros);	//<! output state

		if (this->use_krylov) {
			state = this->krylov_space * state_lanczos;
		} 
		else {
			arma::Col<_type> fi_next(this->N, arma::fill::zeros);
			//if (this->myParams.memory_over_performance)
			//	this->model->hamil_vec_onthefly(rand, fi_next);
			//else
			fi_next = this->H * rand;

			_type alfa = cdot(rand, fi_next);
			fi_next = fi_next - alfa * rand;
			arma::Col<_type> fi_prev = rand;
			for (int j = 1; j < this->params.lanczos_steps; j++) {
				_type beta = arma::norm(fi_next);
				arma::Col<_type> fi = fi_next / beta;

				state += state_lanczos(j) * fi;

				//if (this->myParams.memory_over_performance)
				//	this->hamil_vec_onthefly(fi, fi_next);
				//else
				fi_next = this->H * fi;

				alfa = arma::cdot(fi, fi_next);
				fi_next = fi_next - alfa * fi - beta * fi_prev;

				fi_prev = fi;
			}
		}
		return arma::normalise(state);
	};

	//<! conversion of eigenstate (by index) to original Hilbert space
	template <typename _type>
	inline
	auto
	Lanczos<_type>::conv_to_hilbert_space(
			const arma::Col<_type>& rand,		//<! input random vector
			int state_id						//<! index of state to transform
		) -> arma::Col<_type>
	{
		assert(!this->H_lanczos.is_empty() && "Diagonalize!!");
		return this->conv_to_hilbert_space(rand, this->eigenvectors.col(state_id));
	}

	//! ------------------------------------------------------ from: HILBERT -> to: KRYLOV
	template <typename _type>
	auto Lanczos<_type>::conv_to_krylov_space(
		const arma::Col<_type>& random_vec,	//<! input random vector
		const arma::Col<_type>& input		//<! state to transform from hilbert to krylov
	) -> arma::Col<_type>
	{
		assert(input.size() == this->params.lanczos_steps
			&& "Wrong state dimensions! Required dim = "
			+ std::to_string(this->N)
		);
		arma::Col<_type> transformed_input(
			params.lanczos_steps,
			arma::fill::zeros
			);
		if (this->use_krylov) {
			transformed_input = this->krylov_space.t() * input;
		}
		else {
			transformed_input(0) = arma::cdot(random_vec, input); // =1

			arma::Col<_type> fi_next = H * random_vec;
			arma::Col<_type> fi_prev = random_vec;

			_type alfa = arma::cdot(random_vec, fi_next);
			fi_next = fi_next - alfa * random_vec;
			H_lanczos(0, 0) = alfa;

			//<! lanczos procedure
			for (int j = 1; j < params.lanczos_steps; j++) {
				_type beta = arma::norm(fi_next);
				arma::Col<_type> fi = fi_next / beta;
				transformed_input(j) = arma::cdot(fi, input);
				fi_next = H * fi;

				alfa = arma::cdot(fi, fi_next);
				fi_next = fi_next - alfa * fi - beta * fi_prev;
				fi_prev = fi;
			}
		}
		return arma::normalise(transformed_input);
	}
}