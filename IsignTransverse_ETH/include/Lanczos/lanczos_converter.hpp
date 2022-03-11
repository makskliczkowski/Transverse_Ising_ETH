#pragma once

namespace lanczos {
	//TODO: make separately from the class as standalone and do struct within class?
	//! ------------------------------------------------------ from: KRYLOV -> to: HILBERT
	//<! conversion of input state to original Hilbert space
	template <typename _type>
	template <typename _ty2>
	inline
	auto Lanczos<_type>::conv_to_hilbert_space(
			const arma::Col<std::complex<_ty2>>& state_lanczos	//<! state to transform
		) -> arma::Col<std::complex<_ty2>>
	{
		static_assert(std::is_convertible_v<_ty2, _type>, "No implicit type cast");
		assert(state_lanczos.size() == this->params.lanczos_steps
			&& "Wrong state dimensions! Required dim is the number of lanczos steps "
		);
		arma::Col<std::complex<_ty2>> state(this->N, arma::fill::zeros);	//<! output state

		if (this->use_krylov) {
			state = this->krylov_space * state_lanczos;
		} 
		else {
			arma::Col<_type> fi_next(this->N, arma::fill::zeros);
			//if (this->myParams.memory_over_performance)
			//	this->model->hamil_vec_onthefly(rand, fi_next);
			//else
			fi_next = this->H * this->initial_random_vec;

			_type alfa = cdot(this->initial_random_vec, fi_next);
			fi_next = fi_next - alfa * this->initial_random_vec;
			arma::Col<_type> fi_prev = this->initial_random_vec;
			for (int j = 1; j < this->params.lanczos_steps; j++) {
				_type beta = arma::norm(fi_next);
				arma::Col<_type> fi = fi_next / beta;

				state += state_lanczos(j) * fi;

				//if (this->myParams.memory_over_performance)
				//	this->hamil_vec_onthefly(fi, fi_next);
				//else
				fi_next = this->H * fi;

				alfa = dot_prod(fi, fi_next);
				fi_next = fi_next - alfa * fi - beta * fi_prev;

				fi_prev = fi;
			}
		}
		return arma::normalise(state);
	};
	template <typename _type>
	template <typename _ty2>
	inline
		auto Lanczos<_type>::conv_to_hilbert_space(
			const arma::Col<_ty2>& state_lanczos	//<! state to transform
		) -> arma::Col<_type>
	{
		static_assert(std::is_convertible_v<_ty2, double>, "Not convertible");
		auto state = cpx_real_vec(state_lanczos);
		return arma::real(this->conv_to_hilbert_space(state));
	}
	//<! conversion of eigenstate (by index) to original Hilbert space
	template <typename _type>
	inline
	auto
	Lanczos<_type>::conv_to_hilbert_space(
			int state_id						//<! index of state to transform
		) -> arma::Col<_type>
	{
		assert(!this->H_lanczos.is_empty() && "Diagonalize!!");
		return this->conv_to_hilbert_space(this->eigenvectors.col(state_id));
	}

	//! ------------------------------------------------------ from: HILBERT -> to: KRYLOV
	template <typename _type>
	template <typename _ty2>
	auto Lanczos<_type>::conv_to_krylov_space(
		const arma::Col<std::complex<_ty2>>& input		//<! state to transform from hilbert to krylov
	) -> arma::Col<std::complex<_ty2>>
	{
		static_assert(std::is_convertible_v<_ty2, _type>, "No implicit type cast");
		assert(input.size() == this->N
			&& "Wrong state dimensions! Required dim is the original hilbert space size "
		);
		arma::Col<std::complex<_ty2>> transformed_input(
			params.lanczos_steps,
			arma::fill::zeros
			);
		if (this->use_krylov) {
			transformed_input = this->krylov_space.t() * input;
		}
		else {
			transformed_input(0) = dot_prod(this->initial_random_vec, input); // =1

			arma::Col<_type> fi_next = H * this->initial_random_vec;
			arma::Col<_type> fi_prev = this->initial_random_vec;

			_type alfa = dot_prod(this->initial_random_vec, fi_next);
			fi_next = fi_next - alfa * this->initial_random_vec;
			H_lanczos(0, 0) = alfa;

			//<! lanczos procedure
			for (int j = 1; j < params.lanczos_steps; j++) {
				_type beta = arma::norm(fi_next);
				arma::Col<_type> fi = fi_next / beta;
				transformed_input(j) = dot_prod(fi, input);
				fi_next = H * fi;

				alfa = dot_prod(fi, fi_next);
				fi_next = fi_next - alfa * fi - beta * fi_prev;
				fi_prev = fi;
			}
		}
		return arma::normalise(transformed_input);
	}
	template <typename _type>
	template <typename _ty2>
	inline
		auto Lanczos<_type>::conv_to_krylov_space(
			const arma::Col<_ty2>& input	//<! state to transform
		) -> arma::Col<_type>
	{
		static_assert(std::is_convertible_v<_ty2, double>, "Not convertible");
		auto state = cpx_real_vec(input);
		return arma::real(this->conv_to_krylov_space(state));
	}
}