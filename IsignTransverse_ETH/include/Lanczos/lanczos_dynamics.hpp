#pragma once

namespace lanczos{

	template <typename _type>
	inline
	auto Lanczos<_type>::time_evolution_stationary(
		const arma::cx_vec& input_state,	//<! initial state
		double time							//<! time step
	) -> arma::cx_vec 
	{
		auto state_lanczos = this->conv_to_krylov_space(input_state);
		arma::cx_vec output(this->N, arma::fill::zeros);
		for (int j = 0; j < this->params.lanczos_steps; j++) {
			_type overlap = arma::cdot(this->eigenvectors.col(j), input_state);
			output += std::exp(-im * this->eigenvalues(j) * time) * overlap * this->eigenvectors.col(j);
		}
		return this->conv_to_hilbert_space(output);
	}

	template <typename _type>
	inline
	auto Lanczos<_type>::time_evolution_non_stationary(
		const arma::cx_vec& prev_state,	//<! state at time t: |c(t)>
		double dt,						//<! time step
		int lanczos_steps				//<! number of lanczos steps, here kept small
	) -> arma::cx_vec					//<! output: state at time t+dt: c(t+dt)>
	{ 
		const int M = this->params.lanczos_steps;
		this->params.lanczos_steps = lanczos_steps;
		this->diagonalization();
		this->time_evolution_stationary(prev_state, dt); //<! locally approximated as stationary within (t, t+dt)
		this->params.lanczos_steps = M;
	}
}