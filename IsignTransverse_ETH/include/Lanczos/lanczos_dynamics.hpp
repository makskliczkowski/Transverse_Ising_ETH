#pragma once

namespace lanczos{

	inline
	auto Lanczos::time_evolution_stationary(
		const arma::cx_vec& _state,	//<! initial state
		double time				//<! time step
	) -> arma::cx_vec
	{
		arma::cx_vec state_in_krylov = this->conv_to_krylov_space(_state);
		arma::cx_vec evolved_state(state_in_krylov.size(), arma::fill::zeros);
		for (int j = 0; j < this->params.lanczos_steps; j++) {
			cpx overlap = dot_prod(this->eigenvectors.col(j), state_in_krylov);
			evolved_state += std::exp(-im * this->eigenvalues(j) * time) * overlap * this->eigenvectors.col(j);
		}
		return arma::normalise(this->conv_to_hilbert_space(evolved_state));
	}

	inline
	auto Lanczos::time_evolution_non_stationary(
		arma::cx_vec& _state,	//<! state at time t: |c(t)>
		double dt,				//<! time step
		int lanczos_steps		//<! number of lanczos steps, here kept small
	)
	{ 
		const int M = this->params.lanczos_steps;
		this->params.lanczos_steps = lanczos_steps;
		this->diagonalization(_state);
		_state = this->time_evolution_stationary(_state, dt);  //<! locally approximated as stationary within (t, t+dt)
		this->params.lanczos_steps = M;
	}
}