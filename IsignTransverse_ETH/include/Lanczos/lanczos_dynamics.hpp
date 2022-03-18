#pragma once

namespace lanczos{

	inline
	auto Lanczos::time_evolution_stationary(
		arma::cx_vec& _state,	//<! initial state
		double time				//<! time step
	)
	{
		auto evolve = [this](
			arma::cx_vec& state,
			double time
			) 
		{
			arma::cx_vec temp_state(state.size(), arma::fill::zeros);
			for (int j = 0; j < this->params.lanczos_steps; j++) {
				cpx overlap = dot_prod(this->eigenvectors.col(j), state);
				temp_state += std::exp(-im * this->eigenvalues(j) * time) * overlap * this->eigenvectors.col(j);
			}
			state = temp_state;
		};
		if (false && this->use_krylov) {
			evolve(_state, time);
		}
		else {
			arma::cx_vec state_in_krylov = this->conv_to_krylov_space(_state);
			evolve(state_in_krylov, time);
			_state =  this->conv_to_hilbert_space(state_in_krylov);
		}
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
		this->time_evolution_stationary(_state, dt);  //<! locally approximated as stationary within (t, t+dt)
		this->params.lanczos_steps = M;
	}
}