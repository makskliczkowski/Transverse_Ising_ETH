#pragma once

namespace lanczos{

	inline
	auto Lanczos::time_evolution_stationary(
		const arma::cx_vec& input_state,	//<! initial state
		double time							//<! time step
	) -> arma::cx_vec 
	{
		auto evolve = [this](
			const arma::cx_vec& state,
			arma::cx_vec& output,
			double time
			) {
			for (int j = 0; j < this->params.lanczos_steps; j++) {
				cpx overlap = dot_prod(this->eigenvectors.col(j), state);
				output += std::exp(-im * this->eigenvalues(j) * time) * overlap * this->eigenvectors.col(j);
			}
		};
		if (false && this->use_krylov) {
			arma::cx_vec output(this->N, arma::fill::zeros);
			evolve(input_state, output, time);
			return output;
		}
		else {
			arma::cx_vec output(this->params.lanczos_steps, arma::fill::zeros);
			evolve(this->conv_to_krylov_space(input_state), output, time);
			return this->conv_to_hilbert_space(output);
		}
	}

	inline
	auto Lanczos::time_evolution_non_stationary(
		const arma::cx_vec& prev_state,	//<! state at time t: |c(t)>
		double dt,						//<! time step
		int lanczos_steps				//<! number of lanczos steps, here kept small
	) -> arma::cx_vec					//<! output: state at time t+dt: c(t+dt)>
	{ 
		const int M = this->params.lanczos_steps;
		this->params.lanczos_steps = lanczos_steps;
		this->diagonalization(prev_state);
		auto ret = this->time_evolution_stationary(prev_state, dt); //<! locally approximated as stationary within (t, t+dt)
		this->params.lanczos_steps = M;
		return ret;
	}
}