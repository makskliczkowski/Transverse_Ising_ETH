#pragma once

namespace lanczos {

	//<! testing convergence of lanczos procedure 
	inline
	void Lanczos::convergence(
		int num_state_out				//<! number of states to test convergence
	) 
	{
		arma::vec e_prev(num_state_out, arma::fill::zeros);
		std::ofstream convergence("energy_convergence.txt");
		std::ofstream energy("energy_plot.txt");
		std::ofstream gs_error("state_convergence.txt");
		convergence << std::setprecision(16) << std::fixed;
		gs_error << std::setprecision(16) << std::fixed;
		const int M = this->params.lanczos_steps;
		for (int j = 1; j < M; j++) {
			this->params.lanczos_steps = j;

			this->diagonalization();

			for (int k = 0; k < j; k++) {
				energy << j << "\t" << this->eigenvalues(k) << std::endl;
			}
			if (j >= num_state_out) {
				arma::vec e(num_state_out, arma::fill::zeros);
				for (int id = 0; id < num_state_out; id++)
					e(id) = this->eigenvalues(id);
				convergence << j << "\t";
				gs_error << j << "\t";
				for (int k = 0; k < num_state_out; k++) {
					convergence << fabs((e_prev(k) - e(k)) / e(k)) << "\t";
					arma::cx_vec eigenState = conv_to_hilbert_space(k);
					double error = abs(arma::norm(this->eigenvalues(k) * eigenState - this->H * eigenState));
					gs_error << error << "\t";
				}
				convergence << std::endl;
				gs_error << std::endl;
				e_prev = e;
			}
			std::cout << j << std::endl;
		}
		energy.close();
		convergence.close();
		this->params.lanczos_steps = M;
	}


}