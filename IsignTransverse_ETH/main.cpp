#include "include/user_interface.h"


std::random_device rd;
std::mt19937::result_type seed = rd() ^ (
	(std::mt19937::result_type)
	std::chrono::duration_cast<std::chrono::seconds>(
		std::chrono::system_clock::now().time_since_epoch()
		).count() +
	(std::mt19937::result_type)
	std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::high_resolution_clock::now().time_since_epoch()
		).count());
std::mt19937_64 gen(seed);


int main(const int argc, char* argv[]) {
	std::unique_ptr<user_interface> intface = std::make_unique<isingUI::ui>(argc, argv);
	intface->make_sim();

	
	/*std::ofstream file("results/spec_rep_map_dis.dat");
	for (g = -2.0; g <= 2.0; g += 0.1) {
		for (w = -3.0; w <= 3.0; w += 0.1) {
			std::unique_ptr<IsingModel> Hamil(new IsingModel_disorder(L, J, J0, g, g0, h, w));
			Hamil->diagonalization();
			double av_sigma = 0;
			for (int r = 0; r < realisations; r++) {
				av_sigma += Hamil->spectrum_repulsion(&IsingModel::av_sigma_z, *Hamil, 0) / double(realisations);
			}
			file << g << "\t\t" << w << "\t\t" << av_sigma << std::endl;
		}
	}
	file.close();
	*/

	/*
	std::unique_ptr<IsingModel> Hamil(new IsingModel_sym(L, J, g, h));
	Hamil->diagonalization();
	u64 N = Hamil->get_hilbert_size();
	out << "dim = " << N << std::endl;
	double dw = 0.02;

	std::ofstream fidel("results/Fidelity_L=" + std::to_string(L) + "_g=" + to_string_prec(g, 2)\
		+ "_h=" + to_string_prec(h, 2) + "_BC=" + std::to_string(_BC) + ".dat");
	std::ofstream level("results/level_stat_L=" + std::to_string(L) + "_g=" + to_string_prec(g, 2)\
		+ "_h=" + to_string_prec(h, 2) + "_BC=" + std::to_string(_BC) + ".dat");
	//std::ofstream file_entropy("results/entropy_L=" + std::to_string(L) + "_g=" + to_string_prec(g, 2)\
		+ "_h=" + to_string_prec(h, 2) + "_BC=" + std::to_string(_BC) + ".dat");
	std::ofstream file_ipr("results/ipr_L=" + std::to_string(L) + "_g=" + to_string_prec(g, 2)\
		+ "_h=" + to_string_prec(h, 2) + "_BC=" + std::to_string(_BC) + ".dat");

	for (g = dw; g <= 5.0; g += dw) {
		vec fidelity(bucket_num - 1, fill::zeros);
		vec level_stat(bucket_num - 1, fill::zeros);
		vec entropy(bucket_num - 1, fill::zeros);
		vec ipr(bucket_num - 1, fill::zeros);
		for (int k = 0; k < bucket_num - 1; k++) {
			u64 bucket_left = k * N / double(bucket_num);
			u64 bucket_right = (k + 1) * N / double(bucket_num);
			for (int r = 0; r < realisations; r++) {
				Hamil.reset(new IsingModel_sym(L, J, g, h));
				Hamil->diagonalization();
				fidelity(k) += quantum_fidelity(bucket_left, bucket_right, Hamil, J, g + 1e-4, h, w + 1e-4) / double(realisations);
				level_stat(k) += Hamil->eigenlevel_statistics(bucket_left + 1, bucket_right) / double(realisations);
				for (u64 n = bucket_left; n < bucket_right; n++) {
					//entropy(k) += Hamil->entaglement_entropy(n, L / 2) / double(bucket_right - bucket_left) / double(realisations);
					ipr(k) += Hamil->ipr(n) / double(bucket_right - bucket_left) / double(realisations) / double(N);
				}
			}
		}
		for (int k = 0; k < bucket_num - 1; k++) {
			fidel << g << "\t\t" << k / double(bucket_num - 1) << "\t\t" << fidelity(k) << endl;
			level << g << "\t\t" << k / double(bucket_num - 1) << "\t\t" << level_stat(k) << endl;
			//file_entropy << w << "\t\t" << k / double(bucket_num - 1) << "\t\t" << entropy(k) << endl;
			file_ipr << g << "\t\t" << k / double(bucket_num - 1) << "\t\t" << ipr(k) << endl;
		}
		std::string name = "results/sigma_x_L=" + std::to_string(L) + "_g=" + to_string_prec(g, 2)\
			+ "_h=" + to_string_prec(h, 2) + "_w=" + to_string_prec(w, 2) + "_BC=" + std::to_string(_BC) + ".dat";
		Hamil->operator_av_in_eigenstates(&IsingModel::av_sigma_x, *Hamil, 0, name, "\t\t");
	}
	*/
	return 0;
}