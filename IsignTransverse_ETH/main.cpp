#include "include/IsingModel.h"


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
	auto start = std::chrono::high_resolution_clock::now();

	int L = 8;
	double g = 2.0;
	double h = 0.0;
	double w = 0.0;
	std::vector<double> J(L, 1);

	int realisations = 1000;
	int bucket_num = 10;

	std::unique_ptr<IsingModel> Hamil(new IsingModel_disorder(L, J, g, h, w));
	Hamil->diagonalization();
	u64 N = Hamil->get_hilbert_size();
	double dw = 0.02;

	std::ofstream fidel("results/Fidelity_L=" + std::to_string(L) + "_g=" + to_string_prec(g, 2)\
		+ "_h=" + to_string_prec(h, 2) + "_BC=" + std::to_string(_BC) + ".dat");
	std::ofstream level("results/level_stat_L=" + std::to_string(L) + "_g=" + to_string_prec(g, 2)\
		+ "_h=" + to_string_prec(h, 2) + "_BC=" + std::to_string(_BC) + ".dat");
	std::ofstream file_entropy("results/entropy_L=" + std::to_string(L) + "_g=" + to_string_prec(g, 2)\
		+ "_h=" + to_string_prec(h, 2) + "_BC=" + std::to_string(_BC) + ".dat");
	std::ofstream file_ipr("results/ipr_L=" + std::to_string(L) + "_g=" + to_string_prec(g, 2)\
		+ "_h=" + to_string_prec(h, 2) + "_BC=" + std::to_string(_BC) + ".dat");

	for (w = dw; w <= 5.0; w+=dw) {
		std::string name = "results/sigma_x_L=" + std::to_string(L) + "_g=" + to_string_prec(g, 2)\
			+ "_h=" + to_string_prec(h, 2) + "_w=" + to_string_prec(w, 2) + "_BC=" + std::to_string(_BC) + ".dat";
		vec fidelity(bucket_num - 1, fill::zeros);		
		vec level_stat(bucket_num - 1, fill::zeros);
		vec entropy(bucket_num - 1, fill::zeros);
		vec ipr(bucket_num - 1, fill::zeros);
		for (int k = 0; k < bucket_num - 1; k++) {
			u64 bucket_left = k * N / double(bucket_num);
			u64 bucket_right = (k + 1) * N / double(bucket_num);
			for (int r = 0; r < realisations; r++) {
				Hamil.reset(new IsingModel_disorder(L, J, g, h, w));
				Hamil->diagonalization();
				fidelity(k) += quantum_fidelity(bucket_left, bucket_right, Hamil, J, g, h, w + 1e-4) / double(realisations);
				level_stat(k) += Hamil->eigenlevel_statistics(bucket_left + 1, bucket_right) / double(realisations);
				for (u64 n = bucket_left; n < bucket_right; n++) {
					entropy(k) += Hamil->entaglement_entropy(n, L / 2) / double(bucket_right - bucket_left) / double(realisations);
					ipr(k) += Hamil->ipr(n) / double(bucket_right - bucket_left) / double(realisations);
				}
			}
		}
		for (int k = 0; k < bucket_num - 1; k++) {
			fidel << w << "\t\t" << k / double(bucket_num) << "\t\t" << fidelity(k) << endl;
			level << w << "\t\t" << k / double(bucket_num) << "\t\t" << level_stat(k) << endl;
			file_entropy << w << "\t\t" << k / double(bucket_num) << "\t\t" << entropy(k) << endl;
			file_ipr << w << "\t\t" << k / double(bucket_num) << "\t\t" << ipr(k) << endl;
		}
		Hamil->operator_av_in_eigenstates(&IsingModel::av_sigma_x, *Hamil, 0, name, "\t\t");
	}


	auto stop1 = std::chrono::high_resolution_clock::now();
	out << "Time to finish simulation with symmetry: "<<\
		double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(stop1 - start)).count()) / 1000.0 << " seconds" << endl;
	
	auto stop2 = std::chrono::high_resolution_clock::now();
	out << double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(stop2 - stop1)).count()) / 1000.0 << " seconds" << endl;
	return 0;
}

