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
	std::string dir = "results/disorder/OBC/";
	int L = 8;
	double g = 1.0;
	double g0 = 0;
	double h = 1.0;
	double w = 1.0;
	double J = 1;
	double J0 = 0.2;

	int realisations = 3;
	int bucket_num = 16;
	int mu = 5;
	int site = 0;


	std::ofstream file(dir + "SpectrumRapScaling_J0=" + to_string_prec(J0, 2) + \
		",g=" + to_string_prec(g, 2) + \
		",g0=" + to_string_prec(g0, 2) + \
		",h=" + to_string_prec(h, 2) + \
		",w=" + to_string_prec(w, 2) + ".dat");
	

	std::unique_ptr<IsingModel> Hamil = std::make_unique<IsingModel_sym>(L, J, g, h);
	/*for (L = 6; L <= 14; L++) {
		realisations = 900 - 100 * (L - 6);
		std::unique_ptr<IsingModel> Hamil = std::make_unique<IsingModel_disorder>(L, J, J0, g, g0, h, w);
		u64 N = Hamil->get_hilbert_size();
		arma::vec data(N, arma::fill::zeros);
		vec E(N);
		for (int r = 0; r < realisations; r++) {
			Hamil->hamiltonian();
			Hamil->diagonalization();
			data += Hamil->operator_av_in_eigenstates_return(&IsingModel::av_sigma_z, *Hamil, site);
		}
		data /= double(realisations);
		E = Hamil->get_eigenvalues();

		arma::vec fluct = data_fluctuations(data, mu);
		std::string name = "ProbDistSpectrumRap" + Hamil->get_info();
		probability_distribution(dir, name, data, -0.1, 0.1, 0.001);

		name = dir + "sigma_x" + Hamil->get_info() + ".dat";
		std::ofstream file2(name);
		for (int k = 0; k < data.size(); k++) {
			file2 << E(k) / double(L) << "\t" << data(k) << std::endl;
		}
		file2.close();

		file << N << "\t" << statistics_average(data, 10).t();
	}
	file.close();*/
	
	/*std::ofstream file("results/disorder/ipr_scaling," + Hamil->get_info() + ".dat");
	for (L = 6; L < 14; L++) {
		std::unique_ptr<IsingModel> Hamil = std::make_unique<IsingModel_disorder>(L, J, J0, g, g0, h, w);
		//std::unique_ptr<IsingModel> Hamil = std::make_unique<IsingModel_sym>(L, J, g, h);
		double ipr = 0;
		u64 N = Hamil->get_hilbert_size();
		mu = 0.25 * N;
		for (int r = 0; r < realisations; r++) {
			Hamil->hamiltonian();
			Hamil->diagonalization();
			for (int k = N / 2. - mu; k < N / 2. + mu; k++) {
				ipr += Hamil->ipr(k);
			}
		}
		file << L << "\t" << ipr / double(realisations * 2 * mu) / double(N) << std::endl;
	}
	file.close();*/
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

	auto stop1 = std::chrono::high_resolution_clock::now();
	out << "Time to finish simulation with symmetry: "<<\
		double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(stop1 - start)).count()) / 1000.0 << " seconds" << endl;
	
	auto stop2 = std::chrono::high_resolution_clock::now();
	out << double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(stop2 - stop1)).count()) / 1000.0 << " seconds" << endl;
	return 0;
}

/*for (L = 6; L <= 14; L++) {
		std::unique_ptr<IsingModel> Hamil = std::make_unique<IsingModel_disorder>(L, J, J0, g, g0, h, w);
		u64 N = Hamil->get_hilbert_size();
		Hamil->diagonalization();
		out << L << "\t";
		for (int o = N / 2 - mu; o < N / 2 + mu; o++) {
			const vec state = Hamil->get_eigenvectors().col(o);
			double sum = 0;
#pragma omp parallel for
			for (int k = 0; k < N; k++)
				sum += state(k) * exp(-k * k / 2.);
			out << sum << ", ";
		}
		out << std::endl;
	}*/