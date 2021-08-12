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

	int L = 6;
	double g = 0.5;
	double h = 0.0;
	double w = 0.0;

	std::vector<double> J(L, 1);
		std::unique_ptr<IsingModel> Hamil(new IsingModel_disorder(L, J, g, h, w));
		Hamil->diagonalization();
		u64 N = Hamil->get_hilbert_size();
		std::vector<double> vectttt = quantum_fidelity(0 * N / 20.0, 0 * N / 20.0+1, L, J, g, h, 0.02);
	mat corr_mat = Hamil->correlation_matrix(0);
	out << "<Sz(i)Sz(j)> = \n" << corr_mat << endl << endl;
	out << "S = " << Hamil->total_spin(corr_mat) << endl;
	out << "S^2 = " << arma::accu(corr_mat) << endl;
	vec Sq(L + 1, fill::zeros);
#pragma omp parallel for shared (Sq)
	for (int k = 0; k <= L; k++) {
		double q;
		if (!_BC) q = 2 * k * pi / (double)L;
		else q = k * pi / ((double)L + 1.0);
		cpx Sq_temp;
		for (int i = 0; i < L; i++) {
			for (int j = 0; j < L; j++) {
				Sq_temp += std::exp(1i * q * double(j - i)) * corr_mat(i, j);
			}
		}
		Sq(k) = real(2.0 * Sq_temp / pi / (L + 1.0));
	}
	out << Sq << endl;


	auto stop1 = std::chrono::high_resolution_clock::now();
	out << "Time to finish simulation with symmetry: "<<\
		double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(stop1 - start)).count()) / 1000.0 << " seconds" << endl;
	
	auto stop2 = std::chrono::high_resolution_clock::now();
	out << double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(stop2 - stop1)).count()) / 1000.0 << " seconds" << endl;
	return 0;
}



/*
	int bucket_num = L;
	int average_num = 1;

	ofstream file_Sx;
	ofstream colorMap_r;colorMap_r.open("results" + std::string(kPathSeparator) + "PhaseDiagram_L=" + to_string(L) + ",g=" + to_string_prec(g, 2) + ",h=" + to_string_prec(h, 2) + ".txt");
	if(!colorMap_r.is_open())
		throw "Not open color map file *sadface* \n";

	std::vector<double> J(L);
	std::fill(J.begin(), J.end(), 1.0);
	//for (h = 0.1; h <= 5.0; h += 0.1) {

		file_Sx.open("results" + std::string(kPathSeparator) + "Sx_L="+to_string(L)+",g="+to_string_prec(g,2)+",h="+to_string_prec(h,2)+",w=" + to_string_prec(disorder_strength,2)+ ".txt");
		if(!file_Sx.is_open())
			throw "Not open Sx file *sadface* \n";

		std::unique_ptr<IsingModel> B(new IsingModel_sym(L, J, g, h));
		u64 N = B->get_hilbert_size();
		vec r(bucket_num - 1,fill::zeros);								// save r for each bucket
		vec average(N,arma::fill::zeros);							// for sigma_x average
		for (int av = 0; av < average_num; av++) {
			B->hamiltonian();
			B->diagonalization();
			for(int bucket = 0; bucket < bucket_num - 1; bucket++){
				r(bucket) += B->eigenlevel_statistics(bucket*N/bucket_num + 1, (bucket + 1)* N / bucket_num - 1);
			}
			average += B->operator_av_in_eigenstates_return(&IsingModel::av_sigma_x,*B,1);
		}
		for(int state = 0; state < N; state++)
			file_Sx << B->get_eigenEnergy(state)/(double)L << "\t\t" << average(state) / (double)average_num << endl;
		//for(int bucket = 0; bucket < bucket_num - 1; bucket++)
			//colorMap_r << h << "\t\t" << (double)bucket/(double)bucket_num << "\t\t" << r(bucket)/(double)average_num << endl;
		out << h << "\t\t" << r.t();
		//colorMap_r << endl;
		file_Sx.close();*/
		//}
		//colorMap_r.close();

		//A->operator_av_in_eigenstates(&IsingModel::av_sigma_x, *A, 1, "results/sigma_x_average.txt", "\t\t");

		/*for (int L = 2; L <= 20; L += 2) {
			std::vector<double> J(L);
			std::fill(J.begin(), J.end(), 1.0);
			std::unique_ptr<IsingModel> A(new IsingModel_disorder(L, J, g, h));
			A->diagonalization();
			out << L << "\t\t" << A->spectrum_repulsion(&IsingModel::av_sigma_x, *A, 1) << endl;
		}*/