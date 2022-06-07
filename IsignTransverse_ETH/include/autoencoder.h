//-------------------------------------------------------------------- AUTO-ENCODER
void isingUI::ui::saveDataForAutoEncoder_symmetries(std::initializer_list<op_type> operators, std::initializer_list<std::string> names)
{
	using namespace std::chrono;
	clk::time_point start = std::chrono::system_clock::now();

	stout << "making symmetric model\n";
	auto params = arma::linspace(this->h, this->hs, this->hn);
	for (int system_size = this->L; system_size < this->L + this->Ls * this->Ln; system_size += this->Ls)
	{

		std::ofstream coeffLog;
		openFile(coeffLog, this->saving_dir + "coeffLog" + std::to_string(system_size) + ".dat", ios::out | ios::app);
		for (auto &hx : params)
		{
			stout << "\n\t\t\tSYM : h = " << hx << "\t\t\n";

			const auto start_loop = std::chrono::system_clock::now();
			auto alfa = std::make_unique<IsingModel_sym>(system_size, this->J, this->g, hx,
														 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
			const std::string saving_folder = this->saving_dir + alfa->get_info() + kPSep;
			const std::string saving_folder_wavefun = saving_folder + "wavefunctions" + kPSep;
			fs::create_directories(saving_folder);
			fs::create_directories(saving_folder_wavefun);

			stout << " \t\t--> finished creating model for " << alfa->get_info() << " - in time : "
				  << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
			alfa->diagonalization();
			stout << " \t\t--> finished diagonalizing for " << alfa->get_info() << " - in time : "
				  << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
			const u64 N = alfa->get_hilbert_size();
			stout << " \t\t--> The system size is " << N << "\n";

			std::ofstream wavefunLog;
			openFile(wavefunLog, saving_folder + "wavefun_log.dat", ios::out | ios::app);
			printSeparated(wavefunLog, "\t", 10, true, "filenum", "sigma_x(0)");

			auto ipr = 0.0;
			auto r = 0.0;
			double entropy = 0;
			this->mu = long(0.5 * N);
			long int E_min = alfa->E_av_idx - long(mu / 2);
			long int E_max = alfa->E_av_idx + long(mu / 2);
			int counter = 0;
			for (long int i = E_min; i < E_max; i++)
			{
				ipr += alfa->ipr(i);
				entropy += alfa->information_entropy(i);
				NO_OVERFLOW(r += alfa->eigenlevel_statistics(i, i + 1));
				counter++;
			}
			printSeparated(coeffLog, "\t", 10, true, to_string_prec(this->g), to_string_prec(hx), to_string_prec(ipr / double(N * counter), 8),
						   to_string_prec(r / double(counter), 8), to_string_prec(entropy / double(counter), 8));

			this->mu = long(0.3 * N);
			E_min = alfa->E_av_idx - long(mu / 2);
			E_max = alfa->E_av_idx + long(mu / 2);
			// let's go over that stuff
			int w_c = 0;
			for (long int i = E_min; i < E_max; i++)
			{
				// check sigma_x
				std::ofstream wavefun;
				openFile(wavefun, saving_folder_wavefun + to_string_prec(hx) + "_" + std::to_string(w_c) + "_wavefun_" + alfa->get_info() + ".dat", ios::out);
				const auto sigma_x = alfa->av_sigma_x(i, i, std::vector<int>({0}));
				printSeparated(wavefunLog, "\t", 10, true, to_string_prec(hx) + "_" + std::to_string(w_c), to_string_prec(sigma_x, 8));
				for (u64 j = 0; j < N; j++)
				{
					wavefun << alfa->get_eigenStateValue(i, j) << std::endl;
				}
				wavefun.close();
				w_c += 1;
			}
			wavefunLog.close();
		}
		coeffLog.close();
	}
}
void isingUI::ui::saveDataForAutoEncoder_disorder(std::initializer_list<op_type> operators, std::initializer_list<std::string> names)
{
	if (operators.size() != names.size())
		assert(false && "Set same number of input names as input operators!\n");
	using namespace std::chrono;
	clk::time_point start = std::chrono::system_clock::now();
	// diorder :3
	const int realisations = 50;
	const double delta = 0.025 * this->L; // width of offdiagonal in taking off-diagonal matrix elements -- to be updated maybe
	const int M = 500;					  // number of states taken across the offdiagonal -- for now chosen randomly
	const double w_end = this->w + this->wn * this->ws;
	const double g0_end = this->g0 + this->g0n * this->g0s;
	std::ofstream map;
	openFile(map, this->saving_dir + "map.dat", ios::out);
	printSeparated(map, "\t", 10, true, "real", "w", "g0", "c", "r", "rvar", "g", "h");
	for (double my_g0 = this->g0; my_g0 < g0_end; my_g0 += this->g0s)
	{
		for (double my_w = this->w; my_w < w_end; my_w += this->ws)
		{
			// reset generator
			gen = std::mt19937_64(seed);
			auto Hamil = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, my_g0, this->h, my_w, this->boundary_conditions);
			const u64 N = Hamil->get_hilbert_size();
			stout << "\n\n------------------------------ Doing : " << Hamil->get_info() << "------------------------------\n";

			// create main folder
			const std::string saving_folder = this->saving_dir + Hamil->get_info() + kPSep;

			// make folders for each operator separetely
			std::vector<std::string> opDirDiag, opDirNonDiag;
			std::string saving_folder_wave = saving_folder + "wavefunctions" + kPSep;
			createDirs(saving_folder, saving_folder_wave);

			for (auto &opName : names)
			{
				std::string saving_folder_operator = saving_folder + opName;
				std::string saving_folder_nondiag = saving_folder_operator + kPSep + "NonDiagonalMatrixElements" + kPSep;
				std::string saving_folder_diag = saving_folder_operator + kPSep + "DiagonalMatrixElements" + kPSep;
				createDirs(saving_folder_operator, saving_folder_diag, saving_folder_nondiag);
				opDirDiag.push_back(saving_folder_nondiag);
				opDirNonDiag.push_back(saving_folder_nondiag);
			}

			for (int realis = 0; realis < realisations; realis++)
			{
				Hamil->hamiltonian(); // restart Hamiltonian for new values of disorder
				Hamil->diagonalization();
				stout << " \t\t--> finished diagonalizing for " << Hamil->get_info() << " - in time : " << tim_s(start) << "s. Realisation -> " << realis << "\n";

				// set states from the middle of the spectrum
				this->mu = long(0.5 * N); // to include small system not not exceed Hilbert space
				long int E_min = Hamil->E_av_idx - long(mu / 2);
				long int E_max = Hamil->E_av_idx + long(mu / 2);

				// calculate level statistics
				double r = 0;
				double r_var = 0;
				for (int k = E_min; k < E_max; k++)
				{
					NO_OVERFLOW(double level_stat = Hamil->eigenlevel_statistics(k, k + 1));
					r += level_stat;
					r_var += level_stat * level_stat;
				}
				NO_OVERFLOW(
					r /= -double(E_min - E_max);	 // 1st moment - mean
					r_var /= -double(E_min - E_max); // 2nd moment
				);
				r_var = r_var - r * r; // variance

				// set new from the middle of the spectrum for operators
				this->mu = (M > N) ? long(0.5 * N) : M / 2; // to include small system not not exceed Hilbert space
				E_min = Hamil->E_av_idx - long(mu / 2);
				E_max = Hamil->E_av_idx + long(mu / 2);
				arma::mat H_diag = arma::diagmat(arma::Mat<double>(Hamil->get_hamiltonian()));
				arma::mat H_offdiag = Hamil->get_hamiltonian() - H_diag;
				auto c = matrixVariance(H_diag) / matrixVariance(H_offdiag);
				c = 1. / c;
				// iterate over input lists
				for (int a = 0; a < operators.size(); a++)
				{
					// assign by iterator
					op_type op = *(operators.begin() + a);
					std::string opName = *(names.begin() + a);

					// make file for log
					std::ofstream MatElemDiag;
					std::ofstream MatElemLogNonDiag;
					// save to wavefunctions log
					openFile(MatElemDiag, opDirDiag[a] + "MatrixElements_" + std::to_string(realis) + ".dat", ios::out);
					openFile(MatElemLogNonDiag, opDirNonDiag[a] + "MatrixElements_" + std::to_string(realis) + ".dat", ios::out);

					printSeparated(MatElemDiag, "\t", 10, false, "stateNum");
					printSeparated(MatElemLogNonDiag, "\t", 10, false, "<i|j>");
					for (int i = 0; i < this->L; i++)
					{
						printSeparated(MatElemDiag, "\t", 10, false, opName, "(", i, ")");
						printSeparated(MatElemLogNonDiag, "\t", 10, false, opName, "(", i, ")");
					}
					printSeparated(MatElemDiag, "\t", 10, true, "E_i");
					printSeparated(MatElemLogNonDiag, "\t", 10, true, "E_i - E_j");

					stout << "\n\n\t\t\t------------------------------ Starting operator " + opName + " for: " << Hamil->get_info() << "------------------------------\n";
					// go through the eigenstates
					for (u64 k = E_min; k < E_max; k++)
					{
						std::string wavename = saving_folder_wave +
											   std::to_string(k) + "," + std::to_string(realis) +
											   "," + to_string_prec(my_w, 2) +
											   +"," + to_string_prec(my_g0, 2) + "," + to_string_prec(c, 4);
						std::ofstream wavefunctionsLog;
						openFile(wavefunctionsLog, wavename + ".dat", ios::out);
						const int idx = long(k) - long(E_min);
						// check sigma_x
						// print k state
						printSeparated(MatElemDiag, "\t", 6, false, k);
						for (int i = 0; i < this->L; i++)
						{
							const auto opElem = Hamil->av_operator(k, k, op, std::vector<int>({i}));
							// const auto opElem = Hamil->av_op
							printSeparated(MatElemDiag, "\t", 10, false, to_string_prec(opElem, 8));
						}
						printSeparated(MatElemDiag, "\t", 10, true, to_string_prec(Hamil->get_eigenEnergy(k), 8));

						// give nondiagonal elements
						long int k2 = long(N) - long(k);
						printSeparated(MatElemLogNonDiag, "\t", 10, false, "<", k, "|", k2, ">");
						for (int i = 0; i < this->L; i++)
						{
							const auto opElem = Hamil->av_operator(k, k2, op, std::vector<int>({i})).real();
							printSeparated(MatElemLogNonDiag, "\t", 10, false, to_string_prec(opElem, 8));
						}
						printSeparated(MatElemLogNonDiag, "\t", 10, true, to_string_prec(Hamil->get_eigenEnergy(k) - Hamil->get_eigenEnergy(k2), 8));
						wavefunctionsLog << Hamil->get_eigenState(k);
						wavefunctionsLog.close();
					}
					MatElemDiag.close();
					MatElemLogNonDiag.close();
				}
				printSeparated(map, "\t", 10, true, realis, my_w, my_g0, c, r, r_var, this->g, this->h);
			}
		}
	}
	map.close();
}

