#include "include/user_interface.h"
// set externs
std::uniform_real_distribution<> theta	= std::uniform_real_distribution<>(0.0, pi);
std::uniform_real_distribution<> fi		= std::uniform_real_distribution<>(0.0, pi);
int outer_threads = 1;
//---------------------------------------------------------------------------------------------------------------- UI main
void isingUI::ui::make_sim()
{
	printAllOptions();
	gen = std::mt19937_64(this->seed);

	clk::time_point start = std::chrono::system_clock::now();
	const int Lmin = this->L, Lmax = this->L + this->Ln * this->Ls;
	const double gmin = this->g, gmax = this->g + this->gn * this->gs;
	const double hmin = this->h, hmax = this->h + this->hn * this->hs;

	switch (this->fun)
	{
	case 0: 
		diagonalize(); 
		break;
	case 1:
		calculate_spectrals();
		break;
	case 2:
		entropy_evolution();
		break;
	case 3: 
		spectral_form_factor();
		break;
	case 4:		
		relaxationTimesFromFiles();
		break;
	case 5:
		benchmark();
		break;
	case 6:
		adiabatic_gauge_potential();
		break;
	case 7:
		level_spacing();
		break;
	case 8:
		thouless_times();
		break;
	default:
		printSeparated(std::cout, "\t", 16, true, "L", "J", "g", "h");
		for (int system_size = Lmin; system_size < Lmax; system_size += this->Ls)
		{
			for (double gx = gmin; gx < gmax; gx += this->gs)
			{
				for (double hx = hmin; hx < hmax; hx += this->hs)
				{
					this->L = system_size;
					this->g = gx;
					this->h = hx;
					const auto start_loop = std::chrono::system_clock::now();
//for(this->J = 0.05; this->J <= 1.5; this->J += 0.05)
{
	//if(this->L > 10) this->realisations = 1000;

					printSeparated(std::cout, "\t", 16, true, this->L, this->J, this->g, this->h);
	//auto model = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
	//model->diagonalization();
	//auto E = model->first_interacting_correction();
	//for(int i = 0; i < E.size(); i++)
	//	printSeparated(std::cout, "\t", 16, true, E(i), model->get_eigenEnergy(i), abs(E(i) - model->get_eigenEnergy(i)));
	//stout << "----------------------------------" << std::endl;
	//continue;
					// ----------------------
					//this->diagonalize(); continue;
					for(this->w = 1.0; this->w <= 2.5; this->w += 0.1)
					{
						std::cout << this->w << std::endl;
						arma::vec loc_length(this->L, arma::fill::zeros), energy;
					#pragma omp parallel for
						for(int r = 0; r < this->realisations; r++){
							auto [E, loc] = anderson::get_localisation_length(this->L, this->J, this->w);
							#pragma omp critical
							{
								energy = E;
								loc_length += loc;
							}
						}
						save_to_file(this->saving_dir + "LocLengthDist_" + to_string_prec(this->w, 2) + ".dat", energy, loc_length / double(this->realisations)); 
						
						//for(this->site = 0; this->site < this->L; this->site++)
						//	calculate_spectrals();
						//diagonalize();
						//spectral_form_factor();
						//analyze_spectra();
						average_SFF();
					}
					continue;
					average_SFF(); continue;
					//combine_spectra(); 
					analyze_spectra(); continue;
					thouless_times(); continue;
					spectral_form_factor(); continue;
					std::string info = IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w);
					smoothen_data(this->saving_dir + "SpectralFormFactor" + kPSep, info + ".dat"); continue;
} continue;
					std::string dir = this->saving_dir + "Entropy" + kPSep;// + "Lanczos" + kPSep;

					auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
					////auto alfa = std::make_unique<IsingModel_sym>(this->L, -this->J, -this->g, -this->h,
					////			 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
					const std::string name = alfa->get_info();
					const size_t N = alfa->get_hilbert_size();
					stout << "\n\t\t--> finished creating model for " << name << " - in time : " << tim_s(start) << "s" << std::endl;

					for(this->site = 0; this->site <= this->L / 2.; this->site++){

					auto [opName, dir_suffix] = IsingModel_disorder::opName(this->op, this->site);
					dir = this->saving_dir + "timeEvolution" + kPSep + dir_suffix + "=" + std::to_string(this->site) + kPSep;
					std::string dir2 = dir + "Exponent" + kPSep;
					//std::string dir2 = this->saving_dir + "TimeEvolution" + kPSep + "Quench" + kPSep;
					createDirs(dir, dir2);
					
					smoothen_data(dir, opName + name + ".dat");
					std::ifstream input;
					auto data = readFromFile(input, dir + "smoothed" + kPSep + opName + name + ".dat");
					if(data.empty()) continue;
					arma::vec exponent_exp 		 = non_uniform_derivative(data[0], data[1]);
					arma::vec exponent_power_law = log_derivative(data[0], data[1]);
					std::ofstream output;
					openFile(output, dir2 + opName + name + ".dat", std::ios::out);
					for(int j = 0; j < exponent_exp.size(); j++)
						printSeparated(output, "\t", 16, true, data[0](j), (data[1](j) - data[3](0)) / exponent_exp(j), exponent_power_law(j));
					output.close();
					
					}
					continue;
				}
			}
		}
	}
	stout << " - - - - - - FINISHED CALCULATIONS IN : " << tim_s(start) << " seconds - - - - - - " << std::endl; // simulation end
}


// ----------------------------------------------------------------------------- SIMULATIONS -----------------------------------------------------------------------------
//-------------------------------------------------------------------------- GENERAL ROUTINES
void isingUI::ui::diag_sparse(int num, bool get_eigenvectors, double sigma)
{
	auto start_up = std::chrono::system_clock::now();
	auto save_eigs = [this](
		const auto& eigenvalues,	//<! eigenvalues to save
		std::string info 			//<! info about model parameters
		){
		std::string dir = this->saving_dir + "EIGENVALUES" + kPSep;
		createDirs(dir);
		std::string name = dir + kPSep + info;
		eigenvalues.save(arma::hdf5_name(name + ".h5", "eigenvalues"));
	};
	if(this->m){
		clk::time_point start = std::chrono::system_clock::now();
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		std::string name = alfa->get_info();
		stout << "\n\t\t--> finished creating model for " << name << " - in time : " << tim_s(start) << "s" << std::endl; 
		start = std::chrono::system_clock::now();
		arma::eigs_opts opts;
		opts.maxiter = 5000;
		opts.tol = 1e-16;
		arma::sp_cx_mat H = alfa->get_hamiltonian();
		//arma::cx_vec E = arma::eigs_gen(H, num, sigma, opts);
		//save_eigs(E, name);
		stout << "\t\t	--> finished diagonalizing for " << alfa->get_info()
			  << " - in time : " << tim_s(start) << "\t\nTotal time : " << tim_s(start_up) << "s" << std::endl;
	} else{
		for(int r = 0; r < this->realisations; r++){
			auto start = std::chrono::system_clock::now();
			auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
			std::string name = alfa->get_info();
			stout << "\n\t\t--> finished creating model for " << name << " - in time : " << tim_s(start) << "s" << std::endl;
			start = std::chrono::system_clock::now();
			arma::eigs_opts opts;
			//opts.maxiter = 1000;
			opts.tol = 1e-12;
			arma::sp_mat H = alfa->get_hamiltonian();
			arma::vec E = arma::eigs_sym(H, num, sigma, opts);
			save_eigs(E, name + ((r == 0)? "" : "_real=" + std::to_string(r)));
			stout << "\t\t	--> finished diagonalizing for " << alfa->get_info()
				  << " - in time : " << tim_s(start) << "\t\nTotal time : " << tim_s(start_up) << "s" << std::endl;
		}
	}
}

void isingUI::ui::diagonalize(){
	clk::time_point start = std::chrono::system_clock::now();
	std::string dir = this->saving_dir + "DIAGONALIZATION" + kPSep;
	createDirs(dir);

	auto kernel = [this, &start, &dir](
		auto& alfa, 				//<! model
		std::string _suffix = ""	//<! suffix, like realisation etc
		){
		std::string info = alfa.get_info({});
		stout << "\n\t\t--> finished creating model for " << info + _suffix << " - in time : " << tim_s(start) << "s" << std::endl;
		arma::vec eigenvalues;
		if(alfa.g == 0){
			auto H = alfa.get_hamiltonian();
			const u64 N = alfa.get_hilbert_size();
			arma::cx_vec E(N);
			for(int i = 0; i < N; i++)
				E(i) = H(i,i);
			eigenvalues = real(E);
			sort(eigenvalues.begin(), eigenvalues.end());
		} 
		else if(alfa.J == 0.0){
			eigenvalues = alfa.get_non_interacting_energies();
		} else{
			alfa.diagonalization(this->ch);
			eigenvalues = alfa.get_eigenvalues();
		}
		stout << "\t\t	--> finished diagonalizing for " << info + _suffix<< " - in time : " << tim_s(start) << "s" << std::endl;
		
		std::string name = dir + info + _suffix + ".hdf5";
		eigenvalues.save(arma::hdf5_name(name, "eigenvalues/", arma::hdf5_opts::append));
		stout << "\t\t	--> finished saving eigenvalues for " << info + _suffix << " - in time : " << tim_s(start) << "s" << std::endl;
		//if(this->ch){
		//	auto V = alfa.get_eigenvectors();
		//	V.save(arma::hdf5_name(name, "/eigenvectors/" + _suffix, arma::hdf5_opts::append));
		//	stout << "\t\t	--> finished saving eigenvectors for " << info << " - in time : " << tim_s(start) << "s" << std::endl;
		//}
	};
	//herr_t status;
	// ----------- choose model and run kernel
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		kernel(*alfa);
	} else{
		auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
		std::string info = alfa->get_info({});
		for(int r = 0; r < this->realisations; r++){
			alfa->hamiltonian();
			kernel(*alfa, "_real=" + std::to_string(r + this->jobid));
		}
	}
}
template <typename _type>
auto isingUI::ui::get_eigenvalues(IsingModel<_type>& alfa, std::string _suffix) 
	-> arma::vec
{
	arma::vec eigenvalues;
	std::string dir = this->saving_dir + "DIAGONALIZATION" + kPSep;
	createDirs(dir);
	std::string name = dir + alfa.get_info({});
	bool loaded;
	#pragma omp critical
	{
		loaded = eigenvalues.load(arma::hdf5_name(name + _suffix + ".hdf5", "eigenvalues/"));
		if(!loaded)
			loaded = eigenvalues.load(arma::hdf5_name(name + ".hdf5", "eigenvalues/" + _suffix));
	}
	if(!loaded){
		if(alfa.g == 0){
			auto H = alfa.get_hamiltonian();
			const u64 N = alfa.get_hilbert_size();
			arma::cx_vec E(N);
			for(int i = 0; i < N; i++)
				E(i) = H(i,i);
			eigenvalues = real(E);
			sort(eigenvalues.begin(), eigenvalues.end());
		} 
		else if(alfa.J == 0.0){
			eigenvalues = alfa.get_non_interacting_energies();
		} 
		else {
			#undef MY_MAC
			#if defined(MY_MAC)
				std::cout << "Failed to load energies, returning empty array" << std::endl;
			#else
				alfa.diagonalization(false);
				//stout << "No energies found, diagonalizing matrix now!" << std::endl;
				eigenvalues = alfa.get_eigenvalues();
			#endif
		}
	}
	#ifndef MY_MAC
		// save eigenvalues (yet unsaved)
		if(!loaded)
			eigenvalues.save(arma::hdf5_name(name + _suffix + ".hdf5", "eigenvalues/"));
	#endif
	return eigenvalues;
}

void isingUI::ui::combine_spectra(){
	std::string dir = this->saving_dir + "DIAGONALIZATION" + kPSep;
	auto kernel = [this, &dir](int realis, double x){
		arma::vec eigenvalues;
		std::string _suffix = "_real=" + std::to_string(realis);
		std::string info = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"k", "x", "p"}) 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w);
		std::string name = dir + info + _suffix + ".hdf5";
		bool loaded = eigenvalues.load(arma::hdf5_name(name, "/eigenvalues/dataset"));
		if(loaded){
			#pragma omp critical
			{
				eigenvalues.save(arma::hdf5_name(dir + info + ".hdf5", "/eigenvalues/" + _suffix, arma::hdf5_opts::append));
			}
		}
	};
	outer_threads = 1;
	if(this->m){
		average_over_realisations<Ising_params::h>(false, kernel);
	} else{
		average_over_realisations<Ising_params::h>(false, kernel);
	}
}
//-------------------------------------------------------------------- COMPARING SYMMETRIC TO DISORDERED RESULTS
void isingUI::ui::compare_energies()
{
	// handle disorder
	auto Hamil = std::make_unique<IsingModel_disorder>(L, J, 0, g, 0, h, 0, boundary_conditions);
	Hamil->diagonalization();
	const arma::vec E_dis = Hamil->get_eigenvalues();
	Hamil.release();
	// here we push back the energies from symmetries
	std::vector<double> E_sym = v_1d<double>();
	std::vector<std::string> symms = v_1d<std::string>();
	// go for each symmetry sector
	for (int k = 0; k < L; k++)
	{
		if (k == 0 || k == this->L / 2.)
		{
			for (int p = 0; p <= 1; p++)
			{
				// if the spin flip is unaviable we just use 1
				const int x_max = (this->h != 0) ? 0 : 1;
				for (int x = 0; x <= x_max; x++)
				{
					auto ham = std::make_unique<IsingModel_sym>(L, J, g, h, k, p, x, boundary_conditions);
					ham->diagonalization();
					arma::vec t = ham->get_eigenvalues();
					E_sym.insert(E_sym.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));
					v_1d<std::string> temp_str = v_1d<std::string>(t.size(), "k=" + std::to_string(k) + ",x=" + to_string(x) + ",p=" + to_string(p));
					symms.insert(symms.end(), std::make_move_iterator(temp_str.begin()), std::make_move_iterator(temp_str.end()));
				}
			}
		}
		else
		{
			int x_max = (this->h != 0) ? 0 : 1;
			for (int x = 0; x <= x_max; x++)
			{
				auto ham = std::make_unique<IsingModel_sym>(L, J, g, h, k, 1, x, boundary_conditions);
				ham->diagonalization();
				arma::vec t = ham->get_eigenvalues();
				E_sym.insert(E_sym.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));
				v_1d<std::string> temp_str = v_1d<std::string>(t.size(), "k=" + std::to_string(k) + ",x=" + to_string(x) + ",p=*");
				symms.insert(symms.end(), std::make_move_iterator(temp_str.begin()), std::make_move_iterator(temp_str.end()));
			}
		}
		stout << k << endl;
	}
	// face permutation
	auto permut = sort_permutation(E_sym, [](const double a, const double b)
								   { return a < b; });
	apply_permutation(E_sym, permut);
	apply_permutation(symms, permut);
	stout << endl
		  << E_sym.size() << endl;
	stout << "symmetry sector\t\tEnergy sym\t\tEnergy total\t\tdifference\t\t <n|X|n>" << endl;
	for (int k = 0; k < E_dis.size(); k++)
	{
		stout << symms[k] << "\t\t" << E_sym[k] << "\t\t" << E_dis(k) << "\t\t" << E_sym[k] - E_dis(k) << endl;
	}
}
void isingUI::ui::compare_matrix_elements(op_type op, int k_alfa, int k_beta, int p_alfa, int p_beta, int x_alfa, int x_beta)
{
	string name = IsingModel_sym::set_info(this->L, this->J, this->g, this->h, k_alfa, p_alfa, x_alfa, {"k", "p", "x"});
	// name = name + "k_ab=" + to_string(k_alfa) + ":" + to_string(k_beta) + "p_ab=" + to_string(p_alfa) + ":" + to_string(p_beta) + "x_ab=" + to_string(x_alfa) + ":" + to_string(x_beta);
	std::ofstream file(this->saving_dir + "matrix_comparison" + name + ".dat");
	// disorder
	auto model = std::make_unique<IsingModel_disorder>(L, J, 0, g, 0, h, 0);
	model->diagonalization();
	arma::sp_cx_mat opMatrix = model->create_operator({op}, std::vector<int>({0}));
	file << "WHOLE SIZE = " << model->get_hilbert_size() << endl;
	// symmetries
	std::unique_ptr<IsingModel_sym> alfa = std::make_unique<IsingModel_sym>(L, J, g, h, k_alfa, p_alfa, x_alfa);
	std::unique_ptr<IsingModel_sym> beta = std::make_unique<IsingModel_sym>(L, J, g, h, k_beta, p_beta, x_beta);
	arma::sp_cx_mat opMatrix_alfa = alfa->create_operator({op}, std::vector<int>({0}));
	arma::sp_cx_mat opMatrix_beta = beta->create_operator({op}, std::vector<int>({0}));
	file << "ALFA SECTOR SIZE = " << alfa->get_hilbert_size() << endl;
	file << "BETA SECTOR SIZE = " << beta->get_hilbert_size() << endl;
	alfa->diagonalization();
	beta->diagonalization();

	// maps
	auto map_alfa = mapping_sym_to_original(0, model->get_hilbert_size() - 1, *alfa, *model);
	auto map_beta = mapping_sym_to_original(0, model->get_hilbert_size() - 1, *beta, *model);
	file << "alfa sector size for nondegenerate mapping is : " << map_alfa.size() << endl;
	file << "beta sector size for nondegenerate mapping is : " << map_beta.size() << endl;

	file << "ALPHA SECTOR : " + alfa->get_info({"L", "J", "g", "h"}) + "\n BETA SECTOR :" + beta->get_info({"L", "J", "g", "h"}) + "\n\n";
	file << " - - - - - - SAME SECTORS - - - - - - \n"
		 << " - - - - > FOR ALFA - ALFA: \n"
		 << "ENERGY ALPHA |('.'|)"
		 << "\t\t"
		 << "ENERGY ALFA (/'.')/"
		 << "\t\t"
		 << "<alfa|SIGMA_X|alfa>"
		 << "\t\t"
		 << "<non_sym|SIGMA_X|non_sym>"
		 << "\t\t"
		 << "DIFFERENCE\t\tA/B" << endl;
	for (auto &element : map_alfa)
	{
		for (auto &t : map_alfa)
		{
			cpx A = arma::cdot(alfa->get_eigenState(element.first), opMatrix_alfa * alfa->get_eigenState(t.first));
			cpx B = arma::dot(model->get_eigenState(element.second), opMatrix * model->get_eigenState(t.second));
			if (abs(abs(A) - abs(B)) >= 1e-14)
				file << alfa->get_eigenEnergy(element.first) << "\t\t\t\t\t" << alfa->get_eigenEnergy(t.first) << "\t\t\t\t" << A << "\t\t\t\t" << B << "\t\t\t\t" << abs(A) - abs(B) << "\t\t\t\t" << abs(A) / abs(B) << std::endl;
		}
	}
	file << "\n - - - - > FOR BETA - BETA: \n"
		 << "ENERGY BETA |('.'|)"
		 << "\t\t"
		 << "ENERGY BETA (/'.')/"
		 << "\t\t"
		 << "<beta|SIGMA_X|beta>"
		 << "\t\t"
		 << "<non_sym|SIGMA_X|non_sym>\t\tDIFFERENCE\t\tA/B" << endl;
	for (auto &element : map_beta)
	{
		for (auto &t : map_beta)
		{
			cpx A = arma::cdot(beta->get_eigenState(element.first), opMatrix_beta * beta->get_eigenState(t.first));
			cpx B = arma::dot(model->get_eigenState(element.second), opMatrix * model->get_eigenState(t.second));
			if (abs(abs(A) - abs(B)) >= 1e-14)
				file << beta->get_eigenEnergy(element.first) << "\t\t\t\t\t" << beta->get_eigenEnergy(t.first) << "\t\t\t\t" << A << "\t\t\t\t" << B << "\t\t\t\t" << abs(A) - abs(B) << "\t\t\t\t" << abs(A) / abs(B) << std::endl;
		}
	}
	file << "\n\n - - - - - - DIFFERENT SECTORS - - - - - - \n"
		 << " - - - - > FOR ALFA - BETA: \n"
		 << "ENERGY ALPHA |('.'|)"
		 << "\t\t"
		 << "ENERGY BETA (/'.')/"
		 << "\t\t"
		 << "<alfa|SIGMA_X|beta>"
		 << "\t\t"
		 << "<non_sym_a|SIGMA_X|non_sym_b>"
		 << "\t\t"
		 << "DIFFERENCE\t\tA/B" << endl;
	for (auto &element : map_alfa)
	{
		for (auto &t : map_beta)
		{
			// cpx A = arma::cdot(alfa->get_eigenState(element.first), opMatrix_beta * beta->get_eigenState(t.first));
			cpx A = 0.; // TODO: here add operator between sectors
			cpx B = arma::dot(model->get_eigenState(element.second), opMatrix * model->get_eigenState(t.second));
			if (abs(abs(A) - abs(B)) >= 1e-14)
				file << alfa->get_eigenEnergy(element.first) << "\t\t\t\t\t" << beta->get_eigenEnergy(t.first) << "\t\t\t\t" << A << "\t\t\t\t" << B << "\t\t\t\t" << abs(A) - abs(B) << "\t\t\t\t" << abs(A) / abs(B) << std::endl;
		}
	}
	file.flush();
	file.close();
}
void isingUI::ui::compare_entaglement()
{
	clk::time_point start = std::chrono::system_clock::now();
	auto alfa1 = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, 1e-4, this->boundary_conditions);
	auto alfa2 = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, 1e-3, this->boundary_conditions);
	auto alfa3 = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, 1e-2, this->boundary_conditions);
	stout << "\n\t\t--> finished creating models for w=1e-4,1e-3,1e-2 and  " << alfa1->get_info({"w"}) << " - in time : " << tim_s(start) << "s" << std::endl;
	alfa1->diagonalization();
	alfa2->diagonalization();
	alfa3->diagonalization();
	stout << "\t\t	--> finished diagonalizing for w=1e-4,1e-3,1e-2 and  " << alfa1->get_info({"w"}) << " - in time : " << tim_s(start) << "s" << std::endl;

	auto beta1 = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h, this->L / 2, 1, 1, this->boundary_conditions);
	auto beta2 = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h, 1, 1, 1, this->boundary_conditions);
	stout << "\n\t\t--> finished creating model for k=0,1; p=1,x=1 and " << beta1->get_info({"k", "p", "x"}) << " - in time : " << tim_s(start) << "s" << std::endl;
	beta1->diagonalization();
	beta2->diagonalization();
	stout << "\t\t	--> finished diagonalizing for k=0,1; p=1,x=1 and " << beta1->get_info({"k", "p", "x"}) << " - in time : " << tim_s(start) << "s" << std::endl;

	std::ofstream file;
	std::string dir = this->saving_dir + "Entropy" + kPSep;
	createDirs(dir);
	openFile(file, dir + "compare_to_disorder" + beta1->get_info({}) + ".dat");
	const u64 dim = alfa1->get_hilbert_size();
	std::cout << std::endl;
	printSeparated(std::cout, "\t", 12, true, "L_A", "w = 1e-4", "w = 1e-3", "w = 1e-2", "k = 0", "k = 1");
	for (int i = 3; i < this->L - 2; i++)
	{
		this->mu = dim > 3000 ? 500 : 0.25 * dim;
		u64 E_min = alfa1->E_av_idx - this->mu / 2.;
		u64 E_max = alfa1->E_av_idx + this->mu / 2.;
		double entropy_dis1 = 0.0, entropy_dis2 = 0.0, entropy_dis3 = 0.0;
		for (long k = E_min; k < E_max; k++)
		{
			auto state = alfa1->get_state_in_full_Hilbert(arma::cx_vec(alfa1->get_eigenState(k), arma::vec(dim, arma::fill::zeros)));
			entropy_dis1 += entropy::vonNeumann(state, i, alfa1->L);
			state = alfa2->get_state_in_full_Hilbert(arma::cx_vec(alfa2->get_eigenState(k), arma::vec(dim, arma::fill::zeros)));
			entropy_dis2 += entropy::vonNeumann(state, i, alfa2->L);
			state = alfa3->get_state_in_full_Hilbert(arma::cx_vec(alfa3->get_eigenState(k), arma::vec(dim, arma::fill::zeros)));
			entropy_dis3 += entropy::vonNeumann(state, i, alfa3->L);
		}
		entropy_dis1 /= double(this->mu);
		entropy_dis2 /= double(this->mu);
		entropy_dis3 /= double(this->mu);

		this->mu = beta1->get_hilbert_size() > 3000 ? 500 : 0.25 * beta1->get_hilbert_size();
		E_min = beta1->E_av_idx - this->mu / 2.;
		E_max = beta1->E_av_idx + this->mu / 2.;
		double entropy_sym1 = 0.0;
		for (long k = E_min; k < E_max; k++)
		{
			auto state = beta1->get_state_in_full_Hilbert(beta1->get_eigenState(k));
			entropy_sym1 += entropy::vonNeumann(state, i, this->L);
		}
		entropy_sym1 /= double(this->mu);

		this->mu = beta2->get_hilbert_size() > 3000 ? 500 : 0.25 * beta2->get_hilbert_size();
		E_min = beta2->E_av_idx - this->mu / 2.;
		E_max = beta2->E_av_idx + this->mu / 2.;
		double entropy_sym2 = 0.0;
		for (long k = E_min; k < E_max; k++)
		{
			auto state = beta2->get_state_in_full_Hilbert(beta2->get_eigenState(k));
			entropy_sym2 += entropy::vonNeumann(state, i, this->L);
		}
		entropy_sym2 /= double(this->mu);

		printSeparated(file, "\t", 12, true, i, entropy_dis1, entropy_dis2, entropy_dis3, entropy_sym1, entropy_sym2);
		printSeparated(std::cout, "\t", 12, true, i, entropy_dis1, entropy_dis2, entropy_dis3, entropy_sym1, entropy_sym2);
	}

	std::cout << std::endl;
}
void isingUI::ui::benchmark()
{
	auto save_eigs = [this](auto& alfa){
		auto eigenvalues = alfa.get_eigenvalues();
		std::string info = alfa.get_info({});
		std::string dir = this->saving_dir + "EIGENVALUES" + kPSep;
		createDirs(dir);
		eigenvalues.save(arma::hdf5_name(dir + kPSep + info + ".h5", "eigenvalues"));
		dir = this->saving_dir + "EIGENVECTORS" + kPSep;
		createDirs(dir);
		auto V = alfa.get_eigenvectors();
		V.save(arma::hdf5_name(dir + kPSep + info + ".h5", "eigenvectors"));
	};
	if (!this->ch)
	{
		const int Lmin = this->L, Lmax = this->L + this->Ln * this->Ls;
		int th_max = this->thread_number;
		std::ofstream file;
		std::vector th_list = {1, 2, 4, 8, 16, 24, 32, 40, 48, 64};
		std::string info = this->m ? IsingModel_sym::set_info(L, J, g, h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"L"}, ",")
								   : IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, {"L"}, ",");
		openFile(file, this->saving_dir + "benchmark" + info + ".dat", std::ios::out);
		file << "Maximum number of threads to parallelize diagonalization by OpenMP:\t " << ARMA_OPENMP_THREADS << std::endl;
		printSeparated(file, "\t", 16, true, "#cores", "chain length", "dim", "with eigenvec 'dc'", "only eigenvalues", "time in seconds");
		for (int system_size = Lmin; system_size <= Lmax; system_size += this->Ls)
		{
			for (auto &th : th_list)
			{
				if(th > this->thread_number) continue;
				omp_set_num_threads(th);
				if (this->m)
				{
					auto alfa = std::make_unique<IsingModel_sym>(system_size, this->J, this->g, this->h,
																 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
					auto start = std::chrono::system_clock::now();
					alfa->diagonalization(true, "dc");
					save_eigs(*alfa);
					double tim1 = tim_s(start);
					start = std::chrono::system_clock::now();
					alfa->diagonalization(false);
					double tim2 = tim_s(start);
					printSeparated(file, "\t", 16, true, th, system_size, alfa->get_hilbert_size(), tim1, tim2);
					printSeparated(std::cout, "\t", 16, true, th, system_size, alfa->get_hilbert_size(), tim1, tim2);
				}
				else
				{
					auto alfa = std::make_unique<IsingModel_disorder>(system_size, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
					auto start = std::chrono::system_clock::now();
					alfa->diagonalization(true, "dc");
					save_eigs(*alfa);
					double tim1 = tim_s(start);
					start = std::chrono::system_clock::now();
					alfa->diagonalization(false);
					double tim2 = tim_s(start);
					printSeparated(file, "\t", 16, true, th, system_size, alfa->get_hilbert_size(), tim1, tim2);
					printSeparated(std::cout, "\t", 16, true, th, system_size, alfa->get_hilbert_size(), tim1, tim2);
				}
			}
			file << std::endl;
		}
		file.close();
	}
	else
	{
		if (this->m)
		{
			auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
			auto start = std::chrono::system_clock::now();
			alfa->diagonalization(true, "dc");
			save_eigs(*alfa);
			double tim1 = tim_s(start);
			start = std::chrono::system_clock::now();
			alfa->diagonalization(false);
			double tim2 = tim_s(start);
			printSeparated(std::cout, "\t", 16, true, this->thread_number, this->L, alfa->get_hilbert_size(), tim1, tim2);
			if (this->thread_number == 32)
				std::cout << std::endl;
		}
		else
		{
			auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
			auto start = std::chrono::system_clock::now();
			alfa->diagonalization(true, "dc");
			save_eigs(*alfa);
			double tim1 = tim_s(start);
			start = std::chrono::system_clock::now();
			alfa->diagonalization(false);
			double tim2 = tim_s(start);
			printSeparated(std::cout, "\t", 16, true, this->thread_number, this->L, alfa->get_hilbert_size(), tim1, tim2);
			if (this->thread_number == 32)
				std::cout << std::endl;
		}
	}
}

//--------------------------------------------------------------------- SPECTRAL PROPERTIES and TIME EVOLUTION

//<! calculate all spectral quantities: time evolution, response function,
//<! integrated spectral function and spectral form factor with folded eigenvalues
void isingUI::ui::calculate_spectrals()
{
	const clk::time_point start = std::chrono::system_clock::now();
	// {"k", "x", "p"}
	std::string info = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym) 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w);

	auto [opName_tup, str] = IsingModel_disorder::opName(this->op, this->site);
	std::string opName = opName_tup;
	std::string timeDir = this->saving_dir + "timeEvolution" + kPSep + str + "=" + std::to_string(this->site) + kPSep;
	std::string specDir = this->saving_dir + "ResponseFunction" + kPSep + str + "=" + std::to_string(this->site) + kPSep;
	std::string intDir = this->saving_dir + "IntegratedResponseFunction" + kPSep + str + "=" + std::to_string(this->site) + kPSep;
	createDirs(timeDir, specDir, intDir);

	const double chi = 0.341345;
	size_t N = ULLPOW(this->L);
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		N = alfa->get_hilbert_size();
	}
	const double wH = sqrt(this->L) / (chi * N) * sqrt(this->J * this->J + this->h * this->h + this->g * this->g
												 + ( this->m? 0.0 : (this->w * this->w + this->g0 * this->g0 + this->J0 * this->J0) / 3. ));
	double tH = 1. / wH;
	int num_of_points = 1000;
	int time_end = (int)std::ceil(std::log10(3 * tH));
	time_end = (time_end / std::log10(tH) < 2.0) ? time_end + 1 : time_end;
	auto times = arma::logspace(-2, time_end, num_of_points);
	auto omegas = arma::logspace(-time_end, 2, num_of_points);

	arma::vec opEvol(times.size(), arma::fill::zeros);
	arma::vec opIntSpec(omegas.size(), arma::fill::zeros);
	arma::vec opSpecFun(N, arma::fill::zeros);
	arma::vec omega_spec(N, arma::fill::zeros);
	double Z = 0.0;
	double LTA = 0;
	arma::sp_cx_mat op;
	auto start_loop = std::chrono::system_clock::now();
	auto kernel = [&](
		auto& alfa, int r, 											//<! main parameters required by average_over_realisations
		const spectrals::preset_omega& set_omegas, const int M	//<! additional params for setting omega on logscale
		){
		std::string tdir_realisation = timeDir + "realisation=" + std::to_string(r) + kPSep;
		std::string intdir_realisation = intDir + "realisation=" + std::to_string(r) + kPSep;
		std::string specdir_realisation = specDir + "realisation=" + std::to_string(r) + kPSep;
		createDirs(tdir_realisation, intdir_realisation, specdir_realisation);
		
		stout << "\t\t	--> finished diagonalizing for " << info << " - in time : " << tim_s(start_loop) << "s" << std::endl;
		auto U = alfa.get_eigenvectors();
		
		stout << "\t\t	--> got eigenvectors for " << info << " - in time : " << tim_s(start_loop) << "s" << std::endl;
		arma::cx_mat mat_elem = U.t() * op * U;
		normaliseMat(mat_elem);
		stout << "\t\t	--> set matrix elements for " << info << " - in time : " << tim_s(start_loop) << "s" << std::endl;

		auto [op_tmp, LTA_tmp] = spectrals::autocorrelation_function(alfa, mat_elem, times);
		save_to_file(tdir_realisation + opName + info + ".dat", times, op_tmp, tH, LTA_tmp);
		stout << "\t\t	--> finished time evolution for " << info
			  << " realisation: " << r << " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s\n\tNEXT: spectral function" << std::endl;

		auto [omega_spec_r, specfun_r] = spectrals::spectralFunction(mat_elem, set_omegas, M);
		save_to_file(specdir_realisation + opName + info + ".dat", omega_spec_r, specfun_r, 1. / tH, LTA_tmp);
		stout << "\t\t	--> finished spectral function for " << info
			  << " realisation: " << r << " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s\n\tNEXT: integrated spectral function" << std::endl;
		
		auto res = spectrals::integratedSpectralFunction(alfa, mat_elem, omegas);
		save_to_file(intdir_realisation + opName + info + ".dat", omegas, res, 1. / tH, LTA_tmp);
		stout << "\t\t	--> finished integrated spectral function for " << info
			  << " realisation: " << r << " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s" << std::endl;

		LTA += LTA_tmp;
		opEvol += op_tmp;
		opIntSpec += res;
		opSpecFun += specfun_r;
		omega_spec = omega_spec_r;
		
		start_loop = std::chrono::system_clock::now();
	};
	
	// ----------- choose model and run kernel
	spectrals::preset_omega set_omega;
	int M;
	auto prefix_kernel = [&](auto& alfa){
		alfa.diagonalization();
		auto E = alfa.get_eigenvalues();
		set_omega = spectrals::preset_omega(E, 0.1 * alfa.L, E(alfa.E_av_idx));
		const int size = set_omega.get_size();
		M = int(size / double(num_of_points));
		int num = int(size / double(M));
		opSpecFun.resize(num + 1);
		omega_spec.resize(num + 1);

		op = alfa.chooseOperator(this->op, this->site);
		stout << "\n\t\t--> finished generating operator for " << info << " - in time : " << tim_s(start) << "s" << std::endl;;	
	};
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		prefix_kernel(*alfa);
		average_over_realisations<Ising_params::h>(*alfa, true, kernel, set_omega, M);
	} else{
		auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
		prefix_kernel(*alfa);
		average_over_realisations<Ising_params::h>(*alfa, true, kernel, set_omega, M);
	}

	if(this->jobid > 0) return;
	opEvol /= double(this->realisations);
	LTA /= double(this->realisations);
	opIntSpec /= double(this->realisations);
	opSpecFun /= double(this->realisations);
	save_to_file(timeDir + opName + info + ".dat", times, opEvol, tH, LTA);					smoothen_data(timeDir, opName + info + ".dat");
	save_to_file(intDir + opName + info + ".dat", omegas, opIntSpec, 1. / tH, LTA);			smoothen_data(intDir,  opName + info + ".dat");
	save_to_file(specDir + opName + info + ".dat", omega_spec, opSpecFun, 1. / tH, LTA);	smoothen_data(specDir, opName + info + ".dat");
		
	stout << " - - - - - - FINISHED CALCULATIONS IN : " << tim_s(start) << " seconds - - - - - - \n"
			  << std::endl; // simulation end
}

//<! calculate evolution of entaglement from initial state chosen by -op.
//<! -s sets the subsystem size, if-s=0 the L/2 is assumed 
void isingUI::ui::entropy_evolution(){
	clk::time_point start = std::chrono::system_clock::now();
	// ----------- generate kernel
	auto kernel = [this, start](auto& alfa){
		stout << "\n\t\t--> finished creating model for " << alfa.get_info() << " - in time : " << tim_s(start) << "s" << std::endl;
		const size_t N = alfa.get_hilbert_size();
		const double tH = 1. / alfa.mean_level_spacing_analytical();
		int t_max = (int)std::ceil(std::log10(tH));
		t_max = (t_max / std::log10(tH) < 1.5) ? t_max + 1 : t_max;

		const int division = this->site == 0? this->L / 2 : this->site;
		// ----------- predefinitions
		arma::vec times, entropy, entropy_stddev;
		double dt_new = 1e-2;
		std::string dir = this->saving_dir + "Entropy" + kPSep;
		createDirs(dir);
		alfa.reset_random(this->seed);

		// ----------- diagonalize
		stout << "\t\t	--> start diagonalizing for " << alfa.get_info()
				<< " - in time : " << tim_s(start) << "s" << std::endl;
		if(this->ch){
			dir += "Lanczos" + kPSep;
			createDirs(dir);
			auto H = alfa.get_hamiltonian();
			lanczos::Lanczos lancz(H, lanczosParams(this->mu, 1, true, false));
			//lancz.diagonalization();
			auto range1 = arma::regspace(0.01 * this->dt, 0.01 * this->dt, 0.3 * this->dt);
			auto range2 = arma::regspace(0.4 * this->dt, 0.10 * this->dt, 10.0 * this->dt);
			times = this->scale ? arma::logspace(-2, t_max, 500) : arma::join_cols(range1, range2, arma::regspace(11.0 * this->dt, this->dt, 2e2));
			entropy = arma::vec(times.size(), arma::fill::zeros); entropy_stddev = entropy;

			auto to_ave_time = [&, this, lancz](auto& alfa, int realis) mutable 
				{	// capture lancz by value to access it outsied the if-else scope, i.e. when the lambda is called
					// entropy and times are declared outside the if-else scope so they can be passed by reference
					arma::cx_vec init_state = this->set_init_state(N);
					stout << "\t\t	-->set initial state for " << alfa.get_info()
						<< " - in time : " << tim_s(start) << "s" << std::endl;
					arma::cx_vec state = init_state;
					if(this->scale)
						lancz.diagonalization(init_state);
					for (int i = 0; i < times.size(); i++)
					{
						auto t = times(i);
						if(this->scale)
							state = lancz.time_evolution_stationary(init_state, t);
						else
							lancz.time_evolution_non_stationary(state, t - (i == 0 ? 0.0 : times(i - 1)), this->mu);
						double ent = entropy::vonNeumann(alfa.get_state_in_full_Hilbert(state), division, this->L);
						entropy(i) += ent;
						entropy_stddev(i) += ent * ent;
					}
					stout << "realisation: " << realis << " - in time : " << tim_s(start) << "s" << std::endl;
				};
			average_over_realisations<Ising_params::h>(alfa, false, to_ave_time);
		} else{
			alfa.diagonalization();
			double omega_max = alfa.get_eigenvalues()(N - 1) - alfa.get_eigenvalues()(0);
			dt_new = 10 * this->dt / omega_max;
			times = this->scale ? arma::logspace(-2, t_max, 200) : arma::regspace(dt_new, dt_new, tH);
			entropy = arma::vec(times.size(), arma::fill::zeros); entropy_stddev = entropy;
			auto to_ave_time = [&](auto& alfa, int realis)
				{
					arma::cx_vec init_state = this->set_init_state(N);
					alfa.set_coefficients(init_state);
					for (int i = 0; i < times.size(); i++)
					{
						auto t = times(i);
						arma::cx_vec state = alfa.time_evolve_state(init_state, t);
						double ent = entropy::vonNeumann(alfa.get_state_in_full_Hilbert(state), division, this->L);
						entropy(i) += ent;
						entropy_stddev(i) += ent * ent;
					}
					stout << "realisation: " << realis << " - in time : " << tim_s(start) << "s" << std::endl;
				};
			average_over_realisations<Ising_params::h>(alfa, false, to_ave_time);
		}
		stout << "\t\t	--> finished diagonalizing for " << alfa.get_info()
					<< " - in time : " << tim_s(start) << "s" << std::endl;
		entropy /= double(this->realisations);
		entropy_stddev = arma::sqrt(entropy_stddev / double(this->realisations) - arma::square(entropy));
		std::ofstream file;
		openFile(file, dir + "TimeEvolution" + alfa.get_info({}) + ".dat");
		for (int j = 0; j < times.size(); j++){
			printSeparated(file, "\t", 16, true, times(j), entropy(j), entropy_stddev(j));
			printSeparated(std::cout, "\t", 16, true, times(j), entropy(j), entropy_stddev(j));
		}
		file.close();
	};
	// ----------- choose model and run kernel
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		kernel(*alfa);
	} else{
		auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
		kernel(*alfa);
	}
	
}

//<! loop over all parameters (L, site, g, h) for given disorder
//<! or symmetry sector and find relaxation times as I(w)=1/2 (the later from integrated time evolution)
void isingUI::ui::relaxationTimesFromFiles()
{
	const int Lmin = this->L, Lmax = this->L + this->Ln * this->Ls;
	auto gx_list = arma::linspace(this->g, this->g + this->gs * (this->gn - 1), this->gn);
	auto hx_list = arma::linspace(this->h, this->h + this->hs * (this->hn - 1), this->hn);
	auto kernel = [this](
		int Lx, double gx, double hx, 
		std::ofstream& map, int _site,
		auto... prints
		){
		if(_site < 0) _site = Lx / 2.;
		// read time-evolution data
		std::ifstream file, file2;
		std::string info = this->m? IsingModel_sym::set_info(this->L, this->J,gx, hx, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"J"}) 
						: IsingModel_disorder::set_info(this->L, this->J, this->J0, gx, this->g0, hx, this->w, {"J"});
		std::string opname, s;
		std::tie(opname, s) = IsingModel_disorder::opName(this->op, _site);
		std::string name = opname + info + ".dat";

		auto get_dir_str = [this, &s, &_site](const std::string& inner_dir)
			{ return this->saving_dir + inner_dir + kPSep + s + "=" + std::to_string(_site) + kPSep; };
		std::string filename 	= get_dir_str("IntegratedResponseFunction"						 ) + name;
		std::string dir_out 	= get_dir_str("IntegratedResponseFunction" + kPSep + "DERIVATIVE");
		std::string dir_norm 	= get_dir_str("IntegratedResponseFunction" + kPSep + "NORMALIZED");
		std::string dir_spec 	= get_dir_str("ResponseFunction"								 );
		std::string dir_out_2nd = get_dir_str("ResponseFunction" 		   + kPSep + "DERIVATIVE");
		createDirs(dir_out, dir_norm, dir_out_2nd);

		auto data = readFromFile(file, filename);
		file.close();
		if (data.empty()) return;
		double wH = data[2](0);
		double LTA = data.size() > 3? data[3](0) : data[1](0);
		printSeparated(std::cout, "\t", 16, true, wH, LTA);
		const size_t size = data[1].size();
		
		// take derivative
		arma::vec specFun = non_uniform_derivative(data[0], data[1]);
		arma::vec x = data[0];	x.shed_row(x.size() - 1);
		save_to_file(dir_out + name, x, specFun, wH, LTA);
		smoothen_data(dir_out, name, 10);

		// spectral function
		//smoothen_data(dir_spec, name, 100 * std::pow(this->L / 15.0, 3.0));
		auto data_spec = readFromFile(file, dir_spec + "smoothed" + kPSep + name);
		file.close();
		arma::vec omega_vals;
		arma::vec spectral_function;
		if(!data_spec.empty()){
			omega_vals = data_spec[0];
			spectral_function = data_spec[1];
			//u64 size_spec = data_spec[0].size();
		} else {
			omega_vals = x;
			spectral_function = specFun;
		}
			// find 1st plateau (peak) and normalize
		//	double omega_cut = this->g < 0.7? 0.5 : 1.0;
		//	u64 j = min_element(begin(omega_vals), end(omega_vals), [=](double a, double b) {
		//			return abs(a - omega_cut) < abs(b - omega_cut);
		//			}) - omega_vals.begin();
			u64 i = min_element(begin(omega_vals), end(omega_vals), [=](double a, double b) {
					return abs(a - wH) < abs(b - wH);
					}) - omega_vals.begin();	
//
		//	auto spectral_fun_cut = exctract_vector(spectral_function, i, j);
		//	u64 idx_tmp = i + spectral_fun_cut.index_min();
		//	if(spectral_function(idx_tmp - 1) > 0.5 * spectral_function(i))	
		//		idx_tmp = spectral_function.size() - 1;
		
		std::vector<std::vector<double>> w_cut_vals2 = {std::vector({1.5, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}),		// g=0.05
													    std::vector({0.2, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1}),		// g=0.1
													    std::vector({0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2}),		// g=0.15
													    std::vector({0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2}),		// g=0.2
													    std::vector({0.2, 0.2, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5}),		// g=0.25
													    std::vector({0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 0.5, 0.5}),		// g=0.3
													    std::vector({0.3, 0.3, 0.5, 0.6, 0.8, 0.8, 0.8, 0.8}),		// g=0.35
													    std::vector({0.4, 0.4, 0.5, 0.6, 0.8, 0.8, 0.8, 0.8}),		// g=0.4
													    std::vector({0.4, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0}),		// g=0.45
													    std::vector({0.4, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0}),		// g=0.5
													    std::vector({0.4, 0.4, 0.8, 1.0, 1.2, 1.2, 1.2, 1.2}),		// g=0.55
													    std::vector({0.4, 0.4, 0.8, 1.0, 1.2, 1.2, 1.2, 1.2}),		// g=0.6
													    std::vector({0.4, 0.4, 1.0, 1.2, 1.3, 1.3, 1.3, 1.3}),		// g=0.65
													    std::vector({0.6, 0.6, 1.0, 1.2, 1.3, 1.3, 1.3, 1.3}),		// g=0.7
													    std::vector({0.6, 0.6, 1.2, 1.4, 1.4, 1.4, 1.4, 1.4}),		// g=0.75
													    std::vector({0.8, 0.8, 1.2, 1.4, 1.4, 1.4, 1.4, 1.4}),		// g=0.8
													    std::vector({0.8, 0.8, 1.4, 10., 10., 10., 10., 1.1}),		// g=0.85
													    std::vector({1.0, 1.0, 1.4, 10., 10., 10., 10., 1.1}),		// g=0.9
													    std::vector({1.0, 1.0, 1.6, 10., 10., 10., 10., 1.1}),		// g=0.95
													    std::vector({1.1, 1.1, 1.6, 10., 10., 10., 10., 1.1}),		// g=1.0
													    std::vector({1.1, 1.1, 1.8, 10., 10., 10., 10., 1.1}),		// g=1.05
													    std::vector({1.2, 1.2, 1.8, 10., 10., 10., 10., 1.1}),		// g=1.1
													    std::vector({1.2, 1.2, 1.8, 10., 10., 10., 10., 1.1}),		// g=1.15
													    std::vector({1.3, 1.3, 1.8, 10., 10., 10., 10., 1.1}),		// g=1.2
													    std::vector({1.3, 1.3, 2.0, 10., 10., 10., 10., 1.1}),		// g=1.25
													    std::vector({1.4, 1.4, 2.0, 10., 10., 10., 10., 1.1}),		// g=1.3
													    std::vector({1.4, 1.4, 2.0, 10., 10., 10., 10., 1.1}),		// g=1.35
													    std::vector({1.5, 1.5, 2.0, 10., 10., 10., 10., 1.1}),		// g=1.4
													    std::vector({1.5, 1.5, 2.0, 10., 10., 10., 10., 1.1}),		// g=1.45
													    std::vector({1.5, 1.5, 2.0, 10., 10., 10., 10., 1.1})		// g=1.5
													    };
		
		std::vector<double> w_cut_vals = this->g <= 0.6? std::vector({0.2, 0.4, 0.6, 0.8, 1.1, 1.1, 1.1, 1.1})
													  : std::vector({10.0, 0.4, 1.0, 1.2, 1.2, 1.2, 1.2, 1.2});
		double w_cut = w_cut_vals2[int(100 * gx / 5) - 1][_site];
		auto get_minimum = [&](double w_cut) -> double {
			u64 idx_tmp = min_element(begin(omega_vals), end(omega_vals), [=](double a, double b) {
					return abs(a - w_cut) < abs(b - w_cut);
					}) - omega_vals.begin();
			if(spectral_function(idx_tmp - 1) > 0.25 * spectral_function(i))	
					idx_tmp = spectral_function.size() - 1;
			u64 idx = min_element(begin(data[0]), end(data[0]), [=](double a, double b) {
						return abs(a - omega_vals(idx_tmp)) < abs(b - omega_vals(idx_tmp));
						}) - data[0].begin();
			stout << data[0](idx) << std::endl << std::endl;
			arma::vec renorm_fun = (data[1] - LTA) / (data[0](idx) > 0.5 * wH? (data[1](idx) - LTA) : data[1](idx) );
			// save normalized data
			save_to_file(dir_norm + name, data[0], renorm_fun, wH, LTA);
			double relax = 0.0 / 0.0; // NaN
			for (int k = 0; k < renorm_fun.size(); k++) {	if (renorm_fun(k) 	>= 0.5 && k > 0){ relax = 1. / data[0](k);	break; }}
			return relax;
		};
		
		// find relax rate
		//double relax1 = 0.0 / 0.0; // NaN
		double relax1 = get_minimum(w_cut);
		double relax3 = get_minimum(10.0);
		double relax2 = 0.0 / 0.0; // NaN
		for (int k = 0; k < data[0].size(); k++)	{	if (data[1](k) 		>= 0.5 && k > 0){ relax2 = 1. / data[0](k);	break; }}
		printSeparated(std::cout, "\t", 12, false, opname, prints...);
		//printSeparated(std::cout, "\t", 12, true, relax2, 1. / wH, relax1, idx);
		printSeparated(map, "\t", 12, false, prints...);
		printSeparated(map, "\t", 12, true, relax2, 1. / wH, relax1, relax3);
	};
	std::string dir = this->saving_dir + "RelaxationTimes" + kPSep;
	createDirs(dir);

	for (this->L = Lmin; this->L < Lmax; this->L += this->Ls){
		for (int si = 0; si <= this->L / 2; si++){
				
			std::string op = std::get<0>(IsingModel_disorder::opName(this->op, si));
			std::string info = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"h", "g"}) 
						: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, {"h", "g"});
			std::ofstream map_g, map_h;
			openFile(map_g, dir + "_g" + op + info + ".dat", ios::out);
			for (auto &hx : hx_list)
				for (auto &gx : gx_list)
					kernel(this->L, gx, hx, map_g, si, hx, gx);
			map_g.close();

			openFile(map_h, dir + "_h" + op + info + ".dat", ios::out);
			for (auto &gx : gx_list)
				for (auto &hx : hx_list)
					kernel(this->L, gx, hx, map_h, si, hx, gx);
			map_h.close();
			
		}
	}
	for (auto &gx : gx_list){
		for (auto &hx : hx_list){
			for (int si = -1; si <= 2; si++){
				std::ofstream map_L;
				std::string op = std::get<0>(IsingModel_disorder::opName(this->op, si));
				std::string info = this->m? IsingModel_sym::set_info(this->L, this->J,gx, hx, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"L"}, ",") 
						: IsingModel_disorder::set_info(this->L, this->J, this->J0, gx, this->g0, hx, this->w, {"L"}, ",");
				openFile(map_L, dir + "_L" + op + info + ".dat", ios::out);
				for (this->L = Lmin; this->L < Lmax; this->L += this->Ls)
					kernel(this->L, gx, hx, map_L, si, this->L);
				map_L.close();	
			}
		}
	}
}
void isingUI::ui::intSpecFun_from_timeEvol()
{
	std::string s = (this->op < 3) ? "j" : "q";
	if(this->op == 6) s = "n";
	std::ofstream map_g, map_h;
	openFile(map_g, this->saving_dir + "RelaxationTimes_g" + IsingModel_disorder::set_info(L, 1.0, 0.0, g, 0.0, h, 0.01, {"h", "g"}) + ".dat", ios::out);
	openFile(map_h, this->saving_dir + "RelaxationTimes_h" + IsingModel_disorder::set_info(L, 1.0, 0.0, g, 0.0, h, 0.01, {"h", "g"}) + ".dat", ios::out);
	auto gx_list = arma::linspace(this->g, this->g + this->gs * (this->gn - 1), this->gn);
	auto hx_list = arma::linspace(this->h, this->h + this->hs * (this->hn - 1), this->hn);
	const double eta = 1e-3;
	for (auto &gx : gx_list)
	{
		for (auto &hx : hx_list)
		{

			// read time-evolution data
			std::ifstream file;
			std::string name = "SigmaZ_j=0" + IsingModel_disorder::set_info(L, 1.0, 0.0, gx, 0.0, hx, 0.01) + ".dat";
			std::string filename = this->saving_dir + "timeEvolution" + kPSep + s + "=" + std::to_string(this->site) + kPSep + name;
			auto data = readFromFile(file, filename);
			if (data.empty())
				continue;
			file.close();
			// integrate and find relax rate
			double wH = 1. / data[2](0);
			auto omegas = arma::logspace(std::floor(std::log10(wH)) - 2, std::log10(1. / data[0](0)), data[0].size());
			std::ofstream fileInt;
			std::string dir = this->saving_dir + "IntegratedResponseFunction" + kPSep + "integrated" + kPSep + s + "=" + std::to_string(this->site) + kPSep;
			createDirs(dir); 
			openFile(fileInt, dir + name, ios::out);
			int counter = 0;
			for (auto &w : omegas)
			{
				arma::vec toIntegrate = data[1] % arma::sin(w * data[0]) / (data[0] * pi) % arma::exp(-data[0] * eta);
				double intSpecFun = 2 * simpson_rule(data[0], toIntegrate); // even integral, only t \in (0,inf)
				if (w == omegas(0))
					printSeparated(fileInt, "\t", 12, true, w, intSpecFun, wH);
				else
					printSeparated(fileInt, "\t", 12, true, w, intSpecFun);
				if (counter == 0 && intSpecFun >= 0.5)
				{
					counter++;
					printSeparated(map_h, "\t", 12, true, hx, gx, 1. / w, wH);
					printSeparated(map_g, "\t", 12, true, gx, hx, 1. / w, wH);
				}
			}
			fileInt.close();
		}
	}
	map_g.close();
	map_h.close();
}
		
//<! analyze spectra with unfolding, DOS and level spacing distribution --  all to file
void isingUI::ui::analyze_spectra(){

	//-------PREAMBLE
	const Ising_params par = Ising_params::h;

	std::string info = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym) 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w);

	std::string dir_spacing 	= this->saving_dir + "LevelSpacingDistribution" + kPSep;
	std::string dir_DOS 		= this->saving_dir + "DensityOfStates" + kPSep;
	std::string dir_unfolding 	= this->saving_dir + "Unfolding" + kPSep;
	std::string dir_gap 		= this->saving_dir + "LevelSpacing" + kPSep + "distribution" + kPSep;
	createDirs(dir_DOS, dir_spacing, dir_unfolding, dir_gap);

	size_t N = ULLPOW(this->L);
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		N = alfa->get_hilbert_size();
	}
	const u64 num = 0.6 * N;
	double wH 				= 0.0;
	double wH_typ 			= 0.0;
	double wH_typ_unfolded 	= 0.0;

	arma::vec energies_all, energies_unfolded_all;
	arma::vec spacing, spacing_unfolded, spacing_log, spacing_unfolded_log;
	arma::vec gap_ratio, gap_ratio_unfolded;
	//-------SET KERNEL
	int counter_realis = 0;
	auto lambda_average = [&]
		(int realis, double x){
		arma::vec eigenvalues;
		std::string dir = this->saving_dir + "DIAGONALIZATION" + kPSep;
		std::string _suffix = "_real=" + std::to_string(realis);
		std::string name = dir + info + _suffix + ".hdf5";
		bool loaded;
		#pragma omp critical
		{
			loaded = eigenvalues.load(arma::hdf5_name(name, "eigenvalues/dataset"));
		}
		if(!loaded) return;
		try{
			if(eigenvalues.empty()) return;

			//------------------- Unfolding, cdf and fit	
			//if(realis == 0)
			//{
			//	arma::vec cdf(eigenvalues.size(), arma::fill::zeros);
    		//	std::iota(cdf.begin(), cdf.end(), 0);
    		//	auto p1 = arma::polyfit(eigenvalues, cdf, 3);			arma::vec res1 = arma::polyval(p1, eigenvalues);
			//	auto p2 = arma::polyfit(eigenvalues, cdf, 6);			arma::vec res2 = arma::polyval(p2, eigenvalues);
			//	auto p3 = arma::polyfit(eigenvalues, cdf, this->L);		arma::vec res3 = arma::polyval(p3, eigenvalues);
			//	auto p4 = arma::polyfit(eigenvalues, cdf, this->L + 5);	arma::vec res4 = arma::polyval(p4, eigenvalues);
			//	std::ofstream file;
			//	openFile(file, dir_unfolding + info + ".dat", std::ios::out);
			//	for(int k = 0; k < cdf.size(); k++)
			//		printSeparated(file, "\t", 14, true, eigenvalues(k), cdf(k), res1(k), res2(k), res3(k), res4(k));
			//	file.close();
			//}
			arma::vec energies_unfolded = statistics::unfolding(eigenvalues, this->L);
			//------------------- Get 50% spectrum
			double E_av = arma::trace(eigenvalues) / double(N);
			auto i = min_element(begin(eigenvalues), end(eigenvalues), [=](double x, double y) {
				return abs(x - E_av) < abs(y - E_av);
				});
			u64 E_av_idx = i - eigenvalues.begin();
			const long E_min = E_av_idx - num / 2.;
			const long E_max = E_av_idx + num / 2. + 1;

			double epsilon = sqrt(g*g + (h+w)*(h+w));
			arma::vec energies = this->ch? exctract_vector_between_values(eigenvalues, -0.5 * epsilon, 0.5 * epsilon) : 
											exctract_vector(eigenvalues, E_min, E_max);
			arma::vec energies_unfolded_cut = this->ch? statistics::unfolding(energies, this->L) :
														 exctract_vector(energies_unfolded, E_min, E_max);
			
			//------------------- Level Spacings
			arma::vec level_spacings(energies.size() - 1, arma::fill::zeros);
			arma::vec level_spacings_unfolded(energies.size() - 1, arma::fill::zeros);
			for(int i = 0; i < energies.size() - 1; i++){
				const double delta 			= energies(i+1) 				- energies(i);
				const double delta_unfolded = energies_unfolded_cut(i+1) 	- energies_unfolded_cut(i);

				wH 				 += delta / double(energies.size()-1);
				wH_typ  		 += log(abs(delta)) / double(energies.size()-1);
				wH_typ_unfolded  += log(abs(delta_unfolded)) / double(energies.size()-1);

				level_spacings(i) 			= delta;
				level_spacings_unfolded(i) 	= delta_unfolded;
			}
			arma::vec gap = statistics::eigenlevel_statistics_return(eigenvalues);
			arma::vec gap_unfolded = statistics::eigenlevel_statistics_return(energies_unfolded);

			//------------------- Combine realisations
		#pragma omp critical
			{
				energies_all = arma::join_cols(energies_all, energies);
				energies_unfolded_all = arma::join_cols(energies_unfolded_all, energies_unfolded_cut);
				
				spacing = arma::join_cols(spacing, level_spacings);
				spacing_log = arma::join_cols(spacing_log, arma::log10(level_spacings));
				spacing_unfolded = arma::join_cols(spacing_unfolded, level_spacings_unfolded);
				spacing_unfolded_log = arma::join_cols(spacing_unfolded_log, arma::log10(level_spacings_unfolded));
				
				gap_ratio = arma::join_cols(gap_ratio, gap);
				gap_ratio_unfolded = arma::join_cols(gap_ratio_unfolded, gap_unfolded);
				counter_realis++;
			}
		}
		catch (const std::exception& err) {
			stout << "Exception:\t" << err.what() << "\n";
			exit(1);
			//std::cout << "Caught some error, but current spectrum is:" << std::endl;
			//std::cout << eigenvalues.t() << std::endl;
		}
		catch (...) {
			stout << "Unknown error...!" << "\n";
			assert(false);
		}
	};

	//------CALCULATE FOR MODEL
	double norm;
	if(this->m){
		int counter = 0;
		for(int k = 1; k < this->L; k++)
		{
			if(k == this->L / 2) continue;
			this->symmetries.k_sym = k;
			average_over_realisations<par>(false, lambda_average);
			counter++;
		}
		norm = counter_realis * counter;
	} else{
		average_over_realisations<par>(false, lambda_average);
		norm = counter_realis;
	}
	if(spacing.is_empty() || spacing_log.is_empty() || spacing_unfolded.is_empty() || spacing_unfolded_log.is_empty()
							 || energies_all.is_empty() || energies_unfolded_all.is_empty()) return;
	if(spacing.is_zero() || spacing_log.is_zero() || spacing_unfolded.is_zero() || spacing_unfolded_log.is_zero()
							 || energies_all.is_zero() || energies_unfolded_all.is_zero()) return;

	wH /= norm;	wH_typ /= norm;	wH_typ_unfolded /= norm;
	std:string prefix = this->ch ? "oneband" : "";
	statistics::probability_distribution(dir_spacing, prefix + info, spacing, -1, exp(wH_typ_unfolded), wH, exp(wH_typ));
	statistics::probability_distribution(dir_spacing, prefix + "_log" + info, spacing_log, -1, wH_typ_unfolded, wH, wH_typ);
	statistics::probability_distribution(dir_spacing, prefix + "unfolded" + info, spacing_unfolded, -1, exp(wH_typ_unfolded), wH, exp(wH_typ));
	statistics::probability_distribution(dir_spacing, prefix + "unfolded_log" + info, spacing_unfolded_log, -1, wH_typ_unfolded, wH, wH_typ);
	
	statistics::probability_distribution(dir_DOS, prefix + info, energies_all, -1);
	statistics::probability_distribution(dir_DOS, prefix + "unfolded" + info, energies_unfolded_all, -1);

	statistics::probability_distribution(dir_gap, info, gap_ratio, -1);
	statistics::probability_distribution(dir_gap, info, gap_ratio_unfolded, -1);
}

//--------------------------------------------------------------------- ADIABATIC GAUGE POTENTIAL
void isingUI::ui::adiabatic_gauge_potential(){

	//----- PREAMBLE
	int counter = 0;
	double 	 
			 AGP = 0.0,	//<! adiabatic gauge potential
		typ_susc = 0.0,	//<! typical fidelity susceptibility
			susc = 0.0,	//<! fidelity susceptibility
	   diag_norm = 0.0,	//<! diagonal normalisation factor
	   	 entropy = 0.0,	//<! half-chain entropy
	  		 ipr = 0.0,	//<! inverse participation ratio
	info_entropy = 0.0,	//<! information entropy in eigenstates
	info_ent_rnd = 0.0,	//<! information entropy of random state in eigenbasis
	   gap_ratio = 0.0,	//<! gap ratio
	   		  wH = 0.0,	//<! mean level spacing
	   	  wH_typ = 0.0;	//<! typical level spacing
	std::string info;
	const long num_ent = L >= 10? 100 : 20;
	

	//---- KERNEL LAMBDA
	auto kernel = [&](auto& alfa, int realis, const arma::sp_cx_mat& opMat)
	{
		const u64 N = alfa.get_hilbert_size();
		const double omegaH = alfa.mean_level_spacing_analytical();
		const double rescale = (double)N * omegaH * omegaH / (double)L;
		this->mu = long(0.5 * N);
		info = alfa.get_info();

		long int E_min = alfa.E_av_idx - long(mu / 2);
		long int E_max = alfa.E_av_idx + long(mu / 2);
		const arma::Mat<decltype(alfa.type_var)> U = alfa.get_eigenvectors();
		const arma::vec E = alfa.get_eigenvalues();
		arma::cx_mat mat_elem = U.t() * opMat * U;

		auto [AGP_local, typ_susc_local, susc_local] = adiabatics::gauge_potential(mat_elem, E, this->L);
		typ_susc += typ_susc_local;
		AGP += AGP_local;
		susc += susc_local;

		double wH_typ_local = 0.0;
		for(int i = 0; i < N; i++){
			diag_norm += abs(mat_elem(i,i) * conj(mat_elem(i,i)));
			if(i >= E_min && i < E_max){
				const double gap1 = E(i) - E(i - 1);
				const double gap2 = E(i + 1) - E(i);
				const double min = std::min(gap1, gap2);
				const double max = std::max(gap1, gap2);
				wH += gap2;
				wH_typ_local += std::log(gap2);
        		if (abs(gap1) <= 1e-15 || abs(gap2) <= 1e-15){ 
        		    std::cout << "Index: " << i << std::endl;
        		    assert(false && "Degeneracy!!!\n");
        		}
				gap_ratio += min / max;
	
				const arma::Col<decltype(alfa.type_var)> state = U.col(i);
				ipr += statistics::inverse_participation_ratio(state);
				info_entropy += statistics::information_entropy(state);
				if(i >= alfa.E_av_idx - num_ent / 2. && i <= alfa.E_av_idx + num_ent / 2.)
					entropy += entropy::vonNeumann(cast_cx_vec(state), this->L / 2, this->L);
			}
		}
		wH_typ += std::exp(wH_typ_local / double(this->mu));

		auto state = this->set_init_state(N, 0);
		info_ent_rnd += statistics::information_entropy(state);
		counter++;
	};

	//---- START COMPUTATION
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		arma::sp_cx_mat opMat = alfa->chooseOperator(this->op, this->site);
		average_over_realisations<Ising_params::h>(*alfa, true, kernel, opMat);
	} else{
		auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
		arma::sp_cx_mat opMat = alfa->chooseOperator(this->op, this->site);
		average_over_realisations<Ising_params::h>(*alfa, true, kernel, opMat);
	}

	auto [opName, str] = IsingModel_disorder::opName(this->ch, this->site);
	std::string dir = this->saving_dir + "STATISTICS" + kPSep;
	std::string dir_agp = this->saving_dir + "AGP" + kPSep + opName + kPSep;
	createDirs(dir, dir_agp);
	std::ofstream file;

	openFile(file, dir + info + "_jobid=" + std::to_string(jobid) + ".dat");
	printSeparated(file, "\t", 25, true, "gap ratio", 			   				gap_ratio / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "ipr", 						 			  ipr / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "information entropy", 			 info_entropy / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "information entropy random state", info_ent_rnd / double(counter));
	printSeparated(file, "\t", 25, true, "entropy in ~100 states at E=0", 	  	  entropy / double(num_ent * counter));
	printSeparated(file, "\t", 25, true, "mean level spacing", 			  			   wH / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "typical level spacing", 	  			   wH_typ / double(counter));
	file.close();

	openFile(file, dir_agp + info + "_jobid=" + std::to_string(jobid) + ".dat");
	printSeparated(file, "\t", 25, true, "adiabatic gauge potential", 	 			  AGP / double(counter));
	printSeparated(file, "\t", 25, true, "typical fidelity susceptibility", 	 typ_susc / double(counter));
	printSeparated(file, "\t", 25, true, "fidelity susceptibility", 			 	 susc / double(counter));
	printSeparated(file, "\t", 25, true, "diagonal normalisation factor", 		diag_norm / double(counter));

	file.close();
}


//-------------------------------------------------------------------------- STATISTICS
//<! spectral form factor calculated from eigenvalues in file or diagonalize matrix
void isingUI::ui::spectral_form_factor(){
	// always unfolding!
	this->ch = 1;
	clk::time_point start = std::chrono::system_clock::now();

	std::string dir = this->saving_dir + "SpectralFormFactor" + kPSep;
	std::string dir2 = this->saving_dir + "LevelSpacing" + kPSep + "raw_data" + kPSep;
	createDirs(dir, dir2);
	//------- PREAMBLE
	std::string info = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"k", "x", "p"}) 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w);

	const double chi = 0.341345;
	#ifdef HEISENBERG
		size_t dim = binomial(this->L, this->L / 2.);
	#else
		size_t dim = ULLPOW(this->L);
	#endif
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		dim = alfa->get_hilbert_size();
	}
	const double wH = sqrt(this->L) / (chi * dim) * sqrt(this->J * this->J + this->h * this->h + this->g * this->g
												 + ( this->m? 0.0 : (this->w * this->w + this->g0 * this->g0 + this->J0 * this->J0) / 3. ));
	double tH = 1. / wH;
	double r1 = 0.0, r2 = 0.0;
	int num_times = 5000;
	int time_end = (int)std::ceil(std::log10(5 * tH));
	time_end = (time_end / std::log10(tH) < 1.5) ? time_end + 1 : time_end;

	arma::vec times = arma::logspace(log10(1.0 / (two_pi * dim)), 2, num_times);
	arma::vec times_fold = arma::logspace(-2, time_end, num_times);
	arma::vec sff(num_times, arma::fill::zeros);
	double Z = 0.0;
	double wH_mean = 0.0;
	double wH_typ  = 0.0;
	const Ising_params par = Ising_params::h;
	// ------ SET LAMBDA
	auto lambda_average = [&](
		int realis, double x
		)
	{
		double Jx, gx, hx;
		switch (par)
		{
			case Ising_params::J: Jx = x; 		gx = this->g; 	hx = this->h;	break;
			case Ising_params::g: Jx = this->J; gx = x;			hx = this->h;	break;
			case Ising_params::h: Jx = this->J; gx = this->g;	hx = x;			break;
		default: 				  Jx = 0.0; 	gx = 0.0;		hx = 0.0;		break;
		}
		arma::vec eigenvalues;
		std::string suffix = "_real=" + std::to_string(realis + this->jobid);
		if(this->m){
			auto alfa = std::make_unique<IsingModel_sym>(this->L, Jx, gx, hx,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
			eigenvalues = this->get_eigenvalues(*alfa, suffix);
			info = IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym);
		} else{
			auto alfa = std::make_unique<IsingModel_disorder>(this->L, Jx, this->J0, gx, this->g0, hx, this->w, this->boundary_conditions);
			eigenvalues = this->get_eigenvalues(*alfa, suffix);
		}
		if(this->fun == 3) stout << "\t\t	--> finished loading eigenvalues for " << info + suffix << " - in time : " << tim_s(start) << "s" << std::endl;
		if(eigenvalues.empty()) return;
		const u64 N = eigenvalues.size();
		// ------------------------------------- calculate gap ratio
			double E_av = arma::trace(eigenvalues) / double(N);
			const u64 num = this->L <= 9? 0.25 * N : 0.5 * N;
			const u64 num2 = this->L <= 10? 50 : 500;

			auto i = min_element(begin(eigenvalues), end(eigenvalues), [=](double x, double y) {
				return abs(x - E_av) < abs(y - E_av);
				});
			u64 E_av_idx = i - eigenvalues.begin();

		// ------------------------------------- calculate level statistics
			double r1_tmp = statistics::eigenlevel_statistics((E_av_idx - num / 2) + eigenvalues.begin(), (E_av_idx + num / 2) + eigenvalues.begin());
			double r2_tmp = statistics::eigenlevel_statistics((E_av_idx - num2 / 2) + eigenvalues.begin(), (E_av_idx + num2 / 2) + eigenvalues.begin());
			double wH_mean_r = statistics::mean_level_spacing((E_av_idx - num / 2) + eigenvalues.begin(), (E_av_idx + num / 2) + eigenvalues.begin());
			double wH_typ_r = statistics::typical_level_spacing((E_av_idx - num / 2) + eigenvalues.begin(), (E_av_idx + num / 2) + eigenvalues.begin());
			if(this->fun == 3) stout << "\t\t	--> finished unfolding for " << info + suffix << " - in time : " << tim_s(start) << "s" << std::endl;
			
			auto [sff_r_folded, Z_r_folded] = statistics::spectral_form_factor(eigenvalues, times_fold, 0.5);
			eigenvalues = statistics::unfolding(eigenvalues);

			auto [sff_r, Z_r] = statistics::spectral_form_factor(eigenvalues, times, 0.5);
			#pragma omp critical
			{
				r1 += r1_tmp;
				r2 += r2_tmp;
				sff += sff_r;
				Z += Z_r;
				wH_mean += wH_mean_r;
				wH_typ  += wH_typ_r;
			}
		if(this->fun == 3) stout << "\t\t	--> finished realisation for " << info + suffix << " - in time : " << tim_s(start) << "s" << std::endl;
		
		//--------- SAVE REALISATION TO FILE
		#if !defined(MY_MAC)
			std::string dir_re  = this->saving_dir + "SpectralFormFactor" + kPSep + "realisation=" + std::to_string(this->jobid + realis) + kPSep;
			createDirs(dir_re);
			save_to_file(dir_re + info + ".dat", 			times, sff_r, 		 Z_r, 		 r1_tmp, r2_tmp, wH_mean_r, wH_typ_r);
			save_to_file(dir_re + "folded" + info + ".dat", times, sff_r_folded, Z_r_folded, r1_tmp, r2_tmp, wH_mean_r, wH_typ_r);
		#else
			stout << this->jobid + realis << "\t";
			stout << std::endl;
		#endif
	};
	
	// ----------- choose model and run kernel
	double norm = 0.0;
	if(this->m){
		int counter = 0;
		for(int k = 1; k < this->L; k++)
		{
			if(k == this->L / 2) continue;
			this->symmetries.k_sym = k;
			average_over_realisations<par>(false, lambda_average);
			counter++;
		}
		norm = this->realisations * counter;
	} else{
		average_over_realisations<par>(false, lambda_average);
		norm = this->realisations;
	}
	if(sff.is_empty()) return;
	if(sff.is_zero()) return;
	r1 /= norm;
	r2 /= norm;
	sff = sff / Z;
	wH_mean /= norm;
	wH_typ /= norm;
	
	if(this->jobid > 0) return;
	std::ofstream lvl;
	openFile(lvl, dir2 + info + ".dat", std::ios::out);
	printSeparated(lvl, "\t", 16, true, r1, r2);
	lvl.close();
	// ---------- find Thouless time
	double eps = 8e-2;
	auto K_GOE = [](double t){
		return t < 1? 2 * t - t * log(1+2*t) : 2 - t * log( (2*t+1) / (2*t-1) );
	};
	double thouless_time = 0;
	const double t_max = this->ch? 2.5 : 2.5 * tH;
	double delta_min = 1e6;
	for(int i = 0; i < sff.size(); i++){
		double t = this->ch? times(i) : times(i) * tH;
		double delta = abs(log10( sff(i) / K_GOE(t) )) - eps;
		delta *= delta;
		if(delta < delta_min){
			delta_min = delta;
			thouless_time = times(i); 
		}
		if(times(i) >= t_max) break;
	}
	save_to_file(dir + info + ".dat", times, sff, 1.0 / wH_mean, thouless_time, r1, r2, dim, 1.0 / wH_typ);
	smoothen_data(dir, info + ".dat");
}

void isingUI::ui::average_SFF(){
	std::string dir = this->saving_dir + "SpectralFormFactor" + kPSep;
	std::string dir2 = this->saving_dir + "LevelSpacing" + kPSep + "raw_data" + kPSep;
	
	const Ising_params par = Ising_params::h;
	//------- PREAMBLE
	std::string info = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym) 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w);

	const double chi = 0.341345;
	#ifdef HEISENBERG
		size_t dim = binomial(this->L, this->L / 2.);
	#else
		size_t dim = ULLPOW(this->L);
	#endif
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		dim = alfa->get_hilbert_size();
	}
	const double wH = sqrt(this->L) / (chi * dim) * sqrt(this->J * this->J + this->h * this->h + this->g * this->g
												 + ( this->m? 0.0 : (this->w * this->w + this->g0 * this->g0 + this->J0 * this->J0) / 3. ));
	double tH = 1. / wH;
	arma::vec times, times_fold; // are always the same
	arma::vec sff(5000, arma::fill::zeros);
	arma::vec sff_fold(5000, arma::fill::zeros);
	double Z = 0.0;
	double Z_folded = 0.0;
	double r1 = 0.0;
	double r2 = 0.0;
	// ------ SET LAMBDA
	int counter_realis = 0;
	auto lambda_average = [&](
		int realis, double x
		)
	{
		std::string dir_re  = this->saving_dir + "SpectralFormFactor" + kPSep + "realisation=" + std::to_string(this->jobid + realis) + kPSep;
		std::ifstream file;

		auto data = readFromFile(file, dir_re + info + ".dat");
		if(data.empty()) return;
		if(data[0].size() != sff.size()) {
			std:cout << "Incompatible data dimensions" << std::endl;
			return;
		}
		file.close();
		#pragma omp critical
		{
			times = data[0];
			sff += data[1];
			Z += data[2](0);
			r1 += data[3](0);
			r2 += data[4](0);
			counter_realis++;
		}

		data = readFromFile(file, dir_re + "folded" + info + ".dat");
		if(data.empty()) return;
		if(data[0].size() != sff_fold.size()) {
			std::cout << "Incompatible data dimensions" << std::endl;
			return;
		}
		file.close();
		#pragma omp critical
		{
			times_fold = data[0];
			sff_fold += data[1];
			Z_folded += data[2](0);
		}
	};
	double norm = 0.0;
	if(this->m){
		int counter = 0;
		for(int k = 1; k < this->L; k++)
		{
			if(k == this->L / 2) continue;
			this->symmetries.k_sym = k;
			average_over_realisations<par>(false, lambda_average);
			counter++;
		}
		norm = counter_realis * counter;
	} else{
		average_over_realisations<par>(false, lambda_average);
		norm = counter_realis;
	}

	if(sff.is_empty()) return;
	if(sff.is_zero()) return;
	r1 /= norm;
	r2 /= norm;
	sff = sff / Z;

	std::ofstream lvl;
	openFile(lvl, dir2 + info + ".dat", std::ios::out);
	printSeparated(lvl, "\t", 16, true, r1, r2);
	lvl.close();
	// ---------- find Thouless time
	double eps = 5e-2;
	auto K_GOE = [](double t){
		return t < 1? 2 * t - t * log(1+2*t) : 2 - t * log( (2*t+1) / (2*t-1) );
	};
	double thouless_time = 0;
	double delta_min = 1e6;
	for(int i = 0; i < sff.size(); i++){
		double delta = abs(log10( sff(i) / K_GOE(times(i)) )) - eps;
		delta *= delta;
		if(delta < delta_min){
			delta_min = delta;
			thouless_time = times(i); 
		}
		if(times(i) >= 2.5) break;
	}
	save_to_file(dir + info + ".dat", times, sff, tH, thouless_time, r1, r2, dim);	
	smoothen_data(dir, info + ".dat");

	//---------- FOLDED SFF
	if(sff_fold.is_empty()) return;
	if(sff_fold.is_zero()) return;
	sff_fold = sff_fold / Z_folded;
	double thouless_time_fold = 0;
	delta_min = 1e6;
	for(int i = 0; i < sff_fold.size(); i++){
		double delta = abs(log10( sff_fold(i) / K_GOE(times_fold(i)) )) - eps;
		delta *= delta;
		if(delta < delta_min){
			delta_min = delta;
			thouless_time_fold = times_fold(i); 
		}
		if(times_fold(i) >= 2.5 * tH) break;
	}
	save_to_file(dir + "folded" + info + ".dat", times_fold, sff_fold, tH, thouless_time_fold, r1, r2, dim);	
	smoothen_data(dir, "folded" + info + ".dat");
}

//<! find thouless time with various method as function of h,g,J
void isingUI::ui::thouless_times()
{
	const int Lmin = this->L, Lmax = this->L + this->Ln * this->Ls;
	auto gx_list = arma::linspace(this->g, this->g + this->gs * (this->gn - 1), this->gn);
	auto hx_list = arma::linspace(this->h, this->h + this->hs * (this->hn - 1), this->hn);
	auto wx_list = arma::linspace(this->w, this->w + this->ws * (this->wn - 1), this->wn);
	auto k_list = arma::linspace(0, this->L, this->L);
	auto _list = this->m? k_list : wx_list;
	std::cout << _list << std::endl;
	auto kernel = [this](
		int Lx, double Jx, double gx, double hx, double x,
		std::ofstream& map, auto... prints
		){
		std::ifstream file;
		std::string info = this->m? IsingModel_sym::set_info(Lx, Jx ,gx, hx, x, this->symmetries.p_sym, this->symmetries.x_sym) 
						: IsingModel_disorder::set_info(Lx, Jx, this->J0, gx, this->g0, hx, x);
		std::string filename = this->saving_dir + "SpectralFormFactor/smoothed" + kPSep + info + ".dat";
		auto data = readFromFile(file, filename);
		file.close();
		if (data.empty()) return;
		arma::vec times = data[0];
		arma::vec sff = data[1];
		double tH = data[2](0);
		double r1 = 0.0, r2 = 0.0;
		size_t dim = 0;
		if(data.size() > 4){
			r1 = data[4](0);
			r2 = data[5](0);
		}
		if(data.size() > 6)
			dim = data[6](0);
		// find thouless time
		double eps = 8e-2;
		auto K_GOE = [](double t){
			return t < 1? 2 * t - t * log(1+2*t) : 2 - t * log( (2*t+1) / (2*t-1) );
		};
		double thouless_time = 0;
		const double t_max = this->ch? 2.5 : 2.5 * tH;
		double delta_min = 1e6;
		for(int i = 0; i < sff.size(); i++){
			double t = this->ch? times(i) : times(i) / tH;
			double delta = abs( log10( sff(i) / K_GOE(t) ));
			//if(delta < eps){
			//	thouless_time = t;
			//	break;
			//}
			delta = delta - eps;
			delta *= delta;
			if(delta < delta_min){
				delta_min = delta;
				thouless_time = times(i); 
			}

			if(times(i) >= t_max) break;
		}
		printSeparated(std::cout, "\t", 12, false, prints...);
		printSeparated(std::cout, "\t", 12, true, thouless_time, tH);
		printSeparated(map, "\t", 12, false, prints...);
		printSeparated(map, "\t", 12, true, thouless_time, tH, r1, r2, dim);
	};
	std::string dir = this->saving_dir + "ThoulessTime" + kPSep;
	createDirs(dir);
	std::string info = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"h", "g"}) 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, {"h", "g", "L", "J", "w"}, ",");
	std::ofstream map;
	openFile(map, dir + "_all" + info + ".dat", ios::out);
	for (int size = Lmin; size < Lmax; size += this->Ls){
		for(double Jx = this->J; Jx <= 1.55; Jx += 0.05){
			for (auto &gx : gx_list){
				for (auto &hx : hx_list){
					for(auto& x : _list){	// either disorder w (m=0) or symmetry sector k (m=1)
						kernel(size, Jx, gx, hx, x, map, size, Jx, gx, hx, x);
		}}}}}
	map.close();

	return;
	for (int size = Lmin; size < Lmax; size += this->Ls){
		std::string info = this->m? IsingModel_sym::set_info(size, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"h", "g"}) 
					: IsingModel_disorder::set_info(size, this->J, this->J0, this->g, this->g0, this->h, this->w, {"h", "g"});
		std::ofstream map_g, map_h;
		openFile(map_g, dir + "_g" + info + ".dat", ios::out);
		for (auto &hx : hx_list)
			for (auto &gx : gx_list)
				kernel(size, this->J, gx, hx, this->w, map_g, hx, gx);
		map_g.close();
		openFile(map_h, dir + "_h" + info + ".dat", ios::out);
		for (auto &gx : gx_list)
			for (auto &hx : hx_list)
				kernel(size, this->J, gx, hx, this->w, map_h, hx, gx);
		map_h.close();
	}
	for (auto &gx : gx_list){
		for (auto &hx : hx_list){
			std::ofstream map_L;
			std::string info = this->m? IsingModel_sym::set_info(this->L, this->J,gx, hx, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"L"}, ",") 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, gx, this->g0, hx, this->w, {"L"}, ",");
			openFile(map_L, dir + "_L" + info + ".dat", ios::out);
			for (int size = Lmin; size < Lmax; size += this->Ls)
				kernel(size, this->J, gx, hx, this->w, map_L, size);
			map_L.close();	
		}
	}
}


void isingUI::ui::level_spacing(){
	clk::time_point start = std::chrono::system_clock::now();
	
	const int Lmin = this->L, Lmax = this->L + this->Ln * this->Ls;
	const double gmin = this->g, gmax = this->g + this->gn * this->gs;
	const double hmin = this->h, hmax = this->h + this->hn * this->hs;

	//std::string dir_raw = this->saving_dir + "LevelSpacing" + kPSep + "raw_data" + kPSep;
	std::string dir_raw = this->saving_dir + "SpectralFormFactor" + kPSep;
	std::string dir_ratio = this->saving_dir + "LevelSpacing" + kPSep + "ratio" + kPSep;
	std::string dir_dist = this->saving_dir + "LevelSpacing" + kPSep + "distribution" + kPSep;
	createDirs(dir_dist, dir_ratio);
	
	//----------- SET LAMBDA
	auto kernel = [this, &start, &dir_dist, &dir_raw](
		std::ofstream& file, double par1, double par2,
		int Lx, double gx, double hx)
	{
		std::string info = this->m? IsingModel_sym::set_info(Lx, this->J, gx, hx, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym) 
				: IsingModel_disorder::set_info(Lx, this->J, this->J0, gx, this->g0, hx, this->w);

	#ifdef HEISENBERG
		size_t dim = binomial(Lx, Lx / 2.);
	#else
		size_t dim = ULLPOW(Lx);
	#endif
		if(this->m){
			auto alfa = std::make_unique<IsingModel_sym>(Lx, this->J, gx, hx,
									 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
			dim = alfa->get_hilbert_size();
		}
		double r1 = 0, r2 = 0;
		double _min = 0, _max = 0;
		std::ifstream raw_file;
		auto r_vals = readFromFile(raw_file, dir_raw + info  + ".dat");
		raw_file.close();
		if(!r_vals.empty()){
			r1 = r_vals[4](0);
			r2 = r_vals[5](0);
		} else {
			long n_bins = 1 + long(3.322 * log(dim - 2));
			arma::vec lvl_prob_dist(n_bins, arma::fill::zeros);
			auto lambda_average = [this, &lvl_prob_dist, &r1, &r2, &_min, &_max, &Lx, &gx, &hx, &n_bins](int realis, double Jx)
			{
				arma::vec eigenvalues;
				std::string suffix = (this->m)? "" : "_real=" + std::to_string(realis);
				if(this->m){
					auto alfa = std::make_unique<IsingModel_sym>(Lx, Jx, gx, hx,
											 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
					eigenvalues = this->get_eigenvalues(*alfa, suffix);
				} else{
					auto alfa = std::make_unique<IsingModel_disorder>(Lx, Jx, this->J0, gx, this->g0, hx, this->w, this->boundary_conditions);
					eigenvalues = this->get_eigenvalues(*alfa, suffix);
				};
				if(eigenvalues.empty()) return;

				const u64 N = eigenvalues.size();
				double E_av = arma::trace(eigenvalues) / double(N);
				const u64 num = this->L <= 8? 0.1 * N : 0.25 * N;
				const u64 num2 = (this->L < 10? (this->L <= 8? 10 : 50) : 200);

				auto i = min_element(begin(eigenvalues), end(eigenvalues), [=](double x, double y) {
					return abs(x - E_av) < abs(y - E_av);
					});
				u64 E_av_idx = i - eigenvalues.begin();
				auto lvl_spacing = statistics::eigenlevel_statistics_return(eigenvalues);
				#pragma omp critical
				{
					r1 += statistics::eigenlevel_statistics((E_av_idx - num) + eigenvalues.begin(), (E_av_idx + num) + eigenvalues.begin());
					r2 += statistics::eigenlevel_statistics((E_av_idx - num2) + eigenvalues.begin(), (E_av_idx + num2) + eigenvalues.begin());
					lvl_prob_dist += statistics::probability_distribution(lvl_spacing, n_bins);
					_min = arma::min(lvl_spacing);
					_max = arma::max(lvl_spacing);
				}
			};

			// ----- SET MODEL
			double norm = 0.0;
			if(this->m){
				int counter = 0;
				for(int k = 1; k < this->L; k++)
				{
					if(k == this->L / 2) continue;
					this->symmetries.k_sym = k;
					average_over_realisations<Ising_params::J>(false, lambda_average);
					counter++;
				}
				norm = this->realisations * counter;
			} else{
				average_over_realisations<Ising_params::J>(false, lambda_average);
				norm = this->realisations;
			}
			r1 /= norm;
			r2 /= norm;

			// save distribution of level spacing
			lvl_prob_dist /= norm;
			auto x_vals = arma::linspace(_min, _max, lvl_prob_dist.size());
			save_to_file(dir_dist + info + ".dat", x_vals, lvl_prob_dist, r1);
		}
		printSeparated(file, "\t", 14, true, par1, par2, r1, r2);
		
		stout << "\t\t	--> finished realisations for " << info << " - in time : " << tim_s(start) << "s" << std::endl;
	};

	std::ofstream file;
	for (int system_size = Lmin; system_size < Lmax; system_size += this->Ls)
	{
		std::string info = this->m ? IsingModel_sym::set_info(system_size, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"g", "h"})
						   : IsingModel_disorder::set_info(system_size, this->J, this->J0, this->g, this->g0, this->h, this->w, {"g", "h"});
		openFile(file, dir_ratio + "hMap" + info + ".dat", std::ios::out);
		for (double gx = gmin; gx < gmax; gx += this->gs)
			for (double hx = hmin; hx < hmax; hx += this->hs)
				kernel(file, gx, hx, system_size, gx, hx);
		file.close();	

		openFile(file, dir_ratio + "gMap" + info + ".dat", std::ios::out);
		for (double hx = hmin; hx < hmax; hx += this->hs)
			for (double gx = gmin; gx < gmax; gx += this->gs)
				kernel(file, hx, gx, system_size, gx, hx);
		file.close();	
	}
}

void isingUI::ui::smoothen_data(const std::string& dir, const std::string& name, int bucket_size){
	if(bucket_size < 0)
		bucket_size = this->mu;
	// read data
	std::ifstream input;
	auto data = readFromFile(input, dir + name);
	if(data.empty()) return;
	input.close();
	// smoothed output -- twice
	auto smoothed_data = statistics::remove_fluctuations(data[1], bucket_size);
	//smoothed_data = statistics::remove_fluctuations(smoothed_data, bucket_size);
	// save smoothed
	std::ofstream output;
	std::string new_dir = dir + "smoothed" + kPSep;
	createDirs(new_dir);
	openFile(output, new_dir + name, std::ios::out);
	for(int k = 0; k < data[0].size(); k++){
			printSeparated(output, "\t", 14, false, data[0](k), smoothed_data(k));
		for(int i = 2; i < data.size(); i++){
			printSeparated(output, "\t", 14, false, data[i](k));
		}
		output << std::endl;
	}
	output.close();
}


void isingUI::ui::calculate_statistics(){
	//----- PREAMBLE
	int counter = 0;
	double 	 
	   	 entropy = 0.0,	//<! half-chain entropy
	  		 ipr = 0.0,	//<! inverse participation ratio
	info_entropy = 0.0,	//<! information entropy in eigenstates
	info_ent_rnd = 0.0,	//<! information entropy of random state in eigenbasis
	   gap_ratio = 0.0,	//<! gap ratio
	   		  wH = 0.0,	//<! mean level spacing
	   	  wH_typ = 0.0;	//<! typical level spacing
	std::string info;
	
	const long num_ent = L >= 10? 100 : 20;

	//---- KERNEL LAMBDA
	auto kernel = [&](auto& alfa, int realis)
	{
		const u64 N = alfa.get_hilbert_size();
		const double omegaH = alfa.mean_level_spacing_analytical();
		const double rescale = (double)N * omegaH * omegaH / (double)L;
		this->mu = long(0.5 * N);
	
		info = alfa.get_info();

		long int E_min = alfa.E_av_idx - long(mu / 2);
		long int E_max = alfa.E_av_idx + long(mu / 2);
		const arma::Mat<decltype(alfa.type_var)> U = alfa.get_eigenvectors();
		const arma::vec E = alfa.get_eigenvalues();

		double wH_typ_local = 0.0;
		for(int i = 0; i < N; i++){
			if(i >= E_min && i < E_max){
				const double gap1 = E(i) - E(i - 1);
				const double gap2 = E(i + 1) - E(i);
				const double min = std::min(gap1, gap2);
				const double max = std::max(gap1, gap2);
				wH += gap2;
				wH_typ_local += std::log(gap2);
        		if (abs(gap1) <= 1e-15 || abs(gap2) <= 1e-15){ 
        		    std::cout << "Index: " << i << std::endl;
        		    assert(false && "Degeneracy!!!\n");
        		}
				gap_ratio += min / max;
	
				const arma::Col<decltype(alfa.type_var)> state = U.col(i);
				ipr += statistics::inverse_participation_ratio(state);
				info_entropy += statistics::information_entropy(state);
				if(i >= alfa.E_av_idx - num_ent / 2. && i <= alfa.E_av_idx + num_ent / 2.)
					entropy += entropy::vonNeumann(cast_cx_vec(state), this->L / 2, this->L);
			}
		}
		wH_typ += std::exp(wH_typ_local / double(this->mu));

		auto state = this->set_init_state(N);
		info_ent_rnd += statistics::information_entropy(state);
		counter++;
	};

	//---- START COMPUTATION
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		average_over_realisations<Ising_params::h>(*alfa, true, kernel);
	} else{
		auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
		average_over_realisations<Ising_params::h>(*alfa, true, kernel);
	}

	std::string dir = this->saving_dir + "STATISTICS" + kPSep;
	createDirs(dir);
	std::ofstream file;
	openFile(file, dir + info + "_jobid=" + std::to_string(jobid) + ".dat");

	printSeparated(file, "\t", 25, true, "gap ratio", 			   				gap_ratio / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "ipr", 						 			  ipr / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "information entropy", 			 info_entropy / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "information entropy random state", info_ent_rnd / double(counter));
	printSeparated(file, "\t", 25, true, "entropy in ~100 states at E=0", 	  	  entropy / double(num_ent * counter));
	printSeparated(file, "\t", 25, true, "mean level spacing", 			  			   wH / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "typical level spacing", 	  			   wH_typ / double(counter));

	file.close();
}



//---------------------------------------------------------------------------------------------------------------- IMPLEMENTATION OF UI
//---------------------------------------------------------------------------------------------------------------- FUNCTIONS AND MORE
//---------------------------------------------------------------------------------------------------------------- 

void print_help(){
	printf(
		" Usage: name of executable [options] outputDirectory \n"
		" The input can be both introduced with [options] described below or with giving the input directory(which also is the flag in the options)\n"
		" options:\n"
		"-f input file for all of the options : (default none)\n"
		"-mu bucket size for ergodic coefficients (default 5)\n"
		"-J spin exchange coefficient : (default 1)\n"
		"-J0 random spin exchange set in uniform distribution [-J0,J0]\n"
		"-g transverse magnetic field (x-) constant: (default 1)\n"
		"-gs transverse magnetic field (x-) constant step: (default 0.0)\n"
		"-gn transverse magnetic field (x-) constant number: (default 1)\n"
		"-g0 random transverse field set in uniform distribution [-g0,g0]\n"
		"-g0s transverse disorder strength step: (default 0.0)\n"
		"-g0n transverse disorder strength number: (default 1)\n"
		"-h perpendicular (z-) magnetic field constant: (default 0)\n"
		"-hs perpendicular (z-) magnetic field constant step: (default 0.0)\n"
		"-hn perpendicular (z-) magnetic field constant number: (default 1)\n"
		"-w disorder strength : (default 0 - no disorder introduced)\n"
		"-ws disorder strength step: (default 0.0)\n"
		"-wn disorder strength number: (default 1)\n"
		"-L chain length minimum: bigger than 0 (default 8)\n"
		"-Ls chain length step: bigger equal than 0 (default 0)\n"
		"-Ln chain length number: bigger than 0 (default 1)\n"
		"-b boundary conditions : bigger than 0 (default 0 - PBC)\n"
		"	0 -- PBC\n"
		"	1 -- OBC\n"
		"	2 -- ABC -- none so far implemented\n"
		"-s site to act with local operators (default 0)\n"
		"-op flag to choose operator: \n"
		"	0 -- Sz_i-local\n"
		"	1 -- Sx_i-local\n"
		"	2 -- Hi\n"
		"	3 -- Sz_q\n"
		"	4 -- Sx_q\n"
		"	5 -- Hq\n"
		"	6 -- I+_n - TFIM LIOMs ordered by locality\n"
		"	7 -- I-_n - TFIM LIOMs ordered by locality\n"
		"	8 -- A_i - non-interacting (J=0) LIOMs\n"
		"	  -- to get sum of local Sz or Sx take Sz_q or Sx_q with -s=0\n"
		"	  -- i or q are set to site (flag -s); (default 0)\n"
		""
		"-fun choose function to start calculations: check user_interface.cpp -> make_sim() to find functions\n"
		"\t	0 -- diagonalizing hamiltonian and writing to file eigenvalues. Set -ch=1 to include eigenvector calculation\n"
		"\t	1 -- time evolution (and spectral functions) for any model (disordered is with averaging):\n"
		"\t\t set -op for operator and -s for acting site\n"
		"\t	2 -- evolution of entropy from initial state chosen by the -op flag:\n"
		"\t		for both models (-m flag) and possible to use lanczos iterative method by setting -ch=1\n"
		"\t		use -mu to set number of lanczos steps (<10 is enough) and -ts as time step (divided by E_max - E_min): 0.1 is sufficient\n"
		"\t			* op=0 -- random initial product state averaged over -r realisations\n"
		"\t			* op=1 -- fully ferromagnetically polarised state |111111...>\n"
		"\t			* op=2 -- fully anti-ferromagnetically polarised state |111111...>\n"
		"\t	3 -- spectral form factor calculation (checks if file exists with data, if not then diagonalize and save\n"
		"\t\t is looped over h, g and L, set Ls, Gs, hs, Ln, gn, hn or use defaults and only for specific g, h, L find ssf\n"
		"\t	4 -- relaxation times from integrated spectral function for:\n"
		"\t\t operator -op flag on site -s flag\n"
		"\t\t (also derivative of integrated spectral function is calculated)\n"
		"\t\t looped over system sizes: -L, -Ls, -Ln and sites: from 0 to L/2\n"
		"\t 5 -- benchmark diagonalization routines vs CPU count:\n"
		"\t\t looped over different system sizes set by -L, -Ln, -Ls\n"
		"\t\t for number of threads: 1, 2, 4, 8, 16, 24, 32, 40, 48, 64\n"
		"\t	6 -- AGPs for small disorder (-m=0) as function of h for -ch=1 or as function of g for -ch=0 for input operator from -op flag\n"
		"\t\t SET: -L, -Ln, -Ls, -h, -hn, -hs, -op, -w(default=0.01)\n"
		"\t 7 -- calculate gap ratio <r> either from input file or diagonalize matrix otherwise.\n"
		"\t default -- in make_sim space for user to write function; designed for non-builtin behavior\n"
		""
		"-m model to be choosen : (default 0 - without symmetries)\n"
		"	0 -- nonsymmetric model - only here the disorder is working\n"
		"	1 -- include symmetries - here the parity flag is also working\n"
		"-k translation symetry sector, 0-L, (default 0)\n"
		"-p parity symmetry sector, +-1 (if applicable) (default 1)\n"
		"-x spin flip symmetry sector, +-1 (if applicable) (default 1)\n"
		"-th number of threads to be used for CPU parallelization : depends on the machine specifics, default(1)\n"
		"-ch general boolean flag used in different context (default: 0)\n"
		"-ts time step for evolution (default: 0.1)\n"
		"-scale choose scale for data: either linear-0, log-1 or custom-2 (default: linear)\n"
		"-h quit with help\n");
		std::cout << std::endl;
}

/// <summary>
/// We want to handle files so let's make the c-way input a string
/// </summary>
/// <param name="argc"> number of main input arguments </param>
/// <param name="argv"> main input arguments </param>
/// <returns></returns>
std::vector<std::string> change_input_to_vec_of_str(int argc, char **argv)
{
	// -1 because first is the name of the file
	NO_OVERFLOW(std::vector<std::string> tmp(argc - 1, "");)
	for (int i = 0; i < argc - 1; i++)
	{
		tmp[i] = argv[i + 1];
	}
	return tmp;
}

// ----------------------------------------------------------------------------- USER INTERFACE FOR OPTIONS -----------------------------------------------------------------------------

/// <summary>
/// Find a given option in a vector of string given from cmd parser
/// </summary>
/// <param name="vec">vector of strings from cmd</param>
/// <param name="option">the option that we seek</param>
/// <returns>value for given option if exists, if not an empty string</returns>
string user_interface::getCmdOption(const v_1d<string> &vec, string option) const
{
	if (auto itr = std::find(vec.begin(), vec.end(), option); itr != vec.end() && ++itr != vec.end())
		return *itr;
	return std::string();
}

/// <summary>
///	Sets an option, if it fails or it lower than zero when geq_0 is set, the dafault value is used
/// </summary>
/// <param name="value">value of an option to be set by reference</param>
/// <param name="argv">arguments list</param>
/// <param name="choosen_option">choosen option string</param>
/// <param name="geq_0"></param>
template <typename T>
void user_interface::set_option(T &value, const v_1d<string> &argv, string choosen_option, bool geq_0)
{
	if (std::string option = this->getCmdOption(argv, choosen_option); option != "")
		value = static_cast<T>(stod(option)); // set value to an option
	if (geq_0 && value < 0)					  // if the variable shall be bigger equal 0
		this->set_default_msg(value, choosen_option.substr(1),
							  choosen_option + " cannot be negative\n", isingUI::table);
}

/// <summary>
/// sets the default message for a given option and tells you why
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="value">sets that value to default</param>
/// <param name="option">option corresponding to the value</param>
/// <param name="message">message for the user</param>
/// <param name="map">map with default values</param>
template <typename T>
void user_interface::set_default_msg(T &value, string option, string message, const unordered_map<string, string> &map) const
{
	stout << message;			// print warning
	std::string value_str = ""; // we will set this to value
	if (auto it = map.find(option); it != map.end())
	{
		value_str = it->second; // if in table - we take the enum
	}
	value = static_cast<T>(stod(value_str));
}

// ----------------------------------------------------------------------------- ISING MODEL -----------------------------------------------------------------------------

/// <summary>
/// prints all options from the user interface (if not chosen the default is taken)
/// </summary>
void isingUI::ui::printAllOptions() const
{
	stout << "------------------------------CHOSEN MODEL:" << std::endl;
	#ifdef HEISENBERG
		std::cout << "HEISENBERG:\n\t\t" << "H = \u03A3_i J_i(\u03C3^x_i \u03C3^x_i+1 + \u03C3^y_i \u03C3^y_i+1) + g_i \u03C3^z_i\u03C3^z_i+1 + h_i \u03C3^x_i" << std::endl << std::endl;
	#else
		std::cout << "ISING:\n\t\t" << "H = \u03A3_i J_i \u03C3^z_i \u03C3^z_i+1 + g_i \u03C3^x_i + h_i \u03C3^x_i" << std::endl << std::endl;
	#endif
	std::cout << "J_i \u03B5 [J - J0, J + J0]" << std::endl;
	std::cout << "g_i \u03B5 [g - g0, g + g0]" << std::endl;
	std::cout << "h_i \u03B5 [h - w, h + w]" << std::endl;
	stout << "------------------------------CHOSEN OPTIONS:" << std::endl;
	std::string opName = std::get<0>(IsingModel_disorder::opName(this->op, this->site));
	stout << "DIR = " << this->saving_dir << std::endl
		  << "model = " << (this->m ? "symmetric" : "disordered") << std::endl
		  << "BC = " << (this->boundary_conditions ? "OBC" : "PBC") << std::endl
		  << "L  = " << this->L << std::endl
		  << "Ls = " << this->Ls << std::endl
		  << "Ln = " << this->Ln << std::endl
		  << "J  = " << this->J << std::endl
		  << "h  = " << this->h << std::endl
		  << "hs = " << this->hs << std::endl
		  << "hn = " << this->hn << std::endl
		  << "g  = " << this->g << std::endl
		  << "gs = " << this->gs << std::endl
		  << "gn = " << this->gn << std::endl
		  << "thread_num = " << this->thread_number << std::endl
		  << "site = " << this->site << std::endl
		  << "operator = " << opName << std::endl
		  << "bucket size = " << this->mu << std::endl
		  << "time step = " << this->dt << std::endl
		  << "boolean value = " << this->ch << std::endl
		  << "scale = " << (this->scale == 1 ? "log" : "linear") << std::endl;

	if (this->m == 0)
		stout << "J0  = " << this->J0 << std::endl
			  << "w   = " << this->w << std::endl
			  << "ws  = " << this->ws << std::endl
			  << "wn  = " << this->wn << std::endl
			  << "g0  = " << this->g0 << std::endl
			  << "g0s = " << this->g0s << std::endl
			  << "g0n = " << this->g0n << std::endl
			  << "realisations = " << this->realisations << std::endl;

	if (this->m == 1)
		stout << "k-sector = " << 2 * this->symmetries.k_sym / this->L << "*pi" << std::endl
			  << "p-sector = " << (this->symmetries.p_sym ? 1 : -1) << std::endl
			  << "x-sector = " << (this->symmetries.x_sym ? 1 : -1) << std::endl;
	stout << "---------------------------------------------------------------------------------\n\n";
	#ifdef PRINT_HELP
		print_help();
		stout << "---------------------------------------------------------------------------------\n\n";
	#endif
}
// ----------------------------------------------------------------------------- Connected with the parser
/// <summary>
/// Setting the default parameters for the Ising model
/// </summary>
void isingUI::ui::set_default()
{
	using namespace isingUI;
	this->saving_dir = "." + std::string(kPathSeparator) + "results" + std::string(kPathSeparator); // directory for the result files to be saved into
	this->L = 4;
	this->Ls = 1;
	this->Ln = 1;

	this->J = 1.0;
	this->J0 = 0.0;

	this->h = 0.0;
	this->hs = 0.1;
	this->hn = 1;

	this->w = 0.01;
	this->ws = 0.0;
	this->wn = 1;

	this->g = 1.0;
	this->gs = 0.1;
	this->gn = 1;

	this->g0 = 0;
	this->g0s = 0.0;
	this->g0n = 1;

	this->symmetries.k_sym = 0;
	this->symmetries.p_sym = 1;
	this->symmetries.x_sym = 1;

	this->realisations = 100;
	this->site = 0;
	this->op = 0;
	this->fun = INT_MAX;
	this->mu = 5;

	this->boundary_conditions = 0;
	this->m = 0;
	this->p = true;
	this->thread_number = 1;

	this->ch = false;
	this->dt = 0.1;
	this->scale = 0;

	this->seed = static_cast<long unsigned int>(87178291199L);
	this->jobid = 0;
}

// ------------------------------------- CONSTURCTORS

/// <summary>
/// Ising model user interface constructor, we pass there the command lines arguments from main
/// </summary>
/// <param name="argc"> number of arguments </param>
/// <param name="argv"> arguments list </param>
isingUI::ui::ui(int argc, char **argv)
{
	auto input = change_input_to_vec_of_str(argc, argv);			// change standard input to vec of strings
	input = std::vector<std::string>(input.begin()++, input.end()); // skip the first element which is the name of file
	// plog::init(plog::info, "log.txt");														// initialize logger
	if (std::string option = this->getCmdOption(input, "-f"); option != "")
	{
		input = this->parseInputFile(option); // parse input from file
	}
	this->parseModel((int)input.size(), input); // parse input from CMD directly
}

/// <summary>
/// Function that tells how does the parser work
/// </summary>
void isingUI::ui::exit_with_help() const
{
	print_help();
	std::exit(1);
}

// ------------------------------------- PARSERS

/// <summary>
/// The parser for the Transverse Field Ising model
/// </summary>
/// <param name="argc"> number of arguments </param>
/// <param name="argv"> list of arguments </param>
void isingUI::ui::parseModel(int argc, std::vector<std::string> argv)
{
	using namespace isingUI;
	// SET DEFAULT VALUES
	this->set_default(); // setting default at the very beginning

	std::string choosen_option = "";																// current choosen option
	std::string str_model = "disorder" + std::string(kPathSeparator); // folder for current model
	//---------- SIMULATION PARAMETERS
	// spin coupling
	choosen_option = "-J";
	this->set_option(this->J, argv, choosen_option);

	// spin coupling disorder
	choosen_option = "-J0";
	this->set_option(this->J0, argv, choosen_option);

	// transverse field
	choosen_option = "-g";
	this->set_option(this->g, argv, choosen_option);
	choosen_option = "-gs";
	this->set_option(this->gs, argv, choosen_option, false);
	choosen_option = "-gn";
	this->set_option(this->gn, argv, choosen_option);

	// transverse field disorder
	choosen_option = "-g0";
	this->set_option(this->g0, argv, choosen_option);
	choosen_option = "-g0s";
	this->set_option(this->g0s, argv, choosen_option);
	choosen_option = "-g0n";
	this->set_option(this->g0n, argv, choosen_option);

	// perpendicular field
	choosen_option = "-h";
	this->set_option(this->h, argv, choosen_option);
	choosen_option = "-hs";
	this->set_option(this->hs, argv, choosen_option, false);
	choosen_option = "-hn";
	this->set_option(this->hn, argv, choosen_option);

	// perpendicular field disorder
	choosen_option = "-w";
	this->set_option(this->w, argv, choosen_option);
	choosen_option = "-ws";
	this->set_option(this->ws, argv, choosen_option, false);
	choosen_option = "-wn";
	this->set_option(this->wn, argv, choosen_option);

	// chain length
	choosen_option = "-L";
	this->set_option(this->L, argv, choosen_option);
	choosen_option = "-Ls";
	this->set_option(this->Ls, argv, choosen_option, false);
	choosen_option = "-Ln";
	this->set_option(this->Ln, argv, choosen_option);

	// boundary condition
	choosen_option = "-b";
	this->set_option(this->boundary_conditions, argv, choosen_option);
	if (this->boundary_conditions > 2)
		this->set_default_msg(this->boundary_conditions, choosen_option.substr(1),
							  "max boundary condition is 2", table);

	// choose site
	choosen_option = "-s";
	this->set_option(this->site, argv, choosen_option);

	// choose operator
	choosen_option = "-op";
	this->set_option(this->op, argv, choosen_option);

	// time step and boolean value and scale
	choosen_option = "-dt";
	this->set_option(this->dt, argv, choosen_option);
	choosen_option = "-ch";
	this->set_option(this->ch, argv, choosen_option);
	choosen_option = "-scale";
	this->set_option(this->scale, argv, choosen_option);

	// disorder
	choosen_option = "-jobid";
	this->set_option(this->jobid, argv, choosen_option);
	choosen_option = "-seed";
	this->set_option(this->seed, argv, choosen_option);

	// model
	choosen_option = "-m";
	this->set_option(this->m, argv, choosen_option);
	if (this->m > 1)
		this->set_default_msg(this->m, choosen_option.substr(1),
							  "max model number is 1", table);

	// choose function
	choosen_option = "-fun";
	this->set_option(this->fun, argv, choosen_option, false);
	if (this->fun < 2)
		this->m = 0;

	// buckets
	choosen_option = "-mu";
	this->set_option(this->mu, argv, choosen_option);
	choosen_option = "-r";
	this->set_option(this->realisations, argv, choosen_option);

	// symmetries
	choosen_option = "-k";
	this->set_option(this->symmetries.k_sym, argv, choosen_option);
	if (this->symmetries.k_sym >= this->L)
		this->set_default_msg(this->symmetries.k_sym, choosen_option.substr(1),
							  "max k sector is L = " + std::to_string(this->L), table);
	choosen_option = "-p";
	this->set_option(this->symmetries.p_sym, argv, choosen_option, false);
	choosen_option = "-x";
	this->set_option(this->symmetries.x_sym, argv, choosen_option, false);

	// thread number
	choosen_option = "-th";
	this->set_option(this->thread_number, argv, choosen_option);
	this->thread_number /= outer_threads;
	if (this->thread_number > std::thread::hardware_concurrency())
		this->set_default_msg(this->thread_number, choosen_option.substr(1),
							  "Wrong number of threads\n", table);
	if(this->thread_number < 0) this->thread_number = 1;
	omp_set_num_threads(this->thread_number);
	num_of_threads = this->thread_number;

	// get help
	choosen_option = "-help";
	if (std::string option = this->getCmdOption(argv, choosen_option); option != "")
		exit_with_help();

	// make folder based on a model
	switch (this->m)
	{
	case 0:
		str_model = "disorder" + std::string(kPathSeparator);
		break;
	case 1:
		str_model = "symmetries" + std::string(kPathSeparator);
		break;
	default:
		str_model = "disorder" + std::string(kPathSeparator);
		break;
	}
	// make boundary condition folder
	switch (this->boundary_conditions)
	{
	case 0:
		str_model += "PBC" + std::string(kPathSeparator);
		break;
	case 1:
		str_model += "OBC" + std::string(kPathSeparator);
		break;
	default:
		str_model += "PBC" + std::string(kPathSeparator);
		break;
	}

	#ifdef HEISENBERG
		std::string folder = saving_dir + "HEISENBERG" + kPSep + str_model;
	#else
		std::string folder = saving_dir + str_model;
	#endif
	if (!argv[argc - 1].empty() && argc % 2 != 0) {
		// only if the last command is non-even
		folder = argv[argc - 1] + str_model;
		if (fs::create_directories(folder) || fs::is_directory(folder)) // creating the directory for saving the files with results
			this->saving_dir = folder;									// if can create dir this is is
	} else {
		if (fs::create_directories(folder) || fs::is_directory(folder)) // creating the directory for saving the files with results
			this->saving_dir = folder;									// if can create dir this is is
	}

	std::cout << " - - - - - - MAKING ISING INTERFACE AND USING OUTER THREADS : "
			  << outer_threads << " - - - - - - " << endl; // setting the number of threads to be used with omp

	std::cout << " - - - - - - MAKING ISING INTERFACE AND USING INNER THREADS : "
			  << thread_number << " - - - - - - " << endl; // setting the number of threads to be used with omp

	omp_set_num_threads(this->thread_number);
	return;
}

/// <summary>
/// If the commands are given from file, we must treat them the same as arguments
/// </summary>
/// <param name="filename"> the name of the file that contains the command line </param>
/// <returns></returns>
std::vector<std::string> user_interface::parseInputFile(std::string filename) const
{
	std::vector<std::string> commands(1, "");
	ifstream inputFile(filename);
	std::string line = "";
	if (!inputFile.is_open())
		stout << "Cannot open a file " + filename + " that I could parse. All parameters are default. Sorry :c \n";
	else
	{
		if (std::getline(inputFile, line))
		{
			commands = split_str(line, " "); // saving lines to out vector if it can be done, then the parser shall treat them normally
		}
	}
	return std::vector<std::string>(commands.begin(), commands.end());
}



/*
	double W = 5.5;
	arma::mat energies;
	std::string base = "/Users/rafal.swietek/Downloads/rafal_test_sff/jan_data_3D_Anderson/";
	energies.load(arma::hdf5_name(base + "eigvals_L_16_dim_3_pbc_True_disorder_uniform_ham_type_anderson_t_-1.0_W_0.0_dW_" + to_string_prec(W, 2) + ".hdf5", "Eigenvalues"));
	std::cout << energies.n_cols << "\t" << energies.n_rows << std::endl;
	arma::mat sff_jan1, sff_jan2, sff_jan3;
	bool state1 = sff_jan1.load(arma::hdf5_name(base + "eigvals_L_16_dim_3_pbc_True_disorder_uniform_ham_type_anderson_t_-1.0_W_0.0_dW_" + to_string_prec(W, 2) + ".hdf5", "SFF_spectra_eta_0.1000_filter_gaussian"));
	bool state2 = sff_jan2.load(arma::hdf5_name(base + "eigvals_L_16_dim_3_pbc_True_disorder_uniform_ham_type_anderson_t_-1.0_W_0.0_dW_" + to_string_prec(W, 2) + ".hdf5", "SFF_spectra_eta_0.3000_filter_gaussian"));
	bool state3 = sff_jan3.load(arma::hdf5_name(base + "eigvals_L_16_dim_3_pbc_True_disorder_uniform_ham_type_anderson_t_-1.0_W_0.0_dW_" + to_string_prec(W, 2) + ".hdf5", "SFF_spectra_eta_0.5000_filter_gaussian"));
	stout << state1 << "\t" << state2 << "\t" << state3 << std::endl;
	printSeparated(std::cout, "\t", 12, true, sff_jan1.n_cols, sff_jan1.n_rows, sff_jan2.n_cols, sff_jan2.n_rows, sff_jan3.n_cols, sff_jan3.n_rows);
	std::cout << "Data loaded" << std::endl;

	auto times = arma::logspace(log10(1.0 / 4096.0), log10(5*two_pi), 2000);
	arma::vec sff1(times.size(), arma::fill::zeros);
	arma::vec sff2 = sff1, sff3 = sff1;
	double Z1 = 0, Z2 = 0, Z3 = 0;
	int n_bins = 1 + long(3.322 * log(energies.n_rows));
	arma::vec histE(n_bins, arma::fill::zeros);
	arma::vec histE_unfolded(n_bins, arma::fill::zeros);
	
	for(int r = 0; r < energies.n_cols; r++){
		arma::vec eigenvalues = energies.col(r);
		histE = normalise_dist(arma::conv_to<arma::vec>::from(
            arma::hist(eigenvalues, n_bins)
            ), arma::min(eigenvalues), arma::max(eigenvalues));
		
		statistics::unfolding(eigenvalues);
		histE_unfolded = normalise_dist(arma::conv_to<arma::vec>::from(
            arma::hist(eigenvalues, n_bins)
            ), arma::min(eigenvalues), arma::max(eigenvalues));
		
		auto [sff_r1, Z_r1] = statistics::spectral_form_factor(eigenvalues, times, 0.1);
		auto [sff_r2, Z_r2] = statistics::spectral_form_factor(eigenvalues, times, 0.3);
		auto [sff_r3, Z_r3] = statistics::spectral_form_factor(eigenvalues, times, 0.5);
		sff1 += sff_r1; Z1 += Z_r1;
		sff2 += sff_r2; Z2 += Z_r2;
		sff3 += sff_r3; Z3 += Z_r3;
	}
	//histE_unfolded /= double(energies.n_cols);
	//histE		   /= double(energies.n_cols);
	sff1 /= Z1;
	sff2 /= Z2;
	sff3 /= Z3;

	std::string out_dir = this->saving_dir + "SFF_TESTS" + kPSep;
	createDirs(out_dir);
	std::ofstream sff_file, hist_file;
	openFile(sff_file, out_dir + "sff.dat", std::ios::out);
	openFile(hist_file, out_dir + "hist.dat", std::ios::out);
	arma::vec xvals = arma::linspace(0.0, 1.0, n_bins);
	for(int i = 0; i < n_bins; i++)
		printSeparated(hist_file, "\t", 16, true, xvals(i), histE(i), histE_unfolded(i));
	hist_file.close();
	for(int j = 0; j<times.size(); j++)
		printSeparated(sff_file, "\t", 16, true, times(j), sff1(j), sff2(j), sff3(j));//, sff_jan1(j), sff_jan2(j), sff_jan3(j));
	sff_file.close();

	exit(1);
	*/

/*
	COMMUTATORS
	
	auto h_scale = arma::logspace(-5, 0, 200); arma::vec zeros(1, arma::fill::zeros);
	h_scale = arma::join_cols(zeros, h_scale);
	std::string info =IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, {"h"});
	std::ofstream file; 
	openFile(file, this->saving_dir + "commutator_tfim_lioms" + info + ".dat", std::ios::out);
	for(double hx : h_scale){
		this->h = hx;
		auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
		auto H = alfa->get_hamiltonian();
		printSeparated(file, "\t", 12, false, this->h);
		printSeparated(std::cout, "\t", 12, false, this->h);
		for(this->site = 0; this->site < this->L; this->site++){
			arma::sp_cx_mat op = alfa->chooseOperator(this->op, this->site);
			arma::sp_cx_mat comm = H * op - op * H;
			double norm = 0;
			for(const cpx& x : comm)
				norm += abs(x);
			norm /= double(comm.n_nonzero == 0? 1.0 : comm.n_nonzero);
			printSeparated(file, "\t", 12, false, norm);
			printSeparated(std::cout, "\t", 12, false, norm);
			file.flush();
		}
		file << std::endl;
		std::cout << std::endl;
	}
	exit(123);
 */