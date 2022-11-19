#include "include/user_interface.h"
// set externs
std::uniform_real_distribution<> theta	= std::uniform_real_distribution<>(0.0, pi);
std::uniform_real_distribution<> fi		= std::uniform_real_distribution<>(0.0, pi);
int outer_threads = 1;
int anderson_dim = 3;
//---------------------------------------------------------------------------------------------------------------- UI main
void isingUI::ui::make_sim()
{	
	#if defined(MY_MAC)
		this->seed = static_cast<long unsigned int>(time(0));
	#endif
	my_gen = randomGen(this->seed);
	printAllOptions();

	//check_symmetry_rotation(); return;
	//compare_energies();	return;
	
	//auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
	//auto alfa = std::make_unique<IsingModel_disorder>(this->L, 1, 0, 1, 0, 0, 0, this->boundary_conditions);
	//auto H = alfa->get_hamiltonian();
	//alfa->diagonalization();
	//std::cout << alfa->get_eigenEnergy(0) << std::endl;
	//u64 idx = 286;// 4. / 3. * (ULLPOW(this->L) - 1);
	//std::ifstream filee;
	//auto data = readFromFile(filee, this->saving_dir + "random_vec.dat");
	//arma::cx_vec random = cast_cx_vec(data[0]);
	//lanczos::Lanczos lancz(H, lanczosParams(this->mu, 1, false, false), random);
	//lancz.diagonalization();
	//auto T = lancz.get_lanczos_matrix();
	//std::cout << lancz.get_eigenvalues()(0) << std::endl;
	//std::ofstream file;
	//std::cout << T << std::endl;
	//file.open(this->saving_dir + "lanczos_matrix.dat");
	//file << T;
	//file.close();
	//return;

	auto kernel = [&](int k, int p, int x)
	{
		this->symmetries.k_sym = k;
		this->symmetries.p_sym = p;
		this->symmetries.x_sym = x;
		this->eigenstate_entropy();
	};
	loopSymmetrySectors(kernel);
	return;

	clk::time_point start = std::chrono::system_clock::now();
	auto L_list = this->get_params_array(Ising_params::L);
	auto J_list = this->get_params_array(Ising_params::J);
	auto g_list = this->get_params_array(Ising_params::g);
	auto h_list = this->get_params_array(Ising_params::h);
	auto w_list = this->get_params_array(Ising_params::w);

	Ising_params var;
	switch(this->op){
		case 0: var = Ising_params::L;	break;
		case 1: var = Ising_params::J;	break;
		case 2: var = Ising_params::g;	break;
		case 3: var = Ising_params::h;	break;
		case 4: var = Ising_params::w;	break;
		case 5: var = Ising_params::k;	break;
		default:
			var = Ising_params::g;
	}

	//for (auto& Ll : L_list){
	//	this->L = Ll;
	//	this->site = this->L / 2;
	//	generate_statistic_map(Ising_params::g); 
	//	//thouless_times(var);
	//}; return;

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
		calculate_localisation_length();
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
		thouless_times(var);
		break;
	case 9:
		calculate_statistics();
		break;
	case 10:
		eigenstate_entropy();
		break;
	default:
		for (auto& system_size : L_list){
			for (auto& gx : g_list){
				for (auto& hx : h_list){
					for(auto& Jx : J_list){
						for(auto& wx : w_list){
							this->L = system_size;
							this->g = gx;
							this->h = hx;
							this->J = Jx;
							this->w = wx;
							this->site = this->L / 2.;
							
							const auto start_loop = std::chrono::system_clock::now();
							stout << " - - START NEW ITERATION AT : " << tim_s(start) << " s;\t\t par = "; // simulation end
							printSeparated(std::cout, "\t", 16, true, this->L, this->J, this->g, this->h, this->w);
							//calculate_localisation_length(); continue;
							calculate_statistics(); continue;

							combine_spectrals(); continue;

							average_SFF(); continue;

							analyze_spectra(); continue;
							spectral_form_factor(); continue;
						}}}}}
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
			eigenvalues = arma::real(E);
			std::sort(eigenvalues.begin(), eigenvalues.end());
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
		loaded = eigenvalues.load(arma::hdf5_name(name + _suffix + ".hdf5", "eigenvalues/dataset"));
		if(!loaded)
			loaded = eigenvalues.load(arma::hdf5_name(name + ".hdf5", "eigenvalues/" + _suffix));
	}
	if(!loaded){
		bool status = false;
		#if !defined(ANDERSON) && !defined(HEISENBRG)
			status = true;
		#endif
		if(status && alfa.g == 0){
			auto H = alfa.get_hamiltonian();
			const u64 N = alfa.get_hilbert_size();
			arma::cx_vec E(N);
			for(int i = 0; i < N; i++)
				E(i) = H(i,i);
			eigenvalues = real(E);
			sort(eigenvalues.begin(), eigenvalues.end());
		} 
		else if(status && alfa.J == 0.0){
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

//-------------------------------------------------------------------------- GENERAL 
//<! generate random product state (random orientation of spins on the bloch sphere)
arma::cx_vec isingUI::ui::random_product_state(int system_size)
{
	auto the = my_gen.random_uni<double>(0.0, pi);
	arma::cx_vec init_state = std::cos(the / 2.) * up
		+ std::exp(im * my_gen.random_uni<double>(0.0, pi)) * std::sin(the / 2.) * down;
	for (int j = 1; j < system_size; j++)
	{
		the = my_gen.random_uni<double>(0.0, pi);
		init_state = arma::kron(init_state, std::cos(the / 2.) * up
			+ std::exp(im * my_gen.random_uni<double>(0.0, pi)) * std::sin(the / 2.) * down);
	}
	return init_state;
}

//<! generate initial state given by input (user control) or -op flag (default): random, FM, AFM, ...
arma::cx_vec isingUI::ui::set_init_state(size_t N, int choose)
{
	if(choose < 0)
		choose = this->op;
	arma::cx_vec init_state(N, arma::fill::zeros);
	switch (choose) {
	case 0: // random product state
	{
		init_state = this->random_product_state(this->L); 
		break;
	}
	case 1: // ferromagnetically polarised
	{
		u64 idx = (ULLPOW(this->L)) - 1;
		init_state(idx) = cpx(1.0, 0.0); // 1111111
		break;
	}
	case 2: // anti-ferromagnetically polarised: 1010 + 0101
	{
		u64 idx = ((ULLPOW(this->L)) - 1) / 3;
		init_state(						idx) = cpx(1.0, 0.0); // 10101010
		init_state((ULLPOW(this->L)) -	idx) = cpx(1.0, 0.0); // 01010101
		break;
	}
	default:
		init_state = random_product_state(this->L); 
	}
	return arma::normalise(init_state);
}
//<! get parameter array
arma::vec isingUI::ui::get_params_array(Ising_params par){
	
	// general
	auto J_list = arma::linspace(this->J, this->J + this->Js * (this->Jn - 1), this->Jn);
	auto g_list = arma::linspace(this->g, this->g + this->gs * (this->gn - 1), this->gn);
	auto h_list = arma::linspace(this->h, this->h + this->hs * (this->hn - 1), this->hn);
	auto L_list = arma::linspace(this->L, this->L + this->Ls * (this->Ln - 1), this->Ln);
	
	// disorder
	auto J0_list = arma::linspace(this->J0, this->J0 + this->J0s * (this->J0n - 1), this->J0n);
	auto g0_list = arma::linspace(this->g0, this->g0 + this->g0s * (this->g0n - 1), this->g0n);
	auto w_list = arma::linspace(this->w, this->w + this->ws * (this->wn - 1), this->wn);

	arma::vec result;
	switch(par){
		case Ising_params::L : 	result = L_list;	break;
		case Ising_params::J : 	result = J_list;	break;
		case Ising_params::J0 : result = J0_list;	break;
		case Ising_params::g : 	result = g_list;	break;
		case Ising_params::g0 : result = g0_list;	break;
		case Ising_params::h : 	result = h_list;	break;
		case Ising_params::w : 	result = w_list;	break;
		default:
			std::cout << "No option found, choosing g_array" << std::endl;
			result = g_list;
	}
	return result;
}

//-------------------------------------------------------------------- COMPARING SYMMETRIC TO DISORDERED RESULTS
void isingUI::ui::compare_energies()
{
	// handle disorder
	auto Hamil = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, boundary_conditions);
	Hamil->diagonalization();
	const arma::vec E_dis = Hamil->get_eigenvalues();
	Hamil.release();
	// here we push back the energies from symmetries
	std::vector<double> E_sym = v_1d<double>();
	std::vector<std::string> symms = v_1d<std::string>();
	// go for each symmetry sector
	const int x_max = (this->h != 0) ? 0 : 1;
	const int k_end = (this->boundary_conditions) ? 1 : this->L;
	for (int k = 0; k < k_end; k++)
	{
		if (k == 0 || k == this->L / 2.)
		{
			for (int p = 0; p <= 1; p++)
			{
				// if the spin flip is unaviable we just use 1
				for (int x = 0; x <= x_max; x++)
				{
					auto ham = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h, k, p, x, boundary_conditions);
					ham->diagonalization(false);
					arma::vec t = ham->get_eigenvalues();
					E_sym.insert(E_sym.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));
					v_1d<std::string> temp_str = v_1d<std::string>(t.size(), "k=" + std::to_string(k) + ",x=" + to_string(x) + ",p=" + to_string(p));
					symms.insert(symms.end(), std::make_move_iterator(temp_str.begin()), std::make_move_iterator(temp_str.end()));
				}
			}
		}
		else
		{
			for (int x = 0; x <= x_max; x++)
			{
				auto ham = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h, k, 1, x, boundary_conditions);
				ham->diagonalization(false);
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
	printSeparated(std::cout, "\t", 20, true, "symmetry sector", "Energy sym", "Energy total", "difference");
	for (int k = 0; k < E_dis.size(); k++)
		printSeparated(std::cout, "\t", 20, true, symms[k], E_sym[k], E_dis(k), E_sym[k] - E_dis(k));
	
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
	const u64 N = alfa1->get_hilbert_size();
	std::cout << std::endl;
	printSeparated(std::cout, "\t", 12, true, "L_A", "w = 1e-4", "w = 1e-3", "w = 1e-2", "k = 0", "k = 1");
	for (int i = 3; i < this->L - 2; i++)
	{
		this->mu = N > 3000 ? 500 : 0.25 * N;
		u64 E_min = alfa1->E_av_idx - this->mu / 2.;
		u64 E_max = alfa1->E_av_idx + this->mu / 2.;
		double entropy_dis1 = 0.0, entropy_dis2 = 0.0, entropy_dis3 = 0.0;
		for (long k = E_min; k < E_max; k++)
		{
			auto state = alfa1->get_state_in_full_Hilbert(arma::cx_vec(alfa1->get_eigenState(k), arma::vec(N, arma::fill::zeros)));
			entropy_dis1 += entropy::vonNeumann(state, i, alfa1->L);
			state = alfa2->get_state_in_full_Hilbert(arma::cx_vec(alfa2->get_eigenState(k), arma::vec(N, arma::fill::zeros)));
			entropy_dis2 += entropy::vonNeumann(state, i, alfa2->L);
			state = alfa3->get_state_in_full_Hilbert(arma::cx_vec(alfa3->get_eigenState(k), arma::vec(N, arma::fill::zeros)));
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
		printSeparated(file, "\t", 16, true, "#cores", "chain length", "N", "with eigenvec 'dc'", "only eigenvalues", "time in seconds");
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
void isingUI::ui::check_symmetry_rotation(){
	auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, 0, this->boundary_conditions);
	arma::sp_mat H0 = alfa->get_hamiltonian();

	auto lambda = [this](int k, int p, int x, arma::sp_cx_mat& H){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 k, p, x, this->boundary_conditions);
		auto U = alfa->symmetryRotation();
		arma::sp_cx_mat Hsym = alfa->get_hamiltonian();
		H += U * Hsym * U.t();
	};

	arma::sp_cx_mat H(H0.n_rows, H0.n_cols);

	const int x_max = (abs(this->h) > 0) ? 0 : 1;
	for (int k = 0; k < this->L; k++) {
		if (k == 0 || k == this->L / 2.) {
			for (int p = 0; p <= 1; p++)
				for (int x = 0; x <= x_max; x++)
					lambda(k, p, x, H);
		}
		else {
			for (int x = 0; x <= x_max; x++)
					lambda(k, 0, x, H);
		}
	}
	arma::sp_mat HH = arma::real(H);
	auto N = H0.n_cols;
	arma::sp_cx_mat res = arma::sp_cx_mat(HH - H0, arma::imag(H));
	printSeparated(std::cout, "\t", 32, true, "index i", "index j", "difference", "original hamil", "symmetry hamil");
	cpx x = 0;
	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			cpx val = res(i, j);
			if(abs(H0(i, j)) > 1e-15){
				x += val;
				printSeparated(std::cout, "\t", 32, true, i, j, res(i, j), H0(i, j), H(i, j));
			}
		}
	}
	printSeparated(std::cout, "\t", 32, true, "Sum of suspicious elements: ", x);
}
//--------------------------------------------------------------------- SPECTRAL PROPERTIES and TIME EVOLUTION

//<! calculate all spectral quantities: time evolution, response function,
//<! integrated spectral function and spectral form factor with folded eigenvalues
void isingUI::ui::calculate_spectrals()
{
	const clk::time_point start = std::chrono::system_clock::now();
	//---------------- DIRECTORIES & INFO
	std::string info = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym) 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w);

	auto [opName_tup, subdir] = IsingModel_disorder::opName(this->op, this->site);
	std::string opName = opName_tup;
	std::string timeDir = this->saving_dir + "timeEvolution" + kPSep + subdir + kPSep;
	std::string specDir = this->saving_dir + "ResponseFunction" + kPSep + subdir + kPSep;
	std::string intDir = this->saving_dir + "IntegratedResponseFunction" + kPSep + subdir + kPSep;
	std::string specDir_der = this->saving_dir + "IntegratedResponseFunction" + kPSep + "DERIVATIVE" + kPSep + subdir + kPSep;
	createDirs(timeDir, specDir, intDir, specDir_der);

	//---------------- PREAMBLE
	int counter = 0;
	double 	 
			 AGP = 0.0,	//<! adiabatic gauge potential
		typ_susc = 0.0,	//<! typical fidelity susceptibility
			susc = 0.0,	//<! fidelity susceptibility
	   		  wH = 0.0,	//<! mean level spacing
	   	  wH_typ = 0.0;	//<! typical level spacing
	
	const double chi = 0.341345;
	size_t N = ULLPOW(this->L);
	const double wooH = sqrt(this->L) / (chi * N) * sqrt(this->J * this->J + this->h * this->h + this->g * this->g
												 + ( this->m? 0.0 : (this->w * this->w + this->g0 * this->g0 + this->J0 * this->J0) / 3. ));
	double tH = 1. / wooH;
	int num_of_points = 1000;
	int time_end = (int)std::ceil(std::log10(1.75 * tH));
	time_end = (time_end / std::log10(tH) < 2.0) ? time_end + 1 : time_end;
	auto times = arma::logspace(-2, time_end, num_of_points);
	auto omegas = arma::logspace(-time_end, 2, num_of_points);
	arma::vec omega_spec = omegas; omega_spec.shed_row(omega_spec.size() - 1);

	arma::vec opEvol(times.size(), arma::fill::zeros);
	arma::vec opIntSpec(omegas.size(), arma::fill::zeros);
	std::vector<u64> map;
	
	//---------------- SET KERNEL
	double LTA = 0;
	arma::sp_cx_mat op;
	auto start_loop = std::chrono::system_clock::now();
	auto kernel = [&](
		auto& alfa, int r
		){
		r += this->jobid;
		std::string tdir_realisation = timeDir + "realisation=" + std::to_string(r) + kPSep;
		std::string intdir_realisation = intDir + "realisation=" + std::to_string(r) + kPSep;
		std::string specdir_realisation = specDir + "realisation=" + std::to_string(r) + kPSep;
		std::string specdir_real_mat_elem = specdir_realisation + "mat_elem" + kPSep;
		createDirs(tdir_realisation, intdir_realisation, specdir_realisation, specdir_real_mat_elem);
		
		stout << "\t\t	--> finished diagonalizing for " << info << " - in time : " << tim_s(start_loop) << "s" << std::endl;
		auto U = alfa.get_eigenvectors();
		arma::vec E = alfa.get_eigenvalues();
		
		stout << "\t\t	--> got eigenvectors for " << info << " - in time : " << tim_s(start_loop) << "s" << std::endl;
		op = alfa.chooseOperator(this->op, this->site);
		arma::cx_mat mat_elem = U.t() * op * U;
		normaliseMat(mat_elem);
		stout << "\t\t	--> set matrix elements for " << info << " - in time : " << tim_s(start_loop) << "s" << std::endl;


		N = alfa.get_hilbert_size();
		this->mu = long(0.5 * N);
		long int E_min = alfa.E_av_idx - long(mu / 2);
		long int E_max = alfa.E_av_idx + long(mu / 2);
	
		double wH_local = 0.0;
		double wH_typ_local = 0.0;
		for(int i = 0; i < N; i++){
			if(i >= E_min && i < E_max){
				const double gap = E(i + 1) - E(i);
				wH_local += gap;
				wH_typ_local += std::log(gap);
			}
		}
		wH_local = wH_local / double(this->mu);
		wH += wH_local;
		counter++;

		stout << "\t\t	--> finished statistics for " << info
			  << " realisation: " << r << " - in time : " << tim_s(start_loop) << "s" << std::endl;

		auto [op_tmp, LTA_tmp] = spectrals::autocorrelation_function(mat_elem, E, times);
		save_to_file(tdir_realisation + opName + info + ".dat", times, op_tmp, 1.0 / wH_local, LTA_tmp);
		stout << "\t\t	--> finished time evolution for " << info
			  << " realisation: " << r << " - in time : " << tim_s(start_loop) << "s" << std::endl;
		
		auto res = spectrals::integratedSpectralFunction(mat_elem, E, omegas);
		save_to_file(intdir_realisation + opName + info + ".dat", omegas, res, wH_local, LTA_tmp);
		stout << "\t\t	--> finished integrated spectral function for " << info
			  << " realisation: " << r << " - in time : " << tim_s(start_loop) << "s" << std::endl;

		spectrals::preset_omega set_omega(E, 0.0025 * alfa.L, E(alfa.E_av_idx));
		set_omega.save_matrix_elements(specdir_real_mat_elem + opName + info, mat_elem);

		//auto specfun_r = spectrals::spectralFunction(mat_elem, set_omega, omega_spec);
		//save_to_file(specdir_realisation + opName + info + ".dat", omega_spec, specfun_r, wH_local, LTA_tmp);
		stout << "\t\t	--> finished saving matrix elements for " << info
			  << " realisation: " << r << " - in time : " << tim_s(start_loop) << "s" << std::endl;

		LTA += LTA_tmp;
		opEvol += op_tmp;
		opIntSpec += res;

		auto [AGP_local, typ_susc_local, susc_local] = adiabatics::gauge_potential(mat_elem, E, this->L);
		typ_susc += typ_susc_local;
		AGP += AGP_local;
		susc += susc_local;

		stout << "\t\t	--> finished adiabatic gauge potential for " << info
			  << " realisation: " << r << " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s" << std::endl;

		start_loop = std::chrono::system_clock::now();
	};
	
	//---------------- SIMULATION FOR INPUT MODEL
	
	int M;
	auto prefix_kernel = [&](auto& alfa){
		op = alfa.chooseOperator(this->op, this->site);
		stout << "\n\t\t--> finished generating operator and omega bins for " << info << " - in time : " << tim_s(start) << "s" << std::endl;;	
	};
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		prefix_kernel(*alfa);
		map = alfa->get_mapping();
		average_over_realisations<Ising_params::h>(*alfa, true, kernel);
	} else{
		auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
		prefix_kernel(*alfa);
		#ifdef HEISENBERG
			map = alfa->get_mapping();
		#elif !defined(ANDERSON)
			if(this->g == 0 && this->g0 == 0)
				map = alfa->get_mapping();
		#endif
		average_over_realisations<Ising_params::h>(*alfa, true, kernel);
	}

	std::string dir_agp = this->saving_dir + "AGP" + kPSep + opName + kPSep + "raw_data" + kPSep;
	createDirs(dir_agp);

	opEvol /= double(counter);
	LTA /= double(counter);
	opIntSpec /= double(counter);
	wH /= double(counter);

	std::ofstream file;

	openFile(file, dir_agp + info + "_jobid=" + std::to_string(jobid) + ".dat");
	printSeparated(file, "\t", 25, true, "'adiabatic gauge potential'", 	 	   AGP / double(counter));
	printSeparated(file, "\t", 25, true, "'typical fidelity susceptibility'", typ_susc / double(counter));
	printSeparated(file, "\t", 25, true, "'fidelity susceptibility'", 		 	  susc / double(counter));
	printSeparated(file, "\t", 25, true, "'diagonal normalisation factor'", 	   LTA);

	file.close();
	return;
	std::string filename = opName + info + "_jobid=" + std::to_string(jobid) + ".dat";
	save_to_file(timeDir + filename, times, opEvol, 1.0 / wH, LTA);		//smoothen_data(timeDir, opName + info + ".dat", 10);
	save_to_file(intDir + filename, omegas, opIntSpec, wH, LTA);		//smoothen_data(intDir,  opName + info + ".dat", 10);
	//save_to_file(specDir + opName + info + ".dat", omega_spec, opSpecFun, wH, LTA);		smoothen_data(specDir, opName + info + ".dat");
	
	stout << " - - - - - - FINISHED CALCULATIONS IN : " << tim_s(start) << " seconds - - - - - - \n"
			  << std::endl; // simulation end
}

//<! average separate spectral files for input operator
void isingUI::ui::combine_spectrals(){
	std::string info = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym) 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w);

	auto [opName_tup, subdir] = IsingModel_disorder::opName(this->op, this->site);
	std::string opName = opName_tup;
	std::string timeDir = this->saving_dir + "timeEvolution" + kPSep + subdir + kPSep;
	std::string intDir = this->saving_dir + "IntegratedResponseFunction" + kPSep + subdir + kPSep;
	std::string specDir_der = this->saving_dir + "IntegratedResponseFunction" + kPSep + "DERIVATIVE" + kPSep + subdir + kPSep;
	std::string specDir = this->saving_dir + "ResponseFunction" + kPSep + subdir + kPSep;
	std::string SpecFunDistDir = this->saving_dir + "MatrixElemDistribution" + kPSep + subdir + kPSep;
	std::string HybridDir = this->saving_dir + "Hybrydization" + kPSep + "Distribution" + kPSep;

	std::string dir_agp = this->saving_dir + "AGP" + kPSep + opName + kPSep + "raw_data" + kPSep;
	createDirs(timeDir, intDir, specDir_der, dir_agp);

	int num_of_points = 1000;
	arma::vec opEvol(num_of_points, arma::fill::zeros);
	arma::vec opIntSpec(num_of_points, arma::fill::zeros);
	arma::vec times(num_of_points);
	arma::vec omegas = times;
	arma::vec omegas_spec, mat_elem;
	//---------------- PREAMBLE
	int counter_agp = 0;
	int counter_time = 0;
	int counter_int = 0;
	double wH = 0.0;
	arma::vec agps(4, arma::fill::zeros);
	
	size_t N = 0;
	#ifdef HEISENBERG
		N = binomial(this->L, this->L / 2.);
	#elif defined ANDERSON
		N = this->L * this->L * this->L;
	#else
		N = ULLPOW(this->L);
	#endif

	auto start_loop = std::chrono::system_clock::now();
	auto lambda_average = [&](int realis, double x){
		
		std::string filename = opName + info + ".dat";
		std::ifstream file;

		//-------- AGP
		bool status = openFile(file, dir_agp + info + "_jobid=" + std::to_string(realis) + ".dat");
		if(status){
			int idx = 0;
			std::string line;
			bool isnan = false;
			while(std::getline(file, line)){
				std::istringstream ss(line);
				std::vector<std::string> datarow;
				while(ss){
					std::string element;	ss >> element;
					datarow.push_back(element);
				}
				double value = std::stod(datarow[datarow.size()-2]);
				if(value != value) std::cout << "found nan!!!" << std::endl;
				agps(idx) += value;
				idx++;
			}
			counter_agp++;
			stout << "\t\t	--> finished AGP for " << info
			  << " realisation: " << realis << " - in time : " << tim_s(start_loop) << "s" << std::endl;
		}
		file.close();

		//-------- TIME EVOLUTION
		auto data = readFromFile(file, timeDir + "realisation=" + std::to_string(realis) + kPSep + filename);
		if(!data.empty()){
			times = data[0];
			opEvol += data[1];
			counter_time++;
			stout << "\t\t	--> finished time evolution for " << info
			  << " realisation: " << realis << " - in time : " << tim_s(start_loop) << "s" << std::endl;
		}
		file.close();

		//-------- INTEGRATED SPECTRAL FUNCTION
		data = readFromFile(file, intDir + "realisation=" + std::to_string(realis) + kPSep + filename);
		if(!data.empty()){
			omegas = data[0];
			opIntSpec += data[1];
			wH += data[2][0] / double(0.5 * N);
			counter_int++;
			stout << "\t\t	--> finished integrated spectral function for " << info
			  << " realisation: " << realis << " - in time : " << tim_s(start_loop) << "s" << std::endl;
		}
		file.close();

		//-------- SPECTRAL FUNCTION
		arma::vec omegas, mat_elem_r;
		bool loaded = omegas.load(arma::hdf5_name(specDir + "realisation=" + std::to_string(realis) + kPSep + "mat_elem" + kPSep + opName + info + ".hdf5", "omegas"));
		bool loaded2 = mat_elem_r.load(arma::hdf5_name(specDir + "realisation=" + std::to_string(realis) + kPSep + "mat_elem" + kPSep + opName + info + ".hdf5", "mat_elem"));
		data = readFromFile(file, specDir + "realisation=" + std::to_string(realis) + kPSep + "mat_elem" + kPSep + filename);
		if(loaded && loaded2){
			omegas_spec = arma::join_cols(omegas_spec, omegas);
			mat_elem = arma::join_cols(mat_elem, mat_elem_r);
			stout << "\t\t	--> finished spectral function for " << info
			  << " realisation: " << realis << " - in time : " << tim_s(start_loop) << "s" << std::endl;
		}
	};
	average_over_realisations<Ising_params::J>(false, lambda_average);

	if(counter_agp > 0){
		agps /= double(counter_agp);
		std::ofstream file_out;
		openFile(file_out, dir_agp + info + ".dat");
		printSeparated(file_out, "\t", 25, true, "'adiabatic gauge potential'", 		agps(0));
		printSeparated(file_out, "\t", 25, true, "'typical fidelity susceptibility'",	agps(1));
		printSeparated(file_out, "\t", 25, true, "'fidelity susceptibility'", 			agps(2));
		printSeparated(file_out, "\t", 25, true, "'diagonal normalisation factor'", 	agps(3));
		file_out.close();
	}
	
	std::string filename = opName + info;
	if(counter_int > 0){
		wH /= double(counter_int);
		save_to_file(intDir + filename + ".dat", omegas, opIntSpec / double(counter_int), wH, agps(3));		
		smoothen_data(intDir,  opName + info + ".dat", 10);
		
		std::ifstream file;
		auto data = readFromFile(file, intDir + "smoothed" + kPSep + filename + ".dat");
		auto specFun = non_uniform_derivative(data[0], data[1]);
		arma::vec x = data[0];	x.shed_row(x.size() - 1);
		save_to_file(specDir_der + opName + info + ".dat", x, specFun, wH, agps(3));		
		smoothen_data(specDir_der, opName + info + ".dat");
	}
	if(counter_time > 0){
		save_to_file(timeDir + filename + ".dat", times, opEvol / double(counter_time), 1.0 / wH, agps(3));		
		smoothen_data(timeDir, opName + info + ".dat", 10);
	}

	if(!mat_elem.is_empty() && !omegas_spec.is_empty()){
		//arma::uvec non_zero_elements = arma::find(mat_elem > 1e-34);
		statistics::probability_distribution(SpecFunDistDir, filename, arma::sqrt(mat_elem), -1);
		statistics::probability_distribution(SpecFunDistDir, filename + "_log", 0.5 * arma::log(mat_elem), -1);

		arma::vec hybrid = mat_elem / omegas_spec;
		statistics::probability_distribution(HybridDir, filename, hybrid, -1, arma::mean(hybrid), arma::var(hybrid));
		statistics::probability_distribution(HybridDir, filename + "_log", 0.5 * arma::log(hybrid), -1, arma::mean(hybrid), arma::var(hybrid));

		std::vector<int> numss = {500, 2000, 6000, 10000};
		const long num = (this->realisations == 1)? numss[ (this->L - 12) / 2] : ULLPOW(this->L / 2) * std::sqrt(counter_int);
		spectrals::spectralFunction(omegas_spec, mat_elem, specDir + filename, num);	smoothen_data(specDir, filename + ".dat");
	}
}



//<! calculate evolution of entaglement from initial state chosen by -op.
//<! -s sets the subsystem size, if-s=-1 the L/2 is assumed 
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
		my_gen.reset();
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
	
//<! calculate entaglement in eigenstates and saves to file
//<! -s sets the subsystem size, if-s=0 the L/2 is assumed 
void isingUI::ui::eigenstate_entropy(){

	//---- KERNEL LAMBDA
	std::vector<u64> map;
	std::string dir = this->saving_dir + "Entropy" + kPSep + "Eigenstate" + kPSep;
	createDirs(dir);
	int LA = this->L / 2;
	size_t N = 0;
	arma::vec energies, entropies;
	
	std::string info = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym) 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w);
	int counter = 0;
	std::string filename = info + "_subsize=" + std::to_string(LA);
	auto kernel = [&](auto& alfa, int realis)
	{
		realis += this->jobid;

		std::string dir_realis = this->saving_dir + "Entropy" + kPSep + "Eigenstate" + kPSep + "realisation=" + std::to_string(realis) + kPSep;
		createDirs(dir_realis);

		const arma::vec E = alfa.get_eigenvalues();
		const arma::cx_mat V = alfa.get_eigenvectors_full();
		arma::vec S(N, arma::fill::zeros);
	//#pragma omp parallel for num_threads(outer_threads) schedule(dynamic)
		for(int n = 0; n < N; n++){
			//arma::cx_vec state = alfa.get_state_in_full_Hilbert(n);
			arma::cx_vec state = V.col(n);
			double S_tmp = entropy::vonNeumann(state, LA, this->L, map);
			S(n) = S_tmp;
		}
		E.save(arma::hdf5_name(dir_realis + filename + ".hdf5", "energies"));
		S.save(arma::hdf5_name(dir_realis + filename + ".hdf5", "entropy", arma::hdf5_opts::append));
		//alfa.get_eigenvectors().save(arma::hdf5_name(filename + ".hdf5", "eigenvectors",arma::hdf5_opts::append));
		//save_to_file(filename + ".dat", E, entropies);
		energies += E;
		entropies += S;
		counter++;
	};

	//---- START COMPUTATION
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		if(alfa->using_Sz_symmetry()){
			std::cout << "Using Sz symmetry!!" << std::endl;
			map = generate_full_map(this->L, alfa->using_Sz_symmetry());
		}
		N = alfa->get_hilbert_size();
	 	energies = arma::vec(N, arma::fill::zeros);
		entropies = arma::vec(N, arma::fill::zeros);
		average_over_realisations<Ising_params::J>(*alfa, true, kernel);
	} else{
		auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
		N = alfa->get_hilbert_size();
	 	energies = arma::vec(N, arma::fill::zeros);
		entropies = arma::vec(N, arma::fill::zeros);
		if(alfa->using_Sz_symmetry())
			map = alfa->get_mapping();
		average_over_realisations<Ising_params::J>(*alfa, true, kernel);
	}
	energies /= double(counter);
	entropies /= double(counter);

	energies.save(arma::hdf5_name(dir + filename + "_jobid=" + std::to_string(this->jobid) + ".hdf5", "energies"));
	entropies.save(arma::hdf5_name(dir + filename + "_jobid=" + std::to_string(this->jobid) + ".hdf5", "entropy", arma::hdf5_opts::append));
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

	#ifdef HEISENBERG
		size_t N = binomial(this->L, this->L / 2.);
	#elif defined ANDERSON
		size_t N = this->L * this->L * this->L;
	#else
		size_t N = ULLPOW(this->L);
	#endif
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
			arma::vec energies_unfolded = statistics::unfolding(eigenvalues, std::min(20, this->L));
			//------------------- Get 50% spectrum
			double E_av = arma::trace(eigenvalues) / double(N);
			auto i = min_element(begin(eigenvalues), end(eigenvalues), [=](double x, double y) {
				return abs(x - E_av) < abs(y - E_av);
				});
			u64 E_av_idx = i - eigenvalues.begin();
			const long E_min = E_av_idx - num / 2.;
			const long E_max = E_av_idx + num / 2. + 1;
			const long num_small = (N > 1000)? 500 : 100;
			arma::vec energies = this->ch? exctract_vector(eigenvalues, E_av_idx - num_small / 2., E_av_idx + num_small / 2.) :
											exctract_vector(eigenvalues, E_min, E_max);
			arma::vec energies_unfolded_cut = this->ch? exctract_vector(energies_unfolded, E_av_idx - num_small / 2., E_av_idx + num_small / 2.) :
														 exctract_vector(energies_unfolded, E_min, E_max);
			
			//------------------- Level Spacings
			arma::vec level_spacings(energies.size() - 1, arma::fill::zeros);
			arma::vec level_spacings_unfolded(energies.size() - 1, arma::fill::zeros);
			for(int i = 0; i < energies.size() - 1; i++){
				const double delta 			= energies(i+1) 			 - energies(i);
				const double delta_unfolded = energies_unfolded_cut(i+1) - energies_unfolded_cut(i);

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
	std:string prefix = this->ch ? "_500_states" : "";
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
				ipr += statistics::inverse_participation_ratio(state) / double(N);
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
	std::string dir = this->saving_dir + "STATISTICS" + kPSep + "raw_data" + kPSep;
	std::string dir_agp = this->saving_dir + "AGP" + kPSep + opName + kPSep + "raw_data" + kPSep;
	createDirs(dir, dir_agp);
	std::ofstream file;

	openFile(file, dir + info + "_jobid=" + std::to_string(jobid) + ".dat");
	printSeparated(file, "\t", 25, true, "'gap ratio'", 			   				gap_ratio / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "'ipr'", 						 				  ipr / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "'information entropy'", 			 	 info_entropy / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "'information entropy random state'", 	 info_ent_rnd / double(counter));
	printSeparated(file, "\t", 25, true, "'entropy in ~100 states at E=0'", 	  	  entropy / double(num_ent * counter));
	printSeparated(file, "\t", 25, true, "'mean level spacing'", 			  			   wH / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "'typical level spacing'", 	  			   wH_typ / double(counter));
	file.close();

	openFile(file, dir_agp + info + "_jobid=" + std::to_string(jobid) + ".dat");
	printSeparated(file, "\t", 25, true, "'adiabatic gauge potential'", 	 			  AGP / double(counter));
	printSeparated(file, "\t", 25, true, "'typical fidelity susceptibility'", 	 	 typ_susc / double(counter));
	printSeparated(file, "\t", 25, true, "'fidelity susceptibility'", 			 	 	 susc / double(counter));
	printSeparated(file, "\t", 25, true, "'diagonal normalisation factor'", 		diag_norm / double(counter));

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
		size_t N = binomial(this->L, this->L / 2.);
	#elif defined ANDERSON
		size_t N = this->L * this->L * this->L;
	#else
		size_t N = ULLPOW(this->L);
	#endif
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		N = alfa->get_hilbert_size();
	}
	const double wH = sqrt(this->L) / (chi * N) * sqrt(this->J * this->J + this->h * this->h + this->g * this->g
												 + ( this->m? 0.0 : (this->w * this->w + this->g0 * this->g0 + this->J0 * this->J0) / 3. ));
	double tH = 1. / wH;
	double r1 = 0.0, r2 = 0.0;
	int num_times = 5000;
	int time_end = (int)std::ceil(std::log10(5 * tH));
	time_end = (time_end / std::log10(tH) < 1.5) ? time_end + 1 : time_end;

	arma::vec times = arma::logspace(log10(1.0 / (two_pi * N)), 1, num_times);
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
			u64 E_av_idx = spectrals::get_mean_energy_index(eigenvalues);
			const u64 num = this->L <= 9? 0.25 * N : 0.5 * N;
			const u64 num2 = this->L <= 12? 50 : 500;

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
	save_to_file(dir + info + ".dat", times, sff, 1.0 / wH_mean, thouless_time, r1, r2, N, 1.0 / wH_typ);
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
		size_t N = binomial(this->L, this->L / 2.);
	#elif defined(ANDERSON)
		size_t N = this->L * this->L * this->L;
	#else
		size_t N = ULLPOW(this->L);
	#endif
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		N = alfa->get_hilbert_size();
	}
	const double wH = sqrt(this->L) / (chi * N) * sqrt(this->J * this->J + this->h * this->h + this->g * this->g
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
	sff_fold = sff_fold / Z_folded;

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
	save_to_file(dir + info + ".dat", times, sff, tH, thouless_time, r1, r2, N);	
	smoothen_data(dir, info + ".dat");
	save_to_file(dir + "folded" + info + ".dat", times_fold, sff_fold, tH, thouless_time_fold, r1, r2, N);	
	smoothen_data(dir, "folded" + info + ".dat");
}

//<! find thouless time with various method as function of h,g,J
void isingUI::ui::thouless_times(Ising_params var)
{
	const int Lmin = this->L, Lmax = this->L + this->Ln * this->Ls;
	auto Jx_list = arma::linspace(this->J, this->J + this->Js * (this->Jn - 1), this->Jn);
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
		size_t N = 0;
		if(data.size() > 4){
			r1 = data[4](0);
			r2 = data[5](0);
		}
		if(data.size() > 6)
			N = data[6](0);
		// find thouless time
		double eps = 1e-1;
		auto K_GOE = [](double t){
			return t < 1? 2 * t - t * log(1+2*t) : 2 - t * log( (2*t+1) / (2*t-1) );
		};
		double thouless_time = 0;
		const double t_max = this->ch? 2.5 : 2.5 * tH;
		double delta_min = 1e6;
		for(int i = 0; i < sff.size(); i++){
			double t = times(i);
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
		printSeparated(map, "\t", 12, true, thouless_time, tH, r1, r2, N);
	};
	std::string dir = this->saving_dir + "ThoulessTime" + kPSep;
	createDirs(dir);
	std::string info = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"h", "g"}) 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, {"h", "g", "L", "J", "w"}, ",");
	std::ofstream map;
	openFile(map, dir + "_all" + info + ".dat", ios::out);
	for (int size = Lmin; size < Lmax; size += this->Ls){
		for(auto& Jx : Jx_list){
			for (auto &gx : gx_list){
				for (auto &hx : hx_list){
					for(auto& x : _list){	// either disorder w (m=0) or symmetry sector k (m=1)
						//printSeparated(std::cout, "\t", 16, true, size, Jx, gx, hx, x);
						kernel(size, Jx, gx, hx, x, map, size, Jx, gx, hx, x);
		}}}}}
	map.close();

	return;
	auto xarr = get_params_array(var); 
	std::string var_name = get_params_name(var);
	
	double eps = 1e-1;
	auto K_GOE = [](double t){
		return t < 1? 2 * t - t * log(1 + 2*t) : 2 - t * log( (2*t + 1) / (2*t - 1) );
	};
/*
	std::string dir = this->saving_dir + "ThoulessTime" + kPSep;	createDirs(dir);
	std::ofstream map;
	std::string sep = var_name == "L"? "," : "_";
	std::string info = this->generate_baseinfo({var_name}, sep);
	std::string baseinfo = this->generate_baseinfo();
	openFile(map, dir + info + ".dat", ios::out);
	std::cout << xarr << std::endl;
	for(auto& x : xarr)
	{
		// read SFF file
		std::ifstream file;
		std::string info_loc = this->update_info(baseinfo, var_name, x);
		std::string filename = this->saving_dir + "SpectralFormFactor/smoothed" + kPSep + info_loc + ".dat";
		auto data = readFromFile(file, filename);
		file.close();
		if (data.empty()) continue;

		arma::vec times = data[0];
		arma::vec sff = data[1];
		double tH = data[2](0);
		double r1 = 0.0, r2 = 0.0;
		size_t N = 0;
		if(data.size() > 4){
			r1 = data[4](0);
			r2 = data[5](0);
		}
		if(data.size() > 6)
			N = data[6](0);

		// find thouless time
		double thouless_time = 0;
		const double t_max = this->ch? 2.5 : 2.5 * tH;
		double delta_min = 1e6;
		for(int i = 0; i < sff.size(); i++){
			double t = times(i);
			double delta = abs( log10( sff(i) / K_GOE(t) ));

			delta = delta - eps;
			delta *= delta;
			if(delta < delta_min){
				delta_min = delta;
				thouless_time = times(i); 
			}
			if(times(i) >= t_max) break;
		}
		printSeparated(std::cout, "\t", 12, true, x, thouless_time);
		printSeparated(map, "\t", 12, true, x, thouless_time, tH, r1, r2, N);
	}
	map.close();
	*/
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
		size_t N = binomial(Lx, Lx / 2.);
	#else
		size_t N = ULLPOW(Lx);
	#endif
		if(this->m){
			auto alfa = std::make_unique<IsingModel_sym>(Lx, this->J, gx, hx,
									 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
			N = alfa->get_hilbert_size();
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
			long n_bins = 1 + long(3.322 * log(N - 2));
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
	   	 errS_dis = 0.0,//<! half-chain entropy error with disorder
	   	 var_S = 0.0,	//<! half-chain entropy variance
	  		 ipr = 0.0,	//<! inverse participation ratio
	info_entropy = 0.0,	//<! information entropy in eigenstates
	   gap_ratio = 0.0,	//<! gap ratio
	   		  wH = 0.0,	//<! mean level spacing
	   	  wH_typ = 0.0;	//<! typical level spacing
	std::string info;
	std::vector<u64> map;
	#ifdef HEISENBERG
		const long num_ent = L >= 14? 1000 : (L >= 12? 400 : 100);
	#else
		const long num_ent = L >= 12? 1000 : (L >= 10? 400 : 100);
	#endif
	//---- KERNEL LAMBDA
	auto kernel = [&](auto& alfa, int realis)
	{
		const u64 N = alfa.get_hilbert_size();
		this->mu = long(0.5 * N);
	
		info = alfa.get_info();
		std::cout << info << std::endl;
		long int E_min = alfa.E_av_idx - long(mu / 2);
		long int E_max = alfa.E_av_idx + long(mu / 2);
		const arma::vec E = alfa.get_eigenvalues();

		double wH_typ_local = 0.0;
		double S = 0, S2 = 0;
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
	
				arma::cx_vec state = alfa.get_state_in_full_Hilbert(i);
				ipr += statistics::inverse_participation_ratio(state) / double(N);
				info_entropy += statistics::information_entropy(state);
				if(i >= alfa.E_av_idx - num_ent / 2. && i < alfa.E_av_idx + num_ent / 2.)
				{
					double Stmp = entropy::vonNeumann(state, this->L / 2, this->L, map);
					S += Stmp;
					S2 += Stmp * Stmp;
				}
			}
		}
		S /= double(num_ent);
		S2 /= double(num_ent);
		entropy += S;
		var_S += S2 - S * S;
		errS_dis += S * S;
		wH_typ += std::exp(wH_typ_local / double(this->mu));

		counter++;
		std::cout << counter << std::endl;
	};

	//---- START COMPUTATION
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		if(alfa->using_Sz_symmetry())
			map = generate_full_map(this->L, alfa->using_Sz_symmetry());
			
		average_over_realisations<Ising_params::h>(*alfa, true, kernel);
	} else{
		auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
		if(alfa->using_Sz_symmetry())
			map = alfa->get_mapping();
		average_over_realisations<Ising_params::h>(*alfa, true, kernel);
	}
	errS_dis = errS_dis / double(counter) - entropy / double(counter) * entropy / double(counter);

	std::string dir = this->saving_dir + "STATISTICS" + kPSep + "raw_data" + kPSep;
	createDirs(dir);
	std::ofstream file;
	openFile(file, dir + info + "_jobid=" + std::to_string(jobid) + ".dat");

	printSeparated(file, "\t", 25, true, "'gap ratio'", 			   				gap_ratio / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "'ipr'", 						 				  ipr / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "'information entropy'", 			 	 info_entropy / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "'entropy in ~100 states at E=0'", 	  	  entropy / double(counter));
	printSeparated(file, "\t", 25, true, "'mean level spacing'", 			  			   wH / double(this->mu * counter));
	printSeparated(file, "\t", 25, true, "'typical level spacing'", 	  			   wH_typ / double(counter));
	printSeparated(file, "\t", 25, true, "'entropy var in ~100 states at E=0'", 	    var_S / double(counter));
	printSeparated(file, "\t", 25, true, "'entropy error over realisations'", 	  	 errS_dis);

	file.close();
}

//<! generate map from statistics data
void isingUI::ui::generate_statistic_map(Ising_params varname){
	auto xarr = get_params_array(varname); 
	std::string str = get_params_name(varname);
	std::string sep = str == "L"? "," : "_";
	std::string info = this->generate_baseinfo({str}, sep);
	std::string baseinfo = this->generate_baseinfo();

	std::ofstream map;
	std::ifstream datafile;
	std::string dir = "";
	if(this->ch){ 
		auto [opName, subdir] = IsingModel_disorder::opName(this->op, this->site);
		dir = this->saving_dir + "AGP" + kPSep + opName + kPSep;
		openFile(map, dir + info + ".dat");
		printSeparated(map, "\t", 25, true, "'" + str + "'", "'adiabatic gauge potential'", 	 			
									"'typical fidelity susceptibility'", "'fidelity susceptibility'", "'diagonal normalisation factor'");
	}
	else { 
		dir = this->saving_dir + "STATISTICS" + kPSep;
		openFile(map, dir + info + ".dat");
		printSeparated(map, "\t", 25, true, "'" + str + "'", "'gap ratio'", "'ipr'", "'information entropy'",
											 "'entropy in ~100 states at E=0'", "'mean level spacing'", "'typical level spacing'",
											 "'entropy var in ~100 states at E=0'", "'entropy error over realisations'");
	}
	for(auto& x : xarr){
		std::string info_loc = update_info(baseinfo, str, x);
		bool status = openFile(datafile, dir + "raw_data" + kPSep + info_loc + ".dat");
		if(!status) continue;
		std::string line;
		double value;
		printSeparated(map, "\t", 25, false, x);
		while(std::getline(datafile, line)){
			std::istringstream ss(line);
			std::vector<std::string> datarow;
			while(ss){
				std::string element;
				ss >> element;
				datarow.push_back(element);
			}
			double value = std::stod(datarow[datarow.size()-2]);
			printSeparated(map, "\t", 25, false, value);
		}
		map << std::endl;
		datafile.close();
	}
	map.close();
}
void isingUI::ui::generate_statistic_map(Ising_params varname1, Ising_params varname2){
	auto xarr = get_params_array(varname1);
	auto yarr = get_params_array(varname2);
	std::string str1 = get_params_name(varname1), str2 = get_params_name(varname2);

	std::string info = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {str1, str2}) 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, {str1, str2});
	std::string baseinfo = this->m? IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym) 
					: IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w);

	std::ofstream map;
	std::ifstream datafile;
	std::string dir = "";
	if(this->ch){ 
		auto [opName, subdir] = IsingModel_disorder::opName(this->op, this->site);
		dir = this->saving_dir + "AGP" + kPSep + opName + kPSep;
		openFile(map, dir + info + ".dat");
		printSeparated(map, "\t", 25, true, "'" + str1 + "'", "'" + str2 + "'", "'adiabatic gauge potential'", 	 			
									"'typical fidelity susceptibility'", "'fidelity susceptibility'", "'diagonal normalisation factor'");
	}
	else { 
		dir = this->saving_dir + "STATISTICS" + kPSep;
		openFile(map, dir + info + ".dat");
		printSeparated(map, "\t", 25, true, "'" + str1 + "'", "'" + str2 + "'", "'gap ratio'", "'ipr'", "'information entropy'",
											 "'entropy in ~100 states at E=0'", "'mean level spacing'", "'typical level spacing'");
	}
	for(auto& x : xarr){
		std::string info_loc = update_info(baseinfo, str1, x);
		for(auto& y : yarr){
			info_loc = update_info(info_loc, str2, y);
			bool status = openFile(datafile, dir + "raw_data" + kPSep + info_loc + ".dat");
			if(!status) continue;
			std::string line;
			double value;
			printSeparated(map, "\t", 25, false, x, y);
			while(std::getline(datafile, line)){
				std::istringstream ss(line);
				std::vector<std::string> datarow;
				while(ss){
					std::string element;
					ss >> element;
					datarow.push_back(element);
				}
				double value = std::stod(datarow[datarow.size()-2]);
				printSeparated(map, "\t", 25, false, value);
			}
			map << std::endl;
			datafile.close();
		}
	}
	map.close();
}

//-------------------------------------------------------------------------- ANDERSON


void isingUI::ui::calculate_localisation_length(){
	arma::vec energies(this->L, arma::fill::zeros);
	arma::mat corr_func(this->L / 2, this->L, arma::fill::zeros);

	auto lambda_average = [&](
		int realis, double x
		)
	{
		arma::vec E;
		arma::mat corr_functmp;
		std::tie(E, corr_functmp) = anderson::get_localisation_length1D(this->L, this->J, this->w);
		#pragma omp critical
		{
			energies +=E;
			corr_func += corr_functmp;
		}
	};

	average_over_realisations<Ising_params::g>(false, lambda_average);
	energies /= double(this->realisations);
	corr_func /= double(this->realisations);
	arma::vec loc_length(this->L, arma::fill::zeros);
	
	for(int i = 0; i < this->L; i++){
		arma::vec corr = corr_func.col(i);
        double _min = arma::min(corr);
        if(!std::isfinite(_min) || _min < -20.0)
            _min = -20.0;

		arma::vec r_vals  = arma::linspace(0, this->L / 2., corr.size());
        
		if(this->L < 20 || (i % (this->L / 20) == 0))
			save_to_file("./results/ANDERSON/1D/PBC/CorrelationFunction/_L=" 
                            + std::to_string(this->L) + "_n=" + std::to_string(i) + "_w=" + to_string_prec(this->w, 2) + ".dat", r_vals, corr);
        arma::vec func_to_fit;
		func_to_fit = exctract_vector_between_values(corr, _min, 1.0);
        r_vals  = arma::linspace(0, func_to_fit.size(), func_to_fit.size());
        arma::vec p = arma::polyfit(r_vals, func_to_fit, 1);
        loc_length(i) = -1. / p(0);
	}
	std::string dir = this->saving_dir + "LocalisationLength" + kPSep + "Distribution" + kPSep;
	createDirs(dir);
	save_to_file(dir + IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w) + ".dat", energies / (this->J + this->w), loc_length);
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
		"	8 -- J    - spin current\n"
		"	9 -- Sx    - global X-magnetization\n"
		"	10 -- Sz   - global Z-magnetization\n"
		"	11 -- sum_i Sx_i Sx_i+1   - nearest neighbour X\n"
		"	12 -- sum_i Sz_i Sz_i+1   - nearest neighbour Z\n"
		"	13 -- sum_i Sx_i Sx_i+2   - next nearest neighbour X\n"
		"	14 -- sum_i Sz_i Sz_i+2   - next nearest neighbour Z\n"
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
		"\t	4 -- calculate localisation length for all eigenstates"
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

	#if MODEL == 0
		std::cout << "ISING:\n\t\t" << "H = \u03A3_i J_i \u03C3^z_i \u03C3^z_i+1 + g_i \u03C3^x_i + h_i \u03C3^x_i" << std::endl << std::endl;
	#elif MODEL == 1
		std::cout << "HEISENBERG:\n\t\t" << "H = \u03A3_i J_i(S^x_i S^x_i+1 + S^y_i S^y_i+1) + g_i S^z_iS^z_i+1 + h_i S^x_i" << std::endl << std::endl;
	#elif MODEL == 2
		std::cout << "ANDERSON:\n\t\t" << "H = J/2 \u03A3_i (S^x_i S^x_i+1 + S^y_i S^y_i+1) + h_i S^z_i" << std::endl << std::endl;
		std::cout << "h_i \u03B5 [-w, w]" << std::endl;
		std::cout << "h_c^{3D} ~ 4.1 (rescale by 4 when using 3D in plots)" << std::endl;
	#elif MODEL == 3
		std::cout << "XYZ:\n\t\t" << "H = \u03A3_r J_r\u03A3_i [ (1 - J0)S^x_i S^x_i+r + (1 + J0)S^y_i S^y_i+r) + g S^z_iS^z_i+1 + g0 S^x_i + h S^z_i + w(S^z_{L-1} + S^x_0)" << std::endl << std::endl;
	#elif MODEL == 4
		std::cout << "QUANTUM SUN:\n\t\t" << "H = R_GOE + \u03A3_i g0 * alfa^{u_j} S^x_i S^x_i+1 + h_i S^x_i" << std::endl << std::endl;
		std::cout << "u_j in [j - J0, j + J0]"  << std::endl;
		std::cout << "h_i \u03B5 [h - w, h + w]" << std::endl;
	#else
		std::cout << "ISING:\n\t\t" << "H = \u03A3_i J_i \u03C3^z_i \u03C3^z_i+1 + g_i \u03C3^x_i + h_i \u03C3^x_i" << std::endl << std::endl;
	#endif
	#if MODEL == 0 || MODEL == 1
		std::cout << "J_i \u03B5 [J - J0, J + J0]" << std::endl;
		std::cout << "g_i \u03B5 [g - g0, g + g0]" << std::endl;
		std::cout << "h_i \u03B5 [h - w, h + w]" << std::endl;
	#endif
	stout << "Chosen spin species is S = " << S << std::endl << std::endl;
	stout << "------------------------------CHOSEN OPTIONS:" << std::endl;
	std::string opName = std::get<0>(IsingModel_disorder::opName(this->op, this->site));
	stout << "DIR = " << this->saving_dir << std::endl
		  << "model = " << (this->m ? "symmetric" : "disordered") << std::endl
		  << "BC = " << (this->boundary_conditions ? "OBC" : "PBC") << std::endl
		  << "L  = " << this->L << std::endl
		  << "Ls = " << this->Ls << std::endl
		  << "Ln = " << this->Ln << std::endl
		  << "J  = " << this->J << std::endl
		  << "Jn = " << this->Jn << std::endl
		  << "Js = " << this->Js << std::endl
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
		  << "scale = " << (this->scale == 1 ? "log" : "linear") << std::endl
		  << "realisations = " << this->realisations << std::endl
		  << "seed = " << this->seed << std::endl;

	if (this->m == 0)
		stout << "J0  = " << this->J0 << std::endl
			  << "J0n = " << this->J0n << std::endl
			  << "J0s = " << this->J0s << std::endl
			  << "w   = " << this->w << std::endl
			  << "ws  = " << this->ws << std::endl
			  << "wn  = " << this->wn << std::endl
			  << "g0  = " << this->g0 << std::endl
			  << "g0s = " << this->g0s << std::endl
			  << "g0n = " << this->g0n << std::endl;

	if (this->m == 1){
		if(this->symmetries.k_sym == -1) stout << "k-sector = pi" << std::endl;
		else stout << "k-sector = " << 2 * this->symmetries.k_sym << "/L *pi" << std::endl;
		stout << "p-sector = " << (this->symmetries.p_sym ? 1 : -1) << std::endl
			  << "x-sector = " << (this->symmetries.x_sym ? 1 : -1) << std::endl;
	}
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
	this->Js = 0.0;
	this->Jn = 1;
	this->J0 = 0.0;
	this->J0s = 0.0;
	this->J0n = 1;

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
	this->dim = 3;
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
	choosen_option = "-Js";
	this->set_option(this->Js, argv, choosen_option, false);
	choosen_option = "-Jn";
	this->set_option(this->Jn, argv, choosen_option);

	// spin coupling disorder
	choosen_option = "-J0";
	this->set_option(this->J0, argv, choosen_option);
	choosen_option = "-J0s";
	this->set_option(this->J0s, argv, choosen_option, false);
	choosen_option = "-J0n";
	this->set_option(this->J0n, argv, choosen_option);

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

	//choose dimensionality
	choosen_option = "-dim";
	this->set_option(this->dim, argv, choosen_option);
	anderson_dim = this->dim;

	// model
	choosen_option = "-m";
	this->set_option(this->m, argv, choosen_option);
	if (this->m > 1)
		this->set_default_msg(this->m, choosen_option.substr(1),
							  "max model number is 1", table);

	// choose function
	choosen_option = "-fun";
	this->set_option(this->fun, argv, choosen_option, false);

	// buckets
	choosen_option = "-mu";
	this->set_option(this->mu, argv, choosen_option);
	choosen_option = "-r";
	this->set_option(this->realisations, argv, choosen_option);

	// symmetries
	choosen_option = "-k";
	this->set_option(this->symmetries.k_sym, argv, choosen_option, false);
	if (this->symmetries.k_sym >= this->L)
		this->set_default_msg(this->symmetries.k_sym, choosen_option.substr(1),
							  "max k sector is L = " + std::to_string(this->L), table);
	if(this->symmetries.k_sym < -1)
		this->set_default_msg(this->symmetries.k_sym, choosen_option.substr(1),
							  "negative k sector is only -1 for k=pi for all", table);

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
	#ifdef LOCAL_PERT
		str_model = "local_pert" + std::string(kPathSeparator);
	#else
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
	#endif

	#ifdef ANDERSON
		str_model = std::to_string(this->dim) + "D" + kPSep;
	#endif

	// make boundary condition folder
	#if !defined(QUANTUM_SUN)
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
	#endif

	#if MODEL == 0
		std::string folder = saving_dir + "ISING" + kPSep + str_model;
	#elif MODEL == 1
		std::string folder = saving_dir + "HEISENBERG" + kPSep + str_model;
	#elif MODEL == 2
		std::string folder = saving_dir + "ANDERSON" + kPSep + str_model;
	#elif MODEL == 3
		std::string folder = saving_dir + "XYZ" + kPSep + str_model;
	#elif MODEL == 4
		std::string folder = saving_dir + "QUANTUM_SUN" + kPSep + str_model;
	#else
		std::string folder = saving_dir + "ISING" + kPSep + str_model;
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
	energies.load(arma::hdf5_name(base + "eigvals_L_16_N_3_pbc_True_disorder_uniform_ham_type_anderson_t_-1.0_W_0.0_dW_" + to_string_prec(W, 2) + ".hdf5", "Eigenvalues"));
	std::cout << energies.n_cols << "\t" << energies.n_rows << std::endl;
	arma::mat sff_jan1, sff_jan2, sff_jan3;
	bool state1 = sff_jan1.load(arma::hdf5_name(base + "eigvals_L_16_N_3_pbc_True_disorder_uniform_ham_type_anderson_t_-1.0_W_0.0_dW_" + to_string_prec(W, 2) + ".hdf5", "SFF_spectra_eta_0.1000_filter_gaussian"));
	bool state2 = sff_jan2.load(arma::hdf5_name(base + "eigvals_L_16_N_3_pbc_True_disorder_uniform_ham_type_anderson_t_-1.0_W_0.0_dW_" + to_string_prec(W, 2) + ".hdf5", "SFF_spectra_eta_0.3000_filter_gaussian"));
	bool state3 = sff_jan3.load(arma::hdf5_name(base + "eigvals_L_16_N_3_pbc_True_disorder_uniform_ham_type_anderson_t_-1.0_W_0.0_dW_" + to_string_prec(W, 2) + ".hdf5", "SFF_spectra_eta_0.5000_filter_gaussian"));
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