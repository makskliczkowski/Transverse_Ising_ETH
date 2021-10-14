#include "include/user_interface.h"

/// <summary>
/// We want to handle files so let's make the c-way input a string
/// </summary>
/// <param name="argc"> number of main input arguments </param>
/// <param name="argv"> main input arguments </param>
/// <returns></returns>
std::vector<std::string> change_input_to_vec_of_str(int argc, char** argv)
{
	// -1 because first is the name of the file
	std::vector<std::string> tmp(argc - 1, "");
	for (int i = 0; i < argc - 1; i++) {
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
string user_interface::getCmdOption(const v_1d<string>& vec, string option) const
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
template<typename T>
void user_interface::set_option(T& value, const v_1d<string>& argv, string choosen_option, bool geq_0)
{
	if (std::string option = this->getCmdOption(argv, choosen_option); option != "")
		value = static_cast<T>(stod(option));												// set value to an option
	if (geq_0 && value < 0)																	// if the variable shall be bigger equal 0
		this->set_default_msg(value, choosen_option.substr(1), \
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
template<typename T>
void user_interface::set_default_msg(T& value, string option, string message, const unordered_map<string, string>& map) const
{
	stout << message;																	// print warning
	std::string value_str = "";															// we will set this to value
	if (auto it = map.find(option); it != map.end()) {
		value_str = it->second;															// if in table - we take the enum
	}
	value = stod(value_str);
}

// ----------------------------------------------------------------------------- ISING MODEL -----------------------------------------------------------------------------

// ----------------------------------------------------------------------------- Connected with the parser
/// <summary>
/// Setting the default parameters for the Ising model
/// </summary>
void isingUI::ui::set_default() {
	using namespace isingUI;
	this->saving_dir = "." + std::string(kPathSeparator) + "results" + std::string(kPathSeparator);		// directory for the result files to be saved into
	this->L = 4;
	this->Ls = 0;
	this->Ln = 1;

	this->J = 1.0;
	this->J0 = 0.2;

	this->h = 0.0;
	this->hs = 0.0;
	this->hn = 1;

	this->w = 1.0;
	this->ws = 0.0;
	this->wn = 1;

	this->g = 1.0;
	this->gs = 0.0;
	this->gn = 1;

	this->g0 = 0;

	this->symmetries.k_sym = 0;
	this->symmetries.p_sym = 0;
	this->symmetries.x_sym = 0;

	this->realisations = 100;
	this->site = 0;
	this->mu = 5;

	this->boundary_conditions = 0;
	this->m = 0;
	this->p = true;
	this->thread_number = 1;
}

// ------------------------------------- CONSTURCTORS

/// <summary>
/// Ising model user interface constructor, we pass there the command lines arguments from main
/// </summary>
/// <param name="argc"> number of arguments </param>
/// <param name="argv"> arguments list </param>
isingUI::ui::ui(int argc, char** argv)
{
	auto input = change_input_to_vec_of_str(argc, argv);									// change standard input to vec of strings
	input = std::vector<std::string>(input.begin()++, input.end());							// skip the first element which is the name of file
	// plog::init(plog::info, "log.txt");														// initialize logger
	if (std::string option = this->getCmdOption(input, "-f"); option != "") {
		input = this->parseInputFile(option);												// parse input from file
	}
	this->parseModel(input.size(), input);													// parse input from CMD directly
}

/// <summary>
/// Function that tells how does the parser work
/// </summary>
void isingUI::ui::exit_with_help() const {
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
		"	2 -- ABC\n"
		"-m model to be choosen : (default 0 - without symmetries)\n"
		"	0 -- nonsymmetric model - only here the disorder is working\n"
		"	1 -- include symmetries - here the parity flag is also working\n"
		"-k translation symetry sector, 0-L, (default 0)\n"
		"-p parity symmetry sector, +-1 (if applicable) (default 1)\n"
		"-x spin flip symmetry sector, +-1 (if applicable) (default 1)\n"
		"-th number of threads to be used for CPU parallelization : depends on the machine specifics, default(1)"
		"-h quit with help\n"
	);
	std::exit(1);
}

// ------------------------------------- PARSERS

/// <summary>
/// The parser for the Transverse Field Ising model
/// </summary>
/// <param name="argc"> number of arguments </param>
/// <param name="argv"> list of arguments </param>
void isingUI::ui::parseModel(int argc, std::vector<std::string> argv) {
	using namespace isingUI;
	// SET DEFAULT VALUES
	this->set_default();																				// setting default at the very beginning

	std::string choosen_option = "";																	// current choosen option
	std::string str_model = std::string(kPathSeparator) + "disorder" + std::string(kPathSeparator);		// folder for current model
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
	this->set_option(this->ws, argv, choosen_option);
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
	if (this->boundary_conditions > 2) this->set_default_msg(this->boundary_conditions, choosen_option.substr(1), \
		"max boundary condition is 2", table);

	// model
	choosen_option = "-m";
	this->set_option(this->m, argv, choosen_option);
	if (this->m > 1) this->set_default_msg(this->m, choosen_option.substr(1), \
		"max model number is 1", table);

	// buckets
	choosen_option = "-mu";
	this->set_option(this->mu, argv, choosen_option);

	// symmetries
	choosen_option = "-k";
	this->set_option(this->symmetries.k_sym, argv, choosen_option);
	if (this->symmetries.k_sym >= this->L) this->set_default_msg(this->symmetries.k_sym, choosen_option.substr(1), \
		"max k sector is L = " + std::to_string(this->L), table);
	choosen_option = "-p";
	this->set_option(this->symmetries.p_sym, argv, choosen_option, false);
	choosen_option = "-x";
	this->set_option(this->symmetries.x_sym, argv, choosen_option, false);

	// thread number
	choosen_option = "-th";
	this->set_option(this->thread_number, argv, choosen_option);
	if (this->thread_number > std::thread::hardware_concurrency())
		this->set_default_msg(this->thread_number, choosen_option.substr(1), \
			"Wrong number of threads\n", table);
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
	default:
		str_model += "PBC" + std::string(kPathSeparator);
		break;
	}

	std::string folder = saving_dir + str_model;
	if (!argv[argc - 1].empty() && argc % 2 != 0) {
		// only if the last command is non-even
		folder = argv[argc - 1] + str_model;
		if (fs::create_directories(folder) || fs::is_directory(folder))						// creating the directory for saving the files with results
			this->saving_dir = folder;																// if can create dir this is is
	}
	else {
		if (fs::create_directories(folder) || fs::is_directory(folder))						// creating the directory for saving the files with results
			this->saving_dir = folder;																// if can create dir this is is
	}

	std::cout << " - - - - - - MAKING ISING INTERFACE AND USING OUTER THREADS : " \
		<< thread_number << " - - - - - - " << endl;										// setting the number of threads to be used with omp

	omp_set_num_threads(this->thread_number);
	return;
}

/// <summary>
/// If the commands are given from file, we must treat them the same as arguments
/// </summary>
/// <param name="filename"> the name of the file that contains the command line </param>
/// <returns></returns>
std::vector<std::string> user_interface::parseInputFile(std::string filename) const {
	std::vector<std::string> commands(1, "");
	ifstream inputFile(filename);
	std::string line = "";
	if (!inputFile.is_open())
		stout << "Cannot open a file " + filename + " that I could parse. All parameters are default. Sorry :c \n";
	else {
		if (std::getline(inputFile, line)) {
			commands = split_str(line, " ");									// saving lines to out vector if it can be done, then the parser shall treat them normally
		}
	}
	return std::vector<std::string>(commands.begin(), commands.end());
}

// ----------------------------------------------------------------------------- SIMULATIONS -----------------------------------------------------------------------------

/// <summary>
/// Creates two models, one with symmetries implemented for each symmetry sector and one without and compares their corresponding energies if they agree
/// </summary>
void isingUI::ui::compare_energies() {
	// handle disorder
	auto Hamil = std::make_unique<IsingModel_disorder>(L, J, 0, g, 0, h, 0, boundary_conditions);
	Hamil->diagonalization();
	const vec E_dis = Hamil->get_eigenvalues();
	Hamil.release();
	// here we push back the energies from symmetries
	std::vector<double> E_sym = v_1d<double>();
	std::vector<std::string> symms = v_1d<std::string>();
	// go for each symmetry sector
	for (int k = 0; k < L; k++) {
		if (k == 0 || k == this->L / 2.) {
			for (int p = 0; p <= 1; p++) {
				// if the spin flip is unaviable we just use 1
				const int x_max = (this->h != 0) ? 0 : 1;
				for (int x = 0; x <= x_max; x++) {
					auto ham = std::make_unique<IsingModel_sym>(L, J, g, h, k, p, x, boundary_conditions);
					ham->diagonalization();
					vec t = ham->get_eigenvalues();
					E_sym.insert(E_sym.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));
					v_1d<std::string> temp_str = v_1d<std::string>(t.size(), "k=" + std::to_string(k) + ",x=" + to_string(x) + ",p=" + to_string(p));
					symms.insert(symms.end(), std::make_move_iterator(temp_str.begin()), std::make_move_iterator(temp_str.end()));
				}
			}
		}
		else {
			int x_max = (this->h != 0) ? 0 : 1;
			for (int x = 0; x <= x_max; x++) {
				auto ham = std::make_unique<IsingModel_sym>(L, J, g, h, k, 1, x, boundary_conditions);
				ham->diagonalization();
				vec t = ham->get_eigenvalues();
				E_sym.insert(E_sym.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));
				v_1d<std::string> temp_str = v_1d<std::string>(t.size(), "k=" + std::to_string(k) + ",x=" + to_string(x));
				symms.insert(symms.end(), std::make_move_iterator(temp_str.begin()), std::make_move_iterator(temp_str.end()));
			}
		}
		stout << k << endl;
	}
	// face permutation
	auto permut = sort_permutation(E_sym, [](const double a, const double b) {return a < b; });
	apply_permutation(E_sym, permut);
	apply_permutation(symms, permut);
	stout << endl << E_sym.size() << endl;
	for (int k = 0; k < E_dis.size(); k++) {
		stout << symms[k] << "\t\t\t\t" << E_sym[k] << "\t\t\t\t" << E_dis(k) << "\t\t\t\t" << E_sym[k] - E_dis(k) << endl;
	}
}

/// <summary>
/// Handle disorder simulation
/// </summary>
void isingUI::ui::disorder() {
	auto start = std::chrono::high_resolution_clock::now();
	/*std::ofstream scaling_r_sigmaX(this->saving_dir + "SpectrumRapScalingSigmaX" + \
		",J0=" + to_string_prec(this->J0, 2) + \
		",g=" + to_string_prec(this->g, 2) + \
		",g0=" + to_string_prec(this->g0, 2) + \
		",h=" + to_string_prec(this->h, 2) + \
		",w=" + to_string_prec(this->w, 2) + ".dat", std::ofstream::app);
	std::ofstream scaling_ipr(this->saving_dir + "iprScaling" + \
		",J0=" + to_string_prec(this->J0, 2) + \
		",g=" + to_string_prec(this->g, 2) + \
		",g0=" + to_string_prec(this->g0, 2) + \
		",h=" + to_string_prec(this->h, 2) + \
		",w=" + to_string_prec(this->w, 2) + ".dat", std::ofstream::app);
	for (L = 10; L <= 14; L += 1) {
		realisations = 1000 - L * 50;

		auto Hamil = std::make_unique<IsingModel_disorder>(L, J, J0, g, g0, h, w);
		const u64 N = Hamil->get_hilbert_size();
		Hamil->diagonalization();

		vec av_sigma_x = Hamil->operator_av_in_eigenstates_return(&IsingModel_disorder::av_sigma_x, *Hamil, 0);
		std::ofstream average_sigma_x(this->saving_dir + "SigmaX" + Hamil->get_info() + ".dat");
		for (int k = 0; k < N; k++)
			average_sigma_x << Hamil->get_eigenEnergy(k) / double(L) << "\t\t\t\t" << av_sigma_x(k) << endl;
		average_sigma_x.close();
		stout << "--> finished writing the sigma _x  average for : " << Hamil->get_info() << " <--\n";

		vec r_sigma_x(N - 1);
#pragma omp parallel for
		for (int k = 0; k < N - 1; k++)
			r_sigma_x(k) = av_sigma_x(k + 1) - av_sigma_x(k);
		// spectrum repulsion for < sigma_0^x >

		vec stat_aver = statistics_average(r_sigma_x, 10);
		// stoutliers and scaling

		probability_distribution(this->saving_dir, "rSigmaXDist" + Hamil->get_info(), r_sigma_x, -0.5, 0.5, 0.01);
		stout << "--> finished writing the probability distribution for r_sigma _x repuslion and outliers for : " << Hamil->get_info() << " <--\n";

		double ipr = 0, entropy = 0, r = 0;
		int mu = 100;
		for (int k = 0; k < realisations; k++) {
			Hamil->hamiltonian();
			Hamil->diagonalization();
			vec av_sigma_x = Hamil->operator_av_in_eigenstates_return(&IsingModel::av_sigma_x, *Hamil, 0);
			vec r_sigma_x(N - 1);
#pragma omp parallel for
			for (int k = 0; k < N - 1; k++)
				r_sigma_x(k) = av_sigma_x(k + 1) - av_sigma_x(k);
			stat_aver += statistics_average(r_sigma_x, 10);
			// average in middle spectrum
			for (int f = Hamil->E_av_idx - mu / 2.; f < Hamil->E_av_idx + mu / 2.; f++) {
				ipr += Hamil->ipr(f);
				entropy += Hamil->information_entropy(f);
				r += Hamil->eigenlevel_statistics(f, f + 1);
			}
			if (k % 5 == 0) stout << " \t\t--> " << k << " - in time : " << \
				double(duration_cast<milliseconds>(duration(high_resolution_clock::now() - start)).count()) / 1000.0 << "s" << std::endl;
		}
		stout << "--> finished averaging over realizations for : " << Hamil->get_info() << " <--\n\n\t\t\n\b";
		scaling_r_sigmaX << N << stat_aver.t() / double(realisations);
		double norm = realisations * mu;
		scaling_ipr << L << "\t\t\t\t" << N << "\t\t\t\t" << ipr / norm / (double)N << "\t\t\t\t" << entropy / norm << "\t\t\t\t" << r / norm << endl;
	}
	scaling_r_sigmaX.close();
	scaling_ipr.close();
	*/
	this->m = 0;
	std::string info = "_L=" + std::to_string(this->L) + \
		",J0=" + to_string_prec(this->J0, 2) + \
		",g=" + to_string_prec(this->g, 2) + \
		",g0=" + to_string_prec(this->g0, 2) + \
		",h=" + to_string_prec(this->h, 2);
	std::ofstream scaling_ipr(this->saving_dir + "iprDisorder" + info + ".dat");
	std::ofstream kurt(this->saving_dir + "Moments" + info + ".dat");
	kurt << "w\t\tSigmaX_kurtosis\tSigmaX_binder" << endl;
	realisations = 300;
	for (double w = 0.0; w <= 5.0; w += 0.2) {
		auto Hamil = std::make_unique<IsingModel_disorder>(L, J, J0, g, g0, h, w);
		gen = std::mt19937_64(seed);
		stout << "\n\n------------------------------ Doing : " << Hamil->get_info() << "------------------------------\n";
		const u64 N = Hamil->get_hilbert_size();
		Hamil->diagonalization();
		stout << " \t\t--> finished diagonalizing for " << Hamil->get_info() << " - in time : " << tim_s(start) << "s\n";
		this->mu = 0.5 * N;

		const long int E_min = Hamil->E_av_idx - mu / 2.;
		const long int E_max = Hamil->E_av_idx + mu / 2.;

		vec av_sigma_x(mu, arma::fill::zeros);
		for (u64 k = E_min; k < E_max; k++) {
			const int idx = k - E_min;
			av_sigma_x(idx) = Hamil->av_sigma_x(k, k, { 0 });
		}
		vec fluct = data_fluctuations(av_sigma_x);
		double _min = -2.0, _max = 2.0, step = 2e-3;
		stout << "--> finished writing the sigma _x fluctuations for w = " << w << " <--\n";

		arma::vec prob_dist = probability_distribution_with_return(fluct);
		arma::vec prob_dist_GOE = probability_distribution_with_return(Hamil->eigenlevel_statistics_with_return());
		double ipr = 0, entropy = 0;
		double r = 0;
		for (int k = 0; k < realisations - 1; k++) {
			Hamil->hamiltonian();
			Hamil->diagonalization();
			for (u64 k = E_min; k < E_max; k++) {
				const int idx = k - E_min;
				av_sigma_x(idx) = Hamil->av_sigma_x(k, k, { 0 });
			}
			fluct = data_fluctuations(av_sigma_x);

			prob_dist += probability_distribution_with_return(fluct);
			prob_dist_GOE += probability_distribution_with_return(Hamil->eigenlevel_statistics_with_return());

			// average in middle spectrum
			for (long int f = E_min; f < E_max; f++) {
				ipr += Hamil->ipr(f);
				entropy += Hamil->information_entropy(f);
				r += Hamil->eigenlevel_statistics(f, f + 1);
			}
			if (k % 5 == 0) stout << " \t\t--> " << k << " realisation - in time : " << tim_s(start) << "s" << std::endl;
		}
		stout << "--> finished loop over realisations for w = " << w << " <--\n";

		kurt << w << "\t\t" << binder_cumulant(prob_dist) << "\t\t" << kurtosis_diff(prob_dist) << "\t\t" << endl;
		save_to_file(this->saving_dir, "ProbDistSigmaX" + Hamil->get_info(), arma::linspace(_min, _max, prob_dist.size()), prob_dist);

		save_to_file(this->saving_dir, "ProbDistGap" + Hamil->get_info(), arma::linspace(0, 1, prob_dist_GOE.size()), prob_dist_GOE);

		double norm = realisations * mu;
		scaling_ipr << w << "\t\t\t\t" << ipr / norm / (double)N << "\t\t\t\t" << entropy / norm << "\t\t\t\t" << r / norm << endl;
		stout << " \t\t--> w = " << w << " - in time : " << tim_s(start) << "s" << std::endl;
	}
	scaling_ipr.close();
}

/// <summary>
/// Compares matrix elements between disorder and not and saves it to file
/// </summary>
void isingUI::ui::compare_matrix_elements(op_type op, int k_alfa, int k_beta, int p_alfa, int p_beta, int x_alfa, int x_beta) {
	string name = IsingModel_sym::set_info(this->L, this->J, this->g, this->h, k_alfa, p_alfa, x_alfa, { "k","p","x" });
	name = name + "k_ab=" + to_string(k_alfa) + ":" + to_string(k_beta) + "p_ab=" + to_string(p_alfa) + ":" + to_string(p_beta) + "x_ab=" + to_string(x_alfa) + ":" + to_string(x_beta);
	std::ofstream file("results" + string(kPathSeparator) + "matrix_comparison" + name + ".dat");

	// disorder
	auto model = std::make_unique<IsingModel_disorder>(L, J, 0, g, 0, h, 0);
	model->diagonalization();
	file << "WHOLE SIZE = " << model->get_hilbert_size() << endl;
	// symmetries
	std::unique_ptr<IsingModel_sym> alfa = std::make_unique<IsingModel_sym>(L, J, g, h, 0, 1, 1);
	std::unique_ptr<IsingModel_sym> beta = std::make_unique<IsingModel_sym>(L, J, g, h, L / 2.0, 0, 0);
	file << "ALFA SECTOR SIZE = " << alfa->get_hilbert_size() << endl;
	file << "BETA SECTOR SIZE = " << alfa->get_hilbert_size() << endl;
	alfa->diagonalization();
	beta->diagonalization();

	// maps
	auto map_alfa = mapping_sym_to_original(0, model->get_hilbert_size() - 1, *alfa, *model);
	auto map_beta = mapping_sym_to_original(0, model->get_hilbert_size() - 1, *beta, *model);
	file << "alfa sector size for nondegenerate mapping is : " << map_alfa.size() << endl;
	file << "beta sector size for nondegenerate mapping is : " << map_alfa.size() << endl;

	file << "ALPHA SECTOR : " + alfa->get_info({ "L","J","g","h" }) + ", BETA SECTOR :" + beta->get_info({ "L","J","g","h" }) + "\n\n";
	file << " - - - - - - SAME SECTORS - - - - - - \n" << " - - - - > FOR ALFA - ALFA: \n" << "ENERGY ALPHA |('.'|)" << "\t\t" << \
		"ENERGY ALFA (/'.')/" << "\t\t" << "<alfa|SIGMA_X|alfa>" << "\t\t" << "<non_sym|SIGMA_X|non_sym>" << "\t\t" << "DIFFERENCE" << endl;
	for (auto& element : map_alfa) {
		for (auto& t : map_alfa) {
			cpx A = av_operator(element.first, t.first, *alfa, *alfa, op);
			cpx B = model->av_sigma_x(element.second, t.second);
			stout << alfa->get_eigenEnergy(element.first) << "\t\t" << \
				alfa->get_eigenEnergy(t.first) << "\t\t" << real(A) << "\t\t" << real(B) << "\t\t" << real(abs(A) - abs(B)) << endl;
		}
	}
	file << "\n - - - - > FOR BETA - BETA: \n" << "ENERGY BETA |('.'|)" << "\t\t" << \
		"ENERGY BETA (/'.')/" << "\t\t" << "<beta|SIGMA_X|beta>" << "\t\t" << "<non_sym|SIGMA_X|non_sym>" << "\t\t" << "DIFFERENCE" << endl;
	for (auto& element : map_beta) {
		for (auto& t : map_beta) {
			cpx A = av_operator(element.first, t.first, *beta, *beta, op);
			cpx B = model->av_sigma_x(element.second, t.second);
			file << beta->get_eigenEnergy(element.first) << "\t\t" << \
				beta->get_eigenEnergy(t.first) << "\t\t" << real(A) << "\t\t" << real(B) << "\t\t" << real(abs(A) - abs(B)) << endl;
		}
	}
	file << "\n\n - - - - - - DIFFERENT SECTORS - - - - - - \n" << " - - - - > FOR ALFA - BETA: \n" << "ENERGY ALPHA |('.'|)" << "\t\t" << \
		"ENERGY BETA (/'.')/" << "\t\t" << "<alfa|SIGMA_X|beta>" << "\t\t" << "<non_sym_a|SIGMA_X|non_sym_b>" << "\t\t" << "DIFFERENCE" << endl;
	for (auto& element : map_alfa) {
		for (auto& t : map_beta) {
			cpx A = im * av_operator(element.first, t.first, *alfa, *beta, op);
			A += conj(im * av_operator(t.first, element.first, *beta, *alfa, op));
			A *= 0.5i;
			cpx B = model->av_sigma_x(element.second, t.second);
			stout << alfa->get_eigenEnergy(element.first) << "\t\t" << \
				beta->get_eigenEnergy(t.first) << "\t\t" << real(A) << "\t\t" << real(B) << "\t\t" << real(abs(A) - abs(B)) << endl;
		}
	}
	file.close();
}

/// <summary>
///
/// </summary>
void isingUI::ui::check_dist_other_sector() {
	double step = 1e-3;
	double min = 0.0;
	double max = 0.1;
	const int size = abs(max - min) / step + 1;
	std::vector<double> E_sym = v_1d<double>();
	std::vector<double> av_sig;
	this->mu = 300;
	for (int k = 0; k < L; k++) {
		if (k == 0 || k == this->L / 2.) {
			for (int p = 0; p <= 1; p++) {
				int x_max = (this->h != 0) ? 0 : 1;
				for (int x = 0; x <= x_max; x++) {
					auto Hamil = std::make_unique<IsingModel_sym>(L, J, g, h, k, p, x, boundary_conditions);
					Hamil->diagonalization();
					vec t = Hamil->get_eigenvalues();
					E_sym.insert(E_sym.end(), std::make_move_iterator(t.begin() + Hamil->E_av_idx - u64(this->mu / 2.)), \
						std::make_move_iterator(t.begin() + Hamil->E_av_idx + u64(this->mu / 2.)));
					for (int i = Hamil->E_av_idx - mu / 2.; i < Hamil->E_av_idx + mu / 2.; i++)
						av_sig.push_back(Hamil->av_sigma_x(i, i, { 0 }));
				}
			}
		}
		else {
			int x_max = (this->h != 0) ? 0 : 1;
			for (int x = 0; x <= x_max; x++) {
				auto Hamil = std::make_unique<IsingModel_sym>(L, J, g, h, k, 1, x, boundary_conditions);
				Hamil->diagonalization();
				vec t = Hamil->get_eigenvalues();
				E_sym.insert(E_sym.end(), std::make_move_iterator(t.begin() + Hamil->E_av_idx - u64(this->mu / 2.)), \
					std::make_move_iterator(t.begin() + Hamil->E_av_idx + u64(this->mu / 2.)));
				for (int i = Hamil->E_av_idx - mu / 2.; i < Hamil->E_av_idx + mu / 2.; i++)
					av_sig.push_back(Hamil->av_sigma_x(i, i, { 0 }));
			}
		}
	}
	auto p = sort_permutation(E_sym, [](const double a, const double b) {
		return a < b;
		});
	apply_permutation(E_sym, p);
	apply_permutation(av_sig, p);
	arma::vec dist(size);
	for (u64 k = 0; k < av_sig.size() - 1; k++) {
		double value = abs(av_sig[k + 1] - av_sig[k]);
		setDistElem(dist, min, step, value);
	}
	normalise_dist(dist, min, max);
	save_to_file(this->saving_dir, "ProbDistSpecRapSigmaXAllSectors" + IsingModel_sym::set_info(L, J, g, h, 0, 0, 0, { "k", "p", "x" }), \
		arma::linspace(min, max, size), dist);
}

/// <summary>
///
/// </summary>
void isingUI::ui::fidelity(std::initializer_list<int> symetries) {
	const auto start = std::chrono::high_resolution_clock::now();
	std::vector<int> sym = { 0, 1, 1 };
	for (int i = 0; i < 3; i++)
		if (i < symetries.size())
			sym[i] = *(symetries.begin() + i);
	auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h, sym[0], sym[1], sym[2], this->boundary_conditions);
	alfa->diagonalization();
	stout << " \t\t--> finished diagonalizing 1st for " << alfa->get_info() << " - in time : " << tim_s(start) << "\nTotal time : " << tim_s(start) << "s\n";

	this->mu = 0.5 * alfa->get_hilbert_size();
	std::ofstream file(this->saving_dir + "Fidelity" + alfa->get_info({}) + ".dat");
	std::unique_ptr<IsingModel_sym> beta;
	arma::vec X = arma::logspace(-3, 1, 100);
	stout << alfa->get_hilbert_size() << endl;
	for (auto& de : X) {
		const auto start_loop = std::chrono::high_resolution_clock::now();
		beta.reset(new IsingModel_sym(this->L, this->J, this->g, this->h + de, sym[0], sym[1], sym[2], this->boundary_conditions));
		beta->diagonalization();
		const long int E_min = beta->E_av_idx - mu / 2.;
		const long int E_max = beta->E_av_idx + mu / 2.;
		double fidel = 0, entropy = 0;
		int counter = 0;
		for (int k = E_min; k < E_max; k++) {
			fidel += abs(overlap(*beta, *alfa, k, k));
			entropy += alfa->information_entropy(k, *beta, E_min, E_max);
			counter++;
		}
		fidel /= double(counter);
		entropy /= (double)counter;
		file << de << "\t\t" << fidel << "\t\t" << entropy << endl;
		stout << de << "\t\t" << fidel << "\t\t" << entropy << endl;
		//perturbative_stat_sym((de >=0.1)? de / 25. : de / 10., -1.0, 1.0, de, *alfa, *beta);
		//stout << "\t\t\t\t - - - - - - finished perturbation = " << de << " in : " << tim_s(start_loop) << " s" << "\nTotal time : " << tim_s(start) << "s\n";
	}
	stout << " \t\t--> finished fidelity for " << alfa->get_info() << " - in time : " << tim_s(start) << "s\n";
	file.close();
}

/// <summary>
///
/// </summary>
/// <param name="k"></param>
/// <param name="p"></param>
/// <param name="x"></param>
void isingUI::ui::size_scaling_sym(int k, int p, int x) {
	using namespace std::chrono;
	auto start = std::chrono::high_resolution_clock::now();

	const int L_max = this->L + this->Ln * this->Ls;
	auto beta = std::make_unique<IsingModel_sym>(2, J, g, h, k, p, x, boundary_conditions);
	std::ofstream farante(this->saving_dir + "IprScaling" + beta->get_info({ "L" }) + ".dat");
	std::ofstream fikolo(this->saving_dir + "SpectrumRapScalingSigmaX" + beta->get_info({ "L" }) + ".dat");

	for (int Lx = this->L; Lx <= L_max; Lx += this->Ls) {
		stout << "\n\n------------------------------Doing L = " << Lx << "------------------------------\n";
		auto alfa = std::make_unique<IsingModel_sym>(Lx, J, g, h, k, p, x, boundary_conditions);
		u64 N = alfa->get_hilbert_size();
		if (N <= 0) continue;
		alfa->diagonalization();
		this->mu = 0.4 * N;
		const long int E_min = alfa->E_av_idx - mu / 2.;
		const long int E_max = alfa->E_av_idx + mu / 2.;
		stout << " \t\t--> finished diagonalizing for " << alfa->get_info() << " - in time : " << \
			double(duration_cast<milliseconds>(duration(high_resolution_clock::now() - start)).count()) / 1000.0 << "s" << std::endl;

		// average sigma_x operator at first site
		std::ofstream sigx(this->saving_dir + "SigmaX" + alfa->get_info() + ".dat");
		vec av_sigma_x(mu, fill::zeros);
		for (int i = E_min; i < E_max; i++) {
			const int idx = i - (alfa->E_av_idx - mu / 2.);
			av_sigma_x(idx) = alfa->av_sigma_x(i, i, { 0 });
			sigx << alfa->get_eigenEnergy(idx) / double(Lx) << "\t\t" << av_sigma_x(idx) << endl;
		}
		sigx.close();

		//outliers and prob distribution for sigma_x
		vec r_sigma_x(mu - 1);
#pragma omp parallel for
		for (long int i = E_min; i < E_max - 1; i++) {
			const int idx = i - (alfa->E_av_idx - mu / 2.);
			r_sigma_x(idx) = abs(av_sigma_x(idx + 1) - av_sigma_x(idx));
		}

		vec outliers = statistics_average(r_sigma_x, 4);
		fikolo << Lx << "\t\t" << outliers.t();
		stout << " \t\t--> finished outliers for " << alfa->get_info() << " - in time : " << tim_s(start) << "s" << std::endl;

		probability_distribution(this->saving_dir, "ProbDistSpecRapSigmaX" + alfa->get_info(), r_sigma_x);
		probability_distribution(this->saving_dir, "ProbDistSigmaX" + alfa->get_info(), data_fluctuations(av_sigma_x));
		stout << " \t\t--> finished prob dist for " << alfa->get_info() << " - in time : " << tim_s(start) << "s" << std::endl;

		// eigenlevel statistics and prob distribution
		//vec r = alfa->eigenlevel_statistics_with_return();
		//probability_distribution(this->saving_dir, "ProbDistGap" + alfa->get_info(), r, 0, 1.0, 0.01);

		// ipr & info entropy
		double ipr = 0, ent = 0, r = 0;
		int counter = 0;
		for (int i = E_min; i < E_max; i++) {
			ipr += alfa->ipr(i);
			ent += alfa->information_entropy(i);
			r += alfa->eigenlevel_statistics(i, i + 1);
			counter++;
		}
		farante << Lx << "\t\t" << ipr / double(counter * N) << "\t\t" << ent / double(counter) << "\t\t" << r / double(counter) << endl;
		stout << " \t\t--> finished farante for " << alfa->get_info() << " - in time : " << tim_s(start) << "s" << std::endl;
	}
	fikolo.close();
	farante.close();
	stout << " - - - - - - FINISHED SIZE SCALING for:\nk = " << k << ", p = " << p << ", x = " << x << "\tIN : " << tim_s(start) << " seconds - - - - - - " << endl;
}
/// <summary>
///
/// </summary>
/// <param name="k"></param>
/// <param name="p"></param>
/// <param name="x"></param>
void isingUI::ui::parameter_sweep_sym(int k, int p, int x)
{
	const auto start = std::chrono::high_resolution_clock::now();
	std::string info = IsingModel_sym::set_info(L, J, g, h, k, p, x, { "h","g" });
	//std::ofstream pertGaussMap(this->saving_dir + "PertGaussianityMap" + info + ".dat");
	std::ofstream farante(this->saving_dir + "IprScalingMap" + info + ".dat");
	//std::ofstream kurt(this->saving_dir + "Moments" + info + ".dat");
	//kurt << "g\t\th\t\tSigmaX_kurtosis\tSigmaX_binder\tSigmaX_stddev\tSigmaZ_nnn_kurtosis\tSigmaZ_nnn_binder\tSigmaZ_nnn_stadev" << endl;
	farante << "g\t\th\t\tipr\t\t\tinformation entropy\tr\t\ttypical susceptiblity" << endl;
	//pertGaussMap << "g" << "\t\t" << "h" << "\t\t" << "perturbation" << "\t\t" << "kurtosis sig_x" << "\t\t" << "kurtosis energy" << "\n";

	// parameters for loops
	const double gmax = 0.8 + this->gn * this->gs;
	const double hmax = 0.4 + this->hn * this->hs;
	std::unique_ptr<IsingModel_sym> alfa;
	int counter_g = 0;
	for (double gx = 0.8; gx < gmax; gx += this->gs) {
		int counter_h = 0;
		for(auto& hx : arma::logspace(-4, 1, 100)) {
		//for (double hx = 0.4; hx < hmax; hx += this->hs) {
			//if (hx == 1.0) this->hs = 0.36;
			//else if (hx > 1.8) this->hs = 0.2;
			//else this->hs = 0.025;

			info = IsingModel_sym::set_info(L, J, gx, hx, k, p, x, {});
			//std::ofstream pertGauss(this->saving_dir + "PertGaussianity" + info + ".dat");
			//pertGauss << "perturbation" << "\t\t" << "kurtosis sig_x" << "\t\t" << "kurtosis energy" << "\n";

			const auto start_loop = std::chrono::high_resolution_clock::now();
			stout << "\n\n------------------------------ Doing : g = " << gx << ", h = " << hx << "------------------------------\n";

			// make model
			alfa.reset(new IsingModel_sym(this->L, this->J, gx, hx, k, p, x, this->boundary_conditions));
			stout << " \t\t--> finished creating model for " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
			alfa->diagonalization();
			const u64 N = alfa->get_hilbert_size();
			this->mu = 0.25 * N;

			const long int E_min = alfa->E_av_idx - mu / 2.;
			const long int E_max = alfa->E_av_idx + mu / 2.;
			stout << " \t\t--> finished diagonalizing for " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";

			//average sigma_x operator at first site prob dist
			//vec av_sigma_x(mu, fill::zeros);
			////vec av_sigma_z_nnn(mu, fill::zeros);
			//for (int i = E_min; i < E_max; i++) {
			//	const int idx = i - E_min;
			//	av_sigma_x(idx) = alfa->av_sigma_x(i, i, { 0 });
			//	//av_sigma_z_nnn(idx) = alfa->av_sigma_z(i, i, { 0,2 });
			//}
			//
			//arma::vec distSigmaX = probability_distribution_with_return(data_fluctuations(av_sigma_x));
			////arma::vec distSigmaZ_nnn = probability_distribution_with_return(data_fluctuations(av_sigma_z_nnn));
			//
			//kurt << gx << "\t\t" << hx << "\t\t" << binder_cumulant(distSigmaX) << "\t\t" << kurtosis_diff(distSigmaX) << "\t\t" << endl;// << binder_cumulant(distSigmaZ_nnn) << "\t\t"\
			//	<< kurtosis_diff(distSigmaZ_nnn) << "\t\t" << arma::stddev(distSigmaZ_nnn) << endl;
			////if (counter_h % 5 == 0)
			//save_to_file(this->saving_dir, "ProbDistSigmaX" + alfa->get_info(), arma::linspace(-1.0, 1.0, distSigmaX.size()), distSigmaX);

			//stout << " \t\t--> finished prob dist of sigma_x for " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "s\n";

			// ipr & info entropy
			double ipr = 0;
			double ent = 0;
			double r = 0;
			double typ_susc = 0;
			int counter = 0;
			for (int i = E_min; i < E_max; i++) {
				ipr += alfa->ipr(i);
				ent += alfa->information_entropy(i);
				r += alfa->eigenlevel_statistics(i, i + 1);
				double susc = 0;
				for (int j = E_min; j < E_max; j++) {
					if (j != i)	{
						const cpx overlap = alfa->av_sigma_z(i, j, { 1, 2 });
						const double omega = alfa->get_eigenEnergy(j) - alfa->get_eigenEnergy(i);
						susc += abs(overlap) * abs(overlap) / (omega * omega);
					}
				}
				typ_susc += log((double)L * susc);
				counter++;
			}
			farante << gx << "\t\t" << hx << "\t\t" << ipr / double(counter * N) << "\t\t" << ent / double(counter) << \
				"\t\t" << r / double(counter) << "\t\t" << exp(typ_susc / double(counter)) << endl;

//			if (((abs(hx - 1.0) <= 0.1 || abs(hx - 1.5) <= 0.2 || abs(hx - 2.0) <= 0.2 || abs(hx - 2.8) <= 0.2) && counter_h % 1 == 0) || counter_h % 10 == 0) {
//				double pert_step = 0.05;
//				double pert_max = 0.30;
//				if (hx == 1.0 || (hx >= 1.36 && hx <= 1.64)) {
//					pert_max = 0.3;
//					pert_step = 0.01;
//				}
//				for (double pert = 0.005; pert <= pert_max; pert += pert_step) {
//					auto vec = this->perturbative_stat_sym(pert, *alfa, gx, hx);
//					pertGauss << pert << "\t\t" << vec[0] << "\t\t" << vec[1] << "\n";
//					pertGaussMap << gx << "\t\t" << hx << "\t\t" << pert << "\t\t" << vec[0] << "\t\t" << vec[1] << "\n";
//				}
//				pertGauss.flush();
//				pertGaussMap.flush();
//				//outliers and prob distribution for sigma_x
//				vec r_sigma_x(mu - 1);
//				vec r_sigma_z_nnn(mu - 1);
//#pragma omp parallel for
//				for (long int i = E_min; i < E_max - 1; i++) {
//					const int idx = i - E_min;
//					r_sigma_x(idx) = abs(av_sigma_x(idx + 1) - av_sigma_x(idx));
//					//r_sigma_z_nnn(idx) = abs(av_sigma_z_nnn(idx + 1) - av_sigma_z_nnn(idx));
//				}
				//probability_distribution(this->saving_dir, "ProbDistSpecRapSigmaX" + alfa->get_info(), r_sigma_x);
				//probability_distribution(this->saving_dir, "ProbDistSpecRapSigmaZNNN" + alfa->get_info(), r_sigma_z_nnn);
				//pertGauss.close();
			//}
			counter_h++;
			stout << "\t\t\t--> finished calculating ETH params for " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
		}
		farante << endl << endl;
		counter_g++;
	}
	stout << " - - - - - - FINISHED PARAMETER SCALING for:\nk = " << k << ", p = " << p << ", x = " << x << "IN : " << tim_s(start) << "s\n";
	farante.close();
	//kurt.close();
	//pertGaussMap.close();
}
/// <summary>
///
/// </summary>
/// <param name="alfa_sym"></param>
/// <param name="beta_sym"></param>
/// <param name="omega_dist"></param>
/// <param name="omega_gauss_max"></param>
/// <param name="energy_constraint"></param>
void isingUI::ui::matrix_elements_stat_sym(double min, double max, double step, double omega_dist, int omega_gauss_max, \
	double energy_constraint, int energy_num, std::initializer_list<int> alfa_sym, std::initializer_list<int> beta_sym) const
{
	const auto start = std::chrono::high_resolution_clock::now();
	// in order k, x, p because we can skip p for most cases, x as well but in different order
	std::vector<int> alfa_syms = { 0, 1, 1 };
	std::vector<int> beta_syms = { 0, 1, 1 };
	for (int i = 0; i < 3; i++) {
		if (i < alfa_sym.size())
			alfa_syms[i] = *(alfa_sym.begin() + i);
		if (i < beta_sym.size())
			beta_syms[i] = *(beta_sym.begin() + i);
	}
	auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h, alfa_syms[0], alfa_syms[2], alfa_syms[1], this->boundary_conditions);
	auto beta = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h, beta_syms[0], beta_syms[2], beta_syms[1], this->boundary_conditions);
	alfa->diagonalization();
	beta->diagonalization();
	stout << " \t\t--> finished diagonalizing for " << alfa->get_info() << "\tand\t" << beta->get_info() << " - in time : " << tim_s(start) << "s\n";
	const u64 Na = alfa->get_hilbert_size();
	const u64 Nb = beta->get_hilbert_size();
	const int size = static_cast<int>(abs(max - min) / step + 1);
	const double omega_step = 0.025;
	const int omega_num = static_cast<int>(omega_gauss_max / omega_step);
	// files
	ofstream gaussianity(this->saving_dir + "gaussianity_for_alfa=" + alfa->get_info() + "_beta=" + beta->get_info() + ".dat");
	gaussianity << "omega" << "\t\t" << "sig_x" << "\t\t" << "sig_x_ex\n";
	ofstream op_prob_dist(this->saving_dir + "operatorProbDist_for_alfa=" + alfa->get_info() + "_beta=" + beta->get_info() + ".dat");
	op_prob_dist << "P(|O_ab|)" << "\t\t" << "sig_x" << "\t\t" << "sig_x_ex" << "\t\t" << "sig_z_nnn" << "\t\t" << "sig_z_nnn_ex" << "\t\t" << "spin_flip\n";

	// operators
	auto sig_x = IsingModel_sym::sigma_x;
	auto sig_z = IsingModel_sym::sigma_z;
	auto spin_fl = IsingModel_sym::spin_flip;
	// sigma _z
	vec d_sigma_z_nnn(size, arma::fill::zeros), d_sigma_z_nnn_ex(size, arma::fill::zeros);	// for making the probability distribution of close sigma_zs matrix elements
	//double min_sz_nnn = 0, min_sz_nnn_ex = 0;													// minimas of sigma_z proba distributions
	// sigma _x
	vec d_sigma_x(size, arma::fill::zeros), d_sigma_x_ex(size, arma::fill::zeros);
	v_1d<double> av_sigma_x(omega_num, 0);													// average sigma_x in small omega bucket
	v_1d<double> av_sigma_x_ex(omega_num, 0);												// average sigma_x_ex in small omega bucket
	v_1d<double> av_sigma_x2(omega_num, 0);													// average sigma_x in small omega bucket
	v_1d<double> av_sigma_x2_ex(omega_num, 0);												// average sigma_x_ex in small omega bucket0);
	// spin flip
	//vec d_sigma_flip(size, arma::fill::zeros);

	double E_av = alfa->get_eigenEnergy(alfa->E_av_idx);
	E_av = (abs(E_av) <= 1e-3) ? 1.0 : E_av;
	const int alfa_max_E = static_cast<int>(alfa->E_av_idx + energy_num / 2.0);
	const int beta_max_E = static_cast<int>(beta->E_av_idx + energy_num / 2.0);
	v_1d<int> counter(omega_num, 0);
	for (int a = alfa->E_av_idx - energy_num / 2.; a < alfa_max_E; a++) {
		if (a < 0 || a >= Na) continue;
		const double Ea = alfa->get_eigenEnergy(a);
		for (int b = beta->E_av_idx - energy_num / 2.; b < beta_max_E; b++) {
			if (b < 0 || b >= Nb) continue;
			const double Eb = beta->get_eigenEnergy(b);
			if (abs((Eb + Ea) / (2.0 * E_av) - 1.0) > energy_constraint / 2.0) continue;
			const double omega = abs(Ea - Eb);
			const auto omega_bucket = static_cast<int>(omega / omega_step);

			// sigma_x
			double sigma_x;
			double sigma_x_ex;
			if ((omega_bucket < omega_num) || (omega < omega_dist)) {
				sigma_x = abs(av_operator(a, b, *alfa, *beta, sig_x, { 0 }));
				sigma_x_ex = abs(av_operator(a, b, *alfa, *beta, sig_x));
			}
			// sigma x to averages
			if (omega_bucket < omega_num) {
				av_sigma_x[omega_bucket] += sigma_x;
				av_sigma_x_ex[omega_bucket] += sigma_x_ex;
				av_sigma_x2[omega_bucket] += sigma_x * sigma_x;
				av_sigma_x2_ex[omega_bucket] += sigma_x_ex * sigma_x_ex;
				counter[omega_bucket]++;
			}
			// check omega constraint
			if (omega < omega_dist) {
				// sigma_z
				setDistElem(d_sigma_z_nnn, min, step, abs(av_operator(a, b, *alfa, *beta, sig_z, { 0, 2 })));
				setDistElem(d_sigma_z_nnn_ex, min, step, abs(av_operator(a, b, *alfa, *beta, sig_z, 2)));
				// sigma_x already calculated
				setDistElem(d_sigma_x, min, step, sigma_x);
				setDistElem(d_sigma_x_ex, min, step, sigma_x_ex);
				// spin flip
				//setDistElem(d_sigma_flip, min, step, 0.5 * abs(av_operator(a, b, *alfa, *beta, spin_fl) + conj(av_operator(b, a, *beta, *alfa, spin_fl))));
				//rescale_Distribution(d_sigma_flip,min_sf,step,spin_flip_ex);
			}
		}
		stout << " \t\t--> finished energy: " << alfa->get_eigenEnergy(a) << " in 1st sector -  in time : " << tim_s(start) << "s\n";
	}
	stout << " \t\t--> finished calculating off-diagonal elements - in time : " << tim_s(start) << "s\n";
	d_sigma_x = normalise_dist(d_sigma_x, min, max); d_sigma_x_ex = normalise_dist(d_sigma_x_ex, min, max);
	d_sigma_z_nnn = normalise_dist(d_sigma_z_nnn, min, max); d_sigma_z_nnn_ex = normalise_dist(d_sigma_z_nnn_ex, min, max);
	// d_sigma_flip = normalise_dist(d_sigma_flip, min, max);
	for (int i = 0; i < size; i++)
		op_prob_dist << i * step + min << "\t\t" << d_sigma_x(i) << "\t\t" << d_sigma_x_ex(i) << \
		"\t\t" << d_sigma_z_nnn(i) << "\t\t" << d_sigma_z_nnn_ex(i) << "\t\t\n";// << d_sigma_flip(i) << "\n";
	for (int i = 0; i < omega_num; i++) {
		gaussianity << i * omega_step << "\t\t" << (av_sigma_x2[i] / av_sigma_x[i]) / av_sigma_x[i] * counter[i] << \
			"\t\t" << (av_sigma_x2_ex[i] / av_sigma_x_ex[i]) * av_sigma_x_ex[i] * counter[i] << std::endl;
	}
	gaussianity.close();
	op_prob_dist.close();
}

/// <summary>
///
/// </summary>
/// <param name="dist_step"></param>
/// <param name="pert"></param>
v_1d<double> isingUI::ui::perturbative_stat_sym(double pert, double gx, double hx)
{
	clk::time_point start = std::chrono::high_resolution_clock::now();
	auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, gx, hx, \
		this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
	alfa->diagonalization();
	this->mu = 0.5 * alfa->get_hilbert_size();
	long int size = 1 + 3.322 * log(this->mu);
	int n_av = 20;
	double pert_change = pert / (double)n_av;
	double pert_min = pert - n_av / 2. * pert_change;
	double pert_max = pert + n_av / 2. * pert_change;

	v_1d<double> kurtos;
	double _min_sig_x = 0, _max_sig_x = 0;
	double _min_E = 0, _max_E = 0;
	arma::vec dis_sig_x(size, arma::fill::zeros);
	arma::vec dis_delta_E(size, arma::fill::zeros);
	int counter = 0;
	for (double deps = pert_min; deps <= pert_max; deps += pert_change) {
		auto beta = std::make_unique<IsingModel_sym>(this->L, this->J, gx + deps, hx + deps, \
			this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		beta->diagonalization();
		const long int E_min = beta->E_av_idx - mu / 2.;
		const long int E_max = beta->E_av_idx + mu / 2.;
		// operators
		arma::vec delta_sig_x(mu, arma::fill::zeros);
		arma::vec delta_E(mu, arma::fill::zeros);
		for (long int i = E_min; i < E_max; i++) {
			const long int idx = i - E_min;
			delta_sig_x(idx) = beta->av_sigma_x(i, i, { 0 }) - alfa->av_sigma_x(i, i, { 0 });
			delta_E(idx) = beta->get_eigenEnergy(i) - alfa->get_eigenEnergy(i);
		}
		dis_sig_x += probability_distribution_with_return(delta_sig_x, size);
		dis_delta_E += probability_distribution_with_return(delta_E, size);
		counter++;
		_min_sig_x = arma::min(delta_sig_x);
		_max_sig_x = arma::max(delta_sig_x);
		_min_E = arma::min(delta_E);
		_max_E = arma::max(delta_E);
	}
	dis_sig_x /= double(counter);
	dis_delta_E /= double(counter);
	kurtos.push_back(kurtosis(dis_sig_x));
	kurtos.push_back(kurtosis(dis_delta_E));

	ofstream dis_op(this->saving_dir + "perturbationOperatorsDist" + alfa->get_info() + ",pert=" + to_string_prec(pert, 4) + ".dat");
	dis_op << "P(O_aa)\t\tsig_x\n";
	ofstream dis_E(this->saving_dir + "perturbationEnergyDiffDist" + alfa->get_info() + ",pert=" + to_string_prec(pert, 4) + ".dat");
	const double step_sig_x = abs(_max_sig_x - _min_sig_x) / (double)size;
	const double step_E = abs(_max_E - _min_E) / (double)size;
	for (int i = 0; i < size; i++) {
		dis_op << step_sig_x * i + _min_sig_x << "\t\t" << dis_sig_x(i) << "\t\t" << gaussian(step_sig_x * i + _min_sig_x, 0.0, arma::stddev(dis_sig_x)) << "\n";
		dis_E << step_E * i + _min_E << "\t\t" << dis_delta_E(i) << "\t\t" << gaussian(step_E * i + _min_E, 0.0, arma::stddev(dis_delta_E)) << "\n";
	}
	dis_op.close(); dis_E.close();
	stout << "\t\t\t\t - - - - - - FINISHED perturbation = " << pert << " IN : " << tim_s(start) << " seconds - -----" << endl;
	return kurtos;
}

/// <summary>
///
/// </summary>
/// <param name="dist_step"></param>
/// <param name="min"></param>
/// <param name="max"></param>
/// <param name="pert"></param>
/// <param name="alfa"></param>
/// <param name="beta"></param>
/// <returns></returns>
std::vector<double> isingUI::ui::perturbative_stat_sym(double pert, IsingModel_sym& alfa, double gx, double hx) {
	clk::time_point start = std::chrono::high_resolution_clock::now();
	long int size = 1 + 3.322 * log(this->mu);
	const int n_av = 4;

	double pert_change = pert / (double)n_av;
	if (pert <= 0.005)
		pert_change = 0.000025;
	else if (pert >= 0.01 && pert <= 0.099)
		pert_change = 0.0002;
	else
		pert_change = 0.001;
	double pert_min = pert - n_av / 2. * pert_change;
	double pert_max = pert + n_av / 2. * pert_change;

	// saving sigma_x from alfa not to recalculate it again
	v_1d<double> sigma_x;
	for (long int i = 0; i < alfa.get_hilbert_size(); i++) {
		sigma_x.push_back(alfa.av_sigma_x(i, i, { 0 }));
	}

	v_1d<double> kurtos;
	double _min_sig_x = INT_MAX, _max_sig_x = INT_MIN;
	double _min_E = INT_MAX, _max_E = INT_MIN;
	arma::vec dis_sig_x(size, arma::fill::zeros);
	arma::vec dis_delta_E(size, arma::fill::zeros);
	int counter = 0;
	for (double deps = pert_min; deps <= pert_max; deps += pert_change) {
		auto beta = std::make_unique<IsingModel_sym>(this->L, this->J, gx + deps, hx + deps, \
			this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		beta->diagonalization();
		const long int E_min = beta->E_av_idx - mu / 2.;
		const long int E_max = beta->E_av_idx + mu / 2.;
		// operators
		arma::vec delta_sig_x(mu, arma::fill::zeros);
		arma::vec delta_E(mu, arma::fill::zeros);
		for (long int i = E_min; i < E_max; i++) {
			const long int idx = i - E_min;
			delta_sig_x(idx) = beta->av_sigma_x(i, i, { 0 }) - sigma_x[i];
			delta_E(idx) = beta->get_eigenEnergy(i) - alfa.get_eigenEnergy(i);
			// save mins and maxes
			if (delta_sig_x(idx) > _max_sig_x) _max_sig_x = delta_sig_x(idx);
			else if (delta_sig_x(idx) < _min_sig_x) _min_sig_x = delta_sig_x(idx);
			if (delta_E(idx) > _max_E) _max_E = delta_sig_x(idx);
			else if (delta_E(idx) < _min_E) _min_E = delta_sig_x(idx);
		}
		dis_sig_x += probability_distribution_with_return(delta_sig_x, size);
		dis_delta_E += probability_distribution_with_return(delta_E, size);
		counter++;
	}
	dis_sig_x /= double(counter);
	dis_delta_E /= double(counter);
	kurtos.push_back(kurtosis(dis_sig_x));
	kurtos.push_back(kurtosis(dis_delta_E));

	ofstream dis_op(this->saving_dir + "perturbationOperatorsDist" + alfa.get_info() + ",pert=" + to_string_prec(pert, 4) + ".dat");
	dis_op << "P(O_aa)\t\tsig_x\n";
	ofstream dis_E(this->saving_dir + "perturbationEnergyDiffDist" + alfa.get_info() + ",pert=" + to_string_prec(pert, 4) + ".dat");
	const double step_sig_x = abs(_max_sig_x - _min_sig_x) / (double)size;
	const double step_E = abs(_max_E - _min_E) / (double)size;
	for (int i = 0; i < size; i++) {
		dis_op << step_sig_x * i + _min_sig_x << "\t\t" << dis_sig_x(i) << "\n";
		dis_E << step_E * i + _min_E << "\t\t" << dis_delta_E(i) << "\n";
	}
	dis_op.close(); dis_E.close();
	stout << "\t\t\t\t - - - - - - FINISHED perturbation = " << pert << " IN : " << tim_s(start) << " seconds - -----" << endl;
	return kurtos;
}
/// <summary>
///
/// </summary>
/// <param name="dist_step"></param>
/// <param name="min"></param>
/// <param name="max"></param>
/// <param name="pert"></param>
/// <param name="alfa"></param>
/// <param name="beta"></param>
/// <returns></returns>
std::vector<double> isingUI::ui::perturbative_stat_sym(double dist_step, double min, double max, double pert, IsingModel_sym& alfa, IsingModel_sym& beta) {
	const double E_dist_step = 5 * dist_step;
	const int size = static_cast<int>(abs(max - min) / dist_step);
	const int E_size = static_cast<int>(abs(max - min) / E_dist_step);
	// operators
	v_1d<double> kurtos;
	vec dis_sig_x(size, arma::fill::zeros);
	vec dis_delta_E(E_size, arma::fill::zeros);
	this->mu = 0.5 * alfa.get_hilbert_size();
	for (int i = alfa.E_av_idx - mu / 2.; i <= alfa.E_av_idx + mu / 2.; i++) {
		const double delta_sig_x = abs(beta.av_sigma_x(i, i, { 0 }) - alfa.av_sigma_x(i, i, { 0 }));
		const double delta_E = abs(beta.get_eigenEnergy(i) - alfa.get_eigenEnergy(i));
		setDistElem(dis_sig_x, min, dist_step, delta_sig_x);
		setDistElem(dis_delta_E, min, E_dist_step, delta_E);
		//stout << std::setprecision(4) << i / (double)alfa->get_hilbert_size() * 100 << "%" << endl;
	}
	dis_sig_x = normalise_dist(dis_sig_x, min, max);
	dis_delta_E = normalise_dist(dis_delta_E, min, max);
	kurtos.push_back(kurtosis_diff(dis_sig_x));
	kurtos.push_back(kurtosis_diff(dis_delta_E));

	ofstream dis_op(this->saving_dir + "perturbationOperatorsDist" + alfa.get_info() + ",pert=" + to_string_prec(pert, 4) + ".dat");
	dis_op << "P(O_aa)\t\tsig_x\n";

	ofstream dis_E(this->saving_dir + "perturbationEnergyDiffDist" + alfa.get_info() + ",pert=" + to_string_prec(pert, 4) + ".dat");

	for (int i = 0; i < size; i++) {
		dis_op << dist_step * i + min << "\t\t" << dis_sig_x(i) << "\t\t" << gaussian(dist_step * i + min, 0.0, arma::stddev(dis_sig_x)) << "\n";
		if (i < E_size) dis_E << E_dist_step * i + min << "\t\t" << dis_delta_E(i) << "\t\t" << gaussian(E_dist_step * i + min, 0.0, arma::stddev(dis_sig_x)) << "\n";
	}

	dis_op.close(); dis_E.close();
	return kurtos;
}

/// <summary>
///
/// </summary>
void isingUI::ui::make_sim()
{
	using namespace std::chrono;
	clk::time_point start = std::chrono::high_resolution_clock::now();

	
	//disorder();
	std::string info = IsingModel_sym::set_info(L, J, g, h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"L", "h"});
	std::ofstream farante(this->saving_dir + "AGPtotal_dis" + info + ".dat");
	std::ofstream scaling(this->saving_dir + "AGPsize_DELETE" + info + ".dat");
	farante << std::setprecision(6) << std::scientific;
	scaling << std::setprecision(6) << std::scientific;
	//farante << "hx\t\t\tL=10 susceptibiltiy\t\tAGP,\t\t11 susceptibiltiy\tAGP,\t\t12 susceptibiltiy\tAGP,\t\t13 susceptibiltiy\tAGP,\t\t14 susceptibiltiy\tAGP\n";
	farante << "hx\t\tL=8 AGP,\t\tL=9 AGP,\t\tL=10 AGP,\t\t11 AGP,\t\t12 AGP,\t\t13 AGP,\t\t14 AGP\n";
	//auto params = arma::logspace(-4, 1, 100);
	auto params = arma::linspace(1.0, 3.0, 41);
	//std::vector<double> params = { 0, 1e-4, 2.5e-4, 5e-4, 1e-3, 2.5e-3, 5e-3 , 1e-2, 5e-2, 1e-1};
	//std::vector<double> params = { 1e-0 };
	for (auto& hx : params) {
		farante << hx << "\t\t";
		scaling << "\"h = " + to_string_prec(hx, 2) << "\"" << endl;
		for (L = 10; L <= 12; L++) {
			const auto start_loop = std::chrono::high_resolution_clock::now();
			auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, 0, this->g, 0, hx, 1e-2, 0);
			//auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, hx, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym);
			stout << " \t\t--> finished creating model for " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
			alfa->diagonalization();
			const u64 N = alfa->get_hilbert_size();
			this->mu = 0.5 * N;
			const double mu2 = double(L) / double(N);
			//const double mu2 = std::log(N) / log(2) / double(N);
			const long int E_min = alfa->E_av_idx - mu / 2.;
			const long int E_max = alfa->E_av_idx + mu / 2.;
			stout << " \t\t--> finished diagonalizing for " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
			double typ_susc = 0;
			double AGP = 0;
			int counter = 0;
			const double omega_max = alfa->get_eigenEnergy(N - 1) - alfa->get_eigenEnergy(0);
			arma::mat mat_elem(N, N, arma::fill::zeros);
			for (int i = 0; i < N; i++) {
				double susc = 0;
//#pragma omp parallel for reduction(+: susc, AGP)
				for (int j = i + 1; j < N && j != i; j++) { // if want AGP fast do sum over j > i
					const cpx overlap = (double)L * alfa->av_sigma_x(i, j); // operators defined as 1/L*sum(..) !!!!!!!
					const double nominator = abs(overlap) * abs(overlap);
					const double omega_ij = alfa->get_eigenEnergy(j) - alfa->get_eigenEnergy(i);
					const double denominator = omega_ij * omega_ij + mu2 * mu2;
					susc += nominator / (omega_ij * omega_ij);
					AGP += omega_ij * omega_ij * nominator / (denominator * denominator);
					
					mat_elem(i, j) = nominator;
				}
				if (i < E_max && i > E_min)
					typ_susc += log((double)L * susc);
				counter++;
			}
			AGP *= 2.;
			//farante << exp(typ_susc / double(counter)) << "\t\t" << AGP / double(counter * L) << "\t\t";
			farante << AGP / double(counter * L) << "\t\t";
			scaling << L << "\t\t" << exp(typ_susc / double(counter)) << "\t\t" << AGP / double(counter * L) << "\t\t" << std::log(N) / log(2) << "\n";
			farante.flush();
			scaling.flush();

			// Spectral Function
			std::ofstream reponse_fun(this->saving_dir + "ResponseFunctionSigmaX" + alfa->get_info({ "h" }) + ",h=" + to_string_prec(hx, 5) + ".dat");
			reponse_fun.flush();
			v_1d<long int> idx_beta, idx_alfa;		// indices satysfying first condition in sum
			v_1d<double> energy_diff;				// energy differnece(omega) of the above indices
			const double tol = 0.005 * L;;
			for (int i = 0; i < N; i++)
				for (int j = i + 1; j < N && j != i; j++) 
					if (abs((alfa->get_eigenEnergy(j) + alfa->get_eigenEnergy(i)) / 2. - alfa->get_eigenEnergy(alfa->E_av_idx)) < tol / 2.) {
						idx_alfa.push_back(i);
						idx_beta.push_back(j);
						energy_diff.push_back(abs(alfa->get_eigenEnergy(j) - alfa->get_eigenEnergy(i)));
					}
			auto permut = sort_permutation(energy_diff, [](const double a, const double b) {return a < b; });
			apply_permutation(energy_diff, permut);
			apply_permutation(idx_beta, permut);
			apply_permutation(idx_alfa, permut);
			long int size = energy_diff.size();
			v_1d<int> Mx = { 150, 300, 700 };
			int M = Mx[L - 10];
			long int bucket_num = int(size / (double)M);
			for (int k = 0; k < bucket_num; k++){
				double element = 0;
				double omega = 0;
				for (long int p = k * M; p < (k + 1) * M; p++) {
					element += mat_elem(idx_alfa[p], idx_beta[p]);
					omega += energy_diff[p];
				}
				reponse_fun << omega / (double)M << "\t\t" << 2. * element / double(M) * (double)N << endl;
				reponse_fun.flush();
			}
			double element = 0;
			double omega = 0;
			counter = 0;
			for (long int p = bucket_num * M; p < size; p++) {
				element += mat_elem(idx_alfa[p], idx_beta[p]);
				omega += energy_diff[p];
				counter++;
			}
			reponse_fun << omega / (double)counter << "\t\t" << 2. * element / double(counter) * (double)N << endl;
			reponse_fun.flush();
			reponse_fun.close();
		}
		farante << endl;
		scaling << endl << endl;
	}
	farante.close();
	scaling.close();

	//parameter_sweep_sym(0, 1, 1);
	//this->h = 1.7; this->g = 0.8;
	//check_dist_other_sector();
	// 
	//matrix_elements_stat_sym(0, 1, 5e-4, 0.01, 10, 1e-1, 500, { 1, 1, 1 }, { 1, 1, 1 });
	//matrix_elements_stat_sym(0, 1, 0.001, 0.1, 10, 0.05, 100, { 1, 1, 1 }, { 2, 1, 1 });

	stout << " - - - - - - FINISHED CALCULATIONS IN : " << tim_s(start) << " seconds - - - - - - " << endl;						// simulation end
}