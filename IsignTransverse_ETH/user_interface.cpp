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
	NO_OVERFLOW(std::vector<std::string> tmp(argc - 1, "");)
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
	value = static_cast<T>(stod(value_str));
}

// ----------------------------------------------------------------------------- ISING MODEL -----------------------------------------------------------------------------

/// <summary>
/// prints all options from the user interface (if not chosen the default is taken)
/// </summary>
void isingUI::ui::printAllOptions() const {
	stout << "------------------------------CHOSEN OPTIONS:" << std::endl;
	std::string opName = IsingModel_disorder::opName(this->op, this->site);
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
		<< "bucket size = " << this->mu << std::endl;

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
}
// ----------------------------------------------------------------------------- Connected with the parser
/// <summary>
/// Setting the default parameters for the Ising model
/// </summary>
void isingUI::ui::set_default() {
	using namespace isingUI;
	this->saving_dir = "." + std::string(kPathSeparator) + "results" + std::string(kPathSeparator);		// directory for the result files to be saved into
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
	this->symmetries.p_sym = 0;
	this->symmetries.x_sym = 0;

	this->realisations = 100;
	this->site = 0;
	this->op = 0;
	this->fun = INT_MAX;
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
	this->parseModel((int)input.size(), input);													// parse input from CMD directly
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
		"	2 -- ABC\n"
		"-s site to act with local operators (default 0)"
		"-op flag to choose operator: "
		"	0 -- Sz_i-local"
		"	1 -- Sx_i-local"
		"	2 -- Hi"
		"	3 -- Sz_q"
		"	4 -- Sx_q"
		"	5 -- Hq"
		"	-- i or q are set to site (flag -s); (default 0 -> local Sz)"
		""
		"-fun choose function to start calculations"
		"	0 -- time evolution (and spectral functions) with disorder --> LIOMsdisorder(...)"
		"	1 -- AGPs for small disorder (-m=0) as function of h for input operator from -op flag\n\t\t\t\tSET: -L, -Ln, -Ls, -h, -hn, -hs, -op, -w(default=0.01)"
		"	2 -- AGPs for small disorder (-m=0) as function of g for input operator from -op flag\n\t\t\t\tSET: -L, -Ln, -Ls, -g, -gn, -gs, -op, -w(default=0.01)"
		"	3 -- LIOMs of Tranverse-field-Ising-model --> TFIsingLIOMs()"
		"	4 -- relaxation times from integrated spectral function for operator -op flag on site -s flag"
		" def -- in make_sim space for user to write function; desifned for non-builtin behavior"
		""
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
	if (this->boundary_conditions > 2) this->set_default_msg(this->boundary_conditions, choosen_option.substr(1), \
		"max boundary condition is 2", table);

	//choose site
	choosen_option = "-s";
	this->set_option(this->site, argv, choosen_option);

	//choose operator
	choosen_option = "-op";
	this->set_option(this->op, argv, choosen_option);

	// model
	choosen_option = "-m";
	this->set_option(this->m, argv, choosen_option);
	if (this->m > 1) this->set_default_msg(this->m, choosen_option.substr(1), \
		"max model number is 1", table);

	//choose function
	choosen_option = "-fun";
	this->set_option(this->fun, argv, choosen_option);
	if (this->fun < 2) this->m = 0;

	// buckets
	choosen_option = "-mu";
	this->set_option(this->mu, argv, choosen_option);
	choosen_option = "-r";
	this->set_option(this->realisations, argv, choosen_option);

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
		break;
	default:
		str_model += "PBC" + std::string(kPathSeparator);
		break;
	}

	std::string folder = saving_dir + str_model;
	NO_OVERFLOW(
	if (!argv[argc - 1].empty() && argc % 2 != 0) {
		// only if the last command is non-even
		folder = argv[argc - 1] + str_model;
		if (fs::create_directories(folder) || fs::is_directory(folder))						// creating the directory for saving the files with results
			this->saving_dir = folder;																// if can create dir this is is
	}
	else {
		if (fs::create_directories(folder) || fs::is_directory(folder))						// creating the directory for saving the files with results
			this->saving_dir = folder;																// if can create dir this is is
	})

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
				v_1d<std::string> temp_str = v_1d<std::string>(t.size(), "k=" + std::to_string(k) + ",x=" + to_string(x) + ",p=*");
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
	stout << "symmetry sector\t\tEnergy sym\t\tEnergy total\t\tdifference\t\t <n|X|n>" << endl;
	for (int k = 0; k < E_dis.size(); k++) {
		stout << symms[k] << "\t\t" << E_sym[k] << "\t\t" << E_dis(k) << "\t\t" << E_sym[k] - E_dis(k) << endl;
	}
}

/// <summary>
/// Handle disorder simulation
/// </summary>
void isingUI::ui::disorder() {
	auto start = std::chrono::system_clock::now();
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
				double(duration_cast<milliseconds>(duration(system_clock::now() - start)).count()) / 1000.0 << "s" << std::endl;
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
		this->mu = long(0.5 * N);

		const long int E_min = Hamil->E_av_idx - long(mu / 2);
		const long int E_max = Hamil->E_av_idx + long(mu / 2);

		vec av_sigma_x(mu, arma::fill::zeros);
		for (u64 k = E_min; k < E_max; k++) {
			const long long idx = k - E_min;
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
				const long long idx = k - E_min;
				av_sigma_x(idx) = Hamil->av_sigma_x(k, k, { 0 });
			}
			fluct = data_fluctuations(av_sigma_x);

			prob_dist += probability_distribution_with_return(fluct);
			prob_dist_GOE += probability_distribution_with_return(Hamil->eigenlevel_statistics_with_return());

			// average in middle spectrum
			for (long int f = E_min; f < E_max; f++) {
				ipr += Hamil->ipr(f);
				entropy += Hamil->information_entropy(f);
				NO_OVERFLOW(r += Hamil->eigenlevel_statistics(f, f + 1));
			}
			if (k % 5 == 0) stout << " \t\t--> " << k << " realisation - in time : " << tim_s(start) << "s" << std::endl;
		}
		stout << "--> finished loop over realisations for w = " << w << " <--\n";

		kurt << w << "\t\t" << binder_cumulant(prob_dist) << "\t\t" << kurtosis_diff(prob_dist) << "\t\t" << endl;
		save_to_file(this->saving_dir, "ProbDistSigmaX" + Hamil->get_info(), arma::linspace(_min, _max, prob_dist.size()), prob_dist);

		save_to_file(this->saving_dir, "ProbDistGap" + Hamil->get_info(), arma::linspace(0, 1, prob_dist_GOE.size()), prob_dist_GOE);

		NO_OVERFLOW(double norm = realisations * mu;)
		scaling_ipr << w << "\t\t\t\t" << ipr / norm / (double)N << "\t\t\t\t" << entropy / norm << "\t\t\t\t" << r / norm << endl;
		stout << " \t\t--> w = " << w << " - in time : " << tim_s(start) << "s" << std::endl;
	}
	scaling_ipr.close();
}

/// <summary>
/// Compares matrix elements between disorder and not and saves it to file
/// </summary>
void isingUI::ui::compare_matrix_elements(op_type op, int k_alfa, int k_beta, int p_alfa, int p_beta, int x_alfa, int x_beta) {
	string name = IsingModel_sym::set_info(this->L, this->J, this->g, this->h, k_alfa, p_alfa, x_alfa, { "k","p","x"});
	//name = name + "k_ab=" + to_string(k_alfa) + ":" + to_string(k_beta) + "p_ab=" + to_string(p_alfa) + ":" + to_string(p_beta) + "x_ab=" + to_string(x_alfa) + ":" + to_string(x_beta);
	std::ofstream file(this->saving_dir + "matrix_comparison" + name + ".dat");
	// disorder
	auto model = std::make_unique<IsingModel_disorder>(L, J, 0, g, 0, h, 0);
	model->diagonalization();
	sp_cx_mat opMatrix = model->create_operator({ op }, std::vector<int>({ 0 }));
	file << "WHOLE SIZE = " << model->get_hilbert_size() << endl;
	// symmetries
	std::unique_ptr<IsingModel_sym> alfa = std::make_unique<IsingModel_sym>(L, J, g, h, k_alfa, p_alfa, x_alfa);
	std::unique_ptr<IsingModel_sym> beta = std::make_unique<IsingModel_sym>(L, J, g, h, k_beta, p_beta, x_beta);
	sp_cx_mat opMatrix_alfa = alfa->create_operator({ op }, std::vector<int>({ 0 }));
	sp_cx_mat opMatrix_beta = beta->create_operator({ op }, std::vector<int>({ 0 }));
	file << "ALFA SECTOR SIZE = " << alfa->get_hilbert_size() << endl;
	file << "BETA SECTOR SIZE = " << beta->get_hilbert_size() << endl;
	alfa->diagonalization();
	beta->diagonalization();

	// maps
	auto map_alfa = mapping_sym_to_original(0, model->get_hilbert_size() - 1, *alfa, *model);
	auto map_beta = mapping_sym_to_original(0, model->get_hilbert_size() - 1, *beta, *model);
	file << "alfa sector size for nondegenerate mapping is : " << map_alfa.size() << endl;
	file << "beta sector size for nondegenerate mapping is : " << map_beta.size() << endl;
	
	file << "ALPHA SECTOR : " + alfa->get_info({ "L","J","g","h" }) + "\n BETA SECTOR :" + beta->get_info({ "L","J","g","h" }) + "\n\n";
	file << " - - - - - - SAME SECTORS - - - - - - \n" << " - - - - > FOR ALFA - ALFA: \n" << "ENERGY ALPHA |('.'|)" << "\t\t" << \
		"ENERGY ALFA (/'.')/" << "\t\t" << "<alfa|SIGMA_X|alfa>" << "\t\t" << "<non_sym|SIGMA_X|non_sym>" << "\t\t" << "DIFFERENCE\t\tA/B" << endl;
	for (auto& element : map_alfa) {
		for (auto& t : map_alfa) {
			cpx A = arma::cdot(alfa->get_eigenState(element.first), opMatrix_alfa * alfa->get_eigenState(t.first));
			//cpx A = av_operator(element.first, t.first, *alfa, *alfa, op);
			cpx B = arma::dot(model->get_eigenState(element.second), opMatrix * model->get_eigenState(t.second));
			if (abs(abs(A) - abs(B)) >= 1e-14)
				file << alfa->get_eigenEnergy(element.first) << "\t\t\t\t\t" << \
				alfa->get_eigenEnergy(t.first) << "\t\t\t\t" << A << "\t\t\t\t" << B << \
				"\t\t\t\t" << abs(A) - abs(B) << "\t\t\t\t" << abs(A) / abs(B) << std::endl;
		}
	}
	file << "\n - - - - > FOR BETA - BETA: \n" << "ENERGY BETA |('.'|)" << "\t\t" << \
		"ENERGY BETA (/'.')/" << "\t\t" << "<beta|SIGMA_X|beta>" << "\t\t" << "<non_sym|SIGMA_X|non_sym>\t\tDIFFERENCE\t\tA/B" << endl;
	for (auto& element : map_beta) {
		for (auto& t : map_beta) {
			cpx A = arma::cdot(beta->get_eigenState(element.first), opMatrix_beta * beta->get_eigenState(t.first));
			//cpx A = av_operator(element.first, t.first, *beta, *beta, op);
			cpx B = arma::dot(model->get_eigenState(element.second), opMatrix * model->get_eigenState(t.second));
			if (abs(abs(A) - abs(B)) >= 1e-14)
				file << beta->get_eigenEnergy(element.first) << "\t\t\t\t\t" << \
				beta->get_eigenEnergy(t.first) << "\t\t\t\t" << A << "\t\t\t\t" << B << \
				"\t\t\t\t" << abs(A) - abs(B) << "\t\t\t\t" << abs(A) / abs(B) << std::endl;
		}
	}
	file << "\n\n - - - - - - DIFFERENT SECTORS - - - - - - \n" << " - - - - > FOR ALFA - BETA: \n" << "ENERGY ALPHA |('.'|)" << "\t\t" << \
		"ENERGY BETA (/'.')/" << "\t\t" << "<alfa|SIGMA_X|beta>" << "\t\t" << "<non_sym_a|SIGMA_X|non_sym_b>" << "\t\t" << "DIFFERENCE\t\tA/B" << endl;
	for (auto& element : map_alfa) {
		for (auto& t : map_beta) {
			//cpx A = arma::cdot(alfa->get_eigenState(element.first), opMatrix_beta * beta->get_eigenState(t.first));
			cpx A = 0.;//TODO: here add operator between sectors
			//cpx A = av_operator(element.first, t.first, *alfa, *beta, op);
			cpx B = arma::dot(model->get_eigenState(element.second), opMatrix * model->get_eigenState(t.second));
			if (abs(abs(A) - abs(B)) >= 1e-14)
				file << alfa->get_eigenEnergy(element.first) << "\t\t\t\t\t" << \
				beta->get_eigenEnergy(t.first) << "\t\t\t\t" << A << "\t\t\t\t" << B << \
				"\t\t\t\t" << abs(A) - abs(B) << "\t\t\t\t" << abs(A) / abs(B) << std::endl;
		}
	}
	file.flush();
	file.close();
}

/// <summary>
///
/// </summary>
void isingUI::ui::check_dist_other_sector() {
	double step = 1e-3;
	double min = 0.0;
	double max = 0.1;
	const int size = int(abs(max - min) / step) + 1;
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
					for (int i = Hamil->E_av_idx - mu / 2; i < Hamil->E_av_idx + mu / 2; i++)
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
				for (int i = Hamil->E_av_idx - mu / 2; i < Hamil->E_av_idx + mu / 2; i++)
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
	const auto start = std::chrono::system_clock::now();
	std::vector<int> sym = { 0, 1, 1 };
	for (int i = 0; i < 3; i++)
		if (i < symetries.size())
			sym[i] = *(symetries.begin() + i);
	auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h, sym[0], sym[1], sym[2], this->boundary_conditions);
	alfa->diagonalization();
	stout << " \t\t--> finished diagonalizing 1st for " << alfa->get_info() << " - in time : " << tim_s(start) << "\nTotal time : " << tim_s(start) << "s\n";

	this->mu = long(0.6 * alfa->get_hilbert_size());
	const long int E_min = alfa->E_av_idx - long(mu / 2);
	const long int E_max = alfa->E_av_idx + long(mu / 2);
	std::ofstream file(this->saving_dir + "Fidelity_1_" + alfa->get_info({}) + ".dat");
	std::unique_ptr<IsingModel_sym> beta;
	arma::vec X = arma::linspace(0.02, 2.02, 101);
	//arma::vec X = arma::logspace(-3, 1, 100);
	stout << alfa->get_hilbert_size() << endl;
	for (auto& de : X) {
		const auto start_loop = std::chrono::system_clock::now();
		beta.reset(new IsingModel_sym(this->L, this->J, this->g, this->h + de, sym[0], sym[1], sym[2], this->boundary_conditions));
		beta->diagonalization();
		double fidel = 0, entropy = 0;
		int counter = 0;
		double r = 0;
		for (int k = E_min; k < E_max; k++) {
			fidel += abs(overlap(*beta, *alfa, k, k));
			entropy += alfa->information_entropy(k, *beta, E_min, E_max);
			r += beta->eigenlevel_statistics(E_min, E_max);
			counter++;
		}
		fidel /= double(counter);
		entropy /= (double)counter;
		r /= double(counter);
		file << de << "\t\t" << fidel << "\t\t" << entropy << "\t\t" << r << endl;
		//stout << de << "\t\t" << fidel << "\t\t" << entropy << "\t\t" << r << endl;
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
	auto start = std::chrono::system_clock::now();

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
		this->mu = long(0.4 * N);
		const long int E_min = alfa->E_av_idx - long(mu / 2);
		const long int E_max = alfa->E_av_idx + long(mu / 2);
		stout << " \t\t--> finished diagonalizing for " << alfa->get_info() << " - in time : " << \
			double(duration_cast<milliseconds>(duration(system_clock::now() - start)).count()) / 1000.0 << "s" << std::endl;

		// average sigma_x operator at first site
		std::ofstream sigx(this->saving_dir + "SigmaX" + alfa->get_info() + ".dat");
		vec av_sigma_x(mu, fill::zeros);
		for (int i = E_min; i < E_max; i++) {
			const int idx = i - int(alfa->E_av_idx - mu / 2);
			av_sigma_x(idx) = alfa->av_sigma_x(i, i, { 0 });
			sigx << alfa->get_eigenEnergy(idx) / double(Lx) << "\t\t" << av_sigma_x(idx) << endl;
		}
		sigx.close();

		//outliers and prob distribution for sigma_x
		NO_OVERFLOW(vec r_sigma_x(mu - 1);)
#pragma omp parallel for
		for (long int i = E_min; i < E_max - 1; i++) {
			const int idx = i - int(alfa->E_av_idx - mu / 2);
			NO_OVERFLOW(r_sigma_x(idx) = abs(av_sigma_x(idx + 1) - av_sigma_x(idx)));
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
			NO_OVERFLOW(r += alfa->eigenlevel_statistics(i, i + 1));
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
	const auto start = std::chrono::system_clock::now();
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

			const auto start_loop = std::chrono::system_clock::now();
			stout << "\n\n------------------------------ Doing : g = " << gx << ", h = " << hx << "------------------------------\n";

			// make model
			alfa.reset(new IsingModel_sym(this->L, this->J, gx, hx, k, p, x, this->boundary_conditions));
			stout << " \t\t--> finished creating model for " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
			alfa->diagonalization();
			const u64 N = alfa->get_hilbert_size();
			this->mu = long(0.25 * N);

			const long int E_min = alfa->E_av_idx - long(mu / 2);
			const long int E_max = alfa->E_av_idx + long(mu / 2);
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
				NO_OVERFLOW(r += alfa->eigenlevel_statistics(i, i + 1));
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
	const auto start = std::chrono::system_clock::now();
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
	for (int a = alfa->E_av_idx - energy_num / 2; a < alfa_max_E; a++) {
		if (a < 0 || a >= Na) continue;
		const double Ea = alfa->get_eigenEnergy(a);
		for (int b = beta->E_av_idx - energy_num / 2; b < beta_max_E; b++) {
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
/// <param name="alfa"></param>
/// <param name="gx"></param>
/// <param name="hx"></param>
/// <returns></returns>
v_1d <double> isingUI::ui::perturbative_stat_sym(IsingModel_sym& alfa, double gx, double hx) {
	clk::time_point start = std::chrono::system_clock::now();
	const u64 N = alfa.get_hilbert_size();
	long int size = 1 + 3.322 * log(this->mu);
	long int sizeE = 1 + 3.322 * log(N);
	const int n_av = 0;
	const long int E_min = 0;// alfa.E_av_idx - mu / 2.;
	const long int E_max = (long)N;// alfa.E_av_idx + mu / 2.;
	// ANALITYCAL
	sp_cx_mat pertMatrix = alfa.create_operator({ IsingModel_sym::sigma_x, IsingModel_sym::sigma_z }) * sqrt(alfa.L);
	//pertMatrix /= sqrt(arma::trace(pertMatrix * pertMatrix) / double(N));
	cx_mat U = alfa.get_eigenvectors();
	cx_mat mat_elem = U.t() * pertMatrix * U;
	cpx order_2nd = 0.0, order_3rd = 0.0, order_4th = 0.0, AGP = 0;
	for (long int n = 0; n < N; n++) {
		double temp = 0;
		for (long int m = 0; m < N && m != n; m++) {
			double omega = alfa.get_eigenEnergy(n) - alfa.get_eigenEnergy(m);
			double elem = abs(mat_elem(n, m));
			temp += elem * elem / omega;
			AGP += elem * elem / (omega * omega);
		}
		cpx diag = mat_elem(n, n);
		order_2nd += diag * diag;
		order_3rd += diag * temp;
		order_4th += temp * temp / 4.0;
	}
	order_2nd /= double(N);
	order_3rd /= double(N);
	order_4th /= double(N);
#if defined(OPERATOR)
	arma::sp_cx_mat opMatrix = alfa.create_operator({ IsingModel_sym::sigma_x }, { 0 });
	cpx opVar = arma::trace(opMatrix * opMatrix) / double(N * N);
	arma::cx_vec sigma_x(mu, arma::fill::zeros);// = arma::diagvec(alfa.get_eigenvectors().t() * opMatrix * alfa.get_eigenvectors());
#pragma omp parallel for
	for (long int k = E_min; k < E_max; k++)
		sigma_x(k - E_min) = arma::cdot(alfa.get_eigenState(k), opMatrix * alfa.get_eigenState(k));
	stout << "\t\t\t\t\t - - - - - - FINISHED building operator IN : " << tim_s(start) << " seconds - -----" << endl;
	NO_OVERFLOW(double _min_sig_x = INT_MAX; double _max_sig_x = INT_MIN;)
	arma::vec dis_sig_x(size, arma::fill::zeros);
#endif
	std::ofstream mom_sig(this->saving_dir + "perturbationOperatorsMoments" + alfa.get_info() + ".dat");
	std::ofstream mom_E(this->saving_dir + "perturbationEnergyMoments" + alfa.get_info() + ".dat");
	mom_sig << "de\t\t\t\tmean\t\tvar\t\tkurtosis\t\tAnalitycal formula\n";
	mom_E	<< "de\t\t\t\tmean\t\tvar\t\tkurtosis\t\tAnalitycal formula\n";
	auto perturbation = arma::logspace(-4, 0, 30);
	//auto perturbation = arma::linspace(0.01, 0.39, 20);
	for (auto& pert : perturbation) {
		const double pert_change = pert / 50. / (n_av == 0 ? 1. : (double)n_av);
		const double pert_min = pert - n_av / 2. * pert_change;
		const double pert_max = pert + n_av / 2. * pert_change;
		NO_OVERFLOW(
			double _min_sig_x = INT_MAX; double _max_sig_x = INT_MIN;
			double _min_E = INT_MAX; double _max_E = INT_MIN;
		);
		arma::vec dis_delta_E(sizeE, arma::fill::zeros);
		//arma::vec dis_sig_x(size, arma::fill::zeros);
		//arma::vec dis_delta_E(size, arma::fill::zeros);
		int counter = 0;
		double mean = 0, variance = 0, kurtosis = 0;
		double mean_E = 0, variance_E = 0, kurtosis_E = 0;
		for (double deps = pert_min; deps <= pert_max; deps += pert_change) {
			auto beta = std::make_unique<IsingModel_sym>(alfa.L, alfa.J, gx + deps, hx + deps, \
				this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
			beta->diagonalization();
#if defined(OPERATOR)
			arma::vec delta_sig_x(mu, arma::fill::zeros);
#pragma omp parallel for
			for (long int i = E_min; i < E_max; i++) {
				const long int idx = i - E_min;
				arma::subview_col state = beta->get_eigenvectors().col(i);
				delta_sig_x(idx) = real(arma::cdot(state, opMatrix * state) - sigma_x(idx));
			}
			_max_sig_x = arma::max(delta_sig_x);
			_min_sig_x = arma::min(delta_sig_x);
			//dis_sig_x += probability_distribution_with_return(delta_sig_x, size);
#endif
			arma::vec delta_E(N, arma::fill::zeros);
#pragma omp parallel for
			for (long int i = 0; i < N; i++)
				delta_E(i) = beta->get_eigenEnergy(i) - alfa.get_eigenEnergy(i);
			_max_E = arma::max(delta_E);
			_min_E = arma::min(delta_E);
			//dis_delta_E += probability_distribution_with_return(delta_E, sizeE);

			counter++;
			mean		+= arma::mean(delta_sig_x);
			variance	+= arma::stddev(delta_sig_x);
			kurtosis	+= kurtosis_diff(delta_sig_x);
			mean_E		+= arma::mean(delta_E);
			variance_E	+= arma::stddev(delta_E);
			kurtosis_E	+= kurtosis_diff(delta_E);
			counter++;
		}
#if defined(OPERATOR)
		//dis_sig_x /= double(counter);
		//kurtos[1] /= double(counter);
		//const double step_sig_x = abs(_max_sig_x - _min_sig_x) / (double)size;
		//ofstream dis_op(this->saving_dir + "/perturbation data/OperatorDiff/perturbationOperatorsDist" + alfa.get_info() + ",pert=" + to_string_prec(pert, 4) + ".dat");
		//dis_op << "P(O_aa)\t\tsig_x\n";
		//for (int i = 0; i < size; i++)
		//	dis_op << step_sig_x * i + _min_sig_x << "\t\t" << dis_sig_x(i) << "\n";
		//dis_op.close();
#endif

		//dis_delta_E /= double(counter);
		//kurtos[0] /= double(counter);
		//ofstream dis_E(this->saving_dir + "/perturbation data/EnergyDiff/perturbationEnergyDiffDist" + alfa.get_info() + ",pert=" + to_string_prec(pert, 4) + ".dat");
		//const double step_E = abs(_max_E - _min_E) / (double)sizeE;
		//for (int i = 0; i < sizeE; i++)
		//	dis_E << step_E * i + _min_E << "\t\t" << dis_delta_E(i) << "\n";
		//dis_E.close();
		mom_sig << pert << "\t\t" << mean / double(counter) << "\t\t" << variance / double(counter) << "\t\t"\
			<< kurtosis / double(counter) << "\t\t" << real(2 * pert * pert * opVar * AGP) << std::endl;
		mom_sig.flush();
		cpx vardE_analitycal = order_2nd * std::pow(pert, 2) + order_3rd * std::pow(pert, 3) + order_4th * std::pow(pert, 4);
		mom_E << pert << "\t\t" << mean_E / double(counter) << "\t\t" << variance_E / double(counter) << "\t\t"\
			<< kurtosis_E / double(counter) << "\t\t" << vardE_analitycal.real() << std::endl;
		mom_E.flush();
		stout << "\t\t\t\t - - - - - - FINISHED perturbation = " << pert << " IN : " << tim_s(start) << " seconds - -----" << endl;
	}
	mom_E.close();
	mom_sig.close();
	return { 0,0 };
}
v_1d<double> isingUI::ui::perturbative_stat_sym(double pert, double gx, double hx)
{
	clk::time_point start = std::chrono::system_clock::now();
	auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, gx, hx, \
		this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
	alfa->diagonalization();
	this->mu = long(0.5 * alfa->get_hilbert_size());
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
		const long int E_min = beta->E_av_idx - long(mu / 2);
		const long int E_max = beta->E_av_idx + long(mu / 2);
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
std::vector<double> isingUI::ui::perturbative_stat_sym(double pert, IsingModel_sym& alfa, double gx, double hx) {
	clk::time_point start = std::chrono::system_clock::now();
	const u64 N = alfa.get_hilbert_size();
	long int size = 1 + 3.322 * log(this->mu);
	long int sizeE = 1 + 3.322 * log(N);
	const int n_av = 10;

	double pert_change = pert / 50. / (double)n_av;
	double pert_min = pert - n_av / 2. * pert_change;
	double pert_max = pert + n_av / 2. * pert_change;

	const long int E_min = alfa.E_av_idx - long(mu / 2);
	const long int E_max = alfa.E_av_idx + long(mu / 2);
	// saving sigma_x from alfa not to recalculate it again
	v_1d<double> kurtos(2, 0.0);
#if defined(OPERATOR)
	arma::sp_cx_mat opMatrix = alfa.create_operator({ IsingModel_sym::sigma_x }, { 0 });
	arma::cx_vec sigma_x(mu, arma::fill::zeros);// = arma::diagvec(alfa.get_eigenvectors().t() * opMatrix * alfa.get_eigenvectors());
#pragma omp parallel for
	for (long int k = E_min; k < E_max; k++) {
		NO_OVERFLOW(sigma_x(k - E_min) = arma::cdot(alfa.get_eigenState(k), opMatrix * alfa.get_eigenState(k)););
	}
	stout << "\t\t\t\t\t - - - - - - FINISHED building operator IN : " << tim_s(start) << " seconds - -----" << endl;
	NO_OVERFLOW(double _min_sig_x = INT_MAX; double _max_sig_x = INT_MIN;);
	arma::vec dis_sig_x(size, arma::fill::zeros);
#endif
	NO_OVERFLOW(double _min_E = INT_MAX; double _max_E = INT_MIN;);
	arma::vec dis_delta_E(sizeE, arma::fill::zeros);
	int counter = 0;
	for (double deps = pert_min; deps <= pert_max; deps += pert_change) {
		auto beta = std::make_unique<IsingModel_sym>(this->L, this->J, gx + deps, hx + deps, \
			this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		beta->diagonalization();
		#if defined(OPERATOR)
			arma::vec delta_sig_x(mu, arma::fill::zeros);
#pragma omp parallel for
			for (long int i = E_min; i < E_max; i++) {
				const long int idx = i - E_min;
				delta_sig_x(idx) = real(arma::cdot(beta->get_eigenState(i), opMatrix * beta->get_eigenState(i)) - sigma_x[idx]);
			}
			_max_sig_x = arma::max(delta_sig_x);
			_min_sig_x = arma::min(delta_sig_x);
			dis_sig_x += probability_distribution_with_return(delta_sig_x, size);
			kurtos[1] += kurtosis_diff(delta_sig_x);
		#endif
		arma::vec delta_E(N, arma::fill::zeros);
#pragma omp parallel for
		for (long int i = 0; i < N; i++) 
			delta_E(i) = beta->get_eigenEnergy(i) - alfa.get_eigenEnergy(i);
		_max_E = arma::max(delta_E);
		_min_E = arma::min(delta_E);
		dis_delta_E += probability_distribution_with_return(delta_E, sizeE);
		kurtos[0] += kurtosis_diff(delta_E);
		counter++;
	}
#if defined(OPERATOR)
	dis_sig_x /= double(counter);
	kurtos[1] /= double(counter);
	const double step_sig_x = abs(_max_sig_x - _min_sig_x) / (double)size;
	ofstream dis_op(this->saving_dir + "/perturbation data/OperatorDiff/perturbationOperatorsDist" + alfa.get_info() + ",pert=" + to_string_prec(pert, 4) + ".dat");
	dis_op << "P(O_aa)\t\tsig_x\n";
	for (int i = 0; i < size; i++) 
		dis_op << step_sig_x * i + _min_sig_x << "\t\t" << dis_sig_x(i) << "\n";
	dis_op.close();
#endif

	dis_delta_E /= double(counter);
	kurtos[0] /= double(counter);
	ofstream dis_E(this->saving_dir + "/perturbation data/EnergyDiff/perturbationEnergyDiffDist" + alfa.get_info() + ",pert=" + to_string_prec(pert, 4) + ".dat");
	const double step_E = abs(_max_E - _min_E) / (double)sizeE;
	for (int i = 0; i < sizeE; i++) 
		dis_E << step_E * i + _min_E << "\t\t" << dis_delta_E(i) << "\n";
	dis_E.close();

	stout << "\t\t\t\t - - - - - - FINISHED perturbation = " << pert << " IN : " << tim_s(start) << " seconds - -----" << endl;
	return kurtos;
}
std::vector<double> isingUI::ui::perturbative_stat_sym(double dist_step, double min, double max, double pert, IsingModel_sym& alfa, IsingModel_sym& beta) {
	const double E_dist_step = 5 * dist_step;
	const int size = static_cast<int>(abs(max - min) / dist_step);
	const int E_size = static_cast<int>(abs(max - min) / E_dist_step);
	// operators
	v_1d<double> kurtos;
	vec dis_sig_x(size, arma::fill::zeros);
	vec dis_delta_E(E_size, arma::fill::zeros);
	this->mu = long(0.5 * alfa.get_hilbert_size());
	for (int i = alfa.E_av_idx - mu / 2; i <= alfa.E_av_idx + mu / 2; i++) {
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

//--------------------------------------------------------------------- SPECTRAL PROPERTIES
template <typename _type> void isingUI::ui::spectralFunction(IsingModel<_type>& alfa, const arma::cx_mat& mat_elem, std::string name) {
	const u64 N = alfa.get_hilbert_size();
	v_1d<long int> idx_beta, idx_alfa;		// indices satysfying first condition in sum
	v_1d<double> energy_diff;				// energy differnece(omega) of the above indices
	const double tol = 0.05 * L;
	for (long int i = 0; i < N; i++) {
		for (long int j = 0; j < N && j != i; j++) {
			if (abs((alfa.get_eigenEnergy(j) + alfa.get_eigenEnergy(i)) / 2. - alfa.get_eigenEnergy(alfa.E_av_idx)) < tol / 2.) {
				idx_alfa.push_back(i);
				idx_beta.push_back(j);
				energy_diff.push_back(abs(alfa.get_eigenEnergy(j) - alfa.get_eigenEnergy(i)));
			}
		}
	}
	auto permut = sort_permutation(energy_diff, [](const double a, const double b) {return a < b; });
	apply_permutation(energy_diff, permut);
	apply_permutation(idx_beta, permut);
	apply_permutation(idx_alfa, permut);
	long int size = (int)energy_diff.size();
	//from L=12
	v_1d<int> Mx = v_1d<int>({ 100, 400, 700, 1200, 2000, 5000, 8000, 12000, 15000, 6000, 8000, 10000, 12000, 14000, 16000 });
	NO_OVERFLOW(int M = Mx[alfa.L - 8];);
	std::ofstream reponse_fun;
	openFile(reponse_fun, name + alfa.get_info({}) + "_M=" + std::to_string(M) + ".dat", ios::out);
	//long int M = std::pow(N, 0.75);
	long int bucket_num = int(size / (double)M);
	for (int k = 0; k < bucket_num; k++) {
		double element = 0;
		double omega = 0;
		for (long int p = k * M; p < (k + 1) * M; p++) {
			cpx overlap = mat_elem(idx_alfa[p], idx_beta[p]);
			element += abs(overlap * overlap);
			omega += energy_diff[p];
		}
		reponse_fun << omega / (double)M << "\t\t" << element / double(M) << endl;
		reponse_fun.flush();
	}
	double element = 0;
	double omega = 0;
	int counter = 0;
	for (long int p = bucket_num * M; p < size; p++) {
		cpx overlap = mat_elem(idx_alfa[p], idx_beta[p]);
		element += abs(overlap * overlap);
		omega += energy_diff[p];
		counter++;
	}
	reponse_fun << omega / (double)counter << "\t\t" << element / double(counter) << endl;
	reponse_fun.flush();
	reponse_fun.close();
}
template <typename _type> void isingUI::ui::integratedSpectralFunction(const IsingModel<_type>& alfa, const arma::cx_mat& mat_elem, std::string name) {
	const u64 N = alfa.get_hilbert_size();
	const double wH = alfa.mean_level_spacing_analytical();
	auto omegas = arma::logspace(std::floor(log10(wH)) - 1, 2, 300);
	std::ofstream reponse_fun;
	openFile(reponse_fun, name + alfa.get_info({}) + ".dat", ios::out);
	for (auto& w : omegas) {
		double overlap = 0.;
#pragma omp parallel for reduction(+: overlap)
		for (long int n = 0; n < N; n++) {
			for (long int m = 0; m < N; m++) {
				const double w_nm = abs(alfa.get_eigenEnergy(n) - alfa.get_eigenEnergy(m));
				if(w >= w_nm)
					overlap += abs(mat_elem(n, m) * conj(mat_elem(n, m)));
			}
		}
		overlap *= 1. / double(N);
		if (w == omegas(0)) printSeparated(reponse_fun, "\t", 12, true, w, overlap, wH );
		else				printSeparated(reponse_fun, "\t", 12, true, w, overlap);
	}
	reponse_fun.close();
}
template <typename _type> auto isingUI::ui::integratedSpectralFunction(const IsingModel<_type>& alfa, const arma::cx_mat& mat_elem, const arma::vec& omegas) -> arma::vec{
	const u64 N = alfa.get_hilbert_size();
	arma::vec intSpec(omegas.size());
	for (int i = 0; i < omegas.size(); i++) {
		const double w = omegas(i);
		double overlap = 0.;
#pragma omp parallel for reduction(+: overlap)
		for (long int n = 0; n < N; n++) {
			for (long int m = 0; m < N; m++) {
				const double w_nm = abs(alfa.get_eigenEnergy(n) - alfa.get_eigenEnergy(m));
				if (w >= w_nm)
					overlap += abs(mat_elem(n, m) * conj(mat_elem(n, m)));
			}
		}
		intSpec(i) = overlap / double(N);
	}
	return intSpec;
}

template <typename _type> void isingUI::ui::timeEvolution(const IsingModel<_type>& alfa, const arma::cx_mat& mat_elem, std::string name) {
	const u64 N = alfa.get_hilbert_size();
	const double tH = 1. / alfa.mean_level_spacing_analytical();

	std::ofstream tEvolution;
	openFile(tEvolution, this->saving_dir + "timeEvolution" + name + alfa.get_info({  }) + ".dat", ios::out);
	double norm_diag = 0;
#pragma omp parallel for reduction(+: norm_diag)
	for (long int k = 0; k < N; k++) {
		cpx temp = mat_elem(k, k);
		norm_diag += abs(temp * temp);
	}
	norm_diag /= double(N);
	const int t_max = (int)std::ceil(std::log(tH));
	auto times = arma::logspace(-2, t_max, 500);
	for (auto& t : times) {
		double overlap = 0.;
		if (t > 5 * tH) break;
#pragma omp parallel for reduction(+: overlap)
		for (long int n = 0; n < N; n++) {
			overlap += abs(mat_elem(n, n) * conj(mat_elem(n, n)));
			for (long int m = n + 1; m < N; m++) {
				const double w_nm = alfa.get_eigenEnergy(n) - alfa.get_eigenEnergy(m);
				overlap += 2. * abs(mat_elem(n, m) * conj(mat_elem(n, m))) * std::cos(w_nm * t);
			}
		}
		overlap *= 1. / double(N);
		//tEvolution << t << "\t\t" << overlap << "\t\t" << tH << std::endl;
		if (t == times(0)) printSeparated(tEvolution, "\t", 12, true, t,overlap,tH, norm_diag);
		else printSeparated(tEvolution, "\t", 12, true, t,overlap);
	}
	tEvolution.close();
}
template <typename _type> auto isingUI::ui::timeEvolution(const IsingModel<_type>& alfa, const arma::cx_mat& mat_elem, const arma::vec& times) -> std::pair<arma::vec, double> {
	const u64 N = alfa.get_hilbert_size();
	const double tH = 1. / alfa.mean_level_spacing_analytical();
	double LTA = 0;
#pragma omp parallel for reduction(+: LTA)
	for (long int k = 0; k < N; k++) {
		cpx temp = mat_elem(k, k);
		LTA += abs(temp * temp);
	}
	LTA /= double(N);
	arma::vec timeEv(times.size(), arma::fill::zeros);
	for (long int k = 0; k < times.size(); k++) {
		auto t = times(k);
		double overlap = 0.;
#pragma omp parallel for reduction(+: overlap)
		for (long int n = 0; n < N; n++) {
			overlap += abs(mat_elem(n, n) * conj(mat_elem(n, n)));
			for (long int m = n + 1; m < N; m++) {
				const double w_nm = alfa.get_eigenEnergy(n) - alfa.get_eigenEnergy(m);
				overlap += 2. * abs(mat_elem(n, m) * conj(mat_elem(n, m))) * std::cos(w_nm * t);
			}
		}
		overlap *= 1. / double(N);
		timeEv(k) = overlap;
	}
	return std::make_pair(timeEv, LTA);
}
void isingUI::ui::relaxationTimesFromFiles() {
	std::string s = this->op >= 2 ? "q" : "j";
	std::string dir = this->saving_dir + "RelaxationTimes" + kPSep;
	createDirs(dir);
	std::string op = IsingModel_disorder::opName(this->op, this->site);
	std::ofstream map_g, map_h;
	openFile(map_h, dir + "_h" + op + IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, { "h", "g" }) + ".dat", ios::out);
	//std::vector<double> gx_list = { 0.025, 0.05, 0.1, 0.15, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
	NO_OVERFLOW(auto gx_list = arma::linspace(this->g, this->g + this->gs * (this->gn - 1), this->gn));
	NO_OVERFLOW(auto hx_list = arma::linspace(this->h, this->h + this->hs * (this->hn - 1), this->hn));
	for (auto& gx : gx_list) {
		for (auto& hx : hx_list) {
			// read time-evolution data
			std::ifstream file;
			std::string name = op + IsingModel_disorder::set_info(this->L, this->J, this->J0, gx, this->g0, hx, this->w) + ".dat";
			std::string filename = this->saving_dir + "IntegratedResponseFunction" + kPSep + s + "=" + std::to_string(this->site) + kPSep + name;
			auto data = readFromFile(file, filename);
			file.close();
			if (data.empty()) continue;
			// integrate and find relax rate
			double wH = data[2](0);
			if (data[1](0) <= 0.5) {
				for (int k = 0; k < data[0].size(); k++) {
					if (data[1](k) >= 0.5) {
						printSeparated(map_h, "\t", 12, true, hx, gx, 1. / data[0](k), 1. / wH);
						break;
					}
				}
			}
		}
	}
	map_h.close();
	openFile(map_g, dir + "_g" + op + IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, { "h", "g" }) + ".dat", ios::out);
	for (auto& hx : hx_list) {
		for (auto& gx : gx_list) {
			// read time-evolution data
			std::ifstream file;
			std::string name = op + IsingModel_disorder::set_info(this->L, this->J, this->J0, gx, this->g0, hx, this->w) + ".dat";
			std::string filename = this->saving_dir + "IntegratedResponseFunction" + kPSep + s + "=" + std::to_string(this->site) + kPSep + name;
			auto data = readFromFile(file, filename);
			file.close();
			if (data.empty()) continue;
			// integrate and find relax rate
			double wH = data[2](0);
			if (data[1](0) <= 1) {
				for (int k = 0; k < data[0].size(); k++) {
					if (data[1](k) >= 0.5) {
						printSeparated(map_g, "\t", 12, true, hx, gx, 1. / data[0](k), 1. / wH);
						break;
					}
				}
			}
		}
	}
	map_g.close();
}
void isingUI::ui::intSpecFun_from_timeEvol() {
	std::ofstream map_g, map_h;
	openFile(map_g, this->saving_dir + "RelaxationTimes_g" + IsingModel_disorder::set_info(L, 1.0, 0.0, g, 0.0, h, 0.01, { "h","g" }) + ".dat", ios::out);
	openFile(map_h, this->saving_dir + "RelaxationTimes_h" + IsingModel_disorder::set_info(L, 1.0, 0.0, g, 0.0, h, 0.01, { "h","g" }) + ".dat", ios::out);
	std::vector<double> gx_list = { 0.05, 0.1, 0.15, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
	auto hx_list = arma::linspace(0.1, 1.3, 7);
	const double eta = 1e-3;
	for (auto& gx : gx_list) {
		for (auto& hx : hx_list) {

			// read time-evolution data
			std::ifstream file;
			std::string name = "SigmaZ_j=0" + IsingModel_disorder::set_info(L, 1.0, 0.0, gx, 0.0, hx, 0.01) + ".dat";
			std::string filename = this->saving_dir + "timeEvolution" + kPSep + name;
			auto data = readFromFile(file, filename);
			if (data.empty()) continue;
			file.close();
			// integrate and find relax rate
			double wH = 1. / data[2](0);
			auto omegas = arma::logspace(std::floor(std::log10(wH)) - 2, std::log10(1. / data[0](0)), data[0].size());
			std::ofstream fileInt;
			std::string dir = this->saving_dir + "IntegratedResponseFunction" + kPSep + "integrated" + kPSep;
			createDirs(dir);
			openFile(fileInt, dir + name, ios::out);
			int counter = 0;
			for (auto& w : omegas) {
				arma::vec toIntegrate = data[1] % arma::sin(w * data[0]) / (data[0] * pi) % arma::exp(-data[0] * eta);
				double intSpecFun = 2 * simpson_rule(data[0], toIntegrate); // even integral, only t \in (0,inf)
				if (w == omegas(0)) printSeparated(fileInt, "\t", 12, true, w, intSpecFun, wH);
				else			    printSeparated(fileInt, "\t", 12, true, w, intSpecFun);
				if (counter == 0 && intSpecFun >= 0.5) {
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

void isingUI::ui::adiabaticGaugePotential_dis(bool h_vs_g) {
	clk::time_point start = std::chrono::system_clock::now();
	auto s = h_vs_g ? "h" : "g"; // (now not inversed because separated h sweep) //inversed, cause exclude the other one
	std::string info = IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, { "L" }, ",");
	std::string dir = this->saving_dir + "AGP" + kPSep; createDirs(dir);
	std::ofstream farante;
	openFile(farante, dir + IsingModel_disorder::opName(this->op, this->site) + info + "_" + s + ".dat");
	farante << std::setprecision(6) << std::scientific;
	//std::ofstream scaling(this->saving_dir + "AGPsize_DELETE" + info + ".dat");
	//scaling << std::setprecision(6) << std::scientific;
	NO_OVERFLOW(auto params = h_vs_g ? arma::linspace(this->h, this->h + this->hs * this->hn, this->hn + 1)
		: arma::linspace(this->g, this->g + this->gs * this->gn, this->gn + 1));
	for (auto& x : params) {
		farante << x << "\t\t";
		//scaling << "\"h = " + to_string_prec(hx, 5) << "\"" << endl;
		for (int system_size = this->L; system_size < this->L + this->Ls * this->Ln; system_size += this->Ls) 
		{
			const auto start_loop = std::chrono::system_clock::now();
			std::unique_ptr<IsingModel_disorder> alfa;
			if(h_vs_g) alfa = std::make_unique<IsingModel_disorder>(system_size, this->J, this->J0, this->g, this->g0, x,		this->w);
			else	   alfa = std::make_unique<IsingModel_disorder>(system_size, this->J, this->J0, x	   , this->g0, this->h, this->w);
			stout << " \t\t--> finished creating model for " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
			auto opMat = alfa->chooseOperator(this->op, this->site);
			normaliseOp(opMat);
			std::function
				getValues = [&](double& AGP, double& typ_susc, int& counter)
			{
				stout << " \t\t--> next realisation " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
				const u64 N = alfa->get_hilbert_size();
				const double omegaH = alfa->mean_level_spacing_analytical();
				const double rescale = (double)N * omegaH * omegaH / (double)L;
				this->mu = long(0.5 * N);
				const double mu2 = double(L) / double(N);
				static long int E_min = alfa->E_av_idx - long(mu / 2);
				static long int E_max = alfa->E_av_idx + long(mu / 2);
				double typ_susc_local = 0;
				double AGP_local = 0;
				const mat U = alfa->get_eigenvectors();
				arma::cx_mat mat_elem = U.t() * opMat * U;
			#pragma omp parallel for reduction(+: AGP_local, typ_susc_local)
				for (long int i = 0; i < N; i++) {
					double susc = 0;
					for (long int j = 0; j < N && j != i; j++) {
						const double nominator = abs(mat_elem(i, j) * conj(mat_elem(i, j)));
						const double omega_ij = alfa->get_eigenEnergy(j) - alfa->get_eigenEnergy(i);
						const double denominator = omega_ij * omega_ij + mu2 * mu2;
						AGP_local += omega_ij * omega_ij * nominator / (denominator * denominator);
						susc	  += nominator / (omega_ij * omega_ij);
					}
					if (susc > 0 && (i > E_min && i < E_max))
						typ_susc_local += log(susc);
				}
				typ_susc += exp(typ_susc_local / double(mu));
				AGP += AGP_local / double(N);
				counter++;
			};

			double typ_susc = 0, AGP = 0;
			int counter = 0;
			if (this->realisations > 1)
				average_over_realisations<double&, double&, int&>(*alfa, true, getValues, AGP, typ_susc, counter);
			else {
				alfa->diagonalization();
				getValues(AGP, typ_susc, counter);
			}
			farante << typ_susc / double(counter) << "\t\t" << AGP / double(counter) << "\t\t";
			farante.flush();
			//scaling << L << "\t\t" << typ_susc << "\t\t" << AGP << "\t\t" << std::log(N) / log(2) << "\n";
			//scaling.flush();

		}
		farante << endl;
		//scaling << endl << endl;
	}
	farante.close();
	//scaling.close();
}
void isingUI::ui::adiabaticGaugePotential_sym(bool SigmaZ, bool avSymSectors) {
	clk::time_point start = std::chrono::system_clock::now();
	std::string info = (avSymSectors ? IsingModel_sym::set_info(L, J, g, h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym\
		, { "L", "h", "k", "p", "x" }) : IsingModel_sym::set_info(L, J, g, h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"L", "h"}));
	std::string dir = this->saving_dir + "AGP" + kPSep; createDirs(dir);
	std::ofstream farante(dir + (SigmaZ ? "SigZ" : "SigX") + info + ".dat");
	farante << std::setprecision(6) << std::scientific;
	//std::ofstream scaling(this->saving_dir + "AGPsize_DELETE" + info + ".dat");
	//scaling << std::setprecision(6) << std::scientific;
	farante << "hx\t\t\t susceptibiltiy\t\tAGP,\t\tsusceptibiltiy\tAGP,\t\t susceptibiltiy\tAGP,\t\t susceptibiltiy\tAGP,\t\t susceptibiltiy\tAGP\n";
	auto params = arma::linspace(0.02, 3, 150);
	for (auto& hx : params) {
		farante << hx << "\t\t";
		//scaling << "\"h = " + to_string_prec(hx, 5) << "\"" << endl;
		stout << "\nh = " << hx << "\t\t";
		for (int system_size = this->L; system_size < this->L + this->Ls * this->Ln; system_size += this->Ls) {
			const auto start_loop = std::chrono::system_clock::now();
			//auto alfa = std::make_unique<IsingModel_disorder>(system_size, this->J, 0, this->g, 0, hx, 0, 0);
			std::function
				getValues = [&](int k, int p, int x, double& AGP, double& typ_susc, int& counter)
			{
				auto alfa = std::make_unique<IsingModel_sym>(system_size, this->J, this->g, hx, k, p, x, this->boundary_conditions);
				stout << " \t\t--> finished creating model for " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
				alfa->diagonalization();
				stout << " \t\t--> finished diagonalizing for " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
				const u64 N = alfa->get_hilbert_size();
				const double omegaH = alfa->mean_level_spacing_analytical();
				const double rescale = (double)N * omegaH * omegaH / (double)L;
				this->mu = long(0.5 * N);
				const double mu2 = double(L) / double(N);
				//const double mu2 = std::log2(N) / double(N);
				static long int E_min = 0; // alfa->E_av_idx - mu / 2.;
				static long int E_max = (long)N; // alfa->E_av_idx + mu / 2.;
				double typ_susc_local = 0;
				double AGP_local = 0;
				int counter_tmp = 0;
				const cx_mat U = alfa->get_eigenvectors();
				arma::sp_cx_mat opMatrix = SigmaZ ? alfa->create_operator({ IsingModel_sym::sigma_z }) : \
					alfa->create_operator({ IsingModel_sym::sigma_x });
				cpx norm = arma::trace(opMatrix * opMatrix) / double(N);
				if (abs(norm) > 1e-12) opMatrix /= norm;
				arma::cx_mat mat_elem = U.t() * opMatrix * U;
				for (long int i = 0; i < N; i++) {
					double susc = 0;
					for (long int j = 0; j < N && j != i; j++) {
						const double nominator = abs(mat_elem(i, j) * conj(mat_elem(i, j)));
						const double omega_ij = alfa->get_eigenEnergy(j) - alfa->get_eigenEnergy(i);
						const double denominator = omega_ij * omega_ij + mu2 * mu2;
						AGP_local += omega_ij * omega_ij * nominator / (denominator * denominator);
						susc += nominator / (omega_ij * omega_ij);
					}
					if (susc > 1e-14)
						typ_susc_local += log(susc);
					counter_tmp++;
				}
				typ_susc += exp(typ_susc_local / double(counter_tmp));
				AGP += AGP_local / double(counter_tmp); // AGP /= L is neglected due to 1/L in operator definition
				counter++;
			};

			double typ_susc = 0, AGP = 0;
			int counter = 0;
			if (avSymSectors)
				loopSymmetrySectors<double&, double&, int&>(getValues, hx, system_size, AGP, typ_susc, counter);
			else
				getValues(this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, AGP, typ_susc, counter);
			farante << typ_susc / double(counter) << "\t\t" << AGP / double(counter) << "\t\t";
			farante.flush();
			//scaling << L << "\t\t" << typ_susc << "\t\t" << AGP << "\t\t" << std::log(N) / log(2) << "\n";
			//scaling.flush();

		}
		farante << endl;
		//scaling << endl << endl;
	}
	farante.close();
	//scaling.close();
}
void isingUI::ui::combineAGPfiles() {
	std::ifstream input;
	std::ofstream output;
	std::string line;
	std::string dir = "SPINONsynchro" + kPSep;
	for (double gx = 0.05; gx <= 0.20; gx += 0.05) {
		std::string gstr = "g=" + to_string_prec(gx, 2);
		openFile(output, dir + "SigmaZ_j=" + std::to_string(this->site) + ",J0=0.00," + gstr + ",g0=0.00,w=0.01.dat", std::ios::out);
		for (double h = 0.01; h < 3.01; h += 0.5)
		{
			std::string filename = dir + "SigmaZ_j=" + std::to_string(this->site) + ",J0=0.00," + gstr + ",g0=0.00,h=" + to_string_prec(h, 2) + ",w=0.01_h.dat";
			openFile(input, filename, std::ios::in);
			while (std::getline(input, line))
				output << line << "\n";
			input.close();
		}
		output.close();
	}
}

template <typename _type> std::pair<double, double> isingUI::ui::operator_norm(arma::sp_cx_mat& opMatrix, IsingModel<_type>& alfa) {
	const u64 N = alfa.get_hilbert_size();
	normaliseOp(opMatrix);
	const Mat<_type> U = alfa.get_eigenvectors();
	arma::cx_mat mat_elem = U.t() * opMatrix * U;
	double norm_diag = 0, norm_off = 0;
#pragma omp parallel for reduction(+: norm_diag, norm_off)
	for (long int k = 0; k < N; k++) {
		cpx temp = mat_elem(k, k);
		norm_diag += abs(temp * temp);
		for (long int m = k + 1; m < N; m++) {
			temp = mat_elem(k, m);
			norm_off += 2 * abs(temp * temp);
		}
	}
	return std::make_pair(norm_diag / double(N), norm_off / double(N));
}
template <typename _type> void isingUI::ui::energyEvolution(IsingModel<_type >& alfa) {
	std::ofstream ener(this->saving_dir + "Energies" + alfa.get_info({}) + ".dat");
	sp_cx_mat pertMatrix = alfa.create_operator({ IsingModel_sym::sigma_x, IsingModel_sym::sigma_z }) * sqrt(alfa.L);
	cx_mat U = alfa.get_eigenvectors();
	cx_mat mat_elem = U.t() * pertMatrix * U;
	const u64 N = alfa.get_hilbert_size();
	auto energyEvolution = [&](int i, double pert)->double {
		double second_order = 0;
		for (long int m = 0; m < N && m != i; m++) {
			const double omega = alfa.get_eigenEnergy(i) - alfa.get_eigenEnergy(m);
			const double elem = abs(mat_elem(i, m));
			second_order += elem * elem / omega;
		}
		cpx diag = mat_elem(i, i);
		return real(alfa.get_eigenEnergy(i) + diag * pert + pert * pert / 2. * second_order);
	};
	this->mu = 10;
	const u64 E_min = alfa.E_av_idx - mu / 2.;
	const u64 E_max = alfa.E_av_idx + mu / 2.;
	arma::vec perturb = arma::logspace(-4, -1, 50);
	for (auto& pert : perturb) {
		ener << pert << "\t\t";
		auto beta = std::make_unique<IsingModel_sym>(alfa.L, alfa.J, alfa.g + pert, alfa.h + pert, \
			this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		beta->diagonalization(true);
		for (long int i = E_min; i < E_max; i++) {
			ener << beta->get_eigenEnergy(i) << "\t\t" << energyEvolution(i, pert) << "\t\t";
		}
		ener << std::endl;
	}
	ener.close();
}

template <typename _type> void isingUI::ui::IsingLIOMs(IsingModel<_type>& alfa) {
	const int L = alfa.L;
	std::ofstream ssf;
	bool isopen = openFile(ssf, this->saving_dir + "StaticStructureFactor" + alfa.get_info({}) + ".dat", ios::out);
	crash(!isopen);
	std::vector<arma::sp_cx_mat> Sz;
	for (int j = 0; j < L; j++) {
		auto opMatrix = alfa.create_operator({ IsingModel_sym::sigma_z }, std::vector({ j }));
		Sz.push_back(opMatrix);
		if (j == 0) timeEvolution(alfa, opMatrix, "SigmaZ_j=" + std::to_string(j));
	}
	auto GS = alfa.get_eigenState(0);
	for (int l = 1; l <= L; l++) {
		double q = l * two_pi / (double)L;
		sp_cx_mat Sz_q = Sz[0];
		for (int j = 1; j < L; j++) {
			Sz_q += Sz[j] * std::exp(im * double(j) * q);
		}
		Sz_q /= sqrt(L);
		cpx Sq = arma::cdot(Sz_q * GS, Sz_q * GS);
		printSeparated(ssf, "\t", 10, true, q / (double)two_pi, real(Sq), imag(Sq));
		timeEvolution(alfa, Sz_q, "SigmaZ_q=" + std::to_string(l));
	}
	ssf.close();
}
void isingUI::ui::TFIsingLIOMs() {
	std::string timeDir = this->saving_dir + "timeEvolution" + kPSep + "TFIM" + kPSep + std::to_string(this->L) + kPSep;
	auto body = [&](const auto& alfa) {
		auto opMatrix = alfa.create_LIOMoperator(this->site);
		normaliseOp(opMatrix);
		const auto U = alfa.get_eigenvectors();
		arma::cx_mat mat_elem = U.t() * opMatrix * U;
		timeEvolution(alfa, mat_elem, timeDir + "LIOM" + std::to_string(this->site));
		for (int j = 0; j < alfa.L; j++) {
			auto opMatrix = alfa.create_LIOMoperator_densities(this->site, j);
			timeEvolution(alfa, mat_elem, timeDir + "LIOM" + std::to_string(this->site) + "," + std::to_string(j));
		}
	};
	if (this->m == 1) {
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
			this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		body(*alfa);
	}
	else {
		auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
		body(*alfa);
	}
}

template <typename _type> void isingUI::ui::LevelSpacingDist(IsingModel<_type>& alfa) {
	auto r = alfa.eigenlevel_statistics_with_return();
	auto probDist = probability_distribution_with_return(r);
	double _min = arma::min(r);
	double _max = arma::max(r);
	for (int i = 0; i < 20; i++) {
		double dg = alfa.getRandomValue(-alfa.g / 100., alfa.g / 100.);
		auto beta = std::make_shared<IsingModel<_type>>(alfa);
		beta->g = beta->g + dg;
		beta->hamiltonian();
		beta->diagonalization(true);
		auto level_spacing = beta->eigenlevel_statistics_with_return();
		probDist += probability_distribution_with_return(level_spacing);
		_min += arma::min(level_spacing);
		_max += arma::max(level_spacing);
	}
	_min /= 21.; _max /= 21.; probDist /= 21.;
	fs::create_directory(this->saving_dir + "LevelSpacing" + kPSep);
	save_to_file(this->saving_dir + "LevelSpacing" + kPSep, "LevelSpacing" + alfa->get_info({}),
		arma::linspace(_min, _max, probDist.size()), probDist);
}
//-------------------------------------------------------------------- AUTO-ENCODER
void isingUI::ui::saveDataForAutoEncoder_symmetries(std::initializer_list<op_type> operators, std::initializer_list<std::string> names) {
	using namespace std::chrono;
	clk::time_point start = std::chrono::system_clock::now();

	stout << "making symmetric model\n";
	auto params = arma::linspace(this->h, this->hs, this->hn);
	for (int system_size = this->L; system_size < this->L + this->Ls * this->Ln; system_size += this->Ls) {

		std::ofstream coeffLog;
		openFile(coeffLog, this->saving_dir + "coeffLog" + std::to_string(system_size) + ".dat", ios::out | ios::app);
		for (auto& hx : params) {
			stout << "\n\t\t\tSYM : h = " << hx << "\t\t\n";

			const auto start_loop = std::chrono::system_clock::now();
			auto alfa = std::make_unique<IsingModel_sym>(system_size, this->J, this->g, hx, \
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
			printSeparated(wavefunLog, "\t", 10, true, "filenum","sigma_x(0)");

			auto ipr = 0.0;
			auto r = 0.0;
			double entropy = 0;
			this->mu = long(0.5 * N);
			long int E_min = alfa->E_av_idx - long(mu / 2);
			long int E_max = alfa->E_av_idx + long(mu / 2);
			int counter = 0;
			for (long int i = E_min; i < E_max; i++) {
				ipr += alfa->ipr(i);
				entropy += alfa->information_entropy(i);
				NO_OVERFLOW(r += alfa->eigenlevel_statistics(i, i + 1));
				counter++;
			}
			printSeparated(coeffLog, "\t", 10, true, to_string_prec(this->g), to_string_prec(hx), to_string_prec(ipr / double(N * counter), 8), \
				to_string_prec(r / double(counter), 8), to_string_prec(entropy / double(counter), 8));


			this->mu = long(0.3 * N);
			E_min = alfa->E_av_idx - long(mu / 2);
			E_max = alfa->E_av_idx + long(mu / 2);
			// let's go over that stuff
			int w_c = 0;
			for (long int i = E_min; i < E_max; i++) {
				// check sigma_x
				std::ofstream wavefun;
				openFile(wavefun, saving_folder_wavefun + to_string_prec(hx) + "_" + std::to_string(w_c)   \
					+ "_wavefun_" + alfa->get_info() + ".dat", ios::out);
				const auto sigma_x = alfa->av_sigma_x(i, i, { 0 });
				printSeparated(wavefunLog, "\t", 10, true, to_string_prec(hx) + "_"\
					+ std::to_string(w_c) , to_string_prec(sigma_x,8));
				for (u64 j = 0; j < N; j++) {
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
void isingUI::ui::saveDataForAutoEncoder_disorder(std::initializer_list<op_type> operators, std::initializer_list<std::string> names) {
	if (operators.size() != names.size()) assert(false && "Set same number of input names as input operators!\n");
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
	printSeparated(map, "\t", 10, true, "real", "w", "g0", "c", "r", "rvar", "g", "h" );
	for (double my_g0 = this->g0; my_g0 < g0_end; my_g0 += this->g0s) {
		for (double my_w = this->w; my_w < w_end; my_w += this->ws) {
			// reset generator 
			gen = std::mt19937_64(seed);
			auto Hamil = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, my_g0, this->h, my_w, this->boundary_conditions);
			const u64 N = Hamil->get_hilbert_size();
			stout << "\n\n------------------------------ Doing : " << \
				Hamil->get_info() << "------------------------------\n";

			//create main folder
			const std::string saving_folder = this->saving_dir + Hamil->get_info() + kPSep;

			// make folders for each operator separetely
			std::vector<std::string> opDirDiag, opDirNonDiag;
			std::string saving_folder_wave = saving_folder + "wavefunctions" + kPSep;
			createDirs(saving_folder, saving_folder_wave);

			for (auto& opName : names) {
				std::string saving_folder_operator = saving_folder + opName;
				std::string saving_folder_nondiag = saving_folder_operator + kPSep + "NonDiagonalMatrixElements" + kPSep;
				std::string saving_folder_diag = saving_folder_operator + kPSep + "DiagonalMatrixElements" + kPSep;
				createDirs(saving_folder_operator, saving_folder_diag, saving_folder_nondiag);
				opDirDiag.push_back(saving_folder_nondiag);
				opDirNonDiag.push_back(saving_folder_nondiag);
			}


			for (int realis = 0; realis < realisations; realis++) {
				Hamil->hamiltonian(); // restart Hamiltonian for new values of disorder
				Hamil->diagonalization();
				stout << " \t\t--> finished diagonalizing for " << Hamil->get_info() << \
					" - in time : " << tim_s(start) << "s. Realisation -> " << realis << "\n";

				// set states from the middle of the spectrum
				this->mu = long(0.5 * N); // to include small system not not exceed Hilbert space
				long int E_min = Hamil->E_av_idx - long(mu / 2);
				long int E_max = Hamil->E_av_idx + long(mu / 2);

				// calculate level statistics 
				double r = 0;
				double r_var = 0;
				for (int k = E_min; k < E_max; k++) {
					NO_OVERFLOW(double level_stat = Hamil->eigenlevel_statistics(k, k + 1));
					r += level_stat;
					r_var += level_stat * level_stat;
				}
				NO_OVERFLOW(
					r /= -double(E_min - E_max);		// 1st moment - mean
					r_var /= -double(E_min - E_max); // 2nd moment
				);
				r_var = r_var - r * r;			// variance


				// set new from the middle of the spectrum for operators
				this->mu = (M > N) ? long(0.5 * N) : M / 2; // to include small system not not exceed Hilbert space
				E_min = Hamil->E_av_idx - long(mu / 2);
				E_max = Hamil->E_av_idx + long(mu / 2);
				arma::mat H_diag = arma::diagmat(Mat<double>(Hamil->get_hamiltonian()));
				arma::mat H_offdiag = Hamil->get_hamiltonian() - H_diag;
				auto c = matrixVariance(H_diag) / matrixVariance(H_offdiag);
				c = 1. / c;
				// iterate over input lists
				for (int a = 0; a < operators.size(); a++) {
					// assign by iterator
					op_type op = *(operators.begin() + a);
					std::string opName = *(names.begin() + a);

					// make file for log
					std::ofstream MatElemDiag;
					std::ofstream MatElemLogNonDiag;
					// save to wavefunctions log
					openFile(MatElemDiag, opDirDiag[a] + "MatrixElements_" + std::to_string(realis) + ".dat"\
						, ios::out);
					openFile(MatElemLogNonDiag, opDirNonDiag[a] + "MatrixElements_" + std::to_string(realis) + ".dat"\
						, ios::out);

					printSeparated(MatElemDiag, "\t", 10, false, "stateNum");
					printSeparated(MatElemLogNonDiag, "\t", 10, false, "<i|j>");
					for (int i = 0; i < this->L; i++) {
						printSeparated(MatElemDiag, "\t", 10, false, opName, "(", i, ")");
						printSeparated(MatElemLogNonDiag, "\t", 10, false, opName, "(", i, ")");
					}
					printSeparated(MatElemDiag, "\t", 10, true, "E_i");
					printSeparated(MatElemLogNonDiag, "\t", 10, true, "E_i - E_j");

					stout << "\n\n\t\t\t------------------------------ Starting operator " + opName + " for: " << \
						Hamil->get_info() << "------------------------------\n";
					// go through the eigenstates
					for (u64 k = E_min; k < E_max; k++) {
						std::string wavename = saving_folder_wave + \
							std::to_string(k) + "," + std::to_string(realis) + \
							"," + to_string_prec(my_w, 2) + \
							+ "," + to_string_prec(my_g0, 2) + "," + to_string_prec(c, 4);
						std::ofstream wavefunctionsLog;
						openFile(wavefunctionsLog, wavename + ".dat"\
							, ios::out);
						const int idx = long(k) - long(E_min);
						// check sigma_x
						// print k state
						printSeparated(MatElemDiag, "\t", 6, false, k);
						for (int i = 0; i < this->L; i++) {
							const auto opElem = Hamil->av_operator(k, k, op, { i });
							//const auto opElem = Hamil->av_op
							printSeparated(MatElemDiag, "\t", 10, false, to_string_prec(opElem,8));
						}
						printSeparated(MatElemDiag, "\t", 10, true, to_string_prec(Hamil->get_eigenEnergy(k), 8));

						// give nondiagonal elements
						long int k2 = long(N) - long(k);
						printSeparated(MatElemLogNonDiag, "\t", 10, false, "<", k, "|", k2, ">");
						for (int i = 0; i < this->L; i++) {
							const auto opElem = Hamil->av_operator(k, k2, op, { i }).real();
							printSeparated(MatElemLogNonDiag, "\t", 10, false, to_string_prec(opElem, 8));
						}
						printSeparated(MatElemLogNonDiag, "\t", 10, true, to_string_prec(Hamil->get_eigenEnergy(k) - Hamil->get_eigenEnergy(k2),8));
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

void isingUI::ui::LIOMsdisorder() {
	clk::time_point start = std::chrono::system_clock::now();
	auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
	stout << "\n\t\t--> finished creating model for " << alfa->get_info() << " - in time : "  << tim_s(start) << "s\n";

	//alfa->diagonalization();
	//stout << " \t\t	--> finished diagonalizing for " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
	const double tH = 1. / alfa->mean_level_spacing_analytical();
	int t_max = (int)std::ceil(std::log10(tH));
	t_max = (t_max / std::log10(tH) < 1.5) ? t_max + 1 : t_max;
	auto times = arma::logspace(-2, t_max, 300);
	auto omegas = arma::logspace(std::floor(log10(1. / tH)) - 1, 2, 300);
	std::string opName = IsingModel_disorder::opName(this->op, this->site);
	std::string str = (this->op < 3) ? "j" : "q";
	std::string timeDir = this->saving_dir + "timeEvolution" + kPSep + str + "=" + std::to_string(this->site) + kPSep;
	std::string specDir = this->saving_dir + "ResponseFunction" + kPSep + str + "=" + std::to_string(this->site) + kPSep;
	std::string intDir = this->saving_dir + "IntegratedResponseFunction" + kPSep + str + "=" + std::to_string(this->site) + kPSep;
	createDirs(timeDir, specDir, intDir);

	auto op = alfa->chooseOperator(this->op, this->site);
	//normaliseOp(op);
	// use normaliseMat as below, some operators have zero norm
	arma::cx_mat mat_elem;
	arma::vec opEvol(times.size(), arma::fill::zeros);
	arma::vec opIntSpec(times.size(), arma::fill::zeros);

	double LTA = 0;
	if (this->realisations > 1) {
		alfa->reset_random();
		for (int r = 0; r < this->realisations; r++) {
			const auto start_loop = std::chrono::system_clock::now();
			std::string tdir_realisation = timeDir + "realisation=" + std::to_string(r) + kPSep;
			std::string intdir_realisation = intDir + "realisation=" + std::to_string(r) + kPSep;
			std::string specdir_realisation = specDir + "realisation=" + std::to_string(r) + kPSep;
			createDirs(tdir_realisation, intdir_realisation, specdir_realisation);
			alfa->diagonalization();
			const arma::mat U = alfa->get_eigenvectors();
			mat_elem = U.t() * op * U;
			normaliseMat(mat_elem);
			stout << "\t\t	--> finished diagonalizing for " << alfa->get_info() 
				<< " realisation: " << r << " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s\n";
			
			auto [op_tmp, LTA_tmp] = timeEvolution(*alfa, mat_elem, times);
			save_to_file(tdir_realisation + opName + alfa->get_info({}) + ".dat", times, op_tmp, tH, LTA_tmp);
			stout << "\t\t	--> finished time evolution for " << alfa->get_info() 
				<< " realisation: " << r << " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s\n\tNEXT: integrated spectral function\n";
			
			auto res = integratedSpectralFunction(*alfa, mat_elem, omegas);
			save_to_file(intdir_realisation + opName + alfa->get_info({}) + ".dat", omegas, res, 1. / tH, LTA_tmp);
			stout << "\t\t	--> finished integrated spectral function for " << alfa->get_info()
				<< " realisation: " << r << " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s\n\tNEXT: spectral function\n";
			
			spectralFunction(*alfa, mat_elem, specdir_realisation + opName);
			
			LTA += LTA_tmp;
			opEvol += op_tmp;
			opIntSpec += res;
			alfa->hamiltonian();
		}
		opEvol /= double(this->realisations);
		LTA /= double(this->realisations);
		opIntSpec /= double(this->realisations);
		save_to_file(timeDir + opName + alfa->get_info({}) + ".dat", times, opEvol, tH, LTA);
		save_to_file(intDir + opName + alfa->get_info({}) + ".dat", omegas, opIntSpec, 1. / tH, LTA);
	}
	else {
		alfa->diagonalization();
		const arma::mat U = alfa->get_eigenvectors();
		mat_elem = U.t() * op * U;
		normaliseMat(mat_elem);
		stout << "\t\t	--> finished diagonalizing for " << alfa->get_info() << " - in time : " << tim_s(start) << "s\n";
		auto [opEvol, LTA] = timeEvolution(*alfa, mat_elem, times);
		save_to_file(timeDir + opName + alfa->get_info({}) + ".dat", times, opEvol, tH, LTA);
		stout << "\t\t	--> finished time evolution for " << alfa->get_info() << " - in time : " << tim_s(start) << "s\n\tNEXT: integrated spectral function\n";
		auto res = integratedSpectralFunction(*alfa, mat_elem, omegas);
		save_to_file(intDir + opName + alfa->get_info({}) + ".dat", omegas, res, 1. / tH, LTA);
		stout << "\t\t	--> finished integrated spectral function for " << alfa->get_info() << " - in time : " << tim_s(start) << "s\n\tNEXT: spectral function\n";
		spectralFunction(*alfa, mat_elem, specDir + opName);
	}
	stout << " - - - - - - FINISHED CALCULATIONS IN : " << tim_s(start) << " seconds - - - - - - \n" << endl;						// simulation end
}



//----------------------------------------------------------------------------------------------------------------UI main
void isingUI::ui::make_sim() {
	printAllOptions();
	seed = static_cast<long unsigned int>(time(0));
	clk::time_point start = std::chrono::system_clock::now();
	//compare_energies();
	//disorder();
	//adiabaticGaugePotential(0, 0);
	//relaxationTimesFromFiles();
	switch (this->fun) {
		case 0: LIOMsdisorder();				break;
		case 1: adiabaticGaugePotential_dis(1);	break;
		case 2: adiabaticGaugePotential_dis(0);	break;
		case 3: TFIsingLIOMs();					break;
		case 4: 
			for (this->op = 2; this->op < 5; this->op += 2)
				for (this->L = 10; this->L <= 15; this->L++)
					for (this->site = 0; this->site <= this->L / 2; this->site++)
						relaxationTimesFromFiles();
			break;
	default:
		const int Lmin = this->L, Lmax = this->L + this->Ln * this->Ls;
		const double gmin = this->g, gmax = this->g + this->gn * this->gs;
		const double hmin = this->h, hmax = this->h + this->hn * this->hs;
		for (int system_size = Lmin; system_size < Lmax; system_size += this->Ls) {
			std::ofstream map;
			//openFile(map, this->saving_dir + "FullMap" + IsingModel_sym::set_info(system_size, J, g, h, \
			//	this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, { "h" }) + ".dat", ios::out);
			for (double gx = gmin; gx < gmax; gx += this->gs) {
				for (double hx = hmin; hx < hmax; hx += this->hs) {
					const auto start_loop = std::chrono::system_clock::now();
					stout << "h = " << hx << "\t\t";
					this->L = system_size;
					this->g = gx;
					this->h = hx;
					

					auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
					stout << "\n\t\t--> finished creating model for " << alfa->get_info() << " - in time : " << tim_s(start) << "s\n";
					const long N = alfa->get_hilbert_size();
					const double tH = 1. / alfa->mean_level_spacing_analytical();
					int t_max = (int)std::ceil(std::log10(tH));
					t_max = (t_max / std::log10(tH) < 1.5) ? t_max + 1 : t_max;
					auto times = arma::logspace(-2, t_max, 100);
					std::string dir = this->saving_dir + "Entropy" + kPSep;
					createDirs(dir);
					
					double entropy_1 = 0, entropy_2 = 0, entropy_3 = 0;
					this->mu = N < 1000 ? 0.2 * N : 500;
					std::function to_ave_subsystem = [&]() {
						//alfa->diagonalization();
						stout << "\t\t	--> finished diagonalizing for " << alfa->get_info()
							<< " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s\n";
						const u64 E_min = alfa->E_av_idx - this->mu / 2.;
						const u64 E_max = alfa->E_av_idx + this->mu / 2.;
						double entropy_1_tmp = 0, entropy_2_tmp = 0, entropy_3_tmp = 0;
#pragma omp parallel for reduction(+: entropy_1_tmp, entropy_2_tmp, entropy_3_tmp)
						for (long k = E_min; k < E_max; k++) {
							arma::cx_vec state = cx_vec(alfa->get_eigenState(k), arma::vec(N, arma::fill::zeros));
							entropy_1_tmp += alfa->entaglement_entropy(state, this->L / 2 - 1);
							entropy_2_tmp += alfa->entaglement_entropy(state, this->L / 2);
							entropy_3_tmp += alfa->entaglement_entropy(state, this->L / 2 + 1);
						}
						entropy_1 += entropy_1_tmp / double(this->mu);
						entropy_2 += entropy_2_tmp / double(this->mu);
						entropy_3 += entropy_3_tmp / double(this->mu);
					};
					average_over_realisations(*alfa, true, to_ave_subsystem);
					entropy_1 /= double(this->realisations);
					entropy_2 /= double(this->realisations);
					entropy_3 /= double(this->realisations);
					std::vector S = { entropy_1, entropy_2, entropy_3 };
					std::ofstream file;
					openFile(file, dir + "SubsystemSize" + alfa->get_info({}) + ".dat");
					const double psi_N = digamma(N + 1.0);
					int j0 = this->L / 2 - 1;
					for (int j = j0; j <= j0 + 2; j++) {
						const long d_A = (u64)std::pow(2, j);
						const long d_B = (u64)std::pow(2, this->L - j);
						const double psi_A = digamma(d_A + 1.0);
						const double psi_B = digamma(d_B + 1.0);
						const double S_A = j > this->L / 2 ?
							psi_N - psi_A - (d_B - 1.0) / (2.0 * d_A) :
							psi_N - psi_B - (d_A - 1.0) / (2.0 * d_B);
						printSeparated(file, "\t", 12, true, j / double(this->L), S[j - j0], S_A);
					}
					file.close();

					//arma::vec entropy(times.size(), arma::fill::zeros);					
					//arma::vec down = { 0,1 };
					//arma::vec up = { 1,0 };
					//std::uniform_real_distribution<double> theta(0, pi);
					//std::uniform_real_distribution<double> fi(0, pi);
					//alfa->reset_random();
					//std::function to_ave_time = [&]() {
					//	//alfa->diagonalization();
					//	stout << "\t\t	--> finished diagonalizing for " << alfa->get_info()
					//		<< " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s\n";
					//	auto the = theta(gen);
					//	arma::cx_vec init_state = std::cos(the / 2.) * up + std::exp(im * fi(gen)) * std::sin(the / 2.) * down;
					//	for (int j = 1; j < this->L; j++) {
					//		the = theta(gen);
					//		init_state = arma::kron(init_state, std::cos(the / 2.) * up + std::exp(im * fi(gen)) * std::sin(the / 2.) * down);
					//	}
					//	init_state = arma::normalise(init_state);
					//	for (int i = 0; i < times.size(); i++) {
					//		auto t = times(i);
					//		arma::cx_vec state = arma::normalise(init_state);
					//		alfa->time_evolve_state(state, t);
					//		entropy(i) += alfa->entaglement_entropy(state, this->L / 2);
					//	}
					//
					//};
					//
					//average_over_realisations<>(*alfa, false, to_ave_time);
					//entropy /= double(this->realisations);
					//std::ofstream file;
					//openFile(file, dir + "TimeEvolution" + alfa->get_info({}) + ".dat");
					//for (int j = 0; j < times.size(); j++) 
					//	printSeparated(file, "\t", 12, true, times(j), entropy(j));
					//file.close();


					//this->mu = 0.5 * N;
					//long int E_min = alfa->E_av_idx - mu / 2.;
					//long int E_max = alfa->E_av_idx + mu / 2.;
					//
					//double r = 0, ipr = 0, S = 0;
					//for (long int k = E_min; k < E_max; k++) {
					//	r += alfa->eigenlevel_statistics(k, k + 1);
					//	ipr += alfa->ipr(k);
					//	S += alfa->information_entropy(k);
					//}
					//r /= double(E_max - E_min);
					//ipr /= double(E_max - E_min) * N;
					//S /= double(E_max - E_min);
					//int nums = 10;
					//for (int i = 0; i < nums; i++) {
					//	double dg = alfa->getRandomValue(-gx / 100., gx / 100.);
					//	auto beta = std::make_unique<IsingModel_sym>(system_size, this->J, gx + dg, hx, \
					//		this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
					//	beta->diagonalization();
					//	r += beta->eigenlevel_statistics(E_min, E_max);
					//}
					//r /= double(nums + 1);
					//
					//printSeparated(map, "\t", { gx, hx, r,ipr,S }, 12, true);

					//probability_distribution(this->saving_dir + "LevelSpacing" + kPSep, "LevelSpacing" + alfa->get_info({}), level_spacing);

					//IsingLIOMs(*alfa);
					//TFIsingLIOMs(*alfa);
					//auto opMatrix = alfa->create_operator({ IsingModel_sym::sigma_x });
					//spectralFunction(*alfa, opMatrix, "SigmaX");
					//opMatrix = alfa->create_operator({ IsingModel_sym::sigma_z });
					//spectralFunction(*alfa, opMatrix, "SigmaZ");

				}
			}
			map.close();
		}
	}
	

	

	stout << " - - - - - - FINISHED CALCULATIONS IN : " << tim_s(start) << " seconds - - - - - - " << endl;						// simulation end
}
