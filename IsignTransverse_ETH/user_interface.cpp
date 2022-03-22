#include "include/user_interface.h"
// set externs
std::uniform_real_distribution<> theta	= std::uniform_real_distribution<>(0.0, pi);
std::uniform_real_distribution<> fi		= std::uniform_real_distribution<>(0.0, pi);
//---------------------------------------------------------------------------------------------------------------- UI main
void isingUI::ui::make_sim()
{
	printAllOptions();
	seed = static_cast<long unsigned int>(time(0));
	clk::time_point start = std::chrono::system_clock::now();
	// compare_energies();
	// disorder();
	// adiabaticGaugePotential(0, 0);
	// relaxationTimesFromFiles();
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
		for (this->L = Lmin; this->L < Lmax; this->L += this->Ls)
			for (this->site = 0; this->site <= this->L / 2; this->site++)
				relaxationTimesFromFiles();
		break;
	case 5:
		benchmark();
		break;
	case 6:
		adiabaticGaugePotential_dis();
		break;
	default:
		for (int system_size = Lmin; system_size < Lmax; system_size += this->Ls)
		{
			for (double gx = gmin; gx < gmax; gx += this->gs)
			{
				for (double hx = hmin; hx < hmax; hx += this->hs)
				{
					const auto start_loop = std::chrono::system_clock::now();
					stout << "h = " << hx << "\t\t";
					this->L = system_size;
					this->g = gx;
					this->h = hx;

					// compare_entaglement();
					// continue;
					auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
					stout << "\n\t\t--> finished creating model for " << alfa->get_info() << " - in time : " << tim_s(start) << "s" << std::endl;
					// const long long N = alfa->get_hilbert_size();
					stout << "\t\t	--> start diagonalizing for " << alfa->get_info()
						  << " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s" << std::endl;
					alfa->diagonalization();

					stout << "\t\t	--> finished diagonalizing for " << alfa->get_info()
						  << " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s" << std::endl;
					const double tH = 1. / alfa->mean_level_spacing_analytical();
					int t_max = (int)std::ceil(std::log10(tH));
					t_max = (t_max / std::log10(tH) < 1.5) ? t_max + 1 : t_max;
					// auto times = arma::logspace(-2, t_max, 300);
					std::string dir = this->saving_dir + "Entropy" + kPSep;
					createDirs(dir);
					double omega_max = alfa->get_eigenEnergy(alfa->get_hilbert_size() - 1) - alfa->get_eigenEnergy(0);

					//for (double x = 1e-2; x <= 10; x *= 2)
					double x = this->dt;
					{
						double dt_new = x;//M * 10 * x / omega_max;
						if (this->scale) stout << "WARNING: log only valid for t<10, for larger t linear is resumed with dt = " << dt_new << std::endl;
						auto init_log = arma::logspace(-2, t_max, 4000);
						auto rest_lin = arma::regspace(10.0, dt_new, 10 * tH);
						auto times = this->scale ? arma::join_cols(exctract_vector(init_log, 0.0, 10.0), rest_lin) : arma::regspace(dt_new, dt_new, tH);
						// temporary check
						times = this->scale ? arma::logspace(-3, t_max, 500) : arma::regspace(dt_new, dt_new, tH);
					//	for (int M = 2; M < 10; M++)
					int M = this->mu;
						{
							alfa->reset_random();
							stout << "\t\t	-->set random generators for " << alfa->get_info()
								  << " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s" << std::endl;
							if (M >= alfa->get_hilbert_size())
								continue;
							arma::vec entropy(times.size(), arma::fill::zeros);
							arma::vec entropy_lanczos(times.size(), arma::fill::zeros);

							lanczosParams params(1000, 1, true, false);
							lanczos::Lanczos lancz(alfa->get_hamiltonian(), std::move(params));
							//lancz.diagonalization();
							
							std::function to_ave_time = [&]()
							{
								const arma::cx_vec init_state = random_product_state(this->L);
								arma::cx_vec state2 = init_state;
								lancz.diagonalization(init_state);
								for (int i = 0; i < times.size(); i++)
								{
									auto t = times(i);
									arma::cx_vec state = arma::normalise(init_state);
									state2 = lancz.time_evolution_stationary(state, t);
									alfa->time_evolve_state(state, t);
									//lancz.time_evolution_non_stationary(state2, t - (i == 0 ? 0.0 : times(i - 1)), M);
									//auto state2 =  lancz.time_evolution_stationary(init_state, t);
									entropy(i) += alfa->entaglement_entropy(state, this->L / 2);
									entropy_lanczos(i) += alfa->entaglement_entropy(state2, this->L / 2);
								}
							};
							average_over_realisations(*alfa, false, to_ave_time);
							entropy /= double(this->realisations);
							entropy_lanczos /= double(this->realisations);
							std::ofstream file;
							openFile(file, dir + "TimeEvolution" + alfa->get_info({}) + ".dat");
							for (int j = 0; j < times.size(); j++)
							{
								double diff = entropy(j) - entropy_lanczos(j);
								printSeparated(file, "\t", 16, true, times(j), entropy(j), entropy_lanczos(j), diff);
								if(abs(diff) > 1e-8)
									printSeparated(std::cout, "\t", 16, true, times(j), entropy(j), entropy_lanczos(j), diff);
							}
							file.close();
						}
					}

					{
						// arma::vec entropy_1(N, arma::fill::zeros);
						// arma::vec entropy_2 = entropy_1, entropy_3 = entropy_1;
						// this->mu = N < 1000 ? 0.2 * N : 500;
						// std::function to_ave_subsystem = [&]() {
						//	//alfa->diagonalization();
						//	stout << "\t\t	--> finished diagonalizing for " << alfa->get_info()
						//		<< " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s\n";
						//	const u64 E_min = 0;// alfa->E_av_idx - this->mu / 2.;
						//	const u64 E_max = N;// alfa->E_av_idx + this->mu / 2.;
						//	double entropy_1_tmp = 0, entropy_2_tmp = 0, entropy_3_tmp = 0;
						//#pragma omp parallel for reduction(+: entropy_1_tmp, entropy_2_tmp, entropy_3_tmp)
						//	for (long k = E_min; k < E_max; k++) {
						//		arma::cx_vec state = cx_vec(alfa->get_eigenState(k), arma::vec(N, arma::fill::zeros));
						//		entropy_1(k) = alfa->entaglement_entropy(state, this->L / 2 - 1);
						//		entropy_2(k) = alfa->entaglement_entropy(state, this->L / 2);
						//		entropy_3(k) = alfa->entaglement_entropy(state, this->L / 2 + 1);
						//	}
						//	//entropy_1 += entropy_1_tmp / double(this->mu);
						//	//entropy_2 += entropy_2_tmp / double(this->mu);
						//	//entropy_3 += entropy_3_tmp / double(this->mu);
						// };
						//
						//
						// average_over_realisations(*alfa, true, to_ave_subsystem);
						// entropy_1 /= double(this->realisations);
						// entropy_2 /= double(this->realisations);
						// entropy_3 /= double(this->realisations);
						//
						// std::ofstream file;
						// openFile(file, dir + "Eigenstates" + alfa->get_info({}) + ".dat");
						// for (long i = 0; i < N; i++) {
						//	printSeparated(file, "\t", 12, true, alfa->get_eigenEnergy(i) / double(this->L), entropy_1(i), entropy_2(i), entropy_3(i));
						// }
						// std::vector S = { entropy_1, entropy_2, entropy_3 };
						// std::ofstream file;
						// openFile(file, dir + "SubsystemSize" + alfa->get_info({}) + ".dat");
						// const double psi_N = digamma(N + 1.0);
						// int j0 = this->L / 2 - 1;
						// for (int j = j0; j <= j0 + 2; j++) {
						//	const long d_A = (u64)std::pow(2, j);
						//	const long d_B = (u64)std::pow(2, this->L - j);
						//	const double psi_A = digamma(d_A + 1.0);
						//	const double psi_B = digamma(d_B + 1.0);
						//	const double S_A = j > this->L / 2 ?
						//		psi_N - psi_A - (d_B - 1.0) / (2.0 * d_A) :
						//		psi_N - psi_B - (d_A - 1.0) / (2.0 * d_B);
						//	printSeparated(file, "\t", 12, true, j / double(this->L), S[j - j0], S_A);
						// }
						// file.close();
					}
				}
			}
		}
	}

	stout << " - - - - - - FINISHED CALCULATIONS IN : " << tim_s(start) << " seconds - - - - - - " << std::endl; // simulation end
}

//---------------------------------------------------------------------------------------------------------------- IMPLEMENTATION OF UI
//---------------------------------------------------------------------------------------------------------------- FUNCTIONS AND MORE
//---------------------------------------------------------------------------------------------------------------- 

void print_help(){
	printf(
		"\n Usage: name of executable [options] outputDirectory \n"
		"\n The input can be both introduced with [options] described below or with giving the input directory(which also is the flag in the options)\n"
		"\n options:\n"
		"\n-f input file for all of the options : (default none)\n"
		"\n-mu bucket size for ergodic coefficients (default 5)\n"
		"\n-J spin exchange coefficient : (default 1)\n"
		"\n-J0 random spin exchange set in uniform distribution [-J0,J0]\n"
		"\n-g transverse magnetic field (x-) constant: (default 1)\n"
		"\n-gs transverse magnetic field (x-) constant step: (default 0.0)\n"
		"\n-gn transverse magnetic field (x-) constant number: (default 1)\n"
		"\n-g0 random transverse field set in uniform distribution [-g0,g0]\n"
		"\n-g0s transverse disorder strength step: (default 0.0)\n"
		"\n-g0n transverse disorder strength number: (default 1)\n"
		"\n-h perpendicular (z-) magnetic field constant: (default 0)\n"
		"\n-hs perpendicular (z-) magnetic field constant step: (default 0.0)\n"
		"\n-hn perpendicular (z-) magnetic field constant number: (default 1)\n"
		"\n-w disorder strength : (default 0 - no disorder introduced)\n"
		"\n-ws disorder strength step: (default 0.0)\n"
		"\n-wn disorder strength number: (default 1)\n"
		"\n-L chain length minimum: bigger than 0 (default 8)\n"
		"\n-Ls chain length step: bigger equal than 0 (default 0)\n"
		"\n-Ln chain length number: bigger than 0 (default 1)\n"
		"\n-b boundary conditions : bigger than 0 (default 0 - PBC)\n"
		"\n	0 -- PBC\n"
		"\n	1 -- OBC\n"
		"\n	2 -- ABC -- none so far implemented\n"
		"\n-s site to act with local operators (default 0)"
		"\n-op flag to choose operator: "
		"\n	0 -- Sz_i-local"
		"\n	1 -- Sx_i-local"
		"\n	2 -- Hi"
		"\n	3 -- Sz_q"
		"\n	4 -- Sx_q"
		"\n	5 -- Hq"
		"\n	  -- to get sum of local Sz or Sx take Sz_q or Sx_q with -s=0"
		"\n	  -- i or q are set to site (flag -s); (default 0)"
		"\n"
		"\n-fun choose function to start calculations: check user_interface.cpp -> make_sim() to find functions"
		"\n\t	0 -- diagonalizing hamiltonian and writing to file eigenvalues. Set -ch=1 to include eigenvector calculation"
		"\n\t	1 -- time evolution (and spectral functions) for any model (disordered is with averaging):\n\t\t set -op for operator and -s for acting site"
		"\n\t	2 -- evolution of entropy from initial state chosen by the -op flag:"
		"\n\t		for both models (-m flag) and possible to use lanczos iterative method by setting -ch=1"
		"\n\t		use -mu to set number of lanczos steps (<10 is enough) and -ts as time step (divided by E_max - E_min): 0.1 is sufficient"
		"\n\t			* op=0 -- random initial product state averaged over -r realisations"
		"\n\t			* op=1 -- fully ferromagnetically polarised state |111111...>"
		"\n\t			* op=2 -- fully anti-ferromagnetically polarised state |111111...>"
		"\n\t	3 -- spectral form factor calculation (checks if file exists with data, if not then diagonalize and save"
		"\n\t	4 -- relaxation times from integrated spectral function for:\n\t\t operator -op flag on site -s flag\
		\t\t	(also derivative of integrated spectral function is calculated)\n\t\tlooped over system sizes: -L, -Ls, -Ln and sites: from 0 to L/2"
		"\n\t   5 -- benchmark diagonalization routines vs CPU count:\n\t\tlooped over different system sizes set by -L, -Ln, -Ls\n\t\tfor number of threads: 1, 2, 4, 8, 16, 24, 32, 40, 48, 64"
		"\n\t	6 -- AGPs for small disorder (-m=0) as function of h for -ch=1 or as function of g for -ch=0 for input operator from -op flag"
		"\n\t		SET: -L, -Ln, -Ls, -h, -hn, -hs, -op, -w(default=0.01)"
		"\n def -- in make_sim space for user to write function; designed for non-builtin behavior"
		"\n"
		"\n-m model to be choosen : (default 0 - without symmetries)\n"
		"\n	0 -- nonsymmetric model - only here the disorder is working\n"
		"\n	1 -- include symmetries - here the parity flag is also working\n"
		"\n-k translation symetry sector, 0-L, (default 0)\n"
		"\n-p parity symmetry sector, +-1 (if applicable) (default 1)\n"
		"\n-x spin flip symmetry sector, +-1 (if applicable) (default 1)\n"
		"\n-th number of threads to be used for CPU parallelization : depends on the machine specifics, default(1)"
		"\n-ch general boolean flag used in different context (default: 0)"
		"\n-ts time step for evolution (default: 0.1)"
		"\n-scale choose scale for data: either linear-0 or log-1 (default: linear)"
		"\n-h quit with help\n");
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
	print_help();
	stout << "---------------------------------------------------------------------------------\n\n";
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
	std::string str_model = std::string(kPathSeparator) + "disorder" + std::string(kPathSeparator); // folder for current model
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
	choosen_option = "-ts";
	this->set_option(this->dt, argv, choosen_option);
	choosen_option = "-ch";
	this->set_option(this->ch, argv, choosen_option);
	choosen_option = "-scale";
	this->set_option(this->scale, argv, choosen_option);

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
	if (this->thread_number > std::thread::hardware_concurrency())
		this->set_default_msg(this->thread_number, choosen_option.substr(1),
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
			if (fs::create_directories(folder) || fs::is_directory(folder)) // creating the directory for saving the files with results
				this->saving_dir = folder;									// if can create dir this is is
		} else {
			if (fs::create_directories(folder) || fs::is_directory(folder)) // creating the directory for saving the files with results
				this->saving_dir = folder;									// if can create dir this is is
		})

	std::cout << " - - - - - - MAKING ISING INTERFACE AND USING OUTER THREADS : "
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

// ----------------------------------------------------------------------------- SIMULATIONS -----------------------------------------------------------------------------
//-------------------------------------------------------------------------- GENERAL ROUTINES
void isingUI::ui::diagonalize(){
	auto kernel = [this](auto& alfa){
		std::string info = alfa.get_info({});
		alfa.diagonalization(this->ch);
		auto eigenvalues = alfa.get_eigenvalues();

		std::string dir = this->saving_dir + "EIGENVALUES" + kPSep;
		createDirs(dir);
		std::string name = dir + kPSep + info;
		eigenvalues.save(arma::hdf5_name(name + ".h5", "eigenvalues"));
		if(this->ch){
			dir = this->saving_dir + "EIGENVECTORS" + kPSep;
			name = dir + kPSep + info;
			createDirs(dir);
			auto V = alfa.get_eigenvectors();
			V.save(arma::hdf5_name(name + ".h5", "eigenvectors"));
		}
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
			// cpx A = av_operator(element.first, t.first, *alfa, *alfa, op);
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
			// cpx A = av_operator(element.first, t.first, *beta, *beta, op);
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
			// cpx A = av_operator(element.first, t.first, *alfa, *beta, op);
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

	auto beta1 = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h, 0, 1, 1, this->boundary_conditions);
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
			auto state = arma::cx_vec(alfa1->get_eigenState(k), arma::vec(dim, arma::fill::zeros));
			entropy_dis1 += alfa1->entaglement_entropy(state, i);
			state = arma::cx_vec(alfa2->get_eigenState(k), arma::vec(dim, arma::fill::zeros));
			entropy_dis2 += alfa2->entaglement_entropy(state, i);
			state = arma::cx_vec(alfa3->get_eigenState(k), arma::vec(dim, arma::fill::zeros));
			entropy_dis3 += alfa3->entaglement_entropy(state, i);
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
			auto state = beta1->get_eigenState(k);
			entropy_sym1 += beta1->entaglement_entropy(state, i);
		}
		entropy_sym1 /= double(this->mu);

		this->mu = beta2->get_hilbert_size() > 3000 ? 500 : 0.25 * beta2->get_hilbert_size();
		E_min = beta2->E_av_idx - this->mu / 2.;
		E_max = beta2->E_av_idx + this->mu / 2.;
		double entropy_sym2 = 0.0;
		for (long k = E_min; k < E_max; k++)
		{
			auto state = beta2->get_eigenState(k);
			entropy_sym2 += beta2->entaglement_entropy(state, i);
		}
		entropy_sym2 /= double(this->mu);

		printSeparated(file, "\t", 12, true, i, entropy_dis1, entropy_dis2, entropy_dis3, entropy_sym1, entropy_sym2);
		printSeparated(std::cout, "\t", 12, true, i, entropy_dis1, entropy_dis2, entropy_dis3, entropy_sym1, entropy_sym2);
	}

	std::cout << std::endl;
}
void isingUI::ui::benchmark(bool full)
{
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
		printSeparated(file, "\t", 16, true, "#cores", "chain length", "dim", "with eigenvec 'dc'", "with eigenvec 'std'", "only eigenvalues", "in seconds");
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
					double tim1 = tim_s(start);
					// start = std::chrono::system_clock::now();
					// alfa->diagonalization(false, "std"); double tim2 = tim_s(start);
					start = std::chrono::system_clock::now();
					alfa->diagonalization(false);
					double tim3 = tim_s(start);
					printSeparated(file, "\t", 16, true, th, system_size, alfa->get_hilbert_size(), tim1, "-----------", tim3);
					printSeparated(std::cout, "\t", 16, true, th, system_size, alfa->get_hilbert_size(), tim1, "-----------", tim3);
				}
				else
				{
					auto alfa = std::make_unique<IsingModel_disorder>(system_size, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
					auto start = std::chrono::system_clock::now();
					alfa->diagonalization(true, "dc");
					double tim1 = tim_s(start);
					// start = std::chrono::system_clock::now();
					// alfa->diagonalization(false, "std"); double tim2 = tim_s(start);
					start = std::chrono::system_clock::now();
					alfa->diagonalization(false);
					double tim3 = tim_s(start);
					printSeparated(file, "\t", 16, true, th, system_size, alfa->get_hilbert_size(), tim1, "-----------", tim3);
					printSeparated(std::cout, "\t", 16, true, th, system_size, alfa->get_hilbert_size(), tim1, "-----------", tim3);
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
			double tim1 = tim_s(start);
			start = std::chrono::system_clock::now();
			alfa->diagonalization(true, "std");
			double tim2 = tim_s(start);
			start = std::chrono::system_clock::now();
			alfa->diagonalization(false);
			double tim3 = tim_s(start);
			printSeparated(std::cout, "\t", 16, true, this->thread_number, this->L, alfa->get_hilbert_size(), tim1, tim2, tim3);
			if (this->thread_number == 32)
				std::cout << std::endl;
		}
		else
		{
			auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
			auto start = std::chrono::system_clock::now();
			alfa->diagonalization(true, "dc");
			double tim1 = tim_s(start);
			start = std::chrono::system_clock::now();
			alfa->diagonalization(true, "std");
			double tim2 = tim_s(start);
			start = std::chrono::system_clock::now();
			alfa->diagonalization(false);
			double tim3 = tim_s(start);
			printSeparated(std::cout, "\t", 16, true, this->thread_number, this->L, alfa->get_hilbert_size(), tim1, tim2, tim3);
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
	auto kernel = [this](auto& alfa){
		std::string str = (this->op < 3) ? "j" : "q";
		if(this->op == 6) str = "n";
		std::string timeDir = this->saving_dir + "timeEvolution" + kPSep + str + "=" + std::to_string(this->site) + kPSep;
		std::string specDir = this->saving_dir + "ResponseFunction" + kPSep + str + "=" + std::to_string(this->site) + kPSep;
		std::string intDir = this->saving_dir + "IntegratedResponseFunction" + kPSep + str + "=" + std::to_string(this->site) + kPSep;
		std::string ssfDir = this->saving_dir + "SpectralFormFactor" + kPSep;
		createDirs(timeDir, specDir, intDir, ssfDir);

		const double tH = 1. / alfa.mean_level_spacing_analytical();
		int t_max = (int)std::ceil(std::log10(tH));
		t_max = (t_max / std::log10(tH) < 1.5) ? t_max + 1 : t_max;
		auto times = arma::logspace(-2, t_max, 300);
		auto omegas = arma::logspace(std::floor(log10(1. / tH)) - 1, 2, 300);

		clk::time_point start = std::chrono::system_clock::now();
		auto op = alfa.chooseOperator(this->op, this->site);
		std::string opName = IsingModel_disorder::opName(this->op, this->site);
		stout << "\n\t\t--> finished generating operator for " << alfa.get_info() << " - in time : " << tim_s(start) << "s\n";
		// normaliseOp(op);
		//  use normaliseMat as below, some operators have zero norm

		if (this->realisations > 1 && !this->m)
		{
			arma::vec opEvol(times.size(), arma::fill::zeros);
			arma::vec opIntSpec(times.size(), arma::fill::zeros);
			arma::vec ssf(times.size(), arma::fill::zeros);
			double LTA = 0;
			alfa.reset_random();
			for (int r = 0; r < this->realisations; r++)
			{
				const auto start_loop = std::chrono::system_clock::now();
				std::string tdir_realisation = timeDir + "realisation=" + std::to_string(r) + kPSep;
				std::string intdir_realisation = intDir + "realisation=" + std::to_string(r) + kPSep;
				std::string specdir_realisation = specDir + "realisation=" + std::to_string(r) + kPSep;
				std::string ssfdir_realisation = ssfDir + "realisation=" + std::to_string(r) + kPSep;
				createDirs(tdir_realisation, intdir_realisation, specdir_realisation, ssfdir_realisation);

				stout << "\t\t	--> start diagonalizing for " << alfa.get_info() << " - in time : " << tim_s(start_loop) << " s" << std::endl;
				alfa.diagonalization();
				stout << "\t\t	--> finished diagonalizing for " << alfa.get_info() << " - in time : " << tim_s(start_loop) << "s" << std::endl;
				auto U = alfa.get_eigenvectors();
				stout << "\t\t	--> got eigenvectors for " << alfa.get_info() << " - in time : " << tim_s(start_loop) << "s" << std::endl;
				arma::cx_mat mat_elem = U.t() * op * U;
				stout << "\t\t	--> set matrix elements for " << alfa.get_info() << " - in time : " << tim_s(start_loop) << "s" << std::endl;
				normaliseMat(mat_elem);

				auto [op_tmp, LTA_tmp] = spectrals::timeEvolution(alfa, mat_elem, times);
				save_to_file(tdir_realisation + opName + alfa.get_info({}) + ".dat", times, op_tmp, tH, LTA_tmp);
				stout << "\t\t	--> finished time evolution for " << alfa.get_info()
					  << " realisation: " << r << " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s\n\tNEXT: integrated spectral function" << std::endl;

				auto res = spectrals::integratedSpectralFunction(alfa, mat_elem, omegas);
				save_to_file(intdir_realisation + opName + alfa.get_info({}) + ".dat", omegas, res, 1. / tH, LTA_tmp);
				stout << "\t\t	--> finished integrated spectral function for " << alfa.get_info()
					  << " realisation: " << r << " - in time : " << tim_s(start_loop) << "\t\nTotal time : " << tim_s(start) << "s\n\tNEXT: spectral function" << std::endl;

				auto ssf_temp = spectrals::spectral_form_factor(alfa.get_eigenvalues(), times);
				save_to_file(ssfdir_realisation + alfa.get_info({}) + ".dat", times, ssf_temp, tH, LTA_tmp);
				spectrals::spectralFunction(alfa, mat_elem, specdir_realisation + opName);

				ssf += ssf_temp;
				LTA += LTA_tmp;
				opEvol += op_tmp;
				opIntSpec += res;
				alfa.hamiltonian();
			}
			opEvol /= double(this->realisations);
			LTA /= double(this->realisations);
			opIntSpec /= double(this->realisations);
			ssf /= double(this->realisations);
			save_to_file(timeDir + opName + alfa.get_info({}) + ".dat", times, opEvol, tH, LTA);
			save_to_file(intDir + opName + alfa.get_info({}) + ".dat", omegas, opIntSpec, 1. / tH, LTA);
			save_to_file(ssfDir + alfa.get_info({}) + ".dat", times, ssf, tH);
		}
		else
		{
			alfa.diagonalization();
			stout << "\t\t	--> finished diagonalizing for " << alfa.get_info() << " - in time : " << tim_s(start) << "s" << std::endl;
			auto U = alfa.get_eigenvectors();
			stout << "\t\t	--> got eigenvectors for " << alfa.get_info() << " - in time : " << tim_s(start) << "s" << std::endl;
			arma::cx_mat mat_elem = U.t() * op * U;
			stout << "\t\t	--> set matrix elements for " << alfa.get_info() << " - in time : " << tim_s(start) << "s" << std::endl;
			normaliseMat(mat_elem);
			auto [opEvol, LTA] = spectrals::timeEvolution(alfa, mat_elem, times);
			save_to_file(timeDir + opName + alfa.get_info({}) + ".dat", times, opEvol, tH, LTA);
			stout << "\t\t	--> finished time evolution for " << alfa.get_info() << " - in time : " << tim_s(start) << "s\n\tNEXT: integrated spectral function" << std::endl;
			auto res = spectrals::integratedSpectralFunction(alfa, mat_elem, omegas);
			save_to_file(intDir + opName + alfa.get_info({}) + ".dat", omegas, res, 1. / tH, LTA);
			stout << "\t\t	--> finished integrated spectral function for " << alfa.get_info() << " - in time : " << tim_s(start) << "s\n\tNEXT: spectral function" << std::endl;
			spectrals::spectralFunction(alfa, mat_elem, specDir + opName);
			auto ssf = spectrals::spectral_form_factor(alfa.get_eigenvalues(), times);
			save_to_file(ssfDir + alfa.get_info({}) + ".dat", times, ssf, tH);
		}
		stout << " - - - - - - FINISHED CALCULATIONS IN : " << tim_s(start) << " seconds - - - - - - \n"
			  << std::endl; // simulation end
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

		// ----------- predefinitions
		arma::vec times, entropy;
		double dt_new = 1e-2;
		std::string dir = this->saving_dir + "Entropy" + kPSep;
		std::function<void()> to_ave_time;
		alfa.reset_random();

		auto set_init_state = [N, this]() 
			-> arma::cx_vec 
		{
			arma::cx_vec init_state(N, arma::fill::zeros);
			switch (this->op) {
			case 0: // random product state
			{
				init_state = this->random_product_state(this->L); 
				break;
			}
			case 1: // ferromagnetically polarised
			{
				u64 idx = (ULLPOW(this->L)) - 1;
				init_state(idx) = cpx(1.0, 1.0); // 1111111
				break;
			}
			case 2: // anti-ferromagnetically polarised: 1010 + 0101
			{
				u64 idx = ((ULLPOW(this->L)) - 1) / 3;
				init_state(						idx) = cpx(1.0, 1.0); // 10101010
				init_state((ULLPOW(this->L)) -	idx) = cpx(1.0, 1.0); // 01010101
				break;
			}
			default:
				init_state = random_product_state(this->L); 
			}
			return arma::normalise(init_state);
		};
		
		stout << "\t\t	-->set random generators for " << alfa.get_info()
					<< " - in time : " << tim_s(start) << "s" << std::endl;
		// ----------- diagonalize
		stout << "\t\t	--> start diagonalizing for " << alfa.get_info()
				<< " - in time : " << tim_s(start) << "s" << std::endl;
		int realisation = 0;
		if(this->ch){
			dir += "Lanczos" + kPSep;
			int M = 200;
			auto H = alfa.get_hamiltonian();
			lanczos::Lanczos lancz(H, lanczosParams(M, 1, true, false));
			lancz.diagonalization();
			double omega_max = lancz.get_eigenvalues()(M - 1) - lancz.get_eigenvalues()(0);
			dt_new = 10 * this->dt / omega_max;
			//if (this->scale) stout << "WARNING: log only valid for t<10, for larger t linear is resumed with dt = " << dt_new << std::endl;
			//auto init_log = arma::logspace(-3, t_max, 4000);
			//auto rest_lin = arma::regspace(10.0, dt_new, tH);
			//times = this->scale? arma::join_cols(exctract_vector(init_log, 0.0, 10.0), rest_lin) : arma::regspace(dt_new, dt_new, tH);
			times = this->scale ? arma::logspace(-3, t_max, 500) : arma::regspace(dt_new, dt_new, tH);
			entropy = arma::vec(times.size(), arma::fill::zeros);

			to_ave_time = [&, lancz]() mutable 
				{	// capture lancz by value to access it outsied the if-else scope, i.e. when the lambda is called
					// entropy and times are declared outside the if-else scope so they can be passed by reference
					arma::cx_vec state = set_init_state();
					lancz.diagonalization(state);
					for (int i = 0; i < times.size(); i++)
					{
						auto t = times(i);
						//lancz.time_evolution_non_stationary(state, t - (i == 0 ? 0.0 : times(i - 1)), this->mu);
						auto evolved_state = lancz.time_evolution_stationary(state, t);
						entropy(i) += alfa.entaglement_entropy(evolved_state, this->L / 2);
					}
					stout << "realisation: " << ++realisation << " - in time : " << tim_s(start) << "s" << std::endl;
				};
		} else{
			alfa.diagonalization();
			double omega_max = alfa.get_eigenvalues()(N - 1) - alfa.get_eigenvalues()(0);
			dt_new = 10 * this->dt / omega_max;
			times = this->scale ? arma::logspace(-3, t_max, 500) : arma::regspace(dt_new, dt_new, tH);
			entropy = arma::vec(times.size(), arma::fill::zeros);
			to_ave_time = [&]()
				{
					arma::cx_vec init_state = set_init_state();
					for (int i = 0; i < times.size(); i++)
					{
						auto t = times(i);
						arma::cx_vec state = init_state;
						alfa.time_evolve_state(state, t);
						entropy(i) += alfa.entaglement_entropy(state, this->L / 2);
					}
					stout << "realisation: " << ++realisation << " - in time : " << tim_s(start) << "s" << std::endl;
				};
		}
		stout << "\t\t	--> finished diagonalizing for " << alfa.get_info()
					<< " - in time : " << tim_s(start) << "s" << std::endl;
		createDirs(dir);
		average_over_realisations<>(alfa, false, to_ave_time);
		entropy /= double(this->realisations);
		std::ofstream file;
		openFile(file, dir + "TimeEvolution" + alfa.get_info({}) + ".dat");
		for (int j = 0; j < times.size(); j++)
			printSeparated(file, "\t", 16, true, times(j), entropy(j));
		
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
	std::string s = (this->op < 3) ? "j" : "q";
	if(this->op == 6) s = "n";
	std::string dir = this->saving_dir + "RelaxationTimes" + kPSep;
	createDirs(dir);
	std::string op = IsingModel_disorder::opName(this->op, this->site);
	std::ofstream map_g, map_h;
	openFile(map_h, dir + "_h" + op + IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, {"h", "g"}) + ".dat", ios::out);
	// std::vector<double> gx_list = { 0.025, 0.05, 0.1, 0.15, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
	auto gx_list = arma::linspace(this->g, this->g + this->gs * (this->gn - 1), this->gn);
	auto hx_list = arma::linspace(this->h, this->h + this->hs * (this->hn - 1), this->hn);
	for (auto &gx : gx_list)
	{
		for (auto &hx : hx_list)
		{
			// read time-evolution data
			std::ifstream file;
			std::string name = op + IsingModel_disorder::set_info(this->L, this->J, this->J0, gx, this->g0, hx, this->w) + ".dat";
			std::string filename = this->saving_dir + "IntegratedResponseFunction" + kPSep + s + "=" + std::to_string(this->site) + kPSep + name;
			std::string filename_time = this->saving_dir + "TimeEvolution" + kPSep + s + "=" + std::to_string(this->site) + kPSep + name;
			std::string dir_out = this->saving_dir + "IntegratedResponseFunction" + kPSep + "DERIVATIVE" + kPSep + s + "=" + std::to_string(this->site) + kPSep;
			createDirs(dir_out);
			auto data = readFromFile(file, filename);
			file.close();
			if (data.empty())
				continue;
			// thouless time
			auto timeEvol = readFromFile(file, filename_time);
			auto i = timeEvol[1].index_min();
			double t_Th = timeEvol[0](i);
			file.close();
			// take derivative
			std::ofstream file2;
			openFile(file2, dir_out + name, std::ios::out);
			arma::vec specFun = non_uniform_derivative(data[0], data[1]);
			for (int j = 0; j < specFun.size(); j++)
				printSeparated(file2, "\t", 12, true, data[0](j + 1), specFun(j));
			file2.close();
			// find relax rate
			double wH = data[2](0);
			if (data[1](0) <= 0.5)
			{
				for (int k = 0; k < data[0].size(); k++)
				{
					if (data[1](k) >= 0.5)
					{
						printSeparated(map_h, "\t", 12, true, hx, gx, 1. / data[0](k), 1. / wH, t_Th);
						break;
					}
				}
			}
		}
	}
	map_h.close();
	openFile(map_g, dir + "_g" + op + IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, {"h", "g"}) + ".dat", ios::out);
	for (auto &hx : hx_list)
	{
		for (auto &gx : gx_list)
		{
			// read time-evolution data
			std::ifstream file;
			std::string name = op + IsingModel_disorder::set_info(this->L, this->J, this->J0, gx, this->g0, hx, this->w) + ".dat";
			std::string filename = this->saving_dir + "IntegratedResponseFunction" + kPSep + s + "=" + std::to_string(this->site) + kPSep + name;
			std::string filename_time = this->saving_dir + "TimeEvolution" + kPSep + s + "=" + std::to_string(this->site) + kPSep + name;
			std::string dir_out = this->saving_dir + "IntegratedResponseFunction" + kPSep + "DERIVATIVE" + kPSep + s + "=" + std::to_string(this->site) + kPSep;
			createDirs(dir_out);
			auto data = readFromFile(file, filename);
			file.close();
			if (data.empty())
				continue;
			// thouless time
			auto timeEvol = readFromFile(file, filename_time);
			auto i = timeEvol[1].index_min();
			double t_Th = timeEvol[0](i);
			file.close();
			// take derivative
			std::ofstream file2;
			openFile(file2, dir_out + name, std::ios::out);
			arma::vec specFun = non_uniform_derivative(data[0], data[1]);
			for (int j = 0; j < specFun.size(); j++)
				printSeparated(file2, "\t", 12, true, data[0](j + 1), specFun(j));
			file2.close();
			// find relax rate
			double wH = data[2](0);
			if (data[1](0) <= 0.5)
			{
				for (int k = 0; k < data[0].size(); k++)
				{
					if (data[1](k) >= 0.5)
					{
						printSeparated(map_g, "\t", 12, true, hx, gx, 1. / data[0](k), 1. / wH, t_Th);
						break;
					}
				}
			}
		}
	}
	map_g.close();
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

//<! spectral form factor calculated from eigenvalues in file or diagonalize matrix
void isingUI::ui::spectral_form_factor(){
    // values to be set in lambda
	double tH = 0;
	arma::vec eigenvalues;
	std::string info;

	// lambda 
	auto get_eigenvalues = [this, &tH, &eigenvalues, &info](auto& alfa)
	{
		double tH = 1. / alfa.mean_level_spacing_analytical();
		info = alfa.get_info({});
		std::ifstream energies;
		std::string E_dir = this->saving_dir + "EIGENVALUES" + kPSep;
		createDirs(E_dir);
		std::string name = E_dir + info;
		auto E = readFromFile(energies, name + ".dat");
		if(E.empty()){
			bool loaded = eigenvalues.load(arma::hdf5_name(name + ".h5", "eigenvalues"));
			if(!loaded){
				alfa.diagonalization(false);
				eigenvalues = alfa.get_eigenvalues();
				// save eigenvalues (yet unsaved)
				eigenvalues.save(arma::hdf5_name(name + ".h5", "eigenvalues"));
			}
		} else{
			eigenvalues = E[0];
		}
	};

	// ----------- choose model and run kernel
	if(this->m){
		auto alfa = std::make_unique<IsingModel_sym>(this->L, this->J, this->g, this->h,
								 this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, this->boundary_conditions);
		get_eigenvalues(*alfa);
	} else{
		auto alfa = std::make_unique<IsingModel_disorder>(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, this->boundary_conditions);
		get_eigenvalues(*alfa);
	}
	if(this->ch){
		// unfolding
	}
	int t_max = (int)std::ceil(std::log10(tH));
	t_max = (t_max / std::log10(tH) < 1.5) ? t_max + 1 : t_max;
	auto times = arma::logspace(-2, t_max, 300);
	arma::vec sff = spectrals::spectral_form_factor(eigenvalues, times);
	save_to_file(this->saving_dir + "SpectralFormFactor" + kPSep + info + ".dat", times, sff, tH);
}


//--------------------------------------------------------------------- ADIABATIC GAUGE POTENTIAL
void isingUI::ui::adiabaticGaugePotential_dis()
{
	clk::time_point start = std::chrono::system_clock::now();
	auto s = this->ch ? "h" : "g"; // (now not inversed because separated h sweep) //inversed, cause exclude the other one
	std::string info = IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, {"L"}, ",");
	std::string dir = this->saving_dir + "AGP" + kPSep;
	createDirs(dir);
	std::ofstream farante;
	openFile(farante, dir + IsingModel_disorder::opName(this->op, this->site) + info + "_" + s + ".dat");
	farante << std::setprecision(6) << std::scientific;
	// std::ofstream scaling(this->saving_dir + "AGPsize_DELETE" + info + ".dat");
	// scaling << std::setprecision(6) << std::scientific;
	auto params = this->ch ? arma::linspace(this->h, this->h + this->hs * this->hn, this->hn + 1)
							: arma::linspace(this->g, this->g + this->gs * this->gn, this->gn + 1);
	for (auto &x : params)
	{
		farante << x << "\t\t";
		// scaling << "\"h = " + to_string_prec(hx, 5) << "\"" << endl;
		for (int system_size = this->L; system_size < this->L + this->Ls * this->Ln; system_size += this->Ls)
		{
			const auto start_loop = std::chrono::system_clock::now();
			std::unique_ptr<IsingModel_disorder> alfa;
			if (this->ch)
				alfa = std::make_unique<IsingModel_disorder>(system_size, this->J, this->J0, this->g, this->g0, x, this->w);
			else
				alfa = std::make_unique<IsingModel_disorder>(system_size, this->J, this->J0, x, this->g0, this->h, this->w);
			stout << " \t\t--> finished creating model for " << alfa->get_info() << " - in time : " << tim_s(start_loop) << "\nTotal time : " << tim_s(start) << "s\n";
			auto opMat = alfa->chooseOperator(this->op, this->site);
			normaliseOp(opMat);
			std::function
				getValues = [&](double &AGP, double &typ_susc, int &counter)
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
				const arma::mat U = alfa->get_eigenvectors();
				arma::cx_mat mat_elem = U.t() * opMat * U;
#pragma omp parallel for reduction(+ \
								   : AGP_local, typ_susc_local)
				for (long int i = 0; i < N; i++)
				{
					double susc = 0;
					for (long int j = 0; j < N && j != i; j++)
					{
						const double nominator = abs(mat_elem(i, j) * conj(mat_elem(i, j)));
						const double omega_ij = alfa->get_eigenEnergy(j) - alfa->get_eigenEnergy(i);
						const double denominator = omega_ij * omega_ij + mu2 * mu2;
						AGP_local += omega_ij * omega_ij * nominator / (denominator * denominator);
						susc += nominator / (omega_ij * omega_ij);
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
				average_over_realisations<double, double &, double &, int &>(*alfa, true, getValues, AGP, typ_susc, counter);
			else
			{
				alfa->diagonalization();
				getValues(AGP, typ_susc, counter);
			}
			farante << typ_susc / double(counter) << "\t\t" << AGP / double(counter) << "\t\t";
			farante.flush();
			// scaling << L << "\t\t" << typ_susc << "\t\t" << AGP << "\t\t" << std::log(N) / log(2) << "\n";
			// scaling.flush();
		}
		farante << endl;
		// scaling << endl << endl;
	}
	farante.close();
	// scaling.close();
}
void isingUI::ui::adiabaticGaugePotential_sym(bool SigmaZ, bool avSymSectors)
{
	clk::time_point start = std::chrono::system_clock::now();
	std::string info = (avSymSectors ? IsingModel_sym::set_info(L, J, g, h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"L", "h", "k", "p", "x"}) : IsingModel_sym::set_info(L, J, g, h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, {"L", "h"}));
	std::string dir = this->saving_dir + "AGP" + kPSep;
	createDirs(dir);
	std::ofstream farante(dir + (SigmaZ ? "SigZ" : "SigX") + info + ".dat");
	farante << std::setprecision(6) << std::scientific;
	// std::ofstream scaling(this->saving_dir + "AGPsize_DELETE" + info + ".dat");
	// scaling << std::setprecision(6) << std::scientific;
	farante << "hx\t\t\t susceptibiltiy\t\tAGP,\t\tsusceptibiltiy\tAGP,\t\t susceptibiltiy\tAGP,\t\t susceptibiltiy\tAGP,\t\t susceptibiltiy\tAGP\n";
	auto params = arma::linspace(0.02, 3, 150);
	for (auto &hx : params)
	{
		farante << hx << "\t\t";
		// scaling << "\"h = " + to_string_prec(hx, 5) << "\"" << endl;
		stout << "\nh = " << hx << "\t\t";
		for (int system_size = this->L; system_size < this->L + this->Ls * this->Ln; system_size += this->Ls)
		{
			const auto start_loop = std::chrono::system_clock::now();
			// auto alfa = std::make_unique<IsingModel_disorder>(system_size, this->J, 0, this->g, 0, hx, 0, 0);
			std::function
				getValues = [&](int k, int p, int x, double &AGP, double &typ_susc, int &counter)
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
				// const double mu2 = std::log2(N) / double(N);
				static long int E_min = 0;		 // alfa->E_av_idx - mu / 2.;
				static long int E_max = (long)N; // alfa->E_av_idx + mu / 2.;
				double typ_susc_local = 0;
				double AGP_local = 0;
				int counter_tmp = 0;
				const arma::cx_mat U = alfa->get_eigenvectors();
				arma::sp_cx_mat opMatrix = SigmaZ ? alfa->create_operator({IsingModel_sym::sigma_z}) : alfa->create_operator({IsingModel_sym::sigma_x});
				cpx norm = arma::trace(opMatrix * opMatrix) / double(N);
				if (abs(norm) > 1e-12)
					opMatrix /= norm;
				arma::cx_mat mat_elem = U.t() * opMatrix * U;
				for (long int i = 0; i < N; i++)
				{
					double susc = 0;
					for (long int j = 0; j < N && j != i; j++)
					{
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
				loopSymmetrySectors<double &, double &, int &>(getValues, hx, system_size, AGP, typ_susc, counter);
			else
				getValues(this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, AGP, typ_susc, counter);
			farante << typ_susc / double(counter) << "\t\t" << AGP / double(counter) << "\t\t";
			farante.flush();
			// scaling << L << "\t\t" << typ_susc << "\t\t" << AGP << "\t\t" << std::log(N) / log(2) << "\n";
			// scaling.flush();
		}
		farante << endl;
		// scaling << endl << endl;
	}
	farante.close();
	// scaling.close();
}
void isingUI::ui::combineAGPfiles()
{
	std::ifstream input;
	std::ofstream output;
	std::string line;
	std::string dir = "SPINONsynchro" + kPSep;
	for (double gx = 0.05; gx <= 0.20; gx += 0.05)
	{
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

template <typename _type>
void isingUI::ui::LevelSpacingDist(IsingModel<_type> &alfa)
{
	auto r = alfa.eigenlevel_statistics_with_return();
	auto probDist = probability_distribution_with_return(r);
	double _min = arma::min(r);
	double _max = arma::max(r);
	for (int i = 0; i < 20; i++)
	{
		double dg = alfa.getRandomValue(-alfa.g / 100., alfa.g / 100.);
		auto beta = std::make_shared<IsingModel<_type>>(alfa);
		beta->g = beta->g + dg;
		beta->hamiltonian();
		beta->diagonalization(false);
		auto level_spacing = beta->eigenlevel_statistics_with_return();
		probDist += probability_distribution_with_return(level_spacing);
		_min += arma::min(level_spacing);
		_max += arma::max(level_spacing);
	}
	_min /= 21.;
	_max /= 21.;
	probDist /= 21.;
	fs::create_directory(this->saving_dir + "LevelSpacing" + kPSep);
	save_to_file(this->saving_dir + "LevelSpacing" + kPSep, "LevelSpacing" + alfa->get_info({}),
				 arma::linspace(_min, _max, probDist.size()), probDist);
}


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




