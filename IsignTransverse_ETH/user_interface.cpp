#include "include/user_interface.h"

/// <summary>
/// We want to handle files so let's make the c-way input a string
/// </summary>
/// <param name="argc"> number of main input arguments </param>
/// <param name="argv"> main input arguments </param>
/// <returns></returns>
std::vector<std::string> change_input_to_vec_of_str(int argc, char** argv)
{
	std::vector<std::string> tmp(argc - 1, "");															// -1 because first is the name of the file
	for (int i = 0; i < argc - 1; i++) {
		tmp[i] = argv[i + 1];
	}
	return tmp;
}

// - - - - - - - - - - - - - - - - - - - USER INTERFACE - - - - - - - - - - - - - - - - - - -

/// <summary>
/// Find a given option in a vector of string given from cmd parser
/// </summary>
/// <param name="vec">vector of strings from cmd</param>
/// <param name="option">the option that we seek</param>
/// <returns>value for given option if exists, if not an empty string</returns>
std::string user_interface::getCmdOption(const v_1d<std::string>& vec, std::string option) const
{
	if (auto itr = std::find(vec.begin(), vec.end(), option); itr != vec.end() && ++itr != vec.end())
		return *itr;
	return std::string();
}
/// <summary>
///
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="value"></param>
/// <param name="argv"></param>
/// <param name="choosen_option"></param>
/// <param name="geq_0"></param>
template<typename T>
void user_interface::set_option(T& value, const v_1d<std::string>& argv, std::string choosen_option, bool geq_0)
{
	if (std::string option = this->getCmdOption(argv, choosen_option); option != "")
		value = static_cast<T>(stod(option));												// set value to an option
	if (geq_0 && value < 0)																	// if the variable shall be bigger equal 0
		this->set_default_msg(value, choosen_option.substr(1), \
			choosen_option + " cannot be negative\n", isingUI::table);
}
/// <summary>
///
/// </summary>
/// <typeparam name="T"></typeparam>
/// <param name="value"></param>
/// <param name="option"></param>
/// <param name="message"></param>
/// <param name="map"></param>
template<typename T>
void user_interface::set_default_msg(T& value, std::string option, std::string message, \
	const std::unordered_map<std::string, std::string>& map) const
{
	stout << message;																			// print warning
	std::string value_str = "";																// we will set this to value
	auto it = map.find(option);
	if (it != map.end()) {
		value_str = it->second;															// if in table - we take the enum
	}
	value = stod(value_str);
}
// - - - - - - - - - - - - - - - - - - - ISING MODEL - - - - - - - - - - - - - - - - - - -

/* Connected with the parser */
/// <summary>
/// Setting the default parameters for the Ising model
/// </summary>
void isingUI::ui::set_default() {
	using namespace isingUI;
	this->saving_dir = "." + std::string(kPathSeparator) + "results" + std::string(kPathSeparator);		// directory for the result files to be saved into
	this->L = 4;

	this->J = 1.0;
	this->J0 = 0.2;
	this->h = 0.0;
	this->w = 1.0;
	this->g = 1.0;
	this->g0 = 0;

	this->realisations = 100;
	this->site = 0;
	this->mu = 5;

	this->boundary_conditions = 0;
	this->m = 0;
	this->p = true;
	this->thread_number = 1;
}

// ---- CONSTURCTORS

/// <summary>
///
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
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
		"-J spin exchange coefficient : (default 1)\n"
		"-J0 random spin exchange set in uniform distribution [-J0,J0]\n"
		"-g transverse magnetic field (x-) constant: (default 1)\n"
		"-g0 random transverse field set in uniform distribution [-g0,g0]\n"
		"-h perpendicular (z-) magnetic field constant: (default 0)\n"
		"-w disorder strength : (default 0 - no disorder introduced)\n"
		"-L chain length : bigger than 0 (default 8)\n"
		"-b boundary conditions : bigger than 0 (default 0 - PBC)\n"
		"	0 -- PBC\n"
		"	1 -- OBC\n"
		"	2 -- ABC\n"
		"-m model to be choosen : (default 0 - without symmetries)\n"
		"	0 -- nonsymmetric model - only here the disorder is working\n"
		"	1 -- include symmetries - here the parity flag is also working\n"
		"-th number of threads to be used for CPU parallelization : depends on the machine specifics, default(1)"
		"-h quit with help\n"
	);
	std::exit(1);
}
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
	// transverse field disorder
	choosen_option = "-g0";
	this->set_option(this->g0, argv, choosen_option);
	// perpendicular field
	choosen_option = "-h";
	this->set_option(this->h, argv, choosen_option);
	// perpendicular field disorder
	choosen_option = "-w";
	this->set_option(this->w, argv, choosen_option);

	// chain length
	choosen_option = "-L";
	this->set_option(this->L, argv, choosen_option);
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
		<< thread_number << " - - - - - - " << endl;				// setting the number of threads to be used with omp

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

// ---- SIMULATIONS
void isingUI::ui::compare_energies() {
	auto Hamil = std::make_unique<IsingModel_disorder>(L, J, J0, g, g0, h, w, boundary_conditions);
	u64 N = Hamil->get_hilbert_size();
	Hamil->diagonalization();
	vec E_dis = Hamil->get_eigenvalues();
	std::vector<double> E_sym = v_1d<double>();
	std::vector<std::string> symmetries = v_1d<std::string>();
	for (int k = 0; k < L; k++) {
		if (k == 0 || k == this->L / 2.) {
			for (int p = 0; p <= 1; p++) {
				for (int x = 0; x <= 1; x++) {
					auto Hamil = std::make_unique<IsingModel_sym>(L, J, g, h, k, p, x, boundary_conditions);
					Hamil->diagonalization();
					vec t = Hamil->get_eigenvalues();
					E_sym.insert(E_sym.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));
					v_1d<std::string> temp_str = v_1d<std::string>(t.size(), "k=" + std::to_string(k) + ",x=" + to_string(x) + ",p=" + to_string(p));
					symmetries.insert(symmetries.end(), std::make_move_iterator(temp_str.begin()), std::make_move_iterator(temp_str.end()));
				}
			}
		}
		else {
			for (int x = 0; x <= 1; x++) {
				auto Hamil = std::make_unique<IsingModel_sym>(L, J, g, h, k, 1, x, boundary_conditions);
				Hamil->diagonalization();
				vec t = Hamil->get_eigenvalues();
				E_sym.insert(E_sym.end(), std::make_move_iterator(t.begin()), std::make_move_iterator(t.end()));
				v_1d<std::string> temp_str = v_1d<std::string>(t.size(), "k=" + std::to_string(k) + ",x=" + to_string(x));
				symmetries.insert(symmetries.end(), std::make_move_iterator(temp_str.begin()), std::make_move_iterator(temp_str.end()));
			}
		}
		//out << Hamil->get_eigenvalues().t();
	}
	auto p = sort_permutation(E_sym, [](const double a, const double b) {
		return a < b;
		});
	//sort(E_sym.begin(), E_sym.end());
	apply_permutation(E_sym, p);
	apply_permutation(symmetries, p);
	stout << E_sym.size() << endl;
	for (int k = 0; k < min(E_sym.size(), E_dis.size()); k++) {
		stout << symmetries[k] << "\t\t\t\t" << E_sym[k] << "\t\t\t\t" << E_dis(k) << "\t\t\t\t" << E_sym[k] - E_dis(k) << endl;
	}
}

void isingUI::ui::disorder() {
	/*auto start = std::chrono::high_resolution_clock::now();
	std::ofstream scaling_r_sigmaX(this->saving_dir + "SpectrumRapScalingSigmaX" + \
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

	stout << "\n--> starting loop over disorders <--\n";
	L = 10;
	std::ofstream scaling_ipr(this->saving_dir + "iprDisorder" + \
		"_L=" + std::to_string(this->L) + \
		",J0=" + to_string_prec(this->J0, 2) + \
		",g=" + to_string_prec(this->g, 2) + \
		",g0=" + to_string_prec(this->g0, 2) + \
		",h=" + to_string_prec(this->h, 2) + \
		".dat", std::ofstream::app);
	for (double w = 0.0; w <= 6.0; w += 0.1) {
		realisations = 600;

		auto Hamil = std::make_unique<IsingModel_disorder>(L, J, J0, g, g0, h, w);
		const u64 N = Hamil->get_hilbert_size();
		Hamil->diagonalization();

		vec av_sigma_x = Hamil->operator_av_in_eigenstates_return(&IsingModel::av_sigma_x, *Hamil, 0);
		vec fluct = data_fluctuations(av_sigma_x);
		double _min = -2.0, _max = 2.0, step = 2e-3;
		stout << "--> finished writing the sigma _x fluctuations for w = " << w << " <--\n";

		arma::vec prob_dist = probability_distribution_with_return(fluct, _min, _max, step);
		arma::vec prob_dist_GOE = probability_distribution_with_return(Hamil->eigenlevel_statistics_with_return(), 0, 1, 2 * step);
		double ipr = 0, entropy = 0;
		int mu = 100;
		double r = 0;
		for (int k = 0; k < realisations - 1; k++) {
			Hamil->hamiltonian();
			Hamil->diagonalization();
			av_sigma_x = Hamil->operator_av_in_eigenstates_return(&IsingModel::av_sigma_x, *Hamil, 0);
			fluct = data_fluctuations(av_sigma_x);

			prob_dist += probability_distribution_with_return(fluct, _min, _max, step);
			prob_dist_GOE += probability_distribution_with_return(Hamil->eigenlevel_statistics_with_return(), 0, 1, 2 * step);

			// average in middle spectrum
			for (int f = Hamil->E_av_idx - mu / 2.; f < Hamil->E_av_idx + mu / 2.; f++) {
				ipr += Hamil->ipr(f);
				entropy += Hamil->information_entropy(f);
				r += Hamil->eigenlevel_statistics(f, f + 1);
			}
			//if (k % 5 == 0) stout << " \t\t--> " << k << " - in time : " << \
				double(duration_cast<milliseconds>(duration(high_resolution_clock::now() - start)).count()) / 1000.0 << "s" << std::endl;
		}
		stout << "--> finished loop over realisations for w = " << w << " <--\n";
		std::ofstream ProbDistSigmaX(this->saving_dir + "ProbDistSigmaX" + Hamil->get_info() + ".dat");
		for (int f = 0; f < prob_dist.size(); f++)
			ProbDistSigmaX << _min + f * step << "\t\t\t\t" << prob_dist(f) / double(realisations) << endl;
		ProbDistSigmaX.close();

		std::ofstream ProbDistGap(this->saving_dir + "ProbDistGap" + Hamil->get_info() + ".dat");
		for (int f = 0; f < prob_dist_GOE.size(); f++)
			ProbDistGap << f * 2 * step << "\t\t\t\t" << prob_dist_GOE(f) / double(realisations) << endl;
		ProbDistGap.close();

		double norm = realisations * mu;
		scaling_ipr << w << "\t\t\t\t" << ipr / norm / (double)N << "\t\t\t\t" << entropy / norm << "\t\t\t\t" << r / norm << endl;
		stout << " \t\t--> w = " << w << " - in time : " << \
			double(duration_cast<milliseconds>(duration(high_resolution_clock::now() - start)).count()) / 1000.0 << "s" << std::endl;
	}
	scaling_ipr.close();*/
 }

void isingUI::ui::compare_matrix_elements() {
	std::ofstream file("file.dat");

	auto model = std::make_unique<IsingModel_disorder>(L, J, 0, g, 0, h, 0);
	model->diagonalization();
	//vec E_dis = model->get_eigenvalues();
	//auto E_unique = get_NonDegenerated_Elements(E_dis);
	std::unique_ptr<IsingModel_sym> alfa = std::make_unique<IsingModel_sym>(L, J, g, h, 0, 1, 1);
	std::unique_ptr<IsingModel_sym> beta = std::make_unique<IsingModel_sym>(L, J, g, h, L / 2.0, 0, 0);
	stout << alfa->get_hilbert_size() << endl;
	alfa->diagonalization();
	beta->diagonalization();
	//stout << model->get_eigenvalues() << endl;
	auto map_alfa = mapping_sym_to_original(0, model->get_hilbert_size() - 1, *alfa, *model);
	auto map_beta = mapping_sym_to_original(0, model->get_hilbert_size() - 1, *beta, *model);
	stout << "alfa sector size for nondegenerate mapping is : " << map_alfa.size() << endl;
	stout << "beta sector size for nondegenerate mapping is : " << map_alfa.size() << endl;

	auto sig_x = IsingModel_sym::sigma_x;
	auto sig_z = IsingModel_sym::sigma_z;
	auto spin_cur = IsingModel_sym::spin_flip;
	stout << "ALPHA SECTOR : k=0,x=1,p=1, BETA SECTOR : k=pi,x=1,p=1\n\n";
	stout << " - - - - - - SAME SECTORS - - - - - - \n" << " - - - - > FOR ALFA - ALFA: \n" << "ENERGY ALPHA |('.'|)" << "\t\t" << \
		"ENERGY ALFA (/'.')/" << "\t\t" << "<alfa|SIGMA_X|alfa>" << "\t\t" << "<non_sym|SIGMA_X|non_sym>" << "\t\t" << "DIFFERENCE" << endl;
	for_each(exec::seq, map_alfa.begin(), map_alfa.end(), [&](std::pair<u64, u64> element) {
		for (auto& t : map_alfa) {
			cpx A = alfa->av_spin_current(element.first, t.first, { 1,2 });
			cpx B = model->av_spin_current(element.second, t.second, { 1,2 });
			stout << alfa->get_eigenEnergy(element.first) << "\t\t" << \
				alfa->get_eigenEnergy(t.first) << "\t\t" << real(A) << "\t\t" << real(B) << "\t\t" << real(abs(A) - abs(B)) << endl;
		}
		});
	stout << "\n - - - - > FOR BETA - BETA: \n" << "ENERGY BETA |('.'|)" << "\t\t" << \
		"ENERGY BETA (/'.')/" << "\t\t" << "<beta|SIGMA_X|beta>" << "\t\t" << "<non_sym|SIGMA_X|non_sym>" << "\t\t" << "DIFFERENCE" << endl;
	for_each(exec::seq, map_beta.begin(), map_beta.end(), [&](std::pair<u64, u64> element) {
		for (auto& t : map_beta) {
			cpx A = beta->av_spin_current(element.first, t.first, { 1,2 });
			cpx B = model->av_spin_current(element.second, t.second, { 1,2 });
			file << beta->get_eigenEnergy(element.first) << "\t\t" << \
				beta->get_eigenEnergy(t.first) << "\t\t" << real(A) << "\t\t" << real(B) << "\t\t" << real(abs(A) - abs(B)) << endl;
		}
		});
	stout << "\n\n - - - - - - DIFFERENT SECTORS - - - - - - \n" << " - - - - > FOR ALFA - BETA: \n" << "ENERGY ALPHA |('.'|)" << "\t\t" << \
		"ENERGY BETA (/'.')/" << "\t\t" << "<alfa|SIGMA_X|beta>" << "\t\t" << "<non_sym_a|SIGMA_X|non_sym_b>" << "\t\t" << "DIFFERENCE" << endl;
	for_each(exec::seq, map_alfa.begin(), map_alfa.end(), [&](std::pair<u64, u64> element) {
		for (auto& t : map_beta) {
			cpx A = im * av_operator(element.first, t.first, *alfa, *beta, spin_cur, { 1,2 });
			A += conj(im * av_operator(t.first, element.first, *beta, *alfa, spin_cur, { 1,2 }));
			A *= 0.5i;
			cpx B = model->av_spin_current(element.second, t.second, { 1,2 });
			stout << alfa->get_eigenEnergy(element.first) << "\t\t" << \
				beta->get_eigenEnergy(t.first) << "\t\t" << real(A) << "\t\t" << real(B) << "\t\t" << real(abs(A) - abs(B)) << endl;
		}
		});
	file.close();


}

void isingUI::ui::size_scaling_sym(int k, int p, int x, int L_min, int L_max) {
	using namespace std::chrono;
	auto start = std::chrono::high_resolution_clock::now();

	std::ofstream farante(this->saving_dir + "IprScaling" + \
		",g=" + to_string_prec(this->g, 2) + \
		",h=" + to_string_prec(this->h, 2) + ".dat");
	std::ofstream fikolo(this->saving_dir + "SpectrumRapScalingSigmaX" + \
		",g=" + to_string_prec(this->g, 2) + \
		",h=" + to_string_prec(this->h, 2) + ".dat");
	
	for (L = L_min; L <= L_max; L ++) {
		stout << "\n\n------------------------------Doing L = " << L << "------------------------------\n";
		auto alfa = std::make_unique<IsingModel_sym>(L, J, g, h, k, p, x, boundary_conditions);
		u64 N = alfa->get_hilbert_size();
		alfa->diagonalization();
		stout << " \t\t--> finished diagonalizing for " << alfa->get_info() << " - in time : " << \
			double(duration_cast<milliseconds>(duration(high_resolution_clock::now() - start)).count()) / 1000.0 << "s" << std::endl;

		// average sigma_x operator at first site
		std::ofstream sigx(this->saving_dir + "SigmaX" + alfa->get_info() + ".dat");
		vec av_sigma_x(N, fill::zeros);
		for (int k = 0; k < N; k++) {
			av_sigma_x(k) = alfa->av_sigma_x(k, k, { 0 });
			sigx << alfa->get_eigenEnergy(k) << "\t\t" << av_sigma_x(k) << endl;
		}
		sigx.close();

		//outliers and prob distribution for sigma_x
		vec r_sigma_x(N - 1);
#pragma omp parallel for
		for (int k = 0; k < N - 1; k++)
			r_sigma_x(k) = abs(av_sigma_x(k + 1) - av_sigma_x(k));
		vec outliers = statistics_average(r_sigma_x, 4);
		fikolo << L << "\t\t" << outliers.t();
		stout << " \t\t--> finished outliers for " << alfa->get_info() << " - in time : " << \
			double(duration_cast<milliseconds>(duration(high_resolution_clock::now() - start)).count()) / 1000.0 << "s" << std::endl;

		probability_distribution(this->saving_dir, "ProbDistSpecRapSigmaX" + alfa->get_info(), r_sigma_x, 0, 0.05, 0.0002);
		probability_distribution(this->saving_dir, "ProbDistSigmaX" + alfa->get_info(), data_fluctuations(av_sigma_x), -2.0, 2.0, 0.001);
		stout << " \t\t--> finished prob dist for " << alfa->get_info() << " - in time : " << \
			double(duration_cast<milliseconds>(duration(high_resolution_clock::now() - start)).count()) / 1000.0 << "s" << std::endl;

		// eigenlevel statistics and prob distribution
		double r = alfa->eigenlevel_statistics(alfa->E_av_idx - mu / 2., alfa->E_av_idx + mu / 2.);

		// ipr & info entropy
		double ipr = 0, ent = 0;
		for (int k = alfa->E_av_idx - mu / 2.; k <= alfa->E_av_idx + mu / 2.; k++) {
			ipr += alfa->ipr(k);
			ent += alfa->information_entropy(k);
		}
		farante << L << "\t\t" << ipr / double(mu) / double(N) << "\t\t" << ent / double(mu) << "\t\t" << r << endl;
		stout << " \t\t--> finished farante for " << alfa->get_info() << " - in time : " << \
			double(duration_cast<milliseconds>(duration(high_resolution_clock::now() - start)).count()) / 1000.0 << "s" << std::endl;

	}
	fikolo.close();
	farante.close();
	auto stop = std::chrono::high_resolution_clock::now();
	stout << " - - - - - - FINISHED SIZE SCALING for:\nk = " << k << ", p = " << p << ", x = " << x << "\nIN : " << \
		double(duration_cast<milliseconds>(duration(stop - start)).count()) / 1000.0 << " seconds - - - - - - " << endl;
}

void isingUI::ui::make_sim()
{
	using namespace std::chrono;
	auto start = std::chrono::high_resolution_clock::now();

	//compare_energies();
	this->mu = 10;
	size_scaling_sym(0, 1, 1, 12, 18);
	size_scaling_sym(1, 1, 1, 12, 18);
	//if(h == 0) size_scaling_sym(1, 1, -1, 12, 18, 2);
	

	auto stop = std::chrono::high_resolution_clock::now();
	stout << " - - - - - - FINISHED CALCULATIONS IN : " << \
		double(duration_cast<milliseconds>(duration(stop - start)).count()) / 1000.0 << " seconds - - - - - - " << endl;						// simulation end
}