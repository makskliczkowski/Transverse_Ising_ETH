#include "include/user_interface.h"

/// <summary>
/// We want to handle files so let's make the c-way input a string
/// </summary>
/// <param name="argc"> number of main input arguments </param>
/// <param name="argv"> main input arguments </param>
/// <returns></returns>
std::vector<std::string> change_input_to_vec_of_str(int argc, char** argv)
{
	std::vector<std::string> tmp(argc-1,"");															// -1 because first is the name of the file
	for(int i = 1; i <argc;i++ ){
		tmp[i] = argv[i];
	}
	return tmp;
}



/* - - - - - - - - - - - - - - - - - - - ISING MODEL - - - - - - - - - - - - - - - - - - -  */

/* Connected with spin exchange */
/// <summary>
/// Setting J according to the choosen mode
/// </summary>
/// <param name="mode"> mode - 0 (from file), 1 (random), 2 (constant)</param>
/// <param name="argument"> can be either number, or "r", or the file </param>
/// <returns></returns>
void isingUI::ui::set_J(int mode, std::string argument){
	using namespace isingUI;
	std::vector<double> tempJ;
	std::ifstream file;
	switch (mode)
	{
	case 0:
	// is choosen from file
		file.open(argument);
		if (!file.is_open()){
			out << "Can't open file " + argument + ".\nSetting J to default.\n";
			std::fill(tempJ.begin(), tempJ.end(), 1);
		}
		else{
			// we will read the system size first
			// then the line shoud be separated by spaces
			std::string line;
			// get L
			if(getline(file,line) && isNumber(line)){
				this->L = stoi(line);
			}
			else{
				out << "Don't understand the file: " + argument +  ". Setting J to default\n";
				this->set_J(2,"1");
				break;
			}
			// get J
			if(getline(file,line)){
				// split string to each element
				std::vector<string> Jstring = split_str(line," ");
				for(auto element: Jstring){
					this->J = std::vector<double>();
					// pushback element if it can be pushed
					if(isNumber(element)) this->J.push_back(stod(element));
				}
			}
			correct_J_size(mode);
		}
		break;
	case 1:
	// is choosen to be random
		tempJ = create_random_vec_std(this->L);
		break;
	// is choosen to be constant
	case 2:
		std::fill(tempJ.begin(), tempJ.end(), std::stod(argument));
		break;
	default:
		std::fill(tempJ.begin(), tempJ.end(), 1);
		break;
	}
	this->J = tempJ;
}
/// <summary>
/// Corrects J size according to the length of the Ising chain isingUI::L
/// </summary>
/// <param name="mode"> Can be "constant" 1, "random" 2, "file" 0</param>
void isingUI::ui::correct_J_size(int mode){
	if(this->J.size() == 0){
		std::fill(this->J.begin(),this->J.end(),1.0);			// if is empty set to default value according to isingUI::L
	}
	else{
		std::uniform_real_distribution<double> uni_dist(-1.0,1.0);
		while(this->J.size() < this->L){
		// while L is bigger we add elements
			double temp_value = 1;
			if(mode == 1)
				temp_value = uni_dist(gen);
			else
				temp_value = this->J[0];
		}
		while(this->J.size() > this->L && this->J.size() > 0)
		// while smaller we popback
			this->J.pop_back();
	}
}
/* Connected with the parser */
/// <summary>
/// Setting the default parameters for the Ising model
/// </summary>
void isingUI::ui::set_default(){
	using namespace isingUI;
	saving_dir = "";
	L = 8;
	J = std::vector<double>(L,1.0);
	h = 0;
	g = 1;
	w = 0;
	boundary_conditions = 0;
	m = 0;
	p = true;
	thread_number = 1;
}
/// <summary>
/// Function that tells how does the parser work
/// </summary>
void isingUI::ui::exit_with_help(){
		printf(
		" Usage: name of executable [options] outputDirectory \n"
		" The input can be both introduced with [options] described below or with giving the input directory(which also is the flag in the options)\n"
		" options:\n"
		"-f input file for all of the options : (default none)\n"
		"-J spin exchange coefficient : (default 1)\n"
		"	directory -- (mode 0) -- from file : set coefficients for different lattice sites seperately\n"
		"		-> file should be prepared giving the system size L in the first line and in the second line the J vector separated by spaces\n"
		"	(is numeric) -- (mode 2) -- constant value\n"
		"	r -- (mode 1) -- random value\n"
		"-h perpendicular (z-) magnetic field constant: (default 0)\n"
		"-g transverse magnetic field (x-) constant: (default 1)\n"
		"-w disorder strength : (default 0 - no disorder introduced)\n"
		"-L chain length : bigger than 0 (default 8)\n"
		"-b boundary conditions : bigger than 0 (default 0 - PBC)\n"
		"	0 -- PBC\n"
		"	1 -- OBC\n"
		"	2 -- ABC\n"
		"-m model to be choosen : (default 0 - without symmetries)\n"
		"	0 -- nonsymmetric model - only here the disorder is working\n"
		"	1 -- include symmetries - here the parity flag is also working\n"
		"-p include the parity symmetry : workis only in the model with symmetries (1) (default 1 - with parity symmetry)\n"
		"-th number of threads to be used for CPU parallelization : depends on the machine specifics, default(1)"
		"-q quit with help\n"
	);
	std::exit(1);
}
/// <summary>
/// The parser for the Transverse Field Ising model
/// </summary>
/// <param name="argc"> number of arguments </param>
/// <param name="argv"> list of arguments </param>
void isingUI::ui::parseModel(int argc, std::vector<std::string> argv){
	using namespace isingUI;
	bool reading_from_file = false;																		// if reading_from_file is set to true, the loop on argument breaks instantly, it's better to have -f at the beginning then ;)
	std::vector<std::string> command;																	// if the file for commands is present, then it will be used
	std::string commandFileName;																		
	int Jtype = 0;																						// saving the type of J if L doesn't match the J size. Can be "constant" 2, "random" 1, "file" 0
	this->saving_dir ="." + std::string(kPathSeparator) + "results" + std::string(kPathSeparator);		// directory for the result files to be saved into
	/* SET DEFAULT VALUES */
	this->set_default();																				// setting default at the very beginning
	/* STARTING ARGUMENT LOOP */
	int i = 1;
	for (i = 1; i < argc; i++)
	{
		/* BREAKERS */
		if (argv[i][0] != '-') break;
		if (reading_from_file) break;
		if (++i >= argc)
			this->exit_with_help();
		/* PARSE COMANDS */
		std::string argument = (string(argv[i - 1])).substr(1, string(argv[i - 1]).size() - 1);			// taking the argument to string
		auto it = table.find(argument);																	// looking for argument iterator in map parser
		parsers enum_arg;																				// creating an instance of the enum class for switch-case
		if (it != table.end()) {
			enum_arg = it->second;																		// if in table - we take the enum 
		}
		else {
			enum_arg = parsers::q;
			fprintf(stderr, "Unknown option: -%c\n", argv[i - 1][1]);									// exit if item is not in the parser
			exit_with_help();
		}
		/* PARSE THE ARGUMENT */
		switch (enum_arg)
		{
		default:
			fprintf(stderr, "Unknown option: -%c\n", argv[i - 1][1]);
			//cout << "Setting default training model parameters\n";
			exit_with_help();
			break;
		case isingUI::parsers::f:
			reading_from_file = true;
			commandFileName = argv[i];
			break;
		case isingUI::parsers::J:
			set_J(Jtype, argv[i]);
			break;
		case isingUI::parsers::h:
			
			if (isNumber(argv[i])) {
				h = stod(argv[i]);
			}
			else{
				out << "Bad input for h. Setting default\n";
				h = 0;
			}
			break;
		case isingUI::parsers::L:
			if (isNumber(argv[i])) {
				L = stoi(argv[i]);
			}
			else{
				out << "Bad input for L. Setting default\n";
				L=8;
			}
			// checking if is correct for J vector size
			correct_J_size(Jtype);
			break;
		case isingUI::parsers::g:
			
			if (isNumber(argv[i])) {
				g = stod(argv[i]);
			}
			else{
				out << "Bad input for g. Setting default\n";
				g = 1;
			}
			break;
		case isingUI::parsers::w:
			
			if (isNumber(argv[i])) {
				w = stod(argv[i]);
			}
			else{
				out << "Bad input for h. Setting default\n";
				w = 0;
			}
			break;
		case isingUI::parsers::b:
			if (isNumber(argv[i])) {
				this->boundary_conditions = stoi(argv[i]);
			}
			else
				this->boundary_conditions = 0;
			break;
		case isingUI::parsers::m:
			if (isNumber(argv[i])) {
				m = stoi(argv[i]);
			}
			else
				m = 0;
			break;
		case isingUI::parsers::p:
			if (isNumber(argv[i])) {
				p = bool(stoi(argv[i]));
			}
			else 
				p = true;
			break;
		case isingUI::parsers::th:
			if (isNumber(argv[i])) {
				int thread_num = stoi(argv[i]);
				if(thread_num > 0 && thread_num <= thread::hardware_concurrency())
					thread_number = thread_num;
				else
					thread_number = 1;
			}
			break;
		case isingUI::parsers::q:
			this->exit_with_help();
			break;

		}
	}
	if ((!((argv[i]).size()))==0)
	{
		saving_dir = string(argv[i]);
		std::filesystem::create_directories(saving_dir);													// creating the directory for saving the files with results
		saving_dir = saving_dir + std::string(kPathSeparator);
	}
	if(reading_from_file){
		command = this->parseInputFile(argv[i]);															// turn file into command line - like type
		if(command.size()==0){
			out << "Wrong file. Setting default\n";
			this->set_default();
		}
		this->parseModel(command.size(), command);
	}

	cout << " - - - - - - USING OUTER THREADS : " << thread_number << "- - - - - - " << endl;				// setting the number of threads to be used with omp
	omp_set_num_threads(this->thread_number);
	return;
}
/// <summary>
/// If the commands are given from file, we must treat them the same as arguments
/// </summary>
/// <param name="filename"> the name of the file that contains the command line </param>
/// <returns></returns>
std::vector<std::string> user_interface::parseInputFile(std::string filename){
	std::vector<std::string> commands;
	ifstream inputFile(filename);
	std::string line = "";
	if(!inputFile.is_open()){
		out << "Cannot open a file " + filename + " that I could parse. Setting all parameters to default. Sorry :c \n";
		this->set_default();
	}
	else{
		if(std::getline(inputFile, line)){
			commands = split_str(line, " ");														// saving lines to out vector if it can be done, then the parser shall treat them normally
		}
	}
	return std::vector<std::string>(commands.begin(),commands.end()); 
}
