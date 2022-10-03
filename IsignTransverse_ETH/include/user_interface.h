#pragma once
#ifndef UI
#define UI
//#include "headers.h"
#include "IsingModel.h"
#include "statistics.hpp"
#include "statistics_dist.hpp"
#include "spectrals.hpp"
#include "thermodynamics.hpp"
#include "entaglement.hpp"
#include "adiabatic_gauges.hpp"
#include "ThomasAlgorithm.hpp"


const arma::vec down = { 0, 1 };
const arma::vec up	 = { 1, 0 };
extern std::uniform_real_distribution<> theta;
extern std::uniform_real_distribution<> fi;
extern int outer_threads;
// can't be const cause operator() is for non-const only (on microsoft it can be on const)
using namespace std;
std::vector<std::string> change_input_to_vec_of_str(int argc, char** argv);
// ----------------------------------------------------------------------------- GENERAL CLASS -----------------------------------------------------------------------------
class user_interface {
protected:
	unsigned int thread_number = 1;																								// number of threads
	int boundary_conditions = 1;																						// boundary conditions - 0 - PBC, 1 - OBC, 2 - ABC,...
	std::string saving_dir = "";																								// directory for files to be saved onto

	// ----------------------------------- FUNCTIONS FOR READING THE INPUT

	string getCmdOption(const v_1d<string>& vec, string option) const;					 							// get the option from cmd input

	// ----------------------------------- TEMPLATES

	template <typename T>
	void set_option(T& value, const v_1d<string>& argv, string choosen_option, bool geq_0 = true);					// set an option

	template <typename T>
	void set_default_msg(T& value, string option, string message, const unordered_map <string, string>& map) const;	// setting value to default and sending a message

public:
	virtual ~user_interface() = default;

	// ----------------------------------- HELPING FUNCIONS

	virtual void set_default() = 0;																					// set default parameters
	virtual void exit_with_help() const = 0;
	virtual void printAllOptions() const = 0;
	// ----------------------------------- REAL PARSING

	virtual void parseModel(int argc, std::vector<std::string> argv) = 0;											// the function to parse the command line

	// ----------------------------------- NON-VIRTUALS

	std::vector<std::string> parseInputFile(string filename) const;													// if the input is taken from file we need to make it look the same way as the command line does

	// ----------------------------------- SIMULATIONS
	virtual void make_sim() = 0;
};
// ---------------------------------------------------------------------- ISING MODEL TRANSVERSE UI ----------------------------------------------------------------------
namespace isingUI
{
	// ----------------------------------- MAP FOR STRINGS -----------------------------------
	/// <summary>
	/// The map is used to parse also the two letters
	/// cases and create string variable from enum
	/// </summary>
	std::unordered_map <std::string, std::string> const table{
		{"f",""},						// file to read from directory
		{"J","1.0"},					// spin coupling
		{"Jn","1"},						// spin coupling step
		{"Js","0.0"},					// spin coupling num of points
		{"J0","0.0"},					// spin coupling randomness maximum (-J0 to J0)
		{"J0n","1"},					// spin coupling disorder step
		{"J0s","0.0"},					// spin coupling disorder num of points
		{"h","0.1"},					// perpendicular magnetic field constant
		{"hn","1"},						// longitudal magnetic field sweep number
		{"hs","0.1"},					// longitudal magnetic field sweep step
		{"w","0.01"},					// disorder strength
		{"wn","1"},						// longitudal disorder change number
		{"ws","0.1"},					// longitudal disorder change step
		{"g","1.0"},					// transverse magnetic field constant
		{"gn","1"},						// parameter scaling g number
		{"gs","0.1"},					// parameter scaling g step
		{"g0","0.0"},					// transverse field randomness maximum (-g0 to g0)
		{"g0n","1"},					// transverse disorder change number
		{"g0s","0.0"},					// transverse disorder change step
		{"L","4"},						// chain length
		{"Ln","1"},						// chain length step in size scaling
		{"Ls","1"},						// number of chain lengths in size scaling
		{"k","1"},						// translation symetry sector
		{"p","1"},						// parity symmetry sector
		{"x","1"},						// spin flip symmetry sector
		{"b","0"},						// boundary condition
		{"m","0"},						// choose model
		{"r","100"},					// realisations
		{"mu","5"},						// small bucket for the operator fluctuations to be averaged onto
		{"s","0"},						// site for operator averages
		{"p","0"},						// use parity symmetry?
		{"th","1"},						// number of threads
		{"op","0"},						// choose operator
		{"fun","-1"},					// choose function 
		{"ch", "0"},					// some boolean choose flag
		{"ts", "0.1"},					// time step for 
		{"scale", "0"},					// scale: linear-0 or log-1
		{"seed", "87178291199L"},		// seed foir random generator
		{"jobid", "0"},					// unique job id
		{"dim", "1"}					// dimensionality of anderson model
	};

	// ----------------------------------- UI CLASS SPECIALISATION -----------------------------------

	class ui : public user_interface {
	protected:
		// MODEL PARAMETERS
		double J, h, g;										// external fields
		int Jn, hn, gn, wn, J0n, g0n;						// external fields number of points
		double Js, hs, gs, ws, J0s, g0s;					// external fields step
		double w, g0, J0;									// disorder strengths
		int L, Ls, Ln, m;									// lattice params
		bool p, q, ch;										// boolean values
		int realisations;									// number of realisations to average on for disordered case - symmetries got 1
		size_t seed;										// radnom seed for random generator
		int jobid;											// unique _id given to current job

		int mu;												// small bucket for the operator fluctuations to be averaged onto
		int site;											// site for operator averages
		int op;												// choose operator
		int fun;											// choose function to start calculations
		double dt;											// time step for evolution
		int scale;											// choose scale: either linear or log
		int dim;											// dimension of lattice (mostly for anderson model)
		struct {
			int k_sym;										// translational symmetry generator
			int p_sym;										// parity symmetry generator
			int x_sym;										// spin-flip symmetry generator
		} symmetries;
	public:
		// ----------------------------------- CONSTRUCTORS
		ui() = default;
		ui(int argc, char** argv);														// standard constructor
		// ----------------------------------- PARSER FUNCTION FOR HELP
		void exit_with_help() const override;
		// ----------------------------------- HELPING FUNCIONS
		void set_default() override;													// set default parameters
		void printAllOptions() const override;
		// ----------------------------------- REAL PARSER
		void parseModel(int argc, std::vector<std::string> argv) override;				// the function to parse the command line

		// ----------------------------------- SIMULATION
		void make_sim() override;														// make default simulation
		



		enum class Ising_params{ L, J, J0, g, g0, h, w, k, p, x }; //<! choose which of these parameters
		auto get_params_name(Ising_params par){
			std::string name;
			switch(par){
				case Ising_params::L : 	name = "L"; 	break;
				case Ising_params::J : 	name = "J"; 	break;
				case Ising_params::J0 : name = "J0"; 	break;
				case Ising_params::g : 	name = "g"; 	break;
				case Ising_params::g0 : name = "g0"; 	break;
				case Ising_params::h : 	name = "h"; 	break;
				case Ising_params::w : 	name = "w"; 	break;
				case Ising_params::k : 	name = "k"; 	break;
				case Ising_params::p : 	name = "p"; 	break;
				case Ising_params::x : 	name = "x"; 	break;
				default:
					std::cout << "No option found, choosing g_array" << std::endl;
					name = "";
			}
			return name;
		}
		
		auto update_info(std::string old_info, const std::string& var, double new_value){
			old_info.erase(0,1);
			auto split_array = split_str(old_info, ",");
			std::string new_info = "_";
			for(auto& substr : split_array){
				auto str_tmp = split_str(substr, "=");
				if(str_tmp[0] == var){
					// found position of parameter with name var
					if(str_tmp[1].find(".") != std::string::npos) str_tmp[1] = to_string_prec(new_value, 2);	// print double value
					else str_tmp[1] = std::to_string(int(new_value));											// print int value
				}
				// reconstruct the info
				new_info += str_tmp[0] + "=" + str_tmp[1];
				if(&substr != &split_array.back()) new_info += ",";
			}
			return new_info;
		};


		//<! generate info by input parameters and model
		std::string generate_baseinfo(std::vector<std::string> skip = {}, std::string sep = "_"){
			if(this->m) return IsingModel_sym::set_info(this->L, this->J, this->g, this->h, this->symmetries.k_sym, this->symmetries.p_sym, this->symmetries.x_sym, skip, sep);
			else 		return IsingModel_disorder::set_info(this->L, this->J, this->J0, this->g, this->g0, this->h, this->w, skip, sep);
		}
		//-------------------------------------------------------------------------- GENERAL ROUTINES

		void diagonalize();
		//void diag_sparse(int num, bool get_eigenvectors = false, const char* form = "sa");	// diagonalize for limited number (set as num) of eigevals using form as choosing subset
		void diag_sparse(int num, bool get_eigenvectors = false, double sigma = 0.0);		// diagonalize for limited number (set as num) of eigevals starting at sigma and higher eigvals

		template <typename _type> [[nodiscard]]
		auto get_eigenvalues(IsingModel<_type>& alfa, std::string _suffix = "") -> arma::vec;

		//<! comvbine .hdf5 files seperated
		void combine_spectra();

		// --------------- COMPARISONS
		void compare_energies();
		void compare_matrix_elements(op_type op, int k_alfa, int k_beta, int p_alfa = 1, int p_beta = 1, int x_alfa = 1, int x_beta = 1);
		void compare_entaglement();
		void check_symmetry_rotation();

		void benchmark();
		//-------------------------------------------------------------------------- ANDERSON
		void calculate_localisation_length();

		//-------------------------------------------------------------------------- SPECTRAL PROPERTIES AND HELPERS
		//<! calculate all spectral quantities: time evolution, response function,
		//<! integrated spectral function and spectral form factor with folded eigenvalues
		void calculate_spectrals();

		//<! average separate spectral files for input operator
		void combine_spectrals();

		//<! calculate evolution of entaglement from initial state chosen by -op.
		//<! -s sets the subsystem size, if-s=0 the L/2 is assumed 
		void entropy_evolution();
		
		//<! analyze spectra with unfolding, DOS and level spacing distribution --  all to file
		void analyze_spectra();

		//-------------------------------------------------------------------------- ADIABATIC GAUGE POTENTIALS
		void adiabatic_gauge_potential();



		//-------------------------------------------------------------------------- STATISTICS
		//<! calculate all statistics for given input parameters -- averaged
		void calculate_statistics();
		void generate_statistic_map(Ising_params varname);
		void generate_statistic_map(Ising_params varname1, Ising_params varname2);

		//<! gap ratio map
		void level_spacing();

		//<! spectral form factor calculated from eigenvalues in file or diagonalize matrix
		void spectral_form_factor();
		
		//<! find thouless time with various method as function of h,g,J
		void thouless_times(Ising_params varname);

		void smoothen_data(const std::string& dir, const std::string& name, int mu = -1);

		//-------------------------------------------------------------------------- AVERAGE RAW DATA OVER REALISATIONS
		//<! average data over disorder realisations
		void average_SFF();



		//-------------------------------------------------------------------------- FUNCTIONS TO CALL IN FUN-DEFAULT MODE
		

		//-------------------------------------------------------------------------- INITIALIZERS
		//<! generate random product state (random orientation of spins on the bloch sphere)
		arma::cx_vec random_product_state(int system_size);

		//<! generate initial state given by input (user control) or -op flag (default): random, FM, AFM, ...
		arma::cx_vec set_init_state(size_t N, int choose = -1);

		//<! get parameter array from UI and input option
		arma::vec get_params_array(Ising_params par);

		//-------------------------------------------------------------------------- GENERAL LAMBDA'S
		//<! loop over all symmetry sectors
		template <
			typename callable, 
			typename... _types
			> 
			void loopSymmetrySectors(
			callable& lambda, //!< callable function
			_types... args										   //!< arguments passed to callable interface lambda
		) {
			const int x_max = (this->h != 0) ? 0 : 1;
		#pragma omp parallel for num_threads(outer_threads) schedule(dynamic)
			for (int k = 0; k < this->L; k++) {
				if (k == 0 || k == this->L / 2.) {
					for (int p = 0; p <= 1; p++){
						for (int x = 0; x <= x_max; x++){
							auto dummy_lambda = [&lambda](int k, int p, int x, auto... args){
								lambda(k, p, x, args...);
							};
							dummy_lambda(k, p, x, args...);
						}
					}
				}
				else {
					for (int x = 0; x <= x_max; x++){
							auto dummy_lambda = [&lambda](int k, int p, int x, auto... args){
								lambda(k, p, x, args...);
							};
							dummy_lambda(k, 0, x, args...);
					}
				}
			}
		}

		//<! loop over disorder realisations and call lambda each time
		template <
			Ising_params par, 		//<! parameter to average over in non-disordered case
			typename _ty, 			//<! type of model (double-> disorder, cpx-> sym)
			typename callable, 		//<! type of callable lambda
			typename... _types		//<! input types for lambda
		> 
		void average_over_realisations(
			IsingModel<_ty>& model,		//!< input model (symmetric model has to have average over external random stuff)
			bool with_diagonalization,	//!< checked if each realisation should diagonalize a new matrix
			callable& lambda, 			//!< callable function
			_types... args				//!< arguments passed to callable interface lambda
		) {
			double x = 0.0;
				switch (par)
				{
					case Ising_params::J: x = this->J;	break;
					case Ising_params::h: x = this->h;	break;
					case Ising_params::g: x = this->g;	break;
				default:				  x = 0.0;		break;
				}
			if(this->m){
				arma::vec _vec = x + my_gen.create_random_vec<double>(this->realisations, x / 50.);
				stout << _vec << std::endl;
				for(int r = 0; r < _vec.size(); r++){
					if(this->realisations > 1){
						switch (par)
						{
							case Ising_params::J: model.J = _vec(r);	break;
							case Ising_params::h: model.h = _vec(r);	break;
							case Ising_params::g: model.g = _vec(r);	break;
						default: 
							std::cout << "No default mode, average only performed over J, g, h" << std::endl;	
							break;
						}
					}
					if (with_diagonalization) {
						model.hamiltonian();
						model.diagonalization();
					}
					lambda(model, r, args...);
				}
			} else {
				for (int r = 0; r < this->realisations; r++) {
					if (with_diagonalization) {
						model.hamiltonian();
						model.diagonalization();
					}
					lambda(model, r, args...);
				}
			}
		};

		template <
			Ising_params par,	//<! which parameter to average over for symmetric case
			typename callable,	//<! callable lambda function
			typename... _types	//<! argument-types passed to lambda
		 > 
		void average_over_realisations(
			bool some_switch_I_add_later_cause_fak_compiler,
			callable& lambda, 			//!< callable function
			_types... args				//!< arguments passed to callable interface lambda
		) {
				double x = 0.0;
				switch (par)
				{
					case Ising_params::J: x = this->J;	break;
					case Ising_params::h: x = this->h;	break;
					case Ising_params::g: x = this->g;	break;
				default:				  x = 0.0;		break;
				}
			if(this->m){
				arma::vec _vec = x + my_gen.create_random_vec<double>(this->realisations, x / 50.);
				stout << _vec << std::endl;
			#pragma omp parallel for num_threads(outer_threads) schedule(dynamic)
				for(int r = 0; r < _vec.size(); r++){
					auto dummy_lambda = [&lambda](int real, double x, auto... args){
						lambda(real, x, args...);
					};
					dummy_lambda(r, _vec(r), args...);
				}
			} else{
			#pragma omp parallel for num_threads(outer_threads) schedule(dynamic)
				for (int r = 0; r < this->realisations; r++) {
					auto dummy_lambda = [&lambda](int real, double x, auto... args){
						lambda(real, x, args...);
					};
					dummy_lambda(r, x, args...);
				}
			}
		};
	};
}


#endif
