#pragma once
#ifndef UI
#define UI
//#include "headers.h"
#include "IsingModel.h"
#include "spectrals.hpp"
#include "thermodynamics.hpp"
#include "statistics.hpp"
#include "statistics_dist.hpp"

const arma::vec down = { 0, 1 };
const arma::vec up	 = { 1, 0 };
extern std::uniform_real_distribution<> theta;
extern std::uniform_real_distribution<> fi;
// can't be const cause operator() is for non-const only (on microsoft it can be on const)
using namespace std;
std::vector<std::string> change_input_to_vec_of_str(int argc, char** argv);
// ----------------------------------------------------------------------------- GENERAL CLASS -----------------------------------------------------------------------------
class user_interface {
protected:
	unsigned int thread_number = 1;																								// number of threads
	int boundary_conditions = 1;																						// boundary conditions - 0 - PBC, 1 - OBC, 2 - ABC,...
	string saving_dir = "";																								// directory for files to be saved onto

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
		{"J0","0.0"},					// spin coupling randomness maximum (-J0 to J0)
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
		{"scale", "0"}					// scale: linear-0 or log-1
	};

	// ----------------------------------- UI CLASS SPECIALISATION -----------------------------------

	class ui : public user_interface {
	protected:
		// MODEL PARAMETERS
		double J, h, g;																	// external fields
		int hn, gn, wn, g0n;															// external fields number of points
		double hs, gs, ws, g0s;															// external fields step
		double w, g0, J0;																// disorder strengths
		int L, Ls, Ln, m;																// lattice params
		bool p, q, ch;																	// boolean values
		int realisations;																// number of realisations to average on for disordered case - symmetries got 1
		int mu;																			// small bucket for the operator fluctuations to be averaged onto
		int site;																		// site for operator averages
		int op;																			// choose operator
		int fun;																		// choose function to start calculations
		double dt;																		// time step for evolution
		int scale;																		// choose scale: either linear or log
		struct {
			int k_sym;																	// translational symmetry generator
			int p_sym;																	// parity symmetry generator
			int x_sym;																	// spin-flip symmetry generator
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
		
		//-------------------------------------------------------------------------- GENERAL ROUTINES
		void diagonalize();
		//void diag_sparse(int num, bool get_eigenvectors = false, const char* form = "sa");	// diagonalize for limited number (set as num) of eigevals using form as choosing subset
		void diag_sparse(int num, bool get_eigenvectors = false, double sigma = 0.0);		// diagonalize for limited number (set as num) of eigevals starting at sigma and higher eigvals

		template <typename _type> [[nodiscard]]
		auto get_eigenvalues(IsingModel<_type>& alfa) -> arma::vec;

		// --------------- COMPARISONS
		void compare_energies();
		void compare_matrix_elements(op_type op, int k_alfa, int k_beta, int p_alfa = 1, int p_beta = 1, int x_alfa = 1, int x_beta = 1);
		void compare_entaglement();

		void benchmark();
		
		//-------------------------------------------------------------------------- SPECTRAL PROPERTIES AND HELPERS
		//<! calculate all spectral quantities: time evolution, response function,
		//<! integrated spectral function and spectral form factor with folded eigenvalues
		void calculate_spectrals();

		//<! calculate evolution of entaglement from initial state chosen by -op.
		//<! -s sets the subsystem size, if-s=0 the L/2 is assumed 
		void entropy_evolution();
		
		//<! loop over all parameters (L, site, g, h) for given disorder
		//<! or symmetry sector and find relaxation times as I(w)=1/2 (the later from integrated time evolution)
		void relaxationTimesFromFiles();
		void intSpecFun_from_timeEvol();
	
		//<! spectral form factor calculated from eigenvalues in file or diagonalize matrix
		void spectral_form_factor();
		
		
		//-------------------------------------------------------------------------- ADIABATIC GAUGE POTENTIALS
		template <typename _type> void LevelSpacingDist(IsingModel<_type>& alfa);
		void adiabaticGaugePotential_sym(bool SigmaZ = 0, bool avSymSectors = 0);
		void adiabaticGaugePotential_dis();
		void combineAGPfiles();



		//-------------------------------------------------------------------------- STATISTICS
		void level_spacing_from_distribution();
		void level_spacing();

		void smoothen_data(const std::string& dir, const std::string& name);

		//-------------------------------------------------------------------------- FUNCTIONS TO CALL IN FUN-DEFAULT MODE
		/// <summary>
		/// saves matrix elements for using in the autoencoder
		/// </summary>
		/// <param name="operators"> inializer list of different local operators and also their global friend </param>
		/// <param name="names"> names of operators set in input, must be equal size as operators</param>
		void saveDataForAutoEncoder_disorder(std::initializer_list<op_type> operators, std::initializer_list<std::string> names);
		void saveDataForAutoEncoder_symmetries(std::initializer_list<op_type> operators, std::initializer_list<std::string> names);


		//-------------------------------------------------------------------------- GENERAL LAMBDA'S
		//<! generate random product state (random orientation of spins on the bloch sphere)
		arma::cx_vec random_product_state(int system_size)
		{
			auto the = theta(gen);
			arma::cx_vec init_state = std::cos(the / 2.) * up
				+ std::exp(im * fi(gen)) * std::sin(the / 2.) * down;
			for (int j = 1; j < system_size; j++)
			{
				the = theta(gen);
				init_state = arma::kron(init_state, std::cos(the / 2.) * up
					+ std::exp(im * fi(gen)) * std::sin(the / 2.) * down);
			}
			return init_state;
		};

		//<! generate initial state given by -op flag: random, FM, AFM, ...
		arma::cx_vec set_init_state(size_t N)
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

		//<! loop over all symmetry sectors
		template <typename... _types> void loopSymmetrySectors(
			std::function<void(int,int,int,_types...args)> lambda, //!< callable function
			double hx,											   //!< longitudal field -- whether spin-flip symmetry is allowed
			int Lx,												   //!< system size
			_types... args										   //!< arguments passed to callable interface lambda
		) {
			const int x_max = (hx != 0) ? 0 : 1;
			for (int k = 0; k < Lx; k++) {
				if (k == 0 || k == Lx / 2.) {
					for (int p = 0; p <= 1; p++)
						for (int x = 0; x <= x_max; x++)
							lambda(k, p, x, std::forward<_types>(args)...);
				}
				else {
					for (int x = 0; x <= x_max; x++)
						lambda(k, 0, x, std::forward<_types>(args)...);
				}
			}
		}

		//<! loop over disorder realisations and call lambda each time
		template <typename _ty, typename... _types> 
		void average_over_realisations(
			IsingModel<_ty>& model,				   	   //!< input model (symmetric model has to have average over external random stuff)
			bool with_diagonalization,				   //!< checked if each realisation should diagonalize a new matrix
			std::function<void(_types...args)> lambda, //!< callable function
			_types... args							   //!< arguments passed to callable interface lambda
		) {
			model.reset_random();
			for (int r = 0; r < this->realisations; r++) {
				if (with_diagonalization) {
					model.hamiltonian();
					model.diagonalization();
				}
				lambda(std::forward<_types>(args)...);
			}
		};
	};
}


#endif
