#pragma once
#ifndef UI
#define UI
//#include "headers.h"
#include "IsingModel.h"

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
		{"ts", "0.1"}					// time step for evolution
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
		double ts;																		// time step for evolution
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
		// --------------- DISORDER
		void disorder();
		// --------------- COMPARISONS
		void compare_energies();
		void compare_matrix_elements(op_type op, int k_alfa, int k_beta, int p_alfa = 1, int p_beta = 1, int x_alfa = 1, int x_beta = 1);
		void compare_entaglement();

		// --------------- SYMMETRIES
		void size_scaling_sym(int k, int p, int x);

		void fidelity(std::initializer_list<int> symetries);

		void benchmark(bool full = true);
		//-------------------------------------------------------------------------- RARELY USED SWEEPS OF PARAMETERS
		void parameter_sweep_sym(int k, int p, int x);
		void check_dist_other_sector();
		void matrix_elements_stat_sym(double min, double max, double step, double omega_dist, \
			int omega_gauss_max, double energy_constraint, int energy_num, \
			std::initializer_list<int> alfa_sym = {}, \
			std::initializer_list<int> beta_sym = {}) const;


		//-------------------------------------------------------------------------- SPECIFIC FOR MODELS
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
		template <typename _ty, typename... _types> void average_over_realisations(
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
		//template <typename... _types>
		//-------------------------------------------------------------------------- SPECTRAL PROPERTIES
		/// <summary>
		/// 
		/// </summary>
		/// <param name="alfa"> input model: cpx when with symmetries </param>
		/// <param name="opMatrix"> input operator as sparse matrix </param>
		/// <param name="name"> name for file to store data </param>
		template <typename _type> void spectralFunction(IsingModel<_type>& alfa, const arma::cx_mat& mat_elem, std::string name);
		template <typename _type> void integratedSpectralFunction(const IsingModel<_type>& alfa, const arma::cx_mat& mat_elem, std::string name);
		template <typename _type> auto integratedSpectralFunction(const IsingModel<_type>& alfa, const arma::cx_mat& mat_elem, const arma::vec& omegas) -> arma::vec;

		template <typename _type> void timeEvolution(const IsingModel<_type>& alfa, const arma::cx_mat& mat_elem, std::string name = "STH");
		template <typename _type> auto timeEvolution(const IsingModel<_type>& alfa, const arma::cx_mat& mat_elem, const arma::vec& times) -> std::pair<arma::vec, double>;
		void entropy_evolution();

		void relaxationTimesFromFiles();
		void intSpecFun_from_timeEvol();
		
		template <typename _type> void IsingLIOMs(IsingModel<_type>& alfa);
		void TFIsingLIOMs();
		void LIOMsdisorder();

		template <typename _type> void LevelSpacingDist(IsingModel<_type>& alfa);
		template <typename _type> std::pair<double, double> operator_norm(arma::sp_cx_mat& opMatrix, IsingModel<_type>& alfa);
		void adiabaticGaugePotential_sym(bool SigmaZ = 0, bool avSymSectors = 0);
		void adiabaticGaugePotential_dis(bool h_vs_g = true);
		void combineAGPfiles();
		template <typename _type> void energyEvolution(IsingModel<_type>& model);

		/// <summary>
		/// saves matrix elements for using in the autoencoder
		/// </summary>
		/// <param name="operators"> inializer list of different local operators and also their global friend </param>
		/// <param name="names"> names of operators set in input, must be equal size as operators</param>
		void saveDataForAutoEncoder_disorder(std::initializer_list<op_type> operators, std::initializer_list<std::string> names);
		void saveDataForAutoEncoder_symmetries(std::initializer_list<op_type> operators, std::initializer_list<std::string> names);


		//------------------------------------------------------------------------------ PERTURBATIVE METHODS (DISTRIBUTIONS WITH PERTURBATIONS)
		std::vector<double> perturbative_stat_sym(double pert, double gx, double hx);
		std::vector<double> perturbative_stat_sym(double pert, IsingModel_sym& alfa, double gx, double hx);
		std::vector<double> perturbative_stat_sym(IsingModel_sym& alfa, double gx, double hx);
		std::vector<double> perturbative_stat_sym(double dist_step, double min, double max, double pert, IsingModel_sym& alfa, IsingModel_sym& beta);
		std::vector<double> perturbative_stat_sym(double pert) {
			return perturbative_stat_sym(pert, this->g, this->h);
		}
		std::vector<double> perturbative_stat_sym(double pert, IsingModel_sym& alfa) {
			return perturbative_stat_sym(pert, alfa, this->g, this->h);
		}
		std::vector<double> perturbative_stat_sym(IsingModel_sym& alfa) {
			return perturbative_stat_sym(alfa, this->g, this->h);
		}
	};
}

#endif
