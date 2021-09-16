#pragma once
#ifndef UI
#define UI
//#include "headers.h"
#include "IsingModel.h"

std::vector<std::string> change_input_to_vec_of_str(int argc, char** argv);

class user_interface {
protected:
	int thread_number;																						// number of threads
	int boundary_conditions;																				// boundary conditions - 0 - PBC, 1 - OBC, 2 - ABC,...
	std::string saving_dir;																					// directory for files to be saved onto

	std::string getCmdOption(const v_1d<std::string>& vec, std::string option) const;					 	// get the option from cmd input

	template <typename T>
	void set_option(T& value, const v_1d<std::string>& argv, std::string choosen_option, bool geq_0 = true);	// set an option

	template <typename T>
	void set_default_msg(T& value, std::string option, std::string message, \
		const std::unordered_map <std::string, std::string>& map) const;									// setting value to default and sending a message
	// std::unique_ptr<LatticeModel> model;															 			// a unique pointer to the model used

public:
	virtual ~user_interface() = default;
	/* HELPING FUNCIONS */
	virtual void set_default() = 0;																			// set default parameters

	virtual void exit_with_help() const = 0;
	/* REAL PARSING */
	virtual void parseModel(int argc, std::vector<std::string> argv) = 0;									// the function to parse the command line

/* NON-VIRTUALS */
	std::vector<std::string> parseInputFile(std::string filename) const;									// if the input is taken from file we need to make it look the same way as the command line does

// SIMULATIONS
	virtual void make_sim() = 0;
};

namespace isingUI
{
	/* MAP FOR STRINGS */
	/// <summary>
	/// The map is used to parse also the two letters
	/// cases and create string variable from enum
	/// </summary>
	std::unordered_map <std::string, std::string> const table{
		{"f",""},						// file to read from directory
		{"J","1.0"},					// spin coupling
		{"J0","0.2"},					// spin coupling randomness maximum (-J0 to J0)
		{"h","0.0"},					// perpendicular magnetic field constant
		{"hn","1"},						// longitudal magnetic field sweep number
		{"hs","0.0"},					// longitudal magnetic field sweep step
		{"w","1.0"},					// disorder strength
		{"wn","1"},						// disorder change number
		{"ws","0.0"},					// disorder change step
		{"g","1.0"},					// transverse magnetic field constant
		{"gn","1"},						// parameter scaling g number
		{"gs","0.0"},					// parameter scaling g step
		{"g0","0.0"},					// transverse field randomness maximum (-g0 to g0)
		{"L","4"},						// chain length
		{"Ln","1"},						// chain length step in size scaling
		{"Ls","0"},						// number of chain lengths in size scaling
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
	};

	class ui : public user_interface {
	protected:
		// MODEL PARAMETERS
		double J, h, g;																	// external fields
		int hn, gn, wn;																	// external fields number of points
		double hs, gs, ws;																// external fields step
		double w, g0, J0;																// disorder strengths
		int L, Ls, Ln, m;																// lattice params
		bool p, q;																		//
		int realisations;																// number of realisations to average on for disordered case - symmetries got 1
		int mu;																			// small bucket for the operator fluctuations to be averaged onto
		int site;																		// site for operator averages

		struct {
			double k_sym;																// translational symmetry generator
			int p_sym;																	// parity symmetry generator
			int x_sym;																	// spin-flip symmetry generator
		} symmetries;
		//std::unique_ptr<IsingModel<T>> model;											// pointer to a model
	public:
		// CONSTRUCTORS
		ui() = default;
		ui(int argc, char** argv);														// standard constructor
	// PARSER FUNCTION FOR HELP
		void exit_with_help() const override;
		// HELPING FUNCIONS
		void set_default() override;													// set default parameters
	// REAL PARSER
		void parseModel(int argc, std::vector<std::string> argv) override;				// the function to parse the command line

	// SIMULATION
		void make_sim() override;														// make default simulation
		void compare_energies();
		void disorder();
		
		void compare_matrix_elements();
		void size_scaling_sym(int k, int p, int x);

		void fidelity(std::initializer_list<int> symetries);

		void parameter_sweep_sym(int k, int p, int x);
		void check_dist_other_sector();
		void matrix_elements_stat_sym(double min, double max, double step, double omega_dist,\
			int omega_gauss_max, double energy_constraint, int energy_num,\
			std::initializer_list<int> alfa_sym = {},\
			std::initializer_list<int> beta_sym = {}) const;
		void perturbative_stat_sym(double dist_step, double min, double max, double pert, double gx, double hx);
		void perturbative_stat_sym(double dist_step, double min, double max, double pert, IsingModel_sym& alfa, IsingModel_sym& beta);
	};
}

#endif
