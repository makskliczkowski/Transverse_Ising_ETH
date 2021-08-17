#pragma once
#ifndef UI
#define UI
#include "headers.h"


std::vector<std::string> change_input_to_vec_of_str(int argc, char** argv);

class user_interface{
protected:
	int thread_number;																						// number of threads
	int boundary_conditions;																				// boundary conditions - 0 - PBC, 1 - OBC, 2 - ABC,... 

	std::string saving_dir;																					// directory for files to be saved onto
	std::string getCmdOption(const v_1d<std::string>& vec, std::string option) const;					 	// get the option from cmd input
	
	template <typename T>
	void set_option(T& value,const v_1d<std::string>& argv, std::string choosen_option, bool geq_0 = true);	// set an option

	template <typename T>
	void set_default_msg(T& value,std::string option, std::string message,\
		const std::unordered_map <std::string, std::string>& map) const;									// setting value to default and sending a message							 		
	// std::unique_ptr<LatticeModel> model;															 			// a unique pointer to the model used

public:
	virtual void exit_with_help() = 0;
/* REAL PARSING */
	virtual void parseModel(int argc, std::vector<std::string> argv) = 0;									// the function to parse the command line
/* HELPING FUNCIONS */
	virtual void set_default() = 0;																			// set default parameters
/* NON-VIRTUALS */
	std::vector<std::string> parseInputFile(std::string filename);											// if the input is taken from file we need to make it look the same way as the command line does
};




namespace isingUI
{
	/* MAP FOR STRINGS */
	/// <summary>
	/// The map is used to parse also the two letters
	/// cases and create string variable from enum
	/// </summary>
	std::unordered_map <std::string, std::string> const table{
		{"f",""},								// file to read from directory
		{"J","1.0"},					// spin coupling
		{"h","0.0"},					// perpendicular magnetic field constant
		{"g","0.0"},					// transverse magnetic field constant
		{"w","0.0"},					// disorder strength
		{"L","4"},						// chain length
		{"b","0"},						// boundary condition
		{"m","0"},						// choose model
		{"p","0"},						// use parity symmetry?
		{"th","1"},						// number of threads
		{"q","0"}						// quit with help
	};

	class ui: public user_interface{
	protected:
	// MODEL PARAMETERS
		std::vector<double> J;
		double h,g,w;
		int L,b,m;
		bool p,q;
	public:
	// CONSTRUCTORS 
		ui() = default;
		ui(int argc, char** argv);																	// standard constructor
	// PARSER FUNCTION FOR HELP 
		void exit_with_help() override;
	// REAL PARSER
		void parseModel(int argc, std::vector<std::string> argv) override;							// the function to parse the command line
	// HELPING FUNCIONS 
		void set_default() override;																// set default parameters
		void set_J(int mode, std::string argument);													// setting J according to the given mode
		void correct_J_size(int mode);																// if L is inconsistent with J size, we must add or erease elements according to L
	};
}







#endif

