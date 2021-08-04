#pragma once
#ifndef UI
#define UI
#include "headers.h"

std::vector<std::string> change_input_to_vec_of_str(int argc, char** argv);

class user_interface{
protected:
	int thread_number;																						// number of threads
	int boundary_conditions;																				// boundary conditions - 0 - PBC, 1 - OBC, 2 - ABC,... 

	std::string saving_dir;
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
	enum class parsers {
		f,						// file to read from directory
		J,						// spin coupling
		h,						// perpendicular magnetic field constant
		g,						// transverse magnetic field constant
		w,						// disorder strength
		L,						// chain length
		b,						// boundary condition
		m,						// choose model
		p,						// use parity symmetry?
		th,						// number of threads
		q						// quit with help
	};
	/// <summary>
	/// The map is used to parse also the two letters
	/// cases and create string variable from enum
	/// </summary>
	std::unordered_map <std::string, parsers> const table{
		{"f",parsers::f},
		{"J",parsers::J},
		{"h",parsers::h},
		{"g",parsers::g},
		{"w",parsers::w},
		{"L",parsers::L},
		{"b",parsers::b},
		{"m",parsers::m},
		{"p",parsers::p},
		{"th",parsers::th},
		{"q",parsers::q}
	};

	class ui: public user_interface{
	protected:
		/* MODEL PARAMETERS */
		std::vector<double> J;
		double h,g,w;
		int L,b,m;
		bool p;
	public:
	/* CONSTRUCTORS */
		ui() = default;
	/* PARSER FUNCTION FOR HELP */
		void exit_with_help() override;
	/* REAL PARSER */
		void parseModel(int argc, std::vector<std::string> argv) override;							// the function to parse the command line
	/* HELPING FUNCIONS */
		void set_default() override;																// set default parameters
		void set_J(int mode, std::string argument);													// setting J according to the given mode
		void correct_J_size(int mode);																// if L is inconsistent with J size, we must add or erease elements according to L
	};
}







#endif

