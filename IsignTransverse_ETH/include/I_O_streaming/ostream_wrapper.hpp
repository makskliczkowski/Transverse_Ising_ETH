#pragma once

#ifndef my_IO
    #include "stream_wrapper_base.hpp"
#endif

namespace files{

//----------------------------------------------------------------------
//---------------------------------------------------------------------- OUTPUT STREAM CLASS
//  -- wrapper on std::ostream
//
//<! class for output to stream (write), inheritaed from FILE base
template <std::ios_base::openmode mode>
class oFILE : public FILE<std::ostream, mode> {
private:
	std::string separator = "\t";
	arma::u16 width		  = 10;

public:

	//-------------------------------------------------------------------------- CONSTRUCTORS
	~oFILE() = default;
	oFILE() = default;
	oFILE(
		std::string name	   = "",	//<! name of file used as outstream
		std::string _separator = "\t",	//<! separator between columns in file 
		arma::u16 _width	   = 10,	//<! width of column in file
		std::string _dir	   = ""		//<! directory to file (may be included in name)
	)
		: separator(_separator), width(_width)
	{
		this->filename = name;
		this->dir = _dir;
		this->open(this->filename, this->dir);
	};

	//-------------------------------------------------------------------------- OVERLOADED OPERATORS
	template <has_output_operator _ty> //<! using concept to ensure type has << operator
	friend oFILE& operator<<(oFILE& output, _ty arg) {
		output.file.width(output.width);
		output.file << arg << output.separator;
		return output;
	}

	//-------------------------------------------------------------------------- PRINTERS AND SAVE TO FILE
	//<! variadic printer
	template <has_output_operator... _ty> //<! using concept to ensure type has << operator
	void print(_ty... rest) 
		{ (this << ... << rest) << std::endl; }

	//<! save data to file
	template <has_access_operator _ty, has_access_operator... vec>
	void save(_ty first, vec... data) {
		//if (sizeof...(vec) > 0)
		//	assert([&first](const auto&... args) {
		//	return ((first.size() == args.size()) && ...);
		//		}(data) && "incompatible dimensions");
		for (int i = 0; i < first.size(); i++) {
			this << first;
			([&](const auto& x) {this << x(i); }(data), ...);
		}
	}

	//<! save data to file with additional variables at &x = begin() "x0=0" 
	template <typename... _ty>
	void save_with_x0_vars(
		const arma::vec& x,
		const arma::vec& y,
		_ty... args
	) {
		assert(x.size() == y.size() && "Incompatible dimensions");
		for (int i = 0; i < x.size(); i++) {
			if (i == 0) ((this << x(i) << y(i)) << ... << args) << std::endl;
			else		  this << x(i) << y(i) << std::endl;
		}
	}
};

};