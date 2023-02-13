#pragma once

#ifndef my_IO
    #include "stream_wrapper_base.hpp"
#endif

namespace files{

//----------------------------------------------------------------------
//---------------------------------------------------------------------- INPUT STREAM CLASS
//  -- wrapper on std::istream
// 
//<! class for input from stream (read), inheritaed from FILE base
template <std::ios_base::openmode mode>
class iFILE : public FILE<std::istream, mode> {
private:
	int rows = 0;
	int cols = 0;

	void get_data_dimensions() {
		std::string datarow;
		while (std::getline(this->file, datarow)) {
			std::istringstream ss(datarow);
			if (rows == 0) {
				double value;
				while (ss >> value)	cols++;
			}
			rows++;
		}
		this->file.close();
	}
public:

	//-------------------------------------------------------------------------- CONSTRUCTORS
	~iFILE() = default;
	iFILE() = default;
	iFILE(
		std::string name = "",	//<! name of file used as outstream
		std::string _dir = ""	//<! directory to file (may be included in name)
	) {
		this->filename = name;
		this->dir = _dir;
		this->open(this->filename, this->dir);
	};

	//-------------------------------------------------------------------------- OVERLOADED OPERATORS
	template <has_input_operator _ty>
	friend iFILE& operator<<(iFILE& input, _ty& var) {
		input.file >> var;
		return input;
	}

	//-------------------------------------------------------------------------- READERS (variadic and to dataset)
	template <has_input_operator first, has_input_operator... _ty>
	void variadic_read(first arg1, _ty&... var) {
		if (this->file.peek() == '\n' || this->file.peek() == EOF) {
			std::cout << "Less data in file row than in function input.\nLast" << sizeof...(_ty) + 1 << " variables left unchanged!";
			return;
		}
		this->file >> arg1;
		variadic_read(var...);
	}
	auto read() {
		std::vector<arma::vec> data;
		if (!this->isOpen()) {
			std::cout << "File not open! could not load";
			assert(false);
		}
		get_data_dimensions();
		this->reOpen();
		// set dataset dimensions
		data = std::vector<arma::vec>(cols);
		for (int i = 0; i < cols; i++)
			data[i] = arma::vec(rows, arma::fill::zeros);

		// read from file
		int j = 0;
		std::string datarow;
		while (std::getline(this->file, datarow)) {
			std::istringstream ss(datarow);
			int i = -1; 
			double value;
			while (ss >> value)
				data[++i](j) = value;
			j++;
		}
		return data;
	}
};

};