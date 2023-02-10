#pragma once
#ifndef my_IO
#define my_IO

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#ifndef ARMA
	#include "../armadillo_wrapper.hpp"
#endif

#ifdef __has_include
#  if __has_include(<filesystem>)
#    include <filesystem>
#    define have_filesystem 1
namespace fs = std::filesystem;
#  elif __has_include(<experimental/filesystem>)
#    include <experimental/filesystem>
#    define have_filesystem 1
#    define experimental_filesystem
namespace fs = std::experimental::filesystem;
#  else
#    define have_filesystem 0
#  endif
#endif

inline void createDirs(const std::string& dir) {
	if (!fs::is_directory(dir) || !fs::exists(dir))
		fs::create_directories(dir);
}
template <typename... _Ty>
inline void createDirs(const std::string& dir, const _Ty&... dirs) {
	createDirs(dir);
	createDirs(dirs...);
}

//<! base class to input/output filestream
namespace files {

	template<
		typename _type,							//<! typename of 'file' in class
		std::ios_base::openmode mode			//<! openmode: write/read/append...									
	>
		class FILE {

		//<! only valid if template _type is inheriting from ios_base (base stream class)
		static_check((std::is_base_of_v<std::ios_base, _type>), BAD_INHERITANCE ":: base class is: std::ios_base");

		protected:
			_type file;					//<! instance of file (fstream)
			std::string dir = "";		//<! directory to given file
			std::string filename = "";	//<! name of file (may include directory)

            /// @brief 
            /// @tparam _type 
            void open(std::string filename, std::string dir = "") 
            {
				if (!dir.empty())
					createDirs(dir);
				if (!filename.empty())
					file.open(dir + filename, mode);
			}
			void reOpen() 
                { this->file.open(dir + filename, mode); }
		public:

			//-------------------------------------------------------------------------- CONSTRUCTORS
			FILE() = default;
			virtual ~FILE() 
                { this->file.close(); }

			//-------------------------------------------------------------------------- OVERLOADED OPERATORS
			
            /// @brief 
            /// @tparam _type 
            bool operator()(std::string filename, std::string dir = "") 
            {
				open(filename, dir);
				return this->file.is_open();
			}


			//-------------------------------------------------------------------------- OTHER FUNCTIONS
			
            /// @brief 
            /// @return 
            bool isOpen() 
                { return this->file.is_open(); }

	};
};

#include "istream_wrapper.hpp"
#include "ostream_wrapper.hpp"

#endif
