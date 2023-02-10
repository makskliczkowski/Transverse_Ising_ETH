#pragma once
#ifndef __COMMONS
#define __COMMONS

//-----------------------
#include <cmath>
#include <omp.h>

#include <cassert> // assert terminates program
#include <ctime>
#include <utility> // auto, etc.
#include <memory> // smart ptr
#include <thread>
#include <future>
#include <functional>

#include "compiler_setup_headers/preprocessor_setup.hpp"
#include "compiler_setup_headers/compiler.hpp"
#include "metaprograming/traits.hpp"

#include "miscaleneous/tools.hpp"

#include "armadillo_wrapper.hpp"

#include "I_O_streaming/stream_wrapper_base.hpp"





inline void handle_exception(std::exception_ptr eptr, std::string message) {
	try {
		if (eptr) {
			std::rethrow_exception(eptr);
		}
	}
	catch (const std::runtime_error& err) {
		std::cout << "Runtime error:\t" << err.what() << "\n";
		std::cout << message << std::endl;
		assert(false);
	}
	catch (const std::bad_alloc& err) {
		std::cout << "Bad alloc error:\t" << err.what() << "\n";
		std::cout << message << std::endl;
		assert(false);
	}
	catch (const std::exception& err) {
		std::cout << "Exception:\t" << err.what() << "\n";
		std::cout << message << std::endl;
		assert(false);
	}
	catch (...) {
		std::cout << "Unknown error...!" << "\n";
		std::cout << message << std::endl;
		assert(false);
	}
}
#define BEGIN_CATCH_HANDLER try{
#define END_CATCH_HANDLER(message) } catch(...){ handle_exception(std::current_exception(), message); };

#endif