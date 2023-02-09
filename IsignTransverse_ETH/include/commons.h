#pragma once

#include "compiler_preprocessor/preprocessor_setup.hpp"

//-----------------------
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>
#include <complex>
#include <cmath>
#include <algorithm>
#include <unordered_map>
//#include <hdf5.h>
//#include <mkl.h>
DISABLE_WARNING_PUSH // include <armadillo> and suppress its warnings, cause developers suck
	// armadillo flags:
#define ARMA_64BIT_WORD // enabling 64 integers in armadillo obbjects
#define ARMA_BLAS_LONG_LONG // using long long inside LAPACK call
#define ARMA_USE_OPENMP
#define ARMA_ALLOW_FAKE_GCC
//#define ARMA_USE_SUPERLU
//#define ARMA_USE_HDF5
//#define ARMA_EXTRA_DEBUG
//-------
DISABLE_OVERFLOW;
#if defined(_MSC_VER)
	DISABLE_WARNING(26812); // unscoped enum
	DISABLE_WARNING(26819); // unannotated fallthrough
	DISABLE_WARNING(26439); // may not throw
	DISABLE_WARNING(6011);  // dereferencing NULL ptr 
	DISABLE_WARNING(26495); // unitialized variable
	DISABLE_WARNING(6993);  // ignore OpenMP: use single-thread
	DISABLE_WARNING(4849);  // ignor OpenMP:collapse
#elif defined(__GNUC__) || defined(__clang__)
	DISABLE_WARNING(-Wenum-compare); // unscoped enum
	DISABLE_WARNING(-Wimplicit-fallthrough); // unannotated fallthrough
	DISABLE_WARNING(-Wuninitialized); // unitialized
	DISABLE_WARNING(-Wopenmp);  // ignore OpenMP warning
#else 
	#pragma message ("not recognized compiler to disable armadillo library warnings");
#endif
	#include <armadillo>
	
DISABLE_WARNING_POP

#include <iterator>
#include <cassert> // assert terminates program
#include <omp.h>
#include <ctime>
#include <utility> // auto, etc.
#include <memory> // smart ptr
#include <thread>
#include <future>
//#include <condition_variable>
#include <functional>
#include <concepts>
//#include <execution>
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

#include "compiler_preprocessor/compiler_setup.h"


typedef size_t u64;

// ----------------------------------------------------------------------------- definitions -----------------------------------------------------------------------------
static const char* kPathSeparator =
#ifdef _WIN32
"\\";
#else
"/";
#endif

template<class T>
using v_3d = std::vector<std::vector<std::vector<T>>>;											// 3d double vector
template<class T>
using v_2d = std::vector<std::vector<T>>;														// 2d double vector
template<class T>
using v_1d = std::vector<T>;																	// 1d double vector


#define im cpx(0.0,1.0)
#define stout std::cout << std::setprecision(8) << std::fixed								// standard outstream

// ----------------------------------------------------------------------------- Macros to generate the lookup table (at compile-time) -----------------------------------------------------------------------------
#define R2(n) n, n + 2*64, n + 1*64, n + 3*64
#define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
#define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )
#define REVERSE_BITS R6(0), R6(2), R6(1), R6(3)
#define ULLPOW(k) (1ULL << k)
#define RETURNS(...) -> decltype((__VA_ARGS__)) { return (__VA_ARGS__); }
// ----------------------------------------------------------------------------- lookup table to store the reverse of each index of the table -----------------------------------------------------------------------------
// The macro `REVERSE_BITS` generates the table
const u64 lookup[256] = { REVERSE_BITS };

const v_1d<u64> BinaryPowers = { ULLPOW(0), ULLPOW(1), ULLPOW(2), ULLPOW(3),
								ULLPOW(4), ULLPOW(5), ULLPOW(6), ULLPOW(7),
								ULLPOW(8), ULLPOW(9), ULLPOW(10), ULLPOW(11),
								ULLPOW(12), ULLPOW(13), ULLPOW(14), ULLPOW(15),
								ULLPOW(16), ULLPOW(17), ULLPOW(18), ULLPOW(19),
								ULLPOW(20), ULLPOW(21), ULLPOW(22), ULLPOW(23),
								ULLPOW(24), ULLPOW(25), ULLPOW(26), ULLPOW(27),
								ULLPOW(28), ULLPOW(29), ULLPOW(30), ULLPOW(31),
								ULLPOW(32), ULLPOW(33), ULLPOW(34), ULLPOW(35),
								ULLPOW(36), ULLPOW(37), ULLPOW(38), ULLPOW(39),
								ULLPOW(40), ULLPOW(41), ULLPOW(42), ULLPOW(43),
								ULLPOW(44), ULLPOW(45), ULLPOW(46), ULLPOW(47),
								ULLPOW(48), ULLPOW(49), ULLPOW(50), ULLPOW(51),
								ULLPOW(52), ULLPOW(53), ULLPOW(54), ULLPOW(55),
								ULLPOW(56), ULLPOW(57), ULLPOW(58), ULLPOW(59),
								ULLPOW(60), ULLPOW(61), ULLPOW(62), ULLPOW(63) }; // vector containing powers of 2 from 2^0 to 2^(L-1)


using cpx = std::complex<double>;
using op_type = std::function<std::pair<cpx, u64>(u64, int, std::vector<int>)>;
using clk = std::chrono::system_clock;


//! ------------------------------------------------------ ATTRIBUTES
#define _noreturn			[[noreturn]]
#define _maybe_unused		[[maybe_unused]]
#define _nodiscard			[[nodiscard, gnu::warn_unused_result]]
#define _no_break			[[fallthrough]]
#define _no_unique_address	[[no_unique_address]]

// ----------------------------------------------------------------------------- CONSTANTS -----------------------------------------------------------------------------

constexpr long double pi = 3.141592653589793238462643383279502884L;			// it is me, pi
constexpr long double two_pi = 2 * 3.141592653589793238462643383279502884L;	// it is me, 2pi
//const auto global_seed = std::random_device{}();							// global seed for classes
const std::string kPSep = std::string(kPathSeparator);


// ----------------------------------------------------------------------------- TRY-Catch -----------------------------------------------------------------------------
inline void handle_exception(std::exception_ptr eptr, std::string message) {
	try {
		if (eptr) {
			std::rethrow_exception(eptr);
		}
	}
	catch (const std::runtime_error& err) {
		stout << "Runtime error:\t" << err.what() << "\n";
		stout << message << std::endl;
		assert(false);
	}
	catch (const std::bad_alloc& err) {
		stout << "Bad alloc error:\t" << err.what() << "\n";
		stout << message << std::endl;
		assert(false);
	}
	catch (const std::exception& err) {
		stout << "Exception:\t" << err.what() << "\n";
		stout << message << std::endl;
		assert(false);
	}
	catch (...) {
		stout << "Unknown error...!" << "\n";
		stout << message << std::endl;
		assert(false);
	}
}
#define BEGIN_CATCH_HANDLER try{
#define END_CATCH_HANDLER(message) } catch(...){ handle_exception(std::current_exception(), message); };

// ----------------------------------------------------------------------------- TIME FUNCTIONS -----------------------------------------------------------------------------

inline double tim_s(clk::time_point start) {
	return double(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::duration(\
		std::chrono::system_clock::now() - start)).count()) / 1000.0;
}
// ----------------------------------------------------------------------------- TOOLS -----------------------------------------------------------------------------
enum class coordinate {
	x, 
	y,
	z
};

template <typename _baseTy>
struct types {
	template<typename... _otherTy>
	using convertible = typename std::enable_if<std::conjunction<std::is_convertible<_otherTy, _baseTy>...>::value>::type;
};
namespace _traits {
	//<! check if variadic templates have common types (or are inmplicitly convertible) with _baseTy
	template <class _baseTy, class... _otherTy>
	struct is_same : std::conjunction<std::is_same<_otherTy, _baseTy>...> {};

	template <class _baseTy, class... _otherTy>
	inline constexpr bool is_same_v = is_same<_baseTy, _otherTy...>::value;

	//<! check if typename is among other fixed ones struct (STL only has constexpr)
	template <class _ty, class... fixed>
	struct is_any_of : std::disjunction<std::is_same<_ty, fixed>...> {};

	template <class _ty, class... fixed>
	inline constexpr bool is_any_of_v = is_any_of<_ty, fixed...>::value;
};


template <typename _ty>
inline
std::string type_name(_ty a){
	if(typeid(_ty) == typeid(cpx))			
		return "COMPLEX";
	else if(typeid(_ty) == typeid(double))	
		return "DOUBLE";
	else if(typeid(_ty) == typeid(float))	
		return "FLOAT";
	else 
		return "UNKNOWN";
}