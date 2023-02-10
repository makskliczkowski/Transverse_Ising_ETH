#pragma once
#include <algorithm>
#include <complex>
#include <vector>
#include <string>
#include "../miscaleneous/constants.hpp"

//! -------------------------------------------------------------------------------------------------------------- TYPEDEFs
#if defined(USE_LONG_LONG)
	typedef size_t u64; //TODO: find better name for that
#else
	typedef unsigned long u64;
#endif

#if defined(USE_LONG_DBL)
	typedef std::complex<long double> cpx;
#else
	typedef std::complex<double> cpx;
#endif

static const char* kPathSeparator =
#ifdef _WIN32
	"\\";
#else
	"/";
#endif

typedef std::vector<int> iVec;
typedef std::vector<short int> kernel_arr;
typedef std::vector<cpx> _cx_vec;

using op_type = std::function<std::pair<cpx, u64>(u64, int, std::vector<int>)>;
typedef std::chrono::system_clock clk;
const std::string kPSep = std::string(kPathSeparator);

template<class T> using v_3d = std::vector<std::vector<std::vector<T>>>;		// 3d type T vector
template<class T> using v_2d = std::vector<std::vector<T>>;						// 2d type T vector
template<class T> using v_1d = std::vector<T>;									// 1d type T vector
//<! -------------------------------------------------------------------------------------------------------------- ATTRIBUTES
#define _noreturn			[[noreturn]]
#define _maybe_unused		[[maybe_unused]]
#define _nodiscard			[[nodiscard, gnu::warn_unused_result]]
#define _no_break			[[fallthrough]]
#define _no_unique_address	[[no_unique_address]]


//<! -------------------------------------------------------------------------------------------------------------- MACROs

#define ignore(variable)  ((void)(variable))		// ignore/dicard variable

struct sink { template<typename ..._ty> sink(_ty const && ...) {} };
#define ignore_pack(...)  (sink { std::move(__VA_ARGS__)... })		// ignore/discard variadic input

#define _type_name(name) #name
#define _typeID(name) typeid(name)
#define var_name(name) std::string(_type_name(name));
#define var_name_value(name,prec) std::string(_type_name(name))+ "=" + to_string_prec((name), prec);

// ------- Macros to generate the lookup table (at compile-time)
#define generate_1st_sequence(n) n, n + 2*64, n + 1*64, n + 3*64
#define generate_2nd_sequence(n) generate_1st_sequence(n), generate_1st_sequence(n + 2*16), generate_1st_sequence(n + 1*16), generate_1st_sequence(n + 3*16)
#define generate_3rd_sequence(n) generate_2nd_sequence(n), generate_2nd_sequence(n + 2*4 ), generate_2nd_sequence(n + 1*4 ), generate_2nd_sequence(n + 3*4 )
#define REVERSE_BITS generate_3rd_sequence(0), generate_3rd_sequence(2), generate_3rd_sequence(1), generate_3rd_sequence(3)


//<! -------------------------------------------------------------------------------------------------------------- CONSTANTS
constexpr cpx im = cpx(0.0, 1.0);					// imaginary unit
#if defined(USE_LONG_DBL)
	const long double pi = constants<long double>::pi;
	const long double two_pi = constants<long double>::two_pi;
#else
	const double pi = constants<double>::pi;
	const double two_pi = constants<double>::two_pi;
#endif


//! -------------------------------------------------------------------------------------------------------------- USER MAKRO FUNCTIONS
#if !defined(ULLPOW)
	#define ULLPOW(k) (1ULL << k)
#endif

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

