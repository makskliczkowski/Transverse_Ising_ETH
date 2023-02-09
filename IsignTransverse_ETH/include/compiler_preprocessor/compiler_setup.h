#pragma once

//<! for generic lambda input: has to be a callable
#if defined(HAS_CXX20)
	#define _typename std::invocable
#else
	#define _typename typename
#endif

#if !defined(FUN_SIGNATURE)
	#if defined (__GNUG__)
		#define FUN_SIGNATURE  __PRETTY_FUNCTION__
	#elif defined (_MSC_VER)
		#define FUN_SIGNATURE  __FUNCSIG__ 
	#elif defined(__INTEL_COMPILER)
		#define FUN_SIGNATURE  __FUNCTION__
	#else 
		#define FUN_SIGNATURE  __func__
	#endif
#endif

//! ---------------------------------------------------- define builtin functions
#ifdef _MSC_VER
#include <intrin.h>
#include <nmmintrin.h>
#define __builtin_popcount __popcnt
#define __builtin_popcountll _mm_popcnt_u64
#endif


//! --------------------------------------------------------- STATIC CONDITION CHECK\
// convert to string literals
#define TO_STRING_IMPL(x) #x
#define TO_STRING(x) TO_STRING_IMPL(x)

#define LINE_STR TO_STRING(__LINE__)

#define static_check(condition, str_lit) static_assert(condition, __FILE__"(line=" LINE_STR "): " str_lit)

//! --------------------------------- ERROR MESSAGES
#define NOT_CONSTRUCTIBLE	"ERROR 1: input is not constructible and cannot be stored in std::vector"
#define NOT_ALLOWED_INPUT	"ERROR 2: not allowed input to function"
#define NOT_CONVERTIBLE		"ERROR 3: not convertible to generic type"
#define BAD_INHERITANCE		"ERROR 4: given class is not inheriting from the appropriate base or links to class with 'final' keyword"
#define INCOMPATIBLE_SIZE	"ERROR 5: size of input arrays does not match"

//! --------------------------------------------------------- EXTRA DEBUG SETUP
#if defined(EXTRA_DEBUG)
	#define DESTRUCTOR_CALL std::cout << FUN_SIGNATURE << "->\tdestructor called" << std::endl;
	#define CONSTRUCTOR_CALL std::cout << FUN_SIGNATURE << "->\tconstructor called" << std::endl;
#else 
	#define DESTRUCTOR_CALL 
	#define CONSTRUCTOR_CALL
#endif

