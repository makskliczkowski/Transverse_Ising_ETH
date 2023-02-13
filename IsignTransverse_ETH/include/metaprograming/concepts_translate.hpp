#pragma once

#if defined(HAS_CXX20) && !defined( __clang__ )
    #include "concepts.hpp"
    
	#define has_output_operator         _has_output_operator
    #define has_input_operator          _has_input_operator
    #define has_IO                      _has_IO

    #define has_access_operator         _has_access_operator
    #define has_addition                _has_addition
    #define has_substraction            _has_substraction
    #define has_multiplication          _has_multiplication
    #define has_division                _has_division
    #define has_arithmetic_operators    _has_arithmetic_operators
    #define has_comparison              _has_comparison
    #define has_equality                _has_equality
#else
    #if defined( __clang__ )
        #pragma message ("--> Clang compiler does not recognize concept! Thus is depracated.")
    #endif

	#define has_output_operator         typename
    #define has_input_operator          typename
    #define has_IO                      typename

    #define has_access_operator         typename
    #define has_addition                typename
    #define has_substraction            typename
    #define has_multiplication          typename
    #define has_division                typename
    #define has_arithmetic_operators    typename
    #define has_comparison              typename
    #define has_equality                typename

#endif
