#pragma once
#include <concepts> // has included type_traits, but is available for C++20 or newer

template <class T>
concept _has_output_operator
	= requires(std::ostream & os, T a)
	{ os << a; };

template <class T>
concept _has_input_operator
	= requires(std::istream & os, T a)
	{ os >> a; };

template <class T>
concept _has_IO
	= has_output_operator<T>
	&& has_input_operator<T>;

template <class T>
concept _has_access_operator1
	= requires(T x)
	{ x(); };
	
template <class T>
concept _has_access_operator2
	= requires(T x)
	{ x[]; };

template <class T>
concept _has_access_operator
	= _has_access_operator1<T>
	|| _has_access_operator2<T>;

template <class T>
concept _has_addition
	= requires(T x)
	{ x + x; };

template <class T>
concept _has_substraction
	= requires(T x)
	{ x - x; };

template <class T>
concept _has_multiplication
	= requires(T x, T y)
	{ x * y; };

template <class T>
concept _has_division
	= requires(T x, T y)
	{ x / y; };

template <class T>
concept _has_arithmetic_operators
	= requires(T x, T y)
	{
		x + y; x - y;
		x * y; x / y;
	};

template <class T>
concept _has_comparison
	= requires(T x, T y)
	{
		x > y;  x < y;
		x >= y; x <= y;
	};

template <class T>
concept _has_equality
	= requires(T x, T y)
	{
		x == y;	
		x != y;
	};


// already exist in <concept>: STL standard
/*
	//<!-------------------------------------------------------
// concepts to check integral type; is integral for: bool, char, char#_t, int#_t,...
	template<class T> concept integral = std::is_integral<T>::value;
	template<class T> concept signed_integral = integral<T> && std::is_signed<T>::value;
	template<class T> concept unsigned_integral = integral<T> && !signed_integral<T>;


	//<!-------------------------------------------------------
	// traits for lambda expression passed to functions

	// requires class F to be called with arguments Args... (is invocable with those arguments)
	template< class F, class... Args >
	concept invocable =
		requires(F && f, Args&&... args)
	{
		std::invoke(
			std::forward<F>(f),
			std::forward<Args>(args)...
		);
	};

	// adds to the above
	template< class F, class... Args >
	concept regular_invocable
		= std::invocable<F, Args...>;
 */