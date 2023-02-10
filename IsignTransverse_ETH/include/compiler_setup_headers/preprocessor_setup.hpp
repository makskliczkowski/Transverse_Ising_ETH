#pragma once

//! --------------------------------- SUPPRESS WARNINGS -- portable

#if defined(_MSC_VER) || defined(__INTEL_COMPILER)
// intel also defines _MSC_VEC
#define __PRAGMA(X) __pragma(X)			//<! portable pragma definition to use inside preprocessor commands
#define DISABLE_WARNING_PUSH           __PRAGMA(warning( push ))
#define DISABLE_WARNING_POP            __PRAGMA(warning( pop )) 
#define DISABLE_WARNING(warningNumber) __PRAGMA(warning( disable : warningNumber ))

#define DISABLE_OVERFLOW						 DISABLE_WARNING(26451)
#define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER    DISABLE_WARNING(4100)
#define DISABLE_WARNING_UNREFERENCED_FUNCTION            DISABLE_WARNING(4505)
// other warnings you want to deactivate...

//<! special macro to disable armadillo warnings
#define DISABLE_ARMA_WARNINGS								   \
		DISABLE_OVERFLOW;									   \
		DISABLE_WARNING(26812); /* unscoped enum			*/ \
		DISABLE_WARNING(26819); /* unannotated fallthrough	*/ \
		DISABLE_WARNING(26439); /* may not throw			*/ \
		DISABLE_WARNING(6011);  /* dereferencing NULL ptr 	*/ 

// declare some warnings as errors
#pragma warning (error: 4834) // discard returning value of function with 'nodiscard' attribute --> now as error
//#pragma warning (error: 4458) // declaration of '...' hides class member (temporary variables with same name as class members)

#elif defined(__clang__) || defined(__GNUC__)
#define __PRAGMA(X) _Pragma(#X)
#define DISABLE_WARNING_PUSH           __PRAGMA(GCC diagnostic push)
#define DISABLE_WARNING_POP            __PRAGMA(GCC diagnostic pop) 
#define DISABLE_WARNING(warningName)   __PRAGMA(GCC diagnostic ignored #warningName)

#define DISABLE_OVERFLOW								 DISABLE_WARNING(-Wstrict-overflow)
#define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER    DISABLE_WARNING(-Wunused-parameter)
#define DISABLE_WARNING_UNREFERENCED_FUNCTION            DISABLE_WARNING(-Wunused-function)
// other warnings you want to deactivate... 

//<! special macro to disable armadillo warnings
#define DISABLE_ARMA_WARNINGS													\
	DISABLE_WARNING(-Wenum-compare);			/* unscoped enum*/				\
	DISABLE_WARNING(-Wimplicit-fallthrough);	/* unannotated fallthrough*/	\
	DISABLE_WARNING(-Wuninitialized);			/* unitialized */				\
	DISABLE_WARNING(-Wopenmp-simd);				/* ignore OpenMP warning*/

// declare some warnings as errors
#if defined( __clang__ )
	#pragma clang diagnostic error "-Wunused-result"
#else
	#pragma GCC diagnostic error "-Wunused-result"
#endif

#else
	// another compiler: intel,...
#define DISABLE_WARNING_PUSH
#define DISABLE_WARNING_POP
#define DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER
#define DISABLE_WARNING_UNREFERENCED_FUNCTION
// other warnings you want to deactivate... 

#endif

//<! macro to capture code snippet with disabled overflow
#define NO_OVERFLOW(captured_code)	\
	DISABLE_WARNING_PUSH;			\
	DISABLE_OVERFLOW;				\
	captured_code;					\
	DISABLE_WARNING_POP;






		//DISABLE_WARNING(26495); /* unitialized variable				*/ 
		//DISABLE_WARNING(6993);  /* ignore OpenMP: use single-thread	*/ 
		//DISABLE_WARNING(4849);  /* ignor OpenMP:collapse			*/ 