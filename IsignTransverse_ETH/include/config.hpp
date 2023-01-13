#pragma once
#if defined(MY_MAC)
/*
0 - Ising
1 - Heisenberg
2 - Anderson
3 - XYZ
4 - Quantum Sun
*/
	#define MODEL 3

	#if !defined(LOCAL_PERT)
		//#define LOCAL_PERT
	#endif
#endif


#ifdef LOCAL_PERT
	#pragma message ("Using Heisenberg Model with local perturbation at center")
#endif
#if !defined(DEGENERACIES)
	//#define DEGENERACIES
#endif


#if MODEL == 0
	#define ISING
#elif MODEL == 1
	#if !defined(HEISENBERG)
		#define HEISENBERG
		#pragma message ("Using Heisenberg Model")
	#endif
#elif MODEL == 2
	#if !defined(ANDERSON)
		#define ANDERSON
		#pragma message ("Using Anderson Model (single particle)")
	#endif
#elif MODEL == 3
	#if !defined(XYZ)
		#define XYZ
		#pragma message ("Using XYZ next nearest neighbour model (includes nnn XXZ)")
	#endif
#elif MODEL == 4
	#if !defined(QUANTUM_SUN)
		#define QUANTUM_SUN
		#pragma message ("Using Quantum Sun Model")
	#endif
#else
	#define ISING
#endif