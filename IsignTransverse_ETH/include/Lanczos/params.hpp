#pragma once
#ifndef _LANCZOS_PARAMS
#define _LANCZOS_PARAMS

/// <summary>
/// settings for lanczos procedure
/// </summary>
struct lanczosParams {
	int lanczos_steps = 200;					// number of lanczos iterations
	int random_steps  = 1;						// number of random vectors in FTLM
	bool memory_over_performance = false;		// building hamiltonian as sparse (false) or diagonalizing on-the-fly (true)
	bool use_reorthogonalization = true;		// parameter to define whether use full reorthogonalization

	lanczosParams(
		int M, int R, 
		bool use_reortho = false, 
		bool mem_over_perf = false
	)
		:
		lanczos_steps(M),
		random_steps(R),
		use_reorthogonalization(use_reortho),
		memory_over_performance(mem_over_perf)
	{};
	lanczosParams()											 = default;
	lanczosParams(const lanczosParams& input)				 = default;
	lanczosParams(lanczosParams&& input)			noexcept = default;
	lanczosParams& operator=(lanczosParams&& input) noexcept = default;
};

#endif
