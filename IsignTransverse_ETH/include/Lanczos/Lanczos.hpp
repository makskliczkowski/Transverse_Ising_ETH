#pragma once
#ifndef _LANCZOS
#define _LANCZOS

namespace lanczos {
	template <typename _type> class Lanczos;	//<! forward declaration for FTLM
};
#include "params.hpp"
#include "FTLM.hpp"						//<! Finite-Temperature Lanczos Method

/// <summary>
/// LANCZOS CLASS
/// </summary>
namespace lanczos {

	template <typename _type>
	class Lanczos {

		arma::Mat<_type> krylov_space;		//<! krylov matrix - basis transformation
		arma::Mat<_type> H_lanczos;			//<! Lanczos matrix for tridiagonalization
		arma::Mat<_type> eigenvectors;		//<! eigenvectors from diagonalizing lanczos matrix
		arma::SpMat<_type> H;				//<! hamiltonian matrix -- change to operator instance
		arma::vec eigenvalues;				//<! lanczos eignevalues
		
		arma::Col<_type> randVec_inKrylovSpace;		//<! random vector written in lanczos basis (used in FTLM)

		randomGen ran;						//<! random variable generator -- uniform distribution
		lanczosParams params;				//<! parameters for lanczos procedure: steps, random steps, efficiency etc
		u64 N;								//<! dimension of hilbert space
		bool use_krylov;					//<! boolean value whether useing krylov matrix or not
		//! ----------------------------------------------------- PRIVATE BUILDERS / INITIALISERS
		void initialize();
		void build_lanczos(const arma::Col<_type>& random_vec);
		void build_krylov(const arma::Col<_type>& random_vec);

		void orthogonalize(
			arma::Col<_type>& vec_to_ortho,  //<! vector to orthogonalize
			int j							 //<! current dimension of Krylov space
		);

	public:
		auto get_eigenvalues() const { return this->eigenvalues; }

		friend _returnTy FTLM(Lanczos<_type>&);
		//------------------------------------------------------------------------------------------------ CONSTRUCTOS
		~Lanczos() = default;
		Lanczos() = delete;
		explicit Lanczos(
			const arma::SpMat<_type>& hamiltonian, 
			lanczosParams&& input_params = lanczosParams()
		) : H(hamiltonian), params(std::move(input_params))
		{ initialize(); }

		//------------------------------------------------------------------------------------------------ DIAGONALIZING MATRIX
		void build(const arma::Col<_type>& random_vec);

		void diagonalization(const arma::Col<_type>& random_vec);
		void diagonalization();
		//TODO: some methods with return values

		//------------------------------------------------------------------------------------------------ CAST STATES TO ORIGINAL HILBERT SPACE:
		//enum class base_type {
		//	hilbert,	//<! Hilbert basis, i.e. computational basis
		//	krylov		//<! Krylov basis build from random vector
		//};
		[[nodiscard]] auto conv_to_hilbert_space(const arma::Col<_type>& rand, int state_id)					-> arma::Col<_type>;
		[[nodiscard]] auto conv_to_hilbert_space(const arma::Col<_type>& rand, const arma::Col<_type>& input)	-> arma::Col<_type>;
		[[nodiscard]] auto conv_to_krylov_space(const arma::Col<_type>& rand, const arma::Col<_type>& input)	-> arma::Col<_type>;
		
		//------------------------------------------------------------------------------------------------ TOOLS
		void convergence(const arma::Col<_type>& rand, int num_state_out = 10);
		
		//------------------------------------------------------------------------------------------------ DYNAMICS
		auto time_evolution_stationary(
			const arma::cx_vec& input_state,
			double time
		) -> arma::cx_vec;

		auto time_evolution_non_stationary(
			const arma::cx_vec& prev_state,
			double dt,
			int lanczos_steps = 10
		) -> arma::cx_vec;
	};
};

#include "lanczos_construct_impl.hpp"	//<! constructors etc implemetation
#include "lanczos_tools.hpp"			//<! lanczos tools implementation: i.e. convergence etc.
#include "lanczos_build_impl.hpp"		//<! lanczos implementation with or without krylov space
#include "lanczos_converter.hpp"		//<! casting vectoes inbetween krylov and hilbert spaces
#include "lanczos_eigs.hpp"				//<! diagonalization of lanczos matrix
#include "lanczos_dynamics.hpp"			//<! implementation of dynamic quantities
#endif