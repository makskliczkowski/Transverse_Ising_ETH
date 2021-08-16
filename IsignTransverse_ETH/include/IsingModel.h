#pragma once
#ifndef ISINGMODEL
#define ISINGMODEL
#include "headers.h"

/*-------------------- ISING MODEL WITH TRANSVERSE MAGNETIC FIELD ---------------------*
 * The Ising model with perpendicular magnetic field, known as the quantum Ising model *
 * is a perfect example of a quantum model than can exhibit various quantum phenomena. *
 * There exist plenty of theoretical methods that have tried to test them all but here *
 * we shall focus mainly on the exact diagonalization techniques for the testing of so *
 * called ETH hypothesis. In order to do this, two schemes of work are introduced. The *
 * first one is focusing on the introduction of the symmetries in the model, while the *
 * consideration of the symmetry sector can be controlled. Secondly, because the model *
 * with translational invariance is creating the degeneracy in the eigenstates, an ETH *
 * is impossible to be tested without symmetry sectors. To deal with it, perturbation  *
 * in the perpendicular magnetic field is introduced. For more information, please ref *
 * to:																		 *
 * -> https://github.com/makskliczkowski/Transverse_Ising_ETH						 *
 * The context is partially based on:											 *
 * -> arXiv:2009.09208v1														 *
 * -> 10.1103/PhysRevE.90.052105												 *
 * ------------------------------ All rights reserved -------------------------------- *
 * Authors:																	 *
 * - Rafa³ Œwiêtek, soon to be Phd student										 *
 *	- email: 77swietek77.at.gmail.com											 *	  
 * - Maksymilian Kliczkowski Phd student, Wroc³aw University of Science and Technology *		   
 *	- email: maxgrom97.at.gmail.com												 *
 * ----------------------------------------------------------------------------------- *
 * Special thanks to dr Lev Vidmar at Institute Josef Stefan, with whose support       *
 * the work has been done, while staying in Ljubljana, Slovenia.					 *
 * ---------------------------------------------------------------------------------- */


class IsingModel {
protected:
	std::string info;									// information about the model
	mat H;												// the Hamiltonian
	vec eigenvalues;									// eigenvalues vector
	mat eigenvectors;									// matrix of the eigenvectors in increasing order
	u64 N;												// the Hilbert space size
	std::vector<int> nearest_neighbors;					// vector of nearest neighbors dependend on BC
	std::mutex my_mute_button;							// thread mutex

	std::vector<u64> mapping;							// mapping for the reduced Hilbert space
	std::vector<int> periodicity;						// used for normalization in the symmetry case

	virtual u64 map(u64 index) = 0;						// function returning either the mapping(symmetries) or the input index (no-symmetry: 1to1 correspondance)

	/* SYMMETRY MATRICES */
	sp_mat X;											// sparse construction of Sigma_x Pauli matrix
	sp_mat P;											// sparse construction of parity matrix
	sp_mat T;											// sparse construction of translation matrix

public:
	enum class operators { H, X, P, T };				// implemented operators to be given to a function
	enum class symmetries { T, P, X };					// implemented symmetries

	/* MODEL BASED PARAMETERS */
	int L;												// chain length
	vector<double> J;									// spin exchange
	double g;											// transverse magnetic field
	double h;											// perpendicular magnetic field

	// CONSTRUCTOR 
	virtual ~IsingModel() = 0;

	// Clone functions 
	//virtual std::unique_ptr<IsingModel> clone() const = 0;
	//return unique_ptr<IsingModel>(new IsingModel(*this));
	//virtual std::unique_ptr<IsingModel> move_clone() = 0;

	// GETTERS & SETTERS
	std::string get_info() const;																// get the information about the model params
	u64 get_hilbert_size() const;																// get the Hilbert space size 2^N
	const mat& get_hamiltonian() const;															// get the const reference to a Hamiltonian
	const vec& get_eigenvalues() const;															// get the const reference to eigenvalues
	const mat& get_eigenvectors() const;														// get the const reference to the eigenvectors
	const std::vector<u64>& get_mapping() const;												// constant reference to the mapping

	double get_eigenEnergy(int idx) const;														// get eigenenergy at a given idx
	const vec& get_eigenState(int idx) const;													// get an eigenstate at a given idx

	// PRINTERS
	void print_base_spin_sector(int Sz = 0);													// print basis state with a given total spin (for clarity purposes)
	void print_state(u64 _id);																	// prints the eigenstate at a given idx

	// METHODS
	void set_neighbors();																		// create neighbors list according to the boundary conditions

	virtual void hamiltonian() = 0;																// pure virtual Hamiltonian creator
	virtual void setHamiltonianElem(u64 k, double value, std::vector<bool>&& temp) = 0;

	void diagonalization();																		// diagonalize the Hamiltonian

	// VIRTUALS
	virtual void create_X_matrix() = 0;															// create spin-flip symmetry matrix via Pauli x-matrices
	virtual mat correlation_matrix(u64 state_id) = 0;											// create the spin correlation matrix at a given state
	static double total_spin(const mat& corr_mat);												// the diagonal part of a spin correlation matrix

	// COMMUTATION
	sp_mat choose_operator(IsingModel::operators A);
	bool commutator(IsingModel::operators A, IsingModel::operators B);

	// PHYSICAL QUANTITIES 
	double ipr(int state_idx);																	// calculate the ipr coeffincient
	double eigenlevel_statistics(u64 _min, u64 _max);											// calculate the statistics based on eigenlevels (r coefficient)



	// PHYSICAL OPERATORS (model states dependent) 
	virtual double av_sigma_x(int state_id, int site) = 0;
	virtual double entaglement_entropy(u64 state_id, int subsystem_size) = 0;

	// USING PHYSICAL QUANTITES FOR PARAMTER RANGES, ETC.
	static void operator_av_in_eigenstates(double (IsingModel::* op)(int, int), IsingModel& A, int site, \
		std::string name = "operator_averaged.txt", string separator = "\t\t");
	static vec operator_av_in_eigenstates_return(double (IsingModel::* op)(int, int), IsingModel& A, int site);
	static double spectrum_repulsion(double (IsingModel::* op)(int, int), IsingModel& A, int site);
};

/// <summary>
/// Overlapping of two eigenstates of possibly different matrices A and B
/// </summary>
/// <param name="A">matrix A</param>
/// <param name="B">matrix B</param>
/// <param name="n_a">number of A eigenstate</param>
/// <param name="n_b">number of B eigenstate</param>
/// <returns>A_n_a dot n_b_B</returns>
inline double overlap(const std::unique_ptr<IsingModel>& A, const std::unique_ptr<IsingModel>& B, int n_a, int n_b) {
	return abs(arma::cdot(A->get_eigenState(n_a), B->get_eigenState(n_b)));
}


//-------------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Model with included symmetries and uniform perpendicular magnetic field
/// </summary>
class IsingModel_sym : public IsingModel {
public:
	/* Constructors */
	IsingModel_sym() = default;
	IsingModel_sym(int L, vector<double>& J, double g, double h);
	// IsingModel_sym(const IsingModel_sym& A);
	// IsingModel_sym(IsingModel_sym&& A) noexcept;
	// ~IsingModel_sym();

	/* METHODS */
private:
	void generate_mapping();
	void mapping_kernel(u64 start, u64 stop, std::vector<u64>& map_threaded, int _id);

	void check_periodicity();
	std::vector<u64> find_SEC_representative(const std::vector<bool>& base_vector);
	u64 find_translation_representative(std::vector<bool>& base_vector);

	u64 map(u64 index) override;

public:
	void hamiltonian() override;
	void setHamiltonianElem(u64 k, double value, std::vector<bool>&& temp) override;

	void create_X_matrix() override {};
	double entaglement_entropy(u64 state_id, int subsystem_size) override {
		return 0;
	};

	double av_sigma_x(int state_id, int site) override;
	mat correlation_matrix(u64 state_id) override {
		return mat();
	};
};


//-------------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Model with disorder thus with no symmetries
/// </summary>
class IsingModel_disorder : public IsingModel {
private:
	vec dh;																									// disorder in the system - deviation from a constant h value
	double disorder_strength;																				// the distorder strength to set dh in (-disorder_strength, disorder_strength)
public:
	/* Constructors */
	IsingModel_disorder() = default;
	IsingModel_disorder(int L, vector<double>& J, double g, double h, double disorder_strength);
	// IsingModel_disorder(const IsingModel_disorder& A);
	// IsingModel_disorder(IsingModel_disorder&& A) noexcept;
	// ~IsingModel_disorder();

private:
	void generate_mapping();
	u64 map(u64 index) override;

public:
	// METHODS
	void hamiltonian() override;
	void setHamiltonianElem(u64 k, double value, std::vector<bool>&& temp) override;

	// MATRICES
	void create_X_matrix() override;
	double av_sigma_x(int state_id, int site) override;

	mat correlation_matrix(u64 state_id) override;

	double entaglement_entropy(u64 state_id, int subsystem_size) override;

};


double quantum_fidelity(u64 _min, u64 _max, const std::unique_ptr<IsingModel>& Hamil, const vector<double>& J, double g, double h, double w = 0);




#endif