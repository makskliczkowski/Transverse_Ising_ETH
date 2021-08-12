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
 * to:																				   *
 * -> https://github.com/makskliczkowski/Transverse_Ising_ETH						   *
 * The context is partially based on:												   *
 * -> arXiv:2009.09208v1															   *
 * -> 10.1103/PhysRevE.90.052105													   *
 * ------------------------------ All rights reserved -------------------------------- *
 * Authors:																			   *
 * - Rafa³ Œwiêtek, soon to be Phd student											   *
 *	- email: 77swietek77.at.gmail.com												   *	  
 * - Maksymilian Kliczkowski Phd student, Wroc³aw University of Science and Technology *		   
 *	- email: maxgrom97.at.gmail.com													   *
 * ----------------------------------------------------------------------------------- *
 * Special thanks to dr Lev Vidmar at Institute Josef Stefan, with whose support       *
 * the work has been done, while staying in Ljubljana, Slovenia.					  *
 * ---------------------------------------------------------------------------------- */


class IsingModel {
protected:
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
	sp_mat X;
	sp_mat P;
	sp_mat T;

	
public:
	enum class operators {
		H,
		X,
		P,
		T
	};
	/* MODEL BASED PARAMETERS */
	int L;												// chain length
	vector<double> J;									// spin exchange
	double g;											// transverse magnetic field
	double h;											// perpendicular magnetic field

	/* CONSTRUCTOR */
	virtual ~IsingModel() = 0;

	/* Clone functions */
	//virtual std::unique_ptr<IsingModel> clone() const = 0;
	//return unique_ptr<IsingModel>(new IsingModel(*this));
	//virtual std::unique_ptr<IsingModel> move_clone() = 0;

	/* GETTERS & SETTERS */
	u64 get_hilbert_size() const;								
	mat get_hamiltonian() const;
	vec get_eigenvalues() const;
	mat get_eigenvectors() const;
	std::vector<u64> get_mapping() const;

	double get_eigenEnergy(int idx) const;

	/* METHODS */
	void print_base_spin_sector(int Sz = 0);			// print basis state with a given total spin

	virtual void hamiltonian() = 0;						// pure virtual Hamiltonian creator
	virtual void setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) = 0;

	void set_neighbors();								// create nearest neighbors list
	void diagonalization();								// diagonalize the Hamiltonian

	virtual void create_X_matrix() = 0;					// create spin-flip symmetry matrix via Pauli x-matrices
	bool commutator(IsingModel::operators A, IsingModel::operators B);

	/* PHYSICAL QUANTITIES */
	double ipr(int state_idx);							// calculate the ipr coeffincient
	double eigenlevel_statistics(u64 _min, u64 _max);
	virtual mat correlation_matrix(u64 state_id) = 0;
	static double total_spin(const mat& corr_mat);
	void print_state(u64 _id);

	/* PHYSICAL OPERATORS (model states dependent) */
	virtual double av_sigma_x(int state_id, int site) = 0;

	/* USING PHYSICAL QUANTITES FOR PARAMTER RANGES, ETC.*/
	static void operator_av_in_eigenstates(double (IsingModel::* op)(int, int), IsingModel& A, int site, \
		std::string name = "operator_averaged.txt", string separator = "\t\t");
	static vec operator_av_in_eigenstates_return(double (IsingModel::* op)(int, int), IsingModel& A, int site);

	static double spectrum_repulsion(double (IsingModel::* op)(int, int), IsingModel& A, int site);

	friend double overlap(const IsingModel& A, const IsingModel& B, int n_a, int n_b) {		
		return arma::dot(A.eigenvectors.col(n_a), B.eigenvectors.col(n_b));
	}

};

/// <summary>
/// Model with included symmetries and uniform perpendicular magnetic field
/// </summary>
class IsingModel_sym : public IsingModel {
public:
	/* Constructors */
	IsingModel_sym() = default;
	IsingModel_sym(int L, vector<double>& J, double g, double h);
	IsingModel_sym(const IsingModel_sym& A);
	IsingModel_sym(IsingModel_sym&& A) noexcept;
	~IsingModel_sym();

	/* METHODS */
private:
	void generate_mapping();
	void mapping_kernel(u64 start, u64 stop, std::vector<u64>& map_threaded, int _id);

	void check_periodicity();
	std::vector<u64> find_SEC_representative(std::vector<bool>& base_vector);
	u64 find_translation_representative(std::vector<bool>& base_vector);

	u64 map(u64 index) override;

public:
	void hamiltonian() override;
	void setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) override;

	void create_X_matrix() override {};

	double av_sigma_x(int state_id, int site) override;	
	mat correlation_matrix(u64 state_id) override {
		return mat();
	};
};  


/// <summary>
/// Model with disorder with no symmetries
/// </summary>
class IsingModel_disorder : public IsingModel {
protected:
	vec dh; // disorder in the system - deviation from a constant h value
	double disorder_strength;
public:
	/* Constructors */
	IsingModel_disorder() = default;
	IsingModel_disorder(int L, vector<double>& J, double g, double h, double disorder_strength);
	IsingModel_disorder(const IsingModel_disorder& A);
	IsingModel_disorder(IsingModel_disorder&& A) noexcept;
	~IsingModel_disorder();

private:
	void generate_mapping();
	u64 map(u64 index) override;

public:
	/* METHODS */
	void hamiltonian() override;
	void setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) override;

	void create_X_matrix() override;

	double av_sigma_x(int state_id, int site) override;
	mat correlation_matrix(u64 state_id) override;

	double entaglement_entropy(u64 state_id, int subsystem_size);

};
std::vector<double> quantum_fidelity(u64 _min, u64 _max, int L, std::vector<double> J, double g, double h, double dw = 0.05);




#endif