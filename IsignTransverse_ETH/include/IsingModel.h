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
	randomGen ran;										// consistent quick random number generator

	Mat<cpx> H;											// the Hamiltonian
	Mat<cpx> eigenvectors;								// matrix of the eigenvectors in increasing order
	vec eigenvalues;									// eigenvalues vector

	u64 N;												// the Hilbert space size
	std::vector<int> nearest_neighbors;					// vector of nearest neighbors dependend on BC
	std::mutex my_mute_button;							// thread mutex

	std::vector<u64> mapping;							// mapping for the reduced Hilbert space
	std::vector<cpx> normalisation;						// used for normalization in the symmetry case

	virtual u64 map(u64 index) = 0;						// function returning either the mapping(symmetries) or the input index (no-symmetry: 1to1 correspondance)

public:
	u64 E_av_idx;										// average energy
	/* MODEL BASED PARAMETERS */
	int L;												// chain length
	double J;											// spin exchange
	double g;											// transverse magnetic field
	double h;											// perpendicular magnetic field
	int _BC;											// boundary condition

	// CONSTRUCTOR
	virtual ~IsingModel() = 0;

	// Clone functions
	//virtual std::unique_ptr<IsingModel> clone() const = 0;
	//return unique_ptr<IsingModel>(new IsingModel(*this));
	//virtual std::unique_ptr<IsingModel> move_clone() = 0;

	// GETTERS & SETTERS
	std::string get_info() const;																// get the information about the model params
	u64 get_hilbert_size() const;																// get the Hilbert space size 2^N
	const cx_mat& get_hamiltonian() const;														// get the const reference to a Hamiltonian
	const vec& get_eigenvalues() const;															// get the const reference to eigenvalues
	const cx_mat& get_eigenvectors() const;														// get the const reference to the eigenvectors
	const std::vector<u64>& get_mapping() const;												// constant reference to the mapping

	double get_eigenEnergy(u64 idx) const;														// get eigenenergy at a given idx
	const cx_vec& get_eigenState(u64 idx) const;													// get an eigenstate at a given idx

	// PRINTERS
	void print_base_spin_sector(int Sz = 0);													// print basis state with a given total spin (for clarity purposes)
	void print_state(u64 _id);																	// prints the eigenstate at a given idx

	// METHODS
	void set_neighbors();																		// create neighbors list according to the boundary conditions
	virtual void hamiltonian() = 0;																// pure virtual Hamiltonian creator
	virtual void setHamiltonianElem(u64 k, double value, std::vector<bool>& temp) = 0;

	void diagonalization();																		// diagonalize the Hamiltonian

	// VIRTUALS
	virtual mat correlation_matrix(u64 state_id) = 0;											// create the spin correlation matrix at a given state
	static double total_spin(const mat& corr_mat);												// the diagonal part of a spin correlation matrix

	// PHYSICAL QUANTITIES
	double ipr(int state_idx);																	// calculate the ipr coeffincient (inverse participation ratio)
	double information_entropy(u64 _id);														// calculate the information entropy in a given state (based on the ipr) Von Neuman type
	double eigenlevel_statistics(u64 _min, u64 _max);											// calculate the statistics based on eigenlevels (r coefficient)
	vec eigenlevel_statistics_with_return();													// calculate the eigenlevel statistics and return the vector with the results
	virtual double entaglement_entropy(u64 state_id, int subsystem_size) = 0;					// entanglement entropy based on the density matrices

	// PHYSICAL OPERATORS (model states dependent)
	virtual double av_sigma_z(int site, u64 alfa, u64 beta) = 0;								// check the sigma_z matrix element on a given site
	virtual double av_sigma_z(int site_a, int site_b, u64 alfa, u64 beta) = 0;					// check the matrix element for sigma_z correlations S^z_aS^z_b
	double av_sigma_z_extensive(u64 n, u64 m);
	double av_sigma_z_extensive_corr(const u64 n, const u64 m, int corr_len = 1);

	virtual double av_sigma_x(int site, u64 alfa, u64 beta) = 0;								// check the sigma_x matrix element on a given site
	virtual double av_sigma_x_extensive(const u64 n, const u64 m) = 0;

	// USING PHYSICAL QUANTITES FOR PARAMTER RANGES, ETC.
	static void operator_av_in_eigenstates(double (IsingModel::* op)(int, int), IsingModel& A, int site, \
		std::string name = "operator_averaged.txt", string separator = "\t\t");
	static vec operator_av_in_eigenstates_return(double (IsingModel::* op)(int, int), IsingModel& A, int site);
	static double spectrum_repulsion(double (IsingModel::* op)(int, int), IsingModel& A, int site);

	// TOOLS AND HELPERS
	friend std::unordered_map<u64, u64> mapping_sym_to_original(u64 _min, u64 _max,\
		const IsingModel& symmetry, const IsingModel& original);								// maps the eigenstates of a symmetry sector to a Hamiltonian not involving the symmetries
	friend double overlap(const IsingModel& A, const IsingModel& B, int n_a, int n_b);			// creates the overlap between two eigenstates
};


// ----------------------------------------- SYMMETRIC -----------------------------------------
/// <summary>
/// Model with included symmetries and uniform perpendicular magnetic field
/// </summary>
class IsingModel_sym : public IsingModel {
public:
	/* Constructors */
	IsingModel_sym() = default;
	IsingModel_sym(int L, double J, double g, double h, int k_sym = 0, bool p_sym = 1, bool x_sym = 1, int _BC = 0);

private:
	// REDUCED BASIS AS A SYMMETRY SECTOR
	struct {
		double k_sym;				// translational symmetry generator
		int p_sym;					// parity symmetry generator
		int x_sym;					// spin-flip symmetry generator
	} symmetries;

	std::tuple<u64, int> find_translation_representative(std::vector<bool>& base_vector) const;								
	std::tuple<u64, int> find_SEC_representative(const std::vector<bool>& base_vector) const;

	cpx get_symmetry_normalization(std::vector<bool>& base_vector, u64 k);
	void mapping_kernel(u64 start, u64 stop, std::vector<u64>& map_threaded, std::vector<cpx>& norm_threaded, int _id);							// multithreaded mapping
	void generate_mapping();																													// utilizes the mapping kernel

	u64 map(u64 index) override;																												// finds a map corresponding to index (for inheritance purpose)
public:
	// OVERRIDES OF THE MODEL METHODS
	void hamiltonian() override;
	void setHamiltonianElem(u64 k, double value, std::vector<bool>& temp) override;																

	friend std::pair<u64, cpx> find_rep_and_sym_eigval(v_1d<bool>& base,\
		const IsingModel_sym& sector_alfa, cpx normalisation_beta);																				// returns the index and the value of the minimum representative

	double entaglement_entropy(u64 state_id, int subsystem_size) override {
		return 0;
	};
	// OPERATORS

	double av_sigma_z(int site, u64 alfa, u64 beta) override;
	double av_sigma_z(int site_a, int site_b, u64 alfa, u64 beta) override;

	double av_sigma_x(int site, u64 alfa, u64 beta) override;

	// lambda functions for Sigmas - changes the state and returns the value on the base vector
	static std::pair<cpx, v_1d<bool>> sigma_x(int site, const v_1d<bool>& base_vec) {
		auto tmp = base_vec;
		tmp[site] = !tmp[site];
		return std::make_pair(1.0,tmp);
	};
	static std::pair<cpx, v_1d<bool>> sigma_y(int site, const v_1d<bool>& base_vec) {
		auto tmp = base_vec;
		tmp[site] = !tmp[site];
		return std::make_pair(base_vec[site] ? im : -im, tmp);
	}
	static std::pair<cpx, v_1d<bool>> sigma_z(int site, const v_1d<bool>& base_vec) {
		auto tmp = base_vec;
		return std::make_pair(tmp[site] ? 1.0 : -1.0,tmp);
	}


	friend cpx av_operator(int site, u64 alfa, u64 beta, const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta, std::function<std::pair<cpx,v_1d<bool>>(int, v_1d<bool>&)> op);
	friend cpx apply_sym_overlap(int site, const arma::subview_col<cpx>& alfa, const arma::subview_col<cpx>& beta, const v_1d<bool>& base_vec,\
		const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta, std::function<std::pair<cpx,v_1d<bool>>(int, v_1d<bool>&)> op);

	friend double av_sigma_x_sym_sectors(int site, const u64 beta, const u64 alfa, const IsingModel_sym& sector_alfa, const IsingModel_sym& sector_beta);
	double av_sigma_x_extensive(const u64 n, const u64 m) override;

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

	vec dh;																		// disorder in the system - deviation from a constant h value
	double w;																	// the distorder strength to set dh in (-disorder_strength, disorder_strength)
	vec dJ;																		// disorder in the system - deviation from a constant J0 value
	double J0;																	// spin exchange coefficient
	vec dg;																		// disorder in the system - deviation from a constant g0 value
	double g0;																	// transverse magnetic field
public:
	/* Constructors */
	IsingModel_disorder() = default;
	IsingModel_disorder(int L, double J, double J0, double g, double g0, double h, double w, int _BC = 0);

private:
	void generate_mapping();
	u64 map(u64 index) override;

public:
	// METHODS
	void hamiltonian() override;
	void setHamiltonianElem(u64 k, double value, std::vector<bool>& temp) override;

	// MATRICES

	double av_sigma_z(int site, u64 alfa, u64 beta) override;
	double av_sigma_z(int site_a, int site_b, u64 alfa, u64 beta) override;

	double av_sigma_x(int site, u64 alfa, u64 beta) override;
	double av_sigma_x_extensive(const u64 n, const u64 m) override;

	mat correlation_matrix(u64 state_id) override;

	double entaglement_entropy(u64 state_id, int subsystem_size) override;
};


void probability_distribution(std::string dir, std::string name, const arma::vec& data,\
	double _min, double _max, double step = 0.05);												// creates the probability distribution on a given data and saves it to a directory
arma::vec probability_distribution_with_return(const arma::vec& data,\
	double _min, double _max, double step = 0.05);												// creates the probability distribution on a given data and returns a vector with it
arma::vec data_fluctuations(const arma::vec& data, int mu = 10);								// removes the average from the given data based on small buckets of size mu
arma::vec statistics_average(const arma::vec& data, int num_of_outliers = 3);					// takes the average from the vector and the outliers
double quantum_fidelity(u64 _min, u64 _max, const IsingModel& Hamil,\
	double J, double J0, double g, double g0, double h, double w = 0);							// calculates the quantum fidelity for disorder

#endif