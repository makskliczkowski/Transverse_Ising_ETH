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
 * - Rafa³ Œwiêtek, soon to be Phd student, Josef Stefan Institute					 *
 *	- email: 77swietek77.at.gmail.com											 *
 * - Maksymilian Kliczkowski Phd student, Wroc³aw University of Science and Technology *
 *	- email: maxgrom97.at.gmail.com												 *
 * ----------------------------------------------------------------------------------- *
 * Special thanks to dr Lev Vidmar at Institute Josef Stefan, with whose support       *
 * the work has been done, while staying in Ljubljana, Slovenia.					 *
 * ---------------------------------------------------------------------------------- */

template <typename T>
class IsingModel {
protected:
	std::string info;									// information about the model
	randomGen ran;										// consistent quick random number generator

	Mat<T> H;											// the Hamiltonian
	Mat<T> eigenvectors;								// matrix of the eigenvectors in increasing order
	vec eigenvalues;									// eigenvalues vector

	u64 N;												// the Hilbert space size
	std::vector<int> nearest_neighbors;					// vector of nearest neighbors dependend on BC
	std::vector<int> next_nearest_neighbors;			// vector of next nearest neighbors dependent on BC
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
	// get the information about the model params
	std::string get_info(std::initializer_list<std::string> skip = {}) const {
		auto tmp = split_str(this->info, ",");	
		std::string tmp_str = "";
		for(int i=0; i < tmp.size(); i++){
			bool save = true;
			for(auto & skip_param: skip)
			{
				// skip the element if we don't want it to be included in the info
				if(split_str(tmp[i],"=")[0] == skip_param)
					save = false;
			}
			if(save) tmp_str += tmp[i]+",";
		}
		tmp_str.pop_back();
		return tmp_str; 
	};

	u64 get_hilbert_size() const { return this->N; };											// get the Hilbert space size 2^N
	const Mat<T>& get_hamiltonian() const {	return this->H;	};									// get the const reference to a Hamiltonian
	const vec& get_eigenvalues() const { return this->eigenvalues; };							// get the const reference to eigenvalues
	const Mat<T>& get_eigenvectors() const { return this->eigenvectors; };						// get the const reference to the eigenvectors
	const std::vector<u64>& get_mapping() const { return this->mapping; };						// constant reference to the mapping

	double get_eigenEnergy(u64 idx) const { return this->eigenvalues(idx); };					// get eigenenergy at a given idx
	Col<T> get_eigenState(u64 idx) const { return this->eigenvectors.col(idx);};			// get an eigenstate at a given idx

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
	double information_entropy(u64 _id, const IsingModel<T>& beta, u64 _min, u64 _max);			// calculate the information entropy in basis of other model from input
	double eigenlevel_statistics(u64 _min, u64 _max);											// calculate the statistics based on eigenlevels (r coefficient)
	vec eigenlevel_statistics_with_return();													// calculate the eigenlevel statistics and return the vector with the results
	virtual double entaglement_entropy(u64 state_id, int subsystem_size) = 0;					// entanglement entropy based on the density matrices

	// PHYSICAL OPERATORS (model states dependent)
	virtual double av_sigma_z(u64 alfa, u64 beta) = 0;											// check the sigma_z matrix element extensive
	virtual double av_sigma_z(u64 alfa, u64 beta, int corr_len) = 0;							// check the sigma_z matrix element with correlation length
	virtual double av_sigma_z(u64 alfa, u64 beta, std::initializer_list<int> sites) = 0;		// check the matrix element of sigma_z elements sites correlation

	virtual double av_sigma_x(u64 alfa, u64 beta) = 0;											// check the sigma_z matrix element extensive
	virtual double av_sigma_x(u64 alfa, u64 beta, int corr_len) = 0;							// check the sigma_z matrix element with correlation length
	virtual double av_sigma_x(u64 alfa, u64 beta, std::initializer_list<int> sites) = 0;		// check the matrix element of sigma_x elements sites correlation

	virtual double av_spin_flip(u64 alfa, u64 beta) = 0;										// check the spin flip element extensive
	virtual double av_spin_flip(u64 alfa, u64 beta, std::initializer_list<int> sites) = 0;		// check the spin flip element at input sites (up to 2)
	
	virtual cpx av_spin_current(u64 alfa, u64 beta) = 0;										// check the spin current extensive
	virtual cpx av_spin_current(u64 alfa, u64 beta, std::initializer_list<int> sites) = 0;		// check the spin current at given sites

	// USING PHYSICAL QUANTITES FOR PARAMTER RANGES, ETC.
	static void operator_av_in_eigenstates(double (IsingModel::* op)(int, int), IsingModel& A, int site, \
		std::string name = "operator_averaged.txt", string separator = "\t\t");
	static vec operator_av_in_eigenstates_return(double (IsingModel::* op)(int, int), IsingModel& A, int site);
	static double spectrum_repulsion(double (IsingModel::* op)(int, int), IsingModel& A, int site);

	// TOOLS AND HELPERS
	
};
template <typename T>
T overlap(const IsingModel<T>& A, const IsingModel<T>& B, int n_a, int n_b);			// creates the overlap between two eigenstates


// ----------------------------------------- SYMMETRIC -----------------------------------------
/// <summary>
/// Model with included symmetries and uniform perpendicular magnetic field
/// </summary>
class IsingModel_sym : public IsingModel<cpx> {
public:
	/* Constructors */
	IsingModel_sym() = default;
	IsingModel_sym(int L, double J, double g, double h, int k_sym = 0, bool p_sym = true, bool x_sym = true, int _BC = 0);

private:
	// REDUCED BASIS AS A SYMMETRY SECTOR
	struct {
		double k_sym;				// translational symmetry generator
		int p_sym;					// parity symmetry generator
		int x_sym;					// spin-flip symmetry generator
	} symmetries;
	bool k_sector;					// if the k-sector allows p symmetry
	v_1d<cpx> k_exponents;			// precalculate the symmetry exponents for current k vector

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
	
	// MATRICES & OPERATORS
	double av_sigma_z(u64 alfa, u64 beta) override;											// check the sigma_z matrix element extensive
	double av_sigma_z(u64 alfa, u64 beta, int corr_len) override;							// check the sigma_z matrix element with correlation length extensive
	double av_sigma_z(u64 alfa, u64 beta, std::initializer_list<int> sites) override;		// check the matrix element of sigma_z elements sites correlation

	double av_sigma_x(u64 alfa, u64 beta) override;											// check the sigma_z matrix element extensive
	double av_sigma_x(u64 alfa, u64 beta, int corr_len) override;							// check the sigma_z matrix element with correlation length extensive
	double av_sigma_x(u64 alfa, u64 beta, std::initializer_list<int> sites) override;		// check the matrix element of sigma_x elements sites correlation

	double av_spin_flip(u64 alfa, u64 beta) override;										// check the spin flip element extensive
	double av_spin_flip(u64 alfa, u64 beta, std::initializer_list<int> sites) override;		// check the spin flip element at input sites (up to 2)
	
	cpx av_spin_current(u64 alfa, u64 beta) override;										// check the extensive spin current
	cpx av_spin_current(u64 alfa, u64 beta, std::initializer_list<int> sites) override;		// check the spin current at given sites

	// lambda functions for Sigmas - changes the state and returns the value on the base vector
	static std::pair<cpx, v_1d<bool>> sigma_x(const v_1d<bool>& base_vec, std::initializer_list<int> sites) {
		auto tmp = base_vec;
		for(auto& site: sites)
			tmp[site] = !tmp[site];
		return std::make_pair(1.0,tmp);
	};
	static std::pair<cpx, v_1d<bool>> sigma_y(const v_1d<bool>& base_vec, std::initializer_list<int> sites) {
		auto tmp = base_vec;
		cpx val = 1.0;
		for(auto& site: sites){
			val *= tmp[site] ? im : -im;
			tmp[site] = !tmp[site];
		}
		return std::make_pair(val, tmp);
	};
	static std::pair<cpx, v_1d<bool>> sigma_z(const v_1d<bool>& base_vec, std::initializer_list<int> sites) {
		auto tmp = base_vec;
		double val = 1.0;
		for (auto& site : sites)
			val *= tmp[site] ? 1.0 : -1.0;
		return std::make_pair(val, tmp);
	};
	static std::pair<cpx, v_1d<bool>> spin_flip(const v_1d<bool>& base_vec, std::initializer_list<int> sites) {
		if (sites.size() > 2) throw "Not implemented such exotic operators, choose 1 or 2 sites\n";
		auto tmp = base_vec;
		cpx val = 0.0;
		auto it = sites.begin() + 1;
		auto it2 = sites.begin();
		if (!base_vec[*it]) {
			tmp[*it] = !tmp[*it];
			val = 2.0;
			if (sites.size() > 1) {
				if (base_vec[*it2]) {
					tmp[*it2] = !tmp[*it2];
					val *= 2.0;
				}
				else val = 0.0;
			}
		}
		else val = 0.0;
		return std::make_pair(val, tmp);
	};
	
	friend cpx av_operator(u64 alfa, u64 beta, const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta,\
		std::function<std::pair<cpx,v_1d<bool>>(const v_1d<bool>&, std::initializer_list<int>)> op, std::initializer_list<int> sites);							// calculates the matrix element of operator at given site
	friend cpx av_operator(u64 alfa, u64 beta, const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta,\
		std::function<std::pair<cpx,v_1d<bool>>(const v_1d<bool>&, std::initializer_list<int>)> op);																// calculates the matrix element of operator at given site in extensive form (a sum)
	friend cpx av_operator(u64 alfa, u64 beta, const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta,\
		std::function<std::pair<cpx,v_1d<bool>>(const v_1d<bool>&, std::initializer_list<int>)> op, int corr_len);												// calculates the matrix element of operator at given site in extensive form (a sum) with corr_len

	friend cpx apply_sym_overlap(const arma::subview_col<cpx>& alfa, const arma::subview_col<cpx>& beta, const v_1d<bool>& base_vec,\
		const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta,\
		std::function<std::pair<cpx,v_1d<bool>>(const v_1d<bool>&, std::initializer_list<int>)> op,\
		std::initializer_list<int> sites);

	mat correlation_matrix(u64 state_id) override {
		return mat();
	};
};

//-------------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Model with disorder thus with no symmetries
/// </summary>
class IsingModel_disorder : public IsingModel<double> {
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

	// MATRICES & OPERATORS
	double av_sigma_z(u64 alfa, u64 beta) override;											// check the sigma_z matrix element extensive
	double av_sigma_z(u64 alfa, u64 beta, int corr_len) override;							// check the sigma_z matrix element with correlation length
	double av_sigma_z(u64 alfa, u64 beta, std::initializer_list<int> sites) override;		// check the matrix element of sigma_z elements sites correlation

	double av_sigma_x(u64 alfa, u64 beta) override;											// check the sigma_z matrix element extensive
	double av_sigma_x(u64 alfa, u64 beta, int corr_len) override;							// check the sigma_z matrix element with correlation length
	double av_sigma_x(u64 alfa, u64 beta, std::initializer_list<int> sites) override;		// check the matrix element of sigma_x elements sites correlation

	double av_spin_flip(u64 alfa, u64 beta) override;										// check the spin flip element extensive
	double av_spin_flip(u64 alfa, u64 beta, std::initializer_list<int> sites) override;		// check the spin flip element at input sites (up to 2)
	
	cpx av_spin_current(u64 alfa, u64 beta) override;										// check the extensive spin current
	cpx av_spin_current(u64 alfa, u64 beta, std::initializer_list<int> sites) override;		// check the spin current at given sites

	mat correlation_matrix(u64 state_id) override;

	double entaglement_entropy(u64 state_id, int subsystem_size) override;
};


void probability_distribution(std::string dir, std::string name, const arma::vec& data,\
	double _min, double _max, double step = 0.05);												// creates the probability distribution on a given data and saves it to a directory
arma::vec probability_distribution_with_return(const arma::vec& data,\
	double _min, double _max, double step = 0.05);												// creates the probability distribution on a given data and returns a vector with it
arma::vec data_fluctuations(const arma::vec& data, int mu = 10);								// removes the average from the given data based on small buckets of size mu
arma::vec statistics_average(const arma::vec& data, int num_of_outliers = 3);					// takes the average from the vector and the outliers

/// <summary>
/// nie pisz, póŸniej
/// </summary>
/// <typeparam name="T1"></typeparam>
/// <typeparam name="T2"></typeparam>
/// <param name="_min"></param>
/// <param name="_max"></param>
/// <param name="symmetry"></param>
/// <param name="original"></param>
/// <returns></returns>
template <typename T1, typename T2>
std::unordered_map<u64, u64> mapping_sym_to_original(u64 _min, u64 _max, \
	const IsingModel<T1>& symmetry, const IsingModel<T2>& original) {
	std::unordered_map<u64, u64> map;
	std::vector<double> E_dis = arma::conv_to<std::vector<double>>::from(original.get_eigenvalues());
#pragma omp parallel for
	for (int k = 0; k < symmetry.get_hilbert_size(); k++) {
		double E = symmetry.get_eigenEnergy(k);
		if (E < original.get_eigenEnergy(_min) && E >= original.get_eigenEnergy(_max)) continue;
		auto idx = binary_search(E_dis, _min, _max, E);
		if (idx >= original.get_hilbert_size()) continue;
		double E_prev = (idx == 0) ? (original.get_eigenEnergy(0) - 1.0) : original.get_eigenEnergy(idx - 1);
		double E_next = (idx == original.get_hilbert_size() - 1) ? \
			(original.get_eigenEnergy(original.get_hilbert_size() - 1) + 1.0) : original.get_eigenEnergy(idx + 1);
		if (abs(E - E_prev) > 1e-8 && abs(E - E_next) > 1e-8) {
#pragma omp critical
			map[k] = idx;
		}
	}
	return map;
};



#endif