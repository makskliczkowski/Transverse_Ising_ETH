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
* - Rafa� �wi�tek, soon to be Phd student, Josef Stefan Institute					 *
*	- email: 77swietek77.at.gmail.com											 *
* - Maksymilian Kliczkowski Phd student, Wroc�aw University of Science and Technology *
*	- email: maxgrom97.at.gmail.com												 *
* ----------------------------------------------------------------------------------- *
* Special thanks to dr Lev Vidmar at Institute Josef Stefan, with whose support       *
* the work has been done, while staying in Ljubljana, Slovenia.					 *
* ---------------------------------------------------------------------------------- */

template <typename _type>
class IsingModel {
protected:
	std::string info;									// information about the model
	randomGen ran;										// consistent quick random number generator

	arma::SpMat<_type> H;								// the Hamiltonian
	arma::Mat<_type> eigenvectors;						// matrix of the eigenvectors in increasing order
	arma::vec eigenvalues;								// eigenvalues vector

	arma::cx_vec coeff;									// coefficients of state in eigenbasis

	u64 N;												// the Hilbert space size
	std::vector<int> nearest_neighbors;					// vector of nearest neighbors dependend on BC
	std::vector<int> next_nearest_neighbors;			// vector of next nearest neighbors dependent on BC
	std::mutex my_mute_button;							// thread mutex

	std::vector<u64> mapping;							// mapping for the reduced Hilbert space
	std::vector<cpx> normalisation;						// used for normalization in the symmetry case

	virtual u64 map(u64 index) const = 0;					// function returning either the mapping(symmetries) or the input index (no-symmetry: 1to1 correspondance)

public:
	u64 E_av_idx = -1;										// average energy
	/* MODEL BASED PARAMETERS */
	int L = 8;												// chain length
	double J = 1;											// spin exchange
	double g = 1;											// transverse magnetic field
	double h = 1;											// perpendicular magnetic field
	int _BC = 1;											// boundary condition

	// ---------------------------------- CONSTRUCTOR ----------------------------------
	virtual ~IsingModel() = 0;

	// ---------------------------------- GETTERS ----------------------------------
	//
	// get the information about the model params
	auto get_info(std::vector<std::string> skip = {}) const->std::string {
		auto tmp = split_str(this->info, ",");
		std::string tmp_str = "";
		for (int i = 0; i < tmp.size(); i++) {
			bool save = true;
			for (auto& skip_param : skip)
			{
				// skip the element if we don't want it to be included in the info
				if (split_str(tmp[i], "=")[0] == skip_param)
					save = false;
			}
			if (save) tmp_str += tmp[i] + ",";
		}
		tmp_str.pop_back();
		return tmp_str;
	};

	auto get_hilbert_size()						const { return this->N; }						 // get the Hilbert space size 2^N
	auto get_mapping()							const { return this->mapping; }					 // constant reference to the mapping
	auto& get_hamiltonian()						const { return this->H; }						 // get the const reference to a Hamiltonian
	auto& get_eigenvectors()					const { return this->eigenvectors; }			 // get the const reference to the eigenvectors
	auto get_eigenvalues()						const { return this->eigenvalues; }				 // get the const reference to eigenvalues
	auto get_eigenEnergy(u64 idx)				const { return this->eigenvalues(idx); }		 // get eigenenergy at a given idx
	auto get_eigenState(u64 idx)				const { return this->eigenvectors.col(idx); }	 // get an eigenstate at a given idx
	auto get_eigenStateValue(u64 idx, u64 elem) const { return this->eigenvectors(elem, idx);}	 // get an eigenstate at a given idx
	auto calculate_energy(const arma::cx_vec& input) const{
		return real(arma::cdot(input, this->H * input));
	};
	// ---------------------------------- PRINTERS ----------------------------------
	void print_base_spin_sector(int Sz = 0);													// print basis state with a given total spin (for clarity purposes)
	void print_state(u64 _id);																	// prints the eigenstate at a given idx

	// ---------------------------------- GENERAL METHODS ----------------------------------
	void set_neighbors();																		// create neighbors list according to the boundary conditions
	void diagonalization(bool get_eigenvectors = true, const char* method = "dc");				// diagonalize the Hamiltonian

	virtual void hamiltonian() = 0;																// pure virtual Hamiltonian creator
	virtual void hamiltonian_heisenberg() = 0;													// pure virtual heisenberg hamiltonian creator
	virtual void setHamiltonianElem(u64 k, double value, u64 new_idx) = 0;						// sets the Hamiltonian elements in a virtual way
	void reset_random(size_t new_seed = seed_global) const {
		gen = std::mt19937_64(new_seed);
		//this->ran.reset();
	}
	double getRandomValue(const double _min, const double _max) {
		return this->ran.randomReal_uni(_min, _max);
	}
	int properSite(int site) const {
		if (_BC == 0) {
			if (site < 0)
				return L - (abs(site) % (L));
			if (site >= L)
				return site % L;
			return site % L;
		}
		else {
			if (site < 0 || site >= L)
				assert(false && "exceeding chain with open boundary conditions!\n Please choos PBC or go to hell!");
			return site;
		}
		stout << "RETURNING INT_MAX, which should not be possible" << std::endl;
		return INT_MAX;
	}
	std::vector<int> properSite(std::vector<int> sites) const {
		std::vector<int> proper;
		for(auto& s : sites)
			proper.push_back(this->properSite(s));
		return proper;
	}
	// ---------------------------------- VIRTUALS ----------------------------------
	virtual arma::mat correlation_matrix(u64 state_id) const = 0;											// create the spin correlation matrix at a given state
	static double total_spin(const arma::mat& corr_mat);												// the diagonal part of a spin correlation matrix

	// ---------------------------------- PHYSICAL QUANTITIES ----------------------------------
	double ipr(int state_idx) const;																// calculate the ipr coeffincient (inverse participation ratio)
	double information_entropy(u64 _id) const;														// calculate the information entropy in a given state (based on the ipr) Von Neuman type
	double information_entropy(u64 _id, const IsingModel<_type>& beta, u64 _min, u64 _max) const;	// calculate the information entropy in basis of other model from input
	double eigenlevel_statistics(u64 _min, u64 _max) const;											// calculate the statistics based on eigenlevels (r coefficient)
	arma::vec eigenlevel_statistics_with_return() const;											// calculate the eigenlevel statistics and return the vector with the results
	virtual double mean_level_spacing_analytical() const = 0;										// mean level spacing from analytical formula calcula
	
	// ---------------------------------- THERMODYNAMIC QUANTITIES ----------------------------------
	std::tuple<arma::vec, arma::vec, arma::vec> thermal_quantities(const arma::vec& temperature);	// calculate thermal quantities with input temperature range
																									// heat capacity, entropy, average energy
	// lambda functions for Sigmas - changes the state and returns the value on the base vector
	static std::pair<cpx, u64> sigma_x(u64 base_vec, int L, std::vector<int> sites) {
		for (auto& site : sites) {
			//site = properSite(site);
			base_vec = flip(base_vec, BinaryPowers[L - 1 - site], L - 1 - site);
		}
		return std::make_pair(1.0, base_vec);
	};
	static std::pair<cpx, u64> sigma_y(u64 base_vec, int L, std::vector<int> sites) {
		auto tmp = base_vec;
		cpx val = 1.0;
		for (auto& site : sites) {
			//site = properSite(site);
			val *= checkBit(tmp, L - 1 - site) ? im : -im;
			tmp = flip(tmp, BinaryPowers[L - 1 - site], L - 1 - site);
		}
		return std::make_pair(val, tmp);
	};
	static std::pair<cpx, u64> sigma_z(u64 base_vec, int L, std::vector<int> sites) {
		auto tmp = base_vec;
		double val = 1.0;
		for (auto& site : sites) {
			//site = properSite(site);
			val *= checkBit(tmp, L - 1 - site) ? 1.0 : -1.0;
		}
		return std::make_pair(val, base_vec);
	};
	static std::pair<cpx, u64> spin_flip(u64 base_vec, int L, std::vector<int> sites) {
		if (sites.size() > 2) throw "Not implemented such exotic operators, choose 1 or 2 sites\n";
		auto tmp = base_vec;
		cpx val = 0.0;
		auto it = sites.begin() + 1;
		auto it2 = sites.begin();
		if (!(checkBit(base_vec, L - 1 - *it))) {
			NO_OVERFLOW(tmp = flip(tmp, BinaryPowers[L - 1 - *it], L - 1 - *it);)
			val = 2.0;
			if (sites.size() > 1) {
				if (checkBit(base_vec, L - 1 - *it2)) {
					NO_OVERFLOW(tmp = flip(tmp, BinaryPowers[L - 1 - *it2], L - 1 - *it2);)
					val *= 2.0;
				}
				else val = 0.0;
			}
		}
		else val = 0.0;
		return std::make_pair(val, tmp);
	};


	// ---------------------------------- PHYSICAL OPERATORS (model states dependent) ----------------------------------
	virtual double av_sigma_z(u64 alfa, u64 beta) = 0;											// check the sigma_z matrix element extensive
	virtual double av_sigma_z(u64 alfa, u64 beta, int corr_len) = 0;							// check the sigma_z matrix element with correlation length
	virtual double av_sigma_z(u64 alfa, u64 beta, std::vector<int> sites) = 0;		// check the matrix element of sigma_z elements sites correlation
			
	virtual double av_sigma_x(u64 alfa, u64 beta) = 0;											// check the sigma_z matrix element extensive
	virtual double av_sigma_x(u64 alfa, u64 beta, int corr_len) = 0;							// check the sigma_z matrix element with correlation length
	virtual double av_sigma_x(u64 alfa, u64 beta, std::vector<int> sites) = 0;		// check the matrix element of sigma_x elements sites correlation
			
	virtual double av_spin_flip(u64 alfa, u64 beta) = 0;										// check the spin flip element extensive
	virtual double av_spin_flip(u64 alfa, u64 beta, std::vector<int> sites) = 0;		// check the spin flip element at input sites (up to 2)

	virtual cpx av_spin_current(u64 alfa, u64 beta) = 0;										// check the spin current extensive
	virtual cpx av_spin_current(u64 alfa, u64 beta, std::vector<int> sites) = 0;		// check the spin current at given sites


	// ---------------------------------- USING PHYSICAL QUANTITES FOR PARAMTER RANGES, ETC. ----------------------------------
	static void operator_av_in_eigenstates(double (IsingModel::* op)(int, int), IsingModel& A, int site, \
		std::string name = "operator_averaged.txt", string separator = "\t\t");
	static arma::vec operator_av_in_eigenstates_return(double (IsingModel::* op)(int, int), IsingModel& A, int site);
	static double spectrum_repulsion(double (IsingModel::* op)(int, int), IsingModel& A, int site);

	virtual arma::sp_cx_mat create_operator(std::initializer_list<op_type> operators) const = 0;
	virtual arma::sp_cx_mat create_operator(std::initializer_list<op_type> operators, int corr_len) const = 0;
	virtual arma::sp_cx_mat create_operator(std::initializer_list<op_type> operators, std::vector<int> sites) const = 0;
	// transverse-field Ising LIOMs operator densities
	//arma::sp_cx_mat create_StringOperator(coordinate alfa, coordinate beta, int j, int k) [[expects: k > 0]];
	//arma::sp_cx_mat create_LIOMoperator_densities(int n, int ell) const;
	arma::sp_cx_mat create_tfim_liom_plus(int n) const;
	arma::sp_cx_mat create_tfim_liom_minus(int n) const;
	
	virtual arma::sp_cx_mat createHq(int k) const = 0;
	virtual arma::sp_cx_mat createHlocal(int k) const = 0;
	virtual arma::sp_cx_mat fourierTransform(op_type op, int q) const = 0; // with operator input
	//template <typename op>
	//virtual sp_cx_mat fourierTransform(int q, op&& opGen) const = 0; // with lambda input to create local operator and perform sum 

	arma::sp_cx_mat chooseOperator(int choose, int site) {
		arma::sp_cx_mat op;
		switch (choose) {
			case 0: op = this->create_operator({ IsingModel::sigma_z }, std::vector<int>({ site })); break;
			case 1: op = this->create_operator({ IsingModel::sigma_x }, std::vector<int>({ site })); break;
			case 2: op = this->createHlocal(site); break;
			case 3: op = this->fourierTransform(IsingModel::sigma_z, site); break;
			case 4: op = this->fourierTransform(IsingModel::sigma_x, site); break;
			case 5: op = this->createHq(site); break;
			case 6: op = this->create_tfim_liom_plus(site); break;
			case 7: op = this->create_tfim_liom_minus(site); break;
			default:
				stout << "No operator chosen!\nReturning empty matrix\n\n";
		}
		return op;
	}
	static auto opName(int choose, int site) 
		-> std::pair<std::string, std::string>
	{
		std::string name, dir_suffix;
		switch (choose) {
		case 0: name = "SigmaZ_j=" 	  		+ std::to_string(site);	break;
		case 1: name = "SigmaX_j=" 	  		+ std::to_string(site);	break;
		case 2: name = "H_j="	   	  		+ std::to_string(site);	break;
		case 3: name = "SigmaZ_q=" 	  		+ std::to_string(site);	break;
		case 4: name = "SigmaX_q=" 	  		+ std::to_string(site);	break;
		case 5: name = "H_q="	   	  		+ std::to_string(site);	break;
		case 6: name = "TFIM_LIOM_plus_n="  + std::to_string(site);	break;
		case 7: name = "TFIM_LIOM_minus_n=" + std::to_string(site);	break;
		default:
			stout << "Bad input! Operator -op 0-7 only";
			exit(1);
		}
		if(choose < 3) 						dir_suffix = "j";
		else if(choose >=3 && choose < 6) 	dir_suffix = "q";
		else 								dir_suffix = "n";
		return std::make_pair(name, dir_suffix);
	}


	auto time_evolve_state(const arma::cx_vec& state, double time) -> arma::cx_vec; 	//<! stationary time evolution
	void time_evolve_state_ns(arma::cx_vec& state, double dt, int order = 3); 			// non-stationary time evolution (time-dependent model)
	// entanglement entropy based on the density matrices
	virtual arma::cx_mat reduced_density_matrix(const arma::cx_vec& state, int A_size) const  = 0;

	double entaglement_entropy(	const arma::cx_vec& state, int A_size)						const;
	double reyni_entropy(		const arma::cx_vec& state, int A_size, unsigned alfa = 2)	const;
	double shannon_entropy(		const arma::cx_vec& state, int A_size)						const;

	arma::vec entaglement_entropy(const arma::cx_vec& state) const;

	void set_coefficients(const arma::cx_vec& initial_state);

	//--------------------------------------------------------- dummy functions
	virtual arma::vec get_non_interacting_energies() = 0;
	virtual arma::vec first_interacting_correction() = 0;
};

inline void normaliseOp(arma::sp_cx_mat& op) {
	const u64 N = op.n_cols;
	const cpx operator_normalisation = std::sqrt(arma::trace(op * op) / double(N));
	if (abs(operator_normalisation) >= 1e-15)
		op = op / operator_normalisation; // normalize if non-zero norm
	else
		std::cout << "normalisation <A^2>=0 - somethings wrong, boi.. " << std::scientific << operator_normalisation << std::endl;
}
inline void normaliseMat(arma::cx_mat& _matrix) {
	const u64 N = _matrix.n_cols;
	cpx operator_normalisation = arma::norm(_matrix, "fro"); // has already square root
	operator_normalisation = operator_normalisation / std::sqrt(N);
	_matrix = _matrix / operator_normalisation;
}
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
	struct sym{
		double k_sym;															// translational symmetry generator
		int p_sym;																// parity symmetry generator
		int x_sym;																// spin-flip symmetry generator
		bool operator==(sym other) {
			return (this->k_sym == other.k_sym 
				&& this->p_sym == other.p_sym 
				&& this->x_sym == other.x_sym);
		}
	} symmetries;
	bool k_sector;																// if the k-sector allows p symmetry
	v_1d<cpx> k_exponents;														// precalculate the symmetry exponents for current k vector
	std::vector<std::function<u64(u64, int)>> symmetry_group;					// full symmetry group containing all operations: T^k, P*T^k, ..
	std::vector<cpx> symmetry_eigVal;											// eigenvalues for each symmetry operation
	void createSymmetryGroup();													// create symmetry group elements and their eigenvalues

	//std::tuple<u64, int> find_translation_representative(u64 base_idx) const;
	std::pair<u64, cpx> find_SEC_representative(u64 base_idx) const;

	cpx get_symmetry_normalization(u64 base_idx) const;
	void mapping_kernel(u64 start, u64 stop, std::vector<u64>& map_threaded, std::vector<cpx>& norm_threaded, int _id);							// multithreaded mapping
	void generate_mapping();																													// utilizes the mapping kernel

	u64 map(u64 index) const override;																												// finds a map corresponding to index (for inheritance purpose)
public:
	bool get_k_sector() const { return this->k_sector; }
	v_1d<cpx> get_k_exponents() const { return this->k_exponents; }
	sym get_symmetries() const { return this->symmetries;}
	v_1d<cpx> get_norm() const { return this->normalisation; }
	v_1d<std::function<u64(u64, int)>> get_sym_group() const { return this->symmetry_group; }
	v_1d<cpx> get_sym_eigVal() const { return this->symmetry_eigVal; }
	// OVERRIDES OF THE MODEL METHODS
	void hamiltonian() override;
	void hamiltonian_heisenberg() override;
	void setHamiltonianElem(u64 k, double value, u64 new_idx) override;
	double mean_level_spacing_analytical() const override {
		const double chi = 0.341345;
		return sqrt(L) / (chi * N) * sqrt(J * J + h * h + g * g);
	}

	friend std::pair<u64, cpx> find_rep_and_sym_eigval(u64 base_idx, \
		const IsingModel_sym& sector_alfa, cpx normalisation_beta);																				// returns the index and the value of the minimum representative

	// MATRICES & OPERATORS
	double av_sigma_z(u64 alfa, u64 beta) override;											// check the sigma_z matrix element extensive
	double av_sigma_z(u64 alfa, u64 beta, int corr_len) override;							// check the sigma_z matrix element with correlation length extensive
	double av_sigma_z(u64 alfa, u64 beta, std::vector<int> sites) override;		// check the matrix element of sigma_z elements sites correlation

	double av_sigma_x(u64 alfa, u64 beta) override;											// check the sigma_z matrix element extensive
	double av_sigma_x(u64 alfa, u64 beta, int corr_len) override;							// check the sigma_z matrix element with correlation length extensive
	double av_sigma_x(u64 alfa, u64 beta, std::vector<int> sites) override;		// check the matrix element of sigma_x elements sites correlation

	double av_spin_flip(u64 alfa, u64 beta) override;										// check the spin flip element extensive
	double av_spin_flip(u64 alfa, u64 beta, std::vector<int> sites) override;		// check the spin flip element at input sites (up to 2)

	cpx av_spin_current(u64 alfa, u64 beta) override;										// check the extensive spin current
	cpx av_spin_current(u64 alfa, u64 beta, std::vector<int> sites) override;		// check the spin current at given sites

	
	static std::string set_info(int L, double J, double g, double h, int k_sym, bool p_sym, bool x_sym, std::vector<std::string> skip = {}, std::string sep = "_") {
		std::string name = sep + "L=" + std::to_string(L) + \
			",J=" + to_string_prec(J, 2) + \
			",g=" + to_string_prec(g, 2) + \
			",h=" + to_string_prec(h, 2) + \
			",k=" + std::to_string(k_sym) + \
			",p=" + std::to_string((p_sym) ? 1 : -1) + \
			",x=" + std::to_string(x_sym ? 1 : -1);
		auto tmp = split_str(name, ",");
		std::string tmp_str = "";
		for (int i = 0; i < tmp.size(); i++) {
			bool save = true;
			for (auto& skip_param : skip)
			{
				// skip the element if we don't want it to be included in the info
				if (split_str(tmp[i], "=")[0] == skip_param)
					save = false;
			}
			if (save) tmp_str += tmp[i] + ",";
		}
		tmp_str.pop_back();
		return tmp_str;
	}

	friend cpx av_operator(u64 alfa, u64 beta, const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta, op_type op, std::vector<int> sites);					// calculates the matrix element of operator at given site
	friend cpx av_operator(u64 alfa, u64 beta, const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta, op_type op);											// calculates the matrix element of operator at given site in extensive form (a sum)
	friend cpx av_operator(u64 alfa, u64 beta, const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta, op_type op, int corr_len);							// calculates the matrix element of operator at given site in extensive form (a sum) with corr_len

	friend cpx apply_sym_overlap(const arma::subview_col<cpx>& alfa, const arma::subview_col<cpx>& beta, u64 base_vec, u64 k, \
		const IsingModel_sym& sec_alfa, const IsingModel_sym& sec_beta, op_type op, std::vector<int> sites);

	arma::sp_cx_mat symmetryRotation() const;
	arma::cx_vec symmetryRotation(const arma::cx_vec& state) const;

	//friend sp_cx_mat create_operatorDistinctSectors()
	arma::sp_cx_mat create_operator(std::initializer_list<op_type> operators) const override;													
	arma::sp_cx_mat create_operator(std::initializer_list<op_type> operators, int corr_len) const override;
	arma::sp_cx_mat create_operator(std::initializer_list<op_type> operators, std::vector<int> sites) const override;
	void set_OperatorElem(std::vector<op_type> operators, std::vector<int> sites, arma::sp_cx_mat& operator_matrix, u64 base_vec, u64 cur_idx) const;
	
	arma::sp_cx_mat createHq(int k) const override { stout << "Not implemented yet!!\n\n"; return arma::sp_cx_mat(); };
	arma::sp_cx_mat createHlocal(int k) const override { stout << "Not implemented yet!!\n\n"; return arma::sp_cx_mat(); };
	arma::sp_cx_mat fourierTransform(op_type op, int q) const override;

	arma::cx_mat reduced_density_matrix(const arma::cx_vec& state, int A_size) const override;
	arma::mat correlation_matrix(u64 state_id) const override;

	//--------------------------------------------------------- dummy functions
	virtual arma::vec get_non_interacting_energies() override;
	virtual arma::vec first_interacting_correction() override;
};
//-------------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Model with disorder thus with no symmetries
/// </summary>
class IsingModel_disorder : public IsingModel<double> {
private:

	arma::vec dh;																		// disorder in the system - deviation from a constant h value
	double w;																	// the distorder strength to set dh in (-disorder_strength, disorder_strength)
	arma::vec dJ;																		// disorder in the system - deviation from a constant J0 value
	double J0;																	// spin exchange coefficient
	arma::vec dg;																		// disorder in the system - deviation from a constant g0 value
	double g0;																	// transverse magnetic field
public:
	/* Constructors */
	IsingModel_disorder() = default;
	IsingModel_disorder(int L, double J, double J0, double g, double g0, double h, double w, int _BC = 0);

private:
	void generate_mapping();
	u64 map(u64 index) const override;

public:
	// METHODS
	void hamiltonian() override;
	void hamiltonian_heisenberg() override;
	void setHamiltonianElem(u64 k, double value, u64 new_idx) override;
	double mean_level_spacing_analytical() const override {
		const double chi = 0.341345;
		return sqrt(L) / (chi * N) * sqrt(J * J + h * h + g * g + (w * w + g0 * g0 + J0 * J0) / 3.);
	}

	// MATRICES & OPERATORS
	double av_sigma_z(u64 alfa, u64 beta) override;											// check the sigma_z matrix element extensive
	double av_sigma_z(u64 alfa, u64 beta, int corr_len) override;							// check the sigma_z matrix element with correlation length
	double av_sigma_z(u64 alfa, u64 beta, std::vector<int> sites) override;		// check the matrix element of sigma_z elements sites correlation

	double av_sigma_x(u64 alfa, u64 beta) override;											// check the sigma_z matrix element extensive
	double av_sigma_x(u64 alfa, u64 beta, int corr_len) override;							// check the sigma_z matrix element with correlation length
	double av_sigma_x(u64 alfa, u64 beta, std::vector<int> sites) override;		// check the matrix element of sigma_x elements sites correlation

	double av_spin_flip(u64 alfa, u64 beta) override;										// check the spin flip element extensive
	double av_spin_flip(u64 alfa, u64 beta, std::vector<int> sites) override;		// check the spin flip element at input sites (up to 2)

	cpx av_spin_current(u64 alfa, u64 beta) override;										// check the extensive spin current
	cpx av_spin_current(u64 alfa, u64 beta, std::vector<int> sites) override;		// check the spin current at given sites

	arma::sp_cx_mat create_operator(std::initializer_list<op_type> operators) const override;
	arma::sp_cx_mat create_operator(std::initializer_list<op_type> operators, int corr_len) const override;
	arma::sp_cx_mat create_operator(std::initializer_list<op_type> operators, std::vector<int> sites) const override;

	arma::sp_cx_mat createHq(int k) const override;
	arma::sp_cx_mat createHlocal(int k) const override;
	arma::sp_cx_mat fourierTransform(op_type op, int q) const override;
	arma::mat correlation_matrix(u64 state_id) const override;

	cpx av_operator(u64 alfa, u64 beta, op_type op, std::vector<int> sites);	// calculates the matrix element of operator at given site
	cpx av_operator(u64 alfa, u64 beta, op_type op);							// calculates the matrix element of operator at given site in extensive form (a sum)
	cpx av_operator(u64 alfa, u64 beta, op_type op, int corr_len);

	arma::cx_mat reduced_density_matrix(const arma::cx_vec& state, int A_size) const override;

	static std::string set_info(int L, double J, double J0, double g, double g0, double h, double w, std::vector<std::string> skip = {}, std::string sep = "_") {
		std::string name = sep + "L=" + std::to_string(L) + \
			",J=" + to_string_prec(J, 2) + \
			",J0=" + to_string_prec(J0, 2) + \
			",g=" + to_string_prec(g, 2) + \
			",g0=" + to_string_prec(g0, 2) + \
			",h=" + to_string_prec(h, 2) + \
			",w=" + to_string_prec(w, 2);
		auto tmp = split_str(name, ",");
		std::string tmp_str = "";
		for (int i = 0; i < tmp.size(); i++) {
			bool save = true;
			for (auto& skip_param : skip)
			{
				// skip the element if we don't want it to be included in the info
				if (split_str(tmp[i], "=")[0] == skip_param)
					save = false;
			}
			if (save) tmp_str += tmp[i] + ",";
		}
		tmp_str.pop_back();
		return tmp_str;
	}

	//--------------------------------------------------------- dummy functions
	virtual arma::vec get_non_interacting_energies() override;
	virtual arma::vec first_interacting_correction() override;
};
// ---------------------------------- HELPERS ----------------------------------
template <typename _type>
_type overlap(const IsingModel<_type>& A, const IsingModel<_type>& B, int n_a, int n_b);								// creates the overlap between two eigenstates

// ---------------------------------- PROBABILITY DISTRIBUTORS ----------------------------------

void probability_distribution(std::string dir, std::string name, const arma::vec& data, int n_bins = -1);	// creates the probability distribution on a given data and saves it to a directory
arma::vec probability_distribution_with_return(const arma::vec& data, int n_bins = -1);						// creates the probability distribution on a given data and returns a vector with it
arma::vec data_fluctuations(const arma::vec& data, int mu = 10);											// removes the average from the given data based on small buckets of size mu
arma::vec statistics_average(const arma::vec& data, int num_of_outliers = 3);								// takes the average from the vector and the outliers

/// <summary>
/// Mapping original energies and matches them by indices
/// </summary>
template <typename T1, typename T2>
std::unordered_map<uint64_t, uint64_t> mapping_sym_to_original(uint64_t _min, uint64_t _max, const IsingModel<T1>& symmetry, const IsingModel<T2>& original) {
	std::unordered_map<uint64_t, uint64_t> map;
	std::vector<double> E_dis = arma::conv_to<std::vector<double>>::from(original.get_eigenvalues());
	const u64 N_tot = original.get_hilbert_size();
#pragma omp parallel for
	for (int k = 0; k < symmetry.get_hilbert_size(); k++) {
		double E = symmetry.get_eigenEnergy(k);
		if (E < original.get_eigenEnergy(_min) && E >= original.get_eigenEnergy(_max)) continue;
		auto idx = binary_search(E_dis, _min, _max, E);
		if (idx >= original.get_hilbert_size()) throw "Energy not found";
#if defined (DEGENERACIES)
	#pragma omp critical
		map[k] = idx;
#else 
		double E_prev = (idx == 0)			? (original.get_eigenEnergy(0) - 1.0)			: original.get_eigenEnergy(idx - 1);
		double E_next = (idx == N_tot - 1)	? (original.get_eigenEnergy(N_tot - 1) + 1.0)	: original.get_eigenEnergy(idx + 1);
		if (abs(E - E_prev) > 1e-8 && abs(E - E_next) > 1e-8) {
	#pragma omp critical
			map[k] = idx;
		}
#endif
	}
	return map;
};





using IsingSym 		= std::unique_ptr<IsingModel_sym>;
using IsingDisorder = std::unique_ptr<IsingModel_disorder>;
#endif