#pragma once
#ifndef ISINGMODEL
#define ISINGMODEL
#include "headers.h"
class IsingModel {
protected:
	mat H;
	vec eigenvalues;
	mat eigenvectors;
	u64 N;
	std::vector<int> nearest_neighbours;
	std::mutex my_mute_button;

	bool use_mapping;
	std::vector<u64> mapping;
	std::vector<int> periodicity;
public:

	/* MODEL BASED PARAMETERS */
	int L; //chain length
	vector<double> J; //spin exchange
	double g; //transverse magnetic field
	double h; // perpendicular magnetic field

	/* CONSTRUCTOR */
	~IsingModel() = default;

	/* Clone functions */
	//virtual std::unique_ptr<IsingModel> clone() const = 0;
	//return unique_ptr<IsingModel>(new IsingModel(*this));
	//virtual std::unique_ptr<IsingModel> move_clone() = 0;

	/* getters & setters */
	u64 get_hilbert_size();
	mat get_hamiltonian();
	vec get_eigenvalues();
	mat get_eigenvectors();

	/* METHODS */
	void print_base_spin_sector(int Sz = 0);

	virtual void hamiltonian() = 0;
	virtual void setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) = 0;

	void set_neighbours();
	void diagonalization();

	/* PHYSICAL QUANTITIES */
	double ipr(int state_idx);
	virtual double av_sigma_x(int state_id, int site) = 0;

	static void operator_av_in_eigenstates(double (IsingModel::* op)(int, int), IsingModel& A, int site, \
		std::string name = "operator_averaged.txt", string separator = "\t\t");
	static double spectrum_repulsion(double (IsingModel::* op)(int, int), IsingModel& A, int site);

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

	/* METHODS */
private:
	void generate_mapping();
	void mapping_kernel(u64 start, u64 stop, std::vector<u64>& map_threaded, int _id);

	void check_periodicity();
	std::vector<u64> find_SEC_representative(std::vector<bool>& base_vector, u64 min);
	u64 find_translation_representative(std::vector<bool>& base_vector);

public:
	void hamiltonian() override;
	void setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) override;

	double av_sigma_x(int state_id, int site);	
};  


/// <summary>
/// Model with disorder with no symmetries
/// </summary>
class IsingModel_disorder : public IsingModel {
protected:
	vec dh; // disorder in the system - deviation from a constant h value
public:
	/* Constructors */
	IsingModel_disorder() = default;
	IsingModel_disorder(int L, vector<double>& J, double g, double h);
	IsingModel_disorder(const IsingModel_disorder& A);
	IsingModel_disorder(IsingModel_disorder&& A) noexcept;

	/* METHODS */
	void hamiltonian() override;
	void setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) override;

	double av_sigma_x(int state_id, int site); 
	
};


#endif