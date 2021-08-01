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
	virtual void hamiltonian() = 0;
	virtual void setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) = 0;

	void set_neighbours();
	void diagonalization();
	double ipr(int state_idx);
};

/// <summary>
/// Model with included symmetries and uniform perpendicular magnetic field
/// </summary>
class IsingModel_sym : public IsingModel {
protected: 
	std::vector<u64> mapping;
public:
	/* Constructors */
	IsingModel_sym() = default;
	IsingModel_sym(int L, vector<double>& J, double g, double h);
	IsingModel_sym(const IsingModel_sym& A);
	IsingModel_sym(IsingModel_sym&& A) noexcept;

	/* METHODS */
	void generate_mapping();
	void mapping_kernel(u64 start, u64 stop, std::vector<u64>& map_threaded, int _id);

	void hamiltonian() override;
	void setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) override;
};

/// <summary>
/// Model with disorder with no symmetries
/// </summary>
class IsingModel_disorder : public IsingModel {
protected:
	vec dh;
public:
	/* Constructors */
	IsingModel_disorder() = default;
	IsingModel_disorder(int L, vector<double>& J, double g, double h);
	IsingModel_disorder(const IsingModel_disorder& A);
	IsingModel_disorder(IsingModel_disorder&& A) noexcept;

	/* METHODS */
	void hamiltonian() override;
	void setHamiltonianElem(u64& k, double value, std::vector<bool>&& temp) override;
};

#endif