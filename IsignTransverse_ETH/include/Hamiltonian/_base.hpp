#pragma once

#ifndef __COMMONS
    #include "../headers.h"
    #include "../hilbert_space/hilbert_space_base.hpp"
#endif
#include "../SpinOperators.hpp"

class _hamiltonian {};

/// @brief Geeneral 
/// @tparam _ty 
/// @tparam hilbert 
/// @tparam ...params_types 
template <  typename _ty, 
            class hilbert
            >
class hamiltonian_base : public _hamiltonian, public arma::SpMat<_ty> {
public:
    using element_type =_ty ;

protected:
    hilbert _hilbert_space;
    typedef arma::SpMat<_ty> sparse_matrix;
    typedef arma::Mat<_ty> matrix;
    sparse_matrix H;

    int _boundary_condition;
    u64 dim;
    int system_size;
    //other important params?
    virtual void init() = 0;

public:
    virtual ~hamiltonian_base() = 0;

    //<! GETTERS TO PRIVATE MEMBERS
    auto get_hilbert_space_size()   const { return this->dim; }                             // get dimension of hilbert space
    auto& get_hamiltonian()	        const { return this->H; }			                    // get the const reference to a Hamiltonian
    auto get_mapping()              const { return this->_hilbert_space.get_mapping(); }    // get mapping to reduced basis
    matrix get_dense_hamiltonian()  const { return matrix(this->H); }                       // return the model hamiltonian as dense matrix

    //<! HAMILTONIAN BUILDERS
    virtual void create_hamiltonian() = 0;
    virtual sparse_matrix create_local_hamiltonian(int site) = 0;
    virtual void set_hamiltonian_elements(u64 k, _ty value, u64 new_idx) = 0;

    //<! OVERRIDING OPERATORS
    virtual std::ostream& write(std::ostream&) const = 0;
    virtual std::istream& read(std::istream&) = 0;
};

template <  typename _ty, class hilbert>
hamiltonian_base<_ty, hilbert>::~hamiltonian_base() {}


//<! ------------------------------------------------------------------------------ OVERRIDEN STREAMs IN/OUT
template <  typename _ty, class hilbert>
inline 
std::ostream& operator<<(std::ostream& os, const hamiltonian_base<_ty, hilbert>& obj) 
    { return obj.write(os); }

template <  typename _ty, class hilbert>
inline 
std::istream& operator>>(std::istream& is, hamiltonian_base<_ty, hilbert>& obj) 
    { return obj.read(is); }

// rethink how to include real sectrs

// macro for clean vs disordered model choosing at compiling?
