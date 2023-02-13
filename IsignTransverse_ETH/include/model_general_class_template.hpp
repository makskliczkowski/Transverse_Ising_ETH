#pragma once


// add concept "is constructible"
template <class Hamiltonian, template ..._types>
class model_disordered{
    protected:
        Hamiltonian H;
    ...

    public:
        model_disordered(..., _types... args)
        {
            H = Hamiltonian(args...);
        }
        auto get_eigenvalues() const { return H.get_eigenvalues(); }
}