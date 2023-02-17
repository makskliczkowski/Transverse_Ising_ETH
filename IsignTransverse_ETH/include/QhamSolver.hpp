#pragma once
//#include "Hamiltonian/XXZ.hpp"
#include "Hamiltonian/QuantumSun.hpp"

template <class Hamiltonian>
class QHamSolver{
	static_check((std::is_base_of_v<_hamiltonian, Hamiltonian>), 
                    "\n" BAD_INHERITANCE "\n\t base class is: hamiltonian_base<element_type, hilbert_space>");
    
    protected:
        typedef typename Hamiltonian::element_type _ty;

        Hamiltonian H;
        arma::Mat<_ty> eigenvectors;
        arma::vec eigenvalues;
        u64 dim;
        u64 E_av_idx;
    public:
        //<! ----------------------------------------------- CONSTRUCTORS / DESTRUCTORS
        ~QHamSolver() { DESTRUCTOR_CALL; }

        template <typename... _param_types>
        QHamSolver(_param_types... args);

        //<! ----------------------------------------------- GETTERS
	    auto get_hilbert_size()	 const { return this->dim; }					// get the Hilbert space size
	    auto get_mapping()		 const { return this->H.get_mapping(); }		// constant reference to the mapping
	    auto get_hamiltonian()	 const { return this->H.get_hamiltonian(); }	// get the const reference to a Hamiltonian
	    auto& get_eigenvectors() const { return this->eigenvectors; }			// get the const reference to the eigenvectors
	    auto get_eigenvalues()	 const { return this->eigenvalues; }			// get the const reference to eigenvalues

        //<! ----------------------------------------------- ROUTINES
        void diagonalization(bool get_eigenvectors = true, const char* method = "dc");
};


//<! ---------------------------------------------------------------------------------------------------------------------------------------
//<! ------------------------------------------------------------------------------------------------------------------------ IMPLEMENTATION


//<! ------------------------------------------------------------------------------ CONSTRUCTORS
template <class Hamiltonian>
template <typename... _param_types>
QHamSolver<Hamiltonian>::QHamSolver(_param_types... args)
{
    CONSTRUCTOR_CALL;
    //<! Initlialize model
    H = Hamiltonian(args...);
    this->dim = H.get_hilbert_space_size();

    //<! what else?
    #if defined(EXTRA_DEBUG)
        std::cout << FUN_SIGNATURE << "::\n\t QHamSolver initialized with: "
            << var_name_value(this->dim, 0) << "\n MODEL:\n" << H << std::endl;
    #endif
}


//<! ------------------------------------------------------------------------------ ROUTINES
template <class Hamiltonian>
void QHamSolver<Hamiltonian>::diagonalization(bool get_eigenvectors, const char* method) 
{
    //out << real(H) << endl;
    arma::Mat<_ty> H_temp;
    try {
        H_temp = this->H.get_dense_hamiltonian();
        if (get_eigenvectors)   arma::eig_sym(this->eigenvalues, this->eigenvectors, H_temp, method);
        else                    arma::eig_sym(this->eigenvalues, H_temp);
        #ifndef NO_DEBUG
            std::cout << "\t HAMILTONIAN TYPE: " + type_name<_ty>() + "\n\tsparse - dim(H) = " + matrix_size(this->H.n_nonzero * sizeof(this->H(0, 0)))
                + "\n\tdense - dim(H) = " + matrix_size(H_temp.n_cols * H_temp.n_rows * sizeof(H_temp(0, 0)))
                + "\n\tspectrum size: " + std::to_string(this->dim) << std::endl << std::endl;
        #endif
    }
    catch (...) {
        handle_exception(std::current_exception(), 
            "sparse - dim(H) = " + matrix_size(this->H.n_nonzero * sizeof(this->H(0, 0)))
            + "\ndense - dim(H) = " + matrix_size(H_temp.n_cols * H_temp.n_rows * sizeof(H_temp(0, 0)))
            + "\n spectrum size: " + std::to_string(this->dim)
        );
    }
    double E_av = arma::trace(eigenvalues) / double(this->dim);
    auto i = min_element(begin(eigenvalues), end(eigenvalues), [=](double x, double y) {
        return abs(x - E_av) < abs(y - E_av);
        });
    this->E_av_idx = i - begin(eigenvalues);

    #ifdef EXTRA_DEBUG
        printSeparated(std::cout, "\t", 16, true, "guessed index", "mean energy", "energies close to this value (-1,0,+1) around found index");
        printSeparated(std::cout, "\t", 16, true, this->E_av_idx, E_av, this->eigenvalues(this->E_av_idx - 1), this->eigenvalues(this->E_av_idx),  this->eigenvalues(this->E_av_idx + 1));
    #endif
}