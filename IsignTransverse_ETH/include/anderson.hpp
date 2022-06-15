

//<! anderson model
namespace anderson{

    inline
    auto get_orbitals(
        int system_size, 
        double J, 
        double h
    ) -> std::pair<arma::vec, arma::mat>
    {
        arma::vec energies;
        arma::mat orbitals;
        
        //-- set  anderson single-body matrix
        arma::vec disorder = create_random_vec(system_size, h);
        
        arma::mat anderson_hamiltonian(system_size, system_size, arma::fill::zeros);
        for(int j = 0; j < system_size; j++){
            anderson_hamiltonian(j, j) = disorder(j);
            anderson_hamiltonian(j, (j + 1) % system_size) = J;
            anderson_hamiltonian((j + 1) % system_size, j) = J;
        }
        
        arma::eig_sym(energies, orbitals, anderson_hamiltonian);
        return std::make_pair(energies, orbitals);
    }

    inline
    auto get_localisation_length(
        int system_size,
        double J,
        double h
    ) -> std::pair<arma::vec, arma::vec>
    {
        arma::vec energies;
        arma::mat orbitals;
        std::tie(energies, orbitals) = anderson::get_orbitals(system_size, J, h);
        arma::vec loc_length(system_size, arma::fill::zeros);

    //#pragma omp parallel
        for(int i = 0; i < system_size; i++){
            auto orbital_i = orbitals.col(i);

            arma::vec func_to_fit(system_size / 2, arma::fill::zeros);
            for(int k = 0; k < system_size / 2; k++)
                for(int j = 0; j < system_size; j++)
                    func_to_fit(k) += abs(orbital_i(j) * orbital_i( (j + k) % system_size ) );
            func_to_fit = arma::log(func_to_fit);
            arma::vec p = arma::polyfit(arma::linspace(0, system_size / 2 - 1, system_size / 2), func_to_fit, 1);
            loc_length(i) = -1. / p(0);
        }
        return std::make_pair(energies, loc_length);
    }
    
};

