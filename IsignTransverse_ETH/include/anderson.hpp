

//<! anderson model
namespace anderson{

    //<! main diagonalisation scheme for 1D chain
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
            anderson_hamiltonian(j, (j + 1) % system_size) = J / 2.;
            anderson_hamiltonian((j + 1) % system_size, j) = J / 2.;
            
        }
        
        arma::eig_sym(energies, orbitals, anderson_hamiltonian);
        return std::make_pair(energies, orbitals);
    }

    //<! calculating localisation length for each anderson orbital
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
        const int num_fit = 200;
    //#pragma omp parallel
        for(int i = 0; i < system_size; i++){
            auto orbital_i = orbitals.col(i);

            arma::vec func_to_fit(2 * num_fit + 1, arma::fill::zeros);
            for(int k = -num_fit; k < num_fit; k++){
                double val = 0;
                for(int j = 0; j < system_size; j++){
                    long idx = (j + k) % system_size;
                    if(idx < 0) idx += system_size; // cause modulo in c++ work in negative space
                    val += abs(orbital_i(j) * orbital_i(idx) );
                }
                func_to_fit(k + num_fit) = log(val);
            }
            for(int k = -num_fit; k < num_fit; k++){
                if(k < 0) func_to_fit(k + num_fit) = func_to_fit(0) - func_to_fit(k + num_fit);
            }

            arma::vec r_vals  = arma::linspace(-num_fit, num_fit, 2 * num_fit + 1);
            arma::vec p = arma::polyfit(r_vals, func_to_fit, 1);
            
            loc_length(i) = -1. / p(0);
            
            if(i % 100 > 0) continue;
            func_to_fit.resize(system_size);
            for(int k = -system_size / 2.; k < system_size / 2.; k++){
                double val = 0;
                for(int j = 0; j < system_size; j++){
                    long idx = (j + k) % system_size;
                    if(idx < 0) idx += system_size; // cause modulo in c++ work in negative space
                    val += abs(orbital_i(j) * orbital_i(idx) );
                }
                func_to_fit(k + system_size / 2.) = log(val);
            }
            r_vals  = arma::linspace(-system_size / 2., system_size / 2., system_size);
            save_to_file("./results/disorder/PBC/ObitalCorr_n=" + std::to_string(i) + "_w=" + to_string_prec(h, 2) + ".dat", r_vals, func_to_fit, energies(i), loc_length(i));
        }
        return std::make_pair(energies, loc_length);
    }
    
};

