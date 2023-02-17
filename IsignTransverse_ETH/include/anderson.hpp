

//<! anderson model
namespace anderson{

    //<! build hamiltonian
    inline
    auto hamiltonian(
        const lattice_base& lattice,    //<! lattice class for system
        double t,                       //<! exchange coupling ( 'opping)
        double W                        //<! disorder strength
    ){
        const u64 N = lattice.volume;
        //-- set  anderson single-body matrix
        arma::vec disorder = my_disorder.uniform(N, W);
        
        arma::sp_mat anderson_hamiltonian(N, N);
        for(int j = 0; j < N; j++){
            anderson_hamiltonian(j, j) = disorder(j);
            auto neis = lattice.get_neighbours(j);
            for(auto& nei : neis){
                if(nei > 0){
                    anderson_hamiltonian(j, nei) = t / 2.;
                    anderson_hamiltonian(nei, j) = t / 2.;
                }
            }
        }
        return anderson_hamiltonian;
    }
    //<! main diagonalisation scheme for any dimension lattice
    inline
    auto diagonalize(
        const lattice_base& lattice,    //<! lattice class for system
        double t,                       //<! exchange coupling ( 'opping)
        double W                        //<! disorder stregth
    ){
        arma::vec energies;
        arma::mat orbitals;
        
        arma::sp_mat anderson_hamiltonian = hamiltonian(lattice, t, W);
        
        arma::eig_sym(energies, orbitals, (arma::mat)anderson_hamiltonian);
        return std::make_pair(energies, orbitals);
        
    }

    //<! calculating localisation length for each anderson orbital
    inline
    auto get_localisation_length1D(
        int system_size,
        double t,
        double W
    ) -> std::pair<arma::vec, arma::mat>
    {
        auto lattice = std::make_unique<lattice1D>(system_size);
        arma::vec energies;
        arma::mat orbitals;
        std::tie(energies, orbitals) = anderson::diagonalize(*lattice, t, W);
        arma::vec loc_length(system_size, arma::fill::zeros);
        //return std::make_pair(energies, arma::vec(system_size, arma::fill::zeros));
    //#pragma omp parallel
        arma::mat corr_func(system_size / 2, system_size, arma::fill::zeros);
        for(int i = 0; i < system_size; i++){
            auto orbital_i = orbitals.col(i);

            arma::vec corr(system_size / 2., arma::fill::zeros);
            for(int r = 0; r < system_size / 2.; r++){
                double val = 0;
                for(int j = 0; j < system_size; j++)
                    val += (orbital_i(j) * orbital_i((j + r) % system_size) );
                
                corr(r) = log(abs(val));
            }
            corr_func.col(i) = corr;
        }
        return std::make_pair(energies, corr_func);
    }
    
};

