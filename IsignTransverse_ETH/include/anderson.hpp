

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
        arma::vec disorder = my_gen.create_random_vec<double>(N, W);
        
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
            continue;


            //<! OLD METHOD TO CALCULATE FOR EACH REALISATION SEPERATRELY
            //<! ADD CHOICE BETWEEN OPTIONS!
            arma::vec r_vals  = arma::linspace(0, system_size / 2., corr.size());
            if(i % 100 == 0 || system_size < 20)
                save_to_file("./results/ANDERSON/1D/PBC/CorrelationFunction/_L=" 
                            + std::to_string(system_size) + "_n=" + std::to_string(i) + "_w=" + to_string_prec(W, 2) + ".dat", r_vals, corr);
            //continue;

            double _min = arma::min(corr);
            if(!std::isfinite(_min) || _min < -20.0)
                _min = -20.0;
            arma::vec func_to_fit;
            func_to_fit = exctract_vector_between_values(corr, _min, 1.0);
            r_vals  = arma::linspace(0, func_to_fit.size(), func_to_fit.size());
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
        return std::make_pair(energies, corr_func);
    }
    
};

