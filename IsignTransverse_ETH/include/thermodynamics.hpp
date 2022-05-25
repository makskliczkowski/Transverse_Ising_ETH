#pragma once

namespace thermodynamics{
    using triple_vec = std::tuple<arma::vec, arma::vec, arma::vec>;
    using triple_dbl = std::vector<double>;

    //<! calculate mean energy, entropy and heat capacity for given input temperature range
    inline
    auto
    thermodynamics(
        double T,                       //<! temperature point
        const arma::vec& eigenvalues,   //<! eigenvalues calculated from ED
        int system_size = 0             //<! system size, if set to 0 do division is performed (assumed user does it)
        ) -> triple_dbl 
    {
        int L = system_size == 0? 1.0 : system_size;
        const size_t N = eigenvalues.size();
    	double Z = 0;
    	double E_av = 0, E_av2 = 0;
    	for (long int i = 0; i < N; i++) {
    		double gibbs = std::exp(-(eigenvalues(i) - eigenvalues(0)) / T);
    		Z += gibbs;
    		E_av += eigenvalues(i) * gibbs;
    		E_av2 += eigenvalues(i) * eigenvalues(i) * gibbs;
    	}
    	E_av /= Z;
    	E_av2 /= Z;
    	double S = (std::log(Z) + (E_av - eigenvalues(0)) / T) / (double)L;
    	double Cv = (E_av2 - E_av * E_av) / double(L * T * T);
    	return std::vector( {E_av, S, Cv} );
    };
    //<! calculate mean energy, entropy and heat capacity for given input temperature range
    inline
    auto
    thermodynamics(
        const arma::vec& temperature,   //<! temperature range
        const arma::vec& eigenvalues,   //<! eigenvalues calculated from ED
        int system_size = 0             //<! system size, if set to 0 do division is performed (assumed user does it)
        ) -> triple_vec 
    {
        const size_t N = eigenvalues.size();
    	arma::vec Cv(temperature.size(), arma::fill::zeros);
    	arma::vec S(temperature.size(), arma::fill::zeros);
    	arma::vec E(temperature.size(), arma::fill::zeros);
    #pragma omp parallel for shared(temperature)
    	for (int k = 0; k < temperature.size(); k++) {
    		const double T = temperature(k);
    		auto thermals = thermodynamics(T, eigenvalues, system_size);
    		E(k) = thermals[0];
    		S(k) = thermals[1];
    		Cv(k) = thermals[2];
    	}
    	return std::make_tuple(E, S, Cv);
    }
    
    //<! assign temperature to given energy (negative temperatures mean state above mean energy at T->inf)
    inline
    auto
    assign_temperature(
        const arma::vec& eigenvalues,
        double energy
    ) -> double
    {
        const int size = 5000;
        double T = 0;
        auto temperature = arma::logspace(-5, 2, size);
        int counter = 0;
        double E_prev = eigenvalues(0);
        for(int k = 1; k < size; k++){
            double E = thermodynamics(temperature[k], eigenvalues)[0];
            if(energy >= E_prev && energy < E){ 
                T = temperature[k];
                break;
            }
            E_prev = E;
            counter++;
        }
        if(counter >= size - 1) {
            T = 1.0e6;
            E_prev = eigenvalues(0);
            for(int k = 1; k < size; k++){
                double E = thermodynamics(-temperature[k], eigenvalues)[0];
                if(energy <= E_prev && energy > E){ 
                    T = -temperature[k];
                    break;
                }
                E_prev = E;
            }
        }
        return T;
    }



}
