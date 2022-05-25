#pragma once

namespace lanczos {
    typedef std::tuple<arma::vec, arma::vec, arma::vec, arma::vec> _returnTy;

    inline _returnTy FTLM(Lanczos& lanczos_obj) {
        const arma::vec E = lanczos_obj.eigenvalues;
        const size_t N = E.size();
        const int M = lanczos_obj.params.lanczos_steps;
        const int R = lanczos_obj.params.random_steps;

        const double bandwidth = E(N - 1) - E(0);
        const double T_min = bandwidth / M;
        const arma::vec temperature = arma::logspace(std::floor(std::log10(T_min)) - 1, 1, 200);

        arma::vec E_av(temperature.size(), arma::fill::zeros);
        arma::vec E_av2(temperature.size(), arma::fill::zeros);
        arma::vec Z(temperature.size(), arma::fill::zeros);

        for (int r = 0; r < R; r++) {
            lanczos_obj.diagonalization();
            arma::vec overlap(M);
#pragma omp parallel for num_threads(num_of_threads)
            for (int m = 0; m < M; m++)
            {
                cpx temp = arma::cdot(lanczos_obj.randVec_inKrylovSpace, lanczos_obj.eigenvectors.col(m));
                overlap(m) = abs(temp * conj(temp));
            }

#pragma omp parallel for shared(temperature) num_threads(num_of_threads)
            for (int k = 0; k < temperature.size(); k++)
            {
                const double beta = 1. / temperature[k];
                for (int m = 0; m < M; m++)
                {
                    const double temp = overlap(m);
                    const double delE = E(m) - E(0);
                    Z(k) += temp * std::exp(-beta * delE);
                    E_av(k) += E(m) * temp * std::exp(-beta * delE);
                    E_av2(k) += E(m) * E(m) * temp * std::exp(-beta * delE);
                }
            }
        }
        Z = Z * (double)N / (double)R;
        E_av = E_av / Z * (double)N / (double)R;
        E_av2 = E_av2 / Z * (double)N / (double)R;

        return std::make_tuple(temperature, Z, E_av, E_av2);
    };

};



