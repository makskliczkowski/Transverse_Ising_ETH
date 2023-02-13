#pragma once

#ifndef __COMMONS
    #include "../commons.hpp"
    #include "../hilbert_space/hilbert_space_base.hpp"
    #include "../random_and_disorder/disorder.hpp"
#endif

/// @brief Geeneral 
/// @tparam _ty 
/// @tparam hilbert 
/// @tparam ...params_types 
template <  typename _ty, 
            class hilbert
            >
class hamiltonian_base{

protected:
    hilbert _hilbert_space;

    arma::SpMat<_ty> H;

    int _boundary_condition;
    u64 dim;
    int system_size;
    //other important params?
    virtual void init() = 0;
public:
    virtual ~hamiltonian_base() = 0;
    auto& get_hamiltonian()	    const { return this->H; }			 // get the const reference to a Hamiltonian
    virtual void create_hamiltonian() = 0;
};


class XXZ : 
    public hamiltonian_base<double, su2_space<su2::spin>>
{
    double _J1;
    double _J2;
    double _delta1;
    double _delta2;
    double _hz;
    double _w;

    disorder::uniform<double> _disorder;
    virtual void init() override
    {
        _disorder = disorder::uniform<double>(system_size, this->_w);
    }
public:
    XXZ() = default;
    XXZ(int L, double J1, double delta1, double w, double hz = 0, double J2 = 0, double delta2 = 0)
        : _J1(J1), _delta1(delta1), _w(w), _hz(hz), _J2(J2), _delta2(delta2)
        { this->system_size = L; init(); }

    virtual void create_hamiltonian() override
    {
        v_2d<double> parameters = { {this->_J1, this->_J1, this->_delta1},
                                    {this->_J2, this->_J2, this->_delta2}
                                  };
        for(auto& x : parameters)
            std::cout << x << std::endl;

        std::vector<op_type> _operators = {operators::sigma_x, operators::sigma_y, operators::sigma_z };
        std::vector<int> neighbor_distance = {1, 2};

        for (size_t k = 0; k < this->dim; k++) {
            int base_state = map(k);
            for (int j = 0; j < this->system_size; j++) {
                
                //<! longitudinal field with disorder
                auto [val, op_k] = operators::sigma_z(base_state, this->system_size, { j });
                double fieldZ = this->_hz + this->_disorder(j);
                this->setHamiltonianElem(k, fieldZ * real(val), op_k);

                //<! go over all neighbours
                for(int a = 0; a < neighbor_distance.size(); a++){
                    int r = neighbor_distance[a];
                    int nei = j + r;
                    if(nei >= this->system_size)
                        nei = (this->_boundary_condition)? -1 : nei % this->system_size;

                    if (nei >= 0) {
                        for(int b = 0; b < _operators.size(); b++){
                            op_type op = _operators[b];
                            auto [val1, op_k] = op(base_state, this->system_size, { j });
                            auto [val2, opop_k] = op(op_k, this->system_size, { nei });
                            this->setHamiltonianElem(k, parameters[a][b] * real(val1 * val2), opop_k);
                        }
                    }
                }
            }
        }
    }
};
// rethink how to include real sectrs

// macro for clean vs disordered model choosing at compiling?
