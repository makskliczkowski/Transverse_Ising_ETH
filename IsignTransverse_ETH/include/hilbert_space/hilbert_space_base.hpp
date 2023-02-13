#pragma once

#ifndef __COMMONS
    #include "../commons.hpp"
#endif

#define NOT_ALLOWED_SYM_SECTOR "ERROR 1: Chosen symmetry sector not allowed, exceeds min/max values;"


/// @brief Base class for hilbert space construction
/// @tparam ...constraints 
template <typename... constraints>
class hilbert_space_base {
    
    protected:
        variadic_struct<constraints...> sectors;

        std::unordered_map<u64, u64> mapping;
        int system_size;
        u64 dim;
        virtual void init() = 0;
    public:
        virtual ~hilbert_space_base() = 0;

        auto get_mapping() const { return this->mapping; }
        virtual u64 get_real_index(u64 idx) { return (*this)(idx); };
        virtual void create_basis() = 0;

        virtual u64 operator()(u64 idx) { return this->mapping[idx]; }

        virtual u64 find(u64 idx) = 0;
};

class full_space : public hilbert_space_base<>{
    virtual void init() override {};
    public:
        full_space() = default;
        full_space(int L, int sector)
        { 
            this->system_size = L; 
            this->dim = ULLPOW(L);
            this->init();
        }
        virtual u64 operator()(u64 idx) override { return idx; };
        virtual u64 find(u64 idx)       override { return idx; };
};


//<! Enum for possible SU(2) symmetries: for now charge and spin
enum class su2 {spin, charge};

/// @brief Hilbert space creator with SU(2) symmetry, either spin or charge
/// @tparam boolean value: spinless fermions?  (valid if chosen su2 == charge)
/// @tparam su2_sym choose SU(2) symmetry: spin, charge, ...
template <su2 su2_sym = su2::spin, bool spinless = true>
class su2_space : public hilbert_space_base<int>{
    int su2_sector;
    float max_sector;
    float min_sector;

    virtual void init() override{
        switch(su2_sym){
        case su2::spin :
            max_sector = -this->system_size * S / 2;
            min_sector = this->system_size * S / 2;
        case su2::charge :
            max_sector = (spinless? 1 : 2) * this->system_size;
            min_sector = 0;
        }
        assert((su2_sector < max_sector && su2_sector > min_sector) && NOT_ALLOWED_SYM_SECTOR);
        this->sectors = {su2_sector};
        this->create_basis();
    }
    
public:
    su2_space() = default;
    su2_space(int L, int sector)
    { 
        this->system_size = L; 
        this->su2_sector = sector; 
        this->init();
    }

    //<! -------------------------------------------------------- OVERLOADED OPERATORS
    
};