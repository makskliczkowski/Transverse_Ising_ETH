#pragma once
#include "random.h"

#define enable_if_convertible(type_in, type_base) \
            static_assert(traits::is_convertible_v<type_in, type_base>, __FILE__"(line=" LINE_STR "): " NOT_CONVERTIBLE)


//!< NAMESPACE WITH DISORDER CLASSES
namespace disorder{

    /// @brief BASE CLASS FOR DISORDER LANDSCAPE
    /// @tparam _ty template argument
    template <typename _ty>
    class disorder_base{

        protected:
            arma::Col<_ty> disorder_values;  //<! disorder landscape
            int length = 0;                  //<! length of disorder array
            _ty bound = 0;                   //<! bound for disorder, i.e. the values are distirbuted in [- _bound, _bound]
            randomGen generator;             //<! random generator to distribute numbers

        //<! BASE CLASS MEMBER FUNCTION: PRIVATE
            template <typename index_type>
            void bound_check(index_type idx){ 
                enable_if_convertible(index_type, int);
                #ifdef NO_DEBUG 
                    return;
                #else
                    assert(idx < this->length && idx >= 0 && "ERROR 1: Index exceeds length of array.");  
                #endif
            } 

            /// @brief Initializing function for default constructor
            virtual void init(int L, _ty W) {
                this->length = L;
                this->bound = W;
                this->disorder_values = arma::Col<_ty>(this->length, arma::fill::zeros);
                if(W > 0)
                    this->set();
                CONSTRUCTOR_CALL;
                #if defined(EXTRA_DEBUG)
                    std::cout << FUN_SIGNATURE << "::\n\t disorder initialized with: "
                        << var_name_value(this->L, 0) << "\t" << var_name_value(this->bound, 0) << std::endl;
                #endif
		    }

        public:
            disorder_base() = default;
            disorder_base(int L, double W) { init(L, W); };
            ~disorder_base() { DESTRUCTOR_CALL; };
            
        //<! GETTERS FOR PROTECTED/PRIVATE BASE CLASS MEMBER ------------------------------------------------------------------------
            _nodiscard auto get_disorder_values() const { return this->disorder_values; };
            _nodiscard auto get_length()          const { return this->length; };
            _nodiscard auto get_bound()           const { return this->bound; };

        //<! BASE CLASS MEMBER FUNCTION: PUBLIC -------------------------------------------------------------------------------------
            virtual void set() = 0;

            /// @brief reset values of disorder given by generator
            /// @param L new length (by default old is assumed)
            virtual
            void reset(int L = -1, _ty W = _ty(0)){
                if(L > 0){
                    this->length = L;
                    this->disorder_values = arma::Col<_ty>(this->length, arma::fill::zeros);
                }
                if(W > 0)
                    this->bound = W;
                this->set();
            }

        //<! CONSTRUCTORS -----------------------------------------------------------------------------------------------------------
            
            template <typename other_ty>
            disorder_base(const disorder_base<other_ty>& other)
                { *this = other; }
            
            template <typename other_ty>
            disorder_base(disorder_base<other_ty>&& other)
                { *this = std::move(other); }

            template <typename other_ty>
            auto operator=(const disorder_base<other_ty>& other)
            {       
                enable_if_convertible(other_ty, _ty);
                
                this->disorder_values = other.get_values();
                this->length = other.get_length();
                this->bound = other.get_bound();                       
            }

        //<! ACCESS OPERATORS --------------------------------------------------------------------------------------------------------
            template <typename index_type>
            auto& operator()(index_type idx) 
            { 
                bound_check(idx); 
                return this->disorder_values(idx); 
            }
            
            template <typename index_type>
            auto& operator[](index_type idx) 
            { 
                enable_if_convertible(int, index_type);   
                return this->disorder_values[idx]; 
            }
            
            /// @brief  Get disorder values at given input indices
            /// @tparam ...index_type Type of input indices:    has to be convertible to integer-types
            /// @param ...idx Indices at which disorder will be put out
            /// @return Returns array of values of length same as number of input indices
            template <typename... index_type>
            auto operator()(index_type... idx) 
            { 
                const std::size_t n = sizeof...(index_type);
                arma::Col<_ty> values_at_indices(n);
                
                int i = 0;
                ([&]{   bound_check(idx);    values_at_indices(i++) = this->disorder_values(idx);  } (), ...);
                
                return values_at_indices;
            }

            template <typename _type>
            friend std::ostream& operator<<(std::ostream& os, const disorder_base<_type>& class_to_write)
                { os << class_to_write.values; return os; }
        
        //<! BOLLEAN OPERATOR --------------------------------------------------------------------------------------------------------
            template <typename other_ty>
            bool operator==(const disorder_base<other_ty>& other)
            {
                enable_if_convertible(other_ty, _ty);
                assert(this->length == other.get_length() && "Incompatible sizes of disorder array");
                
                arma::uvec elemets_equality = this->disorder_values == other.get_values();
                bool result = this->bound == other.get_bound();
                for(auto& check_element : elemets_equality){ result = result && check_element; }
                
                return result;
            }

            template <typename other_ty>
            bool operator!=(const disorder_base<other_ty>& other)
            { 
                enable_if_convertible(other_ty, _ty); 
                return !(*this == other); 
            }

        //<! ARITHMETICS ------------------------------------------------------------------------------------------------------------
            #define apply_operation(who, other_ty, _ty, arg, op)    \
                        enable_if_convertible(other_ty, _ty);               \
                        who disorder_values = op(who disorder_values, arg); \
                        who bound           = op(who bound, arg);
            //<!------------------------------------------------------------- *
            template <typename _type> _nodiscard friend              
            auto operator*(const disorder_base& other, _type arg)
                {   disorder_base<_ty> new_inst = other;                     
                    apply_operation(new_inst., _type, _ty, arg, operator*); 
                    return std::move(new_inst); }                    
                                                                     
            template <typename _type> _nodiscard friend              
            auto operator*(_type arg, const disorder_base& other)                   
                { return op(other, arg); }                           

            template <typename _type>  _noreturn                                                     
            auto operator*=(_type arg)                                      
                { apply_operation(this->, _type, _ty, arg, operator*); } 

            //<!------------------------------------------------------------- /
            template <typename _type> _nodiscard friend              
            auto operator/(const disorder_base& other, _type arg)
                {   disorder_base<_ty> new_inst = other;                     
                    apply_operation(new_inst., _type, _ty, arg, operator/); 
                    return std::move(new_inst); }                    
                                                                     
            template <typename _type> _nodiscard friend              
            auto operator/(_type arg, const disorder_base& other)                   
                { return op(other, arg); }                           

            template <typename _type>  _noreturn                                                     
            auto operator/=(_type arg)                                      
                { apply_operation(this->, _type, _ty, arg, operator/); } 
            
            //<!------------------------------------------------------------- +
            template <typename _type> _nodiscard friend              
            auto operator+(const disorder_base& other, _type arg)
                {   disorder_base<_ty> new_inst = other;                     
                    apply_operation(new_inst., _type, _ty, arg, operator+); 
                    return std::move(new_inst); }                    
                                                                     
            template <typename _type> _nodiscard friend              
            auto operator+(_type arg, const disorder_base& other)                   
                { return op(other, arg); }                           

            template <typename _type>  _noreturn                                                     
            auto operator+=(_type arg)                                      
                { apply_operation(this->, _type, _ty, arg, operator+); } 
            
            //<!------------------------------------------------------------- -
            template <typename _type> _nodiscard friend              
            auto operator-(const disorder_base& other, _type arg)
                {   disorder_base<_ty> new_inst = other;                     
                    apply_operation(new_inst., _type, _ty, arg, operator-); 
                    return std::move(new_inst); }                    
                                                                     
            template <typename _type> _nodiscard friend              
            auto operator-(_type arg, const disorder_base& other)                   
                { return op(other, arg); }                           

            template <typename _type>  _noreturn                                                     
            auto operator-=(_type arg)                                      
                { apply_operation(this->, _type, _ty, arg, operator-); } 
    };  

    //<! -----------------------------------------------------------------------------------------------------------------------------------------------------
    //<! -----------------------------------------------------------------------------------------------------------------------------------------------------
    //<! ------------------------------------------------------------ UNIFORM DISTRIBUTION
    template <typename _ty>
    class uniform : public disorder_base<_ty>{
    public:
        uniform() = default;
        uniform(int L, _ty W)
            { this->init(L, W); };

        virtual void set() override
            { this->disorder_values = this->generator.template create_random_vec<_ty>(this->length, this->bound); }
    };


    //<! -----------------------------------------------------------------------------------------------------------------------------------------------------
    //<! -----------------------------------------------------------------------------------------------------------------------------------------------------
    //<! ------------------------------------------------------------ QUASIPERIODIC DISORDER
    template <typename _ty>
    class quasiperiodic : public disorder_base<_ty>{
        double phase = 0;
    public:
        virtual void set() override
            { this->disorder_values = this->generator.template create_random_vec<_ty>(this->length, this->bound); }

    };


    //<! -----------------------------------------------------------------------------------------------------------------------------------------------------
    //<! -----------------------------------------------------------------------------------------------------------------------------------------------------
    //<! ------------------------------------------------------------ IMPURITY POTENTIAL
    template <typename _ty>
    class impurity : public disorder_base<_ty>{
        std::vector<int> positions;         //<! positions for impurities
        std::vector<_ty> impurity_vals;     //<! values of impurities only at given positions
        bool add_edge_impurity = true;      //<! add impurity on edge to break mirror symmetry
        double edge_impurity = 0.1;         //<! valu of impurity at edge

        void check_positions()
        {
            assert(this->positions.size() == this->impurity_vals.size() && "ERROR: size of input arrays does not match");
            auto [min, max] = *std::minmax_element(this->positions.begin(), this->positions.end());
            assert( (max < this->length && min >= 0)
                    && "ERROR: Positions of impurity exceed system length");
        }

    public:
        impurity() = default;
        
        /// @brief Main constructor with arbitrary number of impuriteis
        /// @param L system size (may be larger than values in pos)
        /// @param pos positions of impurities on the lattice
        /// @param values strength of impurity at given positions
        /// @param add_edge whether to add additional impurity on edge (mostly useful for single impurity models)
        /// @param edge_value value of additonal impurity on edge
        impurity(int L, std::vector<int> pos, std::vector<_ty> values, 
                    bool add_edge = true, double edge_value = 0.1)
        {
            this->add_edge_impurity = add_edge;
            this->edge_impurity = edge_value;
            this->positions = pos;
            this->impurity_vals = values;
            check_positions();
            this->init(L, *std::max_element(values.begin(), values.end()));
        }
        
        /// @brief constructor with arbitrary number of impurities of equal strength
        /// @param L system size (may be larger than values in pos)
        /// @param pos positions of impurities on the lattice
        /// @param value strength of impurities
        /// @param add_edge whether to add additional impurity on edge (mostly useful for single impurity models)
        /// @param edge_value value of additonal impurity on edge
        impurity(int L, std::vector<int> pos, _ty value, 
                    bool add_edge = true, double edge_value = 0.1)
        {
            this->add_edge_impurity = add_edge;
            this->edge_impurity = edge_value;
            this->positions = pos;
            this->impurity_vals = std::vector<_ty>(pos.size, value);
            check_positions();
            this->init(L, value);
        }
        
        /// @brief constructor for a single impurity
        /// @param L system size (may be larger than values in pos)
        /// @param pos position of impurity on the lattice
        /// @param value strength of impurity at given positions
        /// @param add_edge whether to add additional impurity on edge (mostly useful for single impurity models)
        /// @param edge_value value of additonal impurity on edge
        impurity(int L, int pos, _ty value, 
                    bool add_edge = true, double edge_value = 0.1)
        {
            this->add_edge_impurity = add_edge;
            this->edge_impurity = edge_value;
            this->positions.append(pos);
            this->impurity_vals.append(value);
            check_positions();
            this->init(L, value);
        }

        virtual void set() override
        {
            if(this->add_edge_impurity)
                this->disorder_values(0) = this->edge_impurity;
            int i = 0;
            for(auto& ell : positions)
                this->disorder_values(ell) = this->impurity_vals[i++];
        }
    };
};