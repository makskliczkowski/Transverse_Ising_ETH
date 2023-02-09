#pragma once

using site_type = signed long long;

class lattice_base{
protected:
    bool boundary_condition = true;    //<! boundary condition: true = PBC
 public:

    site_type volume;
    virtual ~lattice_base() = 0;
    virtual std::vector<site_type> get_neighbours(site_type site) const = 0;
    
};
inline lattice_base::~lattice_base(){}


//-------------------------------------------------------------------------------------------------------------- 1D LATTICES
//--------------------------------------------------------------------------------------------------------------------------
class lattice1D : public lattice_base{
private:
    site_type L;
    
public:

    //---------------------------------------------- CONSTRUCTORS
    lattice1D(site_type size, bool _bound_cond = false)
    {
        this->L = size;
        this->boundary_condition = !_bound_cond;
        this->volume = this->L;
    };

    //---------------------------------------------- NEIGHBOURS
    virtual std::vector<site_type> get_neighbours(site_type site) 
    const 
    override
    {
        std::vector<site_type> neis;
        
        site_type nei = site + 1;
        if(nei >= this->L){
            if(this->boundary_condition) nei = nei % L;
            else nei = -1;
        }
        neis.push_back(nei);
        
        nei = site - 1;
        if(nei < 0){
            if(this->boundary_condition) {
                nei = (nei % L) + this->L;
            }
            else nei = -1;
        }
        neis.push_back(nei);
        return neis;
    }
};



//-------------------------------------------------------------------------------------------------------------- 2D LATTICES
//--------------------------------------------------------------------------------------------------------------------------

class lattice2D : public lattice_base{
private:
    site_type Lx;
    site_type Ly;

    site_type A;    //<! surface area of system
public:

    //---------------------------------------------- CONSTRUCTORS
    //<! general constructor(rectangles)
    lattice2D(site_type _Lx, site_type _Ly, bool _bound_cond = true)
    {
        this->Lx = _Lx;
        this->Ly = _Ly;
        this->boundary_condition = _bound_cond;
        this->A = this->Lx * this->Ly;
        this->volume = this->A;
    };

    //<! square constructor(rectangles)
    lattice2D(site_type _L, bool _bound_cond = true)
    {
        this->Lx = _L;
        this->Ly = _L;
        this->boundary_condition = _bound_cond;
        this->A = this->Lx * this->Ly;
        this->volume = this->A;
    };

    //---------------------------------------------- NEIGHBOURS
    virtual std::vector<site_type> get_neighbours(site_type site) 
    const 
    override
    {
        std::vector<site_type> neis;
        site_type nei;

        // go right along X
        if(site % this->Lx == this->Lx - 1){
            if(this->boundary_condition) nei = site + 1 - this->Lx;
            else nei = -1;
        }
        else 
            nei = site + 1;
        neis.push_back(nei);
        
        // go left along X
        if(site % this->Lx == 0){
            if(this->boundary_condition) nei = site - 1 + this->Lx;
            else nei = -1;
        }
        else 
            nei = site - 1;
        neis.push_back(nei);
        
        // go up along Y
        if(site >= this->A - this->Lx ){
            if(this->boundary_condition) nei = site % this->Lx;
            else nei = -1;
        }
        else 
            nei = site + this->Lx;
        neis.push_back(nei);
        
        // go down along Y
        if(site < this->Lx){
            if(this->boundary_condition) nei = this->A - this->Lx + site;
            else nei = -1;
        }
        else 
            nei = site - this->Lx;
        neis.push_back(nei);
        return neis;
        
    }
};



//-------------------------------------------------------------------------------------------------------------- 3D LATTICES
//--------------------------------------------------------------------------------------------------------------------------
class lattice3D : public lattice_base{
private:
    site_type Lx;
    site_type Ly;
    site_type Lz;

    site_type V;    //<! volume of system
    site_type A;    //<! surface area of base
public:

    //---------------------------------------------- CONSTRUCTORS
    //<! general constructor(rectangles)
    lattice3D(site_type _Lx, site_type _Ly, site_type _Lz, bool _bound_cond = true)
    {
        this->Lx = _Lx;
        this->Ly = _Ly;
        this->Lz = _Lz;
        this->boundary_condition = _bound_cond;
        this->A = this->Lx * this->Ly;
        this->V = this->Lx * this->Ly * this->Lz;
        this->volume = this->V;
    };
    //<! orthombic constructor(rectangles)
    lattice3D(site_type _Lx, site_type _Lz, bool _bound_cond = true)
    {
        this->Lx = _Lx;
        this->Ly = _Lx;
        this->Lz = _Lz;
        this->boundary_condition = _bound_cond;
        this->A = this->Lx * this->Ly;
        this->V = this->Lx * this->Ly * this->Lz;
        this->volume = this->V;
    };

    //<! cube constructor(rectangles)
    lattice3D(site_type _L, bool _bound_cond = true)
    {
        this->Lx = _L;
        this->Ly = _L;
        this->Lz = _L;
        this->boundary_condition = _bound_cond;
        this->A = this->Lx * this->Ly;
        this->V = this->Lx * this->Ly * this->Lz;
        this->volume = this->V;
    };

    //---------------------------------------------- NEIGHBOURS
    virtual std::vector<site_type> get_neighbours(site_type site) 
    const 
    override
    {
        std::vector<site_type> neis;
        site_type nei;

        // go right along X
        if(site % this->Lx == this->Lx - 1){
            if(this->boundary_condition) nei = site + 1 - this->Lx;
            else nei = -1;
        }
        else 
            nei = site + 1;
        neis.push_back(nei);
        
        // go left along X
        if(site % this->Lx == 0){
            if(this->boundary_condition) nei = site - 1 + this->Lx;
            else nei = -1;
        }
        else 
            nei = site - 1;
        neis.push_back(nei);
        
        // go up along Y
        if(site % this->A >= this->A - this->Lx ){
            if(this->boundary_condition) nei = site - (this->A - this->Lx);
            else nei = -1;
        }
        else 
            nei = site + this->Lx;
        neis.push_back(nei);
        
        // go down along Y
        if(site % this->A < this->Lx){
            if(this->boundary_condition) nei = site + (this->A - this->Lx);
            else nei = -1;
        }
        else 
            nei = site - this->Lx;
        neis.push_back(nei);

        // go inside along Z
        if(site >= this->V - this->A ){
            if(this->boundary_condition) nei = site % this->A;
            else nei = -1;
        }
        else 
            nei = site + this->A;
        neis.push_back(nei);
        
        // go outside along Z
        if(site < this->A){
            if(this->boundary_condition) nei = site + this->V - this->A;
            else nei = -1;
        }
        else 
            nei = site - this->A;
        neis.push_back(nei);


        return neis;
    }
};