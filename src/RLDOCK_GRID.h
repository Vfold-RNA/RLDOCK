#ifndef GRID_H
#define GRID_H
#include "RLDOCK_CASE.h"
class GRID
{
    public:
    
    GRID(double igd_in, double disaway_in) : igd(igd_in),disaway(disaway_in){}
    
    int Nnx;
    int Nny;
    int Nnz;
    double x0;
    double y0;
    double z0;
    
    void Set_grid_pocket(const CASE& mol);

    void Set_grid_common(const CASE& mol);
    
    double getigd() const {return igd;}
    // int getNnx(){return Nnx;}
    // int getNny(){return Nny;}
    // int getNnz(){return Nnz;}
    // double getx0(){return x0;}
    // double gety0(){return y0;}
    // double getz0(){return z0;}
  
    vector< vector< vector<int> > > grids_pocket_flag;
    void Set_grid_int();
    
    vector< vector< vector< vector<double> > > >grids_energy;    
    vector< vector< vector< vector<int> > > >lj_flag;
    void Initial_grid_lj(const PARAMETER& P);

    void Find_RNA_grid(const PARAMETER& P, const double& select_radius, const ATOM& rec_atom);
    void Set_LJ(const PARAMETER& P, const double& select_radius, const ATOM& rec_atom, const FORCEFIELD& F);

    private:
    double igd;
    double disaway;
    double dis;
    

    int xi;
    int xa;
    int yi;
    int ya;
    int zi;
    int za;
    
};
#endif