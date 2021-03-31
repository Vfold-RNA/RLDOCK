#ifndef FORCEFIELD_H
#define FORCEFIELD_H
#include "main.h"
#include "RLDOCK_PARAMETER.h"

class FORCEFIELD
{
    public:
    vector< vector< vector< double > > >LJ_potential;
    vector< vector<double> > eql_lj;
    vector< vector< vector< double > > >HB_potential;
    vector< vector<double> > eql_hb;
    vector< vector< vector< double > > >SASA_potential;
    vector< vector<double> > eql_sasa;
    
    
    
    void LJ(const PARAMETER& P);
    void HB(const PARAMETER& P);
    double SASA_Vall(const PARAMETER& P, const double& r1, const double& r2, const double& dis);
    void SASA(const PARAMETER& P);
    void Initialization(const PARAMETER& P);
    protected:
    int M;
    int N;
    vector< vector<double> > eq_length;
    
    
};


#endif