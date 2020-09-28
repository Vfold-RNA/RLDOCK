#ifndef CASE_H
#define CASE_H
#include "RLDOCK_FORCEFIELD.h"
#include "RLDOCK_COORDINATE.h"


class CASE
{
    public:
    std::vector< ATOM > noH;
    std::vector< ATOM > wH;
    std::vector< ATOM > HH;
   

    std::vector< ATOM > center;
    void Move_to_Center_Easy();
    
    void Find_Max_and_Min();
    double x_max() const {return xmax;}
    double y_max() const {return ymax;}
    double z_max() const {return zmax;}
    double x_min() const {return xmin;}
    double y_min() const {return ymin;}
    double z_min() const {return zmin;}

    vector<vector<double> >distance;
    void cal_distance();

    double pol;
    void get_rb0_and_cal_pol(const PARAMETER & P);
    
    //////For receptor
    vector< vector< vector<COORDINATE> > > radius_rw;
    void set_rw_position(const PARAMETER & P);

    //////For ligand conformer
    double score;
    double internal_lj=0;
    
    vector<vector<int> > usepose;
    void cal_internal_lj(const FORCEFIELD & F, const PARAMETER& P);
     
    private:
    double xcen=0;
    double ycen=0;
    double zcen=0;
    double xmax;
    double ymax;
    double zmax;
    double xmin;
    double ymin;
    double zmin;
    
    vector<vector<double> >distance_square;

};

class COMPLEX
{
    public:
    vector<CASE> rec;
    vector<CASE> conformers;
    vector<CASE> lig;//reference

    int lig_atom_num;
    int rec_atom_num;
    int conformers_num;

    void NumInput();
    void Initialization(const PARAMETER& P, const FORCEFIELD& F);
};
void READFILE(const PARAMETER P, vector<CASE>& in_mol, string filename);
#endif