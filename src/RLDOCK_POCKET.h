#ifndef POCKET_H
#define POCKET_H
#include "main.h"

class LABEL
{
    public:
    int lig_idx;
    int initial_idx;
    double lj;
};

bool less_label_lj(const LABEL& l1, const LABEL& l2);

class POCKET
{
    public:
    int index;
    double x;
    double y;
    double z;
    double dis2rec;
    double lj_site_min=1000.0;
    vector< vector<int> > chosen_atom_lig;
    vector< vector<LABEL> > LJ_min;
    vector<double> rb;
    
};

bool less_pocket_lj (const POCKET& p1, const POCKET& p2);

#endif