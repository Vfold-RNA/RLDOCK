#ifndef POSE_H
#define POSE_H
#include "main.h"
#include "RLDOCK_COORDINATE.h"
#include "RLDOCK_POCKET.h"
class POSE
{
    public:
    int num_cluster=1;
    int iconf;
    int iatom;
    int isite;
    int ipose;
    double rmsd;
    double rmsd_show;
    double lj=0;
    double ele=0;
    double pol=0;
    double self=0;
    double self_lig=0;
    double self_rec=0;
    double sasa=0;
    double hb=0;
    double internal_lj=0;
    double eng;
    double eng_rerank;
    vector<COORDINATE> position;
    void calculate_energy(const vector<double>& w){
        eng = lj*w[0]+ele*w[1]+pol*w[2]+self*w[3]+sasa*w[4]+hb*w[5]+internal_lj*w[6];
    }
     void calculate_eng_rerank(const vector<double>& w){
        eng_rerank = lj*w[0]+ele*w[1]+pol*w[2]+self_lig*w[3]+sasa*w[4]+hb*w[5]+internal_lj*w[6]+self_rec*w[7];
    }
    void info_input(int iconf_in, int iatom_in, POCKET& pocket, vector<COORDINATE>& pose);

    
};

bool greater_cluster(const POSE& i, const POSE& j);

bool less_eng (const POSE & p1, const POSE& p2);

bool less_eng_rerank (const POSE & p1, const POSE & p2);
#endif