#ifndef PARAMETER_H
#define PARAMETER_H
#include "main.h"


class PARAMETER
{
    public:
    int num_threads = 30;
    int num_ligand_atom;

    double pi = 3.1415926;
    int num_atom_type = 9;
    double radius_small_sphere = 1.5;
    double radius_large_sphere = 6.0;
    int LJstep = 1000;
    double LJ_cut = 2.50;
    double equ_LJ = 0.8;
    double excludLJ = 0.6;
    double LJ_energy_cutoff = 30 ;//25.91//pow(equ_LJ/excludLJ,12)-pow(equ_LJ/excludLJ,6);
    double Num_in_sphere = 1000;
    int Nats_RNA = 5500;
    int Nats_ligand =200;
    int Nats_RNA_addH = 11000;
    int Nats_ligand_addH = 400;
    double gamma_SASA = 0.005;
    double r_water = 1.4;
    double step_pocket = 0.5;
    int MaxN = 50;
    int Maxpose = 20000;
    double gamma_step1 = 72.0;//10
    double Max_pocket_site = 10000;
    double rmsd_cut_step1 = 0.5;
    double ligand_cut_step1 = 2.0;
    double Hbond_max = 1.3;
    double Hbond_min = 0.8;
    double step_SASA = 0.25;
    int TOP_step1 = 2;
    int TOP_step2_pose = 3;
    int TOP_step2_site = 300;
    int TOP_step2_conformer = 3;
    int TOP_focus1_num_each_site = 10;
    vector<double> weight_focus1={3.3, 1.32, 0.36, 0.58, 0.30, 0.10, 0.66};
    vector<double> weight_focus2={3.3, 1.32, 1.38, 2.78, 1.26, 0.10, 0.66, 4.98};

    double step_grid_energy = 0.2;
    double LJ_cut_new = 10.0;
    double lj_dr = LJ_cut_new/LJstep;
    double hb_cut_new =5.0;
    double hb_dr = hb_cut_new/LJstep;
    double iseed1 = 25478553;
    double iseed2 = 25212112;
    double Tc = 25.0;
    double e1 = 20.0;
    double radius_O = 1.5;
    double radius_P = 1.9;
    double radius_H = 1.0;
    double radius_C = 1.7;
    double radius_N = 1.65;
    double radius_S = 1.80;
    // vector<double> atom_radius = {radius_H,radius_C,radius_N,radius_O,radius_P,radius_S,radius_C,radius_small_sphere,radius_large_sphere};
    vector<double> atom_radius = {radius_H,radius_C,radius_N,radius_O,radius_P,radius_S,radius_C};
    double born_scale_O = 0.85;
    double born_scale_P = 0.86;
    double born_scale_H = 0.85;
    double born_scale_C = 0.72;
    double born_scale_N = 0.79; 
    double born_scale_S = 0.80;
    double e2 = (87.740-0.4008*Tc+9.398*1e-4*Tc*Tc-1.41*1e-6*Tc*Tc*Tc);   // the dielctric consant of water in temperature T 
    double e25 = 87.740-0.4008*25+9.398*1e-4*25*25-1.41*1e-6*25*25*25;    // the dielctric consant of water in temperature 25
    double lB = 7.15*(273+25)*e25/((273+Tc)*e2);   // e^2/(ebs*kB*Tc) 
    double lB0=e2*lB;        // in A lb0=e^2/(4*pi*e0*kB*T) 
    double rate_kcal_to_kt=0.593*(273+Tc)/(273+25);
    double pol_par = lB0*(1./e2-1/e1);
    double ele_par = lB0/e1;

    vector<double> ux;
    vector<double> vy;
    vector<double> wz;
    void SetSphere(string filename);

    vector<double> aa;
    vector<double> bb;
    vector<double> cc;
    vector<double> dd;
    vector<double> ee;
    vector<double> ff;
    vector<double> gg;
    vector<double> hh;
    vector<double> ii;

    vector< vector<double> > poses;

    void SetAngle();
    void UpdatAngle(int new_gamma_step1);
    void InputUpdate(string rf_in, string lf_in, string ref_in,string op_in, int thread_in);

    string receptor_file;
    string ligand_file;
    string output_prefix;
    string reference_file;
    bool reference_flag;

};

class THREAD_PARAMETER
{
    public:
    int ibegin = 0;
    int iend;
};
class RANK_PARAMETER
{
    public:
    int index;
    double val;
};
bool rank_less(const RANK_PARAMETER& i, const RANK_PARAMETER& j);
void Thread_Distribution(vector<THREAD_PARAMETER>& thread_p, int n_total, int n_piece);
#endif