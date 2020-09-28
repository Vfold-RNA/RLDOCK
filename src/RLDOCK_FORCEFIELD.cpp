#include "RLDOCK_FORCEFIELD.h"

void FORCEFIELD::Initialization(const PARAMETER& P)
{
       M = P.atom_radius.size();
       N=P.LJstep;
       LJ_potential = vector< vector< vector< double > > >(M,vector< vector<double>>(M,vector<double>(N,0.0)));
       HB_potential = vector< vector< vector< double > > >(M,vector< vector<double>>(M,vector<double>(N,0.0)));
       SASA_potential = vector< vector< vector< double > > >(M,vector< vector<double>>(M,vector<double>(N,0.0)));
       
       eq_length.resize(M,vector<double>(M,0.0));
       eql_lj.resize(M,vector<double>(M,0.0));
       eql_hb.resize(M,vector<double>(M,0.0));
       eql_sasa.resize(M,vector<double>(M,0.0));
       for(int i=0;i!=P.atom_radius.size();i++)
       for(int j=0;j!=P.atom_radius.size();j++)
       {
           eq_length[i][j] = P.atom_radius[i]+P.atom_radius[j];
           eql_lj[i][j] = eq_length[i][j]*P.LJ_cut;
           eql_hb[i][j] = eq_length[i][j]*P.Hbond_max;
           eql_sasa[i][j] = eq_length[i][j]+P.r_water*2.0;
       }
       LJ(P);
       HB(P);
       SASA(P);
       

}

void FORCEFIELD::LJ(const PARAMETER& P)
{
    ////////LJ index start from 1 in calculation 0 can not be divined.
    for(int i=0;i!=P.atom_radius.size();i++)
    {
        for(int j=0;j!=P.atom_radius.size();j++)
        {
            double rdis = P.equ_LJ*(P.atom_radius[i] + P.atom_radius[j]);
            double dr = eq_length[i][j]*P.LJ_cut/P.LJstep;
            for( int istep=0; istep!=P.LJstep;istep++)
            {
                double dis = (istep+1)*dr;
                double dx = rdis/dis;
                double x6 = dx*dx*dx*dx*dx*dx;
                LJ_potential[i][j][istep] = x6*x6 - x6;
            }
        }
    }
}


void FORCEFIELD::HB(const PARAMETER& P)
{
    ////////LJ index start from 1 in calculation 0 can not be divined.
    // HB_potential=vector< vector< vector<double> > >(P.atom_radius.size(),vector< vector<double> >(P.atom_radius.size(),vector<double>(P.LJstep,0.0)));
    for(int i=0;i!=P.atom_radius.size();i++)
    {
        for(int j=0;j!=P.atom_radius.size();j++)
        {
            double dismax = P.Hbond_max*(P.atom_radius[i] + P.atom_radius[j]);
            double dismin = P.Hbond_min*(P.atom_radius[i] + P.atom_radius[j]);
            double dr=dismax/P.LJstep;
            for( int istep=0; istep!=P.LJstep;istep++)
            {
                double dis = (istep+1)*dr;
                if(dis>dismin)
                {
                    HB_potential[i][j][istep] = 1. - (dis-dismin)/(dismax-dismin);
                }else {
                    HB_potential[i][j][istep] =1.;
                }
            }
        }
    }
}


double FORCEFIELD::SASA_Vall(const PARAMETER& P, const double& r1, const double& r2, const double& dis)
{
    double Vall=0.;
    double Rlarge=max(r1,r2);
    double Rsmall=min(r1,r2);
    if(dis<Rlarge+Rsmall)
    {
        if(Rlarge>dis+Rsmall)
        {
            Vall = 4*P.pi*Rsmall*Rsmall;
        }
        else
        {
            double anglelarge=(Rlarge*Rlarge+dis*dis-Rsmall*Rsmall)/(2*Rlarge*dis);
            double Hlarge = Rlarge*(1-anglelarge);
            double anglesmall=(Rsmall*Rsmall+dis*dis-Rlarge*Rlarge)/(2*Rsmall*dis);
            double Hsmall = Rsmall*(1-anglesmall);
            double Vlarge=2*P.pi*Rlarge*Hlarge;
            double Vsmall=2*P.pi*Rsmall*Hsmall;
            Vall=Vlarge+Vsmall;
        }
    }
    
    return Vall;
}
void FORCEFIELD::SASA(const PARAMETER& P)
{
    for(int i=0;i!=P.atom_radius.size();i++)
    {
        for(int j=0;j!=P.atom_radius.size();j++)
        {
            double rwater_i = P.atom_radius[i] + P.r_water ;
            double rwater_j = P.atom_radius[j] + P.r_water ;
            double dismax = rwater_i + rwater_j;
            double dismax2 = P.atom_radius[i] + rwater_j;
            double dismax3 = rwater_i + P.atom_radius[j];
            double dr = dismax/P.LJstep; 
    
            for( int istep=0; istep!=P.LJstep;istep++)
            {
                double dis = (istep+1)*dr;
                double Vall1=SASA_Vall(P,rwater_i,rwater_j,dis);
                double Vall2=SASA_Vall(P,P.atom_radius[i],rwater_j,dis);
                double Vall3=SASA_Vall(P,rwater_i,P.atom_radius[j],dis);

                if(i!=0 && j!=0)
                {
                    SASA_potential[i][j][istep]=P.gamma_SASA*(Vall1-Vall2+Vall1-Vall3)/P.rate_kcal_to_kt;
                }
                else
                {
                    SASA_potential[i][j][istep]=0.;
                }
            }
        }
    }
}
