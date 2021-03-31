#include "RLDOCK_GRID.h"

void GRID::Set_grid_pocket(const CASE& mol){
    this->Nnx = int((mol.x_max() + this->disaway - mol.x_min() + this->disaway)/this->igd + 0.1);
    this->Nny = int((mol.y_max() + this->disaway - mol.y_min() + this->disaway)/this->igd + 0.1);
    this->Nnz = int((mol.z_max() + this->disaway - mol.z_min() + this->disaway)/this->igd + 0.1);

    this->x0 = mol.x_min() - this->disaway + this->igd;
    this->y0 = mol.y_min() - this->disaway + this->igd;
    this->z0 = mol.z_min() - this->disaway + this->igd;    
}
void GRID::Set_grid_common(const CASE& mol){
    this->Nnx = int((mol.x_max() + this->disaway - mol.x_min() + this->disaway)/this->igd + 0.1)+2;
    this->Nny = int((mol.y_max() + this->disaway - mol.y_min() + this->disaway)/this->igd + 0.1)+2;
    this->Nnz = int((mol.z_max() + this->disaway - mol.z_min() + this->disaway)/this->igd + 0.1)+2;

    this->x0 = mol.x_min() - this->disaway;
    this->y0 = mol.y_min() - this->disaway;
    this->z0 = mol.z_min() - this->disaway;    
}


void GRID::Set_grid_int(){
    this->grids_pocket_flag=vector<vector<vector< int > > >(this->Nnx, vector<vector< int > >(this->Nny, vector< int >(this->Nnz,0)));
}
void GRID::Initial_grid_lj(const PARAMETER& P)
{
    for(int i=0;i!=P.atom_radius.size();i++)
    {
        this->grids_energy.push_back(vector<vector<vector< double > > >(this->Nnx, vector<vector< double > >(this->Nny, vector< double >(this->Nnz,0))) );
        this->lj_flag.push_back(vector<vector<vector< int > > >(this->Nnx, vector<vector< int > >(this->Nny, vector< int >(this->Nnz,0))) );
    }
}
void GRID::Find_RNA_grid(const PARAMETER& P, const double& select_radius, const ATOM& rec_atom)
{
    
    this->xi = (rec_atom.x - this->x0 - select_radius)/this->igd;     if(this->xi<0)this->xi=0;
    this->xa = (rec_atom.x - this->x0 + select_radius)/this->igd + 1; if(this->xa>=this->Nnx)this->xa=this->Nnx-1;
    this->yi = (rec_atom.y - this->y0 - select_radius)/this->igd;     if(this->yi<0)this->yi=0;
    this->ya = (rec_atom.y - this->y0 + select_radius)/this->igd + 1; if(this->ya>=this->Nny)this->ya=this->Nny-1;
    this->zi = (rec_atom.z - this->z0 - select_radius)/this->igd;     if(this->zi<0)this->zi=0;
    this->za = (rec_atom.z - this->z0 + select_radius)/this->igd + 1; if(this->za>=this->Nnz)this->za=this->Nnz-1;
    for(int i=this->xi; i<=this->xa; i++)
    for(int j=this->yi; j<=this->ya; j++)
    for(int k=this->zi; k<=this->za; k++)
    {
        if(this->grids_pocket_flag[i][j][k]==0)
        {
            double xtmp = P.step_pocket*i+this->x0 - rec_atom.x;
            dis  = xtmp*xtmp;
            double ytmp = P.step_pocket*j+this->y0 - rec_atom.y;
            dis += ytmp*ytmp;
            double ztmp = P.step_pocket*k+this->z0 - rec_atom.z;
            dis += ztmp*ztmp;
            dis = sqrt(dis);
            if( dis < rec_atom.rrr + P.radius_small_sphere) this->grids_pocket_flag[i][j][k]=1;
        }
    }
}
void GRID::Set_LJ(const PARAMETER& P, const double& select_radius, const ATOM& rec_atom, const FORCEFIELD& F)
{
    this->xi = (rec_atom.x - this->x0 - select_radius)/this->igd;   //if(xi<0)xi=0;
    this->xa = (rec_atom.x - this->x0 + select_radius)/this->igd+1; //if(xa>=Nnx)xa=Nnx-1;
    this->yi = (rec_atom.y - this->y0 - select_radius)/this->igd;   //if(yi<0)yi=0;
    this->ya = (rec_atom.y - this->y0 + select_radius)/this->igd+1; //if(ya>=Nny)ya=Nny-1;
    this->zi = (rec_atom.z - this->z0 - select_radius)/this->igd;   //if(zi<0)zi=0;
    this->za = (rec_atom.z - this->z0 + select_radius)/this->igd+1; //if(za>=Nnz)za=Nnz-1;
    double disx,disy,disz;
    double xtmp = this->igd*(this->xi-1)+this->x0 - rec_atom.x;
    for(int i=this->xi; i<=this->xa; i++)
    {
        xtmp += this->igd;
        disx  = xtmp * xtmp;
        double ytmp = this->igd*(this->yi-1)+this->y0 - rec_atom.y;
        for(int j=this->yi; j<=this->ya; j++)
        {
            ytmp += this->igd;
            disy = ytmp * ytmp;
            double ztmp = this->igd*(this->zi-1)+this->z0-rec_atom.z;
            for(int k=this->zi; k<=this->za; k++)
            {     
                ztmp += this->igd;
                dis = disx+disy+ztmp*ztmp;
                dis = sqrt(dis);
                // int dis_idx = dis/P.lj_dr;
                for(int l=0;l!=P.atom_radius.size();l++)
                {
                    int recidx=rec_atom.type_idx;
                    //P.excludLJ*(P.atom_radius[l]+rec_atom.rrr)
                    if(dis < F.eql_lj[l][recidx])
                    {
                        int lj_idx = dis/F.eql_lj[l][recidx]*P.LJstep;
                        this->grids_energy[l][i][j][k] += F.LJ_potential[l][recidx][lj_idx];
                    }
                }
            }
        }
    }
}