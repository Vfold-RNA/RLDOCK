#include "main.h"
#include "RLDOCK_PARAMETER.h"
#include "RLDOCK_COORDINATE.h"
#include "RLDOCK_CASE.h"
#include "RLDOCK_POCKET.h"
#include "RLDOCK_POSE.h"
#include "RLDOCK_FORCEFIELD.h"
#include "RLDOCK_GRID.h"


void Find_Pocket(const PARAMETER& P, const GRID& gd, const vector<ATOM>& rec_atom,const THREAD_PARAMETER& thread_p,vector<POCKET>& output)
{
    for (int i = thread_p.ibegin; i<=thread_p.iend;i++)
    for (int j = 0; j != gd.Nny; j++)
    for (int k = 0; k != gd.Nnz; k++)
    {
        double xx = i*P.step_pocket+gd.x0;
        double yy = j*P.step_pocket+gd.y0;
        double zz = k*P.step_pocket+gd.z0;
        int totalblock,blockrx=0,blocklx=0,blockry=0,blockly=0,blockrz=0,blocklz=0;
        if(gd.grids_pocket_flag[i][j][k]==0){
            for(auto&a:rec_atom){
                double disxx = (xx - a.x)*(xx - a.x);
                double disyy = (yy - a.y)*(yy - a.y);
                double diszz = (zz - a.z)*(zz - a.z);

                if(blockrx==0){
                    if(a.x < xx && xx< a.x + P.radius_large_sphere){
                        double dd = sqrt(disyy+diszz);
                        if(dd < a.rrr + P.radius_small_sphere)blockrx=1;
                    }
                }

                if(blocklx==0){
                    if(a.x > xx && xx> a.x - P.radius_large_sphere){
                        double dd = sqrt(disyy+diszz);
                        if(dd < a.rrr + P.radius_small_sphere)blocklx=1;
                    }
                }

                if(blockry==0){
                    if(a.y < yy && yy< a.y + P.radius_large_sphere){
                        double dd = sqrt(disxx+diszz);
                        if(dd < a.rrr + P.radius_small_sphere)blockry=1;
                    }
                }

                if(blockly==0){
                    if(a.y > yy && yy> a.y - P.radius_large_sphere){
                        double dd = sqrt(disxx+diszz);
                        if(dd < a.rrr + P.radius_small_sphere)blockly=1;
                    }
                }

                if(blockrz==0){
                    if(a.z < zz && zz< a.z + P.radius_large_sphere){
                        double dd = sqrt(disyy+disxx);
                        if(dd < a.rrr + P.radius_small_sphere)blockrz=1;
                    }
                }

                if(blocklz==0){
                    if(a.z > zz && zz> a.z - P.radius_large_sphere){
                        double dd = sqrt(disyy+disxx);
                        if(dd < a.rrr + P.radius_small_sphere)blocklz=1;
                    }
                }
                totalblock=blocklx+blockrx+blockly+blockry+blocklz+blockrz;
                if(totalblock==6)break;
            }
             if(totalblock==6){
                POCKET p_tmp;
                p_tmp.x = xx;
                p_tmp.y = yy;
                p_tmp.z = zz;
                output.push_back(p_tmp);
                if(output.size()>P.Max_pocket_site){std::cout<<"Too much pocket sites"<<endl;exit(2);}
             }

        }

    }
}

void Select_Pose(const PARAMETER& P,const THREAD_PARAMETER& thread_p, const vector<COORDINATE>centers, const COMPLEX& cmpl,const GRID& gd, vector<POCKET>& pockets, vector< vector<int> >& chose_pose_idx)
{
    
    for(int ic = thread_p.ibegin; ic <= thread_p.iend;ic++)
    {
        
        int i_conf = ic/cmpl.lig_atom_num;
        int i_atom = ic % cmpl.lig_atom_num; 
       
        int num_pose = P.aa.size();
        COORDINATE void_coordinate;
        void_coordinate.index=0;
        void_coordinate.x=0;
        void_coordinate.y=0;
        void_coordinate.z=0;

        vector< vector<COORDINATE> > usepose;
        vector<COORDINATE> pp(vector<COORDINATE>(cmpl.lig_atom_num,void_coordinate));
        for(int i=0;i!=num_pose;i++)
        {
            for (int j=0;j!=cmpl.lig_atom_num;j++)
            {
                double xtmp = cmpl.conformers[i_conf].center[j].x - centers[ic].x; 
                double ytmp = cmpl.conformers[i_conf].center[j].y - centers[ic].y;
                double ztmp = cmpl.conformers[i_conf].center[j].z - centers[ic].z;
                pp[j].x = P.aa[i]*xtmp + P.bb[i]*ytmp + P.cc[i]*ztmp;
                pp[j].y = P.dd[i]*xtmp + P.ee[i]*ytmp + P.ff[i]*ztmp;
                pp[j].z = P.gg[i]*xtmp + P.hh[i]*ytmp + P.ii[i]*ztmp;
            }
            bool use=true;
            for(int iuse=0;iuse!=usepose.size();iuse++)
            {
                double rmsd = 0;
                for(int j=0;j!=cmpl.lig_atom_num;j++)
                {
                    rmsd += (pp[j].x - usepose[iuse][j].x) * (pp[j].x - usepose[iuse][j].x);
                    rmsd += (pp[j].y - usepose[iuse][j].y) * (pp[j].y - usepose[iuse][j].y);
                    rmsd += (pp[j].z - usepose[iuse][j].z) * (pp[j].z - usepose[iuse][j].z);
                    
                }
                rmsd =sqrt(rmsd/cmpl.lig_atom_num);
                if(rmsd<P.rmsd_cut_step1){use = false;break;}
            }
            if(use){
                usepose.push_back(pp);
                chose_pose_idx[ic].push_back(i);
            }
        }
        for(int isite = 0;isite!=pockets.size();isite++)
        {
            for(int ipose=0;ipose!=usepose.size();ipose++)
            {
                bool record=true;
                double ljenergy = 0.;
                for(int ilig =0; ilig!=cmpl.lig_atom_num;ilig++)
                {
                    double x = pockets[isite].x + usepose[ipose][ilig].x - gd.x0;
                    double y = pockets[isite].y + usepose[ipose][ilig].y - gd.y0;
                    double z = pockets[isite].z + usepose[ipose][ilig].z - gd.z0;

                    int xidx = (x + gd.getigd()*0.5)/gd.getigd();
                    int yidx = (y + gd.getigd()*0.5)/gd.getigd();
                    int zidx = (z + gd.getigd()*0.5)/gd.getigd();

                    
                    int atom_type_idx =  cmpl.conformers[0].noH[ilig].type_idx;
                    double var;
                    if(xidx>=gd.Nnx || yidx>=gd.Nny || zidx>=gd.Nnz || xidx<0 || yidx<0 || zidx<0)var=0;
                    else var = gd.grids_energy[atom_type_idx][xidx][yidx][zidx];
                    if(var>P.LJ_energy_cutoff){record=false;break;}
                    ljenergy += var;
                    if(ljenergy >50.0){record=false;break;}
                }
                if(record){
                    if(ljenergy < pockets[isite].LJ_min[i_conf][i_atom].lj)
                    {
                        pockets[isite].LJ_min[i_conf][i_atom].lj = ljenergy;
                        if(i_conf==0 && ljenergy < pockets[isite].lj_site_min)
                        {
                            pockets[isite].lj_site_min= ljenergy;
                        }
                    } 
                }                
            }
        }
    }
}

void BornRadius(const PARAMETER& P, vector<POCKET>& pockets, CASE& rec)
{
    /////for pockets
    for(int ip = 0; ip!=pockets.size(); ip++)
    {
        pockets[ip].rb.resize(P.atom_radius.size());
        for(int iat=0;iat!=P.atom_radius.size();iat++)
        {
            double intb=0, radi = P.atom_radius[iat];
            for (int ir=0;ir!=rec.noH.size();ir++)
            {
                double dis = (pockets[ip].x - rec.noH[ir].x)*(pockets[ip].x - rec.noH[ir].x);
                dis +=  (pockets[ip].y - rec.noH[ir].y)*(pockets[ip].y - rec.noH[ir].y);
                dis +=  (pockets[ip].z - rec.noH[ir].z)*(pockets[ip].z - rec.noH[ir].z);
                dis = sqrt(dis);
                double intbb =0.;
                // if(dis<P.ele_cut_new)
                // {
                    double radj = rec.noH[ir].rrr;
                    double radj_sb = rec.noH[ir].sb*radj;
                    double radc = dis + radj_sb;
                    double Lij,Uij;
                    if(radi>= radc){Lij=1;Uij=1.;}
                    else
                    {
                        if(radi>(dis-radj_sb)){Lij=radi;}
                        else{Lij = dis - radj_sb;}
                        Uij = dis + radj_sb;
                    }
                    intbb = 0.5*((1./Lij-1./Uij)+(radj_sb*radj_sb/(4.*dis)-dis/4.)*(1./(Lij*Lij)-1./(Uij*Uij))+(1./(2.*dis))*log(Lij/Uij));
                // }
                intb += intbb;
            }
            pockets[ip].rb[iat] = 1./(1./radi-intb);
        }
    }

    ///////
    for(int ir = 0; ir !=rec.noH.size(); ir++)
    {
        double intb=0., radi = rec.noH[ir].rrr;
        for(int jr = 0; jr !=rec.noH.size(); jr++)
        {
            if(jr!=ir)
            {
                double dis=(rec.noH[ir].x - rec.noH[jr].x)*(rec.noH[ir].x - rec.noH[jr].x);
                dis += (rec.noH[ir].y - rec.noH[jr].y)*(rec.noH[ir].y - rec.noH[jr].y);
                dis += (rec.noH[ir].z - rec.noH[jr].z)*(rec.noH[ir].z - rec.noH[jr].z);
                dis = sqrt(dis);
                double radj = rec.noH[jr].rrr;
                double radj_sb = rec.noH[jr].sb*rec.noH[jr].rrr;
                double radc = dis+radj_sb;
                double Lij,Uij;
                if(radi>=radc){Lij=0;Lij=0;}
                else{
                    if(radi>(dis-radj_sb)){Lij=radi;}
                    else{Lij=dis-radj_sb;}
                    Uij=dis+radj_sb;
                }
                intb += 0.5*((1./Lij-1./Uij)+(radj_sb*radj_sb/(4.*dis)-dis/4.)*(1./(Lij*Lij)-1./(Uij*Uij))+(1./(2.*dis))*log(Lij/Uij));
            }
        }
        rec.noH[ir].rb0 = 1./(1./radi-intb);
    }


}

void SF_low(const PARAMETER& P,const THREAD_PARAMETER& thread_p,const COMPLEX& cmpl, vector<POCKET>& pockets,const FORCEFIELD & fc, vector<vector<POSE> >& top_each_site)
{
    // std::cout<<thread_p.ibegin<<" "<<thread_p.iend<<endl;
    
    int num_ligand_atom = cmpl.lig_atom_num;
    int num_rec_atom = cmpl.rec_atom_num;
    int num_conformers = cmpl.conformers_num;
    int num_pose = P.aa.size();
    const CASE & rec = cmpl.rec[0];
    const CASE & lig = cmpl.conformers[0];
    const CASE & ref_lig = cmpl.lig[0];
    COORDINATE void_coordinate;
    void_coordinate.index=0;
    void_coordinate.x=0;
    void_coordinate.y=0;
    void_coordinate.z=0;

    ////////////////
    vector< vector< vector<int> > > pockets_needed(num_conformers ,vector< vector<int> >(num_ligand_atom));
    for(int ip=thread_p.ibegin;ip<=thread_p.iend;ip++)
    {
        // if(pockets[ip].dis2rec<2.0)
        // {
            for(int i_conf=0; i_conf !=num_conformers;i_conf++)
            {
                if(pockets[ip].chosen_atom_lig[i_conf].size()>0)
                {
                    for(int id=0; id !=pockets[ip].chosen_atom_lig[i_conf].size();id++)
                    { 
                        int i_atom = pockets[ip].chosen_atom_lig[i_conf][id];
                        // std::cout<<ip<<" "<<i_conf<<" "<<id<<" "<<i_atom<<" "<<pockets[ip].chosen_atom_lig.size()<<endl;
                        pockets_needed[i_conf][i_atom].push_back(ip);
                    }
                }
            }
        // }
    }
    for(int i_conf=0; i_conf!=num_conformers; i_conf++)
    for(int i_atom=0; i_atom!=num_ligand_atom; i_atom++)
    {
        if(pockets_needed[i_conf][i_atom].size()>0)
        {
            // time_t t00=time(NULL);
            // cout<<i_conf<<" "<<i_atom<<endl;
            vector< vector<COORDINATE> > usepose(vector< vector<COORDINATE> >(cmpl.conformers[i_conf].usepose[i_atom].size(),vector<COORDINATE>(num_ligand_atom)));
            
            //////////////selecting useful poses//////////////
            for(int ip=0;ip!=cmpl.conformers[i_conf].usepose[i_atom].size();ip++)
            {
                int ipose=cmpl.conformers[i_conf].usepose[i_atom][ip];
                for (int j=0;j!=num_ligand_atom;j++)
                {
                    double xtmp = cmpl.conformers[i_conf].center[j].x - cmpl.conformers[i_conf].center[i_atom].x; 
                    double ytmp = cmpl.conformers[i_conf].center[j].y - cmpl.conformers[i_conf].center[i_atom].y;
                    double ztmp = cmpl.conformers[i_conf].center[j].z - cmpl.conformers[i_conf].center[i_atom].z;
                    usepose[ip][j].x = P.aa[ipose]*xtmp + P.bb[ipose]*ytmp + P.cc[ipose]*ztmp;
                    usepose[ip][j].y = P.dd[ipose]*xtmp + P.ee[ipose]*ytmp + P.ff[ipose]*ztmp;
                    usepose[ip][j].z = P.gg[ipose]*xtmp + P.hh[ipose]*ytmp + P.ii[ipose]*ztmp;
                    usepose[ip][j].index = ip;
                }

            }
            
            // time_t t11=time(NULL);
            // t1 += t11-t00;
            /////////////////////////////////////////////////////////
            for(int ip=0;ip!=pockets_needed[i_conf][i_atom].size();ip++)
            {
                int isite = pockets_needed[i_conf][i_atom][ip];
                for(int ipose=0;ipose!=usepose.size();ipose++)
                {
                    // time_t t22=time(NULL);
                    bool looprun=true;
                    vector<COORDINATE> tmp_pose(num_ligand_atom);
                    double lj_rough=0;
                    for(int ilig =0; ilig!=num_ligand_atom && looprun;ilig++)
                    {
                        tmp_pose[ilig].x = pockets[isite].x + usepose[ipose][ilig].x;
                        tmp_pose[ilig].y = pockets[isite].y + usepose[ipose][ilig].y;
                        tmp_pose[ilig].z = pockets[isite].z + usepose[ipose][ilig].z;
                    }
                    // time_t t33=time(NULL);
                    // t2 += t33-t22;
                    vector< vector<double> >distance(num_ligand_atom,vector<double>(num_rec_atom,0.0));
                    double rmsd=0;
                    double lj_tmp=0;
                    for(int ilig =0; ilig!=num_ligand_atom && looprun;ilig++)
                    {
                        double x = pockets[isite].x + usepose[ipose][ilig].x;
                        double y = pockets[isite].y + usepose[ipose][ilig].y;
                        double z = pockets[isite].z + usepose[ipose][ilig].z;

                        if(P.reference_flag)
                        {
                            rmsd+=(tmp_pose[ilig].x-ref_lig.noH[ilig].x)*(tmp_pose[ilig].x-ref_lig.noH[ilig].x)
                                +(tmp_pose[ilig].y-ref_lig.noH[ilig].y)*(tmp_pose[ilig].y-ref_lig.noH[ilig].y)
                                +(tmp_pose[ilig].z-ref_lig.noH[ilig].z)*(tmp_pose[ilig].z-ref_lig.noH[ilig].z);
                        }
                        for(int ir=0;ir!= num_rec_atom && looprun;ir++)
                        {
                            int ligidx =lig.noH[ilig].type_idx;
                            int recidx =rec.noH[ir].type_idx;
                            double dd = (tmp_pose[ilig].x-rec.noH[ir].x)*(tmp_pose[ilig].x-rec.noH[ir].x);
                            dd += (tmp_pose[ilig].y-rec.noH[ir].y)*(tmp_pose[ilig].y-rec.noH[ir].y);
                            dd += (tmp_pose[ilig].z-rec.noH[ir].z)*(tmp_pose[ilig].z-rec.noH[ir].z);
                            dd = sqrt(dd);
                            distance[ilig][ir]=dd;

                            if(dd<fc.eql_lj[ligidx][recidx])
                            {
                                int lj_idx = dd/fc.eql_lj[ligidx][recidx]*P.LJstep;
                                lj_tmp += fc.LJ_potential[ligidx][recidx][lj_idx];
                            }
                            if(lj_tmp>10.0){looprun=false;break;}
                        }
                    }

                    POSE pe;
                    pe.rmsd = sqrt(rmsd/num_ligand_atom);
                    if(looprun)pe.lj=lj_tmp;
                    for(int ilig =0; ilig!=num_ligand_atom && looprun;ilig++)
                    {
                        const ATOM& alig = lig.noH[ilig];
                        // cout<<distance[ilig].size()<<endl;
                        for(int id=0;id!=distance[ilig].size() && looprun;id++)
                        {
                            const ATOM& arec = rec.noH[id];
                            pe.pol += P.pol_par*alig.q*arec.q/sqrt(distance[ilig][id]*distance[ilig][id]+pockets[isite].rb[alig.type_idx]*arec.rb0*exp(-distance[ilig][id]*distance[ilig][id]/(4.0*pockets[isite].rb[alig.type_idx]*arec.rb0)));
                            pe.ele += P.ele_par * alig.q * arec.q/distance[ilig][id];
                            
                            if(distance[ilig][id]<fc.eql_sasa[alig.type_idx][arec.type_idx])
                            {
                                int sasa_idx = distance[ilig][id]/fc.eql_sasa[alig.type_idx][arec.type_idx]*P.LJstep;
                                pe.sasa += fc.SASA_potential[alig.type_idx][arec.type_idx][sasa_idx];
                            }
                            
                            if(alig.carryH || arec.carryH)
                            {   
                                if(distance[ilig][id]<fc.eql_hb[alig.type_idx][arec.type_idx])
                                {
                                    int hb_idx = distance[ilig][id]/fc.eql_hb[alig.type_idx][arec.type_idx]*P.LJstep;
                                    pe.hb -= fc.HB_potential[alig.type_idx][arec.type_idx][hb_idx];
                                }
                            }
                        }
                        pe.self += 0.5*P.pol_par*alig.q*alig.q*(1./pockets[isite].rb[alig.type_idx]-1./alig.rrr);
                       
                    }
                    if(looprun){
                        pe.internal_lj = cmpl.conformers[i_conf].internal_lj;
                        pe.ipose = ipose;
                        pe.isite = isite;
                        pe.calculate_energy(P.weight_focus1);
                    }
                    if(looprun && (pe.eng<top_each_site[isite][P.TOP_focus1_num_each_site].eng))
                    {
                        pe.info_input(i_conf,i_atom,pockets[isite],usepose[ipose]);
                        top_each_site[isite][P.TOP_focus1_num_each_site]=pe;
                        for(int irank=P.TOP_focus1_num_each_site;irank!=0;irank--)
                        {
                            if(top_each_site[isite][irank].eng < top_each_site[isite][irank-1].eng)
                            {swap(top_each_site[isite][irank],top_each_site[isite][irank-1]);}
                            else{break;}
                        }

                    }
                }
            }
        } 
    }
}

void SF_high(const PARAMETER& P,const THREAD_PARAMETER& thread_p, const COMPLEX& cmpl,  vector<POCKET>& pockets,const FORCEFIELD & fc, vector<vector<POSE> >& top_each_site)
{
    // std::cout<<thread_p.ibegin<<" "<<thread_p.iend<<endl;
    int num_ligand_atom = cmpl.lig_atom_num;
    int num_pose = P.aa.size();
    const CASE & rec = cmpl.rec[0];
    const CASE & lig = cmpl.conformers[0];
    COORDINATE void_coordinate;
    void_coordinate.index=0;
    void_coordinate.x=0;
    void_coordinate.y=0;
    void_coordinate.z=0;

    double time_pol=0, time_sasa=0;
    
    // for(auto&isite:top_each_site)
    for(int isite=thread_p.ibegin;isite<=thread_p.iend;isite++)
    for(auto&irank:top_each_site[isite])
    {
        if(irank.eng<10000)
        {
            // cout<<irank.isite<<" "<<irank.iconf<<" "<<irank.iatom<<endl;
            // time_t t_pol;
            // t_pol = time(NULL);
            int iconf=irank.iconf;
            int num_ligand_atom = irank.position.size(), num_receptor_atom=cmpl.rec_atom_num;
            vector<double>rb1_lig(num_ligand_atom,0.);
            vector<double>rb1_rec(num_receptor_atom,0);

            vector<vector<double> > distance(num_receptor_atom,vector<double>(num_ligand_atom,0.0));
            for(int ir=0;ir!=num_receptor_atom;ir++)
            for(int il=0;il!=num_ligand_atom;il++)
            {
                distance[ir][il]=sqrt((rec.noH[ir].x-irank.position[il].x)*(rec.noH[ir].x-irank.position[il].x)+(rec.noH[ir].y-irank.position[il].y)*(rec.noH[ir].y-irank.position[il].y)+(rec.noH[ir].z-irank.position[il].z)*(rec.noH[ir].z-irank.position[il].z));
            }
            //////rbi lig
            for(int il=0;il!=num_ligand_atom;il++)
            {
                double radi=lig.noH[il].rrr;
                // double intb=1./radi - 1./conformers[iconf].coordinate[il].rb0;
                double intb = 0;
                for(int ir=0;ir!=num_receptor_atom;ir++)
                {
                    double radj = rec.noH[ir].rrr, Lij=1, Uij=1;
                    double sb_radj = rec.noH[ir].sb*radj;
                    // if(radi>=(distance[ir][il]+sb_radj)){Lij=1.; Uij=1.;}
                    // else
                     if(radi<(distance[ir][il]+sb_radj))
                    {
                        if(radi>(distance[ir][il]-sb_radj)){Lij=radi;}
                        else{Lij = distance[ir][il]-sb_radj ;}
                        Uij=distance[ir][il]+sb_radj;

                         intb+=0.5*((1./Lij-1./Uij)+(sb_radj*sb_radj/(4.*distance[ir][il])-distance[ir][il]/4.)*(1./(Lij*Lij)-1./(Uij*Uij))+(1./(2.*distance[ir][il]))*log(Lij/Uij));   
                    }
                }
                // rb1_lig[il]=1./(1./radi-intb);
                rb1_lig[il]=1./(1./cmpl.conformers[iconf].noH[il].rb0-intb);
            }
            //////rbi receptor
            for(int ir=0;ir!=num_receptor_atom;ir++)
            {
                double radi=rec.noH[ir].rrr;
                // double intb=1./radi - 1./rec.noH[ir].rb0;
                double intb=0.;
                for(int il=0;il!=num_ligand_atom;il++)
                {
                    double radj = lig.noH[il].rrr, Lij=1, Uij=1;
                    double sb_radj = lig.noH[il].sb*radj;
                    // if(radi>=(distance[ir][il]+sb_radj)){Lij=1.; Uij=1.;}
                    // else
                    if(radi<(distance[ir][il]+sb_radj))
                    {
                        if(radi>(distance[ir][il]-sb_radj)){Lij=radi;}
                        else{Lij = distance[ir][il]-sb_radj ;}
                        Uij=distance[ir][il]+sb_radj;

                        intb+=0.5*((1./Lij-1./Uij)+(sb_radj*sb_radj/(4.*distance[ir][il])-distance[ir][il]/4.)*(1./(Lij*Lij)-1./(Uij*Uij))+(1./(2.*distance[ir][il]))*log(Lij/Uij));    
                    }
                    
                }
                // rb1_rec[ir]=1./(1./radi-intb);
                rb1_rec[ir]=1./(1./rec.noH[ir].rb0-intb);
            }
            /////POL NEW
            double pol_new_lig=0,pol_new_rec=0,pol_new_lig_rec=0;
            for(int il=0;il!=num_ligand_atom;il++)
            {
                for(int jl=0;jl!=il;jl++)
                {
                    double dis=cmpl.conformers[iconf].distance[il][jl];
                    double rr=rb1_lig[il]*rb1_lig[jl];
                    pol_new_lig+=P.pol_par*lig.noH[il].q*lig.noH[jl].q/sqrt(dis*dis+rr*exp(-dis*dis/(4.*rr)));
                }
                irank.self_lig += 0.5*P.pol_par*lig.noH[il].q*lig.noH[il].q*(1./rb1_lig[il]-1./cmpl.conformers[iconf].noH[il].rb0);
            }

            for(int ir=0;ir!=num_receptor_atom;ir++)
            {
                for(int jr=0;jr!=ir;jr++)
                {
                    double dis=rec.distance[ir][jr];
                    double rr=rb1_rec[ir]*rb1_rec[jr];
                    pol_new_rec+=P.pol_par*rec.noH[ir].q*rec.noH[jr].q/sqrt(dis*dis+rr*exp(-dis*dis/(4.*rr)));
                }
                irank.self_rec += 0.5*P.pol_par*rec.noH[ir].q*rec.noH[ir].q*(1./rb1_rec[ir]-1./rec.noH[ir].rb0);
            }
            for(int ir=0;ir!=num_receptor_atom;ir++)
            {
                for(int il=0;il!=num_ligand_atom;il++)
                {
                    double dis=distance[ir][il];
                    double rr=rb1_rec[ir]*rb1_lig[il];
                    pol_new_lig_rec+=P.pol_par*rec.noH[ir].q*lig.noH[il].q/sqrt(dis*dis+rr*exp(-dis*dis/(4.*rr)));
                }
            }
            irank.pol=pol_new_lig+pol_new_rec+pol_new_lig_rec-rec.pol-cmpl.conformers[iconf].pol;
            // time_pol+= time(NULL)-t_pol;
            // time_t t_sasa;
            // t_sasa = time(NULL);
            // cout<<"new pol"<<endl;
            ////near list creation
            vector<vector<int> >near_rec(num_receptor_atom);
            vector<vector<int> >near_lig(num_ligand_atom);
            for(int ir=0;ir!=num_receptor_atom;ir++)
            {
                for(int il=0;il!=num_ligand_atom;il++)
                {
                    if(distance[ir][il]<(2*P.r_water + rec.noH[ir].rrr + lig.noH[il].rrr))
                    {
                        near_rec[ir].push_back(il);
                        near_lig[il].push_back(ir);
                    }
                }
            }
            // cout<<"near rec"<<endl;
            /////sasa receptor
            double sasa_rec=0;
            for(int ir=0;ir!=num_receptor_atom;ir++)
            {
                double rrr_rw=rec.noH[ir].rrr+P.r_water;
                double dsita=P.step_SASA/rrr_rw;
                for(int j=0;j!=rec.radius_rw[ir].size();j++)
                {
                    double sita=(j+1)*dsita;
                    double dphi=dsita*sin(sita);
                    
                    for(int k=0;k!=rec.radius_rw[ir][j].size();k++)
                    {
                        
                        bool calSASA = false;
                        for(int inear=0;inear!=near_rec[ir].size();inear++)
                        {
                            const COORDINATE& ll=irank.position[near_rec[ir][inear]];
                            const ATOM& lla=lig.noH[near_rec[ir][inear]];
                            double dis2lig=(rec.radius_rw[ir][j][k].x-ll.x)*(rec.radius_rw[ir][j][k].x-ll.x)
                                          +(rec.radius_rw[ir][j][k].y-ll.y)*(rec.radius_rw[ir][j][k].y-ll.y)
                                          +(rec.radius_rw[ir][j][k].z-ll.z)*(rec.radius_rw[ir][j][k].z-ll.z);
                            if(dis2lig<=((lla.rrr + P.r_water)*(lla.rrr + P.r_water)))
                            {
                                calSASA=true;break;
                            }
                        }
                        if(calSASA)sasa_rec+=rrr_rw*rrr_rw*sin(sita)*dsita*dphi;
                    }
                }
            }
            // cout<<"sasa rec"<<endl;
            //////sasa ligand
            double sasa_lig=0;
            for(int il=0;il!=num_ligand_atom;il++)
            {
                double rrr_rw=lig.noH[il].rrr+P.r_water;
                double dsita=P.step_SASA/rrr_rw;
                int jmax=P.pi/dsita;
                for(int j=1;j<=jmax;j++)
                {
                    double sita=j*dsita;
                    double dphi=dsita*sin(sita);
                    int kmax=2.*P.pi/dphi+1;
                    for(int k=1;k<=kmax;k++)
                    {
                        double phi=k*dphi;
                        double rx=rrr_rw*sin(sita)*cos(phi)+irank.position[il].x;
                        double ry=rrr_rw*sin(sita)*sin(phi)+irank.position[il].y;
                        double rz=rrr_rw*cos(sita)+irank.position[il].z;
                        bool calSASA = false;
                        for(int inear=0;inear!=near_lig[il].size();inear++)
                        {
                            const ATOM& rra=rec.noH[near_lig[il][inear]];
                            double dis2lig=(rx-rra.x)*(rx-rra.x)+(ry-rra.y)*(ry-rra.y)+(rz-rra.z)*(rz-rra.z);
                            if(dis2lig<=((rra.rrr + P.r_water)*(rra.rrr + P.r_water)))
                            {
                                calSASA=true;break;
                            }
                        }
                        if(calSASA)sasa_lig+=rrr_rw*rrr_rw*sin(sita)*dsita*dphi;
                    }
                }
            }
            irank.sasa = P.gamma_SASA*(sasa_rec+sasa_lig)/P.rate_kcal_to_kt;
            irank.calculate_eng_rerank(P.weight_focus2);
        }
    }

}

int main(int argc, char** argv ) 
{
    int opt = 0 ;
    string receptor_file = "";
    string ligand_file = "";
    string reference_file = "";
    string output_prefix = "";
    int num_threads = 4;


    while((opt = getopt(argc, argv, "hi:l:r:o:n:")) != -1) {
        switch(opt)
        {
            case 'h':
                printf("Usage ./RLDOCK -i <receptor.mol2> -l <ligand.mol2> -o <output prefix> -n <thread number> -r <native_pose.mol2>");
                return 0;
            case 'i':
                receptor_file = string(optarg);
                break;
            case 'l':
                ligand_file = string(optarg);
                break;
            case 'o':
                output_prefix = string(optarg);
                break;
            case 'n':
                num_threads = atoi(optarg);
                break;
            case 'r':
                reference_file = string(optarg);
                break;
        }
    }
    if(receptor_file.length() == 0){cout<<"Receptor mol2 file needed."<<endl;exit(3);}
    else if(ligand_file.length() == 0){cout<<"Ligand mol2 file needed."<<endl;exit(3);}
    else {
        cout<<"Receptor file: "<<receptor_file<<endl;
        cout<<"Ligand   file: "<<ligand_file<<endl;
        if(reference_file.length() == 0){
            cout<<"No reference ligand"<<endl;}
        else{cout<<"Reference ligand file:"<<reference_file<<endl;}
        if(output_prefix.length() == 0){
            output_prefix = "output";
            cout<<"Job Name is set as default: "<<output_prefix<<endl;}
        else{cout<<"Job Name: "<<output_prefix<<endl;}
    }
    std::cout<<"*********start*******"<<endl;
    
    PARAMETER P;
    P.SetSphere("src/sphere.dat");
    P.SetAngle();
    P.InputUpdate(receptor_file,ligand_file,reference_file,output_prefix,num_threads);
    
    FORCEFIELD F;
    F.Initialization(P);

    time_t t_start = time(NULL);
    COMPLEX cmpl;
    READFILE(P, cmpl.rec, P.receptor_file);
    READFILE(P, cmpl.conformers, P.ligand_file);
    if(P.reference_flag){READFILE(P, cmpl.lig, P.reference_file);}
    else{cmpl.lig.push_back(cmpl.conformers[0]);}
    cmpl.NumInput();

    if(cmpl.lig_atom_num>80){P.UpdatAngle(20);}
    else if(cmpl.lig_atom_num>50){P.UpdatAngle(15);}

    cmpl.Initialization(P,F);
    std::cout<<"Input and Initialization Finish "<<endl;
    cout<<cmpl.conformers_num<<endl;
    ////////////////////////////////////////////////
    ///////////step1 pocket generation//////////////
    ////////////////////////////////////////////////
    string pocket_filename = P.output_prefix+"_pocket.dat";
    ifstream pocket_file(pocket_filename.c_str());
    string usepose_filename = P.output_prefix+"_usepose.dat";
    ifstream usepose_file(usepose_filename.c_str());
    vector<POCKET> pockets;
    vector<THREAD_PARAMETER> thread_p;
    if( !pocket_file || !usepose_file)
    {
        GRID grid(P.step_pocket, 2*P.r_water);
        grid.Set_grid_pocket(cmpl.rec[0]);
        grid.Set_grid_int();

        for(auto&arec:cmpl.rec[0].noH)
        {
            double disaway = arec.rrr + P.radius_small_sphere;
            grid.Find_RNA_grid(P, disaway, arec);
        }
        std::cout<<"Grid for binding site Finish"<<endl;

        thread t[P.num_threads];
        Thread_Distribution(thread_p, grid.Nnx, P.num_threads);

        vector<vector<POCKET> > pockets_sub(P.num_threads);

        for(unsigned int i = 0; i != P.num_threads; i++)
        {
            t[i] = thread(Find_Pocket, ref(P), ref(grid), ref(cmpl.rec[0].noH), ref(thread_p[i]), ref(pockets_sub[i]));
        }
        for(unsigned int i = 0; i != P.num_threads; ++i)
        {
            t[i].join();
            pockets.insert(pockets.end(),pockets_sub[i].begin(),pockets_sub[i].end());
        }
        std::cout<<"Find "<<pockets.size()<<" pockets"<<endl;

        GRID grid_energy(P.step_grid_energy, P.LJ_cut_new);
        grid_energy.Set_grid_common(cmpl.rec[0]);
        grid_energy.Initial_grid_lj(P);
        for(auto & arec:cmpl.rec[0].noH)
        {
            grid_energy.Set_LJ(P, P.LJ_cut_new, arec, F);
        }

        vector<COORDINATE>centers;
        for(auto&a:cmpl.conformers)
        {
            centers.insert(centers.end(),a.center.begin(),a.center.end());
        }

        int num_thread_tmp,center_size = centers.size();
        num_thread_tmp=min(P.num_threads,center_size);
        thread t1[num_thread_tmp];
        Thread_Distribution(thread_p, center_size, num_thread_tmp);

        LABEL void_label;
        vector< vector <LABEL> > tmp_lj(cmpl.conformers_num, vector<LABEL>(cmpl.lig_atom_num,void_label));
        for(int i=0; i != cmpl.conformers_num; i++ )
        {
            for(int j=0; j != cmpl.lig_atom_num; j++)
            {
                tmp_lj[i][j].initial_idx=i;
                tmp_lj[i][j].lig_idx=j;
                tmp_lj[i][j].lj=1000.0;

            }
        }
        for(auto& apoc:pockets){apoc.LJ_min = tmp_lj;}
        pockets_sub.resize(num_thread_tmp);
        for(auto& ap:pockets_sub){ap=pockets;}
        vector< vector< vector< int > > >chose_pose_idx_sub(num_thread_tmp, vector< vector<int> >(center_size));

        for(unsigned int i = 0; i != num_thread_tmp; ++i)
        {
            t1[i] = thread(Select_Pose,ref(P),ref(thread_p[i]),ref(centers),ref(cmpl),ref(grid_energy),ref(pockets_sub[i]),ref(chose_pose_idx_sub[i]));
        }
        for(unsigned int i = 0; i != num_thread_tmp; ++i)
        {
            t1[i].join();
            
            for(int isite=0;isite!=pockets.size();isite++)
            {
                for(int ic = thread_p[i].ibegin; ic <= thread_p[i].iend;ic++)
                {
                    int i_conf = ic/cmpl.lig_atom_num;
                    int i_atom = ic % cmpl.lig_atom_num; 
                    pockets[isite].LJ_min[i_conf][i_atom] =  pockets_sub[i][isite].LJ_min[i_conf][i_atom];
                }
                if(pockets_sub[i][isite].lj_site_min < pockets[isite].lj_site_min)
                pockets[isite].lj_site_min = pockets_sub[i][isite].lj_site_min;    
            }
           
            for(int ic = thread_p[i].ibegin; ic <= thread_p[i].iend;ic++)
            {
                int i_conf = ic/cmpl.lig_atom_num;
                int i_atom = ic % cmpl.lig_atom_num;
                cmpl.conformers[i_conf].usepose[i_atom]=chose_pose_idx_sub[i][ic];
            }
        }
        ofstream out_usepose_file(usepose_filename.c_str());
        for(int iconf=0;iconf!=cmpl.conformers_num ;iconf++)
        {
            for(int iatom=0;iatom!=cmpl.lig_atom_num ;iatom++)
            {
                out_usepose_file<<"<pos> "<<iconf<<" "<<iatom<<" "<<cmpl.conformers[iconf].usepose[iatom].size()<<endl;
                for(auto&a:cmpl.conformers[iconf].usepose[iatom])out_usepose_file<<a<<endl;
            }
        }
        out_usepose_file.close();

        double cost_t = time(NULL) - t_start;
        ofstream out_pocket_file(pocket_filename.c_str());
        out_pocket_file<<argv[1]<<" rec_num: "<<cmpl.rec_atom_num<<" lig_num: "<<cmpl.lig_atom_num<<" pocket_num: "<<pockets.size()<<" conformer_num: "<<cmpl.conformers_num<<" "<<cost_t<<" s"<<endl;
        out_pocket_file<<"site_index    x           y           z    distance_to_reference_ligand   lj_eng"<<endl;
        sort(pockets.begin(),pockets.end(),less_pocket_lj);
        for(auto& a:pockets)
        {
            for(auto& p:a.LJ_min)
            {
                sort(p.begin(),p.end(),less_label_lj);
            }
        }
        int num=0;
        for(auto& a:pockets)
        {
            num++;
            double dmin=-1;
            int label_pocket_to_native=-1;
            if(P.reference_flag)
            {
                for(auto&l:cmpl.lig[0].noH)
                {
                    double dis=(l.x-a.x)*(l.x-a.x);
                    dis += (l.y-a.y)*(l.y-a.y);
                    dis += (l.z-a.z)*(l.z-a.z);
                    dis = sqrt(dis);
                    if(dis<dmin){dmin=dis;label_pocket_to_native=l.index;}
                }
            }
            out_pocket_file.setf(ios::fixed);
            out_pocket_file<<"site"<<num<<right<<setw(12)<<setprecision(3)<<a.x<<setw(12)<<a.y<<setw(12)<<a.z<<setw(8)<<label_pocket_to_native<<setw(15)<<dmin<<setw(15)<<a.lj_site_min<<endl;
            for(int i = 0; i!= cmpl.conformers_num; i++)
            {
               out_pocket_file<<setw(4)<<left<<i;
               for(int j = 0; j!= cmpl.lig_atom_num; j++)
               {
                   out_pocket_file<<setw(4)<<a.LJ_min[i][j].lig_idx <<setw(15)<<a.LJ_min[i][j].lj ;
               }
               out_pocket_file<<endl;
            }
        }
        out_pocket_file.close();
        cout<<"_pocket.dat _usepose.dat Finish"<<endl;
        for(int ipck=0;ipck!=pockets.size();ipck++)
        {
            
            vector<int>low_idx(P.TOP_step2_conformer);
            vector<RANK_PARAMETER> rankuse(cmpl.conformers_num);
            for(int ir=0;ir!=cmpl.conformers_num;ir++)
            { 
                rankuse[ir].index = ir;
                rankuse[ir].val = pockets[ipck].LJ_min[ir][0].lj;
            }
            sort(rankuse.begin(),rankuse.end(),rank_less);
            for(int ii=0;ii!=P.TOP_step2_conformer;ii++)
            {
                low_idx[ii]=rankuse[ii].index;
            }
            /////////////
            pockets[ipck].chosen_atom_lig.resize(cmpl.conformers_num);
            for(int top_iconf=0;top_iconf!=low_idx.size();top_iconf++)
            {
                pockets[ipck].chosen_atom_lig[low_idx[top_iconf]].resize(P.TOP_step2_pose);
                for(int i2p=0;i2p!=P.TOP_step2_pose;i2p++)
                {
                    pockets[ipck].chosen_atom_lig[low_idx[top_iconf]][i2p]=pockets[ipck].LJ_min[low_idx[top_iconf]][i2p].lig_idx;
                }
            }
        }
        

        if(pockets.size()>P.TOP_step2_site)pockets.resize(P.TOP_step2_site);
    }
    else
    {
        int num_pockets = 0;
        string sline;
        while(getline(pocket_file,sline))
        {
            if(sline.find("site")!=string::npos)
            {
                POCKET pp;
                std::istringstream ss(sline);
                std::string buf;
                std::vector<std::string> token;
                while(ss >> buf) token.push_back(buf);
                pp.index = stod(sline.substr(4,4));
                pp.x = stod(token[1]);
                pp.y = stod(token[2]);
                pp.z = stod(token[3]);
                pp.dis2rec = stod(token[5]);
               

                vector< vector<int> >pose_tmp(cmpl.conformers_num,vector<int>(P.TOP_step2_pose,0.0));
                vector< double > pose_val(cmpl.conformers_num,0.0);
                vector< int > pose_conformer(cmpl.conformers_num,1);
                for(int i=0; i!=cmpl.conformers_num;i++)
                {
                    getline(pocket_file,sline);
                    std::istringstream ss(sline);
                    std::string buf;
                    std::vector<std::string> token;
                    while(ss >> buf) token.push_back(buf);
                    for(int j=0;j!=P.TOP_step2_pose;j++)
                    {
                        pose_tmp[i][j] = stoi(token[2*j+1]);
                    }
                    pose_val[i] = stod(token[2]);
                }
                vector< double > pose_rank = pose_val;
                sort(pose_rank.begin(),pose_rank.end());
                double biaozhun = pose_rank[P.TOP_step2_conformer];
                pp.chosen_atom_lig.resize(cmpl.conformers_num);
                for(int i=0;i!=cmpl.conformers_num;i++)
                {
                    if(pose_val[i]<biaozhun)
                    {
                        pp.chosen_atom_lig[i]=pose_tmp[i];
                    }
                }
                pockets.push_back(pp);
                num_pockets++;
                if(num_pockets==P.TOP_step2_site)break;
            }
        }
        pocket_file.close();

  
        while(getline(usepose_file,sline))
        {
            if(sline.find("<pos>")!=string::npos)
            {
                std::istringstream ss(sline);
                std::string buf;
                std::vector<std::string> token;
                while(ss >> buf) token.push_back(buf);

                int iconf = stoi(token[1]);
                int iatom = stoi(token[2]);
                int posesize = stoi(token[3]);
                cmpl.conformers[iconf].usepose[iatom].resize(posesize);
                for(int ii=0;ii!=posesize;ii++)
                {
                    getline(usepose_file,sline);
                    cmpl.conformers[iconf].usepose[iatom][ii]=stoi(sline);
                }       
            }
        }
        usepose_file.close();
        std::cout<<pockets.size()<<" pockets input finish"<<endl;
    }

    BornRadius(P,pockets,cmpl.rec[0]);
    ////////////////////////////////////////////////
    ///////////step2 _SF_low.dat       //////////////
    ////////////////////////////////////////////////
    time_t f1_t=time(NULL);
    string sfl_filename = P.output_prefix+"_SF_low.dat";
    ifstream sfl_file(sfl_filename.c_str());
    vector<POSE> per_top_each_site;
    per_top_each_site.resize(P.TOP_focus1_num_each_site+1);
    for(int it=0;it!=per_top_each_site.size();it++)per_top_each_site[it].eng=100000000+it;
    vector<vector<POSE> > top_each_site(pockets.size(),per_top_each_site);
    if(!sfl_file)
    {
        thread t[P.num_threads];
         vector<vector<vector<POSE> > >top_each_site_sub(P.num_threads,top_each_site);
        Thread_Distribution(thread_p,pockets.size(),P.num_threads);

        for(unsigned int i = 0; i != P.num_threads; ++i)
        {             
            t[i] = std::thread(SF_low,std::ref(P), std::ref(thread_p[i]), std::ref(cmpl),std::ref(pockets),std::ref(F), std::ref(top_each_site_sub[i]));
        }    
        for(unsigned int i = 0; i != P.num_threads; ++i)
        {
            t[i].join();
            for(int idx = thread_p[i].ibegin;idx<=thread_p[i].iend;idx++)
            {
                top_each_site[idx]=top_each_site_sub[i][idx];
            }
        } 

        ofstream out_sfl_file(sfl_filename.c_str());
        out_sfl_file<<std::fixed<<setprecision(3);
        double cost_t = time(NULL) - f1_t;
        size_t found = P.output_prefix.find_last_of("/\\");
        out_sfl_file<<"<JobName> "<<P.output_prefix.substr(found+1)<<" receptor_size "<<cmpl.rec_atom_num<<" lig_size "<<cmpl.lig_atom_num<<" "<<cost_t<<"s"<<endl;
        for(int isite=0;isite!=top_each_site.size();isite++)
        {
            for(int irank=0;irank!=P.TOP_focus1_num_each_site;irank++)
            {
                // cout<<isite<<" "<<irank<<" "<<top_each_site[isite][irank].eng<<endl;
                if(top_each_site[isite][irank].eng<100000000)
                {
                    out_sfl_file<<"<site> "<<isite<<" "<<irank<<" /10 iconf "<<top_each_site[isite][irank].iconf<<" iatom "<<top_each_site[isite][irank].iatom<<" ipose "<<top_each_site[isite][irank].ipose<<" rmsd "<<top_each_site[isite][irank].rmsd<<endl;
                    out_sfl_file<<"<energy> "<<top_each_site[isite][irank].eng<<" lj "<<top_each_site[isite][irank].lj<<" ele "<<top_each_site[isite][irank].ele<<" pol "<<top_each_site[isite][irank].pol<<" self "<<top_each_site[isite][irank].self<<" sasa "<<top_each_site[isite][irank].sasa<<" hb "<<top_each_site[isite][irank].hb<<" internal_lj "<<top_each_site[isite][irank].internal_lj<<endl;
                    out_sfl_file<<"<pocket position> "<<pockets[isite].x<<" "<<pockets[isite].y<<" "<<pockets[isite].z<<endl;
                    for(int ilig=0;ilig!=cmpl.lig_atom_num;ilig++)
                    {
                        out_sfl_file<<right<<setw(4)<<ilig<<setw(5)<<cmpl.conformers[0].noH[ilig].name<<setw(3)<<cmpl.conformers[0].noH[ilig].element<<setw(8)<<top_each_site[isite][irank].position[ilig].x<<setw(8)<<top_each_site[isite][irank].position[ilig].y<<setw(8)<<top_each_site[isite][irank].position[ilig].z<<endl;
                    }
                }
            }
        }
        out_sfl_file.close();
         cout<<"_SF_low.dat Finish"<<endl;

    }
    else
    {
        int num_ligand_atom;
        string sline;
        while(getline(sfl_file,sline))
        {
            if(sline.find("<PDBinfo>")!=string::npos)
            {
                std::istringstream ss(sline);
                std::string buf;
                std::vector<std::string> token;
                while(ss >> buf) token.push_back(buf);

                num_ligand_atom = stoi(token[5]);
                if(num_ligand_atom!=cmpl.lig_atom_num){cout<<"focus1 "<<num_ligand_atom<<" != lig.noH "<<cmpl.lig_atom_num<<endl;exit(3);}
            }
            if(sline.find("<site>")!=std::string::npos)
            {
                std::istringstream ss(sline);
                std::string buf;
                std::vector<std::string> token;
                while(ss >> buf) token.push_back(buf);
                int  isite = stoi(token[1]),irank = stoi(token[2]);
                top_each_site[isite][irank].isite=stoi(token[1]);
                top_each_site[isite][irank].iconf=stoi(token[5]);
                top_each_site[isite][irank].iatom=stoi(token[7]);
                top_each_site[isite][irank].ipose=stod(token[9]);
                top_each_site[isite][irank].rmsd=stod(token[11]);

                getline(sfl_file,sline);
                std::istringstream ss2(sline);
                token.clear();
                while(ss2 >> buf) token.push_back(buf);
                top_each_site[isite][irank].eng=stod(token[1]);
                top_each_site[isite][irank].lj=stod(token[3]);
                top_each_site[isite][irank].ele=stod(token[5]);
                top_each_site[isite][irank].pol=stod(token[7]);
                top_each_site[isite][irank].self=stod(token[9]);
                top_each_site[isite][irank].sasa=stod(token[11]);
                top_each_site[isite][irank].hb=stod(token[13]);
                int iconf =  top_each_site[isite][irank].iconf;
                top_each_site[isite][irank].internal_lj=cmpl.conformers[iconf].internal_lj;

                for(int ilig=0;ilig!=num_ligand_atom;ilig++)
                {
                    getline(sfl_file,sline);
                    std::istringstream ss3(sline);
                    token.clear();
                    while(ss3 >> buf) token.push_back(buf);
                    top_each_site[isite][irank].position[ilig].x=stod(token[3]);
                    top_each_site[isite][irank].position[ilig].y=stod(token[4]);
                    top_each_site[isite][irank].position[ilig].z=stod(token[5]);
                }
            }
        }
        sfl_file.close();

    }

    ////////////////////////////////////////////////
    ///////////step3 _SF_high.dat       //////////////
    ////////////////////////////////////////////////
    time_t f2_t=time(NULL);
    thread t_f2[P.num_threads];
    Thread_Distribution(thread_p,top_each_site.size(),P.num_threads);
    
    for(unsigned int i = 0; i != P.num_threads; ++i)
    {             
        t_f2[i] = std::thread(SF_high,std::ref(P), std::ref(thread_p[i]), std::ref(cmpl),std::ref(pockets),std::ref(F), std::ref(top_each_site));
    }    
    for(unsigned int i = 0; i != P.num_threads; ++i)
    {
        t_f2[i].join();
    }
    vector<POSE> all_pose;
    for(auto&a:top_each_site)
    {
        a.resize(P.TOP_focus1_num_each_site);
        for(auto&p:a)
        {
            all_pose.push_back(p);
        }
    }
    sort(all_pose.begin(),all_pose.end(),less_eng_rerank);
    
    string sfh_filename = P.output_prefix+"_SF_high.dat";
    ofstream out_sfhfilename(sfh_filename.c_str());
    out_sfhfilename<<std::fixed;
    double cost_t = time(NULL) - f2_t;
    size_t found = P.output_prefix.find_last_of("/");
    out_sfhfilename<<"<JobName> "<<P.output_prefix.substr(found+1)<<" receptor_size "<<cmpl.rec_atom_num<<" lig_size "<<cmpl.lig_atom_num<<" "<<cost_t<<"s"<<endl;
    
        for(int irank=0;irank!=all_pose.size();irank++)
        {
            if(all_pose[irank].eng<100000000)
            {
                int isite =all_pose[irank].isite;
                // outputfile<<right<<setw(4)<<isite<<setw(3)<<irank<<setw(3)<<top_each_site[isite][irank].iconf<<setw(3)<<top_each_site[isite][irank].iatom<<setw(6)<<top_each_site[isite][irank].ipose<<setw(6)<<setprecision(3)<<top_each_site[isite][irank].rmsd
                // <<setw(8)<<top_each_site[isite][irank].lj<<setw(8)<<top_each_site[isite][irank].ele<<setw(8)<<top_each_site[isite][irank].pol<<setw(8)<<top_each_site[isite][irank].self_lig<<setw(8)<<top_each_site[isite][irank].sasa<<setw(8)<<top_each_site[isite][irank].hb<<setw(8)<<top_each_site[isite][irank].self_rec
                out_sfhfilename<<"<site> "<<all_pose[irank].isite<<" rank "<<irank<<" iconf "<<all_pose[irank].iconf<<" iatom "<<all_pose[irank].iatom<<" ipose "<<all_pose[irank].ipose<<" rmsd "<<all_pose[irank].rmsd<<endl;
                out_sfhfilename<<"<energy> "<<all_pose[irank].eng_rerank<<" lj "<<all_pose[irank].lj<<" ele "<<all_pose[irank].ele<<" pol "<<all_pose[irank].pol<<" self_lig "<<all_pose[irank].self_lig<<" sasa "<<all_pose[irank].sasa<<" hb "<<all_pose[irank].hb<<" self_rec "<<all_pose[irank].self_rec<<" internal_lj "<<all_pose[irank].internal_lj<<endl;
                out_sfhfilename<<"<pocket position> "<<pockets[isite].x<<" "<<pockets[isite].y<<" "<<pockets[isite].z<<endl;
                for(int ilig=0;ilig!=cmpl.lig_atom_num;ilig++)
                {
                    out_sfhfilename<<right<<setw(4)<<ilig<<setw(5)<<cmpl.conformers[0].noH[ilig].name<<setw(3)<<cmpl.conformers[0].noH[ilig].element<<setw(8)<<setprecision(3)<<all_pose[irank].position[ilig].x<<setw(8)<<all_pose[irank].position[ilig].y<<setw(8)<<all_pose[irank].position[ilig].z<<endl;
                }
            }
        }
    
    out_sfhfilename.close();
    cout<<"_SF_high.dat Finish"<<endl;
}