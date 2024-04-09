#include "RLDOCK_CASE.h"


void CASE::Move_to_Center_Easy()
{
    this->center = this->noH;
    xcen=0,ycen=0,zcen=0;
    for(auto&a:this->noH){xcen+=a.x;ycen+=a.y;zcen+=a.z;}
    xcen/=1.0*this->noH.size();
    ycen/=1.0*this->noH.size();
    zcen/=1.0*this->noH.size();
    for(auto&a:this->center){a.x-=xcen;a.y-=ycen;a.z-=zcen;}

}

void CASE::Find_Max_and_Min()
{
    xmax=-1000000.0;
    for(auto&a:this->noH){if(xmax<a.x)xmax=a.x;}

    ymax=-1000000.0;
    for(auto&a:this->noH){if(ymax<a.y)ymax=a.y;}

    zmax=-1000000.0;
    for(auto&a:this->noH){if(zmax<a.z)zmax=a.z;}

    xmin=1000000.0;
    for(auto&a:this->noH){if(xmin>a.x)xmin=a.x;}
    
    ymin=1000000.0;
    for(auto&a:this->noH){if(ymin>a.y)ymin=a.y;}

    zmin=1000000.0;
    for(auto&a:this->noH){if(zmin>a.z)zmin=a.z;}

}

void CASE::cal_distance()
{
    this->distance.resize(this->noH.size(),vector<double>(this->noH.size(),0.0));
    this->distance_square.resize(this->noH.size(),vector<double>(this->noH.size(),0.0));
    for(int ilig=0;ilig!=this->noH.size();ilig++)
    {
        for(int jlig=0;jlig!=this->noH.size();jlig++)
        {
            if(ilig!=jlig)
            {
                this->distance[ilig][jlig]=this->noH[ilig].dis(this->noH[jlig]);
                this->distance_square[ilig][jlig] = this->distance[ilig][jlig]*this->distance[ilig][jlig];
            }
        }
    }
}

void CASE::get_rb0_and_cal_pol(const PARAMETER& P)
{
    for(int ir = 0; ir !=this->noH.size(); ir++)
    {
        double intb=0., radi = this->noH[ir].rrr;
        for(int jr = 0; jr !=this->noH.size(); jr++)
        {
            if(jr!=ir)
            {
                double dis=(this->noH[ir].x - this->noH[jr].x)*(this->noH[ir].x - this->noH[jr].x);
                dis += (this->noH[ir].y - this->noH[jr].y)*(this->noH[ir].y - this->noH[jr].y);
                dis += (this->noH[ir].z - this->noH[jr].z)*(this->noH[ir].z - this->noH[jr].z);
                dis = sqrt(dis);
                double radj = noH[jr].rrr;
                double radj_sb = noH[jr].sb*noH[jr].rrr;
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
        this->noH[ir].rb0 = 1./(1./radi-intb);
    }
    this->pol=0;
    for(int il = 0; il != this->noH.size(); il++)
    {
        for(int jl = 0; jl != il; jl++)
        {
            double rr_square = this->noH[il].rb0*this->noH[jl].rb0;
            this->pol+=P.pol_par * this->noH[il].q * this->noH[jl].q/sqrt(this->distance_square[il][jl]+rr_square*exp(-this->distance_square[il][jl]/4./rr_square)); 
        }
    }
    
}

void CASE::cal_internal_lj(const FORCEFIELD& F, const PARAMETER& P)
{
    this->internal_lj=0;
    for(int ilig =0; ilig!=this->noH.size();ilig++)
    for(int jlig =0; jlig!=ilig;jlig++)
    {
        vector<int>::iterator it = find(this->noH[ilig].bonded_atom.begin(),this->noH[ilig].bonded_atom.end(),jlig);
        if(it==this->noH[ilig].bonded_atom.end() && ilig!=jlig)
        {
        int dis_idx = distance[ilig][jlig]/F.eql_lj[this->noH[ilig].type_idx][this->noH[jlig].type_idx]*P.LJstep;
        if(dis_idx<P.LJstep)this->internal_lj += F.LJ_potential[this->noH[ilig].type_idx][this->noH[jlig].type_idx][dis_idx];
        }
    }
}

void CASE::set_rw_position(const PARAMETER & P)
{
    this->radius_rw.resize(this->noH.size());
    for(int ir=0;ir!=this->noH.size();ir++)
    {
        double rrr_rw=this->noH[ir].rrr+P.r_water;
        double dsita=P.step_SASA/rrr_rw;
        int jmax=P.pi/dsita;
        this->radius_rw[ir].resize(jmax);
        for(int j=1;j<=jmax;j++)
        {
            double sita=j*dsita;
            double dphi=dsita*sin(sita);
            int kmax=2.*P.pi/dphi+1;
            this->radius_rw[ir][j-1].resize(kmax);
            for(int k=1;k<=kmax;k++)
            {
                double phi=k*dphi;
                this->radius_rw[ir][j-1][k-1].x=rrr_rw*sin(sita)*cos(phi)+this->noH[ir].x;
                this->radius_rw[ir][j-1][k-1].y=rrr_rw*sin(sita)*sin(phi)+this->noH[ir].y;
                this->radius_rw[ir][j-1][k-1].z=rrr_rw*cos(sita)+this->noH[ir].z;
            }
        }
    }
}

void COMPLEX::Mol2_Info_Get(string lig_file_name)
{
    for(int i=1;i!=this->conformers.size();i++)
    {
        if(this->conformers[i].noH.size()!=this->conformers[0].noH.size())
        {
            cout<<"The number of heavy atoms of the "<<i+1<<"th conformer is different from that of first conformer."<<endl;
            cout<<"The "<<i+1<<"th conformer has "<<this->conformers[i].noH.size()<<" heavy atoms."<<endl;
            cout<<"The first conformer has "<<this->conformers[0].noH.size()<<" heavy atoms."<<endl;
            exit(3);
        }
    }
    this->lig_atom_num = this->conformers[0].noH.size();
    this->rec_atom_num = this->rec[0].noH.size();
    this->conformers_num = this->conformers.size();

    ifstream infile(lig_file_name.c_str());
    string sline;
    while (getline(infile,sline))
    {
        if(sline.find("@<TRIPOS>MOLECULE")!=std::string::npos)
        {
            if(this->lig_mol2.atom_num != 0)break;
            CASE tmp_mol;
            int num_atom,num_bond,num_sub,num_feat,num_sets;

            getline(infile,sline);
            this->lig_mol2.mol_name = sline;

            getline(infile,sline);
            std::istringstream ss(sline);
            std::string buf;
            std::vector<std::string> token;
            while(ss >> buf) token.push_back(buf);
            if(token.size() != 5) {
                std::cout << "unable to read mol2 file 3rd line! this line should contain 5 integers, please check the file format!" << std::endl;
                exit(2);
            }
            num_atom = stoi(token[0]);
            num_bond = stoi(token[1]);
            num_sub = stoi(token[2]);
            num_feat = stoi(token[3]);
            num_sets = stoi(token[4]);

            getline(infile,sline);
            getline(infile,sline);
            // getline(infile,sline);
            this->lig_mol2.charge_type=sline;
            std::map<int,int> idx_update;

            while(getline(infile,sline))
            {
                if(sline.find("@<TRIPOS>ATOM")!=std::string::npos)
                {
                    for (int ia=0,iheavy=1;ia!=num_atom;ia++)
                    {
                        getline(infile,sline);
                        ATOM a;
                        a.mol2_input(sline);
                        if(a.element != "H" && a.element != "h")
                        {
                            idx_update[ia+1]=iheavy;
                            iheavy++;
                            this->lig_mol2.atoms.push_back(a);
                        }
                        else
                        {
                            idx_update[ia+1]=0;
                        }
                    }
                    this->lig_mol2.atom_num = this->lig_mol2.atoms.size();
                }
                if(sline.find("@<TRIPOS>BOND")!=std::string::npos)
                {
                    for(int ib=0,ib_real=1;ib!=num_bond;ib++)
                    {
                        getline(infile,sline);
                        BOND b;
                        b.mol2_input(sline);
                        int new_left = idx_update[b.left_id];
                        int new_right = idx_update[b.right_id];
                        if(new_left!=0 && new_right != 0)
                        {
                            b.index=ib_real;
                            b.update(ib_real,new_left,new_right);
                            this->lig_mol2.bonds.push_back(b);
                            ib_real++;
                        }
                    }
                    this->lig_mol2.bond_num = this->lig_mol2.bonds.size();
                    break;
                    // if(sline.find("@<TRIPOS>SUBSTRUCTURE")!=std::string::npos)
                    // {
                    //     for(int is=0;is!=num_sub;is++)
                    //     {
                    //         getline(infile,sline);
                    //         SUBSTRUCTURE s;
                    //         s.mol2_input(sline);
                    //         s.root_atom = idx_update[s.root_atom];
                    //         this->lig_mol2.subs.push_back(s);
                    //     }   
                    //     this->lig_mol2.sub_num = this->lig_mol2.subs.size();
                    //     break;             
                    // }
                }
                
            } 
            
        }
    }
    infile.close();
}

void COMPLEX::Initialization(const PARAMETER& P, const FORCEFIELD& F)
{
    for(auto&arec:this->rec)
    {
        arec.Find_Max_and_Min();
        arec.cal_distance();
        arec.get_rb0_and_cal_pol(P);
        arec.set_rw_position(P);
    }

    for(auto&alig:this->conformers)
    {
        alig.cal_distance();
        alig.get_rb0_and_cal_pol(P);
        alig.Move_to_Center_Easy();
        alig.cal_internal_lj(F,P);
        alig.usepose.resize(this->lig_atom_num);
    }

}

void READFILE(const PARAMETER P, vector<CASE>& in_mol, string filename)
{
    ifstream infile(filename.c_str());
    string sline;
    int ii=0;
    while (getline(infile,sline))
    {
        if(sline.find("@<TRIPOS>MOLECULE")!=std::string::npos)
        {
            
            CASE tmp_mol;
            int num_atom,num_bond;
            getline(infile,sline);
            
            infile>>sline;num_atom=stoi(sline);
            infile>>sline;num_bond=stoi(sline);

            while(getline(infile,sline))
            {
                if(sline.find("@<TRIPOS>ATOM")!=std::string::npos)
                {
                    int index_noH=0, index_HH=0;
                    for (int ia=0;ia!=num_atom;ia++)
                    {
                        getline(infile,sline);
                        std::istringstream ss(sline);
                        std::string buf;
                        std::vector<std::string> token;
                        while(ss >> buf) token.push_back(buf);

                        ATOM a;
                        a.index_original = stoi(token[0]);
                        a.name = token[1];
                        a.x = stod(token[2]);
                        a.y = stod(token[3]);
                        a.z = stod(token[4]);
                        a.type = token[5];
                        a.resi_index = stoi(token[6]);
                        a.resi_name = token[7];
                        a.q = stod(token[8]);

                        if(a.type.find(".")!=std::string::npos)
                        {
                            string::size_type position;
                            position = a.type.find(".");
                            a.element = a.type.substr(0,position);
                        }
                        else
                        {
                            a.element = a.type;
                        }

                        if(a.element == "H" || a.element == "h")
                        {
                            a.rrr = P.radius_H;
                            a.sb = P.born_scale_H;
                            a.type_idx = 0;
                        }
                        else if(a.element == "C" || a.element == "c")
                        {
                            a.rrr = P.radius_C;
                            a.sb = P.born_scale_C;
                            a.type_idx = 1;

                        }
                        else if(a.element == "N" || a.element == "n")
                        {
                            a.rrr = P.radius_N;
                            a.sb = P.born_scale_N;
                            a.type_idx = 2;

                        }
                        else if(a.element == "O" || a.element == "o")
                        {
                            a.rrr = P.radius_O;
                            a.sb = P.born_scale_O;
                            a.type_idx = 3;
                        }
                        else if(a.element == "P" || a.element == "p")
                        {
                            a.rrr = P.radius_P;
                            a.sb = P.born_scale_P;
                            a.type_idx = 4;
                        }
                        else if(a.element == "S" || a.element == "s")
                        {
                            a.rrr = P.radius_S;
                            a.sb = P.born_scale_S;
                            a.type_idx = 5;
                        }
                        else
                        {
                            a.rrr = P.radius_C;
                            a.sb = P.born_scale_C;
                            a.type_idx = 6;
                        }

                        if(a.element != "H" && a.element != "h")
                        {
                            a.index = index_noH;
                            index_noH++;
                            tmp_mol.noH.push_back(a);
                            tmp_mol.wH.push_back(a);
                        }
                        else
                        {
                            a.index = index_HH;
                            index_HH++;
                            tmp_mol.HH.push_back(a);
                            tmp_mol.wH.push_back(a);
                        }

                    }

                }
                if(sline.find("@<TRIPOS>BOND")!=std::string::npos)
                {
                    for(int ib=0;ib!=num_bond;ib++)
                    {
                        getline(infile,sline);
                        std::istringstream ss(sline);
                        std::string buf;
                        std::vector<std::string> token;
                        while(ss >> buf) token.push_back(buf);
                        int left_idx = stoi(token[1])-1; int left_actual_index = tmp_mol.wH[left_idx].index;
                        int right_idx = stoi(token[2])-1;int right_actual_index = tmp_mol.wH[right_idx].index;
                        tmp_mol.wH[left_idx].bonded_atom.push_back(right_idx);
                        tmp_mol.wH[right_idx].bonded_atom.push_back(left_idx);

                        if( !tmp_mol.wH[left_idx].if_H() && !tmp_mol.wH[right_idx].if_H()) //both not hydrogen
                        {
                            tmp_mol.noH[left_actual_index].bonded_atom.push_back(right_actual_index);
                            tmp_mol.noH[right_actual_index].bonded_atom.push_back(left_actual_index);
                        }
                        else if(tmp_mol.wH[left_idx].if_H() && !tmp_mol.wH[right_idx].if_H()) //left hydrogen, right not H
                        {
                            tmp_mol.noH[right_actual_index].carryH = true;
                            tmp_mol.noH[right_actual_index].q += tmp_mol.wH[left_idx].q;

                        }
                        else if(!tmp_mol.wH[left_idx].if_H() && tmp_mol.wH[right_idx].if_H()) //left not hydrogen,right hydrogen
                        {
                            tmp_mol.noH[left_actual_index].carryH = true;
                            tmp_mol.noH[left_actual_index].q += tmp_mol.wH[right_idx].q;
                        }
                    }
                    in_mol.push_back(tmp_mol);
                    break;
                   
                }
            } 
            
        }
    }
    infile.close();

}
