#include "RLDOCK_PARAMETER.h"

void PARAMETER::SetSphere(string filename)
{
    ifstream spherefile(filename.c_str());
    string sline;
    while(getline(spherefile,sline))
    {
        std::istringstream ss(sline);
        std::string buf;
        std::vector<std::string> token;
        while(ss >> buf) token.push_back(buf);
        this->ux.push_back(stod(token[1]));
        this->vy.push_back(stod(token[2]));
        this->wz.push_back(stod(token[3]));
    }

}

void PARAMETER::SetAngle()
{
    
    int Ngamma = 360.0/gamma_step1+0.01;
    int index = 0;
    for(int ig =0; ig != Ngamma;ig++)
    {
        double dg = ig*(2*pi/Ngamma);
        double sing = sin(dg);
        double cosg = cos(dg);
        for(int i=0;i!=this->ux.size();i++)
        {
            this->aa.push_back(ux[i]*ux[i]+(vy[i]*vy[i]+wz[i]*wz[i])*cosg);
            this->bb.push_back(ux[i]*vy[i]*(1.-cosg)-wz[i]*sing);
            this->cc.push_back(ux[i]*wz[i]*(1.-cosg)+vy[i]*sing);
            this->dd.push_back(ux[i]*vy[i]*(1.-cosg)+wz[i]*sing);
            this->ee.push_back(vy[i]*vy[i]+(ux[i]*ux[i]+wz[i]*wz[i])*cosg);
            this->ff.push_back(vy[i]*wz[i]*(1.-cosg)-ux[i]*sing);
            this->gg.push_back(ux[i]*wz[i]*(1.-cosg)-vy[i]*sing);
            this->hh.push_back(vy[i]*wz[i]*(1.-cosg)+ux[i]*sing);
            this->ii.push_back(wz[i]*wz[i]+(ux[i]*ux[i]+vy[i]*vy[i])*cosg);
            
            index++;
        }
    }

}

void PARAMETER::UpdatAngle(int new_gamma_step1)
{
    this->gamma_step1 = new_gamma_step1;
    int Ngamma = 360.0/gamma_step1+0.01;
    int index = 0;
    int new_size = Ngamma*ux.size();
    this->aa.resize(new_size);
    this->bb.resize(new_size);
    this->cc.resize(new_size);
    this->dd.resize(new_size);
    this->ee.resize(new_size);
    this->ff.resize(new_size);
    this->gg.resize(new_size);
    this->hh.resize(new_size);
    this->ii.resize(new_size);
    for(int ig =0; ig != Ngamma;ig++)
    {
        double dg = ig*(2*pi/Ngamma);
        double sing = sin(dg);
        double cosg = cos(dg);
        for(int i=0;i!=ux.size();i++)
        {
            this->aa[index]=(ux[i]*ux[i]+(vy[i]*vy[i]+wz[i]*wz[i])*cosg);
            this->bb[index]=(ux[i]*vy[i]*(1.-cosg)-wz[i]*sing);
            this->cc[index]=(ux[i]*wz[i]*(1.-cosg)+vy[i]*sing);
            this->dd[index]=(ux[i]*vy[i]*(1.-cosg)+wz[i]*sing);
            this->ee[index]=(vy[i]*vy[i]+(ux[i]*ux[i]+wz[i]*wz[i])*cosg);
            this->ff[index]=(vy[i]*wz[i]*(1.-cosg)-ux[i]*sing);
            this->gg[index]=(ux[i]*wz[i]*(1.-cosg)-vy[i]*sing);
            this->hh[index]=(vy[i]*wz[i]*(1.-cosg)+ux[i]*sing);
            this->ii[index]=(wz[i]*wz[i]+(ux[i]*ux[i]+vy[i]*vy[i])*cosg);
            
            index++;
        }
    }

}

void PARAMETER::InputUpdate(string rf_in, string lf_in, string ref_in,string op_in, int thread_in)
{
    this->receptor_file = rf_in;
    this->ligand_file = lf_in;
    this->output_prefix = op_in;
    this->num_threads = thread_in;
    this->reference_file = ref_in;
    if(ref_in.length() == 0 ){this->reference_flag = false;}
    else{this->reference_flag = true;}
}

bool rank_less(const RANK_PARAMETER& i, const RANK_PARAMETER& j)
{
    return(i.val<j.val);
}

void Thread_Distribution(vector<THREAD_PARAMETER>& thread_p, int n_total, int n_piece)
{
    thread_p.resize(n_piece);
    int num_per_thread = n_total / n_piece;
    int num_plus = n_total % n_piece;

    for(unsigned int i = 1; i != n_piece; ++i)
    {
        thread_p[i-1].iend = thread_p[i-1].ibegin + num_per_thread - 1;
        if( i < num_plus)thread_p[i-1].iend ++;
        thread_p[i].ibegin = thread_p[i-1].iend + 1;
    }
    thread_p[n_piece-1].iend = n_total -1;
}