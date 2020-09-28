#include "RLDOCK_COORDINATE.h"


double COORDINATE::dis(class COORDINATE& b)
{
    double dis=(this->x-b.x)*(this->x-b.x)+(this->y-b.y)*(this->y-b.y)+(this->z-b.z)*(this->z-b.z);
    dis=sqrt(dis);
    return dis;
}

string ATOM::output(int true_index)
{
    ostringstream ost;
    ost<<fixed<<setprecision(4)<<right<<setw(4)<<true_index<<" "<<left<<setw(8)<<name<<right<<setw(11)<<setprecision(4)<<x<<setw(11)<<y<<setw(11)<<z<<" "<<left<<setw(7)<<type<<" "<<resi_index<<" "<<setw(4)<<resi_name<<right<<setw(11)<<q<<endl;
    return ost.str();

}

void ATOM::mol2_input(string instr)
{
    
        std::istringstream ss(instr);
        std::string buf;
        std::vector<std::string> token;
        while(ss >> buf) token.push_back(buf);
        index = stoi(token[0]);
        name = token[1];
        x = stod(token[2]);
        y = stod(token[3]);
        z = stod(token[4]);
        type = token[5];
        resi_index = stoi(token[6]);
        resi_name = token[7];
        q = stod(token[8]);
}

bool ATOM::if_H()
{
    if(this->element != "H" && this->element != "h")
    {return false;}
    else {return true;}
}