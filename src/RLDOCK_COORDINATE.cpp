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
        string::size_type position;
        position = type.find(".");
        element = type.substr(0,position);
}

bool ATOM::if_H()
{
    if(this->element != "H" && this->element != "h")
    {return false;}
    else {return true;}
}

void BOND::mol2_input(string instr)
{
    
    std::istringstream ss(instr);
    std::string buf;
    std::vector<std::string> token;
    while(ss >> buf) token.push_back(buf);
    this->index = stoi(token[0]);
    this->left_id = stoi(token[1]);
    this->right_id = stoi(token[2]);
    this->atom_type = token[3];
}
void BOND::update(int new_idx, int new_left,int new_right)
{
    this->index = new_idx;
    this->left_id = new_left;
    this->right_id = new_right;
}

void SUBSTRUCTURE::mol2_input(string instr)
{
    
    std::istringstream ss(instr);
    std::string buf;
    std::vector<std::string> token;
    while(ss >> buf) token.push_back(buf);
    this->index = stoi(token[0]);
    this->sub_name = token[1];
    this->root_atom = stoi(token[2]);
    this->sub_type = token[3];
    this->dict_type = stoi(token[4]);
    this->chain = token[5];
    this->chain_type = token[6];
    this->inter_bonds = stoi(token[7]);
    this->status = token[8];
}