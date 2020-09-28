#ifndef COORDINATE_H
#define COORDINATE_H
#include "main.h"
class COORDINATE
{
    public:
    int index;
    double x;
    double y;
    double z; 
    double rb0;

    double dis(class COORDINATE& b);  
    
};

class ATOM : public COORDINATE
{
    public:
    string name;
    string type;
    int resi_index;
    string resi_name;
    double q;

    string element;
    double rrr;
    double sb;
    int type_idx;

    int index_original;
    bool carryH=false;
    
    vector<int> bonded_atom;

    string output(int true_index);
    void mol2_input(string instr);

    bool if_H();

    
};
#endif