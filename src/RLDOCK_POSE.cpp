#include "RLDOCK_POSE.h"

void POSE::info_input(int iconf_in, int iatom_in, POCKET& pocket, vector<COORDINATE>& pose)
{
    iconf = iconf_in;
    iatom = iatom_in;
    position.resize(pose.size());
    for(int ilig=0;ilig!=pose.size();ilig++)
    {
        position[ilig].x = pocket.x + pose[ilig].x;
        position[ilig].y = pocket.y + pose[ilig].y;
        position[ilig].z = pocket.z + pose[ilig].z;
    }
}

bool greater_cluster(const POSE& i, const POSE& j)
{
    return (i.num_cluster > j.num_cluster);
}

bool less_eng (const POSE & p1, const POSE& p2)
{
    return (p1.eng < p2.eng);
}
bool less_eng_rerank (const POSE & p1, const POSE & p2)
{
    return (p1.eng_rerank < p2.eng_rerank);
}