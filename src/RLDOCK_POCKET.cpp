#include "RLDOCK_POCKET.h"
bool less_label_lj(const LABEL& l1, const LABEL& l2)
{
    return (l1.lj < l2.lj);
}

bool less_pocket_lj (const POCKET& p1, const POCKET& p2)
{
    return (p1.lj_site_min < p2.lj_site_min);
}