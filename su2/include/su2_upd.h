#ifndef SU2_UPD_H
#define SU2_UPD_H

#include"su2.h"

void randheat(double k, double *out);
void single_heatbath(Su2 *link, Su2 const * const staple);
void single_overrelaxation(Su2 *link, Su2 const * const staple);

#endif
