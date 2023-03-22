#include "include/traj.h"
void changedt(float pos0x[2][n_partd],float pos0y[2][n_partd], float pos0z[2][n_partd], float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd], int n_part[3])
{
    //recalculate pos0 for time step new time step of half
    for (int p = 0; p < 2; p++)
    {
        for (int n = 0; n < n_part[p]; n++)
        {
            pos0x[p][n] =  (pos1x[p][n] + pos0x[p][n]) * 0.5;
            pos0y[p][n] =  (pos1y[p][n] + pos0y[p][n]) * 0.5;
            pos0z[p][n] =  (pos1z[p][n] + pos0z[p][n]) * 0.5;
        }
    }
}