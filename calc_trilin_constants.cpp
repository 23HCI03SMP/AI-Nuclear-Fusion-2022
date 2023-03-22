#include "include/traj.h"
// calculate constants for each cell for trilinear interpolation
void calc_trilin_constants(float E[3][n_space_divz][n_space_divy][n_space_divx],
                           float Ea[n_space_divz][n_space_divy][n_space_divx][3][ncoeff],
                           float dd[3], float posL[3])
{
    float dV = -dd[0] * dd[1] * dd[2];
    const float dV1 = 1 / dV;
    float *E_flat = &E[0][0][0][0];
    //const int c_skip = n_cells;
    const int i_skip = 1, j_skip = n_space_divx * i_skip, k_skip = n_space_divy * j_skip, c_skip = n_space_divz * k_skip;
    const int ij_skip = i_skip + j_skip, jk_skip = j_skip + k_skip, ik_skip = i_skip + k_skip, ijk_skip = i_skip + j_skip + k_skip;
    int E_idx = 0;
    for (unsigned int k = 0; k < n_space_divz - 1; k++)
    {
        const float z0 = k * dd[2] + posL[2];
        const float z1 = z0 + dd[2];
        for (unsigned int j = 0; j < n_space_divy - 1; j++)
        {
            const float y0 = j * dd[1] + posL[1];
            const float y1 = y0 + dd[1];
            const float y0z0 = y0 * z0;
            const float y0z1 = y0 * z1;
            const float y1z0 = y1 * z0;
            const float y1z1 = y1 * z1;
//#pragma simd
            for (unsigned int i = 0; i < n_space_divx - 1; i++)
            {
                const float x0 = i * dd[0] + posL[0];
                const float x1 = x0 + dd[0];
                const float x0y0z0 = x0 * y0z0;
                const float x0y0z1 = x0 * y0z1;
                const float x0y1z0 = x0 * y1z0;
                const float x0y1z1 = x0 * y1z1;
                const float x1y0z0 = x1 * y0z0;
                const float x1y0z1 = x1 * y0z1;
                const float x1y1z0 = x1 * y1z0;
                const float x1y1z1 = x1 * y1z1;
                const float x0y0 = x0 * y0;
                const float x0y1 = x0 * y1;
                const float x1y0 = x1 * y0;
                const float x1y1 = x1 * y1;
                const float x0z0 = x0 * z0;
                const float x0z1 = x0 * z1;
                const float x1z0 = x1 * z0;
                const float x1z1 = x1 * z1;
                for (int c = 0, offset = E_idx; c < 3; ++c, offset += c_skip)
                {
                    const float c000 = E_flat[offset]; //E[c][k][j][i];
                    const float c100 = E_flat[offset + i_skip]; //E[c][k][j][i1];
                    const float c001 = E_flat[offset + k_skip]; //E[c][k1][j][i];
                    const float c010 = E_flat[offset + j_skip]; //E[c][k][j1][i];
                    const float c011 = E_flat[offset + jk_skip]; //E[c][k1][j1][i];
                    const float c101 = E_flat[offset + ik_skip]; //E[c][k1][j][i1];
                    const float c110 = E_flat[offset + ij_skip]; //E[c][k][j1][i1];
                    const float c111 = E_flat[offset + ijk_skip]; //E[c][k1][j1][i1];
                    //  -c000 x1y1z1 +c001 x1y1z0  +           c010  x1y0z1 -       c011     x1y0z0 +   c100         x0y1z1 - c101           x0y1z0  -    c110        x0y0z1  + c111      x0y0z0
                    Ea[k][j][i][c][0] = (-c000 * x1y1z1 + c001 * x1y1z0 + c010 * x1y0z1 - c011 * x1y0z0 + c100 * x0y1z1 - c101 * x0y1z0 - c110 * x0y0z1 + c111 * x0y0z0) * dV1;
                    // x
                    Ea[k][j][i][c][1] = ((c000 - c100) * y1z1 + (-c001 + c101) * y1z0 + (-c010 + c110) * y0z1 + (c011 - c111) * y0z0) * dV1;
                    // y
                    Ea[k][j][i][c][2] = ((c000 - c010) * x1z1 + (-c001 + c011) * x1z0 + (-c100 + c110) * x0z1 + (c101 - c111) * x0z0) * dV1;
                    // z
                    Ea[k][j][i][c][3] = ((c000 - c001) * x1y1 + (-c010 + c011) * x1y0 + (-c100 + c101) * x0y1 + (c110 - c111) * x0y0) * dV1;
                    // xy
                    Ea[k][j][i][c][4] = ((-c000 + c010 + c100 - c110) * z1 + (c001 - c011 - c101 + c111) * z0) * dV1;
                    // xz
                    Ea[k][j][i][c][5] = ((-c000 + c001 + c100 - c101) * y1 + (c010 - c011 - c110 + c111) * y0) * dV1;
                    // yz
                    Ea[k][j][i][c][6] = ((-c000 + c001 + c010 - c011) * x1 + (c100 - c101 - c110 + c111) * x0) * dV1;
                    Ea[k][j][i][c][7] = (c000 - c001 - c010 + c011 - c100 + c101 + c110 - c111) * dV1;
                }
                ++E_idx;
            }
            ++E_idx; // We only went to n_cellx - 1, so we need to + n_cellx
        }
        E_idx += j_skip; // Why does adding this line instantly make this 20% slower;-;
    }
}
