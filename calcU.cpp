#include "include/traj.h"
namespace {
    struct position {
        float x, y, z;
        unsigned int cell;
        int charge;
    };
    bool sortPositions(position const&lhs, position const&rhs){
        return lhs.cell < rhs.cell;
    }
}
void calcU(float V[n_space_divz][n_space_divy][n_space_divx],
            float E[3][n_space_divz][n_space_divy][n_space_divz], float B[3][n_space_divz][n_space_divy][n_space_divx],
            float posx[2][n_partd], float posy[2][n_partd], float posz[2][n_partd],
            float posL[3], float dd[3], int n_part[2], int q[2][n_partd], float out[2])
{
    const float dd0 = 1/dd[0], dd1 = 1/dd[1], dd2 = 1/dd[2];
    float EUtot = 0.f;
    {
        // Calculate electrical potential energy, from get_densityfields.cpp
        for (int p = 0; p < 2; p++)
        {
            for (int n = 0; n < n_part[p]; ++n)
            {
                float dx = (posx[p][n] - posL[0]) * dd0; // Get the cell positions in decimal
                float dy = (posy[p][n] - posL[1]) * dd1;
                float dz = (posz[p][n] - posL[2]) * dd2;
                unsigned int i = dx, j = dy, k = dz; // Round downwards
                            dx -= i; dy -= k; dz -= k; // and get the "fractional cell" values (ie. located at cell 5.1 -> 0.1)            
                // Perform trilinear interpolation
                float dx1 = 1 - dx;
                float c00 = V[k    ][j    ][i    ] * dx1 + V[k    ][j    ][i + 1] * dx;
                float c01 = V[k + 1][j    ][i    ] * dx1 + V[k + 1][j    ][i + 1] * dx;
                float c10 = V[k    ][j + 1][i    ] * dx1 + V[k    ][j + 1][i + 1] * dx;
                float c11 = V[k + 1][j + 1][i    ] * dx1 + V[k + 1][j + 1][i + 1] * dx;
                float c = (c00 * (1 - dy) + c10 * dy) * (1 - dz) + (c01 * (1 - dy) + c11 * dy) * dz;
                EUtot += c * q[p][n];
            }
        }
        EUtot *= 0.5f * e_charge;
        EUtot *= r_part_spart; // scale to target particles
    }
    // Calculate energy between particles when they are in the same cell
    #ifdef UE_cell
    static auto pos = new position[n_partd * 2];
    const int limit = n_part[0] + n_part[1];
    if (limit > 1)
    {
        // Set the position data
        int idx = 0;
        for (int p = 0; p < 2; ++p)
        {
            for (int n = 0; n < n_part[p]; ++n)
            {
                position* curr = &pos[idx++];
                curr->x = posx[p][n];
                curr->y = posy[p][n];
                curr->z = posz[p][n];
                curr->charge = q[p][n];
                unsigned int dx = (posx[p][n] - posL[0]) * dd0;
                unsigned int dy = (posy[p][n] - posL[1]) * dd1;
                unsigned int dz = (posz[p][n] - posL[2]) * dd2;
                curr->cell = (dz * n_space_divy + dy) * n_space_divx + dx;
            }
        }
        // Sort by cell number
        sort(pos, pos + limit, &sortPositions);
        bool invalid = false; int max_particles_percell = 0;
        float Etot = 0.f;
        int curr_cell = pos[0].cell;
        int start = 0; // First position with .cell == curr_cell
        // TODO: maybe do a mini FFT or something - this is N^2 and could be problematic.
        for (int n = 0; n < limit; ++n)
        {
            position curr = pos[n];
            if (curr_cell != curr.cell) { // Swap to new cell, reset stats
                max_particles_percell = max(max_particles_percell, n - start);
                start = n;
                curr_cell = curr.cell;
            } else {
                float to_add = 0.f;
                invalid |= (n - start > 10000);
                for (int i = start; i < min(n, start + 10000); ++i)
                { // "Add" the particle to the cell, hence add with all other previous particles
                    position other = pos[i];
                    float r = sqrtf(powf(other.x - curr.x, 2) + powf(other.y - curr.y, 2) + powf(other.z - curr.z, 2));
                    if (r == 0.f) continue;
                    to_add += other.charge / r;
                }
                Etot += to_add * curr.charge;
            }
        }
        if (invalid) cout << "Eel invalid (max " << max_particles_percell << " particles in one cell), ";
        Etot *= kc * e_charge * e_charge;
        Etot *= r_part_spart * r_part_spart ;  // as if it were scaled to that many particles
        EUtot += Etot;
    }
    #endif
    // Calculate energy stored in electric field - e0/2 E^2
    #ifdef UE_field
    {
        float E2tot = 0;
        float * E_1d = reinterpret_cast<float *>(E);
        for (int i = 0; i < n_cells * 3; ++i){
            float e = E_1d[i];
            E2tot += e * e; // Why can we do this? Because E^2 = Ex^2 + Ey^2 + Ez^2
        }
        E2tot *= 0.5f * epsilon0 * (dd[0] * dd[1] * dd[2]); // dV
        EUtot += E2tot;
    }
    #endif
    float BUtot = 0;
    #define UB_field
    #ifdef UB_field
    // Energy stored in magnetic field - 1/2 B^2/u0
    {
        float B2tot = 0;
        float * B_1d = reinterpret_cast<float *>(B);
        for (int i = 0; i < n_cells * 3; ++i){
            float b = B_1d[i];
            B2tot += b * b;
        }
        B2tot *= 0.5f / u0 * (dd[0] * dd[1] * dd[2]); // dV
        BUtot += B2tot;
    }
    #endif
    out[0] = EUtot / ev_to_j;
    out[1] = BUtot / ev_to_j;
}