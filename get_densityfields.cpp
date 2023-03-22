#include "include/traj.h"
void get_densityfields(float currentj[2][3][n_space_divz][n_space_divy][n_space_divx],
                       float np[2][n_space_divz][n_space_divy][n_space_divx],
                       float npt[n_space_divz][n_space_divy][n_space_divx],
                       int nt[2], float KEtot[2], float posL[3], float posH[3], float dd[3],
                       float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd],
                       float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd],
                       int q[2][n_partd], float dt[2], int mp[2], int n_part[2],
                       float jc[3][n_space_divz][n_space_divy][n_space_divx])
{
    // find number of particle and current density fields
    // set fields=0 in preparation
    unsigned int n;

    fill(reinterpret_cast<float *>(currentj), reinterpret_cast<float *>(currentj) + n_cells * 6, 0.f);
    fill(reinterpret_cast<float *>(np), reinterpret_cast<float *>(np) + n_cells * 2, 0); // Could split into threads.

    nt[0] = 0;
    nt[1] = 0;
    KEtot[0] = 0;
    KEtot[1] = 0;
    // set limits beyond which particle is considered as "lost"
    static float xld = posL[0] + dd[0] * 1;
    static float yld = posL[1] + dd[1] * 1;
    static float zld = posL[2] + dd[2] * 1;
    static float xhd = posH[0] - dd[0] * 1;
    static float yhd = posH[1] - dd[1] * 1;
    static float zhd = posH[2] - dd[2] * 1;
//#pragma omp parallel for
    for (int p = 0; p < 2; p++)
    {
        //    cout << "n_part[" << p << "]=" << n_part[p] << endl;
        for (n = 0; n < n_part[p]; n++)
        {
            // replace out of bounds particles and reduce the number of particles.
            if ((pos1x[p][n] <= xld) | (pos1y[p][n] <= yld) | (pos1z[p][n] <= zld) |
                (pos1x[p][n] >= xhd) | (pos1y[p][n] >= yhd) | (pos1z[p][n] >= zhd) |
                (pos0x[p][n] <= xld) | (pos0y[p][n] <= yld) | (pos0z[p][n] <= zld) |
                (pos0x[p][n] >= xhd) | (pos0y[p][n] >= yhd) | (pos0z[p][n] >= zhd))
            /* |
            (isnan(pos1x[p][n])) | (isnan(pos1y[p][n])) | (isnan(pos1z[p][n])) |
            (isnan(pos0x[p][n])) | (isnan(pos0y[p][n])) | (isnan(pos0z[p][n])))*/
            {
                int last = --n_part[p];
                //               cout << "n=" << n << ", npart[" << p << "]=" << n_part[p] << ":" << pos0x[p][n] << ", " << pos0y[p][n] << ", " << pos0z[p][n] << endl;
                pos0x[p][n] = pos0x[p][last];
                pos0y[p][n] = pos0y[p][last];
                pos0z[p][n] = pos0z[p][last];
                pos1x[p][n] = pos1x[p][last];
                pos1y[p][n] = pos1y[p][last];
                pos1z[p][n] = pos1z[p][last];
                q[p][n] = q[p][last];
                q[p][last] = 0;
                //                cout << "n=" << n << ", npart[" << p << "]=" << n_part[p] << ":" << pos1x[p][n] << ", " << pos1y[p][n] << ", " << pos1z[p][n] << endl;
                n--;
            }
        }
    }

//#pragma omp parallel for
    for (int p = 0; p < 2; p++)
    {
        //      cout << "npp" << p << n_part[p] << endl;

        //#pragma omp distribute parallel for reduction(+ \
                                              : np[p], nt[p], currentj[p], KEtot[p])
        for (int n = 0; n < n_part[p]; ++n)
        {
            unsigned i, j, k, c;
            float fracx, fracy, fracz; // integer and fractional component of the position
            //           cout << " thread: " << omp_get_thread_num() << endl;
            /*
            if (n == 0)
            {
                cout << "n=" << n << ", p=" << p << endl;
                cout << pos1x[p][n] << ", ";
                cout << pos1y[p][n] << ", ";
                cout << pos1z[p][n] << "\n ";
            }
            */
            float cx = (pos1x[p][n] - posL[0]) / dd[0];
            float cy = (pos1y[p][n] - posL[1]) / dd[1];
            float cz = (pos1z[p][n] - posL[2]) / dd[2];

            i = (unsigned int) cx; fracx = cx - i;
            j = (unsigned int) cy; fracy = cy - j;
            k = (unsigned int) cz; fracz = cz - k;

            // I'm not sure whether this algorithm is correct, but it seems legit
            // Only drawback is that it is not exact (a tiny fraction of charge leaks)
            float fracx1 = 1 - fracx, fracy1 = 1 - fracy, fracz1 = 1 - fracz, charge = q[p][n];
// The problem with any distribution is that it will end up pushing itself, although overall it still leads to slightly better stability
            #ifdef LinearDistribution // Distribute the particle according to 1/r, not defining this will distribute according to 1/r2
            float coeffs[8] = {
                fracx1 * fracy1 * fracz1,
                fracx  * fracy1 * fracz1,
                fracx1 * fracy  * fracz1,
                fracx  * fracy  * fracz1,
                fracx1 * fracy1 * fracz ,
                fracx  * fracy1 * fracz ,
                fracx1 * fracy  * fracz ,
                fracx  * fracy  * fracz  
            };
            #pragma unroll(8)
            for (int i = 0; i < 8; ++i) coeffs[i] *= charge;
            #else
            bool hasZero = (fracx == 0 || fracx1 == 0) && (fracy == 0 || fracy1 == 0) && (fracz == 0 || fracz1 == 0);
            float coeffs[8]; //c000, c100, c010, c110, c001, c101, c011, c111;
            if(hasZero){
                // deduce cell/coefficient index where the particle is precisely at
                unsigned int idx = ((fracx1 == 0) << 0) + ((fracy1 == 0) << 1) + ((fracz1 == 0) << 2);
                fill(coeffs, coeffs + 8, 0.f);
                coeffs[idx] = 1.f;
            } else {
                coeffs[0] = 1.f / (powf(fracx , 2) + powf(fracy , 2) + powf(fracz , 2));
                coeffs[1] = 1.f / (powf(fracx1, 2) + powf(fracy , 2) + powf(fracz , 2));
                coeffs[2] = 1.f / (powf(fracx , 2) + powf(fracy1, 2) + powf(fracz , 2));
                coeffs[3] = 1.f / (powf(fracx1, 2) + powf(fracy1, 2) + powf(fracz , 2));
                coeffs[4] = 1.f / (powf(fracx , 2) + powf(fracy , 2) + powf(fracz1, 2));
                coeffs[5] = 1.f / (powf(fracx1, 2) + powf(fracy , 2) + powf(fracz1, 2));
                coeffs[6] = 1.f / (powf(fracx , 2) + powf(fracy1, 2) + powf(fracz1, 2));
                coeffs[7] = 1.f / (powf(fracx1, 2) + powf(fracy1, 2) + powf(fracz1, 2));
            }
            float total = coeffs[0] + coeffs[1] + coeffs[2] + coeffs[3] + coeffs[4] + coeffs[5] + coeffs[6] + coeffs[7];
            total = charge / total; // Multiply r... by charge, ie total /= charge.
            // Then take the reciprocal for multiplication (faster)
            #pragma unroll(8)
            for (int i = 0; i < 8; ++i) coeffs[i] *= total;
            #endif
            // number of charge (in units of 1.6e-19 C) in each cell
            np[p][k    ][j    ][i    ] += coeffs[0];
            np[p][k    ][j    ][i + 1] += coeffs[1];
            np[p][k    ][j + 1][i    ] += coeffs[2];
            np[p][k    ][j + 1][i + 1] += coeffs[3];
            np[p][k + 1][j    ][i    ] += coeffs[4];
            np[p][k + 1][j    ][i + 1] += coeffs[5];
            np[p][k + 1][j + 1][i    ] += coeffs[6];
            np[p][k + 1][j + 1][i + 1] += coeffs[7];
            nt[p] += q[p][n];
            //cout << coeffs[0] + coeffs[1] + coeffs[2] + coeffs[3] + coeffs[4] + coeffs[5] + coeffs[6] + coeffs[7] << endl;
            // current density p=0 electron j=nev in each cell n in units 1.6e-19 C m/s
            float dx = pos1x[p][n] - pos0x[p][n], dy = pos1y[p][n] - pos0y[p][n], dz = pos1z[p][n] - pos0z[p][n];
            float vx = dx / dt[p], vy = dy / dt[p], vz = dz / dt[p];

            currentj[p][0][k    ][j    ][i    ] += coeffs[0] * vx;
            currentj[p][0][k    ][j    ][i + 1] += coeffs[1] * vx;
            currentj[p][0][k    ][j + 1][i    ] += coeffs[2] * vx;
            currentj[p][0][k    ][j + 1][i + 1] += coeffs[3] * vx;
            currentj[p][0][k + 1][j    ][i    ] += coeffs[4] * vx;
            currentj[p][0][k + 1][j    ][i + 1] += coeffs[5] * vx;
            currentj[p][0][k + 1][j + 1][i    ] += coeffs[6] * vx;
            currentj[p][0][k + 1][j + 1][i + 1] += coeffs[7] * vx;

            currentj[p][1][k    ][j    ][i    ] += coeffs[0] * vy;
            currentj[p][1][k    ][j    ][i + 1] += coeffs[1] * vy;
            currentj[p][1][k    ][j + 1][i    ] += coeffs[2] * vy;
            currentj[p][1][k    ][j + 1][i + 1] += coeffs[3] * vy;
            currentj[p][1][k + 1][j    ][i    ] += coeffs[4] * vy;
            currentj[p][1][k + 1][j    ][i + 1] += coeffs[5] * vy;
            currentj[p][1][k + 1][j + 1][i    ] += coeffs[6] * vy;
            currentj[p][1][k + 1][j + 1][i + 1] += coeffs[7] * vy;
        
            currentj[p][2][k    ][j    ][i    ] += coeffs[0] * vz;
            currentj[p][2][k    ][j    ][i + 1] += coeffs[1] * vz;
            currentj[p][2][k    ][j + 1][i    ] += coeffs[2] * vz;
            currentj[p][2][k    ][j + 1][i + 1] += coeffs[3] * vz;
            currentj[p][2][k + 1][j    ][i    ] += coeffs[4] * vz;
            currentj[p][2][k + 1][j    ][i + 1] += coeffs[5] * vz;
            currentj[p][2][k + 1][j + 1][i    ] += coeffs[6] * vz;
            currentj[p][2][k + 1][j + 1][i + 1] += coeffs[7] * vz;
            /*
            currentj[p][0][k][j][i] += q[p][n] * (pos1x[p][n] - pos0x[p][n]) / dt[p];
            currentj[p][1][k][j][i] += q[p][n] * (pos1y[p][n] - pos0y[p][n]) / dt[p];
            currentj[p][2][k][j][i] += q[p][n] * (pos1z[p][n] - pos0z[p][n]) / dt[p];
            */

            KEtot[p] += dx * dx + dy * dy + dz * dz;
        }
        KEtot[p] *= 0.5 * mp[p] / (e_charge_mass * dt[p] * dt[p]);
        KEtot[p] *= r_part_spart; // as if these particles were actually samples of the greater thing
        //      cout <<maxk <<",";
    }
    //  cout << "get_density";
//#pragma omp parallel for simd
    for (unsigned int i = 0; i < n_cells * 3; i++)
    {
        (reinterpret_cast<float *>(jc))[i] = (reinterpret_cast<float *>(currentj[0]))[i] + (reinterpret_cast<float *>(currentj[1]))[i];
    }
    //#pragma omp  parallel for simd
    for (unsigned int i = 0; i < n_cells; i++)
    {
        (reinterpret_cast<float *>(npt))[i] = (reinterpret_cast<float *>(np[0]))[i] + (reinterpret_cast<float *>(np[1]))[i];
    }
}