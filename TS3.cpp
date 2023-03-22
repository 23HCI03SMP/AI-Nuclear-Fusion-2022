/* TS3.cpp
This contains the main loop for the program. Most of the initialization occurs here, and time steps are iterated through.
For settings (as to what to calculate, eg. E / B field, E / B force) go to the defines in include/traj.h
*/
#include <complex>
#include "include/traj.h"
#include <fftw3.h>
// sphere
int n_part[3] = {n_parte, n_partd, n_parte + n_partd}; // 0,number of "super" electrons, electron +deuteriom ions, total
unsigned int n_space_div[3] = {n_space_divx, n_space_divy, n_space_divz};
unsigned int n_space_div2[3] = {n_space_divx2, n_space_divy2, n_space_divz2};

void log_headers(){
    logger.write("time_large");
    logger.write("time_small");
    logger.write("dtchanged");
    logger.write("ncalc_ele");
    logger.write("ncalc_deut");
    logger.write("dt_ele");
    logger.write("dt_deut");
    logger.write("t_sim");
    logger.write("ne");
    logger.write("ni");
    logger.write("KEtot_ele");
    logger.write("KEtot_deut");
    logger.write("Ele_pot");
    logger.write("Mag_pot");
    logger.write("E_tot");
    logger.newline();
}
void log_entry(int i_time, int ntime, int cdt, int total_ncalc[2], float dt[2], double t, int nt[2], float KEtot[2], float U[2]){
    logger.write(i_time);
    logger.write(ntime);
    logger.write(cdt);
    logger.write(total_ncalc[0]);
    logger.write(total_ncalc[1]);
    logger.write(dt[0]);
    logger.write(dt[1]);
    logger.write(t);
    logger.write(nt[0]);
    logger.write(nt[1]);
    logger.write(KEtot[0]);
    logger.write(KEtot[1]);
    logger.write(U[0]);
    logger.write(U[1]);
    logger.write(KEtot[0] + KEtot[1] + U[0] + U[1]);
    logger.newline();
}
int main()
{
    // Fast printing
    cin.tie(NULL);
    //ios_base::sync_with_stdio(false);
    cout << "Output dir: " << outpath << "\n";
    bool check = std::filesystem::create_directory(outpath);

    timer.mark(); // Yes, 3 time marks. The first is for the overall program dt
    timer.mark(); // The second is for compute_d_time
    timer.mark(); // The third is for start up dt

    omp_set_num_threads(nthreads);

    double t = 0;
    float Bmax = Bmax0;
    float Emax = Emax0;
    const unsigned int n_cells = n_space_divx * n_space_divy * n_space_divz;

    // position of particle and velocity: stored as 2 positions at slightly different times
    /** CL: Ensure that pos0/1.. contain multiple of 64 bytes, ie. multiple of 16 floats **/
    auto *pos0x = reinterpret_cast<float(&)[2][n_partd]>(*((float *)_aligned_malloc(sizeof(float) * n_partd * 2, 4096))); //new float[2][n_partd];
    auto *pos0y = reinterpret_cast<float(&)[2][n_partd]>(*((float *)_aligned_malloc(sizeof(float) * n_partd * 2, 4096))); //new float[2][n_partd];
    auto *pos0z = reinterpret_cast<float(&)[2][n_partd]>(*((float *)_aligned_malloc(sizeof(float) * n_partd * 2, 4096))); //new float[2][n_partd];
    auto *pos1x = reinterpret_cast<float(&)[2][n_partd]>(*((float *)_aligned_malloc(sizeof(float) * n_partd * 2, 4096))); //new float[2][n_partd];
    auto *pos1y = reinterpret_cast<float(&)[2][n_partd]>(*((float *)_aligned_malloc(sizeof(float) * n_partd * 2, 4096))); //new float[2][n_partd];
    auto *pos1z = reinterpret_cast<float(&)[2][n_partd]>(*((float *)_aligned_malloc(sizeof(float) * n_partd * 2, 4096))); //new float[2][n_partd];

    //    charge of particles
    auto *q = new int[2][n_partd]; // charge of each particle +1 for H,D or T or -1 for electron can also be +2 for He for example
    auto *m = new int[2][n_partd]; // mass of of each particle not really useful unless we want to simulate many different types of particles

    // reduced particle position dataset for printing/plotting
    auto *posp = new float[2][n_output_part][3];
    auto *KE = new float[2][n_output_part];

    /** CL: Ensure that Ea/Ba contain multiple of 64 bytes, ie. multiple of 16 floats **/
    auto *E = reinterpret_cast<float(&)[3][n_space_divz][n_space_divy][n_space_divx]>(*fftwf_alloc_real(3*n_cells));//new float[3][n_space_divz][n_space_divy][n_space_divx];
    auto *Ee = new float[3][n_space_divz][n_space_divy][n_space_divx];
    float *Ea1 = (float *)_aligned_malloc(sizeof(float) * n_cells * 3 * ncoeff, 4096); //auto *Ea = new float[3][ncoeff][n_space_divz][n_space_divy][n_space_divx];
    auto *Ea = reinterpret_cast<float(&)[n_space_divz][n_space_divy][n_space_divx][3][ncoeff]>(*Ea1);
    //auto Ex = new float[n_partd];
    //auto Ey = new float[n_partd];
    //auto Ez = new float[n_partd];
    //   auto *Einterpolated= new float[(n_space_divz-1)*4+1][(n_space_divy-1)*4+1][(n_space_divx-1)*4+1][3];
    //   auto *Vfield=new float[n_space_divz][n_space_divy][n_space_divx];

    auto *B = reinterpret_cast<float(&)[3][n_space_divz][n_space_divy][n_space_divx]>(*fftwf_alloc_real(3*n_cells));//new float[3][n_space_divz][n_space_divy][n_space_divx];
    auto *Be = new float[3][n_space_divz][n_space_divy][n_space_divx];
    float *Ba1 = (float *)_aligned_malloc(sizeof(float) * n_cells * 3 * ncoeff, 4096); //auto *Ba = new float[3][ncoeff][n_space_divz][n_space_divy][n_space_divx];
    auto *Ba = reinterpret_cast<float(&)[n_space_divz][n_space_divy][n_space_divx][3][ncoeff]>(*Ba1);
    //auto Bx = new float[n_partd];
    //auto By = new float[n_partd];
    //auto Bz = new float[n_partd];
    //    auto *Afield= new float[n_space_divz][n_space_divy][n_space_divx][3]; // x,y,z components

    auto *V = reinterpret_cast<float(&)[n_space_divz][n_space_divy][n_space_divx]>(*fftwf_alloc_real(n_cells));

    auto *np = new float[2][n_space_divz][n_space_divy][n_space_divx];
    auto *npt = new float[n_space_divz][n_space_divy][n_space_divx];
    int nt[2] = {0, 0};
    float KEtot[2] = {0, 0};
    auto *currentj = new float[2][3][n_space_divz][n_space_divy][n_space_divx];
    auto *jc = new float[3][n_space_divz][n_space_divy][n_space_divx];
    float U[2] = {0, 0};

    ofstream E_file, B_file;

    log_headers();

    cout << std::scientific;
    cout.precision(3);
    cerr << std::scientific;
    cerr.precision(3);
    cout << "float size=" << sizeof(float) << ", "
         << "int32_t size=" << sizeof(int32_t) << ", "
         << "int size=" << sizeof(int) << endl;
    int ncalc[2] = {md_me*1, 1};
    int total_ncalc[2] = {0, 0};
    // particle 0 - electron, particle 1 deuteron
    // set plasma parameters
    int mp[2] = {1, 1835 * 2};
    // float mp[2]= {9.10938356e-31,3.3435837724e-27}; //kg
    int qs[2] = {-1, 1};        // Sign of charge
    float Temp[2] = {1e5, 1e5}; // in K convert to eV divide by 1.160451812e4

    // initial bulk electron, ion velocity
    float v0[2][3] = {{0, 0, 0/*1e6*/}, {0, 0, 0}};

    // typical dimensions of a cell
    float a0 = 2e-3;

    float r0 = (n_space / 8) * a0; // if sphere this is the radius
    float area = 4 * pi * r0 * r0;
    float volume = 4 / 3 * pi * r0 * r0 * r0;
    // float volume=((posHp[0]-posLp[0])*(posHp[1]-posLp[1])*(posHp[2]-posLp[2]))); //cube

    // calculated plasma parameters
    float Density_e = n_partd / volume * r_part_spart;
    // float initial_current=Density_e*e_charge*v0[0][2]*area;
    // float       Bmax=initial_current*2e-7/a0*10;

    float plasma_freq = sqrt(Density_e * e_charge * e_charge_mass / (mp[0] * epsilon0)) / (2 * pi);
    float plasma_period = 1 / plasma_freq;
    float Debye_Length = sqrt(epsilon0 * kb * Temp[0] / (Density_e * e_charge * e_charge));
    float vel_e = sqrt(kb * Temp[0] / (mp[0] * e_mass));
    float Tv = a0 / vel_e; // time for electron to move across 1 cell
    float Tcyclotron = 2.0 * pi * mp[0] / (e_charge_mass * Bmax);
    float TDebye = Debye_Length / vel_e;
    float TE = sqrt(2 * a0 / e_charge_mass / Emax);
    // set time step to allow electrons to gyrate if there is B field or to allow electrons to move slowly throughout the plasma distance


    float dt[2];
    dt[0] = 4*min(min(min(TDebye, min(Tv / ncalc[0], Tcyclotron) / 4), plasma_period / ncalc[0] / 4), TE / ncalc[0]) / 2; // electron should not move more than 1 cell after ncalc*dt and should not make more than 1/4 gyration and must calculate E before the next 1/4 plasma period
    //dt[0] /= 2.f;
    //Bmax *= 2;
    //Emax *= 4;
    dt[1] = dt[0] * (ncalc[0] / ncalc[1]);
    //  float mu0_4pidt[2]= {mu0_4pi/dt[0],mu0_4pi/dt[1]};
    cout << "v0 electron = " << v0[0][0] << "," << v0[0][1] << "," << v0[0][2] << endl;
    //   cout <<"Initial Current = "<<initial_current<<endl;
    //   cout <<"Initial Bmax = "<<initial_current*2e-7/a0<<endl;

    cout << "Start up dt = " << timer.replace() << "s\n";
#define generateRandom
    #ifdef generateRandom
    // set initial positions and velocity
    float sigma[2] = {sqrt(kb * Temp[0] / (mp[0] * e_mass)), sqrt(kb * Temp[1] / (mp[1] * e_mass))};
    long seed;
    gsl_rng *rng;                        // random number generator
    rng = gsl_rng_alloc(gsl_rng_rand48); // pick random number generator
    seed = 1670208073;//time(NULL);
    cout << "seed=" << seed << "\n";
    gsl_rng_set(rng, seed); // set seed

    /*float posLp[3], posHp[3];
    for (int c = 0; c < 3; c++)
    {
        posLp[c] = -a0 * n_space / 4;
        posHp[c] = a0 * n_space / 4;
    }*/

    for (int p = 0; p < 2; p++)
    {
//#pragma omp parallel for reduction(+ \
                                   : nt)
        for (int n = 0; n < n_partd; n++)
        {

            // spherical plasma radius is 1/8 of total extent.
            float r = r0 * pow(gsl_ran_flat(rng, 0, 1), 0.3333333333);
            //if (p == 0) r += n_space / 8 * a0;
            double x, y, z;
            gsl_ran_dir_3d(rng, &x, &y, &z);
            pos0x[p][n] = r * x;
            pos1x[p][n] = pos0x[p][n] + (gsl_ran_gaussian(rng, sigma[p]) + v0[p][0]) * dt[p];
            pos0y[p][n] = r * y;
            pos1y[p][n] = pos0y[p][n] + (gsl_ran_gaussian(rng, sigma[p]) + v0[p][1]) * dt[p];
            pos0z[p][n] = r * z;
            pos1z[p][n] = pos0z[p][n] + (gsl_ran_gaussian(rng, sigma[p]) + v0[p][2]) * dt[p];
            //          if (n==0) cout << "p = " <<p <<", sigma = " <<sigma[p]<<", temp = " << Temp[p] << ",mass of particle = " << mp[p] << dt[p]<<endl;
            q[p][n] = qs[p];
            m[p][n] = mp[p];
            nt[p] += q[p][n];
        }
    }

    gsl_rng_free(rng); // dealloc the rng

    // get limits and spacing of Field cells
    #else
    generateParticles(a0, r0, qs, mp, pos0x, pos0y, pos0z, pos1x, pos1y, pos1z, q, m, nt);
    n_part[0] = abs(nt[0]);
    n_part[1] = abs(nt[1]);
    n_part[2] = n_part[0] + n_part[1];
    #endif
    generateField(Ee, Be);
    cout << "Set initial random positions: " << timer.replace() << "s\n";
    float posL[3], posH[3], posL2[3], dd[3];
    // set spacing between cells
    for (int c = 0; c < 3; c++)
        dd[c] = a0;
    // set position of centers of the cells at extreme ends
    for (int c = 0; c < 3; c++)
    {
        posL[c] = -dd[c] * (n_space_div[c] - 1.0) / 2.0;
        posH[c] = dd[c] * (n_space_div[c] - 1.0) / 2.0;
        posL2[c] = -dd[c] * (n_space_div[c]);
        cout << posL[c] << "," << posH[c] << "," << dd[c] << endl;
    }

    //unsigned int ci[7] = {n_partd, n_cells, n_space_divx, n_space_divy, n_space_divz, 0, n_space * 2 - 1};
    //float cf[11]] = {0, 0, posL[0], posH[0], posL[1], posH[1], posL[2], posH[2], dd[0], dd[1], dd[2]};
    unsigned int ci[2] = {n_partd, 0};
    float cf[2] = {0, 0};
    cl_set_build_options(posL, posH, dd);
    cl_start();

    // print initial conditions
    {
        cout << "electron Temp = " << Temp[0] << " K, electron Density = " << Density_e << " m^-3" << endl;
        cout << "Plasma Frequency(assume cold) = " << plasma_freq << " Hz, Plasma period = " << plasma_period << " s" << endl;
        cout << "Cyclotron period = " << Tcyclotron << " s, Time for electron to move across 1 cell = " << Tv << " s" << endl;
        cout << "Time taken for electron at rest to accelerate across 1 cell due to E = " << TE << " s" << endl;
        cout << "electron thermal velocity = " << vel_e << endl;
        cout << "dt = " << dt[0] << " s, Total time = " << dt[0] * ncalc[0] * ndatapoints * nc << ", s" << endl;
        cout << "Debye Length = " << Debye_Length << " m, initial dimension = " << a0 << " m" << endl;
        cout << "number of particle per cell = " << n_partd / (n_space * n_space * n_space) * 8 << endl;

        E_file.open("info.csv");
        E_file << ",X, Y, Z" << endl;
        E_file << "Data Origin," << posL[0] << "," << posL[1] << "," << posL[0] << endl;
        E_file << "Data Spacing," << dd[0] << "," << dd[1] << "," << dd[2] << endl;
        E_file << "Data extent x, 0," << n_space - 1 << endl;
        E_file << "Data extent y, 0," << n_space - 1 << endl;
        E_file << "Data extent z, 0," << n_space - 1 << endl;

        E_file << "electron Temp = ," << Temp[0] << ",K" << endl;
        E_file << "electron Density =," << Density_e << ",m^-3" << endl;
        E_file << "electron thermal velocity = ," << vel_e << endl;
        E_file << "Maximum expected B = ," << Bmax << endl;
        E_file << "Plasma Frequency(assume cold) = ," << plasma_freq << ", Hz" << endl;
        E_file << "Plasma period =," << plasma_period << ",s" << endl;
        E_file << "Cyclotron period =," << Tcyclotron << ",s" << endl;
        E_file << "Time for electron to move across 1 cell Tv =," << Tv << ",s" << endl;
        E_file << "Time for electron to move across 1 cell TE =," << TE << ",s" << endl;
        E_file << "time step between prints = ," << dt[0] * ncalc[0] * nc << ",s" << endl;
        E_file << "time step between EBcalc = ," << dt[0] * ncalc[0] << ",s" << endl;
        E_file << "dt =," << dt[0] << ",s" << endl;
        E_file << "Debye Length =," << Debye_Length << ",m" << endl;
        E_file << "Larmor radius =," << vel_e/(Bmax*e_charge_mass) << ",m" << endl;
        E_file << "cell size =," << a0 << ",m" << endl;
        E_file << "number of particles per cell = ," << n_partd / (n_space * n_space * n_space) << endl;
        E_file.close();
    }

    int i_time = 0;
    get_densityfields(currentj, np, npt, nt, KEtot, posL, posH, dd, pos1x, pos1y, pos1z, pos0x, pos0y, pos0z, q, dt, mp, n_part, jc);
    calcEBV(V, E, B, Ee, Be, npt, jc, dd, Emax, Bmax);
    calc_trilin_constants(E, Ea, dd, posL);
    calc_trilin_constants(B, Ba, dd, posL);
    #ifdef Uon_
    calcU(V, E, B, pos1x, pos1y, pos1z, posL, dd, n_part, q, U);
    #endif
    cout << i_time << "." << 0 << " (compute_time = " << timer.elapsed() << "s): ";
    cout << "dt = {" << dt[0] << " " << dt[1] << "}, t_sim = " << t << " s" << ", ne = " << nt[0] << ", ni = " << nt[1];
    cout << "\nKEtot e = " << KEtot[0] << ", KEtot i = " << KEtot[1] << ", Eele = " << U[0] << ", Emag = " << U[1] << ", Etot = " << KEtot[0] + KEtot[1] + U[0] + U[1] << " eV\n";
    sel_part_print(n_part, pos1x, pos1y, pos1z, pos0x, pos0y, pos0z, posp, KE, m, dt);
    save_files(i_time, n_space_div, posL, dd, t, np, currentj, V, E, B, KE, posp);    
    cout << "print data: " << timer.elapsed() << "s (no. of electron time steps calculated: " << 0 << ")\n";
    
    // Write everything to log
    log_entry(i_time, 0, 0, total_ncalc, dt, t, nt, KEtot, U);

    for (i_time = 1; i_time < ndatapoints; i_time++)
    {
        for (int ntime = 0; ntime < nc; ntime++)
        {
            /*
            for (int p = 0; p < 2; ++p){
                float dposx = -1.f, dposy = -1.f, dposz = -1.f, dpos = -1.f;
                for (int i = 0; i < n_part[p]; ++i)
                {
                    float dx = pos1x[p][i] - pos0x[p][i];
                    float dy = pos1y[p][i] - pos0y[p][i];
                    float dz = pos1z[p][i] - pos0z[p][i];
                    float dp = sqrt(dx * dx + dy * dy + dz * dz);
                    
                    if (dp > dpos) {
                        dpos = dp;
                        dposx = dx;
                        dposy = dy;
                        dposz = dz;
                    }
                    dpos = max(dpos, dp);
                    dposx = max(dposx, dx);
                    dposy = max(dposy, dy);
                    dposz = max(dposz, dz);
                }
                if(p != 0) continue;
                cout.precision(15);
                cout << "Fastest " << p << ": " << dpos << " (" << dposx/dd[0] << ", " << dposy/dd[1] << ", " << dposz/dd[2] << ")\n";
                cout.precision(3);
            }
            */
            timer.mark(); // For timestep
            // Work out motion
            timer.mark();
            for (int p = 0; p < 2; p++)
            {
                const float coef = (float)qs[p] * e_charge_mass / (float)mp[p] * dt[p] * 0.5f;
                #ifdef BFon_
                cf[0] = coef;
                #else
                cf[0] = 0;
                #endif
                #ifdef EFon_
                cf[1] = coef * dt[p]; // multiply by dt because of the later portion of cl code
                #else
                cf[1] = 0;
                #endif
                ci[1] = ncalc[p];
                ci[0] = n_part[p]; //
                                   //               cout << p << " Bconst=" << cf[0] << ", Econst=" << cf[1] << endl;
                // calculate the next position ncalc[p] times
                tnp(Ea1, Ba1, pos0x[p], pos0y[p], pos0z[p], pos1x[p], pos1y[p], pos1z[p], cf, ci);
                total_ncalc[p] += ncalc[p];
            }
            cout << "motion: " << timer.elapsed() << "s, ";

            t += dt[0] * ncalc[0];

            //  find number of particle and current density fields
            timer.mark();
            get_densityfields(currentj, np, npt, nt, KEtot, posL, posH, dd, pos1x, pos1y, pos1z, pos0x, pos0y, pos0z, q, dt, mp, n_part, jc);
            cout << "density: " << timer.elapsed() << "s, ";

            // find E field must work out every i,j,k depends on charge in every other cell
            timer.mark();
            // set externally applied fields this is inside time loop so we can set time varying E and B field
            // calcEeBe(Ee,Be,t);
            int cdt = calcEBV(V, E, B, Ee, Be, npt, jc, dd, Emax, Bmax);
            /* change time step if E or B too big*/
            if (cdt)
            {
                changedt(pos0x, pos0y, pos0z, pos1x, pos1y, pos1z, n_part);
                dt[0] /= 2;
                dt[1] /= 2;
                ncalc[0] *= 2;
                ncalc[1] *= 2;
                Emax *= 4;
                Bmax *= 2;
                // cout <<"\ndtchanged\n";
            }
            cout << "EBV: " << timer.elapsed() << "s, ";

            #ifdef Uon_
            // calculate the total potential energy U
            timer.mark();
            calcU(V, E, B, pos1x, pos1y, pos1z, posL, dd, n_part, q, U);
            cout << "U: " << timer.elapsed() << "s, ";
            #endif

            // calculate constants for each cell for trilinear interpolation
            timer.mark();
            calc_trilin_constants(E, Ea, dd, posL);
            calc_trilin_constants(B, Ba, dd, posL);
            cout << "trilin const: " << timer.elapsed() << "s";

            cout << "\n\n" << i_time << "." << ntime << " (compute_time = " << timer.elapsed() << "s): ";
            if (cdt) cout << "dtchanged\n";
            cout << "dt = {" << dt[0] << " " << dt[1] << "}, t_sim = " << t << " s" << ", ne = " << nt[0] << ", ni = " << nt[1];
            cout << "\nKEtot e = " << KEtot[0] << ", KEtot i = " << KEtot[1] << ", Eele = " << U[0] << ", Emag = " << U[1] << ", Etot = " << KEtot[0] + KEtot[1] + U[0] + U[1] << " eV\n";
            log_entry(i_time, ntime, cdt?1:0, total_ncalc, dt, t, nt, KEtot, U);
        }

        // print out all files for paraview
        timer.mark();
        sel_part_print(n_part, pos1x, pos1y, pos1z, pos0x, pos0y, pos0z, posp, KE, m, dt);
        save_files(i_time, n_space_div, posL, dd, t, np, currentj, V, E, B, KE, posp);
        cout << "print data: " << timer.elapsed() << "s (no. of electron time steps calculated: " << total_ncalc[0] << ")\n";
    }
    cout << "Overall execution time: " << timer.elapsed() << "s";
    logger.close();
    return 0;
}
