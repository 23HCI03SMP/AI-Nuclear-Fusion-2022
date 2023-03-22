#ifndef TRAJ_H_INCLUDED
#define TRAJ_H_INCLUDED
#define CL_HPP_TARGET_OPENCL_VERSION 300
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <omp.h>
#include <string>
#include <filesystem>
#include <CL/opencl.hpp>
#define trilinon_
#define Uon_ //whether to calculate the electric (V) potential and potential energy (U). Needs Eon to be enabled.
#define Eon_ //whether to calculate the electric (E) field
#define Bon_ //whether to calculate the magnetic (B) field
#define EFon_ //whether to apply electric force
#define BFon_ //whether to apply magnetic force
#define printDensity
#define printParticles
//#define printV //print out V
#define printB //print out B field
#define printE //print out E field
//#define FileIn //whether to load from input file (unused)
#define RamDisk //whether to use RamDisk

using namespace std;
//save file info - initialize filepath
#ifdef RamDisk
const string outpath = "R:\\Temp\\out\\";
#elif !defined(_WIN32) && (defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__)))
/* UNIX-style OS. ------------------------------------------- */
const string outpath = std::filesystem::temp_directory_path().string() + "/out/";
#else
const string outpath = std::filesystem::temp_directory_path().string() + "out/";
#endif

// technical parameters
constexpr int n_space = 64;                             // must be 2 to power of n
constexpr int n_partd = n_space * n_space * n_space * 8; // must be 2 to power of n
constexpr int n_parte = n_partd;
constexpr unsigned int ncoeff = 8;
extern cl::Context context_g;
extern cl::Device default_device_g;
extern cl::Program program_g;
constexpr int n_output_part = min(n_partd, 8192); // maximum number of particles to output to file
// const int nprtd=floor(n_partd/n_output_part);

constexpr int ndatapoints = 300; // total number of time steps to calculate
constexpr int nc = 10;          // number of times to calculate E and B between printouts
constexpr int md_me = 60;       // ratio of electron speed/deuteron speed at the same KE. Used to calculate electron motion more often than deuteron motion
constexpr int ncalc0[2] = {md_me, 1};
constexpr int nthreads = 8; // match with your CPU

// The maximum expected E and B fields. If fields go beyond this, the the time step, cell size etc will be wrong. Should adjust and recalculate.
//  maximum expected magnetic field
constexpr float Bmax0 = 1;
constexpr float Emax0 = 1e6;

constexpr float target_part = 1e10;
constexpr float r_part_spart = target_part / n_partd;//1e12 / n_partd; // ratio of particles per tracked "super" particle
//ie. the field of N particles will be multiplied by (1e12/N), as if there were 1e12 particles

constexpr int n_space_divx = n_space;
constexpr int n_space_divy = n_space;
constexpr int n_space_divz = n_space;
constexpr int n_space_divx2 = n_space_divx * 2;
constexpr int n_space_divy2 = n_space_divy * 2;
constexpr int n_space_divz2 = n_space_divz * 2;
constexpr int n_cells = n_space_divx * n_space_divy * n_space_divz;
constexpr int n_cells8 = n_cells * 8;
// physical "constants"
constexpr float kb = 1.38064852e-23;       // m^2kss^-2K-1
constexpr float e_charge = 1.60217662e-19; // C
constexpr float ev_to_j = e_charge;
constexpr float e_mass = 9.10938356e-31;
constexpr float e_charge_mass = e_charge / e_mass;
constexpr float kc = 8.9875517923e9;         // kg m3 s-2 C-2
constexpr float epsilon0 = 8.8541878128e-12; // F m-1
constexpr float pi = 3.1415926536;
constexpr float u0 = 4e-7 * pi;

class Time {
        private:
                vector<chrono::_V2::system_clock::time_point> marks;
        public:
                void mark();
                float elapsed();
                float replace();

};
class Log {
        private:
                ofstream log_file;
                bool firstEntry = true; // Whether the next item to print is the first item in the line
        public:
                Log();
                template<class T> void write(T text, bool flush = false){
                        if(!firstEntry) log_file << ",";
                        firstEntry = false;
                        log_file << text;
                        if(flush) log_file.flush();
                }
                void newline();
                void close();
};
static Time timer;
static Log logger;
void save_vti_c2(string filename, int i,
                 unsigned int n_space_div[3], float posl[3], float dd[3], uint64_t num, int ncomponents, double t,
                 float data1[3][n_space_divz2][n_space_divy2][n_space_divz2], string typeofdata, int bytesperdata);
void save_vti_c(string filename, int i,
                unsigned int n_space_div[3], float posl[3], float dd[3], uint64_t num, int ncomponents, double t,
                float data1[3][n_space_divz][n_space_divy][n_space_divz], string typeofdata, int bytesperdata);
void save_vti(string filename, int i, unsigned int n_space_div[3], float posl[3], float dd[3], uint64_t num, int ncomponents, double t, const char *data, string typeofdata, int sizeofdata);
void save_pvd(string filename, int ndatapoints);
void save_vtp(string filename, int i, uint64_t num, int ncomponents, double t, const char *data, const char *points);
void set_initial_pos_vel(int n_part_types, int n_particles, float *pos0, float *pos1, float *sigma, int *q, int *m, int *nt);
void cl_start();
void cl_set_build_options(float posL[3], float posH[3], float dd[3]);

void tnp(float *Ea1, float *Ba1, float *pos0x, float *pos0y, float *pos0z, float *pos1x, float *pos1y, float *pos1z,
         float cf[2],
         unsigned int ci[2]
        );
//void get_precalc_r3(float precalc_r3[3][n_space_divz2][n_space_divy2][n_space_divx2], float dd[3]);
int calcEBV(float V[n_space_divz][n_space_divy][n_space_divx],
           float E[3][n_space_divz][n_space_divy][n_space_divx], float B[3][n_space_divz][n_space_divy][n_space_divx],
           float Ee[3][n_space_divz][n_space_divy][n_space_divx], float Be[3][n_space_divz][n_space_divy][n_space_divx],
           float npt[n_space_divz][n_space_divy][n_space_divx], float jc[3][n_space_divz][n_space_divy][n_space_divx],
           float dd[3], float Emax, float Bmax);

void save_files(int i_time, unsigned int n_space_div[3], float posL[3], float dd[3], double t,
                float np[2][n_space_divz][n_space_divy][n_space_divx], float currentj[2][3][n_space_divz][n_space_divy][n_space_divx],
                float V[n_space_divz][n_space_divy][n_space_divx],
                float E[3][n_space_divz][n_space_divy][n_space_divx], float B[3][n_space_divz][n_space_divy][n_space_divx],
                float KE[2][n_output_part], float posp[2][n_output_part][3]);
void sel_part_print(int n_part[3],
                    float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd],
                    float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd],
                    float posp[2][n_output_part][3], float KE[2][n_output_part],
                    int m[2][n_partd], float dt[2]);

void get_densityfields(float currentj[2][3][n_space_divz][n_space_divy][n_space_divx],
                       float np[2][n_space_divz][n_space_divy][n_space_divx],
                       float npt[n_space_divz][n_space_divy][n_space_divx],
                       int nt[2], float KEtot[2], float posL[3], float posH[3], float dd[3],
                       float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd],
                       float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd],
                       int q[2][n_partd], float dt[2], int mp[2], int n_part[2],
                       float jc[3][n_space_divz][n_space_divy][n_space_divz]);
void calc_trilin_constants(float E[3][n_space_divz][n_space_divy][n_space_divx],
                           float Ea[n_space_divz][n_space_divy][n_space_divx][3][ncoeff],
                           float dd[3], float posL[3]);

void changedt(float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd], float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd], int n_part[3]);

void calcU(float V[n_space_divz][n_space_divy][n_space_divx],
                float E[3][n_space_divz][n_space_divy][n_space_divx], float B[3][n_space_divz][n_space_divy][n_space_divx],
                float posx[2][n_partd], float posy[2][n_partd], float posz[2][n_partd],
                float posL[3], float dd[3], int n_part[2], int q[2][n_partd], float out[2]);

void generateParticles(float a0, float r0, int *qs, int *mp, float pos0x[2][n_partd], float pos0y[2][n_partd], float pos0z[2][n_partd],
                        float pos1x[2][n_partd], float pos1y[2][n_partd], float pos1z[2][n_partd], int q[2][n_partd], int m[2][n_partd], int *nt);
void generateField(float Ee[3][n_space_divz][n_space_divy][n_space_divx], float Be[3][n_space_divz][n_space_divy][n_space_divx]);
void id_to_cell(int id, int *x, int *y, int *z);
#endif // TRAJ_H_INCLUDED
