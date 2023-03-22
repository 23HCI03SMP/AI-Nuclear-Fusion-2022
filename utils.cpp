#include "include/traj.h"

void id_to_cell(int id, int*x, int*y, int*z){
    constexpr size_t xy = n_space_divx * n_space_divy;
    *z = id / xy;
    id = id % xy;
    *y = id / n_space_divx;
    *x = id % n_space_divx;
}

void Time::mark(){
    marks.push_back(chrono::high_resolution_clock::now());
}

float Time::elapsed(){
    unsigned long long time = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - marks.back()).count();
    marks.pop_back();
    return (float)time * 1e-6;
}

// Get the same result as elapsed, but also insert the current time point back in
float Time::replace(){
    auto now = chrono::high_resolution_clock::now();
    auto back = marks.back();
    unsigned long long time = chrono::duration_cast<chrono::microseconds>(now - back).count();
    back = now;
    return (float)time * 1e-6;
}

Log::Log() { if(!log_file.is_open()) log_file.open("log.csv"); log_file << setprecision(5); }

void Log::newline(){
    log_file << "\n";
    log_file.flush();
    firstEntry = true;
}
void Log::close(){
    log_file.close();
}