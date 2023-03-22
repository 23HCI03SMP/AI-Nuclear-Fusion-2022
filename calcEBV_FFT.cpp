#include "include/traj.h"
#include <math.h>
#include <complex>
#include <fftw3.h>
void vector_muls(float *A, float Bb, int n)
{
    // Create a command queue
    cl::CommandQueue queue(context_g, default_device_g);
    float B[1] = {Bb};
    //  cout << B[0] << endl;
    // Create memory buffers on the device for each vector
    cl::Buffer buffer_A(context_g, CL_MEM_READ_WRITE, sizeof(float) * n);
    cl::Buffer buffer_B(context_g, CL_MEM_READ_ONLY, sizeof(float));

    // Copy the lists C and B to their respective memory buffers
    queue.enqueueWriteBuffer(buffer_A, CL_TRUE, 0, sizeof(float) * n, A);
    queue.enqueueWriteBuffer(buffer_B, CL_TRUE, 0, sizeof(float), B);

    // Create the OpenCL kernel
    cl::Kernel kernel_add = cl::Kernel(program_g, "vector_muls"); // select the kernel program to run

    // Set the arguments of the kernel
    kernel_add.setArg(0, buffer_A); // the 1st argument to the kernel program
    kernel_add.setArg(1, buffer_B);
    queue.enqueueNDRangeKernel(kernel_add, cl::NullRange, cl::NDRange(n), cl::NullRange);
    queue.finish(); // wait for the end of the kernel program
    // read result arrays from the device to main memory
    queue.enqueueReadBuffer(buffer_A, CL_TRUE, 0, sizeof(float) * n, A);
}

// Vector multiplication for complex numbers. Note that this is not in-place.
void vector_muls(fftwf_complex *dst, fftwf_complex *A, fftwf_complex *B, int n)
{
    // Create a command queue
    cl::CommandQueue queue(context_g, default_device_g);
    // Create memory buffers on the device for each vector
    cl::Buffer buffer_A(context_g, CL_MEM_WRITE_ONLY, sizeof(fftwf_complex) * n);
    cl::Buffer buffer_B(context_g, CL_MEM_READ_ONLY, sizeof(fftwf_complex) * n);
    cl::Buffer buffer_C(context_g, CL_MEM_READ_ONLY, sizeof(fftwf_complex) * n);

    // Copy the lists C and B to their respective memory buffers
    queue.enqueueWriteBuffer(buffer_B, CL_TRUE, 0, sizeof(fftwf_complex) * n, A);
    queue.enqueueWriteBuffer(buffer_C, CL_TRUE, 0, sizeof(fftwf_complex) * n, B);

    // Create the OpenCL kernel
    cl::Kernel kernel_add = cl::Kernel(program_g, "vector_mul_complex"); // select the kernel program to run

    // Set the arguments of the kernel
    kernel_add.setArg(0, buffer_A); // the 1st argument to the kernel program
    kernel_add.setArg(1, buffer_B);
    kernel_add.setArg(2, buffer_C);

    queue.enqueueNDRangeKernel(kernel_add, cl::NullRange, cl::NDRange(n), cl::NullRange);
    queue.finish(); // wait for the end of the kernel program
    // read result arrays from the device to main memory
    queue.enqueueReadBuffer(buffer_A, CL_TRUE, 0, sizeof(fftwf_complex) * n, dst);
}

bool checkInRange(string name, float data[3][n_space_divz][n_space_divy][n_space_divz], float minval, float maxval){
    bool toolow = false, toohigh = false;
    const float *data_1d = reinterpret_cast<float *>(data);
    for (unsigned int i = 0; i < n_cells * 3; ++i){
        toolow |= data_1d[i] < minval;
        toohigh |= data_1d[i] > maxval;
    }
    if (toolow)
    {
        const float *minelement = min_element(data_1d, data_1d + 3 * n_cells);
        size_t pos = minelement - &data[0][0][0][0];
        int count = 0;
        for (unsigned int n = 0; n < n_cells * 3; ++n) count += data_1d[n] < minval;
        int x, y, z;
        id_to_cell(pos, &x, &y, &z);
        cout << "Min " << name << ": " << *minelement<< " (" << x << "," << y << "," << z << ") (" << count << " values below threshold)\n";
    }
    if (toohigh) 
    {
        const float *maxelement = max_element(data_1d, data_1d + 3 * n_cells);
        size_t pos = maxelement - &data[0][0][0][0];
        int count = 0;
        for (unsigned int n = 0; n < n_cells * 3; ++n) count += data_1d[n] > maxval;
        int x, y, z;
        id_to_cell(pos, &x, &y, &z);
        cout << "Max " << name << ": " << *maxelement<< " (" << x << "," << y << "," << z << ") (" << count << " values above threshold)\n";
    }
    return toolow || toohigh;
}

// Shorthand for cleaner code
const size_t N0 = n_space_divx2, N1 = n_space_divy2, N2 = n_space_divz2,
        N0N1 = N0 * N1, N0N1_2 = N0N1 / 2,
        N2_c = N2 / 2 + 1; // Dimension to store the complex data, as required by fftw (from their docs)
const size_t n_cells4 = N0 * N1 * N2_c; // NOTE: This is not actually n_cells * 4, there is an additional buffer that fftw requires.

/*
The flow of this function is as follows, note that the FFT we are using converts real to complex, and the IFFT converts complex to real.
This is because we are making use of the fact that our original data is fully reall, and this saves about half the memory/computation
    Precalculation step on first run:
    precalc_r3_base -> FFT -> precalc_r3

    npt (number of particles per cell) -> copy -> FFT(planforE) -> fft_complex
    fft_complex *= precalc_r3, do this three times for three dimensions
    fft_complex -> IFFT(planbacE) -> fft_real
    E (input) = fft_real + Ee (background constant E also input)

    jc -> copy -> FFT(planforB) -> fft_complex
    fft_complex *= precalc_r3
    fft_complex -> IFFT(planbacB) -> fft_real
    B (input) = fft_real + Be (background constant E also input)

    Bonus: V potential. This is ran in parallel with E.
    npt -> copy -> FFT(planforE) -> fft_complex (reuse)
    fft_complex *= precalc_r2
    fft_complex -> IFFT(planbacV) -> fft_real
    V (input) = fft_real
*/

// Arrays for fft, output is multiplied by 2 because the convolution pattern should be double the size
// (a particle can be both 64 cells up, or 64 cells down, so we need 128 cells to calculate this information)
auto *fft_real = reinterpret_cast<float(&)[4][n_cells8]>(*fftwf_alloc_real(n_cells8 * 4));
auto *fft_complex = reinterpret_cast<fftwf_complex(&)[4][n_cells4]>(*fftwf_alloc_complex(4 * n_cells4));
//  pre-calculate 1/ r3 to make it faster to calculate electric and magnetic fields
auto *precalc_r3 = reinterpret_cast<fftwf_complex(&)[2][3][N2_c][N1][N0]>(*fftwf_alloc_complex(2 * 3 * n_cells4));

#ifdef Uon_ // similar arrays for U, but kept separately in one ifdef
auto *precalc_r2 = reinterpret_cast<fftwf_complex(&)      [N2_c][N1][N0]>(*fftwf_alloc_complex(n_cells4));
#endif

int calcEBV(float V[n_space_divz][n_space_divy][n_space_divx],
           float E[3][n_space_divz][n_space_divy][n_space_divx], float B[3][n_space_divz][n_space_divy][n_space_divx],
           float Ee[3][n_space_divz][n_space_divy][n_space_divx], float Be[3][n_space_divz][n_space_divy][n_space_divx],
           float npt[n_space_divz][n_space_divy][n_space_divx], float jc[3][n_space_divz][n_space_divy][n_space_divx], float dd[3],
           float Emax, float Bmax)
// float precalc_r3[3][n_space_divz2][n_space_divy2][n_space_divx2],  float Aconst, float Vconst,
{
    static int first = 1;
    static fftwf_plan planforE, planforB, planbacE, planbacB;

    static float posL2[3];
    static unsigned int *n_space_div2;
    if (first)
    { // allocate and initialize to 0
        auto precalc_r2_base = new float      [N2][N1][N0];
        auto precalc_r3_base = new float[2][3][N2][N1][N0];

        int dims[3] = {N0, N1, N2};

        // Create fftw plans
        fftwf_plan planfor_k = fftwf_plan_many_dft_r2c(3, dims, 6, reinterpret_cast<float *>(precalc_r3_base[0][0]), NULL, 1, n_cells8, reinterpret_cast<fftwf_complex *>(precalc_r3[0][0]), NULL, 1, n_cells4, FFTW_ESTIMATE);
        planforE = fftwf_plan_dft_r2c_3d(N0, N1, N2, fft_real[0], fft_complex[3], FFTW_MEASURE);
        #ifndef Uon_
        // Perform ifft on the first 3/4 of the array for 3 components of E field
        planbacE = fftwf_plan_many_dft_c2r(3, dims, 3, fft_complex[0], NULL, 1, n_cells4, fft_real[0], NULL, 1, n_cells8, FFTW_MEASURE);
        #else
        // Perform ifft on the entire array; the first 3/4 is used for E while the last 1/4 is used for V
        planbacE = fftwf_plan_many_dft_c2r(3, dims, 4, fft_complex[0], NULL, 1, n_cells4, fft_real[0], NULL, 1, n_cells8, FFTW_MEASURE);
        #endif
        planforB = fftwf_plan_many_dft_r2c(3, dims, 3, fft_real[0], NULL, 1, n_cells8, fft_complex[0], NULL, 1, n_cells4, FFTW_MEASURE);
        planbacB = fftwf_plan_many_dft_c2r(3, dims, 3, fft_complex[0], NULL, 1, n_cells4, fft_real[0], NULL, 1, n_cells8, FFTW_MEASURE);
        
        #ifdef Uon_
        fftwf_plan planfor_k2 = fftwf_plan_dft_r2c_3d(N0, N1, N2, reinterpret_cast<float *>(precalc_r2_base), reinterpret_cast<fftwf_complex *>(precalc_r2), FFTW_ESTIMATE);
        #endif
        cout << "allocate done\n";
        float r3, rx, ry, rz, rx2, ry2, rz2;
        int i, j, k, loc_i, loc_j, loc_k;
        posL2[0] = -dd[0] * ((float)n_space_divx - 0.5);
        posL2[1] = -dd[1] * ((float)n_space_divy - 0.5);
        posL2[2] = -dd[2] * ((float)n_space_divz - 0.5);
        n_space_div2 = new unsigned int[3]{n_space_divx2, n_space_divy2, n_space_divz2};
        // precalculate 1/r^3 (field) and 1/r^2 (energy)
        for (k = -n_space_divz; k < n_space_divz; k++)
        {
            loc_k = k + (k < 0 ? n_space_divz2 : 0); // The "logical" array position
            // We wrap around values smaller than 0 to the other side of the array, since 0, 0, 0 is defined as the center of the convolution pattern an hence rz should be 0
            rz = k * dd[2]; // The change in z coordinate for the k-th cell.
            rz2 = rz * rz;
            for (j = -n_space_divy; j < n_space_divy; j++)
            {
                loc_j = j + (j < 0 ? n_space_divy2 : 0);
                ry = j * dd[1];
                ry2 = ry * ry + rz2;
                for (i = -n_space_divx; i < n_space_divx; i++)
                {
                    loc_i = i + (i < 0 ? n_space_divx2 : 0);
                    rx = i * dd[0];
                    rx2 = rx * rx + ry2;
                    r3 = rx2 == 0? 0.f : powf(rx2, -1.5) / n_cells8;
                    precalc_r3_base[0][0][loc_k][loc_j][loc_i] = r3 * rx;
                    precalc_r3_base[0][1][loc_k][loc_j][loc_i] = r3 * ry;
                    precalc_r3_base[0][2][loc_k][loc_j][loc_i] = r3 * rz;
                    precalc_r3_base[1][0][loc_k][loc_j][loc_i] = precalc_r3_base[0][0][loc_k][loc_j][loc_i];
                    precalc_r3_base[1][1][loc_k][loc_j][loc_i] = precalc_r3_base[0][1][loc_k][loc_j][loc_i];
                    precalc_r3_base[1][2][loc_k][loc_j][loc_i] = precalc_r3_base[0][2][loc_k][loc_j][loc_i];
                    #ifdef Uon_
                    precalc_r2_base[loc_k][loc_j][loc_i] = rx2 == 0? 0.f : powf(rx2, -0.5) / n_cells8;
                    #endif
                }
            }
        }
        // Multiply by the respective constants here, since it is faster to parallelize it
        const float Vconst = kc * e_charge * r_part_spart;
        const float Aconst = 1e-7 * e_charge * r_part_spart;

        vector_muls(reinterpret_cast<float *>(precalc_r3_base[0]), Vconst, n_cells8 * 3);
        vector_muls(reinterpret_cast<float *>(precalc_r3_base[1]), Aconst, n_cells8 * 3);
        #ifdef Uon_
        vector_muls(reinterpret_cast<float *>(precalc_r2_base), Vconst, n_cells8);
        fftwf_execute(planfor_k2);
        fftwf_destroy_plan(planfor_k2);
        #endif

        fftwf_execute(planfor_k); // fft of kernel arr3=fft(arr)
        fftwf_destroy_plan(planfor_k);
        delete[] precalc_r3_base;
        delete[] precalc_r2_base;
    }
    #ifdef Eon_
    //#pragma omp parallel sections
    {
        fill(&fft_real[0][0], &fft_real[3][n_cells8], 0.f);
        //#pragma omp section
        {
            // density field
            size_t i, j, k, jj;
            jj = 0;
            for (k = 0; k < n_space_divz; ++k)
            {
                for (j = 0; j < n_space_divy; ++j)
                {
                    for (i = 0; i < n_space_divx; ++i){
                        fft_real[0][jj + i] = npt[k][j][i];
                    }
                    jj += N0;
                }
                jj += N0N1_2;
            }
            //          cout << "density\n";
            fftwf_execute(planforE); // arrn1 = fft(arrn) multiply fft charge with fft of kernel(i.e field associated with 1 charge)
            for (int c = 0; c < 3; c++){
                //vector_muls(fft_complex[c], fft_complex[3], reinterpret_cast<fftwf_complex *>(precalc_r3[0][c]), n_cells4);
                // ^ read/write is slow, maybe use host pointer iff it is supported
                const auto dst_std = reinterpret_cast<complex<float> *>(fft_complex[c]), src_std = reinterpret_cast<complex<float> *>(fft_complex[3]),
                precalc_r3_std = reinterpret_cast<complex<float> *>(precalc_r3[0][c]);
                for (int i = 0; i < n_cells4; ++i)
                    dst_std[i] = src_std[i] * precalc_r3_std[i];
            }
            #ifdef Uon_
            {
                const auto ptr_std = reinterpret_cast<complex<float> *>(fft_complex[3]),
                precalc_r2_std = reinterpret_cast<complex<float> *>(precalc_r2);
                complex<float> temp;
                for (int i = 0; i < n_cells4; ++i)
                    ptr_std[i] *= precalc_r2_std[i];
            }
            #endif
            fftwf_execute(planbacE); // inverse transform to get convolution
                                     //#pragma omp parallel for
            for (int c = 0; c < 3; c++)
            { // 3 axis
                const float* fft_real_c = fft_real[c];
                size_t i, j, k, jj;
                //               cout << "c " << c << ", thread " << omp_get_thread_num() << ", jj " << jj << endl;
                jj = 0;
                for (k = 0; k < n_space_divz; ++k)
                {
                    for (j = 0; j < n_space_divy; ++j)
                    {
                        for (i = 0; i < n_space_divx; ++i)
                            E[c][k][j][i] = fft_real_c[jj + i] + Ee[c][k][j][i];
                        jj += N0;
                    }
                    jj += N0N1_2;
                }
            }
            #ifdef Uon_
            jj = 0;
            for (k = 0; k < n_space_divz; ++k){
                for (j = 0; j < n_space_divy; ++j)
                    memcpy(V[k][j], &fft_real[3][jj += N0], sizeof(float) * n_space_divx);
                jj += N0N1_2;
            }
            #endif
        }
    }
    #else
    memcpy(reinterpret_cast<float *>(E), reinterpret_cast<float *>(Ee), 3 * n_cells * sizeof(float));
    #endif

    #ifdef Bon_
    {
        fill(&fft_real[0][0], &fft_real[2][n_cells8], 0.f);
        //#pragma omp section
        //#pragma omp parallel for
        for (int c = 0; c < 3; c++)
        { // 3 axis
            size_t i, j, k, jj;
            jj = 0;
            for (k = 0; k < n_space_divz; ++k)
            {
                for (j = 0; j < n_space_divy; ++j)
                    memcpy(&fft_real[c][jj += N0], jc[c][k][j], sizeof(float) * n_space_divx);
                jj += N0N1_2;
            }
        }
        fftwf_execute(planforB); // arrj1 = fft(arrj)

        const auto ptr = reinterpret_cast<complex<float> *>(fft_complex), r3_1d = reinterpret_cast<complex<float> *>(precalc_r3[1]);
        complex<float> temp1, temp2, temp3;
        for (int i = 0, j = n_cells4, k = n_cells4 * 2; i < n_cells4; ++i, ++j, ++k)
        {
            temp1 = ptr[j] * r3_1d[k] - ptr[k] * r3_1d[j];
            temp2 = ptr[k] * r3_1d[i] - ptr[i] * r3_1d[k];
            temp3 = ptr[i] * r3_1d[j] - ptr[j] * r3_1d[i];
            ptr[i] = temp1;
            ptr[j] = temp2;
            ptr[k] = temp3;
        }
        fftwf_execute(planbacB);
        for (int c = 0; c < 3; c++)
        { // 3 axis
            const float* fft_real_c = fft_real[c];
            size_t i, j, k, jj;
            jj = 0;
            for (k = 0; k < n_space_divz; ++k)
            {
                for (j = 0; j < n_space_divy; ++j)
                {
                    for (i = 0; i < n_space_divx; ++i)
                        B[c][k][j][i] = fft_real_c[jj + i] + Be[c][k][j][i];
                    jj += N0;
                }
                jj += N0N1_2;
            }
        }
    }
    #else
    memcpy(reinterpret_cast<float *>(B), reinterpret_cast<float *>(Be), 3 * n_cells * sizeof(float));
    #endif

    #ifdef Uon_
    #ifdef Eon_ // if both Uon and Eon are defined
    {
        // Perform estimate of electric potential energy
        size_t i, j, k, jj = 0;
        float EUtot = 0.f;
        const float* V_1d = reinterpret_cast<float *>(V);
        const float* npt_1d = reinterpret_cast<float *>(npt);
        for (int i = 0; i < n_cells; ++i){
            EUtot += V_1d[i] * npt_1d[i];
        }
        EUtot *= 0.5f;// * e_charge / ev_to_j; <- this is just 1
        cout << "Eele (estimate): " << EUtot << ", ";
    }
    #endif
    #endif
    first = 0;
    bool E_exceeds = checkInRange("E", E, -Emax, Emax), B_exceeds = checkInRange("B", B, -Bmax, Bmax);
    return E_exceeds || B_exceeds;
}