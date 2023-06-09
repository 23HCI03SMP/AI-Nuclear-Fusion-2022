# AI Nuclear Fusion 2022
## Table of Contents
- [AI Nuclear Fusion 2022](#ai-nuclear-fusion-2022)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
  - [Usage](#usage)
  - [TODO](#todo)

## Introduction
This project aims to develop a Particle-in-Cell code to simulate charged particles and their trajectories in β and e fields. [Original Code](https://github.com/23HCI03SMP/AI-Nuclear-Fusion-2021) was developed in 2021 by Hilary, Yin Yue and Chloe, with extensive improvments by Samuel, Ananth and Vishwa.

## Getting Started
### Prerequisites
- MSYS2
- GCC added to PATH in MSYS and Windows
    - In MSYS in the root directory, run `export PATH=$PATH:/mingw64/bin`
- [GDB](https://packages.msys2.org/package/mingw-w64-x86_64-gdb) for debugging
- Required Libraries
    - [GSL](https://packages.msys2.org/package/mingw-w64-x86_64-gsl)
    - [OpenMP](https://packages.msys2.org/package/mingw-w64-x86_64-openmp)
    - [OpenCL-ICD](https://packages.msys2.org/package/mingw-w64-x86_64-opencl-icd)
    - [OpenCL-CLHPP](https://packages.msys2.org/package/mingw-w64-x86_64-opencl-clhpp)
    - [FFTW3](https://packages.msys2.org/package/mingw-w64-x86_64-fftw)
- To get more performance, you may recompile libraries, such as recompiling FFTW3 with OMP enabled
  - ```bash
    > wget https://www.fftw.org/fftw-3.3.10.tar.gz
    > tar xvzf fftw-3.3.10.tar.gz
    > cd fftw-3.3.10
    > ./configure --enable-threads --enable-openmp --enable-avx --enable-avx2 --enable-avx512 --enable-avx-128-fma --enable-float --with-our-malloc --enable-sse2
    > make
    > make install
    ```

## Usage
- Ensure that a `obj` folder is in the root folder
- Include the export directory in `traj.h`
- `make` the project by running `mingw32-make`

## TODO
- [ ] Set up particles for hot rod project 
- [ ] Set up particles for hot plate project
- [ ] Set up particles for MagLIF
  - [ ] Two sets of plamas similar to hot rod project
    - [ ] Low density plasma cylinder
    - [ ] High density thin cylindrical shell around plasma cylinder
  - [ ] Set up external fields
    - [ ] Electric field around the cylinder
- [ ] Add in artificial viscosity to simulate energy loss/gain
    - $\vec{F} = m\frac{d \vec{v}}{dt} = q(\vec{E} - \vec{B} \times \vec{v}) + r \vec{v}$
    - Here, `r` is the artificial viscosity coefficient which is negative when energy is lost and positive when energy is gained
- [ ] Add temperature field Te[x][y][z]
  - Approximate Te as the average kinetic energy of particles
- [ ] viscosity field[p][x][y][z]
  - https://tanimislam.github.io/research/NRL_Formulary_2019.pdf
  - ?
    - $P = Fv = R \times v \times v, \therefore r = P_{perparticle}/v^2$
  - Bremsstrahlung from hydrogen-like plasma
    - $P_{Br} = 1.69 \times 10^{-32}N_eT_e^{1/2} \sum[Z^2N(Z)]$ $watt/cm^3$, where the sum of all ionization states Z
  - Cyclotron radiation loss from NRL Plasma formulary in magnetic field **B**
    - $P_c = 6.21 \times 10^{-38}B^2N_eT_e$ $watt/cm^3$
  - Power density released in the form of charged particles from fusion from NRL Plasma formulary (Assume to be absorbed within the cell for the time being not realistic as fusion products may have very long range. If cell sizes are small, most of energy from fusion products will be lost from the cell) all in in $watts/cm^3$
    - $P_{DD} = 3.3 \times 10^{-13}n_D{}^2(\overline{\sigma v})_{DD}$ $watt/cm^3$
    - $P_{DT} = 5.6 \times 10^{-13}n_Dn_T(\overline{\sigma v})_{DT}$ $watt/cm^3$
    - $P_{DHe^3} = 2.9 \times 10^{-12}n_Dn_{He^3}(\overline{\sigma v})_{DHe^3}$ $watt/cm^3$
  - [Power transfer between charged particles](https://tanimislam.github.io/research/NRL_Formulary_2019.pdf#page=32), where test particle energy $\epsilon$ and field particle temperature $T$ are both in eV; $\mu = m_i/m_p$ where $m_P$ is the proton mass; $Z$ is ion charge state; in electron–electron and ion–ion encounters, field particle quantities are distinguished by a prime.
  - The two expressions given below for each rate hold for very slow ($x^{\alpha \backslash \beta} \ll 1$) and very fast ($x^{\alpha \backslash \beta} \gg 1$) test particles, respectively.
    - Electron-electron
      - $v_s^{e|e}/n_e \lambda_{ee} \approx 5.8 \times 10^{-6}T^{-3/2} \longrightarrow 7.7 \times 10^{-6} \epsilon^{-3/2}$
      - $v_\perp^{e|e}/n_e \lambda_{ee} \approx 5.8 \times 10^{-6}T^{-1/2} \epsilon^{-1} \longrightarrow 7.7 \times 10^{-6} \epsilon^{-3/2}$
      - $v_\parallel^{e|e}/n_e \lambda_{ee} \approx 2.9 \times 10^{-6}T^{-1/2} \epsilon^{-1} \longrightarrow 3.9 \times 10^{-6} T \epsilon^{-5/2}$
    - Electron-ion
      - $v_s^{e|i}/n_iZ^2 \lambda_{ei} \approx 0.23\mu^{3/2}T^{-3/2} \longrightarrow 3.9 \times 10^{-6}\epsilon^{-3/2}$
      - $v_\perp^{e|i}/n_iZ^2\lambda_{ei} \approx 2.5 \times 10^{-4}\mu^{1/2}T^{-1/2}\epsilon^{-1} \longrightarrow 7.7 \times 10^{-6}\epsilon^{-3/2}$
      - $v_\parallel^{e|i}/n_iZ^2 \lambda_{ei} \approx 1.2 \times 10^{-4}\mu^{1/2}T^{-1/2}\epsilon^{-1} \longrightarrow 2.1 \times 10^{-9}\mu^{-1}T\epsilon^{-5/2}$
    - Ion-electron
      - $v_s^{i|e}/n_eZ^2 \lambda_{ie} \approx 1.6 \times 10^{-9}\mu^{-1}T^{-3/2} \longrightarrow 1.7 \times 10^{-4}\mu^{1/2}\epsilon^{-3/2}$
      - $v_\perp^{i|e}/n_eZ^2 \lambda_{ie} \approx 3.2 \times 10^{-9}\mu^{-1}T^{-1/2}\epsilon^{-1} \longrightarrow 1.8 \times 10^{-7}\mu^{-1/2}\epsilon^{-3/2}$
      - $v_\parallel^{i|e}/n_eZ^2 \lambda_{ie} \approx 1.6 \times 10^{-9}\mu^{-1}T^{-1/2}\epsilon^{-1} \longrightarrow 1.7 \times 10^{-4}\mu^{-1/2}\epsilon^{-5/2}$
    - Ion-ion
      - $\frac{v_s^{i|i^\prime}}{n_{i^\prime Z^2Z^{\prime 2}\lambda_{ii^\prime}}} \approx 6.8 \times 10^{-8}\frac{\mu^{\prime 1/2}}{\mu}(1 + \frac{\mu^\prime}{\mu})T^{-3/2} \longrightarrow 9.0 \times 10^{-8}(\frac{1}{\mu} + \frac{1}{\mu^\prime})\frac{\mu^{1/2}}{\epsilon^{3/2}}$
      - $\frac{v_\perp^{i|i^\prime}}{n_{i^\prime Z^2Z^{\prime 2}\lambda_{ii^\prime}}} \approx 1.4 \times 10^{-7}\mu^{\prime 1/2}\mu^{-1}T^{-1/2}\epsilon^{-1} \longrightarrow 1.8 \times 10^{-7}\mu^{-1/2}\epsilon^{-3/2}$
      - $\frac{v_\parallel^{i|i^\prime}}{n_{i^\prime Z^2Z^{\prime 2}\lambda_{ii^\prime}}} \approx 6.8 \times 10^{-8}\mu^{\prime 1/2}\mu^{-1}T^{-1/2}\epsilon^{-1} \longrightarrow 9.0 \times 10^{-8}\mu^{1/2}\mu^\prime\epsilon^{-5/2}$