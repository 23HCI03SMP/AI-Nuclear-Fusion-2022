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
    - ```math
        F = q(E + v \times B) + r \times v \times w
        ```
    - Here, `r` is the artificial viscosity coefficient which is negative when energy is lost and positive when energy is gained
- [ ] Add temperature field Te[x][y][z]
  - Approximate Te as the average kinetic energy of particles
