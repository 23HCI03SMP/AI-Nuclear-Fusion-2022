# AI Nuclear Fusion 2022
## Table of Contents
1. [Introduction](#introduction)
2. [Getting Started](#getting-started)
3. [Usage](#usage)

## Introduction
This project simulates Magnetized Liner Inertial Fusion (MagLIF). [Specifically, it simulates charged particles and their trajectories in β and e fields.]?

## Getting Started
### Prerequisites
- MSYS2
- GCC added to PATH in MSYS
    - In the root directory, run `export PATH=$PATH:/mingw64/bin`
- Required Libraries
    - [GSL](https://packages.msys2.org/package/mingw-w64-x86_64-gsl)
    - [OpenMP](https://packages.msys2.org/package/mingw-w64-x86_64-openmp)
    - [OpenCL-ICD](https://packages.msys2.org/package/mingw-w64-x86_64-opencl-icd)
    - [OpenCL-CLHPP](https://packages.msys2.org/package/mingw-w64-x86_64-opencl-clhpp)
## Usage
- Ensure that a `obj` folder is in the root folder
- Include the export directory in `traj.h`
- `make` the project by running `mingw32-make`