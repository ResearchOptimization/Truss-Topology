# Truss-Topology
In this repository you will find the code to solve truss topology optimization using Bernstein approximation


This directory contains the subroutines used to optimize the truss structures corresponding
to the article "Topology optimization of truss structures under failure probability using the
Bernstein approximation" submitted to Computer and Structures by Alfredo Canelas, Miguel Carrasco, and Julio López.

# The files are:

1. Elaslin.jl    : Julia module with the subroutines to perform the optimization.
2. figuretruss.m : Subroutine for drawing 2D or 3D trusses in Matlab.
3. P2DCrane1.jl  : Subroutine corresponding to the 2D Crane example with a mesh of 13 x 4 nodes.
4. P2DCrane2.jl  : Subroutine corresponding to the 2D Crane example with a mesh of 17 x 5 nodes.
5. P2DCrane3.jl  : Subroutine corresponding to the 2D Crane example with a mesh of 25 x 7 nodes.
6. P3DBeam1.jl   : Subroutine corresponding to the 3D Beam example with a mesh of 3 x 3 x 9 nodes.
7. P3DBeam2.jl   : Subroutine corresponding to the 3D Beam example with a mesh of 5 x 5 x 17 nodes.

# Instalation instructions:

1. Install Julia in your system, see https://julialang.org/
2. Install Mosek in your system, see https://www.mosek.com/
3. Install MATLAB in your system, see https://mathworks.com/

After instalation check if Julia contains the following required packages:
Mosek, SparseArrays, LinearAlgebra, Printf, FinancialToolbox

# Execution instructions:

1. Run the desired example, for example in a Linux terminal:

> julia P2DCrane1.jl

The output for this example is:

Variables:   629

Iterations:  32

Time:        0.3264818589996139


In addition, the provided subroutines write the MATLAB file 'P2DCrane1.m' with the solution.

2. Run the file 'P2DCrane1.m' in MATLAB:

> P2DCrane1

The MATLAB environment will be added with the following variables:

Mnodo:   Matrix of nodal coordinates.

Melem:   Matrix of element conectivities.

PropGeo: Vector of bar volumes (normalized to a total volume = 1).

vtol:    Volume tolerance for drawings.


3. View the optimized structure:

> figuretruss
