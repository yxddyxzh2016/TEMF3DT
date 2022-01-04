# TEMF3DT

This is TEMF3DT, an algorithm that carries out both large-loop and small-loop TEM modeling.

Please cite the following paper if the software is used in your research:

A Fast 3-D finite element modeling algorithm for land transient electromagnetic method with OneAPI acceleration

TEMF3DT is copyright (C) 2021-2022 Xiaodong Yang.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

AUTHOR: Xiaodong Yang

EMAIL: yxddyxzh@mail.ustc.edu.cn

ADDRESS: School of Earth and Space Sciences, University of Science and Technology of China, Hefei 230026, China

AIM: An algorithm that carries out both large-loop and small-loop TEM modeling. The algorithm is based on the 3-D time-domain finite element method with unstructured edge-based meshes, and an unconditional stable time-stepping method is applied to acquire stable solutions. Several models are provided along with the program. The program is compatible with PC, workstation, and cluster usage with OpenMP and MPI parallelization.

# Version history

1.0     2021/03/21      initial release

2.0     2022/01/02      bug fix, open source

# 1 Environment

## 1.1 Windows operating system

Install Visual Studio 2022 with C++ development tool selected, together with Intel OneAPI basekit and hpckit libraries, refer to following website:

https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit/download.html?operatingsystem=window&distributions=webdownload

and

https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit/download.html?operatingsystem=window&distributions=webdownload

## 1.2 Linux operating system

Install Intel OneAPI basekit and hpckit packages, refer to following website:

https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit/download.html?operatingsystem=linux

and

https://software.intel.com/content/www/us/en/develop/tools/oneapi/hpc-toolkit/download.html?operatingsystem=linux

Set environmental variables BEFORE running example codes as:

source /opt/intel/oneapi/setvars.sh --force

Set stacksize as:

ulimit -s unlimited

# 2 Run example models

## 2.1 Large-loop 1-D model

On Windows: ./mpiexec -n 1 -genv OMP_NUM_THREADS=8 -genv I_MPI_PIN_DOMAIN=omp ./FETD -f mod_big_sig_0.001

On Linux: mpirun -n 1 -genv OMP_NUM_THREADS=8 -genv I_MPI_PIN_DOMAIN=omp ./FETD -f mod_big_sig_0.001

## 2.2 Small-loop 1-D model

On Windows: ./mpiexec -n 1 -genv OMP_NUM_THREADS=8 -genv I_MPI_PIN_DOMAIN=omp ./FETD -f mod_small_sig_0.001

On Linux: mpirun -n 1 -genv OMP_NUM_THREADS=8 -genv I_MPI_PIN_DOMAIN=omp ./FETD -f mod_small_sig_0.001

## 2.3 Large-loop 3-D model

On Windows: ./mpiexec -n 1 -genv OMP_NUM_THREADS=8 -genv I_MPI_PIN_DOMAIN=omp ./FETD -f mod3_big

On Linux: mpirun -n 1 -genv OMP_NUM_THREADS=8 -genv I_MPI_PIN_DOMAIN=omp ./FETD -f mod3_big

## 2.4 Small-loop 3-D model

Run the mod3_small.m MatLab script

# 3 Build your own model

## 3.1 Gmsh script

Refer to https://gmsh.info/ for Gmsh scripting info.

Here is a sample Gmsh file mod_big_sig_0.001.geo that could be used by FETD program, you can find it within the package:

// set the min and max mesh size

SetFactory("OpenCASCADE");

Mesh.CharacteristicLengthMin = 1;

Mesh.CharacteristicLengthMax = 20000;

// set the physical volumes with physical ID 11, 12, 13, ..., the physical domains must start with number 11, and no more than 50

Box(1)={-50000,-50000,-50000,100000,100000,50000};

Box(2)={-50000,-50000,     0,100000,100000,50000};

Physical Volume ( "air",11 ) = {1};

Physical Volume ( "layer",12 ) = {2};

// Coherence is needed right after the domain section to avoid duplicate mesh points

Coherence;

// set line sources, each source has its own physical ID 51, 52, 53, ..., the physical number of lines must start with number 51

l1=newreg; Circle(l1) = {0,0,30,500};

Curve{l1} In Volume {2};

Physical Curve ( "lsource1",51 ) = {l1};

// this section is optional, you can set near which receiver points the mesh is going to have a smaller size than the background

p71=newp;Point(p71)={0,0,3};

Coherence;

// field, needs to be changed according to the model

Field[1] = Distance;

Field[1].NNodesByEdge = 160;

Field[1].EdgesList =  {l1};

Field[2] = Threshold;

Field[2].IField = 1;

Field[2].LcMin = 20;

Field[2].LcMax =20000;

Field[2].DistMin = 0;

Field[2].DistMax = 100000;

Field[3] = Distance;

Field[3].NodesList = {p71};

Field[4] = Threshold;

Field[4].IField = 3;

Field[4].LcMin = 1;

Field[4].LcMax =20000;

Field[4].DistMin = 0;

Field[4].DistMax = 50000;

Field[11] = Min;

Field[11].FieldsList = {2,4};

Background Field = 11;

## 3.2 Meshing

The mesh should be saved as 2 seperate files:

1. The gmsh default .msh format file;
3. The gmsh INRIA MEDIT .mesh format file with "physical entities" selected.

## 3.3 .cntl file

This file should have the same prefix as the mesh files. Here is a sample mod_big_sig_0.001.cntl file that could be used by FETD program:

! set the receiver locations with x-, y- and z- coordinates

receivers: 1

0.0      0.0  3.0

! set the media chatacteristics, each record should contain the physical ID as the same as .geo file,

! followed by the resistivity, the relative permittivity and the relative magnetic conductivity.

isotropic: 2

11 1.d-12 1.0 1.0

12 0.001 1.0 1.0

! set source condiguration, the first line defines the number of source lines

! the following lines define the pulse type (1:Step pulse 2:Delta pulse), the pulse width of the source and the maximum current of each source

source: 1

1(should not be changed)    1(TYPE 1 - STEP / TYPE 2 - IMPULSE)   1.0e-7(t_0 as mentioned in the paper)    1.0

! number to divide pulse width, the initial time step is thus calculated by t_0/sp_division

sp_division: 30

! the maximum time of the computation

time_maximum: 0.1

! number of steps to double the time step size

sp_double: 32

! n_multi as mentioned in the paper

sp_dblesize: 2

! define time channels you need a time slice

vtkout: 0

## 3.4 Run and plot

On Windows: ./mpiexec -n 1 -genv OMP_NUM_THREADS=8 -genv I_MPI_PIN_DOMAIN=omp ./FETD -f prefix

On Linux: mpirun -n 1 -genv OMP_NUM_THREADS=8 -genv I_MPI_PIN_DOMAIN=omp ./FETD -f prefix

After running the FETD program, the required field componets at receiver locations will be ouput to .field file.

And if vtkout is required, the field components will be output to seperate .msh files, which could be opened by Gmsh for visualization.
