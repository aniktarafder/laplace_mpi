# laplace_mpi
The code was originally written by **John Urbanic, PSC**. I modified how the arrays are defined and added one subroutine to plot **Temperature** field in **Tecplot** format.

This code solves laplace heat equation in a 2D domain using jacobian iterations. 
* Left and bottom boundary --> 0 Temperature.
* Top boundary   --> Temperature varies from 0 to 100 linearly from left to right.
* Right boundary --> Temperature varies from 0 to 100 linearly from bottom to top.

### The code is hard-coded to run only with 4 processors.
Just clone the **repo** --> **compile** --> **run**

The code will generate 4 **.dat** files of temperature fields. As there are 4 processors involved.

* **To clone** : git clone https://github.com/aniktarafder/laplace_mpi.git
* **To compile** : mpif90 -o laplace_serial.exe laplace_serial.f90
* **To run**     : mpirun -n 4 ./laplace_serial.f90

**Course Link:** https://www.psc.edu/resources/training/xsede-hpc-workshop-may-4-5-2021-mpi/
