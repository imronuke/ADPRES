# ADPRES

Abu Dhabi Polytechnic Reactor Simulator (ADPRES) solves static and transient diffusion equation for two or three dimensional rectangular geometry. The method used is fourth order Nodal Expansion Method (NEM) where the transverse leakage moments are approximated using quadratic leakage fit. ADPRES can solve either forward or adjoint problems as well as fixed source problems. It also can handle problems with Assembly Discontinuity Factors (ADF) to imporove the accuracy.

ADPRES is a great learning tool in the reactor theory class. The input is designed to straightforward. ADPRES' main objective is to make all nuclear engineering students have access on similar nuclear computer code. It is open, so everyone has access to the source code and modify for his/her own purposes.

## How to compile

You need a fortran compiler to compile the ADPRES code. Any fortran fortran compiler should work (I use gnufortran).
To compile using gnufortran, type this command in the working folder




