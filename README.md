# ADPRES

Abu Dhabi Polytechnic Reactor Simulator (ADPRES) solves static and transient diffusion equation for two or three dimensional rectangular geometry. The method used is fourth order Nodal Expansion Method (NEM) where the transverse leakage moments are approximated using quadratic leakage fit. ADPRES can solve either forward or adjoint problems as well as fixed source problems. It also can handle problems with Assembly Discontinuity Factors (ADF) to imporove the accuracy.

ADPRES is a great learning tool in the reactor theory class. The input is designed to straightforward. ADPRES' main objective is to make all nuclear engineering students have access on similar nuclear computer code. It is open, so everyone has access to the source code and modify for his/her own purposes.

## How to compile

You need a fortran compiler to compile the ADPRES code. Any fortran fortran compiler should work (I use gfortran).
To compile using gfortran in an Unix based OS, type this command in the [src](https://github.com/imronuke/ADPRES/tree/master/src) folder

```
gfortran mod_data.f90 mod_io.f90 mod_nodal.f90 mod_th.f90 mod_trans.f90 ADPRES.f90 -O4 -o adpres
```

This command shall produce executable file 'adpres'

## How to run a test

Copy some sample inputs file from folder [smpl](https://github.com/imronuke/ADPRES/tree/master/smpl) into the [src](https://github.com/imronuke/ADPRES/tree/master/src) folder, then execute adpres

```
./adpres
```

then type the input file name. It will produce output file with out extension

## How accurate is ADPRES

ADPRES has been tested for both static and transient reactor problems
* [IAEA 3D PWR](https://engineering.purdue.edu/PARCS/Code/TestSuite/CalculationMode/StandAloneMode/Eigenvalue/IAEA3DPWR) benchmark. A static PWR benchmark. Here are the [input]((https://github.com/imronuke/ADPRES/blob/master/smpl/IAEA3D)) and [output](https://github.com/imronuke/ADPRES/blob/master/smpl/IAEA3D.out)
* [DVP BWR](https://github.com/imronuke/ADPRES/blob/master/smpl/DVP) benchmark. A static BWR bechmark with discontinuity factors
* [LMW](https://github.com/imronuke/ADPRES/blob/master/smpl/LMW) benchmark. Reactor transient benchmark without TH feedbacks.
* [NEACRP 3D PWR Core transient](https://github.com/imronuke/ADPRES/tree/master/smpl/NEACRP_TRANS) bechmark. PWR transient benchmark with TH feedbacks.


#How to give feebacks
Contact me
* muhammad.imron[at]adpoly.ac.ae
* makrus.imron[at]gmail.com

## How to cite

[will be updated]



"The best of people are those that bring most benefit to the rest of mankind." (PROPHET)

