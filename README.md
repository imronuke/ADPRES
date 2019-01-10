# ADPRES

Abu Dhabi Polytechnic Reactor Simulator (ADPRES) solves static and transient diffusion equation for two or three dimensional rectangular geometry. The method used is fourth order Nodal Expansion Method (NEM) where the transverse leakage moments are approximated using quadratic leakage fit. The transient diffusion equation is solved with fully-implicit method. ADPRES can handle problems with Assembly Discontinuity Factors (ADF) to imporove the accuracy. Recently, TH feedbacks module is added to enable ADPRES solving transient problems with TH feedbacks.

ADPRES is a great learning tool in the reactor theory class. The input is designed to straightforward. ADPRES' main objective is to make all nuclear engineering students have access on similar nuclear computer code. It is open, so everyone has access to the source code and modify for his/her own purposes.

## How to compile

You need a fortran compiler to compile the ADPRES code. Any fortran fortran compiler should work (I use gfortran).
First, you need to clone or download this repository to your computer, then to compile using gfortran (in an Unix based OS), type this command in the [src](https://github.com/imronuke/ADPRES/tree/master/src) folder

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

## How to develop my own input?

ADPRES input is designed to self-explanatory. It has 12 input cards, for example: `%mode`, `%geom`, `%xsec`, and so on. Some cards are obligatory for any problems, while some cards are conditional depend on the problem being solved and some cards are optional.

Followings are some step by step tips to develop your own input
* Open and try to understand the easiest input sample, for example [IAEA 3D PWR input](https://github.com/imronuke/ADPRES/blob/master/smpl/IAEA3D). From here you should have some ideas about the input structure.
* Modify the input sample for your own purposes
* In order to get detailed explanation, you may refer to the [ADPRES User Manual](https://github.com/imronuke/ADPRES/blob/reactivity/ADPRES%20USER%20MANUAL.pdf)

## How to read the ouput?

At first don't get confuse. The ouput is divided into three parts
1. Input echo. In this part, you will see your inputs along with the line numbers. This part makes sure that results produced were given from this input. In this [output](https://github.com/imronuke/ADPRES/blob/master/smpl/IAEA3D.out), it is from line 7 to 82
2. Input read. Now ADPRES is reading your input and restate your input to make sure these are the input paramaters that you want.  In this [output](https://github.com/imronuke/ADPRES/blob/master/smpl/IAEA3D.out), it is from line 86 to 236
3. Results. After doing some calculations, ADPRES presents the results in this part. For this [output](https://github.com/imronuke/ADPRES/blob/master/smpl/IAEA3D.out) for example, you can see the evolution of the outer iteration, final multiplication factor, radial power distribution as well as radial flux distribution.

## How accurate is ADPRES

ADPRES has been tested for both static and transient reactor problems
* [IAEA 3D PWR](https://engineering.purdue.edu/PARCS/Code/TestSuite/CalculationMode/StandAloneMode/Eigenvalue/IAEA3DPWR) benchmark. A static PWR benchmark. Here are the [input](https://github.com/imronuke/ADPRES/blob/master/smpl/IAEA3D) and [output](https://github.com/imronuke/ADPRES/blob/master/smpl/IAEA3D.out)
* [DVP BWR](http://li.mit.edu/Stuff/CNSE/Paper/Smith86PNE.pdf) benchmark. A static BWR bechmark with discontinuity factors. Here are the [input](https://github.com/imronuke/ADPRES/blob/master/smpl/DVP) and [output](https://github.com/imronuke/ADPRES/blob/master/smpl/DVP.out)
* [LMW](https://www.sciencedirect.com/science/article/pii/014919709500082U) benchmark. Reactor transient benchmark without TH feedbacks. Here are the [input](https://github.com/imronuke/ADPRES/blob/master/smpl/LMW) and [output](https://github.com/imronuke/ADPRES/blob/master/smpl/LMW.out)
* [NEACRP 3D PWR Core transient](https://www.oecd-nea.org/science/docs/1991/neacrp-l-1991-335.pdf) bechmark.  PWR transient benchmark with TH feedbacks. Inputs and the coressponding outputs are located in [this folder](https://github.com/imronuke/ADPRES/tree/master/smpl/NEACRP_TRANS)

## What to do next
Anyone wants to provide contribution is very welcome. Credits are given no matter how small his/her contribution. Here are the jobs list for future developments
- [ ] Develop subroutine for pin power reconstrurution
- [ ] Develop ADPRES for hexagonal geometry reactor problems (or might be used further for Molten Salt Reactor analysis with new transient module)
- [ ] Develop subroutine for decay heat calculation
- [ ] Develop subroutine for Xenon and Samarium treatment during transient
- [ ] Develop inteface with other TH codes such as RELAP or COBRA for coupled calculations
- [ ] Develop ADPRES for multiphysics calculation


## How to give feedbacks
Contact me
* muhammad.imron[at]adpoly.ac.ae
* makrus.imron[at]gmail.com

## How to cite

[will be updated]



> **"The best of people are those who bring most benefit to the rest of mankind." (THE PROPHET)**

