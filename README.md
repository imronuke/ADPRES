[![Language](https://img.shields.io/badge/made%20with-fortran-blueviolet)](https://gcc.gnu.org/wiki/GFortran) [![Build Status](https://travis-ci.com/imronuke/ADPRES.svg?branch=master)](https://travis-ci.com/imronuke/ADPRES) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/imronuke/ADPRES/blob/master/LICENSE)
# ADPRES

Abu Dhabi Polytechnic Reactor Simulator (ADPRES) is an open nuclear reactor simulator and reactor core analysis tool that solves static and transient neutron diffusion equation for two or three dimensional reactor problems in Cartesian geometry. Currently, ADPRES uses Semi-Analytic Nodal Method (SANM) to spatially discretised the neutron diffusion equation. While theta method is used for the time discretisation.

ADPRES is also a great learning tool for reactor theory classes, and we have been striving hard to make the input is easy to create. Among ADPRES' main objectives is to make all nuclear engineering students have access on reactor simulator code for them to learn. It is open and completely free, so everyone has access to the source code and modify for his/her own purposes.

ADPRES features:
* Input is straightforward, modular and in a free-format form
* Solves both static and transient core problems **with or without TH feedback**
* Performs forward, adjoint and fixed-source calculations
* Performs calculations using branched cross sections data. An example of the library format can be seen [here](https://github.com/imronuke/ADPRES/blob/master/smpl/xsec/SERPENT_CMM/m40.tab)
* Critical boron concentration search
* Rod ejection simulation or Reactivity Initiated Accident (RIA)
* Solves multi-group of neutron energy
* Solves calculations with Assembly Discontinuity Factors (ADFs)
* CMFD accelerated using two-node problem non-linear iteration
* CMFD matrix is solved with the latest linear system solver: BiCGSTAB
* Thermal-hydraulics solutions are obtained by solving mass and energy conservation equations in an enclosed channel
* Three nodal kernels are available:
  * Traditional Finite Difference Method
  * Polynomial Nodal Method (PNM) which is equivalent to Nodal Expansion Method (NEM)
  * Semi-Analytic Nodal Method (SANM)

# User Guides

Here you can find quick and advanced guides on how to use ADPRES. Given you have background on nuclear engineering, **we believe you can create your own ADPRES input within minutes!**
## [Theory and Background](https://imronuke.github.io/ADPRES/method)
## [Source code compilation](https://imronuke.github.io/ADPRES/install)
## [Quick guides](https://imronuke.github.io/ADPRES/quick-guides)
## [Advanced guides](https://imronuke.github.io/ADPRES/card-desc)


# How to give feedbacks
You may raise an issue or contact me at
* muhammad.imron[at]adpoly.ac.ae
* makrus.imron[at]gmail.com

# How to cite
If you find this work helpful and use this work for a publication, you may cite as

**Imron, M. (2019). [Development and verification of open reactor simulator ADPRES](https://doi.org/10.1016/j.anucene.2019.06.049). Annals of Nuclear Energy, 133, 580–588.**


> **"The best of people are those who bring most benefit to the rest of mankind." (THE PROPHET)**
