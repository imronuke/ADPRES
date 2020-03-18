---
title: ADPRES
theme: _config.yml
filename: index
---

# ADPRES

Abu Dhabi Polytechnic Reactor Simulator (ADPRES) is an open reactor core simulator that solves static and transient diffusion equation for two or three dimensional reactor problems in Cartesian geometry. Currently, ADPRES uses Semi-Analytic Nodal Method (SANM) to spatially discretised the neutron diffusion equation. While theta method is used for the time discretisation.

ADPRES is a great learning tool for reactor theory classes, and we have been striving hard to make the input is easy to create. ADPRES' main objective is to make all nuclear engineering students have access on reactor simulator code for them to learn. It is open and completely free, so everyone has access to the source code and modify for his/her own purposes.

ADPRES features:
* Input is straightforward, modular and in a free-format form
* Solves both static and transient core problems with or without TH feedback
* Performs forward, adjoint and fixed-source calculations
* Perform branch calculations. An example of the library format can be seen [here]
* Critical boron concentration search
* Rod ejection simulation or Reactivity Initiated Accident (RIA)
* CMFD accelerated using two-node problem non-linear iteration
* CMFD matrix is solved with the latest matrix solver: BiCGSTAB
* Three nodal kernels are available:
  * Traditional Finite Difference Method
  * Polynomial Nodal Method (PNM)
  * Semi-Analytic Nodal Method (SANM)
* It can handle multi-group of neutron energy
* It handle calculations with Assembly Discontinuity Factors (ADFs)

# User Guides

Here you can find quick guides on how to use ADPRES. Given you have background on nuclear engineering, we believe you can create your own ADPRES input within minutes!

* [Source code compilation](https://imronuke.github.io/ADPRES/install)
* [Quick guides](https://imronuke.github.io/ADPRES/quick-guides)
* [Advanced guides](https://imronuke.github.io/ADPRES/card-desc)
