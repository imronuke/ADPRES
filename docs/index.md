---
title: ADPRES
theme: _config.yml
filename: index
---

# ADPRES

Abu Dhabi Polytechnic Reactor Simulator (ADPRES) is a reactor core simulator that solves static and transient diffusion equation for two or three dimensional reactor problems in Cartesian geometry. Currently, ADPRES uses Semi-Analytic Nodal Method (SANM) to spatially discretized the neutron diffusion equation.

ADPRES is a great learning tool for reactor theory classes, and we have been striving hard to make the input is easy to create. ADPRES' main objective is to make all nuclear engineering students have access on reactor simulator code for them to learn. It is open and completely free, so everyone has access to the source code and modify for his/her own purposes.

ADPRES features:
* Solves both static and transient core problems with or without TH feedback
* CMFD accelerated using two-node problem non-linear iteration
* Perform branch calculations. An example of library format can be seen [here]
* Three nodal kernels are available:
  * Traditional Finite Difference Method
  * Polynomial Nodal Method (PNM)
  * Semi-Analytic Nodal Method
* It can handle multi-group of neutron energy
* Performs forward, adjoint and fixed-source calculations
* It handle calculations with Assembly Discontinuity Factors (ADFs)
* Critical boron concentration search
* Rod ejection simulation
* Theta-method to solve transient problems

# Quick Guides

Here you can find quick guides on how to use ADPRES. We believe you can create your own ADPRES input within minutes, given you have background on nuclear engineering

* [Source code compilation](https://imronuke.github.io/ADPRES/install)
* [Quick guides](https://imronuke.github.io/ADPRES/quick-guides)
* Advanced guides
  * General rules
  * Cards Description
