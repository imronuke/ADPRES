# Background and Methodology   

Current version of ADPRES uses current available methods; therefore, rather than re-explaining the methods, the author would like to just mention the references of the methods implemented in ADPRES.

The ADPRES development began on the early 2017 and was motivated by the hurdles to obtain such similar codes in region where the author was working. The initial version of ADPRES used Nodal Expansion Method (NEM) based on the response matrix formulations [1,2,3]. The transverse integrated leakages were approximated by the Quadratic Transverse Leakage Approximation (QTLA)[4]. The transient and and thermal module were also added to enable ADPRES handle time-dependent problems with thermal-hydraulics (T-H) feedback. The transient diffusion equation was solved with fully implicit method, and T-H solutions were obtained by solving mass and energy conservation equations in an enclosed channel, following NODAL 3 computer code [5] with slight modifications. This version of ADPRES was published in Annals of Nuclear Energy in 2019 [6]. However, this version encountered slow performance compared to other modern nodal simulators.

In the early 2020, the ADPRES was revamped to implement CMFD acceleration with two-node problems non-linear iteration procedures [7] for the sake of rapid calculations. Initially, ADPRES implemented Polynomial Nodal Method [8,5], but then upgraded to Semi-Analytic Nodal Method [9,10] for better accuracy. Also, in the transient calculations, the delayed terms and other terms that do not appear in static calculations are included with transverse leakages in the calculation of transverse moments to save further save memory storage [11]. Users also have option to perform exponential flux [12] transformation for time-dependent problems which might be useful for rod ejection simulation from HZP. By implementing CMFD acceleration, the current ADPRES version is able to run time-dependent problems with running time 10-20 times faster than previous version of ADPRES.

# Acknowledgement

This ADPRES development would be impossible without the God's Mercy and the works done by those mentioned in the references. We would like to thank them and other people who contributed on their works. We also would like to thank to other people who directly or indirectly contributed to this work:

* Dr. Ali Al Naqbi
* Dr. Anthony Hechanova
* Prof. Nam Zin Cho
* Dr. Alexander Agung
* Dr. Andang Widiharto
* Liem Peng Hong, PhD
* All my colleagues and friends.

# References
[1] Bennewitz F., Finnemann H. and Moldaschl H.(1975)  Solution of the multidimensional neutron diffusion equa-tion by nodal expansion. Proc. Conf. Computational Methods in Nuclear Engineering, p. 1-99, CONF-750413.

[2] Lawrence, R.D., (1986) Progress in nodal methods for the solution of the neutron diffusion and transport equations, Progress in Nuclear Energy, Vol. 17, No.3, pp. 271-301.

[3] Okumura, K., (1998) MOSRA-Light: High speed three-dimensional nodal diffusion code for vector computers, JAEA-Research 98-025. (in Japanese)

[4] Finnemann, H., Bennewitz F. and Wagner M. R., (1977) Interface current techniques for multidimensional reactor calculations, Atomkernenergie, Vol. 30, pp. 123-128.

[5] Liem P.H., et al., (2010) NODAL3: A Nodal Neutron Diffusion Code Version 2.0 User Guides (unpublihsed)

[6] Imron, M. (2019). Development and verification of open reactor simulator ADPRES. Annals of Nuclear Energy, 133, 580–588. https://doi.org/10.1016/j.anucene.2019.06.049

[7] Smith, K. S. (1984) Nodal Method Storage Reduction by Nonlinear Iteration. Transactions of American Nuclear Society 44, 265.

[8] Zimin V.G. and Ninokata, H., (1997) Nonlinear Iteration Procedure Based on Legendre Polynomials, Trans. Am. Nucl. Soc., 76, 162.

[9] Zimin, V. G., & Ninokata, H. (1998). Nodal neutron kinetics model based on nonlinear iteration procedure for LWR analysis. Annals of Nuclear Energy, 25(8), 507–528. https://doi.org/10.1016/S0306-4549(97)00078-9

[10] Zimin V.G. and Ninokata, H., (1997) Polynomials ans Semi-Analytic Nodal Methods For Non-Linear Iteration Procedures, Trans. Am. Nucl. Soc., 76, 162.

[11] Engrand, P. R., Maldonado, G. I., A1-Chalabi, R. M. and Turinsky, P. J. (1992) Non-Linear Iterative Strategy for NEM Refinement and Extension. Transactions of American Nuclear Society 65, 221.

[12] Hendricks, J.S., “Finite difference solution of the time dependent neutron group diffusion equations”, Thesis, Department of Nuclear Engineering, Massachusetts Institute of Technology, MITNE-176 (1975).
