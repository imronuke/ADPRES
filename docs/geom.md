---
title: %GEOM
theme: _config.yml
filename: geom
---

# %XSEC Card

This card is used to describe the problem geometry in rectangular coordinate system. This card is mandatory. The coordinate system used in ADPRES is shown in the following figure

![alt text](https://github.com/imronuke/ADPRES/blob/master/docs/images/geom_1.png "ADPRES 3D coordinate system")

The point of origin is located at the corner between west, bottom and south sides. The next figure shows the coordinate system seen from top which typically used for two-dimensional problems.

![alt text](https://github.com/imronuke/ADPRES/blob/master/docs/images/geom_2.png "ADPRES 2D coordinate system")

| `%XSEC` | Variable | Description | Remarks |
| --- |
| LINE 1 | NG | Number of groups |  |
|        | NMAT | Number materials |
| LINE 2 | SIGTR(g) | Transport macroscopic XS for group g | Repeat LINE 2 NG times. And again repeat this input segment NMAT times.(See example in the provided sample inputs) |
| SIGA(g) | Absorption macroscopic XS for group g |
| NUF(g) | Nu \* Fission macroscopic XS for group g |
| SIGF(g) | Fission macroscopic XS for group g |
| CHI(g) | Fission neutron spectrum for group g |
| SIGS(g,1:NG) | Scattering macroscopic XS from group g to other groups |

Example:
```
%XSEC
2  5    ! Number of groups and number of materials
! sigtr   siga   nu*sigf sigf   chi   sigs_g1  sigs_g2
0.222222  0.010  0.000  0.000    1.0   0.1922   0.020
0.833333  0.080  0.135  0.135    0.0   0.000    0.7533   ! MAT1 : Outer Fuel
0.222222  0.010  0.000  0.000    1.0   0.1922   0.020
0.833333  0.085  0.135  0.135    0.0   0.000    0.7483   ! MAT2 : Inner Fuel
0.222222  0.0100 0.000  0.000    1.0   0.1922   0.020
0.833333  0.1300 0.135  0.135    0.0   0.000    0.7033   ! MAT3 : Inner Fuel + Control Rod
0.166667  0.000  0.000  0.000    0.0   0.1267   0.040
1.111111  0.010  0.000  0.000    0.0   0.000    1.1011   ! MAT4 : Reflector
0.166667  0.000  0.000  0.000    0.0   0.000    0.040
1.111111  0.055  0.000  0.000    0.0   0.000    0.000    ! MAT5 : Reflector + Control Rod
```
