---
title: %GEOM
theme: _config.yml
filename: geom
---

# %GEOM Card

This card is used to describe the problem geometry in rectangular coordinate system. This card is mandatory. The coordinate system used in ADPRES is shown in the following figure

![alt text](https://github.com/imronuke/ADPRES/tree/master/docs/images/geom_1.png "ADPRES 3D coordinate system")

The point of origin is located at the corner between west, bottom and south sides. The next figure shows the coordinate system seen from top which typically used for two-dimensional problems.

![alt text](https://github.com/imronuke/ADPRES/tree/master/docs/images/geom_2.png "ADPRES 2D coordinate system")

| `%GEOM` | Variable | Description | Remarks |
| --- |
| LINE 1 | NX | Number of assemblies along X-direction |  |
|   | NY | Number of assemblies along Y-direction |
|   | NZ | Number of assemblies along Z-direction |
| LINE 2 | XSIZE(1:NX) | Assembly size along X-direction(from west to east) |  |
| LINE 3 | XDIV(1:NX) | Assembly division along X-direction(from west to east) |  |
| LINE 4 | YSIZE(1:NY) | Assembly size along Y-direction(from south to north) |  |
| LINE 5 | YDIV(1:NY) | Assembly division along Y-direction(from south to north) |  |
| LINE 6 | ZSIZE(1:NZ) | Assembly size along Z-direction(from bottom to top) |  |
| LINE 7 | ZDIV(1:NZ) | Assembly division along Z-direction(from bottom to top) |  |
| LINE 8 | NP | Number of different core planar with different material composition |  |
| LINE 9 | ZPLN(1:NZ) | Planar assignment along Z-direction from bottom to top |  |
| LINE 10 | DO j = NY, 1, -1     ASM(1:NX)END DO | Planar material map |  |
| LINE 11 | XEAST | East boundary conditions | 0 = Zero flux |
|   | XWEST | West boundary conditions | 1 = Zero incoming current |
|   | YNORTH | North boundary conditions | 2 = Reflective |
|   | YSOUTH | South boundary conditions |
|   | ZBOTT | Bottom boundary conditions |
|   | ZTOP | Top boundary conditions |

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
