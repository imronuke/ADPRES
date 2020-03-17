---
title: General Rules
theme: _config.yml
filename: general-rules
---

# General Rules
Some general rules for the ADPRES inputs:
1.	Comments are marked by `!`. Example:
```
! COMPOSITION 1
0.20574  0.00654  0.00415  0.00415  1.0  0.0  0.01462  !Group 1
0.68866  0.04850  0.06099  0.06099  0.0  0.0  0.00000  !Group 2
```

2.	Currently ADPRES has several input cards. Cards’ name shall be uppercase and marked by `%`. Example:
```
%MODE
FORWARD
%XSEC    ! Cross section card                                                                                                                                  
2 4      ! Number of groups and number of materials
...
...
%GEOM    ! Geometry card
12 12 2  !nx, ny, nz
...
...
```

3.	Numbers can be repeated using `*` mark. For example
```
10.0 8*20.0  is equivalent to  10.0 20.0 20.0 20.0 20.0 20.0 20.0 20.0 20.0
```
