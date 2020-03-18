---
title: %ITER
theme: _config.yml
filename: iter
---

# %ITER Card

This card can be used to control the iterations in ADPRES calculation

| %ITER | Variable | Description | Remarks |
| --- | --- | --- | --- |
| LINE 1 | NOUT  | Maximum number of outer iteration |
|        | NIN   | Maximum number of inner iteration per an outer iteration |
|        | SERC  | Fission source error criteria (relative error) |
|        | FERC  | Flux source error criteria (relative error) |
|        | NAC   | Outer iteration fission source extrapolation interval |
|        | NUPD  | Nodal update interval through outer iteration |
|        | TH_NITER  | Maximum number of T-H iteration |
|        | NTH  | Maximum number of outer iteration per T-H iteration |

Example:
```
! Iteration control card
%ITER
1200 3 1.e-5 1.e-5 20 40 20 80
```
