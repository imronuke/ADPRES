---
title: %ITER
theme: _config.yml
filename: iter
---

# %ITER Card

This card can be used to control the iterations in ADPRES calculation

| %ITER | Variable | Description | Remarks |
| --- | --- | --- | --- |
| LINE 1 | NOUT  | Maximum number of outer iteration | Default = 2 |
|        | NIN   | Maximum number of inner iteration per an outer iteration | Default = 500 |
|        | SERC  | Fission source error criteria (relative error) | Default = 1.e-5 |
|        | FERC  | Flux source error criteria (relative error) | Default = 1.e-5 |
|        | NAC   | Outer iteration fission source extrapolation interval | Default = 5 |
|        | NUPD  | Nodal update interval through outer iteration | Default = (NX+NY+NZ)/2.5. Effective only for ANM and PNM [nodal kernel](https://imronuke.github.io/ADPRES/kern) |
|        | TH_NITER  | Maximum number of T-H iteration | Default = 30. Effective only if [`%THER`](https://imronuke.github.io/ADPRES/ther) present  |
|        | NTH  | Maximum number of outer iteration per T-H iteration | Default = 20. Effective only if [`%THER`](https://imronuke.github.io/ADPRES/ther) present |

Example:
```
! Iteration control card
%ITER
1200 3 1.e-5 1.e-5 20 40 20 80
```
