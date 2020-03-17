---
title: %MODE
theme: _config.yml
filename: mode
---

# %MODE Card

This card is used to describe the calculation mode. It has several options as follow

| `%MODE` | Variable | Description | Available options |
| --- |
| LINE 1 | MODE | Calculation mode | `FORWARD`  : Forward Calculation<br>`ADJOINT`   : Adjoint Calculation<br>`FIXEDSRC` : Fixed Source Calculation<br>`BCSEARCH`: Critical boron concentration search<br>`RODEJECT` : Rod ejection and/or insertion mode (transient problems) |

Example:
```
! Mode card
%MODE
FORWARD
```
