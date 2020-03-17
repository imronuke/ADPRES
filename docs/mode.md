---
title: %MODE
theme: _config.yml
filename: mode
---

# %MODE Card

This card is used to describe the calculation mode. It has several options as follow

| `%MODE` |
| --- |
| LINE 1 | MODE | Calculation mode | FORWARD  : Forward Calculation
ADJOINT   : Adjoint Calculation
FIXEDSRC : Fixed Source Calculation
RODEJECT : Rod ejection and/or insertion mode (transient problems) |

Example:
```
! Mode card
%MODE
FORWARD
```
