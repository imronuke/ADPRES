---
title: %KERN
theme: _config.yml
filename: kern
---

# %MODE Card

This card is used to choose other nodal kernel. By default, ADPRES uses Semi-Analytic Nodal Method.

| `%MODE` | Variable | Description | Available options |
| --- |
| LINE 1 | KERN | Nodal kernel | `FDM`  : Finite Difference Method<br>`PNM`   : Polynomial Nodal Method (equivalent to Nodal Expansion Method)<br>`SANM` : Semi-Analytic Nodal Method |

Example:
```
%KERN
PNM
```
