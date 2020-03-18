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

There are some conditions for some calculation mode.
1. When calculation mode is `FIXEDSRC`, then the [`ESRC`](https://imronuke.github.io/ADPRES/esrc) card must be present
2. When calculation mode is `BCSEARCH`, then the [`cbcs`](https://imronuke.github.io/ADPRES/cbcs) card must be present

Example:
```
! Mode card
%MODE
FORWARD
```
