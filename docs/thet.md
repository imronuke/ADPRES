---
title: %THET
theme: _config.yml
filename: thet
---

# %THET Card

This card is used to set the value of theta. Only effective for `RODEJECT` mode.

| `%THET` | Variable | Description | Available options |
| --- |
| LINE 1 | theta | Theta value | Theta value should between (or equal to) 0.01 to 1.0. Default is 1.0  or correspond to fully-implicit method|

Example:
```
%THET
0.5   !This correspond to Crank-Nicholson method
```
