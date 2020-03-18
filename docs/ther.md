---
title: %THER
theme: _config.yml
filename: ther
---

# %THER Card

This card is used to set the T-H parameters. It also activates T-H feedback.

| `%THER` | Variable | Description | Remarks |
| --- |
| LINE 1 | theta | Theta value | Theta value should between (or equal to) 0.01 to 1.0. Default is 1.0  or correspond to fully-implicit method|

Example:
```
! THERMAL-HYDRAULIC CARD
%THER
100.                                ! Percent power in %
891.25e6                            ! Reactor thermal power for quarter geometry in Watt
560.  82.121243523                  ! Inlet coolant temp. (Kelvin) and Fuel Assembly Mass flow rate (kg/s)
3.951E-03  5.9E-05  5.73E-04  1.26E-2 ! Fuel meat rad., gap thinkness, cladding thinkness and pin picth (m)
264  25                               ! Number of fuel pin and guide tubes
0.0                                   ! FRACTION OF HEAT DEPOSITED IN COOLANT
1
```
