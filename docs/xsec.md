---
title: %XSEC
theme: _config.yml
filename: xsec
---

# `%MODE` Card

This card is used to describe the cross section data. This card is conditional : either `%XSEC` or `XTAB` card must be present

| %XSEC |
| --- |
| LINE 1 | NG | Number of groups | 2  4    ! Number of groups and number of materials |
| NMAT | Number materials |
| LINE 2 | SIGTR(g) | Transport macroscopic XS for group g | Repeat LINE 2 NG times. And again repeat this input segment NMAT times.(See example in the provided sample inputs) |
| SIGA(g) | Absorption macroscopic XS for group g |
| NUF(g) | Nu \* Fission macroscopic XS for group g |
| SIGF(g) | Fission macroscopic XS for group g |
| CHI(g) | Fission neutron spectrum for group g |
| SIGS(g,1:NG) | Scattering macroscopic XS from group g to other groups |
