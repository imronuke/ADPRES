---
title: Quick install guide
theme: jekyll-theme-Cayman
filename: install
---

# Compiling the source codes
To compile ADPRES source codes you just need a Fortran compiler, that's it. We design ADPRES to be very portable, so it can be installed on any machine.

## Compiling in Ubuntu or other GNU-Linux based OS
In UNIX environment such as GNU-LINUX or CYGWIN you can use either gfotran or intel fortran to compile the source codes. You can install gfortran in Ubuntu OS by using command

```
sudo apt install gfortran
```

Then you can download the [ADPRES source codes](https://github.com/imronuke/ADPRES) from Github or clone them

```
git clone https://github.com/imronuke/ADPRES.git
```

In a machine where the gfortran is already installed, go to the [src folder](https://github.com/imronuke/ADPRES/tree/master/src) located inside the ADPRES folder which you had been downloaded or cloned, and compile the source codes by using command:

```
gfortran -O4 -c mod_data.f90
gfortran -O4 -c mod_io.f90
gfortran -O4 -c mod_xsec.f90
gfortran -O4 -c mod_nodal.f90
gfortran -O4 -c mod_cmfd.f90
gfortran -O4 -c mod_th.f90
gfortran -O4 -c mod_trans.f90
gfortran -O4 -c mod_control.f90
gfortran -O4 -c ADPRES.f90
gfortran *.o -o adpres
```

These command will create executable file named `adpres`. For you to be able to execute ADPRES from any folder, you can copy it to `usr/bin`

```
sudo cp adpres /usr/bin
```
