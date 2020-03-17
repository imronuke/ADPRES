---
title: Quick install guide
theme: _config.yml
filename: install
---

# Compiling the source codes
To compile ADPRES source codes you just need a Fortran compiler, that's it. We design ADPRES to be very portable, so it can be installed on any machine.

## Compiling in Ubuntu or other GNU-Linux based OS
In Ubuntu, other GNU-Linux based OS or CYGWIN you can use either gfotran or intel fortran to compile the source codes. You can install gfortran in Ubuntu OS by using command

```
sudo apt install gfortran
```

Then you can download the [ADPRES zip files](https://github.com/imronuke/ADPRES) from Github or clone them (you need to download git first if you don't have one)

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

These command will create an executable file named `adpres`. Now, you can run a test using several examples of inputs file in folder [smpl](https://github.com/imronuke/ADPRES/tree/master/smpl) to see if you had compiled properly. You can run ADPRES using command

```
./adpres [FILE_PATH_NAME]
```

for example, you can run a test by

```
./adpres /home/imronuke/smpl/static/IAEA3Ds
```

For you to be able to execute ADPRES from any folder, you can copy it to `usr/bin` (provided you have an admin privilege)

```
sudo cp adpres /usr/bin
```

And now you can execute ADPRES from the smpl folder directly.

## Compiling in Windows
In Windows you can compile ADPRES using g95 which can be obtained from [here](https://www.fortran.com/wp-content/uploads/2013/05/g95-Mingw_201210.exe).

Then you can download the [ADPRES zip files](https://github.com/imronuke/ADPRES) from Github.

After you install g95 and download ADPRES, go to the [src folder](https://github.com/imronuke/ADPRES/tree/master/src) located inside the ADPRES folder which you had been downloaded, and compile the source codes using g95:

```
g95 -O4 -c mod_data.f90
g95 -O4 -c mod_io.f90
g95 -O4 -c mod_xsec.f90
g95 -O4 -c mod_nodal.f90
g95 -O4 -c mod_cmfd.f90
g95 -O4 -c mod_th.f90
g95 -O4 -c mod_trans.f90
g95 -O4 -c mod_control.f90
g95 -O4 -c ADPRES.f90
g95 *.o -o adpres
```

These command will create an executable file named `adpres`. Now, you can run a test using several examples of inputs file in folder [smpl](https://github.com/imronuke/ADPRES/tree/master/smpl) to see if you had compiled properly. You can run ADPRES using command

```
adpres [FILE_PATH_NAME]
```

for example, you can run a test by

```
adpres C:\Users\imronuke\Downloads\ADPRES-master\smpl\static\IAEA3Ds
```
