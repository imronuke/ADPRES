#!/bin/bash
# start

echo "ADPRES compilation starts"

set FORT=gfotran

echo " "
echo "Compiling mod_data.f90"
gfortran -O4 -c mod_data.f90
echo "Compiling mod_io.f90"
gfortran -O4 -c mod_io.f90
echo "Compiling mod_xsec.f90"
gfortran -O4 -c mod_xsec.f90
echo "Compiling mod_nodal.f90"
gfortran -O4 -c mod_nodal.f90
echo "Compiling mod_cmfd.f90"
gfortran -O4 -c mod_cmfd.f90
echo "Compiling mod_th.f90"
gfortran -O4 -c mod_th.f90
echo "Compiling mod_trans.f90"
gfortran -O4 -c mod_trans.f90
echo "Compiling mod_control.f90"
gfortran -O4 -c mod_control.f90
echo "Compiling ADPRES.f90"
gfortran -O4 -c ADPRES.f90
echo "Combining all together"
gfortran *.o -o adpres

echo " "
echo "Copy adpres to usr/bin"
sudo cp adpres /usr/bin

echo " "
echo "Get rid unnecessary files"
rm *.o *.mod adpres

echo " "
echo "ADPRES successfully compiled"
