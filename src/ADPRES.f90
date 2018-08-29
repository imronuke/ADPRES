PROGRAM main

USE sdata, ONLY: mode
USE InpOutp, ONLY: ounit, inp_read, bwrst, w_rst
USE nodal, ONLY: forward, adjoint, fixedsrc, init
USE trans, ONLY: rod_eject

IMPLICIT NONE

REAL :: st, fn


CALL CPU_TIME(st)

CALL inp_read()

CALL Init()

SELECT CASE(mode)
    CASE('FIXEDSRC')
        CALL fixedsrc()
    CASE('ADJOINT')        
        CALL adjoint()
    CASE('RODEJECT')        
        CALL rod_eject()
    CASE DEFAULT
        CALL forward()
END SELECT

! Write Restart File if required
IF (bwrst == 1)    CALL w_rst()


CALL CPU_TIME(fn)

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) "Total time : ", fn-st, " seconds" 
 
END PROGRAM main