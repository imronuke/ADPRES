PROGRAM main

USE sdata, ONLY: mode
USE InpOutp, ONLY: ounit, inp_read, bwrst, w_rst
USE nodal

IMPLICIT NONE

REAL :: st, fn


CALL CPU_TIME(st)

CALL inp_read()

CALL Init()

SELECT CASE(mode)
    CASE('FORWARD')
	    CALL forward()
    CASE('ADJOINT')		
        CALL adjoint()
	CASE DEFAULT
        CALL fixedsrc()	
END SELECT

! Write Restart File if required
IF (bwrst == 1)	CALL w_rst()


CALL CPU_TIME(fn)

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) "Total time : ", fn-st, " seconds" 
 
END PROGRAM main