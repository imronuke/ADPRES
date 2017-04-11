PROGRAM main

USE sdata, ONLY: mode
USE InpOutp, ONLY: ounit, inp_read
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



CALL CPU_TIME(fn)

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) "Total time : ", fn-st, " seconds" 
 
END PROGRAM main