PROGRAM main

USE sdata, ONLY: mode, negxs
USE InpOutp, ONLY: ounit, inp_read, bwrst, w_rst, bther
USE nodal, ONLY: forward, adjoint, fixedsrc, init
USE trans, ONLY: rod_eject, trod_eject
USE th, ONLY: cbsearch, cbsearcht

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
        IF (bther == 0) THEN
            CALL rod_eject()
  	    ELSE
  		      CALL trod_eject()
  		  END IF
    CASE('BCSEARCH')
  	    IF (bther == 0) THEN
            CALL cbsearch()
  	    ELSE
  		      CALL cbsearcht()
  		  END IF
    CASE DEFAULT
        CALL forward()
END SELECT

! Write Restart File if required
IF (bwrst == 1)    CALL w_rst()

IF (negxs) THEN
    WRITE(ounit,*)
    WRITE(ounit,*) "  WARNING: SOME NEGATIVE CXs (DUE TO CR INSERTION) ARE SUPPRESSED TO ZERO IN SUBROUTINE CROD_UPDT"
END IF


CALL CPU_TIME(fn)

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) "Total time : ", fn-st, " seconds"

END PROGRAM main
