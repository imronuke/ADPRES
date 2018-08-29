MODULE th

!=========================
! Thermal-hydraulics Module to solve Thermal-hydraulics related problems
! =======================

IMPLICIT NONE

SAVE

CONTAINS

SUBROUTINE th_iter(xbcon, xbpos)

!
! Purpose:
!    To calculate adjoint flux
!

USE sdata, ONLY: nnod, ftem, mtem, cden, npow, pow, zdel, node_nf, ix, iy, iz
USE nodal, ONLY: nodal_coup4, outer4, PowDis
USE InpOutp, ONLY: XS_updt, ounit

IMPLICIT NONE

REAL, INTENT(IN) :: xbcon                  ! Boron concentration input
REAL, DIMENSION(:), INTENT(IN) :: xbpos    ! CR Bank position input

REAL, DIMENSION(nnod) :: otem
REAL, DIMENSION(nnod) :: nlpow
REAL :: th_err
REAL :: leng
INTEGER :: n, niter

INTEGER :: nout = 10    ! Maximum number of outer iteration
INTEGER :: nac  = 5     ! Acceleration interval

th_err = 1.
niter = 0
DO
    niter = niter + 1

    ! Save old moderator temp
	otem = ftem

    ! Update XS
    CALL XS_updt(xbcon, ftem, mtem, cden, xbpos)

	! Update nodal couplings
    CALL nodal_coup4()

	! Perform outer inner iteration
    CALL outer4(0, nout, nac)

	! Calculate power density
	CALL PowDis(npow)
	! Calculate linear power density (W/m)
	DO n = 1, nnod
	    npow(n) = npow(n) * pow                             ! Power for each nodes/channel

	    leng = zdel(iz(n)) * 0.01                           ! Convert to meter
	    nlpow(n) = npow(n) / (leng * node_nf(ix(n),iy(n)))  ! Linear power density for each fuel pin
	END DO

	! Update fuel, moderator temp. and coolant density
	CALL th_upd(nlpow)

	th_err = MAXVAL(ABS(ftem - otem))
	WRITE(*,*) th_err
	IF ((th_err < 0.01) .OR. (niter == 5)) EXIT
END DO


END SUBROUTINE th_iter


SUBROUTINE th_upd(xnlpow)

!
! Purpose:
!    To calculate adjoint flux
!

USE sdata, ONLY: mtem, cden, ftem, tin, xyz, npow, mflow, nyy, nzz, ystag, &
                 rf, tg, tc, vdel, yfstn, pi
USE InpOutp, ONLY : ounit

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: xnlpow

INTEGER :: i, j, k
REAL :: cp
REAL :: dTf, dTg, dTc, dTl
REAL :: hs = 2.e4, kc, hg = 1.e4, kf
REAL :: cT, gT, fT, clT
REAL :: alp = 0.7

DO j = 1, nyy
    DO i = yfstn(j)%smin, yfstn(j)%smax

	    ! For the most bottom nodes
	    CALL getcp(tin,cp)                                             ! Get specific heat capacity
	    mtem(xyz(i,j,1)) = tin + npow(xyz(i,j,1)) / (mflow(i,j) * cp)  ! Calculate coolant temp.
		CALL getrho(mtem(xyz(i,j,1)), cden(xyz(i,j,1)))                ! Get node coolant density

		dTl = xnlpow(xyz(i,j,1)) / (hs * 2. * pi * (rf + tc))! Delta T coolant
		cT = mtem(xyz(i,j,1)) + dTl                          ! Temp. on clad surface

		kc = getkc(cT)                                       ! Get thermal conductivity for clad
		dTc = xnlpow(xyz(i,j,1)) * tc / (2. * pi * rf * kc)  ! Delta T clad
		gT = cT + dTc                                        ! Temp. on gap surface

		dTg = xnlpow(xyz(i,j,1)) / (2. * pi * rf * hg)       ! Delta gap
		fT = gT + dTg                                        ! Temp. of fuel surface
		kf = getkf(ftem(xyz(i,j,1)))                         ! Get thermal conductivity for fuel

		dTf = xnlpow(xyz(i,j,1)) / (4. * pi * kf)              ! Delta fuel
		clT = fT + dTf
		ftem(xyz(i,j,1)) = (1. - alp) * clT + alp * fT
		! WRITE(ounit,*) xnlpow(xyz(i,j,1)) * 0.01, mtem(xyz(i,j,1)), fT, clT, ftem(xyz(i,j,1))
	    DO k = 2, nzz
	        ! For the remaining nodes
	        CALL getcp(mtem(xyz(i,j,k-1)),cp)                                          ! Get specific heat capacity
	        mtem(xyz(i,j,k)) = mtem(xyz(i,j,k-1)) + npow(xyz(i,j,k)) / (mflow(i,j)*cp) ! Calculate coolant temp.
		    CALL getrho(mtem(xyz(i,j,k)), cden(xyz(i,j,k)))                            ! Get node coolant density

		    dTl = xnlpow(xyz(i,j,k)) / (hs * 2. * pi * (rf + tc))! Delta T coolant
		    cT = mtem(xyz(i,j,k)) + dTl                          ! Temp. on clad surface

		    kc = getkc(cT)                                       ! Get thermal conductivity for clad
		    dTc = xnlpow(xyz(i,j,k)) * tc / (2. * pi * rf * kc)    ! Delta T clad
		    gT = cT + dTc                                        ! Temp. on gap surface

		    dTg = xnlpow(xyz(i,j,k)) / (2. * pi * rf * hg)       ! Delta gap
		    fT = gT + dTg                                        ! Temp. of fuel surface
		    kf = getkf(ftem(xyz(i,j,k)))                         ! Get thermal conductivity for fuel

		    dTf = xnlpow(xyz(i,j,k)) / (4. * pi * kf)            ! Delta fuel
		    clT = fT + dTf
			ftem(xyz(i,j,k)) = (1.-alp) * clT + alp * fT
			! WRITE(ounit,*) xnlpow(xyz(i,j,k)) * 0.01, mtem(xyz(i,j,k)), fT, clT, ftem(xyz(i,j,k))
		END DO
		! WRITE(ounit,*)
	END DO
END DO

! CALL par_ave(ftem, dTf)
CALL par_ave(mtem, dTg)
! CALL par_ave(cden, dTc)

! WRITE(ounit,*)
! WRITE(ounit, *) dTf, dTg, dTc

! WRITE(*,*) MAXVAL(mtem), dtg

! STOP

END SUBROUTINE th_upd


SUBROUTINE par_ave(par, ave)
!
! Purpose:
!    To calculate average fuel tem, coolant tem, and density
!

USE sdata, ONLY: vdel, nnod, yfstn, xyz, nzz, nyy

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: par
REAL, INTENT(OUT) :: ave
REAL :: dum, dum2
INTEGER :: i, j, k

dum = 0.; dum2 = 0.
DO k = 1, nzz
    DO j = 1, nyy
	    DO i = yfstn(j)%smin, yfstn(j)%smax
		    dum = dum + par(xyz(i,j,k)) * vdel(xyz(i,j,k))
			dum2 = dum2 + vdel(xyz(i,j,k))
        END DO
	END DO
END DO

ave = dum / dum2

END SUBROUTINE par_ave


SUBROUTINE getcp(t,cp)
!
! Purpose:
!    To get specific heat capacity for given coolant temp. from steam table
!

USE sdata, ONLY: thunit
USE InpOutp, ONLY : ounit

IMPLICIT NONE

REAL, INTENT(IN) :: t
REAL, INTENT(OUT) :: cp
REAL :: t1, p1, rho1, cp1
REAL :: t2, p2, rho2, cp2

IF ((t < 295.15) .OR. (t > 641.15)) THEN
    WRITE(ounit,*) '  Coolant temp. : ', t
    WRITE(ounit,*) '  ERROR : MODERATOR TEMP. IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
	WRITE(ounit,*) '  CHECK INPUT MASS FLOW RATE OR POWER'
	STOP
END IF

REWIND(thunit)

READ(thunit,*) t2, p2, rho2, cp2
DO
	t1 = t2
	p1 = p2
	rho1 = rho2
	cp1 = cp2
	READ(thunit,*) t2, p2, rho2, cp2
	IF ((t >= t1) .AND. (t <= t2)) THEN
	    cp = cp1 + (t - t1) / (t2 - t1) * (cp2 - cp1)
	    EXIT
	END IF
END DO


END SUBROUTINE getcp


SUBROUTINE getrho(t,rho)
!
! Purpose:
!    To get density for given coolant temp. from steam table
!

USE sdata, ONLY: thunit
USE InpOutp, ONLY : ounit

IMPLICIT NONE

REAL, INTENT(IN) :: t
REAL, INTENT(OUT) :: rho
REAL :: t1, p1, rho1, cp1
REAL :: t2, p2, rho2, cp2

IF ((t < 295.15) .OR. (t > 641.15)) THEN
    WRITE(ounit,*) '  Coolant temp. : ', t
	WRITE(ounit,*) '  ERROR : MODERATOR TEMP. IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
	WRITE(ounit,*) '  CHECK INPUT MASS FLOW RATE OR POWER'
	STOP
END IF

REWIND(thunit)

READ(thunit,*) t2, p2, rho2, cp2
DO
	t1 = t2
	p1 = p2
	rho1 = rho2
	cp1 = cp2
	READ(thunit,*) t2, p2, rho2, cp2
	IF ((t >= t1) .AND. (t <= t2)) THEN
	    rho = rho1 + (t - t1) / (t2 - t1) * (rho2 - rho1)
	    EXIT
	END IF
END DO


END SUBROUTINE getrho


REAL FUNCTION getkc(t)
!
! Purpose:
!    To calculate thermal conductivity of cladding
!

IMPLICIT NONE

REAL, INTENT(IN) :: t

getkc = 7.51 + 2.09E-2*t - 1.45E-5*t**2 + 7.67E-9*t**3

END FUNCTION getkc


REAL FUNCTION getkf(t)
!
! Purpose:
!    To calculate thermal conductivity of fuel
!

IMPLICIT NONE

REAL, INTENT(IN) :: t

getkf = 1.05 + 2150. / (t - 73.15)

END FUNCTION getkf


REAL FUNCTION Reynold(xmflow, xden)
!
! Purpose:
!    To calculate Reynold Number
!

USE sdata, ONLY: dh, nfpin, farea, ngt

IMPLICIT NONE

REAL, INTENT(IN) :: xden, xmflow
REAL, PARAMETER :: nu = 0.294e-6   ! Water kinematic viscosity at 100 oC (m2/s)
REAL :: dia
REAL :: cvelo, cmflow

WRITE(*,*) ' FLow Rate : ', xmflow
WRITE(*,*) ' Density : ', xden
cmflow = xmflow / (nfpin + ngt)                ! Calculate sub-channel mass flow rate
cvelo = cmflow / (farea * xden * 1000.)        ! Calculate flow velocity (density converted to kg/m3)
Reynold = cvelo * dh / nu                      ! Calculate Reynolds Number


END FUNCTION Reynold


REAL FUNCTION h2o_hs(Rex)
!
! Purpose:
!    To calculate Water thermal conductivity
!

USE sdata, ONLY: dh

IMPLICIT NONE

REAL, INTENT(IN) :: Rex
REAL, PARAMETER :: Pr = 0.81  !Prandtl Number at about 300 oC
REAL, PARAMETER :: tc = 0.68  ! themal conductivity of water at 100 oC (W/m2.K)
REAL :: Nu

Nu = 0.023*(Pr**0.4)*(Rex**0.8)         ! Calculate Nusselt Number
h2o_hs = (tc / dh) * Nu                 ! Calculate heat transfer coefficient

END FUNCTION h2o_hs



END MODULE th
