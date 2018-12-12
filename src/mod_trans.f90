MODULE trans

!=========================
! Transient Module to solve transient diffusion problems
! Using Aadiabatic approach adapated from paper:
! Neutronic Modeling for Modular High Temperature Pebble Bed Reactor during Reactivity Accident
! Author: Peng Hong LIEM and Hiroshi SEKIMOTO
! Journal of NucLEAR SciENCE and TECHNOLOGY, 29(8], pp. 805-812 (August 1992).
! =======================


IMPLICIT NONE

SAVE

REAL, DIMENSION(:,:), ALLOCATABLE :: af      ! Adjoint Flux
REAL, DIMENSION(:), ALLOCATABLE :: beta, C   ! beta (point kinetic), neutron precusor density (point kinetic)
REAL :: ptbet ! Total beta (point kinetic)
REAL, DIMENSION(:,:), ALLOCATABLE :: xvdum

CONTAINS

SUBROUTINE rod_eject2()

!
! Purpose:
!    To perform rod ejection simulation
!

USE sdata, ONLY: ng, nnod, sigr, nuf, sigs, f0, &
                 iBeta, lamb, nf, nout, nac, &
                 ttot, tdiv, tstep1, tstep2, &
                 ix, iy, iz, zdel, vdel, Ke, &
                 bcon, ftem, mtem, cden, &
                 fbpos, bpos, tmove, bspeed, mdir, &
                 nout, nac, nb, xnuf, dnuf, sigr, velo, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, &
                 fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2, &
                 c0, cx1, cy1, cz1, cx2, cy2, cz2, tbeta, &
                 ct, ctx1, cty1, ctz1, ctx2, cty2, ctz2
USE InpOutp, ONLY: XS_updt, bther, ounit
USE nodal, ONLY: nodal_coup4, outer4, outertr, powdis2

IMPLICIT NONE

REAL, DIMENSION(nnod, ng) :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Flux at previous time step

REAL :: t1, t2, lat
REAL :: tpow1, tpow2
INTEGER :: n, i, j, g, imax, step
LOGICAL :: stime

! Allocate fission source
ALLOCATE (fs0(nnod), fsx1(nnod), fsy1(nnod), fsz1(nnod))
ALLOCATE (fsx2(nnod), fsy2(nnod), fsz2(nnod))

! Allocate precusor density
ALLOCATE (c0(nf,nnod),cx1(nf,nnod),cy1(nf,nnod),cz1(nf,nnod))
ALLOCATE (cx2(nf,nnod),cy2(nf,nnod),cz2(nf,nnod))
ALLOCATE (ct(nf,nnod),ctx1(nf,nnod),cty1(nf,nnod),ctz1(nf,nnod))
ALLOCATE (ctx2(nf,nnod),cty2(nf,nnod),ctz2(nf,nnod))

! Calculate forward flux at t=0 and check if keff=1
CALL XS_updt(bcon, ftem, mtem, cden, bpos)
CALL nodal_coup4()
CALL outer4(0)

! If K-EFF NOT EQUAL TO 1.0
IF (ABS(Ke - 1.) > 1.e-5) THEN
   WRITE(ounit, *)
   WRITE(ounit, '(A46,F9.6)') '  INITIAL MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke
   WRITE(ounit, *) '  WARNING: THE STEADY STATE K-EFF IS NOT EQUAL TO 1.0'
   WRITE(ounit, *) '  AND NOW IT IS FORCED TO 1.0 BY MODIFYING THE nu*sigf CROSS SECTIONS '
   WRITE(ounit, *)
   DO i = 1, 10
      xnuf = xnuf / Ke
      dnuf = dnuf / Ke
      CALL XS_updt(bcon, ftem, mtem, cden, bpos)
      CALL outer4(0)
      IF (ABS(Ke-1.0) < 1.e-5) EXIT
   END DO
   IF (i == 10) STOP "K-EFF STILL NOT EQUAL TO ONE. ADPRES IS STOPPING"
END IF

! ReCalculate forward flux at t=0 after keff fixed
CALL outer4(0)

! Calculate power
CALL powdis2(f0,tpow1)

! Calculate Initial precursor density
DO n = 1, nnod
   DO j = 1, nf
      c0(j,n) = iBeta(j) * fs0(n) / lamb(j)
      cx1(j,n) = iBeta(j) * fsx1(n) / lamb(j)
      cy1(j,n) = iBeta(j) * fsy1(n) / lamb(j)
      cz1(j,n) = iBeta(j) * fsz1(n) / lamb(j)
      cx2(j,n) = iBeta(j) * fsx2(n) / lamb(j)
      cy2(j,n) = iBeta(j) * fsy2(n) / lamb(j)
      cz2(j,n) = iBeta(j) * fsz2(n) / lamb(j)
   END DO
END DO

! Total beta
tbeta = 0.
DO j = 1, nf
  tbeta = tbeta + iBeta(j)
END DO

WRITE(ounit, *)
WRITE(ounit, *) " TRANSIENT RESULTS :"
WRITE(ounit, *)
WRITE(ounit, *) " Step  Time(s)  React.($)   Rel. Power   CR Bank Pos. (1-end)"
WRITE(ounit, *) "--------------------------------------------------------------"
WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') 0, 0., 0., 1.0, (bpos(n), n = 1, nb)

! Start transient calculation
step = 0
t2 = 0.
imax = CEILING(tdiv/tstep1)
stime = .FALSE.

! First Time Step
DO i = 1, imax

    step = step + 1
    t1 = t2
    t2 = t1 + tstep1

    IF (t2 > tdiv) THEN
        t2 = tdiv
        stime = .TRUE.
    END IF

    ! Rod bank changes
    DO n = 1, nb
        IF (mdir(n) == 1) THEN   ! If CR moving down
            IF (t2-tmove(n) > 1.d-5 .AND. fbpos(n)-bpos(n) < 1.d-5) THEN
                bpos(n) = bpos(n) - tstep1 *  bspeed(n)
                IF (bpos(n) < fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE IF (mdir(n) == 2) THEN ! If CR moving up
            IF (t2-tmove(n) > 1.d-5 .AND. fbpos(n)-bpos(n) > 1.d-5) THEN
                bpos(n) = bpos(n) + tstep1 *  bspeed(n)
                IF (bpos(n) > fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE
            CONTINUE
        END IF
     END DO

    ! Calculate xsec after pertubation
    CALL XS_updt(bcon, ftem, mtem, cden, bpos)

    ! Save the previous flux and DN precusor density
    ft = f0
    ftx1 = fx1; fty1 = fy1; ftz1 = fz1
    ftx2 = fx2; fty2 = fy2; ftz2 = fz2
    ct = c0
    ctx1 = cx1; cty1 = cy1; ctz1 = cz1
    ctx2 = cx2; cty2 = cy2; ctz2 = cz2

    ! Modify removal xsec
    DO g = 1, ng
       DO n = 1, nnod
          sigr(n,g) = sigr(n,g) + 1. / (velo(g) * tstep1)
       END DO
    END DO

    ! Transient calculation
    CALL nodal_coup4()
    CALL outertr(200, tstep1, ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2)

    ! Update fission source
    fs0 = 0.; fsx1 = 0.; fsy1 = 0.; fsz1 = 0.; fsx2 = 0.; fsy2 = 0.; fsz2 = 0.
    DO g = 1, ng
       DO n = 1, nnod
          fs0(n) = fs0(n) + nuf(n,g) * f0(n,g)
          fsx1(n) = fsx1(n) + nuf(n,g) * fx1(n,g)
          fsy1(n) = fsy1(n) + nuf(n,g) * fy1(n,g)
          fsz1(n) = fsz1(n) + nuf(n,g) * fz1(n,g)
          fsx2(n) = fsx2(n) + nuf(n,g) * fx2(n,g)
          fsy2(n) = fsy2(n) + nuf(n,g) * fy2(n,g)
          fsz2(n) = fsz2(n) + nuf(n,g) * fz2(n,g)
       END DO
    END DO

    ! Update precursor density
    DO n = 1, nnod
       DO j = 1, nf
          lat = (1. + lamb(j) * tstep1)
          c0(j,n) = (ct(j,n) + iBeta(j) * tstep1 * fs0(n)) / lat
          cx1(j,n) = (ctx1(j,n) + iBeta(j) * tstep1 * fsx1(n)) / lat
          cy1(j,n) = (cty1(j,n) + iBeta(j) * tstep1 * fsy1(n)) / lat
          cz1(j,n) = (ctz1(j,n) + iBeta(j) * tstep1 * fsz1(n)) / lat
          cx2(j,n) = (ctx2(j,n) + iBeta(j) * tstep1 * fsx2(n)) / lat
          cy2(j,n) = (cty2(j,n) + iBeta(j) * tstep1 * fsy2(n)) / lat
          cz2(j,n) = (ctz2(j,n) + iBeta(j) * tstep1 * fsz2(n)) / lat
       END DO
    END DO

    ! Calculate power
    CALL powdis2(f0,tpow2)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, 0., tpow2/tpow1, (bpos(n), n = 1, nb)

END DO

! Second Time Step
! imax = CEILING((ttot-tdiv)/tstep2)
! stime = .FALSE.
!
! DO i = 1, imax
!
!
!
! END DO




END SUBROUTINE rod_eject2


SUBROUTINE rod_ejecte()

!
! Purpose:
!    To perform rod ejection simulation
!

USE sdata, ONLY: ng, nnod, sigr, nuf, sigs, f0, &
                 iBeta, lamb, nf, nout, nac, &
                 ttot, tdiv, tstep1, tstep2, &
                 ix, iy, iz, zdel, vdel, Ke, &
                 bcon, ftem, mtem, cden, &
                 fbpos, bpos, tmove, bspeed, mdir, &
                 nout, nac, nb, xnuf, dnuf, sigr, velo, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, &
                 fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2, &
                 c0, cx1, cy1, cz1, cx2, cy2, cz2, tbeta, &
                 ct, ctx1, cty1, ctz1, ctx2, cty2, ctz2, omeg
USE InpOutp, ONLY: XS_updt, bther, ounit
USE nodal, ONLY: nodal_coup4, outer4, outertf, powdis2

IMPLICIT NONE

REAL, DIMENSION(nnod, ng) :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Flux at previous time step

REAL :: t1, t2, lat
REAL :: tpow1, tpow2
INTEGER :: n, i, j, g, imax, step
LOGICAL :: stime

! Allocate fission source
ALLOCATE (fs0(nnod), fsx1(nnod), fsy1(nnod), fsz1(nnod))
ALLOCATE (fsx2(nnod), fsy2(nnod), fsz2(nnod))

! Allocate precusor density
ALLOCATE (c0(nf,nnod),cx1(nf,nnod),cy1(nf,nnod),cz1(nf,nnod))
ALLOCATE (cx2(nf,nnod),cy2(nf,nnod),cz2(nf,nnod))
ALLOCATE (ct(nf,nnod),ctx1(nf,nnod),cty1(nf,nnod),ctz1(nf,nnod))
ALLOCATE (ctx2(nf,nnod),cty2(nf,nnod),ctz2(nf,nnod))

! Allocate Frequency transformation constant
ALLOCATE (omeg(nnod,ng))

! Calculate forward flux at t=0 and check if keff=1
CALL XS_updt(bcon, ftem, mtem, cden, bpos)
CALL nodal_coup4()
CALL outer4(0)

! If K-EFF NOT EQUAL TO 1.0
IF (ABS(Ke - 1.) > 1.e-5) THEN
   WRITE(ounit, *)
   WRITE(ounit, '(A46,F9.6)') '  INITIAL MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke
   WRITE(ounit, *) '  WARNING: THE STEADY STATE K-EFF IS NOT EQUAL TO 1.0'
   WRITE(ounit, *) '  AND NOW IT IS FORCED TO 1.0 BY MODIFYING THE nu*sigf CROSS SECTIONS '
   WRITE(ounit, *)
   DO i = 1, 10
      xnuf = xnuf / Ke
      dnuf = dnuf / Ke
      CALL XS_updt(bcon, ftem, mtem, cden, bpos)
      CALL outer4(0)
      IF (ABS(Ke-1.0) < 1.e-5) EXIT
   END DO
   IF (i == 10) STOP "K-EFF STILL NOT EQUAL TO ONE. ADPRES IS STOPPING"
END IF

! ReCalculate forward flux at t=0 after keff fixed
CALL outer4(0)

! Calculate power
CALL powdis2(f0,tpow1)

! Calculate Initial precursor density
DO n = 1, nnod
   DO j = 1, nf
      c0(j,n) = iBeta(j) * fs0(n) / lamb(j)
      cx1(j,n) = iBeta(j) * fsx1(n) / lamb(j)
      cy1(j,n) = iBeta(j) * fsy1(n) / lamb(j)
      cz1(j,n) = iBeta(j) * fsz1(n) / lamb(j)
      cx2(j,n) = iBeta(j) * fsx2(n) / lamb(j)
      cy2(j,n) = iBeta(j) * fsy2(n) / lamb(j)
      cz2(j,n) = iBeta(j) * fsz2(n) / lamb(j)
   END DO
END DO

! Total beta
tbeta = 0.
DO j = 1, nf
  tbeta = tbeta + iBeta(j)
END DO

WRITE(ounit, *)
WRITE(ounit, *) " TRANSIENT RESULTS :"
WRITE(ounit, *)
WRITE(ounit, *) " Step  Time(s)  React.($)   Rel. Power   CR Bank Pos. (1-end)"
WRITE(ounit, *) "--------------------------------------------------------------"
WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') 0, 0., 0., 1.0, (bpos(n), n = 1, nb)

! Start transient calculation
step = 0
t2 = 0.
imax = CEILING(tdiv/tstep1)
stime = .FALSE.

! First Time Step
DO i = 1, imax

    step = step + 1
    t1 = t2
    t2 = t1 + tstep1

    IF (t2 > tdiv) THEN
        t2 = tdiv
        stime = .TRUE.
    END IF

    ! Rod bank changes
    DO n = 1, nb
        IF (mdir(n) == 1) THEN   ! If CR moving down
            IF (t2-tmove(n) > 1.d-5 .AND. fbpos(n)-bpos(n) < 1.d-5) THEN
                bpos(n) = bpos(n) - tstep1 *  bspeed(n)
                IF (bpos(n) < fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE IF (mdir(n) == 2) THEN ! If CR moving up
            IF (t2-tmove(n) > 1.d-5 .AND. fbpos(n)-bpos(n) > 1.d-5) THEN
                bpos(n) = bpos(n) + tstep1 *  bspeed(n)
                IF (bpos(n) > fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE
            CONTINUE
        END IF
     END DO

    ! Calculate xsec after pertubation
    CALL XS_updt(bcon, ftem, mtem, cden, bpos)

    ! Modify removal xsec
    IF (i > 1) THEN
      DO g = 1, ng
         DO n = 1, nnod
           omeg(n,g) = LOG(f0(n,g) / ft(n,g)) / tstep1
           sigr(n,g) = sigr(n,g) + 1. / (velo(g) * tstep1) + omeg(n,g) / velo(g)
         END DO
      END DO
    ELSE
       DO g = 1, ng
          DO n = 1, nnod
             omeg(n,g) = 0.
             sigr(n,g) = sigr(n,g) + 1. / (velo(g) * tstep1)
          END DO
       END DO
    END IF

    ! Save the previous flux and DN precusor density
    ft = f0
    ftx1 = fx1; fty1 = fy1; ftz1 = fz1
    ftx2 = fx2; fty2 = fy2; ftz2 = fz2
    ct = c0
    ctx1 = cx1; cty1 = cy1; ctz1 = cz1
    ctx2 = cx2; cty2 = cy2; ctz2 = cz2

    ! Transient calculation
    CALL nodal_coup4()
    CALL outertf(200, t2, tstep1, ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2)

    ! Update fission source
    fs0 = 0.; fsx1 = 0.; fsy1 = 0.; fsz1 = 0.; fsx2 = 0.; fsy2 = 0.; fsz2 = 0.
    DO g = 1, ng
       DO n = 1, nnod
          fs0(n) = fs0(n) + nuf(n,g) * f0(n,g)
          fsx1(n) = fsx1(n) + nuf(n,g) * fx1(n,g)
          fsy1(n) = fsy1(n) + nuf(n,g) * fy1(n,g)
          fsz1(n) = fsz1(n) + nuf(n,g) * fz1(n,g)
          fsx2(n) = fsx2(n) + nuf(n,g) * fx2(n,g)
          fsy2(n) = fsy2(n) + nuf(n,g) * fy2(n,g)
          fsz2(n) = fsz2(n) + nuf(n,g) * fz2(n,g)
       END DO
    END DO

    ! Update precursor density
    DO n = 1, nnod
       DO j = 1, nf
          lat = (1. + lamb(j) * tstep1)
          c0(j,n) = (ct(j,n) + iBeta(j) * tstep1 * fs0(n)) / lat
          cx1(j,n) = (ctx1(j,n) + iBeta(j) * tstep1 * fsx1(n)) / lat
          cy1(j,n) = (cty1(j,n) + iBeta(j) * tstep1 * fsy1(n)) / lat
          cz1(j,n) = (ctz1(j,n) + iBeta(j) * tstep1 * fsz1(n)) / lat
          cx2(j,n) = (ctx2(j,n) + iBeta(j) * tstep1 * fsx2(n)) / lat
          cy2(j,n) = (cty2(j,n) + iBeta(j) * tstep1 * fsy2(n)) / lat
          cz2(j,n) = (ctz2(j,n) + iBeta(j) * tstep1 * fsz2(n)) / lat
       END DO
    END DO

    ! Calculate power
    CALL powdis2(f0, tpow2)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, 0., tpow2/tpow1, (bpos(n), n = 1, nb)
    ! CALL powdis2(f0, tpow2)
    !
    ! WRITE(*,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, 0., tpow2/tpow1, (bpos(n), n = 1, nb)

END DO

! Second Time Step
! imax = CEILING((ttot-tdiv)/tstep2)
! stime = .FALSE.
!
! DO i = 1, imax
!
!
!
! END DO




END SUBROUTINE rod_ejecte


SUBROUTINE trod_ejecte()

!
! Purpose:
!    To perform rod ejection simulation
!

USE sdata, ONLY: ng, nnod, sigr, nuf, sigs, f0, &
                 iBeta, lamb, nf, nout, nac, &
                 ttot, tdiv, tstep1, tstep2, &
                 ix, iy, iz, zdel, vdel, Ke, &
                 bcon, ftem, mtem, cden, &
                 fbpos, bpos, tmove, bspeed, mdir, &
                 nout, nac, nb, xnuf, dnuf, sigr, velo, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, &
                 fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2, &
                 c0, cx1, cy1, cz1, cx2, cy2, cz2, tbeta, &
                 ct, ctx1, cty1, ctz1, ctx2, cty2, ctz2, omeg
USE InpOutp, ONLY: XS_updt, bther, ounit
USE nodal, ONLY: nodal_coup4, outer4, outertf, powdis2

IMPLICIT NONE

REAL, DIMENSION(nnod, ng) :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Flux at previous time step

REAL :: t1, t2, lat
REAL :: tpow1, tpow2
INTEGER :: n, i, j, g, imax, step
LOGICAL :: stime

! Allocate fission source
ALLOCATE (fs0(nnod), fsx1(nnod), fsy1(nnod), fsz1(nnod))
ALLOCATE (fsx2(nnod), fsy2(nnod), fsz2(nnod))

! Allocate precusor density
ALLOCATE (c0(nf,nnod),cx1(nf,nnod),cy1(nf,nnod),cz1(nf,nnod))
ALLOCATE (cx2(nf,nnod),cy2(nf,nnod),cz2(nf,nnod))
ALLOCATE (ct(nf,nnod),ctx1(nf,nnod),cty1(nf,nnod),ctz1(nf,nnod))
ALLOCATE (ctx2(nf,nnod),cty2(nf,nnod),ctz2(nf,nnod))

! Allocate Frequency transformation constant
ALLOCATE (omeg(nnod,ng))

! Calculate forward flux at t=0 and check if keff=1
CALL XS_updt(bcon, ftem, mtem, cden, bpos)
CALL nodal_coup4()
CALL outer4(0)

! If K-EFF NOT EQUAL TO 1.0
IF (ABS(Ke - 1.) > 1.e-5) THEN
   WRITE(ounit, *)
   WRITE(ounit, '(A46,F9.6)') '  INITIAL MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke
   WRITE(ounit, *) '  WARNING: THE STEADY STATE K-EFF IS NOT EQUAL TO 1.0'
   WRITE(ounit, *) '  AND NOW IT IS FORCED TO 1.0 BY MODIFYING THE nu*sigf CROSS SECTIONS '
   WRITE(ounit, *)
   DO i = 1, 10
      xnuf = xnuf / Ke
      dnuf = dnuf / Ke
      CALL XS_updt(bcon, ftem, mtem, cden, bpos)
      CALL outer4(0)
      IF (ABS(Ke-1.0) < 1.e-5) EXIT
   END DO
   IF (i == 10) STOP "K-EFF STILL NOT EQUAL TO ONE. ADPRES IS STOPPING"
END IF

! ReCalculate forward flux at t=0 after keff fixed
CALL outer4(0)

! Calculate power
CALL powdis2(f0,tpow1)

! Calculate Initial precursor density
DO n = 1, nnod
   DO j = 1, nf
      c0(j,n) = iBeta(j) * fs0(n) / lamb(j)
      cx1(j,n) = iBeta(j) * fsx1(n) / lamb(j)
      cy1(j,n) = iBeta(j) * fsy1(n) / lamb(j)
      cz1(j,n) = iBeta(j) * fsz1(n) / lamb(j)
      cx2(j,n) = iBeta(j) * fsx2(n) / lamb(j)
      cy2(j,n) = iBeta(j) * fsy2(n) / lamb(j)
      cz2(j,n) = iBeta(j) * fsz2(n) / lamb(j)
   END DO
END DO

! Total beta
tbeta = 0.
DO j = 1, nf
  tbeta = tbeta + iBeta(j)
END DO

WRITE(ounit, *)
WRITE(ounit, *) " TRANSIENT RESULTS :"
WRITE(ounit, *)
WRITE(ounit, *) " Step  Time(s)  React.($)   Rel. Power   CR Bank Pos. (1-end)"
WRITE(ounit, *) "--------------------------------------------------------------"
WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') 0, 0., 0., 1.0, (bpos(n), n = 1, nb)

! Start transient calculation
step = 0
t2 = 0.
imax = CEILING(tdiv/tstep1)
stime = .FALSE.

! First Time Step
DO i = 1, imax

    step = step + 1
    t1 = t2
    t2 = t1 + tstep1

    IF (t2 > tdiv) THEN
        t2 = tdiv
        stime = .TRUE.
    END IF

    ! Rod bank changes
    DO n = 1, nb
        IF (mdir(n) == 1) THEN   ! If CR moving down
            IF (t2-tmove(n) > 1.d-5 .AND. fbpos(n)-bpos(n) < 1.d-5) THEN
                bpos(n) = bpos(n) - tstep1 *  bspeed(n)
                IF (bpos(n) < fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE IF (mdir(n) == 2) THEN ! If CR moving up
            IF (t2-tmove(n) > 1.d-5 .AND. fbpos(n)-bpos(n) > 1.d-5) THEN
                bpos(n) = bpos(n) + tstep1 *  bspeed(n)
                IF (bpos(n) > fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE
            CONTINUE
        END IF
     END DO

    ! Calculate xsec after pertubation
    CALL XS_updt(bcon, ftem, mtem, cden, bpos)

    ! Modify removal xsec
    IF (i > 1) THEN
      DO g = 1, ng
         DO n = 1, nnod
           omeg(n,g) = LOG(f0(n,g) / ft(n,g)) / tstep1
           sigr(n,g) = sigr(n,g) + 1. / (velo(g) * tstep1) + omeg(n,g) / velo(g)
         END DO
      END DO
    ELSE
       DO g = 1, ng
          DO n = 1, nnod
             omeg(n,g) = 0.
             sigr(n,g) = sigr(n,g) + 1. / (velo(g) * tstep1)
          END DO
       END DO
    END IF

    ! Save the previous flux and DN precusor density
    ft = f0
    ftx1 = fx1; fty1 = fy1; ftz1 = fz1
    ftx2 = fx2; fty2 = fy2; ftz2 = fz2
    ct = c0
    ctx1 = cx1; cty1 = cy1; ctz1 = cz1
    ctx2 = cx2; cty2 = cy2; ctz2 = cz2

    ! Transient calculation
    CALL nodal_coup4()
    CALL outertf(200, t2, tstep1, ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2)

    ! Update fission source
    fs0 = 0.; fsx1 = 0.; fsy1 = 0.; fsz1 = 0.; fsx2 = 0.; fsy2 = 0.; fsz2 = 0.
    DO g = 1, ng
       DO n = 1, nnod
          fs0(n) = fs0(n) + nuf(n,g) * f0(n,g)
          fsx1(n) = fsx1(n) + nuf(n,g) * fx1(n,g)
          fsy1(n) = fsy1(n) + nuf(n,g) * fy1(n,g)
          fsz1(n) = fsz1(n) + nuf(n,g) * fz1(n,g)
          fsx2(n) = fsx2(n) + nuf(n,g) * fx2(n,g)
          fsy2(n) = fsy2(n) + nuf(n,g) * fy2(n,g)
          fsz2(n) = fsz2(n) + nuf(n,g) * fz2(n,g)
       END DO
    END DO

    ! Update precursor density
    DO n = 1, nnod
       DO j = 1, nf
          lat = (1. + lamb(j) * tstep1)
          c0(j,n) = (ct(j,n) + iBeta(j) * tstep1 * fs0(n)) / lat
          cx1(j,n) = (ctx1(j,n) + iBeta(j) * tstep1 * fsx1(n)) / lat
          cy1(j,n) = (cty1(j,n) + iBeta(j) * tstep1 * fsy1(n)) / lat
          cz1(j,n) = (ctz1(j,n) + iBeta(j) * tstep1 * fsz1(n)) / lat
          cx2(j,n) = (ctx2(j,n) + iBeta(j) * tstep1 * fsx2(n)) / lat
          cy2(j,n) = (cty2(j,n) + iBeta(j) * tstep1 * fsy2(n)) / lat
          cz2(j,n) = (ctz2(j,n) + iBeta(j) * tstep1 * fsz2(n)) / lat
       END DO
    END DO

    ! Calculate power
    CALL powdis2(f0, tpow2)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, 0., tpow2/tpow1, (bpos(n), n = 1, nb)
    ! CALL powdis2(f0, tpow2)
    !
    ! WRITE(*,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, 0., tpow2/tpow1, (bpos(n), n = 1, nb)

END DO

! Second Time Step
! imax = CEILING((ttot-tdiv)/tstep2)
! stime = .FALSE.
!
! DO i = 1, imax
!
!
!
! END DO




END SUBROUTINE trod_ejecte

END MODULE trans
