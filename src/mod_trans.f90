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
REAL, DIMENSION(:), ALLOCATABLE :: beta, C   ! beta, neutron precusor density
REAL, DIMENSION(:,:), ALLOCATABLE :: xvdum
REAL :: amp                             ! Amplitude function

CONTAINS

SUBROUTINE adj_calc()

!
! Purpose:
!    To calculate adjoint flux
!

USE sdata, ONLY: f0, ng, nnod, nf, bcon, ftem, mtem, cden, bpos
USE nodal, ONLY: nodal_coup4, outer4ad
USE InpOutp, ONLY: XS_updt

IMPLICIT NONE

! Start adjoint flux calculation
CALL XS_updt(bcon, ftem, mtem, cden, bpos)
CALL nodal_coup4()
CALL outer4ad(0)

ALLOCATE(af(nnod,ng))
af = f0       ! Save adjoint flux to af

ALLOCATE(beta(nf), C(nf))
ALLOCATE(xvdum(nnod,ng))

END SUBROUTINE adj_calc



SUBROUTINE kinet_par(dsigr, dnuf, dsigs, sf, xA, xrho)

!
! Purpose:
!    To calculate kinetic parameters
!

USE sdata, ONLY: ng, nnod, chi, velo, nf, iBeta, tbeta, nuf, &
vdel, nod, xdel, ydel, zdel, ix, iy, iz
USE nodal, ONLY: Integrate

IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(IN) :: sf, dsigr, dnuf
REAL, DIMENSION(:,:,:), INTENT(IN) :: dsigs
REAL, INTENT(OUT) :: xA, xrho

INTEGER :: n, i, g, h
REAL, DIMENSION(nnod) :: vdum, vdum2, vdum3, vdum4, vdum5
REAL, DIMENSION(nnod,ng) :: yvdum
REAL :: F2, F, fiss, scat, leak, remo

! Calculate F
vdum = 0.
DO g = 1, ng
    DO n = 1, nnod
        vdum(n) = vdum(n) + nuf(n,g) * sf(n,g)
    END DO
END DO

vdum2 = 0.
DO g = 1, ng
    DO n = 1, nnod
        vdum2(n) = vdum2(n) + chi(n,g) * vdum(n) * af(n,g)
    END DO
END DO

F2 = Integrate(vdum2)


! Calculate neutron gneration time (A)
vdum2 = 0.
DO g = 1, ng
    DO n = 1, nnod
        vdum2(n) = vdum2(n) + af(n,g) * sf(n,g) / velo(g)
    END DO
END DO

xA = Integrate(vdum2) / F2    ! Calculate neutron generation time


! Calculate Delayed neutron fraction (beta)
tbeta = 0.
DO i = 1, nf
    vdum2 = 0.
    DO g = 1, ng
        DO n = 1, nnod
            vdum2(n) = vdum2(n) + chi(n,g) * iBeta(i) * vdum(n) * af(n,g)
        END DO
    END DO
    beta(i) = Integrate(vdum2) / F2
    tbeta = tbeta + beta(i)
END DO


! Calculate reactivity (rho)
vdum2 = 0.; vdum3 = 0.; vdum4 = 0.
DO g = 1, ng
    DO n = 1, nnod
        vdum2(n) = vdum2(n) + dnuf(n,g) * sf(n,g)
        vdum3(n) = vdum3(n) + af(n,g) * dsigr(n,g) * sf(n,g)
    END DO
END DO

DO g = 1, ng
      vdum5 = 0.
    DO h = 1, ng
        DO n = 1, nnod
            vdum5(n) = vdum5(n) + dsigs(n,h,g) * sf(n,h)
        END DO
    END DO
    DO n = 1, nnod
        vdum4(n) = vdum4(n) + af(n,g) * (vdum5(n) + chi(n,g) * vdum2(n))
    END DO
END DO

xrho = Integrate(vdum4 - vdum3) / F2

END SUBROUTINE kinet_par


SUBROUTINE rod_eject()

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
                 nout, nac, nb, xnuf, dnuf, velo, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, &
                 fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2, &
                 c0, cx1, cy1, cz1, cx2, cy2, cz2, tbeta
USE InpOutp, ONLY: XS_updt, bther, ounit
USE nodal, ONLY: nodal_coup4, outer4, outertr, powdis2

IMPLICIT NONE

REAL :: A, rho

REAL, DIMENSION(nnod, ng) :: fh  ! Shape function
REAL, DIMENSION(nnod, ng) :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Flux at previous time step
REAL, DIMENSION(6,nnod) :: ct, ctx1, cty1, ctz1, ctx2, cty2, ctz2
REAL :: gamma1, gamma2

REAL, DIMENSION(nnod,ng) ::  ksigr, knuf
REAL, DIMENSION(nnod,ng,ng) :: ksigs

REAL, DIMENSION(nnod,ng) ::  osigr, onuf
REAL, DIMENSION(nnod,ng,ng) :: osigs

REAL :: t1, t2
REAL :: tpow1, tpow2
INTEGER :: n, i, j, g, imax, step
REAL, PARAMETER :: hp = 0.0001 ! Point Kinetetic Time step
LOGICAL :: stime

! Allocate fission source
ALLOCATE (fs0(nnod), fsx1(nnod), fsy1(nnod), fsz1(nnod))
ALLOCATE (fsx2(nnod), fsy2(nnod), fsz2(nnod))

! Allocate precusor density
ALLOCATE (c0(nf,nnod),cx1(nf,nnod),cy1(nf,nnod),cz1(nf,nnod))
ALLOCATE (cx2(nf,nnod),cy2(nf,nnod),cz2(nf,nnod))

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
      CALL outer4(0)
      IF (ABS(Ke-1.0) < 1.e-5) EXIT
      xnuf = xnuf / Ke
      dnuf = dnuf / Ke
      CALL XS_updt(bcon, ftem, mtem, cden, bpos)
   END DO
   IF (i == 10) STOP "K-EFF STILL NOT EQUAL TO ONE. ADPRES IS STOPPING"
END IF

! Calculate adjoint flux at t=0
CALL adj_calc()

! ReCalculate forward flux at t=0 after keff fixed
CALL outer4(0)

! Calculate power
CALL powdis2(tpow1)

! Calculate Initial precursor density for shape function
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

! Calculate gamma
gamma1 = 0.
DO g = 1, ng
    DO n = 1, nnod
        gamma1 = gamma1 + af(n,g) * f0(n,g) / velo(g) * vdel(n)
    END DO
END DO

! Save old sigr, nuf and sigs
osigr = sigr; onuf = nuf; osigs = sigs
ksigr = 0.; knuf = 0.; ksigs = 0.

! Calculate intgral kinet parameters at t = 0
CALL kinet_par(ksigr, knuf, ksigs, f0, A, rho)

! Calculate Initial precursor density for point kinetic
amp = 1.
DO j = 1, nf
    C(j) = beta(j) / (A * lamb(j))   ! See Eq. 6-32 D&H
END DO

WRITE(ounit, *)
WRITE(ounit, *) " TRANSIENT RESULTS :"
WRITE(ounit, *)
WRITE(ounit, *) " Step  Time(s)  React.($)   Rel. Power   CR Bank Pos. (1-end)"
WRITE(ounit, *) "--------------------------------------------------------------"
WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') 0, 0., 0., amp, (bpos(n), n = 1, nb)

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

    ! Predict the flux
    ft = f0
    ftx1 = fx1; fty1 = fy1; ftz1 = fz1
    ftx2 = fx2; fty2 = fy2; ftz2 = fz2
    ct = c0
    ctx1 = cx1; cty1 = cy1; ctz1 = cz1
    ctx2 = cx2; cty2 = cy2; ctz2 = cz2
    ! CALL nodal_couptr(tstep1)
    ! CALL outertr(200, tstep1, ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2)
    CALL nodal_coup4()
    CALL outer4(0)

    ! Predict precusor density
    ! CALL precusor(tstep1, ct, ctx1, cty1, ctz1, ctx2, cty2, ctz2)

    ! Calculate the shape function
    gamma2 = 0.
    DO g = 1, ng
        DO n = 1, nnod
            gamma2 = gamma2 + af(n,g) * f0(n,g) / velo(g) * vdel(n)
        END DO
    END DO
    gamma2 = gamma2 / gamma1

    ft = f0 / gamma2  ! Note: ft used to store shape function to minimize memory usage
    ftx1 = fx1 / gamma2;  fty1 = fy1 / gamma2;  ftz1 = fz1 / gamma2
    ftx2 = fx2 / gamma2;  fty2 = fy2 / gamma2;  ftz2 = fz2 / gamma2

    ! Calculate xsec changes after rod is ejected
    ksigr = sigr - osigr
    knuf = nuf - onuf
    ksigs = sigs - osigs

    ! Calculate intgral kinet parameters
    ! CALL kinet_par(ksigr, knuf, ksigs, ft, A, rho)
    CALL kinet_par(ksigr, knuf, ksigs, f0, A, rho)

    ! Starts pint kinetic calculation
    CALL point(tstep1, hp, rho, A, beta, amp)

    IF (step>1000) THEN
        WRITE(ounit,*) 'TOO SMALL TIME STEPS. STOPPING'
        STOP
    END IF

    ! Correct the flux
    f0 = ft * amp
    fx1 = ftx1 * amp; fy1 = fty1 * amp; fz1 = ftz1 * amp
    fx2 = ftx2 * amp; fy2 = fty2 * amp; fz2 = ftz2 * amp

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

    ! Correct precusor density
    ! CALL precusor(tstep1, ct, ctx1, cty1, ctz1, ctx2, cty2, ctz2)

    ! Calculate power
    CALL powdis2(tpow2)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho/tbeta, amp, (bpos(n), n = 1, nb)

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




END SUBROUTINE rod_eject


SUBROUTINE precusor(h, ct, ctx1, cty1, ctz1, ctx2, cty2, ctz2)

!
! Purpose:
!    To update neutron precusor density
!

USE sdata, ONLY: ng, nnod, nf, iBeta, lamb, &
fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2, c0, cx1, cy1, cz1, cx2, cy2, cz2

IMPLICIT NONE

REAL, INTENT(IN) :: h
REAL, DIMENSION(:,:), INTENT(IN) :: ct, ctx1, cty1, ctz1, ctx2, cty2, ctz2

INTEGER :: n, j

DO n = 1, nnod
   DO j = 1, nf
      c0(j,n) = (ct(j,n) + iBeta(j) * h * fs0(n)) / (1. + h * lamb(j))
      cx1(j,n) = (ctx1(j,n) + iBeta(j) * h * fsx1(n)) / (1. + h * lamb(j))
      cy1(j,n) = (cty1(j,n) + iBeta(j) * h * fsy1(n)) / (1. + h * lamb(j))
      cz1(j,n) = (ctz1(j,n) + iBeta(j) * h * fsz1(n)) / (1. + h * lamb(j))
      cx2(j,n) = (ctx2(j,n) + iBeta(j) * h * fsx2(n)) / (1. + h * lamb(j))
      cy2(j,n) = (cty2(j,n) + iBeta(j) * h * fsy2(n)) / (1. + h * lamb(j))
      cz2(j,n) = (ctz2(j,n) + iBeta(j) * h * fsz2(n)) / (1. + h * lamb(j))
   END DO
END DO

END SUBROUTINE precusor


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

REAL :: A, rho

REAL, DIMENSION(nnod, ng) :: fh  ! Shape function
REAL, DIMENSION(nnod, ng) :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Flux at previous time step

REAL :: t1, t2
REAL :: tpow1, tpow2
INTEGER :: n, i, j, g, imax, step
REAL, PARAMETER :: hp = 0.0001 ! Point Kinetetic Time step
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


! Calculate power
CALL powdis2(tpow1)

! Calculate Initial precursor density for shape function
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

    ! Predict the flux
    ft = f0
    ftx1 = fx1; fty1 = fy1; ftz1 = fz1
    ftx2 = fx2; fty2 = fy2; ftz2 = fz2
    ct = c0
    ctx1 = cx1; cty1 = cy1; ctz1 = cz1
    ctx2 = cx2; cty2 = cy2; ctz2 = cz2
    DO g = 1, ng
       DO n = 1, nnod
          sigr(n,g) = sigr(n,g) + 1. / (velo(g) * tstep1)
       END DO
    END DO
    CALL nodal_coup4()
    CALL outertr(200, tstep1, ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2)

    ! Predict precusor density
    CALL precusor(tstep1, ct, ctx1, cty1, ctz1, ctx2, cty2, ctz2)

    ! Calculate power
    CALL powdis2(tpow2)

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


SUBROUTINE point(ht, h, xrho, xA, xbet, xamp)

!
! Purpose:
!    To perform transient point reactor calculation using RK-4
!

USE sdata, ONLY: nf, lamb, tbeta

IMPLICIT NONE

REAL, INTENT(IN) :: ht, h                  ! Macro time step, micro time step
REAL, INTENT(IN) :: xrho, xA               ! reactivity and neutron generation time
REAL, DIMENSION(:), INTENT(IN) :: xbet     ! delayed neutron fraction
REAL, INTENT(INOUT) :: xamp                ! amplitude function

INTEGER :: j, i, itot
REAL :: summ
REAL :: k1, k2, k3, k4
REAL :: xtbet

itot = INT(ht / h)
DO i = 1, itot

    !!!Calculate Power for time step = i
    summ = 0.
    DO j = 1, nf
        summ = summ + lamb(j)*C(j)
    END DO

    k1 = (xrho - tbeta) * xamp / xA + summ
    k2 = (xrho - tbeta) * (xamp + k1 * 0.5 * h) / xA + summ
    k3 = (xrho - tbeta) * (xamp + k2 * 0.5 * h) / xA + summ
    k4 = (xrho - tbeta) * (xamp + k3 * h) / xA + summ

    xamp = xamp + h/6. * (k1 + 2.* k2 + 2.* k3 + k4)

    !!!Calculate precursor density for each group j for time step = i
    DO j = 1, nf
        k1 = xbet(j) * xamp / xA - lamb(j) * C(j)
        k2 = xbet(j) * xamp / xA - lamb(j) * (C(j) + k1 * 0.5 * h)
        k3 = xbet(j) * xamp / xA - lamb(j) * (C(j) + k2 * 0.5 * h)
        k4 = xbet(j) * xamp / xA - lamb(j) * (C(j) + k3 * h)

        C(j) =  C(j) + h/6. * (k1 + 2.* k2 + 2.* k3 + k4)
    END DO
END DO


END SUBROUTINE point





END MODULE trans
