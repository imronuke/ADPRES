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
REAL :: tbeta
LOGICAL :: logi = .TRUE.

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

USE sdata, ONLY: ng, nnod, nuf, chi, velo, nf, iBeta, sigs
USE nodal, ONLY: Integrate

IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(IN) :: sf, dsigr, dnuf
REAL, DIMENSION(:,:,:), INTENT(IN) :: dsigs
REAL, INTENT(OUT) :: xA, xrho

INTEGER :: n, i, g, h
REAL, DIMENSION(nnod) :: vdum, vdum2, vdum3, vdum4, vdum5
REAL, DIMENSION(nnod,ng) :: yvdum
REAL :: F2, rho2

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
DO i = 1, nf
    vdum2 = 0.
    DO g = 1, ng
        DO n = 1, nnod
            vdum2(n) = vdum2(n) + chi(n,g) * iBeta(i) * vdum(n) * af(n,g)
        END DO
    END DO
    beta(i) = Integrate(vdum2) / F2
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
                 ix, iy, iz, zdel,Ke, &
                 bcon, ftem, mtem, cden, &
                 fbpos, bpos, tmove, bspeed, mdir, &
                 nout, nac, nb
USE InpOutp, ONLY: XS_updt, bther, ounit
USE nodal, ONLY: nodal_coup4, outer4

IMPLICIT NONE

REAL, DIMENSION(nnod,ng) ::  dsigr, dnuf
REAL, DIMENSION(nnod,ng,ng) :: dsigs

REAL, DIMENSION(nnod,ng) ::  osigr, onuf
REAL, DIMENSION(nnod,ng,ng) :: osigs

REAL :: A, rho

REAL :: t1, t2
REAL :: xppow
INTEGER :: n, i, j, imax, step
REAL, PARAMETER :: hp = 0.0001 ! Point Kinetetic Time step
LOGICAL :: stime


! Calculate forward flux at t=0
CALL XS_updt(bcon, ftem, mtem, cden, bpos)
CALL adj_calc()
CALL nodal_coup4()
CALL outer4(0)


! Save old sigr, nuf and sigs
osigr = sigr; onuf = nuf; osigs = sigs
dsigr = 0.; dnuf = 0.; dsigs = 0.

! Calculate intgral kinet parameters at t = 0
CALL kinet_par(dsigr, dnuf, dsigs, f0, A, rho)

! Calculate Initial precursor density
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

    ! Calculate shape function
    CALL nodal_coup4()
    CALL outer4(0)

    ! Calculate xsec changes after rod is ejected
    dsigr = sigr - osigr
    dnuf = nuf - onuf
    dsigs = sigs - osigs

    ! Calculate intgral kinet parameters
    CALL kinet_par(dsigr, dnuf, dsigs, f0, A, rho)

    tbeta = 0.
    DO j = 1, nf
        tbeta = tbeta + beta(j)
    END DO

    !Calculate amplitude function
    CALL point(t1, t2, hp, rho, A, beta, amp)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho/tbeta, amp, (bpos(n), n = 1, nb)

    IF (stime) EXIT

    IF (step>1000) THEN
        WRITE(ounit,*) 'TOO SMALL TIME STEPS. STOPPING'
        STOP
    END IF

END DO


! Second Time Step
imax = CEILING((ttot-tdiv)/tstep2)
stime = .FALSE.

DO i = 1, imax

    step = step + 1
    t1 = t2
    t2 = t1 + tstep2

    IF (t2 > ttot) THEN
        t2 = ttot
        stime = .TRUE.
    END IF

    ! Rod bank changes
    DO n = 1, nb
        IF (mdir(n) == 1) THEN   ! If CR moving down
            IF (t2-tmove(n) > 1.d-5 .AND. fbpos(n)-bpos(n) < 1.d-5) THEN
                bpos(n) = bpos(n) - tstep2 *  bspeed(n)
                IF (bpos(n) < fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE IF (mdir(n) == 2) THEN ! If CR moving up
            IF (t2-tmove(n) > 1.d-5 .AND. fbpos(n)-bpos(n) > 1.d-5) THEN
                bpos(n) = bpos(n) + tstep2 *  bspeed(n)
                IF (bpos(n) > fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE
            CONTINUE
        END IF
     END DO

    ! Calculate xsec after pertubation
    CALL XS_updt(bcon, ftem, mtem, cden, bpos)

    ! Calculate shape function
    CALL nodal_coup4()
    CALL outer4(0)

    ! Calculate xsec changes after rod is ejected
    dsigr = sigr - osigr
    dnuf = nuf - onuf
    dsigs = sigs - osigs

    ! Calculate intgral kinet parameters
    CALL kinet_par(dsigr, dnuf, dsigs, f0, A, rho)

    tbeta = 0.
    DO j = 1, nf
        tbeta = tbeta + beta(j)
    END DO

    !Calculate amplitude function
    CALL point(t1, t2, hp, rho, A, beta, amp)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho/tbeta, amp, (bpos(n), n = 1, nb)

    IF (stime) EXIT

    IF (step>1000) THEN
        WRITE(ounit,*) 'TOO SMALL TIME STEPS. STOPPING'
        STOP
    END IF

END DO

END SUBROUTINE rod_eject


SUBROUTINE trod_eject()

!
! Purpose:
!    To perform rod ejection simulation
!

USE sdata, ONLY: ng, nnod, sigr, nuf, sigs, f0, &
                 iBeta, lamb, nf, nout, nac, &
                 ttot, tdiv, tstep1, tstep2, &
                 ix, iy, iz, zdel,Ke, &
                 bcon, ftem, mtem, cden, &
                 fbpos, bpos, tmove, bspeed, mdir, &
                 nout, nac, nb, npow, ppow, pow, node_nf
USE InpOutp, ONLY: XS_updt, bther, ounit
USE nodal, ONLY: nodal_coup4, outer4, powdis
USE th, ONLY: th_iter, th_trans3
!, par_ave_f, par_max, par_ave

IMPLICIT NONE

REAL, DIMENSION(nnod,ng) ::  dsigr, dnuf
REAL, DIMENSION(nnod,ng,ng) :: dsigs

REAL, DIMENSION(nnod,ng) ::  osigr, onuf
REAL, DIMENSION(nnod,ng,ng) :: osigs

REAL, DIMENSION(nnod) :: pline       ! Linear power density
REAL, DIMENSION(nnod,ng) :: fl

REAL :: A, rho

REAL :: t1, t2
REAL :: xppow
INTEGER :: n, i, j, imax, step
REAL, PARAMETER :: hp = 0.0001 ! Point Kinetetic Time step
LOGICAL :: stime


! Allocate node power distribution npow
ALLOCATE(npow(nnod))

! Calculate forward flux at t=0
CALL th_iter(0)
fl = f0

!Initial amplitude function
amp = 1.

! Save old sigr, nuf and sigs
osigr = sigr; onuf = nuf; osigs = sigs
dsigr = 0.; dnuf = 0.; dsigs = 0.

WRITE(ounit, *)
WRITE(ounit, *) " TRANSIENT RESULTS :"
WRITE(ounit, *)
WRITE(ounit, *) " Step  Time(s)  React.($)    Power (W)    CR Bank Pos. (1-end)"
WRITE(ounit, *) "--------------------------------------------------------------"
WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') 0, 0., 0., amp*ppow*0.01, (bpos(n), n = 1, nb)

! Adjoint flux at t=0
CALL adj_calc()

! Calculate intgral kinet parameters at t = 0
CALL kinet_par(dsigr, dnuf, dsigs, fl, A, rho)

! Calculate Initial precursor density
DO j = 1, nf
    C(j) = beta(j) / (A * lamb(j))   ! See Eq. 6-32 D&H
END DO

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

    ! Calculate shape function
    CALL nodal_coup4()
    CALL outer4(0)

    ! Calculate xsec changes after rod is ejected
    dsigr = sigr - osigr
    dnuf = nuf - onuf
    dsigs = sigs - osigs

    ! Calculate intgral kinet parameters
    CALL kinet_par(dsigr, dnuf, dsigs, f0, A, rho)

    tbeta = 0.
    DO j = 1, nf
        tbeta = tbeta + beta(j)
    END DO

    !Calculate amplitude function
    CALL point(t1, t2, hp, rho, A, beta, amp)

    ! Calculate node power distribution
    CALL powdis(npow)

    ! Power change
    xppow = ppow * amp * 0.01

    ! Calculate linear power density for each nodes (W/cm)
    DO n = 1, nnod
        pline(n) = npow(n) * pow * xppow &
                 / (node_nf(ix(n),iy(n)) * zdel(iz(n)))     ! Linear power density (W/cm)
    END DO

    ! TH transient
    CALL th_trans3(pline,tstep1)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho/tbeta, xppow, (bpos(n), n = 1, nb)

    IF (stime) EXIT

    IF (step>1000) THEN
        WRITE(ounit,*) 'TOO SMALL TIME STEPS. STOPPING'
        STOP
    END IF

END DO


! Second Time Step
imax = CEILING((ttot-tdiv)/tstep2)
stime = .FALSE.

DO i = 1, imax

    step = step + 1
    t1 = t2
    t2 = t1 + tstep2

    IF (t2 > ttot) THEN
        t2 = ttot
        stime = .TRUE.
    END IF

    ! Rod bank changes
    DO n = 1, nb
        IF (mdir(n) == 1) THEN   ! If CR moving down
            IF (t2-tmove(n) > 1.d-5 .AND. fbpos(n)-bpos(n) < 1.d-5) THEN
                bpos(n) = bpos(n) - tstep2 *  bspeed(n)
                IF (bpos(n) < fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE IF (mdir(n) == 2) THEN ! If CR moving up
            IF (t2-tmove(n) > 1.d-5 .AND. fbpos(n)-bpos(n) > 1.d-5) THEN
                bpos(n) = bpos(n) + tstep2 *  bspeed(n)
                IF (bpos(n) > fbpos(n)) bpos(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE
            CONTINUE
        END IF
     END DO

    ! Calculate xsec after pertubation
    CALL XS_updt(bcon, ftem, mtem, cden, bpos)

    ! Calculate shape function
    CALL nodal_coup4()
    CALL outer4(0)

    ! Calculate xsec changes after rod is ejected
    dsigr = sigr - osigr
    dnuf = nuf - onuf
    dsigs = sigs - osigs

    ! Calculate intgral kinet parameters
    CALL kinet_par(dsigr, dnuf, dsigs, f0, A, rho)

    tbeta = 0.
    DO j = 1, nf
        tbeta = tbeta + beta(j)
    END DO

    !Calculate amplitude function
    CALL point(t1, t2, hp, rho, A, beta, amp)

    ! Calculate node power distribution
    CALL powdis(npow)

    ! Power change
    xppow = ppow * amp * 0.01

    ! Calculate linear power density for each nodes (W/cm)
    DO n = 1, nnod
        pline(n) = npow(n) * pow * xppow  &
                 / (node_nf(ix(n),iy(n)) * zdel(iz(n)))     ! Linear power density (W/cm)
    END DO

    ! TH transient
    CALL th_trans3(pline,tstep2)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho/tbeta, xppow, (bpos(n), n = 1, nb)

    IF (stime) EXIT

    IF (step>1000) THEN
        WRITE(ounit,*) 'TOO SMALL TIME STEPS. STOPPING'
        STOP
    END IF

END DO

END SUBROUTINE trod_eject



SUBROUTINE point(ti, tf, h, xrho, xA, xbet, xamp)

!
! Purpose:
!    To perform transient point reactor calculation using RK-4
!

USE sdata, ONLY: nf, lamb

IMPLICIT NONE

REAL, INTENT(IN) :: ti, tf, h              ! Initial time, final time, time increment
REAL, INTENT(IN) :: xrho, xA               ! reactivity and neutron generation time
REAL, DIMENSION(:), INTENT(IN) :: xbet     ! delayed neutron fraction
REAL, INTENT(INOUT) :: xamp                ! amplitude function

INTEGER :: j, i, itot
REAL :: summ
REAL :: k1, k2, k3, k4
REAL :: xtbet


itot = INT((tf - ti) / h)
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
