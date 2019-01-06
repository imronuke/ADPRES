MODULE trans

!=========================
! Transient Module to solve transient diffusion problems
! Using Fully Implicit method with exponetial transformation
! =======================


IMPLICIT NONE

SAVE

CONTAINS


SUBROUTINE rod_eject()

!
! Purpose:
!    To perform rod ejection simulation
!

USE sdata, ONLY: ng, nnod, sigr, nf, nout, &
                 ttot, tdiv, tstep1, tstep2, Ke, &
                 bcon, ftem, mtem, cden, &
                 fbpos, bpos, tmove, bspeed, mdir, nb, velo, iBeta, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, &
                 fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2, &
                 c0, cx1, cy1, cz1, cx2, cy2, cz2, tbeta, omeg
USE InpOutp, ONLY: XS_updt, ounit
USE nodal, ONLY: nodal_coup4, outer4, outertf, outer4ad, PowTot, Fsrc

IMPLICIT NONE

REAL, DIMENSION(nnod, ng) :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Flux at previous time step
REAL, DIMENSION(nnod, ng) :: af                                      ! adjoint flux
REAL, DIMENSION(nnod, ng) :: sigrp                                   ! Temporary sigr

REAL :: rho
REAL :: t1, t2
REAL :: tpow1, tpow2
INTEGER :: n, i, j, g, imax, step

! Allocate precusor density
ALLOCATE (c0(nf,nnod),cx1(nf,nnod),cy1(nf,nnod),cz1(nf,nnod))
ALLOCATE (cx2(nf,nnod),cy2(nf,nnod),cz2(nf,nnod))

! Allocate Frequency transformation constant
ALLOCATE (omeg(nnod,ng))

! Update xsec
CALL XS_updt(bcon, ftem, mtem, cden, bpos)

! Guess fission source
fs0 = 0.; fsx1 = 0.; fsy1 = 0.; fsz1 = 0.; fsx2 = 0.; fsy2 = 0.; fsz2 = 0.
DO g = 1, ng
   CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
END DO

! Calculate forward flux at t=0 and check if keff=1
CALL nodal_coup4()
CALL outer4(0)

! If K-EFF NOT EQUAL TO 1.0
IF (ABS(Ke - 1.) > 1.e-5) CALL KNE1

! Calculate Adjoint flux
CALL outer4ad(0)
af = f0   ! Save adjoint flux to af

! ReCalculate forward flux
CALL outer4(0)

! Calculate Initial precursor density
CALL iPden()

! Calculate initial power
CALL PowTot(f0,tpow1)

! Total beta
tbeta = 0.
DO j = 1, nf
  tbeta = tbeta + iBeta(j)
END DO

! Calculate reactivity
CALL react(af, sigr, rho)

WRITE(ounit, *)
WRITE(ounit, *) " TRANSIENT RESULTS :"
WRITE(ounit, *)
WRITE(ounit, *) " Step  Time(s)  React.($)   Rel. Power   CR Bank Pos. (1-end)"
WRITE(ounit, *) "--------------------------------------------------------------"
WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') 0, 0., 0., &
1.0, (bpos(n), n = 1, nb)

! Start transient calculation
step = 0
t2 = 0.
imax = INT(tdiv/tstep1)

! First Time Step
DO i = 1, imax

    step = step + 1
    t1 = t2
    t2 = REAL(i)*tstep1

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
    sigrp = sigr    ! Save sigr to sigrp
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

    ! Transient calculation
    CALL nodal_coup4()
    CALL outertf(nout, tstep1, ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2)

    ! Update fission source
    fs0 = 0.; fsx1 = 0.; fsy1 = 0.; fsz1 = 0.; fsx2 = 0.; fsy2 = 0.; fsz2 = 0.
    DO g = 1, ng
       CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
    END DO

    ! Update precursor density
    CALL uPden(tstep1)

    ! Calculate power
    CALL PowTot(f0, tpow2)

    ! Calculate reactivity
    CALL react(af, sigrp, rho)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho/tbeta, &
    tpow2/tpow1, (bpos(n), n = 1, nb)


END DO

! Second Time Step
imax = INT((ttot-tdiv)/tstep2)

DO i = 1, imax

    step = step + 1
    t1 = t2
    t2 = tdiv + REAL(i)*tstep2

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

    ! Modify removal xsec
    sigrp = sigr
    DO g = 1, ng
       DO n = 1, nnod
          omeg(n,g) = 0.
          sigr(n,g) = sigr(n,g) + 1. / (velo(g) * tstep2)
       END DO
    END DO

    ! Save the previous flux and DN precusor density
    ft = f0
    ftx1 = fx1; fty1 = fy1; ftz1 = fz1
    ftx2 = fx2; fty2 = fy2; ftz2 = fz2

    ! Transient calculation
    CALL nodal_coup4()
    CALL outertf(nout, tstep2, ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2)

    ! Update fission source
    fs0 = 0.; fsx1 = 0.; fsy1 = 0.; fsz1 = 0.; fsx2 = 0.; fsy2 = 0.; fsz2 = 0.
    DO g = 1, ng
       CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
    END DO

    ! Update precursor density
    CALL uPden(tstep2)

    ! Calculate power
    CALL PowTot(f0, tpow2)

    ! Calculate reactivity
    CALL react(af, sigrp, rho)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho/tbeta, &
    tpow2/tpow1, (bpos(n), n = 1, nb)

END DO

END SUBROUTINE rod_eject


SUBROUTINE trod_eject()

!
! Purpose:
!    To perform rod ejection simulation with TH feedbacks
!

USE sdata, ONLY: ng, nnod, sigr, nf, nout, &
                 ttot, tdiv, tstep1, tstep2, Ke, &
                 bcon, ftem, mtem, cden, tfm, &
                 fbpos, bpos, tmove, bspeed, mdir, nb, velo, iBeta, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, &
                 fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2, &
                 c0, cx1, cy1, cz1, cx2, cy2, cz2, tbeta, omeg, &
                 npow, pow, ppow, node_nf, ix, iy, iz, zdel
USE InpOutp, ONLY: XS_updt, ounit
USE nodal, ONLY: nodal_coup4, outer4, outertf, outer4ad, PowTot, powdis, Fsrc
USE th, ONLY: th_iter, th_trans3, par_ave, par_max, par_ave_f

IMPLICIT NONE

REAL, DIMENSION(nnod, ng) :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Flux at previous time step
REAL, DIMENSION(nnod, ng) :: af                                      ! adjoint flux
REAL, DIMENSION(nnod, ng) :: sigrp                                   ! Temporary sigr

REAL :: rho
REAL :: t1, t2
REAL :: tpow1, tpow2
INTEGER :: n, i, j, g, imax, step

REAL, DIMENSION(nnod) :: pline       ! Linear power density
REAL :: xppow
REAL :: tf, tm, mtf, mtm

! Allocate node power distribution npow
ALLOCATE(npow(nnod))

! Allocate precusor density
ALLOCATE (c0(nf,nnod),cx1(nf,nnod),cy1(nf,nnod),cz1(nf,nnod))
ALLOCATE (cx2(nf,nnod),cy2(nf,nnod),cz2(nf,nnod))

! Allocate Frequency transformation constant
ALLOCATE (omeg(nnod,ng))

! Determine th paramters distribution
CALL th_iter(0)

! Calculate forward flux
CALL outer4(0)

! If K-EFF NOT EQUAL TO 1.0
IF (ABS(Ke - 1.) > 1.e-5) CALL KNE1

! Calculate Adjoint flux
CALL outer4ad(0)
af = f0   ! Save adjoint flux to af

! ReCalculate forward flux
CALL outer4(0)

! Calculate power
CALL PowTot(f0,tpow1)

! Calculate Initial precursor density
CALL iPden()

! Total beta
tbeta = 0.
DO j = 1, nf
  tbeta = tbeta + iBeta(j)
END DO

! ReCalculate forward flux
! CALL outer4(0)

! Calculate reactivity
CALL react(af, sigr, rho)

CALL par_ave_f(ftem, tf)
CALL par_max(tfm(:,1), mtf)
CALL par_ave(mtem, tm)
CALL par_max(mtem, mtm)

WRITE(ounit, *)
WRITE(ounit, *) " TRANSIENT RESULTS :"
WRITE(ounit, *)
WRITE(ounit, *) " Step  Time(s)  React.($)   Rel. Power   Avg. Tm   Max. Tm   Avg. Tf   Max. Tf"
WRITE(ounit, *) "--------------------------------------------------------------------------------"
WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 4F10.2)') 0, 0., 0., &
ppow*0.01, tm-273.15, mtm-273.15, tf-273.15, mtf-273.15

! Start transient calculation
step = 0
t2 = 0.
imax = INT(tdiv/tstep1)

! First Time Step
DO i = 1, imax

    step = step + 1
    t1 = t2
    t2 = REAL(i)*tstep1

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
    sigrp = sigr
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

    ! Transient calculation
    CALL nodal_coup4()
    CALL outertf(nout, tstep1, ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2)

    ! Update fission source
    fs0 = 0.; fsx1 = 0.; fsy1 = 0.; fsz1 = 0.; fsx2 = 0.; fsy2 = 0.; fsz2 = 0.
    DO g = 1, ng
       CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
    END DO

    ! Update precursor density
    CALL uPden(tstep1)

    ! Calculate power
    CALL PowTot(f0, tpow2)

    ! Calculate node power distribution
    CALL powdis(npow)

    ! Power change
    xppow = ppow * tpow2/tpow1 * 0.01

    ! Calculate linear power density for each nodes (W/cm)
    DO n = 1, nnod
       pline(n) = npow(n) * pow * xppow &
       / (node_nf(ix(n),iy(n)) * zdel(iz(n)))     ! Linear power density (W/cm)
    END DO

    ! TH transient
    CALL th_trans3(pline,tstep1)

    ! Calculate reactivity
    CALL react(af, sigrp, rho)

    CALL par_ave_f(ftem, tf)
    CALL par_max(tfm(:,1), mtf)
    CALL par_ave(mtem, tm)
    CALL par_max(mtem, mtm)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 4F10.2)') step, t2, rho/tbeta, &
    xppow, tm-273.15, mtm-273.15, tf-273.15, mtf-273.15

    IF (step>1000) THEN
        WRITE(ounit,*) 'TOO SMALL TIME STEPS. STOPPING'
        STOP
    END IF

END DO

! Second Time Step
imax = INT((ttot-tdiv)/tstep2)

DO i = 1, imax

    step = step + 1
    t1 = t2
    t2 = tdiv + REAL(i)*tstep2

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

    ! Modify removal xsec
    sigrp = sigr
    DO g = 1, ng
       DO n = 1, nnod
          omeg(n,g) = 0.
          sigr(n,g) = sigr(n,g) + 1. / (velo(g) * tstep2)
       END DO
    END DO

    ! Save the previous flux and DN precusor density
    ft = f0
    ftx1 = fx1; fty1 = fy1; ftz1 = fz1
    ftx2 = fx2; fty2 = fy2; ftz2 = fz2

    ! Transient calculation
    CALL nodal_coup4()
    CALL outertf(nout, tstep2, ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2)

    ! Update fission source
    fs0 = 0.; fsx1 = 0.; fsy1 = 0.; fsz1 = 0.; fsx2 = 0.; fsy2 = 0.; fsz2 = 0.
    DO g = 1, ng
       CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
    END DO

    ! Update precursor density
    CALL uPden(tstep2)

    ! Calculate power
    CALL PowTot(f0, tpow2)

    ! Calculate node power distribution
    CALL powdis(npow)

    ! Power change
    xppow = ppow * tpow2/tpow1 * 0.01

    ! Calculate linear power density for each nodes (W/cm)
    DO n = 1, nnod
       pline(n) = npow(n) * pow * xppow &
       / (node_nf(ix(n),iy(n)) * zdel(iz(n)))     ! Linear power density (W/cm)
    END DO

    ! TH transient
    CALL th_trans3(pline,tstep2)

    ! Calculate reactivity
    CALL react(af, sigrp, rho)

    CALL par_ave_f(ftem, tf)
    CALL par_max(tfm(:,1), mtf)
    CALL par_ave(mtem, tm)
    CALL par_max(mtem, mtm)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 4F10.2)') step, t2, rho/tbeta, &
    xppow, tm-273.15, mtm-273.15, tf-273.15, mtf-273.15


    IF (step>1000) THEN
        WRITE(ounit,*) 'TOO SMALL TIME STEPS. STOPPING'
        STOP
    END IF

END DO


END SUBROUTINE trod_eject


SUBROUTINE KNE1()

!
! Purpose:
!    To adjuts the Keff to 1.0 if it is not equal to 1.0
!

USE sdata, ONLY: Ke, xnuf, dnuf, bcon, ftem, mtem, cden, bpos
USE InpOutp, ONLY: XS_updt, ounit
USE nodal, ONLY: outer4

IMPLICIT NONE

INTEGER :: i

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


END SUBROUTINE KNE1


SUBROUTINE iPden()

!
! Purpose:
!    Calculate Initial precursor density
!

USE sdata, ONLY: nnod, nf, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2, &
c0, cx1, cy1, cz1, cx2, cy2, cz2, iBeta, lamb

IMPLICIT NONE

INTEGER :: n, j

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


END SUBROUTINE iPden


SUBROUTINE uPden(h)

!
! Purpose:
!    To update precursor density
!

USE sdata, ONLY: nnod, nf, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2, &
                 c0, cx1, cy1, cz1, cx2, cy2, cz2, iBeta, lamb

IMPLICIT NONE

REAL, INTENT(IN) :: h

REAL :: lat
INTEGER :: n, j

DO n = 1, nnod
   DO j = 1, nf
      lat = (1. + lamb(j) * h)
      c0(j,n) = (c0(j,n) + iBeta(j) * h * fs0(n)) / lat
      cx1(j,n) = (cx1(j,n) + iBeta(j) * h * fsx1(n)) / lat
      cy1(j,n) = (cy1(j,n) + iBeta(j) * h * fsy1(n)) / lat
      cz1(j,n) = (cz1(j,n) + iBeta(j) * h * fsz1(n)) / lat
      cx2(j,n) = (cx2(j,n) + iBeta(j) * h * fsx2(n)) / lat
      cy2(j,n) = (cy2(j,n) + iBeta(j) * h * fsy2(n)) / lat
      cz2(j,n) = (cz2(j,n) + iBeta(j) * h * fsz2(n)) / lat
   END DO
END DO


END SUBROUTINE uPden


SUBROUTINE react(af,sigrp,rho)

!
! Purpose:
!    To calculate dynamic reactivity
!

USE sdata, ONLY: nnod, ng, f0, sigs, nod, chi, fs0, &
                 vdel, xdel, ydel, zdel, ix, iy, iz

IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(IN) :: af
REAL, DIMENSION(:,:), INTENT(IN) :: sigrp
REAL, INTENT(OUT) :: rho

INTEGER :: n, g, h
REAL, DIMENSION(nnod) :: scg
REAL :: rem, lea, src, fde

src = 0.; rem = 0.; lea = 0.; fde = 0.
DO g = 1, ng
   scg = 0.
   DO h = 1, ng
      DO n = 1, nnod
         IF (g /= h) scg(n) = scg(n) + sigs(n,h,g) * f0(n,h)
      END DO
   END DO
   DO n = 1, nnod
      src = src + af(n,g) * (scg(n) + chi(n,g) * fs0(n)) * vdel(n)
      rem = rem + af(n,g) * sigrp(n,g) * f0(n,g) * vdel(n)
      lea = lea + af(n,g) * &
                      (nod(n,g)%L(1) * ydel(iy(n)) * zdel(iz(n)) + &
                       nod(n,g)%L(2) * xdel(ix(n)) * zdel(iz(n)) + &
                       nod(n,g)%L(3) * xdel(ix(n)) * ydel(iy(n)))
      fde = fde + af(n,g) * chi(n,g) * fs0(n) * vdel(n)
    END DO
END DO

rho = (src - lea - rem) / fde

END SUBROUTINE react


SUBROUTINE kreact()

!
! Purpose:
!    To calculate reactivity from Keff
!    (used only for to verify react SR)
!    shall be used only to calculate control rod worth, not to be used with transient calc.
!

USE sdata, ONLY: ng, nf, ttot, tstep1, Ke, &
                 bcon, ftem, mtem, cden, &
                 fbpos, bpos, tmove, bspeed, mdir, nb, iBeta, &
                 fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2, tbeta
USE InpOutp, ONLY: XS_updt, ounit
USE nodal, ONLY: nodal_coup4, outer4, Fsrc

IMPLICIT NONE

REAL :: t1, t2
INTEGER :: n, i, j, g, imax, step

! Update xsec
CALL XS_updt(bcon, ftem, mtem, cden, bpos)

! Guess fission source
fs0 = 0.; fsx1 = 0.; fsy1 = 0.; fsz1 = 0.; fsx2 = 0.; fsy2 = 0.; fsz2 = 0.
DO g = 1, ng
   CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
END DO

! Calculate forward flux at t=0 and check if keff=1
CALL nodal_coup4()
CALL outer4(0)

! If K-EFF NOT EQUAL TO 1.0
IF (ABS(Ke - 1.) > 1.e-5) CALL KNE1

! ReCalculate forward flux
CALL outer4(0)

! Total beta
tbeta = 0.
DO j = 1, nf
  tbeta = tbeta + iBeta(j)
END DO

WRITE(ounit,'(I4, F12.5, F10.4)') 0, Ke, ((Ke-1.)/Ke)/tbeta

! Start transient calculation
step = 0
t2 = 0.
imax = INT(ttot/tstep1)

! First Time Step
DO i = 1, imax

    step = step + 1
    t1 = t2
    t2 = REAL(i)*tstep1

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


    ! Transient calculation
    CALL nodal_coup4()
    CALL outer4(0)

    WRITE(ounit,'(I4, F12.5, F10.4)') step, Ke, ((Ke-1.)/Ke)/tbeta


END DO

END SUBROUTINE kreact


END MODULE trans
