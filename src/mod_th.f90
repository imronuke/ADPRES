MODULE th

IMPLICIT NONE

SAVE

CONTAINS

SUBROUTINE th_iter(ind)

  !
  ! Purpose:
  !    To do thermal-hydrailics iteration
  !

  USE sdata, ONLY: nnod, ftem, mtem, cden, bcon, bpos, npow, pow, ppow,  &
                   zdel, node_nf, ix, iy, iz, th_err, node_nf, ix, iy, iz, th_niter
  USE nodal, ONLY: nodal_coup4, outer4th, PowDis
  USE InpOutp, ONLY: XS_updt, ounit

  IMPLICIT NONE

  INTEGER, INTENT(IN), OPTIONAL :: ind    ! if iteration reaching th_iter and ind = 0 then STOP
  DOUBLE PRECISION, DIMENSION(nnod) :: pline
  DOUBLE PRECISION, DIMENSION(nnod) :: otem
  INTEGER :: n, l

  th_err = 1.
  DO l = 1, th_niter
      ! Save old fuel temp
      otem = ftem

      ! Update XS
      CALL XS_updt(bcon, ftem, mtem, cden, bpos)

      ! Update nodal couplings
      CALL nodal_coup4()

      ! Perform outer inner iteration
      CALL outer4th(20)

      ! Calculate power density
      CALL PowDis(npow)

      ! Calculate linear power density for each nodes (W/cm)
      DO n = 1, nnod
          pline(n) = npow(n) * pow * ppow * 0.01 &
                   / (node_nf(ix(n),iy(n)) * zdel(iz(n)))     ! Linear power density (W/cm)
      END DO

      ! Update fuel, moderator temp. and coolant density
      CALL th_upd(pline)

      CALL AbsE(ftem, otem, th_err)

      IF (th_err < 0.01) EXIT

  END DO
  IF (PRESENT(ind)) THEN
     IF ((ind == 0) .AND. (l >= 20)) THEN
        WRITE(ounit,*) '  MAXIMUM TH ITERATION REACHED.'
        WRITE(ounit,*) '  CALCULATION MIGHT BE NOT CONVERGED OR CHANGE ITERATION CONTROL'
        STOP
     END IF
  END IF



END SUBROUTINE th_iter


SUBROUTINE AbsE(newF, oldF, rel)

  !
  ! Purpose:
  !    To calculate Max Relative error

USE sdata, ONLY: nnod

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: newF, oldF
DOUBLE PRECISION, INTENT(OUT) :: rel

DOUBLE PRECISION :: error
INTEGER :: n

rel = 0.

DO n= 1, nnod
    IF (ABS(newF(n)) > 1.d-10) THEN
        error = ABS(newF(n) - oldF(n))
        IF (error > rel) rel = error
    END IF
END DO

END SUBROUTINE AbsE


SUBROUTINE par_ave_f(par, ave)
!
! Purpose:
!    To calculate average fuel temp (only for active core)
!

USE sdata, ONLY: vdel, nnod, ng, nuf

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: par
DOUBLE PRECISION, INTENT(OUT) :: ave
DOUBLE PRECISION :: dum, dum2
INTEGER :: n

dum = 0.; dum2 = 0.
DO n = 1, nnod
   IF (nuf(n,ng) > 0.) THEN
      dum = dum + par(n) * vdel(n)
      dum2 = dum2 + vdel(n)
   END IF
END DO

ave = dum / dum2

END SUBROUTINE par_ave_f


SUBROUTINE par_ave(par, ave)
!
! Purpose:
!    To calculate average moderator temp (only for radially active core)
!

USE sdata, ONLY: vdel, nnod

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: par
DOUBLE PRECISION, INTENT(OUT) :: ave
DOUBLE PRECISION :: dum, dum2
INTEGER :: n

dum = 0.; dum2 = 0.
DO n = 1, nnod
   dum = dum + par(n) * vdel(n)
   dum2 = dum2 + vdel(n)
END DO

ave = dum / dum2

END SUBROUTINE par_ave


SUBROUTINE par_max(par, pmax)
!
! Purpose:
!    To calculate maximum fuel tem, coolant tem, and density
!

USE sdata, ONLY: nnod

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: par
DOUBLE PRECISION, INTENT(OUT) :: pmax
INTEGER :: n

pmax = 0.
DO n = 1, nnod
   IF (par(n) > pmax) pmax = par(n)
END DO

END SUBROUTINE par_max


SUBROUTINE getent(t,ent)
!
! Purpose:
!    To get enthalpy for given coolant temp. from steam table
!

USE sdata, ONLY: stab, ntem
USE InpOutp, ONLY : ounit

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: t
DOUBLE PRECISION, INTENT(OUT) :: ent
DOUBLE PRECISION :: t1, ent1
DOUBLE PRECISION :: t2, ent2
INTEGER :: i

IF ((t < 473.15) .OR. (t > 617.91)) THEN
    WRITE(ounit,*) '  Coolant temp. : ', t
    WRITE(ounit,*) '  ERROR : MODERATOR TEMP. IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
    WRITE(ounit,*) '  CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
    WRITE(*,*) '  Coolant temp. : ', t
    WRITE(*,*) '  ERROR : MODERATOR TEMP. IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
    WRITE(*,*) '  CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
    STOP
END IF

t2 = stab(1,1); ent2 = stab(1,3)
DO i = 2, ntem
    t1 = t2
    ent1 = ent2
    t2 = stab(i,1); ent2 = stab(i,3)
    IF ((t >= t1) .AND. (t <= t2)) THEN
        ent = ent1 + (t - t1) / (t2 - t1) * (ent2 - ent1)
        EXIT
    END IF
END DO


END SUBROUTINE getent


SUBROUTINE gettd(ent,t,rho,prx,kvx,tcx)
!
! Purpose:
!    To get enthalpy for given coolant temp. from steam table
!

USE sdata, ONLY: stab, ntem
USE InpOutp, ONLY : ounit

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: ent
DOUBLE PRECISION, INTENT(OUT) :: t, rho, prx, kvx, tcx
DOUBLE PRECISION :: t1, rho1, ent1, kv1, pr1, tc1
DOUBLE PRECISION :: t2, rho2, ent2, kv2, pr2, tc2
DOUBLE PRECISION :: ratx

INTEGER :: i

IF ((ent < 858341.5) .OR. (ent > 1624307.1)) THEN
    WRITE(ounit,*) '  Enthalpy. : ', ent
    WRITE(ounit,*) '  ERROR : ENTHALPY IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
    WRITE(ounit,*) '  CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
    WRITE(*,*) '  Enthalpy. : ', ent
    WRITE(*,*) '  ERROR : ENTHALPY IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
    WRITE(*,*) '  CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
    STOP
END IF




t2 = stab(1,1); rho2 = stab(1,2); ent2 = stab(1,3)
pr2 = stab(1,4); kv2 = stab(1,5); tc2 = stab(1,6)
DO i = 2, ntem
    t1 = t2
    ent1 = ent2
    rho1 = rho2
    pr1 = pr2
    kv1 = kv2
    tc1 = tc2
    t2 = stab(i,1); rho2 = stab(i,2); ent2 = stab(i,3)
    pr2 = stab(i,4); kv2 = stab(i,5); tc2 = stab(i,6)
    IF ((ent >= ent1) .AND. (ent <= ent2)) THEN
        ratx = (ent - ent1) / (ent2 - ent1)
        t   = t1   + ratx * (t2 - t1)
        rho = rho1 + ratx * (rho2 - rho1)
        prx = pr1  + ratx * (pr2 - pr1)
        kvx = kv1  + ratx * (kv2 - kv1)
        tcx = tc1 + ratx * (tc2 - tc1)
        EXIT
    END IF
END DO


END SUBROUTINE gettd


DOUBLE PRECISION FUNCTION getkc(t)
!
! Purpose:
!    To calculate thermal conductivity of cladding
!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: t

getkc = 7.51d0 + 2.09d-2*t - 1.45d-5*t**2 + 7.67d-9*t**3

END FUNCTION getkc


DOUBLE PRECISION FUNCTION getkf(t)
!
! Purpose:
!    To calculate thermal conductivity of fuel
!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: t

getkf = 1.05d0 + 2150.0d0 / (t - 73.15d0)

END FUNCTION getkf


DOUBLE PRECISION FUNCTION getcpc(t)
!
! Purpose:
!    To calculate specific heat capacity of cladding
!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: t

getcpc = 252.54d0 + 0.11474d0*t

END FUNCTION getcpc


DOUBLE PRECISION FUNCTION getcpf(t)
!
! Purpose:
!    To calculate specific heat capacity of fuel
!

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: t

getcpf = 162.3d0 + 0.3038d0*t - 2.391d-4*t**2 + 6.404d-8*t**3

END FUNCTION getcpf


SUBROUTINE TridiaSolve(a,b,c,d,x)
!
! Purpose:
!    To solve tridiagonal matrix
!

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: a, b, c, d
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: x

INTEGER :: i, n

n = SIZE(d)

! Gauss Elimination
c(1) = c(1)/b(1)
d(1) = d(1)/b(1)
DO i = 2, n
    c(i) = c(i) / (b(i) - a(i) * c(i-1))
    d(i) = (d(i) - a(i) * d(i-1)) / (b(i) - a(i) * c(i-1))
END DO

! Back Substitution
x(n) = d(n)
DO i = n-1, 1, -1
    x(i) = d(i) - c(i) * x(i+1)
END DO

END SUBROUTINE TridiaSolve



DOUBLE PRECISION FUNCTION geths(xden, tc, kv, Pr)
!
! Purpose:
!    To calculate heat transfer coef.
!

USE sdata, ONLY: dh, farea, cflow

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: xden  ! coolant densisty
DOUBLE PRECISION, INTENT(IN) :: tc  ! coolant thermal conductivity
DOUBLE PRECISION, INTENT(IN) :: kv  ! kinematic viscosity
DOUBLE PRECISION, INTENT(IN) :: Pr  ! Prandtl Number

DOUBLE PRECISION :: cvelo, Nu, Re

cvelo = cflow / (farea * xden * 1000.d0)        ! Calculate flow velocity (m/s)
Re = cvelo * dh / (kv * 1.d-6)                 ! Calculate Reynolds Number
Nu = 0.023d0*(Pr**0.4d0)*(Re**0.8d0)                ! Calculate Nusselt Number
geths = (tc / dh) * Nu                        ! Calculate heat transfer coefficient


END FUNCTION geths



SUBROUTINE th_trans(xpline, h)

!
! Purpose:
!    To perform fuel pin thermal transient
!

USE sdata, ONLY: mtem, cden, ftem, tin, xyz, cflow, nyy, nzz, nxx, cf, ent, heatf, nnod, &
                 ystag, tfm, nt, rpos, rdel, rf, rg, rc, farea, dia, pi, zdel, ystag

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: xpline    ! Linear Power Density (W/cm)
DOUBLE PRECISION, INTENT(IN) :: h                       ! Time step

INTEGER :: i, j, k, n
DOUBLE PRECISION, DIMENSION(nt+1) :: a, b, c, d
DOUBLE PRECISION :: hs, hg = 1.d4, kt , kt1, kt2          ! coolant heat transfer coef., gap heat transfer coef, and thermal conductivity
DOUBLE PRECISION :: alpha = 0.7d0
DOUBLE PRECISION :: xa, xc
DOUBLE PRECISION :: fdens = 10.412e3            ! UO2 density (kg/m3)
DOUBLE PRECISION :: cdens = 6.6e3               ! Cladding density (kg/m3)
DOUBLE PRECISION :: cp                          ! Specific heat capacity
DOUBLE PRECISION :: eps, eta
DOUBLE PRECISION :: mdens, vol                  ! Coolant density and channel volume
DOUBLE PRECISION, DIMENSION(nnod) :: entp        ! previous enthalpy

DOUBLE PRECISION :: pdens      ! power densisty  (W/m3)
DOUBLE PRECISION :: enti       ! Coolant inlet enthalpy
DOUBLE PRECISION, DIMENSION(nxx, nyy) :: entm
DOUBLE PRECISION :: cpline     ! Coolant Linear power densisty (W/m)
DOUBLE PRECISION :: Pr, kv, tcon ! Coolant Prandtl Number, Kinematic viscosity, and thermal conductivity

CALL getent(tin, enti)
entp = ent

DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax

            mdens = cden(xyz(i,j,k)) * 1000.d0                                    ! Coolant density (kg/m3)
            cpline = heatf(xyz(i,j,k)) * pi * dia  &
                   + cf * xpline(xyz(i,j,k)) * 100.d0       ! Coolant Linear power densisty (W/m)
            vol   = farea * zdel(k) * 0.01d0

            IF (k == 1) THEN                                                    ! Calculate coolant enthalpy
                eps = mdens * vol / h
                ent(xyz(i,j,k)) = (cpline * zdel(k) * 0.01d0 &
                                + 2.d0 * cflow * enti &
                                + eps * entp(xyz(i,j,k))) &
                                / (eps + 2.d0 * cflow)
                CALL gettd(ent(xyz(i,j,k)), mtem(xyz(i,j,k)), &
                          cden(xyz(i,j,k)), Pr, kv, tcon)                       ! Get corresponding temp and density
                entm(i,j) = 2.d0 * ent(xyz(i,j,k)) - enti
            ELSE
                eps = mdens * vol / h
                ent(xyz(i,j,k)) = (cpline * zdel(k) * 0.01d0 &
                                + 2.d0 * cflow * entm(i,j) &
                                + eps * entp(xyz(i,j,k))) &
                                / (eps + 2.d0 * cflow)
                CALL gettd(ent(xyz(i,j,k)), mtem(xyz(i,j,k)), &
                          cden(xyz(i,j,k)), Pr, kv, tcon)                       ! Get corresponding temp and density
                entm(i,j) = 2.d0 * ent(xyz(i,j,k)) - entm(i,j)
            END IF


            hs = geths(cden(xyz(i,j,k)), Pr, kv, tcon)                                               ! Calculate heat transfer coef
            pdens = (1. - cf) * 100.d0 * xpline(xyz(i,j,k)) / (pi * rf**2)                ! Fuel pin Power Density (W/m3)

            ! Calculate tridiagonal matrix: a, b, c and source: d
            ! For nt=1 [FUEL CENTERLINE]
            kt1 = getkf(tfm(xyz(i,j,k),1))                                                     ! Get thermal conductivity
            kt2 = getkf(tfm(xyz(i,j,k),2))
            kt  = 2.d0 * kt1 * kt2 / (kt1 + kt2)
            cp = getcpf(tfm(xyz(i,j,k),1))                                                           ! Get specific heat capacity
            eta = fdens * cp * rpos(1)**2 / (2. * h)
            xc  = kt * rpos(1) / rdel(1)
            b(1) =  xc + eta
            c(1) = -xc
            d(1) = pdens * 0.5d0 * rpos(1)**2 + eta * tfm(xyz(i,j,k),1)

            DO n = 2, nt-2
                kt1 = kt2
                kt2 = getkf(tfm(xyz(i,j,k),n+1))
                kt  = 2.d0 * kt1 * kt2 / (kt1 + kt2)
                cp = getcpf(tfm(xyz(i,j,k),n))
                eta = fdens * cp * (rpos(n)**2 - rpos(n-1)**2) / (2. * h)
                xa = xc
                xc = kt * rpos(n) / rdel(n)
                a(n) = -xa
                b(n) =  xa + xc + eta
                c(n) = -xc
                d(n) = pdens * 0.5d0 * (rpos(n)**2 - rpos(n-1)**2) &
                     + eta * tfm(xyz(i,j,k),n)
            END DO

            ! For nt-1 [FUEL-GAP INTERFACE]
            cp = getcpf(tfm(xyz(i,j,k),nt-1))
            eta = fdens * cp * (rf**2 - rpos(nt-2)**2) / (2. * h)
            xa = xc
            xc = rg * hg
            a(nt-1) = -xa
            b(nt-1) =  xa + xc + eta
            c(nt-1) = -xc
            d(nt-1) = pdens * 0.5d0 * (rf**2 - rpos(nt-2)**2) &
                    + eta * tfm(xyz(i,j,k),nt-1)

            ! For nt [GAP-CLADDING INTERFACE]
            kt1 = getkc(tfm(xyz(i,j,k),nt))
            kt2 = getkc(tfm(xyz(i,j,k),nt+1))
            kt  = 2.d0 * kt1 * kt2 / (kt1 + kt2)     ! For cladding
            cp = getcpc(tfm(xyz(i,j,k),nt))
            eta = cdens * cp * (rpos(nt)**2 - rg**2) / (2. * h)
            xa = xc
            xc = kt * rpos(nt) / rdel(nt)
            a(nt) = -xa
            b(nt) =  xa + xc + eta
            c(nt) = -xc
            d(nt) = eta * tfm(xyz(i,j,k),nt)

            ! For nt+1  [CLADDING-COOLANT INTERFACE]
            cp = getcpc(tfm(xyz(i,j,k),nt+1))
            eta = cdens * cp * (rc**2 - rpos(nt)**2) / (2. * h)
            xa = xc
            xc = rc * hs
            a(nt+1) = -xa
            b(nt+1) =  xa + xc + eta
            d(nt+1) = rc * hs * mtem(xyz(i,j,k)) &
                    + eta * tfm(xyz(i,j,k),nt+1)

            ! Solve tridiagonal matrix
            CALL TridiaSolve(a, b, c, d, tfm(xyz(i,j,k), :))

            ! Get lumped fuel temp
            ftem(xyz(i,j,k)) = (1.-alpha) * tfm(xyz(i,j,k), 1) &
                             + alpha * tfm(xyz(i,j,k), nt-1)

            ! Calculate heat flux
            heatf(xyz(i,j,k)) = hs * (tfm(xyz(i,j,k), nt+1) - mtem(xyz(i,j,k)))

        END DO
    END DO
END DO

END SUBROUTINE th_trans


SUBROUTINE th_upd(xpline)

!
! Purpose:
!    To update thermal parameters
!

USE sdata, ONLY: mtem, cden, ftem, tin, xyz, cflow, nyy, nxx, nzz, cf, ent, heatf, &
                 ystag, tfm, nt, rpos, rdel, rf, rg, rc, pi, zdel, dia, ystag, &
                 farea

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: xpline    ! Linear Power Density (W/cm)

INTEGER :: i, j, k, n
DOUBLE PRECISION, DIMENSION(nt+1) :: a, b, c, d
DOUBLE PRECISION :: hs, Hg = 1.d4, kt, kt1, kt2
DOUBLE PRECISION :: alp = 0.7d0
DOUBLE PRECISION :: xa, xc
DOUBLE PRECISION :: pdens      ! power densisty  (W/m3)
DOUBLE PRECISION :: enti       ! Coolant inlet enthalpy
DOUBLE PRECISION, DIMENSION(nxx, nyy) :: entm
DOUBLE PRECISION :: cpline     ! Coolant Linear power densisty (W/m)
DOUBLE PRECISION :: Pr, kv, tcon ! Coolant Prandtl Number, Kinematic viscosity, and thermal conductivity

CALL getent(tin, enti)

DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax

            cpline = heatf(xyz(i,j,k)) * pi * dia  &
                   + cf * xpline(xyz(i,j,k)) * 100.d0       ! Coolant Linear power densisty (W/m)

            IF (k == 1) THEN                                                    ! Calculate coolant enthalpy and
                ent(xyz(i,j,k)) = enti + 0.5d0 * cpline * zdel(k) * 0.01d0 / cflow
                CALL gettd(ent(xyz(i,j,k)), mtem(xyz(i,j,k)), cden(xyz(i,j,k)), &
                          Pr, kv, tcon)                                         ! Get corresponding temp and density
                entm(i,j) = 2.d0 * ent(xyz(i,j,k)) - enti
            ELSE
                ent(xyz(i,j,k)) = entm(i,j) + 0.5d0 * cpline * zdel(k) * 0.01d0 / cflow
                CALL gettd(ent(xyz(i,j,k)), mtem(xyz(i,j,k)), cden(xyz(i,j,k)), &
                          Pr, kv, tcon)
                entm(i,j) = 2.d0 * ent(xyz(i,j,k)) - entm(i,j)
            END IF

            hs = geths(cden(xyz(i,j,k)), Pr, kv, tcon)
            pdens = (1. - cf) * 100.d0 * xpline(xyz(i,j,k)) / (pi * rf**2)        ! Fuel pin Power Density (W/m3)

            ! Calculate tridiagonal matrix: a, b, c and source: d
            kt1 = getkf(tfm(xyz(i,j,k),1))                                                     ! Get thermal conductivity
            kt2 = getkf(tfm(xyz(i,j,k),2))
            kt  = 2.d0 * kt1 * kt2 / (kt1 + kt2)
            xc  = kt * rpos(1) / rdel(1)
            b(1) =  xc
            c(1) = -xc
            d(1) = pdens * 0.5d0 * rpos(1)**2

            DO n = 2, nt-2
                kt1 = kt2
                kt2 = getkf(tfm(xyz(i,j,k),n+1))
                kt  = 2.d0 * kt1 * kt2 / (kt1 + kt2)
                xa = xc
                xc = kt * rpos(n) / rdel(n)
                a(n) = -xa
                b(n) =  xa + xc
                c(n) = -xc
                d(n) = pdens * 0.5d0 * (rpos(n)**2 - rpos(n-1)**2)
            END DO

            ! For nt-1 [FUEL-GAP INTERFACE]
            xa = xc
            xc = rg * Hg
            a(nt-1) = -xa
            b(nt-1) =  xa + xc
            c(nt-1) = -xc
            d(nt-1) = pdens * 0.5d0 * (rf**2 - rpos(nt-2)**2)

            ! For nt [GAP-CLADDING INTERFACE]
            kt1 = getkc(tfm(xyz(i,j,k),nt))
            kt2 = getkc(tfm(xyz(i,j,k),nt+1))
            kt  = 2.d0 * kt1 * kt2 / (kt1 + kt2)     ! For cladding
            xa = xc
            xc = kt * rpos(nt) / rdel(nt)
            a(nt) = -xa
            b(nt) =  xa + xc
            c(nt) = -xc
            d(nt) = 0.

            ! For nt+1  [CLADDING-COOLANT INTERFACE]
            xa = xc
            a(nt+1) = -xa
            b(nt+1) =  xa + hs * rc
            d(nt+1) = rc * hs * mtem(xyz(i,j,k))

            ! Solve tridiagonal matrix
            CALL TridiaSolve(a, b, c, d, tfm(xyz(i,j,k), :))

            ! Get lumped fuel temp
            ftem(xyz(i,j,k)) = (1.-alp) * tfm(xyz(i,j,k), 1) + alp * tfm(xyz(i,j,k), nt-1)

            ! Calculate heat flux
            heatf(xyz(i,j,k)) = hs * (tfm(xyz(i,j,k), nt+1) - mtem(xyz(i,j,k)))


        END DO
    END DO
END DO


END SUBROUTINE th_upd



SUBROUTINE cbsearch()

!
! Purpose:
!    To search critical boron concentration
!

USE sdata, ONLY: Ke, rbcon, ftem, mtem, cden, bpos, nnod, f0, fer, ser, &
                 aprad, apaxi, afrad, npow
USE InpOutp, ONLY: ounit, XS_updt, AsmFlux, AsmPow, AxiPow
USE nodal, ONLY: nodal_coup4, outer4, powdis

IMPLICIT NONE

DOUBLE PRECISION  :: bc1, bc2, bcon     ! Boron Concentration
DOUBLE PRECISION :: ke1, ke2
INTEGER :: n

! File Output
WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) ' ==============================================' &
            // '=================================================='
WRITE(ounit,*) &
               '                       CRITICAL BORON CONCENTRATION SEARCH'
WRITE(ounit,*) ' ==============================================' &
            // '=================================================='
WRITE(ounit,*)
WRITE(ounit,*) '  Itr  Boron Concentration          K-EFF    FLUX REL. ERROR' &
               //'   FISS. SOURCE REL. ERROR    DOPPLER ERROR'
WRITE(ounit,*) ' -----------------------------------------------------------' &
              // '-------------------------------------------'

! Terminal Output
WRITE(*,*)
WRITE(*,*)
WRITE(*,*) ' ==============================================' &
            // '=========='
WRITE(*,*) &
               '           CRITICAL BORON CONCENTRATION SEARCH'
WRITE(*,*) ' ==============================================' &
            // '=========='
WRITE(*,*)
WRITE(*,*) '  Itr  Boron Concentration          K-EFF    '
WRITE(*,*) ' --------------------------------------------------------'


bcon = rbcon
CALL XS_updt(bcon, ftem, mtem, cden, bpos)
CALL nodal_coup4()
CALL outer4(0)
bc1 = bcon
ke1 = Ke

WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') 1, bc1, Ke1, ser, fer
WRITE(*,'(I5, F15.2, F23.5)') 1, bc1, Ke1

bcon = bcon + (Ke - 1.) * bcon   ! Guess next critical boron concentration
CALL XS_updt(bcon, ftem, mtem, cden, bpos)
CALL nodal_coup4()
CALL outer4(0)
bc2 = bcon
ke2 = Ke

WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') 2, bc2, Ke2, ser, fer
WRITE(*,'(I5, F15.2, F23.5)') 2, bc2, Ke2

n = 3
DO
  bcon = bc2 + (1.0 - ke2) / (ke1 - ke2) * (bc1 - bc2)
  CALL XS_updt(bcon, ftem, mtem, cden, bpos)
  CALL nodal_coup4()
  CALL outer4(0)
  bc1 = bc2
  bc2 = bcon
  ke1 = ke2
  ke2 = ke
  WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') n, bcon, Ke, ser, fer
  WRITE(*,'(I5, F15.2, F23.5)') n, bcon, Ke
    IF ((ABS(Ke - 1.0) < 1.d-5) .AND. (ser < 1.d-5) .AND. (fer < 1.d-5)) EXIT
    n = n + 1
    IF (bcon > 3000.) THEN
        WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
        WRITE(ounit,*) '  ADPRES IS STOPPING'
        WRITE(*,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
        STOP
    END IF
    IF (bcon < 0.) THEN
        WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
        WRITE(ounit,*) '  ADPRES IS STOPPING'
        WRITE(*,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
        STOP
    END IF
    IF (n == 20) THEN
        WRITE(ounit,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
        WRITE(ounit,*) '  ADPRES IS STOPPING'
        WRITE(*,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
        STOP
    END IF
END DO

ALLOCATE(npow(nnod))
IF (aprad == 1 .OR. apaxi == 1) THEN
    CALL PowDis(npow)
END IF

IF (aprad == 1) CALL AsmPow(npow)

IF (apaxi == 1) CALL AxiPow(npow)

IF (afrad == 1) CALL AsmFlux(f0, 1.d0)

END SUBROUTINE cbsearch


SUBROUTINE cbsearcht()

!
! Purpose:
!    To search critical boron concentration with thermal feedback
!

USE sdata, ONLY: Ke, ftem, mtem, bcon, rbcon, npow, nnod, &
                 f0, ser, fer, tfm, aprad, apaxi, afrad, npow, th_err, &
                 serc, ferc
USE InpOutp, ONLY: ounit, XS_updt, AsmFlux, AsmPow, AxiPow, getfq
USE nodal, ONLY: powdis, nodal_coup4, outer4

IMPLICIT NONE

DOUBLE PRECISION  :: bc1, bc2    ! Boron Concentration
DOUBLE PRECISION :: ke1, ke2
INTEGER :: n
DOUBLE PRECISION :: tf, tm, mtm, mtf

! File Output
WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) ' ==============================================' &
            // '=================================================='
WRITE(ounit,*) &
               '                       CRITICAL BORON CONCENTRATION SEARCH'
WRITE(ounit,*) ' ==============================================' &
            // '=================================================='
WRITE(ounit,*)
WRITE(ounit,*) '  Itr  Boron Concentration          K-EFF    FLUX REL. ERROR' &
               //'   FISS. SOURCE REL. ERROR    DOPPLER ERROR'
WRITE(ounit,*) ' -----------------------------------------------------------' &
              // '-------------------------------------------'

! Terminal Output
WRITE(*,*)
WRITE(*,*)
WRITE(*,*) ' ==============================================' &
            // '=========='
WRITE(*,*) &
               '           CRITICAL BORON CONCENTRATION SEARCH'
WRITE(*,*) ' ==============================================' &
            // '=========='
WRITE(*,*)
WRITE(*,*) '  Itr  Boron Concentration          K-EFF    '
WRITE(*,*) ' --------------------------------------------------------'


ALLOCATE(npow(nnod))

bcon = rbcon
CALL th_iter()  ! Start thermal hydarulic iteration with current paramters
bc1 = bcon
ke1 = Ke

WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') 1, bc1, Ke1, ser, fer, th_err
WRITE(*,'(I5, F15.2, F23.5)') 1, bc1, Ke1

bcon = bcon + (Ke - 1.) * bcon   ! Guess next critical boron concentration
CALL th_iter()                 ! Perform second thermal hydarulic iteration with updated parameters
bc2 = bcon
ke2 = Ke

WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') 2, bc2, Ke2, ser, fer, th_err
WRITE(*,'(I5, F15.2, F23.5)') 2, bc2, Ke2

n = 3
DO
    bcon = bc2 + (1.0 - ke2) / (ke1 - ke2) * (bc1 - bc2)
    CALL th_iter()
    bc1 = bc2
    bc2 = bcon
    ke1 = ke2
    ke2 = ke
    WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') n, bcon, Ke, ser, fer, th_err
    WRITE(*,'(I5, F15.2, F23.5)') n, bcon, Ke
    IF ((ABS(Ke - 1.0) < 1.d-5) .AND. (ser < serc) .AND. (fer < ferc)) EXIT
    n = n + 1
    IF (bcon > 3000.) THEN
        WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
        WRITE(ounit,*) '  ADPRES IS STOPPING'
        WRITE(*,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
        STOP
    END IF
    IF (bcon < 0.) THEN
        WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
        WRITE(ounit,*) '  ADPRES IS STOPPING'
        WRITE(*,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
        STOP
    END IF
    IF (n == 30) THEN
        WRITE(ounit,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
        WRITE(ounit,*) '  ADPRES IS STOPPING'
        WRITE(*,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
        STOP
    END IF
END DO

IF (aprad == 1 .OR. apaxi == 1) THEN
    CALL PowDis(npow)
END IF

IF (aprad == 1) CALL AsmPow(npow)

IF (apaxi == 1) CALL AxiPow(npow)

IF (afrad == 1) CALL AsmFlux(f0, 1.d0)

CALL par_ave_f(ftem, tf)
CALL par_ave(mtem, tm)

CALL par_max(tfm(:,1), mtf)
CALL par_max(mtem, mtm)
CALL getfq(npow)

! Write Output
WRITE(ounit,*)
WRITE(ounit, 5001) tf, tf-273.15
WRITE(ounit, 5002)  mtf, mtf-273.15
WRITE(ounit, 5003) tm, tm-273.15
WRITE(ounit, 5004) mtm, mtm-273.15

5001 FORMAT(2X, 'AVERAGE DOPPLER TEMPERATURE     : ', F7.1, ' K (', F7.1, ' C)')
5002 FORMAT(2X, 'MAX FUEL CENTERLINE TEMPERATURE : ', F7.1, ' K (', F7.1, ' C)')
5003 FORMAT(2X, 'AVERAGE MODERATOR TEMPERATURE   : ', F7.1, ' K (', F7.1, ' C)')
5004 FORMAT(2X, 'MAXIMUM MODERATOR TEMPERATURE   : ', F7.1, ' K (', F7.1, ' C)')


END SUBROUTINE cbsearcht


END MODULE th
