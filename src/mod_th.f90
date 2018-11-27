MODULE th

IMPLICIT NONE

SAVE

REAL :: kv   ! Water kinematic viscosity  (10e-6 m2/s)
REAL :: Pr  !Prandtl Number
REAL :: tc  ! Thermal conductivity (W/mK)

INTEGER, PARAMETER:: npres = 2   ! Number of pressure data in steam table
INTEGER, PARAMETER:: npres2 = 3   ! Number of pressure data for Prandtl number
INTEGER, PARAMETER:: ntem = 16   ! Number of temperature in steam table
INTEGER, PARAMETER:: ntem2 = 7   ! Number of temperature in steam table for Prandtl number
TYPE :: WATER_DATA
    REAL :: press   ! Pressure
    REAL, DIMENSION(ntem) :: tem  ! TEMPERATURE
    REAL, DIMENSION(ntem) :: dens  ! density
    REAL, DIMENSION(ntem) :: enth  ! Enthalpy
    REAL, DIMENSION(ntem) :: tcon  ! thermal conductivity
    REAL, DIMENSION(ntem) :: pran  ! Prandtl Numbers
    REAL, DIMENSION(ntem) :: kinvis  ! kinematic viscosity
END TYPE
TYPE(WATER_DATA) :: wdata  ! water thermophysical data

CONTAINS

SUBROUTINE th_iter(ind)

  !
  ! Purpose:
  !    To do thermal-hydrailics iteration
  !

  USE sdata, ONLY: nnod, ftem, mtem, cden, bcon, bpos, npow, pow, ppow, nout, nac,  &
                   zdel, node_nf, ix, iy, iz, th_err, node_nf, ix, iy, iz, th_niter
  USE nodal, ONLY: nodal_coup4, outer4th, PowDis
  USE InpOutp, ONLY: XS_updt, ounit

  IMPLICIT NONE

  INTEGER, INTENT(IN), OPTIONAL :: ind    ! if iteration reaching th_iter and ind = 0 then STOP
  REAL, DIMENSION(nnod) :: pline
  REAL, DIMENSION(nnod) :: otem
  INTEGER :: n, niter

  th_err = 1.
  niter = 0
  DO
      niter = niter + 1

      ! Save old moderator temp
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
      CALL th_upd2(pline)

      th_err = MAXVAL(ABS(ftem - otem))
      IF ((th_err < 0.01) .OR. (niter == th_niter)) EXIT
  END DO

  IF (PRESENT(ind)) THEN
    IF ((niter == th_niter) .AND. (ind == 0)) THEN
       WRITE(ounit,*) '  MAXIMUM TH ITERATION REACHED.'
       WRITE(ounit,*) '  CALCULATION MIGHT BE NOT CONVERGED OR CHANGE ITERATION CONTROL'
       STOP
    END IF
  END IF



END SUBROUTINE th_iter


SUBROUTINE par_ave_f(par, ave)
!
! Purpose:
!    To calculate average fuel temp (only for active core)
!

USE sdata, ONLY: vdel, nnod, ystag, xyz, nzz, nyy

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: par
REAL, INTENT(OUT) :: ave
REAL :: dum, dum2
INTEGER :: i, j, k

dum = 0.; dum2 = 0.
DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax
            dum = dum + par(xyz(i,j,k)) * vdel(xyz(i,j,k))
            dum2 = dum2 + vdel(xyz(i,j,k))
        END DO
    END DO
END DO

ave = dum / dum2

END SUBROUTINE par_ave_f


SUBROUTINE par_ave(par, ave)
!
! Purpose:
!    To calculate average moderator temp (only for radially active core)
!

USE sdata, ONLY: vdel, nnod, ystag, xyz, nzz, nyy

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: par
REAL, INTENT(OUT) :: ave
REAL :: dum, dum2
INTEGER :: i, j, k

dum = 0.; dum2 = 0.
DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax
            dum = dum + par(xyz(i,j,k)) * vdel(xyz(i,j,k))
            dum2 = dum2 + vdel(xyz(i,j,k))
        END DO
    END DO
END DO

ave = dum / dum2

END SUBROUTINE par_ave


SUBROUTINE par_max(par, pmax, ti, tj, tk)
!
! Purpose:
!    To calculate maximum fuel tem, coolant tem, and density
!

USE sdata, ONLY: nnod, ystag, xyz, nzz, nyy

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: par
REAL, INTENT(OUT) :: pmax
INTEGER, INTENT(OUT) :: ti, tj, tk
INTEGER :: i, j, k

pmax = 0.
DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax
            IF (par(xyz(i,j,k)) > pmax) THEN
                pmax = par(xyz(i,j,k))
                ti = i
                tj = j
                tk = k
            END IF
        END DO
    END DO
END DO

END SUBROUTINE par_max


SUBROUTINE getent(t,ent)
!
! Purpose:
!    To get enthalpy for given coolant temp. from steam table
!

USE sdata, ONLY: thunit
USE InpOutp, ONLY : ounit

IMPLICIT NONE

REAL, INTENT(IN) :: t
REAL, INTENT(OUT) :: ent
REAL :: t1, rho1, ent1
REAL :: t2, rho2, ent2

IF ((t < 473.15) .OR. (t > 617.91)) THEN
    WRITE(ounit,*) '  Coolant temp. : ', t
    WRITE(ounit,*) '  ERROR : MODERATOR TEMP. IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
    WRITE(ounit,*) '  CHECK INPUT MASS FLOW RATE OR POWER'
    STOP
END IF

REWIND(thunit)

READ(thunit,*) t2, rho2, ent2
DO
    t1 = t2
    ent1 = ent2
    READ(thunit,*) t2, rho2, ent2
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

USE sdata, ONLY: thunit
USE InpOutp, ONLY : ounit

IMPLICIT NONE

REAL, INTENT(IN) :: ent
REAL, INTENT(OUT) :: t, rho, prx, kvx, tcx
REAL :: t1, rho1, ent1, kv1, pr1, tc1
REAL :: t2, rho2, ent2, kv2, pr2, tc2
REAL :: ratx

IF ((ent < 858341.5) .OR. (ent > 1624307.1)) THEN
    WRITE(ounit,*) '  Enthalpy. : ', ent
    WRITE(ounit,*) '  ERROR : ENTHALPY IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
    WRITE(ounit,*) '  CHECK INPUT MASS FLOW RATE OR POWER'
    STOP
END IF

REWIND(thunit)

READ(thunit,*) t2, rho2, ent2, pr2, kv2, tc2
DO
    t1 = t2
    ent1 = ent2
    rho1 = rho2
    pr1 = pr2
    kv1 = kv2
    tc1 = tc2
    READ(thunit,*) t2, rho2, ent2, pr2, kv2, tc2
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


REAL FUNCTION getcpc(t)
!
! Purpose:
!    To calculate specific heat capacity of cladding
!

IMPLICIT NONE

REAL, INTENT(IN) :: t

getcpc = 252.54 + 0.11474*t

END FUNCTION getcpc


REAL FUNCTION getcpf(t)
!
! Purpose:
!    To calculate specific heat capacity of fuel
!

IMPLICIT NONE

REAL, INTENT(IN) :: t

getcpf = 162.3 + 0.3038*t - 2.391e-4*t**2 + 6.404e-8*t**3

END FUNCTION getcpf


SUBROUTINE TridiaSolve(a,b,c,d,x)
!
! Purpose:
!    To solve tridiagonal matrix
!

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(INOUT) :: a, b, c, d
REAL, DIMENSION(:), INTENT(OUT) :: x

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



REAL FUNCTION geths(xden)
!
! Purpose:
!    To calculate heat transfer coef.
!

USE sdata, ONLY: dh, farea, cflow

IMPLICIT NONE

REAL, INTENT(IN) :: xden

REAL :: cvelo, Nu, Re

cvelo = cflow / (farea * xden * 1000.)        ! Calculate flow velocity (m/s)
Re = cvelo * dh / (kv * 1e-6)                 ! Calculate Reynolds Number
Nu = 0.023*(Pr**0.4)*(Re**0.8)                ! Calculate Nusselt Number
geths = (tc / dh) * Nu                        ! Calculate heat transfer coefficient


END FUNCTION geths


SUBROUTINE th_trans4(xpline, h)

!
! Purpose:
!    To perform fuel pin thermal transient
!

USE sdata, ONLY: mtem, cden, ftem, tin, xyz, cflow, nyy, nzz, cf, ent, heatf, nnod, &
                 ystag, tfm, nt, nf, rpos, rdel, tg, rf, farea, dia, pi, node_nf, zdel, ystag
USE InpOutp, ONLY : ounit

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: xpline    ! Linear Power Density (W/cm)
REAL, INTENT(IN) :: h                       ! Time step

INTEGER :: i, j, k, n
REAL :: hs, hg = 1.e4, kt           ! coolant heat transfer coef., gap heat transfer coef, and thermal conductivity
REAL :: alpha = 0.7
REAL :: xa, xc, tem, tem2
REAL :: pdens                       ! Fuel pin power density channel
REAL :: fdens = 10.412e3            ! UO2 density (kg/m3)
REAL :: cdens = 6.6e3               ! Cladding density (kg/m3)
REAL :: cp                          ! Specific heat capacity
REAL :: eta, alp, beta
REAL :: mdens, cpline, vol
REAL :: enti
REAL, DIMENSION(nnod) :: entp        ! previous enthalpy

CALL getent(tin, enti)
entp = ent

DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax

            mdens = cden(xyz(i,j,k)) * 1000.                                    ! Coolant density (kg/m3)
            cpline = heatf(xyz(i,j,k)) * pi * dia  &
                  + cf * xpline(xyz(i,j,k)) * 100.                              ! Coolant Linear power densisty (W/m)
            vol   = farea * zdel(k) * 0.01
            eta = 0.5 * mdens * vol / h
            IF (k == 1) THEN                                                    ! Calculate coolant enthalpy
                ent(xyz(i,j,k)) = (cpline * zdel(k) * 0.01 &
                                - cflow * (ent(xyz(i,j,k)) - enti) &
                                + eta * ent(xyz(i,j,k))) / eta
                CALL gettd(0.5 * (enti + ent(xyz(i,j,k))), mtem(xyz(i,j,k)), &
                          cden(xyz(i,j,k)), Pr, kv, tc)                             ! Get corresponding temp and density
            ELSE
              ent(xyz(i,j,k)) = (cpline * zdel(k) * 0.01 &
                              - cflow * (ent(xyz(i,j,k)) - entp(xyz(i,j,k-1))) &
                              + eta * (ent(xyz(i,j,k)) + entp(xyz(i,j,k-1)) &
                              - ent(xyz(i,j,k-1)))) / eta
                CALL gettd(0.5 * (ent(xyz(i,j,k-1)) + ent(xyz(i,j,k))), &
                           mtem(xyz(i,j,k)), cden(xyz(i,j,k)), Pr, kv, tc)          ! Get corresponding temp and density
            END IF


            hs = geths(cden(xyz(i,j,k)))                                               ! Calculate heat transfer coef
            pdens = (1. - cf) * 100. * xpline(xyz(i,j,k)) / (pi * rf**2)                ! Fuel pin Power Density (W/m3)

            ! Calculate tridiagonal matrix: a, b, c and source: d
            ! For nt=1 [FUEL CENTERLINE]
            tem = 0.5 * (tfm(xyz(i,j,k),1) + tfm(xyz(i,j,k),2))                        ! Average temp. to get thermal conductivity
            kt = getkf(tem)                                                            ! Get thermal conductivity
            cp = getcpf(tfm(xyz(i,j,k),1))                                                           ! Get specific heat capacity
            eta = fdens * cp * rpos(1)**2 / (2. * h)
            xc  = kt * rpos(1) / rdel(1)
            tfm(xyz(i,j,k),1) = pdens * 0.5 * rpos(1)**2 / eta &
                              + xc * tfm(xyz(i,j,k),2)  / eta &
                              - xc * tfm(xyz(i,j,k),1)  / eta &
                              + tfm(xyz(i,j,k),1)

            DO n = 2, nt-2
                tem = 0.5 * (tfm(xyz(i,j,k),n) + tfm(xyz(i,j,k),n+1))
                kt = getkf(tem)
                cp = getcpf(tfm(xyz(i,j,k),n))
                eta = fdens * cp * (rpos(n)**2 - rpos(n-1)**2) / (2. * h)
                xa = xc
                xc = kt * rpos(n) / rdel(n)
                tfm(xyz(i,j,k),n) = pdens * 0.5 * (rpos(n)**2 - rpos(n-1)**2) / eta &
                                  + xa * tfm(xyz(i,j,k),n-1)  / eta &
                                  + xc * tfm(xyz(i,j,k),n+1)  / eta &
                                  - (xa + xc) * tfm(xyz(i,j,k),n)  / eta &
                                  + tfm(xyz(i,j,k),n)
            END DO

            ! For nt-1 [FUEL SURFACE]
            cp = getcpf(tfm(xyz(i,j,k),nt-1))
            eta = fdens * cp * (rf**2 - rpos(nt-2)**2) / (2. * h)
            xa = xc
            xc = hg * rpos(nt)   ! This is position of inner clad
            tfm(xyz(i,j,k),nt-1) = pdens * 0.5 * (rf**2 - rpos(nt-2)**2) / eta &
                                 + xa * tfm(xyz(i,j,k),nt-2)  / eta &
                                 + xc * tfm(xyz(i,j,k),nt)  / eta &
                                 - (xa + xc) * tfm(xyz(i,j,k),nt-1)  / eta &
                                 + tfm(xyz(i,j,k),nt-1)

            ! For nt [INNER CLAD]
            tem = 0.5 * (tfm(xyz(i,j,k),nt) + tfm(xyz(i,j,k),nt+1))
            kt = getkc(tem)      ! For cladding
            cp = getcpc(tfm(xyz(i,j,k),nt))
            eta = cdens * cp * (rpos(nt+1)**2 - rpos(nt)**2) / (2. * h)
            xa = xc
            xc = kt * rpos(nt+1) / rdel(nt)
            tfm(xyz(i,j,k),nt) = xa * tfm(xyz(i,j,k),nt-1)  / eta &
                               + xc * tfm(xyz(i,j,k),nt+1)  / eta &
                               - (xa + xc) * tfm(xyz(i,j,k),nt)  / eta &
                               + tfm(xyz(i,j,k),nt)

            ! For nt+1  [OUTER CLAD]
            cp = getcpc(tfm(xyz(i,j,k),nt+1))
            eta = cdens * cp * (rpos(nt+2)**2 - rpos(nt+1)**2) / (2. * h)
            xa = xc
            xc = hs * rpos(nt+2)
            tfm(xyz(i,j,k),nt+1) = hs * rpos(nt+2) * mtem(xyz(i,j,k)) / eta &
                                 + xa * tfm(xyz(i,j,k),nt)  / eta &
                                 - (xa + xc) * tfm(xyz(i,j,k),nt+1)  / eta &
                                 + tfm(xyz(i,j,k),nt+1)

            ! Get lumped fuel temp
            ftem(xyz(i,j,k)) = (1.-alpha) * tfm(xyz(i,j,k), 1) &
                             + alpha * tfm(xyz(i,j,k), nt-1)

            ! Calculate heat flux
            heatf(xyz(i,j,k)) = hs * (tfm(xyz(i,j,k), nt+1) - mtem(xyz(i,j,k)))

        END DO
    END DO
END DO

END SUBROUTINE th_trans4


SUBROUTINE th_trans3(xpline, h)

!
! Purpose:
!    To perform fuel pin thermal transient
!

USE sdata, ONLY: mtem, cden, ftem, tin, xyz, cflow, nyy, nzz, cf, ent, heatf, nnod, &
                 ystag, tfm, nt, nf, rpos, rdel, tg, rf, farea, dia, pi, node_nf, zdel, ystag
USE InpOutp, ONLY : ounit

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: xpline    ! Linear Power Density (W/cm)
REAL, INTENT(IN) :: h                       ! Time step

INTEGER :: i, j, k, n
REAL, DIMENSION(nt+1) :: a, b, c, d
REAL :: hs, hg = 1.e4, kt           ! coolant heat transfer coef., gap heat transfer coef, and thermal conductivity
REAL :: alpha = 0.7
REAL :: xa, xc, tem
REAL :: pdens                       ! Fuel pin power density channel
REAL :: fdens = 10.412e3            ! UO2 density (kg/m3)
REAL :: cdens = 6.6e3               ! Cladding density (kg/m3)
REAL :: cp                          ! Specific heat capacity
REAL :: eta, alp, beta
REAL :: mdens, cpline, vol
REAL :: enti
REAL, DIMENSION(nnod) :: entp        ! previous enthalpy

CALL getent(tin, enti)
entp = ent

DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax

            mdens = cden(xyz(i,j,k)) * 1000.                                    ! Coolant density (kg/m3)
            cpline = heatf(xyz(i,j,k)) * pi * dia  &
                  + cf * xpline(xyz(i,j,k)) * 100.                              ! Coolant Linear power densisty (W/m)
            vol   = farea * zdel(k) * 0.01
            IF (k == 1) THEN                                                    ! Calculate coolant enthalpy
                eta = enti + entp(xyz(i,j,k))
                ent(xyz(i,j,k)) = (cpline * zdel(k) * 0.01 * h &
                                + (cflow * h - 0.5 * mdens * vol) * enti &
                                + 0.5 * mdens * vol * eta) &
                                / (0.5 * mdens * vol + cflow * h)
                CALL gettd(0.5 * (enti + ent(xyz(i,j,k))), mtem(xyz(i,j,k)), &
                          cden(xyz(i,j,k)), Pr, kv, tc)                             ! Get corresponding temp and density
            ELSE
                eta = entp(xyz(i,j,k-1)) + entp(xyz(i,j,k))
                ent(xyz(i,j,k)) = (cpline * zdel(k) * 0.01 * h &
                                + (cflow * h - 0.5 * mdens * vol) * ent(xyz(i,j,k-1)) &
                                + 0.5 * mdens * vol * eta) &
                                / (0.5 * mdens * vol + cflow * h)
                CALL gettd(0.5 * (ent(xyz(i,j,k-1)) + ent(xyz(i,j,k))), &
                           mtem(xyz(i,j,k)), cden(xyz(i,j,k)), Pr, kv, tc)          ! Get corresponding temp and density
            END IF


            hs = geths(cden(xyz(i,j,k)))                                               ! Calculate heat transfer coef
            pdens = (1. - cf) * 100. * xpline(xyz(i,j,k)) / (pi * rf**2)                ! Fuel pin Power Density (W/m3)

            ! Calculate tridiagonal matrix: a, b, c and source: d
            ! For nt=1 [FUEL CENTERLINE]
            tem = 0.5 * (tfm(xyz(i,j,k),1) + tfm(xyz(i,j,k),2))                        ! Average temp. to get thermal conductivity
            kt = getkf(tem)                                                            ! Get thermal conductivity
            cp = getcpf(tfm(xyz(i,j,k),1))                                                           ! Get specific heat capacity
            eta = fdens * cp * rpos(1)**2 / (2. * h)
            xc  = kt * rpos(1) / rdel(1)
            b(1) =  xc + eta
            c(1) = -xc
            d(1) = pdens * 0.5 * rpos(1)**2 + eta * tfm(xyz(i,j,k),1)

            DO n = 2, nt-2
                tem = 0.5 * (tfm(xyz(i,j,k),n) + tfm(xyz(i,j,k),n+1))
                kt = getkf(tem)
                cp = getcpf(tfm(xyz(i,j,k),n))
                eta = fdens * cp * (rpos(n)**2 - rpos(n-1)**2) / (2. * h)
                xa = xc
                xc = kt * rpos(n) / rdel(n)
                a(n) = -xa
                b(n) =  xa + xc + eta
                c(n) = -xc
                d(n) = pdens * 0.5 * (rpos(n)**2 - rpos(n-1)**2) &
                     + eta * tfm(xyz(i,j,k),n)
            END DO

            ! For nt-1 [FUEL SURFACE]
            cp = getcpf(tfm(xyz(i,j,k),nt-1))
            eta = fdens * cp * (rf**2 - rpos(nt-2)**2) / (2. * h)
            xa = xc
            xc = hg * rpos(nt)   ! This is position of inner clad
            a(nt-1) = -xa
            b(nt-1) =  xa + xc + eta
            c(nt-1) = -xc
            d(nt-1) = pdens * 0.5 * (rf**2 - rpos(nt-2)**2) &
                    + eta * tfm(xyz(i,j,k),nt-1)

            ! For nt [INNER CLAD]
            tem = 0.5 * (tfm(xyz(i,j,k),nt) + tfm(xyz(i,j,k),nt+1))
            kt = getkc(tem)      ! For cladding
            cp = getcpc(tfm(xyz(i,j,k),nt))
            eta = cdens * cp * (rpos(nt+1)**2 - rpos(nt)**2) / (2. * h)
            xa = xc
            xc = kt * rpos(nt+1) / rdel(nt)
            a(nt) = -xa
            b(nt) =  xa + xc + eta
            c(nt) = -xc
            d(nt) = eta * tfm(xyz(i,j,k),nt)

            ! For nt+1  [OUTER CLAD]
            cp = getcpc(tfm(xyz(i,j,k),nt+1))
            eta = cdens * cp * (rpos(nt+2)**2 - rpos(nt+1)**2) / (2. * h)
            xa = xc
            xc = hs * rpos(nt+2)
            a(nt+1) = -xa
            b(nt+1) =  xa + xc + eta
            d(nt+1) = hs * rpos(nt+2) * mtem(xyz(i,j,k)) &
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

END SUBROUTINE th_trans3



SUBROUTINE th_upd2(xpline)

!
! Purpose:
!    To update thermal parameters
!

USE sdata, ONLY: mtem, cden, ftem, tin, xyz, cflow, nyy, nzz, cf, ent, heatf, &
                 ystag, tfm, nt, nf, rpos, rdel, tg, rf, pi, node_nf, zdel, dia, ystag
USE InpOutp, ONLY : ounit

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: xpline    ! Linear Power Density (W/cm)

INTEGER :: i, j, k, n
REAL :: cp
REAL, DIMENSION(nt+1) :: a, b, c, d
REAL :: hs, hg = 1.e4, kt
REAL :: alp = 0.7
REAL :: xa, xc, tem
REAL :: pdens
REAL :: enti
REAL :: cpline

CALL getent(tin, enti)

DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax
        !WRITE(ounit,*) cflow , node_nf(i,j), i, j
        !IF (node_nf(i,j) < 2) STOP

            cpline = heatf(xyz(i,j,k)) * pi * dia  &
                   + cf * xpline(xyz(i,j,k)) * 100.                             ! Coolant Linear power densisty (W/m)
            !WRITE(*,*) i, j, k, cpline/100.
            IF (k == 1) THEN                                                    ! Calculate coolant enthalpy and
                ent(xyz(i,j,k)) = enti + cpline * zdel(k) * 0.01 / cflow        ! corresponding temp and density
                CALL gettd(0.5 * (enti + ent(xyz(i,j,k))), &
                           mtem(xyz(i,j,k)), cden(xyz(i,j,k)), Pr, kv, tc)          ! Get corresponding temp and density
            ELSE
                ent(xyz(i,j,k)) = ent(xyz(i,j,k-1)) &
                                + cpline * zdel(k) * 0.01 / cflow
                CALL gettd(0.5 * (ent(xyz(i,j,k-1)) + ent(xyz(i,j,k))), &
                          mtem(xyz(i,j,k)), cden(xyz(i,j,k)), Pr, kv, tc)           ! Get corresponding temp and density
            END IF

            hs = geths(cden(xyz(i,j,k)))
            pdens = (1. - cf) * 100. * xpline(xyz(i,j,k)) / (pi * rf**2)        ! Fuel pin Power Density (W/m3)

            ! Calculate tridiagonal matrix: a, b, c and source: d
            tem = 0.5 * (tfm(xyz(i,j,k),1) + tfm(xyz(i,j,k),2))                 ! Average temp. to get thermal conductivity
            kt = getkf(tem)                                                     ! Get thermal conductivity
            xc  = kt * rpos(1) / rdel(1)
            b(1) =  xc
            c(1) = -xc
            d(1) = pdens * 0.5 * rpos(1)**2

            DO n = 2, nt-2
                tem = 0.5 * (tfm(xyz(i,j,k),n) + tfm(xyz(i,j,k),n+1))
                kt = getkf(tem)
                xa = xc
                xc = kt * rpos(n) / rdel(n)
                a(n) = -xa
                b(n) =  xa + xc
                c(n) = -xc
                d(n) = pdens * 0.5 * (rpos(n)**2 - rpos(n-1)**2)
            END DO

            ! For nt-1 [FUEL SURFACE]
            xa = xc
            xc = hg * rpos(nt)   ! This is position of inner clad
            a(nt-1) = -xa
            b(nt-1) =  xa + xc
            c(nt-1) = -xc
            d(nt-1) = pdens * 0.5 * (rf**2 - rpos(nt-2)**2)

            ! For nt [INNER CLAD]
            tem = 0.5 * (tfm(xyz(i,j,k),nt) + tfm(xyz(i,j,k),nt+1))
            kt = getkc(tem)      ! For cladding
            xa = xc
            xc = kt * rpos(nt+1) / rdel(nt)
            a(nt) = -xa
            b(nt) =  xa + xc
            c(nt) = -xc
            d(nt) = 0.

            ! For nt+1  [OUTER CLAD]
            xa = xc
            a(nt+1) = -xa
            b(nt+1) =  xa + hs * rpos(nt+2)
            d(nt+1) = hs * rpos(nt+2) * mtem(xyz(i,j,k))

            ! Solve tridiagonal matrix
            CALL TridiaSolve(a, b, c, d, tfm(xyz(i,j,k), :))

            ! Get lumped fuel temp
            ftem(xyz(i,j,k)) = (1.-alp) * tfm(xyz(i,j,k), 1) + alp * tfm(xyz(i,j,k), nt-1)

            ! Calculate heat flux
            heatf(xyz(i,j,k)) = hs * (tfm(xyz(i,j,k), nt+1) - mtem(xyz(i,j,k)))


        END DO
    END DO


!STOP
END DO


END SUBROUTINE th_upd2


SUBROUTINE cbsearch()

!
! Purpose:
!    To search critical boron concentration
!

USE sdata, ONLY: Ke, rbcon, ftem, mtem, cden, bpos, nnod, ng, f0, fer, ser, &
                 aprad, apaxi, afrad, npow
USE InpOutp, ONLY: ounit, XS_updt, AsmFlux, AsmPow, AxiPow
USE nodal, ONLY: nodal_coup4, outer4, powdis

IMPLICIT NONE

INTEGER :: rKe  ! Rounded Keff
REAL  :: bc, bc1, bc2     ! Boron Concentration
REAL :: ke1, ke2
INTEGER :: n

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
               //'   FISS. SOURCE REL. ERROR'
WRITE(ounit,*) ' -----------------------------------------------------------' &
              // '--------------------------'

CALL XS_updt(rbcon, ftem, mtem, cden, bpos)
CALL nodal_coup4()
CALL outer4(0)
rKe = NINT(Ke * 100000)
bc1 = rbcon
ke1 = Ke

WRITE(ounit,'(I5, F15.2, F23.5, ES15.5, ES17.5)') 1, bc1, Ke1, ser, fer

bc2 = rbcon + (Ke - 1.) * rbcon
CALL XS_updt(bc2, ftem, mtem, cden, bpos)
CALL nodal_coup4()
CALL outer4(0)
rKe = NINT(Ke * 100000)
ke2 = Ke

WRITE(ounit,'(I5, F15.2, F23.5, ES15.5, ES17.5)') 2, bc2, Ke2, ser, fer

n = 3
DO
	bc = bc2 + (1.0 - ke2) / (ke1 - ke2) * (bc1 - bc2)
	CALL XS_updt(bc, ftem, mtem, cden, bpos)
    CALL nodal_coup4()
    CALL outer4(0)
	rKe = NINT(Ke * 100000)
	bc1 = bc2
	bc2 = bc
	ke1 = ke2
	ke2 = ke
    WRITE(ounit,'(I5, F15.1, F23.5, ES15.5, ES17.5)') n, bc, Ke, ser, fer
	IF ((rKe == 100000) .AND. (ser < 1.e-2) .AND. (fer < 1.e-2)) EXIT
	n = n + 1
	IF (bc > 2999. .AND. bc < 3000.) THEN
	    WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
		WRITE(ounit,*) '  ADPRES IS STOPPING'
	    STOP
	END IF
	IF (bc > 0. .AND. bc < 1.) THEN
	    WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
		WRITE(ounit,*) '  ADPRES IS STOPPING'
	    STOP
	END IF
	IF (n == 30) THEN
	    WRITE(ounit,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
		WRITE(ounit,*) '  ADPRES IS STOPPING'
	    STOP
	END IF
END DO

IF (aprad == 1 .OR. apaxi == 1) THEN
    ALLOCATE(npow(nnod))
	CALL PowDis(npow)
END IF

IF (aprad == 1) CALL AsmPow(npow)

IF (apaxi == 1) CALL AxiPow(npow)

IF (afrad == 1) CALL AsmFlux(f0, 1.e0)


END SUBROUTINE cbsearch


SUBROUTINE cbsearcht()

!
! Purpose:
!    To search critical boron concentration with thermal feedback
!

USE sdata, ONLY: Ke, ftem, mtem, cden, bpos, bcon, rbcon, npow, nnod, &
                 f0, ser, fer, tfm, aprad, apaxi, afrad, npow, th_err
USE InpOutp, ONLY: ounit, XS_updt, AsmFlux, AsmPow, AxiPow, getfq
USE nodal, ONLY: powdis, nodal_coup4, outer4

IMPLICIT NONE

INTEGER :: rKe  ! Rounded Keff
REAL  :: bc1, bc2    ! Boron Concentration
REAL :: ke1, ke2
INTEGER :: n
REAL :: tf, tm, maxtf, maxtm, maxfcl, cd, fz

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

ALLOCATE(npow(nnod))

bcon = rbcon
CALL th_iter()  ! Start thermal hydarulic iteration with current paramters
rKe = NINT(Ke * 100000)
bc1 = bcon
ke1 = Ke

WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') 1, bc1, Ke1, ser, fer, th_err

bcon = bcon + (Ke - 1.) * bcon   ! Guess next critical boron concentration
CALL th_iter()                 ! Perform second thermal hydarulic iteration with updated parameters
bc2 = bcon
ke2 = Ke

WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') 2, bc2, Ke2, ser, fer, th_err

n = 3
DO
	bcon = bc2 + (1.0 - ke2) / (ke1 - ke2) * (bc1 - bc2)
	IF ((rKe > 99900) .AND. (rKe < 100100)) THEN
	    CALL th_iter()
	ELSE
	    CALL th_iter()
	END IF
	rKe = NINT(Ke * 100000)
	bc1 = bc2
	bc2 = bcon
	ke1 = ke2
	ke2 = ke
    WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') n, bcon, Ke, ser, fer, th_err
	IF ((rKe == 100000) .AND. (ser < 1.e-5) .AND. (fer < 1.e-5)) EXIT
	n = n + 1
	IF (bcon > 2999. .AND. bcon < 3000.) THEN
	    WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
		WRITE(ounit,*) '  ADPRES IS STOPPING'
	    STOP
	END IF
	IF (bcon > 0. .AND. bcon < 1.) THEN
	    WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
		WRITE(ounit,*) '  ADPRES IS STOPPING'
	    STOP
	END IF
	IF (n == 30) THEN
	    WRITE(ounit,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
		WRITE(ounit,*) '  ADPRES IS STOPPING'
	    STOP
	END IF
END DO

IF (aprad == 1 .OR. apaxi == 1) THEN
	CALL PowDis(npow)
END IF

IF (aprad == 1) CALL AsmPow(npow)

IF (apaxi == 1) CALL AxiPow(npow)

IF (afrad == 1) CALL AsmFlux(f0, 1.e0)

! CALL par_ave_f(ftem, tf)
! CALL par_ave(mtem, tm)

! CALL par_max(ftem, maxtf)
! CALL par_max(tfm(:,1), maxfcl)
! CALL par_max(mtem, maxtm)
! CALL getfq(npow)

! Write Output
! WRITE(ounit,*)
! WRITE(ounit, 5001) tf, tf-273.15
! WRITE(ounit, 5002)  maxfcl, maxfcl-273.15
! WRITE(ounit, 5003) tm, tm-273.15
! WRITE(ounit, 5004) maxtm, maxtm-273.15

5001 FORMAT(2X, 'AVERAGE DOPPLER TEMPERATURE     : ', F7.1, ' K (', F7.1, ' C)')
5002 FORMAT(2X, 'MAX FUEL CENTERLINE TEMPERATURE : ', F7.1, ' K (', F7.1, ' C)')
5003 FORMAT(2X, 'AVERAGE MODERATOR TEMPERATURE   : ', F7.1, ' K (', F7.1, ' C)')
5004 FORMAT(2X, 'MAXIMUM MODERATOR TEMPERATURE   : ', F7.1, ' K (', F7.1, ' C)')


END SUBROUTINE cbsearcht


END MODULE th
