MODULE nodal

!=========================
! Nodal Module to solve diffusion equation using NEM
! For the theoretical background of NEM used here refer to:
! Okumura K (1998) MOSRA-light: high speed three-dimensional nodal diffusion code for vector computers,
! JAERI-Data/Code 98-025. Japan Atomic Energy Research Institute, Tokaimura (in Japanese)
! =======================

IMPLICIT NONE

SAVE

CONTAINS


SUBROUTINE outer4(popt)

!
! Purpose:
!    To perform normal outer iteration


USE sdata, ONLY: ng, nnod, nout, serc, ferc, fer, ser, &
                 Ke, nac, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2
USE InpOutp, ONLY: ounit

IMPLICIT NONE

INTEGER, OPTIONAL, INTENT(IN) :: popt

DOUBLE PRECISION :: Keo                                    !Old Multiplication factor (Keff)
DOUBLE PRECISION, DIMENSION(nnod) :: fs0c                  !old fission source
DOUBLE PRECISION, DIMENSION(nnod,ng,6) :: joc              !Old outgoing currents
DOUBLE PRECISION, DIMENSION(nnod) :: fsx1c, fsy1c, fsz1c
DOUBLE PRECISION, DIMENSION(nnod) :: fsx2c, fsy2c, fsz2c
DOUBLE PRECISION, DIMENSION(nnod) :: ss0                   ! Scattering source
DOUBLE PRECISION, DIMENSION(nnod) :: ssx1, ssy1, ssz1
DOUBLE PRECISION, DIMENSION(nnod) :: ssx2, ssy2, ssz2      ! Scattering source moments
DOUBLE PRECISION :: f, fc                                  ! new and old integrated fission sources
DOUBLE PRECISION :: domiR, e1, e2
INTEGER :: g
INTEGER :: p, npos
LOGICAL :: opt

DOUBLE PRECISION, DIMENSION(nnod) :: errn, erro

IF (PRESENT(popt)) THEN
    opt = .FALSE.
ELSE
    opt = .TRUE.
END IF

! Initialize fission source
fs0  = 0.0
fsx1 = 0.0; fsy1 = 0.0; fsz1 = 0.0
fsx2 = 0.0; fsy2 = 0.0; fsz2 = 0.0

DO g= 1, ng
    CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
END DO

errn = 1.0
f = Integrate(fs0)
e1 = Integrate(errn)

!Start outer iteration
DO p=1, nout
    CALL SaveJo(joc)  ! Save old outgoing currents
    fc = f            ! Save old integrated fission source
    fs0c  = fs0       ! Save old fission source
    fsx1c = fsx1; fsy1c = fsy1; fsz1c = fsz1
    fsx2c = fsx2; fsy2c = fsy2; fsz2c = fsz2
    fs0  = 0.0
    fsx1 = 0.0; fsy1 = 0.0; fsz1 = 0.0
    fsx2 = 0.0; fsy2 = 0.0; fsz2 = 0.0
    Keo = Ke          ! Save old multiplication factor
    erro = errn       ! Save old fission source error/difference
    DO g = 1, ng

        ss0  = 0.0
        ssx1 = 0.0; ssy1 = 0.0; ssz1 = 0.0
        ssx2 = 0.0; ssy2 = 0.0; ssz2 = 0.0

        !!!Calculate Scattering source
        CALL SSrc(g, ss0, ssx1, ssy1, ssz1, ssx2, ssy2, ssz2)

        !!!Calculate total source
        CALL TSrc(g, Keo, fs0c , fsx1c, fsy1c, fsz1c, &
                                 fsx2c, fsy2c, fsz2c, &
                           ss0 , ssx1 , ssy1 , ssz1 , &
                                 ssx2 , ssy2 , ssz2   )

        !!!Inner Iteration
        CALL inner4(g)

        !!!Calculate fission source for next outer iteration
        CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
    END DO

    errn = fs0 - fs0c
    e2 = Integrate(ABS(errn))

    IF (MOD(p,nac) == 0) THEN   ! Fission source extrapolation
        domiR = e2 / e1
        npos = MAXLOC(ABS(erro),1)
        IF (erro(npos) * errn(npos) < 0.0) domiR = -domiR
        fs0 = fs0 + domiR / (1.0 - domiR) * errn
    END IF
    e1 = e2                       ! Save integrated fission source error
    f = Integrate(fs0)            ! Integrate fission source
    Ke = Keo * f / fc             ! Update Keff
    CALL RelE(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
    CALL RelJ(joc, fer)           ! Search maximum point wise outgoing current Error
    IF ((ser < serc) .AND. (fer < ferc)) EXIT  ! If converge, exit.
END DO

IF (p-1 == nout) THEN
    WRITE(ounit,*)
    WRITE(ounit,*) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED.'
    WRITE(ounit,*) '  ITERATION MAY NOT CONVERGE.'
    WRITE(ounit,*) '  CHECK PROBLEM SPECIFICATION OR CHANGE ITERATION CONTROL (%ITER).'
    WRITE(ounit,*) '  ADPRES IS STOPING...'
    STOP
END IF

IF (opt) WRITE(ounit,*)
IF (opt) WRITE(ounit, '(A36,F9.6)') 'MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke


END SUBROUTINE outer4


SUBROUTINE outer4th(maxn)

!
! Purpose:
!    To perform normal outer iteration when %THER card present


USE sdata, ONLY: ng, nnod, serc, ferc, fer, ser,  &
                 Ke, nac, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2

IMPLICIT NONE

INTEGER, INTENT(IN) :: maxn

DOUBLE PRECISION :: Keo                                    !Old Multiplication factor (Keff)
DOUBLE PRECISION, DIMENSION(nnod) :: fs0c                  !Old fission source
DOUBLE PRECISION, DIMENSION(nnod,ng,6) :: joc              !Old outgoing currents
DOUBLE PRECISION, DIMENSION(nnod) :: fsx1c, fsy1c, fsz1c
DOUBLE PRECISION, DIMENSION(nnod) :: fsx2c, fsy2c, fsz2c
DOUBLE PRECISION, DIMENSION(nnod) :: ss0                   ! Scattering source
DOUBLE PRECISION, DIMENSION(nnod) :: ssx1, ssy1, ssz1
DOUBLE PRECISION, DIMENSION(nnod) :: ssx2, ssy2, ssz2      ! Scattering source moments
DOUBLE PRECISION :: f, fc                                  ! new and old integrated fission sources
DOUBLE PRECISION :: domiR, e1, e2
INTEGER :: g, i, n
INTEGER :: p, npos

DOUBLE PRECISION, DIMENSION(nnod) :: errn, erro


! Initialize fission source
fs0  = 0.0
fsx1 = 0.0; fsy1 = 0.0; fsz1 = 0.0
fsx2 = 0.0; fsy2 = 0.0; fsz2 = 0.0

DO g= 1, ng
    CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
END DO

errn = 1.0
f = Integrate(fs0)
e1 = Integrate(errn)

!Start outer iteration
DO p=1, maxn
    CALL SaveJo(joc)  ! Save old outgoing currents
    fc = f            ! Save old integrated fission source
    fs0c  = fs0       ! Save old fission source
    fsx1c = fsx1; fsy1c = fsy1; fsz1c = fsz1
    fsx2c = fsx2; fsy2c = fsy2; fsz2c = fsz2
    fs0  = 0.0
    fsx1 = 0.0; fsy1 = 0.0; fsz1 = 0.0
    fsx2 = 0.0; fsy2 = 0.0; fsz2 = 0.0
    Keo = Ke          ! Save old multiplication factor
    erro = errn       ! Save old fission source error/difference
    DO g = 1, ng

        !!!Calculate Scattering source
        CALL SSrc(g, ss0, ssx1, ssy1, ssz1, ssx2, ssy2, ssz2)

        !!!Calculate total source
        CALL TSrc(g, Keo, fs0c , fsx1c, fsy1c, fsz1c, &
                                 fsx2c, fsy2c, fsz2c, &
                           ss0 , ssx1 , ssy1 , ssz1 , &
                                 ssx2 , ssy2 , ssz2   )

        !!!Inner Iteration
        CALL inner4(g)

        !!!Calculate fission source for next outer iteration
        CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
    END DO

    errn = fs0 - fs0c
    e2 = Integrate(ABS(errn))

    IF (MOD(p,nac) == 0) THEN   ! Fission source extrapolation
        domiR = e2 / e1
        npos = MAXLOC(ABS(erro),1)
        IF (erro(npos) * errn(npos) < 0.0) domiR = -domiR
        fs0 = fs0 + domiR / (1.0 - domiR) * errn
    END IF
    e1 = e2                       ! Save integrated fission source error
    f = Integrate(fs0)            ! Integrate fission source
    Ke = Keo * f / fc             ! Update Keff
    CALL RelE(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
    CALL RelJ(joc, fer)           ! Search maximum point wise outgoing current Error
    IF ((ser < serc) .AND. (fer < ferc)) EXIT  ! If converge, exit.
END DO


END SUBROUTINE outer4th


SUBROUTINE outer4Fx

!
! Purpose:
!    To perform fixed-source outer iteration

USE sdata, ONLY: ng, nnod, nout, serc, ferc,  fer, ser, &
                 Ke, nac, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2
USE InpOutp, ONLY: ounit

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(nnod) :: fs0c                  !Old fission source
DOUBLE PRECISION, DIMENSION(nnod,ng,6) :: joc              !Old outgoing currents
DOUBLE PRECISION, DIMENSION(nnod) :: fsx1c, fsy1c, fsz1c
DOUBLE PRECISION, DIMENSION(nnod) :: fsx2c, fsy2c, fsz2c
DOUBLE PRECISION, DIMENSION(nnod) :: ss0                   ! Scattering source
DOUBLE PRECISION, DIMENSION(nnod) :: ssx1, ssy1, ssz1
DOUBLE PRECISION, DIMENSION(nnod) :: ssx2, ssy2, ssz2      ! Scattering source moments
DOUBLE PRECISION :: domiR, e1, e2
INTEGER :: g
INTEGER :: p, npos

DOUBLE PRECISION, DIMENSION(nnod) :: errn, erro

! Initialize fission source
fs0  = 0.0
fsx1 = 0.0; fsy1 = 0.0; fsz1 = 0.0
fsx2 = 0.0; fsy2 = 0.0; fsz2 = 0.0

DO g= 1, ng
    CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
END DO

errn = 1.0
e1 = Integrate(errn)

!Start outer iteration
DO p=1, nout
    CALL SaveJo(joc)  ! Save old outgoing currents
    fs0c  = fs0       ! Save old fission source
    fsx1c = fsx1; fsy1c = fsy1; fsz1c = fsz1
    fsx2c = fsx2; fsy2c = fsy2; fsz2c = fsz2
    fs0  = 0.0
    fsx1 = 0.0; fsy1 = 0.0; fsz1 = 0.0
    fsx2 = 0.0; fsy2 = 0.0; fsz2 = 0.0
    erro = errn       ! Save old fission source error/difference
    DO g = 1, ng

        !!!Calculate Scattering source
        CALL SSrc(g, ss0, ssx1, ssy1, ssz1, ssx2, ssy2, ssz2)

        !!!Calculate total source
        CALL TSrcFx(g, fs0c , fsx1c, fsy1c, fsz1c, &
                              fsx2c, fsy2c, fsz2c, &
                        ss0 , ssx1 , ssy1 , ssz1 , &
                              ssx2 , ssy2 , ssz2   )

        !!!Inner Iteration
        CALL inner4(g)

        !!!Calculate fission source for next outer iteration
        CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
    END DO

    errn = fs0 - fs0c
    e2 = Integrate(ABS(errn))

    IF (MOD(p,nac) == 0) THEN   ! Fission source extrapolation
        domiR = e2 / e1
        npos = MAXLOC(ABS(erro),1)
        IF (erro(npos) * errn(npos) < 0.0) domiR = -domiR
        fs0 = fs0 + domiR / (1.0 - domiR) * errn
    END IF
    e1 = e2                       ! Save integrated fission source error
    CALL RelE(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
    CALL RelJ(joc, fer)           ! Search maximum point wise outgoing current Error
    IF ((ser < serc) .AND. (fer < ferc)) EXIT  ! If converge, exit.
END DO


IF (p-1 == nout) THEN
    WRITE(ounit,*)
    WRITE(ounit,*) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED.'
    WRITE(ounit,*) '  ITERATION DO NOT CONVERGE.'
    WRITE(ounit,*) '  CHECK PROBLEM SPECIFICATION. ADPRES IS STOPING...'
    STOP
END IF

CALL MultF(Ke)

WRITE(ounit,*)
WRITE(ounit, '(A36,F9.6)') 'MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke

END SUBROUTINE outer4Fx


SUBROUTINE outertf (ht, maxi)

!
! Purpose:
!    To perform outer iteration for transient with flux transformation

USE sdata, ONLY: ng, nnod, serc, ferc, nout,  fer, ser, &
                 nac, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: ht
LOGICAL, INTENT(OUT) :: maxi

DOUBLE PRECISION, DIMENSION(nnod) :: fs0c                  !Old fission source
DOUBLE PRECISION, DIMENSION(nnod,ng,6) :: joc              !Old outgoing currents
DOUBLE PRECISION, DIMENSION(nnod) :: fsx1c, fsy1c, fsz1c
DOUBLE PRECISION, DIMENSION(nnod) :: fsx2c, fsy2c, fsz2c
DOUBLE PRECISION, DIMENSION(nnod) :: ss0                   ! Scattering source
DOUBLE PRECISION, DIMENSION(nnod) :: ssx1, ssy1, ssz1,ssx2, ssy2, ssz2      ! Scattering source moments
DOUBLE PRECISION :: domiR, e1, e2
INTEGER :: g
INTEGER :: p, npos

DOUBLE PRECISION, DIMENSION(nnod) :: errn, erro

errn = 1.0
e1 = Integrate(errn)

!Start outer iteration
DO p=1, nout
    CALL SaveJo(joc)  ! Save old outgoing currents
    fs0c  = fs0       ! Save old fission source
    fsx1c = fsx1; fsy1c = fsy1; fsz1c = fsz1
    fsx2c = fsx2; fsy2c = fsy2; fsz2c = fsz2
    fs0  = 0.0
    fsx1 = 0.0; fsy1 = 0.0; fsz1 = 0.0
    fsx2 = 0.0; fsy2 = 0.0; fsz2 = 0.0
    erro = errn       ! Save old fission source error/difference
    DO g = 1, ng

        !!!Calculate Scattering source
        CALL SSrc(g, ss0, ssx1, ssy1, ssz1, ssx2, ssy2, ssz2)

        !!!Calculate total source
        CALL TSrcT(g,  fs0c, fsx1c, fsy1c, fsz1c,&
                              fsx2c, fsy2c, fsz2c, &
                        ss0 , ssx1 , ssy1 , ssz1 , &
                              ssx2 , ssy2 , ssz2, ht)

        !!!Inner Iteration
        CALL inner4(g)

        !!!Calculate fission source for next outer iteration
        CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
    END DO

    errn = fs0 - fs0c
    e2 = Integrate(ABS(errn))

    IF (MOD(p,nac) == 0) THEN   ! Fission source extrapolation
        domiR = e2 / e1
        npos = MAXLOC(ABS(erro),1)
        IF (erro(npos) * errn(npos) < 0.0) domiR = -domiR
        fs0 = fs0 + domiR / (1.0 - domiR) * errn
    END IF
    e1 = e2                       ! Save integrated fission source error
    CALL RelE(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
    CALL RelJ(joc, fer)           ! Search maximum point wise outgoing current Error
    IF ((ser < serc) .AND. (fer < ferc)) EXIT  ! If converge, exit.
END DO

IF (p==nout+1) THEN
    maxi = .TRUE.
ELSE
    maxi = .FALSE.
END IF

END SUBROUTINE outertf


SUBROUTINE outer4ad(popt)

  !
  ! Purpose:
  !    To perform adjoint outer iteration

USE sdata, ONLY: ng, nnod, nout, serc, ferc,  fer, ser, &
                 Ke, nac, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2
USE InpOutp, ONLY: ounit

IMPLICIT NONE

INTEGER, OPTIONAL, INTENT(IN) :: popt

DOUBLE PRECISION :: Keo                                    !Old Multiplication factor (Keff)
DOUBLE PRECISION, DIMENSION(nnod) :: fs0c                  !Old fission source
DOUBLE PRECISION, DIMENSION(nnod,ng,6) :: joc              !Old outgoing currents
DOUBLE PRECISION, DIMENSION(nnod) :: fsx1c, fsy1c, fsz1c
DOUBLE PRECISION, DIMENSION(nnod) :: fsx2c, fsy2c, fsz2c
DOUBLE PRECISION, DIMENSION(nnod) :: ss0                   ! Scattering source
DOUBLE PRECISION, DIMENSION(nnod) :: ssx1, ssy1, ssz1
DOUBLE PRECISION, DIMENSION(nnod) :: ssx2, ssy2, ssz2      ! Scattering source moments
DOUBLE PRECISION :: f, fc                                  ! new and old integrated fission sources
DOUBLE PRECISION :: domiR, e1, e2
INTEGER :: g
INTEGER :: p, npos

LOGICAL :: opt

DOUBLE PRECISION, DIMENSION(nnod) :: errn, erro

IF (PRESENT(popt)) THEN
    opt = .FALSE.
ELSE
    opt = .TRUE.
END IF

! Initialize fission source
fs0  = 0.0
fsx1 = 0.0; fsy1 = 0.0; fsz1 = 0.0
fsx2 = 0.0; fsy2 = 0.0; fsz2 = 0.0

DO g= 1, ng
    CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
END DO

errn = 1.0
f = Integrate(fs0)
e1 = Integrate(errn)

!Start outer iteration
DO p=1, nout
    CALL SaveJo(joc)  ! Save old outgoing currents
    fc = f            ! Save old integrated fission source
    fs0c  = fs0       ! Save old fission source
    fsx1c = fsx1; fsy1c = fsy1; fsz1c = fsz1
    fsx2c = fsx2; fsy2c = fsy2; fsz2c = fsz2
    fs0  = 0.0
    fsx1 = 0.0; fsy1 = 0.0; fsz1 = 0.0
    fsx2 = 0.0; fsy2 = 0.0; fsz2 = 0.0
    Keo = Ke          ! Save old multiplication factor
    erro = errn       ! Save old fission source error/difference
    DO g = ng,1,-1

        !!!Calculate Scattering source
        CALL SSrcAd(g, ss0, ssx1, ssy1, ssz1, ssx2, ssy2, ssz2)

        !!!Calculate total source
        CALL TSrcAd(g, Keo, fs0c , fsx1c, fsy1c, fsz1c, &
                                 fsx2c, fsy2c, fsz2c, &
                           ss0 , ssx1 , ssy1 , ssz1 , &
                                 ssx2 , ssy2 , ssz2   )

        !!!Inner Iteration
        CALL inner4(g)

        !!!Calculate fission source for next outer iteration
        CALL FSrcAd (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
    END DO

    errn = fs0 - fs0c
    e2 = Integrate(ABS(errn))

    IF (MOD(p,nac) == 0) THEN   ! Fission source extrapolation
        domiR = e2 / e1
        npos = MAXLOC(ABS(erro),1)
        IF (erro(npos) * errn(npos) < 0.0) domiR = -domiR
        fs0 = fs0 + domiR / (1.0 - domiR) * errn
    END IF
    e1 = e2                       ! Save integrated fission source error
    f = Integrate(fs0)            ! Integrate fission source
    Ke = Keo * f / fc             ! Update Keff
    CALL RelE(fs0, fs0c, ser)     ! Search maximum point wise fission source Relative Error
    CALL RelJ(joc, fer)           ! Search maximum point wise outgoing current Error
    IF ((ser < serc) .AND. (fer < ferc)) EXIT  ! If converge, exit.
END DO

IF (p-1 == nout) THEN
    WRITE(ounit,*)
    WRITE(ounit,*) '  MAXIMUM NUMBER OF OUTER ITERATION IS REACHED.'
    WRITE(ounit,*) '  ITERATION MAY NOT CONVERGE.'
    WRITE(ounit,*) '  CHECK PROBLEM SPECIFICATION OR CHANGE ITERATION CONTROL (%ITER).'
    WRITE(ounit,*) '  ADPRES IS STOPING...'
    STOP
END IF

IF (opt) WRITE(ounit,*)
IF (opt) WRITE(ounit, '(A36,F9.6)') 'MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke

END SUBROUTINE outer4ad




SUBROUTINE inner4(g)
!
! Purpose:
!   To perform inner iterations
!

USE sdata, ONLY: nod, nnod, xstag, ystag, ierc, &
                xwest, xeast, ynorth, ysouth, zbott, ztop, &
                f0, ix, iy, iz, xyz, nzz, nin, al

IMPLICIT NONE

INTEGER, INTENT(IN) :: g
INTEGER :: l, n
DOUBLE PRECISION, DIMENSION(6) :: bvec, qvec

! Transverse Leakage Moments(0, Lx1, Ly1, Lz2, Lx2, Ly2, Lz3)
DOUBLE PRECISION, DIMENSION(nnod,7) :: Lm

! Jot Nodals' outgoing currents+flux  (X+, X-, Y+, Y-, Z+, Z-)
! Jin Nodals' ingoing currents+source (X+, X-, Y+, Y-, Z+, Z-)

DO l = 1, nin
    DO n = 1, nnod

            ! Calculate ingoing partial currents from neighborhod nodes
            IF (ix(n) == ystag(iy(n))%smax) THEN                          ! East (X+) BC
                CALL bcond(xeast, n, g, 1)
            ELSE
                nod(n,g)%ji(1) = (nod(xyz( ix(n)+1, iy(n), iz(n) ), g)%jo(2) + &
                                  al(n,g)%dc(1) * nod(n,g)%jo(1)) / &
                                  (1.0 - al(n,g)%dc(1))
            END IF

            IF (ix(n) == ystag(iy(n))%smin) THEN                          ! West (X-) BC
                CALL bcond(xwest, n, g, 2)
            ELSE
                nod(n,g)%ji(2) = (nod(xyz( ix(n)-1, iy(n), iz(n) ), g)%jo(1) + &
                                  al(n,g)%dc(2) * nod(n,g)%jo(2)) / &
                                  (1.0 - al(n,g)%dc(2))
            END IF

            IF (iy(n) == xstag(ix(n))%smax) THEN                          ! North (Y+) BC
                CALL bcond(ynorth, n, g, 3)
            ELSE
                nod(n,g)%ji(3) = (nod(xyz( ix(n), iy(n)+1, iz(n) ), g)%jo(4) + &
                                  al(n,g)%dc(3) * nod(n,g)%jo(3)) / &
                                  (1.0 - al(n,g)%dc(3))
            END IF

            IF (iy(n) == xstag(ix(n))%smin) THEN                          ! South (Y-) BC
                CALL bcond(ysouth, n, g, 4)
            ELSE
                nod(n,g)%ji(4) = (nod(xyz( ix(n), iy(n)-1, iz(n) ), g)%jo(3) + &
                                  al(n,g)%dc(4) * nod(n,g)%jo(4)) / &
                                  (1.0 - al(n,g)%dc(4))
            END IF

            IF (iz(n) == nzz) THEN                                    ! Top (Z+) BC
                CALL bcond(ztop, n, g, 5)
            ELSE
                nod(n,g)%ji(5) = (nod(xyz( ix(n), iy(n), iz(n)+1 ), g)%jo(6) + &
                                  al(n,g)%dc(5) * nod(n,g)%jo(5)) / &
                                  (1.0 - al(n,g)%dc(5))
            END IF

            IF (iz(n) == 1) THEN                                      ! Bottom (Z-)BC
                CALL bcond(zbott, n, g, 6)
            ELSE
                nod(n,g)%ji(6) = (nod(xyz( ix(n), iy(n), iz(n)-1 ), g)%jo(5) + &
                                  al(n,g)%dc(6) * nod(n,g)%jo(6)) / &
                                  (1.0 - al(n,g)%dc(6))
            END IF


            ! Update transverse leakage moments
            CALL TLUpd (n, g, Lm(n,:))
      END DO

      DO n = 1, nnod
            CALL matvec(nod(n,g)%P, nod(n,g)%ji, bvec)

            CALL matvec(nod(n,g)%R, nod(n,g)%Q - Lm(n,:), qvec)

            ! Update outgoing partial currents
            nod(n,g)%jo = qvec+bvec

            ! Update zeroth transverse leakages
            CALL LxyzUpd(n,g)

            ! Update flux and flux moments
            CALL FluxUpd4(n, g, Lm(n,:))
    END DO


    DO n = nnod, 1, -1

            ! Calculate ingoing partial currents from neighborhod nodes
            IF (ix(n) == ystag(iy(n))%smax) THEN                          ! East (X+) BC
                CALL bcond(xeast, n, g, 1)
            ELSE
                nod(n,g)%ji(1) = (nod(xyz( ix(n)+1, iy(n), iz(n) ), g)%jo(2) + &
                                  al(n,g)%dc(1) * nod(n,g)%jo(1)) / &
                                  (1.0 - al(n,g)%dc(1))
            END IF

            IF (ix(n) == ystag(iy(n))%smin) THEN                          ! West (X-) BC
                CALL bcond(xwest, n, g, 2)
            ELSE
                nod(n,g)%ji(2) = (nod(xyz( ix(n)-1, iy(n), iz(n) ), g)%jo(1) + &
                                  al(n,g)%dc(2) * nod(n,g)%jo(2)) / &
                                  (1.0 - al(n,g)%dc(2))
            END IF

            IF (iy(n) == xstag(ix(n))%smax) THEN                          ! North (Y+) BC
                CALL bcond(ynorth, n, g, 3)
            ELSE
                nod(n,g)%ji(3) = (nod(xyz( ix(n), iy(n)+1, iz(n) ), g)%jo(4) + &
                                  al(n,g)%dc(3) * nod(n,g)%jo(3)) / &
                                  (1.0 - al(n,g)%dc(3))
            END IF

            IF (iy(n) == xstag(ix(n))%smin) THEN                          ! South (Y-) BC
                CALL bcond(ysouth, n, g, 4)
            ELSE
                nod(n,g)%ji(4) = (nod(xyz( ix(n), iy(n)-1, iz(n) ), g)%jo(3) + &
                                  al(n,g)%dc(4) * nod(n,g)%jo(4)) / &
                                  (1.0 - al(n,g)%dc(4))
            END IF

            IF (iz(n) == nzz) THEN                                    ! Top (Z+) BC
                CALL bcond(ztop, n, g, 5)
            ELSE
                nod(n,g)%ji(5) = (nod(xyz( ix(n), iy(n), iz(n)+1 ), g)%jo(6) + &
                                  al(n,g)%dc(5) * nod(n,g)%jo(5)) / &
                                  (1.0 - al(n,g)%dc(5))
            END IF

            IF (iz(n) == 1) THEN                                      ! Bottom (Z-)BC
                CALL bcond(zbott, n, g, 6)
            ELSE
                nod(n,g)%ji(6) = (nod(xyz( ix(n), iy(n), iz(n)-1 ), g)%jo(5) + &
                                  al(n,g)%dc(6) * nod(n,g)%jo(6)) / &
                                  (1.0 - al(n,g)%dc(6))
            END IF


            ! Update transverse leakage moments
            CALL TLUpd (n, g, Lm(n,:))
      END DO

      DO n = nnod, 1, -1
            CALL matvec(nod(n,g)%P, nod(n,g)%ji, bvec)

            CALL matvec(nod(n,g)%R, nod(n,g)%Q - Lm(n,:), qvec)

            ! Update outgoing partial currents
            nod(n,g)%jo = qvec+bvec

            ! Update zeroth transverse leakages
            CALL LxyzUpd(n,g)

            ! Update flux and flux moments
            CALL FluxUpd4(n, g, Lm(n,:))
    END DO
END DO




END SUBROUTINE inner4



SUBROUTINE bcond (bc, nt, g, side)

!
! Purpose:
!    To provide proper boundary conditions
!

USE sdata, ONLY: nod

IMPLICIT NONE

INTEGER, INTENT(IN) :: bc, nt, g, side

IF (bc == 0) THEN
    nod(nt,g)%ji(side) = -nod(nt,g)%jo(side)
ELSE IF (bc == 1) THEN
    nod(nt,g)%ji(side) = 0.0
ELSE
    nod(nt,g)%ji(side) = nod(nt,g)%jo(side)
END IF

END SUBROUTINE bcond



SUBROUTINE FluxUpd4 (n, g, L)

USE sdata, ONLY: nod, D, sigr, xdel, ydel, zdel, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, &
                 ix, iy, iz

! Purpose:
   ! To update nod averaged flux and flux moments

IMPLICIT NONE

INTEGER, INTENT(IN) :: g, n
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: L

DOUBLE PRECISION :: Tx, Ty, Tz

! Calculate Zeroth Flux
f0(n,g)  = ( nod(n,g)%Q(1)          &
                 - nod(n,g)%L(1)/xdel(ix(n))   &
                 - nod(n,g)%L(2)/ydel(iy(n))   &
                 - nod(n,g)%L(3)/zdel(iz(n)))  &
                 / sigr(n,g)
IF (f0(n,g) < 0.) f0(n,g) = 0.

! Set parameters Tx, Ty and Tz
Tx = nod(n,g)%jo(1) - nod(n,g)%ji(1) &
   - nod(n,g)%jo(2) + nod(n,g)%ji(2)
Ty = nod(n,g)%jo(3) - nod(n,g)%ji(3) &
   - nod(n,g)%jo(4) + nod(n,g)%ji(4)
Tz = nod(n,g)%jo(5) - nod(n,g)%ji(5) &
   - nod(n,g)%jo(6) + nod(n,g)%ji(6)

! Calculate Flux moments
fx1(n,g) = ( nod(n,g)%Q(2) - L(2)                &
           - 0.5*Tx/xdel(ix(n))                  &
           - 2.*D(n,g)/xdel(ix(n))**2            &
           * (nod(n,g)%jo(1) + nod(n,g)%ji(1)    &
           -  nod(n,g)%jo(2) - nod(n,g)%ji(2)) ) &
           / sigr(n,g)


fy1(n,g) = ( nod(n,g)%Q(3) - L(3)                &
           - 0.5*Ty/ydel(iy(n))                  &
           - 2.*D(n,g)/ydel(iy(n))**2            &
           * (nod(n,g)%jo(3) + nod(n,g)%ji(3)    &
           -  nod(n,g)%jo(4) - nod(n,g)%ji(4)) ) &
           / sigr(n,g)


fz1(n,g) = ( nod(n,g)%Q(4) - L(4)                &
           - 0.5*Tz/zdel(iz(n))                  &
           - 2.*D(n,g)/zdel(iz(n))**2            &
           * (nod(n,g)%jo(5) + nod(n,g)%ji(5)    &
           -  nod(n,g)%jo(6) - nod(n,g)%ji(6)) ) &
           / sigr(n,g)


fx2(n,g) = ( nod(n,g)%Q(5) - L(5)                &
           - 0.5*nod(n,g)%L(1)/xdel(ix(n))       &
           - 6.*D(n,g)/xdel(ix(n))**2            &
           * (nod(n,g)%jo(1) + nod(n,g)%ji(1)    &
           +  nod(n,g)%jo(2) + nod(n,g)%ji(2)    &
           - f0(n,g)) ) / sigr(n,g)


fy2(n,g) = ( nod(n,g)%Q(6) - L(6)                &
           - 0.5*nod(n,g)%L(2)/ydel(iy(n))       &
           - 6.*D(n,g)/ydel(iy(n))**2            &
           * (nod(n,g)%jo(3) + nod(n,g)%ji(3)    &
           +  nod(n,g)%jo(4) + nod(n,g)%ji(4)    &
           - f0(n,g)) ) / sigr(n,g)


fz2(n,g) = ( nod(n,g)%Q(7) - L(7)                &
           - 0.5*nod(n,g)%L(3)/zdel(iz(n))       &
           - 6.*D(n,g)/zdel(iz(n))**2            &
           * (nod(n,g)%jo(5) + nod(n,g)%ji(5)    &
           +  nod(n,g)%jo(6) + nod(n,g)%ji(6)    &
           - f0(n,g)) ) / sigr(n,g)


END SUBROUTINE FluxUpd4


SUBROUTINE FluxUpd2 (nt, g)

USE sdata, ONLY: nod, sigr, xdel, ydel, zdel, &
                 f0, ix, iy, iz

! Purpose:
   ! To update nod averaged flux and flux moments

IMPLICIT NONE

INTEGER, INTENT(IN) :: g, nt

! Calculate Zeroth Flux
f0(nt,g)  = ( nod(nt,g)%Q(1)          &
                 - nod(nt,g)%L(1)/xdel(ix(nt))   &
                 - nod(nt,g)%L(2)/ydel(iy(nt))   &
                 - nod(nt,g)%L(3)/zdel(iz(nt)))  &
                 / sigr(nt,g)

END SUBROUTINE FluxUpd2


SUBROUTINE TLUpd (n, g, L)

USE sdata, ONLY: nod, xdel, ydel, zdel, xstag, ystag, nzz, &
                 xwest, xeast, ynorth, ysouth, zbott, ztop, &
                 ix, iy, iz, xyz

! Purpose:
   ! To calaculate transverse leakage moments

IMPLICIT NONE

INTEGER, INTENT(IN) :: g, n
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: L

DOUBLE PRECISION :: tm, tp
DOUBLE PRECISION :: p1m, p2m, p1p, p2p, pm, pp, p1, p2, p3
DOUBLE PRECISION :: r1xy, r2xy, r1xz, r2xz
DOUBLE PRECISION :: r1yx, r2yx, r1yz, r2yz
DOUBLE PRECISION :: r1zx, r2zx, r1zy, r2zy

! Set paramaters for X-Direction Transverse leakage
IF (ix(n) == ystag(iy(n))%smin) THEN
    IF (xwest == 0 .OR. xwest == 1) THEN
        tp = xdel(ix(n)+1)/xdel(ix(n))
        p1p = tp+1.0
        r1xy = 2. * ( nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(2)   &
             - nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)          &
             ) / p1p
        r2xy = 0.0
        r1xz = 2. * ( nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(3)   &
             - nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)          &
             ) / p1p
        r2xz = 0.0
    ELSE
        tm = 1.0
        tp = xdel(ix(n)+1)/xdel(ix(n))
        p1m = tm+1.0; p2m = 2.*tm+1.0; p1p = tp+1.0
        p2p = 2.*tp+1.0; p2 = tm+tp+2.
        pm = p1m*p2m; pp = p1p*p2p; p1 = pp-pm; p3 = p1m*p1p*(tm+tp+1.0)
        r1xy = (   pm * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(2)   &
                 - pp * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)   &
                 + p1 * nod(n,g)%L(2)                          &
               ) / p3
        r2xy = (   p1m * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(2)  &
                 + p1p * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)  &
                 - p2  * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)  &
               ) / p3
        r1xz = (   pm * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(3)   &
                 - pp * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)   &
                 + p1 * nod(n,g)%L(3)                          &
               ) / p3
        r2xz = (   p1m * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(3)  &
                 + p1p * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)  &
                 - p2  * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)  &
               ) / p3
    END IF
ELSE IF (ix(n) == ystag(iy(n))%smax) THEN
    IF (xeast == 0 .OR. xeast == 1) THEN
        tm = xdel(ix(n)-1)/xdel(ix(n))
        p1m = tm+1.0
        r1xy = 2. * ( nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)   &
             - nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(2)          &
             ) / p1m
        r2xy = 0.0
        r1xz = 2. * ( nod(xyz(ix(n)  ,iy(n) ,iz(n)),g)%L(3)  &
             - nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(3)          &
             ) / p1m
        r2xz = 0.0
    ELSE
        tm = xdel(ix(n)-1)/xdel(ix(n))
        tp = 1.0
        p1m = tm+1.0; p2m = 2.*tm+1.0;p1p = tp+1.0
        p2p = 2.*tp+1.0; p2 = tm+tp+2.
        pm = p1m*p2m; pp = p1p*p2p; p1 = pp-pm; p3 = p1m*p1p*(tm+tp+1.0)
        r1xy = (   pm * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)   &
                 - pp * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(2)   &
                 + p1 * nod(n,g)%L(2)                          &
               ) / p3
        r2xy = (   p1m * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)  &
                 + p1p * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(2)  &
                 - p2  * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)  &
               ) / p3
        r1xz = (   pm * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)   &
                 - pp * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(3)   &
                 + p1 * nod(n,g)%L(3)                          &
                 ) / p3
        r2xz = (   p1m * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)  &
                 + p1p * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(3)  &
                 - p2  * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)  &
               ) / p3
    END IF
ELSE
    tm = xdel(ix(n)-1)/xdel(ix(n))
    tp = xdel(ix(n)+1)/xdel(ix(n))
    p1m = tm+1.0; p2m = 2.*tm+1.0;p1p = tp+1.0
    p2p = 2.*tp+1.0; p2 = tm+tp+2.
    pm = p1m*p2m; pp = p1p*p2p; p1 = pp-pm; p3 = p1m*p1p*(tm+tp+1.0)
    r1xy = (   pm * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(2)   &
             - pp * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(2)   &
             + p1 * nod(n,g)%L(2)                          &
           ) / p3
    r2xy = (   p1m * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(2)  &
             + p1p * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(2)  &
             - p2  * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)  &
           ) / p3
    r1xz = (   pm * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(3)   &
             - pp * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(3)   &
             + p1 * nod(n,g)%L(3)                          &
           ) / p3
    r2xz = (   p1m * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(3)  &
             + p1p * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(3)  &
             - p2  * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)  &
           ) / p3
END IF


! Set paramaters for Y-Direction Transverse leakage
IF (iy(n) == xstag(ix(n))%smin) THEN
    IF (ysouth == 0 .OR. ysouth == 1) THEN
        tp = ydel(iy(n)+1)/ydel(iy(n))
        p1p = tp+1.0
        r1yx = 2. * ( nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(1)   &
             - nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)          &
             ) / p1p
        r2yx = 0.0
        r1yz = 2. * ( nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(3)   &
             - nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)          &
             ) / p1p
        r2yz = 0.0
    ELSE
        tm = 1.0
        tp = ydel(iy(n)+1)/ydel(iy(n))
        p1m = tm+1.0; p2m = 2.*tm+1.0;p1p = tp+1.0
        p2p = 2.*tp+1.0; p2 = tm+tp+2.
        pm = p1m*p2m; pp = p1p*p2p; p1 = pp-pm; p3 = p1m*p1p*(tm+tp+1.0)
        r1yx = (   pm * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(1)   &
                 - pp * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)   &
                 + p1 * nod(n,g)%L(1)                          &
               ) / p3
        r2yx = (   p1m * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(1)  &
                 + p1p * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)  &
                 - p2  * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)  &
               ) / p3
        r1yz = (   pm * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(3)   &
                 - pp * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)   &
                 + p1 * nod(n,g)%L(3)                          &
               ) / p3
        r2yz = (   p1m * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(3)  &
                 + p1p * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)  &
                 - p2  * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)  &
               ) / p3
    END IF
ELSE IF (iy(n) == xstag(ix(n))%smax) THEN
    IF (ynorth == 0 .OR. ynorth == 1) THEN
        tm = ydel(iy(n)-1)/ydel(iy(n))
        p1m = tm+1.0
        r1yx = 2. * ( nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)   &
             - nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(1)          &
             ) / p1m
        r2yx = 0.0
        r1yz = 2. * ( nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)   &
             - nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(3)          &
             ) / p1m
        r2yz = 0.0
    ELSE
        tm = ydel(iy(n)-1)/ydel(iy(n))
        tp = 1.0
        p1m = tm+1.0; p2m = 2.*tm+1.0;p1p = tp+1.0
        p2p = 2.*tp+1.0; p2 = tm+tp+2.
        pm = p1m*p2m; pp = p1p*p2p; p1 = pp-pm; p3 = p1m*p1p*(tm+tp+1.0)
        r1yx = (   pm * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)   &
                 - pp * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(1)   &
                 + p1 * nod(n,g)%L(1)                          &
               ) / p3
        r2yx = (   p1m * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)  &
                 + p1p * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(1)  &
                 - p2  * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)  &
                 ) / p3
        r1yz = (   pm * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)   &
                 - pp * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(3)   &
                 + p1 * nod(n,g)%L(3)                          &
               ) / p3
        r2yz = (   p1m * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)  &
                 + p1p * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(3)  &
                 - p2  * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)  &
               ) / p3
    END IF
ELSE
    tm = ydel(iy(n)-1)/ydel(iy(n))
    tp = ydel(iy(n)+1)/ydel(iy(n))
    p1m = tm+1.0; p2m = 2.*tm+1.0;p1p = tp+1.0
    p2p = 2.*tp+1.0; p2 = tm+tp+2.
    pm = p1m*p2m; pp = p1p*p2p; p1 = pp-pm; p3 = p1m*p1p*(tm+tp+1.0)
    r1yx = (   pm * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(1)   &
             - pp * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(1)   &
             + p1 * nod(n,g)%L(1)                          &
           ) / p3
    r2yx = (   p1m * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(1)  &
             + p1p * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(1)  &
             - p2  * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)  &
           ) / p3
    r1yz = (   pm * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(3)   &
             - pp * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(3)   &
             + p1 * nod(n,g)%L(3)                          &
           ) / p3
    r2yz = (   p1m * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(3)  &
             + p1p * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(3)  &
             - p2  * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)  &
           ) / p3
END IF

! Set paramaters for Z-Direction Transverse leakage
IF (iz(n) == 1 ) THEN
    IF (zbott == 0 .OR. zbott == 1) THEN
        tp = zdel(iz(n)+1)/zdel(iz(n))
        p1p = tp+1.0
        r1zx = 2. * ( nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(1)   &
             - nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)          &
             ) / p1p
        r2zx = 0.0
        r1zy = 2. * ( nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(2)   &
             - nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)          &
             ) / p1p
        r2zy = 0.0
    ELSE
        tm = 1.0
        tp = zdel(iz(n)+1)/zdel(iz(n))
        p1m = tm+1.0; p2m = 2.*tm+1.0; p1p = tp+1.0
        p2p = 2.*tp+1.0; p2 = tm+tp+2.
        pm = p1m*p2m; pp = p1p*p2p; p1 = pp-pm; p3 = p1m*p1p*(tm+tp+1.0)
        r1zx = (   pm * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(1)   &
                 - pp * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)   &
                 + p1 * nod(n,g)%L(1)                          &
               ) / p3
        r2zx = (   p1m * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(1)  &
                 + p1p * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)  &
                 - p2  * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)  &
               ) / p3
        r1zy = (   pm * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(2)   &
                 - pp * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)   &
                 + p1 * nod(n,g)%L(2)                          &
               ) / p3
        r2zy = (   p1m * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(2)  &
                 + p1p * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)  &
                 - p2  * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)  &
               ) / p3
    END IF
ELSE IF (iz(n) == nzz) THEN
    IF (ztop == 0 .OR. ztop == 1) THEN
        tm = zdel(iz(n)-1)/zdel(iz(n))
        p1m = tm+1.0
        r1zx = 2. * ( nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)   &
             - nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(1)          &
             ) / p1m
        r2zx = 0.0
        r1zy = 2. * ( nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)   &
             - nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(2)          &
             ) / p1m
        r2zy = 0.0
    ELSE
        tm = zdel(iz(n)-1)/zdel(iz(n))
        tp = 1.0
        p1m = tm+1.0; p2m = 2.*tm+1.0; p1p = tp+1.0
        p2p = 2.*tp+1.0; p2 = tm+tp+2.
        pm = p1m*p2m; pp = p1p*p2p; p1 = pp-pm; p3 = p1m*p1p*(tm+tp+1.0)
        r1zx = (   pm * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)   &
                 - pp * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(1)   &
                 + p1 * nod(n,g)%L(1)                          &
               ) / p3
        r2zx = (    p1m * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)  &
                  + p1p * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(1)  &
                  - p2  * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)  &
               ) / p3
        r1zy = (   pm * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)   &
                 - pp * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(2)   &
                 + p1 * nod(n,g)%L(2)                          &
               ) / p3
        r2zy = (   p1m * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)  &
                 + p1p * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(2)  &
                 - p2  * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)  &
               ) / p3
    END IF
ELSE
    tm = zdel(iz(n)-1)/zdel(iz(n))
    tp = zdel(iz(n)+1)/zdel(iz(n))
    p1m = tm+1.0; p2m = 2.*tm+1.0; p1p = tp+1.0
    p2p = 2.*tp+1.0; p2 = tm+tp+2.
    pm = p1m*p2m; pp = p1p*p2p; p1 = pp-pm; p3 = p1m*p1p*(tm+tp+1.0)
    r1zx = (   pm * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(1)   &
             - pp * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(1)   &
             + p1 * nod(n,g)%L(1)                          &
           ) / p3
    r2zx = (   p1m * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(1)  &
             + p1p * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(1)  &
             - p2  * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)  &
           ) / p3
    r1zy = (   pm * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(2)   &
             - pp * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(2)   &
             + p1 * nod(n,g)%L(2)                          &
           ) / p3
    r2zy = (   p1m * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(2)  &
             + p1p * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(2)  &
             - p2  * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)  &
           ) / p3
END IF

! Set Transverse leakage Moments
L(1) = 0.0
L(2) = ( r1xy/ydel(iy(n))+r1xz/zdel(iz(n)) ) / 12.
L(3) = ( r1yx/xdel(ix(n))+r1yz/zdel(iz(n)) ) / 12.
L(4) = ( r1zx/xdel(ix(n))+r1zy/ydel(iy(n)) ) / 12.
L(5) = ( r2xy/ydel(iy(n))+r2xz/zdel(iz(n)) ) / 20.
L(6) = ( r2yx/xdel(ix(n))+r2yz/zdel(iz(n)) ) / 20.
L(7) = ( r2zx/xdel(ix(n))+r2zy/ydel(iy(n)) ) / 20.

END SUBROUTINE TLUpd


SUBROUTINE FSrc(g, s, sx1, sy1, sz1, sx2, sy2, sz2)
!
! Purpose:
!   To calculate fission source and fission source moments
!

USE sdata, ONLY: nnod, nuf, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2

IMPLICIT NONE

INTEGER, INTENT(IN) :: g
DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: s, sx1, sy1, sz1
DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: sx2, sy2, sz2

INTEGER :: n

DO n = 1, nnod
    s(n)   = s(n)   + f0 (n,g) * nuf(n,g)
    sx1(n) = sx1(n) + fx1(n,g) * nuf(n,g)
    sy1(n) = sy1(n) + fy1(n,g) * nuf(n,g)
    sz1(n) = sz1(n) + fz1(n,g) * nuf(n,g)
    sx2(n) = sx2(n) + fx2(n,g) * nuf(n,g)
    sy2(n) = sy2(n) + fy2(n,g) * nuf(n,g)
    sz2(n) = sz2(n) + fz2(n,g) * nuf(n,g)
END DO

END SUBROUTINE FSrc


SUBROUTINE FSrcAd(g, s, sx1, sy1, sz1, sx2, sy2, sz2)
!
! Purpose:
!   To calculate fission source and fission source moments for adjoint calc.
!

USE sdata, ONLY: nnod, chi, mat, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2

IMPLICIT NONE

INTEGER, INTENT(IN) :: g
DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: s, sx1, sy1, sz1
DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: sx2, sy2, sz2

INTEGER :: n

DO n = 1, nnod
    s(n)   = s(n)   + f0 (n,g) * chi(mat(n),g)
    sx1(n) = sx1(n) + fx1(n,g) * chi(mat(n),g)
    sy1(n) = sy1(n) + fy1(n,g) * chi(mat(n),g)
    sz1(n) = sz1(n) + fz1(n,g) * chi(mat(n),g)
    sx2(n) = sx2(n) + fx2(n,g) * chi(mat(n),g)
    sy2(n) = sy2(n) + fy2(n,g) * chi(mat(n),g)
    sz2(n) = sz2(n) + fz2(n,g) * chi(mat(n),g)
END DO

END SUBROUTINE FSrcAd



SUBROUTINE SSrc(g, s, sx1, sy1, sz1, sx2, sy2, sz2)
!
! Purpose:
!   To calculate scattering source and scattering source moments
!

USE sdata, ONLY: ng, nnod, sigs, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2

IMPLICIT NONE

INTEGER, INTENT(IN) :: g
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: s, sx1, sy1, sz1
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: sx2, sy2, sz2

INTEGER :: h, n

s = 0.; sx1 = 0.; sy1 = 0.; sz1 = 0.
sx2 = 0.; sy2 = 0.; sz2 = 0.

DO h = 1, ng
    DO n = 1, nnod
        IF (g /= h) THEN
            s(n)   = s(n)   + sigs(n,h,g) * f0(n,h)
            sx1(n) = sx1(n) + sigs(n,h,g) * fx1(n,h)
            sy1(n) = sy1(n) + sigs(n,h,g) * fy1(n,h)
            sz1(n) = sz1(n) + sigs(n,h,g) * fz1(n,h)
            sx2(n) = sx2(n) + sigs(n,h,g) * fx2(n,h)
            sy2(n) = sy2(n) + sigs(n,h,g) * fy2(n,h)
            sz2(n) = sz2(n) + sigs(n,h,g) * fz2(n,h)
        END IF
    END DO
END DO

END SUBROUTINE SSrc



SUBROUTINE SSrcAd(g, s, sx1, sy1, sz1, sx2, sy2, sz2)
!
! Purpose:
!   To calculate scattering source and scattering source moments for adjoint calc.
!

USE sdata, ONLY: ng, nnod, sigs, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2

IMPLICIT NONE

INTEGER, INTENT(IN) :: g
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: s, sx1, sy1, sz1
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: sx2, sy2, sz2

INTEGER :: h, n

s = 0.; sx1 = 0.; sy1 = 0.; sz1 = 0.
sx2 = 0.; sy2 = 0.; sz2 = 0.

DO h = 1, ng
    DO n = 1, nnod
        IF (g /= h) THEN
            s(n)   = s(n)   + sigs(n,g,h) * f0(n,h)
            sx1(n) = sx1(n) + sigs(n,g,h) * fx1(n,h)
            sy1(n) = sy1(n) + sigs(n,g,h) * fy1(n,h)
            sz1(n) = sz1(n) + sigs(n,g,h) * fz1(n,h)
            sx2(n) = sx2(n) + sigs(n,g,h) * fx2(n,h)
            sy2(n) = sy2(n) + sigs(n,g,h) * fy2(n,h)
            sz2(n) = sz2(n) + sigs(n,g,h) * fz2(n,h)
        END IF
    END DO
END DO

END SUBROUTINE SSrcAd



SUBROUTINE TSrc(g, Keff, sf0, sfx1, sfy1, sfz1, sfx2, sfy2, sfz2, &
                          s0,  sx1,  sy1,  sz1,  sx2,  sy2 , sz2   )
!
! Purpose:
!   To update total source
!

USE sdata, ONLY: nod, chi, mat, nnod

IMPLICIT NONE

INTEGER, INTENT(IN) :: g
DOUBLE PRECISION, INTENT(IN)    :: Keff
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: sf0, sfx1, sfy1, sfz1
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: sfx2, sfy2, sfz2
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: s0, sx1, sy1, sz1
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: sx2, sy2, sz2

INTEGER :: n

DO n = 1, nnod
    nod(n,g)%Q(1) = chi(mat(n),g) * sf0(n)/Keff  + s0(n)
    nod(n,g)%Q(2) = chi(mat(n),g) * sfx1(n)/Keff + sx1(n)
    nod(n,g)%Q(3) = chi(mat(n),g) * sfy1(n)/Keff + sy1(n)
    nod(n,g)%Q(4) = chi(mat(n),g) * sfz1(n)/Keff + sz1(n)
    nod(n,g)%Q(5) = chi(mat(n),g) * sfx2(n)/Keff + sx2(n)
    nod(n,g)%Q(6) = chi(mat(n),g) * sfy2(n)/Keff + sy2(n)
    nod(n,g)%Q(7) = chi(mat(n),g) * sfz2(n)/Keff + sz2(n)
END DO

END SUBROUTINE TSrc


SUBROUTINE TSrcFx(g, sf0, sfx1, sfy1, sfz1, sfx2, sfy2, sfz2, &
                      s0,  sx1,  sy1,  sz1,  sx2,  sy2 , sz2   )
!
! Purpose:
!   To update total source for fixed source calcs.
!

USE sdata, ONLY: nod, chi, mat, nnod, exsrc

IMPLICIT NONE

INTEGER, INTENT(IN) :: g
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: sf0, sfx1, sfy1, sfz1
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: sfx2, sfy2, sfz2
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: s0, sx1, sy1, sz1
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: sx2, sy2, sz2

INTEGER :: n

DO n = 1, nnod
    nod(n,g)%Q(1) = chi(mat(n),g) * sf0(n)  + s0(n)  + exsrc(n,g)
    nod(n,g)%Q(2) = chi(mat(n),g) * sfx1(n) + sx1(n)
    nod(n,g)%Q(3) = chi(mat(n),g) * sfy1(n) + sy1(n)
    nod(n,g)%Q(4) = chi(mat(n),g) * sfz1(n) + sz1(n)
    nod(n,g)%Q(5) = chi(mat(n),g) * sfx2(n) + sx2(n)
    nod(n,g)%Q(6) = chi(mat(n),g) * sfy2(n) + sy2(n)
    nod(n,g)%Q(7) = chi(mat(n),g) * sfz2(n) + sz2(n)
END DO

END SUBROUTINE TSrcFx


SUBROUTINE TSrcT(g, sf0, sfx1, sfy1, sfz1, sfx2, sfy2, sfz2, &
                      s0,  sx1,  sy1,  sz1,  sx2,  sy2 , sz2, h)
!
! Purpose:
!   To update total source for transient calcs. with exponetial transformation
!

USE sdata, ONLY: nod, chi, mat, nnod, tbeta, velo, lamb, iBeta, nf, omeg, &
                 c0, cx1, cy1, cz1, cx2, cy2, cz2, &
                 ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2

IMPLICIT NONE

INTEGER, INTENT(IN) :: g
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: sf0, sfx1, sfy1, sfz1, sfx2, sfy2, sfz2
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: s0, sx1, sy1, sz1, sx2, sy2, sz2
DOUBLE PRECISION, INTENT(IN) :: h

DOUBLE PRECISION :: dt, dtx1, dty1, dtz1, dtx2, dty2, dtz2, lat, dfis
INTEGER :: n, i

DO n = 1, nnod
     dt = 0.; dtx1 = 0.; dty1 = 0.; dtz1 = 0.; dtx2 = 0.; dty2 = 0.; dtz2 = 0.
     dfis = 0.
     DO i = 1, nf
        lat = 1. + lamb(i) * h
        dt = dt  + lamb(i) * c0(n,i) / lat
        dtx1 = dtx1 + lamb(i) * cx1(n,i) / lat
        dty1 = dty1 + lamb(i) * cy1(n,i) / lat
        dtz1 = dtz1 + lamb(i) * cz1(n,i) / lat
        dtx2 = dtx2 + lamb(i) * cx2(n,i) / lat
        dty2 = dty2 + lamb(i) * cy2(n,i) / lat
        dtz2 = dtz2 + lamb(i) * cz2(n,i) / lat
        dfis = dfis + chi(mat(n),g) * iBeta(i) * lamb(i) * h / lat
    END DO

    nod(n,g)%Q(1) = ((1. - tbeta) * chi(mat(n),g) + dfis) * sf0(n)  &
    + s0(n) + chi(mat(n),g) * dt + ft(n,g)  * EXP(omeg(n,g) * h) / (velo(g) * h)
    nod(n,g)%Q(2) = ((1. - tbeta) * chi(mat(n),g) + dfis) * sfx1(n)  &
    + sx1(n) + chi(mat(n),g) * dtx1 + ftx1(n,g)  * EXP(omeg(n,g) * h) / (velo(g) * h)
    nod(n,g)%Q(3) = ((1. - tbeta) * chi(mat(n),g) + dfis) * sfy1(n)  &
    + sy1(n) + chi(mat(n),g) * dty1 + fty1(n,g)  * EXP(omeg(n,g) * h) / (velo(g) * h)
    nod(n,g)%Q(4) = ((1. - tbeta) * chi(mat(n),g) + dfis) * sfz1(n)  &
    + sz1(n) + chi(mat(n),g) * dtz1 + ftz1(n,g)  * EXP(omeg(n,g) * h) / (velo(g) * h)
    nod(n,g)%Q(5) = ((1. - tbeta) * chi(mat(n),g) + dfis) * sfx2(n)  &
    + sx2(n) + chi(mat(n),g) * dtx2 + ftx2(n,g)  * EXP(omeg(n,g) * h) / (velo(g) * h)
    nod(n,g)%Q(6) = ((1. - tbeta) * chi(mat(n),g) + dfis) * sfy2(n)  &
    + sy2(n) + chi(mat(n),g) * dty2 + fty2(n,g)  * EXP(omeg(n,g) * h) / (velo(g) * h)
    nod(n,g)%Q(7) = ((1. - tbeta) * chi(mat(n),g) + dfis) * sfz2(n)  &
    + sz2(n) + chi(mat(n),g) * dtz2 + ftz2(n,g)  * EXP(omeg(n,g) * h) / (velo(g) * h)
END DO

END SUBROUTINE TSrcT


SUBROUTINE TSrcAd(g, Keff, sf0, sfx1, sfy1, sfz1, sfx2, sfy2, sfz2, &
                          s0,  sx1,  sy1,  sz1,  sx2,  sy2 , sz2   )
!
! Purpose:
!   To update total source for adjoint calc.
!

USE sdata, ONLY: nod, nuf, nnod

IMPLICIT NONE

INTEGER, INTENT(IN) :: g
DOUBLE PRECISION, INTENT(IN)    :: Keff
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: sf0, sfx1, sfy1, sfz1
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: sfx2, sfy2, sfz2
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: s0, sx1, sy1, sz1
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: sx2, sy2, sz2

INTEGER :: n

DO n = 1, nnod
    nod(n,g)%Q(1) = nuf(n,g) * sf0(n)/Keff  + s0(n)
    nod(n,g)%Q(2) = nuf(n,g) * sfx1(n)/Keff + sx1(n)
    nod(n,g)%Q(3) = nuf(n,g) * sfy1(n)/Keff + sy1(n)
    nod(n,g)%Q(4) = nuf(n,g) * sfz1(n)/Keff + sz1(n)
    nod(n,g)%Q(5) = nuf(n,g) * sfx2(n)/Keff + sx2(n)
    nod(n,g)%Q(6) = nuf(n,g) * sfy2(n)/Keff + sy2(n)
    nod(n,g)%Q(7) = nuf(n,g) * sfz2(n)/Keff + sz2(n)
END DO

END SUBROUTINE TSrcAd


SUBROUTINE LxyzUpd (nt,g)

USE sdata, ONLY: nod

! Purpose:
   ! To update Transverse leakages for group g and nod n

IMPLICIT NONE

INTEGER, INTENT(IN) :: g, nt

nod(nt,g)%L(1) = nod(nt,g)%jo(1) - nod(nt,g)%ji(1) &
                - nod(nt,g)%ji(2) + nod(nt,g)%jo(2)
nod(nt,g)%L(2) = nod(nt,g)%jo(3) - nod(nt,g)%ji(3) &
                - nod(nt,g)%ji(4) + nod(nt,g)%jo(4)
nod(nt,g)%L(3) = nod(nt,g)%jo(5) - nod(nt,g)%ji(5) &
                - nod(nt,g)%ji(6) + nod(nt,g)%jo(6)

END SUBROUTINE LxyzUpd


SUBROUTINE nodal_coup4()
!
! Purpose:
!    To calculate nodal coupling matrix
!

USE sdata, ONLY: ng, nnod, xdel, ydel, zdel, &
                 ix, iy, iz, D, sigr, nod

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(6,6) :: A, B
DOUBLE PRECISION, DIMENSION(6,7) :: C

INTEGER :: n, g

DOUBLE PRECISION:: dx, dy, dz, lx, ly, lz
DOUBLE PRECISION :: ax, ay, az, ax1, ay1, az1, bx, by, bz
DOUBLE PRECISION :: bx1, by1, bz1, xy, xz, yx, yz, zx, zy

DO g= 1, ng
    DO n = 1, nnod

        dx = D(n,g) / xdel(ix(n))
        dy = D(n,g) / ydel(iy(n))
        dz = D(n,g) / zdel(iz(n))

        lx = 1.0 / sigr(n,g) / xdel(ix(n))
        ly = 1.0 / sigr(n,g) / ydel(iy(n))
        lz = 1.0 / sigr(n,g) / zdel(iz(n))

        ax = 1.0+32.0*dx+120.0*dx*lx+960.0*dx*dx*lx+840.0*dx*dx*lx*lx
        ay = 1.0+32.0*dy+120.0*dy*ly+960.0*dy*dy*ly+840.0*dy*dy*ly*ly
        az = 1.0+32.0*dz+120.0*dz*lz+960.0*dz*dz*lz+840.0*dz*dz*lz*lz

        ax1 = 8.0*dx+60.0*dx*lx+720.0*dx*dx*lx+840.0*dx*dx*lx*lx
        ay1 = 8.0*dy+60.0*dy*ly+720.0*dy*dy*ly+840.0*dy*dy*ly*ly
        az1 = 8.0*dz+60.0*dz*lz+720.0*dz*dz*lz+840.0*dz*dz*lz*lz

        bx = 1.0-32.0*dx+120.0*dx*lx-960.0*dx*dx*lx+840.0*dx*dx*lx*lx
        by = 1.0-32.0*dy+120.0*dy*ly-960.0*dy*dy*ly+840.0*dy*dy*ly*ly
        bz = 1.0-32.0*dz+120.0*dz*lz-960.0*dz*dz*lz+840.0*dz*dz*lz*lz

        bx1 = -8.0*dx+60.0*dx*lx-720.0*dx*dx*lx+840.0*dx*dx*lx*lx
        by1 = -8.0*dy+60.0*dy*ly-720.0*dy*dy*ly+840.0*dy*dy*ly*ly
        bz1 = -8.0*dz+60.0*dz*lz-720.0*dz*dz*lz+840.0*dz*dz*lz*lz

        xy = 20.*dx*ly+840.0*dx*dx*lx*ly
        xz = 20.*dx*lz+840.0*dx*dx*lx*lz

        yx = 20.*dy*lx+840.0*dy*dy*ly*lx
        yz = 20.*dy*lz+840.0*dy*dy*ly*lz

        zx = 20.*dz*lx+840.0*dz*dz*lz*lx
        zy = 20.*dz*ly+840.0*dz*dz*lz*ly

        A(1,1) = ax; A(2,2) = ax
        A(3,3) = ay; A(4,4) = ay
        A(5,5) = az; A(6,6) = az

        A(1,2) = ax1; A(2,1) = ax1
        A(3,4) = ay1; A(4,3) = ay1
        A(5,6) = az1; A(6,5) = az1

        A(1,3) = xy; A(1,4) = xy
        A(2,3) = xy; A(2,4) = xy
        A(1,5) = xz; A(1,6) = xz
        A(2,5) = xz; A(2,6) = xz

        A(3,1) = yx; A(3,2) = yx
        A(4,1) = yx; A(4,2) = yx
        A(3,5) = yz; A(3,6) = yz
        A(4,5) = yz; A(4,6) = yz

        A(5,1) = zx; A(5,2) = zx
        A(6,1) = zx; A(6,2) = zx
        A(5,3) = zy; A(5,4) = zy
        A(6,3) = zy; A(6,4) = zy

        B = A

        B(1,1) = bx; B(2,2) = bx
        B(3,3) = by; B(4,4) = by
        B(5,5) = bz; B(6,6) = bz

        B(1,2) = bx1; B(2,1) = bx1
        B(3,4) = by1; B(4,3) = by1
        B(5,6) = bz1; B(6,5) = bz1

        C = 0.0

        ax = 20.*dx*lx*xdel(ix(n))+840.0*dx*dx*lx*lx*xdel(ix(n))
        ay = 20.*dy*ly*ydel(iy(n))+840.0*dy*dy*ly*ly*ydel(iy(n))
        az = 20.*dz*lz*zdel(iz(n))+840.0*dz*dz*lz*lz*zdel(iz(n))

        C(1,1) = ax
        C(2,1) = ax
        C(3,1) = ay
        C(4,1) = ay
        C(5,1) = az
        C(6,1) = az

        ax1 = 60.0*dx*lx*xdel(ix(n))
        ay1 = 60.0*dy*ly*ydel(iy(n))
        az1 = 60.0*dz*lz*zdel(iz(n))

        C(1,2) =  ax1
        C(2,2) = -ax1
        C(3,3) =  ay1
        C(4,3) = -ay1
        C(5,4) =  az1
        C(6,4) = -az1

        ax1 = 140.0*dx*lx*xdel(ix(n))
        ay1 = 140.0*dy*ly*ydel(iy(n))
        az1 = 140.0*dz*lz*zdel(iz(n))

        C(1,5) = ax1
        C(2,5) = ax1
        C(3,6) = ay1
        C(4,6) = ay1
        C(5,7) = az1
        C(6,7) = az1

        CALL inverse(g, n, A)

        nod(n,g)%P = MATMUL(A, B)
        nod(n,g)%R = MATMUL(A, C)

    END DO
END DO

END SUBROUTINE nodal_coup4


SUBROUTINE inverse (g, nt, mat)

!
! Purpose:
!    To perform matrix inverse by LU decomposition
!

USE InpOutp, ONLY: ounit
USE sdata,   ONLY: ix, iy, iz

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: mat
INTEGER, INTENT(IN) :: g, nt

DOUBLE PRECISION, DIMENSION(6,6) :: L, U, imat, pmat
DOUBLE PRECISION, DIMENSION(6) :: y
DOUBLE PRECISION :: piv, isum
INTEGER :: i, j, k

pmat = mat
U = mat
L = 0.0

! Start matrix decomposition
DO i= 1, 6
    IF (ABS(mat(i,i)) < 10e-3) THEN
      WRITE(ounit,*) 'ERROR IN MATRIX DECOMP: DIAGONAL ELEMENTS CLOSE TO ZERO'
      WRITE(ounit,2001) g, ix(nt), iy(nt), iz(nt)
      STOP
    END IF
    L(i,i) = 1.0
    DO j= i+1, 6
        piv = U(j,i)/U(i,i)
        L(j,i) = piv
        DO k= i, 6
            U(j,k) = U(j,k) - piv*U(i,k)
        END DO
        U(j,i) = 0.0
    END DO
END DO

! Check matrix decomposition
DO i = 1,6
    DO j = 1,6
        isum = 0.0
        DO k = 1,6
            isum = isum+L(i,k)*U(k,j)
        END DO
        IF (ABS(mat(i,j)-isum)/ABS(mat(i,j)) > 1.d-3) THEN
            WRITE(ounit,*) 'ERROR IN MATRIX DECOMP: DECOMPOSITION FAILED'
            WRITE(ounit,2001) g, ix(nt), iy(nt), iz(nt)
            STOP
        END IF
    END DO
END DO

!Initialiaze Identity matrix
imat = 0.0
DO i= 1, 6
    imat(i,i) = 1.0
END DO

! Calculate matrix inverse
! Ref: https://www.gamedev.net/resources/_/technical/math-and-physics/matrix-inversion-using-lu-decomposition-r3637
DO j=1,6   ! For each column
    !Solve y in Ly = b (Forward substitution)
    y(1) = imat(1,j)
    DO i=2,6
        isum = 0.0
        DO k =1, i-1
            isum = isum + L(i,k)*y(k)
        END DO
        y(i) = imat(i,j)-isum
    END DO

    ! Solve x in Ux=y(Backward substitution) and store inverse matrix to input matrix 'mat'
    mat(6,j) = y(6)/U(6,6)
    DO i = 5,1,-1
        isum = 0.0
        DO k =i+1,6
            isum = isum + U(i,k)*mat(k,j)
        END DO
        mat(i,j) = (y(i)-isum) / U(i,i)
    END DO
END DO

!Check matrix Inverse
DO i = 1,6
    DO j = 1,6
        isum = 0.0
        DO k = 1,6
            isum = isum+pmat(i,k)*mat(k,j)
        END DO
        IF (ABS(imat(i,j)-isum) > 1.d-4) THEN
            WRITE(ounit,*) 'ERROR IN MATRIX INVERSION'
            WRITE(ounit,2001) g, ix(nt), iy(nt), iz(nt)
            STOP
        END IF
    END DO
END DO

2001 FORMAT(2X, 'Group = ', I2, ', I = ', I2, ', J = ', I2, ', K = ', I2)


END SUBROUTINE inverse


DOUBLE PRECISION FUNCTION Integrate(s)

  !
  ! Purpose:
  !    To perform volume integration

USE sdata, ONLY: nnod, vdel

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION (:), INTENT(IN) :: s
INTEGER :: n

Integrate = 0.
DO n = 1, nnod
    Integrate = Integrate + vdel(n) * s(n)
END DO

END FUNCTION Integrate



SUBROUTINE matvec (mat, vec, rvec)
!
! Purpose:
!    To perform matrix vector multiplication
!

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: mat
DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: vec
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: rvec

INTEGER :: i, j, n, m
DOUBLE PRECISION :: isum

m = SIZE(mat,1)
n = SIZE(vec,1)

DO i= 1, m
    isum = 0.
    DO j = 1, n
        isum = isum + mat(i,j)*vec(j)
    END DO
    rvec(i)     = isum
END DO

END SUBROUTINE matvec



SUBROUTINE SaveJo(jox)

  !
  ! Purpose:
  !    To calculate Max Relative error for outgoing currents

USE sdata, ONLY: nnod, ng, nod

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(OUT) :: jox

INTEGER :: n, g, i

DO i = 1, 6
   DO g = 1, ng
      DO n= 1, nnod
          jox(n,g,i) = nod(n,g)%jo(i)
      END DO
   END DO
END DO

END SUBROUTINE SaveJo



SUBROUTINE RelJ(jox,rel)

  !
  ! Purpose:
  !    To calculate Max Relative error for outgoing currents

USE sdata, ONLY: nnod, ng, nod

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN) :: jox
DOUBLE PRECISION, INTENT(OUT) :: rel

DOUBLE PRECISION :: error
INTEGER :: n, g, i

rel = 0.

DO i = 1, 6
   DO g = 1, ng
      DO n= 1, nnod
         IF (ABS(nod(n,g)%jo(i)) > 1.d-10) THEN
             error = ABS(nod(n,g)%jo(i) - jox(n,g,i)) / ABS(nod(n,g)%jo(i))
             IF (error > rel) rel = error
         END IF
      END DO
  END DO
END DO

END SUBROUTINE RelJ



SUBROUTINE RelE(newF, oldF, rel)

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
        error = ABS(newF(n) - oldF(n)) / ABS(newF(n))
        IF (error > rel) rel = error
    END IF
END DO

END SUBROUTINE RelE


SUBROUTINE MultF(k)

! Purpose: To calculate Keff for fixed source problem

USE sdata, ONLY: ng, nnod, nod, nzz, &
                 ix, iy, iz, nuf, siga, f0, &
                 xstag, ystag, xdel, ydel, zdel

IMPLICIT NONE

DOUBLE PRECISION, INTENT(OUT) :: k

DOUBLE PRECISION :: leak, absp, fiss

INTEGER :: g, n

!! Jot Nodals' outgoing currents  (X+, X-, Y+, Y-, Z+, Z-)
!! Jin Nodals' ingoing currents   (X+, X-, Y+, Y-, Z+, Z-)

leak = 0.0
absp = 0.0
fiss = 0.0
DO g = 1, ng
    DO n = 1, nnod
        !! Get leakages
        IF (ix(n) == ystag(iy(n))%smax) leak = leak + (nod(n,g)%Jo(1) - nod(n,g)%Ji(1)) &
                                             * ydel(iy(n)) * zdel(iz(n))
        IF (ix(n) == ystag(iy(n))%smin) leak = leak + (nod(n,g)%Jo(2) - nod(n,g)%Ji(2)) &
                                             * ydel(iy(n)) * zdel(iz(n))
        IF (iy(n) == xstag(ix(n))%smax) leak = leak + (nod(n,g)%Jo(3) - nod(n,g)%Ji(3)) &
                                             * xdel(ix(n)) * zdel(iz(n))
        IF (iy(n) == xstag(ix(n))%smin) leak = leak + (nod(n,g)%Jo(4) - nod(n,g)%Ji(4)) &
                                             * xdel(ix(n)) * zdel(iz(n))
        IF (iz(n) == nzz) leak = leak + (nod(n,g)%Jo(5) - nod(n,g)%Ji(5)) &
                               * xdel(ix(n)) * ydel(iy(n))
        IF (iz(n) == 1)   leak = leak + (nod(n,g)%Jo(6) - nod(n,g)%Ji(6)) &
                               * xdel(ix(n)) * ydel(iy(n))

        absp = absp + siga(n,g) * f0(n,g) * xdel(ix(n)) * ydel(iy(n)) * zdel(iz(n))
        fiss = fiss + nuf(n,g) * f0(n,g) * xdel(ix(n)) * ydel(iy(n)) * zdel(iz(n))
    END DO
END DO

! WRITE(*,*) fiss, leak , absp
k = fiss / (leak + absp)

END SUBROUTINE MultF



SUBROUTINE Init()

!
! Purpose:
!    To provide initial guess for some variables
!



USE sdata, ONLY: ng, nod, nnod, Ke, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, &
                 fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2
USE InpOutp, ONLY: ounit, brrst, runit, ounit

IMPLICIT NONE

INTEGER :: rnnod, rng
INTEGER :: g, n, istat

ALLOCATE(nod(nnod,ng))

ALLOCATE(f0(nnod,ng), &
         fx1(nnod,ng), fy1(nnod,ng), fz1(nnod,ng), &
         fx2(nnod,ng), fy2(nnod,ng), fz2(nnod,ng), &
         STAT=istat)
IF (istat /= 0) THEN
    WRITE(ounit,*) '[1] NOT ENOUGH MEMORY. PROGRAM STOP'
    STOP
END IF

IF (brrst == 1) THEN
    READ(runit,*) rng, rnnod
    IF (rng /= ng) THEN
        WRITE(ounit,*) '  ERROR: NUMBER OF GROUP IN RESTART DOES NOT MATCH'
        STOP
    END IF
    IF (rnnod /= nnod) THEN
        WRITE(ounit,*) '  ERROR: GEOMETRY (NUMBER OF NODES) IN RESTART DOES NOT MATCH'
        STOP
    END IF

    READ(runit,*) Ke

    DO g= 1, ng
        DO n = 1, nnod
            READ(runit,*) f0(n,g), fx1(n,g), fy1(n,g), &
                          fz1(n,g), fx2(n,g), fy2(n,g), fz2(n,g)
        END DO
    END DO

    DO g= 1, ng
        DO n = 1, nnod
            READ(runit,*) nod(n,g)%ji(1), nod(n,g)%ji(2), nod(n,g)%ji(3), &
                          nod(n,g)%ji(4), nod(n,g)%ji(5), nod(n,g)%ji(6)
            READ(runit,*) nod(n,g)%jo(1), nod(n,g)%jo(2), nod(n,g)%jo(3), &
                          nod(n,g)%jo(4), nod(n,g)%jo(5), nod(n,g)%jo(6)
            CALL LxyzUpd(n,g)
        END DO
    END DO

ELSE
    Ke = 1.0
    DO g= 1, ng
        DO n = 1, nnod
            nod(n,g)%jo = 1.0
            nod(n,g)%ji = 1.0

            CALL LxyzUpd(n,g)

            f0(n,g)  = 1.0
            fx1(n,g) = 1.0
            fy1(n,g) = 1.0
            fz1(n,g) = 1.0
            fx2(n,g) = 1.0
            fy2(n,g) = 1.0
            fz2(n,g) = 1.0
        END DO
    END DO
END IF

! Allocate fission source
ALLOCATE (fs0(nnod), fsx1(nnod), fsy1(nnod), fsz1(nnod))
ALLOCATE (fsx2(nnod), fsy2(nnod), fsz2(nnod))

END SUBROUTINE Init



SUBROUTINE PowDis (p)

!
! Purpose:
!    To calculate power distribution
!


USE sdata, ONLY: ng, nnod, sigf, f0, vdel, mode
USE InpOutp, ONLY: ounit

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: p
INTEGER :: g, n
DOUBLE PRECISION :: tpow, pow

p = 0.0
DO g= 1, ng
    DO n= 1, nnod
       pow = f0(n,g) * sigf(n,g) * vdel(n)
    IF (pow < 0.) pow = 0.
        p(n) = p(n) + pow
    END DO
END DO

! Normalize to 1.0
tpow = 0.
DO n = 1, nnod
    tpow = tpow + p(n)
END DO

IF (tpow <= 0 .AND. mode /= 'FIXEDSRC') THEN
   WRITE(ounit, *) '   ERROR: TOTAL NODES POWER IS ZERO OR LESS'
   WRITE(ounit, *) '   STOP IN SUBROUTINE POWDIS'
   STOP
END IF

DO n = 1, nnod
    p(n) = p(n) / tpow
END DO


END SUBROUTINE PowDis


SUBROUTINE PowTot (fx,tpow)

!
! Purpose:
!    To calculate power distribution
!


USE sdata, ONLY: ng, nnod, sigf, vdel

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: fx
DOUBLE PRECISION, INTENT(OUT) :: tpow

DOUBLE PRECISION, DIMENSION(nnod) :: p
INTEGER :: g, n

p = 0.0
DO g= 1, ng
    DO n= 1, nnod
        p(n) = p(n) + fx(n,g) * sigf(n,g) * vdel(n)
    END DO
END DO


tpow = 0.
DO n = 1, nnod
    tpow = tpow + p(n)
END DO

END SUBROUTINE PowTot


FUNCTION AbsAr (ar)

!
! Purpose:
!    To calculate power distribution
!


USE sdata, ONLY: nnod

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: ar
DOUBLE PRECISION, DIMENSION(nnod) :: AbsAr

INTEGER :: n

DO n= 1, nnod
   AbsAr(n) = ABS(ar(n))
END DO

END FUNCTION AbsAr


SUBROUTINE forward()

!
! Purpose:
!    To solve forward (normal) problems
!

USE sdata, ONLY: nnod, f0, aprad, apaxi, afrad, ftem, mtem, cden, bcon, bpos
USE InpOutp, ONLY: ounit, AsmPow, AxiPow, AsmFlux, XS_updt

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: pow


WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) ' ==============================================' &
// '=================================================='
WRITE(ounit,*) &
'                                       CALCULATION RESULTS'
WRITE(ounit,*) ' ==============================================' &
// '=================================================='
WRITE(ounit,*)
WRITE(ounit,*) '  Itr     k-eff     Fis.Src Error   Inner Error'
WRITE(ounit,*) ' ----------------------------------------------------'

! Update XSEC
CALL XS_updt(bcon, ftem, mtem, cden, bpos)

! Calculate Nodal coupling matrices
CALL nodal_coup4()

CALL outer4()


IF (aprad == 1 .OR. apaxi == 1) THEN
    ALLOCATE(pow(nnod))
    CALL PowDis(pow)
END IF

IF (aprad == 1) CALL AsmPow(pow)

IF (apaxi == 1) CALL AxiPow(pow)

IF (afrad == 1) CALL AsmFlux(f0, 1.d0)



END SUBROUTINE forward



SUBROUTINE adjoint()

  !
  ! Purpose:
  !    To solve adjoint problems
  !

USE sdata, ONLY: nnod, f0, aprad, apaxi, afrad, ftem, mtem, cden, bcon, bpos
USE InpOutp, ONLY: ounit, AsmPow, AxiPow, AsmFlux, XS_updt

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: pow


WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) ' ==============================================' &
// '=================================================='
WRITE(ounit,*) &
'                                       CALCULATION RESULTS'
WRITE(ounit,*) ' ==============================================' &
// '=================================================='
WRITE(ounit,*)
WRITE(ounit,*) '  Itr     k-eff     Fis.Src Error   Inner Error'
WRITE(ounit,*) ' ----------------------------------------------------'

! Update XSEC
CALL XS_updt(bcon, ftem, mtem, cden, bpos)

! Calculate Nodal coupling matrices
CALL nodal_coup4()

CALL outer4ad()

IF (aprad == 1 .OR. apaxi == 1) THEN
    ALLOCATE(pow(nnod))
    CALL PowDis(pow)
END IF

IF (aprad == 1) CALL AsmPow(pow)

IF (apaxi == 1) CALL AxiPow(pow)

IF (afrad == 1) CALL AsmFlux(f0, 1.d0)



END SUBROUTINE adjoint



SUBROUTINE fixedsrc()

  !
  ! Purpose:
  !    To solve fixed source problems
  !

USE sdata, ONLY: nnod, f0, aprad, apaxi, afrad, ftem, mtem, cden, bcon, bpos
USE InpOutp, ONLY: ounit, AsmPow, AxiPow, AsmFlux, XS_updt

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: pow


WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) ' ==============================================' &
// '=================================================='
WRITE(ounit,*) &
'                                       CALCULATION RESULTS'
WRITE(ounit,*) ' ==============================================' &
// '=================================================='
WRITE(ounit,*)
WRITE(ounit,*) '  Itr   Fis.Src Error   Inner Error'
WRITE(ounit,*) ' ----------------------------------------------------'

! Update XSEC
CALL XS_updt(bcon, ftem, mtem, cden, bpos)

! Calculate Nodal coupling matrices
CALL nodal_coup4()

CALL outer4Fx()

IF (aprad == 1 .OR. apaxi == 1) THEN
    ALLOCATE(pow(nnod))
    CALL PowDis(pow)
END IF

IF (aprad == 1) CALL AsmPow(pow)

IF (apaxi == 1) CALL AxiPow(pow)

IF (afrad == 1) CALL AsmFlux(f0)



END SUBROUTINE fixedsrc


END MODULE nodal
