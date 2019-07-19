MODULE InpOutp2

!=========================
! Input output module to read, process and echo input, as well as writing output
! =======================

USE sdata, ONLY: DP

IMPLICIT NONE

SAVE

CHARACTER(LEN=1) :: ind          ! used to read x indicator in input buffer to prevent error
CHARACTER(LEN=100):: iline  ! Input line
CHARACTER(LEN=100) :: message    ! error message


CONTAINS

SUBROUTINE inp_xtab(xbunit)

!
! Purpose:
!    To read tabular xsec file
!

USE sdata, ONLY: ng, nmat, nf, xmat, ncden, nbcon, nftem, nmtem

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln !Line number
INTEGER :: iost  ! IOSTAT status

INTEGER :: i, j, g, s,t,u,v

! XTAB TYPE
TYPE :: XFILE
    CHARACTER(LEN=100) :: fname             ! XTAB File name
    INTEGER :: cnum                        ! Composition number in the XTAB Files
END TYPE
TYPE(XFILE), DIMENSION(:), ALLOCATABLE :: xtab
! INTEGER, DIMENSION(:), ALLOCATABLE :: group

LOGICAL, DIMENSION(:), ALLOCATABLE :: noty   !to check if this buffer was read or not?
INTEGER, PARAMETER :: xunit = 998  !XTAB file unit number
INTEGER, PARAMETER :: tunit = 999  !XTAB Buffer unit number
INTEGER :: comm  ! Position of comment mark

INTEGER :: nskip  !Number of lines to skip

REAL(DP), DIMENSION(:), ALLOCATABLE :: sigtr, siga, nuf, sigf     !Remove this
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sigs, dc                         !Remove this

WRITE(101,*)
WRITE(101,*)
WRITE(101,*) '           >>>> READING TABULAR XSEC FILE <<<<'
WRITE(101,*) '           --------------------------------------------'

READ(xbunit, *, IOSTAT=iost) ind, ln, ng, nmat  !Read numbef of group and material
message = ' error in material number'
CALL er_message(101, iost, ln, message)

ALLOCATE(xtab(nmat), noty(nmat), xmat(nmat))
noty = .TRUE.

! Reading XTAB file names and composition number in input file
DO i= 1, nmat
    READ(xbunit, '(A2,I5,A100)', IOSTAT=iost) ind, ln, iline
    iline = ADJUSTL(iline)         !Adjust to left
    comm = INDEX(iline, ' ')       ! Get space position
    xtab(i)%fname = iline(1:comm-1)  !Get xtab file name
    READ(iline(comm:100),*) xtab(i)%cnum  ! Get composition number (convert to integer)
    message = ' error in reading XTAB files'
    CALL er_message(101, iost, ln, message)
END DO

! Starting to read XTAB File, remove comments and write to buffer file
DO i = 1, nmat
  IF (noty(i)) THEN  !If this composition was not written in buffer
    ! Open XTAB File
    CALL openFile(xunit, xtab(i)%fname, 'XTAB File Open Failed--status')

    ! Start removing comments and rewrite into one input XTAB buffer
    CALL inp_comments(xunit, tunit, '*')

    !This loop to read another composition in the same XTAB file
    DO j = i, nmat
      ! If next material has the same file name
      IF (TRIM(ADJUSTL(xtab(i)%fname)) == TRIM(ADJUSTL(xtab(j)%fname))) THEN

        !Read buffer file and saved the xsec and transient data
        REWIND(tunit)
        READ(tunit, *,IOSTAT=iost) ind, ln, xmat(j)%tadf, xmat(j)%trod     ! Read input control
        message = ' ERROR IN XTAB FILE: CANNOT READ CONTROL PARAMETERS'
        CALL er_message(101, iost, ln, message, j)
        READ(tunit, *, IOSTAT=iost) ind, ln, ncden, nbcon, nftem, nmtem    ! Read number of branch
        message = ' ERROR IN XTAB FILE: CANNOT READ BRANCH DIMENSION'
        CALL er_message(101, iost, ln, message, j)

        ! Check branch dimension
        IF (ncden < 1 .OR. nbcon < 1 .OR. nftem < 1 .OR. nmtem < 1) THEN
          WRITE(101, *) ' ERROR IN XTAB FILE ', TRIM(ADJUSTL(xtab(j)%fname))
          WRITE(101, *) ' ERROR: MINIMUM NUMBER OF BRANCH IS 1'
          WRITE(*, *) ' ERROR IN XTAB FILE ', TRIM(ADJUSTL(xtab(j)%fname))
          WRITE(*, *) ' ERROR: MINIMUM NUMBER OF BRANCH IS 1'
          STOP
        END IF

        ! ALLOCATE XSEC DATA
        ALLOCATE(xmat(j)%chi(ng), xmat(j)%velo(ng))
        ALLOCATE(xmat(j)%xsec(ncden, nbcon, nftem, nmtem))
        ALLOCATE(xmat(j)%rxsec(ncden, nbcon, nftem, nmtem))
        DO s = 1, ncden
          DO t = 1, nbcon
            DO u = 1, nftem
              DO v = 1, nmtem
                ALLOCATE(xmat(j)%xsec(s,t,u,v)%sigtr(ng))
                ALLOCATE(xmat(j)%xsec(s,t,u,v)%siga(ng))
                ALLOCATE(xmat(j)%xsec(s,t,u,v)%sigf(ng))
                ALLOCATE(xmat(j)%xsec(s,t,u,v)%nuf(ng))
                ALLOCATE(xmat(j)%xsec(s,t,u,v)%dc(6,ng))
                ALLOCATE(xmat(j)%xsec(s,t,u,v)%sigs(ng,ng))
                IF (xmat(j)%trod == 1) THEN
                  ALLOCATE(xmat(j)%rxsec(s,t,u,v)%sigtr(ng))
                  ALLOCATE(xmat(j)%rxsec(s,t,u,v)%siga(ng))
                  ALLOCATE(xmat(j)%rxsec(s,t,u,v)%sigf(ng))
                  ALLOCATE(xmat(j)%rxsec(s,t,u,v)%nuf(ng))
                  ALLOCATE(xmat(j)%rxsec(s,t,u,v)%dc(6,ng))
                  ALLOCATE(xmat(j)%rxsec(s,t,u,v)%sigs(ng,ng))
                END IF
              END DO
            END DO
          END DO
        END DO

        ! Allocate and read branch paramaters
        CALL branchPar(tunit, ncden, j, xmat(j)%pcden, xtab(j)%fname, &
        'COOLANT DENSITY')  ! Allocate and read coolant dens. branc paramaters
        CALL branchPar(tunit, nbcon, j, xmat(j)%pbcon, xtab(j)%fname, &
        'BORON CONCENTRATION')  ! Allocate and read Boron conc. branc paramaters
        CALL branchPar(tunit, nftem, j, xmat(j)%pftem, xtab(j)%fname, &
        'FUEL TEMPERATURE')  ! Allocate and read fule temp. branc paramaters
        CALL branchPar(tunit, nmtem, j, xmat(j)%pmtem, xtab(j)%fname, &
        'MODERATOR TEMPERATURE')  ! Allocate and read moderator temp. branc paramaters

        ! Skip lines to read desired composition in the xtab file
        nskip = ng*nbcon*nftem*nmtem
        IF (xmat(j)%tadf  == 1) THEN  ! IF dc present
          IF (xmat(j)%trod  == 1) THEN
            CALL skipRead(tunit, j, (xtab(j)%cnum-1)*(10*nskip+2*ng*nskip+4))
          ELSE
            CALL skipRead(tunit, j, (xtab(j)%cnum-1)*(5*nskip+ng*nskip+4))
          END IF
        ELSE IF (xmat(j)%tadf  == 2) THEN
          IF (xmat(j)%trod  == 1) THEN
            CALL skipRead(tunit, j, (xtab(j)%cnum-1)*(20*nskip+2*ng*nskip+4))
          ELSE
            CALL skipRead(tunit, j, (xtab(j)%cnum-1)*(10*nskip+ng*nskip+4))
          END IF
        ELSE
          IF (xmat(j)%trod  == 1) THEN
            CALL skipRead(tunit, j, (xtab(j)%cnum-1)*(8*nskip+2*ng*nskip+4))
          ELSE
            CALL skipRead(tunit, j, (xtab(j)%cnum-1)*(4*nskip+ng*nskip+4))
          END IF
        END IF

        ! Read unrodded XSEC
        CALL readXS (tunit, j, 0, xmat(j)%xsec)
        ! Read rodded XSEC
        CALL readXS (tunit, j, 1, xmat(j)%rxsec)
        !Read fission spectrum
        READ(tunit, *, IOSTAT=iost) ind, ln, (xmat(j)%chi(g), g = 1, ng)
        message = ' ERROR IN XTAB FILE: CANNOT READ FISSION SPECTRUM'
        CALL er_message(101, iost, ln, message, j)
        !Read neutron Inverse velocity
        READ(tunit, *, IOSTAT=iost) ind, ln, (xmat(j)%velo(g), g = 1, ng)
        message = ' ERROR IN XTAB FILE: CANNOT READ Inverse Velocity'
        CALL er_message(101, iost, ln, message, j)
        DO g = 1, ng
          xmat(j)%velo(g) = 1._DP/xmat(j)%velo(g)  !COnvert to velocity
        END DO
        ! Read decay constant
        READ(tunit, *, IOSTAT=iost) ind, ln, (xmat(j)%lamb(t), t = 1, nf)
        message = ' ERROR IN XTAB FILE: CANNOT READ DECAY CONSTANT'
        CALL er_message(101, iost, ln, message, j)
        ! Read beta
        READ(tunit, *, IOSTAT=iost) ind, ln, (xmat(j)%iBeta(t), t = 1, nf)
        message = ' ERROR IN XTAB FILE: CANNOT READ DELAYED NEUTRON FRACTION'
        CALL er_message(101, iost, ln, message, j)

        ! If read, indicate that it has been read
        noty(j) = .FALSE.
      END IF
    END DO

    ! Close XTAB File and buffer file
    CLOSE(UNIT=xunit); CLOSE(UNIT=tunit)
  END IF
END DO

DEALLOCATE(xtab, noty)


! xmat(1)%xsec = 2.*xmat(1)%xsec
! DO u = 1, nftem
!   DO t = 1, nbcon
!     WRITE(*,'(3ES16.5)') (2*xmat(1)%rxsec(s,t,u,1)%dc(1), s = 1, ncden)
!   END DO
! END DO

ALLOCATE(sigtr(ng), siga(ng), nuf(ng), sigf(ng), sigs(ng,ng), dc(6,ng))
CALL brInterp(0,1, 0.71187_DP,  500._DP, 560._DP, 560._DP, sigtr, siga, nuf, sigf, sigs, dc)

WRITE(*, '(10ES16.5)') sigtr

CALL brInterp(1,1, 0.730_DP,  0._DP, 560._DP, 560._DP, sigtr, siga, nuf, sigf, sigs, dc)

WRITE(*, '(10ES16.5)') sigs

WRITE(*,*) 'OK'
STOP

END SUBROUTINE inp_xtab

!******************************************************************************!

SUBROUTINE readXS (tunit, matnum, rod, xsec)
!Purpose: To read xsec in XTAB file

USE sdata, ONLY: ng, XBRANCH, xmat, ncden, nbcon, nftem, nmtem

IMPLICIT NONE

INTEGER, INTENT(IN) :: tunit, matnum, rod  ! file unit number, material number, and rod indicator
TYPE(XBRANCH), DIMENSION(:,:,:,:), INTENT(INOUT) :: xsec  !Set INOUT, see: http://www.cs.rpi.edu/~szymansk/OOF90/bugs.html#2
INTEGER :: iost, ln
INTEGER :: g, h, s, t, u, v, k

!Read sigtr
DO g = 1, ng
  DO v = 1, nmtem
    DO u = 1, nftem
      DO t = 1, nbcon
        READ(tunit, *, IOSTAT=iost) ind, ln, &
        (xsec(s,t,u,v)%sigtr(g), s = 1, ncden)
        IF (rod == 0) THEN
          message = ' ERROR IN XTAB FILE: CANNOT READ TRANSPORT XSEC'
        ELSE
          message = ' ERROR IN XTAB FILE: CANNOT READ RODDED TRANSPORT XSEC'
        END IF
        CALL er_message(101, iost, ln, message, matnum)
      END DO
    END DO
  END DO
END DO
!Read siga
DO g = 1, ng
  DO v = 1, nmtem
    DO u = 1, nftem
      DO t = 1, nbcon
        READ(tunit, *, IOSTAT=iost) ind, ln, &
        (xsec(s,t,u,v)%siga(g), s = 1, ncden)
        IF (rod == 0) THEN
          message = ' ERROR IN XTAB FILE: CANNOT READ ABSORPTION XSEC'
        ELSE
          message = ' ERROR IN XTAB FILE: CANNOT READ RODDED ABSORPTION XSEC'
        END IF
        CALL er_message(101, iost, ln, message, matnum)
      END DO
    END DO
  END DO
END DO
!Read nu*sigf
DO g = 1, ng
  DO v = 1, nmtem
    DO u = 1, nftem
      DO t = 1, nbcon
        READ(tunit, *, IOSTAT=iost) ind, ln, &
        (xsec(s,t,u,v)%nuf(g), s = 1, ncden)
        IF (rod == 0) THEN
          message = ' ERROR IN XTAB FILE: CANNOT READ NU*SIGF XSEC'
        ELSE
          message = ' ERROR IN XTAB FILE: CANNOT READ RODDED NU*SIGF XSEC'
        END IF
        CALL er_message(101, iost, ln, message, matnum)
      END DO
    END DO
  END DO
END DO
!Read kappa*sigf
DO g = 1, ng
  DO v = 1, nmtem
    DO u = 1, nftem
      DO t = 1, nbcon
        READ(tunit, *, IOSTAT=iost) ind, ln, &
        (xsec(s,t,u,v)%sigf(g), s = 1, ncden)
        IF (rod == 0) THEN
          message = ' ERROR IN XTAB FILE: CANNOT READ KAPPA*SIGF XSEC'
        ELSE
          message = ' ERROR IN XTAB FILE: CANNOT READ RODDED KAPPA*SIGF XSEC'
        END IF
        CALL er_message(101, iost, ln, message, matnum)
      END DO
    END DO
  END DO
END DO
!Read sigs
DO g = 1, ng
  DO h = 1, ng
    DO v = 1, nmtem
      DO u = 1, nftem
        DO t = 1, nbcon
          READ(tunit, *, IOSTAT=iost) ind, ln, &
          (xsec(s,t,u,v)%sigs(g,h), s = 1, ncden)
          IF (rod == 0) THEN
            message = ' ERROR IN XTAB FILE: CANNOT READ SCATTERING XSEC'
          ELSE
            message = ' ERROR IN XTAB FILE: CANNOT READ RODDED SCATTERING XSEC'
          END IF
          CALL er_message(101, iost, ln, message, matnum)
        END DO
      END DO
    END DO
  END DO
END DO
!Read dc
IF (xmat(matnum)%tadf  == 1) THEN  ! IF dc present
  DO g = 1, ng
    DO v = 1, nmtem
      DO u = 1, nftem
        DO t = 1, nbcon
          READ(tunit, *, IOSTAT=iost) ind, ln, &
          (xsec(s,t,u,v)%dc(1,g), s = 1, ncden)
          DO k = 1, 6
            DO s = 1, ncden
              xsec(s,t,u,v)%dc(k,g) = xsec(s,t,u,v)%dc(1,g)
            END DO
          END DO
          IF (rod == 0) THEN
            message = ' ERROR IN XTAB FILE: CANNOT READ ADFs'
          ELSE
            message = ' ERROR IN XTAB FILE: CANNOT READ RODDED ADFs'
          END IF
          CALL er_message(101, iost, ln, message, matnum)
        END DO
      END DO
    END DO
  END DO
ELSE IF (xmat(matnum)%tadf  == 2) THEN
  DO g = 1, ng
    DO k = 1, 6
      DO v = 1, nmtem
        DO u = 1, nftem
          DO t = 1, nbcon
            READ(tunit, *, IOSTAT=iost) ind, ln, &
            (xsec(s,t,u,v)%dc(k,g), s = 1, ncden)
            IF (rod == 0) THEN
              message = ' ERROR IN XTAB FILE: CANNOT READ ADFs'
            ELSE
              message = ' ERROR IN XTAB FILE: CANNOT READ RODDED TRANSPORT ADFs'
            END IF
            CALL er_message(101, iost, ln, message, matnum)
          END DO
        END DO
      END DO
    END DO
  END DO
ELSE
  CONTINUE
END IF


END SUBROUTINE readXS

!******************************************************************************!

SUBROUTINE skipRead (iunit,matnum, nskip)
!Purpose: To allocate and read branch paramaters

IMPLICIT NONE

INTEGER, INTENT(IN) :: iunit, matnum, nskip
INTEGER :: i, eof

DO i = 1, nskip
  READ (iunit, *, IOSTAT=eof)
  IF (eof < 0) THEN              !Check end of file
    WRITE(101,1131) matnum
    WRITE(101,1132)
    WRITE(*,1131) matnum
    WRITE(*,1132)
    STOP
  END IF
END DO

1131 FORMAT(2X, 'ERROR: END OF FILE REACHED FOR XTAB FILE IN MATERIAL NUMBER ', I3)
1132 FORMAT(2X, 'ADPRES IS STOPPING')

END SUBROUTINE skipRead

!******************************************************************************!

SUBROUTINE branchPar (tunit, dim, matnum, par, fname, messPar)
!Purpose: To allocate and read branch paramaters

IMPLICIT NONE

INTEGER, INTENT(IN) :: tunit, dim, matnum
CHARACTER(LEN=*), INTENT(IN) :: messPar, fname
REAL(DP), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: par

INTEGER :: k
INTEGER :: ln, iost

IF (dim > 1) THEN                         ! If branch DIMENSION > 1
  ALLOCATE(par(dim))
  READ(tunit, *, IOSTAT=iost) ind, ln, par(1:dim)
  message = ' ERROR IN XTAB FILE: CANNOT READ BRANCH PARAMETERS ' // messPar
  CALL er_message(101, iost, ln, message, matnum)
  DO k = 2, dim
    IF (par(k-1) > par(k)) THEN
      WRITE(101,*) "  ERROR IN XTAB FILE  ", fname
      WRITE(101,*) "  ", messPar, " PARAMETER SHALL BE IN ORDER, SMALL to BIG"
      STOP
    END IF
  END DO
ELSE
  ALLOCATE(par(1))
  par(1) = 1.0    !Arbitrary
END IF

END SUBROUTINE branchPar

!******************************************************************************!

! SUBROUTINE xxs_updt (xbcon, xftem, xmtem, xcden, xbpos)
! !
! ! Purpose: Update xsec for given paramters
! !
!
! USE sdata, ONLY: ng, nxx, nyy, nzz, xyz, zdel, mat, nod, cusp, f0, &
!                  sigtr, siga, nuf, sigf, sigs, &
!                  dsigtr, dsiga, dnuf, dsigf, dsigs, &
!                  nod, f0, fz1, fz2, coreh, fbmap, pos0, ssize
!
! IMPLICIT NONE
!
! REAL(DP), DIMENSION(:), INTENT(IN) :: xftem  ! Provided fuel temperature
! REAL(DP), DIMENSION(:), INTENT(IN) :: xmtem  ! Provided moderator temperature
! REAL(DP), DIMENSION(:), INTENT(IN) :: xcden  ! Provided coolant density
! REAL(DP), DIMENSION(:), INTENT(IN) :: xbpos  ! Provided control rod bank position
!
! INTEGER ::i, j, k, g, h
! REAL(DP) :: rodh, vfrac
! REAL(DP) :: dum
!
! INTEGER :: n, n1, n2, nmax
! REAL(DP) :: del1, del2, eta1, eta2
! REAL(DP) :: sum1, sum2, sum3, sum4, sumx
! REAL(DP), DIMENSION(ng) :: sum5
! REAL(DP), DIMENSION(:), ALLOCATABLE :: f
! REAL(DP) :: a1, a2, a3, a4, x, tx, f1, f2
!
! DO j = 1, nyy
!   DO i = 1, nxx
!      IF (fbmap(i,j) > 0) THEN
!         !!!(rodh -> posistion the tip of the control rod the top of core)
!          rodh = coreh - pos0  - xbpos(fbmap(i,j))*ssize
!          dum = 0._DP
!          DO k = nzz, 1, -1
!            ! For partially rodded node, get volume weighted homogenized CX (0 < vfrac < 1._DP)
!            IF (rodh >= dum .AND. rodh <= dum+zdel(k)) THEN   ! If this node partially rodded
!               eta1 = rodh - dum
!               eta2 = zdel(k) - rodh + dum
!               IF (cusp == 0 .OR. eta1 < 1. .OR. eta2 < 1.) THEN    ! IF ROD CUSPING NOT ACTIVE OR LESS THAN 1 CM CLOSE TO Boundary
!                  vfrac = (rodh - dum) / zdel(k)
!                  sigtr(xyz(i,j,k),:) = sigtr(xyz(i,j,k),:) + &
!                                     vfrac * dsigtr(mat(xyz(i,j,k)),:)
!                  siga(xyz(i,j,k),:)  = siga(xyz(i,j,k),:) + &
!                                     vfrac * dsiga(mat(xyz(i,j,k)),:)
!                  nuf(xyz(i,j,k),:)   = nuf(xyz(i,j,k),:) + &
!                                     vfrac * dnuf(mat(xyz(i,j,k)),:)
!                  sigf(xyz(i,j,k),:)  = sigf(xyz(i,j,k),:) + &
!                                     vfrac * dsigf(mat(xyz(i,j,k)),:)
!                  sigs(xyz(i,j,k),:,:)  = sigs(xyz(i,j,k),:,:) + &
!                                       vfrac * dsigs(mat(xyz(i,j,k)),:,:)
!               ELSE                    ! IF ROD CUSPING ACTIVE
!                  n1 = CEILING(rodh - dum)        ! Number of mesh in rodded area
!                  del1 = (rodh - dum) / REAL(n1)  ! mesh size in rodded area
!                  n2 = CEILING(zdel(k) - rodh + dum)  ! Number of mesh in non-rodded area
!                  del2 = (zdel(k) - rodh + dum) / REAL(n2)  ! mesh size in non-rodded area
!
!                  nmax = n1 + n2                     ! Total number of mesh
!
!                  ! Calculate vectors a, b, c, d
!                  ALLOCATE(f(nmax))
!
!                  DO g = 1, ng
!                     ! Determine the flux coefficients
!                     a1 = 2. * (nod(xyz(i,j,k),g)%jo(5) + nod(xyz(i,j,k),g)%ji(5) &
!                        - nod(xyz(i,j,k),g)%jo(6) - nod(xyz(i,j,k),g)%ji(6))
!                     a2 = 2. * (nod(xyz(i,j,k),g)%jo(5) + nod(xyz(i,j,k),g)%ji(5) &
!                        + nod(xyz(i,j,k),g)%jo(6) + nod(xyz(i,j,k),g)%ji(6)) &
!                        - 2. * f0(xyz(i,j,k),g)
!                     a3 = 10. * a1 - 120. * fz1(xyz(i,j,k),g)
!                     a4 = 35. * a2 - 700. * fz2(xyz(i,j,k),g)
!
!                     ! Calculate fluxes in rodded area
!                     x = 0.5 * zdel(k)
!                     tx = x / zdel(k)
!                     f1 = f0(xyz(i,j,k),g) + a1 * tx + a2 * (3*tx**2-0.25) &
!                          + a3 * (tx*(tx+0.5)*(tx-0.5)) &
!                          + a4 * ((tx**2-0.05)*(tx+0.5)*(tx-0.5))
!                     DO n = 1, n1
!                        x = x - del1
!                        tx = x / zdel(k)
!                        f2 = f0(xyz(i,j,k),g) + a1 * tx + a2 * (3*tx**2-0.25) &
!                             + a3 * (tx*(tx+0.5)*(tx-0.5)) &
!                             + a4 * ((tx**2-0.05)*(tx+0.5)*(tx-0.5))
!                       f(n) = 0.5 * (f1 + f2)
!                       f1 = f2
!                    END DO
!                    ! Calculate fluxes in non-rodded area
!                    DO n = n1+1, nmax
!                       x = x - del2
!                       tx = x / zdel(k)
!                       f2 = f0(xyz(i,j,k),g) + a1 * tx + a2 * (3*tx**2-0.25) &
!                            + a3 * (tx*(tx+0.5)*(tx-0.5)) &
!                            + a4 * ((tx**2-0.05)*(tx+0.5)*(tx-0.5))
!                       f(n) = 0.5 * (f1 + f2)
!                       f1 = f2
!                    END DO
!
!                    ! Calculate homogenized CXs
!                    ! Rodded area
!                    sumx = 0.
!                    sum1 = 0.; sum2 = 0.; sum3 = 0.; sum4 = 0.; sum5 = 0.
!                    DO n = 1, n1
!                       sumx = sumx + f(n) * del1
!                       sum1 = sum1 + f(n) * (sigtr(xyz(i,j,k),g) &
!                       + dsigtr(mat(xyz(i,j,k)),g)) * del1
!                       sum2 = sum2 + f(n) * (siga(xyz(i,j,k),g) &
!                       + dsiga(mat(xyz(i,j,k)),g)) * del1
!                       sum3 = sum3 + f(n) * (nuf(xyz(i,j,k),g) &
!                       + dnuf(mat(xyz(i,j,k)),g)) * del1
!                       sum4 = sum4 + f(n) * (sigf(xyz(i,j,k),g) &
!                       + dsigf(mat(xyz(i,j,k)),g)) * del1
!                       DO h = 1, ng
!                           sum5(h) = sum5(h) + f(n) * (sigs(xyz(i,j,k),g,h) &
!                           + dsigs(mat(xyz(i,j,k)),g,h)) * del1
!                       END DO
!                    END DO
!                    ! Non-rodded area
!                    DO n = n1+1, nmax
!                       sumx = sumx + f(n) * del2
!                       sum1 = sum1 + f(n) * sigtr(xyz(i,j,k),g) * del2
!                       sum2 = sum2 + f(n) * siga(xyz(i,j,k),g) * del2
!                       sum3 = sum3 + f(n) * nuf(xyz(i,j,k),g) * del2
!                       sum4 = sum4 + f(n) * sigf(xyz(i,j,k),g) * del2
!                       DO h = 1, ng
!                           sum5(h) = sum5(h) + f(n) * sigs(xyz(i,j,k),g,h) * del2
!                       END DO
!                    END DO
!
!                    sigtr(xyz(i,j,k),g) = sum1 / sumx
!                    siga(xyz(i,j,k),g)  = sum2 / sumx
!                    nuf(xyz(i,j,k),g)   = sum3 / sumx
!                    sigf(xyz(i,j,k),g)  = sum4 / sumx
!                    DO h = 1, ng
!                       sigs(xyz(i,j,k),g,h) = sum5(h) / sumx
!                    END DO
!
!                 END DO
!                 DEALLOCATE(f)
!               END IF
!            ELSE IF (rodh > dum+zdel(k))   ! If fully rodded
!              CALL brInterpR(mat(xyz(i,j,k)), xcden(xyz(i,j,k)),  &
!              xbcon(xyz(i,j,k)), xftem(xyz(i,j,k)), xmtem(xyz(i,j,k)), &
!              sigtr, siga, nuf, sigf, sigs, dc)
!            ELSE                            ! If unrodded
!              CALL brInterp(mat(xyz(i,j,k)), xcden(xyz(i,j,k)),  &
!              xbcon(xyz(i,j,k)), xftem(xyz(i,j,k)), xmtem(xyz(i,j,k)), &
!              sigtr, siga, nuf, sigf, sigs, dc)
!            END IF
!          END DO
!      END IF
!   END DO
! END DO
!
!
! END SUBROUTINE xxs_updt

!******************************************************************************!

SUBROUTINE brInterp (rod, matnum, xcden, xbcon, xftem, xmtem, sigtr, siga, nuf, &
  sigf, sigs, dc)
!Purpose: To interpolate the xsec data from xtab file for given bcon,
! ftem, mtem and cden

USE sdata, ONLY: xmat, XBRANCH, ncden, nbcon, nftem, nmtem, ng

IMPLICIT NONE

INTEGER, INTENT(IN) :: rod, matnum  ! CR indicator and material number
REAL(DP), INTENT(IN) :: xbcon, xftem, xmtem, xcden  ! TH Parameters
REAL(DP), DIMENSION(:), INTENT(OUT) :: sigtr, siga, nuf, sigf
REAL(DP), DIMENSION(:,:), INTENT(OUT) :: dc, sigs

INTEGER, PARAMETER :: nx = 8
TYPE(XBRANCH), DIMENSION(nx) :: xs   !Temporary xsec for interpolation
INTEGER :: s, t, u, v
INTEGER :: s1=1, s2=1, t1=1, t2=1, u1=1, u2=1, v1=1, v2=1
LOGICAL :: fnd
INTEGER :: i

!Define + and - operators for XBRANCH type addition and substitution respectively
INTERFACE OPERATOR (+)
  MODULE PROCEDURE brAdd
END INTERFACE
INTERFACE OPERATOR (-)
  MODULE PROCEDURE brSubst
END INTERFACE
INTERFACE OPERATOR (*)
  MODULE PROCEDURE brRealMult
END INTERFACE

! Get 2 closest branch parameters (Given parameters should be in between)
fnd = .TRUE.
DO s = 2, ncden        ! For coolant density
  IF (xcden >= xmat(matnum)%pcden(s-1) .AND. xcden <= xmat(matnum)%pcden(s)) THEN
    s1 = s-1
    s2 = s
    fnd = .FALSE.
  END IF
END DO
IF (fnd .AND. ncden > 1) THEN   ! If not inside branch parameters, STOP
  WRITE(101,*) '  ERROR: COOLANT DENSITY IS OUT OF THE RANGE OF THE BRANCH PARAMETER'
  WRITE(*,*) '  ERROR: COOLANT DENSITY IS OUT OF THE RANGE OF THE BRANCH PARAMETER'
  STOP
END IF

fnd = .TRUE.
DO t = 2, nbcon      ! For Boron concentration
  IF (xbcon >= xmat(matnum)%pbcon(t-1) .AND. xbcon <= xmat(matnum)%pbcon(t)) THEN
    t1 = t-1
    t2 = t
    fnd = .FALSE.
  END IF
END DO
IF (fnd .AND. nbcon > 1) THEN  ! If not inside branch parameters, STOP
  WRITE(101,*) '  ERROR: BORON CONCENTRATION IS OUT OF THE RANGE OF THE BRANCH PARAMETER'
  WRITE(*,*) '  ERROR: BORON CONCENTRATION IS OUT OF THE RANGE OF THE BRANCH PARAMETER'
  STOP
END IF

fnd = .TRUE.
DO u = 2, nftem    ! For fuel temperature
  IF (xftem >= xmat(matnum)%pftem(u-1) .AND. xftem <= xmat(matnum)%pftem(u)) THEN
    u1 = u-1
    u2 = u
    fnd = .FALSE.
  END IF
END DO
IF (fnd .AND. nftem > 1) THEN    ! If not inside branch parameters, STOP
  WRITE(101,*) '  ERROR: FUEL TEMPERATURE IS OUT OF THE RANGE OF THE BRANCH PARAMETER'
  WRITE(*,*) '  ERROR: FUEL TEMPERATURE IS OUT OF THE RANGE OF THE BRANCH PARAMETER'
  STOP
END IF

fnd = .TRUE.
DO v = 2, nmtem    ! For moderator temperature
  IF (xmtem >= xmat(matnum)%pmtem(v-1) .AND. xmtem <= xmat(matnum)%pmtem(v)) THEN
    v1 = v-1
    v2 = v
    fnd = .FALSE.
  END IF
END DO
IF (fnd .AND. nmtem > 1) THEN    ! If not inside branch parameters, STOP
  WRITE(101,*) '  ERROR: MODERATOR TEMPERATURE IS OUT OF THE RANGE OF THE BRANCH PARAMETER'
  WRITE(*,*) '  ERROR: MODERATOR TEMPERATURE IS OUT OF THE RANGE OF THE BRANCH PARAMETER'
  STOP
END IF

!Start doing interpolation
DO i = 1, nx   !Allocate memory to xs
  ALLOCATE(xs(i)%sigtr(ng))
  ALLOCATE(xs(i)%siga(ng))
  ALLOCATE(xs(i)%nuf(ng))
  ALLOCATE(xs(i)%sigf(ng))
  ALLOCATE(xs(i)%dc(6,ng))
  ALLOCATE(xs(i)%sigs(ng,ng))
END DO

IF (rod == 0) THEN   ! For Unrodded XSEC
  !interpolation on Moderator Temperature
  IF (nmtem > 1) THEN
    xs(1) = xmat(matnum)%xsec(s1,t1,u1,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%xsec(s1,t1,u1,v2) - xmat(matnum)%xsec(s1,t1,u1,v1))
    xs(2) = xmat(matnum)%xsec(s1,t1,u2,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%xsec(s1,t1,u2,v2) - xmat(matnum)%xsec(s1,t1,u2,v1))
    xs(3) = xmat(matnum)%xsec(s1,t2,u1,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%xsec(s1,t2,u1,v2) - xmat(matnum)%xsec(s1,t2,u1,v1))
    xs(4) = xmat(matnum)%xsec(s1,t2,u2,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%xsec(s1,t2,u2,v2) - xmat(matnum)%xsec(s1,t2,u2,v1))
    xs(5) = xmat(matnum)%xsec(s2,t1,u1,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%xsec(s2,t1,u1,v2) - xmat(matnum)%xsec(s2,t1,u1,v1))
    xs(6) = xmat(matnum)%xsec(s2,t1,u2,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%xsec(s2,t1,u2,v2) - xmat(matnum)%xsec(s2,t1,u2,v1))
    xs(7) = xmat(matnum)%xsec(s2,t2,u1,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%xsec(s2,t2,u1,v2) - xmat(matnum)%xsec(s2,t2,u1,v1))
    xs(8) = xmat(matnum)%xsec(s2,t2,u2,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%xsec(s2,t2,u2,v2) - xmat(matnum)%xsec(s2,t2,u2,v1))
  ELSE
    xs(1) = xmat(matnum)%xsec(s1,t1,u1,v1)
    xs(2) = xmat(matnum)%xsec(s1,t1,u2,v1)
    xs(3) = xmat(matnum)%xsec(s1,t2,u1,v1)
    xs(4) = xmat(matnum)%xsec(s1,t2,u2,v1)
    xs(5) = xmat(matnum)%xsec(s2,t1,u1,v1)
    xs(6) = xmat(matnum)%xsec(s2,t1,u2,v1)
    xs(7) = xmat(matnum)%xsec(s2,t2,u1,v1)
    xs(8) = xmat(matnum)%xsec(s2,t2,u2,v1)
  END IF

  !interpolation on Fuel Temperature
  IF (nftem > 1) THEN
    xs(1) = xs(1) + (xftem - xmat(matnum)%pftem(u1)) / (xmat(matnum)%pftem(u2) - xmat(matnum)%pftem(u1)) * &
    (xs(2) - xs(1))
    xs(3) = xs(3) + (xftem - xmat(matnum)%pftem(u1)) / (xmat(matnum)%pftem(u2) - xmat(matnum)%pftem(u1)) * &
    (xs(4) - xs(3))
    xs(5) = xs(5) + (xftem - xmat(matnum)%pftem(u1)) / (xmat(matnum)%pftem(u2) - xmat(matnum)%pftem(u1)) * &
    (xs(6) - xs(5))
    xs(7) = xs(7) + (xftem - xmat(matnum)%pftem(u1)) / (xmat(matnum)%pftem(u2) - xmat(matnum)%pftem(u1)) * &
    (xs(8) - xs(7))
  END IF

  !interpolation on Boron concentration
  IF (nbcon > 1) THEN
    xs(1) = xs(1) + (xbcon - xmat(matnum)%pbcon(t1)) / (xmat(matnum)%pbcon(t2) - xmat(matnum)%pbcon(t1)) * &
    (xs(3) - xs(1))
    xs(5) = xs(5) + (xbcon - xmat(matnum)%pbcon(t1)) / (xmat(matnum)%pbcon(t2) - xmat(matnum)%pbcon(t1)) * &
    (xs(7) - xs(3))
  END IF

  !interpolation on coolant density
  IF (ncden > 1) THEN
    xs(1) = xs(1) + (xcden - xmat(matnum)%pcden(s1)) / (xmat(matnum)%pcden(s2) - xmat(matnum)%pcden(s1)) * &
    (xs(5) - xs(1))
  END IF
ELSE   ! For Rodded XSEC
  !interpolation on Moderator Temperature
  IF (nmtem > 1) THEN
    xs(1) = xmat(matnum)%rxsec(s1,t1,u1,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%rxsec(s1,t1,u1,v2) - xmat(matnum)%rxsec(s1,t1,u1,v1))
    xs(2) = xmat(matnum)%rxsec(s1,t1,u2,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%rxsec(s1,t1,u2,v2) - xmat(matnum)%rxsec(s1,t1,u2,v1))
    xs(3) = xmat(matnum)%rxsec(s1,t2,u1,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%rxsec(s1,t2,u1,v2) - xmat(matnum)%rxsec(s1,t2,u1,v1))
    xs(4) = xmat(matnum)%rxsec(s1,t2,u2,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%rxsec(s1,t2,u2,v2) - xmat(matnum)%rxsec(s1,t2,u2,v1))
    xs(5) = xmat(matnum)%rxsec(s2,t1,u1,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%rxsec(s2,t1,u1,v2) - xmat(matnum)%rxsec(s2,t1,u1,v1))
    xs(6) = xmat(matnum)%rxsec(s2,t1,u2,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%rxsec(s2,t1,u2,v2) - xmat(matnum)%rxsec(s2,t1,u2,v1))
    xs(7) = xmat(matnum)%rxsec(s2,t2,u1,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%rxsec(s2,t2,u1,v2) - xmat(matnum)%rxsec(s2,t2,u1,v1))
    xs(8) = xmat(matnum)%rxsec(s2,t2,u2,v1) + &
    (xmtem - xmat(matnum)%pmtem(v1)) / (xmat(matnum)%pmtem(v2) - xmat(matnum)%pmtem(v1)) * &
    (xmat(matnum)%rxsec(s2,t2,u2,v2) - xmat(matnum)%rxsec(s2,t2,u2,v1))
  ELSE
    xs(1) = xmat(matnum)%rxsec(s1,t1,u1,v1)
    xs(2) = xmat(matnum)%rxsec(s1,t1,u2,v1)
    xs(3) = xmat(matnum)%rxsec(s1,t2,u1,v1)
    xs(4) = xmat(matnum)%rxsec(s1,t2,u2,v1)
    xs(5) = xmat(matnum)%rxsec(s2,t1,u1,v1)
    xs(6) = xmat(matnum)%rxsec(s2,t1,u2,v1)
    xs(7) = xmat(matnum)%rxsec(s2,t2,u1,v1)
    xs(8) = xmat(matnum)%rxsec(s2,t2,u2,v1)
  END IF

  !interpolation on Fuel Temperature
  IF (nftem > 1) THEN
    xs(1) = xs(1) + (xftem - xmat(matnum)%pftem(u1)) / (xmat(matnum)%pftem(u2) - xmat(matnum)%pftem(u1)) * &
    (xs(2) - xs(1))
    xs(3) = xs(3) + (xftem - xmat(matnum)%pftem(u1)) / (xmat(matnum)%pftem(u2) - xmat(matnum)%pftem(u1)) * &
    (xs(4) - xs(3))
    xs(5) = xs(5) + (xftem - xmat(matnum)%pftem(u1)) / (xmat(matnum)%pftem(u2) - xmat(matnum)%pftem(u1)) * &
    (xs(6) - xs(5))
    xs(7) = xs(7) + (xftem - xmat(matnum)%pftem(u1)) / (xmat(matnum)%pftem(u2) - xmat(matnum)%pftem(u1)) * &
    (xs(8) - xs(7))
  END IF

  !interpolation on Boron concentration
  IF (nbcon > 1) THEN
    xs(1) = xs(1) + (xbcon - xmat(matnum)%pbcon(t1)) / (xmat(matnum)%pbcon(t2) - xmat(matnum)%pbcon(t1)) * &
    (xs(3) - xs(1))
    xs(5) = xs(5) + (xbcon - xmat(matnum)%pbcon(t1)) / (xmat(matnum)%pbcon(t2) - xmat(matnum)%pbcon(t1)) * &
    (xs(7) - xs(3))
  END IF

  !interpolation on coolant density
  IF (ncden > 1) THEN
    xs(1) = xs(1) + (xcden - xmat(matnum)%pcden(s1)) / (xmat(matnum)%pcden(s2) - xmat(matnum)%pcden(s1)) * &
    (xs(5) - xs(1))
  END IF
END IF

sigtr = xs(1)%sigtr
siga = xs(1)%siga
nuf = xs(1)%nuf
sigf = xs(1)%sigf
sigs = xs(1)%sigs
dc = xs(1)%dc

DO i = 1, nx   !DeAllocate memory to xs
  DEALLOCATE(xs(i)%sigtr)
  DEALLOCATE(xs(i)%siga)
  DEALLOCATE(xs(i)%nuf)
  DEALLOCATE(xs(i)%sigf)
  DEALLOCATE(xs(i)%dc)
  DEALLOCATE(xs(i)%sigs)
END DO


END SUBROUTINE brInterp

!******************************************************************************!

FUNCTION brAdd(A, B) RESULT (C)

  ! To perform XBRANCH data type addition

USE sdata, ONLY: XBRANCH

IMPLICIT NONE

TYPE(XBRANCH), INTENT(IN) :: A, B
TYPE(XBRANCH) :: C

C%sigtr = A%sigtr + B%sigtr
C%siga = A%siga + B%siga
C%nuf = A%nuf + B%nuf
C%sigf = A%sigf + B%sigf
C%sigs = A%sigs + B%sigs
C%dc = A%dc + B%dc

END FUNCTION brAdd

!******************************************************************************!

FUNCTION brSubst(A, B) RESULT (C)

    ! To perform XBRANCH data type substraction

USE sdata, ONLY: XBRANCH

IMPLICIT NONE

TYPE(XBRANCH), INTENT(IN) :: A, B
TYPE(XBRANCH) :: C

C%sigtr = A%sigtr - B%sigtr
C%siga = A%siga - B%siga
C%nuf = A%nuf - B%nuf
C%sigf = A%sigf - B%sigf
C%sigs = A%sigs - B%sigs
C%dc = A%dc - B%dc

END FUNCTION brSubst

!******************************************************************************!

FUNCTION brRealMult(Re, A) RESULT (B)

    ! To perform XBRANCH data type substraction

USE sdata, ONLY: XBRANCH, DP

IMPLICIT NONE

REAL(DP), INTENT(IN) :: Re
TYPE(XBRANCH), INTENT(IN) :: A
TYPE(XBRANCH) :: B

B%sigtr = Re * A%sigtr
B%siga = Re * A%siga
B%nuf = Re * A%nuf
B%sigf = Re * A%sigf
B%sigs = Re * A%sigs
B%dc = Re * A%dc

END FUNCTION brRealMult

!******************************************************************************!

SUBROUTINE openFile(iunit, iname, message)

IMPLICIT NONE

INTEGER :: iunit
CHARACTER(LEN=*) :: iname, message
INTEGER  :: iost

OPEN (UNIT=iunit, FILE=iname, STATUS='OLD', ACTION='READ', &
      IOSTAT = iost)

IF (iost /= 0) THEN
    WRITE(*,1020) message, iost
    WRITE(*,*) '  CANNOT OPEN FILE : ', iname
    1020 FORMAT    (2X, A, I6)
    STOP
END IF

END SUBROUTINE openFile


SUBROUTINE inp_comments (inunit, buffer, mark)
!
! Purpose:
!    To remove the comments in input and rewrite the
!    input into input buffer for each card. Comments marked by !.
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: inunit, buffer
CHARACTER(LEN=*), INTENT(IN) :: mark

INTEGER :: ln                  ! line number
INTEGER :: eof, comm

OPEN (UNIT=buffer, STATUS='SCRATCH', ACTION='READWRITE')

! Start removing comments and rewrite into one input buffer
ln = 0
DO
    ln = ln+1
    READ (inunit, '(A100)', IOSTAT=eof) iline
    IF (eof < 0) EXIT              !Check end of file
    iline = TRIM(ADJUSTL(iline))   ! Remove trailing blanks and adjust to left
    comm = INDEX(iline, mark)       ! Find position '!' if any
    ! If there is no '!' and no first ten blank spaces
    IF (comm == 0 .AND. iline(1:20) /= '                    ')  THEN
        WRITE(buffer,1012)'x ',ln,iline
    END IF
    !If the first character is not '!'
    IF (comm > 1) THEN
        iline = iline(1:comm-1)       ! Take only part of input
        WRITE(buffer,1012)'x ',ln, iline
    END IF
END DO

REWIND(buffer)

1012 FORMAT(A2, I5,' ',A100)

END SUBROUTINE inp_comments


SUBROUTINE er_message (funit, iost, ln, mess, xtab)
!
! Purpose:
!    To provide error message
!

IMPLICIT NONE

INTEGER, INTENT(IN) :: funit, iost, ln
INTEGER, OPTIONAL, INTENT(IN) :: xtab
CHARACTER(LEN=*), INTENT(IN) :: mess

IF (iost < 0) THEN
    IF (PRESENT(xtab)) THEN
      WRITE(funit, 1014) ln, xtab
    ELSE
      WRITE(funit, 1013) ln
    END IF
    WRITE(funit,*) mess
    IF (PRESENT(xtab)) THEN
      WRITE(*, 1014) ln, xtab
    ELSE
      WRITE(*, 1013) ln
    END IF
    WRITE(*,*) mess
    1013 FORMAT(2x, 'ERROR: Line', I4, ' needs more data')
    1014 FORMAT(2x, 'ERROR: Line', I4, &
    'in XTAB file for material number' , I4, '. It needs more data')
    STOP
END IF
IF (iost > 0) THEN
    IF (PRESENT(xtab)) THEN
       WRITE(funit, 1005) ln, xtab
    ELSE
       WRITE(funit, 1004) ln
    END IF
    WRITE(funit,*) mess
    IF (PRESENT(xtab)) THEN
      WRITE(*, 1005) ln, xtab
    ELSE
      WRITE(*, 1004) ln
    END IF
    WRITE(*,*) mess
    1004 FORMAT(2X, 'ERROR: Please check line number', I4)
    1005 FORMAT(2X, 'ERROR: Please check line number', I4, &
    ' in XTAB file for material number ', I4)
    STOP
END IF

END SUBROUTINE er_message


END MODULE InpOutp2
