MODULE nodal

IMPLICIT NONE

SAVE

CONTAINS

SUBROUTINE outer4

USE sdata, ONLY: ng, nnod, mat, ystag, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, Ke
USE InpOutp, ONLY: ounit

IMPLICIT NONE

REAL :: Keo                                ! new and old Multiplication factor (Keff)
REAL, DIMENSION(nnod) :: fs0, fs0c             ! new and old fission source, and scattering source
REAL, DIMENSION(nnod) :: fsx1, fsy1, fsz1
REAL, DIMENSION(nnod) :: fsx2, fsy2, fsz2      ! Fission source moments
REAL, DIMENSION(nnod) :: fsx1c, fsy1c, fsz1c
REAL, DIMENSION(nnod) :: fsx2c, fsy2c, fsz2c
REAL, DIMENSION(nnod) :: ss0                   ! Scattering source
REAL, DIMENSION(nnod) :: ssx1, ssy1, ssz1
REAL, DIMENSION(nnod) :: ssx2, ssy2, ssz2      ! Scattering source moments
REAL :: fer, Ker, ierr                           ! Fission source and Keff error
REAL :: f, fc                                         ! new and old integrated fission sources
REAL :: dr, omega                                            ! dominance ratio
INTEGER :: h, g
INTEGER :: lt
INTEGER :: fneg

REAL, DIMENSION(nnod) :: errn, erro

! K-eff and fission source guess
fs0  = 0.d0
fsx1 = 0.d0; fsy1 = 0.d0; fsz1 = 0.d0
fsx2 = 0.d0; fsy2 = 0.d0; fsz2 = 0.d0

DO g= 1, ng
    CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
END DO

f = Integrate(fs0)
!Start outer iteration
DO lt=1,999
    fc = f
    fs0c  = fs0
	fsx1c = fsx1; fsy1c = fsy1; fsz1c = fsz1 
	fsx2c = fsx2; fsy2c = fsy2; fsz2c = fsz2 
    fs0  = 0.d0
	fsx1 = 0.d0; fsy1 = 0.d0; fsz1 = 0.d0
	fsx2 = 0.d0; fsy2 = 0.d0; fsz2 = 0.d0
    Keo = Ke
    DO g = 1, ng
	
	    ss0  = 0.d0
        ssx1 = 0.d0; ssy1 = 0.d0; ssz1 = 0.d0
		ssx2 = 0.d0; ssy2 = 0.d0; ssz2 = 0.d0

		!!!Calculate Scattering source
        CALL SSrc(g, ss0, ssx1, ssy1, ssz1, ssx2, ssy2, ssz2)
		
		!!!Calculate total source
		CALL TSrc(g, Keo, fs0c , fsx1c, fsy1c, fsz1c, &
		                         fsx2c, fsy2c, fsz2c, &
                           ss0 , ssx1 , ssy1 , ssz1 , &
						         ssx2 , ssy2 , ssz2   )
								 
		!!!Inner Iteration
        CALL inner4(g, ierr)
		
		!!!Calculate fission source for next outer iteration
	    CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2) 
    END DO
	

	f = Integrate(fs0)
    Ke = Keo * f / fc                              ! Update Keff

    CALL RelE(fs0, fs0c, fer)                      ! Search maximum point wise total source Relative Error\
	
    Ker = ABS(Ke - Keo)                            ! Get Keff Abs Error
	
	CALL CNeg(f0, fneg)

    WRITE(ounit,'(I5,F13.6,2ES15.5)') lt, Ke, fer, ierr

    IF ((Ker < 1.d-5) .AND. (fer < 1.d-5) .AND. (ierr < 1.d-5)) EXIT
END DO

WRITE(ounit,*)
WRITE(ounit, '(A36,F9.6)') 'MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke

END SUBROUTINE outer4


SUBROUTINE outer2

USE sdata, ONLY: ng, nnod, mat, nod, ystag, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, Ke, &
				 nuf, sigs, chi
USE InpOutp, ONLY: ounit

IMPLICIT NONE

REAL :: Keo                                ! new and old Multiplication factor (Keff)
REAL, DIMENSION(nnod) :: fs0, fs0c             ! new and old fission source, and scattering source
REAL, DIMENSION(nnod) :: ss0                   ! Scattering source
REAL :: fer, Ker, ierr                         ! Fission source and Keff error
REAL :: f, fc                                  ! new and old integrated fission sources
REAL :: dr, omega                              ! dominance ratio
INTEGER :: h, g
INTEGER :: lt, n
INTEGER :: fneg
REAL :: st, fn

REAL, DIMENSION(nnod) :: errn, erro

! K-eff and fission source guess
Ke = 1.d0
fs0  = 0.d0

DO g = 1,ng
    DO n = 1, nnod
	    fs0(n)   = fs0(n)   + f0 (n,g) * nuf(n,g)
    END DO
END DO	

f = Integrate(fs0)
!Start outer iteration
DO lt=1,500
    CALL CPU_TIME(st)
    fc = f
    fs0c  = fs0
    fs0  = 0.d0
    Keo = Ke
    DO g = 1, ng
	
	    ss0  = 0.d0

		!!!Calculate Scattering source
        DO h = 1, ng
            DO n = 1, nnod
	            IF (g /= h) THEN
			        ss0(n)   = ss0(n)   + sigs(n,h,g) * f0(n,h)
		        END IF
	        END DO
        END DO
		
		!!!Calculate total source
        DO n = 1, nnod
	        nod(n,g)%Q(1) = chi(n,g) * fs0c(n)/Keo  + ss0(n)
        END DO
								 
		!!!Inner Iteration
        CALL inner2(g, ierr)
		
		!!!Calculate fission source for next outer iteration
        DO n = 1, nnod
	        fs0(n)   = fs0(n)   + f0 (n,g) * nuf(n,g)
        END DO
    END DO

	f = Integrate(fs0)
    Ke = Keo * f / fc                              ! Update Keff

    CALL RelE(fs0, fs0c, fer)                      ! Search maximum point wise total source Relative Error\
	
    Ker = ABS(Ke - Keo)                            ! Get Keff Abs Error
	
	CALL CNeg(f0, fneg)
	 
	 CALL CPU_TIME(fn)
	
    WRITE(ounit,'(I5,F15.8,3ES15.5, F8.3, I6)') lt, Ke, Ker, fer, ierr, fn-st, fneg

    IF ((Ker < 1.d-5) .AND. (fer < 1.d-5) .AND. (ierr < 1.d-5)) EXIT
END DO

WRITE(ounit,*)
WRITE(ounit, '(A36,F9.6)') 'MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke

END SUBROUTINE outer2



SUBROUTINE outer4Fx

USE sdata, ONLY: ng, nnod, mat, ystag, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, Ke
USE InpOutp, ONLY: ounit

IMPLICIT NONE

REAL :: Keo                                ! new and old Multiplication factor (Keff)
REAL, DIMENSION(nnod) :: fs0, fs0c             ! new and old fission source, and scattering source
REAL, DIMENSION(nnod) :: fsx1, fsy1, fsz1
REAL, DIMENSION(nnod) :: fsx2, fsy2, fsz2      ! Fission source moments
REAL, DIMENSION(nnod) :: fsx1c, fsy1c, fsz1c
REAL, DIMENSION(nnod) :: fsx2c, fsy2c, fsz2c
REAL, DIMENSION(nnod) :: ss0                   ! Scattering source
REAL, DIMENSION(nnod) :: ssx1, ssy1, ssz1
REAL, DIMENSION(nnod) :: ssx2, ssy2, ssz2      ! Scattering source moments
REAL :: fer, Ker, ierr                         ! Fission source and Keff error
INTEGER :: h, g
INTEGER :: lt
INTEGER :: fneg

REAL, DIMENSION(nnod) :: errn, erro

! K-eff and fission source guess
fs0  = 0.d0
fsx1 = 0.d0; fsy1 = 0.d0; fsz1 = 0.d0
fsx2 = 0.d0; fsy2 = 0.d0; fsz2 = 0.d0

DO g= 1, ng
    CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
END DO

!Start outer iteration
DO lt=1,500
    fs0c  = fs0
	fsx1c = fsx1; fsy1c = fsy1; fsz1c = fsz1 
	fsx2c = fsx2; fsy2c = fsy2; fsz2c = fsz2 
    fs0  = 0.d0
	fsx1 = 0.d0; fsy1 = 0.d0; fsz1 = 0.d0
	fsx2 = 0.d0; fsy2 = 0.d0; fsz2 = 0.d0
    DO g = 1, ng
	
	    ss0  = 0.d0
        ssx1 = 0.d0; ssy1 = 0.d0; ssz1 = 0.d0
		ssx2 = 0.d0; ssy2 = 0.d0; ssz2 = 0.d0

		!!!Calculate Scattering source
        CALL SSrc(g, ss0, ssx1, ssy1, ssz1, ssx2, ssy2, ssz2)
		
		!!!Calculate total source
		CALL TSrcFx(g, Keo, fs0c , fsx1c, fsy1c, fsz1c, &
		                         fsx2c, fsy2c, fsz2c, &
                           ss0 , ssx1 , ssy1 , ssz1 , &
						         ssx2 , ssy2 , ssz2   )
								 
		!!!Inner Iteration
        CALL inner4(g, ierr)
		
		!!!Calculate fission source for next outer iteration
	    CALL FSrc (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2) 
    END DO
	
    CALL RelE(fs0, fs0c, fer)                      ! Search maximum point wise total source Relative Error\
	
	CALL CNeg(f0, fneg)
	
    WRITE(ounit,'(I5,2ES15.5)') lt, fer, ierr

    IF ((fer < 1.d-5) .AND. (ierr < 1.d-5)) EXIT
END DO

CALL MultF(Ke)

WRITE(ounit,*)
WRITE(ounit, '(A36,F9.6)') 'MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke

END SUBROUTINE outer4Fx



SUBROUTINE outer4ad

USE sdata, ONLY: ng, nnod, mat, ystag, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, Ke
USE InpOutp, ONLY: ounit

IMPLICIT NONE

REAL :: Keo                                ! new and old Multiplication factor (Keff)
REAL, DIMENSION(nnod) :: fs0, fs0c             ! new and old fission source, and scattering source
REAL, DIMENSION(nnod) :: fsx1, fsy1, fsz1
REAL, DIMENSION(nnod) :: fsx2, fsy2, fsz2      ! Fission source moments
REAL, DIMENSION(nnod) :: fsx1c, fsy1c, fsz1c
REAL, DIMENSION(nnod) :: fsx2c, fsy2c, fsz2c
REAL, DIMENSION(nnod) :: ss0                   ! Scattering source
REAL, DIMENSION(nnod) :: ssx1, ssy1, ssz1
REAL, DIMENSION(nnod) :: ssx2, ssy2, ssz2      ! Scattering source moments
REAL :: fer, Ker, ierr                           ! Fission source and Keff error
REAL :: f, fc                                         ! new and old integrated fission sources
REAL :: dr, omega                                            ! dominance ratio
INTEGER :: h, g
INTEGER :: lt
INTEGER :: fneg

REAL, DIMENSION(nnod) :: errn, erro

! K-eff and fission source guess
fs0  = 0.d0
fsx1 = 0.d0; fsy1 = 0.d0; fsz1 = 0.d0
fsx2 = 0.d0; fsy2 = 0.d0; fsz2 = 0.d0

DO g= 1, ng
    CALL FSrcAd (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
END DO

f = Integrate(fs0)
!Start outer iteration
DO lt=1,500
    fc = f
    fs0c  = fs0
	fsx1c = fsx1; fsy1c = fsy1; fsz1c = fsz1 
	fsx2c = fsx2; fsy2c = fsy2; fsz2c = fsz2 
    fs0  = 0.d0
	fsx1 = 0.d0; fsy1 = 0.d0; fsz1 = 0.d0
	fsx2 = 0.d0; fsy2 = 0.d0; fsz2 = 0.d0
    Keo = Ke
    DO g = ng,1,-1
	
	    ss0  = 0.d0
        ssx1 = 0.d0; ssy1 = 0.d0; ssz1 = 0.d0
		ssx2 = 0.d0; ssy2 = 0.d0; ssz2 = 0.d0

		!!!Calculate Scattering source
        CALL SSrcAd(g, ss0, ssx1, ssy1, ssz1, ssx2, ssy2, ssz2)
		
		!!!Calculate total source
		CALL TSrcAd(g, Keo, fs0c , fsx1c, fsy1c, fsz1c, &
		                         fsx2c, fsy2c, fsz2c, &
                           ss0 , ssx1 , ssy1 , ssz1 , &
						         ssx2 , ssy2 , ssz2   )
								 
		!!!Inner Iteration
        CALL inner4(g, ierr)
		
		!!!Calculate fission source for next outer iteration
	    CALL FSrcAd (g, fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2) 
    END DO	

	f = Integrate(fs0)
    Ke = Keo * f / fc                              ! Update Keff

    CALL RelE(fs0, fs0c, fer)                      ! Search maximum point wise total source Relative Error\
	
    Ker = ABS(Ke - Keo)                            ! Get Keff Abs Error
	
	CALL CNeg(f0, fneg)
	
    WRITE(ounit,'(I5,F13.6,2ES15.5)') lt, Ke, fer, ierr

    IF ((Ker < 1.d-5) .AND. (fer < 1.d-5) .AND. (ierr < 1.d-5)) EXIT
END DO

WRITE(ounit,*)
WRITE(ounit, '(A36,F9.6)') 'MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke

END SUBROUTINE outer4ad



SUBROUTINE inner4(g,ierr)
!
! Purpose:
!   To perform inner iterations
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nod, nnod, xstag, ystag, &
                xleft, xrigt, yback, yfrnt, zbott, ztop, &
				f0, ix, iy, iz, xyz, nzz

IMPLICIT NONE

INTEGER, INTENT(IN) :: g
REAL, INTENT(OUT) :: ierr
INTEGER :: l, n, i, j
REAL, DIMENSION(6) :: bvec, qvec

! Transverse Leakage Moments(0, Lx1, Ly1, Lz2, Lx2, Ly2, Lz3)
REAL, DIMENSION(7) :: Lm  
               
! To store old fluxes
REAL, DIMENSION(nnod) :: flxc   

! Jot Nodals' outgoing currents+flux  (X+, X-, Y+, Y-, Z+, Z-)
! Jin Nodals' ingoing currents+source (X+, X-, Y+, Y-, Z+, Z-)

DO l = 1, 1   
	flxc = f0(:,g)	
    DO n = 1, nnod
		
	        IF (ix(n) == ystag(iy(n))%smax) THEN                          ! East (X+) BC
		        CALL bcond(xrigt, n, g, 1)
		    ELSE
		        nod(n,g)%ji(1) = nod(xyz( ix(n)+1, iy(n), iz(n) ), g)%jo(2)
            END IF
			
			IF (ix(n) == ystag(iy(n))%smin) THEN                          ! West (X-) BC
				CALL bcond(xleft, n, g, 2)
			ELSE
				nod(n,g)%ji(2) = nod(xyz( ix(n)-1, iy(n), iz(n) ), g)%jo(1)
            END IF

			IF (iy(n) == xstag(ix(n))%smax) THEN                          ! South (Y+) BC
                CALL bcond(yfrnt, n, g, 3)
			ELSE
				nod(n,g)%ji(3) = nod(xyz( ix(n), iy(n)+1, iz(n) ), g)%jo(4)
			END IF

			IF (iy(n) == xstag(ix(n))%smin) THEN                          ! North (Y-) BC
                CALL bcond(yback, n, g, 4)
			ELSE
				nod(n,g)%ji(4) = nod(xyz( ix(n), iy(n)-1, iz(n) ), g)%jo(3)
			END IF

			IF (iz(n) == nzz) THEN                                    ! Top (Z+) BC
                CALL bcond(ztop, n, g, 5)
			ELSE
				nod(n,g)%ji(5) = nod(xyz( ix(n), iy(n), iz(n)+1 ), g)%jo(6)
			END IF

			IF (iz(n) == 1) THEN                                      ! Bottom (Z-)BC
                CALL bcond(zbott, n, g, 6)
			ELSE
   				nod(n,g)%ji(6) = nod(xyz( ix(n), iy(n), iz(n)-1 ), g)%jo(5)
			END IF
		
				
			! Update transverse leakage moments
			CALL TLUpd (n, g, Lm) 				
				
			CALL matvec(nod(n,g)%P, nod(n,g)%ji, bvec)

			CALL matvec(nod(n,g)%R, nod(n,g)%Q - Lm, qvec)
				
            nod(n,g)%jo = qvec+bvec
		
            ! Update zeroth transverse leakages
		    CALL LxyzUpd(n,g)
		
			! Update flux and flux moments
			CALL FluxUpd4(n, g, Lm)	
	END DO
	
    DO n = nnod,1,-1
		
	        IF (ix(n) == ystag(iy(n))%smax) THEN                          ! East (X+) BC
		        CALL bcond(xrigt, n, g, 1)
		    ELSE
		        nod(n,g)%ji(1) = nod(xyz( ix(n)+1, iy(n), iz(n) ), g)%jo(2)
            END IF
			
			IF (ix(n) == ystag(iy(n))%smin) THEN                          ! West (X-) BC
				CALL bcond(xleft, n, g, 2)
			ELSE
				nod(n,g)%ji(2) = nod(xyz( ix(n)-1, iy(n), iz(n) ), g)%jo(1)
            END IF

			IF (iy(n) == xstag(ix(n))%smax) THEN                          ! South (Y+) BC
                CALL bcond(yfrnt, n, g, 3)
			ELSE
				nod(n,g)%ji(3) = nod(xyz( ix(n), iy(n)+1, iz(n) ), g)%jo(4)
			END IF

			IF (iy(n) == xstag(ix(n))%smin) THEN                          ! North (Y-) BC
                CALL bcond(yback, n, g, 4)
			ELSE
				nod(n,g)%ji(4) = nod(xyz( ix(n), iy(n)-1, iz(n) ), g)%jo(3)
			END IF

			IF (iz(n) == nzz) THEN                                    ! Top (Z+) BC
                CALL bcond(ztop, n, g, 5)
			ELSE
				nod(n,g)%ji(5) = nod(xyz( ix(n), iy(n), iz(n)+1 ), g)%jo(6)
			END IF

			IF (iz(n) == 1) THEN                                      ! Bottom (Z-)BC
                CALL bcond(zbott, n, g, 6)
			ELSE
   				nod(n,g)%ji(6) = nod(xyz( ix(n), iy(n), iz(n)-1 ), g)%jo(5)
			END IF
		
				
			! Update transverse leakage moments
			CALL TLUpd (n, g, Lm) 				
				
			CALL matvec(nod(n,g)%P, nod(n,g)%ji, bvec)

			CALL matvec(nod(n,g)%R, nod(n,g)%Q - Lm, qvec)
				
            nod(n,g)%jo = qvec+bvec
		
            ! Update zeroth transverse leakages
		    CALL LxyzUpd(n,g)
		
			! Update flux and flux moments
			CALL FluxUpd4(n, g, Lm)	
	END DO
	
	CALL RelE(f0(:,g), flxc, ierr)
	
	
END DO  

END SUBROUTINE inner4




SUBROUTINE inner4x(g,ierr)
!
! Purpose:
!   To perform inner iterations
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nod, nnod, xstag, ystag, &
                xleft, xrigt, yback, yfrnt, zbott, ztop, &
				f0, ix, iy, iz, xyz, nzz, nxx, nyy

IMPLICIT NONE

INTEGER, INTENT(IN) :: g
REAL, INTENT(OUT) :: ierr
INTEGER :: l, n, i, j, k
REAL, DIMENSION(6) :: bvec, qvec

! Transverse Leakage Moments(0, Lx1, Ly1, Lz2, Lx2, Ly2, Lz3)
REAL, DIMENSION(7) :: Lm  
               
! To store old fluxes
REAL, DIMENSION(nnod) :: flxc   

! Jot Nodals' outgoing currents+flux  (X+, X-, Y+, Y-, Z+, Z-)
! Jin Nodals' ingoing currents+source (X+, X-, Y+, Y-, Z+, Z-)

DO l = 1, 1    
	flxc = f0(:,g)	
    DO k = 1, nzz
	DO j = 1, nyy
	DO i = ystag(j)%smin, ystag(j)%smax
		
	        IF (i == ystag(j)%smax) THEN                          ! East (X+) BC
		        CALL bcond(xrigt, xyz(i,j,k), g, 1)
		    ELSE
		        nod(xyz(i,j,k),g)%ji(1) = nod(xyz( i+1, j, k ), g)%jo(2)
            END IF
			
			IF (i == ystag(j)%smin) THEN                          ! West (X-) BC
				CALL bcond(xleft, xyz(i,j,k), g, 2)
			ELSE
				nod(xyz(i,j,k),g)%ji(2) = nod(xyz( i-1, j, k ), g)%jo(1)
            END IF

			IF (j == xstag(i)%smax) THEN                          ! South (Y+) BC
                CALL bcond(yfrnt, xyz(i,j,k), g, 3)
			ELSE
				nod(xyz(i,j,k),g)%ji(3) = nod(xyz( i , j+1, k ), g)%jo(4)
			END IF

			IF (j == xstag(i)%smin) THEN                          ! North (Y-) BC
                CALL bcond(yback, xyz(i,j,k), g, 4)
			ELSE
				nod(xyz(i,j,k),g)%ji(4) = nod(xyz( i , j-1, k ), g)%jo(3)
			END IF

			IF (k == nzz) THEN                                    ! Top (Z+) BC
                CALL bcond(ztop, xyz(i,j,k), g, 5)
			ELSE
				nod(xyz(i,j,k),g)%ji(5) = nod(xyz( i , j , k+1 ), g)%jo(6)
			END IF

			IF (k == 1) THEN                                      ! Bottom (Z-)BC
                CALL bcond(zbott, xyz(i,j,k), g, 6)
			ELSE
   				nod(xyz(i,j,k),g)%ji(6) = nod(xyz( i , j , k-1 ), g)%jo(5)
			END IF
		
				
			! Update transverse leakage moments
			CALL TLUpd (xyz(i,j,k), g, Lm) 				
				
			CALL matvec(nod(xyz(i,j,k),g)%P, nod(xyz(i,j,k),g)%ji, bvec)

			CALL matvec(nod(xyz(i,j,k),g)%R, nod(xyz(i,j,k),g)%Q - Lm, qvec)
				
            nod(xyz(i,j,k),g)%jo = qvec+bvec
		
            ! Update zeroth transverse leakages
		    CALL LxyzUpd(xyz(i,j,k),g)
		
			! Update flux and flux moments
			CALL FluxUpd4(xyz(i,j,k), g, Lm)	
	END DO
	END DO
	END DO
	
    DO k = 1, nzz
	DO j = nyy, 1, -1
	DO i = ystag(j)%smin, ystag(j)%smax
		
	        IF (i == ystag(j)%smax) THEN                          ! East (X+) BC
		        CALL bcond(xrigt, xyz(i,j,k), g, 1)
		    ELSE
		        nod(xyz(i,j,k),g)%ji(1) = nod(xyz( i+1, j, k ), g)%jo(2)
            END IF
			
			IF (i == ystag(j)%smin) THEN                          ! West (X-) BC
				CALL bcond(xleft, xyz(i,j,k), g, 2)
			ELSE
				nod(xyz(i,j,k),g)%ji(2) = nod(xyz( i-1, j, k ), g)%jo(1)
            END IF

			IF (j == xstag(i)%smax) THEN                          ! South (Y+) BC
                CALL bcond(yfrnt, xyz(i,j,k), g, 3)
			ELSE
				nod(xyz(i,j,k),g)%ji(3) = nod(xyz( i , j+1, k ), g)%jo(4)
			END IF

			IF (j == xstag(i)%smin) THEN                          ! North (Y-) BC
                CALL bcond(yback, xyz(i,j,k), g, 4)
			ELSE
				nod(xyz(i,j,k),g)%ji(4) = nod(xyz( i , j-1, k ), g)%jo(3)
			END IF

			IF (k == nzz) THEN                                    ! Top (Z+) BC
                CALL bcond(ztop, xyz(i,j,k), g, 5)
			ELSE
				nod(xyz(i,j,k),g)%ji(5) = nod(xyz( i , j , k+1 ), g)%jo(6)
			END IF

			IF (k == 1) THEN                                      ! Bottom (Z-)BC
                CALL bcond(zbott, xyz(i,j,k), g, 6)
			ELSE
   				nod(xyz(i,j,k),g)%ji(6) = nod(xyz( i , j , k-1 ), g)%jo(5)
			END IF
		
				
			! Update transverse leakage moments
			CALL TLUpd (xyz(i,j,k), g, Lm) 				
				
			CALL matvec(nod(xyz(i,j,k),g)%P, nod(xyz(i,j,k),g)%ji, bvec)

			CALL matvec(nod(xyz(i,j,k),g)%R, nod(xyz(i,j,k),g)%Q - Lm, qvec)
				
            nod(xyz(i,j,k),g)%jo = qvec+bvec
		
            ! Update zeroth transverse leakages
		    CALL LxyzUpd(xyz(i,j,k),g)
		
			! Update flux and flux moments
			CALL FluxUpd4(xyz(i,j,k), g, Lm)	
	END DO
	END DO
	END DO
	
    DO k = 1, nzz
	DO j = nyy, 1, -1
	DO i = ystag(j)%smax, ystag(j)%smin, -1
		
	        IF (i == ystag(j)%smax) THEN                          ! East (X+) BC
		        CALL bcond(xrigt, xyz(i,j,k), g, 1)
		    ELSE
		        nod(xyz(i,j,k),g)%ji(1) = nod(xyz( i+1, j, k ), g)%jo(2)
            END IF
			
			IF (i == ystag(j)%smin) THEN                          ! West (X-) BC
				CALL bcond(xleft, xyz(i,j,k), g, 2)
			ELSE
				nod(xyz(i,j,k),g)%ji(2) = nod(xyz( i-1, j, k ), g)%jo(1)
            END IF

			IF (j == xstag(i)%smax) THEN                          ! South (Y+) BC
                CALL bcond(yfrnt, xyz(i,j,k), g, 3)
			ELSE
				nod(xyz(i,j,k),g)%ji(3) = nod(xyz( i , j+1, k ), g)%jo(4)
			END IF

			IF (j == xstag(i)%smin) THEN                          ! North (Y-) BC
                CALL bcond(yback, xyz(i,j,k), g, 4)
			ELSE
				nod(xyz(i,j,k),g)%ji(4) = nod(xyz( i , j-1, k ), g)%jo(3)
			END IF

			IF (k == nzz) THEN                                    ! Top (Z+) BC
                CALL bcond(ztop, xyz(i,j,k), g, 5)
			ELSE
				nod(xyz(i,j,k),g)%ji(5) = nod(xyz( i , j , k+1 ), g)%jo(6)
			END IF

			IF (k == 1) THEN                                      ! Bottom (Z-)BC
                CALL bcond(zbott, xyz(i,j,k), g, 6)
			ELSE
   				nod(xyz(i,j,k),g)%ji(6) = nod(xyz( i , j , k-1 ), g)%jo(5)
			END IF
		
				
			! Update transverse leakage moments
			CALL TLUpd (xyz(i,j,k), g, Lm) 				
				
			CALL matvec(nod(xyz(i,j,k),g)%P, nod(xyz(i,j,k),g)%ji, bvec)

			CALL matvec(nod(xyz(i,j,k),g)%R, nod(xyz(i,j,k),g)%Q - Lm, qvec)
				
            nod(xyz(i,j,k),g)%jo = qvec+bvec
		
            ! Update zeroth transverse leakages
		    CALL LxyzUpd(xyz(i,j,k),g)
		
			! Update flux and flux moments
			CALL FluxUpd4(xyz(i,j,k), g, Lm)	
	END DO
	END DO
	END DO
	
    DO k = 1, nzz
	DO j = 1, nyy
	DO i = ystag(j)%smax, ystag(j)%smin, -1
		
	        IF (i == ystag(j)%smax) THEN                          ! East (X+) BC
		        CALL bcond(xrigt, xyz(i,j,k), g, 1)
		    ELSE
		        nod(xyz(i,j,k),g)%ji(1) = nod(xyz( i+1, j, k ), g)%jo(2)
            END IF
			
			IF (i == ystag(j)%smin) THEN                          ! West (X-) BC
				CALL bcond(xleft, xyz(i,j,k), g, 2)
			ELSE
				nod(xyz(i,j,k),g)%ji(2) = nod(xyz( i-1, j, k ), g)%jo(1)
            END IF

			IF (j == xstag(i)%smax) THEN                          ! South (Y+) BC
                CALL bcond(yfrnt, xyz(i,j,k), g, 3)
			ELSE
				nod(xyz(i,j,k),g)%ji(3) = nod(xyz( i , j+1, k ), g)%jo(4)
			END IF

			IF (j == xstag(i)%smin) THEN                          ! North (Y-) BC
                CALL bcond(yback, xyz(i,j,k), g, 4)
			ELSE
				nod(xyz(i,j,k),g)%ji(4) = nod(xyz( i , j-1, k ), g)%jo(3)
			END IF

			IF (k == nzz) THEN                                    ! Top (Z+) BC
                CALL bcond(ztop, xyz(i,j,k), g, 5)
			ELSE
				nod(xyz(i,j,k),g)%ji(5) = nod(xyz( i , j , k+1 ), g)%jo(6)
			END IF

			IF (k == 1) THEN                                      ! Bottom (Z-)BC
                CALL bcond(zbott, xyz(i,j,k), g, 6)
			ELSE
   				nod(xyz(i,j,k),g)%ji(6) = nod(xyz( i , j , k-1 ), g)%jo(5)
			END IF
		
				
			! Update transverse leakage moments
			CALL TLUpd (xyz(i,j,k), g, Lm) 				
				
			CALL matvec(nod(xyz(i,j,k),g)%P, nod(xyz(i,j,k),g)%ji, bvec)

			CALL matvec(nod(xyz(i,j,k),g)%R, nod(xyz(i,j,k),g)%Q - Lm, qvec)
				
            nod(xyz(i,j,k),g)%jo = qvec+bvec
		
            ! Update zeroth transverse leakages
		    CALL LxyzUpd(xyz(i,j,k),g)
		
			! Update flux and flux moments
			CALL FluxUpd4(xyz(i,j,k), g, Lm)	
	END DO
	END DO
	END DO
	
	
	
	CALL RelE(f0(:,g), flxc, ierr)
	
	
END DO  

END SUBROUTINE inner4x



SUBROUTINE inner2(g,ierr)
!
! Purpose:
!   To perform inner iterations
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nod, nnod, xstag, ystag, &
                xleft, xrigt, yback, yfrnt, zbott, ztop, &
				f0, ix, iy, iz, xyz, nzz

IMPLICIT NONE

INTEGER, INTENT(IN) :: g
REAL, INTENT(OUT) :: ierr
INTEGER :: l, n, i, j
REAL, DIMENSION(6) :: bvec, qvec
               
! To store old fluxes
REAL, DIMENSION(nnod) :: flxc   

! Jot Nodals' outgoing currents+flux  (X+, X-, Y+, Y-, Z+, Z-)
! Jin Nodals' ingoing currents+source (X+, X-, Y+, Y-, Z+, Z-)

DO l = 1, 2    ! 10 is maximum inner iteration	
	flxc = f0(:,g)	
    DO n = 1, nnod
		
	        IF (ix(n) == ystag(iy(n))%smax) THEN                          ! Right (X+) BC
		        CALL bcond(xrigt, n, g, 1)
		    ELSE
		        nod(n,g)%ji(1) = nod(xyz( ix(n)+1, iy(n), iz(n) ), g)%jo(2)
            END IF
			
			IF (ix(n) == ystag(iy(n))%smin) THEN                          ! Left (X-) BC
				CALL bcond(xleft, n, g, 2)
			ELSE
				nod(n,g)%ji(2) = nod(xyz( ix(n)-1, iy(n), iz(n) ), g)%jo(1)
            END IF

			IF (iy(n) == xstag(ix(n))%smax) THEN                          ! Front (Y+) BC
                CALL bcond(yfrnt, n, g, 3)
			ELSE
				nod(n,g)%ji(3) = nod(xyz( ix(n), iy(n)+1, iz(n) ), g)%jo(4)
			END IF

			IF (iy(n) == xstag(ix(n))%smin) THEN                          ! Back (Y-) BC
                CALL bcond(yback, n, g, 4)
			ELSE
				nod(n,g)%ji(4) = nod(xyz( ix(n), iy(n)-1, iz(n) ), g)%jo(3)
			END IF

			IF (iz(n) == nzz) THEN                                    ! Top (Z+) BC
                CALL bcond(ztop, n, g, 5)
			ELSE
				nod(n,g)%ji(5) = nod(xyz( ix(n), iy(n), iz(n)+1 ), g)%jo(6)
			END IF

			IF (iz(n) == 1) THEN                                      ! Bottom (Z-)BC
                CALL bcond(zbott, n, g, 6)
			ELSE
   				nod(n,g)%ji(6) = nod(xyz( ix(n), iy(n), iz(n)-1 ), g)%jo(5)
			END IF			
				
			CALL matvec(nod(n,g)%P, nod(n,g)%ji, bvec)

			CALL matvec(nod(n,g)%R, nod(n,g)%Q, qvec)
				
            nod(n,g)%jo = qvec+bvec
		
            ! Update zeroth transverse leakages
		    CALL LxyzUpd(n,g)
		
			! Update flux and flux moments
			CALL FluxUpd2(n, g)	
	END DO
	
	CALL RelE(f0(:,g), flxc, ierr)
	
END DO			    

END SUBROUTINE inner2



SUBROUTINE bcond (bc, nt, gt, side)

!
! Purpose:
!    To provide proper boundary conditions
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nod

IMPLICIT NONE

INTEGER, INTENT(IN) :: bc, nt, gt, side

IF (bc == 0) THEN
    nod(nt,gt)%ji(side) = 0.d0
ELSE
    nod(nt,gt)%ji(side) = nod(nt,gt)%jo(side)
END IF

END SUBROUTINE bcond



SUBROUTINE FluxUpd4 (n, g, L)

USE sdata, ONLY: nod, D, sigr, xdel, ydel, zdel, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, &
				 ix, iy, iz

! Purpose:
   ! To update nod averaged flux and flux moments

  ! Date                Programmer           History
 ! ========================================================
 ! 6 FEB 2017         Muhammad Imron       Original code

IMPLICIT NONE

INTEGER, INTENT(IN) :: g, n
REAL, DIMENSION(:), INTENT(IN) :: L

REAL :: Tx, Ty, Tz

! Calculate Zeroth Flux
f0(n,g)  = ( nod(n,g)%Q(1)          &
                 - nod(n,g)%L(1)/xdel(ix(n))   &
                 - nod(n,g)%L(2)/ydel(iy(n))   &
				 - nod(n,g)%L(3)/zdel(iz(n)))  &
				 / sigr(n,g)

! Set parameters Tx, Ty and Tz				 
Tx = nod(n,g)%jo(1) - nod(n,g)%ji(1) &
   - nod(n,g)%jo(2) + nod(n,g)%ji(2)
Ty = nod(n,g)%jo(3) - nod(n,g)%ji(3) &
   - nod(n,g)%jo(4) + nod(n,g)%ji(4)
Tz = nod(n,g)%jo(5) - nod(n,g)%ji(5) &
   - nod(n,g)%jo(6) + nod(n,g)%ji(6)

! Calculate Flux moments   
fx1(n,g) = ( nod(n,g)%Q(2) - L(2)               &
           - 0.5d0*Tx/xdel(ix(n))                &
		   - 2.d0*D(n,g)/xdel(ix(n))**2          &
		   * (nod(n,g)%jo(1) + nod(n,g)%ji(1)    &
           -  nod(n,g)%jo(2) - nod(n,g)%ji(2)) ) &
		   / sigr(n,g)

				 
fy1(n,g) = ( nod(n,g)%Q(3) - L(3)               &
           - 0.5d0*Ty/ydel(iy(n))                &
		   - 2.d0*D(n,g)/ydel(iy(n))**2          &
		   * (nod(n,g)%jo(3) + nod(n,g)%ji(3)    &
           -  nod(n,g)%jo(4) - nod(n,g)%ji(4)) ) &
		   / sigr(n,g)

				 
fz1(n,g) = ( nod(n,g)%Q(4) - L(4)               &
           - 0.5d0*Tz/zdel(iz(n))                &
		   - 2.d0*D(n,g)/zdel(iz(n))**2          &
		   * (nod(n,g)%jo(5) + nod(n,g)%ji(5)    &
           -  nod(n,g)%jo(6) - nod(n,g)%ji(6)) ) &
		   / sigr(n,g)


fx2(n,g) = ( nod(n,g)%Q(5) - L(5)               &
           - 0.5d0*nod(n,g)%L(1)/xdel(ix(n))     &
		   - 6.d0*D(n,g)/xdel(ix(n))**2          &
		   * (nod(n,g)%jo(1) + nod(n,g)%ji(1)    &
           +  nod(n,g)%jo(2) + nod(n,g)%ji(2)    &
	       - f0(n,g)) ) / sigr(n,g)


fy2(n,g) = ( nod(n,g)%Q(6) - L(6)               &
           - 0.5d0*nod(n,g)%L(2)/ydel(iy(n))     &
		   - 6.d0*D(n,g)/ydel(iy(n))**2          &
		   * (nod(n,g)%jo(3) + nod(n,g)%ji(3)    &
           +  nod(n,g)%jo(4) + nod(n,g)%ji(4)    &
		   - f0(n,g)) ) / sigr(n,g)


fz2(n,g) = ( nod(n,g)%Q(7) - L(7)               &
           - 0.5d0*nod(n,g)%L(3)/zdel(iz(n))     &
		   - 6.d0*D(n,g)/zdel(iz(n))**2          &
		   * (nod(n,g)%jo(5) + nod(n,g)%ji(5)    &
           +  nod(n,g)%jo(6) + nod(n,g)%ji(6)    &
		   - f0(n,g)) ) / sigr(n,g)


END SUBROUTINE FluxUpd4


SUBROUTINE FluxUpd2 (nt, gt)

USE sdata, ONLY: nod, D, sigr, xdel, ydel, zdel, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2, &
				 ix, iy, iz

! Purpose:
   ! To update nod averaged flux and flux moments

  ! Date                Programmer           History
 ! ========================================================
 ! 6 FEB 2017         Muhammad Imron       Original code

IMPLICIT NONE

INTEGER, INTENT(IN) :: gt, nt

REAL :: Tx, Ty, Tz

! Calculate Zeroth Flux
f0(nt,gt)  = ( nod(nt,gt)%Q(1)          &
                 - nod(nt,gt)%L(1)/xdel(ix(nt))   &
                 - nod(nt,gt)%L(2)/ydel(iy(nt))   &
				 - nod(nt,gt)%L(3)/zdel(iz(nt)))  &
				 / sigr(nt,gt)

END SUBROUTINE FluxUpd2



SUBROUTINE TLUpd (n, g, L)

USE sdata, ONLY: nod, xdel, ydel, zdel, xstag, ystag, nzz, &
				 xleft, xrigt, yback, yfrnt, zbott, ztop, &
				 ix, iy, iz, xyz

! Purpose:
   ! To calaculate transverse leakage moments

  ! Date                Programmer           History
 ! ========================================================
 ! 6 FEB 2017         Muhammad Imron       Original code

IMPLICIT NONE

INTEGER, INTENT(IN) :: g, n
REAL, DIMENSION(:), INTENT(OUT) :: L

REAL :: txm, txp, tym, typ, tzm, tzp
REAL :: p1m, p2m, p1p, p2p, p1, p2
REAL :: r1xy, r2xy, r1xz, r2xz
REAL :: r1yx, r2yx, r1yz, r2yz
REAL :: r1zx, r2zx, r1zy, r2zy

! Set paramaters for X-Direction Transverse leakage
IF (ix(n) == ystag(iy(n))%smin) THEN
    IF (xleft == 0) THEN
        txp = xdel(ix(n)+1)/xdel(ix(n))
	    p1p = txp+1.d0
        r1xy = 2.d0 *                     &
		     ( nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(2)   &
	         - nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)   &
			 ) / p1p
        r2xy = 0.d0
        r1xz = 2.d0 *                     &
		     ( nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(3)   &
             - nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)   &
			 ) / p1p
        r2xz = 0.d0
	ELSE
        txm = 1.d0
	    txp = xdel(ix(n)+1)/xdel(ix(n))
	    p1m = txm+1.d0; p2m = 2.d0*txm+1.d0; p1p = txp+1.d0 
	    p2p = 2.d0*txp+1.d0; p1 = txm+txp+1.d0; p2 = txm+txp+2.d0
        r1xy = (   p1m*p2m * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(2)   &
                 - p1p*p2p * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)   &
	             + (p1p*p2p-p1m*p2m)                                &
			     * nod(n,g)%L(2)                                    &
		       ) / (p1m*p1p*p1)
        r2xy = (   p1m * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(2)  &
                 + p1p * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)  &
	             - p2  * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)  & 
		       ) / (p1m*p1p*p1)
        r1xz = (   p1m*p2m * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(3)   &
                 - p1p*p2p * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)   &
	             + (p1p*p2p-p1m*p2m)                                &
			     * nod(n,g)%L(3)                                    &
		       ) / (p1m*p1p*p1)
        r2xz = (   p1m * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(3)  &
                 + p1p * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)  &
	             - p2  * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)  & 
		       ) / (p1m*p1p*p1)
    END IF
ELSE IF (ix(n) == ystag(iy(n))%smax) THEN
    IF (xrigt == 0) THEN
        txm = xdel(ix(n)-1)/xdel(ix(n))
	    p1m = txm+1.d0
        r1xy = 2.d0 *                                 &
		     ( nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)   &
	         - nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(2)   &
			 ) / p1m
        r2xy = 0.d0
        r1xz = 2.d0 *                                 &
		     ( nod(xyz(ix(n)  ,iy(n) ,iz(n)),g)%L(3)  &
             - nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(3)   &
			 ) / p1m
        r2xz = 0.d0
	ELSE
        txm = xdel(ix(n)-1)/xdel(ix(n))
	    txp = 1.d0
	    p1m = txm+1.d0; p2m = 2.d0*txm+1.d0;p1p = txp+1.d0 
	    p2p = 2.d0*txp+1.d0; p1 = txm+txp+1.d0; p2 = txm+txp+2.d0
        r1xy = (   p1m*p2m * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)   &
                 - p1p*p2p * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(2)   &
	             + (p1p*p2p-p1m*p2m)                                &
			     * nod(n,g)%L(2)                                    &
		       ) / (p1m*p1p*p1)
        r2xy = (   p1m * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)  &
                 + p1p * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(2)  &
	             - p2  * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)  & 
		       ) / (p1m*p1p*p1)
        r1xz = (   p1m*p2m * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)   &
                 - p1p*p2p * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(3)   &
	             + (p1p*p2p-p1m*p2m)                                &
			     * nod(n,g)%L(3)                                    &
		         ) / (p1m*p1p*p1)
        r2xz = (   p1m * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)  &
                 + p1p * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(3)  &
	             - p2  * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)  & 
		       ) / (p1m*p1p*p1)
    END IF	
ELSE
    txm = xdel(ix(n)-1)/xdel(ix(n))
	txp = xdel(ix(n)+1)/xdel(ix(n))
	p1m = txm+1.d0; p2m = 2.d0*txm+1.d0;p1p = txp+1.d0 
	p2p = 2.d0*txp+1.d0; p1 = txm+txp+1.d0; p2 = txm+txp+2.d0
    r1xy = (   p1m*p2m * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(2)   &
             - p1p*p2p * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(2)   &
	         + (p1p*p2p-p1m*p2m)                                &
			 * nod(n,g)%L(2)                                    &
		   ) / (p1m*p1p*p1)
    r2xy = (   p1m * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(2)  &
             + p1p * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(2)  &
	         - p2  * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(2)  & 
		   ) / (p1m*p1p*p1)
    r1xz = (   p1m*p2m * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(3)   &
             - p1p*p2p * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(3)   &
	         + (p1p*p2p-p1m*p2m)                                &
			 * nod(n,g)%L(3)                                    &
		   ) / (p1m*p1p*p1)
    r2xz = (   p1m * nod(xyz(ix(n)+1,iy(n),iz(n)),g)%L(3)  &
             + p1p * nod(xyz(ix(n)-1,iy(n),iz(n)),g)%L(3)  &
	         - p2  * nod(xyz(ix(n)  ,iy(n),iz(n)),g)%L(3)  & 
		   ) / (p1m*p1p*p1)
END IF


! Set paramaters for Y-Direction Transverse leakage
IF (iy(n) == xstag(ix(n))%smin) THEN
    IF (yback == 0) THEN
	    typ = ydel(iy(n)+1)/ydel(iy(n))
	    p1p = typ+1.d0
        r1yx = 2.d0 *                                 &
		     ( nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(1)   &
	         - nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)   &
			 ) / p1p
        r2yx = 0.d0
        r1yz = 2.d0 *                                 &
		     ( nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(3)   &
             - nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)   &
			 ) / p1p
        r2yz = 0.d0
	ELSE
        tym = 1.d0
        typ = ydel(iy(n)+1)/ydel(iy(n))
	    p1m = tym+1.d0; p2m = 2.d0*tym+1.d0;p1p = typ+1.d0
	    p2p = 2.d0*typ+1.d0; p1 = tym+typ+1.d0; p2 = tym+typ+2.d0
        r1yx = (   p1m*p2m * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(1)   &
                 - p1p*p2p * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)   &
	             + (p1p*p2p-p1m*p2m)                                &
			     * nod(n,g)%L(1)                                    &
		       ) / (p1m*p1p*p1)
        r2yx = (   p1m * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(1)  &
                 + p1p * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)  &
	             - p2  * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)  & 
		       ) / (p1m*p1p*p1)
        r1yz = (   p1m*p2m * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(3)   &
                 - p1p*p2p * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)   &
	             + (p1p*p2p-p1m*p2m)                                &
			     * nod(n,g)%L(3)                                    &
		       ) / (p1m*p1p*p1)
        r2yz = (   p1m * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(3)  &
                 + p1p * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)  &
	             - p2  * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)  & 
		       ) / (p1m*p1p*p1)	
    END IF
ELSE IF (iy(n) == xstag(ix(n))%smax) THEN
    IF (yfrnt == 0) THEN
        tym = ydel(iy(n)-1)/ydel(iy(n))
	    p1m = tym+1.d0
        r1yx = 2.d0 *                                 &
		     ( nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)   &
	         - nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(1)   &
			 ) / p1m
        r2yx = 0.d0
        r1yz = 2.d0 *                                 &
		     ( nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)   &
             - nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(3)   &
			 ) / p1m
        r2yz = 0.d0
	ELSE
        tym = ydel(iy(n)-1)/ydel(iy(n))
        typ = 1.d0
	    p1m = tym+1.d0; p2m = 2.d0*tym+1.d0;p1p = typ+1.d0
	    p2p = 2.d0*typ+1.d0; p1 = tym+typ+1.d0; p2 = tym+typ+2.d0
        r1yx = (   p1m*p2m * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)   &
                 - p1p*p2p * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(1)   &
	             + (p1p*p2p-p1m*p2m)                                &
			     * nod(n,g)%L(1)                                    &
		       ) / (p1m*p1p*p1)
        r2yx = (   p1m * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)  &
                 + p1p * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(1)  &
	             - p2  * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)  & 
		         ) / (p1m*p1p*p1)
        r1yz = (   p1m*p2m * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)   &
                 - p1p*p2p * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(3)   &
	             + (p1p*p2p-p1m*p2m)                                &
			     * nod(n,g)%L(3)                                    &
		       ) / (p1m*p1p*p1)
        r2yz = (   p1m * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)  &
                 + p1p * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(3)  &
	             - p2  * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)  & 
		       ) / (p1m*p1p*p1)
    END IF
ELSE
    tym = ydel(iy(n)-1)/ydel(iy(n))
    typ = ydel(iy(n)+1)/ydel(iy(n))
	p1m = tym+1.d0; p2m = 2.d0*tym+1.d0;p1p = typ+1.d0
	p2p = 2.d0*typ+1.d0; p1 = tym+typ+1.d0; p2 = tym+typ+2.d0
    r1yx = (   p1m*p2m * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(1)   &
             - p1p*p2p * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(1)   &
	         + (p1p*p2p-p1m*p2m)                                &
			 * nod(n,g)%L(1)                                    &
		   ) / (p1m*p1p*p1)
    r2yx = (   p1m * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(1)  &
             + p1p * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(1)  &
	         - p2  * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(1)  & 
		   ) / (p1m*p1p*p1)
    r1yz = (   p1m*p2m * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(3)   &
             - p1p*p2p * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(3)   &
	         + (p1p*p2p-p1m*p2m)                                &
			 * nod(n,g)%L(3)                                    &
		   ) / (p1m*p1p*p1)
    r2yz = (   p1m * nod(xyz(ix(n),iy(n)+1,iz(n)),g)%L(3)  &
             + p1p * nod(xyz(ix(n),iy(n)-1,iz(n)),g)%L(3)  &
	         - p2  * nod(xyz(ix(n),iy(n)  ,iz(n)),g)%L(3)  & 
		   ) / (p1m*p1p*p1)
END IF
	 
! Set paramaters for Z-Direction Transverse leakage
IF (iz(n) == 1 ) THEN
    IF (zbott == 0) THEN
        tzp = zdel(iz(n)+1)/zdel(iz(n))
	    p1p = tzp+1.d0
        r1zx = 2.d0 *                                 &
		     ( nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(1)   &
	         - nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)   &
			 ) / p1p
        r2zx = 0.d0
        r1zy = 2.d0 *                                 &
		     ( nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(2)   &
             - nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)   &
			 ) / p1p
        r2zy = 0.d0
	ELSE
        tzm = 1.d0
        tzp = zdel(iz(n)+1)/zdel(iz(n))
        p1m = tzm+1.d0; p2m = 2.d0*tzm+1.d0; p1p = tzp+1.d0
	    p2p = 2.d0*tzp+1.d0; p1 = tzm+tzp+1.d0; p2 = tzm+tzp+2.d0
        r1zx = (   p1m*p2m * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(1)   &
                 - p1p*p2p * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)   &
	             + (p1p*p2p-p1m*p2m)                                &
			     * nod(n,g)%L(1)                                    &
		       ) / (p1m*p1p*p1)
        r2zx = (   p1m * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(1)  &
                 + p1p * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)  &
	             - p2  * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)  & 
		       ) / (p1m*p1p*p1)
        r1zy = (   p1m*p2m * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(2)   &
                 - p1p*p2p * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)   &
	             + (p1p*p2p-p1m*p2m)                                &
			     * nod(n,g)%L(2)                                    &
		       ) / (p1m*p1p*p1)
        r2zy = (   p1m * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(2)  &
                 + p1p * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)  &
	             - p2  * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)  & 
		       ) / (p1m*p1p*p1)
    END IF
ELSE IF (iz(n) == nzz) THEN
    IF (ztop == 0) THEN
        tzm = zdel(iz(n)-1)/zdel(iz(n))
	    p1m = tzm+1.d0
        r1zx = 2.d0 *                                 &
		     ( nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)   &
	         - nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(1)   &
			 ) / p1m
        r2zx = 0.d0
        r1zy = 2.d0 *                                 &
		     ( nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)   &
             - nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(2)   &
			 ) / p1m
        r2zy = 0.d0
	ELSE
        tzm = zdel(iz(n)-1)/zdel(iz(n))
        tzp = 1.d0
        p1m = tzm+1.d0; p2m = 2.d0*tzm+1.d0; p1p = tzp+1.d0
	    p2p = 2.d0*tzp+1.d0; p1 = tzm+tzp+1.d0; p2 = tzm+tzp+2.d0
        r1zx = (   p1m*p2m * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)   &
                 - p1p*p2p * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(1)   &
	             + (p1p*p2p-p1m*p2m)                                &
			     * nod(n,g)%L(1)                                    &
		       ) / (p1m*p1p*p1)
        r2zx = (    p1m * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)  &
                  + p1p * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(1)  &
	              - p2  * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)  & 
		       ) / (p1m*p1p*p1)
        r1zy = (   p1m*p2m * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)   &
                 - p1p*p2p * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(2)   &
	             + (p1p*p2p-p1m*p2m)                                &
			     * nod(n,g)%L(2)                                    &
		       ) / (p1m*p1p*p1)
        r2zy = (   p1m * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)  &
                 + p1p * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(2)  &
	             - p2  * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)  & 
		       ) / (p1m*p1p*p1)	
	END IF
ELSE
    tzm = zdel(iz(n)-1)/zdel(iz(n))
    tzp = zdel(iz(n)+1)/zdel(iz(n))
    p1m = tzm+1.d0; p2m = 2.d0*tzm+1.d0; p1p = tzp+1.d0
	p2p = 2.d0*tzp+1.d0; p1 = tzm+tzp+1.d0; p2 = tzm+tzp+2.d0
    r1zx = (   p1m*p2m * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(1)   &
             - p1p*p2p * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(1)   &
	         + (p1p*p2p-p1m*p2m)                                &
			 * nod(n,g)%L(1)                                    &
		   ) / (p1m*p1p*p1)
    r2zx = (   p1m * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(1)  &
             + p1p * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(1)  &
	         - p2  * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(1)  & 
		   ) / (p1m*p1p*p1)
    r1zy = (   p1m*p2m * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(2)   &
             - p1p*p2p * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(2)   &
	         + (p1p*p2p-p1m*p2m)                                &
			 * nod(n,g)%L(2)                                    &
		   ) / (p1m*p1p*p1)
    r2zy = (   p1m * nod(xyz(ix(n),iy(n),iz(n)+1),g)%L(2)  &
             + p1p * nod(xyz(ix(n),iy(n),iz(n)-1),g)%L(2)  &
	         - p2  * nod(xyz(ix(n),iy(n),iz(n)  ),g)%L(2)  & 
		   ) / (p1m*p1p*p1)
END IF

! Set Transverse leakage Moments
L(1) = 0.d0
L(2) = ( r1xy/ydel(iy(n))+r1xz/zdel(iz(n)) ) / 12.d0
L(3) = ( r1yx/xdel(ix(n))+r1yz/zdel(iz(n)) ) / 12.d0
L(4) = ( r1zx/xdel(ix(n))+r1zy/ydel(iy(n)) ) / 12.d0
L(5) = ( r2xy/ydel(iy(n))+r2xz/zdel(iz(n)) ) / 20.d0
L(6) = ( r2yx/xdel(ix(n))+r2yz/zdel(iz(n)) ) / 20.d0
L(7) = ( r2zx/xdel(ix(n))+r2zy/ydel(iy(n)) ) / 20.d0

END SUBROUTINE TLUpd



SUBROUTINE FSrc(gt, s, sx1, sy1, sz1, sx2, sy2, sz2)
!
! Purpose:
!   To calculate fission source and fission source moments
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: mat, nnod, nuf, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2

IMPLICIT NONE

INTEGER, INTENT(IN) :: gt
REAL, DIMENSION(:), INTENT(INOUT) :: s, sx1, sy1, sz1
REAL, DIMENSION(:), INTENT(INOUT) :: sx2, sy2, sz2

INTEGER :: n

DO n = 1, nnod
	s(n)   = s(n)   + f0 (n,gt) * nuf(n,gt)
    sx1(n) = sx1(n) + fx1(n,gt) * nuf(n,gt)
	sy1(n) = sy1(n) + fy1(n,gt) * nuf(n,gt)
	sz1(n) = sz1(n) + fz1(n,gt) * nuf(n,gt)
	sx2(n) = sx2(n) + fx2(n,gt) * nuf(n,gt)
	sy2(n) = sy2(n) + fy2(n,gt) * nuf(n,gt)
	sz2(n) = sz2(n) + fz2(n,gt) * nuf(n,gt)
END DO

END SUBROUTINE FSrc


SUBROUTINE FSrcAd(gt, s, sx1, sy1, sz1, sx2, sy2, sz2)
!
! Purpose:
!   To calculate fission source and fission source moments
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: mat, nnod, chi, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2

IMPLICIT NONE

INTEGER, INTENT(IN) :: gt
REAL, DIMENSION(:), INTENT(INOUT) :: s, sx1, sy1, sz1
REAL, DIMENSION(:), INTENT(INOUT) :: sx2, sy2, sz2

INTEGER :: n

DO n = 1, nnod
	s(n)   = s(n)   + f0 (n,gt) * chi(n,gt)
    sx1(n) = sx1(n) + fx1(n,gt) * chi(n,gt)
	sy1(n) = sy1(n) + fy1(n,gt) * chi(n,gt)
	sz1(n) = sz1(n) + fz1(n,gt) * chi(n,gt)
	sx2(n) = sx2(n) + fx2(n,gt) * chi(n,gt)
	sy2(n) = sy2(n) + fy2(n,gt) * chi(n,gt)
	sz2(n) = sz2(n) + fz2(n,gt) * chi(n,gt)
END DO

END SUBROUTINE FSrcAd



SUBROUTINE SSrc(gt, s, sx1, sy1, sz1, sx2, sy2, sz2)
!
! Purpose:
!   To calculate scattering source and scattering source moments
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: ng, nnod, sigs, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2

IMPLICIT NONE

INTEGER, INTENT(IN) :: gt
REAL, DIMENSION(:), INTENT(INOUT) :: s, sx1, sy1, sz1
REAL, DIMENSION(:), INTENT(INOUT) :: sx2, sy2, sz2

INTEGER :: h, n

DO h = 1, ng
    DO n = 1, nnod
	    IF (gt /= h) THEN
			s(n)   = s(n)   + sigs(n,h,gt) * f0(n,h)
			sx1(n) = sx1(n) + sigs(n,h,gt) * fx1(n,h)
			sy1(n) = sy1(n) + sigs(n,h,gt) * fy1(n,h)
			sz1(n) = sz1(n) + sigs(n,h,gt) * fz1(n,h)
			sx2(n) = sx2(n) + sigs(n,h,gt) * fx2(n,h)
			sy2(n) = sy2(n) + sigs(n,h,gt) * fy2(n,h)
			sz2(n) = sz2(n) + sigs(n,h,gt) * fz2(n,h)
		END IF
	END DO
END DO

END SUBROUTINE SSrc



SUBROUTINE SSrcAd(gt, s, sx1, sy1, sz1, sx2, sy2, sz2)
!
! Purpose:
!   To calculate scattering source and scattering source moments
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: ng, nnod, sigs, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2

IMPLICIT NONE

INTEGER, INTENT(IN) :: gt
REAL, DIMENSION(:), INTENT(INOUT) :: s, sx1, sy1, sz1
REAL, DIMENSION(:), INTENT(INOUT) :: sx2, sy2, sz2

INTEGER :: h, n

DO h = 1, ng
    DO n = 1, nnod
	    IF (gt /= h) THEN
			s(n)   = s(n)   + sigs(n,gt,h) * f0(n,h)
			sx1(n) = sx1(n) + sigs(n,gt,h) * fx1(n,h)
			sy1(n) = sy1(n) + sigs(n,gt,h) * fy1(n,h)
			sz1(n) = sz1(n) + sigs(n,gt,h) * fz1(n,h)
			sx2(n) = sx2(n) + sigs(n,gt,h) * fx2(n,h)
			sy2(n) = sy2(n) + sigs(n,gt,h) * fy2(n,h)
			sz2(n) = sz2(n) + sigs(n,gt,h) * fz2(n,h)
		END IF
	END DO
END DO

END SUBROUTINE SSrcAd



SUBROUTINE TSrc(gt, Keff, sf0, sfx1, sfy1, sfz1, sfx2, sfy2, sfz2, &
                          s0,  sx1,  sy1,  sz1,  sx2,  sy2 , sz2   )
!
! Purpose:
!   To update total source
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nod, chi, nnod

IMPLICIT NONE

INTEGER, INTENT(IN) :: gt
REAL, INTENT(IN)    :: Keff
REAL, DIMENSION(:), INTENT(IN) :: sf0, sfx1, sfy1, sfz1
REAL, DIMENSION(:), INTENT(IN) :: sfx2, sfy2, sfz2
REAL, DIMENSION(:), INTENT(IN) :: s0, sx1, sy1, sz1
REAL, DIMENSION(:), INTENT(IN) :: sx2, sy2, sz2

INTEGER :: n

DO n = 1, nnod
	nod(n,gt)%Q(1) = chi(n,gt) * sf0(n)/Keff  + s0(n)
	nod(n,gt)%Q(2) = chi(n,gt) * sfx1(n)/Keff + sx1(n)
	nod(n,gt)%Q(3) = chi(n,gt) * sfy1(n)/Keff + sy1(n)
	nod(n,gt)%Q(4) = chi(n,gt) * sfz1(n)/Keff + sz1(n)
	nod(n,gt)%Q(5) = chi(n,gt) * sfx2(n)/Keff + sx2(n)
	nod(n,gt)%Q(6) = chi(n,gt) * sfy2(n)/Keff + sy2(n)
	nod(n,gt)%Q(7) = chi(n,gt) * sfz2(n)/Keff + sz2(n)
END DO 

END SUBROUTINE TSrc


SUBROUTINE TSrcFx(gt, Keff, sf0, sfx1, sfy1, sfz1, sfx2, sfy2, sfz2, &
                          s0,  sx1,  sy1,  sz1,  sx2,  sy2 , sz2   )
!
! Purpose:
!   To update total source
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nod, chi, nnod, exsrc

IMPLICIT NONE

INTEGER, INTENT(IN) :: gt
REAL, INTENT(IN)    :: Keff
REAL, DIMENSION(:), INTENT(IN) :: sf0, sfx1, sfy1, sfz1
REAL, DIMENSION(:), INTENT(IN) :: sfx2, sfy2, sfz2
REAL, DIMENSION(:), INTENT(IN) :: s0, sx1, sy1, sz1
REAL, DIMENSION(:), INTENT(IN) :: sx2, sy2, sz2

INTEGER :: n

DO n = 1, nnod
	
	nod(n,gt)%Q(1) = chi(n,gt) * sf0(n)  + s0(n)  + exsrc(n,gt)
	nod(n,gt)%Q(2) = chi(n,gt) * sfx1(n) + sx1(n)
	nod(n,gt)%Q(3) = chi(n,gt) * sfy1(n) + sy1(n)
	nod(n,gt)%Q(4) = chi(n,gt) * sfz1(n) + sz1(n)
	nod(n,gt)%Q(5) = chi(n,gt) * sfx2(n) + sx2(n)
	nod(n,gt)%Q(6) = chi(n,gt) * sfy2(n) + sy2(n)
	nod(n,gt)%Q(7) = chi(n,gt) * sfz2(n) + sz2(n)
	
END DO 

END SUBROUTINE TSrcFx



SUBROUTINE TSrcAd(gt, Keff, sf0, sfx1, sfy1, sfz1, sfx2, sfy2, sfz2, &
                          s0,  sx1,  sy1,  sz1,  sx2,  sy2 , sz2   )
!
! Purpose:
!   To update total source
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nod, nuf, nnod

IMPLICIT NONE

INTEGER, INTENT(IN) :: gt
REAL, INTENT(IN)    :: Keff
REAL, DIMENSION(:), INTENT(IN) :: sf0, sfx1, sfy1, sfz1
REAL, DIMENSION(:), INTENT(IN) :: sfx2, sfy2, sfz2
REAL, DIMENSION(:), INTENT(IN) :: s0, sx1, sy1, sz1
REAL, DIMENSION(:), INTENT(IN) :: sx2, sy2, sz2

INTEGER :: n

DO n = 1, nnod
	nod(n,gt)%Q(1) = nuf(n,gt) * sf0(n)/Keff  + s0(n)
	nod(n,gt)%Q(2) = nuf(n,gt) * sfx1(n)/Keff + sx1(n)
	nod(n,gt)%Q(3) = nuf(n,gt) * sfy1(n)/Keff + sy1(n)
	nod(n,gt)%Q(4) = nuf(n,gt) * sfz1(n)/Keff + sz1(n)
	nod(n,gt)%Q(5) = nuf(n,gt) * sfx2(n)/Keff + sx2(n)
	nod(n,gt)%Q(6) = nuf(n,gt) * sfy2(n)/Keff + sy2(n)
	nod(n,gt)%Q(7) = nuf(n,gt) * sfz2(n)/Keff + sz2(n)
END DO 

END SUBROUTINE TSrcAd


SUBROUTINE LxyzUpd (nt,gt)

USE sdata, ONLY: nod

! Purpose:
   ! To update Transverse leakages for group g and nod n

  ! Date                Programmer           History
 ! ========================================================
 ! 6 FEB 2017         Muhammad Imron       Original code

IMPLICIT NONE

INTEGER, INTENT(IN) :: gt, nt

nod(nt,gt)%L(1) = nod(nt,gt)%jo(1) - nod(nt,gt)%ji(1) &
                - nod(nt,gt)%ji(2) + nod(nt,gt)%jo(2)
nod(nt,gt)%L(2) = nod(nt,gt)%jo(3) - nod(nt,gt)%ji(3) &
                - nod(nt,gt)%ji(4) + nod(nt,gt)%jo(4)
nod(nt,gt)%L(3) = nod(nt,gt)%jo(5) - nod(nt,gt)%ji(5) &
                - nod(nt,gt)%ji(6) + nod(nt,gt)%jo(6)

END SUBROUTINE LxyzUpd


SUBROUTINE nodal_coup4()
!
! Purpose:
!    To calculate nodal coupling matrix
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: ng, nnod, xdel, ydel, zdel, ystag, &
                 ix, iy, iz, D, sigr, nod

IMPLICIT NONE

REAL, DIMENSION(6,6) :: A, B
REAL, DIMENSION(6,7) :: C

INTEGER :: n, g, i,j

REAL:: dx, dy, dz, lx, ly, lz
REAL :: ax, ay, az, ax1, ay1, az1, bx, by, bz
REAL :: bx1, by1, bz1, xy, xz, yx, yz, zx, zy

DO g= 1, ng 
    DO n = 1, nnod  
	
		dx = D(n,g) / xdel(ix(n))
		dy = D(n,g) / ydel(iy(n))
		dz = D(n,g) / zdel(iz(n))
		lx = 1.d0 / sigr(n,g) / xdel(ix(n))
		ly = 1.d0 / sigr(n,g) / ydel(iy(n))
		lz = 1.d0 / sigr(n,g) / zdel(iz(n))
		
        ax = 1.d0+32.d0*dx+120.d0*dx*lx+960.d0*dx*dx*lx+840.d0*dx*dx*lx*lx
		ay = 1.d0+32.d0*dy+120.d0*dy*ly+960.d0*dy*dy*ly+840.d0*dy*dy*ly*ly
		az = 1.d0+32.d0*dz+120.d0*dz*lz+960.d0*dz*dz*lz+840.d0*dz*dz*lz*lz
				
		ax1 = 8.d0*dx+60.d0*dx*lx+720.d0*dx*dx*lx+840.d0*dx*dx*lx*lx
		ay1 = 8.d0*dy+60.d0*dy*ly+720.d0*dy*dy*ly+840.d0*dy*dy*ly*ly
		az1 = 8.d0*dz+60.d0*dz*lz+720.d0*dz*dz*lz+840.d0*dz*dz*lz*lz
				
        bx = 1.d0-32.d0*dx+120.d0*dx*lx-960.d0*dx*dx*lx+840.d0*dx*dx*lx*lx
		by = 1.d0-32.d0*dy+120.d0*dy*ly-960.d0*dy*dy*ly+840.d0*dy*dy*ly*ly
		bz = 1.d0-32.d0*dz+120.d0*dz*lz-960.d0*dz*dz*lz+840.d0*dz*dz*lz*lz
				
		bx1 = -8.d0*dx+60.d0*dx*lx-720.d0*dx*dx*lx+840.d0*dx*dx*lx*lx
		by1 = -8.d0*dy+60.d0*dy*ly-720.d0*dy*dy*ly+840.d0*dy*dy*ly*ly
		bz1 = -8.d0*dz+60.d0*dz*lz-720.d0*dz*dz*lz+840.d0*dz*dz*lz*lz			
				
		xy = 20.d0*dx*ly+840.d0*dx*dx*lx*ly
		xz = 20.d0*dx*lz+840.d0*dx*dx*lx*lz
				
		yx = 20.d0*dy*lx+840.d0*dy*dy*ly*lx
		yz = 20.d0*dy*lz+840.d0*dy*dy*ly*lz
				
		zx = 20.d0*dz*lx+840.d0*dz*dz*lz*lx
		zy = 20.d0*dz*ly+840.d0*dz*dz*lz*ly

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
 
        C = 0.d0
				
		ax = 20.d0*dx*lx*xdel(ix(n))+840.d0*dx*dx*lx*lx*xdel(ix(n))
		ay = 20.d0*dy*ly*ydel(iy(n))+840.d0*dy*dy*ly*ly*ydel(iy(n))
	    az = 20.d0*dz*lz*zdel(iz(n))+840.d0*dz*dz*lz*lz*zdel(iz(n))
				
		C(1,1) = ax
		C(2,1) = ax
		C(3,1) = ay
		C(4,1) = ay
		C(5,1) = az
		C(6,1) = az
				
		ax1 = 60.d0*dx*lx*xdel(ix(n))
		ay1 = 60.d0*dy*ly*ydel(iy(n))
		az1 = 60.d0*dz*lz*zdel(iz(n))
				
		C(1,2) =  ax1
		C(2,2) = -ax1
		C(3,3) =  ay1
		C(4,3) = -ay1
		C(5,4) =  az1
		C(6,4) = -az1
				
		ax1 = 140.d0*dx*lx*xdel(ix(n))
		ay1 = 140.d0*dy*ly*ydel(iy(n))
		az1 = 140.d0*dz*lz*zdel(iz(n))
				
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



SUBROUTINE nodal_coup2()
!
! Purpose:
!    To calculate nodal coupling matrix
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: ng, nnod, xdel, ydel, zdel, ystag, &
                 ix, iy, iz, D, sigr, nod

IMPLICIT NONE

REAL, DIMENSION(6,6) :: A, B
REAL, DIMENSION(6,7) :: C

INTEGER :: n, g, i,j

REAL:: dx, dy, dz, lx, ly, lz
REAL :: ax, ay, az, ax1, ay1, az1, bx, by, bz
REAL :: bx1, by1, bz1, xy, xz, yx, yz, zx, zy


DO g= 1, ng 
    DO n = 1, nnod ! For every nodes from left to right that do not have zero materials
		dx = D(n,g) / xdel(ix(n))
		dy = D(n,g) / ydel(iy(n))
		dz = D(n,g) / zdel(iz(n))
		lx = 1.d0 / sigr(n,g) / xdel(ix(n))
		ly = 1.d0 / sigr(n,g) / ydel(iy(n))
		lz = 1.d0 / sigr(n,g) / zdel(iz(n))

        ax = 1.d0+8.d0*dx+6.d0*dx*lx 
		ay = 1.d0+8.d0*dy+6.d0*dy*ly 
		az = 1.d0+8.d0*dz+6.d0*dz*lz 
				
		ax1 = 4.d0*dx+6.d0*dx*lx
		ay1 = 4.d0*dy+6.d0*dy*ly
		az1 = 4.d0*dz+6.d0*dz*lz
				
        bx = 1.d0-8.d0*dx+6.d0*dx*lx 
		by = 1.d0-8.d0*dy+6.d0*dy*ly 
		bz = 1.d0-8.d0*dz+6.d0*dz*lz 
				
		bx1 = -4.d0*dx+6.d0*dx*lx
		by1 = -4.d0*dy+6.d0*dy*ly
		bz1 = -4.d0*dz+6.d0*dz*lz			
				
		xy = 6.d0*dx*ly
		xz = 6.d0*dx*lz
				
		yx = 6.d0*dy*lx
		yz = 6.d0*dy*lz
				
		zx = 6.d0*dz*lx
		zy = 6.d0*dz*ly

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
 
        C      = 0.d0
		C(1,1) = 6.d0*dx*lx*xdel(ix(n))
		C(2,1) = 6.d0*dx*lx*xdel(ix(n))
		C(3,1) = 6.d0*dy*ly*ydel(iy(n))
		C(4,1) = 6.d0*dy*ly*ydel(iy(n))
		C(5,1) = 6.d0*dz*lz*zdel(iz(n))
		C(6,1) = 6.d0*dz*lz*zdel(iz(n))
				
		CALL inverse(g, n, A) 
				
		nod(n,g)%P = MATMUL(A, B)
		nod(n,g)%R = MATMUL(A, C)
    END DO
END DO 


END SUBROUTINE nodal_coup2



SUBROUTINE inverse (gt, nt, mat)
!
! Purpose:
!    To perform matrix inverse by LU decomposition
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE InpOutp, ONLY: ounit
USE sdata,   ONLY: ix, iy, iz

IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(INOUT) :: mat
INTEGER, INTENT(IN) :: gt, nt

REAL, DIMENSION(6,6) :: L, U, imat, pmat
REAL, DIMENSION(6) :: y
REAL :: piv, isum
INTEGER :: i, j, k

pmat = mat
U = mat
L = 0.d0

! Start matrix decomposition
DO i= 1, 6
    IF (ABS(mat(i,i)) < 10e-3) THEN
      WRITE(ounit,*) 'ERROR IN MATRIX DECOMP: DIAGONAL ELEMENTS CLOSE TO ZERO'
	  WRITE(ounit,2001) gt, ix(nt), iy(nt), iz(nt)
	  STOP
    END IF
	L(i,i) = 1.d0
    DO j= i+1, 6
	    piv = U(j,i)/U(i,i)
		L(j,i) = piv
        DO k= i, 6
            U(j,k) = U(j,k) - piv*U(i,k)		
		END DO  
		U(j,i) = 0.d0
    END DO
END DO

DO i = 1,6
    DO j = 1,6
	    isum = 0.d0
		DO k = 1,6
		    isum = isum+L(i,k)*U(k,j)
		END DO
		IF (ABS(mat(i,j)-isum)/ABS(mat(i,j)) > 1.e-5) THEN
            WRITE(ounit,*) 'ERROR IN MATRIX DECOMP: DECOMPOSITION FAILED'
			WRITE(ounit,2001) gt, ix(nt), iy(nt), iz(nt)
	        STOP		
		END IF
    END DO
END DO	

!Initialiaze Identity matrix
imat = 0.d0
DO i= 1, 6
    imat(i,i) = 1.d0
END DO

! Calculate matrix inverse
DO j=1,6
    !Solve y in Ly = b (Forward substitution)
	y(1) = imat(1,j)
	DO i=2,6
	    isum = 0.d0
	    DO k =1, i-1
	        isum = isum + L(i,k)*y(k)
		END DO
		y(i) = imat(i,j)-isum
	END DO
	
	! Solve x in Ux=y(Backward substitution) and store inverse matrix to input matrix 'mat'
    mat(6,j) = y(6)/U(6,6)
	DO i = 5,1,-1
	    isum = 0.d0
	    DO k =i+1,6
	        isum = isum + U(i,k)*mat(k,j)
		END DO	 
        mat(i,j) = (y(i)-isum) / U(i,i)
	END DO
END DO

!Check matrix Inverse
DO i = 1,6
    DO j = 1,6
	    isum = 0.d0
		DO k = 1,6
		    isum = isum+pmat(i,k)*mat(k,j)
		END DO
		IF (ABS(imat(i,j)-isum) > 1.e-4) THEN
            WRITE(ounit,*) 'ERROR IN MATRIX INVERSION'
			WRITE(ounit,2001) gt, ix(nt), iy(nt), iz(nt)
	        STOP		
		END IF
    END DO
END DO

2001 FORMAT(2X, 'Group = ', I3, ', I = ', I5, ', J = ', I5, ', K = ', I5)


END SUBROUTINE inverse


REAL FUNCTION Integrate(s)

USE sdata, ONLY: nnod, xdel, ydel, zdel, ix, iy, iz

IMPLICIT NONE

REAL, DIMENSION (:), INTENT(IN) :: s
INTEGER :: n

Integrate = 0.
DO n = 1, nnod
    Integrate = Integrate + xdel(ix(n))*ydel(iy(n))*zdel(iz(n)) * s(n)
END DO	 

END FUNCTION Integrate



SUBROUTINE matvec (mat, vec, rvec)
!
! Purpose:
!    To perform matrix vector multiplication
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(IN) :: mat
REAL, DIMENSION(:), INTENT(IN) :: vec
REAL, DIMENSION(:), INTENT(OUT) :: rvec

INTEGER :: i, j, n, m
REAL :: isum

m = SIZE(mat,1)
n = SIZE(vec,1)

DO i= 1, m
    isum = 0.
    DO j = 1, n
        isum = isum + mat(i,j)*vec(j)
	END DO
    rvec(i)	 = isum
END DO		

END SUBROUTINE matvec



SUBROUTINE RelE(newF, oldF, rel)

USE sdata, ONLY: nnod

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: newF, oldF
REAL, INTENT(OUT) :: rel

REAL :: error
INTEGER :: n

rel = 0.

DO n= 1, nnod
    IF (newF(n) /= 0.d0) THEN
	    error = ABS(newF(n) - oldF(n)) / ABS(newF(n))
        IF (error > rel) rel = error
	END IF
END DO

END SUBROUTINE RelE




SUBROUTINE CNeg(f, neg)

USE sdata, ONLY: ng, nnod

IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(IN) :: f
INTEGER, INTENT(OUT) :: neg

INTEGER :: g, n

neg = 0


DO g = 1, ng
    DO n = 1, nnod
            IF (f(n,g) < 0.d0) neg = neg + 1
    END DO
END DO

END SUBROUTINE CNeg



SUBROUTINE MultF(k)

! Purpose: To calculate Keff from fixed source problem

USE sdata, ONLY: ng, nnod, nod, nxx, nyy, nzz, &
                 ix, iy, iz, xyz, nuf, siga, f0, &
				 xstag, ystag, xdel, ydel, zdel

IMPLICIT NONE

REAL, INTENT(OUT) :: k

REAL :: leak, absp, fiss

INTEGER :: g, n

!! Jot Nodals' outgoing currents  (X+, X-, Y+, Y-, Z+, Z-)
!! Jin Nodals' ingoing currents   (X+, X-, Y+, Y-, Z+, Z-)

leak = 0.d0
absp = 0.d0
fiss = 0.d0
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

k = fiss / (leak + absp)

END SUBROUTINE MultF



SUBROUTINE Init()

!
! Purpose:
!    To calculate power distribution
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code


USE sdata, ONLY: ng, nod, nnod, Ke, &
                 f0, fx1, fy1, fz1, fx2, fy2, fz2
USE InpOutp, ONLY: ounit

IMPLICIT NONE

INTEGER :: g, n, istat

Ke = 1.d0

ALLOCATE(nod(nnod,ng))
DO g = 1, ng
    DO n = 1, nnod
        ALLOCATE(nod(n,g)%P(6,6))
        ALLOCATE(nod(n,g)%R(6,7))
    END DO
END DO

CALL nodal_coup4()

ALLOCATE(f0(nnod,ng), &
         fx1(nnod,ng), fy1(nnod,ng), fz1(nnod,ng), &
		 fx2(nnod,ng), fy2(nnod,ng), fz2(nnod,ng), &
		 STAT=istat) 
IF (istat /= 0) THEN 
    WRITE(ounit,*) '[1] NOT ENOUGH MEMORY. PROGRAM STOP'
	STOP
END IF

DO g= 1, ng
    DO n = 1, nnod
		   ALLOCATE(nod(n,g)%Q(7))
		   ALLOCATE(nod(n,g)%L(3))
           ALLOCATE(nod(n,g)%jo(6), nod(n,g)%ji(6))
		   nod(n,g)%jo = 1.d0
		   nod(n,g)%ji = 1.d0
		   nod(n,g)%Q  = 0.d0
		   
		   CALL LxyzUpd(n,g)
		   
		   f0(n,g)  = 1.d0
		   fx1(n,g) = 1.d0
		   fy1(n,g) = 1.d0
		   fz1(n,g) = 1.d0
		   fx2(n,g) = 1.d0
		   fy2(n,g) = 1.d0
		   fz2(n,g) = 1.d0
    END DO
END DO


END SUBROUTINE Init



SUBROUTINE PowDis (p)

!
! Purpose:
!    To calculate power distribution
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code


USE sdata, ONLY: ng, nnod, sigf, f0

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(OUT) :: p
INTEGER :: g, n

p = 0.d0
DO g= 1, ng
    DO n= 1, nnod
	    p(n) = p(n) + f0(n,g) * sigf(n,g)
	END DO					 
END DO	 


END SUBROUTINE PowDis


SUBROUTINE forward()

!
! Purpose:
!    To calculate power distribution
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nnod, f0
USE InpOutp, ONLY: ounit, AsmPow, AxiPow, AsmFlux

IMPLICIT NONE

REAL, DIMENSION(nnod) :: pow


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

CALL outer4()

CALL PowDis(pow)

CALL AsmPow(pow)

CALL AxiPow(pow)

CALL AsmFlux(f0, 1.e0)	
 


END SUBROUTINE forward



SUBROUTINE adjoint()

!
! Purpose:
!    To calculate power distribution
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: f0
USE InpOutp, ONLY: ounit, AsmFlux

IMPLICIT NONE

INTEGER :: g


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

CALL outer4Ad()

CALL AsmFlux(f0, 1.e0)	
 


END SUBROUTINE adjoint



SUBROUTINE fixedsrc()

!
! Purpose:
!    To calculate power distribution
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: f0
USE InpOutp, ONLY: ounit, AsmFlux

IMPLICIT NONE

INTEGER :: g


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

CALL outer4Fx()

CALL AsmFlux(f0)	
 


END SUBROUTINE fixedsrc


END MODULE nodal
