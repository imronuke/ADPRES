MODULE InpOutp

IMPLICIT NONE

SAVE

CHARACTER(LEN=1) :: ind          ! used to read x indicator in input buffer to prevent error
CHARACTER(LEN=100) :: message    ! error message
!Ouput options
LOGICAL, PARAMETER :: ogeom = .TRUE.  ! Geometry output option
LOGICAL, PARAMETER :: oxsec = .TRUE.  ! Macroscopic CXs output option

! Inut, output and buffer input file unit number
INTEGER, PARAMETER :: iunit = 100   !input
INTEGER, PARAMETER :: ounit = 101   !output
INTEGER, PARAMETER :: umode = 111, uxsec = 112, ugeom = 113
INTEGER, PARAMETER :: ucase = 114, uesrc = 115
INTEGER :: bunit

INTEGER :: bmode, bxsec, bgeom, bcase, besrc

CHARACTER(LEN=100):: iline

CONTAINS

SUBROUTINE inp_read()
!
! Purpose:
!    [Main subroutine in this modul] To read input, echo the 
!    input and gives the description
!    to the user about his/her input
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: ng, nnod, chi, nuf, sigf, siga, &
                 D, sigr, sigs, mat, mode

IMPLICIT NONE

INTEGER :: iost
CHARACTER(LEN=20) :: iname, oname
WRITE(*,'(A,A100)',ADVANCE='NO') '  INPUT NAME : '
READ(*,*) iname

iname = TRIM(iname)

OPEN (UNIT=iunit, FILE=iname, STATUS='OLD', ACTION='READ', &
      IOSTAT = iost)
	  


IF (iost /= 0) THEN  
    WRITE(*,1020) iost
	WRITE(*,*) '  NO FILE : ', iname
    1020 FORMAT	(2X, 'File Open Failed--status', I6)
	STOP
END IF	

oname = TRIM(iname) // '.out'
oname = TRIM(oname)

OPEN (UNIT=ounit, FILE=oname, STATUS='REPLACE', ACTION='WRITE')	

OPEN (UNIT=umode, STATUS='SCRATCH', ACTION='READWRITE')	
OPEN (UNIT=uxsec, STATUS='SCRATCH', ACTION='READWRITE')	
OPEN (UNIT=ugeom, STATUS='SCRATCH', ACTION='READWRITE')
OPEN (UNIT=ucase, STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=uesrc, STATUS='SCRATCH', ACTION='READWRITE') 

! OPEN (UNIT=umode, FILE = 'umode', STATUS='REPLACE', ACTION='READWRITE')	
! OPEN (UNIT=uxsec, FILE = 'uxsec', STATUS='REPLACE', ACTION='READWRITE')	
! OPEN (UNIT=ugeom, FILE = 'ugeom', STATUS='REPLACE', ACTION='READWRITE')
! OPEN (UNIT=ucase, FILE = 'ucase', STATUS='REPLACE', ACTION='READWRITE') 
! OPEN (UNIT=uesrc, FILE = 'uesrc', STATUS='REPLACE', ACTION='READWRITE') 

   

CALL inp_echo()
CALL inp_comments ()

REWIND(umode)
REWIND(uxsec)
REWIND(ugeom)
REWIND(ucase)
REWIND(uesrc)

! Start reading buffer files for each card

WRITE(ounit,*) 
WRITE(ounit,*) 
WRITE(ounit,1008) 
WRITE(ounit,*) &
' ***************************************************************************************************'

! Card MODE
IF (bmode > 0) THEN
    CALL inp_mode(umode)
ELSE
    WRITE(ounit,1021) '%MODE'
	STOP
END IF

! Card CASE
IF (bcase == 1) CALL inp_case (ucase)

! Card XSEC
IF (bxsec == 1) THEN
    CALL inp_xsec(uxsec)
ELSE
    WRITE(ounit,1021) '%XSEC'
	STOP
END IF	

! Card GEOM
IF (bgeom == 1) THEN
    CALL inp_geom1(ugeom)
    CALL inp_geom2(ugeom)
ELSE
    WRITE(ounit,1021) '%GEOM'
	STOP
END IF

! Card ESRC
IF (mode == 'FIXEDSRC' .AND. besrc == 1) THEN
    CALL inp_esrc(uesrc)
ELSE IF (mode == 'FIXEDSRC' .AND. besrc /= 1) THEN
    WRITE(ounit,*) '  CALCULATION MODE IS FIXED SOURCE'
    WRITE(ounit,1021) 'ESRC'
	STOP
ELSE IF (mode /= 'FIXEDSRC' .AND. besrc == 1) THEN
    WRITE(ounit,*) '  ESRC CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
	STOP
ELSE
    CONTINUE
END IF	

WRITE(ounit,*) 
WRITE(ounit,*) 
WRITE(ounit,*) &
' **********************************', '  STOP READING INPUT  ', '*******************************************'	


! CXs preparation
ALLOCATE(D(nnod,ng), chi(nnod,ng), nuf(nnod,ng), siga(nnod,ng), &
         sigf(nnod,ng), sigr(nnod,ng), sigs(nnod,ng,ng))
		 
CALL node_cx(mat)	

1008 FORMAT (45X, 'START READING INPUT')
1021 FORMAT(2X, 'CARD', A6, ' DOES NOT PRESENT. THIS CARD IS MANDATORY')	 

END SUBROUTINE inp_read


SUBROUTINE inp_echo()
!
! Purpose:
!    To rewrite the input
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

IMPLICIT NONE

INTEGER :: eof
INTEGER :: nline

WRITE(ounit, *) '                      ###########################################################'
WRITE(ounit, *) '                      #                  ADPRES v.0.1.ALPHA                     #'
WRITE(ounit, *) '                      #        ABU DHABI POLYTECHNIC REACTOR SIMULATOR          #'
WRITE(ounit, *) '                      ###########################################################  '
WRITE(ounit, *)
WRITE(ounit, *)


WRITE(ounit,1002) 'STARTS'

nline = 0
DO
   READ(iunit, '(A100)', IOSTAT=eof) iline
   nline = nline + 1
   IF (eof < 0) THEN
       WRITE(ounit,1002) 'ENDS'
	   WRITE(ounit,*)       
	   EXIT
   END IF
   WRITE(ounit, 1001) nline, iline
END DO

1001 FORMAT (2X, I4, ': ', A100)
1002 FORMAT	(2X, '========================================INPUT DATA',A7,' HERE=====================================')
REWIND (iunit)

END SUBROUTINE inp_echo


SUBROUTINE inp_comments ()
!
! Purpose:
!    To remove the comments in input and rewrite the 
!    input into input buffer. Comments marked by !.
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

IMPLICIT NONE

INTEGER :: ln                  ! line number
INTEGER :: eof, comm, per
INTEGER :: buffer

bmode = 0; bxsec = 0; bgeom = 0
bcase = 0; besrc = 0

buffer = 99
OPEN (UNIT=buffer, STATUS='SCRATCH', ACTION='READWRITE')  

ln = 0   
DO
    ln = ln+1
    READ (iunit, '(A100)', IOSTAT=eof) iline
    IF (eof < 0) EXIT              !Check end of file	
	iline = TRIM(ADJUSTL(iline))   ! Remove trailing blanks 
    comm = INDEX(iline, '!')       ! Find position '!' if any
    ! If there is no '!' and no first ten blank spaces
    IF (comm == 0 .AND. iline(1:10) /= '          ')  THEN
	    WRITE(buffer,1012)'x ',ln,iline
    END IF
    !If the first character is not '!'
    IF (comm > 1) THEN 
	    iline = iline(1:comm-1)       ! Take only part of input
	    WRITE(buffer,1012)'x ',ln, iline
    END IF
END DO  

REWIND(buffer)

DO
    per = 0
    READ (buffer, '(A2,I5,A100)', IOSTAT=eof) ind, ln, iline
    IF (eof < 0) EXIT              !Check end of file
    per = INDEX(iline,'%')
    IF (per > 0) THEN              ! IF %card detected
        iline = iline(per+1:100)
		iline = TRIM(ADJUSTL(iline))
		SELECT CASE (iline)
		    CASE('MODE')
			    bunit = umode
				bmode = 1
		    CASE('XSEC')
			    bunit = uxsec	
				bxsec = 1
		    CASE('GEOM')
			    bunit = ugeom
				bgeom = 1
		    CASE('CASE')
			    bunit = ucase
				bcase = 1
		    CASE('ESRC')
			    bunit = uesrc
				besrc = 1
		    CASE DEFAULT
			    WRITE(ounit,1014) ln, iline
				STOP
		END SELECT
    END IF		 
	IF (per == 0) WRITE(bunit, 1012) 'x ',ln, iline
END DO  

1012 FORMAT(A2, I5,' ',A100)
1014 FORMAT(2X, 'AT LINE', I3, ' : WRONG CARD ', A8)


 
END SUBROUTINE inp_comments


SUBROUTINE inp_mode (xbunit)
!
! Purpose:
!    To read case mode in input
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: mode

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

READ(xbunit, '(A2, I5,A100)', IOSTAT=ios) ind, ln, mode
message = ' error in reading MODE'
CALL er_message(ounit, ios, ln, message)


mode = TRIM(ADJUSTL(mode))  ! ADJUSTL = MOVE PRECEDING BLANK TO TRAILING

SELECT CASE(mode)
    CASE('FORWARD') 
	    WRITE(ounit,1031) 'FORWARD CALCULATION'
		WRITE(ounit,*)
    CASE('ADJOINT') 
	    WRITE(ounit,1031) 'ADJOINT CALCULATION'
		WRITE(ounit,*)
    CASE('FIXEDSRC') 
	    WRITE(ounit,1031) 'FIXED SOURCE CALCULATION'
		WRITE(ounit,*)
	CASE DEFAULT
	    WRITE(ounit,1032) mode
		STOP
END SELECT		

1031 FORMAT(2X, 'CALCULATION MODE : ', A30)	
1032 FORMAT(2X, 'MODE : ', A10, ' UNIDENTIFIED')	

END SUBROUTINE inp_mode


SUBROUTINE inp_case (xbunit)
!
! Purpose:
!    To read case card in input
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

CHARACTER(LEN=100) :: case_id
CHARACTER(LEN=100) :: case_exp

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

READ(xbunit, '(A2, I5,A100)', IOSTAT=ios) ind, ln, case_id
message = ' error in CASE ID'
CALL er_message(ounit, ios, ln, message)
READ(xbunit, '(A2, I5,A100)', IOSTAT=ios) ind, ln, case_exp
message = ' error in case description'
CALL er_message(ounit, ios, ln, message)

WRITE(ounit,1006) case_id
WRITE(ounit,1007) case_exp
1006 FORMAT(2X, 'CASE ID : ', A100)
1007 FORMAT(1X, A100)

END SUBROUTINE inp_case


SUBROUTINE inp_xsec (xbunit)
!
! Purpose:
!    To read CROSS SECTIONS card in input
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nmat, mat, ng

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: i, g, h
INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status
REAL :: dum
INTEGER, DIMENSION(:), ALLOCATABLE :: group

WRITE(ounit,*) 
WRITE(ounit,*) 
WRITE(ounit,*) '           >>>> READING MACROSCOPIC CROSS SECTIONS <<<<'
WRITE(ounit,*) '           --------------------------------------------'

READ(xbunit, *, IOSTAT=ios) ind, ln, ng, nmat
message = ' error in material number'
CALL er_message(ounit, ios, ln, message)


ALLOCATE(group(ng))
DO g = 1,ng
   group(g) = g
END DO

ALLOCATE(mat(nmat))

! Reading MACROSCOPIC CXs
DO i= 1, nmat
	ALLOCATE(mat(i)%sigtr(ng), mat(i)%D(ng), mat(i)%siga(ng), &
    mat(i)%sigr(ng), mat(i)%nuf(ng), mat(i)%sigf(ng), &
	mat(i)%chi(ng), mat(i)%sigs(ng,ng))
    DO g= 1, ng
        READ(xbunit, *, IOSTAT=ios) ind, ln, mat(i)%sigtr(g), &
		mat(i)%siga(g), mat(i)%nuf(g), mat(i)%sigf(g), &
		mat(i)%chi(g), (mat(i)%sigs(g,h), h = 1, ng)
		message = ' error in cross section data'
		CALL er_message(ounit, ios, ln, message)
		mat(i)%D(g) = 1.d0/(3.d0*mat(i)%sigtr(g)) 
		dum = 0.0
		DO h= 1, ng
		    IF (g /= h) dum = dum + mat(i)%sigs(g,h)
		END DO
		mat(i)%sigr(g) =  mat(i)%siga(g)+dum
	END DO
	
	! WRITE(ounit,*) '! MATERIAL ', i
	! DO g = 1, ng
	    ! WRITE(ounit,'(10F6.3)') mat(i)%sigtr(g), mat(i)%siga(g), mat(i)%nuf(g), mat(i)%sigf(g), &
		! mat(i)%chi(g), (mat(i)%sigs(g,h), h = 1, ng)
	! END DO	
END DO	

! Writing output
IF (oxsec) THEN
    DO i= 1, nmat
        WRITE(ounit,*) 
        WRITE(ounit,1009) i
	    WRITE(ounit,1011)'GROUP', 'TRANSPORT', 'DIFFUSION', 'ABSORPTION', &
	    'REMOVAL', 'NU*FISS', 'FISSION','FISS. SPCTR'
        DO g= 1, ng
		    WRITE(ounit,1010) g, mat(i)%sigtr(g), mat(i)%D(g), mat(i)%siga(g), &
		    mat(i)%sigr(g), mat(i)%nuf(g), mat(i)%sigf(g), mat(i)%chi(g)
	    END DO
	    WRITE(ounit,*)'  --SCATTERING MATRIX--'
	    WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
	    DO g= 1, ng
	        WRITE(ounit,1015)g, (mat(i)%sigs(g,h), h=1,ng)
	    END DO
    END DO
END IF	

WRITE(ounit,*)
WRITE(ounit,*) ' ...Macroscopic CXs are sucessfully read...'

1009 FORMAT(5X, 'MATERIAL', I3)
1011 FORMAT(2X, A7, A12, A13, A12, A11, 2A13, A15)
1010 FORMAT(2X, I6, F13.6, 3F12.6, 3F13.6)
1015 FORMAT(4X, I3, F16.6, 20F12.6)

DEALLOCATE(group)

END SUBROUTINE inp_xsec


SUBROUTINE inp_geom1 (xbunit)
!
! Purpose:
!    To read geometry card in input (1st part)
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nx, ny, nz, nxx, nyy, nzz, xdel, ydel, zdel, &
                xdiv, ydiv, zdiv

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln, ios

INTEGER :: i, j, k, lx, ly, lz, xtot, ytot, ztot
REAL, DIMENSION(:), ALLOCATABLE :: xsize, ysize, zsize  !Assembly size
REAL :: div

WRITE(ounit,*) 
WRITE(ounit,*) 
WRITE(ounit,*) '           >>>>>READING CORE GEOMETRY<<<<<'
WRITE(ounit,*) '           -------------------------------'

! Read number of assemblies in x, y and z directions
READ(xbunit, *, IOSTAT=ios) ind, ln, nx, ny, nz
message = ' error in reading number assemblies'
CALL er_message(ounit, ios, ln, message)

IF (nx < 2) THEN
    WRITE(ounit,*) '  Error: nx shall be at least 2'
	STOP
END IF	
IF (ny < 2) THEN
    WRITE(ounit,*) '  Error: ny shall be at least 2'
	STOP
END IF	
IF (nz < 2) THEN
    WRITE(ounit,*) '  Error: nz shall be at least 2'
	STOP
END IF	

ALLOCATE(xsize(nx), ysize(ny), zsize(nz))
ALLOCATE(xdiv(nx), ydiv(ny), zdiv(nz))

! Read assembly sizes and assembly division in x, y and z directions
! x-direction
READ(xbunit, *, IOSTAT=ios) ind, ln, (xsize(i), i=1,nx)
message = ' error in reading x-direction assembly sizes'
CALL er_message(ounit, ios, ln, message)
READ(xbunit, *, IOSTAT=ios) ind, ln, (xdiv(i), i=1,nx)
message = ' error in reading x-direction assembly division'
CALL er_message(ounit, ios, ln, message)
! y-direction
READ(xbunit, *, IOSTAT=ios) ind, ln, (ysize(j), j=1,ny)
message = ' error in reading y-direction assembly sizes'
CALL er_message(ounit, ios, ln, message)
READ(xbunit, *, IOSTAT=ios) ind, ln, (ydiv(j), j=1,ny)
message = ' error in reading y-direction assembly division'
CALL er_message(ounit, ios, ln, message)
! z-direction
READ(xbunit, *, IOSTAT=ios) ind, ln, (zsize(k), k=1,nz)
message = ' error in reading z-direction assembly sizes'
CALL er_message(ounit, ios, ln, message)
READ(xbunit, *, IOSTAT=ios) ind, ln, (zdiv(k), k=1,nz)
message = ' error in reading z-direction assembly division'
CALL er_message(ounit, ios, ln, message)

!Calculate number of nodes in x, y and z directions
!x-direction
nxx=0
DO i= 1,nx
nxx = nxx+xdiv(i)
END DO
!y-direction
nyy=0
DO j= 1,ny
nyy = nyy+ydiv(j)
END DO
!z-direction
nzz=0
DO k= 1,nz
nzz = nzz+zdiv(k)
END DO

ALLOCATE(xdel(nxx), ydel(nyy), zdel(nzz))

!Calculate delta x, y, and z (node sizes)
!Delta x
xtot=0
DO i= 1,nx
    div = xsize(i)/REAL(xdiv(i))
    DO lx= 1, xdiv(i)
	xtot = xtot+1
    xdel(xtot) = div
    END DO	
END DO
!Delta y
ytot=0
DO j= 1,ny
    div = ysize(j)/REAL(ydiv(j))
    DO ly= 1, ydiv(j)
	ytot = ytot+1
    ydel(ytot) = div
    END DO	
END DO
!Delta z
ztot=0
DO k= 1,nz
    div = zsize(k)/REAL(zdiv(k))
    DO lz= 1, zdiv(k)
	ztot = ztot+1
    zdel(ztot) = div
    END DO	
END DO

DEALLOCATE(xsize, ysize, zsize)

END SUBROUTINE inp_geom1




SUBROUTINE inp_geom2 (xbunit)
!
! Purpose:
!    To read geometry card in input (2nd part)
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nx, ny, nz, nxx, nyy, nzz, xdel, ydel, zdel, &
                xleft, xrigt, yback, yfrnt, zbott, ztop, mnum, &
				xstag, ystag, xdiv, ydiv, zdiv, ix, iy, iz, xyz, &
				nnod, nmat

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln, ios

INTEGER :: i, j, k, lx, ly, lz, xtot, ytot, ztot, n
INTEGER :: np                                           ! Numbe rof planars
INTEGER, DIMENSION(:), ALLOCATABLE :: zpln              ! Planar assignment to z direction
TYPE :: MAT_ASGN                                        ! Material assignment
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: asm         ! Material assignment into assembly
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: node        ! Material assignment into nodes
END TYPE
TYPE(MAT_ASGN), DIMENSION(:), ALLOCATABLE :: plnr       ! planar
CHARACTER(LEN=2), DIMENSION(nxx, nyy) :: mmap

! Reading number of planar
READ(xbunit, *, IOSTAT=ios) ind, ln, np
message = ' error in reading number of planars'
CALL er_message(ounit, ios, ln, message)

!Reading planar assignment into z-direction
ALLOCATE(zpln(nz))
READ(xbunit, *, IOSTAT=ios) ind, ln, (zpln(k), k=1,nz)
message = ' error in reading number of planars'
CALL er_message(ounit, ios, ln, message)

! Reading material assignment for each planar
ALLOCATE(plnr(np))
DO k= 1,np
    ALLOCATE(plnr(k)%asm (nx,ny))
	ALLOCATE(plnr(k)%node(nxx,nyy))
	DO j= 1, ny
	    READ(xbunit, *, IOSTAT=ios) ind, ln, (plnr(k)%asm(i,j), i=1,nx)
        message = ' error in reading planar'
        CALL er_message(ounit, ios, ln, message)
	END DO
END DO

! Material assignment into nodes (not part for geom output)
ALLOCATE (mnum(nxx, nyy, nzz))

ztot = 0
DO k= 1, nz
    DO lz= 1, zdiv(k)
        ztot = ztot+1
		ytot = 0
	    DO j= 1, ny
	        DO ly= 1, ydiv(j)
	            ytot = ytot+1
				xtot = 0
		        DO i= 1, nx
			        DO lx= 1, xdiv(i)
				        xtot = xtot+1
				        mnum(xtot, ytot, ztot)        = plnr(zpln(k))%asm(i,j)
						IF (mnum(xtot, ytot, ztot) > nmat) THEN
						    WRITE(ounit,'(2X,A17,I3,A37)') 'ERROR: MATERIAL ', &
							mnum(xtot, ytot, ztot), ' IS GREATER THAN NUMBER OF MATERIAL'
						    STOP
						END IF
						plnr(zpln(k))%node(xtot,ytot) = plnr(zpln(k))%asm(i,j)
					END DO	
				END DO
			END DO
		END DO
	END DO			
END DO

! -Indexing non zero material for staggered mesh-
ALLOCATE(ystag(nyy), xstag(nxx))
!Indexing non zero material for staggered mesh along y direction
DO j= 1, nyy
    ystag(j)%smin = nxx
    DO i = 1, nxx
	    IF (mnum(i,j,1) /= 0) THEN
		    ystag(j)%smin = i
			EXIT
		END IF
	END DO
END DO
DO j= 1, nyy
    ystag(j)%smax = 0
    DO i = nxx, 1, -1
	    IF (mnum(i,j,1) /= 0) THEN
		    ystag(j)%smax = i
			EXIT
		END IF
	END DO
END DO	
!Indexing non zero material for staggered mesh along x direction
DO i= 1, nxx
    xstag(i)%smin = nyy
    DO j = 1, nyy
	    IF (mnum(i,j,1) /= 0) THEN
		    xstag(i)%smin = j
			EXIT
		END IF
	END DO
END DO

DO i= 1, nxx
    xstag(i)%smax = 0
    DO j = nyy, 1, -1
	    IF (mnum(i,j,1) /= 0) THEN
		    xstag(i)%smax = j
			EXIT
		END IF
	END DO
END DO

! Checking zero material between non-zero material
DO k = 1, nzz	
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax	
			IF (mnum(i,j,k) == 0) THEN
				WRITE(ounit,*) 'Zero material found inside core. Check material assignment'
				STOP
		    END IF
        END DO
	END DO
END DO
DO k = 1, nzz	
    DO i = 1, nxx
        DO j = xstag(i)%smin, xstag(i)%smax	
			IF (mnum(i,j,k) == 0) THEN
				WRITE(ounit,*) 'Zero material found inside core. Check material assignment'
				STOP
		    END IF
        END DO
	END DO
END DO

!Reading Boundary Conditions
READ(xbunit, *, IOSTAT=ios) ind, ln, xleft, xrigt, yback, yfrnt, zbott, ztop
message = ' error in reading boundary conditions'
CALL er_message(ounit, ios, ln, message)
	

! Wrting core geometry output
IF (ogeom) THEN
    WRITE(ounit,*)' Number of assembly in x, y and z directions respectively :'
    WRITE(ounit,*) nx, ny, nz
    WRITE(ounit,*)' Number of nodes in x, y and z directions respectively :'
    WRITE(ounit,*) nxx, nyy, nzz
	WRITE(ounit,*) 
	WRITE(ounit,1016) 'x','x'
	WRITE(ounit,'(2X,10F7.2)')(xdel(i), i=1,nxx)
	WRITE(ounit,1016) 'y','y'
	WRITE(ounit,'(2X,10F7.2)')(ydel(j), j=1,nyy)
	WRITE(ounit,*)
	DO k= 1,np
	
	    DO j = 1, nyy
		    DO i = 1, nxx
			    IF (plnr(zpln(k))%node(i,j) == 0) THEN
				    mmap(i,j) = '  '
				ELSE	
			        WRITE (mmap(i,j),'(I2)') plnr(zpln(k))%node(i,j)
					mmap(i,j) = TRIM(ADJUSTL(mmap(i,j)))
				END IF
			END DO
		END DO
		
	    WRITE(ounit,1017) k 
		WRITE(ounit,'(5X,100I3)') (i, i = 1, nxx)
		DO j= 1, nyy
		     WRITE(ounit,'(2X,I4,100A3)') j, (mmap(i,j), i = 1 , nxx)
		END DO
	END DO
	WRITE(ounit,*) 
	WRITE(ounit,1018) 
	WRITE(ounit,*) '--------------------------------------'
    WRITE(ounit,*) '  Plane Number     Planar Region    delta-z'
	ztot = 0
    DO k= 1, nz
        DO lz= 1, zdiv(k)
	        ztot = ztot + 1
	        IF (ztot == nzz) THEN
		        WRITE(ounit,'(I9, A6, I13, F15.2)') ztot, ' (TOP)', zpln(k), zdel(ztot)
	        ELSE IF (ztot == 1) THEN
		         WRITE(ounit,'(I9, A9, I10, F15.2)') ztot, ' (BOTTOM)', zpln(k), zdel(ztot)
		    ELSE	
	            WRITE(ounit,'(I9, I19, F15.2)') ztot, zpln(k), zdel(ztot)
		    END IF
		END DO	
	END DO
	WRITE(ounit,*) 
	WRITE(ounit,*) '  Boundary conditions'
	IF (xleft == 0) THEN
	    WRITE(ounit,*)' X-directed West   (X-): VACUUM'
	ELSE
	    WRITE(ounit,*)' X-directed West   (X-): REFLECTIVE'
	END IF	
	IF (xrigt == 0) THEN
	    WRITE(ounit,*)' X-directed East   (X+): VACUUM'
	ELSE
	    WRITE(ounit,*)' X-directed East   (X+): REFLECTIVE'
	END IF	
	IF (yback == 0) THEN
	    WRITE(ounit,*)' Y-directed North  (Y-): VACUUM'
	ELSE
	    WRITE(ounit,*)' Y-directed North  (Y-): REFLECTIVE'
	END IF	
	IF (yfrnt == 0) THEN
	    WRITE(ounit,*)' Y-directed South  (Y+): VACUUM'
	ELSE
	    WRITE(ounit,*)' Y-directed South  (Y+): REFLECTIVE'
	END IF	
	IF (zbott == 0) THEN
	    WRITE(ounit,*)' Z-directed Bottom (Z-): VACUUM'
	ELSE
	    WRITE(ounit,*)' Z-directed Bottom (Z-): REFLECTIVE'
	END IF	
	IF (ztop == 0) THEN
	    WRITE(ounit,*)' Z-directed Top    (Z+): VACUUM'
	ELSE
	    WRITE(ounit,*)' Z-directed Top    (Z+): REFLECTIVE'
	END IF	
END IF	

1016 FORMAT(2X,A,'-directed nodes divison (delta-',A,')')
1017 FORMAT(2X, 'Planar Region : ', I2)
1018 FORMAT(2X, 'Planar Region Assignment to planes.')

WRITE(ounit,*)
WRITE(ounit,*) ' ...Core geometry is sucessfully read...'

DO k= 1,np
    DEALLOCATE(plnr(k)%asm)
	DEALLOCATE(plnr(k)%node)
END DO
DEALLOCATE(plnr)
DEALLOCATE(zpln)

! Assign array ix, iy, iz and xyz
nnod = 0
DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax	
		    nnod = nnod + 1
		END DO
	END DO
END DO

ALLOCATE(ix(nnod), iy(nnod), iz(nnod))
ALLOCATE(xyz(nxx, nyy, nzz))

! Assign n for boundary nodes
n = 0
DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax
		        n = n + 1
		        ix(n) = i
			    iy(n) = j
			    iz(n) = k
			    xyz(i,j,k) = n
		END DO
	END DO
END DO
	    

END SUBROUTINE inp_geom2



SUBROUTINE inp_esrc (xbunit)
!
! Purpose:
!    To read geometry card in input (1st part)
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: exsrc, ng, nx, ny, nz, nnod, ystag, &
                 xdiv, ydiv, zdiv, xyz, ix, iy, iz, &
				 xdel, ydel, zdel, nxx, nyy, nzz

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln, ios  ! line number and IOSTAT status
INTEGER :: g, i, j, k, n
INTEGER :: xt, yt, zt
INTEGER :: it, jt, kt
INTEGER :: nsrc
REAL    :: sden                                       ! Source density
REAL, DIMENSION(:), ALLOCATABLE :: spec               ! Source Spectrum
CHARACTER(LEN=1), DIMENSION(:,:), ALLOCATABLE :: spos ! Source position
INTEGER :: xpos, ypos, zpos
REAL :: summ

WRITE(ounit,*) 
WRITE(ounit,*) 
WRITE(ounit,*) '           >>>>>READING EXTRA SOURCE<<<<<'
WRITE(ounit,*) '           -------------------------------'

ALLOCATE(exsrc(nnod, ng))
exsrc = 0.d0

! Read number of source
READ(xbunit, *, IOSTAT=ios) ind, ln, nsrc
message = ' error in reading number of extra source'
CALL er_message(ounit, ios, ln, message)

ALLOCATE(spec(ng), spos(nx,ny))

DO n = 1, nsrc
    ! Read source density
    READ(xbunit, *, IOSTAT=ios) ind, ln, sden
    message = ' error in reading source density'
    CALL er_message(ounit, ios, ln, message)
	
	IF (sden <= 0.0) THEN
	    WRITE(ounit,*) '  ERROR: SOURCE DENSITY SHALL BE GREATER THAN ZERO'
		STOP
	END IF

    ! Read source spectrum
    READ(xbunit, *, IOSTAT=ios) ind, ln, (spec(g), g = 1, ng)
    message = ' error in reading source spectrum'
    CALL er_message(ounit, ios, ln, message)
	
	! Is total spectrum = 1.0?
	summ = 0.d0
	DO g = 1, ng
	    summ = summ + spec(g)
	END DO
	IF (summ /= 1.d0) THEN
	    WRITE(ounit,*) 'TOTAL SOURCE SPECTRUM AT LINE', ln, ' IS NOT EQUAL TO 1.0'
		STOP
	END IF
	
	! WRITE OUTPUT
	WRITE(ounit,'(A12,I3)') '     SOURCE ', n
	WRITE(ounit,*)         '-----------------'
	WRITE(ounit,'(A20,ES10.3, A11)') '  Source Density  : ', sden, '  n/(m^3*s)'
	WRITE(ounit,'(A19,100F6.2)') '  Source Spectrum : ', (spec(g), g = 1, ng)
	WRITE(ounit,*) ' Source Position '

	! Read source position
    DO
        READ(xbunit, *, IOSTAT=ios) ind, ln, zpos
        message = ' error in reading axial position (zpos) of extra source'
        CALL er_message(ounit, ios, ln, message)
		IF (zpos < 1) EXIT
		IF (zpos > nz) THEN
		    WRITE(ounit,* ) '  ERROR: WRONG EXTRA SOURCES POSITION (ZPOS)'
			STOP
		END IF
		spos = '0'
		DO 
            READ(xbunit, *, IOSTAT=ios) ind, ln, xpos, ypos
            message = ' error in reading radial position (xpos and ypos) of extra source'
            CALL er_message(ounit, ios, ln, message)
			IF (xpos < 1 .OR. ypos < 1) EXIT
			
			IF (xpos > nx) THEN
		        WRITE(ounit,* ) '  ERROR: WRONG EXTRA SOURCES POSITION (XPOS)'
			    STOP			
			END IF
			
			IF (ypos > ny) THEN
		        WRITE(ounit,* ) '  ERROR: WRONG EXTRA SOURCES POSITION (YPOS)'
			    STOP			
			END IF
			
			spos(xpos,ypos) = 'X'
			
			zt = 0
			kt = 1
			DO k = 1, zpos
			    IF (k > 1) kt = zt + 1
			    zt = zt + zdiv(k)
			END DO

			yt = 0
			jt = 1
			DO j = 1, ypos
			    IF (j > 1) jt = yt + 1
			    yt = yt + ydiv(j)
			END DO
			
			xt = 0
			it = 1
			DO i = 1, xpos
			    IF (i > 1) it = xt + 1
			    xt = xt + xdiv(i)
			END DO

			DO k = kt, zt
			    DO j = jt, yt
				    DO i = it, xt
						DO g = 1, ng
						    exsrc(xyz(i,j,k), g) = exsrc(xyz(i,j,k), g) + &
						    sden * spec(g)
					    END DO
					END DO
				END DO
			END DO	
			
		END DO	
		
		WRITE(ounit,'(A18,I3)') '   Plane number : ', zpos
		WRITE(ounit,'(7X,100I3)') (i, i = 1, nx)
		DO j = 1, ny
		    WRITE(ounit,'(4X,I3, 100A3 )') j, (spos(i,j), i=1, nx)
		END DO
		WRITE(ounit,*)
    END DO
END DO

DEALLOCATE(spec, spos)	

END SUBROUTINE inp_esrc



SUBROUTINE node_cx(matx)
!
! Purpose:
!    To map the CXs to nodes
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: ng, nxx, nyy, nzz, nnod, chi, nuf, sigf, &
                 D, sigr, sigs, siga, ix, iy, iz, mnum, MAT_DATA

IMPLICIT NONE

TYPE(MAT_DATA), DIMENSION(:), INTENT(IN) :: matx

INTEGER :: g, h, n, i, j, k

DO g = 1, ng
    DO n = 1, nnod 	
        D(n,g)    = matx(mnum(ix(n), iy(n), iz(n)))%D(g)
		sigr(n,g) = matx(mnum(ix(n), iy(n), iz(n)))%sigr(g)
		siga(n,g) = matx(mnum(ix(n), iy(n), iz(n)))%siga(g)
		chi(n,g)  = matx(mnum(ix(n), iy(n), iz(n)))%chi(g)
		nuf(n,g)  = matx(mnum(ix(n), iy(n), iz(n)))%nuf(g)
		sigf(n,g)  = matx(mnum(ix(n), iy(n), iz(n)))%sigf(g)
		DO h = 1, ng
		    sigs(n,g,h)  = matx(mnum(ix(n), iy(n), iz(n)))%sigs(g,h)
		END DO
    END DO		
END DO

DEALLOCATE(mnum)

END SUBROUTINE node_cx



SUBROUTINE er_message (funit, ios, ln, mess)
!
! Purpose:
!    To provide error message
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

IMPLICIT NONE

INTEGER, INTENT(IN) :: funit, ios, ln
CHARACTER(LEN=*), INTENT(IN) :: mess

IF (ios < 0) THEN
    WRITE(funit, 1013) ln
	WRITE(funit,*) mess
	1013 FORMAT(2x, 'Line', I4, ' needs more data')
	STOP
END IF
IF (ios > 0) THEN
    WRITE(funit,1004) ln
	WRITE(funit,*) mess
	1004 FORMAT(2X, 'Please check line number', I5)
    STOP
END IF

END SUBROUTINE er_message



SUBROUTINE NodPow (fn)

!
! Purpose:
!    To print axially averaged node-wise distribution
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nxx, nyy, nzz, ystag, nnod, ix, iy, iz

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: fn

REAL, DIMENSION(nxx, nyy, nzz) :: fx
INTEGER :: i, j, k, n
REAL :: summ
REAL, DIMENSION(nxx, nyy) :: fnode


fx = 0.d0
DO n = 1, nnod
    fx(ix(n), iy(n), iz(n)) = fn(n)
END DO

!Calculate radial node-wise distribution
fnode = 0.d0
DO j = 1, nyy
    DO i = ystag(j)%smin, ystag(j)%smax
	    summ = 0.d0
	    DO k = 1, nzz
	        summ = summ + fx(i,j,k)
		END DO
		fnode(i,j)= summ/DBLE(nzz)
	END DO
END DO
	
! Print
DO j = 1, nyy
	WRITE(ounit,'(100F8.3)') (fnode(i,j), i=1,nxx)
END DO

END SUBROUTINE NodPow


SUBROUTINE  AsmPow(fn)

!
! Purpose:
!    To print axially averaged assembly-wise distribution
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nx, ny, nxx, nyy, nzz, &
                xdiv, ydiv, ystag, nnod, ix, iy, iz

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: fn

REAL, DIMENSION(nxx, nyy, nzz) :: fx
INTEGER :: i, j, k, n
INTEGER :: ly, lx, ys, xs, yf, xf
REAL :: summ
REAL, DIMENSION(nxx, nyy) :: fnode
REAL, DIMENSION(nx, ny) :: fasm
REAL :: totp
INTEGER :: nfuel
REAL :: fmax
INTEGER :: xmax, ymax

fx = 0.d0
DO n = 1, nnod
    fx(ix(n), iy(n), iz(n)) = fn(n)
END DO

!Calculate axially averaged node-wise distribution
fnode = 0.d0
DO j = 1, nyy
    DO i = ystag(j)%smin, ystag(j)%smax
	    summ = 0.d0
	    DO k = 1, nzz
	        summ = summ + fx(i,j,k)
		END DO
		fnode(i,j)= summ
	END DO
END DO
      
!Calculate assembly power
nfuel = 0
totp  = 0.d0 
ys = 1
yf = 0
DO j= 1, ny
	yf = yf + ydiv(j)
	xf = 0
	xs = 1
	DO i= 1, nx
		xf = xf + xdiv(i)
		summ = 0.
		DO ly= ys, yf
			DO lx= xs, xf
	            summ = summ + fnode(lx,ly)
		    END DO
		END DO	
		fasm(i,j) = summ / DBLE(ydiv(j)*xdiv(i))
		xs = xs + xdiv(i)
		IF (fasm(i,j) > 0.d0) nfuel = nfuel + 1
		IF (fasm(i,j) > 0.d0) totp  = totp + fasm(i,j)
	END DO
	ys = ys + ydiv(j)
END DO	

! Normalize assembly power to 1.0
xmax = 0; ymax = 0
fmax = 0.d0
DO j = 1, ny
    DO i = 1, nx
	    fasm(i,j) = DBLE(nfuel) / totp * fasm(i,j)
	    IF (fasm(i,j) > fmax) THEN     ! Get max position
		    xmax = i
			ymax = j
			fmax = fasm(i,j)
		END IF
	END DO
END DO


! Print assembly power distribution
WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '    Radial Power Distribution'
WRITE(ounit,*) '  =============================='
WRITE(ounit,'(4X,100I8)') (i, i = 1, nx)
DO j= 1, ny
	WRITE(ounit,'(2X,I4,100F8.3)') j, (fasm(i,j), i=1, nx)
END DO
WRITE(ounit,*)
WRITE(ounit,*) '  MAX POS.       Maximum Value'
WRITE(ounit,1101) ymax, xmax, fasm(xmax, ymax)

1101 FORMAT(2X, '(' , I3, ',', I3,')', F15.3)


END SUBROUTINE AsmPow



SUBROUTINE  AxiPow(fn)

!
! Purpose:
!    To print axially averaged assembly-wise distribution
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nx, ny, nxx, nyy, nzz, &
                xdiv, ydiv, ystag, nnod, ix, iy, iz

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: fn

REAL, DIMENSION(nxx, nyy, nzz) :: fx
INTEGER :: i, j, k, n
INTEGER :: ly, lx, ys, xs, yf, xf
REAL :: summ
REAL, DIMENSION(nzz) :: faxi
REAL :: totp
INTEGER :: nfuel
REAL :: fmax
INTEGER :: amax

fx = 0.d0
DO n = 1, nnod
    fx(ix(n), iy(n), iz(n)) = fn(n)
END DO

! Calculate Axial Power
nfuel = 0
totp  = 0.d0 
DO k = 1, nzz
    summ = 0.d0
    DO j = 1, nyy
	    DO i = 1, nxx
		    summ = summ + fx(i,j,k)
		END DO
	END DO	
	faxi(k) = summ
    IF (faxi(k) > 0.d0) nfuel = nfuel + 1
	IF (faxi(k) > 0.d0) totp  = totp + faxi(k)
END DO

! Normalize Axial power to 1.0
fmax = 0.d0
amax = 0
DO k = 1, nzz
	faxi(k) = DBLE(nfuel) / totp * faxi(k)
	IF (faxi(k) > fmax) THEN
	    amax = k   ! Get max position
		fmax = faxi(k)
	END IF	
END DO

! Print Axial power distribution
WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '    Axial Power Distribution'
WRITE(ounit,*) '  ============================'
WRITE(ounit,*) 
WRITE(ounit,*) '    Plane Number        Power'
WRITE(ounit,*) '   ----------------------------'
DO k= nzz, 1, -1
    IF (k == nzz) THEN
	    WRITE(ounit,'(2X,I8,A7,100F13.3)') k, ' (TOP)', faxi(k)
	ELSE IF (k == 1) THEN
	    WRITE(ounit,'(2X,I8,A10,100F10.3)') k, ' (BOTTOM)', faxi(k)
	ELSE
	    WRITE(ounit,'(2X,I8,100F20.3)') k, faxi(k)
	END IF	
END DO
WRITE(ounit,*)
WRITE(ounit,*) '  MAX POS.       Maximum Value'
WRITE(ounit,1102)  amax, faxi(amax)

1101 FORMAT(2X, '(' , I3, ',', I3,')', F15.3)
1102 FORMAT(4X, '(' , I3, ')', F18.3)


END SUBROUTINE AxiPow


SUBROUTINE  AsmFlux(fn, norm)

!
! Purpose:
!    To print axially averaged assembly-wise distribution
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: ng, nx, ny, nxx, nyy, nzz, &
                xdiv, ydiv, ystag, nnod, ix, iy, iz

IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(IN) :: fn
REAL, OPTIONAL, INTENT(IN) :: norm

REAL, DIMENSION(nxx, nyy, nzz, ng) :: fx
INTEGER :: g, i, j, k, n
INTEGER :: ly, lx, ys, xs, yf, xf
REAL :: summ
REAL, DIMENSION(nxx, nyy, ng) :: fnode
REAL, DIMENSION(nx, ny, ng) :: fasm
REAL, DIMENSION(ng) :: totp
INTEGER, DIMENSION(ng) :: nfuel

fx = 0.d0
DO g = 1, ng
    DO n = 1, nnod
        fx(ix(n), iy(n), iz(n), g) = fn(n,g)
    END DO
END DO

!Calculate axially averaged node-wise distribution
fnode = 0.d0
DO g = 1, ng
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax
	        summ = 0.d0
	        DO k = 1, nzz
	            summ = summ + fx(i,j,k,g)
		    END DO
		    fnode(i,j,g)= summ
	    END DO
	END DO	
END DO
      
!Calculate Radial Flux (assembly wise)
DO g = 1, ng
    nfuel(g) = 0
    totp(g)  = 0.d0 
    ys = 1
    yf = 0
    DO j= 1, ny
	    yf = yf + ydiv(j)
	    xf = 0
	    xs = 1
	    DO i= 1, nx
		    xf = xf + xdiv(i)
		    summ = 0.
		    DO ly= ys, yf
			    DO lx= xs, xf
	                summ = summ + fnode(lx,ly,g)
		        END DO
		    END DO	
		    fasm(i,j,g) = summ / DBLE(ydiv(j)*xdiv(i))
		    xs = xs + xdiv(i)
		    IF (fasm(i,j,g) > 0.d0) nfuel(g) = nfuel(g) + 1
		    IF (fasm(i,j,g) > 0.d0) totp(g)  = totp(g) + fasm(i,j,g)
	    END DO
	    ys = ys + ydiv(j)
	END DO	
END DO	

! Normalize Flux to 1.0
IF (PRESENT(norm)) THEN
    DO g = 1, ng
        DO j = 1, ny
            DO i = 1, nx
	            fasm(i,j,g) = DBLE(nfuel(g)) / totp(g) * fasm(i,j,g) * norm
	        END DO
	    END DO	
    END DO
END IF	


! Print assembly Flux distribution
WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '    Radial Flux Distribution'
WRITE(ounit,*) '  =============================='
DO g = 1, ng
    WRITE(ounit,'(2X, A10,I3)') 'Group = ', g
    WRITE(ounit,'(2X,I13,100I11)') (i, i = 1, nx)
    DO j= 1, ny
	    WRITE(ounit,'(2X,I4,ES13.3,100ES11.3)') j, (fasm(i,j,g), i=1, nx)
    END DO
END DO	


END SUBROUTINE AsmFlux



END MODULE InpOutp
