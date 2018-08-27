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
INTEGER, PARAMETER :: wunit = 201   !write restart file
INTEGER, PARAMETER :: runit = 202   !read restart file
INTEGER, PARAMETER :: umode = 111, uxsec = 112, ugeom = 113
INTEGER, PARAMETER :: ucase = 114, uesrc = 115, uwrst = 116
INTEGER, PARAMETER :: urrst = 117, uiter = 118, uprnt = 119
INTEGER, PARAMETER :: uadf  = 120, ucrod = 121, ubcon = 122
INTEGER, PARAMETER :: uftem = 123, umtem = 124, ucden = 125
INTEGER, PARAMETER :: ucbcs = 126, uejct = 127, uther = 128
INTEGER :: bunit

INTEGER :: bmode = 0, bxsec = 0, bgeom = 0, bcase = 0, besrc = 0
INTEGER :: bwrst = 0, brrst = 0, biter = 0, bprnt = 0, badf  = 0
INTEGER :: bcrod = 0, bbcon = 0, bftem = 0, bmtem = 0, bcden = 0
INTEGER :: bcbcs = 0, bejct = 0, bther = 0

CHARACTER(LEN=100):: iline

! Geometry
INTEGER :: np                                           ! Number of planars
INTEGER, DIMENSION(:), ALLOCATABLE :: zpln              ! Planar assignment to z direction
REAL, DIMENSION(:), ALLOCATABLE :: xsize, ysize, zsize  !Assembly size
TYPE :: MAT_ASGN                                        ! Material assignment
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: asm         ! Material assignment into assembly
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: node        ! Material assignment into nodes
END TYPE
TYPE(MAT_ASGN), DIMENSION(:), ALLOCATABLE :: plnr       ! planar
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: mnum

! CROD CHANGES
INTEGER :: nstep                                                  ! Number of steps
REAL    :: coreh                                                  ! Core Height
INTEGER, DIMENSION(:,:), ALLOCATABLE :: fbmap                     ! Radial control rod bank map (node wise)
REAL :: pos0, ssize                                               ! Zero step position and step size
LOGICAL :: warn = .TRUE.                                          ! To activate warning for first time


CONTAINS

SUBROUTINE inp_read()
!
! Purpose:
!    [Main subroutine in this modul] To read input, echo the 
!    input and gives the description
!    to the user about his/her input
!


USE sdata, ONLY: ng, nnod, mode, al

IMPLICIT NONE

INTEGER :: iost, g, i
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
OPEN (UNIT=uwrst, STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=urrst, STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=uiter, STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=uprnt, STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=uadf,  STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=ucrod, STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=ubcon, STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=uftem, STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=umtem, STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=ucden, STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=ucbcs, STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=uejct, STATUS='SCRATCH', ACTION='READWRITE') 
OPEN (UNIT=uther, STATUS='SCRATCH', ACTION='READWRITE') 


CALL inp_echo()
CALL inp_comments ()

REWIND(umode)
REWIND(uxsec)
REWIND(ugeom)
REWIND(ucase)
REWIND(uesrc)
REWIND(uwrst)
REWIND(urrst)
REWIND(uiter)
REWIND(uprnt)
REWIND(uadf)
REWIND(ucrod)
REWIND(ubcon)
REWIND(uftem)
REWIND(umtem)
REWIND(ucden)
REWIND(ucbcs)
REWIND(uejct)
REWIND(uther)

! Start reading buffer files for each card

WRITE(ounit,*) 
WRITE(ounit,*) 
WRITE(ounit,1008) 
WRITE(ounit,*) &
' ***************************************************************************************************'

! Card MODE
IF (bmode == 1) THEN
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

! Card WRST
IF (bwrst == 1) CALL inp_wrst (uwrst)

! Card RRST
IF (brrst == 1) CALL inp_rrst (urrst)

! Card ITER
IF (biter == 1) CALL inp_iter (uiter)

! Card PRNT
IF (bprnt == 1) CALL inp_prnt (uprnt)

!CARD CROD
IF (bcrod == 1) CALL inp_crod (ucrod)

! Card EJCT (Rod Ejection)
IF (mode == 'RODEJECT' .AND. bejct == 1 .AND. bcrod == 1) THEN
    CALL inp_ejct(uejct)
ELSE IF (mode == 'RODEJECT' .AND. bejct /= 1) THEN
    WRITE(ounit,*) '   CALCULATION MODE ROD EJECTION'
    WRITE(ounit,1041) 'EJCT', 'ROD EJECTION - TRANSIENT'
	STOP
ELSE IF (mode == 'RODEJECT' .AND. bcrod /= 1) THEN
    WRITE(ounit,*) '   CALCULATION MODE ROD EJECTION'
    WRITE(ounit,1041) 'CROD', 'CONTROL ROD'
	STOP
ELSE IF (mode /= 'RODEJECT' .AND. bejct == 1) THEN
    WRITE(ounit,*) '   EJCT CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
	STOP
ELSE
    CONTINUE
END IF

! Miscellaneous things
CALL misc()

!CARD ADF
ALLOCATE(al(nnod,ng))
DO g = 1, ng
    DO i = 1, nnod
		al(i,g)%dc = 0.d0 ! Default alpha = 0.0
	END DO
END DO
IF (badf == 1) CALL inp_adf (uadf)

! Card ESRC
IF (mode == 'FIXEDSRC' .AND. besrc == 1) THEN
    CALL inp_esrc(uesrc)
ELSE IF (mode == 'FIXEDSRC' .AND. besrc /= 1) THEN
    WRITE(ounit,*) '   CALCULATION MODE IS FIXED SOURCE'
    WRITE(ounit,1041) 'ESRC', 'FIXED SOURCE'
	STOP
ELSE IF (mode /= 'FIXEDSRC' .AND. besrc == 1) THEN
    WRITE(ounit,*) '   ESRC CARD IS NOT NECESSARY FOR THIS CALCULATION MODE'
	STOP
ELSE
    CONTINUE
END IF

DEALLOCATE(mnum)
DO i= 1,np
    DEALLOCATE(plnr(i)%asm)
	DEALLOCATE(plnr(i)%node)
END DO
DEALLOCATE(plnr)
DEALLOCATE(zpln)

WRITE(ounit,*) 
WRITE(ounit,*) 
WRITE(ounit,*) &
' **********************************', '  STOP READING INPUT  ', '*******************************************'	

1008 FORMAT (45X, 'START READING INPUT')
1021 FORMAT(2X, 'CARD ', A, ' DOES NOT PRESENT. THIS CARD IS MANDATORY')
1041 FORMAT(2X, 'CARD ', A, ' DOES NOT PRESENT. THIS CARD IS MANDATORY FOR ', A,' CALCULATION MODE')

END SUBROUTINE inp_read


SUBROUTINE inp_echo()
!
! Purpose:
!    To rewrite the input
!

IMPLICIT NONE

INTEGER :: eof
INTEGER :: nline

WRITE(ounit, *) '                      ###########################################################'
WRITE(ounit, *) '                      #                     ADPRES 1.1                          #'
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
!    input into input buffer for each card. Comments marked by !.
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

! Start removing comments and rewrite into one input buffer
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

! Read previous input buffer and rewrite into different buffer for particular cards
! Cards identfied by %
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
		    CASE('WRST')
			    bunit = uwrst
				bwrst = 1
		    CASE('RRST')
			    bunit = urrst
				brrst = 1
		    CASE('ITER')
			    bunit = uiter
				biter = 1
		    CASE('PRNT')
			    bunit = uprnt
				bprnt = 1
		    CASE('ADF')
			    bunit = uadf
				badf = 1
		    CASE('CROD')
			    bunit = ucrod
				bcrod = 1
			CASE('EJCT')
			    bunit = uejct
				bejct = 1
		    CASE DEFAULT
			    WRITE(ounit,1014) ln, iline
				STOP
		END SELECT
    END IF
	! Write input buffer for each card
	IF (per == 0) WRITE(bunit, 1012) 'x ',ln, iline  ! 'x' used to prevent reading next line
END DO  

1012 FORMAT(A2, I5,' ',A100)
1014 FORMAT(2X, 'AT LINE', I3, ' : THIS IS A WRONG INPUT CARD : ', A8)


 
END SUBROUTINE inp_comments


SUBROUTINE inp_mode (xbunit)
!
! Purpose:
!    To read case mode in input
!

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
    CASE('RODEJECT') 
	    WRITE(ounit,1031) 'ROD EJECTION CALCULATION'
		WRITE(ounit,*)
	CASE DEFAULT
	    WRITE(ounit,1032) mode
		STOP
END SELECT		

1031 FORMAT(2X, 'CALCULATION MODE : ', A30)	
1032 FORMAT(2X, 'MODE : ', A10, ' IS UNIDENTIFIED')	

END SUBROUTINE inp_mode


SUBROUTINE inp_case (xbunit)
!
! Purpose:
!    To read case card in input
!

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

USE sdata, ONLY: nmat, ng, xsigtr, xsiga, xnuf, xsigf, &
                 xsigs, xD, xsigr, xchi, mode, ccnuf, ccsigf

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

ALLOCATE(xsigtr(nmat,ng))
ALLOCATE(xsiga (nmat,ng))
ALLOCATE(xnuf  (nmat,ng))
ALLOCATE(xsigf (nmat,ng))
ALLOCATE(xsigs (nmat,ng,ng))
ALLOCATE(xD    (nmat,ng))
ALLOCATE(xsigr (nmat,ng))
ALLOCATE(xchi  (nmat,ng))

! Reading MACROSCOPIC CXs
DO i= 1, nmat
    DO g= 1, ng
        READ(xbunit, *, IOSTAT=ios) ind, ln, xsigtr(i,g), &
		xsiga(i,g), xnuf(i,g), xsigf(i,g), &
		xchi(i,g), (xsigs(i,g,h), h = 1, ng)
		message = ' error in cross section data'
		CALL er_message(ounit, ios, ln, message)
		
		! Check CXs values
		IF (xsigtr(i,g) <= 0.0) THEN
		    WRITE(ounit,1020)i, g
			STOP
		END IF
		IF (xnuf(i,g) > 0.) ccnuf = .FALSE.
		IF (xsigf(i,g) > 0.) ccsigf = .FALSE.
		
		
		xD(i,g) = 1.d0/(3.d0*xsigtr(i,g)) 
		dum = 0.0
		DO h= 1, ng
		    IF (g /= h) dum = dum + xsigs(i,g,h)
		END DO
		xsigr(i,g) =  xsiga(i,g) + dum
	END DO
END DO


! Writing output
IF (oxsec) THEN
    DO i= 1, nmat
        WRITE(ounit,*) 
        WRITE(ounit,1009) i
	    WRITE(ounit,1011)'GROUP', 'TRANSPORT', 'DIFFUSION', 'ABSORPTION', &
	    'REMOVAL', 'NU*FISS', 'FISSION','FISS. SPCTR'
        DO g= 1, ng
		    WRITE(ounit,1010) g, xsigtr(i,g), xD(i,g), xsiga(i,g), &
		    xsigr(i,g), xnuf(i,g), xsigf(i,g), xchi(i,g)
	    END DO
	    WRITE(ounit,*)'  --SCATTERING MATRIX--'
	    WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
	    DO g= 1, ng
	        WRITE(ounit,1015)g, (xsigs(i,g,h), h=1,ng)
	    END DO
    END DO
END IF

IF (ccnuf .AND. mode /= 'FIXEDSRC') THEN
    WRITE(ounit, *) "ERROR: The Problem has no fission material (nu*fission for all materials are zero)" 
	STOP
END IF
IF (ccsigf .AND. mode /= 'FIXEDSRC') THEN
    WRITE(ounit, *) "ERROR: The Problem has no fission material (fission xsec for all materials are zero)"
	STOP
END IF

WRITE(ounit,*)
WRITE(ounit,*) ' ...Macroscopic CX Card is sucessfully read...'

1009 FORMAT(5X, 'MATERIAL', I3)
1011 FORMAT(2X, A7, A12, A13, A12, A11, 2A13, A15)
1010 FORMAT(2X, I6, F13.6, 3F12.6, 3F13.6)
1015 FORMAT(4X, I3, F16.6, 20F12.6)
1020 FORMAT(2X, 'ERROR: Transport cross section (sigtr)is zero or negative in material: ', I3, ' ;group: ', I3)

DEALLOCATE(group)

END SUBROUTINE inp_xsec


SUBROUTINE inp_geom1 (xbunit)
!
! Purpose:
!    To read geometry card in input (1st part)
!

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

! Limit values of nx, ny and nz
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
IF (nx > 33) THEN
    WRITE(ounit,*) '  Error: nx shall be maximum 33'
	STOP
END IF	
IF (ny > 33) THEN
    WRITE(ounit,*) '  Error: ny shall be maximum 33'
	STOP
END IF	
IF (nz > 40) THEN
    WRITE(ounit,*) '  Error: nz shall be maximum 40'
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

! Limit values of nxx, nyy and nzz
IF (nxx > 80) THEN
    WRITE(ounit,*) '  Error: nxx shall be maximum 80'
	STOP
END IF	
IF (nyy > 80) THEN
    WRITE(ounit,*) '  Error: nyy shall be maximum 80'
	STOP
END IF	
IF (nzz > 80) THEN
    WRITE(ounit,*) '  Error: nzz shall be maximum 80'
	STOP
END IF	

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
                xwest, xeast, ynorth, ysouth, zbott, ztop, nnod, &
				xstag, ystag, xdiv, ydiv, zdiv, nmat

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln, ios

INTEGER :: i, j, k, lx, ly, lz, xtot, ytot, ztot, n, g
CHARACTER(LEN=2), DIMENSION(nxx, nyy) :: mmap
INTEGER, PARAMETER :: xm = 36
INTEGER :: ip, ipr, kp
INTEGER :: xs, xf

! Reading number of planar
READ(xbunit, *, IOSTAT=ios) ind, ln, np
message = ' error in reading number of planars'
CALL er_message(ounit, ios, ln, message)

!Reading planar assignment into z-direction
ALLOCATE(zpln(nz))
READ(xbunit, *, IOSTAT=ios) ind, ln, (zpln(k), k=1,nz)
message = ' error in reading planar assignment'
CALL er_message(ounit, ios, ln, message)

DO k = 1, nz
    IF (zpln(k) > np) THEN
	    WRITE(ounit,'(2X,A15,I3,A35)') 'ERROR: PLANAR ', &
		zpln(k), ' IS GREATER THAN NUMBER OF PLANAR'
		STOP
	END IF
    IF (zpln(k) < 1) THEN
	    WRITE(ounit,'(2X,A)') 'ERROR: PLANAR SHOULD BE AT LEAST 1'
		STOP
	END IF
END DO

! Reading material assignment for each planar
ALLOCATE(plnr(np))
DO k= 1,np
    ALLOCATE(plnr(k)%asm (nx,ny))
	ALLOCATE(plnr(k)%node(nxx,nyy))
	DO j= ny, 1, -1
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
				        mnum(xtot, ytot, ztot) = plnr(zpln(k))%asm(i,j)
						IF (mnum(xtot, ytot, ztot) > nmat) THEN
						    WRITE(ounit,'(2X,A17,I3,A37)') 'ERROR: MATERIAL ', &
							mnum(xtot, ytot, ztot), ' IS GREATER THAN NUMBER OF MATERIAL'
						    STOP
						END IF
						IF (mnum(xtot, ytot, ztot) < 0) THEN
						    WRITE(ounit,'(2X,A)') 'ERROR: NEGATIVE MATERIAL FOUND'
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
READ(xbunit, *, IOSTAT=ios) ind, ln, xeast, xwest, ynorth, ysouth, zbott, ztop
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
	
	ip = nxx/xm
    ipr = MOD(nxx,xm) - 1
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
	    xs = 1; xf = xm
        DO kp = 1, ip
            WRITE(ounit,'(6X,100I3)') (i, i = xs, xf)
            DO j= nyy, 1, -1
	            WRITE(ounit,'(2X,I4,1X,100A3)') j, (mmap(i,j), i=xs, xf)
            END DO
	        xs = xs + xm
	        xf = xf + xm
        END DO

        WRITE(ounit,'(6X,100I3)') (i, i = xs, xs+ipr)
		IF (xs+ipr > xs) THEN
            DO j= nyy, 1, -1
	            WRITE(ounit,'(2X,I4,1X,100A3)') j, (mmap(i,j), i=xs, xs+ipr)
            END DO
		END IF
	END DO
	
	WRITE(ounit,*) 
	WRITE(ounit,1018) 
	WRITE(ounit,*) '--------------------------------------'
    WRITE(ounit,*) '  Plane Number     Planar Region    delta-z'
	ztot = nzz
    DO k= nz, 1, -1
        DO lz= 1, zdiv(k)
	        IF (ztot == nzz) THEN
		        WRITE(ounit,'(I9, A6, I13, F15.2)') ztot, ' (TOP)', zpln(k), zdel(ztot)
	        ELSE IF (ztot == 1) THEN
		         WRITE(ounit,'(I9, A9, I10, F15.2)') ztot, ' (BOTTOM)', zpln(k), zdel(ztot)
		    ELSE	
	            WRITE(ounit,'(I9, I19, F15.2)') ztot, zpln(k), zdel(ztot)
		    END IF
			ztot = ztot - 1
		END DO	
	END DO
	WRITE(ounit,*) 
	WRITE(ounit,*) '  Boundary conditions'
	
	IF (xwest == 0) THEN
	    WRITE(ounit,*)' X-directed West   : ZERO FLUX'
	ELSE IF (xwest == 1) THEN
	    WRITE(ounit,*)' X-directed West   : ZERO INCOMIG CURRENT'
	ELSE
	    WRITE(ounit,*)' X-directed West   : REFLECTIVE'
	END IF	
	
	IF (xeast == 0) THEN
	    WRITE(ounit,*)' X-directed East   : ZERO FLUX'
	ELSE IF (xeast == 1) THEN
	    WRITE(ounit,*)' X-directed East   : ZERO INCOMIG CURRENT'
	ELSE
	    WRITE(ounit,*)' X-directed East   : REFLECTIVE'
	END IF	
	
	IF (ynorth == 0) THEN
	    WRITE(ounit,*)' Y-directed North  : ZERO FLUX'
	ELSE IF (ynorth == 1) THEN
	    WRITE(ounit,*)' Y-directed North  : ZERO INCOMIG CURRENT'
	ELSE
	    WRITE(ounit,*)' Y-directed North  : REFLECTIVE'
	END IF	
	
	IF (ysouth == 0) THEN
	    WRITE(ounit,*)' Y-directed South  : ZERO FLUX'
	ELSE IF (ysouth == 1) THEN
	    WRITE(ounit,*)' Y-directed South  : ZERO INCOMIG CURRENT'
	ELSE
	    WRITE(ounit,*)' Y-directed South  : REFLECTIVE'
	END IF	
	
	IF (zbott == 0) THEN
	    WRITE(ounit,*)' Z-directed Bottom : ZERO FLUX'
	ELSE IF (zbott == 1) THEN
	    WRITE(ounit,*)' Z-directed Bottom : ZERO INCOMIG CURRENT'
	ELSE
	    WRITE(ounit,*)' Z-directed Bottom : REFLECTIVE'
	END IF	
	
	IF (ztop == 0) THEN
	    WRITE(ounit,*)' Z-directed Top    : ZERO FLUX'
	ELSE IF (ztop == 1) THEN
	    WRITE(ounit,*)' Z-directed Top    : ZERO INCOMIG CURRENT'
	ELSE
	    WRITE(ounit,*)' Z-directed Top    : REFLECTIVE'
	END IF	
END IF


1016 FORMAT(2X,A,'-directed nodes divison (delta-',A,')')
1017 FORMAT(2X, 'Planar Region : ', I2)
1018 FORMAT(2X, 'Planar Region Assignment to planes.')

WRITE(ounit,*)
WRITE(ounit,*) ' ...Core geometry is sucessfully read...'

! Calculate core height
coreh = 0.d0
DO k = 1, nzz
    coreh = coreh + zdel(k)
END DO

! set Number of nodes (nnod)
nnod = 0
DO k = 1, nzz
    DO j = 1, nyy
        DO i = ystag(j)%smin, ystag(j)%smax	
		    nnod = nnod + 1
		END DO
	END DO
END DO

END SUBROUTINE inp_geom2


SUBROUTINE misc ()
!
! Purpose:
!    To assign material xsec to nodes
!    To arranges the nodes into 1-D array instead of 3-D array
!


USE sdata, ONLY: nxx, nyy, nzz, ix, iy, iz, xyz, &
				 nnod, sigtr, siga, nuf, sigf, &
                 chi, sigs, D, sigr, ng, ystag, &
				 xdel, ydel, zdel, vdel, mode, &
				 mat, xD, xsigr, xchi, &
				 bcon, ftem, mtem, cden, bpos

IMPLICIT NONE

INTEGER :: i, j, k, n

ALLOCATE(ix(nnod), iy(nnod), iz(nnod))
ALLOCATE(xyz(nxx, nyy, nzz))
ALLOCATE(mat(nnod))

! Set ix, iy, iz and xyz
n = 0
DO k = 1, nzz
    DO j = nyy, 1, -1
        DO i = ystag(j)%smin, ystag(j)%smax
		        n = n + 1
		        ix(n) = i
			    iy(n) = j
			    iz(n) = k
			    xyz(i,j,k) = n
				mat(n) = mnum(i,j,k)
		END DO
	END DO
END DO

ALLOCATE(sigtr(nnod,ng))
ALLOCATE(siga (nnod,ng))
ALLOCATE(nuf  (nnod,ng))
ALLOCATE(sigf (nnod,ng))
ALLOCATE(sigs (nnod,ng,ng))
ALLOCATE(D    (nnod,ng))
ALLOCATE(sigr (nnod,ng))
ALLOCATE(chi  (nnod,ng))

DO n = 1, nnod
    chi  (n,:)   = xchi  (mat(n),:)
END DO

IF (mode == 'CBSEARCH' .OR. mode == 'RODEJECT') THEN
    CONTINUE
ELSE
    CALL XS_updt(bcon, ftem, mtem, cden, bpos)
END IF

DEALLOCATE(xD, xsigr, xchi)

! Calculate nodes' volume
ALLOCATE(vdel(nnod))
DO i = 1, nnod
    vdel(i) = xdel(ix(i)) * ydel(iy(i)) * zdel(iz(i)) 
END DO

END SUBROUTINE misc


SUBROUTINE XS_updt (xbcon, xftem, xmtem, xcden, xbpos)
!
! Purpose:
!    To update XS
!


USE sdata, ONLY: 

IMPLICIT NONE

REAL, INTENT(IN) :: xbcon  ! Provided Boron Concentration
REAL, DIMENSION(:), INTENT(IN) :: xftem  ! Provided fuel temperature
REAL, DIMENSION(:), INTENT(IN) :: xmtem  ! Provided moderator temperature
REAL, DIMENSION(:), INTENT(IN) :: xcden  ! Provided coolant density
REAL, DIMENSION(:), INTENT(IN) :: xbpos  ! Provided control rod bank position

CALL base_updt()
! IF (bbcon == 1 .OR. bcbcs == 1) CALL bcon_updt(xbcon)
! IF (bftem == 1) CALL ftem_updt(xftem)
! IF (bmtem == 1) CALL mtem_updt(xmtem)
! IF (bcden == 1) CALL cden_updt(xcden)
IF (bcrod == 1) CALL crod_updt(xbpos)
CALL Dsigr_updt()


END SUBROUTINE XS_updt


SUBROUTINE base_updt ()
!
! Purpose:
!    To update current XS to base XS
!


USE sdata, ONLY: nnod, sigtr, siga, nuf, sigf, sigs, mat, &
                 xsigtr, xsiga, xnuf, xsigf, xsigs

IMPLICIT NONE

INTEGER :: n



DO n = 1, nnod
    sigtr(n,1:)   = xsigtr(mat(n),1:)
    siga (n,1:)   = xsiga (mat(n),1:)
    nuf  (n,1:)   = xnuf  (mat(n),1:)
    sigf (n,1:)   = xsigf (mat(n),1:)
    sigs (n,1:,1:) = xsigs (mat(n),1:,1:)
END DO


END SUBROUTINE base_updt



SUBROUTINE Dsigr_updt ()
!
! Purpose:
!    To update diffusion constant and removal XS
!


USE sdata, ONLY: nnod, ng, sigtr, siga,  &
                 sigs, D, sigr

IMPLICIT NONE

INTEGER :: n, i, g, h
REAL :: dum

DO i = 1, nnod
    DO g = 1, ng
		D(i,g) = 1./(3.*sigtr(i,g))
		dum = 0.
		DO h= 1, ng
		    IF (g /= h) dum = dum + sigs(i,g,h)
		END DO
		sigr(i,g) =  siga(i,g) + dum
	END DO
END DO

END SUBROUTINE Dsigr_updt



SUBROUTINE inp_esrc (xbunit)
!
! Purpose:
!    To read extra sources if any
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
			WRITE(ounit, 2033) ln, zpos
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
				WRITE(ounit, 2033) ln, xpos, ypos
			    STOP			
			END IF
			
			IF (ypos > ny) THEN
		        WRITE(ounit,* ) '  ERROR: WRONG EXTRA SOURCES POSITION (YPOS)'
				WRITE(ounit, 2033) ln, xpos, ypos
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
		DO j = ny, 1, -1
		    WRITE(ounit,'(4X,I3, 100A3 )') j, (spos(i,j), i=1, nx)
		END DO
		WRITE(ounit,*)
    END DO
END DO

2033 FORMAT(2X,'LINE', I4, ' : ', I3, I3)

DEALLOCATE(spec, spos)	

END SUBROUTINE inp_esrc



SUBROUTINE inp_wrst (xbunit)

! Purpose:
!    To read restart file name (WRITE MODE) if any

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

CHARACTER(LEN=20) :: fname

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

READ(xbunit, *, IOSTAT=ios) ind, ln, fname
message = ' error in reading restart file name (WRITE)'
CALL er_message(ounit, ios, ln, message)

fname = TRIM(ADJUSTL(fname))

OPEN (UNIT=wunit, FILE=fname, STATUS='REPLACE', ACTION='WRITE')	

END SUBROUTINE inp_wrst



SUBROUTINE inp_rrst (xbunit)

! Purpose:
!    To read restart file name (READ MODE) if any


IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

CHARACTER(LEN=20) :: fname

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

READ(xbunit, *, IOSTAT=ios) ind, ln, fname
message = ' error in reading restart file name (READ)'
CALL er_message(ounit, ios, ln, message)

fname = TRIM(ADJUSTL(fname))

OPEN (UNIT=runit, FILE=fname, STATUS='REPLACE', ACTION='READ', &
      IOSTAT = ios)
	  
IF (ios /= 0) THEN  
	WRITE(ounit,*) '  ERROR - NO RESTART FILE : ', fname
	STOP
END IF	

END SUBROUTINE inp_rrst



SUBROUTINE inp_iter (xbunit)

!
! Purpose:
!    To read iteration control if any

USE sdata, ONLY: nout, nin, serc, ferc, ierc, nac

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

READ(xbunit, *, IOSTAT=ios) ind, ln, nout, nin, serc, ferc, ierc, nac
message = ' error in reading iteration control'
CALL er_message(ounit, ios, ln, message)

WRITE(ounit,*) 
WRITE(ounit,*) 
WRITE(ounit,*) '           >>>>>READING ITERATION CONTROL<<<<<'
WRITE(ounit,*) '           ------------------------------------'

WRITE(ounit,'(A,I5)') '  MAXIMUM NUMBER OF OUTER ITERATION : ', nout
WRITE(ounit,'(A,I5)') '  MAXIMUM NUMBER OF INNER ITERATION : ', nin
WRITE(ounit,'(A,ES12.3)') '  FISSION SOURCE ERROR CRITERIA     : ', serc
WRITE(ounit,'(A,ES12.3)') '  FLUX ERROR CRITERIA               : ', ferc
WRITE(ounit,'(A,ES12.3)') '  INNER ITERATION ERROR CRITERIA    : ', ierc
WRITE(ounit,'(A,I5)') '  OUTER ITERATION FISSION SOURCE EXTRAPOLATION INTERVAL : ', nac


END SUBROUTINE inp_iter



SUBROUTINE inp_prnt (xbunit)

!
! Purpose:
!    To read output print option if any

USE sdata, ONLY: aprad, apaxi, afrad

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

CHARACTER(LEN=3) :: caprad='YES', capaxi='YES', cafrad='YES'

READ(xbunit, *, IOSTAT=ios) ind, ln, aprad, apaxi, afrad
message = ' error in reading output print option'
CALL er_message(ounit, ios, ln, message)

IF (aprad == 0) caprad='NO'
IF (apaxi == 0) capaxi='NO'
IF (afrad == 0) cafrad='NO'

WRITE(ounit,*) 
WRITE(ounit,*) 
WRITE(ounit,*) '           >>>>>READING OUTPUT PRINT OPTION<<<<<'
WRITE(ounit,*) '           -------------------------------------'

WRITE(ounit,'(A,A)') '  RADIAL ASSEMBLY POWER DISTRIBUTION : ', caprad
WRITE(ounit,'(A,A)') '  AXIAL ASSEMBLY POWER DISTRIBUTION  : ', capaxi
WRITE(ounit,'(A,A)') '  RADIAL FLUX POWER DISTRIBUTION     : ', cafrad


END SUBROUTINE inp_prnt



SUBROUTINE inp_adf (xbunit)

!
! Purpose:
!    To read ADF values if any

USE sdata, ONLY: ADF_TYPE, ng, nmat, nx, ny, nz, nxx, nyy, nzz, &
                 nnod, xdiv, ydiv, zdiv, xyz, ix, iy, iz, ystag, xstag, al

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

TYPE(ADF_TYPE), DIMENSION(nmat,ng) :: mdc
TYPE(ADF_TYPE), DIMENSION(nx,ny,nz,ng) :: xdc
TYPE(ADF_TYPE), DIMENSION(nxx,nyy,nzz,ng) :: xxdc
TYPE(ADF_TYPE), DIMENSION(nxx,nyy,nzz,ng) :: xxal

INTEGER :: g, i, j, k, n
INTEGER :: rot, x1, x2, y1, y2, z1, z2
INTEGER :: xtot, ytot, ztot
INTEGER :: lz, ly, lx
INTEGER, DIMENSION(nx+1) :: tx
INTEGER, DIMENSION(ny+1) :: ty
INTEGER, DIMENSION(nz+1) :: tz
INTEGER :: zp
CHARACTER(LEN=6), DIMENSION(nx, ny) :: cadf
INTEGER, PARAMETER :: xm = 12
INTEGER :: ip, ipr
INTEGER :: xs, xf

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

WRITE(ounit,*) 
WRITE(ounit,*) 
WRITE(ounit,*) '         >>>>>READING ASSEMBLY DISCONTINUITY FACTOR<<<<<'
WRITE(ounit,*) '         -----------------------------------------------'

DO i = 1, nmat
    DO g = 1, ng
        READ(xbunit, *, IOSTAT=ios) ind, ln, (mdc(i,g)%dc(j), j = 1, 6)
        message = ' error in reading Assembly Discontinuity Factor (ADF)'
        CALL er_message(ounit, ios, ln, message)
	END DO
END DO

DO g = 1, ng
    DO k = 1, nz
        DO j = 1, ny
	        DO i = 1, nx
		        IF (plnr(zpln(k))%asm(i,j) /= 0) xdc(i,j,k,g) = mdc(plnr(zpln(k))%asm(i,j),g)
			END DO
		END DO
	END DO
END DO

!!! ADF ROTATION
DO
    READ(xbunit, *, IOSTAT=ios) ind, ln, rot
    message = ' error in reading ADF Rotation'
    CALL er_message(ounit, ios, ln, message)
	IF (rot < 1) EXIT
	IF (rot > 3) THEN
	    WRITE(ounit,*) '  ERROR: MAXIMUM ADF ROTATION IS 3 TIMES'
	    WRITE(ounit,2030) ln, rot
	END IF
	
	DO
	    READ(xbunit, *, IOSTAT=ios) ind, ln, x1, x2, y1, y2, z1, z2
        message = ' error in reading ADF Rotation'
        CALL er_message(ounit, ios, ln, message)
		IF (x1 < 1 .OR. x2 < 1 .OR. y1 < 1 .OR. y2 < 1 .OR. z1 < 1 .OR. z2 < 1) EXIT
		IF (x1 > nx .OR. x2 > nx .OR. y1 > ny .OR. y2 > ny .OR. z1 > nz .OR. z2 > nz) THEN
		    WRITE(ounit,*) '  ERROR: WRONG POSITION FOR ADF ROTATION. ' // &
			               'OUT OF DIMENSION OF THE CORE'
			WRITE(ounit,2032) ln, x1, x2, y1, y2, z1, z2
			STOP
		END IF
		IF (x2 < x1 .OR. y2 < y1 .OR. z2 < z1) THEN
		    WRITE(ounit,*) '  ERROR: WRONG POSITION FOR ADF ROTATION. ' // &
			               'INITIAL POSITION IS SMALLER'
			WRITE(ounit,2032) ln, x1, x2, y1, y2, z1, z2
		END IF
		
		DO g = 1, ng
		    DO k = z1, z2
		        DO j = y1, y2
			        DO i = x1, x2
						CALL rotate(rot, xdc(i,j,k,g)%dc(1), xdc(i,j,k,g)%dc(2), &
							             xdc(i,j,k,g)%dc(3), xdc(i,j,k,g)%dc(4))
					END DO
				END DO
			END DO
		END DO
		
	END DO
END DO

! ADF PRINT OPTION
READ(xbunit, *, IOSTAT=ios) ind, ln, zp
IF (ios == 0 .AND. zp >=1) THEN
    WRITE(ounit,*)
	WRITE(ounit,'(A,I3)') '  ADF VALUES ON PLANAR NUMBER : ', zp
	
	ip = nx/xm
    ipr = MOD(nx,xm) - 1
    DO g = 1, ng
	    WRITE(ounit,*)
	    WRITE(ounit, 1999) g
		WRITE(ounit,*)
		
		WRITE(ounit,*) '  EAST ADF'
        DO j = 1, ny
            DO i = 1, nx
		        !!! If ADF > 0, Convert ADF to character
		        IF (xdc(i,j,zp,g)%dc(1) == 0.) THEN
		            cadf(i,j) = '      '
		        ELSE
                    WRITE (cadf(i,j),'(F6.4)') xdc(i,j,zp,g)%dc(1)
			        cadf(i,j) = TRIM(ADJUSTL(cadf(i,j)))
		        END IF
	        END DO
        END DO
		
        xs = 1; xf = xm
        DO k = 1, ip
            WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
            DO j= ny, 1, -1
	            WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xf)
            END DO
            WRITE(ounit,*)
	        xs = xs + xm
	        xf = xf + xm
        END DO

        WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
		IF (xs+ipr > xs) THEN
            DO j= ny, 1, -1
	            WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xs+ipr)
            END DO
		END IF
		
		
		WRITE(ounit,*) '  WEST ADF'
        DO j = 1, ny
            DO i = 1, nx
		        !!! If ADF > 0, Convert ADF to character
		        IF (xdc(i,j,zp,g)%dc(2) == 0.) THEN
		            cadf(i,j) = '      '
		        ELSE
                    WRITE (cadf(i,j),'(F6.4)') xdc(i,j,zp,g)%dc(2)
			        cadf(i,j) = TRIM(ADJUSTL(cadf(i,j)))
		        END IF
	        END DO
        END DO
		
        xs = 1; xf = xm
        DO k = 1, ip
            WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
            DO j= ny, 1, -1
	            WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xf)
            END DO
            WRITE(ounit,*)
	        xs = xs + xm
	        xf = xf + xm
        END DO

        WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
		IF (xs+ipr > xs) THEN
            DO j= ny, 1, -1
	            WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xs+ipr)
            END DO
		END IF
		
		
		WRITE(ounit,*) '  NORTH ADF'
        DO j = 1, ny
            DO i = 1, nx
		        !!! If ADF > 0, Convert ADF to character
		        IF (xdc(i,j,zp,g)%dc(3) == 0.) THEN
		            cadf(i,j) = '      '
		        ELSE
                    WRITE (cadf(i,j),'(F6.4)') xdc(i,j,zp,g)%dc(3)
			        cadf(i,j) = TRIM(ADJUSTL(cadf(i,j)))
		        END IF
	        END DO
        END DO
		
        xs = 1; xf = xm
        DO k = 1, ip
            WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
            DO j= ny, 1, -1
	            WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xf)
            END DO
            WRITE(ounit,*)
	        xs = xs + xm
	        xf = xf + xm
        END DO

        WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
		IF (xs+ipr > xs) THEN
            DO j= ny, 1, -1
	            WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xs+ipr)
            END DO
		END IF
		
		
		WRITE(ounit,*) '  SOUTH ADF'
        DO j = 1, ny
            DO i = 1, nx
		        !!! If ADF > 0, Convert ADF to character
		        IF (xdc(i,j,zp,g)%dc(4) == 0.) THEN
		            cadf(i,j) = '      '
		        ELSE
                    WRITE (cadf(i,j),'(F6.4)') xdc(i,j,zp,g)%dc(4)
			        cadf(i,j) = TRIM(ADJUSTL(cadf(i,j)))
		        END IF
	        END DO
        END DO
		
        xs = 1; xf = xm
        DO k = 1, ip
            WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
            DO j= ny, 1, -1
	            WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xf)
            END DO
            WRITE(ounit,*)
	        xs = xs + xm
	        xf = xf + xm
        END DO

        WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
		IF (xs+ipr > xs) THEN
            DO j= ny, 1, -1
	            WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xs+ipr)
            END DO
		END IF
		
		
		WRITE(ounit,*) '  TOP ADF'
        DO j = 1, ny
            DO i = 1, nx
		        !!! If ADF > 0, Convert ADF to character
		        IF (xdc(i,j,zp,g)%dc(5) == 0.) THEN
		            cadf(i,j) = '      '
		        ELSE
                    WRITE (cadf(i,j),'(F6.4)') xdc(i,j,zp,g)%dc(5)
			        cadf(i,j) = TRIM(ADJUSTL(cadf(i,j)))
		        END IF
	        END DO
        END DO
		
        xs = 1; xf = xm
        DO k = 1, ip
            WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
            DO j= ny, 1, -1
	            WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xf)
            END DO
            WRITE(ounit,*)
	        xs = xs + xm
	        xf = xf + xm
        END DO

        WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
		IF (xs+ipr > xs) THEN
            DO j= ny, 1, -1
	            WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xs+ipr)
            END DO
		END IF
		
		
		WRITE(ounit,*) '  BOTTOM ADF'
        DO j = 1, ny
            DO i = 1, nx
		        !!! If ADF > 0, Convert ADF to character
		        IF (xdc(i,j,zp,g)%dc(6) == 0.) THEN
		            cadf(i,j) = '      '
		        ELSE
                    WRITE (cadf(i,j),'(F6.4)') xdc(i,j,zp,g)%dc(6)
			        cadf(i,j) = TRIM(ADJUSTL(cadf(i,j)))
		        END IF
	        END DO
        END DO
		
        xs = 1; xf = xm
        DO k = 1, ip
            WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
            DO j= ny, 1, -1
	            WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xf)
            END DO
            WRITE(ounit,*)
	        xs = xs + xm
	        xf = xf + xm
        END DO

        WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
		IF (xs+ipr > xs) THEN
            DO j= ny, 1, -1
	            WRITE(ounit,'(2X,I4,100A8)') j, (cadf(i,j), i=xs, xs+ipr)
            END DO
		END IF
        WRITE(ounit,*)
    END DO
END IF

1999 FORMAT (4X, 'GROUP : ', I3)

WRITE(ounit,*) '  ...Assembly Discontinuity Factors are sucessfully read...'


tx(1) = 1
DO i = 2, nx+1
    tx(i) = tx(i-1) + xdiv(i-1)
END DO
ty(1) = 1
DO j = 2, ny+1
    ty(j) = ty(j-1) + ydiv(j-1)
END DO
tz(1) = 1
DO k = 2, nz+1
    tz(k) = tz(k-1) + zdiv(k-1)
END DO

DO g = 1, ng
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
						    xxdc(xtot,ytot,ztot,g)%dc = 0.d0
							IF (mnum(xtot, ytot, ztot) /= 0) xxdc(xtot,ytot,ztot,g)%dc = 1.d0
						    IF (xtot == tx(i))     xxdc(xtot,ytot,ztot,g)%dc(2) = xdc(i,j,k,g)%dc(2)
							IF (xtot == tx(i+1)-1) xxdc(xtot,ytot,ztot,g)%dc(1) = xdc(i,j,k,g)%dc(1)
						    IF (ytot == ty(j))     xxdc(xtot,ytot,ztot,g)%dc(4) = xdc(i,j,k,g)%dc(4)
							IF (ytot == ty(j+1)-1) xxdc(xtot,ytot,ztot,g)%dc(3) = xdc(i,j,k,g)%dc(3)
							IF (ztot == tz(k))     xxdc(xtot,ytot,ztot,g)%dc(5) = xdc(i,j,k,g)%dc(5)
							IF (ztot == tz(k+1)-1) xxdc(xtot,ytot,ztot,g)%dc(6) = xdc(i,j,k,g)%dc(6)
						END DO
					END DO
				END DO
			END DO
		END DO
	END DO
END DO

! Calculate alpha
DO g = 1, ng
    DO k = 1, nzz
        DO j = 1, nyy
	        DO i = ystag(j)%smin, ystag(j)%smax
		        IF (i /= ystag(j)%smax) xxal(i,j,k,g)%dc(1) = &
				0.5d0 * (1.d0 - xxdc(i,j,k,g)%dc(1) / xxdc(i+1,j,k,g)%dc(2))
				IF (i /= ystag(j)%smin) xxal(i,j,k,g)%dc(2) = &
				0.5d0 * (1.d0 - xxdc(i,j,k,g)%dc(2) / xxdc(i-1,j,k,g)%dc(1))
				IF (j /= xstag(i)%smax) xxal(i,j,k,g)%dc(3) = &
				0.5d0 * (1.d0 - xxdc(i,j,k,g)%dc(3) / xxdc(i,j+1,k,g)%dc(4))
				IF (j /= xstag(i)%smin) xxal(i,j,k,g)%dc(4) = &
				0.5d0 * (1.d0 - xxdc(i,j,k,g)%dc(4) / xxdc(i,j-1,k,g)%dc(3))
				IF (k /= nzz) xxal(i,j,k,g)%dc(5) = &
				0.5d0 * (1.d0 - xxdc(i,j,k,g)%dc(5) / xxdc(i,j,k+1,g)%dc(6))
				IF (k /= 1) xxal(i,j,k,g)%dc(6) = &
				0.5d0 * (1.d0 - xxdc(i,j,k,g)%dc(6) / xxdc(i,j,k-1,g)%dc(5))
			END DO
		END DO
	END DO
END DO

DO g = 1, ng
    DO k = 1, nzz
        DO j = 1, nyy
	        DO i = ystag(j)%smin, ystag(j)%smax
		        al(xyz(i,j,k), g) = xxal(i,j,k,g)
			END DO
		END DO
	END DO
END DO


2030 FORMAT(3X,'LINE', I4, ' : ', I3)
2032 FORMAT(3X,'LINE', I4, ' : ', I3, I3, I3, I3, I3, I3)


END SUBROUTINE inp_adf



SUBROUTINE rotate(rot, a1, a2, a3, a4)

! Purpose:
!           To rotate ADF values (mostly necessary for BWR assemblies)

INTEGER, INTENT(IN) :: rot
REAL, INTENT(INOUT) :: a1, a2, a3, a4
REAL :: x1, x2, x3, x4


x1 = a1
x2 = a2
x3 = a3
x4 = a4

IF (rot == 1) THEN
    a1 = x4
    a2 = x3
    a3 = x1
    a4 = x2
END IF

IF (rot == 2) THEN
    a1 = x2
    a2 = x1
    a3 = x4
    a4 = x3
END IF

IF (rot == 3) THEN
    a1 = x3
    a2 = x4
    a3 = x2
    a4 = x1
END IF

END SUBROUTINE rotate


SUBROUTINE inp_crod (xbunit)

!
! Purpose:
!    To read control rod position

USE sdata, ONLY: nx, ny, nzz, nmat, ng, xdiv, ydiv, &
				 nxx, nyy, nzz, mode, bpos, nb, &
				 dsigtr, dsiga, dnuf, dsigf, dsigs

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

INTEGER :: i, j, k, g, h
REAL :: dum, dumx
INTEGER, DIMENSION(nx, ny) :: bmap       ! Radial control rod bank map (assembly wise)
INTEGER :: popt
INTEGER :: xtot, ytot, kt, ly, lx
INTEGER, DIMENSION(ng) :: group
INTEGER, DIMENSION(:), ALLOCATABLE :: bank

WRITE(ounit,*) 
WRITE(ounit,*) 
WRITE(ounit,*) '           >>>> READING CONTROL RODS INSERTION <<<<'
WRITE(ounit,*) '           --------------------------------------------'

READ(xbunit, *, IOSTAT=ios) ind, ln, nb, nstep
message = ' error in reading number of control rod bank and number of steps'
CALL er_message(ounit, ios, ln, message)

READ(xbunit, *, IOSTAT=ios) ind, ln, pos0, ssize
message = ' error in reading zeroth step rod position and step size'
CALL er_message(ounit, ios, ln, message)

ALLOCATE(bpos(nb))

!!! READ CONTROL ROD BANK POSITIONS
READ(xbunit, *, IOSTAT=ios) ind, ln, (bpos(i), i = 1, nb)
message = ' error in reading bank position'
CALL er_message(ounit, ios, ln, message)


!!! Check Control Rod Bank POSITION
DO i = 1, nb
    IF (bpos(i) > REAL(nstep)) THEN
	    WRITE(ounit,1999) 'ERROR: POSITION OF CONTROL ROD BANK ', i, ' IS ', bpos(i), ' WHICH IS HIGHER THAN NUMBER OF STEPS.'
		STOP
	END IF
    IF (bpos(i) < 0.) THEN
	    WRITE(ounit,1999) 'ERROR: POSITION OF CONTROL ROD BANK ', i, ' IS ', bpos(i), ' WHICH IS LOWER THAN 0.'
		STOP
	END IF
    IF (coreh < bpos(i)*ssize) THEN
	    WRITE(ounit,1998) 'ERROR: CORE HEIGHT ', coreh, ' IS LOWER THAN CONTROL ROD POSITION ', bpos(i)*ssize+pos0
		WRITE(ounit,*) ' BANK NUMBER ', i
		STOP
	END IF
END DO
1999 FORMAT (2X, A, I2, A, F5.1, A)
1998 FORMAT (2X, A, F6.2, A, F6.2)
	
!!! READ CONTROL ROD BANK MAP
DO j = ny, 1, -1
    READ(xbunit, *, IOSTAT=ios) ind, ln, (bmap(i,j), i = 1, nx)
    message = ' error in reading control rod bank map'
    CALL er_message(ounit, ios, ln, message)
    DO i = 1, nx
	    IF (bmap(i,j) > nb) THEN
		    WRITE(ounit,*) '  ERROR: BANK NUMBER ON CR BANK MAP IS GREATER THAN NUMBER OF BANK'
			STOP
		END IF
	END DO
END DO

ALLOCATE(dsigtr(nmat,ng))
ALLOCATE(dsiga (nmat,ng))
ALLOCATE(dnuf  (nmat,ng))
ALLOCATE(dsigf (nmat,ng))
ALLOCATE(dsigs (nmat,ng,ng))

! Reac CX changes due to control rod increment or dcrement
DO i = 1, nmat
	DO g= 1, ng
        READ(xbunit, *, IOSTAT=ios) ind, ln, dsigtr(i,g), &
		dsiga(i,g), dnuf(i,g), dsigf(i,g), (dsigs(i,g,h), h = 1, ng)
		message = ' error in reading macro xs changes due to control rod insertion'
		CALL er_message(ounit, ios, ln, message)
	END DO
END DO

!! CROD PRINT OPTION
READ(xbunit, *, IOSTAT=ios) ind, ln, popt
IF (ios == 0 .AND. popt > 0) THEN
	
	WRITE(ounit,1201) nb
	WRITE(ounit,1216) nstep
	WRITE(ounit,1202) pos0
	WRITE(ounit,1203) ssize
	
	ALLOCATE(bank(nb))
	DO i = 1, nb
	    bank(i) = i
	END DO
	WRITE(ounit,*) ' INITIAL CONTROL ROD BANK POSITION (STEPS) : '
	WRITE(ounit, 1204)(bank(i), i = 1, nb)
	WRITE(ounit, 1205)(bpos(i), i = 1, nb)
	
	WRITE(ounit,*) 
	WRITE(ounit,*) ' CONTROL ROD BANK MAP : '
	WRITE(ounit,*) ' (0 means fully inserted) '
    DO j = ny, 1, -1
        WRITE(ounit,'(100I3)' ) (bmap(i,j), i = 1, nx)
	END DO
	
	WRITE(ounit,*)
	WRITE(ounit,*) ' MATERIAL CX INCREMENT OR DECREMENT DUE TO CR INSERTION : '
    DO i= 1, nmat
       WRITE(ounit,1209) i
	    WRITE(ounit,1211)'GROUP', 'TRANSPORT', 'ABSORPTION', &
	    'NU*FISS', 'FISSION'
        DO g= 1, ng
		    WRITE(ounit,1210) g, dsigtr(i,g), dsiga(i,g), &
		    dnuf(i,g), dsigf(i,g)
			group(g) = g
	    END DO
	    WRITE(ounit,*)'  --SCATTERING MATRIX--'
	    WRITE(ounit,'(4X, A5, 20I9)') "G/G'", (group(g), g=1,ng)
	    DO g= 1, ng
	        WRITE(ounit,1215)g, (dsigs(i,g,h), h=1,ng)
	    END DO
    END DO
	DEALLOCATE(bank)
END IF	

1201 FORMAT(2X, 'NUMBER OF CONTROL ROD BANK  :', I3)
1216 FORMAT(2X, 'NUMBER OF STEPS  :', I3)
1202 FORMAT(2X, 'FULLY INSERTED POSITION (cm): ', F4.1, ' (FROM BOTTOM OF THE CORE)')
1203 FORMAT(2X, 'STEP SIZE (cm)              : ', F8.4)
1204 FORMAT(2X, 10(:, 2X, 'Bank ', I2))
1205 FORMAT(10(:, 2X, F7.1), /)
1209 FORMAT(4X, 'MATERIAL', I3)
1211 FORMAT(2X, A7, A12, A12, 2A13)
1210 FORMAT(2X, I6, F13.6, F12.6, 2F13.6)
1215 FORMAT(4X, I3, F14.6, 20F10.6)
	

!!! Convert assembly wise CR bank map to node wise CR bank map
ALLOCATE(fbmap(nxx,nyy))
ytot = 0
DO j= 1, ny
	DO ly= 1, ydiv(j)
	    ytot = ytot+1
	    xtot = 0
		DO i= 1, nx
			DO lx= 1, xdiv(i)
				 xtot = xtot+1
				 fbmap(xtot, ytot) = bmap(i,j)
			END DO
		END DO
	END DO
END DO


WRITE(ounit,*)
WRITE(ounit,*) ' ...Control Rods Insertion card is sucessfully read...'

END SUBROUTINE inp_crod


SUBROUTINE crod_updt (bpos)
!
! Purpose: TO CALCUALTE VOLUME WEIGHTED HOMOGENIZED CX FOR RODDED NODE
!    

USE sdata, ONLY: ng, nxx, nyy, nzz, xyz, zdel, nnod, mat, &
                 sigtr, siga, nuf, sigf, sigs, D, sigr, &
				 dsigtr, dsiga, dnuf, dsigf, dsigs

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: bpos

INTEGER :: n, i, j, k, g, h, kt
REAL :: rodh, vfrac
REAL :: dum, dumx
LOGICAL :: negxs = .FALSE.


DO j = 1, nyy
    DO i = 1, nxx
	    IF (fbmap(i,j) > 0) THEN
		    !!!(rodh -> posistion the tip of the control from the top of core)
		    rodh = coreh - pos0  - bpos(fbmap(i,j))*ssize
			dum = 0.d0
			DO k = nzz, 1, -1
			    ! For partially rodded node, get volume weighted homogenized CX (0 < vfrac < 1.0)
			    IF (rodh >= dum .AND. rodh < dum+zdel(k)) THEN
				    vfrac = (rodh - dum) / zdel(k)
			        sigtr(xyz(i,j,k),1:) = sigtr(xyz(i,j,k),1:) + &
				                       vfrac * dsigtr(mat(xyz(i,j,k)),1:)
			        siga(xyz(i,j,k),1:)  = siga(xyz(i,j,k),1:) + &
				                       vfrac * dsiga(mat(xyz(i,j,k)),1:)
			        nuf(xyz(i,j,k),1:)   = nuf(xyz(i,j,k),1:) + &
				                       vfrac * dnuf(mat(xyz(i,j,k)),:)
			        sigf(xyz(i,j,k),1:)  = sigf(xyz(i,j,k),1:) + &
				                       vfrac * dsigf(mat(xyz(i,j,k)),1:)
					sigs(xyz(i,j,k),1:,1:)  = sigs(xyz(i,j,k),1:,1:) + &
				                       vfrac * dsigs(mat(xyz(i,j,k)),1:,1:)				   
					
					EXIT
				END IF
				! For fully rodded node, vfrac = 1.
			    sigtr(xyz(i,j,k),:) = sigtr(xyz(i,j,k),:) + &
				                       dsigtr(mat(xyz(i,j,k)),:)
			    siga(xyz(i,j,k),:)  = siga(xyz(i,j,k),:) + &
				                       dsiga(mat(xyz(i,j,k)),:)
			    nuf(xyz(i,j,k),:)   = nuf(xyz(i,j,k),:) + &
				                       dnuf(mat(xyz(i,j,k)),:)
			    sigf(xyz(i,j,k),:)  = sigf(xyz(i,j,k),:) + &
				                       dsigf(mat(xyz(i,j,k)),:)
				sigs(xyz(i,j,k),:,:)  = sigs(xyz(i,j,k),:,:) + &
				                       dsigs(mat(xyz(i,j,k)),:,:)	
		        
			    dum = dum + zdel(k)
			END DO
			! if negative CX found, Surpress CX to zero  and calculate D and sigr
			DO k = nzz, 1, -1
				DO g = 1, ng
					IF (sigtr(xyz(i,j,k),g) < 0.) THEN
						sigtr(xyz(i,j,k),g) = 0.
						negxs = .TRUE.
					END IF
					IF (siga(xyz(i,j,k),g) < 0.) THEN
						siga(xyz(i,j,k),g) = 0.
						negxs = .TRUE.
					END IF
					IF (nuf(xyz(i,j,k),g) < 0.) THEN
						nuf(xyz(i,j,k),g) = 0.
						negxs = .TRUE.
					END IF
					IF (sigf(xyz(i,j,k),g) < 0.) THEN
						sigf(xyz(i,j,k),g) = 0.
						negxs = .TRUE.
					END IF
					DO h = 1, ng
					    IF (sigs(xyz(i,j,k),g,h) < 0.) THEN
						    sigs(xyz(i,j,k),g,h) = 0.
						    negxs = .TRUE.
					    END IF					    
					END DO
				END DO
							
                D(xyz(i,j,k),:) = 1.d0/(3.d0*sigtr(xyz(i,j,k),:)) 
			    DO g = 1, ng
					dumx = 0.0
		            DO h= 1, ng
		                IF (g /= h) dumx = dumx + sigs(xyz(i,j,k),g,h)
		            END DO
		            sigr(xyz(i,j,k),g) =  siga(xyz(i,j,k),g) + dumx
			    END DO
			END DO
		END IF
	END DO
END DO

IF (warn) THEN
    IF (negxs) WRITE(ounit,*)
    IF (negxs) WRITE(ounit,*) "  WARNING: SOME NEGATIVE CXs (DUE TO CR INSERTION) ARE SUPPRESSED TO ZERO IN SUBROUTINE CROD_UPDT"
	warn = .FALSE.
END IF


END SUBROUTINE crod_updt



SUBROUTINE inp_ejct (xbunit)

!
! Purpose:
!    To read rod ejection input

USE sdata, ONLY: nf, ng, lamb, iBeta, velo, nnod, nb, &
                 ttot, tstep1, tdiv, tstep2, bpos, fbpos, tmove, &
				 bspeed, mdir

IMPLICIT NONE

INTEGER, INTENT(IN) :: xbunit

INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

INTEGER :: i, g, n
INTEGER :: popt
INTEGER, DIMENSION(nb) :: bank
CHARACTER(LEN=4) :: cnb         ! number of bank (character type)

WRITE(ounit,*) 
WRITE(ounit,*) 
WRITE(ounit,*) '           >>>>     READING ROD EJECTION DATA      <<<<'
WRITE(ounit,*) '           --------------------------------------------'

ALLOCATE(tmove(nb), bspeed(nb), mdir(nb), fbpos(nb))

! Read Final CR bank position, time to start, and speed
DO i = 1, nb
    READ(xbunit, *, IOSTAT=ios) ind, ln, fbpos(i), tmove(i), bspeed(i)
	WRITE (cnb,'(I4)') nb
	cnb = TRIM(ADJUSTL(cnb))
    message = ' error in reading Final CR Bank Position, time to move and speed for bank : ' // cnb
    CALL er_message(ounit, ios, ln, message)
	IF (ABS(fbpos(i)-bpos(i)) < 1.d-5) THEN
	    mdir(i) = 0
	ELSE IF (fbpos(i)-bpos(i) > 1.d-5) THEN
	    mdir(i) = 2
	ELSE
	    mdir(i) = 1
	END IF
END DO

! Read time for CR to be ejected
READ(xbunit, *, IOSTAT=ios) ind, ln, ttot, tstep1, tdiv, tstep2
message = ' error in time parameters'
CALL er_message(ounit, ios, ln, message)

! ttot must be bigger than tstep1 and tstep2
IF ((ttot < tstep1) .OR. (ttot < tstep2)) THEN
    WRITE(ounit,*) 'TOTAL SIMULATION TIME SHALL BE GREATER THAN TIME STEPS'
	STOP
END IF

! tdiv must be bigger than tstep1
IF (tdiv < tstep1) THEN
    WRITE(ounit,*) 'THEN TIME WHEN SECOND TIME STEP STARTS SHALL BE GREATER THAN FIRST TIME STEP'
	STOP
END IF

! Read beta (delayed neutron fraction)
READ(xbunit, *, IOSTAT=ios) ind, ln, (iBeta(i), i = 1, nf)
message = ' error in reading delayed netron fraction (beta)'
CALL er_message(ounit, ios, ln, message)

! Read precusor decay constant
READ(xbunit, *, IOSTAT=ios) ind, ln, (lamb(i), i = 1, nf)
message = ' error in reading precusor decay constant'
CALL er_message(ounit, ios, ln, message)

! Read neutron velocity
ALLOCATE(velo(ng))
READ(xbunit, *, IOSTAT=ios) ind, ln, (velo(g), g = 1, ng)
message = ' error in reading neutron velocity'
CALL er_message(ounit, ios, ln, message)


!! EJCT PRINT OPTION
READ(xbunit, *, IOSTAT=ios) ind, ln, popt
IF (ios == 0 .AND. popt > 0) THEN
	
	DO i = 1, nb
	    bank(i) = i
	END DO
	WRITE(ounit, 1294)(bank(i), i = 1, nb)
	WRITE(ounit, 1295)(fbpos(i), i = 1, nb)
	WRITE(ounit, 1281)(tmove(i), i = 1, nb)
	WRITE(ounit, 1282)(bspeed(i), i = 1, nb)
	
	WRITE(ounit,*)
	WRITE(ounit,*) ' TIME PARAMETERS IN SECONDS : '
	WRITE(ounit,1297) ttot
	WRITE(ounit,1298) tstep1
	WRITE(ounit,1299) tstep2
	WRITE(ounit,1300) tdiv
	
	WRITE(ounit,*)
	WRITE(ounit,*) ' DELAYED NEUTRON FRACTION : '
	WRITE(ounit,'(100F11.5)') (iBeta(i), i = 1, nf)
	
	WRITE(ounit,*)
	WRITE(ounit,*) ' PRECUSOR DECAY CONSTANT (1/s): '
	WRITE(ounit,'(100F11.5)') (lamb(i), i = 1, nf)
	
	WRITE(ounit,*)
	WRITE(ounit,*) ' NEUTRON VELOCITY (cm/s) : '
	WRITE(ounit,'(100ES15.5)') (velo(g), g = 1, ng)
END IF

WRITE(ounit,*)
WRITE(ounit,*) ' ...Rod Ejection Card is sucessfully read...'

1294 FORMAT(25X, 99(:, 2X, 'Bank ', I2))
1295 FORMAT(2X, 'FINAL BANK POS. (STEP)', 99(:, 2X, F7.1), /)
1281 FORMAT(2X, 'STARTS MOVE (SECOND)  ', 99(:, 2X, F7.1), /)
1282 FORMAT(2X, 'SPEED (STEP/SECOND)   ', 99(:, 2X, F7.1), /)
1297 FORMAT(4X, 'TOTAL SIMULATION TIME         : ', F6.2)
1298 FORMAT(4X, 'FIRST TIME STEP               : ', F6.4)
1299 FORMAT(4X, 'SECOND TIME STEP              : ', F6.4)
1300 FORMAT(4X, 'WHEN SECOND TIME STEP APPLY?  : ', F6.2)

END SUBROUTINE inp_ejct



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
!    To print axially averaged node-wise power distribution
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
DO j = nyy, 1, -1
	WRITE(ounit,'(100F8.3)') (fnode(i,j), i=1,nxx)
END DO

END SUBROUTINE NodPow


SUBROUTINE  AsmPow(fn)

!
! Purpose:
!    To print axially averaged assembly-wise power distribution
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: nx, ny, nxx, nyy, nzz, zdel, &
                xdel, ydel, ystag, nnod, ix, iy, iz, &
				xdiv, ydiv

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: fn

REAL, DIMENSION(nxx, nyy, nzz) :: fx
INTEGER :: i, j, k, n
INTEGER :: ly, lx, ys, xs, yf, xf
REAL :: summ, vsumm
REAL, DIMENSION(nxx, nyy) :: fnode
REAL, DIMENSION(nx, ny) :: fasm
REAL :: totp
INTEGER :: nfuel
REAL :: fmax
INTEGER :: xmax, ymax
CHARACTER(LEN=6), DIMENSION(nx, ny) :: cpow

INTEGER, PARAMETER :: xm = 12
INTEGER :: ip, ipr

fx = 0.d0
DO n = 1, nnod
    fx(ix(n), iy(n), iz(n)) = fn(n)
END DO

!Calculate axially averaged node-wise distribution
fnode = 0.d0
DO j = 1, nyy
    DO i = ystag(j)%smin, ystag(j)%smax
	    summ = 0.d0
		vsumm = 0.d0
	    DO k = 1, nzz
	        summ = summ + fx(i,j,k)*zdel(k)
			vsumm = vsumm + zdel(k)
		END DO
		fnode(i,j)= summ/vsumm
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
		summ = 0.d0
		vsumm = 0.d0
		DO ly= ys, yf
			DO lx= xs, xf
	            summ = summ + fnode(lx,ly)*xdel(lx)*ydel(ly)
				vsumm = vsumm + xdel(lx)*ydel(ly)
		    END DO
		END DO	
		fasm(i,j) = summ / vsumm
		xs = xs + xdiv(i)
		IF (fasm(i,j) > 0.d0) nfuel = nfuel + 1
		IF (fasm(i,j) > 0.d0) totp  = totp + fasm(i,j)
	END DO
	ys = ys + ydiv(j)
END DO
	

! Normalize assembly power to 1.0
xmax = 1; ymax = 1
fmax = 0.d0
DO j = 1, ny
    DO i = 1, nx
	    IF (totp > 0.) fasm(i,j) = DBLE(nfuel) / totp * fasm(i,j) 
	    IF (fasm(i,j) > fmax) THEN     ! Get max position
		    xmax = i
			ymax = j
			fmax = fasm(i,j)
		END IF
		! Convert power to character (If power == 0 convert to blank spaces)
		IF (fasm(i,j) == 0.) THEN
		    cpow(i,j) = '     '
		ELSE
            WRITE (cpow(i,j),'(F6.3)') fasm(i,j)
			cpow(i,j) = TRIM(cpow(i,j))
		END IF
	END DO
END DO


! Print assembly power distribution
WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '    Radial Power Distribution'
WRITE(ounit,*) '  =============================='

ip = nx/xm
ipr = MOD(nx,xm) - 1
xs = 1; xf = xm
DO k = 1, ip
    WRITE(ounit,'(4X,100I8)') (i, i = xs, xf)
    DO j= ny, 1, -1
	    WRITE(ounit,'(2X,I4,100A8)') j, (cpow(i,j), i=xs, xf)
    END DO
    WRITE(ounit,*)
	xs = xs + xm
	xf = xf + xm
END DO


WRITE(ounit,'(4X,100I8)') (i, i = xs, xs+ipr)
IF (xs+ipr > xs) THEN
    DO j= ny, 1, -1
	    WRITE(ounit,'(2X,I4,100A8)') j, (cpow(i,j), i=xs, xs+ipr)
    END DO
END IF



WRITE(ounit,*)

WRITE(ounit,*) '  MAX POS.       Maximum Value'
WRITE(ounit,1101) ymax, xmax, fasm(xmax, ymax)

1101 FORMAT(2X, '(' , I3, ',', I3,')', F15.3)


END SUBROUTINE AsmPow



SUBROUTINE  AxiPow(fn)

!
! Purpose:
!    To print radially averaged  power distribution
!

USE sdata, ONLY: nx, ny, nxx, nyy, nzz, nz, zdiv, &
                vdel, ystag, nnod, ix, iy, iz, xyz, zdel, ystag

IMPLICIT NONE

REAL, DIMENSION(:), INTENT(IN) :: fn

REAL, DIMENSION(nxx, nyy, nzz) :: fx
INTEGER :: i, j, k, n, ztot
INTEGER :: ly, lx, ys, xs, yf, xf, lz
REAL :: summ, vsumm
REAL, DIMENSION(nz) :: faxi
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
ztot = 0
DO k= 1, nz
    summ = 0.d0
	vsumm = 0.d0
    DO lz= 1, zdiv(k)
	    ztot = ztot + 1
        DO j = 1, nyy
	        DO i = ystag(j)%smin, ystag(j)%smax
		        summ = summ + fx(i,j,ztot)
			    vsumm = vsumm + vdel(xyz(i,j,ztot))
		    END DO
	    END DO	
    END DO
	faxi(k) = summ/vsumm
    IF (faxi(k) > 0.d0) nfuel = nfuel + 1
	IF (faxi(k) > 0.d0) totp  = totp + faxi(k)
END DO

! Normalize Axial power to 1.0
fmax = 0.d0
amax = 1
DO k = 1, nz
	faxi(k) = REAL(nfuel) / totp * faxi(k)
	IF (faxi(k) > fmax) THEN
	    amax = k   ! Get max position
		fmax = faxi(k)
	END IF	
END DO

! Print Axial power distribution
WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '    Axial Power Density Distribution'
WRITE(ounit,*) '  ===================================='
WRITE(ounit,*) 
WRITE(ounit,*) '    Plane Number        Power      Height'
WRITE(ounit,*) '   -----------------------------------------'
summ = 0.
ztot = nzz
DO k= nz, 1, -1
    IF (k == nz) THEN
	    WRITE(ounit,'(2X,I8,A7,F13.3, F12.2)') k, ' (TOP)', faxi(k), coreh-summ
	ELSE IF (k == 1) THEN
	    WRITE(ounit,'(2X,I8,A10,F10.3, F12.2)') k, ' (BOTTOM)', faxi(k), coreh-summ
	ELSE
	    WRITE(ounit,'(2X,I8,F20.3, F12.2)') k, faxi(k), coreh-summ
	END IF	
    DO lz = 1, zdiv(k)
        summ = summ + zdel(ztot)
	    ztot = ztot - 1
	END DO
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
!    To print axially averaged assembly-wise flux distribution
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: ng, nx, ny, nxx, nyy, nzz, zdel, &
                xdel, ydel, ystag, nnod, ix, iy, iz, &
				xdiv, ydiv

IMPLICIT NONE

REAL, DIMENSION(:,:), INTENT(IN) :: fn
REAL, OPTIONAL, INTENT(IN) :: norm

REAL, DIMENSION(nxx, nyy, nzz, ng) :: fx
INTEGER :: g, i, j, k, n
INTEGER :: ly, lx, ys, xs, yf, xf
REAL :: summ, vsumm
REAL, DIMENSION(nxx, nyy, ng) :: fnode
REAL, DIMENSION(nx, ny, ng) :: fasm
REAL, DIMENSION(ng) :: totp
INTEGER, DIMENSION(ng) :: nfuel
CHARACTER(LEN=10), DIMENSION(nx, ny) :: cflx

INTEGER, PARAMETER :: xm = 12
INTEGER :: ip, ipr
INTEGER :: negf

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
			vsumm = 0.d0
	        DO k = 1, nzz
	            summ = summ + fx(i,j,k,g)*xdel(i)*ydel(j)*zdel(k)
				vsumm = vsumm + xdel(i)*ydel(j)*zdel(k)
		    END DO
		    fnode(i,j,g)= summ/vsumm
	    END DO
	END DO	
END DO
      
!Calculate Radial Flux (assembly wise)
negf = 0
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
		    summ = 0.d0
			vsumm = 0.d0
		    DO ly= ys, yf
			    DO lx= xs, xf
	                summ = summ + fnode(lx,ly,g)*xdel(lx)*ydel(ly)
					vsumm = vsumm + xdel(lx)*ydel(ly)
		        END DO
		    END DO	
		    fasm(i,j,g) = summ / vsumm
		    xs = xs + xdiv(i)
		    IF (fasm(i,j,g) > 0.d0) THEN
			    nfuel(g) = nfuel(g) + 1
		        totp(g)  = totp(g) + fasm(i,j,g)
			END IF
			! Check if there is negative flux
			IF (fasm(i,j,g) < 0.d0) negf = 1
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
IF (negf > 0) WRITE(ounit,*) '    ....WARNING: NEGATIVE FLUX ENCOUNTERED....'
WRITE(ounit,*)
WRITE(ounit,*) '    Radial Flux Distribution'
WRITE(ounit,*) '  =============================='

ip = nx/xm
ipr = MOD(nx,xm) - 1

DO g = 1, ng
    WRITE(ounit,'(A,I3)') '    Group : ', g
    !!! Convert to character (zero flux convert to blank spaces)
    DO j = 1, ny
        DO i = 1, nx
		    IF (fasm(i,j,g) == 0.) THEN
		        cflx(i,j) = '         '
		    ELSE
                WRITE (cflx(i,j),'(ES10.3)') fasm(i,j,g)
			    cflx(i,j) = TRIM(ADJUSTL(cflx(i,j)))
		    END IF
		END DO
	END DO	

	xs = 1; xf = xm
    DO k = 1, ip
        WRITE(ounit,'(3X,100I11)') (i, i = xs, xf)
        DO j= ny, 1, -1
	        WRITE(ounit,'(2X,I4,2X,100A11)') j, (cflx(i,j), i=xs, xf)
        END DO
        WRITE(ounit,*)
	    xs = xs + xm
	    xf = xf + xm
    END DO

    WRITE(ounit,'(3X,100I11)') (i, i = xs, xs+ipr)
	IF (xs+ipr > xs) THEN
        DO j= ny, 1, -1
	        WRITE(ounit,'(2X,I4,2X,100A11)') j, (cflx(i,j), i=xs, xs+ipr)
        END DO
	END IF
   WRITE(ounit,*)
   
END DO


END SUBROUTINE AsmFlux



SUBROUTINE  w_rst()

!
! Purpose:
!    To Write Restart File
!
!   Date                Programmer           History
!  ========================================================
!  6 FEB 2017         Muhammad Imron       Original code

USE sdata, ONLY: ng, nnod, f0, fx1, fy1, fz1, fx2, fy2, fz2, &
                 nod, Ke

IMPLICIT NONE

INTEGER :: n, g

WRITE(wunit, '(I4, I12)') ng, nnod
WRITE(wunit, '(F13.8)') Ke
DO g = 1, ng
    DO n = 1, nnod
        WRITE(wunit, '(7ES14.5)') f0(n,g), fx1(n,g), fy1(n,g), &
		                          fz1(n,g), fx2(n,g), fy2(n,g), fz2(n,g)
    END DO
END DO
DO g = 1, ng
    DO n = 1, nnod
        WRITE(wunit, '(6ES14.5)') nod(n,g)%ji(1), nod(n,g)%ji(2), nod(n,g)%ji(3), &
		                          nod(n,g)%ji(4), nod(n,g)%ji(5), nod(n,g)%ji(6)
        WRITE(wunit, '(6ES14.5)') nod(n,g)%jo(1), nod(n,g)%jo(2), nod(n,g)%jo(3), &
		                          nod(n,g)%jo(4), nod(n,g)%jo(5), nod(n,g)%jo(6)
    END DO
END DO

END SUBROUTINE w_rst




END MODULE InpOutp
