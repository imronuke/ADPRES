MODULE sdata

SAVE

CHARACTER(LEN=100) :: mode

INTEGER :: ng     ! number of groups
INTEGER :: nmat   ! number of materials
! Material data 
TYPE :: MAT_DATA
	REAL, DIMENSION(:), ALLOCATABLE :: sigtr          ! Transport macroscopic cx
	REAL, DIMENSION(:), ALLOCATABLE :: siga           ! Absorption macroscopic cx
	REAL, DIMENSION(:), ALLOCATABLE :: nuf            ! nu* fission macroscopic cx
	REAL, DIMENSION(:), ALLOCATABLE :: sigf           ! fission macroscopic cx
	REAL, DIMENSION(:), ALLOCATABLE :: chi            ! neutron fission spectrum
	REAL, DIMENSION(:,:), ALLOCATABLE :: sigs         ! Scattering macroscopic cx
	REAL, DIMENSION(:), ALLOCATABLE :: D              ! Diffusion coefficient
	REAL, DIMENSION(:), ALLOCATABLE :: sigr           ! Removal macroscopic cx
END TYPE
TYPE(MAT_DATA), DIMENSION(:), ALLOCATABLE :: mat

! Nodes' CX data
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: mnum
REAL, DIMENSION(:,:), ALLOCATABLE :: chi, nuf, sigf, siga, sigr, D
REAL, DIMENSION(:,:,:), ALLOCATABLE :: sigs

 
INTEGER :: nx, ny, nz                                ! Number of assemblies in x, y, and z directions
INTEGER :: nxx, nyy, nzz                             ! Number of nodes in x, y, and z directions 
INTEGER :: nnod                                      ! Number of nodes
INTEGER, DIMENSION(:), ALLOCATABLE :: ix, iy, iz
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: xyz
INTEGER, DIMENSION(:), ALLOCATABLE :: xdiv, ydiv, zdiv  ! Assembly division
REAL, DIMENSION(:), ALLOCATABLE :: xdel, ydel, zdel  ! Delta x, y and z
INTEGER :: xleft, xrigt, yback, yfrnt, zbott, ztop   ! Boundary conditions

REAL :: Ke
TYPE :: NODE_DATA
	REAL, DIMENSION(:), ALLOCATABLE :: jo             ! Nodals' outgoing currents (X+,X-,Y+, Y-, Z+, Z-)
	REAL, DIMENSION(:), ALLOCATABLE :: ji             ! Nodals' ingoing currents  (X+,X-,Y+, Y-, Z+, Z-)
	REAL, DIMENSION(:), ALLOCATABLE :: L               ! Zeroth transverse leakages (Lx, Ly, Lz)
	REAL, DIMENSION(:), ALLOCATABLE :: Q               ! Nodal's source and source moments (0, x1, y1, z1, x2, y2, z2)
    REAL, DIMENSION(:,:), ALLOCATABLE :: P, R          ! Nodal coupling coefficients matrix (0, x1, y1, z1, x2, y2, z2)
END TYPE
TYPE(NODE_DATA), DIMENSION(:,:), ALLOCATABLE :: nod

REAL, DIMENSION(:,:), ALLOCATABLE :: f0, fx1, fy1, fz1, fx2, fy2, fz2      ! Flux and Flux moments

TYPE :: STAGGERED
	INTEGER :: smax, smin                             ! imax and imin along x and y direction for staggered nodes
END TYPE
TYPE(STAGGERED), DIMENSION(:), ALLOCATABLE :: ystag, xstag 

REAL, DIMENSION(:,:), ALLOCATABLE :: exsrc


END MODULE sdata
