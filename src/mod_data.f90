MODULE sdata

SAVE

CHARACTER(LEN=100) :: mode

INTEGER :: ng     ! number of groups
INTEGER :: nmat   ! number of materials
! MAterial CX data
REAL, DIMENSION(:,:), ALLOCATABLE :: sigtr          ! Transport macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: siga           ! Absorption macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: nuf            ! nu* fission macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: sigf           ! fission macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: chi            ! neutron fission spectrum
REAL, DIMENSION(:,:,:), ALLOCATABLE :: sigs         ! Scattering macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: D              ! Diffusion coefficient
REAL, DIMENSION(:,:), ALLOCATABLE :: sigr           ! Removal macroscopic cx

 
INTEGER :: nx, ny, nz                                ! Number of assemblies in x, y, and z directions
INTEGER :: nxx, nyy, nzz                             ! Number of nodes in x, y, and z directions 
INTEGER :: nnod                                      ! Number of nodes
INTEGER, DIMENSION(:), ALLOCATABLE :: ix, iy, iz
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: xyz
INTEGER, DIMENSION(:), ALLOCATABLE :: xdiv, ydiv, zdiv  ! Assembly division
REAL, DIMENSION(:), ALLOCATABLE :: xdel, ydel, zdel  ! Delta x, y and z
INTEGER :: xwest, xeast, ysouth, ynorth, zbott, ztop   ! Boundary conditions

REAL :: Ke
TYPE :: NODE_DATA
	REAL, DIMENSION(6) :: jo             ! Nodals' outgoing currents (X+,X-,Y+, Y-, Z+, Z-)
	REAL, DIMENSION(6) :: ji             ! Nodals' ingoing currents  (X+,X-,Y+, Y-, Z+, Z-)
	REAL, DIMENSION(3) :: L              ! Zeroth transverse leakages (Lx, Ly, Lz)
	REAL, DIMENSION(7) :: Q              ! Nodal's source and source moments (0, x1, y1, z1, x2, y2, z2)
    REAL, DIMENSION(6,6) :: P            ! Response matrix
    REAL, DIMENSION(6,7) :: R            ! Response matrix	
END TYPE
TYPE(NODE_DATA), DIMENSION(:,:), ALLOCATABLE :: nod

REAL, DIMENSION(:,:), ALLOCATABLE :: f0, fx1, fy1, fz1, fx2, fy2, fz2      ! Flux and Flux moments

TYPE :: STAGGERED
	INTEGER :: smax, smin                             ! imax and imin along x and y direction for staggered nodes
END TYPE
TYPE(STAGGERED), DIMENSION(:), ALLOCATABLE :: ystag, xstag 

! Extra Sources
REAL, DIMENSION(:,:), ALLOCATABLE :: exsrc

! Iteration Control
REAL :: ferc = 1.e-5
REAL :: serc = 1.e-5
REAL :: ierc = 1.e-5
INTEGER :: nin = 2
INTEGER :: nout = 300
INTEGER :: nac = 5

! OUTPUT PRINT OPTION
INTEGER :: aprad=1, apaxi=1, afrad=1

!ADF
TYPE :: ADF_TYPE
	REAL, DIMENSION(6) :: dc            
END TYPE
TYPE(ADF_TYPE), DIMENSION(:,:), ALLOCATABLE :: al


END MODULE sdata
