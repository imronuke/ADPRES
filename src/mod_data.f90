MODULE sdata

SAVE

CHARACTER(LEN=100) :: mode

INTEGER :: ng     ! number of groups
INTEGER :: nmat   ! number of materials
!! CXs Assigned to Nodes
REAL, DIMENSION(:,:), ALLOCATABLE :: sigtr          ! Transport macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: siga           ! Absorption macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: nuf            ! nu* fission macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: sigf           ! fission macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: chi            ! neutron fission spectrum
REAL, DIMENSION(:,:,:), ALLOCATABLE :: sigs         ! Scattering macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: D              ! Diffusion coefficient
REAL, DIMENSION(:,:), ALLOCATABLE :: sigr           ! Removal macroscopic cx

!! CXs Assigned to Materials
REAL, DIMENSION(:,:), ALLOCATABLE :: xsigtr          ! Transport macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: xsiga           ! Absorption macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: xnuf            ! nu* fission macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: xsigf           ! fission macroscopic cx
REAL, DIMENSION(:,:), ALLOCATABLE :: xchi            ! neutron fission spectrum
REAL, DIMENSION(:,:), ALLOCATABLE :: xD              ! Diffusion coefficient
REAL, DIMENSION(:,:), ALLOCATABLE :: xsigr           ! Removal macroscopic cx
REAL, DIMENSION(:,:,:), ALLOCATABLE :: xsigs         ! Scattering macroscopic cx
LOGICAL :: ccnuf = .TRUE.                            ! Logical variable to check the presence of fissile material
LOGICAL :: ccsigf = .TRUE.                           ! Logical variable to check the presence of fissile material

! Geometry 
INTEGER :: nx, ny, nz                                ! Number of assemblies in x, y, and z directions
INTEGER :: nxx, nyy, nzz                             ! Number of nodes in x, y, and z directions 
INTEGER :: nnod                                      ! Number of nodes
INTEGER, DIMENSION(:), ALLOCATABLE :: ix, iy, iz
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: xyz
INTEGER, DIMENSION(:), ALLOCATABLE :: xdiv, ydiv, zdiv     ! Assembly division
REAL, DIMENSION(:), ALLOCATABLE :: xdel, ydel, zdel, vdel  ! Delta x, y and z and nodes' volume in cm3
INTEGER :: xwest, xeast, ysouth, ynorth, zbott, ztop       ! Boundary conditions
INTEGER, DIMENSION(:), ALLOCATABLE :: mat

! Keff, flux and currents
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
INTEGER :: nout = 500
INTEGER :: nac = 5

! OUTPUT PRINT OPTION
INTEGER :: aprad=1, apaxi=1, afrad=1

!ADF
TYPE :: ADF_TYPE
	REAL, DIMENSION(6) :: dc            
END TYPE
TYPE(ADF_TYPE), DIMENSION(:,:), ALLOCATABLE :: al

! FUEL TEMPERATURE
REAL, DIMENSION(:), ALLOCATABLE :: ftem       ! Fuel temperature in Kelvin for each nodes
REAL :: rftem      ! Fuel temperature Reference in Kelvin
REAL, DIMENSION(:,:), ALLOCATABLE :: fsigtr, fsiga, fnuf, fsigf   ! CX changes per fuel temp changes
REAL, DIMENSION(:,:,:), ALLOCATABLE :: fsigs

! MODERATOR TEMPERATURE
REAL, DIMENSION(:), ALLOCATABLE :: mtem       ! Moderator temperature in Kelvin for each nodes
REAL :: rmtem      ! Moderator temperature Reference in Kelvin
REAL, DIMENSION(:,:), ALLOCATABLE :: msigtr, msiga, mnuf, msigf   ! CX changes per Moderator temp changes
REAL, DIMENSION(:,:,:), ALLOCATABLE :: msigs

! COOLANT DENSITY
REAL, DIMENSION(:), ALLOCATABLE :: cden       ! Coolant Density in g/cm3 for each nodes
REAL :: rcden      ! Coolant Density Reference in g/cm3
REAL, DIMENSION(:,:), ALLOCATABLE :: lsigtr, lsiga, lnuf, lsigf   ! CX changes per Coolant density changes
REAL, DIMENSION(:,:,:), ALLOCATABLE :: lsigs

! Crod changes
INTEGER :: nb                                                     ! Number of CR banks
REAL, DIMENSION(:), ALLOCATABLE :: bpos  ! CR bank position
REAL, DIMENSION(:), ALLOCATABLE :: fbpos    ! Final CR bank position
REAL, DIMENSION(:,:), ALLOCATABLE :: dsigtr, dsiga, dnuf, dsigf   ! CX incerement or decrement due to CR insertion
REAL, DIMENSION(:,:,:), ALLOCATABLE :: dsigs
REAL, DIMENSION(:), ALLOCATABLE :: tmove    ! Time when CR bank starts moving
REAL, DIMENSION(:), ALLOCATABLE :: bspeed   ! CR bank movement speed
INTEGER, DIMENSION(:), ALLOCATABLE :: mdir     ! To indicate CR movement direction (0=do not move, 1=down, 2 = up)

! Boron Concentration
REAL :: bcon       ! Boron concentration in ppm
REAL :: rbcon      ! Boron concentration in ppm Reference
REAL, DIMENSION(:,:), ALLOCATABLE :: csigtr, csiga, cnuf, csigf   ! CX changes due to boron concentration
REAL, DIMENSION(:,:,:), ALLOCATABLE :: csigs                      ! Used only for CBCS card

! Transient parameters
INTEGER, PARAMETER :: nf = 6                       ! Number of delaye dneutron precusor family 
REAL, DIMENSION(nf) :: ibeta, lamb                 ! beta and precusor decay constant
! REAL, DIMENSION(:,:), ALLOCATABLE :: iC            ! neutron precusor density
REAL, DIMENSION(:), ALLOCATABLE :: velo            ! Neutron velocity
REAL :: ttot                                       ! TOTAL SIMULATION TIME
REAL :: tstep1                                     ! FIRST TIME STEP
REAL :: tstep2                                     ! SECOND TIME STEP
REAL :: tdiv                                       ! WHEN SECOND TIME STEP APPLY


END MODULE sdata