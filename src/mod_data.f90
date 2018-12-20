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
REAL, DIMENSION(:), ALLOCATABLE :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2      ! Fission source moments
REAL, DIMENSION(:,:), ALLOCATABLE :: c0, cx1, cy1, cz1, cx2, cy2, cz2  ! neutron precusor density
REAL, DIMENSION(:,:), ALLOCATABLE :: ct, ctx1, cty1, ctz1, ctx2, cty2, ctz2 ! previous neutron precusor density

TYPE :: STAGGERED
    INTEGER :: smax, smin                             ! imax and imin along x and y direction for staggered nodes
END TYPE
TYPE(STAGGERED), DIMENSION(:), ALLOCATABLE :: ystag, xstag

! Extra Sources
REAL, DIMENSION(:,:), ALLOCATABLE :: exsrc

! Iteration Control
REAL :: ferc = 1.e-5    ! Flux Error Criteria
REAL :: serc = 1.e-5    ! Fission source Error CRITERIA
REAL :: ierc = 1.e-5    ! Inner Iteration Error Criteria
REAL :: fer, ser        ! Flux and Fission source error in BCSEARCH calcs.
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
LOGICAL :: negxs = .FALSE.                      ! To activate warning for first time
INTEGER :: cusp = 0

! Boron Concentration
REAL :: bcon       ! Boron concentration in ppm
REAL :: rbcon      ! Boron concentration in ppm Reference
REAL, DIMENSION(:,:), ALLOCATABLE :: csigtr, csiga, cnuf, csigf   ! CX changes due to boron concentration
REAL, DIMENSION(:,:,:), ALLOCATABLE :: csigs                      ! Used only for CBCS card

! Transient parameters
INTEGER, PARAMETER :: nf = 6                       ! Number of delaye dneutron precusor family
REAL, DIMENSION(nf) :: ibeta, lamb                 ! beta (delayed neutron fraction) and precusor decay constant
REAL :: tbeta                                      ! total beta
REAL, DIMENSION(:), ALLOCATABLE :: velo            ! Neutron velocity
REAL :: ttot                                       ! TOTAL SIMULATION TIME
REAL :: tstep1                                     ! FIRST TIME STEP
REAL :: tstep2                                     ! SECOND TIME STEP
REAL :: tdiv                                       ! WHEN SECOND TIME STEP APPLY
REAL, DIMENSION(:,:), ALLOCATABLE :: omeg          ! Exponential transformation constant

! Thermal-hydraulics parameters
REAL :: pow                                        ! Reactor power for given geometry (watt)
REAL :: ppow                                       ! Reactor percent power in percent
REAL :: tpow                                       ! Total reactor power
REAL, DIMENSION(:), ALLOCATABLE :: npow            ! nodes power (watt)
REAL :: tin                                        ! coolant inlet temperature (kelvin)
REAL :: cflow                                      ! Sub-channel mass flow rate (kg/s)
REAL :: rf, tg, tc, ppitch                         ! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
REAL :: dia, dh, farea                             ! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
REAL :: cf                                         ! heat fraction deposited into coolant
REAL, DIMENSION(:,:), ALLOCATABLE :: node_nf       ! Number of fuel pin per node
INTEGER :: nm                                      ! number of Fuel meat mesh
INTEGER :: nt                                      ! Number Total mesh
REAL, DIMENSION(:,:), ALLOCATABLE :: tfm           ! Fuel pin mesh temperature for each nodes
REAL, DIMENSION(:), ALLOCATABLE :: rdel            ! mesh delta
REAL, DIMENSION(:), ALLOCATABLE :: rpos            ! mesh position
REAL :: th_err                                     ! Doppler error
REAL, DIMENSION(:), ALLOCATABLE :: ent             ! Coolant Enthalpy (J/Kg)
REAL, DIMENSION(:), ALLOCATABLE :: heatf           ! Heat flux (W/m2)
INTEGER :: th_niter = 20                                ! Maximum number of thermal-hydraulics iteration
INTEGER, PARAMETER :: thunit = 300                 ! Unit number to open steam table file
REAL, PARAMETER :: pi = 3.14159265

! Steam Table data
INTEGER, PARAMETER:: ntem = 16   ! Number of temperature in steam table
REAL, DIMENSION(ntem,6) :: stab


END MODULE sdata
