MODULE sdata

SAVE

CHARACTER(LEN=100) :: mode

INTEGER :: ng     ! number of groups
INTEGER :: nmat   ! number of materials
!! CXs Assigned to Nodes
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: sigtr          ! Transport macroscopic cx
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: siga           ! Absorption macroscopic cx
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: nuf            ! nu* fission macroscopic cx
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: sigf           ! fission macroscopic cx
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: chi            ! neutron fission spectrum
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: sigs         ! Scattering macroscopic cx
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: D              ! Diffusion coefficient
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: sigr           ! Removal macroscopic cx

!! CXs Assigned to Materials
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xsigtr          ! Transport macroscopic cx
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xsiga           ! Absorption macroscopic cx
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xnuf            ! nu* fission macroscopic cx
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xsigf           ! fission macroscopic cx
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xD              ! Diffusion coefficient
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xsigr           ! Removal macroscopic cx
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: xsigs         ! Scattering macroscopic cx
LOGICAL :: ccnuf = .TRUE.                            ! Logical variable to check the presence of fissile material
LOGICAL :: ccsigf = .TRUE.                           ! Logical variable to check the presence of fissile material

! Geometry
INTEGER :: nx, ny, nz                                ! Number of assemblies in x, y, and z directions
INTEGER :: nxx, nyy, nzz                             ! Number of nodes in x, y, and z directions
INTEGER :: nnod                                      ! Number of nodes
INTEGER, DIMENSION(:), ALLOCATABLE :: ix, iy, iz
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: xyz
INTEGER, DIMENSION(:), ALLOCATABLE :: xdiv, ydiv, zdiv     ! Assembly division
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xdel, ydel, zdel, vdel  ! Delta x, y and z and nodes' volume in cm3
INTEGER :: xwest, xeast, ysouth, ynorth, zbott, ztop       ! Boundary conditions
INTEGER, DIMENSION(:), ALLOCATABLE :: mat

! Keff, flux and currents
DOUBLE PRECISION :: Ke
TYPE :: NODE_DATA
    DOUBLE PRECISION, DIMENSION(6) :: jo             ! Nodals' outgoing currents (X+,X-,Y+, Y-, Z+, Z-)
    DOUBLE PRECISION, DIMENSION(6) :: ji             ! Nodals' ingoing currents  (X+,X-,Y+, Y-, Z+, Z-)
    DOUBLE PRECISION, DIMENSION(3) :: L              ! Zeroth transverse leakages (Lx, Ly, Lz)
    DOUBLE PRECISION, DIMENSION(7) :: Q              ! Nodal's source and source moments (0, x1, y1, z1, x2, y2, z2)
    DOUBLE PRECISION, DIMENSION(6,6) :: P            ! Response matrix
    DOUBLE PRECISION, DIMENSION(6,7) :: R            ! Response matrix
END TYPE
TYPE(NODE_DATA), DIMENSION(:,:), ALLOCATABLE :: nod

DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: f0, fx1, fy1, fz1, fx2, fy2, fz2      ! Flux and Flux moments
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2      ! Fission source moments
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: c0, cx1, cy1, cz1, cx2, cy2, cz2  ! neutron precusor density
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Parameters at previous time step

TYPE :: STAGGERED
    INTEGER :: smax, smin                             ! imax and imin along x and y direction for staggered nodes
END TYPE
TYPE(STAGGERED), DIMENSION(:), ALLOCATABLE :: ystag, xstag

! Extra Sources
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: exsrc

! Iteration Control
DOUBLE PRECISION :: ferc = 1.e-5    ! Flux Error Criteria
DOUBLE PRECISION :: serc = 1.e-5    ! Fission source Error CRITERIA
DOUBLE PRECISION :: ierc = 1.e-5    ! Inner Iteration Error Criteria
DOUBLE PRECISION :: fer, ser        ! Flux and Fission source error in BCSEARCH calcs.
INTEGER :: nin = 4      ! Maximum inner iteration
INTEGER :: nout = 500   ! Maximum outer iteration
INTEGER :: nac = 5      ! number of outer iteration before next source EXTRAPOLATION
INTEGER :: th_niter = 20                           ! Maximum number of thermal-hydraulics iteration

! OUTPUT PRINT OPTION
INTEGER :: aprad=1, apaxi=1, afrad=1

!ADF
TYPE :: ADF_TYPE
    DOUBLE PRECISION, DIMENSION(6) :: dc
END TYPE
TYPE(ADF_TYPE), DIMENSION(:,:), ALLOCATABLE :: al

! FUEL TEMPERATURE
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ftem       ! Fuel temperature in Kelvin for each nodes
DOUBLE PRECISION :: rftem      ! Fuel temperature Reference in Kelvin
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: fsigtr, fsiga, fnuf, fsigf   ! CX changes per fuel temp changes
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: fsigs

! MODERATOR TEMPERATURE
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mtem       ! Moderator temperature in Kelvin for each nodes
DOUBLE PRECISION :: rmtem      ! Moderator temperature Reference in Kelvin
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: msigtr, msiga, mnuf, msigf   ! CX changes per Moderator temp changes
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: msigs

! COOLANT DENSITY
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cden       ! Coolant Density in g/cm3 for each nodes
DOUBLE PRECISION :: rcden      ! Coolant Density Reference in g/cm3
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: lsigtr, lsiga, lnuf, lsigf   ! CX changes per Coolant density changes
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: lsigs

! Crod changes
INTEGER :: nb                                                     ! Number of CR banks
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: bpos  ! CR bank position
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fbpos    ! Final CR bank position
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: dsigtr, dsiga, dnuf, dsigf   ! CX incerement or decrement due to CR insertion
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: dsigs
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: tmove    ! Time when CR bank starts moving
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: bspeed   ! CR bank movement speed
INTEGER, DIMENSION(:), ALLOCATABLE :: mdir  ! To indicate CR movement direction (0=do not move, 1=down, 2 = up)
INTEGER :: cusp = 0                         ! Rod cusping option

! Boron Concentration
DOUBLE PRECISION :: bcon       ! Boron concentration in ppm
DOUBLE PRECISION :: rbcon      ! Boron concentration in ppm Reference
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: csigtr, csiga, cnuf, csigf   ! CX changes due to boron concentration
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: csigs                      ! Used only for CBCS card

! Transient parameters
INTEGER, PARAMETER :: nf = 6                       ! Number of delaye dneutron precusor family
DOUBLE PRECISION, DIMENSION(nf) :: ibeta, lamb                 ! beta (delayed neutron fraction) and precusor decay constant
DOUBLE PRECISION :: tbeta                                      ! total beta
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: velo            ! Neutron velocity
DOUBLE PRECISION :: ttot                                       ! TOTAL SIMULATION TIME
DOUBLE PRECISION :: tstep1                                     ! FIRST TIME STEP
DOUBLE PRECISION :: tstep2                                     ! SECOND TIME STEP
DOUBLE PRECISION :: tdiv                                       ! WHEN SECOND TIME STEP APPLY
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: omeg          ! Exponential transformation constant
LOGICAL :: tranw = .FALSE.                        ! To activate unconverged  outer iteration warning

! Thermal-hydraulics parameters
DOUBLE PRECISION :: pow                                        ! Reactor power for given geometry (watt)
DOUBLE PRECISION :: ppow                                       ! Reactor percent power in percent
DOUBLE PRECISION :: tpow                                       ! Total reactor power
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: npow            ! nodes power (watt)
DOUBLE PRECISION :: tin                                        ! coolant inlet temperature (kelvin)
DOUBLE PRECISION :: cflow                                      ! Sub-channel mass flow rate (kg/s)
DOUBLE PRECISION :: rf, tg, tc, ppitch                         ! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
DOUBLE PRECISION :: rg, rc                                     ! Outer radius of gap and cladding
DOUBLE PRECISION :: dia, dh, farea                             ! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
DOUBLE PRECISION :: cf                                         ! heat fraction deposited into coolant
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: node_nf       ! Number of fuel pin per node
INTEGER :: nm                                      ! number of Fuel meat mesh
INTEGER :: nt                                      ! Number Total mesh
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: tfm           ! Fuel pin mesh temperature for each nodes
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rdel            ! mesh delta
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rpos            ! mesh position
DOUBLE PRECISION :: th_err                                     ! Doppler error
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ent             ! Coolant Enthalpy (J/Kg)
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: heatf           ! Heat flux (W/m2
INTEGER, PARAMETER :: thunit = 300                 ! Unit number to open steam table file
DOUBLE PRECISION, PARAMETER :: pi = 3.14159265

! Steam Table data
INTEGER, PARAMETER:: ntem = 16   ! Number of temperature in steam table
DOUBLE PRECISION, DIMENSION(ntem,6) :: stab  ! Steam table matrix


END MODULE sdata
