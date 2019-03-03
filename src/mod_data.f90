MODULE sdata

SAVE

INTEGER, PARAMETER :: DP = KIND(1.E0)

CHARACTER(LEN=100) :: mode

INTEGER :: ng     ! number of groups
INTEGER :: nmat   ! number of materials
!! CXs Assigned to Nodes
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sigtr          ! Transport macroscopic cx
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: siga           ! Absorption macroscopic cx
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: nuf            ! nu* fission macroscopic cx
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sigf           ! fission macroscopic cx
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: chi            ! neutron fission spectrum
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: sigs         ! Scattering macroscopic cx
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: D              ! Diffusion coefficient
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sigr           ! Removal macroscopic cx

!! CXs Assigned to Materials
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xsigtr          ! Transport macroscopic cx
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xsiga           ! Absorption macroscopic cx
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xnuf            ! nu* fission macroscopic cx
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xsigf           ! fission macroscopic cx
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xD              ! Diffusion coefficient
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: xsigr           ! Removal macroscopic cx
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: xsigs         ! Scattering macroscopic cx
LOGICAL :: ccnuf = .TRUE.                            ! Logical variable to check the presence of fissile material
LOGICAL :: ccsigf = .TRUE.                           ! Logical variable to check the presence of fissile material

! Geometry
INTEGER :: nx, ny, nz                                ! Number of assemblies in x, y, and z directions
INTEGER :: nxx, nyy, nzz                             ! Number of nodes in x, y, and z directions
INTEGER :: nnod                                      ! Number of nodes
INTEGER, DIMENSION(:), ALLOCATABLE :: ix, iy, iz
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: xyz
INTEGER, DIMENSION(:), ALLOCATABLE :: xdiv, ydiv, zdiv     ! Assembly division
REAL(DP), DIMENSION(:), ALLOCATABLE :: xdel, ydel, zdel, vdel  ! Delta x, y and z and nodes' volume in cm3
INTEGER :: xwest, xeast, ysouth, ynorth, zbott, ztop       ! Boundary conditions
INTEGER, DIMENSION(:), ALLOCATABLE :: mat

! Keff, flux and currents
REAL(DP) :: Ke
TYPE :: NODE_DATA
    REAL(DP), DIMENSION(6) :: jo             ! Nodals' outgoing currents (X+,X-,Y+, Y-, Z+, Z-)
    REAL(DP), DIMENSION(6) :: ji             ! Nodals' ingoing currents  (X+,X-,Y+, Y-, Z+, Z-)
    REAL(DP), DIMENSION(3) :: L              ! Zeroth transverse leakages (Lx, Ly, Lz)
    REAL(DP), DIMENSION(7) :: Q              ! Nodal's source and source moments (0, x1, y1, z1, x2, y2, z2)
    REAL(DP), DIMENSION(6,6) :: P            ! Response matrix
    REAL(DP), DIMENSION(6,7) :: R            ! Response matrix
END TYPE
TYPE(NODE_DATA), DIMENSION(:,:), ALLOCATABLE :: nod

REAL(DP), DIMENSION(:,:), ALLOCATABLE :: f0, fx1, fy1, fz1, fx2, fy2, fz2      ! Flux and Flux moments
REAL(DP), DIMENSION(:), ALLOCATABLE :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2      ! Fission source moments
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: c0, cx1, cy1, cz1, cx2, cy2, cz2  ! neutron precusor density
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Parameters at previous time step

TYPE :: STAGGERED
    INTEGER :: smax, smin                             ! imax and imin along x and y direction for staggered nodes
END TYPE
TYPE(STAGGERED), DIMENSION(:), ALLOCATABLE :: ystag, xstag

! Extra Sources
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: exsrc

! Iteration Control
REAL(DP) :: ferc = 1.e-5    ! Flux Error Criteria
REAL(DP) :: serc = 1.e-5    ! Fission source Error CRITERIA
REAL(DP) :: ierc = 1.e-5    ! Inner Iteration Error Criteria
REAL(DP) :: fer, ser        ! Flux and Fission source error in BCSEARCH calcs.
INTEGER :: nin = 2      ! Maximum inner iteration
INTEGER :: nout = 500   ! Maximum outer iteration
INTEGER :: nac = 5      ! number of outer iteration before next source EXTRAPOLATION
INTEGER :: th_niter = 20                           ! Maximum number of thermal-hydraulics iteration

! OUTPUT PRINT OPTION
INTEGER :: aprad=1, apaxi=1, afrad=1

!ADF
TYPE :: ADF_TYPE
    REAL(DP), DIMENSION(6) :: dc
END TYPE
TYPE(ADF_TYPE), DIMENSION(:,:), ALLOCATABLE :: al

! FUEL TEMPERATURE
REAL(DP), DIMENSION(:), ALLOCATABLE :: ftem       ! Fuel temperature in Kelvin for each nodes
REAL(DP) :: rftem      ! Fuel temperature Reference in Kelvin
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: fsigtr, fsiga, fnuf, fsigf   ! CX changes per fuel temp changes
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: fsigs

! MODERATOR TEMPERATURE
REAL(DP), DIMENSION(:), ALLOCATABLE :: mtem       ! Moderator temperature in Kelvin for each nodes
REAL(DP) :: rmtem      ! Moderator temperature Reference in Kelvin
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: msigtr, msiga, mnuf, msigf   ! CX changes per Moderator temp changes
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: msigs

! COOLANT DENSITY
REAL(DP), DIMENSION(:), ALLOCATABLE :: cden       ! Coolant Density in g/cm3 for each nodes
REAL(DP) :: rcden      ! Coolant Density Reference in g/cm3
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: lsigtr, lsiga, lnuf, lsigf   ! CX changes per Coolant density changes
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: lsigs

! Crod changes
INTEGER :: nb                                                     ! Number of CR banks
REAL(DP), DIMENSION(:), ALLOCATABLE :: bpos  ! CR bank position
REAL(DP), DIMENSION(:), ALLOCATABLE :: fbpos    ! Final CR bank position
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: dsigtr, dsiga, dnuf, dsigf   ! CX incerement or decrement due to CR insertion
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: dsigs
REAL(DP), DIMENSION(:), ALLOCATABLE :: tmove    ! Time when CR bank starts moving
REAL(DP), DIMENSION(:), ALLOCATABLE :: bspeed   ! CR bank movement speed
INTEGER, DIMENSION(:), ALLOCATABLE :: mdir  ! To indicate CR movement direction (0=do not move, 1=down, 2 = up)
INTEGER :: cusp = 0                         ! Rod cusping option

! Boron Concentration
REAL(DP) :: bcon       ! Boron concentration in ppm
REAL(DP) :: rbcon      ! Boron concentration in ppm Reference
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: csigtr, csiga, cnuf, csigf   ! CX changes due to boron concentration
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: csigs                      ! Used only for CBCS card

! Transient parameters
INTEGER, PARAMETER :: nf = 6                       ! Number of delaye dneutron precusor family
REAL(DP), DIMENSION(nf) :: ibeta, lamb                 ! beta (delayed neutron fraction) and precusor decay constant
REAL(DP) :: tbeta                                      ! total beta
REAL(DP), DIMENSION(:), ALLOCATABLE :: velo            ! Neutron velocity
REAL(DP) :: ttot                                       ! TOTAL SIMULATION TIME
REAL(DP) :: tstep1                                     ! FIRST TIME STEP
REAL(DP) :: tstep2                                     ! SECOND TIME STEP
REAL(DP) :: tdiv                                       ! WHEN SECOND TIME STEP APPLY
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: omeg          ! Exponential transformation constant
LOGICAL :: tranw = .FALSE.                        ! To activate unconverged  outer iteration warning

! Thermal-hydraulics parameters
REAL(DP) :: pow                                        ! Reactor power for given geometry (watt)
REAL(DP) :: ppow                                       ! Reactor percent power in percent
REAL(DP) :: tpow                                       ! Total reactor power
REAL(DP), DIMENSION(:), ALLOCATABLE :: npow            ! nodes power (watt)
REAL(DP) :: tin                                        ! coolant inlet temperature (kelvin)
REAL(DP) :: cflow                                      ! Sub-channel mass flow rate (kg/s)
REAL(DP) :: rf, tg, tc, ppitch                         ! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
REAL(DP) :: rg, rc                                     ! Outer radius of gap and cladding
REAL(DP) :: dia, dh, farea                             ! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
REAL(DP) :: cf                                         ! heat fraction deposited into coolant
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: node_nf       ! Number of fuel pin per node
INTEGER :: nm                                      ! number of Fuel meat mesh
INTEGER :: nt                                      ! Number Total mesh
REAL(DP), DIMENSION(:,:), ALLOCATABLE :: tfm           ! Fuel pin mesh temperature for each nodes
REAL(DP), DIMENSION(:), ALLOCATABLE :: rdel            ! mesh delta
REAL(DP), DIMENSION(:), ALLOCATABLE :: rpos            ! mesh position
REAL(DP) :: th_err                                     ! Doppler error
REAL(DP), DIMENSION(:), ALLOCATABLE :: ent             ! Coolant Enthalpy (J/Kg)
REAL(DP), DIMENSION(:), ALLOCATABLE :: heatf           ! Heat flux (W/m2
INTEGER, PARAMETER :: thunit = 300                 ! Unit number to open steam table file
REAL(DP), PARAMETER :: pi = 3.14159265

! Steam Table data
INTEGER, PARAMETER:: ntem = 16   ! Number of temperature in steam table
REAL(DP), DIMENSION(ntem,6) :: stab  ! Steam table matrix


END MODULE sdata
