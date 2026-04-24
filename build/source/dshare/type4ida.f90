module type4ida

! data types
USE nrtype
USE, intrinsic :: iso_c_binding

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    zLookup,      & ! data vector with variable length dimension (rkind)
                    model_options   ! defines the model decisions
use device_data_types
implicit none

type data4ida
  type(c_ptr)                     :: ida_mem                         ! IDA memory
  real(rkind)                     :: dt                              ! data step
  integer(i4b),device,pointer                    :: nSnow_p(:)                           ! number of snow layers
  integer(i4b)                    :: nSoil                           ! number of soil layers
  integer(i4b),device,pointer                    :: nLayers_p(:)                         ! total number of layers
  integer(i4b) :: nGRU
  integer(i4b)                    :: nState                          ! total number of state variables
  integer(i4b)                    :: ixMatrix                        ! form of matrix (dense or banded)
  logical(lgt)                    :: firstSubStep                    ! flag to indicate if we are processing the first sub-step
  logical(lgt)                    :: firstFluxCall                   ! flag to indicate if we are processing the first flux call
  logical(lgt)                    :: firstSplitOper                  ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt),device,pointer                    :: computeVegFlux_p(:)                  ! flag to indicate if computing fluxes over vegetation
  logical(lgt)                    :: scalarSolution                  ! flag to denote if implementing the scalar solution
  type(model_options),allocatable :: model_decisions(:)              ! model decisions
  type(decisions_device),pointer :: decisions_p
  type(zLookup_device),pointer                   :: lookup_data_p                  ! lookup tables
  type(type_data_device),pointer                     :: type_data_p                       ! type of vegetation and soil
  type(attr_data_device),pointer                     :: attr_data_p                       ! spatial attributes
  type(mpar_data_device),pointer               :: mpar_data_p                       ! model parameters
  type(veg_parameters),pointer :: veg_param_p
  type(forc_data_device),pointer                     :: forc_data_p                       ! model forcing data
  type(bvar_data_device),pointer               :: bvar_data_p                       ! model variables for the local basin
  type(prog_data_device),pointer               :: prog_data_p                       ! prognostic variables for a local HRU
  type(indx_data_device),pointer               :: indx_data_p                       ! indices defining model states and layers
  type(diag_data_device),pointer               :: diag_data_p                       ! diagnostic variables for a local HRU
  type(flux_data_device),pointer               :: flux_data_p                       ! model fluxes for a local HRU
  type(deriv_data_device),pointer               :: deriv_data_p                      ! derivatives in model fluxes w.r.t. relevant state variables
  real(qp), device,pointer           :: sMul_p(:,:)                         ! state vector multiplier (used in the residual calculations)
  real(rkind), device,pointer        :: dMat_p(:,:)                         ! diagonal of the Jacobian matrix
  real(rkind), allocatable,device        :: fluxVec(:,:)                      ! flux vector
  real(qp), allocatable,device           :: resVec(:,:)                       ! residual vector
  real(qp), allocatable,device           :: resSink(:,:)                      ! additional (sink) terms on the RHS of the state equation
  real(rkind),device,allocatable        :: atol(:,:)                         ! vector of absolute tolerances
  real(rkind),device,allocatable        :: rtol(:,:)                         ! vector of relative tolerances
  integer(i4b),device,allocatable                    :: ixSaturation(:)                    ! index of the lowest saturated layer
  real(rkind), device,allocatable        :: dBaseflow_dMatric(:,:,:)          ! derivative in baseflow w.r.t. matric head (s-1)
  integer(i4b)                    :: err                             ! error code
  character(len=50)               :: message                         ! error message
  real(rkind),device,allocatable                     :: scalarCanopyTempPrev(:)            ! previous value for temperature of the vegetation canopy (K)
  real(rkind),device, allocatable        :: mLayerTempPrev(:,:)               ! previous vector of layer temperature (K)
  real(rkind),device, allocatable        :: mLayerMatricHeadPrev(:,:)         ! previous value for total water matric potential (m)
  real(rkind),device,allocatable                     :: scalarCanopyEnthalpyTrial(:)       ! trial value for enthalpy of the vegetation canopy (J m-2)
  real(rkind),device,allocatable                     :: scalarCanopyTempTrial(:)           ! trial value for temperature of the vegetation canopy (K)
  real(rkind),device,allocatable                     :: scalarCanopyWatTrial(:)            ! trial value for mass of total water on the vegetation canopy (kg m-2)
  real(rkind),device, allocatable        :: mLayerTempTrial(:,:)              ! trial vector of layer temperature (K)
  real(rkind),device, allocatable        :: mLayerMatricHeadTrial(:,:)        ! trial value for total water matric potential (m)
  real(rkind),device,allocatable                     :: scalarCanopyTempPrime(:)           ! prime value for temperature of the vegetation canopy (K s-1)
  real(rkind),device,allocatable                     :: scalarCanopyWatPrime(:)            ! prime value for mass of total water on the vegetation canopy (kg m-2 s-1)
  real(rkind),device, allocatable        :: mLayerTempPrime(:,:)              ! prime vector of temperature of each snow and soil layer (K s-1)
  real(rkind),device, allocatable        :: mLayerMatricHeadPrime(:,:)        ! prime vector of matric head of each snow and soil layer (m s-1)
  real(rkind),device, allocatable        :: mLayerVolFracWatPrime(:,:)        ! prime vector of volumetric total water content of each snow and soil layer (s-1)
 end type data4ida


end module type4ida





