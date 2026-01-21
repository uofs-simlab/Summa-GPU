! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE data_types
 ! used to define model data structures
 USE nrtype, integerMissing=>nr_integerMissing
 USE var_lookup,only:maxvarFreq
 USE var_lookup,only:maxvarStat
 USE var_lookup,only:maxvarDecisions  ! maximum number of decisions
 USE var_lookup,only:iLookPARAM       ! lookup indices for model parameters
 USE var_lookup,only:iLookFLUX        ! lookup indices for flux data
 USE var_lookup,only:iLookDERIV       ! lookup indices for derivative data
 USE var_lookup,only:iLookFORCE       ! lookup indices for forcing data 
 USE var_lookup,only:iLookDIAG        ! lookup indices for diagnostic variable data
 USE var_lookup,only:iLookDECISIONS   ! lookup indices for elements of the decision structure
 USE var_lookup,only:iLookPROG        ! lookup indices for prognostic variables
 implicit none
 ! constants necessary for variable defs
 private

 ! ***********************************************************************************************************
 ! Define the model decisions
 ! ***********************************************************************************************************
 ! the model decision structure
 type,public  :: model_options
  character(len=64)                      :: cOption   = 'notPopulatedYet'
  character(len=64)                      :: cDecision = 'notPopulatedYet'
  integer(i4b)                           :: iDecision = integerMissing
 end type model_options

 ! ***********************************************************************************************************
 ! Define metadata for model forcing datafile
 ! ***********************************************************************************************************
 ! define a derived type for the data in the file
 type,public  :: file_info
  character(len=256)                     :: filenmData='notPopulatedYet'  ! name of data file
  integer(i4b)                           :: nVars                         ! number of variables in the file
  integer(i4b)                           :: nTimeSteps                    ! number of variables in the file
  integer(i4b),allocatable               :: var_ix(:)                     ! index of each forcing data variable in the data structure
  integer(i4b),allocatable               :: data_id(:)                    ! netcdf variable id for each forcing data variable
  character(len=256),allocatable         :: varName(:)                    ! netcdf variable name for each forcing data variable
  real(rkind)                            :: firstJulDay                   ! first julian day in forcing file
  real(rkind)                            :: convTime2Days                 ! factor to convert time to days
 end type file_info

 ! ***********************************************************************************************************
 ! Define metadata on model parameters
 ! ***********************************************************************************************************
 ! define a data type to store model parameter information
 type,public  :: par_info
  real(rkind)                            :: default_val                   ! default parameter value
  real(rkind)                            :: lower_limit                   ! lower bound
  real(rkind)                            :: upper_limit                   ! upper bound
 endtype par_info

 ! ***********************************************************************************************************
 ! Define variable metadata
 ! ***********************************************************************************************************
 ! define derived type for model variables, including name, description, and units
 type,public :: var_info
  character(len=64)                      :: varname   = 'empty'           ! variable name
  character(len=128)                     :: vardesc   = 'empty'           ! variable description
  character(len=64)                      :: varunit   = 'empty'           ! variable units
  integer(i4b)                           :: vartype   = integerMissing    ! variable type
  integer(i4b),dimension(maxvarFreq)     :: ncVarID   = integerMissing    ! netcdf variable id (missing if frequency is not desired)
  integer(i4b),dimension(maxvarFreq)     :: statIndex = integerMissing    ! index of desired statistic for temporal aggregation
  logical(lgt)                           :: varDesire = .false.           ! flag to denote if the variable is desired for model output
 endtype var_info

 ! define extended data type (include indices to map onto parent data type)
 type,extends(var_info),public :: extended_info
  integer(i4b)                           :: ixParent                      ! index in the parent data structure
 endtype extended_info

 ! define extended data type (includes named variables for the states affected by each flux)
 type,extends(var_info),public :: flux2state
  integer(i4b)                           :: state1                        ! named variable of the 1st state affected by the flux
  integer(i4b)                           :: state2                        ! named variable of the 2nd state affected by the flux
 endtype flux2state

 ! ***********************************************************************************************************
 ! Define summary of data structures
 ! ***********************************************************************************************************
 ! data structure information
 type,public :: struct_info
  character(len=32)                      :: structName                    ! name of the data structure
  character(len=32)                      :: lookName                      ! name of the look-up variables
  integer(i4b)                           :: nVar                          ! number of variables in each data structure
 end type struct_info

 ! ***********************************************************************************************************
 ! Define data types to map between GRUs and HRUs
 ! ***********************************************************************************************************

 ! hru info data structure
 type, public :: hru_info
  integer(i4b)                           :: hru_nc                        ! index of the hru in the netcdf file
  integer(i4b)                           :: hru_ix                        ! index of the hru in the run domain
  integer(i8b)                           :: hru_id                        ! id (non-sequential number) of the hru
  integer(i4b)                           :: nSnow                         ! number of snow layers
  integer(i4b)                           :: nSoil                         ! number of soil layers
 endtype hru_info

 ! define mapping from GRUs to the HRUs
 type, public :: gru2hru_map
  integer(i8b)                           :: gru_id                        ! id of the gru
  integer(i4b)                           :: hruCount                      ! total number of hrus in the gru
  type(hru_info), allocatable            :: hruInfo(:)                    ! basic information of HRUs within the gru
  integer(i4b)                           :: gru_nc                        ! index of gru in the netcdf file
 endtype gru2hru_map

 ! define the mapping from the HRUs to the GRUs
 type, public :: hru2gru_map
  integer(i4b)                           :: gru_ix                        ! index of gru which the hru belongs to
  integer(i4b)                           :: localHRU_ix                   ! index of a hru within a gru (start from 1 per gru)
 endtype hru2gru_map

 ! ***********************************************************************************************************
 ! Define hierarchal derived data types
 ! ***********************************************************************************************************
  ! define derived types to hold look-up tables for each soil layer
 ! ** double precision type
 type, public :: dLookup
  real(rkind),allocatable                :: lookup(:)                     ! lookup(:)
 endtype dLookup
 ! ** double precision type for a variable number of soil layers; variable length
 type, public :: vLookup
  type(dLookup),allocatable              :: var(:)                        ! var(:)%lookup(:)
 endtype vLookup
 ! ** double precision type for a variable number of soil layers
 type, public :: zLookup
  real(rkind),allocatable :: temperature(:,:)
  real(rkind),allocatable :: psiLiq_int(:,:)
  real(rkind),allocatable :: deriv2(:,:)
  ! type(vLookup),allocatable              :: z(:)                          ! z(:)%var(:)%lookup(:)
 endtype zLookup
 ! ** double precision type for a variable number of soil layers
 type, public :: hru_z_vLookup
  type(zLookup),allocatable              :: hru(:)                        ! hru(:)%z(:)%var(:)%lookup(:)
 endtype hru_z_vLookup
 ! ** double precision type for a variable number of soil layers
 type, public :: gru_hru_z_vLookup
  type(hru_z_vLookup),allocatable        :: gru(:)                        ! gru(:)%hru(:)%z(:)%var(:)%lookup(:)
 endtype gru_hru_z_vLookup
 ! define derived types to hold multivariate data for a single variable (different variables have different length)
 ! NOTE: use derived types here to facilitate adding the "variable" dimension
 ! ** double precision type
 type, public :: dlength
  real(rkind),allocatable                :: dat(:)                        ! dat(:)
 endtype dlength
 ! ** integer type (4 byte)
 type, public :: ilength
  integer(i4b),allocatable               :: dat(:)                        ! dat(:)
 endtype ilength
 ! ** integer type (8 byte)
 type, public :: i8length
  integer(i8b),allocatable               :: dat(:)                        ! dat(:)
 endtype i8length
 ! ** logical type
 type, public :: flagVec
  logical(lgt),allocatable               :: dat(:)                        ! dat(:)
 endtype flagVec

 ! define derived types to hold data for multiple variables
 ! NOTE: use derived types here to facilitate adding extra dimensions (e.g., spatial)

 ! ** double precision type of variable length
 type, public :: var_dlength
  type(dlength),allocatable              :: var(:)                        ! var(:)%dat
 endtype var_dlength
 ! ** integer type of variable length (4 byte)
 type, public :: var_ilength
  type(ilength),allocatable              :: var(:)                        ! var(:)%dat
 endtype var_ilength
 ! ** integer type of variable length (8 byte)
 type, public :: var_i8length
  type(i8length),allocatable             :: var(:)                        ! var(:)%dat
 endtype var_i8length
 ! ** logical type of variable length
 type, public :: var_flagVec
  type(flagVec),allocatable              :: var(:)                        ! var(:)%dat
 endtype var_flagVec

 ! ** double precision type of fixed length
 type, public :: var_d
  real(rkind),allocatable                :: var(:)                        ! var(:)
 endtype var_d
 ! ** integer type of fixed length (4 byte)
 type, public :: var_i
  integer(i4b),allocatable               :: var(:)                        ! var(:)
 endtype var_i
 ! ** integer type of fixed length (8 byte)
 type, public :: var_i8
  integer(i8b),allocatable               :: var(:)                        ! var(:)
 endtype var_i8

 ! ** double precision type of fixed length
 type, public :: hru_d
  real(rkind),allocatable                :: hru(:)                        ! hru(:)
 endtype hru_d
 ! ** integer type of fixed length (4 byte)
 type, public :: hru_i
  integer(i4b),allocatable               :: hru(:)                        ! hru(:)
 endtype hru_i
 ! ** integer type of fixed length (8 byte)
 type, public :: hru_i8
  integer(i8b),allocatable               :: hru(:)                        ! hru(:)
 endtype hru_i8

 ! define derived types to hold JUST the HRU dimension
 ! ** double precision type of variable length
 type, public :: hru_doubleVec
  type(var_dlength),allocatable          :: hru(:)                        ! hru(:)%var(:)%dat
 endtype hru_doubleVec
 ! ** integer type of variable length (4 byte)
 type, public :: hru_intVec
  type(var_ilength),allocatable          :: hru(:)                        ! hru(:)%var(:)%dat
 endtype hru_intVec
 ! ** integer type of variable length (8 byte)
 type, public :: hru_int8Vec
  type(var_i8length),allocatable         :: hru(:)                        ! hru(:)%var(:)%dat
 endtype hru_int8Vec
 ! ** double precision type of fixed length
 type, public :: hru_double
  type(var_d),allocatable                :: hru(:)                        ! hru(:)%var(:)
 endtype hru_double
 ! ** integer type of fixed length (4 byte)
 type, public :: hru_int
  type(var_i),allocatable                :: hru(:)                        ! hru(:)%var(:)
 endtype hru_int
 ! ** integer type of fixed length (8 byte)
 type, public :: hru_int8
  type(var_i8),allocatable               :: hru(:)                        ! hru(:)%var(:)
 endtype hru_int8

 ! define derived types to hold JUST the HRU dimension
 ! ** double precision type of variable length
 type, public :: gru_doubleVec
  type(var_dlength),allocatable          :: gru(:)                        ! gru(:)%var(:)%dat
 endtype gru_doubleVec
 ! ** integer type of variable length (4 byte)
 type, public :: gru_intVec
  type(var_ilength),allocatable          :: gru(:)                        ! gru(:)%var(:)%dat
 endtype gru_intVec
 ! ** integer type of variable length (8 byte)
 type, public :: gru_int8Vec
  type(var_i8length),allocatable         :: gru(:)                        ! gru(:)%var(:)%dat
 endtype gru_int8Vec
 ! ** double precision type of fixed length
 type, public :: gru_double
  type(var_d),allocatable                :: gru(:)                        ! gru(:)%var(:)
 endtype gru_double
 ! ** integer type of variable length (4 byte)
 type, public :: gru_int
  type(var_i),allocatable                :: gru(:)                        ! gru(:)%var(:)
 endtype gru_int
 ! ** integer type of variable length (8 byte)
 type, public :: gru_int8
  type(var_i8),allocatable               :: gru(:)                        ! gru(:)%var(:)
 endtype gru_int8

 ! define derived types to hold BOTH the GRU and HRU dimension
 ! ** double precision type of variable length
 type, public :: gru_hru_doubleVec
  type(hru_doubleVec),allocatable        :: gru(:)                        ! gru(:)%hru(:)%var(:)%dat
 endtype gru_hru_doubleVec
 ! ** integer type of variable length (4 byte)
 type, public :: gru_hru_intVec
  type(hru_intVec),allocatable           :: gru(:)                        ! gru(:)%hru(:)%var(:)%dat
 endtype gru_hru_intVec
 ! ** integer type of variable length (8 byte)
 type, public :: gru_hru_int8Vec
  type(hru_int8Vec),allocatable          :: gru(:)                        ! gru(:)%hru(:)%var(:)%dat
 endtype gru_hru_int8Vec
 ! ** double precision type of fixed length
 type, public :: gru_hru_double
  type(hru_double),allocatable           :: gru(:)                        ! gru(:)%hru(:)%var(:)
 endtype gru_hru_double
 ! ** integer type of variable length (4 byte)
 type, public :: gru_hru_int
  type(hru_int),allocatable              :: gru(:)                        ! gru(:)%hru(:)%var(:)
 endtype gru_hru_int
 ! ** integer type of variable length (8 byte)
 type, public :: gru_hru_int8
  type(hru_int8),allocatable             :: gru(:)                        ! gru(:)%hru(:)%var(:)
 endtype gru_hru_int8
 ! ** double precision type of fixed length
 type, public :: gru_d
  type(hru_d),allocatable                :: gru(:)                        ! gru(:)%hru(:)
 endtype gru_d
 ! ** integer type of fixed length
 type, public :: gru_i
  type(hru_i),allocatable                :: gru(:)                        ! gru(:)%hru(:)
 endtype gru_i

 integer(i4b),parameter :: len_msg=256 ! length of character string used in class definitions

 ! ***********************************************************************************************************
 ! Define classes used to simplify calls to the subrotuines in computFlux
 ! ***********************************************************************************************************
 ! Note: class procedures are located in the contains block of this (data_types) module
 ! ** vegNrgFlux
 type, public :: out_type_vegNrgFlux ! class for intent(out) arguments in vegNrgFlux call
   
   integer(i4b)             :: err                                     ! intent(out): error code
   character(len=len_msg)   :: cmessage                                ! intent(out): error message
  contains
   procedure :: finalize => finalize_out_vegNrgFlux
 end type out_type_vegNrgFlux
 ! ** end vegNrgFlux

 ! ** ssdNrgFlux
 type, public :: out_type_ssdNrgFlux ! class for intent(inout) arguments in ssdNrgFlux call
   integer(i4b)             :: err                               ! intent(out): error code
   character(len=len_msg)   :: cmessage                          ! intent(out): error message
  contains
   procedure :: finalize   => finalize_out_ssdNrgFlux
 end type out_type_ssdNrgFlux
 ! ** end ssdNrgFlux

 ! ** vegLiqFlux
 type, public :: out_type_vegLiqFlux ! class for intent(out) arguments in vegLiqFlux call
   integer(i4b)             :: err                               ! intent(out): error code
   character(len=len_msg)   :: cmessage                          ! intent(out): error message
  contains
   procedure :: finalize => finalize_out_vegLiqFlux
 end type out_type_vegLiqFlux
 ! ** end vegLiqFlux

 ! ** snowLiqFlx
 type, public :: out_type_snowLiqFlx ! class for intent(out) arguments in snowLiqFlx call
   integer(i4b)             :: err                               ! intent(out):   error code
   character(len=len_msg)   :: cmessage                          ! intent(out):   error message
  contains
   procedure :: finalize => finalize_out_snowLiqFlx
 end type out_type_snowLiqFlx
 ! ** end snowLiqFlx

 ! ** soilLiqFlx
 type, public :: out_type_soilLiqFlx ! class for intent(out) arguments in soilLiqFlx call
   integer(i4b)             :: err                               ! intent(out):   error code
   character(len=len_msg)   :: cmessage                          ! intent(out):   error message
  contains
   procedure :: finalize   => finalize_out_soilLiqFlx
 end type out_type_soilLiqFlx
 ! ** end soilLiqFlx

 ! ** groundwatr
 type, public :: out_type_groundwatr ! class for intent(out) arguments in groundwatr call
   integer(i4b)             :: err                               ! intent(out):   error code
   character(len=len_msg)   :: cmessage                          ! intent(out):   error message
  contains
   procedure :: finalize => finalize_out_groundwatr
 end type out_type_groundwatr
 ! ** end groundwatr

 ! ** bigAquifer
 type, public :: out_type_bigAquifer  ! class for intent(out) arguments in bigAquifer call
   integer(i4b)             :: err                               ! intent(out):   error code
   character(len=len_msg)   :: cmessage                          ! intent(out):   error message
  contains
   procedure :: finalize   => finalize_out_bigAquifer
 end type out_type_bigAquifer
 ! ** end bigAquifer

 ! ***********************************************************************************************************
 ! Define classes used to simplify calls to the subrotuines in soilLiqFlx
 ! ***********************************************************************************************************

 ! ** diagv_node
 type, public :: out_type_diagv_node ! intent(out) data
   ! output: error control
   integer(i4b)           :: err                       ! error code
   character(len=len_msg) :: message                   ! error message
  contains
   procedure :: finalize   => finalize_out_diagv_node
 end type out_type_diagv_node
 ! ** end diagv_node

 ! ** surfaceFlx
 type, public :: out_type_surfaceFlx ! intent(out) data
   ! output: error control
   integer(i4b)            :: err     ! error code
   character(len=len_msg)  :: message ! error message
  contains
   procedure :: finalize   => finalize_out_surfaceFlx
 end type out_type_surfaceFlx 
 ! ** end surfaceFlx

 ! ** iLayerFlux
 type, public :: out_type_iLayerFlux ! intent(out) data
   ! output: error control
   integer(i4b)           :: err     ! error code
   character(len=len_msg) :: message ! error message
  contains
   procedure :: finalize => finalize_out_iLayerFlux
 end type out_type_iLayerFlux
 ! ** end iLayerFlux

 ! ** qDrainFlux
 type, public :: out_type_qDrainFlux ! intent(out) data
   ! output: error control
   integer(i4b)           :: err     ! error code
   character(len=len_msg) :: message ! error message
  contains
   procedure :: finalize => finalize_out_qDrainFlux
 end type out_type_qDrainFlux
 ! ** end qDrainFlux

 ! ***********************************************************************************************************
 ! Define classes used to simplify calls to the subrotuines in opSplittin
 ! ***********************************************************************************************************
 ! ** stateFilter
 type, public :: out_type_stateFilter ! class for intent(out) arguments in stateFilter call
   integer(i4b)             :: err                         ! intent(out): error code
   character(len=len_msg)   :: cmessage                    ! intent(out): error message
  contains
   procedure :: finalize   => finalize_out_stateFilter
 end type out_type_stateFilter
 ! ** end stateFilter

 ! ** indexSplit
 type, public :: in_type_indexSplit  ! class for intent(in) arguments in indexSplit call
  !  integer(i4b)             :: nSnow                       ! intent(in): number of snow layers
  !  integer(i4b)             :: nSoil                       ! intent(in): number of soil layers
  !  integer(i4b)             :: nLayers                     ! intent(in): total number of layers
   integer(i4b)             :: nSubset                     ! intent(in): number of states in the subset
  contains
   procedure :: initialize => initialize_in_indexSplit
 end type in_type_indexSplit

 type, public :: out_type_indexSplit ! class for intent(out) arguments in indexSplit call
   integer(i4b)             :: err                         ! intent(out): error code
   character(len=len_msg)   :: cmessage                    ! intent(out): error message
  contains
   procedure :: finalize   => finalize_out_indexSplit
 end type out_type_indexSplit
 ! ** end indexSplit

 ! ** varSubstep
 type, public :: in_type_varSubstep  ! class for intent(in) arguments in varSubstep call
   real(rkind)              :: dt                          ! intent(in): time step (s)
   real(rkind)              :: dtInit                      ! intent(in): initial time step (seconds)
   real(rkind)              :: dt_min                      ! intent(in): minimum time step (seconds)
   real(rkind)              :: whole_step                  ! intent(in): length of whole step for surface drainage and average flux
   integer(i4b)             :: nSubset                     ! intent(in): total number of variables in the state subset
   logical(lgt)             :: doAdjustTemp                ! intent(in): flag to indicate if we adjust the temperature
   logical(lgt)             :: firstSubStep                ! intent(in): flag to denote first sub-step
  !  logical(lgt)             :: computeVegFlux              ! intent(in): flag to denote if computing energy flux over vegetation
   logical(lgt)             :: scalarSolution              ! intent(in): flag to denote computing the scalar solution
   integer(i4b)             :: iStateSplit                 ! intent(in): index of the layer in the splitting operation
  !  type(var_flagVec)        :: fluxMask                    ! intent(in): mask for the fluxes used in this given state subset
  contains
   procedure :: initialize => initialize_in_varSubstep
 end type in_type_varSubstep

 type, public :: io_type_varSubstep  ! class for intent(inout) arguments in varSubstep call
   logical(lgt)             :: firstFluxCall               ! intent(inout): flag to indicate if we are processing the first flux call
  !  type(var_ilength)        :: fluxCount                   ! intent(inout): number of times fluxes are updated (should equal nsubstep)
  !  integer(i4b)             :: ixSaturation                ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
  contains
   procedure :: initialize => initialize_io_varSubstep
   procedure :: finalize   => finalize_io_varSubstep
 end type io_type_varSubstep

 type, public :: out_type_varSubstep  ! class for intent(out) arguments in varSubstep call
   real(rkind)              :: dtMultiplier                ! intent(out): substep multiplier (-)
   integer(i4b)             :: nSubsteps                   ! intent(out): number of substeps taken for a given split
   logical(lgt)             :: failedMinimumStep           ! intent(out): flag for failed substeps
   logical(lgt)             :: reduceCoupledStep           ! intent(out): flag to reduce the length of the coupled step
   logical(lgt)             :: tooMuchMelt                 ! intent(out): flag to denote that ice is insufficient to support melt
   integer(i4b)             :: err                         ! intent(out): error code
   character(len=len_msg)   :: cmessage                    ! intent(out): error message
  contains
   procedure :: finalize   => finalize_out_varSubstep
 end type out_type_varSubstep
 ! ** end varSubstep

 ! ***********************************************************************************************************
 ! Define classes used to simplify calls to the subrotuines in summaSolve4homegrown
 ! ***********************************************************************************************************

 type, public :: in_type_computJacob  ! class for intent(in) arguments in computJacob call
   ! input: model control
   real(rkind)              :: dt                          ! intent(in): length of the time step (seconds)
   integer(i4b)             :: nSnow                       ! intent(in): number of snow layers
   integer(i4b)             :: nSoil                       ! intent(in): number of soil layers
   integer(i4b)             :: nLayers                     ! intent(in): total number of layers in the snow+soil domain
   logical(lgt)             :: computeVegFlux              ! intent(in): flag to indicate if computing fluxes over vegetation
   logical(lgt)             :: computeBaseflow             ! intent(in): flag to indicate if computing baseflow
   integer(i4b)             :: ixMatrix                    ! intent(in): form of the Jacobian matrix
  contains
   procedure :: initialize => initialize_in_computJacob
 end type in_type_computJacob

 type, public :: out_type_computJacob  ! class for intent(out) arguments in computJacob call
   ! output: error control
   integer(i4b)             :: err                         ! intent(out): error code
   character(len=len_msg)   :: cmessage                    ! intent(out): error message
  contains
   procedure :: finalize => finalize_out_computJacob
 end type out_type_computJacob

 type, public :: in_type_lineSearchRefinement  ! class for intent(in) arguments in lineSearchRefinement call
   logical(lgt)             :: doSearch                    ! intent(in): flag to do the line search
   real(rkind)              :: fOld                        ! intent(in): old function value
  contains
   procedure :: initialize => initialize_in_lineSearchRefinement
 end type in_type_lineSearchRefinement

 type, public :: out_type_lineSearchRefinement  ! class for intent(out) arguments in lineSearchRefinement call
   real(rkind)              :: fNew                        ! intent(out): new function evaluation
   logical(lgt)             :: converged                   ! intent(out): convergence flag
   ! output: error control
   integer(i4b)             :: err                         ! intent(out): error code
   character(len=len_msg)   :: message                     ! intent(out): error message
  contains
   procedure :: finalize => finalize_out_lineSearchRefinement
 end type out_type_lineSearchRefinement

 ! ***********************************************************************************************************
 ! Define classes used to simplify calls to the subrotuines in systemSolv
 ! ***********************************************************************************************************

 type, public :: in_type_summaSolve4homegrown  ! class for intent(in) arguments in summaSolve4homegrown call
   real(rkind)              :: dt_cur                   ! intent(in): current stepsize
   real(rkind)              :: dt                       ! intent(in): entire time step for drainage pond rate
   integer(i4b)             :: iter                     ! intent(in): iteration index
   integer(i4b)             :: nSnow                    ! intent(in): number of snow layers
   integer(i4b)             :: nSoil                    ! intent(in): number of soil layers
   integer(i4b)             :: nLayers                  ! intent(in): total number of layers
   integer(i4b)             :: nLeadDim                 ! intent(in): length of the leading dimension of the Jacobian matrix (nBands or nState)
   integer(i4b)             :: nState                   ! intent(in): total number of state variables
   integer(i4b)             :: ixMatrix                 ! intent(in): type of matrix (full or band diagonal)
   logical(lgt)             :: firstSubStep             ! intent(in): flag to indicate if we are processing the first sub-step
   logical(lgt)             :: computeVegFlux           ! intent(in): flag to indicate if computing fluxes over vegetation
   logical(lgt)             :: scalarSolution           ! intent(in): flag to denote if implementing the scalar solution
   real(rkind)              :: fOld                     ! intent(in): old function evaluation
  contains
   procedure :: initialize => initialize_in_summaSolve4homegrown
 end type in_type_summaSolve4homegrown

 type, public :: io_type_summaSolve4homegrown  ! class for intent(inout) arguments in summaSolve4homegrown call
   logical(lgt)             :: firstFluxCall            ! intent(inout): flag to indicate if we are processing the first flux call
   real(rkind)              :: xMin,xMax                ! intent(inout): brackets of the root
   integer(i4b)             :: ixSaturation             ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
  contains
   procedure :: initialize => initialize_io_summaSolve4homegrown
   procedure :: finalize   => finalize_io_summaSolve4homegrown
 end type io_type_summaSolve4homegrown

 type, public :: out_type_summaSolve4homegrown  ! class for intent(out) arguments in summaSolve4homegrown call
   real(rkind)              :: fNew                     ! intent(out): new function evaluation
   logical(lgt)             :: converged                ! intent(out): convergence flag
   integer(i4b)             :: err                      ! intent(out): error code
   character(len=len_msg)   :: message                  ! intent(out): error message
  contains
   procedure :: finalize => finalize_out_summaSolve4homegrown
 end type out_type_summaSolve4homegrown

contains
 
 ! **** vegNrgFlux ****
 subroutine finalize_out_vegNrgFlux(out_vegNrgFlux,err,cmessage)
  class(out_type_vegNrgFlux),intent(in) :: out_vegNrgFlux              ! class object for intent(out) vegNrgFlux arguments
  integer(i4b),intent(out)              :: err                         ! error code
  character(*),intent(out)              :: cmessage                    ! error message from ssdNrgFlux

   ! intent(out) arguments: error control
   err                        =out_vegNrgFlux % err                   ! intent(out): error code
   cmessage                   =out_vegNrgFlux % cmessage              ! intent(out): error message
 end subroutine finalize_out_vegNrgFlux
 ! **** end vegNrgFlux ****

 ! **** ssdNrgFlux ****
 subroutine finalize_out_ssdNrgFlux(out_ssdNrgFlux,err,cmessage)
  class(out_type_ssdNrgFlux),intent(in) :: out_ssdNrgFlux              ! class object for intent(out) ssdNrgFlux arguments
  integer(i4b),intent(out)              :: err                         ! error code
  character(*),intent(out)              :: cmessage                    ! error message from ssdNrgFlux
   err                =out_ssdNrgFlux % err                               ! intent(out): error code
   cmessage           =out_ssdNrgFlux % cmessage                          ! intent(out): error message
 end subroutine finalize_out_ssdNrgFlux
 ! **** end ssdNrgFlux ****
 
 ! **** vegLiqFlux ****
 subroutine finalize_out_vegLiqFlux(out_vegLiqFlux,err,cmessage)
  class(out_type_vegLiqFlux),intent(in)   :: out_vegLiqFlux              ! class object for intent(out) vegLiqFlux arguments
  integer(i4b),intent(out)                :: err                         ! error code
  character(*),intent(out)                :: cmessage                    ! error message from vegLiqFlux
   err                         =out_vegLiqFlux % err                         ! intent(out): error code
   cmessage                    =out_vegLiqFlux % cmessage                    ! intent(out): error control
 end subroutine finalize_out_vegLiqFlux
 ! **** end vegLiqFlux ****

 ! **** snowLiqFlx ****
 subroutine finalize_out_snowLiqFlx(out_snowLiqFlx,err,cmessage)
  class(out_type_snowLiqFlx),intent(in)   :: out_snowLiqFlx              ! class object for intent(out) snowLiqFlx arguments
  integer(i4b),intent(out)                :: err                         ! error code
  character(*),intent(out)                :: cmessage                    ! error message from snowLiqFlx
  ! intent(out) arguments
  err     =out_snowLiqFlx % err                                          ! intent(out):   error code
  cmessage=out_snowLiqFlx % cmessage                                     ! intent(out):   error message
 end subroutine finalize_out_snowLiqFlx
 ! **** end snowLiqFlx ****

 ! **** soilLiqFlx ****
 subroutine finalize_out_soilLiqFlx(out_soilLiqFlx,err,cmessage)
  class(out_type_soilLiqFlx),intent(in) :: out_soilLiqFlx              ! class object for intent(out) soilLiqFlx arguments
  integer(i4b),intent(out)              :: err                         ! error code
  character(*),intent(out)              :: cmessage                    ! error message from groundwatr
  ! intent(out) arguments
  err                     =out_soilLiqFlx % err                        ! intent(out):   error code
  cmessage                =out_soilLiqFlx % cmessage                   ! intent(out):   error message
 end subroutine finalize_out_soilLiqFlx
 ! **** end soilLiqFlx ****

 ! **** groundwatr ****
 subroutine finalize_out_groundwatr(out_groundwatr,err,cmessage)
  class(out_type_groundwatr),intent(in) :: out_groundwatr              ! class object for intent(out) groundwatr arguments
  integer(i4b),intent(out)              :: err                         ! error code
  character(*),intent(out)              :: cmessage                    ! error message from groundwatr
   err               = out_groundwatr % err                                          ! intent(out):   error code
   cmessage          = out_groundwatr % cmessage                                     ! intent(out):   error message
 end subroutine finalize_out_groundwatr
 ! **** end groundwatr ****

 ! **** bigAquifer ****
 subroutine finalize_out_bigAquifer(out_bigAquifer,err,cmessage)
  class(out_type_bigAquifer),intent(in) :: out_bigAquifer ! class object for intent(out) bigAquifer arguments
  integer(i4b),intent(out)              :: err            ! error code
  character(*),intent(out)              :: cmessage       ! error message from bigAquifer
   err                    = out_bigAquifer % err                         ! intent(out):   error code
   cmessage               = out_bigAquifer % cmessage                    ! intent(out):   error message
 end subroutine finalize_out_bigAquifer
 ! **** end bigAquifer ****

 ! **** diagv_node ****
 subroutine finalize_out_diagv_node(out_diagv_node,err,cmessage)
  class(out_type_diagv_node),intent(in)  :: out_diagv_node ! class object for output diagv_node variables
  integer(i4b),intent(out)  :: err                         ! error code
  character(*),intent(out)  :: cmessage                    ! error message

   ! output: error control
   err      = out_diagv_node % err     ! error code
   cmessage = out_diagv_node % message ! error message
 end subroutine finalize_out_diagv_node
 ! **** end diagv_node ****

 ! **** surfaceFlx ****

 subroutine finalize_out_surfaceFlx(out_surfaceFlx,err,cmessage)
  class(out_type_surfaceFlx),intent(in)  :: out_surfaceFlx ! output object for surfaceFlx
  integer(i4b),intent(out)  :: err       ! error code
  character(*),intent(out)  :: cmessage  ! error message

  ! intent(out): error control
  err      = out_surfaceFlx % err     ! error code
  cmessage = out_surfaceFlx % message ! error message
 end subroutine finalize_out_surfaceFlx
 ! **** end surfaceFlx ****

 ! **** iLayerFlux ****
 subroutine finalize_out_iLayerFlux(out_iLayerFlux,err,cmessage)
  class(out_type_iLayerFlux),intent(in)  :: out_iLayerFlux ! class object for output iLayerFlux variables
  integer(i4b),intent(out)  :: err       ! error code
  character(*),intent(out)  :: cmessage  ! error message

   err      = out_iLayerFlux % err     ! error code
   cmessage = out_iLayerFlux % message ! error message
 end subroutine finalize_out_iLayerFlux
 ! **** end iLayerFlux ****

 ! **** qDrainFlux ****
 subroutine finalize_out_qDrainFlux(out_qDrainFlux,err,cmessage)
  class(out_type_qDrainFlux),intent(in) :: out_qDrainFlux  ! class object for output qDrainFlux variables
  integer(i4b),intent(out)  :: err                         ! error code
  character(*),intent(out)  :: cmessage                    ! error message

  ! intent(out): error control
  err      = out_qDrainFlux % err     ! error code 
  cmessage = out_qDrainFlux % message ! error message
 end subroutine finalize_out_qDrainFlux
 ! **** end qDrainFlux ****

 ! **** stateFilter ****
 subroutine finalize_out_stateFilter(out_stateFilter,err,cmessage)
  class(out_type_stateFilter),intent(in) :: out_stateFilter   ! class object for intent(out) stateFilter arguments
  integer(i4b),intent(out)               :: err               ! intent(out): error code
  character(*),intent(out)               :: cmessage          ! intent(out): error message
  err      = out_stateFilter % err                            ! intent(out): error code
  cmessage = out_stateFilter % cmessage                       ! intent(out): error message
 end subroutine finalize_out_stateFilter
 ! **** end stateFilter ****

 ! **** indexSplit ****
 subroutine initialize_in_indexSplit(in_indexSplit,nSubset)
  class(in_type_indexSplit),intent(out) :: in_indexSplit    ! class object for intent(in) indexSplit arguments
  integer(i4b),intent(in)               :: nSubset          ! intent(in): number of states in the subset
  ! in_indexSplit % nSnow   = nSnow                           ! intent(in): number of snow layers          
  ! in_indexSplit % nSoil   = nSoil                           ! intent(in): number of soil layers
  ! in_indexSplit % nLayers = nLayers                         ! intent(in): total number of layers
  in_indexSplit % nSubset = nSubset                         ! intent(in): number of states in the subset
 end subroutine initialize_in_indexSplit

 subroutine finalize_out_indexSplit(out_indexSplit,err,cmessage)
  class(out_type_indexSplit),intent(in) :: out_indexSplit   ! class object for intent(out) indexSplit arguments
  integer(i4b),intent(out)              :: err              ! intent(out): error code
  character(*),intent(out)              :: cmessage         ! intent(out): error message
  err      = out_indexSplit % err                           ! intent(out): error code    
  cmessage = out_indexSplit % cmessage                      ! intent(out): error message
 end subroutine finalize_out_indexSplit
 ! **** end indexSplit ****

 ! **** varSubstep ****
 subroutine initialize_in_varSubstep(in_varSubstep,dt,dtInit,dt_min,whole_step,nSubset,&
                                     doAdjustTemp,firstSubStep,ixSolution,scalar,iStateSplit)
  class(in_type_varSubstep),intent(out) :: in_varSubstep  ! class object for intent(in) varSubstep arguments
  real(rkind),intent(in)                :: dt             ! time step (s)
  real(rkind),intent(in)                :: dtInit         ! initial time step (s)
  real(rkind),intent(in)                :: dt_min         ! minimum time step (s) 
  real(rkind),intent(in)                :: whole_step     ! length of whole step for surface drainage and average flux
  integer(i4b),intent(in)               :: nSubset        ! total number of variables in the state subset
  logical(lgt),intent(in)               :: doAdjustTemp   ! flag to indicate if we adjust the temperature
  logical(lgt),intent(in)               :: firstSubStep   ! flag to denote first sub-step
  integer(i4b),intent(in)               :: ixSolution     ! index of solution method
  integer(i4b),intent(in)               :: scalar         ! scalar solution method
  integer(i4b),intent(in)               :: iStateSplit    ! index of the layer in the splitting operation
 
  ! intent(in) arguments
  in_varSubstep % dt             = dt                     ! intent(in): time step (s)
  in_varSubstep % dtInit         = dtInit                 ! intent(in): initial time step (s)
  in_varSubstep % dt_min         = dt_min                 ! intent(in): minimum time step (s)
  in_varSubstep % whole_step     = whole_step             ! intent(in): length of whole step for surface drainage and average flux
  in_varSubstep % nSubset        = nSubset                ! intent(in): total number of variables in the state subset
  in_varSubstep % doAdjustTemp   = doAdjustTemp           ! intent(in): flag to indicate if we adjust the temperature
  in_varSubstep % firstSubStep   = firstSubStep           ! intent(in): flag to denote first sub-step
  ! in_varSubstep % computeVegFlux = computeVegFlux         ! intent(in): flag to denote if computing energy flux over vegetation
  in_varSubstep % scalarSolution = (ixSolution==scalar)   ! intent(in): flag to denote computing the scalar solution
  in_varSubstep % iStateSplit    = iStateSplit            ! intent(in): index of the layer in the splitting operation
  ! in_varSubstep % fluxMask       = fluxMask               ! intent(in): mask for the fluxes used in this given state subset
 end subroutine initialize_in_varSubstep

 subroutine initialize_io_varSubstep(io_varSubstep,firstFluxCall)
  class(io_type_varSubstep),intent(out) :: io_varSubstep  ! class object for intent(inout) varSubstep arguments
  logical(lgt),intent(in)               :: firstFluxCall  ! flag to indicate if we are processing the first flux call

  ! intent(inout) arguments
  io_varSubstep % firstFluxCall = firstFluxCall           ! intent(inout): flag to indicate if we are processing the first flux call
  ! io_varSubstep % fluxCount     = fluxCount               ! intent(inout): number of times fluxes are updated (should equal nsubstep)
  ! io_varSubstep % ixSaturation  = ixSaturation            ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
 end subroutine initialize_io_varSubstep

 subroutine finalize_io_varSubstep(io_varSubstep,firstFluxCall)
  class(io_type_varSubstep),intent(in)  :: io_varSubstep  ! class object for intent(inout) varSubstep arguments
  logical(lgt),intent(out)              :: firstFluxCall  ! flag to indicate if we are processing the first flux call

  ! intent(inout) arguments
  firstFluxCall = io_varSubstep % firstFluxCall           ! intent(inout): flag to indicate if we are processing the first flux call
  ! fluxCount     = io_varSubstep % fluxCount               ! intent(inout): number of times fluxes are updated (should equal nsubstep)
  ! ixSaturation  = io_varSubstep % ixSaturation            ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
 end subroutine finalize_io_varSubstep

 subroutine finalize_out_varSubstep(out_varSubstep,dtMultiplier,nSubsteps,failedMinimumStep,reduceCoupledStep,tooMuchMelt,err,cmessage)
  class(out_type_varSubstep),intent(in) :: out_varSubstep    ! class object for intent(out) varSubstep arguments
  real(rkind),intent(out)               :: dtMultiplier      ! substep multiplier (-)
  integer(i4b),intent(out)              :: nSubsteps         ! number of substeps taken for a given split
  logical(lgt),intent(out)              :: failedMinimumStep ! flag for failed substeps
  logical(lgt),intent(out)              :: reduceCoupledStep ! flag to reduce the length of the coupled step
  logical(lgt),intent(out)              :: tooMuchMelt       ! flag to denote that ice is insufficient to support melt
  integer(i4b),intent(out)              :: err               ! error code
  character(*),intent(out)              :: cmessage          ! error message                                          

  ! intent(out) arguments
  dtMultiplier      = out_varSubstep % dtMultiplier       ! intent(out): substep multiplier (-)
  nSubsteps         = out_varSubstep % nSubsteps          ! intent(out): number of substeps taken for a given split
  failedMinimumStep = out_varSubstep % failedMinimumStep  ! intent(out): flag for failed substeps
  reduceCoupledStep = out_varSubstep % reduceCoupledStep  ! intent(out): flag to reduce the length of the coupled step
  tooMuchMelt       = out_varSubstep % tooMuchMelt        ! intent(out): flag to denote that ice is insufficient to support melt
  err               = out_varSubstep % err                ! intent(out): error code
  cmessage          = out_varSubstep % cmessage           ! intent(out): error message                                          
 end subroutine finalize_out_varSubstep
 ! **** end varSubstep ****

 ! **** computJacob ****
 subroutine initialize_in_computJacob(in_computJacob,dt,nSnow,nSoil,nLayers,computeVegFlux,computeBaseflow,ixMatrix)
  class(in_type_computJacob),intent(out) :: in_computJacob           ! class object for intent(in) computJacob arguments
  real(rkind),intent(in)              :: dt                          ! intent(in): length of the time step (seconds)
  integer(i4b),intent(in)             :: nSnow                       ! intent(in): number of snow layers
  integer(i4b),intent(in)             :: nSoil                       ! intent(in): number of soil layers
  integer(i4b),intent(in)             :: nLayers                     ! intent(in): total number of layers in the snow+soil domain
  logical(lgt),intent(in)             :: computeVegFlux              ! intent(in): flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)             :: computeBaseflow             ! intent(in): flag to indicate if computing baseflow
  integer(i4b),intent(in)             :: ixMatrix                    ! intent(in): form of the Jacobian matrix                         
 
  ! intent(in) arguments
  in_computJacob % dt               =  dt                            ! intent(in): length of the time step (seconds)                    
  in_computJacob % nSnow            =  nSnow                         ! intent(in): number of snow layers
  in_computJacob % nSoil            =  nSoil                         ! intent(in): number of soil layers
  in_computJacob % nLayers          =  nLayers                       ! intent(in): total number of layers in the snow+soil domain
  in_computJacob % computeVegFlux   =  computeVegFlux                ! intent(in): flag to indicate if computing fluxes over vegetation
  in_computJacob % computeBaseflow  =  computeBaseflow               ! intent(in): flag to indicate if computing baseflow
  in_computJacob % ixMatrix         =  ixMatrix                      ! intent(in): form of the Jacobian matrix                         
 end subroutine initialize_in_computJacob

 subroutine finalize_out_computJacob(out_computJacob,err,cmessage)
  class(out_type_computJacob),intent(in) :: out_computJacob           ! class object for intent(out) computJacob arguments
  integer(i4b),intent(out)               :: err                       ! intent(out): error code
  character(*),intent(out)               :: cmessage                  ! intent(out): error message
  ! intent(out) arguments
  err               = out_computJacob % err                           ! intent(out): error code
  cmessage          = out_computJacob % cmessage                      ! intent(out): error message                                          
 end subroutine finalize_out_computJacob

 ! **** lineSearchRefinement ****
 subroutine initialize_in_lineSearchRefinement(in_lineSearchRefinement,doSearch,fOld)
  class(in_type_lineSearchRefinement),intent(out) :: in_lineSearchRefinement   ! class object for intent(out) arguments
  logical(lgt),intent(in)                         :: doSearch                  ! intent(in): flag to do the line search
  real(rkind) ,intent(in)                         :: fOld                      ! intent(in): old function value
  in_lineSearchRefinement % doSearch     = doSearch                  ! intent(in): flag to do the line search
  in_lineSearchRefinement % fOld         = fOld                      ! intent(in): old function value
 end subroutine initialize_in_lineSearchRefinement
 
 subroutine finalize_out_lineSearchRefinement(out_lineSearchRefinement,fNew,converged,err,message)
  class(out_type_lineSearchRefinement),intent(in) :: out_lineSearchRefinement   ! class object for intent(out) arguments
  real(rkind) ,intent(out)   :: fNew                                  ! intent(out): new function evaluation
  logical(lgt),intent(out)   :: converged                             ! intent(out): convergence flag
  integer(i4b),intent(out)   :: err                                   ! intent(out): error code
  character(*),intent(out)   :: message                               ! intent(out): error message
  fNew      = out_lineSearchRefinement % fNew                         ! intent(out): new function evaluation
  converged = out_lineSearchRefinement % converged                    ! intent(out): convergence flag
  err       = out_lineSearchRefinement % err                          ! intent(out): error code
  message   = out_lineSearchRefinement % message                      ! intent(out): error message
 end subroutine finalize_out_lineSearchRefinement

 ! **** summaSolve4homegrown ****

 subroutine initialize_in_summaSolve4homegrown(in_SS4NR,dt_cur,dt,iter,nSnow,nSoil,nLayers,nLeadDim,nState,ixMatrix,firstSubStep,computeVegFlux,scalarSolution,fOld)
  class(in_type_summaSolve4homegrown),intent(out)    :: in_SS4NR   ! class object for intent(out) arguments
  real(rkind) ,intent(in) :: dt_cur                   ! intent(in): current stepsize
  real(rkind) ,intent(in) :: dt                       ! intent(in): entire time step for drainage pond rate
  integer(i4b),intent(in) :: iter                     ! intent(in): iteration index
  integer(i4b),intent(in) :: nSnow                    ! intent(in): number of snow layers
  integer(i4b),intent(in) :: nSoil                    ! intent(in): number of soil layers
  integer(i4b),intent(in) :: nLayers                  ! intent(in): total number of layers
  integer(i4b),intent(in) :: nLeadDim                 ! intent(in): length of the leading dimension of the Jacobian matrix (nBands or nState)
  integer(i4b),intent(in) :: nState                   ! intent(in): total number of state variables
  integer(i4b),intent(in) :: ixMatrix                 ! intent(in): type of matrix (full or band diagonal)
  logical(lgt),intent(in) :: firstSubStep             ! intent(in): flag to indicate if we are processing the first sub-step
  logical(lgt),intent(in) :: computeVegFlux           ! intent(in): flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in) :: scalarSolution           ! intent(in): flag to denote if implementing the scalar solution
  real(rkind) ,intent(in) :: fOld                     ! intent(in): old function evaluation

  in_SS4NR % dt_cur         = dt_cur        
  in_SS4NR % dt             = dt            
  in_SS4NR % iter           = iter          
  in_SS4NR % nSnow          = nSnow         
  in_SS4NR % nSoil          = nSoil         
  in_SS4NR % nLayers        = nLayers       
  in_SS4NR % nLeadDim       = nLeadDim      
  in_SS4NR % nState         = nState        
  in_SS4NR % ixMatrix       = ixMatrix      
  in_SS4NR % firstSubStep   = firstSubStep  
  in_SS4NR % computeVegFlux = computeVegFlux 
  in_SS4NR % scalarSolution = scalarSolution
  in_SS4NR % fOld           = fOld            
 end subroutine initialize_in_summaSolve4homegrown

 subroutine initialize_io_summaSolve4homegrown(io_SS4NR,firstFluxCall,xMin,xMax,ixSaturation)
  class(io_type_summaSolve4homegrown),intent(out)    :: io_SS4NR   ! class object for intent(inout) arguments
  logical(lgt),intent(in) :: firstFluxCall ! intent(inout): flag to indicate if we are processing the first flux call
  real(rkind) ,intent(in) :: xMin,xMax     ! intent(inout): brackets of the root
  integer(i4b),intent(in) :: ixSaturation  ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)

  io_SS4NR % firstFluxCall = firstFluxCall
  io_SS4NR % xMin          = xMin    
  io_SS4NR % xMax          = xMax    
  io_SS4NR % ixSaturation  = ixSaturation 
 end subroutine initialize_io_summaSolve4homegrown

 subroutine finalize_io_summaSolve4homegrown(io_SS4NR,firstFluxCall,xMin,xMax,ixSaturation)
  class(io_type_summaSolve4homegrown),intent(in)    :: io_SS4NR   ! class object for intent(inout) arguments
  logical(lgt),intent(out) :: firstFluxCall ! intent(inout): flag to indicate if we are processing the first flux call
  real(rkind) ,intent(out) :: xMin,xMax     ! intent(inout): brackets of the root
  integer(i4b),intent(out) :: ixSaturation  ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)

  firstFluxCall = io_SS4NR % firstFluxCall
  xMin          = io_SS4NR % xMin    
  xMax          = io_SS4NR % xMax    
  ixSaturation  = io_SS4NR % ixSaturation 
 end subroutine finalize_io_summaSolve4homegrown

 subroutine finalize_out_summaSolve4homegrown(out_SS4NR,fNew,converged,err,message)
  class(out_type_summaSolve4homegrown),intent(in)    :: out_SS4NR   ! class object for intent(out) arguments
  real(rkind) ,intent(out) :: fNew      ! intent(out): new function evaluation
  logical(lgt),intent(out) :: converged ! intent(out): convergence flag
  integer(i4b),intent(out) :: err       ! intent(out): error code
  character(*),intent(out) :: message   ! intent(out): error message

  fNew      = out_SS4NR % fNew       
  converged = out_SS4NR % converged
  err       = out_SS4NR % err      
  message   = out_SS4NR % message  
 end subroutine finalize_out_summaSolve4homegrown

END MODULE data_types
