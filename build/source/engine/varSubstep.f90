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

module varSubstep_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access the global print flag
USE globalData,only:globalPrintFlag

! domain types
USE globalData,only:iname_cas       ! named variables for the canopy air space
USE globalData,only:iname_veg       ! named variables for vegetation
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! global metadata
USE globalData,only:flux_meta       ! metadata on the model fluxes

! derived types to define the data structures
USE data_types,only:&
                    var_i,              & ! data vector (i4b)
                    var_d,              & ! data vector (rkind)
                    var_flagVec,        & ! data vector with variable length dimension (i4b)
                    var_ilength,        & ! data vector with variable length dimension (i4b)
                    var_dlength,        & ! data vector with variable length dimension (rkind)
                    zLookup,            & ! lookup tables
                    model_options,      & ! defines the model decisions
                    in_type_varSubstep, & ! class for intent(in) arguments
                    io_type_varSubstep, & ! class for intent(inout) arguments
                    out_type_varSubstep   ! class for intent(out) arguments

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure

! look up structure for variable types
USE var_lookup,only:iLookVarType

! constants
USE multiconst,only:&
                    Tfreeze,        & ! freezing temperature                 (K)
                    LH_fus,         & ! latent heat of fusion                (J kg-1)
                    LH_vap,         & ! latent heat of vaporization          (J kg-1)
                    iden_ice,       & ! intrinsic density of ice             (kg m-3)
                    iden_water        ! intrinsic density of liquid water    (kg m-3)

! look-up values for the numerical method
USE mDecisions_module,only:         &
                    homegrown      ,& ! homegrown backward Euler solution using concepts from numerical recipes
                    kinsol         ,& ! SUNDIALS backward Euler solution using Kinsol
                    ida               ! SUNDIALS solution using IDA

! look-up values for the choice of variable in energy equations (BE residual or IDA state variable)
USE mDecisions_module,only:         &
                    closedForm,     & ! use temperature with closed form heat capacity
                    enthalpyFormLU, & ! use enthalpy with soil temperature-enthalpy lookup tables
                    enthalpyForm      ! use enthalpy with soil temperature-enthalpy analytical solution

! safety: set private unless specified otherwise
implicit none
private
public::varSubstep

! algorithmic parameters
real(rkind),parameter     :: verySmall=1.e-6_rkind   ! used as an additive constant to check if substantial difference among real numbers

contains


! **********************************************************************************************************
! public subroutine varSubstep: run the model for a collection of substeps for a given state subset
! **********************************************************************************************************
subroutine varSubstep(&
                      ! input: model control
                      in_varSubstep,     & ! intent(in)    : model control
                      nGRU, &
                      computeVegFlux, &
                      io_varSubstep,     & ! intent(inout) : model control
                      ! input/output: data structures
                      model_decisions,   & ! intent(in)    : model decisions
                      decisions, veg_param, &
                      lookup_data,       & ! intent(in)    : lookup tables
                      type_data,         & ! intent(in)    : type of vegetation and soil
                      attr_data,         & ! intent(in)    : spatial attributes
                      forc_data,         & ! intent(in)    : model forcing data
                      mpar_data,         & ! intent(in)    : model parameters
                      indx_data,         & ! intent(inout) : index data
                      prog_data,         & ! intent(inout) : model prognostic variables for a local HRU
                      diag_data,         & ! intent(inout) : model diagnostic variables for a local HRU
                      flux_data,         & ! intent(inout) : model fluxes for a local HRU
                      flux_mean,         & ! intent(inout) : mean model fluxes for a local HRU
                      deriv_data,        & ! intent(inout) : derivatives in model fluxes w.r.t. relevant state variables
                      bvar_data,         & ! intent(in)    : model variables for the local basin
                      ! output: model control
                      out_varSubstep)      ! intent(out)   : model control
  ! ---------------------------------------------------------------------------------------
  ! structure allocations
  USE allocspace_module,only:allocLocal                ! allocate local data structures
  ! simulation of fluxes and residuals given a trial state vector
  USE getVectorz_module,only:popStateVec                ! populate the state vector
  USE getVectorz_module,only:varExtract                 ! extract variables from the state vector
  USE systemSolv_module,only:systemSolv                 ! solve the system of equations for one time step
  ! identify name of variable type (for error message)
  USE get_ixName_module,only:get_varTypeName           ! to access type strings for error messages
  use device_data_types
  use initialize_device
  implicit none
  ! ---------------------------------------------------------------------------------------
  ! * dummy variables
  ! ---------------------------------------------------------------------------------------
  ! input: model control
  type(in_type_varSubstep),intent(in)    :: in_varSubstep             ! model control
  logical(lgt),device :: computeVegFlux(:)
  type(io_type_varSubstep),intent(inout) :: io_varSubstep             ! model control
  ! input/output: data structures
  type(model_options),intent(in)         :: model_decisions(:)        ! model decisions
  type(zLookup),intent(in)               :: lookup_data               ! lookup tables
  type(type_data_device),intent(in)                 :: type_data                 ! type of vegetation and soil
  type(attr_data_device),intent(in)                 :: attr_data                 ! spatial attributes
  type(forc_data_device),intent(in)                 :: forc_data                 ! model forcing data
  type(mpar_data_device),intent(in)           :: mpar_data                 ! model parameters
  type(indx_data_device),intent(inout)        :: indx_data                 ! indices for a local HRU
  type(prog_data_device),intent(inout)        :: prog_data                 ! prognostic variables for a local HRU
  type(diag_data_device),intent(inout)        :: diag_data                 ! diagnostic variables for a local HRU
  type(flux_data_device),intent(inout)        :: flux_data                 ! model fluxes for a local HRU
  type(flux_data_device),intent(inout)        :: flux_mean                 ! mean model fluxes for a local HRU
  type(deriv_data_device),intent(inout)        :: deriv_data                ! derivatives in model fluxes w.r.t. relevant state variables
  type(bvar_data_device),intent(in)           :: bvar_data                 ! model variables for the local basin
  ! output: model control
  type(out_type_varSubstep),intent(out)  :: out_varSubstep            ! model control
  ! ---------------------------------------------------------------------------------------
  ! * general local variables
  ! ---------------------------------------------------------------------------------------
  ! error control
  character(LEN=256)                 :: cmessage                               ! error message of downwind routine
  ! general local variables
  integer(i4b)                       :: iVar                                   ! index of variables in data structures
  integer(i4b)                       :: iSoil                                  ! index of soil layers
  integer(i4b)                       :: ixLayer                                ! index in a given domain
  integer(i4b),dimension(1)          :: ixMin,ixMax                            ! bounds of a given flux vector
  ! time stepping
  real(rkind)                        :: dtSum                                  ! sum of time from successful steps (seconds)
  real(rkind)                        :: dt_wght                                ! weight given to a given flux calculation
  real(rkind)                        :: dtSubstep                              ! length of a substep (s)
  real(rkind)                        :: maxstep                                ! maximum time step length (seconds)
  integer(i4b)                       :: nSteps                                 ! number of time steps taken in solver
  ! adaptive sub-stepping for the solution
  logical(lgt)                       :: failedSubstep                          ! flag to denote success of substepping for a given split
  integer(i4b)                       :: niter                                  ! number of iterations taken
  integer(i4b),parameter             :: n_inc=5                                ! minimum number of iterations to increase time step
  integer(i4b),parameter             :: n_dec=15                               ! maximum number of iterations to decrease time step
  real(rkind),parameter              :: F_inc = 1.25_rkind                     ! factor used to increase time step
  real(rkind),parameter              :: F_dec = 0.90_rkind                     ! factor used to decrease time step
  ! state and flux vectors (Note: nstate = in_varSubstep % nSubset)
  real(rkind),device                        :: untappedMelt(in_varSubstep % nSubset,nGRU)  ! un-tapped melt energy (J m-3 s-1)
  real(rkind),device                        :: stateVecInit(in_varSubstep % nSubset,nGRU)  ! initial state vector (mixed units)
  real(rkind),device                        :: stateVecTrial(in_varSubstep % nSubset,nGRU) ! trial state vector (mixed units)
  real(rkind),device                        :: stateVecPrime(in_varSubstep % nSubset,nGRU) ! trial state vector (mixed units)
  type(flux_data_device)                  :: flux_temp                              ! temporary model fluxes
  ! flags
  logical(lgt)                       :: firstSplitOper                         ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt)                       :: waterBalanceError                      ! flag to denote that there is a water balance error
  logical(lgt)                       :: nrgFluxModified                        ! flag to denote that the energy fluxes were modified
  ! energy fluxes
  real(rkind),device                        :: sumCanopyEvaporation(nGRU)                   ! sum of canopy evaporation/condensation (kg m-2 s-1)
  real(rkind),device                        :: sumLatHeatCanopyEvap(nGRU)                   ! sum of latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
  real(rkind),device                        :: sumSenHeatCanopy(nGRU)                       ! sum of sensible heat flux from the canopy to the canopy air space (W m-2)
  real(rkind),device                        :: sumSoilCompress(nGRU)                        ! sum of total soil compression
  real(rkind),allocatable,device            :: sumLayerCompress(:,:)                    ! sum of soil compression by layer
  ! balances and residual vectors
  real(rkind),device                        :: fluxVec(in_varSubstep % nSubset,nGRU)       ! substep flux vector (mixed units)
  real(rkind),device                        :: resSink(in_varSubstep % nSubset,nGRU)       ! substep sink terms on the RHS of the state equation
  real(qp),device                           :: resVec(in_varSubstep % nSubset,nGRU)        ! substep residual vector
  real(rkind),device                        :: balance(in_varSubstep % nSubset,nGRU)       ! substep balance per second
  real(rkind),device                        :: sumBalance(in_varSubstep % nSubset,nGRU)    ! sum of substeps balance
  logical(lgt),parameter             :: computMassBalance = .true.             ! flag to compute the mass balance, will affect step length, default true
  logical(lgt),parameter             :: computNrgBalance = .true.              ! flag to compute the energy balance, will not effect solution but will not compute energy balance if false (saves expense)
  logical(lgt)                       :: computeEnthTemp                        ! flag to compute enthalpy regardless of the model decision
  logical(lgt)                       :: enthalpyStateVec                       ! flag if enthalpy is a state variable (ida)
  logical(lgt)                       :: use_lookup                             ! flag to use the lookup table for soil enthalpy, otherwise use analytical solution

  integer(i4b) :: nGRU
  type(decisions_device) :: decisions
  type(veg_parameters) :: veg_param
  integer(i4b) :: iGRU
  integer(i4b) :: nSoil
  ! ---------------------------------------------------------------------------------------
  ! initialize error control
  nSoil = indx_data%nSoil
  ! nLayers = indx_data%nLayers_d
  out_varSubstep % err=0; out_varSubstep % cmessage='varSubstep/'


  ! ---------------------------------------------------------------------------------------
  ! point to variables in the data structures
  ! ---------------------------------------------------------------------------------------
  globalVars: associate(&
    ! input: model control
    dt             => in_varSubstep % dt,             & ! intent(in): time step (seconds)
    dtInit         => in_varSubstep % dtInit,         & ! intent(in): initial time step (seconds)
    dt_min         => in_varSubstep % dt_min,         & ! intent(in): minimum time step (seconds)
    whole_step     => in_varSubstep % whole_step,     & ! intent(in): length of whole step for surface drainage and average flux
    nState         => in_varSubstep % nSubset,        & ! intent(in): total number of state variables
    doAdjustTemp   => in_varSubstep % doAdjustTemp,   & ! intent(in): flag to indicate if we adjust the temperature
    firstSubStep   => in_varSubstep % firstSubStep,   & ! intent(in): flag to indicate if processing the first sub-step
    ! computeVegFlux => in_varSubstep % computeVegFlux, & ! intent(in): flag to indicate if computing fluxes over vegetation (.false. means veg is buried with snow)
    scalarSolution => in_varSubstep % scalarSolution, & ! intent(in): flag to denote implementing the scalar solution
    iStateSplit    => in_varSubstep % iStateSplit,    & ! intent(in): index of the state in the splitting operation
    fluxMask       => in_varSubstep % fluxMask,       & ! intent(in): flags to denote if the flux is calculated in the given state subset
    firstFluxCall  => io_varSubstep % firstFluxCall,  & ! intent(inout): flag to define the first flux call
    fluxCount      => io_varSubstep % fluxCount,      & ! intent(inout): number of times that the flux is updated (should equal nSubsteps)
    ixSaturation   => io_varSubstep % ixSaturation,   & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
    ! model decisions
    ixNumericalMethod       => model_decisions(iLookDECISIONS%num_method)%iDecision   ,& ! intent(in):    [i4b]    choice of numerical solver
    ixNrgConserv            => model_decisions(iLookDECISIONS%nrgConserv)%iDecision   ,& ! intent(in):    [i4b]    choice of variable in either energy backward Euler residual or IDA state variable
    ! number of layers
    ! nSnow                   => indx_data%nSnow                 ,& ! intent(in):    [i4b]    number of snow layers
    nSoil                   => indx_data%nSoil                 ,& ! intent(in):    [i4b]    number of soil layers
    nLayers                 => indx_data%nLayers_d               ,& ! intent(in):    [i4b]    total number of layers
    ! nSoilOnlyHyd            => indx_data%nSoilOnlyHyd         ,& ! intent(in):    [i4b]    number of hydrology variables in the soil domain
    mLayerDepth             => prog_data%mLayerDepth               ,& ! intent(in):    [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! get indices for balances
    ixCasNrg                => indx_data%ixCasNrg              ,& ! intent(in):    [i4b]    index of canopy air space energy state variable
    ixVegNrg                => indx_data%ixVegNrg              ,& ! intent(in):    [i4b]    index of canopy energy state variable
    ixVegHyd                => indx_data%ixVegHyd              ,& ! intent(in):    [i4b]    index of canopy hydrology state variable (mass)
    ixTopNrg                => indx_data%ixTopNrg              ,& ! intent(in):    [i4b]    index of upper-most energy state in the snow+soil subdomain
    ixTopHyd                => indx_data%ixTopHyd              ,& ! intent(in):    [i4b]    index of upper-most hydrology state in the snow+soil subdomain
    ixAqWat                 => indx_data%ixAqWat               ,& ! intent(in):    [i4b]    index of water storage in the aquifer
    ixSoilOnlyHyd           => indx_data%ixSoilOnlyHyd         ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the soil domain
    ixSnowSoilHyd           => indx_data%ixSnowSoilHyd         ,& ! intent(in):    [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    ixSnowSoilNrg           => indx_data%ixSnowSoilNrg         ,& ! intent(in):    [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    nSnowSoilNrg            => indx_data%nLayers_d          ,& ! intent(in):    [i4b]    number of energy state variables in the snow+soil domain
    nSnowSoilHyd            => indx_data%nLayers_d          ,& ! intent(in):    [i4b]    number of hydrology state variables in the snow+soil domain
    ! mapping between state vectors and control volumes
    ! ixLayerActive           => indx_data%ixLayerActive            ,& ! intent(in):    [i4b(:)] list of indices for all active layers (inactive=integerMissing)
    ixMapFull2Subset        => indx_data%ixMapFull2Subset         ,& ! intent(in):    [i4b(:)] mapping of full state vector to the state subset
    ixControlVolume         => indx_data%ixControlVolume          ,& ! intent(in):    [i4b(:)] index of control volume for different domains (veg, snow, soil)
    ! model state variables (vegetation canopy)
    scalarCanairTemp        => prog_data%scalarCanairTemp       ,& ! intent(inout): [dp]     temperature of the canopy air space (K)
    scalarCanopyTemp        => prog_data%scalarCanopyTemp       ,& ! intent(inout): [dp]     temperature of the vegetation canopy (K)
    scalarCanopyIce         => prog_data%scalarCanopyIce        ,& ! intent(inout): [dp]     mass of ice on the vegetation canopy (kg m-2)
    scalarCanopyLiq         => prog_data%scalarCanopyLiq        ,& ! intent(inout): [dp]     mass of liquid water on the vegetation canopy (kg m-2)
    scalarCanopyWat         => prog_data%scalarCanopyWat        ,& ! intent(inout): [dp]     mass of total water on the vegetation canopy (kg m-2)
    ! model state variables (snow and soil domains)
    mLayerTemp              => prog_data%mLayerTemp                ,& ! intent(inout): [dp(:)]  temperature of each snow/soil layer (K)
    mLayerVolFracIce        => prog_data%mLayerVolFracIce          ,& ! intent(inout): [dp(:)]  volumetric fraction of ice (-)
    mLayerVolFracLiq        => prog_data%mLayerVolFracLiq          ,& ! intent(inout): [dp(:)]  volumetric fraction of liquid water (-)
    mLayerVolFracWat        => prog_data%mLayerVolFracWat          ,& ! intent(inout): [dp(:)]  volumetric fraction of total water (-)
    mLayerMatricHead        => prog_data%mLayerMatricHead          ,& ! intent(inout): [dp(:)]  matric head (m)
    mLayerMatricHeadLiq     => diag_data%mLayerMatricHeadLiq       ,& ! intent(inout): [dp(:)]  matric potential of liquid water (m)
    ! model control
    dtMultiplier      => out_varSubstep % dtMultiplier             ,& ! intent(out): substep multiplier (-)
    nSubsteps         => out_varSubstep % nSubsteps                ,& ! intent(out): number of substeps taken for a given split
    failedMinimumStep => out_varSubstep % failedMinimumStep        ,& ! intent(out): flag to denote success of substepping for a given split
    reduceCoupledStep => out_varSubstep % reduceCoupledStep        ,& ! intent(out): flag to denote need to reduce the length of the coupled step
    tooMuchMelt       => out_varSubstep % tooMuchMelt              ,& ! intent(out): flag to denote that ice is insufficient to support melt
    err               => out_varSubstep % err                      ,& ! intent(out): error code
    message           => out_varSubstep % cmessage                  & ! intent(out): error message
    )  ! end association with variables in the data structures
    ! *********************************************************************************************************************************************************

    ! initialize flag for the success of the substepping
    failedMinimumStep=.false.

    ! set the flag to compute enthalpy, may want to have this true always if want to output enthalpy
    computeEnthTemp  = .false.
    enthalpyStateVec = .false.
    use_lookup       = .false.
    if((ixNrgConserv .ne. closedForm .or. computNrgBalance) .and. ixNumericalMethod .ne. ida) computeEnthTemp = .true. ! use enthTemp to conserve energy or compute energy balance
    if(ixNrgConserv .ne. closedForm .and. ixNumericalMethod==ida) enthalpyStateVec = .true. ! enthalpy as state variable
    if(ixNrgConserv==enthalpyFormLU) use_lookup = .true. ! use lookup tables for soil enthalpy instead of analytical solution

    ! initialize the length of the substep
    dtSubstep = dtInit

    ! change maxstep with hard code here to make only the newton step loop in systemSolv* happen more frequently
    !   NOTE: this may just be amplifying the splitting error if maxstep is smaller than the full possible step
    maxstep = mpar_data%maxstep  ! maximum time step (s).

    ! allocate space for the temporary model flux structure
    ! call allocLocal(flux_meta(:),flux_temp,nSnow,nSoil,err,cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! initialize the model fluxes (some model fluxes are not computed in the iterations)
    ! do iVar=1,size(flux_data%var)
      ! flux_temp%var(iVar)%dat(:) = flux_data%var(iVar)%dat(:)
    ! end do
    flux_temp = flux_data

    ! initialize the total energy fluxes (modified in updateProg)
    sumCanopyEvaporation = 0._rkind  ! canopy evaporation/condensation (kg m-2 s-1)
    sumLatHeatCanopyEvap = 0._rkind  ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
    sumSenHeatCanopy     = 0._rkind  ! sensible heat flux from the canopy to the canopy air space (W m-2)
    sumSoilCompress      = 0._rkind  ! total soil compression
    allocate(sumLayerCompress(nSoil,nGRU)); sumLayerCompress = 0._rkind ! soil compression by layer

    ! initialize balances
    sumBalance = 0._rkind

    ! define the first flux call in a splitting operation
    firstSplitOper = (.not.scalarSolution .or. iStateSplit==1)

    ! initialize subStep
    dtSum     = 0._rkind  ! keep track of the portion of the time step that is completed
    nSubsteps = 0

    ! loop through substeps
    ! NOTE: continuous do statement with exit clause
    substeps: do
      dtSubstep = min(dtSubstep,maxstep)

      ! -----
      ! * populate state vectors...
      ! ---------------------------

      ! initialize state vectors
      call popStateVec(&
                      ! input
                      nState,           & ! intent(in):  number of desired state variables
                      nGRU, &
                      enthalpyStateVec, & ! intent(in):  flag to use enthalpy as a state variable
                      prog_data,        & ! intent(in):  model prognostic variables for a local HRU
                      diag_data,        & ! intent(in):  model diagnostic variables for a local HRU
                      indx_data,        & ! intent(in):  indices defining model states and layers
                      ! output
                      stateVecInit,     & ! intent(out): initial model state vector (mixed units)
                      err,cmessage)       ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

      ! -----
      ! * iterative solution...
      ! -----------------------
      ! solve the system of equations for a given state subset
      call systemSolv(&
                      ! input: model control
                      dtSubstep,         & ! intent(in):    time step (s)
                      whole_step,        & ! intent(in):    entire time step (s)
                      nState,            & ! intent(in):    total number of state variables
                      nLayers,           & ! intent(in):    total number of layers
                      nGRU, indx_data%nSnow, &
                      firstSubStep,      & ! intent(in):    flag to denote first sub-step
                      firstFluxCall,     & ! intent(inout): flag to indicate if we are processing the first flux call
                      firstSplitOper,    & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                      computeVegFlux,    & ! intent(in):    flag to denote if computing energy flux over vegetation
                      scalarSolution,    & ! intent(in):    flag to denote if implementing the scalar solution
                      computMassBalance, & ! intent(in):    flag to compute mass balance
                      computNrgBalance,  & ! intent(in):    flag to compute energy balance
                      ! input/output: data structures
                      lookup_data,       & ! intent(in):    lookup tables
                      type_data,         & ! intent(in):    type of vegetation and soil
                      attr_data,         & ! intent(in):    spatial attributes
                      forc_data,         & ! intent(in):    model forcing data
                      mpar_data,         & ! intent(in):    model parameters
                      indx_data,         & ! intent(inout): index data
                      prog_data,         & ! intent(inout): model prognostic variables for a local HRU
                      diag_data,         & ! intent(inout): model diagnostic variables for a local HRU
                      flux_temp,         & ! intent(inout): model fluxes for a local HRU
                      bvar_data,         & ! intent(in):    model variables for the local basin
                      model_decisions,   & ! intent(in):    model decisions
                      decisions, veg_param, &
                      stateVecInit,      & ! intent(in):    initial state vector
                      ! output: model control
                      deriv_data,        & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      ixSaturation,      & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      stateVecTrial,     & ! intent(out):   updated state vector
                      stateVecPrime,     & ! intent(out):   updated state vector if need the prime space (ida)
                      fluxVec,           & ! intent(out):   model flux vector
                      resSink,           & ! intent(out):   additional (sink) terms on the RHS of the state equation
                      resVec,            & ! intent(out):   residual vector
                      untappedMelt,      & ! intent(out):   un-tapped melt energy (J m-3 s-1)
                      ! output: balances (only computed at this level for ida)
                      balance,           & ! intent(out):   balance per state variable
                      ! output  model control
                      niter,             & ! intent(out):   number of iterations taken (homegrown solver)
                      nSteps,            & ! intent(out):   number of time steps taken in solver
                      reduceCoupledStep, & ! intent(out):   flag to reduce the length of the coupled step
                      tooMuchMelt,       & ! intent(out):   flag to denote that ice is insufficient to support melt
                      err,cmessage)        ! intent(out):   error code and error message

      if(err/=0)then ! (check for errors, but do not fail yet)
        message=trim(message)//trim(cmessage)
        if(err>0) return
      endif
 
      ! if too much melt or need to reduce length of the coupled step then return
      ! NOTE: need to go all the way back to coupled_em and merge snow layers, as all splitting operations need to occur with the same layer geometry
      if(tooMuchMelt .or. reduceCoupledStep)then 
        deallocate(sumLayerCompress)
        return
      endif

      ! identify failure
      failedSubstep = (err<0)

      ! check
      if(globalPrintFlag)then
        print*, 'niter, failedSubstep, dtSubstep = ', niter, failedSubstep, dtSubstep
        print*, trim(cmessage)
      endif

      ! reduce step based on failure
      if(failedSubstep)then
        err=0; message='varSubstep/'  ! recover from failed convergence
        dtMultiplier  = 0.5_rkind        ! system failure: step halving
      else
        ! ** implicit Euler: adjust step length based on iteration count
        if(niter<n_inc)then
          dtMultiplier = F_inc
        elseif(niter>n_dec)then
          dtMultiplier = F_dec
        else
          dtMultiplier = 1._rkind
        endif
      endif  ! switch between failure and success

      ! check if we failed the substep
      if(failedSubstep)then

        ! check that the substep is greater than the minimum step
        if(dtSubstep*dtMultiplier<dt_min)then
          ! --> exit, and either (1) try another solution method; or (2) reduce coupled step
          failedMinimumStep=.true.
          exit subSteps

        else ! step is still OK
          dtSubstep = dtSubstep*dtMultiplier
          cycle subSteps
        endif  ! if step is less than the minimum

      endif  ! if failed the substep

      ! -----
      ! * update model fluxes...
      ! ------------------------

      ! NOTE: if we get to here then we are accepting the step of dtSubstep
      if(err/=0)then
        message=trim(message)//'expect err=0 if updating fluxes'
        return
      endif

      ! update prognostic variables, update balances, and check them for possible step reduction if homegrown or kinsol solver
      call updateProg(dtSubstep,indx_data%nSnow,nSoil,nLayers,nGRU,untappedMelt,stateVecTrial,stateVecPrime,                                    & ! input: states
                      doAdjustTemp,computeVegFlux,computMassBalance,computNrgBalance,computeEnthTemp,enthalpyStateVec,use_lookup,& ! input: model control
                      model_decisions,lookup_data,mpar_data,indx_data,flux_temp,prog_data,diag_data,deriv_data,decisions,                 & ! input-output: data structures
                      fluxVec,resVec,balance,waterBalanceError,nrgFluxModified,err,message)                                        ! input-output: balances, flags, and error control
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        if(err>0) return
      endif

      ! if water balance error then reduce the length of the coupled step
      if(waterBalanceError)then
        message=trim(message)//'water balance error'
        reduceCoupledStep=.true.
        deallocate(sumLayerCompress)
        err=-20; return
      endif

      if(globalPrintFlag)&
      print*, trim(cmessage)//': dt = ', dtSubstep

      ! recover from errors in prognostic update
      if(err<0)then

        ! modify step
        err=0  ! error recovery
        dtSubstep = dtSubstep/2._rkind 

        ! check minimum: fail minimum step if there is an error in the update
        if(dtSubstep<dt_min)then
          failedMinimumStep=.true.
          exit subSteps
        ! minimum OK -- try again
        else
          cycle substeps
        endif

      endif  ! if errors in prognostic update

      ! add balances to the total balances
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
      if(ixCasNrg(iGRU)/=integerMissing) sumBalance(ixCasNrg(iGRU),iGRU) = sumBalance(ixCasNrg(iGRU),iGRU) + dtSubstep*balance(ixCasNrg(iGRU),iGRU)
      if(ixVegNrg(iGRU)/=integerMissing) sumBalance(ixVegNrg(iGRU),iGRU) = sumBalance(ixVegNrg(iGRU),iGRU) + dtSubstep*balance(ixVegNrg(iGRU),iGRU)
      if(ixVegHyd(iGRU)/=integerMissing) sumBalance(ixVegHyd(iGRU),iGRU) = sumBalance(ixVegHyd(iGRU),iGRU) + dtSubstep*balance(ixVegHyd(iGRU),iGRU)
      if(ixAqWat(iGRU)/=integerMissing) sumBalance(ixAqWat(iGRU),iGRU) = sumBalance(ixAqWat(iGRU),iGRU) + dtSubstep*balance(ixAqWat(iGRU),iGRU)
    end do

      ! if(nSnowSoilNrg>0) then
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
        do ixLayer=1,nLayers(iGRU)
          if(ixSnowSoilNrg(ixLayer,iGRU)/=integerMissing) sumBalance(ixSnowSoilNrg(ixLayer,iGRU),iGRU) = sumBalance(ixSnowSoilNrg(ixLayer,iGRU),iGRU) + dtSubstep*balance(ixSnowSoilNrg(ixLayer,iGRU),iGRU)
        end do
      end do
      ! endif

      ! if(nSnowSoilHyd>0) then
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
        do ixLayer=1,nLayers(iGRU)
          if(ixSnowSoilHyd(ixLayer,iGRU)/=integerMissing) sumBalance(ixSnowSoilHyd(ixLayer,iGRU),iGRU) = sumBalance(ixSnowSoilHyd(ixLayer,iGRU),iGRU) + dtSubstep*balance(ixSnowSoilHyd(ixLayer,iGRU),iGRU)
        end do
      end do
      ! endif

      ! get the total energy fluxes (modified in updateProg), have to do differently
      if(nrgFluxModified)then
        associate(scalarCanopyEvaporation_temp => flux_temp%scalarCanopyEvaporation, &
          scalarLatHeatCanopyEvap_temp => flux_temp%scalarLatHeatCanopyEvap, &
          scalarSenHeatCanopy_temp => flux_temp%scalarSenHeatCanopy)
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
        sumCanopyEvaporation(iGRU)  = sumCanopyEvaporation(iGRU)  + dtSubstep*scalarCanopyEvaporation_temp(iGRU)  ! canopy evaporation/condensation (kg m-2 s-1)
        sumLatHeatCanopyEvap(iGRU)  = sumLatHeatCanopyEvap(iGRU)  + dtSubstep*scalarLatHeatCanopyEvap_temp(iGRU)  ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
        sumSenHeatCanopy(iGRU)      = sumSenHeatCanopy(iGRU)      + dtSubstep*scalarSenHeatCanopy_temp(iGRU)      ! sensible heat flux from the canopy to the canopy air space (W m-2)
          end do
          end associate
      else
        associate(scalarCanopyEvaporation_data => flux_data%scalarCanopyEvaporation, &
          scalarLatHeatCanopyEvap_data => flux_data%scalarLatHeatCanopyEvap, &
          scalarSenHeatCanopy_data => flux_data%scalarSenHeatCanopy, &
          scalarCanopyEvaporation_temp => flux_temp%scalarCanopyEvaporation, &
          scalarLatHeatCanopyEvap_temp => flux_temp%scalarLatHeatCanopyEvap, &
          scalarSenHeatCanopy_temp => flux_temp%scalarSenHeatCanopy)
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
            if (ixVegNrg(iGRU)/=integerMissing) then
              sumCanopyEvaporation(iGRU)  = sumCanopyEvaporation(iGRU)  + dtSubstep*scalarCanopyEvaporation_temp(iGRU)  ! canopy evaporation/condensation (kg m-2 s-1)
              sumLatHeatCanopyEvap(iGRU)  = sumLatHeatCanopyEvap(iGRU)  + dtSubstep*scalarLatHeatCanopyEvap_temp(iGRU)  ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
              sumSenHeatCanopy(iGRU)      = sumSenHeatCanopy(iGRU)      + dtSubstep*scalarSenHeatCanopy_temp(iGRU)      ! sensible heat flux from the canopy to the canopy air space (W m-2)
            else      
        sumCanopyEvaporation(iGRU)  = sumCanopyEvaporation(iGRU)  + dtSubstep*scalarCanopyEvaporation_data(iGRU)  ! canopy evaporation/condensation (kg m-2 s-1)
        sumLatHeatCanopyEvap(iGRU)  = sumLatHeatCanopyEvap(iGRU)  + dtSubstep*scalarLatHeatCanopyEvap_data(iGRU)  ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
        sumSenHeatCanopy(iGRU)      = sumSenHeatCanopy(iGRU)      + dtSubstep*scalarSenHeatCanopy_data(iGRU)      ! sensible heat flux from the canopy to the canopy air space (W m-2)
            endif
          end do
        end associate
      endif  ! if energy fluxes were modified

      ! get the total soil compression
      if (nSoil>0) then
        associate(scalarSoilCompress => diag_data%scalarSoilCompress, &
          mLayerCompress => diag_data%mLayerCompress_m)
        ! scalar compression
        if(.not.scalarSolution .or. iStateSplit==nSoil) then
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
            sumSoilCompress(iGRU) = sumSoilCompress(iGRU) + dtSubstep*scalarSoilCompress(iGRU) ! total soil compression
          end do
        endif
        ! vector compression
        !$cuf kernel do(2) <<<*,*>>>
          do iGRU=1,nGRU
        do iSoil=1,nSoil
          if(ixSoilOnlyHyd(iSoil,iGRU)/=integerMissing)&
          sumLayerCompress(iSoil,iGRU) = sumLayerCompress(iSoil,iGRU) + dtSubstep*mLayerCompress(iSoil,iGRU) ! soil compression in layers
        end do
      end do
      end associate
      endif

      ! print progress
      if(globalPrintFlag)&
      write(*,'(a,1x,3(f13.2,1x))') 'updating: dtSubstep, dtSum, dt = ', dtSubstep, dtSum, dt

     ! increment fluxes
      dt_wght = dtSubstep/dt ! define weight applied to each sub-step
      call update_flux_mean(flux_mean, flux_temp, dt_wght)

      do iVar=1,size(flux_meta)
        if(count(fluxMask%var(iVar)%dat)>0) then

          ! ** no domain splitting
          ! if(count(ixLayerActive/=integerMissing)==nLayers)then
      !       flux_mean%var(iVar)%dat(:) = flux_mean%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:)*dt_wght
            fluxCount%var(iVar)%dat(:) = fluxCount%var(iVar)%dat(:) + 1

          ! ** domain splitting
          ! else
          !   ixMin=lbound(flux_data%var(iVar)%dat)
          !   ixMax=ubound(flux_data%var(iVar)%dat)
          !   do ixLayer=ixMin(1),ixMax(1)
          !     if(fluxMask%var(iVar)%dat(ixLayer)) then
          !       ! special case of the transpiration sink from soil layers: only computed for the top soil layer
          !       if(iVar==iLookFLUX%mLayerTranspire)then
          !         if(ixLayer==1) flux_mean%var(iVar)%dat(:) = flux_mean%var(iVar)%dat(:) + flux_temp%var(iVar)%dat(:)*dt_wght
          !       ! standard case
          !       else
          !         flux_mean%var(iVar)%dat(ixLayer) = flux_mean%var(iVar)%dat(ixLayer) + flux_temp%var(iVar)%dat(ixLayer)*dt_wght
          !       endif
          !       fluxCount%var(iVar)%dat(ixLayer) = fluxCount%var(iVar)%dat(ixLayer) + 1
          !     endif
          !   end do
          ! endif  ! (domain splitting)

        endif   ! (if the flux is desired)
      end do  ! (loop through fluxes)

      ! increment the number of substeps
      nSubsteps = nSubsteps + nSteps

      ! increment the sub-step legth
      dtSum = dtSum + dtSubstep

      ! check that we have completed the sub-step
      if(dtSum >= dt-verySmall)then
        failedMinimumStep=.false.
        exit subSteps
      endif

      ! adjust length of the sub-step (make sure that we don't exceed the step)
      dtSubstep = min(dt - dtSum, max(dtSubstep*dtMultiplier, dt_min) )

    end do substeps  ! time steps for variable-dependent sub-stepping
    ! NOTE: if we get to here then we are accepting then dtSum should dt

    ! save the fluxes as averages
    do iVar=1,size(flux_meta)

      if(count(fluxMask%var(iVar)%dat)>0) then
        select case(iVar)
        case(iLookFLUX%iLayerConductiveFlux); flux_data%iLayerConductiveFlux_m = flux_mean%iLayerConductiveFlux_m
        case(iLookFLUX%iLayerAdvectiveFlux); flux_data%iLayerAdvectiveFlux_m = flux_mean%iLayerAdvectiveFlux_m
        case(iLookFLUX%mLayerSatHydCond); flux_data%mLayerSatHydCond_m = flux_mean%mLayerSatHydCond_m
        case(iLookFLUX%mLayerSatHydCondMP); flux_data%mLayerSatHydCondMP_m = flux_mean%mLayerSatHydCondMP_m
        case(iLookFLUX%iLayerSatHydCond); flux_data%iLayerSatHydCond_m = flux_mean%iLayerSatHydCond_m
        case(iLookFLUX%scalarExfiltration); flux_data%scalarExfiltration = flux_mean%scalarExfiltration
        case(iLookFLUX%mLayerColumnOutflow); flux_data%mLayerColumnOutflow_m = flux_mean%mLayerColumnOutflow_m
        case(iLookFLUX%scalarSoilBaseflow); flux_data%scalarSoilBaseflow = flux_mean%scalarSoilBaseflow
        case(iLookFLUX%mLayerBaseflow); flux_data%mLayerBaseflow_m = flux_mean%mLayerBaseflow_m
        case(iLookFLUX%scalarTotalRunoff); flux_data%scalarTotalRunoff = flux_mean%scalarTotalRunoff
        case(iLookFLUX%scalarSurfaceRunoff); flux_data%scalarSurfaceRunoff = flux_mean%scalarSurfaceRunoff
        case(iLookFLUX%scalarSoilDrainage); flux_data%scalarSoilDrainage = flux_mean%scalarSoilDrainage
        case(iLookFLUX%scalarSnowfall); flux_data%scalarSnowfall = flux_mean%scalarSnowfall
        case(iLookFLUX%scalarThroughfallSnow); flux_data%scalarThroughfallSnow = flux_mean%scalarThroughfallSnow
        case(iLookFLUX%scalarCanopySunlitPAR); flux_data%scalarCanopySunlitPAR = flux_mean%scalarCanopySunlitPAR
        case(iLookFLUX%scalarCanopyShadedPAR); flux_data%scalarCanopyShadedPAR = flux_mean%scalarCanopyShadedPAR
        case(iLookFLUX%scalarCanopyAbsorbedSolar); flux_data%scalarCanopyAbsorbedSolar = flux_mean%scalarCanopyAbsorbedSolar
        case(iLookFLUX%scalarGroundAbsorbedSolar); flux_data%scalarGroundAbsorbedSolar = flux_mean%scalarGroundAbsorbedSolar
        case(iLookFLUX%scalarLWRadCanopy); flux_data%scalarLWRadCanopy = flux_mean%scalarLWRadCanopy
        case(iLookFLUX%scalarLWRadGround); flux_data%scalarLWRadGround = flux_mean%scalarLWRadGround
        case(iLookFLUX%scalarLWRadUbound2Canopy); flux_data%scalarLWRadUbound2Canopy = flux_mean%scalarLWRadUbound2Canopy
        case(iLookFLUX%scalarLWRadUbound2Ground); flux_data%scalarLWRadUbound2Ground = flux_mean%scalarLWRadUbound2Ground
        case(iLookFLUX%scalarLWRadUbound2Ubound); flux_data%scalarLWRadUbound2Ubound = flux_mean%scalarLWRadUbound2Ubound
        case(iLookFLUX%scalarLWRadCanopy2Ubound); flux_data%scalarLWRadCanopy2Ubound = flux_mean%scalarLWRadCanopy2Ubound
        case(iLookFLUX%scalarLWRadCanopy2Ground); flux_data%scalarLWRadCanopy2Ground = flux_mean%scalarLWRadCanopy2Ground
        case(iLookFLUX%scalarLWRadCanopy2Canopy); flux_data%scalarLWRadCanopy2Canopy = flux_mean%scalarLWRadCanopy2Canopy
        case(iLookFLUX%scalarLWRadGround2Ubound); flux_data%scalarLWRadGround2Ubound = flux_mean%scalarLWRadGround2Ubound
        case(iLookFLUX%scalarLWRadGround2Canopy); flux_data%scalarLWRadGround2Canopy = flux_mean%scalarLWRadGround2Canopy
        case(iLookFLUX%scalarLWNetCanopy); flux_data%scalarLWNetCanopy = flux_mean%scalarLWNetCanopy
        case(iLookFLUX%scalarLWNetGround); flux_data%scalarLWNetGround = flux_mean%scalarLWNetGround
        case(iLookFLUX%scalarLWNetUbound); flux_data%scalarLWNetUbound = flux_mean%scalarLWNetUbound
        case(iLookFLUX%scalarEddyDiffusCanopyTop); flux_data%scalarEddyDiffusCanopyTop = flux_mean%scalarEddyDiffusCanopyTop
        case(iLookFLUX%scalarFrictionVelocity); flux_data%scalarFrictionVelocity = flux_mean%scalarFrictionVelocity
        case(iLookFLUX%scalarWindspdCanopyTop); flux_data%scalarWindspdCanopyTop = flux_mean%scalarWindspdCanopyTop
        case(iLookFLUX%scalarWindspdCanopyBottom); flux_data%scalarWindspdCanopyBottom = flux_mean%scalarWindspdCanopyBottom
        case(iLookFLUX%scalarLeafResistance); flux_data%scalarLeafResistance = flux_mean%scalarLeafResistance
        case(iLookFLUX%scalarGroundResistance); flux_data%scalarGroundResistance = flux_mean%scalarGroundResistance
        case(iLookFLUX%scalarCanopyResistance); flux_data%scalarCanopyResistance = flux_mean%scalarCanopyResistance
        case(iLookFLUX%scalarSoilResistance); flux_data%scalarSoilResistance = flux_mean%scalarSoilResistance
        case(iLookFLUX%scalarStomResistSunlit); flux_data%scalarStomResistSunlit = flux_mean%scalarStomResistSunlit
        case(iLookFLUX%scalarStomResistShaded); flux_data%scalarStomResistShaded = flux_mean%scalarStomResistShaded
        case(iLookFLUX%scalarPhotosynthesisSunlit); flux_data%scalarPhotosynthesisSunlit = flux_mean%scalarPhotosynthesisSunlit
        case(iLookFLUX%scalarPhotosynthesisShaded); flux_data%scalarPhotosynthesisShaded = flux_mean%scalarPhotosynthesisShaded
        case(iLookFLUX%scalarSenHeatTotal); flux_data%scalarSenHeatTotal = flux_mean%scalarSenHeatTotal
        case(iLookFLUX%scalarSenHeatCanopy); flux_data%scalarSenHeatCanopy = flux_mean%scalarSenHeatCanopy
        case(iLookFLUX%scalarSenHeatGround); flux_data%scalarSenHeatGround = flux_mean%scalarSenHeatGround
        case(iLookFLUX%scalarLatHeatTotal); flux_data%scalarLatHeatTotal = flux_mean%scalarLatHeatTotal
        case(iLookFLUX%scalarLatHeatCanopyEvap); flux_data%scalarLatHeatCanopyEvap = flux_mean%scalarLatHeatCanopyEvap
        case(iLookFLUX%scalarLatHeatCanopyTrans); flux_data%scalarLatHeatCanopyTrans = flux_mean%scalarLatHeatCanopyTrans
        case(iLookFLUX%scalarLatHeatGround); flux_data%scalarLatHeatGround = flux_mean%scalarLatHeatGround
        case(iLookFLUX%scalarCanopyAdvectiveHeatFlux); flux_data%scalarCanopyAdvectiveHeatFlux = flux_mean%scalarCanopyAdvectiveHeatFlux
        case(iLookFLUX%scalarGroundAdvectiveHeatFlux); flux_data%scalarGroundAdvectiveHeatFlux = flux_mean%scalarGroundAdvectiveHeatFlux
        case(iLookFLUX%scalarCanopySublimation); flux_data%scalarCanopySublimation = flux_mean%scalarCanopySublimation
        case(iLookFLUX%scalarSnowSublimation); flux_data%scalarSnowSublimation = flux_mean%scalarSnowSublimation
        case(iLookFLUX%scalarCanopyTranspiration); flux_data%scalarCanopyTranspiration = flux_mean%scalarCanopyTranspiration
        case(iLookFLUX%scalarCanopyEvaporation); flux_data%scalarCanopyEvaporation = flux_mean%scalarCanopyEvaporation
        case(iLookFLUX%scalarGroundEvaporation); flux_data%scalarGroundEvaporation = flux_mean%scalarGroundEvaporation
        case(iLookFLUX%scalarTotalET); flux_data%scalarTotalET = flux_mean%scalarTotalET
        case(iLookFLUX%scalarNetRadiation); flux_data%scalarNetRadiation = flux_mean%scalarNetRadiation
        case(iLookFLUX%scalarGroundNetNrgFlux); flux_data%scalarGroundNetNrgFlux = flux_mean%scalarGroundNetNrgFlux
        case(iLookFLUX%iLayerLiqFluxSnow); flux_data%iLayerLiqFluxSnow_m = flux_mean%iLayerLiqFluxSnow_m
        case(iLookFLUX%iLayerLiqFluxSoil); flux_data%iLayerLiqFluxSoil_m = flux_mean%iLayerLiqFluxSoil_m
        case(iLookFLUX%iLayerNrgFLux); flux_data%iLayerNrgFLux_m = flux_mean%iLayerNrgFLux_m
        case(iLookFLUX%scalarThroughfallRain); flux_data%scalarThroughfallRain = flux_mean%scalarThroughfallRain
        case(iLookFLUX%scalarCanopyLiqDrainage); flux_data%scalarCanopyLiqDrainage = flux_mean%scalarCanopyLiqDrainage
        case(iLookFLUX%scalarRainfall); flux_data%scalarRainfall = flux_mean%scalarRainfall
        case(iLookFLUX%scalarCanopyNetLiqFlux); flux_data%scalarCanopyNetLiqFlux = flux_mean%scalarCanopyNetLiqFlux
        case(iLookFLUX%scalarAquiferTranspire); flux_data%scalarAquiferTranspire = flux_mean%scalarAquiferTranspire
        case(iLookFLUX%scalarAquiferRecharge); flux_data%scalarAquiferRecharge = flux_mean%scalarAquiferRecharge
        case(iLookFLUX%scalarAquiferBaseflow); flux_data%scalarAquiferBaseflow = flux_mean%scalarAquiferBaseflow
        case(iLookFLUX%scalarRainPlusMelt); flux_data%scalarRainPlusMelt = flux_mean%scalarRainPlusMelt
        case(iLookFLUX%scalarMaxInfilRate); flux_data%scalarMaxInfilRate = flux_mean%scalarMaxInfilRate
        case(iLookFLUX%mLayerTranspire); flux_data%mLayerTranspire_m = flux_mean%mLayerTranspire_m
        case(iLookFLUX%mLayerHydCond); flux_data%mLayerHydCond_m = flux_mean%mLayerHydCond_m
        case(iLookFLUX%scalarInfiltration); flux_data%scalarInfiltration = flux_mean%scalarInfiltration
        case(iLookFLUX%mLayerLiqFluxSoil); flux_data%mLayerLiqFluxSoil_m = flux_mean%mLayerLiqFluxSoil_m
        case(iLookFLUX%mLayerLiqFluxSnow); flux_data%mLayerLiqFluxSnow_m = flux_mean%mLayerLiqFluxSnow_m
        case(iLookFLUX%scalarSnowDrainage); flux_data%scalarSnowDrainage = flux_mean%scalarSnowDrainage
        case(iLookFLUX%mLayerNrgFlux); flux_data%mLayerNrgFlux_m = flux_mean%mLayerNrgFlux_m
        case(iLookFLUX%scalarCanairNetNrgFlux); flux_data%scalarCanairNetNrgFlux = flux_mean%scalarCanairNetNrgFlux
        case(iLookFLUX%scalarCanopyNetNrgFlux); flux_data%scalarCanopyNetNrgFlux = flux_mean%scalarCanopyNetNrgFlux
        case(iLookFLUX%mLayerColumnInflow); flux_data%mLayerColumnInflow = flux_mean%mLayerColumnInflow
        end select
      endif

    enddo

    ! save the energy fluxes as averages
    associate(flux_scalarCanopyEvaporation => flux_data%scalarCanopyEvaporation, &
      flux_scalarLatHeatCanopyEvap => flux_data%scalarLatHeatCanopyEvap, &
      flux_scalarSenHeatCanopy => flux_data%scalarSenHeatCanopy)
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
    flux_scalarCanopyEvaporation(iGRU) = sumCanopyEvaporation(iGRU) /dtSum      ! canopy evaporation/condensation (kg m-2 s-1)
    flux_scalarLatHeatCanopyEvap(iGRU) = sumLatHeatCanopyEvap(iGRU) /dtSum      ! latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
    flux_scalarSenHeatCanopy(iGRU)     = sumSenHeatCanopy(iGRU)     /dtSum      ! sensible heat flux from the canopy to the canopy air space (W m-2)
      end do
      end associate
      err = cudaDeviceSynchronize()
  
    ! save the soil compression diagnostics as averages
      associate(diag_scalarSoilCompress => diag_data%scalarSoilCompress, &
        diag_mLayerCompress => diag_data%mLayerCompress_m)
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
    diag_scalarSoilCompress(iGRU) = sumSoilCompress(iGRU)/dtSum
        end do
        !$cuf kernel do(2) <<<*,*>>>
    do iGRU=1,nGRU
    do iSoil=1,nSoil
      if(ixSoilOnlyHyd(iSoil,iGRU)/=integerMissing)&
      diag_mLayerCompress(iSoil,iGRU) = sumLayerCompress(iSoil,iGRU)/dtSum
    end do
  end do
  end associate
    deallocate(sumLayerCompress)

    ! save the balance diagnostics as averages
    associate(balanceCasNrg => diag_data%balanceCasNrg, &
      balanceVegNrg => diag_data%balanceVegNrg, &
      balanceVegMass => diag_data%balanceVegMass, &
      balanceAqMass => diag_data%balanceAqMass)
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
    if(ixCasNrg(iGRU)/=integerMissing) balanceCasNrg(iGRU) = sumBalance(ixCasNrg(iGRU),iGRU)/dtSum
    if(ixVegNrg(iGRU)/=integerMissing) balanceVegNrg(iGRU) = sumBalance(ixVegNrg(iGRU),iGRU)/dtSum
    if(ixVegHyd(iGRU)/=integerMissing) balanceVegMass(iGRU) = sumBalance(ixVegHyd(iGRU),iGRU)/dtSum
    if(ixAqWat(iGRU)/=integerMissing) balanceAqMass(iGRU) = sumBalance(ixAqWat(iGRU),iGRU)/dtSum
      end do
    end associate


    ! if(nSnowSoilNrg>0) then
      associate(balanceLayerNrg => diag_data%balanceLayerNrg)
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
      do ixLayer=1,nLayers(iGRU)
        if (ixSnowSoilNrg(ixLayer,iGRU)/=integerMissing) balanceLayerNrg(ixLayer,iGRU) = sumBalance(ixSnowSoilNrg(ixLayer,iGRU),iGRU)/dtSum
      end do
    end do
      end associate
    ! endif
    ! if(nSnowSoilHyd>0) then
      associate(balanceLayerMass => diag_data%balanceLayerMass)
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
      do ixLayer=1,nLayers(iGRU)
        if (ixSnowSoilHyd(ixLayer,iGRU)/=integerMissing) balanceLayerMass(ixLayer,iGRU) = sumBalance(ixSnowSoilHyd(ixLayer,iGRU),iGRU)/dtSum
      end do
    end do
    end associate
    ! endif 

    call deallocate_device_flux_data(flux_temp)


    ! update error codes
    if (failedMinimumStep) then
      err=-20 ! negative = recoverable error
      message=trim(message)//'failed minimum step'
    end if
  ! end associate statements
  end associate globalVars

end subroutine varSubstep


! **********************************************************************************************************
! private subroutine updateProg: update prognostic variables
! **********************************************************************************************************
subroutine updateProg(dt,nSnow,nSoil,nLayers,nGRU, untappedMelt,stateVecTrial,stateVecPrime,                                           & ! input: states
                      doAdjustTemp,computeVegFlux,computMassBalance,computNrgBalance,computeEnthTemp,enthalpyStateVec,use_lookup,& ! input: model control
                      model_decisions,lookup_data,mpar_data,indx_data,flux_data,prog_data,diag_data,deriv_data,decisions,                  & ! input-output: data structures
                      fluxVec,resVec,balance,waterBalanceError,nrgFluxModified,err,message)                                        ! input-output: balances, flags, and error control
USE getVectorz_module,only:varExtract                              ! extract variables from the state vector
#ifdef SUNDIALS_ACTIVE
  USE updateVarsWithPrime_module,only:updateVarsWithPrime          ! update prognostic variables
  use device_data_types
  use initialize_device
#endif
  USE enthalpyTemp_module,only:enthTemp_or_enthalpy                ! add phase change terms to delta temperature component of enthalpy
  implicit none
  ! model control
  real(rkind)      ,intent(in)    :: dt                            ! time step (s)
  integer(i4b),device     ,intent(in)    :: nSnow(:)                         ! number of snow layers
  integer(i4b)     ,intent(in)    :: nSoil                         ! number of soil layers
  integer(i4b),device     ,intent(in)    :: nLayers(:)                       ! total number of layers
  integer(i4b) :: nGRU
  logical(lgt)     ,intent(in)    :: doAdjustTemp                  ! flag to indicate if we adjust the temperature
  logical(lgt)     ,intent(in),device    :: computeVegFlux(:)                ! flag to compute the vegetation flux
  real(rkind),device      ,intent(in)    :: untappedMelt(:,:)               ! un-tapped melt energy (J m-3 s-1)
  real(rkind),device      ,intent(in)    :: stateVecTrial(:,:)              ! trial state vector (mixed units)
  real(rkind),device      ,intent(in)    :: stateVecPrime(:,:)              ! trial state vector (mixed units)
  logical(lgt)     ,intent(in)    :: computMassBalance             ! flag to check the mass balance
  logical(lgt)     ,intent(in)    :: computNrgBalance              ! flag to check the energy balance
  logical(lgt)     ,intent(in)    :: computeEnthTemp               ! flag to compute enthalpy
  logical(lgt)     ,intent(in)    :: enthalpyStateVec              ! flag if enthalpy is a state variable (ida)
  logical(lgt)     ,intent(in)    :: use_lookup                    ! flag to use the lookup table for soil enthalpy, otherwise use analytical solution
  ! data structures
  type(model_options),intent(in)  :: model_decisions(:)            ! model decisions
  type(zLookup),intent(in)        :: lookup_data                   ! lookup tables
  type(mpar_data_device),intent(in)    :: mpar_data                     ! model parameters
  type(indx_data_device),intent(in)    :: indx_data                     ! indices for a local HRU
  type(flux_data_device),intent(inout) :: flux_data                     ! model fluxes for a local HRU
  type(prog_data_device),intent(inout) :: prog_data                     ! prognostic variables for a local HRU
  type(diag_data_device),intent(inout) :: diag_data                     ! diagnostic variables for a local HRU
  type(deriv_data_device),intent(inout) :: deriv_data                    ! derivatives in model fluxes w.r.t. relevant state variables
  type(decisions_device),intent(in) :: decisions
  ! balances, flags, and error control
  real(rkind),device      ,intent(in)    :: fluxVec(:,:)                    ! flux vector (mixed units)
  real(qp),device         ,intent(in)    :: resVec(:,:)    ! NOTE: qp       ! residual vector
  real(rkind),device      ,intent(inout) :: balance(:,:)                    ! balance of energy per domain per second
  logical(lgt)     ,intent(out)   :: waterBalanceError             ! flag to denote that there is a water balance error
  logical(lgt)     ,intent(out)   :: nrgFluxModified               ! flag to denote that the energy fluxes were modified
  integer(i4b)     ,intent(out)   :: err                           ! error code
  character(*)     ,intent(out)   :: message                       ! error message
  ! ==================================================================================================================
  ! general
  integer(i4b)                    :: i                             ! indices
  integer(i4b)                    :: iState                        ! index of model state variable
  integer(i4b)                    :: ixSubset                      ! index within the state subset
  integer(i4b)                    :: ixFullVector                  ! index within full state vector
  integer(i4b)                    :: ixControlIndex                ! index within a given domain
  real(rkind)                     :: volMelt                       ! volumetric melt (kg m-3)
  real(rkind),parameter           :: verySmall=epsilon(1._rkind)   ! a very small number (deal with precision issues)
  real(rkind)                     :: verySmall_veg                 ! precision needs to vary based on set canopy water tolerance for IDA
  real(rkind)                     :: verySmall_snow                ! precision needs to vary based on set snow water tolerance for IDA
  ! mass balance
  real(rkind)                     :: canopyBalance0,canopyBalance1 ! canopy storage at start/end of time step
  real(rkind)                     :: soilBalance0,soilBalance1     ! soil storage at start/end of time step
  real(rkind)                     :: vertFlux                      ! change in storage due to vertical fluxes
  real(rkind)                     :: tranSink,baseSink,compSink    ! change in storage due to sink terms
  real(rkind)                     :: liqError                      ! water balance error
  real(rkind)                     :: fluxNet                       ! net water fluxes (kg m-2 s-1)
  real(rkind)                     :: superflousWat                 ! superflous water used for evaporation (kg m-2 s-1)
  real(rkind)                     :: superflousNrg                 ! superflous energy that cannot be used for evaporation (W m-2 [J m-2 s-1])
  character(LEN=256)              :: cmessage                      ! error message of downwind routine
  ! trial state variables
  real(rkind),device,dimension(nGRU)                     :: scalarCanairTempTrial         ! trial value for temperature of the canopy air space (K)
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyTempTrial         ! trial value for temperature of the vegetation canopy (K)
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyWatTrial          ! trial value for liquid water storage in the canopy (kg m-2)
  real(rkind),device,dimension(maxSnowLayers+nSoil,nGRU)  :: mLayerTempTrial               ! trial vector of temperature of layers in the snow and soil domains (K)
  real(rkind),device,dimension(maxSnowLayers+nSoil,nGRU)  :: mLayerVolFracWatTrial         ! trial vector of volumetric fraction of total water (-)
  real(rkind),device,dimension(nSoil,nGRU)    :: mLayerMatricHeadTrial         ! trial vector of total water matric potential (m)
  real(rkind),device,dimension(nSoil,nGRU)    :: mLayerMatricHeadLiqTrial      ! trial vector of liquid water matric potential (m)
  real(rkind),device,dimension(nGRU)                     :: scalarAquiferStorageTrial     ! trial value for storage of water in the aquifer (m)
  real(rkind),device,dimension(nGRU)                     :: scalarCanairEnthalpyTrial     ! trial value for enthalpy of the canopy air space (J m-3)
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyEnthTempTrial     ! trial value for temperature component of enthalpy of the vegetation canopy (J m-3)
  real(rkind),device,dimension(nSoil+maxSnowLayers,nGRU)  :: mLayerEnthTempTrial           ! trial vector of temperature component of enthalpy of snow + soil (J m-3)
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyEnthalpyTrial     ! trial value for enthalpy of the vegetation canopy (J m-3)
  real(rkind),device,dimension(maxSnowLayers+nSoil,nGRU)  :: mLayerEnthalpyTrial           ! trial vector of enthalpy of each snow and soil layer (J m-3)
  ! diagnostic variables
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyLiqTrial          ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyIceTrial          ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),device,dimension(maxSnowLayers+nSoil,nGRU)  :: mLayerVolFracLiqTrial         ! trial vector of volumetric fraction of liquid water (-)
  real(rkind),device,dimension(maxSnowLayers+nSoil,nGRU)  :: mLayerVolFracIceTrial         ! trial vector of volumetric fraction of ice (-)
  ! prime state variables
  real(rkind),device,dimension(nGRU)                     :: scalarCanairTempPrime         ! trial value for temperature of the canopy air space (K)
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyTempPrime         ! trial value for temperature of the vegetation canopy (K)
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyWatPrime          ! trial value for liquid water storage in the canopy (kg m-2)
  real(rkind),device,dimension(maxSnowLayers+nSoil,nGRU)  :: mLayerTempPrime               ! trial vector of temperature of layers in the snow and soil domains (K)
  real(rkind),device,dimension(maxSnowLayers+nSoil,nGRU)  :: mLayerVolFracWatPrime         ! trial vector of volumetric fraction of total water (-)
  real(rkind),device,dimension(nSoil,nGRU)    :: mLayerMatricHeadPrime         ! trial vector of total water matric potential (m)
  real(rkind),device,dimension(nSoil,nGRU)    :: mLayerMatricHeadLiqPrime      ! trial vector of liquid water matric potential (m)
  real(rkind),device,dimension(nGRU)                     :: scalarAquiferStoragePrime     ! trial value for storage of water in the aquifer (m)
  ! diagnostic prime or delta variables
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyLiqPrime          ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyIcePrime          ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyIceDelta          ! delta value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyHDelta            ! delta value for enthalpy of the vegetation canopy (J m-3)
  real(rkind),device,dimension(maxSnowLayers+nSoil,nGRU)  :: mLayerVolFracLiqPrime         ! trial vector of volumetric fraction of liquid water (-)
  real(rkind),device,dimension(maxSnowLayers+nSoil,nGRU)  :: mLayerVolFracIcePrime         ! trial vector of volumetric fraction of ice (-)
  ! real(rkind),dimension(nLayers)  :: mLayerVolFracIceDelta         ! delta vector volumetric fraction of ice of snow + soil (-)
  ! real(rkind),dimension(nLayers)  :: mLayerHDelta                  ! delta vector of enthalpy of snow+soil (J m-3)
  ! dummy state variables
  real(rkind),device,dimension(nGRU)                     :: scalarCanairNrgTrial        ! trial value for energy of the canopy air space
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyNrgTrial        ! trial value for energy of the vegetation canopy
  real(rkind),device,dimension(maxSnowLayers+nSoil,nGRU)  :: mLayerNrgTrial              ! trial vector of energy of each snow and soil layer
  real(rkind),device,dimension(nGRU)                     :: scalarCanairNrgPrime        ! prime value for energy of the canopy air space
  real(rkind),device,dimension(nGRU)                     :: scalarCanopyNrgPrime        ! prime value for energy of the vegetation canopy
  real(rkind),device,dimension(maxSnowLayers+nSoil,nGRU)  :: mLayerNrgPrime              ! prime vector of energy of each snow and soil layer
  logical(lgt),device :: enthalpyStateVec_d, doAdjustTemp_d


  integer(i4b) :: iGRU

  ! -------------------------------------------------------------------------------------------------------------------
  ! -------------------------------------------------------------------------------------------------------------------
  ! point to flux variables in the data structure
  associate(&
    ! model decisions
    ixNumericalMethod         => model_decisions(iLookDECISIONS%num_method)%iDecision       ,& ! intent(in):  [i4b] choice of numerical solver
    ! get indices for balances
    ixCasNrg                  => indx_data%ixCasNrg                  ,& ! intent(in)   : [i4b]    index of canopy air space energy state variable
    ixVegNrg                  => indx_data%ixVegNrg                  ,& ! intent(in)   : [i4b]    index of canopy energy state variable
    ixVegHyd                  => indx_data%ixVegHyd                  ,& ! intent(in)   : [i4b]    index of canopy hydrology state variable (mass)
    ! ixTopNrg                  => indx_data%ixTopNrg                  ,& ! intent(in)   : [i4b]    index of upper-most energy state in the snow+soil subdomain
    ! ixTopHyd                  => indx_data%ixTopHyd                  ,& ! intent(in)   : [i4b]    index of upper-most hydrology state in the snow+soil subdomain
    ixAqWat                   => indx_data%ixAqWat                   ,& ! intent(in)   : [i4b]    index of water storage in the aquifer
    ixSoilOnlyHyd_m             => indx_data%ixSoilOnlyHyd                ,& ! intent(in)   : [i4b(:)] index in the state subset for hydrology state variables in the soil domain
    ixSnowSoilNrg             => indx_data%ixSnowSoilNrg                ,& ! intent(in)   : [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    ixSnowSoilHyd             => indx_data%ixSnowSoilHyd                ,& ! intent(in)   : [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    nSnowSoilNrg              => indx_data%nLayers_d              ,& ! intent(in)   : [i4b]    number of energy state variables in the snow+soil domain
    nSnowSoilHyd              => indx_data%nLayers_d              ,& ! intent(in)   : [i4b]    number of hydrology state variables in the snow+soil domain
    ! get indices for the un-tapped melt
    ixNrgOnly                 => indx_data%ixNrgOnly                    ,& ! intent(in)   : [i4b(:)] list of indices for all energy states
    ixDomainType              => indx_data%ixDomainType                 ,& ! intent(in)   : [i4b(:)] indices defining the domain of the state (iname_veg, iname_snow, iname_soil)
    ixControlVolume           => indx_data%ixControlVolume              ,& ! intent(in)   : [i4b(:)] index of the control volume for different domains (veg, snow, soil)
    ixMapSubset2Full          => indx_data%ixMapSubset2Full             ,& ! intent(in)   : [i4b(:)] [state subset] list of indices of the full state vector in the state subset
    ! water fluxes
    scalarRainfall            => flux_data%scalarRainfall             ,& ! intent(in)   : [dp]     rainfall rate (kg m-2 s-1)
    scalarThroughfallRain     => flux_data%scalarThroughfallRain      ,& ! intent(in)   : [dp]     rain reaches ground without touching the canopy (kg m-2 s-1)
    scalarCanopyEvaporation   => flux_data%scalarCanopyEvaporation    ,& ! intent(in)   : [dp]     canopy evaporation/condensation (kg m-2 s-1)
    scalarCanopyLiqDrainage   => flux_data%scalarCanopyLiqDrainage    ,& ! intent(in)   : [dp]     drainage liquid water from vegetation canopy (kg m-2 s-1)
    iLayerLiqFluxSoil         => flux_data%iLayerLiqFluxSoil_m             ,& ! intent(in)   : [dp(0:)] vertical liquid water flux at soil layer interfaces (-)
    iLayerNrgFlux             => flux_data%iLayerNrgFlux_m                 ,& ! intent(in)   :
    mLayerNrgFlux             => flux_data%mLayerNrgFlux_m                 ,& ! intent(out)  : [dp]     net energy flux for each layer within the snow+soil domain (J m-3 s-1)
    mLayerTranspire           => flux_data%mLayerTranspire_m               ,& ! intent(in)   : [dp(:)]  transpiration loss from each soil layer (m s-1)
    mLayerBaseflow            => flux_data%mLayerBaseflow_m                ,& ! intent(in)   : [dp(:)]  baseflow from each soil layer (m s-1)
    mLayerCompress            => diag_data%mLayerCompress_m                ,& ! intent(in)   : [dp(:)]  change in storage associated with compression of the soil matrix (-)
    ! energy fluxes
    scalarLatHeatCanopyEvap   => flux_data%scalarLatHeatCanopyEvap    ,& ! intent(in)   : [dp]     latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
    scalarSenHeatCanopy       => flux_data%scalarSenHeatCanopy        ,& ! intent(in)   : [dp]     sensible heat flux from the canopy to the canopy air space (W m-2)
    ! domain depth
    canopyDepth               => diag_data%scalarCanopyDepth          ,& ! intent(in)   : [dp   ]  canopy depth (m)
    mLayerDepth               => prog_data%mLayerDepth                   ,& ! intent(in)   : [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! model state variables (vegetation canopy)
    scalarCanairTemp          => prog_data%scalarCanairTemp           ,& ! intent(inout): [dp]     temperature of the canopy air space (K)
    scalarCanopyTemp          => prog_data%scalarCanopyTemp           ,& ! intent(inout): [dp]     temperature of the vegetation canopy (K)
    scalarCanopyIce           => prog_data%scalarCanopyIce            ,& ! intent(inout): [dp]     mass of ice on the vegetation canopy (kg m-2)
    scalarCanopyLiq           => prog_data%scalarCanopyLiq            ,& ! intent(inout): [dp]     mass of liquid water on the vegetation canopy (kg m-2)
    scalarCanopyWat           => prog_data%scalarCanopyWat            ,& ! intent(inout): [dp]     mass of total water on the vegetation canopy (kg m-2)
    ! model state variables (snow and soil domains)
    mLayerTemp                => prog_data%mLayerTemp                    ,& ! intent(inout): [dp(:)]  temperature of each snow/soil layer (K)
    mLayerVolFracIce          => prog_data%mLayerVolFracIce              ,& ! intent(inout): [dp(:)]  volumetric fraction of ice (-)
    mLayerVolFracLiq          => prog_data%mLayerVolFracLiq              ,& ! intent(inout): [dp(:)]  volumetric fraction of liquid water (-)
    mLayerVolFracWat          => prog_data%mLayerVolFracWat              ,& ! intent(inout): [dp(:)]  volumetric fraction of total water (-)
    mLayerMatricHead          => prog_data%mLayerMatricHead              ,& ! intent(inout): [dp(:)]  matric head (m)
    mLayerMatricHeadLiq       => diag_data%mLayerMatricHeadLiq           ,& ! intent(inout): [dp(:)]  matric potential of liquid water (m)
    ! enthalpy
    scalarCanairEnthalpy      => prog_data%scalarCanairEnthalpy       ,& ! intent(inout): [dp]     enthalpy of the canopy air space (J m-3)
    scalarCanopyEnthalpy      => prog_data%scalarCanopyEnthalpy       ,& ! intent(inout): [dp]     enthalpy of the vegetation canopy (J m-3)
    scalarCanopyEnthTemp      => diag_data%scalarCanopyEnthTemp       ,& ! intent(inout): [dp]     temperature component of enthalpy of the vegetation canopy (J m-3)
    mLayerEnthalpy            => prog_data%mLayerEnthalpy                ,& ! intent(inout): [dp(:)]  enthalpy of the snow+soil layers (J m-3)
    mLayerEnthTemp            => diag_data%mLayerEnthTemp                ,& ! intent(inout): [dp(:)]  temperature component of enthalpy of the snow+soil layers (J m-3)
    ! model state variables (aquifer)
    scalarAquiferStorage      => prog_data%scalarAquiferStorage       & ! intent(inout): [dp(:)]  storage of water in the aquifer (m)
    ! error tolerance
    ! absConvTol_liquid         => mpar_data%absConvTol_liquid          & ! intent(in)   : [dp]     absolute convergence tolerance for vol frac liq water (-)
    ) ! associating flux variables in the data structure
    ! -------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='updateProg/'
    enthalpyStateVec_d = enthalpyStateVec
    doAdjustTemp_d = doAdjustTemp


    ! initialize flags for water balance error and energy flux modification
    waterBalanceError=.false.
    nrgFluxModified = .false.

    ! get storage at the start of the step
    ! canopyBalance0 = merge(scalarCanopyLiq + scalarCanopyIce, realMissing, computeVegFlux)
    ! soilBalance0   = sum( (mLayerVolFracLiq(nSnow+1:nLayers) + mLayerVolFracIce(nSnow+1:nLayers)  )*mLayerDepth(nSnow+1:nLayers) )

    ! -----
    ! * update states...
    ! ------------------

    ! initialize to state variable from the last update
    scalarCanairTempTrial     = scalarCanairTemp
    scalarCanairEnthalpyTrial = scalarCanairEnthalpy
    scalarCanopyTempTrial     = scalarCanopyTemp
    scalarCanopyEnthalpyTrial = scalarCanopyEnthalpy
    scalarCanopyEnthTempTrial = scalarCanopyEnthTemp
    scalarCanopyWatTrial      = scalarCanopyWat
    scalarCanopyLiqTrial      = scalarCanopyLiq
    scalarCanopyIceTrial      = scalarCanopyIce
    mLayerTempTrial           = mLayerTemp
    mLayerEnthalpyTrial       = mLayerEnthalpy
    mLayerEnthTempTrial       = mLayerEnthTemp
    mLayerVolFracWatTrial     = mLayerVolFracWat
    mLayerVolFracLiqTrial     = mLayerVolFracLiq
    mLayerVolFracIceTrial     = mLayerVolFracIce
    mLayerMatricHeadTrial     = mLayerMatricHead
    mLayerMatricHeadLiqTrial  = mLayerMatricHeadLiq
    scalarAquiferStorageTrial = scalarAquiferStorage

    if(enthalpyStateVec)then ! use state variable as enthalpy
      scalarCanairNrgTrial = scalarCanairEnthalpy
      scalarCanopyNrgTrial = realMissing ! currently not splitting in ida so no need to update
      mLayerNrgTrial       = realMissing ! currently not splitting in ida so no need to update
    else
      scalarCanairNrgTrial = scalarCanairTemp
      scalarCanopyNrgTrial = scalarCanopyTemp
      mLayerNrgTrial       = mLayerTemp
    endif
      
    ! extract states from the state vector
    call varExtract(&
                    ! input
                    stateVecTrial,             & ! intent(in):    model state vector (mixed units)
                    indx_data,                 & ! intent(in):    indices defining model states and layers
                    nGRU, &
                    ! output: variables for the vegetation canopy
                    scalarCanairNrgTrial,      & ! intent(inout): trial value of energy of the canopy air space, temperature (K) or enthalpy (J m-3)
                    scalarCanopyNrgTrial,      & ! intent(inout): trial value of energy of the vegetation canopy, temperature (K) or enthalpy (J m-3)
                    scalarCanopyWatTrial,      & ! intent(inout): trial value of canopy total water (kg m-2)
                    scalarCanopyLiqTrial,      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                    ! output: variables for the snow-soil domain
                    mLayerNrgTrial,            & ! intent(inout): trial vector of energy, temperature (K) or enthalpy (J m-3)
                    mLayerVolFracWatTrial,     & ! intent(inout): trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,     & ! intent(inout): trial vector of volumetric liquid water content (-)
                    mLayerMatricHeadTrial,     & ! intent(inout): trial vector of total water matric potential (m)
                    mLayerMatricHeadLiqTrial,  & ! intent(inout): trial vector of liquid water matric potential (m)
                    ! output: variables for the aquifer
                    scalarAquiferStorageTrial, & ! intent(inout): trial value of storage of water in the aquifer (m)
                    ! output: error control
                    err,cmessage)               ! intent(out):   error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)
  
    if(enthalpyStateVec)then ! use state variable as enthalpy
      scalarCanairEnthalpyTrial = scalarCanairNrgTrial
      scalarCanopyEnthalpyTrial = scalarCanopyNrgTrial
      mLayerEnthalpyTrial       = mLayerNrgTrial
    else
      scalarCanairTempTrial = scalarCanairNrgTrial
      scalarCanopyTempTrial = scalarCanopyNrgTrial
      mLayerTempTrial       = mLayerNrgTrial
    endif

    ! Placeholder: if we decide to use splitting, we need to pass all the previous values of the state variables
    scalarCanairNrgPrime      = realMissing
    scalarCanopyNrgPrime      = realMissing
    scalarCanopyWatPrime      = realMissing
    scalarCanopyLiqPrime      = realMissing
    scalarCanopyIcePrime      = realMissing
    mLayerNrgPrime            = realMissing
    mLayerVolFracWatPrime     = realMissing
    mLayerVolFracLiqPrime     = realMissing
    mLayerVolFracIcePrime     = realMissing
    mLayerMatricHeadPrime     = realMissing
    mLayerMatricHeadLiqPrime  = realMissing
    scalarAquiferStoragePrime = realMissing

    ! set the default precision for the very small number
    verySmall_veg  = verySmall*2._rkind
    verySmall_snow = verySmall*2._rkind

        ! IDA precision needs to vary based on set tolerances
        verySmall_veg = mpar_data%absTolWatVeg*2._rkind
        verySmall_snow = mpar_data%absTolWatSnow*2._rkind

        ! extract the derivatives from the state vector
        call varExtract(&
                  ! input
                  stateVecPrime,             & ! intent(in):    derivative of model state vector (mixed units)
                  indx_data,                 & ! intent(in):    indices defining model states and layers
                  nGRU, &
                  ! output: variables for the vegetation canopy
                  scalarCanairNrgPrime,      & ! intent(inout): derivative of energy of the canopy air space, temperature (K s-1) or enthalpy (W m-3)
                  scalarCanopyNrgPrime,      & ! intent(inout): derivative of energy of the vegetation canopy, temperature (K s-1) or enthalpy (W m-3)
                  scalarCanopyWatPrime,      & ! intent(inout): derivative of canopy total water (kg m-2 s-1)
                  scalarCanopyLiqPrime,      & ! intent(inout): derivative of canopy liquid water (kg m-2 s-1)
                  ! output: variables for the snow-soil domain
                  mLayerNrgPrime,            & ! intent(inout): derivative of energy of each snow and soil layer, temperature (K s-1) or enthalpy (W m-3)
                  mLayerVolFracWatPrime,     & ! intent(inout):   derivative of volumetric total water content (-)
                  mLayerVolFracLiqPrime,     & ! intent(inout):   derivative of volumetric liquid water content (-)
                  mLayerMatricHeadPrime,     & ! intent(inout):   derivative of total water matric potential (m)
                  mLayerMatricHeadLiqPrime,  & ! intent(inout):   derivative of liquid water matric potential (m)
                  ! output: variables for the aquifer
                  scalarAquiferStoragePrime, & ! intent(inout):   derivative of storage of water in the aquifer (m)
                  ! output: error control
                  err,cmessage)               ! intent(out):   error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

        if(enthalpyStateVec)then ! use state variable as enthalpy, need to compute temperature
          ! do not use these variables
          scalarCanairTempPrime = realMissing
          scalarCanopyTempPrime = realMissing
          mLayerTempPrime       = realMissing
        else ! use state variable as temperature
          scalarCanairTempPrime = scalarCanairNrgPrime
          scalarCanopyTempPrime = scalarCanopyNrgPrime
          mLayerTempPrime       = mLayerNrgPrime   
        endif !(choice of how conservation of energy is implemented)
    
        ! update diagnostic variables
        call updateVarsWithPrime(&
                    ! input
                    decisions%nrgConserv,                 & ! intent(in):    flag if enthalpy is used as state variable
                    use_lookup,                       & ! intent(in):    flag to use the lookup table for soil enthalpy
                    decisions%false,                          & ! intent(in):    logical flag if computing for Jacobian update
                    doAdjustTemp_d,                     & ! intent(in):    logical flag to adjust temperature to account for the energy used in melt+freeze
                    nGRU, &
                    mpar_data,                        & ! intent(in):    model parameters for a local HRU
                    indx_data,                        & ! intent(in):    indices defining model states and layers
                    prog_data,                        & ! intent(in):    model prognostic variables for a local HRU
                    diag_data,                        & ! intent(inout): model diagnostic variables for a local HRU
                    deriv_data,                       & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                    lookup_data,                      & ! intent(in):    lookup table data structure
                    ! input: enthalpy state variables  
                    scalarCanairEnthalpyTrial,        & ! intent(in):    trial value for enthalpy of the canopy air space (J m-3)
                    scalarCanopyEnthalpyTrial,        & ! intent(in):    trial value for enthalpy of the vegetation canopy (J m-3)
                    mLayerEnthalpyTrial,              & ! intent(in):    trial vector of enthalpy of each snow+soil layer (J m-3)                      
                    ! output: variables for the vegetation canopy
                    scalarCanairTempTrial,            & ! intent(inout): trial value of canopy air space temperature (K)
                    scalarCanopyTempTrial,            & ! intent(inout): trial value of canopy temperature (K)
                    scalarCanopyWatTrial,             & ! intent(inout): trial value of canopy total water (kg m-2)
                    scalarCanopyLiqTrial,             & ! intent(inout): trial value of canopy liquid water (kg m-2)
                    scalarCanopyIceTrial,             & ! intent(inout): trial value of canopy ice content (kg m-2)
                    scalarCanopyTempPrime,            & ! intent(inout): trial value of canopy temperature (K)
                    scalarCanopyWatPrime,             & ! intent(inout): trial value of canopy total water (kg m-2)
                    scalarCanopyLiqPrime,             & ! intent(inout): trial value of canopy liquid water (kg m-2)
                    scalarCanopyIcePrime,             & ! intent(inout): trial value of canopy ice content (kg m-2)
                    ! output: variables for the snow-soil domain
                    mLayerTempTrial,                  & ! intent(inout): trial vector of layer temperature (K)
                    mLayerVolFracWatTrial,            & ! intent(inout): trial vector of volumetric total water content (-)
                    mLayerVolFracLiqTrial,            & ! intent(inout): trial vector of volumetric liquid water content (-)
                    mLayerVolFracIceTrial,            & ! intent(inout): trial vector of volumetric ice water content (-)
                    mLayerMatricHeadTrial,            & ! intent(inout): trial vector of total water matric potential (m)
                    mLayerMatricHeadLiqTrial,         & ! intent(inout): trial vector of liquid water matric potential (m)
                    mLayerTempPrime,                  & ! intent(inout): Prime vector of layer temperature (K)
                    mLayerVolFracWatPrime,            & ! intent(inout): Prime vector of volumetric total water content (-)
                    mLayerVolFracLiqPrime,            & ! intent(inout): Prime vector of volumetric liquid water content (-)
                    mLayerVolFracIcePrime,            & ! intent(inout): Prime vector of volumetric ice water content (-)
                    mLayerMatricHeadPrime,            & ! intent(inout): Prime vector of total water matric potential (m)
                    mLayerMatricHeadLiqPrime,         & ! intent(inout): Prime vector of liquid water matric potential (m)
                    ! output: error control
                    err,cmessage)                       ! intent(out):   error control
                

    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if  ! (check for errors)

    if( .not. computNrgBalance)then
     ! if not checking energy balance set balance to missing
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        if(ixCasNrg(iGRU)/=integerMissing) balance(ixCasNrg(iGRU),iGRU) = realMissing
        if(ixVegNrg(iGRU)/=integerMissing) balance(ixVegNrg(iGRU),iGRU) = realMissing
      end do
      ! if(nSnowSoilNrg>0)then
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
        do i=1,nLayers(iGRU)
          if (ixSnowSoilNrg(i,iGRU)/=integerMissing) then
            balance(ixSnowSoilNrg(i,iGRU),iGRU) = realMissing
          endif
        enddo
      end do
      ! endif
    endif  ! if checking energy balance

    ! -----
    ! * check mass balance...
    ! -----------------------

    ! NOTE: currently this will only fail with kinsol solver, since mass balance is checked in the homegrown solver and not checked for ida solver
    !   Negative error code will mean step will be failed and retried with smaller step size
!     if(computMassBalance)then
! print*, computMassBalance
!       if(ixVegHyd/=integerMissing)then ! check for complete drainage

!         ! handle cases where fluxes empty the canopy
!         fluxNet = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
!         if(-fluxNet*dt > canopyBalance0)then

!           ! --> first add water
!           canopyBalance1 = canopyBalance0 + (scalarRainfall - scalarThroughfallRain)*dt

!           ! --> next, remove canopy evaporation -- put the unsatisfied evap into sensible heat
!           canopyBalance1 = canopyBalance1 + scalarCanopyEvaporation*dt
!           if(canopyBalance1 < 0._rkind)then
!             ! * get superfluous water and energy
!             superflousWat = -canopyBalance1/dt     ! kg m-2 s-1
!             superflousNrg = superflousWat*LH_vap   ! W m-2 (J m-2 s-1)
!             ! * update fluxes and states
!             canopyBalance1          = 0._rkind
!             scalarCanopyEvaporation = scalarCanopyEvaporation + superflousWat
!             scalarLatHeatCanopyEvap = scalarLatHeatCanopyEvap + superflousNrg
!             scalarSenHeatCanopy     = scalarSenHeatCanopy - superflousNrg
!           endif

!           ! --> next, remove canopy drainage
!           canopyBalance1 = canopyBalance1 -scalarCanopyLiqDrainage*dt
!           if(canopyBalance1 < 0._rkind)then
!             superflousWat           = -canopyBalance1/dt     ! kg m-2 s-1
!             canopyBalance1          = 0._rkind
!             scalarCanopyLiqDrainage = scalarCanopyLiqDrainage + superflousWat
!           endif

!           ! update the trial state
!           scalarCanopyWatTrial = canopyBalance1

!           ! set the modification flag
!           nrgFluxModified = .true.

!         else
!           canopyBalance1  = canopyBalance0 + fluxNet*dt
!           nrgFluxModified = .false.
!         endif  ! cases where fluxes empty the canopy
      
!       endif ! check for complete drainage

    ! else ! if not checking mass balance set balance to missing
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
      if(ixVegHyd(iGRU)/=integerMissing) balance(ixVegHyd(iGRU),iGRU) = realMissing
      if(ixAqWat(iGRU)/=integerMissing) balance(ixAqWat(iGRU),iGRU) = realMissing
    end do

      ! if(nSnowSoilHyd>0)then
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
        do i=1,nLayers(iGRU)
          if (ixSnowSoilHyd(i,iGRU)/=integerMissing) balance(ixSnowSoilHyd(i,iGRU),iGRU) = realMissing
        end do
      end do
      ! endif
    ! endif  ! if checking the mass balance
  
    ! -----
    ! * remove untapped melt energy... always 0 at the moment but if use should be in solved as affects state
    ! --------------------------------

    ! only work with energy state variables
    if(size(ixNrgOnly,1)>0)then  ! energy state variables exist

      ! loop through energy state variables
      !$cuf kernel do(2) <<<*,*>>>
      do iGRU=1,nGRU
      do iState=1,size(ixNrgOnly,1)

        ! get index of the control volume within the domain
        ixSubset       = ixNrgOnly(iState,iGRU)             ! index within the state subset
        ixFullVector   = ixMapSubset2Full(ixSubset,iGRU)    ! index within full state vector
        ixControlIndex = ixControlVolume(ixFullVector,iGRU) ! index within a given domain

        ! compute volumetric melt (kg m-3)
        volMelt = dt*untappedMelt(ixSubset,iGRU)/LH_fus  ! (kg m-3)

        ! update ice content
        select case( ixDomainType(ixFullVector,iGRU) )
          case(iname_cas);  cycle ! do nothing, since there is no snow stored in the canopy air space
          case(iname_veg);  scalarCanopyIceTrial(iGRU)                        = scalarCanopyIceTrial(iGRU)                        - volMelt*canopyDepth(iGRU)  ! (kg m-2)
          case(iname_snow); mLayerVolFracIceTrial(ixControlIndex,iGRU)       = mLayerVolFracIceTrial(ixControlIndex,iGRU)       - volMelt/iden_ice     ! (-)
          case(iname_soil); mLayerVolFracIceTrial(ixControlIndex+nSnow(iGRU),iGRU) = mLayerVolFracIceTrial(ixControlIndex+nSnow(iGRU),iGRU) - volMelt/iden_water   ! (-)
          ! case default; err=20; message=trim(message)//'unable to identify domain type [remove untapped melt energy]'; return
        end select

        ! update liquid water content
        select case( ixDomainType(ixFullVector,iGRU) )
          case(iname_cas);  cycle ! do nothing, since there is no snow stored in the canopy air space
          case(iname_veg);  scalarCanopyLiqTrial(iGRU)                        = scalarCanopyLiqTrial(iGRU)                        + volMelt*canopyDepth(iGRU)  ! (kg m-2)
          case(iname_snow); mLayerVolFracLiqTrial(ixControlIndex,iGRU)       = mLayerVolFracLiqTrial(ixControlIndex,iGRU)       + volMelt/iden_water   ! (-)
          case(iname_soil); mLayerVolFracLiqTrial(ixControlIndex+nSnow(iGRU),iGRU) = mLayerVolFracLiqTrial(ixControlIndex+nSnow(iGRU),iGRU) + volMelt/iden_water   ! (-)
          ! case default; err=20; message=trim(message)//'unable to identify domain type [remove untapped melt energy]'; return
        end select

      end do  ! looping through energy variables
    end do

      ! ========================================================================================================

      ! *** ice

      ! --> check if we removed too much water
      ! if(any(scalarCanopyIceTrial < 0._rkind)  .or. any(mLayerVolFracIceTrial_d < 0._rkind) )then

        ! **
        ! canopy within numerical precision
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
        if(scalarCanopyIceTrial(iGRU) < 0._rkind)then

          if(scalarCanopyIceTrial(iGRU) > -verySmall_veg)then
            scalarCanopyLiqTrial(iGRU) = scalarCanopyLiqTrial(iGRU) - scalarCanopyIceTrial(iGRU)
            scalarCanopyIceTrial(iGRU) = 0._rkind

          ! encountered an inconsistency: spit the dummy
          ! else
          !   print*, 'dt = ', dt
          !   print*, 'untappedMelt          = ', untappedMelt
          !   print*, 'untappedMelt*dt       = ', untappedMelt*dt
          !   print*, 'scalarCanopyiceTrial  = ', scalarCanopyIceTrial
          !   message=trim(message)//'melted more than the available water'
          !   err=20; return
          endif  ! (inconsistency)

        endif  ! if checking the canopy
      end do
        ! **
        ! snow+soil within numerical precision
      !$cuf kernel do(2) <<<*,*>>>
      do iGRU=1,nGRU
        do iState=1,size(mLayerVolFracIceTrial,1)

          ! snow layer within numerical precision
          if(mLayerVolFracIceTrial(iState,iGRU) < 0._rkind)then

            if(mLayerVolFracIceTrial(iState,iGRU) > -verySmall_snow)then
              mLayerVolFracLiqTrial(iState,iGRU) = mLayerVolFracLiqTrial(iState,iGRU) - mLayerVolFracIceTrial(iState,iGRU)
              mLayerVolFracIceTrial(iState,iGRU) = 0._rkind

            ! encountered an inconsistency: spit the dummy
            ! else
            !   print*, 'dt = ', dt
            !   print*, 'untappedMelt          = ', untappedMelt
            !   print*, 'untappedMelt*dt       = ', untappedMelt*dt
            !   print*, 'mLayerVolFracIceTrial = ', mLayerVolFracIceTrial
            !   message=trim(message)//'melted more than the available water'
            !   err=20; return
            endif  ! (inconsistency)

          endif  ! if checking a snow layer

        end do ! (looping through state variables)
      end do

      ! endif  ! (if we removed too much water)

      ! ========================================================================================================  
      ! *** liquid water

      ! --> check if we removed too much water
      ! if(scalarCanopyLiqTrial < 0._rkind  .or. any(mLayerVolFracLiqTrial < 0._rkind) )then

        ! **
        ! canopy within numerical precision
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        if(scalarCanopyLiqTrial(iGRU) < 0._rkind)then

          if(scalarCanopyLiqTrial(iGRU) > -verySmall_veg)then
            scalarCanopyIceTrial(iGRU) = scalarCanopyIceTrial(iGRU) - scalarCanopyLiqTrial(iGRU)
            scalarCanopyLiqTrial(iGRU) = 0._rkind

          ! encountered an inconsistency: spit the dummy
          ! else
          !   print*, 'dt = ', dt
          !   print*, 'untappedMelt          = ', untappedMelt
          !   print*, 'untappedMelt*dt       = ', untappedMelt*dt
          !   print*, 'scalarCanopyLiqTrial  = ', scalarCanopyLiqTrial
          !   message=trim(message)//'frozen more than the available water'
          !   err=20; return
          endif  ! (inconsistency)
        endif  ! checking the canopy
      end do

        ! **
        ! snow+soil within numerical precision
      !$cuf kernel do(2) <<<*,*>>>
      do iGRU=1,nGRU
        do iState=1,size(mLayerVolFracLiqTrial,1)

          ! snow layer within numerical precision
          if(mLayerVolFracLiqTrial(iState,iGRU) < 0._rkind)then

            if(mLayerVolFracLiqTrial(iState,iGRU) > -verySmall_snow)then
              mLayerVolFracIceTrial(iState,iGRU) = mLayerVolFracIceTrial(iState,iGRU) - mLayerVolFracLiqTrial(iState,iGRU)
              mLayerVolFracLiqTrial(iState,iGRU) = 0._rkind

            ! encountered an inconsistency: spit the dummy
            ! else
            !   print*, 'dt = ', dt
            !   print*, 'untappedMelt          = ', untappedMelt
            !   print*, 'untappedMelt*dt       = ', untappedMelt*dt
            !   print*, 'mLayerVolFracLiqTrial = ', mLayerVolFracLiqTrial
            !   message=trim(message)//'frozen more than the available water'
            !   err=20; return
            endif  ! (inconsistency)

          endif  ! checking a snow layer

        end do ! (looping through state variables)
      end do

      endif  ! (if we removed too much water)

    ! endif  ! (if energy state variables exist)

    ! -----
    ! * update enthalpy as a diagnostic variable... 
    !   if computeEnthTemp then enthTemp will change, if enthalpyStateVec then enthalpy will change
    ! --------------------------------
    scalarCanairEnthalpy = scalarCanairEnthalpyTrial ! equivalent to scalarCanairEnthTemp
    scalarCanopyEnthTemp = scalarCanopyEnthTempTrial
    scalarCanopyEnthalpy = scalarCanopyEnthalpyTrial
    mLayerEnthTemp       = mLayerEnthTempTrial
    mLayerEnthalpy       = mLayerEnthalpyTrial

    ! -----
    ! * update prognostic variables...
    ! --------------------------------
    ! update state variables for the vegetation canopy
    scalarCanairTemp    = scalarCanairTempTrial    ! trial value of canopy air temperature (K)
    scalarCanopyTemp    = scalarCanopyTempTrial    ! trial value of canopy temperature (K)
    scalarCanopyWat     = scalarCanopyWatTrial     ! trial value of canopy total water (kg m-2)
    scalarCanopyLiq     = scalarCanopyLiqTrial     ! trial value of canopy liquid water (kg m-2)
    scalarCanopyIce     = scalarCanopyIceTrial     ! trial value of canopy ice content (kg m-2)

    ! update state variables for the snow+soil domain
    mLayerTemp          = mLayerTempTrial          ! trial vector of layer temperature (K)
    mLayerVolFracWat    = mLayerVolFracWatTrial    ! trial vector of volumetric total water content (-)
    mLayerVolFracLiq    = mLayerVolFracLiqTrial    ! trial vector of volumetric liquid water content (-)
    mLayerVolFracIce    = mLayerVolFracIceTrial    ! trial vector of volumetric ice water content (-)
    mLayerMatricHead    = mLayerMatricHeadTrial    ! trial vector of matric head (m)
    mLayerMatricHeadLiq = mLayerMatricHeadLiqTrial ! trial vector of matric head (m)

    ! update state variables for the aquifer
    scalarAquiferStorage = scalarAquiferStorageTrial

    ! end associations to info in the data structures
  end associate

end subroutine updateProg
subroutine update_flux_mean(flux_mean,flux_temp,scale)
  use device_data_types
  type(flux_data_device) :: flux_mean, flux_temp
  real(rkind) :: scale
  integer(i4b) :: nGRU
  integer(i4b) :: nLayers,nSnow,nSoil
  nGRU = size(flux_mean%scalarCanopyResistance)
  nLayers = size(flux_mean%mLayerNrgFlux_m,1)
  nSnow = size(flux_mean%mLayerLiqFluxSnow_m,1)
  nSoil = size(flux_mean%mLayerLiqFluxSoil_m,1)
  call update_mean(flux_mean%scalarCanopyResistance, flux_temp%scalarCanopyResistance, scale, nGRU)
  call update_mean(flux_mean%scalarGroundResistance, flux_temp%scalarGroundResistance, scale, nGRU)
  call update_mean(flux_mean%scalarLeafResistance, flux_temp%scalarLeafResistance, scale, nGRU)
  call update_mean(flux_mean%scalarEddyDiffusCanopyTop, flux_temp%scalarEddyDiffusCanopyTop, scale, nGRU)
  call update_mean(flux_mean%scalarFrictionVelocity, flux_temp%scalarFrictionVelocity, scale, nGRU)
  call update_mean(flux_mean%scalarWindspdCanopyTop, flux_temp%scalarWindspdCanopyTop, scale, nGRU)
  call update_mean(flux_mean%scalarWindspdCanopyBottom, flux_temp%scalarWindspdCanopyBottom, scale, nGRU)
  call update_mean(flux_mean%scalarSenHeatTotal, flux_temp%scalarSenHeatTotal, scale, nGRU)
  call update_mean(flux_mean%scalarSenHeatCanopy, flux_temp%scalarSenHeatCanopy, scale, nGRU)
  call update_mean(flux_mean%scalarSenHeatGround, flux_temp%scalarSenHeatGround, scale, nGRU)
  call update_mean(flux_mean%scalarLatHeatTotal, flux_temp%scalarLatHeatTotal, scale, nGRU)
  call update_mean(flux_mean%scalarLatHeatCanopyEvap, flux_temp%scalarLatHeatCanopyEvap, scale, nGRU)
  call update_mean(flux_mean%scalarLatHeatCanopyTrans, flux_temp%scalarLatHeatCanopyTrans, scale, nGRU)
  call update_mean(flux_mean%scalarLatHeatGround, flux_temp%scalarLatHeatGround, scale, nGRU)
  call update_mean(flux_mean%scalarAquiferBaseflow, flux_temp%scalarAquiferBaseflow, scale, nGRU)
  call update_mean(flux_mean%scalarAquiferRecharge, flux_temp%scalarAquiferRecharge, scale, nGRU)
  call update_mean(flux_mean%scalarAquiferTranspire, flux_temp%scalarAquiferTranspire, scale, nGRU)
  call update_mean(flux_mean%scalarCanopyLiqDrainage, flux_temp%scalarCanopyLiqDrainage, scale, nGRU)
  call update_mean(flux_mean%scalarThroughfallRain, flux_temp%scalarThroughfallRain, scale, nGRU)
  call update_mean(flux_mean%scalarRainfall, flux_temp%scalarRainfall, scale, nGRU)
  call update_mean(flux_mean%scalarCanopyNetLiqFlux, flux_temp%scalarCanopyNetLiqFlux, scale, nGRU)
  call update_mean(flux_mean%scalarCanopyEvaporation, flux_temp%scalarCanopyEvaporation, scale, nGRU)
  call update_mean(flux_mean%scalarCanopyAdvectiveHeatFlux, flux_temp%scalarCanopyAdvectiveHeatFlux, scale, nGRU)
  call update_mean(flux_mean%scalarGroundAdvectiveHeatFlux, flux_temp%scalarGroundAdvectiveHeatFlux, scale, nGRU)
  call update_mean(flux_mean%scalarCanopySublimation, flux_temp%scalarCanopySublimation, scale, nGRU)
  call update_mean(flux_mean%scalarSnowSublimation, flux_temp%scalarSnowSublimation, scale, nGRU)
  call update_mean(flux_mean%scalarCanopyTranspiration, flux_temp%scalarCanopyTranspiration, scale, nGRU)
  call update_mean(flux_mean%scalarGroundEvaporation, flux_temp%scalarGroundEvaporation, scale, nGRU)
  call update_mean(flux_mean%scalarLWRadCanopy, flux_temp%scalarLWRadCanopy, scale, nGRU)
  call update_mean(flux_mean%scalarLWRadGround, flux_temp%scalarLWRadGround, scale, nGRU)
  call update_mean(flux_mean%scalarLWRadUbound2Canopy, flux_temp%scalarLWRadUbound2Canopy, scale, nGRU)
  call update_mean(flux_mean%scalarLWRadUbound2Ground, flux_temp%scalarLWRadUbound2Ground, scale, nGRU)
  call update_mean(flux_mean%scalarLWRadUbound2Ubound, flux_temp%scalarLWRadUbound2Ubound, scale, nGRU)
  call update_mean(flux_mean%scalarLWRadCanopy2Ubound, flux_temp%scalarLWRadCanopy2Ubound, scale, nGRU)
  call update_mean(flux_mean%scalarLWRadCanopy2Ground, flux_temp%scalarLWRadCanopy2Ground, scale, nGRU)
  call update_mean(flux_mean%scalarLWRadCanopy2Canopy, flux_temp%scalarLWRadCanopy2Canopy, scale, nGRU)
  call update_mean(flux_mean%scalarLWRadGround2Ubound, flux_temp%scalarLWRadGround2Ubound, scale, nGRU)
  call update_mean(flux_mean%scalarLWRadGround2Canopy, flux_temp%scalarLWRadGround2Canopy, scale, nGRU)
  call update_mean(flux_mean%scalarLWNetCanopy, flux_temp%scalarLWNetCanopy, scale, nGRU)
  call update_mean(flux_mean%scalarLWNetGround, flux_temp%scalarLWNetGround, scale, nGRU)
  call update_mean(flux_mean%scalarLWNetUbound, flux_temp%scalarLWNetUbound, scale, nGRU)
  call update_mean(flux_mean%scalarGroundNetNrgFlux, flux_temp%scalarGroundNetNrgFlux, scale, nGRU)
  call update_mean(flux_mean%scalarCanairNetNrgFlux, flux_temp%scalarCanairNetNrgFlux, scale, nGRU)
  call update_mean(flux_mean%scalarCanopyNetNrgFlux, flux_temp%scalarCanopyNetNrgFlux, scale, nGRU)
  call update_mean(flux_mean%scalarSnowfall, flux_temp%scalarSnowfall, scale, nGRU)
  call update_mean(flux_mean%scalarThroughfallSnow, flux_temp%scalarThroughfallSnow, scale, nGRU)
  call update_mean(flux_mean%scalarCanopyAbsorbedSolar, flux_temp%scalarCanopyAbsorbedSolar, scale, nGRU)
  call update_mean(flux_mean%scalarGroundAbsorbedSolar, flux_temp%scalarGroundAbsorbedSolar, scale, nGRU)
  call update_mean(flux_mean%scalarTotalET, flux_temp%scalarTotalET, scale, nGRU)
  call update_mean(flux_mean%scalarNetRadiation, flux_temp%scalarNetRadiation, scale, nGRU)
  call update_mean_array(flux_mean%iLayerNrgFlux_m, flux_temp%iLayerNrgFlux_m, scale, nGRU,nLayers+1)
  call update_mean_array(flux_mean%iLayerAdvectiveFlux_m, flux_temp%iLayerAdvectiveFlux_m, scale, nGRU,nLayers+1)
  call update_mean_array(flux_mean%iLayerConductiveFlux_m, flux_temp%iLayerConductiveFlux_m, scale, nGRU,nLayers+1)
  call update_mean_array(flux_mean%iLayerLiqFluxSnow_m, flux_temp%iLayerLiqFluxSnow_m, scale, nGRU,nSnow+1)
  call update_mean_array(flux_mean%iLayerLiqFluxSoil_m, flux_temp%iLayerLiqFluxSoil_m, scale, nGRU,nSoil+1)
  call update_mean_array(flux_mean%mLayerNrgFlux_m, flux_temp%mLayerNrgFlux_m, scale, nGRU,nLayers)
  call update_mean(flux_mean%scalarSnowDrainage, flux_temp%scalarSnowDrainage, scale, nGRU)
  call update_mean(flux_mean%scalarRainPlusMelt, flux_temp%scalarRainPlusMelt, scale, nGRU)
  call update_mean(flux_mean%scalarSoilBaseflow, flux_temp%scalarSoilBaseflow, scale, nGRU)
  call update_mean(flux_mean%scalarTotalRunoff, flux_temp%scalarTotalRunoff, scale, nGRU)
  call update_mean(flux_mean%scalarSurfaceRunoff, flux_temp%scalarSurfaceRunoff, scale, nGRU)
  call update_mean(flux_mean%scalarSoilDrainage, flux_temp%scalarSoilDrainage, scale, nGRU)
  call update_mean_array(flux_mean%mLayerTranspire_m, flux_temp%mLayerTranspire_m, scale, nGRU,nSoil)
  call update_mean_array(flux_mean%mLayerHydCond_m, flux_temp%mLayerHydCond_m, scale, nGRU,nSoil)
  call update_mean_array(flux_mean%mLayerLiqFluxSoil_m, flux_temp%mLayerLiqFluxSoil_m, scale, nGRU,nSoil)
  call update_mean_array(flux_mean%mLayerBaseflow_m, flux_temp%mLayerBaseflow_m, scale, nGRU,nSoil)
  call update_mean_array(flux_mean%iLayerSatHydCond_m, flux_temp%iLayerSatHydCond_m, scale, nGRU,nSoil+1)
  call update_mean_array(flux_mean%mLayerSatHydCond_m, flux_temp%mLayerSatHydCond_m, scale, nGRU,nSoil)
  call update_mean_array(flux_mean%mLayerSatHydCondMP_m, flux_temp%mLayerSatHydCondMP_m, scale, nGRU,nSoil)
  call update_mean(flux_mean%scalarExfiltration, flux_temp%scalarExfiltration, scale, nGRU)
  call update_mean(flux_mean%scalarCanopySunlitPAR, flux_temp%scalarCanopySunlitPAR, scale, nGRU)
  call update_mean(flux_mean%scalarCanopyShadedPAR, flux_temp%scalarCanopyShadedPAR, scale, nGRU)
  call update_mean(flux_mean%scalarSoilResistance, flux_temp%scalarSoilResistance, scale, nGRU)
  call update_mean(flux_mean%scalarStomResistSunlit, flux_temp%scalarStomResistSunlit, scale, nGRU)
  call update_mean(flux_mean%scalarStomResistShaded, flux_temp%scalarStomResistShaded, scale, nGRU)
  call update_mean(flux_mean%scalarPhotosynthesisSunlit, flux_temp%scalarPhotosynthesisSunlit, scale, nGRU)
  call update_mean(flux_mean%scalarPhotosynthesisShaded, flux_temp%scalarPhotosynthesisShaded, scale, nGRU)
  call update_mean_array(flux_mean%mLayerColumnOutflow_m, flux_temp%mLayerColumnOutflow_m, scale, nGRU,nSoil)
  call update_mean_array(flux_mean%mLayerColumnInflow, flux_temp%mLayerColumnInflow, scale, nGRU,nSoil)
  if (nSnow /= 0) call update_mean_array(flux_mean%mLayerLiqFluxSnow_m, flux_temp%mLayerLiqFluxSnow_m,  scale, nGRU,nSnow)
end subroutine

subroutine update_mean_array(fmean, ftemp, scale, nGRU,nLayer)
  real(rkind),intent(inout),device :: fmean(:,:), ftemp(:,:)
  real(rkind),intent(in) :: scale
  integer(i4b),intent(in) :: nGRU,nLayer
  integer(i4b) :: iLayer,iGRU

  associate(mean => fmean, temp => ftemp)
  !$cuf kernel do(2) <<<*,*>>>
  do iGRU=1,nGRU
  do iLayer=1,nLayer
    mean(iLayer,iGRU) = mean(iLayer,iGRU) + temp(iLayer,iGRU)*scale
  enddo
end do
end associate
end subroutine
subroutine update_mean(mean, temp, scale,nGRU)
  real(rkind),device :: mean(:), temp(:)
  real(rkind) :: scale
  integer(i4b),intent(in) :: nGRU
  integer(i4b) :: iGRU
  associate(fmean => mean, ftemp => temp)
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
      fmean(iGRU) = fmean(iGRU) + ftemp(iGRU) * scale
    end do
    end associate
end subroutine



end module varSubstep_module
