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

module coupled_em_module

! homegrown solver data types
USE nrtype

! physical constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! data types
USE data_types,only:&
                    var_i,               & ! x%var(:)                (i4b)
                    var_d,               & ! x%var(:)                (rkind)
                    var_ilength,         & ! x%var(:)%dat            (i4b)
                    var_dlength,         & ! x%var(:)%dat            (rkind)
                    zLookup                ! x%z(:)%var(:)%lookup(:) (rkind)

! named variables for parent structures
USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
USE var_lookup,only:iLookPROG              ! named variables for structure elements
USE var_lookup,only:iLookDIAG              ! named variables for structure elements
USE var_lookup,only:iLookFLUX              ! named variables for structure elements
USE var_lookup,only:iLookPARAM             ! named variables for structure elements
USE var_lookup,only:iLookINDEX             ! named variables for structure elements
USE globalData,only:iname_snow             ! named variables for snow
USE globalData,only:iname_soil             ! named variables for soil

! named variables for child structures
! USE var_lookup,only:childFLUX_MEAN
USE var_lookup,only:iLookVarType           ! look up structure for variable typed

! metadata
USE globalData,only:flux_meta              ! metadata on the model fluxes
USE globalData,only:indx_meta              ! metadata on the model index variables
USE globalData,only:diag_meta              ! metadata on the model diagnostic variables
USE globalData,only:prog_meta              ! metadata on the model prognostic variables
! USE globalData,only:averageFlux_meta       ! metadata on the timestep-average model flux structure

! global data
USE globalData,only:data_step              ! time step of forcing data (s)
USE globalData,only:model_decisions        ! model decision structure
USE globalData,only:globalPrintFlag        ! the global print flag
USE globalData,only:realMissing            ! missing double precision number


! look-up values for the maximum interception capacity
USE mDecisions_module,only:         &
                      stickySnow,   &      ! maximum interception capacity an increasing function of temerature
                      lightSnow            ! maximum interception capacity an inverse function of new snow density

! look-up values for the groundwater parameterization
USE mDecisions_module,only:         &
                      qbaseTopmodel,&      ! TOPMODEL-ish baseflow parameterization
                      bigBucket    ,&      ! a big bucket (lumped aquifer model)
                      noExplicit           ! no explicit groundwater parameterization

! look-up values for the spatial representation of groundwater
USE mDecisions_module,only:         &
                      localColumn  ,&      ! separate groundwater representation in each local soil column
                      singleBasin          ! single groundwater store over the entire basin

! look-up values for the numerical method
USE mDecisions_module,only:         &
                      homegrown   ,&      ! homegrown backward Euler solution based on concepts from numerical recipes
                      kinsol       ,&      ! SUNDIALS backward Euler solution using Kinsol
                      ida                  ! SUNDIALS solution using IDA

! look-up values for the choice of variable in energy equations (BE residual or IDA state variable)
USE mDecisions_module,only:         &
                    closedForm,     &      ! use temperature with closed form heat capacity
                    enthalpyFormLU, &      ! use enthalpy with soil temperature-enthalpy lookup tables
                    enthalpyForm           ! use enthalpy with soil temperature-enthalpy analytical solution


! privacy
implicit none
private
public::coupled_em
! algorithmic parameters
real(rkind),parameter     :: verySmall=1.e-6_rkind   ! used as an additive constant to check if substantial difference among real numbers
contains


! ************************************************************************************************
! public subroutine coupled_em: run the coupled energy-mass model for one timestep
! ************************************************************************************************
subroutine coupled_em(&
                      ! model control
                      hruId,             & ! intent(in):    hruId
                      dt_init,           & ! intent(inout): used to initialize the size of the sub-step
                      dt_init_factor,    & ! intent(in):    Used to adjust the length of the timestep in the event of a failure
                      computeVegFlux,    & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                      fracJulDay,        & ! intent(in):    fractional julian days since the start of year
                      yearLength,        & ! intent(in):    number of days in the current year
                      nGRU, &
                      ! data structures (input)
                      type_data,         & ! intent(in):    local classification of soil veg etc. for each HRU
                      attr_data,         & ! intent(in):    local attributes for each HRU
                      forc_data,         & ! intent(in):    model forcing data
                      mpar_data,         & ! intent(in):    model parameters
                      bvar_data,         & ! intent(in):    basin-average variables
                      lookup_data,       & ! intent(in):    lookup tables
                      veg_param, decisions, &
                      ! data structures (input-output)
                      indx_data,         & ! intent(inout): model indices
                      prog_data,         & ! intent(inout): prognostic variables for a local HRU
                      diag_data,         & ! intent(inout): diagnostic variables for a local HRU
                      flux_data,         & ! intent(inout): model fluxes for a local HRU
                      ! error control
                      err,message)         ! intent(out):   error control
  ! structure allocations
  USE allocspace_module,only:allocLocal             ! allocate local data structures
  USE allocspace_module,only:resizeData             ! clone a data structure
  ! simulation of fluxes and residuals given a trial state vector
  USE soil_utils_module,only:liquidHead             ! compute the liquid water matric potential
  ! preliminary subroutines
  USE vegPhenlgy_module,only:vegPhenlgy_d             ! compute vegetation phenology
  USE vegNrgFlux_module,only:wettedFrac             ! compute wetted fraction of the canopy (used in sw radiation fluxes)
  USE snowAlbedo_module,only:snowAlbedo             ! compute snow albedo
  USE vegSWavRad_module,only:vegSWavRad             ! compute canopy sw radiation fluxes
  USE canopySnow_module,only:canopySnow             ! compute interception and unloading of snow from the vegetation canopy
  USE volicePack_module,only:newsnwfall             ! compute change in the top snow layer due to throughfall and unloading
  USE volicePack_module,only:volicePack             ! merge and sub-divide snow layers, if necessary
  USE diagn_evar_module,only:diagn_evar             ! compute diagnostic energy variables -- thermal conductivity and heat capacity
  ! the model solver
  USE indexState_module,only:indexState             ! define indices for all model state variables and layers
  USE opSplittin_module,only:opSplittin             ! solve the system of thermodynamic and hydrology equations for a given substep
  USE time_utils_module,only:elapsedSec             ! calculate the elapsed time
  ! additional subroutines
  USE tempAdjust_module,only:tempAdjust             ! adjust snow temperature associated with new snowfall
  USE var_derive_module,only:calcHeight,calcHeight_d             ! module to calculate height at layer interfaces and layer mid-point
  USE computSnowDepth_module,only:computSnowDepth_d   ! compute snow depth
  USE enthalpyTemp_module,only:T2enthTemp_veg,T2enthTemp_veg_d       ! convert temperature to enthalpy for vegetation
  USE enthalpyTemp_module,only:T2enthTemp_snow,T2enthTemp_snow_d      ! convert temperature to enthalpy for snow
  USE enthalpyTemp_module,only:T2enthTemp_soil,T2enthTemp_soil2      ! convert temperature to enthalpy for soil
  USE enthalpyTemp_module,only:enthTemp_or_enthalpy ! add phase change terms to delta temperature component of enthalpy or vice versa
  use device_data_types
  use initialize_device

  implicit none

  integer(i8b),intent(in)              :: hruId                  ! hruId
  real(rkind),intent(inout)            :: dt_init                ! used to initialize the size of the sub-step
  integer(i4b),intent(in)              :: dt_init_factor         ! Used to adjust the length of the timestep in the event of a failure
  logical(lgt),intent(inout)           :: computeVegFlux         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  ! data structures (input)
  type(type_data_device),intent(in)               :: type_data              ! type of vegetation and soil
  type(attr_data_device),intent(in)               :: attr_data              ! spatial attributes
  type(forc_data_device),intent(in)               :: forc_data              ! model forcing data
  type(mpar_data_device),intent(in)         :: mpar_data              ! model parameters
  type(bvar_data_device),intent(in)         :: bvar_data              ! basin-average model variables
  type(zLookup),intent(in)             :: lookup_data            ! lookup tables
  ! data structures (input-output)
  ! type(var_ilength),intent(inout)      :: indx_data2              ! state vector geometry
  type(prog_data_device),intent(inout)      :: prog_data              ! prognostic variables for a local HRU
  type(diag_data_device),intent(inout)      :: diag_data              ! diagnostic variables for a local HRU
  type(flux_data_device),intent(inout)      :: flux_data              ! model fluxes for a local HRU
  real(rkind),intent(in)               :: fracJulDay             ! fractional julian days since the start of year
  integer(i4b),intent(in)              :: yearLength             ! number of days in the current year
  ! error control
  integer(i4b),intent(out)             :: err                    ! error code
  character(*),intent(out)             :: message                ! error message
  ! =====================================================================================================================================================
  ! =====================================================================================================================================================
  ! local variables
  character(len=256)                   :: cmessage               ! error message
  integer(i4b)                         :: nSnow                  ! number of snow layers
  integer(i4b)                         :: nSoil                  ! number of soil layers
  integer(i4b)                         :: nLayers                ! total number of layers
  integer(i4b)                         :: nState                 ! total number of state variables
  real(rkind)                          :: dtSave                 ! length of last input model whole sub-step (seconds)
  real(rkind)                          :: dt_sub                 ! length of model sub-step (seconds)
  real(rkind)                          :: dt_wght                ! weight applied to model sub-step (dt_sub/data_step)
  real(rkind)                          :: lyr_wght               ! weight applied to domain layer (layer_depth/domain_depth)
  real(rkind)                          :: dt_solv                ! seconds in the data step that have been completed
  real(rkind)                          :: dtMultiplier           ! time step multiplier (-) based on what happenned in "opSplittin"
  real(rkind)                          :: minstep,maxstep        ! minimum and maximum time step length (seconds)
  real(rkind)                          :: maxstep_op             ! maximum time step length (seconds) to run opSplittin over
  real(rkind)                          :: whole_step             ! step the surface pond drainage and sublimation calculated over
  integer(i4b)                         :: nsub                   ! number of substeps
  integer(i4b)                         :: nsub_success           ! number of successful substeps
  logical(lgt)                         :: computeVegFluxOld      ! flag to indicate if we are computing fluxes over vegetation on the previous sub step
  logical(lgt)                         :: includeAquifer         ! flag to denote that an aquifer is included
  logical(lgt)                         :: modifiedLayers         ! flag to denote that snow layers were modified
  logical(lgt)                         :: modifiedVegState       ! flag to denote that vegetation states were modified
  integer(i4b),device                         :: nLayersRoots(nGRU)           ! number of soil layers that contain roots
  real(rkind),device                          :: exposedVAI(nGRU)             ! exposed vegetation area index
  real(rkind),device                          :: dCanopyWetFraction_dWat(nGRU) ! derivative in wetted fraction w.r.t. canopy total water (kg-1 m2)
  real(rkind),device                          :: dCanopyWetFraction_dT(nGRU)   ! derivative in wetted fraction w.r.t. canopy temperature (K-1)
  real(rkind),device                :: varNotUsed1 ! variables used to calculate derivatives (not needed here)
  real(rkind),device                :: varNotUsed2 ! variables used to calculate derivatives (not needed here)
  integer(i4b)                         :: iLayer                 ! index of model layers
  real(rkind)                          :: superflousSub          ! superflous sublimation (kg m-2 s-1)
  real(rkind)                          :: superflousNrg          ! superflous energy that cannot be used for sublimation (W m-2 [J m-2 s-1])
  integer(i4b)                         :: ixSolution             ! solution method used by opSplittin
  logical(lgt)                         :: firstSubStep           ! flag to denote if the first time step
  logical(lgt)                         :: stepFailure            ! flag to denote the need to reduce length of the coupled step and try again
  logical(lgt)                         :: tooMuchMelt            ! flag to denote that there was too much melt in a given time step
  logical(lgt)                         :: tooMuchSublim          ! flag to denote that there was too much sublimation in a given time step
  logical(lgt)                         :: doLayerMerge           ! flag to denote the need to merge snow layers
  logical(lgt)                         :: pauseFlag              ! flag to pause execution
  logical(lgt),parameter               :: backwardsCompatibility=.false.  ! flag to denote a desire to ensure backwards compatibility with previous branches for end of time step flux only 
  logical(lgt)                         :: checkMassBalance_ds    ! flag to check the mass balance over the data step
  type(indx_data_device)                    :: indx_temp              ! temporary model index variables saved only on outer loop
  type(indx_data_device)                    :: indx_temp0             ! temporary model index variables saved every time
  type(prog_data_device)                    :: prog_temp              ! temporary model prognostic variables
  type(diag_data_device)                    :: diag_temp              ! temporary model diagnostic variables
  real(rkind),device,allocatable              :: mLayerVolFracIceInit(:,:)! initial vector for volumetric fraction of ice (-)
  ! check SWE
  real(rkind)                          :: oldSWE                 ! SWE at the start of the substep
  real(rkind),device                          :: innerEffRainfall(nGRU)       ! inner step average effective rainfall into snow (kg m-2 s-1)
  real(rkind),device                          :: effRainfall(nGRU)            ! timestep-average effective rainfall into snow (kg m-2 s-1)
  real(rkind),device                          :: sfcMeltPond(nGRU)            ! surface melt pond (kg m-2)
  ! energy fluxes
  integer(i4b)                         :: iSoil                  ! index of soil layers
  type(flux_data_device)                    :: flux_mean              ! timestep-average model fluxes for a local HRU
  type(flux_data_device)                    :: flux_inner             ! inner step average model fluxes for a local HRU
  real(rkind),allocatable,device              :: meanSoilCompress(:,:)    ! timestep-average soil compression by layer
  real(rkind),allocatable,device              :: innerSoilCompress(:,:)   ! inner step average soil compression by layer
  ! sublimation sums over substep and means over data_step
  real(rkind),device                          :: sumCanopySublimation(nGRU)   ! sum of sublimation from the vegetation canopy (kg m-2 s-1) over substep
  real(rkind),device                          :: sumSnowSublimation(nGRU)     ! sum of sublimation from the snow surface (kg m-2 s-1) over substep
  real(rkind),device                          :: sumLatHeatCanopyEvap(nGRU)   ! sum of latent heat flux for evaporation from the canopy to the canopy air space (W m-2) over substep
  real(rkind),device                          :: sumSenHeatCanopy(nGRU)       ! sum of sensible heat flux from the canopy to the canopy air space (W m-2) over substep
  real(rkind),device                          :: meanCanopySublimation(nGRU)  ! timestep-average sublimation from the vegetation canopy (kg m-2 s-1)
  real(rkind),device                          :: meanLatHeatCanopyEvap(nGRU)  ! timestep-average latent heat flux for evaporation from the canopy to the canopy air space (W m-2)
  real(rkind),device                          :: meanSenHeatCanopy(nGRU)      ! timestep-average sensible heat flux from the canopy to the canopy air space (W m-2)
  ! balance checks
  logical(lgt)                         :: bal_veg                ! flag to denote if computed a vegetation balance
  logical(lgt)                         :: bal_snow               ! flag to denote if computed a snow balance
  logical(lgt)                         :: bal_soil               ! flag to denote if computed a soil balance
  logical(lgt)                         :: bal_aq                 ! flag to denote if computed an aquifer balance
  integer(i4b)                         :: iVar                   ! loop through model variables
  real(rkind)                          :: scalarCanopyWatBalError! water balance error for the vegetation canopy (kg m-2)
  real(rkind)                          :: scalarInitCanopyLiq    ! initial liquid water on the vegetation canopy (kg m-2)
  real(rkind)                          :: scalarInitCanopyIce    ! initial ice          on the vegetation canopy (kg m-2)
  real(rkind)                          :: balanceCanopyWater0    ! total water stored in the vegetation canopy at the start of the step (kg m-2)
  real(rkind)                          :: balanceSoilWater0      ! total soil storage at the start of the step (kg m-2)
  real(rkind)                          :: balanceAquifer0        ! total aquifer storage at the start of the step (kg m-2)
  real(rkind),device                          :: innerBalance(4,nGRU)        ! inner step balances for domain with one layer
  real(rkind),device                          :: meanBalance(8,nGRU)         ! timestep-average balances for domains
  real(rkind),allocatable,device              :: innerBalanceLayerMass(:,:) ! inner step balances for domain with multiple layers
  real(rkind),allocatable,device              :: innerBalanceLayerNrg(:,:)  ! inner step balances for domain with multiple layers
  ! test balance checks
  logical(lgt),parameter               :: printBalance=.false.   ! flag to print the balance checks
  real(rkind),allocatable              :: liqSnowInit(:)         ! volumetric liquid water conetnt of snow at the start of the time step
  real(rkind),allocatable              :: liqSoilInit(:)         ! soil moisture at the start of the time step
  ! timing information
  integer(kind=8)                      :: count_rate 
  integer(kind=8)                      :: i_start, i_end
  real                                 :: elapsed_time
  real(rkind)                          :: mean_step_dt_sub       ! mean solution step for the sub-step
  real(rkind)                          :: sumStepSize            ! sum solution step for the data step
  ! outer loop control
  integer(i4b)                         :: be_steps               ! number of substeps for a BE solver
  logical(lgt)                         :: firstInnerStep         ! flag to denote if the first time step in maxstep subStep
  logical(lgt)                         :: lastInnerStep          ! flag to denote if the last time step in maxstep subStep
  logical(lgt)                         :: do_outer               ! flag to denote if doing the outer steps surrounding the call to opSplittin
  real(rkind)                          :: dt_solvInner           ! seconds in the maxstep subStep that have been completed
  logical(lgt),parameter               :: computNrgBalance_var=.true. ! flag to compute enthalpy, must have computNrgBalance true in varSubStep (will compute enthalpy for BE even if not using enthalpy formulation)
  logical(lgt)                         :: computeEnthalpy        ! flag to compute enthalpy regardless of the model decision
  logical(lgt)                         :: enthalpyStateVec       ! flag if enthalpy is a state variable (IDA)
  logical(lgt)                         :: use_lookup             ! flag to use the lookup table for soil enthalpy, otherwise use analytical solution
  ! ----------------------------------------------------------------------------------------------------------------------------------------------
    type(decisions_device) :: decisions
    type(veg_parameters) :: veg_param
    integer(i4b) :: nGRU
    type(dim3) :: threads,blocks
    integer(i4b) :: iGRU
    type(indx_data_device) :: indx_data
! real(rkind),device :: sumCanopySublimation_d(nGRU)
! real(rkind),device :: sumSnowSublimation_d(nGRU)
! real(rkind),device :: sumLatHeatCanopyEvap_d(nGRU)
! real(rkind),device :: sumSenHeatCanopy_d(nGRU)
! real(rkind),device :: innerEffRainfall_d(nGRU)
real(rkind),device :: whole_step_d
! real(rkind),device :: sfcMeltPond_d(nGRU)
real(rkind),device :: layerDepthSnow(nGRU),layerDepthSoil(nGRU)
! real(rkind),device :: meanBalance_d(8,nGRU)
! real(rkind),device :: meanCanopySublimation_d(nGRU)
! real(rkind),device :: meanLatHeatCanopyEvap_d(nGRU)
! real(rkind),device :: meanSenHeatCanopy_d(nGRU)
! real(rkind),device :: effRainfall_d(nGRU)
real(rkind),device :: sum_depth(nGRU)
logical(lgt),device :: computeVegFlux_d(nGRU)
real(rkind),device :: fracJulDay_d
integer(i4b),device :: yearLength_d

fracJulDay_d = fracJulDay
yearLength_d = yearLength

call allocate_device_decisions(decisions)
! call allocate_device_indx_data(indx_data,indx_data2,nGRU)
varNotUsed1=-9999._rkind
varNotUsed2=-9999._rkind
  ! initialize error control
  err=0; message="coupled_em/"

  call allocate_device_flux_prev(flux_mean,indx_data%nSoil,nGRU)
        call allocate_device_flux_prev(flux_inner,indx_data%nSoil,nGRU)
whole_step_d = whole_step

  ! This is the start of a data step for a local HRU

  ! get the start time
 ! get the start time
  CALL system_clock(count_rate=count_rate)
  CALL system_clock(i_start)

  ! check that the decision is supported
  if(model_decisions(iLookDECISIONS%groundwatr)%iDecision==bigBucket .and. &
      model_decisions(iLookDECISIONS%spatial_gw)%iDecision/=localColumn)then
    message=trim(message)//'expect "spatial_gw" decision to equal localColumn when "groundwatr" decision is bigBucket'
    err=20; return
  endif

  ! check if the aquifer is included
  includeAquifer = (model_decisions(iLookDECISIONS%groundwatr)%iDecision==bigBucket)

  ! initialize variables
  call initialize_coupled_em

  ! link canopy depth to the information in the data structure
  canopy: associate(&
    snowfrz_scale     => mpar_data%snowfrz_scale    ,& ! scaling parameter for the snow freezing curve (K-1)
    canopyDepth       => diag_data%scalarCanopyDepth ,& ! depth of the vegetation canopy (m) 
    specificHeatVeg   => mpar_data%specificHeatVeg  ,& ! specific heat of vegetation (J kg-1 K-1)
    maxMassVegetation => mpar_data%maxMassVegetation & ! maximum mass of vegetation (kg m-2)
    )

    ! start by NOT pausing
    pauseFlag=.false.

    ! start by assuming that the step is successful
    stepFailure  = .false.
    doLayerMerge = .false.

    ! initialize flags to modify the veg layers or modify snow layers
    modifiedLayers    = .false.    ! flag to denote that snow layers were modified
    modifiedVegState  = .false.    ! flag to denote that vegetation states were modified

    ! define the first step and first and last inner steps
    firstSubStep = .true.
    firstInnerStep = .true.
    lastInnerStep = .false.

    ! create temporary data structures for prognostic variables
    ! call resizeData(prog_meta(:),prog_data,prog_temp,err=err,message=cmessage)
    prog_temp = prog_data
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! create temporary data structures for diagnostic variables
    ! call resizeData(diag_meta(:),diag_data,diag_temp,err=err,message=cmessage)
  call allocate_device_diag_temp(diag_temp,nGRU,indx_data%nSoil)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! create temporary data structures for index variables
    indx_temp = indx_data
    indx_temp0 = indx_data
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! allocate space for the local fluxes

    ! initialize surface melt pond
    sfcMeltPond       = 0._rkind  ! change in storage associated with the surface melt pond (kg m-2)

    ! initialize fluxes to average over data_step (averaged over substep in varSubStep)
    call zero_device_flux(flux_mean)

    ! associate local variables with information in the data structures
    associate(&
      nSnow => indx_data%nSnow, &
    ! model decisions
    ixNumericalMethod    => model_decisions(iLookDECISIONS%num_method)%iDecision            ,& ! choice of numerical solver
    ixNrgConserv         => model_decisions(iLookDECISIONS%nrgConserv)%iDecision            ,& ! choice of variable in either energy backward Euler residual or IDA state variable
    ! state variables in the vegetation canopy
    scalarCanopyLiq      => prog_data%scalarCanopyLiq                 ,& ! canopy liquid water (kg m-2)
    scalarCanopyIce      => prog_data%scalarCanopyIce                 ,& ! canopy ice content (kg m-2)
    ! state variables in the soil domain
    mLayerDepth          => prog_data%mLayerDepth       ,& ! depth of each soil layer (m)
    mLayerVolFracIce     => prog_data%mLayerVolFracIce  ,& ! volumetric ice content in each soil layer (-)
    mLayerVolFracLiq     => prog_data%mLayerVolFracLiq  ,& ! volumetric liquid water content in each soil layer (-)
    scalarAquiferStorage => prog_data%scalarAquiferStorage            ,& ! aquifer storage (m)
    scalarTotalSoilIce   => diag_data%scalarTotalSoilIce              ,& ! total ice in the soil column (kg m-2)
    scalarTotalSoilLiq   => diag_data%scalarTotalSoilLiq              ,& ! total liquid water in the soil column (kg m-2)
    scalarTotalSoilWat   => diag_data%scalarTotalSoilWat               & ! total water in the soil column (kg m-2)
    ) ! (association of local variables with information in the data structures


      ! identify the need to check the mass balance, both methods should work if tolerance coarse enough
      checkMassBalance_ds = .false. ! IDA balance agreement levels are controlled by set tolerances
      
      ! set the number of substeps for a BE solver
      be_steps = 1_i4b ! IDA does not use substeps

      ! set the flag to compute enthalpy, may want to have this true always if want to output enthalpy
      computeEnthalpy  = .false.
      enthalpyStateVec = .false.
      use_lookup       = .false.
       if(ixNrgConserv.ne.closedForm) enthalpyStateVec = .true. ! enthalpy as state variable
      if(ixNrgConserv==enthalpyFormLU) use_lookup = .true. ! use lookup tables for soil temperature-enthalpy instead of analytical solution


      ! compute total soil moisture and ice at the *START* of the step (kg m-2)
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        scalarTotalSoilLiq(iGRU) = 0
        scalarTotalSoilIce(iGRU) = 0
        do iLayer=nSnow(iGRU)+1,nLayers
      scalarTotalSoilLiq(iGRU) = scalarTotalSoilLiq(iGRU) + iden_water*mLayerVolFracLiq(iLayer,iGRU)*mLayerDepth(iLayer,iGRU)
      scalarTotalSoilIce(iGRU) = scalarTotalSoilIce(iGRU) + iden_water*mLayerVolFracIce(iLayer,iGRU)*mLayerDepth(iLayer,iGRU)  ! NOTE: no expansion and hence use iden_water
        end do
      end do

    ! end association of local variables with information in the data structures
    end associate
    ! short-cut to the algorithmic control parameters
    ! NOTE - temporary assignment of minstep to foce something reasonable
    ! changing the maxstep parameter will make the outer and inner loop computations here in coupled_em happen more frequently
    ! changing the be_steps parameter will make the inner loop computations in opSplittin happen more frequently (e.g. be_steps = 32.0 give BE32)
    minstep = mpar_data%minstep  ! minimum time step (s)
    maxstep = mpar_data%maxstep  ! maximum time step (s)
    maxstep_op = mpar_data%maxstep/be_steps  ! maximum time step (s) to run opSplittin over
    associate(iLayerHeight => prog_data%iLayerHeight, nSnow => indx_data%nSnow, rootingDepth => mpar_data%rootingDepth)
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        nLayersRoots(iGRU) = 0
        do iLayer=nSnow(iGRU),nLayers-1
           if (iLayerHeight(iLayer,iGRU) < rootingDepth-verySmall) nLayersRoots(iGRU) = nLayersRoots(iGRU) + 1
        end do
      end do
      end associate
    ! if(nLayersRoots == 0)then
    !   message=trim(message)//'no roots within the soil profile'
    !   err=20; return
    ! end if

    ! ! compute the number of layers with roots
    ! associate(nLayersRoots => nLayersRoots_d, iLayerHeight => prog_data_d%iLayerHeight, nSnow => indx_data_d%nSNow, rootingDepth => mpar_data%rootingDepth)
    !   !$cuf kernel do(1) <<<*,*>>>
    !   do iGRU=1,nGRU
    !     nLayersRoots(iGRU) = 0
    !     do iLayer=nSnow(iGRU),nLayers-1
    !       if (iLayerHeight(iLayer,iGRU) < rootingDepth-verySmall) nLayersRoots(iGRU) = nLayersRoots(iGRU) + 1
    !     end do
    !   end do
    !   end associate

    ! define the foliage nitrogen factor
    diag_data%scalarFoliageNitrogenFactor = 1._rkind  ! foliage nitrogen concentration (1.0 = saturated)

    ! save SWE
    ! oldSWE = prog_data%var(iLookPROG%scalarSWE)%dat(1)

    ! *** compute phenology...
    ! ------------------------

    ! compute the temperature of the root zone: used in vegetation phenology
    associate(scalarRootZoneTemp => diag_data%scalarRootZoneTemp, mLayerTemp => prog_data%mLayerTemp, nSnow => indx_data%nSnow, nLayersRoots => nLayersRoots)
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        scalarRootZoneTemp(iGRU) = 0
        do iLayer=nSnow(iGRU) + 1, nSnow(iGRU)+nLayersRoots(iGRU)
          scalarRootZoneTemp(iGRU) = scalarRootZoneTemp(iGRU) + mLayerTemp(iLayer,iGRU) / real(nLayersRoots(iGRU), kind(rkind))
        end do
      end do
    end associate

    ! remember if we compute the vegetation flux on the previous sub-step
    computeVegFluxOld = computeVegFlux
computeVegFlux_d = computeVegFlux

    ! compute the exposed LAI and SAI and whether veg is buried by snow
    call vegPhenlgy_d(&
                    ! model control
                    model_decisions,             & ! intent(in):    model decisions
                    nGRU, &
                    ! input/output: data structures
                    fracJulDay_d,                  & ! intent(in):    fractional julian days since the start of year
                    yearLength_d,                  & ! intent(in):    number of days in the current year
                    type_data,                   & ! intent(in):    type of vegetation and soil
                    attr_data,                   & ! intent(in):    spatial attributes
                    mpar_data,                   & ! intent(in):    model parameters
                    prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                    diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                    decisions, veg_param, &
                    ! output
                    computeVegFlux_d,              & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                    diag_data%scalarcanopyDepth,                 & ! intent(out): canopy depth (m)
                    exposedVAI,                  & ! intent(out): exposed vegetation area index (m2 m-2)
                    err,cmessage)                  ! intent(out): error control
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

    ! check
    ! if(computeVegFlux)then
    !   if(canopyDepth < epsilon(canopyDepth))then
    !     message=trim(message)//'canopy depth is zero when computeVegFlux flag is .true.'
    !     err=20; return
    !   endif
    ! endif

    associate(snowIncept => decisions%snowIncept, &
      scalarCanopyIceMax => diag_data%scalarCanopyIceMax, &
      refInterceptCapSnow => mpar_data%refInterceptCapSnow, &
      scalarCanopyLiqMax => diag_data%scalarCanopyLiqMax, &
      refInterceptCapRain => mpar_data%refInterceptCapRain)

      computeVegFlux = computeVegFlux_d(1)
    ! flag the case where number of vegetation states has changed
    modifiedVegState = (computeVegFlux.neqv.computeVegFluxOld)
    ! *** compute wetted canopy area...
    ! ---------------------------------
    
    ! compute maximum canopy liquid water (kg m-2)
!$cuf kernel do(1) <<<*,*>>>
do iGRU=1,nGRU
    scalarCanopyLiqMax(iGRU) = refInterceptCapRain*exposedVAI(iGRU)
! end do

!       !$cuf kernel do(1) <<<*,*>>>
!       do iGRU=1,nGRU

    ! compute maximum canopy ice content (kg m-2)
    ! NOTE 1: this is used to compute the snow fraction on the canopy, as used in *BOTH* the radiation AND canopy sublimation routines
    ! NOTE 2: this is a different variable than the max ice used in the throughfall (snow interception) calculations
    ! NOTE 3: use maximum per unit leaf area storage capacity for snow (kg m-2)
    select case(snowIncept)
      case(lightSnow);  scalarCanopyIceMax(iGRU) = exposedVAI(iGRU)*refInterceptCapSnow
      case(stickySnow); scalarCanopyIceMax(iGRU) = exposedVAI(iGRU)*refInterceptCapSnow*4._rkind
      ! case default; message=trim(message)//'unable to identify option for maximum branch interception capacity'; err=20; return
    end select ! identifying option for maximum branch interception capacity
  end do
end associate

        ! compute wetted fraction of the canopy
    ! NOTE: assume that the wetted fraction is constant over the substep for the radiation calculations
    if(computeVegFlux)then

      associate(scalarCanopyTemp => prog_data%scalarCanopyTemp, &
        scalarCanopyLiq => prog_data%scalarCanopyLiq, &
        scalarCanopyIce => prog_data%scalarCanopyIce, &
        scalarCanopyLiqMax => diag_data%scalarCanopyLiqMax, &
        scalarCanopyIceMax => diag_data%scalarCanopyIceMax, &
        scalarCanopyWetFraction => diag_data%scalarCanopyWetFraction, &
        false => decisions%false,&
        canopyWettingFactor => mpar_data%canopyWettingFactor, &
        canopyWettingExp => mpar_data%canopyWettingExp)
      ! compute wetted fraction of the canopy
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
      call wettedFrac(&
                      ! input
                      false,                                                      & ! flag to denote if derivatives are required
                      (scalarCanopyTemp(iGRU) < Tfreeze), & ! flag to denote if the canopy is frozen
                      varNotUsed1,                                                  & ! derivative in canopy liquid w.r.t. canopy temperature (kg m-2 K-1)
                      varNotUsed2,                                                  & ! fraction of liquid water on the canopy
                      scalarCanopyLiq(iGRU),              & ! canopy liquid water (kg m-2)
                      scalarCanopyIce(iGRU),              & ! canopy ice (kg m-2)
                      scalarCanopyLiqMax(iGRU),           & ! maximum canopy liquid water (kg m-2)
                      scalarCanopyIceMax(iGRU),           & ! maximum canopy ice content (kg m-2)
                      canopyWettingFactor,         & ! maximum wetted fraction of the canopy (-)
                      canopyWettingExp,            & ! exponent in canopy wetting function (-)
                      ! output
                      scalarCanopyWetFraction(iGRU),      & ! canopy wetted fraction (-)
                      dCanopyWetFraction_dWat(iGRU),                                      & ! derivative in wetted fraction w.r.t. canopy liquid water content (kg-1 m2)
                      dCanopyWetFraction_dT(iGRU)                                        & ! derivative in wetted fraction w.r.t. canopy liquid water content (kg-1 m2)
                      )
        end do
      ! if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
end associate
    ! vegetation is completely buried by snow (or no veg exists at all)
    else
      diag_data%scalarCanopyWetFraction = 0._rkind
      dCanopyWetFraction_dWat                                 = 0._rkind
      dCanopyWetFraction_dT                                   = 0._rkind
    end if

    ! *** compute snow albedo...
    ! --------------------------
    ! NOTE: this should be done before the radiation calculations
    ! NOTE: uses snowfall; should really use canopy throughfall + canopy unloading
    call snowAlbedo(&
                    ! input: model control
                    data_step,                   & ! intent(in): model time step (s)
                    indx_data%nSnow,                 & ! intent(in): logical flag to denote if snow is present
                    nGRU, &
                    ! input/output: data structures
                    model_decisions,             & ! intent(in):    model decisions
                    decisions, &
                    mpar_data,                   & ! intent(in):    model parameters
                    flux_data,                   & ! intent(in):    model flux variables
                    diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                    prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                    ! output: error control
                    err,cmessage)                  ! intent(out): error control
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if


    ! *** compute canopy sw radiation fluxes...
    ! -----------------------------------------
    call vegSWavRad(&
                    data_step,                    & ! intent(in):    time step (s) -- only used in Noah-MP radiation, to compute albedo
                    indx_data%nSnow,                        & ! intent(in):    number of snow layers
                    nSoil,                        & ! intent(in):    number of soil layers
                    nLayers,                      & ! intent(in):    total number of layers
                    nGRU, &
                    computeVegFlux,               & ! intent(in):    logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                    decisions, veg_param, &
                    type_data,                    & ! intent(in):    type of vegetation and soil
                    prog_data,                    & ! intent(inout): model prognostic variables for a local HRU
                    diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                    flux_data,                    & ! intent(inout): model flux variables
                    err,cmessage)                   ! intent(out):   error control
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

    ! *** compute canopy throughfall and unloading...
    ! -----------------------------------------------
    ! NOTE 1: this needs to be done before solving the energy and liquid water equations, to account for the heat advected with precipitation (and throughfall/unloading)
    ! NOTE 2: the unloading flux is computed using canopy drip (scalarCanopyLiqDrainage) from the previous time step
    ! this changes canopy ice
    call canopySnow(&
    nGRU, &
                    ! input: model control
                    decisions%data_step,                   & ! intent(in): time step (seconds)
                    exposedVAI,                  & ! intent(in): exposed vegetation area index (m2 m-2)
                    computeVegFlux_d,              & ! intent(in): flag to denote if computing energy flux over vegetation
                    ! input/output: data structures
                    decisions,             & ! intent(in):    model decisions
                    forc_data,                   & ! intent(in):    model forcing data
                    mpar_data,                   & ! intent(in):    model parameters
                    diag_data,                   & ! intent(in):    model diagnostic variables for a local HRU
                    prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                    flux_data,                   & ! intent(inout): model flux variables
                    ! output: error control
                    err,cmessage)                  ! intent(out): error control
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

    ! adjust canopy temperature to account for new snow
    if(computeVegFlux)then ! logical flag to compute vegetation fluxes (.false. if veg buried by snow)
      call tempAdjust(&
      nGRU, &
                      ! input: derived parameters
                      diag_data%scalarcanopyDepth,                 & ! intent(in):    canopy depth (m)
                      ! input/output: data structures
                      mpar_data,                   & ! intent(in):    model parameters
                      prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                      diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                      ! output: error control
                      err,cmessage)                  ! intent(out): error control
      if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

      ! change enthalpy based on new canopy temperature and water content, only if will need enthalpy for energy balance
      if(enthalpyStateVec .or. computeEnthalpy)then
        ! associate local variables with variables in the data structures
        enthalpyVeg: associate(&
        canopyDepth => diag_data%scalarCanopyDepth, &
        specificHeatVeg => mpar_data%specificHeatVeg, &
        maxMassVegetation => mpar_data%maxMassVegetation, &
        snowfrz_scale => mpar_data%snowfrz_scale, &
          ! state variables in the vegetation canopy
          scalarCanopyTemp     => prog_data%scalarCanopyTemp     ,& ! canopy temperature (K)
          scalarCanopyEnthTemp => diag_data%scalarCanopyEnthTemp ,& ! canopy temperature component of enthalpy (J m-3)
          scalarCanopyEnthalpy => prog_data%scalarCanopyEnthalpy ,& ! enthalpy of the vegetation canopy (J m-3)
          scalarCanopyLiq      => prog_data%scalarCanopyLiq      ,& ! mass of liquid water on the vegetation canopy (kg m-2)
          scalarCanopyIce      => prog_data%scalarCanopyIce       & ! mass of ice on the vegetation canopy (kg m-2)
          )  ! (associate local variables with model parameters)       
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
          call T2enthTemp_veg_d(&
                          canopyDepth(iGRU),            & ! intent(in): canopy depth (m)
                          specificHeatVeg,        & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                          maxMassVegetation,      & ! intent(in): maximum mass of vegetation (kg m-2)
                          snowfrz_scale,          & ! intent(in): scaling parameter for the snow freezing curve  (K-1)
                          scalarCanopyTemp(iGRU),       & ! intent(in): canopy temperature (K)
                          (scalarCanopyLiq(iGRU)+scalarCanopyIce(iGRU)), & ! intent(in): canopy water content (kg m-2)
                          scalarCanopyEnthTemp(iGRU))     ! intent(out): temperature component of enthalpy of the vegetation canopy (J m-3)
          scalarCanopyEnthalpy(iGRU) = scalarCanopyEnthTemp(iGRU)  - LH_fus * scalarCanopyIce(iGRU)/ canopyDepth(iGRU) ! new ice and/or temperature
          end do
        end associate enthalpyVeg
      end if ! (need to recalculate enthalpy state variable)
    end if ! if computing fluxes over vegetation

    ! initialize drainage and throughfall
    ! NOTE 1: this needs to be done before solving the energy and liquid water equations, to account for the heat advected with precipitation
    ! NOTE 2: this initialization needs to be done AFTER the call to canopySnow, since canopySnow uses canopy drip drom the previous time step
    if(.not.computeVegFlux)then
      flux_data%scalarThroughfallRain = flux_data%scalarRainfall
    else
      flux_data%scalarThroughfallRain = 0._rkind
    end if
    flux_data%scalarCanopyLiqDrainage = 0._rkind

    ! ****************************************************************************************************
    ! *** MAIN SOLVER ************************************************************************************
    ! ****************************************************************************************************

    ! initialize the length of the sub-step and counters
    whole_step = maxstep
    dt_solv       = 0._rkind   ! length of time step that has been completed (s)
    dt_solvInner  = 0._rkind   ! length of time step that has been completed (s) in whole_step subStep
    dt_init = min(data_step,whole_step,maxstep_op) / dt_init_factor  ! initial substep length (s)
    dt_sub = dt_init
    dtSave  = whole_step       ! length of whole substep

    ! initialize the number of sub-steps
    nsub = 0
    nsub_success = 0

    ! initialize if used a balance
    bal_veg = .false.
    bal_snow = .false.
    bal_soil = .false.
    bal_aq   = .false.


    ! loop through sub-steps
    substeps: do  ! continuous do statement with exit clause (alternative to "while")

      dt_sub = min(data_step,whole_step,maxstep_op,dt_sub) ! adjust for possible whole_step changes

      ! print progress
      if(globalPrintFlag)then
        write(*,'(a,1x,4(f13.5,1x))') ' start of step: dt_init, dt_sub, dt_solv, data_step: ', dt_init, dt_sub, dt_solv, data_step
        print*, 'stepFailure = ', stepFailure
        print*, 'before resizeData: nSnow, nSoil = ', nSnow, nSoil
      endif

      ! increment the number of sub-steps
      nsub = nsub+1

      ! resize the "indx_data" structure
      ! NOTE: this is necessary because the length of index variables depends on a given split
      !        --> the resize here is overwritten later (in indexSplit)
      !        --> admittedly ugly, and retained for now
      if(stepFailure)then ! resize temp to current data, later in code current data is set to lastInnerStep data
        ! call resizeData(indx_meta(:),indx_temp,indx_data,err=err,message=cmessage)
        if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
      else ! resize current data to temp0, temp0 is saved for next run
        ! call resizeData(indx_meta(:),indx_data,indx_temp0,err=err,message=cmessage)
        ! if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
        ! do iVar=1,size(indx_data%var)
          indx_temp0 = indx_data
        ! end do
      endif
      ! check if on outer loop, always do outer if after failed step and on then on reduced whole_step
      do_outer = .false.
      if(stepFailure) firstInnerStep = .true.
      if(firstInnerStep) do_outer = .true.

      if(do_outer)then

        ! if(.not.stepFailure)then
        !   call resizeData(indx_meta(:),indx_data,indx_temp,err=err,message=cmessage)
        !   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
        ! endif

        ! save/recover copies of index variables, temp saved on lastInnerStep, failed starts at lastInnerStep
        ! do iVar=1,size(indx_data%var)
          select case(stepFailure)
            case(.false.); indx_temp = indx_data
            case(.true.);  indx_data = indx_temp
          end select
        ! end do  ! looping through variables

        ! save/recover copies of prognostic variables
        ! do iVar=1,size(prog_data%var)
          select case(stepFailure)
            case(.false.); prog_temp = prog_data
            case(.true.);  prog_data = prog_temp
          end select
        ! end do  ! looping through variables

        ! save/recover copies of diagnostic variables
        ! do iVar=1,size(diag_data%var)
          select case(stepFailure)
            case(.false.); diag_temp = diag_data
            case(.true.);  diag_data = diag_temp
          end select
        ! end do  ! looping through variables

        ! re-assign dimension lengths
        associate(nSnow => indx_data%nSnow, layerType => indx_data%layerType)
            !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    nSnow(iGRU) = 0
    do iLayer=1,size(layerType,1)
      if (layerType(iLayer,iGRU) == iname_snow) nSnow(iGRU) = nSnow(iGRU) + 1
    end do
  !   if (createLayer(iGRU)) then

  !     print*, 694, iGRU
  ! nSnow(iGRU)   = nSnow(iGRU) + 1
  !   end if
  end do
end associate
        nSnow   = maxval(indx_data%nSnow)
        nSoil   = indx_data%nSoil
        nLayers = nSnow+nSoil

        ! *** merge/sub-divide snow layers...
        ! -----------------------------------
        call volicePack(&
                        ! input/output: model data structures
                        doLayerMerge,               & ! intent(in):    flag to force merge of snow layers
                        nGRU, &
                        decisions,            & ! intent(in):    model decisions
                        mpar_data,                  & ! intent(in):    model parameters
                        indx_data,                  & ! intent(inout): type of each layer
                        prog_data,                  & ! intent(inout): model prognostic variables for a local HRU
                        diag_data,                  & ! intent(inout): model diagnostic variables for a local HRU
                        flux_data,                  & ! intent(inout): model fluxes for a local HRU
                        ! output
                        modifiedLayers,             & ! intent(out): flag to denote that layers were modified
                        err,cmessage)                 ! intent(out): error control
        if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if

        ! save the number of snow and soil layers
        nSnow   = maxval(indx_data%nSnow)
        nSoil   = indx_data%nSoil
        nLayers = nSnow + nSoil

        ! compute the indices for the model state variables
        if(firstSubStep .or. modifiedVegState .or. modifiedLayers)then
          call indexState(computeVegFlux,         & ! intent(in):    flag to denote if computing the vegetation flux
                          includeAquifer,         & ! intent(in):    flag to denote if included the aquifer
                          indx_data%nSnow,nSoil,indx_data%nLayers_d,nGRU,    & ! intent(in):    number of snow and soil layers, and total number of layers
                          indx_data,              & ! intent(inout): indices defining model states and layers
                          err,cmessage)             ! intent(out):   error control
          if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
        end if

        ! get enthalpy from temperature if new layering
        if( (enthalpyStateVec .or. computeEnthalpy) .and. modifiedLayers )then
          ! associate local variables with variables in the data structures
          enthalpySnow: associate(&
          nSnow => indx_data%nSnow, &
            ! variables in the snow and soil domains
            mLayerTemp           => prog_data%mLayerTemp              ,& ! temperature (K)
            mLayerEnthTemp       => diag_data%mLayerEnthTemp          ,& ! temperature component of enthalpy (J m-3)
            mLayerEnthalpy       => prog_data%mLayerEnthalpy          ,& ! enthalpy (J m-3)
            mLayerVolFracWat     => prog_data%mLayerVolFracWat        ,& ! volumetric fraction of total water in each snow layer (-)
            mLayerVolFracLiq     => prog_data%mLayerVolFracLiq        ,& ! volumetric fraction of liquid water (-)
            mLayerVolFracIce     => prog_data%mLayerVolFracIce        ,& ! volumetric fraction of ice in each snow layer (-)
            mLayerMatricHead     => prog_data%mLayerMatricHead        ,& ! matric head (m)
            snowfrz_scale => mpar_data%snowfrz_scale, &
            temperature => lookup_data%temperature, psiLiq_int => lookup_data%psiLiq_int, deriv2 => lookup_data%deriv2, &
            ! depth-varying soil parameters
            soil_dens_intr       => mpar_data%soil_dens_intr         ,& ! surface layer  intrinsic soil density (kg m-3)
            vGn_m                => diag_data%scalarVGn_m_m             ,& ! van Genutchen "m" parameter (-)
            vGn_n                => mpar_data%vGn_n                  ,& ! van Genutchen "n" parameter (-)
            vGn_alpha            => mpar_data%vGn_alpha              ,& ! van Genutchen "alpha" parameter (m-1)
            theta_sat            => mpar_data%theta_sat              ,& ! soil porosity (-)
            theta_res            => mpar_data%theta_res               & ! soil residual volumetric water content (-)
            )  ! (associate local variables with model parameters)    
            !$cuf kernel do(1) <<<*,*>>>
            do iGRU=1,nGRU
              do iLayer=1,nSnow(iGRU)
                print*, iLayer, iGRU
                ! if (iLayer .le. nSnow(iGRU)) then
                mLayerVolFracWat(iLayer,iGRU) = mLayerVolFracLiq(iLayer,iGRU) + mLayerVolFracIce(iLayer,iGRU)*(iden_ice/iden_water)
                ! compute enthalpy for snow layers
                call T2enthTemp_snow_d(&
                             snowfrz_scale,             & ! intent(in):  scaling parameter for the snow freezing curve  (K-1)
                             mLayerTemp(iLayer,iGRU),        & ! intent(in):  layer temperature (K)
                             mLayerVolFracWat(iLayer,iGRU),  & ! intent(in):  volumetric total water content (-)
                             mLayerEnthTemp(iLayer,iGRU))      ! intent(out): temperature component of enthalpy of each snow layer (J m-3)
                mLayerEnthalpy(iLayer,iGRU) = mLayerEnthTemp(iLayer,iGRU) - iden_ice * LH_fus * mLayerVolFracIce(iLayer,iGRU)
                ! else  ! looping through snow layers
              end do
            end do  ! looping through soil layers
            threads = dim3(128,1,1)
            blocks = dim3(nGRU/128+1,1,1)
  
            call enthtemp_soilnSoil<<<blocks,threads>>>(nGRU, nSnow,nSoil, &
  use_lookup, &
  vGn_alpha, vGn_n, theta_sat, theta_res, vGn_m, &
  temperature,psiLiq_int,deriv2, &
  mLayerTemp, mLayerEnthTemp, mLayerEnthalpy, &
  mLayerMatricHead, mLayerVolFracIce, mLayerVolFracWat, mLayerVolFracLiq, soil_dens_intr)

          end associate enthalpySnow
        end if ! (need to recalculate enthalpy state variable)

        ! recreate the temporary data structures
        ! NOTE: resizeData(meta, old, new, ..)
        if(modifiedVegState .or. modifiedLayers)then

          ! create temporary data structures for prognostic variables
          prog_temp = prog_data

          ! create temporary data structures for diagnostic variables
          diag_temp = diag_data

          ! create temporary data structures for index variables          
            select case(stepFailure)
              case(.false.); indx_temp = indx_data
              case(.true.);  indx_data = indx_temp
            end select

        endif  ! if modified the states

        ! define the number of state variables
        nState = maxval(indx_data%nState)

        ! *** compute diagnostic variables for each layer...
        ! --------------------------------------------------
        ! NOTE: this needs to be done AFTER volicePack, since layers may have been sub-divided and/or merged, and need to specifically send in canopy depth
        call diagn_evar(&
                        ! input: control variables
                        computeVegFlux,         & ! intent(in): flag to denote if computing the vegetation flux
                        nGRU, &
                        decisions, &
                        diag_data%scalarCanopyDepth, & ! intent(in): canopy depth (m), send in specific value since diag_data may have changed
                        ! input/output: data structures
                        mpar_data,              & ! intent(in):    model parameters
                        indx_data,              & ! intent(in):    model layer indices
                        prog_data,              & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,              & ! intent(inout): model diagnostic variables for a local HRU
                        ! output: error control
                        err,cmessage)             ! intent(out): error control
        if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; end if

        ! *** compute melt of the "snow without a layer"...
        ! -------------------------------------------------
        ! NOTE: forms a surface melt pond, which drains into the upper-most soil layer through the time step
        ! (check for the special case of "snow without a layer")
        ! this pond melts evenly over entire time of maxstep until it gets recomputed because based on SWE when computed
  
        associate(nSnow => indx_data%nSnow, &
          scalarSWE => prog_data%scalarSWE, &
          scalarSnowDepth => prog_data%scalarSnowDepth, &
          scalarSfcMeltPond => prog_data%scalarSfcMeltPond, &
          mLayerTemp => prog_data%mLayerTemp, &
          mLayerDepth => prog_data%mLayerDepth, &
          mLayerVolHtCapBulk => diag_data%mLayerVolHtCapBulk_m)
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
        if(nSnow(iGRU)==0) then
          call implctMelt(&
                          ! input/output: integrated snowpack properties
                          scalarSWE(iGRU),               & ! intent(inout): snow water equivalent (kg m-2)
                          scalarSnowDepth(iGRU),         & ! intent(inout): snow depth (m)
                          scalarSfcMeltPond(iGRU),       & ! intent(inout): surface melt pond (kg m-2)
                          ! input/output: properties of the upper-most soil layer
                          mLayerTemp(nSnow(iGRU)+1,iGRU),        & ! intent(inout): surface layer temperature (K)
                          mLayerDepth(nSnow(iGRU)+1,iGRU),       & ! intent(inout): surface layer depth (m)
                          mLayerVolHtCapBulk(nSnow(iGRU)+1,iGRU)& ! intent(inout): surface layer volumetric heat capacity (J m-3 K-1)
                          ! output: error control
                                                                  ) ! intent(out): error control
          ! if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
        endif
      end do
      end associate

        ! save volumetric ice content at the start of the step
        ! NOTE: used for volumetric loss due to melt-freeze
        if (allocated(mLayerVolFracIceInit)) deallocate(mLayerVolFracIceInit) ! prep for potential size change
        allocate(mLayerVolFracIceInit(nSoil+maxSnowLayers,nGRU)); mLayerVolFracIceInit = prog_data%mLayerVolFracIce
  

        ! make sure have consistent state variables to start, later done in updateVars
        ! associate local variables with information in the data structures
        init: associate(&
          nSnow => indx_data%nSnow, &
          ! depth-varying soil parameters
          soil_dens_intr          => mpar_data%soil_dens_intr        ,& ! intent(in):    [dp]   surface layer  intrinsic soil density (kg m-3)
          vGn_m                   => diag_data%scalarVGn_m_m               ,& ! intent(in):    [dp(:)]  van Genutchen "m" parameter (-)
          vGn_n                   => mpar_data%vGn_n                    ,& ! intent(in):    [dp(:)]  van Genutchen "n" parameter (-)
          vGn_alpha               => mpar_data%vGn_alpha                ,& ! intent(in):    [dp(:)]  van Genutchen "alpha" parameter (m-1)
          theta_sat               => mpar_data%theta_sat                ,& ! intent(in):    [dp(:)]  soil porosity (-)
          theta_res               => mpar_data%theta_res                ,& ! intent(in):    [dp(:)]  soil residual volumetric water content (-)
          ! variables in the vegetation canopy
          scalarCanopyIce         => prog_data%scalarCanopyIce        ,& ! intent(in):    [dp]     mass of ice on the vegetation canopy (kg m-2)
          scalarCanopyLiq         => prog_data%scalarCanopyLiq        ,& ! intent(in):    [dp]     mass of liquid water on the vegetation canopy (kg m-2)
          scalarCanopyWat         => prog_data%scalarCanopyWat        ,& ! intent(out):   [dp]     mass of total water on the vegetation canopy (kg m-2)
          ! variables in the snow and soil domains
          mLayerVolFracIce        => prog_data%mLayerVolFracIce          ,& ! intent(in):    [dp(:)]  volumetric fraction of ice (-)
          mLayerVolFracLiq        => prog_data%mLayerVolFracLiq          ,& ! intent(in):    [dp(:)]  volumetric fraction of liquid water (-)
          mLayerVolFracWat        => prog_data%mLayerVolFracWat          ,& ! intent(out):   [dp(:)]  volumetric fraction of total water (-)
          mLayerMatricHead        => prog_data%mLayerMatricHead          ,& ! intent(in):    [dp(:)]  matric head (m)
          mLayerMatricHeadLiq     => diag_data%mLayerMatricHeadLiq        & ! intent(out):   [dp(:)]  matric potential of liquid water (m)
          ) ! associations to variables in data structures
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU

          ! compute the total water content in the vegetation canopy
          scalarCanopyWat(iGRU) = scalarCanopyLiq(iGRU) + scalarCanopyIce(iGRU)  ! kg m-2

          ! compute the total water content in snow and soil
          ! NOTE: no ice expansion allowed for soil
            do iLayer=1,nSnow(iGRU)
              mLayerVolFracWat(iLayer,iGRU  ) = mLayerVolFracLiq(      iLayer,iGRU  ) + mLayerVolFracIce(      iLayer,iGRU  )*(iden_ice/iden_water)
            end do
            do iLayer=nSnow(iGRU)+1,nLayers
              mLayerVolFracWat(iLayer,iGRU)   = mLayerVolFracLiq(iLayer,iGRU) + mLayerVolFracIce(iLayer,iGRU)
            end do
          end do
  

          ! compute enthalpy of the top soil layer if changed with surface melt pond
            threads = dim3(128,1,1)
            blocks = dim3(nGRU/128+1,1,1)
  
            call enthtemp_soil1<<<blocks,threads>>>(nGRU, nSnow, &
            enthalpyStateVec,computeEnthalpy,use_lookup, &
            prog_data%scalarSWE,&
            vGn_alpha, vGn_n, theta_sat, theta_res, vGn_m, &
            lookup_data%temperature,lookup_data%psiLiq_int,lookup_data%deriv2, &
            prog_data%mLayerTemp, diag_data%mLayerEnthTemp, prog_data%mLayerEnthalpy, &
            mLayerMatricHead, mLayerVolFracIce, soil_dens_intr(1))
          ! if( (enthalpyStateVec .or. computeEnthalpy) .and. nSnow==0 .and. scalarSWE>0._rkind )then
          !   call T2enthTemp_soil2(&
          !             ! use_lookup_d,                                               & ! intent(in):  flag to use the lookup table for soil enthalpy
          !             soil_dens_intr,                                           & ! intent(in):  intrinsic soil density (kg m-3)
          !             vGn_alpha(1),vGn_n(1),theta_sat(1),theta_res(1),vGn_m(1), & ! intent(in):  van Genutchen soil parameters
          !             ! 1_i4b,                                                    & ! intent(in):  index of the control volume within the domain
          !             ! lookup_data,                                              & ! intent(in):  lookup table data structure
          !             ! lookup_data%temperature(:,1,iGRU),lookup_data%psiLiq_int(:,1,iGRU),lookup_data%deriv2(:,1,iGRU), &
          !             realMissing,                                              & ! intent(in):  lower value of integral (not computed)
          !             mLayerTemp(nSnow+1),         & ! intent(in):  surface layer temperature (K)
          !             mLayerMatricHead(1),                                      & ! intent(in):  surface layer matric head (m)
          !             mLayerEnthTemp(nSnow+1))       ! intent(out): temperature component of enthalpy soil layer (J m-3)
          !   mLayerEnthalpy(nSnow+1) = mLayerEnthTemp(nSnow+1) - iden_water * LH_fus * mLayerVolFracIce(nSnow+1)
          ! end if
        ! end do
        ! end associate

          ! compute the liquid water matric potential (m)
          ! NOTE: include ice content as part of the solid porosity - major effect of ice is to reduce the pore size; ensure that effSat=1 at saturation
          ! (from Zhao et al., J. Hydrol., 1997: Numerical analysis of simultaneous heat and mass transfer...)
          associate(mLayerMatricHead => prog_data%mLayerMatricHead, &
            mLayerVolFracLiq => prog_data%mLayerVolFracLiq,nSNow=>indx_data%nSnow,&
            mLayerVolFracIce => prog_data%mLayerVolFracIce, &
            vGn_alpha => mpar_data%vGn_alpha, vGn_n => mpar_data%vGn_n, &
            theta_sat => mpar_data%theta_sat, theta_res => mpar_data%theta_res, &
            vGn_m => diag_data%scalarVGn_m_m, &
            mLayerMatricHeadLiq => diag_data%mLayerMatricHeadLiq)
            !$cuf kernel do(2) <<<*,*>>>
            do iGRU=1,nGRU
          do iSoil=1,nSoil
            call liquidHead(mLayerMatricHead(iSoil,iGRU),mLayerVolFracLiq(nSnow(iGRU)+iSoil,iGRU),mLayerVolFracIce(nSnow(iGRU)+iSoil,iGRU), & ! input:  state variables
                      vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m(iSoil,iGRU),              & ! input:  parameters
                      matricHeadLiq=mLayerMatricHeadLiq(iSoil,iGRU)                                                  & ! output: liquid water matric potential (m)
                      )                                                                    ! output: error control
            ! if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
          end do  ! looping through soil layers (computing liquid water matric potential)
        end do
          end associate

        end associate init


        ! correct increments (if need to redo inner step) and reset increment
        dt_solv = dt_solv - dt_solvInner
        dt_solvInner = 0._rkind
        lastInnerStep = .false.

        ! initialize sublimation sums to average over whole_step
        sumCanopySublimation = 0._rkind
        sumSnowSublimation   = 0._rkind
        sumLatHeatCanopyEvap = 0._rkind
        sumSenHeatCanopy     = 0._rkind
        ! initialize fluxes to average over whole_step (averaged over substep in varSubStep)
        call zero_device_flux(flux_inner)
        innerEffRainfall  = 0._rkind ! mean total effective rainfall over snow
        innerSoilCompress = 0._rkind ! mean total soil compression
        innerBalance = 0._rkind ! mean total balance array
        if (allocated(innerBalanceLayerNrg))  deallocate(innerBalanceLayerNrg)
        allocate(innerBalanceLayerNrg(nLayers,nGRU)); innerBalanceLayerNrg = 0._rkind ! mean total balance of energy in layers
        if (allocated(innerBalanceLayerMass)) deallocate(innerBalanceLayerMass)    ! deallocate if already allocated to permit size change
        allocate(innerBalanceLayerMass(nLayers,nGRU)); innerBalanceLayerMass = 0._rkind ! mean total balance of mass in layers
        sumStepSize= 0._rkind ! initialize the sum of the step sizes

      endif ! (do_outer loop)

      ! *** solve model equations...
      ! ----------------------------
      ! save input step
      dtSave = whole_step

      ! get the new solution
      call opSplittin(&
                      ! input: model control
                      nSnow,                                  & ! intent(in):    number of snow layers
                      nSoil,                                  & ! intent(in):    number of soil layers
                      nLayers,                                & ! intent(in):    total number of layers
                      nState,                                 & ! intent(in):    total number of layers
                      nGRU, &
                      dt_sub,                                 & ! intent(in):    length of the model sub-step
                      whole_step,                             & ! intent(in):    length of whole step for surface drainage and average flux
                      (dt_solv<whole_step),                   & ! intent(in):    logical flag to denote the first loop of the whole_step in a data_step
                      firstInnerStep,                         & ! intent(in):    flag to denote if the first time step in maxstep subStep
                      computeVegFlux,                         & ! intent(in):    logical flag to compute fluxes within the vegetation canopy
                      ! input/output: data structures
                      type_data,                              & ! intent(in):    type of vegetation and soil
                      attr_data,                              & ! intent(in):    spatial attributes
                      forc_data,                              & ! intent(in):    model forcing data
                      mpar_data,                              & ! intent(in):    model parameters
                      indx_data,                              & ! intent(inout): index data
                      prog_data,                              & ! intent(inout): model prognostic variables for a local HRU
                      diag_data,                              & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,                              & ! intent(inout): model fluxes for a local HRU
                      bvar_data,                              & ! intent(in):    model variables for the local basin
                      lookup_data,                            & ! intent(in):    lookup tables
                      model_decisions,                        & ! intent(in):    model decisions
                      decisions, veg_param, &
                      ! output: model control
                      dtMultiplier,                           & ! intent(out):   substep multiplier (-)
                      tooMuchMelt,                            & ! intent(out):   flag to denote that ice is insufficient to support melt
                      stepFailure,                            & ! intent(out):   flag to denote that the coupled step failed
                      ixSolution,                             & ! intent(out):   solution method used in this iteration
                      mean_step_dt_sub,                       & ! intent(out):   mean solution step for the sub-step
                      err,cmessage)                             ! intent(out):   error code and error message
      ! check for all errors (error recovery within opSplittin)
      if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

      ! process the flag for too much melt
      if(tooMuchMelt)then
        stepFailure  = .true.
        doLayerMerge = .true.
      else
        doLayerMerge = .false.
      endif

      ! handle special case of the step failure
      ! NOTE: need to revert back to the previous state vector that we were happy with and reduce the time step
      ! TODO: ask isn't this what the actors program does without the code block below
      if(stepFailure)then
        ! halve whole_step, for more frequent outer loop updates
        whole_step = dtSave/2._rkind
        ! check that the step is not tiny
        if(whole_step < minstep)then
          print*,ixSolution
          print*, 'dtSave, dt_sub', dtSave, whole_step
          message=trim(message)//'length of the coupled step is below the minimum step length'
          err=20; return
        endif
        ! try again, restart step
        deallocate(mLayerVolFracIceInit)
        deallocate(innerBalanceLayerNrg)
        deallocate(innerBalanceLayerMass)
        cycle substeps
      endif

      associate(&
        scalarCanopySublimation => flux_data%scalarCanopySublimation, &
        scalarSnowSublimation => flux_data%scalarSnowSublimation, &
        scalarLatHeatCanopyEvap => flux_data%scalarLatHeatCanopyEvap, &
        scalarSenHeatCanopy => flux_data%scalarSenHeatCanopy)
      ! increment sublimation sums
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
          sumCanopySublimation(iGRU) = sumCanopySublimation(iGRU) + dt_sub*scalarCanopySublimation(iGRU)
          sumSnowSublimation(iGRU)   = sumSnowSublimation(iGRU)   + dt_sub*scalarSnowSublimation(iGRU)
          sumLatHeatCanopyEvap(iGRU) = sumLatHeatCanopyEvap(iGRU) + dt_sub*scalarLatHeatCanopyEvap(iGRU)
          sumSenHeatCanopy(iGRU)     = sumSenHeatCanopy(iGRU)     + dt_sub*scalarSenHeatCanopy(iGRU)
        end do
      end associate


      ! update the step size sum
      sumStepSize = sumStepSize + mean_step_dt_sub

      ! update first step and first and last inner steps
      firstSubStep = .false.
      firstInnerStep = .false.
      if(dt_solvInner + dt_sub >= whole_step) lastInnerStep = .true.
      if(dt_solv + dt_sub >= data_step-verySmall) lastInnerStep = .true.

      ! check if on outer loop
      do_outer = .false.
      if(lastInnerStep) do_outer = .true.

      if(do_outer)then

        ! ***  remove ice due to sublimation and freeze calculations...
        ! NOTE: In the future this should be moved into the solver, makes a big difference
        ! --------------------------------------------------------------
        sublime: associate(&
          mLayerMeltFreeze        => diag_data%mLayerMeltFreeze,           & ! melt-freeze in each snow and soil layer (kg m-3)
          scalarCanopyLiq         => prog_data%scalarCanopyLiq,         & ! liquid water stored on the vegetation canopy (kg m-2)
          scalarCanopyIce         => prog_data%scalarCanopyIce,         & ! ice          stored on the vegetation canopy (kg m-2)
          scalarCanopyWat         => prog_data%scalarCanopyWat,         & ! canopy ice content (kg m-2)
          mLayerVolFracIce        => prog_data%mLayerVolFracIce,           & ! volumetric fraction of ice in the snow+soil domain (-)
          mLayerVolFracLiq        => prog_data%mLayerVolFracLiq,           & ! volumetric fraction of liquid water in the snow+soil domain (-)
          mLayerVolFracWat        => prog_data%mLayerVolFracWat,           & ! volumetric fraction of total water (-)
          mLayerDepth             => prog_data%mLayerDepth                 & ! depth of each snow+soil layer (m)
          ) ! associations to variables in data structures

          associate(nSnow=>indx_data%nSnow&
            ! nLayers => indx_data_d%nLayers_d, &
            )
          ! compute the melt in each snow and soil layer
            !$cuf kernel do(1) <<<*,*>>>
            do iGRU=1,nGRU
              do iLayer=1,nSnow(iGRU)
                mLayerMeltFreeze(iLayer,iGRU) = -( mLayerVolFracIce(iLayer,iGRU) - mLayerVolFracIceInit(iLayer,iGRU) )*iden_ice
              end do
              do iLayer=nSnow(iGRU)+1,nLayers
                mLayerMeltFreeze(iLayer,iGRU) = -( mLayerVolFracIce(iLayer,iGRU) - mLayerVolFracIceInit(iLayer,iGRU) )*iden_water
                ! print*, iLayer, mLayerMeltFreeze(iLayer,iGRU)
              end do
            end do
            end associate
            deallocate(mLayerVolFracIceInit)

          ! * compute change in canopy ice content due to sublimation...
          ! ------------------------------------------------------------
          if(computeVegFlux)then

            associate( &
              scalarCanopyTemp => prog_data%scalarCanopyTemp, &
              scalarCanopyEnthTemp => diag_data%scalarCanopyEnthTemp, &
              scalarCanopyEnthalpy => prog_data%scalarCanopyEnthalpy, &
              ! scalarCanopyDepth => diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1), &
              canopyDepth => diag_data%scalarCanopyDepth, &
              ! scalarCanopyTemp => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1), &
              ! scalarCanopyEnthTemp => diag_data%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1), &
              ! scalarCanopyEnthalpy => prog_data%var(iLookPROG%scalarCanopyEnthalpy)%dat(1) &
              specificHeatVeg => mpar_data%specificHeatVeg, &
              maxMassVegetation => mpar_data%maxMassVegetation, &
              snowfrz_scale => mpar_data%snowfrz_scale &
              )
              !$cuf kernel do(1) <<<*,*>>>
              do iGRU=1,nGRU
            ! remove mass of ice on the canopy
            scalarCanopyIce(iGRU) = scalarCanopyIce(iGRU) + sumCanopySublimation(iGRU)

            ! if removed all ice, take the remaining sublimation from water
            if(scalarCanopyIce(iGRU) < 0._rkind)then
              scalarCanopyLiq(iGRU) = scalarCanopyLiq(iGRU) + scalarCanopyIce(iGRU)
              scalarCanopyIce(iGRU) = 0._rkind
            endif

            ! modify fluxes and mean fluxes if there is insufficient canopy water to support the converged sublimation rate over the whole time step
            if(scalarCanopyLiq(iGRU) < 0._rkind)then
              ! --> superfluous sublimation flux
              superflousSub = -scalarCanopyLiq(iGRU)/whole_step  ! kg m-2 s-1
              superflousNrg = superflousSub*LH_sub     ! W m-2 (J m-2 s-1)
              ! --> update fluxes and states
              sumCanopySublimation(iGRU) = sumCanopySublimation(iGRU) + superflousSub*whole_step
              sumLatHeatCanopyEvap(iGRU) = sumLatHeatCanopyEvap(iGRU) + superflousNrg*whole_step
              sumSenHeatCanopy(iGRU)     = sumSenHeatCanopy(iGRU)     - superflousNrg*whole_step
              scalarCanopyLiq(iGRU)      = 0._rkind
            endif

            ! update water
            scalarCanopyWat(iGRU) = scalarCanopyLiq(iGRU) + scalarCanopyIce(iGRU)


            ! print*, scalarCanopyWat(iGRU), scalarCanopyLiq(iGRU), scalarCanopyIce(iGRU)

            if(enthalpyStateVec .or. computeEnthalpy)then ! recompute enthalpy of the canopy if changed water and ice content
              call T2enthTemp_veg_d(&
                          canopyDepth(iGRU),    & ! intent(in): canopy depth (m), send in specific value since diag_data may have changed
                          specificHeatVeg,                                      & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                          maxMassVegetation,                                    & ! intent(in): maximum mass of vegetation (kg m-2)
                          snowfrz_scale,                                        & ! intent(in): scaling parameter for the snow freezing curve  (K-1)
                          scalarCanopyTemp(iGRU),     & ! intent(in): canopy temperature (K)
                          scalarCanopyWat(iGRU),                                      & ! intent(in): canopy water content (kg m-2)
                          scalarCanopyEnthTemp(iGRU))   ! intent(out): temperature component of enthalpy of the vegetation canopy (J m-3)
              scalarCanopyEnthalpy(iGRU) = scalarCanopyEnthTemp(iGRU) - LH_fus * scalarCanopyIce(iGRU)/ canopyDepth(iGRU)
            endif
          end do
            end associate
  
          end if  ! (if computing the vegetation flux)

          ! * compute change in ice content of the top snow layer due to sublimation 
          !   and account for compaction and cavitation in the snowpack...
          ! ------------------------------------------------------------------------
          call computSnowDepth_d(&
                    whole_step_d,                               & ! intent(in)
                    indx_data%nSnow,                                    & ! intent(in)
                    nGRU, &
                    sumSnowSublimation,            & ! intent(in)
                    mLayerVolFracLiq,                         & ! intent(inout)
                    mLayerVolFracIce,                         & ! intent(inout)
                    prog_data%mLayerTemp,  & ! intent(in)
                    mLayerMeltFreeze,                         & ! intent(in)
                    mpar_data,                                & ! intent(in)
                    ! output
                    tooMuchSublim,                            & ! intent(out): flag to denote that there was too much sublimation in a given time step
                    mLayerDepth,                              & ! intent(inout)
                    ! error control
                    err,message)                                ! intent(out):   error control
          if(err/=0)then; err=55; return; end if

          ! process the flag for too much sublimation
          if(tooMuchSublim)then
            stepFailure  = .true.
            doLayerMerge = .true.
          else
            doLayerMerge = .false.
          endif

          ! handle special case of the step failure
          ! NOTE: need to revert back to the previous state vector that we were happy with and reduce the time step
          if(stepFailure)then
            ! halve whole_step, for more frequent outer loop updates
            whole_step = dtSave/2._rkind
            ! check that the step is not tiny
            if(whole_step < minstep)then
              print*,ixSolution
              print*, 'dtSave, dt_sub', dtSave, whole_step
              message=trim(message)//'length of the coupled step is below the minimum step length'
              err=20; return
            endif
            ! try again, restart step (at end inner step)
            deallocate(innerBalanceLayerNrg)
            deallocate(innerBalanceLayerMass)
            cycle substeps
          endif

          ! update coordinate variables
          call calcHeight_d(&
          nGRU, &
                     ! input/output: data structures
                     indx_data,   & ! intent(in): layer type
                     prog_data,   & ! intent(inout): model variables for a local HRU
                     ! output: error control
                     err,cmessage)
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if


          ! recompute snow depth, SWE, and layer water
          associate(nSnow => indx_data%nSnow, &
            scalarSnowDepth => prog_data%scalarSnowDepth, &
            scalarSWE => prog_data%scalarSWE, &
            mLayerTemp => prog_data%mLayerTemp, &
            mLayerEnthTemp => diag_data%mLayerEnthTemp, &
            mLayerEnthalpy => prog_data%mLayerEnthalpy, &
            snowfrz_scale => mpar_data%snowfrz_scale)
            !$cuf kernel do(1) <<<*,*>>>
            do iGRU=1,nGRU
          if(nSnow(iGRU) > 0)then
            scalarSnowDepth(iGRU) = 0._rkind
            scalarSWE(iGRU) = 0._rkind
            do iLayer=1,nSnow(iGRU)
              scalarSnowDepth(iGRU) = scalarSnowDepth(iGRU) + mLayerDepth(iLayer,iGRU)
              scalarSWE(iGRU)       = (mLayerVolFracLiq(iLayer,iGRU)*iden_water &
              + mLayerVolFracIce(iLayer,iGRU)*iden_ice) * mLayerDepth(iLayer,iGRU)
              mLayerVolFracWat(iLayer,iGRU) = mLayerVolFracLiq(iLayer,iGRU) + mLayerVolFracIce(iLayer,iGRU)*iden_ice/iden_water
            end do
            if(enthalpyStateVec .or. computeEnthalpy)then ! recompute enthalpy of layers if changed water and ice content
              do iLayer=1,nSnow(iGRU)
                call T2enthTemp_snow_d(&
                             snowfrz_scale,                                       & ! intent(in):  scaling parameter for the snow freezing curve  (K-1)
                             mLayerTemp(iLayer,iGRU),     & ! intent(in):  layer temperature (K)
                             mLayerVolFracWat(iLayer,iGRU),                            & ! intent(in):  volumetric total water content (-)
                             mLayerEnthTemp(iLayer,iGRU))   ! intent(out): temperature component of enthalpy of each snow layer (J m-3)
                mLayerEnthalpy(iLayer,iGRU) = mLayerEnthTemp(iLayer,iGRU) - iden_ice * LH_fus * mLayerVolFracIce(iLayer,iGRU)
              end do  ! looping through snow layers
            endif
          endif
        end do
        end associate
        end associate sublime

        ! increment change in storage associated with the surface melt pond (kg m-2)
        associate(nSnow => indx_data%nSnow, &
          scalarSfcMeltPond => prog_data%scalarSfcMeltPond)
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
        if(nSnow(iGRU)==0) sfcMeltPond(iGRU) = sfcMeltPond(iGRU) + scalarSfcMeltPond(iGRU)
          end do
          end associate

      endif ! (do_outer loop)

      ! ****************************************************************************************************
      ! *** END MAIN SOLVER ********************************************************************************
      ! ****************************************************************************************************

      ! increment mean fluxes, soil compression, and effective rainfall, reset on whole_step
      dt_wght = dt_sub/whole_step ! define weight applied to each sub-step
      ! do iVar=1,size(flux_meta)
      !   if (flux_meta(iVar)%vartype == iLookVarType%scalarv) flux_inner%var(iVar)%dat(:)    = flux_inner%var(iVar)%dat(:) + flux_data%var(iVar)%dat(:)*dt_wght
      ! end do
      call update_flux_inner(flux_inner,flux_data,nGRU,dt_wght)

      associate(mLayerCompress => diag_data%mLayerCompress_m)
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
          do iLayer=1,nSoil
            innerSoilCompress(iLayer,iGRU) = innerSoilCompress(iLayer,iGRU) + mLayerCompress(iLayer,iGRU)*dt_wght
          end do
        end do
      end associate
      associate(nSnow => indx_data%nSnow, scalarThroughfallRain=>flux_data%scalarThroughfallRain, scalarCanopyLiqDrainage=>flux_data%scalarCanopyLiqDrainage)
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
          if (nSnow(iGRU)>0) innerEffRainfall(iGRU) = innerEffRainfall(iGRU) + ( scalarThroughfallRain(iGRU) + scalarCanopyLiqDrainage(iGRU) )*dt_wght
        end do
      end associate


      ! sum the balance of energy and water per state
      associate(&
        balanceCasNrg => diag_data%balanceCasNrg, &
        balanceVegNrg => diag_data%balanceVegNrg, &
        canopyDepth => diag_data%scalarCanopyDepth, &
        balanceVegMass => diag_data%balanceVegMass, &
        balanceAqMass => diag_data%balanceAqMass, &
        balanceLayerNrg => diag_data%balanceLayerNrg, &
        balanceLayerMass => diag_data%balanceLayerMass)
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
          if(computeVegFlux)then
            innerBalance(1,iGRU) = innerBalance(1,iGRU) + balanceCasNrg(iGRU)*dt_wght ! W m-3
            innerBalance(2,iGRU) = innerBalance(2,iGRU) + balanceVegNrg(iGRU)*dt_wght ! W m-3
            innerBalance(3,iGRU) = innerBalance(3,iGRU) + balanceVegMass(iGRU)*dt_wght/canopyDepth(iGRU)  ! kg m-3 s-1
            bal_veg = .true.
          endif
          innerBalance(4,iGRU) = innerBalance(4,iGRU) + balanceAqMass(iGRU)*dt_wght * iden_water  ! kg m-2 s-1 (no depth to aquifer)
        end do
        !$cuf kernel do(2) <<<*,*>>>
        do iGRU=1,nGRU
          do iLayer=1,nLayers
            innerBalanceLayerNrg(iLayer,iGRU) = innerBalanceLayerNrg(iLayer,iGRU) + balanceLayerNrg(iLayer,iGRU)*dt_wght ! W m-3
            innerBalanceLayerMass(iLayer,iGRU) = innerBalanceLayerMass(iLayer,iGRU) + balanceLayerMass(iLayer,iGRU)*dt_wght * iden_water ! kg m-3 s-1
          end do
        end do

      ! save balance of energy and water per snow+soil layer after inner step, since can change nLayers with outer steps
        !$cuf kernel do(2) <<<*,*>>>
        do iGRU=1,nGRU
          do iLayer=1,nLayers
            balanceLayerNrg(iLayer,iGRU) = innerBalanceLayerNrg(iLayer,iGRU)
            balanceLayerMass(iLayer,iGRU) = innerBalanceLayerMass(iLayer,iGRU)
          end do
        end do
      end associate

      ! compute the balance of energy and water per entire snow and soil domain, in W m-3 and kg m-2 s-1 respectively
      associate(balanceSnowMass => diag_data%balanceSnowMass, &
        balanceSnowNrg => diag_data%balanceSnowNrg, &
        balanceSoilMass => diag_data%balanceSoilMass, &
        balanceSoilNrg => diag_data%balanceSoilNrg, &
        layerType => indx_data%layerType, &
        mLayerDepth => prog_data%mLayerDepth, nSnow => indx_data%nSnow)
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
          balanceSnowNrg(iGRU) = 0._rkind
          balanceSoilNrg(iGRU) = 0._rkind
          balanceSnowMass(iGRU) = 0._rkind
          balanceSoilMass(iGRU) = 0._rkind
          layerDepthSnow(iGRU) = 0._rkind
          layerDepthSOil(iGRU) = 0._rkind
          do iLayer=1,nSnow(iGRU)
            layerDepthSnow(iGRU) = layerDepthSnow(iGRU) + mLayerDepth(iLayer,iGRU)
          end do
          do iLayer=1,nSoil
            layerDepthSoil(iGRU) = layerDepthSoil(iGRU) + mLayerDepth(iLayer+nSnow(iGRU),iGRU)
          end do
          do iLayer=1,nLayers
            select case (layerType(iLayer,iGRU))
              case (iname_snow)
                lyr_wght = mLayerDepth(iLayer,iGRU) / layerDepthSnow(iGRU)
                balanceSnowNrg(iGRU)  = balanceSnowNrg(iGRU) + innerBalanceLayerNrg(iLayer,iGRU)*lyr_wght
                balanceSnowMass(iGRU) = balanceSnowMass(iGRU) + innerBalanceLayerMass(iLayer,iGRU)*lyr_wght
                bal_snow = .true.
              case (iname_soil)
                lyr_wght = mLayerDepth(iLayer,iGRU) / layerDepthSoil(iGRU)
                balanceSoilNrg(iGRU)  = balanceSoilNrg(iGRU) + innerBalanceLayerNrg(iLayer,iGRU)*lyr_wght
                balanceSoilMass(iGRU) = balanceSoilMass(iGRU) + innerBalanceLayerMass(iLayer,iGRU)*lyr_wght
                bal_soil = .true.
            end select
          end do
        end do
      end associate
      if (model_decisions(iLookDECISIONS%groundwatr)%iDecision == bigBucket) bal_aq = .true. ! aquifer does not change existance with time steps

      if(do_outer)then
        deallocate(innerBalanceLayerNrg)
        deallocate(innerBalanceLayerMass)
      endif

      ! increment sub-step accepted step
      dt_solvInner = dt_solvInner + dt_sub
      dt_solv = dt_solv + dt_sub

      ! update first and last inner steps if did successful lastInnerStep, increment fluxes and flux variables over data_step
      if (lastInnerStep)then
        firstInnerStep = .true.
        lastInnerStep = .false.
        dt_solvInner = 0._rkind

        dt_wght = whole_step/data_step ! define weight applied to each sub-step
        call update_flux_inner(flux_mean,flux_inner,nGRU,dt_wght)

        associate(&
          balanceSnowNrg => diag_data%balanceSnowNrg, &
          balanceSnowMass => diag_data%balanceSnowMass, &
          balanceSoilNrg => diag_data%balanceSoilNrg, &
          balanceSoilMass => diag_data%balanceSoilMass, &
          scalarCanopySublimation => flux_mean%scalarCanopySublimation, &
          scalarLatHeatCanopyEvap => flux_mean%scalarLatHeatCanopyEvap, &
          scalarSenHeatCanopy => flux_mean%scalarSenHeatCanopy &
          )
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
            meanCanopySublimation(iGRU) = meanCanopySublimation(iGRU) + sumCanopySublimation(iGRU)/data_step
            meanLatHeatCanopyEvap(iGRU) = meanLatHeatCanopyEvap(iGRU) + sumLatHeatCanopyEvap(iGRU)/data_step
            meanSenHeatCanopy(iGRU)     = meanSenHeatCanopy(iGRU)     + sumSenHeatCanopy(iGRU)/data_step
            do iLayer=1,nSoil
              meanSoilCompress(iLayer,iGRU) = meanSoilCompress(iLayer,iGRU) + innerSoilCompress(iLayer,iGRU)*dt_wght
            end do
            meanBalance(1,iGRU) = meanBalance(1,iGRU) + innerBalance(1,iGRU)*dt_wght
            meanBalance(2,iGRU) = meanBalance(2,iGRU) + innerBalance(2,iGRU)*dt_wght
            meanBalance(3,iGRU) = meanBalance(3,iGRU) + innerBalance(3,iGRU)*dt_wght
            meanBalance(4,iGRU) = meanBalance(4,iGRU) + innerBalance(4,iGRU)*dt_wght
            meanBalance(5,iGRU) = meanBalance(5,iGRU) + balanceSnowNrg(iGRU)*dt_wght
            meanBalance(6,iGRU) = meanBalance(6,iGRU) + balanceSoilNrg(iGRU)*dt_wght
            meanBalance(7,iGRU) = meanBalance(7,iGRU) + balanceSnowMass(iGRU)*dt_wght
            meanBalance(8,iGRU) = meanBalance(8,iGRU) + balanceSoilMass(iGRU)*dt_wght

            effRainfall(iGRU) = effRainfall(iGRU) + innerEffRainfall(iGRU)*dt_wght
            scalarCanopySublimation(iGRU) = meanCanopySublimation(iGRU)
            scalarLatHeatCanopyEvap(iGRU) = meanLatHeatCanopyEvap(iGRU)
            scalarSenHeatCanopy(iGRU)     = meanSenHeatCanopy(iGRU)
    
          end do
        end associate  

        ! add mean step size for the data_step to the total step size sum
        diag_data%meanStepSize =  diag_data%meanStepSize + sumStepSize
      endif

      ! save the time step to initialize the subsequent step
      if(dt_solv<data_step .or. nsub==1) dt_init = dt_sub

      ! check
      if(globalPrintFlag)&
      write(*,'(a,1x,3(f18.5,1x))') 'dt_sub, dt_solv, data_step: ', dt_sub, dt_solv, data_step

      nsub_success = nsub_success + 1
      ! check that we have completed the sub-step
      if(dt_solv >= data_step-verySmall) then
        exit substeps
      endif

      ! adjust length of the sub-step (make sure that we don't exceed the step)
      dt_sub = min(data_step - dt_solv, dt_sub)

    end do  substeps ! (sub-step loop)
    diag_data%meanStepSize = diag_data%meanStepSize/nsub_success



    ! *** add snowfall to the snowpack...
    ! -----------------------------------
    ! add new snowfall to the snowpack
    ! NOTE: This needs to be done AFTER the call to canopySnow, since throughfall and unloading are computed in canopySnow
    call newsnwfall(&
                  ! input: model control
                  data_step,                                                 & ! time step (seconds)
                  nGRU, &
                  indx_data%nSnow,                                               & ! logical flag if snow layers exist
                  mpar_data%snowfrz_scale,                                             & ! freeezing curve parameter for snow (K-1)
                  ! input: diagnostic scalar variables
                  diag_data%scalarSnowfallTemp,        & ! computed temperature of fresh snow (K)
                  diag_data%scalarNewSnowDensity,      & ! computed density of new snow (kg m-3)
                  flux_data%scalarThroughfallSnow,     & ! throughfall of snow through the canopy (kg m-2 s-1)
                  flux_data%scalarCanopySnowUnloading, & ! unloading of snow from the canopy (kg m-2 s-1)
                  ! input/output: state variables
                  prog_data%scalarSWE,                 & ! SWE (kg m-2)
                  prog_data%scalarSnowDepth,           & ! total snow depth (m)
                  prog_data%mLayerTemp,                & ! temperature of the top layer (K)
                  prog_data%mLayerDepth,               & ! depth of the top layer (m)
                  prog_data%mLayerVolFracIce,          & ! volumetric fraction of ice of the top layer (-)
                  prog_data%mLayerVolFracLiq,          & ! volumetric fraction of liquid water of the top layer (-)
                  ! output: error control
                  err,cmessage)                                                ! error control
    if(err/=0)then; err=30; message=trim(message)//trim(cmessage); return; end if


   associate(scalarSnowDepth => prog_data%scalarSnowDepth, &
    mLayerDepth => prog_data%mLayerDepth, &
    scalarSWE => prog_data%scalarSWE, &
    nSnow => indx_data%nSnow, &
    mLayerVolFracLiq => prog_data%mLayerVolFracLiq, &
    mLayerVolFracIce => prog_data%mLayerVolFracIce, &
    mLayerVolFracWat => prog_data%mLayerVolFracWat, &
    mLayerTemp => prog_data%mLayerTemp, &
    mlayerEnthTemp => diag_data%mLayerEnthTemp, &
    snowfrz_scale => mpar_data%snowfrz_scale, &
    mLayerEnthalpy => prog_data%mLayerEnthalpy)
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
    ! recompute snow depth, SWE, and top layer water
    if(nSnow(iGRU) > 0)then
      scalarSnowDepth(iGRU) = 0
      scalarSWE(iGRU) = 0
      do iLayer=1,nSnow(iGRU)
        scalarSnowDepth(iGRU) = scalarSnowDepth(iGRU) + mLayerDepth(iLayer,iGRU)
        scalarSWE(iGRU) = scalarSWE(iGRU) + (mLayerVolFracLiq(iLayer,iGRU)*iden_water + &
                                                              mLayerVolFracIce(iLayer,iGRU)*iden_ice) &
                                                            * mLayerDepth(iLayer,iGRU)
      end do
      mLayerVolFracWat(1,iGRU) = mLayerVolFracLiq(1,iGRU) &
                                                        + mLayerVolFracIce(1,iGRU)*iden_ice/iden_water
      if(enthalpyStateVec .or. computeEnthalpy)then ! compute enthalpy of the top snow layer
        call T2enthTemp_snow_d(&
                       snowfrz_scale,                                     & ! intent(in):  scaling parameter for the snow freezing curve  (K-1)
                       mLayerTemp(1,iGRU),        & ! temperature of the top layer (K)
                       mLayerVolFracWat(1,iGRU),  & ! intent(in):  volumetric total water content (-)
                       mLayerEnthTemp(1,iGRU))      ! intent(out): temperature component of enthalpy of each snow layer (J m-3)
        mLayerEnthalpy(1,iGRU) = mLayerEnthTemp(1,iGRU) - iden_ice * LH_fus * mLayerVolFracIce(1,iGRU)
      end if
    end if
  end do
  end associate

      ! re-assign dimension lengths
        associate(nSnow => indx_data%nSnow, layerType => indx_data%layerType)
            !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    nSnow(iGRU) = 0
    do iLayer=1,size(layerType,1)
      if (layerType(iLayer,iGRU) == iname_snow) nSnow(iGRU) = nSnow(iGRU) + 1
    end do
  !   if (createLayer(iGRU)) then

  !     print*, 694, iGRU
  ! nSnow(iGRU)   = nSnow(iGRU) + 1
  !   end if
  end do
end associate
        nSnow   = maxval(indx_data%nSnow)
        nSoil   = indx_data%nSoil
        nLayers = nSnow+nSoil

    ! update coordinate variables
    call calcHeight_d(&
    nGRU, &
                    ! input/output: data structures
                    indx_data,   & ! intent(in): layer type
                    prog_data,   & ! intent(inout): model variables for a local HRU
                    ! output: error control
                    err,cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if



    ! overwrite flux_data and soil compression with the timestep-average value (returns timestep-average fluxes for scalar variables)
    do iVar=1,size(flux_meta)
            if (flux_meta(iVar)%vartype==iLookVarType%scalarv) then
select case(iVar)
    case(iLookFLUX%scalarCanairNetNrgFlux); flux_data%scalarCanairNetNrgFlux = flux_mean%scalarCanairNetNrgFlux
    case(iLookFLUX%scalarCanopyNetNrgFlux); flux_data%scalarCanopyNetNrgFlux = flux_mean%scalarCanopyNetNrgFlux
    case(iLookFLUX%scalarGroundNetNrgFlux); flux_data%scalarGroundNetNrgFlux = flux_mean%scalarGroundNetNrgFlux
    case(iLookFLUX%scalarCanopyNetLiqFlux); flux_data%scalarCanopyNetLiqFlux = flux_mean%scalarCanopyNetLiqFlux
    case(iLookFLUX%scalarRainfall); flux_data%scalarRainfall = flux_mean%scalarRainfall
    case(iLookFLUX%scalarSnowfall); flux_data%scalarSnowfall = flux_mean%scalarSnowfall
    case(iLookFLUX%spectralIncomingDirect); flux_data%spectralIncomingDirect = flux_mean%spectralIncomingDirect
    case(iLookFLUX%spectralIncomingDiffuse); flux_data%spectralIncomingDiffuse = flux_mean%spectralIncomingDiffuse
    case(iLookFLUX%scalarCanopySunlitPAR); flux_data%scalarCanopySunlitPAR = flux_mean%scalarCanopySunlitPAR
    case(iLookFLUX%scalarCanopyShadedPAR); flux_data%scalarCanopyShadedPAR = flux_mean%scalarCanopyShadedPAR
    case(iLookFLUX%spectralBelowCanopyDirect); flux_data%spectralBelowCanopyDirect = flux_mean%spectralBelowCanopyDirect
    case(iLookFLUX%spectralBelowCanopyDiffuse); flux_data%spectralBelowCanopyDiffuse = flux_mean%spectralBelowCanopyDiffuse
    case(iLookFLUX%scalarBelowCanopySolar); flux_data%scalarBelowCanopySolar = flux_mean%scalarBelowCanopySolar
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
    case(iLookFLUX%scalarGroundResistance); flux_data%scalarGroundResistance = flux_mean%scalarGroundResistance
    case(iLookFLUX%scalarCanopyResistance); flux_data%scalarCanopyResistance = flux_mean%scalarCanopyResistance
    case(iLookFLUX%scalarLeafResistance); flux_data%scalarLeafResistance = flux_mean%scalarLeafResistance
    case(iLookFLUX%scalarSoilResistance); flux_data%scalarSoilResistance = flux_mean%scalarSoilResistance
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
    case(iLookFLUX%scalarStomResistSunlit); flux_data%scalarStomResistSunlit = flux_mean%scalarStomResistSunlit
    case(iLookFLUX%scalarStomResistShaded); flux_data%scalarStomResistShaded = flux_mean%scalarStomResistShaded
    case(iLookFLUX%scalarPhotosynthesisSunlit); flux_data%scalarPhotosynthesisSunlit = flux_mean%scalarPhotosynthesisSunlit
    case(iLookFLUX%scalarPhotosynthesisShaded); flux_data%scalarPhotosynthesisShaded = flux_mean%scalarPhotosynthesisShaded
    case(iLookFLUX%scalarCanopyTranspiration); flux_data%scalarCanopyTranspiration = flux_mean%scalarCanopyTranspiration
    case(iLookFLUX%scalarCanopyEvaporation); flux_data%scalarCanopyEvaporation = flux_mean%scalarCanopyEvaporation
    case(iLookFLUX%scalarGroundEvaporation); flux_data%scalarGroundEvaporation = flux_mean%scalarGroundEvaporation
    case(iLookFLUX%mLayerTranspire); flux_data%mLayerTranspire_m = flux_mean%mLayerTranspire_m
    case(iLookFLUX%scalarThroughfallSnow); flux_data%scalarThroughfallSnow = flux_mean%scalarThroughfallSnow
    case(iLookFLUX%scalarThroughfallRain); flux_data%scalarThroughfallRain = flux_mean%scalarThroughfallRain
    case(iLookFLUX%scalarCanopySnowUnloading); flux_data%scalarCanopySnowUnloading = flux_mean%scalarCanopySnowUnloading
    case(iLookFLUX%scalarCanopyLiqDrainage); flux_data%scalarCanopyLiqDrainage = flux_mean%scalarCanopyLiqDrainage
    case(iLookFLUX%scalarCanopyMeltFreeze); flux_data%scalarCanopyMeltFreeze = flux_mean%scalarCanopyMeltFreeze
    case(iLookFLUX%iLayerConductiveFlux); flux_data%iLayerConductiveFlux_m = flux_mean%iLayerConductiveFlux_m
    case(iLookFLUX%iLayerAdvectiveFlux); flux_data%iLayerAdvectiveFlux_m = flux_mean%iLayerAdvectiveFlux_m
    case(iLookFLUX%iLayerNrgFlux); flux_data%iLayerNrgFlux_m = flux_mean%iLayerNrgFlux_m
    case(iLookFLUX%mLayerNrgFlux); flux_data%mLayerNrgFlux_m = flux_mean%mLayerNrgFlux_m
    case(iLookFLUX%scalarSnowDrainage); flux_data%scalarSnowDrainage = flux_mean%scalarSnowDrainage
    case(iLookFLUX%iLayerLiqFluxSnow); flux_data%iLayerLiqFluxSnow_m = flux_mean%iLayerLiqFluxSnow_m
    case(iLookFLUX%mLayerLiqFluxSnow); flux_data%mLayerLiqFluxSnow_m = flux_mean%mLayerLiqFluxSnow_m
    case(iLookFLUX%scalarRainPlusMelt); flux_data%scalarRainPlusMelt = flux_mean%scalarRainPlusMelt
    case(iLookFLUX%scalarMaxInfilRate); flux_data%scalarMaxInfilRate = flux_mean%scalarMaxInfilRate
    case(iLookFLUX%scalarInfiltration); flux_data%scalarInfiltration = flux_mean%scalarInfiltration
    case(iLookFLUX%scalarExfiltration); flux_data%scalarExfiltration = flux_mean%scalarExfiltration
    case(iLookFLUX%scalarSurfaceRunoff); flux_data%scalarSurfaceRunoff = flux_mean%scalarSurfaceRunoff
    case(iLookFLUX%mLayerSatHydCondMP); flux_data%mLayerSatHydCondMP_m = flux_mean%mLayerSatHydCondMP_m
    case(iLookFLUX%mLayerSatHydCond); flux_data%mLayerSatHydCond_m = flux_mean%mLayerSatHydCond_m
    case(iLookFLUX%iLayerSatHydCond); flux_data%iLayerSatHydCond_m = flux_mean%iLayerSatHydCond_m
    case(iLookFLUX%mLayerHydCond); flux_data%mLayerHydCond_m = flux_mean%mLayerHydCond_m
    case(iLookFLUX%iLayerLiqFluxSoil); flux_data%iLayerLiqFluxSoil_m = flux_mean%iLayerLiqFluxSoil_m
    case(iLookFLUX%mLayerLiqFluxSoil); flux_data%mLayerLiqFluxSoil_m = flux_mean%mLayerLiqFluxSoil_m
    case(iLookFLUX%mLayerBaseflow); flux_data%mLayerBaseflow_m = flux_mean%mLayerBaseflow_m
    case(iLookFLUX%mLayerColumnInflow); flux_data%mLayerColumnInflow = flux_mean%mLayerColumnInflow
    case(iLookFLUX%mLayerColumnOutflow); flux_data%mLayerColumnOutflow_m = flux_mean%mLayerColumnOutflow_m
    case(iLookFLUX%scalarSoilBaseflow); flux_data%scalarSoilBaseflow = flux_mean%scalarSoilBaseflow
    case(iLookFLUX%scalarSoilDrainage); flux_data%scalarSoilDrainage = flux_mean%scalarSoilDrainage
    case(iLookFLUX%scalarAquiferRecharge); flux_data%scalarAquiferRecharge = flux_mean%scalarAquiferRecharge
    case(iLookFLUX%scalarAquiferTranspire); flux_data%scalarAquiferTranspire = flux_mean%scalarAquiferTranspire
    case(iLookFLUX%scalarAquiferBaseflow); flux_data%scalarAquiferBaseflow = flux_mean%scalarAquiferBaseflow
    case(iLookFLUX%scalarTotalET); flux_data%scalarTotalET = flux_mean%scalarTotalET
    case(iLookFLUX%scalarTotalRunoff); flux_data%scalarTotalRunoff = flux_mean%scalarTotalRunoff
    case(iLookFLUX%scalarNetRadiation); flux_data%scalarNetRadiation = flux_mean%scalarNetRadiation
    end select
  end if
    end do
    ! keep soil compression as an average like the fluxes, will not want to do this if nSoil can change
    diag_data%mLayerCompress_m = meanSoilCompress
   associate(scalarSoilCompress => diag_data%scalarSoilCompress, &
    mLayerDepth => prog_data%mLayerDepth, &
    nSnow => indx_data%nSnow)
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
      scalarSoilCompress(iGRU) = 0
      do iLayer=1,nSoil
        scalarSoilCompress(iGRU) = meanSoilCompress(iLayer,iGRU)*iden_water &
                                                             * mLayerDepth(nSnow(iGRU)+iLayer,iGRU )
      end do
    end do
    end associate
    deallocate(innerSoilCompress)
    deallocate(meanSoilCompress)

    ! ***********************************************************************************************************************************
    ! ---
    ! *** balance checks and summary variable saving...
    ! ---------------------

    ! save the average compression and melt pond storage in the data structures
    prog_data%scalarSfcMeltPond  = sfcMeltPond

    ! associate local variables with information in the data structures
    associate(&
      ! model forcing
      scalarSnowfall             => flux_mean%scalarSnowfall     ,&  ! computed snowfall rate (kg m-2 s-1)
      scalarRainfall             => flux_mean%scalarRainfall     ,&  ! computed rainfall rate (kg m-2 s-1)
      ! canopy fluxes
      averageThroughfallSnow     => flux_mean%scalarThroughfallSnow     ,&  ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
      averageThroughfallRain     => flux_mean%scalarThroughfallRain     ,&  ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
      averageCanopySnowUnloading => flux_mean%scalarCanopySnowUnloading     ,&  ! unloading of snow from the vegetion canopy (kg m-2 s-1)
      averageCanopyLiqDrainage   => flux_mean%scalarCanopyLiqDrainage     ,&  ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
      averageCanopySublimation   => flux_mean%scalarCanopySublimation     ,&  ! canopy sublimation/frost (kg m-2 s-1)
      averageCanopyEvaporation   => flux_mean%scalarCanopyEvaporation     ,&  ! canopy evaporation/condensation (kg m-2 s-1)
      ! snow fluxes
      averageSnowSublimation     => flux_mean%scalarSnowSublimation     ,&  ! sublimation from the snow surface (kg m-2 s-1)
      averageSnowDrainage        => flux_mean%scalarSnowDrainage     ,&  ! drainage from the bottom of the snowpack (m s-1)
      ! soil fluxes
      averageSoilInflux          => flux_mean%scalarInfiltration     ,&  ! influx of water at the top of the soil profile (m s-1)
      averageSoilDrainage        => flux_mean%scalarSoilDrainage     ,&  ! drainage from the bottom of the soil profile (m s-1)
      averageSoilBaseflow        => flux_mean%scalarSoilBaseflow     ,&  ! total baseflow from throughout the soil profile (m s-1)
      averageSoilCompress        => diag_data%scalarSoilCompress     ,&  ! soil compression (kg m-2 s-1)
      averageGroundEvaporation   => flux_mean%scalarGroundEvaporation    ,&  ! soil evaporation (kg m-2 s-1)
      averageCanopyTranspiration => flux_mean%scalarCanopyTranspiration  ,&  ! canopy transpiration (kg m-2 s-1)
      ! state variables in the vegetation canopy
      scalarCanopyWat            => prog_data%scalarCanopyWat                               ,&  ! canopy ice content (kg m-2)
      scalarCanopyIce            => prog_data%scalarCanopyIce                               ,& ! ice content of the vegetation canopy (kg m-2)
      scalarCanopyEnthTemp       => diag_data%scalarCanopyEnthTemp                          ,& ! temperature component of enthalpy of the vegetation canopy (K)
      scalarCanopyEnthalpy       => prog_data%scalarCanopyEnthalpy                          ,& ! enthalpy of the vegetation canopy (J m-3)
       ! state variables in the snow+soil domains
      scalarSWE                  => prog_data%scalarSWE                                     ,&  ! snow water equivalent (kg m-2)
      mLayerDepth                => prog_data%mLayerDepth                                   ,&  ! depth of each layer (m)
      mLayerVolFracIce           => prog_data%mLayerVolFracIce                              ,&  ! volumetric ice content in each layer (-)
      mLayerVolFracLiq           => prog_data%mLayerVolFracLiq                              ,&  ! volumetric liquid water content in each layer (-)
      scalarTotalSoilWat         => diag_data%scalarTotalSoilWat                            ,&  ! total water in the soil column (kg m-2)
      scalarTotalSoilIce         => diag_data%scalarTotalSoilIce                            ,&  ! total ice in the soil column (kg m-2)
      scalarTotalSoilLiq         => diag_data%scalarTotalSoilLiq                            ,&  ! total liquid water in the soil column (kg m-2)
      mLayerEnthTemp             => diag_data%mLayerEnthTemp                                ,& ! temperature component of enthalpy of each snow+soil layer (K)
      mLayerEnthalpy             => prog_data%mLayerEnthalpy                                ,& ! enthalpy of each snow+soil layer (J m-3)
      scalarTotalSoilEnthalpy    => diag_data%scalarTotalSoilEnthalpy                       ,& ! total enthalpy of the soil column (J m-3)
      scalarTotalSnowEnthalpy    => diag_data%scalarTotalSnowEnthalpy                       ,& ! total enthalpy of the snow column (J m-3)
      ! state variables in the aquifer
      balanceCasNrg => diag_data%balanceCasNrg, &
balanceVegNrg => diag_data%balanceVegNrg, &
balanceVegMass => diag_data%balanceVegMass, &
balanceAqMass => diag_data%balanceAqMass, &
balanceSnowNrg => diag_data%balanceSnowNrg, &
balanceSoilNrg => diag_data%balanceSoilNrg, &
balanceSnowMass => diag_data%balanceSnowMass, &
balanceSoilMass => diag_data%balanceSoilMass, &
      scalarAquiferStorage       => prog_data%scalarAquiferStorage,                          &  ! aquifer storage (m)
      nSnow => indx_data%nSnow &
      ! error tolerance
      ! absConvTol_liquid          => mpar_data_d%absConvTol_liquid                             &  ! absolute convergence tolerance for vol frac liq water (-)
      ) ! (association of local variables with information in the data structures

      ! save balance of energy and water per single layer domain
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU

      balanceCasNrg(iGRU)   = meanBalance(1,iGRU) ! W m-3
      balanceVegNrg(iGRU)   = meanBalance(2,iGRU) ! W m-3      will be realMissing if computeVegFlux is false
      balanceVegMass(iGRU)  = meanBalance(3,iGRU) ! kg m-3 s-1 will be realMissing if computeVegFlux is false
      balanceAqMass(iGRU)   = meanBalance(4,iGRU) ! kg m-2 s-1 will be realMissing if no aquifer
      balanceSnowNrg(iGRU)  = meanBalance(5,iGRU) ! W m-3      will be realMissing if no snow during data step
      balanceSoilNrg(iGRU)  = meanBalance(6,iGRU) ! W m-3       
      balanceSnowMass(iGRU) = meanBalance(7,iGRU) ! kg m-3 s-1 will be realMissing if no snow during data step
      balanceSoilMass(iGRU) = meanBalance(8,iGRU) ! kg m-3 s-1
      end do
      if(.not.bal_veg)then ! will be 0, make realMissing
        diag_data%balanceCasNrg   = realMissing
        diag_data%balanceVegNrg   = realMissing
        diag_data%balanceVegMass  = realMissing
      endif
      if (.not.bal_snow)then ! will be 0, make realMissing
        diag_data%balanceSnowNrg  = realMissing
        diag_data%balanceSnowMass = realMissing
      endif
      if (.not.bal_soil)then ! will be 0, make realMissing
        diag_data%balanceSoilNrg  = realMissing
        diag_data%balanceSoilMass = realMissing
      endif
      if (.not.bal_aq)then ! will be 0, make realMissing
        diag_data%balanceAqMass   = realMissing
      endif



      ! -----
      ! * balance checks for soil...
      ! ----------------------------

      ! compute the liquid water and ice content at the end of the time step
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        scalarTotalSoilLiq(iGRU) = 0
        scalarTotalSoilIce(iGRU) = 0
      do iLayer=1,nSoil
      scalarTotalSoilLiq(iGRU) = scalarTotalSoilLiq(iGRU) + mLayerVolFracLiq(iLayer+nSnow(iGRU),iGRU)*mLayerDepth(iLayer+nSnow(iGRU),iGRU)
      scalarTotalSoilIce(iGRU) = scalarTotalSoilIce(iGRU) + mLayerVolFracIce(iLayer+nSnow(iGRU),iGRU)*mLayerDepth(iLayer+nSnow(iGRU),iGRU)   ! NOTE: no expansion of soil, hence use iden_water
      end do

      ! get the total water in the soil (liquid plus ice) at the end of the time step (kg m-2)
      scalarTotalSoilWat(iGRU) = scalarTotalSoilLiq(iGRU) + scalarTotalSoilIce(iGRU)
    end do
      ! -----
      ! save the enthalpy or temperature component of enthalpy, and total enthalpy
      ! ----------------------------

      if(computeEnthalpy)then ! use enthTemp to conserve energy or compute energy balance  
        ! initialize the enthalpy
        scalarCanopyEnthalpy = scalarCanopyEnthTemp
        mLayerEnthalpy       = mLayerEnthTemp
        ! compute enthalpy for current values
        call enthTemp_or_enthalpy(&
                        ! input: data structures
                        .true.,                & ! intent(in):    flag to convert enthTemp to enthalpy
                        nGRU, &
                        diag_data,             & ! intent(in):    model diagnostic variables for a local HRU
                        indx_data,             & ! intent(in):    model indices
                        ! input: ice content change
                        scalarCanopyIce,       & ! intent(in):    value for canopy ice content (kg m-2)
                        mLayerVolFracIce,      & ! intent(in):    vector of volumetric ice water content (-)
                        ! input/output: enthalpy
                        scalarCanopyEnthalpy,  & ! intent(inout): enthTemp to enthalpy of the vegetation canopy (J m-3)
                        mLayerEnthalpy,        & ! intent(inout): enthTemp to enthalpy of each snow+soil layer (J m-3)
                        ! output: error control    
                        err,cmessage)            ! intent(out): error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if

      if(enthalpyStateVec)then ! enthalpy as state variable
        ! initialize the temperature component of enthalpy
        scalarCanopyEnthTemp = scalarCanopyEnthalpy
        mLayerEnthTemp       = mLayerEnthalpy          
        ! compute temperature component of enthalpy for current values       
        call enthTemp_or_enthalpy(&
                        ! input: data structures
                        .false.,               & ! intent(in):    flag to convert enthalpy to enthTemp
                        nGRU, &
                        diag_data,             & ! intent(in):    model diagnostic variables for a local HRU
                        indx_data,             & ! intent(in):    model indices
                        ! input: ice content change
                        scalarCanopyIce,       & ! intent(in):    value for canopy ice content (kg m-2)
                        mLayerVolFracIce,      & ! intent(in):    vector of volumetric ice water content (-)
                        ! input/output: enthalpy
                        scalarCanopyEnthTemp,  & ! intent(inout): enthalpy to enthTemp of the vegetation canopy (J m-3)
                        mLayerEnthTemp,        & ! intent(inout): enthalpy to enthTemp of each snow+soil layer (J m-3)
                        ! output: error control    
                        err,cmessage)            ! intent(out): error control
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif

      associate(scalarSurfaceTemp => prog_data%scalarSurfaceTemp, &
        mLayerTemp => prog_data%mLayerTemp)
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      ! save the total soil enthalpy
      scalarTotalSoilEnthalpy(iGRU) = 0
      sum_depth(iGRU) = 0
      do iLayer=1,nSoil
        scalarTotalSoilEnthalpy(iGRU) = scalarTotalSoilEnthalpy(iGRU) + mLayerEnthalpy(nSnow(iGRU)+iLayer,iGRU) * mLayerDepth(nSnow(iGRU)+iLayer,iGRU)
        sum_depth(iGRU) = sum_depth(iGRU) + mLayerDepth(nSnow(iGRU)+iLayer,iGRU)
      end do
      scalarTotalSoilEnthalpy(iGRU) = scalarTotalSoilEnthalpy(iGRU) + sum_depth(iGRU)
      if (nSnow(iGRU) > 0) scalarTotalSnowEnthalpy(iGRU) = 0
      sum_depth(iGRU) = 0
      do iLayer=1,nSnow(iGRU)
        scalarTotalSoilEnthalpy(iGRU) = scalarTotalSoilEnthalpy(iGRU) + mLayerEnthalpy(iLayer,iGRU) + mLayerDepth(iLayer,iGRU)
        sum_depth(iGRU) = sum_depth(iGRU) + mLayerDepth(iLayer,iGRU)
      end do
      if (nSnow(iGRU) > 0) scalarTotalSnowEnthalpy(iGRU) = scalarTotalSnowEnthalpy(iGRU) / sum_depth(iGRU)
      ! save the surface temperature (just to make things easier to visualize)
      scalarSurfaceTemp(iGRU) = mLayerTemp(1,iGRU)

    end do
    end associate

    end associate ! end association of local variables with information in the data structures

  end associate canopy  ! end association to canopy parameters

  ! overwrite flux data with timestep-average value for all flux_mean vars, hard-coded to not happen
  if(.not.backwardsCompatibility)then
    do iVar=1,size(flux_meta)
      if (flux_meta(iVar)%vartype==iLookVarType%scalarv) then
        select case(iVar)
        case(iLookFLUX%scalarCanairNetNrgFlux); flux_data%scalarCanairNetNrgFlux = flux_mean%scalarCanairNetNrgFlux
        case(iLookFLUX%scalarCanopyNetNrgFlux); flux_data%scalarCanopyNetNrgFlux = flux_mean%scalarCanopyNetNrgFlux
        case(iLookFLUX%scalarGroundNetNrgFlux); flux_data%scalarGroundNetNrgFlux = flux_mean%scalarGroundNetNrgFlux
        case(iLookFLUX%scalarCanopyNetLiqFlux); flux_data%scalarCanopyNetLiqFlux = flux_mean%scalarCanopyNetLiqFlux
        case(iLookFLUX%scalarRainfall); flux_data%scalarRainfall = flux_mean%scalarRainfall
        case(iLookFLUX%scalarSnowfall); flux_data%scalarSnowfall = flux_mean%scalarSnowfall
        case(iLookFLUX%scalarCanopySunlitPAR); flux_data%scalarCanopySunlitPAR = flux_mean%scalarCanopySunlitPAR
        case(iLookFLUX%scalarCanopyShadedPAR); flux_data%scalarCanopyShadedPAR = flux_mean%scalarCanopyShadedPAR
        case(iLookFLUX%scalarBelowCanopySolar); flux_data%scalarBelowCanopySolar = flux_mean%scalarBelowCanopySolar
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
        case(iLookFLUX%scalarGroundResistance); flux_data%scalarGroundResistance = flux_mean%scalarGroundResistance
        case(iLookFLUX%scalarCanopyResistance); flux_data%scalarCanopyResistance = flux_mean%scalarCanopyResistance
        case(iLookFLUX%scalarLeafResistance); flux_data%scalarLeafResistance = flux_mean%scalarLeafResistance
        case(iLookFLUX%scalarSoilResistance); flux_data%scalarSoilResistance = flux_mean%scalarSoilResistance
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
        case(iLookFLUX%scalarStomResistSunlit); flux_data%scalarStomResistSunlit = flux_mean%scalarStomResistSunlit
        case(iLookFLUX%scalarStomResistShaded); flux_data%scalarStomResistShaded = flux_mean%scalarStomResistShaded
        case(iLookFLUX%scalarPhotosynthesisSunlit); flux_data%scalarPhotosynthesisSunlit = flux_mean%scalarPhotosynthesisSunlit
        case(iLookFLUX%scalarPhotosynthesisShaded); flux_data%scalarPhotosynthesisShaded = flux_mean%scalarPhotosynthesisShaded
        case(iLookFLUX%scalarCanopyTranspiration); flux_data%scalarCanopyTranspiration = flux_mean%scalarCanopyTranspiration
        case(iLookFLUX%scalarCanopyEvaporation); flux_data%scalarCanopyEvaporation = flux_mean%scalarCanopyEvaporation
        case(iLookFLUX%scalarGroundEvaporation); flux_data%scalarGroundEvaporation = flux_mean%scalarGroundEvaporation
        case(iLookFLUX%scalarThroughfallSnow); flux_data%scalarThroughfallSnow = flux_mean%scalarThroughfallSnow
        case(iLookFLUX%scalarThroughfallRain); flux_data%scalarThroughfallRain = flux_mean%scalarThroughfallRain
        case(iLookFLUX%scalarCanopySnowUnloading); flux_data%scalarCanopySnowUnloading = flux_mean%scalarCanopySnowUnloading
        case(iLookFLUX%scalarCanopyLiqDrainage); flux_data%scalarCanopyLiqDrainage = flux_mean%scalarCanopyLiqDrainage
        case(iLookFLUX%scalarCanopyMeltFreeze); flux_data%scalarCanopyMeltFreeze = flux_mean%scalarCanopyMeltFreeze
        case(iLookFLUX%scalarSnowDrainage); flux_data%scalarSnowDrainage = flux_mean%scalarSnowDrainage
        case(iLookFLUX%scalarRainPlusMelt); flux_data%scalarRainPlusMelt = flux_mean%scalarRainPlusMelt
        case(iLookFLUX%scalarMaxInfilRate); flux_data%scalarMaxInfilRate = flux_mean%scalarMaxInfilRate
        case(iLookFLUX%scalarInfiltration); flux_data%scalarInfiltration = flux_mean%scalarInfiltration
        case(iLookFLUX%scalarExfiltration); flux_data%scalarExfiltration = flux_mean%scalarExfiltration
        case(iLookFLUX%scalarSurfaceRunoff); flux_data%scalarSurfaceRunoff = flux_mean%scalarSurfaceRunoff
        case(iLookFLUX%scalarSoilBaseflow); flux_data%scalarSoilBaseflow = flux_mean%scalarSoilBaseflow
        case(iLookFLUX%scalarSoilDrainage); flux_data%scalarSoilDrainage = flux_mean%scalarSoilDrainage
        case(iLookFLUX%scalarAquiferRecharge); flux_data%scalarAquiferRecharge = flux_mean%scalarAquiferRecharge
        case(iLookFLUX%scalarAquiferTranspire); flux_data%scalarAquiferTranspire = flux_mean%scalarAquiferTranspire
        case(iLookFLUX%scalarAquiferBaseflow); flux_data%scalarAquiferBaseflow = flux_mean%scalarAquiferBaseflow
        case(iLookFLUX%scalarTotalET); flux_data%scalarTotalET = flux_mean%scalarTotalET
        case(iLookFLUX%scalarTotalRunoff); flux_data%scalarTotalRunoff = flux_mean%scalarTotalRunoff
        case(iLookFLUX%scalarNetRadiation); flux_data%scalarNetRadiation = flux_mean%scalarNetRadiation
        end select
      end if
    end do
  end if

  if(nsub>50000)then
    write(message,'(a,i0)') trim(cmessage)//'number of sub-steps > 50000 for HRU ', hruId
    err=20; return
  end if

  ! get the end time
  CALL system_clock(i_end)
  elapsed_time = REAL(i_end - i_start) / REAL(count_rate)

  ! get the elapsed time
  diag_data%wallClockTime = elapsed_time
      ! call finalize_device_indx_data(indx_data,indx_data2)
      ! call deallocate_device_indx_data(indx_data)
                call deallocate_device_flux_data(flux_inner)
                                call deallocate_device_flux_data(flux_mean)

        call deallocate_device_indx_data(indx_temp)
        call deallocate_device_diag_data(diag_temp)
        call deallocate_device_prog_data(prog_temp)

contains

 subroutine initialize_coupled_em
  ! *** Initialize steps for coupled_em subroutine ***
  ! Notes: - created to ensure certain variables are initialized prior to use in calculations
  !        - based on warnings from the SUMMA debug build (e.g., -Wall flag)
  !        - additional initial operations may be added here in the future

  ! initialize variables
  innerEffRainfall=0._rkind       ! inner step average effective rainfall into snow (kg m-2 s-1) 
  sumCanopySublimation=0._rkind   ! sum of sublimation from the vegetation canopy (kg m-2 s-1) over substep
  sumLatHeatCanopyEvap=0._rkind   ! sum of latent heat flux for evaporation from the canopy to the canopy air space (W m-2) over substep
  sumSenHeatCanopy=0._rkind       ! sum of sensible heat flux from the canopy to the canopy air space (W m-2) over substep
  sumSnowSublimation=0._rkind     ! sum of sublimation from the snow surface (kg m-2 s-1) over substep
  sumStepSize=0._rkind            ! sum solution step for the data step
  innerBalance = 0._rkind         ! mean total balance array

  ! get initial value of nLayers
        associate(nSnow => indx_data%nSnow, layerType => indx_data%layerType)
            !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    nSnow(iGRU) = 0
    do iLayer=1,size(layerType,1)
      if (layerType(iLayer,iGRU) == iname_snow) nSnow(iGRU) = nSnow(iGRU) + 1
    end do
  !   if (createLayer(iGRU)) then
  end do
  end associate
  nSnow = maxval(indx_data%nSnow)
  nSoil = indx_data%nSoil
  nLayers = nSnow + nSoil

  ! allocate and initialize using the initial value of nLayers
  allocate(innerBalanceLayerMass(nLayers,nGRU)); innerBalanceLayerMass = 0._rkind ! mean total balance of mass in layers
  allocate(innerBalanceLayerNrg(nLayers,nGRU));  innerBalanceLayerNrg = 0._rkind ! mean total balance of energy in layers
  allocate(mLayerVolFracIceInit(nLayers,nGRU));  mLayerVolFracIceInit = prog_data%mLayerVolFracIce ! volume fraction of water ice

   ! initialize the numerix tracking variables
  indx_data%numberFluxCalc        = 0  ! number of flux calculations                     (-)
  ! indx_data%numberStateSplit      = 0  ! number of state splitting solutions             (-)
  ! indx_data%numberDomainSplitNrg  = 0  ! number of domain splitting solutions for energy (-)
  ! indx_data%numberDomainSplitMass = 0  ! number of domain splitting solutions for mass   (-)
  ! indx_data%numberScalarSolutions = 0  ! number of scalar solutions                      (-)

  ! initialize surface melt pond
  sfcMeltPond       = 0._rkind  ! change in storage associated with the surface melt pond (kg m-2)

  ! initialize average over data_step (averaged over substep in varSubStep)
  meanCanopySublimation = 0._rkind ! mean canopy sublimation
  meanLatHeatCanopyEvap = 0._rkind ! mean latent heat flux for evaporation from the canopy
  meanSenHeatCanopy     = 0._rkind ! mean sensible heat flux from the canopy
  effRainfall           = 0._rkind ! mean total effective rainfall over snow

  diag_data%meanStepSize = 0._rkind ! mean step size over data_step

  ! Need mean soil compression for balance checks but it is not in flux structure so handle differently 
  !  This will be a problem if nSoil changes (currently not possible)-- then might need to not keep the average
  allocate(meanSoilCompress(nSoil,nGRU))
  allocate(innerSoilCompress(nSOil,nGRU))
  meanSoilCompress = 0._rkind ! mean total soil compression

  ! initialize the balance checks
  meanBalance = 0._rkind

 end subroutine initialize_coupled_em

end subroutine coupled_em


! *********************************************************************************************************
! private subroutine implctMelt: compute melt of the "snow without a layer"
! *********************************************************************************************************
attributes(device) subroutine implctMelt(&
                      ! input/output: integrated snowpack properties
                      scalarSWE,         & ! intent(inout): snow water equivalent (kg m-2)
                      scalarSnowDepth,   & ! intent(inout): snow depth (m)
                      scalarSfcMeltPond, & ! intent(inout): surface melt pond (kg m-2)
                      ! input/output: properties of the upper-most soil layer
                      soilTemp,          & ! intent(inout): surface layer temperature (K)
                      soilDepth,         & ! intent(inout): surface layer depth (m)
                      soilHeatcap       & ! intent(inout): surface layer volumetric heat capacity (J m-3 K-1)
                      ! output: error control
                      ) ! intent(out): error control
  implicit none
  ! input/output: integrated snowpack properties
  real(rkind),intent(inout)    :: scalarSWE          ! snow water equivalent (kg m-2)
  real(rkind),intent(inout)    :: scalarSnowDepth    ! snow depth (m)
  real(rkind),intent(inout)    :: scalarSfcMeltPond  ! surface melt pond (kg m-2)
  ! input/output: properties of the upper-most soil layer
  real(rkind),intent(inout)    :: soilTemp           ! surface layer temperature (K)
  real(rkind),intent(inout)    :: soilDepth          ! surface layer depth (m)
  real(rkind),intent(inout)    :: soilHeatcap        ! surface layer volumetric heat capacity (J m-3 K-1)
  ! output: error control
  ! integer(i4b),intent(out)  :: err                ! error code
  ! character(*),intent(out)  :: message            ! error message
  ! local variables
  real(rkind)                  :: nrgRequired        ! energy required to melt all the snow (J m-2)
  real(rkind)                  :: nrgAvailable       ! energy available to melt the snow (J m-2)
  real(rkind)                  :: snwDensity         ! snow density (kg m-3)
  ! initialize error control
  ! err=0; message='implctMelt/'

  if(scalarSWE > 0._rkind)then
    ! only melt if temperature of the top soil layer is greater than Tfreeze
    if(soilTemp > Tfreeze)then
      ! compute the energy required to melt all the snow (J m-2)
      nrgRequired     = scalarSWE*LH_fus
      ! compute the energy available to melt the snow (J m-2)
      nrgAvailable    = soilHeatcap*(soilTemp - Tfreeze)*soilDepth
      ! compute the snow density (not saved)
      snwDensity      = scalarSWE/scalarSnowDepth
      ! compute the amount of melt, and update SWE (kg m-2)
      if(nrgAvailable > nrgRequired)then
        scalarSfcMeltPond  = scalarSWE
        scalarSWE          = 0._rkind
      else
        scalarSfcMeltPond  = nrgAvailable/LH_fus
        scalarSWE          = scalarSWE - scalarSfcMeltPond
      end if
      ! update depth
      scalarSnowDepth = scalarSWE/snwDensity
      ! update temperature of the top soil layer (K)
      soilTemp =  soilTemp - (LH_fus*scalarSfcMeltPond/soilDepth)/soilHeatcap
    else  ! melt is zero if the temperature of the top soil layer is less than Tfreeze
      scalarSfcMeltPond = 0._rkind  ! kg m-2
    end if ! (if the temperature of the top soil layer is greater than Tfreeze)
  else  ! melt is zero if the "snow without a layer" does not exist
    scalarSfcMeltPond = 0._rkind  ! kg m-2
  end if ! (if the "snow without a layer" exists)

end subroutine implctMelt


attributes(global) subroutine enthtemp_soil1(nGRU, nSnow, &
  enthalpyStateVec,computeEnthalpy,use_lookup, &
  scalarSWE,&
  vGn_alpha, vGn_n, theta_sat, theta_res, vGn_m, &
  temperature,psiLiq_int,deriv2, &
  mLayerTemp, mLayerEnthTemp, mLayerEnthalpy, &
  mLayerMatricHead, mLayerVolFracIce, soil_dens_intr)
  USE enthalpyTemp_module,only:T2enthTemp_soil      ! convert temperature to enthalpy for soil

  implicit none
  integer(i4b),value :: nGRU
  logical(lgt),value :: enthalpyStateVec,computeEnthalpy,use_lookup
  integer(i4b) :: nSnow(:)
  real(rkind) :: scalarSWE(:)
  real(rkind) :: vGn_alpha(:), vGn_n(:), theta_sat(:), theta_res(:)
  real(rkind) :: vGn_m(:,:)
  real(rkind) :: temperature(:,:,:), psiLiq_int(:,:,:),deriv2(:,:,:)
  real(rkind) :: mLayerTemp(:,:), mLayerEnthTemp(:,:), mLayerEnthalpy(:,:)
  real(rkind) :: mLayerMatricHead(:,:)
  real(rkind) :: mLayerVolFracIce(:,:)
  real(rkind) :: soil_dens_intr
  integer(i4b) :: iGRU
  
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return
  if( (enthalpyStateVec .or. computeEnthalpy) .and. nSnow(iGRU)==0 .and. scalarSWE(iGRU)>0._rkind )then
    call T2enthTemp_soil(&
              use_lookup,                                               & ! intent(in):  flag to use the lookup table for soil enthalpy
              soil_dens_intr,                                           & ! intent(in):  intrinsic soil density (kg m-3)
              vGn_alpha(1),vGn_n(1),theta_sat(1),theta_res(1),vGn_m(1,iGRU), & ! intent(in):  van Genutchen soil parameters
              ! 1_i4b,                                                    & ! intent(in):  index of the control volume within the domain
              ! lookup_data,                                              & ! intent(in):  lookup table data structure
              temperature(:,1,iGRU),psiLiq_int(:,1,iGRU),deriv2(:,1,iGRU), &
              realMissing,                                              & ! intent(in):  lower value of integral (not computed)
              mLayerTemp(nSnow(iGRU)+1,iGRU),         & ! intent(in):  surface layer temperature (K)
              mLayerMatricHead(1,iGRU),                                      & ! intent(in):  surface layer matric head (m)
              mLayerEnthTemp(nSnow(iGRU)+1,iGRU))       ! intent(out): temperature component of enthalpy soil layer (J m-3)
    mLayerEnthalpy(nSnow(iGRU)+1,iGRU) = mLayerEnthTemp(nSnow(iGRU)+1,iGRU) - iden_water * LH_fus * mLayerVolFracIce(nSnow(iGRU)+1,iGRU)
  end if

end subroutine

attributes(global) subroutine enthtemp_soilnSoil(nGRU, nSnow,nSoil, &
  use_lookup, &
  ! scalarSWE,&
  vGn_alpha, vGn_n, theta_sat, theta_res, vGn_m, &
  temperature,psiLiq_int,deriv2, &
  mLayerTemp, mLayerEnthTemp, mLayerEnthalpy, &
  mLayerMatricHead, mLayerVolFracIce, mLayerVolFracWat, mLayerVolFracLiq, soil_dens_intr)
  USE enthalpyTemp_module,only:T2enthTemp_soil      ! convert temperature to enthalpy for soil

  implicit none
  integer(i4b),value :: nGRU
  logical(lgt),value :: use_lookup
  integer(i4b) :: nSnow(:)
  integer(i4b),value :: nSoil
  real(rkind) :: vGn_alpha(:), vGn_n(:), theta_sat(:), theta_res(:)
  real(rkind) :: vGn_m(:,:)
  real(rkind) :: temperature(:,:,:), psiLiq_int(:,:,:),deriv2(:,:,:)
  real(rkind) :: mLayerTemp(:,:), mLayerEnthTemp(:,:), mLayerEnthalpy(:,:)
  real(rkind) :: mLayerMatricHead(:,:)
  real(rkind) :: mLayerVolFracWat(:,:), mLayerVolFracLiq(:,:), mLayerVolFracIce(:,:)
  real(rkind) :: soil_dens_intr(:)
  integer(i4b) :: iGRU, iLayer, iSoil
  
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return
  do iSoil=1,nSoil
    iLayer = iSoil + nSnow(iGRU)
              mLayerVolFracWat(iLayer,iGRU) = mLayerVolFracLiq(iLayer,iGRU) + mLayerVolFracIce(iLayer,iGRU)
              ! compute enthalpy for soil layers
              call T2enthTemp_soil(&
                             use_lookup,                   & ! intent(in):  flag to use the lookup table for soil enthalpy
                             soil_dens_intr(iSoil),        & ! intent(in):  intrinsic soil density (kg m-3)
                             vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m(iSoil,iGRU), & ! intent(in): soil parameters
                             temperature(:,iSoil,iGRU),psiLiq_int(:,iSoil,iGRU),deriv2(:,iSoil,iGRU), &
                            !  iSoil,                        & ! intent(in):  index of the control volume within the domain
                            !  lookup_data,                  & ! intent(in):  lookup table data structure
                             realMissing,                  & ! intent(in):  lower value of integral (not computed)
                             mLayerTemp(iLayer,iGRU),           & ! intent(in):  layer temperature (K)
                             mLayerMatricHead(iSoil,iGRU),      & ! intent(in):  matric head (m)
                             mLayerEnthTemp(iLayer,iGRU))         ! intent(out): temperature component of enthalpy soil layer (J m-3)
              mLayerEnthalpy(iLayer,iGRU) = mLayerEnthTemp(iLayer,iGRU) - iden_water * LH_fus * mLayerVolFracIce(iLayer,iGRU)

  end do
end subroutine

subroutine update_flux_inner_array(inner, data, nGRU, dt_wght)
  real(rkind),device :: inner(:), data(:)
  real(rkind) :: dt_wght
  integer(i4b) :: nGRU
  integer(i4b) :: iGRU

  associate(inner_ => inner, data_ => data)
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
      inner_(iGRU) = inner_(iGRU) + data_(iGRU)*dt_wght
    end do
  end associate
end subroutine

subroutine update_flux_inner(flux_inner, flux_data, nGRU, dt_wght)
  use device_data_types
  type(flux_data_device) :: flux_inner, flux_data
  real(rkind) :: dt_wght
  integer(i4b) :: nGRU

  call update_flux_inner_array(flux_inner%scalarCanairNetNrgFlux, flux_data%scalarCanairNetNrgFlux, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopyNetNrgFlux, flux_data%scalarCanopyNetNrgFlux, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarGroundNetNrgFlux, flux_data%scalarGroundNetNrgFlux, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopyNetLiqFlux, flux_data%scalarCanopyNetLiqFlux, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarRainfall, flux_data%scalarRainfall, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarSnowfall, flux_data%scalarSnowfall, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%spectralIncomingDirect, flux_data%spectralIncomingDirect, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%spectralIncomingDiffuse, flux_data%spectralIncomingDiffuse, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopySunlitPAR, flux_data%scalarCanopySunlitPAR, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopyShadedPAR, flux_data%scalarCanopyShadedPAR, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%spectralBelowCanopyDirect, flux_data%spectralBelowCanopyDirect, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%spectralBelowCanopyDiffuse, flux_data%spectralBelowCanopyDiffuse, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarBelowCanopySolar, flux_data%scalarBelowCanopySolar, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopyAbsorbedSolar, flux_data%scalarCanopyAbsorbedSolar, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarGroundAbsorbedSolar, flux_data%scalarGroundAbsorbedSolar, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWRadCanopy, flux_data%scalarLWRadCanopy, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWRadGround, flux_data%scalarLWRadGround, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWRadUbound2Canopy, flux_data%scalarLWRadUbound2Canopy, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWRadUbound2Ground, flux_data%scalarLWRadUbound2Ground, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWRadUbound2Ubound, flux_data%scalarLWRadUbound2Ubound, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWRadCanopy2Ubound, flux_data%scalarLWRadCanopy2Ubound, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWRadCanopy2Ground, flux_data%scalarLWRadCanopy2Ground, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWRadCanopy2Canopy, flux_data%scalarLWRadCanopy2Canopy, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWRadGround2Ubound, flux_data%scalarLWRadGround2Ubound, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWRadGround2Canopy, flux_data%scalarLWRadGround2Canopy, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWNetCanopy, flux_data%scalarLWNetCanopy, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWNetGround, flux_data%scalarLWNetGround, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLWNetUbound, flux_data%scalarLWNetUbound, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarEddyDiffusCanopyTop, flux_data%scalarEddyDiffusCanopyTop, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarFrictionVelocity, flux_data%scalarFrictionVelocity, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarWindspdCanopyTop, flux_data%scalarWindspdCanopyTop, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarWindspdCanopyBottom, flux_data%scalarWindspdCanopyBottom, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarGroundResistance, flux_data%scalarGroundResistance, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopyResistance, flux_data%scalarCanopyResistance, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLeafResistance, flux_data%scalarLeafResistance, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarSoilResistance, flux_data%scalarSoilResistance, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarSenHeatTotal, flux_data%scalarSenHeatTotal, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarSenHeatCanopy, flux_data%scalarSenHeatCanopy, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarSenHeatGround, flux_data%scalarSenHeatGround, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLatHeatTotal, flux_data%scalarLatHeatTotal, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLatHeatCanopyEvap, flux_data%scalarLatHeatCanopyEvap, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLatHeatCanopyTrans, flux_data%scalarLatHeatCanopyTrans, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarLatHeatGround, flux_data%scalarLatHeatGround, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopyAdvectiveHeatFlux, flux_data%scalarCanopyAdvectiveHeatFlux, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarGroundAdvectiveHeatFlux, flux_data%scalarGroundAdvectiveHeatFlux, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopySublimation, flux_data%scalarCanopySublimation, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarSnowSublimation, flux_data%scalarSnowSublimation, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarStomResistSunlit, flux_data%scalarStomResistSunlit, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarStomResistShaded, flux_data%scalarStomResistShaded, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarPhotosynthesisSunlit, flux_data%scalarPhotosynthesisSunlit, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarPhotosynthesisShaded, flux_data%scalarPhotosynthesisShaded, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopyTranspiration, flux_data%scalarCanopyTranspiration, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopyEvaporation, flux_data%scalarCanopyEvaporation, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarGroundEvaporation, flux_data%scalarGroundEvaporation, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%mLayerTranspire_m, flux_data%mLayerTranspire_m, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarThroughfallSnow, flux_data%scalarThroughfallSnow, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarThroughfallRain, flux_data%scalarThroughfallRain, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopySnowUnloading, flux_data%scalarCanopySnowUnloading, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopyLiqDrainage, flux_data%scalarCanopyLiqDrainage, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarCanopyMeltFreeze, flux_data%scalarCanopyMeltFreeze, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%iLayerConductiveFlux_m, flux_data%iLayerConductiveFlux_m, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%iLayerAdvectiveFlux_m, flux_data%iLayerAdvectiveFlux_m, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%iLayerNrgFlux_m, flux_data%iLayerNrgFlux_m, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%mLayerNrgFlux_m, flux_data%mLayerNrgFlux_m, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarSnowDrainage, flux_data%scalarSnowDrainage, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%iLayerLiqFluxSnow_m, flux_data%iLayerLiqFluxSnow_m, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%mLayerLiqFluxSnow_m, flux_data%mLayerLiqFluxSnow_m, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarRainPlusMelt, flux_data%scalarRainPlusMelt, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarMaxInfilRate, flux_data%scalarMaxInfilRate, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarInfiltration, flux_data%scalarInfiltration, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarExfiltration, flux_data%scalarExfiltration, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarSurfaceRunoff, flux_data%scalarSurfaceRunoff, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%mLayerSatHydCondMP_m, flux_data%mLayerSatHydCondMP_m, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%mLayerSatHydCond_m, flux_data%mLayerSatHydCond_m, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%iLayerSatHydCond_m, flux_data%iLayerSatHydCond_m, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%mLayerHydCond_m, flux_data%mLayerHydCond_m, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%iLayerLiqFluxSoil_m, flux_data%iLayerLiqFluxSoil_m, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%mLayerLiqFluxSoil_m, flux_data%mLayerLiqFluxSoil_m, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%mLayerBaseflow_m, flux_data%mLayerBaseflow_m, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%mLayerColumnInflow, flux_data%mLayerColumnInflow, nGRU, dt_wght)
  ! call update_flux_inner_array(flux_inner%mLayerColumnOutflow_m, flux_data%mLayerColumnOutflow_m, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarSoilBaseflow, flux_data%scalarSoilBaseflow, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarSoilDrainage, flux_data%scalarSoilDrainage, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarAquiferRecharge, flux_data%scalarAquiferRecharge, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarAquiferTranspire, flux_data%scalarAquiferTranspire, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarAquiferBaseflow, flux_data%scalarAquiferBaseflow, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarTotalET, flux_data%scalarTotalET, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarTotalRunoff, flux_data%scalarTotalRunoff, nGRU, dt_wght)
  call update_flux_inner_array(flux_inner%scalarNetRadiation, flux_data%scalarNetRadiation, nGRU, dt_wght)

end subroutine

end module coupled_em_module
