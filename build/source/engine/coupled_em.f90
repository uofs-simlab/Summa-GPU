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
USE globalData,only: verySmall ! a very small number used as an additive constant to check if substantial difference among real numbers

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
USE var_lookup,only:childFLUX_MEAN

! metadata
USE globalData,only:flux_meta              ! metadata on the model fluxes
USE globalData,only:indx_meta              ! metadata on the model index variables
USE globalData,only:diag_meta              ! metadata on the model diagnostic variables
USE globalData,only:prog_meta              ! metadata on the model prognostic variables
USE globalData,only:averageFlux_meta       ! metadata on the timestep-average model flux structure

! global data
USE globalData,only:data_step              ! time step of forcing data (s)
USE globalData,only:model_decisions        ! model decision structure
USE globalData,only:globalPrintFlag        ! the global print flag
USE globalData,only:realMissing            ! missing real number


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
contains


! ************************************************************************************************
! public subroutine coupled_em: run the coupled energy-mass model for one timestep
! ************************************************************************************************
subroutine coupled_em(&
  nGRU, &
                      ! model control
                      hruId,             & ! intent(in):    hruId
                      dt_init,           & ! intent(inout): used to initialize the size of the sub-step
                      dt_init_factor,    & ! intent(in):    Used to adjust the length of the timestep in the event of a failure
                      computeVegFlux,    & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                      fracJulDay,        & ! intent(in):    fractional julian days since the start of year
                      yearLength,        & ! intent(in):    number of days in the current year
                      ! data structures (input)
                      type_data,         & ! intent(in):    local classification of soil veg etc. for each HRU
                      attr_data,         & ! intent(in):    local attributes for each HRU
                      forc_data,         & ! intent(in):    model forcing data
                      mpar_data,         & ! intent(in):    model parameters
                      bvar_data,         & ! intent(in):    basin-average variables
                      lookup_data,       & ! intent(in):    lookup tables
                      decisions, veg_param, &
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
  USE vegPhenlgy_module,only:vegPhenlgy             ! compute vegetation phenology
  USE vegNrgFlux_module,only:wettedFrac             ! compute wetted fraction of the canopy (used in sw radiation fluxes)
  USE snowAlbedo_module,only:snowAlbedo             ! compute snow albedo
  USE vegSWavRad_module,only:vegSWavRad             ! compute canopy sw radiation fluxes
  USE canopySnow_module,only:canopySnow             ! compute interception and unloading of snow from the vegetation canopy
  USE volicePack_module,only:newsnwfall             ! compute change in the top snow layer due to throughfall and unloading
  USE volicePack_module,only:volicePack_device             ! merge and sub-divide snow layers, if necessary
  USE diagn_evar_module,only:diagn_evar             ! compute diagnostic energy variables -- thermal conductivity and heat capacity
  ! the model solver
  USE indexState_module,only:indexState             ! define indices for all model state variables and layers
  USE opSplittin_module,only:opSplittin             ! solve the system of thermodynamic and hydrology equations for a given substep
  USE time_utils_module,only:elapsedSec             ! calculate the elapsed time
  ! additional subroutines
  USE tempAdjust_module,only:tempAdjust             ! adjust snow temperature associated with new snowfall
  USE var_derive_module,only:calcHeight             ! module to calculate height at layer interfaces and layer mid-point
  USE computSnowDepth_module,only:computSnowDepth   ! compute snow depth
  USE enthalpyTemp_module,only:T2enthTemp_veg       ! convert temperature to enthalpy for vegetation
  USE enthalpyTemp_module,only:T2enthTemp_snow      ! convert temperature to enthalpy for snow
  USE enthalpyTemp_module,only:T2enthTemp_soil      ! convert temperature to enthalpy for soil
  USE enthalpyTemp_module,only:enthTemp_or_enthalpy ! add phase change terms to delta temperature component of enthalpy or vice versa
  use device_data_types
  use initialize_device
  use globalData,only:flux2state_orig_d,maxSnowLayers,veryBig
  use enthalpyTemp_module,only:h_lookup,t_lookup
   USE var_lookup,only:iLookVarType                 ! look up structure for variable typed
   use globalData,only:minExpLogHgtFac,urbanVegCategory

  implicit none

  integer(i4b),intent(in) :: nGRU
  integer(i8b),intent(in)              :: hruId                  ! hruId
  real(rkind),intent(inout)            :: dt_init                ! used to initialize the size of the sub-step
  integer(i4b),intent(in)              :: dt_init_factor         ! Used to adjust the length of the timestep in the event of a failure
  logical(lgt),intent(inout),device           :: computeVegFlux(nGRU)         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  ! data structures (input)
  type(type_data_device),intent(in)               :: type_data              ! type of vegetation and soil
  type(attr_data_device),intent(in)               :: attr_data              ! spatial attributes
  type(forc_data_device),intent(in)               :: forc_data              ! model forcing data
  type(mpar_data_device),intent(in)         :: mpar_data              ! model parameters
  type(bvar_data_device),intent(in)         :: bvar_data              ! basin-average model variables
  type(zLookup_device),intent(in)             :: lookup_data            ! lookup tables
  ! data structures (input-output)
  type(indx_data_device),intent(inout)      :: indx_data              ! state vector geometry
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
  logical(lgt),device                         :: computeVegFluxOld(nGRU)      ! flag to indicate if we are computing fluxes over vegetation on the previous sub step
  logical(lgt),device                         :: includeAquifer         ! flag to denote that an aquifer is included
  logical(lgt),device                         :: modifiedLayers(nGRU)         ! flag to denote that snow layers were modified
  logical(lgt),device                         :: modifiedVegState(nGRU)       ! flag to denote that vegetation states were modified
  ! integer(i4b)                         :: nLayersRoots           ! number of soil layers that contain roots
  real(rkind),device                          :: exposedVAI(nGRU)             ! exposed vegetation area index
  real(rkind)                          :: dCanopyWetFraction_dWat ! derivative in wetted fraction w.r.t. canopy total water (kg-1 m2)
  real(rkind)                          :: dCanopyWetFraction_dT   ! derivative in wetted fraction w.r.t. canopy temperature (K-1)
  real(rkind),parameter                :: varNotUsed1=-9999._rkind ! variables used to calculate derivatives (not needed here)
  real(rkind),parameter                :: varNotUsed2=-9999._rkind ! variables used to calculate derivatives (not needed here)
  integer(i4b)                         :: iLayer                 ! index of model layers
  real(rkind),device                          :: superflousSub(nGRU)          ! superflous sublimation (kg m-2 s-1)
  real(rkind),device                          :: superflousNrg(nGRU)          ! superflous energy that cannot be used for sublimation (W m-2 [J m-2 s-1])
  integer(i4b)                         :: ixSolution             ! solution method used by opSplittin
  logical(lgt)                         :: firstSubStep           ! flag to denote if the first time step
  logical(lgt)                         :: stepFailure            ! flag to denote the need to reduce length of the coupled step and try again
  logical(lgt)                         :: tooMuchMelt            ! flag to denote that there was too much melt in a given time step
  logical(lgt)                         :: tooMuchSublim          ! flag to denote that there was too much sublimation in a given time step
  logical(lgt)                         :: doLayerMerge           ! flag to denote the need to merge snow layers
  logical(lgt),parameter               :: backwardsCompatibility=.false.  ! flag to denote a desire to ensure backwards compatibility with previous branches for end of time step flux only 
  logical(lgt)                         :: checkMassBalance_ds    ! flag to check the mass balance over the data step
  type(indx_data_device)                    :: indx_temp              ! temporary model index variables saved only on outer loop
  type(indx_data_device)                    :: indx_temp0             ! temporary model index variables saved every time
  type(prog_data_device)                    :: prog_temp              ! temporary model prognostic variables
  type(diag_data_device)                    :: diag_temp              ! temporary model diagnostic variables
  real(rkind),device,allocatable              :: mLayerVolFracIceInit(:,:)! initial vector for volumetric fraction of ice (-)
  ! check SWE
  real(rkind),device                          :: oldSWE(nGRU)                 ! SWE at the start of the substep
  real(rkind),device                          :: delSWE(nGRU)                 ! change in SWE over the subtep
  real(rkind),device                          :: innerEffRainfall(nGRU)       ! inner step average effective rainfall into snow (kg m-2 s-1)
  real(rkind),device                          :: effRainfall(nGRU)            ! timestep-average effective rainfall into snow (kg m-2 s-1)
  real(rkind),device                          :: effSnowfall(nGRU)            ! effective snowfall (kg m-2 s-1)
  real(rkind),device                          :: sfcMeltPond(nGRU)            ! surface melt pond (kg m-2)
  real(rkind),device                          :: massBalance(nGRU)            ! mass balance error (kg m-2)
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
  logical(lgt),device                         :: bal_veg(nGRU)                ! flag to denote if computed a vegetation balance
  logical(lgt),device                         :: bal_snow(nGRU)               ! flag to denote if computed a snow balance
  logical(lgt),device                         :: bal_soil(nGRU)               ! flag to denote if computed a soil balance
  logical(lgt),device                         :: bal_aq(nGRU)                 ! flag to denote if computed an aquifer balance
  integer(i4b)                         :: iVar                   ! loop through model variables
  real(rkind),device                          :: balanceSoilCompress(nGRU)    ! total soil compression (kg m-2)
  real(rkind),device                          :: scalarCanopyWatBalError(nGRU)! water balance error for the vegetation canopy (kg m-2)
  real(rkind),device                          :: scalarSoilWatBalError(nGRU)  ! water balance error (kg m-2)
  real(rkind),device                          :: scalarInitCanopyLiq(nGRU)    ! initial liquid water on the vegetation canopy (kg m-2)
  real(rkind),device                          :: scalarInitCanopyIce(nGRU)    ! initial ice          on the vegetation canopy (kg m-2)
  real(rkind),device                          :: balanceCanopyWater0(nGRU)    ! total water stored in the vegetation canopy at the start of the step (kg m-2)
  real(rkind),device                          :: balanceSoilWater0(nGRU)      ! total soil storage at the start of the step (kg m-2)
  real(rkind),device                          :: balanceSoilInflux(nGRU)      ! input to the soil zone
  real(rkind),device                          :: balanceSoilBaseflow(nGRU)    ! output from the soil zone
  real(rkind),device                          :: balanceSoilDrainage(nGRU)    ! output from the soil zone
  real(rkind),device                          :: balanceSoilET(nGRU)          ! output from the soil zone
  real(rkind),device                          :: balanceAquifer0(nGRU)        ! total aquifer storage at the start of the step (kg m-2)
  real(rkind),device                          :: balanceAquifer1(nGRU)        ! total aquifer storage at the end of the step (kg m-2)
  real(rkind),device                          :: innerBalance(4,nGRU)        ! inner step balances for domain with one layer
  real(rkind),device                          :: meanBalance(8,nGRU)         ! timestep-average balances for domains
  real(rkind),device,allocatable              :: innerBalanceLayerMass(:,:) ! inner step balances for domain with multiple layers
  real(rkind),device,allocatable              :: innerBalanceLayerNrg(:,:)  ! inner step balances for domain with multiple layers
  ! test balance checks
  logical(lgt),parameter               :: printBalance=.false.   ! flag to print the balance checks
  ! real(rkind),allocatable              :: liqSnowInit(:)         ! volumetric liquid water conetnt of snow at the start of the time step
  ! real(rkind),allocatable              :: liqSoilInit(:)         ! soil moisture at the start of the time step
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
  logical(lgt),device                         :: computeEnthalpy        ! flag to compute enthalpy regardless of the model decision
  logical(lgt),device                         :: enthalpyStateVec       ! flag if enthalpy is a state variable (IDA)
  logical(lgt),device                         :: use_lookup             ! flag to use the lookup table for soil enthalpy, otherwise use analytical solution
  ! ----------------------------------------------------------------------------------------------------------------------------------------------
  type(decisions_device) :: decisions
  type(veg_parameters) :: veg_param
    type(dim3) :: blocks,threads
    logical(lgt),device :: tooMuchSublim_d(nGRU)
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  ! initialize error control
  err=0; message="coupled_em/"

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
    snowfrz_scale     => mpar_data%snowfrz_scale_    ,& ! scaling parameter for the snow freezing curve (K-1)
    canopyDepth       => diag_data%scalarCanopyDepth ,& ! depth of the vegetation canopy (m) 
    specificHeatVeg   => mpar_data%specificHeatVeg_  ,& ! specific heat of vegetation (J kg-1 K-1)
    maxMassVegetation => mpar_data%maxMassVegetation_ & ! maximum mass of vegetation (kg m-2)
    )

    ! define the first step and first and last inner steps
    firstSubStep = .true.
    firstInnerStep = .true.
    lastInnerStep = .false.

    ! create temporary data structures for prognostic variables
    call allocate_device_prog_temp(prog_temp,nGRU,indx_data%nSoil)
    ! call resizeData(prog_meta(:),prog_data,prog_temp,err=err,message=cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! create temporary data structures for diagnostic variables
    call allocate_device_diag_temp(diag_temp,nGRU,indx_data%nSoil)
    ! call resizeData(diag_meta(:),diag_data,diag_temp,err=err,message=cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! create temporary data structures for index variables
    call allocate_device_indx_temp(indx_temp,indx_data,indx_data%nSoil,nGRU)
    call allocate_device_indx_temp(indx_temp0,indx_data,indx_data%nSoil,nGRU)
    ! call resizeData(indx_meta(:),indx_data,indx_temp,err=err,message=cmessage)
    ! call resizeData(indx_meta(:),indx_data,indx_temp0,err=err,message=cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! allocate space for the local fluxes
    call allocate_device_flux_prev(flux_mean,indx_data%nSoil,nGRU)
    ! call allocLocal(flux_meta,flux_mean,nSnow,nSoil,err,cmessage)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
    ! call allocLocal(flux_meta,flux_inner,nSnow,nSoil,err,cmessage)
    call allocate_device_flux_prev(flux_inner,nSoil,nGRU)
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if


tooMuchSublim_d = tooMuchSublim

      ! identify the need to check the mass balance, both methods should work if tolerance coarse enough
      ! select case(ixNumericalMethod)
        ! case(ida);               
    checkMassBalance_ds = .false. ! IDA balance agreement levels are controlled by set tolerances
      !   case(kinsol, homegrown); checkMassBalance_ds = .true.  ! KINSOL or homegrown give finite difference dt_sub fluxes and were summed for an average flux
      !   case default; err=20;    message=trim(message)//'expect num_method to be ida, kinsol, or homegrown (or itertive, which is homegrown)'; return
      ! end select

      ! set the number of substeps for a BE solver
      ! be_steps = NINT(mpar_data%var(iLookPARAM%be_steps)%dat(1)) ! number of substeps for a BE solver
      ! if (be_steps < 1) then
      !   message=trim(message)//'expect be_steps to be greater than 0'
      !   err=20; return
      ! end if
      ! if (ixNumericalMethod == ida) &
      be_steps = 1_i4b ! IDA does not use substeps

    ! short-cut to the algorithmic control parameters
    ! NOTE - temporary assignment of minstep to foce something reasonable
    ! changing the maxstep parameter will make the outer and inner loop computations here in coupled_em happen more frequently
    ! changing the be_steps parameter will make the inner loop computations in opSplittin happen more frequently (e.g. be_steps = 32.0 give BE32)
    minstep = mpar_data%minstep  ! minimum time step (s)
    maxstep = mpar_data%maxstep  ! maximum time step (s)
    maxstep_op = mpar_data%maxstep/be_steps  ! maximum time step (s) to run opSplittin over




    call initialize_coupled_em_kernel<<<blocks,threads>>>(nGRU,flux_data%data,flux_data%ixScalarThroughfallRain,flux_data%ixScalarRainfall,flux_data%ixScalarCanopyLiqDrainage,&
  bal_veg,bal_snow,bal_soil,bal_aq,&
  computeVegFlux,&
  canopyDepth,&
  snowfrz_scale,specificHeatVeg,maxMassVegetation,&
  prog_data%scalarCanopyLiq,prog_data%scalarCanopyIce,prog_data%scalarCanopyTemp,&
  diag_data%scalarBulkVolHeatCapVeg,&
  computeEnthalpy,enthalpyStateVec,&
  diag_data%scalarCanopyEnthTemp,prog_data%scalarCanopyEnthalpy,&
  data_step,exposedVAI,&
  decisions%snowIncept,decisions%snowUnload,&
  forc_data%airtemp_d,&
  mpar_data%refInterceptCapSnow_,mpar_data%ratioDrip2Unloading_,mpar_data%snowUnloadingCoeff_,&
  mpar_data%minTempUnloading_,mpar_data%minWindUnloading_,mpar_data%rateTempUnloading_,mpar_data%rateWindUnloading_,&
  diag_data%scalarNewSnowDensity,prog_data%scalarCanairTemp,&
  flux_data%ixScalarSnowfall,flux_data%ixScalarWindspdCanopyTop,&
  flux_data%ixScalarThroughfallSnow,flux_data%ixScalarCanopySnowUnloading,&
  indx_data%nSnow,indx_data%nSoil,indx_data%nLayers_d,&
  decisions%canopySrad,type_data%vegTypeIndex,&
  prog_data%scalarSWE,prog_data%scalarSnowDepth,&
  prog_data%mLayerVolFracLiq,prog_data%spectralSnowAlbedoDiffuse,&
  prog_data%scalarSnowAlbedo,prog_data%mLayerTemp,&
  diag_data%scalarGroundSnowFraction,diag_data%scalarSnowAge,diag_data%scalarCosZenith,&
  diag_data%spectralSnowAlbedoDirect,diag_data%scalarExposedLAI,diag_data%scalarExposedSAI,&
  diag_data%scalarCanopyWetFraction,diag_data%scalarCanopySunlitFraction,&
  diag_data%scalarCanopySunlitLAI,diag_data%scalarCanopyShadedLAI,&
  diag_data%spectralAlbGndDirect,diag_data%spectralAlbGndDiffuse,diag_data%scalarGroundAlbedo,&
  flux_data%ixSpectralIncomingDirect_start,flux_data%ixSpectralIncomingDirect_end,flux_data%ixSpectralIncomingDiffuse_start,flux_data%ixspectralIncomingDiffuse_end,&
  flux_data%ixScalarCanopySunlitPAR,flux_data%ixscalarCanopyShadedPAR,&
  flux_data%ixSpectralBelowCanopyDirect_start,flux_data%ixSpectralBelowCanopyDirect_end,flux_data%ixSpectralBelowCanopyDiffuse_start,flux_data%ixSpectralBelowCanopyDiffuse_end,&
  flux_data%ixScalarBelowCanopySolar,flux_data%ixScalarCanopyAbsorbedSolar,flux_data%ixScalarGroundAbsorbedSolar,&
  veg_param%ALBSAT,veg_param%ALBDRY,veg_param%RHOL,veg_param%RHOS,veg_param%TAUL,veg_param%TAUS,&
  veg_param%alblak,veg_param%omegas,veg_param%betads,veg_param%betais,veg_param%hvt_,veg_param%hvb_,veg_param%rc,veg_param%opt_rad,veg_param%xl,veg_param%opt_alb,&
  decisions%alb_method,&
  mpar_data%Frad_vis_,mpar_data%Frad_direct_,&
  mpar_data%albedoMax_,mpar_data%albedoMinWinter_,mpar_data%albedoMinSpring_,&
  mpar_data%albedoMaxVisible_,mpar_data%albedoMinVisible_,&
  mpar_data%albedoMaxNearIR_,mpar_data%albedoMinNearIR_,&
  mpar_data%albedoDecayRate_,mpar_data%tempScalGrowth_,&
  mpar_data%albedoSootLoad_,mpar_data%albedoRefresh_,&
  diag_data%scalarCanopyLiqMax,diag_data%scalarCanopyIceMax,&
  mpar_data%canopyWettingFactor_,mpar_data%canopyWettingExp_,&
  computeVegFluxOld,modifiedVegState,&
  mpar_data%refInterceptCapRain_,&
  decisions%bcUpprTdyn,decisions%bcUpprSoiH,&
  fracJulDay,yearLength,&
  attr_data%latitude,&
  mpar_data%z0Snow_,mpar_data%z0Soil_,&
  mpar_data%heightCanopyTop_,mpar_data%heightCanopyBottom_,&
  diag_data%scalarRootZoneTemp,diag_data%scalarLAI,diag_data%scalarSAI,&
  diag_data%scalarGrowingSeasonIndex,&
  veg_param%dveg,veg_param%iswater,veg_param%isbarren,veg_param%issnow,&
  veg_param%laim_,veg_param%saim_,veg_param%tmin,&
  urbanVegCategory,minExpLogHgtFac,&
  prog_data%iLayerHeight,mpar_data%rootingDepth_,diag_data%scalarFoliageNitrogenFactor,oldSWE,&
  prog_data%mLayerDepth,use_lookup,&
  decisions%nrgConserv,&
  scalarInitCanopyLiq,scalarInitCanopyIce,&
  diag_data%scalarTotalSoilLiq,diag_data%scalarTotalSoilIce, &
  prog_data%mLayerVolFracIce, &
  balanceCanopyWater0, balanceSoilWater0, balanceAquifer0, &
  prog_data%scalarAquiferStorage,flux_mean%data,&
  innerEffRainfall, sumCanopySublimation, sumLatHeatCanopyEvap, sumSenHeatCanopy, sumSnowSublimation,&
  innerBalance,innerBalanceLayerMass,innerBalanceLayerNrg,&
  mLayerVolFracIceInit,&
  sfcMeltPond, meanCanopySublimation,meanLatHeatCanopyEvap,meanSenHeatCanopy,&
  effRainfall,meanSoilCompress,meanBalance,modifiedLayers)

  err = cudaDeviceSynchronize()
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

    ! ! initialize if used a balance

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
        if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
        indx_temp0=indx_data!call finalize_device_indx_data(indx_data,indx_temp0)
      endif

      ! check if on outer loop, always do outer if after failed step and on then on reduced whole_step
      do_outer = .false.
      if(stepFailure) firstInnerStep = .true.
      if(firstInnerStep) do_outer = .true.

      if(do_outer)then

        if(.not.stepFailure)then
          ! call resizeData(indx_meta(:),indx_data,indx_temp,err=err,message=cmessage)
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
        endif

        ! save/recover copies of index variables, temp saved on lastInnerStep, failed starts at lastInnerStep
        ! do iVar=1,size(indx_data%var)
          select case(stepFailure)
            case(.false.); indx_temp=indx_data!call finalize_device_indx_data(indx_data,indx_temp)!indx_temp%var(iVar)%dat(:) = indx_data%var(iVar)%dat(:)
            case(.true.);  indx_data=indx_temp!call initialize_device_indx_data(indx_data,indx_temp,indx_data%var(iLookINDEX%nSoil)%dat(1),nGRU)!indx_data%var(iVar)%dat(:) = indx_temp%var(iVar)%dat(:)
          end select
        ! end do  ! looping through variables

        ! save/recover copies of prognostic variables
        ! do iVar=1,size(prog_data%var)
          select case(stepFailure)
            case(.false.); prog_temp = prog_data!call finalize_device_prog_data(prog_data,prog_temp,indx_data%var(iLookINDEX%nLayers)%dat(1),indx_data%var(iLookINDEX%nSoil)%dat(1))!prog_temp%var(iVar)%dat(:) = prog_data%var(iVar)%dat(:)
            case(.true.);  prog_data = prog_temp!call initialize_device_prog_data(prog_data,prog_temp,nGRU,indx_data%var(iLookINDEX%nSoil)%dat(1))!prog_data%var(iVar)%dat(:) = prog_temp%var(iVar)%dat(:)
          end select
        ! end do  ! looping through variables

        ! save/recover copies of diagnostic variables
        ! do iVar=1,size(diag_data%var)
          select case(stepFailure)
            case(.false.); diag_temp=diag_data!call finalize_device_diag_data(diag_data,diag_temp,indx_data%var(iLookINDEX%nSnow)%dat(1),indx_data%var(iLookINDEX%nLayers)%dat(1)) !diag_temp%var(iVar)%dat(:) = diag_data%var(iVar)%dat(:)
            case(.true.);  diag_data=diag_temp!call initialize_device_diag_data(diag_data,diag_temp,nGRU,indx_data%var(iLookINDEX%nSoil)%dat(1)) !diag_data%var(iVar)%dat(:) = diag_temp%var(iVar)%dat(:)
          end select
        ! end do  ! looping through variables

        nSoil = indx_data%nSoil
        ! *** merge/sub-divide snow layers...
        ! -----------------------------------
        ! print*, indx_data%var(iLookINDEX%nSnow)%dat(1), indx_data%var(iLookINDEX%nLayers)%dat(1)
        ! indx_data%var(iLookINDEX%nSnow)%dat(1) = nSnow
        ! indx_data%var(iLookINDEX%nLayers)%dat(1) = nLayers
  err = cudaDeviceSynchronize()

        call resize_snow_kernel<<<blocks,threads>>>(nGRU,decisions%snowLayers,indx_data%nSnow,indx_data%nLayers_d,nSoil, &
 prog_data%mLayerTemp,prog_data%mLayerVolFracIce,prog_data%mLayerVolFracLiq, &
prog_data%mLayerVolFracWat,prog_data%mLayerEnthalpy,prog_data%mLayerDepth, prog_data%mLayerHeight, prog_data%iLayerHeight, &
diag_data%mLayerVolHtCapBulk, diag_data%mLayerCm, diag_data%mLayerThermalC, diag_data%iLayerThermalC, &
diag_data%mLayerEnthTemp, diag_data%mLayerFracLiqSnow, diag_data%mLayerThetaResid, diag_data%mLayerPoreSpace, &
diag_data%mLayerMeltFreeze, diag_data%mLayerVolFracAir, diag_data%balanceLayerNrg, diag_data%balanceLayerMass,&
flux_data%data, &
flux_data%ixILayerConductiveFlux_start,flux_data%ixILayerConductiveFlux_end, &
flux_data%ixILayerAdvectiveFlux_start,flux_data%ixILayerAdvectiveFlux_end, &
flux_data%ixILayerNrgFlux_start,flux_data%ixILayerNrgFlux_end, flux_data%ixMLayerNrgFlux_start,flux_data%ixMLayerNrgFlux_end, &
flux_data%ixILayerLiqFluxSnow_start,flux_data%ixILayerLiqFluxSnow_end,flux_data%ixMLayerLiqFluxSnow_start,flux_data%ixMLayerLiqFluxSnow_end, &
    indx_data%layerType,indx_data%ixHydType, &
    indx_data%ixSnowSoilNrg, indx_data%ixSnowOnlyNrg, &
    indx_data%ixSnowSoilHyd, indx_data%ixSnowOnlyHyd, &
    indx_data%ixNrgLayer, indx_data%ixHydLayer,&
    tooMuchMelt, mpar_data%zMin_,&
    mpar_data%zMinLayer1_,mpar_data%zMinLayer2_,mpar_data%zminLayer3_,mpar_data%zminLayer4_,mpar_data%zminLayer5_, &
    modifiedLayers,snowfrz_scale,lookup_data%h_lookup,lookup_data%t_lookup, &
    prog_data%scalarSnowDepth, prog_data%scalarSWE, &
    mpar_data%zmax_, &
    mpar_data%zmaxLayer1_lower_, mpar_data%zmaxLayer2_lower_, mpar_data%zmaxLayer3_lower_, mpar_data%zmaxLayer4_lower_, &
    mpar_data%zmaxLayer1_upper_, mpar_data%zmaxLayer2_upper_, mpar_data%zmaxLayer3_upper_, mpar_data%zmaxLayer4_upper_, &
    decisions%canopySrad, decisions%alb_method, &
    mpar_data%albedoMax_, mpar_data%albedoMaxVisible_, mpar_data%albedoMaxNearIR_, &
    prog_data%spectralSnowAlbedoDiffuse, diag_data%spectralSnowAlbedoDirect, prog_data%scalarSnowAlbedo, &
    veryBig,mpar_data%Frad_vis_,maxSnowLayers,&
    firstSubStep,modifiedVegState,&
  indx_data%nCasNrg,indx_data%nVegNrg,indx_data%nVegMass,&
  indx_data%nVegState,indx_data%nNrgState,indx_data%nWatState,&
  indx_data%nMatState,indx_data%nMassState,indx_data%nState,&
  indx_data%ixNrgCanair,indx_data%ixNrgCanopy,indx_data%ixHydCanopy,indx_data%ixWatAquifer,&
  indx_data%ixSoilState,&
  indx_data%ixControlVolume,indx_data%ixDomainType,indx_data%ixStateType,&
  computeVegFlux,includeAquifer, &
  mpar_data%vGn_alpha_,mpar_data%vGn_n_,mpar_data%theta_sat_,mpar_data%theta_res_, diag_data%scalarVGn_m, mpar_data%soil_dens_intr_,&
  lookup_data%temperature,lookup_data%psiLiq_int,lookup_data%deriv2,use_lookup,&
  prog_data%mLayerMatricHead,&
  enthalpyStateVec,computeEnthalpy)

  err = cudaDeviceSynchronize()

        ! save the number of snow and soil layers
        nSnow   = maxval(indx_data%nSnow)
        nSoil   = indx_data%nSoil
        nLayers = maxval(indx_data%nLayers_d)

        ! recreate the temporary data structures
        ! NOTE: resizeData(meta, old, new, ..)
        ! if(modifiedVegState .or. modifiedLayers)then

          ! create temporary data structures for prognostic variables
          prog_temp=prog_data!call finalize_device_prog_data(prog_data,prog_temp,nLayers,nSoil)
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

          ! create temporary data structures for diagnostic variables
          diag_temp = diag_data!call finalize_device_diag_data(diag_data,diag_temp,nSnow,nLayers)
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

          ! create temporary data structures for index variables
          indx_temp=indx_data!call finalize_device_indx_data(indx_data,indx_temp)
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
  err = cudaDeviceSynchronize()

          ! do iVar=1,size(indx_data%var)
          !   select case(stepFailure)
          !     case(.false.); indx_temp%var(iVar)%dat(:) = indx_data%var(iVar)%dat(:)
          !     case(.true.);  indx_data%var(iVar)%dat(:) = indx_temp%var(iVar)%dat(:)
          !   end select
          ! end do  ! looping through variables

        ! endif  ! if modified the states

        ! define the number of state variables
        nState = maxval(indx_data%nState)
        nState = (maxSnowLayers + nSoil)*2 + 4

        ! *** compute diagnostic variables for each layer...
        ! --------------------------------------------------
        ! NOTE: this needs to be done AFTER volicePack, since layers may have been sub-divided and/or merged, and need to specifically send in canopy depth
    call update_diagnostics<<<blocks,threads>>>(nGRU,decisions%thCondSnow,decisions%thCondSoil,&
    computeVegFlux,canopyDepth,&
    specificHeatVeg,maxMassVegetation,mpar_data%fixedThermalCond_snow_, &
    mpar_data%soil_dens_intr_,mpar_data%thCond_soil_,mpar_data%theta_sat_,mpar_data%frac_sand_,mpar_data%frac_silt_,mpar_data%frac_clay_,&
    indx_data%nSnow,indx_data%nSoil,indx_data%nLayers_d,indx_data%layerType,&
    prog_data%scalarCanopyIce,prog_data%scalarCanopyLiq,&
    prog_data%mLayerVolFracIce,prog_data%mLayerVolFracLiq,&
    prog_data%mLayerHeight,prog_data%iLayerHeight,&
    diag_data%scalarBulkVolHeatCapVeg, &
    diag_data%mLayerVolHtCapBulk,diag_data%mLayerThermalC,diag_data%iLayerThermalC,diag_data%mLayerVolFracAir,&
    prog_data%scalarSWE, prog_data%scalarSnowDepth, prog_data%scalarSfcMeltPond,&
    prog_data%mLayerTemp,prog_data%mLayerDepth,mLayerVolFracIceInit,&
    prog_data%scalarCanopyWat,prog_data%mLayerVolFracWat,&
    diag_data%mLayerEnthTemp,prog_data%mLayerEnthalpy,&
    mpar_data%soil_dens_intr_,mpar_data%vGn_alpha_,mpar_data%vGn_n_,mpar_data%theta_res_,diag_data%scalarVGn_m,&
    prog_data%mLayerMatricHead,use_lookup,&
    lookup_data%temperature,lookup_data%psiLiq_int,lookup_data%deriv2,&
    computeEnthalpy, enthalpyStateVec,&
    diag_data%mLayerMatricHeadLiq,&
    sumCanopySublimation, sumSnowSublimation, sumLatHeatCanopyEvap, sumSenHeatCanopy,&
    flux_inner%data,&
    innerEffRainfall,innerSoilCompress, &
    innerBalance,innerBalanceLayerNrg,innerBalanceLayerMass)
  err = cudaDeviceSynchronize()


        ! correct increments (if need to redo inner step) and reset increment
        dt_solv = dt_solv - dt_solvInner
        dt_solvInner = 0._rkind
        lastInnerStep = .false.

        ! initialize sublimation sums to average over whole_step
        sumStepSize= 0._rkind ! initialize the sum of the step sizes

      endif ! (do_outer loop)
  err = cudaDeviceSynchronize()

      ! *** solve model equations...
      ! ----------------------------
      ! save input step
      dtSave = whole_step

      ! get the new solution
      call opSplittin(&
      nGRU, &
                      ! input: model control
                      nSnow,                                  & ! intent(in):    number of snow layers
                      nSoil,                                  & ! intent(in):    number of soil layers
                      nLayers,                                & ! intent(in):    total number of layers
                      nState,                                 & ! intent(in):    total number of layers
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
  err = cudaDeviceSynchronize()

      ! process the flag for too much melt
      if(tooMuchMelt)then
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
          print*, 'ixSolution', ixSolution
          print*, 'dtSave, dt_sub', dtSave, whole_step
          message=trim(message)//'length of the coupled step is below the minimum step length'
          err=20; return
        endif
        ! try again, restart step
        ! deallocate(mLayerVolFracIceInit)
        ! deallocate(innerBalanceLayerNrg)
        ! deallocate(innerBalanceLayerMass)
        cycle substeps
      endif
  err = cudaDeviceSynchronize()

      call increment_sums_kernel<<<blocks,threads>>>(nGRU, &
  flux_data%data, &
  sumCanopySublimation, sumSnowSublimation, &
  sumLatHeatCanopyEvap, sumSenHeatCanopy, &
  flux_data%ixScalarCanopySublimation, flux_data%ixScalarSnowSublimation, &
  flux_data%ixScalarLatHeatCanopyEvap, flux_data%ixScalarSenHeatCanopy, &
  dt_sub)
  err = cudaDeviceSynchronize()

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

          call sublime_kernel<<<blocks,threads>>>(nGRU,indx_data%nSnow,indx_data%nLayers_d, &
  mLayerMeltFreeze, mLayerVolFracIce, mLayerVolFracIceInit, &
  computeVegFlux, computeEnthalpy,enthalpyStateVec,&
  scalarCanopyIce, sumCanopySublimation, scalarCanopyLiq, &
  superflousSub, superflousNrg, &
  sumLatHeatCanopyEvap, sumSenHeatCanopy, &
  scalarCanopyWat, canopyDepth, &
  specificHeatVeg, maxMassVegetation, &
  whole_step, snowfrz_scale, &
  prog_data%scalarCanopyTemp, diag_data%scalarCanopyEnthTemp,prog_data%scalarCanopyEnthalpy, &
  sumSnowSublimation, &
  mLayerVolFracLiq, prog_data%mLayerTemp, &
  mpar_data%densScalGrowth_,mpar_data%tempScalGrowth_, &
  mpar_data%grainGrowthRate_,mpar_data%densScalOvrbdn_, &
  mpar_data%tempScalOvrbdn_,mpar_data%baseViscosity_,mLayerDepth,&
  tooMuchSublim_d,prog_data%mLayerHeight,prog_data%iLayerHeight,indx_data%layerType, &
  prog_data%scalarSnowDepth, prog_data%scalarSWE, &
  prog_data%mLayerEnthalpy, diag_data%mLayerEnthTemp, &
  sfcMeltPond, prog_data%scalarSfcMeltPond, &
  mLayerVolFracWat)

          ! deallocate(mLayerVolFracIceInit)

  err = cudaDeviceSynchronize()

          if(err/=0)then; err=55; return; end if
          tooMuchSublim = tooMuchSublim_d(1)

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
              print*, 'ixSolution', ixSolution
              print*, 'dtSave, dt_sub', dtSave, whole_step
              message=trim(message)//'length of the coupled step is below the minimum step length'
              err=20; return
            endif
            ! try again, restart step (at end inner step)
            ! deallocate(innerBalanceLayerNrg)
            ! deallocate(innerBalanceLayerMass)
            cycle substeps
          endif


        end associate sublime


      endif ! (do_outer loop)


      ! ****************************************************************************************************
      ! *** END MAIN SOLVER ********************************************************************************
      ! ****************************************************************************************************

      ! increment mean fluxes, soil compression, and effective rainfall, reset on whole_step
      dt_wght = dt_sub/whole_step ! define weight applied to each sub-step
      call increment_means_kernel<<<blocks,threads>>>(nGRU, &
      flux_inner%data, flux_data%data,&
      innerSoilCompress,diag_data%mLayerCompress,&
      indx_data%nSnow,flux_data%numScalarFluxData,dt_wght,&
      innerEffRainfall,flux_data%ixScalarThroughfallRain,flux_data%ixScalarCanopyLiqDrainage, &
  innerBalance, innerBalanceLayerMass, innerBalanceLayerNrg, &
  computeVegFlux, bal_veg, &
  diag_data%balanceCasNrg, diag_data%balanceVegNrg, diag_data%balanceVegMass, canopyDepth, diag_data%balanceAqMass, &
  diag_data%balanceLayerMass, diag_data%balanceLayerNrg, &
  prog_data%mLayerDepth, indx_data%layerType, &
  diag_data%balanceSnowNrg, diag_data%balanceSnowMass, diag_data%balanceSoilNrg,diag_data%balanceSoilMass, &
  bal_soil, bal_snow, bal_aq, &
  indx_data%nLayers_d, decisions%groundwatr)
  err = cudaDeviceSynchronize()

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
        call mean_inner_step_kernel<<<blocks,threads>>>(nGRU,flux_mean%data,flux_inner%data, &
  meanBalance, innerBalance, flux_inner%numScalarFluxData, dt_wght, &
  meanCanopySublimation, sumCanopySublimation, flux_mean%ixScalarCanopySublimation,&
  meanLatHeatCanopyEvap, sumLatHeatCanopyEvap, flux_mean%ixScalarLatHeatCanopyEvap,&
  meanSenHeatCanopy, sumSenHeatCanopy, flux_mean%ixScalarSenHeatCanopy, &
  meanSoilCompress, innerSoilCompress, &
  diag_data%balanceSnowNrg, diag_data%balanceSnowMass, diag_data%balanceSoilNrg, diag_data%balanceSoilMass, &
  effRainfall, innerEffRainfall,data_step)
  err = cudaDeviceSynchronize()

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
  err = cudaDeviceSynchronize()

      ! adjust length of the sub-step (make sure that we don't exceed the step)
      dt_sub = min(data_step - dt_solv, dt_sub)
    end do  substeps ! (sub-step loop)

    diag_data%meanStepSize = diag_data%meanStepSize/nsub_success

    call finalize_coupled_em<<<blocks,threads>>>(nGRU,diag_data%scalarSnowfallTemp, diag_data%scalarNewSnowDensity, &
  flux_data%data,flux_data%ixScalarThroughfallSnow,flux_data%ixScalarCanopySnowUnloading,&
  prog_data%scalarSWE,prog_data%scalarSnowDepth,&
  prog_data%mLayerTemp,prog_data%mLayerDepth,&
  prog_data%mLayerVolFracIce,prog_data%mLayerVolFracLiq,&
  data_step,indx_data%nSnow,snowfrz_scale,&
  prog_data%mLayerVolFracWat,diag_data%mLayerEnthTemp,prog_data%mLayerEnthalpy,&
  computeEnthalpy,enthalpyStateVec,&
  indx_data%nLayers_d,indx_data%layerType,indx_data%nSoil,&
  prog_data%mLayerHeight,prog_data%iLayerHeight,&
  flux_data%data,flux_data%numScalarFluxData,&
  diag_data%mLayerCompress,meanSoilCompress,diag_data%scalarSoilCompress,&
  prog_data%scalarSfcMeltPond,sfcMeltPond,&
  meanBalance,&
  diag_data%balanceCasNrg,diag_data%balanceVegNrg,diag_data%balanceVegMass,diag_data%balanceAqMass,&
  diag_data%balanceSnowNrg,diag_data%balanceSoilNrg,diag_data%balanceSnowMass,diag_data%balanceSoilMass, &
  bal_veg,bal_soil,bal_snow,bal_aq,&
  computeVegFlux,&
  scalarCanopyWatBalError,prog_data%scalarCanopyWat,balanceCanopyWater0,&
  flux_data%ixScalarSnowfall,flux_mean%ixScalarThroughfallSnow,flux_mean%ixScalarRainfall,flux_mean%ixScalarThroughfallRain,&
  flux_mean%ixScalarCanopySnowUnloading,flux_mean%ixScalarCanopyLiqDrainage,flux_mean%ixScalarCanopySublimation,flux_mean%ixScalarCanopyEvaporation,&
  effSnowfall,delSWE,oldSWE,&
  massBalance,effRainfall,&
  flux_mean%ixScalarSnowSublimation,flux_mean%ixScalarSnowDrainage,&
  diag_data%scalarTotalSoilLiq,diag_data%scalarTotalSoilIce,diag_data%scalarTotalSoilWat,&
  balanceAquifer1,prog_data%scalarAquiferStorage,&
  balanceSoilInflux,flux_mean%ixScalarInfiltration,&
  balanceSoilBaseflow,flux_mean%ixScalarSoilBaseflow,&
  balanceSoilDrainage,flux_mean%ixscalarSoilDrainage, &
  balanceSoilET,flux_mean%ixScalarCanopyTranspiration,flux_mean%ixScalarGroundEvaporation, &
  balanceSoilCompress, &
  scalarSoilWatBalError,balanceSoilWater0,&
  prog_data%scalarCanopyEnthalpy,diag_data%scalarCanopyEnthTemp,canopyDepth,&
  indx_data%ixDomainType,indx_data%ixControlVolume,indx_data%ixStateType,&
  prog_data%scalarCanopyIce,diag_data%scalarTotalSoilEnthalpy,diag_data%scalarTotalSnowEnthalpy,&
  prog_data%scalarSurfaceTemp)
  
  err = cudaDeviceSynchronize()

tooMuchSublim = tooMuchSublim_d(1)

    deallocate(innerSoilCompress)
    deallocate(meanSoilCompress)


  end associate canopy  ! end association to canopy parameters


  if(nsub>50000)then
    write(message,'(a,i0)') trim(cmessage)//'number of sub-steps > 50000 for HRU ', hruId
    err=20; return
  end if

  ! get the end time
  CALL system_clock(i_end)
  elapsed_time = REAL(i_end - i_start) / REAL(count_rate)

  ! get the elapsed time
  diag_data%wallClockTime = elapsed_time
  call deallocate_device_indx_data(indx_temp)
  call deallocate_device_indx_data(indx_temp0)
  call deallocate_device_prog_data(prog_temp)
  call deallocate_device_diag_data(diag_temp)
  call deallocate_device_flux_data(flux_mean)
  call deallocate_device_flux_data(flux_inner)

contains

 subroutine initialize_coupled_em
  ! *** Initialize steps for coupled_em subroutine ***
  ! Notes: - created to ensure certain variables are initialized prior to use in calculations
  !        - based on warnings from the SUMMA debug build (e.g., -Wall flag)
  !        - additional initial operations may be added here in the future

  ! initialize variables
  sumStepSize=0._rkind            ! sum solution step for the data step

  ! get initial value of nLayers
  nSnow = maxval(indx_data%nSnow)
  nSoil = indx_data%nSoil
  nLayers = nSnow + nSoil

  ! allocate and initialize using the initial value of nLayers
  allocate(innerBalanceLayerMass(nSoil+maxSnowLayers+1,nGRU)); innerBalanceLayerMass = 0._rkind ! mean total balance of mass in layers
  allocate(innerBalanceLayerNrg(nSoil+maxSnowLayers+1,nGRU));  innerBalanceLayerNrg = 0._rkind ! mean total balance of energy in layers
  allocate(mLayerVolFracIceInit(nSoil+maxSnowLayers+1,nGRU));  mLayerVolFracIceInit = prog_data%mLayerVolFracIce ! volume fraction of water ice

   ! initialize the numerix tracking variables
  indx_data%numberFluxCalc        = 0  ! number of flux calculations                     (-)
  indx_data%numberStateSplit      = 0  ! number of state splitting solutions             (-)
  indx_data%numberDomainSplitNrg  = 0  ! number of domain splitting solutions for energy (-)
  indx_data%numberDomainSplitMass = 0  ! number of domain splitting solutions for mass   (-)
  indx_data%numberScalarSolutions = 0  ! number of scalar solutions                      (-)

  diag_data%meanStepSize = 0._rkind ! mean step size over data_step

  ! Need mean soil compression for balance checks but it is not in flux structure so handle differently 
  !  This will be a problem if nSoil changes (currently not possible)-- then might need to not keep the average
  allocate(meanSoilCompress(nSoil,nGRU))
  allocate(innerSoilCompress(nSoil,nGRU))

  ! start by assuming that the step is successful
  stepFailure  = .false.
  doLayerMerge = .false.


 end subroutine initialize_coupled_em

end subroutine coupled_em


! *********************************************************************************************************
! private subroutine implctMelt: compute melt of the "snow without a layer"
! *********************************************************************************************************
attributes(device) subroutine implctMelt(&
                      ! input/output: integrated snowpack properties
                      scalarSWE,         & ! intent(inout): snow water equivalent (kg m-2)
                      scalarSnowDepth,   & ! intent(inout): snow depth (m)
                      scalarSfcMeltPond, & ! intent(out):   surface melt pond (kg m-2)
                      ! input/output: properties of the upper-most soil layer
                      soilTemp,          & ! intent(inout): surface layer temperature (K)
                      soilDepth,         & ! intent(inout): surface layer depth (m)
                      soilHeatcap,       & ! intent(in):    surface layer volumetric heat capacity (J m-3 K-1)
                      ! output: error control
                      err        ) ! intent(out): error control
  implicit none
  ! input/output: integrated snowpack properties
  real(rkind),intent(inout)    :: scalarSWE          ! snow water equivalent (kg m-2)
  real(rkind),intent(inout)    :: scalarSnowDepth    ! snow depth (m)
  real(rkind),intent(out)      :: scalarSfcMeltPond  ! surface melt pond (kg m-2)
  ! input/output: properties of the upper-most soil layer
  real(rkind),intent(inout)    :: soilTemp           ! surface layer temperature (K)
  real(rkind),intent(inout)    :: soilDepth          ! surface layer depth (m)
  real(rkind),intent(in)       :: soilHeatcap        ! surface layer volumetric heat capacity (J m-3 K-1)
  ! output: error control
  integer(i4b),intent(out)     :: err                ! error code
  ! character(*),intent(out)     :: message            ! error message
  ! local variables
  real(rkind)                  :: nrgRequired        ! energy required to melt all the snow (J m-2)
  real(rkind)                  :: nrgAvailable       ! energy available to melt the snow (J m-2)
  real(rkind)                  :: snwDensity         ! snow density (kg m-3)
  ! initialize error control
  err=0!; message='implctMelt/'

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

attributes(global) subroutine initialize_coupled_em_kernel(nGRU,flux_data,scalarThroughfallRain,scalarRainfall,scalarCanopyLiqDrainage,&
  bal_veg,bal_snow,bal_soil,bal_aq,&
  computeVegFlux,&
  canopyDepth,&
  snowfrz_scale,specificHeatVeg,maxMassVegetation,&
  scalarCanopyLiq,scalarCanopyIce,scalarCanopyTemp,&
  scalarBulkVolHeatCapVeg,&
  computeEnthalpy,enthalpyStateVec,&
  scalarCanopyEnthTemp,scalarCanopyEnthalpy,&
  data_step,exposedVAI,&
  snowIncept,snowUnload,&
  airtemp,&
  refInterceptCapSnow,ratioDrip2Unloading,snowUnloadingCoeff,&
  minTempUnloading,minWindUnloading,rateTempUnloading,rateWindUnloading,&
  scalarNewSnowDensity,scalarCanairTemp,&
  scalarSnowfall,scalarWindspdCanopyTop,&
  scalarThroughfallSnow,scalarCanopySnowUnloading,&
  nSnow,nSoil,nLayers,&
  ix_canopySrad,vegTypeIndex,&
  scalarSWE,scalarSnowDepth,&
  mLayerVolFracLiq,spectralSnowAlbedoDiffuse,&
  scalarSnowAlbedo,mLayerTemp,&
  scalarGroundSnowFraction,scalarSnowAge,scalarCosZenith,&
  spectralSnowAlbedoDirect,scalarExposedLAI,scalarExposedSAI,&
  scalarCanopyWetFraction,scalarCanopySunlitFraction,&
  scalarCanopySunlitLAI,scalarCanopyShadedLAI,&
  spectralAlbGndDirect,spectralAlbGndDiffuse,scalarGroundAlbedo,&
  spectralIncomingDirect_start,spectralIncomingDirect_end,spectralIncomingDiffuse_start,spectralIncomingDiffuse_end,&
  scalarCanopySunlitPAR,scalarCanopyShadedPAR,&
  spectralBelowCanopyDirect_start,spectralBelowCanopyDirect_end,spectralBelowCanopyDiffuse_start,spectralBelowCanopyDiffuse_end,&
  scalarBelowCanopySolar,scalarCanopyAbsorbedSolar,scalarGroundAbsorbedSolar,&
  ALBSAT,ALBDRY,RHOL,RHOS,TAUL,TAUS,&
  alblak,omegas,betads,betais,hvt,hvb,rc,opt_rad,xl,opt_alb,&
  ixAlbedoMethod,&
  Frad_vis,Frad_direct,&
  albedoMax,albedoMinWinter,albedoMinSpring,&
  albedoMaxVisible,albedoMinVisible,&
  albedoMaxNearIR,albedoMinNearIR,&
  albedoDecayRate,tempScalGrowth,&
  albedoSootLoad,albedoRefresh,&
  scalarCanopyLiqMax,scalarCanopyIceMax,&
  canopyWettingFactor,canopyWettingExp,&
  computeVegFluxOld,modifiedVegState,&
  refInterceptCapRain,&
  bcUpprTdyn,bcUpprSoiH,&
  fracJulDay,yearLength,&
  latitude,&
  z0Snow,z0Soil,&
  heightCanopyTop,heightCanopyBottom,&
  scalarRootZoneTemp,scalarLAI,scalarSAI,&
  scalarGrowingSeasonIndex,&
  dveg,iswater,isbarren,issnow,&
  laim,saim,tmin,&
  urbanVegCategory,minExpLogHgtFac,&
  iLayerHeight,rootingDepth,scalarFoliageNitrogenFactor,oldSWE,&
  mLayerDepth,use_lookup,&
  ixNrgConserv,&
  scalarInitCanopyLiq,scalarInitCanopyIce,&
  scalarTotalSoilLiq,scalarTotalSoilIce, &
  mLayerVolFracIce, &
  balanceCanopyWater0, balanceSoilWater0, balanceAquifer0, &
  scalarAquiferStorage,flux_mean,&
  innerEffRainfall, sumCanopySublimation, sumLatHeatCanopyEvap, sumSenHeatCanopy, sumSnowSublimation,&
  innerBalance,innerBalanceLayerMass,innerBalanceLayerNrg,&
  mLayerVolFracIceInit,&
  sfcMeltPond, meanCanopySublimation,meanLatHeatCanopyEvap,meanSenHeatCanopy,&
  effRainfall,meanSoilCompress,meanBalance,modifiedLayers)
  integer(i4b),intent(in),value :: nGRU
  real(rkind),intent(inout) :: flux_data(:,:)
    integer(i4b),intent(in),value :: scalarThroughfallRain,scalarRainfall,scalarCanopyLiqDrainage
  logical(lgt),intent(inout) :: bal_veg(:),bal_snow(:),bal_soil(:),bal_aq(:)
  logical(lgt),intent(inout) :: computeVegFlux(:)
  real(rkind),intent(inout) :: canopyDepth(:)
  real(rkind),intent(in) :: snowfrz_scale(:),specificHeatVeg(:),maxMassVegetation(:)
  real(rkind),intent(inout) :: scalarCanopyLiq(:),scalarCanopyIce(:),scalarCanopyTemp(:)
  real(rkind),intent(inout) :: scalarBulkVolHeatCapVeg(:)
  logical(lgt),intent(inout) :: computeEnthalpy,enthalpyStateVec
  real(rkind),intent(inout) :: scalarCanopyEnthTemp(:),scalarCanopyEnthalpy(:)
  real(rkind),intent(in),value :: data_step
  real(rkind),intent(inout) :: exposedVAI(:)
  integer(i4b),intent(in) :: snowIncept,snowUnload
  real(rkind),intent(inout) :: airtemp(:)
  real(rkind),intent(inout) :: refInterceptCapSnow(:),ratioDrip2Unloading(:),snowUnloadingCoeff(:)
  real(rkind),intent(inout) :: minTempUnloading(:),minWindUnloading(:),rateTempUnloading(:),rateWindUnloading(:)
  real(rkind),intent(inout) :: scalarNewSnowDensity(:),scalarCanairTemp(:)
  integer(i4b),intent(in),value :: scalarSnowfall,scalarWindspdCanopyTop
  integer(i4b),intent(in),value :: scalarThroughfallSnow,scalarCanopySnowUnloading
  integer(i4b),intent(inout) :: nSnow(:),nLayers(:)
  integer(i4b),intent(in),value :: nSoil
  integer(i4b),intent(in) :: ix_canopySrad
  integer(i4b),intent(in) :: vegTypeIndex(:)
  real(rkind),intent(inout) :: scalarSWE(:),scalarSnowDepth(:)
  real(rkind),intent(inout) :: mLayerVolFracLiq(:,:), spectralSnowAlbedoDiffuse(:,:)
  real(rkind),intent(inout) :: scalarSnowAlbedo(:),mLayerTemp(:,:)
  real(rkind),intent(inout) :: scalarGroundSnowFraction(:),scalarSnowAge(:),scalarCosZenith(:)
  real(rkind),intent(inout) :: spectralSnowAlbedoDirect(:,:),scalarExposedLAI(:),scalarExposedSAI(:)
  real(rkind),intent(inout) :: scalarCanopyWetFraction(:),scalarCanopySunlitFraction(:)
  real(rkind),intent(inout) :: scalarCanopySunlitLAI(:),scalarCanopyShadedLAI(:)
  real(rkind),intent(inout) :: spectralAlbGndDirect(:,:), spectralAlbGndDiffuse(:,:),scalarGroundAlbedo(:)
  integer(i4b),intent(in),value :: spectralIncomingDirect_start,spectralIncomingDirect_end,spectralIncomingDiffuse_start,spectralIncomingDiffuse_end
  integer(i4b),intent(in),value :: scalarCanopySunlitPAR,scalarCanopyShadedPAR
  integer(i4b),intent(in),value :: spectralBelowCanopyDirect_start,spectralBelowCanopyDirect_end, spectralBelowCanopyDiffuse_start,spectralBelowCanopyDiffuse_end
  integer(i4b),intent(in),value :: scalarBelowCanopySolar,scalarCanopyAbsorbedSolar,scalarGroundAbsorbedSolar
 real(rkind),intent(in) :: ALBSAT(:,:),ALBDRY(:,:),RHOL(:,:),RHOS(:,:),TAUL(:,:),TAUS(:,:)
   real(rkind),intent(in) :: alblak(:)
   real(rkind),intent(in) :: omegas(:), betads,betais
   real(rkind),intent(in) :: hvt(:,:),hvb(:,:),rc(:),xl(:)
   integer(i4b),intent(in) :: opt_alb,opt_rad
   integer(i4b),intent(in) :: ixAlbedoMethod
   real(rkind),intent(in) :: Frad_vis(:),Frad_direct(:)
   real(rkind),intent(in) :: albedoMax(:),albedoMinWinter(:),albedoMinSpring(:)
   real(rkind),intent(in) :: albedoMaxVisible(:),albedoMinVisible(:)
   real(rkind),intent(in) :: albedoMaxNearIR(:),albedoMinNearIR(:)
   real(rkind),intent(in) :: albedoDecayRate(:),tempScalGrowth(:)
   real(rkind),intent(in) :: albedoSootLoad(:),albedoRefresh(:)
   real(rkind),intent(inout) :: scalarCanopyLiqMax(:),scalarCanopyIceMax(:)
   real(rkind),intent(in) :: canopyWettingFactor(:),canopyWettingExp(:)
   logical(lgt),intent(inout) :: computeVegFluxOld(:), modifiedVegState(:)
   real(rkind),intent(in) :: refInterceptCapRain(:)
   integer(i4b),intent(in) :: bcUpprTdyn,bcUpprSoiH
   real(rkind),intent(in),value :: fracJulDay
   integer(i4b),intent(in),value :: yearLength
   real(rkind),intent(in) :: latitude(:)
   real(rkind),intent(in) :: z0Snow(:),z0Soil(:)
   real(rkind),intent(in) :: heightCanopyTop(:),heightCanopyBottom(:)
   real(rkind),intent(inout) :: scalarRootZoneTemp(:),scalarLAI(:),scalarSAI(:)
   real(rkind),intent(inout) :: scalarGrowingSeasonIndex(:)
   integer(i4b),intent(in) :: dveg,iswater,isbarren,issnow
   real(rkind),intent(in) :: laim(:,:,:), saim(:,:,:), tmin(:)
   integer(i4b),intent(in),value :: urbanVegCategory
   real(rkind),intent(in),value :: minExpLogHgtFac
   real(rkind),intent(inout) :: iLayerHeight(0:,:)
   real(rkind),intent(inout) :: rootingDepth(:),scalarFoliageNitrogenFactor(:),oldSWE(:)
   real(rkind),intent(inout) :: mLayerDepth(:,:)
   logical(lgt),intent(inout) :: use_lookup
   integer(i4b),intent(in) :: ixNrgConserv
   real(rkind),intent(inout) :: scalarInitCanopyLiq(:),scalarInitCanopyIce(:)
   real(rkind),intent(inout) :: scalarTotalSoilLiq(:),scalarTotalSoilIce(:)
   real(rkind),intent(inout) :: mLayerVolFracIce(:,:)
   real(rkind),intent(inout) :: balanceCanopyWater0(:), balanceSoilWater0(:), balanceAquifer0(:)
   real(rkind),intent(inout) :: scalarAquiferStorage(:)
   real(rkind),intent(inout) :: flux_mean(:,:)
   real(rkind),intent(inout) :: innerEffRainfall(:), sumCanopySublimation(:), sumLatHeatCanopyEvap(:), sumSenHeatCanopy(:), sumSnowSublimation(:)
   real(rkind),intent(inout) :: innerBalance(:,:), innerBalanceLayerMass(:,:), innerBalanceLayerNrg(:,:)
   real(rkind),intent(inout) :: mLayerVolFracIceInit(:,:)
   real(rkind),intent(inout) :: sfcMeltPond(:), meanCanopySublimation(:),meanLatHeatCanopyEvap(:),meanSenHeatCanopy(:)
   real(rkind),intent(inout) :: effRainfall(:)
   real(rkind),intent(inout) :: meanSoilCompress(:,:), meanBalance(:,:)
   logical(lgt),intent(inout) :: modifiedLayers(:)

         integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return
  ! do iGRU=1,nGRU
  call initialize_coupled_em_device(flux_data(scalarThroughfallRain,iGRU),flux_data(scalarRainfall,iGRU),flux_data(scalarCanopyLiqDrainage,iGRU),&
  bal_veg(iGRU),bal_snow(iGRU),bal_soil(iGRU),bal_aq(iGRU),&
  computeVegFlux(iGRU),&
  canopyDepth,&
  snowfrz_scale,specificHeatVeg,maxMassVegetation,&
  scalarCanopyLiq,scalarCanopyIce,scalarCanopyTemp,&
  scalarBulkVolHeatCapVeg,&
  computeEnthalpy,enthalpyStateVec,&
  scalarCanopyEnthTemp(iGRU),scalarCanopyEnthalpy(iGRU),&
  data_step,exposedVAI(iGRU),&
  snowIncept,snowUnload,&
  airtemp(iGRU),&
  refInterceptCapSnow(iGRU),ratioDrip2Unloading(iGRU),snowUnloadingCoeff(iGRU),&
  minTempUnloading(iGRU),minWindUnloading(iGRU),rateTempUnloading(iGRU),rateWindUnloading(iGRU),&
  scalarNewSnowDensity(iGRU),scalarCanairTemp(iGRU),&
  flux_data(scalarSnowfall,iGRU),flux_data(scalarWindspdCanopyTop,iGRU),&
  flux_data(scalarThroughfallSnow,iGRU),flux_data(scalarCanopySnowUnloading,iGRU),&
  nSnow(iGRU),nSoil,nLayers(iGRU),&
  ix_canopySrad,vegTypeIndex(iGRU),&
  scalarSWE(iGRU),scalarSnowDepth(iGRU),&
  mLayerVolFracLiq(:,iGRU),spectralSnowAlbedoDiffuse(:,iGRU),&
  scalarSnowAlbedo(iGRU),mLayerTemp(:,iGRU),&
  scalarGroundSnowFraction(iGRU),scalarSnowAge(iGRU),scalarCosZenith(iGRU),&
  spectralSnowAlbedoDirect(:,iGRU),scalarExposedLAI(iGRU),scalarExposedSAI(iGRU),&
  scalarCanopyWetFraction(iGRU),scalarCanopySunlitFraction(iGRU),&
  scalarCanopySunlitLAI(iGRU),scalarCanopyShadedLAI(iGRU),&
  spectralAlbGndDirect(:,iGRU),spectralAlbGndDiffuse(:,iGRU),scalarGroundAlbedo(iGRU),&
  flux_data(spectralIncomingDirect_start:spectralIncomingDirect_end,iGRU),flux_data(spectralIncomingDiffuse_start:spectralIncomingDiffuse_end,iGRU),&
  flux_data(scalarCanopySunlitPAR,iGRU),flux_data(scalarCanopyShadedPAR,iGRU),&
  flux_data(spectralBelowCanopyDirect_start:spectralBelowCanopyDirect_end,iGRU),flux_data(spectralBelowCanopyDiffuse_start:spectralBelowCanopyDiffuse_end,iGRU),&
  flux_data(scalarBelowCanopySolar,iGRU),flux_data(scalarCanopyAbsorbedSolar,iGRU),flux_data(scalarGroundAbsorbedSolar,iGRU),&
  ALBSAT,ALBDRY,RHOL,RHOS,TAUL,TAUS,&
  alblak,omegas,betads,betais,hvt(:,iGRU),hvb(:,iGRU),rc,opt_rad,xl,opt_alb,&
  ixAlbedoMethod,&
  Frad_vis(iGRU),Frad_direct(iGRU),&
  albedoMax(iGRU),albedoMinWinter(iGRU),albedoMinSpring(iGRU),&
  albedoMaxVisible(iGRU),albedoMinVisible(iGRU),&
  albedoMaxNearIR(iGRU),albedoMinNearIR(iGRU),&
  albedoDecayRate(iGRU),tempScalGrowth(iGRU),&
  albedoSootLoad(iGRU),albedoRefresh(iGRU),&
  scalarCanopyLiqMax(iGRU),scalarCanopyIceMax(iGRU),&
  canopyWettingFactor(iGRU),canopyWettingExp(iGRU),&
  computeVegFluxOld(iGRU),modifiedVegState(iGRU),&
  refInterceptCapRain(iGRU),&
  bcUpprTdyn,bcUpprSoiH,&
  fracJulDay,yearLength,&
  latitude(iGRU),&
  z0Snow(iGRU),z0Soil(iGRU),&
  heightCanopyTop(iGRU),heightCanopyBottom(iGRU),&
  scalarRootZoneTemp(iGRU),scalarLAI(iGRU),scalarSAI(iGRU),&
  scalarGrowingSeasonIndex(iGRU),&
  dveg,iswater,isbarren,issnow,&
  laim(:,:,iGRU),saim(:,:,iGRU),tmin,&
  urbanVegCategory,minExpLogHgtFac,&
  iLayerHeight(:,iGRU),rootingDepth(iGRU),scalarFoliageNitrogenFactor(iGRU),oldSWE(iGRU),&
  mLayerDepth(:,iGRU),use_lookup,&
  ixNrgConserv,&
  scalarInitCanopyLiq(iGRU),scalarInitCanopyIce(iGRU),&
  scalarTotalSoilLiq(iGRU),scalarTotalSoilIce(iGRU), &
  mLayerVolFracIce(:,iGRU), &
  balanceCanopyWater0(iGRU), balanceSoilWater0(iGRU), balanceAquifer0(iGRU), &
  scalarAquiferStorage(iGRU),flux_mean(:,iGRU),&
  innerEffRainfall(iGRU), sumCanopySublimation(iGRU), sumLatHeatCanopyEvap(iGRU), sumSenHeatCanopy(iGRU), sumSnowSublimation(iGRU),&
  innerBalance(:,iGRU),innerBalanceLayerMass(:,iGRU),innerBalanceLayerNrg(:,iGRU),&
  mLayerVolFracIceInit(:,iGRU),&
  sfcMeltPond(iGRU), meanCanopySublimation(iGRU),meanLatHeatCanopyEvap(iGRU),meanSenHeatCanopy(iGRU),&
  effRainfall(iGRU),meanSoilCompress(:,iGRU),meanBalance(:,iGRU),modifiedLayers(iGRU))
  ! end do
end subroutine

attributes(device) subroutine initialize_coupled_em_device(scalarThroughfallRain,scalarRainfall,scalarCanopyLiqDrainage,&
  bal_veg,bal_snow,bal_soil,bal_aq,&
  computeVegFlux,&
  canopyDepth_,&
  snowfrz_scale_,specificHeatVeg_,maxMassVegetation_,&
  scalarCanopyLiq_,scalarCanopyIce_,scalarCanopyTemp_,&
  scalarBulkVolHeatCapVeg_,&
  computeEnthalpy,enthalpyStateVec,&
  scalarCanopyEnthTemp,scalarCanopyEnthalpy,&
  data_step,exposedVAI,&
  snowIncept,snowUnload,&
  airtemp,&
  refInterceptCapSnow,ratioDrip2Unloading,snowUnloadingCoeff,&
  minTempUnloading,minWindUnloading,rateTempUnloading,rateWindUnloading,&
  scalarNewSnowDensity,scalarCanairTemp,&
  scalarSnowfall,scalarWindspdCanopyTop,&
  scalarThroughfallSnow,scalarCanopySnowUnloading,&
  nSnow,nSoil,nLayers,&
  ix_canopySrad,vegTypeIndex,&
  scalarSWE,scalarSnowDepth,&
  mLayerVolFracLiq,spectralSnowAlbedoDiffuse,&
  scalarSnowAlbedo,mLayerTemp,&
  scalarGroundSnowFraction,scalarSnowAge,scalarCosZenith,&
  spectralSnowAlbedoDirect,scalarExposedLAI,scalarExposedSAI,&
  scalarCanopyWetFraction,scalarCanopySunlitFraction,&
  scalarCanopySunlitLAI,scalarCanopyShadedLAI,&
  spectralAlbGndDirect,spectralAlbGndDiffuse,scalarGroundAlbedo,&
  spectralIncomingDirect,spectralIncomingDiffuse,&
  scalarCanopySunlitPAR,scalarCanopyShadedPAR,&
  spectralBelowCanopyDirect,spectralBelowCanopyDiffuse,&
  scalarBelowCanopySolar,scalarCanopyAbsorbedSolar,scalarGroundAbsorbedSolar,&
  ALBSAT,ALBDRY,RHOL,RHOS,TAUL,TAUS,&
  alblak,omegas,betads,betais,hvt,hvb,rc,opt_rad,xl,opt_alb,&
  ixAlbedoMethod,&
  Frad_vis,Frad_direct,&
  albedoMax,albedoMinWinter,albedoMinSpring,&
  albedoMaxVisible,albedoMinVisible,&
  albedoMaxNearIR,albedoMinNearIR,&
  albedoDecayRate,tempScalGrowth,&
  albedoSootLoad,albedoRefresh,&
  scalarCanopyLiqMax,scalarCanopyIceMax,&
  canopyWettingFactor,canopyWettingExp,&
  computeVegFluxOld,modifiedVegState,&
  refInterceptCapRain,&
  bcUpprTdyn,bcUpprSoiH,&
  fracJulDay,yearLength,&
  latitude,&
  z0Snow,z0Soil,&
  heightCanopyTop,heightCanopyBottom,&
  scalarRootZoneTemp,scalarLAI,scalarSAI,&
  scalarGrowingSeasonIndex,&
  dveg,iswater,isbarren,issnow,&
  laim,saim,tmin,&
  urbanVegCategory,minExpLogHgtFac,&
  iLayerHeight,rootingDepth,scalarFoliageNitrogenFactor,oldSWE,&
  mLayerDepth,use_lookup,&
  ixNrgConserv,&
  scalarInitCanopyLiq,scalarInitCanopyIce,&
  scalarTotalSoilLiq,scalarTotalSoilIce, &
  mLayerVolFracIce, &
  balanceCanopyWater0, balanceSoilWater0, balanceAquifer0, &
  scalarAquiferStorage,flux_mean,&
  innerEffRainfall, sumCanopySublimation, sumLatHeatCanopyEvap, sumSenHeatCanopy, sumSnowSublimation,&
  innerBalance,innerBalanceLayerMass,innerBalanceLayerNrg,&
  mLayerVolFracIceInit,&
  sfcMeltPond, meanCanopySublimation,meanLatHeatCanopyEvap,meanSenHeatCanopy,&
  effRainfall,meanSoilCompress,meanBalance,modifiedLayers)
  USE tempAdjust_module,only:tempAdjust             ! adjust snow temperature associated with new snowfall
  use enthalpyTemp_module,only:T2enthTemp_veg
  USE canopySnow_module,only:canopySnow             ! compute interception and unloading of snow from the vegetation canopy
  USE vegSWavRad_module,only:vegSWavRad             ! compute canopy sw radiation fluxes
  USE snowAlbedo_module,only:snowAlbedo             ! compute snow albedo
  USE vegNrgFlux_module,only:wettedFrac             ! compute wetted fraction of the canopy (used in sw radiation fluxes)
  USE vegPhenlgy_module,only:vegPhenlgy_d             ! compute vegetation phenology
  use initialize_device,only:get_iGRU

  real(rkind),intent(inout) :: scalarThroughfallRain,scalarRainfall,scalarCanopyLiqDrainage
  logical(lgt),intent(inout) :: bal_veg,bal_snow,bal_soil,bal_aq
  logical(lgt),intent(inout) :: computeVegFlux
  real(rkind),intent(inout) :: canopyDepth_(:)
  real(rkind),intent(in) :: snowfrz_scale_(:),specificHeatVeg_(:),maxMassVegetation_(:)
  real(rkind),intent(inout) :: scalarCanopyLiq_(:),scalarCanopyIce_(:),scalarCanopyTemp_(:)
  real(rkind),intent(inout) :: scalarBulkVolHeatCapVeg_(:)
  logical(lgt),intent(inout) :: computeEnthalpy,enthalpyStateVec
  real(rkind),intent(inout) :: scalarCanopyEnthTemp,scalarCanopyEnthalpy
  real(rkind),intent(in) :: data_step
  real(rkind),intent(inout) :: exposedVAI
  integer(i4b),intent(in) :: snowIncept,snowUnload
  real(rkind),intent(inout) :: airtemp
  real(rkind),intent(inout) :: refInterceptCapSnow,ratioDrip2Unloading,snowUnloadingCoeff
  real(rkind),intent(inout) :: minTempUnloading,minWindUnloading,rateTempUnloading,rateWindUnloading
  real(rkind),intent(inout) :: scalarNewSnowDensity,scalarCanairTemp
  real(rkind),intent(inout) :: scalarSnowfall,scalarWindspdCanopyTop
  real(rkind),intent(inout) :: scalarThroughfallSnow,scalarCanopySnowUnloading
  integer(i4b),intent(inout) :: nSnow,nSoil,nLayers
  integer(i4b),intent(in) :: ix_canopySrad
  integer(i4b),intent(in) :: vegTypeIndex
  real(rkind),intent(inout) :: scalarSWE,scalarSnowDepth
  real(rkind),intent(inout) :: mLayerVolFracLiq(:), spectralSnowAlbedoDiffuse(:)
  real(rkind),intent(inout) :: scalarSnowAlbedo,mLayerTemp(:)
  real(rkind),intent(inout) :: scalarGroundSnowFraction,scalarSnowAge,scalarCosZenith
  real(rkind),intent(inout) :: spectralSnowAlbedoDirect(:),scalarExposedLAI,scalarExposedSAI
  real(rkind),intent(inout) :: scalarCanopyWetFraction,scalarCanopySunlitFraction
  real(rkind),intent(inout) :: scalarCanopySunlitLAI,scalarCanopyShadedLAI
  real(rkind),intent(inout) :: spectralAlbGndDirect(:), spectralAlbGndDiffuse(:),scalarGroundAlbedo
  real(rkind),intent(inout) :: spectralIncomingDirect(:),spectralIncomingDiffuse(:)
  real(rkind),intent(inout) :: scalarCanopySunlitPAR,scalarCanopyShadedPAR
  real(rkind),intent(inout) :: spectralBelowCanopyDirect(:), spectralBelowCanopyDiffuse(:)
  real(rkind),intent(inout) :: scalarBelowCanopySolar,scalarCanopyAbsorbedSolar,scalarGroundAbsorbedSolar
 real(rkind),intent(in) :: ALBSAT(:,:),ALBDRY(:,:),RHOL(:,:),RHOS(:,:),TAUL(:,:),TAUS(:,:)
   real(rkind),intent(in) :: alblak(:)
   real(rkind),intent(in) :: omegas(:), betads,betais
   real(rkind),intent(in) :: hvt(:),hvb(:),rc(:),xl(:)
   integer(i4b),intent(in) :: opt_alb,opt_rad
   integer(i4b),intent(in) :: ixAlbedoMethod
   real(rkind),intent(in) :: Frad_vis,Frad_direct
   real(rkind),intent(in) :: albedoMax,albedoMinWinter,albedoMinSpring
   real(rkind),intent(in) :: albedoMaxVisible,albedoMinVisible
   real(rkind),intent(in) :: albedoMaxNearIR,albedoMinNearIR
   real(rkind),intent(in) :: albedoDecayRate,tempScalGrowth
   real(rkind),intent(in) :: albedoSootLoad,albedoRefresh
   real(rkind),intent(inout) :: scalarCanopyLiqMax,scalarCanopyIceMax
   real(rkind),intent(in) :: canopyWettingFactor,canopyWettingExp
   logical(lgt),intent(inout) :: computeVegFluxOld,modifiedVegState
   real(rkind),intent(in) :: refInterceptCapRain
   integer(i4b),intent(in) :: bcUpprTdyn,bcUpprSoiH
   real(rkind),intent(in) :: fracJulDay
   integer(i4b),intent(in) :: yearLength
   real(rkind),intent(in) :: latitude
   real(rkind),intent(in) :: z0Snow,z0Soil
   real(rkind),intent(in) :: heightCanopyTop,heightCanopyBottom
   real(rkind),intent(inout) :: scalarRootZoneTemp,scalarLAI,scalarSAI
   real(rkind),intent(inout) :: scalarGrowingSeasonIndex
   integer(i4b),intent(in) :: dveg,iswater,isbarren,issnow
   real(rkind),intent(in) :: laim(:,:), saim(:,:), tmin(:)
   integer(i4b),intent(in) :: urbanVegCategory
   real(rkind),intent(in) :: minExpLogHgtFac
   real(rkind),intent(inout) :: iLayerHeight(0:)
   real(rkind),intent(inout) :: rootingDepth,scalarFoliageNitrogenFactor,oldSWE
   real(rkind),intent(inout) :: mLayerDepth(:)
   logical(lgt),intent(inout) :: use_lookup
   integer(i4b),intent(in) :: ixNrgConserv
   real(rkind),intent(inout) :: scalarInitCanopyLiq,scalarInitCanopyIce
   real(rkind),intent(inout) :: scalarTotalSoilLiq,scalarTotalSoilIce
   real(rkind),intent(inout) :: mLayerVolFracIce(:)
   real(rkind),intent(inout) :: balanceCanopyWater0, balanceSoilWater0, balanceAquifer0
   real(rkind),intent(inout) :: scalarAquiferStorage
   real(rkind),intent(inout) :: flux_mean(:)
   real(rkind),intent(inout) :: innerEffRainfall, sumCanopySublimation, sumLatHeatCanopyEvap, sumSenHeatCanopy, sumSnowSublimation
   real(rkind),intent(inout) :: innerBalance(:), innerBalanceLayerMass(:), innerBalanceLayerNrg(:)
   real(rkind),intent(inout) :: mLayerVolFracIceInit(:)
   real(rkind),intent(inout) :: sfcMeltPond, meanCanopySublimation,meanLatHeatCanopyEvap,meanSenHeatCanopy
   real(rkind),intent(inout) :: effRainfall
   real(rkind),intent(inout) :: meanSoilCompress(:), meanBalance(:)
   logical(lgt),intent(inout) :: modifiedLayers

   integer(i4b) :: nLayersRoots
  integer(i4b) :: err
  real(rkind)                          :: dCanopyWetFraction_dWat ! derivative in wetted fraction w.r.t. canopy total water (kg-1 m2)
  real(rkind)                          :: dCanopyWetFraction_dT   ! derivative in wetted fraction w.r.t. canopy temperature (K-1)
  real(rkind),parameter                :: varNotUsed1=-9999._rkind ! variables used to calculate derivatives (not needed here)
  real(rkind),parameter                :: varNotUsed2=-9999._rkind ! variables used to calculate derivatives (not needed here)
  integer(i4b) :: iGRU
  iGRU = get_iGRU()
  ! initialize variables
  innerEffRainfall=0._rkind       ! inner step average effective rainfall into snow (kg m-2 s-1) 
  sumCanopySublimation=0._rkind   ! sum of sublimation from the vegetation canopy (kg m-2 s-1) over substep
  sumLatHeatCanopyEvap=0._rkind   ! sum of latent heat flux for evaporation from the canopy to the canopy air space (W m-2) over substep
  sumSenHeatCanopy=0._rkind       ! sum of sensible heat flux from the canopy to the canopy air space (W m-2) over substep
  sumSnowSublimation=0._rkind     ! sum of sublimation from the snow surface (kg m-2 s-1) over substep
  innerBalance = 0._rkind         ! mean total balance array

  innerBalanceLayerMass = 0._rkind ! mean total balance of mass in layers
  innerBalanceLayerNrg = 0._rkind ! mean total balance of energy in layers
  mLayerVolFracIceInit = mLayerVolFracIce ! volume fraction of water ice
  ! initialize surface melt pond
  sfcMeltPond       = 0._rkind  ! change in storage associated with the surface melt pond (kg m-2)

  ! initialize average over data_step (averaged over substep in varSubStep)
  meanCanopySublimation = 0._rkind ! mean canopy sublimation
  meanLatHeatCanopyEvap = 0._rkind ! mean latent heat flux for evaporation from the canopy
  meanSenHeatCanopy     = 0._rkind ! mean sensible heat flux from the canopy
  effRainfall           = 0._rkind ! mean total effective rainfall over snow

  meanSoilCompress = 0._rkind ! mean total soil compression

  ! initialize the balance checks
  meanBalance = 0._rkind

  ! initialize flags to modify the veg layers or modify snow layers
  modifiedLayers    = .false.    ! flag to denote that snow layers were modified
  modifiedVegState  = .false.    ! flag to denote that vegetation states were modified

  ! initialize fluxes to average over data_step (averaged over substep in varSubStep)
  flux_mean = 0._rkind

        ! set the flag to compute enthalpy, may want to have this true always if want to output enthalpy
      computeEnthalpy  = .false.
      enthalpyStateVec = .false.
      use_lookup       = .false.
       if(ixNrgConserv.ne.closedForm) enthalpyStateVec = .true. ! enthalpy as state variable
      if(ixNrgConserv==enthalpyFormLU) use_lookup = .true. ! use lookup tables for soil temperature-enthalpy instead of analytical solution

      ! save the liquid water and ice on the vegetation canopy
      scalarInitCanopyLiq = scalarCanopyLiq_(iGRU)    ! initial liquid water on the vegetation canopy (kg m-2)
      scalarInitCanopyIce = scalarCanopyIce_(iGRU)    ! initial ice          on the vegetation canopy (kg m-2)

      ! compute total soil moisture and ice at the *START* of the step (kg m-2)
      scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(1:nSoil)*mLayerDepth(1:nSoil))
      scalarTotalSoilIce = sum(iden_water*mLayerVolFracIce(1:nSoil)*mLayerDepth(1:nSoil))  ! NOTE: no expansion and hence use iden_water

      ! compute storage of water in the canopy and the soil
      balanceCanopyWater0 = scalarCanopyLiq_(iGRU) + scalarCanopyIce_(iGRU)
      balanceSoilWater0   = scalarTotalSoilLiq + scalarTotalSoilIce

      ! get the total aquifer storage at the start of the time step (kg m-2)
      balanceAquifer0 = scalarAquiferStorage*iden_water
      ! save liquid water content
      ! if(printBalance)then
      !   allocate(liqSnowInit(nSnow), liqSoilInit(nSoil), stat=err)
      !   if(err/=0)then
      !     message=trim(message)//'unable to allocate space for the initial vectors'
      !     err=20; return
      !   endif
      !   if(nSnow>0) liqSnowInit = prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)
      !               liqSoilInit =                         mLayerVolFracLiq
      ! endif


    ! compute the number of layers with roots
    nLayersRoots = count(iLayerHeight(nSnow:nLayers-1) < rootingDepth-verySmall)
    if(nLayersRoots == 0)then
      ! message=trim(message)//'no roots within the soil profile'
      err=20; return
    end if

    ! define the foliage nitrogen factor
    scalarFoliageNitrogenFactor = 1._rkind  ! foliage nitrogen concentration (1.0 = saturated)

    ! save SWE
    oldSWE = scalarSWE

      ! *** compute phenology...
    ! ------------------------

    ! compute the temperature of the root zone: used in vegetation phenology
    scalarRootZoneTemp = sum(mLayerTemp(nSnow+1:nSnow+nLayersRoots)) / real(nLayersRoots, kind(rkind))

    ! remember if we compute the vegetation flux on the previous sub-step
    computeVegFluxOld = computeVegFlux

    ! compute the exposed LAI and SAI and whether veg is buried by snow
    call vegPhenlgy_d(&
                    ! model control
                    nSnow,                       & ! intent(in):    number of snow layers
                    ! model_decisions,             & ! intent(in):    model decisions
                    bcUpprTdyn,bcUpprSoiH,&
                    ! input/output: data structures
                    fracJulDay,                  & ! intent(in):    fractional julian days since the start of year
                    yearLength,                  & ! intent(in):    number of days in the current year
                    ! type_data,                   & ! intent(in):    type of vegetation and soil
                    vegTypeIndex, &
                    ! attr_data,                   & ! intent(in):    spatial attributes
                    latitude, &
                    ! mpar_data,                   & ! intent(in):    model parameters
                    z0Snow,z0Soil,&
                    heightCanopyTop,heightCanopyBottom,&
                    ! prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                    scalarSnowDepth,scalarCanopyTemp_(iGRU),scalarCanopyLiq_(iGRU), &
                    ! diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                    scalarRootZoneTemp,scalarLAI,scalarSAI,&
                    scalarExposedLAI,scalarExposedSAI,&
                    scalarGrowingSeasonIndex,scalarGroundSnowFraction,&
                    ! output
                    computeVegFlux,              & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                    canopyDepth_(iGRU),                 & ! intent(out): canopy depth (m)
                    exposedVAI,                  & ! intent(out): exposed vegetation area index (m2 m-2)
                    err,&
                    dveg,iswater,isbarren,issnow,&
                    hvt,hvb,laim,saim,tmin,urbanVegCategory,minExpLogHgtFac)                  ! intent(out): error control
    ! if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

    ! check
    if(computeVegFlux)then
      if(canopyDepth_(iGRU) < epsilon(canopyDepth_(iGRU)))then
        ! message=trim(message)//'canopy depth is zero when computeVegFlux flag is .true.'
        err=20; return
      endif
    endif


    ! flag the case where number of vegetation states has changed
    modifiedVegState = (computeVegFlux.neqv.computeVegFluxOld)

    ! *** compute wetted canopy area...
    ! ---------------------------------

    ! compute maximum canopy liquid water (kg m-2)
    scalarCanopyLiqMax = refInterceptCapRain*exposedVAI

    ! compute maximum canopy ice content (kg m-2)
    ! NOTE 1: this is used to compute the snow fraction on the canopy, as used in *BOTH* the radiation AND canopy sublimation routines
    ! NOTE 2: this is a different variable than the max ice used in the throughfall (snow interception) calculations
    ! NOTE 3: use maximum per unit leaf area storage capacity for snow (kg m-2)
    select case(snowIncept)
      case(lightSnow);  scalarCanopyIceMax = exposedVAI*refInterceptCapSnow
      case(stickySnow); scalarCanopyIceMax = exposedVAI*refInterceptCapSnow*4._rkind
      ! case default; message=trim(message)//'unable to identify option for maximum branch interception capacity'; err=20; return
    end select ! identifying option for maximum branch interception capacity

      ! compute wetted fraction of the canopy
    ! NOTE: assume that the wetted fraction is constant over the substep for the radiation calculations
    if(computeVegFlux)then

      ! compute wetted fraction of the canopy
      call wettedFrac(&
                      ! input
                      .false.,                                                      & ! flag to denote if derivatives are required
                      (scalarCanopyTemp_(iGRU) < Tfreeze), & ! flag to denote if the canopy is frozen
                      varNotUsed1,                                                  & ! derivative in canopy liquid w.r.t. canopy temperature (kg m-2 K-1)
                      varNotUsed2,                                                  & ! fraction of liquid water on the canopy
                      scalarCanopyLiq_(iGRU),              & ! canopy liquid water (kg m-2)
                      scalarCanopyIce_(iGRU),              & ! canopy ice (kg m-2)
                      scalarCanopyLiqMax,           & ! maximum canopy liquid water (kg m-2)
                      scalarCanopyIceMax,           & ! maximum canopy ice content (kg m-2)
                      canopyWettingFactor,         & ! maximum wetted fraction of the canopy (-)
                      canopyWettingExp,            & ! exponent in canopy wetting function (-)
                      ! output
                      scalarCanopyWetFraction,      & ! canopy wetted fraction (-)
                      dCanopyWetFraction_dWat,                                      & ! derivative in wetted fraction w.r.t. canopy liquid water content (kg-1 m2)
                      dCanopyWetFraction_dT,                                        & ! derivative in wetted fraction w.r.t. canopy liquid water content (kg-1 m2)
                      err)
      ! if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

    ! vegetation is completely buried by snow (or no veg exists at all)
    else
      scalarCanopyWetFraction = 0._rkind
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
                    (nSnow > 0),                 & ! intent(in): logical flag to denote if snow is present
                    ! input/output: data structures
                    ! model_decisions,             & ! intent(in):    model decisions
                    ix_canopySrad,ixAlbedoMethod,&
                    ! mpar_data,                   & ! intent(in):    model parameters
                    Frad_vis,Frad_direct,&
                    albedoMax,albedoMinWinter,albedoMinSpring,&
                    albedoMaxVisible,albedoMinVisible,&
                    albedoMaxNearIR,albedoMinNearIR,&
                    albedoDecayRate,tempScalGrowth,&
                    albedoSootLoad,albedoRefresh,snowfrz_scale_(iGRU),&
                    ! flux_data,                   & ! intent(in):    model flux variables
                    scalarSnowfall, &
                    ! diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                    scalarCosZenith,spectralSnowAlbedoDirect,&
                    ! prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                    mLayerTemp(1),scalarSnowAlbedo,spectralSnowAlbedoDiffuse, &
                    ! output: error control
                    err)                  ! intent(out): error control
    ! if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if


    ! *** compute canopy sw radiation fluxes...
    ! -----------------------------------------
    call vegSWavRad(&
                    data_step,                    & ! intent(in):    time step (s) -- only used in Noah-MP radiation, to compute albedo
                    nSnow,                        & ! intent(in):    number of snow layers
                    nSoil,                        & ! intent(in):    number of soil layers
                    nLayers,                      & ! intent(in):    total number of layers
                    computeVegFlux,               & ! intent(in):    logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                       ix_canopySrad, &
                    !    type_data,                    & ! intent(in):    classification of veg, soil etc. for a local HRU
                       vegTypeIndex, &
                    !    prog_data,                    & ! intent(inout): model prognostic variables for a local HRU
                       scalarSWE,scalarSnowDepth,&
                       mLayerVolFracLiq,spectralSnowAlbedoDiffuse,&
                       scalarSnowAlbedo,&
                       mLayerTemp(1),scalarCanopyTemp_(iGRU), &
                    !    diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                       scalarGroundSnowFraction,&
                       scalarSnowAge,scalarCosZenith,spectralSnowAlbedoDirect,&
                       scalarExposedLAI,scalarExposedSAI,scalarCanopyWetFraction, &
                       scalarCanopySunlitFraction,scalarCanopySunlitLAI,scalarCanopyShadedLAI,&
                       spectralAlbGndDirect,spectralAlbGndDiffuse,scalarGroundAlbedo,&
                    !    flux_data,                    & ! intent(inout): model flux variables
                       scalarSnowfall,spectralIncomingDirect,spectralIncomingDiffuse, &
                       scalarCanopySunlitPAR,scalarCanopyShadedPAR,&
                       spectralBelowCanopyDirect,spectralBelowCanopyDiffuse,&
                       scalarBelowCanopySolar,scalarCanopyAbsorbedSolar,scalarGroundAbsorbedSolar,&
                       err,ALBSAT,ALBDRY,RHOL,RHOS,TAUL,TAUS,&
                       alblak,omegas,betads,betais,hvt,hvb,rc,opt_rad,xl,opt_alb)                    ! intent(out): error control
    ! if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if


      ! *** compute canopy throughfall and unloading...
    ! -----------------------------------------------
    ! NOTE 1: this needs to be done before solving the energy and liquid water equations, to account for the heat advected with precipitation (and throughfall/unloading)
    ! NOTE 2: the unloading flux is computed using canopy drip (scalarCanopyLiqDrainage) from the previous time step
    ! this changes canopy ice
    call canopySnow(&
                    ! input: model control
                    data_step,                   & ! intent(in): time step (seconds)
                    exposedVAI,                  & ! intent(in): exposed vegetation area index (m2 m-2)
                    computeVegFlux,              & ! intent(in): flag to denote if computing energy flux over vegetation
                    ! input/output: data structures
                    ! model_decisions,             & ! intent(in):    model decisions
                    snowIncept, &
                    snowUnload, &
                    ! forc_data,                   & ! intent(in):    model forcing data
                    airtemp, &
                    ! mpar_data,                   & ! intent(in):    model parameters
                    refInterceptCapSnow, &
                    ratioDrip2Unloading, &
                    snowUnloadingCoeff, &
                    minTempUnloading, &
                    minWindUnloading, &
                    rateTempUnloading, &
                    rateWindUnloading, &
                    ! diag_data,                   & ! intent(in):    model diagnostic variables for a local HRU &
                    scalarNewSnowDensity, &
                    ! prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                    scalarCanopyIce_(iGRU), &
                    scalarCanairTemp, &
                    ! flux_data,                   & ! intent(inout): model flux variables
                    scalarSnowfall, &
                    scalarCanopyLiqDrainage, &
                    scalarWindspdCanopyTop, &
                    scalarThroughfallSnow, &
                    scalarCanopySnowUnloading, &
                    ! output: error control
                    err)                  ! intent(out): error control
    ! if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if


      ! adjust canopy temperature to account for new snow
    if(computeVegFlux)then ! logical flag to compute vegetation fluxes (.false. if veg buried by snow)
      call tempAdjust(&
                      ! input: derived parameters
                      canopyDepth_,                 & ! intent(in):    canopy depth (m)
                      ! input/output: data structures
                      ! mpar_data,                   & ! intent(in):    model parameters
                      snowfrz_scale_,specificHeatVeg_,maxMassVegetation_,&
                      ! prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                      scalarCanopyLiq_,scalarCanopyIce_,scalarCanopyTemp_,&
                      ! diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                      scalarBulkVolHeatCapVeg_,&
                      ! output: error control
                      err)                  ! intent(out): error control
      ! if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

      ! change enthalpy based on new canopy temperature and water content, only if will need enthalpy for energy balance
      if(enthalpyStateVec .or. computeEnthalpy)then
        ! associate local variables with variables in the data structures
        ! enthalpyVeg: associate(&
        !   ! state variables in the vegetation canopy
        !   scalarCanopyTemp     => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)     ,& ! canopy temperature (K)
        !   scalarCanopyEnthTemp => diag_data%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1) ,& ! canopy temperature component of enthalpy (J m-3)
        !   scalarCanopyEnthalpy => prog_data%var(iLookPROG%scalarCanopyEnthalpy)%dat(1) ,& ! enthalpy of the vegetation canopy (J m-3)
        !   scalarCanopyLiq      => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)      ,& ! mass of liquid water on the vegetation canopy (kg m-2)
        !   scalarCanopyIce      => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)       & ! mass of ice on the vegetation canopy (kg m-2)
        !   )  ! (associate local variables with model parameters)       
          call T2enthTemp_veg(&
                          canopyDepth_(iGRU),            & ! intent(in): canopy depth (m)
                          specificHeatVeg_(iGRU),        & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                          maxMassVegetation_(iGRU),      & ! intent(in): maximum mass of vegetation (kg m-2)
                          snowfrz_scale_(iGRU),          & ! intent(in): scaling parameter for the snow freezing curve  (K-1)
                          scalarCanopyTemp_(iGRU),       & ! intent(in): canopy temperature (K)
                          (scalarCanopyLiq_(iGRU)+scalarCanopyIce_(iGRU)), & ! intent(in): canopy water content (kg m-2)
                          scalarCanopyEnthTemp)     ! intent(out): temperature component of enthalpy of the vegetation canopy (J m-3)
          scalarCanopyEnthalpy = scalarCanopyEnthTemp  - LH_fus * scalarCanopyIce_(iGRU)/ canopyDepth_(iGRU) ! new ice and/or temperature
        ! end associate enthalpyVeg
      end if ! (need to recalculate enthalpy state variable)
    end if ! if computing fluxes over vegetation

    ! initialize drainage and throughfall
    ! NOTE 1: this needs to be done before solving the energy and liquid water equations, to account for the heat advected with precipitation
    ! NOTE 2: this initialization needs to be done AFTER the call to canopySnow, since canopySnow uses canopy drip drom the previous time step
    if(.not.computeVegFlux)then
      scalarThroughfallRain = scalarRainfall
    else
      scalarThroughfallRain = 0._rkind
    end if
    scalarCanopyLiqDrainage = 0._rkind

    ! initialize if used a balance
    bal_veg = .false.
    bal_snow = .false.
    bal_soil = .false.
    bal_aq   = .false.

end subroutine
attributes(global) subroutine finalize_coupled_em(nGRU,scalarSnowfallTemp, scalarNewSnowDensity, &
  flux_data,scalarThroughfallSnow,scalarCanopySnowUnloading,&
  scalarSWE,scalarSnowDepth,&
  mLayerTemp,mLayerDepth,&
  mLayerVolFracIce,mLayerVolFracLiq,&
  data_step,nSnow,snowfrz_scale,&
  mLayerVolFracWat,mLayerEnthTemp,mLayerEnthalpy,&
  computeEnthalpy,enthalpyStateVec,&
  nLayers,layerType,nSoil,&
  mLayerHeight,iLayerHeight,&
  flux_mean,numScalarFluxData,&
  mLayerCompress,meanSoilCompress,scalarSoilCompress,&
  scalarSfcMeltPond,sfcMeltPond,&
  meanBalance,&
  balanceCasNrg,balanceVegNrg,balanceVegMass,balanceAqMass,&
  balanceSnowNrg,balanceSoilNrg,balanceSnowMass,balanceSoilMass, &
  bal_veg,bal_soil,bal_snow,bal_aq,&
  computeVegFlux,&
  scalarCanopyWatBalError,scalarCanopyWat,balanceCanopyWater0,&
  scalarSnowfall,averageThroughfallSnow,scalarRainfall,averageThroughfallRain,&
  averageCanopySnowUnloading,averageCanopyLiqDrainage,averageCanopySublimation,averageCanopyEvaporation,&
  effSnowfall,delSWE,oldSWE,&
  massBalance,effRainfall,&
  averageSnowSublimation,averageSnowDrainage,&
  scalarTotalSoilLiq,scalarTotalSoilIce,scalarTotalSoilWat,&
  balanceAquifer1,scalarAquiferStorage,&
  balanceSoilInflux,averageSoilInflux,&
  balanceSoilBaseflow,averageSoilBaseflow,&
  balanceSoilDrainage,averageSoilDrainage, &
  balanceSoilET,averageCanopyTranspiration,averageGroundEvaporation, &
  balanceSoilCompress, &
  scalarSoilWatBalError,balanceSoilWater0,&
  scalarCanopyEnthalpy,scalarCanopyEnthTemp,canopyDepth,&
  ixDomainType,ixControlVolume,ixStateType,&
  scalarCanopyIce,scalarTotalSoilEnthalpy,scalarTotalSnowEnthalpy,&
  scalarSurfaceTemp)
  integer(i4b),intent(in),value :: nGRU
  real(rkind),intent(inout) :: scalarSnowfallTemp(:), scalarNewSnowDensity(:)
  real(rkind),intent(inout) :: flux_data(:,:)
  integer(i4b),intent(in),value :: scalarThroughfallSnow,scalarCanopySnowUnloading
  real(rkind),intent(inout) :: scalarSWE(:), scalarSnowDepth(:)
  real(rkind),intent(inout) :: mLayerTemp(:,:), mLayerDepth(:,:)
  real(rkind),intent(inout) :: mLayerVolFracIce(:,:), mLayerVolFracLiq(:,:)
  real(rkind),intent(in),value :: data_step
  integer(i4b),intent(inout) :: nSnow(:)
  real(rkind),intent(in) :: snowfrz_scale(:)
  real(rkind),intent(inout) :: mLayerVolFracWat(:,:),mLayerEnthTemp(:,:),mLayerEnthalpy(:,:)
  logical(lgt),intent(in) :: computeEnthalpy,enthalpyStateVec
  integer(i4b),intent(inout) :: nLayers(:),layerType(:,:)
  integer(i4b),intent(in),value :: nSoil
  real(rkind),intent(inout) :: mLayerHeight(0:,:), iLayerHeight(0:,:)
  real(rkind),intent(inout) :: flux_mean(:,:)
  integer(i4b),intent(in),value :: numScalarFluxData
  real(rkind),intent(inout) :: mLayerCompress(:,:), meanSoilCompress(:,:), scalarSoilCompress(:)
real(rkind),intent(inout) :: scalarSfcMeltPond(:),sfcMeltPond(:)
  real(rkind),intent(inout) :: meanBalance(:,:)
  real(rkind),intent(inout) :: balanceCasNrg(:),balanceVegNrg(:),balanceVegMass(:),balanceAqMass(:)
  real(rkind),intent(inout) :: balanceSnowNrg(:),balanceSoilNrg(:),balanceSnowMass(:),balanceSoilMass(:)
  logical(lgt),intent(in) :: bal_veg(:),bal_soil(:),bal_snow(:),bal_aq(:)
  logical(lgt),intent(in) :: computeVegFlux(:)
  real(rkind),intent(inout) :: scalarCanopyWatBalError(:),scalarCanopyWat(:),balanceCanopyWater0(:)
  integer(i4b),intent(in),value :: scalarSnowfall,averageThroughfallSnow,scalarRainfall,averageThroughfallRain
  integer(i4b),intent(in),value :: averageCanopySnowUnloading,averageCanopyLiqDrainage,averageCanopySublimation,averageCanopyEvaporation
  real(rkind),intent(inout) :: effSnowfall(:),delSWE(:),oldSWE(:)
  real(rkind),intent(inout) :: massBalance(:),effRainfall(:)
  integer(i4b),intent(in),value :: averageSnowSublimation,averageSnowDrainage
  real(rkind),intent(inout) :: scalarTotalSoilLiq(:),scalarTotalSoilIce(:),scalarTotalSoilWat(:)
  real(rkind),intent(inout) :: balanceAquifer1(:),scalarAquiferStorage(:)
  real(rkind),intent(inout) :: balanceSoilInflux(:)
  integer(i4b),intent(in),value :: averageSoilInflux
  real(rkind),intent(inout) :: balanceSoilBaseflow(:)
  integer(i4b),intent(in),value :: averageSoilBaseflow
  real(rkind),intent(inout) :: balanceSoilDrainage(:)
  integer(i4b),intent(in),value :: averageSoilDrainage
  real(rkind),intent(inout) :: balanceSoilET(:)
  integer(i4b),intent(in),value :: averageCanopyTranspiration,averageGroundEvaporation
  real(rkind),intent(inout) :: balanceSoilCompress(:)
  real(rkind),intent(inout) :: scalarSoilWatBalError(:),balanceSoilWater0(:)
  real(rkind),intent(inout) :: scalarCanopyEnthalpy(:),scalarCanopyEnthTemp(:)
  real(rkind),intent(inout) :: canopyDepth(:)
  integer(i4b),intent(in) :: ixDomainType(:,:), ixControlVolume(:,:), ixStateType(:,:)
  real(rkind),intent(inout) :: scalarCanopyIce(:),scalarTotalSoilEnthalpy(:),scalarTotalSnowEnthalpy(:)
  real(rkind),intent(inout) :: scalarSurfaceTemp(:)

         integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return
  call finalize_coupled_em_device(scalarSnowfallTemp(iGRU), scalarNewSnowDensity(iGRU), &
  flux_data(scalarThroughfallSnow,iGRU),flux_data(scalarCanopySnowUnloading,iGRU),&
  scalarSWE(iGRU),scalarSnowDepth(iGRU),&
  mLayerTemp(:,iGRU),mLayerDepth(:,iGRU),&
  mLayerVolFracIce(:,iGRU),mLayerVolFracLiq(:,iGRU),&
  data_step,nSnow(iGRU),snowfrz_scale(iGRU),&
  mLayerVolFracWat(:,iGRU),mLayerEnthTemp(:,iGRU),mLayerEnthalpy(:,iGRU),&
  computeEnthalpy,enthalpyStateVec,&
  nLayers(iGRU),layerType(:,iGRU),nSoil,&
  mLayerHeight(:,iGRU),iLayerHeight(:,iGRU),&
  flux_data(:,iGRU),flux_mean(:,iGRU),numScalarFluxData,&
  mLayerCompress(:,iGRU),meanSoilCompress(:,iGRU),scalarSoilCompress(iGRU),&
  scalarSfcMeltPond(iGRU),sfcMeltPond(iGRU),&
  meanBalance(:,iGRU),&
  balanceCasNrg(iGRU),balanceVegNrg(iGRU),balanceVegMass(iGRU),balanceAqMass(iGRU),&
  balanceSnowNrg(iGRU),balanceSoilNrg(iGRU),balanceSnowMass(iGRU),balanceSoilMass(iGRU), &
  bal_veg(iGRU),bal_soil(iGRU),bal_snow(iGRU),bal_aq(iGRU),&
  computeVegFlux(iGRU),&
  scalarCanopyWatBalError(iGRU),scalarCanopyWat(iGRU),balanceCanopyWater0(iGRU),&
  flux_mean(scalarSnowfall,iGRU),flux_mean(averageThroughfallSnow,iGRU),flux_mean(scalarRainfall,iGRU),flux_mean(averageThroughfallRain,iGRU),&
  flux_mean(averageCanopySnowUnloading,iGRU),flux_mean(averageCanopyLiqDrainage,iGRU),flux_mean(averageCanopySublimation,iGRU),flux_mean(averageCanopyEvaporation,iGRU),&
  effSnowfall(iGRU),delSWE(iGRU),oldSWE(iGRU),&
  massBalance(iGRU),effRainfall(iGRU),&
  flux_mean(averageSnowSublimation,iGRU),flux_mean(averageSnowDrainage,iGRU),&
  scalarTotalSoilLiq(iGRU),scalarTotalSoilIce(iGRU),scalarTotalSoilWat(iGRU),&
  balanceAquifer1(iGRU),scalarAquiferStorage(iGRU),&
  balanceSoilInflux(iGRU),flux_mean(averageSoilInflux,iGRU),&
  balanceSoilBaseflow(iGRU),flux_mean(averageSoilBaseflow,iGRU),&
  balanceSoilDrainage(iGRU),flux_mean(averageSoilDrainage,iGRU), &
  balanceSoilET(iGRU),flux_mean(averageCanopyTranspiration,iGRU),flux_mean(averageGroundEvaporation,iGRU), &
  balanceSoilCompress(iGRU), &
  scalarSoilWatBalError(iGRU),balanceSoilWater0(iGRU),&
  scalarCanopyEnthalpy(iGRU),scalarCanopyEnthTemp(iGRU),canopyDepth(iGRU),&
  ixDomainType(:,iGRU),ixControlVolume(:,iGRU),ixStateType(:,iGRU),&
  scalarCanopyIce(iGRU),scalarTotalSoilEnthalpy(iGRU),scalarTotalSnowEnthalpy(iGRU),&
  scalarSurfaceTemp(iGRU))
end subroutine

attributes(device) subroutine finalize_coupled_em_device(scalarSnowfallTemp, scalarNewSnowDensity, &
  scalarThroughfallSnow,scalarCanopySnowUnloading,&
  scalarSWE,scalarSnowDepth,&
  mLayerTemp,mLayerDepth,&
  mLayerVolFracIce,mLayerVolFracLiq,&
  data_step,nSnow,snowfrz_scale,&
  mLayerVolFracWat,mLayerEnthTemp,mLayerEnthalpy,&
  computeEnthalpy,enthalpyStateVec,&
  nLayers,layerType,nSoil,&
  mLayerHeight,iLayerHeight,&
  flux_data,flux_mean,numScalarFluxData,&
  mLayerCompress,meanSoilCompress,scalarSoilCompress,&
  scalarSfcMeltPond,sfcMeltPond,&
  meanBalance,&
  balanceCasNrg,balanceVegNrg,balanceVegMass,balanceAqMass,&
  balanceSnowNrg,balanceSoilNrg,balanceSnowMass,balanceSoilMass, &
  bal_veg,bal_soil,bal_snow,bal_aq,&
  computeVegFlux,&
  scalarCanopyWatBalError,scalarCanopyWat,balanceCanopyWater0,&
  scalarSnowfall,averageThroughfallSnow,scalarRainfall,averageThroughfallRain,&
  averageCanopySnowUnloading,averageCanopyLiqDrainage,averageCanopySublimation,averageCanopyEvaporation,&
  effSnowfall,delSWE,oldSWE,&
  massBalance,effRainfall,&
  averageSnowSublimation,averageSnowDrainage,&
  scalarTotalSoilLiq,scalarTotalSoilIce,scalarTotalSoilWat,&
  balanceAquifer1,scalarAquiferStorage,&
  balanceSoilInflux,averageSoilInflux,&
  balanceSoilBaseflow,averageSoilBaseflow,&
  balanceSoilDrainage,averageSoilDrainage, &
  balanceSoilET,averageCanopyTranspiration,averageGroundEvaporation, &
  balanceSoilCompress, &
  scalarSoilWatBalError,balanceSoilWater0,&
  scalarCanopyEnthalpy,scalarCanopyEnthTemp,canopyDepth,&
  ixDomainType,ixControlVolume,ixStateType,&
  scalarCanopyIce,scalarTotalSoilEnthalpy,scalarTotalSnowEnthalpy,&
  scalarSurfaceTemp)
  USE volicePack_module,only:newsnwfall             ! compute change in the top snow layer due to throughfall and unloading
  USE enthalpyTemp_module,only:T2enthTemp_snow      ! convert temperature to enthalpy for snow
  USE var_derive_module,only:calcHeight_d             ! module to calculate height at layer interfaces and layer mid-point
  use enthalpyTemp_module,only:enthTemp_or_enthalpy

  real(rkind),intent(inout) :: scalarSnowfallTemp, scalarNewSnowDensity
  real(rkind),intent(inout) :: scalarThroughfallSnow,scalarCanopySnowUnloading
  real(rkind),intent(inout) :: scalarSWE, scalarSnowDepth
  real(rkind),intent(inout) :: mLayerTemp(:), mLayerDepth(:)
  real(rkind),intent(inout) :: mLayerVolFracIce(:), mLayerVolFracLiq(:)
  real(rkind),intent(in) :: data_step
  integer(i4b),intent(inout) :: nSnow
  real(rkind),intent(in) :: snowfrz_scale
  real(rkind),intent(inout) :: mLayerVolFracWat(:), mLayerEnthTemp(:), mLayerEnthalpy(:)
  logical(lgt),intent(in) :: computeEnthalpy,enthalpyStateVec
  integer(i4b),intent(inout) :: nLayers,layerType(:)
  integer(i4b),intent(in) :: nSoil
  real(rkind),intent(inout) :: mLayerHeight(0:), iLayerHeight(0:)
  real(rkind),intent(inout) :: flux_data(:), flux_mean(:)
  integer(i4b),intent(in) :: numScalarFluxData
  real(rkind),intent(inout) :: mLayerCompress(:), meanSoilCompress(:), scalarSoilCompress
  real(rkind),intent(inout) :: scalarSfcMeltPond,sfcMeltPond
  real(rkind),intent(inout) :: meanBalance(:)
  real(rkind),intent(inout) :: balanceCasNrg,balanceVegNrg,balanceVegMass,balanceAqMass
  real(rkind),intent(inout) :: balanceSnowNrg,balanceSoilNrg,balanceSnowMass,balanceSoilMass
  logical(lgt),intent(in) :: bal_veg,bal_soil,bal_snow,bal_aq
  logical(lgt),intent(in) :: computeVegFlux
  real(rkind),intent(inout) :: scalarCanopyWatBalError,scalarCanopyWat,balanceCanopyWater0
  real(rkind),intent(inout) :: scalarSnowfall,averageThroughfallSnow,scalarRainfall,averageThroughfallRain
  real(rkind),intent(inout) :: averageCanopySnowUnloading,averageCanopyLiqDrainage,averageCanopySublimation,averageCanopyEvaporation
  real(rkind),intent(inout) :: effSnowfall,delSWE,oldSWE
  real(rkind),intent(inout) :: massBalance,effRainfall
  real(rkind),intent(inout) :: averageSnowSublimation,averageSnowDrainage
  real(rkind),intent(inout) :: scalarTotalSoilLiq,scalarTotalSoilIce,scalarTotalSoilWat
  real(rkind),intent(inout) :: balanceAquifer1,scalarAquiferStorage
  real(rkind),intent(inout) :: balanceSoilInflux,averageSoilInflux
  real(rkind),intent(inout) :: balanceSoilBaseflow,averageSoilBaseflow
  real(rkind),intent(inout) :: balanceSoilDrainage,averageSoilDrainage
  real(rkind),intent(inout) :: balanceSoilET,averageCanopyTranspiration,averageGroundEvaporation
  real(rkind),intent(inout) :: balanceSoilCompress
  real(rkind),intent(inout) :: scalarSoilWatBalError,balanceSoilWater0
  real(rkind),intent(inout) :: scalarCanopyEnthalpy,scalarCanopyEnthTemp
  real(rkind),intent(inout) :: canopyDepth
  integer(i4b),intent(in) :: ixDomainType(:), ixControlVolume(:), ixStateType(:)
  real(rkind),intent(inout) :: scalarCanopyIce,scalarTotalSoilEnthalpy,scalarTotalSnowEnthalpy
  real(rkind),intent(inout) :: scalarSurfaceTemp

  integer(i4b) :: err
  integer(i4b) :: iVar
  ! *** add snowfall to the snowpack...
  ! -----------------------------------
  ! add new snowfall to the snowpack
  ! NOTE: This needs to be done AFTER the call to canopySnow, since throughfall and unloading are computed in canopySnow
  call newsnwfall(&
                ! input: model control
                data_step,                                                 & ! time step (seconds)
                (nSnow > 0),                                               & ! logical flag if snow layers exist
                snowfrz_scale,                                             & ! freeezing curve parameter for snow (K-1)
                ! input: diagnostic scalar variables
                scalarSnowfallTemp,        & ! computed temperature of fresh snow (K)
                scalarNewSnowDensity,      & ! computed density of new snow (kg m-3)
                scalarThroughfallSnow,     & ! throughfall of snow through the canopy (kg m-2 s-1)
                scalarCanopySnowUnloading, & ! unloading of snow from the canopy (kg m-2 s-1)
                ! input/output: state variables
                scalarSWE,                 & ! SWE (kg m-2)
                scalarSnowDepth,           & ! total snow depth (m)
                mLayerTemp(1),                & ! temperature of the top layer (K)
                mLayerDepth(1),               & ! depth of the top layer (m)
                mLayerVolFracIce(1),          & ! volumetric fraction of ice of the top layer (-)
                mLayerVolFracLiq(1),          & ! volumetric fraction of liquid water of the top layer (-)
                ! output: error control
                err)                                                ! error control
  ! if(err/=0)then; err=30; message=trim(message)//trim(cmessage); return; end if

  ! recompute snow depth, SWE, and top layer water
  if(nSnow > 0)then
    scalarSnowDepth = sum(  mLayerDepth(1:nSnow))
    scalarSWE       = sum( (mLayerVolFracLiq(1:nSnow)*iden_water + &
                                                            mLayerVolFracIce(1:nSnow)*iden_ice) &
                                                          * mLayerDepth(1:nSnow) )
    mLayerVolFracWat(1) = mLayerVolFracLiq(1) &
                                                      + mLayerVolFracIce(1)*iden_ice/iden_water
    if(enthalpyStateVec .or. computeEnthalpy)then ! compute enthalpy of the top snow layer
      call T2enthTemp_snow(&
                     snowfrz_scale,                                     & ! intent(in):  scaling parameter for the snow freezing curve  (K-1)
                     mLayerTemp(1),        & ! temperature of the top layer (K)
                     mLayerVolFracWat(1),  & ! intent(in):  volumetric total water content (-)
                     mLayerEnthTemp(1))      ! intent(out): temperature component of enthalpy of each snow layer (J m-3)
      mLayerEnthalpy(1) = mLayerEnthTemp(1) - iden_ice * LH_fus * mLayerVolFracIce(1)
    end if
  end if

    ! re-assign dimension lengths
    nSnow   = count(layerType==iname_snow)
    ! nSoil   = count(layerType==iname_soil)
    nLayers = nSnow+nSoil

    ! update coordinate variables
    call calcHeight_d(&
    nLayers,layerType,&
    mLayerDepth,mLayerHeight,iLayerHeight)

  ! overwrite flux_data and soil compression with the timestep-average value (returns timestep-average fluxes for scalar variables)
  do iVar=1,numScalarFluxData
    flux_data(iVar) = flux_mean(iVar)
  end do
  ! keep soil compression as an average like the fluxes, will not want to do this if nSoil can change
  mLayerCompress = meanSoilCompress
  scalarSoilCompress = sum(  meanSoilCompress(1:nSoil)*iden_water &
                                                             * mLayerDepth(nSnow+1:nLayers) )
  
  ! ***********************************************************************************************************************************
  ! ---
  ! *** balance checks and summary variable saving...
  ! ---------------------

  ! save the average melt pond storage in the data structures
  scalarSfcMeltPond  = sfcMeltPond

  ! save balance of energy and water per single layer domain
  balanceCasNrg   = meanBalance(1) ! W m-3
  balanceVegNrg   = meanBalance(2) ! W m-3      will be realMissing if computeVegFlux is false
  balanceVegMass  = meanBalance(3) ! kg m-3 s-1 will be realMissing if computeVegFlux is false
  balanceAqMass   = meanBalance(4) ! kg m-2 s-1 will be realMissing if no aquifer
  balanceSnowNrg  = meanBalance(5) ! W m-3      will be realMissing if no snow during data step
  balanceSoilNrg  = meanBalance(6) ! W m-3       
  balanceSnowMass = meanBalance(7) ! kg m-3 s-1 will be realMissing if no snow during data step
  balanceSoilMass = meanBalance(8) ! kg m-3 s-1
  if(.not.bal_veg)then ! will be 0, make realMissing
    balanceCasNrg   = realMissing
    balanceVegNrg   = realMissing
    balanceVegMass  = realMissing
  endif
  if (.not.bal_snow)then ! will be 0, make realMissing
    balanceSnowNrg  = realMissing
    balanceSnowMass = realMissing
  endif
  if (.not.bal_soil)then ! will be 0, make realMissing
    balanceSoilNrg  = realMissing
    balanceSoilMass = realMissing
  endif
  if (.not.bal_aq)then ! will be 0, make realMissing
    balanceAqMass   = realMissing
  endif

      ! -----
      ! * balance checks for the canopy...
      ! ----------------------------------

      ! if computing the vegetation flux
      if(computeVegFlux)then
        ! balance checks for the canopy
        ! NOTE: need to put the balance checks in the sub-step loop so that we can recompute if necessary
        scalarCanopyWatBalError = scalarCanopyWat - (balanceCanopyWater0 + (scalarSnowfall - averageThroughfallSnow)*data_step + (scalarRainfall - averageThroughfallRain)*data_step &
                                  - averageCanopySnowUnloading*data_step - averageCanopyLiqDrainage*data_step + averageCanopySublimation*data_step + averageCanopyEvaporation*data_step)
        ! if(abs(scalarCanopyWatBalError) > absConvTol_liquid*iden_water*10._rkind .and. checkMassBalance_ds)then
          ! write(*,'(a,1x,f20.10)') 'data_step                    = ', data_step
          ! write(*,'(a,1x,f20.10)') 'balanceCanopyWater0          = ', balanceCanopyWater0
          ! write(*,'(a,1x,f20.10)') 'balanceCanopyWater1          = ', scalarCanopyWat
          ! write(*,'(a,1x,f20.10)') 'snowfall                     = ', scalarSnowfall*data_step
          ! write(*,'(a,1x,f20.10)') 'rainfall                     = ', scalarRainfall*data_step
          ! write(*,'(a,1x,f20.10)') '(snowfall - throughfallSnow) = ', (scalarSnowfall - averageThroughfallSnow)*data_step
          ! write(*,'(a,1x,f20.10)') '(rainfall - throughfallRain) = ', (scalarRainfall - averageThroughfallRain)*data_step
          ! write(*,'(a,1x,f20.10)') 'canopySnowUnloading          = ', averageCanopySnowUnloading*data_step
          ! write(*,'(a,1x,f20.10)') 'canopyLiqDrainage            = ', averageCanopyLiqDrainage*data_step
          ! write(*,'(a,1x,f20.10)') 'canopySublimation            = ', averageCanopySublimation*data_step
          ! write(*,'(a,1x,f20.10)') 'canopyEvaporation            = ', averageCanopyEvaporation*data_step
          ! write(*,'(a,1x,f20.10)') 'canopyWatBalError            = ', scalarCanopyWatBalError
          ! message=trim(message)//'canopy hydrology does not balance'
          ! err=20; return
        ! end if
      endif  ! if computing the vegetation flux

      ! -----
      ! * balance checks for SWE...
      ! ---------------------------

      ! check the individual layers
      ! if(printBalance .and. nSnow>0)then
        ! write(*,'(a,1x,10(f12.8,1x))') 'liqSnowInit       = ', liqSnowInit
        ! write(*,'(a,1x,10(f12.8,1x))') 'volFracLiq        = ', mLayerVolFracLiq(1:nSnow)
        ! write(*,'(a,1x,10(f12.8,1x))') 'iLayerLiqFluxSnow = ', flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat*iden_water*data_step
        ! write(*,'(a,1x,10(f12.8,1x))') 'mLayerLiqFluxSnow = ', flux_data%var(iLookFLUX%mLayerLiqFluxSnow)%dat*data_step
        ! write(*,'(a,1x,10(f12.8,1x))') 'change volFracLiq = ', mLayerVolFracLiq(1:nSnow) - liqSnowInit
        ! deallocate(liqSnowInit, stat=err)
        ! if(err/=0)then
        !   message=trim(message)//'unable to deallocate space for the initial volumetric liquid water content of snow'
        !   err=20; return
        ! endif
      ! endif

      ! check SWE
      if(nSnow>0)then
        effSnowfall = averageThroughfallSnow + averageCanopySnowUnloading
        ! effRainfall is averageThroughfallRain + averageCanopyLiqDrainage only over snow
        delSWE      = scalarSWE - (oldSWE - sfcMeltPond)
        massBalance = delSWE - (effSnowfall + effRainfall + averageSnowSublimation - averageSnowDrainage*iden_water)*data_step
        ! if(abs(massBalance) > absConvTol_liquid*iden_water*10._rkind .and. checkMassBalance_ds)then
          ! print*,                  'nSnow       = ', nSnow
          ! print*,                  'nSub        = ', nSub
          ! write(*,'(a,1x,f20.10)') 'data_step   = ', data_step
          ! write(*,'(a,1x,f20.10)') 'oldSWE      = ', oldSWE
          ! write(*,'(a,1x,f20.10)') 'newSWE      = ', scalarSWE
          ! write(*,'(a,1x,f20.10)') 'delSWE      = ', delSWE
          ! write(*,'(a,1x,f20.10)') 'effRainfall = ', effRainfall*data_step
          ! write(*,'(a,1x,f20.10)') 'effSnowfall = ', effSnowfall*data_step
          ! write(*,'(a,1x,f20.10)') 'sublimation = ', averageSnowSublimation*data_step
          ! write(*,'(a,1x,f20.10)') 'snwDrainage = ', averageSnowDrainage*iden_water*data_step
          ! write(*,'(a,1x,f20.10)') 'sfcMeltPond = ', sfcMeltPond
          ! write(*,'(a,1x,f20.10)') 'SWE_BalErr  = ', massBalance
          ! message=trim(message)//'SWE does not balance'
          ! err=20; return
        ! endif  ! if failed mass balance check
      endif  ! if snow layers exist

      ! -----
      ! * balance checks for soil...
      ! ----------------------------

      ! compute the liquid water and ice content at the end of the time step
      scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))
      scalarTotalSoilIce = sum(iden_water*mLayerVolFracIce(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))   ! NOTE: no expansion of soil, hence use iden_water

      ! get the total water in the soil (liquid plus ice) at the end of the time step (kg m-2)
      scalarTotalSoilWat = scalarTotalSoilLiq + scalarTotalSoilIce

      ! get the total aquifer storage at the start of the time step (kg m-2)
      balanceAquifer1 = scalarAquiferStorage*iden_water

      ! get the input and output to/from the soil zone (kg m-2)
      balanceSoilInflux        = averageSoilInflux*iden_water*data_step
      balanceSoilBaseflow      = averageSoilBaseflow*iden_water*data_step
      balanceSoilDrainage      = averageSoilDrainage*iden_water*data_step
      balanceSoilET            = (averageCanopyTranspiration + averageGroundEvaporation)*data_step
      balanceSoilCompress      = scalarSoilCompress*data_step

      ! check the individual layers
      ! if(printBalance)then
        ! write(*,'(a,1x,10(f12.8,1x))') 'liqSoilInit       = ', liqSoilInit
        ! write(*,'(a,1x,10(f12.8,1x))') 'volFracLiq        = ', mLayerVolFracLiq(nSnow+1:nLayers)
        ! write(*,'(a,1x,10(f12.8,1x))') 'iLayerLiqFluxSoil = ', flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat*iden_water*data_step
        ! write(*,'(a,1x,10(f12.8,1x))') 'mLayerLiqFluxSoil = ', flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat*data_step
        ! write(*,'(a,1x,10(f12.8,1x))') 'change volFracLiq = ', mLayerVolFracLiq(nSnow+1:nLayers) - liqSoilInit
        ! deallocate(liqSoilInit, stat=err)
        ! if(err/=0)then
        ! message=trim(message)//'unable to deallocate space for the initial soil moisture'
        ! err=20; return
        ! endif
      ! endif

      ! check the soil water balance
      scalarSoilWatBalError  = scalarTotalSoilWat - (balanceSoilWater0 + (balanceSoilInflux + balanceSoilET - balanceSoilBaseflow - balanceSoilDrainage - balanceSoilCompress) )
      ! if(abs(scalarSoilWatBalError) > absConvTol_liquid*iden_water*10._rkind .and. checkMassBalance_ds)then  ! NOTE: kg m-2, so need coarse tolerance to account for precision issues
        ! write(*,*)               'solution method       = ', ixSolution
        ! write(*,'(a,1x,f20.10)') 'data_step             = ', data_step
        ! write(*,'(a,1x,f20.10)') 'balanceSoilCompress   = ', balanceSoilCompress
        ! write(*,'(a,1x,f20.10)') 'scalarTotalSoilLiq    = ', scalarTotalSoilLiq
        ! write(*,'(a,1x,f20.10)') 'scalarTotalSoilIce    = ', scalarTotalSoilIce
        ! write(*,'(a,1x,f20.10)') 'balanceSoilWater0     = ', balanceSoilWater0
        ! write(*,'(a,1x,f20.10)') 'balanceSoilWater1     = ', scalarTotalSoilWat
        ! write(*,'(a,1x,f20.10)') 'balanceSoilInflux     = ', balanceSoilInflux
        ! write(*,'(a,1x,f20.10)') 'balanceSoilBaseflow   = ', balanceSoilBaseflow
        ! write(*,'(a,1x,f20.10)') 'balanceSoilDrainage   = ', balanceSoilDrainage
        ! write(*,'(a,1x,f20.10)') 'balanceSoilET         = ', balanceSoilET
        ! write(*,'(a,1x,f20.10)') 'scalarSoilWatBalError = ', scalarSoilWatBalError
        ! message=trim(message)//'soil hydrology does not balance'
        ! err=20; return
      ! end if

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
                        canopyDepth, &
                        nSnow,ixDomainType,&
                        ixControlVolume,ixStateType, &
                        ! input: ice content change
                        scalarCanopyIce,       & ! intent(in):    value for canopy ice content (kg m-2)
                        mLayerVolFracIce,      & ! intent(in):    vector of volumetric ice water content (-)
                        ! input/output: enthalpy
                        scalarCanopyEnthalpy,  & ! intent(inout): enthTemp to enthalpy of the vegetation canopy (J m-3)
                        mLayerEnthalpy,        & ! intent(inout): enthTemp to enthalpy of each snow+soil layer (J m-3)
                        ! output: error control    
                        err)            ! intent(out): error control
        ! if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      end if

      if(enthalpyStateVec)then ! enthalpy as state variable
        ! initialize the temperature component of enthalpy
        scalarCanopyEnthTemp = scalarCanopyEnthalpy
        mLayerEnthTemp       = mLayerEnthalpy          
        ! compute temperature component of enthalpy for current values       
        call enthTemp_or_enthalpy(&
                        ! input: data structures
                        .false.,               & ! intent(in):    flag to convert enthalpy to enthTemp
                        canopyDepth, &
                        nSnow,ixDomainType,&
                        ixControlVolume,ixStateType, &
                        ! input: ice content change
                        scalarCanopyIce,       & ! intent(in):    value for canopy ice content (kg m-2)
                        mLayerVolFracIce,      & ! intent(in):    vector of volumetric ice water content (-)
                        ! input/output: enthalpy
                        scalarCanopyEnthTemp,  & ! intent(inout): enthalpy to enthTemp of the vegetation canopy (J m-3)
                        mLayerEnthTemp,        & ! intent(inout): enthalpy to enthTemp of each snow+soil layer (J m-3)
                        ! output: error control    
                        err)            ! intent(out): error control
        ! if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif
      ! save the total soil enthalpy
      scalarTotalSoilEnthalpy = sum(mLayerEnthalpy(nSnow+1:nLayers) * mLayerDepth(nSnow+1:nLayers))/sum(mLayerDepth(nSnow+1:nLayers))
      ! save the total snow enthalpy
      if(nSnow>0) scalarTotalSnowEnthalpy = sum(mLayerEnthalpy(1:nSnow) * mLayerDepth(1:nSnow))/sum(mLayerDepth(1:nSnow))

      ! save the surface temperature (just to make things easier to visualize)
      scalarSurfaceTemp = mLayerTemp(1)

  ! overwrite flux data with timestep-average value for all flux_mean vars, hard-coded to not happen
  ! if(.not.backwardsCompatibility)then
  !   do iVar=1,numScalarFluxData
  !     flux_data(iVar) = flux_mean(iVar)
  !   end do
  ! end if


end subroutine
attributes(global) subroutine update_diagnostics(nGRU,ixThCondSnow,ixThCondSoil,&
    computeVegFlux,canopyDepth,&
    specificHeatVeg,maxMassVegetation,fixedThermalCond_snow, &
    iden_soil,thCond_soil,theta_sat,frac_sand,frac_silt,frac_clay,&
    nSnow,nSoil,nLayers,layerType,&
    scalarCanopyIce,scalarCanopyLiquid,&
    mLayerVolFracIce,mLayerVolFracLiq,&
    mLayerHeight,iLayerHeight,&
    scalarBulkVolHeatCapVeg, &
    mLayerVolHtCapBulk,mLayerThermalC,iLayerThermalC,mLayerVolFracAir,&
    scalarSWE, scalarSnowDepth, scalarSfcMeltPond,&
    mLayerTemp,mLayerDepth,mLayerVolFracIceInit,&
    scalarCanopyWat,mLayerVolFracWat,&
    mLayerEnthTemp,mLayerEnthalpy,&
    soil_dens_intr,vGn_alpha,vGn_n,theta_res,vGn_m,&
    mLayerMatricHead,use_lookup,&
    temperature,psiLiq_int,deriv2,&
    computeEnthalpy, enthalpyStateVec,&
    mLayerMatricHeadLiq,&
    sumCanopySublimation, sumSnowSublimation, sumLatHeatCanopyEvap, sumSenHeatCanopy,&
    flux_inner,&
    innerEffRainfall,innerSoilCompress, &
    innerBalance,innerBalanceLayerNrg,innerBalanceLayerMass)
   integer(i4b),intent(in),value :: nGRU
      integer(i4b),intent(in) :: ixThCondSnow,ixThCondSoil
 logical(lgt),intent(in) :: computeVegFlux(:)
 real(rkind),intent(in) :: canopyDepth(:)
 real(rkind),intent(in) :: specificHeatVeg(:),maxMassVegetation(:),fixedThermalCond_snow(:)
 real(rkind),intent(in) :: iden_soil(:,:),thCond_soil(:,:),theta_sat(:,:),frac_sand(:,:),frac_silt(:,:),frac_clay(:,:)
 integer(i4b),intent(in) :: nSnow(:),nLayers(:)
 integer(i4b),intent(in),value :: nSoil
 integer(i4b),intent(in) :: layerType(:,:)
 real(rkind),intent(in) :: scalarCanopyIce(:),scalarCanopyLiquid(:)
 real(rkind),intent(in) :: mLayerVolFracIce(:,:), mLayerVolFracLiq(:,:)
 real(rkind),intent(in) :: mLayerHeight(0:,:),iLayerHeight(0:,:)
 real(rkind),intent(inout) :: scalarBulkVolHeatCapVeg(:)
 real(rkind),intent(inout) :: mLayerVolHtCapBulk(:,:)
 real(rkind),intent(inout) :: mLayerThermalC(:,:)
 real(rkind),intent(inout) :: iLayerThermalC(0:,:)
 real(rkind),intent(inout) :: mLayerVolFracAir(:,:)
 real(rkind),intent(inout) :: scalarSWE(:), scalarSnowDepth(:), scalarSfcMeltPond(:)
 real(rkind),intent(inout) :: mLayerTemp(:,:), mLayerDepth(:,:)
 real(rkind),intent(inout) :: mLayerVolFracIceInit(:,:)
 real(rkind),intent(inout) :: scalarCanopyWat(:), mLayerVolFracWat(:,:)
 real(rkind),intent(inout) :: mLayerEnthTemp(:,:), mLayerEnthalpy(:,:)
 real(rkind),intent(in) :: soil_dens_intr(:,:), vGn_alpha(:,:), vGn_n(:,:),theta_res(:,:), vGn_m(:,:)
 real(rkind),intent(inout) :: mLayerMatricHead(:,:)
 logical(lgt),intent(in) :: use_lookup
 real(rkind),intent(in) :: temperature(:,:,:),psiLiq_int(:,:,:),deriv2(:,:,:)
 logical(lgt),intent(in) :: computeEnthalpy, enthalpyStateVec
 real(rkind),intent(inout) :: mLayerMatricHeadLiq(:,:)
 real(rkind),intent(inout) :: sumCanopySublimation(:), sumSnowSublimation(:), sumLatHeatCanopyEvap(:), sumSenHeatCanopy(:)
 real(rkind),intent(inout) :: flux_inner(:,:)
 real(rkind),intent(inout) :: innerEffRainfall(:),innerSoilCompress(:,:)
 real(rkind),intent(inout) :: innerBalance(:,:), innerBalanceLayerNrg(:,:), innerBalanceLayerMass(:,:)

       integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return

  call update_diagnostics_device(ixThCondSnow,ixThCondSoil,&
    computeVegFlux(iGRU),canopyDepth(iGRU),&
    specificHeatVeg(iGRU),maxMassVegetation(iGRU),fixedThermalCond_snow(iGRU), &
    iden_soil(:,iGRU),thCond_soil(:,iGRU),theta_sat(:,iGRU),frac_sand(:,iGRU),frac_silt(:,iGRU),frac_clay(:,iGRU),&
    nSnow(iGRU),nSoil,nLayers(iGRU),layerType(:,iGRU),&
    scalarCanopyIce(iGRU),scalarCanopyLiquid(iGRU),&
    mLayerVolFracIce(:,iGRU),mLayerVolFracLiq(:,iGRU),&
    mLayerHeight(:,iGRU),iLayerHeight(:,iGRU),&
    scalarBulkVolHeatCapVeg(iGRU), &
    mLayerVolHtCapBulk(:,iGRU),mLayerThermalC(:,iGRU),iLayerThermalC(:,iGRU),mLayerVolFracAir(:,iGRU),&
    scalarSWE(iGRU), scalarSnowDepth(iGRU), scalarSfcMeltPond(iGRU),&
    mLayerTemp(:,iGRU),mLayerDepth(:,iGRU),mLayerVolFracIceInit(:,iGRU),&
    scalarCanopyWat(iGRU),mLayerVolFracWat(:,iGRU),&
    mLayerEnthTemp(:,iGRU),mLayerEnthalpy(:,iGRU),&
    soil_dens_intr(:,iGRU),vGn_alpha(:,iGRU),vGn_n(:,iGRU),theta_res(:,iGRU),vGn_m(:,iGRU),&
    mLayerMatricHead(:,iGRU),use_lookup,&
    temperature(:,:,iGRU),psiLiq_int(:,:,iGRU),deriv2(:,:,iGRU),&
    computeEnthalpy, enthalpyStateVec,&
    mLayerMatricHeadLiq(:,iGRU),&
    sumCanopySublimation(iGRU), sumSnowSublimation(iGRU), sumLatHeatCanopyEvap(iGRU), sumSenHeatCanopy(iGRU),&
    flux_inner(:,iGRU),&
    innerEffRainfall(iGRU),innerSoilCompress(:,iGRU),&
    innerBalance(:,iGRU),innerBalanceLayerNrg(:,iGRU),innerBalanceLayerMass(:,iGRU))
end subroutine

attributes(device) subroutine update_diagnostics_device(ixThCondSnow,ixThCondSoil,&
    computeVegFlux,canopyDepth,&
    specificHeatVeg,maxMassVegetation,fixedThermalCond_snow, &
    iden_soil,thCond_soil,theta_sat,frac_sand,frac_silt,frac_clay,&
    nSnow,nSoil,nLayers,layerType,&
    scalarCanopyIce,scalarCanopyLiq,&
    mLayerVolFracIce,mLayerVolFracLiq,&
    mLayerHeight,iLayerHeight,&
    scalarBulkVolHeatCapVeg, &
    mLayerVolHtCapBulk,mLayerThermalC,iLayerThermalC,mLayerVolFracAir,&
    scalarSWE, scalarSnowDepth, scalarSfcMeltPond,&
    mLayerTemp,mLayerDepth,mLayerVolFracIceInit,&
    scalarCanopyWat,mLayerVolFracWat,&
    mLayerEnthTemp,mLayerEnthalpy,&
    soil_dens_intr,vGn_alpha,vGn_n,theta_res,vGn_m,&
    mLayerMatricHead,use_lookup,&
    temperature,psiLiq_int,deriv2,&
    computeEnthalpy, enthalpyStateVec,&
    mLayerMatricHeadLiq,&
    sumCanopySublimation, sumSnowSublimation, sumLatHeatCanopyEvap, sumSenHeatCanopy,&
    flux_inner,&
    innerEffRainfall,innerSoilCompress,&
    innerBalance,innerBalanceLayerNrg,innerBalanceLayerMass)
  use diagn_evar_module,only:diagn_evar
  use enthalpyTemp_module,only:T2enthTemp_soil
  USE soil_utils_module,only:liquidHead             ! compute the liquid water matric potential

   integer(i4b),intent(in) :: ixThCondSnow,ixThCondSoil
 logical(lgt),intent(in) :: computeVegFlux
 real(rkind),intent(in) :: canopyDepth
 real(rkind),intent(in) :: specificHeatVeg,maxMassVegetation,fixedThermalCond_snow
 real(rkind),intent(in) :: iden_soil(:),thCond_soil(:),theta_sat(:),frac_sand(:),frac_silt(:),frac_clay(:)
 integer(i4b),intent(in) :: nSnow,nSoil,nLayers
 integer(i4b),intent(in) :: layerType(:)
 real(rkind),intent(in) :: scalarCanopyIce,scalarCanopyLiq
 real(rkind),intent(in) :: mLayerVolFracIce(:), mLayerVolFracLiq(:)
 real(rkind),intent(in) :: mLayerHeight(0:),iLayerHeight(0:)
 real(rkind),intent(inout) :: scalarBulkVolHeatCapVeg
 real(rkind),intent(inout) :: mLayerVolHtCapBulk(:)
 real(rkind),intent(inout) :: mLayerThermalC(:)
 real(rkind),intent(inout) :: iLayerThermalC(0:)
 real(rkind),intent(inout) :: mLayerVolFracAir(:)
 real(rkind),intent(inout) :: scalarSWE, scalarSnowDepth, scalarSfcMeltPond
 real(rkind),intent(inout) :: mLayerTemp(:), mLayerDepth(:)
 real(rkind),intent(inout) :: mLayerVolFracIceInit(:)
 real(rkind),intent(inout) :: scalarCanopyWat, mLayerVolFracWat(:)
 real(rkind),intent(inout) :: mLayerEnthTemp(:), mLayerEnthalpy(:)
 real(rkind),intent(in) :: soil_dens_intr(:), vGn_alpha(:), vGn_n(:),theta_res(:), vGn_m(:)
 real(rkind),intent(inout) :: mLayerMatricHead(:)
 logical(lgt),intent(in) :: use_lookup
 real(rkind),intent(in) :: temperature(:,:),psiLiq_int(:,:),deriv2(:,:)
 logical(lgt),intent(in) :: computeEnthalpy, enthalpyStateVec
 real(rkind),intent(inout) :: mLayerMatricHeadLiq(:)
 real(rkind),intent(inout) :: sumCanopySublimation, sumSnowSublimation, sumLatHeatCanopyEvap, sumSenHeatCanopy
 real(rkind),intent(inout) :: flux_inner(:)
 real(rkind),intent(inout) :: innerEffRainfall,innerSoilCompress(:)
 real(rkind),intent(inout) :: innerBalance(:),innerBalanceLayerNrg(:),innerBalanceLayerMass(:)
 integer(i4b) :: err, iSoil

  ! *** compute diagnostic variables for each layer...
  ! --------------------------------------------------
  ! NOTE: this needs to be done AFTER volicePack, since layers may have been sub-divided and/or merged, and need to specifically send in canopy depth
  call diagn_evar(ixThCondSnow,ixThCondSoil,&
    computeVegFlux,canopyDepth,&
    specificHeatVeg,maxMassVegetation,fixedThermalCond_snow, &
    iden_soil,thCond_soil,theta_sat,frac_sand,frac_silt,frac_clay,&
    nSnow,nSoil,nLayers,layerType,&
    scalarCanopyIce,scalarCanopyLiq,&
    mLayerVolFracIce,mLayerVolFracLiq,&
    mLayerHeight,iLayerHeight,&
    scalarBulkVolHeatCapVeg, &
    mLayerVolHtCapBulk,mLayerThermalC,iLayerThermalC,mLayerVolFracAir,err)

  ! *** compute melt of the "snow without a layer"...
  ! -------------------------------------------------
  ! NOTE: forms a surface melt pond, which drains into the upper-most soil layer through the time step
  ! (check for the special case of "snow without a layer")
  ! this pond melts evenly over entire time of maxstep until it gets recomputed because based on SWE when computed
  if(nSnow==0) then
    call implctMelt(&
                    ! input/output: integrated snowpack properties
                    scalarSWE,               & ! intent(inout): snow water equivalent (kg m-2)
                    scalarSnowDepth,         & ! intent(inout): snow depth (m)
                    scalarSfcMeltPond,       & ! intent(out):   surface melt pond (kg m-2)
                    ! input/output: properties of the upper-most soil layer
                    mLayerTemp(nSnow+1),        & ! intent(inout): surface layer temperature (K)
                    mLayerDepth(nSnow+1),       & ! intent(inout): surface layer depth (m)
                    mLayerVolHtCapBulk(nSnow+1),& ! intent(in):    surface layer volumetric heat capacity (J m-3 K-1)
                    ! output: error control
                    err                                        ) ! intent(out): error control
    ! if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
  endif

  ! save volumetric ice content at the start of the step
  ! NOTE: used for volumetric loss due to melt-freeze
  mLayerVolFracIceInit = mLayerVolFracIce

  ! make sure have consistent state variables to start, later done in updateVars
  ! associate local variables with information in the data structures
  ! compute the total water content in the vegetation canopy
  scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce  ! kg m-2

  ! compute the total water content in snow and soil
  ! NOTE: no ice expansion allowed for soil
  if(nSnow>0)&
    mLayerVolFracWat(      1:nSnow  ) = mLayerVolFracLiq(      1:nSnow  ) + mLayerVolFracIce(      1:nSnow  )*(iden_ice/iden_water)
  mLayerVolFracWat(nSnow+1:nLayers)   = mLayerVolFracLiq(nSnow+1:nLayers) + mLayerVolFracIce(nSnow+1:nLayers)

  ! compute enthalpy of the top soil layer if changed with surface melt pond
  if( (enthalpyStateVec .or. computeEnthalpy) .and. nSnow==0 .and. scalarSWE>0._rkind )then
    call T2enthTemp_soil(&
              use_lookup,                                               & ! intent(in):  flag to use the lookup table for soil enthalpy
              soil_dens_intr(1),                                           & ! intent(in):  intrinsic soil density (kg m-3)
              vGn_alpha(1),vGn_n(1),theta_sat(1),theta_res(1),vGn_m(1), & ! intent(in):  van Genutchen soil parameters
              1_i4b,                                                    & ! intent(in):  index of the control volume within the domain
              temperature(:,1),psiLiq_int(:,1),deriv2(:,1),                                              & ! intent(in):  lookup table data structure
              realMissing,                                              & ! intent(in):  lower value of integral (not computed)
              mLayerTemp(nSnow+1),         & ! intent(in):  surface layer temperature (K)
              mLayerMatricHead(1),                                      & ! intent(in):  surface layer matric head (m)
              mLayerEnthTemp(nSnow+1))       ! intent(out): temperature component of enthalpy soil layer (J m-3)
    mLayerEnthalpy(nSnow+1) = mLayerEnthTemp(nSnow+1) - iden_water * LH_fus * mLayerVolFracIce(nSnow+1)
  end if
    
  ! compute the liquid water matric potential (m)
  ! NOTE: include ice content as part of the solid porosity - major effect of ice is to reduce the pore size; ensure that effSat=1 at saturation
  ! (from Zhao et al., J. Hydrol., 1997: Numerical analysis of simultaneous heat and mass transfer...)
  do iSoil=1,nSoil
    call liquidHead(mLayerMatricHead(iSoil),mLayerVolFracLiq(nSnow+iSoil),mLayerVolFracIce(nSnow+iSoil), & ! input:  state variables
              vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m(iSoil),              & ! input:  parameters
              matricHeadLiq=mLayerMatricHeadLiq(iSoil),                                                  & ! output: liquid water matric potential (m)
              err=err)                                                                    ! output: error control
    ! if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  end do  ! looping through soil layers (computing liquid water matric potential)
        
  ! initialize sublimation sums to average over whole_step
  sumCanopySublimation = 0._rkind
  sumSnowSublimation   = 0._rkind
  sumLatHeatCanopyEvap = 0._rkind
  sumSenHeatCanopy     = 0._rkind
  ! initialize fluxes to average over whole_step (averaged over substep in varSubStep)
  ! do iVar=1,size(flux_inner%var)
    flux_inner = 0._rkind
  ! end do
  innerEffRainfall  = 0._rkind ! mean total effective rainfall over snow
  innerSoilCompress = 0._rkind ! mean total soil compression
  innerBalance = 0._rkind ! mean total balance array
  innerBalanceLayerNrg = 0._rkind ! mean total balance of energy in layers
  innerBalanceLayerMass = 0._rkind ! mean total balance of mass in layers


end subroutine

attributes(global) subroutine resize_snow_kernel(nGRU,ix_snowLayers,nSnow,nLayers,nSoil, &
 mLayerTemp,mLayerVolFracIce,mLayerVolFracLiq, &
mLayerVolFracWat,mLayerEnthalpy,mLayerDepth, mLayerHeight, iLayerHeight, &
mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass,&
flux_data, &
iLayerConductiveFlux_start,iLayerConductiveFlux_end, &
iLayerAdvectiveFlux_start,iLayerAdvectiveFlux_end, &
iLayerNrgFlux_start,iLayerNrgFlux_end, mLayerNrgFlux_start,mLayerNrgFlux_end, &
iLayerLiqFluxSnow_start,iLayerLiqFluxSnow_end,mLayerLiqFluxSnow_start,mLayerLiqFluxSnow_end, &
    layerType,ixHydType, &
    ixSnowSoilNrg, ixSnowOnlyNrg, &
    ixSnowSoilHyd, ixSnowOnlyHyd, &
    ixNrgLayer, ixHydLayer,&
    tooMuchMelt, zMin,&
    zMinLayer1,zMinLayer2,zminLayer3,zminLayer4,zminLayer5, &
    modifiedLayers,snowfrz_scale,h_lookup_d,t_lookup_d, &
    scalarSnowDepth, scalarSWE, &
    zmax, &
    zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower, &
    zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper, &
    canopySrad, alb_method, &
    albedoMax, albedoMaxVisible, albedoMaxNearIR, &
    spectralSnowAlbedoDiffuse, spectralSnowAlbedoDirect, scalarSnowAlbedo, &
    veryBig,Frad_vis,maxSnowLayers,&
    firstSubStep,modifiedVegState,&
  nCasNrg,nVegNrg,nVegMass,&
  nVegState,nNrgState,nWatState,&
  nMatState,nMassState,nState,&
  ixNrgCanair,ixNrgCanopy,ixHydCanopy,ixWatAquifer,&
  ixSoilState,&
  ixControlVolume,ixDomainType,ixStateType,&
  computeVegFlux,includeAquifer, &
  vGn_alpha,vGn_n,theta_sat,theta_res, vGn_m, soil_dens_intr,&
  temperature,psiLiq_int,deriv2,use_lookup,&
  mLayerMatricHead,&
  enthalpyStateVec,computeEnthalpy)
   integer(i4b),intent(in),value :: nGRU
        integer(i4b),intent(in) :: ix_snowLayers
    integer(i4b),intent(inout) :: nSnow(:),nLayers(:)
    integer(i4b),intent(in),value :: nSoil
    real(rkind),intent(inout) :: mLayerTemp(:,:)
    real(rkind),intent(inout) :: mLayerVolFracIce(:,:)
    real(rkind),intent(inout) :: mLayerVolFracLiq(:,:)
    real(rkind),intent(inout) :: mLayerVolFracWat(:,:)
    real(rkind),intent(inout) :: mLayerEnthalpy(:,:)
    real(rkind),intent(inout) :: mLayerDepth(:,:)
    real(rkind),intent(inout) :: mLayerHeight(0:,:)
    real(rkind),intent(inout) :: iLayerHeight(0:,:)
  real(rkind),intent(inout) :: mLayerVolHtCapBulk(:,:), mLayerCm(:,:), mLayerThermalC(:,:), iLayerThermalC(0:,:)
  real(rkind),intent(inout) :: mLayerEnthTemp(:,:), mLayerFracLiqSnow(:,:), mLayerThetaResid(:,:), mLayerPoreSpace(:,:)
  real(rkind),intent(inout) :: mLayerMeltFreeze(:,:), mLayerVolFracAir(:,:), balanceLayerNrg(:,:), balanceLayerMass(:,:)
    real(rkind),intent(inout) :: flux_data(:,:)
  integer(i4b),intent(in),value :: iLayerConductiveFlux_start,iLayerConductiveFlux_end
  integer(i4b),intent(in),value :: iLayerAdvectiveFlux_start,iLayerAdvectiveFlux_end
  integer(i4b),intent(in),value :: iLayerNrgFlux_start,iLayerNrgFlux_end, mLayerNrgFlux_start,mLayerNrgFlux_end
  integer(i4b),intent(in),value :: iLayerLiqFluxSnow_start,iLayerLiqFluxSnow_end,mLayerLiqFluxSnow_start,mLayerLiqFluxSnow_end
    integer(i4b),intent(inout) :: layerType(:,:), ixHydType(:,:)
    integer(i4b),intent(inout) :: ixSnowSoilNrg(:,:), ixSnowOnlyNrg(:,:)
    integer(i4b),intent(inout) :: ixSnowSoilHyd(:,:), ixSnowOnlyHyd(:,:)
    integer(i4b),intent(inout) :: ixNrgLayer(:,:), ixHydLayer(:,:)
    logical(lgt),intent(in),value :: tooMuchMelt
    real(rkind),intent(in) :: zMin(:)
    real(rkind),intent(in) :: zMinLayer1(:),zMinLayer2(:),zminLayer3(:),zminLayer4(:),zminLayer5(:)
    logical(lgt),intent(inout) :: modifiedLayers(:)
    real(rkind),intent(in) :: snowfrz_scale(:)
    real(rkind),intent(inout) :: h_lookup_d(:), t_lookup_d(:)
    real(rkind),intent(inout) :: scalarSnowDepth(:), scalarSWE(:)
    real(rkind),intent(in) :: zmax(:)
    real(rkind),intent(in) :: zmaxLayer1_lower(:), zmaxLayer2_lower(:), zmaxLayer3_lower(:), zmaxLayer4_lower(:)
    real(rkind),intent(in) :: zmaxLayer1_upper(:), zmaxLayer2_upper(:), zmaxLayer3_upper(:), zmaxLayer4_upper(:)
    integer(i4b),intent(in) :: canopySrad, alb_method
    real(rkind),intent(in) :: albedoMax(:), albedoMaxVisible(:), albedoMaxNearIR(:)
    real(rkind),intent(inout) :: spectralSnowAlbedoDiffuse(:,:), spectralSnowAlbedoDirect(:,:)
    real(rkind),intent(inout) :: scalarSnowAlbedo(:)
    real(rkind),intent(in),value :: veryBig
    real(rkind),intent(in) :: Frad_vis(:)
    integer(i4b),intent(in),value :: maxSnowLayers
  logical(lgt),intent(in),value :: firstSubStep
  logical(lgt),intent(in) :: modifiedVegState(:)
  integer(i4b),intent(inout) :: nCasNrg(:),nVegNrg(:),nVegMass(:)
  integer(i4b),intent(inout) :: nVegState(:),nNrgState(:),nWatState(:)
  integer(i4b),intent(inout) :: nMatState(:),nMassState(:),nState(:)
  integer(i4b),intent(inout) :: ixNrgCanair(:),ixNrgCanopy(:),ixHydCanopy(:),ixWatAquifer(:)
  integer(i4b),intent(inout) :: ixSoilState(:,:)
  integer(i4b),intent(inout) :: ixControlVolume(:,:), ixDomainType(:,:), ixStateType(:,:)
  logical(lgt),intent(in) :: computeVegFlux(:)
  logical(lgt),intent(in) :: includeAquifer
  real(rkind),intent(in) :: vGn_alpha(:,:), vGn_n(:,:),theta_sat(:,:),theta_res(:,:), vGn_m(:,:), soil_dens_intr(:,:)
  real(rkind),intent(inout) :: temperature(:,:,:), psiLiq_int(:,:,:), deriv2(:,:,:)
  logical(lgt),intent(in) :: use_lookup
  real(rkind),intent(inout) :: mLayerMatricHead(:,:)
  logical(lgt),intent(in) :: enthalpyStateVec,computeEnthalpy

       integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return

  call resize_snow(mLayerTemp(:,iGRU), mLayerVolFracIce(:,iGRU), mLayerVolFracLiq(:,iGRU), mLayerVolFracWat(:,iGRU), &
  mLayerEnthalpy(:,iGRU), mLayerDepth(:,iGRU), mLayerHeight(:,iGRU), iLayerHeight(:,iGRU), &
  nLayers(iGRU), &
  mLayerVolHtCapBulk(:,iGRU), mLayerCm(:,iGRU), mLayerThermalC(:,iGRU), iLayerThermalC(:,iGRU), &
    mLayerEnthTemp(:,iGRU), mLayerFracLiqSnow(:,iGRU), mLayerThetaResid(:,iGRU), mLayerPoreSpace(:,iGRU), &
    mLayerMeltFreeze(:,iGRU), mLayerVolFracAir(:,iGRU), balanceLayerNrg(:,iGRU), balanceLayerMass(:,iGRU), &
flux_data(iLayerConductiveFlux_start:iLayerConductiveFlux_end,iGRU), flux_data(iLayerAdvectiveFlux_start:iLayerAdvectiveFlux_end,iGRU), &
  flux_data(iLayerNrgFlux_start:iLayerNrgFlux_end,iGRU), flux_data(mLayerNrgFlux_start:mLayerNrgFlux_end,iGRU), &
  flux_data(iLayerLiqFluxSnow_start:iLayerLiqFluxSnow_end,iGRU),flux_data(mLayerLiqFluxSnow_start:mLayerLiqFluxSnow_end,iGRU), &
  layerType(:,iGRU), ixHydType(:,iGRU), &
  ixSnowSoilNrg(:,iGRU), ixSnowOnlyNrg(:,iGRU), &
  ixSnowSoilHyd(:,iGRU), ixSnowOnlyHyd(:,iGRU), &
  ixNrgLayer(:,iGRU), ixHydLayer(:,iGRU), nSnow(iGRU), &
  ix_snowLayers, zMax(iGRU), &
  zmaxLayer1_lower(iGRU), zmaxLayer2_lower(iGRU), zmaxLayer3_lower(iGRU), zmaxLayer4_lower(iGRU), &
  zmaxLayer1_upper(iGRU), zmaxLayer2_upper(iGRU), zmaxLayer3_upper(iGRU), zmaxLayer4_upper(iGRU), &
  maxSnowLayers, veryBig, &
  canopySrad, alb_method, &
  albedoMax(iGRU), albedoMaxVisible(iGRU), albedoMaxNearIR(iGRU), &
  spectralSnowAlbedoDiffuse(:,iGRU), spectralSnowAlbedoDirect(:,iGRU), scalarSnowAlbedo(iGRU), &
  scalarSnowDepth(iGRU), scalarSWE(iGRU), &
  Frad_vis(iGRU), nSoil,tooMuchMelt, &
  zMin(iGRU), zMinLayer1(iGRU), zMinLayer2(iGRU), zminLayer3(iGRU), zminLayer4(iGRU), zminLayer5(iGRU), &
  snowfrz_scale(iGRU), h_lookup_d, t_lookup_d, modifiedLayers(iGRU),&
  firstSubStep,modifiedVegState(iGRU),&
  nCasNrg(iGRU),nVegNrg(iGRU),nVegMass(iGRU),&
  nVegState(iGRU),nNrgState(iGRU),nWatState(iGRU),&
  nMatState(iGRU),nMassState(iGRU),nState(iGRU),&
  ixNrgCanair(iGRU),ixNrgCanopy(iGRU),ixHydCanopy(iGRU),ixWatAquifer(iGRU),&
  ixSoilState(:,iGRU),&
  ixControlVolume(:,iGRU),ixDomainType(:,iGRU),ixStateType(:,iGRU),&
  computeVegFlux(iGRU),includeAquifer, &
  vGn_alpha(:,iGRU),vGn_n(:,iGRU),theta_sat(:,iGRU),theta_res(:,iGRU), vGn_m(:,iGRU), soil_dens_intr(:,iGRU),&
  temperature(:,:,iGRU),psiLiq_int(:,:,iGRU),deriv2(:,:,iGRU),use_lookup,&
  mLayerMatricHead(:,iGRU),&
  enthalpyStateVec,computeEnthalpy)

end subroutine

attributes(device) subroutine resize_snow(mLayerTemp, mLayerVolFracIce, mLayerVolFracLiq, mLayerVolFracWat, &
  mLayerEnthalpy, mLayerDepth, mLayerHeight, iLayerHeight, &
  nLayers, &
  mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
    mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
    mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass, &
    iLayerConductiveFlux, iLayerAdvectiveFlux, &
  iLayerNrgFlux, mLayerNrgFlux, &
  iLayerLiqFluxSnow,mLayerLiqFluxSnow, &
  layerType, ixHydType, &
  ixSnowSoilNrg, ixSnowOnlyNrg, &
  ixSnowSoilHyd, ixSnowOnlyHyd, &
  ixNrgLayer, ixHydLayer, nSnow, &
  ix_snowLayers, zMax, &
  zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower, &
  zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper, &
  maxSnowLayers, veryBig, &
  canopySrad, alb_method, &
  albedoMax, albedoMaxVisible, albedoMaxNearIR, &
  spectralSnowAlbedoDiffuse, spectralSnowAlbedoDirect, scalarSnowAlbedo, &
  scalarSnowDepth, scalarSWE, &
  Frad_vis, nSoil,tooMuchMelt, &
  zMin, zMinLayer1, zMinLayer2, zminLayer3, zminLayer4, zminLayer5, &
  snowfrz_scale, h_lookup_d, t_lookup_d, modifiedLayers,&
  firstSubStep,modifiedVegState,&
  nCasNrg,nVegNrg,nVegMass,&
  nVegState,nNrgState,nWatState,&
  nMatState,nMassState,nState,&
  ixNrgCanair,ixNrgCanopy,ixHydCanopy,ixWatAquifer,&
  ixSoilState,&
  ixControlVolume,ixDomainType,ixStateType,&
  computeVegFlux,includeAquifer, &
  vGn_alpha,vGn_n,theta_sat,theta_res, vGn_m, soil_dens_intr,&
  temperature,psiLiq_int,deriv2,use_lookup,&
  mLayerMatricHead,&
  enthalpyStateVec,computeEnthalpy)
  use volicePack_module,only:volicePack_device
  use indexState_module,only:indexState_d
  use enthalpyTemp_module,only:T2enthTemp_snow
  use enthalpyTemp_module,only:T2enthTemp_soil

  real(rkind),intent(inout) :: mLayerTemp(:), mLayerVolFracIce(:), mLayerVolFracLiq(:), mLayerVolFracWat(:)
  real(rkind),intent(inout) :: mLayerEnthalpy(:), mLayerDepth(:), mLayerHeight(0:), iLayerHeight(0:)
  integer(i4b),intent(inout) :: nLayers
  real(rkind),intent(inout) :: mLayerVolHtCapBulk(:), mLayerCm(:), mLayerThermalC(:), iLayerThermalC(0:)
  real(rkind),intent(inout) :: mLayerEnthTemp(:), mLayerFracLiqSnow(:), mLayerThetaResid(:), mLayerPoreSpace(:)
  real(rkind),intent(inout) :: mLayerMeltFreeze(:), mLayerVolFracAir(:), balanceLayerNrg(:), balanceLayerMass(:)
  real(rkind),intent(inout) :: iLayerConductiveFlux(0:)
  real(rkind),intent(inout) :: iLayerAdvectiveFlux(0:)
  real(rkind),intent(inout) :: iLayerNrgFlux(0:)
  real(rkind),intent(inout) :: mLayerNrgFlux(:)
  real(rkind),intent(inout) :: iLayerLiqFluxSnow(0:)
  real(rkind),intent(inout) :: mLayerLiqFluxSnow(:)
  integer(i4b),intent(inout) :: layerType(:), ixHydType(:)
  integer(i4b),intent(inout) :: ixSnowSoilNrg(:), ixSnowOnlyNrg(:)
  integer(i4b),intent(inout) :: ixSnowSoilHyd(:), ixSnowOnlyHyd(:)
  integer(i4b),intent(inout) :: ixNrgLayer(:), ixHydLayer(:)
  integer(i4b),intent(inout) :: nSnow
  integer(i4b),intent(in) :: ix_snowLayers
  real(rkind),intent(in) :: zmax
  real(rkind),intent(in) :: zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower
  real(rkind),intent(in) :: zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper
  real(rkind),intent(in) :: veryBig
  integer(i4b),intent(in) :: maxSnowLayers
  integer(i4b),intent(in) :: canopySrad, alb_method
  real(rkind),intent(in) :: albedoMax, albedoMaxVisible, albedoMaxNearIR
  real(rkind),intent(inout) :: spectralSnowAlbedoDiffuse(:), spectralSnowAlbedoDirect(:), scalarSnowAlbedo
  real(rkind),intent(inout) :: scalarSnowDepth, scalarSWE
  real(rkind),intent(in) :: Frad_vis
  integer(i4b),intent(in) :: nSoil
  logical(lgt),intent(in) :: tooMuchMelt
  real(rkind),intent(in) :: zMin, zMinLayer1, zMinLayer2, zminLayer3, zminLayer4, zminLayer5
  real(rkind),intent(in) :: snowfrz_scale
  real(rkind),intent(in) :: h_lookup_d(:), t_lookup_d(:)
  logical(lgt),intent(inout) :: modifiedLayers
  logical(lgt),intent(in) :: firstSubStep,modifiedVegState
  integer(i4b),intent(inout) :: nCasNrg,nVegNrg,nVegMass
  integer(i4b),intent(inout) :: nVegState,nNrgState,nWatState
  integer(i4b),intent(inout) :: nMatState,nMassState,nState
  integer(i4b),intent(inout) :: ixNrgCanair,ixNrgCanopy,ixHydCanopy,ixWatAquifer
  integer(i4b),intent(inout) :: ixSoilState(:)
  integer(i4b),intent(inout) :: ixControlVolume(:), ixDomainType(:), ixStateType(:)
  logical(lgt),intent(in) :: computeVegFlux,includeAquifer
  real(rkind),intent(in) :: vGn_alpha(:), vGn_n(:),theta_sat(:),theta_res(:), vGn_m(:), soil_dens_intr(:)
  real(rkind),intent(inout) :: temperature(:,:), psiLiq_int(:,:), deriv2(:,:)
  logical(lgt),intent(in) :: use_lookup
  real(rkind),intent(inout) :: mLayerMatricHead(:)
  logical(lgt),intent(in) :: enthalpyStateVec,computeEnthalpy

  integer(i4b) :: err, iLayer,iSoil


  ! re-assign dimension lengths
  nSnow   = count(layerType==iname_snow)
  nLayers = nSnow+nSoil

  call volicePack_device(mLayerTemp, mLayerVolFracIce, mLayerVolFracLiq, mLayerVolFracWat, &
  mLayerEnthalpy, mLayerDepth, mLayerHeight, iLayerHeight, &
  nLayers, &
  mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
    mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
    mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass, &
    iLayerConductiveFlux, iLayerAdvectiveFlux, &
  iLayerNrgFlux, mLayerNrgFlux, &
  iLayerLiqFluxSnow,mLayerLiqFluxSnow, &
  layerType, ixHydType, &
  ixSnowSoilNrg, ixSnowOnlyNrg, &
  ixSnowSoilHyd, ixSnowOnlyHyd, &
  ixNrgLayer, ixHydLayer, nSnow, &
  ix_snowLayers, zMax, &
  zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower, &
  zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper, &
  maxSnowLayers, veryBig, &
  canopySrad, alb_method, &
  albedoMax, albedoMaxVisible, albedoMaxNearIR, &
  spectralSnowAlbedoDiffuse, spectralSnowAlbedoDirect, scalarSnowAlbedo, &
  scalarSnowDepth, scalarSWE, &
  Frad_vis, nSoil,tooMuchMelt, &
  zMin, zMinLayer1, zMinLayer2, zminLayer3, zminLayer4, zminLayer5, &
  snowfrz_scale, h_lookup_d, t_lookup_d, modifiedLayers)

! compute the indices for the model state variables
if(firstSubStep .or. modifiedVegState .or. modifiedLayers)then
  call indexState_d(computeVegFlux,         & ! intent(in):    flag to denote if computing the vegetation flux
                  includeAquifer,         & ! intent(in):    flag to denote if included the aquifer
                  nSnow,nSoil,nLayers,&
                  nCasNrg,nVegNrg,nVegMass,&
                  nVegState,nNrgState,nWatState,&
                  nMatState,nMassState,nState,&
                  ixNrgCanair,ixNrgCanopy,ixHydCanopy,ixWatAquifer,&
                  ixNrgLayer,ixHydLayer,&
                  ixSoilState,&
                  ixControlVolume,ixDomainType,ixStateType,&
                  err)             ! intent(out):   error control
  ! if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
end if
  if( (enthalpyStateVec .or. computeEnthalpy) .and. modifiedLayers )then

    if(nSnow>0)then 
      do iLayer=1,nSnow
        mLayerVolFracWat(iLayer) = mLayerVolFracLiq(iLayer) + mLayerVolFracIce(iLayer)*(iden_ice/iden_water)
        ! compute enthalpy for snow layers
        call T2enthTemp_snow(&
                     snowfrz_scale,             & ! intent(in):  scaling parameter for the snow freezing curve  (K-1)
                     mLayerTemp(iLayer),        & ! intent(in):  layer temperature (K)
                     mLayerVolFracWat(iLayer),  & ! intent(in):  volumetric total water content (-)
                     mLayerEnthTemp(iLayer))      ! intent(out): temperature component of enthalpy of each snow layer (J m-3)
        mLayerEnthalpy(iLayer) = mLayerEnthTemp(iLayer) - iden_ice * LH_fus * mLayerVolFracIce(iLayer)
      end do  ! looping through snow layers
    endif
    do iLayer=nSnow+1,nLayers
      mLayerVolFracWat(iLayer) = mLayerVolFracLiq(iLayer) + mLayerVolFracIce(iLayer)
      ! compute enthalpy for soil layers
      iSoil = iLayer - nSnow
      call T2enthTemp_soil(&
                     use_lookup,                   & ! intent(in):  flag to use the lookup table for soil enthalpy
                     soil_dens_intr(iSoil),        & ! intent(in):  intrinsic soil density (kg m-3)
                     vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m(iSoil), & ! intent(in): soil parameters
                     iSoil,                        & ! intent(in):  index of the control volume within the domain
                     temperature(:,iSoil),psiLiq_int(:,iSoil),deriv2(:,iSoil),                  & ! intent(in):  lookup table data structure
                     realMissing,                  & ! intent(in):  lower value of integral (not computed)
                     mLayerTemp(iLayer),           & ! intent(in):  layer temperature (K)
                     mLayerMatricHead(iSoil),      & ! intent(in):  matric head (m)
                     mLayerEnthTemp(iLayer))         ! intent(out): temperature component of enthalpy soil layer (J m-3)
      mLayerEnthalpy(iLayer) = mLayerEnthTemp(iLayer) - iden_water * LH_fus * mLayerVolFracIce(iLayer)
    end do  ! looping through soil layers
  end if ! (need to recalculate enthalpy state variable)

end subroutine


attributes(global) subroutine mean_inner_step_kernel(nGRU,flux_mean,flux_inner, &
  meanBalance, innerBalance, numScalarFluxData, dt_wght, &
  meanCanopySublimation, sumCanopySublimation, scalarCanopySublimation,&
  meanLatHeatCanopyEvap, sumLatHeatCanopyEvap, scalarLatHeatCanopyEvap,&
  meanSenHeatCanopy, sumSenHeatCanopy, scalarSenHeatCanopy, &
  meanSoilCompress, innerSoilCompress, &
  balanceSnowNrg, balanceSnowMass, balanceSoilNrg, balanceSoilMass, &
  effRainfall, innerEffRainfall,data_step)
  integer(i4b),intent(in),value :: nGRU
    real(rkind) :: flux_mean(:,:), flux_inner(:,:)
  real(rkind) :: meanBalance(:,:), innerBalance(:,:)
  integer(i4b),value :: numScalarFluxData
  real(rkind),value :: dt_wght
  real(rkind) :: meanCanopySublimation(:), sumCanopySublimation(:)
  integer(i4b),intent(in),value :: scalarCanopySublimation
  real(rkind) :: meanLatHeatCanopyEvap(:), sumLatHeatCanopyEvap(:)
  integer(i4b),intent(in),value :: scalarLatHeatCanopyEvap
  real(rkind) :: meanSenHeatCanopy(:), sumSenHeatCanopy(:)
  integer(i4b),intent(in),value :: scalarSenHeatCanopy
  real(rkind) :: meanSoilCompress(:,:), innerSoilCompress(:,:)
  real(rkind) :: balanceSnowNrg(:), balanceSnowMass(:), balanceSoilNrg(:), balanceSoilMass(:)
  real(rkind) :: effRainfall(:), innerEffRainfall(:)
  real(rkind),intent(in),value :: data_step

      integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return

  call means_inner_step(flux_mean(:,iGRU),flux_inner(:,iGRU), &
  meanBalance(:,iGRU), innerBalance(:,iGRU), numScalarFluxData, dt_wght, &
  meanCanopySublimation(iGRU), sumCanopySublimation(iGRU), flux_mean(scalarCanopySublimation,iGRU),&
  meanLatHeatCanopyEvap(iGRU), sumLatHeatCanopyEvap(iGRU), flux_mean(scalarLatHeatCanopyEvap,iGRU),&
  meanSenHeatCanopy(iGRU), sumSenHeatCanopy(iGRU), flux_mean(scalarSenHeatCanopy,iGRU), &
  meanSoilCompress(:,iGRU), innerSoilCompress(:,iGRU), &
  balanceSnowNrg(iGRU), balanceSnowMass(iGRU), balanceSoilNrg(iGRU), balanceSoilMass(iGRU), &
  effRainfall(iGRU), innerEffRainfall(iGRU),data_step)
end subroutine
attributes(device) subroutine means_inner_step(flux_mean,flux_inner, &
  meanBalance, innerBalance, numScalarFluxData, dt_wght, &
  meanCanopySublimation, sumCanopySublimation, scalarCanopySublimation,&
  meanLatHeatCanopyEvap, sumLatHeatCanopyEvap, scalarLatHeatCanopyEvap,&
  meanSenHeatCanopy, sumSenHeatCanopy, scalarSenHeatCanopy, &
  meanSoilCompress, innerSoilCompress, &
  balanceSnowNrg, balanceSnowMass, balanceSoilNrg, balanceSoilMass, &
  effRainfall, innerEffRainfall,data_step)
  real(rkind) :: flux_mean(:), flux_inner(:)
  real(rkind) :: meanBalance(:), innerBalance(:)
  integer(i4b) :: numScalarFluxData
  real(rkind) :: dt_wght
  real(rkind) :: meanCanopySublimation, sumCanopySublimation, scalarCanopySublimation
  real(rkind) :: meanLatHeatCanopyEvap, sumLatHeatCanopyEvap, scalarLatHeatCanopyEvap
  real(rkind) :: meanSenHeatCanopy, sumSenHeatCanopy, scalarSenHeatCanopy
  real(rkind) :: meanSoilCompress(:), innerSoilCompress(:)
  real(rkind) :: balanceSnowNrg, balanceSnowMass, balanceSoilNrg, balanceSoilMass
  real(rkind) :: effRainfall, innerEffRainfall
  real(rkind) :: data_step

  integer(i4b) :: iVar
  do iVar=1,numScalarFluxData
    flux_mean(iVar)    = flux_mean(iVar) + flux_inner(iVar)*dt_wght
  end do
  meanCanopySublimation = meanCanopySublimation + sumCanopySublimation/data_step
  meanLatHeatCanopyEvap = meanLatHeatCanopyEvap + sumLatHeatCanopyEvap/data_step
  meanSenHeatCanopy     = meanSenHeatCanopy     + sumSenHeatCanopy/data_step
  meanSoilCompress(:)   = meanSoilCompress(:)   + innerSoilCompress(:)*dt_wght
  meanBalance(1) = meanBalance(1) + innerBalance(1)*dt_wght
  meanBalance(2) = meanBalance(2) + innerBalance(2)*dt_wght
  meanBalance(3) = meanBalance(3) + innerBalance(3)*dt_wght
  meanBalance(4) = meanBalance(4) + innerBalance(4)*dt_wght
  meanBalance(5) = meanBalance(5) + balanceSnowNrg*dt_wght
  meanBalance(6) = meanBalance(6) + balanceSoilNrg*dt_wght
  meanBalance(7) = meanBalance(7) + balanceSnowMass*dt_wght
  meanBalance(8) = meanBalance(8) + balanceSoilMass*dt_wght

  effRainfall = effRainfall + innerEffRainfall*dt_wght
  scalarCanopySublimation = meanCanopySublimation
  scalarLatHeatCanopyEvap = meanLatHeatCanopyEvap
  scalarSenHeatCanopy     = meanSenHeatCanopy
end subroutine

attributes(global) subroutine increment_means_kernel(nGRU,&
  flux_inner,flux_data, &
  innerSoilCompress, mLayerCompress, &
  nSnow,numScalarFluxData,dt_wght, &
  innerEffRainfall,scalarThroughfallRain,scalarCanopyLiqDrainage, &
  innerBalance, innerBalanceLayerMass, innerBalanceLayerNrg, &
  computeVegFlux, bal_veg, &
  balanceCasNrg, balanceVegNrg, balanceVegMass, canopyDepth, balanceAqMass, &
  balanceLayerMass, balanceLayerNrg, &
  mLayerDepth, layerType, &
  balanceSnowNrg, balanceSnowMass, balanceSoilNrg,balanceSoilMass, &
  bal_soil, bal_snow, bal_aq, &
  nLayers, ixgroundwatr)
  integer(i4b),intent(in),value :: nGRU
  real(rkind) :: flux_inner(:,:), flux_data(:,:)
  real(rkind) :: innerSoilCompress(:,:), mLayerCompress(:,:)
  integer(i4b) :: nSnow(:)
  integer(i4b),value :: numScalarFluxData
  real(rkind),value :: dt_wght
  real(rkind) :: innerEffRainfall(:)
  integer(i4b),value :: scalarThroughfallRain,scalarCanopyLiqDrainage
  real(rkind),intent(inout) :: innerBalance(:,:)
  real(rkind),intent(inout) :: innerBalanceLayerMass(:,:), innerBalanceLayerNrg(:,:)
  logical(lgt),intent(inout) :: computeVegFlux(:), bal_veg(:)
  real(rkind),intent(inout) :: balanceCasNrg(:), balanceVegNrg(:), balanceVegMass(:), canopyDepth(:), balanceAqMass(:)
  real(rkind),intent(inout) :: balanceLayerNrg(:,:), balanceLayerMass(:,:)
  real(rkind),intent(inout) :: mLayerDepth(:,:)
  integer(i4b),intent(inout) :: layerType(:,:)
  real(rkind),intent(inout) :: balanceSnowNrg(:), balanceSnowMass(:), balanceSoilNrg(:),balanceSoilMass(:)
  logical(lgt),intent(inout) :: bal_soil(:), bal_snow(:), bal_aq(:)
  integer(i4b),intent(in) :: nLayers(:), ixgroundwatr

      integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return
  call increment_means_device(flux_inner(:,iGRU),flux_data(:,iGRU),&
  innerSoilCompress(:,iGRU),mLayerCompress(:,iGRU), &
  nSnow(iGRU), numScalarFluxData, dt_wght, &
  innerEffRainfall(iGRU),flux_data(scalarThroughfallRain,iGRU),flux_data(scalarCanopyLiqDrainage,iGRU), &
  innerBalance(:,iGRU), innerBalanceLayerMass(:,iGRU), innerBalanceLayerNrg(:,iGRU), &
  computeVegFlux(iGRU), bal_veg(iGRU), &
  balanceCasNrg(iGRU), balanceVegNrg(iGRU), balanceVegMass(iGRU), canopyDepth(iGRU), balanceAqMass(iGRU), &
  balanceLayerMass(:,iGRU), balanceLayerNrg(:,iGRU), &
  mLayerDepth(:,iGRU), layerType(:,iGRU), &
  balanceSnowNrg(iGRU), balanceSnowMass(iGRU), balanceSoilNrg(iGRU),balanceSoilMass(iGRU), &
  bal_soil(iGRU), bal_snow(iGRU), bal_aq(iGRU), &
  nLayers(iGRU), ixgroundwatr)
end subroutine

attributes(device) subroutine increment_means_device(flux_inner,flux_data,&
  innerSoilCompress,mLayerCompress, &
  nSnow, numScalarFluxData, dt_wght, &
  innerEffRainfall,scalarThroughfallRain,scalarCanopyLiqDrainage, &
  innerBalance, innerBalanceLayerMass, innerBalanceLayerNrg, &
  computeVegFlux, bal_veg, &
  balanceCasNrg, balanceVegNrg, balanceVegMass, canopyDepth, balanceAqMass, &
  balanceLayerMass, balanceLayerNrg, &
  mLayerDepth, layerType, &
  balanceSnowNrg, balanceSnowMass, balanceSoilNrg,balanceSoilMass, &
  bal_soil, bal_snow, bal_aq, &
  nLayers, ixgroundwatr)
  real(rkind),intent(inout) :: flux_inner(:),flux_data(:)
  real(rkind),intent(inout) :: innerSoilCompress(:),mLayerCompress(:)
  integer(i4b),intent(in) :: nSnow, numScalarFluxData
  real(rkind),intent(in) :: dt_wght
  real(rkind),intent(inout) :: innerEffRainfall,scalarThroughfallRain,scalarCanopyLiqDrainage
  real(rkind),intent(inout) :: innerBalance(:)
  real(rkind),intent(inout) :: innerBalanceLayerMass(:), innerBalanceLayerNrg(:)
  logical(lgt),intent(inout) :: computeVegFlux, bal_veg
  real(rkind),intent(inout) :: balanceCasNrg, balanceVegNrg, balanceVegMass, canopyDepth, balanceAqMass
  real(rkind),intent(inout) :: balanceLayerNrg(:), balanceLayerMass(:)
  real(rkind),intent(inout) :: mLayerDepth(:)
  integer(i4b),intent(inout) :: layerType(:)
  real(rkind),intent(inout) :: balanceSnowNrg, balanceSnowMass, balanceSoilNrg,balanceSoilMass
  logical(lgt),intent(inout) :: bal_soil, bal_snow, bal_aq
  integer(i4b),intent(in) :: nLayers, ixgroundwatr


  integer(i4b) :: iVar, iLayer
  real(rkind) :: lyr_wght
  do iVar=1,numScalarFluxData
     flux_inner(iVar)    = flux_inner(iVar) + flux_data(iVar)*dt_wght
  end do
  innerSoilCompress(:) = innerSoilCompress(:) + mLayerCompress*dt_wght
  if (nSnow>0) innerEffRainfall = innerEffRainfall + ( scalarThroughfallRain + scalarCanopyLiqDrainage )*dt_wght

  ! sum the balance of energy and water per state
  if(computeVegFlux)then
    innerBalance(1) = innerBalance(1) + balanceCasNrg*dt_wght ! W m-3
    innerBalance(2) = innerBalance(2) + balanceVegNrg*dt_wght ! W m-3
    innerBalance(3) = innerBalance(3) + balanceVegMass*dt_wght/canopyDepth  ! kg m-3 s-1
    bal_veg = .true.
  endif
  innerBalance(4) = innerBalance(4) + balanceAqMass*dt_wght * iden_water  ! kg m-2 s-1 (no depth to aquifer)
  innerBalanceLayerNrg(:) = innerBalanceLayerNrg(:) + balanceLayerNrg(:)*dt_wght ! W m-3
  innerBalanceLayerMass(:) = innerBalanceLayerMass(:) + balanceLayerMass(:)*dt_wght * iden_water ! kg m-3 s-1

  ! save balance of energy and water per snow+soil layer after inner step, since can change nLayers with outer steps
  balanceLayerNrg(:) = innerBalanceLayerNrg(:)
  balanceLayerMass(:) = innerBalanceLayerMass(:)

  ! compute the balance of energy and water per entire snow and soil domain, in W m-3 and kg m-2 s-1 respectively
  balanceSnowNrg = 0._rkind
  balanceSoilNrg = 0._rkind
  balanceSnowMass = 0._rkind
  balanceSoilMass = 0._rkind
  do iLayer=1,nLayers
    select case (layerType(iLayer))
      case (iname_snow)
        lyr_wght = mLayerDepth(iLayer) / sum( mLayerDepth(1:nSnow) )
        balanceSnowNrg  = balanceSnowNrg + innerBalanceLayerNrg(iLayer)*lyr_wght
        balanceSnowMass = balanceSnowMass + innerBalanceLayerMass(iLayer)*lyr_wght
        bal_snow = .true.
      case (iname_soil)
        lyr_wght = mLayerDepth(iLayer) / sum( mLayerDepth(nSnow+1:nLayers) )
        balanceSoilNrg  = balanceSoilNrg + innerBalanceLayerNrg(iLayer)*lyr_wght
        balanceSoilMass = balanceSoilMass + innerBalanceLayerMass(iLayer)*lyr_wght
        bal_soil = .true.
    end select
  end do
  if (ixgroundwatr == bigBucket) bal_aq = .true. ! aquifer does not change existance with time steps

end subroutine

attributes(global) subroutine sublime_kernel(nGRU,nSnow,nLayers, &
  mLayerMeltFreeze, mLayerVolFracIce, mLayerVolFracIceInit, &
  computeVegFlux, computeEnthalpy,enthalpyStateVec,&
  scalarCanopyIce, sumCanopySublimation, scalarCanopyLiq, &
  superflousSub, superflousNrg, &
  sumLatHeatCanopyEvap, sumSenHeatCanopy, &
  scalarCanopyWat, canopyDepth, &
  specificHeatVeg, maxMassVegetation, &
  whole_step, snowfrz_scale, &
  scalarCanopyTemp, scalarCanopyEnthTemp,scalarCanopyEnthalpy, &
  sumSnowSublimation, &
  mLayerVolFracLiq, mLayerTemp, &
  densScalGrowth,tempScalGrowth, &
  grainGrowthRate,densScalOvrbdn, &
  tempScalOvrbdn,baseViscosity,mLayerDepth,&
  tooMuchSublim,mLayerHeight,iLayerHeight,layerType, &
  scalarSnowDepth, scalarSWE, &
  mLayerEnthalpy, mLayerEnthTemp, &
  sfcMeltPond, scalarSfcMeltPond, &
  mLayerVolFracWat)
  integer(i4b),intent(in),value :: nGRU
  integer(i4b),intent(in) :: nSnow(:),nLayers(:)
  real(rkind),intent(inout) :: mLayerMeltFreeze(:,:), mLayerVolFracIce(:,:), mLayerVolFracIceInit(:,:)
  logical(lgt),intent(in) :: computeVegFlux(:)
  logical(lgt),intent(in) :: computeEnthalpy,enthalpyStateVec
  real(rkind),intent(inout) :: scalarCanopyIce(:), sumCanopySublimation(:), scalarCanopyLiq(:)
  real(rkind),intent(inout) :: superflousSub(:), superflousNrg(:)
  real(rkind),intent(inout) :: sumLatHeatCanopyEvap(:), sumSenHeatCanopy(:)
  real(rkind),intent(inout) :: scalarCanopyWat(:), canopyDepth(:)
  real(rkind),intent(inout) :: specificHeatVeg(:), maxMassVegetation(:)
  real(rkind),intent(in),value :: whole_step
  real(rkind),intent(in) :: snowfrz_scale(:)
  real(rkind),intent(inout) :: scalarCanopyTemp(:), scalarCanopyEnthTemp(:), scalarCanopyEnthalpy(:)
  real(rkind),intent(inout) :: sumSnowSublimation(:)
  real(rkind),intent(inout) :: mLayerVolFracLiq(:,:), mLayerTemp(:,:)
  real(rkind),intent(in) :: densScalGrowth(:),tempScalGrowth(:)
  real(rkind),intent(in) :: grainGrowthRate(:),densScalOvrbdn(:)
  real(rkind),intent(in) :: tempScalOvrbdn(:),baseViscosity(:)
  real(rkind),intent(inout) :: mLayerDepth(:,:)
  logical(lgt),intent(inout) :: tooMuchSublim(:)
  real(rkind),intent(inout) :: mLayerHeight(0:,:), iLayerHeight(0:,:)
  integer(i4b) :: layerType(:,:)
  real(rkind) :: scalarSnowDepth(:), scalarSWE(:)
  real(rkind) :: mLayerEnthalpy(:,:), mLayerEnthTemp(:,:)
  real(rkind) :: sfcMeltPond(:), scalarSfcMeltPond(:)
  real(rkind) :: mLayerVolFracWat(:,:)

    integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return
  call sublime_device(nSnow(iGRU),nLayers(iGRU), &
  mLayerMeltFreeze(:,iGRU), mLayerVolFracIce(:,iGRU), mLayerVolFracIceInit(:,iGRU), &
  computeVegFlux(iGRU), computeEnthalpy,enthalpyStateVec,&
  scalarCanopyIce(iGRU), sumCanopySublimation(iGRU), scalarCanopyLiq(iGRU), &
  superflousSub(iGRU), superflousNrg(iGRU), &
  sumLatHeatCanopyEvap(iGRU), sumSenHeatCanopy(iGRU), &
  scalarCanopyWat(iGRU), canopyDepth(iGRU), &
  specificHeatVeg(iGRU), maxMassVegetation(iGRU), &
  whole_step, snowfrz_scale(iGRU), &
  scalarCanopyTemp(iGRU), scalarCanopyEnthTemp(iGRU),scalarCanopyEnthalpy(iGRU), &
  sumSnowSublimation(iGRU), &
  mLayerVolFracLiq(:,iGRU), mLayerTemp(:,iGRU), &
  densScalGrowth(iGRU),tempScalGrowth(iGRU), &
  grainGrowthRate(iGRU),densScalOvrbdn(iGRU), &
  tempScalOvrbdn(iGRU),baseViscosity(iGRU),mLayerDepth(:,iGRU),&
  tooMuchSublim(iGRU),mLayerHeight(:,iGRU),iLayerHeight(:,iGRU),layerType(:,iGRU), &
  scalarSnowDepth(iGRU), scalarSWE(iGRU), &
  mLayerEnthalpy(:,iGRU), mLayerEnthTemp(:,iGRU), &
  sfcMeltPond(iGRU), scalarSfcMeltPond(iGRU), &
  mLayerVolFracWat(:,iGRU))
end subroutine

attributes(device) subroutine sublime_device(nSnow,nLayers, &
  mLayerMeltFreeze, mLayerVolFracIce, mLayerVolFracIceInit, &
  computeVegFlux, computeEnthalpy,enthalpyStateVec,&
  scalarCanopyIce, sumCanopySublimation, scalarCanopyLiq, &
  superflousSub, superflousNrg, &
  sumLatHeatCanopyEvap, sumSenHeatCanopy, &
  scalarCanopyWat, canopyDepth, &
  specificHeatVeg, maxMassVegetation, &
  whole_step, snowfrz_scale, &
  scalarCanopyTemp, scalarCanopyEnthTemp,scalarCanopyEnthalpy, &
  sumSnowSublimation, &
  mLayerVolFracLiq, mLayerTemp, &
  densScalGrowth,tempScalGrowth, &
  grainGrowthRate,densScalOvrbdn, &
  tempScalOvrbdn,baseViscosity,mLayerDepth,&
  tooMuchSublim,mLayerHeight,iLayerHeight,layerType, &
  scalarSnowDepth, scalarSWE, &
  mLayerEnthalpy, mLayerEnthTemp, &
  sfcMeltPond, scalarSfcMeltPond, &
  mLayerVolFracWat)
  USE enthalpyTemp_module,only:T2enthTemp_veg       ! convert temperature to enthalpy for vegetation
  USE computSnowDepth_module,only:computSnowDepth   ! compute snow depth
  USE var_derive_module,only:calcHeight_d             ! module to calculate height at layer interfaces and layer mid-point
  USE enthalpyTemp_module,only:T2enthTemp_snow      ! convert temperature to enthalpy for snow

  integer(i4b),intent(in) :: nSnow,nLayers
  real(rkind),intent(inout) :: mLayerMeltFreeze(:), mLayerVolFracIce(:), mLayerVolFracIceInit(:)
  logical(lgt),intent(in) :: computeVegFlux,computeEnthalpy,enthalpyStateVec
  real(rkind),intent(inout) :: scalarCanopyIce, sumCanopySublimation, scalarCanopyLiq
  real(rkind),intent(inout) :: superflousSub, superflousNrg
  real(rkind),intent(inout) :: sumLatHeatCanopyEvap, sumSenHeatCanopy
  real(rkind),intent(inout) :: scalarCanopyWat, canopyDepth
  real(rkind),intent(inout) :: specificHeatVeg, maxMassVegetation
  real(rkind),intent(inout) :: whole_step, snowfrz_scale
  real(rkind),intent(inout) :: scalarCanopyTemp, scalarCanopyEnthTemp, scalarCanopyEnthalpy
  real(rkind),intent(inout) :: sumSnowSublimation
  real(rkind),intent(inout) :: mLayerVolFracLiq(:), mLayerTemp(:)
  real(rkind),intent(in) :: densScalGrowth,tempScalGrowth
  real(rkind),intent(in) :: grainGrowthRate,densScalOvrbdn
  real(rkind),intent(in) :: tempScalOvrbdn,baseViscosity
  real(rkind),intent(inout) :: mLayerDepth(:)
  logical(lgt),intent(inout) :: tooMuchSublim
  real(rkind) :: mLayerHeight(0:), iLayerHeight(0:)
  integer(i4b) :: layerType(:)
  real(rkind) :: scalarSnowDepth, scalarSWE
  real(rkind) :: mLayerEnthalpy(:), mLayerEnthTemp(:)
  real(rkind) :: sfcMeltPond, scalarSfcMeltPond
  real(rkind) :: mLayerVolFracWat(:)


  integer(i4b) :: iLayer

  ! compute the melt in each snow and soil layer
  if(nSnow>0)&
    mLayerMeltFreeze(1:nSnow)       = -( mLayerVolFracIce(1:nSnow)         - mLayerVolFracIceInit(1:nSnow) )        *iden_ice
  mLayerMeltFreeze(nSnow+1:nLayers) = -( mLayerVolFracIce(nSnow+1:nLayers) - mLayerVolFracIceInit(nSnow+1:nLayers) )*iden_water
  ! * compute change in canopy ice content due to sublimation...
  ! ------------------------------------------------------------
  if(computeVegFlux)then

    ! remove mass of ice on the canopy
    scalarCanopyIce = scalarCanopyIce + sumCanopySublimation

    ! if removed all ice, take the remaining sublimation from water
    if(scalarCanopyIce < 0._rkind)then
      scalarCanopyLiq = scalarCanopyLiq + scalarCanopyIce
      scalarCanopyIce = 0._rkind
    endif

    ! modify fluxes and mean fluxes if there is insufficient canopy water to support the converged sublimation rate over the whole time step
    if(scalarCanopyLiq < 0._rkind)then
      ! --> superfluous sublimation flux
      superflousSub = -scalarCanopyLiq/whole_step  ! kg m-2 s-1
      superflousNrg = superflousSub*LH_sub     ! W m-2 (J m-2 s-1)
      ! --> update fluxes and states
      sumCanopySublimation = sumCanopySublimation + superflousSub*whole_step
      sumLatHeatCanopyEvap = sumLatHeatCanopyEvap + superflousNrg*whole_step
      sumSenHeatCanopy     = sumSenHeatCanopy     - superflousNrg*whole_step
      scalarCanopyLiq      = 0._rkind
    endif

    ! update water
    scalarCanopyWat = scalarCanopyLiq + scalarCanopyIce
    if(enthalpyStateVec .or. computeEnthalpy)then ! recompute enthalpy of the canopy if changed water and ice content
      call T2enthTemp_veg(&
                  canopyDepth,    & ! intent(in): canopy depth (m), send in specific value since diag_data may have changed
                  specificHeatVeg,                                      & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                  maxMassVegetation,                                    & ! intent(in): maximum mass of vegetation (kg m-2)
                  snowfrz_scale,                                        & ! intent(in): scaling parameter for the snow freezing curve  (K-1)
                  scalarCanopyTemp,     & ! intent(in): canopy temperature (K)
                  scalarCanopyWat,                                      & ! intent(in): canopy water content (kg m-2)
                  scalarCanopyEnthTemp)   ! intent(out): temperature component of enthalpy of the vegetation canopy (J m-3)
      scalarCanopyEnthalpy = scalarCanopyEnthTemp - LH_fus * scalarCanopyIce/ canopyDepth
    endif
  end if  ! (if computing the vegetation flux)

  ! * compute change in ice content of the top snow layer due to sublimation 
  !   and account for compaction and cavitation in the snowpack...
  ! ------------------------------------------------------------------------
  call computSnowDepth(&
            whole_step,                               & ! intent(in)
            nSnow,                                    & ! intent(in)
            sumSnowSublimation/whole_step,            & ! intent(in)
            mLayerVolFracLiq,                         & ! intent(inout)
            mLayerVolFracIce,                         & ! intent(inout)
            mLayerTemp,  & ! intent(in)
            mLayerMeltFreeze,                         & ! intent(in)
            densScalGrowth,tempScalGrowth, &
            grainGrowthRate,densScalOvrbdn, &
            tempScalOvrbdn,baseViscosity, & ! intent(in)
            ! output
            tooMuchSublim,                            & ! intent(out): flag to denote that there was too much sublimation in a given time step
            mLayerDepth)                               ! intent(inout)

  call calcHeight_d(nLayers,layerType,&
            mLayerDepth,mLayerHeight,iLayerHeight)

  ! recompute snow depth, SWE, and layer water
  if(nSnow > 0)then
    scalarSnowDepth = sum( mLayerDepth(1:nSnow) )
    scalarSWE       = sum( (mLayerVolFracLiq(1:nSnow)*iden_water &
                                                      + mLayerVolFracIce(1:nSnow)*iden_ice) * mLayerDepth(1:nSnow) )
    mLayerVolFracWat(1:nSnow) = mLayerVolFracLiq(1:nSnow) + mLayerVolFracIce(1:nSnow)*iden_ice/iden_water
    if(enthalpyStateVec .or. computeEnthalpy)then ! recompute enthalpy of layers if changed water and ice content
      do iLayer=1,nSnow
        call T2enthTemp_snow(&
                     snowfrz_scale,                                       & ! intent(in):  scaling parameter for the snow freezing curve  (K-1)
                     mLayerTemp(iLayer),     & ! intent(in):  layer temperature (K)
                     mLayerVolFracWat(iLayer),                            & ! intent(in):  volumetric total water content (-)
                     mLayerEnthTemp(iLayer))   ! intent(out): temperature component of enthalpy of each snow layer (J m-3)
        mLayerEnthalpy(iLayer) = mLayerEnthTemp(iLayer) - iden_ice * LH_fus * mLayerVolFracIce(iLayer)
      end do  ! looping through snow layers
    endif
  endif

  ! increment change in storage associated with the surface melt pond (kg m-2)
  if(nSnow==0) sfcMeltPond = sfcMeltPond + scalarSfcMeltPond

end subroutine

attributes(global) subroutine increment_sums_kernel(nGRU, &
  flux_data, &
  sumCanopySublimation, sumSnowSublimation, &
  sumLatHeatCanopyEvap, sumSenHeatCanopy, &
  scalarCanopySublimation, scalarSnowSublimation, &
  scalarLatHeatCanopyEvap, scalarSenHeatCanopy, &
  dt_sub)
  integer(i4b),intent(in),value :: nGRU
  real(rkind),intent(inout) :: flux_data(:,:)
  real(rkind),intent(inout) :: sumCanopySublimation(:), sumSnowSublimation(:)
  real(rkind),intent(inout) :: sumLatHeatCanopyEvap(:), sumSenHeatCanopy(:)
  integer(i4b),intent(in),value :: scalarCanopySublimation, scalarSnowSublimation
  integer(i4b),intent(in),value :: scalarLatHeatCanopyEvap, scalarSenHeatCanopy
  real(rkind),intent(in),value :: dt_sub

  integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return
  call increment_sums(sumCanopySublimation(iGRU), sumSnowSublimation(iGRU), &
  sumLatHeatCanopyEvap(iGRU), sumSenHeatCanopy(iGRU), &
  flux_data(scalarCanopySublimation,iGRU), flux_data(scalarSnowSublimation,iGRU), &
  flux_data(scalarLatHeatCanopyEvap,iGRU), flux_data(scalarSenHeatCanopy,iGRU), &
  dt_sub)

end subroutine

attributes(device) subroutine increment_sums(sumCanopySublimation, sumSnowSublimation, &
  sumLatHeatCanopyEvap, sumSenHeatCanopy, &
  scalarCanopySublimation, scalarSnowSublimation, &
  scalarLatHeatCanopyEvap, scalarSenHeatCanopy, &
  dt_sub)
  real(rkind),intent(inout) :: sumCanopySublimation, sumSnowSublimation
  real(rkind),intent(inout) :: sumLatHeatCanopyEvap, sumSenHeatCanopy
  real(rkind),intent(inout) :: scalarCanopySublimation, scalarSnowSublimation
  real(rkind),intent(inout) :: scalarLatHeatCanopyEvap, scalarSenHeatCanopy
  real(rkind),intent(in) :: dt_sub
  sumCanopySublimation = sumCanopySublimation + dt_sub*scalarCanopySublimation
  sumSnowSublimation   = sumSnowSublimation   + dt_sub*scalarSnowSublimation
  sumLatHeatCanopyEvap = sumLatHeatCanopyEvap + dt_sub*scalarLatHeatCanopyEvap
  sumSenHeatCanopy     = sumSenHeatCanopy     + dt_sub*scalarSenHeatCanopy
end subroutine


end module coupled_em_module
