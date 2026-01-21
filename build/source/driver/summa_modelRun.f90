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

module summa_modelRun
! calls the model physics

! access missing values
USE globalData,only:integerMissing   ! missing integer
USE globalData,only:realMissing      ! missing real number

! named variables
USE globalData,only:yes,no           ! .true. and .false.
USE var_lookup,only:iLookTIME        ! named variables for time data structure
USE var_lookup,only:iLookDIAG        ! look-up values for local column model diagnostic variables
USE var_lookup,only:iLookINDEX       ! look-up values for local column index variables
USE summa_util,only:handle_err

! these are needed because we cannot access them in modules locally if we might use those modules with Actors
USE globalData,only:fracJulDay       ! fractional julian days since the start of year
USE globalData,only:yearLength       ! number of days in the current year

! safety: set private unless specified otherwise
implicit none
private
public::summa_runPhysics
contains

 ! calls the model physics
 subroutine summa_runPhysics(modelTimeStep, summa1_struc, err, message)
 ! ---------------------------------------------------------------------------------------
 ! * desired modules
 ! ---------------------------------------------------------------------------------------
 ! data types
 USE nrtype                                                     ! variable types, etc.
 USE summa_type, only:summa1_type_dec                           ! master summa data type
 ! subroutines and functions
 USE nr_utility_module,only:indexx                              ! sort vectors in ascending order
 USE vegPhenlgy_module,only:vegPhenlgy                          ! module to compute vegetation phenology
 USE run_oneGRU_module,only:run_oneGRU                          ! module to run for one GRU
 USE time_utils_module,only:elapsedSec                          ! calculate the elapsed time
 ! global data
 USE globalData,only:gru_struc                                  ! gru-hru mapping structures
 USE globalData,only:model_decisions                            ! model decision structure
 USE globalData,only:startPhysics,endPhysics                    ! date/time for the start and end of the initialization
 USE globalData,only:elapsedPhysics                             ! elapsed time for the initialization
 use initialize_device
 use device_data_types
 use globalData,only:minExpLogHgtFac,urbanVegCategory

 ! ---------------------------------------------------------------------------------------
 ! * variables
 ! ---------------------------------------------------------------------------------------
 implicit none
 ! dummy variables
 integer(i4b),intent(in)               :: modelTimeStep         ! time step index
 type(summa1_type_dec),intent(inout)   :: summa1_struc          ! master summa data structure
 integer(i4b),intent(out)              :: err                   ! error code
 character(*),intent(out)              :: message               ! error message
 ! ---------------------------------------------------------------------------------------
 ! local variables: general
 character(LEN=256)                    :: cmessage              ! error message of downwind routine
 integer(i4b)                          :: iHRU                  ! HRU index
 integer(i4b)                          :: iGRU,jGRU        ! GRU indices
 ! local variables: veg phenology
 logical(lgt)                          :: computeVegFluxFlag    ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(rkind)                           :: notUsed_canopyDepth   ! NOT USED: canopy depth (m)
 real(rkind)                           :: notUsed_exposedVAI    ! NOT USED: exposed vegetation area index (m2 m-2)
!  type(prog_data_device) :: progStruct_d
!  type(flux_data_device) :: fluxStruct_d
!  type(diag_data_device) :: diagStruct_d
!  type(bvar_data_device) :: bvarStruct_d
!  type(attr_data_device) :: attrStruct_d
!  type(mpar_data_device) :: mparStruct_d
!   type(decisions_device) :: decisions
! type(veg_parameters) :: veg_param
! type(veg_param_tables) :: veg_tables
! type(indx_data_device) :: indxStruct_d
! integer(i4b),device,allocatable :: computeVegFlux_d(:)
  type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(summa1_struc%nGRU/128+1,1,1)

 ! ---------------------------------------------------------------------------------------
 ! associate to elements in the data structure
 summaVars: associate(&

  ! primary data structures (scalars)
  timeStruct           => summa1_struc%timeStruct          , & ! x%var(:)                   -- model time data
  forcStruct           => summa1_struc%forcStruct_d          , & ! x%gru(:)%hru(:)%var(:)     -- model forcing data
  attrStruct           => summa1_struc%attrStruct_d          , & ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
  typeStruct           => summa1_struc%typeStruct          , & ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
  idStruct             => summa1_struc%idStruct            , & ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU

  veg_param => summa1_struc%veg_param, &
  decisions => summa1_struc%decisions, &
  ! primary data structures (variable length vectors)
  indxStruct           => summa1_struc%indxStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
  mparStruct           => summa1_struc%mparStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model parameters
  progStruct           => summa1_struc%progStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
  diagStruct           => summa1_struc%diagStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
  fluxStruct           => summa1_struc%fluxStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes

  ! basin-average structures
  bparStruct           => summa1_struc%bparStruct          , & ! x%gru(:)%var(:)            -- basin-average parameters
  bvarStruct           => summa1_struc%bvarStruct_d          , & ! x%gru(:)%var(:)%dat        -- basin-average variables

  ! run time variables
  computeVegFlux       => summa1_struc%computeVegFlux_d      , & ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  dt_init              => summa1_struc%dt_init             , & ! used to initialize the length of the sub-step for each HRU
  nGRU                 => summa1_struc%nGRU                  & ! number of grouped response units

 ) ! assignment to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='summa_runPhysics/'

 ! *******************************************************************************************
 ! *** initialize computeVegFlux (flag to indicate if we are computing fluxes over vegetation)
 ! *******************************************************************************************
 ! if computeVegFlux changes, then the number of state variables changes, and we need to reorganize the data structures
 if(modelTimeStep==1)then
    call initialize_computeVegFlux_kernel<<<blocks,threads>>>(nGRU,indxStruct%nSnow, decisions%bcUpprTdyn,decisions%bcUpprSoiH,&
fracJulDay,yearLength,&
typeStruct%vegTypeIndex,attrStruct%latitude,&
mparStruct%z0Snow_,mparStruct%z0Soil_,mparStruct%heightCanopyTop_,mparStruct%heightCanopyBottom_,&
progStruct%scalarSnowDepth,progStruct%scalarCanopyTemp,progStruct%scalarCanopyLiq,&
diagStruct%scalarRootZoneTemp,diagStruct%scalarLAI,diagStruct%scalarSAI,&
diagStruct%scalarExposedLAI,diagStruct%scalarExposedSAI,&
diagStruct%scalarGrowingSeasonIndex,diagStruct%scalarGroundSnowFraction, &
veg_param%dveg,veg_param%iswater,veg_param%isbarren,veg_param%issnow,&
veg_param%hvt_,veg_param%hvb_,veg_param%laim_,veg_param%saim_,veg_param%tmin,urbanVegCategory,minExpLogHgtFac, computeVegFlux)
 end if  ! if the first time step

 ! ****************************************************************************
 ! *** model simulation
 ! ****************************************************************************

 ! initialize the start of the physics
 call date_and_time(values=startPhysics)

 ! ----- rank the GRUs in terms of their anticipated computational expense -----

 ! estimate computational expense based on persistence
 !  -- assume that that expensive GRUs from a previous time step are also expensive in the current time step


 end associate summaVars

 ! ----- use openMP directives to run GRUs in parallel -------------------------
 ! start of parallel section: define shared and private structure elements
  !$omp parallel  default(none) &
  !$omp          private(iGRU, jGRU) &  ! GRU indices are private for a given thread
  !$omp          shared(summa1_struc, gru_struc) &
  !$omp          private(err, cmessage)
 ! associate to elements in the data structur, gru_struce
 ! need to associate again for the parallelism to work
 summaVars2: associate(&

  ! primary data structures (scalars)
  timeStruct           => summa1_struc%timeStruct          , & ! x%var(:)                   -- model time data
  forcStruct           => summa1_struc%forcStruct_d          , & ! x%gru(:)%hru(:)%var(:)     -- model forcing data
  attrStruct           => summa1_struc%attrStruct_d          , & ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
  typeStruct           => summa1_struc%typeStruct          , & ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
  idStruct             => summa1_struc%idStruct            , & ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU

  ! primary data structures (variable length vectors)
  indxStruct           => summa1_struc%indxStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
  mparStruct           => summa1_struc%mparStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model parameters
  progStruct           => summa1_struc%progStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
  diagStruct           => summa1_struc%diagStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
  fluxStruct           => summa1_struc%fluxStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes

  ! basin-average structures
  bparStruct           => summa1_struc%bparStruct          , & ! x%gru(:)%var(:)            -- basin-average parameters
  bvarStruct           => summa1_struc%bvarStruct_d          , & ! x%gru(:)%var(:)%dat        -- basin-average variables

  ! lookup table structure
  lookupStruct         => summa1_struc%lookupStruct        , & ! x%gru(:)%hru(:)%z(:)%var(:)%lookup    -- lookup-tables
  decisions => summa1_struc%decisions, &
  veg_param => summa1_struc%veg_param, &
  veg_tables => summa1_struc%veg_tables, &

  ! run time variables
  computeVegFlux       => summa1_struc%computeVegFlux_d      , & ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  dt_init              => summa1_struc%dt_init             , & ! used to initialize the length of the sub-step for each HRU
  nGRU                 => summa1_struc%nGRU                  & ! number of grouped response units

 ) ! assignment to variables in the data structures

  !----- run simulation for all GRU ----------------------------------------
  call run_oneGRU(&
  nGRU, &
                  ! model control
                  gru_struc,              & ! intent(inout): HRU information for given GRU (# HRUs, #snow+soil layers)
                  dt_init%gru(1)%hru,        & ! intent(inout): used to initialize the length of the sub-step for each HRU
                  computeVegFlux, & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                  ! data structures (input)
                  timeStruct%var,               & ! intent(in):    model time data
                  typeStruct,         & ! intent(in):    local classification of soil veg etc. for each HRU
                  idStruct%gru(1),           & ! intent(in):    local classification of soil veg etc. for each HRU
                  attrStruct,         & ! intent(in):    local attributes for each HRU
                  lookupStruct,       & ! intent(in):    local lookup tables for each HRU
                  decisions, veg_param, veg_tables, &
                  ! data structures (input-output)
                  mparStruct,         & ! intent(inout): local model parameters
                  indxStruct,         & ! intent(inout): model indices
                  forcStruct,         & ! intent(inout): model forcing data
                  progStruct,         & ! intent(inout): prognostic variables for a local HRU
                  diagStruct,         & ! intent(inout): diagnostic variables for a local HRU
                  fluxStruct,         & ! intent(inout): model fluxes for a local HRU
                  bvarStruct,         & ! intent(inout): basin-average variables
                  ! error control
                  err,cmessage)                   ! intent(out):   error control

  ! check errors
  call handle_err(err, cmessage)



  end associate summaVars2
 !$omp end parallel

 ! identify the end of the physics
 call date_and_time(values=endPhysics)

 ! aggregate the elapsed time for the physics
 elapsedPhysics = elapsedPhysics + elapsedSec(startPhysics, endPhysics)

 ! end associate statements

 end subroutine summa_runPhysics
 attributes(global) subroutine initialize_computeVegFlux_kernel(nGRU,nSnow, bcUpprTdyn,bcUpprSoiH,&
fracJulDay,yearLength,&
vegTypeIndex,latitude,&
z0Snow,z0Soil,heightCanopyTop,heightCanopyBottom,&
scalarSnowDepth,scalarCanopyTemp,scalarCanopyLiq,&
scalarRootZoneTemp,scalarLAI,scalarSAI,&
scalarExposedLAI,scalarExposedSAI,&
scalarGrowingSeasonIndex,scalarGroundSnowFraction, &
dveg,iswater,isbarren,issnow,&
hvt,hvb,laim,saim,tmin,urbanVegCategory,minExpLogHgtFac, computeVegFlux)
  USE nrtype                                                     ! variable types, etc.

  integer(i4b),intent(in),value :: nGRU
 integer(i4b),intent(in) :: nSnow(:), bcUpprTdyn,bcUpprSoiH
 real(rkind),intent(in),value :: fracJulDay
 integer(i4b),intent(in),value :: yearLength
 integer(i4b),intent(in) :: vegTypeIndex(:)
 real(rkind),intent(in) :: latitude(:)
 real(rkind),intent(in) :: z0Snow(:),z0Soil(:),heightCanopyTop(:),heightCanopyBottom(:)
 real(rkind),intent(inout) :: scalarSnowDepth(:),scalarCanopyTemp(:),scalarCanopyLiq(:)
 real(rkind),intent(inout) :: scalarRootZoneTemp(:),scalarLAI(:),scalarSAI(:)
 real(rkind),intent(inout) :: scalarExposedLAI(:),scalarExposedSAI(:)
 real(rkind),intent(inout) :: scalarGrowingSeasonIndex(:),scalarGroundSnowFraction(:)
 integer(i4b),intent(in) :: dveg,iswater,isbarren,issnow
 integer(i4b),intent(in),value :: urbanVegCategory
 real(rkind),intent(in) :: hvt(:,:),hvb(:,:),laim(:,:,:),saim(:,:,:),tmin(:)
 real(rkind),intent(in),value :: minExpLogHgtFac
integer(i4b),intent(inout) :: computeVegFlux(:)

     integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return
  call initialize_computeVegFlux(nSnow(iGRU), bcUpprTdyn,bcUpprSoiH,&
fracJulDay,yearLength,&
vegTypeIndex(iGRU),latitude(iGRU),&
z0Snow(iGRU),z0Soil(iGRU),heightCanopyTop(iGRU),heightCanopyBottom(iGRU),&
scalarSnowDepth(iGRU),scalarCanopyTemp(iGRU),scalarCanopyLiq(iGRU),&
scalarRootZoneTemp(iGRU),scalarLAI(iGRU),scalarSAI(iGRU),&
scalarExposedLAI(iGRU),scalarExposedSAI(iGRU),&
scalarGrowingSeasonIndex(iGRU),scalarGroundSnowFraction(iGRU), &
dveg,iswater,isbarren,issnow,&
hvt(:,iGRU),hvb(:,iGRU),laim(:,:,iGRU),saim(:,:,iGRU),tmin,urbanVegCategory,minExpLogHgtFac, computeVegFlux(iGRU))
print*, iGRU, scalarCanopyTemp(iGRU)
end subroutine

attributes(device) subroutine initialize_computeVegFlux(nSnow, bcUpprTdyn,bcUpprSoiH,&
fracJulDay,yearLength,&
vegTypeIndex,latitude,&
z0Snow,z0Soil,heightCanopyTop,heightCanopyBottom,&
scalarSnowDepth,scalarCanopyTemp,scalarCanopyLiq,&
scalarRootZoneTemp,scalarLAI,scalarSAI,&
scalarExposedLAI,scalarExposedSAI,&
scalarGrowingSeasonIndex,scalarGroundSnowFraction, &
dveg,iswater,isbarren,issnow,&
hvt,hvb,laim,saim,tmin,urbanVegCategory,minExpLogHgtFac, computeVegFlux)
 USE vegPhenlgy_module,only:vegPhenlgy_d                          ! module to compute vegetation phenology
 use nrtype

 integer(i4b),intent(in) :: nSnow, bcUpprTdyn,bcUpprSoiH
 real(rkind),intent(in) :: fracJulDay
 integer(i4b),intent(in) :: yearLength
 integer(i4b),intent(in) :: vegTypeIndex
 real(rkind),intent(in) :: latitude
 real(rkind),intent(in) :: z0Snow,z0Soil,heightCanopyTop,heightCanopyBottom
 real(rkind),intent(inout) :: scalarSnowDepth,scalarCanopyTemp,scalarCanopyLiq
 real(rkind),intent(inout) :: scalarRootZoneTemp,scalarLAI,scalarSAI
 real(rkind),intent(inout) :: scalarExposedLAI,scalarExposedSAI
 real(rkind),intent(inout) :: scalarGrowingSeasonIndex,scalarGroundSnowFraction
 integer(i4b),intent(in) :: dveg,iswater,isbarren,issnow,urbanVegCategory
 real(rkind),intent(in) :: hvt(:),hvb(:),laim(:,:),saim(:,:),tmin(:),minExpLogHgtFac
integer(i4b),intent(inout) :: computeVegFlux
 integer(i4b) :: err
 logical(lgt) :: computeVegFluxFlag
 real(rkind) :: notUsed_canopyDepth, notUsed_exposedVAI
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
                    scalarSnowDepth,scalarCanopyTemp,scalarCanopyLiq, &
                    ! diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                    scalarRootZoneTemp,scalarLAI,scalarSAI,&
                    scalarExposedLAI,scalarExposedSAI,&
                    scalarGrowingSeasonIndex,scalarGroundSnowFraction,&
                    ! output
                    computeVegFluxFlag,                  & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                    notUsed_canopyDepth,                 & ! intent(out): NOT USED: canopy depth (m)
                    notUsed_exposedVAI,                  & ! intent(out): NOT USED: exposed vegetation area index (m2 m-2)
                    err,&
                    dveg,iswater,isbarren,issnow,&
                    hvt,hvb,laim,saim,tmin,urbanVegCategory,minExpLogHgtFac)                  ! intent(out): error control
    ! save the flag for computing the vegetation fluxes
    if(computeVegFluxFlag)      computeVegFlux = yes
    if(.not.computeVegFluxFlag) computeVegFlux = no

end subroutine

end module summa_modelRun
