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

module run_oneHRU_module

! numerical recipes data types
USE nrtype

! data types
USE data_types,only:&
               var_i,                    & ! x%var(:)                (i4b)
               var_d,                    & ! x%var(:)                (rkind)
               var_ilength,              & ! x%var(:)%dat            (i4b)
               var_dlength,              & ! x%var(:)%dat            (rkind)
               zLookup                     ! x%z(:)%var(:)%lookup(:) (rkind)

! access vegetation data
USE globalData,only:greenVegFrac_monthly   ! fraction of green vegetation in each month (0-1)
USE globalData,only:overwriteRSMIN         ! flag to overwrite RSMIN
USE globalData,only:maxSoilLayers          ! Maximum Number of Soil Layers

! provide access to Noah-MP constants
USE module_sf_noahmplsm,only:isWater       ! parameter for water land cover type

! provide access to the named variables that describe elements of parameter structures
USE var_lookup,only:iLookTYPE              ! look-up values for classification of veg, soils etc.
USE var_lookup,only:iLookATTR              ! look-up values for local attributes
USE var_lookup,only:iLookPARAM             ! look-up values for local column model parameters

! provide access to the named variables that describe elements of variable structures
USE var_lookup,only:iLookPROG              ! look-up values for local column model prognostic (state) variables
USE var_lookup,only:iLookDIAG              ! look-up values for local column model diagnostic variables
USE var_lookup,only:iLookINDEX             ! look-up values for local column index variables

! provide access to model decisions
USE globalData,only:model_decisions        ! model decision structure
USE var_lookup,only:iLookDECISIONS         ! look-up values for model decisions

! these are needed because we cannot access them in modules locally if we might use those modules with Actors
USE globalData,only:fracJulDay             ! fractional julian days since the start of year
USE globalData,only:yearLength             ! number of days in the current year
USE globalData,only:tmZoneOffsetFracDay    ! time zone offset in fractional days

! provide access to the named variables that describe model decisions
USE mDecisions_module,only:        &       ! look-up values for LAI decisions
                      monthlyTable,&       ! LAI/SAI taken directly from a monthly table for different vegetation classes
                      specified            ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters

! ----- global variables that are modified ------------------------------------------------------------------------------------------

! Noah-MP parameters
USE NOAHMP_VEG_PARAMETERS,only:SAIM,LAIM   ! 2-d tables for stem area index and leaf area index (vegType,month)
USE NOAHMP_VEG_PARAMETERS,only:HVT,HVB     ! height at the top and bottom of vegetation (vegType)
USE noahmp_globals,only:RSMIN              ! minimum stomatal resistance (vegType)

! urban vegetation category (could be local)
USE globalData,only:urbanVegCategory       ! vegetation category for urban areas

implicit none
private
public::run_oneHRU

contains

 ! ************************************************************************************************
 ! public subroutine run_oneGRU: simulation for a single GRU
 ! ************************************************************************************************

 ! simulation for a single HRU
 subroutine run_oneHRU(&
                       ! model control
                       hru_nc,              & ! intent(in):    hru index in netcdf
                       hruId,               & ! intent(in):    hruId
                       dt_init,             & ! intent(inout): used to initialize the length of the sub-step for each HRU
                       computeVegFlux,      & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                       nSnow,nSoil,nLayers, & ! intent(inout): number of snow and soil layers
                       nGRU, &
                       decisions, veg_param, tables, &
                       ! data structures (input)
                       timeVec,             & ! intent(in):    model time data
                       typeData,            & ! intent(in):    local classification of soil veg etc. for each HRU
                       attrData,            & ! intent(in):    local attributes for each HRU
                       lookupData,          & ! intent(in):    local lookup tables for each HRU
                       bvarData,            & ! intent(in):    basin-average variables
                       ! data structures (input-output)
                       mparData,            & ! intent(inout): local model parameters
                       indxData,            & ! intent(inout): model indices
                       forcData,            & ! intent(inout): model forcing data
                       progData,            & ! intent(inout): prognostic variables for a local HRU
                       diagData,            & ! intent(inout): diagnostic variables for a local HRU
                       fluxData,            & ! intent(inout): model fluxes for a local HRU
                       ! error control
                       err,message)           ! intent(out):   error control

 ! ----- define downstream subroutines -----------------------------------------------------------------------------------

 USE module_sf_noahmplsm,only:redprm          ! module to assign more Noah-MP parameters
 USE derivforce_module,only:derivforce        ! module to compute derived forcing data
 USE coupled_em_module,only:coupled_em        ! module to run the coupled energy and mass model
 use device_data_types
 use initialize_device
 implicit none

 ! ----- define dummy variables ------------------------------------------------------------------------------------------

 ! model control
 integer(i4b)      , intent(in)    :: hru_nc              ! hru index in netcdf
 integer(i8b)      , intent(in)    :: hruId               ! hruId
 real(rkind)       , intent(inout) :: dt_init             ! used to initialize the length of the sub-step for each HRU
 logical(lgt)      , intent(inout),device :: computeVegFlux(:)      ! flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
 integer(i4b)      , intent(inout) :: nSnow,nSoil,nLayers ! number of snow and soil layers
 ! data structures (input)
 integer(i4b)      , intent(in)    :: timeVec(:)          ! int vector               -- model time data
 type(type_data_device)       , intent(in)    :: typeData            ! x%var(:)                 -- local classification of soil veg etc. for each HRU
 type(attr_data_device)       , intent(in)    :: attrData            ! x%var(:)                 -- local attributes for each HRU
 type(zLookup)     , intent(in)    :: lookupData          ! x%z(:)%var(:)%lookup(:)  -- local lookup tables for each HRU
 type(bvar_data_device) , intent(in)    :: bvarData            ! x%var(:)%dat -- basin-average variables
 ! data structures (input-output)
 type(mpar_data_device) , intent(inout) :: mparData            ! x%var(:)%dat -- local (HRU) model parameters
 type(indx_data_device) , intent(inout) :: indxData            ! x%var(:)%dat -- model indices
 type(forc_data_device)       , intent(inout) :: forcData            ! x%var(:)     -- model forcing data
 type(prog_data_device) , intent(inout) :: progData            ! x%var(:)%dat -- model prognostic (state) variables
 type(diag_data_device) , intent(inout) :: diagData            ! x%var(:)%dat -- model diagnostic variables
 type(flux_data_device) , intent(inout) :: fluxData            ! x%var(:)%dat -- model fluxes
 ! error control
 integer(i4b)      , intent(out)   :: err                 ! error code
 character(*)      , intent(out)   :: message             ! error message

 ! ----- define local variables ------------------------------------------------------------------------------------------

 ! local variables
 character(len=256)                   :: cmessage            ! error message
!  real(rkind)          , allocatable   :: zSoilReverseSign(:) ! height at bottom of each soil layer, negative downwards (m)

 real(rkind),device :: tmZoneOffsetFracDay_d(ngRU)
 integer(i4b) :: nGRU
 type(decisions_device) :: decisions
 type(veg_parameters) :: veg_param
 type(veg_param_tables) :: tables
!  integer(i4b) :: nSnow,nSoil,nLayers
 real(rkind),device :: greenVegFrac_monthly_d(size(greenVegFrac_monthly))

 integer(i4b) :: iGRU, iLayer

 greenVegFrac_monthly_d = greenVegFrac_monthly
 nSnow = maxval(indxData%nSnow)
 nSOil = indxData%nSoil
 nLayers = maxval(indxData%nLayers_d)

 tmZoneOffsetFracDay_d = tmZoneOffsetFracDay

 ! initialize error control
 err=0; write(message, '(A21,I0,A10,I0,A2)' ) 'run_oneHRU (hru nc = ',hru_nc ,', hruId = ',hruId,')/'

 ! ----- hru initialization ---------------------------------------------------------------------------------------------
 ! initialize the number of flux calls
 diagData%numFluxCalls = 0._rkind

 ! water pixel: do nothing
!  if (typeData%var(iLookTYPE%vegTypeIndex) == isWater)then
!   ! Set wall_clock time to zero so it does not get a random value
!    diagData%var(iLookDIAG%wallClockTime)%dat(1) = 0._rkind
!    return
!  endif

 ! get height at bottom of each soil layer, negative downwards (used in Noah MP)
!  allocate(zSoilReverseSign(nSoil),stat=err)
!  if(err/=0)then
!   message=trim(message)//'problem allocating space for zSoilReverseSign'
!   err=20; return
!  endif
!  zSoilReverseSign(:) = -progData%var(iLookPROG%iLayerHeight)%dat(nSnow+1:nLayers)

!  ! populate parameters in Noah-MP modules
!  ! Passing a maxSoilLayer in order to pass the check for NROOT, that is done to avoid making any changes to Noah-MP code.
!  !  --> NROOT from Noah-MP veg tables (as read here) is not used in SUMMA
!  call REDPRM(typeData%var(iLookTYPE%vegTypeIndex),      & ! vegetation type index
!              typeData%var(iLookTYPE%soilTypeIndex),     & ! soil type
!              typeData%var(iLookTYPE%slopeTypeIndex),    & ! slope type index
!              zSoilReverseSign,                          & ! * not used: height at bottom of each layer [NOTE: negative] (m)
!              maxSoilLayers,                             & ! number of soil layers
!              urbanVegCategory)                            ! vegetation category for urban areas

 ! deallocate height at bottom of each soil layer(used in Noah MP)
!  deallocate(zSoilReverseSign,stat=err)
!  if(err/=0)then
!   message=trim(message)//'problem deallocating space for zSoilReverseSign'
!   err=20; return
!  endif
call redprm_d(typeData%vegTypeIndex,urbanVegCategory,veg_param,nGRU,tables)

 ! overwrite the minimum resistance
 if(overwriteRSMIN) then
  associate(RSMIN => veg_param%RSMIN_d, minStomatalResistance => mparData%minStomatalResistance)
    !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    RSMIN(iGRU) = minStomatalResistance
  end do
  end associate
end if

 ! overwrite the vegetation height
 associate(HVT => veg_param%HVT, HVB => veg_param%hvb, vegTypeIndex => typeData%vegTypeIndex, &
  heightCanopyTop => mparData%heightCanopyTop, heightCanopyBottom => mparData%heightCanopyBottom)
  !$cuf kernel do(1) <<<*,*>>>
 do iGRU=1,nGRU
  HVT(vegTypeIndex(iGRU)) = heightCanopyTop
  HVB(vegTypeIndex(iGRU)) = heightCanopyBottom
 end do
 end associate

 ! overwrite the tables for LAI and SAI
 associate(SAIM => veg_param%saim, laim=>veg_param%laim, vegTypeIndex => typeData%vegTypeIndex, winterSAI => mparData%winterSAI, summerLAI => mparData%summerLAI)
 if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
  !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    do iLayer=1,size(SAIM,2)
  SAIM(vegTypeIndex(iGRU),iLayer) = winterSAI
  LAIM(vegTypeIndex(iGRU),iLayer) = summerLAI*greenVegFrac_monthly_d(iLayer)
    end do
  end do
 end if
 end associate

 ! ----- hru forcing ----------------------------------------------------------------------------------------------------

 ! compute derived forcing variables
 call derivforce(timeVec,            & ! vector of time information
  decisions, &
 nGRU, &
                 forcData,       & ! vector of model forcing data
                 attrData,       & ! vector of model attributes
                 mparData,           & ! data structure of model parameters
                 progData,           & ! data structure of model prognostic variables
                 diagData,           & ! data structure of model diagnostic variables
                 fluxData,           & ! data structure of model fluxes
                 tmZoneOffsetFracDay_d,& ! time zone offset in fractional days
                 err,cmessage)         ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! ----- run the model --------------------------------------------------------------------------------------------------

 ! run the model for a single HRU
 call coupled_em(&
                 ! model control
                 hruId,            & ! intent(in):    hruId
                 dt_init,          & ! intent(inout): initial time step
                 1,                & ! intent(in):    used to adjust the length of the timestep with failure in Actors (non-Actors here, always 1)
                 computeVegFlux,   & ! intent(inout): flag to indicate if we are computing fluxes over vegetation
                 fracJulDay,       & ! intent(in):    fractional julian days since the start of year
                 yearLength,       & ! intent(in):    number of days in the current year
                 nGRU, &
                 ! data structures (input)
                 typeData,         & ! intent(in):    local classification of soil veg etc. for each HRU
                 attrData,         & ! intent(in):    local attributes for each HRU
                 forcData,         & ! intent(in):    model forcing data
                 mparData,         & ! intent(in):    model parameters
                 bvarData,         & ! intent(in):    basin-average model variables
                 lookupData,       & ! intent(in):    lookup tables
                 veg_param, decisions, &
                 ! data structures (input-output)
                 indxData,         & ! intent(inout): model indices
                 progData,         & ! intent(inout): model prognostic variables for a local HRU
                 diagData,         & ! intent(inout): model diagnostic variables for a local HRU
                 fluxData,         & ! intent(inout): model fluxes for a local HRU
                 ! error control
                 err,cmessage)       ! intent(out):   error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! update the number of layers
 nSnow   = maxval(indxData%nSnow)     ! number of snow layers
 nSoil   = indxData%nSoil     ! number of soil layers
 nLayers = maxval(indxData%nLayers_d)   ! total number of layers
tmZoneOffsetFracDay = tmZoneOffsetFracDay_d(1)

 end subroutine run_oneHRU

 subroutine redprm_d(vegtyp,isurban, veg_param,nGRU, tables)
  use device_data_types
  use module_sf_noahlsm,only:TOPT_DATA,RSMAX_DATA
  integer(i4b),device :: vegtyp(:)
  integer(i4b) :: isurban,nGRU
  type(veg_parameters) :: veg_param
  type(veg_param_tables) :: tables

  integer(i4b) :: iGRU

!       IF (SOILTYP .gt. SLCATS) THEN
!        call wrf_message('SOILTYP must be less than SLCATS:')
!        write(message, '("SOILTYP = ", I6, ";    SLCATS = ", I6)') SOILTYP, SLCATS
!        call wrf_message(trim(message))
!        call wrf_error_fatal ('REDPRM: Error: too many input soil types')
!     END IF
!     IF (VEGTYP .gt. LUCATS) THEN
!        call wrf_message('VEGTYP must be less than LUCATS:')
!        write(message, '("VEGTYP = ", I6, ";    LUCATS = ", I6)') VEGTYP, LUCATS
!        call wrf_message(trim(message))
!        call wrf_error_fatal ('Error: too many input landuse types')
!     END IF


  !    write(*,*) FRZK, FRZX, KDT, SLOPE, SLOPETYP
! ----------------------------------------------------------------------
! SET-UP VEGETATION PARAMETERS
! ----------------------------------------------------------------------
    ! Six redprm_canres variables:
    veg_param%TOPT = TOPT_DATA
    veg_param%RSMAX = RSMAX_DATA

  associate(rgl => veg_param%rgl_d, rgltbl => tables%rgltbl, &
    rsmin => veg_param%rsmin_d, hs => veg_param%hs_d, &
    rstbl => tables%rstbl, hstbl => tables%hstbl, vegtype=>vegtyp )
    !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    RGL(iGRU) = RGLTBL (VEGTYPe(iGRU))
    RSMIN(iGRU) = RSTBL (VEGTYPe(iGRU))
    HS(iGRU) = HSTBL (VEGTYPe(iGRU))
  end do

  end associate

  associate(rsmin => veg_param%rsmin_d, vegtype => vegtyp)
    !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    if (vegtype(iGRU)==isurban) rsmin(iGRU) = 400
  end do
  end associate

 end subroutine

end module run_oneHRU_module
