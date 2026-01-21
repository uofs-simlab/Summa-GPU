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
  nGRU, &
                       ! model control
                       hru_nc,              & ! intent(in):    hru index in netcdf
                       hruId,               & ! intent(in):    hruId
                       dt_init,             & ! intent(inout): used to initialize the length of the sub-step for each HRU
                       computeVegFlux,      & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                       ! data structures (input)
                       timeVec,             & ! intent(in):    model time data
                       typeData,            & ! intent(in):    local classification of soil veg etc. for each HRU
                       attrData,            & ! intent(in):    local attributes for each HRU
                       lookupData,          & ! intent(in):    local lookup tables for each HRU
                       bvarData,            & ! intent(in):    basin-average variables
                       decisions, veg_param, veg_tables, &
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

 USE derivforce_module,only:derivforce        ! module to compute derived forcing data
 USE coupled_em_module,only:coupled_em        ! module to run the coupled energy and mass model
 use initialize_device
 use device_data_types
 implicit none

 ! ----- define dummy variables ------------------------------------------------------------------------------------------

 ! model control
 integer(i4b)      , intent(in)    :: hru_nc              ! hru index in netcdf
 integer(i8b)      , intent(in)    :: hruId               ! hruId
 real(rkind)       , intent(inout) :: dt_init             ! used to initialize the length of the sub-step for each HRU
 logical(lgt),device      , intent(inout) :: computeVegFlux(:)      ! flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
 ! data structures (input)
 integer(i4b)      , intent(in)    :: timeVec(:)          ! int vector               -- model time data
 type(type_data_device)       , intent(in)    :: typeData            ! x%var(:)                 -- local classification of soil veg etc. for each HRU
 type(attr_data_device)       , intent(in)    :: attrData            ! x%var(:)                 -- local attributes for each HRU
 type(zLookup_device)     , intent(in)    :: lookupData          ! x%z(:)%var(:)%lookup(:)  -- local lookup tables for each HRU
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
 real(rkind)          , allocatable   :: zSoilReverseSign(:) ! height at bottom of each soil layer, negative downwards (m)
integer(i4b),intent(in) :: nGRU
type(decisions_device) :: decisions
type(veg_parameters) :: veg_param
type(veg_param_tables) :: veg_tables
real(rkind),device,allocatable :: greenVegFrac_monthly_d(:)
  type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

 ! initialize error control
 err=0; write(message, '(A21,I0,A10,I0,A2)' ) 'run_oneHRU (hru nc = ',hru_nc ,', hruId = ',hruId,')/'
 greenVegFrac_monthly_d = greenVegFrac_monthly

 ! ----- hru initialization ---------------------------------------------------------------------------------------------
 ! initialize the number of flux calls
 diagData%numFluxCalls = 0._rkind

 ! water pixel: do nothing
!  if (typeData%var(iLookTYPE%vegTypeIndex) == isWater)then
  ! Set wall_clock time to zero so it does not get a random value
   diagData%wallClockTime = 0._rkind
  !  return
!  endif


 ! ----- run the model --------------------------------------------------------------------------------------------------
call initialize_run_oneHRU<<<blocks,threads>>>(nGRU,urbanVegCategory,decisions%snowDenNew,decisions%LAI_method,&
typeData%vegTypeIndex,typeData%soilTypeIndex,typeData%slopeTypeIndex,&
fluxData%data,fluxData%ixScalarSnowfall,fluxData%ixScalarRainfall,fluxData%ixSpectralIncomingDirect_start,fluxData%ixSpectralIncomingDirect_end,fluxData%ixSpectralIncomingDiffuse_start,fluxData%ixSpectralIncomingDiffuse_end,&
diagData%scalarNewSnowDensity,diagData%scalarSnowfallTemp,diagData%scalarTwetbulb,diagData%scalarVPair,diagData%scalarFractionDirect,diagData%scalarCosZenith,&
diagData%windspd_x,diagData%windspd_y,diagData%scalarAdjMeasHeight,progData%scalarSnowDepth,diagData%scalarCO2air,diagData%scalarO2air,&
mparData%newSnowDenMin_,mparData%newSnowDenMult_,mparData%newSnowDenScal_,mparData%newSnowDenAdd_,mparData%newSnowDenMultTemp_, &
mparData%newSnowDenMultWind_,mparData%newSnowDenMultAnd_,mparData%newSnowDenBase_,mparData%constSnowDen_,mparData%snowfrz_scale_, &
mparData%tempRangeTimestep_,mparData%tempCritRain_,mparData%frozenPrecipMultip_,mparData%minwind_,mparData%Frad_vis_,mparData%Frad_direct_,mparData%directScale_,&
mparData%heightCanopyTop_,mparData%heightCanopyBottom_,&
forcData%airtemp_d,forcData%windspd_d,forcData%pptrate,forcData%spechum,forcData%airpres_d,forcData%SWRadAtm,forcData%time,&
attrData%longitude,attrData%latitude,attrData%tan_slope,attrData%aspect,attrData%mHeight,&
tmZoneOffsetFracDay,decisions%NC_TIME_ZONE,decisions%data_step,decisions%refJulDay,&
veg_param%RSMIN_d,mparData%minStomatalResistance_,&
veg_param%HVT_,veg_param%HVB_,&
veg_param%SAIM_,veg_param%LAIM_,mparData%winterSAI_,mparData%summerLAI_,greenVegFrac_monthly_d,&
  veg_tables%RSTBL,veg_tables%RGLTBL,veg_tables%HSTBL, &
  veg_tables%TOPT_DATA,veg_tables%RSMAX_DATA, &
  veg_param%RGL_d,veg_param%HS_d,veg_param%RSMAX_d,veg_param%TOPT_d)
 ! run the model for a single HRU
 call coupled_em(&
 nGRU, &
                 ! model control
                 hruId,            & ! intent(in):    hruId
                 dt_init,          & ! intent(inout): initial time step
                 1,                & ! intent(in):    used to adjust the length of the timestep with failure in Actors (non-Actors here, always 1)
                 computeVegFlux,   & ! intent(inout): flag to indicate if we are computing fluxes over vegetation
                 fracJulDay,       & ! intent(in):    fractional julian days since the start of year
                 yearLength,       & ! intent(in):    number of days in the current year
                 ! data structures (input)
                 typeData,         & ! intent(in):    local classification of soil veg etc. for each HRU
                 attrData,         & ! intent(in):    local attributes for each HRU
                 forcData,         & ! intent(in):    model forcing data
                 mparData,         & ! intent(in):    model parameters
                 bvarData,         & ! intent(in):    basin-average model variables
                 lookupData,       & ! intent(in):    lookup tables
                 decisions, veg_param, &
                 ! data structures (input-output)
                 indxData,         & ! intent(inout): model indices
                 progData,         & ! intent(inout): model prognostic variables for a local HRU
                 diagData,         & ! intent(inout): model diagnostic variables for a local HRU
                 fluxData,         & ! intent(inout): model fluxes for a local HRU
                 ! error control
                 err,cmessage)       ! intent(out):   error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif


 end subroutine run_oneHRU

 attributes(global) subroutine initialize_run_oneHRU(nGRU,urbanVegCategory,snowDenNew,LAI_method,&
 vegTypeIndex,soilTypeIndex,slopeTypeIndex,&
flux_data,snowfall,rainfall,spectralIncomingDirect_start,spectralIncomingDirect_end,spectralIncomingDiffuse_start,spectralIncomingDiffuse_end,&
newSnowDensity,snowfallTemp,twetbulb,VPair,scalarFractionDirect,cosZenith,&
windspd_x,windspd_y,adjMeasHeight,scalarSnowDepth,scalarCO2air,scalarO2air,&
newSnowDenMin,newSnowDenMult,newSnowDenScal,newSnowDenAdd,newSnowDenMultTemp, &
newSnowDenMultWind,newSnowDenMultAnd,newSnowDenBase,constSnowDen,snowfrz_scale, &
tempRangeTimestep,tempCritRain,frozenPrecipMultip,minwind,Frad_vis,Frad_direct,directScale, &
heightCanopyTop, heightCanopyBottom,&
airtemp,windspd,pptrate,spechum,airpres,SWRadAtm,time,&
longitude,latitude,tan_slope,aspect,mHeight,&
tmZoneOffsetFracDay,NC_TIME_ZONE,data_step,refJulDay,&
RSMIN,minStomatalResistance,&
HVT,HVB,&
SAIM,LAIM,winterSAI,summerLAI,greenVegFrac_monthly,&
  RSTBL,RGLTBL,HSTBL, &
  TOPT_DATA,RSMAX_DATA, &
  RGL,HS,RSMAX,TOPT)
implicit none
 integer(i4b),intent(in),value :: nGRU,urbanVegCategory
 real(rkind),intent(inout) :: flux_data(:,:)
 integer(i4b),intent(in),value :: snowfall,rainfall,spectralIncomingDirect_start,spectralIncomingDirect_end,spectralIncomingDiffuse_start,spectralIncomingDiffuse_end
integer(i4b),intent(in) :: snowDenNew,LAI_method
real(rkind),intent(inout) :: newSnowDensity(:),snowfallTemp(:),twetbulb(:),VPair(:),scalarFractionDirect(:),cosZenith(:)
real(rkind),intent(inout) :: windspd_x(:), windspd_y(:),adjMeasHeight(:),scalarSnowDepth(:),scalarCO2air(:),scalarO2air(:)
real(rkind),intent(inout) :: newSnowDenMin(:), newSnowDenMult(:),newSnowDenScal(:),newSnowDenAdd(:),newSnowDenMultTemp(:)
real(rkind),intent(inout) :: newSnowDenMultWind(:),newSnowDenMultAnd(:),newSnowDenBase(:),constSnowDen(:),snowfrz_scale(:)
real(rkind),intent(inout) :: tempRangeTimestep(:),tempCritRain(:),frozenPrecipMultip(:),minwind(:),Frad_vis(:),Frad_direct(:),directScale(:)
real(rkind),intent(inout) :: heightCanopyTop(:),heightCanopyBottom(:)
real(rkind),intent(inout) :: airtemp(:),windspd(:),pptrate(:),spechum(:),airpres(:),SWRadAtm(:),time(:)
real(rkind),intent(inout) :: longitude(:), latitude(:), tan_slope(:), aspect(:),mHeight(:)
real(rkind),intent(in),value :: tmZoneOffsetFracDay
integer(i4b) :: NC_TIME_ZONE
real(rkind) :: data_step,refJulDay
real(rkind),intent(inout) :: RSMIN(:), minStomatalResistance(:)
real(rkind),intent(inout) :: HVT(:,:),HVB(:,:)
real(rkind),intent(inout) :: SAIM(:,:,:), LAIM(:,:,:)
integer(i4b),intent(in) :: vegTypeIndex(:),soilTypeIndex(:),slopeTypeIndex(:)
real(rkind),intent(in) :: winterSAI(:),summerLAI(:),greenVegFrac_monthly(:)
    real(rkind),intent(inout) :: RSTBL(:),RGLTBL(:),HSTBL(:)
    real(rkind),intent(inout) :: TOPT_DATA,RSMAX_DATA
    real(rkind),intent(inout) :: rgl(:),hs(:),rsmax(:),topt(:)

   integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return
  call initialize_run_oneHRU_device(snowDenNew,LAI_method,&
  vegTypeIndex(iGRU),soilTypeIndex(iGRU),slopeTypeIndex(iGRU),urbanVegCategory,&
flux_data(snowfall,iGRU),flux_data(rainfall,iGRU),flux_data(spectralIncomingDirect_start:spectralIncomingDirect_end,iGRU),&
flux_data(spectralIncomingDiffuse_start:spectralIncomingDiffuse_end,iGRU),&
newSnowDensity(iGRU),snowfallTemp(iGRU),twetbulb(iGRU),VPair(iGRU),scalarFractionDirect(iGRU),cosZenith(iGRU),&
windspd_x(iGRU),windspd_y(iGRU),adjMeasHeight,scalarSnowDepth(iGRU),scalarCO2air(iGRU),scalarO2air(iGRU),&
newSnowDenMin(iGRU),newSnowDenMult(iGRU),newSnowDenScal(iGRU),newSnowDenAdd(iGRU),newSnowDenMultTemp(iGRU), &
newSnowDenMultWind(iGRU),newSnowDenMultAnd(iGRU),newSnowDenBase(iGRU),constSnowDen(iGRU),snowfrz_scale(iGRU), &
tempRangeTimestep(iGRU),tempCritRain(iGRU),frozenPrecipMultip(iGRU),minwind(iGRU),Frad_vis(iGRU),Frad_direct(iGRU),directScale(iGRU), &
heightCanopyTop(iGRU), heightCanopyBottom(iGRU),&
airtemp(iGRU),windspd(iGRU),pptrate(iGRU),spechum(iGRU),airpres(iGRU),SWRadAtm(iGRU),time(iGRU),&
longitude(iGRU),latitude(iGRU),tan_slope(iGRU),aspect(iGRU),mHeight(iGRU),&
tmZoneOffsetFracDay,NC_TIME_ZONE,data_step,refJulDay,&
RSMIN(iGRU),minStomatalResistance(iGRU),&
HVT(:,iGRU),HVB(:,iGRU),&
SAIM(:,:,iGRU),LAIM(:,:,iGRU),winterSAI(iGRU),summerLAI(iGRU),greenVegFrac_monthly,&
  RSTBL,RGLTBL,HSTBL, &
  TOPT_DATA,RSMAX_DATA, &
  RGL(iGRU),HS(iGRU),RSMAX(iGRU),TOPT(iGRU))
  if (iGRU .eq. 1) then
    do iGRU=1,nGRU
      print*, mHeight(iGRU), adjMeasHeight(iGRU), heightCanopyTop(iGRU),scalarSnowDepth(iGRU)
    end do
  end if

end subroutine

attributes(device) subroutine initialize_run_oneHRU_device(snowDenNew,LAI_method,&
vegTypeIndex,soilTypeIndex,slopeTypeIndex,urbanVegCategory,&
snowfall,rainfall,spectralIncomingDirect,spectralIncomingDiffuse,&
newSnowDensity,snowfallTemp,twetbulb,VPair,scalarFractionDirect,cosZenith,&
windspd_x,windspd_y,adjMeasHeight_,scalarSnowDepth,scalarCO2air,scalarO2air, &
newSnowDenMin,newSnowDenMult,newSnowDenScal,newSnowDenAdd,newSnowDenMultTemp, &
newSnowDenMultWind,newSnowDenMultAnd,newSnowDenBase,constSnowDen,snowfrz_scale, &
tempRangeTimestep,tempCritRain,frozenPrecipMultip,minwind,Frad_vis,Frad_direct,directScale,&
heightCanopyTop, heightCanopyBottom,&
airtemp,windspd,pptrate,spechum,airpres,SWRadAtm,time,&
longitude,latitude,tan_slope,aspect,mHeight,&
tmZoneOffsetFracDay,NC_TIME_ZONE,data_step,refJulDay,&
RSMIN,minStomatalResistance,&
HVT,HVB,&
SAIM,LAIM,winterSAI,summerLAI,greenVegFrac_monthly,&
  RSTBL,RGLTBL,HSTBL, &
  TOPT_DATA,RSMAX_DATA, &
  RGL,HS,RSMAX,TOPT)
use derivforce_module,only:derivforce_d
 USE module_sf_noahmplsm,only:redprm          ! module to assign more Noah-MP parameters

real(rkind),intent(inout) :: snowfall,rainfall,spectralIncomingDirect(:),spectralIncomingDiffuse(:)
integer(i4b),intent(in) :: snowDenNew,LAI_method
real(rkind),intent(inout) :: newSnowDensity,snowfallTemp,twetbulb,VPair,scalarFractionDirect,cosZenith
real(rkind),intent(inout) :: windspd_x,windspd_y,adjMeasHeight_(:),scalarSnowDepth,scalarCO2air,scalarO2air
real(rkind),intent(inout) :: newSnowDenMin, newSnowDenMult,newSnowDenScal,newSnowDenAdd,newSnowDenMultTemp
real(rkind),intent(inout) :: newSnowDenMultWind,newSnowDenMultAnd,newSnowDenBase,constSnowDen,snowfrz_scale
real(rkind),intent(inout) :: tempRangeTimestep,tempCritRain,frozenPrecipMultip,minwind,Frad_vis,Frad_direct,directScale
real(rkind),intent(in),value :: heightCanopyTop,heightCanopyBottom
real(rkind),intent(inout) :: airtemp,windspd,pptrate,spechum,airpres,SWRadAtm,time
real(rkind),intent(in),value :: longitude,latitude,tan_slope,aspect,mHeight
real(rkind),intent(inout) :: tmZoneOffsetFracDay
real(rkind),intent(in) :: data_step,refJulDay
integer(i4b),intent(in) :: NC_TIME_ZONE
real(rkind),intent(inout) :: RSMIN, minStomatalResistance
real(rkind),intent(inout) :: HVT(:),HVB(:)
real(rkind),intent(inout) :: SAIM(:,:), LAIM(:,:)
integer(i4b),intent(in) :: vegTypeIndex,soilTypeIndex,slopeTypeIndex
real(rkind),intent(in) :: winterSAI,summerLAI,greenVegFrac_monthly(:)
    real(rkind),intent(inout) :: RSTBL(:),RGLTBL(:),HSTBL(:)
    real(rkind),intent(inout) :: TOPT_DATA,RSMAX_DATA
    real(rkind),intent(inout) :: rgl,hs,rsmax,topt
    integer(i4b),intent(in) :: urbanVegCategory

 ! populate parameters in Noah-MP modules
 ! Passing a maxSoilLayer in order to pass the check for NROOT, that is done to avoid making any changes to Noah-MP code.
 !  --> NROOT from Noah-MP veg tables (as read here) is not used in SUMMA
call REDPRM (vegTypeIndex,soilTypeIndex,slopeTypeIndex,urbanVegCategory,&
  RSTBL,RGLTBL,HSTBL, &
  TOPT_DATA,RSMAX_DATA, &
  RGL,RSMIN,HS,RSMAX,TOPT)

 ! overwrite the minimum resistance
 if(overwriteRSMIN) RSMIN = minStomatalResistance

 ! overwrite the vegetation height
 HVT(vegTypeIndex) = heightCanopyTop
 HVB(vegTypeIndex) = heightCanopyBottom

 ! overwrite the tables for LAI and SAI
 if(LAI_method == specified)then
  SAIM(vegTypeIndex,:) = winterSAI
  LAIM(vegTypeIndex,:) = summerLAI*greenVegFrac_monthly
 end if

 ! ----- hru forcing ----------------------------------------------------------------------------------------------------

 ! compute derived forcing variables
call derivforce_d(snowDenNew,&
snowfall,rainfall,spectralIncomingDirect,spectralIncomingDiffuse,&
newSnowDensity,snowfallTemp,twetbulb,VPair,scalarFractionDirect,cosZenith,&
windspd_x,windspd_y,adjMeasHeight_,scalarSnowDepth,scalarCO2air,scalarO2air,&
newSnowDenMin,newSnowDenMult,newSnowDenScal,newSnowDenAdd,newSnowDenMultTemp, &
newSnowDenMultWind,newSnowDenMultAnd,newSnowDenBase,constSnowDen,snowfrz_scale, &
tempRangeTimestep,tempCritRain,frozenPrecipMultip,minwind,Frad_vis,Frad_direct,directScale,&
heightCanopyTop,&
airtemp,windspd,pptrate,spechum,airpres,SWRadAtm,time,&
longitude,latitude,tan_slope,aspect,mHeight,&
tmZoneOffsetFracDay,NC_TIME_ZONE,data_step,refJulDay)
end subroutine
end module run_oneHRU_module
