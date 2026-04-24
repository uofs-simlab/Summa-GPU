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

! used to manage output statistics of the model and forcing variables
module output_stats
USE nrtype, realMissing=>nr_realMissing
USE nrtype, integerMissing=>nr_integerMissing
implicit none
private
public :: calcStats
public :: calc_stats_kernel_var,calcStats_d
contains

 ! ******************************************************************************************************
 ! public subroutine calcStats is called at every model timestep to update/store output statistics
 ! from model variables
 ! ******************************************************************************************************
 subroutine calcStats(stat,dat,meta,resetStats,finalizeStats,statCounter,err,message)
 USE nrtype
 USE data_types,only:extended_info,dlength,ilength  ! metadata structure type
 USE var_lookup,only:iLookVarType                   ! named variables for variable types
 USE var_lookup,only:iLookSTAT                      ! named variables for output statistics types
 implicit none

 ! input variables
 type(dlength) ,intent(inout)   :: stat(:)          ! statistics
 class(*)      ,intent(in)      :: dat(:)           ! data
 type(extended_info),intent(in) :: meta(:)          ! metadata
 logical(lgt)  ,intent(in)      :: resetStats(:)    ! vector of flags to reset statistics
 logical(lgt)  ,intent(in)      :: finalizeStats(:) ! vector of flags to reset statistics
 integer(i4b)  ,intent(in)      :: statCounter(:)   ! number of time steps in each output frequency

 ! output variables
 integer(i4b)  ,intent(out)     :: err              ! error code
 character(*)  ,intent(out)     :: message          ! error message

 ! internals
 character(256)                 :: cmessage         ! error message
 integer(i4b)                   :: iVar             ! index for varaiable loop
 integer(i4b)                   :: pVar             ! index into parent structure
 real(rkind)                       :: tdata            ! dummy for pulling info from dat structure

 ! initialize error control
 err=0; message='calcStats/'

 ! loop through variables
 do iVar = 1,size(meta)

  ! don't do anything if var is not requested
  if (.not.meta(iVar)%varDesire) cycle

  ! only treat stats of scalars - all others handled separately
  if (meta(iVar)%varType==iLookVarType%outstat) then

   ! index in parent structure
   pVar = meta(iVar)%ixParent

   ! extract data from the structures
   select type (dat)
    type is (real(rkind));  tdata = dat(pVar)
    class is (dlength) ; tdata = dat(pVar)%dat(1)
    class is (ilength) ; tdata = real(dat(pVar)%dat(1), kind(rkind))
    class default;err=20;message=trim(message)//'dat type not found';return
   end select

   ! calculate statistics
   if (trim(meta(iVar)%varName)=='time') then
    stat(iVar)%dat(iLookSTAT%inst) = tdata
   else
    call calc_stats(meta(iVar),stat(iVar),tdata,resetStats,finalizeStats,statCounter,err,cmessage)
   end if
   if(err/=0)then; message=trim(message)//trim(cmessage);return; end if

  end if  ! if calculating statistics
 end do  ! looping through variables

 return
 end subroutine calcStats

 subroutine calc_stats_bvar_iVar(iVar,nGRU,data_step,stat,dat,resetStats,finalizeStats,statCounter)
    use device_data_types
    use var_lookup,only:iLookBVAR
  integer(i4b) :: iVar
  integer(i4b) :: nGRU
  real(rkind) :: data_step
  type(bvar_stat_device) :: stat
  type(bvar_data_device) :: dat
  logical(lgt),device :: resetStats(:), finalizeStats(:)
  integer(i4b),device :: statCounter(:)
   type(dim3) :: blocks,threads

    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  select case(iVar)
    case(iLookBVAR%basin__totalArea); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_basin__totalArea,stat%stat_basin__totalArea,dat%basin__totalArea,resetStats,finalizeStats,statCounter)
    case(iLookBVAR%basin__SurfaceRunoff); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_basin__SurfaceRunoff,stat%stat_basin__SurfaceRunoff,dat%basin__SurfaceRunoff,resetStats,finalizeStats,statCounter)
    case(iLookBVAR%basin__ColumnOutflow); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_basin__ColumnOutflow,stat%stat_basin__ColumnOutflow,dat%basin__ColumnOutflow,resetStats,finalizeStats,statCounter)
    case(iLookBVAR%basin__AquiferStorage); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_basin__AquiferStorage,stat%stat_basin__AquiferStorage,dat%basin__AquiferStorage,resetStats,finalizeStats,statCounter)
    case(iLookBVAR%basin__AquiferRecharge); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_basin__AquiferRecharge,stat%stat_basin__AquiferRecharge,dat%basin__AquiferRecharge,resetStats,finalizeStats,statCounter)
    case(iLookBVAR%basin__AquiferBaseflow); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_basin__AquiferBaseflow,stat%stat_basin__AquiferBaseflow,dat%basin__AquiferBaseflow,resetStats,finalizeStats,statCounter)
    case(iLookBVAR%basin__AquiferTranspire); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_basin__AquiferTranspire,stat%stat_basin__AquiferTranspire,dat%basin__AquiferTranspire,resetStats,finalizeStats,statCounter)
    case(iLookBVAR%basin__TotalRunoff); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_basin__TotalRunoff,stat%stat_basin__TotalRunoff,dat%basin__TotalRunoff,resetStats,finalizeStats,statCounter)
    case(iLookBVAR%basin__SoilDrainage); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_basin__SoilDrainage,stat%stat_basin__SoilDrainage,dat%basin__SoilDrainage,resetStats,finalizeStats,statCounter)
    case(iLookBVAR%averageInstantRunoff); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_averageInstantRunoff,stat%stat_averageInstantRunoff,dat%averageInstantRunoff,resetStats,finalizeStats,statCounter)
    case(iLookBVAR%averageRoutedRunoff); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_averageRoutedRunoff,stat%stat_averageRoutedRunoff,dat%averageRoutedRunoff,resetStats,finalizeStats,statCounter)
    end select
  end subroutine

 subroutine calc_stats_forc_iVar(iVar,nGRU,data_step,stat,dat,resetStats,finalizeStats,statCounter)
    use device_data_types
    use var_lookup,only:iLookFORCE
  integer(i4b) :: iVar
  integer(i4b) :: nGRU
  real(rkind) :: data_step
  type(forc_stat_device) :: stat
  type(forc_data_device) :: dat
  logical(lgt),device :: resetStats(:), finalizeStats(:)
  integer(i4b),device :: statCounter(:)
   type(dim3) :: blocks,threads

    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  select case(iVar)
  case(iLookFORCE%time); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_time,stat%stat_time,dat%time_h,resetStats,finalizeStats,statCounter)
  case(iLookFORCE%pptrate); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_pptrate,stat%stat_pptrate,dat%pptrate,resetStats,finalizeStats,statCounter)
  case(iLookFORCE%airtemp); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_airtemp,stat%stat_airtemp,dat%airtemp,resetStats,finalizeStats,statCounter)
  case(iLookFORCE%spechum); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_spechum,stat%stat_spechum,dat%spechum,resetStats,finalizeStats,statCounter)
  case(iLookFORCE%windspd); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_windspd,stat%stat_windspd,dat%windspd,resetStats,finalizeStats,statCounter)
  case(iLookFORCE%SWRadAtm); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_SWRadAtm,stat%stat_SWRadAtm,dat%SWRadAtm,resetStats,finalizeStats,statCounter)
  case(iLookFORCE%LWRadAtm); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_LWRadAtm,stat%stat_LWRadAtm,dat%LWRadAtm,resetStats,finalizeStats,statCounter)
  case(iLookFORCE%airpres); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_airpres,stat%stat_airpres,dat%airpres,resetStats,finalizeStats,statCounter)
    end select
  end subroutine


 subroutine calc_stats_prog_iVar(iVar,nGRU,data_step,stat,dat,resetStats,finalizeStats,statCounter)
    use device_data_types
    use var_lookup,only:iLookPROG
  integer(i4b) :: iVar
  integer(i4b) :: nGRU
  real(rkind) :: data_step
  type(prog_stat_device) :: stat
  type(prog_data_device) :: dat
  logical(lgt),device :: resetStats(:), finalizeStats(:)
  integer(i4b),device :: statCounter(:)
   type(dim3) :: blocks,threads

    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  select case(iVar)
  case(iLookPROG%scalarCanopyIce); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyIce,stat%stat_scalarCanopyIce,dat%scalarCanopyIce,resetStats,finalizeStats,statCounter)
  case(iLookPROG%scalarCanopyLiq); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyLiq,stat%stat_scalarCanopyLiq,dat%scalarCanopyLiq,resetStats,finalizeStats,statCounter)
  case(iLookPROG%scalarCanopyWat); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyWat,stat%stat_scalarCanopyWat,dat%scalarCanopyWat,resetStats,finalizeStats,statCounter)
  case(iLookPROG%scalarCanairTemp); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanairTemp,stat%stat_scalarCanairTemp,dat%scalarCanairTemp,resetStats,finalizeStats,statCounter)
  case(iLookPROG%scalarCanopyTemp); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyTemp,stat%stat_scalarCanopyTemp,dat%scalarCanopyTemp,resetStats,finalizeStats,statCounter)
  case(iLookPROG%scalarSnowAlbedo); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSnowAlbedo,stat%stat_scalarSnowAlbedo,dat%scalarSnowAlbedo,resetStats,finalizeStats,statCounter)
  case(iLookPROG%scalarSnowDepth); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSnowDepth,stat%stat_scalarSnowDepth,dat%scalarSnowDepth,resetStats,finalizeStats,statCounter)
  case(iLookPROG%scalarSWE); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSWE,stat%stat_scalarSWE,dat%scalarSWE,resetStats,finalizeStats,statCounter)
  case(iLookPROG%scalarSfcMeltPond); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSfcMeltPond,stat%stat_scalarSfcMeltPond,dat%scalarSfcMeltPond,resetStats,finalizeStats,statCounter)
  case(iLookPROG%scalarCanairEnthalpy); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanairEnthalpy,stat%stat_scalarCanairEnthalpy,dat%scalarCanairEnthalpy,resetStats,finalizeStats,statCounter)
  case(iLookPROG%scalarCanopyEnthalpy); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyEnthalpy,stat%stat_scalarCanopyEnthalpy,dat%scalarCanopyEnthalpy,resetStats,finalizeStats,statCounter)
  case(iLookPROG%scalarAquiferStorage); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarAquiferStorage,stat%stat_scalarAquiferStorage,dat%scalarAquiferStorage,resetStats,finalizeStats,statCounter)
  case(iLookPROG%scalarSurfaceTemp); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSurfaceTemp,stat%stat_scalarSurfaceTemp,dat%scalarSurfaceTemp,resetStats,finalizeStats,statCounter)
    end select
  end subroutine

 subroutine calc_stats_diag_iVar(iVar,nGRU,data_step,stat,dat,resetStats,finalizeStats,statCounter)
    use device_data_types
    use var_lookup,only:iLookDIAG
    implicit none
  integer(i4b) :: iVar
  integer(i4b) :: nGRU
  real(rkind) :: data_step
  type(diag_stat_device) :: stat
  type(diag_data_device) :: dat
  logical(lgt),device :: resetStats(:), finalizeStats(:)
  integer(i4b),device :: statCounter(:)
   type(dim3) :: blocks,threads
   integer(i4b) :: err

    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)
!   print*, 'diag', iVar

  select case(iVar)
  case(iLookDIAG%scalarCanopyDepth); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyDepth,stat%stat_scalarCanopyDepth,dat%scalarCanopyDepth,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarBulkVolHeatCapVeg); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarBulkVolHeatCapVeg,stat%stat_scalarBulkVolHeatCapVeg,dat%scalarBulkVolHeatCapVeg,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarCanopyCm); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyCm,stat%stat_scalarCanopyCm,dat%scalarCanopyCm,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarCanopyEmissivity); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyEmissivity,stat%stat_scalarCanopyEmissivity,dat%scalarCanopyEmissivity,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarRootZoneTemp); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarRootZoneTemp,stat%stat_scalarRootZoneTemp,dat%scalarRootZoneTemp,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarLAI); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLAI,stat%stat_scalarLAI,dat%scalarLAI,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarSAI); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSAI,stat%stat_scalarSAI,dat%scalarSAI,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarExposedLAI); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarExposedLAI,stat%stat_scalarExposedLAI,dat%scalarExposedLAI,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarExposedSAI); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarExposedSAI,stat%stat_scalarExposedSAI,dat%scalarExposedSAI,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarAdjMeasHeight); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarAdjMeasHeight,stat%stat_scalarAdjMeasHeight,dat%scalarAdjMeasHeight,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarCanopyIceMax); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyIceMax,stat%stat_scalarCanopyIceMax,dat%scalarCanopyIceMax,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarCanopyLiqMax); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyLiqMax,stat%stat_scalarCanopyLiqMax,dat%scalarCanopyLiqMax,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarGrowingSeasonIndex); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarGrowingSeasonIndex,stat%stat_scalarGrowingSeasonIndex,dat%scalarGrowingSeasonIndex,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarVolHtCap_air); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarVolHtCap_air,stat%stat_scalarVolHtCap_air,dat%scalarVolHtCap_air,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarVolHtCap_ice); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarVolHtCap_ice,stat%stat_scalarVolHtCap_ice,dat%scalarVolHtCap_ice,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarVolHtCap_soil); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarVolHtCap_soil,stat%stat_scalarVolHtCap_soil,dat%scalarVolHtCap_soil,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarVolHtCap_water); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarVolHtCap_water,stat%stat_scalarVolHtCap_water,dat%scalarVolHtCap_water,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarLambda_drysoil); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLambda_drysoil,stat%stat_scalarLambda_drysoil,dat%scalarLambda_drysoil,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarLambda_wetsoil); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLambda_wetsoil,stat%stat_scalarLambda_wetsoil,dat%scalarLambda_wetsoil,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarCanopyEnthTemp); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyEnthTemp,stat%stat_scalarCanopyEnthTemp,dat%scalarCanopyEnthTemp,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarTotalSoilEnthalpy); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarTotalSoilEnthalpy,stat%stat_scalarTotalSoilEnthalpy,dat%scalarTotalSoilEnthalpy,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarTotalSnowEnthalpy); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarTotalSnowEnthalpy,stat%stat_scalarTotalSnowEnthalpy,dat%scalarTotalSnowEnthalpy,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarVPair); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarVPair,stat%stat_scalarVPair,dat%scalarVPair,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarVP_CanopyAir); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarVP_CanopyAir,stat%stat_scalarVP_CanopyAir,dat%scalarVP_CanopyAir,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarTwetbulb); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarTwetbulb,stat%stat_scalarTwetbulb,dat%scalarTwetbulb,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarSnowfallTemp); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSnowfallTemp,stat%stat_scalarSnowfallTemp,dat%scalarSnowfallTemp,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarNewSnowDensity); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarNewSnowDensity,stat%stat_scalarNewSnowDensity,dat%scalarNewSnowDensity,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarO2air); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarO2air,stat%stat_scalarO2air,dat%scalarO2air,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarCO2air); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCO2air,stat%stat_scalarCO2air,dat%scalarCO2air,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%windspd_x); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_windspd_x,stat%stat_windspd_x,dat%windspd_x,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%windspd_y); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_windspd_y,stat%stat_windspd_y,dat%windspd_y,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarCosZenith); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCosZenith,stat%stat_scalarCosZenith,dat%scalarCosZenith,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarFractionDirect); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarFractionDirect,stat%stat_scalarFractionDirect,dat%scalarFractionDirect,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarCanopySunlitFraction); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopySunlitFraction,stat%stat_scalarCanopySunlitFraction,dat%scalarCanopySunlitFraction,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarCanopySunlitLAI); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopySunlitLAI,stat%stat_scalarCanopySunlitLAI,dat%scalarCanopySunlitLAI,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarCanopyShadedLAI); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyShadedLAI,stat%stat_scalarCanopyShadedLAI,dat%scalarCanopyShadedLAI,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarGroundAlbedo); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarGroundAlbedo,stat%stat_scalarGroundAlbedo,dat%scalarGroundAlbedo,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarLatHeatSubVapCanopy); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLatHeatSubVapCanopy,stat%stat_scalarLatHeatSubVapCanopy,dat%scalarLatHeatSubVapCanopy,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarLatHeatSubVapGround); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLatHeatSubVapGround,stat%stat_scalarLatHeatSubVapGround,dat%scalarLatHeatSubVapGround,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarSatVP_CanopyTemp); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSatVP_CanopyTemp,stat%stat_scalarSatVP_CanopyTemp,dat%scalarSatVP_CanopyTemp,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarSatVP_GroundTemp); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSatVP_GroundTemp,stat%stat_scalarSatVP_GroundTemp,dat%scalarSatVP_GroundTemp,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarZ0Canopy); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarZ0Canopy,stat%stat_scalarZ0Canopy,dat%scalarZ0Canopy,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarWindReductionFactor); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarWindReductionFactor,stat%stat_scalarWindReductionFactor,dat%scalarWindReductionFactor,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarZeroPlaneDisplacement); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarZeroPlaneDisplacement,stat%stat_scalarZeroPlaneDisplacement,dat%scalarZeroPlaneDisplacement,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarRiBulkCanopy); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarRiBulkCanopy,stat%stat_scalarRiBulkCanopy,dat%scalarRiBulkCanopy,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarRiBulkGround); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarRiBulkGround,stat%stat_scalarRiBulkGround,dat%scalarRiBulkGround,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarCanopyStabilityCorrection); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyStabilityCorrection,stat%stat_scalarCanopyStabilityCorrection,dat%scalarCanopyStabilityCorrection,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarGroundStabilityCorrection); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarGroundStabilityCorrection,stat%stat_scalarGroundStabilityCorrection,dat%scalarGroundStabilityCorrection,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarIntercellularCO2Sunlit); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarIntercellularCO2Sunlit,stat%stat_scalarIntercellularCO2Sunlit,dat%scalarIntercellularCO2Sunlit,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarIntercellularCO2Shaded); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarIntercellularCO2Shaded,stat%stat_scalarIntercellularCO2Shaded,dat%scalarIntercellularCO2Shaded,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarTranspireLim); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarTranspireLim,stat%stat_scalarTranspireLim,dat%scalarTranspireLim,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarTranspireLimAqfr); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarTranspireLimAqfr,stat%stat_scalarTranspireLimAqfr,dat%scalarTranspireLimAqfr,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarFoliageNitrogenFactor); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarFoliageNitrogenFactor,stat%stat_scalarFoliageNitrogenFactor,dat%scalarFoliageNitrogenFactor,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarSoilRelHumidity); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSoilRelHumidity,stat%stat_scalarSoilRelHumidity,dat%scalarSoilRelHumidity,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarAquiferRootFrac); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarAquiferRootFrac,stat%stat_scalarAquiferRootFrac,dat%scalarAquiferRootFrac,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarFracLiqVeg); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarFracLiqVeg,stat%stat_scalarFracLiqVeg,dat%scalarFracLiqVeg,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarCanopyWetFraction); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyWetFraction,stat%stat_scalarCanopyWetFraction,dat%scalarCanopyWetFraction,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarSnowAge); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSnowAge,stat%stat_scalarSnowAge,dat%scalarSnowAge,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarGroundSnowFraction); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarGroundSnowFraction,stat%stat_scalarGroundSnowFraction,dat%scalarGroundSnowFraction,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarInfilArea); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarInfilArea,stat%stat_scalarInfilArea,dat%scalarInfilArea,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarFrozenArea); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarFrozenArea,stat%stat_scalarFrozenArea,dat%scalarFrozenArea,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarSoilControl); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSoilControl,stat%stat_scalarSoilControl,dat%scalarSoilControl,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarSoilCompress); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSoilCompress,stat%stat_scalarSoilCompress,dat%scalarSoilCompress,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarTotalSoilLiq); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarTotalSoilLiq,stat%stat_scalarTotalSoilLiq,dat%scalarTotalSoilLiq,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarTotalSoilIce); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarTotalSoilIce,stat%stat_scalarTotalSoilIce,dat%scalarTotalSoilIce,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarTotalSoilWat); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarTotalSoilWat,stat%stat_scalarTotalSoilWat,dat%scalarTotalSoilWat,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarKappa); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarKappa,stat%stat_scalarKappa,dat%scalarKappa,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%scalarVolLatHt_fus); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarVolLatHt_fus,stat%stat_scalarVolLatHt_fus,dat%scalarVolLatHt_fus,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%numFluxCalls); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_numFluxCalls,stat%stat_numFluxCalls,dat%numFluxCalls,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%wallClockTime); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_wallClockTime,stat%stat_wallClockTime,dat%wallClockTime,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%meanStepSize); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_meanStepSize,stat%stat_meanStepSize,dat%meanStepSize,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%balanceCasNrg); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_balanceCasNrg,stat%stat_balanceCasNrg,dat%balanceCasNrg,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%balanceVegNrg); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_balanceVegNrg,stat%stat_balanceVegNrg,dat%balanceVegNrg,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%balanceSnowNrg); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_balanceSnowNrg,stat%stat_balanceSnowNrg,dat%balanceSnowNrg,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%balanceSoilNrg); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_balanceSoilNrg,stat%stat_balanceSoilNrg,dat%balanceSoilNrg,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%balanceVegMass); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_balanceVegMass,stat%stat_balanceVegMass,dat%balanceVegMass,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%balanceSnowMass); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_balanceSnowMass,stat%stat_balanceSnowMass,dat%balanceSnowMass,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%balanceSoilMass); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_balanceSoilMass,stat%stat_balanceSoilMass,dat%balanceSoilMass,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%balanceAqMass); call calc_stats_kernel_var<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_balanceAqMass,stat%stat_balanceAqMass,dat%balanceAqMass,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%numSteps); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_numSteps,stat%stat_numSteps,dat%numSteps,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%numResEvals); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_numResEvals,stat%stat_numResEvals,dat%numResEvals,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%numLinSolvSetups); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_numLinSolvSetups,stat%stat_numLinSolvSetups,dat%numLinSolvSetups,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%numErrTestFails); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_numErrTestFails,stat%stat_numErrTestFails,dat%numErrTestFails,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%kLast); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_kLast,stat%stat_kLast,dat%kLast,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%kCur); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_kCur,stat%stat_kCur,dat%kCur,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%hInitUsed); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_hInitUsed,stat%stat_hInitUsed,dat%hInitUsed,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%hLast); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_hLast,stat%stat_hLast,dat%hLast,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%hCur); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_hCur,stat%stat_hCur,dat%hCur,resetStats,finalizeStats,statCounter)
  case(iLookDIAG%tCur); call calc_stats_kernel_var_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_tCur,stat%stat_tCur,dat%tCur,resetStats,finalizeStats,statCounter)
    end select
  end subroutine

 subroutine calc_stats_indx_iVar(iVar,nGRU,data_step,stat,dat,resetStats,finalizeStats,statCounter)
    use device_data_types
    use var_lookup,only:iLookINDEX
  integer(i4b) :: iVar
  integer(i4b) :: nGRU
  real(rkind) :: data_step
  type(indx_stat_device) :: stat
  type(indx_data_device) :: dat
  logical(lgt),device :: resetStats(:), finalizeStats(:)
  integer(i4b),device :: statCounter(:)
   type(dim3) :: blocks,threads

    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  select case(iVar)
  case(iLookINDEX%nSnow); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_nSnow,stat%stat_nSnow,dat%nSnow,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%nSoil); call calc_stats_kernel_var_i_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_nSoil,stat%stat_nSoil,dat%nSoil,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%nLayers); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_nLayers,stat%stat_nLayers,dat%nLayers,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%nCasNrg); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_nCasNrg,stat%stat_nCasNrg,dat%nCasNrg,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%nVegNrg); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_nVegNrg,stat%stat_nVegNrg,dat%nVegNrg,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%nVegMass); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_nVegMass,stat%stat_nVegMass,dat%nVegMass,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%nVegState); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_nVegState,stat%stat_nVegState,dat%nVegState,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%nNrgState); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_nNrgState,stat%stat_nNrgState,dat%nNrgState,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%nWatState); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_nWatState,stat%stat_nWatState,dat%nWatState,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%nMatState); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_nMatState,stat%stat_nMatState,dat%nMatState,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%nMassState); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_nMassState,stat%stat_nMassState,dat%nMassState,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%nState); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_nState,stat%stat_nState,dat%nState,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%ixCasNrg); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_ixCasNrg,stat%stat_ixCasNrg,dat%ixCasNrg,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%ixVegNrg); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_ixVegNrg,stat%stat_ixVegNrg,dat%ixVegNrg,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%ixVegHyd); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_ixVegHyd,stat%stat_ixVegHyd,dat%ixVegHyd,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%ixTopNrg); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_ixTopNrg,stat%stat_ixTopNrg,dat%ixTopNrg,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%ixTopHyd); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_ixTopHyd,stat%stat_ixTopHyd,dat%ixTopHyd,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%ixAqWat); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_ixAqWat,stat%stat_ixAqWat,dat%ixAqWat,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%ixNrgCanair); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_ixNrgCanair,stat%stat_ixNrgCanair,dat%ixNrgCanair,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%ixNrgCanopy); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_ixNrgCanopy,stat%stat_ixNrgCanopy,dat%ixNrgCanopy,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%ixHydCanopy); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_ixHydCanopy,stat%stat_ixHydCanopy,dat%ixHydCanopy,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%ixWatAquifer); call calc_stats_kernel_var_i<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_ixWatAquifer,stat%stat_ixWatAquifer,dat%ixWatAquifer,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%numberFluxCalc); call calc_stats_kernel_var_i_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_numberFluxCalc,stat%stat_numberFluxCalc,dat%numberFluxCalc,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%numberStateSplit); call calc_stats_kernel_var_i_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_numberStateSplit,stat%stat_numberStateSplit,dat%numberStateSplit,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%numberDomainSplitNrg); call calc_stats_kernel_var_i_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_numberDomainSplitNrg,stat%stat_numberDomainSplitNrg,dat%numberDomainSplitNrg,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%numberDomainSplitMass); call calc_stats_kernel_var_i_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_numberDomainSplitMass,stat%stat_numberDomainSplitMass,dat%numberDomainSplitMass,resetStats,finalizeStats,statCounter)
  case(iLookINDEX%numberScalarSolutions); call calc_stats_kernel_var_i_host<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_numberScalarSolutions,stat%stat_numberScalarSolutions,dat%numberScalarSolutions,resetStats,finalizeStats,statCounter)
    end select
  end subroutine


 subroutine calc_stats_flux_iVar(iVar,nGRU,data_step,stat,dat,resetStats,finalizeStats,statCounter)
    use device_data_types
    use var_lookup,only:iLookFLUX
  integer(i4b) :: iVar
  integer(i4b) :: nGRU
  real(rkind) :: data_step
  type(flux_stat_device) :: stat
  type(flux_data_device) :: dat
  logical(lgt),device :: resetStats(:), finalizeStats(:)
  integer(i4b),device :: statCounter(:)
   type(dim3) :: blocks,threads

    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  select case(iVar)
  case(iLookFLUX%scalarCanairNetNrgFlux); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanairNetNrgFlux,stat%stat_scalarCanairNetNrgFlux,dat%ixscalarCanairNetNrgFlux,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopyNetNrgFlux); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyNetNrgFlux,stat%stat_scalarCanopyNetNrgFlux,dat%ixscalarCanopyNetNrgFlux,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarGroundNetNrgFlux); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarGroundNetNrgFlux,stat%stat_scalarGroundNetNrgFlux,dat%ixscalarGroundNetNrgFlux,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopyNetLiqFlux); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyNetLiqFlux,stat%stat_scalarCanopyNetLiqFlux,dat%ixscalarCanopyNetLiqFlux,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarRainfall); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarRainfall,stat%stat_scalarRainfall,dat%ixscalarRainfall,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarSnowfall); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSnowfall,stat%stat_scalarSnowfall,dat%ixscalarSnowfall,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopySunlitPAR); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopySunlitPAR,stat%stat_scalarCanopySunlitPAR,dat%ixscalarCanopySunlitPAR,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopyShadedPAR); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyShadedPAR,stat%stat_scalarCanopyShadedPAR,dat%ixscalarCanopyShadedPAR,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarBelowCanopySolar); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarBelowCanopySolar,stat%stat_scalarBelowCanopySolar,dat%ixscalarBelowCanopySolar,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopyAbsorbedSolar); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyAbsorbedSolar,stat%stat_scalarCanopyAbsorbedSolar,dat%ixscalarCanopyAbsorbedSolar,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarGroundAbsorbedSolar); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarGroundAbsorbedSolar,stat%stat_scalarGroundAbsorbedSolar,dat%ixscalarGroundAbsorbedSolar,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWRadCanopy); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWRadCanopy,stat%stat_scalarLWRadCanopy,dat%ixscalarLWRadCanopy,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWRadGround); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWRadGround,stat%stat_scalarLWRadGround,dat%ixscalarLWRadGround,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWRadUbound2Canopy); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWRadUbound2Canopy,stat%stat_scalarLWRadUbound2Canopy,dat%ixscalarLWRadUbound2Canopy,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWRadUbound2Ground); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWRadUbound2Ground,stat%stat_scalarLWRadUbound2Ground,dat%ixscalarLWRadUbound2Ground,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWRadUbound2Ubound); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWRadUbound2Ubound,stat%stat_scalarLWRadUbound2Ubound,dat%ixscalarLWRadUbound2Ubound,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWRadCanopy2Ubound); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWRadCanopy2Ubound,stat%stat_scalarLWRadCanopy2Ubound,dat%ixscalarLWRadCanopy2Ubound,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWRadCanopy2Ground); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWRadCanopy2Ground,stat%stat_scalarLWRadCanopy2Ground,dat%ixscalarLWRadCanopy2Ground,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWRadCanopy2Canopy); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWRadCanopy2Canopy,stat%stat_scalarLWRadCanopy2Canopy,dat%ixscalarLWRadCanopy2Canopy,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWRadGround2Ubound); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWRadGround2Ubound,stat%stat_scalarLWRadGround2Ubound,dat%ixscalarLWRadGround2Ubound,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWRadGround2Canopy); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWRadGround2Canopy,stat%stat_scalarLWRadGround2Canopy,dat%ixscalarLWRadGround2Canopy,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWNetCanopy); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWNetCanopy,stat%stat_scalarLWNetCanopy,dat%ixscalarLWNetCanopy,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWNetGround); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWNetGround,stat%stat_scalarLWNetGround,dat%ixscalarLWNetGround,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLWNetUbound); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLWNetUbound,stat%stat_scalarLWNetUbound,dat%ixscalarLWNetUbound,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarEddyDiffusCanopyTop); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarEddyDiffusCanopyTop,stat%stat_scalarEddyDiffusCanopyTop,dat%ixscalarEddyDiffusCanopyTop,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarFrictionVelocity); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarFrictionVelocity,stat%stat_scalarFrictionVelocity,dat%ixscalarFrictionVelocity,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarWindspdCanopyTop); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarWindspdCanopyTop,stat%stat_scalarWindspdCanopyTop,dat%ixscalarWindspdCanopyTop,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarWindspdCanopyBottom); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarWindspdCanopyBottom,stat%stat_scalarWindspdCanopyBottom,dat%ixscalarWindspdCanopyBottom,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarGroundResistance); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarGroundResistance,stat%stat_scalarGroundResistance,dat%ixscalarGroundResistance,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopyResistance); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyResistance,stat%stat_scalarCanopyResistance,dat%ixscalarCanopyResistance,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLeafResistance); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLeafResistance,stat%stat_scalarLeafResistance,dat%ixscalarLeafResistance,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarSoilResistance); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSoilResistance,stat%stat_scalarSoilResistance,dat%ixscalarSoilResistance,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarSenHeatTotal); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSenHeatTotal,stat%stat_scalarSenHeatTotal,dat%ixscalarSenHeatTotal,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarSenHeatCanopy); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSenHeatCanopy,stat%stat_scalarSenHeatCanopy,dat%ixscalarSenHeatCanopy,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarSenHeatGround); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSenHeatGround,stat%stat_scalarSenHeatGround,dat%ixscalarSenHeatGround,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLatHeatTotal); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLatHeatTotal,stat%stat_scalarLatHeatTotal,dat%ixscalarLatHeatTotal,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLatHeatCanopyEvap); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLatHeatCanopyEvap,stat%stat_scalarLatHeatCanopyEvap,dat%ixscalarLatHeatCanopyEvap,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLatHeatCanopyTrans); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLatHeatCanopyTrans,stat%stat_scalarLatHeatCanopyTrans,dat%ixscalarLatHeatCanopyTrans,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarLatHeatGround); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarLatHeatGround,stat%stat_scalarLatHeatGround,dat%ixscalarLatHeatGround,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopyAdvectiveHeatFlux); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyAdvectiveHeatFlux,stat%stat_scalarCanopyAdvectiveHeatFlux,dat%ixscalarCanopyAdvectiveHeatFlux,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarGroundAdvectiveHeatFlux); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarGroundAdvectiveHeatFlux,stat%stat_scalarGroundAdvectiveHeatFlux,dat%ixscalarGroundAdvectiveHeatFlux,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopySublimation); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopySublimation,stat%stat_scalarCanopySublimation,dat%ixscalarCanopySublimation,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarSnowSublimation); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSnowSublimation,stat%stat_scalarSnowSublimation,dat%ixscalarSnowSublimation,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarStomResistSunlit); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarStomResistSunlit,stat%stat_scalarStomResistSunlit,dat%ixscalarStomResistSunlit,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarStomResistShaded); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarStomResistShaded,stat%stat_scalarStomResistShaded,dat%ixscalarStomResistShaded,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarPhotosynthesisSunlit); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarPhotosynthesisSunlit,stat%stat_scalarPhotosynthesisSunlit,dat%ixscalarPhotosynthesisSunlit,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarPhotosynthesisShaded); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarPhotosynthesisShaded,stat%stat_scalarPhotosynthesisShaded,dat%ixscalarPhotosynthesisShaded,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopyTranspiration); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyTranspiration,stat%stat_scalarCanopyTranspiration,dat%ixscalarCanopyTranspiration,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopyEvaporation); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyEvaporation,stat%stat_scalarCanopyEvaporation,dat%ixscalarCanopyEvaporation,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarGroundEvaporation); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarGroundEvaporation,stat%stat_scalarGroundEvaporation,dat%ixscalarGroundEvaporation,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarThroughfallSnow); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarThroughfallSnow,stat%stat_scalarThroughfallSnow,dat%ixscalarThroughfallSnow,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarThroughfallRain); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarThroughfallRain,stat%stat_scalarThroughfallRain,dat%ixscalarThroughfallRain,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopySnowUnloading); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopySnowUnloading,stat%stat_scalarCanopySnowUnloading,dat%ixscalarCanopySnowUnloading,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopyLiqDrainage); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyLiqDrainage,stat%stat_scalarCanopyLiqDrainage,dat%ixscalarCanopyLiqDrainage,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarCanopyMeltFreeze); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarCanopyMeltFreeze,stat%stat_scalarCanopyMeltFreeze,dat%ixscalarCanopyMeltFreeze,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarSnowDrainage); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSnowDrainage,stat%stat_scalarSnowDrainage,dat%ixscalarSnowDrainage,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarRainPlusMelt); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarRainPlusMelt,stat%stat_scalarRainPlusMelt,dat%ixscalarRainPlusMelt,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarMaxInfilRate); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarMaxInfilRate,stat%stat_scalarMaxInfilRate,dat%ixscalarMaxInfilRate,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarInfiltration); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarInfiltration,stat%stat_scalarInfiltration,dat%ixscalarInfiltration,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarExfiltration); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarExfiltration,stat%stat_scalarExfiltration,dat%ixscalarExfiltration,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarSurfaceRunoff); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSurfaceRunoff,stat%stat_scalarSurfaceRunoff,dat%ixscalarSurfaceRunoff,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarSurfaceRunoff_IE); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSurfaceRunoff_IE,stat%stat_scalarSurfaceRunoff_IE,dat%ixscalarSurfaceRunoff_IE,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarSurfaceRunoff_SE); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSurfaceRunoff_SE,stat%stat_scalarSurfaceRunoff_SE,dat%ixscalarSurfaceRunoff_SE,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarSoilBaseflow); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSoilBaseflow,stat%stat_scalarSoilBaseflow,dat%ixscalarSoilBaseflow,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarSoilDrainage); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarSoilDrainage,stat%stat_scalarSoilDrainage,dat%ixscalarSoilDrainage,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarAquiferRecharge); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarAquiferRecharge,stat%stat_scalarAquiferRecharge,dat%ixscalarAquiferRecharge,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarAquiferTranspire); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarAquiferTranspire,stat%stat_scalarAquiferTranspire,dat%ixscalarAquiferTranspire,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarAquiferBaseflow); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarAquiferBaseflow,stat%stat_scalarAquiferBaseflow,dat%ixscalarAquiferBaseflow,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarTotalET); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarTotalET,stat%stat_scalarTotalET,dat%ixscalarTotalET,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarTotalRunoff); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarTotalRunoff,stat%stat_scalarTotalRunoff,dat%ixscalarTotalRunoff,dat%data,resetStats,finalizeStats,statCounter)
  case(iLookFLUX%scalarNetRadiation); call calc_stats_kernel_var_flux<<<blocks,threads>>>(nGRU,data_step,stat%statIndex_scalarNetRadiation,stat%stat_scalarNetRadiation,dat%ixscalarNetRadiation,dat%data,resetStats,finalizeStats,statCounter)
    end select
  end subroutine


 subroutine calcStats_d(nGRU,stat,dat,meta,resetStats,finalizeStats,statCounter,err,message)
 USE nrtype
 USE data_types,only:extended_info,dlength,ilength  ! metadata structure type
 USE var_lookup,only:iLookVarType                   ! named variables for variable types
 USE var_lookup,only:iLookSTAT                      ! named variables for output statistics types
     use device_data_types
 USE globalData,only:data_step        ! forcing timestep

 implicit none

 ! input variables
 integer(i4b),intent(in) :: nGRU
 class(*) ,intent(inout)   :: stat          ! statistics
 class(*)      ,intent(in)      :: dat           ! data
 type(extended_info),intent(in) :: meta(:)          ! metadata
 logical(lgt)  ,intent(in),device      :: resetStats(:)    ! vector of flags to reset statistics
 logical(lgt)  ,intent(in),device      :: finalizeStats(:) ! vector of flags to reset statistics
 integer(i4b)  ,intent(in),device      :: statCounter(:)   ! number of time steps in each output frequency

 ! output variables
 integer(i4b)  ,intent(out)     :: err              ! error code
 character(*)  ,intent(out)     :: message          ! error message

 ! internals
 character(256)                 :: cmessage         ! error message
 integer(i4b)                   :: iVar             ! index for varaiable loop
 integer(i4b)                   :: pVar             ! index into parent structure
 real(rkind)                       :: tdata            ! dummy for pulling info from dat structure

 ! initialize error control
 err=0; message='calcStats/'

 ! loop through variables
 do iVar = 1,size(meta)

  ! don't do anything if var is not requested
  if (.not.meta(iVar)%varDesire) cycle

  ! only treat stats of scalars - all others handled separately
  if (meta(iVar)%varType==iLookVarType%outstat) then

   ! index in parent structure
   pVar = meta(iVar)%ixParent
   
   ! extract data from the structures
      ! calculate statistics

   select type (dat)
    class is (bvar_data_device);
      select type(stat)
       class is (bvar_stat_device); call calc_stats_bvar_iVar(pVar,nGRU,data_step,stat,dat,resetStats,finalizeStats,statCounter)
       class default; err=20;message=trim(message)//'mismatch of dat and stat types[bvar_data_device]';return
      end select
    class is (prog_data_device);
      select type(stat)
       class is (prog_stat_device); call calc_stats_prog_iVar(pVar,nGRU,data_step,stat,dat,resetStats,finalizeStats,statCounter)
       class default; err=20;message=trim(message)//'mismatch of dat and stat types[prog_data_device]';return
      end select
    class is (diag_data_device);
      select type(stat)
       class is (diag_stat_device); call calc_stats_diag_iVar(pVar,nGRU,data_step,stat,dat,resetStats,finalizeStats,statCounter)
       class default; err=20;message=trim(message)//'mismatch of dat and stat types[diag_data_device]';return
      end select
    class is (forc_data_device);
      select type(stat)
       class is (forc_stat_device); call calc_stats_forc_iVar(pVar,nGRU,data_step,stat,dat,resetStats,finalizeStats,statCounter)
       class default; err=20;message=trim(message)//'mismatch of dat and stat types[forc_data_device]';return
      end select
    class is (flux_data_device);
      select type(stat)
       class is (flux_stat_device); call calc_stats_flux_iVar(pVar,nGRU,data_step,stat,dat,resetStats,finalizeStats,statCounter)
       class default; err=20;message=trim(message)//'mismatch of dat and stat types[flux_data_device]';return
      end select
    class is (indx_data_device);
      select type(stat)
       class is (indx_stat_device); call calc_stats_indx_iVar(pVar,nGRU,data_step,stat,dat,resetStats,finalizeStats,statCounter)
       class default; err=20;message=trim(message)//'mismatch of dat and stat types[indx_data_device]';return
      end select

    class default;err=20;message=trim(message)//'dat type not found';return
   end select

   err = cudaDeviceSynchronize()

  end if  ! if calculating statistics
 end do  ! looping through variables

 return
 end subroutine calcStats_d


 ! ***********************************************************************************
 ! Private subroutine calc_stats is a generic function to deal with any variable type.
 ! Called from compile_stats
 ! ***********************************************************************************
 subroutine calc_stats(meta,stat,tdata,resetStats,finalizeStats,statCounter,err,message)
 USE nrtype
 ! data structures
 USE data_types,only:var_info,ilength,dlength ! type dec for meta data structures
 USE var_lookup,only:maxVarFreq       ! # of output frequencies
 ! global variables
 USE globalData,only:data_step        ! forcing timestep
 ! structures of named variables
 USE var_lookup,only:iLookVarType     ! named variables for variable types
 USE var_lookup,only:iLookFREQ        ! named variables for output frequency
 USE var_lookup,only:iLookSTAT        ! named variables for output statistics
 USE var_lookup,only:iLookTIME        ! named variables for time information
 implicit none
 ! input variables
 class(var_info),intent(in)         :: meta              ! meta data structure
 class(*)       ,intent(inout)      :: stat              ! statistics structure
 real(rkind)       ,intent(in)         :: tdata             ! data value
 logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency
 ! output variables
 integer(i4b)   ,intent(out)        :: err               ! error code
 character(*)   ,intent(out)        :: message           ! error message
 ! internals
 real(rkind),dimension(maxvarFreq*2)   :: tstat             ! temporary stats vector
 integer(i4b)                       :: iFreq             ! index of output frequency
 ! initialize error control
 err=0; message='calc_stats/'

 ! extract variable from the data structure
 select type (stat)
  class is (ilength); tstat = real(stat%dat)
  class is (dlength); tstat = stat%dat
  class default;err=20;message=trim(message)//'stat type not found';return
 end select

 ! ---------------------------------------------
 ! reset statistics at new frequency period
 ! ---------------------------------------------
 do iFreq=1,maxVarFreq                              ! loop through output statistics
  if(resetStats(iFreq))then                         ! flag to reset statistics
   if(meta%statIndex(iFreq)==integerMissing) cycle  ! don't bother if output frequency is not desired for a given variable
   if(meta%varType/=iLookVarType%outstat) cycle     ! only calculate stats for scalars
   select case(meta%statIndex(iFreq))               ! act depending on the statistic
    ! -------------------------------------------------------------------------------------
    case (iLookSTAT%totl)                           ! * summation over period                  
     tstat(iFreq) = 0._rkind                           !     - resets stat at beginning of period
    case (iLookSTAT%mean)                           ! * mean over period                       
     tstat(iFreq) = 0._rkind                           !     - resets stat at beginning of period
    case (iLookSTAT%vari)                           ! * variance over period                   
     tstat(iFreq) = 0._rkind                           !     - resets E[X^2] term in var calc    
     tstat(maxVarFreq+iFreq) = 0._rkind                !     - resets E[X]^2 term                 
    case (iLookSTAT%mini)                           ! * minimum over period                    
     tstat(iFreq) = huge(tstat(iFreq))              !     - resets stat at beginning of period 
    case (iLookSTAT%maxi)                           ! * maximum over period                    
     tstat(iFreq) = -huge(tstat(iFreq))             !     - resets stat at beginning of period 
    case (iLookSTAT%mode)                           ! * mode over period      
     tstat(iFreq) = realMissing                     !     - does not work
    case (iLookSTAT%inst)                           ! * instantaneous -- no need to reset
    case default
     message=trim(message)//'unable to identify type of statistic [reset]'
     err=20; return
    ! -------------------------------------------------------------------------------------
   end select
  end if
 end do ! looping through output frequencies

 ! ---------------------------------------------
 ! Calculate each statistic that is requested by user
 ! ---------------------------------------------
 do iFreq=1,maxVarFreq                                ! loop through output statistics
  if(meta%statIndex(iFreq)==integerMissing) cycle     ! don't bother if output frequency is not desired for a given variab;e
  if(meta%varType/=iLookVarType%outstat) cycle        ! only calculate stats for scalars
  select case(meta%statIndex(iFreq))                  ! act depending on the statistic
   ! -------------------------------------------------------------------------------------
   case (iLookSTAT%inst)                              ! * instantaneous value 
    tstat(iFreq) = tdata                              !     - data at a given time
   case (iLookSTAT%totl)                              ! * summation over period                    
    tstat(iFreq) = tstat(iFreq) + tdata*data_step     !     - increment data 
   case (iLookSTAT%mean)                              ! * mean over period                       
    tstat(iFreq) = tstat(iFreq) + tdata               !     -  increment data
   case (iLookSTAT%vari)                              ! * variance over period                   
    tstat(iFreq) = tstat(iFreq) + tdata**2                     ! - E[X^2] term in var calc    
    tstat(maxVarFreq+iFreq) = tstat(maxVarFreq+iFreq) + tdata  ! - E[X]^2 term                 
   case (iLookSTAT%mini)                              ! * minimum over period                    
    if (tdata<tstat(iFreq)) tstat(iFreq) = tdata      !     - check value 
   case (iLookSTAT%maxi)                              ! * maximum over period                    
    if (tdata>tstat(iFreq)) tstat(iFreq) = tdata      !     - check value 
   case (iLookSTAT%mode)                              ! * mode over period (does not workind)       
    tstat(iFreq) = realMissing
   case default
    message=trim(message)//'unable to identify type of statistic [calculating stats]'
    err=20; return
   ! -------------------------------------------------------------------------------------
  end select
 end do ! looping through output frequencies

 ! ---------------------------------------------
 ! finalize statistics at end of frequency period
 ! ---------------------------------------------
 do iFreq=1,maxVarFreq                                ! loop through output statistics
  if(finalizeStats(iFreq))then
   if(meta%statIndex(iFreq)==integerMissing) cycle     ! don't bother if output frequency is not desired for a given variable
   if(meta%varType/=iLookVarType%outstat) cycle        ! only calculate stats for scalars
   select case(meta%statIndex(iFreq))                  ! act depending on the statistic
    ! -------------------------------------------------------------------------------------
    case (iLookSTAT%mean)                              ! * mean over period
     tstat(iFreq) = tstat(iFreq)/statCounter(iFreq)    !     - normalize sum into mean
    case (iLookSTAT%vari)                              ! * variance over period
     tstat(maxVarFreq+iFreq) = tstat(maxVarFreq+1)/statCounter(iFreq)            ! E[X] term
     tstat(iFreq) = tstat(iFreq)/statCounter(iFreq) - tstat(maxVarFreq+iFreq)**2 ! full variance
    case default ! do nothing -- don't need finalization for most stats
    ! -------------------------------------------------------------------------------------
   end select
  end if
 end do ! looping through output frequencies

 ! pack back into struc
 select type (stat)
  class is (ilength); stat%dat = int(tstat)
  class is (dlength); stat%dat = tstat
  class default;err=20;message=trim(message)//'stat type not found';return
 end select

 end subroutine calc_stats

 attributes(global) subroutine calc_stats_kernel_var(nGRU,&
 data_step, &
 statIndex,&
 stat,&
 var, &
resetStats,finalizeStats,statCounter)
! use initialize_device,only:get_iGRU
    integer(i4b),intent(in),value :: nGRU
    real(rkind),intent(in),value :: data_step
    integer(i4b) :: statIndex(:)
    real(rkind) :: stat(:,:)
    real(rkind) :: var(:)
     logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency
 integer(i4b) :: iGRU, err
 iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
 if (iGRU .gt. nGRU) return
    call calc_stats_d_real(data_step,statIndex,stat(:,iGRU),var(iGRU),resetStats,finalizeStats,statCounter,err)

end subroutine

 attributes(global) subroutine calc_stats_kernel_var_i(nGRU,&
 data_step, &
 statIndex,&
 stat,&
 var, &
resetStats,finalizeStats,statCounter)
! use initialize_device,only:get_iGRU
    integer(i4b),intent(in),value :: nGRU
    real(rkind),intent(in),value :: data_step
    integer(i4b) :: statIndex(:)
    real(rkind) :: stat(:,:)
    integer(i4b) :: var(:)
     logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency
 integer(i4b) :: iGRU, err
 iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
 if (iGRU .gt. nGRU) return
    call calc_stats_d_i4b(data_step,statIndex,stat(:,iGRU),var(iGRU),resetStats,finalizeStats,statCounter,err)

end subroutine

 attributes(global) subroutine calc_stats_kernel_var_i_host(nGRU,&
 data_step, &
 statIndex,&
 stat,&
 var, &
resetStats,finalizeStats,statCounter)
! use initialize_device,only:get_iGRU
    integer(i4b),intent(in),value :: nGRU
    real(rkind),intent(in),value :: data_step
    integer(i4b) :: statIndex(:)
    real(rkind) :: stat(:,:)
    integer(i4b),value :: var
     logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency
 integer(i4b) :: iGRU, err
 iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
 if (iGRU .gt. nGRU) return
    call calc_stats_d_i4b(data_step,statIndex,stat(:,iGRU),var,resetStats,finalizeStats,statCounter,err)

end subroutine

attributes(global) subroutine calc_stats_kernel_var_flux(nGRU,&
 data_step, &
 statIndex,&
 stat,&
 var, dat, &
resetStats,finalizeStats,statCounter)
! use initialize_device,only:get_iGRU
    integer(i4b),intent(in),value :: nGRU
    real(rkind),intent(in),value :: data_step
    integer(i4b) :: statIndex(:)
    real(rkind) :: stat(:,:)
    integer(i4b),value :: var
    real(rkind) :: dat(:,:)
     logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency
 integer(i4b) :: iGRU, err
 iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
 if (iGRU .gt. nGRU) return
    call calc_stats_d_real(data_step,statIndex,stat(:,iGRU),dat(var,iGRU),resetStats,finalizeStats,statCounter,err)

end subroutine


 attributes(global) subroutine calc_stats_kernel_var_host(nGRU,&
 data_step, &
 statIndex,&
 stat,&
 var, &
resetStats,finalizeStats,statCounter)
! use initialize_device,only:get_iGRU
    integer(i4b),intent(in),value :: nGRU
    real(rkind),intent(in),value :: data_step
    integer(i4b) :: statIndex(:)
    real(rkind) :: stat(:,:)
    real(rkind),value :: var
     logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency
 integer(i4b) :: iGRU, err
 iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
 if (iGRU .gt. nGRU) return
    call calc_stats_d_real(data_step,statIndex,stat(:,iGRU),var,resetStats,finalizeStats,statCounter,err)

end subroutine


 attributes(global) subroutine calc_stats_kernel_indx(nGRU, data_step, &
 statIndex_nSnow,stat_nSnow,nSnow, &
statIndex_nSoil,stat_nSoil,nSoil, &
statIndex_nLayers,stat_nLayers,nLayers, &
statIndex_nCasNrg,stat_nCasNrg,nCasNrg, &
statIndex_nVegNrg,stat_nVegNrg,nVegNrg, &
statIndex_nVegMass,stat_nVegMass,nVegMass, &
statIndex_nVegState,stat_nVegState,nVegState, &
statIndex_nNrgState,stat_nNrgState,nNrgState, &
statIndex_nWatState,stat_nWatState,nWatState, &
statIndex_nMatState,stat_nMatState,nMatState, &
statIndex_nMassState,stat_nMassState,nMassState, &
statIndex_nState,stat_nState,nState, &
statIndex_ixCasNrg,stat_ixCasNrg,ixCasNrg, &
statIndex_ixVegNrg,stat_ixVegNrg,ixVegNrg, &
statIndex_ixVegHyd,stat_ixVegHyd,ixVegHyd, &
statIndex_ixTopNrg,stat_ixTopNrg,ixTopNrg, &
statIndex_ixTopHyd,stat_ixTopHyd,ixTopHyd, &
statIndex_ixAqWat,stat_ixAqWat,ixAqWat, &
statIndex_ixNrgCanair,stat_ixNrgCanair,ixNrgCanair, &
statIndex_ixNrgCanopy,stat_ixNrgCanopy,ixNrgCanopy, &
statIndex_ixHydCanopy,stat_ixHydCanopy,ixHydCanopy, &
statIndex_ixWatAquifer,stat_ixWatAquifer,ixWatAquifer, &
statIndex_numberFluxCalc,stat_numberFluxCalc,numberFluxCalc, &
statIndex_numberStateSplit,stat_numberStateSplit,numberStateSplit, &
statIndex_numberDomainSplitNrg,stat_numberDomainSplitNrg,numberDomainSplitNrg, &
statIndex_numberDomainSplitMass,stat_numberDomainSplitMass,numberDomainSplitMass, &
statIndex_numberScalarSolutions,stat_numberScalarSolutions,numberScalarSolutions, &
resetStats,finalizeStats,statCounter)
implicit none
integer(i4b),intent(in),value :: nGRU
    real(rkind),intent(in),value :: data_step
    integer(i4b) :: statIndex_nSnow(:)
    real(rkind) :: stat_nSnow(:,:)
    integer(i4b) :: nSnow(:)
integer(i4b) :: statIndex_nSoil(:)
real(rkind) :: stat_nSoil(:,:)
integer(i4b),intent(in),value :: nSoil
integer(i4b) :: statIndex_nLayers(:)
real(rkind) :: stat_nLayers(:,:)
integer(i4b) :: nLayers(:)
integer(i4b) :: statIndex_nCasNrg(:)
real(rkind) :: stat_nCasNrg(:,:)
integer(i4b) :: nCasNrg(:)
integer(i4b) :: statIndex_nVegNrg(:)
real(rkind) :: stat_nVegNrg(:,:)
integer(i4b) :: nVegNrg(:)
integer(i4b) :: statIndex_nVegMass(:)
real(rkind) :: stat_nVegMass(:,:)
integer(i4b) :: nVegMass(:)
integer(i4b) :: statIndex_nVegState(:)
real(rkind) :: stat_nVegState(:,:)
integer(i4b) :: nVegState(:)
integer(i4b) :: statIndex_nNrgState(:)
real(rkind) :: stat_nNrgState(:,:)
integer(i4b) :: nNrgState(:)
integer(i4b) :: statIndex_nWatState(:)
real(rkind) :: stat_nWatState(:,:)
integer(i4b) :: nWatState(:)
integer(i4b) :: statIndex_nMatState(:)
real(rkind) :: stat_nMatState(:,:)
integer(i4b) :: nMatState(:)
integer(i4b) :: statIndex_nMassState(:)
real(rkind) :: stat_nMassState(:,:)
integer(i4b) :: nMassState(:)
integer(i4b) :: statIndex_nState(:)
real(rkind) :: stat_nState(:,:)
integer(i4b) :: nState(:)
integer(i4b) :: statIndex_ixCasNrg(:)
real(rkind) :: stat_ixCasNrg(:,:)
integer(i4b) :: ixCasNrg(:)
integer(i4b) :: statIndex_ixVegNrg(:)
real(rkind) :: stat_ixVegNrg(:,:)
integer(i4b) :: ixVegNrg(:)
integer(i4b) :: statIndex_ixVegHyd(:)
real(rkind) :: stat_ixVegHyd(:,:)
integer(i4b) :: ixVegHyd(:)
integer(i4b) :: statIndex_ixTopNrg(:)
real(rkind) :: stat_ixTopNrg(:,:)
integer(i4b) :: ixTopNrg(:)
integer(i4b) :: statIndex_ixTopHyd(:)
real(rkind) :: stat_ixTopHyd(:,:)
integer(i4b) :: ixTopHyd(:)
integer(i4b) :: statIndex_ixAqWat(:)
real(rkind) :: stat_ixAqWat(:,:)
integer(i4b) :: ixAqWat(:)
integer(i4b) :: statIndex_ixNrgCanair(:)
real(rkind) :: stat_ixNrgCanair(:,:)
integer(i4b) :: ixNrgCanair(:)
integer(i4b) :: statIndex_ixNrgCanopy(:)
real(rkind) :: stat_ixNrgCanopy(:,:)
integer(i4b) :: ixNrgCanopy(:)
integer(i4b) :: statIndex_ixHydCanopy(:)
real(rkind) :: stat_ixHydCanopy(:,:)
integer(i4b) :: ixHydCanopy(:)
integer(i4b) :: statIndex_ixWatAquifer(:)
real(rkind) :: stat_ixWatAquifer(:,:)
integer(i4b) :: ixWatAquifer(:)
integer(i4b) :: statIndex_numberFluxCalc(:)
real(rkind) :: stat_numberFluxCalc(:,:)
integer(i4b),intent(in),value :: numberFluxCalc
integer(i4b) :: statIndex_numberStateSplit(:)
real(rkind) :: stat_numberStateSplit(:,:)
integer(i4b),intent(in),value :: numberStateSplit
integer(i4b) :: statIndex_numberDomainSplitNrg(:)
real(rkind) :: stat_numberDomainSplitNrg(:,:)
integer(i4b),intent(in),value :: numberDomainSplitNrg
integer(i4b) :: statIndex_numberDomainSplitMass(:)
real(rkind) :: stat_numberDomainSplitMass(:,:)
integer(i4b),intent(in),value :: numberDomainSplitMass
integer(i4b) :: statIndex_numberScalarSolutions(:)
real(rkind) :: stat_numberScalarSolutions(:,:)
integer(i4b),intent(in),value :: numberScalarSolutions

     logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency
 integer(i4b) :: iGRU, err
 iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
 if (iGRU .gt. nGRU) return

call calc_stats_d_i4b(data_step, statIndex_nSnow, stat_nSnow(:,iGRU), nSnow(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_nSoil, stat_nSoil(:,iGRU), nSoil, resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_nLayers, stat_nLayers(:,iGRU), nLayers(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_nCasNrg, stat_nCasNrg(:,iGRU), nCasNrg(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_nVegNrg, stat_nVegNrg(:,iGRU), nVegNrg(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_nVegMass, stat_nVegMass(:,iGRU), nVegMass(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_nVegState, stat_nVegState(:,iGRU), nVegState(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_nNrgState, stat_nNrgState(:,iGRU), nNrgState(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_nWatState, stat_nWatState(:,iGRU), nWatState(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_nMatState, stat_nMatState(:,iGRU), nMatState(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_nMassState, stat_nMassState(:,iGRU), nMassState(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_nState, stat_nState(:,iGRU), nState(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_ixCasNrg, stat_ixCasNrg(:,iGRU), ixCasNrg(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_ixVegNrg, stat_ixVegNrg(:,iGRU), ixVegNrg(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_ixVegHyd, stat_ixVegHyd(:,iGRU), ixVegHyd(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_ixTopNrg, stat_ixTopNrg(:,iGRU), ixTopNrg(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_ixTopHyd, stat_ixTopHyd(:,iGRU), ixTopHyd(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_ixAqWat, stat_ixAqWat(:,iGRU), ixAqWat(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_ixNrgCanair, stat_ixNrgCanair(:,iGRU), ixNrgCanair(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_ixNrgCanopy, stat_ixNrgCanopy(:,iGRU), ixNrgCanopy(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_ixHydCanopy, stat_ixHydCanopy(:,iGRU), ixHydCanopy(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_ixWatAquifer, stat_ixWatAquifer(:,iGRU), ixWatAquifer(iGRU), resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_numberFluxCalc, stat_numberFluxCalc(:,iGRU), numberFluxCalc, resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_numberStateSplit, stat_numberStateSplit(:,iGRU), numberStateSplit, resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_numberDomainSplitNrg, stat_numberDomainSplitNrg(:,iGRU), numberDomainSplitNrg, resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_numberDomainSplitMass, stat_numberDomainSplitMass(:,iGRU), numberDomainSplitMass, resetStats, finalizeStats, statCounter, err)
call calc_stats_d_i4b(data_step, statIndex_numberScalarSolutions, stat_numberScalarSolutions(:,iGRU), numberScalarSolutions, resetStats, finalizeStats, statCounter, err)
end subroutine

attributes(global) subroutine calc_stats_kernel_flux(nGRU,data_step,data, &
statIndex_scalarCanairNetNrgFlux,stat_scalarCanairNetNrgFlux,ixscalarCanairNetNrgFlux, &
statIndex_scalarCanopyNetNrgFlux,stat_scalarCanopyNetNrgFlux,ixscalarCanopyNetNrgFlux, &
statIndex_scalarGroundNetNrgFlux,stat_scalarGroundNetNrgFlux,ixscalarGroundNetNrgFlux, &
statIndex_scalarCanopyNetLiqFlux,stat_scalarCanopyNetLiqFlux,ixscalarCanopyNetLiqFlux, &
statIndex_scalarRainfall,stat_scalarRainfall,ixscalarRainfall, &
statIndex_scalarSnowfall,stat_scalarSnowfall,ixscalarSnowfall, &
statIndex_scalarCanopySunlitPAR,stat_scalarCanopySunlitPAR,ixscalarCanopySunlitPAR, &
statIndex_scalarCanopyShadedPAR,stat_scalarCanopyShadedPAR,ixscalarCanopyShadedPAR, &
statIndex_scalarBelowCanopySolar,stat_scalarBelowCanopySolar,ixscalarBelowCanopySolar, &
statIndex_scalarCanopyAbsorbedSolar,stat_scalarCanopyAbsorbedSolar,ixscalarCanopyAbsorbedSolar, &
statIndex_scalarGroundAbsorbedSolar,stat_scalarGroundAbsorbedSolar,ixscalarGroundAbsorbedSolar, &
statIndex_scalarLWRadCanopy,stat_scalarLWRadCanopy,ixscalarLWRadCanopy, &
statIndex_scalarLWRadGround,stat_scalarLWRadGround,ixscalarLWRadGround, &
statIndex_scalarLWRadUbound2Canopy,stat_scalarLWRadUbound2Canopy,ixscalarLWRadUbound2Canopy, &
statIndex_scalarLWRadUbound2Ground,stat_scalarLWRadUbound2Ground,ixscalarLWRadUbound2Ground, &
statIndex_scalarLWRadUbound2Ubound,stat_scalarLWRadUbound2Ubound,ixscalarLWRadUbound2Ubound, &
statIndex_scalarLWRadCanopy2Ubound,stat_scalarLWRadCanopy2Ubound,ixscalarLWRadCanopy2Ubound, &
statIndex_scalarLWRadCanopy2Ground,stat_scalarLWRadCanopy2Ground,ixscalarLWRadCanopy2Ground, &
statIndex_scalarLWRadCanopy2Canopy,stat_scalarLWRadCanopy2Canopy,ixscalarLWRadCanopy2Canopy, &
statIndex_scalarLWRadGround2Ubound,stat_scalarLWRadGround2Ubound,ixscalarLWRadGround2Ubound, &
statIndex_scalarLWRadGround2Canopy,stat_scalarLWRadGround2Canopy,ixscalarLWRadGround2Canopy, &
statIndex_scalarLWNetCanopy,stat_scalarLWNetCanopy,ixscalarLWNetCanopy, &
statIndex_scalarLWNetGround,stat_scalarLWNetGround,ixscalarLWNetGround, &
statIndex_scalarLWNetUbound,stat_scalarLWNetUbound,ixscalarLWNetUbound, &
statIndex_scalarEddyDiffusCanopyTop,stat_scalarEddyDiffusCanopyTop,ixscalarEddyDiffusCanopyTop, &
statIndex_scalarFrictionVelocity,stat_scalarFrictionVelocity,ixscalarFrictionVelocity, &
statIndex_scalarWindspdCanopyTop,stat_scalarWindspdCanopyTop,ixscalarWindspdCanopyTop, &
statIndex_scalarWindspdCanopyBottom,stat_scalarWindspdCanopyBottom,ixscalarWindspdCanopyBottom, &
statIndex_scalarGroundResistance,stat_scalarGroundResistance,ixscalarGroundResistance, &
statIndex_scalarCanopyResistance,stat_scalarCanopyResistance,ixscalarCanopyResistance, &
statIndex_scalarLeafResistance,stat_scalarLeafResistance,ixscalarLeafResistance, &
statIndex_scalarSoilResistance,stat_scalarSoilResistance,ixscalarSoilResistance, &
statIndex_scalarSenHeatTotal,stat_scalarSenHeatTotal,ixscalarSenHeatTotal, &
statIndex_scalarSenHeatCanopy,stat_scalarSenHeatCanopy,ixscalarSenHeatCanopy, &
statIndex_scalarSenHeatGround,stat_scalarSenHeatGround,ixscalarSenHeatGround, &
statIndex_scalarLatHeatTotal,stat_scalarLatHeatTotal,ixscalarLatHeatTotal, &
statIndex_scalarLatHeatCanopyEvap,stat_scalarLatHeatCanopyEvap,ixscalarLatHeatCanopyEvap, &
statIndex_scalarLatHeatCanopyTrans,stat_scalarLatHeatCanopyTrans,ixscalarLatHeatCanopyTrans, &
statIndex_scalarLatHeatGround,stat_scalarLatHeatGround,ixscalarLatHeatGround, &
statIndex_scalarCanopyAdvectiveHeatFlux,stat_scalarCanopyAdvectiveHeatFlux,ixscalarCanopyAdvectiveHeatFlux, &
statIndex_scalarGroundAdvectiveHeatFlux,stat_scalarGroundAdvectiveHeatFlux,ixscalarGroundAdvectiveHeatFlux, &
statIndex_scalarCanopySublimation,stat_scalarCanopySublimation,ixscalarCanopySublimation, &
statIndex_scalarSnowSublimation,stat_scalarSnowSublimation,ixscalarSnowSublimation, &
statIndex_scalarStomResistSunlit,stat_scalarStomResistSunlit,ixscalarStomResistSunlit, &
statIndex_scalarStomResistShaded,stat_scalarStomResistShaded,ixscalarStomResistShaded, &
statIndex_scalarPhotosynthesisSunlit,stat_scalarPhotosynthesisSunlit,ixscalarPhotosynthesisSunlit, &
statIndex_scalarPhotosynthesisShaded,stat_scalarPhotosynthesisShaded,ixscalarPhotosynthesisShaded, &
statIndex_scalarCanopyTranspiration,stat_scalarCanopyTranspiration,ixscalarCanopyTranspiration, &
statIndex_scalarCanopyEvaporation,stat_scalarCanopyEvaporation,ixscalarCanopyEvaporation, &
statIndex_scalarGroundEvaporation,stat_scalarGroundEvaporation,ixscalarGroundEvaporation, &
statIndex_scalarThroughfallSnow,stat_scalarThroughfallSnow,ixscalarThroughfallSnow, &
statIndex_scalarThroughfallRain,stat_scalarThroughfallRain,ixscalarThroughfallRain, &
statIndex_scalarCanopySnowUnloading,stat_scalarCanopySnowUnloading,ixscalarCanopySnowUnloading, &
statIndex_scalarCanopyLiqDrainage,stat_scalarCanopyLiqDrainage,ixscalarCanopyLiqDrainage, &
statIndex_scalarCanopyMeltFreeze,stat_scalarCanopyMeltFreeze,ixscalarCanopyMeltFreeze, &
statIndex_scalarSnowDrainage,stat_scalarSnowDrainage,ixscalarSnowDrainage, &
statIndex_scalarRainPlusMelt,stat_scalarRainPlusMelt,ixscalarRainPlusMelt, &
statIndex_scalarMaxInfilRate,stat_scalarMaxInfilRate,ixscalarMaxInfilRate, &
statIndex_scalarInfiltration,stat_scalarInfiltration,ixscalarInfiltration, &
statIndex_scalarExfiltration,stat_scalarExfiltration,ixscalarExfiltration, &
statIndex_scalarSurfaceRunoff,stat_scalarSurfaceRunoff,ixscalarSurfaceRunoff, &
statIndex_scalarSurfaceRunoff_IE,stat_scalarSurfaceRunoff_IE,ixscalarSurfaceRunoff_IE, &
statIndex_scalarSurfaceRunoff_SE,stat_scalarSurfaceRunoff_SE,ixscalarSurfaceRunoff_SE, &
statIndex_scalarSoilBaseflow,stat_scalarSoilBaseflow,ixscalarSoilBaseflow, &
statIndex_scalarSoilDrainage,stat_scalarSoilDrainage,ixscalarSoilDrainage, &
statIndex_scalarAquiferRecharge,stat_scalarAquiferRecharge,ixscalarAquiferRecharge, &
statIndex_scalarAquiferTranspire,stat_scalarAquiferTranspire,ixscalarAquiferTranspire, &
statIndex_scalarAquiferBaseflow,stat_scalarAquiferBaseflow,ixscalarAquiferBaseflow, &
statIndex_scalarTotalET,stat_scalarTotalET,ixscalarTotalET, &
statIndex_scalarTotalRunoff,stat_scalarTotalRunoff,ixscalarTotalRunoff, &
statIndex_scalarNetRadiation,stat_scalarNetRadiation,ixscalarNetRadiation, &
resetStats,finalizeStats,statCounter)
implicit none
integer(i4b),intent(in),value :: nGRU
real(rkind),intent(in),value :: data_step
real(rkind) :: data(:,:)
integer(i4b),intent(in),value :: ixScalarCanairNetNrgFlux
real(rkind) :: stat_scalarCanairNetNrgFlux(:,:)
integer(i4b) :: statIndex_scalarCanairNetNrgFlux(:)
integer(i4b),intent(in),value :: ixScalarCanopyNetNrgFlux
real(rkind) :: stat_scalarCanopyNetNrgFlux(:,:)
integer(i4b) :: statIndex_scalarCanopyNetNrgFlux(:)
integer(i4b),intent(in),value :: ixScalarGroundNetNrgFlux
real(rkind) :: stat_scalarGroundNetNrgFlux(:,:)
integer(i4b) :: statIndex_scalarGroundNetNrgFlux(:)
integer(i4b),intent(in),value :: ixScalarCanopyNetLiqFlux
real(rkind) :: stat_scalarCanopyNetLiqFlux(:,:)
integer(i4b) :: statIndex_scalarCanopyNetLiqFlux(:)
integer(i4b),intent(in),value :: ixScalarRainfall
real(rkind) :: stat_scalarRainfall(:,:)
integer(i4b) :: statIndex_scalarRainfall(:)
integer(i4b),intent(in),value :: ixScalarSnowfall
real(rkind) :: stat_scalarSnowfall(:,:)
integer(i4b) :: statIndex_scalarSnowfall(:)
integer(i4b),intent(in),value :: ixScalarCanopySunlitPAR
real(rkind) :: stat_scalarCanopySunlitPAR(:,:)
integer(i4b) :: statIndex_scalarCanopySunlitPAR(:)
integer(i4b),intent(in),value :: ixScalarCanopyShadedPAR
real(rkind) :: stat_scalarCanopyShadedPAR(:,:)
integer(i4b) :: statIndex_scalarCanopyShadedPAR(:)
integer(i4b),intent(in),value :: ixScalarBelowCanopySolar
real(rkind) :: stat_scalarBelowCanopySolar(:,:)
integer(i4b) :: statIndex_scalarBelowCanopySolar(:)
integer(i4b),intent(in),value :: ixScalarCanopyAbsorbedSolar
real(rkind) :: stat_scalarCanopyAbsorbedSolar(:,:)
integer(i4b) :: statIndex_scalarCanopyAbsorbedSolar(:)
integer(i4b),intent(in),value :: ixScalarGroundAbsorbedSolar
real(rkind) :: stat_scalarGroundAbsorbedSolar(:,:)
integer(i4b) :: statIndex_scalarGroundAbsorbedSolar(:)
integer(i4b),intent(in),value :: ixScalarLWRadCanopy
real(rkind) :: stat_scalarLWRadCanopy(:,:)
integer(i4b) :: statIndex_scalarLWRadCanopy(:)
integer(i4b),intent(in),value :: ixScalarLWRadGround
real(rkind) :: stat_scalarLWRadGround(:,:)
integer(i4b) :: statIndex_scalarLWRadGround(:)
integer(i4b),intent(in),value :: ixScalarLWRadUbound2Canopy
real(rkind) :: stat_scalarLWRadUbound2Canopy(:,:)
integer(i4b) :: statIndex_scalarLWRadUbound2Canopy(:)
integer(i4b),intent(in),value :: ixScalarLWRadUbound2Ground
real(rkind) :: stat_scalarLWRadUbound2Ground(:,:)
integer(i4b) :: statIndex_scalarLWRadUbound2Ground(:)
integer(i4b),intent(in),value :: ixScalarLWRadUbound2Ubound
real(rkind) :: stat_scalarLWRadUbound2Ubound(:,:)
integer(i4b) :: statIndex_scalarLWRadUbound2Ubound(:)
integer(i4b),intent(in),value :: ixScalarLWRadCanopy2Ubound
real(rkind) :: stat_scalarLWRadCanopy2Ubound(:,:)
integer(i4b) :: statIndex_scalarLWRadCanopy2Ubound(:)
integer(i4b),intent(in),value :: ixScalarLWRadCanopy2Ground
real(rkind) :: stat_scalarLWRadCanopy2Ground(:,:)
integer(i4b) :: statIndex_scalarLWRadCanopy2Ground(:)
integer(i4b),intent(in),value :: ixScalarLWRadCanopy2Canopy
real(rkind) :: stat_scalarLWRadCanopy2Canopy(:,:)
integer(i4b) :: statIndex_scalarLWRadCanopy2Canopy(:)
integer(i4b),intent(in),value :: ixScalarLWRadGround2Ubound
real(rkind) :: stat_scalarLWRadGround2Ubound(:,:)
integer(i4b) :: statIndex_scalarLWRadGround2Ubound(:)
integer(i4b),intent(in),value :: ixScalarLWRadGround2Canopy
real(rkind) :: stat_scalarLWRadGround2Canopy(:,:)
integer(i4b) :: statIndex_scalarLWRadGround2Canopy(:)
integer(i4b),intent(in),value :: ixScalarLWNetCanopy
real(rkind) :: stat_scalarLWNetCanopy(:,:)
integer(i4b) :: statIndex_scalarLWNetCanopy(:)
integer(i4b),intent(in),value :: ixScalarLWNetGround
real(rkind) :: stat_scalarLWNetGround(:,:)
integer(i4b) :: statIndex_scalarLWNetGround(:)
integer(i4b),intent(in),value :: ixScalarLWNetUbound
real(rkind) :: stat_scalarLWNetUbound(:,:)
integer(i4b) :: statIndex_scalarLWNetUbound(:)
integer(i4b),intent(in),value :: ixScalarEddyDiffusCanopyTop
real(rkind) :: stat_scalarEddyDiffusCanopyTop(:,:)
integer(i4b) :: statIndex_scalarEddyDiffusCanopyTop(:)
integer(i4b),intent(in),value :: ixScalarFrictionVelocity
real(rkind) :: stat_scalarFrictionVelocity(:,:)
integer(i4b) :: statIndex_scalarFrictionVelocity(:)
integer(i4b),intent(in),value :: ixScalarWindspdCanopyTop
real(rkind) :: stat_scalarWindspdCanopyTop(:,:)
integer(i4b) :: statIndex_scalarWindspdCanopyTop(:)
integer(i4b),intent(in),value :: ixScalarWindspdCanopyBottom
real(rkind) :: stat_scalarWindspdCanopyBottom(:,:)
integer(i4b) :: statIndex_scalarWindspdCanopyBottom(:)
integer(i4b),intent(in),value :: ixScalarGroundResistance
real(rkind) :: stat_scalarGroundResistance(:,:)
integer(i4b) :: statIndex_scalarGroundResistance(:)
integer(i4b),intent(in),value :: ixScalarCanopyResistance
real(rkind) :: stat_scalarCanopyResistance(:,:)
integer(i4b) :: statIndex_scalarCanopyResistance(:)
integer(i4b),intent(in),value :: ixScalarLeafResistance
real(rkind) :: stat_scalarLeafResistance(:,:)
integer(i4b) :: statIndex_scalarLeafResistance(:)
integer(i4b),intent(in),value :: ixScalarSoilResistance
real(rkind) :: stat_scalarSoilResistance(:,:)
integer(i4b) :: statIndex_scalarSoilResistance(:)
integer(i4b),intent(in),value :: ixScalarSenHeatTotal
real(rkind) :: stat_scalarSenHeatTotal(:,:)
integer(i4b) :: statIndex_scalarSenHeatTotal(:)
integer(i4b),intent(in),value :: ixScalarSenHeatCanopy
real(rkind) :: stat_scalarSenHeatCanopy(:,:)
integer(i4b) :: statIndex_scalarSenHeatCanopy(:)
integer(i4b),intent(in),value :: ixScalarSenHeatGround
real(rkind) :: stat_scalarSenHeatGround(:,:)
integer(i4b) :: statIndex_scalarSenHeatGround(:)
integer(i4b),intent(in),value :: ixScalarLatHeatTotal
real(rkind) :: stat_scalarLatHeatTotal(:,:)
integer(i4b) :: statIndex_scalarLatHeatTotal(:)
integer(i4b),intent(in),value :: ixScalarLatHeatCanopyEvap
real(rkind) :: stat_scalarLatHeatCanopyEvap(:,:)
integer(i4b) :: statIndex_scalarLatHeatCanopyEvap(:)
integer(i4b),intent(in),value :: ixScalarLatHeatCanopyTrans
real(rkind) :: stat_scalarLatHeatCanopyTrans(:,:)
integer(i4b) :: statIndex_scalarLatHeatCanopyTrans(:)
integer(i4b),intent(in),value :: ixScalarLatHeatGround
real(rkind) :: stat_scalarLatHeatGround(:,:)
integer(i4b) :: statIndex_scalarLatHeatGround(:)
integer(i4b),intent(in),value :: ixScalarCanopyAdvectiveHeatFlux
real(rkind) :: stat_scalarCanopyAdvectiveHeatFlux(:,:)
integer(i4b) :: statIndex_scalarCanopyAdvectiveHeatFlux(:)
integer(i4b),intent(in),value :: ixScalarGroundAdvectiveHeatFlux
real(rkind) :: stat_scalarGroundAdvectiveHeatFlux(:,:)
integer(i4b) :: statIndex_scalarGroundAdvectiveHeatFlux(:)
integer(i4b),intent(in),value :: ixScalarCanopySublimation
real(rkind) :: stat_scalarCanopySublimation(:,:)
integer(i4b) :: statIndex_scalarCanopySublimation(:)
integer(i4b),intent(in),value :: ixScalarSnowSublimation
real(rkind) :: stat_scalarSnowSublimation(:,:)
integer(i4b) :: statIndex_scalarSnowSublimation(:)
integer(i4b),intent(in),value :: ixScalarStomResistSunlit
real(rkind) :: stat_scalarStomResistSunlit(:,:)
integer(i4b) :: statIndex_scalarStomResistSunlit(:)
integer(i4b),intent(in),value :: ixScalarStomResistShaded
real(rkind) :: stat_scalarStomResistShaded(:,:)
integer(i4b) :: statIndex_scalarStomResistShaded(:)
integer(i4b),intent(in),value :: ixScalarPhotosynthesisSunlit
real(rkind) :: stat_scalarPhotosynthesisSunlit(:,:)
integer(i4b) :: statIndex_scalarPhotosynthesisSunlit(:)
integer(i4b),intent(in),value :: ixScalarPhotosynthesisShaded
real(rkind) :: stat_scalarPhotosynthesisShaded(:,:)
integer(i4b) :: statIndex_scalarPhotosynthesisShaded(:)
integer(i4b),intent(in),value :: ixScalarCanopyTranspiration
real(rkind) :: stat_scalarCanopyTranspiration(:,:)
integer(i4b) :: statIndex_scalarCanopyTranspiration(:)
integer(i4b),intent(in),value :: ixScalarCanopyEvaporation
real(rkind) :: stat_scalarCanopyEvaporation(:,:)
integer(i4b) :: statIndex_scalarCanopyEvaporation(:)
integer(i4b),intent(in),value :: ixScalarGroundEvaporation
real(rkind) :: stat_scalarGroundEvaporation(:,:)
integer(i4b) :: statIndex_scalarGroundEvaporation(:)
integer(i4b),intent(in),value :: ixScalarThroughfallSnow
real(rkind) :: stat_scalarThroughfallSnow(:,:)
integer(i4b) :: statIndex_scalarThroughfallSnow(:)
integer(i4b),intent(in),value :: ixScalarThroughfallRain
real(rkind) :: stat_scalarThroughfallRain(:,:)
integer(i4b) :: statIndex_scalarThroughfallRain(:)
integer(i4b),intent(in),value :: ixScalarCanopySnowUnloading
real(rkind) :: stat_scalarCanopySnowUnloading(:,:)
integer(i4b) :: statIndex_scalarCanopySnowUnloading(:)
integer(i4b),intent(in),value :: ixScalarCanopyLiqDrainage
real(rkind) :: stat_scalarCanopyLiqDrainage(:,:)
integer(i4b) :: statIndex_scalarCanopyLiqDrainage(:)
integer(i4b),intent(in),value :: ixScalarCanopyMeltFreeze
real(rkind) :: stat_scalarCanopyMeltFreeze(:,:)
integer(i4b) :: statIndex_scalarCanopyMeltFreeze(:)
integer(i4b),intent(in),value :: ixScalarSnowDrainage
real(rkind) :: stat_scalarSnowDrainage(:,:)
integer(i4b) :: statIndex_scalarSnowDrainage(:)
integer(i4b),intent(in),value :: ixScalarRainPlusMelt
real(rkind) :: stat_scalarRainPlusMelt(:,:)
integer(i4b) :: statIndex_scalarRainPlusMelt(:)
integer(i4b),intent(in),value :: ixScalarMaxInfilRate
real(rkind) :: stat_scalarMaxInfilRate(:,:)
integer(i4b) :: statIndex_scalarMaxInfilRate(:)
integer(i4b),intent(in),value :: ixScalarInfiltration
real(rkind) :: stat_scalarInfiltration(:,:)
integer(i4b) :: statIndex_scalarInfiltration(:)
integer(i4b),intent(in),value :: ixScalarExfiltration
real(rkind) :: stat_scalarExfiltration(:,:)
integer(i4b) :: statIndex_scalarExfiltration(:)
integer(i4b),intent(in),value :: ixScalarSurfaceRunoff
real(rkind) :: stat_scalarSurfaceRunoff(:,:)
integer(i4b) :: statIndex_scalarSurfaceRunoff(:)
integer(i4b),intent(in),value :: ixScalarSurfaceRunoff_IE
real(rkind) :: stat_scalarSurfaceRunoff_IE(:,:)
integer(i4b) :: statIndex_scalarSurfaceRunoff_IE(:)
integer(i4b),intent(in),value :: ixScalarSurfaceRunoff_SE
real(rkind) :: stat_scalarSurfaceRunoff_SE(:,:)
integer(i4b) :: statIndex_scalarSurfaceRunoff_SE(:)
integer(i4b),intent(in),value :: ixScalarSoilBaseflow
real(rkind) :: stat_scalarSoilBaseflow(:,:)
integer(i4b) :: statIndex_scalarSoilBaseflow(:)
integer(i4b),intent(in),value :: ixScalarSoilDrainage
real(rkind) :: stat_scalarSoilDrainage(:,:)
integer(i4b) :: statIndex_scalarSoilDrainage(:)
integer(i4b),intent(in),value :: ixScalarAquiferRecharge
real(rkind) :: stat_scalarAquiferRecharge(:,:)
integer(i4b) :: statIndex_scalarAquiferRecharge(:)
integer(i4b),intent(in),value :: ixScalarAquiferTranspire
real(rkind) :: stat_scalarAquiferTranspire(:,:)
integer(i4b) :: statIndex_scalarAquiferTranspire(:)
integer(i4b),intent(in),value :: ixScalarAquiferBaseflow
real(rkind) :: stat_scalarAquiferBaseflow(:,:)
integer(i4b) :: statIndex_scalarAquiferBaseflow(:)
integer(i4b),intent(in),value :: ixScalarTotalET
real(rkind) :: stat_scalarTotalET(:,:)
integer(i4b) :: statIndex_scalarTotalET(:)
integer(i4b),intent(in),value :: ixScalarTotalRunoff
real(rkind) :: stat_scalarTotalRunoff(:,:)
integer(i4b) :: statIndex_scalarTotalRunoff(:)
integer(i4b),intent(in),value :: ixScalarNetRadiation
real(rkind) :: stat_scalarNetRadiation(:,:)
integer(i4b) :: statIndex_scalarNetRadiation(:)
     logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency

integer(i4b) :: iGRU, err
 iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
 if (iGRU .gt. nGRU) return
return
call calc_stats_d_real(data_step,statIndex_scalarCanairNetNrgFlux, stat_scalarCanairNetNrgFlux(:,iGRU), data(ixScalarCanairNetNrgFlux,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopyNetNrgFlux, stat_scalarCanopyNetNrgFlux(:,iGRU), data(ixScalarCanopyNetNrgFlux,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarGroundNetNrgFlux, stat_scalarGroundNetNrgFlux(:,iGRU), data(ixScalarGroundNetNrgFlux,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopyNetLiqFlux, stat_scalarCanopyNetLiqFlux(:,iGRU), data(ixScalarCanopyNetLiqFlux,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarRainfall, stat_scalarRainfall(:,iGRU), data(ixScalarRainfall,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarSnowfall, stat_scalarSnowfall(:,iGRU), data(ixScalarSnowfall,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopySunlitPAR, stat_scalarCanopySunlitPAR(:,iGRU), data(ixScalarCanopySunlitPAR,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopyShadedPAR, stat_scalarCanopyShadedPAR(:,iGRU), data(ixScalarCanopyShadedPAR,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarBelowCanopySolar, stat_scalarBelowCanopySolar(:,iGRU), data(ixScalarBelowCanopySolar,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopyAbsorbedSolar, stat_scalarCanopyAbsorbedSolar(:,iGRU), data(ixScalarCanopyAbsorbedSolar,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarGroundAbsorbedSolar, stat_scalarGroundAbsorbedSolar(:,iGRU), data(ixScalarGroundAbsorbedSolar,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWRadCanopy, stat_scalarLWRadCanopy(:,iGRU), data(ixScalarLWRadCanopy,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWRadGround, stat_scalarLWRadGround(:,iGRU), data(ixScalarLWRadGround,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWRadUbound2Canopy, stat_scalarLWRadUbound2Canopy(:,iGRU), data(ixScalarLWRadUbound2Canopy,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWRadUbound2Ground, stat_scalarLWRadUbound2Ground(:,iGRU), data(ixScalarLWRadUbound2Ground,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWRadUbound2Ubound, stat_scalarLWRadUbound2Ubound(:,iGRU), data(ixScalarLWRadUbound2Ubound,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWRadCanopy2Ubound, stat_scalarLWRadCanopy2Ubound(:,iGRU), data(ixScalarLWRadCanopy2Ubound,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWRadCanopy2Ground, stat_scalarLWRadCanopy2Ground(:,iGRU), data(ixScalarLWRadCanopy2Ground,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWRadCanopy2Canopy, stat_scalarLWRadCanopy2Canopy(:,iGRU), data(ixScalarLWRadCanopy2Canopy,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWRadGround2Ubound, stat_scalarLWRadGround2Ubound(:,iGRU), data(ixScalarLWRadGround2Ubound,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWRadGround2Canopy, stat_scalarLWRadGround2Canopy(:,iGRU), data(ixScalarLWRadGround2Canopy,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWNetCanopy, stat_scalarLWNetCanopy(:,iGRU), data(ixScalarLWNetCanopy,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWNetGround, stat_scalarLWNetGround(:,iGRU), data(ixScalarLWNetGround,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLWNetUbound, stat_scalarLWNetUbound(:,iGRU), data(ixScalarLWNetUbound,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarEddyDiffusCanopyTop, stat_scalarEddyDiffusCanopyTop(:,iGRU), data(ixScalarEddyDiffusCanopyTop,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarFrictionVelocity, stat_scalarFrictionVelocity(:,iGRU), data(ixScalarFrictionVelocity,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarWindspdCanopyTop, stat_scalarWindspdCanopyTop(:,iGRU), data(ixScalarWindspdCanopyTop,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarWindspdCanopyBottom, stat_scalarWindspdCanopyBottom(:,iGRU), data(ixScalarWindspdCanopyBottom,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarGroundResistance, stat_scalarGroundResistance(:,iGRU), data(ixScalarGroundResistance,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopyResistance, stat_scalarCanopyResistance(:,iGRU), data(ixScalarCanopyResistance,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLeafResistance, stat_scalarLeafResistance(:,iGRU), data(ixScalarLeafResistance,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarSoilResistance, stat_scalarSoilResistance(:,iGRU), data(ixScalarSoilResistance,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarSenHeatTotal, stat_scalarSenHeatTotal(:,iGRU), data(ixScalarSenHeatTotal,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarSenHeatCanopy, stat_scalarSenHeatCanopy(:,iGRU), data(ixScalarSenHeatCanopy,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarSenHeatGround, stat_scalarSenHeatGround(:,iGRU), data(ixScalarSenHeatGround,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLatHeatTotal, stat_scalarLatHeatTotal(:,iGRU), data(ixScalarLatHeatTotal,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLatHeatCanopyEvap, stat_scalarLatHeatCanopyEvap(:,iGRU), data(ixScalarLatHeatCanopyEvap,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLatHeatCanopyTrans, stat_scalarLatHeatCanopyTrans(:,iGRU), data(ixScalarLatHeatCanopyTrans,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarLatHeatGround, stat_scalarLatHeatGround(:,iGRU), data(ixScalarLatHeatGround,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopyAdvectiveHeatFlux, stat_scalarCanopyAdvectiveHeatFlux(:,iGRU), data(ixScalarCanopyAdvectiveHeatFlux,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarGroundAdvectiveHeatFlux, stat_scalarGroundAdvectiveHeatFlux(:,iGRU), data(ixScalarGroundAdvectiveHeatFlux,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopySublimation, stat_scalarCanopySublimation(:,iGRU), data(ixScalarCanopySublimation,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarSnowSublimation, stat_scalarSnowSublimation(:,iGRU), data(ixScalarSnowSublimation,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarStomResistSunlit, stat_scalarStomResistSunlit(:,iGRU), data(ixScalarStomResistSunlit,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarStomResistShaded, stat_scalarStomResistShaded(:,iGRU), data(ixScalarStomResistShaded,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarPhotosynthesisSunlit, stat_scalarPhotosynthesisSunlit(:,iGRU), data(ixScalarPhotosynthesisSunlit,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarPhotosynthesisShaded, stat_scalarPhotosynthesisShaded(:,iGRU), data(ixScalarPhotosynthesisShaded,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopyTranspiration, stat_scalarCanopyTranspiration(:,iGRU), data(ixScalarCanopyTranspiration,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopyEvaporation, stat_scalarCanopyEvaporation(:,iGRU), data(ixScalarCanopyEvaporation,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarGroundEvaporation, stat_scalarGroundEvaporation(:,iGRU), data(ixScalarGroundEvaporation,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarThroughfallSnow, stat_scalarThroughfallSnow(:,iGRU), data(ixScalarThroughfallSnow,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarThroughfallRain, stat_scalarThroughfallRain(:,iGRU), data(ixScalarThroughfallRain,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopySnowUnloading, stat_scalarCanopySnowUnloading(:,iGRU), data(ixScalarCanopySnowUnloading,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopyLiqDrainage, stat_scalarCanopyLiqDrainage(:,iGRU), data(ixScalarCanopyLiqDrainage,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarCanopyMeltFreeze, stat_scalarCanopyMeltFreeze(:,iGRU), data(ixScalarCanopyMeltFreeze,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarSnowDrainage, stat_scalarSnowDrainage(:,iGRU), data(ixScalarSnowDrainage,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarRainPlusMelt, stat_scalarRainPlusMelt(:,iGRU), data(ixScalarRainPlusMelt,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarMaxInfilRate, stat_scalarMaxInfilRate(:,iGRU), data(ixScalarMaxInfilRate,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarInfiltration, stat_scalarInfiltration(:,iGRU), data(ixScalarInfiltration,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarExfiltration, stat_scalarExfiltration(:,iGRU), data(ixScalarExfiltration,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarSurfaceRunoff, stat_scalarSurfaceRunoff(:,iGRU), data(ixScalarSurfaceRunoff,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarSurfaceRunoff_IE, stat_scalarSurfaceRunoff_IE(:,iGRU), data(ixScalarSurfaceRunoff_IE,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarSurfaceRunoff_SE, stat_scalarSurfaceRunoff_SE(:,iGRU), data(ixScalarSurfaceRunoff_SE,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarSoilBaseflow, stat_scalarSoilBaseflow(:,iGRU), data(ixScalarSoilBaseflow,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarSoilDrainage, stat_scalarSoilDrainage(:,iGRU), data(ixScalarSoilDrainage,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarAquiferRecharge, stat_scalarAquiferRecharge(:,iGRU), data(ixScalarAquiferRecharge,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarAquiferTranspire, stat_scalarAquiferTranspire(:,iGRU), data(ixScalarAquiferTranspire,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarAquiferBaseflow, stat_scalarAquiferBaseflow(:,iGRU), data(ixScalarAquiferBaseflow,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarTotalET, stat_scalarTotalET(:,iGRU), data(ixScalarTotalET,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarTotalRunoff, stat_scalarTotalRunoff(:,iGRU), data(ixScalarTotalRunoff,iGRU), resetStats,finalizeStats,statCounter,err)
call calc_stats_d_real(data_step,statIndex_scalarNetRadiation, stat_scalarNetRadiation(:,iGRU), data(ixScalarNetRadiation,iGRU), resetStats,finalizeStats,statCounter,err)

end subroutine

 attributes(device) subroutine calc_stats_d_real(data_step,statIndex,stat,tdata,resetStats,finalizeStats,statCounter,err)
 USE nrtype
 ! data structures
 USE data_types,only:var_info,ilength,dlength ! type dec for meta data structures
 USE var_lookup,only:maxVarFreq       ! # of output frequencies
 ! global variables
!  USE globalData,only:data_step        ! forcing timestep
 ! structures of named variables
 USE var_lookup,only:iLookVarType     ! named variables for variable types
 USE var_lookup,only:iLookFREQ        ! named variables for output frequency
 USE var_lookup,only:iLookSTAT        ! named variables for output statistics
 USE var_lookup,only:iLookTIME        ! named variables for time information
 implicit none
 real(rkind) :: data_step
 ! input variables
!  class(var_info),intent(in)         :: meta              ! meta data structure
 integer(i4b), intent(in) :: statIndex(:)
 real(rkind)       ,intent(inout)      :: stat(:)              ! statistics structure
 real(rkind)       ,intent(in)         :: tdata             ! data value
 logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency
 ! output variables
 integer(i4b)   ,intent(out)        :: err               ! error code
!  character(*)   ,intent(out)        :: message           ! error message
 ! internals
!  real(rkind),dimension(maxvarFreq*2)   :: tstat             ! temporary stats vector
 integer(i4b)                       :: iFreq             ! index of output frequency
 ! initialize error control
 err=0


 ! ---------------------------------------------
 ! reset statistics at new frequency period
 ! ---------------------------------------------
 if (size(statIndex).eq.1) return
!  print*, size(resetStats), maxVarFreq
 do iFreq=1,maxVarFreq                              ! loop through output statistics
  if(resetStats(iFreq))then                         ! flag to reset statistics
   if(statIndex(iFreq)==integerMissing) cycle  ! don't bother if output frequency is not desired for a given variable
!    if(meta%varType/=iLookVarType%outstat) cycle     ! only calculate stats for scalars
   select case(statIndex(iFreq))               ! act depending on the statistic
    ! -------------------------------------------------------------------------------------
    case (iLookSTAT%totl)                           ! * summation over period                  
     stat(iFreq) = 0._rkind                           !     - resets stat at beginning of period
    case (iLookSTAT%mean)                           ! * mean over period                       
     stat(iFreq) = 0._rkind                           !     - resets stat at beginning of period
    case (iLookSTAT%vari)                           ! * variance over period                   
     stat(iFreq) = 0._rkind                           !     - resets E[X^2] term in var calc    
     stat(maxVarFreq+iFreq) = 0._rkind                !     - resets E[X]^2 term                 
    case (iLookSTAT%mini)                           ! * minimum over period                    
     stat(iFreq) = huge(stat(iFreq))              !     - resets stat at beginning of period 
    case (iLookSTAT%maxi)                           ! * maximum over period                    
     stat(iFreq) = -huge(stat(iFreq))             !     - resets stat at beginning of period 
    case (iLookSTAT%mode)                           ! * mode over period      
     stat(iFreq) = realMissing                     !     - does not work
    case (iLookSTAT%inst)                           ! * instantaneous -- no need to reset
    case default
    !  message=trim(message)//'unable to identify type of statistic [reset]'
     err=20; return
    ! -------------------------------------------------------------------------------------
   end select
  end if
 end do ! looping through output frequencies

 ! ---------------------------------------------
 ! Calculate each statistic that is requested by user
 ! ---------------------------------------------
 do iFreq=1,maxVarFreq                                ! loop through output statistics
  if(statIndex(iFreq)==integerMissing) cycle     ! don't bother if output frequency is not desired for a given variab;e
!   if(meta%varType/=iLookVarType%outstat) cycle        ! only calculate stats for scalars
  select case(statIndex(iFreq))                  ! act depending on the statistic
   ! -------------------------------------------------------------------------------------
   case (iLookSTAT%inst)                              ! * instantaneous value 
    stat(iFreq) = tdata                              !     - data at a given time
   case (iLookSTAT%totl)                              ! * summation over period                    
    stat(iFreq) = stat(iFreq) + tdata*data_step     !     - increment data 
   case (iLookSTAT%mean)                              ! * mean over period                       
    stat(iFreq) = stat(iFreq) + tdata               !     -  increment data
   case (iLookSTAT%vari)                              ! * variance over period                   
    stat(iFreq) = stat(iFreq) + tdata**2                     ! - E[X^2] term in var calc    
    stat(maxVarFreq+iFreq) = stat(maxVarFreq+iFreq) + tdata  ! - E[X]^2 term                 
   case (iLookSTAT%mini)                              ! * minimum over period                    
    if (tdata<stat(iFreq)) stat(iFreq) = tdata      !     - check value 
   case (iLookSTAT%maxi)                              ! * maximum over period                    
    if (tdata>stat(iFreq)) stat(iFreq) = tdata      !     - check value 
   case (iLookSTAT%mode)                              ! * mode over period (does not workind)       
    stat(iFreq) = realMissing
   case default
    ! message=trim(message)//'unable to identify type of statistic [calculating stats]'
    err=20; return
   ! -------------------------------------------------------------------------------------
  end select
 end do ! looping through output frequencies

 ! ---------------------------------------------
 ! finalize statistics at end of frequency period
 ! ---------------------------------------------
 do iFreq=1,maxVarFreq                                ! loop through output statistics
  if(finalizeStats(iFreq))then
   if(statIndex(iFreq)==integerMissing) cycle     ! don't bother if output frequency is not desired for a given variable
!    if(meta%varType/=iLookVarType%outstat) cycle        ! only calculate stats for scalars
   select case(statIndex(iFreq))                  ! act depending on the statistic
    ! -------------------------------------------------------------------------------------
    case (iLookSTAT%mean)                              ! * mean over period
     stat(iFreq) = stat(iFreq)/statCounter(iFreq)    !     - normalize sum into mean
    case (iLookSTAT%vari)                              ! * variance over period
     stat(maxVarFreq+iFreq) = stat(maxVarFreq+1)/statCounter(iFreq)            ! E[X] term
     stat(iFreq) = stat(iFreq)/statCounter(iFreq) - stat(maxVarFreq+iFreq)**2 ! full variance
    case default ! do nothing -- don't need finalization for most stats
    ! -------------------------------------------------------------------------------------
   end select
  end if
 end do ! looping through output frequencies


 end subroutine calc_stats_d_real

  subroutine calc_stats_h_real(data_step,statIndex,stat,tdata,resetStats,finalizeStats,statCounter,err)
 USE nrtype
 ! data structures
 USE data_types,only:var_info,ilength,dlength ! type dec for meta data structures
 USE var_lookup,only:maxVarFreq       ! # of output frequencies
 ! global variables
!  USE globalData,only:data_step        ! forcing timestep
 ! structures of named variables
 USE var_lookup,only:iLookVarType     ! named variables for variable types
 USE var_lookup,only:iLookFREQ        ! named variables for output frequency
 USE var_lookup,only:iLookSTAT        ! named variables for output statistics
 USE var_lookup,only:iLookTIME        ! named variables for time information
 implicit none
 real(rkind) :: data_step
 ! input variables
!  class(var_info),intent(in)         :: meta              ! meta data structure
 integer(i4b), intent(in) :: statIndex(:)
 real(rkind)       ,intent(inout)      :: stat(:)              ! statistics structure
 real(rkind)       ,intent(in)         :: tdata             ! data value
 logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency
 ! output variables
 integer(i4b)   ,intent(out)        :: err               ! error code
!  character(*)   ,intent(out)        :: message           ! error message
 ! internals
!  real(rkind),dimension(maxvarFreq*2)   :: tstat             ! temporary stats vector
 integer(i4b)                       :: iFreq             ! index of output frequency
 ! initialize error control
 err=0


 ! ---------------------------------------------
 ! reset statistics at new frequency period
 ! ---------------------------------------------
 if (size(statIndex).eq.1) return
 do iFreq=1,maxVarFreq                              ! loop through output statistics
  if(resetStats(iFreq))then                         ! flag to reset statistics
   if(statIndex(iFreq)==integerMissing) cycle  ! don't bother if output frequency is not desired for a given variable
!    if(meta%varType/=iLookVarType%outstat) cycle     ! only calculate stats for scalars
   select case(statIndex(iFreq))               ! act depending on the statistic
    ! -------------------------------------------------------------------------------------
    case (iLookSTAT%totl)                           ! * summation over period                  
     stat(iFreq) = 0._rkind                           !     - resets stat at beginning of period
    case (iLookSTAT%mean)                           ! * mean over period                       
     stat(iFreq) = 0._rkind                           !     - resets stat at beginning of period
    case (iLookSTAT%vari)                           ! * variance over period                   
     stat(iFreq) = 0._rkind                           !     - resets E[X^2] term in var calc    
     stat(maxVarFreq+iFreq) = 0._rkind                !     - resets E[X]^2 term                 
    case (iLookSTAT%mini)                           ! * minimum over period                    
     stat(iFreq) = huge(stat(iFreq))              !     - resets stat at beginning of period 
    case (iLookSTAT%maxi)                           ! * maximum over period                    
     stat(iFreq) = -huge(stat(iFreq))             !     - resets stat at beginning of period 
    case (iLookSTAT%mode)                           ! * mode over period      
     stat(iFreq) = realMissing                     !     - does not work
    case (iLookSTAT%inst)                           ! * instantaneous -- no need to reset
    case default
    !  message=trim(message)//'unable to identify type of statistic [reset]'
     err=20; return
    ! -------------------------------------------------------------------------------------
   end select
  end if
 end do ! looping through output frequencies

 ! ---------------------------------------------
 ! Calculate each statistic that is requested by user
 ! ---------------------------------------------
 do iFreq=1,maxVarFreq                                ! loop through output statistics
  if(statIndex(iFreq)==integerMissing) cycle     ! don't bother if output frequency is not desired for a given variab;e
!   if(meta%varType/=iLookVarType%outstat) cycle        ! only calculate stats for scalars
  select case(statIndex(iFreq))                  ! act depending on the statistic
   ! -------------------------------------------------------------------------------------
   case (iLookSTAT%inst)                              ! * instantaneous value 
    stat(iFreq) = tdata                              !     - data at a given time
   case (iLookSTAT%totl)                              ! * summation over period                    
    stat(iFreq) = stat(iFreq) + tdata*data_step     !     - increment data 
   case (iLookSTAT%mean)                              ! * mean over period                       
    stat(iFreq) = stat(iFreq) + tdata               !     -  increment data
   case (iLookSTAT%vari)                              ! * variance over period                   
    stat(iFreq) = stat(iFreq) + tdata**2                     ! - E[X^2] term in var calc    
    stat(maxVarFreq+iFreq) = stat(maxVarFreq+iFreq) + tdata  ! - E[X]^2 term                 
   case (iLookSTAT%mini)                              ! * minimum over period                    
    if (tdata<stat(iFreq)) stat(iFreq) = tdata      !     - check value 
   case (iLookSTAT%maxi)                              ! * maximum over period                    
    if (tdata>stat(iFreq)) stat(iFreq) = tdata      !     - check value 
   case (iLookSTAT%mode)                              ! * mode over period (does not workind)       
    stat(iFreq) = realMissing
   case default
    ! message=trim(message)//'unable to identify type of statistic [calculating stats]'
    err=20; return
   ! -------------------------------------------------------------------------------------
  end select
 end do ! looping through output frequencies

 ! ---------------------------------------------
 ! finalize statistics at end of frequency period
 ! ---------------------------------------------
 do iFreq=1,maxVarFreq                                ! loop through output statistics
  if(finalizeStats(iFreq))then
   if(statIndex(iFreq)==integerMissing) cycle     ! don't bother if output frequency is not desired for a given variable
!    if(meta%varType/=iLookVarType%outstat) cycle        ! only calculate stats for scalars
   select case(statIndex(iFreq))                  ! act depending on the statistic
    ! -------------------------------------------------------------------------------------
    case (iLookSTAT%mean)                              ! * mean over period
     stat(iFreq) = stat(iFreq)/statCounter(iFreq)    !     - normalize sum into mean
    case (iLookSTAT%vari)                              ! * variance over period
     stat(maxVarFreq+iFreq) = stat(maxVarFreq+1)/statCounter(iFreq)            ! E[X] term
     stat(iFreq) = stat(iFreq)/statCounter(iFreq) - stat(maxVarFreq+iFreq)**2 ! full variance
    case default ! do nothing -- don't need finalization for most stats
    ! -------------------------------------------------------------------------------------
   end select
  end if
 end do ! looping through output frequencies


 end subroutine calc_stats_h_real

 attributes(device) subroutine calc_stats_d_i4b(data_step,statIndex,stat,data,resetStats,finalizeStats,statCounter,err)
 USE nrtype
 ! data structures
 USE data_types,only:var_info,ilength,dlength ! type dec for meta data structures
 USE var_lookup,only:maxVarFreq       ! # of output frequencies
 ! global variables
!  USE globalData,only:data_step        ! forcing timestep
 ! structures of named variables
 USE var_lookup,only:iLookVarType     ! named variables for variable types
 USE var_lookup,only:iLookFREQ        ! named variables for output frequency
 USE var_lookup,only:iLookSTAT        ! named variables for output statistics
 USE var_lookup,only:iLookTIME        ! named variables for time information
 implicit none
 real(rkind) :: data_step
 ! input variables
!  class(var_info),intent(in)         :: meta              ! meta data structure
 integer(i4b), intent(in) :: statIndex(:)
 real(rkind)       ,intent(inout)      :: stat(:)              ! statistics structure
 integer(i4b)       ,intent(in)         :: data             ! data value
 logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency
 ! output variables
 integer(i4b)   ,intent(out)        :: err               ! error code
!  character(*)   ,intent(out)        :: message           ! error message
 ! internals
!  real(rkind),dimension(maxvarFreq*2)   :: tstat             ! temporary stats vector
 integer(i4b)                       :: iFreq             ! index of output frequency
 real(rkind) :: tdata
 ! initialize error control
 err=0

 tdata = real(data)

 ! ---------------------------------------------
 ! reset statistics at new frequency period
 ! ---------------------------------------------
 if (size(statIndex).eq.1) return
 do iFreq=1,maxVarFreq                              ! loop through output statistics
  if(resetStats(iFreq))then                         ! flag to reset statistics
   if(statIndex(iFreq)==integerMissing) cycle  ! don't bother if output frequency is not desired for a given variable
!    if(meta%varType/=iLookVarType%outstat) cycle     ! only calculate stats for scalars
   select case(statIndex(iFreq))               ! act depending on the statistic
    ! -------------------------------------------------------------------------------------
    case (iLookSTAT%totl)                           ! * summation over period                  
     stat(iFreq) = 0._rkind                           !     - resets stat at beginning of period
    case (iLookSTAT%mean)                           ! * mean over period                       
     stat(iFreq) = 0._rkind                           !     - resets stat at beginning of period
    case (iLookSTAT%vari)                           ! * variance over period                   
     stat(iFreq) = 0._rkind                           !     - resets E[X^2] term in var calc    
     stat(maxVarFreq+iFreq) = 0._rkind                !     - resets E[X]^2 term                 
    case (iLookSTAT%mini)                           ! * minimum over period                    
     stat(iFreq) = huge(stat(iFreq))              !     - resets stat at beginning of period 
    case (iLookSTAT%maxi)                           ! * maximum over period                    
     stat(iFreq) = -huge(stat(iFreq))             !     - resets stat at beginning of period 
    case (iLookSTAT%mode)                           ! * mode over period      
     stat(iFreq) = realMissing                     !     - does not work
    case (iLookSTAT%inst)                           ! * instantaneous -- no need to reset
    case default
    !  message=trim(message)//'unable to identify type of statistic [reset]'
     err=20; return
    ! -------------------------------------------------------------------------------------
   end select
  end if
 end do ! looping through output frequencies

 ! ---------------------------------------------
 ! Calculate each statistic that is requested by user
 ! ---------------------------------------------
 do iFreq=1,maxVarFreq                                ! loop through output statistics
  if(statIndex(iFreq)==integerMissing) cycle     ! don't bother if output frequency is not desired for a given variab;e
!   if(meta%varType/=iLookVarType%outstat) cycle        ! only calculate stats for scalars
  select case(statIndex(iFreq))                  ! act depending on the statistic
   ! -------------------------------------------------------------------------------------
   case (iLookSTAT%inst)                              ! * instantaneous value 
    stat(iFreq) = tdata                              !     - data at a given time
   case (iLookSTAT%totl)                              ! * summation over period                    
    stat(iFreq) = stat(iFreq) + tdata*data_step     !     - increment data 
   case (iLookSTAT%mean)                              ! * mean over period                       
    stat(iFreq) = stat(iFreq) + tdata               !     -  increment data
   case (iLookSTAT%vari)                              ! * variance over period                   
    stat(iFreq) = stat(iFreq) + tdata**2                     ! - E[X^2] term in var calc    
    stat(maxVarFreq+iFreq) = stat(maxVarFreq+iFreq) + tdata  ! - E[X]^2 term                 
   case (iLookSTAT%mini)                              ! * minimum over period                    
    if (tdata<stat(iFreq)) stat(iFreq) = tdata      !     - check value 
   case (iLookSTAT%maxi)                              ! * maximum over period                    
    if (tdata>stat(iFreq)) stat(iFreq) = tdata      !     - check value 
   case (iLookSTAT%mode)                              ! * mode over period (does not workind)       
    stat(iFreq) = realMissing
   case default
    ! message=trim(message)//'unable to identify type of statistic [calculating stats]'
    err=20; return
   ! -------------------------------------------------------------------------------------
  end select
 end do ! looping through output frequencies

 ! ---------------------------------------------
 ! finalize statistics at end of frequency period
 ! ---------------------------------------------
 do iFreq=1,maxVarFreq                                ! loop through output statistics
  if(finalizeStats(iFreq))then
   if(statIndex(iFreq)==integerMissing) cycle     ! don't bother if output frequency is not desired for a given variable
!    if(meta%varType/=iLookVarType%outstat) cycle        ! only calculate stats for scalars
   select case(statIndex(iFreq))                  ! act depending on the statistic
    ! -------------------------------------------------------------------------------------
    case (iLookSTAT%mean)                              ! * mean over period
     stat(iFreq) = stat(iFreq)/statCounter(iFreq)    !     - normalize sum into mean
    case (iLookSTAT%vari)                              ! * variance over period
     stat(maxVarFreq+iFreq) = stat(maxVarFreq+1)/statCounter(iFreq)            ! E[X] term
     stat(iFreq) = stat(iFreq)/statCounter(iFreq) - stat(maxVarFreq+iFreq)**2 ! full variance
    case default ! do nothing -- don't need finalization for most stats
    ! -------------------------------------------------------------------------------------
   end select
  end if
 end do ! looping through output frequencies


 end subroutine calc_stats_d_i4b


end module output_stats
