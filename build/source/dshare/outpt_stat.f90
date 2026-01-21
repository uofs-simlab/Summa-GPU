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
public :: calcStats,calc_stats_kernel
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


 attributes(global) subroutine calc_stats_kernel(nGRU,&
 data_step, &
 statIndex_basin__SurfaceRunoff,&
 stat_basin__SurfaceRunoff,&
 basin__SurfaceRunoff, &
statIndex_basin__ColumnOutflow,&
stat_basin__ColumnOutflow,&
basin__ColumnOutflow, &
statIndex_basin__AquiferStorage,&
stat_basin__AquiferStorage,&
basin__AquiferStorage, &
statIndex_basin__AquiferRecharge,&
stat_basin__AquiferRecharge,&
basin__AquiferRecharge, &
statIndex_basin__AquiferBaseflow,&
stat_basin__AquiferBaseflow,&
basin__AquiferBaseflow, &
statIndex_basin__AquiferTranspire,&
stat_basin__AquiferTranspire,&
basin__AquiferTranspire, &
statIndex_basin__TotalRunoff,&
stat_basin__TotalRunoff,&
basin__TotalRunoff, &
statIndex_basin__SoilDrainage,&
stat_basin__SoilDrainage,&
basin__SoilDrainage, &
statIndex_basin__totalArea,&
stat_basin__totalArea,&
basin__totalArea, &
statIndex_averageInstantRunoff,&
stat_averageInstantRunoff,&
averageInstantRunoff, &
statIndex_averageRoutedRunoff,&
stat_averageRoutedRunoff,&
averageRoutedRunoff,&
resetStats,finalizeStats,statCounter)
    integer(i4b),intent(in),value :: nGRU
    real(rkind),intent(in) :: data_step
    integer(i4b) :: statIndex_basin__SurfaceRunoff(:)
    real(rkind) :: stat_basin__SurfaceRunoff(:)
    real(rkind) :: basin__SurfaceRunoff
    integer(i4b) :: statIndex_basin__ColumnOutflow(:)
    real(rkind) :: stat_basin__ColumnOutflow(:)
    real(rkind) :: basin__ColumnOutflow
    integer(i4b) :: statIndex_basin__AquiferStorage(:)
    real(rkind) :: stat_basin__AquiferStorage(:)
    real(rkind) :: basin__AquiferStorage
    integer(i4b) :: statIndex_basin__AquiferRecharge(:)
    real(rkind) :: stat_basin__AquiferRecharge(:)
    real(rkind) :: basin__AquiferRecharge
    integer(i4b) :: statIndex_basin__AquiferBaseflow(:)
    real(rkind) :: stat_basin__AquiferBaseflow(:)
    real(rkind) :: basin__AquiferBaseflow
    integer(i4b) :: statIndex_basin__AquiferTranspire(:)
    real(rkind) :: stat_basin__AquiferTranspire(:)
    real(rkind) :: basin__AquiferTranspire
    integer(i4b) :: statIndex_basin__TotalRunoff(:)
    real(rkind) :: stat_basin__TotalRunoff(:)
    real(rkind) :: basin__TotalRunoff
    integer(i4b) :: statIndex_basin__SoilDrainage(:)
    real(rkind) :: stat_basin__SoilDrainage(:)
    real(rkind) :: basin__SoilDrainage
    integer(i4b) :: statIndex_basin__totalArea(:)
    real(rkind) :: stat_basin__totalArea(:)
    real(rkind) :: basin__totalArea
    integer(i4b) :: statIndex_averageInstantRunoff(:)
    real(rkind) :: stat_averageInstantRunoff(:)
    real(rkind) :: averageInstantRunoff
    integer(i4b) :: statIndex_averageRoutedRunoff(:)
    real(rkind) :: stat_averageRoutedRunoff(:)
    real(rkind) :: averageRoutedRunoff
     logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency

    call calc_stats_bvar(data_step, &
 statIndex_basin__SurfaceRunoff,&
 stat_basin__SurfaceRunoff,&
 basin__SurfaceRunoff, &
statIndex_basin__ColumnOutflow,&
stat_basin__ColumnOutflow,&
basin__ColumnOutflow, &
statIndex_basin__AquiferStorage,&
stat_basin__AquiferStorage,&
basin__AquiferStorage, &
statIndex_basin__AquiferRecharge,&
stat_basin__AquiferRecharge,&
basin__AquiferRecharge, &
statIndex_basin__AquiferBaseflow,&
stat_basin__AquiferBaseflow,&
basin__AquiferBaseflow, &
statIndex_basin__AquiferTranspire,&
stat_basin__AquiferTranspire,&
basin__AquiferTranspire, &
statIndex_basin__TotalRunoff,&
stat_basin__TotalRunoff,&
basin__TotalRunoff, &
statIndex_basin__SoilDrainage,&
stat_basin__SoilDrainage,&
basin__SoilDrainage, &
statIndex_basin__totalArea,&
stat_basin__totalArea,&
basin__totalArea, &
statIndex_averageInstantRunoff,&
stat_averageInstantRunoff,&
averageInstantRunoff, &
statIndex_averageRoutedRunoff,&
stat_averageRoutedRunoff,&
averageRoutedRunoff,&
resetStats,finalizeStats,statCounter)
end subroutine

 attributes(device) subroutine calc_stats_bvar(data_step, &
 statIndex_basin__SurfaceRunoff,&
 stat_basin__SurfaceRunoff,&
 basin__SurfaceRunoff, &
statIndex_basin__ColumnOutflow,&
stat_basin__ColumnOutflow,&
basin__ColumnOutflow, &
statIndex_basin__AquiferStorage,&
stat_basin__AquiferStorage,&
basin__AquiferStorage, &
statIndex_basin__AquiferRecharge,&
stat_basin__AquiferRecharge,&
basin__AquiferRecharge, &
statIndex_basin__AquiferBaseflow,&
stat_basin__AquiferBaseflow,&
basin__AquiferBaseflow, &
statIndex_basin__AquiferTranspire,&
stat_basin__AquiferTranspire,&
basin__AquiferTranspire, &
statIndex_basin__TotalRunoff,&
stat_basin__TotalRunoff,&
basin__TotalRunoff, &
statIndex_basin__SoilDrainage,&
stat_basin__SoilDrainage,&
basin__SoilDrainage, &
statIndex_basin__totalArea,&
stat_basin__totalArea,&
basin__totalArea, &
statIndex_averageInstantRunoff,&
stat_averageInstantRunoff,&
averageInstantRunoff, &
statIndex_averageRoutedRunoff,&
stat_averageRoutedRunoff,&
averageRoutedRunoff,&
resetStats,finalizeStats,statCounter)
    real(rkind),intent(in) :: data_step
    integer(i4b) :: statIndex_basin__SurfaceRunoff(:)
    real(rkind) :: stat_basin__SurfaceRunoff(:)
    real(rkind) :: basin__SurfaceRunoff
    integer(i4b) :: statIndex_basin__ColumnOutflow(:)
    real(rkind) :: stat_basin__ColumnOutflow(:)
    real(rkind) :: basin__ColumnOutflow
    integer(i4b) :: statIndex_basin__AquiferStorage(:)
    real(rkind) :: stat_basin__AquiferStorage(:)
    real(rkind) :: basin__AquiferStorage
    integer(i4b) :: statIndex_basin__AquiferRecharge(:)
    real(rkind) :: stat_basin__AquiferRecharge(:)
    real(rkind) :: basin__AquiferRecharge
    integer(i4b) :: statIndex_basin__AquiferBaseflow(:)
    real(rkind) :: stat_basin__AquiferBaseflow(:)
    real(rkind) :: basin__AquiferBaseflow
    integer(i4b) :: statIndex_basin__AquiferTranspire(:)
    real(rkind) :: stat_basin__AquiferTranspire(:)
    real(rkind) :: basin__AquiferTranspire
    integer(i4b) :: statIndex_basin__TotalRunoff(:)
    real(rkind) :: stat_basin__TotalRunoff(:)
    real(rkind) :: basin__TotalRunoff
    integer(i4b) :: statIndex_basin__SoilDrainage(:)
    real(rkind) :: stat_basin__SoilDrainage(:)
    real(rkind) :: basin__SoilDrainage
    integer(i4b) :: statIndex_basin__totalArea(:)
    real(rkind) :: stat_basin__totalArea(:)
    real(rkind) :: basin__totalArea
    integer(i4b) :: statIndex_averageInstantRunoff(:)
    real(rkind) :: stat_averageInstantRunoff(:)
    real(rkind) :: averageInstantRunoff
    integer(i4b) :: statIndex_averageRoutedRunoff(:)
    real(rkind) :: stat_averageRoutedRunoff(:)
    real(rkind) :: averageRoutedRunoff
     logical(lgt)   ,intent(in)         :: resetStats(:)     ! vector of flags to reset statistics
 logical(lgt)   ,intent(in)         :: finalizeStats(:)  ! vector of flags to reset statistics
 integer(i4b)   ,intent(in)         :: statCounter(:)   ! number of time steps in each output frequency
integer(i4b) :: err
    !basin__SurfaceRunoff
    call calc_stats_d_real(data_step,statIndex_basin__SurfaceRunoff,stat_basin__SurfaceRunoff,basin__SurfaceRunoff,resetStats,finalizeStats,statCounter,err)
    !basin__ColumnOutflow
    call calc_stats_d_real(data_step,statIndex_basin__ColumnOutflow,stat_basin__ColumnOutflow,basin__ColumnOutflow,resetStats,finalizeStats,statCounter,err)
    !basin__AquiferStorage
    call calc_stats_d_real(data_step,statIndex_basin__AquiferStorage,stat_basin__AquiferStorage,basin__AquiferStorage,resetStats,finalizeStats,statCounter,err)
    !basin__AquiferRecharge
    call calc_stats_d_real(data_step,statIndex_basin__AquiferRecharge,stat_basin__AquiferRecharge,basin__AquiferRecharge,resetStats,finalizeStats,statCounter,err)
    !basin__AquiferBaseflow
    call calc_stats_d_real(data_step,statIndex_basin__AquiferBaseflow,stat_basin__AquiferBaseflow,basin__AquiferBaseflow,resetStats,finalizeStats,statCounter,err)
    !basin__AquiferTranspire
    call calc_stats_d_real(data_step,statIndex_basin__AquiferTranspire,stat_basin__AquiferTranspire,basin__AquiferTranspire,resetStats,finalizeStats,statCounter,err)
    !basin__TotalRunoff
    call calc_stats_d_real(data_step,statIndex_basin__TotalRunoff,stat_basin__TotalRunoff,basin__TotalRunoff,resetStats,finalizeStats,statCounter,err)
    !basin__SoilDrainage
    call calc_stats_d_real(data_step,statIndex_basin__SoilDrainage,stat_basin__SoilDrainage,basin__SoilDrainage,resetStats,finalizeStats,statCounter,err)
    !basin__totalArea
    call calc_stats_d_real(data_step,statIndex_basin__totalArea,stat_basin__totalArea,basin__totalArea,resetStats,finalizeStats,statCounter,err)
    !averageInstantRunoff
    call calc_stats_d_real(data_step,statIndex_averageInstantRunoff,stat_averageInstantRunoff,averageInstantRunoff,resetStats,finalizeStats,statCounter,err)
    !averageRoutedRunoff
    call calc_stats_d_real(data_step,statIndex_averageRoutedRunoff,stat_averageRoutedRunoff,averageRoutedRunoff,resetStats,finalizeStats,statCounter,err)
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


end module output_stats
