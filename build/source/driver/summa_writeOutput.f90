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

module summa_writeOutput
! used to define/write output files

! named variables to define new output files
USE globalData, only: noNewFiles              ! no new output files
USE globalData, only: newFileEveryOct1        ! create a new file on Oct 1 every year (start of the USA water year)

! metadata
USE globalData,only:time_meta                 ! metadata on the model time
USE globalData,only:forc_meta                 ! metadata on the model forcing data
USE globalData,only:diag_meta                 ! metadata on the model diagnostic variables
USE globalData,only:prog_meta                 ! metadata on the model prognostic variables
USE globalData,only:flux_meta                 ! metadata on the model fluxes
USE globalData,only:indx_meta                 ! metadata on the model index variables
USE globalData,only:bvar_meta                 ! metadata on basin-average variables

! child metadata for stats
USE globalData,only:statForc_meta             ! child metadata for stats
USE globalData,only:statProg_meta             ! child metadata for stats
USE globalData,only:statDiag_meta             ! child metadata for stats
USE globalData,only:statFlux_meta             ! child metadata for stats
USE globalData,only:statIndx_meta             ! child metadata for stats
USE globalData,only:statBvar_meta             ! child metadata for stats

! index of the child data structure
USE globalData,only:forcChild_map             ! index of the child data structure: stats forc
USE globalData,only:progChild_map             ! index of the child data structure: stats prog
USE globalData,only:diagChild_map             ! index of the child data structure: stats diag
USE globalData,only:fluxChild_map             ! index of the child data structure: stats flux
USE globalData,only:indxChild_map             ! index of the child data structure: stats indx
USE globalData,only:bvarChild_map             ! index of the child data structure: stats bvar

! named variables
USE var_lookup,only:maxvarFreq                ! maximum number of output files
USE var_lookup,only:iLookTIME                 ! named variables for time data structure
USE var_lookup,only:iLookDIAG                 ! named variables for local column model diagnostic variables
USE var_lookup,only:iLookPROG                 ! named variables for local column model prognostic variables
USE var_lookup,only:iLookINDEX                ! named variables for local column index variables
USE var_lookup,only:iLookFREQ                 ! named variables for the frequency structure

! safety: set private unless specified otherwise
implicit none
private
public::summa_writeOutputFiles
contains

 ! used to define/write output files
 subroutine summa_writeOutputFiles(modelTimeStep, summa1_struc, err, message)
 ! ---------------------------------------------------------------------------------------
 ! * desired modules
 ! ---------------------------------------------------------------------------------------
 ! data types
 USE nrtype                                                  ! variable types, etc.
 USE summa_type,only:summa1_type_dec                         ! master summa data type
 ! subroutines and functions
 USE time_utils_module,only:elapsedSec                       ! calculate the elapsed time
 USE summa_alarms,only:summa_setWriteAlarms                  ! set alarms to control model output
 USE summa_defineOutput,only:summa_defineOutputFiles         ! define summa output files
 USE modelwrite_module,only:writeRestart                     ! module to write model Restart
 USE modelwrite_module,only:writeData,writeBasin,writeBasin_device             ! module to write model output
 USE modelwrite_module,only:writeTime                        ! module to write model time
 USE output_stats,only:calcStats                             ! module for compiling output statistics
  USE output_stats,only:calcStats_d                             ! module for compiling output statistics

 ! global data: general
 USE globalData,only:forcingStep                             ! index of current time step in current forcing file
 USE globalData,only:gru_struc                               ! gru-hru mapping structures
 USE globalData,only:structInfo                              ! information on the data structures
 ! global data: time structures
 USE globalData,only:oldTime                                 ! time from the previous time step
 USE globalData,only:finshTime                               ! end time of simulation
 ! global data: decisions for model alarms
 USE globalData,only:ixProgress                              ! define frequency to write progress
 USE globalData,only:ixRestart                               ! define frequency to write restart files
 USE globalData,only:newOutputFile                           ! define option for new output files
 ! controls on statistics output
 USE globalData,only:statCounter                             ! time counter for stats
 USE globalData,only:resetStats                              ! flags to reset statistics
 USE globalData,only:finalizeStats                           ! flags to finalize statistics
 USE globalData,only:outputTimeStep                          ! timestep in output files
 ! output constraints
 USE globalData,only:maxLayers                               ! maximum number of layers
 USE globalData,only:maxSnowLayers                           ! maximum number of snow layers
 ! timing variables
 USE globalData,only:startWrite,endWrite                     ! date/time for the start and end of the model writing
 USE globalData,only:elapsedWrite                            ! elapsed time to write data
 ! file information
 USE summaFileManager,only:OUTPUT_PATH,OUTPUT_PREFIX         ! define output file
 USE summaFileManager,only:STATE_PATH                        ! optional path to state output files (defaults to OUTPUT_PATH)
 USE globalData,only:output_fileSuffix                       ! suffix for the output & state files (optional summa argument)
 USE globalData,only:nHRUrun                                 ! number of HRU in the run
 USE globalData,only:nGRUrun                                 ! number of GRU in the run
 use device_data_types
 use initialize_device
  USE globalData,only:data_step        ! forcing timestep
USE var_lookup,only:iLookBVAR          ! look-up values for basin-average model variables

 ! ---------------------------------------------------------------------------------------
 ! * variables
 ! ---------------------------------------------------------------------------------------
 implicit none
 ! dummy variables
 integer(i4b),intent(in)               :: modelTimeStep      ! time step index
 type(summa1_type_dec),intent(inout)   :: summa1_struc       ! master summa data structure
 integer(i4b),intent(out)              :: err                ! error code
 character(*),intent(out)              :: message            ! error message
 ! local variables
 character(LEN=256)                    :: cmessage           ! error message of downwind routine
 character(len=256)                    :: timeString                 ! portion of restart file name that contains the write-out time
 character(len=256)                    :: restartFile                ! restart file name
 logical(lgt)                          :: printRestart=.false.       ! flag to print a re-start file
 logical(lgt)                          :: printProgress=.false.      ! flag to print simulation progress
 logical(lgt)                          :: defNewOutputFile=.false.   ! flag to define new output files
 integer(i4b)                          :: iGRU,iHRU          ! indices of GRUs and HRUs
 integer(i4b)                          :: iStruct            ! index of model structure
 integer(i4b)                          :: iFreq              ! index of the output frequency
!  type(bvar_stat_device) :: bvarStat_d
 type(dim3) :: blocks,threads
 integer(i4b),device,allocatable :: statCounter_d(:)
 logical(lgt),device,allocatable :: resetStats_d(:),finalizeStats_d(:)
    threads = dim3(128,1,1)
  blocks = dim3(summa1_struc%nGRU/128+1,1,1)


 ! ---------------------------------------------------------------------------------------
 ! associate to elements in the data structure
 summaVars: associate(&

  ! statistics structures
  forcStat             => summa1_struc%forcStat    , & ! x%gru(:)%hru(:)%var(:)%dat -- model forcing data
  progStat             => summa1_struc%progStat    , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
  diagStat             => summa1_struc%diagStat    , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
  fluxStat             => summa1_struc%fluxStat    , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
  indxStat             => summa1_struc%indxStat    , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
  bvarStat             => summa1_struc%bvarStat    , & ! x%gru(:)%var(:)%dat        -- basin-average variabl

  ! primary data structures
  timeStruct           => summa1_struc%timeStruct  , & ! x%var(:)                   -- model time data
  forcStruct           => summa1_struc%forcStruct  , & ! x%gru(:)%hru(:)%var(:)     -- model forcing data
  indxStruct           => summa1_struc%indxStruct  , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
  progStruct           => summa1_struc%progStruct  , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
  diagStruct           => summa1_struc%diagStruct  , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
  fluxStruct           => summa1_struc%fluxStruct  , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
  bvarStruct           => summa1_struc%bvarStruct  , & ! x%gru(:)%var(:)%dat        -- basin-average variables

  ! miscellaneous variables
  nGRU                 => summa1_struc%nGRU        , & ! number of grouped response units
  nHRU                 => summa1_struc%nHRU          & ! number of global hydrologic response units

 ) ! assignment to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='summa_manageOutputFiles/'

 ! identify the start of the writing
 call date_and_time(values=startWrite)

 ! *****************************************************************************
 ! *** initialize statistics
 ! *****************************************************************************

 ! initialize the statistics flags
 if(modelTimeStep==1)then

  ! initialize time step index
  statCounter(1:maxVarFreq) = 1
  outputTimeStep(1:maxVarFreq) = 1

  ! initialize flags to reset/finalize statistics
  resetStats(:)    = .true.   ! start by resetting statistics
  finalizeStats(:) = .false.  ! do not finalize stats on the first time step

  ! set stats flag for the timestep-level output
  finalizeStats(iLookFREQ%timestep)=.true.

  ! initialize number of hru and gru in global data
  nGRUrun = summa1_struc%nGRU
  nHRUrun = summa1_struc%nHRU
 endif  ! if the first time step

 ! *****************************************************************************
 ! *** set alarms for writing data
 ! *****************************************************************************

 ! set alarms to control model output
 call summa_setWriteAlarms(oldTime%var, timeStruct%var, finshTime%var, &   ! time vectors
                           newOutputFile, defNewOutputFile,            &   ! flag to define new output file
                           ixRestart,     printRestart,                &   ! flag to print the restart file
                           ixProgress,    printProgress,               &   ! flag to print simulation progress
                           resetStats,    finalizeStats,               &   ! flags to reset and finalize stats
                           statCounter,                                &   ! statistics counter
                           err, cmessage)                                  ! error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! print progress
 if(printProgress) write(*,'(i4,1x,5(i2,1x))') summa1_struc%timeStruct%var(1:5)

 ! *****************************************************************************
 ! *** define summa output files
 ! *****************************************************************************

 ! check the need to create a new output file
 if(defNewOutputFile .or. modelTimeStep==1)then

  ! define summa output files
  call summa_defineOutputFiles(modelTimeStep, summa1_struc, err, cmessage)
  ! if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! re-initalize the indices for model writing
  outputTimeStep(:)=1

 end if  ! if defining a new file

 ! ****************************************************************************
 ! *** calculate output statistics
 ! ****************************************************************************
statCounter_d = statCounter
 resetStats_d = resetStats
 finalizeStats_d = finalizeStats
 
 ! loop through GRUs and HRUs
!  do iGRU=1,summa1_struc%nGRU
  ! do iHRU=1,gru_struc(iGRU)%hruCount

    !call calcStats(summa1_struc%indxStat%gru(iGRU)%hru(iHRU)%var,summa1_struc%indxStruct%gru(iGRU)%hru(iHRU)%var,statIndx_meta,resetStats,finalizeStats,statCounter,err,cmessage)
   ! calculate output Statistics
   do iStruct=1,size(structInfo)
    select case(trim(structInfo(iStruct)%structName))
     case('forc'); call calcStats_d(nGRU,forcStat,forcStruct,statForc_meta,resetStats_d,finalizeStats_d,statCounter_d,err,cmessage)
     case('prog'); call calcStats_d(nGRU,progStat,progStruct,statProg_meta,resetStats_d,finalizeStats_d,statCounter_d,err,cmessage)
     case('diag'); call calcStats_d(nGRU,diagStat,diagStruct,statDiag_meta,resetStats_d,finalizeStats_d,statCounter_d,err,cmessage)
    !  case('flux'); call calcStats_d(nGRU,fluxStat,fluxStruct,statFlux_meta,resetStats_d,finalizeStats_d,statCounter_d,err,cmessage)
    !  case('indx'); call calcStats_d(nGRU,indxStat,indxStruct,statIndx_meta,resetStats_d,finalizeStats_d,statCounter_d,err,cmessage)
    end select
    if(err/=0)then; message=trim(message)//trim(cmessage)//'['//trim(structInfo(iStruct)%structName)//']'; return; endif
   end do  ! (looping through structures)

  ! end do  ! (looping through HRUs)

  ! calc basin stats
   call calcStats_d(nGRU,bvarStat,bvarStruct,statBvar_meta,resetStats_d,finalizeStats_d,statCounter_d,err,message)
  ! call calcStats(bvarStat%gru(iGRU)%var,bvarStruct%gru(iGRU)%var,statBvar_meta,resetStats,finalizeStats,statCounter,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage)//'[bvar stats]'; return; endif


!  end do  ! (looping through GRUs)

call calc_stats_kernel_indx<<<blocks,threads>>>(nGRU, data_step, &
 indxStat%statIndex_nSnow,indxStat%stat_nSnow,indxStruct%nSnow, &
indxStat%statIndex_nSoil,indxStat%stat_nSoil,indxStruct%nSoil, &
indxStat%statIndex_nLayers,indxStat%stat_nLayers,indxStruct%nLayers, &
indxStat%statIndex_nCasNrg,indxStat%stat_nCasNrg,indxStruct%nCasNrg, &
indxStat%statIndex_nVegNrg,indxStat%stat_nVegNrg,indxStruct%nVegNrg, &
indxStat%statIndex_nVegMass,indxStat%stat_nVegMass,indxStruct%nVegMass, &
indxStat%statIndex_nVegState,indxStat%stat_nVegState,indxStruct%nVegState, &
indxStat%statIndex_nNrgState,indxStat%stat_nNrgState,indxStruct%nNrgState, &
indxStat%statIndex_nWatState,indxStat%stat_nWatState,indxStruct%nWatState, &
indxStat%statIndex_nMatState,indxStat%stat_nMatState,indxStruct%nMatState, &
indxStat%statIndex_nMassState,indxStat%stat_nMassState,indxStruct%nMassState, &
indxStat%statIndex_nState,indxStat%stat_nState,indxStruct%nState, &
indxStat%statIndex_ixCasNrg,indxStat%stat_ixCasNrg,indxStruct%ixCasNrg, &
indxStat%statIndex_ixVegNrg,indxStat%stat_ixVegNrg,indxStruct%ixVegNrg, &
indxStat%statIndex_ixVegHyd,indxStat%stat_ixVegHyd,indxStruct%ixVegHyd, &
indxStat%statIndex_ixTopNrg,indxStat%stat_ixTopNrg,indxStruct%ixTopNrg, &
indxStat%statIndex_ixTopHyd,indxStat%stat_ixTopHyd,indxStruct%ixTopHyd, &
indxStat%statIndex_ixAqWat,indxStat%stat_ixAqWat,indxStruct%ixAqWat, &
indxStat%statIndex_ixNrgCanair,indxStat%stat_ixNrgCanair,indxStruct%ixNrgCanair, &
indxStat%statIndex_ixNrgCanopy,indxStat%stat_ixNrgCanopy,indxStruct%ixNrgCanopy, &
indxStat%statIndex_ixHydCanopy,indxStat%stat_ixHydCanopy,indxStruct%ixHydCanopy, &
indxStat%statIndex_ixWatAquifer,indxStat%stat_ixWatAquifer,indxStruct%ixWatAquifer, &
indxStat%statIndex_numberFluxCalc,indxStat%stat_numberFluxCalc,indxStruct%numberFluxCalc, &
indxStat%statIndex_numberStateSplit,indxStat%stat_numberStateSplit,indxStruct%numberStateSplit, &
indxStat%statIndex_numberDomainSplitNrg,indxStat%stat_numberDomainSplitNrg,indxStruct%numberDomainSplitNrg, &
indxStat%statIndex_numberDomainSplitMass,indxStat%stat_numberDomainSplitMass,indxStruct%numberDomainSplitMass, &
indxStat%statIndex_numberScalarSolutions,indxStat%stat_numberScalarSolutions,indxStruct%numberScalarSolutions, &
resetStats_d,finalizeStats_d,statCounter_d)

call calc_stats_kernel_flux<<<blocks,threads>>>(nGRU,data_step,fluxStruct%data, &
fluxStat%statIndex_scalarCanairNetNrgFlux,fluxStat%stat_scalarCanairNetNrgFlux,fluxStruct%ixscalarCanairNetNrgFlux, &
fluxStat%statIndex_scalarCanopyNetNrgFlux,fluxStat%stat_scalarCanopyNetNrgFlux,fluxStruct%ixscalarCanopyNetNrgFlux, &
fluxStat%statIndex_scalarGroundNetNrgFlux,fluxStat%stat_scalarGroundNetNrgFlux,fluxStruct%ixscalarGroundNetNrgFlux, &
fluxStat%statIndex_scalarCanopyNetLiqFlux,fluxStat%stat_scalarCanopyNetLiqFlux,fluxStruct%ixscalarCanopyNetLiqFlux, &
fluxStat%statIndex_scalarRainfall,fluxStat%stat_scalarRainfall,fluxStruct%ixscalarRainfall, &
fluxStat%statIndex_scalarSnowfall,fluxStat%stat_scalarSnowfall,fluxStruct%ixscalarSnowfall, &
fluxStat%statIndex_scalarCanopySunlitPAR,fluxStat%stat_scalarCanopySunlitPAR,fluxStruct%ixscalarCanopySunlitPAR, &
fluxStat%statIndex_scalarCanopyShadedPAR,fluxStat%stat_scalarCanopyShadedPAR,fluxStruct%ixscalarCanopyShadedPAR, &
fluxStat%statIndex_scalarBelowCanopySolar,fluxStat%stat_scalarBelowCanopySolar,fluxStruct%ixscalarBelowCanopySolar, &
fluxStat%statIndex_scalarCanopyAbsorbedSolar,fluxStat%stat_scalarCanopyAbsorbedSolar,fluxStruct%ixscalarCanopyAbsorbedSolar, &
fluxStat%statIndex_scalarGroundAbsorbedSolar,fluxStat%stat_scalarGroundAbsorbedSolar,fluxStruct%ixscalarGroundAbsorbedSolar, &
fluxStat%statIndex_scalarLWRadCanopy,fluxStat%stat_scalarLWRadCanopy,fluxStruct%ixscalarLWRadCanopy, &
fluxStat%statIndex_scalarLWRadGround,fluxStat%stat_scalarLWRadGround,fluxStruct%ixscalarLWRadGround, &
fluxStat%statIndex_scalarLWRadUbound2Canopy,fluxStat%stat_scalarLWRadUbound2Canopy,fluxStruct%ixscalarLWRadUbound2Canopy, &
fluxStat%statIndex_scalarLWRadUbound2Ground,fluxStat%stat_scalarLWRadUbound2Ground,fluxStruct%ixscalarLWRadUbound2Ground, &
fluxStat%statIndex_scalarLWRadUbound2Ubound,fluxStat%stat_scalarLWRadUbound2Ubound,fluxStruct%ixscalarLWRadUbound2Ubound, &
fluxStat%statIndex_scalarLWRadCanopy2Ubound,fluxStat%stat_scalarLWRadCanopy2Ubound,fluxStruct%ixscalarLWRadCanopy2Ubound, &
fluxStat%statIndex_scalarLWRadCanopy2Ground,fluxStat%stat_scalarLWRadCanopy2Ground,fluxStruct%ixscalarLWRadCanopy2Ground, &
fluxStat%statIndex_scalarLWRadCanopy2Canopy,fluxStat%stat_scalarLWRadCanopy2Canopy,fluxStruct%ixscalarLWRadCanopy2Canopy, &
fluxStat%statIndex_scalarLWRadGround2Ubound,fluxStat%stat_scalarLWRadGround2Ubound,fluxStruct%ixscalarLWRadGround2Ubound, &
fluxStat%statIndex_scalarLWRadGround2Canopy,fluxStat%stat_scalarLWRadGround2Canopy,fluxStruct%ixscalarLWRadGround2Canopy, &
fluxStat%statIndex_scalarLWNetCanopy,fluxStat%stat_scalarLWNetCanopy,fluxStruct%ixscalarLWNetCanopy, &
fluxStat%statIndex_scalarLWNetGround,fluxStat%stat_scalarLWNetGround,fluxStruct%ixscalarLWNetGround, &
fluxStat%statIndex_scalarLWNetUbound,fluxStat%stat_scalarLWNetUbound,fluxStruct%ixscalarLWNetUbound, &
fluxStat%statIndex_scalarEddyDiffusCanopyTop,fluxStat%stat_scalarEddyDiffusCanopyTop,fluxStruct%ixscalarEddyDiffusCanopyTop, &
fluxStat%statIndex_scalarFrictionVelocity,fluxStat%stat_scalarFrictionVelocity,fluxStruct%ixscalarFrictionVelocity, &
fluxStat%statIndex_scalarWindspdCanopyTop,fluxStat%stat_scalarWindspdCanopyTop,fluxStruct%ixscalarWindspdCanopyTop, &
fluxStat%statIndex_scalarWindspdCanopyBottom,fluxStat%stat_scalarWindspdCanopyBottom,fluxStruct%ixscalarWindspdCanopyBottom, &
fluxStat%statIndex_scalarGroundResistance,fluxStat%stat_scalarGroundResistance,fluxStruct%ixscalarGroundResistance, &
fluxStat%statIndex_scalarCanopyResistance,fluxStat%stat_scalarCanopyResistance,fluxStruct%ixscalarCanopyResistance, &
fluxStat%statIndex_scalarLeafResistance,fluxStat%stat_scalarLeafResistance,fluxStruct%ixscalarLeafResistance, &
fluxStat%statIndex_scalarSoilResistance,fluxStat%stat_scalarSoilResistance,fluxStruct%ixscalarSoilResistance, &
fluxStat%statIndex_scalarSenHeatTotal,fluxStat%stat_scalarSenHeatTotal,fluxStruct%ixscalarSenHeatTotal, &
fluxStat%statIndex_scalarSenHeatCanopy,fluxStat%stat_scalarSenHeatCanopy,fluxStruct%ixscalarSenHeatCanopy, &
fluxStat%statIndex_scalarSenHeatGround,fluxStat%stat_scalarSenHeatGround,fluxStruct%ixscalarSenHeatGround, &
fluxStat%statIndex_scalarLatHeatTotal,fluxStat%stat_scalarLatHeatTotal,fluxStruct%ixscalarLatHeatTotal, &
fluxStat%statIndex_scalarLatHeatCanopyEvap,fluxStat%stat_scalarLatHeatCanopyEvap,fluxStruct%ixscalarLatHeatCanopyEvap, &
fluxStat%statIndex_scalarLatHeatCanopyTrans,fluxStat%stat_scalarLatHeatCanopyTrans,fluxStruct%ixscalarLatHeatCanopyTrans, &
fluxStat%statIndex_scalarLatHeatGround,fluxStat%stat_scalarLatHeatGround,fluxStruct%ixscalarLatHeatGround, &
fluxStat%statIndex_scalarCanopyAdvectiveHeatFlux,fluxStat%stat_scalarCanopyAdvectiveHeatFlux,fluxStruct%ixscalarCanopyAdvectiveHeatFlux, &
fluxStat%statIndex_scalarGroundAdvectiveHeatFlux,fluxStat%stat_scalarGroundAdvectiveHeatFlux,fluxStruct%ixscalarGroundAdvectiveHeatFlux, &
fluxStat%statIndex_scalarCanopySublimation,fluxStat%stat_scalarCanopySublimation,fluxStruct%ixscalarCanopySublimation, &
fluxStat%statIndex_scalarSnowSublimation,fluxStat%stat_scalarSnowSublimation,fluxStruct%ixscalarSnowSublimation, &
fluxStat%statIndex_scalarStomResistSunlit,fluxStat%stat_scalarStomResistSunlit,fluxStruct%ixscalarStomResistSunlit, &
fluxStat%statIndex_scalarStomResistShaded,fluxStat%stat_scalarStomResistShaded,fluxStruct%ixscalarStomResistShaded, &
fluxStat%statIndex_scalarPhotosynthesisSunlit,fluxStat%stat_scalarPhotosynthesisSunlit,fluxStruct%ixscalarPhotosynthesisSunlit, &
fluxStat%statIndex_scalarPhotosynthesisShaded,fluxStat%stat_scalarPhotosynthesisShaded,fluxStruct%ixscalarPhotosynthesisShaded, &
fluxStat%statIndex_scalarCanopyTranspiration,fluxStat%stat_scalarCanopyTranspiration,fluxStruct%ixscalarCanopyTranspiration, &
fluxStat%statIndex_scalarCanopyEvaporation,fluxStat%stat_scalarCanopyEvaporation,fluxStruct%ixscalarCanopyEvaporation, &
fluxStat%statIndex_scalarGroundEvaporation,fluxStat%stat_scalarGroundEvaporation,fluxStruct%ixscalarGroundEvaporation, &
fluxStat%statIndex_scalarThroughfallSnow,fluxStat%stat_scalarThroughfallSnow,fluxStruct%ixscalarThroughfallSnow, &
fluxStat%statIndex_scalarThroughfallRain,fluxStat%stat_scalarThroughfallRain,fluxStruct%ixscalarThroughfallRain, &
fluxStat%statIndex_scalarCanopySnowUnloading,fluxStat%stat_scalarCanopySnowUnloading,fluxStruct%ixscalarCanopySnowUnloading, &
fluxStat%statIndex_scalarCanopyLiqDrainage,fluxStat%stat_scalarCanopyLiqDrainage,fluxStruct%ixscalarCanopyLiqDrainage, &
fluxStat%statIndex_scalarCanopyMeltFreeze,fluxStat%stat_scalarCanopyMeltFreeze,fluxStruct%ixscalarCanopyMeltFreeze, &
fluxStat%statIndex_scalarSnowDrainage,fluxStat%stat_scalarSnowDrainage,fluxStruct%ixscalarSnowDrainage, &
fluxStat%statIndex_scalarRainPlusMelt,fluxStat%stat_scalarRainPlusMelt,fluxStruct%ixscalarRainPlusMelt, &
fluxStat%statIndex_scalarMaxInfilRate,fluxStat%stat_scalarMaxInfilRate,fluxStruct%ixscalarMaxInfilRate, &
fluxStat%statIndex_scalarInfiltration,fluxStat%stat_scalarInfiltration,fluxStruct%ixscalarInfiltration, &
fluxStat%statIndex_scalarExfiltration,fluxStat%stat_scalarExfiltration,fluxStruct%ixscalarExfiltration, &
fluxStat%statIndex_scalarSurfaceRunoff,fluxStat%stat_scalarSurfaceRunoff,fluxStruct%ixscalarSurfaceRunoff, &
fluxStat%statIndex_scalarSurfaceRunoff_IE,fluxStat%stat_scalarSurfaceRunoff_IE,fluxStruct%ixscalarSurfaceRunoff_IE, &
fluxStat%statIndex_scalarSurfaceRunoff_SE,fluxStat%stat_scalarSurfaceRunoff_SE,fluxStruct%ixscalarSurfaceRunoff_SE, &
fluxStat%statIndex_scalarSoilBaseflow,fluxStat%stat_scalarSoilBaseflow,fluxStruct%ixscalarSoilBaseflow, &
fluxStat%statIndex_scalarSoilDrainage,fluxStat%stat_scalarSoilDrainage,fluxStruct%ixscalarSoilDrainage, &
fluxStat%statIndex_scalarAquiferRecharge,fluxStat%stat_scalarAquiferRecharge,fluxStruct%ixscalarAquiferRecharge, &
fluxStat%statIndex_scalarAquiferTranspire,fluxStat%stat_scalarAquiferTranspire,fluxStruct%ixscalarAquiferTranspire, &
fluxStat%statIndex_scalarAquiferBaseflow,fluxStat%stat_scalarAquiferBaseflow,fluxStruct%ixscalarAquiferBaseflow, &
fluxStat%statIndex_scalarTotalET,fluxStat%stat_scalarTotalET,fluxStruct%ixscalarTotalET, &
fluxStat%statIndex_scalarTotalRunoff,fluxStat%stat_scalarTotalRunoff,fluxStruct%ixscalarTotalRunoff, &
fluxStat%statIndex_scalarNetRadiation,fluxStat%stat_scalarNetRadiation,fluxStruct%ixscalarNetRadiation, &
resetStats_d,finalizeStats_d,statCounter_d)

! call calc_stats_kernel_diag<<<blocks,threads>>>(nGRU,data_step, &
! diagStat%statIndex_scalarCanopyDepth,diagStat%stat_scalarCanopyDepth,diagStruct%scalarCanopyDepth, &
! diagStat%statIndex_scalarBulkVolHeatCapVeg,diagStat%stat_scalarBulkVolHeatCapVeg,diagStruct%scalarBulkVolHeatCapVeg, &
! diagStat%statIndex_scalarCanopyCm,diagStat%stat_scalarCanopyCm,diagStruct%scalarCanopyCm, &
! diagStat%statIndex_scalarCanopyEmissivity,diagStat%stat_scalarCanopyEmissivity,diagStruct%scalarCanopyEmissivity, &
! diagStat%statIndex_scalarRootZoneTemp,diagStat%stat_scalarRootZoneTemp,diagStruct%scalarRootZoneTemp, &
! diagStat%statIndex_scalarLAI,diagStat%stat_scalarLAI,diagStruct%scalarLAI, &
! diagStat%statIndex_scalarSAI,diagStat%stat_scalarSAI,diagStruct%scalarSAI, &
! diagStat%statIndex_scalarExposedLAI,diagStat%stat_scalarExposedLAI,diagStruct%scalarExposedLAI, &
! diagStat%statIndex_scalarExposedSAI,diagStat%stat_scalarExposedSAI,diagStruct%scalarExposedSAI, &
! diagStat%statIndex_scalarAdjMeasHeight,diagStat%stat_scalarAdjMeasHeight,diagStruct%scalarAdjMeasHeight, &
! diagStat%statIndex_scalarCanopyIceMax,diagStat%stat_scalarCanopyIceMax,diagStruct%scalarCanopyIceMax, &
! diagStat%statIndex_scalarCanopyLiqMax,diagStat%stat_scalarCanopyLiqMax,diagStruct%scalarCanopyLiqMax, &
! diagStat%statIndex_scalarGrowingSeasonIndex,diagStat%stat_scalarGrowingSeasonIndex,diagStruct%scalarGrowingSeasonIndex, &
! diagStat%statIndex_scalarVolHtCap_air,diagStat%stat_scalarVolHtCap_air,diagStruct%scalarVolHtCap_air, &
! diagStat%statIndex_scalarVolHtCap_ice,diagStat%stat_scalarVolHtCap_ice,diagStruct%scalarVolHtCap_ice, &
! diagStat%statIndex_scalarVolHtCap_soil,diagStat%stat_scalarVolHtCap_soil,diagStruct%scalarVolHtCap_soil, &
! diagStat%statIndex_scalarVolHtCap_water,diagStat%stat_scalarVolHtCap_water,diagStruct%scalarVolHtCap_water, &
! diagStat%statIndex_scalarLambda_drysoil,diagStat%stat_scalarLambda_drysoil,diagStruct%scalarLambda_drysoil, &
! diagStat%statIndex_scalarLambda_wetsoil,diagStat%stat_scalarLambda_wetsoil,diagStruct%scalarLambda_wetsoil, &
! diagStat%statIndex_scalarCanopyEnthTemp,diagStat%stat_scalarCanopyEnthTemp,diagStruct%scalarCanopyEnthTemp, &
! diagStat%statIndex_scalarTotalSoilEnthalpy,diagStat%stat_scalarTotalSoilEnthalpy,diagStruct%scalarTotalSoilEnthalpy, &
! diagStat%statIndex_scalarTotalSnowEnthalpy,diagStat%stat_scalarTotalSnowEnthalpy,diagStruct%scalarTotalSnowEnthalpy, &
! diagStat%statIndex_scalarVPair,diagStat%stat_scalarVPair,diagStruct%scalarVPair, &
! diagStat%statIndex_scalarVP_CanopyAir,diagStat%stat_scalarVP_CanopyAir,diagStruct%scalarVP_CanopyAir, &
! diagStat%statIndex_scalarTwetbulb,diagStat%stat_scalarTwetbulb,diagStruct%scalarTwetbulb, &
! diagStat%statIndex_scalarSnowfallTemp,diagStat%stat_scalarSnowfallTemp,diagStruct%scalarSnowfallTemp, &
! diagStat%statIndex_scalarNewSnowDensity,diagStat%stat_scalarNewSnowDensity,diagStruct%scalarNewSnowDensity, &
! diagStat%statIndex_scalarO2air,diagStat%stat_scalarO2air,diagStruct%scalarO2air, &
! diagStat%statIndex_scalarCO2air,diagStat%stat_scalarCO2air,diagStruct%scalarCO2air, &
! diagStat%statIndex_windspd_x,diagStat%stat_windspd_x,diagStruct%windspd_x, &
! diagStat%statIndex_windspd_y,diagStat%stat_windspd_y,diagStruct%windspd_y, &
! diagStat%statIndex_scalarCosZenith,diagStat%stat_scalarCosZenith,diagStruct%scalarCosZenith, &
! diagStat%statIndex_scalarFractionDirect,diagStat%stat_scalarFractionDirect,diagStruct%scalarFractionDirect, &
! diagStat%statIndex_scalarCanopySunlitFraction,diagStat%stat_scalarCanopySunlitFraction,diagStruct%scalarCanopySunlitFraction, &
! diagStat%statIndex_scalarCanopySunlitLAI,diagStat%stat_scalarCanopySunlitLAI,diagStruct%scalarCanopySunlitLAI, &
! diagStat%statIndex_scalarCanopyShadedLAI,diagStat%stat_scalarCanopyShadedLAI,diagStruct%scalarCanopyShadedLAI, &
! diagStat%statIndex_scalarGroundAlbedo,diagStat%stat_scalarGroundAlbedo,diagStruct%scalarGroundAlbedo, &
! diagStat%statIndex_scalarLatHeatSubVapCanopy,diagStat%stat_scalarLatHeatSubVapCanopy,diagStruct%scalarLatHeatSubVapCanopy, &
! diagStat%statIndex_scalarLatHeatSubVapGround,diagStat%stat_scalarLatHeatSubVapGround,diagStruct%scalarLatHeatSubVapGround, &
! diagStat%statIndex_scalarSatVP_CanopyTemp,diagStat%stat_scalarSatVP_CanopyTemp,diagStruct%scalarSatVP_CanopyTemp, &
! diagStat%statIndex_scalarSatVP_GroundTemp,diagStat%stat_scalarSatVP_GroundTemp,diagStruct%scalarSatVP_GroundTemp, &
! diagStat%statIndex_scalarZ0Canopy,diagStat%stat_scalarZ0Canopy,diagStruct%scalarZ0Canopy, &
! diagStat%statIndex_scalarWindReductionFactor,diagStat%stat_scalarWindReductionFactor,diagStruct%scalarWindReductionFactor, &
! diagStat%statIndex_scalarZeroPlaneDisplacement,diagStat%stat_scalarZeroPlaneDisplacement,diagStruct%scalarZeroPlaneDisplacement, &
! diagStat%statIndex_scalarRiBulkCanopy,diagStat%stat_scalarRiBulkCanopy,diagStruct%scalarRiBulkCanopy, &
! diagStat%statIndex_scalarRiBulkGround,diagStat%stat_scalarRiBulkGround,diagStruct%scalarRiBulkGround, &
! diagStat%statIndex_scalarCanopyStabilityCorrection,diagStat%stat_scalarCanopyStabilityCorrection,diagStruct%scalarCanopyStabilityCorrection, &
! diagStat%statIndex_scalarGroundStabilityCorrection,diagStat%stat_scalarGroundStabilityCorrection,diagStruct%scalarGroundStabilityCorrection, &
! diagStat%statIndex_scalarIntercellularCO2Sunlit,diagStat%stat_scalarIntercellularCO2Sunlit,diagStruct%scalarIntercellularCO2Sunlit, &
! diagStat%statIndex_scalarIntercellularCO2Shaded,diagStat%stat_scalarIntercellularCO2Shaded,diagStruct%scalarIntercellularCO2Shaded, &
! diagStat%statIndex_scalarTranspireLim,diagStat%stat_scalarTranspireLim,diagStruct%scalarTranspireLim, &
! diagStat%statIndex_scalarTranspireLimAqfr,diagStat%stat_scalarTranspireLimAqfr,diagStruct%scalarTranspireLimAqfr, &
! diagStat%statIndex_scalarFoliageNitrogenFactor,diagStat%stat_scalarFoliageNitrogenFactor,diagStruct%scalarFoliageNitrogenFactor, &
! diagStat%statIndex_scalarSoilRelHumidity,diagStat%stat_scalarSoilRelHumidity,diagStruct%scalarSoilRelHumidity, &
! diagStat%statIndex_scalarAquiferRootFrac,diagStat%stat_scalarAquiferRootFrac,diagStruct%scalarAquiferRootFrac, &
! diagStat%statIndex_scalarFracLiqVeg,diagStat%stat_scalarFracLiqVeg,diagStruct%scalarFracLiqVeg, &
! diagStat%statIndex_scalarCanopyWetFraction,diagStat%stat_scalarCanopyWetFraction,diagStruct%scalarCanopyWetFraction, &
! diagStat%statIndex_scalarSnowAge,diagStat%stat_scalarSnowAge,diagStruct%scalarSnowAge, &
! diagStat%statIndex_scalarGroundSnowFraction,diagStat%stat_scalarGroundSnowFraction,diagStruct%scalarGroundSnowFraction, &
! diagStat%statIndex_scalarInfilArea,diagStat%stat_scalarInfilArea,diagStruct%scalarInfilArea, &
! diagStat%statIndex_scalarFrozenArea,diagStat%stat_scalarFrozenArea,diagStruct%scalarFrozenArea, &
! diagStat%statIndex_scalarSoilControl,diagStat%stat_scalarSoilControl,diagStruct%scalarSoilControl, &
! diagStat%statIndex_scalarSoilCompress,diagStat%stat_scalarSoilCompress,diagStruct%scalarSoilCompress, &
! diagStat%statIndex_scalarTotalSoilLiq,diagStat%stat_scalarTotalSoilLiq,diagStruct%scalarTotalSoilLiq, &
! diagStat%statIndex_scalarTotalSoilIce,diagStat%stat_scalarTotalSoilIce,diagStruct%scalarTotalSoilIce, &
! diagStat%statIndex_scalarTotalSoilWat,diagStat%stat_scalarTotalSoilWat,diagStruct%scalarTotalSoilWat, &
! diagStat%statIndex_scalarKappa,diagStat%stat_scalarKappa,diagStruct%scalarKappa, &
! diagStat%statIndex_scalarVolLatHt_fus,diagStat%stat_scalarVolLatHt_fus,diagStruct%scalarVolLatHt_fus, &
! diagStat%statIndex_numFluxCalls,diagStat%stat_numFluxCalls,diagStruct%numFluxCalls, &
! diagStat%statIndex_wallClockTime,diagStat%stat_wallClockTime,diagStruct%wallClockTime, &
! diagStat%statIndex_meanStepSize,diagStat%stat_meanStepSize,diagStruct%meanStepSize, &
! diagStat%statIndex_balanceCasNrg,diagStat%stat_balanceCasNrg,diagStruct%balanceCasNrg, &
! diagStat%statIndex_balanceVegNrg,diagStat%stat_balanceVegNrg,diagStruct%balanceVegNrg, &
! diagStat%statIndex_balanceSnowNrg,diagStat%stat_balanceSnowNrg,diagStruct%balanceSnowNrg, &
! diagStat%statIndex_balanceSoilNrg,diagStat%stat_balanceSoilNrg,diagStruct%balanceSoilNrg, &
! diagStat%statIndex_balanceVegMass,diagStat%stat_balanceVegMass,diagStruct%balanceVegMass, &
! diagStat%statIndex_balanceSnowMass,diagStat%stat_balanceSnowMass,diagStruct%balanceSnowMass, &
! diagStat%statIndex_balanceSoilMass,diagStat%stat_balanceSoilMass,diagStruct%balanceSoilMass, &
! diagStat%statIndex_balanceAqMass,diagStat%stat_balanceAqMass,diagStruct%balanceAqMass, &
! diagStat%statIndex_numSteps,diagStat%stat_numSteps,diagStruct%numSteps, &
! diagStat%statIndex_numResEvals,diagStat%stat_numResEvals,diagStruct%numResEvals, &
! diagStat%statIndex_numLinSolvSetups,diagStat%stat_numLinSolvSetups,diagStruct%numLinSolvSetups, &
! diagStat%statIndex_numErrTestFails,diagStat%stat_numErrTestFails,diagStruct%numErrTestFails, &
! resetStats_d,finalizeStats_d,statCounter_d)

! call calc_stats_kernel_diag_part2<<<blocks,threads>>>(nGRU,data_step, &
! diagStat%statIndex_kLast,diagStat%stat_kLast,diagStruct%kLast, &
! diagStat%statIndex_kCur,diagStat%stat_kCur,diagStruct%kCur, &
! diagStat%statIndex_hInitUsed,diagStat%stat_hInitUsed,diagStruct%hInitUsed, &
! diagStat%statIndex_hLast,diagStat%stat_hLast,diagStruct%hLast, &
! diagStat%statIndex_hCur,diagStat%stat_hCur,diagStruct%hCur, &
! diagStat%statIndex_tCur,diagStat%stat_tCur,diagStruct%tCur, &
! resetStats_d,finalizeStats_d,statCounter_d)

! call calc_stats_kernel_forc<<<blocks,threads>>>(nGRU,data_step,&
!     forcStat%statIndex_time,forcStat%stat_time,forcStruct%time_h,&
! forcStat%statIndex_pptrate,forcStat%stat_pptrate,forcStruct%pptrate,&
! forcStat%statIndex_SWRadAtm,forcStat%stat_SWRadAtm,forcStruct%SWRadAtm,&
! forcStat%statIndex_LWRadAtm,forcStat%stat_LWRadAtm,forcStruct%LWRadAtm,&
! forcStat%statIndex_airtemp,forcStat%stat_airtemp,forcStruct%airtemp,&
! forcStat%statIndex_windspd,forcStat%stat_windspd,forcStruct%windspd,&
! forcStat%statIndex_airpres,forcStat%stat_airpres,forcStruct%airpres,&
! forcStat%statIndex_spechum,forcStat%stat_spechum,forcStruct%spechum,&
! resetStats_d,finalizeStats_d,statCounter_d)

!  do iGRU=1,summa1_struc%nGRU

  ! write basin-average variables
  call writeBasin_device(summa1_struc%nGRU,finalizeStats,outputTimeStep,bvar_meta,bvarStat,bvarStruct,bvarChild_map,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage)//'[bvar]'; return; endif
!  end do

 ! ****************************************************************************
 ! *** write data
 ! ****************************************************************************

 ! get the number of HRUs in the run domain
 nHRUrun = sum(gru_struc%hruCount)

 ! write time information
 call writeTime(finalizeStats,outputTimeStep,time_meta,timeStruct%var,err,message)

 ! write the model output to the NetCDF file
 ! Passes the full metadata structure rather than the stats metadata structure because
 !  we have the option to write out data of types other than statistics.
 !  Thus, we must also pass the stats parent->child maps from childStruct.
 do iStruct=1,size(structInfo)
  select case(trim(structInfo(iStruct)%structName))
   case('forc'); call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,forc_meta,forcStat,forcStruct,forcChild_map,indxStruct,err,cmessage)
   case('prog'); call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,prog_meta,progStat,progStruct,progChild_map,indxStruct,err,cmessage)
   case('diag'); call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,diag_meta,diagStat,diagStruct,diagChild_map,indxStruct,err,cmessage)
   case('flux'); call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,flux_meta,fluxStat,fluxStruct,fluxChild_map,indxStruct,err,cmessage)
   case('indx'); call writeData(finalizeStats,outputTimeStep,nHRUrun,maxLayers,indx_meta,indxStat,indxStruct,indxChild_map,indxStruct,err,cmessage)
  end select
  if(err/=0)then; message=trim(message)//trim(cmessage)//'['//trim(structInfo(iStruct)%structName)//']'; return; endif
 end do  ! (looping through structures)
 ! *****************************************************************************
 ! *** write restart file
 ! *****************************************************************************

 ! print a restart file if requested
 if(printRestart)then
  write(timeString,'(i4,3(i2.2))') timeStruct%var(iLookTIME%iyyy),timeStruct%var(iLookTIME%im),timeStruct%var(iLookTIME%id),timeStruct%var(iLookTIME%ih)
  
  if(STATE_PATH == '') then
    restartFile=trim(OUTPUT_PATH)//trim(OUTPUT_PREFIX)//'_restart_'//trim(timeString)//trim(output_fileSuffix)//'.nc'
  else
    restartFile=trim(STATE_PATH)//trim(OUTPUT_PREFIX)//'_restart_'//trim(timeString)//trim(output_fileSuffix)//'.nc'
  endif

  ! call writeRestart(restartFile,nGRU,nHRU,prog_meta,progStruct,bvar_meta,bvarStruct,maxLayers,maxSnowLayers,indx_meta,indxStruct,err,cmessage)  
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end if

 ! *****************************************************************************
 ! *** update counters
 ! *****************************************************************************

 ! increment output file timestep
 do iFreq = 1,maxvarFreq
  statCounter(iFreq) = statCounter(iFreq)+1
  if(finalizeStats(iFreq)) outputTimeStep(iFreq) = outputTimeStep(iFreq) + 1
 end do

 ! increment forcingStep
 forcingStep=forcingStep+1

 ! if finalized stats, then reset stats on the next time step
 resetStats(:) = finalizeStats(:)

 ! save time vector
 oldTime%var(:) = summa1_struc%timeStruct%var(:)

 ! *****************************************************************************
 ! *** finalize
 ! *****************************************************************************

 ! identify the end of the writing
 call date_and_time(values=endWrite)

 ! aggregate the elapsed time for model writing
 elapsedWrite = elapsedWrite + elapsedSec(startWrite, endWrite)

 ! end associate statements
 end associate summaVars

 end subroutine summa_writeOutputFiles
end module summa_writeOutput


