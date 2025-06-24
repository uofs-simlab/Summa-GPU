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

module summa_restart
! read restart data and reset the model state

! access missing values
USE globalData,only:integerMissing   ! missing integer
USE globalData,only:realMissing      ! missing double precision number

! named variables
USE var_lookup,only:iLookPROG                               ! look-up values for local column model prognostic (state) variables
USE var_lookup,only:iLookDIAG                               ! look-up values for local column model diagnostic variables
USE var_lookup,only:iLookFLUX                               ! look-up values for local column model fluxes
USE var_lookup,only:iLookBVAR                               ! look-up values for basin-average model variables
USE var_lookup,only:iLookDECISIONS                          ! look-up values for model decisions
USE var_lookup,only:iLookINDEX                          ! look-up values for model decisions

! safety: set private unless specified otherwise
implicit none
private
public::summa_readRestart
contains

 ! read restart data and reset the model state
 subroutine summa_readRestart(summa1_struc, err, message)
 ! ---------------------------------------------------------------------------------------
 ! * desired modules
 ! ---------------------------------------------------------------------------------------
 ! data types
 USE nrtype                                                  ! variable types, etc.
 USE summa_type, only:summa1_type_dec                        ! master summa data type
 ! functions and subroutines
 USE time_utils_module,only:elapsedSec                       ! calculate the elapsed time
 USE read_icond_module,only:read_icond                       ! module to read initial conditions
 USE check_icond_module,only:check_icond                     ! module to check initial conditions
 USE var_derive_module,only:calcHeight_d                       ! module to calculate height at layer interfaces and layer mid-point
 USE var_derive_module,only:v_shortcut                       ! module to calculate "short-cut" variables
 USE var_derive_module,only:rootDensty                       ! module to calculate the vertical distribution of roots
 USE var_derive_module,only:satHydCond                       ! module to calculate the saturated hydraulic conductivity in each soil layer
 ! global data structures
 USE globalData,only:gru_struc                               ! gru-hru mapping structures
 USE globalData,only:model_decisions                         ! model decision structure
 ! file paths
 USE summaFileManager,only:SETTINGS_PATH                     ! path to settings files (e.g., Noah vegetation tables)
 USE summaFileManager,only:STATE_PATH                        ! optional path to state/init. condition files (defaults to SETTINGS_PATH)
 USE summaFileManager,only:MODEL_INITCOND                    ! name of model initial conditions file
 ! timing variables
 USE globalData,only:startRestart,endRestart                 ! date/time for the start and end of reading model restart files
 USE globalData,only:elapsedRestart                          ! elapsed time to read model restart files
 ! model decisions
 USE mDecisions_module,only:&                                ! look-up values for the choice of method for the spatial representation of groundwater
   localColumn,    & ! separate groundwater representation in each local soil column
   singleBasin       ! single groundwater store over the entire basin
 ! look-up values for the choice of variable in energy equations (BE residual or IDA state variable)
 USE mDecisions_module,only:&
   closedForm,     & ! use temperature with closed form heat capacity
   enthalpyFormLU, & ! use enthalpy with soil temperature-enthalpy lookup tables
   enthalpyForm      ! use enthalpy with soil temperature-enthalpy analytical solution
! look-up values for the choice of full or empty aquifer at start
 USE mDecisions_module,only:&
   fullStart,      & ! start with full aquifer
   emptyStart        ! start with empty aquifer
   use initialize_device
 ! ---------------------------------------------------------------------------------------
 ! * variables
 ! ---------------------------------------------------------------------------------------
 implicit none
 ! dummy variables
 type(summa1_type_dec),intent(inout)   :: summa1_struc       ! master summa data structure
 integer(i4b),intent(out)              :: err                ! error code
 character(*),intent(out)              :: message            ! error message
 ! local variables
 character(LEN=256)                    :: cmessage           ! error message of downwind routine
 character(LEN=256)                    :: restartFile        ! restart file name
 integer(i4b)                          :: iGRU,iHRU          ! looping variables
 logical(lgt)                          :: checkEnthalpy      ! flag if checking enthalpy for consistency
 logical(lgt)                          :: no_icond_enth      ! flag that enthalpy not in initial conditions
 logical(lgt)                          :: use_lookup         ! flag to use the lookup table for soil enthalpy, otherwise use analytical solution
 real(rkind)                           :: aquifer_start      ! initial aquifer storage
 integer(i4b) :: nSnow, nSoil, nLayers
 ! ---------------------------------------------------------------------------------------
 ! associate to elements in the data structure
 summaVars: associate(& 
  ! model decisions
  ixNrgConserv         => model_decisions(iLookDECISIONS%nrgConserv)%iDecision   ,& !choice of variable in either energy backward Euler residual or IDA state variable
  spatial_gw           => model_decisions(iLookDECISIONS%spatial_gw)%iDecision   ,& !choice of method for the spatial representation of groundwater
  aquiferIni           => model_decisions(iLookDECISIONS%aquiferIni)%iDecision   ,& !choice of full or empty aquifer at start
  ! lookup table data structure
  lookupStruct         => summa1_struc%lookupStruct        , & ! x%gru(:)%hru(:)%z(:)%var(:)%lookup(:) -- lookup tables
  ! primary data structures (variable length vectors)
  indxStruct           => summa1_struc%indxStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
  mparStruct           => summa1_struc%mparStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model parameters
  progStruct           => summa1_struc%progStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
  diagStruct           => summa1_struc%diagStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
  fluxStruct           => summa1_struc%fluxStruct          , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
  ! basin-average structures
  bparStruct           => summa1_struc%bparStruct          , & ! x%gru(:)%var(:)            -- basin-average parameters
  bvarStruct           => summa1_struc%bvarStruct          , & ! x%gru(:)%var(:)%dat        -- basin-average variables
  ! miscellaneous variables
  dt_init              => summa1_struc%dt_init             , & ! used to initialize the length of the sub-step for each HRU
  nGRU                 => summa1_struc%nGRU                , & ! number of grouped response units
  nHRU                 => summa1_struc%nHRU                  & ! number of global hydrologic response units
 ) ! assignment to variables in the data structures
 
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='summa_readRestart/'

 ! identify the start of the writing
 call date_and_time(values=startRestart)

 ! *****************************************************************************
 ! *** read/check initial conditions
 ! *****************************************************************************

 ! define restart file path/name
 if(STATE_PATH == '') then
   restartFile = trim(SETTINGS_PATH)//trim(MODEL_INITCOND)
 else
   restartFile = trim(STATE_PATH)//trim(MODEL_INITCOND)
 endif

 ! read initial conditions
 call read_icond(restartFile,                   & ! intent(in):    name of initial conditions file
                 nGRU,                          & ! intent(in):    number of response units
                 mparStruct,                    & ! intent(in):    model parameters
                 progStruct,                    & ! intent(inout): model prognostic variables
                 bvarStruct,                    & ! intent(inout): model basin (GRU) variables
                 indxStruct,                    & ! intent(inout): model indices
                 no_icond_enth,                 & ! intent(out):   flag that enthalpy not in initial conditions
                 err,cmessage)                    ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

! check initial conditions
 checkEnthalpy = .false.
 use_lookup    = .false.
 if(ixNrgConserv .ne. closedForm) checkEnthalpy = .true. ! check enthalpy either for mixed form energy equation or enthalpy state variable
 if(ixNrgConserv==enthalpyFormLU) use_lookup = .true.    ! use lookup tables for soil temperature-enthalpy instead of analytical solution
 call check_icond(nGRU,                         & ! intent(in):    number of response units
                  progStruct,                   & ! intent(inout): model prognostic variables
                  diagStruct,                   & ! intent(inout): model diagnostic variables
                  mparStruct,                   & ! intent(in):    model parameters
                  indxStruct,                   & ! intent(in):    layer indexes
                  lookupStruct,                 & ! intent(in):    lookup tables
                  checkEnthalpy,                & ! intent(in):    flag if need to start with consistent enthalpy
                  no_icond_enth,                & ! intent(in):    flag that enthalpy not in initial conditions
                  use_lookup,                   & ! intent(in):    flag to use the lookup table for soil enthalpy
                  err,cmessage)                   ! intent(out):   error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    nSnow = summa1_struc%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
        nSoil = summa1_struc%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nSoil)%dat(1)
            nLayers = summa1_struc%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nLayers)%dat(1)

    call allocate_device_indx_data(summa1_struc%indxStruct_d,summa1_struc%indxStruct,summa1_struc%nGRU)

  call allocate_device_param_data(summa1_struc%mparStruct_d,summa1_struc%mparStruct%gru(1)%hru(1))
  call allocate_device_attr_data(summa1_struc%attrStruct_d,summa1_struc%attrStruct,summa1_struc%nGRU)
      call allocate_device_type_data(summa1_struc%typeStruct_d,summa1_struc%typeStruct,summa1_struc%nGRU)
 summa1_struc%greenVegFrac_monthly_d = summa1_struc%greenVegFrac_monthly
call allocate_device_diag_data(summa1_struc%diagStruct_d,summa1_struc%diagStruct,nSnow,summa1_struc%nGRU,nLayers,nSoil)
  call allocate_device_forc_data(summa1_struc%forcStruct_d,summa1_struc%forcStruct,summa1_struc%nGRU)
call allocate_device_prog_data(summa1_struc%progStruct_d,summa1_struc%progStruct,summa1_struc%nGRU,nLayers,nSoil)
   call allocate_veg_param_tables(summa1_struc%tables)
  call allocate_veg_parameters(summa1_struc%veg_param,summa1_struc%nGRU)
 call allocate_device_decisions(summa1_struc%decisions)
call allocate_device_flux_data(summa1_struc%fluxStruct_d,summa1_struc%fluxStruct,nSnow,nSoil,summa1_struc%nGRU,nLayers)
  call allocate_device_bvar_data(summa1_struc%bvarStruct_d,summa1_struc%bvarStruct,summa1_struc%nGRU)

  ! *****************************************************************************
  ! *** compute ancillary variables
  ! *****************************************************************************


   ! re-calculate height of each layer
   call calcHeight_d(summa1_struc%nGRU,summa1_struc%indxStruct_d,   & ! layer type
                   summa1_struc%progStruct_d,   & ! model prognostic (state) variables for a local HRU
                   err,cmessage)                       ! error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! calculate vertical distribution of root density
   call rootDensty(summa1_struc%decisions,&
                   summa1_struc%mparStruct_d,   & ! vector of model parameters
                   summa1_struc%indxStruct_d,   & ! data structure of model indices
                   summa1_struc%progStruct_d,   & ! data structure of model prognostic (state) variables
                   summa1_struc%diagStruct_d,   & ! data structure of model diagnostic variables
                   summa1_struc%nGRU, &
                   err,cmessage)                       ! error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


     ! calculate saturated hydraulic conductivity in each soil layer
   call satHydCond(summa1_struc%decisions,&
   summa1_struc%mparStruct_d,   & ! vector of model parameters
                   summa1_struc%indxStruct_d,   & ! data structure of model indices
                   summa1_struc%progStruct_d,   & ! data structure of model prognostic (state) variables
                   summa1_struc%fluxStruct_d,   & ! data structure of model fluxes
                   summa1_struc%nGRU, &
                   err,cmessage)                       ! error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! calculate "short-cut" variables such as volumetric heat capacity
   call v_shortcut(summa1_struc%mparStruct_d,   & ! vector of model parameters
                   summa1_struc%diagStruct_d,   & ! data structure of model diagnostic variables
                   summa1_struc%nGRU, &
                   err,cmessage)                       ! error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! initialize canopy drip
   ! NOTE: canopy drip from the previous time step is used to compute throughfall for the current time step
   summa1_struc%fluxStruct_d%scalarCanopyLiqDrainage = 0._rkind  ! not used

  ! *****************************************************************************
  ! *** initialize time step
  ! *****************************************************************************

    ! *****************************************************************************
  ! *** initialize aquifer storage
  ! *****************************************************************************

  ! initialize aquifer storage
  ! NOTE: this is ugly: need to add capabilities to initialize basin-wide state variables

  ! There are two options for groundwater:
  !  (1) where groundwater is included in the local column (i.e., the HRUs); and
  !  (2) where groundwater is included for the single basin (i.e., the GRUS, where multiple HRUS drain into a GRU).
  ! For water balance calculations it is important to ensure that the local aquifer storage is zero if groundwater is treated as a basin-average state variable (singleBasin);
  !  and ensure that basin-average aquifer storage is zero when groundwater is included in the local columns (localColumn).

  ! select aquifer option
  select case(aquiferIni)
   case(fullStart)
    aquifer_start  = 1._rkind ! Start with full aquifer, since easier to spin up by draining than filling (filling we need to wait for precipitation) 
   case(emptyStart)
    aquifer_start  = 0._rkind ! Start with empty aquifer ! If want to compare model method outputs, empty start leads to quicker equilibrium
   case default
    message=trim(message)//'unable to identify decision for initial aquifer storage'
    return
  end select  ! aquifer option

  ! select groundwater option
  select case(spatial_gw)

   ! the basin-average aquifer storage is not used if the groundwater is included in the local column
   case(localColumn)
    summa1_struc%bvarStruct_d%basin__AquiferStorage = 0._rkind ! set to zero to be clear that there is no basin-average aquifer storage in this configuration
      if(aquiferIni==emptyStart) summa1_struc%progStruct_d%scalarAquiferStorage = aquifer_start ! leave at initialized values if fullStart

   ! the local column aquifer storage is not used if the groundwater is basin-average
   ! (i.e., where multiple HRUs drain to a basin-average aquifer)
   case(singleBasin)
    summa1_struc%bvarStruct_d%basin__AquiferStorage = aquifer_start 
     summa1_struc%progStruct_d%scalarAquiferStorage = 0._rkind  ! set to zero to be clear that there is no local aquifer storage in this configuration

   ! error check
   case default
    message=trim(message)//'unable to identify decision for regional representation of groundwater'
    return

  end select  ! groundwater option

 do iGRU=1,nGRU
  ! initialize time step length for each HRU
  do iHRU=1,gru_struc(iGRU)%hruCount
   dt_init%gru(iGRU)%hru(iHRU) = summa1_struc%progStruct_d%dt_init ! seconds
  end do

 end do  ! end looping through GRUs


  call finalize_device_indx_data(summa1_struc%indxStruct_d,summa1_struc%indxStruct,summa1_struc%nGRU)

         nSnow = summa1_struc%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
        nSoil = summa1_struc%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nSoil)%dat(1)
            nLayers = summa1_struc%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nLayers)%dat(1)

 call finalize_device_forc_data(summa1_struc%forcStruct_d,summa1_struc%forcStruct,summa1_struc%nGRU)
call finalize_device_diag_data(summa1_struc%diagStruct_d,summa1_struc%diagStruct,nSnow,nSoil,nLayers,summa1_struc%nGRU)
call finalize_device_prog_data(summa1_struc%progStruct_d,summa1_struc%progStruct,nLayers,nSoil,summa1_struc%nGRU)
   call finalize_device_flux_data(summa1_struc%fluxStruct_d,summa1_struc%fluxStruct,nSnow,nSoil,summa1_struc%nGRU)
 call finalize_device_bvar_data(summa1_struc%bvarStruct_d,summa1_struc%bvarStruct,summa1_struc%nGRU)


 ! *****************************************************************************
 ! *** finalize
 ! *****************************************************************************

 ! identify the end of the writing
 call date_and_time(values=endRestart)

 ! aggregate the elapsed time for model writing
 elapsedRestart = elapsedSec(startRestart, endRestart)

 ! end associate statements
 end associate summaVars

 end subroutine summa_readRestart
end module summa_restart




