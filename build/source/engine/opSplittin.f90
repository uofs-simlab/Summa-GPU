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

module opSplittin_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing   ! missing integer

! state variable type
USE globalData,only:iname_nrgCanair  ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy  ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy  ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy  ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer   ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer   ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer   ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer   ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer   ! named variable defining the liquid matric potential state variable for soil layers
USE globalData,only:iname_watAquifer ! named variable defining the water storage in the aquifer

! global metadata
USE globalData,only:flux_meta        ! metadata on the model fluxes
USE globalData,only:diag_meta        ! metadata on the model diagnostic variables
USE globalData,only:prog_meta        ! metadata on the model prognostic variables
USE globalData,only:deriv_meta       ! metadata on the model derivatives
USE globalData,only:flux2state_orig_d  ! metadata on flux-to-state mapping (original state variables)
USE globalData,only:flux2state_liq   ! metadata on flux-to-state mapping (liquid water state variables)
  
! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookFLUX        ! named variables for structure elements
USE var_lookup,only:iLookINDEX       ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS   ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM       ! named variables for elements of the parameter structure

! look up structure for variable types
USE var_lookup,only:iLookVarType

! provide access to the number of flux variables
USE var_lookup,only:nFlux=>maxvarFlux ! number of model flux variables

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,                                                     & ! data vector (i4b)
                    var_d,                                                     & ! data vector (rkind)
                    var_flagVec,                                               & ! data vector with variable length dimension (lgt)
                    var_ilength,                                               & ! data vector with variable length dimension (i4b)
                    var_dlength,                                               & ! data vector with variable length dimension (rkind)
                    zLookup,                                                   & ! lookup tables
                    model_options,                                             & ! defines the model decisions
                    out_type_stateFilter,                                      & ! classes for stateFilter objects
                    in_type_indexSplit,out_type_indexSplit,                    & ! classes for indexSplit objects
                    in_type_varSubstep,io_type_varSubstep,out_type_varSubstep    ! classes for varSubstep objects

! look-up values for the numerical method
USE mDecisions_module,only:       &
                    homegrown    ,& ! home-grown backward Euler solution using concepts from numerical recipes
                    kinsol       ,& ! SUNDIALS backward Euler solution using Kinsol
                    ida             ! SUNDIALS solution using IDA

! safety: set private unless specified otherwise
implicit none
private
public::opSplittin

! named variables for the coupling method
integer(i4b),parameter  :: fullyCoupled=1             ! 1st try: fully coupled solution
integer(i4b),parameter  :: stateTypeSplit=2           ! 2nd try: separate solutions for each state type

! named variables for the state variable split
integer(i4b),parameter  :: nrgSplit=1                 ! order in sequence for the energy operation
integer(i4b),parameter  :: massSplit=2                ! order in sequence for the mass operation

! named variables for the domain type split
integer(i4b),parameter  :: vegSplit=1                 ! order in sequence for the vegetation split
integer(i4b),parameter  :: snowSplit=2                ! order in sequence for the snow split
integer(i4b),parameter  :: soilSplit=3                ! order in sequence for the soil split
integer(i4b),parameter  :: aquiferSplit=4             ! order in sequence for the aquifer split

! named variables for the solution method
integer(i4b),parameter  :: vector=1                   ! vector solution method
integer(i4b),parameter  :: scalar=2                   ! scalar solution method
integer(i4b),parameter  :: nSolutions=2               ! number of solution methods

! named variables for the switch between states and domains
integer(i4b),parameter  :: fullDomain=1               ! full domain (veg+snow+soil)
integer(i4b),parameter  :: subDomain=2                ! sub domain (veg, snow, soil, and aquifer separately)

! maximum number of possible splits
integer(i4b),parameter  :: nStateTypes=2              ! number of state types (energy, water)
integer(i4b),parameter  :: nDomains=4                 ! number of domains (vegetation, snow, soil, and aquifer)

! class definitions

type, public :: split_select_type  ! class for selecting operator splitting methods
  ! opSplittin indices (in order)
  integer(i4b)             :: iSplit                      ! iteration counter for split_select_loop
  integer(i4b)             :: ixCoupling
  integer(i4b)             :: iStateTypeSplit
  integer(i4b)             :: ixStateThenDomain           ! 1=state type split; 2=domain split within a given state type 
  integer(i4b)             :: iDomainSplit
  integer(i4b)             :: ixSolution
  integer(i4b)             :: iStateSplit
  ! variables for specifying the split
  integer(i4b)             :: nState                      ! # of state variables
  integer(i4b)             :: nSubset                     ! number of selected state variables for a given split
  type(var_flagVec)        :: fluxMask                    ! integer mask defining model fluxes
  logical(lgt),allocatable :: stateMask(:)                ! mask defining desired state variables
  ! flags for splitting method control
  logical(lgt)             :: stateTypeSplitting,stateThenDomain,domainSplit,solution,stateSplit
 contains
  procedure :: initialize_flags             => split_select_initialize_flags             ! initialize flags that control operations
  procedure :: initialize_ixCoupling        => split_select_initialize_ixCoupling        ! initialize operator splitting indices
  procedure :: initialize_iStateTypeSplit   => split_select_initialize_iStateTypeSplit   ! initialize operator splitting indices
  procedure :: initialize_ixStateThenDomain => split_select_initialize_ixStateThenDomain ! initialize operator splitting indices
  procedure :: initialize_iDomainSplit      => split_select_initialize_iDomainSplit      ! initialize operator splitting indices
  procedure :: initialize_ixSolution        => split_select_initialize_ixSolution        ! initialize operator splitting indices
  procedure :: initialize_iStateSplit       => split_select_initialize_iStateSplit       ! initialize operator splitting indices

  procedure :: get_stateMask                => split_select_compute_stateMask            ! compute stateMask and nSubset and load into class object

  procedure :: advance_iSplit               => split_select_advance_iSplit               ! advance coupling iterator
  procedure :: advance_ixCoupling           => split_select_advance_ixCoupling           ! advance coupling iterator
  procedure :: advance_iStateTypeSplit      => split_select_advance_iStateTypeSplit      ! advance stateTypeSplitting iterator
  procedure :: advance_ixStateThenDomain    => split_select_advance_ixStateThenDomain    ! advance stateThenDomain iterator
  procedure :: advance_iDomainSplit         => split_select_advance_iDomainSplit         ! advance domainSplit iterator
  procedure :: advance_ixSolution           => split_select_advance_ixSolution           ! advance solution iterator
  procedure :: advance_iStateSplit          => split_select_advance_iStateSplit          ! advance stateSplit iterator
  
  procedure :: logic_exit_stateTypeSplitting => split_select_logic_exit_stateTypeSplitting ! get logical for branch
  procedure :: logic_exit_stateThenDomain    => split_select_logic_exit_stateThenDomain    ! get logical for branch
  procedure :: logic_exit_domainSplit        => split_select_logic_exit_domainSplit        ! get logical for branch
  procedure :: logic_exit_solution           => split_select_logic_exit_solution           ! get logical for branch
  procedure :: logic_exit_stateSplit         => split_select_logic_exit_stateSplit         ! get logical for branch

  procedure :: logic_initialize_stateTypeSplitting => split_select_logic_initialize_stateTypeSplitting ! get logical for branch
  procedure :: logic_initialize_stateThenDomain    => split_select_logic_initialize_stateThenDomain    ! get logical for branch
  procedure :: logic_initialize_domainSplit        => split_select_logic_initialize_domainSplit        ! get logical for branch
  procedure :: logic_initialize_solution           => split_select_logic_initialize_solution           ! get logical for branch
  procedure :: logic_initialize_stateSplit         => split_select_logic_initialize_stateSplit         ! get logical for branch

  procedure :: logic_finalize_stateTypeSplitting => split_select_logic_finalize_stateTypeSplitting     ! get logical for branch
  procedure :: logic_finalize_stateThenDomain    => split_select_logic_finalize_stateThenDomain        ! get logical for branch
  procedure :: logic_finalize_domainSplit        => split_select_logic_finalize_domainSplit            ! get logical for branch
  procedure :: logic_finalize_solution           => split_select_logic_finalize_solution               ! get logical for branch
  procedure :: logic_finalize_stateSplit         => split_select_logic_finalize_stateSplit             ! get logical for branch
end type split_select_type

contains

! **********************************************************************************************************
! public subroutine opSplittin: run the coupled energy-mass model for one timestep
!
! The logic of the solver is as follows:
! (1) Attempt different solutions in the following order: (a) fully coupled; (b) split by state type and by
!      domain type for a given energy and mass split (vegetation, snow, and soil); and (c) scalar solution
!      for a given state type and domain subset.
! (2) For a given split, compute a variable number of substeps (in varSubstep).
! **********************************************************************************************************
subroutine opSplittin(&
  nGRU, &
                      ! input: model control
                      nSnow,                & ! intent(in):    number of snow layers
                      nSoil,                & ! intent(in):    number of soil layers
                      nLayers,              & ! intent(in):    total number of layers
                      nState,               & ! intent(in):    total number of state variables
                      dt,                   & ! intent(in):    time step (s)
                      whole_step,           & ! intent(in):    length of whole step for surface drainage and average flux
                      firstSubStep,         & ! intent(in):    flag to denote first sub-step
                      firstInnerStep,       & ! intent(in):    flag to denote if the last time step in maxstep subStep
                      computeVegFlux,       & ! intent(in):    flag to denote if computing energy flux over vegetation
                       ! input/output: data structures
                      type_data,            & ! intent(in):    type of vegetation and soil
                      attr_data,            & ! intent(in):    spatial attributes
                      forc_data,            & ! intent(in):    model forcing data
                      mpar_data,            & ! intent(in):    model parameters
                      indx_data,            & ! intent(inout): index data
                      prog_data,            & ! intent(inout): model prognostic variables for a local HRU
                      diag_data,            & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,            & ! intent(inout): model fluxes for a local HRU
                      bvar_data,            & ! intent(in):    model variables for the local basin
                      lookup_data,          & ! intent(in):    lookup tables
                      model_decisions,      & ! intent(in):    model decisions
                      decisions, &
                      veg_param, &
                      ! output: model control
                      dtMultiplier,         & ! intent(out):   substep multiplier (-)
                      tooMuchMelt,          & ! intent(out):   flag to denote that ice is insufficient to support melt
                      stepFailure,          & ! intent(out):   flag to denote step failure
                      ixSolution,           & ! intent(out):   solution method used in this iteration
                      mean_step_dt,         & ! intent(out):   mean solution step for the time step
                      err,message)            ! intent(out):   error code and error message
  ! ---------------------------------------------------------------------------------------
  ! structure allocations
  USE allocspace_module,only:allocLocal                ! allocate local data structures
  ! population/extraction of state vectors
  USE indexState_module,only:indexSplit                ! get state indices
  USE varSubstep_module,only:varSubstep                ! complete substeps for a given split
  use device_data_types
  use initialize_device,only:allocate_device_deriv_data,deallocate_device_deriv_data
  use initialize_device,only:allocate_device_flux_prev
  use initialize_device,only:allocate_device_flux_count
  use initialize_device,only:allocate_device_diag_temp,deallocate_device_diag_data
  use initialize_device,only:allocate_device_prog_temp,deallocate_device_prog_data
  use initialize_device,only:allocate_device_flux_mask,zero_device_deriv_data
  implicit none
  ! ---------------------------------------------------------------------------------------
  ! * dummy variables
  ! ---------------------------------------------------------------------------------------
  ! input: model control
  integer(i4b),intent(in) :: nGRU
  integer(i4b),intent(in)         :: nSnow                          ! number of snow layers
  integer(i4b),intent(in)         :: nSoil                          ! number of soil layers
  integer(i4b),intent(in)         :: nLayers                        ! total number of layers
  integer(i4b),intent(in)         :: nState                         ! total number of state variables
  real(rkind),intent(in)          :: dt                             ! time step (seconds)
  real(rkind),intent(in)          :: whole_step                     ! length of whole step for surface drainage and average flux
  logical(lgt),intent(in)         :: firstSubStep                   ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(in)         :: firstInnerStep                 ! flag to denote if the first time step in maxstep subStep
  logical(lgt),intent(in),device         :: computeVegFlux(:)                 ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  ! input/output: data structures
  type(type_data_device),intent(in)          :: type_data                      ! type of vegetation and soil
  type(attr_data_device),intent(in)          :: attr_data                      ! spatial attributes
  type(forc_data_device),intent(in)          :: forc_data                      ! model forcing data
  type(mpar_data_device),intent(in)    :: mpar_data                      ! model parameters
  type(indx_data_device),intent(inout) :: indx_data                      ! indices for a local HRU
  type(prog_data_device),intent(inout) :: prog_data                      ! prognostic variables for a local HRU
  type(diag_data_device),intent(inout) :: diag_data                      ! diagnostic variables for a local HRU
  type(flux_data_device),intent(inout) :: flux_data                      ! model fluxes for a local HRU
  type(bvar_data_device),intent(in)    :: bvar_data                      ! model variables for the local basin
  type(zLookup_device),    intent(in)    :: lookup_data                    ! lookup tables
  type(model_options),intent(in)  :: model_decisions(:)             ! model decisions
  type(decisions_device) :: decisions
  type(veg_parameters) :: veg_param
  ! output: model control
  real(rkind),intent(out)         :: dtMultiplier                   ! substep multiplier (-)
  logical(lgt),intent(out)        :: tooMuchMelt                    ! flag to denote that ice is insufficient to support melt
  logical(lgt),intent(out)        :: stepFailure                    ! flag to denote step failure
  integer(i4b),intent(out)        :: ixSolution                     ! index of solution method (1,2)
  real(rkind),intent(out)         :: mean_step_dt                   ! mean solution step for the time step
  integer(i4b),intent(out)        :: err                            ! error code
  character(*),intent(out)        :: message                        ! error message
  ! ---------------------------------------------------------------------------------------
  ! * general local variables
  ! ---------------------------------------------------------------------------------------
  character(LEN=256)              :: cmessage                       ! error message of downwind routine
  integer(i4b)                    :: minLayer                       ! the minimum layer used in assigning flags for flux aggregations
  integer(i4b)                    :: iOffset                        ! offset to account for different indices in the soil domain
  integer(i4b)                    :: iMin(1),iMax(1)                ! bounds of a given vector
  integer(i4b)                    :: iLayer,jLayer                  ! index of model layer
  integer(i4b)                    :: iVar                           ! index of variables in data structures
  logical(lgt)                    :: firstSuccess                   ! flag to define the first success
  logical(lgt)                    :: firstFluxCall                  ! flag to define the first flux call
  logical(lgt)                    :: reduceCoupledStep              ! flag to define the need to reduce the length of the coupled step
  logical(lgt)                    :: return_flag                    ! flag to indicate the execution of a return statement
  type(prog_data_device)               :: prog_temp                      ! temporary model prognostic variables
  type(diag_data_device)               :: diag_temp                      ! temporary model diagnostic variables
  type(flux_data_device)               :: flux_temp                      ! temporary model fluxes
  type(flux_data_device)               :: flux_mean                      ! mean model fluxes
  type(flux_data_device)               :: flux_mntemp                    ! temporary mean model fluxes
  type(deriv_data_device)               :: deriv_data                     ! derivatives in model fluxes w.r.t. relevant state variables
  ! ------------------------------------------------------------------------------------------------------
  ! * operator splitting
  ! ------------------------------------------------------------------------------------------------------
  ! minimum timestep
  real(rkind)                     :: dtmin_coupled                  ! minimum time step for the fully coupled solution (seconds)
  real(rkind)                     :: dtmin_split                    ! minimum time step for the fully split solution (seconds)
  real(rkind)                     :: dtmin_scalar                   ! minimum time step for the scalar solution (seconds)
  real(rkind)                     :: dt_min                         ! minimum time step (seconds)
  real(rkind)                     :: dtInit                         ! initial time step (seconds)
  ! number of substeps taken for a given split
  integer(i4b)                    :: nSubsteps                      ! number of substeps taken for a given split
  ! named variables defining the coupling and solution method
  integer(i4b)                    :: ixCoupling                     ! index of coupling method (1,2)
  integer(i4b)                    :: ixStateThenDomain              ! switch between the state and domain (1,2)
  integer(i4b)                    :: tryDomainSplit                 ! (0,1) - flag to try the domain split
  ! actual number of splits
  integer(i4b)                    :: nStateTypeSplit                ! number of splits for the state type
  integer(i4b)                    :: nDomainSplit                   ! number of splits for the domain
  integer(i4b)                    :: nStateSplit                    ! number of splits for the states within a given domain
  ! indices for the state type and the domain split
  integer(i4b)                    :: iStateTypeSplit                ! index of the state type split
  integer(i4b)                    :: iDomainSplit                   ! index of the domain split
  integer(i4b)                    :: iStateSplit                    ! index of the state split
  ! flux masks
  ! logical(lgt)                    :: neededFlux(nFlux)              ! .true. if flux is needed at all
  logical(lgt)                    :: desiredFlux                    ! .true. if flux is desired for a given split
  type(flux_count_device)               :: fluxCount                      ! number of times each flux is updated (should equal nSubsteps)
  type(flux_mask_device)               :: fluxMask                       ! mask defining model fluxes
  ! state masks
  integer(i4b),dimension(nState)  :: stateCheck                     ! number of times each state variable is updated (should equal 1)
  logical(lgt),dimension(nState)  :: stateMask                      ! mask defining desired state variables
  integer(i4b)                    :: nSubset                        ! number of selected state variables for a given split
  ! flags
  logical(lgt)                    :: failure                        ! flag to denote failure of substepping
  logical(lgt)                    :: doAdjustTemp                   ! flag to adjust temperature after the mass split
  logical(lgt)                    :: failedMinimumStep              ! flag to denote failure of substepping for a given split
  integer(i4b),device                    :: ixSaturation(nGRU)                   ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
  integer(i4b)                    :: nCoupling
  ! mean steps 
  real(rkind)                     :: mean_step_state                ! mean step over the state (with or without domain splits)
  real(rkind)                     :: mean_step_solution             ! mean step for a solution (scalar or vector)
  logical(lgt)                    :: addFirstFlux                   ! flag to add the first flux to the mask
  ! splitting method control variables
  logical(lgt)                    :: exit_split_select,exit_coupling,exit_stateThenDomain,exit_solution
  logical(lgt)                    :: cycle_split_select,cycle_coupling,cycle_stateThenDomain,cycle_domainSplit,cycle_solution
  ! ------------------------ classes for subroutine arguments (classes defined in data_types module) ------------------------
  !      ** intent(in) arguments **         ||       ** intent(inout) arguments **        ||      ** intent(out) arguments **
  type(in_type_indexSplit)  :: in_indexSplit;                                             type(out_type_indexSplit)  :: out_indexSplit;  ! indexSplit arguments
  type(in_type_varSubstep)  :: in_varSubstep;  type(io_type_varSubstep) :: io_varSubstep; type(out_type_varSubstep)  :: out_varSubstep;  ! varSubstep arguments
  ! -------------------------------------------------------------------------------------------------------------------------
  type(split_select_type) :: split_select ! class object for selecting operator splitting methods
  ! -----

  type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)


  ! set up split_select object and prepare for split_select_loop
  call initialize_opSplittin; if (return_flag) return
 
  ! execute split_select_loop
  call update_opSplittin;     if (return_flag) return

  ! validation and error control
  call finalize_opSplittin;   if (return_flag) return  

 contains

  subroutine initialize_opSplittin
   ! *** Initial operations for opSplittin ***
   ! set splitting parameters
   dtmin_coupled = max(1._rkind, mpar_data%maxstep/NINT(mpar_data%be_steps)/10._rkind)
   dtmin_split   = max(1._rkind, mpar_data%maxstep/NINT(mpar_data%be_steps)/40._rkind)
   dtmin_scalar  = max(1._rkind, mpar_data%maxstep/NINT(mpar_data%be_steps)/60._rkind) 
   call initialize_split_select;   if (return_flag) return ! initialize split selector object (split_select)
   call initialize_split_coupling; if (return_flag) return ! prep for first iteration of update_opSplittin 
  end subroutine initialize_opSplittin

  subroutine update_opSplittin
   ! *** Update operations for opSplittin ***
   ! Note: first loop iteration is the fully coupled method, followed by operator splitting methods if required
   split_select_loop: do
     ! initialize variables
     call initialize_split_solve; if (exit_split_select) exit split_select_loop; if (return_flag) return
     ! update the trial solution via fully coupled method or operator splitting
     call update_split_solve;     if (return_flag) return; if (cycle_split_select) cycle split_select_loop
     ! validate the trial solution
     call finalize_split_solve;   if (exit_split_select) exit split_select_loop; if (return_flag) return
   end do split_select_loop
  end subroutine update_opSplittin

  subroutine finalize_opSplittin
   ! *** Final operations for opSplittin ***
   call finalize_split_coupling; if (return_flag) return

   call deallocate_device_deriv_data(deriv_data)
   deallocate(flux_temp%data)
   deallocate(flux_mntemp%data)
   deallocate(flux_mean%data)
   deallocate(fluxCount%data)
   deallocate(fluxMask%data)
   call deallocate_device_diag_data(diag_temp)
   call deallocate_device_prog_data(prog_temp)
  end subroutine finalize_opSplittin

  subroutine initialize_split_solve
   ! *** Initial operations for solving the selected split ***
   call initialize_split_stateTypeSplitting; if (exit_split_select.or.return_flag) return
   cycle_split_select=.false. ! initialize flag for cycle control of split_select_loop
  end subroutine initialize_split_solve

  subroutine update_split_solve
   ! *** Update operations for solving the selected split ***
   ! Note: split_select data components are used to select the splitting method
   if (split_select % stateTypeSplitting) then ! stateTypeSplitting method begins
     call initialize_split_stateThenDomain
     if (split_select % stateThenDomain) then  ! stateThenDomain method begins
       call initialize_split_domainSplit;    if (return_flag) return
       if (split_select % domainSplit) then    ! domainSplit method begins
         call initialize_split_solution
         if (split_select % solution) then     ! solution method begins
           call initialize_split_stateSplit; if (return_flag) return 
           if (split_select % stateSplit) then ! stateSplit method begins
             ! define masks for selected splitting method
             call initialize_split; if (return_flag) return; if (cycle_initialize_split()) then; cycle_split_select=.true.; return; end if
             ! update trial solution for selected splitting method
             call update_split;     if (return_flag) return
             ! validate trial solution
             call finalize_split;   if (return_flag) return; if (cycle_finalize_split())   then; cycle_split_select=.true.; return; end if
           end if                              ! stateSplit method ends
           call finalize_split_stateSplit
         end if                                ! solution method ends
         call finalize_split_solution
       end if                                  ! domainSplit method ends
       call finalize_split_domainSplit
     end if                                    ! stateThenDomain method ends
     call finalize_split_stateThenDomain;    if (return_flag) return
   end if                                      ! stateTypeSplitting method ends
  end subroutine update_split_solve

  subroutine finalize_split_solve
   ! *** Final operations for solving the selected split ***
   integer(i4b),parameter          :: maxSplit=500       ! >= max number of splitting methods (controls upper limit of split_select loop)               
   call finalize_split_stateTypeSplitting; if (exit_split_select.or.return_flag) return
   if (split_select % iSplit.ge.maxSplit) then ! check for errors - execute fail-safe if needed
    err=20; message=trim(message)//'split_select loop exceeded max number of iterations'; return_flag=.true.; return 
   end if
  end subroutine finalize_split_solve

  subroutine initialize_split
   ! *** Initialize logical masks for selected splitting method ***
   call update_stateMask;  if (return_flag)              return ! get the mask for the state subset - return for a non-zero error code
   call validate_split;    if (return_flag)              return ! verify that the split is valid
                           if (cycle_initialize_split()) return ! if needed, proceed to next split
   
   call save_recover                                            ! save/recover copies of variables and fluxes

   call get_split_indices; if (return_flag)              return ! get indices for a given split - return for a non-zero error code
   call update_fluxMask;   if (return_flag)              return ! define the mask for the fluxes used - return for a non-zero error code
  end subroutine initialize_split

  subroutine update_split
   ! *** Update solution for selected splitting method ***
   call solve_subset; if (return_flag) return ! solve variable subset for one time step - return for a positive error code
  end subroutine update_split

  subroutine finalize_split
   ! *** Finalize solution for selected splitting method ***
   call assess_solution;            if (return_flag)            return ! is solution a success or failure? - return for a recovering solution

   call try_other_solution_methods; if (cycle_finalize_split()) return ! if solution failed to converge, try other splitting methods 

   call confirm_variable_updates;   if (return_flag)            return ! check that state variables are updated - return if error 

   call success_check                                                  ! check for success
   call check_exit_stateThenDomain                                     ! check exit criterion for stateThenDomain split
   call check_exit_solution;        if (return_flag)            return ! check exit criterion for solution split - return if error 
  end subroutine finalize_split

  function cycle_initialize_split() result(cycle_flag)
   ! *** Compute loop cycle flag for initialize_split ***
   logical(lgt) :: cycle_flag
   cycle_flag=any([cycle_domainSplit,cycle_solution])
  end function cycle_initialize_split

  function cycle_finalize_split() result(cycle_flag)
   ! *** Compute loop cycle flag for finalize_split ***
   logical(lgt) :: cycle_flag
   cycle_flag=any([cycle_coupling,cycle_stateThenDomain,cycle_solution]) 
  end function cycle_finalize_split

  subroutine initialize_split_select
   ! *** Initialize split_select class object ***
   
   ! initizlaize iteration counter for split_select_loop
   split_select % iSplit = 1 
 
   ! initialize # of state variables
   split_select % nState = nState

   ! allocate data components
   allocate(split_select % stateMask(1:nState),STAT=err) ! allocate split_select components

   ! check for allocation errors
   if (err/=0) then
    message=trim(message)//'allocation error in initialize_split_select routine for split_select % stateMask'
    return_flag=.true.; return
   else
    return_flag=.false.
   end if
   
   ! initialize split_select % stateMask to default initial case
   split_select % stateMask(1:nState) = .true.
 
   ! initialize flags 
   exit_split_select=.false.
   cycle_split_select=.false.
   call split_select % initialize_flags ! initialize control flags 
  end subroutine initialize_split_select

  subroutine initialize_split_coupling
   ! *** Initialize coupling split method ***
   call split_select % initialize_ixCoupling
   call initialize_coupling; if (return_flag) return ! select coupling options and allocate memory - return if error occurs
  end subroutine initialize_split_coupling

  subroutine initialize_split_stateTypeSplitting
   ! *** Initialize stateTypeSplitting split method ***
   if (split_select % logic_initialize_stateTypeSplitting()) then
    ixCoupling=split_select % ixCoupling 
    if (ixCoupling.gt.nCoupling) then; exit_split_select=.true.; return; end if ! exit if all splits are exhausted 
    call initialize_stateTypeSplitting; if (return_flag) return ! setup steps for stateTypeSplitting split method - return if error occurs
    call split_select % initialize_iStateTypeSplit; split_select % stateTypeSplitting=.true.
   end if
   if (split_select % logic_exit_stateTypeSplitting()) then 
    iStateTypeSplit=split_select % iStateTypeSplit; if (iStateTypeSplit.gt.nStateTypeSplit) split_select % stateTypeSplitting=.false.
   end if
  end subroutine initialize_split_stateTypeSplitting

  subroutine initialize_split_stateThenDomain
   ! *** Initialize stateThenDomain split method ***
   if (split_select % logic_initialize_stateThenDomain()) then 
    ! first try the state type split, then try the domain split within a given state type
    call initialize_stateThenDomain ! setup steps for stateThenDomain split method -- identify state-specific variables for a given state split
    call split_select % initialize_ixStateThenDomain; split_select % stateThenDomain=.true.
   end if
   if (split_select % logic_exit_stateThenDomain()) then ! stateThenDomain
    ixStateThenDomain=split_select % ixStateThenDomain 
    if (ixStateThenDomain > (1+tryDomainSplit)) then
     ixStateThenDomain=ixStateThenDomain-1; split_select % ixStateThenDomain = ixStateThenDomain ! correct index needed after exit
     split_select % stateThenDomain=.false. ! eqivalent to exiting the stateThenDomain method
    end if
   end if
  end subroutine initialize_split_stateThenDomain

  subroutine initialize_split_domainSplit
   ! *** Initialize domainSplit split method ***
   if (split_select % logic_initialize_domainSplit()) then
    call initialize_domainSplit; if (return_flag) return ! setup steps for domainSplit split method - return if error occurs
    call split_select % initialize_iDomainSplit; split_select % domainSplit=.true.
   end if
   if (split_select % logic_exit_domainSplit()) then
    iDomainSplit=split_select % iDomainSplit
    if (split_select % iDomainSplit > nDomainSplit) split_select % domainSplit=.false.
   end if
  end subroutine initialize_split_domainSplit

  subroutine initialize_split_solution
   ! *** Initialize solution split method ***
   if (split_select % logic_initialize_solution()) then; call split_select % initialize_ixSolution; split_select % solution=.true.; end if
   if (split_select % logic_exit_solution()) then
    ixSolution=split_select % ixSolution
    if (split_select % ixSolution > nsolutions) split_select % solution=.false.            
   end if
  end subroutine initialize_split_solution

  subroutine initialize_split_stateSplit
   ! *** Initialize stateSplit split method ***
   if (split_select % logic_initialize_stateSplit()) then
    call initialize_stateSplit; if (return_flag) return ! setup steps for stateSplit split method - return if error occurs
    call split_select % initialize_iStateSplit; split_select % stateSplit=.true.; ! loop through layers (NOTE: nStateSplit=1 for the vector solution, hence no looping)
   end if
   if (split_select % logic_exit_stateSplit()) then ! stateSplit begins
    iStateSplit=split_select % iStateSplit
    if (split_select % iStateSplit > nStateSplit) split_select % stateSplit=.false.; !exit stateSplit
   end if
  end subroutine initialize_split_stateSplit

  subroutine check_exit_stateThenDomain
   ! *** check exit criterion for stateThenDomain split ***
   if (exit_stateThenDomain) then ! exit stateThenDomain split if necessary -- deactivate flags for inner splits 
    call split_select % initialize_ixStateThenDomain 
    split_select % stateThenDomain=.false.; split_select % domainSplit=.false.; split_select % solution=.false.; split_select % stateSplit=.false. 
   end if 
  end subroutine check_exit_stateThenDomain

  subroutine check_exit_solution
   ! *** Check exit criterion for solution split - return if needed ***
   if (split_select % stateThenDomain) then
    if (exit_solution) then; split_select % solution=.false.; split_select % stateSplit=.false.; end if
    if (split_select % solution) then
     if (return_flag) return             ! return if error 
     call split_select % advance_iStateSplit
    end if
   end if
  end subroutine check_exit_solution

  subroutine finalize_split_stateSplit
   ! *** Finalize steps for stateSplit split method ***
   if (split_select % logic_finalize_stateSplit()) then
    call split_select % advance_ixSolution
   end if
  end subroutine finalize_split_stateSplit
           
  subroutine finalize_split_solution
   ! *** Finalize steps for solution split method ***
   if (split_select % logic_finalize_solution()) then
    call finalize_solution ! final steps following solution split method
    call split_select % advance_iDomainSplit
   end if
  end subroutine finalize_split_solution

  subroutine finalize_split_domainSplit
   ! *** Finalize steps for domainSplit split method ***
   if (split_select % logic_finalize_domainSplit()) then
    call split_select % advance_ixStateThenDomain
   end if
  end subroutine finalize_split_domainSplit

  subroutine finalize_split_stateThenDomain 
   ! *** Finalize steps for stateThenDomain split method ***
   if (split_select % logic_finalize_stateThenDomain()) then
    call finalize_stateThenDomain; if (return_flag) return ! final steps following the stateThenDomain split method
    call split_select % advance_iStateTypeSplit
   end if
  end subroutine finalize_split_stateThenDomain 

  subroutine finalize_split_stateTypeSplitting
   ! *** Finalize steps for stateTypeSplitting split method ***
   if (split_select % logic_finalize_stateTypeSplitting()) then
    call finalize_stateTypeSplitting 
    if (exit_coupling) then
     call split_select % initialize_ixCoupling; exit_split_select=.true.; return ! success = exit the coupling split method (split_select_loop)
    end if
    call split_select % advance_ixCoupling
   end if
   call split_select % advance_iSplit ! advance iteration counter for split_select_loop
  end subroutine finalize_split_stateTypeSplitting

  subroutine finalize_split_coupling
   ! *** Finalize steps for coupling split method ***
   call finalize_coupling; if (return_flag) return ! check variables and fluxes, and apply step halving if needed
  end subroutine finalize_split_coupling

  subroutine initialize_coupling
   ! *** initial steps for coupling split method ***
   ! initialize error control
   err=0; message="opSplittin/"

   call get_nCoupling; if (return_flag) return ! get nCoupling value -- return if error

   ! initialize the first success call
   firstSuccess=.false.
   if (.not.firstInnerStep) firstSuccess=.true.

   ! initialize the flags
   tooMuchMelt=.false.  ! too much melt (merge snow layers)
   stepFailure=.false.  ! step failure

   ! initialize flag for the success of the substepping
   failure=.false.

   ! initialize the flux check
  !  neededFlux(:) = .false.

   ! initialize the state check
   stateCheck(:) = 0

   ! allocate local structures based on the number of snow and soil layers
   call allocate_memory
   if (return_flag) return ! return if an error occurs during memory allocation 

   ! intialize the flux counter
  !  do iVar=1,size(flux_meta)  ! loop through fluxes
  !    fluxCount%var(iVar)%dat(:) = 0
  !  end do
   fluxCount%data = 0

   ! initialize the model fluxes
   call initialize_coupling_kernel<<<blocks,threads>>>(nGRU,flux2state_orig_d%state1,flux2state_orig_d%state2,&
   flux_data%data,flux_mean%data,firstInnerStep,computeVegFlux)

   ! initialize derivatives
   call zero_device_deriv_data(deriv_data)
  end subroutine initialize_coupling

  subroutine get_nCoupling
   ! *** Get nCoupling value ***
    ! we just solve the fully coupled problem for IDA
    nCoupling = 1
  end subroutine get_nCoupling

  subroutine allocate_memory
   ! *** allocate memory for local structures ***
   return_flag=.false. ! initialize flag

   ! allocate space for the flux mask (used to define when fluxes are updated)
   call allocate_device_flux_mask(fluxMask,nSoil,nGRU)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the flux count (used to check that fluxes are only updated once)
   call allocate_device_flux_count(fluxCount,nSoil,nGRU)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the temporary prognostic variable structure
   call allocate_device_prog_temp(prog_temp,nGRU,nSoil)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the temporary diagnostic variable structure
   call allocate_device_diag_temp(diag_temp,nGRU,nSoil)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the temporary flux variable structure
   call allocate_device_flux_prev(flux_temp,nSoil,nGRU)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the mean flux variable structure
    call allocate_device_flux_prev(flux_mean,nSoil,nGRU)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the temporary mean flux variable structure
   call allocate_device_flux_prev(flux_mntemp,nSoil,nGRU)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if

   ! allocate space for the derivative structure
   call allocate_device_deriv_data(deriv_data,nGRU,nSoil,nLayers,nSnow)
   if (err/=0) then; err=20; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
  end subroutine allocate_memory

  subroutine finalize_coupling
   ! *** final operations for coupling split method ***
   ! check that all state variables were updated
   if (any(stateCheck==0)) then
    message=trim(message)//'some state variables were not updated!'
    err=20; return_flag=.true.; return
   endif

   ! use step halving if unable to complete the fully coupled solution in one substep
   if (ixCoupling/=fullyCoupled .or. nSubsteps>1) dtMultiplier=0.5_rkind
  end subroutine finalize_coupling

  subroutine initialize_stateTypeSplitting
   ! *** Initial steps to prepare for iterations of the stateTypeSplit split method ***
   return_flag=.false. ! initialize flag
   ! initialize the time step
   dtInit = min(merge(dt,            dtmin_coupled, ixCoupling==fullyCoupled), dt) ! initial time step
   dt_min = min(merge(dtmin_coupled, dtmin_split,   ixCoupling==fullyCoupled), dt) ! minimum time step

   ! get nStateTypeSplit and tryDomainSplit values
   call get_nStateTypeSplit_tryDomainSplit(ixCoupling); if (return_flag) return

   mean_step_dt = 0._rkind ! initialize mean step for the time step
   addFirstFlux = .true.     ! flag to add the first flux to the mask
  end subroutine initialize_stateTypeSplitting

  subroutine get_nStateTypeSplit_tryDomainSplit(ixCoupling_value)
   ! *** Get nStateTypeSplit and tryDomainSplit values ***
   integer(i4b),intent(in) :: ixCoupling_value
   ! keep track of the number of state splits
   associate(numberStateSplit => indx_data%numberStateSplit) ! intent(inout): [i4b] number of state splitting solutions
    if (ixCoupling/=fullyCoupled) numberStateSplit = numberStateSplit + 1
   end associate

   ! define the number of operator splits for the state type
   select case(ixCoupling_value)
    case(fullyCoupled); nStateTypeSplit=1
    case(stateTypeSplit); nStateTypeSplit=nStateTypes
    case default; err=20; message=trim(message)//'coupling case not found'; return_flag=.true.; return
   end select  ! operator splitting option

   ! define if we wish to try the domain split
   select case(ixCoupling_value)
    case(fullyCoupled);   tryDomainSplit=0
    case(stateTypeSplit); tryDomainSplit=1
    case default; err=20; message=trim(message)//'coupling case not found'; return_flag=.true.; return
   end select  ! operator splitting option
  end subroutine get_nStateTypeSplit_tryDomainSplit

  subroutine finalize_stateTypeSplitting
   ! *** Final operations subsequent to the stateTypeSplitting split method ***
   exit_coupling=.false. ! initialize flag for control 
   if (ixCoupling==fullyCoupled .and. .not.failure) then; exit_coupling=.true.; return; end if ! success = exit the coupling method in opSplittin
  end subroutine finalize_stateTypeSplitting

  subroutine initialize_stateThenDomain
   ! *** Identify state-specific variables for a given state split ***
   doAdjustTemp = (ixCoupling/=fullyCoupled .and. iStateTypeSplit==massSplit) ! flag to adjust the temperature
  end subroutine initialize_stateThenDomain

  subroutine finalize_stateThenDomain
   ! *** Final steps following the stateThenDomain split method ***
   ! sum the mean steps for the time step over each state type split
   !if (ixStateThenDomain == 2+tryDomainSplit) ixStateThenDomain=1+tryDomainSplit ! correct index value if stateThenDomain method is completed fully 
   select case(ixStateThenDomain) 
     case(fullDomain); mean_step_dt = mean_step_dt + mean_step_solution/nStateTypeSplit
     case(subDomain);  mean_step_dt = mean_step_dt + mean_step_state/nStateTypeSplit
     case default; err=20; message=trim(message)//'ixStateThenDomain case not found'; return_flag=.true.; return
   end select
  end subroutine finalize_stateThenDomain

  subroutine initialize_domainSplit
   ! *** initial operations to set up domainSplit split method ***
   return_flag=.false. ! initialize flag
   associate(numberDomainSplitNrg => indx_data%numberDomainSplitNrg ,& ! intent(inout): [i4b] number of domain splitting solutions for energy (-)
             numberDomainSplitMass => indx_data%numberDomainSplitMass )! intent(inout): [i4b] number of domain splitting solutions for mass   (-)
    ! keep track of the number of domain splits
    if (iStateTypeSplit==nrgSplit  .and. ixStateThenDomain==subDomain) numberDomainSplitNrg  = numberDomainSplitNrg  + 1
    if (iStateTypeSplit==massSplit .and. ixStateThenDomain==subDomain) numberDomainSplitMass = numberDomainSplitMass + 1
   end associate

   call get_nDomainSplit(ixStateThenDomain); if (return_flag) return ! get nDomainSplit value -- return if error occurs

   ! check that we haven't split the domain when we are fully coupled
   if (ixCoupling==fullyCoupled .and. nDomainSplit==nDomains) then
     message=trim(message)//'cannot split domains when fully coupled'
     return_flag=.true. ! return statement required in opSplittin
     err=20; return
   end if

   mean_step_state = 0._rkind ! initialize mean step for state
  end subroutine initialize_domainSplit

  subroutine get_nDomainSplit(ixStateThenDomain_value)
   ! *** Get nDomainSplit value ***
   integer(i4b),intent(in) :: ixStateThenDomain_value
   ! define the number of domain splits for the state type
   select case(ixStateThenDomain_value)
     case(fullDomain); nDomainSplit=1
     case(subDomain);  nDomainSplit=nDomains
     case default; err=20; message=trim(message)//'coupling case not found';
      return_flag=.true. ! return statement required in opSplittin
      return
   end select
  end subroutine get_nDomainSplit

  subroutine finalize_solution
   ! *** final operations following solution split method ***
   ! sum the mean steps for the state over each domain split
   mean_step_state = mean_step_state + mean_step_solution/nDomainSplit
  end subroutine finalize_solution

  subroutine initialize_stateSplit
   ! *** initial operations to set up stateSplit split method ***
   return_flag=.false. ! initialize flag
   mean_step_solution = 0._rkind ! initialize mean step for a solution

   ! initialize error control
   err=0; message="opSplittin/"

   ! refine the time step
   if (ixSolution==scalar) then
    dtInit = min(dtmin_split, dt)    ! initial time step
    dt_min = min(dtmin_scalar, dt)   ! minimum time step
   end if

   ! initialize the first flux call
   firstFluxCall=.true.
   if (.not.firstInnerStep) firstFluxCall=.false.

   call get_nStateSplit(ixSolution); if (return_flag) return ! get nStateSplit value -- return if error occurs
  end subroutine initialize_stateSplit

  subroutine get_nStateSplit(ixSolution_value)
   ! *** Get nStateSplit value ***
   integer(i4b),intent(in) :: ixSolution_value
   ! get the number of split layers
   select case(ixSolution_value)
    case(vector); nStateSplit=1
    case(scalar); nStateSplit=count(stateMask)
    case default; err=20; message=trim(message)//'unknown solution method'; 
     return_flag=.true. ! return statement required in opSplittin
     return
   end select
  end subroutine get_nStateSplit

  ! **** indexSplit ****
  subroutine initialize_indexSplit
   call in_indexSplit % initialize(nSubset)
  end subroutine initialize_indexSplit

  subroutine finalize_indexSplit
   call out_indexSplit % finalize(err,cmessage)
  end subroutine finalize_indexSplit
  ! **** end indexSplit ****

  ! **** varSubstep ****
  subroutine initialize_varSubstep
   call in_varSubstep % initialize(dt,dtInit,dt_min,whole_step,nSubset,doAdjustTemp,firstSubStep,ixSolution,scalar,iStateSplit)
   call io_varSubstep % initialize(firstFluxCall)
  end subroutine initialize_varSubstep

  subroutine finalize_varSubstep
   call io_varSubstep  % finalize(firstFluxCall)
   call out_varSubstep % finalize(dtMultiplier,nSubsteps,failedMinimumStep,reduceCoupledStep,tooMuchMelt,err,cmessage)
  end subroutine finalize_varSubstep

  subroutine solve_subset 

   ! *** Solve variable subset for one time step ***
   return_flag=.false. ! initialize flag
   ! keep track of the number of scalar solutions
   associate(numberScalarSolutions => indx_data%numberScalarSolutions) ! intent(inout): [i4b] number of scalar solutions
    if (ixSolution==scalar) numberScalarSolutions = numberScalarSolutions + 1
   end associate

   ! solve variable subset for one full time step
   call initialize_varSubstep

   call varSubstep(nGRU,in_varSubstep,computeVegFlux,fluxMask,io_varSubstep,ixSaturation,fluxCount,&                                            ! intent(inout): class objects for model control
                   model_decisions,decisions,veg_param,lookup_data,type_data,attr_data,forc_data,mpar_data,&    ! intent(inout): data structures for model properties
                   indx_data,prog_data,diag_data,flux_data,flux_mean,deriv_data,bvar_data,&
                   out_varSubstep)                                                          ! intent(out): class object for model control

   call finalize_varSubstep
   if (err/=0) then 
    message=trim(message)//trim(cmessage) 
    if (err>0) then ! return for positive error codes
     return_flag=.true.; return
    end if 
   end if ! error control

  end subroutine solve_subset 

  subroutine assess_solution
   ! *** determine whether solution is a success or a failure ***
   return_flag=.false. ! initialize flag

   ! reduce coupled step if failed the minimum step for the scalar solution
   if (failedMinimumStep .and. ixSolution==scalar) reduceCoupledStep=.true.

   ! if too much melt (or some other need to reduce the coupled step) then return
   ! NOTE: need to go all the way back to coupled_em and merge snow layers, as all splitting operations need to occur with the same layer geometry
   if (tooMuchMelt .or. reduceCoupledStep) then
     stepFailure=.true.
     err=0 ! recovering
     return_flag=.true. ! return statement required in opSplittin
     return
   end if

   ! define failure
   failure = (failedMinimumStep .or. err<0)
   if (.not.failure) firstSuccess=.true.

  end subroutine assess_solution

  subroutine try_other_solution_methods 
   ! *** if solution failed to converge, try other splitting methods *** 
   ! initialize flags
   cycle_coupling=.false.
   cycle_stateThenDomain=.false.
   cycle_solution=.false.

   ! try the fully split solution if failed to converge with a minimum time step in the coupled solution
   if (ixCoupling==fullyCoupled .and. failure) then
    call split_select % advance_ixCoupling; call split_select % initialize_flags; ! prep for next iteration
    cycle_coupling=.true.; return; ! return required to execute cycle statement in opSplittin
   end if

   ! try the scalar solution if failed to converge with a minimum time step in the split solution
   if (ixCoupling/=fullyCoupled) then
     select case(ixStateThenDomain)
       case(fullDomain)
        if (failure) then
         call split_select % advance_ixStateThenDomain ! prep for next iteration
         split_select % domainSplit=.false.; split_select % solution=.false.; split_select % stateSplit=.false.; 
         cycle_stateThenDomain=.true.; return ! return required to execute cycle statement in opSplittin
        end if
       case(subDomain)
        if (failure) then
         call split_select % advance_ixSolution; split_select % stateSplit=.false.; ! prep for next iteration 
         cycle_solution=.true.; return ! return required to execute cycle statement in opSplittin
        end if
       case default; err=20; message=trim(message)//'unknown ixStateThenDomain case'
     end select
   end if
  end subroutine try_other_solution_methods 

  subroutine update_stateMask
   ! *** Get the mask for the state subset ***
   call split_select % get_stateMask(indx_data,err,cmessage,message,return_flag)
   nSubset = split_select % nSubset; stateMask = split_select % stateMask 
   if (return_flag) return
  end subroutine update_stateMask

  subroutine validate_split 
   ! *** Verify that the split is valid ***
   ! initialize flags
   cycle_domainSplit=.false.
   cycle_solution=.false.
   return_flag=.false.
   nSubset = split_select % nSubset

   ! check that state variables exist
   if (nSubset==0) then
    call split_select % advance_iDomainSplit 
    split_select % solution=.false.; split_select % stateSplit=.false. 
    cycle_domainSplit=.true. 
    return 
   end if

   ! avoid redundant case where vector solution is of length 1
   if (ixSolution==vector .and. count(stateMask)==1) then
    call split_select % advance_ixSolution; 
    split_select % stateSplit=.false.; 
    cycle_solution=.true. 
    return 
   end if

   ! check that we do not attempt the scalar solution for the fully coupled case
   if (ixCoupling==fullyCoupled .and. ixSolution==scalar) then
     message=trim(message)//'only apply the scalar solution to the fully split coupling strategy'
     err=20; return_flag=.true.; return
   end if

   ! reset the flag for the first flux call
   if (.not.firstSuccess) firstFluxCall=.true.
  end subroutine validate_split 

  subroutine save_recover
    real(rkind),allocatable :: scalarCanopyTemp(:)
   ! save/recover copies of prognostic variables
     scalarCanopyTemp = prog_data%scalarCanopyTemp
     print*, 'scT', scalarCanopyTemp

     select case(failure)
       case(.false.); prog_temp = prog_data
       case(.true.);  prog_data = prog_temp
     end select
     scalarCanopyTemp = prog_data%scalarCanopyTemp
     print*, scalarCanopyTemp
     scalarCanopyTemp = prog_temp%scalarCanopyTemp


   ! save/recover copies of diagnostic variables
     select case(failure)
       case(.false.); diag_temp = diag_data
       case(.true.);  diag_data = diag_temp
     end select

   ! save/recover copies of model fluxes and mean fluxes
     select case(failure)
       case(.false.)
         flux_temp   = flux_data
         flux_mntemp = flux_mean
         addFirstFlux = .false.
       case(.true.)
         flux_data   = flux_temp
         flux_mean   = flux_mntemp
         if (addFirstFlux) addFirstFlux = .true.
     end select
  end subroutine save_recover

  subroutine get_split_indices
   ! *** Get indices for a given split ***
   return_flag=.false. ! initialize flag
   call initialize_indexSplit
   call indexSplit(nGRU,in_indexSplit,stateMask,indx_data,out_indexSplit)
   call finalize_indexSplit
   if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
  end subroutine get_split_indices

  subroutine confirm_variable_updates
   ! *** check that state variables updated ***
   return_flag=.false. ! set flag
   ! check that state variables updated
   where(stateMask) stateCheck = stateCheck+1
   if (any(stateCheck>1)) then
     message=trim(message)//'state variable updated more than once!'
     err=20; return_flag=.true.; return
   end if
  end subroutine confirm_variable_updates

  subroutine success_check
   ! initialize flags
   return_flag=.false.
   exit_stateThenDomain=.false.
   exit_solution=.false.
   ! success = exit solution
   if (.not.failure) then
     ! sum the mean steps for the successful solution type
     mean_step_solution = mean_step_solution + (dt/nSubsteps)/nStateSplit
     select case(ixStateThenDomain)
       case(fullDomain); if (iStateSplit==nStateSplit) exit_stateThenDomain=.true. ! exit stateThenDomain
       case(subDomain);  if (iStateSplit==nStateSplit) exit_solution=.true. ! exit solution
       case default; err=20; message=trim(message)//'unknown ixStateThenDomain case'
     end select
   else ! failure
     call check_failure; return_flag=.true.; return ! check reason for failure and return
   end if  ! success check
  end subroutine success_check

  subroutine check_failure
   ! *** Analyze reason for failure ***
   if (ixSolution==scalar) then ! check that we did not fail for the scalar solution (last resort)
     message=trim(message)//'failed the minimum step for the scalar solution'
     err=20; return
   else ! check for an unexpected failure
     message=trim(message)//'unexpected failure'
     err=20; return
   end if
  end subroutine check_failure

  subroutine update_fluxMask
   ! *** update the fluxMask data structure ***
   return_flag=.false. ! initialize flag
 
   call update_fluxMask_kernel<<<blocks,threads>>>(nGRU,fluxMask%data,indx_data%ixStateType,flux2state_orig_d%state1,flux2state_orig_d%state2,&
     fluxMask%ixScalarSoilResistance, fluxMask%ixScalarStomResistSunlit, fluxMask%ixScalarStomResistShaded,&
     fluxMask%ixScalarPhotosynthesisSunlit, fluxMask%ixScalarPhotosynthesisShaded, &
     firstFluxCall, addFirstFlux)

  end subroutine update_fluxMask

end subroutine opSplittin

! ****** Class procedures for split_select_type class ******

subroutine split_select_initialize_flags(split_select)
 ! *** Initialize flags for opSplittin split methods ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % stateTypeSplitting=.false. 
 split_select % stateThenDomain=.false.
 split_select % domainSplit=.false.
 split_select % solution=.false.
 split_select % stateSplit=.false.
end subroutine split_select_initialize_flags

subroutine split_select_advance_iSplit(split_select)
 ! *** Advance index for coupling split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iSplit = split_select % iSplit + 1
end subroutine split_select_advance_iSplit

subroutine split_select_advance_ixCoupling(split_select)
 ! *** Advance index for coupling split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % ixCoupling = split_select % ixCoupling + 1
end subroutine split_select_advance_ixCoupling

subroutine split_select_advance_iStateTypeSplit(split_select)
 ! *** Advance index for stateTypeSplit split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iStateTypeSplit = split_select % iStateTypeSplit + 1
end subroutine split_select_advance_iStateTypeSplit

subroutine split_select_advance_ixStateThenDomain(split_select)
 ! *** Advance index for stateThenDomain split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % ixStateThenDomain = split_select % ixStateThenDomain + 1
end subroutine split_select_advance_ixStateThenDomain

subroutine split_select_advance_iDomainSplit(split_select)
 ! *** Advance index for domainSplit split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iDomainSplit = split_select % iDomainSplit + 1
end subroutine split_select_advance_iDomainSplit

subroutine split_select_advance_ixSolution(split_select)
 ! *** Advance index for solution split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % ixSolution = split_select % ixSolution + 1
end subroutine split_select_advance_ixSolution

subroutine split_select_advance_iStateSplit(split_select)
 ! *** Advance index for stateSplit split method ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iStateSplit = split_select % iStateSplit + 1
end subroutine split_select_advance_iStateSplit

subroutine split_select_initialize_ixCoupling(split_select)
 ! *** initialize operator splitting indices for split_select_type class ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % ixCoupling        = 1       
end subroutine split_select_initialize_ixCoupling

subroutine split_select_initialize_iStateTypeSplit(split_select)
 ! *** initialize operator splitting indices for split_select_type class ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iStateTypeSplit        = 1       
end subroutine split_select_initialize_iStateTypeSplit

subroutine split_select_initialize_ixStateThenDomain(split_select)
 ! *** initialize operator splitting indices for split_select_type class ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % ixStateThenDomain        = 1       
end subroutine split_select_initialize_ixStateThenDomain

subroutine split_select_initialize_iDomainSplit(split_select)
 ! *** initialize operator splitting indices for split_select_type class ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iDomainSplit        = 1       
end subroutine split_select_initialize_iDomainSplit

subroutine split_select_initialize_ixSolution(split_select)
 ! *** initialize operator splitting indices for split_select_type class ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % ixSolution        = 1       
end subroutine split_select_initialize_ixSolution

subroutine split_select_initialize_iStateSplit(split_select)
 ! *** initialize operator splitting indices for split_select_type class ***
 class(split_select_type),intent(inout) :: split_select               ! class object for operator splitting selector
 split_select % iStateSplit        = 1       
end subroutine split_select_initialize_iStateSplit

logical(lgt) function split_select_logic_initialize_stateTypeSplitting(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_initialize_stateTypeSplitting=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..false.).and.(split_select % stateTypeSplitting.eqv..false.)
end function split_select_logic_initialize_stateTypeSplitting

logical(lgt) function split_select_logic_exit_stateTypeSplitting(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_exit_stateTypeSplitting=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..false.).and.(split_select % stateTypeSplitting.eqv..true.)
end function split_select_logic_exit_stateTypeSplitting

logical(lgt) function split_select_logic_initialize_stateThenDomain(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_initialize_stateThenDomain=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..false.)
end function split_select_logic_initialize_stateThenDomain

logical(lgt) function split_select_logic_exit_stateThenDomain(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_exit_stateThenDomain=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..true.)
end function split_select_logic_exit_stateThenDomain

logical(lgt) function split_select_logic_initialize_domainSplit(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_initialize_domainSplit=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.)
end function split_select_logic_initialize_domainSplit

logical(lgt) function split_select_logic_exit_domainSplit(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_exit_domainSplit=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..true.)
end function split_select_logic_exit_domainSplit

logical(lgt) function split_select_logic_initialize_solution(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_initialize_solution=(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.)
end function split_select_logic_initialize_solution

logical(lgt) function split_select_logic_exit_solution(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_exit_solution=(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..true.)
end function split_select_logic_exit_solution

logical(lgt) function split_select_logic_initialize_stateSplit(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_initialize_stateSplit=(split_select % stateSplit.eqv..false.)
end function split_select_logic_initialize_stateSplit

logical(lgt) function split_select_logic_exit_stateSplit(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_exit_stateSplit=(split_select % stateSplit.eqv..true.)
end function split_select_logic_exit_stateSplit

logical(lgt) function split_select_logic_finalize_stateSplit(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_finalize_stateSplit=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..true.).and.(split_select % stateThenDomain.eqv..true.)
end function split_select_logic_finalize_stateSplit

logical(lgt) function split_select_logic_finalize_solution(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_finalize_solution=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % stateThenDomain.eqv..true.)
end function split_select_logic_finalize_solution

logical(lgt) function split_select_logic_finalize_domainSplit(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_finalize_domainSplit=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..true.)
end function split_select_logic_finalize_domainSplit

logical(lgt) function split_select_logic_finalize_stateThenDomain(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_finalize_stateThenDomain=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..false.)
end function split_select_logic_finalize_stateThenDomain

logical(lgt) function split_select_logic_finalize_stateTypeSplitting(split_select)
 ! *** Compute logical for branch in split_select loop ***
 class(split_select_type),intent(in)    :: split_select               ! class object for operator splitting selector
 split_select_logic_finalize_stateTypeSplitting=&
 &(split_select % stateSplit.eqv..false.).and.(split_select % solution.eqv..false.).and.(split_select % domainSplit.eqv..false.).and.(split_select % stateThenDomain.eqv..false.).and.(split_select % stateTypeSplitting.eqv..false.)
end function split_select_logic_finalize_stateTypeSplitting

subroutine split_select_compute_stateMask(split_select,indx_data,err,cmessage,message,return_flag)
  use device_data_types
 ! *** Get the mask for the state subset ***
 class(split_select_type),intent(inout) :: split_select              ! class object for operator splitting selector
 type(indx_data_device),intent(in)           :: indx_data                 ! indices for a local HRU
 integer(i4b),intent(out)               :: err                       ! intent(out): error code
 character(*),intent(out)               :: cmessage                  ! intent(out): error message
 character(*),intent(out)               :: message                   ! error message
 logical(lgt),intent(out)               :: return_flag               ! return flag
 ! local variables
 type(out_type_stateFilter)             :: out_stateFilter           ! number of selected state variables for a given split and error control

 err=0               ! initialize error code 
 return_flag=.false. ! initialize flag
 call stateFilter(indx_data,split_select,out_stateFilter)
 call out_stateFilter % finalize(err,cmessage)
 if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if  ! error control
end subroutine split_select_compute_stateMask


! **********************************************************************************************************
! private subroutine stateFilter: get a mask for the desired state variables
! **********************************************************************************************************
 subroutine stateFilter(indx_data,split_select,out_stateFilter)
 USE indexState_module,only:indxSubset                               ! get state indices
 use device_data_types
 implicit none
 ! input
 type(indx_data_device),intent(in)           :: indx_data                 ! indices for a local HRU
 ! input-output
 type(split_select_type),intent(inout)  :: split_select              ! class object for operator splitting selector
 ! output
 type(out_type_stateFilter),intent(out) :: out_stateFilter           ! number of selected state variables for a given split and error control
 ! local
 integer(i4b),allocatable               :: ixSubset(:)               ! list of indices in the state subset
 character(len=256)                     :: cmessage                  ! error message
 logical(lgt)                           :: return_flag               ! flag to indicate a return 
 ! ----------------------------------------------------------------------------------------------------------------------------------------------------
 ! data structures
 associate(ixCoupling => split_select % ixCoupling  ,&  ! intent(in): [i4b] index of coupling method (1,2)
           err        => out_stateFilter % err      ,&  ! intent(out): error code
           message    => out_stateFilter % cmessage  )  ! intent(out): error message
   
  err=0; message='stateFilter/'; return_flag=.false. ! initialize error control

  ! identify splitting option
  select case(ixCoupling)
   ! *** fully coupled ***
   case(fullyCoupled); call fullyCoupled_stateMask ! get stateMask for fully coupled method 
    ! check
   case default; err=20; message=trim(message)//'unable to identify coupling method'; return_flag=.true.; return
  end select  ! selecting solution method 

  ! initialize ixSubset
  allocate(ixSubset(1_i4b),STAT=err)
  if (err/=0) then; message=trim(message)//'allocation error in stateFilter for ixSubset'; return_flag=.true.; return; end if
  ixSubset = 0._rkind
 end associate

 ! get the number of selected state variables
 split_select % nSubset = count(split_select % stateMask)

contains

 subroutine fullyCoupled_stateMask
  ! *** Get fully coupled stateMask ***
  split_select % stateMask(:) = .true. ! use all state variables
 end subroutine fullyCoupled_stateMask
 
end subroutine stateFilter

attributes(global) subroutine update_fluxMask_kernel(nGRU,fluxMask_data,ixStateType_subset,&
  flux2state_orig_state1,flux2state_orig_state2,&
  scalarSoilResistance, scalarStomResistSunlit, scalarStomResistShaded,&
  scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded, &
  firstFluxCall, addFirstFlux)
  integer(i4b),intent(in),value :: nGRU
  logical(lgt),intent(inout) :: fluxMask_data(:,:)
  integer(i4b),intent(in) :: ixStateType_subset(:,:)
  integer(i4b),intent(in) :: flux2state_orig_state1(:), flux2state_orig_state2(:)
  integer(i4b),intent(in),value :: scalarSoilResistance, scalarStomResistSunlit, scalarStomResistShaded
  integer(i4b),intent(in),value :: scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded
  logical(lgt),intent(in),value :: firstFluxCall, addFirstFlux

  integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return
  call update_fluxMask_device(fluxMask_data(:,iGRU),ixStateType_subset(:,iGRU),&
  flux2state_orig_state1,flux2state_orig_state2,&
  scalarSoilResistance, scalarStomResistSunlit, scalarStomResistShaded,&
  scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded, &
  firstFluxCall, addFirstFlux)
end subroutine

attributes(device) subroutine update_fluxMask_device(fluxMask_data,ixStateType_subset,&
  flux2state_orig_state1,flux2state_orig_state2,&
  scalarSoilResistance, scalarStomResistSunlit, scalarStomResistShaded,&
  scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded, &
  firstFluxCall, addFirstFlux)
  implicit none
  logical(lgt),intent(inout) :: fluxMask_data(:)
  integer(i4b),intent(in) :: ixStateType_subset(:)
  integer(i4b),intent(in) :: flux2state_orig_state1(:), flux2state_orig_state2(:)
  integer(i4b),intent(in) :: scalarSoilResistance, scalarStomResistSunlit, scalarStomResistShaded
  integer(i4b),intent(in) :: scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded
  logical(lgt),intent(in) :: firstFluxCall, addFirstFlux

  logical(lgt) :: desiredFlux
  integer(i4b) :: iVar

  do iVar=1,size(fluxMask_data)
    desiredFlux = any(ixStateType_subset==flux2state_orig_state1(iVar)) .or. any(ixStateType_subset==flux2state_orig_state2(iVar))
    ! print*, iVar, fluxMask_data(iVar), desiredFlux
    fluxMask_data(iVar) = desiredFlux
  end do

  ! make sure firstFluxCall fluxes are included in the mask
  if (firstFluxCall .and. addFirstFlux) then 
    fluxMask_data(scalarSoilResistance) = .true.
    fluxMask_data(scalarStomResistSunlit) = .true.
    fluxMask_data(scalarStomResistShaded) = .true.
    fluxMask_data(scalarPhotosynthesisSunlit) = .true.
    fluxMask_data(scalarPhotosynthesisShaded) = .true.
  end if

end subroutine

attributes(global) subroutine initialize_coupling_kernel(nGRU,flux2state_orig_state1, flux2state_orig_state2, flux_data, flux_mean, firstInnerStep, computeVegFlux)
  integer(i4b),intent(in),value :: nGRU
  integer(i4b),intent(in) :: flux2state_orig_state1(:), flux2state_orig_state2(:)
  real(rkind),intent(inout) :: flux_data(:,:), flux_mean(:,:)
  logical(lgt),intent(in) :: computeVegFlux(:)
  logical(lgt),intent(in),value :: firstInnerStep

    integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return
  call initialize_coupling_fluxes(flux2state_orig_state1, flux2state_orig_state2, flux_data(:,iGRU), flux_mean(:,iGRU), firstInnerStep, computeVegFlux(iGRU))
end subroutine

attributes(device) subroutine initialize_coupling_fluxes(flux2state_orig_state1, flux2state_orig_state2, flux_data, flux_mean, firstInnerStep, computeVegFlux)
  integer(i4b),intent(in) :: flux2state_orig_state1(:), flux2state_orig_state2(:)
  real(rkind),intent(inout) :: flux_data(:), flux_mean(:)
  logical(lgt),intent(in) :: firstInnerStep, computeVegFlux

  integer(i4b) :: iVar
     ! initialize the model fluxes
   do iVar=1,size(flux_mean)  ! loop through fluxes
    if (flux2state_orig_state1(iVar)==integerMissing .and. flux2state_orig_state2(iVar)==integerMissing) cycle ! flux does not depend on state (e.g., input)
    if (flux2state_orig_state1(iVar)==iname_watCanopy .and. .not.computeVegFlux) cycle ! use input fluxes in cases where there is no canopy
    if (firstInnerStep) flux_data(iVar) = 0._rkind
    flux_mean(iVar) = 0._rkind
   end do
  end subroutine


end module opSplittin_module
