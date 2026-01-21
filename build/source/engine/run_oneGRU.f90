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

module run_oneGRU_module

! numerical recipes data types
USE nrtype

! access integers to define "yes" and "no"
USE globalData,only:yes,no             ! .true. and .false.

! define data types
USE data_types,only:&
                    ! GRU-to-HRU mapping
                    gru2hru_map,     & ! HRU info
                    ! no spatial dimension
                    var_ilength,     & ! x%var(:)%dat        (i4b)
                    var_dlength,     & ! x%var(:)%dat        (rkind)
                    ! hru dimension
                    hru_int,         & ! x%hru(:)%var(:)     (i4b)
                    hru_int8,        & ! x%hru(:)%var(:)     (i8b)
                    hru_double,      & ! x%hru(:)%var(:)     (rkind)
                    hru_intVec,      & ! x%hru(:)%var(:)%dat (i4b)
                    hru_doubleVec,   & ! x%hru(:)%var(:)%dat (rkind)
                    ! hru+z dimension
                    hru_z_vLookup      ! x%hru(:)%z(:)%var(:)%lookup(:)

! provide access to the named variables that describe elements of parameter structures
USE var_lookup,only:iLookTYPE          ! look-up values for classification of veg, soils etc.
USE var_lookup,only:iLookID            ! look-up values for hru and gru IDs
USE var_lookup,only:iLookATTR          ! look-up values for local attributes
USE var_lookup,only:iLookINDEX         ! look-up values for local column index variables
USE var_lookup,only:iLookFLUX          ! look-up values for local column model fluxes
USE var_lookup,only:iLookBVAR          ! look-up values for basin-average model variables

! provide access to model decisions
USE globalData,only:model_decisions    ! model decision structure
USE var_lookup,only:iLookDECISIONS     ! look-up values for model decisions

! provide access to the named variables that describe model decisions
USE mDecisions_module,only:&           ! look-up values for the choice of method for the spatial representation of groundwater
 localColumn, &                        ! separate groundwater representation in each local soil column
 singleBasin, &                        ! single groundwater store over the entire basin
 bigBucket                             ! a big bucket (lumped aquifer model)
! -----------------------------------------------------------------------------------------------------------------------------------
implicit none
private
public::run_oneGRU

contains

 ! ************************************************************************************************
 ! public subroutine run_oneGRU: simulation for a single GRU
 ! ************************************************************************************************

 ! simulation for a single GRU
 subroutine run_oneGRU(&
    nGRU, &
                       ! model control
                       gruInfo,            & ! intent(inout): HRU information for given GRU (# HRUs, #snow+soil layers)
                       dt_init,            & ! intent(inout): used to initialize the length of the sub-step for each HRU
                       ixComputeVegFlux,   & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                       ! data structures (input)
                       timeVec,            & ! intent(in):    model time data
                       typeHRU,            & ! intent(in):    local classification of soil veg etc. for each HRU
                       idHRU,              & ! intent(in):    local classification of hru and gru IDs
                       attrHRU,            & ! intent(in):    local attributes for each HRU
                       lookupHRU,          & ! intent(in):    local lookup tables for each HRU
                       decisions, veg_param, veg_tables, &
                       ! data structures (input-output)
                       mparHRU,            & ! intent(inout):    local model parameters
                       indxHRU,            & ! intent(inout): model indices
                       forcHRU,            & ! intent(inout): model forcing data
                       progHRU,            & ! intent(inout): prognostic variables for a local HRU
                       diagHRU,            & ! intent(inout): diagnostic variables for a local HRU
                       fluxHRU,            & ! intent(inout): model fluxes for a local HRU
                       bvarData,           & ! intent(inout): basin-average variables
                       ! error control
                       err,message)          ! intent(out):   error control

 ! ----- define downstream subroutines -----------------------------------------------------------------------------------
 USE run_oneHRU_module,only:run_oneHRU                       ! module to run for one HRU
 use initialize_device
 use device_data_types
 ! ----- define dummy variables ------------------------------------------------------------------------------------------
 implicit none
 ! model control
 type(gru2hru_map)   , intent(inout) :: gruInfo(:)              ! HRU information for given GRU (# HRUs, #snow+soil layers)
 real(rkind)            , intent(inout) :: dt_init(:)           ! used to initialize the length of the sub-step for each HRU
 integer(i4b)        , intent(inout),device :: ixComputeVegFlux(:)  ! flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
 ! data structures (input)
 integer(i4b)        , intent(in)    :: timeVec(:)           ! integer vector      -- model time data
 type(type_data_device)       , intent(in)    :: typeHRU              ! x%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
 type(hru_int8)      , intent(in)    :: idHRU                ! x%hru(:)%var(:)     -- local classification of hru and gru IDs
 type(attr_data_device)    , intent(in)    :: attrHRU              ! x%hru(:)%var(:)     -- local attributes for each HRU
 type(zLookup_device) , intent(in)    :: lookupHRU            ! x%hru(:)%z(:)%var(:)%lookup(:) -- lookup values for each HRU
 ! data structures (input-output)
 type(mpar_data_device) , intent(inout) :: mparHRU              ! x%hru(:)%var(:)%dat -- local (HRU) model parameters
 type(indx_data_device)    , intent(inout) :: indxHRU              ! x%hru(:)%var(:)%dat -- model indices
 type(forc_data_device)    , intent(inout) :: forcHRU              ! x%hru(:)%var(:)     -- model forcing data
 type(prog_data_device) , intent(inout) :: progHRU              ! x%hru(:)%var(:)%dat -- model prognostic (state) variables
 type(diag_data_device) , intent(inout) :: diagHRU              ! x%hru(:)%var(:)%dat -- model diagnostic variables
 type(flux_data_device) , intent(inout) :: fluxHRU              ! x%hru(:)%var(:)%dat -- model fluxes
 type(bvar_data_device)   , intent(inout) :: bvarData             ! x%var(:)%dat        -- basin-average variables
 ! error control
 integer(i4b)        , intent(out)   :: err                  ! error code
 character(*)        , intent(out)   :: message              ! error message
 ! ----- define local variables ------------------------------------------------------------------------------------------
 ! general local variables
 character(len=256)                      :: cmessage               ! error message
 integer(i4b)                            :: iHRU                   ! HRU index
 integer(i4b)                            :: jHRU,kHRU              ! index of the hydrologic response unit
 real(rkind)                             :: fracHRU                ! fractional area of a given HRU (-)
 logical(lgt),device                            :: computeVegFluxFlag(nGRU)     ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)

 type(decisions_device) :: decisions
type(veg_parameters) :: veg_param
type(veg_param_tables) :: veg_tables
 integer(i4b) :: nGRU
 integer(i4b) :: iGRU
   type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

 ! initialize error control
 err=0; !write(message, '(A21,I0,A10,I0,A2)' ) 'run_oneGRU (gru nc = ',gruInfo%gru_nc,', gruId = ',gruInfo%gru_id,')/'

 ! ----- basin initialization --------------------------------------------------------------------------------------------


  call initialize_run_oneGRU<<<blocks,threads>>>(nGRU,bvarData%basin__SurfaceRunoff,bvarData%basin__SoilDrainage,&
bvarData%basin__ColumnOutflow,bvarData%basin__TotalRunoff,&
bvarData%basin__AquiferRecharge,bvarData%basin__AquiferBaseflow,bvarData%basin__AquiferTranspire,&
fluxHRU%data,fluxHRU%ixMLayerColumnInflow_start,fluxHRU%ixMLayerColumnInflow_end,&
computeVegFluxFlag,ixComputeVegFlux)

 ! ********** RUN FOR ONE HRU ********************************************************************************************
 ! loop through HRUs
 do iHRU=1,gruInfo(1)%hruCount

  ! ----- hru initialization ---------------------------------------------------------------------------------------------


  ! ----- run the model --------------------------------------------------------------------------------------------------

  ! simulation for a single HRU
  call run_oneHRU(&
  nGRU, &
                  ! model control
                  gruInfo(1)%hruInfo(iHRU)%hru_nc,    & ! intent(in):    hru count Id
                  gruInfo(1)%hruInfo(iHRU)%hru_id,    & ! intent(in):    hruId
                  dt_init(iHRU),                   & ! intent(inout): initial time step
                  computeVegFluxFlag,              & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                  ! data structures (input)
                  timeVec,                         & ! intent(in):    model time data
                  typeHRU,               & ! intent(in):    local classification of soil veg etc. for each HRU
                  attrHRU,               & ! intent(in):    local attributes for each HRU
                  lookupHRU,             & ! intent(in):    local lookup tables for each HRU
                  bvarData,                        & ! intent(in):    basin-average model variables
                  decisions, veg_param, veg_tables, &
                  ! data structures (input-output)
                  mparHRU,               & ! intent(inout): model parameters
                  indxHRU,               & ! intent(inout): model indices
                  forcHRU,               & ! intent(inout): model forcing data
                  progHRU,               & ! intent(inout): model prognostic variables for a local HRU
                  diagHRU,               & ! intent(inout): model diagnostic variables for a local HRU
                  fluxHRU,               & ! intent(inout): model fluxes for a local HRU
                  ! error control
                  err,cmessage)                      ! intent(out):   error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
  ! update layer numbers that could be changed in run_oneHRU -- needed for model output
  do iGRU=1,nGRU
    gruInfo(iGRU)%hruInfo(iHRU)%nSnow = indxHRU%nSnow(iGRU)
    gruInfo(iGRU)%hruInfo(iHRU)%nSoil = indxHRU%nSoil
  end do

  call finalize_run_oneGRU<<<blocks,threads>>>(nGRU,decisions%spatial_gw,decisions%groundwatr,decisions%subRouting,&
  attrHRU%HRUarea,bvarData%basin__totalArea,&
  bvarData%basin__ColumnOutflow,fluxHRU%data,fluxHRU%ixMLayerColumnOutflow_start,fluxHRU%ixMLayerColumnOutflow_end,&
bvarData%basin__SurfaceRunoff,fluxHRU%ixScalarSurfaceRunoff,&
bvarData%basin__SoilDrainage,fluxHRU%ixScalarSoilDrainage,&
bvarData%basin__AquiferRecharge,fluxHRU%ixScalarAquiferRecharge,&
bvarData%basin__AquiferTranspire,fluxHRU%ixScalarAquiferTranspire,&
bvarData%basin__AquiferBaseflow,fluxHRU%ixScalarAquiferBaseflow,&
bvarData%basin__TotalRunoff,&
bvarData%routingFractionFuture,bvarData%routingRunoffFuture,&
bvarData%averageInstantRunoff,bvarData%averageRoutedRunoff,&
computeVegFluxFlag,ixComputeVegFlux)



 end do  ! (looping through HRUs)
 ! ********** END LOOP THROUGH HRUS **************************************************************************************
 
 end subroutine run_oneGRU


 attributes(global) subroutine finalize_run_oneGRU(nGRU,spatial_gw,groundwatr,subRouting,&
 HRUarea,basin__totalArea,&
 basin__ColumnOutflow,flux_data,mLayerColumnOutflow_start,mLayerColumnOutflow_end,&
 basin__SurfaceRunoff,scalarSurfaceRunoff,&
basin__SoilDrainage,scalarSoilDrainage,&
basin__AquiferRecharge,scalarAquiferRecharge,&
basin__AquiferTranspire,scalarAquiferTranspire,&
basin__AquiferBaseflow,scalarAquiferBaseflow,&
basin__TotalRunoff,&
routingFractionFuture,routingRunoffFuture,&
averageInstantRunoff,averageRoutedRunoff,&
computeVegFluxFlag,ixComputeVegFlux)
 implicit none
 integer(i4b),intent(in),value :: nGRU
 integer(i4b),intent(in) :: spatial_gw,groundwatr,subRouting
   real(rkind),intent(inout) :: HRUarea(:),basin__totalArea(:)
  real(rkind),intent(inout) :: basin__ColumnOutflow(:)
  real(rkind),intent(inout) :: flux_data(:,:)
  integer(i4b),intent(in),value :: mLayerColumnOutflow_start,mLayerColumnOutflow_end
  real(rkind),intent(inout) :: basin__SurfaceRunoff(:)
  integer(i4b),intent(in),value :: scalarSurfaceRunoff
  real(rkind),intent(inout) :: basin__SoilDrainage(:)
  integer(i4b),intent(in),value :: scalarSoilDrainage
  real(rkind),intent(inout) :: basin__AquiferRecharge(:)
  integer(i4b),intent(in),value :: scalarAquiferRecharge
  real(rkind),intent(inout) :: basin__AquiferTranspire(:)
  integer(i4b),intent(in),value :: scalarAquiferTranspire
  real(rkind),intent(inout) :: basin__AquiferBaseflow(:)
  integer(i4b),intent(in),value :: scalarAquiferBaseflow
real(rkind),intent(inout) :: basin__TotalRunoff(:)
real(rkind),intent(inout) :: routingFractionFuture(:,:),routingRunoffFuture(:,:)
real(rkind),intent(inout) :: averageInstantRunoff(:),averageRoutedRunoff(:)
logical(lgt),intent(inout) :: computeVegFluxFlag(:)
integer(i4b),intent(inout) :: ixComputeVegFlux(:)

    integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return
  call finalize_run_oneGRU_device(spatial_gw,groundwatr,subRouting,&
  HRUarea(iGRU),basin__totalArea(iGRU),&
  basin__ColumnOutflow(iGRU),flux_data(mLayerColumnOutflow_start:mLayerColumnOutflow_end,iGRU),&
basin__SurfaceRunoff(iGRU),flux_data(scalarSurfaceRunoff,iGRU),&
basin__SoilDrainage(iGRU),flux_data(scalarSoilDrainage,iGRU),&
basin__AquiferRecharge(iGRU),flux_data(scalarAquiferRecharge,iGRU),&
basin__AquiferTranspire(iGRU),flux_data(scalarAquiferTranspire,iGRU),&
basin__AquiferBaseflow(iGRU),flux_data(scalarAquiferBaseflow,iGRU),&
basin__TotalRunoff(iGRU),&
routingFractionFuture(:,iGRU),routingRunoffFuture(:,iGRU),&
averageInstantRunoff(iGRU),averageRoutedRunoff(iGRU),&
computeVegFluxFlag(iGRU),ixComputeVegFlux(iGRU))
end subroutine

attributes(device) subroutine finalize_run_oneGRU_device(spatial_gw,groundwatr,subRouting,&
HRUarea,basin__totalArea,&
basin__ColumnOutflow,mLayerColumnOutflow,&
basin__SurfaceRunoff,scalarSurfaceRunoff,&
basin__SoilDrainage,scalarSoilDrainage,&
basin__AquiferRecharge,scalarAquiferRecharge,&
basin__AquiferTranspire,scalarAquiferTranspire,&
basin__AquiferBaseflow,scalarAquiferBaseflow,&
basin__TotalRunoff,&
routingFractionFuture,routingRunoffFuture,&
averageInstantRunoff,averageRoutedRunoff,&
computeVegFluxFlag,ixComputeVegFlux)
 USE qTimeDelay_module,only:qOverland                        ! module to route water through an "unresolved" river network

implicit none
  real(rkind),intent(inout) :: HRUarea,basin__totalArea
  real(rkind),intent(inout) :: basin__ColumnOutflow
  real(rkind),intent(inout) :: mLayerColumnOutflow(:)
  real(rkind),intent(inout) :: basin__SurfaceRunoff,scalarSurfaceRunoff
  real(rkind),intent(inout) :: basin__SoilDrainage,scalarSoilDrainage
  integer(i4b),intent(in) :: spatial_gw,groundwatr,subRouting
  real(rkind),intent(inout) :: basin__AquiferRecharge,scalarAquiferRecharge
real(rkind),intent(inout) :: basin__AquiferTranspire,scalarAquiferTranspire
real(rkind),intent(inout) :: basin__AquiferBaseflow,scalarAquiferBaseflow
real(rkind),intent(inout) :: basin__TotalRunoff
real(rkind),intent(inout) :: routingFractionFuture(:), routingRunoffFuture(:)
real(rkind),intent(inout) :: averageInstantRunoff,averageRoutedRunoff
logical(lgt),intent(inout) :: computeVegFluxFlag
integer(i4b),intent(inout) :: ixComputeVegFlux

integer(i4b) :: err

  real(rkind) :: fracHRU

  ! save the flag for computing the vegetation fluxes
  if(computeVegFluxFlag)       ixComputeVegFlux = yes
  if(.not. computeVegFluxFlag) ixComputeVegFlux = no

  ! identify the area covered by the current HRU
  fracHRU = HRUarea / basin__totalArea
  ! (Note:  for efficiency, this could this be done as a setup task, not every timestep)
    
   ! otherwise just increment basin (GRU) column outflow (m3 s-1) with the hru fraction
   basin__ColumnOutflow = basin__ColumnOutflow + sum(mLayerColumnOutflow(:))
  ! ----- calculate weighted basin (GRU) fluxes --------------------------------------------------------------------------------------
  ! increment basin surface runoff (m s-1)
  basin__SurfaceRunoff  = basin__SurfaceRunoff + scalarSurfaceRunoff*fracHRU

  ! increment basin soil drainage (m s-1)
  basin__SoilDrainage   = basin__SoilDrainage  + scalarSoilDrainage *fracHRU

  ! increment aquifer variables -- ONLY if aquifer baseflow is computed individually for each HRU and aquifer is run
  ! NOTE: groundwater computed later for singleBasin
  if(spatial_gw == localColumn .and. groundwatr == bigBucket) then
   basin__AquiferRecharge  = basin__AquiferRecharge  + scalarAquiferRecharge *fracHRU
   basin__AquiferTranspire = basin__AquiferTranspire + scalarAquiferTranspire*fracHRU
   basin__AquiferBaseflow  = basin__AquiferBaseflow  + scalarAquiferBaseflow *fracHRU
  end if

  ! averaging more fluxes (and/or states) can be added to this section as desired

  ! perform the routing

 ! compute water balance for the basin aquifer
!  if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)then
!   message=trim(message)//'multi_driver/bigBucket groundwater code not transferred from old code base yet'
!   err=20; return
!  end if

 ! calculate total runoff depending on whether aquifer is connected
 if(groundwatr == bigBucket) then
  ! aquifer
  basin__TotalRunoff = basin__SurfaceRunoff + basin__ColumnOutflow/basin__totalArea + basin__AquiferBaseflow
 else
  ! no aquifer
  basin__TotalRunoff = basin__SurfaceRunoff + basin__ColumnOutflow/basin__totalArea + basin__SoilDrainage
 endif

 call qOverland(&
                ! input
                subRouting,          &  ! intent(in): index for routing method
                basin__TotalRunoff,             &  ! intent(in): total runoff to the channel from all active components (m s-1)
                routingFractionFuture,             &  ! intent(in): fraction of runoff in future time steps (m s-1)
                routingRunoffFuture,               &  ! intent(in): runoff in future time steps (m s-1)
                ! output
                averageInstantRunoff,           &  ! intent(out): instantaneous runoff (m s-1)
                averageRoutedRunoff,            &  ! intent(out): routed runoff (m s-1)
                err)                                                                  ! intent(out): error control
!  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif


end subroutine
 attributes(global) subroutine initialize_run_oneGRU(nGRU,basin__SurfaceRunoff,basin__SoilDrainage,&
basin__ColumnOutflow,basin__TotalRunoff,&
basin__AquiferRecharge,basin__AquiferBaseflow,basin__AquiferTranspire,&
flux_data,mLayerColumnInflow_start,mLayerColumnInflow_end,&
computeVegFluxFlag,ixComputeVegFlux)
 implicit none
 integer(i4b),intent(in),value :: nGRU
 real(rkind),intent(inout) :: basin__SurfaceRunoff(:), basin__SoilDrainage(:)
 real(rkind),intent(inout) :: basin__ColumnOutflow(:),basin__TotalRunoff(:)
 real(rkind),intent(inout) :: basin__AquiferRecharge(:),basin__AquiferBaseflow(:),basin__AquiferTranspire(:)
   real(rkind),intent(inout) :: flux_data(:,:)
  integer(i4b),intent(in),value :: mLayerColumnInflow_start,mLayerColumnInflow_end
logical(lgt),intent(inout) :: computeVegFluxFlag(:)
integer(i4b),intent(inout) :: ixComputeVegFlux(:)

    integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return
  call initialize_run_oneGRU_device(basin__SurfaceRunoff(iGRU),basin__SoilDrainage(iGRU),&
basin__ColumnOutflow(iGRU),basin__TotalRunoff(iGRU),&
basin__AquiferRecharge(iGRU),basin__AquiferBaseflow(iGRU),basin__AquiferTranspire(iGRU),&
flux_data(mLayerColumnInflow_start:mLayerColumnInflow_end,iGRU),&
computeVegFluxFlag(iGRU),ixComputeVegFlux(iGRU))
end subroutine

attributes(device) subroutine initialize_run_oneGRU_device(basin__SurfaceRunoff,basin__SoilDrainage,&
basin__ColumnOutflow,basin__TotalRunoff,&
basin__AquiferRecharge,basin__AquiferBaseflow,basin__AquiferTranspire,&
mLayerColumnInflow,&
computeVegFluxFlag,ixComputeVegFlux)
implicit none

real(rkind),intent(inout) :: basin__SurfaceRunoff,basin__SoilDrainage
real(rkind),intent(inout) :: basin__ColumnOutflow,basin__TotalRunoff
real(rkind),intent(inout) :: basin__AquiferRecharge,basin__AquiferBaseflow,basin__AquiferTranspire
real(rkind),intent(inout) :: mLayerColumnInflow(:)
logical(lgt),intent(inout) :: computeVegFluxFlag
integer(i4b),intent(inout) :: ixComputeVegFlux

 ! ----- basin initialization --------------------------------------------------------------------------------------------

 ! initialize runoff variables
 basin__SurfaceRunoff    = 0._rkind  ! surface runoff (m s-1)
 basin__SoilDrainage     = 0._rkind  ! soil drainage (m s-1)
 basin__ColumnOutflow    = 0._rkind  ! outflow from all "outlet" HRUs (those with no downstream HRU)
 basin__TotalRunoff      = 0._rkind  ! total runoff to the channel from all active components (m s-1)

 ! initialize baseflow variables
 basin__AquiferRecharge  = 0._rkind ! recharge to the aquifer (m s-1)
 basin__AquiferBaseflow  = 0._rkind ! baseflow from the aquifer (m s-1)
 basin__AquiferTranspire = 0._rkind ! transpiration loss from the aquifer (m s-1)

 ! initialize total inflow for each layer in a soil column
  mLayerColumnInflow(:) = 0._rkind

  ! set the flag to compute the vegetation flux
  computeVegFluxFlag = (ixComputeVegFlux == yes)

end subroutine

end module run_oneGRU_module
