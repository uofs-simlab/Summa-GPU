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

module computFlux_module

! data types
USE nrtype

! provide access to the derived types and classes used to define data structures and class objects
USE data_types,only:&
                    var_i,              & ! data vector (i4b)
                    var_d,              & ! data vector (rkind)
                    var_ilength,        & ! data vector with variable length dimension (i4b)
                    var_dlength,        & ! data vector with variable length dimension (rkind)
                    model_options,      & ! defines the model decisions
                    out_type_vegNrgFlux,                   & ! classes for vegNrgFlux call
                    out_type_ssdNrgFlux,& ! classes for ssdNrgFlux call
                    out_type_vegLiqFlux,                   & ! classes for vegLiqFlux call
                    out_type_snowLiqFlx,& ! classes for snowLiqFlx call                
                    out_type_soilLiqFlx,& ! classes for soilLiqFlx call
                    out_type_groundwatr,& ! classes for groundwatr call
                    out_type_bigAquifer   ! classes for bigAquifer call

! indices that define elements of the data structures
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookFORCE      ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! layer types
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil

! constants
USE multiconst,only:iden_water      ! intrinsic density of liquid water    (kg m-3)

! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:       &
 localColumn,                     & ! separate groundwater representation in each local soil column
 singleBasin                        ! single groundwater store over the entire basin

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:       &
 qbaseTopmodel,                   & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                       & ! a big bucket (lumped aquifer model)
 noExplicit                         ! no explicit groundwater parameterization

! look-up values for the form of Richards' equation
USE mDecisions_module,only:       &
 moisture,                        & ! moisture-based form of Richards' equation
 mixdform                           ! mixed form of Richards' equation

! look-up values for the choice of boundary conditions for hydrology
USE mDecisions_module,only:       &
 prescribedHead,                  & ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
 funcBottomHead,                  & ! function of matric head in the lower-most layer
 freeDrainage,                    & ! free drainage
 liquidFlux,                      & ! liquid water flux
 zeroFlux                           ! zero flux

implicit none
private
public::computFlux
public::soilCmpres
public::soilCmpresPrime,soilCmpresPrime_kernel

contains
! *********************************************************************************************************
! public subroutine computFlux: compute model fluxes
! *********************************************************************************************************
subroutine computFlux(&
                      ! input-output: model control
                      nSnow,                    & ! intent(in):    number of snow layers
                      nSoil,                    & ! intent(in):    number of soil layers
                      nLayers,                  & ! intent(in):    total number of layers
                      nGRU, &
                      firstSubStep,             & ! intent(in):    flag to indicate if we are processing the first sub-step
                      firstFluxCall,            & ! intent(inout): flag to denote the first flux call
                      firstSplitOper,           & ! intent(in):    flag to indicate if we are processing the first flux call in a splitting operation
                      computeVegFlux,           & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      scalarSolution,           & ! intent(in):    flag to indicate the scalar solution
                      checkLWBalance,           & ! intent(in):    flag to check longwave balance
                      drainageMeltPond,         & ! intent(in):    drainage from the surface melt pond (kg m-2 s-1)
                      ! input: state variables
                      scalarCanairTempTrial,    & ! intent(in):    trial value for the temperature of the canopy air space (K)
                      scalarCanopyTempTrial,    & ! intent(in):    trial value for the temperature of the vegetation canopy (K)
                      mLayerTempTrial,          & ! intent(in):    trial value for the temperature of each snow and soil layer (K)
                      mLayerMatricHeadLiqTrial, & ! intent(in):    trial value for the liquid water matric potential in each soil layer (m)
                      mLayerMatricHeadTrial,    & ! intent(in):    trial vector of total water matric potential (m)
                      scalarAquiferStorageTrial,& ! intent(in):    trial value of storage of water in the aquifer (m)
                      ! input: diagnostic variables defining the liquid water and ice content
                      scalarCanopyLiqTrial,     & ! intent(in):    trial value for the liquid water on the vegetation canopy (kg m-2)
                      scalarCanopyIceTrial,     & ! intent(in):    trial value for the ice on the vegetation canopy (kg m-2)
                      mLayerVolFracLiqTrial,    & ! intent(in):    trial value for the volumetric liquid water content in each snow and soil layer (-)
                      mLayerVolFracIceTrial,    & ! intent(in):    trial value for the volumetric ice in each snow and soil layer (-)
                      ! input: data structures
                      decisions,          & ! intent(in):    model decisions
                      model_decisions, &
                      type_data,                & ! intent(in):    type of vegetation and soil
                      attr_data,                & ! intent(in):    spatial attributes
                      mpar_data,                & ! intent(in):    model parameters
                      forc_data,                & ! intent(in):    model forcing data
                      bvar_data,                & ! intent(in):    average model variables for the entire basin
                      prog_data,                & ! intent(in):    model prognostic variables for a local HRU
                      indx_data,                & ! intent(in):    index data
                      veg_param, &
                      ! input-output: data structures
                      diag_data,                & ! intent(inout): model diagnostic variables for a local HRU
                      flux_data,                & ! intent(inout): model fluxes for a local HRU
                      deriv_data,               & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                      ! input-output: flux vector and baseflow derivatives
                      ixSaturation,             & ! intent(inout): index of the lowest saturated layer (NOTE: only computed on the first iteration)
                      dBaseflow_dMatric,        & ! intent(out):   derivative in baseflow w.r.t. matric head (s-1)
                      fluxVec,                  & ! intent(out):   flux vector (mixed units)
                      ! output: error control
                      err,message)                ! intent(out):   error code and error message
  ! provide access to flux subroutines
  USE vegNrgFlux_module,only:vegNrgFlux           ! compute energy fluxes over vegetation
  USE ssdNrgFlux_module,only:ssdNrgFlux           ! compute energy fluxes throughout the snow and soil subdomains
  USE vegLiqFlux_module,only:vegLiqFlux           ! compute liquid water fluxes through vegetation
  USE snowLiqFlx_module,only:snowLiqflx           ! compute liquid water fluxes through snow
  USE soilLiqFlx_module,only:soilLiqflx           ! compute liquid water fluxes through soil
  USE groundwatr_module,only:groundwatr           ! compute the baseflow flux
  USE bigAquifer_module,only:bigAquifer           ! compute fluxes for the big aquifer
  use cudafor
  use device_data_types
  implicit none
  ! -------------------------------------------------------------------------------------------------------------------------
  ! * dummy variables
  ! -------------------------------------------------------------------------------------------------------------------------
  ! input-output: control
  integer(i4b),intent(in),device            :: nSnow(:)                       ! number of snow layers
  integer(i4b),intent(in)            :: nSoil                       ! number of soil layers
  integer(i4b),intent(in),device            :: nLayers(:)                     ! total number of layers
  integer(i4b),intent(in) :: nGRU
  logical(lgt),intent(in)            :: firstSubStep                ! flag to indicate if we are processing the first sub-step
  logical(lgt),intent(inout)         :: firstFluxCall               ! flag to indicate if we are processing the first flux call
  logical(lgt),intent(in)            :: firstSplitOper              ! flag to indicate if we are processing the first flux call in a splitting operation
  logical(lgt),intent(in),device            :: computeVegFlux(:)              ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)            :: scalarSolution              ! flag to denote if implementing the scalar solution
  logical(lgt),intent(in)            :: checkLWBalance              ! flag to check longwave balance
  real(rkind),intent(in),device             :: drainageMeltPond(:)            ! drainage from the surface melt pond (kg m-2 s-1)
  ! input: state variables
  real(rkind),intent(in),device             :: scalarCanairTempTrial(:)       ! trial value for temperature of the canopy air space (K)
  real(rkind),intent(in),device             :: scalarCanopyTempTrial(:)       ! trial value for temperature of the vegetation canopy (K)
  real(rkind),intent(in),device             :: mLayerTempTrial(:,:)          ! trial value for temperature of each snow/soil layer (K)
  real(rkind),intent(in),device             :: mLayerMatricHeadLiqTrial(:,:) ! trial value for the liquid water matric potential (m)
  real(rkind),intent(in),device             :: mLayerMatricHeadTrial(:,:)    ! trial value for the total water matric potential (m)
  real(rkind),intent(in),device             :: scalarAquiferStorageTrial(:)   ! trial value of aquifer storage (m)
  ! input: diagnostic variables
  real(rkind),intent(in),device             :: scalarCanopyLiqTrial(:)        ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
  real(rkind),intent(in),device             :: scalarCanopyIceTrial(:)        ! trial value for mass of ice on the vegetation canopy (kg m-2)
  real(rkind),intent(in),device             :: mLayerVolFracLiqTrial(:,:)    ! trial value for volumetric fraction of liquid water (-)
  real(rkind),intent(in),device             :: mLayerVolFracIceTrial(:,:)    ! trial value for volumetric fraction of ice (-)
  ! input: data structures
  type(decisions_device),intent(in)     :: decisions          ! model decisions
  type(model_options),intent(in)  :: model_decisions(:)          ! model decisions
  type(type_data_device),        intent(in)     :: type_data                   ! type of vegetation and soil
  type(attr_data_device),        intent(in)     :: attr_data                   ! spatial attributes
  type(mpar_data_device),  intent(in)     :: mpar_data                   ! model parameters
  type(forc_data_device),        intent(in)     :: forc_data                   ! model forcing data
  type(bvar_data_device),  intent(in)     :: bvar_data                   ! model variables for the local basin
  type(prog_data_device),  intent(in)     :: prog_data                   ! prognostic variables for a local HRU
  type(indx_data_device),  intent(in)     :: indx_data                   ! indices defining model states and layers
  type(veg_parameters),intent(in) :: veg_param
  ! input-output: data structures
  type(diag_data_device),intent(inout)    :: diag_data                   ! diagnostic variables for a local HRU
  type(flux_data_device),intent(inout)    :: flux_data                   ! model fluxes for a local HRU
  type(deriv_data_device),intent(inout)    :: deriv_data                  ! derivatives in model fluxes w.r.t. relevant state variables
  ! input-output: flux vector and baseflow derivatives
  integer(i4b),intent(inout),device         :: ixSaturation(:)                ! index of the lowest saturated layer (NOTE: only computed on the first iteration)
  real(rkind),intent(out),device            :: dBaseflow_dMatric(:,:,:)      ! derivative in baseflow w.r.t. matric head (s-1)
  real(rkind),intent(out),device            :: fluxVec(:,:)                  ! model flux vector (mixed units)
  ! output: error control
  integer(i4b),intent(out)           :: err                         ! error code
  character(*),intent(out)           :: message                     ! error message
  ! -------------------------------------------------------------------------------------------------------------------------
  ! * local variables
  ! -------------------------------------------------------------------------------------------------------------------------
  integer(i4b)                       :: local_ixGroundwater         ! local index for groundwater representation
  integer(i4b)                       :: iLayer                      ! index of model layers
  logical(lgt)                       :: doVegNrgFlux                ! flag to compute the energy flux over vegetation
  real(rkind),device,dimension(nSoil,nGRU)       :: dHydCond_dMatric            ! derivative in hydraulic conductivity w.r.t matric head (s-1)
  character(LEN=256)                 :: cmessage                    ! error message of downwind routine
  ! ---------------------- classes for flux subroutine arguments (classes defined in data_types module) ----------------------
  !      ** intent(in) arguments **       ||       ** intent(inout) arguments **        ||      ** intent(out) arguments **
  type(out_type_vegNrgFlux) :: out_vegNrgFlux ! vegNrgFlux arguments
  type(out_type_ssdNrgFlux) :: out_ssdNrgFlux ! ssdNrgFlux arguments
  type(out_type_vegLiqFlux) :: out_vegLiqFlux ! vegLiqFlux arguments
  type(out_type_snowLiqFlx) :: out_snowLiqFlx ! snowLiqFlx arguments
  type(out_type_soilLiqFlx) :: out_soilLiqFlx ! soilLiqFlx arguments
  type(out_type_groundwatr) :: out_groundwatr ! groundwatr arguments
  type(out_type_bigAquifer) :: out_bigAquifer ! bigAquifer arguments
  ! -------------------------------------------------------------------------------------------------------------------------

  ! initialize error control
  err=0; message='computFlux/'

  call initialize_computFlux; if(err/=0)then; return; endif ! Preliminary operations to start routine

  ! *** CALCULATE ENERGY FLUXES OVER VEGETATION ***
    ! identify the need to calculate the energy flux over vegetation
    ! if (doVegNrgFlux) then ! if necessary, calculate the energy fluxes over vegetation
      call initialize_vegNrgFlux
      call vegNrgFlux(nGRU,firstSubStep,firstFluxCall,computeVegFlux,checkLWBalance,&
      scalarCanairTempTrial,scalarCanopyTempTrial,mLayerTempTrial,scalarCanopyIceTrial,scalarCanopyLiqTrial,&
      type_data,forc_data,mpar_data,indx_data,prog_data,diag_data,flux_data,bvar_data,deriv_data,decisions,veg_param,out_vegNrgFlux)
      call finalize_vegNrgFlux; if(err/=0)then; return; endif
    ! end if

  ! *** CALCULATE ENERGY FLUXES THROUGH THE SNOW-SOIL DOMAIN ***
      call initialize_ssdNrgFlux
      call ssdNrgFlux(nGRU,mLayerTempTrial,decisions,mpar_data,indx_data,prog_data,diag_data,flux_data,deriv_data,out_ssdNrgFlux)
      call finalize_ssdNrgFlux; if(err/=0)then; return; endif



  ! *** CALCULATE THE LIQUID FLUX THROUGH VEGETATION ***
      call initialize_vegLiqFlux
      call vegLiqFlux(nGRU,computeVegFlux,scalarCanopyLiqTrial,decisions,indx_data,mpar_data,diag_data,flux_data,deriv_data,out_vegLiqFlux)
      call finalize_vegLiqFlux; if(err/=0)then; return; endif

  ! *** CALCULATE THE LIQUID FLUX THROUGH SNOW ***
    ! if (nSnowOnlyHyd>0) then ! if necessary, compute liquid fluxes through snow
      call initialize_snowLiqFlx
      call snowLiqFlx(nGRU,firstFluxCall,nSnow,mLayerVolFracLiqTrial,&
      indx_data,mpar_data,prog_data,diag_data,flux_data,deriv_data,out_snowLiqFlx)
      call finalize_snowLiqFlx; if(err/=0)then; return; endif
    ! else
      call soilForcingNoSnow ! define forcing for the soil domain for the case of no snow layers
    ! end if


  ! *** CALCULATE THE LIQUID FLUX THROUGH SOIL ***
      call initialize_soilLiqFlx
      call soilLiqFlx(nGRU,decisions,mpar_data,indx_data,prog_data,diag_data,flux_data,out_soilLiqFlx,deriv_data,dHydCond_dMatric,nSoil,firstSplitOper,(scalarSolution .and. .not.firstFluxCall),&
      mLayerMatricHeadTrial,mLayerMatricHeadLiqTrial,mLayerVolFracIceTrial,nSnow,mLayerVolFracLiqTrial,nLayers,mLayerTempTrial)
      call finalize_soilLiqFlx; if(err/=0)then; return; endif


  ! *** CALCULATE THE GROUNDWATER FLOW ***
      if (local_ixGroundwater/=qbaseTopmodel) then ! set baseflow fluxes to zero if the topmodel baseflow routine is not used
        call zeroBaseflowFluxes
      else ! compute the baseflow flux for topmodel-ish shallow groundwater
        call initialize_groundwatr; if(err/=0)then; return; endif
        call groundwatr(nGRU,nSnow,nSoil,nLayers,firstFluxCall,mLayerVolFracLiqTrial,mLayerVolFracIceTrial,&
        attr_data,mpar_data,prog_data,diag_data,flux_data,deriv_data,&
        ixSaturation,out_groundwatr,dBaseflow_dMatric)
        call finalize_groundwatr;   if(err/=0)then; return; endif
      end if
      call computeBaseflowRunoff ! compute total baseflow from soil and runoff


  ! *** CALCULATE FLUXES FOR THE DEEP AQUIFER ***
      if (local_ixGroundwater==bigBucket) then ! compute fluxes for the big bucket
        call bigAquifer(nGRU,scalarAquiferStorageTrial,mpar_data,diag_data,indx_data,flux_data,deriv_data,out_bigAquifer)
        call finalize_bigAquifer; if(err/=0)then; return; endif
      else ! if no aquifer, then fluxes are zero
        call zeroAquiferFluxes
      end if ! end check aquifer model decision
          
  call finalize_computFlux; if(err/=0)then; return; endif ! final operations to prep for end of routine

contains

 ! **** Subroutines that handle the absence of model features ****
 subroutine soilForcingNoSnow
  ! define forcing for the soil domain for the case of no snow layers
  ! NOTE: in case where nSnowOnlyHyd==0 AND snow layers exist, then scalarRainPlusMelt is taken from the previous flux evaluation
      type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  call soilForcingNoSnow_kernel<<<blocks,threads>>>(nGRU,nSnow,&
    flux_data%ixScalarRainPlusMelt,flux_data%ixScalarThroughfallRain,flux_data%ixScalarCanopyLiqDrainage,drainageMeltPond,flux_data%data)
 end subroutine soilForcingNoSnow

 subroutine zeroBaseflowFluxes
    type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  call zeroBaseflowFluxes_kernel<<<blocks,threads>>>(nGRU,nSoil,flux_data%data,flux_data%ixScalarExfiltration,flux_data%ixMLayerColumnOutflow_start,flux_data%ixMLayerBaseflow_start)
 end subroutine zeroBaseflowFluxes

 subroutine computeBaseflowRunoff
  ! compute total baseflow from the soil zone (needed for mass balance checks) and total runoff
  ! (Note: scalarSoilBaseflow is zero if topmodel is not used)
  ! (Note: scalarSoilBaseflow may need to re-envisioned in topmodel formulation if parts of it flow into neighboring soil rather than exfiltrate)

      type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  call computeBaseflowRunoff_kernel<<<blocks,threads>>>(nGRU,nSoil,&
  flux_data%data, &
  flux_data%ixScalarSoilBaseflow,flux_data%ixMLayerBaseflow_start,flux_data%ixScalarTotalRunoff,flux_data%ixScalarSurfaceRunoff,flux_data%ixScalarSoilDrainage)
 end subroutine computeBaseflowRunoff  

 subroutine zeroAquiferFluxes
      type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)
call zeroAquiferFluxes_kernel<<<blocks,threads>>>(nGRU,indx_data%ixAqWat,flux_data%data,flux_data%ixScalarAquiferTranspire,flux_data%ixScalarAquiferRecharge,flux_data%ixScalarAquiferBaseflow,deriv_data%dBaseflow_dAquifer)
 end subroutine zeroAquiferFluxes

 ! **** Subroutines for starting/ending operations of computFlux ****
 subroutine initialize_computFlux
    type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
    blocks = dim3(nGRU/128+1,1,1)

  ! operations to prep for the start of computFlux
  associate(&
   numFluxCalls                 => diag_data%numFluxCalls,         & ! intent(out): [dp] number of flux calls (-)
   ixSpatialGroundwater         => model_decisions(iLookDECISIONS%spatial_gw)%iDecision, & ! intent(in): [i4b] spatial representation of groundwater (local-column or single-basin)
   ixGroundwater                => model_decisions(iLookDECISIONS%groundwatr)%iDecision, & ! intent(in): [i4b] groundwater parameterization
   iLayerLiqFluxSnow            => flux_data%ixILayerLiqFluxSnow_start,       & ! intent(out): [dp(0:)] vertical liquid water flux at snow layer interfaces (-)
   iLayerLiqFluxSoil            => flux_data%ixILayerLiqFluxSoil_start        ) ! intent(out): [dp(0:)] vertical liquid water flux at soil layer interfaces (-)

   numFluxCalls = numFluxCalls+1 ! increment the number of flux calls

   ! modify the groundwater representation for this single-column implementation
   select case(ixSpatialGroundwater)
     case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
     case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
     case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
   end select ! end modify the groundwater representation for this single-column implementation

   ! initialize liquid water fluxes throughout the snow and soil domains
   ! NOTE: used in the energy routines, which is called before the hydrology routines
   if (firstFluxCall) then
     call initialize_computeFlux_kernel<<<blocks,threads>>>(nGRU,nSoil,nSnow,flux_data%data,iLayerLiqFluxSnow,iLayerLiqFluxSoil)
   end if
  end associate
 end subroutine initialize_computFlux

 subroutine finalize_computFlux
      type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
    blocks = dim3(nGRU/128+1,1,1)

    call finalize_computFlux_kernel<<<blocks,threads>>>(nGRU,nLayers,&
    indx_data%ixCasNrg,indx_data%ixVegNrg,indx_data%ixVegHyd,&
    flux_data%data,flux_data%ixScalarCanairNetNrgFlux,flux_data%ixScalarCanopyNetNrgFlux,flux_data%ixScalarCanopyNetLiqFlux,diag_data%scalarCanopyDepth,&
    indx_data%ixSnowSoilNrg,flux_data%ixmLayerNrgFlux_start,&
    fluxVec,&
    indx_data%ixAqWat,indx_data%ixSnowSoilHyd,indx_data%layerType,nSnow,&
    flux_data%ixMLayerLiqFluxSnow_start,flux_data%ixMLayerLiqFluxSoil_start,&
    flux_data%ixScalarAquiferTranspire,flux_data%ixScalarAquiferRecharge,flux_data%ixScalarAquiferBaseflow)

   firstFluxCall=.false. ! set the first flux call to false
 end subroutine finalize_computFlux

 ! ----------------------- Initialize and Finalize procedures for the flux routines -----------------------
 ! **** vegNrgFlux ****
 subroutine initialize_vegNrgFlux
        type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  call initialize_vegNrgFlux_kernel<<<blocks,threads>>>(nGRU,indx_data%ixCasNrg,indx_data%ixVegNrg,indx_data%ixTopNrg,&
  deriv_data%dCanLiq_dTcanopy,deriv_data%dTheta_dTkCanopy,diag_data%scalarCanopyDepth)
 end subroutine initialize_vegNrgFlux

 subroutine finalize_vegNrgFlux
  call out_vegNrgFlux%finalize(err,cmessage)
 ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if  ! check for errors
 end subroutine finalize_vegNrgFlux
 ! **** end vegNrgFlux ****

 ! **** ssdNrgFlux ****
 subroutine initialize_ssdNrgFlux
  ! call in_ssdNrgFlux%initialize(scalarSolution,firstFluxCall,mLayerTempTrial,flux_data,deriv_data)
 end subroutine initialize_ssdNrgFlux

 subroutine finalize_ssdNrgFlux
  ! type(dim3) :: blocks,threads
  ! threads = dim3(128,1,1)
  ! blocks = dim3(nGRU/128+1,1,1)
        type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  call out_ssdNrgFlux%finalize(err,cmessage)
   call finalize_ssdNrgFlux_kernel<<<blocks,threads>>>(nGRU,nLayers,&
  flux_data%data,flux_data%ixmLayerNrgFlux_start,flux_data%ixmlayerNrgFlux_end,flux_data%ixiLayerNrgFlux_start,flux_data%ixiLayerNrgFlux_end,prog_data%mLayerDepth)

 end subroutine finalize_ssdNrgFlux
 ! **** end ssdNrgFlux ****

 ! **** vegLiqFlux ****
 subroutine initialize_vegLiqFlux
  ! call in_vegLiqFlux%initialize(computeVegFlux,scalarCanopyLiqTrial,flux_data)
 end subroutine initialize_vegLiqFlux
 
 subroutine finalize_vegLiqFlux
        type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  call out_vegLiqFlux%finalize(err,cmessage)
  associate( &
    ixVegHyd => indx_data%ixVegHyd, &
   scalarThroughfallRain        => flux_data%ixScalarThroughfallRain,         & ! intent(out): [dp] rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
   scalarCanopyLiqDrainage      => flux_data%ixScalarCanopyLiqDrainage,       & ! intent(out): [dp] drainage of liquid water from the vegetation canopy (kg m-2 s-1)
   scalarThroughfallRainDeriv   => deriv_data%scalarThroughfallRainDeriv  ,& ! intent(out): [dp] derivative in throughfall w.r.t. canopy liquid water
   scalarCanopyLiqDrainageDeriv => deriv_data%scalarCanopyLiqDrainageDeriv,& ! intent(out): [dp] derivative in canopy drainage w.r.t. canopy liquid water
   scalarCanopyNetLiqFlux       => flux_data%ixScalarCanopyNetLiqFlux,        & ! intent(out): [dp] net liquid water flux for the vegetation canopy (kg m-2 s-1)
   scalarRainfall               => flux_data%ixScalarRainfall,                & ! intent(in):  [dp] rainfall rate (kg m-2 s-1)
   scalarCanopyEvaporation      => flux_data%ixScalarCanopyEvaporation,       & ! intent(out): [dp] canopy evaporation/condensation (kg m-2 s-1)
   scalarCanopyLiqDeriv         => deriv_data%scalarCanopyLiqDeriv         ) ! intent(out): [dp] derivative in (throughfall + drainage) w.r.t. canopy liquid water
   ! error control
   if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
  !  ! calculate the net liquid water flux for the vegetation canopy
  !  scalarCanopyNetLiqFlux = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
  !  ! calculate the total derivative in the downward liquid flux
  !  scalarCanopyLiqDeriv   = scalarThroughfallRainDeriv + scalarCanopyLiqDrainageDeriv

   call finalize_vegLiqFlux_kernel<<<blocks,threads>>>(nGRU,ixVegHyd,&
  scalarThroughfallRain,scalarCanopyLiqDrainage,scalarThroughfallRainDeriv,&
  scalarCanopyLiqDrainageDeriv,flux_data%data,scalarCanopyNetLiqFlux,&
  scalarRainfall,scalarCanopyEvaporation,scalarCanopyLiqDeriv)

  end associate
 end subroutine finalize_vegLiqFlux
 ! **** end vegLiqFlux ****

 ! **** snowLiqFlx ****
 subroutine initialize_snowLiqFlx
 end subroutine initialize_snowLiqFlx

 subroutine finalize_snowLiqFlx
      type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  call finalize_snowLiqFlx_kernel<<<blocks,threads>>>(nGRU,nSnow,&
  flux_data%data, &
  flux_data%ixScalarRainPlusMelt,flux_data%ixMLayerLiqFluxSnow_start,flux_data%ixILayerLiqFluxSnow_start,prog_data%mLayerDepth,flux_data%ixScalarSnowDrainage)
  call out_snowLiqFlx%finalize(err,cmessage) 
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
 end subroutine finalize_snowLiqFlx
 ! **** end snowLiqFlx ****

 ! **** soilLiqFlx ****
 subroutine initialize_soilLiqFlx
 end subroutine initialize_soilLiqFlx

 subroutine finalize_soilLiqFlx
  ! call io_soilLiqFlx%finalize(nSoil,dHydCond_dMatric,flux_data,diag_data,deriv_data)
  type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)


  call finalize_soilLiqFlx_kernel<<<blocks,threads>>>(nGRU,nSoil,nSnow,&
  flux_data%data, &
  flux_data%ixMLayerLiqFluxSoil_start,flux_data%ixILayerLiqFluxSoil_start,prog_data%mLayerDepth,flux_data%ixScalarSoilDrainage,&
  deriv_data%dq_dHydStateAbove_m,deriv_data%dq_dHydStateBelow_m,deriv_data%dq_dHydStateLayerSurfVec_m,deriv_data%dPsiLiq_dPsi0_m)
 end subroutine finalize_soilLiqFlx
 ! **** end soilLiqFlx ****

 ! **** groundwatr ****
 subroutine initialize_groundwatr
  ! check the derivative matrix is sized appropriately
  if (size(dBaseflow_dMatric,1)/=nSoil .or. size(dBaseflow_dMatric,2)/=nSoil) then
    message=trim(message)//'expect dBaseflow_dMatric to be nSoil x nSoil'
    err=20; return
  end if
  ! call in_groundwatr%initialize(nSnow,nSoil,nLayers,firstFluxCall,mLayerVolFracLiqTrial,mLayerVolFracIceTrial,deriv_data)
 end subroutine initialize_groundwatr

 subroutine finalize_groundwatr
  call out_groundwatr%finalize(err,cmessage)
  ! error control
  if (err/=0) then; message=trim(message)//trim(cmessage); return; end if
 end subroutine finalize_groundwatr
 ! **** end groundwatr ****

 ! **** bigAquifer ****
 subroutine initialize_bigAquifer
 end subroutine initialize_bigAquifer

 subroutine finalize_bigAquifer
    type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)
call finalize_bigAquifer_kernel<<<blocks,threads>>>(nGRU,indx_data%ixAqWat,flux_data%data,flux_data%ixScalarTotalRunoff,flux_data%ixScalarSurfaceRunoff,flux_data%ixScalarAquiferBaseflow)
 end subroutine finalize_bigAquifer
 ! **** end bigAquifer ****

end subroutine computFlux
attributes(global) subroutine initialize_computeFlux_kernel(nGRU,nSoil,nSnow,flux_data,iLayerLiqFluxSnow,iLayerLiqFluxSoil)
  integer(i4b),value :: nGRU,nSoil
  integer(i4b),intent(in) :: nSnow(:)
  real(rkind) :: flux_data(:,:)
  integer(i4b),value :: iLayerLiqFluxSnow, iLayerLiqFluxSoil
  integer(i4b) :: iGRU,iLayer
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

  do iLayer=0,nSnow(iGRU)
    flux_data(iLayerLiqFluxSnow+iLayer,iGRU) = 0._rkind
  end do
  do iLayer=0,nSoil
     flux_data(iLayerLiqFluxSoil+iLayer,iGRU) = 0._rkind
  end do
end subroutine

attributes(global) subroutine finalize_soilLiqFlx_kernel(nGRU,nSoil,nSnow,&
  flux_data, &
  mLayerLiqFluxSoil,iLayerLiqFluxSoil,mLayerDepth,scalarSoilDrainage,&
  dq_dHydStateAbove,dq_dHydStateBelow,dq_dHydStateLayerSurfVec,dPsiLiq_dPsi0)
  integer(i4b),value :: nGRU,nSoil
  integer(i4b),intent(in) :: nSnow(:)
  real(rkind),intent(inout) :: flux_data(:,:)
  integer(i4b),value :: mLayerLiqFluxSoil,iLayerLiqFluxSoil
  real(rkind) :: mLayerDepth(:,:)
  integer(i4b),value :: scalarSoilDrainage
  real(rkind),intent(inout) :: dq_dHydStateAbove(0:,:),dq_dHydStateBelow(0:,:),dq_dHydStateLayerSurfVec(:,:),dPsiLiq_dPsi0(:,:)
  integer(i4b) :: iGRU,iLayer
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

   ! calculate net liquid water fluxes for each soil layer (s-1)
   do iLayer=1,nSoil
     flux_data(mLayerLiqFluxSoil+iLayer-1,iGRU) = -(flux_data(iLayerLiqFluxSoil+iLayer,iGRU) - flux_data(iLayerLiqFluxSoil+iLayer-1,iGRU))/mLayerDepth(iLayer+nSnow(iGRU),iGRU)
   end do
   ! compute drainage from the soil zone (needed for mass balance checks and in aquifer recharge)
   flux_data(scalarSoilDrainage,iGRU) = flux_data(iLayerLiqFluxSoil+nSoil,iGRU)

   do iLayer=1,nSoil
    dq_dHydStateAbove(iLayer,iGRU)   = dq_dHydStateAbove(iLayer,iGRU)  *dPsiLiq_dPsi0(iLayer,iGRU)
    dq_dHydStateBelow(iLayer-1,iGRU) = dq_dHydStateBelow(iLayer-1,iGRU)*dPsiLiq_dPsi0(iLayer,iGRU)
   end do
   if(all(dq_dHydStateLayerSurfVec(:,iGRU)/=realMissing)) dq_dHydStateLayerSurfVec(1:nSoil,iGRU) = dq_dHydStateLayerSurfVec(1:nSoil,iGRU)*dPsiLiq_dPsi0(1:nSoil,iGRU)

  end subroutine

  attributes(global) subroutine zeroBaseflowFluxes_kernel(nGRU,nSoil,flux_data,scalarExfiltration,mLayerColumnOutflow,mLayerBaseflow)
  integer(i4b),value :: nGRU,nSoil
  real(rkind),intent(inout) :: flux_data(:,:)
  integer(i4b),intent(in),value :: scalarExfiltration
  integer(i4b),intent(in),value :: mLayerColumnOutflow, mLayerBaseflow
  integer(i4b) :: iGRU,iLayer
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

   ! diagnostic variables in the data structures
   flux_data(scalarExfiltration,iGRU)     = 0._rkind  ! exfiltration from the soil profile (m s-1)
   do iLayer=1,nSoil
   flux_data(mLayerColumnOutflow+iLayer-1,iGRU) = 0._rkind  ! column outflow from each soil layer (m3 s-1)
   ! variables needed for the numerical solution
   flux_data(mLayerBaseflow+iLayer-1,iGRU)      = 0._rkind  ! baseflow from each soil layer (m s-1)
   end do

  end subroutine

  attributes(global) subroutine computeBaseflowRunoff_kernel(nGRU,nSoil,&
  flux_data,&
  scalarSoilBaseflow,mLayerBaseflow,scalarTotalRunoff,scalarSurfaceRunoff,scalarSoilDrainage)
    integer(i4b),value :: nGRU, nSoil
    real(rkind) :: flux_data(:,:)
    integer(i4b),value :: scalarSoilBaseflow
    integer(i4b),value :: mLayerBaseflow
    integer(i4b),value :: scalarTotalRunoff, scalarSurfaceRunoff, scalarSoilDrainage
  integer(i4b) :: iGRU,iLayer
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

  flux_data(scalarSoilBaseflow,iGRU) = 0._rkind
  do iLayer=1,nSoil
       flux_data(scalarSoilBaseflow,iGRU) = flux_data(scalarSoilBaseflow,iGRU) + flux_data(mLayerBaseflow+iLayer-1,iGRU)                                               ! baseflow from the soil zone 
  end do
   flux_data(scalarTotalRunoff,iGRU)  = flux_data(scalarSurfaceRunoff,iGRU) + flux_data(scalarSoilDrainage,iGRU) + flux_data(scalarSoilBaseflow,iGRU)     ! total runoff
end subroutine

  attributes(global) subroutine finalize_bigAquifer_kernel(nGRU,ixAqWat,&
  flux_data, &
  scalarTotalRunoff,scalarSurfaceRunoff,scalarAquiferBaseflow)
  use initialize_device,only:get_iGRU
  integer(i4b),value :: nGRU
  integer(i4b),intent(in) :: ixAqWat(:)
  real(rkind),intent(inout) :: flux_data(:,:)
  integer(i4b),intent(in),value :: scalarTotalRunoff, scalarSurfaceRunoff, scalarAquiferBaseflow
  integer(i4b) :: iGRU
  iGRU = get_iGRU()
   ! compute total runoff (overwrite previously calculated value before considering aquifer).
   !   (Note:  SoilDrainage goes into aquifer, not runoff)
  if (iGRU .gt. nGRU) return

  if (ixAqWat(iGRU) .eq. integerMissing) return
   flux_data(scalarTotalRunoff,iGRU)  = flux_data(scalarSurfaceRunoff,iGRU) + flux_data(scalarAquiferBaseflow,iGRU)     


end subroutine

  attributes(global) subroutine finalize_vegLiqFlux_kernel(nGRU,ixVegHyd,&
  scalarThroughfallRain,scalarCanopyLiqDrainage,scalarThroughfallRainDeriv,&
  scalarCanopyLiqDrainageDeriv,flux_data,scalarCanopyNetLiqFlux,&
  scalarRainfall,scalarCanopyEvaporation,scalarCanopyLiqDeriv)
  integer(i4b),value :: nGRU
  integer(i4b),intent(in) :: ixVegHyd(:)
  integer(i4b),intent(in),value :: scalarThroughfallRain
  integer(i4b),intent(in),value :: scalarCanopyLiqDrainage
  real(rkind),intent(inout) :: scalarThroughfallRainDeriv(:)
  real(rkind),intent(inout) :: scalarCanopyLiqDrainageDeriv(:)
  real(rkind),intent(inout) :: flux_data(:,:)
  integer(i4b),intent(in),value :: scalarCanopyNetLiqFlux
  integer(i4b),intent(in),value :: scalarRainfall
  integer(i4b),intent(in),value :: scalarCanopyEvaporation
  real(rkind),intent(inout) :: scalarCanopyLiqDeriv(:)
  integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
   ! compute total runoff (overwrite previously calculated value before considering aquifer).
   !   (Note:  SoilDrainage goes into aquifer, not runoff)
  if (iGRU .gt. nGRU) return

  if (ixVegHyd(iGRU) .eq. integerMissing) return
   ! calculate the net liquid water flux for the vegetation canopy
   flux_data(scalarCanopyNetLiqFlux,iGRU) = flux_data(scalarRainfall,iGRU) + flux_data(scalarCanopyEvaporation,iGRU) - flux_data(scalarThroughfallRain,iGRU) - flux_data(scalarCanopyLiqDrainage,iGRU)
   ! calculate the total derivative in the downward liquid flux
   scalarCanopyLiqDeriv(iGRU)   = scalarThroughfallRainDeriv(iGRU) + scalarCanopyLiqDrainageDeriv(iGRU)


end subroutine


  attributes(global) subroutine zeroAquiferFluxes_kernel(nGRU,ixAqWat,&
  flux_data, &
  scalarAquiferTranspire,scalarAquiferRecharge,scalarAquiferBaseflow,dBaseflow_dAquifer)
  use initialize_device,only:get_iGRU
  integer(i4b),value :: nGRU
  integer(i4b),intent(in) :: ixAqWat(:)
  real(rkind),intent(inout) :: flux_data(:,:)
  integer(i4b),intent(in),value :: scalarAquiferTranspire, scalarAquiferRecharge, scalarAquiferBaseflow
  real(rkind),intent(inout) :: dBaseflow_dAquifer(:)
  integer(i4b) :: iGRU
  iGRU = get_iGRU()
   ! compute total runoff (overwrite previously calculated value before considering aquifer).
   !   (Note:  SoilDrainage goes into aquifer, not runoff)
  if (iGRU .gt. nGRU) return

  if (ixAqWat(iGRU) .eq. integerMissing) return
   ! set aquifer fluxes to zero (if no aquifer exists)
   flux_data(scalarAquiferTranspire,iGRU) = 0._rkind  ! transpiration loss from the aquifer (m s-1)
   flux_data(scalarAquiferRecharge,iGRU)  = 0._rkind  ! recharge to the aquifer (m s-1)
   flux_data(scalarAquiferBaseflow,iGRU)  = 0._rkind  ! total baseflow from the aquifer (m s-1)
   dBaseflow_dAquifer(iGRU)     = 0._rkind  ! change in baseflow flux w.r.t. aquifer storage (s-1)


end subroutine

attributes(global) subroutine finalize_snowLiqFlx_kernel(nGRU,nSnow,&
flux_data, &
  scalarRainPlusMelt,mLayerLiqFluxSnow,iLayerLiqFluxSnow,mLayerDepth,scalarSnowDrainage)
  integer(i4b),value :: nGRU
  integer(i4b) :: nSnow(:)
  real(rkind) :: flux_data(:,:)
  integer(i4b),value :: scalarRainPlusMelt
  integer(i4b),value :: mLayerLiqFluxSnow
  integer(i4b),value :: iLayerLiqFluxSnow
  real(rkind) :: mLayerDepth(:,:)
  integer(i4b),value :: scalarSnowDrainage
  integer(i4b) :: iGRU,iLayer
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

  if (nSnow(iGRU).eq.0) return
     ! define forcing for the soil domain
   flux_data(scalarRainPlusMelt,iGRU) = flux_data(iLayerLiqFluxSnow+nSnow(iGRU),iGRU)          ! drainage from the base of the snowpack
   ! calculate net liquid water fluxes for each snow layer (s-1)
   do iLayer=1,nSnow(iGRU)
     flux_data(mLayerLiqFluxSnow+iLayer-1,iGRU) = -(flux_data(iLayerLiqFluxSnow+iLayer,iGRU) - flux_data(iLayerLiqFluxSnow+iLayer-1,iGRU))/mLayerDepth(iLayer,iGRU)
   end do
   ! compute drainage from the soil zone (needed for mass balance checks)
   flux_data(scalarSnowDrainage,iGRU) = flux_data(iLayerLiqFluxSnow+nSnow(iGRU),iGRU)
  end subroutine

  attributes(global) subroutine soilForcingNoSnow_kernel(nGRU,nSnow,&
    scalarRainPlusMelt,scalarThroughfallRain,scalarCanopyLiqDrainage,drainageMeltPond,flux_data)
    integer(i4b),value :: nGRU
    integer(i4b) :: nSnow(:)
    real(rkind) :: drainageMeltPond(:)
    integer(i4b),intent(in),value :: scalarCanopyLiqDrainage,scalarThroughfallRain,scalarRainPlusMelt
    real(rkind),intent(in) :: flux_data(:,:)

    integer(i4b) :: iGRU
    iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

    if (iGRU .gt. nGRU) return

    if (nSnow(iGRU)==0) then !no snow layers
      flux_data(scalarRainPlusMelt,iGRU) = (flux_data(scalarThroughfallRain,iGRU) + flux_data(scalarCanopyLiqDrainage,iGRU))/iden_water &  ! liquid flux from the canopy (m s-1)
                       + drainageMeltPond(iGRU)/iden_water  ! melt of the snow without a layer (m s-1)
    end if ! snow layers or not
  end subroutine

attributes(global) subroutine finalize_ssdNrgFlux_kernel(nGRU,nLayers,&
flux_data, &
  mLayerNrgFlux_start,mLayerNrgFlux_end,iLayerNrgFlux_start,iLayerNrgFlux_end,mLayerDepth)
  integer(i4b),value :: nGRU
  integer(i4b),intent(in) :: nLayers(:)
  real(rkind) :: flux_data(:,:)
  integer(i4b),intent(in),value :: mLayerNrgFlux_start,mLayerNrgFlux_end,iLayerNrgFlux_start,iLayerNrgFlux_end
  real(rkind) :: mLayerDepth(:,:)
  integer(i4b) :: iGRU,iLayer
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

   ! calculate net energy fluxes for each snow and soil layer (J m-3 s-1)
   do iLayer=1,nLayers(iGRU)
     flux_data(mLayerNrgFlux_start+iLayer-1,iGRU) = -(flux_data(iLayerNrgFlux_start+iLayer,iGRU) - flux_data(iLayerNrgFlux_start+iLayer-1,iGRU))/mLayerDepth(iLayer,iGRU)
    !  if (iGRU .eq. 1) print*, iLayer, mLayerNrgFLux(iLayer,iGRU),iLayerNrgFlux(iLayer,iGRU),iLayerNrgFlux(iLayer-1,iGRU),mLayerDepth(iLayer,iGRU)
   end do

  end subroutine

attributes(global) subroutine initialize_vegNrgFlux_kernel(nGRU,&
  ixCasNrg, ixVegNrg, ixTopNrg, &
  dCanLiq_dTcanopy, dTheta_dTkCanopy,canopyDepth)
  integer(i4b),value :: nGRU
  integer(i4b) :: ixCasNrg(:), ixVegNrg(:), ixTopNrg(:)
  real(rkind) :: dCanLiq_dTcanopy(:), dTheta_dTkCanopy(:), canopyDepth(:)
  integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

  if (ixCasNrg(iGRU)/=integerMissing .or. ixVegNrg(iGRU)/=integerMissing .or. ixTopNrg(iGRU)/=integerMissing) then
   dCanLiq_dTcanopy(iGRU) = dTheta_dTkCanopy(iGRU)*iden_water*canopyDepth(iGRU)     ! derivative in canopy liquid storage w.r.t. canopy temperature (kg m-2 K-1)
  end if
  end subroutine


  attributes(global) subroutine finalize_computFlux_kernel(nGRU,nLayers_,&
    ixCasNrg_,ixVegNrg_,ixVegHyd_,&
    flux_data,scalarCanairNetNrgFlux,scalarCanopyNetNrgFlux,scalarCanopyNetLiqFlux,canopyDepth_,&
    ixSnowSoilNrg_,mLayerNrgFlux_start, &
    fluxVec_,&
    ixAqWat_,ixSnowSoilHyd_,layerType_,nSnow_,&
    mLayerLiqFluxSnow_start,mLayerLiqFLuxSoil_start,&
    scalarAquiferTranspire,scalarAquiferRecharge,scalarAquiferBaseflow)
    implicit none
    integer(i4b),value :: nGRU
    integer(i4b),intent(in) :: ixCasNrg_(:), ixVegNrg_(:), ixVegHyd_(:)
    integer(i4b),intent(in) :: nLayers_(:)
    real(rkind),intent(inout) :: flux_data(:,:)
    integer(i4b),intent(in),value :: scalarCanairNetNrgFlux, scalarCanopyNetNrgFlux,scalarCanopyNetLiqFlux
    real(rkind),intent(in) :: canopyDepth_(:)
    integer(i4b),intent(in) :: ixSnowSoilNrg_(:,:)
    integer(i4b),value :: mLayerNrgFlux_start
    real(rkind),intent(inout) :: fluxVec_(:,:)
    integer(i4b),intent(in) :: ixAqWat_(:), ixSnowSoilHyd_(:,:), layerType_(:,:),nSnow_(:)
    integer(i4b),intent(in),value :: mLayerLiqFluxSnow_start, mLayerLiqFLuxSoil_start
    integer(i4b),intent(in),value :: scalarAquiferTranspire, scalarAquiferRecharge,scalarAquiferBaseflow
    integer(i4b) :: iGRU,iLayer
    iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

    if (iGRU .gt. nGRU) return

     ! *** WRAP UP ***
   ! define model flux vector for the vegetation sub-domain
   if (ixCasNrg_(iGRU)/=integerMissing) fluxVec_(ixCasNrg_(iGRU),iGRU) = flux_data(scalarCanairNetNrgFlux,iGRU)/canopyDepth_(iGRU)
   if (ixVegNrg_(iGRU)/=integerMissing) fluxVec_(ixVegNrg_(iGRU),iGRU) = flux_data(scalarCanopyNetNrgFlux,iGRU)/canopyDepth_(iGRU)
   if (ixVegHyd_(iGRU)/=integerMissing) fluxVec_(ixVegHyd_(iGRU),iGRU) = flux_data(scalarCanopyNetLiqFlux,iGRU)   ! NOTE: solid fluxes are handled separately
  !  if (nSnowSoilNrg>0) then ! if necessary, populate the flux vector for energy
     do iLayer=1,nLayers_(iGRU)
      if (ixSnowSoilNrg_(iLayer,iGRU)/=integerMissing) then   ! loop through non-missing energy state variables in the snow+soil domain
       fluxVec_( ixSnowSoilNrg_(iLayer,iGRU),iGRU ) = flux_data(mLayerNrgFlux_start+iLayer-1,iGRU)
      end if
     end do
  !  end if

      ! populate the flux vector for hydrology
   ! NOTE: ixVolFracWat  and ixVolFracLiq can also include states in the soil domain, hence enable primary variable switching
  !  if (nSnowSoilHyd>0) then  ! check if any hydrology states exist
     do iLayer=1,nLayers_(iGRU)     ! loop through non-missing energy state variables in the snow+soil domain
       if (ixSnowSoilHyd_(iLayer,iGRU)/=integerMissing) then   ! check if a given hydrology state exists
         select case(layerType_(iLayer,iGRU))
           case(iname_snow); fluxVec_(ixSnowSoilHyd_(iLayer,iGRU),iGRU) = flux_data(mLayerLiqFluxSnow_start+iLayer-1,iGRU)
           case(iname_soil); fluxVec_(ixSnowSoilHyd_(iLayer,iGRU),iGRU) = flux_data(mLayerLiqFluxSoil_start+iLayer-1-nSnow_(iGRU),iGRU)
          !  case default; err=20; message=trim(message)//'expect layerType to be either iname_snow or iname_soil'; return
         end select
       end if  ! end if a given hydrology state exists
     end do
  !  end if  ! end if any hydrology states exist
   ! compute the flux vector for the aquifer
   if (ixAqWat_(iGRU)/=integerMissing) fluxVec_(ixAqWat_(iGRU),iGRU) = flux_data(scalarAquiferTranspire,iGRU) + flux_data(scalarAquiferRecharge,iGRU) - flux_data(scalarAquiferBaseflow,iGRU)
    end subroutine

! **********************************************************************************************************
! public subroutine soilCmpres: compute soil compressibility (-) and its derivative w.r.t matric head (m-1)
! **********************************************************************************************************
subroutine soilCmpres(&
                      ! input:
                      dt,                                 & ! intent(in):  length of the time step (seconds)
                      ixRichards,                         & ! intent(in):  choice of option for Richards' equation
                      ixBeg,ixEnd,                        & ! intent(in):  start and end indices defining desired layers
                      mLayerMatricHead,                   & ! intent(in):  matric head at the start of the time step (m)
                      mLayerMatricHeadTrial,              & ! intent(in):  trial value of matric head (m)
                      mLayerVolFracLiqTrial,              & ! intent(in):  trial value for the volumetric liquid water content in each soil layer (-)
                      mLayerVolFracIceTrial,              & ! intent(in):  trial value for the volumetric ice content in each soil layer (-)
                      specificStorage,                    & ! intent(in):  specific storage coefficient (m-1)
                      theta_sat,                          & ! intent(in):  soil porosity (-)
                      ! output:
                      compress,                           & ! intent(out): compressibility of the soil matrix (-), per second
                      dCompress_dPsi,                     & ! intent(out): derivative in compressibility w.r.t. matric head (m-1)
                      err,message)                          ! intent(out): error code and error message
  implicit none
  ! input:
  real(rkind),intent(in)         :: dt                        !  length of the time step (seconds)
  integer(i4b),intent(in)        :: ixRichards                ! choice of option for Richards' equation
  integer(i4b),intent(in)        :: ixBeg,ixEnd               ! start and end indices defining desired layers
  real(rkind),intent(in)         :: mLayerMatricHead(:)       ! matric head at the start of the time step (m)
  real(rkind),intent(in)         :: mLayerMatricHeadTrial(:)  ! trial value for matric head (m)
  real(rkind),intent(in)         :: mLayerVolFracLiqTrial(:)  ! trial value for volumetric fraction of liquid water (-)
  real(rkind),intent(in)         :: mLayerVolFracIceTrial(:)  ! trial value for volumetric fraction of ice (-)
  real(rkind),intent(in)         :: specificStorage           ! specific storage coefficient (m-1)
  real(rkind),intent(in)         :: theta_sat(:)              ! soil porosity (-)
  ! output:
  real(rkind),intent(inout)      :: compress(:)               ! soil compressibility (-)
  real(rkind),intent(inout)      :: dCompress_dPsi(:)         ! derivative in soil compressibility w.r.t. matric head (m-1)
  integer(i4b),intent(out)       :: err                       ! error code
  character(*),intent(out)       :: message                   ! error message
  ! local variables
  integer(i4b)                   :: iLayer                    ! index of soil layer
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='soilCmpres/'
  ! (only compute for the mixed form of Richards' equation)
  if (ixRichards==mixdform) then
    do iLayer=1,size(mLayerMatricHead)
      if (iLayer>=ixBeg .and. iLayer<=ixEnd) then
      ! compute the derivative for the compressibility term (m-1), no volume expansion for total water
      dCompress_dPsi(iLayer) = specificStorage*(mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer))/theta_sat(iLayer)
      ! compute the compressibility term (-) per second
      compress(iLayer) = (mLayerMatricHeadTrial(iLayer) - mLayerMatricHead(iLayer))*dCompress_dPsi(iLayer)/dt
      end if
    end do
  else
    compress(:)       = 0._rkind
    dCompress_dPsi(:) = 0._rkind
  end if
end subroutine soilCmpres

! **********************************************************************************************************
! public subroutine soilCmpres: compute soil compressibility (-) and its derivative w.r.t matric head (m-1)
! **********************************************************************************************************

attributes(global) subroutine soilCmpresPrime_kernel(nGRU,ixRichards,&
nSoil,nSnow,nLayers,&
mLayerMatricHeadPrime,mLayerVolFracLiqTrial,mLayerVolFracIceTrial,&
specificStorage,theta_sat,&
mLayerCompress,dCompress_dPsi,&
scalarSoilCompress,&
mLayerDepth)
  integer(i4b),value,intent(in) :: nGRU,nSoil
  integer(i4b),intent(in) :: ixRichards
  integer(i4b),intent(in) :: nSnow(:),nLayers(:)
  real(rkind),intent(in) :: mLayerMatricHeadPrime(:,:)
  real(rkind),intent(in) :: mLayerVolFracLiqTrial(:,:)
  real(rkind),intent(in) :: mLayerVolFracIceTrial(:,:)
  real(rkind),intent(in) :: specificStorage(:)
  real(rkind),intent(in) :: theta_sat(:,:)
  real(rkind),intent(inout) :: mLayerCompress(:,:)
  real(rkind),intent(inout) :: dCompress_dPsi(:,:)
  real(rkind),intent(inout) :: scalarSoilCompress(:)
  real(rkind),intent(in) :: mLayerDepth(:,:)
  integer(i4b) :: err
    integer(i4b) :: iGRU
    iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

    if (iGRU .gt. nGRU) return

    call soilCmpresPrime(ixRichards,1,nSoil,&
    mLayerMatricHeadPrime(1:nSoil,iGRU),&
    mLayerVolFracLiqTrial(nSnow(iGRU)+1:nLayers(iGRU),iGRU),&
                        mLayerVolFracIceTrial(nSnow(iGRU)+1:nLayers(iGRU),iGRU), & ! intent(in):    trial value for the volumetric ice content in each soil layer (-)
                    specificStorage(iGRU),                        & ! intent(in):    specific storage coefficient (m-1)
                    theta_sat(:,iGRU),                              & ! intent(in):    soil porosity (-)
                    ! output:
                    mLayerCompress(:,iGRU),                         & ! intent(inout): compressibility of the soil matrix (-)
                    dCompress_dPsi(:,iGRU),                         & ! intent(inout): derivative in compressibility w.r.t. matric head (m-1)
                    err)

  ! compute the total change in storage associated with compression of the soil matrix (kg m-2 s-1)
  scalarSoilCompress(iGRU) = sum(mLayerCompress(1:nSoil,iGRU)*mLayerDepth(nSnow(iGRU)+1:nLayers(iGRU),iGRU))*iden_water

end subroutine


attributes(device) subroutine soilCmpresPrime(&
                          ! input:
                          ixRichards,                         & ! intent(in):  choice of option for Richards' equation
                          ixBeg,ixEnd,                        & ! intent(in):  start and end indices defining desired layers
                          mLayerMatricHeadPrime,              & ! intent(in):  matric head at the start of the time step (m)
                          mLayerVolFracLiqTrial,              & ! intent(in):  trial value for the volumetric liquid water content in each soil layer (-)
                          mLayerVolFracIceTrial,              & ! intent(in):  trial value for the volumetric ice content in each soil layer (-)
                          specificStorage,                    & ! intent(in):  specific storage coefficient (m-1)
                          theta_sat,                          & ! intent(in):  soil porosity (-)
                          ! output:
                          compress,                           & ! intent(out): compressibility of the soil matrix (-)
                          dCompress_dPsi,                     & ! intent(out): derivative in compressibility w.r.t. matric head (m-1)
                          err)                          ! intent(out): error code and error message
  implicit none
  ! input:
  integer(i4b),intent(in)           :: ixRichards               ! choice of option for Richards' equation
  integer(i4b),intent(in)           :: ixBeg,ixEnd              ! start and end indices defining desired layers
  real(rkind),intent(in)            :: mLayerMatricHeadPrime(:) ! matric head at the start of the time step (m)
  real(rkind),intent(in)            :: mLayerVolFracLiqTrial(:) ! trial value for volumetric fraction of liquid water (-)
  real(rkind),intent(in)            :: mLayerVolFracIceTrial(:) ! trial value for volumetric fraction of ice (-)
  real(rkind),intent(in)            :: specificStorage          ! specific storage coefficient (m-1)
  real(rkind),intent(in)            :: theta_sat(:)             ! soil porosity (-)
  ! output:
  real(rkind),intent(inout)         :: compress(:)              ! soil compressibility (-)
  real(rkind),intent(inout)         :: dCompress_dPsi(:)        ! derivative in soil compressibility w.r.t. matric head (m-1)
  integer(i4b),intent(out)          :: err                      ! error code
  ! local variables
  integer(i4b)                      :: iLayer                   ! index of soil layer
  ! --------------------------------------------------------------
  ! initialize error control
  err=0
  ! (only compute for the mixed form of Richards' equation)
  if (ixRichards==mixdform) then
    do iLayer=1,size(mLayerMatricHeadPrime)
      if (iLayer>=ixBeg .and. iLayer<=ixEnd) then
          ! compute the derivative for the compressibility term (m-1), no volume expansion for total water
          dCompress_dPsi(iLayer) = specificStorage*(mLayerVolFracLiqTrial(iLayer) + mLayerVolFracIceTrial(iLayer))/theta_sat(iLayer)
          ! compute the compressibility term (-) instantaneously
          compress(iLayer) = mLayerMatricHeadPrime(iLayer) * dCompress_dPsi(iLayer)
      end if
    end do
  else
    compress(:)       = 0._rkind
    dCompress_dPsi(:) = 0._rkind
  end if
end subroutine soilCmpresPrime

end module computFlux_module
