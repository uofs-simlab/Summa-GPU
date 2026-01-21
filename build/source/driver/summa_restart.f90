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
USE globalData,only:realMissing      ! missing real number

! named variables
USE var_lookup,only:iLookPROG                               ! look-up values for local column model prognostic (state) variables
USE var_lookup,only:iLookDIAG                               ! look-up values for local column model diagnostic variables
USE var_lookup,only:iLookFLUX                               ! look-up values for local column model fluxes
USE var_lookup,only:iLookBVAR                               ! look-up values for basin-average model variables
USE var_lookup,only:iLookDECISIONS                          ! look-up values for model decisions

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
 USE var_derive_module,only:calcHeight                       ! module to calculate height at layer interfaces and layer mid-point
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
   use device_data_types
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
  type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(summa1_struc%nGRU/128+1,1,1)

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
  indxStruct           => summa1_struc%indxStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model indices
  mparStruct           => summa1_struc%mparStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model parameters
  progStruct           => summa1_struc%progStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
  diagStruct           => summa1_struc%diagStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
  fluxStruct           => summa1_struc%fluxStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
  ! basin-average structures
  bparStruct           => summa1_struc%bparStruct          , & ! x%gru(:)%var(:)            -- basin-average parameters
  bvarStruct           => summa1_struc%bvarStruct_d          , & ! x%gru(:)%var(:)%dat        -- basin-average variables
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


 call readRestart_kernel<<<1,1>>>(nGRU,summa1_struc%decisions%aquiferIni,summa1_struc%decisions%spatial_gw,summa1_struc%decisions%hc_profile,summa1_struc%decisions%rootProfil,summa1_struc%decisions%groundwatr,summa1_struc%decisions%nrgConserv,&
no_icond_enth, &
bvarStruct%basin__AquiferStorage,progStruct%scalarAquiferStorage,&
fluxStruct%ixScalarCanopyLiqDrainage,fluxStruct%data,&
mparStruct%vGn_n_,diagStruct%scalarVGn_m,&
fluxStruct%ixmLayerSatHydCond_start,fluxStruct%ixmLayerSatHydCond_end,fluxStruct%ixmLayerSatHydCondMP_start,fluxStruct%ixmLayerSatHydCondMP_end,fluxStruct%ixiLayerSatHydCond_start,fluxStruct%ixiLayerSatHydCond_end,&
indxStruct%nSnow,indxStruct%nSoil,indxStruct%nLayers_d,&
mparStruct%k_soil_,mparStruct%k_macropore_,mparStruct%compactedDepth_,mparStruct%zScale_TOPMODEL_,&
progStruct%mLayerHeight,progStruct%iLayerHeight,progStruct%mLayerDepth,&
indxStruct%layerType,&
mparStruct%rootScaleFactor1_,mparStruct%rootScaleFactor2_,mparStruct%rootingDepth_,mparStruct%rootDistExp_,&
diagStruct%scalarAquiferRootFrac,diagStruct%mLayerRootDensity,&
 progStruct%scalarSWE,progStruct%mLayerVolFracLiq,progStruct%mLayerVolFracIce,&
 mparStruct%vGn_alpha_,mparStruct%theta_res_,mparStruct%theta_sat_,&
 progStruct%mLayerTemp,diagStruct%mLayerEnthTemp,progStruct%mLayerEnthalpy,&
 progStruct%mLayerMatricHead,mparStruct%soil_dens_intr_,mparStruct%snowfrz_scale_,&
  progStruct%scalarCanairTemp,progStruct%scalarCanairEnthalpy,&
  mparStruct%heightCanopyTop_,mparStruct%heightCanopyBottom_,&
  mparStruct%specificHeatVeg_,mparStruct%maxMassVegetation_,&
  progStruct%scalarCanopyTemp,diagStruct%scalarCanopyEnthTemp,progStruct%scalarCanopyEnthalpy,progStruct%scalarCanopyIce,progStruct%scalarCanopyLiq,&
  progStruct%scalarSnowAlbedo,mparStruct%albedoMax_,mparStruct%albedoMinWinter_,&
  progStruct%spectralSnowAlbedoDiffuse,mparStruct%albedoMaxVisible_,&
  mparStruct%albedoMinVisible_,mparStruct%albedoMaxNearIR_,mparStruct%albedoMinNearIR_)
 ! *****************************************************************************
  ! *** initialize time step
  ! *****************************************************************************
 ! loop through GRUs
 do iGRU=1,nGRU

  ! initialize time step length for each HRU
  do iHRU=1,gru_struc(iGRU)%hruCount
   dt_init%gru(iGRU)%hru(iHRU) = progStruct%dt_init ! seconds
  end do
end do
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

  attributes(global) subroutine readRestart_kernel(nGRU,aquiferIni,spatial_gw,hc_profile,rootProfil,groundwatr,nrgConserv, &
 no_icond_enth, &
 basin__AquiferStorage_,scalarAquiferStorage_, scalarCanopyLiqDrainage_, flux_data,&
 vGn_n,scalarVGn_m_,&
 mLayerSatHydCond_start,mLayerSatHydCond_end,mLayerSatHydCondMP_start,mLayerSatHydCondMP_end,iLayerSatHydCond_start,iLayerSatHydCond_end,&
 nSnow_,nSoil,nLayers_,&
 k_soil,k_macropore,compactedDepth,zScale_TOPMODEL,&
 mLayerHeight_,iLayerHeight_,mLayerDepth_,&
 layerType_,&
 rootScaleFactor1,rootScaleFactor2,rootingDepth,rootDistExp,&
 scalarAquiferRootFrac_,mLayerRootDensity_,&
 scalarSWE_,mLayerVolFracLiq_,mLayerVolFracIce_,&
 vGn_alpha,theta_res,theta_sat,&
 mLayerTemp_,mLayerEnthTemp_,mLayerEnthalpy_,&
 mLayerMatricHead_,soil_dens_intr_,snowfrz_scale_,&
  scalarCanairTemp_,scalarCanairEnthalpy_,&
  heightCanopyTop_,heightCanopyBottom_,&
  specificHeatVeg_,maxMassVegetation_,&
  scalarCanopyTemp_,scalarCanopyEnthTemp_,scalarCanopyEnthalpy_,scalarCanopyIce_,scalarCanopyLiq_,&
  scalarSnowAlbedo_,albedoMax_,albedoMinWinter_,&
  spectralSnowAlbedoDiffuse,albedoMaxVisible_,&
  albedoMinVisible_,albedoMaxNearIR_,albedoMinNearIR_)
 use nrtype
  USE mDecisions_module,only:&
   fullStart,      & ! start with full aquifer
   emptyStart        ! start with empty aquifer
 USE mDecisions_module,only:&                                ! look-up values for the choice of method for the spatial representation of groundwater
   localColumn,    & ! separate groundwater representation in each local soil column
   singleBasin       ! single groundwater store over the entire basin
 USE var_derive_module,only:v_shortcut                       ! module to calculate "short-cut" variables
 USE var_derive_module,only:satHydCond                       ! module to calculate the saturated hydraulic conductivity in each soil layer
 USE var_derive_module,only:rootDensty                       ! module to calculate the vertical distribution of roots
 USE var_derive_module,only:calcHeight_d                       ! module to calculate height at layer interfaces and layer mid-point
 USE check_icond_module,only:check_icond                     ! module to check initial conditions
 USE mDecisions_module,only:&
   closedForm,     & ! use temperature with closed form heat capacity
   enthalpyFormLU, & ! use enthalpy with soil temperature-enthalpy lookup tables
   enthalpyForm      ! use enthalpy with soil temperature-enthalpy analytical solution

 implicit none
  integer(i4b),intent(in),value :: nGRU
  integer(i4b) :: aquiferIni,spatial_gw,hc_profile,rootProfil,groundwatr,nrgConserv
  logical(lgt),intent(in),value :: no_icond_enth
  real(rkind) :: basin__AquiferStorage_(:), scalarAquiferStorage_(:)
  integer(i4b),intent(in),value :: scalarCanopyLiqDrainage_
  real(rkind) :: flux_data(:,:)
  real(rkind) :: vGn_n(:,:), scalarVGn_m_(:,:)
  integer(i4b),intent(in),value :: mLayerSatHydCond_start,mLayerSatHydCond_end,mLayerSatHydCondMP_start,mLayerSatHydCondMP_end,iLayerSatHydCond_start,iLayerSatHydCond_end
  integer(i4b),intent(in) :: nSnow_(:),nLayers_(:)
  integer(i4b),intent(in),value :: nSoil
  real(rkind),intent(in) :: k_soil(:,:), k_macropore(:,:)
  real(rkind),intent(in) :: compactedDepth(:), zScale_TOPMODEL(:)
  real(rkind),intent(inout) :: mLayerHeight_(0:,:), iLayerHeight_(0:,:),mLayerDepth_(:,:)
  integer(i4b),intent(inout) :: layerType_(:,:)
  real(rkind),intent(in) :: rootScaleFactor1(:),rootScaleFactor2(:),rootingDepth(:),rootDistExp(:)
  real(rkind),intent(inout) :: scalarAquiferRootFrac_(:)
  real(rkind),intent(inout) :: mLayerRootDensity_(:,:)
  real(rkind),intent(inout) :: scalarSWE_(:)
  real(rkind),intent(inout) :: mLayerVolFracLiq_(:,:), mLayerVolFracIce_(:,:)
  real(rkind),intent(in) :: vGn_alpha(:,:), theta_res(:,:), theta_sat(:,:)
  real(rkind) :: mLayerTemp_(:,:), mLayerEnthTemp_(:,:),mLayerEnthalpy_(:,:)
  real(rkind) :: temperature_(:,:,:),psiLiq_int_(:,:,:),deriv2_(:,:,:)
  real(rkind) :: mLayerMatricHead_(:,:), soil_dens_intr_(:,:)
  real(rkind),intent(in) :: snowfrz_scale_(:)
  real(rkind) :: scalarCanairTemp_(:),scalarCanairEnthalpy_(:)
    real(rkind) :: heightCanopyTop_(:),heightCanopyBottom_(:)
    real(rkind) :: specificHeatVeg_(:),maxMassVegetation_(:)
    real(rkind) :: scalarCanopyTemp_(:),scalarCanopyEnthTemp_(:),scalarCanopyEnthalpy_(:),scalarCanopyIce_(:),scalarCanopyLiq_(:)
  real(rkind) :: scalarSnowAlbedo_(:), albedoMax_(:), albedoMinWinter_(:)
  real(rkind) :: spectralSnowAlbedoDiffuse(:,:)
  real(rkind) :: albedoMaxVisible_(:), albedoMinVisible_(:), albedoMaxNearIR_(:), albedoMinNearIR_(:)

  logical(lgt) :: checkEnthalpy,use_lookup
  real(rkind) :: aquifer_start
  integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return

  do iGRU=1,nGRU
  ! check initial conditions
 checkEnthalpy = .false.
 use_lookup    = .false.
 if(nrgConserv .ne. closedForm) checkEnthalpy = .true. ! check enthalpy either for mixed form energy equation or enthalpy state variable
 if(nrgConserv==enthalpyFormLU) use_lookup = .true.    ! use lookup tables for soil temperature-enthalpy instead of analytical solution
  call check_icond(checkEnthalpy,no_icond_enth,use_lookup,mLayerDepth_(:,iGRU),iLayerHeight_(:,iGRU),nLayers_(iGRU),nSnow_(iGRU),&
  scalarSWE_,mLayerVolFracLiq_(:,iGRU),mLayerVolFracIce_(:,iGRU),layerType_(:,iGRU),&
  vGn_n(:,iGRU),vGn_alpha(:,iGRU),theta_res(:,iGRU),theta_sat(:,iGRU),&
  mLayerTemp_(:,iGRU),mLayerEnthTemp_(:,iGRU),mLayerEnthalpy_(:,iGRU),&
  temperature_(:,:,iGRU),psiLiq_int_(:,:,iGRU),deriv2_(:,:,iGRU),&
  mLayerMatricHead_(:,iGRU),soil_dens_intr_(:,iGRU),snowfrz_scale_(iGRU),&
  scalarCanairTemp_,scalarCanairEnthalpy_,&
  heightCanopyTop_,heightCanopyBottom_,&
  specificHeatVeg_,maxMassVegetation_,&
  scalarCanopyTemp_,scalarCanopyEnthTemp_,scalarCanopyEnthalpy_,scalarCanopyIce_,scalarCanopyLiq_,&
  scalarSnowAlbedo_,albedoMax_,albedoMinWinter_,&
  spectralSnowAlbedoDiffuse(:,iGRU),albedoMaxVisible_,&
  albedoMinVisible_,albedoMaxNearIR_,albedoMinNearIR_)


  ! *****************************************************************************
  ! *** compute ancillary variables
  ! *****************************************************************************

   ! re-calculate height of each layer
   call calcHeight_d(nLayers_(iGRU),layerType_(:,iGRU), &
   mLayerDepth_(:,iGRU),mLayerHeight_(:,iGRU),iLayerHeight_(:,iGRU))

     ! calculate vertical distribution of root density
   call rootDensty(rootProfil,groundwatr,&
  rootScaleFactor1(iGRU),rootScaleFactor2(iGRU),rootingDepth(iGRU),rootDistExp(iGRU),&
  nSnow_(iGRU),nSoil,nLayers_(iGRU),iLayerHeight_(:,iGRU),&
  scalarAquiferRootFrac_(iGRU),mLayerRootDensity_(:,iGRU))

     ! calculate saturated hydraulic conductivity in each soil layer
   call satHydCond(hc_profile,k_soil(:,iGRU), k_macropore(:,iGRU), compactedDepth(iGRU),zScale_TOPMODEL(iGRU),&
  nSnow_(iGRU),nSoil,nLayers_(iGRU),&
  mLayerHeight_(:,iGRU),iLayerHeight_(:,iGRU),&
  flux_data(mLayerSatHydCondMP_start:mLayerSatHydCondMP_end,iGRU),flux_data(mLayerSatHydCond_start:mLayerSatHydCond_end,iGRU),flux_data(iLayerSatHydCond_start:iLayerSatHydCond_end,iGRU))


   ! calculate "short-cut" variables such as volumetric heat capacity
   call v_shortcut(vGn_n(:,iGRU),   & ! vector of model parameters
                   scalarVGn_m_(:,iGRU))    ! data structure of model diagnostic variables

   ! initialize canopy drip
   ! NOTE: canopy drip from the previous time step is used to compute throughfall for the current time step
  flux_data(scalarCanopyLiqDrainage_,iGRU) = 0._rkind  ! not used

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
  end select  ! aquifer option


    ! select groundwater option
  select case(spatial_gw)

   ! the basin-average aquifer storage is not used if the groundwater is included in the local column
   case(localColumn)
    basin__AquiferStorage_(iGRU) = 0._rkind ! set to zero to be clear that there is no basin-average aquifer storage in this configuration
    if(aquiferIni==emptyStart) scalarAquiferStorage_(iGRU) = aquifer_start ! leave at initialized values if fullStart

   ! the local column aquifer storage is not used if the groundwater is basin-average
   ! (i.e., where multiple HRUs drain to a basin-average aquifer)
   case(singleBasin)
    basin__AquiferStorage_(iGRU) = aquifer_start 
     scalarAquiferStorage_(iGRU) = 0._rkind  ! set to zero to be clear that there is no local aquifer storage in this configuration

  end select  ! groundwater option
end do
do iGRU=1,nGRU
    print*, iGRU, scalarCanopyTemp_(iGRU)
end do
end subroutine

end module summa_restart




