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

module summa_setup
! initializes parameter data structures (e.g. vegetation and soil parameters).

! access missing values
USE globalData,only:integerMissing      ! missing integer
USE globalData,only:realMissing         ! missing real number

! global data on the forcing file
USE globalData,only:data_step           ! length of the data step (s)

! named variables
USE var_lookup,only:iLookATTR           ! look-up values for local attributes
USE var_lookup,only:iLookTYPE           ! look-up values for classification of veg, soils etc.
USE var_lookup,only:iLookPARAM          ! look-up values for local column model parameters
USE var_lookup,only:iLookINDEX          ! look-up values for local column model indices
USE var_lookup,only:iLookLOOKUP         ! look-up values for local column lookup tables
USE var_lookup,only:iLookID             ! look-up values for local column model ids
USE var_lookup,only:iLookBVAR           ! look-up values for basin-average model variables
USE var_lookup,only:iLookDECISIONS      ! look-up values for model decisions
USE globalData,only:urbanVegCategory    ! vegetation category for urban areas

! metadata structures
USE globalData,only:mpar_meta,bpar_meta ! parameter metadata structures

! look-up values for the choice of variable in energy equations (BE residual or IDA state variable)
USE mDecisions_module,only:&
  closedForm,    &                      ! use temperature with closed form heat capacity
  enthalpyFormLU,&                      ! use enthalpy with soil temperature-enthalpy lookup tables
  enthalpyForm                          ! use enthalpy with soil temperature-enthalpy analytical solution

! named variables to define the decisions for snow layers
USE mDecisions_module,only:&
  sameRulesAllLayers,&                  ! SNTHERM option: same combination/sub-dividion rules applied to all layers
  rulesDependLayerIndex                 ! CLM option: combination/sub-dividion rules depend on layer index

! named variables to define LAI decisions
USE mDecisions_module,only:&
 monthlyTable,&                         ! LAI/SAI taken directly from a monthly table for different vegetation classes
 specified                              ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters

! safety: set private unless specified otherwise
implicit none
private
public::summa_paramSetup
contains

 ! initializes parameter data structures (e.g. vegetation and soil parameters).
 subroutine summa_paramSetup(summa1_struc, err, message)
 ! ---------------------------------------------------------------------------------------
 ! * desired modules
 ! ---------------------------------------------------------------------------------------
 USE nrtype                                                  ! variable types, etc.
 USE summa_type, only:summa1_type_dec                        ! master summa data type
 ! subroutines and functions
 use time_utils_module,only:elapsedSec                       ! calculate the elapsed time
 USE mDecisions_module,only:mDecisions                       ! module to read model decisions
 USE ffile_info_module,only:ffile_info                       ! module to read information on forcing datafile
 USE read_attrb_module,only:read_attrb                       ! module to read local attributes
 USE read_pinit_module,only:read_pinit                       ! module to read initial model parameter values
 USE paramCheck_module,only:paramCheck                       ! module to check consistency of model parameters
 USE pOverwrite_module,only:pOverwrite                       ! module to overwrite default parameter values with info from the Noah tables
 USE read_param_module,only:read_param                       ! module to read model parameter sets
 USE enthalpyTemp_module,only:T2H_lookup_snWat               ! module to calculate a look-up table for the snow temperature-enthalpy conversion
 USE enthalpyTemp_module,only:T2L_lookup_soil                ! module to calculate a look-up table for the soil temperature-enthalpy conversion
 USE var_derive_module,only:fracFuture                       ! module to calculate the fraction of runoff in future time steps (time delay histogram)
 USE module_sf_noahmplsm,only:read_mp_veg_parameters         ! module to read NOAH vegetation tables
 ! global data structures
 USE globalData,only:gru_struc                               ! gru-hru mapping structures
 USE globalData,only:localParFallback                        ! local column default parameters
 USE globalData,only:basinParFallback                        ! basin-average default parameters
 USE globalData,only:model_decisions                         ! model decision structure
 USE globalData,only:greenVegFrac_monthly                    ! fraction of green vegetation in each month (0-1)
 ! run time options
 USE globalData,only:startGRU                                ! index of the starting GRU for parallelization run
 USE globalData,only:checkHRU                                ! index of the HRU for a single HRU run
 USE globalData,only:iRunMode                                ! define the current running mode
 ! output constraints
 USE globalData,only:maxLayers                               ! maximum number of layers
 USE globalData,only:maxSnowLayers                           ! maximum number of snow layers
 ! timing variables
 USE globalData,only:startSetup,endSetup                     ! date/time for the start and end of the parameter setup
 USE globalData,only:elapsedSetup                            ! elapsed time for the parameter setup
 ! file paths
 USE summaFileManager,only:SETTINGS_PATH                     ! define path to settings files (e.g., parameters, soil and veg. tables)
 USE summaFileManager,only:LOCAL_ATTRIBUTES                  ! name of model initial attributes file
 USE summaFileManager,only:LOCALPARAM_INFO,BASINPARAM_INFO   ! files defining the default values and constraints for model parameters
 USE summaFileManager,only:GENPARM,VEGPARM,SOILPARM,MPTABLE  ! files defining the noah tables
 ! Noah-MP parameters
 USE NOAHMP_VEG_PARAMETERS,only:SAIM,LAIM                    ! 2-d tables for stem area index and leaf area index (vegType,month)
 USE NOAHMP_VEG_PARAMETERS,only:HVT,HVB                      ! height at the top and bottom of vegetation (vegType)
     USE flxMapping_module,only:flxMapping                       ! module to map fluxes to states
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
 character(len=256)                    :: cmessage           ! error message of downwind routine
 character(len=256)                    :: attrFile           ! attributes file name
 integer(i4b)                          :: jHRU,kHRU          ! HRU indices
 integer(i4b)                          :: iGRU,iHRU          ! looping variables
 integer(i4b)                          :: iVar               ! looping variables
 real(rkind)                           :: absEnergyFac       ! multiplier for absolute value of energy state variable (for enthalpy or temperature)
 logical                               :: needLookup_soil    ! logical to decide if computing soil enthalpy lookup tables
 
 type(attr_data_device) :: attrStruct_d
 type(bvar_data_device) :: bvarStruct_d
 real(rkind),device,allocatable :: greenVegFrac_monthly_d(:)
 integer(i4b) :: soilTypeIndex, vegTypeIndex
   type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(summa1_struc%nGRU/128+1,1,1)

 ! ---------------------------------------------------------------------------------------
 ! associate to elements in the data structure
 summaVars: associate(&

  ! primary data structures (scalars)
  attrStruct           => summa1_struc%attrStruct          , & ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
  typeStruct           => summa1_struc%typeStruct          , & ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
  idStruct             => summa1_struc%idStruct            , & ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU

  ! primary data structures (variable length vectors)
  mparStruct           => summa1_struc%mparStruct_d          , & ! x%gru(:)%hru(:)%var(:)%dat -- model parameters
  dparStruct           => summa1_struc%dparStruct          , & ! x%gru(:)%hru(:)%var(:)     -- default model parameters

  ! basin-average structures
  bparStruct           => summa1_struc%bparStruct          , & ! x%gru(:)%var(:)            -- basin-average parameters
  bvarStruct           => summa1_struc%bvarStruct          , & ! x%gru(:)%var(:)%dat        -- basin-average variables

  ! lookup table structure
  lookupStruct         => summa1_struc%lookupStruct        , & ! x%gru(:)%hru(:)%z(:)%var(:)%lookup    -- lookup-tables

  ! miscellaneous variables
  upArea               => summa1_struc%upArea              , & ! area upslope of each HRU
  nGRU                 => summa1_struc%nGRU                , & ! number of grouped response units
  nHRU                 => summa1_struc%nHRU                  & ! number of global hydrologic response units

 ) ! assignment to variables in the data structures
 ! ---------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='summa_paramSetup/'
 
 ! initialize the start of the initialization
 call date_and_time(values=startSetup)

#ifdef NGEN_FORCING_ACTIVE
 ! *****************************************************************************
 ! if using NGEN forcing only need to set the hourly data_step (fixed)
 ! *****************************************************************************
 data_step = 3600._rkind
#else
 ! *****************************************************************************
 ! *** read description of model forcing datafile used in each HRU
 ! *****************************************************************************
 call ffile_info(nGRU,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
#endif

 ! *****************************************************************************
 ! *** read model decisions
 ! *****************************************************************************
 ! NOTE: Must be after ffile_info because mDecisions uses the data_step
 call mDecisions(err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! decide if computing soil enthalpy lookup tables and vegetation enthalpy lookup tables
 needLookup_soil = .false.
 ! if need enthalpy for either energy backward Euler residual or IDA state variable and not using soil enthalpy hypergeometric function
 if(model_decisions(iLookDECISIONS%nrgConserv)%iDecision == enthalpyFormLU) needLookup_soil = .true. 
 ! if using IDA and enthalpy as a state variable, need temperature-enthalpy lookup tables for soil and vegetation
 
 ! get the maximum number of snow layers
 select case(model_decisions(iLookDECISIONS%snowLayers)%iDecision)
  case(sameRulesAllLayers);    maxSnowLayers = 100
  case(rulesDependLayerIndex); maxSnowLayers = 5
  case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
 end select ! (option to combine/sub-divide snow layers)

 ! get the maximum number of layers
 maxLayers = gru_struc(1)%hruInfo(1)%nSoil + maxSnowLayers
  ! define mapping between fluxes and states
  call flxMapping(err,cmessage, gru_struc(1)%hruInfo(1)%nSoil)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****************************************************************************
 ! *** read local attributes for each HRU
 ! *****************************************************************************

 ! define the attributes file
 attrFile = trim(SETTINGS_PATH)//trim(LOCAL_ATTRIBUTES)

 ! read local attributes for each HRU
 call read_attrb(trim(attrFile),nGRU,attrStruct,typeStruct,idStruct,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****************************************************************************
 ! *** read default model parameters
 ! *****************************************************************************

 ! read default values and constraints for model parameters (local column)
 select case(model_decisions(iLookDECISIONS%nrgConserv)%iDecision)
   case(closedForm) ! ida temperature state variable
     absEnergyFac = 1.e2_rkind ! energy state variable is 2 orders of magnitude larger than mass state variable
   case(enthalpyFormLU,enthalpyForm) ! ida enthalpy state variable
     absEnergyFac = 1.e7_rkind ! energy state variable is 7 orders of magnitude larger than mass state variable
   case default; err=20; message=trim(message)//'unable to identify option for energy conservation'; return
 end select ! (option for energy conservation)

 call read_pinit(LOCALPARAM_INFO,.TRUE., absEnergyFac,mpar_meta,localParFallback,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! read default values and constraints for model parameters (basin-average)
 call read_pinit(BASINPARAM_INFO,.FALSE.,absEnergyFac, bpar_meta,basinParFallback,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****************************************************************************
 ! *** read Noah vegetation and soil tables
 ! *****************************************************************************

 ! define monthly fraction of green vegetation
 greenVegFrac_monthly = (/0.01_rkind, 0.02_rkind, 0.03_rkind, 0.07_rkind, 0.50_rkind, 0.90_rkind, 0.95_rkind, 0.96_rkind, 0.65_rkind, 0.24_rkind, 0.11_rkind, 0.02_rkind/)

 ! read Noah soil and vegetation tables
 call soil_veg_gen_parm(trim(SETTINGS_PATH)//trim(VEGPARM),                            & ! filename for vegetation table
                        trim(SETTINGS_PATH)//trim(SOILPARM),                           & ! filename for soils table
                        trim(SETTINGS_PATH)//trim(GENPARM),                            & ! filename for general table
                        trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision),    & ! classification system used for vegetation
                        trim(model_decisions(iLookDECISIONS%soilCatTbl)%cDecision))      ! classification system used for soils

 ! read Noah-MP vegetation tables
 call read_mp_veg_parameters(trim(SETTINGS_PATH)//trim(MPTABLE),                       & ! filename for Noah-MP table
                             trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision)) ! classification system used for vegetation

 ! define urban vegetation category
 select case(trim(model_decisions(iLookDECISIONS%vegeParTbl)%cDecision))
  case('USGS');                     urbanVegCategory =    1
  case('MODIFIED_IGBP_MODIS_NOAH'); urbanVegCategory =   13
  case('plumberCABLE');             urbanVegCategory = -999
  case('plumberCHTESSEL');          urbanVegCategory = -999
  case('plumberSUMMA');             urbanVegCategory = -999
  case default
   message=trim(message)//'unable to identify vegetation category'
   return
 end select

 ! set default model parameters
 do iGRU=1,nGRU
  do iHRU=1,gru_struc(iGRU)%hruCount

   ! set parameters to their default value
   dparStruct%gru(iGRU)%hru(iHRU)%var(:) = localParFallback(:)%default_val         ! x%hru(:)%var(:)
   vegTypeIndex = typeStruct%vegTypeIndex(iGRU)
   soilTypeIndex = typeStruct%soilTypeIndex(iGRU)

   ! overwrite default model parameters with information from the Noah-MP tables
   call pOverwrite(vegTypeIndex,  &  ! vegetation category
                   soilTypeIndex, &  ! soil category
                   dparStruct%gru(iGRU)%hru(iHRU)%var,                          &  ! default model parameters
                   err,cmessage)                                                   ! error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! copy over to the parameter structure
   ! NOTE: constant for the dat(:) dimension (normally depth)
   do ivar=1,size(localParFallback)
    print*, iVar, iGRU, iHRU
    print*, mpar_meta(ivar)%varname
    ! print*, dparStruct%gru(iGRU)%hru(iHRU)%var(ivar)
    call set_device_param_data_scalar(mparStruct, iGRU, ivar, dparStruct%gru(iGRU)%hru(iHRU)%var(ivar))
    ! mparStruct%gru(iGRU)%hru(iHRU)%var(ivar)%dat(:) = dparStruct%gru(iGRU)%hru(iHRU)%var(ivar)
    print*, dparStruct%gru(iGRU)%hru(iHRU)%var(ivar)
   end do  ! looping through variables

  end do  ! looping through HRUs

  ! set default for basin-average parameters
  bparStruct%gru(iGRU)%var(:) = basinParFallback(:)%default_val

 end do  ! looping through GRUs

 ! *****************************************************************************
 ! *** read trial model parameter values for each HRU, and populate initial data structures
 ! *****************************************************************************
 call read_param(iRunMode,checkHRU,startGRU,nHRU,nGRU,idStruct,mparStruct,bparStruct,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****************************************************************************
 ! *** compute derived model variables that are pretty much constant for the basin as a whole
 ! *****************************************************************************
 ! loop through GRUs
 do iGRU=1,nGRU

  ! calculate the fraction of runoff in future time steps
  call fracFuture(bparStruct%gru(iGRU)%var,    &  ! vector of basin-average model parameters
                  bvarStruct%gru(iGRU),        &  ! data structure of basin-average variables
                  err,cmessage)                   ! error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute total area of the upstream HRUS that flow into each HRU
  do iHRU=1,gru_struc(iGRU)%hruCount
   upArea%gru(iGRU)%hru(iHRU) = 0._rkind
  end do  ! iHRU


 end do ! GRU
   call allocate_veg_parameters(summa1_struc%veg_param,summa1_struc%nGRU)

 call allocate_device_attr_data(attrStruct_d,attrStruct,nGRU)
 call allocate_device_bvar_data(bvarStruct_d,bvarStruct,nGRU)
!  call allocate_device_param_data(mparStruct_d,mparStruct,nGRU)

 greenVegFrac_monthly_d = greenVegFrac_monthly
 call paramSetup_kernel<<<blocks,threads>>>(nGRU,model_decisions(iLookDECISIONS%LAI_method)%iDecision,model_decisions(iLookDECISIONS%snowLayers)%iDecision,needLookup_soil,bvarStruct_d%basin__totalArea,attrStruct_d%HRUarea,&
typeStruct%vegTypeIndex,&
mparStruct%heightCanopyTop_,mparStruct%heightCanopyBottom_,mparStruct%winterSAI_,mparStruct%summerLAI_,greenVegFrac_monthly_d,&
summa1_struc%veg_param%HVT_,summa1_struc%veg_param%HVB_,summa1_struc%veg_param%SAIM_,summa1_struc%veg_param%LAIM_,&
mparStruct%snowfrz_scale_,mparStruct%theta_sat_,mparStruct%theta_res_,mparStruct%vGn_alpha_,mparStruct%vGn_n_,lookupStruct%temperature,lookupStruct%psiLiq_int,lookupStruct%deriv2,gru_struc(1)%hruInfo(1)%nSoil,&
lookupStruct%h_lookup,lookupStruct%t_lookup,&
mparStruct%critSoilTranspire_,mparStruct%critSoilWilting_,mparStruct%fieldCapacity_,mparStruct%zMax_,mparStruct%zMin_,&
 mparStruct%zminLayer1_,mparStruct%zminLayer2_,mparStruct%zminLayer3_,mparStruct%zminLayer4_,mparStruct%zminLayer5_,&
 mparStruct%zmaxLayer1_lower_,mparStruct%zmaxLayer2_lower_,mparStruct%zmaxLayer3_lower_,mparStruct%zmaxLayer4_lower_,&
 mparStruct%zmaxLayer1_upper_,mparStruct%zmaxLayer2_upper_,mparStruct%zmaxLayer3_upper_,mparStruct%zmaxLayer4_upper_)
 call finalize_device_bvar_data(bvarStruct_d,bvarStruct)
 ! identify the end of the initialization
 call date_and_time(values=endSetup)

 ! aggregate the elapsed time for the initialization
 elapsedSetup = elapsedSec(startSetup, endSetup)

 ! end associate statements
 end associate summaVars


 end subroutine summa_paramSetup

attributes(global) subroutine paramSetup_kernel(nGRU,LAI_method,snowLayers,needLookup_soil,totalArea_,HRUarea_,&
vegTypeIndex_,&
heightCanopyTop,heightCanopyBottom,winterSAI,summerLAI,greenVegFrac_monthly,&
HVT,HVB,SAIM,LAIM,&
snowfrz_scale_,theta_sat_,theta_res_,vGn_alpha_,vGn_n_,temperature,psiLiq_int,deriv2,nSoil,h_lookup,t_lookup,&
critSoilTranspire_,critSoilWilting_,fieldCapacity_,zMax,zMin,&
 zminLayer1,zminLayer2,zminLayer3,zminLayer4,zminLayer5,&
 zmaxLayer1_lower,zmaxLayer2_lower,zmaxLayer3_lower,zmaxLayer4_lower,&
 zmaxLayer1_upper,zmaxLayer2_upper,zmaxLayer3_upper,zmaxLayer4_upper)
 USE nrtype                                                  ! variable types, etc.
 USE enthalpyTemp_module,only:T2L_lookup_soil_device                ! module to calculate a look-up table for the soil temperature-enthalpy conversion
 USE enthalpyTemp_module,only:T2H_lookup_snWat               ! module to calculate a look-up table for the snow temperature-enthalpy conversion
 USE paramCheck_module,only:paramCheck_device                       ! module to check consistency of model parameters

 integer(i4b),intent(in),value :: nGRU,LAI_method,snowLayers
 logical(lgt),intent(in),value :: needLookup_soil
 real(rkind),intent(inout) :: totalArea_(:),HRUarea_(:)
 integer(i4b) :: vegTypeIndex_(:)
 real(rkind) :: heightCanopyTop(:),heightCanopyBottom(:),winterSAI(:),summerLAI(:)
 real(rkind) :: greenVegFrac_monthly(:)
 real(rkind) :: hvt(:,:), hvb(:,:)
 real(rkind) :: saim(:,:,:), laim(:,:,:)
 real(rkind) :: snowfrz_scale_(:),theta_sat_(:,:),theta_res_(:,:),vGn_alpha_(:,:),vGn_n_(:,:)
 real(rkind) :: temperature(:,:,:), psiLiq_int(:,:,:), deriv2(:,:,:)
 real(rkind) :: h_lookup(:), t_lookup(:)
 real(rkind) :: critSoilTranspire_(:),critSoilWilting_(:),fieldCapacity_(:)
 integer(i4b),intent(in),value :: nSoil
  real(rkind) :: zMax(:),zMin(:)
 real(rkind) :: zminLayer1(:),zminLayer2(:),zminLayer3(:),zminLayer4(:),zminLayer5(:)
 real(rkind) :: zmaxLayer1_lower(:),zmaxLayer2_lower(:),zmaxLayer3_lower(:),zmaxLayer4_lower(:)
 real(rkind) :: zmaxLayer1_upper(:),zmaxLayer2_upper(:),zmaxLayer3_upper(:),zmaxLayer4_upper(:)

  integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return


  call paramCheck_device(critSoilTranspire_(iGRU),critSoilWilting_(iGRU),theta_sat_(:,iGRU),theta_res_(:,iGRU),&
 heightCanopyTop(iGRU), heightCanopyBottom(iGRU),fieldCapacity_(iGRU),&
 snowLayers,zMax(iGRU),zMin(iGRU),&
 zminLayer1(iGRU),zminLayer2(iGRU),zminLayer3(iGRU),zminLayer4(iGRU),zminLayer5(iGRU),&
 zmaxLayer1_lower(iGRU),zmaxLayer2_lower(iGRU),zmaxLayer3_lower(iGRU),zmaxLayer4_lower(iGRU),&
 zmaxLayer1_upper(iGRU),zmaxLayer2_upper(iGRU),zmaxLayer3_upper(iGRU),zmaxLayer4_upper(iGRU))

  if (iGRU .eq. 1) then
       ! calculate a look-up table for the temperature-enthalpy conversion of snow for future snow layer merging
   ! NOTE1: might be able to make this more efficient by only doing this for the HRUs that have snow
   ! NOTE2: H is the mixture enthalpy of snow liquid and ice
   call T2H_lookup_snWat(snowfrz_scale_(iGRU),                     &  ! intent(in):    parameter data structure
                           h_lookup,t_lookup)
  end if

   ! calculate a lookup table for the temperature-enthalpy conversion of soil 
   ! NOTE: L is the integral of soil Clapeyron equation liquid water matric potential from temperature
   !       multiply by Cp_liq*iden_water to get temperature component of enthalpy
   if(needLookup_soil)then
    call T2L_lookup_soil_device(snowfrz_scale_(iGRU),theta_sat_(:,iGRU),theta_res_(:,iGRU),vGn_alpha_(:,iGRU),vGn_n_(:,iGRU),temperature(:,:,iGRU),psiLiq_int(:,:,iGRU),deriv2(:,:,iGRU),nSoil)
   endif

    ! overwrite the vegetation height
   HVT(vegTypeIndex_(iGRU),iGRU) = heightCanopyTop(iGRU)
   HVB(vegTypeIndex_(iGRU),iGRU) = heightCanopyBottom(iGRU)

   ! overwrite the tables for LAI and SAI
   if(LAI_method == specified)then
    SAIM(vegTypeIndex_(iGRU),:,iGRU) = winterSAI(iGRU)
    LAIM(vegTypeIndex_(iGRU),:,iGRU) = summerLAI(iGRU)*greenVegFrac_monthly
   endif


    ! identify the total basin area for a GRU (m2)
   totalArea_(iGRU) = HRUarea_(iGRU)


end subroutine

 ! **************************************************************************************************
 ! private subroutine SOIL_VEG_GEN_PARM: Read soil, vegetation and other model parameters (from NOAH)
 ! **************************************************************************************************
 !-----------------------------------------------------------------
 SUBROUTINE SOIL_VEG_GEN_PARM(FILENAME_VEGTABLE, FILENAME_SOILTABLE, FILENAME_GENERAL, MMINLU, MMINSL)
 !-----------------------------------------------------------------
  use module_sf_noahlsm, only : shdtbl, nrotbl, rstbl, rgltbl, &
       &                        hstbl, snuptbl, maxalb, laimintbl, &
       &                        bb, drysmc, f11, maxsmc, laimaxtbl, &
       &                        emissmintbl, emissmaxtbl, albedomintbl, &
       &                        albedomaxtbl, wltsmc, qtz, refsmc, &
       &                        z0mintbl, z0maxtbl, &
       &                        satpsi, satdk, satdw, &
       &                        theta_res, theta_sat, vGn_alpha, vGn_n, k_soil, &  ! MPC add van Genutchen parameters
       &                        fxexp_data, lvcoef_data, &
       &                        lutype, maxalb, &
       &                        slope_data, frzk_data, bare, cmcmax_data, &
       &                        cfactr_data, csoil_data, czil_data, &
       &                        refkdt_data, natural, refdk_data, &
       &                        rsmax_data, salp_data, sbeta_data, &
       &                        zbot_data, smhigh_data, smlow_data, &
       &                        lucats, topt_data, slcats, slpcats, sltype

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: FILENAME_VEGTABLE, FILENAME_SOILTABLE, FILENAME_GENERAL
  CHARACTER(LEN=*), INTENT(IN) :: MMINLU, MMINSL
  integer :: LUMATCH, IINDEX, LC, NUM_SLOPE
  integer :: ierr
  INTEGER , PARAMETER :: OPEN_OK = 0

  character*128 :: mess , message

  !-----SPECIFY VEGETATION RELATED CHARACTERISTICS :
  !             ALBBCK: SFC albedo (in percentage)
  !                 Z0: Roughness length (m)
  !             SHDFAC: Green vegetation fraction (in percentage)
  !  Note: The ALBEDO, Z0, and SHDFAC values read from the following table
  !          ALBEDO, amd Z0 are specified in LAND-USE TABLE; and SHDFAC is
  !          the monthly green vegetation data
  !             CMXTBL: MAX CNPY Capacity (m)
  !             NROTBL: Rooting depth (layer)
  !              RSMIN: Mimimum stomatal resistance (s m-1)
  !              RSMAX: Max. stomatal resistance (s m-1)
  !                RGL: Parameters used in radiation stress function
  !                 HS: Parameter used in vapor pressure deficit functio
  !               TOPT: Optimum transpiration air temperature. (K)
  !             CMCMAX: Maximum canopy water capacity
  !             CFACTR: Parameter used in the canopy inteception calculati
  !               SNUP: Threshold snow depth (in water equivalent m) that
  !                     implies 100% snow cover
  !                LAI: Leaf area index (dimensionless)
  !             MAXALB: Upper bound on maximum albedo over deep snow
  !
  !-----READ IN VEGETAION PROPERTIES FROM VEGPARM.TBL
  !

  OPEN(19, FILE=trim(FILENAME_VEGTABLE),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening VEGPARM.TBL'
     CALL wrf_error_fatal ( message )
  END IF

  LUMATCH=0

  FIND_LUTYPE : DO WHILE (LUMATCH == 0)
     READ (19,*,END=2002)
     READ (19,*,END=2002)LUTYPE
     READ (19,*)LUCATS,IINDEX

     IF(LUTYPE.EQ.MMINLU)THEN
        WRITE( mess , * ) 'LANDUSE TYPE = ' // TRIM ( LUTYPE ) // ' FOUND', LUCATS,' CATEGORIES'
        ! CALL wrf_message( mess )
        LUMATCH=1
     ELSE
        call wrf_message ( "Skipping over LUTYPE = " // TRIM ( LUTYPE ) )
        DO LC = 1, LUCATS+12
           read(19,*)
        ENDDO
     ENDIF
  ENDDO FIND_LUTYPE
  ! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
  IF ( SIZE(SHDTBL)       < LUCATS .OR. &
       SIZE(NROTBL)       < LUCATS .OR. &
       SIZE(RSTBL)        < LUCATS .OR. &
       SIZE(RGLTBL)       < LUCATS .OR. &
       SIZE(HSTBL)        < LUCATS .OR. &
       SIZE(SNUPTBL)      < LUCATS .OR. &
       SIZE(MAXALB)       < LUCATS .OR. &
       SIZE(LAIMINTBL)    < LUCATS .OR. &
       SIZE(LAIMAXTBL)    < LUCATS .OR. &
       SIZE(Z0MINTBL)     < LUCATS .OR. &
       SIZE(Z0MAXTBL)     < LUCATS .OR. &
       SIZE(ALBEDOMINTBL) < LUCATS .OR. &
       SIZE(ALBEDOMAXTBL) < LUCATS .OR. &
       SIZE(EMISSMINTBL ) < LUCATS .OR. &
       SIZE(EMISSMAXTBL ) < LUCATS ) THEN
     CALL wrf_error_fatal('Table sizes too small for value of LUCATS in module_sf_noahdrv.F, expand NLUS and MVT parameters to size of vegetation table and recompile')
  ENDIF

  IF(LUTYPE.EQ.MMINLU)THEN
     DO LC=1,LUCATS
        READ (19,*)IINDEX,SHDTBL(LC),                        &
             NROTBL(LC),RSTBL(LC),RGLTBL(LC),HSTBL(LC), &
             SNUPTBL(LC),MAXALB(LC), LAIMINTBL(LC),     &
             LAIMAXTBL(LC),EMISSMINTBL(LC),             &
             EMISSMAXTBL(LC), ALBEDOMINTBL(LC),         &
             ALBEDOMAXTBL(LC), Z0MINTBL(LC), Z0MAXTBL(LC)
     ENDDO

     READ (19,*)
     READ (19,*)TOPT_DATA
     READ (19,*)
     READ (19,*)CMCMAX_DATA
     READ (19,*)
     READ (19,*)CFACTR_DATA
     READ (19,*)
     READ (19,*)RSMAX_DATA
     READ (19,*)
     READ (19,*)BARE
     READ (19,*)
     READ (19,*)NATURAL
  ENDIF

2002 CONTINUE

  CLOSE (19)
  IF (LUMATCH == 0) then
     CALL wrf_error_fatal ("Land Use Dataset '"//MMINLU//"' not found in VEGPARM.TBL.")
  ENDIF

!
!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
!
  OPEN(19, FILE=trim(FILENAME_SOILTABLE),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening SOILPARM.TBL'
     CALL wrf_error_fatal ( message )
  END IF

  WRITE(mess,*) 'INPUT SOIL TEXTURE CLASSIFICATION = ', TRIM ( MMINSL )
  ! CALL wrf_message( mess )

  LUMATCH=0

  ! MPC add a new soil table
  FIND_soilTYPE : DO WHILE (LUMATCH == 0)
   READ (19,*)
   READ (19,*,END=2003)SLTYPE
   READ (19,*)SLCATS,IINDEX
   IF(SLTYPE.EQ.MMINSL)THEN
     WRITE( mess , * ) 'SOIL TEXTURE CLASSIFICATION = ', TRIM ( SLTYPE ) , ' FOUND', &
          SLCATS,' CATEGORIES'
     ! CALL wrf_message ( mess )
     LUMATCH=1
   ELSE
    call wrf_message ( "Skipping over SLTYPE = " // TRIM ( SLTYPE ) )
    DO LC = 1, SLCATS
     read(19,*)
    ENDDO
   ENDIF
  ENDDO FIND_soilTYPE
  ! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
  IF ( SIZE(BB    ) < SLCATS .OR. &
       SIZE(DRYSMC) < SLCATS .OR. &
       SIZE(F11   ) < SLCATS .OR. &
       SIZE(MAXSMC) < SLCATS .OR. &
       SIZE(REFSMC) < SLCATS .OR. &
       SIZE(SATPSI) < SLCATS .OR. &
       SIZE(SATDK ) < SLCATS .OR. &
       SIZE(SATDW ) < SLCATS .OR. &
       SIZE(WLTSMC) < SLCATS .OR. &
       SIZE(QTZ   ) < SLCATS  ) THEN
     CALL wrf_error_fatal('Table sizes too small for value of SLCATS in module_sf_noahdrv.F')
  ENDIF

  ! MPC add new soil table
  select case(trim(SLTYPE))
   case('STAS','STAS-RUC')  ! original soil tables
     DO LC=1,SLCATS
        READ (19,*) IINDEX,BB(LC),DRYSMC(LC),F11(LC),MAXSMC(LC),&
             REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
             WLTSMC(LC), QTZ(LC)
     ENDDO
   case('ROSETTA')          ! new soil table
     DO LC=1,SLCATS
        READ (19,*) IINDEX,&
             ! new soil parameters (from Rosetta)
             theta_res(LC), theta_sat(LC),        &
             vGn_alpha(LC), vGn_n(LC), k_soil(LC), &
             ! original soil parameters
             BB(LC),DRYSMC(LC),F11(LC),MAXSMC(LC),&
             REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
             WLTSMC(LC), QTZ(LC)
     ENDDO
   case default
     CALL wrf_message( 'SOIL TEXTURE IN INPUT FILE DOES NOT ' )
     CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
     CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
  end select

2003 CONTINUE

  CLOSE (19)

  IF(LUMATCH.EQ.0)THEN
     CALL wrf_message( 'SOIL TEXTURE IN INPUT FILE DOES NOT ' )
     CALL wrf_message( 'MATCH SOILPARM TABLE'                 )
     CALL wrf_error_fatal ( 'INCONSISTENT OR MISSING SOILPARM FILE' )
  ENDIF

!
!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL
!
  OPEN(19, FILE=trim(FILENAME_GENERAL),FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
  IF(ierr .NE. OPEN_OK ) THEN
     WRITE(message,FMT='(A)') &
          'module_sf_noahlsm.F: soil_veg_gen_parm: failure opening GENPARM.TBL'
     CALL wrf_error_fatal ( message )
  END IF

  READ (19,*)
  READ (19,*)
  READ (19,*) NUM_SLOPE

  SLPCATS=NUM_SLOPE
! prevent possible array overwrite, Bill Bovermann, IBM, May 6, 2008
  IF ( SIZE(slope_data) < NUM_SLOPE ) THEN
     CALL wrf_error_fatal('NUM_SLOPE too large for slope_data array in module_sf_noahdrv')
  ENDIF

  DO LC=1,SLPCATS
     READ (19,*)SLOPE_DATA(LC)
  ENDDO

  READ (19,*)
  READ (19,*)SBETA_DATA
  READ (19,*)
  READ (19,*)FXEXP_DATA
  READ (19,*)
  READ (19,*)CSOIL_DATA
  READ (19,*)
  READ (19,*)SALP_DATA
  READ (19,*)
  READ (19,*)REFDK_DATA
  READ (19,*)
  READ (19,*)REFKDT_DATA
  READ (19,*)
  READ (19,*)FRZK_DATA
  READ (19,*)
  READ (19,*)ZBOT_DATA
  READ (19,*)
  READ (19,*)CZIL_DATA
  READ (19,*)
  READ (19,*)SMLOW_DATA
  READ (19,*)
  READ (19,*)SMHIGH_DATA
  READ (19,*)
  READ (19,*)LVCOEF_DATA
  CLOSE (19)

!-----------------------------------------------------------------
END SUBROUTINE SOIL_VEG_GEN_PARM
!-----------------------------------------------------------------

end module summa_setup
