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

module layerDivide_module

! variable types
USE nrtype

! physical constants
USE multiconst,only:&
                    iden_ice,       & ! intrinsic density of ice             (kg m-3)
                    iden_water        ! intrinsic density of liquid water    (kg m-3)

! access named variables for snow and soil
USE globalData,only:iname_snow        ! named variables for snow
USE globalData,only:iname_soil        ! named variables for soil

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! access metadata
USE globalData,only:prog_meta,diag_meta,flux_meta,indx_meta   ! metadata

! access the derived types to define the data structures
USE data_types,only:&
                    var_d,            & ! data vector (rkind)
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength,      & ! data vector with variable length dimension (rkind)
                    model_options       ! defines the model decisions

! access named variables defining elements in the data structures
USE var_lookup,only:iLookPROG,iLookDIAG,iLookFLUX,iLookINDEX  ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS                            ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM                                ! named variables for elements of the parameter structure

! define look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers,       & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex       ! CLM option: combination/sub-dividion rules depend on layer index

! define look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:&
 noah_mp,                  & ! full Noah-MP implementation (including albedo)
 CLM_2stream,              & ! CLM 2-stream model (see CLM documentation)
 UEB_2stream,              & ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
 NL_scatter,               & ! Simplified method Nijssen and Lettenmaier (JGR 1999)
 BeersLaw                    ! Beer's Law (as implemented in VIC)

! define look-up values for the choice of albedo method
USE mDecisions_module,only:& ! identify model options for snow albedo
 constantDecay,            & ! constant decay in snow albedo (e.g., VIC, CLASS)
 variableDecay               ! variable decay in snow albedo (e.g., BATS approach, with destructive metamorphism + soot content)

! privacy
implicit none
private
public::layerDivide,layerDivide_device

contains


  ! ***********************************************************************************************************
 ! public subroutine layerDivide: add new snowfall to the system, and increase number of snow layers if needed
 ! ***********************************************************************************************************
 subroutine layerDivide(&
  nGRU, &
  tooMuchMelt, &
                        ! input/output: model data structures
                        model_decisions,                 & ! intent(in):    model decisions
                        mpar_data,                       & ! intent(in):    model parameters
                        indx_data,                       & ! intent(inout): type of each layer
                        prog_data,                       & ! intent(inout): model prognostic variables for a local HRU
                        diag_data,                       & ! intent(inout): model diagnostic variables for a local HRU
                        flux_data,                       & ! intent(inout): model fluxes for a local HRU
                        ! output
                        divideLayer,                     & ! intent(out): flag to denote that a layer was divided
                        err,message)                       ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 ! computational modules
 USE snow_utils_module,only:fracliquid,templiquid          ! functions to compute temperature/liquid water
 use device_data_types
 USE globalData,only:maxSnowLayers, &                      ! maximum number of snow layers
                     veryBig
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input/output: model data structures
 integer(i4b),intent(in) :: nGRU
 logical(lgt),intent(in) :: tooMuchMelt
 type(decisions_device),intent(in)  :: model_decisions     ! model decisions
 type(mpar_data_device),intent(in)    :: mpar_data              ! model parameters
 type(indx_data_device),intent(inout) :: indx_data              ! type of each layer
 type(prog_data_device),intent(inout) :: prog_data              ! model prognostic variables for a local HRU
 type(diag_data_device),intent(inout) :: diag_data              ! model diagnostic variables for a local HRU
 type(flux_data_device),intent(inout) :: flux_data              ! model flux variables
 ! output
 logical(lgt),intent(out),device        :: divideLayer(:)            ! flag to denote that a layer was divided
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! define local variables
 character(LEN=256)              :: cmessage               ! error message of downwind routine
!  integer(i4b)                    :: nSnow                  ! number of snow layers
!  integer(i4b)                    :: nSoil                  ! number of soil layers
!  integer(i4b)                    :: nLayers                ! total number of layers
 integer(i4b)                    :: iLayer                 ! layer index
 integer(i4b)                    :: jLayer                 ! layer index
!  real(rkind),dimension(4)        :: zmax_lower             ! lower value of maximum layer depth
!  real(rkind),dimension(4)        :: zmax_upper             ! upper value of maximum layer depth
 real(rkind)                     :: zmaxCheck              ! value of zmax for a given snow layer
 integer(i4b)                    :: nCheck                 ! number of layers to check to divide
 logical(lgt)                    :: createLayer            ! flag to indicate we are creating a new snow layer
 real(rkind)                     :: depthOriginal          ! original layer depth before sub-division (m)
 real(rkind),parameter           :: fracTop=0.5_rkind      ! fraction of old layer used for the top layer
 real(rkind)                     :: surfaceLayerSoilTemp   ! temperature of the top soil layer (K)
 real(rkind)                     :: maxFrozenSnowTemp      ! maximum temperature when effectively all water is frozen (K)
 real(rkind),parameter           :: unfrozenLiq=0.01_rkind ! unfrozen liquid water used to compute maxFrozenSnowTemp (-)
 real(rkind)                     :: volFracWater           ! volumetric fraction of total water, liquid and ice (-)
 real(rkind)                     :: fracLiq                ! fraction of liquid water (-)
 integer(i4b),parameter          :: ixVisible=1            ! named variable to define index in array of visible part of the spectrum
 integer(i4b),parameter          :: ixNearIR=2             ! named variable to define index in array of near IR part of the spectrum
 real(rkind),parameter           :: snowDepthTol=1.e-10_rkind ! tolerance for the snow depth difference (m)
   type(dim3) :: blocks,threads
  threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)
  

 ! --------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="layerDivide/"

 ! --------------------------------------------------------------------------------------------------------
 ! associate variables in the data structures
 associate(&
 ! model decisions
 ix_snowLayers          => model_decisions%snowLayers, & ! decision for snow combination
 ! model parameters (compute layer temperature)
 fc_param               => mpar_data%snowfrz_scale_,       & ! freezing curve parameter for snow (K-1)
 ! model parameters (control the depth of snow layers)
 zmax                   => mpar_data%zmax_,                & ! maximum layer depth (m)
 zminLayer1             => mpar_data%zminLayer1_,          & ! minimum layer depth for the 1st (top) layer (m)
 zmaxLayer1_lower       => mpar_data%zmaxLayer1_lower_,    & ! maximum layer depth for the 1st (top) layer when only 1 layer (m)
 zmaxLayer2_lower       => mpar_data%zmaxLayer2_lower_,    & ! maximum layer depth for the 2nd layer when only 2 layers (m)
 zmaxLayer3_lower       => mpar_data%zmaxLayer3_lower_,    & ! maximum layer depth for the 3rd layer when only 3 layers (m)
 zmaxLayer4_lower       => mpar_data%zmaxLayer4_lower_,    & ! maximum layer depth for the 4th layer when only 4 layers (m)
 zmaxLayer1_upper       => mpar_data%zmaxLayer1_upper_,    & ! maximum layer depth for the 1st (top) layer when > 1 layer (m)
 zmaxLayer2_upper       => mpar_data%zmaxLayer2_upper_,    & ! maximum layer depth for the 2nd layer when > 2 layers (m)
 zmaxLayer3_upper       => mpar_data%zmaxLayer3_upper_,    & ! maximum layer depth for the 3rd layer when > 3 layers (m)
 zmaxLayer4_upper       => mpar_data%zmaxLayer4_upper_    & ! maximum layer depth for the 4th layer when > 4 layers (m)
 ! diagnostic scalar variables
!  scalarSnowfall         => flux_data%var(iLookFLUX%scalarSnowfall)%dat(1),       & ! snowfall flux (kg m-2 s-1)
!  scalarSnowfallTemp     => diag_data%var(iLookDIAG%scalarSnowfallTemp)%dat(1),   & ! computed temperature of fresh snow (K)
!  scalarSnowDepth        => prog_data%var(iLookPROG%scalarSnowDepth)%dat(1),      & ! total snow depth (m)
!  scalarSWE              => prog_data%var(iLookPROG%scalarSWE)%dat(1)             & ! SWE (kg m-2)
 )  ! end associate statement

 ! ---------------------------------------------------------------------------------------------------

 ! initialize flag to denote that a layer was divided
 divideLayer=.false.

 ! identify algorithmic control parameters to syb-divide and combine snow layers
!  zmax_lower = (/zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower/)
!  zmax_upper = (/zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper/)

 ! initialize the number of snow layers
!  nSnow   = indx_data%var(iLookINDEX%nSnow)%dat(1)
!  nSoil   = indx_data%var(iLookINDEX%nSoil)%dat(1)
!  nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)

    call layerDivide_kernel<<<blocks,threads>>>(nGRU,prog_data%mLayerTemp,prog_data%mLayerVolFracIce, prog_data%mLayerVolFracLiq, prog_data%mLayerVolFracWat, &
  prog_data%mLayerEnthalpy, prog_data%mLayerDepth, prog_data%mLayerHeight, prog_data%iLayerHeight, &
  indx_data%nLayers_d, &
  diag_data%mLayerVolHtCapBulk, diag_data%mLayerCm, diag_data%mLayerThermalC, diag_data%iLayerThermalC, &
    diag_data%mLayerEnthTemp, diag_data%mLayerFracLiqSnow, diag_data%mLayerThetaResid, diag_data%mLayerPoreSpace, &
    diag_data%mLayerMeltFreeze, diag_data%mLayerVolFracAir, diag_data%balanceLayerNrg, diag_data%balanceLayerMass, &
    flux_data%data, &
    flux_data%ixiLayerConductiveFlux_start,flux_data%ixiLayerConductiveFlux_end, flux_data%ixiLayerAdvectiveFlux_start,flux_data%ixiLayerAdvectiveFlux_end, &
  flux_data%ixiLayerNrgFlux_start,flux_data%ixiLayerNrgFlux_end, flux_data%ixmLayerNrgFlux_start,flux_data%ixmLayerNrgFlux_end, &
  flux_data%ixiLayerLiqFluxSnow_start,flux_data%ixiLayerLiqFluxSnow_end,flux_data%ixmLayerLiqFluxSnow_start,flux_data%ixmLayerLiqFluxSnow_end, &
  indx_data%layerType, indx_data%ixHydType, &
  indx_data%ixSnowSoilNrg, indx_data%ixSnowOnlyNrg, &
  indx_data%ixSnowSoilHyd, indx_data%ixSnowOnlyHyd, &
  indx_data%ixNrgLayer, indx_data%ixHydLayer, indx_data%nSnow, &
  divideLayer, ix_snowLayers, zMax, &
  zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower, &
  zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper, &
  maxSnowLayers, veryBig, &
  model_decisions%canopySrad, model_decisions%alb_method,&
  mpar_data%albedoMax_, mpar_data%albedoMaxVisible_, mpar_data%albedoMaxNearIR_, &
  prog_data%spectralSnowAlbedoDiffuse, diag_data%spectralSnowAlbedoDirect, prog_data%scalarSnowAlbedo, &
  prog_data%scalarSnowDepth, prog_data%scalarSWE, &
  fc_param, mpar_data%Frad_vis_,indx_data%nSoil,tooMuchMelt)

 ! end associate variables in data structure
 end associate

 end subroutine layerDivide

 attributes(global) subroutine layerDivide_kernel(nGRU,mLayerTemp, mLayerVolFracIce, mLayerVolFracLiq, mLayerVolFracWat, &
  mLayerEnthalpy, mLayerDepth, mLayerHeight, iLayerHeight, &
  nLayers, &
  mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
    mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
    mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass, &
    flux_data, &
    iLayerConductiveFlux_start,iLayerConductiveFlux_end, iLayerAdvectiveFlux_start,iLayerAdvectiveFlux_end, &
  iLayerNrgFlux_start,iLayerNrgFlux_end, mLayerNrgFlux_start,mLayerNrgFlux_end, &
  iLayerLiqFluxSnow_start,iLayerLiqFluxSnow_end,mLayerLiqFluxSnow_start,mLayerLiqFluxSnow_end, &
  layerType, ixHydType, &
  ixSnowSoilNrg, ixSnowOnlyNrg, &
  ixSnowSoilHyd, ixSnowOnlyHyd, &
  ixNrgLayer, ixHydLayer, nSnow, &
  divideLayer, ix_snowLayers, zMax, &
  zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower, &
  zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper, &
  maxSnowLayers, veryBig, &
  canopySrad, alb_method, &
  albedoMax, albedoMaxVisible, albedoMaxNearIR, &
  spectralSnowAlbedoDiffuse, spectralSnowAlbedoDirect, scalarSnowAlbedo, &
  scalarSnowDepth, scalarSWE, &
  fc_param, Frad_vis, nSoil,tooMuchMelt)
  integer(i4b),value :: nGRU
  real(rkind),intent(inout) :: mLayerTemp(:,:), mLayerVolFracIce(:,:), mLayerVolFracLiq(:,:), mLayerVolFracWat(:,:)
  real(rkind),intent(inout) :: mLayerEnthalpy(:,:), mLayerDepth(:,:), mLayerHeight(0:,:), iLayerHeight(0:,:)
  integer(i4b),intent(inout) :: nLayers(:)
  real(rkind),intent(inout) :: mLayerVolHtCapBulk(:,:), mLayerCm(:,:), mLayerThermalC(:,:), iLayerThermalC(0:,:)
  real(rkind),intent(inout) :: mLayerEnthTemp(:,:), mLayerFracLiqSnow(:,:), mLayerThetaResid(:,:), mLayerPoreSpace(:,:)
  real(rkind),intent(inout) :: mLayerMeltFreeze(:,:), mLayerVolFracAir(:,:), balanceLayerNrg(:,:), balanceLayerMass(:,:)
  real(rkind),intent(inout) :: flux_data(:,:)
  integer(i4b),intent(in),value :: iLayerConductiveFlux_start,iLayerConductiveFlux_end
  integer(i4b),intent(in),value :: iLayerAdvectiveFlux_start,iLayerAdvectiveFlux_end
  integer(i4b),intent(in),value :: iLayerNrgFlux_start,iLayerNrgFlux_end, mLayerNrgFlux_start,mLayerNrgFlux_end
  integer(i4b),intent(in),value :: iLayerLiqFluxSnow_start,iLayerLiqFluxSnow_end,mLayerLiqFluxSnow_start,mLayerLiqFluxSnow_end
  integer(i4b),intent(inout) :: layerType(:,:), ixHydType(:,:)
  integer(i4b),intent(inout) :: ixSnowSoilNrg(:,:), ixSnowOnlyNrg(:,:)
  integer(i4b),intent(inout) :: ixSnowSoilHyd(:,:), ixSnowOnlyHyd(:,:)
  integer(i4b),intent(inout) :: ixNrgLayer(:,:), ixHydLayer(:,:)
  ! integer(i4b),intent(inout) :: ixLayerState(:,:), ixLayerActive(:,:)
  integer(i4b),intent(inout) :: nSnow(:)
  logical(lgt),intent(inout) :: divideLayer(:)
  integer(i4b),intent(in) :: ix_snowLayers
  real(rkind),intent(in) :: zmax(:)
  real(rkind),intent(in) :: zmaxLayer1_lower(:), zmaxLayer2_lower(:), zmaxLayer3_lower(:), zmaxLayer4_lower(:)
  real(rkind),intent(in) :: zmaxLayer1_upper(:), zmaxLayer2_upper(:), zmaxLayer3_upper(:), zmaxLayer4_upper(:)
  real(rkind),intent(in),value :: veryBig
  integer(i4b),intent(in),value :: maxSnowLayers
  integer(i4b),intent(in) :: canopySrad, alb_method
  real(rkind),intent(in) :: albedoMax(:), albedoMaxVisible(:), albedoMaxNearIR(:)
  real(rkind),intent(inout) :: spectralSnowAlbedoDiffuse(:,:), spectralSnowAlbedoDirect(:,:), scalarSnowAlbedo(:)
  real(rkind),intent(inout) :: scalarSnowDepth(:), scalarSWE(:)
  real(rkind),intent(in) :: fc_param(:), Frad_vis(:)
  integer(i4b),intent(in),value :: nSoil
  logical(lgt),intent(in),value :: tooMuchMelt


      integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return
  call layerDivide_device(mLayerTemp(:,iGRU), mLayerVolFracIce(:,iGRU), mLayerVolFracLiq(:,iGRU), mLayerVolFracWat(:,iGRU), &
  mLayerEnthalpy(:,iGRU), mLayerDepth(:,iGRU), mLayerHeight(:,iGRU), iLayerHeight(:,iGRU), &
  nLayers(iGRU), mLayerVolHtCapBulk(:,iGRU), mLayerCm(:,iGRU), mLayerThermalC(:,iGRU), iLayerThermalC(:,iGRU), &
    mLayerEnthTemp(:,iGRU), mLayerFracLiqSnow(:,iGRU), mLayerThetaResid(:,iGRU), mLayerPoreSpace(:,iGRU), &
    mLayerMeltFreeze(:,iGRU), mLayerVolFracAir(:,iGRU), balanceLayerNrg(:,iGRU), balanceLayerMass(:,iGRU), &
    flux_data(iLayerConductiveFlux_start:iLayerConductiveFlux_end,iGRU), flux_data(iLayerAdvectiveFlux_start:iLayerAdvectiveFlux_end,iGRU), &
  flux_data(iLayerNrgFlux_start:iLayerNrgFlux_end,iGRU), flux_data(mLayerNrgFlux_start:mLayerNrgFlux_end,iGRU), &
  flux_data(iLayerLiqFluxSnow_start:iLayerLiqFluxSnow_end,iGRU),flux_data(mLayerLiqFluxSnow_start:mLayerLiqFluxSnow_end,iGRU), &
  layerType(:,iGRU), ixHydType(:,iGRU), &
  ixSnowSoilNrg(:,iGRU), ixSnowOnlyNrg(:,iGRU), &
  ixSnowSoilHyd(:,iGRU), ixSnowOnlyHyd(:,iGRU), &
  ixNrgLayer(:,iGRU), ixHydLayer(:,iGRU), nSnow(iGRU), &
  divideLayer(iGRU), ix_snowLayers, zMax(iGRU), &
  zmaxLayer1_lower(iGRU), zmaxLayer2_lower(iGRU), zmaxLayer3_lower(iGRU), zmaxLayer4_lower(iGRU), &
  zmaxLayer1_upper(iGRU), zmaxLayer2_upper(iGRU), zmaxLayer3_upper(iGRU), zmaxLayer4_upper(iGRU), &
  maxSnowLayers, veryBig, &
  canopySrad, alb_method, &
  albedoMax(iGRU), albedoMaxVisible(iGRU), albedoMaxNearIR(iGRU), &
  spectralSnowAlbedoDiffuse(:,iGRU), spectralSnowAlbedoDirect(:,iGRU), scalarSnowAlbedo(iGRU), &
  scalarSnowDepth(iGRU), scalarSWE(iGRU), &
  fc_param(iGRU), Frad_vis(iGRU),nSoil,tooMuchMelt)
 end subroutine

 attributes(device) subroutine layerDivide_device(mLayerTemp, mLayerVolFracIce, mLayerVolFracLiq, mLayerVolFracWat, &
  mLayerEnthalpy, mLayerDepth, mLayerHeight, iLayerHeight, &
  nLayers, &
  mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
    mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
    mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass, &
    iLayerConductiveFlux, iLayerAdvectiveFlux, &
  iLayerNrgFlux, mLayerNrgFlux, &
  iLayerLiqFluxSnow,mLayerLiqFluxSnow, &
  layerType, ixHydType, &
  ixSnowSoilNrg, ixSnowOnlyNrg, &
  ixSnowSoilHyd, ixSnowOnlyHyd, &
  ixNrgLayer, ixHydLayer, nSnow, &
  divideLayer, ix_snowLayers, zMax, &
  zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower, &
  zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper, &
  maxSnowLayers, veryBig, &
  canopySrad, alb_method, &
  albedoMax, albedoMaxVisible, albedoMaxNearIR, &
  spectralSnowAlbedoDiffuse, spectralSnowAlbedoDirect, scalarSnowAlbedo, &
  scalarSnowDepth, scalarSWE, &
  fc_param, Frad_vis, nSoil,tooMuchMelt)
   USE snow_utils_module,only:fracliquid,templiquid          ! functions to compute temperature/liquid water

   real(rkind),intent(inout) :: mLayerTemp(:), mLayerVolFracIce(:), mLayerVolFracLiq(:), mLayerVolFracWat(:)
  real(rkind),intent(inout) :: mLayerEnthalpy(:), mLayerDepth(:), mLayerHeight(0:), iLayerHeight(0:)
  integer(i4b),intent(inout) :: nLayers
  real(rkind),intent(inout) :: mLayerVolHtCapBulk(:), mLayerCm(:), mLayerThermalC(:), iLayerThermalC(0:)
  real(rkind),intent(inout) :: mLayerEnthTemp(:), mLayerFracLiqSnow(:), mLayerThetaResid(:), mLayerPoreSpace(:)
  real(rkind),intent(inout) :: mLayerMeltFreeze(:), mLayerVolFracAir(:), balanceLayerNrg(:), balanceLayerMass(:)
  real(rkind),intent(inout) :: iLayerConductiveFlux(0:)
  real(rkind),intent(inout) :: iLayerAdvectiveFlux(0:)
  real(rkind),intent(inout) :: iLayerNrgFlux(0:)
  real(rkind),intent(inout) :: mLayerNrgFlux(:)
  real(rkind),intent(inout) :: iLayerLiqFluxSnow(0:)
  real(rkind),intent(inout) :: mLayerLiqFluxSnow(:)
  integer(i4b),intent(inout) :: layerType(:), ixHydType(:)
  integer(i4b),intent(inout) :: ixSnowSoilNrg(:), ixSnowOnlyNrg(:)
  integer(i4b),intent(inout) :: ixSnowSoilHyd(:), ixSnowOnlyHyd(:)
  integer(i4b),intent(inout) :: ixNrgLayer(:), ixHydLayer(:)
  integer(i4b),intent(inout) :: nSnow
  logical(lgt),intent(inout) :: divideLayer
  integer(i4b),intent(in) :: ix_snowLayers
  real(rkind),intent(in) :: zmax
  real(rkind),intent(in) :: zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower
  real(rkind),intent(in) :: zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper
  real(rkind),intent(in) :: veryBig
  integer(i4b),intent(in) :: maxSnowLayers
  integer(i4b),intent(in) :: canopySrad, alb_method
  real(rkind),intent(in) :: albedoMax, albedoMaxVisible, albedoMaxNearIR
  real(rkind),intent(inout) :: spectralSnowAlbedoDiffuse(:), spectralSnowAlbedoDirect(:), scalarSnowAlbedo
  real(rkind),intent(inout) :: scalarSnowDepth, scalarSWE
  real(rkind),intent(in) :: fc_param, Frad_vis
  integer(i4b),intent(in) :: nSoil
  logical(lgt),intent(in) :: tooMuchMelt

   real(rkind)                     :: depthOriginal          ! original layer depth before sub-division (m)
 real(rkind),parameter           :: fracTop=0.5_rkind      ! fraction of old layer used for the top layer
  integer(i4b),parameter          :: ixVisible=1            ! named variable to define index in array of visible part of the spectrum
 integer(i4b),parameter          :: ixNearIR=2             ! named variable to define index in array of near IR part of the spectrum
 real(rkind)                     :: fracLiq                ! fraction of liquid water (-)
 real(rkind)                     :: volFracWater           ! volumetric fraction of total water, liquid and ice (-)
 real(rkind)                     :: surfaceLayerSoilTemp   ! temperature of the top soil layer (K)
 real(rkind)                     :: maxFrozenSnowTemp      ! maximum temperature when effectively all water is frozen (K)
 real(rkind),parameter           :: unfrozenLiq=0.01_rkind ! unfrozen liquid water used to compute maxFrozenSnowTemp (-)
 logical(lgt)                    :: createLayer            ! flag to indicate we are creating a new snow layer

 integer(i4b) :: iLayer, jLayer
 integer(i4b) :: nCheck
 real(rkind) :: zmaxCheck

 if (tooMuchMelt) return
 ! ***** special case of no snow layers
 if(nSnow==0)then
  divideLayer=.false.

  ! check if create the first snow layer
  select case(ix_snowLayers)
   case(sameRulesAllLayers);    createLayer = (scalarSnowDepth > zmax)
   case(rulesDependLayerIndex); createLayer = (scalarSnowDepth > zmaxLayer1_lower)
  !  case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
  end select ! (option to combine/sub-divide snow layers)

  ! ** create a new snow layer
  if(createLayer)then

   ! flag that the layers have changed
   divideLayer=.true.

   ! add a layer to all model variables
   iLayer=0 ! (layer to divide: 0 is the special case of "snow without a layer")
   call addModelLayer_prog(mLayerTemp, mLayerVolFracIce, mLayerVolFracLiq, mLayerVolFracWat, &
  mLayerEnthalpy, mLayerDepth, mLayerHeight, iLayerHeight, &
  nLayers, iLayer)
    call addModelLayer_diag(mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
    mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
    mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass)
  call addModelLayer_flux(iLayerConductiveFlux, iLayerAdvectiveFlux, &
  iLayerNrgFlux, mLayerNrgFlux, &
  iLayerLiqFluxSnow,mLayerLiqFluxSnow)
  call addModelLayer_indx(layerType, ixHydType, &
  ixSnowSoilNrg, ixSnowOnlyNrg, &
  ixSnowSoilHyd, ixSnowOnlyHyd, &
  ixNrgLayer, ixHydLayer)


   ! get the layer depth
   mLayerDepth(1) = scalarSnowDepth

   ! compute surface layer temperature
   surfaceLayerSoilTemp = mLayerTemp(2)    ! temperature of the top soil layer (K)
   maxFrozenSnowTemp    = templiquid(unfrozenLiq,fc_param)               ! snow temperature at fraction "unfrozenLiq" (K)
   mLayerTemp(1)        = min(maxFrozenSnowTemp,surfaceLayerSoilTemp)    ! snow temperature  (K)

   ! compute the fraction of liquid water associated with the layer temperature
   fracLiq = fracliquid(mLayerTemp(1),fc_param)

   ! compute volumeteric fraction of liquid water and ice
   volFracWater = (scalarSWE/scalarSnowDepth)/iden_water  ! volumetric fraction of total water (liquid and ice)
   mLayerVolFracIce(1) = (1._rkind - fracLiq)*volFracWater*(iden_water/iden_ice)   ! volumetric fraction of ice (-)
   mLayerVolFracLiq(1) =             fracLiq *volFracWater                         ! volumetric fraction of liquid water (-)

   ! initialize albedo
   ! NOTE: albedo is computed within the Noah-MP radiation routine
   if(canopySrad /= noah_mp)then
    select case(alb_method)
     ! (constant decay rate -- albedo the same for all spectral bands)
     case(constantDecay)
      scalarSnowAlbedo          = albedoMax
      spectralSnowAlbedoDiffuse = albedoMax
     ! (variable decay rate)
     case(variableDecay)
      spectralSnowAlbedoDiffuse(ixVisible) = albedoMaxVisible
      spectralSnowAlbedoDiffuse(ixNearIR)  = albedoMaxNearIR
      scalarSnowAlbedo                  = (        Frad_vis)*albedoMaxVisible + &
                                          (1._rkind - Frad_vis)*albedoMaxNearIR
    !  case default; err=20; message=trim(message)//'unable to identify option for snow albedo'; return
    end select  ! identify option for snow albedo
    ! set direct albedo to diffuse albedo
    spectralSnowAlbedoDirect = spectralSnowAlbedoDiffuse
   end if  ! (if NOT using the Noah-MP radiation routine)

  end if  ! if creating a new layer

 ! end special case of nSnow=0
 ! ********************************************************************************************************************
 ! ********************************************************************************************************************

 ! ***** sub-divide snow layers, if necessary
 else ! if nSnow>0

  ! identify the number of layers to check for need for sub-division
  nCheck = min(nSnow, maxSnowLayers-1) ! the depth of the last layer, if it exists, does not have a maximum value
  ! loop through all layers, and sub-divide a given layer, if necessary
  do iLayer=1,nCheck
   divideLayer=.false.

   ! identify the maximum depth of the layer
   select case(ix_snowLayers)
    case(sameRulesAllLayers)
     if (nCheck >= maxSnowLayers-1) then
      ! make sure we don't divide so make very big
      zmaxCheck = veryBig
     else
      zmaxCheck = zmax
     end if
    case(rulesDependLayerIndex)
     if(iLayer == nSnow)then
      select case(iLayer)
      case(1); zmaxCheck = zmaxLayer1_lower
      case(2); zmaxCheck = zmaxLayer2_lower
      case(3); zmaxCheck = zmaxLayer3_lower
      case(4); zmaxCheck = zmaxLayer4_lower
      end select
     else
      select case(iLayer)
      case(1); zmaxCheck = zmaxLayer1_upper
      case(2); zmaxCheck = zmaxLayer2_upper
      case(3); zmaxCheck = zmaxLayer3_upper
      case(4); zmaxCheck = zmaxLayer4_upper
      end select
     end if
    ! case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
   end select ! (option to combine/sub-divide snow layers)

   ! check the need to sub-divide
   if(mLayerDepth(iLayer) > zmaxCheck)then

    ! flag that layers were divided
    divideLayer=.true.

        ! add a layer to all model variables
   call addModelLayer_prog(mLayerTemp, mLayerVolFracIce, mLayerVolFracLiq, mLayerVolFracWat, &
  mLayerEnthalpy, mLayerDepth, mLayerHeight, iLayerHeight, &
  nLayers, iLayer)
  call addModelLayer_diag(mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
    mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
    mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass)
  call addModelLayer_flux(iLayerConductiveFlux, iLayerAdvectiveFlux, &
  iLayerNrgFlux, mLayerNrgFlux, &
  iLayerLiqFluxSnow,mLayerLiqFluxSnow)
  call addModelLayer_indx(layerType, ixHydType, &
  ixSnowSoilNrg, ixSnowOnlyNrg, &
  ixSnowSoilHyd, ixSnowOnlyHyd, &
  ixNrgLayer, ixHydLayer)

      ! define the layer depth
    depthOriginal = mLayerDepth(iLayer)
    mLayerDepth(iLayer)   = fracTop*depthOriginal
    mLayerDepth(iLayer+1) = (1._rkind - fracTop)*depthOriginal
    exit  ! NOTE: only sub-divide one layer per substep

   end if   ! (if sub-dividing layer)

  end do  ! (looping through layers)

end if

 ! update coordinates
 if(divideLayer)then
  ! update the layer type
  layerType(1:nSnow+1)         = iname_snow
  layerType(nSnow+2:nLayers+1) = iname_soil

  ! identify the number of snow and soil layers, and check all is a-OK
  nSnow   = count(layerType(1:nLayers+1)==iname_snow)
  ! nSoil   = count(layerType(1:nLayers+1)==iname_soil)
  nLayers = nSnow + nSoil

  ! re-set coordinate variables
  iLayerHeight(0) = -scalarSnowDepth
  do jLayer=1,nLayers
   iLayerHeight(jLayer) = iLayerHeight(jLayer-1) + mLayerDepth(jLayer)
   mLayerHeight(jLayer) = (iLayerHeight(jLayer-1) + iLayerHeight(jLayer))/2._rkind
  end do

    ! check
  ! if(abs(sum(mLayerDepth(1:nSnow)) - scalarSnowDepth) > snowDepthTol)then
  !  print*, 'nSnow = ', nSnow
  !  write(*,'(a,1x,f30.25,1x)') 'sum(mLayerDepth(1:nSnow)) = ', sum(mLayerDepth(1:nSnow))
  !  write(*,'(a,1x,f30.25,1x)') 'scalarSnowDepth           = ', scalarSnowDepth
  !  write(*,'(a,1x,f30.25,1x)') 'snowDepthTol              = ', snowDepthTol
  !  message=trim(message)//'sum of layer depths does not equal snow depth'
  !  err=20; return
  ! end if

 end if
end subroutine


 attributes(device) subroutine addModelLayer_prog(mLayerTemp, mLayerVolFracIce, mLayerVolFracLiq, mLayerVolFracWat, &
  mLayerEnthalpy, mLayerDepth, mLayerHeight, iLayerHeight, &
  nLayers, ix_divide)
  real(rkind),intent(inout) :: mLayerTemp(:), mLayerVolFracIce(:), mLayerVolFracLiq(:), mLayerVolFracWat(:)
  real(rkind),intent(inout) :: mLayerEnthalpy(:), mLayerDepth(:), mLayerHeight(0:), iLayerHeight(0:)
  integer(i4b),intent(inout) :: nLayers, ix_divide
  integer(i4b) :: iLayer
    mLayerEnthalpy = realMissing
    mLayerHeight = realMissing
    iLayerHeight = realMissing
    mLayerVolFracWat = realMissing
    do iLayer=nLayers+1,ix_divide+1,-1
      mLayerDepth(iLayer) = mLayerDepth(iLayer-1)
      mLayerTemp(iLayer) = mLayerTemp(iLayer-1)
      mLayerVolFracIce(iLayer) = mLayerVolFracIce(iLayer-1)
      mLayerVolFracLiq(iLayer) = mLayerVolFracLiq(iLayer-1)
    end do
 end subroutine

 attributes(device) subroutine addModelLayer_diag(mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
    mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
    mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass)
  real(rkind),intent(inout) :: mLayerVolHtCapBulk(:), mLayerCm(:), mLayerThermalC(:), iLayerThermalC(0:)
  real(rkind),intent(inout) :: mLayerEnthTemp(:), mLayerFracLiqSnow(:), mLayerThetaResid(:), mLayerPoreSpace(:)
  real(rkind),intent(inout) :: mLayerMeltFreeze(:), mLayerVolFracAir(:), balanceLayerNrg(:), balanceLayerMass(:)

  mLayerVolHtCapBulk = realMissing
  mLayerCm = realMissing
  mLayerThermalC = realMissing
  iLayerThermalC = realMissing
  mLayerEnthTemp = realMissing
  mLayerFracLiqSnow = realMissing
  mLayerThetaResid = realMissing
  mLayerPoreSpace = realMissing
  mLayerMeltFreeze = realMissing
  mLayerVolFracAir = realMissing
  balanceLayerNrg = realMissing
  balanceLayerMass = realMissing
  
end subroutine

attributes(device) subroutine addModelLayer_flux(iLayerConductiveFlux, iLayerAdvectiveFlux, &
  iLayerNrgFlux, mLayerNrgFlux, &
  iLayerLiqFluxSnow,mLayerLiqFluxSnow)
  real(rkind),intent(inout) :: iLayerConductiveFlux(0:)
  real(rkind),intent(inout) :: iLayerAdvectiveFlux(0:)
  real(rkind),intent(inout) :: iLayerNrgFlux(0:)
  real(rkind),intent(inout) :: mLayerNrgFlux(:)
  real(rkind),intent(inout) :: iLayerLiqFluxSnow(0:)
  real(rkind),intent(inout) :: mLayerLiqFluxSnow(:)

  iLayerConductiveFlux = realMissing
  iLayerAdvectiveFlux = realMissing
  iLayerNrgFlux = realMissing
  mLayerNrgFlux = realMissing
  iLayerLiqFluxSnow = realMissing
  mLayerLiqFluxSnow = realMissing


end subroutine

attributes(device) subroutine addModelLayer_indx(layerType, ixHydType, &
  ixSnowSoilNrg, ixSnowOnlyNrg, &
  ixSnowSoilHyd, ixSnowOnlyHyd, &
  ixNrgLayer, ixHydLayer)
  integer(i4b),intent(inout) :: layerType(:), ixHydType(:)
  integer(i4b),intent(inout) :: ixSnowSoilNrg(:), ixSnowOnlyNrg(:)
  integer(i4b),intent(inout) :: ixSnowSoilHyd(:), ixSnowOnlyHyd(:)
  integer(i4b),intent(inout) :: ixNrgLayer(:), ixHydLayer(:)
  ! integer(i4b),intent(inout) :: ixLayerState(:), ixLayerActive(:)

  layerType = integerMissing
  ixHydType = integerMissing
  ixSnowSoilNrg = integerMissing
  ixSnowOnlyNrg = integerMissing
  ixSnowSoilHyd = integerMissing
  ixSnowOnlyHyd = integerMissing
  ixNrgLayer = integerMissing
  ixHydLayer = integerMissing
  ! ixLayerState = integerMissing
  ! ixLayerActive = integerMissing
end subroutine

 ! ************************************************************************************************
 ! private subroutine addModelLayer: add an additional layer to all model vectors
 ! ************************************************************************************************
 subroutine addModelLayer(dataStruct,metaStruct,ix_divide,nSnow,nLayers,err,message)
 USE var_lookup,only:iLookVarType                     ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName           ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc                ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength          ! data vectors with variable length dimension
 USE data_types,only:var_info                         ! metadata structure
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 class(*),intent(inout)          :: dataStruct        ! data structure
 type(var_info),intent(in)       :: metaStruct(:)     ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in)         :: ix_divide         ! index of the layer to divide
 integer(i4b),intent(in)         :: nSnow,nLayers     ! number of snow layers, total number of layers
 ! output: error control
 integer(i4b),intent(out)        :: err               ! error code
 character(*),intent(out)        :: message           ! error message
 ! ---------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: ivar              ! index of model variable
 integer(i4b)                    :: ix_lower          ! lower bound of the vector
 integer(i4b)                    :: ix_upper          ! upper bound of the vector
 logical(lgt)                    :: stateVariable     ! .true. if variable is a state variable
 real(rkind),allocatable         :: tempVec_rkind(:)  ! temporary vector (double precision)
 integer(i4b),allocatable        :: tempVec_i4b(:)    ! temporary vector (integer)
 character(LEN=256)              :: cmessage          ! error message of downwind routine
 ! ---------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='addModelLayer/'

 ! ***** add a layer to each model variable
 do ivar=1,size(metaStruct)

  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=nSnow
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=nSnow
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nLayers
   case default; cycle
  end select

  ! identify whether it is a state variable
  select case(trim(metaStruct(ivar)%varname))
   case('mLayerDepth','mLayerTemp','mLayerVolFracIce','mLayerVolFracLiq'); stateVariable=.true.
   case default; stateVariable=.false.
  end select

  ! divide layers
  select type(dataStruct)

   ! ** double precision
   type is (var_dlength)
    ! check allocated
    if(.not.allocated(dataStruct%var(ivar)%dat))then; err=20; message='data vector is not allocated'; return; end if
    ! assign the data vector to the temporary vector
    call cloneStruc(tempVec_rkind, ix_lower, source=dataStruct%var(ivar)%dat, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! reallocate space for the new vector
    deallocate(dataStruct%var(ivar)%dat,stat=err)
    if(err/=0)then; err=20; message='problem in attempt to deallocate memory for data vector'; return; end if
    allocate(dataStruct%var(ivar)%dat(ix_lower:ix_upper+1),stat=err)
    if(err/=0)then; err=20; message='problem in attempt to reallocate memory for data vector'; return; end if
    ! populate the state vector
    if(stateVariable)then
     if(ix_upper > 0)then  ! (only copy data if the vector exists -- can be a variable for snow, with no layers)
      if(ix_divide > 0)then
       dataStruct%var(ivar)%dat(1:ix_divide) = tempVec_rkind(1:ix_divide)  ! copy data
       dataStruct%var(ivar)%dat(ix_divide+1) = tempVec_rkind(ix_divide)    ! repeat data for the sub-divided layer
      end if
      if(ix_upper > ix_divide) &
       dataStruct%var(ivar)%dat(ix_divide+2:ix_upper+1) = tempVec_rkind(ix_divide+1:ix_upper)  ! copy data
     end if  ! if the vector exists
    ! not a state variable
    else
     dataStruct%var(ivar)%dat(:) = realMissing
    end if
    ! deallocate the temporary vector: strictly not necessary, but include to be safe
    deallocate(tempVec_rkind,stat=err)
    if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

   ! ** integer
   type is (var_ilength)
    ! check allocated
    if(.not.allocated(dataStruct%var(ivar)%dat))then; err=20; message='data vector is not allocated'; return; end if
    ! assign the data vector to the temporary vector
    call cloneStruc(tempVec_i4b, ix_lower, source=dataStruct%var(ivar)%dat, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! reallocate space for the new vector
    deallocate(dataStruct%var(ivar)%dat,stat=err)
    if(err/=0)then; err=20; message='problem in attempt to deallocate memory for data vector'; return; end if
    allocate(dataStruct%var(ivar)%dat(ix_lower:ix_upper+1),stat=err)
    if(err/=0)then; err=20; message='problem in attempt to reallocate memory for data vector'; return; end if
    ! populate the state vector
    if(stateVariable)then
     if(ix_upper > 0)then  ! (only copy data if the vector exists -- can be a variable for snow, with no layers)
      if(ix_divide > 0)then
       dataStruct%var(ivar)%dat(1:ix_divide) = tempVec_i4b(1:ix_divide)  ! copy data
       dataStruct%var(ivar)%dat(ix_divide+1) = tempVec_i4b(ix_divide)    ! repeat data for the sub-divided layer
      end if
      if(ix_upper > ix_divide) &
       dataStruct%var(ivar)%dat(ix_divide+2:ix_upper+1) = tempVec_i4b(ix_divide+1:ix_upper)  ! copy data
     end if  ! if the vector exists
    ! not a state variable
    else
     dataStruct%var(ivar)%dat(:) = integerMissing
    end if
    ! deallocate the temporary vector: strictly not necessary, but include to be safe
    deallocate(tempVec_i4b,stat=err)
    if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

   ! check that we found the data type
   class default; err=20; message=trim(message)//'unable to identify the data type'; return

  end select ! dependence on data types

 end do  ! looping through variables

 end subroutine addModelLayer

end module layerDivide_module
