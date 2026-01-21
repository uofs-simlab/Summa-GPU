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

module volicePack_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_d,            & ! data vector (rkind)
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength,      & ! data vector with variable length dimension (rkind)
                    model_options       ! defines the model decisions

! named variables for snow and soil
USE globalData,only:iname_snow          ! named variables for snow
USE globalData,only:iname_soil          ! named variables for soil

! named variables for parent structures
USE var_lookup,only:iLookINDEX          ! named variables for structure elements

! physical constants
USE multiconst,only:&
                    iden_ice, & ! intrinsic density of ice    (kg m-3)
                    iden_water  ! intrinsic density of water  (kg m-3)

! privacy
implicit none
private
public::volicePack_device
public::newsnwfall

contains


 ! ************************************************************************************************
 ! public subroutine volicePack: combine and sub-divide layers if necessary)
 ! ************************************************************************************************
 attributes(global) subroutine volicePack_kernel(nGRU,ix_snowLayers,nSnow,nLayers,nSoil, &
 mLayerTemp,mLayerVolFracIce,mLayerVolFracLiq, &
mLayerVolFracWat,mLayerEnthalpy,mLayerDepth, mLayerHeight, iLayerHeight, &
mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass,&
flux_data, &
iLayerConductiveFlux_start,iLayerConductiveFlux_end, &
iLayerAdvectiveFlux_start,iLayerAdvectiveFlux_end, &
iLayerNrgFlux_start,iLayerNrgFlux_end, mLayerNrgFlux_start,mLayerNrgFlux_end, &
iLayerLiqFluxSnow_start,iLayerLiqFluxSnow_end,mLayerLiqFluxSnow_start,mLayerLiqFluxSnow_end, &
    layerType,ixHydType, &
    ixSnowSoilNrg, ixSnowOnlyNrg, &
    ixSnowSoilHyd, ixSnowOnlyHyd, &
    ixNrgLayer, ixHydLayer,&
    tooMuchMelt, zMin,&
    zMinLayer1,zMinLayer2,zminLayer3,zminLayer4,zminLayer5, &
    modifiedLayers,snowfrz_scale,h_lookup_d,t_lookup_d, &
    scalarSnowDepth, scalarSWE, &
    zmax, &
    zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower, &
    zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper, &
    canopySrad, alb_method, &
    albedoMax, albedoMaxVisible, albedoMaxNearIR, &
    spectralSnowAlbedoDiffuse, spectralSnowAlbedoDirect, scalarSnowAlbedo, &
    veryBig,Frad_vis,maxSnowLayers)
   integer(i4b),intent(in),value :: nGRU
        integer(i4b),intent(in) :: ix_snowLayers
    integer(i4b),intent(inout) :: nSnow(:),nLayers(:)
    integer(i4b),intent(in),value :: nSoil
    real(rkind),intent(inout) :: mLayerTemp(:,:)
    real(rkind),intent(inout) :: mLayerVolFracIce(:,:)
    real(rkind),intent(inout) :: mLayerVolFracLiq(:,:)
    real(rkind),intent(inout) :: mLayerVolFracWat(:,:)
    real(rkind),intent(inout) :: mLayerEnthalpy(:,:)
    real(rkind),intent(inout) :: mLayerDepth(:,:)
    real(rkind),intent(inout) :: mLayerHeight(0:,:)
    real(rkind),intent(inout) :: iLayerHeight(0:,:)
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
    logical(lgt),intent(in),value :: tooMuchMelt
    real(rkind),intent(in) :: zMin
    real(rkind),intent(in) :: zMinLayer1,zMinLayer2,zminLayer3,zminLayer4,zminLayer5
    logical(lgt),intent(inout) :: modifiedLayers(:)
    real(rkind),intent(in) :: snowfrz_scale
    real(rkind),intent(inout) :: h_lookup_d(:), t_lookup_d(:)
    real(rkind),intent(inout) :: scalarSnowDepth(:), scalarSWE(:)
    real(rkind),intent(in) :: zmax
    real(rkind),intent(in) :: zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower
    real(rkind),intent(in) :: zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper
    integer(i4b),intent(in) :: canopySrad, alb_method
    real(rkind),intent(in) :: albedoMax(:), albedoMaxVisible(:), albedoMaxNearIR(:)
    real(rkind),intent(inout) :: spectralSnowAlbedoDiffuse(:,:), spectralSnowAlbedoDirect(:,:)
    real(rkind),intent(inout) :: scalarSnowAlbedo(:)
    real(rkind),intent(in),value :: veryBig
    real(rkind),intent(in) :: Frad_vis
    integer(i4b),intent(in),value :: maxSnowLayers

       integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return

  call volicePack_device(mLayerTemp(:,iGRU), mLayerVolFracIce(:,iGRU), mLayerVolFracLiq(:,iGRU), mLayerVolFracWat(:,iGRU), &
  mLayerEnthalpy(:,iGRU), mLayerDepth(:,iGRU), mLayerHeight(:,iGRU), iLayerHeight(:,iGRU), &
  nLayers(iGRU), &
  mLayerVolHtCapBulk(:,iGRU), mLayerCm(:,iGRU), mLayerThermalC(:,iGRU), iLayerThermalC(:,iGRU), &
    mLayerEnthTemp(:,iGRU), mLayerFracLiqSnow(:,iGRU), mLayerThetaResid(:,iGRU), mLayerPoreSpace(:,iGRU), &
    mLayerMeltFreeze(:,iGRU), mLayerVolFracAir(:,iGRU), balanceLayerNrg(:,iGRU), balanceLayerMass(:,iGRU), &
flux_data(iLayerConductiveFlux_start:iLayerConductiveFlux_end,iGRU), flux_data(iLayerAdvectiveFlux_start:iLayerAdvectiveFlux_end,iGRU), &
  flux_data(iLayerNrgFlux_start:iLayerNrgFlux_end,iGRU), flux_data(mLayerNrgFlux_start:mLayerNrgFlux_end,iGRU), &
  flux_data(iLayerLiqFluxSnow_start:iLayerLiqFluxSnow_end,iGRU),flux_data(mLayerLiqFluxSnow_start:mLayerLiqFluxSnow_end,iGRU), &
  layerType(:,iGRU), ixHydType(:,iGRU), &
  ixSnowSoilNrg(:,iGRU), ixSnowOnlyNrg(:,iGRU), &
  ixSnowSoilHyd(:,iGRU), ixSnowOnlyHyd(:,iGRU), &
  ixNrgLayer(:,iGRU), ixHydLayer(:,iGRU), nSnow(iGRU), &
  ix_snowLayers, zMax, &
  zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower, &
  zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper, &
  maxSnowLayers, veryBig, &
  canopySrad, alb_method, &
  albedoMax(iGRU), albedoMaxVisible(iGRU), albedoMaxNearIR(iGRU), &
  spectralSnowAlbedoDiffuse(:,iGRU), spectralSnowAlbedoDirect(:,iGRU), scalarSnowAlbedo(iGRU), &
  scalarSnowDepth(iGRU), scalarSWE(iGRU), &
  Frad_vis, nSoil,tooMuchMelt, &
  zMin, zMinLayer1, zMinLayer2, zminLayer3, zminLayer4, zminLayer5, &
  snowfrz_scale, h_lookup_d, t_lookup_d, modifiedLayers(iGRU))

end subroutine

 attributes(device) subroutine volicePack_device(mLayerTemp, mLayerVolFracIce, mLayerVolFracLiq, mLayerVolFracWat, &
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
  ix_snowLayers, zMax, &
  zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower, &
  zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper, &
  maxSnowLayers, veryBig, &
  canopySrad, alb_method, &
  albedoMax, albedoMaxVisible, albedoMaxNearIR, &
  spectralSnowAlbedoDiffuse, spectralSnowAlbedoDirect, scalarSnowAlbedo, &
  scalarSnowDepth, scalarSWE, &
  Frad_vis, nSoil,tooMuchMelt, &
  zMin, zMinLayer1, zMinLayer2, zminLayer3, zminLayer4, zminLayer5, &
  snowfrz_scale, h_lookup_d, t_lookup_d, modifiedLayers)
   USE layerMerge_module,only:layerMerge_device   ! merge snow layers if they are too thin
 USE layerDivide_module,only:layerDivide_device ! sub-divide layers if they are too thick
 integer(i4b) :: iLayer

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
  real(rkind),intent(in) :: Frad_vis
  integer(i4b),intent(in) :: nSoil
  logical(lgt),intent(in) :: tooMuchMelt
  real(rkind),intent(in) :: zMin, zMinLayer1, zMinLayer2, zminLayer3, zminLayer4, zminLayer5
  real(rkind),intent(in) :: snowfrz_scale
  real(rkind),intent(in) :: h_lookup_d(:), t_lookup_d(:)
  logical(lgt),intent(inout) :: modifiedLayers

  logical(lgt) :: divideLayer, mergedLayers
  
   call layerDivide_device(mLayerTemp, mLayerVolFracIce, mLayerVolFracLiq, mLayerVolFracWat, &
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
  snowfrz_scale, Frad_vis, nSoil,tooMuchMelt)
   
   call layerMerge_device(ix_snowLayers,nSnow,nLayers,nSoil, &
 mLayerTemp,mLayerVolFracIce,mLayerVolFracLiq, &
mLayerVolFracWat,mLayerEnthalpy,mLayerDepth, mLayerHeight, iLayerHeight, &
mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass,&
iLayerConductiveFlux,iLayerAdvectiveFlux,&
iLayerNrgFlux,mLayerNrgFlux,&
iLayerLiqFluxSnow,mLayerLiqFluxSnow,&
    layerType,ixHydType, &
    ixSnowSoilNrg, ixSnowOnlyNrg, &
    ixSnowSoilHyd, ixSnowOnlyHyd, &
    ixNrgLayer, ixHydLayer,&
    tooMuchMelt, zMin,&
    zMinLayer1,zMinLayer2,zminLayer3,zminLayer4,zminLayer5, &
    mergedLayers,snowfrz_scale,h_lookup_d,t_lookup_d, &
    scalarSnowDepth, scalarSWE)

    ! update the number of layers
   nSnow   = count(layerType==iname_snow)
   nLayers = nSnow + nSoil

  modifiedLayers = (mergedLayers .or. divideLayer)

end subroutine



 ! ************************************************************************************************
 ! public subroutine newsnwfall: add new snowfall to the system
 ! ************************************************************************************************
 attributes(device) subroutine newsnwfall(&
                       ! input: model control
                       dt,                        & ! time step (seconds)
                       snowLayers,                & ! logical flag if snow layers exist
                       fc_param,                  & ! freeezing curve parameter for snow (K-1)
                       ! input: diagnostic scalar variables
                       scalarSnowfallTemp,        & ! computed temperature of fresh snow (K)
                       scalarNewSnowDensity,      & ! computed density of new snow (kg m-3)
                       scalarThroughfallSnow,     & ! throughfall of snow through the canopy (kg m-2 s-1)
                       scalarCanopySnowUnloading, & ! unloading of snow from the canopy (kg m-2 s-1)
                       ! input/output: state variables
                       scalarSWE,                 & ! SWE (kg m-2)
                       scalarSnowDepth,           & ! total snow depth (m)
                       surfaceLayerTemp,          & ! temperature of surface layer (K)
                       surfaceLayerDepth,         & ! depth of surface layer (m)
                       surfaceLayerVolFracIce,    & ! volumetric fraction of ice in surface layer (-)
                       surfaceLayerVolFracLiq,    & ! volumetric fraction of liquid water in surface layer (-)
                       ! output: error control
                       err                ) ! error control
 ! computational modules
 USE snow_utils_module,only:fracliquid,templiquid                  ! functions to compute temperature/liquid water
 ! add new snowfall to the system
 implicit none
 ! input: model control
 real(rkind),intent(in)                 :: dt                         ! time step (seconds)
 logical(lgt),intent(in)                :: snowLayers                 ! logical flag if snow layers exist
 real(rkind),intent(in)                 :: fc_param                   ! freeezing curve parameter for snow (K-1)
 ! input: diagnostic scalar variables
 real(rkind),intent(in)                 :: scalarSnowfallTemp         ! computed temperature of fresh snow (K)
 real(rkind),intent(in)                 :: scalarNewSnowDensity       ! computed density of new snow (kg m-3)
 real(rkind),intent(in)                 :: scalarThroughfallSnow      ! throughfall of snow through the canopy (kg m-2 s-1)
 real(rkind),intent(in)                 :: scalarCanopySnowUnloading  ! unloading of snow from the canopy (kg m-2 s-1)
 ! input/output: state variables
 real(rkind),intent(inout)              :: scalarSWE                  ! SWE (kg m-2)
 real(rkind),intent(inout)              :: scalarSnowDepth            ! total snow depth (m)
 real(rkind),intent(inout)              :: surfaceLayerTemp           ! temperature of surface layer (K)
 real(rkind),intent(inout)              :: surfaceLayerDepth          ! depth of each layer (m)
 real(rkind),intent(inout)              :: surfaceLayerVolFracIce     ! volumetric fraction of ice in surface layer (-)
 real(rkind),intent(inout)              :: surfaceLayerVolFracLiq     ! volumetric fraction of liquid water in surface layer (-)
 ! output: error control
 integer(i4b),intent(out)               :: err                        ! error code
!  character(*),intent(out)               :: message                    ! error message
 ! define local variables
 real(rkind)                            :: newSnowfall                ! new snowfall -- throughfall and unloading (kg m-2 s-1)
 real(rkind)                            :: newSnowDepth               ! new snow depth (m)
 real(rkind),parameter                  :: densityCanopySnow=200._rkind  ! density of snow on the vegetation canopy (kg m-3)
 real(rkind)                            :: totalMassIceSurfLayer      ! total mass of ice in the surface layer (kg m-2)
 real(rkind)                            :: totalDepthSurfLayer        ! total depth of the surface layer (m)
 real(rkind)                            :: volFracWater               ! volumetric fraction of total water, liquid and ice (-)
 real(rkind)                            :: fracLiq                    ! fraction of liquid water (-)
 real(rkind)                            :: SWE                        ! snow water equivalent after snowfall (kg m-2)
 real(rkind)                            :: tempSWE0                   ! temporary SWE before snowfall, used to check mass balance (kg m-2)
 real(rkind)                            :: tempSWE1                   ! temporary SWE after snowfall, used to check mass balance (kg m-2)
 real(rkind)                            :: xMassBalance               ! mass balance check (kg m-2)
 real(rkind),parameter                  :: massBalTol=1.e-8_rkind     ! tolerance for mass balance check (kg m-2)
 ! initialize error control
 err=0;! message="newsnwfall/"

 ! compute the new snowfall (kg m-2 s-1)
 newSnowfall = scalarThroughfallSnow + scalarCanopySnowUnloading

 ! early return if there is no snowfall
 if(newSnowfall < tiny(dt)) return

 ! compute depth of new snow
 newSnowDepth     = dt*(scalarThroughfallSnow/scalarNewSnowDensity + scalarCanopySnowUnloading/densityCanopySnow)  ! new snow depth (m)

 ! process special case of "snow without a layer"
 if(.not.snowLayers)then
  ! increment depth and water equivalent
  scalarSnowDepth = scalarSnowDepth + newSnowDepth
  scalarSWE       = scalarSWE + dt*newSnowfall

 ! add snow to the top layer (more typical case where snow layers already exist)
 else

  ! get SWE in the upper layer (used to check mass balance)
  tempSWE0 = (surfaceLayerVolFracIce*iden_ice + surfaceLayerVolFracLiq*iden_water)*surfaceLayerDepth

  ! get the total mass of liquid water and ice (kg m-2)
  totalMassIceSurfLayer  = iden_ice*surfaceLayerVolFracIce*surfaceLayerDepth + newSnowfall*dt
  ! get the total snow depth
  totalDepthSurfLayer    = surfaceLayerDepth + newSnowDepth
  ! compute the new temperature
  surfaceLayerTemp       = (surfaceLayerTemp*surfaceLayerDepth + scalarSnowfallTemp*newSnowDepth) / totalDepthSurfLayer
  ! compute new SWE for the upper layer (kg m-2)
  SWE = totalMassIceSurfLayer + iden_water*surfaceLayerVolFracLiq*surfaceLayerDepth
  ! compute new volumetric fraction of liquid water and ice (-)
  volFracWater = (SWE/totalDepthSurfLayer)/iden_water
  fracLiq      = fracliquid(surfaceLayerTemp,fc_param)                           ! fraction of liquid water
  surfaceLayerVolFracIce = (1._rkind - fracLiq)*volFracWater*(iden_water/iden_ice)  ! volumetric fraction of ice (-)
  surfaceLayerVolFracLiq =          fracLiq *volFracWater                        ! volumetric fraction of liquid water (-)
  ! update new layer depth (m)
  surfaceLayerDepth      = totalDepthSurfLayer

  ! get SWE in the upper layer (used to check mass balance)
  tempSWE1 = (surfaceLayerVolFracIce*iden_ice + surfaceLayerVolFracLiq*iden_water)*surfaceLayerDepth

  ! check SWE
  xMassBalance = tempSWE1 - (tempSWE0 + newSnowfall*dt)
  if (abs(xMassBalance) > massBalTol)then
  !  write(*,'(a,1x,f20.10)') 'SWE mass balance = ', xMassBalance
  !  message=trim(message)//'mass balance problem'
   err=20; return
  end if

 end if  ! if snow layers already exist

 end subroutine newsnwfall


end module volicePack_module
