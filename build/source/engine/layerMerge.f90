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

module layerMerge_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! access named variables for snow and soil
USE globalData,only:iname_snow        ! named variables for snow
USE globalData,only:iname_soil        ! named variables for soil

! access metadata
USE globalData,only:prog_meta,diag_meta,flux_meta,indx_meta   ! metadata

! physical constants
USE multiconst,only:&
                    iden_ice,       & ! intrinsic density of ice             (kg m-3)
                    iden_water        ! intrinsic density of liquid water    (kg m-3)

! access the derived types to define the data structures
USE data_types,only:&
                    var_d,            & ! data vector (rkind)
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength,      & ! data vector with variable length dimension (rkind)
                    model_options       ! defines the model decisions

! access named variables defining elements in the data structures
USE var_lookup,only:iLookPARAM,iLookPROG,iLookINDEX  ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS                   ! named variables for elements of the decision structure

! look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers, & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index

! provide access to external modules
USE var_derive_module,only:calcHeight,calcHeight_d ! module to calculate height at layer interfaces and layer mid-point

! privacy
implicit none
private
public::layerMerge,layerMerge_device

contains


 ! *****************************************************************************************************************
 ! public subroutine layerMerge: merge layers if the thickness is less than zmin
 ! *****************************************************************************************************************
 subroutine layerMerge(&
    nGRU, &
                       ! input/output: model data structures
                       tooMuchMelt,                 & ! intent(in):    flag to force merge of snow layers
                       model_decisions,             & ! intent(in):    model decisions
                       mpar_data,                   & ! intent(in):    model parameters
                       indx_data,                   & ! intent(inout): type of each layer
                       prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,                   & ! intent(inout): model fluxes for a local HRU
                       ! output
                       mergedLayers,                & ! intent(out): flag to denote that layers were merged
                       err,message)                   ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 use device_data_types
  USE enthalpyTemp_module,only:enthalpy2T_snwWat,T2enthalpy_snwWat,h_lookup,t_lookup ! convert temperature to liq+ice enthalpy for a snow layer

 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input/output: model data structures
 integer(i4b) :: nGRU
 logical(lgt),intent(in)         :: tooMuchMelt         ! flag to denote that ice is insufficient to support melt
 type(decisions_device),intent(in)  :: model_decisions  ! model decisions
 type(mpar_data_device),intent(in)    :: mpar_data           ! model parameters
 type(indx_data_device),intent(inout) :: indx_data           ! type of each layer
 type(prog_data_device),intent(inout) :: prog_data           ! model prognostic variables for a local HRU
 type(diag_data_device),intent(inout) :: diag_data           ! model diagnostic variables for a local HRU
 type(flux_data_device),intent(inout) :: flux_data           ! model flux variables
 ! output
 logical(lgt),intent(out),device        :: mergedLayers(:)        ! flag to denote that layers were merged
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! define local variables
 character(LEN=256)              :: cmessage            ! error message of downwind routine
!  real(rkind),dimension(5)           :: zminLayer           ! minimum layer depth in each layer (m)
 logical(lgt)                    :: removeLayer         ! flag to indicate need to remove a layer
 integer(i4b)                    :: nCheck              ! number of layers to check for combination
 integer(i4b)                    :: iSnow               ! index of snow layers (looping)
 integer(i4b)                    :: jSnow               ! index of snow layer identified for combination with iSnow
 integer(i4b)                    :: kSnow               ! index of the upper layer of the two layers identified for combination
  real(rkind),device,allocatable :: H_lookup_d(:), T_lookup_d(:)
 
    type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)
  H_lookup_d = h_lookup
  t_lookup_d = t_lookup

 ! initialize error control
 err=0; message="layerMerge/"
 ! --------------------------------------------------------------------------------------------------------
 ! associate variables to the data structures
 associate(&

 ! model decisions
 ix_snowLayers    => model_decisions%snowLayers, & ! decision for snow combination

 ! model parameters (control the depth of snow layers)
 zmin             => mpar_data%zmin_,                & ! minimum layer depth (m)
 zminLayer1       => mpar_data%zminLayer1_,          & ! minimum layer depth for the 1st (top) layer (m)
 zminLayer2       => mpar_data%zminLayer2_,          & ! minimum layer depth for the 2nd layer (m)
 zminLayer3       => mpar_data%zminLayer3_,          & ! minimum layer depth for the 3rd layer (m)
 zminLayer4       => mpar_data%zminLayer4_,          & ! minimum layer depth for the 4th layer (m)
 zminLayer5       => mpar_data%zminLayer5_,          & ! minimum layer depth for the 5th (bottom) layer (m)

 ! diagnostic scalar variables
 scalarSnowDepth  => prog_data%scalarSnowDepth,      & ! total snow depth (m)
 scalarSWE        => prog_data%scalarSWE             & ! SWE (kg m-2)

 ) ! end associate statement
 ! --------------------------------------------------------------------------------------------------------

 ! identify algorithmic control parameters to sub-divide and combine snow layers
!  zminLayer = (/zminLayer1, zminLayer2, zminLayer3, zminLayer4, zminLayer5/)



    call layerMerge_kernel<<<blocks,threads>>>(nGRU,ix_snowLayers,indx_data%nSnow,indx_data%nLayers_d,indx_data%nSoil, &
 prog_data%mLayerTemp,prog_data%mLayerVolFracIce,prog_data%mLayerVolFracLiq, &
prog_data%mLayerVolFracWat,prog_data%mLayerEnthalpy,prog_data%mLayerDepth, prog_data%mLayerHeight, prog_data%iLayerHeight, &
diag_data%mLayerVolHtCapBulk, diag_data%mLayerCm, diag_data%mLayerThermalC, diag_data%iLayerThermalC, &
diag_data%mLayerEnthTemp, diag_data%mLayerFracLiqSnow, diag_data%mLayerThetaResid, diag_data%mLayerPoreSpace, &
diag_data%mLayerMeltFreeze, diag_data%mLayerVolFracAir, diag_data%balanceLayerNrg, diag_data%balanceLayerMass,&
flux_data%data, &
flux_data%ixiLayerConductiveFlux_start,flux_data%ixiLayerConductiveFlux_end, &
flux_data%ixiLayerAdvectiveFlux_start,flux_data%ixiLayerAdvectiveFlux_end, &
flux_data%ixiLayerNrgFlux_start,flux_data%ixiLayerNrgFlux_end, flux_data%ixmLayerNrgFlux_start,flux_data%ixmLayerNrgFlux_end, &
flux_data%ixiLayerLiqFluxSnow_start,flux_data%ixiLayerLiqFluxSnow_end,flux_data%ixmLayerLiqFluxSnow_start,flux_data%ixmLayerLiqFluxSnow_end, &
    indx_data%layerType,indx_data%ixHydType, &
    indx_data%ixSnowSoilNrg, indx_data%ixSnowOnlyNrg, &
    indx_data%ixSnowSoilHyd, indx_data%ixSnowOnlyHyd, &
    indx_data%ixNrgLayer, indx_data%ixHydLayer,&
    tooMuchMelt, zMin,&
    zMinLayer1,zMinLayer2,zminLayer3,zminLayer4,zminLayer5, &
    mergedLayers,mpar_data%snowfrz_scale_,h_lookup_d,t_lookup_d, &
    prog_data%scalarSnowDepth, prog_data%scalarSWE)

 ! end association to variables in the data structure
 end associate

 end subroutine layerMerge

 attributes(global) subroutine layerMerge_kernel(nGRU,ix_snowLayers,nSnow,nLayers,nSoil, &
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
    mergedLayers,snowfrz_scale,h_lookup_d,t_lookup_d, &
    scalarSnowDepth, scalarSWE)
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
    real(rkind),intent(in) :: zMin(:)
    real(rkind),intent(in) :: zMinLayer1(:),zMinLayer2(:),zminLayer3(:),zminLayer4(:),zminLayer5(:)
    logical(lgt),intent(inout) :: mergedLayers(:)
    real(rkind),intent(in) :: snowfrz_scale(:)
    real(rkind),intent(inout) :: h_lookup_d(:), t_lookup_d(:)
    real(rkind),intent(inout) :: scalarSnowDepth(:), scalarSWE(:)

    integer(i4b),intent(in),value :: nGRU
       integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return
  call layerMerge_device(ix_snowLayers,nSnow(iGRU),nLayers(iGRU),nSoil, &
 mLayerTemp(:,iGRU),mLayerVolFracIce(:,iGRU),mLayerVolFracLiq(:,iGRU), &
mLayerVolFracWat(:,iGRU),mLayerEnthalpy(:,iGRU),mLayerDepth(:,iGRU), mLayerHeight(:,iGRU), iLayerHeight(:,iGRU), &
mLayerVolHtCapBulk(:,iGRU), mLayerCm(:,iGRU), mLayerThermalC(:,iGRU), iLayerThermalC(:,iGRU), &
mLayerEnthTemp(:,iGRU), mLayerFracLiqSnow(:,iGRU), mLayerThetaResid(:,iGRU), mLayerPoreSpace(:,iGRU), &
mLayerMeltFreeze(:,iGRU), mLayerVolFracAir(:,iGRU), balanceLayerNrg(:,iGRU), balanceLayerMass(:,iGRU),&
flux_data(iLayerConductiveFlux_start:iLayerConductiveFlux_end,iGRU), flux_data(iLayerAdvectiveFlux_start:iLayerAdvectiveFlux_end,iGRU), &
  flux_data(iLayerNrgFlux_start:iLayerNrgFlux_end,iGRU), flux_data(mLayerNrgFlux_start:mLayerNrgFlux_end,iGRU), &
  flux_data(iLayerLiqFluxSnow_start:iLayerLiqFluxSnow_end,iGRU),flux_data(mLayerLiqFluxSnow_start:mLayerLiqFluxSnow_end,iGRU),&
    layerType(:,iGRU),ixHydType(:,iGRU), &
    ixSnowSoilNrg(:,iGRU), ixSnowOnlyNrg(:,iGRU), &
    ixSnowSoilHyd(:,iGRU), ixSnowOnlyHyd(:,iGRU), &
    ixNrgLayer(:,iGRU), ixHydLayer(:,iGRU),&
    tooMuchMelt, zMin(iGRU),&
    zMinLayer1(iGRU),zMinLayer2(iGRU),zminLayer3(iGRU),zminLayer4(iGRU),zminLayer5(iGRU), &
    mergedLayers(iGRU),snowfrz_scale(iGRU),h_lookup_d,t_lookup_d, &
    scalarSnowDepth(iGRU), scalarSWE(iGRU))
end subroutine

attributes(device) subroutine layerMerge_device(ix_snowLayers,nSnow,nLayers,nSoil, &
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
    integer(i4b),intent(in) :: ix_snowLayers
    integer(i4b),intent(inout) :: nSnow,nLayers
    integer(i4b),intent(in) :: nSoil
    real(rkind),intent(inout) :: mLayerTemp(:)
    real(rkind),intent(inout) :: mLayerVolFracIce(:)
    real(rkind),intent(inout) :: mLayerVolFracLiq(:)
    real(rkind),intent(inout) :: mLayerVolFracWat(:)
    real(rkind),intent(inout) :: mLayerEnthalpy(:)
    real(rkind),intent(inout) :: mLayerDepth(:)
    real(rkind),intent(inout) :: mLayerHeight(0:)
    real(rkind),intent(inout) :: iLayerHeight(0:)
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
    logical(lgt),intent(in) :: tooMuchMelt
    real(rkind),intent(in) :: zMin
    real(rkind),intent(in) :: zMinLayer1,zMinLayer2,zminLayer3,zminLayer4,zminLayer5
    logical(lgt),intent(inout) :: mergedLayers
    real(rkind),intent(in) :: snowfrz_scale
    real(rkind),intent(inout) :: h_lookup_d(:), t_lookup_d(:)
    real(rkind),intent(inout) :: scalarSnowDepth, scalarSWE

integer(i4b) :: iSnow,jSnow,kSnow,nCheck
logical(lgt) :: removeLayer
mergedLayers = .false.
 kSnow=0 ! initialize first layer to test (top layer)
 do ! attempt to remove multiple layers in a single time step (continuous do loop with exit clause)

  ! special case of >5 layers: add an offset to use maximum threshold from layer above
  if(ix_snowLayers == rulesDependLayerIndex .and. nSnow > 5)then
   nCheck=5
  else
   nCheck=nSnow
  end if

  ! loop through snow layers
  do iSnow=kSnow+1,nCheck

   ! associate local variables with the information in the data structures
   ! NOTE: do this here, since the layer variables are re-defined
   ! check if the layer depth is less than the depth threshold
   select case(ix_snowLayers)
    case(sameRulesAllLayers);    removeLayer = (mLayerDepth(iSnow) < zmin)
    case(rulesDependLayerIndex); 
        select case(iSnow)
        case(1); removeLayer = (mLayerDepth(iSnow) < zminLayer1)
        case(2); removeLayer = (mLayerDepth(iSnow) < zminLayer2)
        case(3); removeLayer = (mLayerDepth(iSnow) < zminLayer3)
        case(4); removeLayer = (mLayerDepth(iSnow) < zminLayer4)
        case(5); removeLayer = (mLayerDepth(iSnow) < zminLayer5)
    end select
    ! case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
   end select ! (option to combine/sub-divide snow layers)

   ! check if we have too much melt
   ! NOTE: assume that this is the top snow layer; need more trickery to relax this assumption
   if(tooMuchMelt .and. iSnow==1) removeLayer=.true.

   ! check if need to remove a layer
   if(removeLayer)then

    ! flag that we modified a layer
    mergedLayers=.true.

    ! ***** handle special case of a single layer
    if(nSnow==1)then
     ! set the variables defining "snow without a layer"
     ! NOTE: ignoring cold content!!! Need to fix later...
     scalarSnowDepth = mLayerDepth(1)
     scalarSWE       = (mLayerVolFracIce(1)*iden_ice + mLayerVolFracLiq(1)*iden_water)*mLayerDepth(1)
     ! remove the top layer from all model variable vectors
     ! NOTE: nSnow-1 = 0, so routine removes layer #1
           call rmLyAllVars_prog(nSnow-1,nLayers, &
 mLayerTemp,mLayerVolFracIce,mLayerVolFracLiq, &
mLayerVolFracWat,mLayerEnthalpy,mLayerDepth, mLayerHeight, iLayerHeight)
call rmLyAllVars_diag(nSnow-1,nLayers,nSnow, &
mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass)
call rmLyAllVars_flux(nSnow-1,nLayers,nSnow,&
iLayerConductiveFlux,iLayerAdvectiveFlux,&
iLayerNrgFlux,mLayerNrgFlux,&
iLayerLiqFluxSnow,mLayerLiqFluxSnow)
  call rmLyAllVars_indx(nSnow-1,nSnow,nLayers,&
    layerType,ixHydType, &
    ixSnowSoilNrg, ixSnowOnlyNrg, &
    ixSnowSoilHyd, ixSnowOnlyHyd, &
    ixNrgLayer, ixHydLayer)
     ! update the total number of layers
     nSnow   = count(layerType==iname_snow)
    !  nSoil   = count(layerType==iname_soil)
     nLayers = nSnow + nSoil
     ! update coordinate variables
     call calcHeight_d(&
                 nLayers,layerType, &
                 mLayerDepth,mLayerHeight,iLayerHeight)
    !  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
     ! exit the do loop (no more snow layers to remove)
     return
    end if  ! (special case of 1 layer --> snow without a layer)

    ! ***** identify the layer to combine
    if(iSnow==1)then
     jSnow = iSnow+1  ! upper-most layer, combine with its lower neighbor
    elseif(iSnow==nSnow)then
     jSnow = nSnow-1  ! lower-most layer, combine with its upper neighbor
    else
     if(mLayerDepth(iSnow-1)<mLayerDepth(iSnow+1))then; jSnow = iSnow-1; else; jSnow = iSnow+1; end if
    end if

    ! ***** combine layers
    ! identify the layer closest to the surface
    kSnow=min(iSnow,jSnow)
    ! combine layer with identified neighbor
    call layer_combine_device(kSnow,nSoil,nLayers,nSnow, &
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
    ixNrgLayer, ixHydLayer,snowfrz_scale, &
    H_lookup_d, T_lookup_d)


    ! exit the loop to try again
    exit

   end if  ! (if layer is below the mass threshold)

   kSnow=iSnow ! ksnow is used for completion test, so include here

  end do ! (looping through snow layers)

  ! exit if finished
  if(kSnow==nCheck)exit

 end do ! continuous do

  ! handle special case of > 5 layers in the CLM option
 if(nSnow > 5 .and. ix_snowLayers == rulesDependLayerIndex)then
  ! flag that layers were merged
  mergedLayers=.true.
  ! initial check to ensure everything is wonderful in the universe
!   if(nSnow /= 6)then; err=5; message=trim(message)//'special case of > 5 layers: expect only six layers'; return; end if
  ! combine 5th layer with layer below
    call layer_combine_device(5,nSoil,nLayers,nSnow, &
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
    ixNrgLayer, ixHydLayer,snowfrz_scale, &
    H_lookup_d, T_lookup_d)
 end if

 ! check that there are no more than 5 layers in the CLM option
 if(ix_snowLayers == rulesDependLayerIndex)then
  if(nSnow > 5)then
!    message=trim(message)//'expect no more than 5 layers when combination/sub-division rules depend on the layer index (CLM option)'
!    err=20; return
  end if
 end if

end subroutine

 ! ***********************************************************************************************************
 ! private subroutine layer_combine: combine snow layers and re-compute model state variables
 ! ***********************************************************************************************************
 ! combines layer iSnow with iSnow+1
 ! ***********************************************************************************************************

attributes(device) subroutine layer_combine_device(iSnow,nSoil,nLayers,nSnow, &
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
    ixNrgLayer, ixHydLayer,snowfrz_scale, &
    H_lookup_d, T_lookup_d)
     USE snow_utils_module,only:fracliquid                            ! compute fraction of liquid water
 USE enthalpyTemp_module,only:enthalpy2T_snwWat,T2enthalpy_snwWat ! convert temperature to liq+ice enthalpy for a snow layer

    integer(i4b),intent(in) :: iSnow,nSoil
    integer(i4b),intent(inout) :: nLayers,nSnow
    real(rkind),intent(inout) :: mLayerTemp(:)
    real(rkind),intent(inout) :: mLayerVolFracIce(:)
    real(rkind),intent(inout) :: mLayerVolFracLiq(:)
    real(rkind),intent(inout) :: mLayerVolFracWat(:)
    real(rkind),intent(inout) :: mLayerEnthalpy(:)
    real(rkind),intent(inout) :: mLayerDepth(:)
    real(rkind),intent(inout) :: mLayerHeight(0:)
    real(rkind),intent(inout) :: iLayerHeight(0:)
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
    real(rkind),intent(inout) :: snowfrz_scale
    real(rkind),intent(in) :: H_lookup_d(:), T_lookup_d(:)

 real(rkind)                     :: massIce(2)               ! mass of ice in the two layers identified for combination (kg m-2)
 real(rkind)                     :: massLiq(2)               ! mass of liquid water in the two layers identified for combination (kg m-2)
 real(rkind)                     :: bulkDenWat(2)            ! bulk density if total water (liquid water plus ice) in the two layers identified for combination (kg m-3)
 real(rkind)                     :: cBulkDenWat              ! combined bulk density of total water (liquid water plus ice) in the two layers identified for combination (kg m-3)
 real(rkind)                     :: cTemp                    ! combined layer temperature
 real(rkind)                     :: cDepth                   ! combined layer depth
 real(rkind)                     :: cVolFracIce              ! combined layer volumetric fraction of ice
 real(rkind)                     :: cVolFracLiq              ! combined layer volumetric fraction of liquid water
 real(rkind)                     :: l1Enthalpy,l2Enthalpy    ! enthalpy in the two layers identified for combination (J m-3)
 real(rkind)                     :: cEnthalpy                ! combined layer enthalpy (J m-3)
 real(rkind)                     :: fLiq                     ! fraction of liquid water at the combined temperature cTemp
 integer(i4b) :: iLayer

 ! compute combined depth
 cDepth       = mLayerDepth(isnow) + mLayerDepth(isnow+1)

 ! compute mass of each layer (kg m-2)
 massIce(1:2) = iden_ice*mLayerVolFracIce(iSnow:iSnow+1)*mLayerDepth(iSnow:iSnow+1)
 massLiq(1:2) = iden_water*mLayerVolFracLiq(iSnow:iSnow+1)*mLayerDepth(iSnow:iSnow+1)

 ! compute bulk density of water (kg m-3)
 bulkDenWat(1:2) = (massIce(1:2) + massLiq(1:2))/mLayerDepth(iSnow:iSnow+1)
 cBulkDenWat     = (mLayerDepth(isnow)*bulkDenWat(1) + mLayerDepth(isnow+1)*bulkDenWat(2))/cDepth

 ! compute enthalpy for each layer (J m-3)
 l1Enthalpy = T2enthalpy_snwWat(mLayerTemp(iSnow),  BulkDenWat(1),snowfrz_scale)
 l2Enthalpy = T2enthalpy_snwWat(mLayerTemp(iSnow+1),BulkDenWat(2),snowfrz_scale)

 ! compute combined enthalpy (J m-3)
 cEnthalpy = (mLayerDepth(isnow)*l1Enthalpy + mLayerDepth(isnow+1)*l2Enthalpy)/cDepth

 ! convert enthalpy (J m-3) to temperature (K)
 call enthalpy2T_snwWat(cEnthalpy,cBulkDenWat,snowfrz_scale,cTemp,H_lookup_d,T_lookup_d)
!  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if

 ! test enthalpy conversion
!  if(abs(T2enthalpy_snwWat(cTemp,cBulkDenWat,snowfrz_scale)/cBulkDenWat - cEnthalpy/cBulkDenWat) > eTol)then
!   write(*,'(a,1x,f12.5,1x,2(e20.10,1x))') 'enthalpy test', cBulkDenWat, T2enthalpy_snwWat(cTemp,cBulkDenWat,snowfrz_scale)/cBulkDenWat, cEnthalpy/cBulkDenWat
!   message=trim(message)//'problem with enthalpy-->temperature conversion'
!   err=20; return
!  end if

 ! check temperature is within the two temperatures
 ! NOTE: use tolerance, for cases of merging a layer that has just been split
!  if(cTemp > max(mLayerTemp(iSnow),mLayerTemp(iSnow+1))+eTol)then; err=20; message=trim(message)//'merged temperature > max(temp1,temp2)'; return; end if
!  if(cTemp < min(mLayerTemp(iSnow),mLayerTemp(iSnow+1))-eTol)then; err=20; message=trim(message)//'merged temperature < min(temp1,temp2)'; return; end if

 ! compute volumetric fraction of liquid water
 fLiq = fracLiquid(cTemp,snowfrz_scale)

 ! compute volumetric fraction of ice and liquid water
 cVolFracLiq =          fLiq *cBulkDenWat/iden_water
 cVolFracIce = (1._rkind - fLiq)*cBulkDenWat/iden_ice

      call rmLyAllVars_prog(iSnow,nLayers, &
 mLayerTemp,mLayerVolFracIce,mLayerVolFracLiq, &
mLayerVolFracWat,mLayerEnthalpy,mLayerDepth, mLayerHeight, iLayerHeight)
call rmLyAllVars_diag(iSnow,nLayers,nSnow, &
mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass)
call rmLyAllVars_flux(iSnow,nLayers,nSnow,&
iLayerConductiveFlux,iLayerAdvectiveFlux,&
iLayerNrgFlux,mLayerNrgFlux,&
iLayerLiqFluxSnow,mLayerLiqFluxSnow)
  call rmLyAllVars_indx(iSnow,nSnow,nLayers,&
    layerType,ixHydType, &
    ixSnowSoilNrg, ixSnowOnlyNrg, &
    ixSnowSoilHyd, ixSnowOnlyHyd, &
    ixNrgLayer, ixHydLayer)

 ! define the combined layer as snow
 layerType(iSnow) = iname_snow

 ! save the number of layers in the data structures
 nSnow   = count(layerType==iname_snow)
 nLayers = nSnow + nSoil

  ! ***** put state variables for the combined layer in the appropriate place
 mLayerTemp(iSnow)       = cTemp
 mLayerDepth(iSnow)      = cDepth
 mLayerVolFracIce(iSnow) = cVolFracIce
 mLayerVolFracLiq(iSnow) = cVolFracLiq

  ! ***** adjust coordinate variables
 call calcHeight_d(&
                 nLayers,layerType, &
                 mLayerDepth,mLayerHeight,iLayerHeight)


end subroutine
 attributes(device) subroutine rmLyAllVars_prog(iSnow,nLayers, &
 mLayerTemp,mLayerVolFracIce,mLayerVolFracLiq, &
mLayerVolFracWat,mLayerEnthalpy,mLayerDepth, mLayerHeight, iLayerHeight)
    integer(i4b),intent(in) :: iSnow,nLayers
    real(rkind),intent(inout) :: mLayerTemp(:)
    real(rkind),intent(inout) :: mLayerVolFracIce(:)
    real(rkind),intent(inout) :: mLayerVolFracLiq(:)
    real(rkind),intent(inout) :: mLayerVolFracWat(:)
    real(rkind),intent(inout) :: mLayerEnthalpy(:)
    real(rkind),intent(inout) :: mLayerDepth(:)
    real(rkind),intent(inout) :: mLayerHeight(0:)
    real(rkind),intent(inout) :: iLayerHeight(0:)
    call rmLyAllVars_real(iSnow,1,nLayers,mLayerTemp)
    call rmLyAllVars_real(iSnow,1,nLayers,mLayerVolFracIce)
    call rmLyAllVars_real(iSnow,1,nLayers,mLayerVolFracLiq)
    call rmLyAllVars_real(iSnow,1,nLayers,mLayerVolFracWat)
    call rmLyAllVars_real(iSnow,1,nLayers,mLayerEnthalpy)
    call rmLyAllVars_real(iSnow,1,nLayers,mLayerDepth)
    call rmLyAllVars_real(iSnow,0,nLayers,mLayerHeight)
    call rmLyAllVars_real(iSnow,0,nLayers,iLayerHeight)
end subroutine

attributes(device) subroutine rmLyAllVars_diag(iSnow,nLayers,nSnow, &
mLayerVolHtCapBulk, mLayerCm, mLayerThermalC, iLayerThermalC, &
mLayerEnthTemp, mLayerFracLiqSnow, mLayerThetaResid, mLayerPoreSpace, &
mLayerMeltFreeze, mLayerVolFracAir, balanceLayerNrg, balanceLayerMass)
integer(i4b),intent(in) :: iSnow, nLayers, nSnow
  real(rkind),intent(inout) :: mLayerVolHtCapBulk(:), mLayerCm(:), mLayerThermalC(:), iLayerThermalC(0:)
  real(rkind),intent(inout) :: mLayerEnthTemp(:), mLayerFracLiqSnow(:), mLayerThetaResid(:), mLayerPoreSpace(:)
  real(rkind),intent(inout) :: mLayerMeltFreeze(:), mLayerVolFracAir(:), balanceLayerNrg(:), balanceLayerMass(:)
  call rmLyAllVars_real(iSnow,1,nLayers,mLayerVolHtCapBulk)
  call rmLyAllVars_real(iSnow,1,nLayers,mLayerCm)
  call rmLyAllVars_real(iSnow,1,nLayers,mLayerThermalC)
  call rmLyAllVars_real(iSnow,0,nLayers,iLayerThermalC)
  call rmLyAllVars_real(iSnow,1,nLayers,mLayerEnthTemp)
  call rmLyAllVars_real(iSnow,1,nSnow,mLayerFracLiqSnow)
  call rmLyAllVars_real(iSnow,1,nSnow,mLayerThetaResid)
  call rmLyAllVars_real(iSnow,1,nSnow,mLayerPoreSpace)
  call rmLyAllVars_real(iSnow,1,nLayers,mLayerMeltFreeze)
  call rmLyAllVars_real(iSnow,1,nLayers,mLayerVolFracAir)
  call rmLyAllVars_real(iSnow,1,nLayers,balanceLayerNrg)
  call rmLyAllVars_real(iSnow,1,nLayers,balanceLayerMass)
end subroutine

attributes(device) subroutine rmLyAllVars_flux(iSnow,nLayers,nSnow,&
iLayerConductiveFlux,iLayerAdvectiveFlux,&
iLayerNrgFlux,mLayerNrgFlux,&
iLayerLiqFluxSnow,mLayerLiqFluxSnow)
    integer(i4b),intent(in) :: iSnow,nLayers,nSnow
    real(rkind),intent(inout) :: iLayerConductiveFlux(0:)
    real(rkind),intent(inout) :: iLayerAdvectiveFlux(0:)
    real(rkind),intent(inout) :: iLayerNrgFlux(0:)
    real(rkind),intent(inout) :: mLayerNrgFlux(:)
    real(rkind),intent(inout) :: iLayerLiqFluxSnow(0:)
    real(rkind),intent(inout) :: mLayerLiqFluxSnow(:)

    call rmLyAllVars_real(iSnow,0,nLayers,iLayerConductiveFlux)
    call rmLyAllVars_real(iSnow,0,nLayers,iLayerAdvectiveFlux)
    call rmLyAllVars_real(iSnow,0,nLayers,iLayerNrgFlux)
    call rmLyAllVars_real(iSnow,1,nLayers,mLayerNrgFlux)
    call rmLyAllVars_real(iSnow,0,nSnow,iLayerLiqFluxSnow)
    call rmLyAllVars_real(iSnow,1,nSnow,mLayerLiqFluxSnow)
end subroutine

attributes(device) subroutine rmLyAllVars_indx(iSnow,nSnow,nLayers,&
    layerType,ixHydType, &
    ixSnowSoilNrg, ixSnowOnlyNrg, &
    ixSnowSoilHyd, ixSnowOnlyHyd, &
    ixNrgLayer, ixHydLayer)
    integer(i4b),intent(in) :: iSnow,nSnow,nLayers
    integer(i4b),intent(inout) :: layerType(:), ixHydType(:)
    integer(i4b),intent(inout) :: ixSnowSoilNrg(:), ixSnowOnlyNrg(:)
    integer(i4b),intent(inout) :: ixSnowSoilHyd(:), ixSnowOnlyHyd(:)
    integer(i4b),intent(inout) :: ixNrgLayer(:), ixHydLayer(:)
    
  call rmLyAllVars_int(iSnow,1,nLayers,layerType)
  call rmLyAllVars_int(iSnow,1,nLayers,ixHydType)
  call rmLyAllVars_int(iSnow,1,nLayers,ixSnowSoilNrg)
  call rmLyAllVars_int(iSnow,1,nSnow,ixSnowOnlyNrg)
  call rmLyAllVars_int(iSnow,1,nLayers,ixSnowSoilHyd)
  call rmLyAllVars_int(iSnow,1,nSnow,ixSnowOnlyHyd)
  call rmLyAllVars_int(iSnow,1,nLayers,ixNrgLayer)
  call rmLyAllVars_int(iSnow,1,nLayers,ixHydLayer)
!   call rmLyAllVars_int(iSnow,1,nLayers,ixLayerState)
!   call rmLyAllVars_int(iSnow,1,nLayers,ixLayerActive)
end subroutine


 attributes(device) subroutine rmLyAllVars_real(iSnow, ix_lower,ix_upper, vec_rkind)
    real(rkind),intent(inout) :: vec_rkind(ix_lower:)
    integer(i4b),intent(in) :: iSnow, ix_upper, ix_lower
    integer(i4b) :: iLayer

    if(iSnow>=ix_lower)  Vec_rkind(iSnow)              = realMissing ! set merged layer to missing (fill in later)
    if(iSnow>ix_lower)   Vec_rkind(ix_lower:iSnow-1)   = Vec_rkind(ix_lower:iSnow-1)
    if(iSnow+1<ix_upper) then
        do iLayer=iSnow+1,ix_upper-1 ! skip iSnow+1
            Vec_rkind(iLayer) = Vec_rkind(iLayer+1)  
        end do
    end if
    vec_rkind(ix_upper) = realMissing
 end subroutine
 attributes(device) subroutine rmLyAllVars_int(iSnow, ix_lower,ix_upper, vec)
    integer(i4b),intent(inout) :: vec(ix_lower:)
    integer(i4b),intent(in) :: iSnow, ix_upper, ix_lower
    integer(i4b) :: iLayer

    if(iSnow>=ix_lower)  Vec(iSnow)              = integerMissing ! set merged layer to missing (fill in later)
    if(iSnow>ix_lower)   Vec(ix_lower:iSnow-1)   = Vec(ix_lower:iSnow-1)
    if(iSnow+1<ix_upper) then
        do iLayer=iSnow+1,ix_upper-1 ! skip iSnow+1
            Vec(iLayer) = Vec(iLayer+1)  
        end do
    end if
    Vec(ix_upper) = integerMissing
 end subroutine

 ! ***********************************************************************************************************
 ! private subroutine rmLyAllVars: reduce the length of the vectors in data structures
 ! ***********************************************************************************************************
 ! removes layer "iSnow+1" and sets layer "iSnow" to a missing value
 ! (layer "iSnow" will be filled with a combined layer later)
 ! ***********************************************************************************************************
 subroutine rmLyAllVars(dataStruct,metaStruct,iSnow,nSnow,nLayers,err,message)
 USE var_lookup,only:iLookVarType                 ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName       ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc            ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength      ! data vectors with variable length dimension
 USE data_types,only:var_info                     ! metadata structure
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 class(*),intent(inout)          :: dataStruct     ! data structure
 type(var_info),intent(in)       :: metaStruct(:)  ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in)         :: iSnow          ! new layer
 integer(i4b),intent(in)         :: nSnow,nLayers  ! number of snow layers, total number of layers
 ! output: error control
 integer(i4b),intent(out)        :: err            ! error code
 character(*),intent(out)        :: message        ! error message
 ! locals
 integer(i4b)                    :: ivar           ! variable index
 integer(i4b)                    :: ix_lower       ! lower bound of the vector
 integer(i4b)                    :: ix_upper       ! upper bound of the vector
 real(rkind),allocatable            :: tempVec_rkind(:)  ! temporary vector (double precision)
 integer(i4b),allocatable        :: tempVec_i4b(:) ! temporary vector (integer)
 character(LEN=256)              :: cmessage       ! error message of downwind routine
 ! initialize error control
 err=0; message="rmLyAllVars/"

 ! check dimensions
 select type(dataStruct)
  type is (var_dlength); if(size(dataStruct%var) /= size(metaStruct)) err=20
  type is (var_ilength); if(size(dataStruct%var) /= size(metaStruct)) err=20
  class default; err=20; message=trim(message)//'unable to identify the data type'; return
 end select
 if(err/=0)then; message=trim(message)//'dimensions of data structure and metadata structures do not match'; return; end if

 ! ***** loop through model variables and remove one layer
 do ivar=1,size(metaStruct)

  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=nSnow
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=nSnow
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nLayers
   case default; cycle  ! no need to remove soil layers or scalar variables
  end select

  ! remove layers
  select type(dataStruct)

   ! ** double precision
   type is (var_dlength)
    ! check allocated
    if(.not.allocated(dataStruct%var(ivar)%dat))then; err=20; message='data vector is not allocated'; return; end if
    ! allocate the temporary vector
    allocate(tempVec_rkind(ix_lower:ix_upper-1), stat=err)
    if(err/=0)then; err=20; message=trim(message)//'unable to allocate temporary vector'; return; end if
    ! copy elements across to the temporary vector
    if(iSnow>=ix_lower)  tempVec_rkind(iSnow)              = realMissing ! set merged layer to missing (fill in later)
    if(iSnow>ix_lower)   tempVec_rkind(ix_lower:iSnow-1)   = dataStruct%var(ivar)%dat(ix_lower:iSnow-1)
    if(iSnow+1<ix_upper) tempVec_rkind(iSnow+1:ix_upper-1) = dataStruct%var(ivar)%dat(iSnow+2:ix_upper)  ! skip iSnow+1
    ! deallocate the data vector: strictly not necessary, but include to be safe
    deallocate(dataStruct%var(ivar)%dat,stat=err)
    if(err/=0)then; err=20; message='problem deallocating data vector'; return; end if
    ! create the new data structure using the temporary vector as the source
    call cloneStruc(dataStruct%var(ivar)%dat, ix_lower, source=tempVec_rkind, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! deallocate the temporary data vector: strictly not necessary, but include to be safe
    deallocate(tempVec_rkind,stat=err)
    if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

   ! ** integer
   type is (var_ilength)
    ! check allocated
    if(.not.allocated(dataStruct%var(ivar)%dat))then; err=20; message='data vector is not allocated'; return; end if
    ! allocate the temporary vector
    allocate(tempVec_i4b(ix_lower:ix_upper-1), stat=err)
    if(err/=0)then; err=20; message=trim(message)//'unable to allocate temporary vector'; return; end if
    ! copy elements across to the temporary vector
    if(iSnow>=ix_lower)  tempVec_i4b(iSnow)              = integerMissing ! set merged layer to missing (fill in later)
    if(iSnow>ix_lower)   tempVec_i4b(ix_lower:iSnow-1)   = dataStruct%var(ivar)%dat(ix_lower:iSnow-1)
    if(iSnow+1<ix_upper) tempVec_i4b(iSnow+1:ix_upper-1) = dataStruct%var(ivar)%dat(iSnow+2:ix_upper)  ! skip iSnow+1
    ! deallocate the data vector: strictly not necessary, but include to be safe
    deallocate(dataStruct%var(ivar)%dat,stat=err)
    if(err/=0)then; err=20; message='problem deallocating data vector'; return; end if
    ! create the new data structure using the temporary vector as the source
    call cloneStruc(dataStruct%var(ivar)%dat, ix_lower, source=tempVec_i4b, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! deallocate the temporary data vector: strictly not necessary, but include to be safe
    deallocate(tempVec_i4b,stat=err)
    if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

   ! check that we found the data type
   class default; err=20; message=trim(message)//'unable to identify the data type'; return

  end select ! dependence on data types

 end do  ! looping through variables

 end subroutine rmLyAllVars

end module layerMerge_module
