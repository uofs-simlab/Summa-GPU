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

module computJacob_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_ilength,         & ! data vector with variable length dimension (i4b)
                    var_dlength,         & ! data vector with variable length dimension (rkind)
                    model_options,       & ! defines the model decisions
                    in_type_computJacob, & ! class for computJacob arguments
                    out_type_computJacob   ! class for computJacob arguments

! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! named variables to describe the state variable type
USE globalData,only:iname_watLayer   ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:maxVolIceContent ! snow maximum volumetric ice content to store water (-)

! access named variables to describe the form and structure of the matrices used in the numerical solver
USE globalData,only: ku             ! number of super-diagonal bands, assume ku>=3
USE globalData,only: kl             ! number of sub-diagonal bands, assume kl>=4
USE globalData,only: ixDiag         ! index for the diagonal band
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:       &
 qbaseTopmodel,                   & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                       & ! a big bucket (lumped aquifer model)
 noExplicit                         ! no explicit groundwater parameterization

implicit none
private
public::fluxJacAdd
logical::fullMatrix
contains


! ***********************************************************************************************************
! public subroutine to compute flux parts of the Jacobian that are shared between IDA and BE
! ***********************************************************************************************************
attributes(device) subroutine fluxJacAdd(&
                      ! input: model control
                      ! passed_fullMatrix,          & ! intent(in):    flag to indicate if the matrix is full (true) or banded (false)
                      dt,                         & ! intent(in):    length of the time step (seconds)
                      nSnow_,                      & ! intent(in):    number of snow layers
                      nSoil,                      & ! intent(in):    number of soil layers
                      nLayers_,                    & ! intent(in):    total number of layers
                      computeVegFlux_,             & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      computeBaseflow,            & ! intent(in):    flag to indicate if we need to compute baseflow
                      ! input: data structures
                      ! indx_data,                  & ! intent(in):    index data
                      ixCasNrg_, ixVegNrg_, ixVegHyd_, ixTopNrg_, ixTopHyd_, ixAqWat_, &
                      ixSnowSoilNrg_, ixSnowOnlyNrg_, ixSoilOnlyNrg_, &
                      ixSnowSoilHyd_, ixSnowOnlyHyd_, ixSoilOnlyHyd_, &
                      ixHydType_, &
                      ! prog_data,                  & ! intent(in):    model prognostic variables for a local HRU
                      mLayerDepth_, mLayerVolFracIce_, &
                      ! diag_data,                  & ! intent(in):    model diagnostic variables for a local HRU
                      mLayerFracLiqSnow_, &
                      scalarFracLiqVeg_, scalarSoilControl_, canopyDepth_, &
                      ! deriv_data,                 & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                      dCanairNetFlux_dCanairTemp_,dCanairNetFlux_dCanopyTemp_,dCanairNetFlux_dGroundTemp_,dCanopyNetFlux_dCanairTemp_,dCanopyNetFlux_dCanopyTemp_,dCanopyNetFlux_dGroundTemp_,dGroundNetFlux_dCanairTemp_,dGroundNetFlux_dCanopyTemp_,dGroundNetFlux_dCanWat_, &
                      dCanopyEvaporation_dTCanair_,dCanopyEvaporation_dTCanopy_,dCanopyEvaporation_dTGround_,dCanopyEvaporation_dCanWat_,dGroundEvaporation_dTCanair_,dGroundEvaporation_dTCanopy_,dGroundEvaporation_dTGround_,dGroundEvaporation_dCanWat_, &
                      scalarCanopyLiqDeriv_, dCanLiq_dTcanopy_, &
                      dNrgFlux_dTempAbove_, dNrgFlux_dTempBelow_, &
                      dNrgFlux_dWatAbove_, dNrgFlux_dWatBelow_, &
                      mLayerdTrans_dTCanair_,mLayerdTrans_dTCanopy_,mLayerdTrans_dTGround_,mLayerdTrans_dCanWat_,&
                      dAquiferTrans_dTCanair_,dAquiferTrans_dTCanopy_,dAquiferTrans_dTGround_,dAquiferTrans_dCanWat_, &
                      iLayerLiqFluxSnowDeriv_, &
                      dq_dHydStateAbove_,dq_dHydStateBelow_,dq_dHydStateLayerSurfVec_,&
                      dBaseflow_dAquifer_, &
                      dq_dNrgStateAbove_,dq_dNrgStateBelow_,dq_dNrgStateLayerSurfVec_,&
                      mLayerdTheta_dTk_, &
                      dBaseflow_dMatric_,          & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
                      ! input-output: Jacobian and its diagonal
                      dMat,                       & ! intent(in):    diagonal of the Jacobian matrix
                      aJac_,                       & ! intent(inout): Jacobian matrix with flux terms added
                      ! output: error control
                      err)                  ! intent(out):   error code and error message
  ! -----------------------------------------------------------------------------------------------------------------
  use initialize_device,only:get_iGRU
  implicit none
  ! input: model control
  ! logical(lgt),intent(in)              :: passed_fullMatrix          ! flag to indicate if the matrix is full (true) or banded (false)
  real(rkind),intent(in)               :: dt                         ! length of the time step (seconds)
  integer(i4b),intent(in)              :: nSnow_(:)                      ! number of snow layers
  integer(i4b),intent(in)              :: nSoil                      ! number of soil layers
  integer(i4b),intent(in)              :: nLayers_(:)                    ! total number of layers in the snow and soil domains
  logical(lgt),intent(in)              :: computeVegFlux_(:)             ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)              :: computeBaseflow            ! flag to indicate if computing baseflow
  ! input: data structures
  ! type(var_ilength),intent(in)         :: indx_data                  ! indices defining model states and layers
  integer(i4b),intent(in) :: ixCasNrg_(:), ixVegNrg_(:), ixVegHyd_(:), ixTopNrg_(:), ixTopHyd_(:), ixAqWat_(:)
  integer(i4b),intent(in) :: ixSnowSoilNrg_(:,:), ixSnowOnlyNrg_(:,:), ixSoilOnlyNrg_(:,:)
  integer(i4b),intent(in) :: ixSnowSoilHyd_(:,:), ixSnowOnlyHyd_(:,:), ixSoilOnlyHyd_(:,:)
  integer(i4b),intent(in) :: ixHydType_(:,:)
  ! type(var_dlength),intent(in)         :: prog_data                  ! prognostic variables for a local HRU
  real(rkind),intent(in) :: mLayerDepth_(:,:), mLayerVolFracIce_(:,:)
  ! type(var_dlength),intent(in)         :: diag_data                  ! diagnostic variables for a local HRU
  real(rkind),intent(in) :: mLayerFracLiqSnow_(:,:)
  real(rkind),intent(in) :: scalarFracLiqVeg_(:), scalarSoilControl_(:), canopyDepth_(:)
  ! type(var_dlength),intent(in)         :: deriv_data                 ! derivatives in model fluxes w.r.t. relevant state variables
  real(rkind),intent(in) :: dCanairNetFlux_dCanairTemp_(:),dCanairNetFlux_dCanopyTemp_(:),dCanairNetFlux_dGroundTemp_(:),dCanopyNetFlux_dCanairTemp_(:),dCanopyNetFlux_dCanopyTemp_(:),dCanopyNetFlux_dGroundTemp_(:),dGroundNetFlux_dCanairTemp_(:),dGroundNetFlux_dCanopyTemp_(:),dGroundNetFlux_dCanWat_(:)
  real(rkind),intent(in) :: dCanopyEvaporation_dTCanair_(:),dCanopyEvaporation_dTCanopy_(:),dCanopyEvaporation_dTGround_(:),dCanopyEvaporation_dCanWat_(:),dGroundEvaporation_dTCanair_(:),dGroundEvaporation_dTCanopy_(:),dGroundEvaporation_dTGround_(:),dGroundEvaporation_dCanWat_(:)
  real(rkind),intent(in) :: scalarCanopyLiqDeriv_(:), dCanLiq_dTcanopy_(:)
  real(rkind),intent(in) :: dNrgFlux_dTempAbove_(0:,:), dNrgFlux_dTempBelow_(0:,:)
  real(rkind),intent(in) :: dNrgFlux_dWatAbove_(0:,:), dNrgFlux_dWatBelow_(0:,:)
  real(rkind),intent(in) :: mLayerdTrans_dTCanair_(:,:),mLayerdTrans_dTCanopy_(:,:),mLayerdTrans_dTGround_(:,:),mLayerdTrans_dCanWat_(:,:)
  real(rkind),intent(in) :: dAquiferTrans_dTCanair_(:),dAquiferTrans_dTCanopy_(:),dAquiferTrans_dTGround_(:),dAquiferTrans_dCanWat_(:)
  real(rkind),intent(in) :: iLayerLiqFluxSnowDeriv_(0:,:)
  real(rkind),intent(in) :: dq_dHydStateAbove_(0:,:),dq_dHydStateBelow_(0:,:),dq_dHydStateLayerSurfVec_(0:,:)
  real(rkind),intent(in) :: dBaseflow_dAquifer_(:)
  real(rkind),intent(in) :: dq_dNrgStateAbove_(0:,:),dq_dNrgStateBelow_(0:,:),dq_dNrgStateLayerSurfVec_(0:,:)
  real(rkind),intent(in) :: mLayerdTheta_dTk_(:,:)
  real(rkind),intent(in)               :: dBaseflow_dMatric_(:,:,:)     ! derivative in baseflow w.r.t. matric head (s-1)
  ! input-output: Jacobian and its diagonal
  real(rkind),intent(in)               :: dMat(:,:)                    ! diagonal of the Jacobian matrix
  real(rkind),intent(inout)            :: aJac_(:,:,:)                  ! Jacobian matrix with flux terms added
  ! output variables
  integer(i4b),intent(out)             :: err                        ! error code
  ! character(*),intent(out)             :: message                    ! error message
  ! --------------------------------------------------------------
  ! * local variables
  ! --------------------------------------------------------------
  ! indices of model state variables
  integer(i4b)                         :: qState                     ! index of cross-derivative state variable for baseflow
  integer(i4b)                         :: nrgState                   ! energy state variable
  integer(i4b)                         :: watState                   ! hydrology state variable
  ! indices of model layers
  integer(i4b)                         :: iLayer,pLayer              ! index of model layer
  integer(i4b)                         :: jLayer                     ! index of model layer within the full state vector (hydrology)
  integer(i4b)                         :: denseLimit                 ! index of the limiting dense layer
  ! conversion factors
  real(rkind)                          :: convLiq2tot                ! factor to convert liquid water derivative to total water derivative
  integer(i4b) :: iGRU
  iGRU = get_iGRU()
  ! --------------------------------------------------------------
  ! associate variables from data structures
    ! --------------------------------------------------------------
    ! initialize error control
    err=0!; message='fluxJacAdd/'
    ! fullMatrix = passed_fullMatrix ! local copy of the flag to indicate if the matrix is full (true) or banded (false)
    ! -----
    ! * energy and liquid fluxes over vegetation...
    ! ---------------------------------------------
    if(computeVegFlux_(iGRU))then ! (derivatives only defined when vegetation protrudes over the surface)

      ! * energy fluxes with the canopy water
      if(ixVegHyd_(iGRU)/=integerMissing)then

        ! * cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
        if(ixCasNrg_(iGRU)/=integerMissing) aJac_(ixVegHyd_(iGRU),ixCasNrg_(iGRU),iGRU) = -dCanopyEvaporation_dTCanair_(iGRU)*dt
        ! dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy is the derivative in throughfall and canopy drainage with canopy temperature
        if(ixVegNrg_(iGRU)/=integerMissing) aJac_(ixVegHyd_(iGRU),ixVegNrg_(iGRU),iGRU) = -dCanopyEvaporation_dTCanopy_(iGRU)*dt + dt*scalarCanopyLiqDeriv_(iGRU)*dCanLiq_dTcanopy_(iGRU)
        ! * liquid water fluxes for vegetation canopy (-), dt*scalarFracLiqVeg*scalarCanopyLiqDeriv is the derivative in throughfall and canopy drainage with canopy water
                                     aJac_(ixVegHyd_(iGRU),ixVegHyd_(iGRU),iGRU) = -scalarFracLiqVeg_(iGRU)*(dCanopyEvaporation_dCanWat_(iGRU) - scalarCanopyLiqDeriv_(iGRU))*dt + aJac_(ixVegHyd_(iGRU),ixVegHyd_(iGRU),iGRU)
        if(ixTopNrg_(iGRU)/=integerMissing) aJac_(ixVegHyd_(iGRU),ixTopNrg_(iGRU),iGRU) = -dCanopyEvaporation_dTGround_(iGRU)*dt

        ! * cross-derivative terms w.r.t. canopy water (kg-1 m2)
        if(nSnow_(iGRU)>0)then
          if(ixTopHyd_(iGRU)/=integerMissing) aJac_(ixTopHyd_(iGRU),ixVegHyd_(iGRU),iGRU) = (dt/mLayerDepth_(1,iGRU))*(-scalarFracLiqVeg_(iGRU)*scalarCanopyLiqDeriv_(iGRU))/iden_water
        else
          if(ixTopHyd_(iGRU)/=integerMissing) aJac_(ixTopHyd_(iGRU),ixVegHyd_(iGRU),iGRU) = (dt/mLayerDepth_(1,iGRU))*(-scalarSoilControl_(iGRU)*scalarFracLiqVeg_(iGRU)*scalarCanopyLiqDeriv_(iGRU))/iden_water
        endif

        ! * cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
        if(ixTopNrg_(iGRU)/=integerMissing) aJac_(ixTopNrg_(iGRU),ixVegHyd_(iGRU),iGRU) = (dt/mLayerDepth_(1,iGRU))*(-dGroundNetFlux_dCanWat_(iGRU))
      endif

      ! * cross-derivative terms w.r.t. canopy temperature (K-1)
      if(ixVegNrg_(iGRU)/=integerMissing)then
        if(nSnow_(iGRU)>0)then
          if(ixTopHyd_(iGRU)/=integerMissing) aJac_(ixTopHyd_(iGRU),ixVegNrg_(iGRU),iGRU) = (dt/mLayerDepth_(1,iGRU))*(-scalarCanopyLiqDeriv_(iGRU)*dCanLiq_dTcanopy_(iGRU))/iden_water
        else
          if(ixTopHyd_(iGRU)/=integerMissing) aJac_(ixTopHyd_(iGRU),ixVegNrg_(iGRU),iGRU) = (dt/mLayerDepth_(1,iGRU))*(-scalarSoilControl_(iGRU)*scalarCanopyLiqDeriv_(iGRU)*dCanLiq_dTcanopy_(iGRU))/iden_water
        endif
      endif

      ! * energy fluxes with the canopy air space (J m-3 K-1)
      if(ixCasNrg_(iGRU)/=integerMissing)then
                                     aJac_(ixCasNrg_(iGRU),ixCasNrg_(iGRU),iGRU) = (dt/canopyDepth_(iGRU))*(-dCanairNetFlux_dCanairTemp_(iGRU)) + aJac_(ixCasNrg_(iGRU),ixCasNrg_(iGRU),iGRU)
        if(ixVegNrg_(iGRU)/=integerMissing) aJac_(ixCasNrg_(iGRU),ixVegNrg_(iGRU),iGRU) = (dt/canopyDepth_(iGRU))*(-dCanairNetFlux_dCanopyTemp_(iGRU))
        if(ixTopNrg_(iGRU)/=integerMissing) aJac_(ixCasNrg_(iGRU),ixTopNrg_(iGRU),iGRU) = (dt/canopyDepth_(iGRU))*(-dCanairNetFlux_dGroundTemp_(iGRU))
      endif

      ! * energy fluxes with the vegetation canopy (J m-3 K-1)
      if(ixVegNrg_(iGRU)/=integerMissing)then
        if(ixCasNrg_(iGRU)/=integerMissing) aJac_(ixVegNrg_(iGRU),ixCasNrg_(iGRU),iGRU) = (dt/canopyDepth_(iGRU))*(-dCanopyNetFlux_dCanairTemp_(iGRU))
                                     aJac_(ixVegNrg_(iGRU),ixVegNrg_(iGRU),iGRU) = (dt/canopyDepth_(iGRU))*(-dCanopyNetFlux_dCanopyTemp_(iGRU)) + aJac_(ixVegNrg_(iGRU),ixVegNrg_(iGRU),iGRU)
        if(ixTopNrg_(iGRU)/=integerMissing) aJac_(ixVegNrg_(iGRU),ixTopNrg_(iGRU),iGRU) = (dt/canopyDepth_(iGRU))*(-dCanopyNetFlux_dGroundTemp_(iGRU))
      endif

      ! * energy fluxes with the surface (J m-3 K-1)
      if(ixTopNrg_(iGRU)/=integerMissing)then
        if(ixCasNrg_(iGRU)/=integerMissing) aJac_(ixTopNrg_(iGRU),ixCasNrg_(iGRU),iGRU) = (dt/mLayerDepth_(1,iGRU))*(-dGroundNetFlux_dCanairTemp_(iGRU))
        if(ixVegNrg_(iGRU)/=integerMissing) aJac_(ixTopNrg_(iGRU),ixVegNrg_(iGRU),iGRU) = (dt/mLayerDepth_(1,iGRU))*(-dGroundNetFlux_dCanopyTemp_(iGRU))
      endif

    endif  ! if there is a need to compute energy fluxes within vegetation

    ! -----
    ! * energy fluxes for the snow and soil domains...
    ! -------------------------------------------
    ! if(nLayers_(iGRU)>0)then
      do iLayer=1,nLayers_(iGRU) ! loop through all layers in the snow and soil domains

        ! check if the state is in the subset
        if(ixSnowSoilNrg_(iLayer,iGRU)==integerMissing) cycle
        ! - define index within the state subset and the full state vector
        nrgState = ixSnowSoilNrg_(iLayer,iGRU)        ! index within the state subset

        ! - diagonal elements
        aJac_(nrgState,nrgState,iGRU) = (dt/mLayerDepth_(iLayer,iGRU))*(-dNrgFlux_dTempBelow_(iLayer-1,iGRU) + dNrgFlux_dTempAbove_(iLayer,iGRU)) + aJac_(nrgState,nrgState,iGRU)

        ! - super-diagonal elements
        if(iLayer>1)then
          if(ixSnowSoilNrg_(iLayer-1,iGRU)/=integerMissing) aJac_(ixSnowSoilNrg_(iLayer-1,iGRU),nrgState,iGRU) = (dt/mLayerDepth_(iLayer-1,iGRU))*( dNrgFlux_dTempBelow_(iLayer-1,iGRU) )
        endif

        ! - sub-diagonal elements
        if(iLayer<nLayers_(iGRU))then
          if(ixSnowSoilNrg_(iLayer+1,iGRU)/=integerMissing) aJac_(ixSnowSoilNrg_(iLayer+1,iGRU),nrgState,iGRU) = (dt/mLayerDepth_(iLayer+1,iGRU))*(-dNrgFlux_dTempAbove_(iLayer,iGRU  ) )
        endif

      end do ! (looping through energy states in the snow and soil domains)
    ! endif ! (if the subset includes energy state variables in the snow and soil domains)

    ! -----
    ! * liquid water fluxes for the snow domain...
    ! --------------------------------------------
    ! if(nSnow_(iGRU)>0)then
      do iLayer=1,nSnow_(iGRU) ! loop through layers in the snow domain

        ! - check that the snow layer is desired
        if(ixSnowOnlyHyd_(iLayer,iGRU)==integerMissing) cycle
        ! - define state indices for the current layer
        watState = ixSnowOnlyHyd_(iLayer,iGRU)   ! hydrology state index within the state subset

        ! compute factor to convert liquid water derivative to total water derivative
        select case( ixHydType_(iLayer,iGRU) )
          case(iname_watLayer); convLiq2tot = mLayerFracLiqSnow_(iLayer,iGRU)
          case default;         convLiq2tot = 1._rkind
        end select

        ! - diagonal elements, water does not move upwards in snow
        aJac_(watState,watState,iGRU) = (dt/mLayerDepth_(iLayer,iGRU))*iLayerLiqFluxSnowDeriv_(iLayer,iGRU)*convLiq2tot + aJac_(watState,watState,iGRU)

        ! - sub-diagonal elements for snow, sub-diagonal only (water does not move upwards in snow)
        if(iLayer<nSnow_(iGRU) .and. mLayerVolFracIce_(iLayer+1,iGRU)<=maxVolIceContent)then
          aJac_(ixSnowOnlyHyd_(iLayer+1,iGRU),watState,iGRU) = -(dt/mLayerDepth_(iLayer+1,iGRU))*iLayerLiqFluxSnowDeriv_(iLayer,iGRU)*convLiq2tot  ! dVol(below)/dLiq(above)
        endif

      end do ! (looping through liquid water states in the snow domain)
    ! endif ! (if the subset includes hydrology state variables in the snow domain)

    ! -----
    ! * cross derivatives in the snow domain...
    ! ----------------------------------------
    ! if(nSnow_(iGRU)>0)then
      do iLayer=1,nSnow_(iGRU)  ! loop through layers in the snow domain

        ! (define the energy state)
        nrgState = ixSnowOnlyNrg_(iLayer,iGRU)       ! index within the full state vector
        ! - define state indices for the current layer
        watState = ixSnowOnlyHyd_(iLayer,iGRU)   ! hydrology state index within the state subset

        if(nrgState/=integerMissing .and. watState/=integerMissing)then
          ! - include derivatives of water fluxes w.r.t energy fluxes for current layer, water does not move upwards in snow
          aJac_(watState,nrgState,iGRU) = (dt/mLayerDepth_(iLayer,iGRU))*iLayerLiqFluxSnowDeriv_(iLayer,iGRU)*mLayerdTheta_dTk_(iLayer,iGRU)  ! (dVol/dT)             
        endif ! (if both the energy and water states for the current layer are within the state subset)

        if(watState/=integerMissing)then
          ! - include derivatives of heat capacity w.r.t water for layer above
          if(iLayer>1)then ! have layer above
            if(ixSnowSoilNrg_(iLayer-1,iGRU)/=integerMissing) aJac_(ixSnowSoilNrg_(iLayer-1,iGRU),watState,iGRU) = (dt/mLayerDepth_(iLayer-1,iGRU))*( dNrgFlux_dWatBelow_(iLayer-1,iGRU) )
          endif

          ! - include derivatives of heat capacity w.r.t water for layer below
          if(iLayer<nSnow_(iGRU) .or. (iLayer==nSnow_(iGRU) .and. nSoil>0))then ! have layer below
            if(ixSnowSoilNrg_(iLayer+1,iGRU)/=integerMissing) aJac_(ixSnowSoilNrg_(iLayer+1,iGRU),watState,iGRU) = (dt/mLayerDepth_(iLayer+1,iGRU))*(-dNrgFlux_dWatAbove_(iLayer,iGRU) )
          endif
        endif ! (if the water state for the current layer is within the state subset)

        if(nrgState/=integerMissing)then
          ! - sub-diagonal elements for snow, sub-diagonal only (water does not move upwards in snow)
          if(iLayer<nSnow_(iGRU) .and. mLayerVolFracIce_(iLayer+1,iGRU)<=maxVolIceContent)then
            aJac_(ixSnowOnlyHyd_(iLayer+1,iGRU),nrgState,iGRU) = -(dt/mLayerDepth_(iLayer+1,iGRU))*iLayerLiqFluxSnowDeriv_(iLayer,iGRU)*mLayerdTheta_dTk_(iLayer,iGRU)  ! dVol(below)/dT(above)
          endif 
        endif ! (if the energy state for the current layer is within the state subset)

      end do ! (looping through snow layers)
    ! endif ! (if there are state variables for both water and energy in the snow domain)

    ! -----
    ! * liquid water fluxes for the soil domain...
    ! --------------------------------------------
    if(nSoil>0)then
      do iLayer=1,nSoil

        ! - check that the soil layer is desired
        if(ixSoilOnlyHyd_(iLayer,iGRU)==integerMissing) cycle
        ! - define state indices
        watState = ixSoilOnlyHyd_(iLayer,iGRU)         ! hydrology state index within the state subset
        ! - define indices of the soil layers
        jLayer   = iLayer+nSnow_(iGRU)                  ! index of layer in the snow+soil vector
        ! - compute the diagonal elements
        ! all terms *excluding* baseflow
        aJac_(watState,watState,iGRU) = (dt/mLayerDepth_(jLayer,iGRU))*(-dq_dHydStateBelow_(iLayer-1,iGRU) + dq_dHydStateAbove_(iLayer,iGRU)) + aJac_(watState,watState,iGRU)

        ! - compute the super-diagonal elements
        if(iLayer>1)then
          if(ixSoilOnlyHyd_(iLayer-1,iGRU)/=integerMissing) aJac_(ixSoilOnlyHyd_(iLayer-1,iGRU),watState,iGRU) = (dt/mLayerDepth_(jLayer-1,iGRU))*( dq_dHydStateBelow_(iLayer-1,iGRU))
        endif

        ! - compute the sub-diagonal elements
        if(iLayer<nSoil)then
          if(ixSoilOnlyHyd_(iLayer+1,iGRU)/=integerMissing) aJac_(ixSoilOnlyHyd_(iLayer+1,iGRU),watState,iGRU) = (dt/mLayerDepth_(jLayer+1,iGRU))*(-dq_dHydStateAbove_(iLayer,iGRU))
        endif

        ! - include baseflow derivatives
        if(computeBaseflow)then
          do pLayer=1,nSoil
            qState = ixSoilOnlyHyd_(pLayer,iGRU)  ! hydrology state index within the state subset
            if(qState/=integerMissing)then
              aJac_(watState,qState,iGRU) = (dt/mLayerDepth_(jLayer,iGRU))*dBaseflow_dMatric_(iLayer,pLayer,iGRU) + aJac_(watState,qState,iGRU)
            endif
          end do
        endif ! (if computed baseflow)

        ! - include derivatives for surface infiltration below surface
        if(ixSoilOnlyHyd_(1,iGRU)/=integerMissing .and. all(dq_dHydStateLayerSurfVec_(:,iGRU)/=realMissing))then
          aJac_(ixSoilOnlyHyd_(1,iGRU),watState,iGRU) = -(dt/mLayerDepth_(nSnow_(iGRU)+1,iGRU))*dq_dHydStateLayerSurfVec_(iLayer,iGRU) + aJac_(ixSoilOnlyHyd_(1,iGRU),watState,iGRU)
        endif
      end do ! (looping through hydrology states in the soil domain)

      ! - include derivatives for surface infiltration above surface if there is snow (vegetation handled already)
      if(nSnow_(iGRU)>0 .and. ixSoilOnlyHyd_(1,iGRU)/=integerMissing .and. all(dq_dHydStateLayerSurfVec_(:,iGRU)/=realMissing))then ! have snow above first soil layer
        denseLimit = nSnow_(iGRU) ! if passed through a too dense snowpack, need to find top dense layer (bottom layer always included, dense or not)
        do pLayer=nSnow_(iGRU),1,-1
          if(mLayerVolFracIce_(pLayer,iGRU)<=maxVolIceContent) exit
          denseLimit = pLayer
        end do
        do pLayer=denseLimit,nSnow_(iGRU)
          if(ixSnowOnlyHyd_(pLayer,iGRU)/=integerMissing)then
            ! compute factor to convert liquid water derivative to total water derivative
            select case( ixHydType_(pLayer,iGRU) )
              case(iname_watLayer); convLiq2tot = mLayerFracLiqSnow_(pLayer,iGRU)
              case default;         convLiq2tot = 1._rkind
            end select
            aJac_(ixSoilOnlyHyd_(1,iGRU),ixSnowOnlyHyd_(pLayer,iGRU),iGRU) = -(dt/mLayerDepth_(nSnow_(iGRU)+1,iGRU))*scalarSoilControl_(iGRU)*iLayerLiqFluxSnowDeriv_(pLayer,iGRU)*convLiq2tot + aJac_(ixSoilOnlyHyd_(1,iGRU),ixSnowOnlyHyd_(pLayer,iGRU),iGRU)
          endif
        end do ! (looping through snow layers above soil until non-dense layer)
      endif ! (if snow present above soil)
    endif ! (if the subset includes hydrology state variables in the soil domain)

    ! -----
    ! * liquid water fluxes for the aquifer...
    ! ----------------------------------------
    if(ixAqWat_(iGRU)/=integerMissing) then
      aJac_(ixAqWat_(iGRU),ixAqWat_(iGRU),iGRU) = -dBaseflow_dAquifer_(iGRU)*dt + aJac_(ixAqWat_(iGRU),ixAqWat_(iGRU),iGRU)
      if(ixSoilOnlyNrg_(nSoil,iGRU)/=integerMissing) aJac_(ixAqWat_(iGRU),ixSoilOnlyNrg_(nSoil,iGRU),iGRU) = -dq_dNrgStateAbove_(nSoil,iGRU)*dt ! dAquiferRecharge_dTk  = d_iLayerLiqFluxSoil(nSoil)_dTk
      if(ixSoilOnlyHyd_(nSoil,iGRU)/=integerMissing) aJac_(ixAqWat_(iGRU),ixSoilOnlyHyd_(nSoil,iGRU),iGRU) = -dq_dHydStateAbove_(nSoil,iGRU)*dt ! dAquiferRecharge_dWat = d_iLayerLiqFluxSoil(nSoil)_dWat
      ! - include derivatives of energy and water w.r.t soil transpiration (dependent on canopy transpiration)
      if(computeVegFlux_(iGRU))then
        if(ixCasNrg_(iGRU)/=integerMissing)then
          aJac_(ixAqWat_(iGRU),ixCasNrg_(iGRU),iGRU) = -dAquiferTrans_dTCanair_(iGRU)*dt ! dVol/dT (K-1)
        endif
        if(ixVegNrg_(iGRU)/=integerMissing)then
          aJac_(ixAqWat_(iGRU),ixVegNrg_(iGRU),iGRU) = -dAquiferTrans_dTCanopy_(iGRU)*dt ! dVol/dT (K-1)
        endif
        if(ixVegHyd_(iGRU)/=integerMissing)then
          aJac_(ixAqWat_(iGRU),ixVegHyd_(iGRU),iGRU) = -dAquiferTrans_dCanWat_(iGRU)*dt  ! dVol/dLiq (kg m-2)-1
        endif
        if(ixTopNrg_(iGRU)/=integerMissing)then
          aJac_(ixAqWat_(iGRU),ixTopNrg_(iGRU),iGRU) = -dAquiferTrans_dTGround_(iGRU)*dt ! dVol/dT (K-1)
        endif
      endif
    endif ! (if aquifer water state is in the subset)

    ! -----
    ! * cross derivatives in the soil domain...
    ! ----------------------------------------
    if(nSoil>0)then
      do iLayer=1,nSoil

        ! - define indices of the soil layers
        jLayer   = iLayer+nSnow_(iGRU)                  ! index of layer in the snow+soil vector
        ! - define the energy state variable
        nrgState = ixSoilOnlyNrg_(iLayer,iGRU)         ! index within the full state vector
        ! - define index of hydrology state variable within the state subset
        watState = ixSoilOnlyHyd_(iLayer,iGRU)

        if(watState/=integerMissing .and. nrgState/=integerMissing)then
          ! - include derivatives in liquid water fluxes w.r.t. temperature for current layer
          aJac_(watState,nrgState,iGRU) = (dt/mLayerDepth_(jLayer,iGRU))*(-dq_dNrgStateBelow_(iLayer-1,iGRU) + dq_dNrgStateAbove_(iLayer,iGRU))   ! dVol/dT (K-1) -- flux depends on ice impedance

          ! - include derivatives w.r.t. ground evaporation
          if(nSnow_(iGRU)==0 .and. iLayer==1)then 
            aJac_(ixTopHyd_(iGRU),ixTopNrg_(iGRU),iGRU) = (dt/mLayerDepth_(jLayer,iGRU))*(-dGroundEvaporation_dTGround_(iGRU)/iden_water) + aJac_(ixTopHyd_(iGRU),ixTopNrg_(iGRU),iGRU) ! dVol/dT (K-1)
          endif
        endif   !(if both the energy and water states for the current layer are within the state subset)

        if(watState/=integerMissing)then
          ! - include derivatives w.r.t. ground evaporation
          if(nSnow_(iGRU)==0 .and. iLayer==1)then 
            if(computeVegFlux_(iGRU))then ! surface soil layer, assume here that kl>=4
              if(ixCasNrg_(iGRU)/=integerMissing) aJac_(ixTopHyd_(iGRU),ixCasNrg_(iGRU),iGRU) = (dt/mLayerDepth_(jLayer,iGRU))*(-dGroundEvaporation_dTCanair_(iGRU)/iden_water) ! dVol/dT (K-1)
              if(ixVegNrg_(iGRU)/=integerMissing) aJac_(ixTopHyd_(iGRU),ixVegNrg_(iGRU),iGRU) = (dt/mLayerDepth_(jLayer,iGRU))*(-dGroundEvaporation_dTCanopy_(iGRU)/iden_water) + aJac_(ixTopHyd_(iGRU),ixVegNrg_(iGRU),iGRU) ! dVol/dT (K-1)
              if(ixVegHyd_(iGRU)/=integerMissing) aJac_(ixTopHyd_(iGRU),ixVegHyd_(iGRU),iGRU) = (dt/mLayerDepth_(jLayer,iGRU))*(-dGroundEvaporation_dCanWat_(iGRU)/iden_water)  + aJac_(ixTopHyd_(iGRU),ixVegHyd_(iGRU),iGRU) ! dVol/dLiq (kg m-2)-1
            endif
          endif

          ! - include derivatives of energy and water w.r.t soil transpiration (dependent on canopy transpiration)
          if(computeVegFlux_(iGRU))then
            if(ixCasNrg_(iGRU)/=integerMissing)then
              aJac_(watState,ixCasNrg_(iGRU),iGRU) = (dt/mLayerDepth_(jLayer,iGRU))*(-mLayerdTrans_dTCanair_(iLayer,iGRU)) + aJac_(watState,ixCasNrg_(iGRU),iGRU) ! dVol/dT (K-1)
            endif
            if(ixVegNrg_(iGRU)/=integerMissing)then
              aJac_(watState,ixVegNrg_(iGRU),iGRU) = (dt/mLayerDepth_(jLayer,iGRU))*(-mLayerdTrans_dTCanopy_(iLayer,iGRU)) + aJac_(watState,ixVegNrg_(iGRU),iGRU) ! dVol/dT (K-1)
            endif
            if(ixVegHyd_(iGRU)/=integerMissing)then
              aJac_(watState,ixVegHyd_(iGRU),iGRU) = (dt/mLayerDepth_(jLayer,iGRU))*(-mLayerdTrans_dCanWat_(iLayer,iGRU))  + aJac_(watState,ixVegHyd_(iGRU),iGRU) ! dVol/dLiq (kg m-2)-1
            endif
            if(ixTopNrg_(iGRU)/=integerMissing)then
              aJac_(watState,ixTopNrg_(iGRU),iGRU) = (dt/mLayerDepth_(jLayer,iGRU))*(-mLayerdTrans_dTGround_(iLayer,iGRU)) + aJac_(watState,ixTopNrg_(iGRU),iGRU) ! dVol/dT (K-1)
            endif
          endif

          ! - include derivatives of heat capacity w.r.t water fluxes for layer above
          if(iLayer>1 .or. (iLayer==1 .and. nSnow_(iGRU)>0))then ! have layer above
            if(ixSnowSoilNrg_(jLayer-1,iGRU)/=integerMissing) aJac_(ixSnowSoilNrg_(jLayer-1,iGRU),watState,iGRU) = (dt/mLayerDepth_(jLayer-1,iGRU))*( dNrgFlux_dWatBelow_(jLayer-1,iGRU) )
          endif

          ! include derivatives of heat capacity w.r.t water fluxes for layer below
          if(iLayer<nSoil)then ! have layer below
            if(ixSnowSoilNrg_(jLayer+1,iGRU)/=integerMissing) aJac_(ixSnowSoilNrg_(jLayer+1,iGRU),watState,iGRU) = (dt/mLayerDepth_(jLayer+1,iGRU))*(-dNrgFlux_dWatAbove_(jLayer,iGRU) )
          endif
        endif ! (if the water state for the current layer is within the state subset)

        if(nrgState/=integerMissing)then
          ! - compute super-diagonal elements
          if(iLayer>1)then
            if(ixSoilOnlyHyd_(iLayer-1,iGRU)/=integerMissing) aJac_(ixSoilOnlyHyd_(iLayer-1,iGRU),nrgState,iGRU) = (dt/mLayerDepth_(jLayer-1,iGRU))*( dq_dNrgStateBelow_(iLayer-1,iGRU))   ! K-1
          endif

          ! compute sub-diagonal elements
          if(iLayer<nSoil)then
            if(ixSoilOnlyHyd_(iLayer+1,iGRU)/=integerMissing) aJac_(ixSoilOnlyHyd_(iLayer+1,iGRU),nrgState,iGRU) = (dt/mLayerDepth_(jLayer+1,iGRU))*(-dq_dNrgStateAbove_(iLayer,iGRU))     ! K-1
          endif

          ! - include derivatives for surface infiltration below surface
          if(ixSoilOnlyHyd_(1,iGRU)/=integerMissing .and. all(dq_dNrgStateLayerSurfVec_(:,iGRU)/=realMissing))then
            aJac_(ixSoilOnlyHyd_(1,iGRU),nrgState,iGRU) = -(dt/mLayerDepth_(nSnow_(iGRU)+1,iGRU))*dq_dNrgStateLayerSurfVec_(iLayer,iGRU) + aJac_(ixSoilOnlyHyd_(1,iGRU),nrgState,iGRU)
          endif
        endif ! (if the energy state for the current layer is within the state subset)
      end do ! (looping through energy states in the soil domain)

      ! - include derivatives for surface infiltration above surface if there is snow (vegetation handled already)
      if(nSnow_(iGRU)>0 .and. ixSoilOnlyHyd_(1,iGRU)/=integerMissing .and. all(dq_dNrgStateLayerSurfVec_(:,iGRU)/=realMissing))then ! have snow above first soil layer
        denseLimit = nSnow_(iGRU) ! if passed through a too dense snowpack, need to find top dense layer (bottom layer always included, dense or not)
        do pLayer=nSnow_(iGRU),1,-1
          if(mLayerVolFracIce_(pLayer,iGRU)<=maxVolIceContent) exit
          denseLimit = pLayer
        end do
        do pLayer=denseLimit,nSnow_(iGRU)
          if(ixSnowOnlyNrg_(pLayer,iGRU)/=integerMissing)then
            aJac_(ixSoilOnlyHyd_(1,iGRU),ixSnowOnlyNrg_(pLayer,iGRU),iGRU) = -(dt/mLayerDepth_(nSnow_(iGRU)+1,iGRU))*scalarSoilControl_(iGRU)*iLayerLiqFluxSnowDeriv_(pLayer,iGRU)*mLayerdTheta_dTk_(pLayer,iGRU) + aJac_(ixSoilOnlyHyd_(1,iGRU),ixSnowOnlyNrg_(pLayer,iGRU),iGRU)
          endif
        end do ! (looping through snow layers above soil until non-dense layer)
      endif ! (if snow present above soil)

    endif ! (if there are state variables for both water and energy in the soil domain)
   

end subroutine fluxJacAdd

! **********************************************************************************************************
! private function: get the index in the band-diagonal matrix or full matrix
! **********************************************************************************************************
function ixInd(jState,iState)
  implicit none
  integer(i4b),intent(in)  :: jState ! off-diagonal state
  integer(i4b),intent(in)  :: iState ! diagonal state
  integer(i4b)             :: ixInd  ! index in the band-diagonal matrix or full matrix

  if(fullMatrix) then
    ixInd = jState
  else
    ixInd = ixDiag + jState - iState
  endif
end function ixInd

end module computJacob_module
