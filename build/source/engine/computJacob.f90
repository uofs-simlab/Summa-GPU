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
                      nSnow,                      & ! intent(in):    number of snow layers
                      nSoil,                      & ! intent(in):    number of soil layers
                      nLayers,                    & ! intent(in):    total number of layers
                      computeVegFlux,             & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      computeBaseflow,            & ! intent(in):    flag to indicate if we need to compute baseflow
                      ! input: data structures
                      ! indx_data,                  & ! intent(in):    index data
                      ixCasNrg, ixVegNrg, ixVegHyd, ixTopNrg, ixTopHyd, ixAqWat, &
                      ixSnowSoilNrg, ixSnowOnlyNrg, ixSoilOnlyNrg, &
                      ixSnowSoilHyd, ixSnowOnlyHyd, ixSoilOnlyHyd, &
                      ixHydType, &
                      ! prog_data,                  & ! intent(in):    model prognostic variables for a local HRU
                      mLayerDepth, mLayerVolFracIce, &
                      ! diag_data,                  & ! intent(in):    model diagnostic variables for a local HRU
                      mLayerFracLiqSnow, &
                      scalarFracLiqVeg, scalarSoilControl, canopyDepth, &
                      ! deriv_data,                 & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                      dCanairNetFlux_dCanairTemp,dCanairNetFlux_dCanopyTemp,dCanairNetFlux_dGroundTemp,dCanopyNetFlux_dCanairTemp,dCanopyNetFlux_dCanopyTemp,dCanopyNetFlux_dGroundTemp,dGroundNetFlux_dCanairTemp,dGroundNetFlux_dCanopyTemp,dGroundNetFlux_dCanWat, &
                      dCanopyEvaporation_dTCanair,dCanopyEvaporation_dTCanopy,dCanopyEvaporation_dTGround,dCanopyEvaporation_dCanWat,dGroundEvaporation_dTCanair,dGroundEvaporation_dTCanopy,dGroundEvaporation_dTGround,dGroundEvaporation_dCanWat, &
                      scalarCanopyLiqDeriv, dCanLiq_dTcanopy, &
                      dNrgFlux_dTempAbove, dNrgFlux_dTempBelow, &
                      dNrgFlux_dWatAbove, dNrgFlux_dWatBelow, &
                      mLayerdTrans_dTCanair,mLayerdTrans_dTCanopy,mLayerdTrans_dTGround,mLayerdTrans_dCanWat,&
                      dAquiferTrans_dTCanair,dAquiferTrans_dTCanopy,dAquiferTrans_dTGround,dAquiferTrans_dCanWat, &
                      iLayerLiqFluxSnowDeriv, &
                      dq_dHydStateAbove,dq_dHydStateBelow,dq_dHydStateLayerSurfVec,&
                      dBaseflow_dAquifer, &
                      dq_dNrgStateAbove,dq_dNrgStateBelow,dq_dNrgStateLayerSurfVec,&
                      mLayerdTheta_dTk, &
                      dBaseflow_dMatric,          & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
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
  integer(i4b),intent(in)              :: nSnow                      ! number of snow layers
  integer(i4b),intent(in)              :: nSoil                      ! number of soil layers
  integer(i4b),intent(in)              :: nLayers                    ! total number of layers in the snow and soil domains
  logical(lgt),intent(in)              :: computeVegFlux             ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)              :: computeBaseflow            ! flag to indicate if computing baseflow
  ! input: data structures
  ! type(var_ilength),intent(in)         :: indx_data                  ! indices defining model states and layers
  integer(i4b),intent(in) :: ixCasNrg, ixVegNrg, ixVegHyd, ixTopNrg, ixTopHyd, ixAqWat
  integer(i4b),intent(in) :: ixSnowSoilNrg(:), ixSnowOnlyNrg(:), ixSoilOnlyNrg(:)
  integer(i4b),intent(in) :: ixSnowSoilHyd(:), ixSnowOnlyHyd(:), ixSoilOnlyHyd(:)
  integer(i4b),intent(in) :: ixHydType(:)
  ! type(var_dlength),intent(in)         :: prog_data                  ! prognostic variables for a local HRU
  real(rkind),intent(in) :: mLayerDepth(:), mLayerVolFracIce(:)
  ! type(var_dlength),intent(in)         :: diag_data                  ! diagnostic variables for a local HRU
  real(rkind),intent(in) :: mLayerFracLiqSnow(:)
  real(rkind),intent(in) :: scalarFracLiqVeg, scalarSoilControl, canopyDepth
  ! type(var_dlength),intent(in)         :: deriv_data                 ! derivatives in model fluxes w.r.t. relevant state variables
  real(rkind),intent(in) :: dCanairNetFlux_dCanairTemp,dCanairNetFlux_dCanopyTemp,dCanairNetFlux_dGroundTemp,dCanopyNetFlux_dCanairTemp,dCanopyNetFlux_dCanopyTemp,dCanopyNetFlux_dGroundTemp,dGroundNetFlux_dCanairTemp,dGroundNetFlux_dCanopyTemp,dGroundNetFlux_dCanWat
  real(rkind),intent(in) :: dCanopyEvaporation_dTCanair,dCanopyEvaporation_dTCanopy,dCanopyEvaporation_dTGround,dCanopyEvaporation_dCanWat,dGroundEvaporation_dTCanair,dGroundEvaporation_dTCanopy,dGroundEvaporation_dTGround,dGroundEvaporation_dCanWat
  real(rkind),intent(in) :: scalarCanopyLiqDeriv, dCanLiq_dTcanopy
  real(rkind),intent(in) :: dNrgFlux_dTempAbove(0:), dNrgFlux_dTempBelow(0:)
  real(rkind),intent(in) :: dNrgFlux_dWatAbove(0:), dNrgFlux_dWatBelow(0:)
  real(rkind),intent(in) :: mLayerdTrans_dTCanair(:),mLayerdTrans_dTCanopy(:),mLayerdTrans_dTGround(:),mLayerdTrans_dCanWat(:)
  real(rkind),intent(in) :: dAquiferTrans_dTCanair,dAquiferTrans_dTCanopy,dAquiferTrans_dTGround,dAquiferTrans_dCanWat
  real(rkind),intent(in) :: iLayerLiqFluxSnowDeriv(0:)
  real(rkind),intent(in) :: dq_dHydStateAbove(0:),dq_dHydStateBelow(0:),dq_dHydStateLayerSurfVec(0:)
  real(rkind),intent(in) :: dBaseflow_dAquifer
  real(rkind),intent(in) :: dq_dNrgStateAbove(0:),dq_dNrgStateBelow(0:),dq_dNrgStateLayerSurfVec(0:)
  real(rkind),intent(in) :: mLayerdTheta_dTk(:)
  real(rkind),intent(in)               :: dBaseflow_dMatric(:,:)     ! derivative in baseflow w.r.t. matric head (s-1)
  ! input-output: Jacobian and its diagonal
  real(rkind),intent(in)               :: dMat(:)                    ! diagonal of the Jacobian matrix
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
    if(computeVegFlux)then ! (derivatives only defined when vegetation protrudes over the surface)

      ! * energy fluxes with the canopy water
      if(ixVegHyd/=integerMissing)then

        ! * cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
        if(ixCasNrg/=integerMissing) aJac_(ixVegHyd,ixCasNrg,iGRU) = -dCanopyEvaporation_dTCanair*dt
        ! dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy is the derivative in throughfall and canopy drainage with canopy temperature
        if(ixVegNrg/=integerMissing) aJac_(ixVegHyd,ixVegNrg,iGRU) = -dCanopyEvaporation_dTCanopy*dt + dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy
        ! * liquid water fluxes for vegetation canopy (-), dt*scalarFracLiqVeg*scalarCanopyLiqDeriv is the derivative in throughfall and canopy drainage with canopy water
                                     aJac_(ixVegHyd,ixVegHyd,iGRU) = -scalarFracLiqVeg*(dCanopyEvaporation_dCanWat - scalarCanopyLiqDeriv)*dt + dMat(ixVegHyd)
        if(ixTopNrg/=integerMissing) aJac_(ixVegHyd,ixTopNrg,iGRU) = -dCanopyEvaporation_dTGround*dt

        ! * cross-derivative terms w.r.t. canopy water (kg-1 m2)
        if(nSnow>0)then
          if(ixTopHyd/=integerMissing) aJac_(ixTopHyd,ixVegHyd,iGRU) = (dt/mLayerDepth(1))*(-scalarFracLiqVeg*scalarCanopyLiqDeriv)/iden_water
        else
          if(ixTopHyd/=integerMissing) aJac_(ixTopHyd,ixVegHyd,iGRU) = (dt/mLayerDepth(1))*(-scalarSoilControl*scalarFracLiqVeg*scalarCanopyLiqDeriv)/iden_water
        endif

        ! * cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
        if(ixTopNrg/=integerMissing) aJac_(ixTopNrg,ixVegHyd,iGRU) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanWat)
      endif

      ! * cross-derivative terms w.r.t. canopy temperature (K-1)
      if(ixVegNrg/=integerMissing)then
        if(nSnow>0)then
          if(ixTopHyd/=integerMissing) aJac_(ixTopHyd,ixVegNrg,iGRU) = (dt/mLayerDepth(1))*(-scalarCanopyLiqDeriv*dCanLiq_dTcanopy)/iden_water
        else
          if(ixTopHyd/=integerMissing) aJac_(ixTopHyd,ixVegNrg,iGRU) = (dt/mLayerDepth(1))*(-scalarSoilControl*scalarCanopyLiqDeriv*dCanLiq_dTcanopy)/iden_water
        endif
      endif

      ! * energy fluxes with the canopy air space (J m-3 K-1)
      if(ixCasNrg/=integerMissing)then
                                     aJac_(ixCasNrg,ixCasNrg,iGRU) = (dt/canopyDepth)*(-dCanairNetFlux_dCanairTemp) + dMat(ixCasNrg)
        if(ixVegNrg/=integerMissing) aJac_(ixCasNrg,ixVegNrg,iGRU) = (dt/canopyDepth)*(-dCanairNetFlux_dCanopyTemp)
        if(ixTopNrg/=integerMissing) aJac_(ixCasNrg,ixTopNrg,iGRU) = (dt/canopyDepth)*(-dCanairNetFlux_dGroundTemp)
      endif

      ! * energy fluxes with the vegetation canopy (J m-3 K-1)
      if(ixVegNrg/=integerMissing)then
        if(ixCasNrg/=integerMissing) aJac_(ixVegNrg,ixCasNrg,iGRU) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanairTemp)
                                     aJac_(ixVegNrg,ixVegNrg,iGRU) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanopyTemp) + dMat(ixVegNrg)
        if(ixTopNrg/=integerMissing) aJac_(ixVegNrg,ixTopNrg,iGRU) = (dt/canopyDepth)*(-dCanopyNetFlux_dGroundTemp)
      endif

      ! * energy fluxes with the surface (J m-3 K-1)
      if(ixTopNrg/=integerMissing)then
        if(ixCasNrg/=integerMissing) aJac_(ixTopNrg,ixCasNrg,iGRU) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanairTemp)
        if(ixVegNrg/=integerMissing) aJac_(ixTopNrg,ixVegNrg,iGRU) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanopyTemp)
      endif

    endif  ! if there is a need to compute energy fluxes within vegetation

    ! -----
    ! * energy fluxes for the snow and soil domains...
    ! -------------------------------------------
    if(nLayers>0)then
      do iLayer=1,nLayers ! loop through all layers in the snow and soil domains

        ! check if the state is in the subset
        if(ixSnowSoilNrg(iLayer)==integerMissing) cycle
        ! - define index within the state subset and the full state vector
        nrgState = ixSnowSoilNrg(iLayer)        ! index within the state subset

        ! - diagonal elements
        aJac_(nrgState,nrgState,iGRU) = (dt/mLayerDepth(iLayer))*(-dNrgFlux_dTempBelow(iLayer-1) + dNrgFlux_dTempAbove(iLayer)) + dMat(nrgState)

        ! - super-diagonal elements
        if(iLayer>1)then
          if(ixSnowSoilNrg(iLayer-1)/=integerMissing) aJac_(ixSnowSoilNrg(iLayer-1),nrgState,iGRU) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dTempBelow(iLayer-1) )
        endif

        ! - sub-diagonal elements
        if(iLayer<nLayers)then
          if(ixSnowSoilNrg(iLayer+1)/=integerMissing) aJac_(ixSnowSoilNrg(iLayer+1),nrgState,iGRU) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dTempAbove(iLayer  ) )
        endif

      end do ! (looping through energy states in the snow and soil domains)
    endif ! (if the subset includes energy state variables in the snow and soil domains)

    ! -----
    ! * liquid water fluxes for the snow domain...
    ! --------------------------------------------
    if(nSnow>0)then
      do iLayer=1,nSnow ! loop through layers in the snow domain

        ! - check that the snow layer is desired
        if(ixSnowOnlyHyd(iLayer)==integerMissing) cycle
        ! - define state indices for the current layer
        watState = ixSnowOnlyHyd(iLayer)   ! hydrology state index within the state subset

        ! compute factor to convert liquid water derivative to total water derivative
        select case( ixHydType(iLayer) )
          case(iname_watLayer); convLiq2tot = mLayerFracLiqSnow(iLayer)
          case default;         convLiq2tot = 1._rkind
        end select

        ! - diagonal elements, water does not move upwards in snow
        aJac_(watState,watState,iGRU) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*convLiq2tot + dMat(watState)

        ! - sub-diagonal elements for snow, sub-diagonal only (water does not move upwards in snow)
        if(iLayer<nSnow .and. mLayerVolFracIce(iLayer+1)<=maxVolIceContent)then
          aJac_(ixSnowOnlyHyd(iLayer+1),watState,iGRU) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*convLiq2tot  ! dVol(below)/dLiq(above)
        endif

      end do ! (looping through liquid water states in the snow domain)
    endif ! (if the subset includes hydrology state variables in the snow domain)

    ! -----
    ! * cross derivatives in the snow domain...
    ! ----------------------------------------
    if(nSnow>0)then
      do iLayer=1,nSnow  ! loop through layers in the snow domain

        ! (define the energy state)
        nrgState = ixSnowOnlyNrg(iLayer)       ! index within the full state vector
        ! - define state indices for the current layer
        watState = ixSnowOnlyHyd(iLayer)   ! hydrology state index within the state subset

        if(nrgState/=integerMissing .and. watState/=integerMissing)then
          ! - include derivatives of water fluxes w.r.t energy fluxes for current layer, water does not move upwards in snow
          aJac_(watState,nrgState,iGRU) = (dt/mLayerDepth(iLayer))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)  ! (dVol/dT)             
        endif ! (if both the energy and water states for the current layer are within the state subset)

        if(watState/=integerMissing)then
          ! - include derivatives of heat capacity w.r.t water for layer above
          if(iLayer>1)then ! have layer above
            if(ixSnowSoilNrg(iLayer-1)/=integerMissing) aJac_(ixSnowSoilNrg(iLayer-1),watState,iGRU) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dWatBelow(iLayer-1) )
          endif

          ! - include derivatives of heat capacity w.r.t water for layer below
          if(iLayer<nSnow .or. (iLayer==nSnow .and. nSoil>0))then ! have layer below
            if(ixSnowSoilNrg(iLayer+1)/=integerMissing) aJac_(ixSnowSoilNrg(iLayer+1),watState,iGRU) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dWatAbove(iLayer) )
          endif
        endif ! (if the water state for the current layer is within the state subset)

        if(nrgState/=integerMissing)then
          ! - sub-diagonal elements for snow, sub-diagonal only (water does not move upwards in snow)
          if(iLayer<nSnow .and. mLayerVolFracIce(iLayer+1)<=maxVolIceContent)then
            aJac_(ixSnowOnlyHyd(iLayer+1),nrgState,iGRU) = -(dt/mLayerDepth(iLayer+1))*iLayerLiqFluxSnowDeriv(iLayer)*mLayerdTheta_dTk(iLayer)  ! dVol(below)/dT(above)
          endif 
        endif ! (if the energy state for the current layer is within the state subset)

      end do ! (looping through snow layers)
    endif ! (if there are state variables for both water and energy in the snow domain)

    ! -----
    ! * liquid water fluxes for the soil domain...
    ! --------------------------------------------
    if(nSoil>0)then
      do iLayer=1,nSoil

        ! - check that the soil layer is desired
        if(ixSoilOnlyHyd(iLayer)==integerMissing) cycle
        ! - define state indices
        watState = ixSoilOnlyHyd(iLayer)         ! hydrology state index within the state subset
        ! - define indices of the soil layers
        jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector
        ! - compute the diagonal elements
        ! all terms *excluding* baseflow
        aJac_(watState,watState,iGRU) = (dt/mLayerDepth(jLayer))*(-dq_dHydStateBelow(iLayer-1) + dq_dHydStateAbove(iLayer)) + dMat(watState)

        ! - compute the super-diagonal elements
        if(iLayer>1)then
          if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac_(ixSoilOnlyHyd(iLayer-1),watState,iGRU) = (dt/mLayerDepth(jLayer-1))*( dq_dHydStateBelow(iLayer-1))
        endif

        ! - compute the sub-diagonal elements
        if(iLayer<nSoil)then
          if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac_(ixSoilOnlyHyd(iLayer+1),watState,iGRU) = (dt/mLayerDepth(jLayer+1))*(-dq_dHydStateAbove(iLayer))
        endif

        ! - include baseflow derivatives
        if(computeBaseflow)then
          do pLayer=1,nSoil
            qState = ixSoilOnlyHyd(pLayer)  ! hydrology state index within the state subset
            if(qState/=integerMissing)then
              aJac_(watState,qState,iGRU) = (dt/mLayerDepth(jLayer))*dBaseflow_dMatric(iLayer,pLayer) + aJac_(watState,qState,iGRU)
            endif
          end do
        endif ! (if computed baseflow)

        ! - include derivatives for surface infiltration below surface
        if(ixSoilOnlyHyd(1)/=integerMissing .and. all(dq_dHydStateLayerSurfVec/=realMissing))then
          aJac_(ixSoilOnlyHyd(1),watState,iGRU) = -(dt/mLayerDepth(nSnow+1))*dq_dHydStateLayerSurfVec(iLayer) + aJac_(ixSoilOnlyHyd(1),watState,iGRU)
        endif
      end do ! (looping through hydrology states in the soil domain)

      ! - include derivatives for surface infiltration above surface if there is snow (vegetation handled already)
      if(nSnow>0 .and. ixSoilOnlyHyd(1)/=integerMissing .and. all(dq_dHydStateLayerSurfVec/=realMissing))then ! have snow above first soil layer
        denseLimit = nSnow ! if passed through a too dense snowpack, need to find top dense layer (bottom layer always included, dense or not)
        do pLayer=nSnow,1,-1
          if(mLayerVolFracIce(pLayer)<=maxVolIceContent) exit
          denseLimit = pLayer
        end do
        do pLayer=denseLimit,nSnow
          if(ixSnowOnlyHyd(pLayer)/=integerMissing)then
            ! compute factor to convert liquid water derivative to total water derivative
            select case( ixHydType(pLayer) )
              case(iname_watLayer); convLiq2tot = mLayerFracLiqSnow(pLayer)
              case default;         convLiq2tot = 1._rkind
            end select
            aJac_(ixSoilOnlyHyd(1),ixSnowOnlyHyd(pLayer),iGRU) = -(dt/mLayerDepth(nSnow+1))*scalarSoilControl*iLayerLiqFluxSnowDeriv(pLayer)*convLiq2tot + aJac_(ixSoilOnlyHyd(1),ixSnowOnlyHyd(pLayer),iGRU)
          endif
        end do ! (looping through snow layers above soil until non-dense layer)
      endif ! (if snow present above soil)
    endif ! (if the subset includes hydrology state variables in the soil domain)

    ! -----
    ! * liquid water fluxes for the aquifer...
    ! ----------------------------------------
    if(ixAqWat/=integerMissing) then
      aJac_(ixAqWat,ixAqWat,iGRU) = -dBaseflow_dAquifer*dt + dMat(ixAqWat)
      if(ixSoilOnlyNrg(nSoil)/=integerMissing) aJac_(ixAqWat,ixSoilOnlyNrg(nSoil),iGRU) = -dq_dNrgStateAbove(nSoil)*dt ! dAquiferRecharge_dTk  = d_iLayerLiqFluxSoil(nSoil)_dTk
      if(ixSoilOnlyHyd(nSoil)/=integerMissing) aJac_(ixAqWat,ixSoilOnlyHyd(nSoil),iGRU) = -dq_dHydStateAbove(nSoil)*dt ! dAquiferRecharge_dWat = d_iLayerLiqFluxSoil(nSoil)_dWat
      ! - include derivatives of energy and water w.r.t soil transpiration (dependent on canopy transpiration)
      if(computeVegFlux)then
        if(ixCasNrg/=integerMissing)then
          aJac_(ixAqWat,ixCasNrg,iGRU) = -dAquiferTrans_dTCanair*dt ! dVol/dT (K-1)
        endif
        if(ixVegNrg/=integerMissing)then
          aJac_(ixAqWat,ixVegNrg,iGRU) = -dAquiferTrans_dTCanopy*dt ! dVol/dT (K-1)
        endif
        if(ixVegHyd/=integerMissing)then
          aJac_(ixAqWat,ixVegHyd,iGRU) = -dAquiferTrans_dCanWat*dt  ! dVol/dLiq (kg m-2)-1
        endif
        if(ixTopNrg/=integerMissing)then
          aJac_(ixAqWat,ixTopNrg,iGRU) = -dAquiferTrans_dTGround*dt ! dVol/dT (K-1)
        endif
      endif
    endif ! (if aquifer water state is in the subset)

    ! -----
    ! * cross derivatives in the soil domain...
    ! ----------------------------------------
    if(nSoil>0)then
      do iLayer=1,nSoil

        ! - define indices of the soil layers
        jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector
        ! - define the energy state variable
        nrgState = ixSoilOnlyNrg(iLayer)         ! index within the full state vector
        ! - define index of hydrology state variable within the state subset
        watState = ixSoilOnlyHyd(iLayer)

        if(watState/=integerMissing .and. nrgState/=integerMissing)then
          ! - include derivatives in liquid water fluxes w.r.t. temperature for current layer
          aJac_(watState,nrgState,iGRU) = (dt/mLayerDepth(jLayer))*(-dq_dNrgStateBelow(iLayer-1) + dq_dNrgStateAbove(iLayer))   ! dVol/dT (K-1) -- flux depends on ice impedance

          ! - include derivatives w.r.t. ground evaporation
          if(nSnow==0 .and. iLayer==1)then 
            aJac_(ixTopHyd,ixTopNrg,iGRU) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTGround/iden_water) + aJac_(ixTopHyd,ixTopNrg,iGRU) ! dVol/dT (K-1)
          endif
        endif   !(if both the energy and water states for the current layer are within the state subset)

        if(watState/=integerMissing)then
          ! - include derivatives w.r.t. ground evaporation
          if(nSnow==0 .and. iLayer==1)then 
            if(computeVegFlux)then ! surface soil layer, assume here that kl>=4
              if(ixCasNrg/=integerMissing) aJac_(ixTopHyd,ixCasNrg,iGRU) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTCanair/iden_water) ! dVol/dT (K-1)
              if(ixVegNrg/=integerMissing) aJac_(ixTopHyd,ixVegNrg,iGRU) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dTCanopy/iden_water) + aJac_(ixTopHyd,ixVegNrg,iGRU) ! dVol/dT (K-1)
              if(ixVegHyd/=integerMissing) aJac_(ixTopHyd,ixVegHyd,iGRU) = (dt/mLayerDepth(jLayer))*(-dGroundEvaporation_dCanWat/iden_water)  + aJac_(ixTopHyd,ixVegHyd,iGRU) ! dVol/dLiq (kg m-2)-1
            endif
          endif

          ! - include derivatives of energy and water w.r.t soil transpiration (dependent on canopy transpiration)
          if(computeVegFlux)then
            if(ixCasNrg/=integerMissing)then
              aJac_(watState,ixCasNrg,iGRU) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dTCanair(iLayer)) + aJac_(watState,ixCasNrg,iGRU) ! dVol/dT (K-1)
            endif
            if(ixVegNrg/=integerMissing)then
              aJac_(watState,ixVegNrg,iGRU) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dTCanopy(iLayer)) + aJac_(watState,ixVegNrg,iGRU) ! dVol/dT (K-1)
            endif
            if(ixVegHyd/=integerMissing)then
              aJac_(watState,ixVegHyd,iGRU) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dCanWat(iLayer))  + aJac_(watState,ixVegHyd,iGRU) ! dVol/dLiq (kg m-2)-1
            endif
            if(ixTopNrg/=integerMissing)then
              aJac_(watState,ixTopNrg,iGRU) = (dt/mLayerDepth(jLayer))*(-mLayerdTrans_dTGround(iLayer)) + aJac_(watState,ixTopNrg,iGRU) ! dVol/dT (K-1)
            endif
          endif

          ! - include derivatives of heat capacity w.r.t water fluxes for layer above
          if(iLayer>1 .or. (iLayer==1 .and. nSnow>0))then ! have layer above
            if(ixSnowSoilNrg(jLayer-1)/=integerMissing) aJac_(ixSnowSoilNrg(jLayer-1),watState,iGRU) = (dt/mLayerDepth(jLayer-1))*( dNrgFlux_dWatBelow(jLayer-1) )
          endif

          ! include derivatives of heat capacity w.r.t water fluxes for layer below
          if(iLayer<nSoil)then ! have layer below
            if(ixSnowSoilNrg(jLayer+1)/=integerMissing) aJac_(ixSnowSoilNrg(jLayer+1),watState,iGRU) = (dt/mLayerDepth(jLayer+1))*(-dNrgFlux_dWatAbove(jLayer) )
          endif
        endif ! (if the water state for the current layer is within the state subset)

        if(nrgState/=integerMissing)then
          ! - compute super-diagonal elements
          if(iLayer>1)then
            if(ixSoilOnlyHyd(iLayer-1)/=integerMissing) aJac_(ixSoilOnlyHyd(iLayer-1),nrgState,iGRU) = (dt/mLayerDepth(jLayer-1))*( dq_dNrgStateBelow(iLayer-1))   ! K-1
          endif

          ! compute sub-diagonal elements
          if(iLayer<nSoil)then
            if(ixSoilOnlyHyd(iLayer+1)/=integerMissing) aJac_(ixSoilOnlyHyd(iLayer+1),nrgState,iGRU) = (dt/mLayerDepth(jLayer+1))*(-dq_dNrgStateAbove(iLayer))     ! K-1
          endif

          ! - include derivatives for surface infiltration below surface
          if(ixSoilOnlyHyd(1)/=integerMissing .and. all(dq_dNrgStateLayerSurfVec/=realMissing))then
            aJac_(ixSoilOnlyHyd(1),nrgState,iGRU) = -(dt/mLayerDepth(nSnow+1))*dq_dNrgStateLayerSurfVec(iLayer) + aJac_(ixSoilOnlyHyd(1),nrgState,iGRU)
          endif
        endif ! (if the energy state for the current layer is within the state subset)
      end do ! (looping through energy states in the soil domain)

      ! - include derivatives for surface infiltration above surface if there is snow (vegetation handled already)
      if(nSnow>0 .and. ixSoilOnlyHyd(1)/=integerMissing .and. all(dq_dNrgStateLayerSurfVec/=realMissing))then ! have snow above first soil layer
        denseLimit = nSnow ! if passed through a too dense snowpack, need to find top dense layer (bottom layer always included, dense or not)
        do pLayer=nSnow,1,-1
          if(mLayerVolFracIce(pLayer)<=maxVolIceContent) exit
          denseLimit = pLayer
        end do
        do pLayer=denseLimit,nSnow
          if(ixSnowOnlyNrg(pLayer)/=integerMissing)then
            aJac_(ixSoilOnlyHyd(1),ixSnowOnlyNrg(pLayer),iGRU) = -(dt/mLayerDepth(nSnow+1))*scalarSoilControl*iLayerLiqFluxSnowDeriv(pLayer)*mLayerdTheta_dTk(pLayer) + aJac_(ixSoilOnlyHyd(1),ixSnowOnlyNrg(pLayer),iGRU)
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
