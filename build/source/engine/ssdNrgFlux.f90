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

module ssdNrgFlux_module

  ! data types
  USE nrtype
  
  ! data types
  USE data_types,only:var_d               ! x%var(:)     [rkind]
  USE data_types,only:var_dlength         ! x%var(:)%dat [rkind]
  USE data_types,only:var_ilength         ! x%var(:)%dat [i4b]
  USE data_types,only:in_type_ssdNrgFlux  ! intent(in) arguments for ssdNrgFlux
  USE data_types,only:io_type_ssdNrgFlux  ! intent(inout) arguments for ssdNrgFlux
  USE data_types,only:out_type_ssdNrgFlux ! intent(out) arguments for ssdNrgFlux
  
  ! physical constants
  USE multiconst,only:&
                      iden_water,  &  ! intrinsic density of water    (kg m-3)
                      Cp_water        ! specific heat of liquid water (J kg-1 K-1)
  
  ! missing values
  USE globalData,only:integerMissing  ! missing integer
  USE globalData,only:realMissing     ! missing real number
  
  ! named variables for snow and soil
  USE globalData,only:iname_snow      ! named variables for snow
  USE globalData,only:iname_soil      ! named variables for soil
  
  ! named variables
  USE var_lookup,only:iLookPROG       ! named variables for structure elements
  USE var_lookup,only:iLookDIAG       ! named variables for structure elements
  USE var_lookup,only:iLookFLUX       ! named variables for structure elements
  USE var_lookup,only:iLookPARAM      ! named variables for structure elements
  USE var_lookup,only:iLookINDEX      ! named variables for structure elements
  
  ! model decisions
  USE globalData,only:model_decisions ! model decision structure
  USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
  
  ! provide access to look-up values for model decisions
  USE mDecisions_module,only:      &
   ! look-up values for method used to compute derivative
   numerical,                      &  ! numerical solution
   analytical,                     &  ! analytical solution
   ! look-up values for choice of boundary conditions for thermodynamics
   prescribedTemp,                 &  ! prescribed temperature
   energyFlux,                     &  ! energy flux
   zeroFlux                           ! zero flux
  ! -------------------------------------------------------------------------------------------------
  implicit none
  private
  public :: ssdNrgFlux
  ! global parameters
  real(rkind),parameter            :: dx=1.e-10_rkind         ! finite difference increment (K)
  contains
  ! **********************************************************************************************************
  ! public subroutine ssdNrgFlux: compute energy fluxes and derivatives at layer interfaces
  ! **********************************************************************************************************
  subroutine ssdNrgFlux(&
                        ! input: model control, fluxes, trial variables, and  derivatives
                        ! in_ssdNrgFlux,                      & ! intent(in):     model control, fluxes, trial variables, and  derivatives
                        mLayerTempTrial, &
                        ! input-output: data structures and derivatives
                        mpar_data,                          & ! intent(in):    model parameters
                        indx_data,                          & ! intent(in):    model indices
                        nGRU, &
                        prog_data,                          & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,                          & ! intent(in):    model diagnostic variables for a local HRU
                        flux_data,                          & ! intent(inout): model fluxes for a local HRU
                        deriv_data, &
                        decisions, &
                        ! io_ssdNrgFlux,                      & ! intent(inout): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                        ! output: fluxes and derivatives at all layer interfaces and error control
                        err, message)                       ! intent(out):   derivatives and error control
    ! -------------------------------------------------------------------------------------------------------------------------------------------------
                        use device_data_types
                        use cudafor
    implicit none
    ! input: model control, fluxes, trial variables, and  derivatives
    ! type(in_type_ssdNrgFlux),intent(in)     :: in_ssdNrgFlux          ! input ssdNrgFlux arguments
    real(rkind),intent(in),device :: mLayerTempTrial(:,:)
    ! input-output: data structures
    type(mpar_data_device),intent(in)            :: mpar_data              ! model parameters
    type(indx_data_device),intent(in)            :: indx_data              ! state vector geometry
    integer(i4b) :: nGRU
    type(prog_data_device),intent(in)            :: prog_data              ! prognostic variables for a local HRU
    type(diag_data_device),intent(in)            :: diag_data              ! diagnostic variables for a local HRU
    type(flux_data_device),intent(inout)         :: flux_data              ! model fluxes for a local HRU
    type(deriv_data_device),intent(inout) :: deriv_data
    type(decisions_device),intent(in) :: decisions
    ! input-output: derivatives
    ! type(io_type_ssdNrgFlux),intent(inout)  :: io_ssdNrgFlux          ! input-output ssdNrgFlux arguments
    ! output: fluxes and derivatives at all layer interfaces
    integer(i4b) :: err
    character(*) :: message         ! output ssdNrgFlux arguments
    ! ------------------------------------------------------------------------------------------------------------------------------------------------------
    ! local variables
    !character(LEN=256)                  :: cmessage                   ! error message of downwind routine
    ! integer(i4b)                        :: nLayers                    ! number of model layers
    integer(i4b)                        :: iLayer                     ! index of model layers
    integer(i4b)                        :: ixLayerDesired(1)          ! layer desired (scalar solution)
    integer(i4b)                        :: ixTop                      ! top layer in subroutine call
    ! integer(i4b)                        :: ixBot                      ! bottom layer in subroutine call
    real(rkind)                         :: qFlux                      ! liquid flux at layer interfaces (m s-1)
    real(rkind)                         :: dz                         ! height difference (m)
    integer(i4b) :: iGRU
    ! ------------------------------------------------------------------------------------------------------------------------------------------------------
    ! allocate intent(out) data structure components
    ! nLayers=indx_data%nLayers

    ! allocate(&
    !   out_ssdNrgFlux % iLayerNrgFlux(0:nLayers),                          & ! energy flux at the layer interfaces (W m-2)
    !   out_ssdNrgFlux % dNrgFlux_dTempAbove(0:nLayers),                    & ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
    !   out_ssdNrgFlux % dNrgFlux_dTempBelow(0:nLayers),                    & ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
    !   out_ssdNrgFlux % dNrgFlux_dWatAbove(0:nLayers),                     & ! derivatives in the flux w.r.t. water state in the layer above (J m-2 s-1 K-1)
    !   out_ssdNrgFlux % dNrgFlux_dWatBelow(0:nLayers))                       ! derivatives in the flux w.r.t. water state in the layer below (J m-2 s-1 K-1)
    ! make association of local variables with information in the data structures
    associate(&
      nLayers => indx_data%nLayers_d, &
      ! input: model control
      ! scalarSolution             => in_ssdNrgFlux % scalarSolution,             & ! intent(in):    flag to denote if implementing the scalar solution
      ! input: fluxes and derivatives at the upper boundary
      groundNetFlux              => flux_data % scalarGroundNetNrgFlux,     & ! intent(in):    net energy flux for the ground surface (W m-2)
      dGroundNetFlux_dGroundTemp => deriv_data % dGroundNetFlux_dGroundTemp, & ! intent(inout): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
      ! input: liquid water fluxes
      iLayerLiqFluxSnow          => flux_data % iLayerLiqFluxSnow_m,          & ! intent(in):    liquid flux at the interface of each snow layer (m s-1)
      iLayerLiqFluxSoil          => flux_data % iLayerLiqFluxSoil_m,          & ! intent(in):    liquid flux at the interface of each soil layer (m s-1)
      ! input: trial model state variables
      ! mLayerTempTrial            => in_ssdNrgFlux % mLayerTempTrial,            & ! intent(in):     temperature in each layer at the current iteration (m)
      ! input: derivatives
      dThermalC_dWatAbove        => deriv_data % dThermalC_dWatAbove_m,  & ! intent(in): derivative in the thermal conductivity w.r.t. water state in the layer above
      dThermalC_dWatBelow        => deriv_data % dThermalC_dWatBelow_m,  & ! intent(in): derivative in the thermal conductivity w.r.t. water state in the layer above
      dThermalC_dTempAbove       => deriv_data % dThermalC_dTempAbove_m, & ! intent(in): derivative in the thermal conductivity w.r.t. energy state in the layer above
      dThermalC_dTempBelow       => deriv_data % dThermalC_dTempBelow_m, & ! intent(in): derivative in the thermal conductivity w.r.t. energy state in the layer above
      ! input: boundary conditions
      ix_bcUpprTdyn           => decisions%bcUpprTdyn, & ! intent(in):  method used to calculate the upper boundary condition for thermodynamics
      ix_bcLowrTdyn           => decisions%bcLowrTdyn, & ! intent(in):  method used to calculate the lower boundary condition for thermodynamics
      ! input: coordinate variables
      nSnow                   => indx_data%nSnow,               & ! intent(in):  number of snow layers
      layerType_m               => indx_data%layerType,              & ! intent(in):  layer type (iname_soil or iname_snow)
      ! ixLayerState            => indx_data%ixLayerState,           & ! intent(in):  list of indices for all model layers
      ixSnowSoilNrg_m           => indx_data%ixSnowSoilNrg,          & ! intent(in):  index in the state subset for energy state variables in the snow+soil domain
      mLayerDepth             => prog_data%mLayerDepth,             & ! intent(in):  depth of each layer (m)
      mLayerHeight            => prog_data%mLayerHeight,            & ! intent(in):  height at the mid-point of each layer (m)
      ! input: thermal properties
      upperBoundTemp          => mpar_data%upperBoundTemp,      & ! intent(in):  temperature of the upper boundary (K)
      lowerBoundTemp          => mpar_data%lowerBoundTemp,      & ! intent(in):  temperature of the lower boundary (K)
      iLayerThermalC          => diag_data%iLayerThermalC_m,          & ! intent(in):  thermal conductivity at the interface of each layer (W m-1 K-1)
       ! output: diagnostic fluxes
      iLayerConductiveFlux => flux_data%iLayerConductiveFlux_m,       & ! intent(out): conductive energy flux at layer interfaces at end of time step (W m-2)
      iLayerAdvectiveFlux  => flux_data%iLayerAdvectiveFlux_m,        & ! intent(out): advective energy flux at layer interfaces at end of time step (W m-2)
      ! output: fluxes and derivatives at all layer interfaces
      iLayerNrgFlux        => flux_data % iLayerNrgFlux_m,          & ! intent(out): energy flux at the layer interfaces (W m-2)
      dFlux_dTempAbove     => deriv_data % dNrgFlux_dTempAbove_m,    & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
      dFlux_dTempBelow_m     => deriv_data % dNrgFlux_dTempBelow_m,    & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
      dFlux_dWatAbove_m      => deriv_data % dNrgFlux_dWatAbove_m,     & ! intent(out): derivatives in the flux w.r.t. water state in the layer above (J m-2 s-1 K-1)
      dFlux_dWatBelow_m      => deriv_data % dNrgFlux_dWatBelow_m     & ! intent(out): derivatives in the flux w.r.t. water state in the layer below (J m-2 s-1 K-1)
      ! output: error control
      ! err                  => out_ssdNrgFlux % err,                    & ! intent(out): error code
      ! message              => out_ssdNrgFlux % cmessage                & ! intent(out): error message
      )  ! end association of local variables with information in the data structures
      ! ------------------------------------------------------------------------------------------------------------------------------------------------------
      ! initialize error control
      err=0; message='ssdNrgFlux/'
  
      ! set conductive and advective fluxes to missing in the upper boundary
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      iLayerConductiveFlux(0,iGRU) = realMissing
      iLayerAdvectiveFlux(0,iGRU)  = realMissing  !included in the ground heat flux
      end do
  
      ! get the indices for the snow+soil layers
      ! if (scalarSolution) then
      !   ixLayerDesired = pack(ixLayerState, ixSnowSoilNrg/=integerMissing)
      !   ixTop = ixLayerDesired(1)
      !   ixBot = ixLayerDesired(1)
      ! else
        ixTop = 1
        ! ixBot = nLayers
      ! end if
  
      ! -------------------------------------------------------------------------------------------------------------------------
      ! ***** compute the conductive fluxes at layer interfaces *****
      ! -------------------------------------------------------------------------------------------------------------------------
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      do iLayer=ixTop,nLayers(iGRU)
        if (iLayer==nLayers(iGRU)) then ! lower boundary fluxes -- positive downwards
          ! flux depends on the type of lower boundary condition
          select case(ix_bcLowrTdyn) ! identify the lower boundary condition for thermodynamics
            case(prescribedTemp); iLayerConductiveFlux(iLayer,iGRU) = -iLayerThermalC(iLayer,iGRU)*(lowerBoundTemp - mLayerTempTrial(iLayer,iGRU))/(mLayerDepth(iLayer,iGRU)*0.5_rkind)
            case(zeroFlux);       iLayerConductiveFlux(iLayer,iGRU) = 0._rkind
          end select  ! identifying the lower boundary condition for thermodynamics
        else ! domain boundary fluxes -- positive downwards
          iLayerConductiveFlux(iLayer,iGRU)  = -iLayerThermalC(iLayer,iGRU)*(mLayerTempTrial(iLayer+1,iGRU) - mLayerTempTrial(iLayer,iGRU)) / &
                                          (mLayerHeight(iLayer+1,iGRU) - mLayerHeight(iLayer,iGRU))
        end if ! the type of layer
      end do  ! end looping through layers
      end do
  
      ! -------------------------------------------------------------------------------------------------------------------------
      ! ***** compute the advective fluxes at layer interfaces *****
      ! -------------------------------------------------------------------------------------------------------------------------
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      do iLayer=ixTop,nLayers(iGRU)
        select case(layerType_m(iLayer,iGRU)) ! get the liquid flux at layer interfaces
          case(iname_snow); qFlux = iLayerLiqFluxSnow(iLayer,iGRU)
          case(iname_soil); qFlux = iLayerLiqFluxSoil(iLayer-nSnow(iGRU),iGRU)
          ! case default; err=20; 
            ! message=trim(message)//'unable to identify layer type'; return
        end select
        if (iLayer==nLayers(iGRU)) then ! compute fluxes at the lower boundary -- positive downwards
          iLayerAdvectiveFlux(iLayer,iGRU) = -Cp_water*iden_water*qFlux*(lowerBoundTemp - mLayerTempTrial(iLayer,iGRU))
        else ! compute fluxes within the domain -- positive downwards
          iLayerAdvectiveFlux(iLayer,iGRU) = -Cp_water*iden_water*qFlux*(mLayerTempTrial(iLayer+1,iGRU) - mLayerTempTrial(iLayer,iGRU))
        end if
      end do  ! end looping through layers
      end do
  
      ! -------------------------------------------------------------------------------------------------------------------------
      ! ***** compute the total fluxes at layer interfaces *****
      ! -------------------------------------------------------------------------------------------------------------------------
      ! NOTE: ignore advective fluxes for now
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        do iLayer=ixTop,nLayers(iGRU)
          iLayerNrgFlux(iLayer,iGRU) = iLayerConductiveFlux(iLayer,iGRU)
        end do
      end do
  
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      iLayerNrgFlux(0,iGRU)           = groundNetFlux(iGRU) ! from vegNrgFlux module
      ! iLayerNrgFlux(ixTop:ixBot,iGRU) = iLayerConductiveFlux(ixTop:ixBot,iGRU)
  
      ! -------------------------------------------------------------------------------------------------------------------
      ! ***** compute the derivative in fluxes at layer interfaces w.r.t state in the layer above and the layer below *****
      ! -------------------------------------------------------------------------------------------------------------------
  
      ! initialize un-used elements
      ! ***** the upper boundary
      dFlux_dTempAbove(0,iGRU) = 0._rkind ! this will be in canopy
      dFlux_dWatAbove_m(0,iGRU) = 0._rkind ! this will be in canopy
  
      ! ***** the lower boundary
      dFlux_dTempBelow_m(nLayers(iGRU),iGRU) = -huge(lowerBoundTemp)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems
      dFlux_dWatBelow_m(nLayers(iGRU),iGRU) = -huge(lowerBoundTemp)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems
  
      ! ***** the upper boundary, always do
      select case(ix_bcUpprTdyn)
  
        ! * prescribed temperature at the upper boundary
        case(prescribedTemp)
          dz = mLayerHeight(1,iGRU)*0.5_rkind
          dFlux_dWatBelow_m(0,iGRU)  = -dThermalC_dWatBelow(0,iGRU) * ( mLayerTempTrial(1,iGRU) - upperBoundTemp )/dz
          dFlux_dTempBelow_m(0,iGRU) = -dThermalC_dTempBelow(0,iGRU) * ( mLayerTempTrial(1,iGRU) - upperBoundTemp )/dz - iLayerThermalC(0,iGRU)/dz
  
        ! * zero flux at the upper boundary
        case(zeroFlux)
          dFlux_dWatBelow_m(0,iGRU) = 0._rkind
          dFlux_dTempBelow_m(0,iGRU) = 0._rkind
  
        ! * compute flux inside vegetation energy flux routine, use here
        case(energyFlux)
          dFlux_dWatBelow_m(0,iGRU) = 0._rkind !dGroundNetFlux_dGroundWat, does not exist in vegNrgFlux
          dFlux_dTempBelow_m(0,iGRU) = dGroundNetFlux_dGroundTemp(iGRU)
  
        ! case default; err=20; 
          ! message=trim(message)//'unable to identify upper boundary condition for thermodynamics'; return
  
      end select  ! end identifying the upper boundary condition for thermodynamics
      !dGroundNetFlux_dGroundWat  = dFlux_dWatBelow(0) ! this is true, but since not used in vegNrgFlux do not define
      ! print*, iGRU
      ! print*, dGroundNetFlux_dGroundTemp(iGRU)
      ! print*, dFlux_dTempBelow_m(0,iGRU)
      dGroundNetFlux_dGroundTemp(iGRU) = dFlux_dTempBelow_m(0,iGRU) ! need this in vegNrgFlux
      end do
  
      ! loop through INTERFACES...
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      do iLayer=ixTop,nLayers(iGRU)
          ! ***** the lower boundary
          if (iLayer==nLayers(iGRU)) then  ! if lower boundary
          ! identify the lower boundary condition
          select case(ix_bcLowrTdyn) ! prescribed temperature at the lower boundary
            case(prescribedTemp)
              dz = mLayerDepth(iLayer,iGRU)*0.5_rkind
              dFlux_dWatAbove_m(iLayer,iGRU)  = -dThermalC_dWatAbove(iLayer,iGRU) * ( lowerBoundTemp - mLayerTempTrial(iLayer,iGRU) )/dz
              dFlux_dTempAbove(iLayer,iGRU) = -dThermalC_dTempAbove(iLayer,iGRU) * ( lowerBoundTemp - mLayerTempTrial(iLayer,iGRU) )/dz + iLayerThermalC(iLayer,iGRU)/dz
            case(zeroFlux)  ! zero flux at the lower boundary
              dFlux_dWatAbove_m(iLayer,iGRU) = 0._rkind
              dFlux_dTempAbove(iLayer,iGRU) = 0._rkind
            ! case default; err=20; 
              ! message=trim(message)//'unable to identify lower boundary condition for thermodynamics'; return
            end select  ! end identifying the lower boundary condition for thermodynamics
        ! ***** internal layers
        else
          dz = (mLayerHeight(iLayer+1,iGRU) - mLayerHeight(iLayer,iGRU))
          dFlux_dWatAbove_m(iLayer,iGRU)  = -dThermalC_dWatAbove(iLayer,iGRU) * ( mLayerTempTrial(iLayer+1,iGRU) - mLayerTempTrial(iLayer,iGRU) )/dz
          dFlux_dWatBelow_m(iLayer,iGRU)  = -dThermalC_dWatBelow(iLayer,iGRU) * ( mLayerTempTrial(iLayer+1,iGRU) - mLayerTempTrial(iLayer,iGRU) )/dz
          dFlux_dTempAbove(iLayer,iGRU) = -dThermalC_dTempAbove(iLayer,iGRU) * ( mLayerTempTrial(iLayer+1,iGRU) - mLayerTempTrial(iLayer,iGRU) )/dz + iLayerThermalC(iLayer,iGRU)/dz
          dFlux_dTempBelow_m(iLayer,iGRU) = -dThermalC_dTempBelow(iLayer,iGRU) * ( mLayerTempTrial(iLayer+1,iGRU) - mLayerTempTrial(iLayer,iGRU) )/dz - iLayerThermalC(iLayer,iGRU)/dz
        end if  ! type of layer (upper, internal, or lower)
      end do  ! end looping through layers
      end do
  
  
    end associate ! end association of local variables with information in the data structures
  
  end subroutine ssdNrgFlux
  
  end module ssdNrgFlux_module
  
  