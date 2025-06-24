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

module vegLiqFlux_module

  ! data types
  USE nrtype
  
  ! data types
  USE data_types,only:var_d                ! x%var(:)       (rkind)
  USE data_types,only:var_dlength          ! x%var(:)%dat   (rkind)
  USE data_types,only:in_type_vegLiqFlux   ! class type for intent(in) arguments
  USE data_types,only:out_type_vegLiqFlux  ! class type for intent(out) arguments
  
  ! named variables
  USE var_lookup,only:iLookPARAM,iLookDIAG ! named variables for structure elements
  
  ! model decisions
  USE globalData,only:model_decisions      ! model decision structure
  USE var_lookup,only:iLookDECISIONS       ! named variables for elements of the decision structure
  USE globalData,only:integerMissing  ! missing integer
  
  ! decisions on canopy interception parameterization
  USE mDecisions_module,only:         &
                        unDefined,    &    ! original model (no flexibility in canopy interception): 100% of rainfall is intercepted by the vegetation canopy
                        sparseCanopy, &    ! fraction of rainfall that never hits the canopy (throughfall); drainage above threshold
                        storageFunc        ! throughfall a function of canopy storage; 100% throughfall when canopy is at capacity
  
  ! privacy
  implicit none
  private
  public :: vegLiqFlux
  contains
  ! ************************************************************************************************
  ! public subroutine vegLiqFlux: compute water balance for the vegetation canopy
  ! ************************************************************************************************
  subroutine vegLiqFlux(&
                      ! input:
                      computeVegFlux,                & ! intent(in): model control, trial value, and rainfall rate
                      scalarCanopyLiqTrial, &
                      ! input-output: data structures
                      indx_data, &
                      mpar_data,                    & ! intent(in): model parameters
                      diag_data,                    & ! intent(in): local HRU model diagnostic variables
                      flux_data, &
                      deriv_data, &
                      decisions, &
                      nGRU, &
                      ! output
                      err,message)                 ! intent(out): output rates, derivatives, and error control
                      use device_data_types
    implicit none
    ! input
    logical(lgt)        :: computeVegFlux                  ! model control, trial value, and rainfall rate
    real(rkind),intent(in),device :: scalarCanopyLiqTrial(:)
    ! input-output: data structures
    type(indx_data_device) :: indx_data
    type(mpar_data_device),intent(in)    :: mpar_data                      ! model parameters
    type(diag_data_device),intent(inout) :: diag_data                      ! model diagnostic variables for the local basin
    type(flux_data_device),intent(inout) :: flux_data
    type(deriv_data_device),intent(inout) :: deriv_data
    type(decisions_device),intent(in) :: decisions
    integer(i4b) :: nGRU
    ! output
    integer(i4b) :: err
    character(*) :: message                 ! output rates, derivatives, and error control
    integer(i4b) :: iGRU
    ! ------------------------------------------------------------------------------------------------------------------------------------------------------
    ! make association of local variables with information in the data structures
    associate(&
      ! computeVegFlux       => in_vegLiqFlux % computeVegFlux,       & ! intent(in): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
      ixVegHyd => indx_data%ixVegHyd, &
      ! scalarCanopyLiqTrial => in_vegLiqFlux % scalarCanopyLiqTrial, & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
      scalarRainfall       => flux_data % scalarRainfall,       & ! intent(in): rainfall (kg m-2 s-1)
      ixCanopyInterception       => decisions%cIntercept, & ! intent(in): index defining choice of parameterization for canopy interception
      scalarCanopyLiqMax         => diag_data%scalarCanopyLiqMax,   & ! intent(in): maximum storage before canopy drainage begins (kg m-2 s-1)
      scalarThroughfallScaleRain => mpar_data%throughfallScaleRain,& ! intent(in): fraction of rain that hits the ground without touching the canopy (-)
      scalarCanopyDrainageCoeff  => mpar_data%canopyDrainageCoeff,  & ! intent(in): canopy drainage coefficient (s-1)
      scalarThroughfallRain        => flux_data % scalarThroughfallRain,        & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
      scalarCanopyLiqDrainage      => flux_data % scalarCanopyLiqDrainage,      & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
      scalarThroughfallRainDeriv   => deriv_data % scalarThroughfallRainDeriv,   & ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
      scalarCanopyLiqDrainageDeriv => deriv_data % scalarCanopyLiqDrainageDeriv & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
      ! err                          => out_vegLiqFlux % err,                          & ! intent(out): error code
      ! message                      => out_vegLiqFlux % cmessage                      & ! intent(out): error message
      )
      ! ------------------------------------------------------------------------------------------------------------------------------------------------------
      ! initialize error control
      err=0; message="vegLiqFlux/"
  
      ! set throughfall to inputs if vegetation is completely buried with snow
      if (.not.computeVegFlux) then
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
          if (ixVegHyd(iGRU)/=integerMissing) then
          scalarThroughfallRain(iGRU)        = scalarRainfall(iGRU)
          scalarCanopyLiqDrainage(iGRU)      = 0._rkind
          scalarThroughfallRainDeriv(iGRU)   = 0._rkind
          scalarCanopyLiqDrainageDeriv(iGRU) = 0._rkind
          end if
        end do
        return
      end if
  
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        if (ixVegHyd(iGRU)/=integerMissing) then
  
      ! compute throughfall
      select case(ixCanopyInterception)
        ! original model (no flexibility in canopy interception): 100% of rainfall is intercepted by the vegetation canopy
        ! NOTE: this could be done with scalarThroughfallScaleRain=0, though requires setting scalarThroughfallScaleRain in all test cases
        case(unDefined)
          scalarThroughfallRain(iGRU)      = 0._rkind
          scalarThroughfallRainDeriv(iGRU) = 0._rkind
        ! fraction of rainfall hits the ground without ever touching the canopy
        case(sparseCanopy)
          scalarThroughfallRain(iGRU)      = scalarThroughfallScaleRain*scalarRainfall(iGRU)
          scalarThroughfallRainDeriv(iGRU) = 0._rkind
        ! throughfall a function of canopy storage
        case(storageFunc)
          ! throughfall during wetting-up phase
          if(scalarCanopyLiqTrial(iGRU) < scalarCanopyLiqMax(iGRU))then
            scalarThroughfallRain(iGRU)      = scalarRainfall(iGRU)*(scalarCanopyLiqTrial(iGRU)/scalarCanopyLiqMax(iGRU))
            scalarThroughfallRainDeriv(iGRU) = scalarRainfall(iGRU)/scalarCanopyLiqMax(iGRU)
          ! all rain falls through the canopy when the canopy is at capacity
          else
            scalarThroughfallRain(iGRU)      = scalarRainfall(iGRU)
            scalarThroughfallRainDeriv(iGRU) = 0._rkind
          end if
        ! case default; err=20; message=trim(message)//'unable to identify option for canopy interception'; return
      end select ! (option for canopy interception)
  
      ! compute canopy drainage
      if(scalarCanopyLiqTrial(iGRU) > scalarCanopyLiqMax(iGRU))then
        scalarCanopyLiqDrainage(iGRU)       = scalarCanopyDrainageCoeff*(scalarCanopyLiqTrial(iGRU) - scalarCanopyLiqMax(iGRU))
        scalarCanopyLiqDrainageDeriv(iGRU)  = scalarCanopyDrainageCoeff
      else
        scalarCanopyLiqDrainage(iGRU)       = 0._rkind
        scalarCanopyLiqDrainageDeriv(iGRU)  = 0._rkind
      end if
    endif
      enddo
  
    end associate ! end association of local variables with information in the data structures
  
  end subroutine vegLiqFlux
  
  end module vegLiqFlux_module
  