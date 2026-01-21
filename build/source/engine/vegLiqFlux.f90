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
  nGRU, &
                    computeVegFlux,scalarCanopyLiqTrial,&
                    decisions,&
                    ! input-output: data structures
                    indx_data, &
                    mpar_data,                    & ! intent(in): model parameters
                    diag_data,                    & ! intent(in): local HRU model diagnostic variables
                    flux_data,deriv_data,&
                    ! output
                    out_vegLiqFlux)                 ! intent(out): output rates, derivatives, and error control
  use device_data_types
  implicit none
  ! input
  integer(i4b) :: nGRU
  logical(lgt),device :: computeVegFlux(:)
  real(rkind),device :: scalarCanopyLiqTrial(:)
  type(decisions_device) :: decisions
  ! input-output: data structures
  type(mpar_data_device),intent(in)    :: mpar_data                      ! model parameters
  type(diag_data_device),intent(inout) :: diag_data                      ! model diagnostic variables for the local basin
  type(flux_data_device),intent(inout) :: flux_data
  type(deriv_data_device),intent(inout) :: deriv_data
  type(indx_data_device),intent(in) :: indx_data
  ! output
  type(out_type_vegLiqFlux)       :: out_vegLiqFlux                 ! output rates, derivatives, and error control
  type(dim3) :: blocks,threads
  threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  ! make association of local variables with information in the data structures
  associate(&
    scalarRainfall       => flux_data%ixScalarRainfall,       & ! intent(in): rainfall (kg m-2 s-1)
    ixCanopyInterception       => decisions%cIntercept, & ! intent(in): index defining choice of parameterization for canopy interception
    scalarCanopyLiqMax         => diag_data%scalarCanopyLiqMax,   & ! intent(in): maximum storage before canopy drainage begins (kg m-2 s-1)
    scalarThroughfallScaleRain => mpar_data%throughfallScaleRain_,& ! intent(in): fraction of rain that hits the ground without touching the canopy (-)
    scalarCanopyDrainageCoeff  => mpar_data%canopyDrainageCoeff_,  & ! intent(in): canopy drainage coefficient (s-1)
    scalarThroughfallRain        => flux_data%ixScalarThroughfallRain,        & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
    scalarCanopyLiqDrainage      => flux_data%ixScalarCanopyLiqDrainage,      & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
    scalarThroughfallRainDeriv   => deriv_data%scalarThroughfallRainDeriv,   & ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
    scalarCanopyLiqDrainageDeriv => deriv_data%scalarCanopyLiqDrainageDeriv, & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
    err                          => out_vegLiqFlux % err,                          & ! intent(out): error code
    message                      => out_vegLiqFlux % cmessage                      & ! intent(out): error message
    )
    ! ------------------------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="vegLiqFlux/"

    call vegLiqFlux_kernel<<<blocks,threads>>>(nGRU,flux_data%numFluxData,indx_data%ixVegHyd,computeVegFlux,ixCanopyInterception,&
  scalarCanopyLiqTrial,flux_data%data,scalarRainfall,&
  scalarCanopyLiqMax,scalarThroughfallScaleRain,scalarCanopyDrainageCoeff,&
  scalarThroughfallRain,scalarCanopyLiqDrainage,scalarThroughfallRainDeriv,scalarCanopyLiqDrainageDeriv)

  end associate ! end association of local variables with information in the data structures

end subroutine vegLiqFlux

attributes(global) subroutine vegLiqFlux_kernel(nGRU,numFluxData,ixVegHyd,computeVegFlux,ixCanopyInterception,&
  scalarCanopyLiqTrial,flux_data,scalarRainfall,&
  scalarCanopyLiqMax,scalarThroughfallScaleRain,scalarCanopyDrainageCoeff,&
  scalarThroughfallRain,scalarCanopyLiqDrainage,scalarThroughfallRainDeriv,scalarCanopyLiqDrainageDeriv)
  integer(i4b),value :: nGRU,numFluxData
  integer(i4b) :: ixVegHyd(nGRU)
      logical(lgt),intent(in) :: computeVegFlux(nGRU) ! intent(in): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
    real(rkind),intent(in) :: scalarCanopyLiqTrial(nGRU)  ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
    real(rkind),intent(in) :: flux_data(numFluxData,nGRU)
    integer(i4b),intent(in),value :: scalarRainfall        ! intent(in): rainfall (kg m-2 s-1)
    integer(i4b),intent(in) :: ixCanopyInterception        ! intent(in): index defining choice of parameterization for canopy interception
    real(rkind),intent(in) :: scalarCanopyLiqMax(nGRU)           ! intent(in): maximum storage before canopy drainage begins (kg m-2 s-1)
    real(rkind),intent(in) :: scalarThroughfallScaleRain(nGRU)   ! intent(in): fraction of rain that hits the ground without touching the canopy (-)
    real(rkind),intent(in) :: scalarCanopyDrainageCoeff(nGRU)    ! intent(in): canopy drainage coefficient (s-1)
    integer(i4b),intent(in),value :: scalarThroughfallRain         ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
    integer(i4b),intent(in),value :: scalarCanopyLiqDrainage       ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
    real(rkind),intent(inout) :: scalarThroughfallRainDeriv(nGRU)    ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
    real(rkind),intent(inout) :: scalarCanopyLiqDrainageDeriv(nGRU)  ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)

    integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return

  if (ixVegHyd(iGRU)==integerMissing) return
  call vegLiqFlux_device(computeVegFlux,ixCanopyInterception,&
  scalarCanopyLiqTrial,flux_data,scalarRainfall,&
  scalarCanopyLiqMax,scalarThroughfallScaleRain,scalarCanopyDrainageCoeff,&
  scalarThroughfallRain,scalarCanopyLiqDrainage,scalarThroughfallRainDeriv,scalarCanopyLiqDrainageDeriv)
end subroutine

attributes(device) subroutine vegLiqFlux_device(computeVegFlux_,ixCanopyInterception,&
  scalarCanopyLiqTrial_,flux_data,scalarRainfall_,&
  scalarCanopyLiqMax_,scalarThroughfallScaleRain_,scalarCanopyDrainageCoeff_,&
  scalarThroughfallRain_,scalarCanopyLiqDrainage_,scalarThroughfallRainDeriv_,scalarCanopyLiqDrainageDeriv_)
  use initialize_device,only:get_iGRU
    logical(lgt),intent(in) :: computeVegFlux_(:) ! intent(in): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
    real(rkind),intent(in) :: scalarCanopyLiqTrial_(:)  ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
    real(rkind) :: flux_data(:,:)
    integer(i4b),intent(in) :: scalarRainfall_        ! intent(in): rainfall (kg m-2 s-1)
    integer(i4b),intent(in) :: ixCanopyInterception        ! intent(in): index defining choice of parameterization for canopy interception
    real(rkind),intent(in) :: scalarCanopyLiqMax_(:)           ! intent(in): maximum storage before canopy drainage begins (kg m-2 s-1)
    real(rkind),intent(in) :: scalarThroughfallScaleRain_(:)   ! intent(in): fraction of rain that hits the ground without touching the canopy (-)
    real(rkind),intent(in) :: scalarCanopyDrainageCoeff_(:)    ! intent(in): canopy drainage coefficient (s-1)
    integer(i4b),intent(in) :: scalarThroughfallRain_         ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
    integer(i4b),intent(in) :: scalarCanopyLiqDrainage_       ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
    real(rkind),intent(inout) :: scalarThroughfallRainDeriv_(:)    ! intent(out): derivative in throughfall w.r.t. canopy liquid water (s-1)
    real(rkind),intent(inout) :: scalarCanopyLiqDrainageDeriv_(:)  ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
    integer(i4b) :: iGRU
    iGRU = get_iGRU()

    ! set throughfall to inputs if vegetation is completely buried with snow
    if (.not.computeVegFlux_(iGRU)) then
      flux_data(scalarThroughfallRain_,iGRU)        = flux_data(scalarRainfall_,iGRU)
      flux_data(scalarCanopyLiqDrainage_,iGRU)      = 0._rkind
      scalarThroughfallRainDeriv_(iGRU)   = 0._rkind
      scalarCanopyLiqDrainageDeriv_(iGRU) = 0._rkind
      return
    end if

    ! compute throughfall
    select case(ixCanopyInterception)
      ! original model (no flexibility in canopy interception): 100% of rainfall is intercepted by the vegetation canopy
      ! NOTE: this could be done with scalarThroughfallScaleRain=0, though requires setting scalarThroughfallScaleRain in all test cases
      case(unDefined)
        flux_data(scalarThroughfallRain_,iGRU)      = 0._rkind
        scalarThroughfallRainDeriv_(iGRU) = 0._rkind
      ! fraction of rainfall hits the ground without ever touching the canopy
      case(sparseCanopy)
        flux_data(scalarThroughfallRain_,iGRU)      = scalarThroughfallScaleRain_(iGRU)*flux_data(scalarRainfall_,iGRU)
        scalarThroughfallRainDeriv_(iGRU) = 0._rkind
      ! throughfall a function of canopy storage
      case(storageFunc)
        ! throughfall during wetting-up phase
        if(scalarCanopyLiqTrial_(iGRU) < scalarCanopyLiqMax_(iGRU))then
          flux_data(scalarThroughfallRain_,iGRU)      = flux_data(scalarRainfall_,iGRU)*(scalarCanopyLiqTrial_(iGRU)/scalarCanopyLiqMax_(iGRU))
          scalarThroughfallRainDeriv_(iGRU) = flux_data(scalarRainfall_,iGRU)/scalarCanopyLiqMax_(iGRU)
        ! all rain falls through the canopy when the canopy is at capacity
        else
          flux_data(scalarThroughfallRain_,iGRU)      = flux_data(scalarRainfall_,iGRU)
          scalarThroughfallRainDeriv_(iGRU) = 0._rkind
        end if
      ! case default; err=20; message=trim(message)//'unable to identify option for canopy interception'; return
    end select ! (option for canopy interception)

    ! compute canopy drainage
    if(scalarCanopyLiqTrial_(iGRU) > scalarCanopyLiqMax_(iGRU))then
      flux_data(scalarCanopyLiqDrainage_,iGRU)       = scalarCanopyDrainageCoeff_(iGRU)*(scalarCanopyLiqTrial_(iGRU) - scalarCanopyLiqMax_(iGRU))
      scalarCanopyLiqDrainageDeriv_(iGRU)  = scalarCanopyDrainageCoeff_(iGRU)
    else
      flux_data(scalarCanopyLiqDrainage_,iGRU)       = 0._rkind
      scalarCanopyLiqDrainageDeriv_(iGRU)  = 0._rkind
    end if


end subroutine

end module vegLiqFlux_module
