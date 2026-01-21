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

module snowLiqFlx_module

! access modules
USE nrtype                                 ! numerical recipes data types
USE multiconst,only:iden_ice,iden_water    ! intrinsic density of ice and water (kg m-3)

! access missing values
USE globalData,only:integerMissing         ! missing integer
USE globalData,only:realMissing            ! missing real number
USE globalData,only:maxVolIceContent       ! snow maximum volumetric ice content to store water (-)

! named variables
USE var_lookup,only:iLookINDEX             ! named variables for structure elements
USE var_lookup,only:iLookPARAM             ! named variables for structure elements
USE var_lookup,only:iLookPROG              ! named variables for structure elements
USE var_lookup,only:iLookDIAG              ! named variables for structure elements
USE var_lookup,only:iLookDERIV,iLookFLUX

! data types
USE data_types,only:var_d                  ! x%var(:)     [rkind]
USE data_types,only:var_dlength            ! x%var(:)%dat [rkind]
USE data_types,only:var_ilength            ! x%var(:)%dat [i4b]
USE data_types,only:out_type_snowLiqFlx    ! data type for intent(out) arguments

! privacy
implicit none
private
public :: snowLiqFlx
contains
! ************************************************************************************************
! public subroutine snowLiqFlx: compute liquid water flux through the snowpack
! ************************************************************************************************
subroutine snowLiqFlx(&
                      ! input: model control, forcing, and model state vector
  nGRU, firstFluxCall,nSnow,&
                      mLayerVolFracLiqTrial,&
                      ! input-output: data structures
                      indx_data,               & ! intent(in):    model indices
                      mpar_data,               & ! intent(in):    model parameters
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                      ! input-output: fluxes and derivatives
                      flux_data, deriv_data, &
                      ! output: error control
                      out_snowLiqFlx)            ! intent(out):   error control
  use device_data_types
  use cudafor
  implicit none

  ! input: model control, forcing, and model state vector
  integer(i4b) :: nGRU
  logical(lgt) :: firstFluxCall
  integer(i4b),device                      :: nSnow(:)                      ! number of snow layers
  real(rkind),device :: mLayerVolFracLiqTrial(:,:)
  ! input-output: data structures
  type(indx_data_device),intent(in)      :: indx_data                  ! model indices
  type(mpar_data_device),intent(in)      :: mpar_data                  ! model parameters
  type(prog_data_device),intent(in)      :: prog_data                  ! prognostic variables for a local HRU
  type(diag_data_device),intent(inout)   :: diag_data                  ! diagnostic variables for a local HRU
  type(flux_data_device) :: flux_data
  type(deriv_data_device) :: deriv_data
  ! input-output: fluxes and derivatives
  ! output: error control
  type(out_type_snowLiqFlx)         :: out_snowLiqFlx             ! error control
  ! ------------------------------  ------------------------------------------------------------------------------------------------------------
  ! local variables
  integer(i4b)                      :: i                          ! search index for scalar solution
  integer(i4b)                      :: iLayer                     ! layer index
  integer(i4b)                      :: ixTop                      ! top layer in subroutine call
  integer(i4b)                      :: ixBot                      ! bottom layer in subroutine call
  real(rkind)                       :: multResid                  ! multiplier for the residual water content (-)
  real(rkind),parameter             :: residThrs=550._rkind       ! ice density threshold to reduce residual liquid water content (kg m-3)
  real(rkind),parameter             :: residScal=10._rkind        ! scaling factor for residual liquid water content reduction factor (kg m-3)
  real(rkind)                       :: availCap                   ! available storage capacity [0,1] (-)
  real(rkind)                       :: relSaturn                  ! relative saturation [0,1] (-)
  type(dim3) :: blocks,threads
  threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  ! ------------------------------------------------------------------------------------------------------------------------------------------
  ! make association of local variables with information in the data structures
  associate(&
    ! input: model control
    ! input: forcing for the snow domain
    scalarThroughfallRain   => flux_data%ixScalarThroughfallRain,   & ! intent(in): computed throughfall rate (kg m-2 s-1)
    scalarCanopyLiqDrainage => flux_data%ixScalarCanopyLiqDrainage, & ! intent(in): computed drainage of liquid water (kg m-2 s-1)
    ! input: model state vector
    ! input: snow properties and parameters
    mLayerVolFracIce => prog_data%mLayerVolFracIce, & ! intent(in):    volumetric ice content at the start of the time step (-)
    Fcapil           => mpar_data%Fcapil_,                & ! intent(in):    capillary retention as a fraction of the total pore volume (-)
    k_snow           => mpar_data%k_snow_,                & ! intent(in):    hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
    mw_exp           => mpar_data%mw_exp_,                & ! intent(in):    exponent for meltwater flow (-)
    ! input-output: diagnostic variables -- only computed for the first iteration
    mLayerPoreSpace  => diag_data%mLayerPoreSpace,           & ! intent(inout): pore space in each snow layer (-)
    mLayerThetaResid => diag_data%mLayerThetaResid,          & ! intent(inout): esidual volumetric liquid water content in each snow layer (-)
    ! input-output: fluxes and derivatives
    iLayerLiqFluxSnow_start      => flux_data%ixiLayerLiqFluxSnow_start,                & ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
    iLayerLiqFluxSnow_end => flux_data%ixiLayerLiqFluxSnow_end, &
    iLayerLiqFluxSnowDeriv => deriv_data%iLayerLiqFluxSnowDeriv_m,           & ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
    ! output: error control
    err                    => out_snowLiqFlx % err,                             & ! intent(out):   error code
    message                => out_snowLiqFlx % cmessage                         & ! intent(out):   error message
    ) ! end association of local variables with information in the data structures
    ! ------------------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='snowLiqFlx/'

    ! check that the input vectors match nSnow
    ! if (size(mLayerVolFracLiqTrial)/=nSnow .or. size(mLayerVolFracIce)/=nSnow .or. &
    !     size(iLayerLiqFluxSnow)/=nSnow+1 .or. size(iLayerLiqFluxSnowDeriv)/=nSnow+1) then
    !   err=20; message=trim(message)//'size mismatch of input/output vectors'; return
    ! end if

    ! check the meltwater exponent is >=1
    ! if (mw_exp<1._rkind) then; err=20; message=trim(message)//'meltwater exponent < 1'; return; end if

    ! get the indices for the snow+soil layers
    ! ixTop = 1
    ! ixBot = nSnow

    call snowLiqFlx_kernel<<<blocks,threads>>>(nGRU,firstFluxCall,&
  nSnow,&
  flux_data%data, &
  scalarThroughfallRain,scalarCanopyLiqDrainage,&
  mLayerVolFracLiqTrial, &
  mLayerVolFracIce,Fcapil,k_snow,mw_exp,&
  mLayerPoreSpace,mLayerThetaResid,&
  iLayerLiqFluxSnow_start,iLayerLiqFluxSnow_end,iLayerLiqFluxSnowDeriv)


  end associate ! end association of local variables with information in the data structures

end subroutine snowLiqFlx



attributes(global) subroutine snowLiqFlx_kernel(nGRU,firstFluxCall,&
  nSnow,&
  flux_data, &
  scalarThroughfallRain,scalarCanopyLiqDrainage,&
  mLayerVolFracLiqTrial, &
  mLayerVolFracIce,Fcapil,k_snow,mw_exp,&
  mLayerPoreSpace,mLayerThetaResid,&
  iLayerLiqFluxSnow_start,iLayerLiqFluxSnow_end,iLayerLiqFluxSnowDeriv)
  integer(i4b),value :: nGRU
  logical(lgt),value :: firstFluxCall
  integer(i4b) :: nSnow(:)
  real(rkind),intent(inout) :: flux_data(:,:)
      ! input: forcing for the snow domain
    integer(i4b),intent(in),value :: scalarThroughfallRain    ! intent(in): computed throughfall rate (kg m-2 s-1)
    integer(i4b),intent(in),value :: scalarCanopyLiqDrainage  ! intent(in): computed drainage of liquid water (kg m-2 s-1)
    ! input: model state vector
    real(rkind) :: mLayerVolFracLiqTrial(:,:) ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
    ! input: snow properties and parameters
    real(rkind) :: mLayerVolFracIce(:,:) ! intent(in):    volumetric ice content at the start of the time step (-)
    real(rkind) :: Fcapil(:)              ! intent(in):    capillary retention as a fraction of the total pore volume (-)
    real(rkind) :: k_snow(:)              ! intent(in):    hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
    real(rkind) :: mw_exp(:)              ! intent(in):    exponent for meltwater flow (-)
    ! input-output: diagnostic variables -- only computed for the first iteration
    real(rkind) :: mLayerPoreSpace(:,:)   ! intent(inout): pore space in each snow layer (-)
    real(rkind) :: mLayerThetaResid(:,:)  ! intent(inout): esidual volumetric liquid water content in each snow layer (-)
    ! input-output: fluxes and derivatives
    integer(i4b),intent(in),value :: iLayerLiqFluxSnow_start,iLayerLiqFluxSnow_end       ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
    real(rkind) :: iLayerLiqFluxSnowDeriv(0:,:)  ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)

    integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return

  if (nSnow(iGRU) .eq. 0) return

  call snowLiqFlx_device(firstFluxCall,nSnow(iGRU),&
  flux_data(scalarThroughfallRain,iGRU),flux_data(scalarCanopyLiqDrainage,iGRU),&
  mLayerVolFracLiqTrial(:,iGRU), &
  mLayerVolFracIce(:,iGRU),Fcapil(iGRU),k_snow(iGRU),mw_exp(iGRU),&
  mLayerPoreSpace(:,iGRU),mLayerThetaResid(:,iGRU),&
  flux_data(iLayerLiqFluxSnow_start:iLayerLiqFluxSnow_end,iGRU),iLayerLiqFluxSnowDeriv(:,iGRU))

end subroutine

attributes(device) subroutine snowLiqFlx_device(firstFluxCall,nSnow,&
  scalarThroughfallRain,scalarCanopyLiqDrainage,&
  mLayerVolFracLiqTrial, &
  mLayerVolFracIce,Fcapil,k_snow,mw_exp,&
  mLayerPoreSpace,mLayerThetaResid,&
  iLayerLiqFluxSnow,iLayerLiqFluxSnowDeriv)
  implicit none
    ! input: model control
    logical(lgt) :: firstFluxCall           ! intent(in): the first flux call
    ! logical(lgt) :: scalarSolution          ! intent(in): flag to denote if implementing the scalar solution
    integer(i4b) :: nSnow
    ! input: forcing for the snow domain
    real(rkind) :: scalarThroughfallRain    ! intent(in): computed throughfall rate (kg m-2 s-1)
    real(rkind) :: scalarCanopyLiqDrainage  ! intent(in): computed drainage of liquid water (kg m-2 s-1)
    ! input: model state vector
    real(rkind) :: mLayerVolFracLiqTrial(:) ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
    ! input: snow properties and parameters
    real(rkind) :: mLayerVolFracIce(:) ! intent(in):    volumetric ice content at the start of the time step (-)
    real(rkind) :: Fcapil              ! intent(in):    capillary retention as a fraction of the total pore volume (-)
    real(rkind) :: k_snow              ! intent(in):    hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
    real(rkind) :: mw_exp              ! intent(in):    exponent for meltwater flow (-)
    ! input-output: diagnostic variables -- only computed for the first iteration
    real(rkind) :: mLayerPoreSpace(:)   ! intent(inout): pore space in each snow layer (-)
    real(rkind) :: mLayerThetaResid(:)  ! intent(inout): esidual volumetric liquid water content in each snow layer (-)
    ! input-output: fluxes and derivatives
    real(rkind) :: iLayerLiqFluxSnow(0:)       ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
    real(rkind) :: iLayerLiqFluxSnowDeriv(0:)  ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)

    integer(i4b)                      :: iLayer                     ! layer index
    real(rkind)                       :: multResid                  ! multiplier for the residual water content (-)
    real(rkind),parameter             :: residThrs=550._rkind       ! ice density threshold to reduce residual liquid water content (kg m-3)
    real(rkind),parameter             :: residScal=10._rkind        ! scaling factor for residual liquid water content reduction factor (kg m-3)
    real(rkind)                       :: availCap                   ! available storage capacity [0,1] (-)
    real(rkind)                       :: relSaturn                  ! relative saturation [0,1] (-)

      ! define the liquid flux at the upper boundary (m s-1)
    iLayerLiqFluxSnow(0)      = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water
    iLayerLiqFluxSnowDeriv(0) = 0._rkind !computed inside computJacob

    ! compute properties fixed over the time step
    if (firstFluxCall) then
      ! loop through snow layers
      do iLayer=1,nSnow ! loop through snow layers
        multResid = 1._rkind/(1._rkind + exp((mLayerVolFracIce(iLayer)*iden_ice - residThrs)/residScal)) ! compute the reduction in liquid water holding capacity at high snow density (-)
        mLayerPoreSpace(iLayer)  = 1._rkind - mLayerVolFracIce(iLayer) ! compute the pore space (-)
        mLayerThetaResid(iLayer) = Fcapil*mLayerPoreSpace(iLayer)*multResid ! compute the residual volumetric liquid water content (-)
      end do  ! end looping through snow layers
    end if  ! end if the first flux call
     
    ! compute fluxes
    do iLayer=1,nSnow  ! loop through snow layers
      if (mLayerVolFracLiqTrial(iLayer) > mLayerThetaResid(iLayer)) then ! check that flow occurs
        ! compute the relative saturation (-)
        availCap  = mLayerPoreSpace(iLayer) - mLayerThetaResid(iLayer)                 ! available capacity
        relSaturn = (mLayerVolFracLiqTrial(iLayer) - mLayerThetaResid(iLayer)) / availCap    ! relative saturation
        iLayerLiqFluxSnow(iLayer)      = k_snow*relSaturn**mw_exp
        iLayerLiqFluxSnowDeriv(iLayer) = ( (k_snow*mw_exp)/availCap ) * relSaturn**(mw_exp - 1._rkind)
        if (mLayerVolFracIce(iLayer) > maxVolIceContent) then ! NOTE: use start-of-step ice content, to avoid convergence problems
          ! ** allow liquid water to pass through under very high ice density
          iLayerLiqFluxSnow(iLayer) = iLayerLiqFluxSnow(iLayer) + iLayerLiqFluxSnow(iLayer-1)
        end if
      else  ! flow does not occur
        iLayerLiqFluxSnow(iLayer)      = 0._rkind
        iLayerLiqFluxSnowDeriv(iLayer) = 0._rkind
      end if  ! storage above residual content
    end do  ! end loop through snow layers
  end subroutine

end module snowLiqFlx_module
