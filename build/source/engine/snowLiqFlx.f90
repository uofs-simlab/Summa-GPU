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
  
  ! named variables
  USE var_lookup,only:iLookINDEX             ! named variables for structure elements
  USE var_lookup,only:iLookPARAM             ! named variables for structure elements
  USE var_lookup,only:iLookPROG              ! named variables for structure elements
  USE var_lookup,only:iLookDIAG              ! named variables for structure elements
  
  ! data types
  USE data_types,only:var_d                  ! x%var(:)     [rkind]
  USE data_types,only:var_dlength            ! x%var(:)%dat [rkind]
  USE data_types,only:var_ilength            ! x%var(:)%dat [i4b]
  USE data_types,only:in_type_snowLiqFlx     ! data type for intent(in) arguments
  USE data_types,only:io_type_snowLiqFlx     ! data type for intent(inout) arguments
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
                        firstFluxCall,           & ! intent(in):    model control, forcing, and model state vector
                        mLayerVolFracLiqTrial, &
                        ! input-output: data structures
                        indx_data,               & ! intent(in):    model indices
                        nGRU, &
                        mpar_data,               & ! intent(in):    model parameters
                        prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                        flux_data, &
                        deriv_data, &
                        ! input-output: fluxes and derivatives
                        ! io_snowLiqFlx,           & ! intent(inout): fluxes and derivatives
                        ! output: error control
                        err,message)            ! intent(out):   error control
                        use device_data_types
    implicit none
    ! input: model control, forcing, and model state vector
    logical(lgt) :: firstFluxCall              ! model control, forcing, and model state vector
    real(rkind),intent(in),device :: mLayerVolFracLiqTrial(:,:)
    ! input-output: data structures
    type(indx_data_device),intent(in)      :: indx_data                  ! model indices
    integer(i4b) :: nGRU
    type(mpar_data_device),intent(in)      :: mpar_data                  ! model parameters
    type(prog_data_device),intent(in)      :: prog_data                  ! prognostic variables for a local HRU
    type(diag_data_device),intent(inout)   :: diag_data                  ! diagnostic variables for a local HRU
    type(flux_data_device),intent(inout) :: flux_data
    type(deriv_data_device),intent(inout) :: deriv_data
    ! input-output: fluxes and derivatives
    ! type(io_type_snowLiqFlx)          :: io_snowLiqFlx              ! fluxes and derivatives
    ! output: error control
    ! type(out_type_snowLiqFlx)         :: out_snowLiqFlx             ! error control
    integer(i4b) :: err
    character(*) :: message
    ! ------------------------------  ------------------------------------------------------------------------------------------------------------
    ! local variables
    ! integer(i4b)                      :: nSnow                      ! number of snow layers
    integer(i4b)                      :: i                          ! search index for scalar solution
    integer(i4b)                      :: iLayer                     ! layer index
    integer(i4b)                      :: ixTop                      ! top layer in subroutine call
    ! integer(i4b)                      :: ixBot                      ! bottom layer in subroutine call
    real(rkind)                       :: multResid                  ! multiplier for the residual water content (-)
    real(rkind),parameter             :: residThrs=550._rkind       ! ice density threshold to reduce residual liquid water content (kg m-3)
    real(rkind),parameter             :: residScal=10._rkind        ! scaling factor for residual liquid water content reduction factor (kg m-3)
    real(rkind),parameter             :: maxVolIceContent=0.7_rkind ! maximum volumetric ice content to store water (-)
    real(rkind)                       :: availCap                   ! available storage capacity [0,1] (-)
    real(rkind)                       :: relSaturn                  ! relative saturation [0,1] (-)
    integer(i4b) :: iGRU
    ! ------------------------------------------------------------------------------------------------------------------------------------------
    ! make association of local variables with information in the data structures
    ! nSnow=in_snowLiqFlx % nSnow ! get number of snow layers
    associate(&
      ! input: model control
      ! firstFluxCall           => in_snowLiqFlx % firstFluxCall,           & ! intent(in): the first flux call
      ! scalarSolution          => in_snowLiqFlx % scalarSolution,          & ! intent(in): flag to denote if implementing the scalar solution
      ! input: forcing for the snow domain
      scalarThroughfallRain   => flux_data % scalarThroughfallRain,   & ! intent(in): computed throughfall rate (kg m-2 s-1)
      scalarCanopyLiqDrainage => flux_data % scalarCanopyLiqDrainage, & ! intent(in): computed drainage of liquid water (kg m-2 s-1)
      ! input: model state vector
      ! mLayerVolFracLiqTrial   => in_snowLiqFlx % mLayerVolFracLiqTrial,   & ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
      ! input: layer indices
      ! ixLayerState     => indx_data%ixLayerState,             & ! intent(in):    list of indices for all model layers
      ixSnowOnlyHyd_m    => indx_data%ixSnowOnlyHyd,            & ! intent(in):    index in the state subset for hydrology state variables in the snow domain
      ! nSnow => indx_data%max_nSnow, &
      nSnow => indx_data%nSnow, &
      ! input: snow properties and parameters
      mLayerVolFracIce => prog_data%mLayerVolFracIce, & ! intent(in):    volumetric ice content at the start of the time step (-)
      Fcapil           => mpar_data%Fcapil,                & ! intent(in):    capillary retention as a fraction of the total pore volume (-)
      k_snow           => mpar_data%k_snow,                & ! intent(in):    hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
      mw_exp           => mpar_data%mw_exp,                & ! intent(in):    exponent for meltwater flow (-)
      ! input-output: diagnostic variables -- only computed for the first iteration
      mLayerPoreSpace  => diag_data%mLayerPoreSpace_m,           & ! intent(inout): pore space in each snow layer (-)
      mLayerThetaResid => diag_data%mLayerThetaResid_m,          & ! intent(inout): esidual volumetric liquid water content in each snow layer (-)
      ! input-output: fluxes and derivatives
      iLayerLiqFluxSnow      => flux_data % iLayerLiqFluxSnow_m,                & ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
      iLayerLiqFluxSnowDeriv => deriv_data % iLayerLiqFluxSnowDeriv_m           & ! intent(inout): derivative in vertical liquid water flux at layer interfaces (m s-1)
      ! output: error control
      ! err                    => out_snowLiqFlx % err,                             & ! intent(out):   error code
      ! message                => out_snowLiqFlx % cmessage                         & ! intent(out):   error message
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
      ixTop = integerMissing
      ! if (scalarSolution) then
      !   do i=1,size(ixSnowOnlyHyd)
      !     if (ixSnowOnlyHyd(i) /= integerMissing) then
      !       ixTop=ixLayerState(i)
      !       ixBot=ixTop
      !       exit  ! break out of loop once found
      !     end if
      !   end do
      !   if (ixTop == integerMissing) then
      !     err=20; message=trim(message)//'Unable to identify snow layer for scalar solution!'; return
      !   end if
      ! else
        ixTop = 1
        ! ixBot = nSnow
      ! end if
  
      ! define the liquid flux at the upper boundary (m s-1)
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      iLayerLiqFluxSnow(0,iGRU)      = (scalarThroughfallRain(iGRU) + scalarCanopyLiqDrainage(iGRU))/iden_water
      iLayerLiqFluxSnowDeriv(0,iGRU) = 0._rkind !computed inside computJacob
      enddo
  
      ! compute properties fixed over the time step
      if (firstFluxCall) then
        ! loop through snow layers
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
        do iLayer=1,nSnow(iGRU) ! loop through snow layers
          ! if (iLayer .gt. nSnow_d) cycle
          multResid = 1._rkind/(1._rkind + exp((mLayerVolFracIce(iLayer,iGRU)*iden_ice - residThrs)/residScal)) ! compute the reduction in liquid water holding capacity at high snow density (-)
          mLayerPoreSpace(iLayer,iGRU)  = 1._rkind - mLayerVolFracIce(iLayer,iGRU) ! compute the pore space (-)
          mLayerThetaResid(iLayer,iGRU) = Fcapil*mLayerPoreSpace(iLayer,iGRU)*multResid ! compute the residual volumetric liquid water content (-)
        end do  ! end looping through snow layers
      enddo
      end if  ! end if the first flux call
       
      ! compute fluxes
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      do iLayer=ixTop,nSnow(iGRU)  ! loop through snow layers
        ! if (iLayer .gt. nSnow_d) cycle
  
        if (mLayerVolFracLiqTrial(iLayer,iGRU) > mLayerThetaResid(iLayer,iGRU)) then ! check that flow occurs
          ! compute the relative saturation (-)
          availCap  = mLayerPoreSpace(iLayer,iGRU) - mLayerThetaResid(iLayer,iGRU)                 ! available capacity
          relSaturn = (mLayerVolFracLiqTrial(iLayer,iGRU) - mLayerThetaResid(iLayer,iGRU)) / availCap    ! relative saturation
          iLayerLiqFluxSnow(iLayer,iGRU)      = k_snow*relSaturn**mw_exp
          iLayerLiqFluxSnowDeriv(iLayer,iGRU) = ( (k_snow*mw_exp)/availCap ) * relSaturn**(mw_exp - 1._rkind)
          if (mLayerVolFracIce(iLayer,iGRU) > maxVolIceContent) then ! NOTE: use start-of-step ice content, to avoid convergence problems
            ! ** allow liquid water to pass through under very high ice density
            iLayerLiqFluxSnow(iLayer,iGRU) = iLayerLiqFluxSnow(iLayer,iGRU) + iLayerLiqFluxSnow(iLayer-1,iGRU) !NOTE: derivative may need to be updated in future.
          end if
        else  ! flow does not occur
          iLayerLiqFluxSnow(iLayer,iGRU)      = 0._rkind
          iLayerLiqFluxSnowDeriv(iLayer,iGRU) = 0._rkind
        end if  ! storage above residual content
      end do  ! end loop through snow layers
    end do
  
    end associate ! end association of local variables with information in the data structures
  
  end subroutine snowLiqFlx
  
  end module snowLiqFlx_module
  