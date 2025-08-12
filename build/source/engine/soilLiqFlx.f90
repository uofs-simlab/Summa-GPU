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

module soilLiqFlx_module
  ! -----------------------------------------------------------------------------------------------------------
  
  ! data types
  USE nrtype
  USE data_types,only:var_d                  ! x%var(:)       (rkind)
  USE data_types,only:var_ilength            ! x%var(:)%dat   (i4b)
  USE data_types,only:var_dlength            ! x%var(:)%dat   (rkind)
  USE data_types,only:in_type_soilLiqFlx     ! derived type for intent(in) arguments
  USE data_types,only:io_type_soilLiqFlx     ! derived type for intent(inout) arguments
  USE data_types,only:out_type_soilLiqFlx    ! derived type for intent(out) arguments
  USE data_types,only:in_type_diagv_node     ! derived type for intent(in) arguments 
  USE data_types,only:out_type_diagv_node    ! derived type for intent(out) arguments 
  USE data_types,only:in_type_surfaceFlx     ! derived type for intent(in) arguments
  USE data_types,only:io_type_surfaceFlx     ! derived type for intent(inout) arguments
  USE data_types,only:out_type_surfaceFlx    ! derived type for intent(out) arguments
  USE data_types,only:in_type_iLayerFlux     ! derived type for intent(in) arguments 
  USE data_types,only:out_type_iLayerFlux    ! derived type for intent(out) arguments 
  USE data_types,only:in_type_qDrainFlux     ! derived type for intent(in) arguments 
  USE data_types,only:out_type_qDrainFlux    ! derived type for intent(out) arguments 
  
  ! missing values
  USE globalData,only:integerMissing         ! missing integer
  USE globalData,only:realMissing            ! missing real number
  
  ! physical constants
  USE multiconst,only:iden_water             ! intrinsic density of water    (kg m-3)
  
  ! named variables
  USE var_lookup,only:iLookPROG              ! named variables for structure elements
  USE var_lookup,only:iLookDIAG              ! named variables for structure elements
  USE var_lookup,only:iLookFLUX              ! named variables for structure elements
  USE var_lookup,only:iLookPARAM             ! named variables for structure elements
  USE var_lookup,only:iLookINDEX,iLookDERIV             ! named variables for structure elements
  
  ! model decisions
  USE globalData,only:model_decisions        ! model decision structure
  USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
  
  ! provide access to look-up values for model decisions
  USE mDecisions_module,only:   &
    ! look-up values for method used to compute derivative
    numerical,                  & ! numerical solution
    analytical,                 & ! analytical solution
    ! look-up values for the form of Richards' equation
    moisture,                   & ! moisture-based form of Richards' equation
    mixdform,                   & ! mixed form of Richards' equation
    ! look-up values for the type of hydraulic conductivity profile
    constant,                   & ! constant hydraulic conductivity with depth
    powerLaw_profile,           & ! power-law profile
    ! look-up values for the choice of groundwater parameterization
    qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
    bigBucket,                  & ! a big bucket (lumped aquifer model)
    noExplicit,                 & ! no explicit groundwater parameterization
    ! look-up values for the choice of boundary conditions for hydrology
    prescribedHead,             & ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
    funcBottomHead,             & ! function of matric head in the lower-most layer
    freeDrainage,               & ! free drainage
    liquidFlux,                 & ! liquid water flux
    zeroFlux                      ! zero flux
  
  ! -----------------------------------------------------------------------------------------------------------
  implicit none
  private
  public :: soilLiqFlx
  ! constant parameters
  real(rkind),parameter     :: verySmall=1.e-12_rkind       ! a very small number (used to avoid divide by zero)
  real(rkind),parameter     :: dx=1.e-8_rkind               ! finite difference increment
  contains
  ! ***************************************************************************************************************
  ! public subroutine soilLiqFlx: compute liquid water fluxes and their derivatives
  ! ***************************************************************************************************************
  subroutine soilLiqFlx(&
                        ! input: model control, trial state variables, derivatives, and fluxes
                        firstSplitOper,                & ! intent(in): model control, trial state variables, derivatives, and fluxes
    nSoil, &
                        decisions, &
                        nGRU, &
                        above_soilLiqFluxDeriv,above_soildLiq_dTk,above_soilFracLiq, &
                        mLayerTempTrial, mLayerVolFracLiqTrial,mLayerVolFracIceTrial, &
                        mLayerMatricHeadLiqTrial,mLayerMatricHeadTrial, &
                        ! input-output: data structures
                        mpar_data,                    & ! intent(in):    model parameters
                        indx_data,                    & ! intent(in):    model indices
                        prog_data,                    & ! intent(in):    model prognostic variables for a local HRU
                        diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                        flux_data,                    & ! intent(inout): model fluxes for a local HRU
                        deriv_data, &
                        ! input-output: diagnostic variables, fluxes, and derivatives
                        ! io_soilLiqFlx,                & ! intent(inout): diagnostic variables, fluxes, and derivatives
                        dHydCond_dMatric, &
                        ! output: error control
                        err,message)                 ! intent(out): error control
    ! utility modules
    USE soil_utils_module,only:volFracLiq               ! compute volumetric fraction of liquid water
    USE soil_utils_module,only:matricHead               ! compute matric head (m)
    USE soil_utils_module,only:dTheta_dPsi              ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
    USE soil_utils_module,only:dPsi_dTheta              ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
    USE soil_utils_module,only:hydCond_psi              ! compute hydraulic conductivity as a function of matric head
    USE soil_utils_module,only:hydCond_liq              ! compute hydraulic conductivity as a function of volumetric liquid water content
    USE soil_utils_module,only:hydCondMP_liq            ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
    use device_data_types
    ! -------------------------------------------------------------------------------------------------------------------------------------------------
    implicit none
    ! input: model control, trial state variables, derivatives, and fluxes
    logical(lgt) :: firstSplitOper
    type(decisions_device),intent(in) :: decisions
    integer(i4b) :: nGRU
    real(rkind),intent(in),device :: above_soilLiqFluxDeriv(:),above_soildLiq_dTk(:),above_soilFracLiq(:)
    real(rkind),intent(in),device :: mLayerTempTrial(:,:), mLayerVolFracLiqTrial(:,:), mLayerVOlFracIceTrial(:,:)
    real(rkind),intent(in),device :: mLayerMatricHeadLiqTrial(:,:), mLayerMatricHeadTrial(:,:)
    ! input-output: data structures
    type(mpar_data_device),intent(in)           :: mpar_data                  ! model parameters
    type(indx_data_device),intent(in)           :: indx_data                  ! state vector geometry
    type(prog_data_device),intent(in)           :: prog_data                  ! prognostic variables for a local HRU
    type(diag_data_device),intent(inout)        :: diag_data                  ! diagnostic variables for a local HRU
    type(flux_data_device),intent(inout)        :: flux_data                  ! model fluxes for a local HRU
    type(deriv_data_device),intent(inout) :: deriv_data
    ! input-output: diagnostic variables, fluxes, and derivatives
    ! type(io_type_soilLiqFlx),intent(inout) :: io_soilLiqFlx              ! diagnostic variables, fluxes, and derivatives
    real(rkind),intent(inout),device :: dHydCond_dMatric(:,:)
    ! output: error control
    integer(i4b) :: err
    character(*) :: message
    ! -----------------------------------------------------------------------------------------------------------------------------------------------------
    ! local variables: general
    character(LEN=256)                  :: cmessage                      ! error message of downwind routine
    integer(i4b)                        :: nSoil                         ! number of soil layers
    ! integer(i4b)                        :: ibeg,iend                     ! start and end indices of the soil layers in concatanated snow-soil vector
    integer(i4b)                        :: iLayer,iSoil                  ! index of soil layer
    integer(i4b)                        :: ixLayerDesired(1)             ! layer desired (scalar solution)
    integer(i4b)                        :: ixTop                         ! top layer in subroutine call
    integer(i4b)                        :: ixBot                         ! bottom layer in subroutine call
    ! transpiration sink term
    real(rkind),dimension(nSoil,nGRU),device    :: mLayerTranspireFrac     ! fraction of transpiration allocated to each soil layer (-)
    ! diagnostic variables
    real(rkind),dimension(nSoil,nGRU),device    :: iceImpedeFac            ! ice impedence factor at layer mid-points (-)
    real(rkind),dimension(nSoil,nGRU),device    :: mLayerDiffuse           ! diffusivity at layer mid-point (m2 s-1)
    real(rkind),dimension(nSoil,nGRU),device    :: dHydCond_dVolLiq        ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
    real(rkind),dimension(nSoil,nGRU),device    :: dDiffuse_dVolLiq        ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
    real(rkind),dimension(nSoil,nGRU),device    :: dHydCond_dTemp          ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
    real(rkind),dimension(0:nSoil,nGRU),device  :: iLayerHydCond           ! hydraulic conductivity at layer interface (m s-1)
    real(rkind),dimension(0:nSoil,nGRU),device  :: iLayerDiffuse           ! diffusivity at layer interface (m2 s-1)
    ! compute surface flux
    integer(i4b),dimension(nGRU),device                                    :: nRoots                  ! number of soil layers with roots
    integer(i4b),dimension(nGRU),device                                    :: ixIce                   ! index of the lowest soil layer that contains ice
    real(rkind),dimension(0:nSoil,nGRU),device  :: iLayerHeight            ! height of the layer interfaces (m)
    ! error control
    logical(lgt)                                    :: return_flag             ! flag for return statements
    ! real(rkind),device :: mLayerDiffuse_d(size(mLayerDiffuse),1)
    ! real(rkind),device :: dHydCond_dTemp_d(size(dHydCond_dTemp),1)
    ! real(rkind),device :: dHydCond_dVolLiq_d(size(dHydCond_dVolLiq),1)
    ! real(rkind),device :: dDiffuse_dVolLiq_d(size(dDiffuse_dVolLiq),1)
    ! real(rkind),device :: dHydCond_dMatric_d(size(dHydCond_dMatric),1)
    ! real(rkind),device :: mLayerVolFracLiqTrial_d(size(mLayerVolFracLiqTrial),1)
    ! real(rkind),device :: mLayerMatricHeadLiqTrial_d(size(mLayerMatricHeadLiqTrial),1)
    ! real(rkind),device :: iLayerHydCond_d(0:in_soilLiqFlx%nSoil,1)
    ! real(rkind),device :: iLayerDiffuse_d(0:in_soilLiqFlx%nSoil,1)
    ! real(rkind),device :: iceImpedeFac_d(size(iceImpedeFac),1)
  
    ! real(rkind),device :: above_soildLiq_dTk_d(1)
    ! real(rkind),device :: above_soilLiqFluxDeriv_d(1)
    ! real(rkind),device :: above_soilFracLiq_d(1)
  
    ! real(rkind),device :: mLayerTempTrial_d(size(mLayerTempTrial),1)
    ! real(rkind),device :: mLayerVolFracIceTrial_d(size(mLayerVolFracIceTrial),1)
    ! real(rkind),device :: mLayerMatricHeadTrial_d(size(mLayerMatricHeadTrial),1)
    ! real(rkind),device :: iLayerHeight_d(0:in_soilLiqFlx%nSoil,1)
  
    ! integer(i4b),device :: nRoots_d(1)
    ! integer(i4b),device :: ixIce_d(1)
    
    ! real(rkind),device :: mLayerTranspireFrac_d(in_soilLiqFlx%nSoil,1)
  
    ! type(mpar_data_device) :: mpar_data_d
    ! type(flux_data_device) :: flux_data_d
    ! type(prog_data_device) :: prog_data_d
    ! type(decisions_device) :: decisions
    ! type(diag_data_device) :: diag_data_d
    ! type(deriv_data_device) :: deriv_data_d
    ! integer(i4b),device :: nSnow_d,nSoil_d
    ! -------------------------------------------------------------------------------------------------------------------------------------------------
  
  
    ! ** Initialize indices, error control, and get layer information ** 
    call initialize_soilLiqFlx; if (return_flag) return 
    ! ** Compute transpiration, diagnostic variables, infiltration, and interface fluxes **
    call update_soilLiqFlx(decisions,mpar_data,prog_data,diag_data,flux_data,deriv_data,indx_data%nSnow,nGRU,indx_data%nSoil);     if (return_flag) return
  
    ! ** Final error control **
    call finalize_soilLiqFlx;   if (return_flag) return
  
   
  contains
  
   subroutine initialize_soilLiqFlx
    ! **** Initial operations for soilLiqFlx module subroutine ****
  
    ! ** assign variables used in main associate block **
  
    ! get indices for the data structures
    ! ibeg = indx_data%nSnow + 1
    ! iend = indx_data%nSnow + indx_data%nSoil
  
    ! get a copy of iLayerHeight
    ! NOTE: performance hit, though cannot define the shape (0:) with the associate construct
    associate(iLayerH => prog_data%iLayerHeight,&
      nSnow_d => indx_data%nSnow)
    !$cuf kernel do(2) <<<*,*>>>
    do iGRU=1,nGRU
      do iLayer=0,nSoil
        iLayerHeight(iLayer,iGRU) = iLayerH(iLayer+nSnow_d(iGRU),iGRU)
      end do
    end do
    end associate
  
    ! ** initialize error control **
    return_flag=.false.
    ! associate(&
    !   err                   => out_soilLiqFlx % err,                  & ! intent(out): error code
    !   message               => out_soilLiqFlx % cmessage              & ! intent(out): error message
    ! &)
     err=0; message='soilLiqFlx/' ! initialize error control
    ! end associate
  
    ! ** get the indices for the soil layers **
    ! associate(&
    !  scalarSolution => in_soilLiqFlx % scalarSolution,             & ! intent(in): flag to denote if implementing the scalar solution
    !  ixMatricHead   => indx_data%var(iLookINDEX%ixMatricHead)%dat, & ! intent(in): indices of soil layers where matric head is the state variable
    !  ixSoilOnlyHyd  => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat & ! intent(in): index in the state subset for hydrology state variables in the soil domain
    ! &)
    !  if (scalarSolution) then
    !    ixLayerDesired = pack(ixMatricHead, ixSoilOnlyHyd/=integerMissing)
    !    ixTop = ixLayerDesired(1)
    !    ixBot = ixLayerDesired(1)
    !  else
       ixTop = 1
       ixBot = nSoil
    !  end if
    ! end associate
  
    ! ** identify the number of layers that contain roots **
    associate(&
     rootingDepth => mpar_data%rootingDepth& ! intent(in): rooting depth (m)
    !  err          => out_soilLiqFlx % err,                         & ! intent(out): error code
    !  message      => out_soilLiqFlx % cmessage                     & ! intent(out): error message
    &) 
    nRoots = 0
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
      do iLayer=0,nSoil-1
        if (iLayerHeight(iLayer,iGRU) < rootingDepth-verySmall) nRoots(iGRU) = nRoots(iGRU) + 1
      end do
    end do
  
    ! nRoots = 0
    !  nRoots = count(iLayerHeight(0:nSoil-1) < rootingDepth-verySmall)
    !  if (nRoots==0) then
    !    message=trim(message)//'no layers with roots'
    !    err=20; return_flag=.true.; return
    !  end if
    end associate
  
    ! ** identify lowest soil layer with ice **
    ! NOTE: cannot use count because there may be an unfrozen wedge
    associate(&
      nSnow => indx_data%nSnow &
    !   mLayerVolFracIceTrial => in_soilLiqFlx % mLayerVolFracIceTrial & ! intent(in): volumetric fraction of ice at the current iteration (-)
    &)
     ixIce = 0  ! initialize the index of the ice layer (0 means no ice in the soil profile)
     !$cuf kernel do(1) <<<*,*>>>
     do iGRU=1,nGRU
     do iLayer=1,nSoil ! (loop through soil layers)
       if (mLayerVolFracIceTrial(iLayer+nSnow(iGRU),iGRU) > verySmall) ixIce(iGRU) = iLayer
     end do
    end do
  
    end associate
   end subroutine initialize_soilLiqFlx
  
   subroutine update_soilLiqFlx(decisions,mpar_data,prog_data,diag_data,flux_data,deriv_data,nSnow_d,nGRU,nSoil_d)
    ! **** Main computations for soilLiqFlx module subroutine ****
    use device_data_types
    type(prog_data_device),intent(in) :: prog_data
    type(decisions_device),intent(inout) :: decisions
    type(flux_data_device),intent(inout) :: flux_data
    type(deriv_data_device),intent(inout) :: deriv_data
    type(mpar_data_device),intent(in) :: mpar_data
    type(diag_data_device),intent(inout) :: diag_data
    integer(i4b),device,intent(inout) :: nSnow_d(:)
    integer(i4b),intent(inout) :: nGRU
    integer(i4b),intent(in) :: nSoil_d
  
    ! if ( .not. (in_soilLiqFlx % scalarSolution .and. ixTop>1) ) then ! check the need to compute transpiration
     call compute_transpiration_sink(nGRU,decisions,diag_data,flux_data,deriv_data); if (return_flag) return
    ! end if  
  
    call compute_diagnostic_variables(nSnow_d,nGRU,decisions,mpar_data,flux_data,deriv_data,diag_data); if (return_flag) return
  
    call compute_surface_infiltration(nSnow_d,nSoil_d,nGRU,decisions,mpar_data,prog_data,flux_data,deriv_data,diag_data); if (return_flag) return
  
    call compute_interface_fluxes_derivatives(nGRU,nSnow_d,prog_data,decisions,flux_data,deriv_data); if (return_flag) return
  
    ! if ( .not. (in_soilLiqFlx % scalarSolution .and. ixTop<nSoil) ) then ! define the need to compute drainage
     call compute_drainage_flux(nGRU,nSoil_d,nSnow_d,decisions,mpar_data,prog_data,diag_data,flux_data,deriv_data); if (return_flag) return
    ! end if
  
   end subroutine update_soilLiqFlx
  
   subroutine finalize_soilLiqFlx
    ! **** Final operations for soilLiqFlx module subroutine ****
   
    ! final error control check for robustness
    ! associate(&
    !  err          => out_soilLiqFlx % err,                         & ! intent(out): error code
    !  message      => out_soilLiqFlx % cmessage                     & ! intent(out): error message
    ! &)
     if (err/=0) then; message=trim(message)//trim("finalize_soilLiqFlx: final error check failed"); return_flag=.true.; return; end if
    ! end associate
   end subroutine finalize_soilLiqFlx
  
   subroutine compute_transpiration_sink(nGRU,decisions,diag_data_d,flux_data_d,deriv_data_d)
    integer(i4b) :: nGRU
    type(decisions_device) :: decisions
    type(diag_data_device) :: diag_data_d
    type(flux_data_device) :: flux_data_d
    type(deriv_data_device) :: deriv_data_d
    real(rkind) :: sum_rootDens
    integer(i4b) :: iLayer2
    ! **** Compute the transpiration sink term ****
    ! **** Update the fraction of transpiration loss from each soil layer *****
    associate(&
      scalarTranspireLim => diag_data_d%scalarTranspireLim, & ! intent(in): weighted average of the transpiration limiting factor (-)
      mLayerRootDensity  => diag_data_d%mLayerRootDensity_m,     & ! intent(in): root density in each layer (-)
      mLayerTranspireLim => diag_data_d%mLayerTranspireLim_m,     & ! intent(in): transpiration limiting factor in each layer (-)
      scalarCanopyTranspiration => flux_data_d%scalarCanopyTranspiration, & ! canopy transpiration (kg m-2 s-1)
      mLayerTranspire           => flux_data_d%mLayerTranspire_m,   & ! transpiration loss from each soil layer (m s-1)
      ! intent(inout): derivatives in the soil layer transpiration flux ...
      mLayerdTrans_dCanWat  => deriv_data_d%mLayerdTrans_dCanWat_m,  & ! ... w.r.t. canopy total water
      mLayerdTrans_dTCanair => deriv_data_d%mLayerdTrans_dTCanair_m, & ! ... w.r.t. canopy air temperature
      mLayerdTrans_dTCanopy => deriv_data_d%mLayerdTrans_dTCanopy_m, & ! ... w.r.t. canopy temperature
      mLayerdTrans_dTGround => deriv_data_d%mLayerdTrans_dTGround_m, & ! ... w.r.t. ground temperature
      ! intent(in): derivative in canopy transpiration ...
      dCanopyTrans_dCanWat  => deriv_data_d%dCanopyTrans_dCanWat,  & ! ... w.r.t. canopy total water content (s-1)
      dCanopyTrans_dTCanair => deriv_data_d%dCanopyTrans_dTCanair, & ! ... w.r.t. canopy air temperature (kg m-2 s-1 K-1)
      dCanopyTrans_dTCanopy => deriv_data_d%dCanopyTrans_dTCanopy, & ! ... w.r.t. canopy temperature (kg m-2 s-1 K-1)
      dCanopyTrans_dTGround => deriv_data_d%dCanopyTrans_dTGround, & ! ... w.r.t. ground temperature (kg m-2 s-1 K-1)
      ! intent(in): index of the upper boundary conditions for soil hydrology
      ixBcUpperSoilHydrology => decisions%bcUpprSoiH & 
     &)
     !$cuf kernel do(2) <<<*,*>>>
     do iGRU=1,nGRU
      do iLayer=1,nSoil
        sum_rootDens = 0
        ! transpiration may be non-zero even if the soil moisture limiting factor is zero
        if (scalarTranspireLim(iGRU) > tiny(scalarTranspireLim(iGRU))) then 
          mLayerTranspireFrac(iLayer,iGRU) = mLayerRootDensity(iLayer,iGRU)*mLayerTranspireLim(iLayer,iGRU)/scalarTranspireLim(iGRU)
        else ! possibility of non-zero conductance and therefore transpiration in this case
          do iLayer2 = 1,nSoil
            sum_rootDens = sum_rootDens + mlayerRootDensity(iLayer2,iGRU)
          end do
          mLayerTranspireFrac(iLayer,iGRU) = mLayerRootDensity(iLayer,iGRU) / sum_rootDens
        end if
      if (ixBcUpperSoilHydrology==prescribedHead) then ! special case of prescribed head -- no transpiration
       mLayerTranspire(iLayer,iGRU)      = 0._rkind
       ! derivatives in transpiration w.r.t. canopy state variables
       mLayerdTrans_dCanWat(iLayer,iGRU) = 0._rkind
       mLayerdTrans_dTCanair(iLayer,iGRU)= 0._rkind
       mLayerdTrans_dTCanopy(iLayer,iGRU)= 0._rkind
       mLayerdTrans_dTGround(iLayer,iGRU)= 0._rkind
      else
       mLayerTranspire(iLayer,iGRU) = mLayerTranspireFrac(iLayer,iGRU)*scalarCanopyTranspiration(iGRU)/iden_water
       ! * derivatives in transpiration w.r.t. canopy state variables *
       mLayerdTrans_dCanWat(iLayer,iGRU)  = mLayerTranspireFrac(iLayer,iGRU)*dCanopyTrans_dCanWat(iGRU) /iden_water
       mLayerdTrans_dTCanair(iLayer,iGRU) = mLayerTranspireFrac(iLayer,iGRU)*dCanopyTrans_dTCanair(iGRU)/iden_water
       mLayerdTrans_dTCanopy(iLayer,iGRU) = mLayerTranspireFrac(iLayer,iGRU)*dCanopyTrans_dTCanopy(iGRU)/iden_water
       mLayerdTrans_dTGround(iLayer,iGRU) = mLayerTranspireFrac(iLayer,iGRU)*dCanopyTrans_dTGround(iGRU)/iden_water
      end if
    end do
  end do
     end associate
  
         ! **** Finalize operations for the fraction of transpiration loss from each soil layer *****
    ! associate(&
    !   err          => out_soilLiqFlx % err,     & ! intent(out): error code
    !   message      => out_soilLiqFlx % cmessage & ! intent(out): error message
    !  &)
      ! check fractions sum to one
      if (abs(sum(mLayerTranspireFrac) - nGRU) > verySmall) then
       message=trim(message)//'fraction transpiration in soil layers does not sum to one'
       err=20; return_flag=.true.; return
      end if
    !  end associate
   
   
   end subroutine compute_transpiration_sink
  
   subroutine compute_diagnostic_variables(nSnow_d,nGRU,decisions,mpar_data_d,flux_data_d,deriv_data_d,diag_data_d)
    ! **** compute diagnostic variables at the nodes throughout the soil profile ****
    ! type(in_type_diagv_node)  :: in_diagv_node  ! input data object for diagv_node
    ! type(out_type_diagv_node) :: out_diagv_node ! output data object for diagv_node
    integer(i4b),device :: nSnow_d(:)
    integer(i4b) :: nGRU
    type(decisions_device) :: decisions
    type(mpar_data_device),intent(in) :: mpar_data_d
    type(flux_data_device),intent(inout) :: flux_data_d
    type(deriv_data_device),intent(inout) :: deriv_data_d
    type(diag_data_device),intent(inout) :: diag_data_d
  
     call initialize_compute_diagnostic_variables()
  
     call update_compute_diagnostic_variables(nSnow_d,nGRU,decisions,mpar_data_d,flux_data_d,deriv_data_d,diag_data_d)
  
     call finalize_compute_diagnostic_variables(); if (return_flag) return
   end subroutine compute_diagnostic_variables
  
   subroutine initialize_compute_diagnostic_variables()
    ! **** Initialize operations for the compute_diagnostic_variables subroutine ****
    ! type(in_type_diagv_node),intent(out) :: in_diagv_node  ! input data object for diagv_node
    ! interface local name space to input data object for diagv_node
    ! call in_diagv_node % initialize(in_soilLiqFlx,model_decisions)
   end subroutine initialize_compute_diagnostic_variables
  
   subroutine update_compute_diagnostic_variables(nSnow_d,nGRU,decisions,mpar_data_d,flux_data_d,deriv_data_d,diag_data_d)
    ! **** Update operations for the compute_diagnostic_variables subroutine ****
    ! type(in_type_diagv_node) ,intent(in)  :: in_diagv_node  ! input data object for diagv_node
    ! type(out_type_diagv_node),intent(out) :: out_diagv_node ! output data object for diagv_node
    integer(i4b),device :: nSnow_d(:)
    integer(i4b) :: nGRU
    type(decisions_device) :: decisions
    type(mpar_data_device),intent(in) :: mpar_data_d
    type(flux_data_device),intent(inout) :: flux_data_d
    type(deriv_data_device),intent(inout) :: deriv_data_d
    type(diag_data_device),intent(inout) :: diag_data_d
  
    ! compute diagnostic variables
    call diagv_node(model_decisions(iLookDECISIONS%f_Richards)%iDecision,ixTop,min(ixBot+1,nSoil),nGRU,indx_data%nSnow,decisions,mpar_data_d,diag_data_d,flux_data_d,deriv_data_d,&
    mLayerMatricHeadLiqTrial,mLayerVolFracLiqTrial,mLayerVolFracIceTrial,&
    mLayerDiffuse,iceImpedeFac,&
    dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dTemp,dHydCond_dMatric,err,cmessage)
   end subroutine update_compute_diagnostic_variables
  
   subroutine finalize_compute_diagnostic_variables()
    ! **** Finalize operations for the compute_diagnostic_variables subroutine ****
    ! type(out_type_diagv_node),intent(in) :: out_diagv_node ! output data object for diagv_node
    ! interface output data object for diagv_node to local name space
    ! associate(&
    !  err          => out_soilLiqFlx % err,     & ! error code
    !  message      => out_soilLiqFlx % cmessage & ! error message
    ! &)
    !  call out_diagv_node % finalize(err,cmessage)
     if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
    ! end associate
   end subroutine finalize_compute_diagnostic_variables
  
   subroutine compute_surface_infiltration(nSnow_d,nSoil_d,nGRU,decisions,mpar_data_d,prog_data_d,flux_data_d,deriv_data_d,diag_data_d)
    ! **** compute infiltration at the surface and its derivative w.r.t. mass in the upper soil layer ****
    ! type(in_type_surfaceFlx)  ::  in_surfaceFlx
    ! type(out_type_surfaceFlx) :: out_surfaceFlx
    integer(i4b),device :: nSnow_d(:)
    integer(i4b) :: nSoil_d
    integer(i4b) :: nGRU
    type(decisions_device) :: decisions
    type(mpar_data_device),intent(in) :: mpar_data_d
    type(prog_data_device),intent(in) :: prog_data_d
    type(flux_data_device),intent(inout) :: flux_data_d
    type(deriv_data_device),intent(inout) :: deriv_data_d
    type(diag_data_device),intent(inout) :: diag_data_d
  
  
    call initialize_compute_surface_infiltration(nGRU,deriv_data_d)
  
    call update_compute_surface_infiltration(nSnow_d,nSoil_d,nGRU,decisions,mpar_data_d,prog_data_d,flux_data_d,deriv_data_d,diag_data_d)
  
    call finalize_compute_surface_infiltration(nGRU,flux_data_d,deriv_data_d); if (return_flag) return
  
   end subroutine compute_surface_infiltration
  
   subroutine initialize_compute_surface_infiltration(nGRU,deriv_data_d)
    ! **** Initialize operations for compute_surface_infiltration ****
    ! type(in_type_surfaceFlx),intent(out) :: in_surfaceFlx
    integer(i4b) :: nGRU
    type(deriv_data_device) :: deriv_data_d
    ! set derivative w.r.t. state above to zero (does not exist)
    associate(&
     ! intent(inout): flux derivatives ... 
     dq_dHydStateAbove => deriv_data_d%dq_dHydStateAbove_m,& ! ... in layer interfaces w.r.t. state variables in the layer above
     dq_dNrgStateAbove => deriv_data_d%dq_dNrgStateAbove_m & ! ... w.r.t. temperature in the layer above (m s-1 K-1)
    &)
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
     dq_dHydStateAbove(0,iGRU) = 0._rkind
     dq_dNrgStateAbove(0,iGRU) = 0._rkind
    end do
    end associate
  
    ! compute surface flux and its derivative...
    ! call in_surfaceFlx % initialize(nSoil,in_soilLiqFlx,model_decisions)
   end subroutine initialize_compute_surface_infiltration
  
   subroutine update_compute_surface_infiltration(nSnow_d,nSoil_d,nGRU,decisions,mpar_data_d,prog_data_d,flux_data_d,deriv_data_d,diag_data_d)
    ! **** Update operations for compute_surface_infiltration ****
    ! type(in_type_surfaceFlx) ,intent(in)    ::  in_surfaceFlx
    ! type(out_type_surfaceFlx),intent(out)   :: out_surfaceFlx
    integer(i4b),device :: nSnow_d(:)
    integer(i4b) :: nSoil_d
    integer(i4b) :: nGRU
    type(decisions_device) :: decisions
    type(mpar_data_device),intent(in) :: mpar_data_d
    type(prog_data_device),intent(in) :: prog_data_d
    type(flux_data_device),intent(inout) :: flux_data_d
    type(deriv_data_device),intent(inout) :: deriv_data_d
    type(diag_data_device),intent(inout) :: diag_data_d
  
    call surfaceFlx(firstSplitOper,nSnow_d,nSoil_d,nSoil,nRoots,ixIce,nGRU,decisions,mpar_data_d,prog_data_d,flux_data_d,deriv_data_d,diag_data_d,err,cmessage,&
    iLayerHydCond,iLayerDiffuse,mLayerMatricHeadLiqTrial,&
    mLayerVolFracLiqTrial,mLayerTempTrial,mLayerMatricHeadTrial,mLayerVolFracIceTrial,&
    dHydCond_dTemp,iceImpedeFac,iLayerHeight,&
    above_soilLiqFluxDeriv,above_soildLiq_dTk,above_soilFracLiq)
   end subroutine update_compute_surface_infiltration
  
   subroutine finalize_compute_surface_infiltration(nGRU,flux_data_d,deriv_data_d)
    ! **** Finalize operations for compute_surface_infiltration ****
    ! type(out_type_surfaceFlx),intent(in) :: out_surfaceFlx
    type(flux_data_device) :: flux_data_d
    type(deriv_data_device) :: deriv_data_d
    integer(i4b) :: nGRU
  
    ! interface object data components with local name space
    ! associate(&
    !  err     => out_soilLiqFlx % err,     & ! error code
    !  message => out_soilLiqFlx % cmessage & ! error message
    ! &)
    !  call out_surfaceFlx % finalize(err,cmessage)
     if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
    ! end associate
  
    ! include base soil evaporation as the upper boundary flux
    associate(&
     iLayerLiqFluxSoil         => flux_data_d%iLayerLiqFluxSoil_m,      & ! liquid flux at soil layer interfaces (m s-1)
     scalarGroundEvaporation   => flux_data_d%scalarGroundEvaporation,& ! ground evaporation (kg m-2 s-1)
     scalarSurfaceInfiltration => flux_data_d%scalarInfiltration,     & ! surface infiltration rate (m s-1)
     dq_dHydStateBelow         => deriv_data_d%dq_dHydStateBelow_m,      & ! derivative in the flux in layer interfaces w.r.t. state variables in the layer below
     dq_dNrgStateBelow         => deriv_data_d%dq_dNrgStateBelow_m       & ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
    &)
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
     iLayerLiqFluxSoil(0,iGRU) = scalarGroundEvaporation(iGRU)/iden_water + scalarSurfaceInfiltration(iGRU)
  
     dq_dHydStateBelow(0,iGRU) = 0._rkind ! contribution will be in dq_dHydStateLayerSurfVec(1)
     dq_dNrgStateBelow(0,iGRU) = 0._rkind ! contribution will be in dq_dNrgStateLayerSurfVec(1)
    end do
    end associate
   end subroutine finalize_compute_surface_infiltration
  
   subroutine compute_interface_fluxes_derivatives(nGRU,nSnow_d,prog_data_d,decisions,flux_data_d,deriv_data_d)
    ! **** compute fluxes and derivatives at layer interfaces ****
    integer(i4b),intent(in) :: nGRU
    integer(i4b),intent(in),device :: nSnow_d(:)
    type(prog_data_device) :: prog_data_d
    type(decisions_device) :: decisions
    type(flux_data_device) :: flux_data_d
    type(deriv_data_device) :: deriv_data_d
  
    ! computing flux at the bottom of the layer
  
     call update_compute_interface_fluxes_derivatives(nGRU,nSnow_d,prog_data_d,decisions,flux_data_d,deriv_data_d)
   
    ! end do 
   end subroutine compute_interface_fluxes_derivatives
  
   subroutine update_compute_interface_fluxes_derivatives(nGRU,nSnow_d,prog_data_d,decisions,flux_data_d,deriv_data_d)
  
  
  
    ! **** Update operations for compute_interface_fluxes_derivatives subroutine ****
    integer(i4b),intent(in) :: nGRU
    integer(i4b),intent(in),device :: nSnow_d(:)
    type(prog_data_device) :: prog_data_d
    type(decisions_device) :: decisions
    type(flux_data_device) :: flux_data_d
    type(deriv_data_device) :: deriv_data_d
  
  
    ! compute fluxes at layer interface
  
  
  
    call iLayerFlux(ixTop,min(ixBot,nSoil-1),nSnow_d,nGRU,decisions,prog_data_d,flux_data_d,deriv_data_d,&
    mLayerDiffuse,dHydCond_dTemp,dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dMatric,&
    mLayerVolFracLiqTrial,mLayerMatricHeadLiqTrial,iLayerHydCond,iLayerDiffuse)
   end subroutine update_compute_interface_fluxes_derivatives
  
   subroutine compute_drainage_flux(nGRU,nSoil_d,nSnow_d,decisions,mpar_data_d,prog_data_d,diag_data_d,flux_data_d,deriv_data_d)
    ! **** Compute the drainage flux from the bottom of the soil profile and its derivative ****
    ! type(in_type_qDrainFlux)  :: in_qDrainFlux
    ! type(out_type_qDrainFlux) :: out_qDrainFlux
    integer(i4b),intent(in) :: nGRU
    integer(i4b),device :: nSnow_d(:)
    integer(i4b) :: nSoil_d
    type(decisions_device) :: decisions
    type(mpar_data_device) :: mpar_data_d
    type(prog_data_device) :: prog_data_d
    type(diag_data_device) :: diag_data_d
    type(flux_data_device) :: flux_data_d
    type(deriv_data_device) :: deriv_data_d
  
  
    call initialize_compute_drainage_flux()
  
    call update_compute_drainage_flux(nGRU,nSoil_d,nSnow_d,decisions,mpar_data_d,prog_data_d,diag_data_d,flux_data_d,deriv_data_d)
  
    call finalize_compute_drainage_flux(nGRU,deriv_data_d,nSoil_d); if (return_flag) return
  
   end subroutine compute_drainage_flux
  
   subroutine initialize_compute_drainage_flux()
    ! **** Initialize operations for compute_drainage_flux ****
    ! type(in_type_qDrainFlux),intent(out) :: in_qDrainFlux
   end subroutine initialize_compute_drainage_flux
   
   subroutine update_compute_drainage_flux(nGRU,nSoil_d,nSnow_d,decisions,mpar_data_d,prog_data_d,diag_data_d,flux_data_d,deriv_data_d)
    ! **** Update operations for compute_drainage_flux ****
    integer(i4b),intent(in) :: nGRU
    integer(i4b),device :: nSnow_d(:)
    integer(i4b) :: nSoil_d
    type(decisions_device) :: decisions
    type(mpar_data_device) :: mpar_data_d
    type(prog_data_device) :: prog_data_d
    type(diag_data_device) :: diag_data_d
    type(flux_data_device) :: flux_data_d
    type(deriv_data_device) :: deriv_data_d
    call qDrainFlux(model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision,nSoil,nSnow_d,nGRU,decisions,mpar_data_d,prog_data_d,diag_data_d,flux_data_d,deriv_data_d,err,cmessage, iLayerHydCond, iLayerDiffuse, &
    iceImpedeFac,dHydCond_dVolLiq,dHydCond_dTemp,dHydCond_dMatric,mLayerMatricHeadLiqTrial,mLayerVolFracLiqTrial)
   end subroutine update_compute_drainage_flux
  
   subroutine finalize_compute_drainage_flux(nGRU,deriv_data_d,nSoil_d)
    ! **** finalize operations for compute_drainage_flux ****
    ! type(out_type_qDrainFlux),intent(in) :: out_qDrainFlux
    integer(i4b) :: nGRU
    type(deriv_data_device) :: deriv_data_d
    integer(i4b) :: nSoil_d
    integer(i4b) :: iGRU
    ! associate(&
    !  err     => out_soilLiqFlx % err,                       & ! error code
    !  message => out_soilLiqFlx % cmessage                   & ! error message
    ! &)
    !  call out_qDrainFlux % finalize(err,cmessage)
     if (err/=0) then; message=trim(message)//trim(cmessage); return_flag=.true.; return; end if
    ! end associate
  
    ! no dependence on the aquifer for drainage
    associate(&
     ! derivatives in flux w.r.t. ...
     dq_dHydStateBelow => deriv_data_d%dq_dHydStateBelow_m,& ! ... hydrology state variables in the layer below
     dq_dNrgStateBelow => deriv_data_d%dq_dNrgStateBelow_m & ! ... temperature in the layer below (m s-1 K-1)
    &)
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
     dq_dHydStateBelow(nSoil_d,iGRU) = 0._rkind  ! keep this here in case we want to couple some day....
     dq_dNrgStateBelow(nSoil_d,iGRU) = 0._rkind  ! keep this here in case we want to couple some day....
    end do
    end associate
   end subroutine finalize_compute_drainage_flux
  end subroutine soilLiqFlx
  
  ! ***************************************************************************************************************
  ! private subroutine diagv_node: compute transmittance and derivatives for model nodes
  ! ***************************************************************************************************************
  subroutine diagv_node(ixRichards,ixTop,ixBot,nGRU,nSnow_d,decisions,mpar_data,diag_data,flux_data,deriv_data,&
    scalarMatricHeadLiqTrial,scalarVolFracLiqTrial,scalarVolFracIceTrial,&
    scalarDiffuse,iceImpedeFac,&
    dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dTemp,dHydCond_dMatric,err,message) 
    use device_data_types
    USE soil_utils_module,only:iceImpede            ! compute the ice impedence factor
    USE soil_utils_module,only:volFracLiq           ! compute volumetric fraction of liquid water as a function of matric head
    USE soil_utils_module,only:matricHead           ! compute matric head (m)
    USE soil_utils_module,only:hydCond_psi          ! compute hydraulic conductivity as a function of matric head
    USE soil_utils_module,only:hydCond_liq          ! compute hydraulic conductivity as a function of volumetric liquid water content
    USE soil_utils_module,only:hydCondMP_liq        ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
    USE soil_utils_module,only:dTheta_dPsi_d          ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
    USE soil_utils_module,only:dPsi_dTheta_d          ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
    USE soil_utils_module,only:dPsi_dTheta2         ! compute derivative in dPsi_dTheta (m)
    USE soil_utils_module,only:dHydCond_dLiq        ! compute derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
    USE soil_utils_module,only:dHydCond_dPsi        ! compute derivative in hydraulic conductivity w.r.t. matric head (s-1)
    USE soil_utils_module,only:dIceImpede_dTemp     ! compute the derivative in the ice impedance factor w.r.t. temperature (K-1)
    ! compute hydraulic transmittance and derivatives for all layers
    implicit none
    ! input: model control, variables, derivatives, and parameters
    integer(i4b) :: ixRichards
    integer(i4b),intent(in) :: ixTop,ixBot,nGRU
    integer(i4b),intent(in),device :: nSnow_d(:)
    type(decisions_device),intent(in) :: decisions
    type(mpar_data_device),intent(in) :: mpar_data
    type(diag_data_device),intent(in) :: diag_data
    type(flux_data_device),intent(inout) :: flux_data
    type(deriv_data_device),intent(inout) :: deriv_data
    real(rkind),intent(in),device :: scalarMatricHeadLiqTrial(:,:),scalarVolFracLiqTrial(:,:),scalarVolFracIceTrial(:,:)
    real(rkind),intent(inout),device :: scalarDiffuse(:,:), iceImpedeFac(:,:)
    real(rkind),intent(inout),device :: dHydCond_dVolLiq(:,:),dDiffuse_dVolLiq(:,:),dHydCond_dTemp(:,:),dHydCond_dMatric(:,:)
    integer(i4b) :: iSoil,iGRU
    ! output: characteristic derivatives, transmittance variables, and error control
    integer(i4b) :: err
    character(*) :: message
    ! local variables
    ! real(rkind)                      :: localVolFracLiq           ! local volumetric fraction of liquid water
    ! real(rkind)                      :: scalarHydCondMP           ! hydraulic conductivity of macropores at layer mid-points (m s-1)
    ! real(rkind)                      :: dIceImpede_dT             ! derivative in ice impedance factor w.r.t. temperature (K-1)
    ! real(rkind)                      :: dHydCondMacro_dVolLiq     ! derivative in hydraulic conductivity of macropores w.r.t volumetric liquid water content (m s-1)
    ! real(rkind)                      :: dHydCondMacro_dMatric     ! derivative in hydraulic conductivity of macropores w.r.t matric head (s-1)
    ! real(rkind)                      :: dHydCondMicro_dMatric     ! derivative in hydraulic conductivity of micropores w.r.t matric head (s-1)
    ! real(rkind)                      :: dHydCondMicro_dTemp       ! derivative in hydraulic conductivity of micropores w.r.t temperature (m s-1 K-1)
    ! real(rkind)                      :: dPsi_dTheta2a             ! derivative in dPsi_dTheta (analytical)
    ! real(rkind)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
    ! real(rkind)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
    ! real(rkind)                      :: dK_dLiq__noIce            ! derivative in hydraulic conductivity w.r.t volumetric liquid water content, in the absence of ice (m s-1)
    ! real(rkind)                      :: dK_dPsi__noIce            ! derivative in hydraulic conductivity w.r.t matric head, in the absence of ice (s-1)
    ! real(rkind)                      :: relSatMP                  ! relative saturation of macropores (-)
    logical(lgt)                     :: return_flag               ! flag for return statements
  
      call initialize_diagv_node
  
      call update_diagv_node;   if (return_flag) return
  
      call finalize_diagv_node; if (return_flag) return
  
  contains
  
   subroutine initialize_diagv_node
    ! **** Initialize operations for diagv_node ****
    ! initialize error control
    return_flag=.false. 
    ! associate(&
    !   ! ixRichards => in_diagv_node % ixRichards, &
    !  err     => out_diagv_node % err    , & ! error code
    !  message => out_diagv_node % message  & ! error message
    ! &)
     err=0; message="diagv_node/"
     if (ixRichards == moisture) then
      ! haven't included macropores yet -- return with error for now
      err=20; message=trim(message)//'still need to include macropores for the moisture-based form of Richards eqn'
      return_flag=.true.; return
     end if
    ! end associate
   end subroutine initialize_diagv_node
  
   subroutine update_diagv_node
    ! **** Update operations for diagv_node ****
  
    associate(&
      ! input: model control
      ixRichards    => decisions % f_Richards, & ! index defining the option for Richards' equation (moisture or mixdform)
      ! input: state and diagnostic variables
      ! scalarMatricHeadLiqTrial => in_diagv_node % scalarMatricHeadLiqTrial, & ! liquid matric head in each layer (m)
      ! scalarVolFracLiqTrial    => in_diagv_node % scalarVolFracLiqTrial   , & ! volumetric fraction of liquid water in a given layer (-)
      ! scalarVolFracIceTrial    => in_diagv_node % scalarVolFracIceTrial   , & ! volumetric fraction of ice in a given layer (-)
      ! input: pre-computed deriavatives
      dTheta_dTk    => deriv_data%mLayerdTheta_dTk_m   , & ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
      dPsiLiq_dTemp => deriv_data%dPsiLiq_dTemp_m, & ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
      ! input: soil parameters
      vGn_alpha => mpar_data%vGn_alpha, & ! van Genuchten "alpha" parameter (m-1)
      vGn_n     => mpar_data%vGn_n    , & ! van Genuchten "n" parameter (-)
      vGn_m     => diag_data%scalarVGn_m_m    , & ! van Genuchten "m" parameter (-)
      mpExp     => mpar_data%mpExp    , & ! empirical exponent in macropore flow equation (-)
      theta_sat => mpar_data%theta_sat, & ! soil porosity (-)
      theta_res => mpar_data%theta_res, & ! soil residual volumetric water content (-)
      theta_mp  => mpar_data%theta_mp , & ! volumetric liquid water content when macropore flow begins (-)
      f_impede  => mpar_data%f_impede , & ! ice impedence factor (-) 
      ! input: saturated hydraulic conductivity ...
      scalarSatHydCond   => flux_data%mLayerSatHydCond_m,  & ! ... at the mid-point of a given layer (m s-1)
      scalarSatHydCondMP => flux_data%mLayerSatHydCondMP_m,& ! ... of macropores at the mid-point of a given layer (m s-1)
      ! output: derivative in the soil water characteristic
      scalardPsi_dTheta => deriv_data%mLayerdPsi_dTheta_m, & ! derivative in the soil water characteristic
      scalardTheta_dPsi => deriv_data%mLayerdTheta_dPsi_m, & ! derivative in the soil water characteristic
      ! output: transmittance
      scalarHydCond => flux_data%mLayerHydCond_m & ! hydraulic conductivity at layer mid-points (m s-1)
      ! scalarDiffuse => out_diagv_node % scalarDiffuse, & ! diffusivity at layer mid-points (m2 s-1)
      ! iceImpedeFac  => out_diagv_node % iceImpedeFac , & ! ice impedence factor in each layer (-)
      ! output: transmittance derivatives in ...
      ! dHydCond_dVolLiq => out_diagv_node % dHydCond_dVolLiq, & ! ... hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
      ! dDiffuse_dVolLiq => out_diagv_node % dDiffuse_dVolLiq, & ! ... hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
      ! dHydCond_dMatric => out_diagv_node % dHydCond_dMatric, & ! ... hydraulic conductivity w.r.t matric head (s-1)
      ! dHydCond_dTemp   => out_diagv_node % dHydCond_dTemp,   & ! ... hydraulic conductivity w.r.t temperature (m s-1 K-1)
  
      ! output: error control
      ! err     => out_diagv_node % err    , & ! error code
      ! message => out_diagv_node % message  & ! error message
     &)
   
     !$cuf kernel do(2) <<<*,*>>> 
     do iGRU=1,nGRU
     do iSoil=ixTop,ixBot
    call calc_diagv_node_update(ixRichards,&
    scalarMatricHeadLiqTrial(iSoil,iGRU),scalarVolFracLiqTrial(iSoil+nSnow_d(iGRU),iGRU),scalarVolFracIceTrial(iSoil+nSnow_d(iGRU),iGRU),&
    dTheta_dTk(iSoil+nSnow_d(iGRU),iGRU),dPsiLiq_dTemp(iSoil,iGRU),&
    vGn_alpha(iSoil),vGn_n(iSoil),vGn_m(iSoil,iGRU),mpExp,theta_sat(iSoil),theta_res(iSoil),theta_mp,f_impede,&
    scalarSatHydCond(iSoil,iGRU),scalarSatHydCondMP(iSoil,iGRU),&
    scalardPsi_dTheta(iSoil,iGRU),scalardTheta_dPsi(iSoil,iGRU),&
    scalarHydCond(iSoil,iGRU),scalarDiffuse(iSoil,iGRU),iceImpedeFac(iSoil,iGRU),&
    dHydCond_dVolLiq(iSoil,iGRU),dDiffuse_dVolLiq(iSoil,iGRU),dHydCond_dMatric(iSoil,iGRU),dHydCond_dTemp(iSoil,iGRU))
    end do
  end do
  
    end associate
  
   end subroutine update_diagv_node
  
   subroutine finalize_diagv_node
    ! **** Finalize operations for diagv_node ****
    ! associate(&
    !  deriv_desired => in_diagv_node % deriv_desired        & ! flag indicating if derivatives are desired
    !  dHydCond_dVolLiq => out_diagv_node % dHydCond_dVolLiq, & ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
    !  dDiffuse_dVolLiq => out_diagv_node % dDiffuse_dVolLiq, & ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
    !  dHydCond_dMatric => out_diagv_node % dHydCond_dMatric  & ! derivative in hydraulic conductivity w.r.t matric head (s-1)
    ! &)
     ! if derivatives are not desired, then set values to missing
    !  if (.not.deriv_desired) then
    !    dHydCond_dVolLiq   = realMissing ! not used, so cause problems
    !    dDiffuse_dVolLiq   = realMissing ! not used, so cause problems
    !    dHydCond_dMatric   = realMissing ! not used, so cause problems
    !  end if
    ! end associate
  
    ! associate(&
    !  err     => out_diagv_node % err    , & ! error code
    !  message => out_diagv_node % message  & ! error message
    ! &)
     ! final error check
     if (err /= 0_i4b) then
      message=trim(message)//'unanticipated error in diagv_node'
      return_flag=.true.; return
     end if
    ! end associate
   end subroutine finalize_diagv_node
  
  end subroutine diagv_node
  
  attributes(device) subroutine calc_diagv_node_update(ixRichards,&
    scalarMatricHeadLiqTrial,scalarVolFracLiqTrial,scalarVolFracIceTrial,&
    dTheta_dTk,dPsiLiq_dTemp,&
    vGn_alpha,vGn_n,vGn_m,mpExp,theta_sat,theta_res,theta_mp,f_impede,&
    scalarSatHydCond,scalarSatHydCondMP,&
    scalardPsi_dTheta,scalardTheta_dPsi,&
    scalarHydCond,scalarDiffuse,iceImpedeFac,&
    dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dMatric,dHydCond_dTemp)
    USE soil_utils_module,only:iceImpede            ! compute the ice impedence factor
    USE soil_utils_module,only:volFracLiq           ! compute volumetric fraction of liquid water as a function of matric head
    USE soil_utils_module,only:matricHead           ! compute matric head (m)
    USE soil_utils_module,only:hydCond_psi          ! compute hydraulic conductivity as a function of matric head
    USE soil_utils_module,only:hydCond_liq          ! compute hydraulic conductivity as a function of volumetric liquid water content
    USE soil_utils_module,only:hydCondMP_liq        ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
    USE soil_utils_module,only:dTheta_dPsi          ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
    USE soil_utils_module,only:dPsi_dTheta          ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
    USE soil_utils_module,only:dPsi_dTheta2         ! compute derivative in dPsi_dTheta (m)
    USE soil_utils_module,only:dHydCond_dLiq        ! compute derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
    USE soil_utils_module,only:dHydCond_dPsi        ! compute derivative in hydraulic conductivity w.r.t. matric head (s-1)
    USE soil_utils_module,only:dIceImpede_dTemp     ! compute the derivative in the ice impedance factor w.r.t. temperature (K-1)
  
    integer(i4b) :: ixRichards
    real(rkind) :: scalarMatricHeadLiqTrial,scalarVolFracLiqTrial,scalarVolFracIceTrial
    real(rkind) :: dTheta_dTk,dPsiLiq_dTemp
    real(rkind) :: vGn_alpha,vGn_n,vGn_m,mpExp,theta_sat,theta_res,theta_mp,f_impede
    real(rkind) :: scalarSatHydCond,scalarSatHydCondMP
    real(rkind) :: scalardPsi_dTheta,scalardTheta_dPsi
    real(rkind) :: scalarHydCond,scalarDiffuse,iceImpedeFac
    real(rkind) :: dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dMatric,dHydCond_dTemp
  
    real(rkind)                      :: localVolFracLiq           ! local volumetric fraction of liquid water
    real(rkind)                      :: scalarHydCondMP           ! hydraulic conductivity of macropores at layer mid-points (m s-1)
    real(rkind)                      :: dIceImpede_dT             ! derivative in ice impedance factor w.r.t. temperature (K-1)
    real(rkind)                      :: dHydCondMacro_dVolLiq     ! derivative in hydraulic conductivity of macropores w.r.t volumetric liquid water content (m s-1)
    real(rkind)                      :: dHydCondMacro_dMatric     ! derivative in hydraulic conductivity of macropores w.r.t matric head (s-1)
    real(rkind)                      :: dHydCondMicro_dMatric     ! derivative in hydraulic conductivity of micropores w.r.t matric head (s-1)
    real(rkind)                      :: dHydCondMicro_dTemp       ! derivative in hydraulic conductivity of micropores w.r.t temperature (m s-1 K-1)
    real(rkind)                      :: dPsi_dTheta2a             ! derivative in dPsi_dTheta (analytical)
    real(rkind)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
    real(rkind)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
    real(rkind)                      :: dK_dLiq__noIce            ! derivative in hydraulic conductivity w.r.t volumetric liquid water content, in the absence of ice (m s-1)
    real(rkind)                      :: dK_dPsi__noIce            ! derivative in hydraulic conductivity w.r.t matric head, in the absence of ice (s-1)
    real(rkind)                      :: relSatMP                  ! relative saturation of macropores (-)
  
    select case(ixRichards)
      case(moisture)
        scalardPsi_dTheta = dPsi_dTheta(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
        scalardTheta_dPsi = realMissing  ! deliberately cause problems if this is ever used
      case(mixdform)
        scalardTheta_dPsi = dTheta_dPsi(scalarMatricHeadLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
        scalardPsi_dTheta = dPsi_dTheta(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
    end select
    ! compute the ice impedence factor and its derivative w.r.t. volumetric liquid water content (-)
    call iceImpede(scalarVolFracIceTrial,f_impede, &  ! input
                    iceImpedeFac,dIceImpede_dLiq)     ! output
    select case(ixRichards)
      case(moisture) ! moisture-based form of Richards' equation
        call hydraulic_conductivity_moisture(scalarVolFracLiqTrial,scalarVolFracIceTrial,&
        vGn_alpha,vGn_n,vGn_m,theta_res,theta_sat,&
        scalarSatHydCond,scalardPsi_dTheta,&
        scalarHydCond,scalarDiffuse,iceImpedeFac,&
        dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dMatric,&
        hydCond_noIce,dK_dLiq__noIce,dIceImpede_dLiq,dPsi_dTheta2a); 
      case(mixdform) ! mixed form of Richards' equation -- just compute hydraulic condictivity
        call hydraulic_conductivity_mixed_form(scalarMatricHeadLiqTrial,scalarVolFracIceTrial,&
        dTheta_dTk,dPsiLiq_dTemp,&
        vGn_alpha,vGn_n,vGn_m,mpExp,theta_sat,theta_res,theta_mp,f_impede,&
        scalarSatHydCond,scalarSatHydCondMP,&
        scalardTheta_dPsi,&
        scalarHydCond,scalarDiffuse,iceImpedeFac,&
        dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dMatric,dHydCond_dTemp,&
        hydCond_noIce,localVolFracLiq,scalarHydCondMP,&
        relSatMP,dHydCondMacro_dVolLiq,dHydCondMacro_dMatric,&
        dK_dPsi__noIce,dHydCondMicro_dTemp,dHydCondMicro_dMatric,&
        dIceImpede_dLiq,dIceImpede_dT);    
      ! case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
    end select 
  
  end subroutine calc_diagv_node_update
  
  attributes(device) subroutine hydraulic_conductivity_moisture(scalarVolFracLiqTrial,scalarVolFracIceTrial,&
    vGn_alpha,vGn_n,vGn_m,theta_res,theta_sat,&
    scalarSatHydCond,scalardPsi_dTheta,&
    scalarHydCond,scalarDiffuse,iceImpedeFac,&
    dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dMatric,&
    hydCond_noIce,dK_dLiq__noIce,dIceImpede_dLiq,dPsi_dTheta2a)
  
    USE soil_utils_module,only:iceImpede            ! compute the ice impedence factor
    USE soil_utils_module,only:volFracLiq           ! compute volumetric fraction of liquid water as a function of matric head
    USE soil_utils_module,only:matricHead           ! compute matric head (m)
    USE soil_utils_module,only:hydCond_psi          ! compute hydraulic conductivity as a function of matric head
    USE soil_utils_module,only:hydCond_liq          ! compute hydraulic conductivity as a function of volumetric liquid water content
    USE soil_utils_module,only:hydCondMP_liq        ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
    USE soil_utils_module,only:dTheta_dPsi          ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
    USE soil_utils_module,only:dPsi_dTheta          ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
    USE soil_utils_module,only:dPsi_dTheta2         ! compute derivative in dPsi_dTheta (m)
    USE soil_utils_module,only:dHydCond_dLiq        ! compute derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
    USE soil_utils_module,only:dHydCond_dPsi        ! compute derivative in hydraulic conductivity w.r.t. matric head (s-1)
    USE soil_utils_module,only:dIceImpede_dTemp     ! compute the derivative in the ice impedance factor w.r.t. temperature (K-1)
  
    real(rkind),intent(in) :: scalarVolFracLiqTrial,scalarVolFracIceTrial
  real(rkind),intent(in) :: vGn_alpha,vGn_n,vGn_m,theta_sat,theta_res
  real(rkind),intent(in) :: scalarSatHydCond,scalardPsi_dTheta
  real(rkind),intent(inout) :: scalarHydCond,scalarDiffuse,iceImpedeFac
  real(rkind),intent(inout) :: dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dMatric
  ! real(rkind)                      :: localVolFracLiq           ! local volumetric fraction of liquid water
  ! real(rkind)                      :: scalarHydCondMP           ! hydraulic conductivity of macropores at layer mid-points (m s-1)
  ! real(rkind)                      :: dIceImpede_dT             ! derivative in ice impedance factor w.r.t. temperature (K-1)
  ! real(rkind)                      :: dHydCondMacro_dVolLiq     ! derivative in hydraulic conductivity of macropores w.r.t volumetric liquid water content (m s-1)
  ! real(rkind)                      :: dHydCondMacro_dMatric     ! derivative in hydraulic conductivity of macropores w.r.t matric head (s-1)
  ! real(rkind)                      :: dHydCondMicro_dMatric     ! derivative in hydraulic conductivity of micropores w.r.t matric head (s-1)
  ! real(rkind)                      :: dHydCondMicro_dTemp       ! derivative in hydraulic conductivity of micropores w.r.t temperature (m s-1 K-1)
  real(rkind)                      :: dPsi_dTheta2a             ! derivative in dPsi_dTheta (analytical)
  real(rkind)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
  real(rkind)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
  real(rkind)                      :: dK_dLiq__noIce            ! derivative in hydraulic conductivity w.r.t volumetric liquid water content, in the absence of ice (m s-1)
  ! real(rkind)                      :: dK_dPsi__noIce            ! derivative in hydraulic conductivity w.r.t matric head, in the absence of ice (s-1)
  ! real(rkind)                      :: relSatMP                  ! relative saturation of macropores (-)
     ! compute the hydraulic conductivity (m s-1) and diffusivity (m2 s-1) for a given layer
    hydCond_noIce = hydCond_liq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m)
    scalarHydCond = hydCond_noIce*iceImpedeFac
    scalarDiffuse = scalardPsi_dTheta * scalarHydCond
    ! compute derivative in hydraulic conductivity (m s-1) and hydraulic diffusivity (m2 s-1)
      if (scalarVolFracIceTrial > epsilon(iceImpedeFac)) then
        dK_dLiq__noIce   = dHydCond_dLiq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m,.true.)  ! [.true. = analytical]
        dHydCond_dVolLiq = hydCond_noIce*dIceImpede_dLiq + dK_dLiq__noIce*iceImpedeFac
      else
        dHydCond_dVolLiq = dHydCond_dLiq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m,.true.)
      end if
        dPsi_dTheta2a    = dPsi_dTheta2(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m,.true.)   ! [.true. = analytical] compute derivative in dPsi_dTheta (m)
        dDiffuse_dVolLiq = dHydCond_dVolLiq*scalardPsi_dTheta + scalarHydCond*dPsi_dTheta2a
        dHydCond_dMatric = realMissing ! not used, so cause problems
  end subroutine
  
  attributes(device) subroutine hydraulic_conductivity_mixed_form(scalarMatricHeadLiqTrial,scalarVolFracIceTrial,&
        dTheta_dTk,dPsiLiq_dTemp,&
        vGn_alpha,vGn_n,vGn_m,mpExp,theta_sat,theta_res,theta_mp,f_impede,&
        scalarSatHydCond,scalarSatHydCondMP,&
        scalardTheta_dPsi,&
        scalarHydCond,scalarDiffuse,iceImpedeFac,&
        dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dMatric,dHydCond_dTemp,&
        hydCond_noIce,localVolFracLiq,scalarHydCondMP,&
        relSatMP,dHydCondMacro_dVolLiq,dHydCondMacro_dMatric,&
        dK_dPsi__noIce,dHydCondMicro_dTemp,dHydCondMicro_dMatric,&
        dIceImpede_dLiq,dIceImpede_dT)
        USE soil_utils_module,only:iceImpede            ! compute the ice impedence factor
        USE soil_utils_module,only:volFracLiq           ! compute volumetric fraction of liquid water as a function of matric head
        USE soil_utils_module,only:matricHead           ! compute matric head (m)
        USE soil_utils_module,only:hydCond_psi          ! compute hydraulic conductivity as a function of matric head
        USE soil_utils_module,only:hydCond_liq          ! compute hydraulic conductivity as a function of volumetric liquid water content
        USE soil_utils_module,only:hydCondMP_liq        ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
        USE soil_utils_module,only:dTheta_dPsi          ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
        USE soil_utils_module,only:dPsi_dTheta          ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
        USE soil_utils_module,only:dPsi_dTheta2         ! compute derivative in dPsi_dTheta (m)
        USE soil_utils_module,only:dHydCond_dLiq        ! compute derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
        USE soil_utils_module,only:dHydCond_dPsi        ! compute derivative in hydraulic conductivity w.r.t. matric head (s-1)
        USE soil_utils_module,only:dIceImpede_dTemp     ! compute the derivative in the ice impedance factor w.r.t. temperature (K-1)
      
    real(rkind),intent(in) :: scalarMatricHeadLiqTrial,scalarVolFracIceTrial
    real(rkind),intent(in) :: dTheta_dTk,dPsiLiq_dTemp
    real(rkind),intent(in) :: vGn_alpha,vGn_n,vGn_m,mpExp,theta_sat,theta_res,theta_mp,f_impede
    real(rkind),intent(in) :: scalarSatHydCond,scalarSatHydCondMP
    real(rkind),intent(inout) :: scalardTheta_dPsi
    real(rkind),intent(inout) :: scalarHydCond,scalarDiffuse,iceImpedeFac
    real(rkind),intent(inout) :: dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dMatric,dHydCond_dTemp
    real(rkind)                      :: localVolFracLiq           ! local volumetric fraction of liquid water
    real(rkind)                      :: scalarHydCondMP           ! hydraulic conductivity of macropores at layer mid-points (m s-1)
    real(rkind)                      :: dIceImpede_dT             ! derivative in ice impedance factor w.r.t. temperature (K-1)
    real(rkind)                      :: dHydCondMacro_dVolLiq     ! derivative in hydraulic conductivity of macropores w.r.t volumetric liquid water content (m s-1)
    real(rkind)                      :: dHydCondMacro_dMatric     ! derivative in hydraulic conductivity of macropores w.r.t matric head (s-1)
    real(rkind)                      :: dHydCondMicro_dMatric     ! derivative in hydraulic conductivity of micropores w.r.t matric head (s-1)
    real(rkind)                      :: dHydCondMicro_dTemp       ! derivative in hydraulic conductivity of micropores w.r.t temperature (m s-1 K-1)
    ! real(rkind)                      :: dPsi_dTheta2a             ! derivative in dPsi_dTheta (analytical)
    real(rkind)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
    real(rkind)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
    ! real(rkind)                      :: dK_dLiq__noIce            ! derivative in hydraulic conductivity w.r.t volumetric liquid water content, in the absence of ice (m s-1)
    real(rkind)                      :: dK_dPsi__noIce            ! derivative in hydraulic conductivity w.r.t matric head, in the absence of ice (s-1)
    real(rkind)                      :: relSatMP                  ! relative saturation of macropores (-)
    
    ! compute the hydraulic conductivity (m s-1) and diffusivity (m2 s-1) for a given layer
   hydCond_noIce = hydCond_psi(scalarMatricHeadLiqTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m)
   scalarDiffuse = realMissing ! not used, so cause problems
   ! compute the hydraulic conductivity of macropores (m s-1)
   localVolFracLiq = volFracLiq(scalarMatricHeadLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   scalarHydCondMP = hydCondMP_liq(localVolFracLiq,theta_sat,theta_mp,mpExp,scalarSatHydCondMP,scalarSatHydCond)
       

   scalarHydCond   = hydCond_noIce*iceImpedeFac + scalarHydCondMP
   ! compute derivative in hydraulic conductivity (m s-1)
     ! compute derivative for macropores
     if (localVolFracLiq > theta_mp) then
       relSatMP              = (localVolFracLiq - theta_mp)/(theta_sat - theta_mp)
       dHydCondMacro_dVolLiq = ((scalarSatHydCondMP - scalarSatHydCond)/(theta_sat - theta_mp))*mpExp*(relSatMP**(mpExp - 1._rkind))
       dHydCondMacro_dMatric = scalardTheta_dPsi*dHydCondMacro_dVolLiq
     else
       dHydCondMacro_dVolLiq = 0._rkind
       dHydCondMacro_dMatric = 0._rkind
     end if
     ! compute derivatives for micropores
     if (scalarVolFracIceTrial > verySmall) then
       dK_dPsi__noIce        = dHydCond_dPsi(scalarMatricHeadLiqTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m,.true.)  ! analytical
       dHydCondMicro_dTemp   = dPsiLiq_dTemp*dK_dPsi__noIce  ! m s-1 K-1
       dHydCondMicro_dMatric = hydCond_noIce*dIceImpede_dLiq*scalardTheta_dPsi + dK_dPsi__noIce*iceImpedeFac
     else
       dHydCondMicro_dTemp   = 0._rkind
       dHydCondMicro_dMatric = dHydCond_dPsi(scalarMatricHeadLiqTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m,.true.)
     end if
     ! combine derivatives
     dHydCond_dMatric = dHydCondMicro_dMatric + dHydCondMacro_dMatric
  
     ! compute analytical derivative for change in ice impedance factor w.r.t. temperature
     call dIceImpede_dTemp(scalarVolFracIceTrial, & ! intent(in):  trial value of volumetric ice content (-)
                           dTheta_dTk,            & ! intent(in):  derivative in volumetric liquid water content w.r.t. temperature (K-1)
                           f_impede,              & ! intent(in):  ice impedance parameter (-)
                           dIceImpede_dT          ) ! intent(out): derivative in ice impedance factor w.r.t. temperature (K-1)
     ! compute derivative in hydraulic conductivity w.r.t. temperature
     dHydCond_dTemp = hydCond_noIce*dIceImpede_dT + dHydCondMicro_dTemp*iceImpedeFac
     ! set values that are not used to missing
     dHydCond_dVolLiq = realMissing ! not used, so cause problems
     dDiffuse_dVolLiq = realMissing ! not used, so cause problems
  end subroutine
  
  ! ***************************************************************************************************************
  ! private subroutine surfaceFlx: compute the surface flux and its derivative
  ! ***************************************************************************************************************
  subroutine surfaceFlx(firstSplitOper,nSnow,nSoil_d,nSoil,nRoots,ixIce,nGRU,decisions,mpar_data,prog_data,flux_data,deriv_data,diag_data,err,message,surfaceHydCond,surfaceDiffuse,scalarMatricHeadLiq,mLayerVolFracLiq,&
    mLayerTemp,mLayerMatricHead,mLayerVolFracIce,&
    dHydCond_dTemp,iceImpedeFac,iLayerHeight,&
    above_soilLiqFluxDeriv,above_soildLiq_dTk,above_soilFracLiq)
    USE soil_utils_module,only:volFracLiq            ! compute volumetric fraction of liquid water as a function of matric head (-)
    USE soil_utils_module,only:hydCond_psi           ! compute hydraulic conductivity as a function of matric head (m s-1)
    USE soil_utils_module,only:hydCond_liq           ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
    USE soil_utils_module,only:dPsi_dTheta           ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
    USE soil_utils_module,only:crit_soilT            ! compute critical temperature below which ice exists
    USE soil_utils_module,only:gammp                 ! compute the cumulative probabilty based on the Gamma distribution
    use device_data_types
    ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
    implicit none
    ! -----------------------------------------------------------------------------------------------------------------------------
    ! input: model control, variables, derivatives, soil layer depth, boundary conditions, fluxes, and transmittance and soil parameters
    logical(lgt) :: firstSplitOper          ! input object for surfaceFlx
    integer(i4b),device :: nSnow(:)
    integer(i4b) :: nSoil,nSoil_d
    integer(i4b),device :: nRoots(:),ixIce(:)
    integer(i4b) :: nGRU
    ! input-output: hydraulic conductivity and diffusivity, and infiltration parameters
    type(decisions_device) :: decisions
    type(mpar_data_device),intent(in) :: mpar_data
    type(prog_data_device),intent(in) :: prog_data
    type(flux_data_device),intent(inout) :: flux_data
    type(deriv_data_device),intent(inout) :: deriv_data
    type(diag_data_device),intent(inout) :: diag_data
    ! output: runoff, infiltration, derivatives, and error control
    integer(i4b) :: err
    character(*) :: message
    real(rkind),device :: surfaceHydCond(0:,:),surfaceDiffuse(0:,:)
    real(rkind),device :: scalarMatricHeadLiq(:,:), mLayerVolFracLiq(:,:)
    real(rkind),device :: mLayerTemp(:,:),mLayerMatricHead(:,:),mLayerVolFracIce(:,:)
    real(rkind),device :: dHydCond_dTemp(:,:), iceImpedeFac(:,:), iLayerHeight(0:,:)
    real(rkind),device :: above_soilLiqFluxDeriv(:),above_soildLiq_dTk(:),above_soilFracLiq(:)
    ! -----------------------------------------------------------------------------------------------------------------------------
    ! local variables
    ! general
    integer(i4b),device                     :: iLayer                              ! index of soil layer
    ! real(rkind)                      :: Tcrit                               ! temperature where all water is unfrozen (K)
    ! real(rkind)                      :: fPart1,fPart2                       ! different parts of a function
    ! real(rkind)                      :: dPart1(1:in_surfaceFlx % nSoil)     ! derivatives for different parts of a function
    ! real(rkind)                      :: dPart2(1:in_surfaceFlx % nSoil)     ! derivatives for different parts of a function
    ! real(rkind)                      :: dfracCap(1:in_surfaceFlx % nSoil)   ! derivatives for different parts of a function
    ! real(rkind)                      :: dfInfRaw(1:in_surfaceFlx % nSoil)   ! derivatives for different parts of a function
    ! real(rkind)                      :: total_soil_depth                    ! total depth of soil (m)
    ! head boundary condition
    integer(i4b) :: iGRU
    real(rkind)                      :: cFlux                               ! capillary flux (m s-1)
    real(rkind)                      :: dNum                                ! numerical derivative
    ! simplified Green-Ampt infiltration
    ! real(rkind)                      :: rootZoneLiq                         ! depth of liquid water in the root zone (m)
    ! real(rkind)                      :: rootZoneIce                         ! depth of ice in the root zone (m)
    ! real(rkind)                      :: availCapacity                       ! available storage capacity in the root zone (m)
    ! real(rkind)                      :: depthWettingFront                   ! depth to the wetting front (m)
    ! real(rkind)                      :: hydCondWettingFront                 ! hydraulic conductivity at the wetting front (m s-1)
    ! saturated area associated with variable storage capacity
    ! real(rkind)                      :: fracCap                             ! fraction of pore space filled with liquid water and ice (-)
    ! real(rkind)                      :: fInfRaw                             ! infiltrating area before imposing solution constraints (-)
    ! real(rkind),parameter            :: maxFracCap=0.995_rkind              ! maximum fraction capacity -- used to avoid numerical problems associated with an enormous derivative
    ! real(rkind),parameter            :: scaleFactor=0.000001_rkind          ! scale factor for the smoothing function (-)
    ! real(rkind),parameter            :: qSurfScaleMax=1000._rkind           ! maximum surface runoff scaling factor (-)
    ! fraction of impermeable area associated with frozen ground
    ! real(rkind)                      :: alpha                               ! shape parameter in the Gamma distribution
    ! real(rkind)                      :: xLimg                               ! upper limit of the integral
    ! derivatives in ...
    ! real(rkind) :: dVolFracLiq_dWat(1:in_surfaceFlx % nSoil)  ! ... vol fraction of liquid w.r.t. water state variable in root layers
    ! real(rkind) :: dVolFracIce_dWat(1:in_surfaceFlx % nSoil)  ! ... vol fraction of ice w.r.t. water state variable in root layers
    ! real(rkind) :: dVolFracLiq_dTk(1:in_surfaceFlx % nSoil)   ! ... vol fraction of liquid w.r.t. temperature in root layers
    ! real(rkind) :: dVolFracIce_dTk(1:in_surfaceFlx % nSoil)   ! ... vol fraction of ice w.r.t. temperature in root layers
    ! real(rkind) :: dRootZoneLiq_dWat(1:in_surfaceFlx % nSoil) ! ... vol fraction of scalar root zone liquid w.r.t. water state variable in root layers
    ! real(rkind) :: dRootZoneIce_dWat(1:in_surfaceFlx % nSoil) ! ... vol fraction of scalar root zone ice w.r.t. water state variable in root layers
    ! real(rkind) :: dRootZoneLiq_dTk(1:in_surfaceFlx % nSoil)  ! ... vol fraction of scalar root zone liquid w.r.t. temperature in root layers
    ! real(rkind) :: dRootZoneIce_dTk(1:in_surfaceFlx % nSoil)  ! ... vol fraction of scalar root zone ice w.r.t. temperature in root layers
    ! real(rkind) :: dDepthWettingFront_dWat(1:in_surfaceFlx % nSoil) ! ... scalar depth of wetting front w.r.t. water state variable in root layers
    ! real(rkind) :: dDepthWettingFront_dTk(1:in_surfaceFlx % nSoil)  ! ... scalar depth of wetting front w.r.t. temperature in root layers
    ! real(rkind) :: dxMaxInfilRate_dWat(1:in_surfaceFlx % nSoil) ! ... scalar max infiltration rate w.r.t. water state variable in root layers
    ! real(rkind) :: dxMaxInfilRate_dTk(1:in_surfaceFlx % nSoil)  ! ... scalar max infiltration rate w.r.t. temperature in root layers
    ! real(rkind) :: dInfilArea_dWat(0:in_surfaceFlx % nSoil)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
    ! real(rkind) :: dInfilArea_dTk(0:in_surfaceFlx % nSoil)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
    ! real(rkind) :: dFrozenArea_dWat(0:in_surfaceFlx % nSoil) ! ... scalar frozen area w.r.t. water state variable in canopy or snow and root layers
    ! real(rkind) :: dFrozenArea_dTk(0:in_surfaceFlx % nSoil)  ! ... scalar frozen area w.r.t. temperature in canopy or snow and root layers
    ! real(rkind) :: dInfilRate_dWat(0:in_surfaceFlx % nSoil)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
    ! real(rkind) :: dInfilRate_dTk(0:in_surfaceFlx % nSoil)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
    ! error control
    logical(lgt) :: return_flag ! logical flag for return statements
  
    call initialize_surfaceFlx
  
    call update_surfaceFlx;   if (return_flag) return
  
    call finalize_surfaceFlx; if (return_flag) return
  
  contains
  
   subroutine initialize_surfaceFlx
    ! **** Initialize operations for surfaceFlx ****
   
    ! allocate output object array components
    ! out_surfaceFlx % dq_dHydStateVec = deriv_data%var(iLookDERIV%dq_dHydStateLayerSurfVec)%dat
    ! out_surfaceFlx % dq_dNrgStateVec = deriv_data%var(iLookDERIV%dq_dNrgStateLayerSurfVec)%dat
  
    ! initialize error control
    return_flag=.false.
    ! associate(&
    !  err     => out_surfaceFlx % err    , & ! error code
    !  message => out_surfaceFlx % message  & ! error message
    ! &)
     err=0; message="surfaceFlx/"
    ! end associate
   
    ! initialize derivatives
    associate(&
     ! output: derivatives in surface infiltration w.r.t. ...
     dq_dHydStateVec => deriv_data%dq_dHydStateLayerSurfVec_m , & ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
     dq_dNrgStateVec => deriv_data%dq_dNrgStateLayerSurfVec_m   & ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
    &)
     dq_dHydStateVec = 0._rkind
     dq_dNrgStateVec = 0._rkind
    end associate
   end subroutine initialize_surfaceFlx
  
   subroutine update_surfaceFlx
    ! **** Update operations for surfaceFlx ****
    associate(&
     ! input: model control
     bc_upper => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision & ! index defining the type of boundary conditions
     ! output: error control
    !  err      => out_surfaceFlx % err    , & ! error code
    !  message  => out_surfaceFlx % message  & ! error message
    &)
  
     ! compute the surface flux and its derivative
     select case(bc_upper)
   
       case(prescribedHead) ! head condition
         call update_surfaceFlx_prescribedHead; if (return_flag) return 
   
       case(liquidFlux) ! flux condition
         call update_surfaceFlx_liquidFlux;     if (return_flag) return 
   
       case default; err=20; message=trim(message)//'unknown upper boundary condition for soil hydrology'; return_flag=.true.; return
   
     end select 
  
    end associate
   end subroutine update_surfaceFlx
  
   subroutine update_surfaceFlx_prescribedHead
    ! **** Update operations for surfaceFlx: prescribed pressure head condition ****
    associate(&
     ! input: model control
    !  deriv_desired  => in_surfaceFlx % deriv_desired  , & ! flag to indicate if derivatives are desired
     ixRichards     => decisions % f_Richards     , & ! index defining the option for Richards' equation (moisture or mixdform)
     ! input: state and diagnostic variables
    !  scalarMatricHeadLiq => in_surfaceFlx % scalarMatricHeadLiq , & ! liquid matric head in the upper-most soil layer (m)
    !  scalarVolFracLiq    => in_surfaceFlx % scalarVolFracLiq    , & ! volumetric liquid water content in the upper-most soil layer (-)
     ! input: depth of upper-most soil layer (m)
     mLayerDepth  => prog_data%mLayerDepth  , & ! depth of upper-most soil layer (m)
     ! input: diriclet boundary conditions
     upperBoundHead   => mpar_data%upperBoundHead  , & ! upper boundary condition for matric head (m)
     upperBoundTheta  => mpar_data%upperBoundTheta , & ! upper boundary condition for volumetric liquid water content (-)
     ! input: transmittance
     surfaceSatHydCond => flux_data%iLayerSatHydCond_m , & ! saturated hydraulic conductivity at the surface (m s-1)
    !  dHydCond_dTemp    => in_surfaceFlx % dHydCond_dTemp    , & ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
    !  iceImpedeFac      => in_surfaceFlx % iceImpedeFac      , & ! ice impedence factor in the upper-most soil layer (-)
     ! input: soil parameters
     vGn_alpha           => mpar_data%vGn_alpha           , & ! van Genuchten "alpha" parameter (m-1)
     vGn_n               => mpar_data%vGn_n               , & ! van Genuchten "n" parameter (-)
     vGn_m               => diag_data%scalarVGn_m_m               , & ! van Genuchten "m" parameter (-)
     theta_sat           => mpar_data%theta_sat           , & ! soil porosity (-)
     theta_res           => mpar_data%theta_res           , & ! soil residual volumetric water content (-)
     ! input-output: hydraulic conductivity and diffusivity at the surface
     ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
    !  surfaceHydCond => io_surfaceFlx % surfaceHydCond , & ! hydraulic conductivity (m s-1)
    !  surfaceDiffuse => io_surfaceFlx % surfaceDiffuse , & ! hydraulic diffusivity at the surface (m
     ! output: runoff and infiltration
     scalarSurfaceRunoff       => flux_data%scalarSurfaceRunoff       , & ! surface runoff (m s-1)
     scalarSurfaceInfiltration => flux_data%scalarInfiltration , & ! surface infiltration (m s-1)
     ! output: derivatives in surface infiltration w.r.t. ...
     dq_dHydStateVec => deriv_data%dq_dHydStateLayerSurfVec_m , & ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
     dq_dNrgStateVec => deriv_data%dq_dNrgStateLayerSurfVec_m  & ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
     ! output: error control
    !  err     => out_surfaceFlx % err    , & ! error code
    !  message => out_surfaceFlx % message  & ! error message
    &)
  
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
    call surfaceFlx_prescribedHead(scalarSurfaceRunoff(iGRU),surfaceHydCond(0,iGRU),surfaceDiffuse(0,iGRU),&
    scalarMatricHeadLiq(1,iGRU), mLayerVolFracLiq(nSnow(iGRU)+1,iGRU),&
    mLayerDepth(1+nSnow(iGRU),iGRU),dq_dHydStateVec(1,iGRU),dq_dNrgStateVec(1,iGRU),&
    upperBoundHead,upperBoundTheta,&
    surfaceSatHydCond(0,iGRU),dHydCond_dTemp(1,iGRU),iceImpedeFac(1,iGRU),&
    vGn_alpha(1),vGn_n(1),vGn_m(1,iGRU),theta_sat(1),theta_res(1),&
    scalarSurfaceInfiltration(iGRU),&
    ixRichards)
    end do
  
    end associate
   end subroutine update_surfaceFlx_prescribedHead
  
   subroutine update_surfaceFlx_liquidFlux 
    ! **** Update operations for surfaceFlx: flux condition ****
  
    ! force infiltration to be constant over the iterations
    associate(&
     ! input: model control
    !  firstSplitOper => in_surfaceFlx % firstSplitOper , & ! flag indicating if desire to compute infiltration
     ! output: derivatives in surface infiltration w.r.t. ...
     dq_dHydStateVec => deriv_data%dq_dHydStateLayerSurfVec_m , & ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
     dq_dNrgStateVec => deriv_data%dq_dNrgStateLayerSurfVec_m   & ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
    &)
     if (firstSplitOper) then
       call update_surfaceFlx_liquidFlux_computation; if (return_flag) return 
     else ! do not compute infiltration after first flux call in a splitting operation
       dq_dHydStateVec = 0._rkind
       dq_dNrgStateVec = 0._rkind
     end if 
    end associate
  
    call update_surfaceFlx_liquidFlux_infiltration ! final computations for infiltration and runoff
   end subroutine update_surfaceFlx_liquidFlux
  
   subroutine update_surfaceFlx_liquidFlux_computation 
    integer(i4b),device :: iGRU_d
    ! **** Update operations for surfaceFlx: flux condition -- main computations ****
    associate(&
     ! input: model control
    !  deriv_desired  => in_surfaceFlx % deriv_desired  , & ! flag to indicate if derivatives are desired
     ixRichards     => decisions % f_Richards     , & ! index defining the option for Richards' equation (moisture or mixdform)
    !  nRoots         => in_surfaceFlx % nRoots         , & ! number of layers that contain roots
     ! input: state and diagnostic variables
    !  mLayerTemp          => in_surfaceFlx % mLayerTemp          , & ! temperature (K)
    !  mLayerMatricHead    => in_surfaceFlx % mLayerMatricHead    , & ! matric head in each soil layer (m)
    !  mLayerVolFracLiq    => in_surfaceFlx % mLayerVolFracLiq    , & ! volumetric liquid water content in each soil layer (-)
    !  mLayerVolFracIce    => in_surfaceFlx % mLayerVolFracIce    , & ! volumetric ice content in each soil layer (-)
     ! input: pre-computed derivatives in ...
     ! note: all of these would need to be recomputed if wanted a numerical derivative
     dTheta_dTk             => deriv_data%mLayerdTheta_dTk_m             , & ! ... volumetric liquid water content w.r.t. temperature (K-1)
     dTheta_dPsi            => deriv_data%mLayerdTheta_dPsi_m            , & ! ... the soil water characteristic w.r.t. psi (m-1)
     mLayerdPsi_dTheta      => deriv_data%mLayerdPsi_dTheta_m      , & ! ... the soil water characteristic w.r.t. theta (m)
     ! input: depth of upper-most soil layer (m)
     mLayerDepth  => prog_data%mLayerDepth  , & ! depth of upper-most soil layer (m)
    !  iLayerHeight => in_surfaceFlx % iLayerHeight , & ! height at the interface of each layer (m)
     ! input: soil parameters
     rootingDepth        => mpar_data%rootingDepth, & ! rooting depth (m)
     theta_sat           => mpar_data%theta_sat,    & ! soil porosity (-)
     ! output: error control
    !  err     => out_surfaceFlx % err    , & ! error code
    !  message => out_surfaceFlx % message,  & ! error message
     ! input: transmittance
     surfaceSatHydCond => flux_data%iLayerSatHydCond_m, & ! saturated hydraulic conductivity at the surface (m s-1)
     ! input: soil parameters
     zScale_TOPMODEL     => mpar_data%zScale_TOPMODEL     , & ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
     wettingFrontSuction => mpar_data%wettingFrontSuction , & ! Green-Ampt wetting front suction (m)
     ! input-output: surface runoff and infiltration flux (m s-1)
     xMaxInfilRate    => flux_data%scalarMaxInfilRate,  & ! maximum infiltration rate (m s-1)
     ! input: model control
    !  nSoil          => in_surfaceFlx % nSoil           , & ! number of soil layers
     ! input: soil parameters
     qSurfScale       => mpar_data%qSurfScale    , & ! scaling factor in the surface runoff parameterization (-)
     ! input-output: surface runoff and infiltration flux (m s-1)
     scalarInfilArea  => diag_data%scalarInfilArea, & ! fraction of unfrozen area where water can infiltrate (-)
     ! input: model control
    !  ixIce          => in_surfaceFlx % ixIce , & ! index of lowest ice layer
     ! input: pre-computed derivatives in ...
     ! note: all of these would need to be recomputed if wanted a numerical derivative
    !  above_soilLiqFluxDeriv => in_surfaceFlx % above_soilLiqFluxDeriv , & ! ... layer above soil (canopy or snow) liquid flux w.r.t. liquid water
    !  above_soildLiq_dTk     => in_surfaceFlx % above_soildLiq_dTk     , & ! ... layer above soil (canopy or snow) liquid flux w.r.t. temperature
    !  above_soilFracLiq      => in_surfaceFlx % above_soilFracLiq      , & ! ... liquid water layer above soil (canopy or snow) (-)
     ! input: flux at the upper boundary
     scalarRainPlusMelt => flux_data%scalarRainPlusMelt , & ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
     ! input: soil parameters
     soilIceScale        => mpar_data%soilIceScale        , & ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
     soilIceCV           => mpar_data%soilIceCV           , & ! soil ice CV in Gamma distribution used to define frozen area (-)
     ! input-output: surface runoff and infiltration flux (m s-1)
     scalarFrozenArea => diag_data%scalarFrozenArea,   & ! fraction of area that is considered impermeable due to soil ice (-)
     ! output: derivatives in surface infiltration w.r.t. ...
     dq_dHydStateVec => deriv_data%dq_dHydStateLayerSurfVec_m  , & ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
     dq_dNrgStateVec => deriv_data%dq_dNrgStateLayerSurfVec_m    & ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
    &)
    call surfaceFlx_liquidFlux_computation(nSoil,nSoil_d,nSnow,ixIce,nGRU,nRoots,ixRichards,mLayerTemp,&
    mLayerMatricHead,mLayerVolFracLiq,mLayerVolFracIce,&
    dTheta_dTk, dTheta_dPsi, mLayerdPsi_dTheta,&
    mLayerDepth,iLayerHeight,&
    rootingDepth, theta_sat, surfaceSatHydCond,zScale_TOPMODEL,wettingFrontSuction,&
    xMaxInfilRate, qSurfScale,scalarInfilArea,scalarFrozenArea, &
    above_soilLiqFluxDeriv, above_soildLiq_dTk, above_soilFracLiq, &
    scalarRainPlusMelt, soilIceScale, soilIceCV,&
    dq_dHydStateVec, dq_dNrgStateVec)
  
  
    end associate
   end subroutine update_surfaceFlx_liquidFlux_computation
  
   subroutine update_surfaceFlx_liquidFlux_infiltration
    ! **** Update operations for surfaceFlx: flux condition -- final infiltration and runoff calculations ****
    associate(&
     ! input: flux at the upper boundary
     scalarRainPlusMelt => flux_data%scalarRainPlusMelt , & ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
     ! input-output: hydraulic conductivity and diffusivity at the surface
     ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
    !  surfaceHydCond => io_surfaceFlx % surfaceHydCond , & ! hydraulic conductivity (m s-1)
    !  surfaceDiffuse => io_surfaceFlx % surfaceDiffuse , & ! hydraulic diffusivity at the surface (m
     ! input-output: surface runoff and infiltration flux (m s-1)
     xMaxInfilRate    => flux_data%scalarMaxInfilRate    , & ! maximum infiltration rate (m s-1)
     scalarInfilArea  => diag_data%scalarInfilArea  , & ! fraction of unfrozen area where water can infiltrate (-)
     scalarFrozenArea => diag_data%scalarFrozenArea , & ! fraction of area that is considered impermeable due to soil ice (-)
     ! output: runoff and infiltration
     scalarSurfaceRunoff       => flux_data%scalarSurfaceRunoff       , & ! surface runoff (m s-1)
     scalarSurfaceInfiltration => flux_data%scalarInfiltration   & ! surface infiltration (m s-1)
    &)
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
     ! compute infiltration (m s-1), if after first flux call in a splitting operation does not change
     scalarSurfaceInfiltration(iGRU) = (1._rkind - scalarFrozenArea(iGRU))*scalarInfilArea(iGRU)*min(scalarRainPlusMelt(iGRU),xMaxInfilRate(iGRU))
   
     ! compute surface runoff (m s-1)
     scalarSurfaceRunoff(iGRU) = scalarRainPlusMelt(iGRU) - scalarSurfaceInfiltration(iGRU)
   
     ! set surface hydraulic conductivity and diffusivity to missing (not used for flux condition)
     surfaceHydCond(0,iGRU) = realMissing
     surfaceDiffuse(0,iGRU) = realMissing
    end do
    end associate
  
   end subroutine update_surfaceFlx_liquidFlux_infiltration
  
   subroutine finalize_surfaceFlx
    ! **** Finalize operations for surfaceFlx ****
    ! final error check
    ! associate(&
    !  err     => out_surfaceFlx % err    , & ! error code
    !  message => out_surfaceFlx % message  & ! error message
    ! &)
     if (err /= 0_i4b) then
      message=trim(message)//'unanticipated error in surfaceFlx subroutine'; return_flag=.true.; return
     end if
    ! end associate
   end subroutine finalize_surfaceFlx
  
  end subroutine surfaceFlx
  
  attributes(device) subroutine surfaceFlx_prescribedHead(scalarSurfaceRunoff,surfaceHydCond,surfaceDiffuse,&
    scalarMatricHeadLiq, scalarVolFracLiq,&
    mLayerDepth,dq_dHydStateVec,dq_dNrgStateVec,&
    upperBoundHead,upperBoundTheta,&
    surfaceSatHydCond,dHydCond_dTemp,iceImpedeFac,&
    vGn_alpha,vGn_n,vGn_m,theta_sat,theta_res,&
    scalarSurfaceInfiltration,&
    ixRichards)
    USE soil_utils_module,only:volFracLiq            ! compute volumetric fraction of liquid water as a function of matric head (-)
    USE soil_utils_module,only:hydCond_psi           ! compute hydraulic conductivity as a function of matric head (m s-1)
    USE soil_utils_module,only:hydCond_liq           ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
    USE soil_utils_module,only:dPsi_dTheta           ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
    USE soil_utils_module,only:crit_soilT            ! compute critical temperature below which ice exists
    USE soil_utils_module,only:gammp                 ! compute the cumulative probabilty based on the Gamma distribution
  
    real(rkind),intent(inout) :: scalarSurfaceRunoff,surfaceHydCond,surfaceDiffuse
    real(rkind),intent(in) :: scalarMatricHeadLiq, scalarVolFracLiq
    real(rkind),intent(inout) :: mLayerDepth,dq_dHydStateVec,dq_dNrgStateVec
    real(rkind),intent(in) :: upperBoundHead,upperBoundTheta
    real(rkind),intent(in) :: surfaceSatHydCond,dHydCond_dTemp,iceImpedeFac
    real(rkind),intent(in) :: vGn_alpha,vGn_n,vGn_m,theta_sat,theta_res
    real(rkind),intent(inout) :: scalarSurfaceInfiltration
    integer(i4b),intent(in) :: ixRichards
  
    real(rkind)                      :: cFlux                               ! capillary flux (m s-1)
    real(rkind)                      :: dNum                                ! numerical derivative
  
       ! surface runoff iz zero for the head condition
    scalarSurfaceRunoff = 0._rkind
  
    ! compute transmission and the capillary flux
    select case(ixRichards)  ! select form of Richards' equation
      case(moisture)
        ! compute the hydraulic conductivity and diffusivity at the boundary
        surfaceHydCond = hydCond_liq(upperBoundTheta,surfaceSatHydCond,theta_res,theta_sat,vGn_m) * iceImpedeFac
        surfaceDiffuse = dPsi_dTheta(upperBoundTheta,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * surfaceHydCond
        ! compute the capillary flux
        cflux = -surfaceDiffuse*(scalarVolFracLiq - upperBoundTheta) / (mLayerDepth*0.5_rkind)
      case(mixdform)
        ! compute the hydraulic conductivity and diffusivity at the boundary
        surfaceHydCond = hydCond_psi(upperBoundHead,surfaceSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
        surfaceDiffuse = realMissing
        ! compute the capillary flux
        cflux = -surfaceHydCond*(scalarMatricHeadLiq - upperBoundHead) / (mLayerDepth*0.5_rkind)
    end select  ! end select form of Richards' eqn
  
    ! compute the total flux
    scalarSurfaceInfiltration = cflux + surfaceHydCond
  
    ! compute the derivative
      ! compute the hydrology derivative at the surface
      select case(ixRichards)  ! select form of Richards' equation
        case(moisture); dq_dHydStateVec = -surfaceDiffuse/(mLayerDepth/2._rkind)
        case(mixdform); dq_dHydStateVec = -surfaceHydCond/(mLayerDepth/2._rkind)
        ! case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
      end select
      ! compute the energy derivative at the surface
      dq_dNrgStateVec = -(dHydCond_dTemp/2._rkind)*(scalarMatricHeadLiq - upperBoundHead)/(mLayerDepth*0.5_rkind) + dHydCond_dTemp/2._rkind
  end subroutine surfaceFlx_prescribedHead
  subroutine surfaceFlx_liquidFlux_computation(nSoil,nSoil_d,nSnow,ixIce,nGRU,nRoots,ixRichards,mLayerTemp,&
    mLayerMatricHead,mLayerVolFracLiq,mLayerVolFracIce,&
    dTheta_dTk, dTheta_dPsi, mLayerdPsi_dTheta,&
    mLayerDepth,iLayerHeight,&
    rootingDepth, theta_sat, surfaceSatHydCond,zScale_TOPMODEL,wettingFrontSuction,&
    xMaxInfilRate, qSurfScale,scalarInfilArea,scalarFrozenArea, &
    above_soilLiqFluxDeriv, above_soildLiq_dTk, above_soilFracLiq, &
    scalarRainPlusMelt, soilIceScale, soilIceCV,&
    dq_dHydStateVec, dq_dNrgStateVec)
    ! USE soil_utils_module,only:volFracLiq            ! compute volumetric fraction of liquid water as a function of matric head (-)
    ! USE soil_utils_module,only:hydCond_psi           ! compute hydraulic conductivity as a function of matric head (m s-1)
    ! USE soil_utils_module,only:hydCond_liq           ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
    ! USE soil_utils_module,only:dPsi_dTheta           ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
    USE soil_utils_module,only:crit_soilT_d            ! compute critical temperature below which ice exists
    ! USE soil_utils_module,only:gammp                 ! compute the cumulative probabilty based on the Gamma distribution
  
    integer(i4b),intent(in),device :: ixIce(:), nRoots(:), ixRichards,nSnow(:)
    integer(i4b) :: iGRU,nGRU,nSoil,nSoil_d
    real(rkind),intent(in),device :: mLayerTemp(:,:)
    real(rkind),intent(in),device :: mLayerMatricHead(:,:), mLayerVolFracLiq(:,:), mLayerVolFracIce(:,:)
    real(rkind),intent(in),device :: dTheta_dTk(:,:), dTheta_dPsi(:,:), mLayerdPsi_dTheta(:,:)
    real(rkind),intent(in),device :: mLayerDepth(:,:)
    real(rkind),intent(in),device :: iLayerHeight(0:,:)
    real(rkind),intent(in),device :: rootingDepth
    real(rkind),intent(in),device :: theta_sat(:), surfaceSatHydCond(:,:),zScale_TOPMODEL,wettingFrontSuction
    real(rkind),intent(inout),device :: xMaxInfilRate(:),scalarInfilArea(:),scalarFrozenArea(:)
    real(rkind),intent(in),device :: qSurfScale
    real(rkind),intent(in),device :: above_soilLiqFluxDeriv(:), above_soildLiq_dTk(:), above_soilFracLiq(:)
    real(rkind),intent(in),device :: scalarRainPlusMelt(:), soilIceScale, soilIceCV
    real(rkind),intent(inout),device :: dq_dHydStateVec(:,:), dq_dNrgStateVec(:,:)
  
      ! -----------------------------------------------------------------------------------------------------------------------------
    ! local variables
    ! general
    real(rkind) :: satCond1, satCond2
    integer(i4b)                     :: iLayer,iLayer2                              ! index of soil layer
    real(rkind)                      :: Tcrit                               ! temperature where all water is unfrozen (K)
    real(rkind)                      :: fPart1,fPart2                       ! different parts of a function
    real(rkind)                      :: dPart1     ! derivatives for different parts of a function
    real(rkind)                      :: dPart2     ! derivatives for different parts of a function
    real(rkind)                      :: dfracCap   ! derivatives for different parts of a function
    real(rkind)                      :: dfInfRaw   ! derivatives for different parts of a function
    real(rkind)                      :: total_soil_depth                    ! total depth of soil (m)
    ! head boundary condition
    real(rkind)                      :: cFlux                               ! capillary flux (m s-1)
    real(rkind)                      :: dNum                                ! numerical derivative
    ! simplified Green-Ampt infiltration
    real(rkind)                      :: rootZoneLiq                         ! depth of liquid water in the root zone (m)
    real(rkind)                      :: rootZoneIce                         ! depth of ice in the root zone (m)
    real(rkind)                      :: availCapacity                       ! available storage capacity in the root zone (m)
    real(rkind)                      :: depthWettingFront                   ! depth to the wetting front (m)
    real(rkind)                      :: hydCondWettingFront                 ! hydraulic conductivity at the wetting front (m s-1)
    ! saturated area associated with variable storage capacity
    real(rkind)                      :: fracCap                             ! fraction of pore space filled with liquid water and ice (-)
    real(rkind)                      :: fInfRaw                             ! infiltrating area before imposing solution constraints (-)
    real(rkind),parameter            :: maxFracCap=0.995_rkind              ! maximum fraction capacity -- used to avoid numerical problems associated with an enormous derivative
    real(rkind),parameter            :: scaleFactor=0.000001_rkind          ! scale factor for the smoothing function (-)
    real(rkind),parameter            :: qSurfScaleMax=1000._rkind           ! maximum surface runoff scaling factor (-)
    ! fraction of impermeable area associated with frozen ground
    real(rkind)                      :: alpha                               ! shape parameter in the Gamma distribution
    real(rkind)                      :: xLimg                               ! upper limit of the integral
    ! derivatives in ...
    real(rkind) :: dVolFracLiq_dWat  ! ... vol fraction of liquid w.r.t. water state variable in root layers
    real(rkind) :: dVolFracIce_dWat  ! ... vol fraction of ice w.r.t. water state variable in root layers
    real(rkind) :: dVolFracLiq_dTk   ! ... vol fraction of liquid w.r.t. temperature in root layers
    real(rkind) :: dVolFracIce_dTk   ! ... vol fraction of ice w.r.t. temperature in root layers
    real(rkind) :: dRootZoneLiq_dWat ! ... vol fraction of scalar root zone liquid w.r.t. water state variable in root layers
    real(rkind) :: dRootZoneIce_dWat ! ... vol fraction of scalar root zone ice w.r.t. water state variable in root layers
    real(rkind) :: dRootZoneLiq_dTk  ! ... vol fraction of scalar root zone liquid w.r.t. temperature in root layers
    real(rkind) :: dRootZoneIce_dTk  ! ... vol fraction of scalar root zone ice w.r.t. temperature in root layers
    real(rkind) :: dDepthWettingFront_dWat ! ... scalar depth of wetting front w.r.t. water state variable in root layers
    real(rkind) :: dDepthWettingFront_dTk  ! ... scalar depth of wetting front w.r.t. temperature in root layers
    real(rkind) :: dxMaxInfilRate_dWat ! ... scalar max infiltration rate w.r.t. water state variable in root layers
    real(rkind) :: dxMaxInfilRate_dTk  ! ... scalar max infiltration rate w.r.t. temperature in root layers
    real(rkind) :: dInfilArea_dWat,dInfilArea_dWat0  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
    real(rkind) :: dInfilArea_dTk,dInfilArea_dTk0   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
    real(rkind) :: dFrozenArea_dWat,dFrozenArea_dWat0 ! ... scalar frozen area w.r.t. water state variable in canopy or snow and root layers
    real(rkind) :: dFrozenArea_dTk,dFrozenArea_dTk0  ! ... scalar frozen area w.r.t. temperature in canopy or snow and root layers
    real(rkind) :: dInfilRate_dWat,dInfilRate_dWat0  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
    real(rkind) :: dInfilRate_dTk,dInfilRate_dTk0   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  
    !$cuf kernel do(2) <<<*,*>>>
    do iGRU=1,nGRU
    do iLayer=1,nSoil
        ! process root layers only liquid and ice derivatives
      dVolFracLiq_dWat = 0._rkind
      dVolFracIce_dWat = 0._rkind
      dVolFracLiq_dTk  = 0._rkind
      dVolFracIce_dTk  = 0._rkind
      if (nRoots(iGRU) > 0) then
        select case(ixRichards)  ! form of Richards' equation
          case(moisture)
              dVolFracLiq_dWat = 1._rkind
              dVolFracIce_dWat = mLayerdPsi_dTheta(iLayer,iGRU) - 1._rkind
          case(mixdform)
            if (iLayer .le. nRoots(iGRU)) then
              Tcrit = crit_soilT_d( mLayerMatricHead(iLayer,iGRU) )
              if (mLayerTemp(iLayer+nSnow(iGRU),iGRU) < Tcrit) then
                dVolFracLiq_dWat = 0._rkind
                dVolFracIce_dWat = dTheta_dPsi(iLayer,iGRU)
              else
                dVolFracLiq_dWat = dTheta_dPsi(iLayer,iGRU)
                dVolFracIce_dWat = 0._rkind
              end if
            end if
        end select 
          dVolFracLiq_dTk = dTheta_dTk(iLayer+nSnow(iGRU),iGRU) !already zeroed out if not below critical temperature
          dVolFracIce_dTk = -dVolFracLiq_dTk !often can and will simplify one of these terms out
      end if
  
      ! define the storage in the root zone (m) and derivatives
      rootZoneLiq = 0._rkind
      rootZoneIce = 0._rkind
      dRootZoneLiq_dWat = 0._rkind
      dRootZoneIce_dWat = 0._rkind
      dRootZoneLiq_dTk  = 0._rkind
      dRootZoneIce_dTk  = 0._rkind
      ! process layers where the roots extend to the bottom of the layer
      if (nRoots(iGRU) > 1) then
        do iLayer2=1,nRoots(iGRU)-1
          rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(iLayer2+nSnow(iGRU),iGRU)*mLayerDepth(iLayer2+nSnow(iGRU),iGRU)
          rootZoneIce = rootZoneIce + mLayerVolFracIce(iLayer2+nSnow(iGRU),iGRU)*mLayerDepth(iLayer2+nSnow(iGRU),iGRU)
        end do
        if (iLayer .le. nRoots(iGRU) - 1) then
          dRootZoneLiq_dWat = dVolFracLiq_dWat*mLayerDepth(iLayer+nSnow(iGRU),iGRU)
          dRootZoneIce_dWat = dVolFracIce_dWat*mLayerDepth(iLayer+nSnow(iGRU),iGRU)
          dRootZoneLiq_dTk  = dVolFracLiq_dTk *mLayerDepth(iLayer+nSnow(iGRU),iGRU)
          dRootZoneIce_dTk  = dVolFracIce_dTk *mLayerDepth(iLayer+nSnow(iGRU),iGRU)
        end if
      end if
  
      ! process layers where the roots end in the current layer
      rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(nRoots(iGRU)+nSnow(iGRU),iGRU)*(rootingDepth - iLayerHeight(nRoots(iGRU)-1,iGRU))
      rootZoneIce = rootZoneIce + mLayerVolFracIce(nRoots(iGRU)+nSnow(iGRU),iGRU)*(rootingDepth - iLayerHeight(nRoots(iGRU)-1,iGRU))
      if (iLayer .eq. nRoots(iGRU)) then
        dRootZoneLiq_dWat = dVolFracLiq_dWat*(rootingDepth - iLayerHeight(nRoots(iGRU)-1,iGRU))
        dRootZoneIce_dWat = dVolFracIce_dWat*(rootingDepth - iLayerHeight(nRoots(iGRU)-1,iGRU))
        dRootZoneLiq_dTk  = dVolFracLiq_dTk* (rootingDepth - iLayerHeight(nRoots(iGRU)-1,iGRU))
        dRootZoneIce_dTk  = dVolFracIce_dTk* (rootingDepth - iLayerHeight(nRoots(iGRU)-1,iGRU))
      end if
      availCapacity = theta_sat(1)*rootingDepth - rootZoneIce
      ! if (rootZoneLiq > availCapacity+verySmall) then
      !   ! message=trim(message)//'liquid water in the root zone exceeds capacity'
      !   ! return
      ! end if
       ! define the depth to the wetting front (m) and derivatives
      total_soil_depth = 0
      ! not cure why the standard loop setup doesn't work here, but it seems to cause seg faults
      iLayer2 = 1
      do
        if (iLayer2 .gt. nSoil_d) exit
        total_soil_depth = total_soil_depth+mLayerDepth(iLayer2+nSnow(iGRU),iGRU)
        iLayer2 = iLayer2 + 1
      end do
      ! do iLayer2=1,nSoil_d
      !   ! total_soil_depth = total_soil_depth+mLayerDepth(iLayer2+nSnow,iGRU)
      ! end do
      depthWettingFront = (rootZoneLiq/availCapacity)*min(rootingDepth, total_soil_depth)
      dDepthWettingFront_dWat=( dRootZoneLiq_dWat*min(rootingDepth, total_soil_depth) + dRootZoneIce_dWat*depthWettingFront )/availCapacity
      dDepthWettingFront_dTk =( dRootZoneLiq_dTk *min(rootingDepth, total_soil_depth) + dRootZoneIce_dTk*depthWettingFront  )/availCapacity
  
      ! define the hydraulic conductivity at depth=depthWettingFront (m s-1)
      hydCondWettingFront =  surfaceSatHydCond(0,iGRU) * ( (1._rkind - depthWettingFront/total_soil_depth)**(zScale_TOPMODEL - 1._rkind) )
  
      ! define the maximum infiltration rate (m s-1) and derivatives
      xMaxInfilRate(iGRU) = hydCondWettingFront*( (wettingFrontSuction + depthWettingFront)/depthWettingFront )  ! maximum infiltration rate (m s-1)
  
      fPart1    = hydCondWettingFront
      fPart2    = (wettingFrontSuction + depthWettingFront)/depthWettingFront
      dPart1 = surfaceSatHydCond(0,iGRU)*(zScale_TOPMODEL - 1._rkind) * ( (1._rkind - depthWettingFront/total_soil_depth)**(zScale_TOPMODEL - 2._rkind) ) * (-dDepthWettingFront_dWat)/total_soil_depth
      dPart2 = -dDepthWettingFront_dWat*wettingFrontSuction / (depthWettingFront**2_i4b)
      dxMaxInfilRate_dWat = fPart1*dPart2 + fPart2*dPart1
      dPart1 = surfaceSatHydCond(0,iGRU)*(zScale_TOPMODEL - 1._rkind) * ( (1._rkind - depthWettingFront/total_soil_depth)**(zScale_TOPMODEL - 2._rkind) ) * (-dDepthWettingFront_dTk)/total_soil_depth
      dPart2 = -dDepthWettingFront_dTk*wettingFrontSuction / (depthWettingFront**2_i4b)
      dxMaxInfilRate_dTk  = fPart1*dPart2 + fPart2*dPart1
  
  
      ! define the infiltrating area and derivatives for the non-frozen part of the cell/basin
      if (qSurfScale < qSurfScaleMax) then
        fracCap         = rootZoneLiq/(maxFracCap*availCapacity)                              ! fraction of available root zone filled with water
        fInfRaw         = 1._rkind - exp(-qSurfScale*(1._rkind - fracCap))                          ! infiltrating area -- allowed to violate solution constraints
        scalarInfilArea(iGRU) = min(0.5_rkind*(fInfRaw + sqrt(fInfRaw**2_i4b + scaleFactor)), 1._rkind)   ! infiltrating area -- constrained
        if (0.5_rkind*(fInfRaw + sqrt(fInfRaw**2_i4b + scaleFactor))< 1._rkind) then
          dfracCap = ( dRootZoneLiq_dWat/maxFracCap + dRootZoneIce_dWat*fracCap )/availCapacity
          dfInfRaw = -qSurfScale*dfracCap * exp(-qSurfScale*(1._rkind - fracCap))
          dInfilArea_dWat = 0.5_rkind*dfInfRaw * (1._rkind + fInfRaw/sqrt(fInfRaw**2_i4b + scaleFactor))
          dfracCap = ( dRootZoneLiq_dTk/maxFracCap + dRootZoneIce_dTk*fracCap )/availCapacity
          dfInfRaw = -qSurfScale*dfracCap * exp(-qSurfScale*(1._rkind - fracCap))
          dInfilArea_dTk  = 0.5_rkind*dfInfRaw * (1._rkind + fInfRaw/sqrt(fInfRaw**2_i4b + scaleFactor))
        else ! scalarInfilArea = 1._rkind
          dInfilArea_dWat = 0._rkind
          dInfilArea_dTk  = 0._rkind
        end if
      else
        scalarInfilArea(iGRU) = 1._rkind
        dInfilArea_dWat = 0._rkind
        dInfilArea_dTk  = 0._rkind
      end if
  
       ! check to ensure we are not infiltrating into a fully saturated column
      if (ixIce(iGRU)<nRoots(iGRU)) then
        satCond1 = 0
        satCond2 = 0
        do iLayer2=ixIce(iGRU)+1,nRoots(iGRU)
          satCond1 = satCond1 + mLayerVolFracLiq(iLayer2+nSnow(iGRU),iGRU)*mLayerDepth(iLayer2+nSnow(iGRU),iGRU)
          satCond2 = satCond2 + mLayerDepth(iLayer2+nSnow(iGRU),iGRU)
        end do
  
        if (satCond1 > 0.9999_rkind*theta_sat(1)*satCond2) scalarInfilArea(iGRU)=0._rkind
      end if
  
          ! define the impermeable area and derivatives due to frozen ground
      if (rootZoneIce > tiny(rootZoneIce)) then  ! (avoid divide by zero)
        alpha            = 1._rkind/(soilIceCV**2_i4b)        ! shape parameter in the Gamma distribution
        xLimg            = alpha*soilIceScale/rootZoneIce  ! upper limit of the integral
  
        !if we use this, we will have a derivative of scalarFrozenArea w.r.t. water and temperature in each layer (through mLayerVolFracIce)
        scalarFrozenArea(iGRU) = 0._rkind
        dFrozenArea_dWat = 0._rkind
        dFrozenArea_dTk  = 0._rkind
      else
        scalarFrozenArea(iGRU) = 0._rkind
        dFrozenArea_dWat = 0._rkind
        dFrozenArea_dTk  = 0._rkind
      end if
      dInfilArea_dWat0 = 0._rkind
      dInfilArea_dTk0  = 0._rkind
      dFrozenArea_dWat0 = 0._rkind
      dFrozenArea_dTk0  = 0._rkind
  
      if (xMaxInfilRate(iGRU) < scalarRainPlusMelt(iGRU)) then ! = dxMaxInfilRate_d, dependent on layers not at surface
        dInfilRate_dWat0 = 0._rkind
        dInfilRate_dTk0  = 0._rkind
        dInfilRate_dWat = dxMaxInfilRate_dWat
        dInfilRate_dTk  = dxMaxInfilRate_dTk
      else ! = dRainPlusMelt_d, dependent on above layer (canopy or snow) water and temp
        dInfilRate_dWat0 = above_soilLiqFluxDeriv(iGRU)*above_soilFracLiq(iGRU)
        dInfilRate_dTk0  = above_soilLiqFluxDeriv(iGRU)*above_soildLiq_dTk(iGRU)
        dInfilRate_dWat = 0._rkind
        dInfilRate_dTk  = 0._rkind
      end if
    dq_dHydStateVec(0,iGRU) = (1._rkind - scalarFrozenArea(iGRU))&
    & * ( scalarInfilArea(iGRU)*dInfilRate_dWat0 )
  dq_dNrgStateVec(0,iGRU) = (1._rkind - scalarFrozenArea(iGRU))&
    & * ( scalarInfilArea(iGRU)*dInfilRate_dTk0  )
  
       ! dq w.r.t. infiltration only, scalarRainPlusMelt accounted for in computJacob module
    dq_dHydStateVec(iLayer,iGRU) = (1._rkind - scalarFrozenArea(iGRU))&
    & * ( dInfilArea_dWat*min(scalarRainPlusMelt(iGRU),xMaxInfilRate(iGRU)) + scalarInfilArea(iGRU)*dInfilRate_dWat )&
    & + (-dFrozenArea_dWat)*scalarInfilArea(iGRU)*min(scalarRainPlusMelt(iGRU),xMaxInfilRate(iGRU))
  dq_dNrgStateVec(iLayer,iGRU) = (1._rkind - scalarFrozenArea(iGRU))&
    & * ( dInfilArea_dTk *min(scalarRainPlusMelt(iGRU),xMaxInfilRate(iGRU)) + scalarInfilArea(iGRU)*dInfilRate_dTk  )&
    & + (-dFrozenArea_dTk) *scalarInfilArea(iGRU)*min(scalarRainPlusMelt(iGRU),xMaxInfilRate(iGRU))
  
    end do
  end do
  
  end subroutine
  ! ***************************************************************************************************************
  ! private subroutine iLayerFlux: compute the fluxes and derivatives at layer interfaces
  ! ***************************************************************************************************************
  subroutine iLayerFlux(ixTop,ixBot,nSnow,nGRU,decisions,prog_data,flux_data,deriv_data,&
    nodeDiffuseTrial,dHydCond_dTemp,dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dMatric,&
                                      nodeVolFracLiqTrial, nodeMatricHeadLiqTrial, iLayerHydCond, iLayerDiffuse)
    ! ---------------------------------------------------------------------------------------------------------------------------
    ! input: model control, state variables, coordinate variables, temperature derivatives, transmittance variables
                                      use device_data_types
                                      use data_types
    ! type(in_type_iLayerFlux),intent(in)   :: in_iLayerFlux   ! class object for input data
    integer(i4b),intent(in) :: ixTop,ixBot,nGRU
    integer(i4b),intent(in),device :: nSnow(:)
    ! output: transmittance variables and vertical flux at layer interface, derivatives, and error control
    type(decisions_device),intent(in) :: decisions
    type(prog_data_device),intent(in) :: prog_data
    type(flux_data_device),intent(inout) :: flux_data
    type(deriv_data_device),intent(inout) :: deriv_data
    ! type(out_type_iLayerFlux),intent(out) :: out_iLayerFlux  ! class object for output data
    real(rkind),intent(in),device                :: nodeDiffuseTrial(:,:)  ! diffusivity at layer mid-point (m2 s-1)
    real(rkind),intent(in),device :: dHydCond_dTemp(:,:)   ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
    real(rkind),intent(in),device :: dHydCond_dVolLiq(:,:) ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
    real(rkind),intent(in),device :: dDiffuse_dVolLiq(:,:) ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
    real(rkind),intent(in),device :: dHydCond_dMatric(:,:)
    real(rkind),intent(in),device :: nodeVolFracLiqTrial(:,:), nodeMatricHeadLiqTrial(:,:)
    real(rkind),intent(inout),device :: iLayerHydCond(0:,:), iLayerDiffuse(0:,:)
  
    ! ---------------------------------------------------------------------------------------------------------------------------
    ! local variables (named variables to provide index of 2-element vectors)
    integer(i4b),parameter           :: ixUpper=1            ! index of upper node in the 2-element vectors
    integer(i4b),parameter           :: ixLower=2            ! index of lower node in the 2-element vectors
    integer(i4b) :: iLayer,iGRU
    logical(lgt),parameter           :: useGeometric=.false. ! switch between the arithmetic and geometric mean
    ! local variables (Darcy flux)
    real(rkind)                      :: dPsi                 ! spatial difference in matric head (m)
    real(rkind)                      :: dLiq                 ! spatial difference in volumetric liquid water (-)
    real(rkind)                      :: dz                   ! spatial difference in layer mid-points (m)
    real(rkind)                      :: cflux                ! capillary flux (m s-1)
    ! error control
    logical(lgt)                     :: return_flag          ! flag for return statements
    ! ---------------------------------------------------------------------------------------------------------------------------
  
    ! **** Update operations for iLayerFlux: compute fluxes ****
    associate(&
     ! input: model control
     ixRichards    => decisions%f_Richards   , & ! index defining the option for Richards' equation (moisture or mixdform)
     ! input: state variables
    !  nodeMatricHeadLiqTrial => in_iLayerFlux % nodeMatricHeadLiqTrial, & ! liquid matric head at the soil nodes (m)
    !  nodeVolFracLiqTrial    => in_iLayerFlux % nodeVolFracLiqTrial   , & ! volumetric fraction of liquid water at the soil nodes (-)
     ! input: model coordinate variables
     nodeHeight => prog_data%mLayerHeight, & ! height at the mid-point of the lower layer (m)
     ! input: transmittance
     nodeHydCondTrial => flux_data%mLayerHydCond_m, & ! hydraulic conductivity at layer mid-points (m s-1)
    !  nodeDiffuseTrial => in_iLayerFlux % nodeDiffuseTrial, & ! diffusivity at layer mid-points (m2 s-1)
     ! output: tranmsmittance at the layer interface (scalars)
    !  iLayerHydCond => out_iLayerFlux % iLayerHydCond, & ! hydraulic conductivity at the interface between layers (m s-1)
    !  iLayerDiffuse => out_iLayerFlux % iLayerDiffuse, & ! hydraulic diffusivity at the interface between layers (m2 s-1)
     ! output: vertical flux at the layer interface (scalars)
     iLayerLiqFluxSoil => flux_data%iLayerLiqFluxSoil_m, & ! vertical flux of liquid water at the layer interface (m s-1)
     ! output: error control
     ! input: temperature derivatives
     dPsiLiq_dTemp   => deriv_data%dPsiLiq_dTemp_m , & ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
    !  dHydCond_dTemp  => in_iLayerFlux % dHydCond_dTemp, & ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
     ! input: transmittance derivatives
    !  dHydCond_dVolLiq => in_iLayerFlux % dHydCond_dVolLiq, & ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
    !  dDiffuse_dVolLiq => in_iLayerFlux % dDiffuse_dVolLiq, & ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
    !  dHydCond_dMatric => in_iLayerFlux % dHydCond_dMatric, & ! derivative in hydraulic conductivity w.r.t matric head (m s-1)
     ! output: derivatives in fluxes w.r.t. ...  
     dq_dHydStateAbove => deriv_data%dq_dHydStateAbove_m, & ! ... matric head or volumetric lquid water in the layer above (m s-1 or s-1)
     dq_dHydStateBelow => deriv_data%dq_dHydStateBelow_m, & ! ... matric head or volumetric lquid water in the layer below (m s-1 or s-1)
     ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
     dq_dNrgStateAbove => deriv_data%dq_dNrgStateAbove_m, & ! derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
     dq_dNrgStateBelow => deriv_data%dq_dNrgStateBelow_m & ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
     ! output: error control
    !  err     => out_iLayerFlux % err    , & ! error code
    !  message => out_iLayerFlux % message  & ! error message
    &)
  
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
    do iLayer=ixTop,ixBot
    call iLayerFlux_flux_deriv(ixRichards,&
    nodeHydCondTrial(iLayer+1,iGRU), nodeHydCondTrial(iLayer,iGRU), &
    nodeHeight(iLayer+1+nSnow(iGRU),iGRU), nodeHeight(iLayer+nSnow(iGRU),iGRU), &
    nodeDiffuseTrial(iLayer+1,iGRU), nodeDiffuseTrial(iLayer,iGRU), &
    nodeVolFracLiqTrial(iLayer+1+nSnow(iGRU),iGRU), nodeVolFracLiqTrial(iLayer+nSnow(iGRU),iGRU), &
    nodeMatricHeadLiqTrial(iLayer+1,iGRU), nodeMatricHeadLiqTrial(iLayer,iGRU), &
    dPsiLiq_dTemp(iLayer+1,iGRU), dPsiLiq_dTemp(iLayer,iGRU), &
    dHydCond_dTemp(iLayer+1,iGRU), dHydCond_dTemp(iLayer,iGRU), &
    dHydCond_dVolLiq(iLayer+1,iGRU), dHydCond_dVolLiq(iLayer,iGRU), &
    dDiffuse_dVolLiq(iLayer+1,iGRU), dDiffuse_dVolLiq(iLayer,iGRU), &
    dHydCond_dMatric(iLayer+1,iGRU), dHydCond_dMatric(iLayer,iGRU), &
    iLayerHydCond(iLayer,iGRU), iLayerDiffuse(iLayer,iGRU), iLayerLiqFluxSoil(iLayer,iGRU), &
    dq_dHydStateAbove(iLayer,iGRU), dq_dHydStateBelow(iLayer,iGRU), dq_dNrgStateAbove(iLayer,iGRU), dq_dNrgStateBelow(iLayer,iGRU))
    end do
  end do
    end associate
   
  end subroutine iLayerFlux
  
  attributes(device) subroutine iLayerFlux_flux_deriv(ixRichards,&
    nodeHydCondTrial_Lower, nodeHydCondTrial_Upper, &
    nodeHeight_Lower, nodeHeight_Upper, &
    nodeDiffuseTrial_Lower, nodeDiffuseTrial_Upper, &
    nodeVolFracLiqTrial_Lower, nodeVolFracLiqTrial_Upper, &
    nodeMatricHeadLiqTrial_Lower, nodeMatricHeadLiqTrial_Upper, &
    dPsiLiq_dTemp_Lower, dPsiLiq_dTemp_Upper, &
    dHydCond_dTemp_Lower, dHydCond_dTemp_Upper, &
    dHydCond_dVolLiq_Lower, dHydCond_dVolLiq_Upper, &
    dDiffuse_dVolLiq_Lower, dDiffuse_dVolLiq_Upper, &
    dHydCond_dMatric_Lower, dHydCond_dMatric_Upper, &
    iLayerHydCond, iLayerDiffuse, iLayerLiqFluxSoil, &
    dq_dHydStateAbove, dq_dHydStateBelow, dq_dNrgStateAbove, dq_dNrgStateBelow)
    integer(i4b),intent(in) :: ixRichards
    real(rkind),intent(in) :: nodeHydCondTrial_Lower, nodeHydCondTrial_Upper
    real(rkind),intent(in) :: nodeHeight_Lower, nodeHeight_Upper
    real(rkind),intent(in) :: nodeDiffuseTrial_Lower, nodeDiffuseTrial_Upper
    real(rkind),intent(in) :: nodeVolFracLiqTrial_Lower, nodeVolFracLiqTrial_Upper
    real(rkind),intent(in) :: nodeMatricHeadLiqTrial_Lower, nodeMatricHeadLiqTrial_Upper
    real(rkind),intent(in) :: dPsiLiq_dTemp_Lower, dPsiLiq_dTemp_Upper
    real(rkind),intent(in) :: dHydCond_dTemp_Lower, dHydCond_dTemp_Upper
    real(rkind),intent(in) :: dHydCond_dVolLiq_Lower, dHydCond_dVolLiq_Upper
    real(rkind),intent(in) :: dDiffuse_dVolLiq_Lower, dDiffuse_dVolLiq_Upper
    real(rkind),intent(in) :: dHydCond_dMatric_Lower, dHydCond_dMatric_Upper
  
    real(rkind),intent(inout) :: iLayerHydCond, iLayerDiffuse, iLayerLiqFluxSoil
    real(rkind),intent(inout) :: dq_dHydStateAbove, dq_dHydStateBelow, dq_dNrgStateAbove, dq_dNrgStateBelow
    logical(lgt),parameter           :: useGeometric=.false. ! switch between the arithmetic and geometric mean
  
    real(rkind) :: dz, dPsi, dLiq, cflux
     ! * local variables (derivative in Darcy's flux) *
    ! deriviatives at the layer interface
    real(rkind) :: dHydCondIface_dVolLiqAbove  ! hydraulic conductivity w.r.t. volumetric liquid water content in layer above
    real(rkind) :: dHydCondIface_dVolLiqBelow  ! hydraulic conductivity w.r.t. volumetric liquid water content in layer below
    real(rkind) :: dDiffuseIface_dVolLiqAbove  ! hydraulic diffusivity  w.r.t. volumetric liquid water content in layer above
    real(rkind) :: dDiffuseIface_dVolLiqBelow  ! hydraulic diffusivity  w.r.t. volumetric liquid water content in layer below
    real(rkind) :: dHydCondIface_dMatricAbove  ! hydraulic conductivity w.r.t. matric head in layer above
    real(rkind) :: dHydCondIface_dMatricBelow  ! hydraulic conductivity w.r.t. matric head in layer below
   
     ! compute the vertical flux of liquid water
     ! compute the hydraulic conductivity at the interface
     if (useGeometric) then
       iLayerHydCond   = sqrt(nodeHydCondTrial_Lower   * nodeHydCondTrial_Upper)
     else
       iLayerHydCond   = (nodeHydCondTrial_Lower   + nodeHydCondTrial_Upper)*0.5_rkind
     end if
     
     dz = nodeHeight_Lower - nodeHeight_Upper
     ! compute the capillary flux
     select case(ixRichards)  ! select form of Richards' equation
       case(moisture)
        iLayerDiffuse = sqrt(nodeDiffuseTrial_Lower * nodeDiffuseTrial_Upper)
        dLiq          = nodeVolFracLiqTrial_Lower - nodeVolFracLiqTrial_Upper
        cflux         = -iLayerDiffuse * dLiq/dz
       case(mixdform)
        iLayerDiffuse = realMissing
        dPsi          = nodeMatricHeadLiqTrial_Lower - nodeMatricHeadLiqTrial_Upper
        cflux         = -iLayerHydCond * dPsi/dz
       case default; 
        !err=10; message=trim(message)//"unable to identify option for Richards' equation"; 
        ! return_flag=.true.; 
        return
     end select
     ! compute the total flux (add gravity flux, positive downwards)
     iLayerLiqFluxSoil = cflux + iLayerHydCond
  
     select case(ixRichards)  ! select form of Richards' equation
       case(moisture)
         ! still need to implement arithmetric mean for the moisture-based form
         if (.not.useGeometric) then
          !  message=trim(message)//'only currently implemented for geometric mean -- change local flag'
          !  err=20; return_flag=.true.; 
          return
         end if
         ! derivatives in hydraulic conductivity at the layer interface (m s-1)
         dHydCondIface_dVolLiqAbove = dHydCond_dVolLiq_Upper*nodeHydCondTrial_Lower * 0.5_rkind/max(iLayerHydCond,verySmall)
         dHydCondIface_dVolLiqBelow = dHydCond_dVolLiq_Lower*nodeHydCondTrial_Upper * 0.5_rkind/max(iLayerHydCond,verySmall)
         ! derivatives in hydraulic diffusivity at the layer interface (m2 s-1)
         dDiffuseIface_dVolLiqAbove = dDiffuse_dVolLiq_Upper*nodeDiffuseTrial_Lower * 0.5_rkind/max(iLayerDiffuse,verySmall)
         dDiffuseIface_dVolLiqBelow = dDiffuse_dVolLiq_Lower*nodeDiffuseTrial_Upper * 0.5_rkind/max(iLayerDiffuse,verySmall)
         ! derivatives in the flux w.r.t. volumetric liquid water content
         dq_dHydStateAbove = -dDiffuseIface_dVolLiqAbove*dLiq/dz + iLayerDiffuse/dz + dHydCondIface_dVolLiqAbove
         dq_dHydStateBelow = -dDiffuseIface_dVolLiqBelow*dLiq/dz - iLayerDiffuse/dz + dHydCondIface_dVolLiqBelow
       case(mixdform)
         ! derivatives in hydraulic conductivity
         if (useGeometric) then
           dHydCondIface_dMatricAbove = dHydCond_dMatric_Upper*nodeHydCondTrial_Lower * 0.5_rkind/max(iLayerHydCond,verySmall)
           dHydCondIface_dMatricBelow = dHydCond_dMatric_Lower*nodeHydCondTrial_Upper * 0.5_rkind/max(iLayerHydCond,verySmall)
         else
           dHydCondIface_dMatricAbove = dHydCond_dMatric_Upper/2._rkind
           dHydCondIface_dMatricBelow = dHydCond_dMatric_Lower/2._rkind
         end if
         ! derivatives in the flux w.r.t. matric head
         dq_dHydStateAbove = -dHydCondIface_dMatricAbove*dPsi/dz + iLayerHydCond/dz + dHydCondIface_dMatricAbove
         dq_dHydStateBelow = -dHydCondIface_dMatricBelow*dPsi/dz - iLayerHydCond/dz + dHydCondIface_dMatricBelow
         ! derivative in the flux w.r.t. temperature
         dq_dNrgStateAbove = -(dHydCond_dTemp_Upper/2._rkind)*dPsi/dz + iLayerHydCond*dPsiLiq_dTemp_Upper/dz + dHydCond_dTemp_Upper/2._rkind
         dq_dNrgStateBelow = -(dHydCond_dTemp_Lower/2._rkind)*dPsi/dz - iLayerHydCond*dPsiLiq_dTemp_Lower/dz + dHydCond_dTemp_Lower/2._rkind
       case default; 
        ! err=10; message=trim(message)//"unknown form of Richards' equation"; 
        ! return_flag=.true.; 
        return
     end select
  
   end subroutine iLayerFlux_flux_deriv
  
  ! ***************************************************************************************************************
  ! private subroutine qDrainFlux: compute the drainage flux from the bottom of the soil profile and its derivative
  ! ***************************************************************************************************************
  subroutine qDrainFlux(bc_lower,nSoil,nSnow,nGRU,decisions,mpar_data,prog_data,diag_data,flux_data,deriv_data,err,message, bottomHydCond, bottomDiffuse, &
    iceImpedeFac, dHydCond_dVolLiq, dHydCond_dTemp, dHydCond_dMatric, nodeMatricHeadLiq, nodeVolFracLiq)
    USE soil_utils_module,only:volFracLiq  ! compute volumetric fraction of liquid water as a function of matric head (-)
    USE soil_utils_module,only:matricHead  ! compute matric head as a function of volumetric fraction of liquid water (m)
    USE soil_utils_module,only:hydCond_psi ! compute hydraulic conductivity as a function of matric head (m s-1)
    USE soil_utils_module,only:hydCond_liq ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
    USE soil_utils_module,only:dPsi_dTheta ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
    ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
    use device_data_types
    implicit none
    ! -----------------------------------------------------------------------------------------------------------------------------
    ! input: model control, variables, boundary conditions, transmittance variables, and soil parameters
    integer(i4b) :: bc_lower      
    integer(i4b),intent(in) :: nSoil
    integer(i4b),intent(in),device :: nSnow(:)
    integer(i4b),intent(in) :: nGRU
    type(decisions_device),intent(in) :: decisions
    type(mpar_data_device),intent(inout) :: mpar_data
    type(prog_data_device),intent(inout) :: prog_data
    type(diag_data_device),intent(inout) :: diag_data
    type(flux_data_device),intent(inout) :: flux_data
    type(deriv_data_device),intent(inout) :: deriv_data
    ! output: hydraulic conductivity and diffusivity, drainage fluxes and derivatives, and error control
    integer(i4b) :: err
    character(*),intent(out) :: message
    real(rkind),intent(inout),device :: bottomHydCond(0:,:), bottomDiffuse(0:,:)
    real(rkind),intent(in),device :: iceImpedeFac(:,:)     ! ice impedence factor at layer mid-points (-)
    real(rkind),intent(in),device :: dHydCond_dVolLiq(:,:) ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
    real(rkind),intent(in),device :: dHydCond_dTemp(:,:)   ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
    real(rkind),intent(in),device :: dHydCond_dMatric(:,:)
    real(rkind),intent(in),device :: nodeMatricHeadLiq(:,:), nodeVolFracLiq(:,:)
  
    ! -----------------------------------------------------------------------------------------------------------------------------
    ! local variables
    integer(i4b) :: iGRU,iLayer
    ! real(rkind)                      :: zWater                  ! effective water table depth (m)
    ! real(rkind)                      :: nodePsi                 ! matric head in the lowest unsaturated node (m)
    ! real(rkind)                      :: cflux                   ! capillary flux (m s-1)
    ! error control
    logical(lgt)                     :: return_flag             ! flag for return statements
    ! -----------------------------------------------------------------------------------------------------------------------------
  
     call initialize_qDrainFlux
  
     call update_qDrainFlux;   if (return_flag) return
  
     call finalize_qDrainFlux; if (return_flag) return
  
  contains
  
   subroutine initialize_qDrainFlux
    ! ** Initialize operations for qDrainFlux **
    return_flag=.false. ! initialize return flag
    ! associate(&
     ! output: error control
    !  err     => out_qDrainFlux % err    , & ! error code
    !  message => out_qDrainFlux % message  & ! error message
    ! &)
     ! initialize error control
     err=0; message="qDrainFlux/"
    ! end associate
   end subroutine initialize_qDrainFlux
  
   subroutine update_qDrainFlux
    ! ** Update operations for qDrainFlux **
    ! associate(&
     ! input: model control
    !  bc_lower      => in_qDrainFlux % bc_lower, & ! index defining the type of boundary conditions
     ! output: error control
    !  err     => out_qDrainFlux % err    , &       ! error code
    !  message => out_qDrainFlux % message  &       ! error message
    ! &)
  
     ! determine lower boundary condition
     select case(bc_lower)
       case(prescribedHead) ! specified matric head value
         call update_qDrainFlux_prescribedHead; if (return_flag) return
       case(funcBottomHead) ! specified matric head function
         call update_qDrainFlux_funcBottomHead; if (return_flag) return
       case(freeDrainage)   ! free drainage 
         call update_qDrainFlux_freeDrainage;   if (return_flag) return
       case(zeroFlux)       ! zero flux
         call update_qDrainFlux_zeroFlux;       if (return_flag) return
       case default;
          err=20; message=trim(message)//'unknown lower boundary condition for soil hydrology'; return_flag=.true.; return
     end select 
  
    ! end associate
   end subroutine update_qDrainFlux
  
   subroutine update_qDrainFlux_prescribedHead
    ! ** Update operations for qDrainFlux: prescribed pressure head value at bottom boundary **
    associate(&
     ! input: model control
    !  deriv_desired => in_qDrainFlux % deriv_desired, &          ! flag to indicate if derivatives are desired
     ixRichards    => decisions % f_Richards   , &          ! index defining the option for Richards' equation (moisture or mixdform)
     ! input: state and diagnostic variables
    !  nodeMatricHeadLiq => in_qDrainFlux % nodeMatricHeadLiq, &  ! liquid matric head in the lowest unsaturated node (m)
    !  nodeVolFracLiq    => in_qDrainFlux % nodeVolFracLiq   , &  ! volumetric liquid water content in the lowest unsaturated node (-)
     ! input: model coordinate variables
     nodeDepth  => prog_data%mLayerDepth , &                ! depth of the lowest unsaturated soil layer (m)
     ! input: diriclet boundary conditions
     lowerBoundHead  => mpar_data%lowerBoundHead , &      ! lower boundary condition for matric head (m)
     lowerBoundTheta => mpar_data%lowerBoundTheta, &      ! lower boundary condition for volumetric liquid water content (-)
     ! input: transmittance
     bottomSatHydCond  => flux_data%iLayerSatHydCond_m, &  ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
    !  iceImpedeFac      => in_qDrainFlux % iceImpedeFac     , &  ! ice impedence factor in the upper-most soil layer (-)
     ! input: transmittance derivatives
    !  dHydCond_dTemp   => in_qDrainFlux % dHydCond_dTemp  , &    ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
     ! input: soil parameters
     vGn_alpha       => mpar_data%vGn_alpha      , &      ! van Genuchten "alpha" parameter (m-1)
     vGn_n           => mpar_data%vGn_n          , &      ! van Genuchten "n" parameter (-)
     vGn_m           => diag_data%scalarVGn_m_m          , &      ! van Genuchten "m" parameter (-)
     theta_sat       => mpar_data%theta_sat      , &      ! soil porosity (-)
     theta_res       => mpar_data%theta_res      , &      ! soil residual volumetric water content (-)
     ! output: hydraulic conductivity at the bottom of the unsaturated zone
    !  bottomHydCond => out_qDrainFlux % bottomHydCond, &         ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
    !  bottomDiffuse => out_qDrainFlux % bottomDiffuse, &         ! hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
     ! output: drainage flux from the bottom of the soil profile
     scalarDrainage => flux_data%iLayerLiqFluxSoil_m, &       ! drainage flux from the bottom of the soil profile (m s-1)
     ! output: derivatives in drainage flux w.r.t. ...
     dq_dHydStateUnsat => deriv_data%dq_dHydStateAbove_m, & ! ... state variable in lowest unsaturated node (m s-1 or s-1)
     dq_dNrgStateUnsat => deriv_data%dq_dNrgStateAbove_m & ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
     ! output: error control
    !  err     => out_qDrainFlux % err    , &                     ! error code
    !  message => out_qDrainFlux % message  &                     ! error message
    &)
  
    !$cuf kernel do(2) <<<*,*>>>
    do iGRU=1,nGRU
      do iLayer=nSoil,nSoil
      call qDrainFlux_prescribedHead(bottomHydCond(iLayer,iGRU), bottomDiffuse(iLayer,iGRU), &
      scalarDrainage(iLayer,iGRU), dq_dHydStateUnsat(iLayer,iGRU), dq_dNrgStateUnsat(iLayer,iGRU),&
      nodeMatricHeadLiq(iLayer,iGRU), nodeVolFracLiq(iLayer+nSnow(iGRU),iGRU), nodeDepth(iLayer+nSnow(iGRU),iGRU), &
      lowerBoundHead, lowerBoundTheta, &
      bottomSatHydCond(iLayer,iGRU), iceImpedeFac(iLayer,iGRU), &
      dHydCond_dTemp(iLayer,iGRU), &
      vGn_alpha(iLayer), vGn_n(iLayer), vGn_m(iLayer,iGRU), theta_sat(iLayer), theta_res(iLayer), &
      ixRichards)
    end do
  end do
   
    end associate
   end subroutine update_qDrainFlux_prescribedHead
  
   subroutine update_qDrainFlux_funcBottomHead
    ! ** Update operations for qDrainFlux: prescribed pressure head function at bottom boundary **
    associate(&
     ! input: model control
    !  deriv_desired => in_qDrainFlux % deriv_desired, &          ! flag to indicate if derivatives are desired
     ixRichards    => decisions % f_Richards   , &          ! index defining the option for Richards' equation (moisture or mixdform)
     ! input: state and diagnostic variables
    !  nodeMatricHeadLiq => in_qDrainFlux % nodeMatricHeadLiq, &  ! liquid matric head in the lowest unsaturated node (m)
    !  nodeVolFracLiq    => in_qDrainFlux % nodeVolFracLiq   , &  ! volumetric liquid water content in the lowest unsaturated node (-)
     ! input: model coordinate variables
     nodeHeight => prog_data%mLayerHeight, &                ! height of the lowest unsaturated soil node (m)
     ! input: derivative in soil water characteristic
     node_dPsi_dTheta => deriv_data%mLayerdPsi_dTheta_m, &    ! derivative of the soil moisture characteristic w.r.t. theta (m)
     ! input: transmittance
     surfaceSatHydCond => flux_data%iLayerSatHydCond_m, &  ! saturated hydraulic conductivity at the surface (m s-1)
     ! input: soil parameters
     vGn_alpha       => mpar_data%vGn_alpha      , &      ! van Genuchten "alpha" parameter (m-1)
     vGn_n           => mpar_data%vGn_n          , &      ! van Genuchten "n" parameter (-)
     vGn_m           => diag_data%scalarVGn_m_m          , &      ! van Genuchten "m" parameter (-)
     theta_sat       => mpar_data%theta_sat      , &      ! soil porosity (-)
     theta_res       => mpar_data%theta_res      , &      ! soil residual volumetric water content (-)
     kAnisotropic    => mpar_data%kAnisotropic   , &      ! anisotropy factor for lateral hydraulic conductivity (-)
     zScale_TOPMODEL => mpar_data%zScale_TOPMODEL, &      ! scale factor for TOPMODEL-ish baseflow parameterization (m)
     ! output: drainage flux from the bottom of the soil profile
     scalarDrainage => flux_data%iLayerLiqFluxSoil_m, &       ! drainage flux from the bottom of the soil profile (m s-1)
     ! output: derivatives in drainage flux w.r.t. ...
     dq_dHydStateUnsat => deriv_data%dq_dHydStateAbove_m, & ! ... state variable in lowest unsaturated node (m s-1 or s-1)
     dq_dNrgStateUnsat => deriv_data%dq_dNrgStateAbove_m & ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
     ! output: error control
    !  err     => out_qDrainFlux % err    , &                     ! error code
    !  message => out_qDrainFlux % message  &                     ! error message
    &)
    !$cuf kernel do(2) <<<*,*>>>
    do iGRU=1,nGRU
      do iLayer=nSoil,nSoil
  
    call qDrainFlux_funcBottomHead(scalarDrainage(iLayer,iGRU), dq_dHydStateUnsat(iLayer,iGRU), dq_dNrgStateUnsat(iLayer,iGRU),&
    nodeMatricHeadLiq(iLayer,iGRU),nodeVolFracLiq(iLayer+nSnow(iGRU),iGRU),nodeHeight(iLayer+nSnow(iGRU),iGRU),&
    node_dPsi_dTheta(iLayer,iGRU),surfaceSatHydCond(0,iGRU), &
    vGn_alpha(iLayer), vGn_n(iLayer), vGn_m(iLayer,iGRU), theta_sat(iLayer), theta_res(iLayer), kAnisotropic, zScale_TOPMODEL, &
    ixRichards)
    end do
  end do
    ! energy derivatives
    err=20; message=trim(message)//"not yet implemented energy derivatives"; return_flag=.true.; return
  
    end associate
   end subroutine update_qDrainFlux_funcBottomHead
  
   subroutine update_qDrainFlux_freeDrainage
    ! ** Update operations for qDrainFlux: free drainage at bottom boundary **
    associate(&
     ! input: model control
    !  deriv_desired => in_qDrainFlux % deriv_desired, &          ! flag to indicate if derivatives are desired
     ixRichards    => decisions % f_Richards   , &          ! index defining the option for Richards' equation (moisture or mixdform)
     ! input: transmittance
     nodeHydCond       => flux_data%mLayerHydCond_m    , &    ! hydraulic conductivity at the node itself (m s-1)
     ! input: transmittance derivatives
    !  dHydCond_dVolLiq => in_qDrainFlux % dHydCond_dVolLiq, &    ! derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
    !  dHydCond_dMatric => in_qDrainFlux % dHydCond_dMatric, &    ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
    !  dHydCond_dTemp   => in_qDrainFlux % dHydCond_dTemp  , &    ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
     ! input: soil parameters
     kAnisotropic    => mpar_data%kAnisotropic  , &       ! anisotropy factor for lateral hydraulic conductivity (-)
     ! output: drainage flux from the bottom of the soil profile
     scalarDrainage => flux_data%iLayerLiqFluxSoil_m, &       ! drainage flux from the bottom of the soil profile (m s-1)
     ! output: derivatives in drainage flux w.r.t. ...
     dq_dHydStateUnsat => deriv_data%dq_dHydStateAbove_m, & ! ... state variable in lowest unsaturated node (m s-1 or s-1)
     dq_dNrgStateUnsat => deriv_data%dq_dNrgStateAbove_m & ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
     ! output: error control
    !  err     => out_qDrainFlux % err    , &                     ! error code
    !  message => out_qDrainFlux % message  &                     ! error message
    &)
      !$cuf kernel do(2) <<<*,*>>>
    do iGRU=1,nGRU
      do iLayer=nSoil,nSoil
  
    call qDrainFlux_freeDrainage(scalarDrainage(iLayer,iGRU),dq_dHydStateUnsat(iLayer,iGRU),dq_dNrgStateUnsat(iLayer,iGRU),&
    nodeHydCond(iLayer,iGRU), dHydCond_dVolLiq(iLayer,iGRU),dHydCond_dMatric(iLayer,iGRU),dHydCond_dTemp(iLayer,iGRU),&
    kAnisotropic,ixRichards)
    end do
  end do
    end associate
   end subroutine update_qDrainFlux_freeDrainage
  
   subroutine update_qDrainFlux_zeroFlux
    ! ** Update operations for qDrainFlux: zero flux condition at bottom boundary **
    associate(&
     ! input: model control
    !  deriv_desired => in_qDrainFlux % deriv_desired, &          ! flag to indicate if derivatives are desired
     ! output: drainage flux from the bottom of the soil profile
     scalarDrainage => flux_data%iLayerLiqFluxSoil_m, &       ! drainage flux from the bottom of the soil profile (m s-1)
     ! output: derivatives in drainage flux w.r.t. ...
     dq_dHydStateUnsat => deriv_data%dq_dHydStateAbove_m, & ! ... state variable in lowest unsaturated node (m s-1 or s-1)
     dq_dNrgStateUnsat => deriv_data%dq_dNrgStateAbove_m  & ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
    &)
    !$cuf kernel do(2) <<<*,*>>>
    do iGRU=1,nGRU
      do iLayer=nSoil,nSoil
  
    call qDrainFlux_zeroFlux(scalarDrainage(iLayer,iGRU),dq_dHydStateUnsat(iLayer,iGRU),dq_dNrgStateUnsat(iLayer,iGRU))
      end do
    end do
    end associate
   end subroutine update_qDrainFlux_zeroFlux
  
   subroutine finalize_qDrainFlux
    ! ** Finalize operations for qDrainFlux **
    ! associate(&
     ! output: error control
    !  err     => out_qDrainFlux % err    , & ! error code
    !  message => out_qDrainFlux % message  & ! error message
    ! &)
     ! final error check
     if (err /= 0_i4b) then
      message=trim(message)//'unanticipated error in qDrainFlux'
      return_flag=.true.; return
     end if
    ! end associate
   end subroutine finalize_qDrainFlux
  
  end subroutine qDrainFlux
  
  attributes(device) subroutine qDrainFlux_prescribedHead(bottomHydCond, bottomDiffuse, &
    scalarDrainage, dq_dHydStateUnsat, dq_dNrgStateUnsat, &
    nodeMatricHeadLiq, nodeVolFracLiq, nodeDepth, &
    lowerBoundHead, lowerBoundTheta, &
    bottomSatHydCond, iceImpedeFac, &
    dHydCond_dTemp, &
    vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, &
    ixRichards)
    ! USE soil_utils_module,only:volFracLiq  ! compute volumetric fraction of liquid water as a function of matric head (-)
    ! USE soil_utils_module,only:matricHead  ! compute matric head as a function of volumetric fraction of liquid water (m)
    USE soil_utils_module,only:hydCond_psi ! compute hydraulic conductivity as a function of matric head (m s-1)
    USE soil_utils_module,only:hydCond_liq ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
    USE soil_utils_module,only:dPsi_dTheta ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
  
    real(rkind),intent(inout) :: bottomHydCond, bottomDiffuse
    real(rkind),intent(inout) :: scalarDrainage, dq_dHydStateUnsat, dq_dNrgStateUnsat
    real(rkind),intent(in) :: nodeMatricHeadLiq, nodeVolFracLiq, nodeDepth
    real(rkind),intent(in) :: lowerBoundHead, lowerBoundTheta
    real(rkind),intent(in) :: bottomSatHydCond, iceImpedeFac
    real(rkind),intent(in) :: dHydCond_dTemp
    real(rkind),intent(in) :: vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res
    integer(i4b),intent(in) :: ixRichards
  
    ! local variables
    real(rkind)                      :: zWater                  ! effective water table depth (m)
    real(rkind)                      :: nodePsi                 ! matric head in the lowest unsaturated node (m)
    real(rkind)                      :: cflux                   ! capillary flux (m s-1)
    ! compute flux
    select case(ixRichards)
      case(moisture) ! moisture-based form of Richards' equation
        ! compute the hydraulic conductivity and diffusivity at the boundary
        bottomHydCond = hydCond_liq(lowerBoundTheta,bottomSatHydCond,theta_res,theta_sat,vGn_m) * iceImpedeFac
        bottomDiffuse = dPsi_dTheta(lowerBoundTheta,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * bottomHydCond
        ! compute the capillary flux
        cflux = -bottomDiffuse*(lowerBoundTheta - nodeVolFracLiq) / (nodeDepth*0.5_rkind)
      case(mixdform) ! mixed form of Richards' equation
        ! compute the hydraulic conductivity and diffusivity at the boundary
        bottomHydCond = hydCond_psi(lowerBoundHead,bottomSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
        bottomDiffuse = realMissing
        ! compute the capillary flux
        cflux = -bottomHydCond*(lowerBoundHead  - nodeMatricHeadLiq) / (nodeDepth*0.5_rkind)
      ! case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
    end select 
    scalarDrainage = cflux + bottomHydCond
  
    ! hydrology derivatives
    select case(ixRichards)  ! select form of Richards' equation
      case(moisture); dq_dHydStateUnsat = bottomDiffuse/(nodeDepth/2._rkind)
      case(mixdform); dq_dHydStateUnsat = bottomHydCond/(nodeDepth/2._rkind)
      ! case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
    end select
    ! energy derivatives
    dq_dNrgStateUnsat = -(dHydCond_dTemp/2._rkind)*(lowerBoundHead  - nodeMatricHeadLiq)/(nodeDepth*0.5_rkind)&
                      & + dHydCond_dTemp/2._rkind
  end subroutine
  
  
  attributes(device) subroutine qDrainFlux_funcBottomHead(scalarDrainage,dq_dHydStateUnsat,dq_dNrgStateUnsat, &
    nodeMatricHeadLiq, nodeVolFracLiq, nodeHeight, &
    node_dPsi_dTheta, surfaceSatHydCond, &
    vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, kAnisotropic, zScale_TOPMODEL, &
    ixRichards)
    ! USE soil_utils_module,only:volFracLiq  ! compute volumetric fraction of liquid water as a function of matric head (-)
    USE soil_utils_module,only:matricHead  ! compute matric head as a function of volumetric fraction of liquid water (m)
    ! USE soil_utils_module,only:hydCond_psi ! compute hydraulic conductivity as a function of matric head (m s-1)
    ! USE soil_utils_module,only:hydCond_liq ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
    ! USE soil_utils_module,only:dPsi_dTheta ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
    real(rkind),intent(inout) :: scalarDrainage, dq_dHydStateUnsat, dq_dNrgStateUnsat
    real(rkind),intent(in) :: nodeMatricHeadLiq, nodeVolFracLiq, nodeHeight
    real(rkind),intent(in) :: node_dPsi_dTheta, surfaceSatHydCond
    real(rkind),intent(in) :: vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, kAnisotropic, zScale_TOPMODEL
    integer(i4b),intent(in) :: ixRichards
    ! local variables
    real(rkind)                      :: zWater                  ! effective water table depth (m)
    real(rkind)                      :: nodePsi                 ! matric head in the lowest unsaturated node (m)
    real(rkind)                      :: cflux                   ! capillary flux (m s-1)
    
     ! compute flux
    select case(ixRichards) ! select form of Richards' equation
      case(moisture); nodePsi = matricHead(nodeVolFracLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
      case(mixdform); nodePsi = nodeMatricHeadLiq
    end select
    zWater = nodeHeight - nodePsi
    scalarDrainage = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)
  
    ! hydrology derivatives
    select case(ixRichards)  ! select form of Richards' equation
      case(moisture); dq_dHydStateUnsat = kAnisotropic*surfaceSatHydCond * node_dPsi_dTheta*exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
      case(mixdform); dq_dHydStateUnsat = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
      ! case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
    end select
    ! energy derivatives
    ! err=20; message=trim(message)//"not yet implemented energy derivatives"; return_flag=.true.; return
  
  end subroutine
  
  attributes(device) subroutine qDrainFlux_freeDrainage(scalarDrainage, dq_dHydStateUnsat, dq_dNrgStateUnsat, &
    nodeHydCond, &
    dHydCond_dVolLiq, dHydCond_dMatric, dHydCond_dTemp, &
    kAnisotropic, &
    ixRichards)
    ! USE soil_utils_module,only:volFracLiq  ! compute volumetric fraction of liquid water as a function of matric head (-)
    ! USE soil_utils_module,only:matricHead  ! compute matric head as a function of volumetric fraction of liquid water (m)
    ! USE soil_utils_module,only:hydCond_psi ! compute hydraulic conductivity as a function of matric head (m s-1)
    ! USE soil_utils_module,only:hydCond_liq ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
    ! USE soil_utils_module,only:dPsi_dTheta ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
    real(rkind),intent(inout) :: scalarDrainage, dq_dHydStateUnsat, dq_dNrgStateUnsat
    real(rkind),intent(in) :: nodeHydCond
    real(rkind),intent(in) :: dHydCond_dVolLiq, dHydCond_dMatric, dHydCond_dTemp
    real(rkind),intent(in) :: kAnisotropic
    integer(i4b),intent(in) :: ixRichards
    ! local variables
    ! real(rkind)                      :: zWater                  ! effective water table depth (m)
    ! real(rkind)                      :: nodePsi                 ! matric head in the lowest unsaturated node (m)
    ! real(rkind)                      :: cflux                   ! capillary flux (m s-1)
  
    scalarDrainage = nodeHydCond*kAnisotropic ! compute flux
  
    ! hydrology derivatives
    select case(ixRichards)  ! select form of Richards' equation
      case(moisture); dq_dHydStateUnsat = dHydCond_dVolLiq*kAnisotropic
      case(mixdform); dq_dHydStateUnsat = dHydCond_dMatric*kAnisotropic
      ! case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
    end select
    ! energy derivatives
    dq_dNrgStateUnsat = dHydCond_dTemp*kAnisotropic
  end subroutine
  
  attributes(device) subroutine qDrainFlux_zeroFlux(scalarDrainage, dq_dHydStateUnsat, dq_dNrgStateUnsat)
    real(rkind),intent(inout) :: scalarDrainage, dq_dHydStateUnsat, dq_dNrgStateUnsat
    scalarDrainage = 0._rkind
    dq_dHydStateUnsat = 0._rkind
    dq_dNrgStateUnsat = 0._rkind
  end subroutine
  end module soilLiqFlx_module
  