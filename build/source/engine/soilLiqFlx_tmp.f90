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
  USE globalData,only:veryBig                ! a very big number
  USE globalData,only:verySmall              ! a small number used as an additive constant to check if substantial difference among real numbers
  USE globalData,only:verySmaller            ! a smaller number used as an additive constant to check if substantial difference among real numbers

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
    zeroFlux,                   & ! zero flux
    zero_IE,                    & ! zero infiltration excess surface runoff parameterization 
    zero_SE,                    & ! zero saturation excess surface runoff parameterization 
    homegrown_IE,               & ! homegrown infiltration excess surface runoff parameterization 
    homegrown_SE,               & ! homegrown saturation excess surface runoff parameterization 
    FUSEPRMS,                   & ! FUSE PRMS     surface runoff parameterization 
    FUSEAVIC,                   & ! FUSE ARNO/VIC surface runoff parameterization
    FUSETOPM,                   & ! FUSE TOPMODEL surface runoff parameterization 
    ! look-up values for the maximum infiltration rate parameterization
    GreenAmpt,                  & ! Green-Ampt parameterization
    topmodel_GA,                & ! Green-Ampt parameterization with conductivity profile from TOPMODEL-ish parameterization  
    noInfiltrationExcess          ! no infiltration excess runoff
  
  ! -----------------------------------------------------------------------------------------------------------
  implicit none
  private
  public :: soilLiqFlx
  ! constant parameters
  ! real(rkind),parameter     :: dx=1.e-8_rkind               ! finite difference increment
  ! flag to denote if updating infiltration during iterations for testing purposes
  logical(lgt),parameter :: updateInfil=.true. 

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
                        ! above_soilLiqFluxDeriv,above_soildLiq_dTk,above_soilFracLiq, &
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
    ! real(rkind),intent(in),device :: above_soilLiqFluxDeriv(:),above_soildLiq_dTk(:),above_soilFracLiq(:)
    real(rkind),intent(in),device :: mLayerTempTrial(:,:), mLayerVolFracLiqTrial(:,:), mLayerVolFracIceTrial(:,:)
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
    ! real(rkind),dimension(0:nSoil,nGRU),device  :: iLayerHeight            ! height of the layer interfaces (m)
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
    call update_soilLiqFlx(decisions,mpar_data,prog_data,diag_data,flux_data,deriv_data,indx_data%nSnow,nGRU,indx_data%nSoil,iceImpedeFac,mLayerDiffuse,dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dTemp);     if (return_flag) return
  
    ! ** Final error control **
    call finalize_soilLiqFlx;   if (return_flag) return
  
   
  contains
  
   subroutine initialize_soilLiqFlx
    integer(i4b) :: iGRU
    ! **** Initial operations for soilLiqFlx module subroutine ****
  
    ! ** initialize error control **
    return_flag=.false.
    err=0; message='soilLiqFlx/' ! initialize error control
  
    ! ** get the indices for the soil layers **
    ixTop = 1
    ixBot = nSoil
  
    ! ** identify the number of layers that contain roots **
    associate(&
     rootingDepth => mpar_data%rootingDepth,& ! intent(in): rooting depth (m)
     iLayerHeight => prog_data%iLayerHeight, &
     nSnow => indx_data%nSnow &
    &) 
    nRoots = 0
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
      do iLayer=0,nSoil-1
        if (iLayerHeight(iLayer+nSnow(iGRU),iGRU) < rootingDepth-verySmall) nRoots(iGRU) = nRoots(iGRU) + 1
      end do
    end do
  
    end associate
  
    ! ** identify lowest soil layer with ice **
    ! NOTE: cannot use count because there may be an unfrozen wedge
    associate(&
      nSnow => indx_data%nSnow &
    &)
     ixIce = 0  ! initialize the index of the ice layer (0 means no ice in the soil profile)
     !$cuf kernel do(1) <<<*,*>>>
     do iGRU=1,nGRU
     do iLayer=1,nSoil ! (loop through soil layers)
       if (mLayerVolFracIceTrial(iLayer+nSnow(iGRU),iGRU) > verySmaller) ixIce(iGRU) = iLayer
     end do
    end do
  
    end associate
   end subroutine initialize_soilLiqFlx
  
   subroutine update_soilLiqFlx(decisions,mpar_data,prog_data,diag_data,flux_data,deriv_data,nSnow_d,nGRU,nSoil_d,iceImpedeFac,mLayerDiffuse,dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dTemp)
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
    real(rkind),device :: iceImpedeFac(:,:)
    real(rkind),device :: mLayerDiffuse(:,:)
    real(rkind),device :: dHydCond_dVolLiq(:,:)
    real(rkind),device :: dDiffuse_dVolLiq(:,:)
    real(rkind),device :: dHydCond_dTemp(:,:)
    type(dim3) :: blocks_2d, threads_2d
    type(dim3) :: blocks, threads
    threads = dim3(128,1,1)
    blocks = dim3(nGRU/128+1,1,1)
    threads_2d = dim3(4,8,1)
    blocks_2d = dim3(nGRU/4+1,nSoil/8+1,1)

    call compute_transpiration_sink_kernel<<<blocks_2d,threads_2d>>>(nGRU, nSoil, &
      diag_data%scalarTranspireLim, diag_data%mLayerRootDensity_m, diag_data%mLayerTranspireLim_m, mLayerTranspireFrac, &
      decisions%bcUpprSoiH, &
      flux_data%scalarCanopyTranspiration, flux_data%mLayerTranspire_m, &
      deriv_data%mLayerdTrans_dCanWat_m, deriv_data%mLayerdTrans_dTCanair_m, deriv_data%mLayerdTrans_dTCanopy_m, deriv_data%mLayerdTrans_dTGround_m, &
      deriv_data%dCanopyTrans_dCanWat, deriv_data%dCanopyTrans_dTCanair, deriv_data%dCanopyTrans_dTCanopy, deriv_data%dCanopyTrans_dTGround)

    call diagv_node_kernel<<<blocks_2d,threads_2d>>>(nGRU, nSoil, nSnow_d, &
      decisions%f_Richards,&
      mLayerVolFracLiqTrial,mLayerMatricHeadLiqTrial,mLayerVolFracIceTrial, &
        deriv_data%mLayerdTheta_dTk_m, deriv_data%dPsiLiq_dTemp_m, &
        mpar_data%vGn_alpha, mpar_data%vGn_n, diag_data%scalarVGn_m_m, mpar_data%mpExp, mpar_data%theta_sat, mpar_data%theta_res, mpar_data%theta_mp, mpar_data%f_impede,&
        flux_data%mLayerSatHydCond_m, flux_data%mLayerSatHydCondMP_m, &
        deriv_data%mLayerdPsi_dTheta_m, deriv_data%mLayerdTheta_dPsi_m, &
        flux_data%mLayerHydCond_m,mLayerDiffuse,iceImpedeFac, &
        dHydCond_dVolLiq, dDiffuse_dVolLiq,dHydCond_dMatric,dHydCond_dTemp)
  
    call compute_surface_infiltration_kernel<<<blocks,threads>>>(nGRU,firstSplitOper, decisions%bcUpprSoiH, decisions%surfRun_IE, decisions%surfRun_SE,decisions%f_Richards,decisions%infRateMax,nRoots,nSoil,ixIce,nSnow_d,indx_data%nLayers_d,&
      mpar_data%FUSE_Ac_max, mpar_data%FUSE_phi_tens, mpar_data%FUSE_b, mpar_data%FUSE_lambda, mpar_data%FUSE_chi, mpar_data%FUSE_mu, mpar_data%FUSE_n, &
      flux_data%scalarRainPlusMelt, &
      mLayerVolFracLiqTrial, mLayerTempTrial, mLayerMatricHeadTrial, mLayerVolFracIceTrial, &
      prog_data%mLayerDepth, prog_data%iLayerHeight, &
      mpar_data%upperBoundHead, mpar_data%upperBoundTheta, &
      flux_data%iLayerSatHydCond_m, dHydCond_dTemp, iceImpedeFac, &
      mpar_data%vGn_alpha, mpar_data%vGn_n, diag_data%scalarVGn_m_m, mpar_data%theta_sat, mpar_data%theta_res, mpar_data%rootingDepth, &
      mpar_data%zScale_TOPMODEL, mpar_data%wettingFrontSuction, mpar_data%qSurfScale, mpar_data%soilIceScale, mpar_data%soilIceCV, &
      deriv_data%mLayerdTheta_dTk_m, deriv_data%mLayerdTheta_dPsi_m, deriv_data%mLayerdPsi_dTheta_m, &
      iLayerHydCond, iLayerDiffuse, &
      flux_data%scalarSurfaceRunoff, flux_data%scalarSurfaceRunoff_IE, flux_data%scalarSurfaceRunoff_SE, flux_data%scalarInfiltration, &
      flux_data%scalarMaxInfilRate, diag_data%scalarInfilArea, diag_data%scalarFrozenArea, &
      diag_data%scalarSoilControl, deriv_data%dq_dHydStateLayerSurfVec_m, deriv_data%dq_dNrgStateLayerSurfVec_m, &
      flux_data%iLayerLiqFluxSoil_m,deriv_data%dq_dHydStateAbove_m,deriv_data%dq_dNrgStateAbove_m, &
      flux_data%scalarGroundEvaporation, deriv_data%dq_dHydStateBelow_m, deriv_data%dq_dNrgStateBelow_m)

    call iLayerFlux_kernel<<<blocks_2d,threads_2d>>>(decisions%f_Richards, nSnow_d, nGRU, nSoil, &
      mLayerMatricHeadLiqTrial, mLayerVolFracLiqTrial, &
      prog_data%mLayerHeight, &
      deriv_data%dPsiLiq_dTemp_m, dHydCond_dTemp, &
      flux_data%mLayerHydCond_m, mLayerDiffuse, &
      dHydCond_dVolLiq, dDiffuse_dVolLiq, dHydCond_dMatric, &
      iLayerHydCond, iLayerDiffuse, &
      flux_data%iLayerLiqFluxSoil_m, &
      deriv_data%dq_dHydStateAbove_m, deriv_data%dq_dHydStateBelow_m, &
      deriv_data%dq_dNrgStateAbove_m, deriv_data%dq_dNrgStateBelow_m)

    call compute_drainage_flux_kernel<<<blocks,threads>>>(nGRU, nSoil, nSnow_d, decisions%bcLowrSoiH, decisions%f_Richards, &
      mLayerMatricHeadLiqTrial, mLayerVolFracLiqTrial, &
      prog_data%mLayerDepth, prog_data%mLayerHeight, &
      mpar_data%lowerBoundHead, mpar_data%lowerBoundTheta, &
      deriv_data%mLayerdPsi_dTheta_m, deriv_data%dPsiLiq_dTemp_m, &
      flux_data%iLayerSatHydCond_m, iceImpedeFac, &
      dHydCond_dVolLiq, dHydCond_dMatric, dHydCond_dTemp, &
      mpar_data%vGn_alpha, mpar_data%vGn_n, diag_data%scalarVGn_m_m, mpar_data%theta_sat, mpar_data%theta_res, mpar_data%kAnisotropic, mpar_data%zScale_TOPMODEL, &
      flux_data%mLayerHydCond_m, mLayerDiffuse, &
      flux_data%iLayerLiqFluxSoil_m, &
      deriv_data%dq_dHydStateAbove_m, deriv_data%dq_dNrgStateAbove_m, &
      deriv_data%dq_dHydStateBelow_m, deriv_data%dq_dNrgStateBelow_m)
  
   end subroutine update_soilLiqFlx
  
   subroutine finalize_soilLiqFlx
    ! **** Final operations for soilLiqFlx module subroutine ****
   
    ! final error control check for robustness
     if (err/=0) then; message=trim(message)//trim("finalize_soilLiqFlx: final error check failed"); return_flag=.true.; return; end if
   end subroutine finalize_soilLiqFlx
  
    
  end subroutine soilLiqFlx

  attributes(global) subroutine compute_transpiration_sink_kernel(nGRU, nSoil, &
    scalarTranspireLim, mLayerRootDensity, mLayerTranspireLim, mLayerTranspireFrac, &
    ixBcUpperSoilHydrology, &
    scalarCanopyTranspiration, mLayerTranspire, &
    mLayerdTrans_dCanWat, mLayerdTrans_dTCanair, mLayerdTrans_dTCanopy, mLayerdTrans_dTGround, &
    dCanopyTrans_dCanWat, dCanopyTrans_dTCanair, dCanopyTrans_dTCanopy, dCanopyTrans_dTGround)
    integer(i4b),value :: nGRU, nSoil
    real(rkind),intent(in) :: scalarTranspireLim(:)  ! intent(in): weighted average of the transpiration limiting factor (-)
    real(rkind),intent(in) :: mLayerRootDensity(:,:)   ! intent(in): root density in each layer (-)
    real(rkind),intent(in) :: mLayerTranspireLim(:,:)  ! intent(in): transpiration limiting factor in each layer (-)
    real(rkind),intent(inout) :: mLayerTranspireFrac(:,:)
    ! intent(in): index of the upper boundary conditions for soil hydrology
    integer(i4b),intent(in) :: ixBcUpperSoilHydrology
    ! intent(inout)
    real(rkind),intent(in) :: scalarCanopyTranspiration(:)  ! canopy transpiration (kg m-2 s-1)
    real(rkind),intent(inout) :: mLayerTranspire(:,:)         ! transpiration loss from each soil layer (m s-1)
    ! intent(inout): derivatives in the soil layer transpiration flux ...
    real(rkind),intent(inout) :: mLayerdTrans_dCanWat(:,:)    ! ... w.r.t. canopy total water
    real(rkind),intent(inout) :: mLayerdTrans_dTCanair(:,:)   ! ... w.r.t. canopy air temperature
    real(rkind),intent(inout) :: mLayerdTrans_dTCanopy(:,:)   ! ... w.r.t. canopy temperature
    real(rkind),intent(inout) :: mLayerdTrans_dTGround(:,:)   ! ... w.r.t. ground temperature
    ! intent(in): derivative in canopy transpiration ...
    real(rkind),intent(in) :: dCanopyTrans_dCanWat(:)       ! ... w.r.t. canopy total water content (s-1)
    real(rkind),intent(in) :: dCanopyTrans_dTCanair(:)      ! ... w.r.t. canopy air temperature (kg m-2 s-1 K-1)
    real(rkind),intent(in) :: dCanopyTrans_dTCanopy(:)      ! ... w.r.t. canopy temperature (kg m-2 s-1 K-1)
    real(rkind),intent(in) :: dCanopyTrans_dTGround(:)      ! ... w.r.t. ground temperature (kg m-2 s-1 K-1)

    integer(i4b) :: iGRU, iSoil

    iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

    if (iGRU .gt. nGRU) return
    do iSoil=1,nSoil
    call update_transpiration_loss_fraction(scalarTranspireLim(iGRU), &
    mLayerRootDensity(iSoil,iGRU), mLayerRootDensity(:,iGRU), mLayerTranspireLim(iSoil,iGRU), mLayerTranspireFrac(iSoil,iGRU))

    call update_transpiration_loss(ixBcUpperSoilHydrology, &
    scalarCanopyTranspiration(iGRU), mLayerTranspire(iSoil,iGRU), mLayerTranspireFrac(iSoil,iGRU), &
    mLayerdTrans_dCanWat(iSoil,iGRU), mLayerdTrans_dTCanair(iSoil,iGRU), mLayerdTrans_dTCanopy(iSoil,iGRU), mLayerdTrans_dTGround(iSoil,iGRU), &
    dCanopyTrans_dCanWat(iGRU), dCanopyTrans_dTCanair(iGRU), dCanopyTrans_dTCanopy(iGRU), dCanopyTrans_dTGround(iGRU))
    end do

  end subroutine compute_transpiration_sink_kernel

  attributes(device) subroutine update_transpiration_loss_fraction(scalarTranspireLim, &
    mLayerRootDensity, mLayerRootDensity_GRU, mLayerTranspireLim, mLayerTranspireFrac)
    implicit none
    real(rkind),intent(in) :: scalarTranspireLim  ! intent(in): weighted average of the transpiration limiting factor (-)
    real(rkind),intent(in) :: mLayerRootDensity   ! intent(in): root density in each layer (-)
    real(rkind),intent(in) :: mLayerRootDensity_GRU(:)   ! intent(in): root density in each layer (-)
    real(rkind),intent(in) :: mLayerTranspireLim  ! intent(in): transpiration limiting factor in each layer (-)
    real(rkind),intent(inout) :: mLayerTranspireFrac

    ! transpiration may be non-zero even if the soil moisture limiting factor is zero
    if (scalarTranspireLim > tiny(scalarTranspireLim)) then 
      mLayerTranspireFrac = mLayerRootDensity*mLayerTranspireLim/scalarTranspireLim
    else ! possibility of non-zero conductance and therefore transpiration in this case
      mLayerTranspireFrac = mLayerRootDensity / sum(mLayerRootDensity_GRU)
    end if
  end subroutine update_transpiration_loss_fraction

  attributes(device) subroutine update_transpiration_loss(ixBcUpperSoilHydrology, &
    scalarCanopyTranspiration, mLayerTranspire, mLayerTranspireFrac, &
    mLayerdTrans_dCanWat, mLayerdTrans_dTCanair, mLayerdTrans_dTCanopy, mLayerdTrans_dTGround, &
    dCanopyTrans_dCanWat, dCanopyTrans_dTCanair, dCanopyTrans_dTCanopy, dCanopyTrans_dTGround)
   implicit none
   ! intent(in): index of the upper boundary conditions for soil hydrology
   integer(i4b),intent(in) :: ixBcUpperSoilHydrology
   ! intent(inout)
   real(rkind),intent(in) :: scalarCanopyTranspiration  ! canopy transpiration (kg m-2 s-1)
   real(rkind),intent(inout) :: mLayerTranspire         ! transpiration loss from each soil layer (m s-1)
   real(rkind),intent(in) :: mLayerTranspireFrac
   ! intent(inout): derivatives in the soil layer transpiration flux ...
   real(rkind),intent(inout) :: mLayerdTrans_dCanWat    ! ... w.r.t. canopy total water
   real(rkind),intent(inout) :: mLayerdTrans_dTCanair   ! ... w.r.t. canopy air temperature
   real(rkind),intent(inout) :: mLayerdTrans_dTCanopy   ! ... w.r.t. canopy temperature
   real(rkind),intent(inout) :: mLayerdTrans_dTGround   ! ... w.r.t. ground temperature
   ! intent(in): derivative in canopy transpiration ...
   real(rkind),intent(in) :: dCanopyTrans_dCanWat       ! ... w.r.t. canopy total water content (s-1)
   real(rkind),intent(in) :: dCanopyTrans_dTCanair      ! ... w.r.t. canopy air temperature (kg m-2 s-1 K-1)
   real(rkind),intent(in) :: dCanopyTrans_dTCanopy      ! ... w.r.t. canopy temperature (kg m-2 s-1 K-1)
   real(rkind),intent(in) :: dCanopyTrans_dTGround      ! ... w.r.t. ground temperature (kg m-2 s-1 K-1)

   if (ixBcUpperSoilHydrology==prescribedHead) then ! special case of prescribed head -- no transpiration
    mLayerTranspire      = 0._rkind
    ! derivatives in transpiration w.r.t. canopy state variables
    mLayerdTrans_dCanWat = 0._rkind
    mLayerdTrans_dTCanair= 0._rkind
    mLayerdTrans_dTCanopy= 0._rkind
    mLayerdTrans_dTGround= 0._rkind
   else
    mLayerTranspire = mLayerTranspireFrac*scalarCanopyTranspiration/iden_water
    ! * derivatives in transpiration w.r.t. canopy state variables *
    mLayerdTrans_dCanWat  = mLayerTranspireFrac*dCanopyTrans_dCanWat /iden_water
    mLayerdTrans_dTCanair = mLayerTranspireFrac*dCanopyTrans_dTCanair/iden_water
    mLayerdTrans_dTCanopy = mLayerTranspireFrac*dCanopyTrans_dTCanopy/iden_water
    mLayerdTrans_dTGround = mLayerTranspireFrac*dCanopyTrans_dTGround/iden_water
   end if
  end subroutine update_transpiration_loss


  attributes(global) subroutine compute_surface_infiltration_kernel(nGRU,firstSplitOper, bc_upper, surfRun_IE, surfRun_SE,ixRichards,ixInfRateMax,nRoots,nSoil,ixIce,nSnow,nLayers,&
  FUSE_Ac_max, FUSE_phi_tens, FUSE_b, FUSE_lambda, FUSE_chi, FUSE_mu, FUSE_n, &
  scalarRainPlusMelt, &
  mLayerVolFracLiq, mLayerTemp, mLayerMatricHead, mLayerVolFracIce, &
  mLayerDepth, iLayerHeight, &
  upperBoundHead, upperBoundTheta, &
  surfaceSatHydCond, dHydCond_dTemp, iceImpedeFac, &
  vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, rootingDepth, &
  zScale_TOPMODEL, wettingFrontSuction, qSurfScale, soilIceScale, soilIceCV, &
  dTheta_dTk, dTheta_dPsi, mLayerdPsi_dTheta, &
  surfaceHydCond, surfaceDiffuse, &
  scalarSurfaceRunoff, scalarSurfaceRunoff_IE, scalarSurfaceRunoff_SE, scalarSurfaceInfiltration, &
  xMaxInfilRate, scalarInfilArea, scalarFrozenArea, &
  scalarSoilControl, dq_dHydStateVec, dq_dNrgStateVec, &
  iLayerLiqFluxSoil,dq_dHydStateAbove,dq_dNrgStateAbove, &
  scalarGroundEvaporation, dq_dHydStateBelow, dq_dNrgStateBelow)
    integer(i4b),value :: nGRU
      logical(lgt),intent(in),value :: firstSplitOper ! flag indicating if desire to compute infiltration
    integer(i4b),intent(in) :: bc_upper ! index defining the type of boundary conditions
    integer(i4b),intent(in) :: surfRun_IE ! index defining the infiltration excess surface runoff method
    integer(i4b),intent(in) :: surfRun_SE ! index defining the saturation excess surface runoff method
   integer(i4b),intent(in) :: ixRichards ! index defining the option for Richards' equation (moisture or mixdform)
           integer(i4b),intent(in) :: ixInfRateMax ! index defining the maximum infiltration rate method (GreenAmpt, topmodel_GA, noInfiltrationExcess)
   integer(i4b),intent(in) :: nRoots(:)          ! number of layers that contain roots
      integer(i4b),intent(in),value :: nSoil           ! number of soil layers
   integer(i4b),intent(in) :: ixIce(:)           ! index of lowest ice layer
   integer(i4b),intent(in) :: nSnow(:),nLayers(:)
   ! input: parameters
   real(rkind),intent(in) :: FUSE_Ac_max, FUSE_phi_tens
   real(rkind),intent(in) :: FUSE_b
   real(rkind),intent(in) :: FUSE_lambda
   real(rkind),intent(in) :: FUSE_chi
   real(rkind),intent(in) :: FUSE_mu
   real(rkind),intent(in) :: FUSE_n
   ! input: rain plus melt
   real(rkind),intent(in) :: scalarRainPlusMelt(:) ! rain plus melt  (m s-1)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: mLayerVolFracLiq(:,:)  ! volumetric liquid water content in each soil layer (-)
   real(rkind),intent(in) :: mLayerTemp(:,:)          ! temperature (K)
   real(rkind),intent(in) :: mLayerMatricHead(:,:)    ! matric head in each soil layer (m)
      real(rkind),intent(in) :: mLayerVolFracIce(:,:)     ! volumetric ice content in each soil layer (-)

   ! input: depth of upper-most soil layer (m)
   real(rkind),intent(in) :: mLayerDepth(:,:)   ! depth of upper-most soil layer (m)
     real(rkind),intent(in) :: iLayerHeight(0:,:)      ! height at the interface of each layer for soil layers only (m)

   ! input: diriclet boundary conditions
   real(rkind),intent(in) :: upperBoundHead    ! upper boundary condition for matric head (m)
   real(rkind),intent(in) :: upperBoundTheta   ! upper boundary condition for volumetric liquid water content (-)
   ! input: transmittance
   real(rkind),intent(in) :: surfaceSatHydCond(0:,:)  ! saturated hydraulic conductivity at the surface (m s-1)
   real(rkind),intent(in) :: dHydCond_dTemp(:,:)     ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   real(rkind),intent(in) :: iceImpedeFac(:,:)       ! ice impedence factor in the upper-most soil layer (-)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha(:)            ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n(:)                ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m(:,:)                ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: theta_sat(:)            ! soil porosity (-)
   real(rkind),intent(in) :: theta_res(:)            ! soil residual volumetric water content (-)
      real(rkind),intent(in) :: rootingDepth         ! rooting depth (m)
     real(rkind),intent(in) :: zScale_TOPMODEL      ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
  real(rkind),intent(in) :: wettingFrontSuction  ! Green-Ampt wetting front suction (m)
     real(rkind),intent(in) :: qSurfScale ! scaling factor in the surface runoff parameterization (-)
   real(rkind),intent(in) :: soilIceScale         ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
   real(rkind),intent(in) :: soilIceCV            ! soil ice CV in Gamma distribution used to define frozen area (-)

     ! input: pre-computed derivatives in ...
  real(rkind),intent(in) :: dTheta_dTk(:,:)          ! ... volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in) :: dTheta_dPsi(:,:)         ! ... the soil water characteristic w.r.t. psi (m-1)
   real(rkind),intent(in) :: mLayerdPsi_dTheta(:,:)      ! ... the soil water characteristic w.r.t. theta (m)

   ! input-output: hydraulic conductivity and diffusivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   real(rkind),intent(inout) :: surfaceHydCond(0:,:)  ! hydraulic conductivity (m s-1)
   real(rkind),intent(inout) :: surfaceDiffuse(0:,:)  ! hydraulic diffusivity at the surface (m2 s-1)
   ! output: runoff and infiltration
   real(rkind),intent(inout) :: scalarSurfaceRunoff(:)        ! surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceRunoff_IE(:)     ! infiltration excess surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceRunoff_SE(:)     ! saturation excess surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceInfiltration(:)  ! surface infiltration (m s-1)
   ! input-output
   real(rkind),intent(inout) :: xMaxInfilRate(:), scalarInfilArea(:), scalarFrozenArea(:)
      ! output: derivatives in surface infiltration w.r.t. ...
   real(rkind),intent(inout) :: scalarSoilControl(:)   ! soil control on infiltration for derivative
   real(rkind),intent(inout) :: dq_dHydStateVec(:,:)     ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateVec(:,:)     ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
   real(rkind),intent(inout) :: iLayerLiqFluxSoil(0:,:), dq_dHydStateAbove(0:,:), dq_dNrgStateAbove(0:,:)
   real(rkind),intent(inout) :: scalarGroundEvaporation(:)
   real(rkind),intent(inout) :: dq_dHydStateBelow(0:,:), dq_dNrgStateBelow(0:,:)

    integer(i4b) :: iGRU
    iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
    if (iGRU .gt. nGRU) return

    call initialize_compute_surface_infiltration(dq_dHydStateAbove(:,iGRU), dq_dNrgStateAbove(:,iGRU))
    call surfaceFlx_kernel(firstSplitOper, bc_upper, surfRun_IE, surfRun_SE,ixRichards,ixInfRateMax,nRoots,nSoil,ixIce,nSnow,nLayers,nGRU,&
  FUSE_Ac_max, FUSE_phi_tens, FUSE_b, FUSE_lambda, FUSE_chi, FUSE_mu, FUSE_n, &
  scalarRainPlusMelt, &
  mLayerVolFracLiq, mLayerTemp, mLayerMatricHead, mLayerVolFracIce, &
  mLayerDepth, iLayerHeight, &
  upperBoundHead, upperBoundTheta, &
  surfaceSatHydCond, dHydCond_dTemp, iceImpedeFac, &
  vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, rootingDepth, &
  zScale_TOPMODEL, wettingFrontSuction, qSurfScale, soilIceScale, soilIceCV, &
  dTheta_dTk, dTheta_dPsi, mLayerdPsi_dTheta, &
  surfaceHydCond, surfaceDiffuse, &
  scalarSurfaceRunoff, scalarSurfaceRunoff_IE, scalarSurfaceRunoff_SE, scalarSurfaceInfiltration, &
  xMaxInfilRate, scalarInfilArea, scalarFrozenArea, &
  scalarSoilControl, dq_dHydStateVec, dq_dNrgStateVec)
    call finalize_compute_surface_infiltration(iLayerLiqFluxSoil(:,iGRU), scalarGroundEvaporation(iGRU), scalarSurfaceInfiltration(iGRU), dq_dHydStateBelow(:,iGRU), dq_dNrgStateBelow(:,iGRU))
  end subroutine
  attributes(device) subroutine initialize_compute_surface_infiltration(dq_dHydStateAbove, dq_dNrgStateAbove)
    real(rkind) :: dq_dHydStateAbove(0:), dq_dNrgStateAbove(0:)

    dq_dHydStateAbove(0) = 0._rkind
    dq_dNrgStateAbove(0) = 0._rkind
  end subroutine
  attributes(device) subroutine finalize_compute_surface_infiltration(iLayerLiqFluxSoil, scalarGroundEvaporation, scalarSurfaceInfiltration, dq_dHydStateBelow, dq_dNrgStateBelow)
    real(rkind) :: iLayerLiqFluxSoil(0:), dq_dHydStateBelow(0:), dq_dNrgStateBelow(0:)
    real(rkind) :: scalarGroundEvaporation, scalarSurfaceInfiltration

    iLayerLiqFluxSoil(0) = scalarGroundEvaporation/iden_water + scalarSurfaceInfiltration

    dq_dHydStateBelow(0) = 0._rkind ! contribution will be in dq_dHydStateLayerSurfVec(1)
    dq_dNrgStateBelow(0) = 0._rkind ! contribution will be in dq_dNrgStateLayerSurfVec(1)
  end subroutine


  attributes(global) subroutine compute_drainage_flux_kernel(nGRU, nSoil, nSnow, bc_lower, ixRichards, &
    mLayerMatricHeadLiqTrial, mLayerVolFracLiqTrial, &
    nodeDepth, nodeHeight, &
    lowerBoundHead, lowerBoundTheta, &
    node_dPsi_dTheta, dPsiLiq_dTemp, &
    bottomSatHydCond, iceImpedeFac, &
    dHydCond_dVolLiq, dHydCond_dMatric, dHydCond_dTemp, &
    vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, kAnisotropic, zScale_TOPMODEL, &
    mLayerHydCond, mLayerDiffuse, &
    scalarDrainage, &
    dq_dHydStateUnsat, dq_dNrgStateUnsat, &
    dq_dHydStateBelow, dq_dNrgStateBelow)
    ! input: model control
    integer(i4b),intent(in),value :: nGRU, nSoil
    integer(i4b),intent(in) :: nSnow(:)
   integer(i4b),intent(in) :: bc_lower ! index defining the type of boundary conditions
   integer(i4b),intent(in) :: ixRichards          ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: mLayerMatricHeadLiqTrial(:,:)   ! liquid matric head in the lowest unsaturated node (m)
   real(rkind),intent(in) :: mLayerVolFracLiqTrial(:,:)      ! volumetric liquid water content in the lowest unsaturated node (-)
   ! input: model coordinate variables
   real(rkind),intent(in) :: nodeDepth(:,:)           ! depth of the lowest unsaturated soil layer (m)
   real(rkind),intent(in) :: nodeHeight(0:,:)                ! height of the lowest unsaturated soil node (m)
   ! input: diriclet boundary conditions
   real(rkind),intent(in) :: lowerBoundHead      ! lower boundary condition for matric head (m)
   real(rkind),intent(in) :: lowerBoundTheta     ! lower boundary condition for volumetric liquid water content (-)
   ! input: derivative in soil water characteristic
   real(rkind),intent(in) :: node_dPsi_dTheta(:,:)      ! derivative of the soil moisture characteristic w.r.t. theta (m)
   real(rkind),intent(in) :: dPsiLiq_dTemp(:,:)    ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: transmittance
   real(rkind),intent(in) :: bottomSatHydCond(0:,:)    ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   real(rkind),intent(in) :: iceImpedeFac(:,:)        ! ice impedence factor in the upper-most soil layer (-)
   ! input: transmittance derivatives
   real(rkind),intent(in) :: dHydCond_dVolLiq(:,:)     ! derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
   real(rkind),intent(in) :: dHydCond_dMatric(:,:)     ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
   real(rkind),intent(in) :: dHydCond_dTemp(:,:)      ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha(:)           ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n(:)               ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m(:,:)               ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: theta_sat(:)           ! soil porosity (-)
   real(rkind),intent(in) :: theta_res(:)           ! soil residual volumetric water content (-)
   real(rkind),intent(in) :: kAnisotropic           ! anisotropy factor for lateral hydraulic conductivity (-)
   real(rkind),intent(in) :: zScale_TOPMODEL       ! scale factor for TOPMODEL-ish baseflow parameterization (m)
   ! output: hydraulic conductivity at the bottom of the unsaturated zone
   real(rkind),intent(inout) :: mLayerHydCond(:,:)       ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   real(rkind),intent(inout) :: mLayerDiffuse(:,:)       ! hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
   ! output: drainage flux from the bottom of the soil profile
   real(rkind),intent(inout) :: scalarDrainage(0:,:)      ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   real(rkind),intent(inout) :: dq_dHydStateUnsat(0:,:)   ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateUnsat(0:,:)   ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
   real(rkind) :: dq_dHydStateBelow(:,:), dq_dNrgStateBelow(:,:)


    integer(i4b) :: iGRU
    iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
    if (iGRU .gt. nGRU) return

    call qDrainFlux_kernel(nGRU, nSoil, nSnow,bc_lower,ixRichards, &
 mLayerMatricHeadLiqTrial, mLayerVolFracLiqTrial, &
 nodeDepth, nodeHeight, &
 lowerBoundHead, lowerBoundTheta, &
 node_dPsi_dTheta, dPsiLiq_dTemp, &
 bottomSatHydCond, iceImpedeFac, bottomSatHydCond, &
 dHydCond_dVolLiq, dHydCond_dMatric, dHydCond_dTemp, &
 vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, kAnisotropic, zScale_TOPMODEL, &
 mLayerHydCond, mLayerDiffuse, &
 scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat)

    call finalize_compute_drainage_flux(dq_dHydStateBelow(:,iGRU), dq_dNrgStateBelow(:,iGRU), nSoil)

  end subroutine
  attributes(device) subroutine finalize_compute_drainage_flux(dq_dHydStateBelow, dq_dNrgStateBelow, nSoil)
    integer(i4b) :: nSoil
    real(rkind) :: dq_dHydStateBelow(0:), dq_dNrgStateBelow(0:)

    dq_dHydStateBelow(nSoil) = 0._rkind  ! keep this here in case we want to couple some day....
    dq_dNrgStateBelow(nSoil) = 0._rkind  ! keep this here in case we want to couple some day....
  end subroutine



  

  ! ***************************************************************************************************************
  ! private subroutine diagv_node: compute transmittance and derivatives for model nodes
  ! ***************************************************************************************************************

  attributes(global) subroutine diagv_node_kernel(nGRU, nSoil, nSnow, &
      ixRichards,&
      mLayerVolFracLiqTrial,mLayerMatricHeadLiqTrial,mLayerVolFracIceTrial, &
        dTheta_dTk, dPsiLiq_dTemp, &
        vGn_alpha, vGn_n, vGn_m, mpExp, theta_sat, theta_res, theta_mp, f_impede,&
        mLayerSatHydCond, mLayerSatHydCondMP, &
        mLayerdPsi_dTheta, mLayerdTheta_dPsi, &
        mLayerHydCond,mLayerDiffuse,iceImpedeFac, &
        dHydCond_dVolLiq, dDiffuse_dVolLiq,dHydCond_dMatric,dHydCond_dTemp)
    integer(i4b),value :: nGRU, nSoil
    integer(i4b) :: nSnow(:)
    ! input: model control
    integer(i4b) :: ixRichards
    ! input: state and diagnostic variables
    real(rkind),intent(in) :: mLayerVolFracLiqTrial(:,:) ! volumetric fraction of liquid water in a given layer (-)
    real(rkind),intent(in) :: mLayerMatricHeadLiqTrial(:,:) ! liquid matric head in each layer (m)
    real(rkind),intent(in) :: mLayerVolFracIceTrial(:,:) ! volumetric fraction of ice in a given layer (-)
    ! input: pre-computed derivatives
    real(rkind),intent(in) :: dTheta_dTk(:,:) ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
    real(rkind),intent(in) :: dPsiLiq_dTemp(:,:) ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
    ! input: soil parameters
    real(rkind),intent(in) :: vGn_alpha(:) ! van Genuchten "alpha" parameter (m-1)
    real(rkind),intent(in) :: vGn_n(:) ! van Genuchten "n" parameter (-)
    real(rkind),intent(in) :: vGn_m(:,:) ! van Genuchten "m" parameter (-)
    real(rkind),intent(in) :: mpExp ! empirical exponent in macropore flow equation (-)
    real(rkind),intent(in) :: theta_sat(:) ! soil porosity (-)
    real(rkind),intent(in) :: theta_res(:) ! soil residual volumetric water content (-)
    real(rkind),intent(in) :: theta_mp ! volumetric liquid water content when macropore flow begins (-)
    real(rkind),intent(in) :: f_impede ! ice impedence factor (-)
    ! input: saturated hydraulic conductivity ...
    real(rkind),intent(in) :: mLayerSatHydCond(:,:) ! ... at the mid-point of a given layer (m s-1)
    real(rkind),intent(in) :: mLayerSatHydCondMP(:,:) ! ... of macropores at the mid-point of a given layer (m s-1)
    ! output: derivative in the soil water characteristic
    real(rkind),intent(inout) :: mLayerdPsi_dTheta(:,:) ! derivative in the soil water characteristic
    real(rkind),intent(inout) :: mLayerdTheta_dPsi(:,:) ! derivative in the soil water characteristic
    ! output: transmittance
    real(rkind),intent(inout) :: mLayerHydCond(:,:) ! hydraulic conductivity at layer mid-points (m s-1)
    real(rkind),intent(inout) :: mLayerDiffuse(:,:) ! diffusivity at layer mid-points (m2 s-1)
    real(rkind),intent(inout) :: iceImpedeFac(:,:) ! ice impedence factor in each layer (-)
    ! output: transmittance derivatives in ...
    real(rkind),intent(inout) :: dHydCond_dVolLiq(:,:) ! ... hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
    real(rkind),intent(inout) :: dDiffuse_dVolLiq(:,:) ! ... hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
    real(rkind),intent(inout) :: dHydCond_dMatric(:,:) ! ... hydraulic conductivity w.r.t matric head (s-1)
    real(rkind),intent(inout) :: dHydCond_dTemp(:,:) ! ... hydraulic conductivity w.r.t temperature (m s-1 K-1)
    integer(i4b) :: iGRU, iSoil

     iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
     iSoil = (blockidx%y-1) * blockdim%y + threadidx%y

     if (iGRU .gt. nGRU .or. iSoil .gt. nSoil) return
     call diagv_node(ixRichards,&
        mLayerVolFracLiqTrial(iSoil+nSnow(iGRU),iGRU),mLayerMatricHeadLiqTrial(iSoil,iGRU),mLayerVOlFracIceTrial(iSoil+nSnow(iGRU),iGRU), &
        dTheta_dTk(iSoil+nSnow(iGRU),iGRU), dPsiLiq_dTemp(iSoil,iGRU), &
        vGn_alpha(iSoil), vGn_n(iSoil), vGn_m(iSoil,iGRU), mpExp, theta_sat(iSoil), theta_res(iSoil), theta_mp, f_impede,&
        mLayerSatHydCond(iSoil,iGRU), mLayerSatHydCondMP(iSoil,iGRU), &
        mLayerdPsi_dTheta(iSoil,iGRU), mLayerdTheta_dPsi(iSoil,iGRU), &
        mLayerHydCond(iSoil,iGRU),mLayerDiffuse(iSoil,iGRU),iceImpedeFac(iSoil,iGRU), &
        dHydCond_dVolLiq(iSoil,iGRU), dDiffuse_dVolLiq(iSoil,iGRU),dHydCond_dMatric(iSoil,iGRU),dHydCond_dTemp(iSoil,iGRU))
  end subroutine diagv_node_kernel
  attributes(device) subroutine diagv_node(ixRichards,&
  scalarVolFracLiqTrial,scalarMatricHeadLiqTrial,scalarVolFracIceTrial,&
  dTheta_dTk, dPsiLiq_dTemp, &
  vGn_alpha, vGn_n, vGn_m, mpExp, theta_sat, theta_res, theta_mp, f_impede, &
  scalarSatHydCond, scalarSatHydCondMP, &
  scalardPsi_dTheta, scalardTheta_dPsi, &
  scalarHydCond, scalarDiffuse, iceImpedeFac, &
  dHydCond_dVolLiq, dDiffuse_dVolLiq, dHydCond_dMatric, dHydCond_dTemp &
  ) 
    use device_data_types
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
    ! compute hydraulic transmittance and derivatives for all layers
    implicit none
    ! input: model control
    integer(i4b) :: ixRichards
    ! input: state and diagnostic variables
    real(rkind),intent(in) :: scalarVolFracLiqTrial ! volumetric fraction of liquid water in a given layer (-)
    real(rkind),intent(in) :: scalarMatricHeadLiqTrial ! liquid matric head in each layer (m)
    real(rkind),intent(in) :: scalarVolFracIceTrial ! volumetric fraction of ice in a given layer (-)
    ! input: pre-computed derivatives
    real(rkind),intent(in) :: dTheta_dTk ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
    real(rkind),intent(in) :: dPsiLiq_dTemp ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
    ! input: soil parameters
    real(rkind),intent(in) :: vGn_alpha ! van Genuchten "alpha" parameter (m-1)
    real(rkind),intent(in) :: vGn_n ! van Genuchten "n" parameter (-)
    real(rkind),intent(in) :: vGn_m ! van Genuchten "m" parameter (-)
    real(rkind),intent(in) :: mpExp ! empirical exponent in macropore flow equation (-)
    real(rkind),intent(in) :: theta_sat ! soil porosity (-)
    real(rkind),intent(in) :: theta_res ! soil residual volumetric water content (-)
    real(rkind),intent(in) :: theta_mp ! volumetric liquid water content when macropore flow begins (-)
    real(rkind),intent(in) :: f_impede ! ice impedence factor (-)
    ! input: saturated hydraulic conductivity ...
    real(rkind),intent(in) :: scalarSatHydCond ! ... at the mid-point of a given layer (m s-1)
    real(rkind),intent(in) :: scalarSatHydCondMP ! ... of macropores at the mid-point of a given layer (m s-1)
    ! output: derivative in the soil water characteristic
    real(rkind),intent(inout) :: scalardPsi_dTheta ! derivative in the soil water characteristic
    real(rkind),intent(inout) :: scalardTheta_dPsi ! derivative in the soil water characteristic
    ! output: transmittance
    real(rkind),intent(inout) :: scalarHydCond ! hydraulic conductivity at layer mid-points (m s-1)
    real(rkind),intent(inout) :: scalarDiffuse ! diffusivity at layer mid-points (m2 s-1)
    real(rkind),intent(inout) :: iceImpedeFac ! ice impedence factor in each layer (-)
    ! output: transmittance derivatives in ...
    real(rkind),intent(inout) :: dHydCond_dVolLiq ! ... hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
    real(rkind),intent(inout) :: dDiffuse_dVolLiq ! ... hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
    real(rkind),intent(inout) :: dHydCond_dMatric ! ... hydraulic conductivity w.r.t matric head (s-1)
    real(rkind),intent(inout) :: dHydCond_dTemp ! ... hydraulic conductivity w.r.t temperature (m s-1 K-1)
    ! local variables
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
    logical(lgt)                     :: return_flag               ! flag for return statements
  
  
      call update_diagv_node(&
    ixRichards, &
    scalarVolFracLiqTrial, &
    scalarMatricHeadLiqTrial, &
    scalarVolFracIceTrial, &
    dTheta_dTk, &
    dPsiLiq_dTemp, &
    vGn_alpha, &
    vGn_n, &
    vGn_m, &
    mpExp, &
    theta_sat, &
    theta_res, &
    theta_mp, &
    f_impede, &
    scalarSatHydCond, &
    scalarSatHydCondMP, &
    scalardPsi_dTheta, &
    scalardTheta_dPsi, &
    scalarHydCond, &
    scalarDiffuse, &
    iceImpedeFac, &
    dHydCond_dVolLiq, &
    dDiffuse_dVolLiq, &
    dHydCond_dMatric, &
    dHydCond_dTemp, &
    localVolFracLiq, &
    scalarHydCondMP, &
    dIceImpede_dT, &
    dHydCondMacro_dVolLiq, &
    dHydCondMacro_dMatric, &
    dHydCondMicro_dMatric, &
    dHydCondMicro_dTemp, &
    dPsi_dTheta2a, &
    dIceImpede_dLiq, &
    hydCond_noIce, &
    dK_dPsi__noIce, &
    dK_dLiq__noIce, &
    relSatMP)
        
  end subroutine diagv_node
  
    attributes(device) subroutine update_diagv_node(&
    ixRichards, &
    scalarVolFracLiqTrial, &
    scalarMatricHeadLiqTrial, &
    scalarVolFracIceTrial, &
    dTheta_dTk, &
    dPsiLiq_dTemp, &
    vGn_alpha, &
    vGn_n, &
    vGn_m, &
    mpExp, &
    theta_sat, &
    theta_res, &
    theta_mp, &
    f_impede, &
    scalarSatHydCond, &
    scalarSatHydCondMP, &
    scalardPsi_dTheta, &
    scalardTheta_dPsi, &
    scalarHydCond, &
    scalarDiffuse, &
    iceImpedeFac, &
    dHydCond_dVolLiq, &
    dDiffuse_dVolLiq, &
    dHydCond_dMatric, &
    dHydCond_dTemp, &
    localVolFracLiq, &
    scalarHydCondMP, &
    dIceImpede_dT, &
    dHydCondMacro_dVolLiq, &
    dHydCondMacro_dMatric, &
    dHydCondMicro_dMatric, &
    dHydCondMicro_dTemp, &
    dPsi_dTheta2a, &
    dIceImpede_dLiq, &
    hydCond_noIce, &
    dK_dPsi__noIce, &
    dK_dLiq__noIce, &
    relSatMP)
   ! input: model control
    integer(i4b),intent(in) :: ixRichards ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: scalarVolFracLiqTrial ! volumetric fraction of liquid water in a given layer (-)
   real(rkind),intent(in) :: scalarMatricHeadLiqTrial ! liquid matric head in each layer (m)
   real(rkind),intent(in) :: scalarVolFracIceTrial ! volumetric fraction of ice in a given layer (-)
   ! input: pre-computed derivatives
   real(rkind),intent(in) :: dTheta_dTk ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
   real(rkind),intent(in) :: dPsiLiq_dTemp ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: mpExp ! empirical exponent in macropore flow equation (-)
   real(rkind),intent(in) :: theta_sat ! soil porosity (-)
   real(rkind),intent(in) :: theta_res ! soil residual volumetric water content (-)
   real(rkind),intent(in) :: theta_mp ! volumetric liquid water content when macropore flow begins (-)
   real(rkind),intent(in) :: f_impede ! ice impedence factor (-)
   ! input: saturated hydraulic conductivity ...
   real(rkind),intent(in) :: scalarSatHydCond ! ... at the mid-point of a given layer (m s-1)
   real(rkind),intent(in) :: scalarSatHydCondMP ! ... of macropores at the mid-point of a given layer (m s-1)
   ! output: derivative in the soil water characteristic
   real(rkind),intent(inout) :: scalardPsi_dTheta ! derivative in the soil water characteristic
   real(rkind),intent(inout) :: scalardTheta_dPsi ! derivative in the soil water characteristic
   ! output: transmittance
   real(rkind),intent(inout) :: scalarHydCond ! hydraulic conductivity at layer mid-points (m s-1)
   real(rkind),intent(inout) :: scalarDiffuse ! diffusivity at layer mid-points (m2 s-1)
   real(rkind),intent(inout) :: iceImpedeFac ! ice impedence factor in each layer (-)
   ! output: transmittance derivatives in ...
   real(rkind),intent(inout) :: dHydCond_dVolLiq ! ... hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
   real(rkind),intent(inout) :: dDiffuse_dVolLiq ! ... hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
   real(rkind),intent(inout) :: dHydCond_dMatric ! ... hydraulic conductivity w.r.t matric head (s-1)
   real(rkind),intent(inout) :: dHydCond_dTemp ! ... hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! inout: local to diagv_node
   real(rkind),intent(inout)                      :: localVolFracLiq           ! local volumetric fraction of liquid water
   real(rkind),intent(inout)                      :: scalarHydCondMP           ! hydraulic conductivity of macropores at layer mid-points (m s-1)
   real(rkind),intent(inout)                      :: dIceImpede_dT             ! derivative in ice impedance factor w.r.t. temperature (K-1)
   real(rkind),intent(inout)                      :: dHydCondMacro_dVolLiq     ! derivative in hydraulic conductivity of macropores w.r.t volumetric liquid water content (m s-1)
   real(rkind),intent(inout)                      :: dHydCondMacro_dMatric     ! derivative in hydraulic conductivity of macropores w.r.t matric head (s-1)
   real(rkind),intent(inout)                      :: dHydCondMicro_dMatric     ! derivative in hydraulic conductivity of micropores w.r.t matric head (s-1)
   real(rkind),intent(inout)                      :: dHydCondMicro_dTemp       ! derivative in hydraulic conductivity of micropores w.r.t temperature (m s-1 K-1)
   real(rkind),intent(inout)                      :: dPsi_dTheta2a             ! derivative in dPsi_dTheta (analytical)
   real(rkind),intent(inout)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
   real(rkind),intent(inout)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
   real(rkind),intent(inout)                      :: dK_dPsi__noIce            ! derivative in hydraulic conductivity w.r.t matric head, in the absence of ice (s-1)
   real(rkind),intent(inout)                      :: dK_dLiq__noIce            ! derivative in hydraulic conductivity w.r.t volumetric liquid water content, in the absence of ice (m s-1)
   real(rkind),intent(inout)                      :: relSatMP                  ! relative saturation of macropores (-)

  ! **** Update operations for diagv_node ****

   call update_diagv_node_characteristic_derivatives(&
    ixRichards, &
    scalarMatricHeadLiqTrial, &
    scalarVolFracLiqTrial, &
    vGn_alpha, &
    vGn_n, &
    vGn_m, &
    mpExp, &
    theta_sat, &
    theta_res, &
    scalardPsi_dTheta, &
    scalardTheta_dPsi &
  )

   call update_diagv_node_hydraulic_conductivity(&
    ixRichards, &
    scalarVolFracLiqTrial, &
    scalarMatricHeadLiqTrial, &
    scalarVolFracIceTrial, &
    dTheta_dTk, &
    dPsiLiq_dTemp, &
    vGn_alpha, &
    vGn_n, &
    vGn_m, &
    mpExp, &
    theta_sat, &
    theta_res, &
    theta_mp, &
    f_impede, &
    scalarSatHydCond, &
    scalarSatHydCondMP, &
    scalardPsi_dTheta, &
    scalardTheta_dPsi, &
    scalarHydCond, &
    scalarDiffuse, &
    iceImpedeFac, &
    dHydCond_dVolLiq, &
    dDiffuse_dVolLiq, &
    dHydCond_dMatric, &
    dHydCond_dTemp, &
    localVolFracLiq, &
    scalarHydCondMP, &
    dIceImpede_dT, &
    dHydCondMacro_dVolLiq, &
    dHydCondMacro_dMatric, &
    dHydCondMicro_dMatric, &
    dHydCondMicro_dTemp, &
    dPsi_dTheta2a, &
    dIceImpede_dLiq, &
    hydCond_noIce, &
    dK_dPsi__noIce, &
    dK_dLiq__noIce, &
    relSatMP)

  end subroutine update_diagv_node


  attributes(device) subroutine update_diagv_node_characteristic_derivatives(&
    ixRichards, &
    scalarMatricHeadLiqTrial, &
    scalarVolFracLiqTrial, &
    vGn_alpha, &
    vGn_n, &
    vGn_m, &
    mpExp, &
    theta_sat, &
    theta_res, &
    scalardPsi_dTheta, &
    scalardTheta_dPsi &
  )
    USE soil_utils_module,only:dTheta_dPsi          ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
    USE soil_utils_module,only:dPsi_dTheta          ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
    USE soil_utils_module,only:dPsi_dTheta2         ! compute derivative in dPsi_dTheta (m)
    implicit none

   ! **** Update operations for diagv_node: compute characteristic derivatives ****
   ! compute the derivative in the soil water characteristic
   ! input: model control
   integer(i4b),intent(in) :: ixRichards     ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: scalarMatricHeadLiqTrial ! liquid matric head in each layer (m)
   real(rkind),intent(in) :: scalarVolFracLiqTrial ! volumetric fraction of liquid water in a given layer (-)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha  ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n      ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m      ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: mpExp      ! empirical exponent in macropore flow equation (-)
   real(rkind),intent(in) :: theta_sat  ! soil porosity (-)
   real(rkind),intent(in) :: theta_res  ! soil residual volumetric water content (-)
   ! output: derivative in the soil water characteristic
   real(rkind),intent(inout) :: scalardPsi_dTheta ! derivative in the soil water characteristic
   real(rkind),intent(inout) :: scalardTheta_dPsi ! derivative in the soil water characteristic
   select case(ixRichards)
     case(moisture)
       scalardPsi_dTheta = dPsi_dTheta(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
       scalardTheta_dPsi = realMissing  ! deliberately cause problems if this is ever used
     case(mixdform)
       scalardTheta_dPsi = dTheta_dPsi(scalarMatricHeadLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
       scalardPsi_dTheta = dPsi_dTheta(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
    !  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select
  end subroutine update_diagv_node_characteristic_derivatives

  attributes(device) subroutine update_diagv_node_hydraulic_conductivity(&
    ixRichards, &
    scalarVolFracLiqTrial, &
    scalarMatricHeadLiqTrial, &
    scalarVolFracIceTrial, &
    dTheta_dTk, &
    dPsiLiq_dTemp, &
    vGn_alpha, &
    vGn_n, &
    vGn_m, &
    mpExp, &
    theta_sat, &
    theta_res, &
    theta_mp, &
    f_impede, &
    scalarSatHydCond, &
    scalarSatHydCondMP, &
    scalardPsi_dTheta, &
    scalardTheta_dPsi, &
    scalarHydCond, &
    scalarDiffuse, &
    iceImpedeFac, &
    dHydCond_dVolLiq, &
    dDiffuse_dVolLiq, &
    dHydCond_dMatric, &
    dHydCond_dTemp, &
    localVolFracLiq, &
    scalarHydCondMP, &
    dIceImpede_dT, &
    dHydCondMacro_dVolLiq, &
    dHydCondMacro_dMatric, &
    dHydCondMicro_dMatric, &
    dHydCondMicro_dTemp, &
    dPsi_dTheta2a, &
    dIceImpede_dLiq, &
    hydCond_noIce, &
    dK_dPsi__noIce, &
    dK_dLiq__noIce, &
    relSatMP)
    USE soil_utils_module,only:iceImpede            ! compute the ice impedence factor
    integer(i4b),intent(in) :: ixRichards ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: scalarVolFracLiqTrial ! volumetric fraction of liquid water in a given layer (-)
   real(rkind),intent(in) :: scalarMatricHeadLiqTrial ! liquid matric head in each layer (m)
   real(rkind),intent(in) :: scalarVolFracIceTrial ! volumetric fraction of ice in a given layer (-)
   ! input: pre-computed derivatives
   real(rkind),intent(in) :: dTheta_dTk ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
   real(rkind),intent(in) :: dPsiLiq_dTemp ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: mpExp ! empirical exponent in macropore flow equation (-)
   real(rkind),intent(in) :: theta_sat ! soil porosity (-)
   real(rkind),intent(in) :: theta_res ! soil residual volumetric water content (-)
   real(rkind),intent(in) :: theta_mp ! volumetric liquid water content when macropore flow begins (-)
   real(rkind),intent(in) :: f_impede ! ice impedence factor (-)
   ! input: saturated hydraulic conductivity ...
   real(rkind),intent(in) :: scalarSatHydCond ! ... at the mid-point of a given layer (m s-1)
   real(rkind),intent(in) :: scalarSatHydCondMP ! ... of macropores at the mid-point of a given layer (m s-1)
   ! output: derivative in the soil water characteristic
   real(rkind),intent(inout) :: scalardPsi_dTheta ! derivative in the soil water characteristic
   real(rkind),intent(inout) :: scalardTheta_dPsi ! derivative in the soil water characteristic
   ! output: transmittance
   real(rkind),intent(inout) :: scalarHydCond ! hydraulic conductivity at layer mid-points (m s-1)
   real(rkind),intent(inout) :: scalarDiffuse ! diffusivity at layer mid-points (m2 s-1)
   real(rkind),intent(inout) :: iceImpedeFac ! ice impedence factor in each layer (-)
   ! output: transmittance derivatives in ...
   real(rkind),intent(inout) :: dHydCond_dVolLiq ! ... hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
   real(rkind),intent(inout) :: dDiffuse_dVolLiq ! ... hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
   real(rkind),intent(inout) :: dHydCond_dMatric ! ... hydraulic conductivity w.r.t matric head (s-1)
   real(rkind),intent(inout) :: dHydCond_dTemp ! ... hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! inout: local to diagv_node
   real(rkind),intent(inout)                      :: localVolFracLiq           ! local volumetric fraction of liquid water
   real(rkind),intent(inout)                      :: scalarHydCondMP           ! hydraulic conductivity of macropores at layer mid-points (m s-1)
   real(rkind),intent(inout)                      :: dIceImpede_dT             ! derivative in ice impedance factor w.r.t. temperature (K-1)
   real(rkind),intent(inout)                      :: dHydCondMacro_dVolLiq     ! derivative in hydraulic conductivity of macropores w.r.t volumetric liquid water content (m s-1)
   real(rkind),intent(inout)                      :: dHydCondMacro_dMatric     ! derivative in hydraulic conductivity of macropores w.r.t matric head (s-1)
   real(rkind),intent(inout)                      :: dHydCondMicro_dMatric     ! derivative in hydraulic conductivity of micropores w.r.t matric head (s-1)
   real(rkind),intent(inout)                      :: dHydCondMicro_dTemp       ! derivative in hydraulic conductivity of micropores w.r.t temperature (m s-1 K-1)
   real(rkind),intent(inout)                      :: dPsi_dTheta2a             ! derivative in dPsi_dTheta (analytical)
   real(rkind),intent(inout)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
   real(rkind),intent(inout)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
   real(rkind),intent(inout)                      :: dK_dPsi__noIce            ! derivative in hydraulic conductivity w.r.t matric head, in the absence of ice (s-1)
   real(rkind),intent(inout)                      :: dK_dLiq__noIce            ! derivative in hydraulic conductivity w.r.t volumetric liquid water content, in the absence of ice (m s-1)
   real(rkind),intent(inout)                      :: relSatMP                  ! relative saturation of macropores (-)

    ! **** Update operations for diagv_node: compute hydraulic conductivity and derivatives ****
    ! compute hydraulic conductivity and its derivative in each soil layer
    ! compute the ice impedence factor and its derivative w.r.t. volumetric liquid water content (-)
    call iceImpede(scalarVolFracIceTrial,f_impede, &  ! input
                   iceImpedeFac,dIceImpede_dLiq)     ! output
   select case(ixRichards)
     case(moisture) ! moisture-based form of Richards' equation
       call update_diagv_node_hydraulic_conductivity_moisture_form(&
    scalarVolFracLiqTrial, &
    scalarVolFracIceTrial, &
    vGn_alpha, &
    vGn_n, &
    vGn_m, &
    theta_sat, &
    theta_res, &
    scalarSatHydCond, &
    scalardPsi_dTheta, &
    scalarHydCond, &
    scalarDiffuse, &
    iceImpedeFac, &
    dHydCond_dVolLiq, &
    dDiffuse_dVolLiq, &
    dHydCond_dMatric, &
    hydCond_noIce, &
    dK_dLiq__noIce, &
    dIceImpede_dLiq, &
    dPsi_dTheta2a)
     case(mixdform) ! mixed form of Richards' equation -- just compute hydraulic condictivity
       call update_diagv_node_hydraulic_conductivity_mixed_form(&
    scalarMatricHeadLiqTrial, &
    scalarVolFracIceTrial, &
    dTheta_dTk, &
    dPsiLiq_dTemp, &
    vGn_alpha, &
    vGn_n, &
    vGn_m, &
    mpExp, &
    theta_sat, &
    theta_res, &
    theta_mp, &
    f_impede, &
    scalarSatHydCond, &
    scalarSatHydCondMP, &
    scalardTheta_dPsi, &
    scalarHydCond, &
    scalarDiffuse, &
    iceImpedeFac, &
    dHydCond_dVolLiq, &
    dDiffuse_dVolLiq, &
    dHydCond_dMatric, &
    dHydCond_dTemp, &
    localVolFracLiq, &
    scalarHydCondMP, &
    dIceImpede_dT, &
    dHydCondMacro_dVolLiq, &
    dHydCondMacro_dMatric, &
    dHydCondMicro_dMatric, &
    dHydCondMicro_dTemp, &
    dIceImpede_dLiq, &
    hydCond_noIce, &
    dK_dPsi__noIce, &
    relSatMP)
   end select 
 end subroutine update_diagv_node_hydraulic_conductivity

 attributes(device) subroutine update_diagv_node_hydraulic_conductivity_moisture_form(&
    scalarVolFracLiqTrial, &
    scalarVolFracIceTrial, &
    vGn_alpha, &
    vGn_n, &
    vGn_m, &
    theta_sat, &
    theta_res, &
    scalarSatHydCond, &
    scalardPsi_dTheta, &
    scalarHydCond, &
    scalarDiffuse, &
    iceImpedeFac, &
    dHydCond_dVolLiq, &
    dDiffuse_dVolLiq, &
    dHydCond_dMatric, &
    hydCond_noIce, &
    dK_dLiq__noIce, &
    dIceImpede_dLiq, &
    dPsi_dTheta2a)
        USE soil_utils_module,only:hydCond_liq          ! compute hydraulic conductivity as a function of volumetric liquid water content
    USE soil_utils_module,only:dHydCond_dLiq        ! compute derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
    USE soil_utils_module,only:dPsi_dTheta2         ! compute derivative in dPsi_dTheta (m)

   ! input: state and diagnostic variables
   real(rkind),intent(in) :: scalarVolFracLiqTrial ! volumetric fraction of liquid water in a given layer (-)
   real(rkind),intent(in) :: scalarVolFracIceTrial ! volumetric fraction of ice in a given layer (-)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: theta_sat ! soil porosity (-)
   real(rkind),intent(in) :: theta_res ! soil residual volumetric water content (-)
   ! input: saturated hydraulic conductivity ...
   real(rkind),intent(in) :: scalarSatHydCond ! ... at the mid-point of a given layer (m s-1)
   ! output: derivative in the soil water characteristic
   real(rkind),intent(inout) :: scalardPsi_dTheta ! derivative in the soil water characteristic
   ! output: transmittance
   real(rkind),intent(inout) :: scalarHydCond ! hydraulic conductivity at layer mid-points (m s-1)
   real(rkind),intent(inout) :: scalarDiffuse ! diffusivity at layer mid-points (m2 s-1)
   real(rkind),intent(inout) :: iceImpedeFac ! ice impedence factor in each layer (-)
   ! output: transmittance derivatives in ...
   real(rkind),intent(inout) :: dHydCond_dVolLiq ! ... hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
   real(rkind),intent(inout) :: dDiffuse_dVolLiq ! ... hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
   real(rkind),intent(inout) :: dHydCond_dMatric ! ... hydraulic conductivity w.r.t matric head (s-1)

   real(rkind),intent(inout)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
   real(rkind),intent(inout)                      :: dK_dLiq__noIce            ! derivative in hydraulic conductivity w.r.t volumetric liquid water content, in the absence of ice (m s-1)
   real(rkind),intent(inout)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
   real(rkind),intent(inout)                      :: dPsi_dTheta2a             ! derivative in dPsi_dTheta (analytical)


  ! **** Update operations for diagv_node: compute hydraulic conductivity and derivatives for moisture form of Richards' equation ****

  ! validation
   ! haven't included macropores yet -- return with error for now
   ! err=20; message=trim(message)//'still need to include macropores for the moisture-based form of Richards eqn'
   return

  ! computation
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

 end subroutine update_diagv_node_hydraulic_conductivity_moisture_form

 attributes(device) subroutine update_diagv_node_hydraulic_conductivity_mixed_form(&
    scalarMatricHeadLiqTrial, &
    scalarVolFracIceTrial, &
    dTheta_dTk, &
    dPsiLiq_dTemp, &
    vGn_alpha, &
    vGn_n, &
    vGn_m, &
    mpExp, &
    theta_sat, &
    theta_res, &
    theta_mp, &
    f_impede, &
    scalarSatHydCond, &
    scalarSatHydCondMP, &
    scalardTheta_dPsi, &
    scalarHydCond, &
    scalarDiffuse, &
    iceImpedeFac, &
    dHydCond_dVolLiq, &
    dDiffuse_dVolLiq, &
    dHydCond_dMatric, &
    dHydCond_dTemp, &
    localVolFracLiq, &
    scalarHydCondMP, &
    dIceImpede_dT, &
    dHydCondMacro_dVolLiq, &
    dHydCondMacro_dMatric, &
    dHydCondMicro_dMatric, &
    dHydCondMicro_dTemp, &
    dIceImpede_dLiq, &
    hydCond_noIce, &
    dK_dPsi__noIce, &
    relSatMP)
        USE soil_utils_module,only:hydCond_psi          ! compute hydraulic conductivity as a function of matric head
    USE soil_utils_module,only:volFracLiq           ! compute volumetric fraction of liquid water as a function of matric head
    USE soil_utils_module,only:hydCondMP_liq        ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
    USE soil_utils_module,only:dHydCond_dPsi        ! compute derivative in hydraulic conductivity w.r.t. matric head (s-1)
    USE soil_utils_module,only:dIceImpede_dTemp     ! compute the derivative in the ice impedance factor w.r.t. temperature (K-1)

   ! input: state and diagnostic variables
   real(rkind),intent(in) :: scalarMatricHeadLiqTrial ! liquid matric head in each layer (m)
   real(rkind),intent(in) :: scalarVolFracIceTrial ! volumetric fraction of ice in a given layer (-)
   ! input: pre-computed derivatives
   real(rkind),intent(in) :: dTheta_dTk ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
   real(rkind),intent(in) :: dPsiLiq_dTemp ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: mpExp ! empirical exponent in macropore flow equation (-)
   real(rkind),intent(in) :: theta_sat ! soil porosity (-)
   real(rkind),intent(in) :: theta_res ! soil residual volumetric water content (-)
   real(rkind),intent(in) :: theta_mp ! volumetric liquid water content when macropore flow begins (-)
   real(rkind),intent(in) :: f_impede ! ice impedence factor (-)
   ! input: saturated hydraulic conductivity ...
   real(rkind),intent(in) :: scalarSatHydCond ! ... at the mid-point of a given layer (m s-1)
   real(rkind),intent(in) :: scalarSatHydCondMP ! ... of macropores at the mid-point of a given layer (m s-1)
   ! output: derivative in the soil water characteristic
   real(rkind),intent(inout) :: scalardTheta_dPsi ! derivative in the soil water characteristic
   ! output: transmittance
   real(rkind),intent(inout) :: scalarHydCond ! hydraulic conductivity at layer mid-points (m s-1)
   real(rkind),intent(inout) :: scalarDiffuse ! diffusivity at layer mid-points (m2 s-1)
   real(rkind),intent(inout) :: iceImpedeFac ! ice impedence factor in each layer (-)
   ! output: transmittance derivatives in ...
   real(rkind),intent(inout) :: dHydCond_dVolLiq ! ... hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
   real(rkind),intent(inout) :: dDiffuse_dVolLiq ! ... hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
   real(rkind),intent(inout) :: dHydCond_dMatric ! ... hydraulic conductivity w.r.t matric head (s-1)
   real(rkind),intent(inout) :: dHydCond_dTemp ! ... hydraulic conductivity w.r.t temperature (m s-1 K-1)

   real(rkind),intent(inout)                      :: localVolFracLiq           ! local volumetric fraction of liquid water
   real(rkind),intent(inout)                      :: scalarHydCondMP           ! hydraulic conductivity of macropores at layer mid-points (m s-1)
   real(rkind),intent(inout)                      :: dIceImpede_dT             ! derivative in ice impedance factor w.r.t. temperature (K-1)
   real(rkind),intent(inout)                      :: dHydCondMacro_dVolLiq     ! derivative in hydraulic conductivity of macropores w.r.t volumetric liquid water content (m s-1)
   real(rkind),intent(inout)                      :: dHydCondMacro_dMatric     ! derivative in hydraulic conductivity of macropores w.r.t matric head (s-1)
   real(rkind),intent(inout)                      :: dHydCondMicro_dMatric     ! derivative in hydraulic conductivity of micropores w.r.t matric head (s-1)
   real(rkind),intent(inout)                      :: dHydCondMicro_dTemp       ! derivative in hydraulic conductivity of micropores w.r.t temperature (m s-1 K-1)
   real(rkind),intent(inout)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
   real(rkind),intent(inout)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
   real(rkind),intent(inout)                      :: dK_dPsi__noIce            ! derivative in hydraulic conductivity w.r.t matric head, in the absence of ice (s-1)
   real(rkind),intent(inout)                      :: relSatMP                  ! relative saturation of macropores (-)

  ! **** Update operations for diagv_node: compute hydraulic conductivity and derivatives for mixed form of Richards' equation ****

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
   if (scalarVolFracIceTrial > verySmaller) then
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

 end subroutine update_diagv_node_hydraulic_conductivity_mixed_form
  
  ! ***************************************************************************************************************
  ! private subroutine surfaceFlx: compute the surface flux and its derivative
  ! ***************************************************************************************************************
 attributes(device) subroutine surfaceFlx_kernel(firstSplitOper, bc_upper, surfRun_IE, surfRun_SE,ixRichards,ixInfRateMax,nRoots,nSoil,ixIce,nSnow,nLayers,nGRU,&
  FUSE_Ac_max, FUSE_phi_tens, FUSE_b, FUSE_lambda, FUSE_chi, FUSE_mu, FUSE_n, &
  scalarRainPlusMelt, &
  mLayerVolFracLiq, mLayerTemp, mLayerMatricHead, mLayerVolFracIce, &
  mLayerDepth, iLayerHeight, &
  upperBoundHead, upperBoundTheta, &
  surfaceSatHydCond, dHydCond_dTemp, iceImpedeFac, &
  vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, rootingDepth, &
  zScale_TOPMODEL, wettingFrontSuction, qSurfScale, soilIceScale, soilIceCV, &
  dTheta_dTk, dTheta_dPsi, mLayerdPsi_dTheta, &
  surfaceHydCond, surfaceDiffuse, &
  scalarSurfaceRunoff, scalarSurfaceRunoff_IE, scalarSurfaceRunoff_SE, scalarSurfaceInfiltration, &
  xMaxInfilRate, scalarInfilArea, scalarFrozenArea, &
  scalarSoilControl, dq_dHydStateVec, dq_dNrgStateVec)
      logical(lgt),intent(in),value :: firstSplitOper ! flag indicating if desire to compute infiltration
    integer(i4b),intent(in) :: bc_upper ! index defining the type of boundary conditions
    integer(i4b),intent(in) :: surfRun_IE ! index defining the infiltration excess surface runoff method
    integer(i4b),intent(in) :: surfRun_SE ! index defining the saturation excess surface runoff method
   integer(i4b),intent(in) :: ixRichards ! index defining the option for Richards' equation (moisture or mixdform)
           integer(i4b),intent(in) :: ixInfRateMax ! index defining the maximum infiltration rate method (GreenAmpt, topmodel_GA, noInfiltrationExcess)
   integer(i4b),intent(in) :: nRoots(:)          ! number of layers that contain roots
      integer(i4b),intent(in),value :: nSoil           ! number of soil layers
   integer(i4b),intent(in) :: ixIce(:)           ! index of lowest ice layer
   integer(i4b),intent(in) :: nSnow(:),nLayers(:)
   integer(i4b),intent(in),value :: nGRU
   ! input: parameters
   real(rkind),intent(in) :: FUSE_Ac_max, FUSE_phi_tens
   real(rkind),intent(in) :: FUSE_b
   real(rkind),intent(in) :: FUSE_lambda
   real(rkind),intent(in) :: FUSE_chi
   real(rkind),intent(in) :: FUSE_mu
   real(rkind),intent(in) :: FUSE_n
   ! input: rain plus melt
   real(rkind),intent(in) :: scalarRainPlusMelt(:) ! rain plus melt  (m s-1)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: mLayerVolFracLiq(:,:)  ! volumetric liquid water content in each soil layer (-)
   real(rkind),intent(in) :: mLayerTemp(:,:)          ! temperature (K)
   real(rkind),intent(in) :: mLayerMatricHead(:,:)    ! matric head in each soil layer (m)
      real(rkind),intent(in) :: mLayerVolFracIce(:,:)     ! volumetric ice content in each soil layer (-)

   ! input: depth of upper-most soil layer (m)
   real(rkind),intent(in) :: mLayerDepth(:,:)   ! depth of upper-most soil layer (m)
     real(rkind),intent(in) :: iLayerHeight(0:,:)      ! height at the interface of each layer for soil layers only (m)

   ! input: diriclet boundary conditions
   real(rkind),intent(in) :: upperBoundHead    ! upper boundary condition for matric head (m)
   real(rkind),intent(in) :: upperBoundTheta   ! upper boundary condition for volumetric liquid water content (-)
   ! input: transmittance
   real(rkind),intent(in) :: surfaceSatHydCond(0:,:)  ! saturated hydraulic conductivity at the surface (m s-1)
   real(rkind),intent(in) :: dHydCond_dTemp(:,:)     ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   real(rkind),intent(in) :: iceImpedeFac(:,:)       ! ice impedence factor in the upper-most soil layer (-)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha(:)            ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n(:)                ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m(:,:)                ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: theta_sat(:)            ! soil porosity (-)
   real(rkind),intent(in) :: theta_res(:)            ! soil residual volumetric water content (-)
      real(rkind),intent(in) :: rootingDepth         ! rooting depth (m)
     real(rkind),intent(in) :: zScale_TOPMODEL      ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
  real(rkind),intent(in) :: wettingFrontSuction  ! Green-Ampt wetting front suction (m)
     real(rkind),intent(in) :: qSurfScale ! scaling factor in the surface runoff parameterization (-)
   real(rkind),intent(in) :: soilIceScale         ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
   real(rkind),intent(in) :: soilIceCV            ! soil ice CV in Gamma distribution used to define frozen area (-)

     ! input: pre-computed derivatives in ...
  real(rkind),intent(in) :: dTheta_dTk(:,:)          ! ... volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in) :: dTheta_dPsi(:,:)         ! ... the soil water characteristic w.r.t. psi (m-1)
   real(rkind),intent(in) :: mLayerdPsi_dTheta(:,:)      ! ... the soil water characteristic w.r.t. theta (m)

   ! input-output: hydraulic conductivity and diffusivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   real(rkind),intent(inout) :: surfaceHydCond(0:,:)  ! hydraulic conductivity (m s-1)
   real(rkind),intent(inout) :: surfaceDiffuse(0:,:)  ! hydraulic diffusivity at the surface (m2 s-1)
   ! output: runoff and infiltration
   real(rkind),intent(inout) :: scalarSurfaceRunoff(:)        ! surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceRunoff_IE(:)     ! infiltration excess surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceRunoff_SE(:)     ! saturation excess surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceInfiltration(:)  ! surface infiltration (m s-1)
   ! input-output
   real(rkind),intent(inout) :: xMaxInfilRate(:), scalarInfilArea(:), scalarFrozenArea(:)
      ! output: derivatives in surface infiltration w.r.t. ...
   real(rkind),intent(inout) :: scalarSoilControl(:)   ! soil control on infiltration for derivative
   real(rkind),intent(inout) :: dq_dHydStateVec(:,:)     ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateVec(:,:)     ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)

   integer(i4b) :: iGRU

   iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return

  call surfaceFlx(firstSplitOper,bc_upper, surfRun_IE, surfRun_SE,ixRichards,ixInfRateMax,nRoots(iGRU),nSoil,ixIce(iGRU),nSnow(iGRU),nLayers(iGRU),&
  FUSE_Ac_max, FUSE_phi_tens, FUSE_b, FUSE_lambda, FUSE_chi, FUSE_mu, FUSE_n, &
  scalarRainPlusMelt(iGRU), &
  mLayerMatricHead(1,iGRU), mLayerVolFracLiq(nSnow(iGRU)+1,iGRU), mLayerVolFracLiq(:,iGRU), mLayerTemp(:,iGRU), mLayerMatricHead(:,iGRU), mLayerVolFracIce(:,iGRU), &
  mLayerDepth(:,iGRU), iLayerHeight(nSnow(iGRU):nLayers(iGRU),iGRU), &
  upperBoundHead, upperBoundTheta, &
  surfaceSatHydCond(0,iGRU), dHydCond_dTemp(1,iGRU), iceImpedeFac(1,iGRU), &
  vGn_alpha(1), vGn_n(1), vGn_m(1,iGRU), theta_sat(1), theta_res(1), rootingDepth, &
  zScale_TOPMODEL, wettingFrontSuction, qSurfScale, soilIceScale, soilIceCV, &
  dTheta_dTk(:,iGRU), dTheta_dPsi(:,iGRU), mLayerdPsi_dTheta(:,iGRU), &
  surfaceHydCond(0,iGRU), surfaceDiffuse(0,iGRU), &
  scalarSurfaceRunoff(iGRU), scalarSurfaceRunoff_IE(iGRU), scalarSurfaceRunoff_SE(iGRU), scalarSurfaceInfiltration(iGRU), &
  xMaxInfilRate(iGRU), scalarInfilArea(iGRU), scalarFrozenArea(iGRU), &
  scalarSoilControl(iGRU), dq_dHydStateVec(:,iGRU), dq_dNrgStateVec(:,iGRU))

  end subroutine


  attributes(device) subroutine surfaceFlx(firstSplitOper,bc_upper, surfRun_IE, surfRun_SE,ixRichards,ixInfRateMax,nRoots,nSoil,ixIce,nSnow,nLayers,&
  FUSE_Ac_max, FUSE_phi_tens, FUSE_b, FUSE_lambda, FUSE_chi, FUSE_mu, FUSE_n, &
  scalarRainPlusMelt, &
  scalarMatricHeadLiq, scalarVolFracLiq, mLayerVolFracLiq, mLayerTemp, mLayerMatricHead, mLayerVolFracIce, &
  mLayerDepth, iLayerHeight, &
  upperBoundHead, upperBoundTheta, &
  surfaceSatHydCond, dHydCond_dTemp, iceImpedeFac, &
  vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, rootingDepth, &
  zScale_TOPMODEL, wettingFrontSuction, qSurfScale, soilIceScale, soilIceCV, &
  dTheta_dTk, dTheta_dPsi, mLayerdPsi_dTheta, &
  surfaceHydCond, surfaceDiffuse, &
  scalarSurfaceRunoff, scalarSurfaceRunoff_IE, scalarSurfaceRunoff_SE, scalarSurfaceInfiltration, &
  xMaxInfilRate, scalarInfilArea, scalarFrozenArea, &
  scalarSoilControl, dq_dHydStateVec, dq_dNrgStateVec)
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
        ! input: model control
    logical(lgt),intent(in) :: firstSplitOper ! flag indicating if desire to compute infiltration
    integer(i4b),intent(in) :: bc_upper ! index defining the type of boundary conditions
    integer(i4b),intent(in) :: surfRun_IE ! index defining the infiltration excess surface runoff method
    integer(i4b),intent(in) :: surfRun_SE ! index defining the saturation excess surface runoff method
   integer(i4b),intent(in) :: ixRichards ! index defining the option for Richards' equation (moisture or mixdform)
           integer(i4b),intent(in) :: ixInfRateMax ! index defining the maximum infiltration rate method (GreenAmpt, topmodel_GA, noInfiltrationExcess)
   integer(i4b),intent(in) :: nRoots          ! number of layers that contain roots
      integer(i4b),intent(in) :: nSoil           ! number of soil layers
   integer(i4b),intent(in) :: ixIce           ! index of lowest ice layer
   integer(i4b),intent(in) :: nSnow,nLayers
   ! input: parameters
   real(rkind),intent(in) :: FUSE_Ac_max, FUSE_phi_tens
   real(rkind),intent(in) :: FUSE_b
   real(rkind),intent(in) :: FUSE_lambda
   real(rkind),intent(in) :: FUSE_chi
   real(rkind),intent(in) :: FUSE_mu
   real(rkind),intent(in) :: FUSE_n
   ! input: rain plus melt
   real(rkind),intent(in) :: scalarRainPlusMelt ! rain plus melt  (m s-1)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: scalarMatricHeadLiq  ! liquid matric head in the upper-most soil layer (m)
   real(rkind),intent(in) :: scalarVolFracLiq     ! volumetric liquid water content in the upper-most soil layer (-)
   real(rkind),intent(in) :: mLayerVolFracLiq(:)  ! volumetric liquid water content in each soil layer (-)
   real(rkind),intent(in) :: mLayerTemp(:)          ! temperature (K)
   real(rkind),intent(in) :: mLayerMatricHead(:)    ! matric head in each soil layer (m)
      real(rkind),intent(in) :: mLayerVolFracIce(:)     ! volumetric ice content in each soil layer (-)

   ! input: depth of upper-most soil layer (m)
   real(rkind),intent(in) :: mLayerDepth(:)   ! depth of upper-most soil layer (m)
     real(rkind),intent(in) :: iLayerHeight(0:)      ! height at the interface of each layer for soil layers only (m)

   ! input: diriclet boundary conditions
   real(rkind),intent(in) :: upperBoundHead    ! upper boundary condition for matric head (m)
   real(rkind),intent(in) :: upperBoundTheta   ! upper boundary condition for volumetric liquid water content (-)
   ! input: transmittance
   real(rkind),intent(in) :: surfaceSatHydCond  ! saturated hydraulic conductivity at the surface (m s-1)
   real(rkind),intent(in) :: dHydCond_dTemp     ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   real(rkind),intent(in) :: iceImpedeFac       ! ice impedence factor in the upper-most soil layer (-)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha            ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n                ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m                ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: theta_sat            ! soil porosity (-)
   real(rkind),intent(in) :: theta_res            ! soil residual volumetric water content (-)
      real(rkind),intent(in) :: rootingDepth         ! rooting depth (m)
     real(rkind),intent(in) :: zScale_TOPMODEL      ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
  real(rkind),intent(in) :: wettingFrontSuction  ! Green-Ampt wetting front suction (m)
     real(rkind),intent(in) :: qSurfScale ! scaling factor in the surface runoff parameterization (-)
   real(rkind),intent(in) :: soilIceScale         ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
   real(rkind),intent(in) :: soilIceCV            ! soil ice CV in Gamma distribution used to define frozen area (-)

     ! input: pre-computed derivatives in ...
  real(rkind),intent(in) :: dTheta_dTk(:)          ! ... volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in) :: dTheta_dPsi(:)         ! ... the soil water characteristic w.r.t. psi (m-1)
   real(rkind),intent(in) :: mLayerdPsi_dTheta(:)      ! ... the soil water characteristic w.r.t. theta (m)

   ! input-output: hydraulic conductivity and diffusivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   real(rkind),intent(inout) :: surfaceHydCond  ! hydraulic conductivity (m s-1)
   real(rkind),intent(inout) :: surfaceDiffuse  ! hydraulic diffusivity at the surface (m2 s-1)
   ! output: runoff and infiltration
   real(rkind),intent(inout) :: scalarSurfaceRunoff        ! surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceRunoff_IE     ! infiltration excess surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceRunoff_SE     ! saturation excess surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceInfiltration  ! surface infiltration (m s-1)
   ! input-output
   real(rkind),intent(inout) :: xMaxInfilRate, scalarInfilArea, scalarFrozenArea
      ! output: derivatives in surface infiltration w.r.t. ...
   real(rkind),intent(inout) :: scalarSoilControl   ! soil control on infiltration for derivative
   real(rkind),intent(inout) :: dq_dHydStateVec(0:)     ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateVec(0:)     ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
    ! -----------------------------------------------------------------------------------------------------------------------------
    ! local variables
  ! local variables
  ! general
  ! integer(i4b)                     :: iLayer                              ! index of soil layer
  real(rkind)                      :: Tcrit                               ! temperature where all water is unfrozen (K)
  real(rkind)                      :: fPart1,fPart2                       ! different parts of a function
  real(rkind)                      :: dPart1(1:nSoil)     ! derivatives for different parts of a function
  ! integer(i4b),parameter :: dPart1 = 0
  ! integer(i4b),parameter :: dPart2 = 1
  real(rkind)                      :: dPart2(1:nSoil)     ! derivatives for different parts of a function
  real(rkind)                      :: dfracCap(1:nSoil)   ! derivatives for different parts of a function
  ! integer(i4b),parameter :: dfracCap = 2
  real(rkind)                      :: dfInfRaw(1:nSoil)   ! derivatives for different parts of a function
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
  real(rkind) :: dVolFracLiq_dWat(1:nSoil)  ! ... vol fraction of liquid w.r.t. water state variable in root layers
  real(rkind) :: dVolFracIce_dWat(1:nSoil)  ! ... vol fraction of ice w.r.t. water state variable in root layers
  real(rkind) :: dVolFracLiq_dTk(1:nSoil)   ! ... vol fraction of liquid w.r.t. temperature in root layers
  real(rkind) :: dVolFracIce_dTk(1:nSoil)   ! ... vol fraction of ice w.r.t. temperature in root layers
  real(rkind) :: dRootZoneLiq_dWat(1:nSoil) ! ... vol fraction of scalar root zone liquid w.r.t. water state variable in root layers
  real(rkind) :: dRootZoneIce_dWat(1:nSoil) ! ... vol fraction of scalar root zone ice w.r.t. water state variable in root layers
  real(rkind) :: dRootZoneLiq_dTk(1:nSoil)  ! ... vol fraction of scalar root zone liquid w.r.t. temperature in root layers
  real(rkind) :: dRootZoneIce_dTk(1:nSoil)  ! ... vol fraction of scalar root zone ice w.r.t. temperature in root layers
  real(rkind) :: dDepthWettingFront_dWat(1:nSoil) ! ... scalar depth of wetting front w.r.t. water state variable in root layers
  real(rkind) :: dDepthWettingFront_dTk(1:nSoil)  ! ... scalar depth of wetting front w.r.t. temperature in root layers
  real(rkind) :: dxMaxInfilRate_dWat(1:nSoil) ! ... scalar max infiltration rate w.r.t. water state variable in root layers
  real(rkind) :: dxMaxInfilRate_dTk(1:nSoil)  ! ... scalar max infiltration rate w.r.t. temperature in root layers
  real(rkind) :: dInfilArea_dWat(1:nSoil)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
  real(rkind) :: dInfilArea_dTk(1:nSoil)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  real(rkind) :: dFrozenArea_dWat(1:nSoil) ! ... scalar frozen area w.r.t. water state variable in canopy or snow and root layers
  real(rkind) :: dFrozenArea_dTk(1:nSoil)  ! ... scalar frozen area w.r.t. temperature in canopy or snow and root layers
  real(rkind) :: dInfilRate_dWat(1:nSoil)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
  real(rkind) :: dInfilRate_dTk(1:nSoil)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  ! component variables for infiltration excess (IE) and saturation excess (SE) surface runoff
  real(rkind) :: SR_IE ! infiltration excess surface runoff component
  real(rkind) :: SR_SE ! saturation excess surface runoff component
  real(rkind) :: dq_dHydStateVec_IE(0:nSoil) ! derivative of infiltration w.r.t hydrology state variable (infiltration excess component)
  ! integer(i4b),parameter :: dq_dHydStateVec_IE=0
  real(rkind) :: dq_dHydStateVec_SE(0:nSoil) ! derivative of infiltration w.r.t hydrology state variable (saturation excess component)
  ! integer(i4b),parameter :: dq_dHydStateVec_SE=1
  real(rkind) :: dq_dNrgStateVec_IE(0:nSoil) ! derivative of infiltration w.r.t energy state variable (infiltration excess component)
  ! integer(i4b),parameter :: dq_dNrgStateVec_IE=2
  real(rkind) :: dq_dNrgStateVec_SE(0:nSoil) ! derivative of infiltration w.r.t energy state variable (saturation excess component)
  ! integer(i4b),parameter :: dq_dNrgStateVec_SE=3
  
    call initialize_surfaceFlx(&
      dq_dHydStateVec, &
      dq_dNrgStateVec, &
      dq_dHydStateVec_IE, &
      dq_dHydStateVec_SE, &
      dq_dNrgStateVec_IE, &
      dq_dNrgStateVec_SE, &
      scalarSurfaceRunoff, &
      scalarSurfaceRunoff_IE, &
      scalarSurfaceRunoff_SE, &
      scalarSurfaceInfiltration)
      print*, firstSplitOper, bc_upper, surfRun_IE, surfRun_SE
      print*, nSnow+1, nLayers
!     call update_surfaceFlx(firstSplitOper, &
!       bc_upper, &
!       surfRun_IE, &
!       surfRun_SE, &
!       ixRichards, ixInfRateMax, &
!       nRoots, nSoil, ixIce, &
!       FUSE_Ac_max, FUSE_phi_tens, FUSE_b, FUSE_lambda, FUSE_chi, FUSE_mu, FUSE_n, &
!       scalarRainPlusMelt, &
!       scalarMatricHeadLiq, scalarVolFracLiq, mLayerVolFracLiq(nSnow+1:nLayers), mLayerTemp(nSnow+1:nLayers), mLayerMatricHead, mLayerVolFracIce(nSnow+1:nLayers), &
!       mLayerDepth(nSnow+1:nLayers),iLayerHeight, &
!       upperBoundHead, upperBoundTheta, &
!       surfaceSatHydCond, dHydCond_dTemp, iceImpedeFac, &
!       vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, &
!       rootingDepth, zScale_TOPMODEL, wettingFrontSuction, qSurfScale, soilIceScale, soilIceCV, &
!       dTheta_dTk(nSnow+1:nLayers), dTheta_dPsi, mLayerdPsi_dTheta, &
!       surfaceHydCond, surfaceDiffuse, &
!       scalarSurfaceRunoff, scalarSurfaceRunoff_IE, scalarSurfaceRunoff_SE, scalarSurfaceInfiltration, &
!       xMaxInfilRate, scalarInfilArea, scalarFrozenArea, &
!       scalarSoilControl, dq_dHydStateVec, dq_dNrgStateVec, &
!       Tcrit, fPart1, fPart2, dPart1, dPart2, dfracCap, dfInfRaw, &
!  total_soil_depth, cFlux, rootZoneLiq, rootZoneIce, availCapacity, depthWettingFront, hydCondWettingFront, &
!  fracCap, fInfRaw, alpha, xLimg, &
!   dVolFracLiq_dWat, dVolFracIce_dWat, dVolFracLiq_dTk, dVolFracIce_dTk, &
!   dRootZoneLiq_dWat, dRootZoneIce_dWat, dRootZoneLiq_dTk, dRootZoneIce_dTk, &
!    dDepthWettingFront_dWat, dDepthWettingFront_dTk, &
!  dxMaxInfilRate_dWat, dxMaxInfilRate_dTk, &
!  dInfilArea_dWat, dInfilArea_dTk, dFrozenArea_dWat, dFrozenArea_dTk, &
!  dInfilRate_dWat, dInfilRate_dTk, &
! SR_IE, SR_SE, dq_dHydStateVec_IE, dq_dHydStateVec_SE, dq_dNrgStateVec_IE, dq_dNrgStateVec_SE)
    ! call finalize_surfaceFlx
    
  end subroutine surfaceFlx

  attributes(device) subroutine initialize_surfaceFlx(&
      dq_dHydStateVec, &
      dq_dNrgStateVec, &
      dq_dHydStateVec_IE, &
      dq_dHydStateVec_SE, &
      dq_dNrgStateVec_IE, &
      dq_dNrgStateVec_SE, &
      scalarSurfaceRunoff, &
      scalarSurfaceRunoff_IE, &
      scalarSurfaceRunoff_SE, &
      scalarSurfaceInfiltration)
    real(rkind),intent(inout) :: dq_dHydStateVec(0:)
    real(rkind),intent(inout) :: dq_dNrgStateVec(0:)
    real(rkind),intent(inout) :: dq_dHydStateVec_IE(0:)
    real(rkind),intent(inout) :: dq_dHydStateVec_SE(0:)
    real(rkind),intent(inout) :: dq_dNrgStateVec_IE(0:)
    real(rkind),intent(inout) :: dq_dNrgStateVec_SE(0:)
    real(rkind),intent(inout) :: scalarSurfaceRunoff
    real(rkind),intent(inout) :: scalarSurfaceRunoff_IE
    real(rkind),intent(inout) :: scalarSurfaceRunoff_SE
    real(rkind),intent(inout) :: scalarSurfaceInfiltration

    integer(i4b) :: iLayer
    ! **** Initialize operations for surfaceFlx ****
 
    ! initialize derivatives
    do iLayer=0,size(dq_dHydStateVec)-1
      dq_dHydStateVec(iLayer) = 0._rkind
      dq_dNrgStateVec(iLayer) = 0._rkind ! energy state variable is temperature (transformed outside soilLiqFlx_module if needed)
    end do

    ! allocate and initialize (to zero) surface runoff component arrays for infiltration derivatives ...
    do iLayer=0,size(dq_dHydStateVec)-1
      dq_dHydStateVec_IE(iLayer) = dq_dHydStateVec(iLayer) ! ... w.r.t hydrology state variable (infiltration excess component)  
      dq_dHydStateVec_SE(iLayer) = dq_dHydStateVec(iLayer) ! ... w.r.t hydrology state variable (saturation excess component)
      dq_dNrgStateVec_IE(iLayer) = dq_dNrgStateVec(iLayer) ! ... w.r.t energy state variable (infiltration excess component)
      dq_dNrgStateVec_SE(iLayer) = dq_dNrgStateVec(iLayer) ! ... w.r.t energy state variable (saturation excess component)
    end do

  ! initialize runoff and infiltration values
   scalarSurfaceRunoff       = 0._rkind 
   scalarSurfaceRunoff_IE    = 0._rkind  
   scalarSurfaceRunoff_SE    = 0._rkind  
   scalarSurfaceInfiltration = 0._rkind 
  end subroutine initialize_surfaceFlx

  attributes(device) subroutine update_surfaceFlx(firstSplitOper, &
      bc_upper, &
      surfRun_IE, &
      surfRun_SE, &
      ixRichards, ixInfRateMax, &
      nRoots, nSoil, ixIce, &
      FUSE_Ac_max, FUSE_phi_tens, FUSE_b, FUSE_lambda, FUSE_chi, FUSE_mu, FUSE_n, &
      scalarRainPlusMelt, &
      scalarMatricHeadLiq, scalarVolFracLiq, mLayerVolFracLiq, mLayerTemp, mLayerMatricHead, mLayerVolFracIce, &
      mLayerDepth,iLayerHeight, &
      upperBoundHead, upperBoundTheta, &
      surfaceSatHydCond, dHydCond_dTemp, iceImpedeFac, &
      vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, &
      rootingDepth, zScale_TOPMODEL, wettingFrontSuction, qSurfScale, soilIceScale, soilIceCV, &
      dTheta_dTk, dTheta_dPsi, mLayerdPsi_dTheta, &
      surfaceHydCond, surfaceDiffuse, &
      scalarSurfaceRunoff, scalarSurfaceRunoff_IE, scalarSurfaceRunoff_SE, scalarSurfaceInfiltration, &
      xMaxInfilRate, scalarInfilArea, scalarFrozenArea, &
      scalarSoilControl, dq_dHydStateVec, dq_dNrgStateVec, &
      Tcrit, fPart1, fPart2, dPart1, dPart2, dfracCap, dfInfRaw, &
 total_soil_depth, cFlux, rootZoneLiq, rootZoneIce, availCapacity, depthWettingFront, hydCondWettingFront, &
 fracCap, fInfRaw, alpha, xLimg, &
  dVolFracLiq_dWat, dVolFracIce_dWat, dVolFracLiq_dTk, dVolFracIce_dTk, &
  dRootZoneLiq_dWat, dRootZoneIce_dWat, dRootZoneLiq_dTk, dRootZoneIce_dTk, &
   dDepthWettingFront_dWat, dDepthWettingFront_dTk, &
 dxMaxInfilRate_dWat, dxMaxInfilRate_dTk, &
 dInfilArea_dWat, dInfilArea_dTk, dFrozenArea_dWat, dFrozenArea_dTk, &
 dInfilRate_dWat, dInfilRate_dTk, &
SR_IE, SR_SE, dq_dHydStateVec_IE, dq_dHydStateVec_SE, dq_dNrgStateVec_IE, dq_dNrgStateVec_SE)
    ! input: model control
    logical(lgt),intent(in) :: firstSplitOper ! flag indicating if desire to compute infiltration
    integer(i4b),intent(in) :: bc_upper ! index defining the type of boundary conditions
    integer(i4b),intent(in) :: surfRun_IE ! index defining the infiltration excess surface runoff method
    integer(i4b),intent(in) :: surfRun_SE ! index defining the saturation excess surface runoff method
   integer(i4b),intent(in) :: ixRichards ! index defining the option for Richards' equation (moisture or mixdform)
           integer(i4b),intent(in) :: ixInfRateMax ! index defining the maximum infiltration rate method (GreenAmpt, topmodel_GA, noInfiltrationExcess)
   integer(i4b),intent(in) :: nRoots          ! number of layers that contain roots
      integer(i4b),intent(in) :: nSoil           ! number of soil layers
   integer(i4b),intent(in) :: ixIce           ! index of lowest ice layer

   ! input: parameters
   real(rkind),intent(in) :: FUSE_Ac_max, FUSE_phi_tens
   real(rkind),intent(in) :: FUSE_b
   real(rkind),intent(in) :: FUSE_lambda
   real(rkind),intent(in) :: FUSE_chi
   real(rkind),intent(in) :: FUSE_mu
   real(rkind),intent(in) :: FUSE_n
   ! input: rain plus melt
   real(rkind),intent(in) :: scalarRainPlusMelt ! rain plus melt  (m s-1)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: scalarMatricHeadLiq  ! liquid matric head in the upper-most soil layer (m)
   real(rkind),intent(in) :: scalarVolFracLiq     ! volumetric liquid water content in the upper-most soil layer (-)
   real(rkind),intent(in) :: mLayerVolFracLiq(:)  ! volumetric liquid water content in each soil layer (-)
   real(rkind),intent(in) :: mLayerTemp(:)          ! temperature (K)
   real(rkind),intent(in) :: mLayerMatricHead(:)    ! matric head in each soil layer (m)
      real(rkind),intent(in) :: mLayerVolFracIce(:)     ! volumetric ice content in each soil layer (-)

   ! input: depth of upper-most soil layer (m)
   real(rkind),intent(in) :: mLayerDepth(:)   ! depth of upper-most soil layer (m)
     real(rkind),intent(in) :: iLayerHeight(0:)      ! height at the interface of each layer for soil layers only (m)

   ! input: diriclet boundary conditions
   real(rkind),intent(in) :: upperBoundHead    ! upper boundary condition for matric head (m)
   real(rkind),intent(in) :: upperBoundTheta   ! upper boundary condition for volumetric liquid water content (-)
   ! input: transmittance
   real(rkind),intent(in) :: surfaceSatHydCond  ! saturated hydraulic conductivity at the surface (m s-1)
   real(rkind),intent(in) :: dHydCond_dTemp     ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   real(rkind),intent(in) :: iceImpedeFac       ! ice impedence factor in the upper-most soil layer (-)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha            ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n                ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m                ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: theta_sat            ! soil porosity (-)
   real(rkind),intent(in) :: theta_res            ! soil residual volumetric water content (-)
      real(rkind),intent(in) :: rootingDepth         ! rooting depth (m)
     real(rkind),intent(in) :: zScale_TOPMODEL      ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
  real(rkind),intent(in) :: wettingFrontSuction  ! Green-Ampt wetting front suction (m)
     real(rkind),intent(in) :: qSurfScale ! scaling factor in the surface runoff parameterization (-)
   real(rkind),intent(in) :: soilIceScale         ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
   real(rkind),intent(in) :: soilIceCV            ! soil ice CV in Gamma distribution used to define frozen area (-)

     ! input: pre-computed derivatives in ...
  real(rkind),intent(in) :: dTheta_dTk(:)          ! ... volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in) :: dTheta_dPsi(:)         ! ... the soil water characteristic w.r.t. psi (m-1)
   real(rkind),intent(in) :: mLayerdPsi_dTheta(:)      ! ... the soil water characteristic w.r.t. theta (m)

   ! input-output: hydraulic conductivity and diffusivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   real(rkind),intent(inout) :: surfaceHydCond  ! hydraulic conductivity (m s-1)
   real(rkind),intent(inout) :: surfaceDiffuse  ! hydraulic diffusivity at the surface (m2 s-1)
   ! output: runoff and infiltration
   real(rkind),intent(inout) :: scalarSurfaceRunoff        ! surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceRunoff_IE     ! infiltration excess surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceRunoff_SE     ! saturation excess surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceInfiltration  ! surface infiltration (m s-1)
   ! input-output
   real(rkind),intent(inout) :: xMaxInfilRate, scalarInfilArea, scalarFrozenArea
      ! output: derivatives in surface infiltration w.r.t. ...
   real(rkind),intent(inout) :: scalarSoilControl   ! soil control on infiltration for derivative
   real(rkind),intent(inout) :: dq_dHydStateVec(0:)     ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateVec(0:)     ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)

   ! inout: local to surfaceFlux
     real(rkind),intent(inout)                      :: Tcrit                               ! temperature where all water is unfrozen (K)
      real(rkind),intent(inout)                      :: fPart1,fPart2                       ! different parts of a function
    real(rkind),intent(inout)                      :: dPart1(:)     ! derivatives for different parts of a function
    real(rkind),intent(inout)                      :: dPart2(:)     ! derivatives for different parts of a function
       real(rkind),intent(inout)                      :: dfracCap(:)   ! derivatives for different parts of a function
   real(rkind),intent(inout)                      :: dfInfRaw(:)   ! derivatives for different parts of a function
    real(rkind),intent(inout)                      :: total_soil_depth                    ! total depth of soil (m)
       real(rkind),intent(inout)                      :: cFlux                               ! capillary flux (m s-1)

    real(rkind),intent(inout)                      :: rootZoneLiq                         ! depth of liquid water in the root zone (m)
      real(rkind),intent(inout)                      :: rootZoneIce                         ! depth of ice in the root zone (m)

    real(rkind),intent(inout)                      :: availCapacity                       ! available storage capacity in the root zone (m)
    real(rkind),intent(inout)                      :: depthWettingFront                   ! depth to the wetting front (m)
    real(rkind),intent(inout)                      :: hydCondWettingFront                 ! hydraulic conductivity at the wetting front (m s-1)
   real(rkind),intent(inout)                      :: fracCap                             ! fraction of pore space filled with liquid water and ice (-)
   real(rkind),intent(inout)                      :: fInfRaw                             ! infiltrating area before imposing solution constraints (-)
      real(rkind),intent(inout)                      :: alpha                               ! shape parameter in the Gamma distribution
   real(rkind),intent(inout)                      :: xLimg                               ! upper limit of the integral

   real(rkind),intent(inout) :: dVolFracLiq_dWat(:)  ! ... vol fraction of liquid w.r.t. water state variable in root layers
   real(rkind),intent(inout) :: dVolFracIce_dWat(:)  ! ... vol fraction of ice w.r.t. water state variable in root layers
  real(rkind),intent(inout) :: dVolFracLiq_dTk(:)   ! ... vol fraction of liquid w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dVolFracIce_dTk(:)   ! ... vol fraction of ice w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dRootZoneLiq_dWat(:) ! ... vol fraction of scalar root zone liquid w.r.t. water state variable in root layers
  real(rkind),intent(inout) :: dRootZoneLiq_dTk(:)  ! ... vol fraction of scalar root zone liquid w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dRootZoneIce_dTk(:)  ! ... vol fraction of scalar root zone ice w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dRootZoneIce_dWat(:) ! ... vol fraction of scalar root zone ice w.r.t. water state variable in root layers
      real(rkind),intent(inout) :: dDepthWettingFront_dWat(:) ! ... scalar depth of wetting front w.r.t. water state variable in root layers
    real(rkind),intent(inout) :: dDepthWettingFront_dTk(:)  ! ... scalar depth of wetting front w.r.t. temperature in root layers
    real(rkind),intent(inout) :: dxMaxInfilRate_dWat(:) ! ... scalar max infiltration rate w.r.t. water state variable in root layers
    real(rkind),intent(inout) :: dxMaxInfilRate_dTk(:)  ! ... scalar max infiltration rate w.r.t. temperature in root layers
   real(rkind),intent(inout) :: dInfilArea_dWat(:)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
   real(rkind),intent(inout) :: dInfilArea_dTk(:)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  real(rkind),intent(inout) :: dFrozenArea_dWat(:) ! ... scalar frozen area w.r.t. water state variable in canopy or snow and root layers
  real(rkind),intent(inout) :: dFrozenArea_dTk(:)  ! ... scalar frozen area w.r.t. temperature in canopy or snow and root layers
  real(rkind),intent(inout) :: dInfilRate_dWat(:)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
  real(rkind),intent(inout) :: dInfilRate_dTk(:)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
     real(rkind),intent(inout) :: SR_IE, SR_SE
  real(rkind),intent(inout) :: dq_dHydStateVec_IE(0:) ! derivative of infiltration w.r.t hydrology state variable (infiltration excess component)
  real(rkind),intent(inout) :: dq_dHydStateVec_SE(0:) ! derivative of infiltration w.r.t hydrology state variable (saturation excess component)
  real(rkind),intent(inout) :: dq_dNrgStateVec_IE(0:) ! derivative of infiltration w.r.t energy state variable (infiltration excess component)
  real(rkind),intent(inout) :: dq_dNrgStateVec_SE(0:) ! derivative of infiltration w.r.t energy state variable (saturation excess component)


    ! **** Update operations for surfaceFlx ****

   ! output: derivatives in surface infiltration w.r.t. ...
  !  dq_dHydStateVec => out_surfaceFlx % dq_dHydStateVec, & ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
  !  dq_dNrgStateVec => out_surfaceFlx % dq_dNrgStateVec, & ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
  !  ! output: error control
  !  err      => out_surfaceFlx % err    , & ! error code
  !  message  => out_surfaceFlx % message  & ! error message
  ! &)
 
  !  ! compute the surface flux and its derivative
   if (firstSplitOper .or. updateInfil) then
     select case(bc_upper)
       case(prescribedHead) ! head condition
         call update_surfaceFlx_prescribedHead(ixRichards, &
    scalarMatricHeadLiq, scalarVolFracLiq, &
    mLayerDepth, &
    upperBoundHead, upperBoundTheta, &
    surfaceSatHydCond, dHydCond_dTemp, iceImpedeFac, &
    vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, &
    surfaceHydCond, surfaceDiffuse, &
    scalarSurfaceRunoff, scalarSurfaceRunoff_IE, scalarSurfaceRunoff_SE, scalarSurfaceInfiltration, &
    scalarSoilControl, dq_dHydStateVec, dq_dNrgStateVec, &
    xMaxInfilRate, scalarInfilArea, scalarFrozenArea, &
    cFlux)
        !  if (return_flag) return 
 
       case(liquidFlux)     ! flux condition

         select case(surfRun_IE) ! infiltration excess surface runoff
           case(zero_IE)         ! zero infiltration excess surface runoff
             call update_surfaceFlx_zero_IE(SR_IE, dq_dHydStateVec_IE, dq_dNrgStateVec_IE);        
             ! if (return_flag) return 

           case(homegrown_IE)    ! homegrown infiltration excess surface runoff
             call update_surfaceFlx_liquidFlux(ixRichards, ixInfRateMax, surfRun_IE, surfRun_SE, nRoots, nSoil, ixIce, &
 scalarRainPlusMelt, &
  mLayerTemp, mLayerMatricHead, mLayerVolFracLiq, mLayerVolFracIce, &
  dTheta_dTk, dTheta_dPsi, mLayerdPsi_dTheta, &
  mLayerDepth, iLayerHeight, &
  surfaceSatHydCond, &
  theta_sat, rootingDepth, zScale_TOPMODEL, wettingFrontSuction, qSurfScale, soilIceScale, soilIceCV, &
  xMaxInfilRate, scalarInfilArea, scalarSoilControl, scalarFrozenArea, &
  surfaceHydCond, surfaceDiffuse, &
  dq_dHydStateVec, dq_dNrgStateVec, &
  Tcrit, fPart1, fPart2, dPart1, dPart2, dfracCap, dfInfRaw, &
 total_soil_depth, rootZoneLiq, rootZoneIce, availCapacity, depthWettingFront, hydCondWettingFront, &
 fracCap, fInfRaw, alpha, xLimg, &
  dVolFracLiq_dWat, dVolFracIce_dWat, dVolFracLiq_dTk, dVolFracIce_dTk, &
  dRootZoneLiq_dWat, dRootZoneIce_dWat, dRootZoneLiq_dTk, dRootZoneIce_dTk, &
   dDepthWettingFront_dWat, dDepthWettingFront_dTk, &
 dxMaxInfilRate_dWat, dxMaxInfilRate_dTk, &
 dInfilArea_dWat, dInfilArea_dTk, dFrozenArea_dWat, dFrozenArea_dTk, &
 dInfilRate_dWat, dInfilRate_dTk, &
SR_IE, SR_SE, dq_dHydStateVec_IE, dq_dHydStateVec_SE, dq_dNrgStateVec_IE, dq_dNrgStateVec_SE);     
             ! if (return_flag) return 

           case default
            !  err=20; message=trim(message)//'unknown infiltration excess surface runoff method';
            !  return_flag=.true.; 
            return
         end select

         select case(surfRun_SE) ! saturation excess surface runoff
           case(zero_SE)         ! zero saturation excess surface runoff
             call update_surfaceFlx_zero_SE(SR_SE, dq_dHydStateVec_SE, dq_dNrgStateVec_SE);        
             ! if (return_flag) return 

           case(homegrown_SE)    ! homegrown saturation excess surface runoff
             if (surfRun_IE /= homegrown_IE) then ! avoid repeating computations
              call update_surfaceFlx_liquidFlux(ixRichards, ixInfRateMax, surfRun_IE, surfRun_SE, nRoots, nSoil, ixIce, &
 scalarRainPlusMelt, &
  mLayerTemp, mLayerMatricHead, mLayerVolFracLiq, mLayerVolFracIce, &
  dTheta_dTk, dTheta_dPsi, mLayerdPsi_dTheta, &
  mLayerDepth, iLayerHeight, &
  surfaceSatHydCond, &
  theta_sat, rootingDepth, zScale_TOPMODEL, wettingFrontSuction, qSurfScale, soilIceScale, soilIceCV, &
  xMaxInfilRate, scalarInfilArea, scalarSoilControl, scalarFrozenArea, &
  surfaceHydCond, surfaceDiffuse, &
  dq_dHydStateVec, dq_dNrgStateVec, &
  Tcrit, fPart1, fPart2, dPart1, dPart2, dfracCap, dfInfRaw, &
 total_soil_depth, rootZoneLiq, rootZoneIce, availCapacity, depthWettingFront, hydCondWettingFront, &
 fracCap, fInfRaw, alpha, xLimg, &
  dVolFracLiq_dWat, dVolFracIce_dWat, dVolFracLiq_dTk, dVolFracIce_dTk, &
  dRootZoneLiq_dWat, dRootZoneIce_dWat, dRootZoneLiq_dTk, dRootZoneIce_dTk, &
   dDepthWettingFront_dWat, dDepthWettingFront_dTk, &
 dxMaxInfilRate_dWat, dxMaxInfilRate_dTk, &
 dInfilArea_dWat, dInfilArea_dTk, dFrozenArea_dWat, dFrozenArea_dTk, &
 dInfilRate_dWat, dInfilRate_dTk, &
SR_IE, SR_SE, dq_dHydStateVec_IE, dq_dHydStateVec_SE, dq_dNrgStateVec_IE, dq_dNrgStateVec_SE);     
              ! if (return_flag) return 
             end if

           case(FUSEPRMS)        ! FUSE PRMS surface runoff
             call update_surfaceFlx_FUSE_PRMS(ixRichards, &
 FUSE_Ac_max, FUSE_phi_tens, theta_sat, &
 scalarRainPlusMelt, &
 mLayerVolFracLiq, mLayerTemp, mLayerMatricHead, &
 dTheta_dTk, dTheta_dPsi, &
 nSoil, &
 mLayerDepth, iLayerHeight, &
 dVolFracLiq_dWat, dVolFracLiq_dTk, dq_dHydStateVec_SE, dq_dNrgStateVec_SE, SR_SE, Tcrit);      
            !  if (return_flag) return 

           case(FUSEAVIC)        ! FUSE ARNO/VIC surface runoff
             call update_surfaceFlx_FUSE_ARNO_VIC(ixRichards,&
 FUSE_b, theta_sat, &
 scalarRainPlusMelt, &
 mLayerTemp, mLayerMatricHead, mLayerVolFracLiq, &
 dTheta_dTk, dTheta_dPsi, &
 nSoil, mLayerDepth, iLayerHeight, &
 dVolFracLiq_dWat, dVolFracLiq_dTk, SR_SE, dq_dHydStateVec_SE, dq_dNrgStateVec_SE, Tcrit);  
            !  if (return_flag) return

           case(FUSETOPM)        ! FUSE TOPMODEL surface runoff
             call update_surfaceFlx_FUSE_TOPMODEL(ixRichards, &
    FUSE_lambda, FUSE_chi, FUSE_mu, FUSE_n, theta_sat, &
    scalarRainPlusMelt, &
    mLayerTemp, mLayerMatricHead, mLayerVolFracLiq, &
    dTheta_dTk, dTheta_dPsi, &
    nSoil, mLayerDepth, iLayerHeight, &
    dVolFracLiq_dWat, dVolFracLiq_dTk, SR_SE, &
    dq_dHydStateVec_SE, dq_dNrgStateVec_SE, Tcrit);  
            !  if (return_flag) return
           case default
            !  err=20; message=trim(message)//'unknown saturation excess surface runoff method';
            !  return_flag=.true.; 
            return
         end select

         call update_gather_runoff_components(scalarRainPlusMelt,&
      scalarSurfaceRunoff_IE, &
      scalarSurfaceRunoff_SE, &
      scalarSurfaceRunoff, &
      scalarSurfaceInfiltration, &
      scalarSoilControl, &
      dq_dHydStateVec, &
      dq_dNrgStateVec, &
      dq_dHydStateVec_IE, dq_dHydStateVec_SE, dq_dNrgStateVec_IE, dq_dNrgStateVec_SE, SR_IE, SR_SE);  
        !  if (return_flag) return 

       case default; 
        ! err=20; message=trim(message)//'unknown upper boundary condition for soil hydrology'; return_flag=.true.; 
        return
 
     end select 
   else ! do not compute infiltration after first flux call in a splitting operation unless updateInfil is true
     dq_dHydStateVec(:) = 0._rkind
     dq_dNrgStateVec(:) = 0._rkind ! energy state variable is temperature (transformed outside soilLiqFlx_module if needed)
   end if 

  ! end associate

  end subroutine update_surfaceFlx  

  attributes(device) subroutine update_gather_runoff_components(scalarRainPlusMelt,&
      scalarSurfaceRunoff_IE, &
      scalarSurfaceRunoff_SE, &
      scalarSurfaceRunoff, &
      scalarSurfaceInfiltration, &
      scalarSoilControl, &
      dq_dHydStateVec, &
      dq_dNrgStateVec, &
      dq_dHydStateVec_IE, dq_dHydStateVec_SE, dq_dNrgStateVec_IE, dq_dNrgStateVec_SE, SR_IE, SR_SE)
      implicit none
    real(rkind),intent(in) :: scalarRainPlusMelt ! rain plus melt  (m s-1)
    real(rkind),intent(inout) :: scalarSurfaceRunoff_IE ! infiltration excess surface runoff (m s-1)
    real(rkind),intent(inout) :: scalarSurfaceRunoff_SE ! saturation excess surface runoff (m s-1)
    real(rkind),intent(inout) :: scalarSurfaceRunoff ! surface runoff (m s-1)
    real(rkind),intent(inout) :: scalarSurfaceInfiltration ! surface infiltration (m s-1)
    real(rkind),intent(inout) :: scalarSoilControl ! soil control on infiltration for derivative
    real(rkind),intent(inout) :: dq_dHydStateVec(0:) ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
    real(rkind),intent(inout) :: dq_dNrgStateVec(0:) ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
    real(rkind),intent(in) :: dq_dHydStateVec_IE(0:), dq_dHydStateVec_SE(0:), dq_dNrgStateVec_IE(0:), dq_dNrgStateVec_SE(0:)
    real(rkind),intent(inout) :: SR_IE, SR_SE

    integer(i4b) :: iLayer

  ! **** Gather surface runoff components for the liquid flux upper hydrology boundary condition ****
  real(rkind) :: roundoff_tolerance   ! tolerance for round-off error

   ! validate against rain plus melt
   if (SR_IE > scalarRainPlusMelt+roundoff_tolerance) then ! runoff above RPM outside tolerance
    ! err=20; message=trim(message)//'infiltration excess surface runoff greater than rain plus melt'; return_flag=.true.; 
    return
   else if (SR_IE > scalarRainPlusMelt) then               ! runoff slightly above RPM but within tolerance
    SR_IE = scalarRainPlusMelt
   end if
   if (SR_SE > scalarRainPlusMelt+roundoff_tolerance) then ! runoff above RPM outside tolerance
    ! err=20; message=trim(message)//'saturation excess surface runoff greater than rain plus melt'; return_flag=.true.; 
    return
   else if (SR_SE > scalarRainPlusMelt) then               ! runoff slightly above RPM but within tolerance
    SR_SE = scalarRainPlusMelt
   end if

   ! validate against zero
   if (SR_IE < (-roundoff_tolerance)) then ! runoff below zero outside tolerance
    ! err=20; message=trim(message)//'infiltration excess surface runoff below zero'; return_flag=.true.; return
   else if (SR_IE < 0._rkind) then      ! runoff slightly below zero but within tolerance
    SR_IE = 0._rkind
   end if
   if (SR_SE < (-roundoff_tolerance)) then ! runoff below zero outside tolerance
    ! err=20; message=trim(message)//'saturation excess surface runoff below zero'; return_flag=.true.; return
   else if (SR_SE < 0._rkind) then      ! runoff slightly below zero but within tolerance
    SR_SE = 0._rkind
   end if

   ! validate against zero
   if (SR_IE < (-roundoff_tolerance)) then ! runoff below zero outside tolerance
    !err=20; message=trim(message)//'infiltration excess surface runoff below zero'; return_flag=.true.; 
    return
   else if (SR_IE < 0._rkind) then      ! runoff slightly below zero but within tolerance
    SR_IE = 0._rkind
   end if
   if (SR_SE < (-roundoff_tolerance)) then ! runoff below zero outside tolerance
    !err=20; message=trim(message)//'saturation excess surface runoff below zero'; return_flag=.true.; 
    return
   else if (SR_SE < 0._rkind) then      ! runoff slightly below zero but within tolerance
    SR_SE = 0._rkind
   end if

  ! interface surface runoff variables to surfaceFlx output object
   scalarSurfaceRunoff_IE    = SR_IE       ! infiltration excess runoff 
   scalarSurfaceRunoff_SE    = SR_SE       ! saturation excess runoff   
   scalarSurfaceRunoff       = SR_IE+SR_SE ! total surface runoff

  ! check total surface runoff
  ! note: - due to combinations of independent methods, it is possible that runoff exceeds RPM in extreme cases
   if (scalarSurfaceRunoff > scalarRainPlusMelt) then
    ! err=10;
    ! message=trim(message)//&
    ! &"update_gather_runoff_components: sum of infiltration and saturation excess surface runoff components exceeds rain plus melt";
    ! return_flag=.true.; 
    return
   end if

  ! interface surface infiltration variable to surfaceFlx output object
   scalarSurfaceInfiltration = scalarRainPlusMelt - scalarSurfaceRunoff ! surface infiltration  

  ! compute soil control factor
  ! note: infiltration = scalarSoilControl * p
   if (scalarRainPlusMelt > 0._rkind) then
    scalarSoilControl = scalarSurfaceInfiltration/scalarRainPlusMelt
   else
    scalarSoilControl = 0._rkind
   end if

  ! interface infiltration derivatives w.r.t state variables to surfaceFlx output object
  ! note: rain plus melt (RPM) derivatives are zero within soil layers
  if (updateInfil) then
    do iLayer=0,size(dq_dHydStateVec)-1
      dq_dHydStateVec(iLayer) = dq_dHydStateVec_IE(iLayer) + dq_dHydStateVec_SE(iLayer) ! infiltration derivative w.r.t hydrology state variable 
      dq_dNrgStateVec(iLayer) = dq_dNrgStateVec_IE(iLayer) + dq_dNrgStateVec_SE(iLayer) ! infiltration derivative w.r.t energy state variable 
    end do
  else
    do iLayer=0,size(dq_dHydStateVec)-1
      dq_dHydStateVec(iLayer) = realMissing ! not used, so cause problems
      dq_dNrgStateVec(iLayer) = realMissing ! not used, so cause problems
    end do
  end if



 end subroutine update_gather_runoff_components

  attributes(device) subroutine update_surfaceFlx_zero_IE(SR_IE, dq_dHydStateVec_IE, dq_dNrgStateVec_IE)
  real(rkind),intent(inout) :: SR_IE
  real(rkind),intent(inout) :: dq_dHydStateVec_IE(0:)
  real(rkind),intent(inout) :: dq_dNrgStateVec_IE(0:)
  ! **** Update operations for surfaceFlx: zero infiltration excess surface runoff ****
  ! set infiltration excess components
  ! note: it is assumed that rain plus melt does not depend on state variables for infiltration derivatives
  SR_IE              = 0._rkind ! surface runoff
  dq_dHydStateVec_IE = 0._rkind ! surface infiltration derivative w.r.t hydrology state variable
  dq_dNrgStateVec_IE = 0._rkind ! surface infiltration derivative w.r.t energy state variable
 end subroutine update_surfaceFlx_zero_IE 

 attributes(device) subroutine update_surfaceFlx_zero_SE(SR_SE, dq_dHydStateVec_SE, dq_dNrgStateVec_SE)
 real(rkind),intent(inout) :: SR_SE
 real(rkind),intent(inout) :: dq_dHydStateVec_SE(0:)
 real(rkind),intent(inout) :: dq_dNrgStateVec_SE(0:)
  ! **** Update operations for surfaceFlx: zero saturation excess surface runoff ****
  ! set saturation excess components
  ! note: it is assumed that rain plus melt does not depend on state variables for infiltration derivatives
  SR_SE              = 0._rkind ! surface runoff
  dq_dHydStateVec_SE = 0._rkind ! surface infiltration derivative w.r.t hydrology state variable
  dq_dNrgStateVec_SE = 0._rkind ! surface infiltration derivative w.r.t energy state variable
 end subroutine update_surfaceFlx_zero_SE 


 attributes(device) subroutine update_surfaceFlx_FUSE_PRMS(ixRichards, &
 FUSE_Ac_max, FUSE_phi_tens, theta_sat, &
 scalarRainPlusMelt, &
 mLayerVolFracLiq, mLayerTemp, mLayerMatricHead, &
 dTheta_dTk, dTheta_dPsi, &
 nSoil, &
 mLayerDepth, iLayerHeight, &
 dVolFracLiq_dWat, dVolFracLiq_dTk, dq_dHydStateVec_SE, dq_dNrgStateVec_SE, SR_SE, Tcrit)
  !  USE soil_utils_module,only:volFracLiq            ! compute volumetric fraction of liquid water as a function of matric head (-)
  ! USE soil_utils_module,only:hydCond_psi           ! compute hydraulic conductivity as a function of matric head (m s-1)
  ! USE soil_utils_module,only:hydCond_liq           ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
  ! USE soil_utils_module,only:dPsi_dTheta           ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
  ! USE soil_utils_module,only:crit_soilT            ! compute critical temperature below which ice exists
  ! USE soil_utils_module,only:gammp,gammp_complex   ! compute the regularized lower incomplete Gamma function

  ! **** Update operations for surfaceFlx: surface runoff from Clark et al. (2008, doi:10.1029/2007WR006735) -- PRMS ****
  ! note: this parameterization utilizes saturation excess surface runoff only
  use soil_utils_module,only:LogSumExp  ! smooth max/min
  use soil_utils_module,only:SoftArgMax ! smooth arg max/min (for derivatives of LogSumExp)
  ! input: model control
  integer(i4b),intent(in) :: ixRichards         ! index defining the option for Richards' equation (moisture or mixdform)
  ! input: parameters
  real(rkind),intent(in) :: FUSE_Ac_max, FUSE_phi_tens
  real(rkind),intent(in) :: theta_sat         ! soil porosity (-)
  ! input: rain plus melt
  real(rkind),intent(in) :: scalarRainPlusMelt ! rain plus melt  (m s-1)
  ! input: state and diagnostic variables
  real(rkind),intent(in) :: mLayerVolFracLiq(:)  ! volumetric liquid water content in each soil layer (-)
  real(rkind),intent(in) :: mLayerTemp(:)          ! temperature (K)
  real(rkind),intent(in) :: mLayerMatricHead(:)    ! matric head in each soil layer (m)
  ! input: pre-computed derivatives in ...
  real(rkind),intent(in) :: dTheta_dTk(:)          ! ... volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in) :: dTheta_dPsi(:)         ! ... the soil water characteristic w.r.t. psi (m-1)
  ! input: soil layers
  integer(i4b),intent(in) :: nSoil ! number of soil layers
  real(rkind),intent(in) :: mLayerDepth(:)         ! depth of upper-most soil layer (m)
  real(rkind),intent(in) :: iLayerHeight(0:)      ! height at the interface of each layer for soil layers only (m)
  ! inout: local to surfaceFlux
  real(rkind),intent(inout) :: dVolFracLiq_dWat(:)  ! ... vol fraction of liquid w.r.t. water state variable in root layers
  real(rkind),intent(inout) :: dVolFracLiq_dTk(1:nSoil)   ! ... vol fraction of liquid w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dq_dHydStateVec_SE(0:) ! derivative of infiltration w.r.t hydrology state variable (saturation excess component)
  real(rkind),intent(inout) :: dq_dNrgStateVec_SE(0:) ! derivative of infiltration w.r.t energy state variable (saturation excess component)
  real(rkind),intent(inout) :: SR_SE ! saturation excess surface runoff component
  real(rkind),intent(inout)                      :: Tcrit                               ! temperature where all water is unfrozen (K)

  ! local variables
  integer(i4b)                     :: iLayer                              ! index of soil layer
  real(rkind),parameter :: alpha_LSE=1.e3_rkind                ! smoothness parameter for LSE smoother function
  real(rkind)           :: Ac_max                              ! maximum saturated area (-)
  real(rkind)           :: phi_tens                            ! fraction of total storage as tension storage (m)
  real(rkind)           :: Ac                                  ! saturated area (-)
  real(rkind)           :: S1                                  ! total water content in upper soil layer (m)
  real(rkind)           :: S1_max                              ! Maximum storage in the upper layer (m)
  real(rkind)           :: S1_T                                ! tension water content in upper soil layer (m)
  real(rkind)           :: S1_T_max                            ! maximum tension water content in upper soil layer (m)
  real(rkind)           :: dS1_dWat(1:nSoil)   ! derivative of S1 w.r.t. water content
  real(rkind)           :: S1_T_derivatives(1:2)               ! array of derivatives for S1_T
  real(rkind)           :: dS1_T_dS1                           ! derivative of S1_T w.r.t S1
  real(rkind)           :: dS1_T_dWat(1:nSoil) ! derivative of S1_T w.r.t water content
  real(rkind)           :: dAc_dWat(1:nSoil)   ! derivative of Ac w.r.t water content 
  integer(i4b) :: err

  ! validation of parameters
   ! interface input parameters
   Ac_max   = FUSE_Ac_max 
   phi_tens = FUSE_phi_tens
   ! validate input parameters 
   if ((Ac_max<0._rkind).or.(Ac_max>1._rkind)) then
    ! err=10; message=trim(message)//"FUSE PRMS surface runoff error: invalid Ac_max (max saturated area) value"; return_flag=.true.; 
    return
   end if
   if ((phi_tens<0._rkind).or.(phi_tens>1._rkind)) then
    ! err=10; message=trim(message)//"FUSE PRMS surface runoff error: invalid phi_tens (tension storage fraction) value"; return_flag=.true.; 
    return
   end if

  ! compute water content in upper FUSE layer
   S1     = sum( mLayerDepth(:) * mLayerVolFracLiq(:) ) ! total water content in upper FUSE layer (m)
   S1_max = iLayerHeight(nSoil) * theta_sat             ! max water storage for upper FUSE layer (m)

  ! compute tension water content
   S1_T_max = phi_tens * S1_max
   S1_T     = LogSumExp(-alpha_LSE,[S1,S1_T_max],err) ! smooth approximation to S1_T=min(S1,S1_T_max)
   if (err /= 0) then
    ! err=10; message=trim(message)//"FUSE PRMS surface runoff: error in LogSumExp"; return_flag=.true.; 
    return
   end if
   if (S1_T < 0._rkind) then ! check for errors
    ! err=10; message=trim(message)//"FUSE PRMS surface runoff: S1_T is negative (may need to adjust magnitude of alpha_LSE)"
    ! return_flag=.true.; 
    return
   end if

  ! compute saturated area
  Ac = (S1_T/S1_T_max)*Ac_max

  ! compute surface runoff
   SR_SE = Ac * scalarRainPlusMelt ! saturation excess surface runoff component

  ! * compute the derivatives for infiltration *
   if (updateInfil) then

    ! compute derivatives needed for infiltration derivative
    dS1_dWat          = mLayerDepth(:)                       ! derivative of S1 w.r.t. water content
    call SoftArgMax(-alpha_LSE,[S1,S1_T_max], S1_T_derivatives) ! compute vector of derivatives for S1_T
    dS1_T_dS1         = S1_T_derivatives(1)                  ! extract S1_T derivative w.r.t S1
    dS1_T_dWat        = dS1_T_dS1 * dS1_dWat(:)              ! derivative of S1_T w.r.t water content
    dAc_dWat          = (dS1_T_dWat(:)/S1_T_max)*Ac_max      ! derivative of Ac w.r.t water content 

    ! process liquid derivatives
    dVolFracLiq_dWat(:) = 0._rkind ! w.r.t hydrology state variable (depends on form of Richards' equation)
    dVolFracLiq_dTk(:)  = 0._rkind ! w.r.t to energy (temperature) state variable
    select case(ixRichards) ! form of Richards' equation
     case(moisture) ! water content state variable
       dVolFracLiq_dWat(:) = 1._rkind
     case(mixdform) ! pressure head state variable (also take freezing into account)
       do iLayer = 1,nSoil
         Tcrit = crit_soilT_wrapper( mLayerMatricHead(iLayer) )
         if (mLayerTemp(iLayer) < Tcrit) then ! water is frozen in the soil layer
           dVolFracLiq_dWat(iLayer) = 0._rkind
         else                                 ! water is unfrozen -- use water retention curve
           dVolFracLiq_dWat(iLayer) = dTheta_dPsi(iLayer)
         end if
       end do
    !  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
    end select 
    dVolFracLiq_dTk(:) = dTheta_dTk(:) ! already zeroed out if not below critical temperature

    ! * compute the hydrology derivatives (only saturation excess components for FUSE) *
    ! scalarSurfaceInfiltration = scalarRainPlusMelt - scalarRainPlusMelt*Ac
    ! note: rain plus melt derivatives are zero in soil layers
    dq_dHydStateVec_SE(:) = -scalarRainPlusMelt * dAc_dWat(:) * dVolFracLiq_dWat(:) 

    ! * compute the energy derivatives (only saturation excess components for FUSE) *
    ! energy state variable is temperature (transformed outside soilLiqFlx_module if needed)
    dq_dNrgStateVec_SE(:) = -scalarRainPlusMelt * dAc_dWat(:) * dVolFracLiq_dTk(:) 

   end if


 end subroutine update_surfaceFlx_FUSE_PRMS

 attributes(device) subroutine update_surfaceFlx_FUSE_ARNO_VIC(ixRichards,&
 FUSE_b, theta_sat, &
 scalarRainPlusMelt, &
 mLayerTemp, mLayerMatricHead, mLayerVolFracLiq, &
 dTheta_dTk, dTheta_dPsi, &
 nSoil, mLayerDepth, iLayerHeight, &
 dVolFracLiq_dWat, dVolFracLiq_dTk, SR_SE, dq_dHydStateVec_SE, dq_dNrgStateVec_SE, Tcrit)
  ! **** Update operations for surfaceFlx: surface runoff from Clark et al. (2008, doi:10.1029/2007WR006735) -- ARNO/VIC ****
  ! note: this parameterization utilizes saturation excess surface runoff only
  use soil_utils_module,only:LogSumExp  ! smooth max/min
  use soil_utils_module,only:SoftArgMax ! smooth arg max/min (for derivatives of LogSumExp)
  implicit none
  ! input: model control
  integer(i4b),intent(in) :: ixRichards ! index defining the option for Richards' equation (moisture or mixdform)
  ! input: parameters
  real(rkind),intent(in) :: FUSE_b
  real(rkind),intent(in) :: theta_sat ! soil porosity (-)
  ! input: rain plus melt
  real(rkind),intent(in) :: scalarRainPlusMelt ! rain plus melt  (m s-1)
  ! input: state and diagnostic variables
  real(rkind),intent(in) :: mLayerTemp(:) ! temperature (K)
  real(rkind),intent(in) :: mLayerMatricHead(:) ! matric head in each soil layer (m)
  real(rkind),intent(in) :: mLayerVolFracLiq(:) ! volumetric liquid water content in each soil layer (-)
  ! input: pre-computed derivatives in ...
  real(rkind),intent(in) :: dTheta_dTk(:) ! ... volumetric liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in) :: dTheta_dPsi(:) ! ... the soil water characteristic w.r.t. psi (m-1)
  ! input: soil layers
  integer(i4b),intent(in) :: nSoil ! number of soil layers
  real(rkind),intent(in) :: mLayerDepth(:) ! depth of upper-most soil layer (m)
  real(rkind),intent(in) :: iLayerHeight(0:)      ! height at the interface of each layer for soil layers only (m)
  ! inout: local to surfaceFlux
  real(rkind),intent(inout) :: dVolFracLiq_dWat(:)  ! ... vol fraction of liquid w.r.t. water state variable in root layers
  real(rkind),intent(inout) :: dVolFracLiq_dTk(:)   ! ... vol fraction of liquid w.r.t. temperature in root layers
  real(rkind),intent(inout) :: SR_SE ! saturation excess surface runoff component
  real(rkind),intent(inout) :: dq_dHydStateVec_SE(0:) ! derivative of infiltration w.r.t hydrology state variable (saturation excess component)
  real(rkind),intent(inout) :: dq_dNrgStateVec_SE(0:) ! derivative of infiltration w.r.t energy state variable (saturation excess component)
  real(rkind),intent(inout)                      :: Tcrit                               ! temperature where all water is unfrozen (K)


  ! local variables
  integer(i4b)                     :: iLayer                              ! index of soil layer
  logical(lgt),parameter :: smoother = .false.                ! control for optional smoothing in base variable  
  real(rkind) ,parameter :: alpha_LSE=1.e3_rkind              ! smoothness parameter for LSE smoother function
  real(rkind)            :: b                                 ! ARNO/VIC exponent (-) 
  real(rkind)            :: S1                                ! total water content in upper FUSE layer (m)
  real(rkind)            :: dS1_dWat(1:nSoil) ! derivative of S1 w.r.t. water content
  real(rkind)            :: S1_max                            ! Maximum storage in the FUSE layer (m)
  real(rkind)            :: S1_star                           ! total water content in upper FUSE layer computed with a smoothed min (m)
  real(rkind)            :: dS1_star_dS1                      ! derivative in S1_star w.r.t S1
  real(rkind)            :: base                              ! base used in saturated area formula
  real(rkind)            :: dbase_dS1                         ! derivative of base w.r.t S1
  real(rkind)            :: Ac                                ! saturated area (-)
  real(rkind)            :: dAc_dWat(1:nSoil) ! derivative of Ac w.r.t water content 
  real(rkind)            :: S1_star_derivatives(1:2)          ! array of derivatives for S1_star from SoftArgMax function
  real(rkind)            :: roundoff_tolerance                ! tolerance for round-off error
  integer(i4b) :: err

  ! validation of input parameters
  b = FUSE_b ! interface ARNO/VIC exponent
   if ((b < 0.001_rkind).or.(b > 3._rkind)) then
    !err=10; message=trim(message)//"FUSE ARNO/VIC exponent must be between 0.001 and 3"; return_flag=.true.; 
    return
   end if

  ! compute water content in upper FUSE layer
   S1     = sum( mLayerDepth(:) * mLayerVolFracLiq(:) ) ! total water content in upper FUSE layer (m)
   S1_max = iLayerHeight(nSoil) * theta_sat             ! max water storage for upper FUSE layer (m)

  ! compute saturated area
  ! Original FUSE: Ac = 1 - (1-S1/S1_max)**b
  ! Optional: - smoothed to prevent negative bases using a smooth approximation of S1_star = min(S1,S1_max)
  !           - (Smoothed Ac) = 1 - (1-S1_star/S1_max)**b 
   ! compute S1_star (smooth approximation of min(S1,S1_max))
   if (smoother) then ! with smooth approximation of min(S1,S1_max)
    S1_star = LogSumExp(-alpha_LSE,[S1,S1_max],err) ! smooth approximation of min(S1,S1_max) to prevent negative bases
    if (err /= 0) then
    !  err=10; message=trim(message)//"FUSE ARNO/VIC surface runoff: error in LogSumExp"; return_flag=.true.; 
      return
    end if
   else               ! no smoothing
    S1_star = S1
   end if
   if (S1_star < 0._rkind) then ! check for errors
    ! err=10; message=trim(message)//&
    ! &"FUSE ARNO/VIC surface runoff: S1_star is negative (may need to apply smoothing or increase magnitude of alpha_LSE)"
    ! return_flag=.true.; 
    return
   end if

   ! compute base value
   base = 1._rkind - S1_star/S1_max

   ! validate base value and add tolerance for round-off error
   roundoff_tolerance = 1.e2_rkind * epsilon(1._rkind) ! tolerance for round-off error is near machine epsilon 
   if (base < -roundoff_tolerance) then ! if below zero outside of tolerance
    ! err=10; message=trim(message)//"FUSE ARNO/VIC base value is negative"; return_flag=.true.; 
    return
   else if (base < 0._rkind) then       ! if below zero within tolerance
    base = 0._rkind
   end if
 
  ! compute saturated area
  Ac = 1._rkind - base**b  

  ! compute surface runoff
   SR_SE = Ac * scalarRainPlusMelt ! saturation excess surface runoff component

  ! ** compute the derivatives for infiltration **
   if(updateInfil)then
    ! compute derivatives needed for infiltration derivative
    ! Ac   = 1._rkind - base**b 
    dS1_dWat  = mLayerDepth(:)                                 ! derivative of S1 w.r.t. water content
    if (smoother) then ! with smooth approximation of min(S1,S1_max)
     call SoftArgMax(-alpha_LSE,[S1,S1_max],S1_star_derivatives)
     dS1_star_dS1 = S1_star_derivatives(1)                     ! extract S1_star derivative w.r.t S1
    else               ! no smoothing
     dS1_star_dS1 = 1._rkind                                   ! S1_star = S1 if no smoothing
    end if
    dbase_dS1 = -1._rkind/S1_max * dS1_star_dS1                ! derivative of base w.r.t S1
    dAc_dWat  = -b*base**(b-1._rkind)*dbase_dS1*dS1_dWat(:)    ! derivative of Ac w.r.t water content 

    ! process liquid derivatives
    dVolFracLiq_dWat(:) = 0._rkind
    dVolFracLiq_dTk(:)  = 0._rkind
    select case(ixRichards) ! form of Richards' equation
     case(moisture) ! state variable is water content
       dVolFracLiq_dWat(:) = 1._rkind
     case(mixdform) ! state variable is pressure head
       do iLayer = 1,nSoil
         Tcrit = crit_soilT_wrapper( mLayerMatricHead(iLayer) )
         if (mLayerTemp(iLayer) < Tcrit) then ! frozen layer
           dVolFracLiq_dWat(iLayer) = 0._rkind
         else                                 ! unfrozen layer
           dVolFracLiq_dWat(iLayer) = dTheta_dPsi(iLayer)
         end if
       end do
    !  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
    end select 
    dVolFracLiq_dTk(:) = dTheta_dTk(:) ! already zeroed out if not below critical temperature

    ! * compute the hydrology derivatives (only saturation excess components for FUSE) *
    ! scalarSurfaceInfiltration = scalarRainPlusMelt - scalarRainPlusMelt*Ac
    ! note: rain plus melt derivatives are zero in soil layers
    dq_dHydStateVec_SE(:) = -scalarRainPlusMelt * dAc_dWat(:) * dVolFracLiq_dWat(:)

    ! * compute the energy derivatives components (only saturation excess components for FUSE) *
    ! note: energy state variable is temperature (transformed outside soilLiqFlx_module if needed)
    dq_dNrgStateVec_SE(:) = -scalarRainPlusMelt * dAc_dWat(:) * dVolFracLiq_dTk(:)
   end if


 end subroutine update_surfaceFlx_FUSE_ARNO_VIC

 attributes(device) subroutine update_surfaceFlx_FUSE_TOPMODEL(ixRichards, &
    FUSE_lambda, FUSE_chi, FUSE_mu, FUSE_n, theta_sat, &
    scalarRainPlusMelt, &
    mLayerTemp, mLayerMatricHead, mLayerVolFracLiq, &
    dTheta_dTk, dTheta_dPsi, &
    nSoil, mLayerDepth, iLayerHeight, &
    dVolFracLiq_dWat, dVolFracLiq_dTk, SR_SE, &
    dq_dHydStateVec_SE, dq_dNrgStateVec_SE, Tcrit)
      USE soil_utils_module,only:gammp,gammp_complex   ! compute the regularized lower incomplete Gamma function
  ! **** Update operations for surfaceFlx: surface runoff from Clark et al. (2008, doi:10.1029/2007WR006735) -- TOPMODEL ****
  ! note: this parameterization utilizes saturation excess surface runoff only
    implicit none
   ! input: model control
   integer(i4b),intent(in) :: ixRichards ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: parameters
   real(rkind),intent(in) :: FUSE_lambda
   real(rkind),intent(in) :: FUSE_chi
   real(rkind),intent(in) :: FUSE_mu
   real(rkind),intent(in) :: FUSE_n
   real(rkind),intent(in) :: theta_sat          ! soil porosity (-)
   ! input: rain plus melt
   real(rkind),intent(in) :: scalarRainPlusMelt ! rain plus melt  (m s-1)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: mLayerTemp(:)          ! temperature (K)
   real(rkind),intent(in) :: mLayerMatricHead(:)    ! matric head in each soil layer (m)
   real(rkind),intent(in) :: mLayerVolFracLiq(:)    ! volumetric liquid water content in each soil layer (-)
   ! input: pre-computed derivatives in ...
   real(rkind),intent(in) :: dTheta_dTk(:)          ! ... volumetric liquid water content w.r.t. temperature (K-1)
   real(rkind),intent(in) :: dTheta_dPsi(:)         ! ... the soil water characteristic w.r.t. psi (m-1)
   ! input: soil layers
   integer(i4b),intent(in) :: nSoil               ! number of soil layers
   real(rkind),intent(in) :: mLayerDepth(:)         ! depth of upper-most soil layer (m)
   real(rkind),intent(in) :: iLayerHeight(0:) ! height at the interface of each layer for soil layers only (m)
  ! inout: local to surfaceFlux
  real(rkind),intent(inout) :: dVolFracLiq_dWat(:)  ! ... vol fraction of liquid w.r.t. water state variable in root layers
  real(rkind),intent(inout) :: dVolFracLiq_dTk(:)   ! ... vol fraction of liquid w.r.t. temperature in root layers
  real(rkind),intent(inout) :: SR_SE ! saturation excess surface runoff component
  real(rkind),intent(inout) :: dq_dHydStateVec_SE(0:) ! derivative of infiltration w.r.t hydrology state variable (saturation excess component)
  real(rkind),intent(inout) :: dq_dNrgStateVec_SE(0:) ! derivative of infiltration w.r.t energy state variable (saturation excess component)
  real(rkind),intent(inout)                      :: Tcrit                               ! temperature where all water is unfrozen (K)


  ! * local variables *
  integer(i4b) :: iLayer
  ! runoff and infiltration variables
  real(rkind) :: Ac     ! saturated area (-)
  ! FUSE parameters and variables
  real(rkind) :: lambda ! mean
  real(rkind) :: chi    ! scale
  real(rkind) :: mu     ! offset
  real(rkind) :: phi    ! shape (computed from other parameters)
  ! Gamma distribution parameters and variables
  real(rkind) :: alpha  ! shape
  real(rkind) :: theta  ! scale
  real(rkind) :: x_crit ! critical x (random variable) value
  ! topographic index variables
  real(rkind),parameter :: zeta_upper=1.e3_rkind ! upper limit of integral (approaches infinity, but ~1000 provides an accurate result) 
  real(rkind) :: zeta_crit_n ! critical topographic index value (power-transfomred)
  real(rkind) :: zeta_crit   ! critical topographic index value (log space)
  complex(rkind) :: F1,F2    ! temporary storage for regularized lower incomplete gamma function values
  complex(rkind) :: lambda_n ! mean of the power-transformed topographic index
  ! lower FUSE layer variables
  real(rkind) :: S2_max ! max storage in lower FUSE layer (m)
  real(rkind) :: S2     ! total water content in lower FUSE layer (m)
  real(rkind) :: n      ! TOPMODEL exponent exponent (must be sufficiently large to avoid divergence of lambda_n -- n>=3.5 or so)
  ! derivative variables
  real(rkind) :: dS2_dWat(1:nSoil) ! derivative in S2 w.r.t water content 
  real(rkind) :: dAc_dWat(1:nSoil) ! derivative of Ac w.r.t water content 
  real(rkind) :: dzeta_crit_n_dS2                  ! derivative of zeta_crit_n w.r.t S2
  real(rkind) :: dzeta_crit_dzeta_crit_n           ! derivative of zeta_crit w.r.t zeta_crit_n
  real(rkind) :: dx_crit_dzeta_crit                ! derivative of x_crit w.r.t zeta_crit
  real(rkind) :: dx_crit_dS2                       ! derivative of x_crit w.r.t S2
  real(rkind) :: dgammp_dx_crit                    ! derivative of gammp function in Ac w.r.t x_crit
  ! tolerance variables
  real(rkind) :: roundoff_tolerance                ! tolerance for round-off error

  ! interface FUSE input parameters
  lambda = FUSE_lambda
  chi    = FUSE_chi
  mu     = FUSE_mu
  n      = FUSE_n

  ! compute water content in lower FUSE layer
  ! note: the entire soil column is used
   S2     = sum( mLayerDepth(:) * mLayerVolFracLiq(:) ) ! total water content in upper FUSE layer (m)
   S2_max = iLayerHeight(nSoil) * theta_sat             ! max water storage for upper FUSE layer (m)

  ! validation of parameters
   ! validate gamma distribution parameters
   if ((lambda < 5._rkind ).or.(lambda > 10._rkind)) then
    !err=10; message=trim(message)//"FUSE TOPMODEL lambda value must be between 5 and 10"; return_flag=.true.; 
    return
   end if
   if (lambda <= mu) then
    !err=10; message=trim(message)//"FUSE TOPMODEL lambda value must be greater than mu value"; return_flag=.true.; 
    return
   end if
   if ((chi < 2._rkind ).or.(chi > 5._rkind)) then
    !err=10; message=trim(message)//"FUSE TOPMODEL chi value must be between 2 and 5"; return_flag=.true.; 
    return
   end if
   if ((mu < 2.5_rkind ).or.(mu > 3.5_rkind)) then
    !err=10; message=trim(message)//"FUSE TOPMODEL mu value must be between 2.5 and 3.5"; return_flag=.true.; 
    return
   end if

   if ((n < 3.5_rkind).or.(n > 10._rkind)) then ! validate TOPMODEL exponent to avoid divergence of lambda_n
    !err=10; message=trim(message)//"FUSE TOPMODEL exponent must be between 3.5 and 10"; return_flag=.true.; 
    return
   end if
   if (S2 < 0._rkind) then ! check for negative water content values in the lower FUSE layer
    !err=10; message=trim(message)//"negative water content value detected in lower FUSE layer"; return_flag=.true.; 
    return
   end if
   if (S2 > S2_max) then   ! check if water content in lower FUSE layer exceeds the maximum storage
    !err=10; message=trim(message)//"water content in lower FUSE layer exceeds max storage"; return_flag=.true.; 
    return
   end if

  if (S2 > 0._rkind) then ! if some water is stored in lower FUSE layer

   ! set FUSE parameters - input parameters are lambda, chi, and mu
   phi=(lambda-mu)/chi

   ! set Gamma distribution parameters
   alpha=phi
   theta=chi

   ! * compute the mean power-transformed topographic index *
   ! compute regularized lower incomplete Gamma function values
   F1=gammp_complex(alpha,(-(mu*n - mu*theta - (n - theta)*zeta_upper)/n)/theta)
   F2=gammp_complex(alpha,(-(mu*n - mu*theta)/n)/theta)

   ! mean power-transformed topographic index (translated to Fortran from SageMath)
   lambda_n=(cmplx(-mu + zeta_upper,0._rkind,rkind)**alpha*(F1 - 1)*exp(mu/n)*gamma(alpha)/cmplx(-(mu*n - mu*theta - &
           &(n - theta)*zeta_upper)/(n*theta),0._rkind,rkind)**alpha - cmplx(-mu,0._rkind,rkind)**alpha*(F2 - 1)*exp(mu/n)*gamma(alpha)/&
           &cmplx(-(mu*n - mu*theta)/(n*theta),0._rkind,rkind)**alpha)/(cmplx(theta,0._rkind,rkind)**alpha*gamma(alpha))

   ! compute critical zeta value
   ! note: to obtain physical topography values, only the real part of lambda_n is used 
   zeta_crit_n=lambda_n%re*S2_max/S2 ! power-transformed critical topographic index
   if (zeta_crit_n <= 0._rkind) then
    !  err=10; message=trim(message)//"FUSE TOPMODEL zeta_crit_n <= 0"
    !  return_flag=.true.; 
    return
   end if

   zeta_crit=n*log(zeta_crit_n) ! critical topographic index in log space

   ! transform to x random variable and validate result
   x_crit=zeta_crit-mu
   roundoff_tolerance = 1.e2_rkind * epsilon(1._rkind) ! tolerance for round-off error is near machine epsilon 
   if (x_crit < -roundoff_tolerance) then ! less than zero outside tolerance
    !  err=10; message=trim(message)//"FUSE TOPMODEL zeta_crit must be greater or equal to mu -- &
    !                 &try increasing lambda or decreasing mu";
    !  return_flag=.true.; 
     return
   else if (x_crit < 0._rkind) then       ! less than zero but within tolerance
    x_crit = 0._rkind
   end if

   ! compute saturated area
   Ac = 1._rkind-gammp(alpha,x_crit/theta)

  else ! if no water is stored in lower FUSE layer (based on asymptotic behaviour of integral in eq. 9c of Clark et al. (2008))
   Ac = 0._rkind
  end if

  ! compute surface runoff
   SR_SE = Ac * scalarRainPlusMelt ! saturation excess surface runoff component

  ! ** compute the derivatives for infiltration **

   if(updateInfil)then
    ! compute derivatives needed for infiltration derivative
    if (S2 > 0._rkind) then ! for S2 > 0: Ac = 1._rkind-gammp(alpha,x_crit/theta)
     dS2_dWat  = mLayerDepth(:)                       ! derivative of S2 w.r.t. water content      
     dzeta_crit_n_dS2 = -lambda_n%re*S2_max/S2**2_i4b ! derivative of zeta_crit_n=lambda_n%re*S2_max/S2 w.r.t S2     
     dzeta_crit_dzeta_crit_n = ( n*zeta_crit_n**(n-1._rkind) ) / zeta_crit_n**n    ! derivative of zeta_crit=log(zeta_crit_n**n) w.r.t zeta_crit_n
     dx_crit_dzeta_crit = 1._rkind                                                 ! derivative of x_crit=zeta_crit-mu w.r.t zeta_crit
     dx_crit_dS2 = dx_crit_dzeta_crit * dzeta_crit_dzeta_crit_n * dzeta_crit_n_dS2 ! derivative of x_crit w.r.t S2 via chain rule
     dgammp_dx_crit = ( (x_crit/theta)**(alpha-1._rkind) * exp(-x_crit/theta) )/theta/gamma(alpha) ! derivative of gammp function in Ac w.r.t x_crit 
     dAc_dWat = -dgammp_dx_crit * dx_crit_dS2 * dS2_dWat(:) ! derivative of Ac w.r.t water content via chain rule 
    else ! for S2 = 0: Ac = 0
     dAc_dWat(:) = 0._rkind
    end if

    ! process liquid derivatives
    dVolFracLiq_dWat(:) = 0._rkind
    dVolFracLiq_dTk(:)  = 0._rkind
    select case(ixRichards) ! form of Richards' equation
     case(moisture) ! state variable is water content
       dVolFracLiq_dWat(:) = 1._rkind
     case(mixdform) ! state variable is pressure head
       do iLayer = 1,nSoil
         Tcrit = crit_soilT_wrapper( mLayerMatricHead(iLayer) )
         if (mLayerTemp(iLayer) < Tcrit) then ! frozen layer
           dVolFracLiq_dWat(iLayer) = 0._rkind
         else                                 ! unfrozen layer
           dVolFracLiq_dWat(iLayer) = dTheta_dPsi(iLayer)
         end if
       end do
    !  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
    end select 
    dVolFracLiq_dTk(:) = dTheta_dTk(:) ! already zeroed out if not below critical temperature

    ! * compute the hydrology derivatives (only saturation excess components for FUSE) *
    ! scalarSurfaceInfiltration = scalarRainPlusMelt - scalarRainPlusMelt*Ac
    ! note: rain plus melt derivatives are zero in soil layers
    dq_dHydStateVec_SE(:) = -scalarRainPlusMelt * dAc_dWat(:) * dVolFracLiq_dWat(:)

    ! * compute the energy derivatives components (only saturation excess components for FUSE) *
    ! note: energy state variable is temperature (transformed outside soilLiqFlx_module if needed)
    dq_dNrgStateVec_SE(:) = -scalarRainPlusMelt * dAc_dWat(:) * dVolFracLiq_dTk(:)
   end if


 end subroutine update_surfaceFlx_FUSE_TOPMODEL


 attributes(device)subroutine update_surfaceFlx_prescribedHead(ixRichards, &
    scalarMatricHeadLiq, scalarVolFracLiq, &
    mLayerDepth, &
    upperBoundHead, upperBoundTheta, &
    surfaceSatHydCond, dHydCond_dTemp, iceImpedeFac, &
    vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, &
    surfaceHydCond, surfaceDiffuse, &
    scalarSurfaceRunoff, scalarSurfaceRunoff_IE, scalarSurfaceRunoff_SE, scalarSurfaceInfiltration, &
    scalarSoilControl, dq_dHydStateVec, dq_dNrgStateVec, &
    xMaxInfilRate, scalarInfilArea, scalarFrozenArea, &
    cFlux)
    USE soil_utils_module,only:hydCond_psi           ! compute hydraulic conductivity as a function of matric head (m s-1)
    USE soil_utils_module,only:hydCond_liq           ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
    USE soil_utils_module,only:dPsi_dTheta           ! compute derivative of the soil moisture characteristic w.r.t. theta (m)

  ! **** Update operations for surfaceFlx: prescribed pressure head condition ****
 implicit none
   ! input: model control
   integer(i4b),intent(in) :: ixRichards ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: scalarMatricHeadLiq  ! liquid matric head in the upper-most soil layer (m)
   real(rkind),intent(in) :: scalarVolFracLiq     ! volumetric liquid water content in the upper-most soil layer (-)
   ! input: depth of upper-most soil layer (m)
   real(rkind),intent(in) :: mLayerDepth(:)   ! depth of upper-most soil layer (m)
   ! input: diriclet boundary conditions
   real(rkind),intent(in) :: upperBoundHead    ! upper boundary condition for matric head (m)
   real(rkind),intent(in) :: upperBoundTheta   ! upper boundary condition for volumetric liquid water content (-)
   ! input: transmittance
   real(rkind),intent(in) :: surfaceSatHydCond  ! saturated hydraulic conductivity at the surface (m s-1)
   real(rkind),intent(in) :: dHydCond_dTemp     ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   real(rkind),intent(in) :: iceImpedeFac       ! ice impedence factor in the upper-most soil layer (-)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha            ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n                ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m                ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: theta_sat            ! soil porosity (-)
   real(rkind),intent(in) :: theta_res            ! soil residual volumetric water content (-)
   ! input-output: hydraulic conductivity and diffusivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   real(rkind),intent(inout) :: surfaceHydCond  ! hydraulic conductivity (m s-1)
   real(rkind),intent(inout) :: surfaceDiffuse  ! hydraulic diffusivity at the surface (m2 s-1)
   ! output: runoff and infiltration
   real(rkind),intent(inout) :: scalarSurfaceRunoff        ! surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceRunoff_IE     ! infiltration excess surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceRunoff_SE     ! saturation excess surface runoff (m s-1)
   real(rkind),intent(inout) :: scalarSurfaceInfiltration  ! surface infiltration (m s-1)
   ! output: derivatives in surface infiltration w.r.t. ...
   real(rkind),intent(inout) :: scalarSoilControl   ! soil control on infiltration for derivative
   real(rkind),intent(inout) :: dq_dHydStateVec(0:)     ! ... hydrology state in above soil snow or canopy and every soil layer (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateVec(0:)     ! ... energy state in above soil snow or canopy and every soil layer  (m s-1 K-1)
   ! input-output
   real(rkind),intent(inout) :: xMaxInfilRate, scalarInfilArea, scalarFrozenArea
   ! inout - local to surfaceFlux
   real(rkind),intent(inout)                      :: cFlux                               ! capillary flux (m s-1)

   ! surface runoff iz zero for the head condition
   scalarSurfaceRunoff_IE = 0._rkind ! infiltration excess runoff 
   scalarSurfaceRunoff_SE = 0._rkind ! saturation excess runoff 
   scalarSurfaceRunoff    = 0._rkind 

   ! compute transmission and the capillary flux
   select case(ixRichards)  ! select form of Richards' equation
     case(moisture)
       ! compute the hydraulic conductivity and diffusivity at the boundary
       surfaceHydCond = hydCond_liq(upperBoundTheta,surfaceSatHydCond,theta_res,theta_sat,vGn_m) * iceImpedeFac
       surfaceDiffuse = dPsi_dTheta(upperBoundTheta,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * surfaceHydCond
       ! compute the capillary flux
       cflux = -surfaceDiffuse*(scalarVolFracLiq - upperBoundTheta) / (mLayerDepth(1)*0.5_rkind)
     case(mixdform)
       ! compute the hydraulic conductivity and diffusivity at the boundary
       surfaceHydCond = hydCond_psi(upperBoundHead,surfaceSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
       surfaceDiffuse = realMissing
       ! compute the capillary flux
       cflux = -surfaceHydCond*(scalarMatricHeadLiq - upperBoundHead) / (mLayerDepth(1)*0.5_rkind)
    !  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select  ! end select form of Richards' eqn

   ! compute the total flux
   scalarSurfaceInfiltration = cflux + surfaceHydCond
   scalarSoilControl = 0._rkind 

   ! compute the derivatives at the surface, only has a non-zero value for the upper-most soil layer
   dq_dHydStateVec(:) = 0._rkind
   dq_dNrgStateVec(:) = 0._rkind
   if(updateInfil)then
     select case(ixRichards)  ! select form of Richards' equation
       case(moisture); dq_dHydStateVec(1) = -surfaceDiffuse/(mLayerDepth(1)/2._rkind)
       case(mixdform); dq_dHydStateVec(1) = -surfaceHydCond/(mLayerDepth(1)/2._rkind)
      !  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
     end select
     ! note: energy state variable is temperature (transformed outside soilLiqFlx_module if needed)
     dq_dNrgStateVec(1) = -(dHydCond_dTemp/2._rkind)*(scalarMatricHeadLiq - upperBoundHead)/(mLayerDepth(1)*0.5_rkind) + dHydCond_dTemp/2._rkind
   end if

   ! * additional assignment statements for surfaceFlx input-output object based on presribed head values *
   ! the infiltration is always constrained by the prescribed head so the maximum infiltration rate is set to missing
   xMaxInfilRate    = realMissing ! maximum infiltration rate (m s-1)
   ! no soil ice assumed for FUSE PRMS
   scalarInfilArea  = 1._rkind ! fraction of unfrozen area where water can infiltrate (-)
   scalarFrozenArea = 0._rkind      ! fraction of area that is considered impermeable due to soil ice (-)

 end subroutine update_surfaceFlx_prescribedHead

 attributes(device) subroutine update_surfaceFlx_liquidFlux(ixRichards, ixInfRateMax, surfRun_IE, surfRun_SE, nRoots, nSoil, ixIce, &
 scalarRainPlusMelt, &
  mLayerTemp, mLayerMatricHead, mLayerVolFracLiq, mLayerVolFracIce, &
  dTheta_dTk, dTheta_dPsi, mLayerdPsi_dTheta, &
  mLayerDepth, iLayerHeight, &
  surfaceSatHydCond, &
  theta_sat, rootingDepth, zScale_TOPMODEL, wettingFrontSuction, qSurfScale, soilIceScale, soilIceCV, &
  xMaxInfilRate, scalarInfilArea, scalarSoilControl, scalarFrozenArea, &
  surfaceHydCond, surfaceDiffuse, &
  dq_dHydStateVec, dq_dNrgStateVec, &
  Tcrit, fPart1, fPart2, dPart1, dPart2, dfracCap, dfInfRaw, &
 total_soil_depth, rootZoneLiq, rootZoneIce, availCapacity, depthWettingFront, hydCondWettingFront, &
 fracCap, fInfRaw, alpha, xLimg, &
  dVolFracLiq_dWat, dVolFracIce_dWat, dVolFracLiq_dTk, dVolFracIce_dTk, &
  dRootZoneLiq_dWat, dRootZoneIce_dWat, dRootZoneLiq_dTk, dRootZoneIce_dTk, &
   dDepthWettingFront_dWat, dDepthWettingFront_dTk, &
 dxMaxInfilRate_dWat, dxMaxInfilRate_dTk, &
 dInfilArea_dWat, dInfilArea_dTk, dFrozenArea_dWat, dFrozenArea_dTk, &
 dInfilRate_dWat, dInfilRate_dTk, &
SR_IE, SR_SE, dq_dHydStateVec_IE, dq_dHydStateVec_SE, dq_dNrgStateVec_IE, dq_dNrgStateVec_SE)
  ! **** Update operations for surfaceFlx: flux condition ****
    ! input: model control
   integer(i4b),intent(in) :: ixRichards      ! index defining the option for Richards' equation (moisture or mixdform)
        integer(i4b),intent(in) :: ixInfRateMax ! index defining the maximum infiltration rate method (GreenAmpt, topmodel_GA, noInfiltrationExcess)
   integer(i4b),intent(in) :: surfRun_IE ! index defining the infiltration excess surface runoff method
   integer(i4b),intent(in) :: surfRun_SE ! index defining the saturation excess surface runoff method

   integer(i4b),intent(in) :: nRoots          ! number of layers that contain roots
      integer(i4b),intent(in) :: nSoil           ! number of soil layers
   integer(i4b),intent(in) :: ixIce           ! index of lowest ice layer
    ! input: flux at the upper boundary
   real(rkind),intent(in) :: scalarRainPlusMelt ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)

   ! input: state and diagnostic variables
   real(rkind),intent(in) :: mLayerTemp(:)           ! temperature (K)
   real(rkind),intent(in) :: mLayerMatricHead(:)     ! matric head in each soil layer (m)
   real(rkind),intent(in) :: mLayerVolFracLiq(:)     ! volumetric liquid water content in each soil layer (-)
   real(rkind),intent(in) :: mLayerVolFracIce(:)     ! volumetric ice content in each soil layer (-)
   ! input: pre-computed derivatives in ...
   real(rkind),intent(in) :: dTheta_dTk(:)             ! ... volumetric liquid water content w.r.t. temperature (K-1)
   real(rkind),intent(in) :: dTheta_dPsi(:)            ! ... the soil water characteristic w.r.t. psi (m-1)
   real(rkind),intent(in) :: mLayerdPsi_dTheta(:)      ! ... the soil water characteristic w.r.t. theta (m)
   ! input: depth of soil layers (m)
   real(rkind),intent(in) :: mLayerDepth(:)   ! depth of upper-most soil layer (m)
   real(rkind),intent(in) :: iLayerHeight(0:)  ! height at the interface of each layer for soil layers only (m)
     ! input: transmittance
  real(rkind),intent(in) :: surfaceSatHydCond ! saturated hydraulic conductivity at the surface (m s-1)

   ! input: soil parameters
      real(rkind),intent(in) :: theta_sat ! soil porosity (-)
   real(rkind),intent(in) :: rootingDepth         ! rooting depth (m)
     real(rkind),intent(in) :: zScale_TOPMODEL      ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
  real(rkind),intent(in) :: wettingFrontSuction  ! Green-Ampt wetting front suction (m)
     real(rkind),intent(in) :: qSurfScale ! scaling factor in the surface runoff parameterization (-)
   real(rkind),intent(in) :: soilIceScale         ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
   real(rkind),intent(in) :: soilIceCV            ! soil ice CV in Gamma distribution used to define frozen area (-)

  ! input-output: surface runoff and infiltration flux (m s-1)
  real(rkind),intent(inout) :: xMaxInfilRate ! maximum infiltration rate (m s-1)
   real(rkind),intent(inout) :: scalarInfilArea   ! fraction of unfrozen area where water can infiltrate (-)
      real(rkind),intent(inout) :: scalarSoilControl  ! soil control on infiltration for derivative
   real(rkind),intent(inout) :: scalarFrozenArea  ! fraction of area that is considered impermeable due to soil ice (-)
         ! input-output: hydraulic conductivity and diffusivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   real(rkind),intent(inout) :: surfaceHydCond ! hydraulic conductivity (m s-1)
   real(rkind),intent(inout) :: surfaceDiffuse ! hydraulic diffusivity at the surface (m2 s-1)

   ! output
   real(rkind),intent(inout) :: dq_dHydStateVec(0:)
   real(rkind),intent(inout) :: dq_dNrgStateVec(0:)

   ! inout: local to surfaceFlux
     real(rkind),intent(inout)                      :: Tcrit                               ! temperature where all water is unfrozen (K)
      real(rkind),intent(inout)                      :: fPart1,fPart2                       ! different parts of a function
    real(rkind),intent(inout)                      :: dPart1(:)     ! derivatives for different parts of a function
    real(rkind),intent(inout)                      :: dPart2(:)     ! derivatives for different parts of a function
       real(rkind),intent(inout)                      :: dfracCap(:)   ! derivatives for different parts of a function
   real(rkind),intent(inout)                      :: dfInfRaw(:)   ! derivatives for different parts of a function
    real(rkind),intent(inout)                      :: total_soil_depth                    ! total depth of soil (m)
    real(rkind),intent(inout)                      :: rootZoneLiq                         ! depth of liquid water in the root zone (m)
      real(rkind),intent(inout)                      :: rootZoneIce                         ! depth of ice in the root zone (m)

    real(rkind),intent(inout)                      :: availCapacity                       ! available storage capacity in the root zone (m)
    real(rkind),intent(inout)                      :: depthWettingFront                   ! depth to the wetting front (m)
    real(rkind),intent(inout)                      :: hydCondWettingFront                 ! hydraulic conductivity at the wetting front (m s-1)
   real(rkind),intent(inout)                      :: fracCap                             ! fraction of pore space filled with liquid water and ice (-)
   real(rkind),intent(inout)                      :: fInfRaw                             ! infiltrating area before imposing solution constraints (-)
      real(rkind),intent(inout)                      :: alpha                               ! shape parameter in the Gamma distribution
   real(rkind),intent(inout)                      :: xLimg                               ! upper limit of the integral

   real(rkind),intent(inout) :: dVolFracLiq_dWat(:)  ! ... vol fraction of liquid w.r.t. water state variable in root layers
   real(rkind),intent(inout) :: dVolFracIce_dWat(:)  ! ... vol fraction of ice w.r.t. water state variable in root layers
  real(rkind),intent(inout) :: dVolFracLiq_dTk(:)   ! ... vol fraction of liquid w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dVolFracIce_dTk(:)   ! ... vol fraction of ice w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dRootZoneLiq_dWat(:) ! ... vol fraction of scalar root zone liquid w.r.t. water state variable in root layers
  real(rkind),intent(inout) :: dRootZoneLiq_dTk(:)  ! ... vol fraction of scalar root zone liquid w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dRootZoneIce_dTk(:)  ! ... vol fraction of scalar root zone ice w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dRootZoneIce_dWat(:) ! ... vol fraction of scalar root zone ice w.r.t. water state variable in root layers
      real(rkind),intent(inout) :: dDepthWettingFront_dWat(:) ! ... scalar depth of wetting front w.r.t. water state variable in root layers
    real(rkind),intent(inout) :: dDepthWettingFront_dTk(:)  ! ... scalar depth of wetting front w.r.t. temperature in root layers
    real(rkind),intent(inout) :: dxMaxInfilRate_dWat(:) ! ... scalar max infiltration rate w.r.t. water state variable in root layers
    real(rkind),intent(inout) :: dxMaxInfilRate_dTk(:)  ! ... scalar max infiltration rate w.r.t. temperature in root layers
   real(rkind),intent(inout) :: dInfilArea_dWat(:)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
   real(rkind),intent(inout) :: dInfilArea_dTk(:)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  real(rkind),intent(inout) :: dFrozenArea_dWat(:) ! ... scalar frozen area w.r.t. water state variable in canopy or snow and root layers
  real(rkind),intent(inout) :: dFrozenArea_dTk(:)  ! ... scalar frozen area w.r.t. temperature in canopy or snow and root layers
  real(rkind),intent(inout) :: dInfilRate_dWat(:)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
  real(rkind),intent(inout) :: dInfilRate_dTk(:)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
     real(rkind),intent(inout) :: SR_IE, SR_SE
  real(rkind),intent(inout) :: dq_dHydStateVec_IE(0:) ! derivative of infiltration w.r.t hydrology state variable (infiltration excess component)
  real(rkind),intent(inout) :: dq_dHydStateVec_SE(0:) ! derivative of infiltration w.r.t hydrology state variable (saturation excess component)
  real(rkind),intent(inout) :: dq_dNrgStateVec_IE(0:) ! derivative of infiltration w.r.t energy state variable (infiltration excess component)
  real(rkind),intent(inout) :: dq_dNrgStateVec_SE(0:) ! derivative of infiltration w.r.t energy state variable (saturation excess component)


  ! THIS WOULD BE A LOT CLEANER IF IT WAS ALL IN ONE SUBROUTINE JUST LIKE THE OTHERS, FIX
  ! note: the routine may be too long if we combine - this routine is more complicated than all the others
  ! -- main computations
  call update_surfaceFlx_liquidFlux_computation_root_layers(ixRichards, nRoots, &
  mLayerTemp, mLayerMatricHead, mLayerVolFracLiq, mLayerVolFracIce, &
  dTheta_dTk, dTheta_dPsi, mLayerdPsi_dTheta, &
  mLayerDepth, iLayerHeight, &
  rootingDepth, &
  dVolFracLiq_dWat, dVolFracIce_dWat, dVolFracLiq_dTk, dVolFracIce_dTk, &
  dRootZoneLiq_dWat, dRootZoneIce_dWat, dRootZoneLiq_dTk, dRootZoneIce_dTk, &
  Tcrit, rootZoneLiq, rootZoneIce)
  call update_surfaceFlx_liquidFlux_computation_available_capacity(theta_sat, rootingDepth, &
  rootZoneLiq, rootZoneIce, availCapacity); 
  ! if (return_flag) return 
  call update_surfaceFlx_liquidFlux_computation_wetting_front(ixInfRateMax, &
 mLayerDepth, &
 surfaceSatHydCond, &
 zScale_TOPMODEL, rootingDepth, wettingFrontSuction, &
 xMaxInfilRate, &
 fPart1, fPart2, dPart1, dPart2, &
 total_soil_depth, rootZoneLiq, availCapacity, depthWettingFront, hydCondWettingFront, &
 dRootZoneLiq_dWat, dRootZoneIce_dWat, dRootZoneLiq_dTk, dRootZoneIce_dTk, &
 dDepthWettingFront_dWat, dDepthWettingFront_dTk, &
 dxMaxInfilRate_dWat, dxMaxInfilRate_dTk)
  call update_surfaceFlx_liquidFlux_computation_infiltrating_area(nSoil, &
 qSurfScale, &
 scalarInfilArea, &
 dfracCap, dfInfRaw, rootZoneLiq, availCapacity, fracCap, fInfRaw, &
 dRootZoneLiq_dWat, dRootZoneIce_dWat, dRootZoneLiq_dTk, dRootZoneIce_dTk, &
 dInfilArea_dWat, dInfilArea_dTk)
  call update_surfaceFlx_liquidFlux_computation_validate_infiltration(nRoots,ixIce, &
 mLayerVolFracLiq, &
 mLayerDepth, &
 theta_sat, &
 scalarInfilArea)
  call update_surfaceFlx_liquidFlux_computation_impermeable_area(nSoil, &
 scalarRainPlusMelt, &
 soilIceScale, soilIceCV, &
 xMaxInfilRate, scalarFrozenArea, &
 rootZoneIce, alpha, xLimg, &
 dxMaxInfilRate_dWat, dxMaxInfilRate_dTk, dFrozenArea_dWat, dFrozenArea_dTk, dInfilRate_dWat, dInfilRate_dTk)
  if (updateInfil) call update_surfaceFlx_liquidFlux_computation_flux_derivatives(surfRun_IE, surfRun_SE, &
 scalarRainPlusMelt, &
 xMaxInfilRate, scalarInfilArea, scalarFrozenArea, &
 dq_dHydStateVec, dq_dNrgStateVec, &
 dInfilArea_dWat, dInfilArea_dTk, dFrozenArea_dWat, dFrozenArea_dTk, &
 dInfilRate_dWat, dInfilRate_dTk, &
 dq_dHydStateVec_IE, dq_dHydStateVec_SE, dq_dNrgStateVec_IE, dq_dNrgStateVec_SE)
  ! -- put it all together
  call update_surfaceFlx_liquidFlux_infiltration(surfRun_IE, surfRun_SE, &
 scalarRainPlusMelt, &
 xMaxInfilRate, scalarInfilArea, scalarSoilControl, scalarFrozenArea, &
 surfaceHydCond, surfaceDiffuse, &
 SR_IE, SR_SE)

 end subroutine update_surfaceFlx_liquidFlux

 attributes(device) subroutine update_surfaceFlx_liquidFlux_computation_root_layers(ixRichards, nRoots, &
  mLayerTemp, mLayerMatricHead, mLayerVolFracLiq, mLayerVolFracIce, &
  dTheta_dTk, dTheta_dPsi, mLayerdPsi_dTheta, &
  mLayerDepth, iLayerHeight, &
  rootingDepth, &
  dVolFracLiq_dWat, dVolFracIce_dWat, dVolFracLiq_dTk, dVolFracIce_dTk, &
  dRootZoneLiq_dWat, dRootZoneIce_dWat, dRootZoneLiq_dTk, dRootZoneIce_dTk, &
  Tcrit, rootZoneLiq, rootZoneIce)
  ! **** Update operations for surfaceFlx: flux condition -- main computations (root layers) ****
   ! input: model control
   integer(i4b),intent(in) :: ixRichards      ! index defining the option for Richards' equation (moisture or mixdform)
   integer(i4b),intent(in) :: nRoots          ! number of layers that contain roots
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: mLayerTemp(:)           ! temperature (K)
   real(rkind),intent(in) :: mLayerMatricHead(:)     ! matric head in each soil layer (m)
   real(rkind),intent(in) :: mLayerVolFracLiq(:)     ! volumetric liquid water content in each soil layer (-)
   real(rkind),intent(in) :: mLayerVolFracIce(:)     ! volumetric ice content in each soil layer (-)
   ! input: pre-computed derivatives in ...
   real(rkind),intent(in) :: dTheta_dTk(:)             ! ... volumetric liquid water content w.r.t. temperature (K-1)
   real(rkind),intent(in) :: dTheta_dPsi(:)            ! ... the soil water characteristic w.r.t. psi (m-1)
   real(rkind),intent(in) :: mLayerdPsi_dTheta(:)      ! ... the soil water characteristic w.r.t. theta (m)
   ! input: depth of soil layers (m)
   real(rkind),intent(in) :: mLayerDepth(:)   ! depth of upper-most soil layer (m)
   real(rkind),intent(in) :: iLayerHeight(0:)  ! height at the interface of each layer for soil layers only (m)
   ! input: soil parameters
   real(rkind),intent(in) :: rootingDepth         ! rooting depth (m)
   ! inout: local to surfaceFlux
   real(rkind),intent(inout) :: dVolFracLiq_dWat(:)  ! ... vol fraction of liquid w.r.t. water state variable in root layers
   real(rkind),intent(inout) :: dVolFracIce_dWat(:)  ! ... vol fraction of ice w.r.t. water state variable in root layers
  real(rkind),intent(inout) :: dVolFracLiq_dTk(:)   ! ... vol fraction of liquid w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dVolFracIce_dTk(:)   ! ... vol fraction of ice w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dRootZoneLiq_dWat(:) ! ... vol fraction of scalar root zone liquid w.r.t. water state variable in root layers
  real(rkind),intent(inout) :: dRootZoneLiq_dTk(:)  ! ... vol fraction of scalar root zone liquid w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dRootZoneIce_dTk(:)  ! ... vol fraction of scalar root zone ice w.r.t. temperature in root layers
  real(rkind),intent(inout) :: dRootZoneIce_dWat(:) ! ... vol fraction of scalar root zone ice w.r.t. water state variable in root layers
  real(rkind),intent(inout)                      :: Tcrit                               ! temperature where all water is unfrozen (K)
  real(rkind),intent(inout)                      :: rootZoneLiq                         ! depth of liquid water in the root zone (m)
  real(rkind),intent(inout)                      :: rootZoneIce                         ! depth of ice in the root zone (m)

  integer(i4b) :: iLayer

   ! process root layers only liquid and ice derivatives, first initialize
   dVolFracLiq_dWat(:) = 0._rkind
   dVolFracIce_dWat(:) = 0._rkind
   dVolFracLiq_dTk(:)  = 0._rkind
   dVolFracIce_dTk(:)  = 0._rkind
   if(updateInfil)then
     if (nRoots > 0) then
       select case(ixRichards)  ! form of Richards' equation
         case(moisture)
           dVolFracLiq_dWat(:) = 1._rkind
           dVolFracIce_dWat(:) = mLayerdPsi_dTheta(:) - 1._rkind
         case(mixdform)
           do iLayer=1,nRoots
             Tcrit = crit_soilT_wrapper( mLayerMatricHead(iLayer) )
             if (mLayerTemp(iLayer) < Tcrit) then
               dVolFracLiq_dWat(iLayer) = 0._rkind
               dVolFracIce_dWat(iLayer) = dTheta_dPsi(iLayer)
             else
               dVolFracLiq_dWat(iLayer) = dTheta_dPsi(iLayer)
               dVolFracIce_dWat(iLayer) = 0._rkind
             end if
           end do
       end select 
       dVolFracLiq_dTk(:) = dTheta_dTk(:) !already zeroed out if not below critical temperature
       dVolFracIce_dTk(:) = -dVolFracLiq_dTk(:) !often can and will simplify one of these terms out
     end if
   endif
 
   ! define the storage in the root zone (m) and derivatives, first initialize
   rootZoneLiq = 0._rkind
   rootZoneIce = 0._rkind
   dRootZoneLiq_dWat(:) = 0._rkind
   dRootZoneIce_dWat(:) = 0._rkind
   dRootZoneLiq_dTk(:)  = 0._rkind
   dRootZoneIce_dTk(:)  = 0._rkind
 
   ! process layers where the roots extend to the bottom of the layer
   if (nRoots > 1) then
     do iLayer=1,nRoots-1
       rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(iLayer)*mLayerDepth(iLayer)
       rootZoneIce = rootZoneIce + mLayerVolFracIce(iLayer)*mLayerDepth(iLayer)
       if(updateInfil)then
         dRootZoneLiq_dWat(iLayer) = dVolFracLiq_dWat(iLayer)*mLayerDepth(iLayer)
         dRootZoneIce_dWat(iLayer) = dVolFracIce_dWat(iLayer)*mLayerDepth(iLayer)
         dRootZoneLiq_dTk(iLayer)  = dVolFracLiq_dTk(iLayer) *mLayerDepth(iLayer)
         dRootZoneIce_dTk(iLayer)  = dVolFracIce_dTk(iLayer) *mLayerDepth(iLayer)
       end if
     end do
   end if
   ! process layers where the roots end in the current layer
   rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
   rootZoneIce = rootZoneIce + mLayerVolFracIce(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
   if(updateInfil)then
     dRootZoneLiq_dWat(nRoots) = dVolFracLiq_dWat(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
     dRootZoneIce_dWat(nRoots) = dVolFracIce_dWat(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
     dRootZoneLiq_dTk(nRoots)  = dVolFracLiq_dTk(nRoots)* (rootingDepth - iLayerHeight(nRoots-1))
     dRootZoneIce_dTk(nRoots)  = dVolFracIce_dTk(nRoots)* (rootingDepth - iLayerHeight(nRoots-1))
   endif

 end subroutine update_surfaceFlx_liquidFlux_computation_root_layers 


  attributes(device) subroutine update_surfaceFlx_liquidFlux_computation_available_capacity(theta_sat, rootingDepth, &
  rootZoneLiq, rootZoneIce, availCapacity)
  ! **** Update operations for surfaceFlx: flux condition -- main computations (check available capacity) ****
  implicit none
   ! input: soil parameters
   real(rkind),intent(in) :: theta_sat ! soil porosity (-)
   real(rkind),intent(in) :: rootingDepth ! rooting depth (m)
   ! inout: local to surfaceFlux
  real(rkind),intent(inout)                      :: rootZoneLiq                         ! depth of liquid water in the root zone (m)
  real(rkind),intent(inout)                      :: rootZoneIce                         ! depth of ice in the root zone (m)
  real(rkind),intent(inout)                      :: availCapacity                       ! available storage capacity in the root zone (m)

  ! compute and check available capacity to hold water (m)

   availCapacity = theta_sat*rootingDepth - rootZoneIce
   if (rootZoneLiq > availCapacity+verySmaller) then
     !message=trim(message)//'liquid water in the root zone exceeds capacity'
     !err=20; return_flag=.true.; 
    return
   end if

 end subroutine update_surfaceFlx_liquidFlux_computation_available_capacity 

 attributes(device) subroutine update_surfaceFlx_liquidFlux_computation_wetting_front(ixInfRateMax, &
 mLayerDepth, &
 surfaceSatHydCond, &
 zScale_TOPMODEL, rootingDepth, wettingFrontSuction, &
 xMaxInfilRate, &
 fPart1, fPart2, dPart1, dPart2, &
 total_soil_depth, rootZoneLiq, availCapacity, depthWettingFront, hydCondWettingFront, &
 dRootZoneLiq_dWat, dRootZoneIce_dWat, dRootZoneLiq_dTk, dRootZoneIce_dTk, &
 dDepthWettingFront_dWat, dDepthWettingFront_dTk, &
 dxMaxInfilRate_dWat, dxMaxInfilRate_dTk)
  ! **** Update operations for surfaceFlx: flux condition -- main computations (wetting front and derivatives) ****
  ! input: model control
  integer(i4b),intent(in) :: ixInfRateMax ! index defining the maximum infiltration rate method (GreenAmpt, topmodel_GA, noInfiltrationExcess)
  ! input: depth of upper-most soil layer (m)
  real(rkind),intent(in) :: mLayerDepth(:) ! depth of upper-most soil layer (m)
  ! input: transmittance
  real(rkind),intent(in) :: surfaceSatHydCond ! saturated hydraulic conductivity at the surface (m s-1)
  ! input: soil parameters
  real(rkind),intent(in) :: zScale_TOPMODEL      ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
  real(rkind),intent(in) :: rootingDepth         ! rooting depth (m)
  real(rkind),intent(in) :: wettingFrontSuction  ! Green-Ampt wetting front suction (m)
  ! input-output: surface runoff and infiltration flux (m s-1)
  real(rkind),intent(inout) :: xMaxInfilRate ! maximum infiltration rate (m s-1)
  ! inout: local to surfaceFlux
    real(rkind),intent(inout)                      :: fPart1,fPart2                       ! different parts of a function
    real(rkind),intent(inout)                      :: dPart1(:)     ! derivatives for different parts of a function
    real(rkind),intent(inout)                      :: dPart2(:)     ! derivatives for different parts of a function
    real(rkind),intent(inout)                      :: total_soil_depth                    ! total depth of soil (m)
    real(rkind),intent(inout)                      :: rootZoneLiq                         ! depth of liquid water in the root zone (m)
    real(rkind),intent(inout)                      :: availCapacity                       ! available storage capacity in the root zone (m)
    real(rkind),intent(inout)                      :: depthWettingFront                   ! depth to the wetting front (m)
    real(rkind),intent(inout)                      :: hydCondWettingFront                 ! hydraulic conductivity at the wetting front (m s-1)
    real(rkind),intent(inout) :: dRootZoneLiq_dWat(:) ! ... vol fraction of scalar root zone liquid w.r.t. water state variable in root layers
    real(rkind),intent(inout) :: dRootZoneIce_dWat(:) ! ... vol fraction of scalar root zone ice w.r.t. water state variable in root layers
    real(rkind),intent(inout) :: dRootZoneLiq_dTk(:)  ! ... vol fraction of scalar root zone liquid w.r.t. temperature in root layers
    real(rkind),intent(inout) :: dRootZoneIce_dTk(:)  ! ... vol fraction of scalar root zone ice w.r.t. temperature in root layers
    real(rkind),intent(inout) :: dDepthWettingFront_dWat(:) ! ... scalar depth of wetting front w.r.t. water state variable in root layers
    real(rkind),intent(inout) :: dDepthWettingFront_dTk(:)  ! ... scalar depth of wetting front w.r.t. temperature in root layers
    real(rkind),intent(inout) :: dxMaxInfilRate_dWat(:) ! ... scalar max infiltration rate w.r.t. water state variable in root layers
    real(rkind),intent(inout) :: dxMaxInfilRate_dTk(:)  ! ... scalar max infiltration rate w.r.t. temperature in root layers

   ! define the depth to the wetting front (m) and derivatives
   total_soil_depth = sum(mLayerDepth)
   depthWettingFront = (rootZoneLiq/availCapacity)*min(rootingDepth, total_soil_depth)
   if(updateInfil)then
     dDepthWettingFront_dWat(:)=( dRootZoneLiq_dWat(:)*min(rootingDepth, total_soil_depth) + dRootZoneIce_dWat(:)*depthWettingFront )/availCapacity
     dDepthWettingFront_dTk(:) =( dRootZoneLiq_dTk(:) *min(rootingDepth, total_soil_depth) + dRootZoneIce_dTk(:)*depthWettingFront  )/availCapacity
    end if

   ! process hydraulic conductivity-controlled infiltration rate
   select case(ixInfRateMax)  ! maximum infiltration rate parameterization
    case(topmodel_GA)
     ! define the hydraulic conductivity at depth=depthWettingFront (m s-1)
     hydCondWettingFront = surfaceSatHydCond * ( (1._rkind - depthWettingFront/total_soil_depth)**(zScale_TOPMODEL - 1._rkind) )
     ! define the maximum infiltration rate (m s-1)
     xMaxInfilRate = hydCondWettingFront*( (wettingFrontSuction + depthWettingFront)/depthWettingFront )  ! maximum infiltration rate (m s-1)
     ! initialize the derivatives
     dxMaxInfilRate_dWat(:) = 0._rkind
     dxMaxInfilRate_dTk(:)  = 0._rkind
     ! define the derivatives
     if(updateInfil)then
       fPart1    = hydCondWettingFront
       fPart2    = (wettingFrontSuction + depthWettingFront)/depthWettingFront
       dPart1(:) = surfaceSatHydCond*(zScale_TOPMODEL - 1._rkind) * ( (1._rkind - depthWettingFront/total_soil_depth)**(zScale_TOPMODEL - 2._rkind) ) * (-dDepthWettingFront_dWat(:))/total_soil_depth
       dPart2(:) = -dDepthWettingFront_dWat(:)*wettingFrontSuction / (depthWettingFront**2_i4b)
       dxMaxInfilRate_dWat(:) = fPart1*dPart2(:) + fPart2*dPart1(:)
       dPart1(:) = surfaceSatHydCond*(zScale_TOPMODEL - 1._rkind) * ( (1._rkind - depthWettingFront/total_soil_depth)**(zScale_TOPMODEL - 2._rkind) ) * (-dDepthWettingFront_dTk(:))/total_soil_depth
       dPart2(:) = -dDepthWettingFront_dTk(:)*wettingFrontSuction / (depthWettingFront**2_i4b)
       dxMaxInfilRate_dTk(:)  = fPart1*dPart2(:) + fPart2*dPart1(:)
     endif
    case(GreenAmpt)
      ! define the hydraulic conductivity at depth=depthWettingFront (m s-1)
      hydCondWettingFront = surfaceSatHydCond ! Green-Ampt assumes homogeneous soil, therefore the whole soil column has the same hydraulic conductivity
      ! define the maximum infiltration rate (m s-1)
      xMaxInfilRate = hydCondWettingFront * (1._rkind + (1._rkind - depthWettingFront/total_soil_depth) * wettingFrontSuction/depthWettingFront) ! Ks * (1 + (Md) * S/F)
      ! define the derivatives
      if(updateInfil)then
        dxMaxInfilRate_dWat(:) = -hydCondWettingFront*wettingFrontSuction*dDepthWettingFront_dWat(:)/depthWettingFront**2_i4b
        dxMaxInfilRate_dTk(:)  = -hydCondWettingFront*wettingFrontSuction*dDepthWettingFront_dTk(:)/depthWettingFront**2_i4b
      endif
    case(noInfiltrationExcess)
      ! define the hydraulic conductivity at depth=depthWettingFront (m s-1)
      !hydCondWettingFront =  surfaceSatHydCond ! this is not needed for this calculation, but keeping it here in case not setting this will cause unanticipated problems down the line
      ! define the maximum infiltration rate (m s-1), derivatives are zero
      xMaxInfilRate = veryBig ! If maximum infiltration is very big we'll never have a rainfall rate that exceeds it, so no infiltration excess
   end select
 end subroutine update_surfaceFlx_liquidFlux_computation_wetting_front

 attributes(device) subroutine update_surfaceFlx_liquidFlux_computation_infiltrating_area(nSoil, &
 qSurfScale, &
 scalarInfilArea, &
 dfracCap, dfInfRaw, rootZoneLiq, availCapacity, fracCap, fInfRaw, &
 dRootZoneLiq_dWat, dRootZoneIce_dWat, dRootZoneLiq_dTk, dRootZoneIce_dTk, &
 dInfilArea_dWat, dInfilArea_dTk)
  ! **** Update operations for surfaceFlx: flux condition -- main computations (infiltrating area) ****
   ! input: model control
   integer(i4b),intent(in) :: nSoil           ! number of soil layers
   ! input: soil parameters
   real(rkind),intent(in) :: qSurfScale ! scaling factor in the surface runoff parameterization (-)
   ! input-output: surface runoff and infiltration flux (m s-1)
   real(rkind),intent(inout) :: scalarInfilArea   ! fraction of unfrozen area where water can infiltrate (-)
   ! inout: local to surfaceFlux
   real(rkind),intent(inout)                      :: dfracCap(:)   ! derivatives for different parts of a function
   real(rkind),intent(inout)                      :: dfInfRaw(:)   ! derivatives for different parts of a function
   real(rkind),intent(inout)                      :: rootZoneLiq                         ! depth of liquid water in the root zone (m)
   real(rkind),intent(inout)                      :: availCapacity                       ! available storage capacity in the root zone (m)
   real(rkind),intent(inout)                      :: fracCap                             ! fraction of pore space filled with liquid water and ice (-)
   real(rkind),intent(inout)                      :: fInfRaw                             ! infiltrating area before imposing solution constraints (-)
   real(rkind),intent(inout) :: dRootZoneLiq_dWat(:) ! ... vol fraction of scalar root zone liquid w.r.t. water state variable in root layers
   real(rkind),intent(inout) :: dRootZoneIce_dWat(:) ! ... vol fraction of scalar root zone ice w.r.t. water state variable in root layers
   real(rkind),intent(inout) :: dRootZoneLiq_dTk(:)  ! ... vol fraction of scalar root zone liquid w.r.t. temperature in root layers
   real(rkind),intent(inout) :: dRootZoneIce_dTk(:)  ! ... vol fraction of scalar root zone ice w.r.t. temperature in root layers
   real(rkind),intent(inout) :: dInfilArea_dWat(:)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
   real(rkind),intent(inout) :: dInfilArea_dTk(:)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  real(rkind),parameter            :: qSurfScaleMax=1000._rkind           ! maximum surface runoff scaling factor (-)
  real(rkind),parameter            :: maxFracCap=0.995_rkind              ! maximum fraction capacity -- used to avoid numerical problems associated with an enormous derivative
  real(rkind),parameter            :: scaleFactor=0.000001_rkind          ! scale factor for the smoothing function (-)


   ! define the infiltrating area and derivatives for the non-frozen part of the cell/basin, first initialize
   dInfilArea_dWat(:) = 0._rkind
   dInfilArea_dTk(:)  = 0._rkind
   if (qSurfScale < qSurfScaleMax) then
     fracCap         = rootZoneLiq/(maxFracCap*availCapacity)                              ! fraction of available root zone filled with water
     fInfRaw         = 1._rkind - exp(-qSurfScale*(1._rkind - fracCap))                          ! infiltrating area -- allowed to violate solution constraints
     scalarInfilArea = min(0.5_rkind*(fInfRaw + sqrt(fInfRaw**2_i4b + scaleFactor)), 1._rkind)   ! infiltrating area -- constrained
     ! define the derivatives
     if(updateInfil)then
       if (0.5_rkind*(fInfRaw + sqrt(fInfRaw**2_i4b + scaleFactor))< 1._rkind) then
         dfracCap(:) = ( dRootZoneLiq_dWat(:)/maxFracCap + dRootZoneIce_dWat(:)*fracCap )/availCapacity
         dfInfRaw(:) = -qSurfScale*dfracCap(:) * exp(-qSurfScale*(1._rkind - fracCap))
         dInfilArea_dWat(:) = 0.5_rkind*dfInfRaw(:) * (1._rkind + fInfRaw/sqrt(fInfRaw**2_i4b + scaleFactor))
         dfracCap(:) = ( dRootZoneLiq_dTk(:)/maxFracCap + dRootZoneIce_dTk(:)*fracCap )/availCapacity
         dfInfRaw(:) = -qSurfScale*dfracCap(:) * exp(-qSurfScale*(1._rkind - fracCap))
         dInfilArea_dTk(:)  = 0.5_rkind*dfInfRaw(:) * (1._rkind + fInfRaw/sqrt(fInfRaw**2_i4b + scaleFactor))
       endif ! else derivatives are zero
     endif
   else
     scalarInfilArea = 1._rkind ! derivatives are zero
   end if
 end subroutine update_surfaceFlx_liquidFlux_computation_infiltrating_area

 attributes(device) subroutine update_surfaceFlx_liquidFlux_computation_validate_infiltration(nRoots,ixIce, &
 mLayerVolFracLiq, &
 mLayerDepth, &
 theta_sat, &
 scalarInfilArea)
  ! **** Update operations for surfaceFlx: flux condition -- main computations (validate infiltration) ****
   ! input: model control
   integer(i4b),intent(in) :: nRoots          ! number of layers that contain roots
   integer(i4b),intent(in) :: ixIce           ! index of lowest ice layer
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: mLayerVolFracLiq(:) ! volumetric liquid water content in each soil layer (-)
   ! input: depth of upper-most soil layer (m)
   real(rkind),intent(in) :: mLayerDepth(:)   ! depth of upper-most soil layer (m)
   ! input: soil parameters
   real(rkind),intent(in) :: theta_sat            ! soil porosity (-)
   ! input-output: surface runoff and infiltration flux (m s-1)
   real(rkind),intent(inout) :: scalarInfilArea   ! fraction of unfrozen area where water can infiltrate (-)
   ! check to ensure we are not infiltrating into a fully saturated column
   if (ixIce<nRoots) then
     if (sum(mLayerVolFracLiq(ixIce+1:nRoots)*mLayerDepth(ixIce+1:nRoots)) > 0.9999_rkind*theta_sat*sum(mLayerDepth(ixIce+1:nRoots))) scalarInfilArea=0._rkind
   end if
 end subroutine update_surfaceFlx_liquidFlux_computation_validate_infiltration

 attributes(device) subroutine update_surfaceFlx_liquidFlux_computation_impermeable_area(nSoil, &
 scalarRainPlusMelt, &
 soilIceScale, soilIceCV, &
 xMaxInfilRate, scalarFrozenArea, &
 rootZoneIce, alpha, xLimg, &
 dxMaxInfilRate_dWat, dxMaxInfilRate_dTk, dFrozenArea_dWat, dFrozenArea_dTk, dInfilRate_dWat, dInfilRate_dTk)
  ! **** Update operations for surfaceFlx: flux condition -- main computations (impermeable area) ****
   ! input: model control
   integer(i4b),intent(in) :: nSoil           ! number of soil layers
   ! input: flux at the upper boundary
   real(rkind),intent(in) :: scalarRainPlusMelt  ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
   ! input: soil parameters
   real(rkind),intent(in) :: soilIceScale         ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
   real(rkind),intent(in) :: soilIceCV            ! soil ice CV in Gamma distribution used to define frozen area (-)
   ! input-output: surface runoff and infiltration flux (m s-1)
   real(rkind),intent(inout) :: xMaxInfilRate     ! maximum infiltration rate (m s-1)
   real(rkind),intent(inout) :: scalarFrozenArea  ! fraction of area that is considered impermeable due to soil ice (-)
   ! inout: local to surfaceFlux
   real(rkind),intent(inout)                      :: rootZoneIce                         ! depth of ice in the root zone (m)
   real(rkind),intent(inout)                      :: alpha                               ! shape parameter in the Gamma distribution
   real(rkind),intent(inout)                      :: xLimg                               ! upper limit of the integral
   real(rkind),intent(inout) :: dxMaxInfilRate_dWat(:) ! ... scalar max infiltration rate w.r.t. water state variable in root layers
   real(rkind),intent(inout) :: dxMaxInfilRate_dTk(:)  ! ... scalar max infiltration rate w.r.t. temperature in root layers
   real(rkind),intent(inout) :: dFrozenArea_dWat(:) ! ... scalar frozen area w.r.t. water state variable in canopy or snow and root layers
   real(rkind),intent(inout) :: dFrozenArea_dTk(:)  ! ... scalar frozen area w.r.t. temperature in canopy or snow and root layers
   real(rkind),intent(inout) :: dInfilRate_dWat(:)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
   real(rkind),intent(inout) :: dInfilRate_dTk(:)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  
   ! define the impermeable area and derivatives due to frozen ground, first initialize
    dFrozenArea_dWat(:) = 0._rkind
    dFrozenArea_dTk(:)  = 0._rkind
   if (rootZoneIce > tiny(rootZoneIce)) then  ! (avoid divide by zero)
     alpha            = 1._rkind/(soilIceCV**2_i4b)        ! shape parameter in the Gamma distribution
     xLimg            = alpha*soilIceScale/rootZoneIce  ! upper limit of the integral
     !if we use this, we will have a derivative of scalarFrozenArea w.r.t. water and temperature in each layer (through mLayerVolFracIce)
     ! Should fix to deal with frozen area in the root zone, calculations would be expensive
     !scalarFrozenArea = 1._rkind - gammp(alpha,xLimg)      ! fraction of frozen area
     !if(updateInfil)then
     !  dFrozenArea_dWat(:) = -dgammp_dx(alpha,xLimg)*(-alpha*soilIceScale/rootZoneIce**2_i4b)*dRootZoneIce_dWat(:)
     !  dFrozenArea_dTk(:)  = -dgammp_dx(alpha,xLimg)*(-alpha*soilIceScale/rootZoneIce**2_i4b)*dRootZoneIce_dTk(:)
     !end if
     scalarFrozenArea = 0._rkind
   else
     scalarFrozenArea = 0._rkind
   end if
   
   ! infiltration rate derivatives, first initialize
    dInfilRate_dWat(:) = 0._rkind
    dInfilRate_dTk(:)  = 0._rkind
   if(updateInfil)then
     if (xMaxInfilRate < scalarRainPlusMelt) then ! = dxMaxInfilRate_d, dependent on layers not at surface
       dInfilRate_dWat(:) = dxMaxInfilRate_dWat(:)
       dInfilRate_dTk(:)  = dxMaxInfilRate_dTk(:)
     end if
   end if
 end subroutine update_surfaceFlx_liquidFlux_computation_impermeable_area

 attributes(device) subroutine update_surfaceFlx_liquidFlux_computation_flux_derivatives(surfRun_IE, surfRun_SE, &
 scalarRainPlusMelt, &
 xMaxInfilRate, scalarInfilArea, scalarFrozenArea, &
 dq_dHydStateVec, dq_dNrgStateVec, &
 dInfilArea_dWat, dInfilArea_dTk, dFrozenArea_dWat, dFrozenArea_dTk, &
 dInfilRate_dWat, dInfilRate_dTk, &
 dq_dHydStateVec_IE, dq_dHydStateVec_SE, dq_dNrgStateVec_IE, dq_dNrgStateVec_SE)
  ! **** Update operations for surfaceFlx: flux condition -- main computations (flux derivatives, nonzero only if updateInfil) ****
    ! input: model control
   integer(i4b),intent(in) :: surfRun_IE ! index defining the infiltration excess surface runoff method
   integer(i4b),intent(in) :: surfRun_SE ! index defining the saturation excess surface runoff method
    ! input: flux at the upper boundary
   real(rkind),intent(in) :: scalarRainPlusMelt ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
   ! input-output: surface runoff and infiltration flux (m s-1)
   real(rkind),intent(inout) :: xMaxInfilRate     ! maximum infiltration rate (m s-1)
   real(rkind),intent(inout) :: scalarInfilArea   ! fraction of unfrozen area where water can infiltrate (-)
   real(rkind),intent(inout) :: scalarFrozenArea  ! fraction of area that is considered impermeable due to soil ice (-)
   ! output
   real(rkind),intent(inout) :: dq_dHydStateVec(0:)
   real(rkind),intent(inout) :: dq_dNrgStateVec(0:)
   ! inout: local to surfaceFlux
  real(rkind),intent(inout) :: dInfilArea_dWat(:)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
  real(rkind),intent(inout) :: dInfilArea_dTk(:)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  real(rkind),intent(inout) :: dFrozenArea_dWat(:) ! ... scalar frozen area w.r.t. water state variable in canopy or snow and root layers
  real(rkind),intent(inout) :: dFrozenArea_dTk(:)  ! ... scalar frozen area w.r.t. temperature in canopy or snow and root layers
  real(rkind),intent(inout) :: dInfilRate_dWat(:)  ! ... scalar infiltration rate w.r.t. water state variable in canopy or snow and root layers
  real(rkind),intent(inout) :: dInfilRate_dTk(:)   ! ... scalar infiltration rate w.r.t. temperature in canopy or snow and root layers
  real(rkind),intent(inout) :: dq_dHydStateVec_IE(0:) ! derivative of infiltration w.r.t hydrology state variable (infiltration excess component)
  real(rkind),intent(inout) :: dq_dHydStateVec_SE(0:) ! derivative of infiltration w.r.t hydrology state variable (saturation excess component)
  real(rkind),intent(inout) :: dq_dNrgStateVec_IE(0:) ! derivative of infiltration w.r.t energy state variable (infiltration excess component)
  real(rkind),intent(inout) :: dq_dNrgStateVec_SE(0:) ! derivative of infiltration w.r.t energy state variable (saturation excess component)

  ! * local variables *
  ! surface runoff component arrays for infiltration derivatives ...
  real(rkind) :: dq_dHyd(0:size(dq_dHydStateVec)-1)    ! ... w.r.t hydrology state variable
  real(rkind) :: dq_dHyd_IE(0:size(dq_dHydStateVec)-1) ! ... w.r.t hydrology state variable (infiltration excess component)
  real(rkind) :: dq_dHyd_SE(0:size(dq_dHydStateVec)-1) ! ... w.r.t hydrology state variable (saturation excess component)
  real(rkind) :: dq_dNrg(0:size(dq_dHydStateVec)-1)    ! ... w.r.t energy state variable
  real(rkind) :: dq_dNrg_IE(0:size(dq_dHydStateVec)-1) ! ... w.r.t energy state variable (infiltration excess component)
  real(rkind) :: dq_dNrg_SE(0:size(dq_dHydStateVec)-1) ! ... w.r.t energy state variable (saturation excess component)     

  ! allocate and initialize surface runoff component arrays for infiltration derivatives ...
  dq_dHyd    = dq_dHydStateVec ! ... w.r.t hydrology state variable
  dq_dHyd_IE = dq_dHydStateVec ! ... w.r.t hydrology state variable (infiltration excess component)  
  dq_dHyd_SE = dq_dHydStateVec ! ... w.r.t hydrology state variable (saturation excess component)
  dq_dNrg    = dq_dNrgStateVec ! ... w.r.t energy state variable
  dq_dNrg_IE = dq_dNrgStateVec ! ... w.r.t energy state variable (infiltration excess component)
  dq_dNrg_SE = dq_dNrgStateVec ! ... w.r.t energy state variable (saturation excess component) 

   ! dq w.r.t. infiltration only, scalarRainPlusMelt accounted for in computJacob module
   dq_dHyd(:) = (1._rkind - scalarFrozenArea)&
              & * ( dInfilArea_dWat(:)*min(scalarRainPlusMelt,xMaxInfilRate) + scalarInfilArea*dInfilRate_dWat(:) )&
              & + (-dFrozenArea_dWat(:))*scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)
   ! energy state variable is temperature (transformed outside soilLiqFlx_module if needed)
   dq_dNrg(:) = (1._rkind - scalarFrozenArea)&
              & * ( dInfilArea_dTk(:) *min(scalarRainPlusMelt,xMaxInfilRate) + scalarInfilArea*dInfilRate_dTk(:)  )&
              & + (-dFrozenArea_dTk(:)) *scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)
   ! compute infiltration excess (IE) and saturation excess (SE) components
   if (scalarRainPlusMelt.gt.xMaxInfilRate) then ! infiltration excess surface runoff (SR) occurs
    ! * saturation excess surface runoff *
    ! SR_SE    = RPM * (1 - InfilArea_unfrozen) ! (rain plus melt) * (saturated area)
    ! Infil_SE = RPM - SR_SE = RPM * (1 - A_frozen) * InfilArea ! infiltration if SE occurs alone
    ! SE infiltration derivatives (RPM derivative is zero in soil layers): 
    dq_dHyd_SE(:) = (1._rkind - scalarFrozenArea)&
                  & * ( dInfilArea_dWat(:)*scalarRainPlusMelt )&
                  & + (-dFrozenArea_dWat(:))*scalarInfilArea*scalarRainPlusMelt
    dq_dNrg_SE(:) = (1._rkind - scalarFrozenArea)&
                  & * ( dInfilArea_dTk(:) *scalarRainPlusMelt )&
                  & + (-dFrozenArea_dTk(:)) *scalarInfilArea*scalarRainPlusMelt
    ! * infiltration excess surface runoff *
    ! SR_IE = SR - SR_SE ! infiltration excess surface runoff 
    ! Infil_IE = RPM - SR_IE = RPM - SR - SR_SE    ! infiltration if IE occurs alone
    ! IE infiltration derivatives (RPM derivative is zero in soil layers): 
    dq_dHyd_IE = dq_dHyd(:) - dq_dHyd_SE(:)
    dq_dNrg_IE = dq_dNrg(:) - dq_dNrg_SE(:)
   else ! infiltration excess runoff does not occur
    ! SR_SE = SR ! saturation excess surface runoff 
    ! Infil_SE = RPM - SR ! infiltration if SE occurs alone
    ! SE infiltration derivatives (RPM derivative is zero in soil layers): 
    dq_dHyd_SE = dq_dHyd(:)
    dq_dNrg_SE = dq_dNrg(:)
    ! SR_IE = 0._rkind ! infiltration excess surface runoff 
    ! Infil_IE = RPM - SR_IE = RPM ! infiltration if IE occurs alone
    ! IE infiltration derivatives (RPM derivative is zero in soil layers): 
    dq_dHyd_IE(:) = 0._rkind
    dq_dNrg_IE(:) = 0._rkind
   end if

  ! interface derivative arrays with surface runoff component variables from surfaceFlx name space
  ! note: model decisions determine which surface runoff components are used 
    select case(surfRun_IE) ! infiltration excess surface runoff
      case(homegrown_IE) ! homegrown infiltration excess surface runoff
        dq_dHydStateVec_IE = dq_dHyd_IE(:) 
        dq_dNrgStateVec_IE = dq_dNrg_IE(:) 
    end select
    select case(surfRun_SE) ! saturation excess surface runoff
      case(homegrown_SE) ! homegrown saturation excess surface runoff
        dq_dHydStateVec_SE = dq_dHyd_SE(:) 
        dq_dNrgStateVec_SE = dq_dNrg_SE(:) 
    end select
 end subroutine update_surfaceFlx_liquidFlux_computation_flux_derivatives

 attributes(device) subroutine update_surfaceFlx_liquidFlux_infiltration(surfRun_IE, surfRun_SE, &
 scalarRainPlusMelt, &
 xMaxInfilRate, scalarInfilArea, scalarSoilControl, scalarFrozenArea, &
 surfaceHydCond, surfaceDiffuse, &
 SR_IE, SR_SE)
  ! **** Update operations for surfaceFlx: flux condition -- final infiltration and runoff calculations ****
  ! input: model control
  integer(i4b),intent(in) :: surfRun_IE ! index defining the infiltration excess surface runoff method
  integer(i4b),intent(in) :: surfRun_SE ! index defining the saturation excess surface runoff method

  ! input: flux at the upper boundary
  real(rkind),intent(in) :: scalarRainPlusMelt ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
   ! input-output: surface runoff and infiltration flux (m s-1)
   real(rkind),intent(inout) :: xMaxInfilRate      ! maximum infiltration rate (m s-1)
   real(rkind),intent(inout) :: scalarInfilArea    ! fraction of unfrozen area where water can infiltrate (-)
   real(rkind),intent(inout) :: scalarSoilControl  ! soil control on infiltration for derivative
   real(rkind),intent(inout) :: scalarFrozenArea   ! fraction of area that is considered impermeable due to soil ice (-)
      ! input-output: hydraulic conductivity and diffusivity at the surface
   ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
   real(rkind),intent(inout) :: surfaceHydCond ! hydraulic conductivity (m s-1)
   real(rkind),intent(inout) :: surfaceDiffuse ! hydraulic diffusivity at the surface (m2 s-1)
   ! inout: local to surfaceFlux
   real(rkind),intent(inout) :: SR_IE, SR_SE


  ! local variables
  real(rkind) :: surfaceInfiltration      ! surface infiltration
  real(rkind) :: surfaceRunoff            ! surface runoff 
  real(rkind) :: surfaceRunoff_IE         ! infiltration excess component of surface runoff 
  real(rkind) :: surfaceRunoff_SE         ! saturation excess component of surface runoff 
  real(rkind) :: scalarInfilArea_unfrozen ! infiltration area that is not frozen

  ! compute infiltration and runoff
   ! unfrozen infiltration area
   scalarInfilArea_unfrozen=(1._rkind - scalarFrozenArea)*scalarInfilArea
   if (xMaxInfilRate > scalarRainPlusMelt) then
     scalarSoilControl = scalarInfilArea_unfrozen
   else
     scalarSoilControl = 0._rkind
   end if
   if(.not.updateInfil) then
     scalarSoilControl = 0._rkind
   end if

   ! compute infiltration (m s-1)
   surfaceInfiltration = scalarInfilArea_unfrozen*min(scalarRainPlusMelt,xMaxInfilRate)
 
   ! compute surface runoff (m s-1)
   surfaceRunoff = scalarRainPlusMelt - surfaceInfiltration
   if (scalarRainPlusMelt.gt.xMaxInfilRate) then ! infiltration excess surface runoff occurs
    ! saturation excess surface runoff
    surfaceRunoff_SE = scalarRainPlusMelt * (1._rkind - scalarInfilArea_unfrozen) ! (rain plus melt) * (saturated area) 
    ! remaining surface runoff is infiltration excess
    surfaceRunoff_IE = surfaceRunoff - surfaceRunoff_SE ! infiltration excess surface runoff     
   else ! infiltration excess runoff does not occur
    surfaceRunoff_SE = surfaceRunoff ! saturation excess surface runoff 
    surfaceRunoff_IE = 0._rkind      ! infiltration excess surface runoff 
   end if

  ! interface with infiltration excess and saturation excess component variables from surfaceFlx name space
  ! note: model decisions determine which surface runoff components are used 
    select case(surfRun_IE) ! infiltration excess surface runoff
      case(homegrown_IE); SR_IE = surfaceRunoff_IE ! homegrown infiltration excess surface runoff
    end select
    select case(surfRun_SE) ! saturation excess surface runoff
      case(homegrown_SE); SR_SE = surfaceRunoff_SE ! homegrown saturation excess surface runoff
    end select

  ! set surface hydraulic conductivity and diffusivity to missing (not used for flux condition)
   surfaceHydCond = realMissing
   surfaceDiffuse = realMissing

 end subroutine update_surfaceFlx_liquidFlux_infiltration


  attributes(global) subroutine iLayerFlux_kernel(ixRichards, nSnow, nGRU, nSoil, &
  mLayerMatricHeadLiqTrial, mLayerVolFracLiqTrial, &
  mLayerHeight, &
  dPsiLiq_dTemp, dHydCond_dTemp, &
  mLayerHydCond, mLayerDiffuse, &
  dHydCond_dVolLiq, dDiffuse_dVolLiq, dHydCond_dMatric, &
  iLayerHydCond, iLayerDiffuse, &
  iLayerLiqFluxSoil, &
  dq_dHydStateAbove, dq_dHydStateBelow, &
  dq_dNrgStateAbove, dq_dNrgStateBelow)
  integer(i4b) :: ixRichards
  integer(i4b) :: nSnow(:) 
  integer(i4b),value :: nSoil, nGRU
  real(rkind),intent(in) :: mLayerMatricHeadLiqTrial(:,:), mLayerVolFracLiqTrial(:,:)
  real(rkind),intent(in) :: mLayerHeight(0:,:)
  real(rkind),intent(in) :: dPsiLiq_dTemp(:,:), dHydCond_dTemp(:,:)
  real(rkind),intent(in) :: mLayerHydCond(:,:), mLayerDiffuse(:,:)
  real(rkind),intent(in) :: dHydCond_dVolLiq(:,:), dDiffuse_dVolLiq(:,:), dHydCond_dMatric(:,:)
  real(rkind),intent(inout) :: iLayerHydCond(0:,:), iLayerDiffuse(0:,:)
  real(rkind),intent(inout) :: iLayerLiqFluxSoil(0:,:)
  real(rkind),intent(inout) :: dq_dHydStateAbove(0:,:), dq_dHydStateBelow(0:,:)
  real(rkind),intent(inout) :: dq_dNrgStateAbove(0:,:), dq_dNrgStateBelow(0:,:)

  integer(i4b) :: iGRU, iLayer

     iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
     iLayer = (blockidx%y-1) * blockdim%y + threadidx%y

     if (iGRU .gt. nGRU .or. iLayer .ge. nSoil) return
     call iLayerFlux(ixRichards, &
    mLayerMatricHeadLiqTrial(iLayer:iLayer+1,iGRU), mLayerVolFracLiqTrial(nSnow(iGRU)+iLayer:nSnow(iGRU)+iLayer+1,iGRU), &
    mLayerHeight(nSnow(iGRU)+iLayer:nSnow(iGRU)+iLayer+1,iGRU), &
    dPsiLiq_dTemp(iLayer:iLayer+1,iGRU), dHydCond_dTemp(iLayer:iLayer+1,iGRU), &
    mLayerHydCond(iLayer:iLayer+1,iGRU), mLayerDiffuse(iLayer:iLayer+1,iGRU), &
    dHydCond_dVolLiq(iLayer:iLayer+1,iGRU), dDiffuse_dVolLiq(iLayer:iLayer+1,iGRU), dHydCond_dMatric(iLayer:iLayer+1,iGRU), &
    iLayerHydCond(iLayer,iGRU), iLayerDiffuse(iLayer,iGRU), &
    iLayerLiqFluxSoil(iLayer,iGRU), &
    dq_dHydStateAbove(iLayer,iGRU), dq_dHydStateBelow(iLayer,iGRU), &
    dq_dNrgStateAbove(iLayer,iGRU), dq_dNrgStateBelow(iLayer,iGRU), iGRU, iLayer)



end subroutine

  ! ***************************************************************************************************************
  ! private subroutine iLayerFlux: compute the fluxes and derivatives at layer interfaces
  ! ***************************************************************************************************************
  attributes(device) subroutine iLayerFlux(ixRichards, &
  nodeMatricHeadLiqTrial, nodeVolFracLiqTrial, &
  nodeHeight, &
  dPsiLiq_dTemp, dHydCond_dTemp, &
  nodeHydCondTrial, nodeDiffuseTrial, &
  dHydCond_dVolLiq, dDiffuse_dVolLiq, dHydCond_dMatric, &
  iLayerHydCond, iLayerDiffuse, &
  iLayerLiqFluxSoil, &
  dq_dHydStateAbove, dq_dHydStateBelow, &
  dq_dNrgStateAbove, dq_dNrgStateBelow, iGRU, iLayer)
  ! (ixTop,ixBot,nSnow,nGRU,decisions,prog_data,flux_data,deriv_data,&
  !   nodeDiffuseTrial,dHydCond_dTemp,dHydCond_dVolLiq,dDiffuse_dVolLiq,dHydCond_dMatric,&
  !                                     nodeVolFracLiqTrial, nodeMatricHeadLiqTrial, iLayerHydCond, iLayerDiffuse)
    ! ---------------------------------------------------------------------------------------------------------------------------
    ! input: model control, state variables, coordinate variables, temperature derivatives, transmittance variables
                                      ! use device_data_types
                                      ! use data_types
  integer(i4b) :: ixRichards
  real(rkind),intent(in) :: nodeMatricHeadLiqTrial(:), nodeVolFracLiqTrial(:)
  real(rkind),intent(in) :: nodeHeight(:)
  real(rkind),intent(in) :: dPsiLiq_dTemp(:), dHydCond_dTemp(:)
  real(rkind),intent(in) :: nodeHydCondTrial(:), nodeDiffuseTrial(:)
  real(rkind),intent(in) :: dHydCond_dVolLiq(:), dDiffuse_dVolLiq(:), dHydCond_dMatric(:)
  real(rkind),intent(inout) :: iLayerHydCond, iLayerDiffuse
  real(rkind),intent(inout) :: iLayerLiqFluxSoil
  real(rkind),intent(inout) :: dq_dHydStateAbove, dq_dHydStateBelow
  real(rkind),intent(inout) :: dq_dNrgStateAbove, dq_dNrgStateBelow
  integer(i4b) :: iGRU, iLayer
    ! type(in_type_iLayerFlux),intent(in)   :: in_iLayerFlux   ! class object for input data
    ! integer(i4b),intent(in) :: ixTop,ixBot,nGRU
    ! integer(i4b),intent(in),device :: nSnow(:)
    ! ! output: transmittance variables and vertical flux at layer interface, derivatives, and error control
    ! type(decisions_device),intent(in) :: decisions
    ! type(prog_data_device),intent(in) :: prog_data
    ! type(flux_data_device),intent(inout) :: flux_data
    ! type(deriv_data_device),intent(inout) :: deriv_data
    ! ! type(out_type_iLayerFlux),intent(out) :: out_iLayerFlux  ! class object for output data
    ! real(rkind),intent(in),device                :: nodeDiffuseTrial(:,:)  ! diffusivity at layer mid-point (m2 s-1)
    ! real(rkind),intent(in),device :: dHydCond_dTemp(:,:)   ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
    ! real(rkind),intent(in),device :: dHydCond_dVolLiq(:,:) ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
    ! real(rkind),intent(in),device :: dDiffuse_dVolLiq(:,:) ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
    ! real(rkind),intent(in),device :: dHydCond_dMatric(:,:)
    ! real(rkind),intent(in),device :: nodeVolFracLiqTrial(:,:), nodeMatricHeadLiqTrial(:,:)
    ! real(rkind),intent(inout),device :: iLayerHydCond(0:,:), iLayerDiffuse(0:,:)
  
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
  call update_iLayerFlux(ixRichards, &
    nodeMatricHeadLiqTrial, nodeVolFracLiqTrial, &
    nodeHeight, &
    dPsiLiq_dTemp, dHydCond_dTemp, &
    nodeHydCondTrial, nodeDiffuseTrial, &
    dHydCond_dVolLiq, dDiffuse_dVolLiq, dHydCond_dMatric, &
    iLayerHydCond, iLayerDiffuse, &
    iLayerLiqFluxSoil, &
    dq_dHydStateAbove, dq_dHydStateBelow, &
    dq_dNrgStateAbove, dq_dNrgStateBelow, &
    useGeometric, &
    ixLower, ixUpper, &
    dPsi, dLiq, dz, cflux, iGRU, iLayer)
   
  end subroutine iLayerFlux

  attributes(device)  subroutine update_iLayerFlux(ixRichards, &
  nodeMatricHeadLiqTrial, nodeVolFracLiqTrial, &
  nodeHeight, &
  dPsiLiq_dTemp, dHydCond_dTemp, &
  nodeHydCondTrial, nodeDiffuseTrial, &
  dHydCond_dVolLiq, dDiffuse_dVolLiq, dHydCond_dMatric, &
  iLayerHydCond, iLayerDiffuse, &
  iLayerLiqFluxSoil, &
  dq_dHydStateAbove, dq_dHydStateBelow, &
  dq_dNrgStateAbove, dq_dNrgStateBelow, &
  useGeometric, &
  ixLower, ixUpper, &
  dPsi, dLiq, dz, cflux, iGRU, iLayer)
  ! **** Update operations for iLayerFlux ****
   ! input: model control
   integer(i4b),intent(in) :: ixRichards ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state variables
   real(rkind),intent(in) :: nodeMatricHeadLiqTrial(:)  ! liquid matric head at the soil nodes (m)
   real(rkind),intent(in) :: nodeVolFracLiqTrial(:)     ! volumetric fraction of liquid water at the soil nodes (-)
   ! input: model coordinate variables
   real(rkind),intent(in) :: nodeHeight(:)  ! height at the mid-point of the lower layer (m)
      ! input: temperature derivatives
   real(rkind),intent(in) :: dPsiLiq_dTemp(:)     ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   real(rkind),intent(in) :: dHydCond_dTemp(:)    ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: transmittance
   real(rkind),intent(in) :: nodeHydCondTrial(:)  ! hydraulic conductivity at layer mid-points (m s-1)
   real(rkind),intent(in) :: nodeDiffuseTrial(:)  ! diffusivity at layer mid-points (m2 s-1)
   ! input: transmittance derivatives
   real(rkind),intent(in) :: dHydCond_dVolLiq(:) ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
   real(rkind),intent(in) :: dDiffuse_dVolLiq(:) ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
   real(rkind),intent(in) :: dHydCond_dMatric(:) ! derivative in hydraulic conductivity w.r.t matric head (m s-1)
   ! output: tranmsmittance at the layer interface (scalars)
   real(rkind),intent(inout) :: iLayerHydCond  ! hydraulic conductivity at the interface between layers (m s-1)
   real(rkind),intent(inout) :: iLayerDiffuse  ! hydraulic diffusivity at the interface between layers (m2 s-1)
   ! output: vertical flux at the layer interface (scalars)
   real(rkind),intent(inout) :: iLayerLiqFluxSoil  ! vertical flux of liquid water at the layer interface (m s-1)
   ! output: derivatives in fluxes w.r.t. ...  
   real(rkind),intent(inout) :: dq_dHydStateAbove ! ... matric head or volumetric liquid water in the layer above (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dHydStateBelow ! ... matric head or volumetric liquid water in the layer below (m s-1 or s-1)
   ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
   real(rkind),intent(inout) :: dq_dNrgStateAbove ! derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
   real(rkind),intent(inout) :: dq_dNrgStateBelow ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
   ! inout: local to iLayerFlux
   logical(lgt),intent(in) :: useGeometric
   integer(i4b),intent(in) :: ixLower, ixUpper
   real(rkind),intent(inout)                      :: dPsi                 ! spatial difference in matric head (m)
   real(rkind),intent(inout)                      :: dLiq                 ! spatial difference in volumetric liquid water (-)
   real(rkind),intent(inout)                      :: dz                   ! spatial difference in layer mid-points (m)
   real(rkind),intent(inout)                      :: cflux                ! capillary flux (m s-1)
   integer(i4b) :: iGRU, iLayer
 
  ! ** compute the fluxes
  call update_iLayerFlux_fluxes(ixRichards, &
  nodeMatricHeadLiqTrial, nodeVolFracLiqTrial, &
  nodeHeight, &
  nodeHydCondTrial, nodeDiffuseTrial, &
  iLayerHydCond, iLayerDiffuse, &
  iLayerLiqFluxSoil, &
  useGeometric, &
  ixLower, ixUpper, &
  dPsi, dLiq, dz, cflux, iGRU, iLayer)

  ! ** compute the derivatives
  call update_iLayerFlux_derivatives(ixRichards, &
  dPsiLiq_dTemp, dHydCond_dTemp, &
  nodeHydCondTrial, nodeDiffuseTrial, &
  dHydCond_dVolLiq, dDiffuse_dVolLiq, dHydCond_dMatric, &
  iLayerHydCond, iLayerDiffuse, &
  dq_dHydStateAbove, dq_dHydStateBelow, &
  dq_dNrgStateAbove, dq_dNrgStateBelow, &
  useGeometric, &
  ixLower, ixUpper, &
  dPsi, dLiq, dz)

 end subroutine update_iLayerFlux


  attributes(device) subroutine update_iLayerFlux_fluxes(ixRichards, &
  nodeMatricHeadLiqTrial, nodeVolFracLiqTrial, &
  nodeHeight, &
  nodeHydCondTrial, nodeDiffuseTrial, &
  iLayerHydCond, iLayerDiffuse, &
  iLayerLiqFluxSoil, &
  useGeometric, &
  ixLower, ixUpper, &
  dPsi, dLiq, dz, cflux, iGRU, iLayer)
  ! **** Update operations for iLayerFlux: compute fluxes ****
   ! input: model control
   integer(i4b),intent(in) :: ixRichards ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state variables
   real(rkind),intent(in) :: nodeMatricHeadLiqTrial(:)  ! liquid matric head at the soil nodes (m)
   real(rkind),intent(in) :: nodeVolFracLiqTrial(:)     ! volumetric fraction of liquid water at the soil nodes (-)
   ! input: model coordinate variables
   real(rkind),intent(in) :: nodeHeight(:)  ! height at the mid-point of the lower layer (m)
   ! input: transmittance
   real(rkind),intent(in) :: nodeHydCondTrial(:)  ! hydraulic conductivity at layer mid-points (m s-1)
   real(rkind),intent(in) :: nodeDiffuseTrial(:)  ! diffusivity at layer mid-points (m2 s-1)
   ! output: tranmsmittance at the layer interface (scalars)
   real(rkind),intent(inout) :: iLayerHydCond  ! hydraulic conductivity at the interface between layers (m s-1)
   real(rkind),intent(inout) :: iLayerDiffuse  ! hydraulic diffusivity at the interface between layers (m2 s-1)
   ! output: vertical flux at the layer interface (scalars)
   real(rkind),intent(inout) :: iLayerLiqFluxSoil  ! vertical flux of liquid water at the layer interface (m s-1)
   ! inout: local to iLayerFlux
   logical(lgt),intent(in) :: useGeometric
   integer(i4b),intent(in) :: ixLower, ixUpper
   real(rkind),intent(inout)                      :: dPsi                 ! spatial difference in matric head (m)
   real(rkind),intent(inout)                      :: dLiq                 ! spatial difference in volumetric liquid water (-)
   real(rkind),intent(inout)                      :: dz                   ! spatial difference in layer mid-points (m)
   real(rkind),intent(inout)                      :: cflux                ! capillary flux (m s-1)
   integer(i4b) :: iGRU, iLayer

   ! compute the vertical flux of liquid water
   ! compute the hydraulic conductivity at the interface
   if (useGeometric) then
     iLayerHydCond   = sqrt(nodeHydCondTrial(ixLower)   * nodeHydCondTrial(ixUpper))
   else
     iLayerHydCond   = (nodeHydCondTrial(ixLower)   + nodeHydCondTrial(ixUpper))*0.5_rkind
   end if
   
   dz = nodeHeight(ixLower) - nodeHeight(ixUpper)
   ! compute the capillary flux
   select case(ixRichards)  ! select form of Richards' equation
     case(moisture)
      iLayerDiffuse = sqrt(nodeDiffuseTrial(ixLower) * nodeDiffuseTrial(ixUpper))
      dLiq          = nodeVolFracLiqTrial(ixLower) - nodeVolFracLiqTrial(ixUpper)
      cflux         = -iLayerDiffuse * dLiq/dz
     case(mixdform)
      iLayerDiffuse = realMissing
      dPsi          = nodeMatricHeadLiqTrial(ixLower) - nodeMatricHeadLiqTrial(ixUpper)
      cflux         = -iLayerHydCond * dPsi/dz
    !  case default; err=10; message=trim(message)//"unable to identify option for Richards' equation"; return_flag=.true.; return
   end select
   ! compute the total flux (add gravity flux, positive downwards)
   iLayerLiqFluxSoil = cflux + iLayerHydCond

 end subroutine update_iLayerFlux_fluxes

  attributes(device) subroutine update_iLayerFlux_derivatives(ixRichards, &
  dPsiLiq_dTemp, dHydCond_dTemp, &
  nodeHydCondTrial, nodeDiffuseTrial, &
  dHydCond_dVolLiq, dDiffuse_dVolLiq, dHydCond_dMatric, &
  iLayerHydCond, iLayerDiffuse, &
  dq_dHydStateAbove, dq_dHydStateBelow, &
  dq_dNrgStateAbove, dq_dNrgStateBelow, &
  useGeometric, &
  ixLower, ixUpper, &
  dPsi, dLiq, dz)
  ! **** Update operations for iLayerFlux: compute derivatives ****
   ! input: model control
   integer(i4b),intent(in) :: ixRichards    ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: temperature derivatives
   real(rkind),intent(in) :: dPsiLiq_dTemp(:)     ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   real(rkind),intent(in) :: dHydCond_dTemp(:)    ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: transmittance
   real(rkind),intent(in) :: nodeHydCondTrial(:)  ! hydraulic conductivity at layer mid-points (m s-1)
   real(rkind),intent(in) :: nodeDiffuseTrial(:)  ! diffusivity at layer mid-points (m2 s-1)
   ! input: transmittance derivatives
   real(rkind),intent(in) :: dHydCond_dVolLiq(:) ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
   real(rkind),intent(in) :: dDiffuse_dVolLiq(:) ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
   real(rkind),intent(in) :: dHydCond_dMatric(:) ! derivative in hydraulic conductivity w.r.t matric head (m s-1)
   ! output: tranmsmittance at the layer interface (scalars)
   real(rkind),intent(inout) :: iLayerHydCond ! hydraulic conductivity at the interface between layers (m s-1)
   real(rkind),intent(inout) :: iLayerDiffuse ! hydraulic diffusivity at the interface between layers (m2 s-1)
   ! output: derivatives in fluxes w.r.t. ...  
   real(rkind),intent(inout) :: dq_dHydStateAbove ! ... matric head or volumetric liquid water in the layer above (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dHydStateBelow ! ... matric head or volumetric liquid water in the layer below (m s-1 or s-1)
   ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
   real(rkind),intent(inout) :: dq_dNrgStateAbove ! derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
   real(rkind),intent(inout) :: dq_dNrgStateBelow ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
   ! inout: local to iLayerFlux
   logical(lgt),intent(in) :: useGeometric
   integer(i4b),intent(in) :: ixLower, ixUpper
   real(rkind),intent(inout)                      :: dPsi                 ! spatial difference in matric head (m)
   real(rkind),intent(inout)                      :: dLiq                 ! spatial difference in volumetric liquid water (-)
   real(rkind),intent(inout)                      :: dz                   ! spatial difference in layer mid-points (m)

  ! * local variables (derivative in Darcy's flux) *
  ! deriviatives at the layer interface
  real(rkind) :: dHydCondIface_dVolLiqAbove  ! hydraulic conductivity w.r.t. volumetric liquid water content in layer above
  real(rkind) :: dHydCondIface_dVolLiqBelow  ! hydraulic conductivity w.r.t. volumetric liquid water content in layer below
  real(rkind) :: dDiffuseIface_dVolLiqAbove  ! hydraulic diffusivity  w.r.t. volumetric liquid water content in layer above
  real(rkind) :: dDiffuseIface_dVolLiqBelow  ! hydraulic diffusivity  w.r.t. volumetric liquid water content in layer below
  real(rkind) :: dHydCondIface_dMatricAbove  ! hydraulic conductivity w.r.t. matric head in layer above
  real(rkind) :: dHydCondIface_dMatricBelow  ! hydraulic conductivity w.r.t. matric head in layer below

   select case(ixRichards)  ! select form of Richards' equation
     case(moisture)
       ! still need to implement arithmetric mean for the moisture-based form
       if (.not.useGeometric) then
        !  message=trim(message)//'only currently implemented for geometric mean -- change local flag'
        !  err=20; return_flag=.true.; return
       end if
       ! derivatives in hydraulic conductivity at the layer interface (m s-1)
       dHydCondIface_dVolLiqAbove = dHydCond_dVolLiq(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_rkind/max(iLayerHydCond,verySmaller)
       dHydCondIface_dVolLiqBelow = dHydCond_dVolLiq(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_rkind/max(iLayerHydCond,verySmaller)
       ! derivatives in hydraulic diffusivity at the layer interface (m2 s-1)
       dDiffuseIface_dVolLiqAbove = dDiffuse_dVolLiq(ixUpper)*nodeDiffuseTrial(ixLower) * 0.5_rkind/max(iLayerDiffuse,verySmaller)
       dDiffuseIface_dVolLiqBelow = dDiffuse_dVolLiq(ixLower)*nodeDiffuseTrial(ixUpper) * 0.5_rkind/max(iLayerDiffuse,verySmaller)
       ! derivatives in the flux w.r.t. volumetric liquid water content
       dq_dHydStateAbove = -dDiffuseIface_dVolLiqAbove*dLiq/dz + iLayerDiffuse/dz + dHydCondIface_dVolLiqAbove
       dq_dHydStateBelow = -dDiffuseIface_dVolLiqBelow*dLiq/dz - iLayerDiffuse/dz + dHydCondIface_dVolLiqBelow
     case(mixdform)
       ! derivatives in hydraulic conductivity
       if (useGeometric) then
         dHydCondIface_dMatricAbove = dHydCond_dMatric(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_rkind/max(iLayerHydCond,verySmaller)
         dHydCondIface_dMatricBelow = dHydCond_dMatric(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_rkind/max(iLayerHydCond,verySmaller)
       else
         dHydCondIface_dMatricAbove = dHydCond_dMatric(ixUpper)/2._rkind
         dHydCondIface_dMatricBelow = dHydCond_dMatric(ixLower)/2._rkind
       end if
       ! derivatives in the flux w.r.t. matric head
       dq_dHydStateAbove = -dHydCondIface_dMatricAbove*dPsi/dz + iLayerHydCond/dz + dHydCondIface_dMatricAbove
       dq_dHydStateBelow = -dHydCondIface_dMatricBelow*dPsi/dz - iLayerHydCond/dz + dHydCondIface_dMatricBelow
       ! derivative in the flux w.r.t. temperature
       dq_dNrgStateAbove = -(dHydCond_dTemp(ixUpper)/2._rkind)*dPsi/dz + iLayerHydCond*dPsiLiq_dTemp(ixUpper)/dz + dHydCond_dTemp(ixUpper)/2._rkind
       dq_dNrgStateBelow = -(dHydCond_dTemp(ixLower)/2._rkind)*dPsi/dz - iLayerHydCond*dPsiLiq_dTemp(ixLower)/dz + dHydCond_dTemp(ixLower)/2._rkind
    !  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select


 end subroutine update_iLayerFlux_derivatives
  
attributes(device) subroutine qDrainFlux_kernel(nGRU, nSoil, nSnow,bc_lower,ixRichards, &
 mLayerMatricHeadLiqTrial, mLayerVolFracLiqTrial, &
 nodeDepth, nodeHeight, &
 lowerBoundHead, lowerBoundTheta, &
 node_dPsi_dTheta, dPsiLiq_dTemp, &
 bottomSatHydCond, iceImpedeFac, surfaceSatHydCond, &
 dHydCond_dVolLiq, dHydCond_dMatric, dHydCond_dTemp, &
 vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, kAnisotropic, zScale_TOPMODEL, &
 mLayerHydCond, mLayerDiffuse, &
 scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat)
    implicit none
   ! input: model control
    integer(i4b),intent(in),value :: nGRU, nSoil
    integer(i4b),intent(in) :: nSnow(:)
   integer(i4b),intent(in) :: bc_lower ! index defining the type of boundary conditions
   integer(i4b),intent(in) :: ixRichards          ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: mLayerMatricHeadLiqTrial(:,:)   ! liquid matric head in the lowest unsaturated node (m)
   real(rkind),intent(in) :: mLayerVolFracLiqTrial(:,:)      ! volumetric liquid water content in the lowest unsaturated node (-)
   ! input: model coordinate variables
   real(rkind),intent(in) :: nodeDepth(:,:)           ! depth of the lowest unsaturated soil layer (m)
   real(rkind),intent(in) :: nodeHeight(0:,:)                ! height of the lowest unsaturated soil node (m)
   ! input: diriclet boundary conditions
   real(rkind),intent(in) :: lowerBoundHead      ! lower boundary condition for matric head (m)
   real(rkind),intent(in) :: lowerBoundTheta     ! lower boundary condition for volumetric liquid water content (-)
   ! input: derivative in soil water characteristic
   real(rkind),intent(in) :: node_dPsi_dTheta(:,:)      ! derivative of the soil moisture characteristic w.r.t. theta (m)
   real(rkind),intent(in) :: dPsiLiq_dTemp(:,:)    ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: transmittance
   real(rkind),intent(in) :: bottomSatHydCond(0:,:)    ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   real(rkind),intent(in) :: iceImpedeFac(:,:)        ! ice impedence factor in the upper-most soil layer (-)
   real(rkind),intent(in) :: surfaceSatHydCond(0:,:)   ! saturated hydraulic conductivity at the surface (m s-1)
   ! input: transmittance derivatives
   real(rkind),intent(in) :: dHydCond_dVolLiq(:,:)     ! derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
   real(rkind),intent(in) :: dHydCond_dMatric(:,:)     ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
   real(rkind),intent(in) :: dHydCond_dTemp(:,:)      ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha(:)           ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n(:)               ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m(:,:)               ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: theta_sat(:)           ! soil porosity (-)
   real(rkind),intent(in) :: theta_res(:)           ! soil residual volumetric water content (-)
   real(rkind),intent(in) :: kAnisotropic           ! anisotropy factor for lateral hydraulic conductivity (-)
   real(rkind),intent(in) :: zScale_TOPMODEL       ! scale factor for TOPMODEL-ish baseflow parameterization (m)
   ! output: hydraulic conductivity at the bottom of the unsaturated zone
   real(rkind),intent(inout) :: mLayerHydCond(:,:)       ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   real(rkind),intent(inout) :: mLayerDiffuse(:,:)       ! hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
   ! output: drainage flux from the bottom of the soil profile
   real(rkind),intent(inout) :: scalarDrainage(0:,:)      ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   real(rkind),intent(inout) :: dq_dHydStateUnsat(0:,:)   ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateUnsat(0:,:)   ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
   integer(i4b) :: iGRU

   iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
   if (iGRU .gt. nGRU) return
       call qDrainFlux(bc_lower,ixRichards, &
 mLayerMatricHeadLiqTrial(nSoil,iGRU), mLayerVolFracLiqTrial(nSoil,iGRU), &
 nodeDepth(nSnow(iGRU)+nSoil,iGRU), nodeHeight(nSnow(iGRU)+nSoil,iGRU), &
 lowerBoundHead, lowerBoundTheta, &
 node_dPsi_dTheta(nSoil,iGRU), dPsiLiq_dTemp(nSoil,iGRU), &
 bottomSatHydCond(nSoil,iGRU), iceImpedeFac(nSoil,iGRU), mLayerHydCond(nSoil,iGRU), surfaceSatHydCond(0,iGRU), &
 dHydCond_dVolLiq(nSoil,iGRU), dHydCond_dMatric(nSoil,iGRU), dHydCond_dTemp(nSoil,iGRU), &
 vGn_alpha(nSoil), vGn_n(nSoil), vGn_m(nSoil,iGRU), theta_sat(nSOil), theta_res(nSoil), kAnisotropic, zScale_TOPMODEL, &
 mLayerHydCond(nSoil,iGRU), mLayerDiffuse(nSoil,iGRU), &
 scalarDrainage(nSoil,iGRU), &
 dq_dHydStateUnsat(nSoil,iGRU), dq_dNrgStateUnsat(nSoil,iGRU))


end subroutine

  ! ***************************************************************************************************************
  ! private subroutine qDrainFlux: compute the drainage flux from the bottom of the soil profile and its derivative
  ! ***************************************************************************************************************
  attributes(device) subroutine qDrainFlux(bc_lower,ixRichards, &
 nodeMatricHeadLiq, nodeVolFracLiq, &
 nodeDepth, nodeHeight, &
 lowerBoundHead, lowerBoundTheta, &
 node_dPsi_dTheta, node_dPsiLiq_dTemp, &
 bottomSatHydCond, iceImpedeFac, nodeHydCond, surfaceSatHydCond, &
 dHydCond_dVolLiq, dHydCond_dMatric, dHydCond_dTemp, &
 vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, kAnisotropic, zScale_TOPMODEL, &
 bottomHydCond, bottomDiffuse, &
 scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat)
    implicit none
   ! input: model control
   integer(i4b),intent(in) :: bc_lower ! index defining the type of boundary conditions
   integer(i4b),intent(in) :: ixRichards          ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: nodeMatricHeadLiq   ! liquid matric head in the lowest unsaturated node (m)
   real(rkind),intent(in) :: nodeVolFracLiq      ! volumetric liquid water content in the lowest unsaturated node (-)
   ! input: model coordinate variables
   real(rkind),intent(in) :: nodeDepth           ! depth of the lowest unsaturated soil layer (m)
   real(rkind),intent(in) :: nodeHeight                ! height of the lowest unsaturated soil node (m)
   ! input: diriclet boundary conditions
   real(rkind),intent(in) :: lowerBoundHead      ! lower boundary condition for matric head (m)
   real(rkind),intent(in) :: lowerBoundTheta     ! lower boundary condition for volumetric liquid water content (-)
   ! input: derivative in soil water characteristic
   real(rkind),intent(in) :: node_dPsi_dTheta      ! derivative of the soil moisture characteristic w.r.t. theta (m)
   real(rkind),intent(in) :: node_dPsiLiq_dTemp    ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: transmittance
   real(rkind),intent(in) :: bottomSatHydCond    ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   real(rkind),intent(in) :: iceImpedeFac        ! ice impedence factor in the upper-most soil layer (-)
   real(rkind),intent(in) :: nodeHydCond           ! hydraulic conductivity at the node itself (m s-1)
   real(rkind),intent(in) :: surfaceSatHydCond   ! saturated hydraulic conductivity at the surface (m s-1)
   ! input: transmittance derivatives
   real(rkind),intent(in) :: dHydCond_dVolLiq     ! derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
   real(rkind),intent(in) :: dHydCond_dMatric     ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
   real(rkind),intent(in) :: dHydCond_dTemp      ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha           ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n               ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m               ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: theta_sat           ! soil porosity (-)
   real(rkind),intent(in) :: theta_res           ! soil residual volumetric water content (-)
   real(rkind),intent(in) :: kAnisotropic           ! anisotropy factor for lateral hydraulic conductivity (-)
   real(rkind),intent(in) :: zScale_TOPMODEL       ! scale factor for TOPMODEL-ish baseflow parameterization (m)
   ! output: hydraulic conductivity at the bottom of the unsaturated zone
   real(rkind),intent(inout) :: bottomHydCond       ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   real(rkind),intent(inout) :: bottomDiffuse       ! hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
   ! output: drainage flux from the bottom of the soil profile
   real(rkind),intent(inout) :: scalarDrainage      ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   real(rkind),intent(inout) :: dq_dHydStateUnsat   ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateUnsat   ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
    ! -----------------------------------------------------------------------------------------------------------------------------
    ! local variables
    real(rkind)                      :: zWater                  ! effective water table depth (m)
    real(rkind)                      :: nodePsi                 ! matric head in the lowest unsaturated node (m)
    real(rkind)                      :: cflux                   ! capillary flux (m s-1)
    ! error control
    logical(lgt)                     :: return_flag             ! flag for return statements
    ! -----------------------------------------------------------------------------------------------------------------------------
  
    !  call initialize_qDrainFlux
  
     call update_qDrainFlux(bc_lower,ixRichards, &
 nodeMatricHeadLiq, nodeVolFracLiq, &
 nodeDepth, nodeHeight, &
 lowerBoundHead, lowerBoundTheta, &
 node_dPsi_dTheta, node_dPsiLiq_dTemp, &
 bottomSatHydCond, iceImpedeFac, nodeHydCond, surfaceSatHydCond, &
 dHydCond_dVolLiq, dHydCond_dMatric, dHydCond_dTemp, &
 vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, kAnisotropic, zScale_TOPMODEL, &
 bottomHydCond, bottomDiffuse, &
 scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat, &
 cflux, zWater, nodePsi)!;   if (return_flag) return
  
    !  call finalize_qDrainFlux; if (return_flag) return
  
  end subroutine qDrainFlux
   
 
 attributes(device) subroutine update_qDrainFlux(bc_lower,ixRichards, &
 nodeMatricHeadLiq, nodeVolFracLiq, &
 nodeDepth, nodeHeight, &
 lowerBoundHead, lowerBoundTheta, &
 node_dPsi_dTheta, node_dPsiLiq_dTemp, &
 bottomSatHydCond, iceImpedeFac, nodeHydCond, surfaceSatHydCond, &
 dHydCond_dVolLiq, dHydCond_dMatric, dHydCond_dTemp, &
 vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, kAnisotropic, zScale_TOPMODEL, &
 bottomHydCond, bottomDiffuse, &
 scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat, &
 cflux, zWater, nodePsi)
  ! ** Update operations for qDrainFlux **
   ! input: model control
   integer(i4b),intent(in) :: bc_lower ! index defining the type of boundary conditions
   integer(i4b),intent(in) :: ixRichards          ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: nodeMatricHeadLiq   ! liquid matric head in the lowest unsaturated node (m)
   real(rkind),intent(in) :: nodeVolFracLiq      ! volumetric liquid water content in the lowest unsaturated node (-)
   ! input: model coordinate variables
   real(rkind),intent(in) :: nodeDepth           ! depth of the lowest unsaturated soil layer (m)
   real(rkind),intent(in) :: nodeHeight                ! height of the lowest unsaturated soil node (m)
   ! input: diriclet boundary conditions
   real(rkind),intent(in) :: lowerBoundHead      ! lower boundary condition for matric head (m)
   real(rkind),intent(in) :: lowerBoundTheta     ! lower boundary condition for volumetric liquid water content (-)
   ! input: derivative in soil water characteristic
   real(rkind),intent(in) :: node_dPsi_dTheta      ! derivative of the soil moisture characteristic w.r.t. theta (m)
   real(rkind),intent(in) :: node_dPsiLiq_dTemp    ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: transmittance
   real(rkind),intent(in) :: bottomSatHydCond    ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   real(rkind),intent(in) :: iceImpedeFac        ! ice impedence factor in the upper-most soil layer (-)
   real(rkind),intent(in) :: nodeHydCond           ! hydraulic conductivity at the node itself (m s-1)
   real(rkind),intent(in) :: surfaceSatHydCond   ! saturated hydraulic conductivity at the surface (m s-1)
   ! input: transmittance derivatives
   real(rkind),intent(in) :: dHydCond_dVolLiq     ! derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
   real(rkind),intent(in) :: dHydCond_dMatric     ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
   real(rkind),intent(in) :: dHydCond_dTemp      ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha           ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n               ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m               ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: theta_sat           ! soil porosity (-)
   real(rkind),intent(in) :: theta_res           ! soil residual volumetric water content (-)
   real(rkind),intent(in) :: kAnisotropic           ! anisotropy factor for lateral hydraulic conductivity (-)
   real(rkind),intent(in) :: zScale_TOPMODEL       ! scale factor for TOPMODEL-ish baseflow parameterization (m)
   ! output: hydraulic conductivity at the bottom of the unsaturated zone
   real(rkind),intent(inout) :: bottomHydCond       ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   real(rkind),intent(inout) :: bottomDiffuse       ! hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
   ! output: drainage flux from the bottom of the soil profile
   real(rkind),intent(inout) :: scalarDrainage      ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   real(rkind),intent(inout) :: dq_dHydStateUnsat   ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateUnsat   ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
  ! local variables (to qDrainFlux)
  real(rkind),intent(inout)                      :: cflux                   ! capillary flux (m s-1)
  real(rkind),intent(inout)                      :: zWater                  ! effective water table depth (m)
  real(rkind),intent(inout)                      :: nodePsi                 ! matric head in the lowest unsaturated node (m)

   ! determine lower boundary condition
   select case(bc_lower)
     case(prescribedHead) ! specified matric head value
       call update_qDrainFlux_prescribedHead(ixRichards, &
 nodeMatricHeadLiq, nodeVolFracLiq, &
 nodeDepth, &
 lowerBoundHead, lowerBoundTheta, &
 bottomSatHydCond, iceImpedeFac, &
 dHydCond_dTemp, &
 vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, &
 bottomHydCond, bottomDiffuse, &
 scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat, &
 cflux)
     case(funcBottomHead) ! specified matric head function
       call update_qDrainFlux_funcBottomHead(ixRichards, &
 nodeMatricHeadLiq, nodeVolFracLiq, &
 nodeHeight, &
 node_dPsi_dTheta, node_dPsiLiq_dTemp, &
 surfaceSatHydCond, &
 vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, kAnisotropic, zScale_TOPMODEL, &
 scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat, &
 zWater, nodePsi)
     case(freeDrainage)   ! free drainage 
       call update_qDrainFlux_freeDrainage(ixRichards, &
 nodeHydCond, &
 dHydCond_dVolLiq, dHydCond_dMatric, dHydCond_dTemp, &
 kAnisotropic, &
 scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat)
     case(zeroFlux)       ! zero flux
       call update_qDrainFlux_zeroFlux(scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat)
     case default;
        ! err=20; message=trim(message)//'unknown lower boundary condition for soil hydrology'; return_flag=.true.; return
   end select 

  end subroutine update_qDrainFlux

 attributes(device) subroutine update_qDrainFlux_prescribedHead(ixRichards, &
 nodeMatricHeadLiq, nodeVolFracLiq, &
 nodeDepth, &
 lowerBoundHead, lowerBoundTheta, &
 bottomSatHydCond, iceImpedeFac, &
 dHydCond_dTemp, &
 vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, &
 bottomHydCond, bottomDiffuse, &
 scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat, &
 cflux)
  USE soil_utils_module,only:hydCond_psi ! compute hydraulic conductivity as a function of matric head (m s-1)
  USE soil_utils_module,only:hydCond_liq ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
  USE soil_utils_module,only:dPsi_dTheta ! compute derivative of the soil moisture characteristic w.r.t. theta (m)

 
  ! ** Update operations for qDrainFlux: prescribed pressure head value at bottom boundary **
   ! input: model control
   integer(i4b),intent(in) :: ixRichards          ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: nodeMatricHeadLiq   ! liquid matric head in the lowest unsaturated node (m)
   real(rkind),intent(in) :: nodeVolFracLiq      ! volumetric liquid water content in the lowest unsaturated node (-)
   ! input: model coordinate variables
   real(rkind),intent(in) :: nodeDepth           ! depth of the lowest unsaturated soil layer (m)
   ! input: diriclet boundary conditions
   real(rkind),intent(in) :: lowerBoundHead      ! lower boundary condition for matric head (m)
   real(rkind),intent(in) :: lowerBoundTheta     ! lower boundary condition for volumetric liquid water content (-)
   ! input: transmittance
   real(rkind),intent(in) :: bottomSatHydCond    ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   real(rkind),intent(in) :: iceImpedeFac        ! ice impedence factor in the upper-most soil layer (-)
   ! input: transmittance derivatives
   real(rkind),intent(in) :: dHydCond_dTemp      ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha           ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n               ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m               ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: theta_sat           ! soil porosity (-)
   real(rkind),intent(in) :: theta_res           ! soil residual volumetric water content (-)
   ! output: hydraulic conductivity at the bottom of the unsaturated zone
   real(rkind),intent(inout) :: bottomHydCond       ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
   real(rkind),intent(inout) :: bottomDiffuse       ! hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
   ! output: drainage flux from the bottom of the soil profile
   real(rkind),intent(inout) :: scalarDrainage      ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   real(rkind),intent(inout) :: dq_dHydStateUnsat   ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateUnsat   ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
  ! local variables (to qDrainFlux)
  real(rkind),intent(inout)                      :: cflux                   ! capillary flux (m s-1)

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
    !  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select 
   scalarDrainage = cflux + bottomHydCond

   ! hydrology derivatives
   select case(ixRichards)  ! select form of Richards' equation
     case(moisture); dq_dHydStateUnsat = bottomDiffuse/(nodeDepth/2._rkind)
     case(mixdform); dq_dHydStateUnsat = bottomHydCond/(nodeDepth/2._rkind)
    !  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select
   ! energy derivatives
   dq_dNrgStateUnsat = -(dHydCond_dTemp/2._rkind)*(lowerBoundHead  - nodeMatricHeadLiq)/(nodeDepth*0.5_rkind)&
                     & + dHydCond_dTemp/2._rkind
 
 end subroutine update_qDrainFlux_prescribedHead

 attributes(device) subroutine update_qDrainFlux_funcBottomHead(ixRichards, &
 nodeMatricHeadLiq, nodeVolFracLiq, &
 nodeHeight, &
 node_dPsi_dTheta, node_dPsiLiq_dTemp, &
 surfaceSatHydCond, &
 vGn_alpha, vGn_n, vGn_m, theta_sat, theta_res, kAnisotropic, zScale_TOPMODEL, &
 scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat, &
 zWater, nodePsi)
  USE soil_utils_module,only:matricHead  ! compute matric head as a function of volumetric fraction of liquid water (m)

  ! ** Update operations for qDrainFlux: prescribed pressure head function at bottom boundary **
   ! input: model control
   integer(i4b),intent(in) :: ixRichards              ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: state and diagnostic variables
   real(rkind),intent(in) :: nodeMatricHeadLiq   ! liquid matric head in the lowest unsaturated node (m)
   real(rkind),intent(in) :: nodeVolFracLiq      ! volumetric liquid water content in the lowest unsaturated node (-)
   ! input: model coordinate variables
   real(rkind),intent(in) :: nodeHeight                ! height of the lowest unsaturated soil node (m)
   ! input: derivative in soil water characteristic
   real(rkind),intent(in) :: node_dPsi_dTheta      ! derivative of the soil moisture characteristic w.r.t. theta (m)
   real(rkind),intent(in) :: node_dPsiLiq_dTemp    ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
   ! input: transmittance
   real(rkind),intent(in) :: surfaceSatHydCond   ! saturated hydraulic conductivity at the surface (m s-1)
   ! input: soil parameters
   real(rkind),intent(in) :: vGn_alpha             ! van Genuchten "alpha" parameter (m-1)
   real(rkind),intent(in) :: vGn_n                 ! van Genuchten "n" parameter (-)
   real(rkind),intent(in) :: vGn_m                 ! van Genuchten "m" parameter (-)
   real(rkind),intent(in) :: theta_sat             ! soil porosity (-)
   real(rkind),intent(in) :: theta_res             ! soil residual volumetric water content (-)
   real(rkind),intent(in) :: kAnisotropic          ! anisotropy factor for lateral hydraulic conductivity (-)
   real(rkind),intent(in) :: zScale_TOPMODEL       ! scale factor for TOPMODEL-ish baseflow parameterization (m)
   ! output: drainage flux from the bottom of the soil profile
   real(rkind),intent(inout) :: scalarDrainage        ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   real(rkind),intent(inout) :: dq_dHydStateUnsat  ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateUnsat  ! ... energy state variable in lowest unsaturated node (m s-1 K-1)
   ! local variables
  real(rkind),intent(inout)                      :: zWater                  ! effective water table depth (m)
  real(rkind),intent(inout)                      :: nodePsi                 ! matric head in the lowest unsaturated node (m)
  

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
    !  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select
   ! energy derivatives
   dq_dNrgStateUnsat = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)*node_dPsiLiq_dTemp/zScale_TOPMODEL

 end subroutine update_qDrainFlux_funcBottomHead

 attributes(device) subroutine update_qDrainFlux_freeDrainage(ixRichards, &
 nodeHydCond, &
 dHydCond_dVolLiq, dHydCond_dMatric, dHydCond_dTemp, &
 kAnisotropic, &
 scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat)
  ! ** Update operations for qDrainFlux: free drainage at bottom boundary **
   ! input: model control
   integer(i4b),intent(in) :: ixRichards              ! index defining the option for Richards' equation (moisture or mixdform)
   ! input: transmittance
   real(rkind),intent(in) :: nodeHydCond           ! hydraulic conductivity at the node itself (m s-1)
   ! input: transmittance derivatives
   real(rkind),intent(in) :: dHydCond_dVolLiq     ! derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
   real(rkind),intent(in) :: dHydCond_dMatric     ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
   real(rkind),intent(in) :: dHydCond_dTemp       ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
   ! input: soil parameters
   real(rkind),intent(in) :: kAnisotropic           ! anisotropy factor for lateral hydraulic conductivity (-)
   ! output: drainage flux from the bottom of the soil profile
   real(rkind),intent(inout) :: scalarDrainage        ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   real(rkind),intent(inout) :: dq_dHydStateUnsat  ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateUnsat  ! ... energy state variable in lowest unsaturated node (m s-1 K-1)

   scalarDrainage = nodeHydCond*kAnisotropic ! compute flux

   ! hydrology derivatives
   select case(ixRichards)  ! select form of Richards' equation
     case(moisture); dq_dHydStateUnsat = dHydCond_dVolLiq*kAnisotropic
     case(mixdform); dq_dHydStateUnsat = dHydCond_dMatric*kAnisotropic
    !  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return_flag=.true.; return
   end select
   ! energy derivatives
   dq_dNrgStateUnsat = dHydCond_dTemp*kAnisotropic

 end subroutine update_qDrainFlux_freeDrainage

 attributes(device) subroutine update_qDrainFlux_zeroFlux(scalarDrainage, &
 dq_dHydStateUnsat, dq_dNrgStateUnsat)
  ! ** Update operations for qDrainFlux: zero flux condition at bottom boundary **
   ! output: drainage flux from the bottom of the soil profile
   real(rkind),intent(inout) :: scalarDrainage        ! drainage flux from the bottom of the soil profile (m s-1)
   ! output: derivatives in drainage flux w.r.t. ...
   real(rkind),intent(inout) :: dq_dHydStateUnsat  ! ... state variable in lowest unsaturated node (m s-1 or s-1)
   real(rkind),intent(inout) :: dq_dNrgStateUnsat  ! ... energy state variable in lowest unsaturated node (m s-1 K-1)

   scalarDrainage = 0._rkind
   dq_dHydStateUnsat = 0._rkind
   dq_dNrgStateUnsat = 0._rkind

 end subroutine update_qDrainFlux_zeroFlux

   attributes(device) function crit_soilT_wrapper(psi)
      USE soil_utils_module,only:crit_soilT            ! compute critical temperature below which ice exists

    real(rkind) :: psi, crit_soilT_wrapper
    crit_soilT_wrapper = crit_soilT(psi)
  end function


  end module soilLiqFlx_module
  