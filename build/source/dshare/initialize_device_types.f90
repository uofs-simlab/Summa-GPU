module initialize_device
  use data_types
  use device_data_types
  USE var_lookup,only:iLookFLUX        ! lookup indices for flux data
  USE var_lookup,only:iLookDERIV       ! lookup indices for derivative data
  USE var_lookup,only:iLookFORCE       ! lookup indices for forcing data 
  USE var_lookup,only:iLookDIAG        ! lookup indices for diagnostic variable data
  USE var_lookup,only:iLookPROG,iLookPARAM, iLookINDEX, iLookBVAR, iLookTYPE, iLookATTR
  ! USE globalData,only: nSpecBands                 ! number of spectral bands

  implicit none
  private
  integer(i4b),parameter :: nSpecBands = 2
  public::get_iGRU
  public::allocate_device_flux2state
  public::deallocate_device_flux_data,deallocate_device_deriv_data
  public::allocate_device_flux_prev,allocate_device_deriv_data,allocate_device_flux_count
  public::update_flux_sum
  public::allocate_device_diag_temp,deallocate_device_diag_data
  public::allocate_device_prog_temp,deallocate_device_prog_data
  public::allocate_device_flux_mask,zero_device_deriv_data
  public::allocate_device_decisions
  public::allocate_veg_parameters,allocate_veg_table,allocate_device_indx_data,allocate_device_attr_data
  public::deallocate_device_forc_data,allocate_device_forc_data
public::deallocate_device_indx_data,allocate_device_indx_temp

public::get_device_param_data
public::set_device_param_data_scalar
public::set_device_param_data
public::get_device_type_data
public::allocate_device_bvar_data,allocate_device_prog_data,allocate_device_diag_data
public::deallocate_device_prog_data
public::deallocate_device_bvar_data
public::deallocate_device_diag_data
public::deallocate_device_flux_data,allocate_device_flux_data
public::deallocate_device_param_data
public::deallocate_device_attr_data
public::deallocate_device_decisions
public::finalize_device_forc_data
public::finalize_device_indx_data
public::finalize_device_prog_data
public::finalize_device_bvar_data,allocate_device_forc_data
public::finalize_device_diag_data
public::finalize_device_flux_data
public::allocate_device_param_data

  contains

  subroutine allocate_device_param_data(mpar_data_device,mpar_data,nGRU)
    type(mpar_data_device), intent(inout) :: mpar_data_device
    type(gru_hru_doubleVec),intent(in) :: mpar_data
    integer(i4b) :: nGRU,iGRU
    iGRU = 1
    allocate(mpar_data_device%k_soil_(size(mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%k_soil)%dat),nGRU))
    allocate(mpar_data_device%k_macropore_(size(mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%k_macropore)%dat),nGRU))
    allocate(mpar_data_device%compactedDepth_(nGRU))
  allocate(mpar_data_device%soil_dens_intr_(size(mpar_data%gru(1)%hru(1)%var(iLookPARAM%soil_dens_intr)%dat),nGRU))
  allocate(mpar_data_device%thCond_soil_(size(mpar_data%gru(1)%hru(1)%var(iLookPARAM%thCond_soil)%dat),nGRU))
  allocate(mpar_data_device%frac_sand_(size(mpar_data%gru(1)%hru(1)%var(iLookPARAM%frac_sand)%dat),nGRU))
  allocate(mpar_data_device%frac_clay_(size(mpar_data%gru(1)%hru(1)%var(iLookPARAM%frac_clay)%dat),nGRU))
  allocate(mpar_data_device%frac_silt_(size(mpar_data%gru(1)%hru(1)%var(iLookPARAM%frac_silt)%dat),nGRU))
  allocate(mpar_data_device%snowfrz_scale_(nGRU))
  allocate(mpar_data_device%specificHeatVeg_(nGRU))
  allocate(mpar_data_device%maxMassVegetation_(nGRU))
  allocate(mpar_data_device%Fcapil_(nGRU))
  allocate(mpar_data_device%k_snow_(nGRU))
  allocate(mpar_data_device%mw_exp_(nGRU))
  allocate(mpar_data_device%fixedThermalCond_snow_(nGRU))
  allocate(mpar_data_device%rootScaleFactor1_(nGRU))
allocate(mpar_data_device%rootScaleFactor2_(nGRU))
allocate(mpar_data_device%rootDistExp_(nGRU))
    allocate(mpar_data_device%Louis79_bparam_(nGRU))
allocate(mpar_data_device%Louis79_cStar_(nGRU))
allocate(mpar_data_device%Mahrt87_eScale_(nGRU))
    allocate(mpar_data_device%aquiferBaseflowExp_(nGRU))
    allocate(mpar_data_device%aquiferBaseflowRate_(nGRU))
    allocate(mpar_data_device%aquiferScaleFactor_(nGRU))
    allocate(mpar_data_device%canopyDrainageCoeff_(nGRU))
    allocate(mpar_data_device%canopyWettingExp_(nGRU))
    allocate(mpar_data_device%canopyWettingFactor_(nGRU))
    allocate(mpar_data_device%critAquiferTranspire_(nGRU))
allocate(mpar_data_device%critRichNumber_(nGRU))
allocate(mpar_data_device%critSoilTranspire_(nGRU))
allocate(mpar_data_device%critSoilWilting_(nGRU))
allocate(mpar_data_device%f_impede_(nGRU))
    allocate(mpar_data_device%heightCanopyBottom_(nGRU))
    allocate(mpar_data_device%heightCanopyTop_(nGRU))
    allocate(mpar_data_device%kAnisotropic_(nGRU))
    allocate(mpar_data_device%leafDimension_(nGRU))
    allocate(mpar_data_device%leafExchangeCoeff_(nGRU))
    allocate(mpar_data_device%lowerBoundHead_(nGRU))
allocate(mpar_data_device%lowerBoundTemp_(nGRU))
allocate(mpar_data_device%lowerBoundTheta_(nGRU))
    allocate(mpar_data_device%minStomatalResistance_(nGRU))
    allocate(mpar_data_device%mpExp_(nGRU))
    allocate(mpar_data_device%plantWiltPsi_(nGRU))
    allocate(mpar_data_device%qSurfScale_(nGRU))
allocate(mpar_data_device%rootingDepth_(nGRU))
allocate(mpar_data_device%soilIceCV_(nGRU))
allocate(mpar_data_device%soilIceScale_(nGRU))
  do iGRU=1,nGRU
    mpar_data_device%k_soil_(:,iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%k_macropore)%dat
    mpar_data_device%k_macropore_(:,iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%k_macropore)%dat
    mpar_data_device%compactedDepth_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%compactedDepth)%dat(1)
    mpar_data_device%soil_dens_intr_(:,iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%soil_dens_intr)%dat
mpar_data_device%thCond_soil_(:,iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%thCond_soil)%dat
mpar_data_device%frac_sand_(:,iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%frac_sand)%dat
mpar_data_device%frac_clay_(:,iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%frac_clay)%dat
mpar_data_device%frac_silt_(:,iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%frac_silt)%dat
mpar_data_device%snowfrz_scale_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%snowfrz_scale)%dat(1)
mpar_data_device%specificHeatVeg_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%specificHeatVeg)%dat(1)
mpar_data_device%maxMassVegetation_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%maxMassVegetation)%dat(1)
mpar_data_device%Fcapil_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%Fcapil)%dat(1)
mpar_data_device%k_snow_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%k_snow)%dat(1)
mpar_data_device%mw_exp_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%mw_exp)%dat(1)
mpar_data_device%fixedThermalCond_snow_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%fixedThermalCond_snow)%dat(1)
mpar_data_device%rootScaleFactor1_(iGRU) = mpar_data%gru(1)%hru(1)%var(iLookPARAM%rootScaleFactor1)%dat(1)
mpar_data_device%rootScaleFactor2_(iGRU) = mpar_data%gru(1)%hru(1)%var(iLookPARAM%rootScaleFactor2)%dat(1)
mpar_data_device%rootDistExp_(iGRU) = mpar_data%gru(1)%hru(1)%var(iLookPARAM%rootDistExp)%dat(1)
    mpar_data_device%Louis79_bparam_(iGRU) = mpar_data%gru(1)%hru(1)%var(iLookPARAM%Louis79_bparam)%dat(1)
    mpar_data_device%Louis79_cStar_(iGRU) = mpar_data%gru(1)%hru(1)%var(iLookPARAM%Louis79_cStar)%dat(1)
    mpar_data_device%Mahrt87_eScale_(iGRU) = mpar_data%gru(1)%hru(1)%var(iLookPARAM%Mahrt87_eScale)%dat(1)
    mpar_data_device%aquiferBaseflowExp_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%aquiferBaseflowExp)%dat(1)
    mpar_data_device%aquiferBaseflowRate_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%aquiferBaseflowRate)%dat(1)
    mpar_data_device%aquiferScaleFactor_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%aquiferScaleFactor)%dat(1)
    mpar_data_device%canopyDrainageCoeff_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%canopyDrainageCoeff)%dat(1)
    mpar_data_device%canopyWettingExp_ = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%canopyWettingExp)%dat(1)
    mpar_data_device%canopyWettingFactor_ = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%canopyWettingFactor)%dat(1)
    mpar_data_device%critAquiferTranspire_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%critAquiferTranspire)%dat(1)
    mpar_data_device%critRichNumber_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%critRichNumber)%dat(1)
    mpar_data_device%critSoilTranspire_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%critSoilTranspire)%dat(1)
    mpar_data_device%critSoilWilting_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%critSoilWilting)%dat(1)
    mpar_data_device%f_impede_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%f_impede)%dat(1)
      mpar_data_device%heightCanopyBottom_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%heightCanopyBottom)%dat(1)
      mpar_data_device%heightCanopyTop_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%heightCanopyTop)%dat(1)
      mpar_data_device%kAnisotropic_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%kAnisotropic)%dat(1)
      mpar_data_device%leafDimension_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%leafDimension)%dat(1)
      mpar_data_device%leafExchangeCoeff_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%leafExchangeCoeff)%dat(1)
      mpar_data_device%lowerBoundHead_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%lowerBoundHead)%dat(1)
      mpar_data_device%lowerBoundTemp_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%lowerBoundTemp)%dat(1)
      mpar_data_device%lowerBoundTheta_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%lowerBoundTheta)%dat(1)
      mpar_data_device%minStomatalResistance_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%minStomatalResistance)%dat(1)
      mpar_data_device%mpExp_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%mpExp)%dat(1)
      mpar_data_device%plantWiltPsi_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%plantWiltPsi)%dat(1)
      mpar_data_device%qSurfScale_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%qSurfScale)%dat(1)
      mpar_data_device%rootingDepth_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%rootingDepth)%dat(1)
    mpar_data_device%soilIceCV_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%soilIceCV)%dat(1)
    mpar_data_device%soilIceScale_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%soilIceScale)%dat(1)
    end do
    iGRU = 1
    allocate(mpar_data_device%soilStressParam_(nGRU))
    allocate(mpar_data_device%specificStorage_(nGRU))
    allocate(mpar_data_device%theta_mp_(nGRU))
    allocate(mpar_data_device%theta_res_(size(mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%theta_res)%dat),nGRU))
    allocate(mpar_data_device%theta_sat_(size(mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%theta_res)%dat),nGRU))
    allocate(mpar_data_device%throughfallScaleRain_(nGRU))
    allocate(mpar_data_device%upperBoundHead_(nGRU))
    allocate(mpar_data_device%upperBoundTemp_(nGRU))
    allocate(mpar_data_device%upperBoundTheta_(nGRU))
    allocate(mpar_data_device%vGn_alpha_(size(mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%vGn_alpha)%dat),nGRU))
    allocate(mpar_data_device%vGn_n_(size(mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%vGn_alpha)%dat),nGRU))
    allocate(mpar_data_device%wettingFrontSuction_(nGRU))
    allocate(mpar_data_device%windReductionParam_(nGRU))
    allocate(mpar_data_device%z0Canopy_(nGRU))
    allocate(mpar_data_device%z0Snow_(nGRU))
    allocate(mpar_data_device%z0Soil_(nGRU))
    allocate(mpar_data_device%zScale_TOPMODEL_(nGRU))
    allocate(mpar_data_device%zpdFraction_(nGRU))
    allocate(mpar_data_device%Kc25_(nGRU))
allocate(mpar_data_device%Ko25_(nGRU))
allocate(mpar_data_device%Kc_qFac_(nGRU))
allocate(mpar_data_device%Ko_qFac_(nGRU))
allocate(mpar_data_device%kc_Ha_(nGRU))
allocate(mpar_data_device%ko_Ha_(nGRU))
allocate(mpar_data_device%vcmax25_canopyTop_(nGRU))
allocate(mpar_data_device%vcmax_qFac_(nGRU))
allocate(mpar_data_device%vcmax_Ha_(nGRU))
allocate(mpar_data_device%vcmax_Hd_(nGRU))
allocate(mpar_data_device%vcmax_Sv_(nGRU))
allocate(mpar_data_device%vcmax_Kn_(nGRU))
allocate(mpar_data_device%jmax25_scale_(nGRU))
allocate(mpar_data_device%jmax_Ha_(nGRU))
allocate(mpar_data_device%jmax_Hd_(nGRU))
allocate(mpar_data_device%jmax_Sv_(nGRU))
allocate(mpar_data_device%fractionJ_(nGRU))
allocate(mpar_data_device%quantamYield_(nGRU))
allocate(mpar_data_device%vpScaleFactor_(nGRU))
allocate(mpar_data_device%cond2photo_slope_(nGRU))
allocate(mpar_data_device%minStomatalConductance_(nGRU))
    do iGRU=1,nGRU
      mpar_data_device%soilStressParam_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%soilStressParam)%dat(1)
      mpar_data_device%specificStorage_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%specificStorage)%dat(1)
      mpar_data_device%theta_mp_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%theta_mp)%dat(1)
      mpar_data_device%theta_res_(:,iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%theta_res)%dat
      mpar_data_device%theta_sat_(:,iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%theta_sat)%dat
      mpar_data_device%throughfallScaleRain_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%throughfallScaleRain)%dat(1)
      mpar_data_device%upperBoundHead_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%upperBoundHead)%dat(1)
      mpar_data_device%upperBoundTemp_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%upperBoundTemp)%dat(1)
      mpar_data_device%upperBoundTheta_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%upperBoundTheta)%dat(1)
      mpar_data_device%vGn_alpha_(:,iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%vGn_alpha)%dat
      mpar_data_device%vGn_n_(:,iGRU)       = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%vGn_n      )%dat
      mpar_data_device%wettingFrontSuction_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%wettingFrontSuction)%dat(1)
      mpar_data_device%windReductionParam_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%windReductionParam)%dat(1)
    mpar_data_device%z0Canopy_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%z0Canopy)%dat(1)
    mpar_data_device%z0Snow_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%z0Snow)%dat(1)
    mpar_data_device%z0Soil_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%z0Soil)%dat(1)
    mpar_data_device%zScale_TOPMODEL_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zScale_TOPMODEL)%dat(1)
      mpar_data_device%zpdFraction_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zpdFraction)%dat(1)
    mpar_data_device%Kc25_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%Kc25)%dat(1)
    mpar_data_device%Ko25_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%Ko25)%dat(1)
    mpar_data_device%Kc_qFac_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%Kc_qFac)%dat(1)
    mpar_data_device%Ko_qFac_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%Ko_qFac)%dat(1)
    mpar_data_device%kc_Ha_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%kc_Ha)%dat(1)
    mpar_data_device%ko_Ha_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%ko_Ha)%dat(1)
    mpar_data_device%vcmax25_canopyTop_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%vcmax25_canopyTop)%dat(1)
    mpar_data_device%vcmax_qFac_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%vcmax_qFac)%dat(1)
    mpar_data_device%vcmax_Ha_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%vcmax_Ha)%dat(1)
    mpar_data_device%vcmax_Hd_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%vcmax_Hd)%dat(1)
    mpar_data_device%vcmax_Sv_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%vcmax_Sv)%dat(1)
    mpar_data_device%vcmax_Kn_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%vcmax_Kn)%dat(1)
    mpar_data_device%jmax25_scale_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%jmax25_scale)%dat(1)
    mpar_data_device%jmax_Ha_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%jmax_Ha)%dat(1)
    mpar_data_device%jmax_Hd_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%jmax_Hd)%dat(1)
    mpar_data_device%jmax_Sv_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%jmax_Sv)%dat(1)
    mpar_data_device%fractionJ_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%fractionJ)%dat(1)
    mpar_data_device%quantamYield_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%quantamYield)%dat(1)
    mpar_data_device%vpScaleFactor_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%vpScaleFactor)%dat(1)
    mpar_data_device%cond2photo_slope_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%cond2photo_slope)%dat(1)
    mpar_data_device%minStomatalConductance_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%minStomatalConductance)%dat(1)
    end do
    iGRU = 1
    mpar_data_device%idaMaxOrder = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%idaMaxOrder)%dat(1)
    mpar_data_device%idaMaxInternalSteps = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%idaMaxInternalSteps)%dat(1)
    mpar_data_device%idaInitStepSize = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%idaInitStepSize)%dat(1)
    mpar_data_device%idaMinStepSize = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%idaMinStepSize)%dat(1)
    mpar_data_device%idaMaxErrTestFail = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%idaMaxErrTestFail)%dat(1)
    allocate(mpar_data_device%fieldCapacity_(nGRU))
    do iGRU=1,nGRU
    mpar_data_device%fieldCapacity_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%fieldCapacity)%dat(1)
    end do
    iGRU = 1
    mpar_data_device%absTolTempCas = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%absTolTempCas)%dat(1)
    mpar_data_device%relTolTempCas = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%relTolTempCas)%dat(1)
    mpar_data_device%absTolTempVeg = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%absTolTempVeg)%dat(1)
    mpar_data_device%relTolTempVeg = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%relTolTempVeg)%dat(1)
    mpar_data_device%absTolWatVeg = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%absTolWatVeg)%dat(1)
    mpar_data_device%relTolWatVeg = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%relTolWatVeg)%dat(1)
    mpar_data_device%absTolTempSoilSnow = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%absTolTempSoilSnow)%dat(1)
    mpar_data_device%relTolTempSoilSnow = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%relTolTempSoilSnow)%dat(1)
    mpar_data_device%absTolWatSnow = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%absTolWatSnow)%dat(1)
    mpar_data_device%relTolWatSnow = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%relTolWatSnow)%dat(1)
    mpar_data_device%absTolMatric = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%absTolMatric)%dat(1)
    mpar_data_device%relTolMatric = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%relTolMatric)%dat(1)
    mpar_data_device%absTolAquifr = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%absTolAquifr)%dat(1)
    mpar_data_device%relTolAquifr = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%relTolAquifr)%dat(1)
    mpar_data_device%maxstep = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%maxstep)%dat(1)
    mpar_data_device%be_steps = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%be_steps)%dat(1)
    mpar_data_device%minstep = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%minstep)%dat(1)
    mpar_data_device%maxstep = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%maxstep)%dat(1)
    allocate(mpar_data_device%refInterceptCapRain_(nGRU))
allocate(mpar_data_device%refInterceptCapSnow_(nGRU))
allocate(mpar_data_device%Frad_vis_(nGRU))
allocate(mpar_data_device%Frad_direct_(nGRU))
allocate(mpar_data_device%directScale_(nGRU))
    allocate(mpar_data_device%albedoMax_(nGRU))
    allocate(mpar_data_device%albedoMinWinter_(nGRU))
    allocate(mpar_data_device%albedoMinSpring_(nGRU))
    do iGRU=1,nGRU
    mpar_data_device%refInterceptCapRain_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%refInterceptCapRain)%dat(1)
    mpar_data_device%refInterceptCapSnow_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%refInterceptCapSnow)%dat(1)
    mpar_data_device%Frad_vis_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%Frad_vis)%dat(1)
    mpar_data_device%Frad_direct_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%Frad_direct)%dat(1)
    mpar_data_device%directScale_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%directScale)%dat(1)
    mpar_data_device%albedoMax_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%albedoMax)%dat(1)
    mpar_data_device%albedoMinWinter_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%albedoMinWinter)%dat(1)
    mpar_data_device%albedoMinSpring_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%albedoMinSpring)%dat(1)
    end do
    iGRU = 1
    allocate(mpar_data_device%albedoMaxVisible_(nGRU))
    allocate(mpar_data_device%albedoMinVisible_(nGRU))
    allocate(mpar_data_device%albedoMaxNearIR_(nGRU))
    allocate(mpar_data_device%albedoMinNearIR_(nGRU))
    do iGRU=1,nGRU
    mpar_data_device%albedoMaxVisible_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%albedoMaxVisible)%dat(1)
    mpar_data_device%albedoMinVisible_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%albedoMinVisible)%dat(1)
    mpar_data_device%albedoMaxNearIR_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%albedoMaxNearIR)%dat(1)
    mpar_data_device%albedoMinNearIR_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%albedoMinNearIR)%dat(1)
    end do
    allocate(mpar_data_device%albedoDecayRate_(nGRU))
allocate(mpar_data_device%tempScalGrowth_(nGRU))
allocate(mpar_data_device%albedoSootLoad_(nGRU))
allocate(mpar_data_device%albedoRefresh_(nGRU))
allocate(mpar_data_device%snowUnloadingCoeff_(nGRU))
allocate(mpar_data_device%ratioDrip2Unloading_(nGRU))
allocate(mpar_data_device%minTempUnloading_(nGRU))
allocate(mpar_data_device%rateTempUnloading_(nGRU))
allocate(mpar_data_device%minWindUnloading_(nGRU))
allocate(mpar_data_device%rateWindUnloading_(nGRU))
allocate(mpar_data_device%densScalGrowth_(nGRU))
allocate(mpar_data_device%grainGrowthRate_(nGRU))
allocate(mpar_data_device%densScalOvrbdn_(nGRU))
allocate(mpar_data_device%tempScalOvrbdn_(nGRU))
allocate(mpar_data_device%baseViscosity_(nGRU))
allocate(mpar_data_device%zmax_(nGRU))
allocate(mpar_data_device%zmin_(nGRU))
allocate(mpar_data_device%zminLayer1_(nGRU))
allocate(mpar_data_device%zminLayer2_(nGRU))
allocate(mpar_data_device%zminLayer3_(nGRU))
allocate(mpar_data_device%zminLayer4_(nGRU))
allocate(mpar_data_device%zminLayer5_(nGRU))
allocate(mpar_data_device%zmaxLayer1_lower_(nGRU))
allocate(mpar_data_device%zmaxLayer2_lower_(nGRU))
allocate(mpar_data_device%zmaxLayer3_lower_(nGRU))
allocate(mpar_data_device%zmaxLayer4_lower_(nGRU))
allocate(mpar_data_device%zmaxLayer1_upper_(nGRU))
allocate(mpar_data_device%zmaxLayer2_upper_(nGRU))
allocate(mpar_data_device%zmaxLayer3_upper_(nGRU))
allocate(mpar_data_device%zmaxLayer4_upper_(nGRU))
    iGRU = 1
    do iGRU=1,nGRU
    mpar_data_device%albedoDecayRate_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%albedoDecayRate)%dat(1)
    mpar_data_device%tempScalGrowth_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%tempScalGrowth)%dat(1)
    mpar_data_device%albedoSootLoad_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%albedoSootLoad)%dat(1)
    mpar_data_device%albedoRefresh_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%albedoRefresh)%dat(1)
    mpar_data_device%snowUnloadingCoeff_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%snowUnloadingCoeff)%dat(1)
    mpar_data_device%ratioDrip2Unloading_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%ratioDrip2Unloading)%dat(1)
    mpar_data_device%minTempUnloading_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%minTempUnloading)%dat(1)
    mpar_data_device%rateTempUnloading_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%rateTempUnloading)%dat(1)
    mpar_data_device%minWindUnloading_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%minWindUnloading)%dat(1)
    mpar_data_device%rateWindUnloading_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%rateWindUnloading)%dat(1)
    mpar_data_device%densScalGrowth_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%densScalGrowth)%dat(1)
    mpar_data_device%grainGrowthRate_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%grainGrowthRate)%dat(1)
    mpar_data_device%densScalOvrbdn_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%densScalOvrbdn)%dat(1)
    mpar_data_device%tempScalOvrbdn_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%tempScalOvrbdn)%dat(1)
    mpar_data_device%baseViscosity_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%baseViscosity)%dat(1)
    mpar_data_device%zmax_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zmax)%dat(1)
    mpar_data_device%zmin_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zmin)%dat(1)
mpar_data_device%zminLayer1_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zminLayer1)%dat(1)
mpar_data_device%zminLayer2_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zminLayer2)%dat(1)
mpar_data_device%zminLayer3_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zminLayer3)%dat(1)
mpar_data_device%zminLayer4_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zminLayer4)%dat(1)
mpar_data_device%zminLayer5_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zminLayer5)%dat(1)
mpar_data_device%zmaxLayer1_lower_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zmaxLayer1_lower)%dat(1)
mpar_data_device%zmaxLayer2_lower_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zmaxLayer2_lower)%dat(1)
mpar_data_device%zmaxLayer3_lower_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zmaxLayer3_lower)%dat(1)
mpar_data_device%zmaxLayer4_lower_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zmaxLayer4_lower)%dat(1)
mpar_data_device%zmaxLayer1_upper_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zmaxLayer1_upper)%dat(1)
mpar_data_device%zmaxLayer2_upper_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zmaxLayer2_upper)%dat(1)
mpar_data_device%zmaxLayer3_upper_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zmaxLayer3_upper)%dat(1)
mpar_data_device%zmaxLayer4_upper_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%zmaxLayer4_upper)%dat(1)
    end do
    allocate(mpar_data_device%newSnowDenMin_(nGRU))
allocate(mpar_data_device%newSnowDenMult_(nGRU))
allocate(mpar_data_device%newSnowDenScal_(nGRU))
allocate(mpar_data_device%constSnowDen_(nGRU))
allocate(mpar_data_device%newSnowDenAdd_(nGRU))
allocate(mpar_data_device%newSnowDenMultTemp_(nGRU))
allocate(mpar_data_device%newSnowDenMultWind_(nGRU))
allocate(mpar_data_device%newSnowDenMultAnd_(nGRU))
allocate(mpar_data_device%newSnowDenBase_(nGRU))
allocate(mpar_data_device%tempRangeTimestep_(nGRU))
allocate(mpar_data_device%tempCritRain_(nGRU))
allocate(mpar_data_device%frozenPrecipMultip_(nGRU))
allocate(mpar_data_device%minwind_(nGRU))
allocate(mpar_data_device%winterSAI_(nGRU))
allocate(mpar_data_device%summerLAI_(nGRU))
do iGRU=1,nGRU
mpar_data_device%newSnowDenMin_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%newSnowDenMin)%dat(1)
mpar_data_device%newSnowDenMult_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%newSnowDenMult)%dat(1)
mpar_data_device%newSnowDenScal_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%newSnowDenScal)%dat(1)
mpar_data_device%constSnowDen_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%constSnowDen)%dat(1)
mpar_data_device%newSnowDenAdd_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%newSnowDenAdd)%dat(1)
mpar_data_device%newSnowDenMultTemp_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%newSnowDenMultTemp)%dat(1)
mpar_data_device%newSnowDenMultWind_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%newSnowDenMultWind)%dat(1)
mpar_data_device%newSnowDenMultAnd_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%newSnowDenMultAnd)%dat(1)
mpar_data_device%newSnowDenBase_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%newSnowDenBase)%dat(1)
mpar_data_device%tempRangeTimestep_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%tempRangeTimestep)%dat(1)
mpar_data_device%tempCritRain_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%tempCritRain)%dat(1)
mpar_data_device%frozenPrecipMultip_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%frozenPrecipMultip)%dat(1)
mpar_data_device%minwind_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%minwind)%dat(1)
mpar_data_device%winterSAI_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%winterSAI)%dat(1)
mpar_data_device%summerLAI_(iGRU) = mpar_data%gru(iGRU)%hru(1)%var(iLookPARAM%summerLAI)%dat(1)
end do
  end subroutine allocate_device_param_data

  subroutine set_device_param_data_scalar(mparStruct, iGRU, ixParam, parVector)
    implicit none
    type(mpar_data_device) :: mparStruct
    integer(i4b) :: iGRU, ixParam
    real(rkind) :: parVector

    select case(ixParam)
    case(iLookPARAM%upperBoundHead); print*, size(mparStruct%upperBoundHead_); mparStruct%upperBoundHead_(iGRU) = parVector
    case(iLookPARAM%lowerBoundHead); mparStruct%lowerBoundHead_(iGRU) = parVector
    case(iLookPARAM%upperBoundTheta); mparStruct%upperBoundTheta_(iGRU) = parVector
    case(iLookPARAM%lowerBoundTheta); mparStruct%lowerBoundTheta_(iGRU) = parVector
    case(iLookPARAM%upperBoundTemp); mparStruct%upperBoundTemp_(iGRU) = parVector
    case(iLookPARAM%lowerBoundTemp); mparStruct%lowerBoundTemp_(iGRU) = parVector
    case(iLookPARAM%tempCritRain); mparStruct%tempCritRain_(iGRU) = parVector
    case(iLookPARAM%tempRangeTimestep); mparStruct%tempRangeTimestep_(iGRU) = parVector
    case(iLookPARAM%frozenPrecipMultip); mparStruct%frozenPrecipMultip_(iGRU) = parVector
    case(iLookPARAM%snowfrz_scale); mparStruct%snowfrz_scale_(iGRU) = parVector
    case(iLookPARAM%fixedThermalCond_snow); mparStruct%fixedThermalCond_snow_(iGRU) = parVector
    case(iLookPARAM%albedoMax); mparStruct%albedoMax_(iGRU) = parVector
    case(iLookPARAM%albedoMinWinter); mparStruct%albedoMinWinter_(iGRU) = parVector
    case(iLookPARAM%albedoMinSpring); mparStruct%albedoMinSpring_(iGRU) = parVector
    case(iLookPARAM%albedoMaxVisible); mparStruct%albedoMaxVisible_(iGRU) = parVector
    case(iLookPARAM%albedoMinVisible); mparStruct%albedoMinVisible_(iGRU) = parVector
    case(iLookPARAM%albedoMaxNearIR); mparStruct%albedoMaxNearIR_(iGRU) = parVector
    case(iLookPARAM%albedoMinNearIR); mparStruct%albedoMinNearIR_(iGRU) = parVector
    case(iLookPARAM%albedoDecayRate); mparStruct%albedoDecayRate_(iGRU) = parVector
    case(iLookPARAM%albedoSootLoad); mparStruct%albedoSootLoad_(iGRU) = parVector
    case(iLookPARAM%albedoRefresh); mparStruct%albedoRefresh_(iGRU) = parVector
    case(iLookPARAM%directScale); mparStruct%directScale_(iGRU) = parVector
    case(iLookPARAM%Frad_direct); mparStruct%Frad_direct_(iGRU) = parVector
    case(iLookPARAM%Frad_vis); mparStruct%Frad_vis_(iGRU) = parVector
    case(iLookPARAM%newSnowDenMin); mparStruct%newSnowDenMin_(iGRU) = parVector
    case(iLookPARAM%newSnowDenMult); mparStruct%newSnowDenMult_(iGRU) = parVector
    case(iLookPARAM%newSnowDenScal); mparStruct%newSnowDenScal_(iGRU) = parVector
    case(iLookPARAM%constSnowDen); mparStruct%constSnowDen_(iGRU) = parVector
    case(iLookPARAM%newSnowDenAdd); mparStruct%newSnowDenAdd_(iGRU) = parVector
    case(iLookPARAM%newSnowDenMultTemp); mparStruct%newSnowDenMultTemp_(iGRU) = parVector
    case(iLookPARAM%newSnowDenMultWind); mparStruct%newSnowDenMultWind_(iGRU) = parVector
    case(iLookPARAM%newSnowDenMultAnd); mparStruct%newSnowDenMultAnd_(iGRU) = parVector
    case(iLookPARAM%newSnowDenBase); mparStruct%newSnowDenBase_(iGRU) = parVector
    case(iLookPARAM%densScalGrowth); mparStruct%densScalGrowth_(iGRU) = parVector
    case(iLookPARAM%tempScalGrowth); mparStruct%tempScalGrowth_(iGRU) = parVector
    case(iLookPARAM%grainGrowthRate); mparStruct%grainGrowthRate_(iGRU) = parVector
    case(iLookPARAM%densScalOvrbdn); mparStruct%densScalOvrbdn_(iGRU) = parVector
    case(iLookPARAM%tempScalOvrbdn); mparStruct%tempScalOvrbdn_(iGRU) = parVector
    case(iLookPARAM%baseViscosity); mparStruct%baseViscosity_(iGRU) = parVector
    case(iLookPARAM%Fcapil); mparStruct%Fcapil_(iGRU) = parVector
    case(iLookPARAM%k_snow); mparStruct%k_snow_(iGRU) = parVector
    case(iLookPARAM%mw_exp); mparStruct%mw_exp_(iGRU) = parVector
    case(iLookPARAM%z0Snow); mparStruct%z0Snow_(iGRU) = parVector
    case(iLookPARAM%z0Soil); mparStruct%z0Soil_(iGRU) = parVector
    case(iLookPARAM%z0Canopy); mparStruct%z0Canopy_(iGRU) = parVector
    case(iLookPARAM%zpdFraction); mparStruct%zpdFraction_(iGRU) = parVector
    case(iLookPARAM%critRichNumber); mparStruct%critRichNumber_(iGRU) = parVector
    case(iLookPARAM%Louis79_bparam); mparStruct%Louis79_bparam_(iGRU) = parVector
    case(iLookPARAM%Louis79_cStar); mparStruct%Louis79_cStar_(iGRU) = parVector
    case(iLookPARAM%Mahrt87_eScale); mparStruct%Mahrt87_eScale_(iGRU) = parVector
    case(iLookPARAM%leafExchangeCoeff); mparStruct%leafExchangeCoeff_(iGRU) = parVector
    case(iLookPARAM%windReductionParam); mparStruct%windReductionParam_(iGRU) = parVector
    case(iLookPARAM%Kc25); mparStruct%Kc25_(iGRU) = parVector
    case(iLookPARAM%Ko25); mparStruct%Ko25_(iGRU) = parVector
    case(iLookPARAM%Kc_qFac); mparStruct%Kc_qFac_(iGRU) = parVector
    case(iLookPARAM%Ko_qFac); mparStruct%Ko_qFac_(iGRU) = parVector
    case(iLookPARAM%kc_Ha); mparStruct%kc_Ha_(iGRU) = parVector
    case(iLookPARAM%ko_Ha); mparStruct%ko_Ha_(iGRU) = parVector
    case(iLookPARAM%vcmax25_canopyTop); mparStruct%vcmax25_canopyTop_(iGRU) = parVector
    case(iLookPARAM%vcmax_qFac); mparStruct%vcmax_qFac_(iGRU) = parVector
    case(iLookPARAM%vcmax_Ha); mparStruct%vcmax_Ha_(iGRU) = parVector
    case(iLookPARAM%vcmax_Hd); mparStruct%vcmax_Hd_(iGRU) = parVector
    case(iLookPARAM%vcmax_Sv); mparStruct%vcmax_Sv_(iGRU) = parVector
    case(iLookPARAM%vcmax_Kn); mparStruct%vcmax_Kn_(iGRU) = parVector
    case(iLookPARAM%jmax25_scale); mparStruct%jmax25_scale_(iGRU) = parVector
    case(iLookPARAM%jmax_Ha); mparStruct%jmax_Ha_(iGRU) = parVector
    case(iLookPARAM%jmax_Hd); mparStruct%jmax_Hd_(iGRU) = parVector
    case(iLookPARAM%jmax_Sv); mparStruct%jmax_Sv_(iGRU) = parVector
    case(iLookPARAM%fractionJ); mparStruct%fractionJ_(iGRU) = parVector
    case(iLookPARAM%quantamYield); mparStruct%quantamYield_(iGRU) = parVector
    case(iLookPARAM%vpScaleFactor); mparStruct%vpScaleFactor_(iGRU) = parVector
    case(iLookPARAM%cond2photo_slope); mparStruct%cond2photo_slope_(iGRU) = parVector
    case(iLookPARAM%minStomatalConductance); mparStruct%minStomatalConductance_(iGRU) = parVector
    case(iLookPARAM%winterSAI); mparStruct%winterSAI_(iGRU) = parVector
    case(iLookPARAM%summerLAI); mparStruct%summerLAI_(iGRU) = parVector
    case(iLookPARAM%rootScaleFactor1); mparStruct%rootScaleFactor1_(iGRU) = parVector
    case(iLookPARAM%rootScaleFactor2); mparStruct%rootScaleFactor2_(iGRU) = parVector
    case(iLookPARAM%rootingDepth); mparStruct%rootingDepth_(iGRU) = parVector
    case(iLookPARAM%rootDistExp); mparStruct%rootDistExp_(iGRU) = parVector
    case(iLookPARAM%plantWiltPsi); mparStruct%plantWiltPsi_(iGRU) = parVector
    case(iLookPARAM%soilStressParam); mparStruct%soilStressParam_(iGRU) = parVector
    case(iLookPARAM%critSoilWilting); mparStruct%critSoilWilting_(iGRU) = parVector
    case(iLookPARAM%critSoilTranspire); mparStruct%critSoilTranspire_(iGRU) = parVector
    case(iLookPARAM%critAquiferTranspire); mparStruct%critAquiferTranspire_(iGRU) = parVector
    case(iLookPARAM%minStomatalResistance); mparStruct%minStomatalResistance_(iGRU) = parVector
    case(iLookPARAM%leafDimension); mparStruct%leafDimension_(iGRU) = parVector
    case(iLookPARAM%heightCanopyTop); mparStruct%heightCanopyTop_(iGRU) = parVector
    case(iLookPARAM%heightCanopyBottom); mparStruct%heightCanopyBottom_(iGRU) = parVector
    case(iLookPARAM%specificHeatVeg); mparStruct%specificHeatVeg_(iGRU) = parVector
    case(iLookPARAM%maxMassVegetation); mparStruct%maxMassVegetation_(iGRU) = parVector
    case(iLookPARAM%throughfallScaleRain); mparStruct%throughfallScaleRain_(iGRU) = parVector
    case(iLookPARAM%refInterceptCapSnow); mparStruct%refInterceptCapSnow_(iGRU) = parVector
    case(iLookPARAM%refInterceptCapRain); mparStruct%refInterceptCapRain_(iGRU) = parVector
    case(iLookPARAM%snowUnloadingCoeff); mparStruct%snowUnloadingCoeff_(iGRU) = parVector
    case(iLookPARAM%canopyDrainageCoeff); mparStruct%canopyDrainageCoeff_(iGRU) = parVector
    case(iLookPARAM%ratioDrip2Unloading); mparStruct%ratioDrip2Unloading_(iGRU) = parVector
    case(iLookPARAM%canopyWettingFactor); mparStruct%canopyWettingFactor_(iGRU) = parVector
    case(iLookPARAM%canopyWettingExp); mparStruct%canopyWettingExp_(iGRU) = parVector
    case(iLookPARAM%minTempUnloading); mparStruct%minTempUnloading_(iGRU) = parVector
    case(iLookPARAM%rateTempUnloading); mparStruct%rateTempUnloading_(iGRU) = parVector
    case(iLookPARAM%minWindUnloading); mparStruct%minWindUnloading_(iGRU) = parVector
    case(iLookPARAM%rateWindUnloading); mparStruct%rateWindUnloading_(iGRU) = parVector
    case(iLookPARAM%soil_dens_intr); mparStruct%soil_dens_intr_(:,iGRU) = parVector
    case(iLookPARAM%thCond_soil); mparStruct%thCond_soil_(:,iGRU) = parVector
    case(iLookPARAM%frac_sand); mparStruct%frac_sand_(:,iGRU) = parVector
    case(iLookPARAM%frac_silt); mparStruct%frac_silt_(:,iGRU) = parVector
    case(iLookPARAM%frac_clay); mparStruct%frac_clay_(:,iGRU) = parVector
    case(iLookPARAM%theta_sat); mparStruct%theta_sat_(:,iGRU) = parVector
    case(iLookPARAM%theta_res); mparStruct%theta_res_(:,iGRU) = parVector
    case(iLookPARAM%vGn_alpha); mparStruct%vGn_alpha_(:,iGRU) = parVector
    case(iLookPARAM%vGn_n); mparStruct%vGn_n_(:,iGRU) = parVector
    case(iLookPARAM%k_soil); mparStruct%k_soil_(:,iGRU) = parVector
    case(iLookPARAM%k_macropore); mparStruct%k_macropore_(:,iGRU) = parVector
    case(iLookPARAM%fieldCapacity); mparStruct%fieldCapacity_(iGRU) = parVector
    case(iLookPARAM%wettingFrontSuction); mparStruct%wettingFrontSuction_(iGRU) = parVector
    case(iLookPARAM%theta_mp); mparStruct%theta_mp_(iGRU) = parVector
    case(iLookPARAM%mpExp); mparStruct%mpExp_(iGRU) = parVector
    case(iLookPARAM%kAnisotropic); mparStruct%kAnisotropic_(iGRU) = parVector
    case(iLookPARAM%zScale_TOPMODEL); mparStruct%zScale_TOPMODEL_(iGRU) = parVector
    case(iLookPARAM%compactedDepth); mparStruct%compactedDepth_(iGRU) = parVector
    case(iLookPARAM%aquiferBaseflowRate); mparStruct%aquiferBaseflowRate_(iGRU) = parVector
    case(iLookPARAM%aquiferScaleFactor); mparStruct%aquiferScaleFactor_(iGRU) = parVector
    case(iLookPARAM%aquiferBaseflowExp); mparStruct%aquiferBaseflowExp_(iGRU) = parVector
    case(iLookPARAM%qSurfScale); mparStruct%qSurfScale_(iGRU) = parVector
    ! case(iLookPARAM%specificYield); mparStruct%specificYield_(iGRU) = parVector
    case(iLookPARAM%specificStorage); mparStruct%specificStorage_(iGRU) = parVector
    case(iLookPARAM%f_impede); mparStruct%f_impede_(iGRU) = parVector
    case(iLookPARAM%soilIceScale); mparStruct%soilIceScale_(iGRU) = parVector
    case(iLookPARAM%soilIceCV); mparStruct%soilIceCV_(iGRU) = parVector
    case(iLookPARAM%minwind); mparStruct%minwind_(iGRU) = parVector
    case(iLookPARAM%minstep); mparStruct%minstep = parVector
    case(iLookPARAM%maxstep); mparStruct%maxstep = parVector
    case(iLookPARAM%be_steps); mparStruct%be_steps = parVector
    case(iLookPARAM%relTolTempCas); mparStruct%relTolTempCas = parVector
    case(iLookPARAM%absTolTempCas); mparStruct%absTolTempCas = parVector
    case(iLookPARAM%relTolTempVeg); mparStruct%relTolTempVeg = parVector
    case(iLookPARAM%absTolTempVeg); mparStruct%absTolTempVeg = parVector
    case(iLookPARAM%relTolWatVeg); mparStruct%relTolWatVeg = parVector
    case(iLookPARAM%absTolWatVeg); mparStruct%absTolWatVeg = parVector
    case(iLookPARAM%relTolTempSoilSnow); mparStruct%relTolTempSoilSnow = parVector
    case(iLookPARAM%absTolTempSoilSnow); mparStruct%absTolTempSoilSnow = parVector
    case(iLookPARAM%relTolWatSnow); mparStruct%relTolWatSnow = parVector
    case(iLookPARAM%absTolWatSnow); mparStruct%absTolWatSnow = parVector
    case(iLookPARAM%relTolMatric); mparStruct%relTolMatric = parVector
    case(iLookPARAM%absTolMatric); mparStruct%absTolMatric = parVector
    case(iLookPARAM%relTolAquifr); mparStruct%relTolAquifr = parVector
    case(iLookPARAM%absTolAquifr); mparStruct%absTolAquifr = parVector
    case(iLookPARAM%idaMaxOrder); mparStruct%idaMaxOrder = parVector
    case(iLookPARAM%idaMaxInternalSteps); mparStruct%idaMaxInternalSteps = parVector
    case(iLookPARAM%idaInitStepSize); mparStruct%idaInitStepSize = parVector
    case(iLookPARAM%idaMinStepSize); mparStruct%idaMinStepSize = parVector
    case(iLookPARAM%idaMaxErrTestFail); mparStruct%idaMaxErrTestFail = parVector
    case(iLookPARAM%zmin); mparStruct%zmin_(iGRU) = parVector
    case(iLookPARAM%zmax); mparStruct%zmax_(iGRU) = parVector
    case(iLookPARAM%zminLayer1); mparStruct%zminLayer1_(iGRU) = parVector
    case(iLookPARAM%zminLayer2); mparStruct%zminLayer2_(iGRU) = parVector
    case(iLookPARAM%zminLayer3); mparStruct%zminLayer3_(iGRU) = parVector
    case(iLookPARAM%zminLayer4); mparStruct%zminLayer4_(iGRU) = parVector
    case(iLookPARAM%zminLayer5); mparStruct%zminLayer5_(iGRU) = parVector
    case(iLookPARAM%zmaxLayer1_lower); mparStruct%zmaxLayer1_lower_(iGRU) = parVector
    case(iLookPARAM%zmaxLayer2_lower); mparStruct%zmaxLayer2_lower_(iGRU) = parVector
    case(iLookPARAM%zmaxLayer3_lower); mparStruct%zmaxLayer3_lower_(iGRU) = parVector
    case(iLookPARAM%zmaxLayer4_lower); mparStruct%zmaxLayer4_lower_(iGRU) = parVector
    case(iLookPARAM%zmaxLayer1_upper); mparStruct%zmaxLayer1_upper_(iGRU) = parVector
    case(iLookPARAM%zmaxLayer2_upper); mparStruct%zmaxLayer2_upper_(iGRU) = parVector
    case(iLookPARAM%zmaxLayer3_upper); mparStruct%zmaxLayer3_upper_(iGRU) = parVector
    case(iLookPARAM%zmaxLayer4_upper); mparStruct%zmaxLayer4_upper_(iGRU) = parVector
    case(iLookPARAM%FUSE_Ac_max); mparStruct%FUSE_Ac_max_(iGRU) = parVector
    case(iLookPARAM%FUSE_phi_tens); mparStruct%FUSE_phi_tens_(iGRU) = parVector
    case(iLookPARAM%FUSE_b); mparStruct%FUSE_b_(iGRU) = parVector
    case(iLookPARAM%FUSE_lambda); mparStruct%FUSE_lambda_(iGRU) = parVector
    case(iLookPARAM%FUSE_chi); mparStruct%FUSE_chi_(iGRU) = parVector
    case(iLookPARAM%FUSE_mu); mparStruct%FUSE_mu_(iGRU) = parVector
    case(iLookPARAM%FUSE_n); mparStruct%FUSE_n_(iGRU) = parVector
    end select
  end subroutine

  subroutine get_device_param_data(mparStruct, iVar, iGRU, parVector)
    implicit none
    type(mpar_data_device) :: mparStruct
    integer(i4b) :: iGRU, iVar
    real(rkind),allocatable :: parVector(:)

    if (allocated(parVector)) deallocate(parVector)
    select case(iVar)
    case(iLookPARAM%upperBoundHead); allocate(parVector(1)); parVector = mparStruct%upperBoundHead_(iGRU)
    case(iLookPARAM%lowerBoundHead); allocate(parVector(1)); parVector = mparStruct%lowerBoundHead_(iGRU)
    case(iLookPARAM%upperBoundTheta); allocate(parVector(1)); parVector = mparStruct%upperBoundTheta_(iGRU)
    case(iLookPARAM%lowerBoundTheta); allocate(parVector(1)); parVector = mparStruct%lowerBoundTheta_(iGRU)
    case(iLookPARAM%upperBoundTemp); allocate(parVector(1)); parVector = mparStruct%upperBoundTemp_(iGRU)
    case(iLookPARAM%lowerBoundTemp); allocate(parVector(1)); parVector = mparStruct%lowerBoundTemp_(iGRU)
    case(iLookPARAM%tempCritRain); allocate(parVector(1)); parVector = mparStruct%tempCritRain_(iGRU)
    case(iLookPARAM%tempRangeTimestep); allocate(parVector(1)); parVector = mparStruct%tempRangeTimestep_(iGRU)
    case(iLookPARAM%frozenPrecipMultip); allocate(parVector(1)); parVector = mparStruct%frozenPrecipMultip_(iGRU)
    case(iLookPARAM%snowfrz_scale); allocate(parVector(1)); parVector = mparStruct%snowfrz_scale_(iGRU)
    case(iLookPARAM%fixedThermalCond_snow); allocate(parVector(1)); parVector = mparStruct%fixedThermalCond_snow_(iGRU)
    case(iLookPARAM%albedoMax); allocate(parVector(1)); parVector = mparStruct%albedoMax_(iGRU)
    case(iLookPARAM%albedoMinWinter); allocate(parVector(1)); parVector = mparStruct%albedoMinWinter_(iGRU)
    case(iLookPARAM%albedoMinSpring); allocate(parVector(1)); parVector = mparStruct%albedoMinSpring_(iGRU)
    case(iLookPARAM%albedoMaxVisible); allocate(parVector(1)); parVector = mparStruct%albedoMaxVisible_(iGRU)
    case(iLookPARAM%albedoMinVisible); allocate(parVector(1)); parVector = mparStruct%albedoMinVisible_(iGRU)
    case(iLookPARAM%albedoMaxNearIR); allocate(parVector(1)); parVector = mparStruct%albedoMaxNearIR_(iGRU)
    case(iLookPARAM%albedoMinNearIR); allocate(parVector(1)); parVector = mparStruct%albedoMinNearIR_(iGRU)
    case(iLookPARAM%albedoDecayRate); allocate(parVector(1)); parVector = mparStruct%albedoDecayRate_(iGRU)
    case(iLookPARAM%albedoSootLoad); allocate(parVector(1)); parVector = mparStruct%albedoSootLoad_(iGRU)
    case(iLookPARAM%albedoRefresh); allocate(parVector(1)); parVector = mparStruct%albedoRefresh_(iGRU)
    case(iLookPARAM%directScale); allocate(parVector(1)); parVector = mparStruct%directScale_(iGRU)
    case(iLookPARAM%Frad_direct); allocate(parVector(1)); parVector = mparStruct%Frad_direct_(iGRU)
    case(iLookPARAM%Frad_vis); allocate(parVector(1)); parVector = mparStruct%Frad_vis_(iGRU)
    case(iLookPARAM%newSnowDenMin); allocate(parVector(1)); parVector = mparStruct%newSnowDenMin_(iGRU)
    case(iLookPARAM%newSnowDenMult); allocate(parVector(1)); parVector = mparStruct%newSnowDenMult_(iGRU)
    case(iLookPARAM%newSnowDenScal); allocate(parVector(1)); parVector = mparStruct%newSnowDenScal_(iGRU)
    case(iLookPARAM%constSnowDen); allocate(parVector(1)); parVector = mparStruct%constSnowDen_(iGRU)
    case(iLookPARAM%newSnowDenAdd); allocate(parVector(1)); parVector = mparStruct%newSnowDenAdd_(iGRU)
    case(iLookPARAM%newSnowDenMultTemp); allocate(parVector(1)); parVector = mparStruct%newSnowDenMultTemp_(iGRU)
    case(iLookPARAM%newSnowDenMultWind); allocate(parVector(1)); parVector = mparStruct%newSnowDenMultWind_(iGRU)
    case(iLookPARAM%newSnowDenMultAnd); allocate(parVector(1)); parVector = mparStruct%newSnowDenMultAnd_(iGRU)
    case(iLookPARAM%newSnowDenBase); allocate(parVector(1)); parVector = mparStruct%newSnowDenBase_(iGRU)
    case(iLookPARAM%densScalGrowth); allocate(parVector(1)); parVector = mparStruct%densScalGrowth_(iGRU)
    case(iLookPARAM%tempScalGrowth); allocate(parVector(1)); parVector = mparStruct%tempScalGrowth_(iGRU)
    case(iLookPARAM%grainGrowthRate); allocate(parVector(1)); parVector = mparStruct%grainGrowthRate_(iGRU)
    case(iLookPARAM%densScalOvrbdn); allocate(parVector(1)); parVector = mparStruct%densScalOvrbdn_(iGRU)
    case(iLookPARAM%tempScalOvrbdn); allocate(parVector(1)); parVector = mparStruct%tempScalOvrbdn_(iGRU)
    case(iLookPARAM%baseViscosity); allocate(parVector(1)); parVector = mparStruct%baseViscosity_(iGRU)
    case(iLookPARAM%Fcapil); allocate(parVector(1)); parVector = mparStruct%Fcapil_(iGRU)
    case(iLookPARAM%k_snow); allocate(parVector(1)); parVector = mparStruct%k_snow_(iGRU)
    case(iLookPARAM%mw_exp); allocate(parVector(1)); parVector = mparStruct%mw_exp_(iGRU)
    case(iLookPARAM%z0Snow); allocate(parVector(1)); parVector = mparStruct%z0Snow_(iGRU)
    case(iLookPARAM%z0Soil); allocate(parVector(1)); parVector = mparStruct%z0Soil_(iGRU)
    case(iLookPARAM%z0Canopy); allocate(parVector(1)); parVector = mparStruct%z0Canopy_(iGRU)
    case(iLookPARAM%zpdFraction); allocate(parVector(1)); parVector = mparStruct%zpdFraction_(iGRU)
    case(iLookPARAM%critRichNumber); allocate(parVector(1)); parVector = mparStruct%critRichNumber_(iGRU)
    case(iLookPARAM%Louis79_bparam); allocate(parVector(1)); parVector = mparStruct%Louis79_bparam_(iGRU)
    case(iLookPARAM%Louis79_cStar); allocate(parVector(1)); parVector = mparStruct%Louis79_cStar_(iGRU)
    case(iLookPARAM%Mahrt87_eScale); allocate(parVector(1)); parVector = mparStruct%Mahrt87_eScale_(iGRU)
    case(iLookPARAM%leafExchangeCoeff); allocate(parVector(1)); parVector = mparStruct%leafExchangeCoeff_(iGRU)
    case(iLookPARAM%windReductionParam); allocate(parVector(1)); parVector = mparStruct%windReductionParam_(iGRU)
    case(iLookPARAM%Kc25); allocate(parVector(1)); parVector = mparStruct%Kc25_(iGRU)
    case(iLookPARAM%Ko25); allocate(parVector(1)); parVector = mparStruct%Ko25_(iGRU)
    case(iLookPARAM%Kc_qFac); allocate(parVector(1)); parVector = mparStruct%Kc_qFac_(iGRU)
    case(iLookPARAM%Ko_qFac); allocate(parVector(1)); parVector = mparStruct%Ko_qFac_(iGRU)
    case(iLookPARAM%kc_Ha); allocate(parVector(1)); parVector = mparStruct%kc_Ha_(iGRU)
    case(iLookPARAM%ko_Ha); allocate(parVector(1)); parVector = mparStruct%ko_Ha_(iGRU)
    case(iLookPARAM%vcmax25_canopyTop); allocate(parVector(1)); parVector = mparStruct%vcmax25_canopyTop_(iGRU)
    case(iLookPARAM%vcmax_qFac); allocate(parVector(1)); parVector = mparStruct%vcmax_qFac_(iGRU)
    case(iLookPARAM%vcmax_Ha); allocate(parVector(1)); parVector = mparStruct%vcmax_Ha_(iGRU)
    case(iLookPARAM%vcmax_Hd); allocate(parVector(1)); parVector = mparStruct%vcmax_Hd_(iGRU)
    case(iLookPARAM%vcmax_Sv); allocate(parVector(1)); parVector = mparStruct%vcmax_Sv_(iGRU)
    case(iLookPARAM%vcmax_Kn); allocate(parVector(1)); parVector = mparStruct%vcmax_Kn_(iGRU)
    case(iLookPARAM%jmax25_scale); allocate(parVector(1)); parVector = mparStruct%jmax25_scale_(iGRU)
    case(iLookPARAM%jmax_Ha); allocate(parVector(1)); parVector = mparStruct%jmax_Ha_(iGRU)
    case(iLookPARAM%jmax_Hd); allocate(parVector(1)); parVector = mparStruct%jmax_Hd_(iGRU)
    case(iLookPARAM%jmax_Sv); allocate(parVector(1)); parVector = mparStruct%jmax_Sv_(iGRU)
    case(iLookPARAM%fractionJ); allocate(parVector(1)); parVector = mparStruct%fractionJ_(iGRU)
    case(iLookPARAM%quantamYield); allocate(parVector(1)); parVector = mparStruct%quantamYield_(iGRU)
    case(iLookPARAM%vpScaleFactor); allocate(parVector(1)); parVector = mparStruct%vpScaleFactor_(iGRU)
    case(iLookPARAM%cond2photo_slope); allocate(parVector(1)); parVector = mparStruct%cond2photo_slope_(iGRU)
    case(iLookPARAM%minStomatalConductance); allocate(parVector(1)); parVector = mparStruct%minStomatalConductance_(iGRU)
    case(iLookPARAM%winterSAI); allocate(parVector(1)); parVector = mparStruct%winterSAI_(iGRU)
    case(iLookPARAM%summerLAI); allocate(parVector(1)); parVector = mparStruct%summerLAI_(iGRU)
    case(iLookPARAM%rootScaleFactor1); allocate(parVector(1)); parVector = mparStruct%rootScaleFactor1_(iGRU)
    case(iLookPARAM%rootScaleFactor2); allocate(parVector(1)); parVector = mparStruct%rootScaleFactor2_(iGRU)
    case(iLookPARAM%rootingDepth); allocate(parVector(1)); parVector = mparStruct%rootingDepth_(iGRU)
    case(iLookPARAM%rootDistExp); allocate(parVector(1)); parVector = mparStruct%rootDistExp_(iGRU)
    case(iLookPARAM%plantWiltPsi); allocate(parVector(1)); parVector = mparStruct%plantWiltPsi_(iGRU)
    case(iLookPARAM%soilStressParam); allocate(parVector(1)); parVector = mparStruct%soilStressParam_(iGRU)
    case(iLookPARAM%critSoilWilting); allocate(parVector(1)); parVector = mparStruct%critSoilWilting_(iGRU)
    case(iLookPARAM%critSoilTranspire); allocate(parVector(1)); parVector = mparStruct%critSoilTranspire_(iGRU)
    case(iLookPARAM%critAquiferTranspire); allocate(parVector(1)); parVector = mparStruct%critAquiferTranspire_(iGRU)
    case(iLookPARAM%minStomatalResistance); allocate(parVector(1)); parVector = mparStruct%minStomatalResistance_(iGRU)
    case(iLookPARAM%leafDimension); allocate(parVector(1)); parVector = mparStruct%leafDimension_(iGRU)
    case(iLookPARAM%heightCanopyTop); allocate(parVector(1)); parVector = mparStruct%heightCanopyTop_(iGRU)
    case(iLookPARAM%heightCanopyBottom); allocate(parVector(1)); parVector = mparStruct%heightCanopyBottom_(iGRU)
    case(iLookPARAM%specificHeatVeg); allocate(parVector(1)); parVector = mparStruct%specificHeatVeg_(iGRU)
    case(iLookPARAM%maxMassVegetation); allocate(parVector(1)); parVector = mparStruct%maxMassVegetation_(iGRU)
    case(iLookPARAM%throughfallScaleRain); allocate(parVector(1)); parVector = mparStruct%throughfallScaleRain_(iGRU)
    case(iLookPARAM%refInterceptCapSnow); allocate(parVector(1)); parVector = mparStruct%refInterceptCapSnow_(iGRU)
    case(iLookPARAM%refInterceptCapRain); allocate(parVector(1)); parVector = mparStruct%refInterceptCapRain_(iGRU)
    case(iLookPARAM%snowUnloadingCoeff); allocate(parVector(1)); parVector = mparStruct%snowUnloadingCoeff_(iGRU)
    case(iLookPARAM%canopyDrainageCoeff); allocate(parVector(1)); parVector = mparStruct%canopyDrainageCoeff_(iGRU)
    case(iLookPARAM%ratioDrip2Unloading); allocate(parVector(1)); parVector = mparStruct%ratioDrip2Unloading_(iGRU)
    case(iLookPARAM%canopyWettingFactor); allocate(parVector(1)); parVector = mparStruct%canopyWettingFactor_(iGRU)
    case(iLookPARAM%canopyWettingExp); allocate(parVector(1)); parVector = mparStruct%canopyWettingExp_(iGRU)
    case(iLookPARAM%minTempUnloading); allocate(parVector(1)); parVector = mparStruct%minTempUnloading_(iGRU)
    case(iLookPARAM%rateTempUnloading); allocate(parVector(1)); parVector = mparStruct%rateTempUnloading_(iGRU)
    case(iLookPARAM%minWindUnloading); allocate(parVector(1)); parVector = mparStruct%minWindUnloading_(iGRU)
    case(iLookPARAM%rateWindUnloading); allocate(parVector(1)); parVector = mparStruct%rateWindUnloading_(iGRU)
    case(iLookPARAM%soil_dens_intr); parVector = mparStruct%soil_dens_intr_(:,iGRU)
    case(iLookPARAM%thCond_soil); parVector = mparStruct%thCond_soil_(:,iGRU)
    case(iLookPARAM%frac_sand); parVector = mparStruct%frac_sand_(:,iGRU)
    case(iLookPARAM%frac_silt); parVector = mparStruct%frac_silt_(:,iGRU)
    case(iLookPARAM%frac_clay); parVector = mparStruct%frac_clay_(:,iGRU)
    case(iLookPARAM%theta_sat); parVector = mparStruct%theta_sat_(:,iGRU)
    case(iLookPARAM%theta_res); parVector = mparStruct%theta_res_(:,iGRU)
    case(iLookPARAM%vGn_alpha); parVector = mparStruct%vGn_alpha_(:,iGRU)
    case(iLookPARAM%vGn_n); parVector = mparStruct%vGn_n_(:,iGRU)
    case(iLookPARAM%k_soil); parVector = mparStruct%k_soil_(:,iGRU)
    case(iLookPARAM%k_macropore); parVector = mparStruct%k_macropore_(:,iGRU)
    case(iLookPARAM%fieldCapacity); allocate(parVector(1)); parVector = mparStruct%fieldCapacity_(iGRU)
    case(iLookPARAM%wettingFrontSuction); allocate(parVector(1)); parVector = mparStruct%wettingFrontSuction_(iGRU)
    case(iLookPARAM%theta_mp); allocate(parVector(1)); parVector = mparStruct%theta_mp_(iGRU)
    case(iLookPARAM%mpExp); allocate(parVector(1)); parVector = mparStruct%mpExp_(iGRU)
    case(iLookPARAM%kAnisotropic); allocate(parVector(1)); parVector = mparStruct%kAnisotropic_(iGRU)
    case(iLookPARAM%zScale_TOPMODEL); allocate(parVector(1)); parVector = mparStruct%zScale_TOPMODEL_(iGRU)
    case(iLookPARAM%compactedDepth); allocate(parVector(1)); parVector = mparStruct%compactedDepth_(iGRU)
    case(iLookPARAM%aquiferBaseflowRate); allocate(parVector(1)); parVector = mparStruct%aquiferBaseflowRate_(iGRU)
    case(iLookPARAM%aquiferScaleFactor); allocate(parVector(1)); parVector = mparStruct%aquiferScaleFactor_(iGRU)
    case(iLookPARAM%aquiferBaseflowExp); allocate(parVector(1)); parVector = mparStruct%aquiferBaseflowExp_(iGRU)
    case(iLookPARAM%qSurfScale); allocate(parVector(1)); parVector = mparStruct%qSurfScale_(iGRU)
    case(iLookPARAM%specificStorage); allocate(parVector(1)); parVector = mparStruct%specificStorage_(iGRU)
    case(iLookPARAM%f_impede); allocate(parVector(1)); parVector = mparStruct%f_impede_(iGRU)
    case(iLookPARAM%soilIceScale); allocate(parVector(1)); parVector = mparStruct%soilIceScale_(iGRU)
    case(iLookPARAM%soilIceCV); allocate(parVector(1)); parVector = mparStruct%soilIceCV_(iGRU)
    case(iLookPARAM%minwind); allocate(parVector(1)); parVector = mparStruct%minwind_(iGRU)
    case(iLookPARAM%minstep); allocate(parVector(1)); parVector = mparStruct%minstep
    case(iLookPARAM%maxstep); allocate(parVector(1)); parVector = mparStruct%maxstep
    case(iLookPARAM%be_steps); allocate(parVector(1)); parVector = mparStruct%be_steps
    case(iLookPARAM%relTolTempCas); allocate(parVector(1)); parVector = mparStruct%relTolTempCas
    case(iLookPARAM%absTolTempCas); allocate(parVector(1)); parVector = mparStruct%absTolTempCas
    case(iLookPARAM%relTolTempVeg); allocate(parVector(1)); parVector = mparStruct%relTolTempVeg
    case(iLookPARAM%absTolTempVeg); allocate(parVector(1)); parVector = mparStruct%absTolTempVeg
    case(iLookPARAM%relTolWatVeg); allocate(parVector(1)); parVector = mparStruct%relTolWatVeg
    case(iLookPARAM%absTolWatVeg); allocate(parVector(1)); parVector = mparStruct%absTolWatVeg
    case(iLookPARAM%relTolTempSoilSnow); allocate(parVector(1)); parVector = mparStruct%relTolTempSoilSnow
    case(iLookPARAM%absTolTempSoilSnow); allocate(parVector(1)); parVector = mparStruct%absTolTempSoilSnow
    case(iLookPARAM%relTolWatSnow); allocate(parVector(1)); parVector = mparStruct%relTolWatSnow
    case(iLookPARAM%absTolWatSnow); allocate(parVector(1)); parVector = mparStruct%absTolWatSnow
    case(iLookPARAM%relTolMatric); allocate(parVector(1)); parVector = mparStruct%relTolMatric
    case(iLookPARAM%absTolMatric); allocate(parVector(1)); parVector = mparStruct%absTolMatric
    case(iLookPARAM%relTolAquifr); allocate(parVector(1)); parVector = mparStruct%relTolAquifr
    case(iLookPARAM%absTolAquifr); allocate(parVector(1)); parVector = mparStruct%absTolAquifr
    case(iLookPARAM%idaMaxOrder); allocate(parVector(1)); parVector = mparStruct%idaMaxOrder
    case(iLookPARAM%idaMaxInternalSteps); allocate(parVector(1)); parVector = mparStruct%idaMaxInternalSteps
    case(iLookPARAM%idaInitStepSize); allocate(parVector(1)); parVector = mparStruct%idaInitStepSize
    case(iLookPARAM%idaMinStepSize); allocate(parVector(1)); parVector = mparStruct%idaMinStepSize
    case(iLookPARAM%idaMaxErrTestFail); allocate(parVector(1)); parVector = mparStruct%idaMaxErrTestFail
    case(iLookPARAM%zmin); allocate(parVector(1)); parVector = mparStruct%zmin_(iGRU)
    case(iLookPARAM%zmax); allocate(parVector(1)); parVector = mparStruct%zmax_(iGRU)
    case(iLookPARAM%zminLayer1); allocate(parVector(1)); parVector = mparStruct%zminLayer1_(iGRU)
    case(iLookPARAM%zminLayer2); allocate(parVector(1)); parVector = mparStruct%zminLayer2_(iGRU)
    case(iLookPARAM%zminLayer3); allocate(parVector(1)); parVector = mparStruct%zminLayer3_(iGRU)
    case(iLookPARAM%zminLayer4); allocate(parVector(1)); parVector = mparStruct%zminLayer4_(iGRU)
    case(iLookPARAM%zminLayer5); allocate(parVector(1)); parVector = mparStruct%zminLayer5_(iGRU)
    case(iLookPARAM%zmaxLayer1_lower); allocate(parVector(1)); parVector = mparStruct%zmaxLayer1_lower_(iGRU)
    case(iLookPARAM%zmaxLayer2_lower); allocate(parVector(1)); parVector = mparStruct%zmaxLayer2_lower_(iGRU)
    case(iLookPARAM%zmaxLayer3_lower); allocate(parVector(1)); parVector = mparStruct%zmaxLayer3_lower_(iGRU)
    case(iLookPARAM%zmaxLayer4_lower); allocate(parVector(1)); parVector = mparStruct%zmaxLayer4_lower_(iGRU)
    case(iLookPARAM%zmaxLayer1_upper); allocate(parVector(1)); parVector = mparStruct%zmaxLayer1_upper_(iGRU)
    case(iLookPARAM%zmaxLayer2_upper); allocate(parVector(1)); parVector = mparStruct%zmaxLayer2_upper_(iGRU)
    case(iLookPARAM%zmaxLayer3_upper); allocate(parVector(1)); parVector = mparStruct%zmaxLayer3_upper_(iGRU)
    case(iLookPARAM%zmaxLayer4_upper); allocate(parVector(1)); parVector = mparStruct%zmaxLayer4_upper_(iGRU)
    case(iLookPARAM%FUSE_Ac_max); allocate(parVector(1)); parVector = mparStruct%FUSE_Ac_max_(iGRU)
    case(iLookPARAM%FUSE_phi_tens); allocate(parVector(1)); parVector = mparStruct%FUSE_phi_tens_(iGRU)
    case(iLookPARAM%FUSE_b); allocate(parVector(1)); parVector = mparStruct%FUSE_b_(iGRU)
    case(iLookPARAM%FUSE_lambda); allocate(parVector(1)); parVector = mparStruct%FUSE_lambda_(iGRU)
    case(iLookPARAM%FUSE_chi); allocate(parVector(1)); parVector = mparStruct%FUSE_chi_(iGRU)
    case(iLookPARAM%FUSE_mu); allocate(parVector(1)); parVector = mparStruct%FUSE_mu_(iGRU)
    case(iLookPARAM%FUSE_n); allocate(parVector(1)); parVector = mparStruct%FUSE_n_(iGRU)
    end select
    if (.not.allocated(parVector)) then
      allocate(parVector(1)); parVector = -9999_rkind
    end if
  end subroutine


  subroutine set_device_param_data(mparStruct, iGRU, ixParam, parVector)
    implicit none
    type(mpar_data_device) :: mparStruct
    integer(i4b) :: iGRU, ixParam
    real(rkind) :: parVector(:)

    select case(ixParam)
    case(iLookPARAM%soil_dens_intr); mparStruct%soil_dens_intr_(:,iGRU) = parVector
    case(iLookPARAM%thCond_soil); mparStruct%thCond_soil_(:,iGRU) = parVector
    case(iLookPARAM%frac_sand); mparStruct%frac_sand_(:,iGRU) = parVector
    case(iLookPARAM%frac_silt); mparStruct%frac_silt_(:,iGRU) = parVector
    case(iLookPARAM%frac_clay); mparStruct%frac_clay_(:,iGRU) = parVector
    case(iLookPARAM%theta_sat); mparStruct%theta_sat_(:,iGRU) = parVector
    case(iLookPARAM%theta_res); mparStruct%theta_res_(:,iGRU) = parVector
    case(iLookPARAM%vGn_alpha); mparStruct%vGn_alpha_(:,iGRU) = parVector
    case(iLookPARAM%vGn_n); mparStruct%vGn_n_(:,iGRU) = parVector
    case(iLookPARAM%k_soil); mparStruct%k_soil_(:,iGRU) = parVector
    case(iLookPARAM%k_macropore); mparStruct%k_macropore_(:,iGRU) = parVector
    end select
  end subroutine


  subroutine deallocate_device_param_data(mpar_data_device)
    type(mpar_data_device), intent(inout) :: mpar_data_device
    deallocate(mpar_data_device%Fcapil_)
    deallocate(mpar_data_device%Louis79_bparam_)
    deallocate(mpar_data_device%Louis79_cStar_)
    deallocate(mpar_data_device%Mahrt87_eScale_)
    deallocate(mpar_data_device%aquiferBaseflowExp_)
    deallocate(mpar_data_device%aquiferBaseflowRate_)
    deallocate(mpar_data_device%aquiferScaleFactor_)
    deallocate(mpar_data_device%canopyDrainageCoeff_)
    deallocate(mpar_data_device%canopyWettingExp_)
    deallocate(mpar_data_device%canopyWettingFactor_)
    deallocate(mpar_data_device%critAquiferTranspire_)
    deallocate(mpar_data_device%critRichNumber_)
    deallocate(mpar_data_device%critSoilTranspire_)
    deallocate(mpar_data_device%critSoilWilting_)
    deallocate(mpar_data_device%f_impede_)
    deallocate(mpar_data_device%fixedThermalCond_snow_)
    deallocate(mpar_data_device%frac_clay_)
    deallocate(mpar_data_device%frac_sand_)
    deallocate(mpar_data_device%frac_silt_)
    deallocate(mpar_data_device%heightCanopyBottom_)
    deallocate(mpar_data_device%heightCanopyTop_)
    deallocate(mpar_data_device%kAnisotropic_)
    deallocate(mpar_data_device%k_snow_)
    deallocate(mpar_data_device%leafDimension_)
    deallocate(mpar_data_device%leafExchangeCoeff_)
    deallocate(mpar_data_device%lowerBoundHead_)
    deallocate(mpar_data_device%lowerBoundTemp_)
    deallocate(mpar_data_device%lowerBoundTheta_)
    deallocate(mpar_data_device%maxMassVegetation_)
    deallocate(mpar_data_device%minStomatalResistance_)
    deallocate(mpar_data_device%mpExp_)
    deallocate(mpar_data_device%mw_exp_)
    deallocate(mpar_data_device%plantWiltPsi_)
    deallocate(mpar_data_device%qSurfScale_)
    deallocate(mpar_data_device%rootingDepth_)
    deallocate(mpar_data_device%snowfrz_scale_)
    deallocate(mpar_data_device%soilIceCV_)
    deallocate(mpar_data_device%soilIceScale_)
    deallocate(mpar_data_device%soilStressParam_)
    deallocate(mpar_data_device%soil_dens_intr_)
    deallocate(mpar_data_device%specificHeatVeg_)
    deallocate(mpar_data_device%specificStorage_)
    deallocate(mpar_data_device%thCond_soil_)
    deallocate(mpar_data_device%theta_mp_)
    deallocate(mpar_data_device%theta_res_)
    deallocate(mpar_data_device%theta_sat_)
    deallocate(mpar_data_device%throughfallScaleRain_)
    deallocate(mpar_data_device%upperBoundHead_)
    deallocate(mpar_data_device%upperBoundTemp_)
    deallocate(mpar_data_device%upperBoundTheta_)
    deallocate(mpar_data_device%vGn_alpha_)
    deallocate(mpar_data_device%vGn_n_)
    deallocate(mpar_data_device%wettingFrontSuction_)
    deallocate(mpar_data_device%windReductionParam_)
    deallocate(mpar_data_device%z0Canopy_)
    deallocate(mpar_data_device%z0Snow_)
    deallocate(mpar_data_device%z0Soil_)
    deallocate(mpar_data_device%zScale_TOPMODEL_)
    deallocate(mpar_data_device%zpdFraction_)
    deallocate(mpar_data_device%Kc25_)
    deallocate(mpar_data_device%Ko25_)
    deallocate(mpar_data_device%Kc_qFac_)
    deallocate(mpar_data_device%Ko_qFac_)
    deallocate(mpar_data_device%kc_Ha_)
    deallocate(mpar_data_device%ko_Ha_)
    deallocate(mpar_data_device%vcmax25_canopyTop_)
    deallocate(mpar_data_device%vcmax_qFac_)
    deallocate(mpar_data_device%vcmax_Ha_)
    deallocate(mpar_data_device%vcmax_Hd_)
    deallocate(mpar_data_device%vcmax_Sv_)
    deallocate(mpar_data_device%vcmax_Kn_)
    deallocate(mpar_data_device%jmax25_scale_)
    deallocate(mpar_data_device%jmax_Ha_)
    deallocate(mpar_data_device%jmax_Hd_)
    deallocate(mpar_data_device%jmax_Sv_)
    deallocate(mpar_data_device%fractionJ_)
    deallocate(mpar_data_device%quantamYield_)
    deallocate(mpar_data_device%vpScaleFactor_)
    deallocate(mpar_data_device%cond2photo_slope_)
    deallocate(mpar_data_device%minStomatalConductance_)
      
  end subroutine deallocate_device_param_data


  subroutine allocate_device_indx_data(indx_data_device, indx_data, nSoil,nGRU)
    use globalData,only:maxSnowLayers,integerMissing
    implicit none
    type(indx_data_device), intent(inout) :: indx_data_device
    type(gru_hru_intVec),intent(in) :: indx_data              ! indices defining model states and layers
    integer(i4b), intent(in) :: nSoil,nGRU
    integer(i4b) :: iGRU
      
    integer(i4b) :: nSnow, nLayers,nState,maxSize
    nSnow = maxSnowLayers+1
    nLayers = nSoil + maxSnowLayers+1
    nState = 4+2*nLayers

    allocate(indx_data_device%nSnow(nGRU))
    do iGRU=1,nGRU
    indx_data_device%nSnow(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
    end do
    indx_data_device%nSoil = indx_data%gru(1)%hru(1)%var(iLookINDEX%nSoil)%dat(1)
    allocate(indx_data_device%nLayers_d(nGRU))
    do iGRU=1,nGRU
    indx_data_device%nLayers_d(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nLayers)%dat(1)
    end do
  allocate(indx_data_device%nState(nGRU))
  do iGRU=1,nGRU
  indx_data_device%nState(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nState)%dat(1)
  end do
  
indx_data_device%numberDomainSplitMass = indx_data%gru(1)%hru(1)%var(iLookINDEX%numberDomainSplitMass)%dat(1)
indx_data_device%numberDomainSplitNrg = indx_data%gru(1)%hru(1)%var(iLookINDEX%numberDomainSplitNrg)%dat(1)
indx_data_device%numberScalarSolutions = indx_data%gru(1)%hru(1)%var(iLookINDEX%numberScalarSolutions)%dat(1)
indx_data_device%numberStateSplit = indx_data%gru(1)%hru(1)%var(iLookINDEX%numberStateSplit)%dat(1)
    allocate(indx_data_device%ixNrgLayer(nLayers,nGRU))
    allocate(indx_data_device%ixHydLayer(nLayers,nGRU))
    ! indx_data_device%maxSnow = nSnow

    ! indx_data_device%ixCasNrg = indx_data%var(iLookINDEX%ixCasNrg)%dat(1)
    ! indx_data_device%ixVegNrg = indx_data%var(iLookINDEX%ixVegNrg)%dat(1)
    ! indx_data_device%ixVegHyd = indx_data%var(iLookINDEX%ixVegHyd)%dat(1)
    ! indx_data_device%ixAqWat = indx_data%var(iLookINDEX%ixAqWat)%dat(1)
    allocate(indx_data_device%ixCasNrg(nGRU))
    do iGRU=1,nGRU
    indx_data_device%ixCasNrg(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixCasNrg)%dat(1)
    end do
    allocate(indx_data_device%ixVegNrg(nGRU))
    do iGRU=1,nGRU
    indx_data_device%ixVegNrg(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixVegNrg)%dat(1)
    end do
    allocate(indx_data_device%ixVegHyd(nGRU))
    do iGRU=1,nGRU
    indx_data_device%ixVegHyd(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixVegHyd)%dat(1)
    end do
    allocate(indx_data_device%ixAqWat(nGRU))
    do iGRU=1,nGRU
    indx_data_device%ixAqWat(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixAqWat)%dat(1)
    end do
    ! indx_data_device%ixTopNrg = indx_data%var(iLookINDEX%ixTopNrg)%dat(1)
    ! indx_data_device%ixTopHyd = indx_data%var(iLookINDEX%ixTopHyd)%dat(1)
    allocate(indx_data_device%ixTopNrg(nGRU))
    do iGRU=1,nGRU
    indx_data_device%ixTopNrg(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixTopNrg)%dat(1)
    end do
    allocate(indx_data_device%ixTopHyd(nGRU))
    do iGRU=1,nGRU
    indx_data_device%ixTopHyd(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixTopHyd)%dat(1)
    end do
    allocate(indx_data_device%ixSnowSoilNrg(nLayers,nGRU))
    allocate(indx_data_device%ixSnowSoilHyd(nLayers,nGRU))
    allocate(indx_data_device%ixHydType(nLayers,nGRU))
    allocate(indx_data_device%ixStateType(nState,nGRU))
    allocate(indx_data_device%ixHydCanopy(nGRU))
    allocate(indx_data_device%ixNrgCanopy(nGRU))
    allocate(indx_data_device%ixNrgCanair(nGRU))
    do iGRU=1,nGRU
    indx_data_device%ixNrgCanair(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgCanair)%dat(1)
    end do
    allocate(indx_data_device%ixControlVolume(nState,nGRU))
    allocate(indx_data_device%layerType(nLayers,nGRU))
    allocate(indx_data_device%ixSoilOnlyHyd(nSoil,nGRU))
    allocate(indx_data_device%ixSnowOnlyHyd(nSnow,nGRU))
    allocate(indx_data_device%ixSoilOnlyNrg(nSoil,nGRU))
    allocate(indx_data_device%ixSnowOnlyNrg(nSnow,nGRU))
    allocate(indx_data_device%nCasNrg(nGRU))
    do iGRU=1,nGRU
    indx_data_device%nCasNrg(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nCasNrg)%dat(1)
    end do
    allocate(indx_data_device%nVegNrg(nGRU))
    do iGRU=1,nGRU
    indx_data_device%nVegNrg(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nVegNrg)%dat(1)
    end do
    allocate(indx_data_device%nVegMass(nGRU))
    do iGRU=1,nGRU
    indx_data_device%nVegMass(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nVegMass)%dat(1)
    end do
    allocate(indx_data_device%nVegState(nGRU))
    do iGRU=1,nGRU
    indx_data_device%nVegState(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nVegState)%dat(1)
    end do
    allocate(indx_data_device%nNrgState(nGRU))
    do iGRU=1,nGRU
    indx_data_device%nNrgState(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nNrgState)%dat(1)
    end do
    allocate(indx_data_device%nWatState(nGRU))
    do iGRU=1,nGRU
    indx_data_device%nWatState(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nWatState)%dat(1)
    end do
    allocate(indx_data_device%nMatState(nGRU))
    do iGRU=1,nGRU
    indx_data_device%nMatState(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nMatState)%dat(1)
    end do
    allocate(indx_data_device%nMassState(nGRU))
    do iGRU=1,nGRU
    indx_data_device%nMassState(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nMassState)%dat(1)
    end do
    allocate(indx_data_device%ixWatAquifer(nGRU))
    do iGRU=1,nGRU
    indx_data_device%ixWatAquifer(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixWatAquifer)%dat(1)
    end do
        
    allocate(indx_data_device%ixSoilState(nSoil,nGRU))

    indx_data_device%numberFluxCalc = indx_data%gru(1)%hru(1)%var(iLookINDEX%numberFluxCalc)%dat(1)
    maxSize = 0
    do iGRU=1,nGRU
      maxSize = max(size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgOnly)%dat), maxSize)
    end do
    if (maxSize.ne.0) allocate(indx_data_device%ixNrgOnly(maxSize,nGRU))
    allocate(indx_data_device%ixDomainType(nState,nGRU))
    do iGRU=1,nGRU
    nLayers = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nLayers)%dat(1)
    nSnow = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
    nState = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nState)%dat(1)
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixStateType)%dat).ne.0) indx_data_device%ixStateType(1:nState,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixStateType)%dat
      if (allocated(indx_data_device%ixHydCanopy)) indx_data_device%ixHydCanopy(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydCanopy)%dat(1)
      if (allocated(indx_data_device%ixNrgCanopy)) indx_data_device%ixNrgCanopy(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgCanopy)%dat(1)
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixControlVolume)%dat).ne.0) indx_data_device%ixControlVolume(1:nState,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixControlVolume)%dat
      indx_data_device%layerType(1:nLayers,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%layerType)%dat(1:nLayers)

      indx_data_device%ixHydType(1:nLayers,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydType)%dat(1:nLayers)
      if (allocated(indx_data_device%ixNrgLayer)) indx_data_device%ixNrgLayer(1:nLayers,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgLayer)%dat(1:nLayers)
      if (allocated(indx_data_device%ixHydLayer)) indx_data_device%ixHydLayer(1:nLayers,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydLayer)%dat(1:nLayers)
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilHyd)%dat).ne.0) then
        indx_data_device%ixSnowSoilHyd(1:nLayers,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilHyd)%dat(1:nLayers)
      else
        indx_data_device%ixSnowSoilHyd(:,iGRU) = integerMissing
      endif
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilNrg)%dat).ne.0) then
        indx_data_device%ixSnowSoilNrg(1:nLayers,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilNrg)%dat(1:nLayers)
      else
        indx_data_device%ixSnowSoilNrg(:,iGRU) = integerMissing
      endif
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilOnlyHyd)%dat).ne.0) then
        indx_data_device%ixSoilOnlyHyd(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilOnlyHyd)%dat
      else
        indx_data_device%ixSoilOnlyHyd(:,iGRU) = integerMissing
      endif
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyHyd)%dat).ne.0) then
        indx_data_device%ixSnowOnlyHyd(1:nSnow,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyHyd)%dat(1:nSnow)
      else
        indx_data_device%ixSnowOnlyHyd(:,iGRU) = integerMissing
      endif
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilOnlyNrg)%dat).ne.0) then
        indx_data_device%ixSoilOnlyNrg(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilOnlyNrg)%dat
      else
        indx_data_device%ixSoilOnlyNrg(:,iGRU) = integerMissing
      endif
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyNrg)%dat).ne.0) then
        indx_data_device%ixSnowOnlyNrg(1:nSnow,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyNrg)%dat(1:nSnow)
      else
        indx_data_device%ixSnowOnlyNrg(:,iGRU) = integerMissing
      endif
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixDomainType)%dat).ne.0) indx_data_device%ixDomainType(1:nState,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixDomainType)%dat
      maxSize = size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgOnly)%dat)
      if (maxSize.ne.0) indx_data_device%ixNrgOnly(1:maxSize,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgOnly)%dat
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilState)%dat).ne.0) indx_data_device%ixSoilState(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilState)%dat

    end do
  
  end subroutine allocate_device_indx_data
  subroutine allocate_device_indx_temp(indx_data_device, indx_data,nSoil,nGRU)
    use globalData,only:maxSnowLayers
    implicit none
    type(indx_data_device), intent(inout) :: indx_data_device
    type(indx_data_device) :: indx_data
    integer(i4b), intent(in) :: nSoil,nGRU
    integer(i4b) :: iGRU
      
    integer(i4b) :: nSnow, nLayers,nState
    nSnow = maxSnowLayers+1
    nLayers = nSoil + maxSnowLayers+1
    nState = 4+2*nLayers

    allocate(indx_data_device%nSnow(nGRU))
    allocate(indx_data_device%nLayers_d(nGRU))
  allocate(indx_data_device%nState(nGRU))
  
    allocate(indx_data_device%ixNrgLayer(nLayers,nGRU))
    allocate(indx_data_device%ixHydLayer(nLayers,nGRU))
    ! indx_data_device%maxSnow = nSnow

    ! indx_data_device%ixCasNrg = indx_data%var(iLookINDEX%ixCasNrg)%dat(1)
    ! indx_data_device%ixVegNrg = indx_data%var(iLookINDEX%ixVegNrg)%dat(1)
    ! indx_data_device%ixVegHyd = indx_data%var(iLookINDEX%ixVegHyd)%dat(1)
    ! indx_data_device%ixAqWat = indx_data%var(iLookINDEX%ixAqWat)%dat(1)
    allocate(indx_data_device%ixCasNrg(nGRU))
    allocate(indx_data_device%ixVegNrg(nGRU))
    allocate(indx_data_device%ixVegHyd(nGRU))
    allocate(indx_data_device%ixAqWat(nGRU))
    allocate(indx_data_device%ixTopNrg(nGRU))
    allocate(indx_data_device%ixTopHyd(nGRU))
    allocate(indx_data_device%ixSnowSoilNrg(nLayers,nGRU))
    allocate(indx_data_device%ixSnowSoilHyd(nLayers,nGRU))
    allocate(indx_data_device%ixHydType(nLayers,nGRU))
    allocate(indx_data_device%ixStateType(nState,nGRU))
    allocate(indx_data_device%ixHydCanopy(nGRU))
    allocate(indx_data_device%ixNrgCanopy(nGRU))
    allocate(indx_data_device%ixNrgCanair(nGRU))
    allocate(indx_data_device%ixControlVolume(nState,nGRU))
    allocate(indx_data_device%layerType(nLayers,nGRU))
    allocate(indx_data_device%ixSoilOnlyHyd(nSoil,nGRU))
    allocate(indx_data_device%ixSnowOnlyHyd(nSnow,nGRU))
    allocate(indx_data_device%ixSoilOnlyNrg(nSoil,nGRU))
    allocate(indx_data_device%ixSnowOnlyNrg(nSnow,nGRU))
        
    allocate(indx_data_device%ixSoilState(nSoil,nGRU))
    allocate(indx_data_device%nCasNrg(nGRU))
    allocate(indx_data_device%nVegNrg(nGRU))
    allocate(indx_data_device%nVegMass(nGRU))
    allocate(indx_data_device%nVegState(nGRU))
    allocate(indx_data_device%nNrgState(nGRU))
    allocate(indx_data_device%nWatState(nGRU))
    allocate(indx_data_device%nMatState(nGRU))
    allocate(indx_data_device%nMassState(nGRU))
    allocate(indx_data_device%ixWatAquifer(nGRU))

    indx_data_device%numberFluxCalc = indx_data%numberFluxCalc
    if (size(indx_data%ixNrgOnly).ne.0) allocate(indx_data_device%ixNrgOnly(size(indx_data%ixNrgOnly,1),nGRU))
    allocate(indx_data_device%ixDomainType(nState,nGRU))
  
  end subroutine allocate_device_indx_temp
  subroutine initialize_device_indx_data(indx_data_device, indx_data, nSoil,nGRU)
    use globalData,only:maxSnowLayers,integerMissing
    implicit none
    type(indx_data_device), intent(inout) :: indx_data_device
    type(var_ilength),intent(in) :: indx_data              ! indices defining model states and layers
    integer(i4b), intent(in) :: nSoil,nGRU
    integer(i4b) :: iGRU
      
    integer(i4b) :: nSnow, nLayers,nState
    nSnow = maxSnowLayers+1
    nLayers = nSoil + maxSnowLayers+1

    indx_data_device%nCasNrg = indx_data%var(iLookINDEX%nCasNrg)%dat(1)
indx_data_device%nVegNrg = indx_data%var(iLookINDEX%nVegNrg)%dat(1)
indx_data_device%nVegMass = indx_data%var(iLookINDEX%nVegMass)%dat(1)
indx_data_device%nVegState = indx_data%var(iLookINDEX%nVegState)%dat(1)
indx_data_device%nNrgState = indx_data%var(iLookINDEX%nNrgState)%dat(1)
indx_data_device%nWatState = indx_data%var(iLookINDEX%nWatState)%dat(1)
indx_data_device%nMatState = indx_data%var(iLookINDEX%nMatState)%dat(1)
indx_data_device%nMassState = indx_data%var(iLookINDEX%nMassState)%dat(1)
indx_data_device%ixWatAquifer = indx_data%var(iLookINDEX%ixWatAquifer)%dat(1)
    indx_data_device%nSnow = indx_data%var(iLookINDEX%nSnow)%dat(1)
    indx_data_device%nSoil = indx_data%var(iLookINDEX%nSoil)%dat(1)
    indx_data_device%nLayers_d = indx_data%var(iLookINDEX%nLayers)%dat(1)
  indx_data_device%nState = indx_data%var(iLookINDEX%nState)%dat(1)
  
indx_data_device%numberDomainSplitMass = indx_data%var(iLookINDEX%numberDomainSplitMass)%dat(1)
indx_data_device%numberDomainSplitNrg = indx_data%var(iLookINDEX%numberDomainSplitNrg)%dat(1)
indx_data_device%numberScalarSolutions = indx_data%var(iLookINDEX%numberScalarSolutions)%dat(1)
indx_data_device%numberStateSplit = indx_data%var(iLookINDEX%numberStateSplit)%dat(1)
    ! indx_data_device%maxSnow = nSnow

    ! indx_data_device%ixCasNrg = indx_data%var(iLookINDEX%ixCasNrg)%dat(1)
    ! indx_data_device%ixVegNrg = indx_data%var(iLookINDEX%ixVegNrg)%dat(1)
    ! indx_data_device%ixVegHyd = indx_data%var(iLookINDEX%ixVegHyd)%dat(1)
    ! indx_data_device%ixAqWat = indx_data%var(iLookINDEX%ixAqWat)%dat(1)
    indx_data_device%ixCasNrg = indx_data%var(iLookINDEX%ixCasNrg)%dat(1)
    indx_data_device%ixVegNrg = indx_data%var(iLookINDEX%ixVegNrg)%dat(1)
    indx_data_device%ixVegHyd = indx_data%var(iLookINDEX%ixVegHyd)%dat(1)
    indx_data_device%ixAqWat = indx_data%var(iLookINDEX%ixAqWat)%dat(1)
    ! indx_data_device%ixTopNrg = indx_data%var(iLookINDEX%ixTopNrg)%dat(1)
    ! indx_data_device%ixTopHyd = indx_data%var(iLookINDEX%ixTopHyd)%dat(1)
    indx_data_device%ixTopNrg = indx_data%var(iLookINDEX%ixTopNrg)%dat(1)
    indx_data_device%ixTopHyd = indx_data%var(iLookINDEX%ixTopHyd)%dat(1)
    indx_data_device%ixNrgCanair = indx_data%var(iLookINDEX%ixNrgCanair)%dat(1)
        

    indx_data_device%numberFluxCalc = indx_data%var(iLookINDEX%numberFluxCalc)%dat(1)
    nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)
    nSnow = indx_data%var(iLookINDEX%nSnow)%dat(1)
    nState = indx_data%var(iLookINDEX%nState)%dat(1)
    do iGRU=1,nGRU
      if (size(indx_data%var(iLookINDEX%ixStateType)%dat).ne.0) indx_data_device%ixStateType(1:nState,iGRU) = indx_data%var(iLookINDEX%ixStateType)%dat
      if (allocated(indx_data_device%ixHydCanopy)) indx_data_device%ixHydCanopy(iGRU) = indx_data%var(iLookINDEX%ixHydCanopy)%dat(1)
      if (allocated(indx_data_device%ixNrgCanopy)) indx_data_device%ixNrgCanopy(iGRU) = indx_data%var(iLookINDEX%ixNrgCanopy)%dat(1)
      if (allocated(indx_data_device%ixControlVolume)) indx_data_device%ixControlVolume(1:nState,iGRU) = indx_data%var(iLookINDEX%ixControlVolume)%dat
      indx_data_device%layerType(1:nLayers,iGRU) = indx_data%var(iLookINDEX%layerType)%dat(1:nLayers)

      indx_data_device%ixHydType(1:nLayers,iGRU) = indx_data%var(iLookINDEX%ixHydType)%dat(1:nLayers)
      if (allocated(indx_data_device%ixNrgLayer)) indx_data_device%ixNrgLayer(1:nLayers,iGRU) = indx_data%var(iLookINDEX%ixNrgLayer)%dat(1:nLayers)
      if (allocated(indx_data_device%ixHydLayer)) indx_data_device%ixHydLayer(1:nLayers,iGRU) = indx_data%var(iLookINDEX%ixHydLayer)%dat(1:nLayers)
      if (size(indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat).ne.0) then
        indx_data_device%ixSnowSoilHyd(1:nLayers,iGRU) = indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat(1:nLayers)
      else
        indx_data_device%ixSnowSoilHyd(:,iGRU) = integerMissing
      endif
      if (size(indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat).ne.0) then
        indx_data_device%ixSnowSoilNrg(1:nLayers,iGRU) = indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat(1:nLayers)
      else
        indx_data_device%ixSnowSoilNrg(:,iGRU) = integerMissing
      endif
      if (size(indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat).ne.0) then
        indx_data_device%ixSoilOnlyHyd(:,iGRU) = indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat
      else
        indx_data_device%ixSoilOnlyHyd(:,iGRU) = integerMissing
      endif
      if (size(indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat).ne.0) then
        indx_data_device%ixSnowOnlyHyd(1:nSnow,iGRU) = indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat(1:nSnow)
      else
        indx_data_device%ixSnowOnlyHyd(:,iGRU) = integerMissing
      endif
      if (size(indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat).ne.0) then
        indx_data_device%ixSoilOnlyNrg(:,iGRU) = indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat
      else
        indx_data_device%ixSoilOnlyNrg(:,iGRU) = integerMissing
      endif
      if (size(indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat).ne.0) then
        indx_data_device%ixSnowOnlyNrg(1:nSnow,iGRU) = indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat(1:nSnow)
      else
        indx_data_device%ixSnowOnlyNrg(:,iGRU) = integerMissing
      endif
      if (size(indx_data%var(iLookINDEX%ixDomainType)%dat).ne.0) indx_data_device%ixDomainType(1:nState,iGRU) = indx_data%var(iLookINDEX%ixDomainType)%dat
      if (allocated(indx_data_device%ixNrgOnly)) indx_data_device%ixNrgOnly(:,iGRU) = indx_data%var(iLookINDEX%ixNrgOnly)%dat
      if (size(indx_data%var(iLookINDEX%ixSoilState)%dat).ne.0) indx_data_device%ixSoilState(:,iGRU) = indx_data%var(iLookINDEX%ixSoilState)%dat

    end do
  
  end subroutine 

  subroutine finalize_device_indx_data(indx_data_device, indx_data)
    implicit none
    type(indx_data_device), intent(in) :: indx_data_device
    type(gru_hru_intVec),intent(inout) :: indx_data              ! indices defining model states and layers

    integer(i4b) :: nLayers,nSnow,nState,iGRU,nGRU
    nGRU = size(indx_data%gru)

    do iGRU=1,nGRU
        nLayers = indx_data_device%nLayers_d(iGRU)
    nSnow = indx_data_device%nSnow(iGRU)
    nState = indx_data_device%nState(iGRU)

    
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nCasNrg)%dat(1) = indx_data_device%nCasNrg(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nVegNrg)%dat(1) = indx_data_device%nVegNrg(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nVegMass)%dat(1) = indx_data_device%nVegMass(iGRU)
indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nVegState)%dat(1) = indx_data_device%nVegState(iGRU)
indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nNrgState)%dat(1) = indx_data_device%nNrgState(iGRU)
indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nWatState)%dat(1) = indx_data_device%nWatState(iGRU)
indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nMatState)%dat(1) = indx_data_device%nMatState(iGRU)
indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nMassState)%dat(1) = indx_data_device%nMassState(iGRU)
indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nState)%dat(1) = indx_data_device%nState(iGRU)
indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixWatAquifer)%dat(1) = indx_data_device%ixWatAquifer(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%numberFluxCalc)%dat(1) =indx_data_device%numberFluxCalc
    if (allocated(indx_data_device%ixNrgLayer)) then
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgLayer)%dat).ne.nLayers) then
        deallocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgLayer)%dat)
        allocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgLayer)%dat(1:nLayers))
      end if
      indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgLayer)%dat = indx_data_device%ixNrgLayer(1:nLayers,iGRU)
    end if
    if (allocated(indx_data_device%ixHydLayer)) then
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydLayer)%dat).ne.nLayers) then
        deallocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydLayer)%dat)
        allocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydLayer)%dat(1:nLayers))
      end if
      indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydLayer)%dat = indx_data_device%ixHydLayer(1:nLayers,iGRU)
    end if
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSoil)%dat(1) = indx_data_device%nSoil
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nLayers)%dat(1) = indx_data_device%nLayers_d(iGRU)
    if (allocated(indx_data_device%nSnow)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnow)%dat(1) = indx_data_device%nSnow(iGRU)
    if (allocated(indx_data_device%ixCasNrg)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixCasNrg)%dat(1) = indx_data_device%ixCasNrg(iGRU)
    if (allocated(indx_data_device%ixVegNrg)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixVegNrg)%dat(1) = indx_data_device%ixVegNrg(iGRU)
    if (allocated(indx_data_device%ixVegHyd)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixVegHyd)%dat(1) = indx_data_device%ixVegHyd(iGRU)
    if (allocated(indx_data_device%ixAqWat)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixAqWat)%dat(1) = indx_data_device%ixAqWat(iGRU)
    if (allocated(indx_data_device%ixTopHyd)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixTopHyd)%dat(1) = indx_data_device%ixTopHyd(iGRU)
    if (allocated(indx_data_device%ixTopNrg)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixTopNrg)%dat(1) = indx_data_device%ixTopNrg(iGRU)
    if (allocated(indx_data_device%ixSnowSoilNrg)) then
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilNrg)%dat).ne.nLayers) deallocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilNrg)%dat)
      indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilNrg)%dat = indx_data_device%ixSnowSoilNrg(1:nLayers,iGRU)
    end if
    
    if (allocated(indx_data_device%ixSnowSoilHyd )) then
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilHyd)%dat).ne.nLayers) deallocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilHyd)%dat)
      indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilHyd)%dat = indx_data_device%ixSnowSoilHyd(1:nLayers,iGRU)
    end if
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnowSoilNrg)%dat(1) = indx_data_device%nLayers_d(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnowSoilHyd)%dat(1) = indx_data_device%nLayers_d(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnowOnlyNrg)%dat(1) = indx_data_device%nSnow(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnowOnlyHyd)%dat(1) = indx_data_device%nSnow(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSoilOnlyNrg)%dat(1) = indx_data_device%nSoil
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSoilOnlyHyd)%dat(1) = indx_data_device%nSoil
    if (allocated(indx_data_device%ixHydType)) then
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydType)%dat).ne.nLayers) deallocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydType)%dat)
      indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydType)%dat = indx_data_device%ixHydType(1:nLayers,iGRU)
    end if
    if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixStateType)%dat).ne.nState) deallocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixStateType)%dat)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixStateType)%dat = indx_data_device%ixStateType(1:nState,iGRU)
    if (allocated(indx_data_device%ixHydCanopy)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydCanopy)%dat(1) = indx_data_device%ixHydCanopy(iGRU)
    if (allocated(indx_data_device%ixNrgCanopy)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgCanopy)%dat(1) = indx_data_device%ixNrgCanopy(iGRU)
    if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixControlVolume)%dat).ne.nState) deallocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixControlVolume)%dat)
    if (allocated(indx_data_device%ixControlVolume)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixControlVolume)%dat = indx_data_device%ixControlVolume(1:nState,iGRU)
    if (allocated(indx_data_device%layerType)) then
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%layerType)%dat).ne.nLayers) deallocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%layerType)%dat)
      indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%layerType)%dat = indx_data_device%layerType(1:nLayers,iGRU)
    end if
    if (allocated(indx_data_device%ixSoilOnlyHyd)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilOnlyHyd)%dat = indx_data_device%ixSoilOnlyHyd(:,iGRU)
    if (nSnow .ne. 0) then
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyHyd)%dat).ne.nSnow) deallocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyHyd)%dat)
      indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyHyd)%dat = indx_data_device%ixSnowOnlyHyd(1:nSnow,iGRU)
    end if
    if (allocated(indx_data_device%ixSoilOnlyNrg)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilOnlyNrg)%dat = indx_data_device%ixSoilOnlyNrg(:,iGRU)
    if (nSnow .ne. 0) then
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyNrg)%dat).ne.nSnow) deallocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyNrg)%dat)
      indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyNrg)%dat = indx_data_device%ixSnowOnlyNrg(1:nSnow,iGRU)
    end if
      
    if (allocated(indx_data_device%ixSoilState)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilState)%dat = indx_data_device%ixSoilState(:,iGRU)
    if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixDomainType)%dat).ne.nState) deallocate(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixDomainType)%dat)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixDomainType)%dat = indx_data_device%ixDomainType(1:nState,iGRU)
    if (allocated(indx_data_device%ixNrgOnly)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgOnly)%dat = indx_data_device%ixNrgOnly(:,iGRU)
      indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%numberDomainSplitMass)%dat(1) = indx_data_device%numberDomainSplitMass
indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%numberDomainSplitNrg)%dat(1) = indx_data_device%numberDomainSplitNrg
indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%numberScalarSolutions)%dat(1) = indx_data_device%numberScalarSolutions
indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%numberStateSplit)%dat(1) = indx_data_device%numberStateSplit

  end do
  end subroutine

      
  subroutine deallocate_device_indx_data(indx_data_device)
    type(indx_data_device), intent(inout) :: indx_data_device
    if (allocated(indx_data_device%ixNrgLayer)) deallocate(indx_data_device%ixNrgLayer)
    if (allocated(indx_data_device%ixHydLayer)) deallocate(indx_data_device%ixHydLayer)
    if (allocated(indx_data_device%nSnow)) deallocate(indx_data_device%nSnow)
    if (allocated(indx_data_device%ixCasNrg)) deallocate(indx_data_device%ixCasNrg)
    if (allocated(indx_data_device%ixVegNrg)) deallocate(indx_data_device%ixVegNrg)
    if (allocated(indx_data_device%ixVegHyd)) deallocate(indx_data_device%ixVegHyd)
    if (allocated(indx_data_device%ixAqWat)) deallocate(indx_data_device%ixAqWat)
    if (allocated(indx_data_device%ixTopHyd)) deallocate(indx_data_device%ixTopHyd)
    if (allocated(indx_data_device%ixTopNrg)) deallocate(indx_data_device%ixTopNrg)
    if (allocated(indx_data_device%ixHydType)) deallocate(indx_data_device%ixHydType)
    if (allocated(indx_data_device%ixStateType)) deallocate(indx_data_device%ixStateType)
    if (allocated(indx_data_device%ixHydCanopy)) deallocate(indx_data_device%ixHydCanopy)
    if (allocated(indx_data_device%ixNrgCanopy)) deallocate(indx_data_device%ixNrgCanopy)
    if (allocated(indx_data_device%ixControlVolume)) deallocate(indx_data_device%ixControlVolume)
    if (allocated(indx_data_device%layerType)) deallocate(indx_data_device%layerType)
    if (allocated(indx_data_device%ixSoilOnlyHyd)) deallocate(indx_data_device%ixSoilOnlyHyd)
    if (allocated(indx_data_device%ixSnowOnlyHyd)) deallocate(indx_data_device%ixSnowOnlyHyd)
    if (allocated(indx_data_device%ixSoilOnlyNrg)) deallocate(indx_data_device%ixSoilOnlyNrg)
    if (allocated(indx_data_device%ixSnowOnlyNrg)) deallocate(indx_data_device%ixSnowOnlyNrg)
    if (allocated(indx_data_device%ixSnowSoilNrg)) deallocate(indx_data_device%ixSnowSoilNrg)
    if (allocated(indx_data_device%ixSnowSoilHyd)) deallocate(indx_data_device%ixSnowSoilHyd)
    if (allocated(indx_data_device%ixNrgOnly)) deallocate(indx_data_device%ixNrgOnly)
    if (allocated(indx_data_device%ixDomainType)) deallocate(indx_data_device%ixDomainType)
    if (allocated(indx_data_device%ixSoilState)) deallocate(indx_data_device%ixSoilState)
      
  end subroutine deallocate_device_indx_data
      

  subroutine allocate_device_prog_data(prog_data_device, prog_data, nGRU, nSoil)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    type(gru_hru_doubleVec),intent(in) :: prog_data
    integer(i4b),intent(in) :: nGRU,nSoil

    integer(i4b) :: nLayers
    integer(i4b) :: iGRU
    nLayers = nSoil + maxSnowLayers+1
    ! state variables for vegetation
    allocate(prog_data_device%scalarCanopyIce(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyLiq(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyWat(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanairTemp(nGRU))                ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyTemp(nGRU))                ! get_ixVarType('scalarv')
    ! state variables for snow
    allocate(prog_data_device%spectralSnowAlbedoDiffuse(nSpecBands,nGRU))     ! get_ixVarType('wLength')
    allocate(prog_data_device%scalarSnowAlbedo(nGRU))                ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSnowDepth(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSWE(nGRU))                       ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSfcMeltPond(nGRU))               ! get_ixVarType('scalarv')
    ! define state variables for the snow+soil domain
    allocate(prog_data_device%mLayerTemp(nLayers,nGRU))                    ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerVolFracIce(nLayers,nGRU))              ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerVolFracLiq(nLayers,nGRU))              ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerVolFracWat(nLayers,nGRU))              ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerMatricHead(nSoil,nGRU))              ! get_ixVarType('midSoil')
    ! enthalpy
    allocate(prog_data_device%scalarCanairEnthalpy(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyEnthalpy(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%mLayerEnthalpy(nLayers,nGRU))                ! get_ixVarType('midToto')
    ! other state variables
    allocate(prog_data_device%scalarAquiferStorage(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSurfaceTemp(nGRU))               ! get_ixVarType('scalarv')
    ! define coordinate variables
    allocate(prog_data_device%mLayerDepth(nLayers,nGRU))                   ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerHeight(0:nLayers,nGRU))                  ! get_ixVarType('ifcToto')
    allocate(prog_data_device%iLayerHeight(0:nLayers,nGRU))                  ! get_ixVarType('ifcToto')
    do iGRU=1,nGRU

      print*, iGRU
    nLayers = size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerTemp)%dat)
      prog_data_device%spectralSnowAlbedoDiffuse(1:nSpecBands,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nSpecBands)
      prog_data_device%mLayerTemp(1:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerTemp)%dat(1:nLayers)
      prog_data_device%mLayerVolFracIce(1:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracIce)%dat(1:nLayers)
      prog_data_device%mLayerVolFracLiq(1:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracLiq)%dat(1:nLayers)
      prog_data_device%mLayerVolFracWat(1:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracWat)%dat(1:nLayers)
      prog_data_device%mLayerMatricHead(1:nSoil,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerMatricHead)%dat(1:nSoil)
      prog_data_device%mLayerEnthalpy(1:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerEnthalpy)%dat(1:nLayers)
      prog_data_device%mLayerDepth(1:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerDepth)%dat(1:nLayers)
      prog_data_device%mLayerHeight(0:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerHeight)%dat(0:nLayers)
      prog_data_device%iLayerHeight(0:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%iLayerHeight)%dat(0:nLayers)
    prog_data_device%scalarCanopyIce(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyIce)%dat(1)
    prog_data_device%scalarCanopyLiq(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyLiq)%dat(1)
    prog_data_device%scalarCanopyWat(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyWat)%dat(1)
    prog_data_device%scalarCanairTemp(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanairTemp)%dat(1)
    prog_data_device%scalarCanopyTemp(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyTemp)%dat(1)
    prog_data_device%scalarSnowAlbedo(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSnowAlbedo)%dat(1)
    prog_data_device%scalarSnowDepth(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSnowDepth)%dat(1)
    prog_data_device%scalarSWE(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSWE)%dat(1)
    prog_data_device%scalarSfcMeltPond(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSfcMeltPond)%dat(1)
    prog_data_device%scalarCanairEnthalpy(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanairEnthalpy)%dat(1)
    prog_data_device%scalarCanopyEnthalpy(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyEnthalpy)%dat(1)
    prog_data_device%scalarAquiferStorage(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarAquiferStorage)%dat(1)
    prog_data_device%scalarSurfaceTemp(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSurfaceTemp)%dat(1)

    enddo
    prog_data_device%dt_init = prog_data%gru(1)%hru(1)%var(iLookPROG%dt_init)%dat(1)

  end subroutine allocate_device_prog_data  
    subroutine initialize_device_prog_data(prog_data_device, prog_data, nGRU, nSoil)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    type(var_dlength),intent(in) :: prog_data
    integer(i4b),intent(in) :: nGRU,nSoil

    integer(i4b) :: nLayers
    integer(i4b) :: iGRU
    nLayers = nSoil + maxSnowLayers+1

    nLayers = size(prog_data%var(iLookPROG%mLayerTemp)%dat)
    do iGRU=1,nGRU
      prog_data_device%spectralSnowAlbedoDiffuse(1:nSpecBands,iGRU) = prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nSpecBands)
      prog_data_device%mLayerTemp(1:nLayers,iGRU) = prog_data%var(iLookPROG%mLayerTemp)%dat(1:nLayers)
      prog_data_device%mLayerVolFracIce(1:nLayers,iGRU) = prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nLayers)
      prog_data_device%mLayerVolFracLiq(1:nLayers,iGRU) = prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nLayers)
      prog_data_device%mLayerVolFracWat(1:nLayers,iGRU) = prog_data%var(iLookPROG%mLayerVolFracWat)%dat(1:nLayers)
      prog_data_device%mLayerMatricHead(1:nSoil,iGRU) = prog_data%var(iLookPROG%mLayerMatricHead)%dat(1:nSoil)
      prog_data_device%mLayerEnthalpy(1:nLayers,iGRU) = prog_data%var(iLookPROG%mLayerEnthalpy)%dat(1:nLayers)
      prog_data_device%mLayerDepth(1:nLayers,iGRU) = prog_data%var(iLookPROG%mLayerDepth)%dat(1:nLayers)
      prog_data_device%mLayerHeight(0:nLayers,iGRU) = prog_data%var(iLookPROG%mLayerHeight)%dat(0:nLayers)
      prog_data_device%iLayerHeight(0:nLayers,iGRU) = prog_data%var(iLookPROG%iLayerHeight)%dat(0:nLayers)

    enddo
    prog_data_device%dt_init = prog_data%var(iLookPROG%dt_init)%dat(1)
    prog_data_device%scalarCanopyIce = prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)
    prog_data_device%scalarCanopyLiq = prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)
    prog_data_device%scalarCanopyWat = prog_data%var(iLookPROG%scalarCanopyWat)%dat(1)
    prog_data_device%scalarCanairTemp = prog_data%var(iLookPROG%scalarCanairTemp)%dat(1)
    prog_data_device%scalarCanopyTemp = prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1)
    prog_data_device%scalarSnowAlbedo = prog_data%var(iLookPROG%scalarSnowAlbedo)%dat(1)
    prog_data_device%scalarSnowDepth = prog_data%var(iLookPROG%scalarSnowDepth)%dat(1)
    prog_data_device%scalarSWE = prog_data%var(iLookPROG%scalarSWE)%dat(1)
    prog_data_device%scalarSfcMeltPond = prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1)
    prog_data_device%scalarCanairEnthalpy = prog_data%var(iLookPROG%scalarCanairEnthalpy)%dat(1)
    prog_data_device%scalarCanopyEnthalpy = prog_data%var(iLookPROG%scalarCanopyEnthalpy)%dat(1)
    prog_data_device%scalarAquiferStorage = prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1)
    prog_data_device%scalarSurfaceTemp = prog_data%var(iLookPROG%scalarSurfaceTemp)%dat(1)

  end subroutine   

  subroutine allocate_device_prog_temp(prog_data_device, nGRU, nSoil)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    integer(i4b),intent(in) :: nGRU,nSoil

    integer(i4b) :: nLayers
    nLayers = nSoil + maxSnowLayers+1
    ! state variables for vegetation
    allocate(prog_data_device%scalarCanopyIce(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyLiq(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyWat(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanairTemp(nGRU))                ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyTemp(nGRU))                ! get_ixVarType('scalarv')
    ! state variables for snow
    allocate(prog_data_device%spectralSnowAlbedoDiffuse(nSpecBands,nGRU))     ! get_ixVarType('wLength')
    allocate(prog_data_device%scalarSnowAlbedo(nGRU))                ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSnowDepth(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSWE(nGRU))                       ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSfcMeltPond(nGRU))               ! get_ixVarType('scalarv')
    ! define state variables for the snow+soil domain
    allocate(prog_data_device%mLayerTemp(nLayers,nGRU))                    ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerVolFracIce(nLayers,nGRU))              ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerVolFracLiq(nLayers,nGRU))              ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerVolFracWat(nLayers,nGRU))              ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerMatricHead(nSoil,nGRU))              ! get_ixVarType('midSoil')
    ! enthalpy
    allocate(prog_data_device%scalarCanairEnthalpy(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyEnthalpy(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%mLayerEnthalpy(nLayers,nGRU))                ! get_ixVarType('midToto')
    ! other state variables
    allocate(prog_data_device%scalarAquiferStorage(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSurfaceTemp(nGRU))               ! get_ixVarType('scalarv')
    ! define coordinate variables
    allocate(prog_data_device%mLayerDepth(nLayers,nGRU))                   ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerHeight(0:nLayers,nGRU))                  ! get_ixVarType('ifcToto')
    allocate(prog_data_device%iLayerHeight(0:nLayers,nGRU))                  ! get_ixVarType('ifcToto')
  end subroutine allocate_device_prog_temp  


  subroutine finalize_device_prog_data(prog_data_device, prog_data,indxStruct)
    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    type(gru_hru_doubleVec),intent(inout) :: prog_data
    type(gru_hru_intVec),intent(inout) :: indxStruct
    integer(i4b) :: nLayers,nSoil
    ! integer(i4b),intent(in) :: nGRU,nLayers,nSoil

    integer(i4b) :: iGRU,nGRU
    nGRU = size(prog_data%gru)

    do iGRU=1,nGRU
      nLayers = indxStruct%gru(iGRU)%hru(1)%var(iLookINDEX%nLayers)%dat(1)
      nSoil = indxStruct%gru(iGRU)%hru(1)%var(iLookINDEX%nSoil)%dat(1)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nSpecBands) = prog_data_device%spectralSnowAlbedoDiffuse(1:nSpecBands,iGRU)

    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerTemp)%dat).ne.nLayers) deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerTemp)%dat)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerTemp)%dat = prog_data_device%mLayerTemp(1:nLayers,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracIce)%dat).ne.nLayers) deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracIce)%dat)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracIce)%dat = prog_data_device%mLayerVolFracIce(1:nLayers,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracLiq)%dat).ne.nLayers) deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracLiq)%dat)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracLiq)%dat = prog_data_device%mLayerVolFracLiq(1:nLayers,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracWat)%dat).ne.nLayers) deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracWat)%dat)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracWat)%dat = prog_data_device%mLayerVolFracWat(1:nLayers,iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerMatricHead)%dat(1:nSoil) = prog_data_device%mLayerMatricHead(1:nSoil,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerEnthalpy)%dat).ne.nLayers) deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerEnthalpy)%dat)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerEnthalpy)%dat = prog_data_device%mLayerEnthalpy(1:nLayers,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerDepth)%dat).ne.nLayers) deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerDepth)%dat)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerDepth)%dat = prog_data_device%mLayerDepth(1:nLayers,iGRU)
    print*, nLayers, prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerDepth)%dat
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerHeight)%dat).ne.nLayers+1) then
      deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerHeight)%dat)
      allocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerHeight)%dat(0:nLayers))
    end if
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerHeight)%dat(0:nLayers) = prog_data_device%mLayerHeight(0:nLayers,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%iLayerHeight)%dat).ne.nLayers+1) then
      deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%iLayerHeight)%dat)
      allocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%iLayerHeight)%dat(0:nLayers))
    end if
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%iLayerHeight)%dat(0:nLayers) = prog_data_device%iLayerHeight(0:nLayers,iGRU)

    
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%dt_init)%dat(1) = prog_data_device%dt_init
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyIce)%dat(1) = prog_data_device%scalarCanopyIce(iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyLiq)%dat(1) = prog_data_device%scalarCanopyLiq(iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyWat)%dat(1) = prog_data_device%scalarCanopyWat(iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanairTemp)%dat(1) = prog_data_device%scalarCanairTemp(iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyTemp)%dat(1) = prog_data_device%scalarCanopyTemp(iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = prog_data_device%scalarSnowAlbedo(iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSnowDepth)%dat(1) = prog_data_device%scalarSnowDepth(iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSWE)%dat(1) = prog_data_device%scalarSWE(iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSfcMeltPond)%dat(1) = prog_data_device%scalarSfcMeltPond(iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanairEnthalpy)%dat(1) = prog_data_device%scalarCanairEnthalpy(iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyEnthalpy)%dat(1) = prog_data_device%scalarCanopyEnthalpy(iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarAquiferStorage)%dat(1) = prog_data_device%scalarAquiferStorage(iGRU)
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSurfaceTemp)%dat(1) = prog_data_device%scalarSurfaceTemp(iGRU)
  end do

  end subroutine finalize_device_prog_data  

  subroutine deallocate_device_prog_data(prog_data_device)
    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    ! integer(i4b),intent(in) :: nGRU,nLayers,nSoil

    deallocate(prog_data_device%scalarCanopyIce)
    deallocate(prog_data_device%scalarCanopyLiq)
    deallocate(prog_data_device%scalarCanopyWat)
    deallocate(prog_data_device%scalarCanairTemp)
    deallocate(prog_data_device%scalarCanopyTemp)
    deallocate(prog_data_device%spectralSnowAlbedoDiffuse)
    deallocate(prog_data_device%scalarSnowAlbedo)
    deallocate(prog_data_device%scalarSnowDepth)
    deallocate(prog_data_device%scalarSWE)
    deallocate(prog_data_device%scalarSfcMeltPond)
    deallocate(prog_data_device%mLayerTemp)
    deallocate(prog_data_device%mLayerVolFracIce)
    deallocate(prog_data_device%mLayerVolFracLiq)
    deallocate(prog_data_device%mLayerVolFracWat)
    deallocate(prog_data_device%mLayerMatricHead)
    deallocate(prog_data_device%scalarCanairEnthalpy)
    deallocate(prog_data_device%scalarCanopyEnthalpy)
    deallocate(prog_data_device%mLayerEnthalpy)
    deallocate(prog_data_device%scalarAquiferStorage)
    deallocate(prog_data_device%scalarSurfaceTemp)
    deallocate(prog_data_device%mLayerDepth)
    deallocate(prog_data_device%mLayerHeight)
    deallocate(prog_data_device%iLayerHeight)
  
  end subroutine deallocate_device_prog_data
      

  subroutine allocate_device_flux_data(flux_data_device, flux_data,nSoil,nGRU)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    type(flux_data_device), intent(inout) :: flux_data_device
    type(gru_hru_doubleVec), intent(in) :: flux_data
    integer(i4b),intent(in) :: nSoil,nGRU

    integer(i4b) :: iGRU, iLayer

    integer(i4b) :: nSnow, nLayers
    nSnow = maxSnowLayers+1
    nLayers = nSoil + nSnow
      
    flux_data_device%ixSpectralIncomingDirect_start = flux_data_device%numScalarFluxData + 1
    flux_data_device%ixSpectralIncomingDirect_end = flux_data_device%ixSpectralIncomingDirect_start + nSpecBands - 1
    flux_data_device%ixSpectralIncomingDiffuse_start = flux_data_device%ixSpectralIncomingDirect_end + 1
    flux_data_device%ixSpectralIncomingDiffuse_end = flux_data_device%ixSpectralIncomingDiffuse_start + nSpecBands - 1
    flux_data_device%ixSpectralBelowCanopyDirect_start = flux_data_device%ixSpectralIncomingDiffuse_end + 1
    flux_data_device%ixSpectralBelowCanopyDirect_end = flux_data_device%ixSpectralBelowCanopyDirect_start + nSpecBands - 1
    flux_data_device%ixSpectralBelowCanopyDiffuse_start = flux_data_device%ixSpectralBelowCanopyDirect_end + 1
    flux_data_device%ixSpectralBelowCanopyDiffuse_end = flux_data_device%ixSpectralBelowCanopyDiffuse_start + nSpecBands - 1
    flux_data_device%ixMLayerTranspire_start = flux_data_device%ixSpectralBelowCanopyDiffuse_end + 1
    flux_data_device%ixMLayerTranspire_end = flux_data_device%ixMLayerTranspire_start + nSoil - 1
    flux_data_device%ixILayerConductiveFlux_start = flux_data_device%ixMLayerTranspire_end + 1
    flux_data_device%ixILayerConductiveFlux_end = flux_data_device%ixILayerConductiveFlux_start + nLayers
    flux_data_device%ixILayerAdvectiveFlux_start = flux_data_device%ixILayerConductiveFlux_end + 1
    flux_data_device%ixILayerAdvectiveFlux_end = flux_data_device%ixILayerAdvectiveFlux_start + nLayers
    flux_data_device%ixILayerNrgFlux_start = flux_data_device%ixILayerAdvectiveFlux_end + 1
    flux_data_device%ixILayerNrgFlux_end = flux_data_device%ixILayerNrgFlux_start + nLayers
    flux_data_device%ixMLayerNrgFlux_start = flux_data_device%ixILayerNrgFlux_end + 1
    flux_data_device%ixMLayerNrgFlux_end = flux_data_device%ixMLayerNrgFlux_start + nLayers-1
    flux_data_device%ixILayerLiqFluxSnow_start = flux_data_device%ixMLayerNrgFlux_end + 1
    flux_data_device%ixILayerLiqFluxSnow_end = flux_data_device%ixILayerLiqFluxSnow_start + nSnow
    flux_data_device%ixMLayerLiqFluxSnow_start = flux_data_device%ixILayerLiqFluxSnow_end + 1
    flux_data_device%ixMLayerLiqFluxSnow_end = flux_data_device%ixMLayerLiqFluxSnow_start + nSnow - 1
    flux_data_device%ixMLayerSatHydCondMP_start = flux_data_device%ixMLayerLiqFluxSnow_end + 1
    flux_data_device%ixMLayerSatHydCondMP_end = flux_data_device%ixMLayerSatHydCondMP_start + nSoil - 1
    flux_data_device%ixMLayerSatHydCond_start = flux_data_device%ixMLayerSatHydCondMP_end + 1
    flux_data_device%ixMLayerSatHydCond_end = flux_data_device%ixMLayerSatHydCond_start + nSoil - 1
    flux_data_device%ixILayerSatHydCond_start = flux_data_device%ixMLayerSatHydCond_end + 1
    flux_data_device%ixILayerSatHydCond_end = flux_data_device%ixILayerSatHydCond_start + nSoil
    flux_data_device%ixMLayerHydCond_start = flux_data_device%ixILayerSatHydCond_end + 1
    flux_data_device%ixMLayerHydCond_end = flux_data_device%ixMLayerHydCond_start + nSoil - 1
    flux_data_device%ixILayerLiqFluxSoil_start = flux_data_device%ixMLayerHydCond_end + 1
    flux_data_device%ixILayerLiqFluxSoil_end = flux_data_device%ixILayerLiqFluxSoil_start + nSoil
    flux_data_device%ixMLayerLiqFluxSoil_start = flux_data_device%ixILayerLiqFluxSoil_end + 1
    flux_data_device%ixMLayerLiqFluxSoil_end = flux_data_device%ixMLayerLiqFluxSoil_start + nSoil - 1
    flux_data_device%ixMLayerBaseflow_start = flux_data_device%ixMLayerLiqFluxSoil_end + 1
    flux_data_device%ixMLayerBaseflow_end = flux_data_device%ixMLayerBaseflow_start + nSoil - 1
    flux_data_device%ixMLayerColumnInflow_start = flux_data_device%ixMLayerBaseflow_end + 1
    flux_data_device%ixMLayerColumnInflow_end = flux_data_device%ixMLayerColumnInflow_start + nSoil - 1
    flux_data_device%ixMLayerColumnOutflow_start = flux_data_device%ixMLayerColumnInflow_end + 1
    flux_data_device%ixMLayerColumnOutflow_end = flux_data_device%ixMLayerColumnOutflow_start + nSoil - 1
    flux_data_device%numFluxData = flux_data_device%ixMLayerColumnOutflow_end


    allocate(flux_data_device%data(flux_data_device%numFluxData,nGRU))
    do iGRU=1,nGRU
    nLayers = size(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerNrgFlux)%dat)
    nSnow = size(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerLiqFluxSnow)%dat)
      flux_data_device%data(flux_data_device%ixScalarCanairNetNrgFlux,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarCanopyNetNrgFlux,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarGroundNetNrgFlux,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarCanopyNetLiqFlux,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarRainfall,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarRainfall)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarSnowfall,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSnowfall)%dat(1)
      flux_data_device%data(flux_data_device%ixSpectralIncomingDirect_start:flux_data_device%ixSpectralIncomingDirect_end,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%spectralIncomingDirect)%dat
      flux_data_device%data(flux_data_device%ixSpectralIncomingDiffuse_start:flux_data_device%ixSpectralIncomingDiffuse_end,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%spectralIncomingDiffuse)%dat
      flux_data_device%data(flux_data_device%ixScalarCanopySunlitPAR,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopySunlitPAR)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarCanopyShadedPAR,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyShadedPAR)%dat(1)
      flux_data_device%data(flux_data_device%ixSpectralBelowCanopyDirect_start:flux_data_device%ixSpectralBelowCanopyDirect_end,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%spectralBelowCanopyDirect)%dat
      flux_data_device%data(flux_data_device%ixSpectralBelowCanopyDiffuse_start:flux_data_device%ixSpectralBelowCanopyDiffuse_end,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%spectralBelowCanopyDiffuse)%dat
      flux_data_device%data(flux_data_device%ixScalarBelowCanopySolar,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarBelowCanopySolar)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarCanopyAbsorbedSolar,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyAbsorbedSolar)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarGroundAbsorbedSolar,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWRadCanopy,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadCanopy)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWRadGround,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadGround)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWRadUbound2Canopy,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadUbound2Canopy)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWRadUbound2Ground,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadUbound2Ground)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWRadUbound2Ubound,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadUbound2Ubound)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWRadCanopy2Ubound,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadCanopy2Ubound)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWRadCanopy2Ground,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadCanopy2Ground)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWRadCanopy2Canopy,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadCanopy2Canopy)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWRadGround2Ubound,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadGround2Ubound)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWRadGround2Canopy,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadGround2Canopy)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWNetCanopy,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWNetCanopy)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWNetGround,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWNetGround)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLWNetUbound,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWNetUbound)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarEddyDiffusCanopyTop,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarEddyDiffusCanopyTop)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarFrictionVelocity,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarFrictionVelocity)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarWindspdCanopyTop,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarWindspdCanopyTop)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarWindspdCanopyBottom,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarWindspdCanopyBottom)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarGroundResistance,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundResistance)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarCanopyResistance,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyResistance)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLeafResistance,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLeafResistance)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarSoilResistance,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSoilResistance)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarSenHeatTotal,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSenHeatTotal)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarSenHeatCanopy,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarSenHeatGround,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSenHeatGround)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLatHeatTotal,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLatHeatTotal)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLatHeatCanopyEvap,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLatHeatCanopyTrans,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLatHeatCanopyTrans)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarLatHeatGround,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLatHeatGround)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarCanopyAdvectiveHeatFlux,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyAdvectiveHeatFlux)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarGroundAdvectiveHeatFlux,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundAdvectiveHeatFlux)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarCanopySublimation,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopySublimation)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarSnowSublimation,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSnowSublimation)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarStomResistSunlit,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarStomResistSunlit)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarStomResistShaded,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarStomResistShaded)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarPhotosynthesisSunlit,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarPhotosynthesisSunlit)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarPhotosynthesisShaded,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarPhotosynthesisShaded)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarCanopyTranspiration,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyTranspiration)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarCanopyEvaporation,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarGroundEvaporation,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundEvaporation)%dat(1)
      flux_data_device%data(flux_data_device%ixMLayerTranspire_start:flux_data_device%ixMLayerTranspire_end,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerTranspire)%dat
      flux_data_device%data(flux_data_device%ixScalarThroughfallSnow,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarThroughfallSnow)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarThroughfallRain,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarThroughfallRain)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarCanopySnowUnloading,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarCanopyLiqDrainage,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)
      flux_data_device%data(flux_data_device%ixScalarCanopyMeltFreeze,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyMeltFreeze)%dat(1)
do iLayer=0,nLayers
  flux_data_device%data(flux_data_device%ixILayerConductiveFlux_start+iLayer,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerConductiveFlux)%dat(iLayer)
  flux_data_device%data(flux_data_device%ixILayerAdvectiveFlux_start+iLayer,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerAdvectiveFlux)%dat(iLayer)
  flux_data_device%data(flux_data_device%ixILayerNrgFlux_start+iLayer,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerNrgFlux)%dat(iLayer)
end do
do iLayer=1,nLayers
  flux_data_device%data(flux_data_device%ixMLayerNrgFlux_start+iLayer-1,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerNrgFlux)%dat(iLayer)
end do
flux_data_device%data(flux_data_device%ixScalarSnowDrainage,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSnowDrainage)%dat(1)
do iLayer=0,nSnow
  flux_data_device%data(flux_data_device%ixILayerLiqFluxSnow_start + iLayer,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerLiqFluxSnow)%dat(iLayer)
end do
do iLayer=1,nSnow
  flux_data_device%data(flux_data_device%ixMLayerLiqFluxSnow_start+iLayer-1,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerLiqFluxSnow)%dat(iLayer)
end do

flux_data_device%data(flux_data_device%ixScalarRainPlusMelt,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarRainPlusMelt)%dat(1)
flux_data_device%data(flux_data_device%ixScalarMaxInfilRate,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarMaxInfilRate)%dat(1)
flux_data_device%data(flux_data_device%ixScalarInfiltration,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarInfiltration)%dat(1)
flux_data_device%data(flux_data_device%ixScalarExfiltration,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarExfiltration)%dat(1)
flux_data_device%data(flux_data_device%ixScalarSurfaceRunoff,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)
flux_data_device%data(flux_data_device%ixScalarSurfaceRunoff_IE,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSurfaceRunoff_IE)%dat(1)
flux_data_device%data(flux_data_device%ixScalarSurfaceRunoff_SE,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSurfaceRunoff_SE)%dat(1)

do iLayer=1,nSoil
  flux_data_device%data(flux_data_device%ixMLayerSatHydCondMP_start+iLayer-1,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerSatHydCondMP)%dat(iLayer)
  flux_data_device%data(flux_data_device%ixMLayerSatHydCond_start+iLayer-1,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerSatHydCond)%dat(iLayer)
end do
do iLayer=0,nSoil
  flux_data_device%data(flux_data_device%ixILayerSatHydCond_start+iLayer,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerSatHydCond)%dat(iLayer)
end do
do iLayer=1,nSoil
flux_data_device%data(flux_data_device%ixmLayerHydCond_start + iLayer-1,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerHydCond)%dat(iLayer)
end do
do iLayer=0,nSoil
flux_data_device%data(flux_data_device%ixiLayerLiqFluxSoil_start + iLayer,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerLiqFluxSoil)%dat(iLayer)
end do
do iLayer=1,nSoil
flux_data_device%data(flux_data_device%ixmLayerLiqFluxSoil_start + iLayer-1,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerLiqFluxSoil)%dat(iLayer)
flux_data_device%data(flux_data_device%ixmLayerBaseflow_start + iLayer-1,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerBaseflow)%dat(iLayer)
flux_data_device%data(flux_data_device%ixmLayerColumnInflow_start + iLayer-1,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerColumnInflow)%dat(iLayer)
flux_data_device%data(flux_data_device%ixmLayerColumnOutflow_start + iLayer-1,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerColumnOutflow)%dat(iLayer)
end do
flux_data_device%data(flux_data_device%ixScalarSoilBaseflow,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSoilBaseflow)%dat(1)
flux_data_device%data(flux_data_device%ixScalarSoilDrainage,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSoilDrainage)%dat(1)
flux_data_device%data(flux_data_device%ixScalarAquiferRecharge,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarAquiferRecharge)%dat(1)
flux_data_device%data(flux_data_device%ixScalarAquiferTranspire,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarAquiferTranspire)%dat(1)
flux_data_device%data(flux_data_device%ixScalarAquiferBaseflow,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarAquiferBaseflow)%dat(1)
flux_data_device%data(flux_data_device%ixScalarTotalET,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarTotalET)%dat(1)
flux_data_device%data(flux_data_device%ixScalarTotalRunoff,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarTotalRunoff)%dat(1)
flux_data_device%data(flux_data_device%ixScalarNetRadiation,iGRU) = flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarNetRadiation)%dat(1)
    end do

  end subroutine allocate_device_flux_data

  subroutine allocate_device_flux_Mask(flux_data_device,nSoil,nGRU)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    type(flux_mask_device), intent(inout) :: flux_data_device
    integer(i4b),intent(in) :: nSoil,nGRU
    integer(i4b) :: nSnow, nLayers

    integer(i4b) :: iGRU, iLayer
    nSnow = maxSnowLayers+1
    nLayers = nSoil + nSnow
    flux_data_device%ixSpectralIncomingDirect_start = flux_data_device%numScalarFluxData + 1
    flux_data_device%ixSpectralIncomingDirect_end = flux_data_device%ixSpectralIncomingDirect_start + nSpecBands - 1
    flux_data_device%ixSpectralIncomingDiffuse_start = flux_data_device%ixSpectralIncomingDirect_end + 1
    flux_data_device%ixSpectralIncomingDiffuse_end = flux_data_device%ixSpectralIncomingDiffuse_start + nSpecBands - 1
    flux_data_device%ixSpectralBelowCanopyDirect_start = flux_data_device%ixSpectralIncomingDiffuse_end + 1
    flux_data_device%ixSpectralBelowCanopyDirect_end = flux_data_device%ixSpectralBelowCanopyDirect_start + nSpecBands - 1
    flux_data_device%ixSpectralBelowCanopyDiffuse_start = flux_data_device%ixSpectralBelowCanopyDirect_end + 1
    flux_data_device%ixSpectralBelowCanopyDiffuse_end = flux_data_device%ixSpectralBelowCanopyDiffuse_start + nSpecBands - 1
    flux_data_device%ixMLayerTranspire_start = flux_data_device%ixSpectralBelowCanopyDiffuse_end + 1
    flux_data_device%ixMLayerTranspire_end = flux_data_device%ixMLayerTranspire_start + nSoil - 1
    flux_data_device%ixILayerConductiveFlux_start = flux_data_device%ixMLayerTranspire_end + 1
    flux_data_device%ixILayerConductiveFlux_end = flux_data_device%ixILayerConductiveFlux_start + nLayers
    flux_data_device%ixILayerAdvectiveFlux_start = flux_data_device%ixILayerConductiveFlux_end + 1
    flux_data_device%ixILayerAdvectiveFlux_end = flux_data_device%ixILayerAdvectiveFlux_start + nLayers
    flux_data_device%ixILayerNrgFlux_start = flux_data_device%ixILayerAdvectiveFlux_end + 1
    flux_data_device%ixILayerNrgFlux_end = flux_data_device%ixILayerNrgFlux_start + nLayers
    flux_data_device%ixMLayerNrgFlux_start = flux_data_device%ixILayerNrgFlux_end + 1
    flux_data_device%ixMLayerNrgFlux_end = flux_data_device%ixMLayerNrgFlux_start + nLayers-1
    flux_data_device%ixILayerLiqFluxSnow_start = flux_data_device%ixMLayerNrgFlux_end + 1
    flux_data_device%ixILayerLiqFluxSnow_end = flux_data_device%ixILayerLiqFluxSnow_start + nSnow
    flux_data_device%ixMLayerLiqFluxSnow_start = flux_data_device%ixILayerLiqFluxSnow_end + 1
    flux_data_device%ixMLayerLiqFluxSnow_end = flux_data_device%ixMLayerLiqFluxSnow_start + nSnow - 1
    flux_data_device%ixMLayerSatHydCondMP_start = flux_data_device%ixMLayerLiqFluxSnow_end + 1
    flux_data_device%ixMLayerSatHydCondMP_end = flux_data_device%ixMLayerSatHydCondMP_start + nSoil - 1
    flux_data_device%ixMLayerSatHydCond_start = flux_data_device%ixMLayerSatHydCondMP_end + 1
    flux_data_device%ixMLayerSatHydCond_end = flux_data_device%ixMLayerSatHydCond_start + nSoil - 1
    flux_data_device%ixILayerSatHydCond_start = flux_data_device%ixMLayerSatHydCond_end + 1
    flux_data_device%ixILayerSatHydCond_end = flux_data_device%ixILayerSatHydCond_start + nSoil
    flux_data_device%ixMLayerHydCond_start = flux_data_device%ixILayerSatHydCond_end + 1
    flux_data_device%ixMLayerHydCond_end = flux_data_device%ixMLayerHydCond_start + nSoil - 1
    flux_data_device%ixILayerLiqFluxSoil_start = flux_data_device%ixMLayerHydCond_end + 1
    flux_data_device%ixILayerLiqFluxSoil_end = flux_data_device%ixILayerLiqFluxSoil_start + nSoil
    flux_data_device%ixMLayerLiqFluxSoil_start = flux_data_device%ixILayerLiqFluxSoil_end + 1
    flux_data_device%ixMLayerLiqFluxSoil_end = flux_data_device%ixMLayerLiqFluxSoil_start + nSoil - 1
    flux_data_device%ixMLayerBaseflow_start = flux_data_device%ixMLayerLiqFluxSoil_end + 1
    flux_data_device%ixMLayerBaseflow_end = flux_data_device%ixMLayerBaseflow_start + nSoil - 1
    flux_data_device%ixMLayerColumnInflow_start = flux_data_device%ixMLayerBaseflow_end + 1
    flux_data_device%ixMLayerColumnInflow_end = flux_data_device%ixMLayerColumnInflow_start + nSoil - 1
    flux_data_device%ixMLayerColumnOutflow_start = flux_data_device%ixMLayerColumnInflow_end + 1
    flux_data_device%ixMLayerColumnOutflow_end = flux_data_device%ixMLayerColumnOutflow_start + nSoil - 1
    flux_data_device%numFluxData = flux_data_device%ixMLayerColumnOutflow_end


    allocate(flux_data_device%data(flux_data_device%numFluxData,nGRU))

  end subroutine allocate_device_flux_Mask

    subroutine allocate_device_flux_Count(flux_data_device,nSoil,nGRU)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    type(flux_count_device), intent(inout) :: flux_data_device
    integer(i4b),intent(in) :: nSoil,nGRU

    integer(i4b) :: iGRU, iLayer
    integer(i4b) :: nSnow, nLayers
    nSnow = maxSnowLayers+1
    nLayers = nSnow + nSoil
    flux_data_device%ixSpectralIncomingDirect_start = flux_data_device%numScalarFluxData + 1
    flux_data_device%ixSpectralIncomingDirect_end = flux_data_device%ixSpectralIncomingDirect_start + nSpecBands - 1
    flux_data_device%ixSpectralIncomingDiffuse_start = flux_data_device%ixSpectralIncomingDirect_end + 1
    flux_data_device%ixSpectralIncomingDiffuse_end = flux_data_device%ixSpectralIncomingDiffuse_start + nSpecBands - 1
    flux_data_device%ixSpectralBelowCanopyDirect_start = flux_data_device%ixSpectralIncomingDiffuse_end + 1
    flux_data_device%ixSpectralBelowCanopyDirect_end = flux_data_device%ixSpectralBelowCanopyDirect_start + nSpecBands - 1
    flux_data_device%ixSpectralBelowCanopyDiffuse_start = flux_data_device%ixSpectralBelowCanopyDirect_end + 1
    flux_data_device%ixSpectralBelowCanopyDiffuse_end = flux_data_device%ixSpectralBelowCanopyDiffuse_start + nSpecBands - 1
    flux_data_device%ixMLayerTranspire_start = flux_data_device%ixSpectralBelowCanopyDiffuse_end + 1
    flux_data_device%ixMLayerTranspire_end = flux_data_device%ixMLayerTranspire_start + nSoil - 1
    flux_data_device%ixILayerConductiveFlux_start = flux_data_device%ixMLayerTranspire_end + 1
    flux_data_device%ixILayerConductiveFlux_end = flux_data_device%ixILayerConductiveFlux_start + nLayers
    flux_data_device%ixILayerAdvectiveFlux_start = flux_data_device%ixILayerConductiveFlux_end + 1
    flux_data_device%ixILayerAdvectiveFlux_end = flux_data_device%ixILayerAdvectiveFlux_start + nLayers
    flux_data_device%ixILayerNrgFlux_start = flux_data_device%ixILayerAdvectiveFlux_end + 1
    flux_data_device%ixILayerNrgFlux_end = flux_data_device%ixILayerNrgFlux_start + nLayers
    flux_data_device%ixMLayerNrgFlux_start = flux_data_device%ixILayerNrgFlux_end + 1
    flux_data_device%ixMLayerNrgFlux_end = flux_data_device%ixMLayerNrgFlux_start + nLayers-1
    flux_data_device%ixILayerLiqFluxSnow_start = flux_data_device%ixMLayerNrgFlux_end + 1
    flux_data_device%ixILayerLiqFluxSnow_end = flux_data_device%ixILayerLiqFluxSnow_start + nSnow
    flux_data_device%ixMLayerLiqFluxSnow_start = flux_data_device%ixILayerLiqFluxSnow_end + 1
    flux_data_device%ixMLayerLiqFluxSnow_end = flux_data_device%ixMLayerLiqFluxSnow_start + nSnow - 1
    flux_data_device%ixMLayerSatHydCondMP_start = flux_data_device%ixMLayerLiqFluxSnow_end + 1
    flux_data_device%ixMLayerSatHydCondMP_end = flux_data_device%ixMLayerSatHydCondMP_start + nSoil - 1
    flux_data_device%ixMLayerSatHydCond_start = flux_data_device%ixMLayerSatHydCondMP_end + 1
    flux_data_device%ixMLayerSatHydCond_end = flux_data_device%ixMLayerSatHydCond_start + nSoil - 1
    flux_data_device%ixILayerSatHydCond_start = flux_data_device%ixMLayerSatHydCond_end + 1
    flux_data_device%ixILayerSatHydCond_end = flux_data_device%ixILayerSatHydCond_start + nSoil
    flux_data_device%ixMLayerHydCond_start = flux_data_device%ixILayerSatHydCond_end + 1
    flux_data_device%ixMLayerHydCond_end = flux_data_device%ixMLayerHydCond_start + nSoil - 1
    flux_data_device%ixILayerLiqFluxSoil_start = flux_data_device%ixMLayerHydCond_end + 1
    flux_data_device%ixILayerLiqFluxSoil_end = flux_data_device%ixILayerLiqFluxSoil_start + nSoil
    flux_data_device%ixMLayerLiqFluxSoil_start = flux_data_device%ixILayerLiqFluxSoil_end + 1
    flux_data_device%ixMLayerLiqFluxSoil_end = flux_data_device%ixMLayerLiqFluxSoil_start + nSoil - 1
    flux_data_device%ixMLayerBaseflow_start = flux_data_device%ixMLayerLiqFluxSoil_end + 1
    flux_data_device%ixMLayerBaseflow_end = flux_data_device%ixMLayerBaseflow_start + nSoil - 1
    flux_data_device%ixMLayerColumnInflow_start = flux_data_device%ixMLayerBaseflow_end + 1
    flux_data_device%ixMLayerColumnInflow_end = flux_data_device%ixMLayerColumnInflow_start + nSoil - 1
    flux_data_device%ixMLayerColumnOutflow_start = flux_data_device%ixMLayerColumnInflow_end + 1
    flux_data_device%ixMLayerColumnOutflow_end = flux_data_device%ixMLayerColumnOutflow_start + nSoil - 1
    flux_data_device%numFluxData = flux_data_device%ixMLayerColumnOutflow_end


    allocate(flux_data_device%data(flux_data_device%numFluxData,nGRU))

  end subroutine allocate_device_flux_Count

    subroutine allocate_device_flux2state(flux_data_device, flux_data,nSoil)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    type(flux2state_device), intent(inout) :: flux_data_device
    type(flux2state), intent(in) :: flux_data(:)
    integer(i4b),intent(in) :: nSoil

    integer(i4b) :: nSnow, nLayers
    nSnow = maxSnowLayers+1
    nLayers = nSnow + nSoil
    flux_data_device%ixSpectralIncomingDirect_start = flux_data_device%numScalarFluxData + 1
    flux_data_device%ixSpectralIncomingDirect_end = flux_data_device%ixSpectralIncomingDirect_start + nSpecBands - 1
    flux_data_device%ixSpectralIncomingDiffuse_start = flux_data_device%ixSpectralIncomingDirect_end + 1
    flux_data_device%ixSpectralIncomingDiffuse_end = flux_data_device%ixSpectralIncomingDiffuse_start + nSpecBands - 1
    flux_data_device%ixSpectralBelowCanopyDirect_start = flux_data_device%ixSpectralIncomingDiffuse_end + 1
    flux_data_device%ixSpectralBelowCanopyDirect_end = flux_data_device%ixSpectralBelowCanopyDirect_start + nSpecBands - 1
    flux_data_device%ixSpectralBelowCanopyDiffuse_start = flux_data_device%ixSpectralBelowCanopyDirect_end + 1
    flux_data_device%ixSpectralBelowCanopyDiffuse_end = flux_data_device%ixSpectralBelowCanopyDiffuse_start + nSpecBands - 1
    flux_data_device%ixMLayerTranspire_start = flux_data_device%ixSpectralBelowCanopyDiffuse_end + 1
    flux_data_device%ixMLayerTranspire_end = flux_data_device%ixMLayerTranspire_start + nSoil - 1
    flux_data_device%ixILayerConductiveFlux_start = flux_data_device%ixMLayerTranspire_end + 1
    flux_data_device%ixILayerConductiveFlux_end = flux_data_device%ixILayerConductiveFlux_start + nLayers
    flux_data_device%ixILayerAdvectiveFlux_start = flux_data_device%ixILayerConductiveFlux_end + 1
    flux_data_device%ixILayerAdvectiveFlux_end = flux_data_device%ixILayerAdvectiveFlux_start + nLayers
    flux_data_device%ixILayerNrgFlux_start = flux_data_device%ixILayerAdvectiveFlux_end + 1
    flux_data_device%ixILayerNrgFlux_end = flux_data_device%ixILayerNrgFlux_start + nLayers
    flux_data_device%ixMLayerNrgFlux_start = flux_data_device%ixILayerNrgFlux_end + 1
    flux_data_device%ixMLayerNrgFlux_end = flux_data_device%ixMLayerNrgFlux_start + nLayers-1
    flux_data_device%ixILayerLiqFluxSnow_start = flux_data_device%ixMLayerNrgFlux_end + 1
    flux_data_device%ixILayerLiqFluxSnow_end = flux_data_device%ixILayerLiqFluxSnow_start + nSnow
    flux_data_device%ixMLayerLiqFluxSnow_start = flux_data_device%ixILayerLiqFluxSnow_end + 1
    flux_data_device%ixMLayerLiqFluxSnow_end = flux_data_device%ixMLayerLiqFluxSnow_start + nSnow - 1
    flux_data_device%ixMLayerSatHydCondMP_start = flux_data_device%ixMLayerLiqFluxSnow_end + 1
    flux_data_device%ixMLayerSatHydCondMP_end = flux_data_device%ixMLayerSatHydCondMP_start + nSoil - 1
    flux_data_device%ixMLayerSatHydCond_start = flux_data_device%ixMLayerSatHydCondMP_end + 1
    flux_data_device%ixMLayerSatHydCond_end = flux_data_device%ixMLayerSatHydCond_start + nSoil - 1
    flux_data_device%ixILayerSatHydCond_start = flux_data_device%ixMLayerSatHydCond_end + 1
    flux_data_device%ixILayerSatHydCond_end = flux_data_device%ixILayerSatHydCond_start + nSoil
    flux_data_device%ixMLayerHydCond_start = flux_data_device%ixILayerSatHydCond_end + 1
    flux_data_device%ixMLayerHydCond_end = flux_data_device%ixMLayerHydCond_start + nSoil - 1
    flux_data_device%ixILayerLiqFluxSoil_start = flux_data_device%ixMLayerHydCond_end + 1
    flux_data_device%ixILayerLiqFluxSoil_end = flux_data_device%ixILayerLiqFluxSoil_start + nSoil
    flux_data_device%ixMLayerLiqFluxSoil_start = flux_data_device%ixILayerLiqFluxSoil_end + 1
    flux_data_device%ixMLayerLiqFluxSoil_end = flux_data_device%ixMLayerLiqFluxSoil_start + nSoil - 1
    flux_data_device%ixMLayerBaseflow_start = flux_data_device%ixMLayerLiqFluxSoil_end + 1
    flux_data_device%ixMLayerBaseflow_end = flux_data_device%ixMLayerBaseflow_start + nSoil - 1
    flux_data_device%ixMLayerColumnInflow_start = flux_data_device%ixMLayerBaseflow_end + 1
    flux_data_device%ixMLayerColumnInflow_end = flux_data_device%ixMLayerColumnInflow_start + nSoil - 1
    flux_data_device%ixMLayerColumnOutflow_start = flux_data_device%ixMLayerColumnInflow_end + 1
    flux_data_device%ixMLayerColumnOutflow_end = flux_data_device%ixMLayerColumnOutflow_start + nSoil - 1
    flux_data_device%numFluxData = flux_data_device%ixMLayerColumnOutflow_end


    allocate(flux_data_device%state1(flux_data_device%numFluxData))
    allocate(flux_data_device%state2(flux_data_device%numFluxData))

      flux_data_device%state1(flux_data_device%ixScalarCanairNetNrgFlux) = flux_data(iLookFLUX%scalarCanairNetNrgFlux)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopyNetNrgFlux) = flux_data(iLookFLUX%scalarCanopyNetNrgFlux)%state1
      flux_data_device%state1(flux_data_device%ixScalarGroundNetNrgFlux) = flux_data(iLookFLUX%scalarGroundNetNrgFlux)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopyNetLiqFlux) = flux_data(iLookFLUX%scalarCanopyNetLiqFlux)%state1
      flux_data_device%state1(flux_data_device%ixScalarRainfall) = flux_data(iLookFLUX%scalarRainfall)%state1
      flux_data_device%state1(flux_data_device%ixScalarSnowfall) = flux_data(iLookFLUX%scalarSnowfall)%state1
      flux_data_device%state1(flux_data_device%ixSpectralIncomingDirect_start:flux_data_device%ixSpectralIncomingDirect_end) = flux_data(iLookFLUX%spectralIncomingDirect)%state1
      flux_data_device%state1(flux_data_device%ixSpectralIncomingDiffuse_start:flux_data_device%ixSpectralIncomingDiffuse_end) = flux_data(iLookFLUX%spectralIncomingDiffuse)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopySunlitPAR) = flux_data(iLookFLUX%scalarCanopySunlitPAR)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopyShadedPAR) = flux_data(iLookFLUX%scalarCanopyShadedPAR)%state1
      flux_data_device%state1(flux_data_device%ixSpectralBelowCanopyDirect_start:flux_data_device%ixSpectralBelowCanopyDirect_end) = flux_data(iLookFLUX%spectralBelowCanopyDirect)%state1
      flux_data_device%state1(flux_data_device%ixSpectralBelowCanopyDiffuse_start:flux_data_device%ixSpectralBelowCanopyDiffuse_end) = flux_data(iLookFLUX%spectralBelowCanopyDiffuse)%state1
      flux_data_device%state1(flux_data_device%ixScalarBelowCanopySolar) = flux_data(iLookFLUX%scalarBelowCanopySolar)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopyAbsorbedSolar) = flux_data(iLookFLUX%scalarCanopyAbsorbedSolar)%state1
      flux_data_device%state1(flux_data_device%ixScalarGroundAbsorbedSolar) = flux_data(iLookFLUX%scalarGroundAbsorbedSolar)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWRadCanopy) = flux_data(iLookFLUX%scalarLWRadCanopy)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWRadGround) = flux_data(iLookFLUX%scalarLWRadGround)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWRadUbound2Canopy) = flux_data(iLookFLUX%scalarLWRadUbound2Canopy)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWRadUbound2Ground) = flux_data(iLookFLUX%scalarLWRadUbound2Ground)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWRadUbound2Ubound) = flux_data(iLookFLUX%scalarLWRadUbound2Ubound)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWRadCanopy2Ubound) = flux_data(iLookFLUX%scalarLWRadCanopy2Ubound)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWRadCanopy2Ground) = flux_data(iLookFLUX%scalarLWRadCanopy2Ground)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWRadCanopy2Canopy) = flux_data(iLookFLUX%scalarLWRadCanopy2Canopy)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWRadGround2Ubound) = flux_data(iLookFLUX%scalarLWRadGround2Ubound)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWRadGround2Canopy) = flux_data(iLookFLUX%scalarLWRadGround2Canopy)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWNetCanopy) = flux_data(iLookFLUX%scalarLWNetCanopy)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWNetGround) = flux_data(iLookFLUX%scalarLWNetGround)%state1
      flux_data_device%state1(flux_data_device%ixScalarLWNetUbound) = flux_data(iLookFLUX%scalarLWNetUbound)%state1
      flux_data_device%state1(flux_data_device%ixScalarEddyDiffusCanopyTop) = flux_data(iLookFLUX%scalarEddyDiffusCanopyTop)%state1
      flux_data_device%state1(flux_data_device%ixScalarFrictionVelocity) = flux_data(iLookFLUX%scalarFrictionVelocity)%state1
      flux_data_device%state1(flux_data_device%ixScalarWindspdCanopyTop) = flux_data(iLookFLUX%scalarWindspdCanopyTop)%state1
      flux_data_device%state1(flux_data_device%ixScalarWindspdCanopyBottom) = flux_data(iLookFLUX%scalarWindspdCanopyBottom)%state1
      flux_data_device%state1(flux_data_device%ixScalarGroundResistance) = flux_data(iLookFLUX%scalarGroundResistance)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopyResistance) = flux_data(iLookFLUX%scalarCanopyResistance)%state1
      flux_data_device%state1(flux_data_device%ixScalarLeafResistance) = flux_data(iLookFLUX%scalarLeafResistance)%state1
      flux_data_device%state1(flux_data_device%ixScalarSoilResistance) = flux_data(iLookFLUX%scalarSoilResistance)%state1
      flux_data_device%state1(flux_data_device%ixScalarSenHeatTotal) = flux_data(iLookFLUX%scalarSenHeatTotal)%state1
      flux_data_device%state1(flux_data_device%ixScalarSenHeatCanopy) = flux_data(iLookFLUX%scalarSenHeatCanopy)%state1
      flux_data_device%state1(flux_data_device%ixScalarSenHeatGround) = flux_data(iLookFLUX%scalarSenHeatGround)%state1
      flux_data_device%state1(flux_data_device%ixScalarLatHeatTotal) = flux_data(iLookFLUX%scalarLatHeatTotal)%state1
      flux_data_device%state1(flux_data_device%ixScalarLatHeatCanopyEvap) = flux_data(iLookFLUX%scalarLatHeatCanopyEvap)%state1
      flux_data_device%state1(flux_data_device%ixScalarLatHeatCanopyTrans) = flux_data(iLookFLUX%scalarLatHeatCanopyTrans)%state1
      flux_data_device%state1(flux_data_device%ixScalarLatHeatGround) = flux_data(iLookFLUX%scalarLatHeatGround)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopyAdvectiveHeatFlux) = flux_data(iLookFLUX%scalarCanopyAdvectiveHeatFlux)%state1
      flux_data_device%state1(flux_data_device%ixScalarGroundAdvectiveHeatFlux) = flux_data(iLookFLUX%scalarGroundAdvectiveHeatFlux)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopySublimation) = flux_data(iLookFLUX%scalarCanopySublimation)%state1
      flux_data_device%state1(flux_data_device%ixScalarSnowSublimation) = flux_data(iLookFLUX%scalarSnowSublimation)%state1
      flux_data_device%state1(flux_data_device%ixScalarStomResistSunlit) = flux_data(iLookFLUX%scalarStomResistSunlit)%state1
      flux_data_device%state1(flux_data_device%ixScalarStomResistShaded) = flux_data(iLookFLUX%scalarStomResistShaded)%state1
      flux_data_device%state1(flux_data_device%ixScalarPhotosynthesisSunlit) = flux_data(iLookFLUX%scalarPhotosynthesisSunlit)%state1
      flux_data_device%state1(flux_data_device%ixScalarPhotosynthesisShaded) = flux_data(iLookFLUX%scalarPhotosynthesisShaded)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopyTranspiration) = flux_data(iLookFLUX%scalarCanopyTranspiration)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopyEvaporation) = flux_data(iLookFLUX%scalarCanopyEvaporation)%state1
      flux_data_device%state1(flux_data_device%ixScalarGroundEvaporation) = flux_data(iLookFLUX%scalarGroundEvaporation)%state1
      flux_data_device%state1(flux_data_device%ixMLayerTranspire_start:flux_data_device%ixMLayerTranspire_end) = flux_data(iLookFLUX%mLayerTranspire)%state1
      flux_data_device%state1(flux_data_device%ixScalarThroughfallSnow) = flux_data(iLookFLUX%scalarThroughfallSnow)%state1
      flux_data_device%state1(flux_data_device%ixScalarThroughfallRain) = flux_data(iLookFLUX%scalarThroughfallRain)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopySnowUnloading) = flux_data(iLookFLUX%scalarCanopySnowUnloading)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopyLiqDrainage) = flux_data(iLookFLUX%scalarCanopyLiqDrainage)%state1
      flux_data_device%state1(flux_data_device%ixScalarCanopyMeltFreeze) = flux_data(iLookFLUX%scalarCanopyMeltFreeze)%state1
      flux_data_device%state1(flux_data_device%ixILayerConductiveFlux_start:flux_data_device%ixILayerConductiveFlux_end) = flux_data(iLookFLUX%iLayerConductiveFlux)%state1
      flux_data_device%state1(flux_data_device%ixILayerAdvectiveFlux_start:flux_data_device%ixILayerAdvectiveFlux_end) = flux_data(iLookFLUX%iLayerAdvectiveFlux)%state1
      flux_data_device%state1(flux_data_device%ixILayerNrgFlux_start:flux_data_device%ixILayerNrgFlux_end) = flux_data(iLookFLUX%iLayerNrgFlux)%state1
      flux_data_device%state1(flux_data_device%ixMLayerNrgFlux_start:flux_data_device%ixMLayerNrgFlux_end) = flux_data(iLookFLUX%mLayerNrgFlux)%state1
      flux_data_device%state1(flux_data_device%ixScalarSnowDrainage) = flux_data(iLookFLUX%scalarSnowDrainage)%state1
      flux_data_device%state1(flux_data_device%ixILayerLiqFluxSnow_start:flux_data_device%ixILayerLiqFluxSnow_end) = flux_data(iLookFLUX%iLayerLiqFluxSnow)%state1
      flux_data_device%state1(flux_data_device%ixMLayerLiqFluxSnow_start:flux_data_device%ixMLayerLiqFluxSnow_end) = flux_data(iLookFLUX%mLayerLiqFluxSnow)%state1
      flux_data_device%state1(flux_data_device%ixScalarRainPlusMelt) = flux_data(iLookFLUX%scalarRainPlusMelt)%state1
      flux_data_device%state1(flux_data_device%ixScalarMaxInfilRate) = flux_data(iLookFLUX%scalarMaxInfilRate)%state1
      flux_data_device%state1(flux_data_device%ixScalarInfiltration) = flux_data(iLookFLUX%scalarInfiltration)%state1
      flux_data_device%state1(flux_data_device%ixScalarExfiltration) = flux_data(iLookFLUX%scalarExfiltration)%state1
      flux_data_device%state1(flux_data_device%ixScalarSurfaceRunoff) = flux_data(iLookFLUX%scalarSurfaceRunoff)%state1
      flux_data_device%state1(flux_data_device%ixScalarSurfaceRunoff_IE) = flux_data(iLookFLUX%scalarSurfaceRunoff_IE)%state1
      flux_data_device%state1(flux_data_device%ixScalarSurfaceRunoff_SE) = flux_data(iLookFLUX%scalarSurfaceRunoff_SE)%state1
      flux_data_device%state1(flux_data_device%ixMLayerSatHydCondMP_start:flux_data_device%ixMLayerSatHydCondMP_end) = flux_data(iLookFLUX%mLayerSatHydCondMP)%state1
      flux_data_device%state1(flux_data_device%ixMLayerSatHydCond_start:flux_data_device%ixMLayerSatHydCond_end) = flux_data(iLookFLUX%mLayerSatHydCond)%state1
      flux_data_device%state1(flux_data_device%ixILayerSatHydCond_start:flux_data_device%ixILayerSatHydCond_end) = flux_data(iLookFLUX%iLayerSatHydCond)%state1
      flux_data_device%state1(flux_data_device%ixmLayerHydCond_start:flux_data_device%ixmLayerHydCond_end) = flux_data(iLookFLUX%mLayerHydCond)%state1
      flux_data_device%state1(flux_data_device%ixiLayerLiqFluxSoil_start:flux_data_device%ixiLayerLiqFluxSoil_end) = flux_data(iLookFLUX%iLayerLiqFluxSoil)%state1
      flux_data_device%state1(flux_data_device%ixmLayerLiqFluxSoil_start:flux_data_device%ixmLayerLiqFluxSoil_end) = flux_data(iLookFLUX%mLayerLiqFluxSoil)%state1
      flux_data_device%state1(flux_data_device%ixmLayerBaseflow_start:flux_data_device%ixmLayerBaseflow_end) = flux_data(iLookFLUX%mLayerBaseflow)%state1
      flux_data_device%state1(flux_data_device%ixmLayerColumnInflow_start:flux_data_device%ixmLayerColumnInflow_end) = flux_data(iLookFLUX%mLayerColumnInflow)%state1
      flux_data_device%state1(flux_data_device%ixmLayerColumnOutflow_start:flux_data_device%ixmLayerColumnOutflow_end) = flux_data(iLookFLUX%mLayerColumnOutflow)%state1
      flux_data_device%state1(flux_data_device%ixScalarSoilBaseflow) = flux_data(iLookFLUX%scalarSoilBaseflow)%state1
      flux_data_device%state1(flux_data_device%ixScalarSoilDrainage) = flux_data(iLookFLUX%scalarSoilDrainage)%state1
      flux_data_device%state1(flux_data_device%ixScalarAquiferRecharge) = flux_data(iLookFLUX%scalarAquiferRecharge)%state1
      flux_data_device%state1(flux_data_device%ixScalarAquiferTranspire) = flux_data(iLookFLUX%scalarAquiferTranspire)%state1
      flux_data_device%state1(flux_data_device%ixScalarAquiferBaseflow) = flux_data(iLookFLUX%scalarAquiferBaseflow)%state1
      flux_data_device%state1(flux_data_device%ixScalarTotalET) = flux_data(iLookFLUX%scalarTotalET)%state1
      flux_data_device%state1(flux_data_device%ixScalarTotalRunoff) = flux_data(iLookFLUX%scalarTotalRunoff)%state1
      flux_data_device%state1(flux_data_device%ixScalarNetRadiation) = flux_data(iLookFLUX%scalarNetRadiation)%state1

!!! state2

      flux_data_device%state2(flux_data_device%ixScalarCanairNetNrgFlux) = flux_data(iLookFLUX%scalarCanairNetNrgFlux)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopyNetNrgFlux) = flux_data(iLookFLUX%scalarCanopyNetNrgFlux)%state2
      flux_data_device%state2(flux_data_device%ixScalarGroundNetNrgFlux) = flux_data(iLookFLUX%scalarGroundNetNrgFlux)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopyNetLiqFlux) = flux_data(iLookFLUX%scalarCanopyNetLiqFlux)%state2
      flux_data_device%state2(flux_data_device%ixScalarRainfall) = flux_data(iLookFLUX%scalarRainfall)%state2
      flux_data_device%state2(flux_data_device%ixScalarSnowfall) = flux_data(iLookFLUX%scalarSnowfall)%state2
      flux_data_device%state2(flux_data_device%ixSpectralIncomingDirect_start:flux_data_device%ixSpectralIncomingDirect_end) = flux_data(iLookFLUX%spectralIncomingDirect)%state2
      flux_data_device%state2(flux_data_device%ixSpectralIncomingDiffuse_start:flux_data_device%ixSpectralIncomingDiffuse_end) = flux_data(iLookFLUX%spectralIncomingDiffuse)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopySunlitPAR) = flux_data(iLookFLUX%scalarCanopySunlitPAR)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopyShadedPAR) = flux_data(iLookFLUX%scalarCanopyShadedPAR)%state2
      flux_data_device%state2(flux_data_device%ixSpectralBelowCanopyDirect_start:flux_data_device%ixSpectralBelowCanopyDirect_end) = flux_data(iLookFLUX%spectralBelowCanopyDirect)%state2
      flux_data_device%state2(flux_data_device%ixSpectralBelowCanopyDiffuse_start:flux_data_device%ixSpectralBelowCanopyDiffuse_end) = flux_data(iLookFLUX%spectralBelowCanopyDiffuse)%state2
      flux_data_device%state2(flux_data_device%ixScalarBelowCanopySolar) = flux_data(iLookFLUX%scalarBelowCanopySolar)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopyAbsorbedSolar) = flux_data(iLookFLUX%scalarCanopyAbsorbedSolar)%state2
      flux_data_device%state2(flux_data_device%ixScalarGroundAbsorbedSolar) = flux_data(iLookFLUX%scalarGroundAbsorbedSolar)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWRadCanopy) = flux_data(iLookFLUX%scalarLWRadCanopy)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWRadGround) = flux_data(iLookFLUX%scalarLWRadGround)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWRadUbound2Canopy) = flux_data(iLookFLUX%scalarLWRadUbound2Canopy)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWRadUbound2Ground) = flux_data(iLookFLUX%scalarLWRadUbound2Ground)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWRadUbound2Ubound) = flux_data(iLookFLUX%scalarLWRadUbound2Ubound)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWRadCanopy2Ubound) = flux_data(iLookFLUX%scalarLWRadCanopy2Ubound)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWRadCanopy2Ground) = flux_data(iLookFLUX%scalarLWRadCanopy2Ground)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWRadCanopy2Canopy) = flux_data(iLookFLUX%scalarLWRadCanopy2Canopy)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWRadGround2Ubound) = flux_data(iLookFLUX%scalarLWRadGround2Ubound)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWRadGround2Canopy) = flux_data(iLookFLUX%scalarLWRadGround2Canopy)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWNetCanopy) = flux_data(iLookFLUX%scalarLWNetCanopy)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWNetGround) = flux_data(iLookFLUX%scalarLWNetGround)%state2
      flux_data_device%state2(flux_data_device%ixScalarLWNetUbound) = flux_data(iLookFLUX%scalarLWNetUbound)%state2
      flux_data_device%state2(flux_data_device%ixScalarEddyDiffusCanopyTop) = flux_data(iLookFLUX%scalarEddyDiffusCanopyTop)%state2
      flux_data_device%state2(flux_data_device%ixScalarFrictionVelocity) = flux_data(iLookFLUX%scalarFrictionVelocity)%state2
      flux_data_device%state2(flux_data_device%ixScalarWindspdCanopyTop) = flux_data(iLookFLUX%scalarWindspdCanopyTop)%state2
      flux_data_device%state2(flux_data_device%ixScalarWindspdCanopyBottom) = flux_data(iLookFLUX%scalarWindspdCanopyBottom)%state2
      flux_data_device%state2(flux_data_device%ixScalarGroundResistance) = flux_data(iLookFLUX%scalarGroundResistance)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopyResistance) = flux_data(iLookFLUX%scalarCanopyResistance)%state2
      flux_data_device%state2(flux_data_device%ixScalarLeafResistance) = flux_data(iLookFLUX%scalarLeafResistance)%state2
      flux_data_device%state2(flux_data_device%ixScalarSoilResistance) = flux_data(iLookFLUX%scalarSoilResistance)%state2
      flux_data_device%state2(flux_data_device%ixScalarSenHeatTotal) = flux_data(iLookFLUX%scalarSenHeatTotal)%state2
      flux_data_device%state2(flux_data_device%ixScalarSenHeatCanopy) = flux_data(iLookFLUX%scalarSenHeatCanopy)%state2
      flux_data_device%state2(flux_data_device%ixScalarSenHeatGround) = flux_data(iLookFLUX%scalarSenHeatGround)%state2
      flux_data_device%state2(flux_data_device%ixScalarLatHeatTotal) = flux_data(iLookFLUX%scalarLatHeatTotal)%state2
      flux_data_device%state2(flux_data_device%ixScalarLatHeatCanopyEvap) = flux_data(iLookFLUX%scalarLatHeatCanopyEvap)%state2
      flux_data_device%state2(flux_data_device%ixScalarLatHeatCanopyTrans) = flux_data(iLookFLUX%scalarLatHeatCanopyTrans)%state2
      flux_data_device%state2(flux_data_device%ixScalarLatHeatGround) = flux_data(iLookFLUX%scalarLatHeatGround)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopyAdvectiveHeatFlux) = flux_data(iLookFLUX%scalarCanopyAdvectiveHeatFlux)%state2
      flux_data_device%state2(flux_data_device%ixScalarGroundAdvectiveHeatFlux) = flux_data(iLookFLUX%scalarGroundAdvectiveHeatFlux)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopySublimation) = flux_data(iLookFLUX%scalarCanopySublimation)%state2
      flux_data_device%state2(flux_data_device%ixScalarSnowSublimation) = flux_data(iLookFLUX%scalarSnowSublimation)%state2
      flux_data_device%state2(flux_data_device%ixScalarStomResistSunlit) = flux_data(iLookFLUX%scalarStomResistSunlit)%state2
      flux_data_device%state2(flux_data_device%ixScalarStomResistShaded) = flux_data(iLookFLUX%scalarStomResistShaded)%state2
      flux_data_device%state2(flux_data_device%ixScalarPhotosynthesisSunlit) = flux_data(iLookFLUX%scalarPhotosynthesisSunlit)%state2
      flux_data_device%state2(flux_data_device%ixScalarPhotosynthesisShaded) = flux_data(iLookFLUX%scalarPhotosynthesisShaded)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopyTranspiration) = flux_data(iLookFLUX%scalarCanopyTranspiration)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopyEvaporation) = flux_data(iLookFLUX%scalarCanopyEvaporation)%state2
      flux_data_device%state2(flux_data_device%ixScalarGroundEvaporation) = flux_data(iLookFLUX%scalarGroundEvaporation)%state2
      flux_data_device%state2(flux_data_device%ixMLayerTranspire_start:flux_data_device%ixMLayerTranspire_end) = flux_data(iLookFLUX%mLayerTranspire)%state2
      flux_data_device%state2(flux_data_device%ixScalarThroughfallSnow) = flux_data(iLookFLUX%scalarThroughfallSnow)%state2
      flux_data_device%state2(flux_data_device%ixScalarThroughfallRain) = flux_data(iLookFLUX%scalarThroughfallRain)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopySnowUnloading) = flux_data(iLookFLUX%scalarCanopySnowUnloading)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopyLiqDrainage) = flux_data(iLookFLUX%scalarCanopyLiqDrainage)%state2
      flux_data_device%state2(flux_data_device%ixScalarCanopyMeltFreeze) = flux_data(iLookFLUX%scalarCanopyMeltFreeze)%state2
      flux_data_device%state2(flux_data_device%ixILayerConductiveFlux_start:flux_data_device%ixILayerConductiveFlux_end) = flux_data(iLookFLUX%iLayerConductiveFlux)%state2
      flux_data_device%state2(flux_data_device%ixILayerAdvectiveFlux_start:flux_data_device%ixILayerAdvectiveFlux_end) = flux_data(iLookFLUX%iLayerAdvectiveFlux)%state2
      flux_data_device%state2(flux_data_device%ixILayerNrgFlux_start:flux_data_device%ixILayerNrgFlux_end) = flux_data(iLookFLUX%iLayerNrgFlux)%state2
      flux_data_device%state2(flux_data_device%ixMLayerNrgFlux_start:flux_data_device%ixMLayerNrgFlux_end) = flux_data(iLookFLUX%mLayerNrgFlux)%state2
      flux_data_device%state2(flux_data_device%ixScalarSnowDrainage) = flux_data(iLookFLUX%scalarSnowDrainage)%state2
      flux_data_device%state2(flux_data_device%ixILayerLiqFluxSnow_start:flux_data_device%ixILayerLiqFluxSnow_end) = flux_data(iLookFLUX%iLayerLiqFluxSnow)%state2
      flux_data_device%state2(flux_data_device%ixMLayerLiqFluxSnow_start:flux_data_device%ixMLayerLiqFluxSnow_end) = flux_data(iLookFLUX%mLayerLiqFluxSnow)%state2
      flux_data_device%state2(flux_data_device%ixScalarRainPlusMelt) = flux_data(iLookFLUX%scalarRainPlusMelt)%state2
      flux_data_device%state2(flux_data_device%ixScalarMaxInfilRate) = flux_data(iLookFLUX%scalarMaxInfilRate)%state2
      flux_data_device%state2(flux_data_device%ixScalarInfiltration) = flux_data(iLookFLUX%scalarInfiltration)%state2
      flux_data_device%state2(flux_data_device%ixScalarExfiltration) = flux_data(iLookFLUX%scalarExfiltration)%state2
      flux_data_device%state2(flux_data_device%ixScalarSurfaceRunoff) = flux_data(iLookFLUX%scalarSurfaceRunoff)%state2
      flux_data_device%state2(flux_data_device%ixScalarSurfaceRunoff_IE) = flux_data(iLookFLUX%scalarSurfaceRunoff_IE)%state2
      flux_data_device%state2(flux_data_device%ixScalarSurfaceRunoff_SE) = flux_data(iLookFLUX%scalarSurfaceRunoff_SE)%state2
      flux_data_device%state2(flux_data_device%ixMLayerSatHydCondMP_start:flux_data_device%ixMLayerSatHydCondMP_end) = flux_data(iLookFLUX%mLayerSatHydCondMP)%state2
      flux_data_device%state2(flux_data_device%ixMLayerSatHydCond_start:flux_data_device%ixMLayerSatHydCond_end) = flux_data(iLookFLUX%mLayerSatHydCond)%state2
      flux_data_device%state2(flux_data_device%ixILayerSatHydCond_start:flux_data_device%ixILayerSatHydCond_end) = flux_data(iLookFLUX%iLayerSatHydCond)%state2
      flux_data_device%state2(flux_data_device%ixmLayerHydCond_start:flux_data_device%ixmLayerHydCond_end) = flux_data(iLookFLUX%mLayerHydCond)%state2
      flux_data_device%state2(flux_data_device%ixiLayerLiqFluxSoil_start:flux_data_device%ixiLayerLiqFluxSoil_end) = flux_data(iLookFLUX%iLayerLiqFluxSoil)%state2
      flux_data_device%state2(flux_data_device%ixmLayerLiqFluxSoil_start:flux_data_device%ixmLayerLiqFluxSoil_end) = flux_data(iLookFLUX%mLayerLiqFluxSoil)%state2
      flux_data_device%state2(flux_data_device%ixmLayerBaseflow_start:flux_data_device%ixmLayerBaseflow_end) = flux_data(iLookFLUX%mLayerBaseflow)%state2
      flux_data_device%state2(flux_data_device%ixmLayerColumnInflow_start:flux_data_device%ixmLayerColumnInflow_end) = flux_data(iLookFLUX%mLayerColumnInflow)%state2
      flux_data_device%state2(flux_data_device%ixmLayerColumnOutflow_start:flux_data_device%ixmLayerColumnOutflow_end) = flux_data(iLookFLUX%mLayerColumnOutflow)%state2
      flux_data_device%state2(flux_data_device%ixScalarSoilBaseflow) = flux_data(iLookFLUX%scalarSoilBaseflow)%state2
      flux_data_device%state2(flux_data_device%ixScalarSoilDrainage) = flux_data(iLookFLUX%scalarSoilDrainage)%state2
      flux_data_device%state2(flux_data_device%ixScalarAquiferRecharge) = flux_data(iLookFLUX%scalarAquiferRecharge)%state2
      flux_data_device%state2(flux_data_device%ixScalarAquiferTranspire) = flux_data(iLookFLUX%scalarAquiferTranspire)%state2
      flux_data_device%state2(flux_data_device%ixScalarAquiferBaseflow) = flux_data(iLookFLUX%scalarAquiferBaseflow)%state2
      flux_data_device%state2(flux_data_device%ixScalarTotalET) = flux_data(iLookFLUX%scalarTotalET)%state2
      flux_data_device%state2(flux_data_device%ixScalarTotalRunoff) = flux_data(iLookFLUX%scalarTotalRunoff)%state2
      flux_data_device%state2(flux_data_device%ixScalarNetRadiation) = flux_data(iLookFLUX%scalarNetRadiation)%state2

  end subroutine allocate_device_flux2state


  subroutine allocate_device_flux_prev(flux_data_device,nSoil,nGRU)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    type(flux_data_device), intent(inout) :: flux_data_device
    integer(i4b),intent(in) :: nSoil,nGRU

    integer(i4b) :: iGRU, iLayer
    integer(i4b) :: nSnow
    integer(i4b) :: nLayers
    nSnow = maxSnowLayers+1
    nLayers = maxSnowLayers+1 + nSoil
      
    flux_data_device%ixSpectralIncomingDirect_start = flux_data_device%numScalarFluxData + 1
    flux_data_device%ixSpectralIncomingDirect_end = flux_data_device%ixSpectralIncomingDirect_start + nSpecBands - 1
    flux_data_device%ixSpectralIncomingDiffuse_start = flux_data_device%ixSpectralIncomingDirect_end + 1
    flux_data_device%ixSpectralIncomingDiffuse_end = flux_data_device%ixSpectralIncomingDiffuse_start + nSpecBands - 1
    flux_data_device%ixSpectralBelowCanopyDirect_start = flux_data_device%ixSpectralIncomingDiffuse_end + 1
    flux_data_device%ixSpectralBelowCanopyDirect_end = flux_data_device%ixSpectralBelowCanopyDirect_start + nSpecBands - 1
    flux_data_device%ixSpectralBelowCanopyDiffuse_start = flux_data_device%ixSpectralBelowCanopyDirect_end + 1
    flux_data_device%ixSpectralBelowCanopyDiffuse_end = flux_data_device%ixSpectralBelowCanopyDiffuse_start + nSpecBands - 1
    flux_data_device%ixMLayerTranspire_start = flux_data_device%ixSpectralBelowCanopyDiffuse_end + 1
    flux_data_device%ixMLayerTranspire_end = flux_data_device%ixMLayerTranspire_start + nSoil - 1
    flux_data_device%ixILayerConductiveFlux_start = flux_data_device%ixMLayerTranspire_end + 1
    flux_data_device%ixILayerConductiveFlux_end = flux_data_device%ixILayerConductiveFlux_start + nLayers
    flux_data_device%ixILayerAdvectiveFlux_start = flux_data_device%ixILayerConductiveFlux_end + 1
    flux_data_device%ixILayerAdvectiveFlux_end = flux_data_device%ixILayerAdvectiveFlux_start + nLayers
    flux_data_device%ixILayerNrgFlux_start = flux_data_device%ixILayerAdvectiveFlux_end + 1
    flux_data_device%ixILayerNrgFlux_end = flux_data_device%ixILayerNrgFlux_start + nLayers
    flux_data_device%ixMLayerNrgFlux_start = flux_data_device%ixILayerNrgFlux_end + 1
    flux_data_device%ixMLayerNrgFlux_end = flux_data_device%ixMLayerNrgFlux_start + nLayers-1
    flux_data_device%ixILayerLiqFluxSnow_start = flux_data_device%ixMLayerNrgFlux_end + 1
    flux_data_device%ixILayerLiqFluxSnow_end = flux_data_device%ixILayerLiqFluxSnow_start + nSnow
    flux_data_device%ixMLayerLiqFluxSnow_start = flux_data_device%ixILayerLiqFluxSnow_end + 1
    flux_data_device%ixMLayerLiqFluxSnow_end = flux_data_device%ixMLayerLiqFluxSnow_start + nSnow - 1
    flux_data_device%ixMLayerSatHydCondMP_start = flux_data_device%ixMLayerLiqFluxSnow_end + 1
    flux_data_device%ixMLayerSatHydCondMP_end = flux_data_device%ixMLayerSatHydCondMP_start + nSoil - 1
    flux_data_device%ixMLayerSatHydCond_start = flux_data_device%ixMLayerSatHydCondMP_end + 1
    flux_data_device%ixMLayerSatHydCond_end = flux_data_device%ixMLayerSatHydCond_start + nSoil - 1
    flux_data_device%ixILayerSatHydCond_start = flux_data_device%ixMLayerSatHydCond_end + 1
    flux_data_device%ixILayerSatHydCond_end = flux_data_device%ixILayerSatHydCond_start + nSoil
    flux_data_device%ixMLayerHydCond_start = flux_data_device%ixILayerSatHydCond_end + 1
    flux_data_device%ixMLayerHydCond_end = flux_data_device%ixMLayerHydCond_start + nSoil - 1
    flux_data_device%ixILayerLiqFluxSoil_start = flux_data_device%ixMLayerHydCond_end + 1
    flux_data_device%ixILayerLiqFluxSoil_end = flux_data_device%ixILayerLiqFluxSoil_start + nSoil
    flux_data_device%ixMLayerLiqFluxSoil_start = flux_data_device%ixILayerLiqFluxSoil_end + 1
    flux_data_device%ixMLayerLiqFluxSoil_end = flux_data_device%ixMLayerLiqFluxSoil_start + nSoil - 1
    flux_data_device%ixMLayerBaseflow_start = flux_data_device%ixMLayerLiqFluxSoil_end + 1
    flux_data_device%ixMLayerBaseflow_end = flux_data_device%ixMLayerBaseflow_start + nSoil - 1
    flux_data_device%ixMLayerColumnInflow_start = flux_data_device%ixMLayerBaseflow_end + 1
    flux_data_device%ixMLayerColumnInflow_end = flux_data_device%ixMLayerColumnInflow_start + nSoil - 1
    flux_data_device%ixMLayerColumnOutflow_start = flux_data_device%ixMLayerColumnInflow_end + 1
    flux_data_device%ixMLayerColumnOutflow_end = flux_data_device%ixMLayerColumnOutflow_start + nSoil - 1
    flux_data_device%numFluxData = flux_data_device%ixMLayerColumnOutflow_end


    allocate(flux_data_device%data(flux_data_device%numFluxData,nGRU))
    flux_data_device%data = 0._rkind

  end subroutine allocate_device_flux_prev

  subroutine finalize_device_flux_data(flux_data_device, flux_data,indxStruct)
    implicit none
    type(flux_data_device), intent(inout) :: flux_data_device
    type(gru_hru_doubleVec), intent(inout) :: flux_data
    type(gru_hru_intVec),intent(in) :: indxStruct
    integer(i4b) :: nSnow,nSoil
    integer(i4b) :: nLayers,iGRU,iLayer,nGRU

    nGRU = size(flux_data%gru)
    do iGRU=1,nGRU
      nSnow = indxStruct%gru(iGRU)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
      nSoil = indxStruct%gru(iGRU)%hru(1)%var(iLookINDEX%nSoil)%dat(1)
    nLayers = nSnow + nSoil

    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanairNetNrgFlux,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyNetNrgFlux,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1) = flux_data_device%data(flux_data_device%ixScalarGroundNetNrgFlux,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyNetLiqFlux,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarRainfall)%dat(1) = flux_data_device%data(flux_data_device%ixScalarRainfall,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSnowfall)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSnowfall,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%spectralIncomingDirect)%dat = flux_data_device%data(flux_data_device%ixSpectralIncomingDirect_start:flux_data_device%ixSpectralIncomingDirect_end,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%spectralIncomingDiffuse)%dat = flux_data_device%data(flux_data_device%ixSpectralIncomingDiffuse_start:flux_data_device%ixSpectralIncomingDiffuse_end,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopySunlitPAR)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopySunlitPAR,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyShadedPAR)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyShadedPAR,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%spectralBelowCanopyDirect)%dat = flux_data_device%data(flux_data_device%ixSpectralBelowCanopyDirect_start:flux_data_device%ixSpectralBelowCanopyDirect_end,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%spectralBelowCanopyDiffuse)%dat = flux_data_device%data(flux_data_device%ixSpectralBelowCanopyDiffuse_start:flux_data_device%ixSpectralBelowCanopyDiffuse_end,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarBelowCanopySolar)%dat(1) = flux_data_device%data(flux_data_device%ixScalarBelowCanopySolar,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyAbsorbedSolar)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyAbsorbedSolar,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1) = flux_data_device%data(flux_data_device%ixScalarGroundAbsorbedSolar,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadCanopy)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadCanopy,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadGround)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadGround,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadUbound2Canopy)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadUbound2Canopy,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadUbound2Ground)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadUbound2Ground,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadUbound2Ubound)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadUbound2Ubound,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadCanopy2Ubound)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadCanopy2Ubound,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadCanopy2Ground)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadCanopy2Ground,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadCanopy2Canopy)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadCanopy2Canopy,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadGround2Ubound)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadGround2Ubound,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadGround2Canopy)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadGround2Canopy,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWNetCanopy)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWNetCanopy,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWNetGround)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWNetGround,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWNetUbound)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWNetUbound,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarEddyDiffusCanopyTop)%dat(1) = flux_data_device%data(flux_data_device%ixScalarEddyDiffusCanopyTop,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarFrictionVelocity)%dat(1) = flux_data_device%data(flux_data_device%ixScalarFrictionVelocity,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarWindspdCanopyTop)%dat(1) = flux_data_device%data(flux_data_device%ixScalarWindspdCanopyTop,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarWindspdCanopyBottom)%dat(1) = flux_data_device%data(flux_data_device%ixScalarWindspdCanopyBottom,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundResistance)%dat(1) = flux_data_device%data(flux_data_device%ixScalarGroundResistance,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyResistance)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyResistance,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLeafResistance)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLeafResistance,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSoilResistance)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSoilResistance,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSenHeatTotal)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSenHeatTotal,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSenHeatCanopy)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSenHeatCanopy,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSenHeatGround)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSenHeatGround,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLatHeatTotal)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLatHeatTotal,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLatHeatCanopyEvap,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLatHeatCanopyTrans)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLatHeatCanopyTrans,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLatHeatGround)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLatHeatGround,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyAdvectiveHeatFlux)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyAdvectiveHeatFlux,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundAdvectiveHeatFlux)%dat(1) = flux_data_device%data(flux_data_device%ixScalarGroundAdvectiveHeatFlux,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopySublimation)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopySublimation,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSnowSublimation)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSnowSublimation,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarStomResistSunlit)%dat(1) = flux_data_device%data(flux_data_device%ixScalarStomResistSunlit,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarStomResistShaded)%dat(1) = flux_data_device%data(flux_data_device%ixScalarStomResistShaded,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarPhotosynthesisSunlit)%dat(1) = flux_data_device%data(flux_data_device%ixScalarPhotosynthesisSunlit,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarPhotosynthesisShaded)%dat(1) = flux_data_device%data(flux_data_device%ixScalarPhotosynthesisShaded,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyTranspiration)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyTranspiration,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyEvaporation)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyEvaporation,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundEvaporation)%dat(1) = flux_data_device%data(flux_data_device%ixScalarGroundEvaporation,iGRU)
if (size(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerTranspire)%dat).ne.nLayers) then
  deallocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerTranspire)%dat)
  allocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerTranspire)%dat(1:nLayers))
end if
do iLayer=1,nLayers
  flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerTranspire)%dat(iLayer) = flux_data_device%data(flux_data_device%ixMLayerTranspire_start+iLayer-1,iGRU)
end do
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarThroughfallSnow)%dat(1) = flux_data_device%data(flux_data_device%ixScalarThroughfallSnow,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarThroughfallRain)%dat(1) = flux_data_device%data(flux_data_device%ixScalarThroughfallRain,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopySnowUnloading,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyLiqDrainage,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyMeltFreeze)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyMeltFreeze,iGRU)
if (size(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerConductiveFlux)%dat).ne.nLayers+1) then
  deallocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerConductiveFlux)%dat)
  allocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerConductiveFlux)%dat(0:nLayers))
end if
if (size(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerAdvectiveFlux)%dat).ne.nLayers+1) then
  deallocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerAdvectiveFlux)%dat)
  allocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerAdvectiveFlux)%dat(0:nLayers))
end if
if (size(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerNrgFlux)%dat).ne.nLayers+1) then
  deallocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerNrgFlux)%dat)
  allocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerNrgFlux)%dat(0:nLayers))
end if
do iLayer=0,nLayers
  flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerConductiveFlux)%dat(iLayer) = flux_data_device%data(flux_data_device%ixILayerConductiveFlux_start+iLayer,iGRU)
  flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerAdvectiveFlux)%dat(iLayer) = flux_data_device%data(flux_data_device%ixILayerAdvectiveFlux_start+iLayer,iGRU)
  flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerNrgFlux)%dat(iLayer) = flux_data_device%data(flux_data_device%ixILayerNrgFlux_start+iLayer,iGRU)
end do
if (size(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerNrgFlux)%dat).ne.nLayers) then
  deallocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerNrgFlux)%dat)
  allocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerNrgFlux)%dat(1:nLayers))
end if
do iLayer=1,nLayers
  flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerNrgFlux)%dat(iLayer) = flux_data_device%data(flux_data_device%ixMLayerNrgFlux_start+iLayer-1,iGRU)
end do
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSnowDrainage)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSnowDrainage,iGRU)
if (size(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerLiqFluxSnow)%dat).ne.nSnow+1) then
  deallocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerLiqFluxSnow)%dat)
  allocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerLiqFluxSnow)%dat(0:nSnow))
end if
do iLayer=0,nSnow
  flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerLiqFluxSnow)%dat(iLayer) = flux_data_device%data(flux_data_device%ixILayerLiqFluxSnow_start + iLayer,iGRU)
end do
if (size(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerLiqFluxSnow)%dat).ne.nSnow) then
  deallocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerLiqFluxSnow)%dat)
  allocate(flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerLiqFluxSnow)%dat(1:nSnow))
end if
do iLayer=1,nSnow
  flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerLiqFluxSnow)%dat(iLayer) = flux_data_device%data(flux_data_device%ixMLayerLiqFluxSnow_start+iLayer-1,iGRU)
end do

flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarRainPlusMelt)%dat(1) = flux_data_device%data(flux_data_device%ixScalarRainPlusMelt,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarMaxInfilRate)%dat(1) = flux_data_device%data(flux_data_device%ixScalarMaxInfilRate,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarInfiltration)%dat(1) = flux_data_device%data(flux_data_device%ixScalarInfiltration,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarExfiltration)%dat(1) = flux_data_device%data(flux_data_device%ixScalarExfiltration,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSurfaceRunoff)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSurfaceRunoff,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSurfaceRunoff_IE)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSurfaceRunoff_IE,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSurfaceRunoff_SE)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSurfaceRunoff_SE,iGRU)

do iLayer=1,nSoil
  flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerSatHydCondMP)%dat(iLayer) = flux_data_device%data(flux_data_device%ixMLayerSatHydCondMP_start+iLayer-1,iGRU)
  flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerSatHydCond)%dat(iLayer) = flux_data_device%data(flux_data_device%ixMLayerSatHydCond_start+iLayer-1,iGRU)
end do
do iLayer=0,nSoil
  flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerSatHydCond)%dat(iLayer) = flux_data_device%data(flux_data_device%ixILayerSatHydCond_start+iLayer,iGRU)
end do
do iLayer=1,nSoil
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerHydCond)%dat(iLayer) = flux_data_device%data(flux_data_device%ixmLayerHydCond_start + iLayer-1,iGRU)
end do
do iLayer=0,nSoil
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerLiqFluxSoil)%dat(iLayer) = flux_data_device%data(flux_data_device%ixiLayerLiqFluxSoil_start + iLayer,iGRU)
end do
do iLayer=1,nSoil
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerLiqFluxSoil)%dat(iLayer) = flux_data_device%data(flux_data_device%ixmLayerLiqFluxSoil_start + iLayer-1,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerBaseflow)%dat(iLayer) = flux_data_device%data(flux_data_device%ixmLayerBaseflow_start + iLayer-1,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerColumnInflow)%dat(iLayer) = flux_data_device%data(flux_data_device%ixmLayerColumnInflow_start + iLayer-1,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerColumnOutflow)%dat(iLayer) = flux_data_device%data(flux_data_device%ixmLayerColumnOutflow_start + iLayer-1,iGRU)
end do
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSoilBaseflow)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSoilBaseflow,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSoilDrainage)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSoilDrainage,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarAquiferRecharge)%dat(1) = flux_data_device%data(flux_data_device%ixScalarAquiferRecharge,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarAquiferTranspire)%dat(1) = flux_data_device%data(flux_data_device%ixScalarAquiferTranspire,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarAquiferBaseflow)%dat(1) = flux_data_device%data(flux_data_device%ixScalarAquiferBaseflow,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarTotalET)%dat(1) = flux_data_device%data(flux_data_device%ixScalarTotalET,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarTotalRunoff)%dat(1) = flux_data_device%data(flux_data_device%ixScalarTotalRunoff,iGRU)
flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarNetRadiation)%dat(1) = flux_data_device%data(flux_data_device%ixScalarNetRadiation,iGRU)
end do

  end subroutine finalize_device_flux_data
  subroutine finalize_device_flux_count(flux_data_device, flux_data,nSnow,nSoil)
    implicit none
    type(flux_count_device), intent(inout) :: flux_data_device
    type(var_ilength), intent(inout) :: flux_data
    integer(i4b),intent(in) :: nSnow,nSoil
    integer(i4b) :: nLayers,iGRU,iLayer

    nLayers = nSnow + nSoil
    iGRU = 1

    flux_data%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanairNetNrgFlux,iGRU)
    flux_data%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyNetNrgFlux,iGRU)
    flux_data%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1) = flux_data_device%data(flux_data_device%ixScalarGroundNetNrgFlux,iGRU)
    flux_data%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyNetLiqFlux,iGRU)
    flux_data%var(iLookFLUX%scalarRainfall)%dat(1) = flux_data_device%data(flux_data_device%ixScalarRainfall,iGRU)
    flux_data%var(iLookFLUX%scalarSnowfall)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSnowfall,iGRU)
    flux_data%var(iLookFLUX%spectralIncomingDirect)%dat = flux_data_device%data(flux_data_device%ixSpectralIncomingDirect_start:flux_data_device%ixSpectralIncomingDirect_end,iGRU)
    flux_data%var(iLookFLUX%spectralIncomingDiffuse)%dat = flux_data_device%data(flux_data_device%ixSpectralIncomingDiffuse_start:flux_data_device%ixSpectralIncomingDiffuse_end,iGRU)
    flux_data%var(iLookFLUX%scalarCanopySunlitPAR)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopySunlitPAR,iGRU)
    flux_data%var(iLookFLUX%scalarCanopyShadedPAR)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyShadedPAR,iGRU)
    flux_data%var(iLookFLUX%spectralBelowCanopyDirect)%dat = flux_data_device%data(flux_data_device%ixSpectralBelowCanopyDirect_start:flux_data_device%ixSpectralBelowCanopyDirect_end,iGRU)
    flux_data%var(iLookFLUX%spectralBelowCanopyDiffuse)%dat = flux_data_device%data(flux_data_device%ixSpectralBelowCanopyDiffuse_start:flux_data_device%ixSpectralBelowCanopyDiffuse_end,iGRU)
    flux_data%var(iLookFLUX%scalarBelowCanopySolar)%dat(1) = flux_data_device%data(flux_data_device%ixScalarBelowCanopySolar,iGRU)
    flux_data%var(iLookFLUX%scalarCanopyAbsorbedSolar)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyAbsorbedSolar,iGRU)
    flux_data%var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1) = flux_data_device%data(flux_data_device%ixScalarGroundAbsorbedSolar,iGRU)
flux_data%var(iLookFLUX%scalarLWRadCanopy)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadCanopy,iGRU)
flux_data%var(iLookFLUX%scalarLWRadGround)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadGround,iGRU)
flux_data%var(iLookFLUX%scalarLWRadUbound2Canopy)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadUbound2Canopy,iGRU)
flux_data%var(iLookFLUX%scalarLWRadUbound2Ground)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadUbound2Ground,iGRU)
flux_data%var(iLookFLUX%scalarLWRadUbound2Ubound)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadUbound2Ubound,iGRU)
flux_data%var(iLookFLUX%scalarLWRadCanopy2Ubound)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadCanopy2Ubound,iGRU)
flux_data%var(iLookFLUX%scalarLWRadCanopy2Ground)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadCanopy2Ground,iGRU)
flux_data%var(iLookFLUX%scalarLWRadCanopy2Canopy)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadCanopy2Canopy,iGRU)
flux_data%var(iLookFLUX%scalarLWRadGround2Ubound)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadGround2Ubound,iGRU)
flux_data%var(iLookFLUX%scalarLWRadGround2Canopy)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWRadGround2Canopy,iGRU)
flux_data%var(iLookFLUX%scalarLWNetCanopy)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWNetCanopy,iGRU)
flux_data%var(iLookFLUX%scalarLWNetGround)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWNetGround,iGRU)
flux_data%var(iLookFLUX%scalarLWNetUbound)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLWNetUbound,iGRU)
flux_data%var(iLookFLUX%scalarEddyDiffusCanopyTop)%dat(1) = flux_data_device%data(flux_data_device%ixScalarEddyDiffusCanopyTop,iGRU)
flux_data%var(iLookFLUX%scalarFrictionVelocity)%dat(1) = flux_data_device%data(flux_data_device%ixScalarFrictionVelocity,iGRU)
flux_data%var(iLookFLUX%scalarWindspdCanopyTop)%dat(1) = flux_data_device%data(flux_data_device%ixScalarWindspdCanopyTop,iGRU)
flux_data%var(iLookFLUX%scalarWindspdCanopyBottom)%dat(1) = flux_data_device%data(flux_data_device%ixScalarWindspdCanopyBottom,iGRU)
flux_data%var(iLookFLUX%scalarGroundResistance)%dat(1) = flux_data_device%data(flux_data_device%ixScalarGroundResistance,iGRU)
flux_data%var(iLookFLUX%scalarCanopyResistance)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyResistance,iGRU)
flux_data%var(iLookFLUX%scalarLeafResistance)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLeafResistance,iGRU)
flux_data%var(iLookFLUX%scalarSoilResistance)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSoilResistance,iGRU)
flux_data%var(iLookFLUX%scalarSenHeatTotal)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSenHeatTotal,iGRU)
flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSenHeatCanopy,iGRU)
flux_data%var(iLookFLUX%scalarSenHeatGround)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSenHeatGround,iGRU)
flux_data%var(iLookFLUX%scalarLatHeatTotal)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLatHeatTotal,iGRU)
flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLatHeatCanopyEvap,iGRU)
flux_data%var(iLookFLUX%scalarLatHeatCanopyTrans)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLatHeatCanopyTrans,iGRU)
flux_data%var(iLookFLUX%scalarLatHeatGround)%dat(1) = flux_data_device%data(flux_data_device%ixScalarLatHeatGround,iGRU)
flux_data%var(iLookFLUX%scalarCanopyAdvectiveHeatFlux)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyAdvectiveHeatFlux,iGRU)
flux_data%var(iLookFLUX%scalarGroundAdvectiveHeatFlux)%dat(1) = flux_data_device%data(flux_data_device%ixScalarGroundAdvectiveHeatFlux,iGRU)
flux_data%var(iLookFLUX%scalarCanopySublimation)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopySublimation,iGRU)
flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSnowSublimation,iGRU)
flux_data%var(iLookFLUX%scalarStomResistSunlit)%dat(1) = flux_data_device%data(flux_data_device%ixScalarStomResistSunlit,iGRU)
flux_data%var(iLookFLUX%scalarStomResistShaded)%dat(1) = flux_data_device%data(flux_data_device%ixScalarStomResistShaded,iGRU)
flux_data%var(iLookFLUX%scalarPhotosynthesisSunlit)%dat(1) = flux_data_device%data(flux_data_device%ixScalarPhotosynthesisSunlit,iGRU)
flux_data%var(iLookFLUX%scalarPhotosynthesisShaded)%dat(1) = flux_data_device%data(flux_data_device%ixScalarPhotosynthesisShaded,iGRU)
flux_data%var(iLookFLUX%scalarCanopyTranspiration)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyTranspiration,iGRU)
flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyEvaporation,iGRU)
flux_data%var(iLookFLUX%scalarGroundEvaporation)%dat(1) = flux_data_device%data(flux_data_device%ixScalarGroundEvaporation,iGRU)
do iLayer=1,nLayers
  flux_data%var(iLookFLUX%mLayerTranspire)%dat(iLayer) = flux_data_device%data(flux_data_device%ixMLayerTranspire_start+iLayer-1,iGRU)
end do
    flux_data%var(iLookFLUX%scalarThroughfallSnow)%dat(1) = flux_data_device%data(flux_data_device%ixScalarThroughfallSnow,iGRU)
flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1) = flux_data_device%data(flux_data_device%ixScalarThroughfallRain,iGRU)
flux_data%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopySnowUnloading,iGRU)
flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyLiqDrainage,iGRU)
flux_data%var(iLookFLUX%scalarCanopyMeltFreeze)%dat(1) = flux_data_device%data(flux_data_device%ixScalarCanopyMeltFreeze,iGRU)

do iLayer=0,nLayers
  flux_data%var(iLookFLUX%iLayerConductiveFlux)%dat(iLayer) = flux_data_device%data(flux_data_device%ixILayerConductiveFlux_start+iLayer,iGRU)
  flux_data%var(iLookFLUX%iLayerAdvectiveFlux)%dat(iLayer) = flux_data_device%data(flux_data_device%ixILayerAdvectiveFlux_start+iLayer,iGRU)
  flux_data%var(iLookFLUX%iLayerNrgFlux)%dat(iLayer) = flux_data_device%data(flux_data_device%ixILayerNrgFlux_start+iLayer,iGRU)
end do
do iLayer=1,nLayers
  flux_data%var(iLookFLUX%mLayerNrgFlux)%dat(iLayer) = flux_data_device%data(flux_data_device%ixMLayerNrgFlux_start+iLayer-1,iGRU)
end do
flux_data%var(iLookFLUX%scalarSnowDrainage)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSnowDrainage,iGRU)
do iLayer=0,nSnow
  flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat(iLayer) = flux_data_device%data(flux_data_device%ixILayerLiqFluxSnow_start + iLayer,iGRU)
end do
do iLayer=1,nSnow
  flux_data%var(iLookFLUX%mLayerLiqFluxSnow)%dat(iLayer) = flux_data_device%data(flux_data_device%ixMLayerLiqFluxSnow_start+iLayer-1,iGRU)
end do

flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1) = flux_data_device%data(flux_data_device%ixScalarRainPlusMelt,iGRU)
flux_data%var(iLookFLUX%scalarMaxInfilRate)%dat(1) = flux_data_device%data(flux_data_device%ixScalarMaxInfilRate,iGRU)
flux_data%var(iLookFLUX%scalarInfiltration)%dat(1) = flux_data_device%data(flux_data_device%ixScalarInfiltration,iGRU)
flux_data%var(iLookFLUX%scalarExfiltration)%dat(1) = flux_data_device%data(flux_data_device%ixScalarExfiltration,iGRU)
flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSurfaceRunoff,iGRU)
flux_data%var(iLookFLUX%scalarSurfaceRunoff_IE)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSurfaceRunoff_IE,iGRU)
flux_data%var(iLookFLUX%scalarSurfaceRunoff_SE)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSurfaceRunoff_SE,iGRU)

do iLayer=1,nSoil
  flux_data%var(iLookFLUX%mLayerSatHydCondMP)%dat(iLayer) = flux_data_device%data(flux_data_device%ixMLayerSatHydCondMP_start+iLayer-1,iGRU)
  flux_data%var(iLookFLUX%mLayerSatHydCond)%dat(iLayer) = flux_data_device%data(flux_data_device%ixMLayerSatHydCond_start+iLayer-1,iGRU)
end do
do iLayer=0,nSoil
  flux_data%var(iLookFLUX%iLayerSatHydCond)%dat(iLayer) = flux_data_device%data(flux_data_device%ixILayerSatHydCond_start+iLayer,iGRU)
end do
do iLayer=1,nSoil
flux_data%var(iLookFLUX%mLayerHydCond)%dat(iLayer) = flux_data_device%data(flux_data_device%ixmLayerHydCond_start + iLayer-1,iGRU)
end do
do iLayer=0,nSoil
flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat(iLayer) = flux_data_device%data(flux_data_device%ixiLayerLiqFluxSoil_start + iLayer,iGRU)
end do
do iLayer=1,nSoil
flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat(iLayer) = flux_data_device%data(flux_data_device%ixmLayerLiqFluxSoil_start + iLayer-1,iGRU)
flux_data%var(iLookFLUX%mLayerBaseflow)%dat(iLayer) = flux_data_device%data(flux_data_device%ixmLayerBaseflow_start + iLayer-1,iGRU)
flux_data%var(iLookFLUX%mLayerColumnInflow)%dat(iLayer) = flux_data_device%data(flux_data_device%ixmLayerColumnInflow_start + iLayer-1,iGRU)
flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat(iLayer) = flux_data_device%data(flux_data_device%ixmLayerColumnOutflow_start + iLayer-1,iGRU)
end do
flux_data%var(iLookFLUX%scalarSoilBaseflow)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSoilBaseflow,iGRU)
flux_data%var(iLookFLUX%scalarSoilDrainage)%dat(1) = flux_data_device%data(flux_data_device%ixScalarSoilDrainage,iGRU)
flux_data%var(iLookFLUX%scalarAquiferRecharge)%dat(1) = flux_data_device%data(flux_data_device%ixScalarAquiferRecharge,iGRU)
flux_data%var(iLookFLUX%scalarAquiferTranspire)%dat(1) = flux_data_device%data(flux_data_device%ixScalarAquiferTranspire,iGRU)
flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1) = flux_data_device%data(flux_data_device%ixScalarAquiferBaseflow,iGRU)
flux_data%var(iLookFLUX%scalarTotalET)%dat(1) = flux_data_device%data(flux_data_device%ixScalarTotalET,iGRU)
flux_data%var(iLookFLUX%scalarTotalRunoff)%dat(1) = flux_data_device%data(flux_data_device%ixScalarTotalRunoff,iGRU)
flux_data%var(iLookFLUX%scalarNetRadiation)%dat(1) = flux_data_device%data(flux_data_device%ixScalarNetRadiation,iGRU)


  end subroutine finalize_device_flux_count

  subroutine deallocate_device_flux_data(flux_data_device)
    type(flux_data_device), intent(inout) :: flux_data_device
    deallocate(flux_data_device%data)

  end subroutine deallocate_device_flux_data
    
  subroutine allocate_device_deriv_data(deriv_data_device,nGRU,nSoil,nLayers,nSnow)
    type(deriv_data_device), intent(inout) :: deriv_data_device
    integer(i4b),intent(in) :: nGRU,nSoil,nLayers,nSnow
      
    integer(i4b) :: iGRU
    allocate(deriv_data_device%dAquiferTrans_dCanWat(nGRU))
    allocate(deriv_data_device%dAquiferTrans_dTCanair(nGRU))
    allocate(deriv_data_device%dAquiferTrans_dTCanopy(nGRU))
    allocate(deriv_data_device%dAquiferTrans_dTGround(nGRU))
    allocate(deriv_data_device%dBaseflow_dAquifer(nGRU))
    allocate(deriv_data_device%scalarCanopyLiqDrainageDeriv(nGRU))
    allocate(deriv_data_device%scalarThroughfallRainDeriv(nGRU))
    allocate(deriv_data_device%scalarCanopyLiqDeriv(nGRU))
    allocate(deriv_data_device%d2Theta_dTkCanopy2(nGRU))
    allocate(deriv_data_device%dTheta_dTkCanopy(nGRU))
    allocate(deriv_data_device%dFracLiqVeg_dTkCanopy(nGRU))
    allocate(deriv_data_device%dVolTot_dPsi0_m(nSoil,nGRU))
    allocate(deriv_data_device%dPsiLiq_dPsi0_m(nSoil,nGRU))
    allocate(deriv_data_device%dPsiLiq_dTemp_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerdTheta_dTk_m(nLayers,nGRU))
    allocate(deriv_data_device%dFracLiqWat_dTk_m(nLayers,nGRU))
    allocate(deriv_data_device%d2VolTot_dPsi02_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerd2Theta_dTk2_m(nSoil,nGRU))
    allocate(deriv_data_device%dCanairTemp_dEnthalpy(nGRU))
    allocate(deriv_data_device%dCanopyTemp_dEnthalpy(nGRU))
    allocate(deriv_data_device%dCanopyTemp_dCanWat(nGRU))
    allocate(deriv_data_device%dTemp_dEnthalpy_m(nLayers,nGRU))
    allocate(deriv_data_device%dTemp_dTheta_m(nLayers,nGRU))
    allocate(deriv_data_device%dTemp_dPsi0_m(nSoil,nGRU))
    allocate(deriv_data_device%dVolHtCapBulk_dTkCanopy(nGRU))
    allocate(deriv_data_device%dVolHtCapBulk_dCanWat(nGRU))
    allocate(deriv_data_device%dVolHtCapBulk_dPsi0_m(nSoil,nGRU))
    allocate(deriv_data_device%dVolHtCapBulk_dTheta_m(nLayers,nGRU))
    allocate(deriv_data_device%dVolHtCapBulk_dTk_m(nLayers,nGRU))
    allocate(deriv_data_device%dThermalC_dTempAbove_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dThermalC_dTempBelow_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dThermalC_dWatAbove_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dThermalC_dWatBelow_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dCm_dTk_m(nLayers,nGRU))
    allocate(deriv_data_device%dCm_dPsi0_m(nSoil,nGRU))
    allocate(deriv_data_device%dCm_dTkCanopy(nGRU))
    allocate(deriv_data_device%dCanLiq_dTcanopy(nGRU))
    allocate(deriv_data_device%dGroundNetFlux_dGroundTemp(nGRU))
    allocate(deriv_data_device%dGroundNetFlux_dCanairTemp(nGRU))
    allocate(deriv_data_device%dGroundNetFlux_dCanopyTemp(nGRU))
    allocate(deriv_data_device%dCanopyNetFlux_dCanWat(nGRU))
    allocate(deriv_data_device%dGroundNetFlux_dCanWat(nGRU))
    allocate(deriv_data_device%dCanopyEvaporation_dCanWat(nGRU))
    allocate(deriv_data_device%dCanopyEvaporation_dTCanair(nGRU))
    allocate(deriv_data_device%dCanopyEvaporation_dTCanopy(nGRU))
    allocate(deriv_data_device%dCanopyEvaporation_dTGround(nGRU))
    allocate(deriv_data_device%dGroundEvaporation_dCanWat(nGRU))
    allocate(deriv_data_device%dGroundEvaporation_dTCanair(nGRU))
    allocate(deriv_data_device%dGroundEvaporation_dTCanopy(nGRU))
    allocate(deriv_data_device%dGroundEvaporation_dTGround(nGRU))
    allocate(deriv_data_device%dCanopyNetFlux_dCanairTemp(nGRU))
    allocate(deriv_data_device%dCanopyNetFlux_dCanopyTemp(nGRU))
    allocate(deriv_data_device%dCanopyNetFlux_dGroundTemp(nGRU))
    allocate(deriv_data_device%dCanopyTrans_dCanWat(nGRU))
    allocate(deriv_data_device%dCanopyTrans_dTCanair(nGRU))
    allocate(deriv_data_device%dCanopyTrans_dTCanopy(nGRU))
    allocate(deriv_data_device%dCanopyTrans_dTGround(nGRU))
    allocate(deriv_data_device%dCanairNetFlux_dCanairTemp(nGRU))
    allocate(deriv_data_device%dCanairNetFlux_dCanopyTemp(nGRU))
    allocate(deriv_data_device%dCanairNetFlux_dGroundTemp(nGRU))
    allocate(deriv_data_device%dNrgFlux_dTempAbove_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dNrgFlux_dTempBelow_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dNrgFlux_dWatAbove_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dNrgFlux_dWatBelow_m(0:nLayers,nGRU))
    allocate(deriv_data_device%iLayerLiqFluxSnowDeriv_m(0:nSnow,nGRU))
    allocate(deriv_data_device%dq_dHydStateAbove_m(0:nSoil,nGRU))
    allocate(deriv_data_device%dq_dHydStateBelow_m(0:nSoil,nGRU))
    allocate(deriv_data_device%dq_dHydStateLayerSurfVec_m(0:nSoil,nGRU))
    allocate(deriv_data_device%dq_dNrgStateAbove_m(0:nSoil,nGRU))
    allocate(deriv_data_device%dq_dNrgStateBelow_m(0:nSoil,nGRU))
    allocate(deriv_data_device%dq_dNrgStateLayerSurfVec_m(0:nSoil,nGRU))
    allocate(deriv_data_device%mLayerdPsi_dTheta_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerdTheta_dPsi_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerdTrans_dTCanair_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerdTrans_dTCanopy_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerdTrans_dTGround_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerdTrans_dCanWat_m(nSoil,nGRU))
    allocate(deriv_data_device%dCompress_dPsi_m(nSoil,nGRU))
  end subroutine allocate_device_deriv_data
  subroutine zero_device_deriv_data(deriv_data_device)
    type(deriv_data_device), intent(inout) :: deriv_data_device
      
    deriv_data_device%d2Theta_dTkCanopy2 = 0._rkind
      deriv_data_device%d2VolTot_dPsi02_m = 0._rkind
      deriv_data_device%dCm_dTk_m = 0._rkind
      deriv_data_device%dCm_dPsi0_m = 0._rkind
      deriv_data_device%dCompress_dPsi_m = 0._rkind
      deriv_data_device%dFracLiqWat_dTk_m = 0._rkind
      deriv_data_device%dNrgFlux_dTempAbove_m = 0._rkind
      deriv_data_device%dNrgFlux_dTempBelow_m = 0._rkind
      deriv_data_device%dNrgFlux_dWatAbove_m = 0._rkind
      deriv_data_device%dNrgFlux_dWatBelow_m = 0._rkind
      deriv_data_device%dPsiLiq_dPsi0_m = 0._rkind
      deriv_data_device%dPsiLiq_dTemp_m = 0._rkind
      deriv_data_device%dTemp_dEnthalpy_m = 0._rkind
      deriv_data_device%dTemp_dPsi0_m = 0._rkind
      deriv_data_device%dTemp_dTheta_m = 0._rkind
      deriv_data_device%dThermalC_dTempAbove_m = 0._rkind
      deriv_data_device%dThermalC_dTempBelow_m = 0._rkind
      deriv_data_device%dThermalC_dWatAbove_m = 0._rkind
      deriv_data_device%dThermalC_dWatBelow_m = 0._rkind
      deriv_data_device%dVolHtCapBulk_dPsi0_m = 0._rkind
      deriv_data_device%dVolHtCapBulk_dTheta_m = 0._rkind
      deriv_data_device%dVolHtCapBulk_dTk_m = 0._rkind
      deriv_data_device%dVolTot_dPsi0_m = 0._rkind
      deriv_data_device%dq_dHydStateAbove_m = 0._rkind
      deriv_data_device%dq_dHydStateBelow_m = 0._rkind
      deriv_data_device%dq_dHydStateLayerSurfVec_m = 0._rkind
      deriv_data_device%dq_dNrgStateAbove_m = 0._rkind
      deriv_data_device%dq_dNrgStateBelow_m = 0._rkind
      deriv_data_device%dq_dNrgStateLayerSurfVec_m = 0._rkind
      deriv_data_device%iLayerLiqFluxSnowDeriv_m = 0._rkind
      deriv_data_device%mLayerd2Theta_dTk2_m = 0._rkind
      deriv_data_device%mLayerdPsi_dTheta_m = 0._rkind
      deriv_data_device%mLayerdTheta_dPsi_m = 0._rkind
      deriv_data_device%mLayerdTheta_dTk_m = 0._rkind
      deriv_data_device%mLayerdTrans_dCanWat_m = 0._rkind
      deriv_data_device%mLayerdTrans_dTCanair_m = 0._rkind
      deriv_data_device%mLayerdTrans_dTCanopy_m = 0._rkind
      deriv_data_device%mLayerdTrans_dTGround_m = 0._rkind
    deriv_data_device%dAquiferTrans_dCanWat = 0._rkind
    deriv_data_device%dAquiferTrans_dTCanair = 0._rkind
    deriv_data_device%dAquiferTrans_dTCanopy = 0._rkind
    deriv_data_device%dAquiferTrans_dTGround = 0._rkind
    deriv_data_device%dBaseflow_dAquifer = 0._rkind
    deriv_data_device%dCanLiq_dTcanopy = 0._rkind
    deriv_data_device%dCanairNetFlux_dCanairTemp = 0._rkind
    deriv_data_device%dCanairNetFlux_dCanopyTemp = 0._rkind
    deriv_data_device%dCanairNetFlux_dGroundTemp = 0._rkind
    deriv_data_device%dCanairTemp_dEnthalpy = 0._rkind
    deriv_data_device%dCanopyEvaporation_dCanWat = 0._rkind
    deriv_data_device%dCanopyEvaporation_dTCanair = 0._rkind
    deriv_data_device%dCanopyEvaporation_dTCanopy = 0._rkind
    deriv_data_device%dCanopyEvaporation_dTGround = 0._rkind
    deriv_data_device%dCanopyNetFlux_dCanWat = 0._rkind
    deriv_data_device%dCanopyNetFlux_dCanairTemp = 0._rkind
    deriv_data_device%dCanopyNetFlux_dCanopyTemp = 0._rkind
    deriv_data_device%dCanopyNetFlux_dGroundTemp = 0._rkind
    deriv_data_device%dCanopyTemp_dCanWat = 0._rkind
    deriv_data_device%dCanopyTemp_dEnthalpy = 0._rkind
    deriv_data_device%dCanopyTrans_dCanWat = 0._rkind
    deriv_data_device%dCanopyTrans_dTCanair = 0._rkind
    deriv_data_device%dCanopyTrans_dTCanopy = 0._rkind
    deriv_data_device%dCanopyTrans_dTGround = 0._rkind
    deriv_data_device%dCm_dTkCanopy = 0._rkind
    deriv_data_device%dFracLiqVeg_dTkCanopy = 0._rkind
    deriv_data_device%dGroundEvaporation_dCanWat = 0._rkind
    deriv_data_device%dGroundEvaporation_dTCanair = 0._rkind
    deriv_data_device%dGroundEvaporation_dTCanopy = 0._rkind
    deriv_data_device%dGroundEvaporation_dTGround = 0._rkind
    deriv_data_device%dGroundNetFlux_dCanWat = 0._rkind
    deriv_data_device%dGroundNetFlux_dCanairTemp = 0._rkind
    deriv_data_device%dGroundNetFlux_dCanopyTemp = 0._rkind
    deriv_data_device%dGroundNetFlux_dGroundTemp = 0._rkind
    deriv_data_device%dTheta_dTkCanopy = 0._rkind
    deriv_data_device%dVolHtCapBulk_dCanWat = 0._rkind
    deriv_data_device%dVolHtCapBulk_dTkCanopy = 0._rkind
    deriv_data_device%scalarCanopyLiqDeriv = 0._rkind
    deriv_data_device%scalarCanopyLiqDrainageDeriv = 0._rkind
    deriv_data_device%scalarThroughfallRainDeriv = 0._rkind
  end subroutine zero_device_deriv_data

  subroutine allocate_device_deriv_prev(deriv_data_device,nGRU,nSoil,nLayers,nSnow)
    type(deriv_data_device), intent(inout) :: deriv_data_device
    integer(i4b),intent(in) :: nGRU,nSoil,nLayers,nSnow
      
    integer(i4b) :: iGRU

    allocate(deriv_data_device%dAquiferTrans_dCanWat(nGRU))
    allocate(deriv_data_device%dAquiferTrans_dTCanair(nGRU))
    allocate(deriv_data_device%dAquiferTrans_dTCanopy(nGRU))
    allocate(deriv_data_device%dAquiferTrans_dTGround(nGRU))
    allocate(deriv_data_device%dBaseflow_dAquifer(nGRU))
    allocate(deriv_data_device%scalarCanopyLiqDrainageDeriv(nGRU))
    allocate(deriv_data_device%scalarThroughfallRainDeriv(nGRU))
    allocate(deriv_data_device%scalarCanopyLiqDeriv(nGRU))
    allocate(deriv_data_device%d2Theta_dTkCanopy2(nGRU))
    allocate(deriv_data_device%dTheta_dTkCanopy(nGRU))
    allocate(deriv_data_device%dFracLiqVeg_dTkCanopy(nGRU))
    allocate(deriv_data_device%dVolTot_dPsi0_m(nSoil,nGRU))
    allocate(deriv_data_device%dPsiLiq_dPsi0_m(nSoil,nGRU))
    allocate(deriv_data_device%dPsiLiq_dTemp_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerdTheta_dTk_m(nLayers,nGRU))
    allocate(deriv_data_device%dFracLiqWat_dTk_m(nLayers,nGRU))
    allocate(deriv_data_device%d2VolTot_dPsi02_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerd2Theta_dTk2_m(nSoil,nGRU))
    allocate(deriv_data_device%dCanairTemp_dEnthalpy(nGRU))
    allocate(deriv_data_device%dCanopyTemp_dEnthalpy(nGRU))
    allocate(deriv_data_device%dCanopyTemp_dCanWat(nGRU))
    allocate(deriv_data_device%dTemp_dEnthalpy_m(nLayers,nGRU))
    allocate(deriv_data_device%dTemp_dTheta_m(nLayers,nGRU))
    allocate(deriv_data_device%dTemp_dPsi0_m(nSoil,nGRU))
    allocate(deriv_data_device%dVolHtCapBulk_dTkCanopy(nGRU))
    allocate(deriv_data_device%dVolHtCapBulk_dCanWat(nGRU))
    allocate(deriv_data_device%dVolHtCapBulk_dPsi0_m(nSoil,nGRU))
    allocate(deriv_data_device%dVolHtCapBulk_dTheta_m(nLayers,nGRU))
    allocate(deriv_data_device%dVolHtCapBulk_dTk_m(nLayers,nGRU))
    allocate(deriv_data_device%dThermalC_dTempAbove_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dThermalC_dTempBelow_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dThermalC_dWatAbove_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dThermalC_dWatBelow_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dCm_dTk_m(nLayers,nGRU))
    allocate(deriv_data_device%dCm_dPsi0_m(nSoil,nGRU))
    allocate(deriv_data_device%dCm_dTkCanopy(nGRU))
    allocate(deriv_data_device%dCanLiq_dTcanopy(nGRU))
    allocate(deriv_data_device%dGroundNetFlux_dGroundTemp(nGRU))
    allocate(deriv_data_device%dGroundNetFlux_dCanairTemp(nGRU))
    allocate(deriv_data_device%dGroundNetFlux_dCanopyTemp(nGRU))
    allocate(deriv_data_device%dCanopyNetFlux_dCanWat(nGRU))
    allocate(deriv_data_device%dGroundNetFlux_dCanWat(nGRU))
    allocate(deriv_data_device%dCanopyEvaporation_dCanWat(nGRU))
    allocate(deriv_data_device%dCanopyEvaporation_dTCanair(nGRU))
    allocate(deriv_data_device%dCanopyEvaporation_dTCanopy(nGRU))
    allocate(deriv_data_device%dCanopyEvaporation_dTGround(nGRU))
    allocate(deriv_data_device%dGroundEvaporation_dCanWat(nGRU))
    allocate(deriv_data_device%dGroundEvaporation_dTCanair(nGRU))
    allocate(deriv_data_device%dGroundEvaporation_dTCanopy(nGRU))
    allocate(deriv_data_device%dGroundEvaporation_dTGround(nGRU))
    allocate(deriv_data_device%dCanopyNetFlux_dCanairTemp(nGRU))
    allocate(deriv_data_device%dCanopyNetFlux_dCanopyTemp(nGRU))
    allocate(deriv_data_device%dCanopyNetFlux_dGroundTemp(nGRU))
    allocate(deriv_data_device%dCanopyTrans_dCanWat(nGRU))
    allocate(deriv_data_device%dCanopyTrans_dTCanair(nGRU))
    allocate(deriv_data_device%dCanopyTrans_dTCanopy(nGRU))
    allocate(deriv_data_device%dCanopyTrans_dTGround(nGRU))
    allocate(deriv_data_device%dCanairNetFlux_dCanairTemp(nGRU))
    allocate(deriv_data_device%dCanairNetFlux_dCanopyTemp(nGRU))
    allocate(deriv_data_device%dCanairNetFlux_dGroundTemp(nGRU))
    allocate(deriv_data_device%dNrgFlux_dTempAbove_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dNrgFlux_dTempBelow_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dNrgFlux_dWatAbove_m(0:nLayers,nGRU))
    allocate(deriv_data_device%dNrgFlux_dWatBelow_m(0:nLayers,nGRU))
    allocate(deriv_data_device%iLayerLiqFluxSnowDeriv_m(0:nSnow,nGRU))
    allocate(deriv_data_device%dq_dHydStateAbove_m(0:nSoil,nGRU))
    allocate(deriv_data_device%dq_dHydStateBelow_m(0:nSoil,nGRU))
    allocate(deriv_data_device%dq_dHydStateLayerSurfVec_m(0:nSoil,nGRU))
    allocate(deriv_data_device%dq_dNrgStateAbove_m(0:nSoil,nGRU))
    allocate(deriv_data_device%dq_dNrgStateBelow_m(0:nSoil,nGRU))
    allocate(deriv_data_device%dq_dNrgStateLayerSurfVec_m(0:nSoil,nGRU))
    allocate(deriv_data_device%mLayerdPsi_dTheta_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerdTheta_dPsi_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerdTrans_dTCanair_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerdTrans_dTCanopy_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerdTrans_dTGround_m(nSoil,nGRU))
    allocate(deriv_data_device%mLayerdTrans_dCanWat_m(nSoil,nGRU))
    allocate(deriv_data_device%dCompress_dPsi_m(nSoil,nGRU))
    deriv_data_device%d2Theta_dTkCanopy2 = 0._rkind
    deriv_data_device%d2VolTot_dPsi02_m = 0._rkind
    deriv_data_device%dCm_dTk_m = 0._rkind
    deriv_data_device%dCm_dPsi0_m = 0._rkind
    deriv_data_device%dCompress_dPsi_m = 0._rkind
    deriv_data_device%dFracLiqWat_dTk_m = 0._rkind
    deriv_data_device%dNrgFlux_dTempAbove_m = 0._rkind
    deriv_data_device%dNrgFlux_dTempBelow_m = 0._rkind
    deriv_data_device%dNrgFlux_dWatAbove_m = 0._rkind
    deriv_data_device%dNrgFlux_dWatBelow_m = 0._rkind
    deriv_data_device%dPsiLiq_dPsi0_m = 0._rkind
    deriv_data_device%dPsiLiq_dTemp_m = 0._rkind
    deriv_data_device%dTemp_dEnthalpy_m = 0._rkind
    deriv_data_device%dTemp_dPsi0_m = 0._rkind
    deriv_data_device%dTemp_dTheta_m = 0._rkind
    deriv_data_device%dThermalC_dTempAbove_m = 0._rkind
    deriv_data_device%dThermalC_dTempBelow_m = 0._rkind
    deriv_data_device%dThermalC_dWatAbove_m = 0._rkind
    deriv_data_device%dThermalC_dWatBelow_m = 0._rkind
    deriv_data_device%dVolHtCapBulk_dPsi0_m = 0._rkind
    deriv_data_device%dVolHtCapBulk_dTheta_m = 0._rkind
    deriv_data_device%dVolHtCapBulk_dTk_m = 0._rkind
    deriv_data_device%dVolTot_dPsi0_m = 0._rkind
    deriv_data_device%dq_dHydStateAbove_m = 0._rkind
    deriv_data_device%dq_dHydStateBelow_m = 0._rkind
    deriv_data_device%dq_dHydStateLayerSurfVec_m = 0._rkind
    deriv_data_device%dq_dNrgStateAbove_m = 0._rkind
    deriv_data_device%dq_dNrgStateBelow_m = 0._rkind
    deriv_data_device%dq_dNrgStateLayerSurfVec_m = 0._rkind
    deriv_data_device%iLayerLiqFluxSnowDeriv_m = 0._rkind
    deriv_data_device%mLayerd2Theta_dTk2_m = 0._rkind
    deriv_data_device%mLayerdPsi_dTheta_m = 0._rkind
    deriv_data_device%mLayerdTheta_dPsi_m = 0._rkind
    deriv_data_device%mLayerdTheta_dTk_m = 0._rkind
    deriv_data_device%mLayerdTrans_dCanWat_m = 0._rkind
    deriv_data_device%mLayerdTrans_dTCanair_m = 0._rkind
    deriv_data_device%mLayerdTrans_dTCanopy_m = 0._rkind
    deriv_data_device%mLayerdTrans_dTGround_m = 0._rkind

    deriv_data_device%dAquiferTrans_dCanWat = 0._rkind
    deriv_data_device%dAquiferTrans_dTCanair = 0._rkind
    deriv_data_device%dAquiferTrans_dTCanopy = 0._rkind
    deriv_data_device%dAquiferTrans_dTGround = 0._rkind
    deriv_data_device%dBaseflow_dAquifer = 0._rkind
    deriv_data_device%dCanLiq_dTcanopy = 0._rkind
    deriv_data_device%dCanairNetFlux_dCanairTemp = 0._rkind
    deriv_data_device%dCanairNetFlux_dCanopyTemp = 0._rkind
    deriv_data_device%dCanairNetFlux_dGroundTemp = 0._rkind
    deriv_data_device%dCanairTemp_dEnthalpy = 0._rkind
    deriv_data_device%dCanopyEvaporation_dCanWat = 0._rkind
    deriv_data_device%dCanopyEvaporation_dTCanair = 0._rkind
    deriv_data_device%dCanopyEvaporation_dTCanopy = 0._rkind
    deriv_data_device%dCanopyEvaporation_dTGround = 0._rkind
    deriv_data_device%dCanopyNetFlux_dCanWat = 0._rkind
    deriv_data_device%dCanopyNetFlux_dCanairTemp = 0._rkind
    deriv_data_device%dCanopyNetFlux_dCanopyTemp = 0._rkind
    deriv_data_device%dCanopyNetFlux_dGroundTemp = 0._rkind
    deriv_data_device%dCanopyTemp_dCanWat = 0._rkind
    deriv_data_device%dCanopyTemp_dEnthalpy = 0._rkind
    deriv_data_device%dCanopyTrans_dCanWat = 0._rkind
    deriv_data_device%dCanopyTrans_dTCanair = 0._rkind
    deriv_data_device%dCanopyTrans_dTCanopy = 0._rkind
    deriv_data_device%dCanopyTrans_dTGround = 0._rkind
    deriv_data_device%dCm_dTkCanopy = 0._rkind
    deriv_data_device%dFracLiqVeg_dTkCanopy = 0._rkind
    deriv_data_device%dGroundEvaporation_dCanWat = 0._rkind
    deriv_data_device%dGroundEvaporation_dTCanair = 0._rkind
    deriv_data_device%dGroundEvaporation_dTCanopy = 0._rkind
    deriv_data_device%dGroundEvaporation_dTGround = 0._rkind
    deriv_data_device%dGroundNetFlux_dCanWat = 0._rkind
    deriv_data_device%dGroundNetFlux_dCanairTemp = 0._rkind
    deriv_data_device%dGroundNetFlux_dCanopyTemp = 0._rkind
    deriv_data_device%dGroundNetFlux_dGroundTemp = 0._rkind
    deriv_data_device%dTheta_dTkCanopy = 0._rkind
    deriv_data_device%dVolHtCapBulk_dCanWat = 0._rkind
    deriv_data_device%dVolHtCapBulk_dTkCanopy = 0._rkind
    deriv_data_device%scalarCanopyLiqDeriv = 0._rkind
    deriv_data_device%scalarCanopyLiqDrainageDeriv = 0._rkind
    deriv_data_device%scalarThroughfallRainDeriv = 0._rkind
  end subroutine allocate_device_deriv_prev

  subroutine finalize_device_deriv_data(deriv_data_device, deriv_data)
    type(deriv_data_device), intent(inout) :: deriv_data_device
    type(var_dlength), intent(inout) :: deriv_data
      
    deriv_data%var(iLookDERIV%d2Theta_dTkCanopy2)%dat(1) = deriv_data_device%d2Theta_dTkCanopy2(1)
    deriv_data%var(iLookDERIV%d2VolTot_dPsi02)%dat = deriv_data_device%d2VolTot_dPsi02_m(:,1)
    deriv_data%var(iLookDERIV%dAquiferTrans_dCanWat)%dat(1) = deriv_data_device%dAquiferTrans_dCanWat(1)
    deriv_data%var(iLookDERIV%dAquiferTrans_dTCanair)%dat(1) = deriv_data_device%dAquiferTrans_dTCanair(1)
    deriv_data%var(iLookDERIV%dAquiferTrans_dTCanopy)%dat(1) = deriv_data_device%dAquiferTrans_dTCanopy(1)
    deriv_data%var(iLookDERIV%dAquiferTrans_dTGround)%dat(1) = deriv_data_device%dAquiferTrans_dTGround(1)
    deriv_data%var(iLookDERIV%dBaseflow_dAquifer)%dat(1) = deriv_data_device%dBaseflow_dAquifer(1)
    deriv_data%var(iLookDERIV%dCanLiq_dTcanopy)%dat(1) = deriv_data_device%dCanLiq_dTcanopy(1)
    deriv_data%var(iLookDERIV%dCanairNetFlux_dCanairTemp)%dat(1) = deriv_data_device%dCanairNetFlux_dCanairTemp(1)
    deriv_data%var(iLookDERIV%dCanairNetFlux_dCanopyTemp)%dat(1) = deriv_data_device%dCanairNetFlux_dCanopyTemp(1)
    deriv_data%var(iLookDERIV%dCanairNetFlux_dGroundTemp)%dat(1) = deriv_data_device%dCanairNetFlux_dGroundTemp(1)
    deriv_data%var(iLookDERIV%dCanairTemp_dEnthalpy)%dat(1) = deriv_data_device%dCanairTemp_dEnthalpy(1)
    deriv_data%var(iLookDERIV%dCanopyEvaporation_dCanWat)%dat(1) = deriv_data_device%dCanopyEvaporation_dCanWat(1)
    deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanair)%dat(1) = deriv_data_device%dCanopyEvaporation_dTCanair(1)
    deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanopy)%dat(1) = deriv_data_device%dCanopyEvaporation_dTCanopy(1)
    deriv_data%var(iLookDERIV%dCanopyEvaporation_dTGround)%dat(1) = deriv_data_device%dCanopyEvaporation_dTGround(1)
    deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanWat)%dat(1) = deriv_data_device%dCanopyNetFlux_dCanWat(1)
    deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanairTemp)%dat(1) = deriv_data_device%dCanopyNetFlux_dCanairTemp(1)
    deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanopyTemp)%dat(1) = deriv_data_device%dCanopyNetFlux_dCanopyTemp(1)
    deriv_data%var(iLookDERIV%dCanopyNetFlux_dGroundTemp)%dat(1) = deriv_data_device%dCanopyNetFlux_dGroundTemp(1)
    deriv_data%var(iLookDERIV%dCanopyTemp_dCanWat)%dat(1) = deriv_data_device%dCanopyTemp_dCanWat(1)
    deriv_data%var(iLookDERIV%dCanopyTemp_dEnthalpy)%dat(1) = deriv_data_device%dCanopyTemp_dEnthalpy(1)
    deriv_data%var(iLookDERIV%dCanopyTrans_dCanWat)%dat(1) = deriv_data_device%dCanopyTrans_dCanWat(1)
    deriv_data%var(iLookDERIV%dCanopyTrans_dTCanair)%dat(1) = deriv_data_device%dCanopyTrans_dTCanair(1)
    deriv_data%var(iLookDERIV%dCanopyTrans_dTCanopy)%dat(1) = deriv_data_device%dCanopyTrans_dTCanopy(1)
    deriv_data%var(iLookDERIV%dCanopyTrans_dTGround)%dat(1) = deriv_data_device%dCanopyTrans_dTGround(1)
    deriv_data%var(iLookDERIV%dCm_dTk)%dat = deriv_data_device%dCm_dTk_m(:,1)
    deriv_data%var(iLookDERIV%dCm_dPsi0)%dat = deriv_data_device%dCm_dPsi0_m(:,1)
    deriv_data%var(iLookDERIV%dCm_dTkCanopy)%dat(1) = deriv_data_device%dCm_dTkCanopy(1)
    deriv_data%var(iLookDERIV%dCompress_dPsi)%dat = deriv_data_device%dCompress_dPsi_m(:,1)
    deriv_data%var(iLookDERIV%dFracLiqWat_dTk)%dat = deriv_data_device%dFracLiqWat_dTk_m(:,1)
    deriv_data%var(iLookDERIV%dFracLiqVeg_dTkCanopy)%dat(1) = deriv_data_device%dFracLiqVeg_dTkCanopy(1)
    deriv_data%var(iLookDERIV%dGroundEvaporation_dCanWat)%dat(1) = deriv_data_device%dGroundEvaporation_dCanWat(1)
    deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanair)%dat(1) = deriv_data_device%dGroundEvaporation_dTCanair(1)
    deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanopy)%dat(1) = deriv_data_device%dGroundEvaporation_dTCanopy(1)
    deriv_data%var(iLookDERIV%dGroundEvaporation_dTGround)%dat(1) = deriv_data_device%dGroundEvaporation_dTGround(1)
    deriv_data%var(iLookDERIV%dGroundNetFlux_dCanWat)%dat(1) = deriv_data_device%dGroundNetFlux_dCanWat(1)
    deriv_data%var(iLookDERIV%dGroundNetFlux_dCanairTemp)%dat(1) = deriv_data_device%dGroundNetFlux_dCanairTemp(1)
    deriv_data%var(iLookDERIV%dGroundNetFlux_dCanopyTemp)%dat(1) = deriv_data_device%dGroundNetFlux_dCanopyTemp(1)
    deriv_data%var(iLookDERIV%dGroundNetFlux_dGroundTemp)%dat(1) = deriv_data_device%dGroundNetFlux_dGroundTemp(1)
    deriv_data%var(iLookDERIV%dNrgFlux_dTempAbove)%dat = deriv_data_device%dNrgFlux_dTempAbove_m(:,1)
    deriv_data%var(iLookDERIV%dNrgFlux_dTempBelow)%dat = deriv_data_device%dNrgFlux_dTempBelow_m(:,1)
    deriv_data%var(iLookDERIV%dNrgFlux_dWatAbove)%dat = deriv_data_device%dNrgFlux_dWatAbove_m(:,1)
    deriv_data%var(iLookDERIV%dNrgFlux_dWatBelow)%dat = deriv_data_device%dNrgFlux_dWatBelow_m(:,1)
    deriv_data%var(iLookDERIV%dPsiLiq_dPsi0)%dat = deriv_data_device%dPsiLiq_dPsi0_m(:,1)
    deriv_data%var(iLookDERIV%dPsiLiq_dTemp)%dat = deriv_data_device%dPsiLiq_dTemp_m(:,1)
    deriv_data%var(iLookDERIV%dTemp_dEnthalpy)%dat = deriv_data_device%dTemp_dEnthalpy_m(:,1)
    deriv_data%var(iLookDERIV%dTemp_dPsi0)%dat = deriv_data_device%dTemp_dPsi0_m(:,1)
    deriv_data%var(iLookDERIV%dTemp_dTheta)%dat = deriv_data_device%dTemp_dTheta_m(:,1)
    deriv_data%var(iLookDERIV%dThermalC_dTempAbove)%dat = deriv_data_device%dThermalC_dTempAbove_m(:,1)
    deriv_data%var(iLookDERIV%dThermalC_dTempBelow)%dat = deriv_data_device%dThermalC_dTempBelow_m(:,1)
    deriv_data%var(iLookDERIV%dThermalC_dWatAbove)%dat = deriv_data_device%dThermalC_dWatAbove_m(:,1)
    deriv_data%var(iLookDERIV%dThermalC_dWatBelow)%dat = deriv_data_device%dThermalC_dWatBelow_m(:,1)
    deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1) = deriv_data_device%dTheta_dTkCanopy(1)
    deriv_data%var(iLookDERIV%dVolHtCapBulk_dCanWat)%dat(1) = deriv_data_device%dVolHtCapBulk_dCanWat(1)
    deriv_data%var(iLookDERIV%dVolHtCapBulk_dPsi0)%dat = deriv_data_device%dVolHtCapBulk_dPsi0_m(:,1)
    deriv_data%var(iLookDERIV%dVolHtCapBulk_dTheta)%dat = deriv_data_device%dVolHtCapBulk_dTheta_m(:,1)
    deriv_data%var(iLookDERIV%dVolHtCapBulk_dTk)%dat = deriv_data_device%dVolHtCapBulk_dTk_m(:,1)
    deriv_data%var(iLookDERIV%dVolHtCapBulk_dTkCanopy)%dat(1) = deriv_data_device%dVolHtCapBulk_dTkCanopy(1)
    deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat = deriv_data_device%dVolTot_dPsi0_m(:,1)
    deriv_data%var(iLookDERIV%dq_dHydStateAbove)%dat = deriv_data_device%dq_dHydStateAbove_m(:,1)
    deriv_data%var(iLookDERIV%dq_dHydStateBelow)%dat = deriv_data_device%dq_dHydStateBelow_m(:,1)
    deriv_data%var(iLookDERIV%dq_dHydStateLayerSurfVec)%dat = deriv_data_device%dq_dHydStateLayerSurfVec_m(:,1)
    deriv_data%var(iLookDERIV%dq_dNrgStateAbove)%dat = deriv_data_device%dq_dNrgStateAbove_m(:,1)
    deriv_data%var(iLookDERIV%dq_dNrgStateBelow)%dat = deriv_data_device%dq_dNrgStateBelow_m(:,1)
    deriv_data%var(iLookDERIV%dq_dNrgStateLayerSurfVec)%dat = deriv_data_device%dq_dNrgStateLayerSurfVec_m(:,1)
    deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv)%dat = deriv_data_device%iLayerLiqFluxSnowDeriv_m(:,1)
    deriv_data%var(iLookDERIV%mLayerd2Theta_dTk2)%dat = deriv_data_device%mLayerd2Theta_dTk2_m(:,1)
    deriv_data%var(iLookDERIV%mLayerdPsi_dTheta)%dat = deriv_data_device%mLayerdPsi_dTheta_m(:,1)
    deriv_data%var(iLookDERIV%mLayerdTheta_dPsi)%dat = deriv_data_device%mLayerdTheta_dPsi_m(:,1)
    deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat = deriv_data_device%mLayerdTheta_dTk_m(:,1)
    deriv_data%var(iLookDERIV%mLayerdTrans_dCanWat)%dat = deriv_data_device%mLayerdTrans_dCanWat_m(:,1)
    deriv_data%var(iLookDERIV%mLayerdTrans_dTCanair)%dat = deriv_data_device%mLayerdTrans_dTCanair_m(:,1)
    deriv_data%var(iLookDERIV%mLayerdTrans_dTCanopy)%dat = deriv_data_device%mLayerdTrans_dTCanopy_m(:,1)
    deriv_data%var(iLookDERIV%mLayerdTrans_dTGround)%dat = deriv_data_device%mLayerdTrans_dTGround_m(:,1)
    deriv_data%var(iLookDERIV%scalarCanopyLiqDeriv)%dat(1) = deriv_data_device%scalarCanopyLiqDeriv(1)
    deriv_data%var(iLookDERIV%scalarCanopyLiqDrainageDeriv)%dat(1) = deriv_data_device%scalarCanopyLiqDrainageDeriv(1)
    deriv_data%var(iLookDERIV%scalarThroughfallRainDeriv)%dat(1) = deriv_data_device%scalarThroughfallRainDeriv(1)
  end subroutine finalize_device_deriv_data
      
  subroutine deallocate_device_deriv_data(deriv_data_device)
    type(deriv_data_device), intent(inout) :: deriv_data_device
      
    deallocate(deriv_data_device%d2Theta_dTkCanopy2)
    deallocate(deriv_data_device%d2VolTot_dPsi02_m)
    deallocate(deriv_data_device%dAquiferTrans_dCanWat)
    deallocate(deriv_data_device%dAquiferTrans_dTCanair)
    deallocate(deriv_data_device%dAquiferTrans_dTCanopy)
    deallocate(deriv_data_device%dAquiferTrans_dTGround)
    deallocate(deriv_data_device%dBaseflow_dAquifer)
    deallocate(deriv_data_device%dCanLiq_dTcanopy)
    deallocate(deriv_data_device%dCanairNetFlux_dCanairTemp)
    deallocate(deriv_data_device%dCanairNetFlux_dCanopyTemp)
    deallocate(deriv_data_device%dCanairNetFlux_dGroundTemp)
    deallocate(deriv_data_device%dCanairTemp_dEnthalpy)
    deallocate(deriv_data_device%dCanopyEvaporation_dCanWat)
    deallocate(deriv_data_device%dCanopyEvaporation_dTCanair)
    deallocate(deriv_data_device%dCanopyEvaporation_dTCanopy)
    deallocate(deriv_data_device%dCanopyEvaporation_dTGround)
    deallocate(deriv_data_device%dCanopyNetFlux_dCanWat)
    deallocate(deriv_data_device%dCanopyNetFlux_dCanairTemp)
    deallocate(deriv_data_device%dCanopyNetFlux_dCanopyTemp)
    deallocate(deriv_data_device%dCanopyNetFlux_dGroundTemp)
    deallocate(deriv_data_device%dCanopyTemp_dCanWat)
    deallocate(deriv_data_device%dCanopyTemp_dEnthalpy)
    deallocate(deriv_data_device%dCanopyTrans_dCanWat)
    deallocate(deriv_data_device%dCanopyTrans_dTCanair)
    deallocate(deriv_data_device%dCanopyTrans_dTCanopy)
    deallocate(deriv_data_device%dCanopyTrans_dTGround)
    deallocate(deriv_data_device%dCm_dTk_m)
    deallocate(deriv_data_device%dCm_dPsi0_m)
    deallocate(deriv_data_device%dCm_dTkCanopy)
    deallocate(deriv_data_device%dCompress_dPsi_m)
    deallocate(deriv_data_device%dFracLiqWat_dTk_m)
    deallocate(deriv_data_device%dFracLiqVeg_dTkCanopy)
    deallocate(deriv_data_device%dGroundEvaporation_dCanWat)
    deallocate(deriv_data_device%dGroundEvaporation_dTCanair)
    deallocate(deriv_data_device%dGroundEvaporation_dTCanopy)
    deallocate(deriv_data_device%dGroundEvaporation_dTGround)
    deallocate(deriv_data_device%dGroundNetFlux_dCanWat)
    deallocate(deriv_data_device%dGroundNetFlux_dCanairTemp)
    deallocate(deriv_data_device%dGroundNetFlux_dCanopyTemp)
    deallocate(deriv_data_device%dGroundNetFlux_dGroundTemp)
    deallocate(deriv_data_device%dNrgFlux_dTempAbove_m)
    deallocate(deriv_data_device%dNrgFlux_dTempBelow_m)
    deallocate(deriv_data_device%dNrgFlux_dWatAbove_m)
    deallocate(deriv_data_device%dNrgFlux_dWatBelow_m)
    deallocate(deriv_data_device%dPsiLiq_dPsi0_m)
    deallocate(deriv_data_device%dPsiLiq_dTemp_m)
    deallocate(deriv_data_device%dTemp_dEnthalpy_m)
    deallocate(deriv_data_device%dTemp_dPsi0_m)
    deallocate(deriv_data_device%dTemp_dTheta_m)
    deallocate(deriv_data_device%dThermalC_dTempAbove_m)
    deallocate(deriv_data_device%dThermalC_dTempBelow_m)
    deallocate(deriv_data_device%dThermalC_dWatAbove_m)
    deallocate(deriv_data_device%dThermalC_dWatBelow_m)
    deallocate(deriv_data_device%dTheta_dTkCanopy)
    deallocate(deriv_data_device%dVolHtCapBulk_dCanWat)
    deallocate(deriv_data_device%dVolHtCapBulk_dPsi0_m)
    deallocate(deriv_data_device%dVolHtCapBulk_dTheta_m)
    deallocate(deriv_data_device%dVolHtCapBulk_dTk_m)
    deallocate(deriv_data_device%dVolHtCapBulk_dTkCanopy)
    deallocate(deriv_data_device%dVolTot_dPsi0_m)
    deallocate(deriv_data_device%dq_dHydStateAbove_m)
    deallocate(deriv_data_device%dq_dHydStateBelow_m)
    deallocate(deriv_data_device%dq_dHydStateLayerSurfVec_m)
    deallocate(deriv_data_device%dq_dNrgStateAbove_m)
    deallocate(deriv_data_device%dq_dNrgStateBelow_m)
    deallocate(deriv_data_device%dq_dNrgStateLayerSurfVec_m)
    deallocate(deriv_data_device%iLayerLiqFluxSnowDeriv_m)
    deallocate(deriv_data_device%mLayerd2Theta_dTk2_m)
    deallocate(deriv_data_device%mLayerdPsi_dTheta_m)
    deallocate(deriv_data_device%mLayerdTheta_dPsi_m)
    deallocate(deriv_data_device%mLayerdTheta_dTk_m)
    deallocate(deriv_data_device%mLayerdTrans_dCanWat_m)
    deallocate(deriv_data_device%mLayerdTrans_dTCanair_m)
    deallocate(deriv_data_device%mLayerdTrans_dTCanopy_m)
    deallocate(deriv_data_device%mLayerdTrans_dTGround_m)
    deallocate(deriv_data_device%scalarCanopyLiqDeriv)
    deallocate(deriv_data_device%scalarCanopyLiqDrainageDeriv)
    deallocate(deriv_data_device%scalarThroughfallRainDeriv)
  end subroutine deallocate_device_deriv_data

  subroutine allocate_device_diag_data(diag_data_device,diag_data, nGRU, nSoil)
    use globalData,only:maxSnowLayers
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_device
    type(gru_hru_doubleVec),intent(in) :: diag_data
    integer(i4b),intent(in) :: nGRU,nSoil

    integer(i4b) :: nLayers, nSnow
    integer(i4b) :: iGRU

    nLayers = nSoil + maxSnowLayers+1
    nSnow = maxSnowLayers+1

    allocate(diag_data_device%scalarCanopyDepth(nGRU))
    allocate(diag_data_device%scalarBulkVolHeatCapVeg(nGRU))
    allocate(diag_data_device%scalarCanopyCm(nGRU))
    allocate(diag_data_device%scalarCanopyEmissivity(nGRU))
    allocate(diag_data_device%scalarRootZoneTemp(nGRU))
    allocate(diag_data_device%scalarLAI(nGRU))
    allocate(diag_data_device%scalarSAI(nGRU))
    allocate(diag_data_device%scalarExposedLAI(nGRU))
    allocate(diag_data_device%scalarExposedSAI(nGRU))
    allocate(diag_data_device%scalarAdjMeasHeight(nGRU))
    allocate(diag_data_device%scalarCanopyIceMax(nGRU))
    allocate(diag_data_device%scalarCanopyLiqMax(nGRU))
    allocate(diag_data_device%scalarGrowingSeasonIndex(nGRU))
    allocate(diag_data_device%scalarVolHtCap_air(nGRU))
    allocate(diag_data_device%scalarVolHtCap_ice(nGRU))
    allocate(diag_data_device%scalarVolHtCap_soil(nGRU))
    allocate(diag_data_device%scalarVolHtCap_water(nGRU))
    allocate(diag_data_device%mLayerVolHtCapBulk(nLayers,nGRU))
    allocate(diag_data_device%mLayerCm(nLayers,nGRU))
    allocate(diag_data_device%scalarLambda_drysoil(nGRU))
    allocate(diag_data_device%scalarLambda_wetsoil(nGRU))
    allocate(diag_data_device%mLayerThermalC(nLayers,nGRU))
    allocate(diag_data_device%iLayerThermalC(0:nLayers,nGRU))
    allocate(diag_data_device%scalarCanopyEnthTemp(nGRU))
    allocate(diag_data_device%mLayerEnthTemp(nLayers,nGRU))
    allocate(diag_data_device%scalarTotalSoilEnthalpy(nGRU))
    allocate(diag_data_device%scalarTotalSnowEnthalpy(nGRU))
    allocate(diag_data_device%scalarVPair(nGRU))
    allocate(diag_data_device%scalarVP_CanopyAir(nGRU))
    allocate(diag_data_device%scalarTwetbulb(nGRU))
    allocate(diag_data_device%scalarSnowfallTemp(nGRU))
    allocate(diag_data_device%scalarNewSnowDensity(nGRU))
    allocate(diag_data_device%scalarO2air(nGRU))
    allocate(diag_data_device%scalarCO2air(nGRU))
    allocate(diag_data_device%windspd_x(nGRU))
    allocate(diag_data_device%windspd_y(nGRU))
    allocate(diag_data_device%scalarCosZenith(nGRU))
    allocate(diag_data_device%scalarFractionDirect(nGRU))
    allocate(diag_data_device%scalarCanopySunlitFraction(nGRU))
    allocate(diag_data_device%scalarCanopySunlitLAI(nGRU))
    allocate(diag_data_device%scalarCanopyShadedLAI(nGRU))
    allocate(diag_data_device%spectralAlbGndDirect(nSpecBands,nGRU))
    allocate(diag_data_device%spectralAlbGndDiffuse(nSpecBands,nGRU))
    allocate(diag_data_device%scalarGroundAlbedo(nGRU))
    allocate(diag_data_device%scalarLatHeatSubVapCanopy(nGRU))
    allocate(diag_data_device%scalarLatHeatSubVapGround(nGRU))
    allocate(diag_data_device%scalarSatVP_CanopyTemp(nGRU))
    allocate(diag_data_device%scalarSatVP_GroundTemp(nGRU))
    allocate(diag_data_device%scalarZ0Canopy(nGRU))
    allocate(diag_data_device%scalarWindReductionFactor(nGRU))
    allocate(diag_data_device%scalarZeroPlaneDisplacement(nGRU))
    allocate(diag_data_device%scalarRiBulkCanopy(nGRU))
    allocate(diag_data_device%scalarRiBulkGround(nGRU))
    allocate(diag_data_device%scalarCanopyStabilityCorrection(nGRU))
    allocate(diag_data_device%scalarGroundStabilityCorrection(nGRU))
    allocate(diag_data_device%scalarIntercellularCO2Sunlit(nGRU))
    allocate(diag_data_device%scalarIntercellularCO2Shaded(nGRU))
    allocate(diag_data_device%scalarTranspireLim(nGRU))
    allocate(diag_data_device%scalarTranspireLimAqfr(nGRU))
    allocate(diag_data_device%scalarFoliageNitrogenFactor(nGRU))
    allocate(diag_data_device%scalarSoilRelHumidity(nGRU))
    allocate(diag_data_device%mLayerTranspireLim(nSoil,nGRU))
    allocate(diag_data_device%mLayerRootDensity(nSoil,nGRU))
    allocate(diag_data_device%scalarAquiferRootFrac(nGRU))
    allocate(diag_data_device%scalarFracLiqVeg(nGRU))
    allocate(diag_data_device%scalarCanopyWetFraction(nGRU))
    allocate(diag_data_device%scalarSnowAge(nGRU))
    allocate(diag_data_device%scalarGroundSnowFraction(nGRU))
    allocate(diag_data_device%spectralSnowAlbedoDirect(nSpecBands,nGRU))
    allocate(diag_data_device%mLayerFracLiqSnow(nSnow,nGRU))
    allocate(diag_data_device%mLayerThetaResid(nSnow,nGRU))
    allocate(diag_data_device%mLayerPoreSpace(nSnow,nGRU))
    allocate(diag_data_device%mLayerMeltFreeze(nLayers,nGRU))
    allocate(diag_data_device%scalarInfilArea(nGRU))
    allocate(diag_data_device%scalarFrozenArea(nGRU))
    allocate(diag_data_device%scalarSoilControl(nGRU))
    allocate(diag_data_device%mLayerVolFracAir(nLayers,nGRU))
    allocate(diag_data_device%mLayerTcrit(nSoil,nGRU))
    allocate(diag_data_device%mLayerCompress(nSoil,nGRU))
    allocate(diag_data_device%scalarSoilCompress(nGRU))
    allocate(diag_data_device%mLayerMatricHeadLiq(nSoil,nGRU))
    allocate(diag_data_device%scalarTotalSoilLiq(nGRU))
    allocate(diag_data_device%scalarTotalSoilIce(nGRU))
    allocate(diag_data_device%scalarTotalSoilWat(nGRU))

    allocate(diag_data_device%scalarVGn_m(nSoil,nGRU))
    allocate(diag_data_device%scalarKappa(nGRU))
    allocate(diag_data_device%scalarVolLatHt_fus(nGRU))

    allocate(diag_data_device%balanceCasNrg(nGRU))
    allocate(diag_data_device%balanceVegNrg(nGRU))
    allocate(diag_data_device%balanceLayerNrg(nLayers,nGRU))
    allocate(diag_data_device%balanceSnowNrg(nGRU))
    allocate(diag_data_device%balanceSoilNrg(nGRU))
    allocate(diag_data_device%balanceVegMass(nGRU))
    allocate(diag_data_device%balanceLayerMass(nLayers,nGRU))
    allocate(diag_data_device%balanceSnowMass(nGRU))
    allocate(diag_data_device%balanceSoilMass(nGRU))
    allocate(diag_data_device%balanceAqMass(nGRU))

    do iGRU=1,nGRU
    nLayers = size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerEnthTemp)%dat)
    nSnow = size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerFracLiqSnow)%dat)

    diag_data_device%scalarCanopyDepth(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyDepth)%dat(1)
    diag_data_device%scalarBulkVolHeatCapVeg(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1)
    diag_data_device%scalarCanopyCm(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyCm)%dat(1)
    diag_data_device%scalarCanopyEmissivity(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyEmissivity)%dat(1)
    diag_data_device%scalarRootZoneTemp(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarRootZoneTemp)%dat(1)
    diag_data_device%scalarLAI(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLAI)%dat(1)
    diag_data_device%scalarSAI(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSAI)%dat(1)
    diag_data_device%scalarExposedLAI(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarExposedLAI)%dat(1)
    diag_data_device%scalarExposedSAI(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarExposedSAI)%dat(1)
    diag_data_device%scalarAdjMeasHeight(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarAdjMeasHeight)%dat(1)
    diag_data_device%scalarCanopyIceMax(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyIceMax)%dat(1)
    diag_data_device%scalarCanopyLiqMax(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyLiqMax)%dat(1)
    diag_data_device%scalarGrowingSeasonIndex(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGrowingSeasonIndex)%dat(1)
    diag_data_device%scalarVolHtCap_air(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolHtCap_air)%dat(1)
    diag_data_device%scalarVolHtCap_ice(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolHtCap_ice)%dat(1)
    diag_data_device%scalarVolHtCap_soil(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolHtCap_soil)%dat(1)
    diag_data_device%scalarVolHtCap_water(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolHtCap_water)%dat(1)
      diag_data_device%mLayerVolHtCapBulk(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolHtCapBulk)%dat(1:nLayers)
      diag_data_device%mLayerCm(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerCm)%dat(1:nLayers)
    diag_data_device%scalarLambda_drysoil(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLambda_drysoil)%dat(1)
    diag_data_device%scalarLambda_wetsoil(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLambda_wetsoil)%dat(1)
      diag_data_device%mLayerThermalC(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThermalC)%dat(1:nLayers)
      diag_data_device%iLayerThermalC(0:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%iLayerThermalC)%dat(0:nLayers)
    diag_data_device%scalarCanopyEnthTemp(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1)
      diag_data_device%mLayerEnthTemp(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerEnthTemp)%dat(1:nLayers)
    diag_data_device%scalarTotalSoilEnthalpy(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilEnthalpy)%dat(1)
    diag_data_device%scalarTotalSnowEnthalpy(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSnowEnthalpy)%dat(1)
    diag_data_device%scalarVPair(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVPair)%dat(1)
    diag_data_device%scalarVP_CanopyAir(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVP_CanopyAir)%dat(1)
    diag_data_device%scalarTwetbulb(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTwetbulb)%dat(1)
    diag_data_device%scalarSnowfallTemp(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSnowfallTemp)%dat(1)
    diag_data_device%scalarNewSnowDensity(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarNewSnowDensity)%dat(1)
    diag_data_device%scalarO2air(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarO2air)%dat(1)
    diag_data_device%scalarCO2air(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCO2air)%dat(1)
    diag_data_device%windspd_x(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%windspd_x)%dat(1)
    diag_data_device%windspd_y(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%windspd_y)%dat(1)
    diag_data_device%scalarCosZenith(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCosZenith)%dat(1)
    diag_data_device%scalarFractionDirect(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFractionDirect)%dat(1)
    diag_data_device%scalarCanopySunlitFraction(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopySunlitFraction)%dat(1)
    diag_data_device%scalarCanopySunlitLAI(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopySunlitLAI)%dat(1)
    diag_data_device%scalarCanopyShadedLAI(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyShadedLAI)%dat(1)
      diag_data_device%spectralAlbGndDirect(:,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%spectralAlbGndDirect)%dat
      diag_data_device%spectralAlbGndDiffuse(:,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%spectralAlbGndDiffuse)%dat
    diag_data_device%scalarGroundAlbedo(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGroundAlbedo)%dat(1)
    diag_data_device%scalarLatHeatSubVapCanopy(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLatHeatSubVapCanopy)%dat(1)
    diag_data_device%scalarLatHeatSubVapGround(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLatHeatSubVapGround)%dat(1)
    diag_data_device%scalarSatVP_CanopyTemp(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSatVP_CanopyTemp)%dat(1)
    diag_data_device%scalarSatVP_GroundTemp(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSatVP_GroundTemp)%dat(1)
    diag_data_device%scalarZ0Canopy(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarZ0Canopy)%dat(1)
    diag_data_device%scalarWindReductionFactor(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarWindReductionFactor)%dat(1)
    diag_data_device%scalarZeroPlaneDisplacement(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarZeroPlaneDisplacement)%dat(1)
    diag_data_device%scalarRiBulkCanopy(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarRiBulkCanopy)%dat(1)
    diag_data_device%scalarRiBulkGround(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarRiBulkGround)%dat(1)
    diag_data_device%scalarCanopyStabilityCorrection(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyStabilityCorrection)%dat(1)
    diag_data_device%scalarGroundStabilityCorrection(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGroundStabilityCorrection)%dat(1)

    diag_data_device%scalarIntercellularCO2Sunlit(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarIntercellularCO2Sunlit)%dat(1)
    diag_data_device%scalarIntercellularCO2Shaded(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarIntercellularCO2Shaded)%dat(1)
    diag_data_device%scalarTranspireLim(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTranspireLim)%dat(1)
    diag_data_device%scalarTranspireLimAqfr(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTranspireLimAqfr)%dat(1)
    diag_data_device%scalarFoliageNitrogenFactor(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFoliageNitrogenFactor)%dat(1)
    diag_data_device%scalarSoilRelHumidity(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSoilRelHumidity)%dat(1)
      diag_data_device%mLayerTranspireLim(:,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerTranspireLim)%dat
      diag_data_device%mLayerRootDensity(:,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerRootDensity)%dat
    diag_data_device%scalarAquiferRootFrac(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarAquiferRootFrac)%dat(1)
    diag_data_device%scalarFracLiqVeg(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFracLiqVeg)%dat(1)
    diag_data_device%scalarCanopyWetFraction(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyWetFraction)%dat(1)

    diag_data_device%scalarSnowAge(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSnowAge)%dat(1)
diag_data_device%scalarGroundSnowFraction(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGroundSnowFraction)%dat(1)
  diag_data_device%spectralSnowAlbedoDirect(:,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%spectralSnowAlbedoDirect)%dat
  if (nSnow.ne.0) diag_data_device%mLayerFracLiqSnow(1:nSnow,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerFracLiqSnow)%dat(1:nSnow)
  if(nSnow.ne.0) diag_data_device%mLayerThetaResid(1:nSnow,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThetaResid)%dat(1:nSnow)
  if(nSnow.ne.0) diag_data_device%mLayerPoreSpace(1:nSnow,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerPoreSpace)%dat(1:nSnow)
  diag_data_device%mLayerMeltFreeze(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerMeltFreeze)%dat(1:nLayers)

  
  diag_data_device%scalarInfilArea(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarInfilArea)%dat(1)
  diag_data_device%scalarFrozenArea(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFrozenArea)%dat(1)
  diag_data_device%scalarSoilControl(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSoilControl)%dat(1)
    diag_data_device%mLayerVolFracAir(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolFracAir)%dat(1:nLayers)
  diag_data_device%mLayerTcrit(:,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerTcrit)%dat
  diag_data_device%mLayerCompress(:,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerCompress)%dat
  diag_data_device%scalarSoilCompress(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSoilCompress)%dat(1)
  diag_data_device%mLayerMatricHeadLiq(:,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerMatricHeadLiq)%dat
  diag_data_device%scalarTotalSoilLiq(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilLiq)%dat(1)
  diag_data_device%scalarTotalSoilIce(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilIce)%dat(1)
  diag_data_device%scalarTotalSoilWat(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilWat)%dat(1)
      diag_data_device%scalarVGn_m(1:nSoil,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVGn_m)%dat(1:nSoil)

    diag_data_device%scalarKappa(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarKappa)%dat(1)
    diag_data_device%scalarVolLatHt_fus(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolLatHt_fus)%dat(1)

    diag_data_device%numFluxCalls = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%numFluxCalls)%dat(1)
    diag_data_device%wallClockTime = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%wallClockTime)%dat(1)
    diag_data_device%meanStepSize = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%meanStepSize)%dat(1)


    diag_data_device%balanceCasNrg(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceCasNrg)%dat(1)
    diag_data_device%balanceVegNrg(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceVegNrg)%dat(1)
      diag_data_device%balanceLayerNrg(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerNrg)%dat(1:nLayers)
    diag_data_device%balanceSnowNrg(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSnowNrg)%dat(1)
    diag_data_device%balanceSoilNrg(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSoilNrg)%dat(1)
    diag_data_device%balanceVegMass(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceVegMass)%dat(1)
      diag_data_device%balanceLayerMass(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerMass)%dat(1:nLayers)
    diag_data_device%balanceSnowMass(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSnowMass)%dat(1)
    diag_data_device%balanceSoilMass(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSoilMass)%dat(1)
    diag_data_device%balanceAqMass(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceAqMass)%dat(1)

    diag_data_device%numSteps = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%numSteps)%dat(1)
    diag_data_device%numResEvals = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%numResEvals)%dat(1)
    diag_data_device%numLinSolvSetups = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%numLinSolvSetups)%dat(1)
    diag_data_device%numErrTestFails = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%numErrTestFails)%dat(1)
    diag_data_device%kLast = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%kLast)%dat(1)
    diag_data_device%kCur = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%kCur)%dat(1)
    diag_data_device%hInitUsed = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%hInitUsed)%dat(1)
    diag_data_device%hLast = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%hLast)%dat(1)
    diag_data_device%hCur = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%hCur)%dat(1)
    diag_data_device%tCur = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%tCur)%dat(1)
  end do
  end subroutine allocate_device_diag_data
    subroutine initialize_device_diag_data(diag_data_device,diag_data, nGRU, nSoil)
    use globalData,only:maxSnowLayers
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_device
    type(var_dlength),intent(in) :: diag_data
    integer(i4b),intent(in) :: nGRU,nSoil

    integer(i4b) :: nLayers, nSnow
    integer(i4b) :: iGRU

    nLayers = nSoil + maxSnowLayers+1
    nSnow = maxSnowLayers+1


    nLayers = size(diag_data%var(iLookDIAG%mLayerEnthTemp)%dat)
    nSnow = size(diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat)

    diag_data_device%scalarCanopyDepth = diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1)
    diag_data_device%scalarBulkVolHeatCapVeg = diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1)
    diag_data_device%scalarCanopyCm = diag_data%var(iLookDIAG%scalarCanopyCm)%dat(1)
    diag_data_device%scalarCanopyEmissivity = diag_data%var(iLookDIAG%scalarCanopyEmissivity)%dat(1)
    diag_data_device%scalarRootZoneTemp = diag_data%var(iLookDIAG%scalarRootZoneTemp)%dat(1)
    diag_data_device%scalarLAI = diag_data%var(iLookDIAG%scalarLAI)%dat(1)
    diag_data_device%scalarSAI = diag_data%var(iLookDIAG%scalarSAI)%dat(1)
    diag_data_device%scalarExposedLAI = diag_data%var(iLookDIAG%scalarExposedLAI)%dat(1)
    diag_data_device%scalarExposedSAI = diag_data%var(iLookDIAG%scalarExposedSAI)%dat(1)
    diag_data_device%scalarAdjMeasHeight = diag_data%var(iLookDIAG%scalarAdjMeasHeight)%dat(1)
    diag_data_device%scalarCanopyIceMax = diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1)
    diag_data_device%scalarCanopyLiqMax = diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1)
    diag_data_device%scalarGrowingSeasonIndex = diag_data%var(iLookDIAG%scalarGrowingSeasonIndex)%dat(1)
    diag_data_device%scalarVolHtCap_air = diag_data%var(iLookDIAG%scalarVolHtCap_air)%dat(1)
    diag_data_device%scalarVolHtCap_ice = diag_data%var(iLookDIAG%scalarVolHtCap_ice)%dat(1)
    diag_data_device%scalarVolHtCap_soil = diag_data%var(iLookDIAG%scalarVolHtCap_soil)%dat(1)
    diag_data_device%scalarVolHtCap_water = diag_data%var(iLookDIAG%scalarVolHtCap_water)%dat(1)
    do iGRU=1,nGRU
      diag_data_device%mLayerVolHtCapBulk(1:nLayers,iGRU) = diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat(1:nLayers)
      diag_data_device%mLayerCm(1:nLayers,iGRU) = diag_data%var(iLookDIAG%mLayerCm)%dat(1:nLayers)
    end do
    diag_data_device%scalarLambda_drysoil = diag_data%var(iLookDIAG%scalarLambda_drysoil)%dat(1)
    diag_data_device%scalarLambda_wetsoil = diag_data%var(iLookDIAG%scalarLambda_wetsoil)%dat(1)
    do iGRU=1,nGRU
      diag_data_device%mLayerThermalC(1:nLayers,iGRU) = diag_data%var(iLookDIAG%mLayerThermalC)%dat(1:nLayers)
      diag_data_device%iLayerThermalC(0:nLayers,iGRU) = diag_data%var(iLookDIAG%iLayerThermalC)%dat(0:nLayers)
    end do
    diag_data_device%scalarCanopyEnthTemp = diag_data%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1)
    do iGRU=1,nGRU
      diag_data_device%mLayerEnthTemp(1:nLayers,iGRU) = diag_data%var(iLookDIAG%mLayerEnthTemp)%dat(1:nLayers)
    end do
    diag_data_device%scalarTotalSoilEnthalpy = diag_data%var(iLookDIAG%scalarTotalSoilEnthalpy)%dat(1)
    diag_data_device%scalarTotalSnowEnthalpy = diag_data%var(iLookDIAG%scalarTotalSnowEnthalpy)%dat(1)
    diag_data_device%scalarVPair = diag_data%var(iLookDIAG%scalarVPair)%dat(1)
    diag_data_device%scalarVP_CanopyAir = diag_data%var(iLookDIAG%scalarVP_CanopyAir)%dat(1)
    diag_data_device%scalarTwetbulb = diag_data%var(iLookDIAG%scalarTwetbulb)%dat(1)
    diag_data_device%scalarSnowfallTemp = diag_data%var(iLookDIAG%scalarSnowfallTemp)%dat(1)
    diag_data_device%scalarNewSnowDensity = diag_data%var(iLookDIAG%scalarNewSnowDensity)%dat(1)
    diag_data_device%scalarO2air = diag_data%var(iLookDIAG%scalarO2air)%dat(1)
    diag_data_device%scalarCO2air = diag_data%var(iLookDIAG%scalarCO2air)%dat(1)
    diag_data_device%windspd_x = diag_data%var(iLookDIAG%windspd_x)%dat(1)
    diag_data_device%windspd_y = diag_data%var(iLookDIAG%windspd_y)%dat(1)
    diag_data_device%scalarCosZenith = diag_data%var(iLookDIAG%scalarCosZenith)%dat(1)
    diag_data_device%scalarFractionDirect = diag_data%var(iLookDIAG%scalarFractionDirect)%dat(1)
    diag_data_device%scalarCanopySunlitFraction = diag_data%var(iLookDIAG%scalarCanopySunlitFraction)%dat(1)
    diag_data_device%scalarCanopySunlitLAI = diag_data%var(iLookDIAG%scalarCanopySunlitLAI)%dat(1)
    diag_data_device%scalarCanopyShadedLAI = diag_data%var(iLookDIAG%scalarCanopyShadedLAI)%dat(1)
    do iGRU=1,nGRU
      diag_data_device%spectralAlbGndDirect(:,iGRU) = diag_data%var(iLookDIAG%spectralAlbGndDirect)%dat
      diag_data_device%spectralAlbGndDiffuse(:,iGRU) = diag_data%var(iLookDIAG%spectralAlbGndDiffuse)%dat
    end do
    diag_data_device%scalarGroundAlbedo = diag_data%var(iLookDIAG%scalarGroundAlbedo)%dat(1)
    diag_data_device%scalarLatHeatSubVapCanopy = diag_data%var(iLookDIAG%scalarLatHeatSubVapCanopy)%dat(1)
    diag_data_device%scalarLatHeatSubVapGround = diag_data%var(iLookDIAG%scalarLatHeatSubVapGround)%dat(1)
    diag_data_device%scalarSatVP_CanopyTemp = diag_data%var(iLookDIAG%scalarSatVP_CanopyTemp)%dat(1)
    diag_data_device%scalarSatVP_GroundTemp = diag_data%var(iLookDIAG%scalarSatVP_GroundTemp)%dat(1)
    diag_data_device%scalarZ0Canopy = diag_data%var(iLookDIAG%scalarZ0Canopy)%dat(1)
    diag_data_device%scalarWindReductionFactor = diag_data%var(iLookDIAG%scalarWindReductionFactor)%dat(1)
    diag_data_device%scalarZeroPlaneDisplacement = diag_data%var(iLookDIAG%scalarZeroPlaneDisplacement)%dat(1)
    diag_data_device%scalarRiBulkCanopy = diag_data%var(iLookDIAG%scalarRiBulkCanopy)%dat(1)
    diag_data_device%scalarRiBulkGround = diag_data%var(iLookDIAG%scalarRiBulkGround)%dat(1)
    diag_data_device%scalarCanopyStabilityCorrection = diag_data%var(iLookDIAG%scalarCanopyStabilityCorrection)%dat(1)
    diag_data_device%scalarGroundStabilityCorrection = diag_data%var(iLookDIAG%scalarGroundStabilityCorrection)%dat(1)

    diag_data_device%scalarIntercellularCO2Sunlit = diag_data%var(iLookDIAG%scalarIntercellularCO2Sunlit)%dat(1)
    diag_data_device%scalarIntercellularCO2Shaded = diag_data%var(iLookDIAG%scalarIntercellularCO2Shaded)%dat(1)
    diag_data_device%scalarTranspireLim = diag_data%var(iLookDIAG%scalarTranspireLim)%dat(1)
    diag_data_device%scalarTranspireLimAqfr = diag_data%var(iLookDIAG%scalarTranspireLimAqfr)%dat(1)
    diag_data_device%scalarFoliageNitrogenFactor = diag_data%var(iLookDIAG%scalarFoliageNitrogenFactor)%dat(1)
    diag_data_device%scalarSoilRelHumidity = diag_data%var(iLookDIAG%scalarSoilRelHumidity)%dat(1)
    do iGRU=1,nGRU
      diag_data_device%mLayerTranspireLim(:,iGRU) = diag_data%var(iLookDIAG%mLayerTranspireLim)%dat
      diag_data_device%mLayerRootDensity(:,iGRU) = diag_data%var(iLookDIAG%mLayerRootDensity)%dat
    end do
    diag_data_device%scalarAquiferRootFrac = diag_data%var(iLookDIAG%scalarAquiferRootFrac)%dat(1)
    diag_data_device%scalarFracLiqVeg = diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)
    diag_data_device%scalarCanopyWetFraction = diag_data%var(iLookDIAG%scalarCanopyWetFraction)%dat(1)

    diag_data_device%scalarSnowAge = diag_data%var(iLookDIAG%scalarSnowAge)%dat(1)
diag_data_device%scalarGroundSnowFraction = diag_data%var(iLookDIAG%scalarGroundSnowFraction)%dat(1)
do iGRU=1,nGRU
  diag_data_device%spectralSnowAlbedoDirect(:,iGRU) = diag_data%var(iLookDIAG%spectralSnowAlbedoDirect)%dat
  if (nSnow.ne.0) diag_data_device%mLayerFracLiqSnow(1:nSnow,iGRU) = diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat(1:nSnow)
  if(nSnow.ne.0) diag_data_device%mLayerThetaResid(1:nSnow,iGRU) = diag_data%var(iLookDIAG%mLayerThetaResid)%dat(1:nSnow)
  if(nSnow.ne.0) diag_data_device%mLayerPoreSpace(1:nSnow,iGRU) = diag_data%var(iLookDIAG%mLayerPoreSpace)%dat(1:nSnow)
  diag_data_device%mLayerMeltFreeze(1:nLayers,iGRU) = diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat(1:nLayers)

  
end do
  diag_data_device%scalarInfilArea = diag_data%var(iLookDIAG%scalarInfilArea)%dat(1)
  diag_data_device%scalarFrozenArea = diag_data%var(iLookDIAG%scalarFrozenArea)%dat(1)
  diag_data_device%scalarSoilControl = diag_data%var(iLookDIAG%scalarSoilControl)%dat(1)
  do iGRU=1,nGRU
    diag_data_device%mLayerVolFracAir(1:nLayers,iGRU) = diag_data%var(iLookDIAG%mLayerVolFracAir)%dat(1:nLayers)
  diag_data_device%mLayerTcrit(:,iGRU) = diag_data%var(iLookDIAG%mLayerTcrit)%dat
  diag_data_device%mLayerCompress(:,iGRU) = diag_data%var(iLookDIAG%mLayerCompress)%dat
  end do
  diag_data_device%scalarSoilCompress = diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1)
  do iGRU=1,nGRU
  diag_data_device%mLayerMatricHeadLiq(:,iGRU) = diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat
  end do
  diag_data_device%scalarTotalSoilLiq = diag_data%var(iLookDIAG%scalarTotalSoilLiq)%dat(1)
  diag_data_device%scalarTotalSoilIce = diag_data%var(iLookDIAG%scalarTotalSoilIce)%dat(1)
  diag_data_device%scalarTotalSoilWat = diag_data%var(iLookDIAG%scalarTotalSoilWat)%dat(1)
    do iGRU=1,nGRU
      diag_data_device%scalarVGn_m(1:nSoil,iGRU) = diag_data%var(iLookDIAG%scalarVGn_m)%dat(1:nSoil)

    end do

    diag_data_device%scalarKappa = diag_data%var(iLookDIAG%scalarKappa)%dat(1)
    diag_data_device%scalarVolLatHt_fus = diag_data%var(iLookDIAG%scalarVolLatHt_fus)%dat(1)

    diag_data_device%numFluxCalls = diag_data%var(iLookDIAG%numFluxCalls)%dat(1)
    diag_data_device%wallClockTime = diag_data%var(iLookDIAG%wallClockTime)%dat(1)
    diag_data_device%meanStepSize = diag_data%var(iLookDIAG%meanStepSize)%dat(1)


    diag_data_device%balanceCasNrg = diag_data%var(iLookDIAG%balanceCasNrg)%dat(1)
    diag_data_device%balanceVegNrg = diag_data%var(iLookDIAG%balanceVegNrg)%dat(1)
    do iGRU=1,nGRU
      diag_data_device%balanceLayerNrg(1:nLayers,iGRU) = diag_data%var(iLookDIAG%balanceLayerNrg)%dat(1:nLayers)
    end do
    diag_data_device%balanceSnowNrg = diag_data%var(iLookDIAG%balanceSnowNrg)%dat(1)
    diag_data_device%balanceSoilNrg = diag_data%var(iLookDIAG%balanceSoilNrg)%dat(1)
    diag_data_device%balanceVegMass = diag_data%var(iLookDIAG%balanceVegMass)%dat(1)
    do iGRU=1,nGRU
      diag_data_device%balanceLayerMass(1:nLayers,iGRU) = diag_data%var(iLookDIAG%balanceLayerMass)%dat(1:nLayers)
    end do
    diag_data_device%balanceSnowMass = diag_data%var(iLookDIAG%balanceSnowMass)%dat(1)
    diag_data_device%balanceSoilMass = diag_data%var(iLookDIAG%balanceSoilMass)%dat(1)
    diag_data_device%balanceAqMass = diag_data%var(iLookDIAG%balanceAqMass)%dat(1)

    diag_data_device%numSteps = diag_data%var(iLookDIAG%numSteps)%dat(1)
    diag_data_device%numResEvals = diag_data%var(iLookDIAG%numResEvals)%dat(1)
    diag_data_device%numLinSolvSetups = diag_data%var(iLookDIAG%numLinSolvSetups)%dat(1)
    diag_data_device%numErrTestFails = diag_data%var(iLookDIAG%numErrTestFails)%dat(1)
    diag_data_device%kLast = diag_data%var(iLookDIAG%kLast)%dat(1)
    diag_data_device%kCur = diag_data%var(iLookDIAG%kCur)%dat(1)
    diag_data_device%hInitUsed = diag_data%var(iLookDIAG%hInitUsed)%dat(1)
    diag_data_device%hLast = diag_data%var(iLookDIAG%hLast)%dat(1)
    diag_data_device%hCur = diag_data%var(iLookDIAG%hCur)%dat(1)
    diag_data_device%tCur = diag_data%var(iLookDIAG%tCur)%dat(1)
  end subroutine initialize_device_diag_data

      subroutine allocate_device_diag_temp(diag_data_device, nGRU, nSoil)
    use globalData,only:maxSnowLayers
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_device
    integer(i4b),intent(in) :: nGRU,nSoil

    integer(i4b) :: nLayers, nSnow
    integer(i4b) :: iGRU

    nLayers = nSoil + maxSnowLayers+1
    nSnow = maxSnowLayers+1

    allocate(diag_data_device%scalarCanopyDepth(nGRU))
    allocate(diag_data_device%scalarBulkVolHeatCapVeg(nGRU))
    allocate(diag_data_device%scalarCanopyCm(nGRU))
    allocate(diag_data_device%scalarCanopyEmissivity(nGRU))
    allocate(diag_data_device%scalarRootZoneTemp(nGRU))
    allocate(diag_data_device%scalarLAI(nGRU))
    allocate(diag_data_device%scalarSAI(nGRU))
    allocate(diag_data_device%scalarExposedLAI(nGRU))
    allocate(diag_data_device%scalarExposedSAI(nGRU))
    allocate(diag_data_device%scalarAdjMeasHeight(nGRU))
    allocate(diag_data_device%scalarCanopyIceMax(nGRU))
    allocate(diag_data_device%scalarCanopyLiqMax(nGRU))
    allocate(diag_data_device%scalarGrowingSeasonIndex(nGRU))
    allocate(diag_data_device%scalarVolHtCap_air(nGRU))
    allocate(diag_data_device%scalarVolHtCap_ice(nGRU))
    allocate(diag_data_device%scalarVolHtCap_soil(nGRU))
    allocate(diag_data_device%scalarVolHtCap_water(nGRU))
    allocate(diag_data_device%mLayerVolHtCapBulk(nLayers,nGRU))
    allocate(diag_data_device%mLayerCm(nLayers,nGRU))
    allocate(diag_data_device%scalarLambda_drysoil(nGRU))
    allocate(diag_data_device%scalarLambda_wetsoil(nGRU))
    allocate(diag_data_device%mLayerThermalC(nLayers,nGRU))
    allocate(diag_data_device%iLayerThermalC(0:nLayers,nGRU))
    allocate(diag_data_device%scalarCanopyEnthTemp(nGRU))
    allocate(diag_data_device%mLayerEnthTemp(nLayers,nGRU))
    allocate(diag_data_device%scalarTotalSoilEnthalpy(nGRU))
    allocate(diag_data_device%scalarTotalSnowEnthalpy(nGRU))
    allocate(diag_data_device%scalarVPair(nGRU))
    allocate(diag_data_device%scalarVP_CanopyAir(nGRU))
    allocate(diag_data_device%scalarTwetbulb(nGRU))
    allocate(diag_data_device%scalarSnowfallTemp(nGRU))
    allocate(diag_data_device%scalarNewSnowDensity(nGRU))
    allocate(diag_data_device%scalarO2air(nGRU))
    allocate(diag_data_device%scalarCO2air(nGRU))
    allocate(diag_data_device%windspd_x(nGRU))
    allocate(diag_data_device%windspd_y(nGRU))
    allocate(diag_data_device%scalarCosZenith(nGRU))
    allocate(diag_data_device%scalarFractionDirect(nGRU))
    allocate(diag_data_device%scalarCanopySunlitFraction(nGRU))
    allocate(diag_data_device%scalarCanopySunlitLAI(nGRU))
    allocate(diag_data_device%scalarCanopyShadedLAI(nGRU))
    allocate(diag_data_device%spectralAlbGndDirect(nSpecBands,nGRU))
    allocate(diag_data_device%spectralAlbGndDiffuse(nSpecBands,nGRU))
    allocate(diag_data_device%scalarGroundAlbedo(nGRU))
    allocate(diag_data_device%scalarLatHeatSubVapCanopy(nGRU))
    allocate(diag_data_device%scalarLatHeatSubVapGround(nGRU))
    allocate(diag_data_device%scalarSatVP_CanopyTemp(nGRU))
    allocate(diag_data_device%scalarSatVP_GroundTemp(nGRU))
    allocate(diag_data_device%scalarZ0Canopy(nGRU))
    allocate(diag_data_device%scalarWindReductionFactor(nGRU))
    allocate(diag_data_device%scalarZeroPlaneDisplacement(nGRU))
    allocate(diag_data_device%scalarRiBulkCanopy(nGRU))
    allocate(diag_data_device%scalarRiBulkGround(nGRU))
    allocate(diag_data_device%scalarCanopyStabilityCorrection(nGRU))
    allocate(diag_data_device%scalarGroundStabilityCorrection(nGRU))
    allocate(diag_data_device%scalarIntercellularCO2Sunlit(nGRU))
    allocate(diag_data_device%scalarIntercellularCO2Shaded(nGRU))
    allocate(diag_data_device%scalarTranspireLim(nGRU))
    allocate(diag_data_device%scalarTranspireLimAqfr(nGRU))
    allocate(diag_data_device%scalarFoliageNitrogenFactor(nGRU))
    allocate(diag_data_device%scalarSoilRelHumidity(nGRU))
    allocate(diag_data_device%mLayerTranspireLim(nSoil,nGRU))
    allocate(diag_data_device%mLayerRootDensity(nSoil,nGRU))
    allocate(diag_data_device%scalarAquiferRootFrac(nGRU))
    allocate(diag_data_device%scalarFracLiqVeg(nGRU))
    allocate(diag_data_device%scalarCanopyWetFraction(nGRU))
    allocate(diag_data_device%scalarSnowAge(nGRU))
    allocate(diag_data_device%scalarGroundSnowFraction(nGRU))
    allocate(diag_data_device%spectralSnowAlbedoDirect(nSpecBands,nGRU))
    allocate(diag_data_device%mLayerFracLiqSnow(nSnow,nGRU))
    allocate(diag_data_device%mLayerThetaResid(nSnow,nGRU))
    allocate(diag_data_device%mLayerPoreSpace(nSnow,nGRU))
    allocate(diag_data_device%mLayerMeltFreeze(nLayers,nGRU))
    allocate(diag_data_device%scalarInfilArea(nGRU))
    allocate(diag_data_device%scalarFrozenArea(nGRU))
    allocate(diag_data_device%scalarSoilControl(nGRU))
    allocate(diag_data_device%mLayerVolFracAir(nLayers,nGRU))
    allocate(diag_data_device%mLayerTcrit(nSoil,nGRU))
    allocate(diag_data_device%mLayerCompress(nSoil,nGRU))
    allocate(diag_data_device%scalarSoilCompress(nGRU))
    allocate(diag_data_device%mLayerMatricHeadLiq(nSoil,nGRU))
    allocate(diag_data_device%scalarTotalSoilLiq(nGRU))
    allocate(diag_data_device%scalarTotalSoilIce(nGRU))
    allocate(diag_data_device%scalarTotalSoilWat(nGRU))

    allocate(diag_data_device%scalarVGn_m(nSoil,nGRU))
    allocate(diag_data_device%scalarKappa(nGRU))
    allocate(diag_data_device%scalarVolLatHt_fus(nGRU))

    allocate(diag_data_device%balanceCasNrg(nGRU))
    allocate(diag_data_device%balanceVegNrg(nGRU))
    allocate(diag_data_device%balanceLayerNrg(nLayers,nGRU))
    allocate(diag_data_device%balanceSnowNrg(nGRU))
    allocate(diag_data_device%balanceSoilNrg(nGRU))
    allocate(diag_data_device%balanceVegMass(nGRU))
    allocate(diag_data_device%balanceLayerMass(nLayers,nGRU))
    allocate(diag_data_device%balanceSnowMass(nGRU))
    allocate(diag_data_device%balanceSoilMass(nGRU))
    allocate(diag_data_device%balanceAqMass(nGRU))

  end subroutine allocate_device_diag_temp
  subroutine finalize_device_diag_data(diag_data_device,diag_data,indxStruct)
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_device
    type(gru_hru_doubleVec),intent(inout) :: diag_data
    type(gru_hru_intVec),intent(in) :: indxStruct

    integer(i4b) :: nSnow,nSoil,nLayers,iGRU,nGRU
    nGRU = size(diag_data%gru)
    do iGRU=1,nGRU
      nSnow=indxStruct%gru(iGRU)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
      nSoil = indxStruct%gru(iGRU)%hru(1)%var(iLookINDEX%nSoil)%dat(1)
      nLayers=indxStruct%gru(iGRU)%hru(1)%var(iLookINDEX%nLayers)%dat(1)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyDepth)%dat(1) = diag_data_device%scalarCanopyDepth(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1) = diag_data_device%scalarBulkVolHeatCapVeg(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyCm)%dat(1) = diag_data_device%scalarCanopyCm(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyEmissivity)%dat(1) = diag_data_device%scalarCanopyEmissivity(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarRootZoneTemp)%dat(1) = diag_data_device%scalarRootZoneTemp(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLAI)%dat(1) = diag_data_device%scalarLAI(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSAI)%dat(1) = diag_data_device%scalarSAI(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarExposedLAI)%dat(1) = diag_data_device%scalarExposedLAI(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarExposedSAI)%dat(1) = diag_data_device%scalarExposedSAI(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarAdjMeasHeight)%dat(1) = diag_data_device%scalarAdjMeasHeight(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyIceMax)%dat(1) = diag_data_device%scalarCanopyIceMax(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyLiqMax)%dat(1) = diag_data_device%scalarCanopyLiqMax(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGrowingSeasonIndex)%dat(1) = diag_data_device%scalarGrowingSeasonIndex(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolHtCap_air)%dat(1) = diag_data_device%scalarVolHtCap_air(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolHtCap_ice)%dat(1) = diag_data_device%scalarVolHtCap_ice(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolHtCap_soil)%dat(1) = diag_data_device%scalarVolHtCap_soil(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolHtCap_water)%dat(1) = diag_data_device%scalarVolHtCap_water(iGRU)
    if (size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolHtCapBulk)%dat).ne.nLayers) then
      deallocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolHtCapBulk)%dat)
      allocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolHtCapBulk)%dat(1:nLayers))
    end if
if (size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerCm)%dat).ne.nLayers) then
  deallocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerCm)%dat)
  allocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerCm)%dat(1:nLayers))
end if
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolHtCapBulk)%dat(1:nLayers) = diag_data_device%mLayerVolHtCapBulk(1:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerCm)%dat(1:nLayers) = diag_data_device%mLayerCm(1:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLambda_drysoil)%dat(1) = diag_data_device%scalarLambda_drysoil(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLambda_wetsoil)%dat(1) = diag_data_device%scalarLambda_wetsoil(iGRU)
    if (size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThermalC)%dat).ne.nLayers) then
      deallocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThermalC)%dat)
      allocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThermalC)%dat(1:nLayers))
    end if
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThermalC)%dat(1:nLayers) = diag_data_device%mLayerThermalC(1:nLayers,iGRU)
    if (size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%iLayerThermalC)%dat).ne.nLayers+1) then
      deallocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%iLayerThermalC)%dat)
      allocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%iLayerThermalC)%dat(0:nLayers))
    end if
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%iLayerThermalC)%dat(0:nLayers) = diag_data_device%iLayerThermalC(0:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1) = diag_data_device%scalarCanopyEnthTemp(iGRU)
    if (size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerEnthTemp)%dat).ne.nLayers) then
      deallocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerEnthTemp)%dat)
      allocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerEnthTemp)%dat(1:nLayers))
    end if
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerEnthTemp)%dat(1:nLayers) = diag_data_device%mLayerEnthTemp(1:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilEnthalpy)%dat(1) = diag_data_device%scalarTotalSoilEnthalpy(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSnowEnthalpy)%dat(1) = diag_data_device%scalarTotalSnowEnthalpy(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVPair)%dat(1) = diag_data_device%scalarVPair(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVP_CanopyAir)%dat(1) = diag_data_device%scalarVP_CanopyAir(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTwetbulb)%dat(1) = diag_data_device%scalarTwetbulb(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSnowfallTemp)%dat(1) = diag_data_device%scalarSnowfallTemp(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarNewSnowDensity)%dat(1) = diag_data_device%scalarNewSnowDensity(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarO2air)%dat(1) = diag_data_device%scalarO2air(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCO2air)%dat(1) = diag_data_device%scalarCO2air(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%windspd_x)%dat(1) = diag_data_device%windspd_x(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%windspd_y)%dat(1) = diag_data_device%windspd_y(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCosZenith)%dat(1) = diag_data_device%scalarCosZenith(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFractionDirect)%dat(1) = diag_data_device%scalarFractionDirect(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopySunlitFraction)%dat(1) = diag_data_device%scalarCanopySunlitFraction(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopySunlitLAI)%dat(1) = diag_data_device%scalarCanopySunlitLAI(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyShadedLAI)%dat(1) = diag_data_device%scalarCanopyShadedLAI(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%spectralAlbGndDirect)%dat = diag_data_device%spectralAlbGndDirect(:,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%spectralAlbGndDiffuse)%dat = diag_data_device%spectralAlbGndDiffuse(:,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGroundAlbedo)%dat(1) = diag_data_device%scalarGroundAlbedo(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLatHeatSubVapCanopy)%dat(1) = diag_data_device%scalarLatHeatSubVapCanopy(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLatHeatSubVapGround)%dat(1) = diag_data_device%scalarLatHeatSubVapGround(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSatVP_CanopyTemp)%dat(1) = diag_data_device%scalarSatVP_CanopyTemp(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSatVP_GroundTemp)%dat(1) = diag_data_device%scalarSatVP_GroundTemp(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarZ0Canopy)%dat(1) = diag_data_device%scalarZ0Canopy(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarWindReductionFactor)%dat(1) = diag_data_device%scalarWindReductionFactor(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarZeroPlaneDisplacement)%dat(1) = diag_data_device%scalarZeroPlaneDisplacement(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarRiBulkCanopy)%dat(1) = diag_data_device%scalarRiBulkCanopy(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarRiBulkGround)%dat(1) = diag_data_device%scalarRiBulkGround(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyStabilityCorrection)%dat(1) = diag_data_device%scalarCanopyStabilityCorrection(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGroundStabilityCorrection)%dat(1) = diag_data_device%scalarGroundStabilityCorrection(iGRU)

    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarIntercellularCO2Sunlit)%dat(1) = diag_data_device%scalarIntercellularCO2Sunlit(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarIntercellularCO2Shaded)%dat(1) = diag_data_device%scalarIntercellularCO2Shaded(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTranspireLim)%dat(1) = diag_data_device%scalarTranspireLim(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTranspireLimAqfr)%dat(1) = diag_data_device%scalarTranspireLimAqfr(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFoliageNitrogenFactor)%dat(1) = diag_data_device%scalarFoliageNitrogenFactor(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSoilRelHumidity)%dat(1) = diag_data_device%scalarSoilRelHumidity(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerTranspireLim)%dat = diag_data_device%mLayerTranspireLim(:,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerRootDensity)%dat = diag_data_device%mLayerRootDensity(:,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarAquiferRootFrac)%dat(1) = diag_data_device%scalarAquiferRootFrac(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFracLiqVeg)%dat(1) = diag_data_device%scalarFracLiqVeg(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyWetFraction)%dat(1) = diag_data_device%scalarCanopyWetFraction(iGRU)

    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSnowAge)%dat(1) = diag_data_device%scalarSnowAge(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGroundSnowFraction)%dat(1) = diag_data_device%scalarGroundSnowFraction(iGRU)
  diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%spectralSnowAlbedoDirect)%dat = diag_data_device%spectralSnowAlbedoDirect(:,iGRU)
  if (size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerFracLiqSnow)%dat).ne.nSnow) then
    deallocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerFracLiqSnow)%dat)
    allocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerFracLiqSnow)%dat(1:nSnow))
  end if
if (size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThetaResid)%dat).ne.nSnow) then
  deallocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThetaResid)%dat)
  allocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThetaResid)%dat(1:nSnow))
end if
if (size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerPoreSpace)%dat).ne.nSnow) then
  deallocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerPoreSpace)%dat)
  allocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerPoreSpace)%dat(1:nSnow))
end if
  if (nSnow.ne.0) diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerFracLiqSnow)%dat(1:nSnow) = diag_data_device%mLayerFracLiqSnow(1:nSnow,iGRU)
  if(nSnow.ne.0) diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThetaResid)%dat(1:nSnow) = diag_data_device%mLayerThetaResid(1:nSnow,iGRU)
  if(nSnow.ne.0) diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerPoreSpace)%dat(1:nSnow) = diag_data_device%mLayerPoreSpace(1:nSnow,iGRU)
  if (size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerMeltFreeze)%dat).ne.nLayers) then
    deallocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerMeltFreeze)%dat)
    allocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerMeltFreeze)%dat(1:nLayers))
  end if
  diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerMeltFreeze)%dat(1:nLayers) = diag_data_device%mLayerMeltFreeze(1:nLayers,iGRU)


  diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarInfilArea)%dat(1) = diag_data_device%scalarInfilArea(iGRU)
  diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFrozenArea)%dat(1) = diag_data_device%scalarFrozenArea(iGRU)
  diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSoilControl)%dat(1) = diag_data_device%scalarSoilControl(iGRU)
  if (size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolFracAir)%dat).ne.nLayers) then
    deallocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolFracAir)%dat)
    allocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolFracAir)%dat(1:nLayers))
  end if
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolFracAir)%dat(1:nLayers) = diag_data_device%mLayerVolFracAir(1:nLayers,iGRU)
  diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerTcrit)%dat = diag_data_device%mLayerTcrit(:,iGRU)
  diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerCompress)%dat = diag_data_device%mLayerCompress(:,iGRU)
  diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSoilCompress)%dat(1) = diag_data_device%scalarSoilCompress(iGRU)
  diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerMatricHeadLiq)%dat = diag_data_device%mLayerMatricHeadLiq(:,iGRU)
  diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilLiq)%dat(1) = diag_data_device%scalarTotalSoilLiq(iGRU)
  diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilIce)%dat(1) = diag_data_device%scalarTotalSoilIce(iGRU)
  diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilWat)%dat(1) = diag_data_device%scalarTotalSoilWat(iGRU)
      diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVGn_m)%dat = diag_data_device%scalarVGn_m(:,iGRU)

    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarKappa)%dat(1) = diag_data_device%scalarKappa(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolLatHt_fus)%dat(1) = diag_data_device%scalarVolLatHt_fus(iGRU)

    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%numFluxCalls)%dat(1) = diag_data_device%numFluxCalls
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%wallClockTime)%dat(1) = diag_data_device%wallClockTime
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%meanStepSize)%dat(1) = diag_data_device%meanStepSize


    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceCasNrg)%dat(1) = diag_data_device%balanceCasNrg(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceVegNrg)%dat(1) = diag_data_device%balanceVegNrg(iGRU)
    if (size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerNrg)%dat).ne.nLayers) then
      deallocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerNrg)%dat)
      allocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerNrg)%dat(1:nLayers))
    end if
      diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerNrg)%dat(1:nLayers) = diag_data_device%balanceLayerNrg(1:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSnowNrg)%dat(1) = diag_data_device%balanceSnowNrg(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSoilNrg)%dat(1) = diag_data_device%balanceSoilNrg(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceVegMass)%dat(1) = diag_data_device%balanceVegMass(iGRU)
    if (size(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerMass)%dat).ne.nLayers) then
      deallocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerMass)%dat)
      allocate(diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerMass)%dat(1:nLayers))
    end if
      diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerMass)%dat(1:nLayers) = diag_data_device%balanceLayerMass(1:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSnowMass)%dat(1) = diag_data_device%balanceSnowMass(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSoilMass)%dat(1) = diag_data_device%balanceSoilMass(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceAqMass)%dat(1) = diag_data_device%balanceAqMass(iGRU)

    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%numSteps)%dat(1) = diag_data_device%numSteps
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%numResEvals)%dat(1) = diag_data_device%numResEvals
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%numLinSolvSetups)%dat(1) = diag_data_device%numLinSolvSetups
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%numErrTestFails)%dat(1) = diag_data_device%numErrTestFails
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%kLast)%dat(1) = diag_data_device%kLast
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%kCur)%dat(1) = diag_data_device%kCur
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%hInitUsed)%dat(1) = diag_data_device%hInitUsed
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%hLast)%dat(1) = diag_data_device%hLast
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%hCur)%dat(1) = diag_data_device%hCur
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%tCur)%dat(1) = diag_data_device%tCur
  end do
  end subroutine finalize_device_diag_data

  subroutine deallocate_device_diag_data(diag_data_device)
    type(diag_data_device), intent(inout) :: diag_data_device
    deallocate(diag_data_device%scalarCanopyDepth)
    deallocate(diag_data_device%scalarBulkVolHeatCapVeg)
    deallocate(diag_data_device%scalarCanopyCm)
    deallocate(diag_data_device%scalarCanopyEmissivity)
    deallocate(diag_data_device%scalarRootZoneTemp)
    deallocate(diag_data_device%scalarLAI)
    deallocate(diag_data_device%scalarSAI)
    deallocate(diag_data_device%scalarExposedLAI)
    deallocate(diag_data_device%scalarExposedSAI)
    deallocate(diag_data_device%scalarAdjMeasHeight)
    deallocate(diag_data_device%scalarCanopyIceMax)
    deallocate(diag_data_device%scalarCanopyLiqMax)
    deallocate(diag_data_device%scalarGrowingSeasonIndex)
    deallocate(diag_data_device%scalarVolHtCap_air)
    deallocate(diag_data_device%scalarVolHtCap_ice)
    deallocate(diag_data_device%scalarVolHtCap_soil)
    deallocate(diag_data_device%scalarVolHtCap_water)
    deallocate(diag_data_device%mLayerVolHtCapBulk)
    deallocate(diag_data_device%mLayerCm)
    deallocate(diag_data_device%scalarLambda_drysoil)
    deallocate(diag_data_device%scalarLambda_wetsoil)
    deallocate(diag_data_device%mLayerThermalC)
    deallocate(diag_data_device%iLayerThermalC)
    deallocate(diag_data_device%scalarCanopyEnthTemp)
    deallocate(diag_data_device%mLayerEnthTemp)
    deallocate(diag_data_device%scalarTotalSoilEnthalpy)
    deallocate(diag_data_device%scalarTotalSnowEnthalpy)
    deallocate(diag_data_device%scalarVPair)
    deallocate(diag_data_device%scalarVP_CanopyAir)
    deallocate(diag_data_device%scalarTwetbulb)
    deallocate(diag_data_device%scalarSnowfallTemp)
    deallocate(diag_data_device%scalarNewSnowDensity)
    deallocate(diag_data_device%scalarO2air)
    deallocate(diag_data_device%scalarCO2air)
    deallocate(diag_data_device%windspd_x)
    deallocate(diag_data_device%windspd_y)
    deallocate(diag_data_device%scalarCosZenith)
    deallocate(diag_data_device%scalarFractionDirect)
    deallocate(diag_data_device%scalarCanopySunlitFraction)
    deallocate(diag_data_device%scalarCanopySunlitLAI)
    deallocate(diag_data_device%scalarCanopyShadedLAI)
    deallocate(diag_data_device%spectralAlbGndDirect)
    deallocate(diag_data_device%spectralAlbGndDiffuse)
    deallocate(diag_data_device%scalarGroundAlbedo)
    deallocate(diag_data_device%scalarLatHeatSubVapCanopy)
    deallocate(diag_data_device%scalarLatHeatSubVapGround)
    deallocate(diag_data_device%scalarSatVP_CanopyTemp)
    deallocate(diag_data_device%scalarSatVP_GroundTemp)
    deallocate(diag_data_device%scalarZ0Canopy)
    deallocate(diag_data_device%scalarWindReductionFactor)
    deallocate(diag_data_device%scalarZeroPlaneDisplacement)
    deallocate(diag_data_device%scalarRiBulkCanopy)
    deallocate(diag_data_device%scalarRiBulkGround)
    deallocate(diag_data_device%scalarCanopyStabilityCorrection)
    deallocate(diag_data_device%scalarGroundStabilityCorrection)
    deallocate(diag_data_device%scalarIntercellularCO2Sunlit)
    deallocate(diag_data_device%scalarIntercellularCO2Shaded)
    deallocate(diag_data_device%scalarTranspireLim)
    deallocate(diag_data_device%scalarTranspireLimAqfr)
    deallocate(diag_data_device%scalarFoliageNitrogenFactor)
    deallocate(diag_data_device%scalarSoilRelHumidity)
    deallocate(diag_data_device%mLayerTranspireLim)
    deallocate(diag_data_device%mLayerRootDensity)
    deallocate(diag_data_device%scalarAquiferRootFrac)
    deallocate(diag_data_device%scalarFracLiqVeg)
    deallocate(diag_data_device%scalarCanopyWetFraction)
    deallocate(diag_data_device%scalarSnowAge)
    deallocate(diag_data_device%scalarGroundSnowFraction)
    deallocate(diag_data_device%spectralSnowAlbedoDirect)
    deallocate(diag_data_device%mLayerFracLiqSnow)
    deallocate(diag_data_device%mLayerThetaResid)
    deallocate(diag_data_device%mLayerPoreSpace)
    deallocate(diag_data_device%mLayerMeltFreeze)
    deallocate(diag_data_device%scalarInfilArea)
    deallocate(diag_data_device%scalarFrozenArea)
    deallocate(diag_data_device%scalarSoilControl)
    deallocate(diag_data_device%mLayerVolFracAir)
    deallocate(diag_data_device%mLayerTcrit)
    deallocate(diag_data_device%mLayerCompress)
    deallocate(diag_data_device%scalarSoilCompress)
    deallocate(diag_data_device%mLayerMatricHeadLiq)
    deallocate(diag_data_device%scalarTotalSoilLiq)
    deallocate(diag_data_device%scalarTotalSoilIce)
    deallocate(diag_data_device%scalarTotalSoilWat)

    deallocate(diag_data_device%scalarVGn_m)
    deallocate(diag_data_device%scalarKappa)
    deallocate(diag_data_device%scalarVolLatHt_fus)

    deallocate(diag_data_device%balanceCasNrg)
    deallocate(diag_data_device%balanceVegNrg)
    deallocate(diag_data_device%balanceLayerNrg)
    deallocate(diag_data_device%balanceSnowNrg)
    deallocate(diag_data_device%balanceSoilNrg)
    deallocate(diag_data_device%balanceVegMass)
    deallocate(diag_data_device%balanceLayerMass)
    deallocate(diag_data_device%balanceSnowMass)
    deallocate(diag_data_device%balanceSoilMass)
    deallocate(diag_data_device%balanceAqMass)
       
  end subroutine deallocate_device_diag_data
      

  subroutine allocate_device_bvar_data(bvar_data_device, bvar_data,nGRU)
    implicit none
    type(bvar_data_device),intent(inout) :: bvar_data_device
    type(gru_doubleVec),intent(in) :: bvar_data
    integer(i4b),intent(in) :: nGRU
    integer(i4b) :: iGRU

    allocate(bvar_data_device%basin__AquiferStorage(nGRU))
    allocate(bvar_data_device%basin__SurfaceRunoff(nGRU))
allocate(bvar_data_device%basin__ColumnOutflow(nGRU))
allocate(bvar_data_device%basin__AquiferRecharge(nGRU))
allocate(bvar_data_device%basin__AquiferBaseflow(nGRU))
allocate(bvar_data_device%basin__AquiferTranspire(nGRU))
allocate(bvar_data_device%basin__TotalRunoff(nGRU))
allocate(bvar_data_device%basin__SoilDrainage(nGRU))
allocate(bvar_data_device%basin__totalArea(nGRU))
allocate(bvar_data_device%routingRunoffFuture(size(bvar_data%gru(1)%var(iLookBVAR%routingRunoffFuture)%dat),nGRU))
allocate(bvar_data_device%routingFractionFuture(size(bvar_data%gru(1)%var(iLookBVAR%routingFractionFuture)%dat),nGRU))
allocate(bvar_data_device%averageInstantRunoff(nGRU))
allocate(bvar_data_device%averageRoutedRunoff(nGRU))
do iGRU=1,nGRU
    bvar_data_device%basin__SurfaceRunoff(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)
    bvar_data_device%basin__AquiferStorage(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferStorage)%dat(1)
bvar_data_device%basin__ColumnOutflow(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1)
bvar_data_device%basin__AquiferRecharge(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferRecharge)%dat(1)
bvar_data_device%basin__AquiferBaseflow(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)
bvar_data_device%basin__AquiferTranspire(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferTranspire)%dat(1)
bvar_data_device%basin__TotalRunoff(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__TotalRunoff)%dat(1)
bvar_data_device%basin__SoilDrainage(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__SoilDrainage)%dat(1)
bvar_data_device%basin__totalArea(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)
  bvar_data_device%routingRunoffFuture(:,iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%routingRunoffFuture)%dat
  bvar_data_device%routingFractionFuture(:,iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%routingFractionFuture)%dat
  bvar_data_device%averageInstantRunoff(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%averageInstantRunoff)%dat(1)
  bvar_data_device%averageRoutedRunoff(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(1)
end do
  end subroutine allocate_device_bvar_data

  subroutine finalize_device_bvar_data(bvar_data_device, bvar_data)
    type(bvar_data_device),intent(inout) :: bvar_data_device
    type(gru_doubleVec),intent(in) :: bvar_data
    integer(i4b) :: iGRU,nGRU
    nGRU = size(bvar_data%gru)
    do iGRU=1,nGRU
        
    bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferStorage)%dat(1) = bvar_data_device%basin__AquiferStorage(iGRU)
    bvar_data%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1) = bvar_data_device%basin__SurfaceRunoff(iGRU)
bvar_data%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1) = bvar_data_device%basin__ColumnOutflow(iGRU)
bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferRecharge)%dat(1) = bvar_data_device%basin__AquiferRecharge(iGRU)
bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1) = bvar_data_device%basin__AquiferBaseflow(iGRU)
bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = bvar_data_device%basin__AquiferTranspire(iGRU)
bvar_data%gru(iGRU)%var(iLookBVAR%basin__TotalRunoff)%dat(1) = bvar_data_device%basin__TotalRunoff(iGRU)
bvar_data%gru(iGRU)%var(iLookBVAR%basin__SoilDrainage)%dat(1) = bvar_data_device%basin__SoilDrainage(iGRU)
bvar_data%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1) = bvar_data_device%basin__totalArea(iGRU)
  bvar_data%gru(iGRU)%var(iLookBVAR%routingRunoffFuture)%dat = bvar_data_device%routingRunoffFuture(:,iGRU)
  bvar_data%gru(iGRU)%var(iLookBVAR%routingFractionFuture)%dat = bvar_data_device%routingFractionFuture(:,iGRU)
bvar_data%gru(iGRU)%var(iLookBVAR%averageInstantRunoff)%dat(1) = bvar_data_device%averageInstantRunoff(iGRU)
bvar_data%gru(iGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(1) = bvar_data_device%averageRoutedRunoff(iGRU)
    end do
  end subroutine finalize_device_bvar_data

  subroutine deallocate_device_bvar_data(bvar_data_device)
    type(bvar_data_device),intent(inout) :: bvar_data_device
        
    deallocate(bvar_data_device%basin__AquiferStorage)
  end subroutine deallocate_device_bvar_data

  subroutine allocate_device_forc_data(forc_data_device, forc_data,nGRU)
    type(forc_data_device),intent(inout) :: forc_data_device
    type(gru_hru_double),intent(in) :: forc_data
    integer(i4b) :: nGRU,iGRU

    allocate(forc_data_device%time(nGRU))
    allocate(forc_data_device%pptrate(nGRU))
    allocate(forc_data_device%SWRadAtm(nGRU))
    allocate(forc_data_device%LWRadAtm_d(nGRU))
    allocate(forc_data_device%airtemp_d(nGRU))
    allocate(forc_data_device%windspd_d(nGRU))
    allocate(forc_data_device%airpres_d(nGRU))
    allocate(forc_data_device%spechum(nGRU))

    do iGRU=1,nGRU
      forc_data_device%time(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%time)
      forc_data_device%pptrate(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%pptrate)
      forc_data_device%SWRadAtm(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%SWRadAtm)
      forc_data_device%LWRadAtm_d(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%LWRadAtm)
      forc_data_device%airtemp_d(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%airtemp)
      forc_data_device%windspd_d(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%windspd)
      forc_data_device%airpres_d(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%airpres)
      forc_data_device%spechum(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%spechum)
    end do

  end subroutine allocate_device_forc_data
          
  subroutine finalize_device_forc_data(forc_data_device, forc_data)
    type(forc_data_device),intent(inout) :: forc_data_device
    type(gru_hru_double),intent(in) :: forc_data
    integer(i4b) :: nGRU,iGRU

    nGRU = size(forc_data%gru)

    do iGRU=1,nGRU
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%time) = forc_data_device%time(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%pptrate) = forc_data_device%pptrate(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%SWRadAtm) = forc_data_device%SWRadAtm(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%LWRadAtm) = forc_data_device%LWRadAtm_d(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%airtemp) = forc_data_device%airtemp_d(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%windspd) = forc_data_device%windspd_d(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%airpres) = forc_data_device%airpres_d(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%spechum) = forc_data_device%spechum(iGRU)
    end do

  end subroutine finalize_device_forc_data

  subroutine deallocate_device_forc_data(forc_data_device)
    type(forc_data_device),intent(inout) :: forc_data_device

    deallocate(forc_data_device%time)
    deallocate(forc_data_device%pptrate)
    deallocate(forc_data_device%SWRadAtm)
    deallocate(forc_data_device%LWRadAtm_d)
    deallocate(forc_data_device%airtemp_d)
    deallocate(forc_data_device%windspd_d)
    deallocate(forc_data_device%airpres_d)
    deallocate(forc_data_device%spechum)
  end subroutine deallocate_device_forc_data
    
  subroutine allocate_device_decisions(device_decisions)
    USE globalData,only:model_decisions ! model decision structure
    USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
    use globalData,only:urbanVegCategory
    use globalData,only:data_step
    use globalData,only:refJulDay
    use summaFileManager,only:NC_TIME_ZONE
    implicit none

    type(decisions_device) :: device_decisions
    device_decisions%true = .true.
    device_decisions%false = .false.
    device_decisions%one = 1
        
    device_decisions%soilCatTbl = model_decisions(iLookDECISIONS%soilCatTbl)%iDecision
    device_decisions%vegeParTbl = model_decisions(iLookDECISIONS%vegeParTbl)%iDecision
    device_decisions%soilStress = model_decisions(iLookDECISIONS%soilStress)%iDecision
    device_decisions%stomResist = model_decisions(iLookDECISIONS%stomResist)%iDecision
    device_decisions%bbTempFunc = model_decisions(iLookDECISIONS%bbTempFunc)%iDecision
    device_decisions%bbHumdFunc = model_decisions(iLookDECISIONS%bbHumdFunc)%iDecision
    device_decisions%bbElecFunc = model_decisions(iLookDECISIONS%bbElecFunc)%iDecision
    device_decisions%bbCO2point = model_decisions(iLookDECISIONS%bbCO2point)%iDecision
    device_decisions%bbNumerics = model_decisions(iLookDECISIONS%bbNumerics)%iDecision
    device_decisions%bbAssimFnc = model_decisions(iLookDECISIONS%bbAssimFnc)%iDecision
    device_decisions%bbCanIntg8 = model_decisions(iLookDECISIONS%bbCanIntg8)%iDecision
    device_decisions%num_method = model_decisions(iLookDECISIONS%num_method)%iDecision
    device_decisions%fDerivMeth = model_decisions(iLookDECISIONS%fDerivMeth)%iDecision
    device_decisions%LAI_method = model_decisions(iLookDECISIONS%LAI_method)%iDecision
    device_decisions%cIntercept = model_decisions(iLookDECISIONS%cIntercept)%iDecision
    device_decisions%f_Richards = model_decisions(iLookDECISIONS%f_Richards)%iDecision
    device_decisions%groundwatr = model_decisions(iLookDECISIONS%groundwatr)%iDecision
    device_decisions%hc_profile = model_decisions(iLookDECISIONS%hc_profile)%iDecision
    device_decisions%bcUpprTdyn = model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision
    device_decisions%bcLowrTdyn = model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision
    device_decisions%bcUpprSoiH = model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision
    device_decisions%bcLowrSoiH = model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision
    device_decisions%veg_traits = model_decisions(iLookDECISIONS%veg_traits)%iDecision
    device_decisions%rootProfil = model_decisions(iLookDECISIONS%rootProfil)%iDecision
    device_decisions%canopyEmis = model_decisions(iLookDECISIONS%canopyEmis)%iDecision
    device_decisions%snowIncept = model_decisions(iLookDECISIONS%snowIncept)%iDecision
    device_decisions%snowUnload = model_decisions(iLookDECISIONS%snowUnload)%iDecision
    device_decisions%windPrfile = model_decisions(iLookDECISIONS%windPrfile)%iDecision
    device_decisions%astability = model_decisions(iLookDECISIONS%astability)%iDecision
    device_decisions%canopySrad = model_decisions(iLookDECISIONS%canopySrad)%iDecision
    device_decisions%alb_method = model_decisions(iLookDECISIONS%alb_method)%iDecision
    device_decisions%snowLayers = model_decisions(iLookDECISIONS%snowLayers)%iDecision
    device_decisions%compaction = model_decisions(iLookDECISIONS%compaction)%iDecision
    device_decisions%thCondSnow = model_decisions(iLookDECISIONS%thCondSnow)%iDecision
    device_decisions%thCondSoil = model_decisions(iLookDECISIONS%thCondSoil)%iDecision
    device_decisions%spatial_gw = model_decisions(iLookDECISIONS%spatial_gw)%iDecision
    device_decisions%subRouting = model_decisions(iLookDECISIONS%subRouting)%iDecision
    device_decisions%snowDenNew = model_decisions(iLookDECISIONS%snowDenNew)%iDecision
    device_decisions%nrgConserv = model_decisions(iLookDECISIONS%nrgConserv)%iDecision
    device_decisions%aquiferIni = model_decisions(iLookDECISIONS%aquiferIni)%iDecision
    device_decisions%infRateMax = model_decisions(iLookDECISIONS%infRateMax)%iDecision
    device_decisions%surfRun_IE = model_decisions(iLookDECISIONS%surfRun_IE)%iDecision
    device_decisions%surfRun_SE = model_decisions(iLookDECISIONS%surfRun_SE)%iDecision
    device_decisions%urbanVegCategory = urbanVegCategory

    device_decisions%data_step = data_step
    device_decisions%refJulDay = refJulDay
    device_decisions%NC_TIME_ZONE = NC_TIME_ZONE
  end subroutine

  subroutine deallocate_device_decisions(device_decisions)
    type(decisions_device) :: device_decisions
    deallocate(device_decisions%bcUpprTdyn)
    deallocate(device_decisions%bcLowrTdyn)
    deallocate(device_decisions%cIntercept)
    deallocate(device_decisions%veg_traits)
    deallocate(device_decisions%canopyEmis)
    deallocate(device_decisions%windPrfile)
    deallocate(device_decisions%astability)
    deallocate(device_decisions%soilStress)
    deallocate(device_decisions%groundwatr)
    deallocate(device_decisions%stomResist)
    deallocate(device_decisions%spatial_gw)
    deallocate(device_decisions%bbTempFunc)
    deallocate(device_decisions%bbHumdFunc)
    deallocate(device_decisions%bbElecFunc)
    deallocate(device_decisions%bbCO2point)
    deallocate(device_decisions%bbNumerics)
    deallocate(device_decisions%bbAssimFnc)
    deallocate(device_decisions%bbCanIntg8)
    deallocate(device_decisions%f_Richards)
    deallocate(device_decisions%true)
    deallocate(device_decisions%one)
    deallocate(device_decisions%false)
    deallocate(device_decisions%bcUpprSoiH)
    deallocate(device_decisions%nrgConserv)
    deallocate(device_decisions%thCondSnow)
    deallocate(device_decisions%thCondSoil)

  end subroutine

  subroutine get_device_type_data(type_data, iVar, iGRU,output)
    type(type_data_device) :: type_data
    integer(i4b),intent(in) :: iVar, iGRU
    integer(i4b),allocatable,intent(inout) :: output(:)
    if (allocated(output)) deallocate(output)
    allocate(output(1))
    select case (iVar)
    case(iLookTYPE%vegTypeIndex); output = type_data%vegTypeIndex(iGRU)
    case(iLookTYPE%soilTypeIndex); output = type_data%soilTypeIndex(iGRU)
    case(iLookTYPE%slopeTypeIndex); output = type_data%slopeTypeIndex(iGRU)
    case(iLookTYPE%downHRUindex); output = type_data%downHRUindex(iGRU)
    end select
  end subroutine

subroutine allocate_veg_table(veg_param_tables)
  use module_sf_noahlsm
  type(veg_param_tables) :: veg_param_tables
  veg_param_tables%rgltbl = rgltbl
  veg_param_tables%rstbl = rstbl
  veg_param_tables%hstbl = hstbl
  veg_param_tables%TOPT_DATA = TOPT_DATA
  veg_param_tables%RSMAX_DATA = RSMAX_DATA
end subroutine

  subroutine allocate_veg_parameters(veg_param,nGRU)
    use NOAHMP_VEG_PARAMETERS
    use noahmp_globals
    use NOAHMP_RAD_PARAMETERS
    implicit none
    type(veg_parameters) :: veg_param
    integer(i4b) :: nGRU
    integer(i4b) :: iGRU

    veg_param%bp = bp
    veg_param%mp = mp
    veg_param%c3psn = c3psn
    veg_param%avcmx = avcmx
    veg_param%vcmx25 = vcmx25
    veg_param%ako = ako
    veg_param%ko25 = ko25
    veg_param%akc = akc
    veg_param%kc25 = kc25
    veg_param%qe25 = qe25
    veg_param%folnmx = folnmx
    allocate(veg_param%rgl_d(nGRU))
    ! veg_param%rgl_d = rgl
    allocate(veg_param%rsmin_d(nGRU))
    ! veg_param%rsmin_d = rsmin
    allocate(veg_param%rsmax_d(nGRU))
    ! veg_param%rsmax = rsmax
    allocate(veg_param%topt_d(nGRU))
    ! veg_param%topt = topt
    allocate(veg_param%hs_d(nGRU))
    ! veg_param%hs_d = hs
    veg_param%rhol = rhol
    veg_param%rhos = rhos
    veg_param%taul = taul
    veg_param%taus = taus
    veg_param%opt_alb = opt_alb
    veg_param%omegas = omegas
    veg_param%opt_rad = opt_rad
    veg_param%tfrz = tfrz
    veg_param%betais = betais
    veg_param%betads = betads
    veg_param%xl = xl
    allocate(veg_param%hvt_(size(hvt),nGRU))
    allocate(veg_param%hvb_(size(hvb),nGRU))

    do iGRU=1,nGRU
      veg_param%hvt_(:,iGRU) = hvt
      veg_param%hvb_(:,iGRU) = hvb
    end do
    veg_param%rc = rc
    veg_param%swemx = swemx
    veg_param%albsat = albsat
    veg_param%albdry = albdry
    veg_param%alblak = alblak
    veg_param%isbarren = isbarren
    veg_param%ISSNOW = ISSNOW
    veg_param%ISWATER = ISWATER
    allocate(veg_param%saim_(size(saim,1),size(saim,2),nGRU))
    allocate(veg_param%laim_(size(laim,1),size(laim,2),nGRU))

    do iGRU=1,nGRU
    veg_param%saim_(:,:,iGRU) = saim
    veg_param%laim_(:,:,iGRU) = laim
    end do
    veg_param%dveg = dveg
    veg_param%tmin = tmin

  end subroutine
  subroutine allocate_device_attr_data(attr_data_device,attr_data,nGRU)
    type(attr_data_device) :: attr_data_device
    type(gru_hru_double) :: attr_data
    integer(i4b) :: nGRU,iGRU
    allocate(attr_data_device%latitude(nGRU))
    allocate(attr_data_device%longitude(nGRU))
    allocate(attr_data_device%elevation(nGRU))
    allocate(attr_data_device%tan_slope(nGRU))
    allocate(attr_data_device%contourLength(nGRU))
    allocate(attr_data_device%HRUarea(nGRU))
    allocate(attr_data_device%mHeight(nGRU))
    allocate(attr_data_device%aspect(nGRU))
    do iGRU=1,nGRU
      attr_data_device%latitude(iGRU) = attr_data%gru(iGRU)%hru(1)%var(iLookATTR%latitude)
      attr_data_device%longitude(iGRU) = attr_data%gru(iGRU)%hru(1)%var(iLookATTR%longitude)
      attr_data_device%elevation(iGRU) = attr_data%gru(iGRU)%hru(1)%var(iLookATTR%elevation)
      attr_data_device%tan_slope(iGRU) = attr_data%gru(iGRU)%hru(1)%var(iLookATTR%tan_slope)
      attr_data_device%contourLength(iGRU) = attr_data%gru(iGRU)%hru(1)%var(iLookATTR%contourLength)
      attr_data_device%HRUarea(iGRU) = attr_data%gru(iGRU)%hru(1)%var(iLookATTR%HRUarea)
      attr_data_device%mHeight(iGRU) = attr_data%gru(iGRU)%hru(1)%var(iLookATTR%mHeight)
      print*, iGRU, attr_data%gru(iGRU)%hru(1)%var(iLookATTR%mHeight)
      attr_data_device%aspect(iGRU) = attr_data%gru(iGRU)%hru(1)%var(iLookATTR%aspect)
    end do
  end subroutine

  subroutine deallocate_device_attr_data(attr_data_device)
    type(attr_data_device) :: attr_data_device
    deallocate(attr_data_device%latitude)
    deallocate(attr_data_device%longitude)
    deallocate(attr_data_device%elevation)
    deallocate(attr_data_device%tan_slope)
    deallocate(attr_data_device%contourLength)
    deallocate(attr_data_device%HRUarea)
    deallocate(attr_data_device%mHeight)
    deallocate(attr_data_device%aspect)
  end subroutine

  subroutine update_flux_sum(flux_sum, flux_data, flux_prev, weight,&
    nGRU)
    real(rkind) :: weight
    type(flux_data_device) :: flux_sum, flux_data, flux_prev
    integer(i4b),intent(in) :: nGRU

    integer(i4b) :: numFluxData, iGRU, ix
    numFluxData = flux_data%numFluxData
    associate(flux_sum_ => flux_sum%data, flux_data_ => flux_data%data, flux_prev_ => flux_prev%data)
      !$cuf kernel do(2) <<<*,*>>>
    do iGRU=1,nGRU
      do ix=1,numFluxData
        flux_sum_(ix,iGRU) = flux_sum_(ix,iGRU) + (flux_data_(ix,iGRU) + flux_prev_(ix,iGRU)) * weight
      end do
    end do
    end associate
    

  end subroutine

  attributes(device) function get_iGRU()
    integer(i4b) :: get_iGRU
    get_iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  end function

end module initialize_device      