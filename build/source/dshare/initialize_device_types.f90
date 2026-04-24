module initialize_device
  use data_types
  use device_data_types
  USE var_lookup,only:iLookFLUX        ! lookup indices for flux data
  USE var_lookup,only:iLookDERIV       ! lookup indices for derivative data
  USE var_lookup,only:iLookFORCE       ! lookup indices for forcing data 
  USE var_lookup,only:iLookDIAG        ! lookup indices for diagnostic variable data
  USE var_lookup,only:iLookPROG,iLookPARAM, iLookINDEX, iLookBVAR, iLookBPAR, iLookTYPE, iLookATTR
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
public::allocate_device_bvar_stat,finalize_device_bvar_stat,get_device_bvar_stat
public::get_device_param_data
public::get_device_attr_data,set_device_attr_data
public::set_device_param_data_scalar
public::set_device_param_data
public::get_device_type_data
public::allocate_device_bvar_data,allocate_device_prog_data,allocate_device_prog_data_snow,allocate_device_diag_data
public::deallocate_device_prog_data
public::deallocate_device_bvar_data,get_device_bvar_data
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
public::get_device_flux_data,get_device_diag_data,get_device_prog_data
public::allocate_device_indx_data_snow,allocate_device_diag_data_snow
public::get_device_forc_data,get_device_indx_data
public::get_device_prog_stat,allocate_device_prog_stat
public::get_device_indx_stat,allocate_device_indx_stat
public::get_device_flux_stat,allocate_device_flux_stat
public::get_device_diag_stat,allocate_device_diag_stat
public::get_device_forc_stat,allocate_device_forc_stat


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
    case(iLookPARAM%upperBoundHead); mparStruct%upperBoundHead_(iGRU) = parVector
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

  subroutine get_device_flux_data(fluxStruct, iVar, iGRU, parVector)
    implicit none
    type(flux_data_device) :: fluxStruct
    integer(i4b) :: iGRU, iVar
    real(rkind),allocatable :: parVector(:)


    if (allocated(parVector)) deallocate(parVector)
    select case(iVar)
    ! net energy and mass fluxes for the vegetation domain
    case(iLookFLUX%scalarCanairNetNrgFlux); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixScalarCanairNetNrgFlux,iGRU)
    case(iLookFLUX%scalarCanopyNetNrgFlux); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixScalarCanopyNetNrgFlux,iGRU)
    case(iLookFLUX%scalarGroundNetNrgFlux); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixScalarGroundNetNrgFlux,iGRU)
    case(iLookFLUX%scalarCanopyNetLiqFlux); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixScalarCanopyNetLiqFlux,iGRU)
    ! forcing
    case(iLookFLUX%scalarRainfall); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixScalarRainfall,iGRU)
    case(iLookFLUX%scalarSnowfall); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixScalarSnowfall,iGRU)
    ! shortwave radiation
    case(iLookFLUX%spectralIncomingDirect); parVector = fluxStruct%data(fluxStruct%ixspectralIncomingDirect_start:fluxStruct%ixspectralIncomingDirect_end,iGRU)
    case(iLookFLUX%spectralIncomingDiffuse); parVector = fluxStruct%data(fluxStruct%ixspectralIncomingDiffuse_start:fluxStruct%ixspectralIncomingDiffuse_end,iGRU)
    case(iLookFLUX%scalarCanopySunlitPAR); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarCanopySunlitPAR,iGRU)
    case(iLookFLUX%scalarCanopyShadedPAR); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarCanopyShadedPAR,iGRU)
    case(iLookFLUX%spectralBelowCanopyDirect); parVector = fluxStruct%data(fluxStruct%ixspectralBelowCanopyDirect_start:fluxStruct%ixspectralBelowCanopyDirect_end,iGRU)
    case(iLookFLUX%spectralBelowCanopyDiffuse); parVector = fluxStruct%data(fluxStruct%ixspectralBelowCanopyDiffuse_start:fluxStruct%ixspectralBelowCanopyDiffuse_end,iGRU)
    case(iLookFLUX%scalarBelowCanopySolar); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarBelowCanopySolar,iGRU)
    case(iLookFLUX%scalarCanopyAbsorbedSolar); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarCanopyAbsorbedSolar,iGRU)
    case(iLookFLUX%scalarGroundAbsorbedSolar); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarGroundAbsorbedSolar,iGRU)
    ! longwave radiation
    case(iLookFLUX%scalarLWRadCanopy); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWRadCanopy,iGRU)
    case(iLookFLUX%scalarLWRadGround); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWRadGround,iGRU)
    case(iLookFLUX%scalarLWRadUbound2Canopy); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWRadUbound2Canopy,iGRU)
    case(iLookFLUX%scalarLWRadUbound2Ground); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWRadUbound2Ground,iGRU)
    case(iLookFLUX%scalarLWRadUbound2Ubound); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWRadUbound2Ubound,iGRU)
    case(iLookFLUX%scalarLWRadCanopy2Ubound); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWRadCanopy2Ubound,iGRU)
    case(iLookFLUX%scalarLWRadCanopy2Ground); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWRadCanopy2Ground,iGRU)
    case(iLookFLUX%scalarLWRadCanopy2Canopy); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWRadCanopy2Canopy,iGRU)
    case(iLookFLUX%scalarLWRadGround2Ubound); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWRadGround2Ubound,iGRU)
    case(iLookFLUX%scalarLWRadGround2Canopy); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWRadGround2Canopy,iGRU)
    case(iLookFLUX%scalarLWNetCanopy); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWNetCanopy,iGRU)
    case(iLookFLUX%scalarLWNetGround); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWNetGround,iGRU)
    case(iLookFLUX%scalarLWNetUbound); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLWNetUbound,iGRU)
    ! turbulent heat transfer
    case(iLookFLUX%scalarEddyDiffusCanopyTop); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarEddyDiffusCanopyTop,iGRU)
    case(iLookFLUX%scalarFrictionVelocity); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarFrictionVelocity,iGRU)
    case(iLookFLUX%scalarWindspdCanopyTop); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarWindspdCanopyTop,iGRU)
    case(iLookFLUX%scalarWindspdCanopyBottom); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarWindspdCanopyBottom,iGRU)
    case(iLookFLUX%scalarGroundResistance); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarGroundResistance,iGRU)
    case(iLookFLUX%scalarCanopyResistance); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarCanopyResistance,iGRU)
    case(iLookFLUX%scalarLeafResistance); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLeafResistance,iGRU)
    case(iLookFLUX%scalarSoilResistance); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarSoilResistance,iGRU)
    case(iLookFLUX%scalarSenHeatTotal); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarSenHeatTotal,iGRU)
    case(iLookFLUX%scalarSenHeatCanopy); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarSenHeatCanopy,iGRU)
    case(iLookFLUX%scalarSenHeatGround); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarSenHeatGround,iGRU)
    case(iLookFLUX%scalarLatHeatTotal); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLatHeatTotal,iGRU)
    case(iLookFLUX%scalarLatHeatCanopyEvap); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLatHeatCanopyEvap,iGRU)
    case(iLookFLUX%scalarLatHeatCanopyTrans); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLatHeatCanopyTrans,iGRU)
    case(iLookFLUX%scalarLatHeatGround); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarLatHeatGround,iGRU)
    case(iLookFLUX%scalarCanopyAdvectiveHeatFlux); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarCanopyAdvectiveHeatFlux,iGRU)
    case(iLookFLUX%scalarGroundAdvectiveHeatFlux); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarGroundAdvectiveHeatFlux,iGRU)
    case(iLookFLUX%scalarCanopySublimation); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarCanopySublimation,iGRU)
    case(iLookFLUX%scalarSnowSublimation); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarSnowSublimation,iGRU)
    ! liquid water fluxes associated with evapotranspiration
    case(iLookFLUX%scalarStomResistSunlit); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarStomResistSunlit,iGRU)
    case(iLookFLUX%scalarStomResistShaded); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarStomResistShaded,iGRU)
    case(iLookFLUX%scalarPhotosynthesisSunlit); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarPhotosynthesisSunlit,iGRU)
    case(iLookFLUX%scalarPhotosynthesisShaded); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarPhotosynthesisShaded,iGRU)
    case(iLookFLUX%scalarCanopyTranspiration); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarCanopyTranspiration,iGRU)
    case(iLookFLUX%scalarCanopyEvaporation); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarCanopyEvaporation,iGRU)
    case(iLookFLUX%scalarGroundEvaporation); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarGroundEvaporation,iGRU)
    case(iLookFLUX%mLayerTranspire); parVector = fluxStruct%data(fluxStruct%ixmLayerTranspire_start:fluxStruct%ixmLayerTranspire_end,iGRU)
    ! liquid and solid water fluxes through the canopy
    case(iLookFLUX%scalarThroughfallSnow); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarThroughfallSnow,iGRU)
    case(iLookFLUX%scalarThroughfallRain); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarThroughfallRain,iGRU)
    case(iLookFLUX%scalarCanopySnowUnloading); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarCanopySnowUnloading,iGRU)
    case(iLookFLUX%scalarCanopyLiqDrainage); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarCanopyLiqDrainage,iGRU)
    case(iLookFLUX%scalarCanopyMeltFreeze); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarCanopyMeltFreeze,iGRU)
    ! energy fluxes and for the snow and soil domains
    case(iLookFLUX%iLayerConductiveFlux); parVector = fluxStruct%data(fluxStruct%ixiLayerConductiveFlux_start:fluxStruct%ixiLayerConductiveFlux_end,iGRU)
    case(iLookFLUX%iLayerAdvectiveFlux); parVector = fluxStruct%data(fluxStruct%ixiLayerAdvectiveFlux_start:fluxStruct%ixiLayerAdvectiveFlux_end,iGRU)
    case(iLookFLUX%iLayerNrgFlux); parVector = fluxStruct%data(fluxStruct%ixiLayerNrgFlux_start:fluxStruct%ixiLayerNrgFlux_end,iGRU)
    case(iLookFLUX%mLayerNrgFlux); parVector = fluxStruct%data(fluxStruct%ixmLayerNrgFlux_start:fluxStruct%ixmLayerNrgFlux_end,iGRU)
    ! liquid water fluxes for the snow domain
    case(iLookFLUX%scalarSnowDrainage); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarSnowDrainage,iGRU)
    case(iLookFLUX%iLayerLiqFluxSnow); parVector = fluxStruct%data(fluxStruct%ixiLayerLiqFluxSnow_start:fluxStruct%ixiLayerLiqFluxSnow_end,iGRU)
    case(iLookFLUX%mLayerLiqFluxSnow); parVector = fluxStruct%data(fluxStruct%ixmLayerLiqFluxSnow_start:fluxStruct%ixmLayerLiqFluxSnow_end,iGRU)
    ! liquid water fluxes for the soil domain
    case(iLookFLUX%scalarRainPlusMelt); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarRainPlusMelt,iGRU)
    case(iLookFLUX%scalarMaxInfilRate); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarMaxInfilRate,iGRU)
    case(iLookFLUX%scalarInfiltration); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarInfiltration,iGRU)
    case(iLookFLUX%scalarExfiltration); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarExfiltration,iGRU)
    case(iLookFLUX%scalarSurfaceRunoff); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarSurfaceRunoff,iGRU)
    case(iLookFLUX%scalarSurfaceRunoff_IE); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarSurfaceRunoff_IE,iGRU)
    case(iLookFLUX%scalarSurfaceRunoff_SE); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarSurfaceRunoff_SE,iGRU)
    case(iLookFLUX%mLayerSatHydCondMP); parVector = fluxStruct%data(fluxStruct%ixmLayerSatHydCondMP_start:fluxStruct%ixmLayerSatHydCondMP_end,iGRU)
    case(iLookFLUX%mLayerSatHydCond); parVector = fluxStruct%data(fluxStruct%ixmLayerSatHydCond_start:fluxStruct%ixmLayerSatHydCond_end,iGRU)
    case(iLookFLUX%iLayerSatHydCond); parVector = fluxStruct%data(fluxStruct%ixiLayerSatHydCond_start:fluxStruct%ixiLayerSatHydCond_end,iGRU)
    case(iLookFLUX%mLayerHydCond); parVector = fluxStruct%data(fluxStruct%ixmLayerHydCond_start:fluxStruct%ixmLayerHydCond_end,iGRU)
    case(iLookFLUX%iLayerLiqFluxSoil); parVector = fluxStruct%data(fluxStruct%ixiLayerLiqFluxSoil_start:fluxStruct%ixiLayerLiqFluxSoil_end,iGRU)
    case(iLookFLUX%mLayerLiqFluxSoil); parVector = fluxStruct%data(fluxStruct%ixmLayerLiqFluxSoil_start:fluxStruct%ixmLayerLiqFluxSoil_end,iGRU)
    case(iLookFLUX%mLayerBaseflow); parVector = fluxStruct%data(fluxStruct%ixmLayerBaseflow_start:fluxStruct%ixmLayerBaseflow_end,iGRU)
    case(iLookFLUX%mLayerColumnInflow); parVector = fluxStruct%data(fluxStruct%ixmLayerColumnInflow_start:fluxStruct%ixmLayerColumnInflow_end,iGRU)
    case(iLookFLUX%mLayerColumnOutflow); parVector = fluxStruct%data(fluxStruct%ixmLayerColumnOutflow_start:fluxStruct%ixmLayerColumnOutflow_end,iGRU)
    case(iLookFLUX%scalarSoilBaseflow); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarSoilBaseflow,iGRU)
    case(iLookFLUX%scalarSoilDrainage); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarSoilDrainage,iGRU)
    case(iLookFLUX%scalarAquiferRecharge); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarAquiferRecharge,iGRU)
    case(iLookFLUX%scalarAquiferTranspire); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarAquiferTranspire,iGRU)
    case(iLookFLUX%scalarAquiferBaseflow); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarAquiferBaseflow,iGRU)
    ! derived variables
    case(iLookFLUX%scalarTotalET); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarTotalET,iGRU)
    case(iLookFLUX%scalarTotalRunoff); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarTotalRunoff,iGRU)
    case(iLookFLUX%scalarNetRadiation); allocate(parVector(1)); parVector = fluxStruct%data(fluxStruct%ixscalarNetRadiation,iGRU)

    end select
    if (.not.allocated(parVector)) then
      allocate(parVector(1)); parVector = -9999_rkind
    end if
  end subroutine

  subroutine get_device_forc_data(forcStruct,iVar,iGRU,parVector)
    implicit none
    type(forc_data_device) :: forcStruct
    integer(i4b) :: iGRU, iVar
    real(rkind),allocatable :: parVector(:)
    if (allocated(parVector)) deallocate(parVector)
    select case(iVar)
    case(iLookFORCE%time); allocate(parVector(1)); parVector = forcStruct%time_h
    case(iLookFORCE%pptrate); allocate(parVector(1)); parVector = forcStruct%pptrate(iGRU)
    case(iLookFORCE%SWRadAtm); allocate(parVector(1)); parVector = forcStruct%SWRadAtm(iGRU)
    case(iLookFORCE%LWRadAtm); allocate(parVector(1)); parVector = forcStruct%LWRadAtm(iGRU)
    case(iLookFORCE%airtemp); allocate(parVector(1)); parVector = forcStruct%airtemp(iGRU)
    case(iLookFORCE%windspd); allocate(parVector(1)); parVector = forcStruct%windspd(iGRU)
    case(iLookFORCE%airpres); allocate(parVector(1)); parVector = forcStruct%airpres(iGRU)
    case(iLookFORCE%spechum); allocate(parVector(1)); parVector = forcStruct%spechum(iGRU)
    end select
    if (.not.allocated(parVector)) then
      allocate(parVector(1)); parVector = -9999_rkind
    end if

  end subroutine
  subroutine get_device_indx_data(indxStruct,iVar,iGRU,parVector)
    implicit none
    type(indx_data_device) :: indxStruct
    integer(i4b) :: iGRU, iVar
    integer(i4b),allocatable :: parVector(:)
    if (allocated(parVector)) deallocate(parVector)
    select case(iVar)
    case(iLookINDEX%nSnow); allocate(parVector(1)); parVector = indxStruct%nSnow(iGRU)
    case(iLookINDEX%nSoil); allocate(parVector(1)); parVector = indxStruct%nSoil
    case(iLookINDEX%nLayers); allocate(parVector(1)); parVector = indxStruct%nLayers(iGRU)
    case(iLookINDEX%layerType); parVector = indxStruct%layerType(:,iGRU)
    case(iLookINDEX%nCasNrg); allocate(parVector(1)); parVector = indxStruct%nCasNrg(iGRU)
    case(iLookINDEX%nVegNrg); allocate(parVector(1)); parVector = indxStruct%nVegNrg(iGRU)
    case(iLookINDEX%nVegMass); allocate(parVector(1)); parVector = indxStruct%nVegMass(iGRU)
    case(iLookINDEX%nVegState); allocate(parVector(1)); parVector = indxStruct%nVegState(iGRU)
    case(iLookINDEX%nNrgState); allocate(parVector(1)); parVector = indxStruct%nNrgState(iGRU)
    case(iLookINDEX%nWatState); allocate(parVector(1)); parVector = indxStruct%nWatState(iGRU)
    case(iLookINDEX%nMatState); allocate(parVector(1)); parVector = indxStruct%nMatState(iGRU)
    case(iLookINDEX%nMassState); allocate(parVector(1)); parVector = indxStruct%nMassState(iGRU)
    case(iLookINDEX%nState); allocate(parVector(1)); parVector = indxStruct%nState(iGRU)
    ! case(iLookINDEX%nSnowSoilNrg); allocate(parVector(1)); parVector = indxStruct%nSnowSoilNrg(iGRU)
    ! case(iLookINDEX%nSnowOnlyNrg); allocate(parVector(1)); parVector = indxStruct%nSnowOnlyNrg(iGRU)
    ! case(iLookINDEX%nSoilOnlyNrg); allocate(parVector(1)); parVector = indxStruct%nSoilOnlyNrg(iGRU)
    ! case(iLookINDEX%nSnowSoilHyd); allocate(parVector(1)); parVector = indxStruct%nSnowSoilHyd(iGRU)
    ! case(iLookINDEX%nSnowOnlyHyd); allocate(parVector(1)); parVector = indxStruct%nSnowOnlyHyd(iGRU)
    ! case(iLookINDEX%nSoilOnlyHyd); allocate(parVector(1)); parVector = indxStruct%nSoilOnlyHyd(iGRU)
    case(iLookINDEX%ixControlVolume); parVector = indxStruct%ixControlVolume(:,iGRU)
    case(iLookINDEX%ixDomainType); parVector = indxStruct%ixDomainType(:,iGRU)
    case(iLookINDEX%ixStateType); parVector = indxStruct%ixStateType(:,iGRU)
    case(iLookINDEX%ixHydType); parVector = indxStruct%ixHydType(:,iGRU)
    ! case(iLookINDEX%ixDomainType_subset); parVector = indxStruct%ixDomainType_subset(:,iGRU)
    ! case(iLookINDEX%ixStateType_subset); parVector = indxStruct%ixStateType_subset(:,iGRU)
    ! case(iLookINDEX%ixMapFull2Subset); parVector = indxStruct%ixMapFull2Subset(:,iGRU)
    ! case(iLookINDEX%ixMapSubset2Full); parVector = indxStruct%ixMapSubset2Full(:,iGRU)
    case(iLookINDEX%ixCasNrg); allocate(parVector(1)); parVector = indxStruct%ixCasNrg(iGRU)
    case(iLookINDEX%ixVegNrg); allocate(parVector(1)); parVector = indxStruct%ixVegNrg(iGRU)
    case(iLookINDEX%ixVegHyd); allocate(parVector(1)); parVector = indxStruct%ixVegHyd(iGRU)
    case(iLookINDEX%ixTopNrg); allocate(parVector(1)); parVector = indxStruct%ixTopNrg(iGRU)
    case(iLookINDEX%ixTopHyd); allocate(parVector(1)); parVector = indxStruct%ixTopHyd(iGRU)
    case(iLookINDEX%ixAqWat); allocate(parVector(1)); parVector = indxStruct%ixAqWat(iGRU)
    case(iLookINDEX%ixNrgOnly); parVector = indxStruct%ixNrgOnly(:,iGRU)
    ! case(iLookINDEX%ixHydOnly); parVector = indxStruct%ixHydOnly(:,iGRU)
    ! case(iLookINDEX%ixMatOnly); parVector = indxStruct%ixMatOnly(:,iGRU)
    ! case(iLookINDEX%ixMassOnly); parVector = indxStruct%ixMassOnly(:,iGRU)
    case(iLookINDEX%ixSnowSoilNrg); parVector = indxStruct%ixSnowSoilNrg(:,iGRU)
    case(iLookINDEX%ixSnowOnlyNrg); parVector = indxStruct%ixSnowOnlyNrg(:,iGRU)
    case(iLookINDEX%ixSoilOnlyNrg); parVector = indxStruct%ixSoilOnlyNrg(:,iGRU)
    case(iLookINDEX%ixSnowSoilHyd); parVector = indxStruct%ixSnowSoilHyd(:,iGRU)
    case(iLookINDEX%ixSnowOnlyHyd); parVector = indxStruct%ixSnowOnlyHyd(:,iGRU)
    case(iLookINDEX%ixSoilOnlyHyd); parVector = indxStruct%ixSoilOnlyHyd(:,iGRU)
    case(iLookINDEX%ixNrgCanair); allocate(parVector(1)); parVector = indxStruct%ixNrgCanair(iGRU)
    case(iLookINDEX%ixNrgCanopy); allocate(parVector(1)); parVector = indxStruct%ixNrgCanopy(iGRU)
    case(iLookINDEX%ixHydCanopy); allocate(parVector(1)); parVector = indxStruct%ixHydCanopy(iGRU)
    case(iLookINDEX%ixNrgLayer); parVector = indxStruct%ixNrgLayer(:,iGRU)
    case(iLookINDEX%ixHydLayer); parVector = indxStruct%ixHydLayer(:,iGRU)
    case(iLookINDEX%ixWatAquifer); allocate(parVector(1)); parVector = indxStruct%ixWatAquifer(iGRU)
    ! case(iLookINDEX%ixVolFracWat); parVector = indxStruct%ixVolFracWat(:,iGRU)
    ! case(iLookINDEX%ixMatricHead); parVector = indxStruct%ixMatricHead(:,iGRU)
    case(iLookINDEX%ixAllState); parVector = indxStruct%ixAllState(:,iGRU)
    case(iLookINDEX%ixSoilState); parVector = indxStruct%ixSoilState(:,iGRU)
    case(iLookINDEX%ixLayerState); parVector = indxStruct%ixLayerState(:,iGRU)
    ! case(iLookINDEX%ixLayerActive); parVector = indxStruct%ixLayerActive(:,iGRU)
    case(iLookINDEX%numberFluxCalc); allocate(parVector(1)); parVector = indxStruct%numberFluxCalc
    case(iLookINDEX%numberStateSplit); allocate(parVector(1)); parVector = indxStruct%numberStateSplit
    case(iLookINDEX%numberDomainSplitNrg); allocate(parVector(1)); parVector = indxStruct%numberDomainSplitNrg
    case(iLookINDEX%numberDomainSplitMass); allocate(parVector(1)); parVector = indxStruct%numberDomainSplitMass
    case(iLookINDEX%numberScalarSolutions); allocate(parVector(1)); parVector = indxStruct%numberScalarSolutions
    end select
    if (.not.allocated(parVector)) then
      allocate(parVector(1)); parVector = -9999
    end if

  end subroutine

  subroutine get_device_diag_data(diagStruct, iVar, iGRU, parVector)
    implicit none
    type(diag_data_device) :: diagStruct
    integer(i4b) :: iGRU, iVar
    real(rkind),allocatable :: parVector(:)


    if (allocated(parVector)) deallocate(parVector)
    select case(iVar)
    ! local properties
    case(iLookDIAG%scalarCanopyDepth); allocate(parVector(1)); parVector = diagStruct%scalarCanopyDepth(iGRU)
    case(iLookDIAG%scalarBulkVolHeatCapVeg); allocate(parVector(1)); parVector = diagStruct%scalarBulkVolHeatCapVeg(iGRU)
    case(iLookDIAG%scalarCanopyCm); allocate(parVector(1)); parVector = diagStruct%scalarCanopyCm(iGRU)
    case(iLookDIAG%scalarCanopyEmissivity); allocate(parVector(1)); parVector = diagStruct%scalarCanopyEmissivity(iGRU)
    case(iLookDIAG%scalarRootZoneTemp); allocate(parVector(1)); parVector = diagStruct%scalarRootZoneTemp(iGRU)
    case(iLookDIAG%scalarLAI); allocate(parVector(1)); parVector = diagStruct%scalarLAI(iGRU)
    case(iLookDIAG%scalarSAI); allocate(parVector(1)); parVector = diagStruct%scalarSAI(iGRU)
    case(iLookDIAG%scalarExposedLAI); allocate(parVector(1)); parVector = diagStruct%scalarExposedLAI(iGRU)
    case(iLookDIAG%scalarExposedSAI); allocate(parVector(1)); parVector = diagStruct%scalarExposedSAI(iGRU)
    case(iLookDIAG%scalarAdjMeasHeight); allocate(parVector(1)); parVector = diagStruct%scalarAdjMeasHeight(iGRU)
    case(iLookDIAG%scalarCanopyIceMax); allocate(parVector(1)); parVector = diagStruct%scalarCanopyIceMax(iGRU)
    case(iLookDIAG%scalarCanopyLiqMax); allocate(parVector(1)); parVector = diagStruct%scalarCanopyLiqMax(iGRU)
    case(iLookDIAG%scalarGrowingSeasonIndex); allocate(parVector(1)); parVector = diagStruct%scalarGrowingSeasonIndex(iGRU)
    case(iLookDIAG%scalarVolHtCap_air); allocate(parVector(1)); parVector = diagStruct%scalarVolHtCap_air(iGRU)
    case(iLookDIAG%scalarVolHtCap_ice); allocate(parVector(1)); parVector = diagStruct%scalarVolHtCap_ice(iGRU)
    case(iLookDIAG%scalarVolHtCap_soil); allocate(parVector(1)); parVector = diagStruct%scalarVolHtCap_soil(iGRU)
    case(iLookDIAG%scalarVolHtCap_water); allocate(parVector(1)); parVector = diagStruct%scalarVolHtCap_water(iGRU)
    case(iLookDIAG%mLayerVolHtCapBulk); parVector = diagStruct%mLayerVolHtCapBulk(:,iGRU)
    case(iLookDIAG%mLayerCm); parVector = diagStruct%mLayerCm(:,iGRU)
    case(iLookDIAG%scalarLambda_drysoil); allocate(parVector(1)); parVector = diagStruct%scalarLambda_drysoil(iGRU)
    case(iLookDIAG%scalarLambda_wetsoil); allocate(parVector(1)); parVector = diagStruct%scalarLambda_wetsoil(iGRU)
    case(iLookDIAG%mLayerThermalC); parVector = diagStruct%mLayerThermalC(:,iGRU)
    case(iLookDIAG%iLayerThermalC); parVector = diagStruct%iLayerThermalC(:,iGRU)
    ! enthalpy
    case(iLookDIAG%scalarCanopyEnthTemp); allocate(parVector(1)); parVector = diagStruct%scalarCanopyEnthTemp(iGRU)
    case(iLookDIAG%mLayerEnthTemp); parVector = diagStruct%mLayerEnthTemp(:,iGRU)
    case(iLookDIAG%scalarTotalSoilEnthalpy); allocate(parVector(1)); parVector = diagStruct%scalarTotalSoilEnthalpy(iGRU)
    case(iLookDIAG%scalarTotalSnowEnthalpy); allocate(parVector(1)); parVector = diagStruct%scalarTotalSnowEnthalpy(iGRU)
    ! forcing
    case(iLookDIAG%scalarVPair); allocate(parVector(1)); parVector = diagStruct%scalarVPair(iGRU)
    case(iLookDIAG%scalarVP_CanopyAir); allocate(parVector(1)); parVector = diagStruct%scalarVP_CanopyAir(iGRU)
    case(iLookDIAG%scalarTwetbulb); allocate(parVector(1)); parVector = diagStruct%scalarTwetbulb(iGRU)
    case(iLookDIAG%scalarSnowfallTemp); allocate(parVector(1)); parVector = diagStruct%scalarSnowfallTemp(iGRU)
    case(iLookDIAG%scalarNewSnowDensity); allocate(parVector(1)); parVector = diagStruct%scalarNewSnowDensity(iGRU)
    case(iLookDIAG%scalarO2air); allocate(parVector(1)); parVector = diagStruct%scalarO2air(iGRU)
    case(iLookDIAG%scalarCO2air); allocate(parVector(1)); parVector = diagStruct%scalarCO2air(iGRU)
    case(iLookDIAG%windspd_x); allocate(parVector(1)); parVector = diagStruct%windspd_x(iGRU)
    case(iLookDIAG%windspd_y); allocate(parVector(1)); parVector = diagStruct%windspd_y(iGRU)
    ! shortwave radiation
    case(iLookDIAG%scalarCosZenith); allocate(parVector(1)); parVector = diagStruct%scalarCosZenith(iGRU)
    case(iLookDIAG%scalarFractionDirect); allocate(parVector(1)); parVector = diagStruct%scalarFractionDirect(iGRU)
    case(iLookDIAG%scalarCanopySunlitFraction); allocate(parVector(1)); parVector = diagStruct%scalarCanopySunlitFraction(iGRU)
    case(iLookDIAG%scalarCanopySunlitLAI); allocate(parVector(1)); parVector = diagStruct%scalarCanopySunlitLAI(iGRU)
    case(iLookDIAG%scalarCanopyShadedLAI); allocate(parVector(1)); parVector = diagStruct%scalarCanopyShadedLAI(iGRU)
    case(iLookDIAG%spectralAlbGndDirect); parVector = diagStruct%spectralAlbGndDirect(:,iGRU)
    case(iLookDIAG%spectralAlbGndDiffuse); parVector = diagStruct%spectralAlbGndDiffuse(:,iGRU)
    case(iLookDIAG%scalarGroundAlbedo); allocate(parVector(1)); parVector = diagStruct%scalarGroundAlbedo(iGRU)
    ! turbulent heat transfer
    case(iLookDIAG%scalarLatHeatSubVapCanopy); allocate(parVector(1)); parVector = diagStruct%scalarLatHeatSubVapCanopy(iGRU)
    case(iLookDIAG%scalarLatHeatSubVapGround); allocate(parVector(1)); parVector = diagStruct%scalarLatHeatSubVapGround(iGRU)
    case(iLookDIAG%scalarSatVP_CanopyTemp); allocate(parVector(1)); parVector = diagStruct%scalarSatVP_CanopyTemp(iGRU)
    case(iLookDIAG%scalarSatVP_GroundTemp); allocate(parVector(1)); parVector = diagStruct%scalarSatVP_GroundTemp(iGRU)
    case(iLookDIAG%scalarZ0Canopy); allocate(parVector(1)); parVector = diagStruct%scalarZ0Canopy(iGRU)
    case(iLookDIAG%scalarWindReductionFactor); allocate(parVector(1)); parVector = diagStruct%scalarWindReductionFactor(iGRU)
    case(iLookDIAG%scalarZeroPlaneDisplacement); allocate(parVector(1)); parVector = diagStruct%scalarZeroPlaneDisplacement(iGRU)
    case(iLookDIAG%scalarRiBulkCanopy); allocate(parVector(1)); parVector = diagStruct%scalarRiBulkCanopy(iGRU)
    case(iLookDIAG%scalarRiBulkGround); allocate(parVector(1)); parVector = diagStruct%scalarRiBulkGround(iGRU)
    case(iLookDIAG%scalarCanopyStabilityCorrection); allocate(parVector(1)); parVector = diagStruct%scalarCanopyStabilityCorrection(iGRU)
    case(iLookDIAG%scalarGroundStabilityCorrection); allocate(parVector(1)); parVector = diagStruct%scalarGroundStabilityCorrection(iGRU)
  ! evapotranspiration
    case(iLookDIAG%scalarIntercellularCO2Sunlit); allocate(parVector(1)); parVector = diagStruct%scalarIntercellularCO2Sunlit(iGRU)
    case(iLookDIAG%scalarIntercellularCO2Shaded); allocate(parVector(1)); parVector = diagStruct%scalarIntercellularCO2Shaded(iGRU)
    case(iLookDIAG%scalarTranspireLim); allocate(parVector(1)); parVector = diagStruct%scalarTranspireLim(iGRU)
    case(iLookDIAG%scalarTranspireLimAqfr); allocate(parVector(1)); parVector = diagStruct%scalarTranspireLimAqfr(iGRU)
    case(iLookDIAG%scalarFoliageNitrogenFactor); allocate(parVector(1)); parVector = diagStruct%scalarFoliageNitrogenFactor(iGRU)
    case(iLookDIAG%scalarSoilRelHumidity); allocate(parVector(1)); parVector = diagStruct%scalarSoilRelHumidity(iGRU)
    case(iLookDIAG%mLayerTranspireLim); parVector = diagStruct%mLayerTranspireLim(:,iGRU)
    case(iLookDIAG%mLayerRootDensity); parVector = diagStruct%mLayerRootDensity(:,iGRU)
    case(iLookDIAG%scalarAquiferRootFrac); allocate(parVector(1)); parVector = diagStruct%scalarAquiferRootFrac(iGRU)
  ! canopy hydrology
    case(iLookDIAG%scalarFracLiqVeg); allocate(parVector(1)); parVector = diagStruct%scalarFracLiqVeg(iGRU)
    case(iLookDIAG%scalarCanopyWetFraction); allocate(parVector(1)); parVector = diagStruct%scalarCanopyWetFraction(iGRU)
  ! snow hydrology
    case(iLookDIAG%scalarSnowAge); allocate(parVector(1)); parVector = diagStruct%scalarSnowAge(iGRU)
    case(iLookDIAG%scalarGroundSnowFraction); allocate(parVector(1)); parVector = diagStruct%scalarGroundSnowFraction(iGRU)
    case(iLookDIAG%spectralSnowAlbedoDirect); parVector = diagStruct%spectralSnowAlbedoDirect(:,iGRU)
    case(iLookDIAG%mLayerFracLiqSnow); parVector = diagStruct%mLayerFracLiqSnow(:,iGRU)
    case(iLookDIAG%mLayerThetaResid); parVector = diagStruct%mLayerThetaResid(:,iGRU)
    case(iLookDIAG%mLayerPoreSpace); parVector = diagStruct%mLayerPoreSpace(:,iGRU)
    case(iLookDIAG%mLayerMeltFreeze); parVector = diagStruct%mLayerMeltFreeze(:,iGRU)
  ! soil hydrology
    case(iLookDIAG%scalarInfilArea); allocate(parVector(1)); parVector = diagStruct%scalarInfilArea(iGRU)
    case(iLookDIAG%scalarFrozenArea); allocate(parVector(1)); parVector = diagStruct%scalarFrozenArea(iGRU)
    case(iLookDIAG%scalarSoilControl); allocate(parVector(1)); parVector = diagStruct%scalarSoilControl(iGRU)
    case(iLookDIAG%mLayerVolFracAir); parVector = diagStruct%mLayerVolFracAir(:,iGRU)
    case(iLookDIAG%mLayerTcrit); parVector = diagStruct%mLayerTcrit(:,iGRU)
    case(iLookDIAG%mLayerCompress); parVector = diagStruct%mLayerCompress(:,iGRU)
    case(iLookDIAG%scalarSoilCompress); allocate(parVector(1)); parVector = diagStruct%scalarSoilCompress(iGRU)
    case(iLookDIAG%mLayerMatricHeadLiq); parVector = diagStruct%mLayerMatricHeadLiq(:,iGRU)
    case(iLookDIAG%scalarTotalSoilLiq); allocate(parVector(1)); parVector = diagStruct%scalarTotalSoilLiq(iGRU)
    case(iLookDIAG%scalarTotalSoilIce); allocate(parVector(1)); parVector = diagStruct%scalarTotalSoilIce(iGRU)
    case(iLookDIAG%scalarTotalSoilWat); allocate(parVector(1)); parVector = diagStruct%scalarTotalSoilWat(iGRU)
  ! variable shortcuts
    case(iLookDIAG%scalarVGn_m); parVector = diagStruct%scalarVGn_m(:,iGRU)
    case(iLookDIAG%scalarKappa); allocate(parVector(1)); parVector = diagStruct%scalarKappa(iGRU)
    case(iLookDIAG%scalarVolLatHt_fus); allocate(parVector(1)); parVector = diagStruct%scalarVolLatHt_fus(iGRU)
  ! timing information
    case(iLookDIAG%numFluxCalls); allocate(parVector(1)); parVector = diagStruct%numFluxCalls
    case(iLookDIAG%wallClockTime); allocate(parVector(1)); parVector = diagStruct%wallClockTime
    case(iLookDIAG%meanStepSize); allocate(parVector(1)); parVector = diagStruct%meanStepSize
  ! balances
    case(iLookDIAG%balanceCasNrg); allocate(parVector(1)); parVector = diagStruct%balanceCasNrg(iGRU)
    case(iLookDIAG%balanceVegNrg); allocate(parVector(1)); parVector = diagStruct%balanceVegNrg(iGRU)
    case(iLookDIAG%balanceLayerNrg); parVector = diagStruct%balanceLayerNrg(:,iGRU)
    case(iLookDIAG%balanceSnowNrg); allocate(parVector(1)); parVector = diagStruct%balanceSnowNrg(iGRU)
    case(iLookDIAG%balanceSoilNrg); allocate(parVector(1)); parVector = diagStruct%balanceSoilNrg(iGRU)
    case(iLookDIAG%balanceVegMass); allocate(parVector(1)); parVector = diagStruct%balanceVegMass(iGRU)
    case(iLookDIAG%balanceLayerMass); parVector = diagStruct%balanceLayerMass(:,iGRU)
    case(iLookDIAG%balanceSnowMass); allocate(parVector(1)); parVector = diagStruct%balanceSnowMass(iGRU)
    case(iLookDIAG%balanceSoilMass); allocate(parVector(1)); parVector = diagStruct%balanceSoilMass(iGRU)
    case(iLookDIAG%balanceAqMass); allocate(parVector(1)); parVector = diagStruct%balanceAqMass(iGRU)
  ! sundials integrator stats
    case(iLookDIAG%numSteps); allocate(parVector(1)); parVector = diagStruct%numSteps
    case(iLookDIAG%numResEvals); allocate(parVector(1)); parVector = diagStruct%numResEvals
    case(iLookDIAG%numLinSolvSetups); allocate(parVector(1)); parVector = diagStruct%numLinSolvSetups
    case(iLookDIAG%numErrTestFails); allocate(parVector(1)); parVector = diagStruct%numErrTestFails
    case(iLookDIAG%kLast); allocate(parVector(1)); parVector = diagStruct%kLast
    case(iLookDIAG%kCur); allocate(parVector(1)); parVector = diagStruct%kCur
    case(iLookDIAG%hInitUsed); allocate(parVector(1)); parVector = diagStruct%hInitUsed
    case(iLookDIAG%hLast); allocate(parVector(1)); parVector = diagStruct%hLast
    case(iLookDIAG%hCur); allocate(parVector(1)); parVector = diagStruct%hCur
    case(iLookDIAG%tCur); allocate(parVector(1)); parVector = diagStruct%tCur

    end select
    if (.not.allocated(parVector)) then
      allocate(parVector(1)); parVector = -9999_rkind
    end if
  end subroutine

    subroutine get_device_prog_data(progStruct, iVar, iGRU, parVector)
    implicit none
    type(prog_data_device) :: progStruct
    integer(i4b) :: iGRU, iVar
    real(rkind),allocatable :: parVector(:)

    if (allocated(parVector)) deallocate(parVector)
    select case(iVar)
        ! define variables for time stepping
    case(iLookPROG%dt_init); allocate(parVector(1)); parVector = progStruct%dt_init
  ! state variables for vegetation
    case(iLookPROG%scalarCanopyIce); allocate(parVector(1)); parVector = progStruct%scalarCanopyIce(iGRU)
    case(iLookPROG%scalarCanopyLiq); allocate(parVector(1)); parVector = progStruct%scalarCanopyLiq(iGRU)
    case(iLookPROG%scalarCanopyWat); allocate(parVector(1)); parVector = progStruct%scalarCanopyWat(iGRU)
    case(iLookPROG%scalarCanairTemp); allocate(parVector(1)); parVector = progStruct%scalarCanairTemp(iGRU)
    case(iLookPROG%scalarCanopyTemp); allocate(parVector(1)); parVector = progStruct%scalarCanopyTemp(iGRU)
  ! state variables for snow
    case(iLookPROG%spectralSnowAlbedoDiffuse); parVector = progStruct%spectralSnowAlbedoDiffuse(:,iGRU)
    case(iLookPROG%scalarSnowAlbedo); allocate(parVector(1)); parVector = progStruct%scalarSnowAlbedo(iGRU)
    case(iLookPROG%scalarSnowDepth); allocate(parVector(1)); parVector = progStruct%scalarSnowDepth(iGRU)
    case(iLookPROG%scalarSWE); allocate(parVector(1)); parVector = progStruct%scalarSWE(iGRU)
    case(iLookPROG%scalarSfcMeltPond); allocate(parVector(1)); parVector = progStruct%scalarSfcMeltPond(iGRU)
  ! define state variables for the snow+soil domain
    case(iLookPROG%mLayerTemp); parVector = progStruct%mLayerTemp(:,iGRU)
    case(iLookPROG%mLayerVolFracIce); parVector = progStruct%mLayerVolFracIce(:,iGRU)
    case(iLookPROG%mLayerVolFracLiq); parVector = progStruct%mLayerVolFracLiq(:,iGRU)
    case(iLookPROG%mLayerVolFracWat); parVector = progStruct%mLayerVolFracWat(:,iGRU)
    case(iLookPROG%mLayerMatricHead); parVector = progStruct%mLayerMatricHead(:,iGRU)
  ! enthalpy
    case(iLookPROG%scalarCanairEnthalpy); allocate(parVector(1)); parVector = progStruct%scalarCanairEnthalpy(iGRU)
    case(iLookPROG%scalarCanopyEnthalpy); allocate(parVector(1)); parVector = progStruct%scalarCanopyEnthalpy(iGRU)
    case(iLookPROG%mLayerEnthalpy); parVector = progStruct%mLayerEnthalpy(:,iGRU)
  ! other state variables
    case(iLookPROG%scalarAquiferStorage); allocate(parVector(1)); parVector = progStruct%scalarAquiferStorage(iGRU)
    case(iLookPROG%scalarSurfaceTemp); allocate(parVector(1)); parVector = progStruct%scalarSurfaceTemp(iGRU)
  ! define coordinate variables
    case(iLookPROG%mLayerDepth); parVector = progStruct%mLayerDepth(:,iGRU)
    case(iLookPROG%mLayerHeight); parVector = progStruct%mLayerHeight(:,iGRU)
    case(iLookPROG%iLayerHeight); parVector = progStruct%iLayerHeight(:,iGRU)

    end select
    if (.not.allocated(parVector)) then
      allocate(parVector(1)); parVector = -9999_rkind
    end if
  end subroutine



    subroutine get_device_bvar_data(bvarStruct, iVar, iGRU, parVector)
    implicit none
    type(bvar_data_device) :: bvarStruct
    integer(i4b) :: iGRU, iVar
    real(rkind),allocatable :: parVector(:)

    if (allocated(parVector)) deallocate(parVector)
    select case(iVar)
    case(iLookBVAR%basin__totalArea); allocate(parVector(1)); parVector = bvarStruct%basin__totalArea(iGRU)
    case(iLookBVAR%basin__SurfaceRunoff); allocate(parVector(1)); parVector = bvarStruct%basin__SurfaceRunoff(iGRU)
    case(iLookBVAR%basin__ColumnOutflow); allocate(parVector(1)); parVector = bvarStruct%basin__ColumnOutflow(iGRU)
    case(iLookBVAR%basin__AquiferStorage); allocate(parVector(1)); parVector = bvarStruct%basin__AquiferStorage(iGRU)
    case(iLookBVAR%basin__AquiferRecharge); allocate(parVector(1)); parVector = bvarStruct%basin__AquiferRecharge(iGRU)
    case(iLookBVAR%basin__AquiferBaseflow); allocate(parVector(1)); parVector = bvarStruct%basin__AquiferBaseflow(iGRU)
    case(iLookBVAR%basin__AquiferTranspire); allocate(parVector(1)); parVector = bvarStruct%basin__AquiferTranspire(iGRU)
    case(iLookBVAR%basin__TotalRunoff); allocate(parVector(1)); parVector = bvarStruct%basin__TotalRunoff(iGRU)
    case(iLookBVAR%basin__SoilDrainage); allocate(parVector(1)); parVector = bvarStruct%basin__SoilDrainage(iGRU)
    case(iLookBVAR%routingRunoffFuture); parVector = bvarStruct%routingRunoffFuture(:,iGRU)
    case(iLookBVAR%routingFractionFuture); parVector = bvarStruct%routingFractionFuture(:,iGRU)
    case(iLookBVAR%averageInstantRunoff); allocate(parVector(1)); parVector = bvarStruct%averageInstantRunoff(iGRU)
    case(iLookBVAR%averageRoutedRunoff); allocate(parVector(1)); parVector = bvarStruct%averageRoutedRunoff(iGRU)
    end select
    if (.not.allocated(parVector)) then
      allocate(parVector(1)); parVector = -9999_rkind
    end if
  end subroutine
    subroutine get_device_bvar_stat(bvarStruct, iVar, nGRU,iFreq, parVector)
    implicit none
    type(bvar_stat_device) :: bvarStruct
    integer(i4b) :: nGRU, iVar,iFreq
    real(rkind),allocatable :: parVector(:)

    if (allocated(parVector)) deallocate(parVector)
    select case(iVar)
    case(iLookBVAR%basin__totalArea); allocate(parVector(nGRU)); parVector = bvarStruct%stat_basin__totalArea(iFreq,:)
    case(iLookBVAR%basin__SurfaceRunoff); allocate(parVector(nGRU)); parVector = bvarStruct%stat_basin__SurfaceRunoff(iFreq,:)
    case(iLookBVAR%basin__ColumnOutflow); allocate(parVector(nGRU)); parVector = bvarStruct%stat_basin__ColumnOutflow(iFreq,:)
    case(iLookBVAR%basin__AquiferStorage); allocate(parVector(nGRU)); parVector = bvarStruct%stat_basin__AquiferStorage(iFreq,:)
    case(iLookBVAR%basin__AquiferRecharge); allocate(parVector(nGRU)); parVector = bvarStruct%stat_basin__AquiferRecharge(iFreq,:)
    case(iLookBVAR%basin__AquiferBaseflow); allocate(parVector(nGRU)); parVector = bvarStruct%stat_basin__AquiferBaseflow(iFreq,:)
    case(iLookBVAR%basin__AquiferTranspire); allocate(parVector(nGRU)); parVector = bvarStruct%stat_basin__AquiferTranspire(iFreq,:)
    case(iLookBVAR%basin__TotalRunoff); allocate(parVector(nGRU)); parVector = bvarStruct%stat_basin__TotalRunoff(iFreq,:)
    case(iLookBVAR%basin__SoilDrainage); allocate(parVector(nGRU)); parVector = bvarStruct%stat_basin__SoilDrainage(iFreq,:)
    case(iLookBVAR%averageInstantRunoff); allocate(parVector(nGRU)); parVector = bvarStruct%stat_averageInstantRunoff(iFreq,:)
    case(iLookBVAR%averageRoutedRunoff); allocate(parVector(nGRU)); parVector = bvarStruct%stat_averageRoutedRunoff(iFreq,:)
    end select
    if (.not.allocated(parVector)) then
      allocate(parVector(nGRU)); parVector = -9999_rkind
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



  subroutine allocate_device_bvar_stat(bvarStat_d,nGRU)
    USE globalData,only:statBvar_meta             ! child metadata for stats
    USE var_lookup,only:maxvarFreq             ! allocation dimension (output frequency)
    implicit none
    type(bvar_stat_device) :: bvarStat_d
    ! type(gru_doubleVec) :: bvarStat

    integer(i4b) :: nGRU, iGRU, iVar, pVar
    integer :: err

    err = cudaDeviceSynchronize()
    allocate(bvarStat_d%stat_basin__SurfaceRunoff(maxvarfreq*2,nGRU))
     allocate(bvarStat_d%stat_basin__ColumnOutflow(maxvarFreq*2,nGRU))
 allocate(bvarStat_d%stat_basin__AquiferStorage(maxvarFreq*2,nGRU))
 allocate(bvarStat_d%stat_basin__AquiferRecharge(maxvarFreq*2,nGRU))
 allocate(bvarStat_d%stat_basin__AquiferBaseflow(maxvarFreq*2,nGRU))
 allocate(bvarStat_d%stat_basin__AquiferTranspire(maxvarFreq*2,nGRU))
 allocate(bvarStat_d%stat_basin__TotalRunoff(maxvarFreq*2,nGRU))
 allocate(bvarStat_d%stat_basin__SoilDrainage(maxvarFreq*2,nGRU))
 allocate(bvarStat_d%stat_basin__totalArea(maxvarFreq*2,nGRU))
 allocate(bvarStat_d%stat_averageInstantRunoff(maxvarFreq*2,nGRU))
 allocate(bvarStat_d%stat_averageRoutedRunoff(maxvarFreq*2,nGRU))
 do iVar=1,size(statBvar_meta)
  pVar = statBvar_meta(iVar)%ixParent
  select case(pVar)
    case(iLookBVAR%basin__SurfaceRunoff); bvarStat_d%statIndex_basin__SurfaceRunoff = statBvar_meta(iVar)%statIndex
    case(iLookBVAR%basin__ColumnOutflow); bvarStat_d%statIndex_basin__ColumnOutflow = statBvar_meta(iVar)%statIndex
    case(iLookBVAR%basin__AquiferStorage); bvarStat_d%statIndex_basin__AquiferStorage = statBvar_meta(iVar)%statIndex
    case(iLookBVAR%basin__AquiferRecharge); bvarStat_d%statIndex_basin__AquiferRecharge = statBvar_meta(iVar)%statIndex
    case(iLookBVAR%basin__AquiferBaseflow); bvarStat_d%statIndex_basin__AquiferBaseflow = statBvar_meta(iVar)%statIndex
    case(iLookBVAR%basin__AquiferTranspire); bvarStat_d%statIndex_basin__AquiferTranspire = statBvar_meta(iVar)%statIndex
    case(iLookBVAR%basin__TotalRunoff); bvarStat_d%statIndex_basin__TotalRunoff = statBvar_meta(iVar)%statIndex
    case(iLookBVAR%basin__SoilDrainage); bvarStat_d%statIndex_basin__SoilDrainage = statBvar_meta(iVar)%statIndex
    case(iLookBVAR%basin__totalArea); bvarStat_d%statIndex_basin__totalArea = statBvar_meta(iVar)%statIndex
    case(iLookBVAR%averageInstantRunoff); bvarStat_d%statIndex_averageInstantRunoff = statBvar_meta(iVar)%statIndex
    case(iLookBVAR%averageRoutedRunoff); bvarStat_d%statIndex_averageRoutedRunoff = statBvar_meta(iVar)%statIndex
    end select
  end do
  if (.not.allocated(bvarStat_d%statIndex_basin__SurfaceRunoff)) allocate(bvarStat_d%statIndex_basin__SurfaceRunoff(1))
if (.not.allocated(bvarStat_d%statIndex_basin__ColumnOutflow)) allocate(bvarStat_d%statIndex_basin__ColumnOutflow(1))
if (.not.allocated(bvarStat_d%statIndex_basin__AquiferStorage)) allocate(bvarStat_d%statIndex_basin__AquiferStorage(1))
if (.not.allocated(bvarStat_d%statIndex_basin__AquiferRecharge)) allocate(bvarStat_d%statIndex_basin__AquiferRecharge(1))
if (.not.allocated(bvarStat_d%statIndex_basin__AquiferBaseflow)) allocate(bvarStat_d%statIndex_basin__AquiferBaseflow(1))
if (.not.allocated(bvarStat_d%statIndex_basin__AquiferTranspire)) allocate(bvarStat_d%statIndex_basin__AquiferTranspire(1))
if (.not.allocated(bvarStat_d%statIndex_basin__TotalRunoff)) allocate(bvarStat_d%statIndex_basin__TotalRunoff(1))
if (.not.allocated(bvarStat_d%statIndex_basin__SoilDrainage)) allocate(bvarStat_d%statIndex_basin__SoilDrainage(1))
if (.not.allocated(bvarStat_d%statIndex_basin__totalArea)) allocate(bvarStat_d%statIndex_basin__totalArea(1))
if (.not.allocated(bvarStat_d%statIndex_averageInstantRunoff)) allocate(bvarStat_d%statIndex_averageInstantRunoff(1))
if (.not.allocated(bvarStat_d%statIndex_averageRoutedRunoff)) allocate(bvarStat_d%statIndex_averageRoutedRunoff(1))

  end subroutine

  subroutine finalize_device_bvar_stat(bvarStat_d, bvarStat)
    USE globalData,only:statBvar_meta             ! child metadata for stats

    type(bvar_stat_device) :: bvarStat_d
    type(gru_doubleVec) :: bvarStat

    integer(i4b) :: nGRU, iGRU, iVar, pVar
    nGRU = size(bvarStat%gru)

    do iGRU=1,nGRU
      do iVar=1,size(bvarStat%gru(iGRU)%var)
    pVar = statBvar_meta(iVar)%ixParent
    if (allocated(bvarStat%gru(iGRU)%var(iVar)%dat)) then
      select case (pVar)
      case(iLookBVAR%basin__SurfaceRunoff); bvarStat%gru(iGRU)%var(iVar)%dat = bvarStat_d%stat_basin__SurfaceRunoff(:,iGRU)
      case(iLookBVAR%basin__ColumnOutflow); bvarStat%gru(iGRU)%var(iVar)%dat = bvarStat_d%stat_basin__ColumnOutflow(:,iGRU)
      case(iLookBVAR%basin__AquiferStorage); bvarStat%gru(iGRU)%var(iVar)%dat = bvarStat_d%stat_basin__AquiferStorage(:,iGRU)
      case(iLookBVAR%basin__AquiferRecharge); bvarStat%gru(iGRU)%var(iVar)%dat = bvarStat_d%stat_basin__AquiferRecharge(:,iGRU)
      case(iLookBVAR%basin__AquiferBaseflow); bvarStat%gru(iGRU)%var(iVar)%dat = bvarStat_d%stat_basin__AquiferBaseflow(:,iGRU)
      case(iLookBVAR%basin__AquiferTranspire); bvarStat%gru(iGRU)%var(iVar)%dat = bvarStat_d%stat_basin__AquiferTranspire(:,iGRU)
      case(iLookBVAR%basin__TotalRunoff); bvarStat%gru(iGRU)%var(iVar)%dat = bvarStat_d%stat_basin__TotalRunoff(:,iGRU)
      case(iLookBVAR%basin__SoilDrainage); bvarStat%gru(iGRU)%var(iVar)%dat = bvarStat_d%stat_basin__SoilDrainage(:,iGRU)
      case(iLookBVAR%basin__totalArea); bvarStat%gru(iGRU)%var(iVar)%dat = bvarStat_d%stat_basin__totalArea(:,iGRU)
      case(iLookBVAR%averageInstantRunoff); bvarStat%gru(iGRU)%var(iVar)%dat = bvarStat_d%stat_averageInstantRunoff(:,iGRU)
      case(iLookBVAR%averageRoutedRunoff); bvarStat%gru(iGRU)%var(iVar)%dat = bvarStat_d%stat_averageRoutedRunoff(:,iGRU)
      end select
    end if
    end do
 end do
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


  subroutine allocate_device_indx_data_snow(indx_data_device, nSoil,nGRU)
    use globalData,only:maxSnowLayers,integerMissing
    implicit none
    type(indx_data_device), intent(inout) :: indx_data_device
    integer(i4b), intent(in) :: nSoil,nGRU
      
    integer(i4b) :: nSnow, nLayers,nState,maxSize
    nSnow = maxSnowLayers+1
    nLayers = nSoil + maxSnowLayers+1
    nState = 4+2*nLayers
  
    allocate(indx_data_device%ixNrgLayer(nLayers,nGRU))
    allocate(indx_data_device%ixHydLayer(nLayers,nGRU))

    allocate(indx_data_device%ixSnowSoilNrg(nLayers,nGRU))
    allocate(indx_data_device%ixSnowSoilHyd(nLayers,nGRU))
    allocate(indx_data_device%ixHydType(nLayers,nGRU))
    allocate(indx_data_device%ixStateType(nState,nGRU))
    allocate(indx_data_device%ixControlVolume(nState,nGRU))
    allocate(indx_data_device%layerType(nLayers,nGRU))
    allocate(indx_data_device%ixSnowOnlyHyd(nSnow,nGRU))
    allocate(indx_data_device%ixSnowOnlyNrg(nSnow,nGRU))
        

    maxSize = 0
    if (maxSize.ne.0) allocate(indx_data_device%ixNrgOnly(maxSize,nGRU))
    allocate(indx_data_device%ixDomainType(nState,nGRU))
      indx_data_device%ixStateType = integerMissing
      indx_data_device%ixControlVolume = integerMissing
      indx_data_device%layerType = integerMissing

      indx_data_device%ixHydType = integerMissing
      indx_data_device%ixNrgLayer = integerMissing
      indx_data_device%ixHydLayer = integerMissing
      indx_data_device%ixSnowSoilHyd = integerMissing
      indx_data_device%ixSnowSoilNrg = integerMissing
      indx_data_device%ixSnowOnlyHyd = integerMissing
      indx_data_device%ixSnowOnlyNrg = integerMissing
      indx_data_device%ixDomainType = integerMissing
      if (allocated(indx_data_device%ixNrgOnly)) indx_data_device%ixNrgOnly = integerMissing
  
  end subroutine allocate_device_indx_data_snow

    subroutine allocate_device_indx_data(indx_data_device, nSoil,nGRU)
    use globalData,only:maxSnowLayers,integerMissing
    implicit none
    type(indx_data_device), intent(inout) :: indx_data_device
    integer(i4b), intent(in) :: nSoil,nGRU
          allocate(indx_data_device%nSnow(nGRU))
              allocate(indx_data_device%ixSoilState(nSoil,nGRU))

    indx_data_device%nSnow = integerMissing
    indx_data_device%nSoil = integerMissing
    allocate(indx_data_device%nLayers(nGRU))
    indx_data_device%nLayers = integerMissing
  allocate(indx_data_device%nState(nGRU))
  indx_data_device%nState = integerMissing
  indx_data_device%numberDomainSplitMass = integerMissing
indx_data_device%numberDomainSplitNrg = integerMissing
indx_data_device%numberScalarSolutions = integerMissing
indx_data_device%numberStateSplit = integerMissing
    indx_data_device%numberFluxCalc = integerMissing

    allocate(indx_data_device%ixCasNrg(nGRU))
    indx_data_device%ixCasNrg = integerMissing
    allocate(indx_data_device%ixVegNrg(nGRU))
    indx_data_device%ixVegNrg = integerMissing
    allocate(indx_data_device%ixVegHyd(nGRU))
    indx_data_device%ixVegHyd = integerMissing
    allocate(indx_data_device%ixAqWat(nGRU))
    indx_data_device%ixAqWat = integerMissing
    allocate(indx_data_device%ixTopNrg(nGRU))
    indx_data_device%ixTopNrg = integerMissing
    allocate(indx_data_device%ixTopHyd(nGRU))
    indx_data_device%ixTopHyd = integerMissing

    allocate(indx_data_device%ixHydCanopy(nGRU))
    allocate(indx_data_device%ixNrgCanopy(nGRU))
    allocate(indx_data_device%ixNrgCanair(nGRU))
        allocate(indx_data_device%ixSoilOnlyHyd(nSoil,nGRU))
    allocate(indx_data_device%ixSoilOnlyNrg(nSoil,nGRU))

          indx_data_device%ixHydCanopy = integerMissing
          indx_data_device%ixNrgCanopy = integerMissing

    indx_data_device%ixNrgCanair = integerMissing
        indx_data_device%ixSoilOnlyHyd = integerMissing
        indx_data_device%ixSoilOnlyNrg = integerMissing
            indx_data_device%ixSoilState = integerMissing


    allocate(indx_data_device%nCasNrg(nGRU))
    indx_data_device%nCasNrg = integerMissing
    allocate(indx_data_device%nVegNrg(nGRU))
    indx_data_device%nVegNrg = integerMissing
    allocate(indx_data_device%nVegMass(nGRU))
    indx_data_device%nVegMass = integerMissing
    allocate(indx_data_device%nVegState(nGRU))
    indx_data_device%nVegState = integerMissing
    allocate(indx_data_device%nNrgState(nGRU))
    indx_data_device%nNrgState = integerMissing
    allocate(indx_data_device%nWatState(nGRU))
    indx_data_device%nWatState = integerMissing
    allocate(indx_data_device%nMatState(nGRU))
    indx_data_device%nMatState = integerMissing
    allocate(indx_data_device%nMassState(nGRU))
    indx_data_device%nMassState = integerMissing
    allocate(indx_data_device%ixWatAquifer(nGRU))
    indx_data_device%ixWatAquifer = integerMissing

  
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
    allocate(indx_data_device%nLayers(nGRU))
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

  subroutine finalize_device_indx_data(indx_data_device, indx_data)
    implicit none
    type(indx_data_device), intent(in) :: indx_data_device
    type(gru_hru_intVec),intent(inout) :: indx_data              ! indices defining model states and layers

    integer(i4b) :: nLayers,nSnow,nState,iGRU,nGRU
    nGRU = size(indx_data%gru)

    do iGRU=1,nGRU
        nLayers = indx_data_device%nLayers(iGRU)
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
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nLayers)%dat(1) = indx_data_device%nLayers(iGRU)
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
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnowSoilNrg)%dat(1) = indx_data_device%nLayers(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnowSoilHyd)%dat(1) = indx_data_device%nLayers(iGRU)
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
      

  subroutine allocate_device_prog_data(prog_data_device, nGRU, nSoil)
    use globalData,only:realMissing
    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    integer(i4b),intent(in) :: nGRU,nSoil

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
    allocate(prog_data_device%mLayerMatricHead(nSoil,nGRU))              ! get_ixVarType('midSoil')
    ! enthalpy
    allocate(prog_data_device%scalarCanairEnthalpy(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyEnthalpy(nGRU))            ! get_ixVarType('scalarv')
    ! other state variables
    allocate(prog_data_device%scalarAquiferStorage(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSurfaceTemp(nGRU))               ! get_ixVarType('scalarv')
      prog_data_device%spectralSnowAlbedoDiffuse = realMissing
      prog_data_device%mLayerMatricHead = realMissing
    prog_data_device%scalarCanopyIce = realMissing
    prog_data_device%scalarCanopyLiq = realMissing
    prog_data_device%scalarCanopyWat = realMissing
    prog_data_device%scalarCanairTemp = realMissing
    prog_data_device%scalarCanopyTemp = realMissing
    prog_data_device%scalarSnowAlbedo = realMissing
    prog_data_device%scalarSnowDepth = realMissing
    prog_data_device%scalarSWE = realMissing
    prog_data_device%scalarSfcMeltPond = realMissing
    prog_data_device%scalarCanairEnthalpy = realMissing
    prog_data_device%scalarCanopyEnthalpy = realMissing
    prog_data_device%scalarAquiferStorage = realMissing
    prog_data_device%scalarSurfaceTemp = realMissing

    prog_data_device%dt_init = realMissing

  end subroutine allocate_device_prog_data  

    subroutine allocate_device_prog_data_snow(prog_data_device, nGRU, nSoil)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers
    USE globalData,only:realMissing            ! missing real number
    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    integer(i4b),intent(in) :: nGRU,nSoil

    integer(i4b) :: nLayers
    nLayers = nSoil + maxSnowLayers+1
    allocate(prog_data_device%mLayerTemp(nLayers,nGRU))
    allocate(prog_data_device%mLayerVolFracIce(nLayers,nGRU))
    allocate(prog_data_device%mLayerVolFracLiq(nLayers,nGRU))
    allocate(prog_data_device%mLayerVolFracWat(nLayers,nGRU))
    allocate(prog_data_device%mLayerEnthalpy(nLayers,nGRU))
    allocate(prog_data_device%mLayerDepth(nLayers,nGRU))
    allocate(prog_data_device%mLayerHeight(0:nLayers,nGRU))
    allocate(prog_data_device%iLayerHeight(0:nLayers,nGRU))

    prog_data_device%mLayerTemp = realMissing
    prog_data_device%mLayerVolFracIce = realMissing
    prog_data_device%mLayerVolFracLiq = realMissing
    prog_data_device%mLayerVolFracWat = realMissing
    prog_data_device%mLayerEnthalpy = realMissing
    prog_data_device%mLayerDepth = realMissing
    prog_data_device%mLayerHeight = realMissing
    prog_data_device%iLayerHeight = realMissing

  end subroutine allocate_device_prog_data_snow  


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
      

  subroutine allocate_device_flux_data(flux_data_device,nSoil,nGRU)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers
    use globalData,only:realMissing

    type(flux_data_device), intent(inout) :: flux_data_device
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
    flux_data_device%data = realMissing
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

  subroutine allocate_device_diag_data_snow(diag_data_device, nGRU, nSoil)
    use globalData,only:maxSnowLayers,realMissing
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_device
    integer(i4b),intent(in) :: nGRU,nSoil

    integer(i4b) :: nLayers, nSnow
    integer(i4b) :: iGRU

    nLayers = nSoil + maxSnowLayers+1
    nSnow = maxSnowLayers+1

    allocate(diag_data_device%mLayerVolHtCapBulk(nLayers,nGRU))
    allocate(diag_data_device%mLayerCm(nLayers,nGRU))
    allocate(diag_data_device%mLayerThermalC(nLayers,nGRU))
    allocate(diag_data_device%iLayerThermalC(0:nLayers,nGRU))
    allocate(diag_data_device%mLayerEnthTemp(nLayers,nGRU))
    allocate(diag_data_device%mLayerFracLiqSnow(nSnow,nGRU))
    allocate(diag_data_device%mLayerThetaResid(nSnow,nGRU))
    allocate(diag_data_device%mLayerPoreSpace(nSnow,nGRU))
    allocate(diag_data_device%mLayerMeltFreeze(nLayers,nGRU))
    allocate(diag_data_device%mLayerVolFracAir(nLayers,nGRU))
    allocate(diag_data_device%balanceLayerNrg(nLayers,nGRU))
    allocate(diag_data_device%balanceLayerMass(nLayers,nGRU))

      diag_data_device%mLayerVolHtCapBulk = realMissing
      diag_data_device%mLayerCm = realMissing
      diag_data_device%mLayerThermalC = realMissing
      diag_data_device%iLayerThermalC = realMissing
      diag_data_device%mLayerEnthTemp = realMissing
  diag_data_device%mLayerFracLiqSnow = realMissing
  diag_data_device%mLayerThetaResid = realMissing
  diag_data_device%mLayerPoreSpace = realMissing
  diag_data_device%mLayerMeltFreeze = realMissing

    diag_data_device%mLayerVolFracAir = realMissing
      diag_data_device%balanceLayerNrg = realMissing
      diag_data_device%balanceLayerMass = realMissing

  end subroutine allocate_device_diag_data_snow


    subroutine allocate_device_diag_data(diag_data_device, nGRU, nSoil)
      use globalData,only:realMissing
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_device
    integer(i4b),intent(in) :: nGRU,nSoil



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
    allocate(diag_data_device%scalarLambda_drysoil(nGRU))
    allocate(diag_data_device%scalarLambda_wetsoil(nGRU))
    allocate(diag_data_device%scalarCanopyEnthTemp(nGRU))
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
    allocate(diag_data_device%scalarInfilArea(nGRU))
    allocate(diag_data_device%scalarFrozenArea(nGRU))
    allocate(diag_data_device%scalarSoilControl(nGRU))
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
    allocate(diag_data_device%balanceSnowNrg(nGRU))
    allocate(diag_data_device%balanceSoilNrg(nGRU))
    allocate(diag_data_device%balanceVegMass(nGRU))
    allocate(diag_data_device%balanceSnowMass(nGRU))
    allocate(diag_data_device%balanceSoilMass(nGRU))
    allocate(diag_data_device%balanceAqMass(nGRU))

    diag_data_device%scalarCanopyDepth = realMissing
    diag_data_device%scalarBulkVolHeatCapVeg = realMissing
    diag_data_device%scalarCanopyCm = realMissing
    diag_data_device%scalarCanopyEmissivity = realMissing
    diag_data_device%scalarRootZoneTemp = realMissing
    diag_data_device%scalarLAI = realMissing
    diag_data_device%scalarSAI = realMissing
    diag_data_device%scalarExposedLAI = realMissing
    diag_data_device%scalarExposedSAI = realMissing
    diag_data_device%scalarAdjMeasHeight = realMissing
    diag_data_device%scalarCanopyIceMax = realMissing
    diag_data_device%scalarCanopyLiqMax = realMissing
    diag_data_device%scalarGrowingSeasonIndex = realMissing
    diag_data_device%scalarVolHtCap_air = realMissing
    diag_data_device%scalarVolHtCap_ice = realMissing
    diag_data_device%scalarVolHtCap_soil = realMissing
    diag_data_device%scalarVolHtCap_water = realMissing
    diag_data_device%scalarLambda_drysoil = realMissing
    diag_data_device%scalarLambda_wetsoil = realMissing
    diag_data_device%scalarCanopyEnthTemp = realMissing
    diag_data_device%scalarTotalSoilEnthalpy = realMissing
    diag_data_device%scalarTotalSnowEnthalpy = realMissing
    diag_data_device%scalarVPair = realMissing
    diag_data_device%scalarVP_CanopyAir = realMissing
    diag_data_device%scalarTwetbulb = realMissing
    diag_data_device%scalarSnowfallTemp = realMissing
    diag_data_device%scalarNewSnowDensity = realMissing
    diag_data_device%scalarO2air = realMissing
    diag_data_device%scalarCO2air = realMissing
    diag_data_device%windspd_x = realMissing
    diag_data_device%windspd_y = realMissing
    diag_data_device%scalarCosZenith = realMissing
    diag_data_device%scalarFractionDirect = realMissing
    diag_data_device%scalarCanopySunlitFraction = realMissing
    diag_data_device%scalarCanopySunlitLAI = realMissing
    diag_data_device%scalarCanopyShadedLAI = realMissing
      diag_data_device%spectralAlbGndDirect = realMissing
      diag_data_device%spectralAlbGndDiffuse = realMissing
    diag_data_device%scalarGroundAlbedo = realMissing
    diag_data_device%scalarLatHeatSubVapCanopy = realMissing
    diag_data_device%scalarLatHeatSubVapGround = realMissing
    diag_data_device%scalarSatVP_CanopyTemp = realMissing
    diag_data_device%scalarSatVP_GroundTemp = realMissing
    diag_data_device%scalarZ0Canopy = realMissing
    diag_data_device%scalarWindReductionFactor = realMissing
    diag_data_device%scalarZeroPlaneDisplacement = realMissing
    diag_data_device%scalarRiBulkCanopy = realMissing
    diag_data_device%scalarRiBulkGround = realMissing
    diag_data_device%scalarCanopyStabilityCorrection = realMissing
    diag_data_device%scalarGroundStabilityCorrection = realMissing

    diag_data_device%scalarIntercellularCO2Sunlit = realMissing
    diag_data_device%scalarIntercellularCO2Shaded = realMissing
    diag_data_device%scalarTranspireLim = realMissing
    diag_data_device%scalarTranspireLimAqfr = realMissing
    diag_data_device%scalarFoliageNitrogenFactor = realMissing
    diag_data_device%scalarSoilRelHumidity = realMissing
      diag_data_device%mLayerTranspireLim = realMissing
      diag_data_device%mLayerRootDensity = realMissing
    diag_data_device%scalarAquiferRootFrac = realMissing
    diag_data_device%scalarFracLiqVeg = realMissing
    diag_data_device%scalarCanopyWetFraction = realMissing

    diag_data_device%scalarSnowAge = realMissing
diag_data_device%scalarGroundSnowFraction = realMissing
  diag_data_device%spectralSnowAlbedoDirect = realMissing

  
  diag_data_device%scalarInfilArea = realMissing
  diag_data_device%scalarFrozenArea = realMissing
  diag_data_device%scalarSoilControl = realMissing
  diag_data_device%mLayerTcrit = realMissing
  diag_data_device%mLayerCompress = realMissing
  diag_data_device%scalarSoilCompress = realMissing
  diag_data_device%mLayerMatricHeadLiq = realMissing
  diag_data_device%scalarTotalSoilLiq = realMissing
  diag_data_device%scalarTotalSoilIce = realMissing
  diag_data_device%scalarTotalSoilWat = realMissing
      diag_data_device%scalarVGn_m = realMissing

    diag_data_device%scalarKappa = realMissing
    diag_data_device%scalarVolLatHt_fus = realMissing

    diag_data_device%numFluxCalls = realMissing
    diag_data_device%wallClockTime = realMissing
    diag_data_device%meanStepSize = realMissing


    diag_data_device%balanceCasNrg = realMissing
    diag_data_device%balanceVegNrg = realMissing
    diag_data_device%balanceSnowNrg = realMissing
    diag_data_device%balanceSoilNrg = realMissing
    diag_data_device%balanceVegMass = realMissing
    diag_data_device%balanceSnowMass = realMissing
    diag_data_device%balanceSoilMass = realMissing
    diag_data_device%balanceAqMass = realMissing

    diag_data_device%numSteps = realMissing
    diag_data_device%numResEvals = realMissing
    diag_data_device%numLinSolvSetups = realMissing
    diag_data_device%numErrTestFails = realMissing
    diag_data_device%kLast = realMissing
    diag_data_device%kCur = realMissing
    diag_data_device%hInitUsed = realMissing
    diag_data_device%hLast = realMissing
    diag_data_device%hCur = realMissing
    diag_data_device%tCur = realMissing
  end subroutine allocate_device_diag_data


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

  subroutine allocate_device_bvar_data(bvar_data_device,nGRU)
    use globalData,only:nTimeDelay
    USE globalData,only:bvar_meta                     ! metadata structures

    implicit none
    type(bvar_data_device),intent(inout) :: bvar_data_device
    integer(i4b),intent(in) :: nGRU
    integer(i4b) :: iGRU,iVar

    allocate(bvar_data_device%basin__AquiferStorage(nGRU))
    allocate(bvar_data_device%basin__SurfaceRunoff(nGRU))
allocate(bvar_data_device%basin__ColumnOutflow(nGRU))
allocate(bvar_data_device%basin__AquiferRecharge(nGRU))
allocate(bvar_data_device%basin__AquiferBaseflow(nGRU))
allocate(bvar_data_device%basin__AquiferTranspire(nGRU))
allocate(bvar_data_device%basin__TotalRunoff(nGRU))
allocate(bvar_data_device%basin__SoilDrainage(nGRU))
allocate(bvar_data_device%basin__totalArea(nGRU))
allocate(bvar_data_device%routingRunoffFuture(nTimeDelay,nGRU))
allocate(bvar_data_device%routingFractionFuture(nTimeDelay,nGRU))
allocate(bvar_data_device%averageInstantRunoff(nGRU))
allocate(bvar_data_device%averageRoutedRunoff(nGRU))

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
  ! bvar_data%gru(iGRU)%var(iLookBVAR%routingRunoffFuture)%dat = bvar_data_device%routingRunoffFuture(:,iGRU)
  ! bvar_data%gru(iGRU)%var(iLookBVAR%routingFractionFuture)%dat = bvar_data_device%routingFractionFuture(:,iGRU)
bvar_data%gru(iGRU)%var(iLookBVAR%averageInstantRunoff)%dat(1) = bvar_data_device%averageInstantRunoff(iGRU)
bvar_data%gru(iGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(1) = bvar_data_device%averageRoutedRunoff(iGRU)
    end do
  end subroutine finalize_device_bvar_data

  subroutine deallocate_device_bvar_data(bvar_data_device)
    type(bvar_data_device),intent(inout) :: bvar_data_device
        
    deallocate(bvar_data_device%basin__AquiferStorage)
  end subroutine deallocate_device_bvar_data

  subroutine allocate_device_forc_data(forc_data_device,nGRU)
    use globalData,only:realMissing
    type(forc_data_device),intent(inout) :: forc_data_device
    integer(i4b) :: nGRU

    allocate(forc_data_device%pptrate(nGRU))
    allocate(forc_data_device%SWRadAtm(nGRU))
    allocate(forc_data_device%LWRadAtm(nGRU))
    allocate(forc_data_device%airtemp(nGRU))
    allocate(forc_data_device%windspd(nGRU))
    allocate(forc_data_device%airpres(nGRU))
    allocate(forc_data_device%spechum(nGRU))

      forc_data_device%time_h = realMissing
      forc_data_device%pptrate = realMissing
      forc_data_device%SWRadAtm = realMissing
      forc_data_device%LWRadAtm = realMissing
      forc_data_device%airtemp = realMissing
      forc_data_device%windspd = realMissing
      forc_data_device%airpres = realMissing
      forc_data_device%spechum = realMissing

  end subroutine allocate_device_forc_data
          
  subroutine finalize_device_forc_data(forc_data_device, forc_data)
    type(forc_data_device),intent(inout) :: forc_data_device
    type(gru_hru_double),intent(in) :: forc_data
    integer(i4b) :: nGRU,iGRU

    nGRU = size(forc_data%gru)

    do iGRU=1,nGRU
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%time) = forc_data_device%time_h
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%pptrate) = forc_data_device%pptrate(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%SWRadAtm) = forc_data_device%SWRadAtm(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%LWRadAtm) = forc_data_device%LWRadAtm(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%airtemp) = forc_data_device%airtemp(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%windspd) = forc_data_device%windspd(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%airpres) = forc_data_device%airpres(iGRU)
      forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%spechum) = forc_data_device%spechum(iGRU)
    end do

  end subroutine finalize_device_forc_data

  subroutine deallocate_device_forc_data(forc_data_device)
    type(forc_data_device),intent(inout) :: forc_data_device

    deallocate(forc_data_device%pptrate)
    deallocate(forc_data_device%SWRadAtm)
    deallocate(forc_data_device%LWRadAtm)
    deallocate(forc_data_device%airtemp)
    deallocate(forc_data_device%windspd)
    deallocate(forc_data_device%airpres)
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
    use NOAHMP_RAD_PARAMETERS,only:MBAND,msc
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
    allocate(veg_param%rsmin_d(nGRU))
    allocate(veg_param%rsmax_d(nGRU))
    allocate(veg_param%topt_d(nGRU))
    allocate(veg_param%hs_d(nGRU))
    veg_param%rhol = rhol
    veg_param%rhos = rhos
    veg_param%taul = taul
    veg_param%taus = taus
    veg_param%opt_alb = opt_alb
    veg_param%omegas = (/0.8, 0.4/)
    veg_param%opt_rad = opt_rad
    veg_param%tfrz = tfrz
    veg_param%betais = 0.5
    veg_param%betads = 0.5
    veg_param%xl = xl
    allocate(veg_param%hvt_(size(hvt),nGRU))
    allocate(veg_param%hvb_(size(hvb),nGRU))

    do iGRU=1,nGRU
      veg_param%hvt_(:,iGRU) = hvt
      veg_param%hvb_(:,iGRU) = hvb
    end do
    veg_param%rc = rc
    veg_param%swemx = swemx
    allocate(veg_param%albsat(msc,mband))
    veg_param%albsat(:,1)=(/0.15,0.11,0.10,0.09,0.08,0.07,0.06,0.05/)
    veg_param%albsat(:,2)=(/0.30,0.22,0.20,0.18,0.16,0.14,0.12,0.10/)
    allocate(veg_param%albdry(msc,mband))
    veg_param%albdry(:,1) = (/0.27,0.22,0.20,0.18,0.16,0.14,0.12,0.10/)
    veg_param%albdry(:,2) = (/0.54,0.44,0.40,0.36,0.32,0.28,0.24,0.20/)
    veg_param%alblak=(/0.60, 0.40/)
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
      attr_data_device%aspect(iGRU) = attr_data%gru(iGRU)%hru(1)%var(iLookATTR%aspect)
    end do
  end subroutine

  subroutine get_device_attr_data(attrStruct, iVar, iGRU, parVector)
    implicit none
    type(attr_data_device) :: attrStruct
    integer(i4b) :: iGRU, iVar
    real(rkind),allocatable :: parVector(:)

    if (allocated(parVector)) deallocate(parVector)
    allocate(parVector(1))
    parVector = -9999_rkind
    select case(iVar)
      case(iLookATTR%latitude); parVector = attrStruct%latitude(iGRU)
      case(iLookATTR%longitude); parVector = attrStruct%longitude(iGRU)
      case(iLookATTR%elevation); parVector = attrStruct%elevation(iGRU)
      case(iLookATTR%tan_slope); parVector = attrStruct%tan_slope(iGRU)
      case(iLookATTR%contourLength); parVector = attrStruct%contourLength(iGRU)
      case(iLookATTR%HRUarea); parVector = attrStruct%HRUarea(iGRU)
      case(iLookATTR%mHeight); parVector = attrStruct%mHeight(iGRU)
      case(iLookATTR%aspect); parVector = attrStruct%aspect(iGRU)
    end select
    if (.not.allocated(parVector)) then
      allocate(parVector(1)); parVector = -9999_rkind
    end if
  end subroutine

  subroutine set_device_attr_data(attrStruct, iVar, iGRU, parVector)
    implicit none
    type(attr_data_device) :: attrStruct
    integer(i4b) :: iGRU, iVar
    real(rkind) :: parVector

    select case(iVar)
      case(iLookATTR%latitude); attrStruct%latitude(iGRU) = parVector
      case(iLookATTR%longitude); attrStruct%longitude(iGRU) = parVector
      case(iLookATTR%elevation); attrStruct%elevation(iGRU) = parVector
      case(iLookATTR%tan_slope); attrStruct%tan_slope(iGRU) = parVector
      case(iLookATTR%contourLength); attrStruct%contourLength(iGRU) = parVector
      case(iLookATTR%HRUarea); attrStruct%HRUarea(iGRU) = parVector
      case(iLookATTR%mHeight); attrStruct%mHeight(iGRU) = parVector
      case(iLookATTR%aspect); attrStruct%aspect(iGRU) = parVector
    end select
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



    subroutine get_device_prog_stat(progStruct, iVar, nGRU,iFreq, parVector)
    implicit none
    type(prog_stat_device) :: progStruct
    integer(i4b) :: nGRU, iVar,iFreq
    real(rkind) :: parVector(:)

    parVector = -9999_rkind

    select case(iVar)
    case(iLookPROG%scalarCanopyIce); parVector = progStruct%stat_scalarCanopyIce(iFreq,:)
    case(iLookPROG%scalarCanopyLiq); parVector = progStruct%stat_scalarCanopyLiq(iFreq,:)
    case(iLookPROG%scalarCanopyWat); parVector = progStruct%stat_scalarCanopyWat(iFreq,:)
    case(iLookPROG%scalarCanairTemp); parVector = progStruct%stat_scalarCanairTemp(iFreq,:)
    case(iLookPROG%scalarCanopyTemp); parVector = progStruct%stat_scalarCanopyTemp(iFreq,:)
    case(iLookPROG%scalarSnowAlbedo); parVector = progStruct%stat_scalarSnowAlbedo(iFreq,:)
    case(iLookPROG%scalarSnowDepth); parVector = progStruct%stat_scalarSnowDepth(iFreq,:)
    case(iLookPROG%scalarSWE); parVector = progStruct%stat_scalarSWE(iFreq,:)
    case(iLookPROG%scalarSfcMeltPond); parVector = progStruct%stat_scalarSfcMeltPond(iFreq,:)
    case(iLookPROG%scalarCanairEnthalpy); parVector = progStruct%stat_scalarCanairEnthalpy(iFreq,:)
    case(iLookPROG%scalarCanopyEnthalpy); parVector = progStruct%stat_scalarCanopyEnthalpy(iFreq,:)
    case(iLookPROG%scalarAquiferStorage); parVector = progStruct%stat_scalarAquiferStorage(iFreq,:)
    case(iLookPROG%scalarSurfaceTemp); parVector = progStruct%stat_scalarSurfaceTemp(iFreq,:)
    end select
  end subroutine
subroutine get_device_forc_stat(forcStat,iVar,nGRU,iFreq,parVector)
  implicit none
  type(forc_stat_device) :: forcStat
  integer(i4b) :: nGRU, iVar, iFreq
  real(rkind) :: parVector(:)
  select case(iVar)
  case(iLookFORCE%time); parVector = forcStat%stat_time(iFreq,:)
  case(iLookFORCE%pptrate); parVector = forcStat%stat_pptrate(iFreq,:)
  case(iLookFORCE%SWRadAtm); parVector = forcStat%stat_SWRadAtm(iFreq,:)
  case(iLookFORCE%LWRadAtm); parVector = forcStat%stat_LWRadAtm(iFreq,:)
  case(iLookFORCE%airtemp); parVector = forcStat%stat_airtemp(iFreq,:)
  case(iLookFORCE%windspd); parVector = forcStat%stat_windspd(iFreq,:)
  case(iLookFORCE%airpres); parVector = forcStat%stat_airpres(iFreq,:)
  case(iLookFORCE%spechum); parVector = forcStat%stat_spechum(iFreq,:)
  case default; parVector = -9999_rkind
  end select
end subroutine

subroutine allocate_device_forc_stat(forcStat,nGRU)
  use globalData,only:statForc_meta
  use var_lookup,only:maxVarFreq
  implicit none
  type(forc_stat_device) :: forcStat
  integer(i4b) :: nGRU, iGRU, iVar, pVar

  do iVar=1,size(statForc_meta)
    pVar = statForc_meta(iVar)%ixParent
    select case(pVar)
    case(iLookFORCE%time); forcStat%statIndex_time = statForc_meta(iVar)%statIndex
      allocate(forcStat%stat_time(2*maxVarFreq,nGRU))
    case(iLookFORCE%pptrate); forcStat%statIndex_pptrate = statForc_meta(iVar)%statIndex
      allocate(forcStat%stat_pptrate(2*maxVarFreq,nGRU))
    case(iLookFORCE%SWRadAtm); forcStat%statIndex_SWRadAtm = statForc_meta(iVar)%statIndex
      allocate(forcStat%stat_SWRadAtm(2*maxVarFreq,nGRU))
    case(iLookFORCE%LWRadAtm); forcStat%statIndex_LWRadAtm = statForc_meta(iVar)%statIndex
      allocate(forcStat%stat_LWRadAtm(2*maxVarFreq,nGRU))
    case(iLookFORCE%airtemp); forcStat%statIndex_airtemp = statForc_meta(iVar)%statIndex
      allocate(forcStat%stat_airtemp(2*maxVarFreq,nGRU))
    case(iLookFORCE%windspd); forcStat%statIndex_windspd = statForc_meta(iVar)%statIndex
      allocate(forcStat%stat_windspd(2*maxVarFreq,nGRU))
    case(iLookFORCE%airpres); forcStat%statIndex_airpres = statForc_meta(iVar)%statIndex
      allocate(forcStat%stat_airpres(2*maxVarFreq,nGRU))
    case(iLookFORCE%spechum); forcStat%statIndex_spechum = statForc_meta(iVar)%statIndex
      allocate(forcStat%stat_spechum(2*maxVarFreq,nGRU))
    end select
  end do
  if (.not.allocated(forcStat%statIndex_time)) allocate(forcStat%statIndex_time(1))
  if (.not.allocated(forcStat%stat_time)) allocate(forcStat%stat_time(1,1))
if (.not.allocated(forcStat%statIndex_pptrate)) allocate(forcStat%statIndex_pptrate(1))
if (.not.allocated(forcStat%stat_pptrate)) allocate(forcStat%stat_pptrate(1,1))
if (.not.allocated(forcStat%statIndex_SWRadAtm)) allocate(forcStat%statIndex_SWRadAtm(1))
if (.not.allocated(forcStat%stat_SWRadAtm)) allocate(forcStat%stat_SWRadAtm(1,1))
if (.not.allocated(forcStat%statIndex_LWRadAtm)) allocate(forcStat%statIndex_LWRadAtm(1))
if (.not.allocated(forcStat%stat_LWRadAtm)) allocate(forcStat%stat_LWRadAtm(1,1))
if (.not.allocated(forcStat%statIndex_airtemp)) allocate(forcStat%statIndex_airtemp(1))
if (.not.allocated(forcStat%stat_airtemp)) allocate(forcStat%stat_airtemp(1,1))
if (.not.allocated(forcStat%statIndex_windspd)) allocate(forcStat%statIndex_windspd(1))
if (.not.allocated(forcStat%stat_windspd)) allocate(forcStat%stat_windspd(1,1))
if (.not.allocated(forcStat%statIndex_airpres)) allocate(forcStat%statIndex_airpres(1))
if (.not.allocated(forcStat%stat_airpres)) allocate(forcStat%stat_airpres(1,1))
if (.not.allocated(forcStat%statIndex_spechum)) allocate(forcStat%statIndex_spechum(1))
if (.not.allocated(forcStat%stat_spechum)) allocate(forcStat%stat_spechum(1,1))
end subroutine



  subroutine allocate_device_prog_stat(progStat_d,nGRU)
    USE globalData,only:statProg_meta             ! child metadata for stats
    USE var_lookup,only:maxvarFreq             ! allocation dimension (output frequency)
    implicit none
    type(prog_stat_device) :: progStat_d
    ! type(gru_doubleVec) :: bvarStat

    integer(i4b) :: nGRU, iGRU, iVar, pVar
    integer :: err

    err = cudaDeviceSynchronize()
 do iVar=1,size(statProg_meta)
  pVar = statProg_meta(iVar)%ixParent
  select case(pVar)
  case(iLookPROG%scalarCanopyIce); progStat_d%statIndex_scalarCanopyIce = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarCanopyIce(maxVarFreq*2,nGRU))
  case(iLookPROG%scalarCanopyLiq); progStat_d%statIndex_scalarCanopyLiq = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarCanopyLiq(maxVarFreq*2,nGRU))
  case(iLookPROG%scalarCanopyWat); progStat_d%statIndex_scalarCanopyWat = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarCanopyWat(maxVarFreq*2,nGRU))
  case(iLookPROG%scalarCanairTemp); progStat_d%statIndex_scalarCanairTemp = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarCanairTemp(maxVarFreq*2,nGRU))
  case(iLookPROG%scalarCanopyTemp); progStat_d%statIndex_scalarCanopyTemp = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarCanopyTemp(maxVarFreq*2,nGRU))
  case(iLookPROG%scalarSnowAlbedo); progStat_d%statIndex_scalarSnowAlbedo = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarSnowAlbedo(maxVarFreq*2,nGRU))
  case(iLookPROG%scalarSnowDepth); progStat_d%statIndex_scalarSnowDepth = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarSnowDepth(maxVarFreq*2,nGRU))
  case(iLookPROG%scalarSWE); progStat_d%statIndex_scalarSWE = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarSWE(maxVarFreq*2,nGRU))
  case(iLookPROG%scalarSfcMeltPond); progStat_d%statIndex_scalarSfcMeltPond = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarSfcMeltPond(maxVarFreq*2,nGRU))
  case(iLookPROG%scalarCanairEnthalpy); progStat_d%statIndex_scalarCanairEnthalpy = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarCanairEnthalpy(maxVarFreq*2,nGRU))
  case(iLookPROG%scalarCanopyEnthalpy); progStat_d%statIndex_scalarCanopyEnthalpy = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarCanopyEnthalpy(maxVarFreq*2,nGRU))
  case(iLookPROG%scalarAquiferStorage); progStat_d%statIndex_scalarAquiferStorage = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarAquiferStorage(maxVarFreq*2,nGRU))
  case(iLookPROG%scalarSurfaceTemp); progStat_d%statIndex_scalarSurfaceTemp = statProg_meta(iVar)%statIndex
    allocate(progStat_d%stat_scalarSurfaceTemp(maxVarFreq*2,nGRU))
    end select
  end do
  if (.not. allocated(progStat_d%statIndex_scalarCanopyIce)) allocate(progStat_d%statIndex_scalarCanopyIce(1))
  if (.not.allocated(progStat_d%stat_scalarCanopyIce)) allocate(progStat_d%stat_scalarCanopyIce(1,1))
if (.not. allocated(progStat_d%statIndex_scalarCanopyLiq)) allocate(progStat_d%statIndex_scalarCanopyLiq(1))
if (.not.allocated(progStat_d%stat_scalarCanopyLiq)) allocate(progStat_d%stat_scalarCanopyLiq(1,1))
if (.not. allocated(progStat_d%statIndex_scalarCanopyWat)) allocate(progStat_d%statIndex_scalarCanopyWat(1))
if (.not.allocated(progStat_d%stat_scalarCanopyWat)) allocate(progStat_d%stat_scalarCanopyWat(1,1))
if (.not. allocated(progStat_d%statIndex_scalarCanairTemp)) allocate(progStat_d%statIndex_scalarCanairTemp(1))
if (.not.allocated(progStat_d%stat_scalarCanairTemp)) allocate(progStat_d%stat_scalarCanairTemp(1,1))
if (.not. allocated(progStat_d%statIndex_scalarCanopyTemp)) allocate(progStat_d%statIndex_scalarCanopyTemp(1))
if (.not.allocated(progStat_d%stat_scalarCanopyTemp)) allocate(progStat_d%stat_scalarCanopyTemp(1,1))
if (.not. allocated(progStat_d%statIndex_scalarSnowAlbedo)) allocate(progStat_d%statIndex_scalarSnowAlbedo(1))
if (.not.allocated(progStat_d%stat_scalarSnowAlbedo)) allocate(progStat_d%stat_scalarSnowAlbedo(1,1))
if (.not. allocated(progStat_d%statIndex_scalarSnowDepth)) allocate(progStat_d%statIndex_scalarSnowDepth(1))
if (.not.allocated(progStat_d%stat_scalarSnowDepth)) allocate(progStat_d%stat_scalarSnowDepth(1,1))
if (.not. allocated(progStat_d%statIndex_scalarSWE)) allocate(progStat_d%statIndex_scalarSWE(1))
if (.not.allocated(progStat_d%stat_scalarSWE)) allocate(progStat_d%stat_scalarSWE(1,1))
if (.not. allocated(progStat_d%statIndex_scalarSfcMeltPond)) allocate(progStat_d%statIndex_scalarSfcMeltPond(1))
if (.not.allocated(progStat_d%stat_scalarSfcMeltPond)) allocate(progStat_d%stat_scalarSfcMeltPond(1,1))
if (.not. allocated(progStat_d%statIndex_scalarCanairEnthalpy)) allocate(progStat_d%statIndex_scalarCanairEnthalpy(1))
if (.not.allocated(progStat_d%stat_scalarCanairEnthalpy)) allocate(progStat_d%stat_scalarCanairEnthalpy(1,1))
if (.not. allocated(progStat_d%statIndex_scalarCanopyEnthalpy)) allocate(progStat_d%statIndex_scalarCanopyEnthalpy(1))
if (.not.allocated(progStat_d%stat_scalarCanopyEnthalpy)) allocate(progStat_d%stat_scalarCanopyEnthalpy(1,1))
if (.not. allocated(progStat_d%statIndex_scalarAquiferStorage)) allocate(progStat_d%statIndex_scalarAquiferStorage(1))
if (.not.allocated(progStat_d%stat_scalarAquiferStorage)) allocate(progStat_d%stat_scalarAquiferStorage(1,1))
if (.not. allocated(progStat_d%statIndex_scalarSurfaceTemp)) allocate(progStat_d%statIndex_scalarSurfaceTemp(1))
if (.not.allocated(progStat_d%stat_scalarSurfaceTemp)) allocate(progStat_d%stat_scalarSurfaceTemp(1,1))
!   if (.not. allocated(bvarStat_d%statIndex_basin__SurfaceRunoff)) allocate(bvarStat_d%statIndex_basin__SurfaceRunoff(1))
! if (.not. allocated(bvarStat_d%statIndex_basin__ColumnOutflow)) allocate(bvarStat_d%statIndex_basin__ColumnOutflow(1))
! if (.not. allocated(bvarStat_d%statIndex_basin__AquiferStorage)) allocate(bvarStat_d%statIndex_basin__AquiferStorage(1))
! if (.not. allocated(bvarStat_d%statIndex_basin__AquiferRecharge)) allocate(bvarStat_d%statIndex_basin__AquiferRecharge(1))
! if (.not. allocated(bvarStat_d%statIndex_basin__AquiferBaseflow)) allocate(bvarStat_d%statIndex_basin__AquiferBaseflow(1))
! if (.not. allocated(bvarStat_d%statIndex_basin__AquiferTranspire)) allocate(bvarStat_d%statIndex_basin__AquiferTranspire(1))
! if (.not. allocated(bvarStat_d%statIndex_basin__TotalRunoff)) allocate(bvarStat_d%statIndex_basin__TotalRunoff(1))
! if (.not. allocated(bvarStat_d%statIndex_basin__SoilDrainage)) allocate(bvarStat_d%statIndex_basin__SoilDrainage(1))
! if (.not. allocated(bvarStat_d%statIndex_basin__totalArea)) allocate(bvarStat_d%statIndex_basin__totalArea(1))
! if (.not. allocated(bvarStat_d%statIndex_averageInstantRunoff)) allocate(bvarStat_d%statIndex_averageInstantRunoff(1))
! if (.not. allocated(bvarStat_d%statIndex_averageRoutedRunoff)) allocate(bvarStat_d%statIndex_averageRoutedRunoff(1))

  end subroutine


        subroutine get_device_indx_stat(indxStruct, iVar, nGRU,iFreq, parVector)
    implicit none
    type(indx_stat_device) :: indxStruct
    integer(i4b) :: nGRU, iVar,iFreq
    real(rkind) :: parVector(:)

    parVector = -9999_rkind
    select case(iVar)
    case(iLookINDEX%nSnow); parVector = indxStruct%stat_nSnow(iFreq,:)
    case(iLookINDEX%nSoil); parVector = indxStruct%stat_nSoil(iFreq,:)
    case(iLookINDEX%nLayers); parVector = indxStruct%stat_nLayers(iFreq,:)
    case(iLookINDEX%nCasNrg); parVector = indxStruct%stat_nCasNrg(iFreq,:)
    case(iLookINDEX%nVegNrg); parVector = indxStruct%stat_nVegNrg(iFreq,:)
    case(iLookINDEX%nVegMass); parVector = indxStruct%stat_nVegMass(iFreq,:)
    case(iLookINDEX%nVegState); parVector = indxStruct%stat_nVegState(iFreq,:)
    case(iLookINDEX%nNrgState); parVector = indxStruct%stat_nNrgState(iFreq,:)
    case(iLookINDEX%nWatState); parVector = indxStruct%stat_nWatState(iFreq,:)
    case(iLookINDEX%nMatState); parVector = indxStruct%stat_nMatState(iFreq,:)
    case(iLookINDEX%nMassState); parVector = indxStruct%stat_nMassState(iFreq,:)
    case(iLookINDEX%nState); parVector = indxStruct%stat_nState(iFreq,:)
    ! case(iLookINDEX%nSnowSoilNrg); allocate(parVector(nGRU)); parVector = indxStruct%stat_nSnowSoilNrg(iFreq,:)
    ! case(iLookINDEX%nSnowOnlyNrg); allocate(parVector(nGRU)); parVector = indxStruct%stat_nSnowOnlyNrg(iFreq,:)
    ! case(iLookINDEX%nSoilOnlyNrg); allocate(parVector(nGRU)); parVector = indxStruct%stat_nSoilOnlyNrg(iFreq,:)
    ! case(iLookINDEX%nSnowSoilHyd); allocate(parVector(nGRU)); parVector = indxStruct%stat_nSnowSoilHyd(iFreq,:)
    ! case(iLookINDEX%nSnowOnlyHyd); allocate(parVector(nGRU)); parVector = indxStruct%stat_nSnowOnlyHyd(iFreq,:)
    ! case(iLookINDEX%nSoilOnlyHyd); allocate(parVector(nGRU)); parVector = indxStruct%stat_nSoilOnlyHyd(iFreq,:)
    case(iLookINDEX%ixCasNrg); parVector = indxStruct%stat_ixCasNrg(iFreq,:)
    case(iLookINDEX%ixVegNrg); parVector = indxStruct%stat_ixVegNrg(iFreq,:)
    case(iLookINDEX%ixVegHyd); parVector = indxStruct%stat_ixVegHyd(iFreq,:)
    case(iLookINDEX%ixTopNrg); parVector = indxStruct%stat_ixTopNrg(iFreq,:)
    case(iLookINDEX%ixTopHyd); parVector = indxStruct%stat_ixTopHyd(iFreq,:)
    case(iLookINDEX%ixAqWat); parVector = indxStruct%stat_ixAqWat(iFreq,:)
    case(iLookINDEX%ixNrgCanair); parVector = indxStruct%stat_ixNrgCanair(iFreq,:)
    case(iLookINDEX%ixNrgCanopy); parVector = indxStruct%stat_ixNrgCanopy(iFreq,:)
    case(iLookINDEX%ixHydCanopy); parVector = indxStruct%stat_ixHydCanopy(iFreq,:)
    case(iLookINDEX%ixWatAquifer); parVector = indxStruct%stat_ixWatAquifer(iFreq,:)
    case(iLookINDEX%numberFluxCalc); parVector = indxStruct%stat_numberFluxCalc(iFreq,:)
    case(iLookINDEX%numberStateSplit); parVector = indxStruct%stat_numberStateSplit(iFreq,:)
    case(iLookINDEX%numberDomainSplitNrg); parVector = indxStruct%stat_numberDomainSplitNrg(iFreq,:)
    case(iLookINDEX%numberDomainSplitMass); parVector = indxStruct%stat_numberDomainSplitMass(iFreq,:)
    case(iLookINDEX%numberScalarSolutions); parVector = indxStruct%stat_numberScalarSolutions(iFreq,:)
    end select
  end subroutine
  subroutine allocate_device_indx_stat(indxStat_d,nGRU)
        USE globalData,only:statIndx_meta             ! child metadata for stats
    USE var_lookup,only:maxvarFreq             ! allocation dimension (output frequency)
    implicit none
    type(indx_stat_device) :: indxStat_d
    ! type(gru_doubleVec) :: bvarStat

    integer(i4b) :: nGRU, iGRU, iVar, pVar
    integer :: err

 do iVar=1,size(statIndx_meta)
  pVar = statIndx_meta(iVar)%ixParent
  select case(pVar)
  case(iLookINDEX%nSnow); indxStat_d%statIndex_nSnow = statIndx_meta(iVar)%statIndex
        allocate(indxStat_d%stat_nSnow(maxVarFreq*2,nGRU))
  case(iLookINDEX%nSoil); indxStat_d%statIndex_nSoil = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_nSoil(maxVarFreq*2,nGRU))
  case(iLookINDEX%nLayers); indxStat_d%statIndex_nLayers = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_nLayers(maxVarFreq*2,nGRU))
  case(iLookINDEX%nCasNrg); indxStat_d%statIndex_nCasNrg = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_nCasNrg(maxVarFreq*2,nGRU))
  case(iLookINDEX%nVegNrg); indxStat_d%statIndex_nVegNrg = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_nVegNrg(maxVarFreq*2,nGRU))
  case(iLookINDEX%nVegMass); indxStat_d%statIndex_nVegMass = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_nVegMass(maxVarFreq*2,nGRU))
  case(iLookINDEX%nVegState); indxStat_d%statIndex_nVegState = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_nVegState(maxVarFreq*2,nGRU))
  case(iLookINDEX%nNrgState); indxStat_d%statIndex_nNrgState = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_nNrgState(maxVarFreq*2,nGRU))
  case(iLookINDEX%nWatState); indxStat_d%statIndex_nWatState = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_nWatState(maxVarFreq*2,nGRU))
  case(iLookINDEX%nMatState); indxStat_d%statIndex_nMatState = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_nMatState(maxVarFreq*2,nGRU))
  case(iLookINDEX%nMassState); indxStat_d%statIndex_nMassState = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_nMassState(maxVarFreq*2,nGRU))
  case(iLookINDEX%nState); indxStat_d%statIndex_nState = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_nState(maxVarFreq*2,nGRU))
  case(iLookINDEX%ixCasNrg); indxStat_d%statIndex_ixCasNrg = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_ixCasNrg(maxVarFreq*2,nGRU))
  case(iLookINDEX%ixVegNrg); indxStat_d%statIndex_ixVegNrg = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_ixVegNrg(maxVarFreq*2,nGRU))
  case(iLookINDEX%ixVegHyd); indxStat_d%statIndex_ixVegHyd = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_ixVegHyd(maxVarFreq*2,nGRU))
  case(iLookINDEX%ixTopNrg); indxStat_d%statIndex_ixTopNrg = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_ixTopNrg(maxVarFreq*2,nGRU))
  case(iLookINDEX%ixTopHyd); indxStat_d%statIndex_ixTopHyd = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_ixTopHyd(maxVarFreq*2,nGRU))
  case(iLookINDEX%ixAqWat); indxStat_d%statIndex_ixAqWat = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_ixAqWat(maxVarFreq*2,nGRU))
  case(iLookINDEX%ixNrgCanair); indxStat_d%statIndex_ixNrgCanair = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_ixNrgCanair(maxVarFreq*2,nGRU))
  case(iLookINDEX%ixNrgCanopy); indxStat_d%statIndex_ixNrgCanopy = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_ixNrgCanopy(maxVarFreq*2,nGRU))
  case(iLookINDEX%ixHydCanopy); indxStat_d%statIndex_ixHydCanopy = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_ixHydCanopy(maxVarFreq*2,nGRU))
  case(iLookINDEX%ixWatAquifer); indxStat_d%statIndex_ixWatAquifer = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_ixWatAquifer(maxVarFreq*2,nGRU))
  case(iLookINDEX%numberFluxCalc); indxStat_d%statIndex_numberFluxCalc = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_numberFluxCalc(maxVarFreq*2,nGRU))
  case(iLookINDEX%numberStateSplit); indxStat_d%statIndex_numberStateSplit = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_numberStateSplit(maxVarFreq*2,nGRU))
  case(iLookINDEX%numberDomainSplitNrg); indxStat_d%statIndex_numberDomainSplitNrg = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_numberDomainSplitNrg(maxVarFreq*2,nGRU))
  case(iLookINDEX%numberDomainSplitMass); indxStat_d%statIndex_numberDomainSplitMass = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_numberDomainSplitMass(maxVarFreq*2,nGRU))
  case(iLookINDEX%numberScalarSolutions); indxStat_d%statIndex_numberScalarSolutions = statIndx_meta(iVar)%statIndex
    allocate(indxStat_d%stat_numberScalarSolutions(maxVarFreq*2,nGRU))

  end select
end do

if (.not.allocated(indxStat_d%statIndex_nSnow)) allocate(indxStat_d%statIndex_nSnow(1))
if (.not.allocated(indxStat_d%stat_nSnow)) allocate(indxStat_d%stat_nSnow(1,1))
if (.not.allocated(indxStat_d%statIndex_nSoil)) allocate(indxStat_d%statIndex_nSoil(1))
if (.not.allocated(indxStat_d%stat_nSoil)) allocate(indxStat_d%stat_nSoil(1,1))
if (.not.allocated(indxStat_d%statIndex_nLayers)) allocate(indxStat_d%statIndex_nLayers(1))
if (.not.allocated(indxStat_d%stat_nLayers)) allocate(indxStat_d%stat_nLayers(1,1))
if (.not.allocated(indxStat_d%statIndex_nCasNrg)) allocate(indxStat_d%statIndex_nCasNrg(1))
if (.not.allocated(indxStat_d%stat_nCasNrg)) allocate(indxStat_d%stat_nCasNrg(1,1))
if (.not.allocated(indxStat_d%statIndex_nVegNrg)) allocate(indxStat_d%statIndex_nVegNrg(1))
if (.not.allocated(indxStat_d%stat_nVegNrg)) allocate(indxStat_d%stat_nVegNrg(1,1))
if (.not.allocated(indxStat_d%statIndex_nVegMass)) allocate(indxStat_d%statIndex_nVegMass(1))
if (.not.allocated(indxStat_d%stat_nVegMass)) allocate(indxStat_d%stat_nVegMass(1,1))
if (.not.allocated(indxStat_d%statIndex_nVegState)) allocate(indxStat_d%statIndex_nVegState(1))
if (.not.allocated(indxStat_d%stat_nVegState)) allocate(indxStat_d%stat_nVegState(1,1))
if (.not.allocated(indxStat_d%statIndex_nNrgState)) allocate(indxStat_d%statIndex_nNrgState(1))
if (.not.allocated(indxStat_d%stat_nNrgState)) allocate(indxStat_d%stat_nNrgState(1,1))
if (.not.allocated(indxStat_d%statIndex_nWatState)) allocate(indxStat_d%statIndex_nWatState(1))
if (.not.allocated(indxStat_d%stat_nWatState)) allocate(indxStat_d%stat_nWatState(1,1))
if (.not.allocated(indxStat_d%statIndex_nMatState)) allocate(indxStat_d%statIndex_nMatState(1))
if (.not.allocated(indxStat_d%stat_nMatState)) allocate(indxStat_d%stat_nMatState(1,1))
if (.not.allocated(indxStat_d%statIndex_nMassState)) allocate(indxStat_d%statIndex_nMassState(1))
if (.not.allocated(indxStat_d%stat_nMassState)) allocate(indxStat_d%stat_nMassState(1,1))
if (.not.allocated(indxStat_d%statIndex_nState)) allocate(indxStat_d%statIndex_nState(1))
if (.not.allocated(indxStat_d%stat_nState)) allocate(indxStat_d%stat_nState(1,1))
if (.not.allocated(indxStat_d%statIndex_ixCasNrg)) allocate(indxStat_d%statIndex_ixCasNrg(1))
if (.not.allocated(indxStat_d%stat_ixCasNrg)) allocate(indxStat_d%stat_ixCasNrg(1,1))
if (.not.allocated(indxStat_d%statIndex_ixVegNrg)) allocate(indxStat_d%statIndex_ixVegNrg(1))
if (.not.allocated(indxStat_d%stat_ixVegNrg)) allocate(indxStat_d%stat_ixVegNrg(1,1))
if (.not.allocated(indxStat_d%statIndex_ixVegHyd)) allocate(indxStat_d%statIndex_ixVegHyd(1))
if (.not.allocated(indxStat_d%stat_ixVegHyd)) allocate(indxStat_d%stat_ixVegHyd(1,1))
if (.not.allocated(indxStat_d%statIndex_ixTopNrg)) allocate(indxStat_d%statIndex_ixTopNrg(1))
if (.not.allocated(indxStat_d%stat_ixTopNrg)) allocate(indxStat_d%stat_ixTopNrg(1,1))
if (.not.allocated(indxStat_d%statIndex_ixTopHyd)) allocate(indxStat_d%statIndex_ixTopHyd(1))
if (.not.allocated(indxStat_d%stat_ixTopHyd)) allocate(indxStat_d%stat_ixTopHyd(1,1))
if (.not.allocated(indxStat_d%statIndex_ixAqWat)) allocate(indxStat_d%statIndex_ixAqWat(1))
if (.not.allocated(indxStat_d%stat_ixAqWat)) allocate(indxStat_d%stat_ixAqWat(1,1))
if (.not.allocated(indxStat_d%statIndex_ixNrgCanair)) allocate(indxStat_d%statIndex_ixNrgCanair(1))
if (.not.allocated(indxStat_d%stat_ixNrgCanair)) allocate(indxStat_d%stat_ixNrgCanair(1,1))
if (.not.allocated(indxStat_d%statIndex_ixNrgCanopy)) allocate(indxStat_d%statIndex_ixNrgCanopy(1))
if (.not.allocated(indxStat_d%stat_ixNrgCanopy)) allocate(indxStat_d%stat_ixNrgCanopy(1,1))
if (.not.allocated(indxStat_d%statIndex_ixHydCanopy)) allocate(indxStat_d%statIndex_ixHydCanopy(1))
if (.not.allocated(indxStat_d%stat_ixHydCanopy)) allocate(indxStat_d%stat_ixHydCanopy(1,1))
if (.not.allocated(indxStat_d%statIndex_ixWatAquifer)) allocate(indxStat_d%statIndex_ixWatAquifer(1))
if (.not.allocated(indxStat_d%stat_ixWatAquifer)) allocate(indxStat_d%stat_ixWatAquifer(1,1))
if (.not.allocated(indxStat_d%statIndex_numberFluxCalc)) allocate(indxStat_d%statIndex_numberFluxCalc(1))
if (.not.allocated(indxStat_d%stat_numberFluxCalc)) allocate(indxStat_d%stat_numberFluxCalc(1,1))
if (.not.allocated(indxStat_d%statIndex_numberStateSplit)) allocate(indxStat_d%statIndex_numberStateSplit(1))
if (.not.allocated(indxStat_d%stat_numberStateSplit)) allocate(indxStat_d%stat_numberStateSplit(1,1))
if (.not.allocated(indxStat_d%statIndex_numberDomainSplitNrg)) allocate(indxStat_d%statIndex_numberDomainSplitNrg(1))
if (.not.allocated(indxStat_d%stat_numberDomainSplitNrg)) allocate(indxStat_d%stat_numberDomainSplitNrg(1,1))
if (.not.allocated(indxStat_d%statIndex_numberDomainSplitMass)) allocate(indxStat_d%statIndex_numberDomainSplitMass(1))
if (.not.allocated(indxStat_d%stat_numberDomainSplitMass)) allocate(indxStat_d%stat_numberDomainSplitMass(1,1))
if (.not.allocated(indxStat_d%statIndex_numberScalarSolutions)) allocate(indxStat_d%statIndex_numberScalarSolutions(1))
if (.not.allocated(indxStat_d%stat_numberScalarSolutions)) allocate(indxStat_d%stat_numberScalarSolutions(1,1))


  end subroutine

  subroutine get_device_flux_stat(fluxStruct, iVar, nGRU, iFreq, parVector)
    implicit none
    type(flux_stat_device) :: fluxStruct
    integer(i4b) :: nGRU, iVar, iFreq
    real(rkind) :: parVector(:)

select case(iVar)
  case(iLookFLUX%scalarCanairNetNrgFlux); parVector = fluxStruct%stat_scalarCanairNetNrgFlux(iFreq,:)
case(iLookFLUX%scalarCanopyNetNrgFlux); parVector = fluxStruct%stat_scalarCanopyNetNrgFlux(iFreq,:)
case(iLookFLUX%scalarGroundNetNrgFlux); parVector = fluxStruct%stat_scalarGroundNetNrgFlux(iFreq,:)
case(iLookFLUX%scalarCanopyNetLiqFlux); parVector = fluxStruct%stat_scalarCanopyNetLiqFlux(iFreq,:)
case(iLookFLUX%scalarRainfall); parVector = fluxStruct%stat_scalarRainfall(iFreq,:)
case(iLookFLUX%scalarSnowfall); parVector = fluxStruct%stat_scalarSnowfall(iFreq,:)
case(iLookFLUX%scalarCanopySunlitPAR); parVector = fluxStruct%stat_scalarCanopySunlitPAR(iFreq,:)
case(iLookFLUX%scalarCanopyShadedPAR); parVector = fluxStruct%stat_scalarCanopyShadedPAR(iFreq,:)
case(iLookFLUX%scalarBelowCanopySolar); parVector = fluxStruct%stat_scalarBelowCanopySolar(iFreq,:)
case(iLookFLUX%scalarCanopyAbsorbedSolar); parVector = fluxStruct%stat_scalarCanopyAbsorbedSolar(iFreq,:)
case(iLookFLUX%scalarGroundAbsorbedSolar); parVector = fluxStruct%stat_scalarGroundAbsorbedSolar(iFreq,:)
case(iLookFLUX%scalarLWRadCanopy); parVector = fluxStruct%stat_scalarLWRadCanopy(iFreq,:)
case(iLookFLUX%scalarLWRadGround); parVector = fluxStruct%stat_scalarLWRadGround(iFreq,:)
case(iLookFLUX%scalarLWRadUbound2Canopy); parVector = fluxStruct%stat_scalarLWRadUbound2Canopy(iFreq,:)
case(iLookFLUX%scalarLWRadUbound2Ground); parVector = fluxStruct%stat_scalarLWRadUbound2Ground(iFreq,:)
case(iLookFLUX%scalarLWRadUbound2Ubound); parVector = fluxStruct%stat_scalarLWRadUbound2Ubound(iFreq,:)
case(iLookFLUX%scalarLWRadCanopy2Ubound); parVector = fluxStruct%stat_scalarLWRadCanopy2Ubound(iFreq,:)
case(iLookFLUX%scalarLWRadCanopy2Ground); parVector = fluxStruct%stat_scalarLWRadCanopy2Ground(iFreq,:)
case(iLookFLUX%scalarLWRadCanopy2Canopy); parVector = fluxStruct%stat_scalarLWRadCanopy2Canopy(iFreq,:)
case(iLookFLUX%scalarLWRadGround2Ubound); parVector = fluxStruct%stat_scalarLWRadGround2Ubound(iFreq,:)
case(iLookFLUX%scalarLWRadGround2Canopy); parVector = fluxStruct%stat_scalarLWRadGround2Canopy(iFreq,:)
case(iLookFLUX%scalarLWNetCanopy); parVector = fluxStruct%stat_scalarLWNetCanopy(iFreq,:)
case(iLookFLUX%scalarLWNetGround); parVector = fluxStruct%stat_scalarLWNetGround(iFreq,:)
case(iLookFLUX%scalarLWNetUbound); parVector = fluxStruct%stat_scalarLWNetUbound(iFreq,:)
case(iLookFLUX%scalarEddyDiffusCanopyTop); parVector = fluxStruct%stat_scalarEddyDiffusCanopyTop(iFreq,:)
case(iLookFLUX%scalarFrictionVelocity); parVector = fluxStruct%stat_scalarFrictionVelocity(iFreq,:)
case(iLookFLUX%scalarWindspdCanopyTop); parVector = fluxStruct%stat_scalarWindspdCanopyTop(iFreq,:)
case(iLookFLUX%scalarWindspdCanopyBottom); parVector = fluxStruct%stat_scalarWindspdCanopyBottom(iFreq,:)
case(iLookFLUX%scalarGroundResistance); parVector = fluxStruct%stat_scalarGroundResistance(iFreq,:)
case(iLookFLUX%scalarCanopyResistance); parVector = fluxStruct%stat_scalarCanopyResistance(iFreq,:)
case(iLookFLUX%scalarLeafResistance); parVector = fluxStruct%stat_scalarLeafResistance(iFreq,:)
case(iLookFLUX%scalarSoilResistance); parVector = fluxStruct%stat_scalarSoilResistance(iFreq,:)
case(iLookFLUX%scalarSenHeatTotal); parVector = fluxStruct%stat_scalarSenHeatTotal(iFreq,:)
case(iLookFLUX%scalarSenHeatCanopy); parVector = fluxStruct%stat_scalarSenHeatCanopy(iFreq,:)
case(iLookFLUX%scalarSenHeatGround); parVector = fluxStruct%stat_scalarSenHeatGround(iFreq,:)
case(iLookFLUX%scalarLatHeatTotal); parVector = fluxStruct%stat_scalarLatHeatTotal(iFreq,:)
case(iLookFLUX%scalarLatHeatCanopyEvap); parVector = fluxStruct%stat_scalarLatHeatCanopyEvap(iFreq,:)
case(iLookFLUX%scalarLatHeatCanopyTrans); parVector = fluxStruct%stat_scalarLatHeatCanopyTrans(iFreq,:)
case(iLookFLUX%scalarLatHeatGround); parVector = fluxStruct%stat_scalarLatHeatGround(iFreq,:)
case(iLookFLUX%scalarCanopyAdvectiveHeatFlux); parVector = fluxStruct%stat_scalarCanopyAdvectiveHeatFlux(iFreq,:)
case(iLookFLUX%scalarGroundAdvectiveHeatFlux); parVector = fluxStruct%stat_scalarGroundAdvectiveHeatFlux(iFreq,:)
case(iLookFLUX%scalarCanopySublimation); parVector = fluxStruct%stat_scalarCanopySublimation(iFreq,:)
case(iLookFLUX%scalarSnowSublimation); parVector = fluxStruct%stat_scalarSnowSublimation(iFreq,:)
case(iLookFLUX%scalarStomResistSunlit); parVector = fluxStruct%stat_scalarStomResistSunlit(iFreq,:)
case(iLookFLUX%scalarStomResistShaded); parVector = fluxStruct%stat_scalarStomResistShaded(iFreq,:)
case(iLookFLUX%scalarPhotosynthesisSunlit); parVector = fluxStruct%stat_scalarPhotosynthesisSunlit(iFreq,:)
case(iLookFLUX%scalarPhotosynthesisShaded); parVector = fluxStruct%stat_scalarPhotosynthesisShaded(iFreq,:)
case(iLookFLUX%scalarCanopyTranspiration); parVector = fluxStruct%stat_scalarCanopyTranspiration(iFreq,:)
case(iLookFLUX%scalarCanopyEvaporation); parVector = fluxStruct%stat_scalarCanopyEvaporation(iFreq,:)
case(iLookFLUX%scalarGroundEvaporation); parVector = fluxStruct%stat_scalarGroundEvaporation(iFreq,:)
case(iLookFLUX%scalarThroughfallSnow); parVector = fluxStruct%stat_scalarThroughfallSnow(iFreq,:)
case(iLookFLUX%scalarThroughfallRain); parVector = fluxStruct%stat_scalarThroughfallRain(iFreq,:)
case(iLookFLUX%scalarCanopySnowUnloading); parVector = fluxStruct%stat_scalarCanopySnowUnloading(iFreq,:)
case(iLookFLUX%scalarCanopyLiqDrainage);  parVector = fluxStruct%stat_scalarCanopyLiqDrainage(iFreq,:)
case(iLookFLUX%scalarCanopyMeltFreeze);  parVector = fluxStruct%stat_scalarCanopyMeltFreeze(iFreq,:)
case(iLookFLUX%scalarSnowDrainage);  parVector = fluxStruct%stat_scalarSnowDrainage(iFreq,:)
case(iLookFLUX%scalarRainPlusMelt);  parVector = fluxStruct%stat_scalarRainPlusMelt(iFreq,:)
case(iLookFLUX%scalarMaxInfilRate);  parVector = fluxStruct%stat_scalarMaxInfilRate(iFreq,:)
case(iLookFLUX%scalarInfiltration);  parVector = fluxStruct%stat_scalarInfiltration(iFreq,:)
case(iLookFLUX%scalarExfiltration);  parVector = fluxStruct%stat_scalarExfiltration(iFreq,:)
case(iLookFLUX%scalarSurfaceRunoff);  parVector = fluxStruct%stat_scalarSurfaceRunoff(iFreq,:)
case(iLookFLUX%scalarSurfaceRunoff_IE);  parVector = fluxStruct%stat_scalarSurfaceRunoff_IE(iFreq,:)
case(iLookFLUX%scalarSurfaceRunoff_SE);  parVector = fluxStruct%stat_scalarSurfaceRunoff_SE(iFreq,:)
case(iLookFLUX%scalarSoilBaseflow);  parVector = fluxStruct%stat_scalarSoilBaseflow(iFreq,:)
case(iLookFLUX%scalarSoilDrainage);  parVector = fluxStruct%stat_scalarSoilDrainage(iFreq,:)
case(iLookFLUX%scalarAquiferRecharge); parVector = fluxStruct%stat_scalarAquiferRecharge(iFreq,:)
case(iLookFLUX%scalarAquiferTranspire); parVector = fluxStruct%stat_scalarAquiferTranspire(iFreq,:)
case(iLookFLUX%scalarAquiferBaseflow); parVector = fluxStruct%stat_scalarAquiferBaseflow(iFreq,:)
case(iLookFLUX%scalarTotalET); parVector = fluxStruct%stat_scalarTotalET(iFreq,:)
case(iLookFLUX%scalarTotalRunoff); parVector = fluxStruct%stat_scalarTotalRunoff(iFreq,:)
case(iLookFLUX%scalarNetRadiation); parVector = fluxStruct%stat_scalarNetRadiation(iFreq,:)
case default; parVector = -9999_rkind
end select
end subroutine
subroutine allocate_device_flux_stat(fluxStat,nGRU)
          USE globalData,only:statFlux_meta             ! child metadata for stats
    USE var_lookup,only:maxvarFreq             ! allocation dimension (output frequency)
    implicit none
    type(flux_stat_device) :: fluxStat
    ! type(gru_doubleVec) :: bvarStat

    integer(i4b) :: nGRU, iGRU, iVar, pVar

do iVar = 1,size(statFlux_meta)
  pVar = statFlux_meta(iVar)%ixParent
  select case (pVar)
  case(iLookFLUX%scalarCanairNetNrgFlux); fluxStat%statIndex_scalarCanairNetNrgFlux = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanairNetNrgFlux(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopyNetNrgFlux); fluxStat%statIndex_scalarCanopyNetNrgFlux = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopyNetNrgFlux(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarGroundNetNrgFlux); fluxStat%statIndex_scalarGroundNetNrgFlux = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarGroundNetNrgFlux(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopyNetLiqFlux); fluxStat%statIndex_scalarCanopyNetLiqFlux = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopyNetLiqFlux(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarRainfall); fluxStat%statIndex_scalarRainfall = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarRainfall(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarSnowfall); fluxStat%statIndex_scalarSnowfall = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarSnowfall(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopySunlitPAR); fluxStat%statIndex_scalarCanopySunlitPAR = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopySunlitPAR(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopyShadedPAR); fluxStat%statIndex_scalarCanopyShadedPAR = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopyShadedPAR(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarBelowCanopySolar); fluxStat%statIndex_scalarBelowCanopySolar = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarBelowCanopySolar(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopyAbsorbedSolar); fluxStat%statIndex_scalarCanopyAbsorbedSolar = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopyAbsorbedSolar(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarGroundAbsorbedSolar); fluxStat%statIndex_scalarGroundAbsorbedSolar = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarGroundAbsorbedSolar(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWRadCanopy); fluxStat%statIndex_scalarLWRadCanopy = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWRadCanopy(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWRadGround); fluxStat%statIndex_scalarLWRadGround = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWRadGround(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWRadUbound2Canopy); fluxStat%statIndex_scalarLWRadUbound2Canopy = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWRadUbound2Canopy(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWRadUbound2Ground); fluxStat%statIndex_scalarLWRadUbound2Ground = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWRadUbound2Ground(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWRadUbound2Ubound); fluxStat%statIndex_scalarLWRadUbound2Ubound = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWRadUbound2Ubound(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWRadCanopy2Ubound); fluxStat%statIndex_scalarLWRadCanopy2Ubound = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWRadCanopy2Ubound(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWRadCanopy2Ground); fluxStat%statIndex_scalarLWRadCanopy2Ground = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWRadCanopy2Ground(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWRadCanopy2Canopy); fluxStat%statIndex_scalarLWRadCanopy2Canopy = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWRadCanopy2Canopy(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWRadGround2Ubound); fluxStat%statIndex_scalarLWRadGround2Ubound = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWRadGround2Ubound(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWRadGround2Canopy); fluxStat%statIndex_scalarLWRadGround2Canopy = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWRadGround2Canopy(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWNetCanopy); fluxStat%statIndex_scalarLWNetCanopy = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWNetCanopy(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWNetGround); fluxStat%statIndex_scalarLWNetGround = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWNetGround(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLWNetUbound); fluxStat%statIndex_scalarLWNetUbound = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLWNetUbound(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarEddyDiffusCanopyTop); fluxStat%statIndex_scalarEddyDiffusCanopyTop = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarEddyDiffusCanopyTop(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarFrictionVelocity); fluxStat%statIndex_scalarFrictionVelocity = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarFrictionVelocity(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarWindspdCanopyTop); fluxStat%statIndex_scalarWindspdCanopyTop = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarWindspdCanopyTop(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarWindspdCanopyBottom); fluxStat%statIndex_scalarWindspdCanopyBottom = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarWindspdCanopyBottom(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarGroundResistance); fluxStat%statIndex_scalarGroundResistance = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarGroundResistance(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopyResistance); fluxStat%statIndex_scalarCanopyResistance = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopyResistance(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLeafResistance); fluxStat%statIndex_scalarLeafResistance = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLeafResistance(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarSoilResistance); fluxStat%statIndex_scalarSoilResistance = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarSoilResistance(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarSenHeatTotal); fluxStat%statIndex_scalarSenHeatTotal = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarSenHeatTotal(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarSenHeatCanopy); fluxStat%statIndex_scalarSenHeatCanopy = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarSenHeatCanopy(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarSenHeatGround); fluxStat%statIndex_scalarSenHeatGround = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarSenHeatGround(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLatHeatTotal); fluxStat%statIndex_scalarLatHeatTotal = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLatHeatTotal(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLatHeatCanopyEvap); fluxStat%statIndex_scalarLatHeatCanopyEvap = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLatHeatCanopyEvap(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLatHeatCanopyTrans); fluxStat%statIndex_scalarLatHeatCanopyTrans = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLatHeatCanopyTrans(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarLatHeatGround); fluxStat%statIndex_scalarLatHeatGround = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarLatHeatGround(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopyAdvectiveHeatFlux); fluxStat%statIndex_scalarCanopyAdvectiveHeatFlux = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopyAdvectiveHeatFlux(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarGroundAdvectiveHeatFlux); fluxStat%statIndex_scalarGroundAdvectiveHeatFlux = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarGroundAdvectiveHeatFlux(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopySublimation); fluxStat%statIndex_scalarCanopySublimation = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopySublimation(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarSnowSublimation); fluxStat%statIndex_scalarSnowSublimation = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarSnowSublimation(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarStomResistSunlit); fluxStat%statIndex_scalarStomResistSunlit = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarStomResistSunlit(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarStomResistShaded); fluxStat%statIndex_scalarStomResistShaded = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarStomResistShaded(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarPhotosynthesisSunlit); fluxStat%statIndex_scalarPhotosynthesisSunlit = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarPhotosynthesisSunlit(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarPhotosynthesisShaded); fluxStat%statIndex_scalarPhotosynthesisShaded = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarPhotosynthesisShaded(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopyTranspiration); fluxStat%statIndex_scalarCanopyTranspiration = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopyTranspiration(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopyEvaporation); fluxStat%statIndex_scalarCanopyEvaporation = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopyEvaporation(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarGroundEvaporation); fluxStat%statIndex_scalarGroundEvaporation = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarGroundEvaporation(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarThroughfallSnow); fluxStat%statIndex_scalarThroughfallSnow = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarThroughfallSnow(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarThroughfallRain); fluxStat%statIndex_scalarThroughfallRain = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarThroughfallRain(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopySnowUnloading); fluxStat%statIndex_scalarCanopySnowUnloading = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopySnowUnloading(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopyLiqDrainage); fluxStat%statIndex_scalarCanopyLiqDrainage = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopyLiqDrainage(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarCanopyMeltFreeze); fluxStat%statIndex_scalarCanopyMeltFreeze = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarCanopyMeltFreeze(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarSnowDrainage); fluxStat%statIndex_scalarSnowDrainage = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarSnowDrainage(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarRainPlusMelt); fluxStat%statIndex_scalarRainPlusMelt = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarRainPlusMelt(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarMaxInfilRate); fluxStat%statIndex_scalarMaxInfilRate = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarMaxInfilRate(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarInfiltration); fluxStat%statIndex_scalarInfiltration = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarInfiltration(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarExfiltration); fluxStat%statIndex_scalarExfiltration = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarExfiltration(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarSurfaceRunoff); fluxStat%statIndex_scalarSurfaceRunoff = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarSurfaceRunoff(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarSurfaceRunoff_IE); fluxStat%statIndex_scalarSurfaceRunoff_IE = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarSurfaceRunoff_IE(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarSurfaceRunoff_SE); fluxStat%statIndex_scalarSurfaceRunoff_SE = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarSurfaceRunoff_SE(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarSoilBaseflow); fluxStat%statIndex_scalarSoilBaseflow = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarSoilBaseflow(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarSoilDrainage); fluxStat%statIndex_scalarSoilDrainage = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarSoilDrainage(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarAquiferRecharge); fluxStat%statIndex_scalarAquiferRecharge = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarAquiferRecharge(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarAquiferTranspire); fluxStat%statIndex_scalarAquiferTranspire = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarAquiferTranspire(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarAquiferBaseflow); fluxStat%statIndex_scalarAquiferBaseflow = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarAquiferBaseflow(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarTotalET); fluxStat%statIndex_scalarTotalET = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarTotalET(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarTotalRunoff); fluxStat%statIndex_scalarTotalRunoff = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarTotalRunoff(2*maxVarFreq,nGRU))
  case(iLookFLUX%scalarNetRadiation); fluxStat%statIndex_scalarNetRadiation = statFlux_meta(iVar)%statIndex
    allocate(fluxStat%stat_scalarNetRadiation(2*maxVarFreq,nGRU))
  end select
end do
if (.not. allocated(fluxStat%statIndex_scalarCanairNetNrgFlux)) allocate(fluxStat%statIndex_scalarCanairNetNrgFlux(1))
if (.not. allocated(fluxStat%stat_scalarCanairNetNrgFlux)) allocate(fluxStat%stat_scalarCanairNetNrgFlux(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopyNetNrgFlux)) allocate(fluxStat%statIndex_scalarCanopyNetNrgFlux(1))
if (.not. allocated(fluxStat%stat_scalarCanopyNetNrgFlux)) allocate(fluxStat%stat_scalarCanopyNetNrgFlux(1,1))
if (.not. allocated(fluxStat%statIndex_scalarGroundNetNrgFlux)) allocate(fluxStat%statIndex_scalarGroundNetNrgFlux(1))
if (.not. allocated(fluxStat%stat_scalarGroundNetNrgFlux)) allocate(fluxStat%stat_scalarGroundNetNrgFlux(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopyNetLiqFlux)) allocate(fluxStat%statIndex_scalarCanopyNetLiqFlux(1))
if (.not. allocated(fluxStat%stat_scalarCanopyNetLiqFlux)) allocate(fluxStat%stat_scalarCanopyNetLiqFlux(1,1))
if (.not. allocated(fluxStat%statIndex_scalarRainfall)) allocate(fluxStat%statIndex_scalarRainfall(1))
if (.not. allocated(fluxStat%stat_scalarRainfall)) allocate(fluxStat%stat_scalarRainfall(1,1))
if (.not. allocated(fluxStat%statIndex_scalarSnowfall)) allocate(fluxStat%statIndex_scalarSnowfall(1))
if (.not. allocated(fluxStat%stat_scalarSnowfall)) allocate(fluxStat%stat_scalarSnowfall(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopySunlitPAR)) allocate(fluxStat%statIndex_scalarCanopySunlitPAR(1))
if (.not. allocated(fluxStat%stat_scalarCanopySunlitPAR)) allocate(fluxStat%stat_scalarCanopySunlitPAR(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopyShadedPAR)) allocate(fluxStat%statIndex_scalarCanopyShadedPAR(1))
if (.not. allocated(fluxStat%stat_scalarCanopyShadedPAR)) allocate(fluxStat%stat_scalarCanopyShadedPAR(1,1))
if (.not. allocated(fluxStat%statIndex_scalarBelowCanopySolar)) allocate(fluxStat%statIndex_scalarBelowCanopySolar(1))
if (.not. allocated(fluxStat%stat_scalarBelowCanopySolar)) allocate(fluxStat%stat_scalarBelowCanopySolar(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopyAbsorbedSolar)) allocate(fluxStat%statIndex_scalarCanopyAbsorbedSolar(1))
if (.not. allocated(fluxStat%stat_scalarCanopyAbsorbedSolar)) allocate(fluxStat%stat_scalarCanopyAbsorbedSolar(1,1))
if (.not. allocated(fluxStat%statIndex_scalarGroundAbsorbedSolar)) allocate(fluxStat%statIndex_scalarGroundAbsorbedSolar(1))
if (.not. allocated(fluxStat%stat_scalarGroundAbsorbedSolar)) allocate(fluxStat%stat_scalarGroundAbsorbedSolar(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWRadCanopy)) allocate(fluxStat%statIndex_scalarLWRadCanopy(1))
if (.not. allocated(fluxStat%stat_scalarLWRadCanopy)) allocate(fluxStat%stat_scalarLWRadCanopy(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWRadGround)) allocate(fluxStat%statIndex_scalarLWRadGround(1))
if (.not. allocated(fluxStat%stat_scalarLWRadGround)) allocate(fluxStat%stat_scalarLWRadGround(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWRadUbound2Canopy)) allocate(fluxStat%statIndex_scalarLWRadUbound2Canopy(1))
if (.not. allocated(fluxStat%stat_scalarLWRadUbound2Canopy)) allocate(fluxStat%stat_scalarLWRadUbound2Canopy(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWRadUbound2Ground)) allocate(fluxStat%statIndex_scalarLWRadUbound2Ground(1))
if (.not. allocated(fluxStat%stat_scalarLWRadUbound2Ground)) allocate(fluxStat%stat_scalarLWRadUbound2Ground(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWRadUbound2Ubound)) allocate(fluxStat%statIndex_scalarLWRadUbound2Ubound(1))
if (.not. allocated(fluxStat%stat_scalarLWRadUbound2Ubound)) allocate(fluxStat%stat_scalarLWRadUbound2Ubound(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWRadCanopy2Ubound)) allocate(fluxStat%statIndex_scalarLWRadCanopy2Ubound(1))
if (.not. allocated(fluxStat%stat_scalarLWRadCanopy2Ubound)) allocate(fluxStat%stat_scalarLWRadCanopy2Ubound(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWRadCanopy2Ground)) allocate(fluxStat%statIndex_scalarLWRadCanopy2Ground(1))
if (.not. allocated(fluxStat%stat_scalarLWRadCanopy2Ground)) allocate(fluxStat%stat_scalarLWRadCanopy2Ground(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWRadCanopy2Canopy)) allocate(fluxStat%statIndex_scalarLWRadCanopy2Canopy(1))
if (.not. allocated(fluxStat%stat_scalarLWRadCanopy2Canopy)) allocate(fluxStat%stat_scalarLWRadCanopy2Canopy(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWRadGround2Ubound)) allocate(fluxStat%statIndex_scalarLWRadGround2Ubound(1))
if (.not. allocated(fluxStat%stat_scalarLWRadGround2Ubound)) allocate(fluxStat%stat_scalarLWRadGround2Ubound(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWRadGround2Canopy)) allocate(fluxStat%statIndex_scalarLWRadGround2Canopy(1))
if (.not. allocated(fluxStat%stat_scalarLWRadGround2Canopy)) allocate(fluxStat%stat_scalarLWRadGround2Canopy(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWNetCanopy)) allocate(fluxStat%statIndex_scalarLWNetCanopy(1))
if (.not. allocated(fluxStat%stat_scalarLWNetCanopy)) allocate(fluxStat%stat_scalarLWNetCanopy(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWNetGround)) allocate(fluxStat%statIndex_scalarLWNetGround(1))
if (.not. allocated(fluxStat%stat_scalarLWNetGround)) allocate(fluxStat%stat_scalarLWNetGround(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLWNetUbound)) allocate(fluxStat%statIndex_scalarLWNetUbound(1))
if (.not. allocated(fluxStat%stat_scalarLWNetUbound)) allocate(fluxStat%stat_scalarLWNetUbound(1,1))
if (.not. allocated(fluxStat%statIndex_scalarEddyDiffusCanopyTop)) allocate(fluxStat%statIndex_scalarEddyDiffusCanopyTop(1))
if (.not. allocated(fluxStat%stat_scalarEddyDiffusCanopyTop)) allocate(fluxStat%stat_scalarEddyDiffusCanopyTop(1,1))
if (.not. allocated(fluxStat%statIndex_scalarFrictionVelocity)) allocate(fluxStat%statIndex_scalarFrictionVelocity(1))
if (.not. allocated(fluxStat%stat_scalarFrictionVelocity)) allocate(fluxStat%stat_scalarFrictionVelocity(1,1))
if (.not. allocated(fluxStat%statIndex_scalarWindspdCanopyTop)) allocate(fluxStat%statIndex_scalarWindspdCanopyTop(1))
if (.not. allocated(fluxStat%stat_scalarWindspdCanopyTop)) allocate(fluxStat%stat_scalarWindspdCanopyTop(1,1))
if (.not. allocated(fluxStat%statIndex_scalarWindspdCanopyBottom)) allocate(fluxStat%statIndex_scalarWindspdCanopyBottom(1))
if (.not. allocated(fluxStat%stat_scalarWindspdCanopyBottom)) allocate(fluxStat%stat_scalarWindspdCanopyBottom(1,1))
if (.not. allocated(fluxStat%statIndex_scalarGroundResistance)) allocate(fluxStat%statIndex_scalarGroundResistance(1))
if (.not. allocated(fluxStat%stat_scalarGroundResistance)) allocate(fluxStat%stat_scalarGroundResistance(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopyResistance)) allocate(fluxStat%statIndex_scalarCanopyResistance(1))
if (.not. allocated(fluxStat%stat_scalarCanopyResistance)) allocate(fluxStat%stat_scalarCanopyResistance(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLeafResistance)) allocate(fluxStat%statIndex_scalarLeafResistance(1))
if (.not. allocated(fluxStat%stat_scalarLeafResistance)) allocate(fluxStat%stat_scalarLeafResistance(1,1))
if (.not. allocated(fluxStat%statIndex_scalarSoilResistance)) allocate(fluxStat%statIndex_scalarSoilResistance(1))
if (.not. allocated(fluxStat%stat_scalarSoilResistance)) allocate(fluxStat%stat_scalarSoilResistance(1,1))
if (.not. allocated(fluxStat%statIndex_scalarSenHeatTotal)) allocate(fluxStat%statIndex_scalarSenHeatTotal(1))
if (.not. allocated(fluxStat%stat_scalarSenHeatTotal)) allocate(fluxStat%stat_scalarSenHeatTotal(1,1))
if (.not. allocated(fluxStat%statIndex_scalarSenHeatCanopy)) allocate(fluxStat%statIndex_scalarSenHeatCanopy(1))
if (.not. allocated(fluxStat%stat_scalarSenHeatCanopy)) allocate(fluxStat%stat_scalarSenHeatCanopy(1,1))
if (.not. allocated(fluxStat%statIndex_scalarSenHeatGround)) allocate(fluxStat%statIndex_scalarSenHeatGround(1))
if (.not. allocated(fluxStat%stat_scalarSenHeatGround)) allocate(fluxStat%stat_scalarSenHeatGround(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLatHeatTotal)) allocate(fluxStat%statIndex_scalarLatHeatTotal(1))
if (.not. allocated(fluxStat%stat_scalarLatHeatTotal)) allocate(fluxStat%stat_scalarLatHeatTotal(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLatHeatCanopyEvap)) allocate(fluxStat%statIndex_scalarLatHeatCanopyEvap(1))
if (.not. allocated(fluxStat%stat_scalarLatHeatCanopyEvap)) allocate(fluxStat%stat_scalarLatHeatCanopyEvap(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLatHeatCanopyTrans)) allocate(fluxStat%statIndex_scalarLatHeatCanopyTrans(1))
if (.not. allocated(fluxStat%stat_scalarLatHeatCanopyTrans)) allocate(fluxStat%stat_scalarLatHeatCanopyTrans(1,1))
if (.not. allocated(fluxStat%statIndex_scalarLatHeatGround)) allocate(fluxStat%statIndex_scalarLatHeatGround(1))
if (.not. allocated(fluxStat%stat_scalarLatHeatGround)) allocate(fluxStat%stat_scalarLatHeatGround(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopyAdvectiveHeatFlux)) allocate(fluxStat%statIndex_scalarCanopyAdvectiveHeatFlux(1))
if (.not. allocated(fluxStat%stat_scalarCanopyAdvectiveHeatFlux)) allocate(fluxStat%stat_scalarCanopyAdvectiveHeatFlux(1,1))
if (.not. allocated(fluxStat%statIndex_scalarGroundAdvectiveHeatFlux)) allocate(fluxStat%statIndex_scalarGroundAdvectiveHeatFlux(1))
if (.not. allocated(fluxStat%stat_scalarGroundAdvectiveHeatFlux)) allocate(fluxStat%stat_scalarGroundAdvectiveHeatFlux(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopySublimation)) allocate(fluxStat%statIndex_scalarCanopySublimation(1))
if (.not. allocated(fluxStat%stat_scalarCanopySublimation)) allocate(fluxStat%stat_scalarCanopySublimation(1,1))
if (.not. allocated(fluxStat%statIndex_scalarSnowSublimation)) allocate(fluxStat%statIndex_scalarSnowSublimation(1))
if (.not. allocated(fluxStat%stat_scalarSnowSublimation)) allocate(fluxStat%stat_scalarSnowSublimation(1,1))
if (.not. allocated(fluxStat%statIndex_scalarStomResistSunlit)) allocate(fluxStat%statIndex_scalarStomResistSunlit(1))
if (.not. allocated(fluxStat%stat_scalarStomResistSunlit)) allocate(fluxStat%stat_scalarStomResistSunlit(1,1))
if (.not. allocated(fluxStat%statIndex_scalarStomResistShaded)) allocate(fluxStat%statIndex_scalarStomResistShaded(1))
if (.not. allocated(fluxStat%stat_scalarStomResistShaded)) allocate(fluxStat%stat_scalarStomResistShaded(1,1))
if (.not. allocated(fluxStat%statIndex_scalarPhotosynthesisSunlit)) allocate(fluxStat%statIndex_scalarPhotosynthesisSunlit(1))
if (.not. allocated(fluxStat%stat_scalarPhotosynthesisSunlit)) allocate(fluxStat%stat_scalarPhotosynthesisSunlit(1,1))
if (.not. allocated(fluxStat%statIndex_scalarPhotosynthesisShaded)) allocate(fluxStat%statIndex_scalarPhotosynthesisShaded(1))
if (.not. allocated(fluxStat%stat_scalarPhotosynthesisShaded)) allocate(fluxStat%stat_scalarPhotosynthesisShaded(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopyTranspiration)) allocate(fluxStat%statIndex_scalarCanopyTranspiration(1))
if (.not. allocated(fluxStat%stat_scalarCanopyTranspiration)) allocate(fluxStat%stat_scalarCanopyTranspiration(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopyEvaporation)) allocate(fluxStat%statIndex_scalarCanopyEvaporation(1))
if (.not. allocated(fluxStat%stat_scalarCanopyEvaporation)) allocate(fluxStat%stat_scalarCanopyEvaporation(1,1))
if (.not. allocated(fluxStat%statIndex_scalarGroundEvaporation)) allocate(fluxStat%statIndex_scalarGroundEvaporation(1))
if (.not. allocated(fluxStat%stat_scalarGroundEvaporation)) allocate(fluxStat%stat_scalarGroundEvaporation(1,1))
if (.not. allocated(fluxStat%statIndex_scalarThroughfallSnow)) allocate(fluxStat%statIndex_scalarThroughfallSnow(1))
if (.not. allocated(fluxStat%stat_scalarThroughfallSnow)) allocate(fluxStat%stat_scalarThroughfallSnow(1,1))
if (.not. allocated(fluxStat%statIndex_scalarThroughfallRain)) allocate(fluxStat%statIndex_scalarThroughfallRain(1))
if (.not. allocated(fluxStat%stat_scalarThroughfallRain)) allocate(fluxStat%stat_scalarThroughfallRain(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopySnowUnloading)) allocate(fluxStat%statIndex_scalarCanopySnowUnloading(1))
if (.not. allocated(fluxStat%stat_scalarCanopySnowUnloading)) allocate(fluxStat%stat_scalarCanopySnowUnloading(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopyLiqDrainage)) allocate(fluxStat%statIndex_scalarCanopyLiqDrainage(1))
if (.not. allocated(fluxStat%stat_scalarCanopyLiqDrainage)) allocate(fluxStat%stat_scalarCanopyLiqDrainage(1,1))
if (.not. allocated(fluxStat%statIndex_scalarCanopyMeltFreeze)) allocate(fluxStat%statIndex_scalarCanopyMeltFreeze(1))
if (.not. allocated(fluxStat%stat_scalarCanopyMeltFreeze)) allocate(fluxStat%stat_scalarCanopyMeltFreeze(1,1))
if (.not. allocated(fluxStat%statIndex_scalarSnowDrainage)) allocate(fluxStat%statIndex_scalarSnowDrainage(1))
if (.not. allocated(fluxStat%stat_scalarSnowDrainage)) allocate(fluxStat%stat_scalarSnowDrainage(1,1))
if (.not. allocated(fluxStat%statIndex_scalarRainPlusMelt)) allocate(fluxStat%statIndex_scalarRainPlusMelt(1))
if (.not. allocated(fluxStat%stat_scalarRainPlusMelt)) allocate(fluxStat%stat_scalarRainPlusMelt(1,1))
if (.not. allocated(fluxStat%statIndex_scalarMaxInfilRate)) allocate(fluxStat%statIndex_scalarMaxInfilRate(1))
if (.not. allocated(fluxStat%stat_scalarMaxInfilRate)) allocate(fluxStat%stat_scalarMaxInfilRate(1,1))
if (.not. allocated(fluxStat%statIndex_scalarInfiltration)) allocate(fluxStat%statIndex_scalarInfiltration(1))
if (.not. allocated(fluxStat%stat_scalarInfiltration)) allocate(fluxStat%stat_scalarInfiltration(1,1))
if (.not. allocated(fluxStat%statIndex_scalarExfiltration)) allocate(fluxStat%statIndex_scalarExfiltration(1))
if (.not. allocated(fluxStat%stat_scalarExfiltration)) allocate(fluxStat%stat_scalarExfiltration(1,1))
if (.not. allocated(fluxStat%statIndex_scalarSurfaceRunoff)) allocate(fluxStat%statIndex_scalarSurfaceRunoff(1))
if (.not. allocated(fluxStat%stat_scalarSurfaceRunoff)) allocate(fluxStat%stat_scalarSurfaceRunoff(1,1))
if (.not. allocated(fluxStat%statIndex_scalarSurfaceRunoff_IE)) allocate(fluxStat%statIndex_scalarSurfaceRunoff_IE(1))
if (.not. allocated(fluxStat%stat_scalarSurfaceRunoff_IE)) allocate(fluxStat%stat_scalarSurfaceRunoff_IE(1,1))
if (.not. allocated(fluxStat%statIndex_scalarSurfaceRunoff_SE)) allocate(fluxStat%statIndex_scalarSurfaceRunoff_SE(1))
if (.not. allocated(fluxStat%stat_scalarSurfaceRunoff_SE)) allocate(fluxStat%stat_scalarSurfaceRunoff_SE(1,1))
if (.not. allocated(fluxStat%statIndex_scalarSoilBaseflow)) allocate(fluxStat%statIndex_scalarSoilBaseflow(1))
if (.not. allocated(fluxStat%stat_scalarSoilBaseflow)) allocate(fluxStat%stat_scalarSoilBaseflow(1,1))
if (.not. allocated(fluxStat%statIndex_scalarSoilDrainage)) allocate(fluxStat%statIndex_scalarSoilDrainage(1))
if (.not. allocated(fluxStat%stat_scalarSoilDrainage)) allocate(fluxStat%stat_scalarSoilDrainage(1,1))
if (.not. allocated(fluxStat%statIndex_scalarAquiferRecharge)) allocate(fluxStat%statIndex_scalarAquiferRecharge(1))
if (.not. allocated(fluxStat%stat_scalarAquiferRecharge)) allocate(fluxStat%stat_scalarAquiferRecharge(1,1))
if (.not. allocated(fluxStat%statIndex_scalarAquiferTranspire)) allocate(fluxStat%statIndex_scalarAquiferTranspire(1))
if (.not. allocated(fluxStat%stat_scalarAquiferTranspire)) allocate(fluxStat%stat_scalarAquiferTranspire(1,1))
if (.not. allocated(fluxStat%statIndex_scalarAquiferBaseflow)) allocate(fluxStat%statIndex_scalarAquiferBaseflow(1))
if (.not. allocated(fluxStat%stat_scalarAquiferBaseflow)) allocate(fluxStat%stat_scalarAquiferBaseflow(1,1))
if (.not. allocated(fluxStat%statIndex_scalarTotalET)) allocate(fluxStat%statIndex_scalarTotalET(1))
if (.not. allocated(fluxStat%stat_scalarTotalET)) allocate(fluxStat%stat_scalarTotalET(1,1))
if (.not. allocated(fluxStat%statIndex_scalarTotalRunoff)) allocate(fluxStat%statIndex_scalarTotalRunoff(1))
if (.not. allocated(fluxStat%stat_scalarTotalRunoff)) allocate(fluxStat%stat_scalarTotalRunoff(1,1))
if (.not. allocated(fluxStat%statIndex_scalarNetRadiation)) allocate(fluxStat%statIndex_scalarNetRadiation(1))
if (.not. allocated(fluxStat%stat_scalarNetRadiation)) allocate(fluxStat%stat_scalarNetRadiation(1,1))
end subroutine

subroutine allocate_device_diag_stat(diagStat,nGRU)
            USE globalData,only:statDiag_meta             ! child metadata for stats
    USE var_lookup,only:maxvarFreq             ! allocation dimension (output frequency)

  implicit none
  type(diag_stat_device) :: diagStat
      integer(i4b) :: nGRU, iGRU, iVar, pVar

do iVar = 1,size(statDiag_meta)
  pVar = statDiag_meta(iVar)%ixParent
  select case (pVar)

  case(iLookDIAG%scalarCanopyDepth); diagStat%statIndex_scalarCanopyDepth = statDiag_meta(iVar)%statIndex
    allocate(diagStat%stat_scalarCanopyDepth(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarBulkVolHeatCapVeg); diagStat%statIndex_scalarBulkVolHeatCapVeg = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarBulkVolHeatCapVeg(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarCanopyCm); diagStat%statIndex_scalarCanopyCm = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarCanopyCm(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarCanopyEmissivity); diagStat%statIndex_scalarCanopyEmissivity = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarCanopyEmissivity(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarRootZoneTemp); diagStat%statIndex_scalarRootZoneTemp = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarRootZoneTemp(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarLAI); diagStat%statIndex_scalarLAI = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarLAI(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarSAI); diagStat%statIndex_scalarSAI = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarSAI(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarExposedLAI); diagStat%statIndex_scalarExposedLAI = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarExposedLAI(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarExposedSAI); diagStat%statIndex_scalarExposedSAI = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarExposedSAI(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarAdjMeasHeight); diagStat%statIndex_scalarAdjMeasHeight = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarAdjMeasHeight(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarCanopyIceMax); diagStat%statIndex_scalarCanopyIceMax = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarCanopyIceMax(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarCanopyLiqMax); diagStat%statIndex_scalarCanopyLiqMax = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarCanopyLiqMax(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarGrowingSeasonIndex); diagStat%statIndex_scalarGrowingSeasonIndex = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarGrowingSeasonIndex(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarVolHtCap_air); diagStat%statIndex_scalarVolHtCap_air = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarVolHtCap_air(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarVolHtCap_ice); diagStat%statIndex_scalarVolHtCap_ice = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarVolHtCap_ice(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarVolHtCap_soil); diagStat%statIndex_scalarVolHtCap_soil = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarVolHtCap_soil(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarVolHtCap_water); diagStat%statIndex_scalarVolHtCap_water = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarVolHtCap_water(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarLambda_drysoil); diagStat%statIndex_scalarLambda_drysoil = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarLambda_drysoil(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarLambda_wetsoil); diagStat%statIndex_scalarLambda_wetsoil = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarLambda_wetsoil(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarCanopyEnthTemp); diagStat%statIndex_scalarCanopyEnthTemp = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarCanopyEnthTemp(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarTotalSoilEnthalpy); diagStat%statIndex_scalarTotalSoilEnthalpy = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarTotalSoilEnthalpy(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarTotalSnowEnthalpy); diagStat%statIndex_scalarTotalSnowEnthalpy = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarTotalSnowEnthalpy(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarVPair); diagStat%statIndex_scalarVPair = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarVPair(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarVP_CanopyAir); diagStat%statIndex_scalarVP_CanopyAir = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarVP_CanopyAir(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarTwetbulb); diagStat%statIndex_scalarTwetbulb = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarTwetbulb(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarSnowfallTemp); diagStat%statIndex_scalarSnowfallTemp = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarSnowfallTemp(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarNewSnowDensity); diagStat%statIndex_scalarNewSnowDensity = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarNewSnowDensity(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarO2air); diagStat%statIndex_scalarO2air = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarO2air(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarCO2air); diagStat%statIndex_scalarCO2air = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarCO2air(2*maxVarFreq,nGRU))
case(iLookDIAG%windspd_x); diagStat%statIndex_windspd_x = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_windspd_x(2*maxVarFreq,nGRU))
case(iLookDIAG%windspd_y); diagStat%statIndex_windspd_y = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_windspd_y(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarCosZenith); diagStat%statIndex_scalarCosZenith = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarCosZenith(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarFractionDirect); diagStat%statIndex_scalarFractionDirect = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarFractionDirect(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarCanopySunlitFraction); diagStat%statIndex_scalarCanopySunlitFraction = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarCanopySunlitFraction(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarCanopySunlitLAI); diagStat%statIndex_scalarCanopySunlitLAI = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarCanopySunlitLAI(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarCanopyShadedLAI); diagStat%statIndex_scalarCanopyShadedLAI = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarCanopyShadedLAI(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarGroundAlbedo); diagStat%statIndex_scalarGroundAlbedo = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarGroundAlbedo(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarLatHeatSubVapCanopy); diagStat%statIndex_scalarLatHeatSubVapCanopy = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarLatHeatSubVapCanopy(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarLatHeatSubVapGround); diagStat%statIndex_scalarLatHeatSubVapGround = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarLatHeatSubVapGround(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarSatVP_CanopyTemp); diagStat%statIndex_scalarSatVP_CanopyTemp = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarSatVP_CanopyTemp(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarSatVP_GroundTemp); diagStat%statIndex_scalarSatVP_GroundTemp = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarSatVP_GroundTemp(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarZ0Canopy); diagStat%statIndex_scalarZ0Canopy = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarZ0Canopy(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarWindReductionFactor); diagStat%statIndex_scalarWindReductionFactor = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarWindReductionFactor(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarZeroPlaneDisplacement); diagStat%statIndex_scalarZeroPlaneDisplacement = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarZeroPlaneDisplacement(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarRiBulkCanopy); diagStat%statIndex_scalarRiBulkCanopy = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarRiBulkCanopy(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarRiBulkGround); diagStat%statIndex_scalarRiBulkGround = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarRiBulkGround(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarCanopyStabilityCorrection); diagStat%statIndex_scalarCanopyStabilityCorrection = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarCanopyStabilityCorrection(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarGroundStabilityCorrection); diagStat%statIndex_scalarGroundStabilityCorrection = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarGroundStabilityCorrection(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarIntercellularCO2Sunlit); diagStat%statIndex_scalarIntercellularCO2Sunlit = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarIntercellularCO2Sunlit(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarIntercellularCO2Shaded); diagStat%statIndex_scalarIntercellularCO2Shaded = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarIntercellularCO2Shaded(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarTranspireLim); diagStat%statIndex_scalarTranspireLim = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarTranspireLim(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarTranspireLimAqfr); diagStat%statIndex_scalarTranspireLimAqfr = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarTranspireLimAqfr(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarFoliageNitrogenFactor); diagStat%statIndex_scalarFoliageNitrogenFactor = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarFoliageNitrogenFactor(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarSoilRelHumidity); diagStat%statIndex_scalarSoilRelHumidity = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarSoilRelHumidity(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarAquiferRootFrac); diagStat%statIndex_scalarAquiferRootFrac = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarAquiferRootFrac(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarFracLiqVeg); diagStat%statIndex_scalarFracLiqVeg = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarFracLiqVeg(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarCanopyWetFraction); diagStat%statIndex_scalarCanopyWetFraction = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarCanopyWetFraction(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarSnowAge); diagStat%statIndex_scalarSnowAge = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarSnowAge(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarGroundSnowFraction); diagStat%statIndex_scalarGroundSnowFraction = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarGroundSnowFraction(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarInfilArea); diagStat%statIndex_scalarInfilArea = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarInfilArea(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarFrozenArea); diagStat%statIndex_scalarFrozenArea = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarFrozenArea(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarSoilControl); diagStat%statIndex_scalarSoilControl = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarSoilControl(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarSoilCompress); diagStat%statIndex_scalarSoilCompress = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarSoilCompress(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarTotalSoilLiq); diagStat%statIndex_scalarTotalSoilLiq = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarTotalSoilLiq(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarTotalSoilIce); diagStat%statIndex_scalarTotalSoilIce = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarTotalSoilIce(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarTotalSoilWat); diagStat%statIndex_scalarTotalSoilWat = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarTotalSoilWat(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarKappa); diagStat%statIndex_scalarKappa = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarKappa(2*maxVarFreq,nGRU))
case(iLookDIAG%scalarVolLatHt_fus); diagStat%statIndex_scalarVolLatHt_fus = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_scalarVolLatHt_fus(2*maxVarFreq,nGRU))
case(iLookDIAG%numFluxCalls); diagStat%statIndex_numFluxCalls = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_numFluxCalls(2*maxVarFreq,nGRU))
case(iLookDIAG%wallClockTime); diagStat%statIndex_wallClockTime = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_wallClockTime(2*maxVarFreq,nGRU))
case(iLookDIAG%meanStepSize); diagStat%statIndex_meanStepSize = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_meanStepSize(2*maxVarFreq,nGRU))
case(iLookDIAG%balanceCasNrg); diagStat%statIndex_balanceCasNrg = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_balanceCasNrg(2*maxVarFreq,nGRU))
case(iLookDIAG%balanceVegNrg); diagStat%statIndex_balanceVegNrg = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_balanceVegNrg(2*maxVarFreq,nGRU))
case(iLookDIAG%balanceSnowNrg); diagStat%statIndex_balanceSnowNrg = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_balanceSnowNrg(2*maxVarFreq,nGRU))
case(iLookDIAG%balanceSoilNrg); diagStat%statIndex_balanceSoilNrg = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_balanceSoilNrg(2*maxVarFreq,nGRU))
case(iLookDIAG%balanceVegMass); diagStat%statIndex_balanceVegMass = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_balanceVegMass(2*maxVarFreq,nGRU))
case(iLookDIAG%balanceSnowMass); diagStat%statIndex_balanceSnowMass = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_balanceSnowMass(2*maxVarFreq,nGRU))
case(iLookDIAG%balanceSoilMass); diagStat%statIndex_balanceSoilMass = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_balanceSoilMass(2*maxVarFreq,nGRU))
case(iLookDIAG%balanceAqMass); diagStat%statIndex_balanceAqMass = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_balanceAqMass(2*maxVarFreq,nGRU))
case(iLookDIAG%numSteps); diagStat%statIndex_numSteps = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_numSteps(2*maxVarFreq,nGRU))
case(iLookDIAG%numResEvals); diagStat%statIndex_numResEvals = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_numResEvals(2*maxVarFreq,nGRU))
case(iLookDIAG%numLinSolvSetups); diagStat%statIndex_numLinSolvSetups = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_numLinSolvSetups(2*maxVarFreq,nGRU))
case(iLookDIAG%numErrTestFails); diagStat%statIndex_numErrTestFails = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_numErrTestFails(2*maxVarFreq,nGRU))
case(iLookDIAG%kLast); diagStat%statIndex_kLast = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_kLast(2*maxVarFreq,nGRU))
case(iLookDIAG%kCur); diagStat%statIndex_kCur = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_kCur(2*maxVarFreq,nGRU))
case(iLookDIAG%hInitUsed); diagStat%statIndex_hInitUsed = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_hInitUsed(2*maxVarFreq,nGRU))
case(iLookDIAG%hLast); diagStat%statIndex_hLast = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_hLast(2*maxVarFreq,nGRU))
case(iLookDIAG%hCur); diagStat%statIndex_hCur = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_hCur(2*maxVarFreq,nGRU))
case(iLookDIAG%tCur); diagStat%statIndex_tCur = statDiag_meta(iVar)%statIndex
  allocate(diagStat%stat_tCur(2*maxVarFreq,nGRU))
end select
end do

if (.not.allocated(diagStat%statIndex_scalarCanopyDepth)) allocate(diagStat%statIndex_scalarCanopyDepth(1))
if (.not.allocated(diagStat%stat_scalarCanopyDepth)) allocate(diagStat%stat_scalarCanopyDepth(1,1))
if (.not.allocated(diagStat%statIndex_scalarBulkVolHeatCapVeg)) allocate(diagStat%statIndex_scalarBulkVolHeatCapVeg(1))
if (.not.allocated(diagStat%stat_scalarBulkVolHeatCapVeg)) allocate(diagStat%stat_scalarBulkVolHeatCapVeg(1,1))
if (.not.allocated(diagStat%statIndex_scalarCanopyCm)) allocate(diagStat%statIndex_scalarCanopyCm(1))
if (.not.allocated(diagStat%stat_scalarCanopyCm)) allocate(diagStat%stat_scalarCanopyCm(1,1))
if (.not.allocated(diagStat%statIndex_scalarCanopyEmissivity)) allocate(diagStat%statIndex_scalarCanopyEmissivity(1))
if (.not.allocated(diagStat%stat_scalarCanopyEmissivity)) allocate(diagStat%stat_scalarCanopyEmissivity(1,1))
if (.not.allocated(diagStat%statIndex_scalarRootZoneTemp)) allocate(diagStat%statIndex_scalarRootZoneTemp(1))
if (.not.allocated(diagStat%stat_scalarRootZoneTemp)) allocate(diagStat%stat_scalarRootZoneTemp(1,1))
if (.not.allocated(diagStat%statIndex_scalarLAI)) allocate(diagStat%statIndex_scalarLAI(1))
if (.not.allocated(diagStat%stat_scalarLAI)) allocate(diagStat%stat_scalarLAI(1,1))
if (.not.allocated(diagStat%statIndex_scalarSAI)) allocate(diagStat%statIndex_scalarSAI(1))
if (.not.allocated(diagStat%stat_scalarSAI)) allocate(diagStat%stat_scalarSAI(1,1))
if (.not.allocated(diagStat%statIndex_scalarExposedLAI)) allocate(diagStat%statIndex_scalarExposedLAI(1))
if (.not.allocated(diagStat%stat_scalarExposedLAI)) allocate(diagStat%stat_scalarExposedLAI(1,1))
if (.not.allocated(diagStat%statIndex_scalarExposedSAI)) allocate(diagStat%statIndex_scalarExposedSAI(1))
if (.not.allocated(diagStat%stat_scalarExposedSAI)) allocate(diagStat%stat_scalarExposedSAI(1,1))
if (.not.allocated(diagStat%statIndex_scalarAdjMeasHeight)) allocate(diagStat%statIndex_scalarAdjMeasHeight(1))
if (.not.allocated(diagStat%stat_scalarAdjMeasHeight)) allocate(diagStat%stat_scalarAdjMeasHeight(1,1))
if (.not.allocated(diagStat%statIndex_scalarCanopyIceMax)) allocate(diagStat%statIndex_scalarCanopyIceMax(1))
if (.not.allocated(diagStat%stat_scalarCanopyIceMax)) allocate(diagStat%stat_scalarCanopyIceMax(1,1))
if (.not.allocated(diagStat%statIndex_scalarCanopyLiqMax)) allocate(diagStat%statIndex_scalarCanopyLiqMax(1))
if (.not.allocated(diagStat%stat_scalarCanopyLiqMax)) allocate(diagStat%stat_scalarCanopyLiqMax(1,1))
if (.not.allocated(diagStat%statIndex_scalarGrowingSeasonIndex)) allocate(diagStat%statIndex_scalarGrowingSeasonIndex(1))
if (.not.allocated(diagStat%stat_scalarGrowingSeasonIndex)) allocate(diagStat%stat_scalarGrowingSeasonIndex(1,1))
if (.not.allocated(diagStat%statIndex_scalarVolHtCap_air)) allocate(diagStat%statIndex_scalarVolHtCap_air(1))
if (.not.allocated(diagStat%stat_scalarVolHtCap_air)) allocate(diagStat%stat_scalarVolHtCap_air(1,1))
if (.not.allocated(diagStat%statIndex_scalarVolHtCap_ice)) allocate(diagStat%statIndex_scalarVolHtCap_ice(1))
if (.not.allocated(diagStat%stat_scalarVolHtCap_ice)) allocate(diagStat%stat_scalarVolHtCap_ice(1,1))
if (.not.allocated(diagStat%statIndex_scalarVolHtCap_soil)) allocate(diagStat%statIndex_scalarVolHtCap_soil(1))
if (.not.allocated(diagStat%stat_scalarVolHtCap_soil)) allocate(diagStat%stat_scalarVolHtCap_soil(1,1))
if (.not.allocated(diagStat%statIndex_scalarVolHtCap_water)) allocate(diagStat%statIndex_scalarVolHtCap_water(1))
if (.not.allocated(diagStat%stat_scalarVolHtCap_water)) allocate(diagStat%stat_scalarVolHtCap_water(1,1))
if (.not.allocated(diagStat%statIndex_scalarLambda_drysoil)) allocate(diagStat%statIndex_scalarLambda_drysoil(1))
if (.not.allocated(diagStat%stat_scalarLambda_drysoil)) allocate(diagStat%stat_scalarLambda_drysoil(1,1))
if (.not.allocated(diagStat%statIndex_scalarLambda_wetsoil)) allocate(diagStat%statIndex_scalarLambda_wetsoil(1))
if (.not.allocated(diagStat%stat_scalarLambda_wetsoil)) allocate(diagStat%stat_scalarLambda_wetsoil(1,1))
if (.not.allocated(diagStat%statIndex_scalarCanopyEnthTemp)) allocate(diagStat%statIndex_scalarCanopyEnthTemp(1))
if (.not.allocated(diagStat%stat_scalarCanopyEnthTemp)) allocate(diagStat%stat_scalarCanopyEnthTemp(1,1))
if (.not.allocated(diagStat%statIndex_scalarTotalSoilEnthalpy)) allocate(diagStat%statIndex_scalarTotalSoilEnthalpy(1))
if (.not.allocated(diagStat%stat_scalarTotalSoilEnthalpy)) allocate(diagStat%stat_scalarTotalSoilEnthalpy(1,1))
if (.not.allocated(diagStat%statIndex_scalarTotalSnowEnthalpy)) allocate(diagStat%statIndex_scalarTotalSnowEnthalpy(1))
if (.not.allocated(diagStat%stat_scalarTotalSnowEnthalpy)) allocate(diagStat%stat_scalarTotalSnowEnthalpy(1,1))
if (.not.allocated(diagStat%statIndex_scalarVPair)) allocate(diagStat%statIndex_scalarVPair(1))
if (.not.allocated(diagStat%stat_scalarVPair)) allocate(diagStat%stat_scalarVPair(1,1))
if (.not.allocated(diagStat%statIndex_scalarVP_CanopyAir)) allocate(diagStat%statIndex_scalarVP_CanopyAir(1))
if (.not.allocated(diagStat%stat_scalarVP_CanopyAir)) allocate(diagStat%stat_scalarVP_CanopyAir(1,1))
if (.not.allocated(diagStat%statIndex_scalarTwetbulb)) allocate(diagStat%statIndex_scalarTwetbulb(1))
if (.not.allocated(diagStat%stat_scalarTwetbulb)) allocate(diagStat%stat_scalarTwetbulb(1,1))
if (.not.allocated(diagStat%statIndex_scalarSnowfallTemp)) allocate(diagStat%statIndex_scalarSnowfallTemp(1))
if (.not.allocated(diagStat%stat_scalarSnowfallTemp)) allocate(diagStat%stat_scalarSnowfallTemp(1,1))
if (.not.allocated(diagStat%statIndex_scalarNewSnowDensity)) allocate(diagStat%statIndex_scalarNewSnowDensity(1))
if (.not.allocated(diagStat%stat_scalarNewSnowDensity)) allocate(diagStat%stat_scalarNewSnowDensity(1,1))
if (.not.allocated(diagStat%statIndex_scalarO2air)) allocate(diagStat%statIndex_scalarO2air(1))
if (.not.allocated(diagStat%stat_scalarO2air)) allocate(diagStat%stat_scalarO2air(1,1))
if (.not.allocated(diagStat%statIndex_scalarCO2air)) allocate(diagStat%statIndex_scalarCO2air(1))
if (.not.allocated(diagStat%stat_scalarCO2air)) allocate(diagStat%stat_scalarCO2air(1,1))
if (.not.allocated(diagStat%statIndex_windspd_x)) allocate(diagStat%statIndex_windspd_x(1))
if (.not.allocated(diagStat%stat_windspd_x)) allocate(diagStat%stat_windspd_x(1,1))
if (.not.allocated(diagStat%statIndex_windspd_y)) allocate(diagStat%statIndex_windspd_y(1))
if (.not.allocated(diagStat%stat_windspd_y)) allocate(diagStat%stat_windspd_y(1,1))
if (.not.allocated(diagStat%statIndex_scalarCosZenith)) allocate(diagStat%statIndex_scalarCosZenith(1))
if (.not.allocated(diagStat%stat_scalarCosZenith)) allocate(diagStat%stat_scalarCosZenith(1,1))
if (.not.allocated(diagStat%statIndex_scalarFractionDirect)) allocate(diagStat%statIndex_scalarFractionDirect(1))
if (.not.allocated(diagStat%stat_scalarFractionDirect)) allocate(diagStat%stat_scalarFractionDirect(1,1))
if (.not.allocated(diagStat%statIndex_scalarCanopySunlitFraction)) allocate(diagStat%statIndex_scalarCanopySunlitFraction(1))
if (.not.allocated(diagStat%stat_scalarCanopySunlitFraction)) allocate(diagStat%stat_scalarCanopySunlitFraction(1,1))
if (.not.allocated(diagStat%statIndex_scalarCanopySunlitLAI)) allocate(diagStat%statIndex_scalarCanopySunlitLAI(1))
if (.not.allocated(diagStat%stat_scalarCanopySunlitLAI)) allocate(diagStat%stat_scalarCanopySunlitLAI(1,1))
if (.not.allocated(diagStat%statIndex_scalarCanopyShadedLAI)) allocate(diagStat%statIndex_scalarCanopyShadedLAI(1))
if (.not.allocated(diagStat%stat_scalarCanopyShadedLAI)) allocate(diagStat%stat_scalarCanopyShadedLAI(1,1))
if (.not.allocated(diagStat%statIndex_scalarGroundAlbedo)) allocate(diagStat%statIndex_scalarGroundAlbedo(1))
if (.not.allocated(diagStat%stat_scalarGroundAlbedo)) allocate(diagStat%stat_scalarGroundAlbedo(1,1))
if (.not.allocated(diagStat%statIndex_scalarLatHeatSubVapCanopy)) allocate(diagStat%statIndex_scalarLatHeatSubVapCanopy(1))
if (.not.allocated(diagStat%stat_scalarLatHeatSubVapCanopy)) allocate(diagStat%stat_scalarLatHeatSubVapCanopy(1,1))
if (.not.allocated(diagStat%statIndex_scalarLatHeatSubVapGround)) allocate(diagStat%statIndex_scalarLatHeatSubVapGround(1))
if (.not.allocated(diagStat%stat_scalarLatHeatSubVapGround)) allocate(diagStat%stat_scalarLatHeatSubVapGround(1,1))
if (.not.allocated(diagStat%statIndex_scalarSatVP_CanopyTemp)) allocate(diagStat%statIndex_scalarSatVP_CanopyTemp(1))
if (.not.allocated(diagStat%stat_scalarSatVP_CanopyTemp)) allocate(diagStat%stat_scalarSatVP_CanopyTemp(1,1))
if (.not.allocated(diagStat%statIndex_scalarSatVP_GroundTemp)) allocate(diagStat%statIndex_scalarSatVP_GroundTemp(1))
if (.not.allocated(diagStat%stat_scalarSatVP_GroundTemp)) allocate(diagStat%stat_scalarSatVP_GroundTemp(1,1))
if (.not.allocated(diagStat%statIndex_scalarZ0Canopy)) allocate(diagStat%statIndex_scalarZ0Canopy(1))
if (.not.allocated(diagStat%stat_scalarZ0Canopy)) allocate(diagStat%stat_scalarZ0Canopy(1,1))
if (.not.allocated(diagStat%statIndex_scalarWindReductionFactor)) allocate(diagStat%statIndex_scalarWindReductionFactor(1))
if (.not.allocated(diagStat%stat_scalarWindReductionFactor)) allocate(diagStat%stat_scalarWindReductionFactor(1,1))
if (.not.allocated(diagStat%statIndex_scalarZeroPlaneDisplacement)) allocate(diagStat%statIndex_scalarZeroPlaneDisplacement(1))
if (.not.allocated(diagStat%stat_scalarZeroPlaneDisplacement)) allocate(diagStat%stat_scalarZeroPlaneDisplacement(1,1))
if (.not.allocated(diagStat%statIndex_scalarRiBulkCanopy)) allocate(diagStat%statIndex_scalarRiBulkCanopy(1))
if (.not.allocated(diagStat%stat_scalarRiBulkCanopy)) allocate(diagStat%stat_scalarRiBulkCanopy(1,1))
if (.not.allocated(diagStat%statIndex_scalarRiBulkGround)) allocate(diagStat%statIndex_scalarRiBulkGround(1))
if (.not.allocated(diagStat%stat_scalarRiBulkGround)) allocate(diagStat%stat_scalarRiBulkGround(1,1))
if (.not.allocated(diagStat%statIndex_scalarCanopyStabilityCorrection)) allocate(diagStat%statIndex_scalarCanopyStabilityCorrection(1))
if (.not.allocated(diagStat%stat_scalarCanopyStabilityCorrection)) allocate(diagStat%stat_scalarCanopyStabilityCorrection(1,1))
if (.not.allocated(diagStat%statIndex_scalarGroundStabilityCorrection)) allocate(diagStat%statIndex_scalarGroundStabilityCorrection(1))
if (.not.allocated(diagStat%stat_scalarGroundStabilityCorrection)) allocate(diagStat%stat_scalarGroundStabilityCorrection(1,1))
if (.not.allocated(diagStat%statIndex_scalarIntercellularCO2Sunlit)) allocate(diagStat%statIndex_scalarIntercellularCO2Sunlit(1))
if (.not.allocated(diagStat%stat_scalarIntercellularCO2Sunlit)) allocate(diagStat%stat_scalarIntercellularCO2Sunlit(1,1))
if (.not.allocated(diagStat%statIndex_scalarIntercellularCO2Shaded)) allocate(diagStat%statIndex_scalarIntercellularCO2Shaded(1))
if (.not.allocated(diagStat%stat_scalarIntercellularCO2Shaded)) allocate(diagStat%stat_scalarIntercellularCO2Shaded(1,1))
if (.not.allocated(diagStat%statIndex_scalarTranspireLim)) allocate(diagStat%statIndex_scalarTranspireLim(1))
if (.not.allocated(diagStat%stat_scalarTranspireLim)) allocate(diagStat%stat_scalarTranspireLim(1,1))
if (.not.allocated(diagStat%statIndex_scalarTranspireLimAqfr)) allocate(diagStat%statIndex_scalarTranspireLimAqfr(1))
if (.not.allocated(diagStat%stat_scalarTranspireLimAqfr)) allocate(diagStat%stat_scalarTranspireLimAqfr(1,1))
if (.not.allocated(diagStat%statIndex_scalarFoliageNitrogenFactor)) allocate(diagStat%statIndex_scalarFoliageNitrogenFactor(1))
if (.not.allocated(diagStat%stat_scalarFoliageNitrogenFactor)) allocate(diagStat%stat_scalarFoliageNitrogenFactor(1,1))
if (.not.allocated(diagStat%statIndex_scalarSoilRelHumidity)) allocate(diagStat%statIndex_scalarSoilRelHumidity(1))
if (.not.allocated(diagStat%stat_scalarSoilRelHumidity)) allocate(diagStat%stat_scalarSoilRelHumidity(1,1))
if (.not.allocated(diagStat%statIndex_scalarAquiferRootFrac)) allocate(diagStat%statIndex_scalarAquiferRootFrac(1))
if (.not.allocated(diagStat%stat_scalarAquiferRootFrac)) allocate(diagStat%stat_scalarAquiferRootFrac(1,1))
if (.not.allocated(diagStat%statIndex_scalarFracLiqVeg)) allocate(diagStat%statIndex_scalarFracLiqVeg(1))
if (.not.allocated(diagStat%stat_scalarFracLiqVeg)) allocate(diagStat%stat_scalarFracLiqVeg(1,1))
if (.not.allocated(diagStat%statIndex_scalarCanopyWetFraction)) allocate(diagStat%statIndex_scalarCanopyWetFraction(1))
if (.not.allocated(diagStat%stat_scalarCanopyWetFraction)) allocate(diagStat%stat_scalarCanopyWetFraction(1,1))
if (.not.allocated(diagStat%statIndex_scalarSnowAge)) allocate(diagStat%statIndex_scalarSnowAge(1))
if (.not.allocated(diagStat%stat_scalarSnowAge)) allocate(diagStat%stat_scalarSnowAge(1,1))
if (.not.allocated(diagStat%statIndex_scalarGroundSnowFraction)) allocate(diagStat%statIndex_scalarGroundSnowFraction(1))
if (.not.allocated(diagStat%stat_scalarGroundSnowFraction)) allocate(diagStat%stat_scalarGroundSnowFraction(1,1))
if (.not.allocated(diagStat%statIndex_scalarInfilArea)) allocate(diagStat%statIndex_scalarInfilArea(1))
if (.not.allocated(diagStat%stat_scalarInfilArea)) allocate(diagStat%stat_scalarInfilArea(1,1))
if (.not.allocated(diagStat%statIndex_scalarFrozenArea)) allocate(diagStat%statIndex_scalarFrozenArea(1))
if (.not.allocated(diagStat%stat_scalarFrozenArea)) allocate(diagStat%stat_scalarFrozenArea(1,1))
if (.not.allocated(diagStat%statIndex_scalarSoilControl)) allocate(diagStat%statIndex_scalarSoilControl(1))
if (.not.allocated(diagStat%stat_scalarSoilControl)) allocate(diagStat%stat_scalarSoilControl(1,1))
if (.not.allocated(diagStat%statIndex_scalarSoilCompress)) allocate(diagStat%statIndex_scalarSoilCompress(1))
if (.not.allocated(diagStat%stat_scalarSoilCompress)) allocate(diagStat%stat_scalarSoilCompress(1,1))
if (.not.allocated(diagStat%statIndex_scalarTotalSoilLiq)) allocate(diagStat%statIndex_scalarTotalSoilLiq(1))
if (.not.allocated(diagStat%stat_scalarTotalSoilLiq)) allocate(diagStat%stat_scalarTotalSoilLiq(1,1))
if (.not.allocated(diagStat%statIndex_scalarTotalSoilIce)) allocate(diagStat%statIndex_scalarTotalSoilIce(1))
if (.not.allocated(diagStat%stat_scalarTotalSoilIce)) allocate(diagStat%stat_scalarTotalSoilIce(1,1))
if (.not.allocated(diagStat%statIndex_scalarTotalSoilWat)) allocate(diagStat%statIndex_scalarTotalSoilWat(1))
if (.not.allocated(diagStat%stat_scalarTotalSoilWat)) allocate(diagStat%stat_scalarTotalSoilWat(1,1))
if (.not.allocated(diagStat%statIndex_scalarKappa)) allocate(diagStat%statIndex_scalarKappa(1))
if (.not.allocated(diagStat%stat_scalarKappa)) allocate(diagStat%stat_scalarKappa(1,1))
if (.not.allocated(diagStat%statIndex_scalarVolLatHt_fus)) allocate(diagStat%statIndex_scalarVolLatHt_fus(1))
if (.not.allocated(diagStat%stat_scalarVolLatHt_fus)) allocate(diagStat%stat_scalarVolLatHt_fus(1,1))
if (.not.allocated(diagStat%statIndex_numFluxCalls)) allocate(diagStat%statIndex_numFluxCalls(1))
if (.not.allocated(diagStat%stat_numFluxCalls)) allocate(diagStat%stat_numFluxCalls(1,1))
if (.not.allocated(diagStat%statIndex_wallClockTime)) allocate(diagStat%statIndex_wallClockTime(1))
if (.not.allocated(diagStat%stat_wallClockTime)) allocate(diagStat%stat_wallClockTime(1,1))
if (.not.allocated(diagStat%statIndex_meanStepSize)) allocate(diagStat%statIndex_meanStepSize(1))
if (.not.allocated(diagStat%stat_meanStepSize)) allocate(diagStat%stat_meanStepSize(1,1))
if (.not.allocated(diagStat%statIndex_balanceCasNrg)) allocate(diagStat%statIndex_balanceCasNrg(1))
if (.not.allocated(diagStat%stat_balanceCasNrg)) allocate(diagStat%stat_balanceCasNrg(1,1))
if (.not.allocated(diagStat%statIndex_balanceVegNrg)) allocate(diagStat%statIndex_balanceVegNrg(1))
if (.not.allocated(diagStat%stat_balanceVegNrg)) allocate(diagStat%stat_balanceVegNrg(1,1))
if (.not.allocated(diagStat%statIndex_balanceSnowNrg)) allocate(diagStat%statIndex_balanceSnowNrg(1))
if (.not.allocated(diagStat%stat_balanceSnowNrg)) allocate(diagStat%stat_balanceSnowNrg(1,1))
if (.not.allocated(diagStat%statIndex_balanceSoilNrg)) allocate(diagStat%statIndex_balanceSoilNrg(1))
if (.not.allocated(diagStat%stat_balanceSoilNrg)) allocate(diagStat%stat_balanceSoilNrg(1,1))
if (.not.allocated(diagStat%statIndex_balanceVegMass)) allocate(diagStat%statIndex_balanceVegMass(1))
if (.not.allocated(diagStat%stat_balanceVegMass)) allocate(diagStat%stat_balanceVegMass(1,1))
if (.not.allocated(diagStat%statIndex_balanceSnowMass)) allocate(diagStat%statIndex_balanceSnowMass(1))
if (.not.allocated(diagStat%stat_balanceSnowMass)) allocate(diagStat%stat_balanceSnowMass(1,1))
if (.not.allocated(diagStat%statIndex_balanceSoilMass)) allocate(diagStat%statIndex_balanceSoilMass(1))
if (.not.allocated(diagStat%stat_balanceSoilMass)) allocate(diagStat%stat_balanceSoilMass(1,1))
if (.not.allocated(diagStat%statIndex_balanceAqMass)) allocate(diagStat%statIndex_balanceAqMass(1))
if (.not.allocated(diagStat%stat_balanceAqMass)) allocate(diagStat%stat_balanceAqMass(1,1))
if (.not.allocated(diagStat%statIndex_numSteps)) allocate(diagStat%statIndex_numSteps(1))
if (.not.allocated(diagStat%stat_numSteps)) allocate(diagStat%stat_numSteps(1,1))
if (.not.allocated(diagStat%statIndex_numResEvals)) allocate(diagStat%statIndex_numResEvals(1))
if (.not.allocated(diagStat%stat_numResEvals)) allocate(diagStat%stat_numResEvals(1,1))
if (.not.allocated(diagStat%statIndex_numLinSolvSetups)) allocate(diagStat%statIndex_numLinSolvSetups(1))
if (.not.allocated(diagStat%stat_numLinSolvSetups)) allocate(diagStat%stat_numLinSolvSetups(1,1))
if (.not.allocated(diagStat%statIndex_numErrTestFails)) allocate(diagStat%statIndex_numErrTestFails(1))
if (.not.allocated(diagStat%stat_numErrTestFails)) allocate(diagStat%stat_numErrTestFails(1,1))
if (.not.allocated(diagStat%statIndex_kLast)) allocate(diagStat%statIndex_kLast(1))
if (.not.allocated(diagStat%stat_kLast)) allocate(diagStat%stat_kLast(1,1))
if (.not.allocated(diagStat%statIndex_kCur)) allocate(diagStat%statIndex_kCur(1))
if (.not.allocated(diagStat%stat_kCur)) allocate(diagStat%stat_kCur(1,1))
if (.not.allocated(diagStat%statIndex_hInitUsed)) allocate(diagStat%statIndex_hInitUsed(1))
if (.not.allocated(diagStat%stat_hInitUsed)) allocate(diagStat%stat_hInitUsed(1,1))
if (.not.allocated(diagStat%statIndex_hLast)) allocate(diagStat%statIndex_hLast(1))
if (.not.allocated(diagStat%stat_hLast)) allocate(diagStat%stat_hLast(1,1))
if (.not.allocated(diagStat%statIndex_hCur)) allocate(diagStat%statIndex_hCur(1))
if (.not.allocated(diagStat%stat_hCur)) allocate(diagStat%stat_hCur(1,1))
if (.not.allocated(diagStat%statIndex_tCur)) allocate(diagStat%statIndex_tCur(1))
if (.not.allocated(diagStat%stat_tCur)) allocate(diagStat%stat_tCur(1,1))
end subroutine

subroutine get_device_diag_stat(diagStat, iVar, nGRU, iFreq, parVector)
  implicit none
  type(diag_stat_device) :: diagStat
  integer(i4b) :: nGRU, iVar, iFreq
  real(rkind) :: parVector(:)
  select case(iVar)
case(iLookDIAG%scalarCanopyDepth); parVector = diagStat%stat_scalarCanopyDepth(iFreq,:)
case(iLookDIAG%scalarBulkVolHeatCapVeg); parVector = diagStat%stat_scalarBulkVolHeatCapVeg(iFreq,:)
case(iLookDIAG%scalarCanopyCm); parVector = diagStat%stat_scalarCanopyCm(iFreq,:)
case(iLookDIAG%scalarCanopyEmissivity); parVector = diagStat%stat_scalarCanopyEmissivity(iFreq,:)
case(iLookDIAG%scalarRootZoneTemp); parVector = diagStat%stat_scalarRootZoneTemp(iFreq,:)
case(iLookDIAG%scalarLAI); parVector = diagStat%stat_scalarLAI(iFreq,:)
case(iLookDIAG%scalarSAI); parVector = diagStat%stat_scalarSAI(iFreq,:)
case(iLookDIAG%scalarExposedLAI); parVector = diagStat%stat_scalarExposedLAI(iFreq,:)
case(iLookDIAG%scalarExposedSAI); parVector = diagStat%stat_scalarExposedSAI(iFreq,:)
case(iLookDIAG%scalarAdjMeasHeight); parVector = diagStat%stat_scalarAdjMeasHeight(iFreq,:)
case(iLookDIAG%scalarCanopyIceMax); parVector = diagStat%stat_scalarCanopyIceMax(iFreq,:)
case(iLookDIAG%scalarCanopyLiqMax); parVector = diagStat%stat_scalarCanopyLiqMax(iFreq,:)
case(iLookDIAG%scalarGrowingSeasonIndex); parVector = diagStat%stat_scalarGrowingSeasonIndex(iFreq,:)
case(iLookDIAG%scalarVolHtCap_air); parVector = diagStat%stat_scalarVolHtCap_air(iFreq,:)
case(iLookDIAG%scalarVolHtCap_ice); parVector = diagStat%stat_scalarVolHtCap_ice(iFreq,:)
case(iLookDIAG%scalarVolHtCap_soil); parVector = diagStat%stat_scalarVolHtCap_soil(iFreq,:)
case(iLookDIAG%scalarVolHtCap_water); parVector = diagStat%stat_scalarVolHtCap_water(iFreq,:)
case(iLookDIAG%scalarLambda_drysoil); parVector = diagStat%stat_scalarLambda_drysoil(iFreq,:)
case(iLookDIAG%scalarLambda_wetsoil); parVector = diagStat%stat_scalarLambda_wetsoil(iFreq,:)
case(iLookDIAG%scalarCanopyEnthTemp); parVector = diagStat%stat_scalarCanopyEnthTemp(iFreq,:)
case(iLookDIAG%scalarTotalSoilEnthalpy); parVector = diagStat%stat_scalarTotalSoilEnthalpy(iFreq,:)
case(iLookDIAG%scalarTotalSnowEnthalpy); parVector = diagStat%stat_scalarTotalSnowEnthalpy(iFreq,:)
case(iLookDIAG%scalarVPair); parVector = diagStat%stat_scalarVPair(iFreq,:)
case(iLookDIAG%scalarVP_CanopyAir); parVector = diagStat%stat_scalarVP_CanopyAir(iFreq,:)
case(iLookDIAG%scalarTwetbulb); parVector = diagStat%stat_scalarTwetbulb(iFreq,:)
case(iLookDIAG%scalarSnowfallTemp); parVector = diagStat%stat_scalarSnowfallTemp(iFreq,:)
case(iLookDIAG%scalarNewSnowDensity); parVector = diagStat%stat_scalarNewSnowDensity(iFreq,:)
case(iLookDIAG%scalarO2air); parVector = diagStat%stat_scalarO2air(iFreq,:)
case(iLookDIAG%scalarCO2air); parVector = diagStat%stat_scalarCO2air(iFreq,:)
case(iLookDIAG%windspd_x); parVector = diagStat%stat_windspd_x(iFreq,:)
case(iLookDIAG%windspd_y); parVector = diagStat%stat_windspd_y(iFreq,:)
case(iLookDIAG%scalarCosZenith); parVector = diagStat%stat_scalarCosZenith(iFreq,:)
case(iLookDIAG%scalarFractionDirect); parVector = diagStat%stat_scalarFractionDirect(iFreq,:)
case(iLookDIAG%scalarCanopySunlitFraction); parVector = diagStat%stat_scalarCanopySunlitFraction(iFreq,:)
case(iLookDIAG%scalarCanopySunlitLAI); parVector = diagStat%stat_scalarCanopySunlitLAI(iFreq,:)
case(iLookDIAG%scalarCanopyShadedLAI); parVector = diagStat%stat_scalarCanopyShadedLAI(iFreq,:)
case(iLookDIAG%scalarGroundAlbedo); parVector = diagStat%stat_scalarGroundAlbedo(iFreq,:)
case(iLookDIAG%scalarLatHeatSubVapCanopy); parVector = diagStat%stat_scalarLatHeatSubVapCanopy(iFreq,:)
case(iLookDIAG%scalarLatHeatSubVapGround); parVector = diagStat%stat_scalarLatHeatSubVapGround(iFreq,:)
case(iLookDIAG%scalarSatVP_CanopyTemp); parVector = diagStat%stat_scalarSatVP_CanopyTemp(iFreq,:)
case(iLookDIAG%scalarSatVP_GroundTemp); parVector = diagStat%stat_scalarSatVP_GroundTemp(iFreq,:)
case(iLookDIAG%scalarZ0Canopy); parVector = diagStat%stat_scalarZ0Canopy(iFreq,:)
case(iLookDIAG%scalarWindReductionFactor); parVector = diagStat%stat_scalarWindReductionFactor(iFreq,:)
case(iLookDIAG%scalarZeroPlaneDisplacement); parVector = diagStat%stat_scalarZeroPlaneDisplacement(iFreq,:)
case(iLookDIAG%scalarRiBulkCanopy); parVector = diagStat%stat_scalarRiBulkCanopy(iFreq,:)
case(iLookDIAG%scalarRiBulkGround); parVector = diagStat%stat_scalarRiBulkGround(iFreq,:)
case(iLookDIAG%scalarCanopyStabilityCorrection); parVector = diagStat%stat_scalarCanopyStabilityCorrection(iFreq,:)
case(iLookDIAG%scalarGroundStabilityCorrection); parVector = diagStat%stat_scalarGroundStabilityCorrection(iFreq,:)
case(iLookDIAG%scalarIntercellularCO2Sunlit); parVector = diagStat%stat_scalarIntercellularCO2Sunlit(iFreq,:)
case(iLookDIAG%scalarIntercellularCO2Shaded); parVector = diagStat%stat_scalarIntercellularCO2Shaded(iFreq,:)
case(iLookDIAG%scalarTranspireLim); parVector = diagStat%stat_scalarTranspireLim(iFreq,:)
case(iLookDIAG%scalarTranspireLimAqfr); parVector = diagStat%stat_scalarTranspireLimAqfr(iFreq,:)
case(iLookDIAG%scalarFoliageNitrogenFactor); parVector = diagStat%stat_scalarFoliageNitrogenFactor(iFreq,:)
case(iLookDIAG%scalarSoilRelHumidity); parVector = diagStat%stat_scalarSoilRelHumidity(iFreq,:)
case(iLookDIAG%scalarAquiferRootFrac); parVector = diagStat%stat_scalarAquiferRootFrac(iFreq,:)
case(iLookDIAG%scalarFracLiqVeg); parVector = diagStat%stat_scalarFracLiqVeg(iFreq,:)
case(iLookDIAG%scalarCanopyWetFraction); parVector = diagStat%stat_scalarCanopyWetFraction(iFreq,:)
case(iLookDIAG%scalarSnowAge); parVector = diagStat%stat_scalarSnowAge(iFreq,:)
case(iLookDIAG%scalarGroundSnowFraction); parVector = diagStat%stat_scalarGroundSnowFraction(iFreq,:)
case(iLookDIAG%scalarInfilArea); parVector = diagStat%stat_scalarInfilArea(iFreq,:)
case(iLookDIAG%scalarFrozenArea); parVector = diagStat%stat_scalarFrozenArea(iFreq,:)
case(iLookDIAG%scalarSoilControl); parVector = diagStat%stat_scalarSoilControl(iFreq,:)
case(iLookDIAG%scalarSoilCompress); parVector = diagStat%stat_scalarSoilCompress(iFreq,:)
case(iLookDIAG%scalarTotalSoilLiq); parVector = diagStat%stat_scalarTotalSoilLiq(iFreq,:)
case(iLookDIAG%scalarTotalSoilIce); parVector = diagStat%stat_scalarTotalSoilIce(iFreq,:)
case(iLookDIAG%scalarTotalSoilWat); parVector = diagStat%stat_scalarTotalSoilWat(iFreq,:)
case(iLookDIAG%scalarKappa); parVector = diagStat%stat_scalarKappa(iFreq,:)
case(iLookDIAG%scalarVolLatHt_fus); parVector = diagStat%stat_scalarVolLatHt_fus(iFreq,:)
case(iLookDIAG%numFluxCalls); parVector = diagStat%stat_numFluxCalls(iFreq,:)
case(iLookDIAG%wallClockTime); parVector = diagStat%stat_wallClockTime(iFreq,:)
case(iLookDIAG%meanStepSize); parVector = diagStat%stat_meanStepSize(iFreq,:)
case(iLookDIAG%balanceCasNrg); parVector = diagStat%stat_balanceCasNrg(iFreq,:)
case(iLookDIAG%balanceVegNrg); parVector = diagStat%stat_balanceVegNrg(iFreq,:)
case(iLookDIAG%balanceSnowNrg); parVector = diagStat%stat_balanceSnowNrg(iFreq,:)
case(iLookDIAG%balanceSoilNrg); parVector = diagStat%stat_balanceSoilNrg(iFreq,:)
case(iLookDIAG%balanceVegMass); parVector = diagStat%stat_balanceVegMass(iFreq,:)
case(iLookDIAG%balanceSnowMass); parVector = diagStat%stat_balanceSnowMass(iFreq,:)
case(iLookDIAG%balanceSoilMass); parVector = diagStat%stat_balanceSoilMass(iFreq,:)
case(iLookDIAG%balanceAqMass); parVector = diagStat%stat_balanceAqMass(iFreq,:)
case(iLookDIAG%numSteps); parVector = diagStat%stat_numSteps(iFreq,:)
case(iLookDIAG%numResEvals); parVector = diagStat%stat_numResEvals(iFreq,:)
case(iLookDIAG%numLinSolvSetups); parVector = diagStat%stat_numLinSolvSetups(iFreq,:)
case(iLookDIAG%numErrTestFails); parVector = diagStat%stat_numErrTestFails(iFreq,:)
case(iLookDIAG%kLast); parVector = diagStat%stat_kLast(iFreq,:)
case(iLookDIAG%kCur); parVector = diagStat%stat_kCur(iFreq,:)
case(iLookDIAG%hInitUsed); parVector = diagStat%stat_hInitUsed(iFreq,:)
case(iLookDIAG%hLast); parVector = diagStat%stat_hLast(iFreq,:)
case(iLookDIAG%hCur); parVector = diagStat%stat_hCur(iFreq,:)
case(iLookDIAG%tCur); parVector = diagStat%stat_tCur(iFreq,:)
case default; parVector = -9999_rkind
end select
end subroutine

end module initialize_device      