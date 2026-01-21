module device_data_types
use cudafor 
USE nrtype

  integer(i4b),parameter        :: nLook=500            ! number of elements in the lookup table

type,public :: mpar_data_device

  real(rkind),device,allocatable :: soil_dens_intr_(:,:)
  real(rkind),device,allocatable :: thCond_soil_(:,:)
  real(rkind),device,allocatable :: frac_sand_(:,:)
  real(rkind),device,allocatable :: frac_clay_(:,:)
  real(rkind),device,allocatable :: frac_silt_(:,:)
  real(rkind),device,allocatable :: snowfrz_scale_(:)
  real(rkind),device,allocatable :: specificHeatVeg_(:)
  real(rkind),device,allocatable :: maxMassVegetation_(:)
  real(rkind),device,allocatable :: Fcapil_(:), k_snow_(:), mw_exp_(:)
  real(rkind),device,allocatable :: fixedThermalCond_snow_(:)

  real(rkind),device,allocatable :: heightCanopyTop_(:)
  real(rkind),device,allocatable :: heightCanopyBottom_(:)
  real(rkind),device,allocatable :: canopyWettingFactor_(:)
  real(rkind),device,allocatable :: canopyWettingExp_(:)

  real(rkind),device,allocatable :: z0Snow_(:)
  real(rkind),device,allocatable :: z0Soil_(:)
  real(rkind),device,allocatable :: z0Canopy_(:)
  real(rkind),device,allocatable :: zpdFraction_(:)
  real(rkind),device,allocatable :: critRichNumber_(:)
  real(rkind),device,allocatable :: Louis79_bparam_(:)
  real(rkind),device,allocatable :: Louis79_cStar_(:)
  real(rkind),device,allocatable :: Mahrt87_eScale_(:)
  real(rkind),device,allocatable :: windReductionParam_(:)
  real(rkind),device,allocatable :: leafExchangeCoeff_(:)
  real(rkind),device,allocatable :: leafDimension_(:)
  real(rkind),device,allocatable :: vGn_n_(:,:), vGn_alpha_(:,:)
  real(rkind),device,allocatable :: theta_sat_(:,:)
  real(rkind),device,allocatable :: theta_res_(:,:)
  real(rkind),device,allocatable :: plantWiltPsi_(:)
  real(rkind),device,allocatable :: soilStressParam_(:)
  real(rkind),device,allocatable :: critSoilWilting_(:)
  real(rkind),device,allocatable :: critSoilTranspire_(:)
  real(rkind),device,allocatable :: critAquiferTranspire_(:)
  real(rkind),device,allocatable :: minStomatalResistance_(:)

  real(rkind),device,allocatable :: upperBoundTemp_(:),lowerBoundTemp_(:)

  real(rkind),device,allocatable :: throughfallScaleRain_(:), canopyDrainageCoeff_(:)

  real(rkind),device,allocatable :: aquiferBaseflowRate_(:), aquiferScaleFactor_(:), aquiferBaseflowExp_(:)
  real(rkind),device,allocatable :: Kc25_(:)
  real(rkind),device,allocatable :: Ko25_(:)
  real(rkind),device,allocatable :: Kc_qFac_(:)
  real(rkind),device,allocatable :: Ko_qFac_(:)
  real(rkind),device,allocatable :: kc_Ha_(:)
  real(rkind),device,allocatable :: ko_Ha_(:)
  real(rkind),device,allocatable :: vcmax25_canopyTop_(:)
  real(rkind),device,allocatable :: vcmax_qFac_(:)
  real(rkind),device,allocatable :: vcmax_Ha_(:)
  real(rkind),device,allocatable :: vcmax_Hd_(:)
  real(rkind),device,allocatable :: vcmax_Sv_(:)
  real(rkind),device,allocatable :: vcmax_Kn_(:)
  real(rkind),device,allocatable :: jmax25_scale_(:)
  real(rkind),device,allocatable :: jmax_Ha_(:)
  real(rkind),device,allocatable :: jmax_Hd_(:)
  real(rkind),device,allocatable :: jmax_Sv_(:)
  real(rkind),device,allocatable :: fractionJ_(:)
  real(rkind),device,allocatable :: quantamYield_(:)
  real(rkind),device,allocatable :: vpScaleFactor_(:)
  real(rkind),device,allocatable :: cond2photo_slope_(:)
  real(rkind),device,allocatable :: minStomatalConductance_(:)
  real(rkind),device,allocatable :: upperBoundHead_(:), upperBoundTheta_(:), lowerBoundHead_(:), lowerBoundTheta_(:)
  real(rkind),device,allocatable :: wettingFrontSuction_(:),rootingDepth_(:),kAnisotropic_(:),zScale_TOPMODEL_(:),qSurfScale_(:),f_impede_(:),soilIceScale_(:),soilIceCV_(:)
  real(rkind),device,allocatable :: theta_mp_(:), mpExp_(:)
  real(rkind),device,allocatable :: specificStorage_(:)
  real(rkind),device,allocatable :: fieldCapacity_(:)
  real(rkind) :: idaMaxOrder
  real(rkind) :: idaMaxInternalSteps
  real(rkind) :: idaInitStepSize
  real(rkind) :: idaMinStepSize
  real(rkind) :: idaMaxErrTestFail

  real(rkind),device,allocatable :: absTolTempCas
  real(rkind),device,allocatable :: relTolTempCas
  real(rkind),device,allocatable :: absTolTempVeg
  real(rkind),device,allocatable :: relTolTempVeg
  real(rkind),device,allocatable :: absTolWatVeg
  real(rkind),device,allocatable :: relTolWatVeg
  real(rkind),device,allocatable :: absTolTempSoilSnow
  real(rkind),device,allocatable :: relTolTempSoilSnow
  real(rkind),device,allocatable :: absTolWatSnow
  real(rkind),device,allocatable :: relTolWatSnow
  real(rkind),device,allocatable :: absTolMatric
  real(rkind),device,allocatable :: relTolMatric
  real(rkind),device,allocatable :: absTolAquifr
  real(rkind),device,allocatable :: relTolAquifr
  real(rkind) :: maxstep
  real(rkind) :: be_steps

  real(rkind),device,allocatable :: tempCritRain_(:)
  real(rkind),device,allocatable :: tempRangeTimestep_(:)
  real(rkind),device,allocatable :: frozenPrecipMultip_(:)
  real(rkind),device,allocatable :: albedoMinWinter_(:)
  real(rkind),device,allocatable :: albedoMinSpring_(:)
  real(rkind),device,allocatable :: albedoMinVisible_(:)
  real(rkind),device,allocatable :: albedoMinNearIR_(:)
  real(rkind),device,allocatable :: albedoDecayRate_(:)
  real(rkind),device,allocatable :: albedoSootLoad_(:)
  real(rkind),device,allocatable :: albedoRefresh_(:)
  ! real(rkind) :: radExt_snow
  real(rkind),device,allocatable :: directScale_(:)
  real(rkind),device,allocatable :: Frad_direct_(:)
  real(rkind),device,allocatable :: newSnowDenMin_(:)
  real(rkind),device,allocatable :: newSnowDenMult_(:)
  real(rkind),device,allocatable :: newSnowDenScal_(:)
  real(rkind),device,allocatable :: constSnowDen_(:)
  real(rkind),device,allocatable :: newSnowDenAdd_(:)
  real(rkind),device,allocatable :: newSnowDenMultTemp_(:)
  real(rkind),device,allocatable :: newSnowDenMultWind_(:)
  real(rkind),device,allocatable :: newSnowDenMultAnd_(:)
  real(rkind),device,allocatable :: newSnowDenBase_(:)
  ! real(rkind) :: densScalGrowth
  ! real(rkind) :: tempScalGrowth
  ! real(rkind) :: grainGrowthRate
  ! real(rkind) :: densScalOvrbdn
  ! real(rkind) :: tempScalOvrbdn
  ! real(rkind) :: baseViscosity
  real(rkind),device,allocatable :: winterSAI_(:)
  real(rkind),device,allocatable :: summerLAI_(:)
  real(rkind),device,allocatable :: rootScaleFactor1_(:)
  real(rkind),device,allocatable :: rootScaleFactor2_(:)
  real(rkind),device,allocatable :: rootDistExp_(:)
  ! real(rkind) :: throughfallScaleSnow
  real(rkind),device,allocatable :: refInterceptCapSnow_(:)
  real(rkind),device,allocatable :: refInterceptCapRain_(:)
  real(rkind),device,allocatable :: snowUnloadingCoeff_(:)
  real(rkind),device,allocatable :: ratioDrip2Unloading_(:)
  real(rkind),device,allocatable :: minTempUnloading_(:)
  real(rkind),device,allocatable :: rateTempUnloading_(:)
  real(rkind),device,allocatable :: minWindUnloading_(:)
  real(rkind),device,allocatable :: rateWindUnloading_(:)

  real(rkind),device,allocatable :: compactedDepth_(:)
  ! real(rkind) :: specificYield
  real(rkind),device,allocatable :: minwind_(:)
  real(rkind) :: minstep
  ! real(rkind) :: maxstep
  ! real(rkind) :: be_steps
  ! real(rkind) :: wimplicit
  ! real(rkind) :: maxiter
  ! real(rkind) :: relConvTol_liquid
  ! real(rkind) :: absConvTol_liquid
  ! real(rkind) :: relConvTol_matric
  ! real(rkind) :: absConvTol_matric
  ! real(rkind) :: relConvTol_energy
  ! real(rkind) :: absConvTol_energy
  ! real(rkind) :: relConvTol_aquifr
  ! real(rkind) :: absConvTol_aquifr
  ! real(rkind) :: relTolTempCas
  ! real(rkind) :: absTolTempCas
  ! real(rkind) :: relTolTempVeg
  ! real(rkind) :: absTolTempVeg
  ! real(rkind) :: relTolWatVeg
  ! real(rkind) :: absTolWatVeg
  ! real(rkind) :: relTolTempSoilSnow
  ! real(rkind) :: absTolTempSoilSnow
  ! real(rkind) :: relTolWatSnow
  ! real(rkind) :: absTolWatSnow
  ! real(rkind) :: relTolMatric
  ! real(rkind) :: absTolMatric
  ! real(rkind) :: relTolAquifr
  ! real(rkind) :: absTolAquifr
  real(rkind),device,allocatable :: zmin_(:)
  real(rkind),device,allocatable :: zminLayer1_(:)
  real(rkind),device,allocatable :: zminLayer2_(:)
  real(rkind),device,allocatable :: zminLayer3_(:)
  real(rkind),device,allocatable :: zminLayer4_(:)
  real(rkind),device,allocatable :: zminLayer5_(:)

  real(rkind),device,allocatable :: k_soil_(:,:)
  real(rkind),device,allocatable :: k_macropore_(:,:)
  real(rkind),allocatable,device :: densScalGrowth_(:)
  real(rkind),allocatable,device :: tempScalGrowth_(:)
  real(rkind),allocatable,device :: densScalOvrbdn_(:)
  real(rkind),allocatable,device :: tempScalOvrbdn_(:)
  real(rkind),allocatable,device :: baseViscosity_(:)
  real(rkind),allocatable,device :: grainGrowthRate_(:)
  real(rkind),device,allocatable :: zmaxLayer1_lower_(:)
  real(rkind),device,allocatable :: zmaxLayer2_lower_(:)
  real(rkind),device,allocatable :: zmaxLayer3_lower_(:)
  real(rkind),device,allocatable :: zmaxLayer4_lower_(:)
  real(rkind),device,allocatable :: zmaxLayer1_upper_(:)
  real(rkind),device,allocatable :: zmaxLayer2_upper_(:)
  real(rkind),device,allocatable :: zmaxLayer3_upper_(:)
  real(rkind),device,allocatable :: zmaxLayer4_upper_(:)
  real(rkind),device,allocatable :: zmax_(:)
  real(rkind),device,allocatable :: Frad_vis_(:)
  real(rkind),device,allocatable :: albedoMaxNearIR_(:)
  real(rkind),device,allocatable :: albedoMaxVisible_(:)
  real(rkind),device,allocatable :: albedoMax_(:)
  real(rkind),device,allocatable :: FUSE_Ac_max_(:)
  real(rkind),device,allocatable :: FUSE_phi_tens_(:)
  real(rkind),device,allocatable :: FUSE_b_(:)
  real(rkind),device,allocatable :: FUSE_lambda_(:)
  real(rkind),device,allocatable :: FUSE_chi_(:)
  real(rkind),device,allocatable :: FUSE_mu_(:)
  real(rkind),device,allocatable :: FUSE_n_(:)

end type mpar_data_device

type, public :: diag_data_device !!!!

  ! local properties
  real(rkind),device,allocatable :: scalarCanopyDepth(:)               ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarBulkVolHeatCapVeg(:)         ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyCm(:)                  ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyEmissivity(:)          ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarRootZoneTemp(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLAI(:)                       ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSAI(:)                       ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarExposedLAI(:)                ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarExposedSAI(:)                ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarAdjMeasHeight(:)             ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyIceMax(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyLiqMax(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarGrowingSeasonIndex(:)        ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarVolHtCap_air(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarVolHtCap_ice(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarVolHtCap_soil(:)             ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarVolHtCap_water(:)            ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: mLayerVolHtCapBulk(:,:)              ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: mLayerCm(:,:)                        ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: scalarLambda_drysoil(:)            ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLambda_wetsoil(:)            ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: mLayerThermalC(:,:)                  ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: iLayerThermalC(:,:)                  ! get_ixVarType('ifcToto')
  ! enthalpy
  real(rkind),device,allocatable :: scalarCanopyEnthTemp(:)            ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: mLayerEnthTemp(:,:)                  ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: scalarTotalSoilEnthalpy(:)         ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarTotalSnowEnthalpy(:)         ! get_ixVarType('scalarv')
  ! forcing
  real(rkind),device,allocatable :: scalarVPair(:)                     ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarVP_CanopyAir(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarTwetbulb(:)                  ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSnowfallTemp(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarNewSnowDensity(:)            ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarO2air(:)                     ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCO2air(:)                    ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: windspd_x(:)                       ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: windspd_y(:)                       ! get_ixVarType('scalarv')
  ! shortwave radiation
  real(rkind),device,allocatable :: scalarCosZenith(:)                 ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarFractionDirect(:)            ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopySunlitFraction(:)      ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopySunlitLAI(:)           ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyShadedLAI(:)           ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: spectralAlbGndDirect(:,:)            ! get_ixVarType('wLength')
  real(rkind),device,allocatable :: spectralAlbGndDiffuse(:,:)           ! get_ixVarType('wLength')
  real(rkind),device,allocatable :: scalarGroundAlbedo(:)              ! get_ixVarType('scalarv')
  ! turbulent heat transfer
  real(rkind),device,allocatable :: scalarLatHeatSubVapCanopy(:)       ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLatHeatSubVapGround(:)       ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSatVP_CanopyTemp(:)          ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSatVP_GroundTemp(:)          ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarZ0Canopy(:)                  ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarWindReductionFactor(:)       ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarZeroPlaneDisplacement(:)     ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarRiBulkCanopy(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarRiBulkGround(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyStabilityCorrection(:) ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarGroundStabilityCorrection(:) ! get_ixVarType('scalarv')
  ! evapotranspiration
  real(rkind),device,allocatable :: scalarIntercellularCO2Sunlit(:)    ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarIntercellularCO2Shaded(:)    ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarTranspireLim(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarTranspireLimAqfr(:)          ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarFoliageNitrogenFactor(:)     ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSoilRelHumidity(:)           ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: mLayerTranspireLim(:,:)              ! get_ixVarType('midSoil')
  real(rkind),device,allocatable :: mLayerRootDensity(:,:)               ! get_ixVarType('midSoil')
  real(rkind),device,allocatable :: scalarAquiferRootFrac(:)           ! get_ixVarType('scalarv')
  ! canopy hydrology
  real(rkind),device,allocatable :: scalarFracLiqVeg(:)                ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyWetFraction(:)         ! get_ixVarType('scalarv')
  ! snow hydrology
  real(rkind),device,allocatable :: scalarSnowAge(:)                   ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarGroundSnowFraction(:)        ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: spectralSnowAlbedoDirect(:,:)        ! get_ixVarType('wLength')
  real(rkind),device,allocatable :: mLayerFracLiqSnow(:,:)               ! get_ixVarType('midSnow')
  real(rkind),device,allocatable :: mLayerThetaResid(:,:)                ! get_ixVarType('midSnow')
  real(rkind),device,allocatable :: mLayerPoreSpace(:,:)                 ! get_ixVarType('midSnow')
  real(rkind),device,allocatable :: mLayerMeltFreeze(:,:)                ! get_ixVarType('midToto')
  ! soil hydrology
  real(rkind),device,allocatable :: scalarInfilArea(:)                 ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarFrozenArea(:)                ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSoilControl(:)               ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: mLayerVolFracAir(:,:)                ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: mLayerTcrit(:,:)                     ! get_ixVarType('midSoil')
  real(rkind),device,allocatable :: mLayerCompress(:,:)                  ! get_ixVarType('midSoil')
  real(rkind),device,allocatable :: scalarSoilCompress(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: mLayerMatricHeadLiq(:,:)             ! get_ixVarType('midSoil')
  real(rkind),device,allocatable :: scalarTotalSoilLiq(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarTotalSoilIce(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarTotalSoilWat(:)              ! get_ixVarType('scalarv')
  ! variable shortcuts
  real(rkind),device,allocatable :: scalarVGn_m(:,:)                     ! get_ixVarType('midSoil')
  real(rkind),device,allocatable :: scalarKappa(:)                     ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarVolLatHt_fus(:)              ! get_ixVarType('scalarv')
  ! timing information
  real(rkind) :: numFluxCalls                    ! get_ixVarType('scalarv')
  real(rkind) :: wallClockTime                   ! get_ixVarType('scalarv')
  real(rkind) :: meanStepSize                    ! get_ixVarType('scalarv')
  ! balances
  real(rkind),device,allocatable :: balanceCasNrg(:)                   ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: balanceVegNrg(:)                   ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: balanceLayerNrg(:,:)                 ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: balanceSnowNrg(:)                  ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: balanceSoilNrg(:)                  ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: balanceVegMass(:)                  ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: balanceLayerMass(:,:)                ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: balanceSnowMass(:)                 ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: balanceSoilMass(:)                 ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: balanceAqMass(:)                   ! get_ixVarType('scalarv')
  ! sundials integrator stats
  real(rkind) :: numSteps
  real(rkind) :: numResEvals
  real(rkind) :: numLinSolvSetups
  real(rkind) :: numErrTestFails
  real(rkind) :: kLast
  real(rkind) :: kCur
  real(rkind) :: hInitUsed
  real(rkind) :: hLast
  real(rkind) :: hCur
  real(rkind) :: tCur
  
  
end type diag_data_device



type,public:: flux_data_device !!!!
  real(rkind),device,allocatable :: data(:,:)
  integer(i4b) :: ixScalarCanairNetNrgFlux = 1          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyNetNrgFlux = 2          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundNetNrgFlux = 3          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyNetLiqFlux = 4          !get_ixVarType('scalarv')
  ! forcing
  integer(i4b) :: ixScalarRainfall = 5                  ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSnowfall = 6                  ! get_ixVarType('scalarv')
  ! shortwave radiation
  integer(i4b) :: ixSpectralIncomingDirect_start, ixSpectralIncomingDirect_end          ! get_ixVarType('wLength')
  integer(i4b) :: ixSpectralIncomingDiffuse_start, ixSpectralIncomingDiffuse_end         ! get_ixVarType('wLength')
  integer(i4b) :: ixScalarCanopySunlitPAR = 7           ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyShadedPAR = 8           ! get_ixVarType('scalarv')
  integer(i4b) :: ixSpectralBelowCanopyDirect_start,ixSpectralBelowCanopyDirect_end       ! get_ixVarType('wLength')
  integer(i4b) :: ixSpectralBelowCanopyDiffuse_start,ixSpectralBelowCanopyDiffuse_end      ! get_ixVarType('wLength')
  integer(i4b) :: ixScalarBelowCanopySolar = 9          ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyAbsorbedSolar = 10       ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundAbsorbedSolar = 11       ! get_ixVarType('scalarv')
  ! longwave radiation
  integer(i4b) :: ixScalarLWRadCanopy = 12               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadGround = 13               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadUbound2Canopy = 14        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadUbound2Ground = 15        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadUbound2Ubound = 16        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadCanopy2Ubound = 17        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadCanopy2Ground = 18        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadCanopy2Canopy = 19        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadGround2Ubound = 20        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadGround2Canopy = 21        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWNetCanopy = 22               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWNetGround = 23               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWNetUbound = 24               !get_ixVarType('scalarv')
  ! turbulent heat transfer
  integer(i4b) :: ixScalarEddyDiffusCanopyTop = 25       !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarFrictionVelocity = 26          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarWindspdCanopyTop = 27          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarWindspdCanopyBottom = 28       !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundResistance = 29          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyResistance = 30          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLeafResistance = 31            !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSoilResistance = 32            !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSenHeatTotal = 33              !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSenHeatCanopy = 34             !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSenHeatGround = 35             !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatTotal = 36              !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatCanopyEvap = 37         !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatCanopyTrans = 38        !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatGround = 39             !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyAdvectiveHeatFlux = 40   !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundAdvectiveHeatFlux = 41   !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopySublimation = 42         !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSnowSublimation = 43           !, get_ixVarType('scalarv')
  ! liquid water fluxes associated with evapotranspiration
  integer(i4b) :: ixScalarStomResistSunlit = 44          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarStomResistShaded = 45          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarPhotosynthesisSunlit = 46      !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarPhotosynthesisShaded = 47      !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyTranspiration = 48       !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyEvaporation = 49         !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundEvaporation = 50         !get_ixVarType('scalarv')
  integer(i4b) :: ixMLayerTranspire_start, ixMLayerTranspire_end                 !get_ixVarType('midSoil')
  ! liquid and solid water fluxes through the canopy
  integer(i4b) :: ixScalarThroughfallSnow = 51 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarThroughfallRain = 52 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopySnowUnloading = 53 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyLiqDrainage = 54 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyMeltFreeze = 55 ! get_ixVarType('scalarv')
  ! energy fluxes and for the snow and soil domains
  integer(i4b) :: ixILayerConductiveFlux_start,ixILayerConductiveFlux_end ! get_ixVarType('ifcToto')
  integer(i4b) :: ixILayerAdvectiveFlux_start, ixILayerAdvectiveFlux_end ! get_ixVarType('ifcToto')
  integer(i4b) :: ixILayerNrgFlux_start, ixILayerNrgFlux_end ! get_ixVarType('ifcToto')
  integer(i4b) :: ixMLayerNrgFlux_start, ixMLayerNrgFlux_end ! get_ixVarType('midToto')
  ! liquid water fluxes for the snow domain
  integer(i4b) :: ixScalarSnowDrainage = 56! get_ixVarType('scalarv')
  integer(i4b) :: ixILayerLiqFluxSnow_start, ixILayerLiqFluxSnow_end! get_ixVarType('ifcSnow')
  integer(i4b) :: ixMLayerLiqFluxSnow_start, ixMLayerLiqFluxSnow_end! get_ixVarType('midSnow')
  ! liquid water fluxes for the soil domain
  integer(i4b) ::ixScalarRainPlusMelt = 57!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarMaxInfilRate = 58!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarInfiltration = 59!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarExfiltration = 60!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarSurfaceRunoff = 61!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarSurfaceRunoff_IE = 62
  integer(i4b) :: ixScalarSurfaceRunoff_SE = 63
  integer(i4b) :: ixMLayerSatHydCondMP_start, ixMLayerSatHydCondMP_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerSatHydCond_start, ixMLayerSatHydCond_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixILayerSatHydCond_start, ixILayerSatHydCond_end!, get_ixVarType('ifcSoil'), 
  integer(i4b) :: ixMLayerHydCond_start, ixMLayerHydCond_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixILayerLiqFluxSoil_start, ixILayerLiqFluxSoil_end!, get_ixVarType('ifcSoil'), 
  integer(i4b) :: ixMLayerLiqFluxSoil_start, ixMLayerLiqFluxSoil_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerBaseflow_start, ixMLayerBaseflow_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerColumnInflow_start, ixMLayerColumnInflow_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerColumnOutflow_start, ixMLayerColumnOutflow_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixScalarSoilBaseflow = 64!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarSoilDrainage = 65!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarAquiferRecharge = 66!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarAquiferTranspire = 67!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarAquiferBaseflow = 68!, get_ixVarType('scalarv'), 
  ! derived variables
  integer(i4b) :: ixScalarTotalET = 69
  integer(i4b) :: ixScalarTotalRunoff = 70
  integer(i4b) :: ixScalarNetRadiation = 71
  integer(i4b) :: numScalarFluxData = 71
  integer(i4b) :: numFluxData


end type flux_data_device
type,public:: flux_mask_device !!!!
  logical(lgt),device,allocatable :: data(:,:)
  integer(i4b) :: ixScalarCanairNetNrgFlux = 1          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyNetNrgFlux = 2          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundNetNrgFlux = 3          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyNetLiqFlux = 4          !get_ixVarType('scalarv')
  ! forcing
  integer(i4b) :: ixScalarRainfall = 5                  ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSnowfall = 6                  ! get_ixVarType('scalarv')
  ! shortwave radiation
  integer(i4b) :: ixSpectralIncomingDirect_start, ixSpectralIncomingDirect_end          ! get_ixVarType('wLength')
  integer(i4b) :: ixSpectralIncomingDiffuse_start, ixSpectralIncomingDiffuse_end         ! get_ixVarType('wLength')
  integer(i4b) :: ixScalarCanopySunlitPAR = 7           ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyShadedPAR = 8           ! get_ixVarType('scalarv')
  integer(i4b) :: ixSpectralBelowCanopyDirect_start,ixSpectralBelowCanopyDirect_end       ! get_ixVarType('wLength')
  integer(i4b) :: ixSpectralBelowCanopyDiffuse_start,ixSpectralBelowCanopyDiffuse_end      ! get_ixVarType('wLength')
  integer(i4b) :: ixScalarBelowCanopySolar = 9          ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyAbsorbedSolar = 10       ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundAbsorbedSolar = 11       ! get_ixVarType('scalarv')
  ! longwave radiation
  integer(i4b) :: ixScalarLWRadCanopy = 12               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadGround = 13               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadUbound2Canopy = 14        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadUbound2Ground = 15        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadUbound2Ubound = 16        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadCanopy2Ubound = 17        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadCanopy2Ground = 18        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadCanopy2Canopy = 19        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadGround2Ubound = 20        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadGround2Canopy = 21        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWNetCanopy = 22               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWNetGround = 23               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWNetUbound = 24               !get_ixVarType('scalarv')
  ! turbulent heat transfer
  integer(i4b) :: ixScalarEddyDiffusCanopyTop = 25       !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarFrictionVelocity = 26          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarWindspdCanopyTop = 27          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarWindspdCanopyBottom = 28       !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundResistance = 29          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyResistance = 30          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLeafResistance = 31            !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSoilResistance = 32            !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSenHeatTotal = 33              !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSenHeatCanopy = 34             !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSenHeatGround = 35             !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatTotal = 36              !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatCanopyEvap = 37         !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatCanopyTrans = 38        !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatGround = 39             !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyAdvectiveHeatFlux = 40   !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundAdvectiveHeatFlux = 41   !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopySublimation = 42         !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSnowSublimation = 43           !, get_ixVarType('scalarv')
  ! liquid water fluxes associated with evapotranspiration
  integer(i4b) :: ixScalarStomResistSunlit = 44          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarStomResistShaded = 45          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarPhotosynthesisSunlit = 46      !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarPhotosynthesisShaded = 47      !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyTranspiration = 48       !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyEvaporation = 49         !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundEvaporation = 50         !get_ixVarType('scalarv')
  integer(i4b) :: ixMLayerTranspire_start, ixMLayerTranspire_end                 !get_ixVarType('midSoil')
  ! liquid and solid water fluxes through the canopy
  integer(i4b) :: ixScalarThroughfallSnow = 51 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarThroughfallRain = 52 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopySnowUnloading = 53 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyLiqDrainage = 54 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyMeltFreeze = 55 ! get_ixVarType('scalarv')
  ! energy fluxes and for the snow and soil domains
  integer(i4b) :: ixILayerConductiveFlux_start,ixILayerConductiveFlux_end ! get_ixVarType('ifcToto')
  integer(i4b) :: ixILayerAdvectiveFlux_start, ixILayerAdvectiveFlux_end ! get_ixVarType('ifcToto')
  integer(i4b) :: ixILayerNrgFlux_start, ixILayerNrgFlux_end ! get_ixVarType('ifcToto')
  integer(i4b) :: ixMLayerNrgFlux_start, ixMLayerNrgFlux_end ! get_ixVarType('midToto')
  ! liquid water fluxes for the snow domain
  integer(i4b) :: ixScalarSnowDrainage = 56! get_ixVarType('scalarv')
  integer(i4b) :: ixILayerLiqFluxSnow_start, ixILayerLiqFluxSnow_end! get_ixVarType('ifcSnow')
  integer(i4b) :: ixMLayerLiqFluxSnow_start, ixMLayerLiqFluxSnow_end! get_ixVarType('midSnow')
  ! liquid water fluxes for the soil domain
  integer(i4b) ::ixScalarRainPlusMelt = 57!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarMaxInfilRate = 58!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarInfiltration = 59!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarExfiltration = 60!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarSurfaceRunoff = 61!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarSurfaceRunoff_IE = 62
  integer(i4b) :: ixScalarSurfaceRunoff_SE = 63
  integer(i4b) :: ixMLayerSatHydCondMP_start, ixMLayerSatHydCondMP_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerSatHydCond_start, ixMLayerSatHydCond_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixILayerSatHydCond_start, ixILayerSatHydCond_end!, get_ixVarType('ifcSoil'), 
  integer(i4b) :: ixMLayerHydCond_start, ixMLayerHydCond_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixILayerLiqFluxSoil_start, ixILayerLiqFluxSoil_end!, get_ixVarType('ifcSoil'), 
  integer(i4b) :: ixMLayerLiqFluxSoil_start, ixMLayerLiqFluxSoil_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerBaseflow_start, ixMLayerBaseflow_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerColumnInflow_start, ixMLayerColumnInflow_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerColumnOutflow_start, ixMLayerColumnOutflow_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixScalarSoilBaseflow = 64!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarSoilDrainage = 65!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarAquiferRecharge = 66!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarAquiferTranspire = 67!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarAquiferBaseflow = 68!, get_ixVarType('scalarv'), 
  ! derived variables
  integer(i4b) :: ixScalarTotalET = 69
  integer(i4b) :: ixScalarTotalRunoff = 70
  integer(i4b) :: ixScalarNetRadiation = 71
  integer(i4b) :: numScalarFluxData = 71
  integer(i4b) :: numFluxData


end type flux_mask_device
type,public:: flux_count_device !!!!
  integer(i4b),device,allocatable :: data(:,:)
  integer(i4b) :: ixScalarCanairNetNrgFlux = 1          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyNetNrgFlux = 2          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundNetNrgFlux = 3          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyNetLiqFlux = 4          !get_ixVarType('scalarv')
  ! forcing
  integer(i4b) :: ixScalarRainfall = 5                  ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSnowfall = 6                  ! get_ixVarType('scalarv')
  ! shortwave radiation
  integer(i4b) :: ixSpectralIncomingDirect_start, ixSpectralIncomingDirect_end          ! get_ixVarType('wLength')
  integer(i4b) :: ixSpectralIncomingDiffuse_start, ixSpectralIncomingDiffuse_end         ! get_ixVarType('wLength')
  integer(i4b) :: ixScalarCanopySunlitPAR = 7           ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyShadedPAR = 8           ! get_ixVarType('scalarv')
  integer(i4b) :: ixSpectralBelowCanopyDirect_start,ixSpectralBelowCanopyDirect_end       ! get_ixVarType('wLength')
  integer(i4b) :: ixSpectralBelowCanopyDiffuse_start,ixSpectralBelowCanopyDiffuse_end      ! get_ixVarType('wLength')
  integer(i4b) :: ixScalarBelowCanopySolar = 9          ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyAbsorbedSolar = 10       ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundAbsorbedSolar = 11       ! get_ixVarType('scalarv')
  ! longwave radiation
  integer(i4b) :: ixScalarLWRadCanopy = 12               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadGround = 13               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadUbound2Canopy = 14        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadUbound2Ground = 15        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadUbound2Ubound = 16        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadCanopy2Ubound = 17        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadCanopy2Ground = 18        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadCanopy2Canopy = 19        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadGround2Ubound = 20        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadGround2Canopy = 21        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWNetCanopy = 22               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWNetGround = 23               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWNetUbound = 24               !get_ixVarType('scalarv')
  ! turbulent heat transfer
  integer(i4b) :: ixScalarEddyDiffusCanopyTop = 25       !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarFrictionVelocity = 26          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarWindspdCanopyTop = 27          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarWindspdCanopyBottom = 28       !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundResistance = 29          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyResistance = 30          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLeafResistance = 31            !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSoilResistance = 32            !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSenHeatTotal = 33              !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSenHeatCanopy = 34             !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSenHeatGround = 35             !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatTotal = 36              !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatCanopyEvap = 37         !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatCanopyTrans = 38        !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatGround = 39             !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyAdvectiveHeatFlux = 40   !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundAdvectiveHeatFlux = 41   !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopySublimation = 42         !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSnowSublimation = 43           !, get_ixVarType('scalarv')
  ! liquid water fluxes associated with evapotranspiration
  integer(i4b) :: ixScalarStomResistSunlit = 44          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarStomResistShaded = 45          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarPhotosynthesisSunlit = 46      !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarPhotosynthesisShaded = 47      !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyTranspiration = 48       !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyEvaporation = 49         !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundEvaporation = 50         !get_ixVarType('scalarv')
  integer(i4b) :: ixMLayerTranspire_start, ixMLayerTranspire_end                 !get_ixVarType('midSoil')
  ! liquid and solid water fluxes through the canopy
  integer(i4b) :: ixScalarThroughfallSnow = 51 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarThroughfallRain = 52 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopySnowUnloading = 53 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyLiqDrainage = 54 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyMeltFreeze = 55 ! get_ixVarType('scalarv')
  ! energy fluxes and for the snow and soil domains
  integer(i4b) :: ixILayerConductiveFlux_start,ixILayerConductiveFlux_end ! get_ixVarType('ifcToto')
  integer(i4b) :: ixILayerAdvectiveFlux_start, ixILayerAdvectiveFlux_end ! get_ixVarType('ifcToto')
  integer(i4b) :: ixILayerNrgFlux_start, ixILayerNrgFlux_end ! get_ixVarType('ifcToto')
  integer(i4b) :: ixMLayerNrgFlux_start, ixMLayerNrgFlux_end ! get_ixVarType('midToto')
  ! liquid water fluxes for the snow domain
  integer(i4b) :: ixScalarSnowDrainage = 56! get_ixVarType('scalarv')
  integer(i4b) :: ixILayerLiqFluxSnow_start, ixILayerLiqFluxSnow_end! get_ixVarType('ifcSnow')
  integer(i4b) :: ixMLayerLiqFluxSnow_start, ixMLayerLiqFluxSnow_end! get_ixVarType('midSnow')
  ! liquid water fluxes for the soil domain
  integer(i4b) ::ixScalarRainPlusMelt = 57!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarMaxInfilRate = 58!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarInfiltration = 59!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarExfiltration = 60!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarSurfaceRunoff = 61!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarSurfaceRunoff_IE = 62
  integer(i4b) :: ixScalarSurfaceRunoff_SE = 63
  integer(i4b) :: ixMLayerSatHydCondMP_start, ixMLayerSatHydCondMP_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerSatHydCond_start, ixMLayerSatHydCond_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixILayerSatHydCond_start, ixILayerSatHydCond_end!, get_ixVarType('ifcSoil'), 
  integer(i4b) :: ixMLayerHydCond_start, ixMLayerHydCond_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixILayerLiqFluxSoil_start, ixILayerLiqFluxSoil_end!, get_ixVarType('ifcSoil'), 
  integer(i4b) :: ixMLayerLiqFluxSoil_start, ixMLayerLiqFluxSoil_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerBaseflow_start, ixMLayerBaseflow_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerColumnInflow_start, ixMLayerColumnInflow_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerColumnOutflow_start, ixMLayerColumnOutflow_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixScalarSoilBaseflow = 64!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarSoilDrainage = 65!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarAquiferRecharge = 66!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarAquiferTranspire = 67!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarAquiferBaseflow = 68!, get_ixVarType('scalarv'), 
  ! derived variables
  integer(i4b) :: ixScalarTotalET = 69
  integer(i4b) :: ixScalarTotalRunoff = 70
  integer(i4b) :: ixScalarNetRadiation = 71
  integer(i4b) :: numScalarFluxData = 71
  integer(i4b) :: numFluxData


end type flux_count_device
type,public:: flux2state_device !!!!
  integer(i4b),device,allocatable :: state1(:)
  integer(i4b),device,allocatable :: state2(:)
  integer(i4b) :: ixScalarCanairNetNrgFlux = 1          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyNetNrgFlux = 2          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundNetNrgFlux = 3          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyNetLiqFlux = 4          !get_ixVarType('scalarv')
  ! forcing
  integer(i4b) :: ixScalarRainfall = 5                  ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSnowfall = 6                  ! get_ixVarType('scalarv')
  ! shortwave radiation
  integer(i4b) :: ixSpectralIncomingDirect_start, ixSpectralIncomingDirect_end          ! get_ixVarType('wLength')
  integer(i4b) :: ixSpectralIncomingDiffuse_start, ixSpectralIncomingDiffuse_end         ! get_ixVarType('wLength')
  integer(i4b) :: ixScalarCanopySunlitPAR = 7           ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyShadedPAR = 8           ! get_ixVarType('scalarv')
  integer(i4b) :: ixSpectralBelowCanopyDirect_start,ixSpectralBelowCanopyDirect_end       ! get_ixVarType('wLength')
  integer(i4b) :: ixSpectralBelowCanopyDiffuse_start,ixSpectralBelowCanopyDiffuse_end      ! get_ixVarType('wLength')
  integer(i4b) :: ixScalarBelowCanopySolar = 9          ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyAbsorbedSolar = 10       ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundAbsorbedSolar = 11       ! get_ixVarType('scalarv')
  ! longwave radiation
  integer(i4b) :: ixScalarLWRadCanopy = 12               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadGround = 13               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadUbound2Canopy = 14        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadUbound2Ground = 15        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadUbound2Ubound = 16        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadCanopy2Ubound = 17        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadCanopy2Ground = 18        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadCanopy2Canopy = 19        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadGround2Ubound = 20        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWRadGround2Canopy = 21        !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWNetCanopy = 22               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWNetGround = 23               !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLWNetUbound = 24               !get_ixVarType('scalarv')
  ! turbulent heat transfer
  integer(i4b) :: ixScalarEddyDiffusCanopyTop = 25       !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarFrictionVelocity = 26          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarWindspdCanopyTop = 27          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarWindspdCanopyBottom = 28       !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundResistance = 29          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyResistance = 30          !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLeafResistance = 31            !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSoilResistance = 32            !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSenHeatTotal = 33              !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSenHeatCanopy = 34             !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSenHeatGround = 35             !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatTotal = 36              !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatCanopyEvap = 37         !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatCanopyTrans = 38        !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarLatHeatGround = 39             !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyAdvectiveHeatFlux = 40   !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundAdvectiveHeatFlux = 41   !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopySublimation = 42         !, get_ixVarType('scalarv')
  integer(i4b) :: ixScalarSnowSublimation = 43           !, get_ixVarType('scalarv')
  ! liquid water fluxes associated with evapotranspiration
  integer(i4b) :: ixScalarStomResistSunlit = 44          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarStomResistShaded = 45          !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarPhotosynthesisSunlit = 46      !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarPhotosynthesisShaded = 47      !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyTranspiration = 48       !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyEvaporation = 49         !get_ixVarType('scalarv')
  integer(i4b) :: ixScalarGroundEvaporation = 50         !get_ixVarType('scalarv')
  integer(i4b) :: ixMLayerTranspire_start, ixMLayerTranspire_end                 !get_ixVarType('midSoil')
  ! liquid and solid water fluxes through the canopy
  integer(i4b) :: ixScalarThroughfallSnow = 51 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarThroughfallRain = 52 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopySnowUnloading = 53 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyLiqDrainage = 54 ! get_ixVarType('scalarv')
  integer(i4b) :: ixScalarCanopyMeltFreeze = 55 ! get_ixVarType('scalarv')
  ! energy fluxes and for the snow and soil domains
  integer(i4b) :: ixILayerConductiveFlux_start,ixILayerConductiveFlux_end ! get_ixVarType('ifcToto')
  integer(i4b) :: ixILayerAdvectiveFlux_start, ixILayerAdvectiveFlux_end ! get_ixVarType('ifcToto')
  integer(i4b) :: ixILayerNrgFlux_start, ixILayerNrgFlux_end ! get_ixVarType('ifcToto')
  integer(i4b) :: ixMLayerNrgFlux_start, ixMLayerNrgFlux_end ! get_ixVarType('midToto')
  ! liquid water fluxes for the snow domain
  integer(i4b) :: ixScalarSnowDrainage = 56! get_ixVarType('scalarv')
  integer(i4b) :: ixILayerLiqFluxSnow_start, ixILayerLiqFluxSnow_end! get_ixVarType('ifcSnow')
  integer(i4b) :: ixMLayerLiqFluxSnow_start, ixMLayerLiqFluxSnow_end! get_ixVarType('midSnow')
  ! liquid water fluxes for the soil domain
  integer(i4b) ::ixScalarRainPlusMelt = 57!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarMaxInfilRate = 58!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarInfiltration = 59!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarExfiltration = 60!, get_ixVarType('scalarv'), 
  integer(i4b) ::ixScalarSurfaceRunoff = 61!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarSurfaceRunoff_IE = 62
  integer(i4b) :: ixScalarSurfaceRunoff_SE = 63
  integer(i4b) :: ixMLayerSatHydCondMP_start, ixMLayerSatHydCondMP_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerSatHydCond_start, ixMLayerSatHydCond_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixILayerSatHydCond_start, ixILayerSatHydCond_end!, get_ixVarType('ifcSoil'), 
  integer(i4b) :: ixMLayerHydCond_start, ixMLayerHydCond_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixILayerLiqFluxSoil_start, ixILayerLiqFluxSoil_end!, get_ixVarType('ifcSoil'), 
  integer(i4b) :: ixMLayerLiqFluxSoil_start, ixMLayerLiqFluxSoil_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerBaseflow_start, ixMLayerBaseflow_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerColumnInflow_start, ixMLayerColumnInflow_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixMLayerColumnOutflow_start, ixMLayerColumnOutflow_end!, get_ixVarType('midSoil'), 
  integer(i4b) :: ixScalarSoilBaseflow = 64!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarSoilDrainage = 65!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarAquiferRecharge = 66!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarAquiferTranspire = 67!, get_ixVarType('scalarv'), 
  integer(i4b) :: ixScalarAquiferBaseflow = 68!, get_ixVarType('scalarv'), 
  ! derived variables
  integer(i4b) :: ixScalarTotalET = 69
  integer(i4b) :: ixScalarTotalRunoff = 70
  integer(i4b) :: ixScalarNetRadiation = 71
  integer(i4b) :: numScalarFluxData = 71
  integer(i4b) :: numFluxData


end type flux2state_device

type,public:: prog_data_device
  real(rkind) :: dt_init               ! get_ixVarType('scalarv')
  ! state variables for vegetation
  real(rkind),device,allocatable :: scalarCanopyIce(:)                 ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyLiq(:)                 ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyWat(:)                 ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanairTemp(:)                ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyTemp(:)                ! get_ixVarType('scalarv')
  ! state variables for snow
  real(rkind),device,allocatable :: spectralSnowAlbedoDiffuse(:,:)     ! get_ixVarType('wLength')
  real(rkind),device,allocatable :: scalarSnowAlbedo(:)                ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSnowDepth(:)                 ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSWE(:)                       ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSfcMeltPond(:)               ! get_ixVarType('scalarv')
  ! define state variables for the snow+soil domain
  real(rkind),device,allocatable :: mLayerTemp(:,:)                    ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: mLayerVolFracIce(:,:)              ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: mLayerVolFracLiq(:,:)              ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: mLayerVolFracWat(:,:)              ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: mLayerMatricHead(:,:)              ! get_ixVarType('midSoil')
  ! enthalpy
  real(rkind),device,allocatable :: scalarCanairEnthalpy(:)            ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyEnthalpy(:)            ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: mLayerEnthalpy(:,:)                ! get_ixVarType('midToto')
  ! other state variables
  real(rkind),device,allocatable :: scalarAquiferStorage(:)            ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSurfaceTemp(:)               ! get_ixVarType('scalarv')
  ! define coordinate variables
  real(rkind),device,allocatable :: mLayerDepth(:,:)                   ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: mLayerHeight(:,:)                  ! get_ixVarType('ifcToto')
  real(rkind),device,allocatable :: iLayerHeight(:,:)                  ! get_ixVarType('ifcToto')


end type prog_data_device

type,public:: deriv_data_device !!!
  real(rkind),device,allocatable :: dThermalC_dWatAbove_m(:,:)
  real(rkind),device,allocatable :: dThermalC_dWatBelow_m(:,:)
  real(rkind),device,allocatable :: dThermalC_dTempAbove_m(:,:)
  real(rkind),device,allocatable :: dThermalC_dTempBelow_m(:,:)
  real(rkind),device,allocatable :: d2Theta_dTkCanopy2(:)
  real(rkind),device,allocatable :: dVolTot_dPsi0_m(:,:)
  real(rkind),device,allocatable :: dPsiLiq_dPsi0_m(:,:), dPsiLiq_dTemp_m(:,:)
  real(rkind),device,allocatable :: mLayerdTheta_dTk_m(:,:)
  real(rkind),device,allocatable :: dTheta_dTkCanopy(:)
  real(rkind),device,allocatable :: dFracLiqWat_dTk_m(:,:)
  real(rkind),device,allocatable :: dFracLiqVeg_dTkCanopy(:)
  real(rkind),device,allocatable :: d2VolTot_dPsi02_m(:,:), mLayerd2Theta_dTk2_m(:,:)
  real(rkind),device,allocatable :: dCanairTemp_dEnthalpy(:), dCanopyTemp_dEnthalpy(:)
  real(rkind),device,allocatable :: dTemp_dEnthalpy_m(:,:)
  real(rkind),device,allocatable :: dCanopyTemp_dCanWat(:)
  real(rkind),device,allocatable :: dTemp_dTheta_m(:,:), dTemp_dPsi0_m(:,:)

  real(rkind),device,allocatable :: dGroundNetFlux_dGroundTemp(:)
  real(rkind),device,allocatable :: dNrgFlux_dTempAbove_m(:,:), dNrgFlux_dTempBelow_m(:,:), dNrgFlux_dWatAbove_m(:,:), dNrgFlux_dWatBelow_m(:,:)
  real(rkind),device,allocatable :: scalarThroughfallRainDeriv(:), scalarCanopyLiqDrainageDeriv(:), scalarCanopyLiqDeriv(:)

  real(rkind),device,allocatable :: dBaseflow_dAquifer(:)
  real(rkind),device,allocatable :: dAquiferTrans_dTCanair(:)
  real(rkind),device,allocatable :: dAquiferTrans_dTCanopy(:), dAquiferTrans_dTGround(:), dAquiferTrans_dCanWat(:)
  real(rkind),device,allocatable :: dCanopyTrans_dCanWat(:), dCanopyTrans_dTCanair(:), dCanopyTrans_dTCanopy(:), dCanopyTrans_dTGround(:)
  real(rkind),device,allocatable :: mLayerdTheta_dPsi_m(:,:)
  real(rkind),device,allocatable :: dq_dHydStateAbove_m(:,:), dq_dHydStateBelow_m(:,:),dq_dHydStateLayerSurfVec_m(:,:), dq_dNrgStateAbove_m(:,:), dq_dNrgStateBelow_m(:,:), dq_dNrgStateLayerSurfVec_m(:,:)
  real(rkind),device,allocatable :: mLayerdPsi_dTheta_m(:,:)
  real(rkind),device,allocatable :: mLayerdTrans_dTCanair_m(:,:), mLayerdTrans_dTCanopy_m(:,:), mLayerdTrans_dTGround_m(:,:), mLayerdTrans_dCanWat_m(:,:)
  real(rkind),device,allocatable :: dCanLiq_dTcanopy(:)
  real(rkind),device,allocatable :: iLayerLiqFluxSnowDeriv_m(:,:)
  real(rkind),device,allocatable :: dVolHtCapBulk_dPsi0_m(:,:)
  real(rkind),device,allocatable :: dVolHtCapBulk_dTheta_m(:,:)
  real(rkind),device,allocatable :: dVolHtCapBulk_dTk_m(:,:)
  real(rkind),device,allocatable :: dCm_dTk_m(:,:)
  real(rkind),device,allocatable :: dCm_dPsi0_m(:,:)
  real(rkind),device,allocatable :: dCanairNetFlux_dCanairTemp(:), dCanairNetFlux_dCanopyTemp(:), dCanairNetFlux_dGroundTemp(:)
  real(rkind),device,allocatable :: dCanopyNetFlux_dCanairTemp(:), dCanopyNetFlux_dCanopyTemp(:), dCanopyNetFlux_dGroundTemp(:)
  real(rkind),device,allocatable :: dGroundNetFlux_dCanairTemp(:), dGroundNetFlux_dCanopyTemp(:)
  real(rkind),device,allocatable :: dCanopyEvaporation_dCanWat(:), dCanopyEvaporation_dTCanair(:), dCanopyEvaporation_dTCanopy(:), dCanopyEvaporation_dTGround(:)
  real(rkind),device,allocatable :: dGroundEvaporation_dCanWat(:), dGroundEvaporation_dTCanair(:), dGroundEvaporation_dTCanopy(:), dGroundEvaporation_dTGround(:)
  real(rkind),device,allocatable :: dCanopyNetFlux_dCanWat(:), dGroundNetFlux_dCanWat(:)
  real(rkind),device,allocatable :: dCm_dTkCanopy(:)
  real(rkind),device,allocatable :: dVolHtCapBulk_dTkCanopy(:),  dVolHtCapBulk_dCanWat(:)
  real(rkind),device,allocatable :: dCompress_dPsi_m(:,:)
end type deriv_data_device

type,public:: indx_data_device
  integer(i4b),device,allocatable :: nSnow(:)
  integer(i4b) :: nSoil 
  integer(i4b),device,allocatable :: nLayers_d(:) !!!
  integer(i4b),device,allocatable :: layerType(:,:)
  ! ! number of state variables of different type
  integer(i4b),device,allocatable :: nCasNrg(:)
  integer(i4b),device,allocatable :: nVegNrg(:)
  integer(i4b),device,allocatable :: nVegMass(:)
  integer(i4b),device,allocatable :: nVegState(:)
  integer(i4b),device,allocatable :: nNrgState(:)
  integer(i4b),device,allocatable :: nWatState(:)
  integer(i4b),device,allocatable :: nMatState(:)
  integer(i4b),device,allocatable :: nMassState(:)
  integer(i4b),device,allocatable :: nState(:)
  ! ! number of state variables within different domains in the snow+soil system
  ! integer(i4b),device,allocatable :: nSnowSoilNrg_d(:)
  ! integer(i4b),device,allocatable :: nSnowOnlyNrg_d(:)
  ! integer(i4b),device,allocatable :: nSoilOnlyNrg_d(:)
  ! integer(i4b),device,allocatable :: nSnowSoilHyd_d(:)
  ! integer(i4b),device,allocatable :: nSnowOnlyHyd_d(:)
  ! integer(i4b),device,allocatable :: nSoilOnlyHyd_d(:)
  ! type of model state variables
  integer(i4b),device,allocatable :: ixControlVolume(:,:)
  integer(i4b),device,allocatable :: ixDomainType(:,:)
  integer(i4b),device,allocatable :: ixStateType(:,:)
  integer(i4b),device,allocatable :: ixHydType(:,:)
  ! ! ! type of model state variables (state subset)
  ! integer(i4b),device,allocatable :: ixDomainType_subset(:,:)
  ! integer(i4b),device,allocatable :: ixStateType_subset(:,:)
  ! ! ! mapping between state subset and the full state vector
  ! integer(i4b),device,allocatable :: ixMapFull2Subset(:,:)
  ! integer(i4b),device,allocatable :: ixMapSubset2Full(:,:)
  ! indices of model specific state variables
  integer(i4b),device,allocatable :: ixCasNrg(:)
  integer(i4b),device,allocatable :: ixVegNrg(:)
  integer(i4b),device,allocatable :: ixVegHyd(:)
  integer(i4b),device,allocatable :: ixTopNrg(:)
  integer(i4b),device,allocatable :: ixTopHyd(:)
  integer(i4b),device,allocatable :: ixAqWat(:)
  ! vectors of indices for specific state types
  integer(i4b),device,allocatable :: ixNrgOnly(:,:)
  ! integer(i4b),device,allocatable :: ixHydOnly(:,:)
  ! integer(i4b),device,allocatable :: ixMatOnly(:,:)
  ! integer(i4b),device,allocatable :: ixMassOnly(:,:)
  ! vectors of indices for specific state types within specific sub-domains
  integer(i4b),device,allocatable :: ixSnowSoilNrg(:,:)
  integer(i4b),device,allocatable :: ixSnowOnlyNrg(:,:)
  integer(i4b),device,allocatable :: ixSoilOnlyNrg(:,:)
  integer(i4b),device,allocatable :: ixSnowSoilHyd(:,:)
  integer(i4b),device,allocatable :: ixSnowOnlyHyd(:,:)
  integer(i4b),device,allocatable :: ixSoilOnlyHyd(:,:)
  ! vectors of indices for specfic state types within specific sub-domains
  integer(i4b),device,allocatable :: ixNrgCanair(:)
  integer(i4b),device,allocatable :: ixNrgCanopy(:)
  integer(i4b),device,allocatable :: ixHydCanopy(:)
  integer(i4b),device,allocatable :: ixNrgLayer(:,:)
  integer(i4b),device,allocatable :: ixHydLayer(:,:)
  integer(i4b),device,allocatable :: ixWatAquifer(:)

  ! vectors of indices for specific state types IN SPECIFIC SUB-DOMAINS
  ! ! indices within state vectors
  integer(i4b),device,allocatable :: ixAllState(:,:)
  integer(i4b),device,allocatable :: ixSoilState(:,:)
  integer(i4b),device,allocatable :: ixLayerState(:,:)

  ! integer(i4b),device,allocatable :: ixNrgLayer(:,:), ixHydLayer(:,:)
  ! integer(i4b) :: nSoil!, nLayers
  ! integer(i4b),device,allocatable :: nSnow(:)!,nSoil_d
  ! integer(i4b),device,allocatable :: nLayers_d(:)
  ! integer(i4b),device,allocatable :: ixCasNrg(:), ixVegNrg(:), ixVegHyd(:), ixAqWat(:)
  ! integer(i4b),device,allocatable :: ixTopHyd(:), ixTopNrg(:)
  ! integer(i4b),device,allocatable :: ixSnowSoilNrg(:,:)
  ! integer(i4b),device,allocatable :: ixSnowSoilHyd(:,:) 
  ! ! integer(i4b) :: nSnowSoilNrg, nSnowSoilHyd
  ! integer(i4b),device,allocatable :: ixHydType(:,:)
  ! integer(i4b),device,allocatable :: ixStateType(:,:)
  ! ! integer(i4b),device,allocatable :: ixStateType_subset(:,:)
  ! integer(i4b),device,allocatable :: ixHydCanopy(:)
  ! integer(i4b),device,allocatable :: ixNrgCanopy(:)
  ! ! integer(i4b) :: ixNrgCanair
  ! ! integer(i4b),device,allocatable :: ixMapSubset2Full(:,:)!, ixMapFull2Subset(:,:)
  ! !integer(i4b),device,allocatable :: ixDomainType_subset(:,:)
  ! integer(i4b),device,allocatable :: ixControlVolume(:,:)
  ! integer(i4b),device,allocatable :: layerType(:,:)
  ! integer(i4b),device,allocatable :: ixSoilOnlyHyd(:,:)
  ! integer(i4b),device,allocatable :: ixSnowOnlyHyd(:,:)
  ! integer(i4b),device,allocatable :: ixSoilOnlyNrg(:,:)
  ! integer(i4b),device,allocatable :: ixSnowOnlyNrg(:,:)

  ! integer(i4b),device,allocatable :: ixSoilState(:,:)
  ! ! integer(i4b),device,allocatable :: ixLayerState(:)
  integer(i4b) :: numberFluxCalc
  ! integer(i4b),device,allocatable :: ixDomainType(:,:)
  ! integer(i4b),device,allocatable :: ixNrgOnly(:,:)

  ! Not used:
  ! ixMatricHead(:,:)
  ! ixVolFracWat(:,:)
  ! ixLayerActive(:,:)
  integer(i4b) :: numberDomainSplitMass, numberDomainSplitNrg, numberScalarSolutions, numberStateSplit

end type indx_data_device


type,public:: forc_data_device
  real(rkind),device,allocatable :: airtemp_d(:)
  real(rkind),device,allocatable :: windspd_d(:)
  real(rkind),device,allocatable :: airpres_d(:)
  real(rkind),device,allocatable :: LWRadAtm_d(:)

  !!! UNUSED (so far)
  real(rkind),device,allocatable :: time(:)
  real(rkind),device,allocatable :: pptrate(:)
  real(rkind),device,allocatable :: SWRadAtm(:)
  real(rkind),device,allocatable :: spechum(:)
end type forc_data_device

type,public :: bvar_data_device
  real(rkind),device,allocatable :: basin__SurfaceRunoff(:)
real(rkind),device,allocatable :: basin__ColumnOutflow(:)
real(rkind),device,allocatable :: basin__AquiferStorage(:)
real(rkind),device,allocatable :: basin__AquiferRecharge(:)
real(rkind),device,allocatable :: basin__AquiferBaseflow(:)
real(rkind),device,allocatable :: basin__AquiferTranspire(:)
real(rkind),device,allocatable :: basin__TotalRunoff(:)
real(rkind),device,allocatable :: basin__SoilDrainage(:)
real(rkind),device,allocatable :: basin__totalArea(:)
real(rkind),device,allocatable :: routingRunoffFuture(:,:)
real(rkind),device,allocatable :: routingFractionFuture(:,:)
real(rkind),device,allocatable :: averageInstantRunoff(:)
real(rkind),device,allocatable :: averageRoutedRunoff(:)
end type bvar_data_device

type,public :: decisions_device
  logical(lgt),device,allocatable :: true
  logical(lgt),device,allocatable :: false
  integer(i4b),device,allocatable :: one
  integer(i4b),device,allocatable :: soilCatTbl
integer(i4b),device,allocatable :: vegeParTbl
integer(i4b),device,allocatable :: soilStress
integer(i4b),device,allocatable :: stomResist
integer(i4b),device,allocatable :: bbTempFunc
integer(i4b),device,allocatable :: bbHumdFunc
integer(i4b),device,allocatable :: bbElecFunc
integer(i4b),device,allocatable :: bbCO2point
integer(i4b),device,allocatable :: bbNumerics
integer(i4b),device,allocatable :: bbAssimFnc
integer(i4b),device,allocatable :: bbCanIntg8
integer(i4b),device,allocatable :: num_method
integer(i4b),device,allocatable :: fDerivMeth
integer(i4b),device,allocatable :: LAI_method
integer(i4b),device,allocatable :: cIntercept
integer(i4b),device,allocatable :: f_Richards
integer(i4b),device,allocatable :: groundwatr
integer(i4b),device,allocatable :: hc_profile
integer(i4b),device,allocatable :: bcUpprTdyn
integer(i4b),device,allocatable :: bcLowrTdyn
integer(i4b),device,allocatable :: bcUpprSoiH
integer(i4b),device,allocatable :: bcLowrSoiH
integer(i4b),device,allocatable :: veg_traits
integer(i4b),device,allocatable :: rootProfil
integer(i4b),device,allocatable :: canopyEmis
integer(i4b),device,allocatable :: snowIncept
integer(i4b),device,allocatable :: snowUnload
integer(i4b),device,allocatable :: windPrfile
integer(i4b),device,allocatable :: astability
integer(i4b),device,allocatable :: canopySrad
integer(i4b),device,allocatable :: alb_method
integer(i4b),device,allocatable :: snowLayers
integer(i4b),device,allocatable :: compaction
integer(i4b),device,allocatable :: thCondSnow
integer(i4b),device,allocatable :: thCondSoil
integer(i4b),device,allocatable :: spatial_gw
integer(i4b),device,allocatable :: subRouting
integer(i4b),device,allocatable :: snowDenNew
integer(i4b),device,allocatable :: nrgConserv
integer(i4b),device,allocatable :: aquiferIni
integer(i4b),device,allocatable :: urbanVegCategory
integer(i4b),device,allocatable :: infRateMax
integer(i4b),device,allocatable :: surfRun_IE
integer(i4b),device,allocatable :: surfRun_SE

real(rkind),device,allocatable :: data_step
real(rkind),device,allocatable :: refJulDay
integer(i4b),device,allocatable :: NC_TIME_ZONE

end type decisions_device

type,public :: type_data_device
  integer(i4b),device,allocatable :: vegTypeIndex(:)
  integer(i4b),device,allocatable :: soilTypeIndex(:)
  integer(i4b),device,allocatable :: slopeTypeIndex(:)
  integer(i4b),device,allocatable :: downHRUindex(:)
end type type_data_device

type,public :: veg_parameters
  real(rkind),device,allocatable :: bp(:)
  real(rkind),device,allocatable :: mp(:)
  real(rkind),device,allocatable :: c3psn(:)
  real(rkind),device,allocatable :: avcmx(:)
  real(rkind),device,allocatable :: vcmx25(:)
  real(rkind),device,allocatable :: ako(:)
  real(rkind),device,allocatable :: ko25(:)
  real(rkind),device,allocatable :: akc(:)
  real(rkind),device,allocatable :: kc25(:)
  real(rkind),device,allocatable :: qe25(:)
  real(rkind),device,allocatable :: folnmx(:)
  real(rkind),device,allocatable :: rgl_d(:)
  real(rkind),device,allocatable :: rsmin_d(:)
  real(rkind),device,allocatable :: rsmax_d(:)
  real(rkind),device,allocatable :: topt_d(:)
  real(rkind),device,allocatable :: hs_d(:)
  real(rkind),device,allocatable :: rhol(:,:)
  real(rkind),device,allocatable :: rhos(:,:)
  real(rkind),device,allocatable :: taul(:,:)
  real(rkind),device,allocatable :: taus(:,:)
  integer(i4b),device,allocatable :: opt_alb
  real(rkind),device,allocatable :: omegas(:)
  integer(i4b),device,allocatable :: opt_rad
  real(rkind),device,allocatable :: tfrz
  real(rkind),device,allocatable :: betais
  real(rkind),device,allocatable :: betads
  real(rkind),device,allocatable :: xl(:)
  real(rkind),device,allocatable :: hvt_(:,:)
  real(rkind),device,allocatable :: hvb_(:,:)
  real(rkind),device,allocatable :: rc(:)
  real(rkind),device,allocatable :: swemx
  real(rkind),device,allocatable :: albsat(:,:)
  real(rkind),device,allocatable :: albdry(:,:)
  real(rkind),device,allocatable :: alblak(:)
  integer(i4b),device,allocatable :: isbarren
integer(i4b),device,allocatable :: ISSNOW
integer(i4b),device,allocatable :: ISWATER
real(rkind),device,allocatable :: saim_(:,:,:)
real(rkind),device,allocatable :: laim_(:,:,:)
integer(i4b),device,allocatable :: dveg
real(rkind),device,allocatable :: tmin(:)
end type
type, public :: veg_param_tables
real(rkind),device,allocatable :: rgltbl(:), rstbl(:), hstbl(:)
real(rkind),device,allocatable :: TOPT_DATA,RSMAX_DATA
end type

type,public :: attr_data_device
  real(rkind),device,allocatable :: latitude(:)
  real(rkind),device,allocatable :: longitude(:)
  real(rkind),device,allocatable :: elevation(:)
  real(rkind),device,allocatable :: tan_slope(:)
  real(rkind),device,allocatable :: contourLength(:)
  real(rkind),device,allocatable :: HRUarea(:)
  real(rkind),device,allocatable :: mHeight(:)
  real(rkind),device,allocatable :: aspect(:)
end type

type,public :: zLookup_device
real(rkind),device,allocatable :: temperature(:,:,:), psiLiq_int(:,:,:), deriv2(:,:,:)
real(rkind),device,allocatable :: h_lookup(:), t_lookup(:)

end type

end module device_data_types