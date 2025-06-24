module device_data_types
use cudafor 
USE nrtype

type,public :: mpar_data_device

  real(rkind),device,allocatable :: soil_dens_intr(:)
  real(rkind),device,allocatable :: thCond_soil(:)
  real(rkind),device,allocatable :: frac_sand(:)
  real(rkind),device,allocatable :: frac_clay(:)
  real(rkind),device,allocatable :: frac_silt(:)
  real(rkind),device,allocatable :: snowfrz_scale
  real(rkind),device,allocatable :: specificHeatVeg
  real(rkind),device,allocatable :: maxMassVegetation
  real(rkind),device,allocatable :: Fcapil, k_snow, mw_exp
  real(rkind),device,allocatable :: fixedThermalCond_snow

  real(rkind),device,allocatable :: heightCanopyTop
  real(rkind),device,allocatable :: heightCanopyBottom
  real(rkind),device,allocatable :: canopyWettingFactor
  real(rkind),device,allocatable :: canopyWettingExp

  real(rkind),device,allocatable :: z0Snow
  real(rkind),device,allocatable :: z0Soil
  real(rkind),device,allocatable :: z0Canopy
  real(rkind),device,allocatable :: zpdFraction
  real(rkind),device,allocatable :: critRichNumber
  real(rkind),device,allocatable :: Louis79_bparam
  real(rkind),device,allocatable :: Louis79_cStar
  real(rkind),device,allocatable :: Mahrt87_eScale
  real(rkind),device,allocatable :: windReductionParam
  real(rkind),device,allocatable :: leafExchangeCoeff
  real(rkind),device,allocatable :: leafDimension
  real(rkind),device,allocatable :: vGn_n(:), vGn_alpha(:)
  real(rkind),device,allocatable :: theta_sat(:)
  real(rkind),device,allocatable :: theta_res(:)
  real(rkind),device,allocatable :: plantWiltPsi
  real(rkind),device,allocatable :: soilStressParam
  real(rkind),device,allocatable :: critSoilWilting
  real(rkind),device,allocatable :: critSoilTranspire
  real(rkind),device,allocatable :: critAquiferTranspire
  real(rkind),device,allocatable :: minStomatalResistance

  real(rkind),device,allocatable :: upperBoundTemp,lowerBoundTemp

  real(rkind),device,allocatable :: throughfallScaleRain, canopyDrainageCoeff

  real(rkind),device,allocatable :: aquiferBaseflowRate, aquiferScaleFactor, aquiferBaseflowExp
  real(rkind),device,allocatable :: Kc25
  real(rkind),device,allocatable :: Ko25
  real(rkind),device,allocatable :: Kc_qFac
  real(rkind),device,allocatable :: Ko_qFac
  real(rkind),device,allocatable :: kc_Ha
  real(rkind),device,allocatable :: ko_Ha
  real(rkind),device,allocatable :: vcmax25_canopyTop
  real(rkind),device,allocatable :: vcmax_qFac
  real(rkind),device,allocatable :: vcmax_Ha
  real(rkind),device,allocatable :: vcmax_Hd
  real(rkind),device,allocatable :: vcmax_Sv
  real(rkind),device,allocatable :: vcmax_Kn
  real(rkind),device,allocatable :: jmax25_scale
  real(rkind),device,allocatable :: jmax_Ha
  real(rkind),device,allocatable :: jmax_Hd
  real(rkind),device,allocatable :: jmax_Sv
  real(rkind),device,allocatable :: fractionJ
  real(rkind),device,allocatable :: quantamYield
  real(rkind),device,allocatable :: vpScaleFactor
  real(rkind),device,allocatable :: cond2photo_slope
  real(rkind),device,allocatable :: minStomatalConductance
  real(rkind),device,allocatable :: upperBoundHead, upperBoundTheta, lowerBoundHead, lowerBoundTheta
  real(rkind),device,allocatable :: wettingFrontSuction,rootingDepth,kAnisotropic,zScale_TOPMODEL,qSurfScale,f_impede,soilIceScale,soilIceCV
  real(rkind),device,allocatable :: theta_mp, mpExp
  real(rkind),device,allocatable :: specificStorage
  real(rkind),device,allocatable :: fieldCapacity
  real(rkind) :: idaMaxOrder
  real(rkind) :: idaMaxInternalSteps
  real(rkind) :: idaInitStepSize
  real(rkind) :: idaMinStepSize
  real(rkind) :: idaMaxErrTestFail

  real(rkind) :: absTolTempCas
  real(rkind) :: relTolTempCas
  real(rkind) :: absTolTempVeg
  real(rkind) :: relTolTempVeg
  real(rkind) :: absTolWatVeg
  real(rkind) :: relTolWatVeg
  real(rkind) :: absTolTempSoilSnow
  real(rkind) :: relTolTempSoilSnow
  real(rkind) :: absTolWatSnow
  real(rkind) :: relTolWatSnow
  real(rkind) :: absTolMatric
  real(rkind) :: relTolMatric
  real(rkind) :: absTolAquifr
  real(rkind) :: relTolAquifr
  real(rkind) :: maxstep
  real(rkind) :: be_steps

  !!! UNUSED (so far)
  real(rkind),device,allocatable :: tempCritRain
  real(rkind),device,allocatable :: tempRangeTimestep
  real(rkind),device,allocatable :: frozenPrecipMultip
  real(rkind),device,allocatable :: albedoMinWinter
  real(rkind),device,allocatable :: albedoMinSpring
  real(rkind),device,allocatable :: albedoMinVisible
  real(rkind),device,allocatable :: albedoMinNearIR
  real(rkind),device,allocatable :: albedoDecayRate
  real(rkind),device,allocatable :: albedoSootLoad
  real(rkind),device,allocatable :: albedoRefresh
  ! real(rkind) :: radExt_snow
  real(rkind),device,allocatable :: directScale
  real(rkind),device,allocatable :: Frad_direct
  real(rkind),device,allocatable :: newSnowDenMin
  real(rkind),device,allocatable :: newSnowDenMult
  real(rkind),device,allocatable :: newSnowDenScal
  real(rkind),device,allocatable :: constSnowDen
  real(rkind),device,allocatable :: newSnowDenAdd
  real(rkind),device,allocatable :: newSnowDenMultTemp
  real(rkind),device,allocatable :: newSnowDenMultWind
  real(rkind),device,allocatable :: newSnowDenMultAnd
  real(rkind),device,allocatable :: newSnowDenBase
  ! real(rkind) :: densScalGrowth
  ! real(rkind) :: tempScalGrowth
  ! real(rkind) :: grainGrowthRate
  ! real(rkind) :: densScalOvrbdn
  ! real(rkind) :: tempScalOvrbdn
  ! real(rkind) :: baseViscosity
  real(rkind),device,allocatable :: winterSAI
  real(rkind),device,allocatable :: summerLAI
  real(rkind),device,allocatable :: rootScaleFactor1
  real(rkind),device,allocatable :: rootScaleFactor2
  real(rkind),device,allocatable :: rootDistExp
  ! real(rkind) :: throughfallScaleSnow
  real(rkind),device,allocatable :: refInterceptCapSnow
  real(rkind),device,allocatable :: refInterceptCapRain
  real(rkind),device,allocatable :: snowUnloadingCoeff
  real(rkind),device,allocatable :: ratioDrip2Unloading
  real(rkind),device,allocatable :: minTempUnloading
  real(rkind),device,allocatable :: rateTempUnloading
  real(rkind),device,allocatable :: minWindUnloading
  real(rkind),device,allocatable :: rateWindUnloading

  real(rkind),device,allocatable :: compactedDepth
  ! real(rkind) :: specificYield
  real(rkind),device,allocatable :: minwind
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
  real(rkind),device,allocatable :: zmin
  real(rkind) :: zminLayer1
  real(rkind) :: zminLayer2
  real(rkind) :: zminLayer3
  real(rkind) :: zminLayer4
  real(rkind) :: zminLayer5

  real(rkind),device,allocatable :: k_soil(:)
  real(rkind),device,allocatable :: k_macropore(:)
  real(rkind),allocatable,device :: densScalGrowth
  real(rkind),allocatable,device :: tempScalGrowth
  real(rkind),allocatable,device :: densScalOvrbdn
  real(rkind),allocatable,device :: tempScalOvrbdn
  real(rkind),allocatable,device :: baseViscosity
  real(rkind),allocatable,device :: grainGrowthRate
  real(rkind) :: zmaxLayer1_lower
  real(rkind) :: zmaxLayer2_lower
  real(rkind) :: zmaxLayer3_lower
  real(rkind) :: zmaxLayer4_lower
  real(rkind) :: zmaxLayer1_upper
  real(rkind) :: zmaxLayer2_upper
  real(rkind) :: zmaxLayer3_upper
  real(rkind) :: zmaxLayer4_upper
  real(rkind),device,allocatable :: zmax
  real(rkind),device,allocatable :: Frad_vis
  real(rkind),device,allocatable :: albedoMaxNearIR
  real(rkind),device,allocatable :: albedoMaxVisible
  real(rkind),device,allocatable :: albedoMax

end type mpar_data_device

type, public :: diag_data_device !!!!

  ! local properties
  real(rkind),device,allocatable :: scalarCanopyDepth(:)               ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarGreenVegFraction(:)          ! get_ixVarType('scalarv')
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
  real(rkind),device,allocatable :: mLayerVolHtCapBulk_m(:,:)              ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: mLayerCm_m(:,:)                        ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: scalarLambda_drysoil(:)            ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLambda_wetsoil(:)            ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: mLayerThermalC_m(:,:)                  ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: iLayerThermalC_m(:,:)                  ! get_ixVarType('ifcToto')
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
  real(rkind),device,allocatable :: mLayerTranspireLim_m(:,:)              ! get_ixVarType('midSoil')
  real(rkind),device,allocatable :: mLayerRootDensity_m(:,:)               ! get_ixVarType('midSoil')
  real(rkind),device,allocatable :: scalarAquiferRootFrac(:)           ! get_ixVarType('scalarv')
  ! canopy hydrology
  real(rkind),device,allocatable :: scalarFracLiqVeg(:)                ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyWetFraction(:)         ! get_ixVarType('scalarv')
  ! snow hydrology
  real(rkind),device,allocatable :: scalarSnowAge(:)                   ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarGroundSnowFraction(:)        ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: spectralSnowAlbedoDirect(:,:)        ! get_ixVarType('wLength')
  real(rkind),device,allocatable :: mLayerFracLiqSnow_m(:,:)               ! get_ixVarType('midSnow')
  real(rkind),device,allocatable :: mLayerThetaResid_m(:,:)                ! get_ixVarType('midSnow')
  real(rkind),device,allocatable :: mLayerPoreSpace_m(:,:)                 ! get_ixVarType('midSnow')
  real(rkind),device,allocatable :: mLayerMeltFreeze(:,:)                ! get_ixVarType('midToto')
  ! soil hydrology
  real(rkind),device,allocatable :: scalarInfilArea(:)                 ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarFrozenArea(:)                ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSoilControl(:)               ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: mLayerVolFracAir_m(:,:)                ! get_ixVarType('midToto')
  real(rkind),device,allocatable :: mLayerTcrit(:,:)                     ! get_ixVarType('midSoil')
  real(rkind),device,allocatable :: mLayerCompress_m(:,:)                  ! get_ixVarType('midSoil')
  real(rkind),device,allocatable :: scalarSoilCompress(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: mLayerMatricHeadLiq(:,:)             ! get_ixVarType('midSoil')
  real(rkind),device,allocatable :: scalarTotalSoilLiq(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarTotalSoilIce(:)              ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarTotalSoilWat(:)              ! get_ixVarType('scalarv')
  ! variable shortcuts
  real(rkind),device,allocatable :: scalarVGn_m_m(:,:)                     ! get_ixVarType('midSoil')
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
  real(rkind),device,allocatable :: scalarCanairNetNrgFlux(:)          !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyNetNrgFlux(:)          !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarGroundNetNrgFlux(:)          !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyNetLiqFlux(:)          !get_ixVarType('scalarv')
  ! forcing
  real(rkind),device,allocatable :: scalarRainfall(:)                  ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSnowfall(:)                  ! get_ixVarType('scalarv')
  ! shortwave radiation
  real(rkind),device,allocatable :: spectralIncomingDirect(:,:)          ! get_ixVarType('wLength')
  real(rkind),device,allocatable :: spectralIncomingDiffuse(:,:)         ! get_ixVarType('wLength')
  real(rkind),device,allocatable :: scalarCanopySunlitPAR(:)           ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyShadedPAR(:)           ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: spectralBelowCanopyDirect(:,:)       ! get_ixVarType('wLength')
  real(rkind),device,allocatable :: spectralBelowCanopyDiffuse(:,:)      ! get_ixVarType('wLength')
  real(rkind),device,allocatable :: scalarBelowCanopySolar(:)          ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyAbsorbedSolar(:)       ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarGroundAbsorbedSolar(:)       ! get_ixVarType('scalarv')
  ! longwave radiation
  real(rkind),device,allocatable :: scalarLWRadCanopy(:)               !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLWRadGround(:)               !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLWRadUbound2Canopy(:)        !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLWRadUbound2Ground(:)        !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLWRadUbound2Ubound(:)        !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLWRadCanopy2Ubound(:)        !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLWRadCanopy2Ground(:)        !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLWRadCanopy2Canopy(:)        !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLWRadGround2Ubound(:)        !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLWRadGround2Canopy(:)        !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLWNetCanopy(:)               !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLWNetGround(:)               !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLWNetUbound(:)               !get_ixVarType('scalarv')
  ! turbulent heat transfer
  real(rkind),device,allocatable :: scalarEddyDiffusCanopyTop(:)       !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarFrictionVelocity(:)          !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarWindspdCanopyTop(:)          !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarWindspdCanopyBottom(:)       !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarGroundResistance(:)          !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyResistance(:)          !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLeafResistance(:)            !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSoilResistance(:)            !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSenHeatTotal(:)              !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSenHeatCanopy(:)             !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSenHeatGround(:)             !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLatHeatTotal(:)              !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLatHeatCanopyEvap(:)         !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLatHeatCanopyTrans(:)        !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarLatHeatGround(:)             !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyAdvectiveHeatFlux(:)   !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarGroundAdvectiveHeatFlux(:)   !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopySublimation(:)         !, get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarSnowSublimation(:)           !, get_ixVarType('scalarv')
  ! liquid water fluxes associated with evapotranspiration
  real(rkind),device,allocatable :: scalarStomResistSunlit(:)          !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarStomResistShaded(:)          !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarPhotosynthesisSunlit(:)      !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarPhotosynthesisShaded(:)      !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyTranspiration(:)       !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyEvaporation(:)         !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarGroundEvaporation(:)         !get_ixVarType('scalarv')
  real(rkind),device,allocatable :: mLayerTranspire_m(:,:)                 !get_ixVarType('midSoil')
  ! liquid and solid water fluxes through the canopy
  real(rkind),device,allocatable :: scalarThroughfallSnow(:) ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarThroughfallRain(:) ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopySnowUnloading(:) ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyLiqDrainage(:) ! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: scalarCanopyMeltFreeze(:) ! get_ixVarType('scalarv')
  ! energy fluxes and for the snow and soil domains
  real(rkind),device,allocatable :: iLayerConductiveFlux_m(:,:) ! get_ixVarType('ifcToto')
  real(rkind),device,allocatable :: iLayerAdvectiveFlux_m(:,:) ! get_ixVarType('ifcToto')
  real(rkind),device,allocatable :: iLayerNrgFlux_m(:,:) ! get_ixVarType('ifcToto')
  real(rkind),device,allocatable :: mLayerNrgFlux_m(:,:) ! get_ixVarType('midToto')
  ! liquid water fluxes for the snow domain
  real(rkind),device,allocatable :: scalarSnowDrainage(:)! get_ixVarType('scalarv')
  real(rkind),device,allocatable :: iLayerLiqFluxSnow_m(:,:)! get_ixVarType('ifcSnow')
  real(rkind),device,allocatable :: mLayerLiqFluxSnow_m(:,:)! get_ixVarType('midSnow')
  ! liquid water fluxes for the soil domain
  real(rkind),device,allocatable ::scalarRainPlusMelt(:)!, get_ixVarType('scalarv'), 
  real(rkind),device,allocatable ::scalarMaxInfilRate(:)!, get_ixVarType('scalarv'), 
  real(rkind),device,allocatable ::scalarInfiltration(:)!, get_ixVarType('scalarv'), 
  real(rkind),device,allocatable ::scalarExfiltration(:)!, get_ixVarType('scalarv'), 
  real(rkind),device,allocatable ::scalarSurfaceRunoff(:)!, get_ixVarType('scalarv'), 
  real(rkind),device,allocatable ::mLayerSatHydCondMP_m(:,:)!, get_ixVarType('midSoil'), 
  real(rkind),device,allocatable ::mLayerSatHydCond_m(:,:)!, get_ixVarType('midSoil'), 
  real(rkind),device,allocatable ::iLayerSatHydCond_m(:,:)!, get_ixVarType('ifcSoil'), 
  real(rkind),device,allocatable ::mLayerHydCond_m(:,:)!, get_ixVarType('midSoil'), 
  real(rkind),device,allocatable ::iLayerLiqFluxSoil_m(:,:)!, get_ixVarType('ifcSoil'), 
  real(rkind),device,allocatable ::mLayerLiqFluxSoil_m(:,:)!, get_ixVarType('midSoil'), 
  real(rkind),device,allocatable ::mLayerBaseflow_m(:,:)!, get_ixVarType('midSoil'), 
  real(rkind),device,allocatable ::mLayerColumnInflow(:,:)!, get_ixVarType('midSoil'), 
  real(rkind),device,allocatable ::mLayerColumnOutflow_m(:,:)!, get_ixVarType('midSoil'), 
  real(rkind),device,allocatable ::scalarSoilBaseflow(:)!, get_ixVarType('scalarv'), 
  real(rkind),device,allocatable ::scalarSoilDrainage(:)!, get_ixVarType('scalarv'), 
  real(rkind),device,allocatable ::scalarAquiferRecharge(:)!, get_ixVarType('scalarv'), 
  real(rkind),device,allocatable ::scalarAquiferTranspire(:)!, get_ixVarType('scalarv'), 
  real(rkind),device,allocatable ::scalarAquiferBaseflow(:)!, get_ixVarType('scalarv'), 
  ! derived variables
  real(rkind),device,allocatable :: scalarTotalET(:)
  real(rkind),device,allocatable :: scalarTotalRunoff(:)
  real(rkind),device,allocatable :: scalarNetRadiation(:)


end type flux_data_device

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
  real(rkind),device,allocatable :: dFracLiqSnow_dTk_m(:,:)
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
  ! integer(i4b),device,allocatable :: nSnow(:)
  ! integer(i4b) :: nSoil 
  ! integer(i4b) :: nLayers !!!
  ! integer(i4b),device,allocatable :: layerType(:,:)
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
  ! ! ! number of state variables within different domains in the snow+soil system
  ! ! integer(i4b),device,allocatable :: nSnowSoilNrg_d(:)
  ! ! integer(i4b),device,allocatable :: nSnowOnlyNrg_d(:)
  ! ! integer(i4b),device,allocatable :: nSoilOnlyNrg_d(:)
  ! ! integer(i4b),device,allocatable :: nSnowSoilHyd_d(:)
  ! ! integer(i4b),device,allocatable :: nSnowOnlyHyd_d(:)
  ! ! integer(i4b),device,allocatable :: nSoilOnlyHyd_d(:)
  ! ! ! type of model state variables
  ! integer(i4b),device,allocatable :: ixControlVolume(:,:)
  ! integer(i4b),device,allocatable :: ixDomainType(:,:)
  ! integer(i4b),device,allocatable :: ixStateType(:,:)
  ! integer(i4b),device,allocatable :: ixHydType(:,:)
  ! ! ! type of model state variables (state subset)
  ! integer(i4b),device,allocatable :: ixDomainType_subset(:,:)
  ! integer(i4b),device,allocatable :: ixStateType_subset(:,:)
  ! ! ! mapping between state subset and the full state vector
  ! integer(i4b),device,allocatable :: ixMapFull2Subset(:,:)
  ! integer(i4b),device,allocatable :: ixMapSubset2Full(:,:)
  ! ! ! indices of model specific state variables
  ! integer(i4b),device,allocatable :: ixCasNrg(:)
  ! integer(i4b),device,allocatable :: ixVegNrg(:)
  ! integer(i4b),device,allocatable :: ixVegHyd(:)
  ! integer(i4b),device,allocatable :: ixTopNrg(:)
  ! integer(i4b),device,allocatable :: ixTopHyd(:)
  ! integer(i4b),device,allocatable :: ixAqWat(:)
  ! ! ! vectors of indices for specific state types
  ! integer(i4b),device,allocatable :: ixNrgOnly(:,:)
  ! integer(i4b),device,allocatable :: ixHydOnly(:,:)
  ! integer(i4b),device,allocatable :: ixMatOnly(:,:)
  ! integer(i4b),device,allocatable :: ixMassOnly(:,:)
  ! ! vectors of indices for specific state types within specific sub-domains
  ! integer(i4b),device,allocatable :: ixSnowSoilNrg(:,:)
  ! integer(i4b),device,allocatable :: ixSnowOnlyNrg(:,:)
  ! integer(i4b),device,allocatable :: ixSoilOnlyNrg(:,:)
  ! integer(i4b),device,allocatable :: ixSnowSoilHyd(:,:)
  ! integer(i4b),device,allocatable :: ixSnowOnlyHyd(:,:)
  ! integer(i4b),device,allocatable :: ixSoilOnlyHyd(:,:)
  ! ! ! vectors of indices for specfic state types within specific sub-domains
  integer(i4b),device,allocatable :: ixNrgCanair(:)
  ! integer(i4b),device,allocatable :: ixNrgCanopy(:)
  ! integer(i4b),device,allocatable :: ixHydCanopy(:)
  ! integer(i4b),device,allocatable :: ixNrgLayer(:,:)
  ! integer(i4b),device,allocatable :: ixHydLayer(:,:)
  integer(i4b),device,allocatable :: ixWatAquifer(:)

  ! ! ! vectors of indices for specific state types IN SPECIFIC SUB-DOMAINS
  ! ! ! indices within state vectors
  integer(i4b),device,allocatable :: ixAllState(:,:)
  ! integer(i4b),device,allocatable :: ixSoilState(:,:)
  integer(i4b),device,allocatable :: ixLayerState(:,:)

  integer(i4b),device,allocatable :: ixNrgLayer(:,:), ixHydLayer(:,:)
  integer(i4b) :: nSoil!, nLayers
  integer(i4b),device,allocatable :: nSnow(:)!,nSoil_d
  integer(i4b),device,allocatable :: nLayers_d(:)
  integer(i4b),device,allocatable :: ixCasNrg(:), ixVegNrg(:), ixVegHyd(:), ixAqWat(:)
  integer(i4b),device,allocatable :: ixTopHyd(:), ixTopNrg(:)
  integer(i4b),device,allocatable :: ixSnowSoilNrg(:,:)
  integer(i4b),device,allocatable :: ixSnowSoilHyd(:,:) 
  ! integer(i4b) :: nSnowSoilNrg, nSnowSoilHyd
  integer(i4b),device,allocatable :: ixHydType(:,:)
  integer(i4b),device,allocatable :: ixStateType(:,:)
  integer(i4b),device,allocatable :: ixStateType_subset(:,:)
  integer(i4b),device,allocatable :: ixHydCanopy(:)
  integer(i4b),device,allocatable :: ixNrgCanopy(:)
  ! integer(i4b) :: ixNrgCanair
  integer(i4b),device,allocatable :: ixMapSubset2Full(:,:), ixMapFull2Subset(:,:)
  integer(i4b),device,allocatable :: ixDomainType_subset(:,:)
  integer(i4b),device,allocatable :: ixControlVolume(:,:)
  integer(i4b),device,allocatable :: layerType(:,:)
  integer(i4b),device,allocatable :: ixSoilOnlyHyd(:,:)
  integer(i4b),device,allocatable :: ixSnowOnlyHyd(:,:)
  integer(i4b),device,allocatable :: ixSoilOnlyNrg(:,:)
  integer(i4b),device,allocatable :: ixSnowOnlyNrg(:,:)

  integer(i4b),device,allocatable :: ixSoilState(:,:)
  ! integer(i4b),device,allocatable :: ixLayerState(:)
  integer(i4b) :: numberFluxCalc
  integer(i4b),device,allocatable :: ixDomainType(:,:)
  integer(i4b),device,allocatable :: ixNrgOnly(:,:)

  ! Not used:
  ! ixMatricHead(:,:)
  ! ixVolFracWat(:,:)
  ! ixLayerActive(:,:)

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

real(rkind),device,allocatable :: data_step
real(rkind),device,allocatable :: refJulDay

real(rkind),device,allocatable :: h_lookup(:), t_lookup(:)
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
  real(rkind),device,allocatable :: rsmax
  real(rkind),device,allocatable :: topt
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
  real(rkind),device,allocatable :: hvt(:)
  real(rkind),device,allocatable :: hvb(:)
  real(rkind),device,allocatable :: rc(:)
  real(rkind),device,allocatable :: swemx
  real(rkind),device,allocatable :: albsat(:,:)
  real(rkind),device,allocatable :: albdry(:,:)
  real(rkind),device,allocatable :: alblak(:)
  integer(i4b),device,allocatable :: isbarren
integer(i4b),device,allocatable :: ISSNOW
integer(i4b),device,allocatable :: ISWATER
real(rkind),device,allocatable :: saim(:,:)
real(rkind),device,allocatable :: laim(:,:)
integer(i4b),device,allocatable :: dveg
real(rkind),device,allocatable :: tmin(:)
end type
type, public :: veg_param_tables
real(rkind),device,allocatable :: rgltbl(:), rstbl(:), hstbl(:)
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

end module device_data_types