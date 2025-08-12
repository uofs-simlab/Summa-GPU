module initialize_device
    use data_types
    use device_data_types
    USE var_lookup,only:iLookFLUX        ! lookup indices for flux data
    USE var_lookup,only:iLookDERIV       ! lookup indices for derivative data
    USE var_lookup,only:iLookFORCE       ! lookup indices for forcing data 
    USE var_lookup,only:iLookDIAG        ! lookup indices for diagnostic variable data
    USE var_lookup,only:iLookPROG,iLookPARAM, iLookINDEX, iLookBVAR, iLookTYPE, iLookATTR
    USE globalData,only: nBand,maxSnowLayers                 ! number of spectral bands
    USE globalData,only:integerMissing            ! missing integer

    contains

    subroutine allocate_device_param_data(mpar_data_device,mpar_data)
        type(mpar_data_device), intent(inout) :: mpar_data_device
        type(var_dlength),intent(in) :: mpar_data
        mpar_data_device%Fcapil = mpar_data%var(iLookPARAM%Fcapil)%dat(1)
        mpar_data_device%Louis79_bparam = mpar_data%var(iLookPARAM%Louis79_bparam)%dat(1)
        mpar_data_device%Louis79_cStar = mpar_data%var(iLookPARAM%Louis79_cStar)%dat(1)
        mpar_data_device%Mahrt87_eScale = mpar_data%var(iLookPARAM%Mahrt87_eScale)%dat(1)
        mpar_data_device%aquiferBaseflowExp = mpar_data%var(iLookPARAM%aquiferBaseflowExp)%dat(1)
        mpar_data_device%aquiferBaseflowRate = mpar_data%var(iLookPARAM%aquiferBaseflowRate)%dat(1)
        mpar_data_device%aquiferScaleFactor = mpar_data%var(iLookPARAM%aquiferScaleFactor)%dat(1)
        mpar_data_device%canopyDrainageCoeff = mpar_data%var(iLookPARAM%canopyDrainageCoeff)%dat(1)
        mpar_data_device%canopyWettingExp = mpar_data%var(iLookPARAM%canopyWettingExp)%dat(1)
        mpar_data_device%canopyWettingFactor = mpar_data%var(iLookPARAM%canopyWettingFactor)%dat(1)
        mpar_data_device%critAquiferTranspire = mpar_data%var(iLookPARAM%critAquiferTranspire)%dat(1)
        mpar_data_device%critRichNumber = mpar_data%var(iLookPARAM%critRichNumber)%dat(1)
        mpar_data_device%critSoilTranspire = mpar_data%var(iLookPARAM%critSoilTranspire)%dat(1)
        mpar_data_device%critSoilWilting = mpar_data%var(iLookPARAM%critSoilWilting)%dat(1)
        mpar_data_device%f_impede = mpar_data%var(iLookPARAM%f_impede)%dat(1)
        mpar_data_device%fixedThermalCond_snow = mpar_data%var(iLookPARAM%fixedThermalCond_snow)%dat(1)
        mpar_data_device%frac_clay = mpar_data%var(iLookPARAM%frac_clay)%dat
        mpar_data_device%frac_sand = mpar_data%var(iLookPARAM%frac_sand)%dat
        mpar_data_device%frac_silt = mpar_data%var(iLookPARAM%frac_silt)%dat
        mpar_data_device%heightCanopyBottom = mpar_data%var(iLookPARAM%heightCanopyBottom)%dat(1)
        mpar_data_device%heightCanopyTop = mpar_data%var(iLookPARAM%heightCanopyTop)%dat(1)
        mpar_data_device%kAnisotropic = mpar_data%var(iLookPARAM%kAnisotropic)%dat(1)
        mpar_data_device%k_snow = mpar_data%var(iLookPARAM%k_snow)%dat(1)
        mpar_data_device%leafDimension = mpar_data%var(iLookPARAM%leafDimension)%dat(1)
        mpar_data_device%leafExchangeCoeff = mpar_data%var(iLookPARAM%leafExchangeCoeff)%dat(1)
        mpar_data_device%lowerBoundHead = mpar_data%var(iLookPARAM%lowerBoundHead)%dat(1)
        mpar_data_device%lowerBoundTemp = mpar_data%var(iLookPARAM%lowerBoundTemp)%dat(1)
        mpar_data_device%lowerBoundTheta = mpar_data%var(iLookPARAM%lowerBoundTheta)%dat(1)
        mpar_data_device%maxMassVegetation = mpar_data%var(iLookPARAM%maxMassVegetation)%dat(1)
        mpar_data_device%minStomatalResistance = mpar_data%var(iLookPARAM%minStomatalResistance)%dat(1)
        mpar_data_device%mpExp = mpar_data%var(iLookPARAM%mpExp)%dat(1)
        mpar_data_device%mw_exp = mpar_data%var(iLookPARAM%mw_exp)%dat(1)
        mpar_data_device%plantWiltPsi = mpar_data%var(iLookPARAM%plantWiltPsi)%dat(1)
        mpar_data_device%qSurfScale = mpar_data%var(iLookPARAM%qSurfScale)%dat(1)
        mpar_data_device%rootingDepth = mpar_data%var(iLookPARAM%rootingDepth)%dat(1)
        mpar_data_device%snowfrz_scale = mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)
        mpar_data_device%soilIceCV = mpar_data%var(iLookPARAM%soilIceCV)%dat(1)
        mpar_data_device%soilIceScale = mpar_data%var(iLookPARAM%soilIceScale)%dat(1)
        mpar_data_device%soilStressParam = mpar_data%var(iLookPARAM%soilStressParam)%dat(1)
        mpar_data_device%soil_dens_intr = mpar_data%var(iLookPARAM%soil_dens_intr)%dat
        mpar_data_device%specificHeatVeg = mpar_data%var(iLookPARAM%specificHeatVeg)%dat(1)
        mpar_data_device%specificStorage = mpar_data%var(iLookPARAM%specificStorage)%dat(1)
        mpar_data_device%thCond_soil = mpar_data%var(iLookPARAM%thCond_soil)%dat
        mpar_data_device%theta_mp = mpar_data%var(iLookPARAM%theta_mp)%dat(1)
        mpar_data_device%theta_res = mpar_data%var(iLookPARAM%theta_res)%dat
        mpar_data_device%theta_sat = mpar_data%var(iLookPARAM%theta_sat)%dat
        mpar_data_device%throughfallScaleRain = mpar_data%var(iLookPARAM%throughfallScaleRain)%dat(1)
        mpar_data_device%upperBoundHead = mpar_data%var(iLookPARAM%upperBoundHead)%dat(1)
        mpar_data_device%upperBoundTemp = mpar_data%var(iLookPARAM%upperBoundTemp)%dat(1)
        mpar_data_device%upperBoundTheta = mpar_data%var(iLookPARAM%upperBoundTheta)%dat(1)
        mpar_data_device%vGn_alpha = mpar_data%var(iLookPARAM%vGn_alpha)%dat
        mpar_data_device%vGn_n       = mpar_data%var(iLookPARAM%vGn_n      )%dat
        mpar_data_device%wettingFrontSuction = mpar_data%var(iLookPARAM%wettingFrontSuction)%dat(1)
        mpar_data_device%windReductionParam = mpar_data%var(iLookPARAM%windReductionParam)%dat(1)
        mpar_data_device%z0Canopy = mpar_data%var(iLookPARAM%z0Canopy)%dat(1)
        mpar_data_device%z0Snow = mpar_data%var(iLookPARAM%z0Snow)%dat(1)
        mpar_data_device%z0Soil = mpar_data%var(iLookPARAM%z0Soil)%dat(1)
        mpar_data_device%zScale_TOPMODEL = mpar_data%var(iLookPARAM%zScale_TOPMODEL)%dat(1)
        mpar_data_device%zpdFraction = mpar_data%var(iLookPARAM%zpdFraction)%dat(1)
        mpar_data_device%Kc25 = mpar_data%var(iLookPARAM%Kc25)%dat(1)
        mpar_data_device%Ko25 = mpar_data%var(iLookPARAM%Ko25)%dat(1)
        mpar_data_device%Kc_qFac = mpar_data%var(iLookPARAM%Kc_qFac)%dat(1)
        mpar_data_device%Ko_qFac = mpar_data%var(iLookPARAM%Ko_qFac)%dat(1)
        mpar_data_device%kc_Ha = mpar_data%var(iLookPARAM%kc_Ha)%dat(1)
        mpar_data_device%ko_Ha = mpar_data%var(iLookPARAM%ko_Ha)%dat(1)
        mpar_data_device%vcmax25_canopyTop = mpar_data%var(iLookPARAM%vcmax25_canopyTop)%dat(1)
        mpar_data_device%vcmax_qFac = mpar_data%var(iLookPARAM%vcmax_qFac)%dat(1)
        mpar_data_device%vcmax_Ha = mpar_data%var(iLookPARAM%vcmax_Ha)%dat(1)
        mpar_data_device%vcmax_Hd = mpar_data%var(iLookPARAM%vcmax_Hd)%dat(1)
        mpar_data_device%vcmax_Sv = mpar_data%var(iLookPARAM%vcmax_Sv)%dat(1)
        mpar_data_device%vcmax_Kn = mpar_data%var(iLookPARAM%vcmax_Kn)%dat(1)
        mpar_data_device%jmax25_scale = mpar_data%var(iLookPARAM%jmax25_scale)%dat(1)
        mpar_data_device%jmax_Ha = mpar_data%var(iLookPARAM%jmax_Ha)%dat(1)
        mpar_data_device%jmax_Hd = mpar_data%var(iLookPARAM%jmax_Hd)%dat(1)
        mpar_data_device%jmax_Sv = mpar_data%var(iLookPARAM%jmax_Sv)%dat(1)
        mpar_data_device%fractionJ = mpar_data%var(iLookPARAM%fractionJ)%dat(1)
        mpar_data_device%quantamYield = mpar_data%var(iLookPARAM%quantamYield)%dat(1)
        mpar_data_device%vpScaleFactor = mpar_data%var(iLookPARAM%vpScaleFactor)%dat(1)
        mpar_data_device%cond2photo_slope = mpar_data%var(iLookPARAM%cond2photo_slope)%dat(1)
        mpar_data_device%minStomatalConductance = mpar_data%var(iLookPARAM%minStomatalConductance)%dat(1)

        mpar_data_device%idaMaxOrder = mpar_data%var(iLookPARAM%idaMaxOrder)%dat(1)
        mpar_data_device%idaMaxInternalSteps = mpar_data%var(iLookPARAM%idaMaxInternalSteps)%dat(1)
        mpar_data_device%idaInitStepSize = mpar_data%var(iLookPARAM%idaInitStepSize)%dat(1)
        mpar_data_device%idaMinStepSize = mpar_data%var(iLookPARAM%idaMinStepSize)%dat(1)
        mpar_data_device%idaMaxErrTestFail = mpar_data%var(iLookPARAM%idaMaxErrTestFail)%dat(1)
        mpar_data_device%fieldCapacity = mpar_data%var(iLookPARAM%fieldCapacity)%dat(1)
        mpar_data_device%absTolTempCas = mpar_data%var(iLookPARAM%absTolTempCas)%dat(1)
        mpar_data_device%relTolTempCas = mpar_data%var(iLookPARAM%relTolTempCas)%dat(1)
        mpar_data_device%absTolTempVeg = mpar_data%var(iLookPARAM%absTolTempVeg)%dat(1)
        mpar_data_device%relTolTempVeg = mpar_data%var(iLookPARAM%relTolTempVeg)%dat(1)
        mpar_data_device%absTolWatVeg = mpar_data%var(iLookPARAM%absTolWatVeg)%dat(1)
        mpar_data_device%relTolWatVeg = mpar_data%var(iLookPARAM%relTolWatVeg)%dat(1)
        mpar_data_device%absTolTempSoilSnow = mpar_data%var(iLookPARAM%absTolTempSoilSnow)%dat(1)
        mpar_data_device%relTolTempSoilSnow = mpar_data%var(iLookPARAM%relTolTempSoilSnow)%dat(1)
        mpar_data_device%absTolWatSnow = mpar_data%var(iLookPARAM%absTolWatSnow)%dat(1)
        mpar_data_device%relTolWatSnow = mpar_data%var(iLookPARAM%relTolWatSnow)%dat(1)
        mpar_data_device%absTolMatric = mpar_data%var(iLookPARAM%absTolMatric)%dat(1)
        mpar_data_device%relTolMatric = mpar_data%var(iLookPARAM%relTolMatric)%dat(1)
        mpar_data_device%absTolAquifr = mpar_data%var(iLookPARAM%absTolAquifr)%dat(1)
        mpar_data_device%relTolAquifr = mpar_data%var(iLookPARAM%relTolAquifr)%dat(1)
        mpar_data_device%maxstep = mpar_data%var(iLookPARAM%maxstep)%dat(1)
        mpar_data_device%be_steps = mpar_data%var(iLookPARAM%be_steps)%dat(1)
        mpar_data_device%densScalGrowth = mpar_data%var(iLookPARAM%densScalGrowth)%dat(1)
        mpar_data_device%tempScalGrowth = mpar_data%var(iLookPARAM%tempScalGrowth)%dat(1)
        mpar_data_device%densScalOvrbdn = mpar_data%var(iLookPARAM%densScalOvrbdn)%dat(1)
        mpar_data_device%tempScalOvrbdn = mpar_data%var(iLookPARAM%tempScalOvrbdn)%dat(1)
        mpar_data_device%baseViscosity = mpar_data%var(iLookPARAM%baseViscosity)%dat(1)
        mpar_data_device%grainGrowthRate = mpar_data%var(iLookPARAM%grainGrowthRate)%dat(1)
        mpar_data_device%zmaxLayer1_lower = mpar_data%var(iLookPARAM%zmaxLayer1_lower)%dat(1)
        mpar_data_device%zmaxLayer2_lower = mpar_data%var(iLookPARAM%zmaxLayer2_lower)%dat(1)
        mpar_data_device%zmaxLayer3_lower = mpar_data%var(iLookPARAM%zmaxLayer3_lower)%dat(1)
        mpar_data_device%zmaxLayer4_lower = mpar_data%var(iLookPARAM%zmaxLayer4_lower)%dat(1)
        mpar_data_device%zmaxLayer1_upper = mpar_data%var(iLookPARAM%zmaxLayer1_upper)%dat(1)
        mpar_data_device%zmaxLayer2_upper = mpar_data%var(iLookPARAM%zmaxLayer2_upper)%dat(1)
        mpar_data_device%zmaxLayer3_upper = mpar_data%var(iLookPARAM%zmaxLayer3_upper)%dat(1)
        mpar_data_device%zmaxLayer4_upper = mpar_data%var(iLookPARAM%zmaxLayer4_upper)%dat(1)
        mpar_data_device%zmax = mpar_data%var(iLookPARAM%zmax)%dat(1)
        mpar_data_device%Frad_vis = mpar_data%var(iLookPARAM%Frad_vis)%dat(1)
        mpar_data_device%albedoMaxNearIR = mpar_data%var(iLookPARAM%albedoMaxNearIR)%dat(1)
        mpar_data_device%albedoMaxVisible = mpar_data%var(iLookPARAM%albedoMaxVisible)%dat(1)
        mpar_data_device%albedoMax = mpar_data%var(iLookPARAM%albedoMax)%dat(1)
        mpar_data_device%zmin = mpar_data%var(iLookPARAM%zmin)%dat(1)
        mpar_data_device%zminLayer1 = mpar_data%var(iLookPARAM%zminLayer1)%dat(1)
        mpar_data_device%zminLayer2 = mpar_data%var(iLookPARAM%zminLayer2)%dat(1)
        mpar_data_device%zminLayer3 = mpar_data%var(iLookPARAM%zminLayer3)%dat(1)
        mpar_data_device%zminLayer4 = mpar_data%var(iLookPARAM%zminLayer4)%dat(1)
        mpar_data_device%zminLayer5 = mpar_data%var(iLookPARAM%zminLayer5)%dat(1)
        mpar_data_device%refInterceptCapSnow = mpar_data%var(iLookPARAM%refInterceptCapSnow)%dat(1)
                mpar_data_device%refInterceptCapRain = mpar_data%var(iLookPARAM%refInterceptCapRain)%dat(1)

mpar_data_device%ratioDrip2Unloading = mpar_data%var(iLookPARAM%ratioDrip2Unloading)%dat(1)
mpar_data_device%snowUnloadingCoeff = mpar_data%var(iLookPARAM%snowUnloadingCoeff)%dat(1)
mpar_data_device%minTempUnloading = mpar_data%var(iLookPARAM%minTempUnloading)%dat(1)
mpar_data_device%minWindUnloading = mpar_data%var(iLookPARAM%minWindUnloading)%dat(1)
mpar_data_device%rateTempUnloading = mpar_data%var(iLookPARAM%rateTempUnloading)%dat(1)
mpar_data_device%rateWindUnloading = mpar_data%var(iLookPARAM%rateWindUnloading)%dat(1)
mpar_data_device%albedoMinWinter = mpar_data%var(iLookPARAM%albedoMinWinter)%dat(1)
mpar_data_device%albedoMinSpring = mpar_data%var(iLookPARAM%albedoMinSpring)%dat(1)
mpar_data_device%albedoMinVisible = mpar_data%var(iLookPARAM%albedoMinVisible)%dat(1)
mpar_data_device%albedoMinNearIR = mpar_data%var(iLookPARAM%albedoMinNearIR)%dat(1)
mpar_data_device%albedoDecayRate = mpar_data%var(iLookPARAM%albedoDecayRate)%dat(1)
mpar_data_device%albedoSootLoad = mpar_data%var(iLookPARAM%albedoSootLoad)%dat(1)
mpar_data_device%albedoRefresh = mpar_data%var(iLookPARAM%albedoRefresh)%dat(1)
mpar_data_device%Frad_direct = mpar_data%var(iLookPARAM%Frad_direct)%dat(1)
mpar_data_device%minstep = mpar_data%var(iLookPARAM%minstep)%dat(1)
mpar_data_device%directScale = mpar_data%var(iLookPARAM%directScale)%dat(1)
mpar_data_device%newSnowDenMin = mpar_data%var(iLookPARAM%newSnowDenMin)%dat(1)
mpar_data_device%newSnowDenMult = mpar_data%var(iLookPARAM%newSnowDenMult)%dat(1)
mpar_data_device%newSnowDenScal = mpar_data%var(iLookPARAM%newSnowDenScal)%dat(1)
mpar_data_device%constSnowDen = mpar_data%var(iLookPARAM%constSnowDen)%dat(1)
mpar_data_device%newSnowDenAdd = mpar_data%var(iLookPARAM%newSnowDenAdd)%dat(1)
mpar_data_device%newSnowDenMultTemp = mpar_data%var(iLookPARAM%newSnowDenMultTemp)%dat(1)
mpar_data_device%newSnowDenMultWind = mpar_data%var(iLookPARAM%newSnowDenMultWind)%dat(1)
mpar_data_device%newSnowDenMultAnd = mpar_data%var(iLookPARAM%newSnowDenMultAnd)%dat(1)
mpar_data_device%newSnowDenBase = mpar_data%var(iLookPARAM%newSnowDenBase)%dat(1)
mpar_data_device%tempCritRain = mpar_data%var(iLookPARAM%tempCritRain)%dat(1)
mpar_data_device%tempRangeTimestep = mpar_data%var(iLookPARAM%tempRangeTimestep)%dat(1)
mpar_data_device%frozenPrecipMultip = mpar_data%var(iLookPARAM%frozenPrecipMultip)%dat(1)
mpar_data_device%minwind = mpar_data%var(iLookPARAM%minwind)%dat(1)
mpar_data_device%winterSAI = mpar_data%var(iLookPARAM%winterSAI)%dat(1)
mpar_data_device%summerLAI = mpar_data%var(iLookPARAM%summerLAI)%dat(1)
mpar_data_device%k_soil = mpar_data%var(iLookPARAM%k_soil)%dat
mpar_data_device%k_macropore = mpar_data%var(iLookPARAM%k_macropore)%dat
mpar_data_device%compactedDepth = mpar_data%var(iLookPARAM%compactedDepth)%dat(1)
mpar_data_device%rootScaleFactor1 = mpar_data%var(iLookPARAM%rootScaleFactor1)%dat(1)
mpar_data_device%rootScaleFactor2 = mpar_data%var(iLookPARAM%rootScaleFactor2)%dat(1)
mpar_data_device%rootDistExp = mpar_data%var(iLookPARAM%rootDistExp)%dat(1)
      end subroutine allocate_device_param_data

      subroutine deallocate_device_param_data(mpar_data_device)
        type(mpar_data_device), intent(inout) :: mpar_data_device
        deallocate(mpar_data_device%Fcapil)
        deallocate(mpar_data_device%Louis79_bparam)
        deallocate(mpar_data_device%Louis79_cStar)
        deallocate(mpar_data_device%Mahrt87_eScale)
        deallocate(mpar_data_device%aquiferBaseflowExp)
        deallocate(mpar_data_device%aquiferBaseflowRate)
        deallocate(mpar_data_device%aquiferScaleFactor)
        deallocate(mpar_data_device%canopyDrainageCoeff)
        deallocate(mpar_data_device%canopyWettingExp)
        deallocate(mpar_data_device%canopyWettingFactor)
        deallocate(mpar_data_device%critAquiferTranspire)
        deallocate(mpar_data_device%critRichNumber)
        deallocate(mpar_data_device%critSoilTranspire)
        deallocate(mpar_data_device%critSoilWilting)
        deallocate(mpar_data_device%f_impede)
        deallocate(mpar_data_device%fixedThermalCond_snow)
        deallocate(mpar_data_device%frac_clay)
        deallocate(mpar_data_device%frac_sand)
        deallocate(mpar_data_device%frac_silt)
        deallocate(mpar_data_device%heightCanopyBottom)
        deallocate(mpar_data_device%heightCanopyTop)
        deallocate(mpar_data_device%kAnisotropic)
        deallocate(mpar_data_device%k_snow)
        deallocate(mpar_data_device%leafDimension)
        deallocate(mpar_data_device%leafExchangeCoeff)
        deallocate(mpar_data_device%lowerBoundHead)
        deallocate(mpar_data_device%lowerBoundTemp)
        deallocate(mpar_data_device%lowerBoundTheta)
        deallocate(mpar_data_device%maxMassVegetation)
        deallocate(mpar_data_device%minStomatalResistance)
        deallocate(mpar_data_device%mpExp)
        deallocate(mpar_data_device%mw_exp)
        deallocate(mpar_data_device%plantWiltPsi)
        deallocate(mpar_data_device%qSurfScale)
        deallocate(mpar_data_device%rootingDepth)
        deallocate(mpar_data_device%snowfrz_scale)
        deallocate(mpar_data_device%soilIceCV)
        deallocate(mpar_data_device%soilIceScale)
        deallocate(mpar_data_device%soilStressParam)
        deallocate(mpar_data_device%soil_dens_intr)
        deallocate(mpar_data_device%specificHeatVeg)
        deallocate(mpar_data_device%specificStorage)
        deallocate(mpar_data_device%thCond_soil)
        deallocate(mpar_data_device%theta_mp)
        deallocate(mpar_data_device%theta_res)
        deallocate(mpar_data_device%theta_sat)
        deallocate(mpar_data_device%throughfallScaleRain)
        deallocate(mpar_data_device%upperBoundHead)
        deallocate(mpar_data_device%upperBoundTemp)
        deallocate(mpar_data_device%upperBoundTheta)
        deallocate(mpar_data_device%vGn_alpha)
        deallocate(mpar_data_device%vGn_n)
        deallocate(mpar_data_device%wettingFrontSuction)
        deallocate(mpar_data_device%windReductionParam)
        deallocate(mpar_data_device%z0Canopy)
        deallocate(mpar_data_device%z0Snow)
        deallocate(mpar_data_device%z0Soil)
        deallocate(mpar_data_device%zScale_TOPMODEL)
        deallocate(mpar_data_device%zpdFraction)
        deallocate(mpar_data_device%Kc25)
        deallocate(mpar_data_device%Ko25)
        deallocate(mpar_data_device%Kc_qFac)
        deallocate(mpar_data_device%Ko_qFac)
        deallocate(mpar_data_device%kc_Ha)
        deallocate(mpar_data_device%ko_Ha)
        deallocate(mpar_data_device%vcmax25_canopyTop)
        deallocate(mpar_data_device%vcmax_qFac)
        deallocate(mpar_data_device%vcmax_Ha)
        deallocate(mpar_data_device%vcmax_Hd)
        deallocate(mpar_data_device%vcmax_Sv)
        deallocate(mpar_data_device%vcmax_Kn)
        deallocate(mpar_data_device%jmax25_scale)
        deallocate(mpar_data_device%jmax_Ha)
        deallocate(mpar_data_device%jmax_Hd)
        deallocate(mpar_data_device%jmax_Sv)
        deallocate(mpar_data_device%fractionJ)
        deallocate(mpar_data_device%quantamYield)
        deallocate(mpar_data_device%vpScaleFactor)
        deallocate(mpar_data_device%cond2photo_slope)
        deallocate(mpar_data_device%minStomatalConductance)
      
      end subroutine deallocate_device_param_data


  subroutine allocate_device_indx_data(indx_data_device, indx_data, nGRU)
    implicit none
    type(indx_data_device), intent(inout) :: indx_data_device
    type(gru_hru_intVec),intent(in) :: indx_data              ! indices defining model states and layers
    integer(i4b), intent(in) :: nGRU
    integer(i4b) :: iGRU
    integer(i4b) :: nSnow,nSoil,nLayers
        allocate(indx_data_device%nSnow(nGRU))
allocate(indx_data_device%nLayers_d(nGRU))

    do iGRU=1,nGRU
      indx_data_device%nLayers_d(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nLayers)%dat(1)
      indx_data_device%nSnow(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
    end do
    nSnow = maxval(indx_data_device%nSnow)
    nSoil = indx_data%gru(1)%hru(1)%var(iLookINDEX%nSoil)%dat(1)
    nLayers = maxval(indx_data_device%nLayers_d)
    indx_data_device%nSoil = nSoil
    allocate(indx_data_device%layerType(nSoil+maxSnowLayers,nGRU))
    allocate(indx_data_device%nCasNrg(nGRU))
    allocate(indx_data_device%nVegNrg(nGRU))
    allocate(indx_data_device%nVegMass(nGRU))
    allocate(indx_data_device%nVegState(nGRU))
    allocate(indx_data_device%nNrgState(nGRU))
    allocate(indx_data_device%nWatState(nGRU))
    allocate(indx_data_device%nMatState(nGRU))
    allocate(indx_data_device%nMassState(nGRU))
    allocate(indx_data_device%nState(nGRU))
    ! allocate(indx_data_device%nSnowSoilNrg_d(nGRU))
    ! indx_data_device%nSnowSoilNrg = indx_data%var(iLookINDEX%nSnowSoilNrg)%dat(1)
    ! allocate(indx_data_device%nSnowOnlyNrg_d(nGRU))
    ! indx_data_device%nSnowOnlyNrg = indx_data%var(iLookINDEX%nSnowOnlyNrg)%dat(1)
    ! allocate(indx_data_device%nSoilOnlyNrg_d(nGRU))
    ! indx_data_device%nSoilOnlyNrg = indx_data%var(iLookINDEX%nSoilOnlyNrg)%dat(1)
    ! allocate(indx_data_device%nSnowSoilHyd_d(nGRU))
    ! indx_data_device%nSnowSoilHyd = indx_data%var(iLookINDEX%nSnowSoilHyd)%dat(1)
    ! allocate(indx_data_device%nSnowOnlyHyd_d(nGRU))
    ! indx_data_device%nSnowOnlyHyd = indx_data%var(iLookINDEX%nSnowOnlyHyd)%dat(1)
    ! allocate(indx_data_device%nSoilOnlyHyd_d(nGRU))
    ! indx_data_device%nSoilOnlyHyd = indx_data%var(iLookINDEX%nSoilOnlyHyd)%dat(1)
    ! print*, nSoil, nSnow, maxSnowLayers, nLayers
    if (size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixControlVolume)%dat).ne.0) allocate(indx_data_device%ixControlVolume(size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixControlVolume)%dat),nGRU))
    if (size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixDomainType)%dat).ne.0) allocate(indx_data_device%ixDomainType(size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixDomainType)%dat),nGRU))
    allocate(indx_data_device%ixStateType(MAX(1,size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixStateType)%dat)),nGRU))
    allocate(indx_data_device%ixHydType(nSoil+maxSnowLayers,nGRU))
    !if (size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixDomainType_subset)%dat).ne.0) allocate(indx_data_device%ixDomainType_subset(size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixDomainType_subset)%dat),nGRU))
    ! if (size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixStateType_subset)%dat).ne.0) allocate(indx_data_device%ixStateType_subset(size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixStateType_subset)%dat),nGRU))
    ! allocate(indx_data_device%ixMapFull2Subset(MAX(1,size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixMapFull2Subset)%dat)),nGRU))
    ! print*, size(indx_data%var(iLookINDEX%ixMapFull2Subset)%dat)
    ! allocate(indx_data_device%ixMapSubset2Full(MAX(1,size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixMapSubset2Full)%dat)),nGRU))
    ! print*, size(indx_data%var(iLookINDEX%ixMapSubset2Full)%dat)
    allocate(indx_data_device%ixCasNrg(nGRU))
    allocate(indx_data_device%ixVegNrg(nGRU))
    allocate(indx_data_device%ixVegHyd(nGRU))
    allocate(indx_data_device%ixTopNrg(nGRU))
    allocate(indx_data_device%ixTopHyd(nGRU))
    allocate(indx_data_device%ixAqWat(nGRU))
    if (size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixNrgOnly)%dat).ne.0) allocate(indx_data_device%ixNrgOnly(size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixNrgOnly)%dat),nGRU))
    ! if (size(indx_data%var(iLookINDEX%ixHydOnly)%dat).ne.0) allocate(indx_data_device%ixHydOnly(size(indx_data%var(iLookINDEX%ixHydOnly)%dat),nGRU))
    ! if (size(indx_data%var(iLookINDEX%ixMatOnly)%dat).ne.0) allocate(indx_data_device%ixMatOnly(size(indx_data%var(iLookINDEX%ixMatOnly)%dat),nGRU))
    ! if (size(indx_data%var(iLookINDEX%ixMassOnly)%dat).ne.0) allocate(indx_data_device%ixMassOnly(size(indx_data%var(iLookINDEX%ixMassOnly)%dat),nGRU))
    allocate(indx_data_device%ixSnowSoilNrg(nSoil+maxSnowLayers,nGRU))
    allocate(indx_data_device%ixSnowOnlyNrg(maxSnowLayers,nGRU))
    allocate(indx_data_device%ixSoilOnlyNrg(nSoil,nGRU))
    allocate(indx_data_device%ixSnowSoilHyd(nSoil+maxSnowLayers,nGRU))
    allocate(indx_data_device%ixSnowOnlyHyd(maxSnowLayers,nGRU))
    allocate(indx_data_device%ixSoilOnlyHyd(nSoil,nGRU))
    allocate(indx_data_device%ixNrgCanair(nGRU))
    allocate(indx_data_device%ixNrgCanopy(nGRU))
    allocate(indx_data_device%ixHydCanopy(nGRU))
    allocate(indx_data_device%ixNrgLayer(nSoil+maxSnowLayers,nGRU))
    allocate(indx_data_device%ixHydLayer(nSoil+maxSnowLayers,nGRU))
    allocate(indx_data_device%ixWatAquifer(nGRU))
    ! if (size(indx_data%var(iLookINDEX%ixVolFracWat)%dat).ne.0) allocate(indx_data_device%ixVolFracWat(size(indx_data%var(iLookINDEX%ixVolFracWat)%dat),nGRU))
    ! if (size(indx_data%var(iLookINDEX%ixMatricHead)%dat).ne.0) allocate(indx_data_device%ixMatricHead(size(indx_data%var(iLookINDEX%ixMatricHead)%dat),nGRU))
    if (size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixAllState)%dat).ne.0) allocate(indx_data_device%ixAllState(size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixAllState)%dat),nGRU))
    if (size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixSoilState)%dat).ne.0) allocate(indx_data_device%ixSoilState(size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixSoilState)%dat),nGRU))
    if (size(indx_data%gru(1)%hru(1)%var(iLookINDEX%ixLayerState)%dat).ne.0) allocate(indx_data_device%ixLayerState(nSoil+maxSnowLayers,nGRU))
    ! if (size(indx_data%var(iLookINDEX%ixLayerActive)%dat).ne.0) allocate(indx_data_device%ixLayerActive(nSoil+maxSnowLayers,nGRU))
    indx_data_device%numberFluxCalc = indx_data%gru(1)%hru(1)%var(iLookINDEX%numberFluxCalc)%dat(1)
    ! indx_data_device%numberStateSplit = indx_data%var(iLookINDEX%numberStateSplit)%dat(1)
    ! indx_data_device%numberDomainSplitNrg = indx_data%var(iLookINDEX%numberDomainSplitNrg)%dat(1)
    ! indx_data_device%numberDomainSplitMass = indx_data%var(iLookINDEX%numberDomainSplitMass)%dat(1)
    ! indx_data_device%numberScalarSolutions = indx_data%var(iLookINDEX%numberScalarSolutions)%dat(1)

    do iGRU=1,nGRU
        indx_data_device%nCasNrg(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nCasNrg)%dat(1)
    indx_data_device%nVegNrg(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nVegNrg)%dat(1)
    indx_data_device%nVegMass(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nVegMass)%dat(1)
    indx_data_device%nVegState(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nVegState)%dat(1)
    indx_data_device%nNrgState(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nNrgState)%dat(1)
    indx_data_device%nWatState(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nWatState)%dat(1)
    indx_data_device%nMatState(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nMatState)%dat(1)
    indx_data_device%nMassState(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nMassState)%dat(1)
    indx_data_device%nState(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nState)%dat(1)
    indx_data_device%ixCasNrg(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixCasNrg)%dat(1)
    indx_data_device%ixVegNrg(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixVegNrg)%dat(1)
    indx_data_device%ixVegHyd(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixVegHyd)%dat(1)
    indx_data_device%ixTopNrg(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixTopNrg)%dat(1)
    indx_data_device%ixTopHyd(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixTopHyd)%dat(1)
    indx_data_device%ixAqWat(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixAqWat)%dat(1)
    indx_data_device%ixNrgCanair(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgCanair)%dat(1)
    indx_data_device%ixNrgCanopy(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgCanopy)%dat(1)
    indx_data_device%ixHydCanopy(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydCanopy)%dat(1)
    indx_data_device%ixWatAquifer(iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixWatAquifer)%dat(1)
    end do

    do iGRU=1,nGRU
      if (allocated(indx_data_device%layerType)) indx_data_device%layerType(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%layerType)%dat
      indx_data_device%layerType(nLayers+1:,iGRU) = integerMissing
      if (allocated(indx_data_device%ixControlVolume)) indx_data_device%ixControlVolume(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixControlVolume)%dat
      if (allocated(indx_data_device%ixDomainType)) indx_data_device%ixDomainType(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixDomainType)%dat
      if (allocated(indx_data_device%ixStateType)) indx_data_device%ixStateType(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixStateType)%dat
      if (allocated(indx_data_device%ixHydType)) indx_data_device%ixHydType(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydType)%dat
      !if (allocated(indx_data_device%ixDomainType_subset)) indx_data_device%ixDomainType_subset(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixDomainType_subset)%dat
      !if (allocated(indx_data_device%ixStateType_subset)) indx_data_device%ixStateType_subset(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixStateType_subset)%dat
      ! if (allocated(indx_data_device%ixMapFull2Subset)) indx_data_device%ixMapFull2Subset(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixMapFull2Subset)%dat
      ! if (allocated(indx_data_device%ixMapSubset2Full)) indx_data_device%ixMapSubset2Full(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixMapSubset2Full)%dat
      if (allocated(indx_data_device%ixNrgOnly)) indx_data_device%ixNrgOnly(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgOnly)%dat
      ! if (allocated(indx_data_device%ixHydOnly)) indx_data_device%ixHydOnly(:,iGRU) = indx_data%var(iLookINDEX%ixHydOnly)%dat
      ! if (allocated(indx_data_device%ixMatOnly)) indx_data_device%ixMatOnly(:,iGRU) = indx_data%var(iLookINDEX%ixMatOnly)%dat
      ! if (allocated(indx_data_device%ixMassOnly)) indx_data_device%ixMassOnly(:,iGRU) = indx_data%var(iLookINDEX%ixMassOnly)%dat
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilNrg)%dat).ne.0) then
        indx_data_device%ixSnowSoilNrg(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilNrg)%dat
        indx_data_device%ixSnowSOilNrg(nLayers+1:,iGRU) = integerMissing
      else
        indx_data_device%ixSnowSoilNrg(:,iGRU) = integerMissing
      endif
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyNrg)%dat).ne.0) then
        indx_data_device%ixSnowOnlyNrg(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyNrg)%dat
        indx_data_device%ixSnowOnlyNrg(nSnow+1:,iGRU) = integerMissing
      else
        indx_data_device%ixSnowOnlyNrg(:,iGRU) = integerMissing
      endif
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilOnlyNrg)%dat).ne.0) then
        indx_data_device%ixSoilOnlyNrg(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilOnlyNrg)%dat
      else
        indx_data_device%ixSoilOnlyNrg(:,iGRU) = integerMissing
      endif
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilHyd)%dat).ne.0) then
        indx_data_device%ixSnowSoilHyd(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilHyd)%dat
        indx_data_device%ixSnowSoilHyd(nLayers+1:,iGRU) = integerMissing
      else
        indx_data_device%ixSnowSoilHyd(:,iGRU) = integerMissing
      endif
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyHyd)%dat).ne.0) then
        indx_data_device%ixSnowOnlyHyd(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyHyd)%dat
        indx_data_device%ixSnowOnlyHyd(nSnow+1:,iGRU) = integerMissing
      else
        indx_data_device%ixSnowOnlyHyd(:,iGRU) = integerMissing
      endif
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilOnlyHyd)%dat).ne.0) then
        indx_data_device%ixSoilOnlyHyd(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilOnlyHyd)%dat
      else
        indx_data_device%ixSoilOnlyHyd(:,iGRU) = integerMissing
      endif
      if (allocated(indx_data_device%ixNrgLayer)) indx_data_device%ixNrgLayer(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgLayer)%dat
      indx_data_device%ixNrgLayer(nLayers+1:,iGRU) = integerMissing
      if (allocated(indx_data_device%ixHydLayer)) indx_data_device%ixHydLayer(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydLayer)%dat
      indx_data_device%ixHydLayer(nLayers+1:,iGRU) = integerMissing
      ! if (size(indx_data%var(iLookINDEX%ixVolFracWat)%dat).ne.0) indx_data_device%ixVolFracWat(:,iGRU) = indx_data%var(iLookINDEX%ixVolFracWat)%dat
      ! if (size(indx_data%var(iLookINDEX%ixMatricHead)%dat).ne.0) indx_data_device%ixMatricHead(:,iGRU) = indx_data%var(iLookINDEX%ixMatricHead)%dat
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixAllState)%dat).ne.0) indx_data_device%ixAllState(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixAllState)%dat
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilState)%dat).ne.0) indx_data_device%ixSoilState(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilState)%dat
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixLayerState)%dat).ne.0) indx_data_device%ixLayerState(:,iGRU) = indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixLayerState)%dat
      indx_data_device%ixLayerState(nLayers+1:,iGRU) = integerMissing
      ! if (size(indx_data%var(iLookINDEX%ixLayerActive)%dat).ne.0) indx_data_device%ixLayerActive(:,iGRU) = indx_data%var(iLookINDEX%ixLayerActive)%dat


    end do
  
  end subroutine allocate_device_indx_data

      subroutine finalize_device_indx_data(indx_data_device, indx_data,nGRU)
        implicit none
        type(indx_data_device), intent(in) :: indx_data_device
        type(gru_hru_intVec),intent(inout) :: indx_data              ! indices defining model states and layers
        integer(i4b) :: nLayers
        ! integer(i4b) :: nSnow, nSoil
        integer(i4b) :: nGRU
        integer(i4b) :: iGRU

        nLayers = maxval(indx_data_device%nLayers_d)
    !     nSnow = indx_data_device%nSnow(1)
    !     nSoil = indx_data_device%nSoil
    ! indx_data%var(iLookINDEX%nSnow)%dat(1) = indx_data_device%nSnow(1)
    ! indx_data%var(iLookINDEX%nSoil)%dat(1) = indx_data_device%nSoil
    ! indx_data%var(iLookINDEX%nLayers)%dat(1) = indx_data_device%nLayers
        do iGRU=1,nGRU
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nCasNrg)%dat(1) = indx_data_device%nCasNrg(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nVegNrg)%dat(1) = indx_data_device%nVegNrg(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nVegMass)%dat(1) = indx_data_device%nVegMass(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nVegState)%dat(1) = indx_data_device%nVegState(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nNrgState)%dat(1) = indx_data_device%nNrgState(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nWatState)%dat(1) = indx_data_device%nWatState(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nMatState)%dat(1) = indx_data_device%nMatState(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nMassState)%dat(1) = indx_data_device%nMassState(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nState)%dat(1) = indx_data_device%nState(iGRU)
    ! indx_data%var(iLookINDEX%nSnowSoilNrg)%dat(1) = indx_data_device%nSnowSoilNrg_d(1)
    ! indx_data%var(iLookINDEX%nSnowOnlyNrg)%dat(1) = indx_data_device%nSnowOnlyNrg_d(1)
    ! indx_data%var(iLookINDEX%nSoilOnlyNrg)%dat(1) = indx_data_device%nSoilOnlyNrg_d(1)
    ! indx_data%var(iLookINDEX%nSnowSoilHyd)%dat(1) = indx_data_device%nSnowSoilHyd_d(1)
    ! indx_data%var(iLookINDEX%nSnowOnlyHyd)%dat(1) = indx_data_device%nSnowOnlyHyd_d(1)
    ! indx_data%var(iLookINDEX%nSoilOnlyHyd)%dat(1) = indx_data_device%nSoilOnlyHyd_d(1)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixCasNrg)%dat(1) = indx_data_device%ixCasNrg(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixVegNrg)%dat(1) = indx_data_device%ixVegNrg(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixVegHyd)%dat(1) = indx_data_device%ixVegHyd(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixTopNrg)%dat(1) = indx_data_device%ixTopNrg(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixTopHyd)%dat(1) = indx_data_device%ixTopHyd(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixAqWat)%dat(1) = indx_data_device%ixAqWat(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgCanair)%dat(1) = indx_data_device%ixNrgCanair(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgCanopy)%dat(1) = indx_data_device%ixNrgCanopy(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydCanopy)%dat(1) = indx_data_device%ixHydCanopy(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixWatAquifer)%dat(1) = indx_data_device%ixWatAquifer(iGRU)
    indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%numberFluxCalc)%dat(1) = indx_data_device%numberFluxCalc
    ! indx_data%var(iLookINDEX%numberStateSplit)%dat(1) = indx_data_device%numberStateSplit
    ! indx_data%var(iLookINDEX%numberDomainSplitNrg)%dat(1) = indx_data_device%numberDomainSplitNrg
    ! indx_data%var(iLookINDEX%numberDomainSplitMass)%dat(1) = indx_data_device%numberDomainSplitMass
    ! indx_data%var(iLookINDEX%numberScalarSolutions)%dat(1) = indx_data_device%numberScalarSolutions

    ! if (allocated(indx_data_device%ixHydOnly)) indx_data%var(iLookINDEX%ixHydOnly)%dat = indx_data_device%ixHydOnly(:,1)
    ! if (allocated(indx_data_device%ixMatOnly)) indx_data%var(iLookINDEX%ixMatOnly)%dat = indx_data_device%ixMatOnly(:,1)
    ! if (allocated(indx_data_device%ixMassOnly)) indx_data%var(iLookINDEX%ixMassOnly)%dat = indx_data_device%ixMassOnly(:,1)
      ! if (size(indx_data%var(iLookINDEX%ixVolFracWat)%dat).ne.0) indx_data%var(iLookINDEX%ixVolFracWat)%dat = indx_data_device%ixVolFracWat(:,1)
      if (size(indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixAllState)%dat).ne.0) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixAllState)%dat = indx_data_device%ixAllState(:,iGRU)

      if (allocated(indx_data_device%ixLayerState)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixLayerState)%dat = indx_data_device%ixLayerState(1:nLayers,iGRU)

      ! if (size(indx_data%var(iLookINDEX%ixLayerActive)%dat).ne.0) indx_data%var(iLookINDEX%ixLayerActive)%dat = indx_data_device%ixLayerActive(1:nLayers,1)



        indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%numberFluxCalc)%dat(1) = indx_data_device%numberFluxCalc
        if (allocated(indx_data_device%ixNrgLayer)) then
                    ! print*, indx_data%var(iLookINDEX%ixNrgLayer)%dat
                    ! print*, 'nLayers', nLayers, size(indx_data%var(iLookINDEX%ixNrgLayer)%dat)

          indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgLayer)%dat = indx_data_device%ixNrgLayer(1:nLayers,iGRU)
          ! print*, indx_data%var(iLookINDEX%ixNrgLayer)%dat
        end if
        if (allocated(indx_data_device%ixHydLayer)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydLayer)%dat = indx_data_device%ixHydLayer(1:nLayers,iGRU)
        indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSoil)%dat(1) = indx_data_device%nSoil
        indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nLayers)%dat(1) = indx_data_device%nLayers_d(iGRU)

        if (allocated(indx_data_device%nSnow)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnow)%dat(1) = indx_data_device%nSnow(iGRU)
        if (allocated(indx_data_device%ixCasNrg)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixCasNrg)%dat(1) = indx_data_device%ixCasNrg(iGRU)
        if (allocated(indx_data_device%ixVegNrg)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixVegNrg)%dat(1) = indx_data_device%ixVegNrg(iGRU)
        if (allocated(indx_data_device%ixVegHyd)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixVegHyd)%dat(1) = indx_data_device%ixVegHyd(iGRU)
        if (allocated(indx_data_device%ixAqWat)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixAqWat)%dat(1) = indx_data_device%ixAqWat(iGRU)
        if (allocated(indx_data_device%ixTopHyd)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixTopHyd)%dat(1) = indx_data_device%ixTopHyd(iGRU)
        if (allocated(indx_data_device%ixTopNrg)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixTopNrg)%dat(1) = indx_data_device%ixTopNrg(iGRU)
        if (allocated(indx_data_device%ixSnowSoilNrg)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilNrg)%dat = indx_data_device%ixSnowSoilNrg(:,iGRU)
        if (allocated(indx_data_device%ixSnowSoilHyd )) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowSoilHyd)%dat = indx_data_device%ixSnowSoilHyd(:,iGRU)
        indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnowSoilNrg)%dat(1) = indx_data_device%nLayers_d(iGRU)
        indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnowSoilHyd)%dat(1) = indx_data_device%nLayers_d(iGRU)
        indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnowOnlyNrg)%dat(1) = indx_data_device%nSnow(iGRU)
        indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSnowOnlyHyd)%dat(1) = indx_data_device%nSnow(iGRU)
        indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSoilOnlyNrg)%dat(1) = indx_data_device%nSoil
        indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%nSoilOnlyHyd)%dat(1) = indx_data_device%nSoil
        if (allocated(indx_data_device%ixHydType)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydType)%dat = indx_data_device%ixHydType(:,iGRU)
        if (allocated(indx_data_device%ixStateType)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixStateType)%dat = indx_data_device%ixStateType(:,iGRU)
        !if (allocated(indx_data_device%ixStateType_subset)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixStateType_subset)%dat = indx_data_device%ixStateType_subset(:,iGRU)
        if (allocated(indx_data_device%ixHydCanopy)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixHydCanopy)%dat = indx_data_device%ixHydCanopy(iGRU)
        if (allocated(indx_data_device%ixNrgCanopy)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgCanopy)%dat(1) = indx_data_device%ixNrgCanopy(iGRU)
        ! if (allocated(indx_data_device%ixMapSubset2Full)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixMapSubset2Full)%dat = indx_data_device%ixMapSubset2Full(:,iGRU)
        ! if (allocated(indx_data_device%ixMapFull2Subset)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixMapFull2Subset)%dat = indx_data_device%ixMapFull2Subset(:,iGRU)
        !if (allocated(indx_data_device%ixDomainType_subset)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixDomainType_subset)%dat = indx_data_device%ixDomainType_subset(:,iGRU)
        if (allocated(indx_data_device%ixControlVolume)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixControlVolume)%dat = indx_data_device%ixControlVolume(:,iGRU)
        if (allocated(indx_data_device%layerType)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%layerType)%dat = indx_data_device%layerType(:,iGRU)
        if (allocated(indx_data_device%ixSoilOnlyHyd)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilOnlyHyd)%dat = indx_data_device%ixSoilOnlyHyd(:,iGRU)
        if (allocated(indx_data_device%ixSnowOnlyHyd)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyHyd)%dat = indx_data_device%ixSnowOnlyHyd(:,iGRU)
        if (allocated(indx_data_device%ixSoilOnlyNrg)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilOnlyNrg)%dat = indx_data_device%ixSoilOnlyNrg(:,iGRU)
        if (allocated(indx_data_device%ixSnowOnlyNrg)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSnowOnlyNrg)%dat = indx_data_device%ixSnowOnlyNrg(:,iGRU)
      
        ! if (allocated(indx_data_device%ixMatricHead)) indx_data%var(iLookINDEX%ixMatricHead)%dat = indx_data_device%ixMatricHead(:,1)
        if (allocated(indx_data_device%ixSoilState)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixSoilState)%dat = indx_data_device%ixSoilState(:,iGRU)
        if (allocated(indx_data_device%ixDomainType)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixDomainType)%dat = indx_data_device%ixDomainType(:,iGRU)
        if (allocated(indx_data_device%ixNrgOnly)) indx_data%gru(iGRU)%hru(1)%var(iLookINDEX%ixNrgOnly)%dat = indx_data_device%ixNrgOnly(:,iGRU)
      end do

      end subroutine

      
      subroutine deallocate_device_indx_data(indx_data_device)
        type(indx_data_device), intent(inout) :: indx_data_device

  if (allocated(indx_data_device%nCasNrg)) deallocate(indx_data_device%nCasNrg)
  if (allocated(indx_data_device%nVegNrg)) deallocate(indx_data_device%nVegNrg)
  if (allocated(indx_data_device%nVegMass)) deallocate(indx_data_device%nVegMass)
  if (allocated(indx_data_device%nVegState)) deallocate(indx_data_device%nVegState)
  if (allocated(indx_data_device%nNrgState)) deallocate(indx_data_device%nNrgState)
  if (allocated(indx_data_device%nWatState)) deallocate(indx_data_device%nWatState)
  if (allocated(indx_data_device%nMatState)) deallocate(indx_data_device%nMatState)
  if (allocated(indx_data_device%nMassState)) deallocate(indx_data_device%nMassState)
  if (allocated(indx_data_device%nState)) deallocate(indx_data_device%nState)
  if (allocated(indx_data_device%ixNrgCanair)) deallocate(indx_data_device%ixNrgCanair)
  if (allocated(indx_data_device%ixWatAquifer)) deallocate(indx_data_device%ixWatAquifer)
  if (allocated(indx_data_device%ixAllState)) deallocate(indx_data_device%ixAllState)
  if (allocated(indx_data_device%ixLayerState)) deallocate(indx_data_device%ixLayerState)

  if (allocated(indx_data_device%ixNrgLayer)) deallocate(indx_data_device%ixNrgLayer)
  if (allocated(indx_data_device%ixHydLayer)) deallocate(indx_data_device%ixHydLayer)
  if (allocated(indx_data_device%nSnow)) deallocate(indx_data_device%nSnow)
  if (allocated(indx_data_device%nLayers_d)) deallocate(indx_data_device%nLayers_d)
  if (allocated(indx_data_device%ixCasNrg)) deallocate(indx_data_device%ixCasNrg)
  if (allocated(indx_data_device%ixVegNrg)) deallocate(indx_data_device%ixVegNrg)
  if (allocated(indx_data_device%ixVegHyd)) deallocate(indx_data_device%ixVegHyd)
  if (allocated(indx_data_device%ixAqWat)) deallocate(indx_data_device%ixAqWat)
  if (allocated(indx_data_device%ixTopHyd)) deallocate(indx_data_device%ixTopHyd)
  if (allocated(indx_data_device%ixTopNrg)) deallocate(indx_data_device%ixTopNrg)
  if (allocated(indx_data_device%ixSnowSoilNrg)) deallocate(indx_data_device%ixSnowSoilNrg)
  if (allocated(indx_data_device%ixSnowSoilHyd)) deallocate(indx_data_device%ixSnowSoilHyd)
  if (allocated(indx_data_device%ixHydType)) deallocate(indx_data_device%ixHydType)
  if (allocated(indx_data_device%ixStateType)) deallocate(indx_data_device%ixStateType)
  !if (allocated(indx_data_device%ixStateType_subset)) deallocate(indx_data_device%ixStateType_subset)
  if (allocated(indx_data_device%ixHydCanopy)) deallocate(indx_data_device%ixHydCanopy)
  if (allocated(indx_data_device%ixNrgCanopy)) deallocate(indx_data_device%ixNrgCanopy)
  ! if (allocated(indx_data_device%ixMapSubset2Full)) deallocate(indx_data_device%ixMapSubset2Full)
  ! if (allocated(indx_data_device%ixMapFull2Subset)) deallocate(indx_data_device%ixMapFull2Subset)
  !if (allocated(indx_data_device%ixDomainType_subset)) deallocate(indx_data_device%ixDomainType_subset)
  if (allocated(indx_data_device%ixControlVolume)) deallocate(indx_data_device%ixControlVolume)
  if (allocated(indx_data_device%layerType)) deallocate(indx_data_device%layerType)
  if (allocated(indx_data_device%ixSoilOnlyHyd)) deallocate(indx_data_device%ixSoilOnlyHyd)
  if (allocated(indx_data_device%ixSnowOnlyHyd)) deallocate(indx_data_device%ixSnowOnlyHyd)
  if (allocated(indx_data_device%ixSoilOnlyNrg)) deallocate(indx_data_device%ixSoilOnlyNrg)
  if (allocated(indx_data_device%ixSnowOnlyNrg)) deallocate(indx_data_device%ixSnowOnlyNrg)

  if (allocated(indx_data_device%ixSoilState)) deallocate(indx_data_device%ixSoilState)
  if (allocated(indx_data_device%ixDomainType)) deallocate(indx_data_device%ixDomainType)
  if (allocated(indx_data_device%ixNrgOnly)) deallocate(indx_data_device%ixNrgOnly)

        ! print*, 'start deallocate', maxval(indx_data_device%nSnow)
  ! if(allocated(indx_data_device%nSnow)) deallocate(indx_data_device%nSnow)
  if(allocated(indx_data_device%nCasNrg)) deallocate(indx_data_device%nCasNrg)
  if(allocated(indx_data_device%nVegNrg)) deallocate(indx_data_device%nVegNrg)
  if(allocated(indx_data_device%nVegMass)) deallocate(indx_data_device%nVegMass)
  if(allocated(indx_data_device%nVegState)) deallocate(indx_data_device%nVegState)
  if(allocated(indx_data_device%nNrgState)) deallocate(indx_data_device%nNrgState)
  if(allocated(indx_data_device%nWatState)) deallocate(indx_data_device%nWatState)
  if(allocated(indx_data_device%nMatState)) deallocate(indx_data_device%nMatState)
  if(allocated(indx_data_device%nMassState)) deallocate(indx_data_device%nMassState)
  if(allocated(indx_data_device%nState)) deallocate(indx_data_device%nState)
  ! if(allocated(indx_data_device%nSnowSoilNrg_d)) deallocate(indx_data_device%nSnowSoilNrg_d)
  ! if(allocated(indx_data_device%nSnowOnlyNrg_d)) deallocate(indx_data_device%nSnowOnlyNrg_d)
  ! if(allocated(indx_data_device%nSoilOnlyNrg_d)) deallocate(indx_data_device%nSoilOnlyNrg_d)
  ! if(allocated(indx_data_device%nSnowSoilHyd_d)) deallocate(indx_data_device%nSnowSoilHyd_d)
  ! if(allocated(indx_data_device%nSnowOnlyHyd_d)) deallocate(indx_data_device%nSnowOnlyHyd_d)
  ! if(allocated(indx_data_device%nSoilOnlyHyd_d)) deallocate(indx_data_device%nSoilOnlyHyd_d)
  ! if(allocated(indx_data_device%ixHydOnly)) deallocate(indx_data_device%ixHydOnly)
  ! if(allocated(indx_data_device%ixMatOnly)) deallocate(indx_data_device%ixMatOnly)
  ! if(allocated(indx_data_device%ixMassOnly)) deallocate(indx_data_device%ixMassOnly)
  if(allocated(indx_data_device%ixNrgCanair)) deallocate(indx_data_device%ixNrgCanair)
  if(allocated(indx_data_device%ixWatAquifer)) deallocate(indx_data_device%ixWatAquifer)
  ! if(allocated(indx_data_device%ixVolFracWat)) deallocate(indx_data_device%ixVolFracWat)
  if(allocated(indx_data_device%ixAllState)) deallocate(indx_data_device%ixAllState)
  if(allocated(indx_data_device%ixLayerState)) deallocate(indx_data_device%ixLayerState)
  ! if(allocated(indx_data_device%ixLayerActive)) deallocate(indx_data_device%ixLayerActive)

  
        if (allocated(indx_data_device%ixNrgLayer)) deallocate(indx_data_device%ixNrgLayer)
        if (allocated(indx_data_device%ixHydLayer)) deallocate(indx_data_device%ixHydLayer)
        if (allocated(indx_data_device%nSnow)) deallocate(indx_data_device%nSnow)
        ! if (allocated(indx_data_device%nSoil_d)) deallocate(indx_data_device%nSoil_d)
        if (allocated(indx_data_device%ixCasNrg)) deallocate(indx_data_device%ixCasNrg)
        if (allocated(indx_data_device%ixVegNrg)) deallocate(indx_data_device%ixVegNrg)
        if (allocated(indx_data_device%ixVegHyd)) deallocate(indx_data_device%ixVegHyd)
        if (allocated(indx_data_device%ixAqWat)) deallocate(indx_data_device%ixAqWat)
        if (allocated(indx_data_device%ixTopHyd)) deallocate(indx_data_device%ixTopHyd)
        if (allocated(indx_data_device%ixTopNrg)) deallocate(indx_data_device%ixTopNrg)
        if (allocated(indx_data_device%ixHydType)) deallocate(indx_data_device%ixHydType)
        if (allocated(indx_data_device%ixStateType)) deallocate(indx_data_device%ixStateType)
        !if (allocated(indx_data_device%ixStateType_subset)) deallocate(indx_data_device%ixStateType_subset)
        if (allocated(indx_data_device%ixHydCanopy)) deallocate(indx_data_device%ixHydCanopy)
        if (allocated(indx_data_device%ixNrgCanopy)) deallocate(indx_data_device%ixNrgCanopy)
        ! if (allocated(indx_data_device%ixMapSubset2Full)) deallocate(indx_data_device%ixMapSubset2Full)
        ! print*, size(indx_data_device%ixMapSubset2Full,1), 's2f'
        ! if (allocated(indx_data_device%ixMapFull2Subset)) deallocate(indx_data_device%ixMapFull2Subset)
        ! print*, size(indx_data_device%ixMapFull2Subset,1), 'f2s'
        !if (allocated(indx_data_device%ixDomainType_subset)) deallocate(indx_data_device%ixDomainType_subset)
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
        ! if (allocated(indx_data_device%ixMatricHead)) deallocate(indx_data_device%ixMatricHead)
        if (allocated(indx_data_device%ixSoilState)) deallocate(indx_data_device%ixSoilState)
          
        ! print*, 'done deallocate'
      end subroutine deallocate_device_indx_data
      

  subroutine allocate_device_prog_data(prog_data_device, prog_data, nGRU, nLayers, nSoil)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    type(gru_hru_doubleVec),intent(in) :: prog_data
    integer(i4b),intent(in) :: nGRU,nLayers,nSoil

    integer(i4b) :: iGRU
    ! print*, 'prog_data', nGRU
    ! state variables for vegetation
    allocate(prog_data_device%scalarCanopyIce(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyLiq(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyWat(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanairTemp(nGRU))                ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyTemp(nGRU))                ! get_ixVarType('scalarv')
    ! state variables for snow
    allocate(prog_data_device%spectralSnowAlbedoDiffuse(nBand,nGRU))     ! get_ixVarType('wLength')
    allocate(prog_data_device%scalarSnowAlbedo(nGRU))                ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSnowDepth(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSWE(nGRU))                       ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSfcMeltPond(nGRU))               ! get_ixVarType('scalarv')
    ! define state variables for the snow+soil domain
    allocate(prog_data_device%mLayerTemp(nSoil+maxSnowLayers,nGRU))                    ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerVolFracIce(nSoil+maxSnowLayers,nGRU))              ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerVolFracLiq(nSoil+maxSnowLayers,nGRU))              ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerVolFracWat(nSoil+maxSnowLayers,nGRU))              ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerMatricHead(nSoil,nGRU))              ! get_ixVarType('midSoil')
    ! enthalpy
    allocate(prog_data_device%scalarCanairEnthalpy(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyEnthalpy(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%mLayerEnthalpy(nSoil+maxSnowLayers,nGRU))                ! get_ixVarType('midToto')
    ! other state variables
    allocate(prog_data_device%scalarAquiferStorage(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSurfaceTemp(nGRU))               ! get_ixVarType('scalarv')
    ! define coordinate variables
    allocate(prog_data_device%mLayerDepth(nSoil+maxSnowLayers,nGRU))                   ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerHeight(0:nSoil+maxSnowLayers,nGRU))                  ! get_ixVarType('ifcToto')
    allocate(prog_data_device%iLayerHeight(0:nSoil+maxSnowLayers,nGRU))                  ! get_ixVarType('ifcToto')

    do iGRU=1,nGRU
      prog_data_device%spectralSnowAlbedoDiffuse(1:nBand,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nBand)
      prog_data_device%mLayerTemp(1:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerTemp)%dat(1:nLayers)
      prog_data_device%mLayerVolFracIce(1:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracIce)%dat(1:nLayers)
      prog_data_device%mLayerVolFracLiq(1:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracLiq)%dat(1:nLayers)
      prog_data_device%mLayerVolFracWat(1:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracWat)%dat(1:nLayers)
      prog_data_device%mLayerMatricHead(1:nSoil,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerMatricHead)%dat(1:nSoil)
      prog_data_device%mLayerEnthalpy(1:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerEnthalpy)%dat(1:nLayers)
      prog_data_device%mLayerDepth(1:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerDepth)%dat(1:nLayers)
      prog_data_device%mLayerHeight(0:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerHeight)%dat(0:nLayers)
      prog_data_device%iLayerHeight(0:nLayers,iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%iLayerHeight)%dat(0:nLayers)

    enddo
    prog_data_device%dt_init = prog_data%gru(1)%hru(1)%var(iLookPROG%dt_init)%dat(1)

    do iGRU=1,nGRU
    prog_data_device%scalarCanopyIce(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyIce)%dat(1)
    prog_data_device%scalarCanopyLiq(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyLiq)%dat(1)
    prog_data_device%scalarCanopyWat(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyWat)%dat(1)
    prog_data_device%scalarCanairTemp(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanairTemp)%dat(1)
    prog_data_device%scalarCanopyTemp(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyTemp)%dat(1)
    prog_data_device%scalarSnowAlbedo(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSnowAlbedo)%dat(1)
    prog_data_device%scalarSnowDepth(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSnowDepth)%dat(1)
    ! prog_data_device%scalarSWE(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSWE)%dat(1)
    prog_data_device%scalarSfcMeltPond(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSfcMeltPond)%dat(1)
    prog_data_device%scalarCanairEnthalpy(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanairEnthalpy)%dat(1)
    prog_data_device%scalarCanopyEnthalpy(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarCanopyEnthalpy)%dat(1)
    prog_data_device%scalarAquiferStorage(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarAquiferStorage)%dat(1)
    prog_data_device%scalarSurfaceTemp(iGRU) = prog_data%gru(iGRU)%hru(1)%var(iLookPROG%scalarSurfaceTemp)%dat(1)
    end do
  end subroutine allocate_device_prog_data  

  subroutine allocate_device_prog_temp(prog_data_device, nGRU, nLayers, nSoil)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    integer(i4b),intent(in) :: nGRU,nLayers,nSoil

        ! print*, 'prog_temp', nGRU

    ! state variables for vegetation
    allocate(prog_data_device%scalarCanopyIce(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyLiq(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyWat(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanairTemp(nGRU))                ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyTemp(nGRU))                ! get_ixVarType('scalarv')
    ! state variables for snow
    allocate(prog_data_device%spectralSnowAlbedoDiffuse(nBand,nGRU))     ! get_ixVarType('wLength')
    allocate(prog_data_device%scalarSnowAlbedo(nGRU))                ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSnowDepth(nGRU))                 ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSWE(nGRU))                       ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSfcMeltPond(nGRU))               ! get_ixVarType('scalarv')
    ! define state variables for the snow+soil domain
    allocate(prog_data_device%mLayerTemp(nSoil+maxSnowLayers,nGRU))                    ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerVolFracIce(nSoil+maxSnowLayers,nGRU))              ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerVolFracLiq(nSoil+maxSnowLayers,nGRU))              ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerVolFracWat(nSoil+maxSnowLayers,nGRU))              ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerMatricHead(nSoil,nGRU))              ! get_ixVarType('midSoil')
    ! enthalpy
    allocate(prog_data_device%scalarCanairEnthalpy(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarCanopyEnthalpy(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%mLayerEnthalpy(nSoil+maxSnowLayers,nGRU))                ! get_ixVarType('midToto')
    ! other state variables
    allocate(prog_data_device%scalarAquiferStorage(nGRU))            ! get_ixVarType('scalarv')
    allocate(prog_data_device%scalarSurfaceTemp(nGRU))               ! get_ixVarType('scalarv')
    ! define coordinate variables
    allocate(prog_data_device%mLayerDepth(nSoil+maxSnowLayers,nGRU))                   ! get_ixVarType('midToto')
    allocate(prog_data_device%mLayerHeight(0:nSoil+maxSnowLayers,nGRU))                  ! get_ixVarType('ifcToto')
    allocate(prog_data_device%iLayerHeight(0:nSoil+maxSnowLayers,nGRU))                  ! get_ixVarType('ifcToto')
  
  end subroutine allocate_device_prog_temp 

  subroutine finalize_device_prog_data(prog_data_device, prog_data,nLayers,nSoil,nGRU)
    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    type(gru_hru_doubleVec),intent(inout) :: prog_data
    integer(i4b),intent(in) :: nLayers,nSoil,nGRU
    ! integer(i4b),intent(in) :: nGRU,nLayers,nSoil

    integer(i4b) :: iGRU

    do iGRU=1,nGRU
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat,1)/=size(prog_data_device%spectralSnowAlbedoDiffuse,1)) then; deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat); allocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nBand)); endif;
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nBand) = prog_data_device%spectralSnowAlbedoDiffuse(1:nBand,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerTemp)%dat,1)/=size(prog_data_device%mLayerTemp,1)) then; deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerTemp)%dat); allocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerTemp)%dat(1:nLayers)); endif;
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerTemp)%dat(1:nLayers) = prog_data_device%mLayerTemp(1:nLayers,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracIce)%dat,1)/=size(prog_data_device%mLayerVolFracIce,1)) then; deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracIce)%dat); allocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracIce)%dat(1:nLayers)); endif;
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracIce)%dat(1:nLayers) = prog_data_device%mLayerVolFracIce(1:nLayers,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracLiq)%dat,1)/=size(prog_data_device%mLayerVolFracLiq,1)) then; deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracLiq)%dat); allocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracLiq)%dat(1:nLayers)); endif;
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracLiq)%dat(1:nLayers) = prog_data_device%mLayerVolFracLiq(1:nLayers,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracWat)%dat,1)/=size(prog_data_device%mLayerVolFracWat,1)) then; deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracWat)%dat); allocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracWat)%dat(1:nLayers)); endif;
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerVolFracWat)%dat(1:nLayers) = prog_data_device%mLayerVolFracWat(1:nLayers,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerMatricHead)%dat,1)/=size(prog_data_device%mLayerMatricHead,1)) then; deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerMatricHead)%dat); allocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerMatricHead)%dat(1:nSoil)); endif;
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerMatricHead)%dat(1:nSoil) = prog_data_device%mLayerMatricHead(1:nSoil,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerEnthalpy)%dat,1)/=size(prog_data_device%mLayerEnthalpy,1)) then; deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerEnthalpy)%dat); allocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerEnthalpy)%dat(1:nLayers)); endif;
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerEnthalpy)%dat(1:nLayers) = prog_data_device%mLayerEnthalpy(1:nLayers,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerDepth)%dat,1)/=size(prog_data_device%mLayerDepth,1)) then; deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerDepth)%dat); allocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerDepth)%dat(1:nLayers)); endif;
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerDepth)%dat(1:nLayers) = prog_data_device%mLayerDepth(1:nLayers,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerHeight)%dat)/=nLayers+1) then; deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerHeight)%dat); allocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerHeight)%dat(0:nLayers)); endif
    prog_data%gru(iGRU)%hru(1)%var(iLookPROG%mLayerHeight)%dat(0:nLayers) = prog_data_device%mLayerHeight(0:nLayers,iGRU)
    if (size(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%iLayerHeight)%dat)/=nLayers+1) then; deallocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%iLayerHeight)%dat); allocate(prog_data%gru(iGRU)%hru(1)%var(iLookPROG%iLayerHeight)%dat(0:nLayers)); endif
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
      

  subroutine allocate_device_flux_data(flux_data_device, flux_data,nSnow,nSoil,nGRU,nLayers)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    type(flux_data_device), intent(inout) :: flux_data_device
    type(gru_hru_doubleVec), intent(in) :: flux_data
    integer(i4b),intent(in) :: nSnow,nSoil,nGRU,nLayers

    integer(i4b) :: iGRU
    ! print*, maxSnowLayers
  
    allocate(flux_data_device%scalarCanairNetNrgFlux(nGRU))
    allocate(flux_data_device%scalarCanopyNetNrgFlux(nGRU))
    allocate(flux_data_device%scalarGroundNetNrgFlux(nGRU))
    allocate(flux_data_device%scalarCanopyNetLiqFlux(nGRU))
    allocate(flux_data_device%scalarRainfall(nGRU))
    allocate(flux_data_device%scalarSnowfall(nGRU))
    allocate(flux_data_device%spectralIncomingDirect(nBand,nGRU))
    allocate(flux_data_device%spectralIncomingDiffuse(nBand,nGRU))
    allocate(flux_data_device%scalarCanopySunlitPAR(nGRU))
    allocate(flux_data_device%scalarCanopyShadedPAR(nGRU))
    allocate(flux_data_device%spectralBelowCanopyDirect(nBand,nGRU))
    allocate(flux_data_device%spectralBelowCanopyDiffuse(nBand,nGRU))
    allocate(flux_data_device%scalarBelowCanopySolar(nGRU))
    allocate(flux_data_device%scalarCanopyAbsorbedSolar(nGRU))
    allocate(flux_data_device%scalarGroundAbsorbedSolar(nGRU))
    allocate(flux_data_device%scalarLWRadCanopy(nGRU))
    allocate(flux_data_device%scalarLWRadGround(nGRU))
    allocate(flux_data_device%scalarLWRadUbound2Canopy(nGRU))
    allocate(flux_data_device%scalarLWRadUbound2Ground(nGRU))
    allocate(flux_data_device%scalarLWRadUbound2Ubound(nGRU))
    allocate(flux_data_device%scalarLWRadCanopy2Ubound(nGRU))
    allocate(flux_data_device%scalarLWRadCanopy2Ground(nGRU))
    allocate(flux_data_device%scalarLWRadCanopy2Canopy(nGRU))
    allocate(flux_data_device%scalarLWRadGround2Ubound(nGRU))
    allocate(flux_data_device%scalarLWRadGround2Canopy(nGRU))
    allocate(flux_data_device%scalarLWNetCanopy(nGRU))
    allocate(flux_data_device%scalarLWNetGround(nGRU))
    allocate(flux_data_device%scalarLWNetUbound(nGRU))
    allocate(flux_data_device%scalarEddyDiffusCanopyTop(nGRU))
    allocate(flux_data_device%scalarFrictionVelocity(nGRU))
    allocate(flux_data_device%scalarWindspdCanopyTop(nGRU))
    allocate(flux_data_device%scalarWindspdCanopyBottom(nGRU))
    allocate(flux_data_device%scalarGroundResistance(nGRU))
    allocate(flux_data_device%scalarCanopyResistance(nGRU))
    allocate(flux_data_device%scalarLeafResistance(nGRU))
    allocate(flux_data_device%scalarSoilResistance(nGRU))
    allocate(flux_data_device%scalarSenHeatTotal(nGRU))
    allocate(flux_data_device%scalarSenHeatCanopy(nGRU))
    allocate(flux_data_device%scalarSenHeatGround(nGRU))
    allocate(flux_data_device%scalarLatHeatTotal(nGRU))
    allocate(flux_data_device%scalarLatHeatCanopyEvap(nGRU))
    allocate(flux_data_device%scalarLatHeatCanopyTrans(nGRU))
    allocate(flux_data_device%scalarLatHeatGround(nGRU))
    allocate(flux_data_device%scalarCanopyAdvectiveHeatFlux(nGRU))
    allocate(flux_data_device%scalarGroundAdvectiveHeatFlux(nGRU))
    allocate(flux_data_device%scalarCanopySublimation(nGRU))
    allocate(flux_data_device%scalarSnowSublimation(nGRU))
    allocate(flux_data_device%scalarStomResistSunlit(nGRU))
    allocate(flux_data_device%scalarStomResistShaded(nGRU))
    allocate(flux_data_device%scalarPhotosynthesisSunlit(nGRU))
    allocate(flux_data_device%scalarPhotosynthesisShaded(nGRU))
    allocate(flux_data_device%scalarCanopyTranspiration(nGRU))
    allocate(flux_data_device%scalarCanopyEvaporation(nGRU))
    allocate(flux_data_device%scalarGroundEvaporation(nGRU))
    allocate(flux_data_device%mLayerTranspire_m(nSoil,nGRU))
    allocate(flux_data_device%scalarThroughfallSnow(nGRU))
    allocate(flux_data_device%scalarThroughfallRain(nGRU))
    allocate(flux_data_device%scalarCanopySnowUnloading(nGRU))
    allocate(flux_data_device%scalarCanopyLiqDrainage(nGRU))
    allocate(flux_data_device%scalarCanopyMeltFreeze(nGRU))
    allocate(flux_data_device%iLayerConductiveFlux_m(0:nSoil+maxSnowLayers,nGRU))
    allocate(flux_data_device%iLayerAdvectiveFlux_m(0:nSoil+maxSnowLayers,nGRU))
    allocate(flux_data_device%iLayerNrgFlux_m(0:nSoil+maxSnowLayers,nGRU))
    allocate(flux_data_device%mLayerNrgFlux_m(nSoil+maxSnowLayers,nGRU))
    allocate(flux_data_device%scalarSnowDrainage(nGRU))
    allocate(flux_data_device%iLayerLiqFluxSnow_m(0:maxSnowLayers,nGRU))
    ! if (nSnow.ne.0) 
    allocate(flux_data_device%mLayerLiqFluxSnow_m(maxSnowLayers,nGRU))
    allocate(flux_data_device%scalarRainPlusMelt(nGRU))
    allocate(flux_data_device%scalarMaxInfilRate(nGRU))
    allocate(flux_data_device%scalarInfiltration(nGRU))
    allocate(flux_data_device%scalarExfiltration(nGRU))
    allocate(flux_data_device%scalarSurfaceRunoff(nGRU))
    allocate(flux_data_device%mLayerSatHydCondMP_m(nSoil,nGRU))
    allocate(flux_data_device%mLayerSatHydCond_m(nSoil,nGRU))
    allocate(flux_data_device%iLayerSatHydCond_m(0:nSoil,nGRU))
    allocate(flux_data_device%mLayerHydCond_m(nSoil,nGRU))
    allocate(flux_data_device%iLayerLiqFluxSoil_m(0:nSoil,nGRU))
    allocate(flux_data_device%mLayerLiqFluxSoil_m(nSoil,nGRU))
    allocate(flux_data_device%mLayerBaseflow_m(nSoil,nGRU))
    allocate(flux_data_device%mLayerColumnInflow(nSoil,nGRU))
    allocate(flux_data_device%mLayerColumnOutflow_m(nSoil,nGRU))
    allocate(flux_data_device%scalarSoilBaseflow(nGRU))
    allocate(flux_data_device%scalarSoilDrainage(nGRU))
    allocate(flux_data_device%scalarAquiferRecharge(nGRU))
    allocate(flux_data_device%scalarAquiferTranspire(nGRU))
    allocate(flux_data_device%scalarAquiferBaseflow(nGRU))
    allocate(flux_data_device%scalarTotalET(nGRU))
    allocate(flux_data_device%scalarTotalRunoff(nGRU))
    allocate(flux_data_device%scalarNetRadiation(nGRU))
  end subroutine allocate_device_flux_data

  subroutine allocate_device_flux_prev(flux_data_device,nSoil,nGRU)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    implicit none
    type(flux_data_device), intent(inout) :: flux_data_device
    integer(i4b),intent(in) :: nSoil,nGRU

    integer(i4b) :: iGRU
     
  
    allocate(flux_data_device%scalarCanairNetNrgFlux(nGRU))
    allocate(flux_data_device%scalarCanopyNetNrgFlux(nGRU))
    allocate(flux_data_device%scalarGroundNetNrgFlux(nGRU))
    allocate(flux_data_device%scalarCanopyNetLiqFlux(nGRU))
    allocate(flux_data_device%scalarRainfall(nGRU))
    allocate(flux_data_device%scalarSnowfall(nGRU))
    allocate(flux_data_device%spectralIncomingDirect(nBand,nGRU))
    allocate(flux_data_device%spectralIncomingDiffuse(nBand,nGRU))
    allocate(flux_data_device%scalarCanopySunlitPAR(nGRU))
    allocate(flux_data_device%scalarCanopyShadedPAR(nGRU))
    allocate(flux_data_device%spectralBelowCanopyDirect(nBand,nGRU))
    allocate(flux_data_device%spectralBelowCanopyDiffuse(nBand,nGRU))
    allocate(flux_data_device%scalarBelowCanopySolar(nGRU))
    allocate(flux_data_device%scalarCanopyAbsorbedSolar(nGRU))
    allocate(flux_data_device%scalarGroundAbsorbedSolar(nGRU))
    allocate(flux_data_device%scalarLWRadCanopy(nGRU))
    allocate(flux_data_device%scalarLWRadGround(nGRU))
    allocate(flux_data_device%scalarLWRadUbound2Canopy(nGRU))
    allocate(flux_data_device%scalarLWRadUbound2Ground(nGRU))
    allocate(flux_data_device%scalarLWRadUbound2Ubound(nGRU))
    allocate(flux_data_device%scalarLWRadCanopy2Ubound(nGRU))
    allocate(flux_data_device%scalarLWRadCanopy2Ground(nGRU))
    allocate(flux_data_device%scalarLWRadCanopy2Canopy(nGRU))
    allocate(flux_data_device%scalarLWRadGround2Ubound(nGRU))
    allocate(flux_data_device%scalarLWRadGround2Canopy(nGRU))
    allocate(flux_data_device%scalarLWNetCanopy(nGRU))
    allocate(flux_data_device%scalarLWNetGround(nGRU))
    allocate(flux_data_device%scalarLWNetUbound(nGRU))
    allocate(flux_data_device%scalarEddyDiffusCanopyTop(nGRU))
    allocate(flux_data_device%scalarFrictionVelocity(nGRU))
    allocate(flux_data_device%scalarWindspdCanopyTop(nGRU))
    allocate(flux_data_device%scalarWindspdCanopyBottom(nGRU))
    allocate(flux_data_device%scalarGroundResistance(nGRU))
    allocate(flux_data_device%scalarCanopyResistance(nGRU))
    allocate(flux_data_device%scalarLeafResistance(nGRU))
    allocate(flux_data_device%scalarSoilResistance(nGRU))
    allocate(flux_data_device%scalarSenHeatTotal(nGRU))
    allocate(flux_data_device%scalarSenHeatCanopy(nGRU))
    allocate(flux_data_device%scalarSenHeatGround(nGRU))
    allocate(flux_data_device%scalarLatHeatTotal(nGRU))
    allocate(flux_data_device%scalarLatHeatCanopyEvap(nGRU))
    allocate(flux_data_device%scalarLatHeatCanopyTrans(nGRU))
    allocate(flux_data_device%scalarLatHeatGround(nGRU))
    allocate(flux_data_device%scalarCanopyAdvectiveHeatFlux(nGRU))
    allocate(flux_data_device%scalarGroundAdvectiveHeatFlux(nGRU))
    allocate(flux_data_device%scalarCanopySublimation(nGRU))
    allocate(flux_data_device%scalarSnowSublimation(nGRU))
    allocate(flux_data_device%scalarStomResistSunlit(nGRU))
    allocate(flux_data_device%scalarStomResistShaded(nGRU))
    allocate(flux_data_device%scalarPhotosynthesisSunlit(nGRU))
    allocate(flux_data_device%scalarPhotosynthesisShaded(nGRU))
    allocate(flux_data_device%scalarCanopyTranspiration(nGRU))
    allocate(flux_data_device%scalarCanopyEvaporation(nGRU))
    allocate(flux_data_device%scalarGroundEvaporation(nGRU))
    allocate(flux_data_device%mLayerTranspire_m(nSoil,nGRU))
    allocate(flux_data_device%scalarThroughfallSnow(nGRU))
    allocate(flux_data_device%scalarThroughfallRain(nGRU))
    allocate(flux_data_device%scalarCanopySnowUnloading(nGRU))
    allocate(flux_data_device%scalarCanopyLiqDrainage(nGRU))
    allocate(flux_data_device%scalarCanopyMeltFreeze(nGRU))
    allocate(flux_data_device%iLayerConductiveFlux_m(0:nSoil+maxSnowLayers,nGRU))
    allocate(flux_data_device%iLayerAdvectiveFlux_m(0:nSoil+maxSnowLayers,nGRU))
    allocate(flux_data_device%iLayerNrgFlux_m(0:nSoil+maxSnowLayers,nGRU))
    allocate(flux_data_device%mLayerNrgFlux_m(nSoil+maxSnowLayers,nGRU))
    allocate(flux_data_device%scalarSnowDrainage(nGRU))
    allocate(flux_data_device%iLayerLiqFluxSnow_m(0:maxSnowLayers,nGRU))
    allocate(flux_data_device%mLayerLiqFluxSnow_m(maxSnowLayers,nGRU))
    allocate(flux_data_device%scalarRainPlusMelt(nGRU))
    allocate(flux_data_device%scalarMaxInfilRate(nGRU))
    allocate(flux_data_device%scalarInfiltration(nGRU))
    allocate(flux_data_device%scalarExfiltration(nGRU))
    allocate(flux_data_device%scalarSurfaceRunoff(nGRU))
    allocate(flux_data_device%mLayerSatHydCondMP_m(nSoil,nGRU))
    allocate(flux_data_device%mLayerSatHydCond_m(nSoil,nGRU))
    allocate(flux_data_device%iLayerSatHydCond_m(0:nSoil,nGRU))
    allocate(flux_data_device%mLayerHydCond_m(nSoil,nGRU))
    allocate(flux_data_device%iLayerLiqFluxSoil_m(0:nSoil,nGRU))
    allocate(flux_data_device%mLayerLiqFluxSoil_m(nSoil,nGRU))
    allocate(flux_data_device%mLayerBaseflow_m(nSoil,nGRU))
    allocate(flux_data_device%mLayerColumnInflow(nSoil,nGRU))
    allocate(flux_data_device%mLayerColumnOutflow_m(nSoil,nGRU))
    allocate(flux_data_device%scalarSoilBaseflow(nGRU))
    allocate(flux_data_device%scalarSoilDrainage(nGRU))
    allocate(flux_data_device%scalarAquiferRecharge(nGRU))
    allocate(flux_data_device%scalarAquiferTranspire(nGRU))
    allocate(flux_data_device%scalarAquiferBaseflow(nGRU))
    allocate(flux_data_device%scalarTotalET(nGRU))
    allocate(flux_data_device%scalarTotalRunoff(nGRU))
    allocate(flux_data_device%scalarNetRadiation(nGRU))
    flux_data_device%spectralIncomingDirect = 0._rkind
    flux_data_device%spectralIncomingDiffuse = 0._rkind
    flux_data_device%spectralBelowCanopyDirect = 0._rkind
    flux_data_device%spectralBelowCanopyDiffuse = 0._rkind
    flux_data_device%iLayerConductiveFlux_m = 0._rkind
    flux_data_device%iLayerAdvectiveFlux_m = 0._rkind
    flux_data_device%iLayerNrgFlux_m = 0._rkind
    flux_data_device%mLayerNrgFlux_m = 0._rkind
    flux_data_device%iLayerLiqFluxSnow_m = 0._rkind
    flux_data_device%mLayerLiqFluxSnow_m = 0._rkind
    flux_data_device%mLayerSatHydCondMP_m = 0._rkind
    flux_data_device%mLayerSatHydCond_m = 0._rkind
    flux_data_device%iLayerSatHydCond_m = 0._rkind
    flux_data_device%mLayerHydCond_m = 0._rkind
    flux_data_device%iLayerLiqFluxSoil_m = 0._rkind
    flux_data_device%mLayerLiqFluxSoil_m = 0._rkind
    flux_data_device%mLayerBaseflow_m = 0._rkind
    flux_data_device%mLayerColumnInflow = 0._rkind
    flux_data_device%mLayerColumnOutflow_m = 0._rkind
    flux_data_device%scalarCanairNetNrgFlux = 0._rkind
    flux_data_device%scalarCanopyNetNrgFlux = 0._rkind
    flux_data_device%scalarGroundNetNrgFlux = 0._rkind
    flux_data_device%scalarCanopyNetLiqFlux = 0._rkind
    flux_data_device%scalarRainfall = 0._rkind
    flux_data_device%scalarSnowfall = 0._rkind
    flux_data_device%scalarCanopySunlitPAR = 0._rkind
    flux_data_device%scalarCanopyShadedPAR = 0._rkind
    flux_data_device%scalarBelowCanopySolar = 0._rkind
    flux_data_device%scalarCanopyAbsorbedSolar = 0._rkind
    flux_data_device%scalarGroundAbsorbedSolar = 0._rkind
    flux_data_device%scalarLWRadCanopy = 0._rkind
    flux_data_device%scalarLWRadGround = 0._rkind
    flux_data_device%scalarLWRadUbound2Canopy = 0._rkind
    flux_data_device%scalarLWRadUbound2Ground = 0._rkind
    flux_data_device%scalarLWRadUbound2Ubound = 0._rkind
    flux_data_device%scalarLWRadCanopy2Ubound = 0._rkind
    flux_data_device%scalarLWRadCanopy2Ground = 0._rkind
    flux_data_device%scalarLWRadCanopy2Canopy = 0._rkind
    flux_data_device%scalarLWRadGround2Ubound = 0._rkind
    flux_data_device%scalarLWRadGround2Canopy = 0._rkind
    flux_data_device%scalarLWNetCanopy = 0._rkind
    flux_data_device%scalarLWNetGround = 0._rkind
    flux_data_device%scalarLWNetUbound = 0._rkind
    flux_data_device%scalarEddyDiffusCanopyTop = 0._rkind
    flux_data_device%scalarFrictionVelocity = 0._rkind
    flux_data_device%scalarWindspdCanopyTop = 0._rkind
    flux_data_device%scalarWindspdCanopyBottom = 0._rkind
    flux_data_device%scalarGroundResistance = 0._rkind
    flux_data_device%scalarCanopyResistance = 0._rkind
    flux_data_device%scalarLeafResistance = 0._rkind
    flux_data_device%scalarSoilResistance = 0._rkind
    flux_data_device%scalarSenHeatTotal = 0._rkind
    flux_data_device%scalarSenHeatCanopy = 0._rkind
    flux_data_device%scalarSenHeatGround = 0._rkind
    flux_data_device%scalarLatHeatTotal = 0._rkind
    flux_data_device%scalarLatHeatCanopyEvap = 0._rkind
    flux_data_device%scalarLatHeatCanopyTrans = 0._rkind
    flux_data_device%scalarLatHeatGround = 0._rkind
    flux_data_device%scalarCanopyAdvectiveHeatFlux = 0._rkind
    flux_data_device%scalarGroundAdvectiveHeatFlux = 0._rkind
    flux_data_device%scalarCanopySublimation = 0._rkind
    flux_data_device%scalarSnowSublimation = 0._rkind
    flux_data_device%scalarStomResistSunlit = 0._rkind
    flux_data_device%scalarStomResistShaded = 0._rkind
    flux_data_device%scalarPhotosynthesisSunlit = 0._rkind
    flux_data_device%scalarPhotosynthesisShaded = 0._rkind
    flux_data_device%scalarCanopyTranspiration = 0._rkind
    flux_data_device%scalarCanopyEvaporation = 0._rkind
    flux_data_device%scalarGroundEvaporation = 0._rkind
    flux_data_device%scalarThroughfallSnow = 0._rkind
    flux_data_device%scalarThroughfallRain = 0._rkind
    flux_data_device%scalarCanopySnowUnloading = 0._rkind
    flux_data_device%scalarCanopyLiqDrainage = 0._rkind
    flux_data_device%scalarCanopyMeltFreeze = 0._rkind
    flux_data_device%scalarSnowDrainage = 0._rkind
    flux_data_device%scalarRainPlusMelt = 0._rkind
    flux_data_device%scalarMaxInfilRate = 0._rkind
    flux_data_device%scalarInfiltration = 0._rkind
    flux_data_device%scalarExfiltration = 0._rkind
    flux_data_device%scalarSurfaceRunoff = 0._rkind
    flux_data_device%scalarSoilBaseflow = 0._rkind
    flux_data_device%scalarSoilDrainage = 0._rkind
    flux_data_device%scalarAquiferRecharge = 0._rkind
    flux_data_device%scalarAquiferTranspire = 0._rkind
    flux_data_device%scalarAquiferBaseflow = 0._rkind
    flux_data_device%scalarTotalET = 0._rkind
    flux_data_device%scalarTotalRunoff = 0._rkind
    flux_data_device%scalarNetRadiation = 0._rkind
      
  end subroutine allocate_device_flux_prev
  subroutine zero_device_flux(flux_data_device)
    ! USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    implicit none
    type(flux_data_device), intent(inout) :: flux_data_device
    ! integer(i4b),intent(in) :: nSnow,nSoil,nGRU,nLayers

    ! integer(i4b) :: iGRU
      
  
    flux_data_device%spectralIncomingDirect = 0._rkind
    flux_data_device%spectralIncomingDiffuse = 0._rkind
    flux_data_device%spectralBelowCanopyDirect = 0._rkind
    flux_data_device%spectralBelowCanopyDiffuse = 0._rkind
    flux_data_device%iLayerConductiveFlux_m = 0._rkind
    flux_data_device%iLayerAdvectiveFlux_m = 0._rkind
    flux_data_device%iLayerNrgFlux_m = 0._rkind
    flux_data_device%mLayerNrgFlux_m = 0._rkind
    flux_data_device%iLayerLiqFluxSnow_m = 0._rkind
    flux_data_device%mLayerLiqFluxSnow_m = 0._rkind
    flux_data_device%mLayerSatHydCondMP_m = 0._rkind
    flux_data_device%mLayerSatHydCond_m = 0._rkind
    flux_data_device%iLayerSatHydCond_m = 0._rkind
    flux_data_device%mLayerHydCond_m = 0._rkind
    flux_data_device%iLayerLiqFluxSoil_m = 0._rkind
    flux_data_device%mLayerLiqFluxSoil_m = 0._rkind
    flux_data_device%mLayerBaseflow_m = 0._rkind
    flux_data_device%mLayerColumnInflow = 0._rkind
    flux_data_device%mLayerColumnOutflow_m = 0._rkind
    flux_data_device%scalarCanairNetNrgFlux = 0._rkind
    flux_data_device%scalarCanopyNetNrgFlux = 0._rkind
    flux_data_device%scalarGroundNetNrgFlux = 0._rkind
    flux_data_device%scalarCanopyNetLiqFlux = 0._rkind
    flux_data_device%scalarRainfall = 0._rkind
    flux_data_device%scalarSnowfall = 0._rkind
    flux_data_device%scalarCanopySunlitPAR = 0._rkind
    flux_data_device%scalarCanopyShadedPAR = 0._rkind
    flux_data_device%scalarBelowCanopySolar = 0._rkind
    flux_data_device%scalarCanopyAbsorbedSolar = 0._rkind
    flux_data_device%scalarGroundAbsorbedSolar = 0._rkind
    flux_data_device%scalarLWRadCanopy = 0._rkind
    flux_data_device%scalarLWRadGround = 0._rkind
    flux_data_device%scalarLWRadUbound2Canopy = 0._rkind
    flux_data_device%scalarLWRadUbound2Ground = 0._rkind
    flux_data_device%scalarLWRadUbound2Ubound = 0._rkind
    flux_data_device%scalarLWRadCanopy2Ubound = 0._rkind
    flux_data_device%scalarLWRadCanopy2Ground = 0._rkind
    flux_data_device%scalarLWRadCanopy2Canopy = 0._rkind
    flux_data_device%scalarLWRadGround2Ubound = 0._rkind
    flux_data_device%scalarLWRadGround2Canopy = 0._rkind
    flux_data_device%scalarLWNetCanopy = 0._rkind
    flux_data_device%scalarLWNetGround = 0._rkind
    flux_data_device%scalarLWNetUbound = 0._rkind
    flux_data_device%scalarEddyDiffusCanopyTop = 0._rkind
    flux_data_device%scalarFrictionVelocity = 0._rkind
    flux_data_device%scalarWindspdCanopyTop = 0._rkind
    flux_data_device%scalarWindspdCanopyBottom = 0._rkind
    flux_data_device%scalarGroundResistance = 0._rkind
    flux_data_device%scalarCanopyResistance = 0._rkind
    flux_data_device%scalarLeafResistance = 0._rkind
    flux_data_device%scalarSoilResistance = 0._rkind
    flux_data_device%scalarSenHeatTotal = 0._rkind
    flux_data_device%scalarSenHeatCanopy = 0._rkind
    flux_data_device%scalarSenHeatGround = 0._rkind
    flux_data_device%scalarLatHeatTotal = 0._rkind
    flux_data_device%scalarLatHeatCanopyEvap = 0._rkind
    flux_data_device%scalarLatHeatCanopyTrans = 0._rkind
    flux_data_device%scalarLatHeatGround = 0._rkind
    flux_data_device%scalarCanopyAdvectiveHeatFlux = 0._rkind
    flux_data_device%scalarGroundAdvectiveHeatFlux = 0._rkind
    flux_data_device%scalarCanopySublimation = 0._rkind
    flux_data_device%scalarSnowSublimation = 0._rkind
    flux_data_device%scalarStomResistSunlit = 0._rkind
    flux_data_device%scalarStomResistShaded = 0._rkind
    flux_data_device%scalarPhotosynthesisSunlit = 0._rkind
    flux_data_device%scalarPhotosynthesisShaded = 0._rkind
    flux_data_device%scalarCanopyTranspiration = 0._rkind
    flux_data_device%scalarCanopyEvaporation = 0._rkind
    flux_data_device%scalarGroundEvaporation = 0._rkind
    flux_data_device%scalarThroughfallSnow = 0._rkind
    flux_data_device%scalarThroughfallRain = 0._rkind
    flux_data_device%scalarCanopySnowUnloading = 0._rkind
    flux_data_device%scalarCanopyLiqDrainage = 0._rkind
    flux_data_device%scalarCanopyMeltFreeze = 0._rkind
    flux_data_device%scalarSnowDrainage = 0._rkind
    flux_data_device%scalarRainPlusMelt = 0._rkind
    flux_data_device%scalarMaxInfilRate = 0._rkind
    flux_data_device%scalarInfiltration = 0._rkind
    flux_data_device%scalarExfiltration = 0._rkind
    flux_data_device%scalarSurfaceRunoff = 0._rkind
    flux_data_device%scalarSoilBaseflow = 0._rkind
    flux_data_device%scalarSoilDrainage = 0._rkind
    flux_data_device%scalarAquiferRecharge = 0._rkind
    flux_data_device%scalarAquiferTranspire = 0._rkind
    flux_data_device%scalarAquiferBaseflow = 0._rkind
    flux_data_device%scalarTotalET = 0._rkind
    flux_data_device%scalarTotalRunoff = 0._rkind
    flux_data_device%scalarNetRadiation = 0._rkind
      
  end subroutine zero_device_flux

  subroutine finalize_device_flux_data(flux_data_device, flux_data,nSnow,nSoil,nGRU)
    implicit none
    type(flux_data_device), intent(inout) :: flux_data_device
    type(gru_hru_doubleVec), intent(inout) :: flux_data
    integer(i4b),intent(in) :: nSnow,nSoil,nGRU
    integer(i4b) :: iGRU

    do iGRU=1,nGRU
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%spectralIncomingDirect)%dat = flux_data_device%spectralIncomingDirect(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%spectralIncomingDiffuse)%dat = flux_data_device%spectralIncomingDiffuse(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%spectralBelowCanopyDirect)%dat = flux_data_device%spectralBelowCanopyDirect(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%spectralBelowCanopyDiffuse)%dat = flux_data_device%spectralBelowCanopyDiffuse(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerConductiveFlux)%dat = flux_data_device%iLayerConductiveFlux_m(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerAdvectiveFlux)%dat = flux_data_device%iLayerAdvectiveFlux_m(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerNrgFlux)%dat = flux_data_device%iLayerNrgFlux_m(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerNrgFlux)%dat = flux_data_device%mLayerNrgFlux_m(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerLiqFluxSnow)%dat = flux_data_device%iLayerLiqFluxSnow_m(0:nSnow,iGRU)
    if (nSnow.ne.0) flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerLiqFluxSnow)%dat = flux_data_device%mLayerLiqFluxSnow_m(1:nSnow,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerSatHydCondMP)%dat = flux_data_device%mLayerSatHydCondMP_m(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerSatHydCond)%dat = flux_data_device%mLayerSatHydCond_m(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerSatHydCond)%dat = flux_data_device%iLayerSatHydCond_m(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerHydCond)%dat = flux_data_device%mLayerHydCond_m(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%iLayerLiqFluxSoil)%dat = flux_data_device%iLayerLiqFluxSoil_m(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerLiqFluxSoil)%dat = flux_data_device%mLayerLiqFluxSoil_m(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerBaseflow)%dat = flux_data_device%mLayerBaseflow_m(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerColumnInflow)%dat = flux_data_device%mLayerColumnInflow(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%mLayerColumnOutflow)%dat = flux_data_device%mLayerColumnOutflow_m(:,iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1) = flux_data_device%scalarCanairNetNrgFlux(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1) = flux_data_device%scalarCanopyNetNrgFlux(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1) = flux_data_device%scalarGroundNetNrgFlux(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1) = flux_data_device%scalarCanopyNetLiqFlux(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarRainfall)%dat(1) = flux_data_device%scalarRainfall(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSnowfall)%dat(1) = flux_data_device%scalarSnowfall(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopySunlitPAR)%dat(1) = flux_data_device%scalarCanopySunlitPAR(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyShadedPAR)%dat(1) = flux_data_device%scalarCanopyShadedPAR(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarBelowCanopySolar)%dat(1) = flux_data_device%scalarBelowCanopySolar(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyAbsorbedSolar)%dat(1) = flux_data_device%scalarCanopyAbsorbedSolar(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1) = flux_data_device%scalarGroundAbsorbedSolar(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadCanopy)%dat(1) = flux_data_device%scalarLWRadCanopy(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadGround)%dat(1) = flux_data_device%scalarLWRadGround(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadUbound2Canopy)%dat(1) = flux_data_device%scalarLWRadUbound2Canopy(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadUbound2Ground)%dat(1) = flux_data_device%scalarLWRadUbound2Ground(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadUbound2Ubound)%dat(1) = flux_data_device%scalarLWRadUbound2Ubound(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadCanopy2Ubound)%dat(1) = flux_data_device%scalarLWRadCanopy2Ubound(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadCanopy2Ground)%dat(1) = flux_data_device%scalarLWRadCanopy2Ground(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadCanopy2Canopy)%dat(1) = flux_data_device%scalarLWRadCanopy2Canopy(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadGround2Ubound)%dat(1) = flux_data_device%scalarLWRadGround2Ubound(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWRadGround2Canopy)%dat(1) = flux_data_device%scalarLWRadGround2Canopy(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWNetCanopy)%dat(1) = flux_data_device%scalarLWNetCanopy(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWNetGround)%dat(1) = flux_data_device%scalarLWNetGround(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLWNetUbound)%dat(1) = flux_data_device%scalarLWNetUbound(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarEddyDiffusCanopyTop)%dat(1) = flux_data_device%scalarEddyDiffusCanopyTop(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarFrictionVelocity)%dat(1) = flux_data_device%scalarFrictionVelocity(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarWindspdCanopyTop)%dat(1) = flux_data_device%scalarWindspdCanopyTop(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarWindspdCanopyBottom)%dat(1) = flux_data_device%scalarWindspdCanopyBottom(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundResistance)%dat(1) = flux_data_device%scalarGroundResistance(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyResistance)%dat(1) = flux_data_device%scalarCanopyResistance(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLeafResistance)%dat(1) = flux_data_device%scalarLeafResistance(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSoilResistance)%dat(1) = flux_data_device%scalarSoilResistance(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSenHeatTotal)%dat(1) = flux_data_device%scalarSenHeatTotal(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSenHeatCanopy)%dat(1) = flux_data_device%scalarSenHeatCanopy(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSenHeatGround)%dat(1) = flux_data_device%scalarSenHeatGround(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLatHeatTotal)%dat(1) = flux_data_device%scalarLatHeatTotal(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1) = flux_data_device%scalarLatHeatCanopyEvap(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLatHeatCanopyTrans)%dat(1) = flux_data_device%scalarLatHeatCanopyTrans(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarLatHeatGround)%dat(1) = flux_data_device%scalarLatHeatGround(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyAdvectiveHeatFlux)%dat(1) = flux_data_device%scalarCanopyAdvectiveHeatFlux(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundAdvectiveHeatFlux)%dat(1) = flux_data_device%scalarGroundAdvectiveHeatFlux(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopySublimation)%dat(1) = flux_data_device%scalarCanopySublimation(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSnowSublimation)%dat(1) = flux_data_device%scalarSnowSublimation(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarStomResistSunlit)%dat(1) = flux_data_device%scalarStomResistSunlit(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarStomResistShaded)%dat(1) = flux_data_device%scalarStomResistShaded(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarPhotosynthesisSunlit)%dat(1) = flux_data_device%scalarPhotosynthesisSunlit(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarPhotosynthesisShaded)%dat(1) = flux_data_device%scalarPhotosynthesisShaded(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyTranspiration)%dat(1) = flux_data_device%scalarCanopyTranspiration(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyEvaporation)%dat(1) = flux_data_device%scalarCanopyEvaporation(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarGroundEvaporation)%dat(1) = flux_data_device%scalarGroundEvaporation(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarThroughfallSnow)%dat(1) = flux_data_device%scalarThroughfallSnow(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarThroughfallRain)%dat(1) = flux_data_device%scalarThroughfallRain(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1) = flux_data_device%scalarCanopySnowUnloading(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) = flux_data_device%scalarCanopyLiqDrainage(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarCanopyMeltFreeze)%dat(1) = flux_data_device%scalarCanopyMeltFreeze(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSnowDrainage)%dat(1) = flux_data_device%scalarSnowDrainage(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarRainPlusMelt)%dat(1) = flux_data_device%scalarRainPlusMelt(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarMaxInfilRate)%dat(1) = flux_data_device%scalarMaxInfilRate(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarInfiltration)%dat(1) = flux_data_device%scalarInfiltration(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarExfiltration)%dat(1) = flux_data_device%scalarExfiltration(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSurfaceRunoff)%dat(1) = flux_data_device%scalarSurfaceRunoff(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSoilBaseflow)%dat(1) = flux_data_device%scalarSoilBaseflow(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarSoilDrainage)%dat(1) = flux_data_device%scalarSoilDrainage(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarAquiferRecharge)%dat(1) = flux_data_device%scalarAquiferRecharge(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarAquiferTranspire)%dat(1) = flux_data_device%scalarAquiferTranspire(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarAquiferBaseflow)%dat(1) = flux_data_device%scalarAquiferBaseflow(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarTotalET)%dat(1) = flux_data_device%scalarTotalET(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarTotalRunoff)%dat(1) = flux_data_device%scalarTotalRunoff(iGRU)
    flux_data%gru(iGRU)%hru(1)%var(iLookFLUX%scalarNetRadiation)%dat(1) = flux_data_device%scalarNetRadiation(iGRU)
    end do
  end subroutine finalize_device_flux_data

  subroutine deallocate_device_flux_data(flux_data_device)
    implicit none
    type(flux_data_device), intent(inout) :: flux_data_device

    integer(i4b) :: iGRU

    deallocate(flux_data_device%scalarCanairNetNrgFlux)
    deallocate(flux_data_device%scalarCanopyNetNrgFlux)
    deallocate(flux_data_device%scalarGroundNetNrgFlux)
    deallocate(flux_data_device%scalarCanopyNetLiqFlux)
    deallocate(flux_data_device%scalarRainfall)
    deallocate(flux_data_device%scalarSnowfall)
    deallocate(flux_data_device%spectralIncomingDirect)
    deallocate(flux_data_device%spectralIncomingDiffuse)
    deallocate(flux_data_device%scalarCanopySunlitPAR)
    deallocate(flux_data_device%scalarCanopyShadedPAR)
    deallocate(flux_data_device%spectralBelowCanopyDirect)
    deallocate(flux_data_device%spectralBelowCanopyDiffuse)
    deallocate(flux_data_device%scalarBelowCanopySolar)
    deallocate(flux_data_device%scalarCanopyAbsorbedSolar)
    deallocate(flux_data_device%scalarGroundAbsorbedSolar)
    deallocate(flux_data_device%scalarLWRadCanopy)
    deallocate(flux_data_device%scalarLWRadGround)
    deallocate(flux_data_device%scalarLWRadUbound2Canopy)
    deallocate(flux_data_device%scalarLWRadUbound2Ground)
    deallocate(flux_data_device%scalarLWRadUbound2Ubound)
    deallocate(flux_data_device%scalarLWRadCanopy2Ubound)
    deallocate(flux_data_device%scalarLWRadCanopy2Ground)
    deallocate(flux_data_device%scalarLWRadCanopy2Canopy)
    deallocate(flux_data_device%scalarLWRadGround2Ubound)
    deallocate(flux_data_device%scalarLWRadGround2Canopy)
    deallocate(flux_data_device%scalarLWNetCanopy)
    deallocate(flux_data_device%scalarLWNetGround)
    deallocate(flux_data_device%scalarLWNetUbound)
    deallocate(flux_data_device%scalarEddyDiffusCanopyTop)
    deallocate(flux_data_device%scalarFrictionVelocity)
    deallocate(flux_data_device%scalarWindspdCanopyTop)
    deallocate(flux_data_device%scalarWindspdCanopyBottom)
    deallocate(flux_data_device%scalarGroundResistance)
    deallocate(flux_data_device%scalarCanopyResistance)
    deallocate(flux_data_device%scalarLeafResistance)
    deallocate(flux_data_device%scalarSoilResistance)
    deallocate(flux_data_device%scalarSenHeatTotal)
    deallocate(flux_data_device%scalarSenHeatCanopy)
    deallocate(flux_data_device%scalarSenHeatGround)
    deallocate(flux_data_device%scalarLatHeatTotal)
    deallocate(flux_data_device%scalarLatHeatCanopyEvap)
    deallocate(flux_data_device%scalarLatHeatCanopyTrans)
    deallocate(flux_data_device%scalarLatHeatGround)
    deallocate(flux_data_device%scalarCanopyAdvectiveHeatFlux)
    deallocate(flux_data_device%scalarGroundAdvectiveHeatFlux)
    deallocate(flux_data_device%scalarCanopySublimation)
    deallocate(flux_data_device%scalarSnowSublimation)
    deallocate(flux_data_device%scalarStomResistSunlit)
    deallocate(flux_data_device%scalarStomResistShaded)
    deallocate(flux_data_device%scalarPhotosynthesisSunlit)
    deallocate(flux_data_device%scalarPhotosynthesisShaded)
    deallocate(flux_data_device%scalarCanopyTranspiration)
    deallocate(flux_data_device%scalarCanopyEvaporation)
    deallocate(flux_data_device%scalarGroundEvaporation)
    deallocate(flux_data_device%mLayerTranspire_m)
    deallocate(flux_data_device%scalarThroughfallSnow)
    deallocate(flux_data_device%scalarThroughfallRain)
    deallocate(flux_data_device%scalarCanopySnowUnloading)
    deallocate(flux_data_device%scalarCanopyLiqDrainage)
    deallocate(flux_data_device%scalarCanopyMeltFreeze)
    deallocate(flux_data_device%iLayerConductiveFlux_m)
    deallocate(flux_data_device%iLayerAdvectiveFlux_m)
    deallocate(flux_data_device%iLayerNrgFlux_m)
    deallocate(flux_data_device%mLayerNrgFlux_m)
    deallocate(flux_data_device%scalarSnowDrainage)
    deallocate(flux_data_device%iLayerLiqFluxSnow_m)
    deallocate(flux_data_device%mLayerLiqFluxSnow_m)
    deallocate(flux_data_device%scalarRainPlusMelt)
    deallocate(flux_data_device%scalarMaxInfilRate)
    deallocate(flux_data_device%scalarInfiltration)
    deallocate(flux_data_device%scalarExfiltration)
    deallocate(flux_data_device%scalarSurfaceRunoff)
    deallocate(flux_data_device%mLayerSatHydCondMP_m)
    deallocate(flux_data_device%mLayerSatHydCond_m)
    deallocate(flux_data_device%iLayerSatHydCond_m)
    deallocate(flux_data_device%mLayerHydCond_m)
    deallocate(flux_data_device%iLayerLiqFluxSoil_m)
    deallocate(flux_data_device%mLayerLiqFluxSoil_m)
    deallocate(flux_data_device%mLayerBaseflow_m)
    deallocate(flux_data_device%mLayerColumnInflow)
    deallocate(flux_data_device%mLayerColumnOutflow_m)
    deallocate(flux_data_device%scalarSoilBaseflow)
    deallocate(flux_data_device%scalarSoilDrainage)
    deallocate(flux_data_device%scalarAquiferRecharge)
    deallocate(flux_data_device%scalarAquiferTranspire)
    deallocate(flux_data_device%scalarAquiferBaseflow)
    deallocate(flux_data_device%scalarTotalET)
    deallocate(flux_data_device%scalarTotalRunoff)
    deallocate(flux_data_device%scalarNetRadiation)
  end subroutine deallocate_device_flux_data
      
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
        allocate(deriv_data_device%dFracLiqSnow_dTk_m(nLayers,nGRU))
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
          deriv_data_device%dFracLiqSnow_dTk_m = 0._rkind
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
        deallocate(deriv_data_device%dFracLiqSnow_dTk_m)
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
      
subroutine allocate_device_diag_data(diag_data_device,diag_data,nSnow, nGRU, nLayers, nSoil)
    use globalData,only:maxSnowLayers
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_device
    type(gru_hru_doubleVec),intent(in) :: diag_data
    integer(i4b),intent(in) :: nSnow
    integer(i4b),intent(in) :: nGRU,nLayers,nSoil

    integer(i4b) :: iGRU

    allocate(diag_data_device%scalarCanopyDepth(nGRU))
    allocate(diag_data_device%scalarGreenVegFraction(nGRU))
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
    allocate(diag_data_device%mLayerVolHtCapBulk_m(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%mLayerCm_m(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%scalarLambda_drysoil(nGRU))
    allocate(diag_data_device%scalarLambda_wetsoil(nGRU))
    allocate(diag_data_device%mLayerThermalC_m(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%iLayerThermalC_m(0:nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%scalarCanopyEnthTemp(nGRU))
    allocate(diag_data_device%mLayerEnthTemp(nSoil+maxSnowLayers,nGRU))
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
    allocate(diag_data_device%spectralAlbGndDirect(nBand,nGRU))
    allocate(diag_data_device%spectralAlbGndDiffuse(nBand,nGRU))
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
    allocate(diag_data_device%mLayerTranspireLim_m(nSoil,nGRU))
    allocate(diag_data_device%mLayerRootDensity_m(nSoil,nGRU))
    allocate(diag_data_device%scalarAquiferRootFrac(nGRU))
    allocate(diag_data_device%scalarFracLiqVeg(nGRU))
    allocate(diag_data_device%scalarCanopyWetFraction(nGRU))
    allocate(diag_data_device%scalarSnowAge(nGRU))
    allocate(diag_data_device%scalarGroundSnowFraction(nGRU))
    allocate(diag_data_device%spectralSnowAlbedoDirect(nBand,nGRU))
    allocate(diag_data_device%mLayerFracLiqSnow_m(maxSnowLayers,nGRU))
    allocate(diag_data_device%mLayerThetaResid_m(maxSnowLayers,nGRU))
    allocate(diag_data_device%mLayerPoreSpace_m(maxSnowLayers,nGRU))
    allocate(diag_data_device%mLayerMeltFreeze(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%scalarInfilArea(nGRU))
    allocate(diag_data_device%scalarFrozenArea(nGRU))
    allocate(diag_data_device%scalarSoilControl(nGRU))
    allocate(diag_data_device%mLayerVolFracAir_m(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%mLayerTcrit(nSoil,nGRU))
    allocate(diag_data_device%mLayerCompress_m(nSoil,nGRU))
    allocate(diag_data_device%scalarSoilCompress(nGRU))
    allocate(diag_data_device%mLayerMatricHeadLiq(nSoil,nGRU))
    allocate(diag_data_device%scalarTotalSoilLiq(nGRU))
    allocate(diag_data_device%scalarTotalSoilIce(nGRU))
    allocate(diag_data_device%scalarTotalSoilWat(nGRU))
    allocate(diag_data_device%scalarVGn_m_m(nSoil,nGRU))
    allocate(diag_data_device%scalarKappa(nGRU))
    allocate(diag_data_device%scalarVolLatHt_fus(nGRU))
    allocate(diag_data_device%balanceCasNrg(nGRU))
    allocate(diag_data_device%balanceVegNrg(nGRU))
    allocate(diag_data_device%balanceLayerNrg(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%balanceSnowNrg(nGRU))
    allocate(diag_data_device%balanceSoilNrg(nGRU))
    allocate(diag_data_device%balanceVegMass(nGRU))
    allocate(diag_data_device%balanceLayerMass(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%balanceSnowMass(nGRU))
    allocate(diag_data_device%balanceSoilMass(nGRU))
    allocate(diag_data_device%balanceAqMass(nGRU))
  
    diag_data_device%numFluxCalls = diag_data%gru(1)%hru(1)%var(iLookDIAG%numFluxCalls)%dat(1)
    diag_data_device%wallClockTime = diag_data%gru(1)%hru(1)%var(iLookDIAG%wallClockTime)%dat(1)
    diag_data_device%meanStepSize = diag_data%gru(1)%hru(1)%var(iLookDIAG%meanStepSize)%dat(1)

    diag_data_device%numSteps = diag_data%gru(1)%hru(1)%var(iLookDIAG%numSteps)%dat(1)
    diag_data_device%numResEvals = diag_data%gru(1)%hru(1)%var(iLookDIAG%numResEvals)%dat(1)
    diag_data_device%numLinSolvSetups = diag_data%gru(1)%hru(1)%var(iLookDIAG%numLinSolvSetups)%dat(1)
    diag_data_device%numErrTestFails = diag_data%gru(1)%hru(1)%var(iLookDIAG%numErrTestFails)%dat(1)
    diag_data_device%kLast = diag_data%gru(1)%hru(1)%var(iLookDIAG%kLast)%dat(1)
    diag_data_device%kCur = diag_data%gru(1)%hru(1)%var(iLookDIAG%kCur)%dat(1)
    diag_data_device%hInitUsed = diag_data%gru(1)%hru(1)%var(iLookDIAG%hInitUsed)%dat(1)
    diag_data_device%hLast = diag_data%gru(1)%hru(1)%var(iLookDIAG%hLast)%dat(1)
    diag_data_device%hCur = diag_data%gru(1)%hru(1)%var(iLookDIAG%hCur)%dat(1)
    diag_data_device%tCur = diag_data%gru(1)%hru(1)%var(iLookDIAG%tCur)%dat(1)

    do iGRU=1,nGRU
    ! diag_data_device%scalarCanopyDepth(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyDepth)%dat(1)
    ! diag_data_device%scalarGreenVegFraction(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGreenVegFraction)%dat(1)
    ! diag_data_device%scalarBulkVolHeatCapVeg(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1)
    ! diag_data_device%scalarCanopyCm(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyCm)%dat(1)
    ! diag_data_device%scalarCanopyEmissivity(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyEmissivity)%dat(1)
    ! diag_data_device%scalarRootZoneTemp(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarRootZoneTemp)%dat(1)
    ! diag_data_device%scalarLAI(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLAI)%dat(1)
    ! diag_data_device%scalarSAI(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSAI)%dat(1)
    ! diag_data_device%scalarExposedLAI(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarExposedLAI)%dat(1)
    ! diag_data_device%scalarExposedSAI(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarExposedSAI)%dat(1)
    ! diag_data_device%scalarAdjMeasHeight(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarAdjMeasHeight)%dat(1)
    ! diag_data_device%scalarCanopyIceMax(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyIceMax)%dat(1)
    ! diag_data_device%scalarCanopyLiqMax(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyLiqMax)%dat(1)
    ! diag_data_device%scalarGrowingSeasonIndex(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGrowingSeasonIndex)%dat(1)
    ! diag_data_device%scalarVolHtCap_air(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolHtCap_air)%dat(1)
    ! diag_data_device%scalarVolHtCap_ice(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolHtCap_ice)%dat(1)
    ! diag_data_device%scalarVolHtCap_soil(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolHtCap_soil)%dat(1)
    ! diag_data_device%scalarVolHtCap_water(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolHtCap_water)%dat(1)
    ! diag_data_device%scalarLambda_drysoil(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLambda_drysoil)%dat(1)
    ! diag_data_device%scalarLambda_wetsoil(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLambda_wetsoil)%dat(1)
    ! diag_data_device%scalarCanopyEnthTemp(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1)
    ! diag_data_device%scalarTotalSoilEnthalpy(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilEnthalpy)%dat(1)
    ! diag_data_device%scalarTotalSnowEnthalpy(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSnowEnthalpy)%dat(1)
    ! diag_data_device%scalarVPair(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVPair)%dat(1)
    ! diag_data_device%scalarVP_CanopyAir(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVP_CanopyAir)%dat(1)
    ! diag_data_device%scalarTwetbulb(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTwetbulb)%dat(1)
    ! diag_data_device%scalarSnowfallTemp(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSnowfallTemp)%dat(1)
    ! diag_data_device%scalarNewSnowDensity(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarNewSnowDensity)%dat(1)
    ! diag_data_device%scalarO2air(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarO2air)%dat(1)
    ! diag_data_device%scalarCO2air(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCO2air)%dat(1)
    ! diag_data_device%windspd_x(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%windspd_x)%dat(1)
    ! diag_data_device%windspd_y(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%windspd_y)%dat(1)
    ! diag_data_device%scalarCosZenith(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCosZenith)%dat(1)
    ! diag_data_device%scalarFractionDirect(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFractionDirect)%dat(1)
    ! diag_data_device%scalarCanopySunlitFraction(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopySunlitFraction)%dat(1)
    ! diag_data_device%scalarCanopySunlitLAI(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopySunlitLAI)%dat(1)
    ! diag_data_device%scalarCanopyShadedLAI(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyShadedLAI)%dat(1)
    ! diag_data_device%scalarGroundAlbedo(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGroundAlbedo)%dat(1)
    ! diag_data_device%scalarLatHeatSubVapCanopy(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLatHeatSubVapCanopy)%dat(1)
    ! diag_data_device%scalarLatHeatSubVapGround(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLatHeatSubVapGround)%dat(1)
    ! diag_data_device%scalarSatVP_CanopyTemp(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSatVP_CanopyTemp)%dat(1)
    ! diag_data_device%scalarSatVP_GroundTemp(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSatVP_GroundTemp)%dat(1)
    ! diag_data_device%scalarZ0Canopy(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarZ0Canopy)%dat(1)
    ! diag_data_device%scalarWindReductionFactor(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarWindReductionFactor)%dat(1)
    ! diag_data_device%scalarZeroPlaneDisplacement(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarZeroPlaneDisplacement)%dat(1)
    ! diag_data_device%scalarRiBulkCanopy(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarRiBulkCanopy)%dat(1)
    ! diag_data_device%scalarRiBulkGround(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarRiBulkGround)%dat(1)
    ! diag_data_device%scalarCanopyStabilityCorrection(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyStabilityCorrection)%dat(1)
    ! diag_data_device%scalarGroundStabilityCorrection(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGroundStabilityCorrection)%dat(1)
    ! diag_data_device%scalarIntercellularCO2Sunlit(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarIntercellularCO2Sunlit)%dat(1)
    ! diag_data_device%scalarIntercellularCO2Shaded(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarIntercellularCO2Shaded)%dat(1)
    ! diag_data_device%scalarTranspireLim(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTranspireLim)%dat(1)
    ! diag_data_device%scalarTranspireLimAqfr(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTranspireLimAqfr)%dat(1)
    ! diag_data_device%scalarFoliageNitrogenFactor(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFoliageNitrogenFactor)%dat(1)
    ! diag_data_device%scalarSoilRelHumidity(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSoilRelHumidity)%dat(1)
    ! diag_data_device%scalarAquiferRootFrac(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarAquiferRootFrac)%dat(1)
    ! diag_data_device%scalarFracLiqVeg(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFracLiqVeg)%dat(1)
    ! diag_data_device%scalarCanopyWetFraction(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyWetFraction)%dat(1)
    ! diag_data_device%scalarSnowAge(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSnowAge)%dat(1)
    ! diag_data_device%scalarGroundSnowFraction(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGroundSnowFraction)%dat(1)
    ! diag_data_device%scalarInfilArea(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarInfilArea)%dat(1)
    ! diag_data_device%scalarFrozenArea(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFrozenArea)%dat(1)
    ! diag_data_device%scalarSoilControl(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSoilControl)%dat(1)
    ! diag_data_device%scalarSoilCompress(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSoilCompress)%dat(1)
    ! diag_data_device%scalarTotalSoilLiq(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilLiq)%dat(1)
    ! diag_data_device%scalarTotalSoilIce(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilIce)%dat(1)
    ! diag_data_device%scalarTotalSoilWat(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilWat)%dat(1)
    ! diag_data_device%scalarKappa(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarKappa)%dat(1)
    ! diag_data_device%scalarVolLatHt_fus(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolLatHt_fus)%dat(1)
    ! diag_data_device%balanceCasNrg(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceCasNrg)%dat(1)
    ! diag_data_device%balanceVegNrg(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceVegNrg)%dat(1)
    ! diag_data_device%balanceSnowNrg(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSnowNrg)%dat(1)
    ! diag_data_device%balanceSoilNrg(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSoilNrg)%dat(1)
    ! diag_data_device%balanceVegMass(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceVegMass)%dat(1)
    ! diag_data_device%balanceSnowMass(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSnowMass)%dat(1)
    ! diag_data_device%balanceSoilMass(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSoilMass)%dat(1)
    ! diag_data_device%balanceAqMass(iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceAqMass)%dat(1)
    end do
    do iGRU=1,nGRU
      ! diag_data_device%mLayerEnthTemp(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerEnthTemp)%dat(1:nLayers)
      ! diag_data_device%mLayerVolHtCapBulk_m(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolHtCapBulk)%dat(1:nLayers)
      ! diag_data_device%mLayerCm_m(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerCm)%dat(1:nLayers)
      ! diag_data_device%mLayerThermalC_m(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThermalC)%dat(1:nLayers)
      ! diag_data_device%iLayerThermalC_m(0:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%iLayerThermalC)%dat(0:nLayers)
      ! diag_data_device%spectralAlbGndDirect(1:nBand,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%spectralAlbGndDirect)%dat(1:nBand)
      ! diag_data_device%spectralAlbGndDiffuse(1:nBand,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%spectralAlbGndDiffuse)%dat(1:nBand)
      ! diag_data_device%mLayerTranspireLim_m(1:nSoil,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerTranspireLim)%dat(1:nSoil)
      ! ! diag_data_device%mLayerRootDensity_m(1:nSoil,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerRootDensity)%dat(1:nSoil)
      ! diag_data_device%spectralSnowAlbedoDirect(:,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%spectralSnowAlbedoDirect)%dat
      ! if (nSnow.ne.0) diag_data_device%mLayerFracLiqSnow_m(1:nSnow,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerFracLiqSnow)%dat(1:nSnow)
      ! if (nSnow.ne.0) diag_data_device%mLayerThetaResid_m(1:nSnow,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThetaResid)%dat(1:nSnow)
      ! if (nSnow.ne.0) diag_data_device%mLayerPoreSpace_m(1:nSnow,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerPoreSpace)%dat(1:nSnow)
      ! diag_data_device%mLayerMeltFreeze(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerMeltFreeze)%dat(1:nLayers)
      ! diag_data_device%mLayerVolFracAir_m(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolFracAir)%dat(1:nLayers)
      ! diag_data_device%mLayerTcrit(1:nSoil,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerTcrit)%dat(1:nSoil)
      ! diag_data_device%mLayerCompress_m(1:nSoil,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerCompress)%dat(1:nSoil)
      ! diag_data_device%mLayerMatricHeadLiq(1:nSoil,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerMatricHeadLiq)%dat(1:nSoil)
      ! ! diag_data_device%scalarVGn_m_m(1:nSoil,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVGn_m)%dat(1:nSoil)
      ! diag_data_device%balanceLayerNrg(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerNrg)%dat(1:nLayers)
      ! diag_data_device%balanceLayerMass(1:nLayers,iGRU) = diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerMass)%dat(1:nLayers)
  
    end do
  end subroutine allocate_device_diag_data
      
  subroutine allocate_device_diag_temp(diag_data_device,nGRU,nSoil)
    use globaldata,only:maxSnowLayers
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_device
    integer(i4b),intent(in) :: nGRU,nSoil

    allocate(diag_data_device%scalarCanopyDepth(nGRU))
    allocate(diag_data_device%scalarGreenVegFraction(nGRU))
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
    allocate(diag_data_device%mLayerVolHtCapBulk_m(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%mLayerCm_m(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%scalarLambda_drysoil(nGRU))
    allocate(diag_data_device%scalarLambda_wetsoil(nGRU))
    allocate(diag_data_device%mLayerThermalC_m(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%iLayerThermalC_m(0:nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%scalarCanopyEnthTemp(nGRU))
    allocate(diag_data_device%mLayerEnthTemp(nSoil+maxSnowLayers,nGRU))
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
    allocate(diag_data_device%spectralAlbGndDirect(nBand,nGRU))
    allocate(diag_data_device%spectralAlbGndDiffuse(nBand,nGRU))
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
    allocate(diag_data_device%mLayerTranspireLim_m(nSoil,nGRU))
    allocate(diag_data_device%mLayerRootDensity_m(nSoil,nGRU))
    allocate(diag_data_device%scalarAquiferRootFrac(nGRU))
    allocate(diag_data_device%scalarFracLiqVeg(nGRU))
    allocate(diag_data_device%scalarCanopyWetFraction(nGRU))
    allocate(diag_data_device%scalarSnowAge(nGRU))
    allocate(diag_data_device%scalarGroundSnowFraction(nGRU))
    allocate(diag_data_device%spectralSnowAlbedoDirect(nBand,nGRU))
    allocate(diag_data_device%mLayerFracLiqSnow_m(maxSnowLayers,nGRU))
    allocate(diag_data_device%mLayerThetaResid_m(maxSnowLayers,nGRU))
    allocate(diag_data_device%mLayerPoreSpace_m(maxSnowLayers,nGRU))
    allocate(diag_data_device%mLayerMeltFreeze(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%scalarInfilArea(nGRU))
    allocate(diag_data_device%scalarFrozenArea(nGRU))
    allocate(diag_data_device%scalarSoilControl(nGRU))
    allocate(diag_data_device%mLayerVolFracAir_m(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%mLayerTcrit(nSoil,nGRU))
    allocate(diag_data_device%mLayerCompress_m(nSoil,nGRU))
    allocate(diag_data_device%scalarSoilCompress(nGRU))
    allocate(diag_data_device%mLayerMatricHeadLiq(nSoil,nGRU))
    allocate(diag_data_device%scalarTotalSoilLiq(nGRU))
    allocate(diag_data_device%scalarTotalSoilIce(nGRU))
    allocate(diag_data_device%scalarTotalSoilWat(nGRU))
    allocate(diag_data_device%scalarVGn_m_m(nSoil,nGRU))
    allocate(diag_data_device%scalarKappa(nGRU))
    allocate(diag_data_device%scalarVolLatHt_fus(nGRU))
    allocate(diag_data_device%balanceCasNrg(nGRU))
    allocate(diag_data_device%balanceVegNrg(nGRU))
    allocate(diag_data_device%balanceLayerNrg(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%balanceSnowNrg(nGRU))
    allocate(diag_data_device%balanceSoilNrg(nGRU))
    allocate(diag_data_device%balanceVegMass(nGRU))
    allocate(diag_data_device%balanceLayerMass(nSoil+maxSnowLayers,nGRU))
    allocate(diag_data_device%balanceSnowMass(nGRU))
    allocate(diag_data_device%balanceSoilMass(nGRU))
    allocate(diag_data_device%balanceAqMass(nGRU))

  end subroutine allocate_device_diag_temp

  subroutine copy_device_diag(diag_data_out,diag_data_in)
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_in,diag_data_out

    diag_data_out%scalarCanopyDepth = diag_data_in%scalarCanopyDepth
    diag_data_out%scalarGreenVegFraction = diag_data_in%scalarGreenVegFraction
    diag_data_out%scalarBulkVolHeatCapVeg = diag_data_in%scalarBulkVolHeatCapVeg
    diag_data_out%scalarCanopyCm = diag_data_in%scalarCanopyCm
    diag_data_out%scalarCanopyEmissivity = diag_data_in%scalarCanopyEmissivity
    diag_data_out%scalarRootZoneTemp = diag_data_in%scalarRootZoneTemp
    diag_data_out%scalarLAI = diag_data_in%scalarLAI
    diag_data_out%scalarSAI = diag_data_in%scalarSAI
    diag_data_out%scalarExposedLAI = diag_data_in%scalarExposedLAI
    diag_data_out%scalarExposedSAI = diag_data_in%scalarExposedSAI
    diag_data_out%scalarAdjMeasHeight = diag_data_in%scalarAdjMeasHeight
    diag_data_out%scalarCanopyIceMax = diag_data_in%scalarCanopyIceMax
    diag_data_out%scalarCanopyLiqMax = diag_data_in%scalarCanopyLiqMax
    diag_data_out%scalarGrowingSeasonIndex = diag_data_in%scalarGrowingSeasonIndex
    diag_data_out%scalarVolHtCap_air = diag_data_in%scalarVolHtCap_air
    diag_data_out%scalarVolHtCap_ice = diag_data_in%scalarVolHtCap_ice
    diag_data_out%scalarVolHtCap_soil = diag_data_in%scalarVolHtCap_soil
    diag_data_out%scalarVolHtCap_water = diag_data_in%scalarVolHtCap_water
    diag_data_out%mLayerVolHtCapBulk_m = diag_data_in%mLayerVolHtCapBulk_m
    diag_data_out%mLayerCm_m = diag_data_in%mLayerCm_m
    diag_data_out%scalarLambda_drysoil = diag_data_in%scalarLambda_drysoil
    diag_data_out%scalarLambda_wetsoil = diag_data_in%scalarLambda_wetsoil
    diag_data_out%mLayerThermalC_m = diag_data_in%mLayerThermalC_m
    diag_data_out%iLayerThermalC_m = diag_data_in%iLayerThermalC_m
    diag_data_out%scalarCanopyEnthTemp = diag_data_in%scalarCanopyEnthTemp
    diag_data_out%mLayerEnthTemp = diag_data_in%mLayerEnthTemp
    diag_data_out%scalarTotalSoilEnthalpy = diag_data_in%scalarTotalSoilEnthalpy
    diag_data_out%scalarTotalSnowEnthalpy = diag_data_in%scalarTotalSnowEnthalpy
    diag_data_out%scalarVPair = diag_data_in%scalarVPair
    diag_data_out%scalarVP_CanopyAir = diag_data_in%scalarVP_CanopyAir
    diag_data_out%scalarTwetbulb = diag_data_in%scalarTwetbulb
    diag_data_out%scalarSnowfallTemp = diag_data_in%scalarSnowfallTemp
    diag_data_out%scalarNewSnowDensity = diag_data_in%scalarNewSnowDensity
    diag_data_out%scalarO2air = diag_data_in%scalarO2air
    diag_data_out%scalarCO2air = diag_data_in%scalarCO2air
    diag_data_out%windspd_x = diag_data_in%windspd_x
    diag_data_out%windspd_y = diag_data_in%windspd_y
    diag_data_out%scalarCosZenith = diag_data_in%scalarCosZenith
    diag_data_out%scalarFractionDirect = diag_data_in%scalarFractionDirect
    diag_data_out%scalarCanopySunlitFraction = diag_data_in%scalarCanopySunlitFraction
    diag_data_out%scalarCanopySunlitLAI = diag_data_in%scalarCanopySunlitLAI
    diag_data_out%scalarCanopyShadedLAI = diag_data_in%scalarCanopyShadedLAI
    diag_data_out%spectralAlbGndDirect = diag_data_in%spectralAlbGndDirect
    diag_data_out%spectralAlbGndDiffuse = diag_data_in%spectralAlbGndDiffuse
    diag_data_out%scalarGroundAlbedo = diag_data_in%scalarGroundAlbedo
    diag_data_out%scalarLatHeatSubVapCanopy = diag_data_in%scalarLatHeatSubVapCanopy
    diag_data_out%scalarLatHeatSubVapGround = diag_data_in%scalarLatHeatSubVapGround
    diag_data_out%scalarSatVP_CanopyTemp = diag_data_in%scalarSatVP_CanopyTemp
    diag_data_out%scalarSatVP_GroundTemp = diag_data_in%scalarSatVP_GroundTemp
    diag_data_out%scalarZ0Canopy = diag_data_in%scalarZ0Canopy
    diag_data_out%scalarWindReductionFactor = diag_data_in%scalarWindReductionFactor
    diag_data_out%scalarZeroPlaneDisplacement = diag_data_in%scalarZeroPlaneDisplacement
    diag_data_out%scalarRiBulkCanopy = diag_data_in%scalarRiBulkCanopy
    diag_data_out%scalarRiBulkGround = diag_data_in%scalarRiBulkGround
    diag_data_out%scalarCanopyStabilityCorrection = diag_data_in%scalarCanopyStabilityCorrection
    diag_data_out%scalarGroundStabilityCorrection = diag_data_in%scalarGroundStabilityCorrection
    diag_data_out%scalarIntercellularCO2Sunlit = diag_data_in%scalarIntercellularCO2Sunlit
    diag_data_out%scalarIntercellularCO2Shaded = diag_data_in%scalarIntercellularCO2Shaded
    diag_data_out%scalarTranspireLim = diag_data_in%scalarTranspireLim
    diag_data_out%scalarTranspireLimAqfr = diag_data_in%scalarTranspireLimAqfr
    diag_data_out%scalarFoliageNitrogenFactor = diag_data_in%scalarFoliageNitrogenFactor
    diag_data_out%scalarSoilRelHumidity = diag_data_in%scalarSoilRelHumidity
    diag_data_out%mLayerTranspireLim_m = diag_data_in%mLayerTranspireLim_m
    diag_data_out%mLayerRootDensity_m = diag_data_in%mLayerRootDensity_m
    diag_data_out%scalarAquiferRootFrac = diag_data_in%scalarAquiferRootFrac
    diag_data_out%scalarFracLiqVeg = diag_data_in%scalarFracLiqVeg
    diag_data_out%scalarCanopyWetFraction = diag_data_in%scalarCanopyWetFraction
    diag_data_out%scalarSnowAge = diag_data_in%scalarSnowAge
    diag_data_out%scalarGroundSnowFraction = diag_data_in%scalarGroundSnowFraction
    diag_data_out%spectralSnowAlbedoDirect = diag_data_in%spectralSnowAlbedoDirect
    diag_data_out%mLayerFracLiqSnow_m = diag_data_in%mLayerFracLiqSnow_m
    diag_data_out%mLayerThetaResid_m = diag_data_in%mLayerThetaResid_m
    diag_data_out%mLayerPoreSpace_m = diag_data_in%mLayerPoreSpace_m
    diag_data_out%mLayerMeltFreeze = diag_data_in%mLayerMeltFreeze
    diag_data_out%scalarInfilArea = diag_data_in%scalarInfilArea
    diag_data_out%scalarFrozenArea = diag_data_in%scalarFrozenArea
    diag_data_out%scalarSoilControl = diag_data_in%scalarSoilControl
    diag_data_out%mLayerVolFracAir_m = diag_data_in%mLayerVolFracAir_m
    diag_data_out%mLayerTcrit = diag_data_in%mLayerTcrit
    diag_data_out%mLayerCompress_m = diag_data_in%mLayerCompress_m
    diag_data_out%scalarSoilCompress = diag_data_in%scalarSoilCompress
    diag_data_out%mLayerMatricHeadLiq = diag_data_in%mLayerMatricHeadLiq
    diag_data_out%scalarTotalSoilLiq = diag_data_in%scalarTotalSoilLiq
    diag_data_out%scalarTotalSoilIce = diag_data_in%scalarTotalSoilIce
    diag_data_out%scalarTotalSoilWat = diag_data_in%scalarTotalSoilWat
    diag_data_out%scalarVGn_m_m = diag_data_in%scalarVGn_m_m
    diag_data_out%scalarKappa = diag_data_in%scalarKappa
    diag_data_out%scalarVolLatHt_fus = diag_data_in%scalarVolLatHt_fus
    diag_data_out%balanceCasNrg = diag_data_in%balanceCasNrg
    diag_data_out%balanceVegNrg = diag_data_in%balanceVegNrg
    diag_data_out%balanceLayerNrg = diag_data_in%balanceLayerNrg
    diag_data_out%balanceSnowNrg = diag_data_in%balanceSnowNrg
    diag_data_out%balanceSoilNrg = diag_data_in%balanceSoilNrg
    diag_data_out%balanceVegMass = diag_data_in%balanceVegMass
    diag_data_out%balanceLayerMass = diag_data_in%balanceLayerMass
    diag_data_out%balanceSnowMass = diag_data_in%balanceSnowMass
    diag_data_out%balanceSoilMass = diag_data_in%balanceSoilMass
    diag_data_out%balanceAqMass = diag_data_in%balanceAqMass

  end subroutine copy_device_diag


  subroutine finalize_device_diag_data(diag_data_device,diag_data,nSnow,nSoil,nLayers,nGRU)
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_device
    type(gru_hru_doubleVec),intent(inout) :: diag_data
    integer(i4b),intent(in) :: nSnow,nSoil,nLayers,nGRU
    integer(i4b) :: iGRU

    do iGRU=1,nGRU
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%numFluxCalls)%dat(1) = diag_data_device%numFluxCalls
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%wallClockTime)%dat(1) = diag_data_device%wallClockTime
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%meanStepSize)%dat(1) = diag_data_device%meanStepSize

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

    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyDepth)%dat(1) = diag_data_device%scalarCanopyDepth(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGreenVegFraction)%dat(1) = diag_data_device%scalarGreenVegFraction(iGRU)
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
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLambda_drysoil)%dat(1) = diag_data_device%scalarLambda_drysoil(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarLambda_wetsoil)%dat(1) = diag_data_device%scalarLambda_wetsoil(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1) = diag_data_device%scalarCanopyEnthTemp(iGRU)
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
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarAquiferRootFrac)%dat(1) = diag_data_device%scalarAquiferRootFrac(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFracLiqVeg)%dat(1) = diag_data_device%scalarFracLiqVeg(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarCanopyWetFraction)%dat(1) = diag_data_device%scalarCanopyWetFraction(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSnowAge)%dat(1) = diag_data_device%scalarSnowAge(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarGroundSnowFraction)%dat(1) = diag_data_device%scalarGroundSnowFraction(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarInfilArea)%dat(1) = diag_data_device%scalarInfilArea(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarFrozenArea)%dat(1) = diag_data_device%scalarFrozenArea(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSoilControl)%dat(1) = diag_data_device%scalarSoilControl(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarSoilCompress)%dat(1) = diag_data_device%scalarSoilCompress(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilLiq)%dat(1) = diag_data_device%scalarTotalSoilLiq(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilIce)%dat(1) = diag_data_device%scalarTotalSoilIce(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarTotalSoilWat)%dat(1) = diag_data_device%scalarTotalSoilWat(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarKappa)%dat(1) = diag_data_device%scalarKappa(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVolLatHt_fus)%dat(1) = diag_data_device%scalarVolLatHt_fus(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceCasNrg)%dat(1) = diag_data_device%balanceCasNrg(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceVegNrg)%dat(1) = diag_data_device%balanceVegNrg(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSnowNrg)%dat(1) = diag_data_device%balanceSnowNrg(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSoilNrg)%dat(1) = diag_data_device%balanceSoilNrg(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceVegMass)%dat(1) = diag_data_device%balanceVegMass(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSnowMass)%dat(1) = diag_data_device%balanceSnowMass(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceSoilMass)%dat(1) = diag_data_device%balanceSoilMass(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceAqMass)%dat(1) = diag_data_device%balanceAqMass(iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerEnthTemp)%dat = diag_data_device%mLayerEnthTemp(:,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolHtCapBulk)%dat(1:nLayers) = diag_data_device%mLayerVolHtCapBulk_m(1:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerCm)%dat(1:nLayers) = diag_data_device%mLayerCm_m(1:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThermalC)%dat(1:nLayers) = diag_data_device%mLayerThermalC_m(1:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%iLayerThermalC)%dat(0:nLayers) = diag_data_device%iLayerThermalC_m(0:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%spectralAlbGndDirect)%dat(1:nBand) = diag_data_device%spectralAlbGndDirect(1:nBand,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%spectralAlbGndDiffuse)%dat(1:nBand) = diag_data_device%spectralAlbGndDiffuse(1:nBand,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerTranspireLim)%dat(1:nSoil) = diag_data_device%mLayerTranspireLim_m(1:nSoil,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerRootDensity)%dat(1:nSoil) = diag_data_device%mLayerRootDensity_m(1:nSoil,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%spectralSnowAlbedoDirect)%dat(1:nBand) = diag_data_device%spectralSnowAlbedoDirect(1:nBand,iGRU)
    if (nSnow.ne.0) diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerFracLiqSnow)%dat(1:nSnow) = diag_data_device%mLayerFracLiqSnow_m(1:nSnow,iGRU)
    if (nSnow.ne.0) diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerThetaResid)%dat(1:nSnow) = diag_data_device%mLayerThetaResid_m(1:nSnow,iGRU)
    if (nSnow.ne.0) diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerPoreSpace)%dat(1:nSnow) = diag_data_device%mLayerPoreSpace_m(1:nSnow,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerMeltFreeze)%dat(1:nLayers) = diag_data_device%mLayerMeltFreeze(1:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerVolFracAir)%dat(1:nLayers) = diag_data_device%mLayerVolFracAir_m(1:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerTcrit)%dat(1:nSoil) = diag_data_device%mLayerTcrit(1:nSoil,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerCompress)%dat(1:nSoil) = diag_data_device%mLayerCompress_m(1:nSoil,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%mLayerMatricHeadLiq)%dat(1:nSoil) = diag_data_device%mLayerMatricHeadLiq(1:nSoil,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%scalarVGn_m)%dat(1:nSoil) = diag_data_device%scalarVGn_m_m(1:nSoil,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerNrg)%dat(1:nLayers) = diag_data_device%balanceLayerNrg(1:nLayers,iGRU)
    diag_data%gru(iGRU)%hru(1)%var(iLookDIAG%balanceLayerMass)%dat(1:nLayers)   = diag_data_device%balanceLayerMass(1:nLayers,iGRU)
    end do
  end subroutine finalize_device_diag_data

  subroutine deallocate_device_diag_data(diag_data_device)
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_device
    deallocate(diag_data_device%scalarCanopyDepth)
    deallocate(diag_data_device%scalarGreenVegFraction)
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
    deallocate(diag_data_device%mLayerVolHtCapBulk_m)
    deallocate(diag_data_device%mLayerCm_m)
    deallocate(diag_data_device%scalarLambda_drysoil)
    deallocate(diag_data_device%scalarLambda_wetsoil)
    deallocate(diag_data_device%mLayerThermalC_m)
    deallocate(diag_data_device%iLayerThermalC_m)
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
    deallocate(diag_data_device%mLayerTranspireLim_m)
    deallocate(diag_data_device%mLayerRootDensity_m)
    deallocate(diag_data_device%scalarAquiferRootFrac)
    deallocate(diag_data_device%scalarFracLiqVeg)
    deallocate(diag_data_device%scalarCanopyWetFraction)
    deallocate(diag_data_device%scalarSnowAge)
    deallocate(diag_data_device%scalarGroundSnowFraction)
    deallocate(diag_data_device%spectralSnowAlbedoDirect)
    deallocate(diag_data_device%mLayerFracLiqSnow_m)
    deallocate(diag_data_device%mLayerThetaResid_m)
    deallocate(diag_data_device%mLayerPoreSpace_m)
    deallocate(diag_data_device%mLayerMeltFreeze)
    deallocate(diag_data_device%scalarInfilArea)
    deallocate(diag_data_device%scalarFrozenArea)
    deallocate(diag_data_device%scalarSoilControl)
    deallocate(diag_data_device%mLayerVolFracAir_m)
    deallocate(diag_data_device%mLayerTcrit)
    deallocate(diag_data_device%mLayerCompress_m)
    deallocate(diag_data_device%scalarSoilCompress)
    deallocate(diag_data_device%mLayerMatricHeadLiq)
    deallocate(diag_data_device%scalarTotalSoilLiq)
    deallocate(diag_data_device%scalarTotalSoilIce)
    deallocate(diag_data_device%scalarTotalSoilWat)
    deallocate(diag_data_device%scalarVGn_m_m)
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
        integer(i4b) :: nGRU
        integer(i4b) :: iGRU
        
        allocate(bvar_data_device%basin__SurfaceRunoff(nGRU))
        allocate(bvar_data_device%basin__ColumnOutflow(nGRU))
        allocate(bvar_data_device%basin__AquiferStorage(nGRU))
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
  bvar_data_device%routingRunoffFuture(:,iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%routingRunoffFuture)%dat
bvar_data_device%routingFractionFuture(:,iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%routingFractionFuture)%dat
        bvar_data_device%basin__SurfaceRunoff(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)
        bvar_data_device%basin__ColumnOutflow(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1)
        ! bvar_data_device%basin__AquiferStorage(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferStorage)%dat(1)
        bvar_data_device%basin__AquiferRecharge(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferRecharge)%dat(1)
        bvar_data_device%basin__AquiferBaseflow(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)
        bvar_data_device%basin__AquiferTranspire(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferTranspire)%dat(1)
        bvar_data_device%basin__TotalRunoff(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__TotalRunoff)%dat(1)
        bvar_data_device%basin__SoilDrainage(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__SoilDrainage)%dat(1)
        bvar_data_device%basin__totalArea(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)
bvar_data_device%averageInstantRunoff(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%averageInstantRunoff)%dat(1)
bvar_data_device%averageRoutedRunoff(iGRU) = bvar_data%gru(iGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(1)
end do
      end subroutine allocate_device_bvar_data

      subroutine finalize_device_bvar_data(bvar_data_device,bvar_data,nGRU)
        implicit none
        type(bvar_data_device),intent(inout) :: bvar_data_device
        type(gru_doubleVec),intent(inout) :: bvar_data
        integer(i4b) :: nGRU
        integer(i4b) :: iGRU
        
        do iGRU=1,nGRU
        bvar_data%gru(iGRU)%var(iLookBVAR%basin__SurfaceRunoff)%dat(1) = bvar_data_device%basin__SurfaceRunoff(iGRU)
        bvar_data%gru(iGRU)%var(iLookBVAR%basin__ColumnOutflow)%dat(1) = bvar_data_device%basin__ColumnOutflow(iGRU)
        bvar_data%gru(iGRU)%var(iLookBVAR%basin__AquiferStorage)%dat(1) = bvar_data_device%basin__AquiferStorage(iGRU)
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
      end subroutine
      subroutine deallocate_device_bvar_data(bvar_data_device)
        type(bvar_data_device),intent(inout) :: bvar_data_device
        
        deallocate(bvar_data_device%basin__AquiferStorage)
      end subroutine deallocate_device_bvar_data

      subroutine allocate_device_forc_data(forc_data_device, forc_data,nGRU)
        type(forc_data_device),intent(inout) :: forc_data_device
        type(gru_hru_double),intent(in) :: forc_data
        integer(i4b) :: nGRU
        integer(i4b) :: iGRU

        allocate(forc_data_device%LWRadAtm_d(nGRU))
        allocate(forc_data_device%airpres_d(nGRU))
        allocate(forc_data_device%airtemp_d(nGRU))
        allocate(forc_data_device%windspd_d(nGRU))
          allocate(forc_data_device%time(nGRU))
  allocate(forc_data_device%pptrate(nGRU))
  allocate(forc_data_device%SWRadAtm(nGRU))
  allocate(forc_data_device%spechum(nGRU))
  do iGRU=1,nGRU
          forc_data_device%LWRadAtm_d(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%LWRadAtm)
        forc_data_device%airpres_d(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%airpres)
        forc_data_device%airtemp_d(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%airtemp)
        forc_data_device%windspd_d(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%windspd)
          forc_data_device%time(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%time)
  forc_data_device%pptrate(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%pptrate)
  forc_data_device%SWRadAtm(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%SWRadAtm)
  forc_data_device%spechum(iGRU) = forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%spechum)
  end do
      end subroutine allocate_device_forc_data
    
      subroutine finalize_device_forc_data(forc_data_device,forc_data,nGRU)
                type(forc_data_device),intent(inout) :: forc_data_device
        type(gru_hru_double),intent(inout) :: forc_data
        integer(i4b) :: nGRU
        integer(i4b) :: iGRU

        do iGRU=1,nGRU
        forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%LWRadAtm) = forc_data_device%LWRadAtm_d(iGRU)
        forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%airpres) = forc_data_device%airpres_d(iGRU)
        forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%airtemp) = forc_data_device%airtemp_d(iGRU)
        forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%windspd) = forc_data_device%windspd_d(iGRU)
        forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%time) = forc_data_device%time(iGRU)
        forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%pptrate) = forc_data_device%pptrate(iGRU)
        forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%SWRadAtm) = forc_data_device%SWRadAtm(iGRU)
        forc_data%gru(iGRU)%hru(1)%var(iLookFORCE%spechum) = forc_data_device%spechum(iGRU)
        end do
      end subroutine
      subroutine deallocate_device_forc_data(forc_data_device)
        type(forc_data_device),intent(inout) :: forc_data_device

        deallocate(forc_data_device%LWRadAtm_d)
        deallocate(forc_data_device%airpres_d)
        deallocate(forc_data_device%airtemp_d)
        deallocate(forc_data_device%windspd_d)
      end subroutine deallocate_device_forc_data
    
      subroutine allocate_device_decisions(device_decisions)
        USE globalData,only:model_decisions ! model decision structure
        USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
        use globalData,only:urbanVegCategory
        use globalData,only:data_step
        use globalData,only:refJulDay
        use enthalpyTemp_module,only:h_lookup,t_lookup
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
        device_decisions%urbanVegCategory = urbanVegCategory

        device_decisions%data_step = data_step
device_decisions%refJulDay = refJulDay
device_decisions%h_lookup = h_lookup
device_decisions%t_lookup = t_lookup
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

      subroutine allocate_device_type_data(type_data_d, type_data, nGRU)
        type(type_data_device) :: type_data_d
        type(gru_hru_int) :: type_data
        integer(i4b) :: nGRU
        integer(i4b) :: iGRU
        allocate(type_data_d%vegTypeIndex(nGRU))
        allocate(type_data_d%soilTypeIndex(nGRU))
        allocate(type_data_d%slopeTypeIndex(nGRU))
        allocate(type_data_d%downHRUindex(nGRU))

        do iGRU=1,nGRU
        type_data_d%vegTypeIndex(iGRU) = type_data%gru(iGRU)%hru(1)%var(iLookTYPE%vegTypeIndex)
        type_data_d%soilTypeIndex(iGRU) = type_data%gru(iGRU)%hru(1)%var(iLookTYPE%soilTypeIndex)
        type_data_d%slopeTypeIndex(iGRU) = type_data%gru(iGRU)%hru(1)%var(iLookTYPE%slopeTypeIndex)
        type_data_d%downHRUindex(iGRU) = type_data%gru(iGRU)%hru(1)%var(iLookTYPE%downHRUindex)
        end do
      end subroutine
      subroutine deallocate_device_type_data(type_data_d)
        type(type_data_device) :: type_data_d
        deallocate(type_data_d%vegTypeIndex)
        deallocate(type_data_d%soilTypeIndex)
        deallocate(type_data_d%slopeTypeIndex)
        deallocate(type_data_d%downHRUindex)

      end subroutine

      subroutine allocate_veg_parameters(veg_param,nGRU)
        use NOAHMP_VEG_PARAMETERS
        use noahmp_globals
        use NOAHMP_RAD_PARAMETERS
        implicit none
        type(veg_parameters) :: veg_param
        integer(i4b) :: nGRU

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
        veg_param%rgl_d = rgl
        allocate(veg_param%rsmin_d(nGRU))
        veg_param%rsmin_d = rsmin
        allocate(veg_param%rsmax)
        veg_param%rsmax = rsmax
        allocate(veg_param%topt)
        veg_param%topt = topt
        allocate(veg_param%hs_d(nGRU))
        veg_param%hs_d = hs
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
veg_param%hvt = hvt
veg_param%hvb = hvb
veg_param%rc = rc
veg_param%swemx = swemx
veg_param%albsat = albsat
veg_param%albdry = albdry
veg_param%alblak = alblak
veg_param%isbarren = isbarren
veg_param%ISSNOW = ISSNOW
veg_param%ISWATER = ISWATER
veg_param%saim = saim
veg_param%laim = laim
veg_param%dveg = dveg
veg_param%tmin = tmin

      end subroutine

      subroutine allocate_veg_param_tables(tables)
        use module_sf_noahlsm
        implicit none
        type(veg_param_tables) :: tables

        tables%hstbl = hstbl
tables%rstbl = rstbl
tables%rgltbl = rgltbl
      end subroutine
        
  subroutine allocate_device_attr_data(attr_data_device,attr_data,nGRU)
    type(attr_data_device) :: attr_data_device
    type(gru_hru_double) :: attr_data
    integer(i4b) :: nGRU
    integer(i4b) :: iGRU
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
    end module initialize_device      