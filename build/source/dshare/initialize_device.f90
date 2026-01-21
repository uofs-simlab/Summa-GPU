module initialize_device
    use data_types
    use device_data_types
    USE var_lookup,only:iLookFLUX        ! lookup indices for flux data
    USE var_lookup,only:iLookDERIV       ! lookup indices for derivative data
    USE var_lookup,only:iLookFORCE       ! lookup indices for forcing data 
    USE var_lookup,only:iLookDIAG        ! lookup indices for diagnostic variable data
    USE var_lookup,only:iLookPROG,iLookPARAM, iLookINDEX, iLookBVAR, iLookTYPE, iLookATTR
    ! USE globalData,only: nSpecBands                 ! number of spectral bands

    integer(i4b),parameter :: nSpecBands = 2
    contains

    subroutine allocate_device_param_data(mpar_data_device,mpar_data,nGRU)
        type(mpar_data_device), intent(inout) :: mpar_data_device
        type(var_dlength),intent(in) :: mpar_data
        integer(i4b) :: nGRU
        mpar_data_device%Fcapil = mpar_data%var(iLookPARAM%Fcapil)%dat(1)
        mpar_data_device%Louis79_bparam = mpar_data%var(iLookPARAM%Louis79_bparam)%dat(1)
        mpar_data_device%Louis79_cStar = mpar_data%var(iLookPARAM%Louis79_cStar)%dat(1)
        mpar_data_device%Mahrt87_eScale = mpar_data%var(iLookPARAM%Mahrt87_eScale)%dat(1)
        allocate(mpar_data_device%aquiferBaseflowExp_(nGRU))
        mpar_data_device%aquiferBaseflowExp_ = mpar_data%var(iLookPARAM%aquiferBaseflowExp)%dat(1)
        allocate(mpar_data_device%aquiferBaseflowRate_(nGRU))
        mpar_data_device%aquiferBaseflowRate_ = mpar_data%var(iLookPARAM%aquiferBaseflowRate)%dat(1)
        allocate(mpar_data_device%aquiferScaleFactor_(nGRU))
        mpar_data_device%aquiferScaleFactor_ = mpar_data%var(iLookPARAM%aquiferScaleFactor)%dat(1)
        mpar_data_device%canopyDrainageCoeff = mpar_data%var(iLookPARAM%canopyDrainageCoeff)%dat(1)
        allocate(mpar_data_device%canopyWettingExp_(nGRU))
        mpar_data_device%canopyWettingExp_ = mpar_data%var(iLookPARAM%canopyWettingExp)%dat(1)
        allocate(mpar_data_device%canopyWettingFactor_(nGRU))
        mpar_data_device%canopyWettingFactor_ = mpar_data%var(iLookPARAM%canopyWettingFactor)%dat(1)
        mpar_data_device%critAquiferTranspire = mpar_data%var(iLookPARAM%critAquiferTranspire)%dat(1)
        mpar_data_device%critRichNumber = mpar_data%var(iLookPARAM%critRichNumber)%dat(1)
        mpar_data_device%critSoilTranspire = mpar_data%var(iLookPARAM%critSoilTranspire)%dat(1)
        mpar_data_device%critSoilWilting = mpar_data%var(iLookPARAM%critSoilWilting)%dat(1)
        mpar_data_device%f_impede = mpar_data%var(iLookPARAM%f_impede)%dat(1)
        mpar_data_device%fixedThermalCond_snow = mpar_data%var(iLookPARAM%fixedThermalCond_snow)%dat(1)
        mpar_data_device%frac_clay = mpar_data%var(iLookPARAM%frac_clay)%dat
        mpar_data_device%frac_sand = mpar_data%var(iLookPARAM%frac_sand)%dat
        mpar_data_device%frac_silt = mpar_data%var(iLookPARAM%frac_silt)%dat
        allocate(mpar_data_device%heightCanopyBottom_(nGRU))
        allocate(mpar_data_device%heightCanopyTop_(nGRU))
        mpar_data_device%heightCanopyBottom_ = mpar_data%var(iLookPARAM%heightCanopyBottom)%dat(1)
        mpar_data_device%heightCanopyTop_ = mpar_data%var(iLookPARAM%heightCanopyTop)%dat(1)
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
        mpar_data_device%minstep = mpar_data%var(iLookPARAM%minstep)%dat(1)
        mpar_data_device%maxstep = mpar_data%var(iLookPARAM%maxstep)%dat(1)
        mpar_data_device%refInterceptCapRain = mpar_data%var(iLookPARAM%refInterceptCapRain)%dat(1)
        mpar_data_device%refInterceptCapSnow = mpar_data%var(iLookPARAM%refInterceptCapSnow)%dat(1)
        mpar_data_device%Frad_vis = mpar_data%var(iLookPARAM%Frad_vis)%dat(1)
        mpar_data_device%Frad_direct = mpar_data%var(iLookPARAM%Frad_direct)%dat(1)
        allocate(mpar_data_device%albedoMax_(nGRU))
        mpar_data_device%albedoMax_ = mpar_data%var(iLookPARAM%albedoMax)%dat(1)
        allocate(mpar_data_device%albedoMinWinter_(nGRU))
        mpar_data_device%albedoMinWinter_ = mpar_data%var(iLookPARAM%albedoMinWinter)%dat(1)
        mpar_data_device%albedoMinSpring = mpar_data%var(iLookPARAM%albedoMinSpring)%dat(1)
        allocate(mpar_data_device%albedoMaxVisible_(nGRU))
        mpar_data_device%albedoMaxVisible_ = mpar_data%var(iLookPARAM%albedoMaxVisible)%dat(1)
        allocate(mpar_data_device%albedoMinVisible_(nGRU))
        mpar_data_device%albedoMinVisible_ = mpar_data%var(iLookPARAM%albedoMinVisible)%dat(1)
        allocate(mpar_data_device%albedoMaxNearIR_(nGRU))
        mpar_data_device%albedoMaxNearIR_ = mpar_data%var(iLookPARAM%albedoMaxNearIR)%dat(1)
        allocate(mpar_data_device%albedoMinNearIR_(nGRU))
        mpar_data_device%albedoMinNearIR_ = mpar_data%var(iLookPARAM%albedoMinNearIR)%dat(1)
        mpar_data_device%albedoDecayRate = mpar_data%var(iLookPARAM%albedoDecayRate)%dat(1)
        mpar_data_device%tempScalGrowth = mpar_data%var(iLookPARAM%tempScalGrowth)%dat(1)
        mpar_data_device%albedoSootLoad = mpar_data%var(iLookPARAM%albedoSootLoad)%dat(1)
        mpar_data_device%albedoRefresh = mpar_data%var(iLookPARAM%albedoRefresh)%dat(1)
        mpar_data_device%snowUnloadingCoeff = mpar_data%var(iLookPARAM%snowUnloadingCoeff)%dat(1)
        mpar_data_device%ratioDrip2Unloading = mpar_data%var(iLookPARAM%ratioDrip2Unloading)%dat(1)
        mpar_data_device%minTempUnloading = mpar_data%var(iLookPARAM%minTempUnloading)%dat(1)
        mpar_data_device%rateTempUnloading = mpar_data%var(iLookPARAM%rateTempUnloading)%dat(1)
        mpar_data_device%minWindUnloading = mpar_data%var(iLookPARAM%minWindUnloading)%dat(1)
        mpar_data_device%rateWindUnloading = mpar_data%var(iLookPARAM%rateWindUnloading)%dat(1)
        mpar_data_device%densScalGrowth = mpar_data%var(iLookPARAM%densScalGrowth)%dat(1)
        mpar_data_device%grainGrowthRate = mpar_data%var(iLookPARAM%grainGrowthRate)%dat(1)
        mpar_data_device%densScalOvrbdn = mpar_data%var(iLookPARAM%densScalOvrbdn)%dat(1)
        mpar_data_device%tempScalOvrbdn = mpar_data%var(iLookPARAM%tempScalOvrbdn)%dat(1)
        mpar_data_device%baseViscosity = mpar_data%var(iLookPARAM%baseViscosity)%dat(1)

      end subroutine allocate_device_param_data

      subroutine deallocate_device_param_data(mpar_data_device)
        type(mpar_data_device), intent(inout) :: mpar_data_device
        deallocate(mpar_data_device%Fcapil)
        deallocate(mpar_data_device%Louis79_bparam)
        deallocate(mpar_data_device%Louis79_cStar)
        deallocate(mpar_data_device%Mahrt87_eScale)
        deallocate(mpar_data_device%aquiferBaseflowExp_)
        deallocate(mpar_data_device%aquiferBaseflowRate_)
        deallocate(mpar_data_device%aquiferScaleFactor_)
        deallocate(mpar_data_device%canopyDrainageCoeff)
        deallocate(mpar_data_device%canopyWettingExp_)
        deallocate(mpar_data_device%canopyWettingFactor_)
        deallocate(mpar_data_device%critAquiferTranspire)
        deallocate(mpar_data_device%critRichNumber)
        deallocate(mpar_data_device%critSoilTranspire)
        deallocate(mpar_data_device%critSoilWilting)
        deallocate(mpar_data_device%f_impede)
        deallocate(mpar_data_device%fixedThermalCond_snow)
        deallocate(mpar_data_device%frac_clay)
        deallocate(mpar_data_device%frac_sand)
        deallocate(mpar_data_device%frac_silt)
        deallocate(mpar_data_device%heightCanopyBottom_)
        deallocate(mpar_data_device%heightCanopyTop_)
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


      subroutine allocate_device_indx_data(indx_data_device, indx_data, nSnow, nSoil, nLayers,nGRU)
        type(indx_data_device), intent(inout) :: indx_data_device
        type(var_ilength),intent(in) :: indx_data              ! indices defining model states and layers
        integer(i4b), intent(in) :: nSnow, nSoil, nLayers,nGRU
        integer(i4b) :: iGRU
      
        allocate(indx_data_device%ixNrgLayer(MAX(1,size(indx_data%var(iLookINDEX%ixNrgLayer)%dat)),nGRU))
        allocate(indx_data_device%ixHydLayer(MAX(1,size(indx_data%var(iLookINDEX%ixHydLayer)%dat)),nGRU))
        indx_data_device%nSoil = nSoil
        allocate(indx_data_device%nLayers_d(nGRU))
        indx_data_device%nLayers_d = nLayers
        allocate(indx_data_device%nSnow(nGRU))
        indx_data_device%nSnow = indx_data%var(iLookINDEX%nSnow)%dat(1)
        ! indx_data_device%maxSnow = nSnow

        ! indx_data_device%ixCasNrg = indx_data%var(iLookINDEX%ixCasNrg)%dat(1)
        ! indx_data_device%ixVegNrg = indx_data%var(iLookINDEX%ixVegNrg)%dat(1)
        ! indx_data_device%ixVegHyd = indx_data%var(iLookINDEX%ixVegHyd)%dat(1)
        ! indx_data_device%ixAqWat = indx_data%var(iLookINDEX%ixAqWat)%dat(1)
        allocate(indx_data_device%ixCasNrg(nGRU))
        indx_data_device%ixCasNrg = indx_data%var(iLookINDEX%ixCasNrg)%dat(1)
        allocate(indx_data_device%ixVegNrg(nGRU))
        indx_data_device%ixVegNrg = indx_data%var(iLookINDEX%ixVegNrg)%dat(1)
        allocate(indx_data_device%ixVegHyd(nGRU))
        indx_data_device%ixVegHyd = indx_data%var(iLookINDEX%ixVegHyd)%dat(1)
        allocate(indx_data_device%ixAqWat(nGRU))
        indx_data_device%ixAqWat = indx_data%var(iLookINDEX%ixAqWat)%dat(1)
        ! indx_data_device%ixTopNrg = indx_data%var(iLookINDEX%ixTopNrg)%dat(1)
        ! indx_data_device%ixTopHyd = indx_data%var(iLookINDEX%ixTopHyd)%dat(1)
        allocate(indx_data_device%ixTopNrg(nGRU))
        indx_data_device%ixTopNrg = indx_data%var(iLookINDEX%ixTopNrg)%dat(1)
        allocate(indx_data_device%ixTopHyd(nGRU))
        indx_data_device%ixTopHyd = indx_data%var(iLookINDEX%ixTopHyd)%dat(1)
        allocate(indx_data_device%ixSnowSoilNrg(MAX(1,size(indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat)),nGRU))
        allocate(indx_data_device%ixSnowSoilHyd(MAX(1,size(indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat)),nGRU))
        allocate(indx_data_device%ixHydType(MAX(1,size(indx_data%var(iLookINDEX%ixHydType)%dat)),nGRU))
        allocate(indx_data_device%ixStateType(MAX(1,size(indx_data%var(iLookINDEX%ixStateType)%dat)),nGRU))
        allocate(indx_data_device%ixHydCanopy(nGRU))
        allocate(indx_data_device%ixNrgCanopy(nGRU))
        ! indx_data_device%ixNrgCanair = indx_data%var(iLookINDEX%ixNrgCanair)%dat(1)
        if (size(indx_data%var(iLookINDEX%ixControlVolume)%dat).ne.0) allocate(indx_data_device%ixControlVolume(size(indx_data%var(iLookINDEX%ixControlVolume)%dat),nGRU))
        allocate(indx_data_device%layerType(MAX(1,size(indx_data%var(iLookINDEX%layerType)%dat)),nGRU))
        allocate(indx_data_device%ixSoilOnlyHyd(MAX(1,size(indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat)),nGRU))
        allocate(indx_data_device%ixSnowOnlyHyd(MAX(1,size(indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat)),nGRU))
        allocate(indx_data_device%ixSoilOnlyNrg(MAX(1,size(indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat)),nGRU))
        allocate(indx_data_device%ixSnowOnlyNrg(MAX(1,size(indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat)),nGRU))
        
        if (size(indx_data%var(iLookINDEX%ixSoilState)%dat).ne.0) allocate(indx_data_device%ixSoilState(size(indx_data%var(iLookINDEX%ixSoilState)%dat),nGRU))

        ! indx_data_device%ixLayerState = indx_data%var(iLookINDEX%ixLayerState)%dat
        indx_data_device%numberFluxCalc = indx_data%var(iLookINDEX%numberFluxCalc)%dat(1)
        if (size(indx_data%var(iLookINDEX%ixNrgOnly)%dat).ne.0) allocate(indx_data_device%ixNrgOnly(size(indx_data%var(iLookINDEX%ixNrgOnly)%dat),nGRU))
        if (size(indx_data%var(iLookINDEX%ixDomainType)%dat).ne.0) allocate(indx_data_device%ixDomainType(size(indx_data%var(iLookINDEX%ixDomainType)%dat),nGRU))
        do iGRU=1,nGRU
          if (allocated(indx_data_device%ixStateType)) indx_data_device%ixStateType(:,iGRU) = indx_data%var(iLookINDEX%ixStateType)%dat
          if (allocated(indx_data_device%ixHydCanopy)) indx_data_device%ixHydCanopy(iGRU) = indx_data%var(iLookINDEX%ixHydCanopy)%dat(1)
          if (allocated(indx_data_device%ixNrgCanopy)) indx_data_device%ixNrgCanopy(iGRU) = indx_data%var(iLookINDEX%ixNrgCanopy)%dat(1)
          if (allocated(indx_data_device%ixControlVolume)) indx_data_device%ixControlVolume(:,iGRU) = indx_data%var(iLookINDEX%ixControlVolume)%dat
          if (allocated(indx_data_device%layerType)) indx_data_device%layerType(:,iGRU) = indx_data%var(iLookINDEX%layerType)%dat

          if (allocated(indx_data_device%ixHydType)) indx_data_device%ixHydType(:,iGRU) = indx_data%var(iLookINDEX%ixHydType)%dat
          if (allocated(indx_data_device%ixNrgLayer)) indx_data_device%ixNrgLayer(:,iGRU) = indx_data%var(iLookINDEX%ixNrgLayer)%dat
          if (allocated(indx_data_device%ixHydLayer)) indx_data_device%ixHydLayer(:,iGRU) = indx_data%var(iLookINDEX%ixHydLayer)%dat
          if (size(indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat).ne.0) then
            indx_data_device%ixSnowSoilHyd(:,iGRU) = indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat
          else
            indx_data_device%ixSnowSoilHyd(:,iGRU) = integerMissing
          endif
          if (size(indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat).ne.0) then
            indx_data_device%ixSnowSoilNrg(:,iGRU) = indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat
          else
            indx_data_device%ixSnowSoilNrg(:,iGRU) = integerMissing
          endif
          if (size(indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat).ne.0) then
            indx_data_device%ixSoilOnlyHyd(:,iGRU) = indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat
          else
            indx_data_device%ixSoilOnlyHyd(:,iGRU) = integerMissing
          endif
          if (size(indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat).ne.0) then
            indx_data_device%ixSnowOnlyHyd(:,iGRU) = indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat
          else
            indx_data_device%ixSnowOnlyHyd(:,iGRU) = integerMissing
          endif
          if (size(indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat).ne.0) then
            indx_data_device%ixSoilOnlyNrg(:,iGRU) = indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat
          else
            indx_data_device%ixSoilOnlyNrg(:,iGRU) = integerMissing
          endif
          if (size(indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat).ne.0) then
            indx_data_device%ixSnowOnlyNrg(:,iGRU) = indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat
          else
            indx_data_device%ixSnowOnlyNrg(:,iGRU) = integerMissing
          endif
          if (allocated(indx_data_device%ixDomainType)) indx_data_device%ixDomainType(:,iGRU) = indx_data%var(iLookINDEX%ixDomainType)%dat
          if (allocated(indx_data_device%ixNrgOnly)) indx_data_device%ixNrgOnly(:,iGRU) = indx_data%var(iLookINDEX%ixNrgOnly)%dat
          if (size(indx_data%var(iLookINDEX%ixSoilState)%dat).ne.0) indx_data_device%ixSoilState(:,iGRU) = indx_data%var(iLookINDEX%ixSoilState)%dat

        end do
  
      end subroutine allocate_device_indx_data

      subroutine finalize_device_indx_data(indx_data_device, indx_data)
        type(indx_data_device), intent(in) :: indx_data_device
        type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
        indx_data%var(iLookINDEX%numberFluxCalc)%dat(1) = indx_data_device%numberFluxCalc
        if (allocated(indx_data_device%ixNrgLayer)) then
          ! print*, indx_data%var(iLookINDEX%ixNrgLayer)%dat
          indx_data%var(iLookINDEX%ixNrgLayer)%dat = indx_data_device%ixNrgLayer(:,1)
          ! print*, indx_data%var(iLookINDEX%ixNrgLayer)%dat
        end if
        if (allocated(indx_data_device%ixHydLayer)) indx_data%var(iLookINDEX%ixHydLayer)%dat = indx_data_device%ixHydLayer(:,1)
        indx_data%var(iLookINDEX%nSoil)%dat(1) = indx_data_device%nSoil
        indx_data%var(iLookINDEX%nLayers)%dat(1) = indx_data_device%nLayers_d(1)
        if (allocated(indx_data_device%nSnow)) indx_data%var(iLookINDEX%nSnow)%dat(1) = indx_data_device%nSnow(1)
        if (allocated(indx_data_device%ixCasNrg)) indx_data%var(iLookINDEX%ixCasNrg)%dat(1) = indx_data_device%ixCasNrg(1)
        if (allocated(indx_data_device%ixVegNrg)) indx_data%var(iLookINDEX%ixVegNrg)%dat(1) = indx_data_device%ixVegNrg(1)
        if (allocated(indx_data_device%ixVegHyd)) indx_data%var(iLookINDEX%ixVegHyd)%dat(1) = indx_data_device%ixVegHyd(1)
        if (allocated(indx_data_device%ixAqWat)) indx_data%var(iLookINDEX%ixAqWat)%dat(1) = indx_data_device%ixAqWat(1)
        if (allocated(indx_data_device%ixTopHyd)) indx_data%var(iLookINDEX%ixTopHyd)%dat(1) = indx_data_device%ixTopHyd(1)
        if (allocated(indx_data_device%ixTopNrg)) indx_data%var(iLookINDEX%ixTopNrg)%dat(1) = indx_data_device%ixTopNrg(1)
        if (allocated(indx_data_device%ixSnowSoilNrg)) indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat = indx_data_device%ixSnowSoilNrg(:,1)
        if (allocated(indx_data_device%ixSnowSoilHyd )) indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat = indx_data_device%ixSnowSoilHyd(:,1)
        indx_data%var(iLookINDEX%nSnowSoilNrg)%dat(1) = indx_data_device%nLayers_d(1)
        indx_data%var(iLookINDEX%nSnowSoilHyd)%dat(1) = indx_data_device%nLayers_d(1)
        indx_data%var(iLookINDEX%nSnowOnlyNrg)%dat(1) = indx_data_device%nSnow(1)
        indx_data%var(iLookINDEX%nSnowOnlyHyd)%dat(1) = indx_data_device%nSnow(1)
        indx_data%var(iLookINDEX%nSoilOnlyNrg)%dat(1) = indx_data_device%nSoil
        indx_data%var(iLookINDEX%nSoilOnlyHyd)%dat(1) = indx_data_device%nSoil
        if (allocated(indx_data_device%ixHydType)) indx_data%var(iLookINDEX%ixHydType)%dat = indx_data_device%ixHydType(:,1)
        if (allocated(indx_data_device%ixStateType)) indx_data%var(iLookINDEX%ixStateType)%dat = indx_data_device%ixStateType(:,1)
        if (allocated(indx_data_device%ixHydCanopy)) indx_data%var(iLookINDEX%ixHydCanopy)%dat(1) = indx_data_device%ixHydCanopy(1)
        if (allocated(indx_data_device%ixNrgCanopy)) indx_data%var(iLookINDEX%ixNrgCanopy)%dat(1) = indx_data_device%ixNrgCanopy(1)
        if (allocated(indx_data_device%ixControlVolume)) indx_data%var(iLookINDEX%ixControlVolume)%dat = indx_data_device%ixControlVolume(:,1)
        if (allocated(indx_data_device%layerType)) indx_data%var(iLookINDEX%layerType)%dat = indx_data_device%layerType(:,1)
        if (allocated(indx_data_device%ixSoilOnlyHyd)) indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat = indx_data_device%ixSoilOnlyHyd(:,1)
        if (allocated(indx_data_device%ixSnowOnlyHyd)) indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat = indx_data_device%ixSnowOnlyHyd(:,1)
        if (allocated(indx_data_device%ixSoilOnlyNrg)) indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat = indx_data_device%ixSoilOnlyNrg(:,1)
        if (allocated(indx_data_device%ixSnowOnlyNrg)) indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat = indx_data_device%ixSnowOnlyNrg(:,1)
      
        if (allocated(indx_data_device%ixSoilState)) indx_data%var(iLookINDEX%ixSoilState)%dat = indx_data_device%ixSoilState(:,1)
        if (allocated(indx_data_device%ixDomainType)) indx_data%var(iLookINDEX%ixDomainType)%dat = indx_data_device%ixDomainType(:,1)
        if (allocated(indx_data_device%ixNrgOnly)) indx_data%var(iLookINDEX%ixNrgOnly)%dat = indx_data_device%ixNrgOnly(:,1)
      

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
      

  subroutine allocate_device_prog_data(prog_data_device, prog_data, nGRU, nLayers, nSoil)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    type(var_dlength),intent(in) :: prog_data
    integer(i4b),intent(in) :: nGRU,nLayers,nSoil

    integer(i4b) :: iGRU
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

  end subroutine allocate_device_prog_data  

  subroutine allocate_device_prog_temp(prog_data_device, nGRU, nLayers, nSoil)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    integer(i4b),intent(in) :: nGRU,nLayers,nSoil

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

  subroutine finalize_device_prog_data(prog_data_device, prog_data,nLayers,nSoil)
    implicit none
    type(prog_data_device), intent(inout) :: prog_data_device
    type(var_dlength),intent(inout) :: prog_data
    integer(i4b),intent(in) :: nLayers,nSoil
    ! integer(i4b),intent(in) :: nGRU,nLayers,nSoil

    ! integer(i4b) :: iGRU

    prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nSpecBands) = prog_data_device%spectralSnowAlbedoDiffuse(1:nSpecBands,1)
    prog_data%var(iLookPROG%mLayerTemp)%dat(1:nLayers) = prog_data_device%mLayerTemp(1:nLayers,1)
    prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nLayers) = prog_data_device%mLayerVolFracIce(1:nLayers,1)
    prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nLayers) = prog_data_device%mLayerVolFracLiq(1:nLayers,1)
    prog_data%var(iLookPROG%mLayerVolFracWat)%dat(1:nLayers) = prog_data_device%mLayerVolFracWat(1:nLayers,1)
    prog_data%var(iLookPROG%mLayerMatricHead)%dat(1:nSoil) = prog_data_device%mLayerMatricHead(1:nSoil,1)
    prog_data%var(iLookPROG%mLayerEnthalpy)%dat(1:nLayers) = prog_data_device%mLayerEnthalpy(1:nLayers,1)
    prog_data%var(iLookPROG%mLayerDepth)%dat(1:nLayers) = prog_data_device%mLayerDepth(1:nLayers,1)
    prog_data%var(iLookPROG%mLayerHeight)%dat(0:nLayers) = prog_data_device%mLayerHeight(0:nLayers,1)
    prog_data%var(iLookPROG%iLayerHeight)%dat(0:nLayers) = prog_data_device%iLayerHeight(0:nLayers,1)

    
    prog_data%var(iLookPROG%dt_init)%dat(1) = prog_data_device%dt_init
    prog_data%var(iLookPROG%scalarCanopyIce)%dat(1) = prog_data_device%scalarCanopyIce(1)
    prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1) = prog_data_device%scalarCanopyLiq(1)
    prog_data%var(iLookPROG%scalarCanopyWat)%dat(1) = prog_data_device%scalarCanopyWat(1)
    prog_data%var(iLookPROG%scalarCanairTemp)%dat(1) = prog_data_device%scalarCanairTemp(1)
    prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1) = prog_data_device%scalarCanopyTemp(1)
    prog_data%var(iLookPROG%scalarSnowAlbedo)%dat(1) = prog_data_device%scalarSnowAlbedo(1)
    prog_data%var(iLookPROG%scalarSnowDepth)%dat(1) = prog_data_device%scalarSnowDepth(1)
    prog_data%var(iLookPROG%scalarSWE)%dat(1) = prog_data_device%scalarSWE(1)
    prog_data%var(iLookPROG%scalarSfcMeltPond)%dat(1) = prog_data_device%scalarSfcMeltPond(1)
    prog_data%var(iLookPROG%scalarCanairEnthalpy)%dat(1) = prog_data_device%scalarCanairEnthalpy(1)
    prog_data%var(iLookPROG%scalarCanopyEnthalpy)%dat(1) = prog_data_device%scalarCanopyEnthalpy(1)
    prog_data%var(iLookPROG%scalarAquiferStorage)%dat(1) = prog_data_device%scalarAquiferStorage(1)
    prog_data%var(iLookPROG%scalarSurfaceTemp)%dat(1) = prog_data_device%scalarSurfaceTemp(1)

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
    type(var_dlength), intent(in) :: flux_data
    integer(i4b),intent(in) :: nSnow,nSoil,nGRU,nLayers

    integer(i4b) :: iGRU
      
    ! print*, maxSnowLayers
  
    allocate(flux_data_device%scalarCanairNetNrgFlux(nGRU))
    allocate(flux_data_device%scalarCanopyNetNrgFlux(nGRU))
    allocate(flux_data_device%scalarGroundNetNrgFlux(nGRU))
    allocate(flux_data_device%scalarCanopyNetLiqFlux(nGRU))
    allocate(flux_data_device%scalarRainfall(nGRU))
    allocate(flux_data_device%scalarSnowfall(nGRU))
    allocate(flux_data_device%spectralIncomingDirect(nSpecBands,nGRU))
    allocate(flux_data_device%spectralIncomingDiffuse(nSpecBands,nGRU))
    allocate(flux_data_device%scalarCanopySunlitPAR(nGRU))
    allocate(flux_data_device%scalarCanopyShadedPAR(nGRU))
    allocate(flux_data_device%spectralBelowCanopyDirect(nSpecBands,nGRU))
    allocate(flux_data_device%spectralBelowCanopyDiffuse(nSpecBands,nGRU))
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
    allocate(flux_data_device%scalarSurfaceRunoff_IE(nGRU))
    allocate(flux_data_device%scalarSurfaceRunoff_SE(nGRU))


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
    do iGRU=1,nGRU
      flux_data_device%spectralIncomingDirect(:,iGRU) = flux_data%var(iLookFLUX%spectralIncomingDirect)%dat
      flux_data_device%spectralIncomingDiffuse(:,iGRU) = flux_data%var(iLookFLUX%spectralIncomingDiffuse)%dat
      flux_data_device%spectralBelowCanopyDirect(:,iGRU) = flux_data%var(iLookFLUX%spectralBelowCanopyDirect)%dat
      flux_data_device%spectralBelowCanopyDiffuse(:,iGRU) = flux_data%var(iLookFLUX%spectralBelowCanopyDiffuse)%dat
      flux_data_device%iLayerConductiveFlux_m(1:nSoil+maxSnowLayers,iGRU) = flux_data%var(iLookFLUX%iLayerConductiveFlux)%dat
      flux_data_device%iLayerAdvectiveFlux_m(1:nSoil+maxSnowLayers,iGRU) = flux_data%var(iLookFLUX%iLayerAdvectiveFlux)%dat
      flux_data_device%iLayerNrgFlux_m(0:nSoil+maxSnowLayers,iGRU) = flux_data%var(iLookFLUX%iLayerNrgFlux)%dat
      flux_data_device%mLayerNrgFlux_m(1:nSoil+maxSnowLayers,iGRU) = flux_data%var(iLookFLUX%mLayerNrgFlux)%dat
      flux_data_device%iLayerLiqFluxSnow_m(0:nSnow,iGRU) = flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat
      if (nSnow.ne.0) flux_data_device%mLayerLiqFluxSnow_m(1:nSnow,iGRU) = flux_data%var(iLookFLUX%mLayerLiqFluxSnow)%dat
      flux_data_device%mLayerSatHydCondMP_m(:,iGRU) = flux_data%var(iLookFLUX%mLayerSatHydCondMP)%dat
      flux_data_device%mLayerSatHydCond_m(:,iGRU) = flux_data%var(iLookFLUX%mLayerSatHydCond)%dat
      flux_data_device%iLayerSatHydCond_m(:,iGRU) = flux_data%var(iLookFLUX%iLayerSatHydCond)%dat
      flux_data_device%mLayerHydCond_m(:,iGRU) = flux_data%var(iLookFLUX%mLayerHydCond)%dat
      flux_data_device%iLayerLiqFluxSoil_m(:,iGRU) = flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat
      flux_data_device%mLayerLiqFluxSoil_m(:,iGRU) = flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat
      flux_data_device%mLayerBaseflow_m(:,iGRU) = flux_data%var(iLookFLUX%mLayerBaseflow)%dat
      flux_data_device%mLayerColumnInflow(:,iGRU) = flux_data%var(iLookFLUX%mLayerColumnInflow)%dat
      flux_data_device%mLayerColumnOutflow_m(:,iGRU) = flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat
    end do
    flux_data_device%scalarCanairNetNrgFlux = flux_data%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1)
    flux_data_device%scalarCanopyNetNrgFlux = flux_data%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1)
    flux_data_device%scalarGroundNetNrgFlux = flux_data%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1)
    flux_data_device%scalarCanopyNetLiqFlux = flux_data%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1)
    flux_data_device%scalarRainfall = flux_data%var(iLookFLUX%scalarRainfall)%dat(1)
    flux_data_device%scalarSnowfall = flux_data%var(iLookFLUX%scalarSnowfall)%dat(1)
    flux_data_device%scalarCanopySunlitPAR = flux_data%var(iLookFLUX%scalarCanopySunlitPAR)%dat(1)
    flux_data_device%scalarCanopyShadedPAR = flux_data%var(iLookFLUX%scalarCanopyShadedPAR)%dat(1)
    flux_data_device%scalarBelowCanopySolar = flux_data%var(iLookFLUX%scalarBelowCanopySolar)%dat(1)
    flux_data_device%scalarCanopyAbsorbedSolar = flux_data%var(iLookFLUX%scalarCanopyAbsorbedSolar)%dat(1)
    flux_data_device%scalarGroundAbsorbedSolar = flux_data%var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1)
    flux_data_device%scalarLWRadCanopy = flux_data%var(iLookFLUX%scalarLWRadCanopy)%dat(1)
    flux_data_device%scalarLWRadGround = flux_data%var(iLookFLUX%scalarLWRadGround)%dat(1)
    flux_data_device%scalarLWRadUbound2Canopy = flux_data%var(iLookFLUX%scalarLWRadUbound2Canopy)%dat(1)
    flux_data_device%scalarLWRadUbound2Ground = flux_data%var(iLookFLUX%scalarLWRadUbound2Ground)%dat(1)
    flux_data_device%scalarLWRadUbound2Ubound = flux_data%var(iLookFLUX%scalarLWRadUbound2Ubound)%dat(1)
    flux_data_device%scalarLWRadCanopy2Ubound = flux_data%var(iLookFLUX%scalarLWRadCanopy2Ubound)%dat(1)
    flux_data_device%scalarLWRadCanopy2Ground = flux_data%var(iLookFLUX%scalarLWRadCanopy2Ground)%dat(1)
    flux_data_device%scalarLWRadCanopy2Canopy = flux_data%var(iLookFLUX%scalarLWRadCanopy2Canopy)%dat(1)
    flux_data_device%scalarLWRadGround2Ubound = flux_data%var(iLookFLUX%scalarLWRadGround2Ubound)%dat(1)
    flux_data_device%scalarLWRadGround2Canopy = flux_data%var(iLookFLUX%scalarLWRadGround2Canopy)%dat(1)
    flux_data_device%scalarLWNetCanopy = flux_data%var(iLookFLUX%scalarLWNetCanopy)%dat(1)
    flux_data_device%scalarLWNetGround = flux_data%var(iLookFLUX%scalarLWNetGround)%dat(1)
    flux_data_device%scalarLWNetUbound = flux_data%var(iLookFLUX%scalarLWNetUbound)%dat(1)
    flux_data_device%scalarEddyDiffusCanopyTop = flux_data%var(iLookFLUX%scalarEddyDiffusCanopyTop)%dat(1)
    flux_data_device%scalarFrictionVelocity = flux_data%var(iLookFLUX%scalarFrictionVelocity)%dat(1)
    flux_data_device%scalarWindspdCanopyTop = flux_data%var(iLookFLUX%scalarWindspdCanopyTop)%dat(1)
    flux_data_device%scalarWindspdCanopyBottom = flux_data%var(iLookFLUX%scalarWindspdCanopyBottom)%dat(1)
    flux_data_device%scalarGroundResistance = flux_data%var(iLookFLUX%scalarGroundResistance)%dat(1)
    flux_data_device%scalarCanopyResistance = flux_data%var(iLookFLUX%scalarCanopyResistance)%dat(1)
    flux_data_device%scalarLeafResistance = flux_data%var(iLookFLUX%scalarLeafResistance)%dat(1)
    flux_data_device%scalarSoilResistance = flux_data%var(iLookFLUX%scalarSoilResistance)%dat(1)
    flux_data_device%scalarSenHeatTotal = flux_data%var(iLookFLUX%scalarSenHeatTotal)%dat(1)
    flux_data_device%scalarSenHeatCanopy = flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1)
    flux_data_device%scalarSenHeatGround = flux_data%var(iLookFLUX%scalarSenHeatGround)%dat(1)
    flux_data_device%scalarLatHeatTotal = flux_data%var(iLookFLUX%scalarLatHeatTotal)%dat(1)
    flux_data_device%scalarLatHeatCanopyEvap = flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1)
    flux_data_device%scalarLatHeatCanopyTrans = flux_data%var(iLookFLUX%scalarLatHeatCanopyTrans)%dat(1)
    flux_data_device%scalarLatHeatGround = flux_data%var(iLookFLUX%scalarLatHeatGround)%dat(1)
    flux_data_device%scalarCanopyAdvectiveHeatFlux = flux_data%var(iLookFLUX%scalarCanopyAdvectiveHeatFlux)%dat(1)
    flux_data_device%scalarGroundAdvectiveHeatFlux = flux_data%var(iLookFLUX%scalarGroundAdvectiveHeatFlux)%dat(1)
    flux_data_device%scalarCanopySublimation = flux_data%var(iLookFLUX%scalarCanopySublimation)%dat(1)
    flux_data_device%scalarSnowSublimation = flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1)
    flux_data_device%scalarStomResistSunlit = flux_data%var(iLookFLUX%scalarStomResistSunlit)%dat(1)
    flux_data_device%scalarStomResistShaded = flux_data%var(iLookFLUX%scalarStomResistShaded)%dat(1)
    flux_data_device%scalarPhotosynthesisSunlit = flux_data%var(iLookFLUX%scalarPhotosynthesisSunlit)%dat(1)
    flux_data_device%scalarPhotosynthesisShaded = flux_data%var(iLookFLUX%scalarPhotosynthesisShaded)%dat(1)
    flux_data_device%scalarCanopyTranspiration = flux_data%var(iLookFLUX%scalarCanopyTranspiration)%dat(1)
    flux_data_device%scalarCanopyEvaporation = flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1)
    flux_data_device%scalarGroundEvaporation = flux_data%var(iLookFLUX%scalarGroundEvaporation)%dat(1)
    flux_data_device%scalarThroughfallSnow = flux_data%var(iLookFLUX%scalarThroughfallSnow)%dat(1)
    flux_data_device%scalarThroughfallRain = flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1)
    flux_data_device%scalarCanopySnowUnloading = flux_data%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1)
    flux_data_device%scalarCanopyLiqDrainage = flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)
    flux_data_device%scalarCanopyMeltFreeze = flux_data%var(iLookFLUX%scalarCanopyMeltFreeze)%dat(1)
    flux_data_device%scalarSnowDrainage = flux_data%var(iLookFLUX%scalarSnowDrainage)%dat(1)
    flux_data_device%scalarRainPlusMelt = flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1)
    flux_data_device%scalarMaxInfilRate = flux_data%var(iLookFLUX%scalarMaxInfilRate)%dat(1)
    flux_data_device%scalarInfiltration = flux_data%var(iLookFLUX%scalarInfiltration)%dat(1)
    flux_data_device%scalarExfiltration = flux_data%var(iLookFLUX%scalarExfiltration)%dat(1)
    flux_data_device%scalarSurfaceRunoff = flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)
    flux_data_device%scalarSurfaceRunoff_IE = flux_data%var(iLookFLUX%scalarSurfaceRunoff_IE)%dat(1)
    flux_data_device%scalarSurfaceRunoff_SE = flux_data%var(iLookFLUX%scalarSurfaceRunoff_SE)%dat(1)

    flux_data_device%scalarSoilBaseflow = flux_data%var(iLookFLUX%scalarSoilBaseflow)%dat(1)
    flux_data_device%scalarSoilDrainage = flux_data%var(iLookFLUX%scalarSoilDrainage)%dat(1)
    flux_data_device%scalarAquiferRecharge = flux_data%var(iLookFLUX%scalarAquiferRecharge)%dat(1)
    flux_data_device%scalarAquiferTranspire = flux_data%var(iLookFLUX%scalarAquiferTranspire)%dat(1)
    flux_data_device%scalarAquiferBaseflow = flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1)
    flux_data_device%scalarTotalET = flux_data%var(iLookFLUX%scalarTotalET)%dat(1)
    flux_data_device%scalarTotalRunoff = flux_data%var(iLookFLUX%scalarTotalRunoff)%dat(1)
    flux_data_device%scalarNetRadiation = flux_data%var(iLookFLUX%scalarNetRadiation)%dat(1)
      
  end subroutine allocate_device_flux_data

  subroutine allocate_device_flux_prev(flux_data_device,nSnow,nSoil,nGRU,nLayers)
    USE globalData,only:maxSnowLayers                           ! maximum number of snow layers

    type(flux_data_device), intent(inout) :: flux_data_device
    integer(i4b),intent(in) :: nSnow,nSoil,nGRU,nLayers

    integer(i4b) :: iGRU
      
  
    allocate(flux_data_device%scalarCanairNetNrgFlux(nGRU))
    allocate(flux_data_device%scalarCanopyNetNrgFlux(nGRU))
    allocate(flux_data_device%scalarGroundNetNrgFlux(nGRU))
    allocate(flux_data_device%scalarCanopyNetLiqFlux(nGRU))
    allocate(flux_data_device%scalarRainfall(nGRU))
    allocate(flux_data_device%scalarSnowfall(nGRU))
    allocate(flux_data_device%spectralIncomingDirect(nSpecBands,nGRU))
    allocate(flux_data_device%spectralIncomingDiffuse(nSpecBands,nGRU))
    allocate(flux_data_device%scalarCanopySunlitPAR(nGRU))
    allocate(flux_data_device%scalarCanopyShadedPAR(nGRU))
    allocate(flux_data_device%spectralBelowCanopyDirect(nSpecBands,nGRU))
    allocate(flux_data_device%spectralBelowCanopyDiffuse(nSpecBands,nGRU))
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
    allocate(flux_data_device%scalarSurfaceRunoff_IE(nGRU))
    allocate(flux_data_device%scalarSurfaceRunoff_SE(nGRU))
    
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
    flux_data_device%scalarSurfaceRunoff_IE = 0._rkind
    flux_data_device%scalarSurfaceRunoff_SE = 0._rkind

    flux_data_device%scalarSoilBaseflow = 0._rkind
    flux_data_device%scalarSoilDrainage = 0._rkind
    flux_data_device%scalarAquiferRecharge = 0._rkind
    flux_data_device%scalarAquiferTranspire = 0._rkind
    flux_data_device%scalarAquiferBaseflow = 0._rkind
    flux_data_device%scalarTotalET = 0._rkind
    flux_data_device%scalarTotalRunoff = 0._rkind
    flux_data_device%scalarNetRadiation = 0._rkind
      
  end subroutine allocate_device_flux_prev

  subroutine finalize_device_flux_data(flux_data_device, flux_data,nSnow,nSoil)
    type(flux_data_device), intent(inout) :: flux_data_device
    type(var_dlength), intent(inout) :: flux_data
    integer(i4b),intent(in) :: nSnow,nSoil

    flux_data%var(iLookFLUX%spectralIncomingDirect)%dat = flux_data_device%spectralIncomingDirect(:,1)
    flux_data%var(iLookFLUX%spectralIncomingDiffuse)%dat = flux_data_device%spectralIncomingDiffuse(:,1)
    flux_data%var(iLookFLUX%spectralBelowCanopyDirect)%dat = flux_data_device%spectralBelowCanopyDirect(:,1)
    flux_data%var(iLookFLUX%spectralBelowCanopyDiffuse)%dat = flux_data_device%spectralBelowCanopyDiffuse(:,1)
    flux_data%var(iLookFLUX%iLayerConductiveFlux)%dat = flux_data_device%iLayerConductiveFlux_m(:,1)
    flux_data%var(iLookFLUX%iLayerAdvectiveFlux)%dat = flux_data_device%iLayerAdvectiveFlux_m(:,1)
    flux_data%var(iLookFLUX%iLayerNrgFlux)%dat = flux_data_device%iLayerNrgFlux_m(:,1)
    flux_data%var(iLookFLUX%mLayerNrgFlux)%dat = flux_data_device%mLayerNrgFlux_m(:,1)
    flux_data%var(iLookFLUX%iLayerLiqFluxSnow)%dat = flux_data_device%iLayerLiqFluxSnow_m(0:nSnow,1)
    if (nSnow.ne.0) flux_data%var(iLookFLUX%mLayerLiqFluxSnow)%dat = flux_data_device%mLayerLiqFluxSnow_m(1:nSnow,1)
    flux_data%var(iLookFLUX%mLayerSatHydCondMP)%dat = flux_data_device%mLayerSatHydCondMP_m(:,1)
    flux_data%var(iLookFLUX%mLayerSatHydCond)%dat = flux_data_device%mLayerSatHydCond_m(:,1)
    flux_data%var(iLookFLUX%iLayerSatHydCond)%dat = flux_data_device%iLayerSatHydCond_m(:,1)
    flux_data%var(iLookFLUX%mLayerHydCond)%dat = flux_data_device%mLayerHydCond_m(:,1)
    flux_data%var(iLookFLUX%iLayerLiqFluxSoil)%dat = flux_data_device%iLayerLiqFluxSoil_m(:,1)
    flux_data%var(iLookFLUX%mLayerLiqFluxSoil)%dat = flux_data_device%mLayerLiqFluxSoil_m(:,1)
    flux_data%var(iLookFLUX%mLayerBaseflow)%dat = flux_data_device%mLayerBaseflow_m(:,1)
    flux_data%var(iLookFLUX%mLayerColumnInflow)%dat = flux_data_device%mLayerColumnInflow(:,1)
    flux_data%var(iLookFLUX%mLayerColumnOutflow)%dat = flux_data_device%mLayerColumnOutflow_m(:,1)
    flux_data%var(iLookFLUX%scalarCanairNetNrgFlux)%dat(1) = flux_data_device%scalarCanairNetNrgFlux(1)
    flux_data%var(iLookFLUX%scalarCanopyNetNrgFlux)%dat(1) = flux_data_device%scalarCanopyNetNrgFlux(1)
    flux_data%var(iLookFLUX%scalarGroundNetNrgFlux)%dat(1) = flux_data_device%scalarGroundNetNrgFlux(1)
    flux_data%var(iLookFLUX%scalarCanopyNetLiqFlux)%dat(1) = flux_data_device%scalarCanopyNetLiqFlux(1)
    flux_data%var(iLookFLUX%scalarRainfall)%dat(1) = flux_data_device%scalarRainfall(1)
    flux_data%var(iLookFLUX%scalarSnowfall)%dat(1) = flux_data_device%scalarSnowfall(1)
    flux_data%var(iLookFLUX%scalarCanopySunlitPAR)%dat(1) = flux_data_device%scalarCanopySunlitPAR(1)
    flux_data%var(iLookFLUX%scalarCanopyShadedPAR)%dat(1) = flux_data_device%scalarCanopyShadedPAR(1)
    flux_data%var(iLookFLUX%scalarBelowCanopySolar)%dat(1) = flux_data_device%scalarBelowCanopySolar(1)
    flux_data%var(iLookFLUX%scalarCanopyAbsorbedSolar)%dat(1) = flux_data_device%scalarCanopyAbsorbedSolar(1)
    flux_data%var(iLookFLUX%scalarGroundAbsorbedSolar)%dat(1) = flux_data_device%scalarGroundAbsorbedSolar(1)
    flux_data%var(iLookFLUX%scalarLWRadCanopy)%dat(1) = flux_data_device%scalarLWRadCanopy(1)
    flux_data%var(iLookFLUX%scalarLWRadGround)%dat(1) = flux_data_device%scalarLWRadGround(1)
    flux_data%var(iLookFLUX%scalarLWRadUbound2Canopy)%dat(1) = flux_data_device%scalarLWRadUbound2Canopy(1)
    flux_data%var(iLookFLUX%scalarLWRadUbound2Ground)%dat(1) = flux_data_device%scalarLWRadUbound2Ground(1)
    flux_data%var(iLookFLUX%scalarLWRadUbound2Ubound)%dat(1) = flux_data_device%scalarLWRadUbound2Ubound(1)
    flux_data%var(iLookFLUX%scalarLWRadCanopy2Ubound)%dat(1) = flux_data_device%scalarLWRadCanopy2Ubound(1)
    flux_data%var(iLookFLUX%scalarLWRadCanopy2Ground)%dat(1) = flux_data_device%scalarLWRadCanopy2Ground(1)
    flux_data%var(iLookFLUX%scalarLWRadCanopy2Canopy)%dat(1) = flux_data_device%scalarLWRadCanopy2Canopy(1)
    flux_data%var(iLookFLUX%scalarLWRadGround2Ubound)%dat(1) = flux_data_device%scalarLWRadGround2Ubound(1)
    flux_data%var(iLookFLUX%scalarLWRadGround2Canopy)%dat(1) = flux_data_device%scalarLWRadGround2Canopy(1)
    flux_data%var(iLookFLUX%scalarLWNetCanopy)%dat(1) = flux_data_device%scalarLWNetCanopy(1)
    flux_data%var(iLookFLUX%scalarLWNetGround)%dat(1) = flux_data_device%scalarLWNetGround(1)
    flux_data%var(iLookFLUX%scalarLWNetUbound)%dat(1) = flux_data_device%scalarLWNetUbound(1)
    flux_data%var(iLookFLUX%scalarEddyDiffusCanopyTop)%dat(1) = flux_data_device%scalarEddyDiffusCanopyTop(1)
    flux_data%var(iLookFLUX%scalarFrictionVelocity)%dat(1) = flux_data_device%scalarFrictionVelocity(1)
    flux_data%var(iLookFLUX%scalarWindspdCanopyTop)%dat(1) = flux_data_device%scalarWindspdCanopyTop(1)
    flux_data%var(iLookFLUX%scalarWindspdCanopyBottom)%dat(1) = flux_data_device%scalarWindspdCanopyBottom(1)
    flux_data%var(iLookFLUX%scalarGroundResistance)%dat(1) = flux_data_device%scalarGroundResistance(1)
    flux_data%var(iLookFLUX%scalarCanopyResistance)%dat(1) = flux_data_device%scalarCanopyResistance(1)
    flux_data%var(iLookFLUX%scalarLeafResistance)%dat(1) = flux_data_device%scalarLeafResistance(1)
    flux_data%var(iLookFLUX%scalarSoilResistance)%dat(1) = flux_data_device%scalarSoilResistance(1)
    flux_data%var(iLookFLUX%scalarSenHeatTotal)%dat(1) = flux_data_device%scalarSenHeatTotal(1)
    flux_data%var(iLookFLUX%scalarSenHeatCanopy)%dat(1) = flux_data_device%scalarSenHeatCanopy(1)
    flux_data%var(iLookFLUX%scalarSenHeatGround)%dat(1) = flux_data_device%scalarSenHeatGround(1)
    flux_data%var(iLookFLUX%scalarLatHeatTotal)%dat(1) = flux_data_device%scalarLatHeatTotal(1)
    flux_data%var(iLookFLUX%scalarLatHeatCanopyEvap)%dat(1) = flux_data_device%scalarLatHeatCanopyEvap(1)
    flux_data%var(iLookFLUX%scalarLatHeatCanopyTrans)%dat(1) = flux_data_device%scalarLatHeatCanopyTrans(1)
    flux_data%var(iLookFLUX%scalarLatHeatGround)%dat(1) = flux_data_device%scalarLatHeatGround(1)
    flux_data%var(iLookFLUX%scalarCanopyAdvectiveHeatFlux)%dat(1) = flux_data_device%scalarCanopyAdvectiveHeatFlux(1)
    flux_data%var(iLookFLUX%scalarGroundAdvectiveHeatFlux)%dat(1) = flux_data_device%scalarGroundAdvectiveHeatFlux(1)
    flux_data%var(iLookFLUX%scalarCanopySublimation)%dat(1) = flux_data_device%scalarCanopySublimation(1)
    flux_data%var(iLookFLUX%scalarSnowSublimation)%dat(1) = flux_data_device%scalarSnowSublimation(1)
    flux_data%var(iLookFLUX%scalarStomResistSunlit)%dat(1) = flux_data_device%scalarStomResistSunlit(1)
    flux_data%var(iLookFLUX%scalarStomResistShaded)%dat(1) = flux_data_device%scalarStomResistShaded(1)
    flux_data%var(iLookFLUX%scalarPhotosynthesisSunlit)%dat(1) = flux_data_device%scalarPhotosynthesisSunlit(1)
    flux_data%var(iLookFLUX%scalarPhotosynthesisShaded)%dat(1) = flux_data_device%scalarPhotosynthesisShaded(1)
    flux_data%var(iLookFLUX%scalarCanopyTranspiration)%dat(1) = flux_data_device%scalarCanopyTranspiration(1)
    flux_data%var(iLookFLUX%scalarCanopyEvaporation)%dat(1) = flux_data_device%scalarCanopyEvaporation(1)
    flux_data%var(iLookFLUX%scalarGroundEvaporation)%dat(1) = flux_data_device%scalarGroundEvaporation(1)
    flux_data%var(iLookFLUX%scalarThroughfallSnow)%dat(1) = flux_data_device%scalarThroughfallSnow(1)
    flux_data%var(iLookFLUX%scalarThroughfallRain)%dat(1) = flux_data_device%scalarThroughfallRain(1)
    flux_data%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1) = flux_data_device%scalarCanopySnowUnloading(1)
    flux_data%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1) = flux_data_device%scalarCanopyLiqDrainage(1)
    flux_data%var(iLookFLUX%scalarCanopyMeltFreeze)%dat(1) = flux_data_device%scalarCanopyMeltFreeze(1)
    flux_data%var(iLookFLUX%scalarSnowDrainage)%dat(1) = flux_data_device%scalarSnowDrainage(1)
    flux_data%var(iLookFLUX%scalarRainPlusMelt)%dat(1) = flux_data_device%scalarRainPlusMelt(1)
    flux_data%var(iLookFLUX%scalarMaxInfilRate)%dat(1) = flux_data_device%scalarMaxInfilRate(1)
    flux_data%var(iLookFLUX%scalarInfiltration)%dat(1) = flux_data_device%scalarInfiltration(1)
    flux_data%var(iLookFLUX%scalarExfiltration)%dat(1) = flux_data_device%scalarExfiltration(1)
    flux_data%var(iLookFLUX%scalarSurfaceRunoff)%dat(1) = flux_data_device%scalarSurfaceRunoff(1)
    flux_data%var(iLookFLUX%scalarSurfaceRunoff_IE)%dat(1) = flux_data_device%scalarSurfaceRunoff_IE(1)
    flux_data%var(iLookFLUX%scalarSurfaceRunoff_SE)%dat(1) = flux_data_device%scalarSurfaceRunoff_SE(1)

    flux_data%var(iLookFLUX%scalarSoilBaseflow)%dat(1) = flux_data_device%scalarSoilBaseflow(1)
    flux_data%var(iLookFLUX%scalarSoilDrainage)%dat(1) = flux_data_device%scalarSoilDrainage(1)
    flux_data%var(iLookFLUX%scalarAquiferRecharge)%dat(1) = flux_data_device%scalarAquiferRecharge(1)
    flux_data%var(iLookFLUX%scalarAquiferTranspire)%dat(1) = flux_data_device%scalarAquiferTranspire(1)
    flux_data%var(iLookFLUX%scalarAquiferBaseflow)%dat(1) = flux_data_device%scalarAquiferBaseflow(1)
    flux_data%var(iLookFLUX%scalarTotalET)%dat(1) = flux_data_device%scalarTotalET(1)
    flux_data%var(iLookFLUX%scalarTotalRunoff)%dat(1) = flux_data_device%scalarTotalRunoff(1)
    flux_data%var(iLookFLUX%scalarNetRadiation)%dat(1) = flux_data_device%scalarNetRadiation(1)
      
  end subroutine finalize_device_flux_data

  subroutine deallocate_device_flux_data(flux_data_device,nSnow,nSoil)
    type(flux_data_device), intent(inout) :: flux_data_device
    integer(i4b),intent(in) :: nSnow,nSoil

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
    deallocate(flux_data_device%scalarSurfaceRunoff_IE)
    deallocate(flux_data_device%scalarSurfaceRunoff_SE)

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
    
      subroutine allocate_device_deriv_data(deriv_data_device, deriv_data,nGRU,nSoil,nLayers,nSnow)
        type(deriv_data_device), intent(inout) :: deriv_data_device
        type(var_dlength), intent(inout) :: deriv_data
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
        deriv_data_device%d2Theta_dTkCanopy2 = deriv_data%var(iLookDERIV%d2Theta_dTkCanopy2)%dat(1)
        do iGRU=1,nGRU
          deriv_data_device%d2VolTot_dPsi02_m(:,iGRU) = deriv_data%var(iLookDERIV%d2VolTot_dPsi02)%dat
          deriv_data_device%dCm_dTk_m(:,iGRU) = deriv_data%var(iLookDERIV%dCm_dTk)%dat
          deriv_data_device%dCm_dPsi0_m(:,iGRU) = deriv_data%var(iLookDERIV%dCm_dPsi0)%dat
          deriv_data_device%dCompress_dPsi_m(:,iGRU) = deriv_data%var(iLookDERIV%dCompress_dPsi)%dat
          deriv_data_device%dFracLiqWat_dTk_m(:,iGRU) = deriv_data%var(iLookDERIV%dFracLiqWat_dTk)%dat
          deriv_data_device%dNrgFlux_dTempAbove_m(:,iGRU) = deriv_data%var(iLookDERIV%dNrgFlux_dTempAbove)%dat
          deriv_data_device%dNrgFlux_dTempBelow_m(:,iGRU) = deriv_data%var(iLookDERIV%dNrgFlux_dTempBelow)%dat
          deriv_data_device%dNrgFlux_dWatAbove_m(:,iGRU) = deriv_data%var(iLookDERIV%dNrgFlux_dWatAbove)%dat
          deriv_data_device%dNrgFlux_dWatBelow_m(:,iGRU) = deriv_data%var(iLookDERIV%dNrgFlux_dWatBelow)%dat
          deriv_data_device%dPsiLiq_dPsi0_m(:,iGRU) = deriv_data%var(iLookDERIV%dPsiLiq_dPsi0)%dat
          deriv_data_device%dPsiLiq_dTemp_m(:,iGRU) = deriv_data%var(iLookDERIV%dPsiLiq_dTemp)%dat
          deriv_data_device%dTemp_dEnthalpy_m(:,iGRU) = deriv_data%var(iLookDERIV%dTemp_dEnthalpy)%dat
          deriv_data_device%dTemp_dPsi0_m(:,iGRU) = deriv_data%var(iLookDERIV%dTemp_dPsi0)%dat
          deriv_data_device%dTemp_dTheta_m(:,iGRU) = deriv_data%var(iLookDERIV%dTemp_dTheta)%dat
          deriv_data_device%dThermalC_dTempAbove_m(:,iGRU) = deriv_data%var(iLookDERIV%dThermalC_dTempAbove)%dat
          deriv_data_device%dThermalC_dTempBelow_m(:,iGRU) = deriv_data%var(iLookDERIV%dThermalC_dTempBelow)%dat
          deriv_data_device%dThermalC_dWatAbove_m(:,iGRU) = deriv_data%var(iLookDERIV%dThermalC_dWatAbove)%dat
          deriv_data_device%dThermalC_dWatBelow_m(:,iGRU) = deriv_data%var(iLookDERIV%dThermalC_dWatBelow)%dat
          deriv_data_device%dVolHtCapBulk_dPsi0_m(:,iGRU) = deriv_data%var(iLookDERIV%dVolHtCapBulk_dPsi0)%dat
          deriv_data_device%dVolHtCapBulk_dTheta_m(:,iGRU) = deriv_data%var(iLookDERIV%dVolHtCapBulk_dTheta)%dat
          deriv_data_device%dVolHtCapBulk_dTk_m(:,iGRU) = deriv_data%var(iLookDERIV%dVolHtCapBulk_dTk)%dat
          deriv_data_device%dVolTot_dPsi0_m(:,iGRU) = deriv_data%var(iLookDERIV%dVolTot_dPsi0)%dat
          deriv_data_device%dq_dHydStateAbove_m(:,iGRU) = deriv_data%var(iLookDERIV%dq_dHydStateAbove)%dat
          deriv_data_device%dq_dHydStateBelow_m(:,iGRU) = deriv_data%var(iLookDERIV%dq_dHydStateBelow)%dat
          deriv_data_device%dq_dHydStateLayerSurfVec_m(:,iGRU) = deriv_data%var(iLookDERIV%dq_dHydStateLayerSurfVec)%dat
          deriv_data_device%dq_dNrgStateAbove_m(:,iGRU) = deriv_data%var(iLookDERIV%dq_dNrgStateAbove)%dat
          deriv_data_device%dq_dNrgStateBelow_m(:,iGRU) = deriv_data%var(iLookDERIV%dq_dNrgStateBelow)%dat
          deriv_data_device%dq_dNrgStateLayerSurfVec_m(:,iGRU) = deriv_data%var(iLookDERIV%dq_dNrgStateLayerSurfVec)%dat
          deriv_data_device%iLayerLiqFluxSnowDeriv_m(:,iGRU) = deriv_data%var(iLookDERIV%iLayerLiqFluxSnowDeriv)%dat
          deriv_data_device%mLayerd2Theta_dTk2_m(:,iGRU) = deriv_data%var(iLookDERIV%mLayerd2Theta_dTk2)%dat
          deriv_data_device%mLayerdPsi_dTheta_m(:,iGRU) = deriv_data%var(iLookDERIV%mLayerdPsi_dTheta)%dat
          deriv_data_device%mLayerdTheta_dPsi_m(:,iGRU) = deriv_data%var(iLookDERIV%mLayerdTheta_dPsi)%dat
          deriv_data_device%mLayerdTheta_dTk_m(:,iGRU) = deriv_data%var(iLookDERIV%mLayerdTheta_dTk)%dat
          deriv_data_device%mLayerdTrans_dCanWat_m(:,iGRU) = deriv_data%var(iLookDERIV%mLayerdTrans_dCanWat)%dat
          deriv_data_device%mLayerdTrans_dTCanair_m(:,iGRU) = deriv_data%var(iLookDERIV%mLayerdTrans_dTCanair)%dat
          deriv_data_device%mLayerdTrans_dTCanopy_m(:,iGRU) = deriv_data%var(iLookDERIV%mLayerdTrans_dTCanopy)%dat
          deriv_data_device%mLayerdTrans_dTGround_m(:,iGRU) = deriv_data%var(iLookDERIV%mLayerdTrans_dTGround)%dat

        enddo
        deriv_data_device%dAquiferTrans_dCanWat = deriv_data%var(iLookDERIV%dAquiferTrans_dCanWat)%dat(1)
        deriv_data_device%dAquiferTrans_dTCanair = deriv_data%var(iLookDERIV%dAquiferTrans_dTCanair)%dat(1)
        deriv_data_device%dAquiferTrans_dTCanopy = deriv_data%var(iLookDERIV%dAquiferTrans_dTCanopy)%dat(1)
        deriv_data_device%dAquiferTrans_dTGround = deriv_data%var(iLookDERIV%dAquiferTrans_dTGround)%dat(1)
        deriv_data_device%dBaseflow_dAquifer = deriv_data%var(iLookDERIV%dBaseflow_dAquifer)%dat(1)
        deriv_data_device%dCanLiq_dTcanopy = deriv_data%var(iLookDERIV%dCanLiq_dTcanopy)%dat(1)
        deriv_data_device%dCanairNetFlux_dCanairTemp = deriv_data%var(iLookDERIV%dCanairNetFlux_dCanairTemp)%dat(1)
        deriv_data_device%dCanairNetFlux_dCanopyTemp = deriv_data%var(iLookDERIV%dCanairNetFlux_dCanopyTemp)%dat(1)
        deriv_data_device%dCanairNetFlux_dGroundTemp = deriv_data%var(iLookDERIV%dCanairNetFlux_dGroundTemp)%dat(1)
        deriv_data_device%dCanairTemp_dEnthalpy = deriv_data%var(iLookDERIV%dCanairTemp_dEnthalpy)%dat(1)
        deriv_data_device%dCanopyEvaporation_dCanWat = deriv_data%var(iLookDERIV%dCanopyEvaporation_dCanWat)%dat(1)
        deriv_data_device%dCanopyEvaporation_dTCanair = deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanair)%dat(1)
        deriv_data_device%dCanopyEvaporation_dTCanopy = deriv_data%var(iLookDERIV%dCanopyEvaporation_dTCanopy)%dat(1)
        deriv_data_device%dCanopyEvaporation_dTGround = deriv_data%var(iLookDERIV%dCanopyEvaporation_dTGround)%dat(1)
        deriv_data_device%dCanopyNetFlux_dCanWat = deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanWat)%dat(1)
        deriv_data_device%dCanopyNetFlux_dCanairTemp = deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanairTemp)%dat(1)
        deriv_data_device%dCanopyNetFlux_dCanopyTemp = deriv_data%var(iLookDERIV%dCanopyNetFlux_dCanopyTemp)%dat(1)
        deriv_data_device%dCanopyNetFlux_dGroundTemp = deriv_data%var(iLookDERIV%dCanopyNetFlux_dGroundTemp)%dat(1)
        deriv_data_device%dCanopyTemp_dCanWat = deriv_data%var(iLookDERIV%dCanopyTemp_dCanWat)%dat(1)
        deriv_data_device%dCanopyTemp_dEnthalpy = deriv_data%var(iLookDERIV%dCanopyTemp_dEnthalpy)%dat(1)
        deriv_data_device%dCanopyTrans_dCanWat = deriv_data%var(iLookDERIV%dCanopyTrans_dCanWat)%dat(1)
        deriv_data_device%dCanopyTrans_dTCanair = deriv_data%var(iLookDERIV%dCanopyTrans_dTCanair)%dat(1)
        deriv_data_device%dCanopyTrans_dTCanopy = deriv_data%var(iLookDERIV%dCanopyTrans_dTCanopy)%dat(1)
        deriv_data_device%dCanopyTrans_dTGround = deriv_data%var(iLookDERIV%dCanopyTrans_dTGround)%dat(1)
        deriv_data_device%dCm_dTkCanopy = deriv_data%var(iLookDERIV%dCm_dTkCanopy)%dat(1)
        deriv_data_device%dFracLiqVeg_dTkCanopy = deriv_data%var(iLookDERIV%dFracLiqVeg_dTkCanopy)%dat(1)
        deriv_data_device%dGroundEvaporation_dCanWat = deriv_data%var(iLookDERIV%dGroundEvaporation_dCanWat)%dat(1)
        deriv_data_device%dGroundEvaporation_dTCanair = deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanair)%dat(1)
        deriv_data_device%dGroundEvaporation_dTCanopy = deriv_data%var(iLookDERIV%dGroundEvaporation_dTCanopy)%dat(1)
        deriv_data_device%dGroundEvaporation_dTGround = deriv_data%var(iLookDERIV%dGroundEvaporation_dTGround)%dat(1)
        deriv_data_device%dGroundNetFlux_dCanWat = deriv_data%var(iLookDERIV%dGroundNetFlux_dCanWat)%dat(1)
        deriv_data_device%dGroundNetFlux_dCanairTemp = deriv_data%var(iLookDERIV%dGroundNetFlux_dCanairTemp)%dat(1)
        deriv_data_device%dGroundNetFlux_dCanopyTemp = deriv_data%var(iLookDERIV%dGroundNetFlux_dCanopyTemp)%dat(1)
        deriv_data_device%dGroundNetFlux_dGroundTemp = deriv_data%var(iLookDERIV%dGroundNetFlux_dGroundTemp)%dat(1)
        deriv_data_device%dTheta_dTkCanopy = deriv_data%var(iLookDERIV%dTheta_dTkCanopy)%dat(1)
        deriv_data_device%dVolHtCapBulk_dCanWat = deriv_data%var(iLookDERIV%dVolHtCapBulk_dCanWat)%dat(1)
        deriv_data_device%dVolHtCapBulk_dTkCanopy = deriv_data%var(iLookDERIV%dVolHtCapBulk_dTkCanopy)%dat(1)
        deriv_data_device%scalarCanopyLiqDeriv = deriv_data%var(iLookDERIV%scalarCanopyLiqDeriv)%dat(1)
        deriv_data_device%scalarCanopyLiqDrainageDeriv = deriv_data%var(iLookDERIV%scalarCanopyLiqDrainageDeriv)%dat(1)
        deriv_data_device%scalarThroughfallRainDeriv = deriv_data%var(iLookDERIV%scalarThroughfallRainDeriv)%dat(1)
      end subroutine allocate_device_deriv_data
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
  subroutine allocate_device_diag_data(diag_data_device,diag_data,nSnow, nGRU, nLayers, nSoil)
    use globalData,only:maxSnowLayers
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_device
    type(var_dlength),intent(in) :: diag_data
    integer(i4b),intent(in) :: nSnow
    integer(i4b),intent(in) :: nGRU,nLayers,nSoil

    integer(i4b) :: iGRU

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
    allocate(diag_data_device%mLayerTranspireLim_m(nSoil,nGRU))
    allocate(diag_data_device%mLayerRootDensity_m(nSoil,nGRU))
    allocate(diag_data_device%scalarAquiferRootFrac(nGRU))
    allocate(diag_data_device%scalarFracLiqVeg(nGRU))
    allocate(diag_data_device%scalarCanopyWetFraction(nGRU))
    allocate(diag_data_device%scalarSnowAge(nGRU))
    allocate(diag_data_device%scalarGroundSnowFraction(nGRU))
    allocate(diag_data_device%spectralSnowAlbedoDirect(nSpecBands,nGRU))
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
  
    diag_data_device%numFluxCalls = diag_data%var(iLookDIAG%numFluxCalls)%dat(1)
    diag_data_device%wallClockTime = diag_data%var(iLookDIAG%wallClockTime)%dat(1)
    diag_data_device%meanStepSize = diag_data%var(iLookDIAG%meanStepSize)%dat(1)

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
    diag_data_device%scalarLambda_drysoil = diag_data%var(iLookDIAG%scalarLambda_drysoil)%dat(1)
    diag_data_device%scalarLambda_wetsoil = diag_data%var(iLookDIAG%scalarLambda_wetsoil)%dat(1)
    diag_data_device%scalarCanopyEnthTemp = diag_data%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1)
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
    diag_data_device%scalarAquiferRootFrac = diag_data%var(iLookDIAG%scalarAquiferRootFrac)%dat(1)
    diag_data_device%scalarFracLiqVeg = diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1)
    diag_data_device%scalarCanopyWetFraction = diag_data%var(iLookDIAG%scalarCanopyWetFraction)%dat(1)
    diag_data_device%scalarSnowAge = diag_data%var(iLookDIAG%scalarSnowAge)%dat(1)
    diag_data_device%scalarGroundSnowFraction = diag_data%var(iLookDIAG%scalarGroundSnowFraction)%dat(1)
    diag_data_device%scalarInfilArea = diag_data%var(iLookDIAG%scalarInfilArea)%dat(1)
    diag_data_device%scalarFrozenArea = diag_data%var(iLookDIAG%scalarFrozenArea)%dat(1)
    diag_data_device%scalarSoilControl = diag_data%var(iLookDIAG%scalarSoilControl)%dat(1)
    diag_data_device%scalarSoilCompress = diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1)
    diag_data_device%scalarTotalSoilLiq = diag_data%var(iLookDIAG%scalarTotalSoilLiq)%dat(1)
    diag_data_device%scalarTotalSoilIce = diag_data%var(iLookDIAG%scalarTotalSoilIce)%dat(1)
    diag_data_device%scalarTotalSoilWat = diag_data%var(iLookDIAG%scalarTotalSoilWat)%dat(1)
    diag_data_device%scalarKappa = diag_data%var(iLookDIAG%scalarKappa)%dat(1)
    diag_data_device%scalarVolLatHt_fus = diag_data%var(iLookDIAG%scalarVolLatHt_fus)%dat(1)
    diag_data_device%balanceCasNrg = diag_data%var(iLookDIAG%balanceCasNrg)%dat(1)
    diag_data_device%balanceVegNrg = diag_data%var(iLookDIAG%balanceVegNrg)%dat(1)
    diag_data_device%balanceSnowNrg = diag_data%var(iLookDIAG%balanceSnowNrg)%dat(1)
    diag_data_device%balanceSoilNrg = diag_data%var(iLookDIAG%balanceSoilNrg)%dat(1)
    diag_data_device%balanceVegMass = diag_data%var(iLookDIAG%balanceVegMass)%dat(1)
    diag_data_device%balanceSnowMass = diag_data%var(iLookDIAG%balanceSnowMass)%dat(1)
    diag_data_device%balanceSoilMass = diag_data%var(iLookDIAG%balanceSoilMass)%dat(1)
    diag_data_device%balanceAqMass = diag_data%var(iLookDIAG%balanceAqMass)%dat(1)
    do iGRU=1,nGRU
      diag_data_device%mLayerEnthTemp(1:nLayers,iGRU) = diag_data%var(iLookDIAG%mLayerEnthTemp)%dat(1:nLayers)
      diag_data_device%mLayerVolHtCapBulk_m(1:nLayers,iGRU) = diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat(1:nLayers)
      diag_data_device%mLayerCm_m(1:nLayers,iGRU) = diag_data%var(iLookDIAG%mLayerCm)%dat(1:nLayers)
      diag_data_device%mLayerThermalC_m(1:nLayers,iGRU) = diag_data%var(iLookDIAG%mLayerThermalC)%dat(1:nLayers)
      diag_data_device%iLayerThermalC_m(0:nLayers,iGRU) = diag_data%var(iLookDIAG%iLayerThermalC)%dat(0:nLayers)
      diag_data_device%spectralAlbGndDirect(1:nSpecBands,iGRU) = diag_data%var(iLookDIAG%spectralAlbGndDirect)%dat(1:nSpecBands)
      diag_data_device%spectralAlbGndDiffuse(1:nSpecBands,iGRU) = diag_data%var(iLookDIAG%spectralAlbGndDiffuse)%dat(1:nSpecBands)
      diag_data_device%mLayerTranspireLim_m(1:nSoil,iGRU) = diag_data%var(iLookDIAG%mLayerTranspireLim)%dat(1:nSoil)
      diag_data_device%mLayerRootDensity_m(1:nSoil,iGRU) = diag_data%var(iLookDIAG%mLayerRootDensity)%dat(1:nSoil)
      diag_data_device%spectralSnowAlbedoDirect(:,iGRU) = diag_data%var(iLookDIAG%spectralSnowAlbedoDirect)%dat
      if (nSnow.ne.0) diag_data_device%mLayerFracLiqSnow_m(1:nSnow,iGRU) = diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat(1:nSnow)
      if (nSnow.ne.0) diag_data_device%mLayerThetaResid_m(1:nSnow,iGRU) = diag_data%var(iLookDIAG%mLayerThetaResid)%dat(1:nSnow)
      if (nSnow.ne.0) diag_data_device%mLayerPoreSpace_m(1:nSnow,iGRU) = diag_data%var(iLookDIAG%mLayerPoreSpace)%dat(1:nSnow)
      diag_data_device%mLayerMeltFreeze(1:nLayers,iGRU) = diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat(1:nLayers)
      diag_data_device%mLayerVolFracAir_m(1:nLayers,iGRU) = diag_data%var(iLookDIAG%mLayerVolFracAir)%dat(1:nLayers)
      diag_data_device%mLayerTcrit(1:nSoil,iGRU) = diag_data%var(iLookDIAG%mLayerTcrit)%dat(1:nSoil)
      diag_data_device%mLayerCompress_m(1:nSoil,iGRU) = diag_data%var(iLookDIAG%mLayerCompress)%dat(1:nSoil)
      diag_data_device%mLayerMatricHeadLiq(1:nSoil,iGRU) = diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat(1:nSoil)
      diag_data_device%scalarVGn_m_m(1:nSoil,iGRU) = diag_data%var(iLookDIAG%scalarVGn_m)%dat(1:nSoil)
      diag_data_device%balanceLayerNrg(1:nLayers,iGRU) = diag_data%var(iLookDIAG%balanceLayerNrg)%dat(1:nLayers)
      diag_data_device%balanceLayerMass(1:nLayers,iGRU) = diag_data%var(iLookDIAG%balanceLayerMass)%dat(1:nLayers)
  
    end do
  end subroutine allocate_device_diag_data
      
  subroutine allocate_device_diag_temp(diag_data_device,nSnow, nGRU, nLayers, nSoil)
    use globaldata,only:maxSnowLayers
    implicit none
    type(diag_data_device), intent(inout) :: diag_data_device
    integer(i4b),intent(in) :: nSnow
    integer(i4b),intent(in) :: nGRU,nLayers,nSoil

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
    allocate(diag_data_device%mLayerTranspireLim_m(nSoil,nGRU))
    allocate(diag_data_device%mLayerRootDensity_m(nSoil,nGRU))
    allocate(diag_data_device%scalarAquiferRootFrac(nGRU))
    allocate(diag_data_device%scalarFracLiqVeg(nGRU))
    allocate(diag_data_device%scalarCanopyWetFraction(nGRU))
    allocate(diag_data_device%scalarSnowAge(nGRU))
    allocate(diag_data_device%scalarGroundSnowFraction(nGRU))
    allocate(diag_data_device%spectralSnowAlbedoDirect(nSpecBands,nGRU))
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


  subroutine finalize_device_diag_data(diag_data_device,diag_data,nSnow)
    type(diag_data_device), intent(inout) :: diag_data_device
    type(var_dlength),intent(inout) :: diag_data
    integer(i4b),intent(in) :: nSnow

    diag_data%var(iLookDIAG%numFluxCalls)%dat(1) = diag_data_device%numFluxCalls
    diag_data%var(iLookDIAG%wallClockTime)%dat(1) = diag_data_device%wallClockTime
    diag_data%var(iLookDIAG%meanStepSize)%dat(1) = diag_data_device%meanStepSize

    diag_data%var(iLookDIAG%numSteps)%dat(1) = diag_data_device%numSteps
    diag_data%var(iLookDIAG%numResEvals)%dat(1) = diag_data_device%numResEvals
    diag_data%var(iLookDIAG%numLinSolvSetups)%dat(1) = diag_data_device%numLinSolvSetups
    diag_data%var(iLookDIAG%numErrTestFails)%dat(1) = diag_data_device%numErrTestFails
    diag_data%var(iLookDIAG%kLast)%dat(1) = diag_data_device%kLast
    diag_data%var(iLookDIAG%kCur)%dat(1) = diag_data_device%kCur
    diag_data%var(iLookDIAG%hInitUsed)%dat(1) = diag_data_device%hInitUsed
    diag_data%var(iLookDIAG%hLast)%dat(1) = diag_data_device%hLast
    diag_data%var(iLookDIAG%hCur)%dat(1) = diag_data_device%hCur
    diag_data%var(iLookDIAG%tCur)%dat(1) = diag_data_device%tCur

    diag_data%var(iLookDIAG%scalarCanopyDepth)%dat(1) = diag_data_device%scalarCanopyDepth(1)
    diag_data%var(iLookDIAG%scalarBulkVolHeatCapVeg)%dat(1) = diag_data_device%scalarBulkVolHeatCapVeg(1)
    diag_data%var(iLookDIAG%scalarCanopyCm)%dat(1) = diag_data_device%scalarCanopyCm(1)
    diag_data%var(iLookDIAG%scalarCanopyEmissivity)%dat(1) = diag_data_device%scalarCanopyEmissivity(1)
    diag_data%var(iLookDIAG%scalarRootZoneTemp)%dat(1) = diag_data_device%scalarRootZoneTemp(1)
    diag_data%var(iLookDIAG%scalarLAI)%dat(1) = diag_data_device%scalarLAI(1)
    diag_data%var(iLookDIAG%scalarSAI)%dat(1) = diag_data_device%scalarSAI(1)
    diag_data%var(iLookDIAG%scalarExposedLAI)%dat(1) = diag_data_device%scalarExposedLAI(1)
    diag_data%var(iLookDIAG%scalarExposedSAI)%dat(1) = diag_data_device%scalarExposedSAI(1)
    diag_data%var(iLookDIAG%scalarAdjMeasHeight)%dat(1) = diag_data_device%scalarAdjMeasHeight(1)
    diag_data%var(iLookDIAG%scalarCanopyIceMax)%dat(1) = diag_data_device%scalarCanopyIceMax(1)
    diag_data%var(iLookDIAG%scalarCanopyLiqMax)%dat(1) = diag_data_device%scalarCanopyLiqMax(1)
    diag_data%var(iLookDIAG%scalarGrowingSeasonIndex)%dat(1) = diag_data_device%scalarGrowingSeasonIndex(1)
    diag_data%var(iLookDIAG%scalarVolHtCap_air)%dat(1) = diag_data_device%scalarVolHtCap_air(1)
    diag_data%var(iLookDIAG%scalarVolHtCap_ice)%dat(1) = diag_data_device%scalarVolHtCap_ice(1)
    diag_data%var(iLookDIAG%scalarVolHtCap_soil)%dat(1) = diag_data_device%scalarVolHtCap_soil(1)
    diag_data%var(iLookDIAG%scalarVolHtCap_water)%dat(1) = diag_data_device%scalarVolHtCap_water(1)
    diag_data%var(iLookDIAG%scalarLambda_drysoil)%dat(1) = diag_data_device%scalarLambda_drysoil(1)
    diag_data%var(iLookDIAG%scalarLambda_wetsoil)%dat(1) = diag_data_device%scalarLambda_wetsoil(1)
    diag_data%var(iLookDIAG%scalarCanopyEnthTemp)%dat(1) = diag_data_device%scalarCanopyEnthTemp(1)
    diag_data%var(iLookDIAG%scalarTotalSoilEnthalpy)%dat(1) = diag_data_device%scalarTotalSoilEnthalpy(1)
    diag_data%var(iLookDIAG%scalarTotalSnowEnthalpy)%dat(1) = diag_data_device%scalarTotalSnowEnthalpy(1)
    diag_data%var(iLookDIAG%scalarVPair)%dat(1) = diag_data_device%scalarVPair(1)
    diag_data%var(iLookDIAG%scalarVP_CanopyAir)%dat(1) = diag_data_device%scalarVP_CanopyAir(1)
    diag_data%var(iLookDIAG%scalarTwetbulb)%dat(1) = diag_data_device%scalarTwetbulb(1)
    diag_data%var(iLookDIAG%scalarSnowfallTemp)%dat(1) = diag_data_device%scalarSnowfallTemp(1)
    diag_data%var(iLookDIAG%scalarNewSnowDensity)%dat(1) = diag_data_device%scalarNewSnowDensity(1)
    diag_data%var(iLookDIAG%scalarO2air)%dat(1) = diag_data_device%scalarO2air(1)
    diag_data%var(iLookDIAG%scalarCO2air)%dat(1) = diag_data_device%scalarCO2air(1)
    diag_data%var(iLookDIAG%windspd_x)%dat(1) = diag_data_device%windspd_x(1)
    diag_data%var(iLookDIAG%windspd_y)%dat(1) = diag_data_device%windspd_y(1)
    diag_data%var(iLookDIAG%scalarCosZenith)%dat(1) = diag_data_device%scalarCosZenith(1)
    diag_data%var(iLookDIAG%scalarFractionDirect)%dat(1) = diag_data_device%scalarFractionDirect(1)
    diag_data%var(iLookDIAG%scalarCanopySunlitFraction)%dat(1) = diag_data_device%scalarCanopySunlitFraction(1)
    diag_data%var(iLookDIAG%scalarCanopySunlitLAI)%dat(1) = diag_data_device%scalarCanopySunlitLAI(1)
    diag_data%var(iLookDIAG%scalarCanopyShadedLAI)%dat(1) = diag_data_device%scalarCanopyShadedLAI(1)
    diag_data%var(iLookDIAG%scalarGroundAlbedo)%dat(1) = diag_data_device%scalarGroundAlbedo(1)
    diag_data%var(iLookDIAG%scalarLatHeatSubVapCanopy)%dat(1) = diag_data_device%scalarLatHeatSubVapCanopy(1)
    diag_data%var(iLookDIAG%scalarLatHeatSubVapGround)%dat(1) = diag_data_device%scalarLatHeatSubVapGround(1)
    diag_data%var(iLookDIAG%scalarSatVP_CanopyTemp)%dat(1) = diag_data_device%scalarSatVP_CanopyTemp(1)
    diag_data%var(iLookDIAG%scalarSatVP_GroundTemp)%dat(1) = diag_data_device%scalarSatVP_GroundTemp(1)
    diag_data%var(iLookDIAG%scalarZ0Canopy)%dat(1) = diag_data_device%scalarZ0Canopy(1)
    diag_data%var(iLookDIAG%scalarWindReductionFactor)%dat(1) = diag_data_device%scalarWindReductionFactor(1)
    diag_data%var(iLookDIAG%scalarZeroPlaneDisplacement)%dat(1) = diag_data_device%scalarZeroPlaneDisplacement(1)
    diag_data%var(iLookDIAG%scalarRiBulkCanopy)%dat(1) = diag_data_device%scalarRiBulkCanopy(1)
    diag_data%var(iLookDIAG%scalarRiBulkGround)%dat(1) = diag_data_device%scalarRiBulkGround(1)
    diag_data%var(iLookDIAG%scalarCanopyStabilityCorrection)%dat(1) = diag_data_device%scalarCanopyStabilityCorrection(1)
    diag_data%var(iLookDIAG%scalarGroundStabilityCorrection)%dat(1) = diag_data_device%scalarGroundStabilityCorrection(1)
    diag_data%var(iLookDIAG%scalarIntercellularCO2Sunlit)%dat(1) = diag_data_device%scalarIntercellularCO2Sunlit(1)
    diag_data%var(iLookDIAG%scalarIntercellularCO2Shaded)%dat(1) = diag_data_device%scalarIntercellularCO2Shaded(1)
    diag_data%var(iLookDIAG%scalarTranspireLim)%dat(1) = diag_data_device%scalarTranspireLim(1)
    diag_data%var(iLookDIAG%scalarTranspireLimAqfr)%dat(1) = diag_data_device%scalarTranspireLimAqfr(1)
    diag_data%var(iLookDIAG%scalarFoliageNitrogenFactor)%dat(1) = diag_data_device%scalarFoliageNitrogenFactor(1)
    diag_data%var(iLookDIAG%scalarSoilRelHumidity)%dat(1) = diag_data_device%scalarSoilRelHumidity(1)
    diag_data%var(iLookDIAG%scalarAquiferRootFrac)%dat(1) = diag_data_device%scalarAquiferRootFrac(1)
    diag_data%var(iLookDIAG%scalarFracLiqVeg)%dat(1) = diag_data_device%scalarFracLiqVeg(1)
    diag_data%var(iLookDIAG%scalarCanopyWetFraction)%dat(1) = diag_data_device%scalarCanopyWetFraction(1)
    diag_data%var(iLookDIAG%scalarSnowAge)%dat(1) = diag_data_device%scalarSnowAge(1)
    diag_data%var(iLookDIAG%scalarGroundSnowFraction)%dat(1) = diag_data_device%scalarGroundSnowFraction(1)
    diag_data%var(iLookDIAG%scalarInfilArea)%dat(1) = diag_data_device%scalarInfilArea(1)
    diag_data%var(iLookDIAG%scalarFrozenArea)%dat(1) = diag_data_device%scalarFrozenArea(1)
    diag_data%var(iLookDIAG%scalarSoilControl)%dat(1) = diag_data_device%scalarSoilControl(1)
    diag_data%var(iLookDIAG%scalarSoilCompress)%dat(1) = diag_data_device%scalarSoilCompress(1)
    diag_data%var(iLookDIAG%scalarTotalSoilLiq)%dat(1) = diag_data_device%scalarTotalSoilLiq(1)
    diag_data%var(iLookDIAG%scalarTotalSoilIce)%dat(1) = diag_data_device%scalarTotalSoilIce(1)
    diag_data%var(iLookDIAG%scalarTotalSoilWat)%dat(1) = diag_data_device%scalarTotalSoilWat(1)
    diag_data%var(iLookDIAG%scalarKappa)%dat(1) = diag_data_device%scalarKappa(1)
    diag_data%var(iLookDIAG%scalarVolLatHt_fus)%dat(1) = diag_data_device%scalarVolLatHt_fus(1)
    diag_data%var(iLookDIAG%balanceCasNrg)%dat(1) = diag_data_device%balanceCasNrg(1)
    diag_data%var(iLookDIAG%balanceVegNrg)%dat(1) = diag_data_device%balanceVegNrg(1)
    diag_data%var(iLookDIAG%balanceSnowNrg)%dat(1) = diag_data_device%balanceSnowNrg(1)
    diag_data%var(iLookDIAG%balanceSoilNrg)%dat(1) = diag_data_device%balanceSoilNrg(1)
    diag_data%var(iLookDIAG%balanceVegMass)%dat(1) = diag_data_device%balanceVegMass(1)
    diag_data%var(iLookDIAG%balanceSnowMass)%dat(1) = diag_data_device%balanceSnowMass(1)
    diag_data%var(iLookDIAG%balanceSoilMass)%dat(1) = diag_data_device%balanceSoilMass(1)
    diag_data%var(iLookDIAG%balanceAqMass)%dat(1) = diag_data_device%balanceAqMass(1)
    diag_data%var(iLookDIAG%mLayerEnthTemp)%dat = diag_data_device%mLayerEnthTemp(:,1)
    diag_data%var(iLookDIAG%mLayerVolHtCapBulk)%dat = diag_data_device%mLayerVolHtCapBulk_m(:,1)
    diag_data%var(iLookDIAG%mLayerCm)%dat = diag_data_device%mLayerCm_m(:,1)
    diag_data%var(iLookDIAG%mLayerThermalC)%dat = diag_data_device%mLayerThermalC_m(:,1)
    diag_data%var(iLookDIAG%iLayerThermalC)%dat = diag_data_device%iLayerThermalC_m(:,1)
    diag_data%var(iLookDIAG%spectralAlbGndDirect)%dat = diag_data_device%spectralAlbGndDirect(:,1)
    diag_data%var(iLookDIAG%spectralAlbGndDiffuse)%dat = diag_data_device%spectralAlbGndDiffuse(:,1)
    diag_data%var(iLookDIAG%mLayerTranspireLim)%dat = diag_data_device%mLayerTranspireLim_m(:,1)
    diag_data%var(iLookDIAG%mLayerRootDensity)%dat = diag_data_device%mLayerRootDensity_m(:,1)
    diag_data%var(iLookDIAG%spectralSnowAlbedoDirect)%dat = diag_data_device%spectralSnowAlbedoDirect(:,1)
    if (nSnow.ne.0) diag_data%var(iLookDIAG%mLayerFracLiqSnow)%dat = diag_data_device%mLayerFracLiqSnow_m(:,1)
    if (nSnow.ne.0) diag_data%var(iLookDIAG%mLayerThetaResid)%dat = diag_data_device%mLayerThetaResid_m(:,1)
    if (nSnow.ne.0) diag_data%var(iLookDIAG%mLayerPoreSpace)%dat = diag_data_device%mLayerPoreSpace_m(:,1)
    diag_data%var(iLookDIAG%mLayerMeltFreeze)%dat = diag_data_device%mLayerMeltFreeze(:,1)
    diag_data%var(iLookDIAG%mLayerVolFracAir)%dat = diag_data_device%mLayerVolFracAir_m(:,1)
    diag_data%var(iLookDIAG%mLayerTcrit)%dat = diag_data_device%mLayerTcrit(:,1)
    diag_data%var(iLookDIAG%mLayerCompress)%dat = diag_data_device%mLayerCompress_m(:,1)
    diag_data%var(iLookDIAG%mLayerMatricHeadLiq)%dat = diag_data_device%mLayerMatricHeadLiq(:,1)
    diag_data%var(iLookDIAG%scalarVGn_m)%dat = diag_data_device%scalarVGn_m_m(:,1)
    diag_data%var(iLookDIAG%balanceLayerNrg)%dat = diag_data_device%balanceLayerNrg(:,1)
    diag_data%var(iLookDIAG%balanceLayerMass)%dat   = diag_data_device%balanceLayerMass(:,1)

  end subroutine finalize_device_diag_data

  subroutine deallocate_device_diag_data(diag_data_device,nSnow)
    type(diag_data_device), intent(inout) :: diag_data_device
    integer(i4b),intent(in) :: nSnow
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
      

      subroutine allocate_device_bvar_data(bvar_data_device, bvar_data)
        type(bvar_data_device),intent(inout) :: bvar_data_device
        type(var_dlength),intent(in) :: bvar_data
        
        bvar_data_device%basin__AquiferStorage = bvar_data%var(iLookBVAR%basin__AquiferStorage)%dat(1)
      end subroutine allocate_device_bvar_data

      subroutine deallocate_device_bvar_data(bvar_data_device)
        type(bvar_data_device),intent(inout) :: bvar_data_device
        
        deallocate(bvar_data_device%basin__AquiferStorage)
      end subroutine deallocate_device_bvar_data

  subroutine allocate_device_forc_data(forc_data_device, forc_data,nGRU)
    type(forc_data_device),intent(inout) :: forc_data_device
    type(var_d),intent(in) :: forc_data
    integer(i4b) :: nGRU

    allocate(forc_data_device%time(nGRU))
    allocate(forc_data_device%pptrate(nGRU))
    allocate(forc_data_device%SWRadAtm(nGRU))
    allocate(forc_data_device%LWRadAtm_d(nGRU))
    allocate(forc_data_device%airtemp_d(nGRU))
    allocate(forc_data_device%windspd_d(nGRU))
    allocate(forc_data_device%airpres_d(nGRU))
    allocate(forc_data_device%spechum(nGRU))

    forc_data_device%time = forc_data%var(iLookFORCE%time)
    forc_data_device%pptrate = forc_data%var(iLookFORCE%pptrate)
    forc_data_device%SWRadAtm = forc_data%var(iLookFORCE%SWRadAtm)
    forc_data_device%LWRadAtm_d = forc_data%var(iLookFORCE%LWRadAtm)
    forc_data_device%airtemp_d = forc_data%var(iLookFORCE%airtemp)
    forc_data_device%windspd_d = forc_data%var(iLookFORCE%windspd)
    forc_data_device%airpres_d = forc_data%var(iLookFORCE%airpres)
    forc_data_device%spechum = forc_data%var(iLookFORCE%spechum)

  end subroutine allocate_device_forc_data
    
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
        type(var_i) :: type_data
        allocate(type_data_d%vegTypeIndex(nGRU))
        allocate(type_data_d%soilTypeIndex(nGRU))
        allocate(type_data_d%slopeTypeIndex(nGRU))
        allocate(type_data_d%downHRUindex(nGRU))

        type_data_d%vegTypeIndex = type_data%var(iLookTYPE%vegTypeIndex)
        type_data_d%soilTypeIndex = type_data%var(iLookTYPE%soilTypeIndex)
        type_data_d%slopeTypeIndex = type_data%var(iLookTYPE%slopeTypeIndex)
        type_data_d%downHRUindex = type_data%var(iLookTYPE%downHRUindex)
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
veg_param%saim = saim
veg_param%laim = laim
veg_param%dveg = dveg
veg_param%tmin = tmin

      end subroutine
subroutine allocate_device_attr_data(attr_data_device,attr_data,nGRU)
  type(attr_data_device) :: attr_data_device
  type(var_d) :: attr_data
  integer(i4b) :: nGRU
  allocate(attr_data_device%latitude(nGRU))
  allocate(attr_data_device%longitude(nGRU))
  allocate(attr_data_device%elevation(nGRU))
  allocate(attr_data_device%tan_slope(nGRU))
  allocate(attr_data_device%contourLength(nGRU))
  allocate(attr_data_device%HRUarea(nGRU))
  allocate(attr_data_device%mHeight(nGRU))
  allocate(attr_data_device%aspect(nGRU))
  attr_data_device%latitude = attr_data%var(iLookATTR%latitude)
  attr_data_device%longitude = attr_data%var(iLookATTR%longitude)
  attr_data_device%elevation = attr_data%var(iLookATTR%elevation)
  attr_data_device%tan_slope = attr_data%var(iLookATTR%tan_slope)
  attr_data_device%contourLength = attr_data%var(iLookATTR%contourLength)
  attr_data_device%HRUarea = attr_data%var(iLookATTR%HRUarea)
  attr_data_device%mHeight = attr_data%var(iLookATTR%mHeight)
  attr_data_device%aspect = attr_data%var(iLookATTR%aspect)
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
    end module initialize_device      