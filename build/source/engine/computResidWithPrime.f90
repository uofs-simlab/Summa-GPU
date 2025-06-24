

module computResidWithPrime_module

  ! data types
  USE nrtype
  
  ! derived types to define the data structures
  USE data_types,only:&
                      var_ilength,  & ! data vector with variable length dimension (i4b)
                      var_dlength     ! data vector with variable length dimension (rkind)
  
  ! named variables
  USE var_lookup,only:iLookPROG       ! named variables for structure elements
  USE var_lookup,only:iLookDIAG       ! named variables for structure elements
  USE var_lookup,only:iLookFLUX       ! named variables for structure elements
  USE var_lookup,only:iLookINDEX      ! named variables for structure elements
  
  ! access the global print flag
  USE globalData,only:globalPrintFlag
  
  ! access missing values
  USE globalData,only:integerMissing  ! missing integer
  USE globalData,only:realMissing     ! missing real number
  
  ! define access to state variables to print
  USE globalData,only: iJac1          ! first layer of the Jacobian to print
  USE globalData,only: iJac2          ! last layer of the Jacobian to print
  
  ! domain types
  USE globalData,only:iname_veg       ! named variables for vegetation
  USE globalData,only:iname_snow      ! named variables for snow
  USE globalData,only:iname_soil      ! named variables for soil
  
  ! named variables to describe the state variable type
  USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
  USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
  USE globalData,only:iname_watCanopy ! named variable defining the mass of water on the vegetation canopy
  USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
  USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
  USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
  USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
  USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers
  
  ! constants
  USE multiconst,only:&
                      LH_fus,       & ! latent heat of fusion                (J kg-1)
                      iden_ice,     & ! intrinsic density of ice             (kg m-3)
                      iden_water      ! intrinsic density of liquid water    (kg m-3)
  ! privacy
  implicit none
  public::computResidWithPrime
  contains
  
  ! **********************************************************************************************************
  ! public subroutine computResidWithPrime: compute the residual vector
  ! **********************************************************************************************************
  subroutine computResidWithPrime(&
                        ! input: model control
                        dt,                        & ! intent(in):  length of the time step (seconds)
                        nSnow,                     & ! intent(in):  number of snow layers
                        nSoil,                     & ! intent(in):  number of soil layers
                        nLayers,                   & ! intent(in):  total number of layers
                        nGRU, &
                        enthalpyStateVec,               & ! intent(in):  flag if enthalpy is state variable
                        ! input: flux vectors
                        sMul,                      & ! intent(in):  state vector multiplier (used in the residual calculations)
                        fVec,                      & ! intent(in):  flux vector
                        ! input: state variables (already disaggregated into scalars and vectors)
                        scalarCanairTempPrime,     & ! intent(in):  prime value for the temperature of the canopy air space (K s-1)
                        scalarCanopyTempPrime,     & ! intent(in):  prime value for the temperature of the vegetation canopy (K s-1)
                        scalarCanopyWatPrime,      & ! intent(in):  prime value for the water on the vegetation canopy (kg m-2 s-1)
                        mLayerTempPrime,           & ! intent(in):  prime vector of the temperature of each snow and soil layer (K s-1)
                        scalarAquiferStoragePrime, & ! intent(in):  prime value for storage of water in the aquifer (m s-1)
                        ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
                        scalarCanopyIcePrime,      & ! intent(in):  prime value for the ice on the vegetation canopy (kg m-2 s-1)
                        scalarCanopyLiqPrime,      & ! intent(in):  prime value for the liq on the vegetation canopy (kg m-2 s-1)
                        mLayerVolFracIcePrime,     & ! intent(in):  prime vector of the volumetric ice in each snow and soil layer (s-1)
                        mLayerVolFracWatPrime,     & ! intent(in):  prime vector of the volumetric water in each snow and soil layer (s-1)
                        mLayerVolFracLiqPrime,     & ! intent(in):  prime vector of the volumetric liq in each snow and soil layer (s-1)
                        ! input: enthalpy terms
                        scalarCanopyCm_noLHTrial,  & ! intent(in):  Cm without latent heat part for vegetation canopy (J kg K-1)
                        mLayerCm_noLHTrial,        & ! intent(in):  Cm without latent heat part for each snow and soil layer (J kg K-1)
                        scalarCanairEnthalpyPrime, & ! intent(in):  prime value for the enthalpy of the canopy air space (W m-3)
                        scalarCanopyEnthalpyPrime, & ! intent(in):  prime value for the of enthalpy of the vegetation canopy (W m-3)
                        mLayerEnthalpyPrime,       & ! intent(in):  prime vector of the of enthalpy of each snow and soil layer (W m-3)
                        ! input: data structures
                        prog_data,                 & ! intent(in):  model prognostic variables for a local HRU
                        diag_data,                 & ! intent(in):  model diagnostic variables for a local HRU
                        flux_data,                 & ! intent(in):  model fluxes for a local HRU
                        indx_data,                 & ! intent(in):  index data
                        ! output
                        rAdd,                      & ! intent(out): additional (sink) terms on the RHS of the state equation
                        rVec,                      & ! intent(out): residual vector
                        err,message)                 ! intent(out): error control
                        use ieee_arithmetic
                        use device_data_types
    ! --------------------------------------------------------------------------------------------------------------------------------
    implicit none
    ! input: model control
    real(rkind),intent(in)          :: dt                        ! length of the time step (seconds)
    integer(i4b),intent(in),device         :: nSnow(:)                     ! number of snow layers
    integer(i4b),intent(in)         :: nSoil                     ! number of soil layers
    integer(i4b),intent(in),device         :: nLayers(:)                   ! total number of layers in the snow+soil domain
    integer(i4b),intent(in) :: nGRU
    logical(lgt),intent(in)         :: enthalpyStateVec          ! flag if enthalpy is state variable
    ! input: flux vectors
    real(qp),intent(in),device             :: sMul(:,:)   ! NOTE: qp      ! state vector multiplier (used in the residual calculations)
    real(rkind),intent(in),device          :: fVec(:,:)                   ! flux vector
    ! input: state variables (already disaggregated into scalars and vectors)
    real(rkind),intent(in),device          :: scalarCanairTempPrime(:)     ! prime value for temperature of the canopy air space (K s-1)
    real(rkind),intent(in),device          :: scalarCanopyTempPrime(:)     ! prime value for temperature of the vegetation canopy (K s-1)
    real(rkind),intent(in),device          :: scalarCanopyWatPrime(:)      ! prime value for canopy total water content (kg m-2 s-1)
    real(rkind),intent(in),device          :: mLayerTempPrime(:,:)        ! prime vector of temperature of each snow/soil layer (K s-1) content
    real(rkind),intent(in),device          :: scalarAquiferStoragePrime(:) ! prime value of aquifer storage (m s-1)
    ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
    real(rkind),intent(in),device          :: scalarCanopyIcePrime(:)      ! prime value for mass of ice on the vegetation canopy (kg m-2 s-1)
    real(rkind),intent(in),device          :: scalarCanopyLiqPrime(:)      ! prime value for the liq on the vegetation canopy (kg m-2 s-1)
    real(rkind),intent(in),device          :: mLayerVolFracIcePrime(:,:)  ! prime vector of volumetric fraction of ice (s-1)
    real(rkind),intent(in),device          :: mLayerVolFracWatPrime(:,:)  ! prime vector of the volumetric water in each snow and soil layer (s-1)
    real(rkind),intent(in),device          :: mLayerVolFracLiqPrime(:,:)  ! prime vector of the volumetric water in each snow and soil layer (s-1)
    ! input: enthalpy terms
    real(rkind),intent(in),device          :: scalarCanopyCm_noLHTrial(:)  ! Cm without latent heat part for vegetation canopy (-)
    real(rkind),intent(in),device          :: mLayerCm_noLHTrial(:,:)     ! Cm without latent heat part for each snow and soil layer (-)
    real(rkind),intent(in),device          :: scalarCanairEnthalpyPrime(:) ! prime value for enthalpy of the canopy air space (W m-3)
    real(rkind),intent(in),device          :: scalarCanopyEnthalpyPrime(:) ! prime value for enthalpy of the vegetation canopy (W m-3)
    real(rkind),intent(in),device          :: mLayerEnthalpyPrime(:,:)    ! prime vector of enthalpy of each snow and soil layer (W m-3)
    ! input: data structures
    type(prog_data_device),intent(in)    :: prog_data                 ! prognostic variables for a local HRU
    type(diag_data_device),intent(in)    :: diag_data                 ! diagnostic variables for a local HRU
    type(flux_data_device),intent(in)    :: flux_data                 ! model fluxes for a local HRU
    type(indx_data_device),intent(in)    :: indx_data                 ! indices defining model states and layers
    ! output
    real(rkind),intent(out),device         :: rAdd(:,:)                   ! additional (sink) terms on the RHS of the state equation
    real(qp),intent(out),device            :: rVec(:,:)   ! NOTE: qp      ! residual vector
    integer(i4b),intent(out)        :: err                       ! error code
    character(*),intent(out)        :: message                   ! error message
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! local variables
    ! --------------------------------------------------------------------------------------------------------------------------------
    integer(i4b)                     :: iLayer                   ! index of layer within the snow+soil domain
    integer(i4b),parameter           :: ixVegVolume=1            ! index of the desired vegetation control volumne (currently only one veg layer)
    real(rkind)                      :: scalarCanopyHydPrime     ! trial value for canopy water (kg m-2), either liquid water content or total water content
    real(rkind)   :: mLayerVolFracHydPrime    ! vector of volumetric water content (-), either liquid water content or total water content
    integer(i4b) :: iGRU
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! link to the necessary variables for the residual computations
    associate(&
      ! canopy and layer depth
      canopyDepth             => diag_data%scalarCanopyDepth      ,& ! intent(in): [dp]      canopy depth (m)
      mLayerDepth             => prog_data%mLayerDepth               ,& ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
      ! model fluxes (sink terms in the soil domain)
      mLayerTranspire         => flux_data%mLayerTranspire_m           ,& ! intent(in): [dp]     transpiration loss from each soil layer (m s-1)
      mLayerBaseflow          => flux_data%mLayerBaseflow_m            ,& ! intent(in): [dp(:)]  baseflow from each soil layer (m s-1)
      mLayerCompress          => diag_data%mLayerCompress_m            ,& ! intent(in): [dp(:)]  change in storage associated with compression of the soil matrix (-)
      ! number of state variables of a specific type
      nSnowSoilNrg            => indx_data%nLayers_d         ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
      nSnowSoilHyd            => indx_data%nLayers_d         ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
      nSoilOnlyHyd            => indx_data%nSoil         ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
      ! model indices
      ixCasNrg                => indx_data%ixCasNrg              ,& ! intent(in): [i4b]    index of canopy air space energy state variable
      ixVegNrg                => indx_data%ixVegNrg              ,& ! intent(in): [i4b]    index of canopy energy state variable
      ixVegHyd                => indx_data%ixVegHyd              ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
      ixAqWat                 => indx_data%ixAqWat               ,& ! intent(in): [i4b]    index of water storage in the aquifer
      ixSnowSoilNrg_m           => indx_data%ixSnowSoilNrg            ,& ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
      ixSnowSoilHyd           => indx_data%ixSnowSoilHyd            ,& ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
      ixSoilOnlyHyd_m           => indx_data%ixSoilOnlyHyd            ,& ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
      ixStateType_m             => indx_data%ixStateType              ,& ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
      ixHydCanopy_m             => indx_data%ixHydCanopy              ,& ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
      ixHydType_m               => indx_data%ixHydType                ,& ! intent(in): [i4b(:)] named variables defining the type of hydrology states in snow+soil domain
      layerType_m               => indx_data%layerType                 & ! intent(in): [i4b(:)] named variables defining the type of layer in snow+soil domain
      ) ! association to necessary variables for the residual computations
      ! --------------------------------------------------------------------------------------------------------------------------------
      ! initialize error control
      err=0; message="computResidWithPrime/"
  
      ! ---
      ! * compute sink terms...
      ! -----------------------
  
      ! intialize additional terms on the RHS as zero
      rAdd = 0._rkind
  
      ! add melt freeze terms only if not using enthalpy terms 
      ! NOTE: would need to use these if were using enthTemp terms
      if(.not.enthalpyStateVec)then
        ! compute energy associated with melt freeze for the vegetation canopy
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
        if(ixVegNrg(iGRU)/=integerMissing) rAdd(ixVegNrg(iGRU),iGRU) = rAdd(ixVegNrg(iGRU),iGRU) + LH_fus*scalarCanopyIcePrime(iGRU)/canopyDepth(iGRU)   ! energy associated with melt/freeze (J m-3)
        end do
   
        ! compute energy associated with melt/freeze for snow
        ! NOTE: allow expansion of ice during melt-freeze for snow; deny expansion of ice during melt-freeze for soil
        ! if(nSnowSoilNrg>0)then
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
          do iLayer=1,nLayers(iGRU)
            if(ixSnowSoilNrg_m(iLayer,iGRU)/=integerMissing) then   ! (loop through non-missing energy state variables in the snow+soil domain)
            select case( layerType_m(iLayer,iGRU) )
              case(iname_snow); rAdd( ixSnowSoilNrg_m(iLayer,iGRU),iGRU ) = rAdd( ixSnowSoilNrg_m(iLayer,iGRU),iGRU ) + LH_fus*iden_ice * mLayerVolFracIcePrime(iLayer,iGRU)
              case(iname_soil); rAdd( ixSnowSoilNrg_m(iLayer,iGRU),iGRU ) = rAdd( ixSnowSoilNrg_m(iLayer,iGRU),iGRU ) + LH_fus*iden_water * mLayerVolFracIcePrime(iLayer,iGRU)
            end select
          end if
          end do  ! looping through non-missing energy state variables in the snow+soil domain
        end do
        ! endif
  
      endif
  
      ! sink terms soil hydrology (-)
      ! NOTE 1: state variable is volumetric water content, so melt-freeze is not included
      ! NOTE 2: ground evaporation was already included in the flux at the upper boundary
      ! NOTE 3: rAdd(ixSnowOnlyWat)=0, and is defined in the initialization above
      ! NOTE 4: same sink terms for matric head and liquid matric potential
      ! if(nSoilOnlyHyd>0)then
        !$cuf kernel do(2) <<<*,*>>>
        do iGRU=1,nGRU
        do iLayer=1,nSoil
          if(ixSoilOnlyHyd_m(iLayer,iGRU)/=integerMissing) then   ! (loop through non-missing hydrology state variables in the snow+soil domain)
         rAdd( ixSoilOnlyHyd_m(iLayer,iGRU) ,iGRU) = rAdd( ixSoilOnlyHyd_m(iLayer,iGRU),iGRU ) + ( ( mLayerTranspire(iLayer,iGRU) - mLayerBaseflow(iLayer,iGRU) )/mLayerDepth(iLayer+nSnow(iGRU),iGRU) - mLayerCompress(iLayer,iGRU) )*dt
          end if
        end do  ! looping through non-missing energy state variables in the snow+soil domain
      end do
      ! endif
  
      ! ---
      ! * compute the residual vector...
      ! --------------------------------
  
      ! compute the residual vector for the vegetation canopy
      ! NOTE: sMul(ixVegHyd) = 1, but include as it converts all variables to quadruple precision
      ! --> energy balance
      if(enthalpyStateVec)then
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
        if(ixCasNrg(iGRU)/=integerMissing) rVec(ixCasNrg(iGRU),iGRU) = scalarCanairEnthalpyPrime(iGRU) - ( fVec(ixCasNrg(iGRU),iGRU)*dt + rAdd(ixCasNrg(iGRU),iGRU) )
        if(ixVegNrg(iGRU)/=integerMissing) rVec(ixVegNrg(iGRU),iGRU) = scalarCanopyEnthalpyPrime(iGRU) - ( fVec(ixVegNrg(iGRU),iGRU)*dt + rAdd(ixVegNrg(iGRU),iGRU) )
        end do
      else
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
        if(ixCasNrg(iGRU)/=integerMissing) rVec(ixCasNrg(iGRU),iGRU) = sMul(ixCasNrg(iGRU),iGRU) * scalarCanairTempPrime(iGRU) - ( fVec(ixCasNrg(iGRU),iGRU)*dt + rAdd(ixCasNrg(iGRU),iGRU) )
        if(ixVegNrg(iGRU)/=integerMissing) rVec(ixVegNrg(iGRU),iGRU) = sMul(ixVegNrg(iGRU),iGRU) * scalarCanopyTempPrime(iGRU) + scalarCanopyCm_noLHTrial(iGRU) * scalarCanopyWatPrime(iGRU)/canopyDepth(iGRU) &
                                                     - ( fVec(ixVegNrg(iGRU),iGRU)*dt + rAdd(ixVegNrg(iGRU),iGRU) )
        end do
      endif   
      !$cuf kernel do(1) <<<*,*>>>                                         
      do iGRU=1,nGRU
      ! --> mass balance
      if(ixVegHyd(iGRU)/=integerMissing)then
        scalarCanopyHydPrime = merge(scalarCanopyWatPrime(iGRU), scalarCanopyLiqPrime(iGRU), (ixStateType_m( ixHydCanopy_m(iGRU),iGRU )==iname_watCanopy) )
        rVec(ixVegHyd(iGRU),iGRU) = sMul(ixVegHyd(iGRU),iGRU)*scalarCanopyHydPrime - ( fVec(ixVegHyd(iGRU),iGRU)*dt + rAdd(ixVegHyd(iGRU),iGRU) )
      endif
      ! compute the residual vector for the aquifer
      if(ixAqWat(iGRU)/=integerMissing)  rVec(ixAqWat(iGRU),iGRU) = sMul(ixAqWat(iGRU),iGRU)*scalarAquiferStoragePrime(iGRU) - ( fVec(ixAqWat(iGRU),iGRU)*dt + rAdd(ixAqWat(iGRU),iGRU) )
    end do
  
      ! compute the residual vector for the snow and soil sub-domains for energy
      ! if(nSnowSoilNrg>0)then
        if (enthalpyStateVec) then
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
            do iLayer=1,nLayers(iGRU)
              if(ixSnowSoilNrg_m(iLayer,iGRU)/=integerMissing) then   ! (loop through non-missing energy state variables in the snow+soil domain)
                rVec( ixSnowSoilNrg_m(iLayer,iGRU),iGRU ) = mLayerEnthalpyPrime(iLayer,iGRU) - ( fVec( ixSnowSoilNrg_m(iLayer,iGRU) ,iGRU)*dt + rAdd( ixSnowSoilNrg_m(iLayer,iGRU) ,iGRU) )
              end if
            end do  ! looping through non-missing energy state variables in the snow+soil domain
          end do
        else
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
            do iLayer=1,nLayers(iGRU)
              if(ixSnowSoilNrg_m(iLayer,iGRU)/=integerMissing) then   ! (loop through non-missing energy state variables in the snow+soil domain)
                rVec( ixSnowSoilNrg_m(iLayer,iGRU) ,iGRU) = sMul( ixSnowSoilNrg_m(iLayer,iGRU) ,iGRU) * mLayerTempPrime(iLayer,iGRU) + mLayerCm_noLHTrial(iLayer,iGRU) * mLayerVolFracWatPrime(iLayer,iGRU) &
                                           - ( fVec( ixSnowSoilNrg_m(iLayer,iGRU),iGRU )*dt + rAdd( ixSnowSoilNrg_m(iLayer,iGRU) ,iGRU) )
              endif
            end do  ! looping through non-missing energy state variables in the snow+soil domain
          end do
        end if
  
      ! endif
  
      ! compute the residual vector for the snow and soil sub-domains for hydrology
      ! NOTE: residual depends on choice of state variable
      ! if(nSnowSoilHyd>0)then
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
          do  iLayer=1,nLayers(iGRU)
            if (ixSnowSoilHyd(iLayer,iGRU)/=integerMissing) then   ! (loop through non-missing hydrology state variables in the snow+soil domain)
              ! (get the correct state variable)
              mLayerVolFracHydPrime = merge(mLayerVolFracWatPrime(iLayer,iGRU), mLayerVolFracLiqPrime(iLayer,iGRU), (ixHydType_m(iLayer,iGRU)==iname_watLayer .or. ixHydType_m(iLayer,iGRU)==iname_matLayer) )
              ! (compute the residual)
              rVec( ixSnowSoilHyd(iLayer,iGRU) ,iGRU) = mLayerVolFracHydPrime - ( fVec( ixSnowSoilHyd(iLayer,iGRU),iGRU )*dt + rAdd( ixSnowSoilHyd(iLayer,iGRU),iGRU ) )
            end if
          end do  ! looping through non-missing energy state variables in the snow+soil domain
        end do
      ! endif
  
      ! print result
      ! if(globalPrintFlag)then
      !   write(*,'(a,1x,100(e12.5,1x))') 'rVec = ', rVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
      !   write(*,'(a,1x,100(e12.5,1x))') 'fVec = ', fVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
      ! endif
  
      ! check
      ! if(any(ieee_is_Nan(rVec)))then
        ! message=trim(message)//'vector of residuals contains NaN value(s) ' ! formerly known as the Indian bread error
        ! write(*,'(a,1x,100(e12.5,1x))') 'rVec = ', rVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
        ! write(*,'(a,1x,100(e12.5,1x))') 'fVec = ', fVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
        ! err=20; return
      ! endif
      
    end associate
  
  end subroutine computResidWithPrime
  
  end module computResidWithPrime_module
  