! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
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

module updateVarsWithPrime_module

  ! data types
  USE nrtype
  
  ! missing values
  USE globalData,only:integerMissing  ! missing integer
  USE globalData,only:realMissing     ! missing real number
  
  ! access the global print flag
  USE globalData,only:globalPrintFlag
  
  ! domain types
  USE globalData,only:iname_cas       ! named variables for canopy air space
  USE globalData,only:iname_veg       ! named variables for vegetation canopy
  USE globalData,only:iname_snow      ! named variables for snow
  USE globalData,only:iname_soil      ! named variables for soil
  USE globalData,only:iname_aquifer   ! named variables for the aquifer
  
  ! named variables to describe the state variable type
  USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
  USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
  USE globalData,only:iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
  USE globalData,only:iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
  USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
  USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
  USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
  USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
  USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers
  
  ! metadata for information in the data structures
  USE globalData,only:indx_meta       ! metadata for the variables in the index structure
  
  ! constants
  USE multiconst,only:&
                      Tfreeze,      & ! temperature at freezing              (K)
                      LH_fus,       & ! latent heat of fusion                (J kg-1)
                      iden_ice,     & ! intrinsic density of ice             (kg m-3)
                      iden_water      ! intrinsic density of liquid water    (kg m-3)
  
  ! provide access to the derived types to define the data structures
  USE data_types,only:&
                      var_i,        & ! data vector (i4b)
                      var_d,        & ! data vector (rkind)
                      var_ilength,  & ! data vector with variable length dimension (i4b)
                      zLookup,      & ! data vector with variable length dimension (rkind)
                      var_dlength     ! data vector with variable length dimension (rkind)
  
  ! provide access to indices that define elements of the data structures
  USE var_lookup,only:iLookDIAG             ! named variables for structure elements
  USE var_lookup,only:iLookPROG             ! named variables for structure elements
  USE var_lookup,only:iLookDERIV            ! named variables for structure elements
  USE var_lookup,only:iLookPARAM            ! named variables for structure elements
  USE var_lookup,only:iLookINDEX            ! named variables for structure elements
  
  ! provide access to routines to update states
  USE updatStateWithPrime_module,only:updateSnowPrime     ! update snow states
  USE updatStateWithPrime_module,only:updateSoilPrime     ! update soil states
  
  ! provide access to functions for the constitutive functions and derivatives
  USE snow_utils_module,only:fracliquid              ! compute the fraction of liquid water (snow)
  USE snow_utils_module,only:dFracLiq_dTk            ! differentiate the freezing curve w.r.t. temperature (snow)
  USE soil_utils_module,only:dTheta_dTk              ! differentiate the freezing curve w.r.t. temperature (soil)
  USE soil_utils_module,only:dTheta_dPsi             ! derivative in the soil water characteristic (soil)
  USE soil_utils_module,only:dPsi_dTheta             ! derivative in the soil water characteristic (soil)
  USE soil_utils_module,only:matricHead              ! compute the matric head based on volumetric water content
  USE soil_utils_module,only:volFracLiq              ! compute volumetric fraction of liquid water
  USE soil_utils_module,only:crit_soilT              ! compute critical temperature below which ice exists
  USE soil_utilsAddPrime_module,only:liquidHeadPrime ! compute the liquid water matric potential
  USE soil_utilsAddPrime_module,only:d2Theta_dPsi2   ! second derivative in the soil water characteristic (soil)
  USE soil_utilsAddPrime_module,only:d2Theta_dTk2    ! second derivative in the freezing curve w.r.t. temperature (soil)
  USE enthalpyTemp_module,only:enthalpy2T_cas        ! compute canopy air space temperature from enthalpy
  USE enthalpyTemp_module,only:enthalpy2T_veg        ! compute canopy temperature from enthalpy and water content
  USE enthalpyTemp_module,only:enthalpy2T_snow       ! compute snow layer temperature from enthalpy and water content
  USE enthalpyTemp_module,only:enthalpy2T_soil       ! compute soil layer temperature from enthalpy and matric potential
  
  ! IEEE checks
  USE, intrinsic :: ieee_arithmetic            ! check values (NaN, etc.)
  USE mDecisions_module,only:  &
   closedForm,                 & ! use temperature with closed form heat capacity
   enthalpyFormLU,             & ! use enthalpy with soil temperature-enthalpy lookup tables
   enthalpyForm                  ! use enthalpy with soil temperature-enthalpy analytical solution
  
  implicit none
  private
  public::updateVarsWithPrime
  
  contains
  
  ! **********************************************************************************************************
  ! public subroutine updateVarsWithPrime: compute diagnostic variables and derivatives for Prime Jacobian
  ! **********************************************************************************************************
  subroutine updateVarsWithPrime(&
                       ! input
                       ixNrgConserv,                          & ! intent(in):    flag if enthalpy is the state variable
                       use_lookup,                                & ! intent(in):    flag to use the lookup table for soil enthalpy
                       computJac,                                 & ! intent(in):    flag if computing for Jacobian update
                       do_adjustTemp,                             & ! intent(in):    flag to adjust temperature to account for the energy used in melt+freeze
                       nGRU, &
                       mpar_data,                                 & ! intent(in):    model parameters for a local HRU
                       indx_data,                                 & ! intent(in):    indices defining model states and layers
                       prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,                                 & ! intent(inout): model diagnostic variables for a local HRU
                       deriv_data,                                & ! intent(inout): derivatives in model fluxes w.r.t. relevant state variables
                       lookup_data,                               & ! intent(in):    lookup table data structure
                       ! input: enthalpy state variables  
                       scalarCanairEnthalpyTrial,                 & ! intent(in):    trial value for enthalpy of the canopy air space (J m-3)
                       scalarCanopyEnthalpyTrial,                 & ! intent(in):    trial value for enthalpy of the vegetation canopy (J m-3)
                       mLayerEnthalpyTrial,                       & ! intent(in):    trial vector of enthalpy of each snow+soil layer (J m-3)                      
                       ! output: variables for the vegetation canopy
                       scalarCanairTempTrial,                     & ! intent(inout): trial value of canopy air space temperature (K)
                       scalarCanopyTempTrial,                     & ! intent(inout): trial value of canopy temperature (K)
                       scalarCanopyWatTrial,                      & ! intent(inout): trial value of canopy total water (kg m-2)
                       scalarCanopyLiqTrial,                      & ! intent(inout): trial value of canopy liquid water (kg m-2)
                       scalarCanopyIceTrial,                      & ! intent(inout): trial value of canopy ice content (kg m-2)
                       scalarCanopyTempPrime,                     & ! intent(inout): trial value of time derivative canopy temperature (K)
                       scalarCanopyWatPrime,                      & ! intent(inout): trial value of time derivative canopy total water (kg m-2)
                       scalarCanopyLiqPrime,                      & ! intent(inout): trial value of time derivative canopy liquid water (kg m-2)
                       scalarCanopyIcePrime,                      & ! intent(inout): trial value of time derivative canopy ice content (kg m-2)
                       ! output: variables for the snow-soil domain
                       mLayerTempTrial,                           & ! intent(inout): trial vector of layer temperature (K)
                       mLayerVolFracWatTrial,                     & ! intent(inout): trial vector of volumetric total water content (-)
                       mLayerVolFracLiqTrial,                     & ! intent(inout): trial vector of volumetric liquid water content (-)
                       mLayerVolFracIceTrial,                     & ! intent(inout): trial vector of volumetric ice water content (-)
                       mLayerMatricHeadTrial,                     & ! intent(inout): trial vector of total water matric potential (m)
                       mLayerMatricHeadLiqTrial,                  & ! intent(inout): trial vector of liquid water matric potential (m)
                       mLayerTempPrime,                           & ! intent(inout): trial value of time derivative layer temperature (K)
                       mLayerVolFracWatPrime,                     & ! intent(inout): trial value of time derivative volumetric total water content (-)
                       mLayerVolFracLiqPrime,                     & ! intent(inout): trial value of time derivative volumetric liquid water content (-)
                       mLayerVolFracIcePrime,                     & ! intent(inout): trial value of time derivative volumetric ice water content (-)
                       mLayerMatricHeadPrime,                     & ! intent(inout): trial value of time derivative total water matric potential (m)
                       mLayerMatricHeadLiqPrime,                  & ! intent(inout): trial value of time derivative liquid water matric potential (m)
                       ! output: error control
                       err,message)                                 ! intent(out):   error control
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------------------------------------------------------
    use device_data_types
    implicit none
    ! input
    integer(i4b)     ,intent(in),device       :: ixNrgConserv                ! flag if enthalpy is the state variable
    logical(lgt)     ,intent(in)       :: use_lookup                      ! flag to use the lookup table for soil enthalpy, otherwise use hypergeometric function
    logical(lgt)     ,intent(in),device       :: computJac                       ! flag if computing for Jacobian update
    logical(lgt)     ,intent(in),device       :: do_adjustTemp                   ! flag to adjust temperature to account for the energy used in melt+freeze
    integer(i4b), intent(in) :: nGRU
    type(mpar_data_device),intent(in)       :: mpar_data                       ! model parameters for a local HRU
    type(indx_data_device),intent(in)       :: indx_data                       ! indices defining model states and layers
    type(prog_data_device),intent(in)       :: prog_data                       ! prognostic variables for a local HRU
    type(diag_data_device),intent(inout)    :: diag_data                       ! diagnostic variables for a local HRU
    type(deriv_data_device),intent(inout)    :: deriv_data                      ! derivatives in model fluxes w.r.t. relevant state variables
    type(zLookup)    ,intent(in)       :: lookup_data                     ! lookup tables
    ! input: enthalpy state variables  
    real(rkind),intent(in),device             :: scalarCanairEnthalpyTrial(:)       ! trial value for enthalpy of the canopy air space (J m-3)
    real(rkind),intent(in),device             :: scalarCanopyEnthalpyTrial(:)       ! trial value for enthalpy of the vegetation canopy (J m-3)
    real(rkind),intent(in),device             :: mLayerEnthalpyTrial(:,:)          ! trial vector of enthalpy of each snow+soil layer (J m-3)                      
    ! output: variables for the vegetation canopy
    real(rkind),intent(inout),device          :: scalarCanairTempTrial(:)           ! trial value of canopy air space temperature (K)
    real(rkind),intent(inout),device          :: scalarCanopyTempTrial(:)           ! trial value of canopy temperature (K)
    real(rkind),intent(inout),device          :: scalarCanopyWatTrial(:)            ! trial value of canopy total water (kg m-2)
    real(rkind),intent(inout),device          :: scalarCanopyLiqTrial(:)            ! trial value of canopy liquid water (kg m-2)
    real(rkind),intent(inout),device          :: scalarCanopyIceTrial(:)            ! trial value of canopy ice content (kg m-2)
    real(rkind),intent(inout),device          :: scalarCanopyTempPrime(:)           ! trial value of time derivative canopy temperature (K)
    real(rkind),intent(inout),device          :: scalarCanopyWatPrime(:)            ! trial value of time derivative canopy total water (kg m-2)
    real(rkind),intent(inout),device          :: scalarCanopyLiqPrime(:)            ! trial value of time derivative canopy liquid water (kg m-2)
    real(rkind),intent(inout),device          :: scalarCanopyIcePrime(:)            ! trial value of time derivative canopy ice content (kg m-2)
    ! output: variables for the snow-soil domain
    real(rkind),intent(inout),device          :: mLayerTempTrial(:,:)              ! trial vector of layer temperature (K)
    real(rkind),intent(inout),device          :: mLayerVolFracWatTrial(:,:)        ! trial vector of volumetric total water content (-)
    real(rkind),intent(inout),device          :: mLayerVolFracLiqTrial(:,:)        ! trial vector of volumetric liquid water content (-)
    real(rkind),intent(inout),device          :: mLayerVolFracIceTrial(:,:)        ! trial vector of volumetric ice water content (-)
    real(rkind),intent(inout),device          :: mLayerMatricHeadTrial(:,:)        ! trial vector of total water matric potential (m)
    real(rkind),intent(inout),device          :: mLayerMatricHeadLiqTrial(:,:)     ! trial vector of liquid water matric potential (m)
    real(rkind),intent(inout),device          :: mLayerTempPrime(:,:)              ! trial value of time derivative layer temperature (K)
    real(rkind),intent(inout),device          :: mLayerVolFracWatPrime(:,:)        ! trial value of time derivative volumetric total water content (-)
    real(rkind),intent(inout),device          :: mLayerVolFracLiqPrime(:,:)        ! trial value of time derivative volumetric liquid water content (-)
    real(rkind),intent(inout),device          :: mLayerVolFracIcePrime(:,:)        ! trial value of time derivative volumetric ice water content (-)
    real(rkind),intent(inout),device          :: mLayerMatricHeadPrime(:,:)        ! trial value of time derivative total water matric potential (m)
    real(rkind),intent(inout),device          :: mLayerMatricHeadLiqPrime(:,:)     ! trial value of time derivative liquid water matric potential (m)
    ! output: error control
    integer(i4b),intent(out)           :: err                             ! error code
    character(*),intent(out)           :: message                         ! error message
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! general local variables
    ! integer(i4b)                       :: iState                          ! index of model state variable
    ! integer(i4b)                       :: iLayer                          ! index of layer within the snow+soil domain
    ! integer(i4b) :: iGRU
    ! integer(i4b)                       :: ixFullVector                    ! index within full state vector
    ! integer(i4b)                       :: ixDomainType                    ! name of a given model domain
    ! integer(i4b)                       :: ixControlIndex                  ! index within a given model domain
    ! integer(i4b),device                       :: ixOther,ixOtherLocal            ! index of the coupled state variable within the (full, local) vector
    ! logical(lgt),device,dimension(size(indx_data%ixMapSubset2Full),nGRU)                       :: isCoupled                       ! .true. if a given variable shared another state variable in the same control volume
    ! logical(lgt),device,dimension(size(indx_data%ixMapSubset2Full),nGRU)                       :: isNrgState                      ! .true. if a given variable is an energy state
    ! logical(lgt),allocatable           :: computedCoupling(:)             ! .true. if computed the coupling for a given state variable
    real(rkind)                        :: scalarVolFracLiq                ! volumetric fraction of liquid water (-)
    real(rkind)                        :: scalarVolFracIce                ! volumetric fraction of ice (-)
    real(rkind)                        :: scalarVolFracLiqPrime           ! time derivative volumetric fraction of liquid water (-)
    real(rkind)                        :: scalarVolFracIcePrime           ! time derivative volumetric fraction of ice (-)
    real(rkind)                        :: Tcrit                           ! critical soil temperature below which ice exists (K)
    real(rkind)                        :: xTemp                           ! temporary temperature (K)
    real(rkind)                        :: fLiq                            ! fraction of liquid water (-)
    real(rkind)                        :: effSat                          ! effective saturation (-)
    real(rkind)                        :: avPore                          ! available pore space (-)
    character(len=256)                 :: cMessage                        ! error message of downwind routine
    logical(lgt),parameter             :: printFlag=.false.               ! flag to turn on printing
    ! iterative solution for temperature
    real(rkind)                        :: meltNrg                         ! energy for melt+freeze (J m-3)
    real(rkind)                        :: residual                        ! residual in the energy equation (J m-3)
    real(rkind)                        :: derivative                      ! derivative in the energy equation (J m-3 K-1)
    real(rkind)                        :: tempInc                         ! iteration increment (K)
    integer(i4b)                       :: iter                            ! iteration index
    integer(i4b)                       :: niter                           ! number of iterations
    integer(i4b),parameter             :: maxiter=100                     ! maximum number of iterations
    real(rkind),parameter              :: nrgConvTol=1.e-4_rkind          ! convergence tolerance for energy (J m-3)
    real(rkind),parameter              :: tempConvTol=1.e-6_rkind         ! convergence tolerance for temperature (K)
    real(rkind)                        :: critDiff                        ! temperature difference from critical (K)
    real(rkind)                        :: tempMin                         ! minimum bracket for temperature (K)
    real(rkind)                        :: tempMax                         ! maximum bracket for temperature (K)
    logical(lgt)                       :: bFlag                           ! flag to denote that iteration increment was constrained using bi-section
    real(rkind),parameter              :: epsT=1.e-7_rkind                ! small interval above/below critical temperature (K)
    type(dim3) :: blocks, threads
    ! iGRU = 1
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! make association with variables in the data structures
    associate(&
      ! number of model layers, and layer type
      nSnow                   => indx_data%nSnow                    ,& ! intent(in):  [i4b]    total number of snow layers
      temperature => lookup_data%temperature, psiLiq_int => lookup_data%psiLiq_int, deriv2 => lookup_data%deriv2, &
      nSoil                   => indx_data%nSoil                    ,& ! intent(in):  [i4b]    total number of soil layers
      ! nLayers                 => indx_data%nLayers                  ,& ! intent(in):  [i4b]    total number of snow and soil layers
      mLayerDepth             => prog_data%mLayerDepth                  ,& ! intent(in):  [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
      ! indices defining model states and layers
      ixVegNrg                => indx_data%ixVegNrg                 ,& ! intent(in):  [i4b]    index of canopy energy state variable
      ixVegHyd                => indx_data%ixVegHyd                 ,& ! intent(in):  [i4b]    index of canopy hydrology state variable (mass)
      ! indices in the full vector for specific domains
      ! ixNrgCanair             => indx_data%ixNrgCanair                 ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for energy states in canopy air space domain
      ixNrgCanopy_m             => indx_data%ixNrgCanopy                 ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for energy states in the canopy domain
      ixHydCanopy_m             => indx_data%ixHydCanopy                 ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the canopy domain
      ixNrgLayer              => indx_data%ixNrgLayer                  ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
      ixHydLayer              => indx_data%ixHydLayer                  ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
      ! mapping between the full state vector and the state subset  
      ! ixMapFull2Subset_m        => indx_data%ixMapFull2Subset            ,& ! intent(in):  [i4b(:)] list of indices in the state subset for each state in the full state vector
      ! ixMapSubset2Full_m        => indx_data%ixMapSubset2Full            ,& ! intent(in):  [i4b(:)] [state subset] list of indices of the full state vector in the state subset
      ! type of domain, type of state variable, and index of control volume within domain
      ixDomainType_subset_m     => indx_data%ixDomainType         ,& ! intent(in):  [i4b(:)] [state subset] id of domain for desired model state variables
      ixControlVolume_m         => indx_data%ixControlVolume             ,& ! intent(in):  [i4b(:)] index of the control volume for different domains (veg, snow, soil)
      ixStateType_m             => indx_data%ixStateType                 ,& ! intent(in):  [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
      ! snow parameters
      snowfrz_scale           => mpar_data%snowfrz_scale            ,& ! intent(in):  [dp    ] scaling parameter for the snow freezing curve (K-1)
      ! depth-varying model parameters
      soil_dens_intr          => mpar_data%soil_dens_intr              ,& ! intent(in):  [dp(:)]  intrinsic soil density (kg m-3)
      vGn_m                   => diag_data%scalarVGn_m_m                  ,& ! intent(in):  [dp(:)]  van Genutchen "m" parameter (-)
      vGn_n                   => mpar_data%vGn_n                       ,& ! intent(in):  [dp(:)]  van Genutchen "n" parameter (-)
      vGn_alpha               => mpar_data%vGn_alpha                   ,& ! intent(in):  [dp(:)]  van Genutchen "alpha" parameter (m-1)
      theta_sat               => mpar_data%theta_sat                   ,& ! intent(in):  [dp(:)]  soil porosity (-)
      theta_res               => mpar_data%theta_res                   ,& ! intent(in):  [dp(:)]  soil residual volumetric water content (-)
      ! model diagnostic variables (heat capacity, enthalpy)
      specificHeatVeg         => mpar_data%specificHeatVeg          ,& ! intent(in):  [dp   ]  specific heat of vegetation (J kg-1 K-1)
      maxMassVegetation       => mpar_data%maxMassVegetation        ,& ! intent(in):  [dp   ]  maximum mass of vegetation (kg m-2)
      canopyDepth             => diag_data%scalarCanopyDepth         ,& ! intent(in):  [dp   ]  canopy depth (m)
      scalarBulkVolHeatCapVeg => diag_data%scalarBulkVolHeatCapVeg   ,& ! intent(in):  [dp   ]  volumetric heat capacity of the vegetation (J m-3 K-1)
      mLayerVolHtCapBulk      => diag_data%mLayerVolHtCapBulk_m           ,& ! intent(in):  [dp(:)]  volumetric heat capacity in each layer (J m-3 K-1)
      ! model diagnostic variables (fraction of liquid water)
      scalarFracLiqVeg        => diag_data%scalarFracLiqVeg          ,& ! intent(out): [dp]     fraction of liquid water on vegetation (-)
      mLayerFracLiqSnow       => diag_data%mLayerFracLiqSnow_m            ,& ! intent(out): [dp(:)]  fraction of liquid water in each snow layer (-)
      ! model states from a previous solution
      scalarCanopyTemp        => prog_data%scalarCanopyTemp          ,& ! intent(in):  [dp]     temperature of the vegetation canopy (K)
      mLayerTemp              => prog_data%mLayerTemp                   ,& ! intent(in):  [dp(:)]  temperature of each snow/soil layer (K)
      ! model diagnostic variables from a previous solution
      scalarCanopyIce         => prog_data%scalarCanopyIce           ,& ! intent(in):  [dp(:)]  mass of ice on the vegetation canopy (kg m-2)
      mLayerVolFracIce        => prog_data%mLayerVolFracIce             ,& ! intent(in):  [dp(:)]  volumetric fraction of ice (-)
      ! derivatives
      dVolTot_dPsi0           => deriv_data%dVolTot_dPsi0_m              ,& ! intent(out): [dp(:)]  derivative in total water content w.r.t. total water matric potential
      dPsiLiq_dPsi0           => deriv_data%dPsiLiq_dPsi0_m              ,& ! intent(out): [dp(:)]  derivative in liquid water matric pot w.r.t. the total water matric pot (-)
      dPsiLiq_dTemp           => deriv_data%dPsiLiq_dTemp_m              ,& ! intent(out): [dp(:)]  derivative in the liquid water matric potential w.r.t. temperature
      mLayerdTheta_dTk        => deriv_data%mLayerdTheta_dTk_m           ,& ! intent(out): [dp(:)]  derivative of volumetric liquid water content w.r.t. temperature
      dTheta_dTkCanopy        => deriv_data%dTheta_dTkCanopy        ,& ! intent(out): [dp]     derivative of volumetric liquid water content w.r.t. temperature
      dFracLiqSnow_dTk        => deriv_data%dFracLiqSnow_dTk_m           ,& ! intent(out): [dp(:)]  derivative in fraction of liquid snow w.r.t. temperature
      dFracLiqVeg_dTkCanopy   => deriv_data%dFracLiqVeg_dTkCanopy   ,& ! intent(out): [dp   ]  derivative in fraction of (throughfall + drainage) w.r.t. temperature
      ! derivatives inside solver for Jacobian only
      d2VolTot_dPsi02         => deriv_data%d2VolTot_dPsi02_m            ,& ! intent(out): [dp(:)]  second derivative in total water content w.r.t. total water matric potential
      mLayerd2Theta_dTk2      => deriv_data%mLayerd2Theta_dTk2_m         ,& ! intent(out): [dp(:)]  second derivative of volumetric liquid water content w.r.t. temperature
      d2Theta_dTkCanopy2      => deriv_data%d2Theta_dTkCanopy2      ,& ! intent(out): [dp   ]  second derivative of volumetric liquid water content w.r.t. temperature
      ! derivatives of temperature if enthalpy is the state variable
      dCanairTemp_dEnthalpy     => deriv_data%dCanairTemp_dEnthalpy ,& ! intent(out): [dp]     derivative of canopy air temperature w.r.t. enthalpy
      dCanopyTemp_dEnthalpy     => deriv_data%dCanopyTemp_dEnthalpy ,& ! intent(out): [dp]     derivative of canopy temperature w.r.t. enthalpy 
      dTemp_dEnthalpy           => deriv_data%dTemp_dEnthalpy_m          ,& ! intent(out): [dp(:)]  derivative of temperature w.r.t. enthalpy
      dCanopyTemp_dCanWat       => deriv_data%dCanopyTemp_dCanWat   ,& ! intent(out): [dp]     derivative of canopy temperature w.r.t. volumetric water content
      dTemp_dTheta              => deriv_data%dTemp_dTheta_m             ,& ! intent(out): [dp(:)]  derivative of temperature w.r.t. volumetric water content
      dTemp_dPsi0               => deriv_data%dTemp_dPsi0_m               & ! intent(out): [dp(:)]  derivative of temperature w.r.t. total water matric potential
      ) ! association with variables in the data structures
  
      ! --------------------------------------------------------------------------------------------------------------------------------
      ! --------------------------------------------------------------------------------------------------------------------------------
  
      ! initialize error control
      err=0; message='updateVarsWithPrime/'
      ! print*, shape(temperature)
  
      ! allocate space and assign values to the flag vector
      ! allocate(computedCoupling(size(ixMapSubset2Full)),stat=err)        ! .true. if computed the coupling for a given state variable
      if(err/=0)then; message=trim(message)//'problem allocating computedCoupling'; return; endif
      ! computedCoupling(:)=.false.
  
      ! loop through model state variables
      ! !$cuf kernel do(1) <<<1,1>>>
      ! do iGRU=1,nGRU
      ! do iState=1,size(ixMapSubset2Full)
  
        ! check the need for the computations
        ! if(computedCoupling(iState)) cycle
  
        ! -----
        ! - compute indices...
        ! --------------------
  
        ! get the layer index
        ! select case(ixDomainType(iState,iGRU))
        !   case(iname_cas);     iLayer = 0
        !   case(iname_veg);     iLayer = 0
        !   case(iname_snow);    iLayer = ixControlIndex
        !   case(iname_soil);    iLayer = ixControlIndex + nSnow
        !   case(iname_aquifer); cycle ! aquifer: do nothing
        !   ! case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
        ! end select
  
        ! if(ixOtherLocal/=integerMissing) computedCoupling(ixOtherLocal)=.true.
  
  
  
        ! if(printFlag)then
        !   print*, 'iState         = ', iState, size(ixMapSubset2Full)
        !   print*, 'ixFullVector   = ', ixFullVector
        !   print*, 'ixDomainType   = ', ixDomainType(iState)
        !   print*, 'ixControlIndex = ', ixControlIndex
        !   print*, 'ixOther        = ', ixOther
        !   print*, 'ixOtherLocal   = ', ixOtherLocal
        !   print*, 'do_adjustTemp  = ', do_adjustTemp
        !   print*, 'isCoupled      = ', isCoupled
        !   print*, 'isNrgState     = ', isNrgState
        ! endif
        ! if (iState < ixOtherLocal .or. ixOtherLocal==integerMissing) then
          ! iGRU_d = iGRU
          ! ixControlIndex_d = ixControlIndex
          ! print*, iState, iGRU, ixControlIndex
          ! print*, temperature(1,ixControlIndex,iGRU)
          threads = dim3(8,16,1)
          blocks = dim3(size(ixDomainType_subset_m,1)/8+1,nGRU/16+1,1)
        call update<<<blocks,threads>>>(mLayerVolFracWatPrime, mLayerVolFracWatTrial, &
        mLayerVolFracLiqPrime, mLayerVolFracLiqTrial, &
        mLayerVolFracIcePrime, mLayerVolFracIceTrial, &
        mLayerMatricHeadPrime, mLayerMatricHeadTrial, &
        mLayerMatricHeadLiqPrime, mLayerMatricHeadLiqTrial, &
        computJac, do_adjustTemp, ixNrgConserv, &
        scalarCanairEnthalpyTrial, &
        scalarCanairTempTrial, &
        dCanairTemp_dEnthalpy, &
        vGn_alpha, vGn_m, vGn_n, &
        theta_res, theta_sat, &
        canopyDepth, specificHeatVeg, &
        maxMassVegetation, snowfrz_scale, &
        scalarCanopyEnthalpyTrial, scalarCanopyWatTrial, scalarCanopyTempTrial, &
        dCanopyTemp_dEnthalpy, dCanopyTemp_dCanWat, &
        dTemp_dEnthalpy, dTemp_dTheta, &
        mLayerTempPrime, mLayerTempTrial, mLayerEnthalpyTrial, &
        soil_dens_intr, dTemp_dPsi0, &
        dVolTot_dPsi0, d2VolTot_dPsi02, &
        dFracLiqSnow_dTk, mLayerdTheta_dTk, mLayerd2Theta_dTk2, &
        nGRU, &
        ixStateType_m, &
        scalarCanopyLiqTrial, scalarCanopyIceTrial, &
        scalarCanopyWatPrime, scalarCanopyLiqPrime, scalarCanopyIcePrime, &
        dPsiLiq_dPsi0, dPsiLiq_dTemp, &
        dFracLiqVeg_dTkCanopy, dTheta_dTkCanopy, d2Theta_dTkCanopy2, &
        scalarFracLiqVeg, &
        mLayerFracLiqSnow, &
        scalarCanopyTempPrime, scalarCanopyTemp, scalarCanopyIce, &
        scalarBulkVolHeatCapVeg, &
        mLayerVolHtCapBulk, mLayerTemp, mLayerVolFracIce, &
        ! ixMapFull2Subset_m, &
        ixHydCanopy_m, ixNrgCanopy_m, ixHydLayer, ixNrgLayer, nSnow,ixControlVolume_m, &
        ! ixMapSubset2Full_m, &
        ixDomainType_subset_m,temperature,&
        psiLiq_int, deriv2)
        ! endif
      ! end do ! looping through state variables
    ! enddo
  
      ! deallocate space
      ! deallocate(computedCoupling,stat=err)        ! .true. if computed the coupling for a given state variable
      if(err/=0)then; message=trim(message)//'problem deallocating computedCoupling'; return; endif
  
      ! end association to the variables in the data structures
  end associate
  
  end subroutine updateVarsWithPrime
  
  
  attributes(global) subroutine update(mLayerVolFracWatPrime_in, mLayerVolFracWatTrial_, &
                    mLayerVolFracLiqPrime_, mLayerVolFracLiqTrial_, &
                    mLayerVolFracIcePrime_, mLayerVolFracIceTrial_, &
                    mLayerMatricHeadPrime_, mLayerMatricHeadTrial_, &
                    mLayerMatricHeadLiqPrime_, mLayerMatricHeadLiqTrial_, &
                    computJac, do_adjustTemp, ixNrgConserv, &
                    ! isCoupled, &
                    scalarCanairEnthalpyTrial_, &
                    scalarCanairTempTrial_, &
                    dCanairTemp_dEnthalpy_, &
                    vGn_alpha_, vGn_m_, vGn_n_, &
                    theta_res_, theta_sat_, &
                    canopyDepth_, specificHeatVeg, &
                    maxMassVegetation, snowfrz_scale, &
                    scalarCanopyEnthalpyTrial, scalarCanopyWatTrial, scalarCanopyTempTrial_, &
                    dCanopyTemp_dEnthalpy_, dCanopyTemp_dCanWat_, &
                    dTemp_dEnthalpy_, dTemp_dTheta_, &
                    mLayerTempPrime_, mLayerTempTrial_, mLayerEnthalpyTrial_, &
                    soil_dens_intr_, dTemp_dPsi0_, &
                    dVolTot_dPsi0_, d2VolTot_dPsi02_, &
                    dFracLiqSnow_dTk_, mLayerdTheta_dTk_, mLayerd2Theta_dTk2_, &
                    nGRU,&
                    ixStateType_, &
                    scalarCanopyLiqTrial_, scalarCanopyIceTrial_, &
                    scalarCanopyWatPrime_, scalarCanopyLiqPrime_, scalarCanopyIcePrime_, &
                    dPsiLiq_dPsi0_, dPsiLiq_dTemp_, &
                    dFracLiqVeg_dTkCanopy_, dTheta_dTkCanopy_, d2Theta_dTkCanopy2_, &
                    scalarFracLiqVeg_, &
                    mLayerFracLiqSnow_, &
                    scalarCanopyTempPrime_, scalarCanopyTemp_, scalarCanopyIce_, &
                    scalarBulkVolHeatCapVeg_, &
                    mLayerVolHtCapBulk_, mLayerTemp_, mLayerVolFracIce_, &
                    ! ixMapFull2Subset, &
                    ixHydCanopy, ixNrgCanopy, ixHydLayer, ixNrgLayer, nSnow,ixControlVolume, &
                    ! ixMapSubset2Full, &
                    ixDomainType_subset, Tk, Ly, L2)
    real(rkind),intent(inout) :: mLayerVolFracWatPrime_in(:,:), mLayerVolFracWatTrial_(:,:) ! *(iLayer)
    real(rkind),intent(inout) :: mLayerVolFracLiqPrime_(:,:), mLayerVolFracLiqTrial_(:,:) ! *(iLayer)
    real(rkind),intent(inout) :: mLayerVolFracIcePrime_(:,:), mLayerVolFracIceTrial_(:,:) ! *(iLayer)
    real(rkind),intent(inout) :: mLayerMatricHeadPrime_(:,:), mLayerMatricHeadTrial_(:,:) ! *(ixControlIndex)
    real(rkind),intent(inout) :: mLayerMatricHeadLiqPrime_(:,:), mLayerMatricHeadLiqTrial_(:,:) ! *(ixControlIndex)
    logical(lgt),intent(in) :: computJac, do_adjustTemp
    integer(i4b),intent(in) :: ixNrgConserv
    logical(lgt) :: enthalpyStateVec
    ! logical(lgt),intent(in) :: isCoupled
    real(rkind),intent(in) :: scalarCanairEnthalpyTrial_(:)
    real(rkind),intent(inout) :: scalarCanairTempTrial_(:)
    real(rkind),intent(inout) :: dCanairTemp_dEnthalpy_(:)
    real(rkind),intent(in) :: vGn_alpha_(:), vGn_m_(:,:), vGn_n_(:) ! *(ixControlIndex)
    real(rkind),intent(in) :: theta_res_(:), theta_sat_(:) ! *(ixControlIndex)
    real(rkind),intent(in) :: canopyDepth_(:), specificHeatVeg
    real(rkind),intent(in) :: maxMassVegetation, snowfrz_scale
    real(rkind),intent(in) :: scalarCanopyEnthalpyTrial(:)
    real(rkind),intent(inout) ::  scalarCanopyWatTrial(:), scalarCanopyTempTrial_(:)
    real(rkind),intent(inout) :: dCanopyTemp_dEnthalpy_(:), dCanopyTemp_dCanWat_(:)
    real(rkind),intent(inout) :: dTemp_dEnthalpy_(:,:), dTemp_dTheta_(:,:) ! *(iLayer)
    real(rkind),intent(inout) :: mLayerTempPrime_(:,:), mLayerTempTrial_(:,:), mLayerEnthalpyTrial_(:,:) ! *(iLayer)
    real(rkind),intent(inout) :: soil_dens_intr_(:), dTemp_dPsi0_(:,:) ! *(ixControlIndex)
    real(rkind),intent(inout) :: dVolTot_dPsi0_(:,:), d2VolTot_dPsi02_(:,:) ! *(ixControlIndex)
    real(rkind),intent(inout) :: dFracLiqSnow_dTk_(:,:), mLayerdTheta_dTk_(:,:), mLayerd2Theta_dTk2_(:,:) ! *(iLayer)
    integer(i4b),value :: nGRU
    integer(i4b),intent(in) :: ixStateType_(:,:) ! (ixFullVector)
    real(rkind),intent(inout) :: scalarCanopyLiqTrial_(:), scalarCanopyIceTrial_(:)
    real(rkind),intent(inout) :: scalarCanopyWatPrime_(:), scalarCanopyLiqPrime_(:), scalarCanopyIcePrime_(:)
    real(rkind),intent(inout) :: dPsiLiq_dPsi0_(:,:), dPsiLiq_dTemp_(:,:) !*(ixControlIndex)
    real(rkind),intent(inout) :: dFracLiqVeg_dTkCanopy_(:), dTheta_dTkCanopy_(:), d2Theta_dTkCanopy2_(:)
    real(rkind),intent(inout) :: scalarFracLiqVeg_(:)
    real(rkind),intent(inout) :: mLayerFracLiqSnow_(:,:) ! (iLayer)
    real(rkind),intent(inout) :: scalarCanopyTempPrime_(:), scalarCanopyTemp_(:), scalarCanopyIce_(:)
    real(rkind),intent(inout) :: scalarBulkVolHeatCapVeg_(:)
    real(rkind),intent(inout) :: mLayerVolHtCapBulk_(:,:), mLayerTemp_(:,:), mLayerVolFracIce_(:,:) ! (iLayer)
    real(rkind) :: Tk(:,:,:)
    real(rkind) ::  Ly(:,:,:), L2(:,:,:)
    integer(i4b) :: nrgConserv
    ! integer(i4b) :: ixMapFull2Subset(:,:)
    integer(i4b) :: ixHydCanopy(:), ixNrgCanopy(:), ixHydLayer(:,:), ixNrgLayer(:,:)
    integer(i4b) :: nSnow(:)
    integer(i4b) :: ixControlVolume(:,:)
    ! integer(i4b) :: ixMapSubset2Full(:,:), 
    integer(i4b) :: ixDomainType_subset(:,:)
    ! integer(i4b) :: iLayer
    ! integer(i4b) :: ixOther
  
    logical(lgt) :: use_lookup
  
    integer(i4b) :: err
    ! integer(i4b)                       :: iState                          ! index of model state variable
    integer(i4b)                       :: iLayer                          ! index of layer within the snow+soil domain
    integer(i4b)                       :: ixFullVector                    ! index within full state vector
    integer(i4b)                       :: ixDomainType                    ! name of a given model domain
    integer(i4b) :: ixStateType
    integer(i4b)                       :: ixControlIndex                  ! index within a given model domain
    integer(i4b)                       :: ixOther, ixOtherLocal            ! index of the coupled state variable within the (full, local) vector
    logical(lgt)                       :: isCoupled                       ! .true. if a given variable shared another state variable in the same control volume
    logical(lgt)                       :: isNrgState                      ! .true. if a given variable is an energy state
    ! logical(lgt),allocatable           :: computedCoupling(:)             ! .true. if computed the coupling for a given state variable
    real(rkind)                        :: scalarVolFracLiq                ! volumetric fraction of liquid water (-)
    real(rkind)                        :: scalarVolFracIce                ! volumetric fraction of ice (-)
    real(rkind)                        :: scalarVolFracLiqPrime           ! time derivative volumetric fraction of liquid water (-)
    real(rkind)                        :: scalarVolFracIcePrime           ! time derivative volumetric fraction of ice (-)
    real(rkind)                        :: Tcrit                           ! critical soil temperature below which ice exists (K)
    real(rkind)                        :: xTemp                           ! temporary temperature (K)
    real(rkind)                        :: fLiq                            ! fraction of liquid water (-)
    real(rkind)                        :: effSat                          ! effective saturation (-)
    real(rkind)                        :: avPore                          ! available pore space (-)
    logical(lgt),parameter             :: printFlag=.false.               ! flag to turn on printing
    ! iterative solution for temperature
    real(rkind)                        :: meltNrg                         ! energy for melt+freeze (J m-3)
    real(rkind)                        :: residual                        ! residual in the energy equation (J m-3)
    real(rkind)                        :: derivative                      ! derivative in the energy equation (J m-3 K-1)
    real(rkind)                        :: tempInc                         ! iteration increment (K)
    integer(i4b)                       :: iter                            ! iteration index
    integer(i4b)                       :: niter                           ! number of iterations
    integer(i4b),parameter             :: maxiter=100                     ! maximum number of iterations
    real(rkind),parameter              :: nrgConvTol=1.e-4_rkind          ! convergence tolerance for energy (J m-3)
    real(rkind),parameter              :: tempConvTol=1.e-6_rkind         ! convergence tolerance for temperature (K)
    real(rkind)                        :: critDiff                        ! temperature difference from critical (K)
    real(rkind)                        :: tempMin                         ! minimum bracket for temperature (K)
    real(rkind)                        :: tempMax                         ! maximum bracket for temperature (K)
    logical(lgt)                       :: bFlag                           ! flag to denote that iteration increment was constrained using bi-section
    real(rkind),parameter              :: epsT=1.e-7_rkind                ! small interval above/below critical temperature (K)
    integer(i4b) :: iState, iGRU
    ! real(rkind) :: mLayerVolFracWatPrime
    ! print*, 538
  
    iState = (blockidx%x-1) * blockdim%x + threadidx%x
    iGRU = (blockidx%y-1) * blockdim%y + threadidx%y
    ! print*, iState
    if (iState .gt. size(ixDomainType_subset,1) .or. iGRU .gt. nGRU) return
    enthalpyStateVec = ixNrgConserv.ne.closedForm
        ! get domain type, and index of the control volume within the domain
    ixFullVector   = iState       ! index within full state vector
    ixDomainType   = ixDomainType_subset(iState,iGRU)    ! named variables defining the domain (iname_cas, iname_veg, etc.)
    ixControlIndex = ixControlVolume(ixFullVector,iGRU)  ! index within a given domain
    ixStateType = ixStateType_(ixFullVector,iGRU)
    ! print*, ixDomainType, iGRU, iState
  
    ! print*, ixControlIndex, ixFullVector, ixDomainType, iname_cas
    ! get the layer index
    select case(ixDomainType)
    case(iname_cas);     iLayer = 0
    case(iname_veg);     iLayer = 0
    case(iname_snow);    iLayer = ixControlIndex
    case(iname_soil);    iLayer = ixControlIndex + nSnow(iGRU)
    case(iname_aquifer); return ! aquifer: do nothing
    case default; return
    ! case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil, iname_aquifer'; return
  end select

  ! print*, iGRU, iState, ixControlIndex
  
          ! get the index of the other (energy or mass) state variable within the full state vector
    select case(ixDomainType)
      case(iname_cas)             ; ixOther = integerMissing
      case(iname_veg)             ; ixOther = merge(ixHydCanopy(iGRU),    ixNrgCanopy(iGRU),    ixStateType==iname_nrgCanopy)
      case(iname_snow, iname_soil); ixOther = merge(ixHydLayer(iLayer,iGRU),ixNrgLayer(iLayer,iGRU),ixStateType==iname_nrgLayer)
      ! case default; err=20; message=trim(message)//'expect case to be iname_cas, iname_veg, iname_snow, iname_soil'; return
    end select
  
    ! print*, iLayer, ixOther
    ! get the index in the local state vector
    if(ixDomainType==iname_cas)then
      ixOtherLocal = integerMissing
    else
      ixOtherLocal = ixOther  ! ixOtherLocal could equal integerMissing
    endif
    ! print*, ixOtherLocal, iLayer, iGRU
  
    ! check if we have a coupled solution
    isCoupled    = (ixOtherLocal/=integerMissing)
  
    ! check if we are an energy state
    isNrgState   = (ixStateType==iname_nrgCanair .or. ixStateType==iname_nrgCanopy .or. ixStateType==iname_nrgLayer)
  
  
    ! print*, iState, ixOtherLocal, ixOther, iLayer
    if (.not. (iState < ixOtherLocal .or. ixOtherLocal==integerMissing)) return
    ! print*, iState, iGRU
    ! mLayerVolFracWatPrime = mLayerVolFracWatPrime_in(iLayer,iGRU)
  
    use_lookup = ixnrgConserv .eq. enthalpyFormLU
    ! print*, Tk(1,ixControlIndex,iGRU)
    ! compute temperature from enthalpy for canopy air space
    if(ixDomainType==iname_cas)then
      if(enthalpyStateVec)then
        call enthalpy2T_cas(&
                 computJac,                  & ! intent(in):  flag if computing for Jacobian update
                 scalarCanairEnthalpyTrial_(iGRU),  & ! intent(in):  trial value for enthalpy of the canopy air space (J m-3)
                 scalarCanairTempTrial_(iGRU),      & ! intent(out): trial value for canopy air temperature (K)
                 dCanairTemp_dEnthalpy_(iGRU),      & ! intent(out): derivative of canopy air temperature with enthalpy
                 err)                 ! intent(out): error control
        if(err/=0)then;  return; endif
      else
          dCanairTemp_dEnthalpy_(iGRU) = 0._rkind
      endif
      return ! no more to do on canopy air space
    end if
  
    ! update hydrology state variables for the uncoupled solution
    if(.not.isNrgState .and. .not.isCoupled)then
  
      ! update the total water from volumetric liquid water
      if(ixStateType==iname_liqCanopy .or. ixStateType==iname_liqLayer)then
        select case(ixDomainType)
          case(iname_veg)
              scalarCanopyWatTrial(iGRU)          = scalarCanopyLiqTrial_(iGRU)          + scalarCanopyIceTrial_(iGRU)
              scalarCanopyWatPrime_(iGRU)          = scalarCanopyLiqPrime_(iGRU)          + scalarCanopyIcePrime_(iGRU)
          case(iname_snow)
              mLayerVolFracWatTrial_(iLayer,iGRU) = mLayerVolFracLiqTrial_(iLayer,iGRU) + mLayerVolFracIceTrial_(iLayer,iGRU)*iden_ice/iden_water
              mLayerVolFracWatPrime_in(iLayer,iGRU) = mLayerVolFracLiqPrime_(iLayer,iGRU) + mLayerVolFracIcePrime_(iLayer,iGRU)*iden_ice/iden_water
          case(iname_soil)
              mLayerVolFracWatTrial_(iLayer,iGRU) = mLayerVolFracLiqTrial_(iLayer,iGRU) + mLayerVolFracIceTrial_(iLayer,iGRU) ! no volume expansion
              mLayerVolFracWatPrime_in(iLayer,iGRU) = mLayerVolFracLiqPrime_(iLayer,iGRU) + mLayerVolFracIcePrime_(iLayer,iGRU)
          case default; err=20; return
        end select
      endif
  
      ! update the total water and the total water matric potential
      if(ixDomainType==iname_soil)then
        select case( ixStateType )
          ! --> update the total water from the liquid water matric potential
          case(iname_lmpLayer)
            effSat = volFracLiq(mLayerMatricHeadLiqTrial_(ixControlIndex,iGRU),vGn_alpha_(ixControlIndex),0._rkind,1._rkind,vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU))  ! effective saturation
            avPore = theta_sat_(ixControlIndex) - mLayerVolFracIceTrial_(iLayer,iGRU) - theta_res_(ixControlIndex)  ! available pore space
            mLayerVolFracLiqTrial_(iLayer,iGRU) = effSat*avPore + theta_res_(ixControlIndex)
            mLayerVolFracWatTrial_(iLayer,iGRU) = mLayerVolFracLiqTrial_(iLayer,iGRU) + mLayerVolFracIceTrial_(iLayer,iGRU) ! no volume expansion
            mLayerVolFracWatPrime_in(iLayer,iGRU) = mLayerVolFracLiqPrime_(iLayer,iGRU) + mLayerVolFracIcePrime_(iLayer,iGRU)
            mLayerMatricHeadTrial_(ixControlIndex,iGRU) = matricHead(mLayerVolFracWatTrial_(iLayer,iGRU),vGn_alpha_(ixControlIndex),theta_res_(ixControlIndex),theta_sat_(ixControlIndex),vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU))
            mLayerMatricHeadPrime_(ixControlIndex,iGRU) =  dPsi_dTheta(mLayerVolFracWatTrial_(iLayer,iGRU),vGn_alpha_(ixControlIndex),theta_res_(ixControlIndex),theta_sat_(ixControlIndex),vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU)) * mLayerVolFracWatPrime_in(iLayer,iGRU)
            ! --> update the total water from the total water matric potential
          case(iname_matLayer)
            mLayerVolFracWatTrial_(iLayer,iGRU) = volFracLiq(mLayerMatricHeadTrial_(ixControlIndex,iGRU),vGn_alpha_(ixControlIndex),theta_res_(ixControlIndex),theta_sat_(ixControlIndex),vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU))
            mLayerVolFracWatPrime_in(iLayer,iGRU) = dTheta_dPsi(mLayerMatricHeadTrial_(ixControlIndex,iGRU),vGn_alpha_(ixControlIndex),theta_res_(ixControlIndex),theta_sat_(ixControlIndex),vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU)) *mLayerMatricHeadPrime_(ixControlIndex,iGRU)
            ! --> update the total water matric potential (assume already have mLayerVolFracWatTrial given block above)
          case(iname_liqLayer, iname_watLayer)
            mLayerMatricHeadTrial_(ixControlIndex,iGRU) = matricHead(mLayerVolFracWatTrial_(iLayer,iGRU),vGn_alpha_(ixControlIndex),theta_res_(ixControlIndex),theta_sat_(ixControlIndex),vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU))
            mLayerMatricHeadPrime_(ixControlIndex,iGRU) = dPsi_dTheta(mLayerVolFracWatTrial_(iLayer,iGRU),vGn_alpha_(ixControlIndex),theta_res_(ixControlIndex),theta_sat_(ixControlIndex),vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU)) * mLayerVolFracWatPrime_in(iLayer,iGRU)
          case default; err=20; return
        end select
      endif  ! if in the soil domain
  
    endif  ! if hydrology state variable or uncoupled solution
  
    ! compute temperature from enthalpy and water content for remaining domains 
    if(ixDomainType==iname_veg)then
      if(enthalpyStateVec)then
        call enthalpy2T_veg(&
                 computJac,                  & ! intent(in):    flag if computing for Jacobian update          
                 canopyDepth_(iGRU),                & ! intent(in):    canopy depth (m)
                 specificHeatVeg,            & ! intent(in):    specific heat of vegetation (J kg-1 K-1)
                 maxMassVegetation,          & ! intent(in):    maximum mass of vegetation (kg m-2)
                 snowfrz_scale,              & ! intent(in):    scaling parameter for the snow freezing curve  (K-1)
                 scalarCanopyEnthalpyTrial(iGRU),  & ! intent(in):    trial value for enthalpy of the vegetation canopy (J m-3)
                 scalarCanopyWatTrial(iGRU),       & ! intent(in):    trial value for canopy total water (kg m-2)
                 scalarCanopyTempTrial_(iGRU),      & ! intent(inout): trial value for canopy temperature (K)
                 dCanopyTemp_dEnthalpy_(iGRU),      & ! intent(inout): derivative of canopy temperature with enthalpy
                 dCanopyTemp_dCanWat_(iGRU),        & ! intent(inout): derivative of canopy temperature with canopy water
                 err)                 ! intent(out):   error control
        if(err/=0)then; return; endif
      else
        dCanopyTemp_dEnthalpy_(iGRU) = 0._rkind
        dCanopyTemp_dCanWat_(iGRU)   = 0._rkind
      endif
    elseif(ixDomainType==iname_snow)then
      if(enthalpyStateVec)then
        call enthalpy2T_snow(&
                 computJac,                      & ! intent(in):    flag if computing for Jacobian update       
                 snowfrz_scale,                  & ! intent(in):    scaling parameter for the snow freezing curve (K-1)
                 mLayerEnthalpyTrial_(iLayer,iGRU),    & ! intent(in):    enthalpy of snow+soil layer (J m-3)
                 mLayerVolFracWatTrial_(iLayer,iGRU),  & ! intent(in):    volumetric total water content (-)
                 mLayerTempTrial_(iLayer,iGRU),        & ! intent(inout): layer temperature (K)
                 dTemp_dEnthalpy_(iLayer,iGRU),        & ! intent(inout): derivative of layer temperature with enthalpy
                 dTemp_dTheta_(iLayer,iGRU),           & ! intent(inout): derivative of layer temperature with volumetric total water content
                 err)                     ! intent(out):   error control
        if(err/=0)then; return; endif
      else
        dTemp_dEnthalpy_(iLayer,iGRU) = 0._rkind
        dTemp_dTheta_(iLayer,iGRU)    = 0._rkind
      endif
    elseif(ixDomainType==iname_soil)then
      if(enthalpyStateVec)then
        call enthalpy2T_soil(&
                 computJac,                              & ! intent(in):    flag if computing for Jacobian update
                 use_lookup,                             & ! intent(in):    flag to use the lookup table for soil enthalpy
                 soil_dens_intr_(ixControlIndex),         & ! intent(in):    intrinsic soil density (kg m-3)
                 vGn_alpha_(ixControlIndex),              & ! intent(in):    van Genutchen "alpha" parameter
                 vGn_n_(ixControlIndex),                  & ! intent(in):    van Genutchen "n" parameter
                 theta_sat_(ixControlIndex),              & ! intent(in):    soil porosity (-)
                 theta_res_(ixControlIndex),              & ! intent(in):    soil residual volumetric water content (-)
                 vGn_m_(ixControlIndex,iGRU),                  & ! intent(in):    van Genutchen "m" parameter (-)
                !  ixControlIndex,                         & ! intent(in):    index of the control volume within the domain
                !  lookup_data,                            & ! intent(in):    lookup table data structure
                 Tk(:,ixControlIndex,iGRU), Ly(:,ixControlIndex,iGRU), L2(:,ixControlIndex,iGRU), &
                 mLayerEnthalpyTrial_(iLayer,iGRU),            & ! intent(in):    trial vector of enthalpy of each snow+soil layer (J m-3)
                 mLayerMatricHeadTrial_(ixControlIndex,iGRU),  & ! intent(in):    trial vector of total water matric potential (m)
                 mLayerTempTrial_(iLayer,iGRU),                & ! intent(inout): trial vector of layer temperature (K)
                 dTemp_dEnthalpy_(iLayer,iGRU),                & ! intent(inout): derivative of layer temperature with enthalpy
                 dTemp_dTheta_(iLayer,iGRU),                   & ! intent(inout): derivative of layer temperature with volumetric total water content
                 dTemp_dPsi0_(ixControlIndex,iGRU),            & ! intent(inout): derivative of layer temperature with total water matric potential
                 err)                             ! intent(out):   error control
         if(err/=0)then; return; endif
      else
        dTemp_dEnthalpy_(iLayer,iGRU)     = 0._rkind
        dTemp_dTheta_(iLayer,iGRU)        = 0._rkind
        dTemp_dPsi0_(ixControlIndex,iGRU) = 0._rkind 
      endif
    endif 
  
    ! compute the critical soil temperature below which ice exists
    select case(ixDomainType)
      case(iname_veg, iname_snow); Tcrit = Tfreeze
      case(iname_soil);            Tcrit = crit_soilT( mLayerMatricHeadTrial_(ixControlIndex,iGRU) )
      case default; err=20; return
    end select
  
    ! initialize temperature
    select case(ixDomainType)
      case(iname_veg);              xTemp = scalarCanopyTempTrial_(iGRU)
      case(iname_snow, iname_soil); xTemp = mLayerTempTrial_(iLayer,iGRU)
      case default; err=20; return
    end select
  
    ! define brackets for the root
    ! NOTE: start with an enormous range; updated quickly in the iterations
    tempMin = xTemp - 10._rkind
    tempMax = xTemp + 10._rkind
  
    ! get iterations (set to maximum iterations if adjusting the temperature)
    niter = merge(maxiter, 1, do_adjustTemp)
  
    ! iterate
    iterations: do iter=1,niter
  
      ! restrict temperature
      if(xTemp <= tempMin .or. xTemp >= tempMax)then
        xTemp = 0.5_rkind*(tempMin + tempMax)  ! new value
        bFlag = .true.
      else
        bFlag = .false.
      endif
  
      ! -----
      ! - compute derivatives...
      ! ------------------------
  
      ! compute the derivative in total water content w.r.t. total water matric potential (m-1)
      ! NOTE 1: valid for frozen and unfrozen conditions
      ! NOTE 2: for case "iname_lmpLayer", dVolTot_dPsi0 = dVolLiq_dPsi
      if(ixDomainType==iname_soil)then
        select case( ixStateType )
          case(iname_lmpLayer)
            dVolTot_dPsi0_(ixControlIndex,iGRU) = dTheta_dPsi(mLayerMatricHeadLiqTrial_(ixControlIndex,iGRU),vGn_alpha_(ixControlIndex),0._rkind,1._rkind,vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU))*avPore
            if(computJac) d2VolTot_dPsi02_(ixControlIndex,iGRU) = d2Theta_dPsi2(mLayerMatricHeadLiqTrial_(ixControlIndex,iGRU),vGn_alpha_(ixControlIndex),0._rkind,1._rkind,vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU))*avPore
          case default
            dVolTot_dPsi0_(ixControlIndex,iGRU) = dTheta_dPsi(mLayerMatricHeadTrial_(ixControlIndex,iGRU),vGn_alpha_(ixControlIndex),theta_res_(ixControlIndex),theta_sat_(ixControlIndex),vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU))
            if(computJac) d2VolTot_dPsi02_(ixControlIndex,iGRU) = d2Theta_dPsi2(mLayerMatricHeadTrial_(ixControlIndex,iGRU),vGn_alpha_(ixControlIndex),theta_res_(ixControlIndex),theta_sat_(ixControlIndex),&
                                              vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU))
          end select
      endif
  
      ! compute the derivative in liquid water content w.r.t. temperature
      ! --> partially frozen: dependence of liquid water on temperature
      if(xTemp<Tcrit)then
        select case(ixDomainType)
          case(iname_veg)
            dFracLiqVeg_dTkCanopy_(iGRU) = dFracLiq_dTk(xTemp,snowfrz_scale)
            dTheta_dTkCanopy_(iGRU) = dFracLiqVeg_dTkCanopy_(iGRU) * scalarCanopyWatTrial(iGRU)/(iden_water*canopyDepth_(iGRU))
            if(computJac)then
              fLiq = fracLiquid(xTemp,snowfrz_scale)
              d2Theta_dTkCanopy2_(iGRU) = 2._rkind * snowfrz_scale**2_i4b * ( (Tfreeze - xTemp) * 2._rkind * fLiq * dFracLiqVeg_dTkCanopy_(iGRU) - fLiq**2_i4b ) * scalarCanopyWatTrial(iGRU)/(iden_water*canopyDepth_(iGRU))
            endif
          case(iname_snow)
            dFracLiqSnow_dTk_(iLayer,iGRU) = dFracLiq_dTk(xTemp,snowfrz_scale)
            mLayerdTheta_dTk_(iLayer,iGRU) = dFracLiqSnow_dTk_(iLayer,iGRU) * mLayerVolFracWatTrial_(iLayer,iGRU)
            if(computJac)then
              fLiq = fracLiquid(xTemp,snowfrz_scale)
              mLayerd2Theta_dTk2_(iLayer,iGRU) = 2._rkind * snowfrz_scale**2_i4b * ( (Tfreeze - xTemp) * 2._rkind * fLiq * dFracLiqSnow_dTk_(iLayer,iGRU) - fLiq**2_i4b ) * mLayerVolFracWatTrial_(iLayer,iGRU)
            endif
          case(iname_soil)
            dFracLiqSnow_dTk_(iLayer,iGRU) = 0._rkind !dTheta_dTk(xTemp,theta_res(ixControlIndex),theta_sat(ixControlIndex),vGn_alpha(ixControlIndex),vGn_n(ixControlIndex),vGn_m(ixControlIndex))/ mLayerVolFracWatTrial(iLayer)
            mLayerdTheta_dTk_(iLayer,iGRU) = dTheta_dTk(xTemp,theta_res_(ixControlIndex),theta_sat_(ixControlIndex),vGn_alpha_(ixControlIndex),vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU))
            if(computJac)then
              mLayerd2Theta_dTk2_(iLayer,iGRU) = d2Theta_dTk2(xTemp,theta_res_(ixControlIndex),theta_sat_(ixControlIndex),vGn_alpha_(ixControlIndex),vGn_n_(ixControlIndex),vGn_m_(ixControlIndex,iGRU))
            endif
          case default; err=20; return
        end select  ! domain type
  
      ! --> unfrozen: no dependence of liquid water on temperature
      else
        select case(ixDomainType)
          case(iname_veg);              dTheta_dTkCanopy_(iGRU)         = 0._rkind; d2Theta_dTkCanopy2_(iGRU)         = 0._rkind; dFracLiqVeg_dTkCanopy_(iGRU)    = 0._rkind
          case(iname_snow, iname_soil); mLayerdTheta_dTk_(iLayer,iGRU) = 0._rkind; mLayerd2Theta_dTk2_(iLayer,iGRU) = 0._rkind; dFracLiqSnow_dTk_(iLayer,iGRU) = 0._rkind
        end select  ! domain type
      endif
  
      ! -----
      ! - update volumetric fraction of liquid water and ice...
      !    => case of hydrology state uncoupled with energy (and when not adjusting the temperature)...
      ! -----------------------------------------------------------------------------------------------
  
      ! case of hydrology state uncoupled with energy (and when not adjusting the temperature)
      if(.not.do_adjustTemp .and. .not.isNrgState .and. .not.isCoupled)then
  
        ! compute the fraction of snow
        select case(ixDomainType)
          case(iname_veg);   scalarFracLiqVeg_(iGRU)          = fracliquid(xTemp,snowfrz_scale)
          case(iname_snow);  mLayerFracLiqSnow_(iLayer,iGRU) = fracliquid(xTemp,snowfrz_scale)
          case(iname_soil)  ! do nothing
          case default; err=20; return
        end select  ! domain type
  
        ! -----
        ! - update volumetric fraction of liquid water and ice...
        !    => case of energy state or coupled solution (or adjusting the temperature)...
        ! --------------------------------------------------------------------------------
  
      ! case of energy state OR coupled solution (or adjusting the temperature)
      elseif(do_adjustTemp .or. ( (isNrgState .or. isCoupled) ) )then
  
        ! identify domain type
        select case(ixDomainType)
  
          ! *** vegetation canopy
          case(iname_veg)
  
            ! compute volumetric fraction of liquid water and ice
            call updateSnowPrime(&
                            xTemp,                                        & ! intent(in):  temperature (K)
                            scalarCanopyWatTrial(iGRU)/(iden_water*canopyDepth_(iGRU)),& ! intent(in):  volumetric fraction of total water (-)
                            snowfrz_scale,                                & ! intent(in):  scaling parameter for the snow freezing curve (K-1)
                            scalarCanopyTempPrime_(iGRU),                        & ! intent(in):  canopy temperature time derivative (K/s)
                            scalarCanopyWatPrime_(iGRU)/(iden_water*canopyDepth_(iGRU)),& ! intent(in):  volumetric fraction of total water time derivative (-)
                            scalarVolFracLiq,                             & ! intent(out): trial canopy liquid water (-)
                            scalarVolFracIce,                             & ! intent(out): trial volumetric canopy ice (-)
                            scalarVolFracLiqPrime,                        & ! intent(out): trial volumetric canopy liquid water (-)
                            scalarVolFracIcePrime,                        & ! intent(out): trial volumetric canopy ice (-)
                            scalarFracLiqVeg_(iGRU),                             & ! intent(out): fraction of liquid water (-)
                            err)                                   ! intent(out): error control
            if(err/=0)then; return; endif
  
            ! compute mass of water on the canopy
            ! NOTE: possibilities for speed-up here
            scalarCanopyLiqTrial_(iGRU) =             scalarFracLiqVeg_(iGRU) *scalarCanopyWatTrial(iGRU) !(kg m-2), scalarVolFracLiq*iden_water*canopyDepth
            scalarCanopyLiqPrime_(iGRU) =             scalarVolFracLiqPrime*iden_water*canopyDepth_(iGRU)
            scalarCanopyIceTrial_(iGRU) = (1._rkind - scalarFracLiqVeg_(iGRU))*scalarCanopyWatTrial(iGRU) !(kg m-2), scalarVolFracIce* iden_ice *canopyDepth
            scalarCanopyIcePrime_(iGRU) =             scalarVolFracIcePrime* iden_ice *canopyDepth_(iGRU)
  
          ! *** snow layers
          case(iname_snow)
  
            ! compute volumetric fraction of liquid water and ice
            call updateSnowPrime(&
                            xTemp,                          & ! intent(in):  temperature (K)
                            mLayerVolFracWatTrial_(iLayer,iGRU),  & ! intent(in):  mass state variable = trial volumetric fraction of water (-)
                            snowfrz_scale,                  & ! intent(in):  scaling parameter for the snow freezing curve (K-1)
                            mLayerTempPrime_(iLayer,iGRU),        & ! intent(in):  temperature time derivative (K/s)
                            mLayerVolFracWatPrime_in(iLayer,iGRU),  & ! intent(in):  volumetric fraction of total water time derivative (-)
                            mLayerVolFracLiqTrial_(iLayer,iGRU),  & ! intent(out): trial volumetric fraction of liquid water (-)
                            mLayerVolFracIceTrial_(iLayer,iGRU),  & ! intent(out): trial volumetric fraction if ice (-)
                            mLayerVolFracLiqPrime_(iLayer,iGRU),  & ! intent(out): volumetric fraction of liquid water time derivative (-)
                            mLayerVolFracIcePrime_(iLayer,iGRU),  & ! intent(out): volumetric fraction of ice time derivative (-)
                            mLayerFracLiqSnow_(iLayer,iGRU),      & ! intent(out): fraction of liquid water (-)
                            err)                     ! intent(out): error control
            if(err/=0)then; return; endif
  
          ! *** soil layers
          case(iname_soil)
          ! print*, 'update', iLayer, iGRU
  
            ! compute volumetric fraction of liquid water and ice
            call updateSoilPrime(&
                            xTemp,                                  & ! intent(in):  temperature (K)
                            mLayerMatricHeadTrial_(ixControlIndex,iGRU),  & ! intent(in):  total water matric potential (m)
                            mLayerTempPrime_(iLayer,iGRU),                & ! intent(in):  temperature time derivative (K/s)
                            mLayerMatricHeadPrime_(ixControlIndex,iGRU),  & ! intent(in):  total water matric potential time derivative (m/s)
                            vGn_alpha_(ixControlIndex),vGn_n_(ixControlIndex),theta_sat_(ixControlIndex),theta_res_(ixControlIndex),vGn_m_(ixControlIndex,iGRU), & ! intent(in): soil parameters
                            mLayerVolFracWatTrial_(iLayer,iGRU),          & ! intent(in):  mass state variable = trial volumetric fraction of water (-)
                            mLayerVolFracLiqTrial_(iLayer,iGRU),          & ! intent(out): trial volumetric fraction of liquid water (-)
                            mLayerVolFracIceTrial_(iLayer,iGRU),          & ! intent(out): trial volumetric fraction if ice (-)
                            mLayerVolFracWatPrime_in(iLayer,iGRU),          & ! intent(out): volumetric fraction of total water time derivative (-)
                            mLayerVolFracLiqPrime_(iLayer,iGRU),          & ! intent(out): volumetric fraction of liquid water time derivative (-)
                            mLayerVolFracIcePrime_(iLayer,iGRU),          & ! intent(out): volumetric fraction of ice time derivative (-)
                            err)                             ! intent(out): error control
            if(err/=0)then; return; endif
  
          ! check
          case default; err=20; return
  
        end select  ! domain type
  
      ! final check
      else
  
        ! do nothing (input = output) -- and check that we got here correctly
        if( (isNrgState .or. isCoupled) )then
          scalarVolFracLiq = realMissing
          scalarVolFracIce = realMissing
        else
          err=20; return
        endif
  
      endif  ! if energy state or solution is coupled
  
      ! -----
      ! ------------------------
  
      ! check the need to adjust temperature (will always be false if inside solver)
      !  can be true if inside varSubstep, outside solver, if in a splitting case
      !  NOTE: should be adjusting enthalpy if that is the state variabls
      if(do_adjustTemp)then
  
        ! get the melt energy
        meltNrg = merge(LH_fus*iden_ice, LH_fus*iden_water, ixDomainType==iname_snow)
  
        ! compute the residual and the derivative
        select case(ixDomainType)
  
          ! * vegetation
          case(iname_veg)
            call xTempSolve(&
                            ! constant over iterations
                            meltNrg         = meltNrg                                 ,&  ! intent(in):    energy for melt+freeze (J m-3)
                            heatCap         = scalarBulkVolHeatCapVeg_(iGRU)                 ,&  ! intent(in):    volumetric heat capacity (J m-3 K-1)
                            tempInit        = scalarCanopyTemp_(iGRU)                        ,&  ! intent(in):    initial temperature (K)
                            volFracIceInit  = scalarCanopyIce_(iGRU)/(iden_water*canopyDepth_(iGRU)),&  ! intent(in):    initial volumetric fraction of ice (-)
                            ! trial values
                            xTemp           = xTemp                                   ,&  ! intent(inout): trial value of temperature
                            dLiq_dT         = dTheta_dTkCanopy_(iGRU)                        ,&  ! intent(in):    derivative in liquid water content w.r.t. temperature (K-1)
                            volFracIceTrial = scalarVolFracIce                        ,&  ! intent(in):    trial value for volumetric fraction of ice
                            ! residual and derivative
                            residual        = residual                                ,&  ! intent(out):   residual (J m-3)
                            derivative      = derivative                               )  ! intent(out):   derivative (J m-3 K-1)
  
                ! * snow and soil
          case(iname_snow, iname_soil)
            call xTempSolve(&
                            ! constant over iterations
                            meltNrg         = meltNrg                        ,&  ! intent(in):    energy for melt+freeze (J m-3)
                            heatCap         = mLayerVolHtCapBulk_(iLayer,iGRU)     ,&  ! intent(in):    volumetric heat capacity (J m-3 K-1)
                            tempInit        = mLayerTemp_(iLayer,iGRU)             ,&  ! intent(in):    initial temperature (K)
                            volFracIceInit  = mLayerVolFracIce_(iLayer,iGRU)       ,&  ! intent(in):    initial volumetric fraction of ice (-)
                            ! trial values
                            xTemp           = xTemp                          ,&  ! intent(inout): trial value of temperature
                            dLiq_dT         = mLayerdTheta_dTk_(iLayer,iGRU)       ,&  ! intent(in):    derivative in liquid water content w.r.t. temperature (K-1)
                            volFracIceTrial = mLayerVolFracIceTrial_(iLayer,iGRU)  ,&  ! intent(in):    trial value for volumetric fraction of ice
                            ! residual and derivative
                            residual        = residual                       ,&  ! intent(out):   residual (J m-3)
                            derivative      = derivative                      )  ! intent(out):   derivative (J m-3 K-1)
  
          ! * check
          case default; err=20; return
  
        end select  ! domain type
  
        ! check validity of residual
        if( ieee_is_nan(residual) )then
          err=20; return
        endif
  
        ! update bracket
        if(residual < 0._rkind)then
          tempMax = min(xTemp,tempMax)
        else
          tempMin = max(tempMin,xTemp)
        end if
  
        ! compute iteration increment
        tempInc    = residual/derivative  ! K
  
        ! check
        ! if(globalPrintFlag)&
        ! write(*,'(i4,1x,e20.10,1x,5(f20.10,1x),L1)') iter, residual, xTemp-Tcrit, tempInc, Tcrit, tempMin, tempMax, bFlag
  
        ! check convergence
        if(abs(residual) < nrgConvTol .or. abs(tempInc) < tempConvTol) exit iterations
  
        ! add constraints for snow temperature
        if(ixDomainType==iname_veg .or. ixDomainType==iname_snow)then
          if(tempInc > Tcrit - xTemp) tempInc=(Tcrit - xTemp)*0.5_rkind  ! simple bi-section method
        endif  ! if the domain is vegetation or snow
  
        ! deal with the discontinuity between partially frozen and unfrozen soil
        if(ixDomainType==iname_soil)then
        ! difference from the temperature below which ice exists
          critDiff = Tcrit - xTemp
          ! --> initially frozen (T < Tcrit)
          if(critDiff > 0._rkind)then
            if(tempInc > critDiff) tempInc = critDiff + epsT  ! set iteration increment to slightly above critical temperature
          ! --> initially unfrozen (T > Tcrit)
          else
            if(tempInc < critDiff) tempInc = critDiff - epsT  ! set iteration increment to slightly below critical temperature
          endif
        endif  ! if the domain is soil
  
        ! update the temperature trial
        xTemp = xTemp + tempInc
  
        ! check failed convergence
        if(iter==maxiter)then
          err=-20; return ! negative error code = try to recover
        endif
  
      endif   ! if adjusting the temperature
  
    end do iterations ! iterating
  
    ! save temperature
    select case(ixDomainType)
      case(iname_veg);              scalarCanopyTempTrial_(iGRU)   = xTemp
      case(iname_snow, iname_soil); mLayerTempTrial_(iLayer,iGRU) = xTemp
    end select
  
    ! =======================================================================================================================================
    ! =======================================================================================================================================
  
    ! -----
    ! - compute the liquid water matric potential (and necessary derivatives)...
    ! -------------------------------------------------------------------------
  
    ! only for soil
    if(ixDomainType==iname_soil)then
  
      ! check liquid water
      if(mLayerVolFracLiqTrial_(iLayer,iGRU) > theta_sat_(ixControlIndex) )then
        err=20; return
      endif
  
      ! case of hydrology state uncoupled with energy
      if(.not.isNrgState .and. .not.isCoupled)then
  
        ! derivatives relating liquid water matric potential to total water matric potential and temperature
        dPsiLiq_dPsi0_(ixControlIndex,iGRU) = 1._rkind  ! exact correspondence (psiLiq=psi0)
        dPsiLiq_dTemp_(ixControlIndex,iGRU) = 0._rkind  ! no relationship between liquid water matric potential and temperature
  
      ! case of energy state or coupled solution
      else
        ! compute the liquid matric potential (and the derivatives w.r.t. total matric potential and temperature)
        call liquidHeadPrime(&
                        ! input
                        mLayerMatricHeadTrial_(ixControlIndex,iGRU)     ,& ! intent(in):  total water matric potential (m)
                        mLayerMatricHeadPrime_(ixControlIndex,iGRU)     ,& ! intent(in):  total water matric potential time derivative (m s-1)
                        mLayerVolFracLiqTrial_(iLayer,iGRU)             ,& ! intent(in):  volumetric fraction of liquid water (-)
                        mLayerVolFracIceTrial_(iLayer,iGRU)             ,& ! intent(in):  volumetric fraction of ice (-)
                        vGn_alpha_(ixControlIndex),vGn_n_(ixControlIndex),theta_sat_(ixControlIndex),theta_res_(ixControlIndex),vGn_m_(ixControlIndex,iGRU) ,& ! intent(in): soil parameters
                        dVolTot_dPsi0_(ixControlIndex,iGRU)             ,& ! intent(in):  derivative in the soil water characteristic (m-1)
                        mLayerdTheta_dTk_(iLayer,iGRU)                  ,& ! intent(in):  derivative in volumetric total water w.r.t. temperature (K-1)
                        mLayerVolFracLiqPrime_(iLayer,iGRU)             ,& ! intent(in):  volumetric fraction of liquid water time derivative (-)
                        mLayerVolFracIcePrime_(iLayer,iGRU)             ,& ! intent(in):  volumetric fraction of ice time derivative (-)
                        ! output
                        mLayerMatricHeadLiqTrial_(ixControlIndex,iGRU)  ,& ! intent(out): liquid water matric potential (m)
                        mLayerMatricHeadLiqPrime_(ixControlIndex,iGRU)  ,& ! intent(out): liquid water matric potential time derivative (m s-1)
                        dPsiLiq_dPsi0_(ixControlIndex,iGRU)             ,& ! intent(out): derivative in the liquid water matric potential w.r.t. the total water matric potential (-)
                        dPsiLiq_dTemp_(ixControlIndex,iGRU)             ,& ! intent(out): derivative in the liquid water matric potential w.r.t. temperature (m K-1)
                        err)                                ! intent(out): error control
        if(err/=0)then; return; endif
  
      endif  ! switch between hydrology and energy state
  
    endif  ! if domain is soil

  end subroutine
  
  ! **********************************************************************************************************
  ! private subroutine xTempSolve: compute residual and derivative for temperature
  ! **********************************************************************************************************
  attributes(device) subroutine xTempSolve(&
                        ! input: constant over iterations
                        meltNrg          ,&  ! intent(in):    energy for melt+freeze (J m-3)
                        heatCap          ,&  ! intent(in):    volumetric heat capacity (J m-3 K-1)
                        tempInit         ,&  ! intent(in):    initial temperature (K)
                        volFracIceInit   ,&  ! intent(in):    initial volumetric fraction of ice (-)
                        ! input-output: trial values
                        xTemp            ,&  ! intent(inout): trial value of temperature
                        dLiq_dT          ,&  ! intent(in):    derivative in liquid water content w.r.t. temperature (K-1)
                        volFracIceTrial  ,&  ! intent(in):    trial value for volumetric fraction of ice
                        ! output: residual and derivative
                        residual         ,&  ! intent(out):   residual (J m-3)
                        derivative        )  ! intent(out):   derivative (J m-3 K-1)
    implicit none
    ! input: constant over iterations
    real(rkind),intent(in)             :: meltNrg          ! energy for melt+freeze (J m-3)
    real(rkind),intent(in)             :: heatCap          ! volumetric heat capacity (J m-3 K-1)
    real(rkind),intent(in)             :: tempInit         ! initial temperature (K)
    real(rkind),intent(in)             :: volFracIceInit   ! initial volumetric fraction of ice (-)
    ! input-output: trial values
    real(rkind),intent(inout)          :: xTemp            ! trial value for temperature
    real(rkind),intent(in)             :: dLiq_dT          ! derivative in liquid water content w.r.t. temperature (K-1)
    real(rkind),intent(in)             :: volFracIceTrial  ! trial value for the volumetric fraction of ice (-)
    ! output: residual and derivative
    real(rkind),intent(out)            :: residual         ! residual (J m-3)
    real(rkind),intent(out)            :: derivative       ! derivative (J m-3 K-1)
    ! subroutine starts here
    residual   = -heatCap*(xTemp - tempInit) + meltNrg*(volFracIceTrial - volFracIceInit)  ! J m-3
    derivative = heatCap + LH_fus*iden_water*dLiq_dT  ! J m-3 K-1
  end subroutine xTempSolve
  
  end module updateVarsWithPrime_module
  