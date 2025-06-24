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

module computJacobWithPrime_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    model_options   ! defines the model decisions

! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! named variables to describe the state variable type
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers

! access named variables to describe the form and structure of the matrices used in the numerical solver
USE globalData,only: ku             ! number of super-diagonal bands, assume ku>=3
USE globalData,only: kl             ! number of sub-diagonal bands, assume kl>=4
USE globalData,only: ixDiag         ! index for the diagonal band
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! constants
USE multiconst,only:&
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:  &
 qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit                    ! no explicit groundwater parameterization

! look-up values for the form of Richards' equation
USE mDecisions_module,only:  &
 moisture,                   & ! moisture-based form of Richards' equation
 mixdform                      ! mixed form of Richards' equation

! look-up values for the choice of variable in energy equations (BE residual or IDA state variable)
USE mDecisions_module,only:  &
 closedForm,                 & ! use temperature with closed form heat capacity
 enthalpyFormLU,             & ! use enthalpy with soil temperature-enthalpy lookup tables
 enthalpyForm                  ! use enthalpy with soil temperature-enthalpy analytical solution

implicit none
! define constants
real(rkind),parameter     :: verySmall=tiny(1.0_rkind)     ! a very small number

private
public::computJacobWithPrime
public::computJacob4ida

contains


! **********************************************************************************************************
! public subroutine computJacobWithPrime: compute the Jacobian matrix
! **********************************************************************************************************
subroutine computJacobWithPrime(&
                      ! input: model control
                      cj,                         & ! intent(in):    this scalar changes whenever the step size or method order changes
                      dt,                         & ! intent(in):    length of the time step (seconds)
                      nSnow,                      & ! intent(in):    number of snow layers
                      nSoil,                      & ! intent(in):    number of soil layers
                      nLayers,                    & ! intent(in):    total number of layers
                      nGRU, &
                      computeVegFlux,             & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      computeBaseflow,            & ! intent(in):    flag to indicate if we need to compute baseflow
                      ixMatrix,                   & ! intent(in):    form of the Jacobian matrix
                      specificStorage,            & ! intent(in):    specific storage coefficient (m-1)
                      theta_sat,                  & ! intent(in):    soil porosity (-)
                      ixRichards,                 & ! intent(in):    choice of option for Richards' equation
                      enthalpyStateVec,           & ! intent(in):    flag if enthalpy is state variable
                      ! input: data structures
                      indx_data,                  & ! intent(in):    index data
                      prog_data,                  & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,                  & ! intent(in):    model diagnostic variables for a local HRU
                      deriv_data,                 & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                      dBaseflow_dMatric,          & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
                      ! input: state variables
                      mLayerTempPrime,            & ! intent(in):    vector of derivative value for layer temperature (K)
                      mLayerMatricHeadPrime,      & ! intent(in):    vector of derivative value for layer matric head
                      mLayerVolFracWatPrime,      & ! intent(in):    vector of derivative value for layer water volume fraction
                      scalarCanopyTempPrime,      & ! intent(in):    derivative value for temperature of the vegetation canopy (K)
                      scalarCanopyWatPrime,       & ! intent(in):    derivative value for water content of the vegetation canopy
                      ! input-output: Jacobian and its diagonal
                      dMat,                       & ! intent(inout): diagonal of the Jacobian matrix
                      aJac,                       & ! intent(out):   Jacobian matrix
                      ! output: error control
                      err,message)                  ! intent(out):   error code and error message
  ! -----------------------------------------------------------------------------------------------------------------
                      use ieee_arithmetic
                      use device_data_types
  implicit none
  ! input: model control
  real(rkind),intent(in)               :: cj
  real(rkind),intent(in)               :: dt                         ! length of the time step (seconds)
  integer(i4b),intent(in),device              :: nSnow(:)                      ! number of snow layers
  integer(i4b),intent(in)              :: nSoil                      ! number of soil layers
  integer(i4b),intent(in),device              :: nLayers(:)                    ! total number of layers in the snow+soil domain
  integer(i4b),intent(in) :: nGRU
  logical(lgt),intent(in)              :: computeVegFlux             ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)              :: computeBaseflow            ! flag to indicate if computing baseflow
  integer(i4b),intent(in)              :: ixMatrix                   ! form of the Jacobian matrix
  real(rkind),intent(in),device               :: specificStorage            ! specific storage coefficient (m-1)
  real(rkind),intent(in),device               :: theta_sat(:)               ! soil porosity (-)
  integer(i4b),intent(in),device              :: ixRichards                 ! choice of option for Richards' equation
  logical(lgt),intent(in)              :: enthalpyStateVec           ! flag if enthalpy is state variable
  ! input: data structures
  type(indx_data_device),intent(in)         :: indx_data                  ! indices defining model states and layers
  type(prog_data_device),intent(in)         :: prog_data                  ! prognostic variables for a local HRU
  type(diag_data_device),intent(in)         :: diag_data                  ! diagnostic variables for a local HRU
  type(deriv_data_device),intent(in)         :: deriv_data                 ! derivatives in model fluxes w.r.t. relevant state variables
  real(rkind),intent(in),device               :: dBaseflow_dMatric(:,:,:)     ! derivative in baseflow w.r.t. matric head (s-1)
  ! input: state variables
  real(rkind),intent(in),device               :: mLayerTempPrime(:,:)         ! vector of derivative value for layer temperature
  real(rkind),intent(in),device               :: mLayerMatricHeadPrime(:,:)   ! vector of derivative value for layer matric head
  real(rkind),intent(in),device               :: mLayerVolFracWatPrime(:,:)   ! vector of derivative value for layer water volume fraction
  real(rkind),intent(in),device               :: scalarCanopyTempPrime(:)      ! derivative value for temperature of the vegetation canopy (K)
  real(rkind),intent(in),device               :: scalarCanopyWatPrime(:)       ! derivative value for water content of the vegetation canopy
  ! input-output: Jacobian and its diagonal
  real(rkind),intent(inout),device            :: dMat(:,:)                    ! diagonal of the Jacobian matrix
  real(rkind),intent(out),device              :: aJac(:,:,:)                  ! Jacobian matrix
  ! output variables
  integer(i4b),intent(out)             :: err                        ! error code
  character(*),intent(out)             :: message                    ! error message
  ! --------------------------------------------------------------
  ! * local variables
  ! --------------------------------------------------------------
  ! indices of model state variables
  integer(i4b)                         :: jState          ! index of state within the state subset
  integer(i4b)                         :: qState          ! index of cross-derivative state variable for baseflow
  integer(i4b)                         :: nrgState        ! energy state variable
  integer(i4b)                         :: watState        ! hydrology state variable
  integer(i4b)                         :: nState          ! number of state variables
  integer(i4b),allocatable             :: nrgRows(:)      ! indices of rows for energy column in banded matrix
  integer(i4b),allocatable             :: watRows(:)      ! indices of rows for hydrology column in banded matrix
  ! indices of model layers
  integer(i4b)                         :: iLayer          ! index of model layer
  integer(i4b)                         :: jLayer          ! index of model layer within the full state vector (hydrology)
  integer(i4b)                         :: pLayer          ! indices of soil layers (used for the baseflow derivatives)
  ! conversion factors
  real(rkind)                          :: LH_fu0          ! latent heat of fusion, modified to be 0 if using enthalpy formulation and not using
  real(rkind)                          :: convLiq2tot     ! factor to convert liquid water derivative to total water derivative
  type(decisions_device) :: decisions
  integer(i4b) :: iGRU

  ! --------------------------------------------------------------
  ! associate variables from data structures
  associate(&
    ! indices of model state variables
    ixCasNrg                     => indx_data%ixCasNrg                       ,& ! intent(in): [i4b]    index of canopy air space energy state variable
    ixVegNrg                     => indx_data%ixVegNrg                       ,& ! intent(in): [i4b]    index of canopy energy state variable
    ixVegHyd                     => indx_data%ixVegHyd                       ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixTopNrg                     => indx_data%ixTopNrg                       ,& ! intent(in): [i4b]    index of upper-most energy state in the snow+soil subdomain
    ixTopHyd                     => indx_data%ixTopHyd                       ,& ! intent(in): [i4b]    index of upper-most hydrology state in the snow+soil subdomain
    ixAqWat                      => indx_data%ixAqWat                        ,& ! intent(in): [i4b]    index of water storage in the aquifer
    ! vectors of indices for specfic state types within specific sub-domains IN THE FULL STATE VECTOR
    ixNrgLayer                   => indx_data%ixNrgLayer                        ,& ! intent(in): [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
    ! vector of energy indices for the snow and soil domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnowSoilNrg_m                => indx_data%ixSnowSoilNrg                     ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    ixSnowOnlyNrg_m                => indx_data%ixSnowOnlyNrg                     ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow domain
    ixSoilOnlyNrg_m                => indx_data%ixSoilOnlyNrg                     ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the soil domain
    ! vector of hydrology indices for the snow and soil domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnowSoilHyd_m                => indx_data%ixSnowSoilHyd                     ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    ixSnowOnlyHyd_m                => indx_data%ixSnowOnlyHyd                     ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow domain
    ixSoilOnlyHyd_m                => indx_data%ixSoilOnlyHyd                     ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the soil domain
    ! number of state variables of a specific type
    ! nSnowSoilNrg                 => indx_data%nLayers                   ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
    ! nSnowOnlyNrg                 => indx_data%maxSnow                   ,& ! intent(in): [i4b]    number of energy state variables in the snow domain
    nSoilOnlyNrg                 => indx_data%nSoil                   ,& ! intent(in): [i4b]    number of energy state variables in the soil domain
    ! nSnowSoilHyd                 => indx_data%nLayers                   ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
    ! nSnowOnlyHyd                 => indx_data%maxSnow                   ,& ! intent(in): [i4b]    number of hydrology variables in the snow domain
    nSoilOnlyHyd                 => indx_data%nSoil                   ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
    ! type and index of model control volume
    ixHydType                    => indx_data%ixHydType                         ,& ! intent(in): [i4b(:)] index of the type of hydrology states in snow+soil domain
    ! mapping between states and model layers
    ixMapSubset2Full_m             => indx_data%ixMapSubset2Full                  ,& ! intent(in): [i4b(:)] list of indices in the full state vector that are in the state subset
    ixMapFull2Subset_m             => indx_data%ixMapFull2Subset                  ,& ! intent(in): [i4b(:)] list of indices in the state subset in each element of the full state vector
    ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
    dCanairNetFlux_dCanairTemp   => deriv_data%dCanairNetFlux_dCanairTemp    ,& ! intent(in): [dp]     derivative in net canopy air space flux w.r.t. canopy air temperature
    dCanairNetFlux_dCanopyTemp   => deriv_data%dCanairNetFlux_dCanopyTemp    ,& ! intent(in): [dp]     derivative in net canopy air space flux w.r.t. canopy temperature
    dCanairNetFlux_dGroundTemp   => deriv_data%dCanairNetFlux_dGroundTemp    ,& ! intent(in): [dp]     derivative in net canopy air space flux w.r.t. ground temperature
    dCanopyNetFlux_dCanairTemp   => deriv_data%dCanopyNetFlux_dCanairTemp    ,& ! intent(in): [dp]     derivative in net canopy flux w.r.t. canopy air temperature
    dCanopyNetFlux_dCanopyTemp   => deriv_data%dCanopyNetFlux_dCanopyTemp    ,& ! intent(in): [dp]     derivative in net canopy flux w.r.t. canopy temperature
    dCanopyNetFlux_dGroundTemp   => deriv_data%dCanopyNetFlux_dGroundTemp    ,& ! intent(in): [dp]     derivative in net canopy flux w.r.t. ground temperature
    dCanopyNetFlux_dCanWat       => deriv_data%dCanopyNetFlux_dCanWat        ,& ! intent(in): [dp]     derivative in net canopy fluxes w.r.t. canopy total water content
    dGroundNetFlux_dCanairTemp   => deriv_data%dGroundNetFlux_dCanairTemp    ,& ! intent(in): [dp]     derivative in net ground flux w.r.t. canopy air temperature
    dGroundNetFlux_dCanopyTemp   => deriv_data%dGroundNetFlux_dCanopyTemp    ,& ! intent(in): [dp]     derivative in net ground flux w.r.t. canopy temperature
    dGroundNetFlux_dCanWat       => deriv_data%dGroundNetFlux_dCanWat        ,& ! intent(in): [dp]     derivative in net ground fluxes w.r.t. canopy total water content
    ! derivatives in evaporative fluxes w.r.t. relevant state variables
    dCanopyEvaporation_dTCanair  => deriv_data%dCanopyEvaporation_dTCanair   ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. canopy air temperature
    dCanopyEvaporation_dTCanopy  => deriv_data%dCanopyEvaporation_dTCanopy   ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. canopy temperature
    dCanopyEvaporation_dTGround  => deriv_data%dCanopyEvaporation_dTGround   ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. ground temperature
    dCanopyEvaporation_dCanWat   => deriv_data%dCanopyEvaporation_dCanWat    ,& ! intent(in): [dp]     derivative in canopy evaporation w.r.t. canopy total water content
    dGroundEvaporation_dTCanair  => deriv_data%dGroundEvaporation_dTCanair   ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. canopy air temperature
    dGroundEvaporation_dTCanopy  => deriv_data%dGroundEvaporation_dTCanopy   ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. canopy temperature
    dGroundEvaporation_dTGround  => deriv_data%dGroundEvaporation_dTGround   ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. ground temperature
    dGroundEvaporation_dCanWat   => deriv_data%dGroundEvaporation_dCanWat    ,& ! intent(in): [dp]     derivative in ground evaporation w.r.t. canopy total water content
    ! derivatives in canopy water w.r.t canopy temperature
    dCanLiq_dTcanopy             => deriv_data%dCanLiq_dTcanopy              ,& ! intent(in): [dp]     derivative in canopy liquid storage w.r.t. temperature
    dTheta_dTkCanopy             => deriv_data%dTheta_dTkCanopy              ,& ! intent(in): [dp]     derivative in volumetric liquid water content w.r.t. temperature
    d2Theta_dTkCanopy2           => deriv_data%d2Theta_dTkCanopy2            ,& ! intent(in): [dp]     second derivative of volumetric liquid water content w.r.t. temperature
    dFracLiqVeg_dTkCanopy        => deriv_data%dFracLiqVeg_dTkCanopy         ,& ! intent(in): [dp]     derivative in fraction of (throughfall + drainage)  w.r.t. temperature
    ! derivatives in canopy liquid fluxes w.r.t. canopy water
    scalarCanopyLiqDeriv         => deriv_data%scalarCanopyLiqDeriv          ,& ! intent(in): [dp]     derivative in (throughfall + drainage) w.r.t. canopy liquid water
    ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. temperature in layers above and below
    dNrgFlux_dTempAbove          => deriv_data%dNrgFlux_dTempAbove_m              ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. temperature in the layer above
    dNrgFlux_dTempBelow          => deriv_data%dNrgFlux_dTempBelow_m              ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. temperature in the layer below
    ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. water state in layers above and below
    dNrgFlux_dWatAbove           => deriv_data%dNrgFlux_dWatAbove_m               ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. water state in the layer above
    dNrgFlux_dWatBelow           => deriv_data%dNrgFlux_dWatBelow_m               ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. water state in the layer below
    ! derivatives in soil transpiration w.r.t. canopy state variables
    mLayerdTrans_dTCanair        => deriv_data%mLayerdTrans_dTCanair_m            ,& ! intent(in): [dp(:)]  derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
    mLayerdTrans_dTCanopy        => deriv_data%mLayerdTrans_dTCanopy_m            ,& ! intent(in): [dp(:)]  derivatives in the soil layer transpiration flux w.r.t. canopy temperature
    mLayerdTrans_dTGround        => deriv_data%mLayerdTrans_dTGround_m            ,& ! intent(in): [dp(:)]  derivatives in the soil layer transpiration flux w.r.t. ground temperature
    mLayerdTrans_dCanWat         => deriv_data%mLayerdTrans_dCanWat_m             ,& ! intent(in): [dp(:)]  derivatives in the soil layer transpiration flux w.r.t. canopy total water
    ! derivatives in aquifer transpiration w.r.t. canopy state variables
    dAquiferTrans_dTCanair       => deriv_data%dAquiferTrans_dTCanair        ,& ! intent(in): [dp]     derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
    dAquiferTrans_dTCanopy       => deriv_data%dAquiferTrans_dTCanopy        ,& ! intent(in): [dp]     derivatives in the aquifer transpiration flux w.r.t. canopy temperature
    dAquiferTrans_dTGround       => deriv_data%dAquiferTrans_dTGround        ,& ! intent(in): [dp]     derivatives in the aquifer transpiration flux w.r.t. ground temperature
    dAquiferTrans_dCanWat        => deriv_data%dAquiferTrans_dCanWat         ,& ! intent(in): [dp]     derivatives in the aquifer transpiration flux w.r.t. canopy total water
    ! derivative in liquid water fluxes at the interface of snow layers w.r.t. volumetric liquid water content in the layer above
    iLayerLiqFluxSnowDeriv       => deriv_data%iLayerLiqFluxSnowDeriv_m           ,& ! intent(in): [dp(:)]  derivative in vertical liquid water flux at layer interfaces
    ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
    dVolTot_dPsi0                => deriv_data%dVolTot_dPsi0_m                    ,& ! intent(in): [dp(:)]  derivative in total water content w.r.t. total water matric potential
    d2VolTot_dPsi02              => deriv_data%d2VolTot_dPsi02_m                  ,& ! intent(in): [dp(:)]  second derivative in total water content w.r.t. total water matric potential
    dCompress_dPsi               => deriv_data%dCompress_dPsi_m                   ,& ! intent(in): [dp(:)]  derivative in compressibility w.r.t matric head
    dq_dHydStateAbove            => deriv_data%dq_dHydStateAbove_m                ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer above
    dq_dHydStateBelow            => deriv_data%dq_dHydStateBelow_m                ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer below
    dq_dHydStateLayerSurfVec     => deriv_data%dq_dHydStateLayerSurfVec_m         ,& ! intent(in): [dp(:)]  change in the flux in soil surface interface w.r.t. state variables in layers
    ! derivative in baseflow flux w.r.t. aquifer storage
    dBaseflow_dAquifer           => deriv_data%dBaseflow_dAquifer            ,& ! intent(in): [dp(:)]  derivative in baseflow flux w.r.t. aquifer storage (s-1)
    ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
    dq_dNrgStateAbove            => deriv_data%dq_dNrgStateAbove_m                ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer above
    dq_dNrgStateBelow            => deriv_data%dq_dNrgStateBelow_m                ,& ! intent(in): [dp(:)]  change in flux at layer interfaces w.r.t. states in the layer below
    dq_dNrgStateLayerSurfVec     => deriv_data%dq_dNrgStateLayerSurfVec_m         ,& ! intent(in): [dp(:)]  change in the flux in soil surface interface w.r.t. state variables in layers
    ! derivative in liquid water fluxes for the soil and snow domain w.r.t temperature
    dFracLiqSnow_dTk             => deriv_data%dFracLiqSnow_dTk_m                 ,& ! intent(in): [dp(:)]  derivative in fraction of liquid snow w.r.t. temperature
    mLayerdTheta_dTk             => deriv_data%mLayerdTheta_dTk_m                 ,& ! intent(in): [dp(:)]  derivative in volumetric liquid water content w.r.t. temperature
    mLayerd2Theta_dTk2           => deriv_data%mLayerd2Theta_dTk2_m               ,& ! intent(in): [dp(:)]  second derivative of volumetric liquid water content w.r.t. temperature
    ! derivative in bulk heat capacity w.r.t. relevant state variables
    dVolHtCapBulk_dPsi0          => deriv_data%dVolHtCapBulk_dPsi0_m           ,& ! intent(in): [dp(:)]  derivative in bulk heat capacity w.r.t. matric potential
    dVolHtCapBulk_dTheta         => deriv_data%dVolHtCapBulk_dTheta_m          ,& ! intent(in): [dp(:)]  derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dCanWat        => deriv_data%dVolHtCapBulk_dCanWat         ,& ! intent(in): [dp   ]  derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dTk            => deriv_data%dVolHtCapBulk_dTk_m             ,& ! intent(in): [dp(:)]  derivative in bulk heat capacity w.r.t. temperature
    dVolHtCapBulk_dTkCanopy      => deriv_data%dVolHtCapBulk_dTkCanopy       ,& ! intent(in): [dp   ]  derivative in bulk heat capacity w.r.t. temperature
    ! derivative in Cm w.r.t. relevant state variables
    dCm_dPsi0                    => deriv_data%dCm_dPsi0_m                     ,& ! intent(in): [dp(:)]  derivative in heat capacity w.r.t. matric potential (J kg-1)
    dCm_dTk                      => deriv_data%dCm_dTk_m                       ,& ! intent(in): [dp(:)]  derivative in heat capacity w.r.t. temperature (J kg-1 K-2)
    dCm_dTkCanopy                => deriv_data%dCm_dTkCanopy                 ,& ! intent(in): [dp   ]  derivative in heat capacity w.r.t. canopy temperature (J kg-1 K-2)
    ! derivatives of temperature if enthalpy is the state variable
    dCanairTemp_dEnthalpy        => deriv_data%dCanairTemp_dEnthalpy         ,& ! intent(in): [dp]     derivative of canopy air temperature w.r.t. enthalpy
    dCanopyTemp_dEnthalpy        => deriv_data%dCanopyTemp_dEnthalpy         ,& ! intent(in): [dp]     derivative of canopy temperature w.r.t. enthalpy 
    dTemp_dEnthalpy              => deriv_data%dTemp_dEnthalpy_m               ,& ! intent(in): [dp(:)]  derivative of temperature w.r.t. enthalpy
    dCanopyTemp_dCanWat          => deriv_data%dCanopyTemp_dCanWat           ,& ! intent(in): [dp]     derivative of canopy temperature w.r.t. volumetric water content
    dTemp_dTheta                 => deriv_data%dTemp_dTheta_m                  ,& ! intent(in): [dp(:)]  derivative of temperature w.r.t. volumetric water content
    dTemp_dPsi0                  => deriv_data%dTemp_dPsi0_m                   ,& ! intent(in): [dp(:)]  derivative of temperature w.r.t. total water matric potential
    ! diagnostic variables
    scalarFracLiqVeg             => diag_data%scalarFracLiqVeg                ,& ! intent(in): [dp]     fraction of liquid water on vegetation (-)
    scalarBulkVolHeatCapVeg      => diag_data%scalarBulkVolHeatCapVeg         ,& ! intent(in): [dp]     bulk volumetric heat capacity of vegetation (J m-3 K-1)
    scalarCanopyCm               => diag_data%scalarCanopyCm                  ,& ! intent(in): [dp]     Cm of canopy (J kg-1 K-1)
    mLayerFracLiqSnow            => diag_data%mLayerFracLiqSnow_m               ,& ! intent(in): [dp(:)]  fraction of liquid water in each snow layer (-)
    mLayerVolHtCapBulk           => diag_data%mLayerVolHtCapBulk_m              ,& ! intent(in): [dp(:)]  bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
    mLayerCm                     => diag_data%mLayerCm_m                        ,& ! intent(in): [dp(:)]  Cm in each snow and soil layer (J kg-1 K-1)
    scalarSoilControl            => diag_data%scalarSoilControl               ,& ! intent(in): [dp]     soil control on infiltration, zero or one
    ! canopy and layer depth
    canopyDepth                  => diag_data%scalarCanopyDepth               ,& ! intent(in): [dp   ]  canopy depth (m)
    mLayerDepth                  => prog_data%mLayerDepth                        ,& ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    layerType_m                    => indx_data%layerType                          & ! intent(in): [i4b(:)] named variables defining the type of layer in snow+soil domain
    ) ! making association with data in structures
    ! --------------------------------------------------------------
    ! initialize error control


    err=0; message='computJacobWithPrime/'

    ! *********************************************************************************************************************************************************
    ! * PART 0: PRELIMINARIES (INITIALIZE JACOBIAN AND COMPUTE TIME-VARIABLE DIAGONAL TERMS)
    ! *********************************************************************************************************************************************************

    ! get the number of state variables
    nState = size(dMat,1)

    ! initialize the Jacobian
    ! NOTE: this needs to be done every time, since Jacobian matrix is modified in the solver
    aJac = 0._rkind  ! analytical Jacobian matrix

    ! compute terms in the Jacobian for vegetation (excluding fluxes)
    ! NOTE: energy for vegetation is computed *within* the iteration loop as it includes phase change
      ! print*, nSnow, nSoil
     !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
    ! compute terms in the Jacobian for vegetation (excluding fluxes)
    ! NOTE: energy for vegetation is computed *within* the iteration loop as it includes phase change
    if(ixVegNrg(iGRU)/=integerMissing)then
      dMat(ixVegNrg(iGRU),iGRU) = ( scalarBulkVolHeatCapVeg(iGRU) + LH_fus*iden_water*dTheta_dTkCanopy(iGRU) ) * cj &
                      + dVolHtCapBulk_dTkCanopy(iGRU) * scalarCanopyTempPrime(iGRU) &
                      + dCm_dTkCanopy(iGRU) * scalarCanopyWatPrime(iGRU) / canopyDepth(iGRU) &
                      + LH_fus*iden_water * scalarCanopyTempPrime(iGRU) * d2Theta_dTkCanopy2(iGRU) &
                      + LH_fus            * dFracLiqVeg_dTkCanopy(iGRU) * scalarCanopyWatPrime(iGRU) / canopyDepth(iGRU)
    endif
  end do
    ! compute additional terms for the Jacobian for the snow-soil domain (excluding fluxes)
    ! NOTE: energy for snow+soil is computed *within* the iteration loop as it includes phase change
  !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    do iLayer=1,nLayers(iGRU)
      if(ixSnowSoilNrg_m(iLayer,iGRU)/=integerMissing)then
        dMat(ixSnowSoilNrg_m(iLayer,iGRU),iGRU) = ( mLayerVolHtCapBulk(iLayer,iGRU) + LH_fus*iden_water*mLayerdTheta_dTk(iLayer,iGRU) ) * cj &
                                    + dVolHtCapBulk_dTk(iLayer,iGRU) * mLayerTempPrime(iLayer,iGRU) &
                                    + dCm_dTk(iLayer,iGRU) * mLayerVolFracWatPrime(iLayer,iGRU) &
                                    + LH_fus*iden_water * mLayerTempPrime(iLayer,iGRU)  * mLayerd2Theta_dTk2(iLayer,iGRU) &
                                    + LH_fus*iden_water * dFracLiqSnow_dTk(iLayer,iGRU) * mLayerVolFracWatPrime(iLayer,iGRU)
      endif
    end do
  end do

    ! compute additional terms for the Jacobian for the soil domain (excluding fluxes)
  !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    do iLayer=1,nSoil
      if(ixSoilOnlyHyd_m(iLayer,iGRU)/=integerMissing)then
        dMat(ixSoilOnlyHyd_m(iLayer,iGRU),iGRU) = ( dVolTot_dPsi0(iLayer,iGRU) + dCompress_dPsi(iLayer,iGRU) ) * cj + d2VolTot_dPsi02(iLayer,iGRU) * mLayerMatricHeadPrime(iLayer,iGRU)

      if(ixRichards==mixdform)then
        dMat(ixSoilOnlyHyd_m(iLayer,iGRU),iGRU) = dMat(ixSoilOnlyHyd_m(iLayer,iGRU),iGRU) + specificStorage * dVolTot_dPsi0(iLayer,iGRU) * mLayerMatricHeadPrime(iLayer,iGRU) / theta_sat(iLayer)
      endif

      endif
    end do
  end do

    ! if using enthalpy as a state variable, zero out usual RHS terms and add them end of the iteration loop
    if(enthalpyStateVec)then 
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      if(ixCasNrg(iGRU)/=integerMissing) dMat(ixCasNrg(iGRU),iGRU) = 0._rkind
      if(ixVegNrg(iGRU)/=integerMissing) dMat(ixVegNrg(iGRU),iGRU) = 0._rkind
      end do
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      do iLayer=1,nLayers(iGRU)
        if(ixSnowSoilNrg_m(iLayer,iGRU)/=integerMissing) dMat(ixSnowSoilNrg_m(iLayer,iGRU),iGRU) = 0._rkind
      end do
    end do
      LH_fu0 = 0._rkind ! set to 0 to not use RHS terms
    else
      LH_fu0 = LH_fus ! use regular value
    endif

      ! *********************************************************************************************************************************************************
      ! * PART 2: FULL MATRIX
      ! *********************************************************************************************************************************************************

        ! -----
        ! * energy and liquid fluxes over vegetation...
        ! ---------------------------------------------
        if(computeVegFlux)then  ! (derivatives only defined when vegetation protrudes over the surface)

          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
          ! * energy fluxes with the canopy water
          if(ixVegHyd(iGRU)/=integerMissing)then

            ! * cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
            if(ixCasNrg(iGRU)/=integerMissing) aJac(ixVegHyd(iGRU),ixCasNrg(iGRU),iGRU) = -dCanopyEvaporation_dTCanair(iGRU)*dt
            ! dt*scalarCanopyLiqDeriv*dCanLiq_dTcanopy is the derivative in throughfall and canopy drainage with canopy temperature
            if(ixVegNrg(iGRU)/=integerMissing) aJac(ixVegHyd(iGRU),ixVegNrg(iGRU),iGRU) = -dCanopyEvaporation_dTCanopy(iGRU)*dt + dt*scalarCanopyLiqDeriv(iGRU)*dCanLiq_dTcanopy(iGRU)
            ! * liquid water fluxes for vegetation canopy (-), dt*scalarFracLiqVeg*scalarCanopyLiqDeriv is the derivative in throughfall and canopy drainage with canopy water
                                          aJac(ixVegHyd(iGRU),ixVegHyd(iGRU),iGRU) = -scalarFracLiqVeg(iGRU)*(dCanopyEvaporation_dCanWat(iGRU) - scalarCanopyLiqDeriv(iGRU))*dt + 1._rkind * cj
            if(ixTopNrg(iGRU)/=integerMissing) aJac(ixVegHyd(iGRU),ixTopNrg(iGRU),iGRU) = -dCanopyEvaporation_dTGround(iGRU)*dt

            ! * cross-derivative terms w.r.t. canopy water (kg-1 m2)
            if(ixTopHyd(iGRU)/=integerMissing) aJac(ixTopHyd(iGRU),ixVegHyd(iGRU),iGRU) = (dt/mLayerDepth(1,iGRU))*(-scalarSoilControl(iGRU)*scalarFracLiqVeg(iGRU)*scalarCanopyLiqDeriv(iGRU))/iden_water

            ! * cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
            ! NOTE: dIce/dLiq = (1 - scalarFracLiqVeg); dIce*LH_fu0/canopyDepth = J m-3; dLiq = kg m-2
            if(ixVegNrg(iGRU)/=integerMissing) aJac(ixVegNrg(iGRU),ixVegHyd(iGRU),iGRU) = (-1._rkind + scalarFracLiqVeg(iGRU))*LH_fu0/canopyDepth(iGRU) * cj &
                                                                  + dVolHtCapBulk_dCanWat(iGRU) * scalarCanopyTempPrime(iGRU) + scalarCanopyCm(iGRU)/canopyDepth(iGRU) * cj&
                                                                  - (dt/canopyDepth(iGRU)) * dCanopyNetFlux_dCanWat(iGRU) &
                                                                  + LH_fu0 * scalarCanopyTempPrime(iGRU) * dFracLiqVeg_dTkCanopy(iGRU) / canopyDepth(iGRU)
            if(ixTopNrg(iGRU)/=integerMissing) aJac(ixTopNrg(iGRU),ixVegHyd(iGRU),iGRU) = (dt/mLayerDepth(1,iGRU))*(-dGroundNetFlux_dCanWat(iGRU))
          endif

          ! * cross-derivative terms w.r.t. canopy temperature (K-1)
          if(ixVegNrg(iGRU)/=integerMissing)then
            if(ixTopHyd(iGRU)/=integerMissing) aJac(ixTopHyd(iGRU),ixVegNrg(iGRU),iGRU) = (dt/mLayerDepth(1,iGRU))*(-scalarSoilControl(iGRU)*scalarCanopyLiqDeriv(iGRU)*dCanLiq_dTcanopy(iGRU))/iden_water
          endif

          ! * energy fluxes with the canopy air space (J m-3 K-1)
          if(ixCasNrg(iGRU)/=integerMissing)then
                                         aJac(ixCasNrg(iGRU),ixCasNrg(iGRU),iGRU) = (dt/canopyDepth(iGRU))*(-dCanairNetFlux_dCanairTemp(iGRU)) + dMat(ixCasNrg(iGRU),iGRU) * cj
            if(ixVegNrg(iGRU)/=integerMissing) aJac(ixCasNrg(iGRU),ixVegNrg(iGRU),iGRU) = (dt/canopyDepth(iGRU))*(-dCanairNetFlux_dCanopyTemp(iGRU))
            if(ixTopNrg(iGRU)/=integerMissing) aJac(ixCasNrg(iGRU),ixTopNrg(iGRU),iGRU) = (dt/canopyDepth(iGRU))*(-dCanairNetFlux_dGroundTemp(iGRU))
          endif

          ! * energy fluxes with the vegetation canopy (J m-3 K-1)
          if(ixVegNrg(iGRU)/=integerMissing)then
            if(ixCasNrg(iGRU)/=integerMissing) aJac(ixVegNrg(iGRU),ixCasNrg(iGRU),iGRU) = (dt/canopyDepth(iGRU))*(-dCanopyNetFlux_dCanairTemp(iGRU))
            aJac(ixVegNrg(iGRU),ixVegNrg(iGRU),iGRU) = (dt/canopyDepth(iGRU))*(-dCanopyNetFlux_dCanopyTemp(iGRU)) + dMat(ixVegNrg(iGRU),iGRU)
            if(ixTopNrg(iGRU)/=integerMissing) aJac(ixVegNrg(iGRU),ixTopNrg(iGRU),iGRU) = (dt/canopyDepth(iGRU))*(-dCanopyNetFlux_dGroundTemp(iGRU))
          endif

          ! * energy fluxes with the surface (J m-3 K-1)
          if(ixTopNrg(iGRU)/=integerMissing)then
            if(ixCasNrg(iGRU)/=integerMissing) aJac(ixTopNrg(iGRU),ixCasNrg(iGRU),iGRU) = (dt/mLayerDepth(1,iGRU))*(-dGroundNetFlux_dCanairTemp(iGRU))
            if(ixVegNrg(iGRU)/=integerMissing) aJac(ixTopNrg(iGRU),ixVegNrg(iGRU),iGRU) = (dt/mLayerDepth(1,iGRU))*(-dGroundNetFlux_dCanopyTemp(iGRU))
          endif
        end do

        endif  ! if there is a need to compute energy fluxes within vegetation
    
        ! -----
        ! * energy fluxes for the snow+soil domain...
        ! -------------------------------------------
        ! if(nSnowSoilNrg>0)then
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
          do iLayer=1,nLayers(iGRU)  ! loop through all layers in the snow+soil domain

            ! check if the state is in the subset
            if(ixSnowSoilNrg_m(iLayer,iGRU)==integerMissing) cycle

            ! - define index within the state subset and the full state vector
            jState = ixSnowSoilNrg_m(iLayer,iGRU)        ! index within the state subset

            ! - diagonal elements
            aJac(jState,jState,iGRU)   = (dt/mLayerDepth(iLayer,iGRU))*(-dNrgFlux_dTempBelow(iLayer-1,iGRU) + dNrgFlux_dTempAbove(iLayer,iGRU)) + dMat(jState,iGRU)

            ! - lower-diagonal elements
            if(iLayer>1)then
              if(ixSnowSoilNrg_m(iLayer-1,iGRU)/=integerMissing) aJac(ixSnowSoilNrg_m(iLayer-1,igRU),jState,iGRU) = (dt/mLayerDepth(iLayer-1,iGRU))*( dNrgFlux_dTempBelow(iLayer-1,iGRU) )
            endif

            ! - upper diagonal elements
            if(iLayer<nLayers(iGRU))then
              if(ixSnowSoilNrg_m(iLayer+1,iGRU)/=integerMissing) aJac(ixSnowSoilNrg_m(iLayer+1,iGRU),jState,iGRU) = (dt/mLayerDepth(iLayer+1,iGRU))*(-dNrgFlux_dTempAbove(iLayer,iGRU  ) )
            endif

          end do  ! (looping through energy states in the snow+soil domain)
        end do
        ! endif   ! (if the subset includes energy state variables in the snow+soil domain)
        ! -----
        ! * liquid water fluxes for the snow domain...
        ! --------------------------------------------
        ! if(nSnowOnlyHyd>0)then

          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
          do iLayer=1,nSnow(iGRU)  ! loop through layers in the snow domain
            ! - check that the snow layer is desired
            if(ixSnowOnlyHyd_m(iLayer,iGRU)==integerMissing) cycle

            ! - define state indices for the current layer
            watState = ixSnowOnlyHyd_m(iLayer,iGRU)   ! hydrology state index within the state subset

            ! compute factor to convert liquid water derivative to total water derivative
            select case( ixHydType(iLayer,iGRU) )
              case(iname_watLayer); convLiq2tot = mLayerFracLiqSnow(iLayer,iGRU)
              case default;         convLiq2tot = 1._rkind
            end select

            ! - diagonal elements
            aJac(watState,watState,iGRU) = (dt/mLayerDepth(iLayer,iGRU))*iLayerLiqFluxSnowDeriv(iLayer,iGRU)*convLiq2tot + dMat(watState,iGRU) * cj
            ! print*, 673, iLayer,iGRU,aJac(watState,watState,iGRU)

            ! - lower-diagonal elements
            if(iLayer>1)then
              if(ixSnowOnlyHyd_m(iLayer-1,iGRU)/=integerMissing) aJac(ixSnowOnlyHyd_m(iLayer-1,iGRU),watState,iGRU) = 0._rkind  ! sub-diagonal: no dependence on other layers
              ! print*, 678, iLayer,iGRU,aJac(ixSnowOnlyHyd(iLayer-1),watState,iGRU)
            endif

            ! - upper diagonal elements
            if(iLayer<nSnow(iGRU))then
              if(ixSnowOnlyHyd_m(iLayer+1,iGRU)/=integerMissing) aJac(ixSnowOnlyHyd_m(iLayer+1,iGRU),watState,iGRU) = -(dt/mLayerDepth(iLayer+1,iGRU))*iLayerLiqFluxSnowDeriv(iLayer,iGRU)*convLiq2tot       ! dVol(below)/dLiq(above) -- (-)
              ! print*, 684, iLayer,iGRU, aJac(ixSnowOnlyHyd(iLayer+1),watState,iGRU)
            endif

          end do  ! (looping through liquid water states in the snow domain)
        end do
        ! endif   ! (if the subset includes hydrology state variables in the snow domain)
          ! print*, 685, maxval(aJac), minval(aJac)
                  ! -----
        ! * cross derivatives in the snow domain...
        ! ----------------------------------------
        ! if(nSnowOnlyHyd>0 .and. nSnowOnlyNrg>0)then
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
          do iLayer=1,nSnow(iGRU)  ! loop through layers in the snow domain

            ! - check that the snow layer is desired
            if(ixSnowOnlyNrg_m(iLayer,iGRU)==integerMissing) cycle

            ! (define the energy state)
            nrgState = ixSnowOnlyNrg_m(iLayer,iGRU)       ! index within the full state vector

            ! - define state indices for the current layer
            watState = ixSnowOnlyHyd_m(iLayer,iGRU)   ! hydrology state index within the state subset

            if(watstate/=integerMissing)then       ! (energy state for the current layer is within the state subset)

              ! - include derivatives of energy fluxes w.r.t water fluxes for current layer
              aJac(nrgState,watState,iGRU) = (-1._rkind + mLayerFracLiqSnow(iLayer,iGRU))*LH_fu0*iden_water * cj &
                                          + dVolHtCapBulk_dTheta(iLayer,iGRU) * mLayerTempPrime(iLayer,iGRU) + mLayerCm(iLayer,iGRU) * cj &
                                          + (dt/mLayerDepth(iLayer,iGRU))*(-dNrgFlux_dWatBelow(iLayer-1,iGRU) + dNrgFlux_dWatAbove(iLayer,iGRU)) &
                                          + LH_fu0*iden_water * mLayerTempPrime(iLayer,iGRU) * dFracLiqSnow_dTk(iLayer,iGRU)    ! (dF/dLiq)
                                          ! print*, 715, iLayer, iGRU, aJac(nrgState,watState,iGRU)

              ! - include derivatives of water fluxes w.r.t energy fluxes for current layer
              aJac(watState,nrgState,iGRU) = (dt/mLayerDepth(iLayer,iGRU))*iLayerLiqFluxSnowDeriv(iLayer,iGRU)*mLayerdTheta_dTk(iLayer,iGRU)  ! (dVol/dT)
              ! print*, 719, iLayer, iGRU, aJac(watState,nrgState,iGRU)
              ! (cross-derivative terms for the layer below)
              if(iLayer<nSnow(iGRU))then
                if(ixSnowOnlyHyd_m(iLayer+1,iGRU)/=integerMissing) aJac(ixSnowOnlyHyd_m(iLayer+1,iGRU),nrgState,iGRU) = -(dt/mLayerDepth(iLayer+1,iGRU))*iLayerLiqFluxSnowDeriv(iLayer,iGRU)*mLayerdTheta_dTk(iLayer,iGRU)    ! dVol(below)/dT(above) -- K-1
                ! print*, 723, iLayer,iGRU, aJac(ixSnowOnlyHyd(iLayer+1),nrgState,iGRU)
              endif ! (if there is a water state in the layer below the current layer in the given state subset)

              ! - include derivatives of heat capacity w.r.t water fluxes for surrounding layers starting with layer above
              if(iLayer>1)then
                if(ixSnowOnlyNrg_m(iLayer-1,iGRU)/=integerMissing) aJac(ixSnowOnlyNrg_m(iLayer-1,iGRU),watState,iGRU) = (dt/mLayerDepth(iLayer-1,iGRU))*( dNrgFlux_dWatBelow(iLayer-1,iGRU) )
                ! print*, iLayer,iGRU, aJac(ixSnowOnlyNrg(iLayer-1),watState,iGRU)
              endif

              ! (cross-derivative terms for the layer below)
              if(iLayer<nSnow(iGRU))then
                if(ixSnowOnlyNrg_m(iLayer+1,iGRU)/=integerMissing) aJac(ixSnowOnlyNrg_m(iLayer+1,iGRU),watState,iGRU) = (dt/mLayerDepth(iLayer+1,iGRU))*(-dNrgFlux_dWatAbove(iLayer  ,iGRU) )
                ! print*, iLayer,iGRU,aJac(ixSnowOnlyNrg(iLayer+1),watState,iGRU)
              elseif(iLayer==nSnow(iGRU) .and. nSoil>0)then !bottom snow layer and there is soil below
                if(ixSoilOnlyNrg_m(1,iGRU)/=integerMissing) aJac(ixSoilOnlyNrg_m(1,iGRU),watState,iGRU) = (dt/mLayerDepth(nSnow(iGRU)+1,iGRU))*(-dNrgFlux_dWatAbove(nSnow(iGRU),iGRU) )
                ! print*, iLayer,iGRU, aJac(ixSoilOnlyNrg(1),watState,iGRU)
              endif

            endif   ! (if the energy state for the current layer is within the state subset)

          end do  ! (looping through snow layers)
        end do
        ! endif   ! (if there are state variables for both water and energy in the snow domain)

        ! -----
        ! * liquid water fluxes for the soil domain...
        ! --------------------------------------------
        if(nSoilOnlyHyd>0)then
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
          do iLayer=1,nSoil

            ! - check that the soil layer is desired
            if(ixSoilOnlyHyd_m(iLayer,iGRU)==integerMissing) cycle

            ! - define state indices
            watState = ixSoilOnlyHyd_m(iLayer,iGRU)         ! hydrology state index within the state subset

            ! - define indices of the soil layers
            jLayer   = iLayer+nSnow(iGRU)                  ! index of layer in the snow+soil vector

            ! - compute the diagonal elements
            ! all terms *excluding* baseflow
            aJac(watState,watState,iGRU) = (dt/mLayerDepth(jLayer,iGRU))*(-dq_dHydStateBelow(iLayer-1,iGRU) + dq_dHydStateAbove(iLayer,iGRU)) + dMat(watState,iGRU)

            ! - compute the lower-diagonal elements
            if(iLayer>1)then
              if(ixSoilOnlyHyd_m(iLayer-1,iGRU)/=integerMissing) aJac(ixSoilOnlyHyd_m(iLayer-1,iGRU),watState,iGRU) = (dt/mLayerDepth(jLayer-1,iGRU))*( dq_dHydStateBelow(iLayer-1,iGRU))
            endif

            ! - compute the upper-diagonal elements
            if(iLayer<nSoil)then
              if(ixSoilOnlyHyd_m(iLayer+1,iGRU)/=integerMissing) aJac(ixSoilOnlyHyd_m(iLayer+1,iGRU),watState,iGRU) = (dt/mLayerDepth(jLayer+1,iGRU))*(-dq_dHydStateAbove(iLayer,iGRU))
            endif
            ! print*, iLayer, mLayerDepth(1+nSnow)

            ! if(computeBaseflow) then ! .and. nSoilOnlyHyd==nSoil)then
            !   do pLayer=1,nSoil
            !   qState = ixSoilOnlyHyd_m(pLayer,iGRU)  ! hydrology state index within the state subset
            !   if(qState/=integerMissing) aJac(watState,qState,iGRU) = (dt/mLayerDepth(jLayer,iGRU))*dBaseflow_dMatric(iLayer,pLayer,iGRU) + aJac(watState,qState,iGRU)
            !   end do
            ! endif

            ! - include terms for surface infiltration below surface
            ! if(ixSoilOnlyHyd_m(1,iGRU)/=integerMissing) aJac(ixSoilOnlyHyd_m(1,iGRU),watState,iGRU) = -(dt/mLayerDepth(1+nSnow,iGRU))*dq_dHydStateLayerSurfVec(iLayer,iGRU) + aJac(ixSoilOnlyHyd_m(1,iGRU),watState,iGRU)

          end do  ! (looping through hydrology states in the soil domain)
        end do


                  !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
          do iLayer=1,nSoil
            jLayer=nSnow(iGRU)

            ! - check that the soil layer is desired
            if(ixSoilOnlyHyd_m(iLayer,iGRU)==integerMissing) cycle

            ! - define state indices
            watState = ixSoilOnlyHyd_m(iLayer,iGRU)         ! hydrology state index within the state subset

            ! - define indices of the soil layers
            ! jLayer   = iLayer+nSnow                  ! index of layer in the snow+soil vector

            ! - include terms for surface infiltration below surface
            if(ixSoilOnlyHyd_m(1,iGRU)/=integerMissing) then
              ! print*, iLayer, iGRU, mLayerDepth(1+jLayer,iGRU)
              aJac(ixSoilOnlyHyd_m(1,iGRU),watState,iGRU) = -(dt/mLayerDepth(1+nSnow(iGRU),iGRU))*dq_dHydStateLayerSurfVec(iLayer,iGRU) + aJac(ixSoilOnlyHyd_m(1,iGRU),watState,iGRU)
            end if

          end do  ! (looping through hydrology states in the soil domain)
      end do

        ! err = cudaDeviceSynchronize()
        ! print*, 805, err

        ! - include terms for baseflow
        if (computeBaseflow) then! .and. nSoilOnlyHyd==nSoil) then
          ! print*, 807
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
            do iLayer=1,nSoil
              do pLayer=1,nSoil
            ! - check that the soil layer is desired
              if(ixSoilOnlyHyd_m(iLayer,iGRU)==integerMissing) cycle

              ! - define state indices
              watState = ixSoilOnlyHyd_m(iLayer,iGRU)         ! hydrology state index within the state subset
  
              ! - define indices of the soil layers
              jLayer   = iLayer+nSnow(iGRU)                  ! index of layer in the snow+soil vector
  
                qState = ixSoilOnlyHyd_m(pLayer,iGRU)  ! hydrology state index within the state subset
                if(qState/=integerMissing) aJac(watState,qState,iGRU) = (dt/mLayerDepth(jLayer,iGRU))*dBaseflow_dMatric(iLayer,pLayer,iGRU) + aJac(watState,qState,iGRU)
              end do
            end do
          end do
        end if

        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
          ! - include terms for surface infiltration above surface
          if(nSnow(iGRU)>0)then !have snow above first soil layer
            if(ixSnowOnlyHyd_m(nSnow(iGRU),iGRU)/=integerMissing .and. ixSoilOnlyHyd_m(1,iGRU)/=integerMissing) aJac(ixSoilOnlyHyd_m(1,iGRU),ixSnowOnlyHyd_m(nSnow(iGRU),iGRU),iGRU) = -(dt/mLayerDepth(1+nSnow(iGRU),iGRU))*dq_dHydStateLayerSurfVec(0,iGRU)
            ! print*, iGRU, aJac(ixSoilOnlyHyd_m(1,iGRU),ixSnowOnlyHyd(nSnow),iGRU)
          elseif(computeVegFlux)then !have vegetation above first soil layer, ixTopHyd = ixSoilOnlyHyd_m(1,iGRU)
            if(ixVegHyd(iGRU)/=integerMissing .and. ixTopHyd(iGRU)/=integerMissing) aJac(ixTopHyd(iGRU),ixVegHyd(iGRU),iGRU) = -(dt/mLayerDepth(1+nSnow(iGRU),iGRU))*dq_dHydStateLayerSurfVec(0,iGRU) + aJac(ixTopHyd(iGRU),ixVegHyd(iGRU),iGRU)
          endif
        end do

        endif   ! (if the subset includes hydrology state variables in the soil domain)
        
        ! -----
        ! * liquid water fluxes for the aquifer...
        ! ----------------------------------------
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU

          if(ixAqWat(iGRU)/=integerMissing) then
            aJac(ixAqWat(iGRU),ixAqWat(iGRU),iGRU) = -dBaseflow_dAquifer(iGRU)*dt + dMat(ixAqWat(iGRU),iGRU) * cj
            if(ixSoilOnlyNrg_m(nSoil,iGRU)/=integerMissing) aJac(ixAqWat(iGRU),ixSoilOnlyNrg_m(nSoil,iGRU),iGRU) = -dq_dNrgStateAbove(nSoil,iGRU)*dt ! dAquiferRecharge_dTk  = d_iLayerLiqFluxSoil(nSoil)_dTk
            if(ixSoilOnlyHyd_m(nSoil,iGRU)/=integerMissing) aJac(ixAqWat(iGRU),ixSoilOnlyHyd_m(nSoil,iGRU),iGRU) = -dq_dHydStateAbove(nSoil,iGRU)*dt ! dAquiferRecharge_dWat = d_iLayerLiqFluxSoil(nSoil)_dWat
            ! - include derivatives of energy and water w.r.t soil transpiration (dependent on canopy transpiration)
            if(computeVegFlux)then
              if(ixCasNrg(iGRU)/=integerMissing) aJac(ixAqWat(iGRU),ixCasNrg(iGRU),iGRU) = -dAquiferTrans_dTCanair(iGRU)*dt ! dVol/dT (K-1)
              if(ixVegNrg(iGRU)/=integerMissing) aJac(ixAqWat(iGRU),ixVegNrg(iGRU),iGRU) = -dAquiferTrans_dTCanopy(iGRU)*dt ! dVol/dT (K-1)
              if(ixVegHyd(iGRU)/=integerMissing) aJac(ixAqWat(iGRU),ixVegHyd(iGRU),iGRU) = -dAquiferTrans_dCanWat(iGRU)*dt  ! dVol/dLiq (kg m-2)-1
              if(ixTopNrg(iGRU)/=integerMissing) aJac(ixAqWat(iGRU),ixTopNrg(iGRU),iGRU) = -dAquiferTrans_dTGround(iGRU)*dt ! dVol/dT (K-1)
            endif
          endif
        end do
  
        ! -----
        ! * cross derivatives in the soil domain...
        ! ----------------------------------------
        if(nSoilOnlyHyd>0 .and. nSoilOnlyNrg>0)then
          !$cuf kernel do(2) <<<*,*>>>
          do iGRU=1,nGRU
          do iLayer=1,nSoilOnlyNrg

            ! - check that the soil layer is desired
            if(ixSoilOnlyNrg_m(iLayer,iGRU)==integerMissing) cycle

            ! - define indices of the soil layers
            jLayer   = iLayer+nSnow(iGRU)                  ! index of layer in the snow+soil vector

            ! - define the energy state variable
            nrgState = ixNrgLayer(jLayer,iGRU)       ! index within the full state vector

            ! - define index of hydrology state variable within the state subset
            watState = ixSoilOnlyHyd_m(iLayer,iGRU)

            ! only compute derivatives if the water state for the current layer is within the state subset
            if(watstate/=integerMissing)then

              ! - include derivatives in liquid water fluxes w.r.t. temperature for current layer
              aJac(watState,nrgState,iGRU) = (dt/mLayerDepth(jLayer,iGRU))*(-dq_dNrgStateBelow(iLayer-1,iGRU) + dq_dNrgStateAbove(iLayer,iGRU))   ! dVol/dT (K-1) -- flux depends on ice impedance

              ! - compute lower diagonal elements
              if(iLayer>1)then
                if(ixSoilOnlyHyd_m(iLayer-1,iGRU)/=integerMissing) aJac(ixSoilOnlyHyd_m(iLayer-1,iGRU),nrgState,iGRU) = (dt/mLayerDepth(jLayer-1,iGRU))*( dq_dNrgStateBelow(iLayer-1,iGRU))   ! K-1
              endif

              ! compute upper-diagonal elements
              if(iLayer<nSoil)then
                if(ixSoilOnlyHyd_m(iLayer+1,iGRU)/=integerMissing) aJac(ixSoilOnlyHyd_m(iLayer+1,iGRU),nrgState,iGRU) = (dt/mLayerDepth(jLayer+1,iGRU))*(-dq_dNrgStateAbove(iLayer,iGRU))     ! K-1
              endif

              ! - include derivatives of energy w.r.t. ground evaporation
              if(nSnow(iGRU)==0 .and. iLayer==1)then  ! upper-most soil layer
                if(computeVegFlux)then
                  if(ixCasNrg(iGRU)/=integerMissing) aJac(ixTopHyd(iGRU),ixCasNrg(iGRU),iGRU) = (dt/mLayerDepth(jLayer,iGRU))*(-dGroundEvaporation_dTCanair(iGRU)/iden_water) ! dVol/dT (K-1)
                  if(ixVegNrg(iGRU)/=integerMissing) aJac(ixTopHyd(iGRU),ixVegNrg(iGRU),iGRU) = (dt/mLayerDepth(jLayer,iGRU))*(-dGroundEvaporation_dTCanopy(iGRU)/iden_water) + aJac(ixTopHyd(iGRU),ixVegNrg(iGRU),iGRU) ! dVol/dT (K-1)
                  if(ixVegHyd(iGRU)/=integerMissing) aJac(ixTopHyd(iGRU),ixVegHyd(iGRU),iGRU) = (dt/mLayerDepth(jLayer,iGRU))*(-dGroundEvaporation_dCanWat(iGRU)/iden_water)  + aJac(ixTopHyd(iGRU),ixVegHyd(iGRU),iGRU) ! dVol/dLiq (kg m-2)-1
                endif
                if(ixTopNrg(iGRU)/=integerMissing) aJac(ixTopHyd(iGRU),ixTopNrg(iGRU),iGRU) = (dt/mLayerDepth(jLayer,iGRU))*(-dGroundEvaporation_dTGround(iGRU)/iden_water) + aJac(ixTopHyd(iGRU),ixTopNrg(iGRU),iGRU) ! dVol/dT (K-1)
              endif

              ! - include derivatives of energy and water w.r.t soil transpiration (dependent on canopy transpiration)
              if(computeVegFlux)then
                if(ixCasNrg(iGRU)/=integerMissing) aJac(watState,ixCasNrg(iGRU),iGRU) = (dt/mLayerDepth(jLayer,iGRU))*(-mLayerdTrans_dTCanair(iLayer,iGRU)) + aJac(watState,ixCasNrg(iGRU),iGRU) ! dVol/dT (K-1)
                if(ixVegNrg(iGRU)/=integerMissing) aJac(watState,ixVegNrg(iGRU),iGRU) = (dt/mLayerDepth(jLayer,iGRU))*(-mLayerdTrans_dTCanopy(iLayer,iGRU)) + aJac(watState,ixVegNrg(iGRU),iGRU) ! dVol/dT (K-1)
                if(ixVegHyd(iGRU)/=integerMissing) aJac(watState,ixVegHyd(iGRU),iGRU) = (dt/mLayerDepth(jLayer,iGRU))*(-mLayerdTrans_dCanWat(iLayer,iGRU))  + aJac(watState,ixVegHyd(iGRU),iGRU) ! dVol/dLiq (kg m-2)-1
                if(ixTopNrg(iGRU)/=integerMissing) aJac(watState,ixTopNrg(iGRU),iGRU) = (dt/mLayerDepth(jLayer,iGRU))*(-mLayerdTrans_dTGround(iLayer,iGRU)) + aJac(watState,ixTopNrg(iGRU),iGRU) ! dVol/dT (K-1)
              endif

              ! - include derivatives in energy fluxes w.r.t. with respect to water for current layer
              aJac(nrgState,watState,iGRU) = dVolHtCapBulk_dPsi0(iLayer,iGRU) * mLayerTempPrime(jLayer,iGRU) &
                                       + mLayerCm(jLayer,iGRU) * dVolTot_dPsi0(iLayer,iGRU) * cj + dCm_dPsi0(iLayer,iGRU) * mLayerVolFracWatPrime(jLayer,iGRU) &
                                       + (dt/mLayerDepth(jLayer,iGRU))*(-dNrgFlux_dWatBelow(jLayer-1,iGRU) + dNrgFlux_dWatAbove(jLayer,iGRU)) + mLayerCm(jLayer,iGRU) * d2VolTot_dPsi02(iLayer,iGRU) * mLayerMatricHeadPrime(iLayer,iGRU)
              if(mLayerdTheta_dTk(jLayer,iGRU) > verySmall)then  ! ice is present
                aJac(nrgState,watState,iGRU) = -LH_fu0*iden_water * dVolTot_dPsi0(iLayer,iGRU) * cj &
                                         - LH_fu0*iden_water * mLayerMatricHeadPrime(iLayer,iGRU) * d2VolTot_dPsi02(iLayer,iGRU) + aJac(nrgState,watState,iGRU) ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
              endif

              ! - include derivatives of heat capacity w.r.t water fluxes for surrounding layers starting with layer above
              if(iLayer>1)then
                if(ixSoilOnlyNrg_m(iLayer-1,iGRU)/=integerMissing) aJac(ixSoilOnlyNrg_m(iLayer-1,iGRU),watState,iGRU) = (dt/mLayerDepth(jLayer-1,iGRU))*( dNrgFlux_dWatBelow(jLayer-1,iGRU) )
              elseif(iLayer==1 .and. nSnow(iGRU)>0)then !top soil layer and there is snow above
                if(ixSnowOnlyNrg_m(nSnow(iGRU),iGRU)/=integerMissing) aJac(ixSnowOnlyNrg_m(nSnow(iGRU),iGRU),watState,iGRU) = (dt/mLayerDepth(nSnow(iGRU),iGRU))*( dNrgFlux_dWatBelow(nSnow(iGRU),iGRU) )
              endif

              ! (cross-derivative terms for the layer below)
              if(iLayer<nSoil)then
                if(ixSoilOnlyHyd_m(iLayer+1,iGRU)/=integerMissing) aJac(ixSoilOnlyNrg_m(iLayer+1,iGRU),watState,iGRU) = (dt/mLayerDepth(jLayer+1,iGRU))*(-dNrgFlux_dWatAbove(jLayer,iGRU  ) )
              endif

            endif   ! (if the water state for the current layer is within the state subset)

            ! - include terms for surface infiltration below surface
            ! if(ixSoilOnlyHyd_m(1,iGRU)/=integerMissing) aJac(ixSoilOnlyHyd_m(1,iGRU),nrgState,iGRU) = -(dt/mLayerDepth(1+nSnow,iGRU))*dq_dNrgStateLayerSurfVec(iLayer,iGRU) + aJac(ixSoilOnlyHyd_m(1,iGRU),nrgState,iGRU)

          end do  ! (looping through soil layers)
        end do

                  !$cuf kernel do(2) <<<*,*>>>
        do iGRU=1,nGRU
          do iLayer=1,nSoilOnlyNrg
            ! - check that the soil layer is desired
            if(ixSoilOnlyNrg_m(iLayer,iGRU)==integerMissing) cycle

            ! - define indices of the soil layers
            jLayer   = iLayer+nSnow(iGRU)                  ! index of layer in the snow+soil vector

            ! - define the energy state variable
            nrgState = ixNrgLayer(jLayer,iGRU)       ! index within the full state vector

            ! - include terms for surface infiltration below surface
            if(ixSoilOnlyHyd_m(1,iGRU)/=integerMissing) aJac(ixSoilOnlyHyd_m(1,iGRU),nrgState,iGRU) = -(dt/mLayerDepth(1+nSnow(iGRU),iGRU))*dq_dNrgStateLayerSurfVec(iLayer,iGRU) + aJac(ixSoilOnlyHyd_m(1,iGRU),nrgState,iGRU)
          end do  ! (looping through soil layers)
        end do


        ! err = cudaDeviceSynchronize()
        ! print*, 963, err
          ! - include terms for surface infiltration above surface
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
          if(nSnow(iGRU)>0)then !have snow above first soil layer
            if(ixSnowOnlyNrg_m(nSnow(iGRU),iGRU)/=integerMissing .and. ixSoilOnlyHyd_m(1,iGRU)/=integerMissing) aJac(ixSoilOnlyHyd_m(1,iGRU),ixSnowOnlyNrg_m(nSnow(iGRU),iGRU),iGRU) = -(dt/mLayerDepth(1+nSnow(iGRU),iGRU))*dq_dNrgStateLayerSurfVec(0,iGRU)
          elseif(computeVegFlux)then !have vegetation above first soil layer, ixTopHyd = ixSoilOnlyHyd_m(1,iGRU)
            if(ixVegNrg(iGRU)/=integerMissing .and. ixTopHyd(iGRU)/=integerMissing) aJac(ixTopHyd(iGRU),ixVegNrg(iGRU),iGRU) = -(dt/mLayerDepth(1+nSnow(iGRU),iGRU))*dq_dNrgStateLayerSurfVec(0,iGRU) + aJac(ixTopHyd(iGRU),ixVegNrg(iGRU),iGRU)
          endif
        end do
        endif   ! (if there are state variables for both water and energy in the soil domain)

    ! *********************************************************************************************************************************************************
    ! -----
    ! * if desired, modify to use enthalpy as a state variable instead of temperature 
    ! NOTE, dMat has been set to 0 and now 1._rkind * cj is added instead 
    ! ----------------------------------------
        if(enthalpyStateVec)then 

          ! allocate(watRows(nBands),nrgRows(nBands))
          ! do jLayer=1,nBands-1
            ! watRows(jLayer) = jLayer
            ! nrgRows(jLayer) = jLayer + 1
          ! end do
          ! watRows(nBands) = nBands
          ! nrgRows(nBands) = nBands
    
          !$cuf kernel do(2) <<<*,*>>>
          do iGRU=1,nGRU
          do iLayer=1,nState
            if(ixCasNrg(iGRU)/=integerMissing) aJac(iLayer,ixCasNrg(iGRU),iGRU) = aJac(iLayer,ixCasNrg(iGRU),iGRU) * dCanairTemp_dEnthalpy(iGRU)
            if (ixVegNrg(iGRU)/=integerMissing .and. ixVegHyd(iGRU)/=integerMissing) aJac(iLayer,ixVegHyd(iGRU),iGRU) = aJac(iLayer,ixVegHyd(iGRU),iGRU) + aJac(iLayer,ixVegNrg(iGRU),iGRU) * dCanopyTemp_dCanWat(iGRU)
            if (ixVegNrg(iGRU)/=integerMissing) aJac(iLayer,ixVegNrg(iGRU),iGRU) = aJac(iLayer,ixVegNrg(iGRU),iGRU) * dCanopyTemp_dEnthalpy(iGRU)
          end do
        end do
    
        !$cuf kernel do(1) <<<*,*>>>
        do iGRU=1,nGRU
          if(ixCasNrg(iGRU)/=integerMissing)then
            aJac(ixCasNrg(iGRU), ixCasNrg(iGRU),iGRU) = aJac(ixCasNrg(iGRU), ixCasNrg(iGRU),iGRU) + 1._rkind * cj
          endif
          
          if(ixVegNrg(iGRU)/=integerMissing)then
            aJac(ixVegNrg(iGRU), ixVegNrg(iGRU),iGRU) = aJac(ixVegNrg(iGRU), ixVegNrg(iGRU),iGRU) + 1._rkind * cj
          endif
        end do
          
          ! if(nSnowSoilNrg>0)then
            !$cuf kernel do(1) <<<*,*>>>
            do iGRU=1,nGRU
            do iLayer=1,nLayers(iGRU)
              do jLayer=1,nState
                nrgState = ixSnowSoilNrg_m(iLayer,iGRU)       
                ! if(nrgState==integerMissing) cycle
                watState = ixSnowSoilHyd_m(iLayer,iGRU)
                if(nrgState/=integerMissing .and. watstate/=integerMissing)then 
                  if(iLayer<=nSnow(iGRU)) then
                      aJac(jlayer,watState,iGRU) = aJac(jLayer,watState,iGRU) + aJac(jLayer,nrgState,iGRU) * dTemp_dTheta(iLayer,iGRU)
                  end if
                  if(iLayer>nSnow(iGRU))  then
                      aJac(jLayer,watState,iGRU) = aJac(jLayer,watState,iGRU) + aJac(jLayer,nrgState,iGRU) * dTemp_dPsi0(iLayer-nSnow(iGRU),iGRU)
                  end if
              endif
                if (nrgState/=integerMissing) aJac(jLayer,nrgState,iGRU) = aJac(jLayer,nrgState,iGRU) * dTemp_dEnthalpy(iLayer,iGRU)
            end do
          end do
        end do
          !$cuf kernel do(1) <<<*,*>>>
          do iGRU=1,nGRU
            do iLayer=1,nLayers(iGRU)
              nrgState = ixSnowSoilNrg_m(iLayer,iGRU)       
              if(nrgState==integerMissing) cycle
              aJac(nrgState, nrgState,iGRU) = aJac(nrgState, nrgState,iGRU) + 1._rkind * cj
            enddo
          end do
          endif
        ! endif
    
    ! if(any(ieee_is_nan(aJac)))then
      ! message=trim(message)//'NaN in Jacobian'
      ! err=20; return
    ! endif

  ! end association to variables in the data structures
  end associate

  associate(nState => indx_data%nState)
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
      do iLayer=nState(iGRU)+1,size(aJac,1)
        aJac(iLayer,iLayer,iGRU) = 1._rkind
      end do
    end do
  end associate

end subroutine computJacobWithPrime

! **********************************************************************************************************
! public function computJacob4ida: the interface to compute the Jacobian matrix dF/dy + c dF/dy' for IDA solver
! **********************************************************************************************************
! Return values:
!    0 = success,
!    1 = recoverable error,
!   -1 = non-recoverable error
! ----------------------------------------------------------------
integer(c_int) function computJacob4ida(t, cj, sunvec_y, sunvec_yp, sunvec_r, &
                    sunmat_J, user_data, sunvec_temp1, sunvec_temp2, sunvec_temp3) &
                    result(ierr) bind(C,name='computJacob4ida')

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fnvector_serial_mod
  use fsunmatrix_band_mod
  use fsunmatrix_dense_mod
  use fsunmatrix_magmadense_mod
  use type4ida

  !======= Declarations =========
  implicit none

  ! calling variables
  real(rkind), value            :: t              ! current time
  real(rkind), value            :: cj             ! step size scaling factor
  type(N_Vector)                :: sunvec_y       ! solution N_Vector
  type(N_Vector)                :: sunvec_yp      ! derivative N_Vector
  type(N_Vector)                :: sunvec_r       ! residual N_Vector
  type(SUNMatrix)               :: sunmat_J       ! Jacobian SUNMatrix
  type(c_ptr), value            :: user_data      ! user-defined data
  type(N_Vector)                :: sunvec_temp1   ! temporary N_Vector
  type(N_Vector)                :: sunvec_temp2   ! temporary N_Vector
  type(N_Vector)                :: sunvec_temp3   ! temporary N_Vector

  ! pointers to data in SUNDIALS vectors
  real(rkind), pointer,device          :: Jac(:,:,:)       ! Jacobian matrix
  type(data4ida), pointer       :: eqns_data      ! equations data
  ! ----------------------------------------------------------------

  ! get equations data from user-defined data
  call c_f_pointer(user_data, eqns_data)

  ! get data arrays from SUNDIALS vectors
   Jac => FSUNMatrix_MagmaDense_Data(sunmat_J)

  ! compute the analytical Jacobian matrix
  ! NOTE: The derivatives were computed in the previous call to computFlux
  !       This occurred either at the call to eval8summaWithPrime at the start of systemSolv
  !        or in the call to eval8summaWithPrime in the previous iteration
  call computJacobWithPrime(&
                ! input: model control
                cj,                                       & ! intent(in):    this scalar changes whenever the step size or method order changes
                1._qp,                                    & ! intent(in):    length of the time step (seconds)
                eqns_data%indx_data%nSnow,                          & ! intent(in):    number of snow layers
                eqns_data%nSoil,                          & ! intent(in):    number of soil layers
                eqns_data%indx_data%nLayers_d,                        & ! intent(in):    total number of layers
                eqns_data%nGRU, &
                eqns_data%computeVegFlux,                 & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                eqns_data%model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel, & ! intent(in): flag to indicate if we need to compute baseflow
                eqns_data%ixMatrix,                                                            & ! intent(in): form of the Jacobian matrix
                eqns_data%mpar_data%specificStorage,                    & ! intent(in): specific storage coefficient (m-1)
                eqns_data%mpar_data%theta_sat,                             & ! intent(in): soil porosity (-)
                eqns_data%decisions%f_Richards,                & ! intent(in): choice of option for Richards' equation
                eqns_data%model_decisions(iLookDECISIONS%nrgConserv)%iDecision.ne.closedForm,  & ! intent(in): flag if enthalpy is state variable
                ! input: data structures
                eqns_data%indx_data,                      & ! intent(in):    index data
                eqns_data%prog_data,                      & ! intent(in):    model prognostic variables for a local HRU
                eqns_data%diag_data,                      & ! intent(in):    model diagnostic variables for a local HRU
                eqns_data%deriv_data,                     & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                eqns_data%dBaseflow_dMatric,              & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
                ! input: state variables
                eqns_data%mLayerTempPrime,                & ! intent(in):    derivative value for temperature of each snow and soil layer (K)
                eqns_data%mLayerMatricHeadPrime,          & ! intent(in):    derivative value for matric head of each snow and soil layer (m)
                eqns_data%mLayerVolFracWatPrime,          & ! intent(in):    derivative value for volumetric total water content of each snow and soil layer (-)
                eqns_data%scalarCanopyTempPrime,          & ! intent(in):    derivative value for temperature of the vegetation canopy (K)
                eqns_data%scalarCanopyWatPrime,           & ! intent(in):    derivative value for total water content of the vegetation canopy (kg m-2)
                ! input-output: Jacobian and its diagonal
                eqns_data%dMat,                           & ! intent(inout): diagonal of the Jacobian matrix
                Jac,                                      & ! intent(out):   Jacobian matrix
                ! output: error control
                eqns_data%err,eqns_data%message)            ! intent(out):   error code and error message
  if(eqns_data%err > 0)then; eqns_data%message=trim(eqns_data%message); ierr=-1; return; endif
  if(eqns_data%err < 0)then; eqns_data%message=trim(eqns_data%message); ierr=1; return; endif

  ! return success
  ierr = 0
  return

end function computJacob4ida

! **********************************************************************************************************
! private function: get the off-diagonal index in the band-diagonal matrix
! **********************************************************************************************************
function ixOffDiag(jState,iState)
  implicit none
  integer(i4b),intent(in)  :: jState    ! off-diagonal state
  integer(i4b),intent(in)  :: iState    ! diagonal state
  integer(i4b)             :: ixOffDiag ! off-diagonal index in the band-diagonal matrix

  ixOffDiag = ixDiag + jState - iState
end function ixOffDiag

end module computJacobWithPrime_module
