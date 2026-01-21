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
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    model_options   ! defines the model decisions

! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDIAG       ! named variables for structure elements
USE var_lookup,only:iLookDERIV      ! named variables for structure elements

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! named variables to describe the state variable type
USE globalData,only:iname_watLayer   ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:maxVolIceContent ! snow maximum volumetric ice content to store water (-)

! access named variables to describe the form and structure of the matrices used in the numerical solver
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
USE mDecisions_module,only:       &
 qbaseTopmodel,                   & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                       & ! a big bucket (lumped aquifer model)
 noExplicit                         ! no explicit groundwater parameterization

! look-up values for the form of Richards' equation
USE mDecisions_module,only:       &
 moisture,                        & ! moisture-based form of Richards' equation
 mixdform                           ! mixed form of Richards' equation

! look-up values for the choice of variable in energy equations (BE residual or IDA state variable)
USE mDecisions_module,only:       &
 closedForm,                      & ! use temperature with closed form heat capacity
 enthalpyFormLU,                  & ! use enthalpy with soil temperature-enthalpy lookup tables
 enthalpyForm                       ! use enthalpy with soil temperature-enthalpy analytical solution

implicit none
private
public::computJacobWithPrime
public::computJacob4ida
logical::fullMatrix
contains


! **********************************************************************************************************
! public subroutine computJacobWithPrime: compute the Jacobian matrix
! **********************************************************************************************************
subroutine computJacobWithPrime(&
  nGRU, &
                      ! input: model control
                      cj,                         & ! intent(in):    this scalar changes whenever the step size or method order changes
                      dt,                         & ! intent(in):    length of the time step (seconds)
                      nSnow,                      & ! intent(in):    number of snow layers
                      nSoil,                      & ! intent(in):    number of soil layers
                      nLayers,                    & ! intent(in):    total number of layers
                      computeVegFlux,             & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                      computeBaseflow,            & ! intent(in):    flag to indicate if we need to compute baseflow
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
                      dMat0,                      & ! intent(in):    diagonal of the Jacobian matrix excluding fluxes, not depending on the state vector
                      aJac,                       & ! intent(out):   Jacobian matrix
                      ! output: error control
                      err,message)                  ! intent(out):   error code and error message
  ! -----------------------------------------------------------------------------------------------------------------
  ! provide access to subroutines
  use computJacob_module,only:fluxJacAdd
  use ieee_arithmetic
  use cudafor
  use device_data_types
  ! -----------------------------------------------------------------------------------------------------------------
  implicit none
  ! input: model control
  integer(i4b),intent(in) :: nGRU
  real(rkind),intent(in)               :: cj                         ! this scalar changes whenever the step size or method order changes
  real(rkind),intent(in)               :: dt                         ! length of the time step (seconds)
  integer(i4b),intent(in),device              :: nSnow(:)                      ! number of snow layers
  integer(i4b),intent(in)              :: nSoil                      ! number of soil layers
  integer(i4b),intent(in),device              :: nLayers(:)                    ! total number of layers in the snow and soil domains
  logical(lgt),intent(in),device              :: computeVegFlux(:)             ! flag to indicate if computing fluxes over vegetation
  logical(lgt),intent(in)              :: computeBaseflow            ! flag to indicate if computing baseflow
  real(rkind),intent(in),device               :: specificStorage(:)            ! specific storage coefficient (m-1)
  real(rkind),intent(in),device               :: theta_sat(:,:)               ! soil porosity (-)
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
  real(rkind),intent(in),device               :: dMat0(:,:)                   ! diagonal of the Jacobian matrix excluding fluxes, not depending on the state vector
  real(rkind),intent(out),device              :: aJac(:,:,:)                  ! Jacobian matrix
  ! output variables
  integer(i4b),intent(out)             :: err                        ! error code
  character(*),intent(out)             :: message                    ! error message
  ! --------------------------------------------------------------
  ! * local variables
  ! --------------------------------------------------------------
  real(rkind),allocatable,device              :: dMat(:,:)                    ! diagonal of the Jacobian matrix excluding fluxes, depending on the state vector
  ! indices of model state variables
  integer(i4b)                         :: nrgState                   ! energy state variable
  integer(i4b)                         :: watState                   ! hydrology state variable
  integer(i4b)                         :: nState                     ! number of state variables
  ! integer(i4b),allocatable             :: nrgRows(:)                 ! indices of rows for energy column
  ! integer(i4b),allocatable             :: watRows(:)                 ! indices of rows for hydrology column
  ! indices of model layers
  integer(i4b)                         :: iLayer                     ! index of model layer
  integer(i4b)                         :: jLayer                     ! index of model layer within the full state vector (hydrology)
  ! conversion factors
  real(rkind)                          :: LH_fu0                     ! latent heat of fusion, modified to be 0 if using enthalpy formulation and not using
  character(LEN=256)                   :: cmessage                   ! error message of downwind routine

  integer(i4b) :: iGRU
        type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  ! --------------------------------------------------------------
  ! associate variables from data structures
  associate(&
    ! indices of model state variables
    ixCasNrg                     => indx_data%ixCasNrg                  ,& ! intent(in): [i4b]    index of canopy air space energy state variable
    ixVegNrg                     => indx_data%ixVegNrg                  ,& ! intent(in): [i4b]    index of canopy energy state variable
    ixVegHyd                     => indx_data%ixVegHyd                  ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixAqWat                      => indx_data%ixAqWat                   ,& ! intent(in): [i4b]    index of water storage in the aquifer
    ! vector of energy indices for the snow and soil domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnowSoilNrg                => indx_data%ixSnowSoilNrg                ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow and soil domains
    ixSnowOnlyNrg                => indx_data%ixSnowOnlyNrg                ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow domain
    ixSoilOnlyNrg                => indx_data%ixSoilOnlyNrg                ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the soil domain
    ! vector of hydrology indices for the snow and soil domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnowSoilHyd                => indx_data%ixSnowSoilHyd                ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow and soil domains
    ixSnowOnlyHyd                => indx_data%ixSnowOnlyHyd                ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow domain
    ixSoilOnlyHyd                => indx_data%ixSoilOnlyHyd                ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the soil domain
    ! number of state variables of a specific type
    ! nSnowSoilNrg                 => indx_data%var(iLookINDEX%nSnowSoilNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the snow and soil domains
    ! nSnowOnlyNrg                 => indx_data%var(iLookINDEX%nSnowOnlyNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the snow domain
    ! nSoilOnlyNrg                 => indx_data%var(iLookINDEX%nSoilOnlyNrg)%dat(1)              ,& ! intent(in): [i4b]    number of energy state variables in the soil domain
    ! nSnowOnlyHyd                 => indx_data%var(iLookINDEX%nSnowOnlyHyd)%dat(1)              ,& ! intent(in): [i4b]    number of hydrology variables in the snow domain
    ! nSoilOnlyHyd                 => indx_data%var(iLookINDEX%nSoilOnlyHyd)%dat(1)              ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
    ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
    dCanopyNetFlux_dCanWat       => deriv_data%dCanopyNetFlux_dCanWat   ,& ! intent(in): [dp]     derivative in net canopy fluxes w.r.t. canopy total water content
    ! derivatives in canopy water w.r.t canopy temperature
    dTheta_dTkCanopy             => deriv_data%dTheta_dTkCanopy         ,& ! intent(in): [dp]     derivative in volumetric liquid water content w.r.t. temperature
    d2Theta_dTkCanopy2           => deriv_data%d2Theta_dTkCanopy2       ,& ! intent(in): [dp]     second derivative in volumetric liquid water content w.r.t. temperature
    dFracLiqVeg_dTkCanopy        => deriv_data%dFracLiqVeg_dTkCanopy    ,& ! intent(in): [dp]     derivative in fraction of (throughfall + drainage)  w.r.t. temperature
    ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. water state in layers above and below
    dNrgFlux_dWatAbove           => deriv_data%dNrgFlux_dWatAbove_m          ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. water state in the layer above
    dNrgFlux_dWatBelow           => deriv_data%dNrgFlux_dWatBelow_m          ,& ! intent(in): [dp(:)]  derivatives in the flux w.r.t. water state in the layer below
    ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
    dVolTot_dPsi0                => deriv_data%dVolTot_dPsi0_m               ,& ! intent(in): [dp(:)]  derivatives in total water content w.r.t. total water matric potential
    d2VolTot_dPsi02              => deriv_data%d2VolTot_dPsi02_m             ,& ! intent(in): [dp(:)]  second derivatives in total water content w.r.t. total water matric potential
    dCompress_dPsi               => deriv_data%dCompress_dPsi_m              ,& ! intent(in): [dp(:)]  derivatives in compressibility w.r.t matric head
    ! derivative in liquid water fluxes for the soil and snow domain w.r.t temperature
    dFracLiqWat_dTk              => deriv_data%dFracLiqWat_dTk_m             ,& ! intent(in): [dp(:)]  derivatives in fraction of liquid water w.r.t. temperature
    mLayerdTheta_dTk             => deriv_data%mLayerdTheta_dTk_m            ,& ! intent(in): [dp(:)]  derivatives in volumetric liquid water content w.r.t. temperature
    mLayerd2Theta_dTk2           => deriv_data%mLayerd2Theta_dTk2_m          ,& ! intent(in): [dp(:)]  second derivatives of volumetric liquid water content w.r.t. temperature
    ! derivative in bulk heat capacity w.r.t. relevant state variables
    dVolHtCapBulk_dPsi0          => deriv_data%dVolHtCapBulk_dPsi0_m         ,& ! intent(in): [dp(:)]  derivatives in bulk heat capacity w.r.t. matric potential
    dVolHtCapBulk_dTheta         => deriv_data%dVolHtCapBulk_dTheta_m        ,& ! intent(in): [dp(:)]  derivatives in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dCanWat        => deriv_data%dVolHtCapBulk_dCanWat    ,& ! intent(in): [dp   ]  derivative in bulk heat capacity w.r.t. volumetric water content
    dVolHtCapBulk_dTk            => deriv_data%dVolHtCapBulk_dTk_m           ,& ! intent(in): [dp(:)]  derivatives in bulk heat capacity w.r.t. temperature
    dVolHtCapBulk_dTkCanopy      => deriv_data%dVolHtCapBulk_dTkCanopy  ,& ! intent(in): [dp   ]  derivative in bulk heat capacity w.r.t. temperature
    ! derivative in Cm w.r.t. relevant state variables
    dCm_dPsi0                    => deriv_data%dCm_dPsi0_m                   ,& ! intent(in): [dp(:)]  derivative sin heat capacity w.r.t. matric potential (J kg-1)
    dCm_dTk                      => deriv_data%dCm_dTk_m                     ,& ! intent(in): [dp(:)]  derivatives in heat capacity w.r.t. temperature (J kg-1 K-2)
    dCm_dTkCanopy                => deriv_data%dCm_dTkCanopy            ,& ! intent(in): [dp   ]  derivative in heat capacity w.r.t. canopy temperature (J kg-1 K-2)
    ! derivatives of temperature if enthalpy is the state variable
    dCanairTemp_dEnthalpy        => deriv_data%dCanairTemp_dEnthalpy    ,& ! intent(in): [dp]     derivative incanopy air temperature w.r.t. enthalpy
    dCanopyTemp_dEnthalpy        => deriv_data%dCanopyTemp_dEnthalpy    ,& ! intent(in): [dp]     derivative incanopy temperature w.r.t. enthalpy 
    dTemp_dEnthalpy              => deriv_data%dTemp_dEnthalpy_m             ,& ! intent(in): [dp(:)]  derivatives in temperature w.r.t. enthalpy
    dCanopyTemp_dCanWat          => deriv_data%dCanopyTemp_dCanWat      ,& ! intent(in): [dp]     derivative incanopy temperature w.r.t. volumetric water content
    dTemp_dTheta                 => deriv_data%dTemp_dTheta_m                ,& ! intent(in): [dp(:)]  derivatives in temperature w.r.t. volumetric water content
    dTemp_dPsi0                  => deriv_data%dTemp_dPsi0_m                 ,& ! intent(in): [dp(:)]  derivatives in temperature w.r.t. total water matric potential
    ! diagnostic variables
    scalarFracLiqVeg             => diag_data%scalarFracLiqVeg           ,& ! intent(in): [dp]     fraction of liquid water on vegetation (-)
    scalarBulkVolHeatCapVeg      => diag_data%scalarBulkVolHeatCapVeg    ,& ! intent(in): [dp]     bulk volumetric heat capacity of vegetation (J m-3 K-1)
    scalarCanopyCm               => diag_data%scalarCanopyCm             ,& ! intent(in): [dp]     Cm of canopy (J kg-1 K-1)
    mLayerFracLiqSnow            => diag_data%mLayerFracLiqSnow          ,& ! intent(in): [dp(:)]  fraction of liquid water in each snow layer (-)
    mLayerVolHtCapBulk           => diag_data%mLayerVolHtCapBulk         ,& ! intent(in): [dp(:)]  bulk volumetric heat capacity in each snow+soil layer (J m-3 K-1)
    mLayerCm                     => diag_data%mLayerCm                   ,& ! intent(in): [dp(:)]  Cm in each snow+soil layer (J kg-1 K-1)
    ! canopy and layer depth
    canopyDepth                  => diag_data%scalarCanopyDepth          ,& ! intent(in): [dp   ]  canopy depth (m)
    mLayerDepth                  => prog_data%mLayerDepth                    & ! intent(in): [dp(:)]  depth of each layer in the snow+soil sub-domain (m)
    ) ! making association with data in structures
    ! --------------------------------------------------------------
    ! initialize error control
    err=0; message='computJacobWithPrime/'

    ! *********************************************************************************************************************************************************
    ! * PART 0: PRELIMINARIES (INITIALIZE JACOBIAN AND COMPUTE TIME-VARIABLE DIAGONAL TERMS)
    ! *********************************************************************************************************************************************************
    ! get the number of state variables
    nState = size(dMat0,1)
    allocate(dMat(nState,nGRU)) 

    ! call allocate_device_indx_data(indx_data_d,indx_data,indx_data%var(iLookINDEX%nSnow)%dat(1),nSoil,indx_data%var(iLookINDEX%nLayers)%dat(1),nGRU)
    ! call allocate_device_deriv_data(deriv_data_d,deriv_data,nGRU,nSoil,indx_data%var(iLookINDEX%nLayers)%dat(1),indx_data%var(iLookINDEX%nSnow)%dat(1))
    ! initialize the Jacobian and diagonal
    ! NOTE: this needs to be done every time, since Jacobian matrix is modified in the solver and dMat is modified below
    ! aJac(:,:) = 0._rkind  ! analytical Jacobian matrix
    ! dMat = dMat0 * cj ! dMat0(ixCasNrg) = Cp_air*iden_air and dMat0(Wat states) = 1.0

    call computJacobWithPrime_kernel<<<blocks,threads>>>(nGRU,aJac, dMat, dMat0, cj, &
  computeVegFlux, enthalpyStateVec, ixRichards, &
  nSnow, nSoil, nLayers, indx_data%nState,&
  ixCasNrg, ixVegNrg, ixVegHyd, ixAqWat, &
  ixSnowSoilNrg, ixSnowOnlyNrg, ixSoilOnlyNrg, &
  ixSnowSoilHyd, ixSnowOnlyHyd, ixSoilOnlyHyd,&
  diag_data%scalarBulkVolHeatCapVeg, dTheta_dTkCanopy, dVolHtCapBulk_dTkCanopy,&
  scalarCanopyTempPrime, dCm_dTkCanopy, scalarCanopyWatPrime,&
  diag_data%scalarCanopyDepth, d2Theta_dTkCanopy2, dFracLiqVeg_dTkCanopy, &
  diag_data%mLayerVolHtCapBulk, mLayerdTheta_dTk,&
  dVolHtCapBulk_dTk, mLayerTempPrime, &
  dCm_dTk, dFracLiqWat_dTk, &
  mLayerVolFracWatPrime, mLayerd2Theta_dTk2,&
  dVolTot_dPsi0, d2VolTot_dPsi02, dCompress_dPsi, &
  mLayerMatricHeadPrime, theta_sat, &
  specificStorage, &
  diag_data%scalarFracLiqVeg, dVolHtCapBulk_dCanWat, diag_data%scalarCanopyCm, &
  dt, dCanopyNetFlux_dCanWat,&
  diag_data%mLayerFracLiqSnow, dVolHtCapBulk_dTheta, diag_data%mLayerCm, &
  prog_data%mLayerDepth, &
  dNrgFlux_dWatAbove, dNrgFlux_dWatBelow,&
  dCm_dPsi0, dVolHtCapBulk_dPsi0,&
  computeBaseflow, &
  indx_data%ixTopNrg, indx_data%ixTopHyd, indx_data%ixHydType,&
  prog_data%mLayerVolFracIce, diag_data%scalarSoilControl,&
  deriv_data%dCanairNetFlux_dCanairTemp, deriv_data%dCanairNetFlux_dCanopyTemp, deriv_data%dCanairNetFlux_dGroundTemp, &
  deriv_data%dCanopyNetFlux_dCanairTemp, deriv_data%dCanopyNetFlux_dCanopyTemp, deriv_data%dCanopyNetFlux_dGroundTemp, &
  deriv_data%dGroundNetFlux_dCanairTemp, deriv_data%dGroundNetFlux_dCanopyTemp, deriv_data%dGroundNetFlux_dCanWat,&
  deriv_data%dCanopyEvaporation_dTCanair,deriv_data%dCanopyEvaporation_dTCanopy,deriv_data%dCanopyEvaporation_dTGround,deriv_data%dCanopyEvaporation_dCanWat,&
  deriv_data%dGroundEvaporation_dTCanair,deriv_data%dGroundEvaporation_dTCanopy,deriv_data%dGroundEvaporation_dTGround,deriv_data%dGroundEvaporation_dCanWat,&
  deriv_data%scalarCanopyLiqDeriv, deriv_data%dCanLiq_dTcanopy,&
  deriv_data%dNrgFlux_dTempAbove_m, deriv_data%dNrgFlux_dTempBelow_m,&
  deriv_data%mLayerdTrans_dTCanair_m, deriv_data%mLayerdTrans_dTCanopy_m, deriv_data%mLayerdTrans_dTGround_m, deriv_data%mLayerdTrans_dCanWat_m, &
  deriv_data%dAquiferTrans_dTCanair, deriv_data%dAquiferTrans_dTCanopy, deriv_data%dAquiferTrans_dTGround, deriv_data%dAquiferTrans_dCanWat,&
  deriv_data%iLayerLiqFluxSnowDeriv_m,&
  deriv_data%dq_dHydStateAbove_m, deriv_data%dq_dHydStateBelow_m, deriv_data%dq_dHydStateLayerSurfVec_m, &
  deriv_data%dBaseflow_dAquifer, &
  deriv_data%dq_dNrgStateAbove_m, deriv_data%dq_dNrgStateBelow_m, deriv_data%dq_dNrgStateLayerSurfVec_m,&
  dBaseflow_dMatric,&
  dTemp_dTheta,dTemp_dPsi0,dTemp_dEnthalpy,&
  dCanairTemp_dEnthalpy, dCanopyTemp_dCanWat, dCanopyTemp_dEnthalpy)


    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    deallocate(dMat)

    ! !$cuf kernel do(1) <<<1,1>>>
    ! do iGRU=1,nGRU
    !   do iLayer=1,size(aJac,1)
    !     do jLayer=1,size(aJac,2)
    !       print*, iLayer,jLayer,iGRU,aJac(iLayer,jLayer,iGRU)
    !     end do
    !   end do
    ! end do
    
    ! print the Jacobian
    ! if(any(ieee_is_nan(aJac)))then; message=trim(message)//'NaN in Jacobian';err=20; return; endif

  end associate ! end association to variables in the data structures  

end subroutine computJacobWithPrime

attributes(global)subroutine computJacobWithPrime_kernel(nGRU,aJac, dMat, dMat0, cj, &
  computeVegFlux, enthalpyStateVec, ixRichards, &
  nSnow, nSoil, nLayers, nState,&
  ixCasNrg, ixVegNrg, ixVegHyd, ixAqWat, &
  ixSnowSoilNrg, ixSnowOnlyNrg, ixSoilOnlyNrg, &
  ixSnowSoilHyd, ixSnowOnlyHyd, ixSoilOnlyHyd,&
  scalarBulkVolHeatCapVeg, dTheta_dTkCanopy, dVolHtCapBulk_dTkCanopy,&
  scalarCanopyTempPrime, dCm_dTkCanopy, scalarCanopyWatPrime,&
  canopyDepth, d2Theta_dTkCanopy2, dFracLiqVeg_dTkCanopy, &
  mLayerVolHtCapBulk, mLayerdTheta_dTk,&
  dVolHtCapBulk_dTk, mLayerTempPrime, &
  dCm_dTk, dFracLiqWat_dTk, &
  mLayerVolFracWatPrime, mLayerd2Theta_dTk2, &
  dVolTot_dPsi0, d2VolTot_dPsi02, dCompress_dPsi, &
  mLayerMatricHeadPrime, theta_sat, &
  specificStorage, &
  scalarFracLiqVeg, dVolHtCapBulk_dCanWat, scalarCanopyCm, &
  dt, dCanopyNetFlux_dCanWat, &
  mLayerFracLiqSnow, dVolHtCapBulk_dTheta, mLayerCm, &
  mLayerDepth, &
  dNrgFlux_dWatAbove, dNrgFlux_dWatBelow,&
  dCm_dPsi0, dVolHtCapBulk_dPsi0,&
  computeBaseflow, &
  ixTopNrg, ixTopHyd, ixHydType,&
  mLayerVolFracIce, scalarSoilControl,&
  dCanairNetFlux_dCanairTemp, dCanairNetFlux_dCanopyTemp, dCanairNetFlux_dGroundTemp, &
  dCanopyNetFlux_dCanairTemp, dCanopyNetFlux_dCanopyTemp, dCanopyNetFlux_dGroundTemp, &
  dGroundNetFlux_dCanairTemp, dGroundNetFlux_dCanopyTemp, dGroundNetFlux_dCanWat,&
  dCanopyEvaporation_dTCanair,dCanopyEvaporation_dTCanopy,dCanopyEvaporation_dTGround,dCanopyEvaporation_dCanWat,&
  dGroundEvaporation_dTCanair,dGroundEvaporation_dTCanopy,dGroundEvaporation_dTGround,dGroundEvaporation_dCanWat,&
  scalarCanopyLiqDeriv, dCanLiq_dTcanopy,&
  dNrgFlux_dTempAbove, dNrgFlux_dTempBelow,&
  mLayerdTrans_dTCanair, mLayerdTrans_dTCanopy, mLayerdTrans_dTGround, mLayerdTrans_dCanWat, &
  dAquiferTrans_dTCanair, dAquiferTrans_dTCanopy, dAquiferTrans_dTGround, dAquiferTrans_dCanWat,&
  iLayerLiqFluxSnowDeriv,&
  dq_dHydStateAbove, dq_dHydStateBelow, dq_dHydStateLayerSurfVec, &
  dBaseflow_dAquifer, &
  dq_dNrgStateAbove, dq_dNrgStateBelow, dq_dNrgStateLayerSurfVec,&
  dBaseflow_dMatric,&
  dTemp_dTheta,dTemp_dPsi0,dTemp_dEnthalpy,&
  dCanairTemp_dEnthalpy, dCanopyTemp_dCanWat, dCanopyTemp_dEnthalpy)
  integer(i4b), value :: nGRU, nSoil
  real(rkind),intent(inout) :: aJac(:,:,:)
  real(rkind),intent(inout) :: dMat(:,:)
  real(rkind),intent(in) :: dMat0(:,:)
  real(rkind),intent(in),value :: cj
  logical(lgt),intent(in) :: computeVegFlux(:)
  logical(lgt),intent(in),value :: enthalpyStateVec
  integer(i4b),intent(in) :: ixRichards
  integer(i4b),intent(in) :: nSnow(:), nLayers(:),nState(:)
  integer(i4b),intent(in) :: ixCasNrg(:), ixVegNrg(:), ixVegHyd(:), ixAqWat(:)
  integer(i4b),intent(in) :: ixSnowSoilNrg(:,:), ixSnowOnlyNrg(:,:), ixSoilOnlyNrg(:,:)
  integer(i4b),intent(in) :: ixSnowSoilHyd(:,:), ixSnowOnlyHyd(:,:), ixSoilOnlyHyd(:,:)
  real(rkind),intent(in) :: scalarBulkVolHeatCapVeg(:), dTheta_dTkCanopy(:), dVolHtCapBulk_dTkCanopy(:)
  real(rkind),intent(in) :: scalarCanopyTempPrime(:), dCm_dTkCanopy(:), scalarCanopyWatPrime(:)
  real(rkind),intent(in) :: canopyDepth(:), d2Theta_dTkCanopy2(:), dFracLiqVeg_dTkCanopy(:)
  real(rkind),intent(in) :: mLayerVolHtCapBulk(:,:), mLayerdTheta_dTk(:,:)
  real(rkind),intent(in) :: dVolHtCapBulk_dTk(:,:), mLayerTempPrime(:,:)
  real(rkind),intent(in) :: dCm_dTk(:,:), dFracLiqWat_dTk(:,:)
  real(rkind),intent(in) :: mLayerVolFracWatPrime(:,:), mLayerd2Theta_dTk2(:,:)
  real(rkind),intent(in) :: dVolTot_dPsi0(:,:), d2VolTot_dPsi02(:,:), dCompress_dPsi(:,:)
  real(rkind),intent(in) :: mLayerMatricHeadPrime(:,:), theta_sat(:,:)
  real(rkind),intent(in) :: specificStorage(:)
  real(rkind),intent(in) :: scalarFracLiqVeg(:), dVolHtCapBulk_dCanWat(:), scalarCanopyCm(:)
  real(rkind),intent(in),value :: dt
  real(rkind),intent(in) :: dCanopyNetFlux_dCanWat(:)
  real(rkind),intent(in) :: mLayerFracLiqSnow(:,:), dVolHtCapBulk_dTheta(:,:), mLayerCm(:,:)
  real(rkind),intent(in) :: mLayerDepth(:,:)
  real(rkind),intent(in) :: dNrgFlux_dWatAbove(0:,:), dNrgFlux_dWatBelow(0:,:)
  real(rkind),intent(in) :: dCm_dPsi0(:,:), dVolHtCapBulk_dPsi0(:,:)
  logical(lgt),intent(in),value :: computeBaseflow
  integer(i4b),intent(in) :: ixTopNrg(:), ixTopHyd(:)
  integer(i4b),intent(in) :: ixHydType(:,:)
  real(rkind),intent(in) :: mLayerVolFracIce(:,:)
  real(rkind),intent(in) :: scalarSoilControl(:)
  real(rkind),intent(in) :: dCanairNetFlux_dCanairTemp(:), dCanairNetFlux_dCanopyTemp(:), dCanairNetFlux_dGroundTemp(:)
  real(rkind),intent(in) :: dCanopyNetFlux_dCanairTemp(:),dCanopyNetFlux_dCanopyTemp(:),dCanopyNetFlux_dGroundTemp(:)
  real(rkind),intent(in) :: dGroundNetFlux_dCanairTemp(:),dGroundNetFlux_dCanopyTemp(:),dGroundNetFlux_dCanWat(:)
  real(rkind),intent(in) :: dCanopyEvaporation_dTCanair(:),dCanopyEvaporation_dTCanopy(:),dCanopyEvaporation_dTGround(:),dCanopyEvaporation_dCanWat(:)
  real(rkind),intent(in) :: dGroundEvaporation_dTCanair(:),dGroundEvaporation_dTCanopy(:),dGroundEvaporation_dTGround(:),dGroundEvaporation_dCanWat(:)
  real(rkind),intent(in) :: scalarCanopyLiqDeriv(:), dCanLiq_dTcanopy(:)
  real(rkind),intent(in) :: dNrgFlux_dTempAbove(0:,:), dNrgFlux_dTempBelow(0:,:)
  real(rkind),intent(in) :: mLayerdTrans_dTCanair(:,:),mLayerdTrans_dTCanopy(:,:),mLayerdTrans_dTGround(:,:),mLayerdTrans_dCanWat(:,:)
  real(rkind),intent(in) :: dAquiferTrans_dTCanair(:),dAquiferTrans_dTCanopy(:),dAquiferTrans_dTGround(:),dAquiferTrans_dCanWat(:)
  real(rkind),intent(in) :: iLayerLiqFluxSnowDeriv(0:,:)
  real(rkind),intent(in) :: dq_dHydStateAbove(0:,:),dq_dHydStateBelow(0:,:),dq_dHydStateLayerSurfVec(0:,:)
  real(rkind),intent(in) :: dBaseflow_dAquifer(:)
  real(rkind),intent(in) :: dq_dNrgStateAbove(0:,:),dq_dNrgStateBelow(0:,:),dq_dNrgStateLayerSurfVec(0:,:)
  real(rkind),intent(in) :: dBaseflow_dMatric(:,:,:)
  real(rkind),intent(in) :: dTemp_dTheta(:,:), dTemp_dPsi0(:,:), dTemp_dEnthalpy(:,:)
  real(rkind),intent(in) :: dCanairTemp_dEnthalpy(:), dCanopyTemp_dCanWat(:), dCanopyTemp_dEnthalpy(:)

    integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

  call computJacobWithPrime_device(aJac, dMat(:,iGRU), dMat0(:,iGRU), cj, &
  computeVegFlux(iGRU), enthalpyStateVec, ixRichards, &
  nSnow(iGRU), nSoil, nLayers(iGRU), nState(iGRU),&
  ixCasNrg(iGRU), ixVegNrg(iGRU), ixVegHyd(iGRU), ixAqWat(iGRU), &
  ixSnowSoilNrg(:,iGRU), ixSnowOnlyNrg(:,iGRU), ixSoilOnlyNrg(:,iGRU), &
  ixSnowSoilHyd(:,iGRU), ixSnowOnlyHyd(:,iGRU), ixSoilOnlyHyd(:,iGRU),&
  scalarBulkVolHeatCapVeg(iGRU), dTheta_dTkCanopy(iGRU), dVolHtCapBulk_dTkCanopy(iGRU),&
  scalarCanopyTempPrime(iGRU), dCm_dTkCanopy(iGRU), scalarCanopyWatPrime(iGRU),&
  canopyDepth(iGRU), d2Theta_dTkCanopy2(iGRU), dFracLiqVeg_dTkCanopy(iGRU),&
  mLayerVolHtCapBulk(:,iGRU), mLayerdTheta_dTk(:,iGRU),&
  dVolHtCapBulk_dTk(:,iGRU), mLayerTempPrime(:,iGRU), &
  dCm_dTk(:,iGRU), dFracLiqWat_dTk(:,iGRU), &
  mLayerVolFracWatPrime(:,iGRU), mLayerd2Theta_dTk2(:,iGRU), &
  dVolTot_dPsi0(:,iGRU), d2VolTot_dPsi02(:,iGRU), dCompress_dPsi(:,iGRU), &
  mLayerMatricHeadPrime(:,iGRU), theta_sat(:,iGRU), &
  specificStorage(iGRU), &
  scalarFracLiqVeg(iGRU), dVolHtCapBulk_dCanWat(iGRU), scalarCanopyCm(iGRU), &
  dt, dCanopyNetFlux_dCanWat(iGRU),&
  mLayerFracLiqSnow(:,iGRU), dVolHtCapBulk_dTheta(:,iGRU), mLayerCm(:,iGRU), &
  mLayerDepth(:,iGRU), &
  dNrgFlux_dWatAbove(:,iGRU), dNrgFlux_dWatBelow(:,iGRU),&
  dCm_dPsi0(:,iGRU), dVolHtCapBulk_dPsi0(:,iGRU),&
  computeBaseflow, &
  ixTopNrg(iGRU), ixTopHyd(iGRU), ixHydType(:,iGRU),&
  mLayerVolFracIce(:,iGRU), scalarSoilControl(iGRU),&
  dCanairNetFlux_dCanairTemp(iGRU), dCanairNetFlux_dCanopyTemp(iGRU), dCanairNetFlux_dGroundTemp(iGRU), &
  dCanopyNetFlux_dCanairTemp(iGRU), dCanopyNetFlux_dCanopyTemp(iGRU), dCanopyNetFlux_dGroundTemp(iGRU), &
  dGroundNetFlux_dCanairTemp(iGRU), dGroundNetFlux_dCanopyTemp(iGRU), dGroundNetFlux_dCanWat(iGRU),&
  dCanopyEvaporation_dTCanair(iGRU),dCanopyEvaporation_dTCanopy(iGRU),dCanopyEvaporation_dTGround(iGRU),dCanopyEvaporation_dCanWat(iGRU),&
  dGroundEvaporation_dTCanair(iGRU),dGroundEvaporation_dTCanopy(iGRU),dGroundEvaporation_dTGround(iGRU),dGroundEvaporation_dCanWat(iGRU),&
  scalarCanopyLiqDeriv(iGRU), dCanLiq_dTcanopy(iGRU),&
  dNrgFlux_dTempAbove(:,iGRU), dNrgFlux_dTempBelow(:,iGRU),&
  mLayerdTrans_dTCanair(:,iGRU), mLayerdTrans_dTCanopy(:,iGRU), mLayerdTrans_dTGround(:,iGRU), mLayerdTrans_dCanWat(:,iGRU), &
  dAquiferTrans_dTCanair(iGRU), dAquiferTrans_dTCanopy(iGRU), dAquiferTrans_dTGround(iGRU), dAquiferTrans_dCanWat(iGRU),&
  iLayerLiqFluxSnowDeriv(:,iGRU),&
  dq_dHydStateAbove(:,iGRU), dq_dHydStateBelow(:,iGRU), dq_dHydStateLayerSurfVec(:,iGRU), &
  dBaseflow_dAquifer(iGRU), &
  dq_dNrgStateAbove(:,iGRU), dq_dNrgStateBelow(:,iGRU), dq_dNrgStateLayerSurfVec(:,iGRU),&
  dBaseflow_dMatric(:,:,iGRU),&
  dTemp_dTheta(:,iGRU),dTemp_dPsi0(:,iGRU),dTemp_dEnthalpy(:,iGRU),&
  dCanairTemp_dEnthalpy(iGRU), dCanopyTemp_dCanWat(iGRU), dCanopyTemp_dEnthalpy(iGRU))

end subroutine

attributes(device) subroutine computJacobWithPrime_device(aJac, dMat, dMat0, cj, &
  computeVegFlux, enthalpyStateVec, ixRichards, &
  nSnow, nSoil, nLayers, nState,&
  ixCasNrg, ixVegNrg, ixVegHyd, ixAqWat, &
  ixSnowSoilNrg, ixSnowOnlyNrg, ixSoilOnlyNrg, &
  ixSnowSoilHyd, ixSnowOnlyHyd, ixSoilOnlyHyd,&
  scalarBulkVolHeatCapVeg, dTheta_dTkCanopy, dVolHtCapBulk_dTkCanopy,&
  scalarCanopyTempPrime, dCm_dTkCanopy, scalarCanopyWatPrime,&
  canopyDepth, d2Theta_dTkCanopy2, dFracLiqVeg_dTkCanopy, &
  mLayerVolHtCapBulk, mLayerdTheta_dTk,&
  dVolHtCapBulk_dTk, mLayerTempPrime, &
  dCm_dTk, dFracLiqWat_dTk, &
  mLayerVolFracWatPrime, mLayerd2Theta_dTk2,&
  dVolTot_dPsi0, d2VolTot_dPsi02, dCompress_dPsi, &
  mLayerMatricHeadPrime, theta_sat, &
  specificStorage, &
  scalarFracLiqVeg, dVolHtCapBulk_dCanWat, scalarCanopyCm, &
  dt, dCanopyNetFlux_dCanWat, &
  mLayerFracLiqSnow, dVolHtCapBulk_dTheta, mLayerCm, &
  mLayerDepth, &
  dNrgFlux_dWatAbove, dNrgFlux_dWatBelow,&
  dCm_dPsi0, dVolHtCapBulk_dPsi0,&
  computeBaseflow, &
  ixTopNrg, ixTopHyd, ixHydType,&
  mLayerVolFracIce, scalarSoilControl,&
  dCanairNetFlux_dCanairTemp, dCanairNetFlux_dCanopyTemp, dCanairNetFlux_dGroundTemp, &
  dCanopyNetFlux_dCanairTemp, dCanopyNetFlux_dCanopyTemp, dCanopyNetFlux_dGroundTemp, &
  dGroundNetFlux_dCanairTemp, dGroundNetFlux_dCanopyTemp, dGroundNetFlux_dCanWat,&
  dCanopyEvaporation_dTCanair,dCanopyEvaporation_dTCanopy,dCanopyEvaporation_dTGround,dCanopyEvaporation_dCanWat,&
  dGroundEvaporation_dTCanair,dGroundEvaporation_dTCanopy,dGroundEvaporation_dTGround,dGroundEvaporation_dCanWat,&
  scalarCanopyLiqDeriv, dCanLiq_dTcanopy,&
  dNrgFlux_dTempAbove, dNrgFlux_dTempBelow,&
  mLayerdTrans_dTCanair, mLayerdTrans_dTCanopy, mLayerdTrans_dTGround, mLayerdTrans_dCanWat, &
  dAquiferTrans_dTCanair, dAquiferTrans_dTCanopy, dAquiferTrans_dTGround, dAquiferTrans_dCanWat,&
  iLayerLiqFluxSnowDeriv,&
  dq_dHydStateAbove, dq_dHydStateBelow, dq_dHydStateLayerSurfVec, &
  dBaseflow_dAquifer, &
  dq_dNrgStateAbove, dq_dNrgStateBelow, dq_dNrgStateLayerSurfVec,&
  dBaseflow_dMatric,&
  dTemp_dTheta,dTemp_dPsi0,dTemp_dEnthalpy,&
  dCanairTemp_dEnthalpy, dCanopyTemp_dCanWat, dCanopyTemp_dEnthalpy)
  use computJacob_module,only:fluxJacAdd
  use initialize_device,only:get_iGRU
  implicit none
  real(rkind),intent(inout) :: aJac(:,:,:)
  real(rkind),intent(inout) :: dMat(:)
  real(rkind),intent(in) :: dMat0(:)
  real(rkind),intent(in) :: cj
  logical(lgt),intent(in) :: computeVegFlux, enthalpyStateVec
  integer(i4b),intent(in) :: ixRichards
  integer(i4b),intent(in) :: nSnow, nSoil, nLayers,nState
  integer(i4b),intent(in) :: ixCasNrg, ixVegNrg, ixVegHyd, ixAqWat
  integer(i4b),intent(in) :: ixSnowSoilNrg(:), ixSnowOnlyNrg(:), ixSoilOnlyNrg(:)
  integer(i4b),intent(in) :: ixSnowSoilHyd(:), ixSnowOnlyHyd(:), ixSoilOnlyHyd(:)
  real(rkind),intent(in) :: scalarBulkVolHeatCapVeg, dTheta_dTkCanopy, dVolHtCapBulk_dTkCanopy
  real(rkind),intent(in) :: scalarCanopyTempPrime, dCm_dTkCanopy, scalarCanopyWatPrime
  real(rkind),intent(in) :: canopyDepth, d2Theta_dTkCanopy2, dFracLiqVeg_dTkCanopy
  real(rkind),intent(in) :: mLayerVolHtCapBulk(:), mLayerdTheta_dTk(:)
  real(rkind),intent(in) :: dVolHtCapBulk_dTk(:), mLayerTempPrime(:)
  real(rkind),intent(in) :: dCm_dTk(:), dFracLiqWat_dTk(:)
  real(rkind),intent(in) :: mLayerVolFracWatPrime(:), mLayerd2Theta_dTk2(:)
  real(rkind),intent(in) :: dVolTot_dPsi0(:), d2VolTot_dPsi02(:), dCompress_dPsi(:)
  real(rkind),intent(in) :: mLayerMatricHeadPrime(:), theta_sat(:)
  real(rkind),intent(in) :: specificStorage
  real(rkind),intent(in) :: scalarFracLiqVeg, dVolHtCapBulk_dCanWat, scalarCanopyCm
  real(rkind),intent(in) :: dt
  real(rkind),intent(in) :: dCanopyNetFlux_dCanWat
  real(rkind),intent(in) :: mLayerFracLiqSnow(:), dVolHtCapBulk_dTheta(:), mLayerCm(:)
  real(rkind),intent(in) :: mLayerDepth(:)
  real(rkind),intent(in) :: dNrgFlux_dWatAbove(0:), dNrgFlux_dWatBelow(0:)
  real(rkind),intent(in) :: dCm_dPsi0(:), dVolHtCapBulk_dPsi0(:)
  logical(lgt),intent(in) :: computeBaseflow
  integer(i4b),intent(in) :: ixTopNrg, ixTopHyd
  integer(i4b),intent(in) :: ixHydType(:)
  real(rkind),intent(in) :: mLayerVolFracIce(:)
  real(rkind),intent(in) :: scalarSoilControl
  real(rkind),intent(in) :: dCanairNetFlux_dCanairTemp, dCanairNetFlux_dCanopyTemp, dCanairNetFlux_dGroundTemp
  real(rkind),intent(in) :: dCanopyNetFlux_dCanairTemp,dCanopyNetFlux_dCanopyTemp,dCanopyNetFlux_dGroundTemp
  real(rkind),intent(in) :: dGroundNetFlux_dCanairTemp,dGroundNetFlux_dCanopyTemp,dGroundNetFlux_dCanWat
  real(rkind),intent(in) :: dCanopyEvaporation_dTCanair,dCanopyEvaporation_dTCanopy,dCanopyEvaporation_dTGround,dCanopyEvaporation_dCanWat
  real(rkind),intent(in) :: dGroundEvaporation_dTCanair,dGroundEvaporation_dTCanopy,dGroundEvaporation_dTGround,dGroundEvaporation_dCanWat
  real(rkind),intent(in) :: scalarCanopyLiqDeriv, dCanLiq_dTcanopy
  real(rkind),intent(in) :: dNrgFlux_dTempAbove(0:), dNrgFlux_dTempBelow(0:)
  real(rkind),intent(in) :: mLayerdTrans_dTCanair(:),mLayerdTrans_dTCanopy(:),mLayerdTrans_dTGround(:),mLayerdTrans_dCanWat(:)
  real(rkind),intent(in) :: dAquiferTrans_dTCanair,dAquiferTrans_dTCanopy,dAquiferTrans_dTGround,dAquiferTrans_dCanWat
  real(rkind),intent(in) :: iLayerLiqFluxSnowDeriv(0:)
  real(rkind),intent(in) :: dq_dHydStateAbove(0:),dq_dHydStateBelow(0:),dq_dHydStateLayerSurfVec(0:)
  real(rkind),intent(in) :: dBaseflow_dAquifer
  real(rkind),intent(in) :: dq_dNrgStateAbove(0:),dq_dNrgStateBelow(0:),dq_dNrgStateLayerSurfVec(0:)
  real(rkind),intent(in) :: dBaseflow_dMatric(:,:)
  real(rkind),intent(in) :: dTemp_dTheta(:), dTemp_dPsi0(:), dTemp_dEnthalpy(:)
  real(rkind),intent(in) :: dCanairTemp_dEnthalpy, dCanopyTemp_dCanWat, dCanopyTemp_dEnthalpy

  integer(i4b) :: iLayer, jLayer
  real(rkind) :: LH_fu0
  integer(i4b) :: nrgState, watState
  integer(i4b) :: err
  integer(i4b) :: iGRU
  iGRU = get_iGRU()

  aJac(:,:,iGRU) = 0._rkind
  dMat = dMat0 * cj ! dMat0(ixCasNrg) = Cp_air*iden_air and dMat0(Wat states) = 1.0

  if(computeVegFlux)then
    ! compute terms in the Jacobian for vegetation (excluding fluxes)
    if(ixVegNrg/=integerMissing)&
      dMat(ixVegNrg) = ( scalarBulkVolHeatCapVeg + LH_fus*iden_water*dTheta_dTkCanopy ) * cj &
                          + dVolHtCapBulk_dTkCanopy * scalarCanopyTempPrime &
                          + dCm_dTkCanopy * scalarCanopyWatPrime/canopyDepth &
                          + LH_fus*iden_water * scalarCanopyTempPrime * d2Theta_dTkCanopy2 &
                          + LH_fus * dFracLiqVeg_dTkCanopy * scalarCanopyWatPrime/canopyDepth   
  endif

  ! compute terms for the Jacobian for the snow and soil domains (excluding fluxes)
    do iLayer=1,nLayers
      if(ixSnowSoilNrg(iLayer)/=integerMissing)&
          dMat(ixSnowSoilNrg(iLayer)) = ( mLayerVolHtCapBulk(iLayer) + LH_fus*iden_water*mLayerdTheta_dTk(iLayer) ) * cj &
                                       + dVolHtCapBulk_dTk(iLayer) * mLayerTempPrime(iLayer) &
                                       + dCm_dTk(iLayer) * mLayerVolFracWatPrime(iLayer) &
                                       + LH_fus*iden_water * mLayerTempPrime(iLayer) * mLayerd2Theta_dTk2(iLayer) &
                                       + LH_fus*iden_water * dFracLiqWat_dTk(iLayer) * mLayerVolFracWatPrime(iLayer)
    end do

    ! compute terms for the Jacobian for the soil domain (excluding fluxes)
    do iLayer=1,nSoil
      if(ixSoilOnlyHyd(iLayer)/=integerMissing)then ! writes over dMat(ixSoilOnlyHyd(iLayer) = 1.0 * cj above
        dMat(ixSoilOnlyHyd(iLayer)) = ( dVolTot_dPsi0(iLayer) + dCompress_dPsi(iLayer) ) * cj + d2VolTot_dPsi02(iLayer) * mLayerMatricHeadPrime(iLayer)
        if(ixRichards==mixdform)&
            dMat(ixSoilOnlyHyd(iLayer)) = dMat(ixSoilOnlyHyd(iLayer)) + specificStorage * dVolTot_dPsi0(iLayer) * mLayerMatricHeadPrime(iLayer)/theta_sat(iLayer)
      endif
    end do

    ! if using enthalpy as a state variable, zero out usual RHS terms and add them end of the iteration loop 
    ! NOTE: other terms on RHS that are not fluxes are zeroed out by not computing heat capacity and Cm and their derivatives
    if(enthalpyStateVec)then 
      if(ixCasNrg/=integerMissing) dMat(ixCasNrg) = 0._rkind
      if(ixVegNrg/=integerMissing) dMat(ixVegNrg) = 0._rkind
      do iLayer=1,nLayers
        if(ixSnowSoilNrg(iLayer)/=integerMissing) dMat(ixSnowSoilNrg(iLayer)) = 0._rkind
      end do
      LH_fu0 = 0._rkind ! set to 0 to not use RHS terms
    else
      LH_fu0 = LH_fus ! use regular value
    endif

    ! *********************************************************************************************************************************************************
    ! * PART 1: COMPUTE CROSS-DERIVATIVE JACOBIAN TERMS 
    ! *********************************************************************************************************************************************************
    ! -----
    ! * cross derivatives in the vegetation...
    ! ---------------------------------------------
    if(computeVegFlux)then ! (derivatives only defined when vegetation protrudes over the surface)
      if(ixVegHyd/=integerMissing .and. ixVegNrg/=integerMissing)&
          ! NOTE: dIce/dLiq = (1 - scalarFracLiqVeg); dIce*LH_fu0/canopyDepth = J m-3; dLiq = kg m-2
          aJac(ixVegNrg,ixVegHyd,iGRU) = (-1._rkind + scalarFracLiqVeg)*LH_fu0/canopyDepth * cj &
                                                   + dVolHtCapBulk_dCanWat * scalarCanopyTempPrime + scalarCanopyCm/canopyDepth * cj &
                                                   - (dt/canopyDepth) * dCanopyNetFlux_dCanWat &
                                                   + LH_fu0 * scalarCanopyTempPrime * dFracLiqVeg_dTkCanopy/canopyDepth
    endif  ! if there is a need to compute energy fluxes within vegetation

        ! -----
    ! * cross derivatives in the snow domain...
    ! ----------------------------------------
    if(nSnow>0)then
      do iLayer=1,nSnow  ! loop through layers in the snow domain

        ! - check that the snow layer is desired
        if(ixSnowOnlyNrg(iLayer)==integerMissing) cycle
        ! (define the energy state)
        nrgState = ixSnowOnlyNrg(iLayer)       ! index within the full state vector
        ! - define state indices for the current layer
        watState = ixSnowOnlyHyd(iLayer)   ! hydrology state index within the state subset

        if(watState/=integerMissing)then       ! (water state for the current layer is within the state subset)
          ! - include derivatives of energy fluxes w.r.t water fluxes for current layer
          aJac(nrgState,watState,iGRU) = (-1._rkind + mLayerFracLiqSnow(iLayer))*LH_fu0*iden_water * cj &
                                      + dVolHtCapBulk_dTheta(iLayer) * mLayerTempPrime(iLayer) + mLayerCm(iLayer) * cj &
                                      + (dt/mLayerDepth(iLayer))*(-dNrgFlux_dWatBelow(iLayer-1) + dNrgFlux_dWatAbove(iLayer)) &
                                      + LH_fu0*iden_water * mLayerTempPrime(iLayer) * dFracLiqWat_dTk(iLayer)    ! (dF/dLiq)
        endif ! (if the water state for the current layer is within the state subset)

      end do ! (looping through snow layers)
    endif ! (if there are state variables for both water and energy in the snow domain)


        ! -----
    ! * cross derivatives in the soil domain...
    ! ----------------------------------------
    if(nSoil>0)then
      do iLayer=1,nSoil

        ! - check that the soil layer is desired
        if(ixSoilOnlyNrg(iLayer)==integerMissing) cycle
        ! - define indices of the soil layers
        jLayer   = iLayer+nSnow                ! index of layer in the snow+soil vector
        ! - define the energy state variable
        nrgState = ixSoilOnlyNrg(iLayer)       ! index within the full state vector
        ! - define index of hydrology state variable within the state subset
        watState = ixSoilOnlyHyd(iLayer)

        ! only compute derivatives if the water state for the current layer is within the state subset
        if(watState/=integerMissing)then
          ! - include derivatives in energy fluxes w.r.t. with respect to water for current layer
          aJac(nrgState,watState,iGRU) = dVolHtCapBulk_dPsi0(iLayer) * mLayerTempPrime(jLayer) &
                                                       + mLayerCm(jLayer) * dVolTot_dPsi0(iLayer) * cj + dCm_dPsi0(iLayer) * mLayerVolFracWatPrime(jLayer) &
                                                       + (dt/mLayerDepth(jLayer))*(-dNrgFlux_dWatBelow(jLayer-1) + dNrgFlux_dWatAbove(jLayer)) &
                                                       + mLayerCm(jLayer) * d2VolTot_dPsi02(iLayer) * mLayerMatricHeadPrime(iLayer)
          if(mLayerdTheta_dTk(jLayer) > tiny(1.0_rkind))&  ! ice is present
              aJac(nrgState,watState,iGRU) = -LH_fu0*iden_water * dVolTot_dPsi0(iLayer) * cj &
                                                       - LH_fu0*iden_water * mLayerMatricHeadPrime(iLayer) * d2VolTot_dPsi02(iLayer) + aJac(nrgState,watState,iGRU) ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
        endif ! (if the water state for the current layer is within the state subset)

      end do ! (looping through energy states in the soil domain)
    endif ! (if there are state variables for both water and energy in the soil domain)

    ! *********************************************************************************************************************************************************
    ! * PART 2: COMPUTE FLUX JACOBIAN TERMS 
    ! *********************************************************************************************************************************************************
    call fluxJacAdd(dt,nSnow,nSoil,nLayers,computeVegFlux,computeBaseflow,&
                    ixCasNrg, ixVegNrg, ixVegHyd, ixTopNrg, ixTopHyd, ixAqWat, &
                    ixSnowSoilNrg, ixSnowOnlyNrg, ixSoilOnlyNrg, &
                    ixSnowSoilHyd, ixSnowOnlyHyd, ixSoilOnlyHyd, &
                    ixHydType, &
                    mLayerDepth, mLayerVolFracIce, &
                    ! diag_data,&
                    mLayerFracLiqSnow, &
                    scalarFracLiqVeg, scalarSoilControl, canopyDepth, &
                    ! deriv_data,&
                    dCanairNetFlux_dCanairTemp,dCanairNetFlux_dCanopyTemp,dCanairNetFlux_dGroundTemp,&
                    dCanopyNetFlux_dCanairTemp,dCanopyNetFlux_dCanopyTemp,dCanopyNetFlux_dGroundTemp,&
                    dGroundNetFlux_dCanairTemp,dGroundNetFlux_dCanopyTemp,dGroundNetFlux_dCanWat, &
                    dCanopyEvaporation_dTCanair,dCanopyEvaporation_dTCanopy,dCanopyEvaporation_dTGround,dCanopyEvaporation_dCanWat,&
                    dGroundEvaporation_dTCanair,dGroundEvaporation_dTCanopy,dGroundEvaporation_dTGround,dGroundEvaporation_dCanWat, &
                    scalarCanopyLiqDeriv, dCanLiq_dTcanopy, &
                    dNrgFlux_dTempAbove, dNrgFlux_dTempBelow, &
                    dNrgFlux_dWatAbove, dNrgFlux_dWatBelow, &
                    mLayerdTrans_dTCanair,mLayerdTrans_dTCanopy,mLayerdTrans_dTGround,mLayerdTrans_dCanWat, &
                    dAquiferTrans_dTCanair,dAquiferTrans_dTCanopy,dAquiferTrans_dTGround,dAquiferTrans_dCanWat, &
                    iLayerLiqFluxSnowDeriv, &
                    dq_dHydStateAbove,dq_dHydStateBelow,dq_dHydStateLayerSurfVec, &
                    dBaseflow_dAquifer, &
                    dq_dNrgStateAbove,dq_dNrgStateBelow,dq_dNrgStateLayerSurfVec,&
                    mLayerdTheta_dTk, &
                    dBaseflow_dMatric,&
                    dMat,aJac,err)

    ! *********************************************************************************************************************************************************
    ! * PART 3: CLEAN UP JACOBIAN (IF USING ENTHALPY AS A STATE VARIABLE) AND PRINT (IF DESIRED)
    ! *********************************************************************************************************************************************************
    ! * if desired, modify to use enthalpy as a state variable instead of temperature 
    ! NOTE, dMat(Nrg states) was used as 0 and now 1._rkind * cj is added instead 
    ! ----------------------------------------
    if(enthalpyStateVec)then 

      if(ixCasNrg/=integerMissing)then
        aJac(:,ixCasNrg,iGRU) = aJac(:,ixCasNrg,iGRU) * dCanairTemp_dEnthalpy
        aJac(ixCasNrg,ixCasNrg,iGRU) = aJac(ixCasNrg,ixCasNrg,iGRU) + 1._rkind * cj
      endif
      
      if(ixVegNrg/=integerMissing)then
        if(ixVegHyd/=integerMissing) aJac(:,ixVegHyd,iGRU) = aJac(:,ixVegHyd,iGRU) + aJac(:,ixVegNrg,iGRU) * dCanopyTemp_dCanWat
        aJac(:,ixVegNrg,iGRU) = aJac(:,ixVegNrg,iGRU) * dCanopyTemp_dEnthalpy
        aJac(ixVegNrg,ixVegNrg,iGRU) = aJac(ixVegNrg,ixVegNrg,iGRU) + 1._rkind * cj
      endif
      
      if(nLayers>0)then
        do iLayer=1,nLayers
          nrgState = ixSnowSoilNrg(iLayer)       
          if(nrgState==integerMissing) cycle
          watState = ixSnowSoilHyd(iLayer)
          if(watState/=integerMissing)then 
            if(iLayer<=nSnow)then
              aJac(:,watState,iGRU) = aJac(:,watState,iGRU) + aJac(:,nrgState,iGRU) * dTemp_dTheta(iLayer)
            else
              aJac(:,watState,iGRU) = aJac(:,watState,iGRU) + aJac(:,nrgState,iGRU) * dTemp_dPsi0(iLayer-nSnow)
            endif
          endif
          aJac(:,nrgState,iGRU) = aJac(:,nrgState,iGRU) * dTemp_dEnthalpy(iLayer)
          aJac(nrgState,nrgState,iGRU) = aJac(nrgState,nrgState,iGRU) + 1._rkind * cj
        enddo
      endif
    endif

  do iLayer=nState+1,size(aJac,1)
    do jLayer=nState+1,size(aJac,2)
      aJac(iLayer,jLayer,iGRU) = 0._rkind
    end do
    aJac(iLayer,iLayer,iGRU) = 1._rkind
  end do


end subroutine

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
  eqns_data%nGRU, &
                ! input: model control
                cj,                                       & ! intent(in):    this scalar changes whenever the step size or method order changes
                1._qp,                                    & ! intent(in):    length of the time step (seconds)
                eqns_data%nSnow,                          & ! intent(in):    number of snow layers
                eqns_data%nSoil,                          & ! intent(in):    number of soil layers
                eqns_data%nLayers,                        & ! intent(in):    total number of layers
                eqns_data%computeVegFlux,                 & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                eqns_data%model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel, & ! intent(in): flag to indicate if we need to compute baseflow
                eqns_data%mpar_data%specificStorage_,                    & ! intent(in): specific storage coefficient (m-1)
                eqns_data%mpar_data%theta_sat_,                             & ! intent(in): soil porosity (-)
                eqns_data%decisions%f_Richards,                & ! intent(in): choice of option for Richards' equation
                eqns_data%model_decisions(iLookDECISIONS%nrgConserv)%iDecision.ne.closedForm,  & ! intent(in): flag if enthalpy is state variable
                ! input: data structures
                eqns_data%indx_data,                      & ! intent(in):    index data
                eqns_data%prog_data,                      & ! intent(in):    model prognostic variables for a local HRU
                eqns_data%diag_data,                      & ! intent(in):    model diagnostic variables for a local HRU
                eqns_data%deriv_data,                     & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                eqns_data%dBaseflow_dMatric,              & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
                ! input: state variables
                eqns_data%mLayerTempPrime,                & ! intent(in):    derivative value for temperature of each snow+soil layer (K)
                eqns_data%mLayerMatricHeadPrime,          & ! intent(in):    derivative value for matric head of each snow+soil layer (m)
                eqns_data%mLayerVolFracWatPrime,          & ! intent(in):    derivative value for volumetric total water content of each snow+soil layer (-)
                eqns_data%scalarCanopyTempPrime,          & ! intent(in):    derivative value for temperature of the vegetation canopy (K)
                eqns_data%scalarCanopyWatPrime,           & ! intent(in):    derivative value for total water content of the vegetation canopy (kg m-2)
                ! input-output: Jacobian and its diagonal
                eqns_data%dMat,                           & ! intent(in):    diagonal of the Jacobian matrix excluding fluxes, not depending on the state vector
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
! private function: get the index in the band-diagonal matrix or full matrix
! **********************************************************************************************************
function ixInd(jState,iState)
  implicit none
  integer(i4b),intent(in)  :: jState ! off-diagonal state
  integer(i4b),intent(in)  :: iState ! diagonal state
  integer(i4b)             :: ixInd  ! index in the band-diagonal matrix or full matrix

  if(fullMatrix) then
    ixInd = jState
  else
    ixInd = ixDiag + jState - iState
  endif
end function ixInd

end module computJacobWithPrime_module
