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
  use globalData,only:maxSnowLayers
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
  ! real(rkind),allocatable,device              :: dMat(:,:)                    ! diagonal of the Jacobian matrix excluding fluxes, depending on the state vector
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

  integer(i4b) :: iGRU, nBaseflow
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
    ! if (.not. allocated(dMat)) allocate(dMat(nState,nGRU))
    ! if (size(dMat,1).ne.nState .or. size(dMat,2).ne.nGRU) then
    !   deallocate(dMat)
    !   allocate(dMat(nState,nGRU))
    ! end if
    ! allocate(dMat(nState,nGRU)) 
    nBaseflow = 1
    if (computeBaseflow) nBaseflow = nSoil

    call computJacobWithPrime_kernel<<<blocks,threads>>>(nGRU,maxSnowLayers,nBaseflow,aJac, dMat0, cj, &
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
    ! deallocate(dMat)

    ! print the Jacobian
    ! if(any(ieee_is_nan(aJac)))then; message=trim(message)//'NaN in Jacobian';err=20; return; endif

  end associate ! end association to variables in the data structures  

end subroutine computJacobWithPrime

attributes(global)subroutine computJacobWithPrime_kernel(nGRU,maxSnowLayers,nBaseflow,aJac, dMat0, cj, &
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
  integer(i4b), value :: nGRU, nSoil,maxSnowLayers,nBaseflow
  real(rkind),intent(inout) :: aJac(2*(nSoil+maxSnowLayers)+4,2*(nSoil+maxSnowLayers)+4,nGRU)
  ! real(rkind),intent(inout) :: dMat(2*(nSoil+maxSnowLayers)+4,nGRU)
  real(rkind),intent(in) :: dMat0(2*(nSoil+maxSnowLayers)+4,nGRU)
  real(rkind),intent(in),value :: cj
  logical(lgt),intent(in) :: computeVegFlux(nGRU)
  logical(lgt),intent(in),value :: enthalpyStateVec
  integer(i4b),intent(in) :: ixRichards
  integer(i4b),intent(in) :: nSnow(nGRU), nLayers(nGRU),nState(nGRU)
  integer(i4b),intent(in) :: ixCasNrg(nGRU), ixVegNrg(nGRU), ixVegHyd(nGRU), ixAqWat(nGRU)
  integer(i4b),intent(in) :: ixSnowSoilNrg(nSoil+maxSnowLayers+1,nGRU), ixSnowOnlyNrg(maxSnowLayers+1,nGRU), ixSoilOnlyNrg(nSoil,nGRU)
  integer(i4b),intent(in) :: ixSnowSoilHyd(nSoil+maxSnowLayers+1,nGRU), ixSnowOnlyHyd(maxSnowLayers+1,nGRU), ixSoilOnlyHyd(nSoil,nGRU)
  real(rkind),intent(in) :: scalarBulkVolHeatCapVeg(nGRU), dTheta_dTkCanopy(nGRU), dVolHtCapBulk_dTkCanopy(nGRU)
  real(rkind),intent(in) :: scalarCanopyTempPrime(nGRU), dCm_dTkCanopy(nGRU), scalarCanopyWatPrime(nGRU)
  real(rkind),intent(in) :: canopyDepth(nGRU), d2Theta_dTkCanopy2(nGRU), dFracLiqVeg_dTkCanopy(nGRU)
  real(rkind),intent(in) :: mLayerVolHtCapBulk(nSoil+maxSnowLayers+1,nGRU), mLayerdTheta_dTk(nSoil+maxSnowLayers,nGRU)
  real(rkind),intent(in) :: dVolHtCapBulk_dTk(nSoil+maxSnowLayers+1,nGRU), mLayerTempPrime(nSoil+maxSnowLayers+1,nGRU)
  real(rkind),intent(in) :: dCm_dTk(nSoil+maxSnowLayers,nGRU), dFracLiqWat_dTk(nSoil+maxSnowLayers,nGRU)
  real(rkind),intent(in) :: mLayerVolFracWatPrime(nSoil+maxSnowLayers+1,nGRU), mLayerd2Theta_dTk2(nSoil,nGRU)
  real(rkind),intent(in) :: dVolTot_dPsi0(nSoil,nGRU), d2VolTot_dPsi02(nSoil,nGRU), dCompress_dPsi(nSoil,nGRU)
  real(rkind),intent(in) :: mLayerMatricHeadPrime(nSoil,nGRU), theta_sat(nSoil,nGRU)
  real(rkind),intent(in) :: specificStorage(nGRU)
  real(rkind),intent(in) :: scalarFracLiqVeg(nGRU), dVolHtCapBulk_dCanWat(nGRU), scalarCanopyCm(nGRU)
  real(rkind),intent(in),value :: dt
  real(rkind),intent(in) :: dCanopyNetFlux_dCanWat(nGRU)
  real(rkind),intent(in) :: mLayerFracLiqSnow(maxSnowLayers+1,nGRU), dVolHtCapBulk_dTheta(nSoil+maxSnowLayers,nGRU), mLayerCm(nSoil+maxSnowLayers+1,nGRU)
  real(rkind),intent(in) :: mLayerDepth(nSoil+maxSnowLayers+1,nGRU)
  real(rkind),intent(in) :: dNrgFlux_dWatAbove(0:nSoil+maxSnowLayers,nGRU), dNrgFlux_dWatBelow(0:nSoil+maxSnowLayers,nGRU)
  real(rkind),intent(in) :: dCm_dPsi0(nSoil,nGRU), dVolHtCapBulk_dPsi0(nSoil,nGRU)
  logical(lgt),intent(in),value :: computeBaseflow
  integer(i4b),intent(in) :: ixTopNrg(nGRU), ixTopHyd(nGRU)
  integer(i4b),intent(in) :: ixHydType(nSoil+maxSnowLayers+1,nGRU)
  real(rkind),intent(in) :: mLayerVolFracIce(nSoil+maxSnowLayers+1,nGRU)
  real(rkind),intent(in) :: scalarSoilControl(nGRU)
  real(rkind),intent(in) :: dCanairNetFlux_dCanairTemp(nGRU), dCanairNetFlux_dCanopyTemp(nGRU), dCanairNetFlux_dGroundTemp(nGRU)
  real(rkind),intent(in) :: dCanopyNetFlux_dCanairTemp(nGRU),dCanopyNetFlux_dCanopyTemp(nGRU),dCanopyNetFlux_dGroundTemp(nGRU)
  real(rkind),intent(in) :: dGroundNetFlux_dCanairTemp(nGRU),dGroundNetFlux_dCanopyTemp(nGRU),dGroundNetFlux_dCanWat(nGRU)
  real(rkind),intent(in) :: dCanopyEvaporation_dTCanair(nGRU),dCanopyEvaporation_dTCanopy(nGRU),dCanopyEvaporation_dTGround(nGRU),dCanopyEvaporation_dCanWat(nGRU)
  real(rkind),intent(in) :: dGroundEvaporation_dTCanair(nGRU),dGroundEvaporation_dTCanopy(nGRU),dGroundEvaporation_dTGround(nGRU),dGroundEvaporation_dCanWat(nGRU)
  real(rkind),intent(in) :: scalarCanopyLiqDeriv(nGRU), dCanLiq_dTcanopy(nGRU)
  real(rkind),intent(in) :: dNrgFlux_dTempAbove(0:nSoil+maxSnowLayers,nGRU), dNrgFlux_dTempBelow(0:nSoil+maxSnowLayers,nGRU)
  real(rkind),intent(in) :: mLayerdTrans_dTCanair(nSoil,nGRU),mLayerdTrans_dTCanopy(nSoil,nGRU),mLayerdTrans_dTGround(nSoil,nGRU),mLayerdTrans_dCanWat(nSoil,nGRU)
  real(rkind),intent(in) :: dAquiferTrans_dTCanair(nGRU),dAquiferTrans_dTCanopy(nGRU),dAquiferTrans_dTGround(nGRU),dAquiferTrans_dCanWat(nGRU)
  real(rkind),intent(in) :: iLayerLiqFluxSnowDeriv(0:maxSnowLayers,nGRU)
  real(rkind),intent(in) :: dq_dHydStateAbove(0:nSoil,nGRU),dq_dHydStateBelow(0:nSoil,nGRU),dq_dHydStateLayerSurfVec(0:nSoil,nGRU)
  real(rkind),intent(in) :: dBaseflow_dAquifer(nGRU)

  real(rkind),intent(in) :: dq_dNrgStateAbove(0:nSoil,nGRU),dq_dNrgStateBelow(0:nSoil,nGRU),dq_dNrgStateLayerSurfVec(0:nSoil,nGRU)
  real(rkind),intent(in) :: dBaseflow_dMatric(nBaseflow,nBaseflow,nGRU)
  real(rkind),intent(in) :: dTemp_dTheta(nSoil+maxSnowLayers,nGRU), dTemp_dPsi0(nSoil,nGRU), dTemp_dEnthalpy(nSoil+maxSnowLayers,nGRU)
  real(rkind),intent(in) :: dCanairTemp_dEnthalpy(nGRU), dCanopyTemp_dCanWat(nGRU), dCanopyTemp_dEnthalpy(nGRU)

    integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

  call computJacobWithPrime_device(aJac, dMat0, cj, &
  computeVegFlux, enthalpyStateVec, ixRichards, &
  nSnow, nSoil, nLayers, nState,&
  ixCasNrg, ixVegNrg, ixVegHyd, ixAqWat, &
  ixSnowSoilNrg, ixSnowOnlyNrg, ixSoilOnlyNrg, &
  ixSnowSoilHyd, ixSnowOnlyHyd, ixSoilOnlyHyd,&
  scalarBulkVolHeatCapVeg, dTheta_dTkCanopy, dVolHtCapBulk_dTkCanopy,&
  scalarCanopyTempPrime, dCm_dTkCanopy, scalarCanopyWatPrime,&
  canopyDepth, d2Theta_dTkCanopy2, dFracLiqVeg_dTkCanopy,&
  mLayerVolHtCapBulk, mLayerdTheta_dTk,&
  dVolHtCapBulk_dTk, mLayerTempPrime, &
  dCm_dTk, dFracLiqWat_dTk, &
  mLayerVolFracWatPrime, mLayerd2Theta_dTk2, &
  dVolTot_dPsi0, d2VolTot_dPsi02, dCompress_dPsi, &
  mLayerMatricHeadPrime, theta_sat, &
  specificStorage, &
  scalarFracLiqVeg, dVolHtCapBulk_dCanWat, scalarCanopyCm, &
  dt, dCanopyNetFlux_dCanWat,&
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

end subroutine

attributes(device) subroutine computJacobWithPrime_device(aJac, dMat0, cj, &
  computeVegFlux_, enthalpyStateVec, ixRichards, &
  nSnow_, nSoil, nLayers_, nState_,&
  ixCasNrg_, ixVegNrg_, ixVegHyd_, ixAqWat_, &
  ixSnowSoilNrg_, ixSnowOnlyNrg_, ixSoilOnlyNrg_, &
  ixSnowSoilHyd_, ixSnowOnlyHyd_, ixSoilOnlyHyd_,&
  scalarBulkVolHeatCapVeg_, dTheta_dTkCanopy_, dVolHtCapBulk_dTkCanopy_,&
  scalarCanopyTempPrime_, dCm_dTkCanopy_, scalarCanopyWatPrime_,&
  canopyDepth_, d2Theta_dTkCanopy2_, dFracLiqVeg_dTkCanopy_, &
  mLayerVolHtCapBulk_, mLayerdTheta_dTk_,&
  dVolHtCapBulk_dTk_, mLayerTempPrime_, &
  dCm_dTk_, dFracLiqWat_dTk_, &
  mLayerVolFracWatPrime_, mLayerd2Theta_dTk2_,&
  dVolTot_dPsi0_, d2VolTot_dPsi02_, dCompress_dPsi_, &
  mLayerMatricHeadPrime_, theta_sat_, &
  specificStorage_, &
  scalarFracLiqVeg_, dVolHtCapBulk_dCanWat_, scalarCanopyCm_, &
  dt, dCanopyNetFlux_dCanWat_, &
  mLayerFracLiqSnow_, dVolHtCapBulk_dTheta_, mLayerCm_, &
  mLayerDepth_, &
  dNrgFlux_dWatAbove_, dNrgFlux_dWatBelow_,&
  dCm_dPsi0_, dVolHtCapBulk_dPsi0_,&
  computeBaseflow, &
  ixTopNrg_, ixTopHyd_, ixHydType_,&
  mLayerVolFracIce_, scalarSoilControl_,&
  dCanairNetFlux_dCanairTemp_, dCanairNetFlux_dCanopyTemp_, dCanairNetFlux_dGroundTemp_, &
  dCanopyNetFlux_dCanairTemp_, dCanopyNetFlux_dCanopyTemp_, dCanopyNetFlux_dGroundTemp_, &
  dGroundNetFlux_dCanairTemp_, dGroundNetFlux_dCanopyTemp_, dGroundNetFlux_dCanWat_,&
  dCanopyEvaporation_dTCanair_,dCanopyEvaporation_dTCanopy_,dCanopyEvaporation_dTGround_,dCanopyEvaporation_dCanWat_,&
  dGroundEvaporation_dTCanair_,dGroundEvaporation_dTCanopy_,dGroundEvaporation_dTGround_,dGroundEvaporation_dCanWat_,&
  scalarCanopyLiqDeriv_, dCanLiq_dTcanopy_,&
  dNrgFlux_dTempAbove_, dNrgFlux_dTempBelow_,&
  mLayerdTrans_dTCanair_, mLayerdTrans_dTCanopy_, mLayerdTrans_dTGround_, mLayerdTrans_dCanWat_, &
  dAquiferTrans_dTCanair_, dAquiferTrans_dTCanopy_, dAquiferTrans_dTGround_, dAquiferTrans_dCanWat_,&
  iLayerLiqFluxSnowDeriv_,&
  dq_dHydStateAbove_, dq_dHydStateBelow_, dq_dHydStateLayerSurfVec_, &
  dBaseflow_dAquifer_, &
  dq_dNrgStateAbove_, dq_dNrgStateBelow_, dq_dNrgStateLayerSurfVec_,&
  dBaseflow_dMatric_,&
  dTemp_dTheta_,dTemp_dPsi0_,dTemp_dEnthalpy_,&
  dCanairTemp_dEnthalpy_, dCanopyTemp_dCanWat_, dCanopyTemp_dEnthalpy_)
  use computJacob_module,only:fluxJacAdd
  use initialize_device,only:get_iGRU
  implicit none
  real(rkind),intent(inout) :: aJac(:,:,:)
  ! real(rkind),intent(inout) :: dMat(:,:)
  real(rkind),intent(in) :: dMat0(:,:)
  real(rkind),intent(in) :: cj
  logical(lgt),intent(in) :: computeVegFlux_(:), enthalpyStateVec
  integer(i4b),intent(in) :: ixRichards
  integer(i4b),intent(in) :: nSnow_(:), nSoil, nLayers_(:),nState_(:)
  integer(i4b),intent(in) :: ixCasNrg_(:), ixVegNrg_(:), ixVegHyd_(:), ixAqWat_(:)
  integer(i4b),intent(in) :: ixSnowSoilNrg_(:,:), ixSnowOnlyNrg_(:,:), ixSoilOnlyNrg_(:,:)
  integer(i4b),intent(in) :: ixSnowSoilHyd_(:,:), ixSnowOnlyHyd_(:,:), ixSoilOnlyHyd_(:,:)
  real(rkind),intent(in) :: scalarBulkVolHeatCapVeg_(:), dTheta_dTkCanopy_(:), dVolHtCapBulk_dTkCanopy_(:)
  real(rkind),intent(in) :: scalarCanopyTempPrime_(:), dCm_dTkCanopy_(:), scalarCanopyWatPrime_(:)
  real(rkind),intent(in) :: canopyDepth_(:), d2Theta_dTkCanopy2_(:), dFracLiqVeg_dTkCanopy_(:)
  real(rkind),intent(in) :: mLayerVolHtCapBulk_(:,:), mLayerdTheta_dTk_(:,:)
  real(rkind),intent(in) :: dVolHtCapBulk_dTk_(:,:), mLayerTempPrime_(:,:)
  real(rkind),intent(in) :: dCm_dTk_(:,:), dFracLiqWat_dTk_(:,:)
  real(rkind),intent(in) :: mLayerVolFracWatPrime_(:,:), mLayerd2Theta_dTk2_(:,:)
  real(rkind),intent(in) :: dVolTot_dPsi0_(:,:), d2VolTot_dPsi02_(:,:), dCompress_dPsi_(:,:)
  real(rkind),intent(in) :: mLayerMatricHeadPrime_(:,:), theta_sat_(:,:)
  real(rkind),intent(in) :: specificStorage_(:)
  real(rkind),intent(in) :: scalarFracLiqVeg_(:), dVolHtCapBulk_dCanWat_(:), scalarCanopyCm_(:)
  real(rkind),intent(in) :: dt
  real(rkind),intent(in) :: dCanopyNetFlux_dCanWat_(:)
  real(rkind),intent(in) :: mLayerFracLiqSnow_(:,:), dVolHtCapBulk_dTheta_(:,:), mLayerCm_(:,:)
  real(rkind),intent(in) :: mLayerDepth_(:,:)
  real(rkind),intent(in) :: dNrgFlux_dWatAbove_(0:,:), dNrgFlux_dWatBelow_(0:,:)
  real(rkind),intent(in) :: dCm_dPsi0_(:,:), dVolHtCapBulk_dPsi0_(:,:)
  logical(lgt),intent(in) :: computeBaseflow
  integer(i4b),intent(in) :: ixTopNrg_(:), ixTopHyd_(:)
  integer(i4b),intent(in) :: ixHydType_(:,:)
  real(rkind),intent(in) :: mLayerVolFracIce_(:,:)
  real(rkind),intent(in) :: scalarSoilControl_(:)
  real(rkind),intent(in) :: dCanairNetFlux_dCanairTemp_(:), dCanairNetFlux_dCanopyTemp_(:), dCanairNetFlux_dGroundTemp_(:)
  real(rkind),intent(in) :: dCanopyNetFlux_dCanairTemp_(:),dCanopyNetFlux_dCanopyTemp_(:),dCanopyNetFlux_dGroundTemp_(:)
  real(rkind),intent(in) :: dGroundNetFlux_dCanairTemp_(:),dGroundNetFlux_dCanopyTemp_(:),dGroundNetFlux_dCanWat_(:)
  real(rkind),intent(in) :: dCanopyEvaporation_dTCanair_(:),dCanopyEvaporation_dTCanopy_(:),dCanopyEvaporation_dTGround_(:),dCanopyEvaporation_dCanWat_(:)
  real(rkind),intent(in) :: dGroundEvaporation_dTCanair_(:),dGroundEvaporation_dTCanopy_(:),dGroundEvaporation_dTGround_(:),dGroundEvaporation_dCanWat_(:)
  real(rkind),intent(in) :: scalarCanopyLiqDeriv_(:), dCanLiq_dTcanopy_(:)
  real(rkind),intent(in) :: dNrgFlux_dTempAbove_(0:,:), dNrgFlux_dTempBelow_(0:,:)
  real(rkind),intent(in) :: mLayerdTrans_dTCanair_(:,:),mLayerdTrans_dTCanopy_(:,:),mLayerdTrans_dTGround_(:,:),mLayerdTrans_dCanWat_(:,:)
  real(rkind),intent(in) :: dAquiferTrans_dTCanair_(:),dAquiferTrans_dTCanopy_(:),dAquiferTrans_dTGround_(:),dAquiferTrans_dCanWat_(:)
  real(rkind),intent(in) :: iLayerLiqFluxSnowDeriv_(0:,:)
  real(rkind),intent(in) :: dq_dHydStateAbove_(0:,:),dq_dHydStateBelow_(0:,:),dq_dHydStateLayerSurfVec_(0:,:)
  real(rkind),intent(in) :: dBaseflow_dAquifer_(:)
  real(rkind),intent(in) :: dq_dNrgStateAbove_(0:,:),dq_dNrgStateBelow_(0:,:),dq_dNrgStateLayerSurfVec_(0:,:)
  real(rkind),intent(in) :: dBaseflow_dMatric_(:,:,:)
  real(rkind),intent(in) :: dTemp_dTheta_(:,:), dTemp_dPsi0_(:,:), dTemp_dEnthalpy_(:,:)
  real(rkind),intent(in) :: dCanairTemp_dEnthalpy_(:), dCanopyTemp_dCanWat_(:), dCanopyTemp_dEnthalpy_(:)

  integer(i4b) :: iLayer, jLayer
  real(rkind) :: LH_fu0
  integer(i4b) :: nrgState, watState
  integer(i4b) :: err
  integer(i4b) :: iGRU
  iGRU = get_iGRU()

  ! initialize the Jacobian and diagonal
  ! NOTE: this needs to be done every time, since Jacobian matrix is modified in the solver and dMat is modified below
  aJac(:,:,iGRU) = 0._rkind
  ! dMat(:,iGRU) = dMat0(:,iGRU) * cj ! dMat0(ixCasNrg) = Cp_air*iden_air and dMat0(Wat states) = 1.0
  do iLayer=1,nState_(iGRU)
    aJac(iLayer,iLayer,iGRU) = dMat0(iLayer,iGRU) * cj
  end do

  if(computeVegFlux_(iGRU))then
    ! compute terms in the Jacobian for vegetation (excluding fluxes)
    if(ixVegNrg_(iGRU)/=integerMissing)&
      aJac(ixVegNrg_(iGRU),ixVegNrg_(iGRU),iGRU) = ( scalarBulkVolHeatCapVeg_(iGRU) + LH_fus*iden_water*dTheta_dTkCanopy_(iGRU) ) * cj &
                          + dVolHtCapBulk_dTkCanopy_(iGRU) * scalarCanopyTempPrime_(iGRU) &
                          + dCm_dTkCanopy_(iGRU) * scalarCanopyWatPrime_(iGRU)/canopyDepth_(iGRU) &
                          + LH_fus*iden_water * scalarCanopyTempPrime_(iGRU) * d2Theta_dTkCanopy2_(iGRU) &
                          + LH_fus * dFracLiqVeg_dTkCanopy_(iGRU) * scalarCanopyWatPrime_(iGRU)/canopyDepth_(iGRU)   
  endif

  ! compute terms for the Jacobian for the snow and soil domains (excluding fluxes)
    do iLayer=1,nLayers_(iGRU)
      if(ixSnowSoilNrg_(iLayer,iGRU)/=integerMissing)&
          aJac(ixSnowSoilNrg_(iLayer,iGRU),ixSnowSoilNrg_(iLayer,iGRU),iGRU) = ( mLayerVolHtCapBulk_(iLayer,iGRU) + LH_fus*iden_water*mLayerdTheta_dTk_(iLayer,iGRU) ) * cj &
                                       + dVolHtCapBulk_dTk_(iLayer,iGRU) * mLayerTempPrime_(iLayer,iGRU) &
                                       + dCm_dTk_(iLayer,iGRU) * mLayerVolFracWatPrime_(iLayer,iGRU) &
                                       + LH_fus*iden_water * mLayerTempPrime_(iLayer,iGRU) * mLayerd2Theta_dTk2_(iLayer,iGRU) &
                                       + LH_fus*iden_water * dFracLiqWat_dTk_(iLayer,iGRU) * mLayerVolFracWatPrime_(iLayer,iGRU)
    end do

    ! compute terms for the Jacobian for the soil domain (excluding fluxes)
    do iLayer=1,nSoil
      if(ixSoilOnlyHyd_(iLayer,iGRU)/=integerMissing)then ! writes over dMat(ixSoilOnlyHyd(iLayer) = 1.0 * cj above
        aJac(ixSoilOnlyHyd_(iLayer,iGRU),ixSoilOnlyHyd_(iLayer,iGRU),iGRU) = ( dVolTot_dPsi0_(iLayer,iGRU) + dCompress_dPsi_(iLayer,iGRU) ) * cj + d2VolTot_dPsi02_(iLayer,iGRU) * mLayerMatricHeadPrime_(iLayer,iGRU)
        if(ixRichards==mixdform)&
            aJac(ixSoilOnlyHyd_(iLayer,iGRU),ixSoilOnlyHyd_(iLayer,iGRU),iGRU) = aJac(ixSoilOnlyHyd_(iLayer,iGRU),ixSoilOnlyHyd_(iLayer,iGRU),iGRU) + specificStorage_(iGRU) * dVolTot_dPsi0_(iLayer,iGRU) * mLayerMatricHeadPrime_(iLayer,iGRU)/theta_sat_(iLayer,iGRU)
      endif
    end do

    ! if using enthalpy as a state variable, zero out usual RHS terms and add them end of the iteration loop 
    ! NOTE: other terms on RHS that are not fluxes are zeroed out by not computing heat capacity and Cm and their derivatives
    ! if(enthalpyStateVec)then 
      if(ixCasNrg_(iGRU)/=integerMissing) aJac(ixCasNrg_(iGRU),ixCasNrg_(iGRU),iGRU) = 0._rkind
      if(ixVegNrg_(iGRU)/=integerMissing) aJac(ixVegNrg_(iGRU),ixVegNrg_(iGRU),iGRU) = 0._rkind
      do iLayer=1,nLayers_(iGRU)
        if(ixSnowSoilNrg_(iLayer,iGRU)/=integerMissing) aJac(ixSnowSoilNrg_(iLayer,iGRU),ixSnowSoilNrg_(iLayer,iGRU),iGRU) = 0._rkind
      end do
      LH_fu0 = 0._rkind ! set to 0 to not use RHS terms
    ! else
    !   LH_fu0 = LH_fus ! use regular value
    ! endif

    ! *********************************************************************************************************************************************************
    ! * PART 1: COMPUTE CROSS-DERIVATIVE JACOBIAN TERMS 
    ! *********************************************************************************************************************************************************
    ! -----
    ! * cross derivatives in the vegetation...
    ! ---------------------------------------------
    if(computeVegFlux_(iGRU))then ! (derivatives only defined when vegetation protrudes over the surface)
      if(ixVegHyd_(iGRU)/=integerMissing .and. ixVegNrg_(iGRU)/=integerMissing)&
          ! NOTE: dIce/dLiq = (1 - scalarFracLiqVeg); dIce*LH_fu0/canopyDepth = J m-3; dLiq = kg m-2
          aJac(ixVegNrg_(iGRU),ixVegHyd_(iGRU),iGRU) = (-1._rkind + scalarFracLiqVeg_(iGRU))*LH_fu0/canopyDepth_(iGRU) * cj &
                                                   + dVolHtCapBulk_dCanWat_(iGRU) * scalarCanopyTempPrime_(iGRU) + scalarCanopyCm_(iGRU)/canopyDepth_(iGRU) * cj &
                                                   - (dt/canopyDepth_(iGRU)) * dCanopyNetFlux_dCanWat_(iGRU) &
                                                   + LH_fu0 * scalarCanopyTempPrime_(iGRU) * dFracLiqVeg_dTkCanopy_(iGRU)/canopyDepth_(iGRU)
    endif  ! if there is a need to compute energy fluxes within vegetation

        ! -----
    ! * cross derivatives in the snow domain...
    ! ----------------------------------------
    ! if(nSnow_(iGRU)>0)then
      do iLayer=1,nSnow_(iGRU)  ! loop through layers in the snow domain

        ! - check that the snow layer is desired
        if(ixSnowOnlyNrg_(iLayer,iGRU)==integerMissing) cycle
        ! (define the energy state)
        nrgState = ixSnowOnlyNrg_(iLayer,iGRU)       ! index within the full state vector
        ! - define state indices for the current layer
        watState = ixSnowOnlyHyd_(iLayer,iGRU)   ! hydrology state index within the state subset

        if(watState/=integerMissing)then       ! (water state for the current layer is within the state subset)
          ! - include derivatives of energy fluxes w.r.t water fluxes for current layer
          aJac(nrgState,watState,iGRU) = (-1._rkind + mLayerFracLiqSnow_(iLayer,iGRU))*LH_fu0*iden_water * cj &
                                      + dVolHtCapBulk_dTheta_(iLayer,iGRU) * mLayerTempPrime_(iLayer,iGRU) + mLayerCm_(iLayer,iGRU) * cj &
                                      + (dt/mLayerDepth_(iLayer,iGRU))*(-dNrgFlux_dWatBelow_(iLayer-1,iGRU) + dNrgFlux_dWatAbove_(iLayer,iGRU)) &
                                      + LH_fu0*iden_water * mLayerTempPrime_(iLayer,iGRU) * dFracLiqWat_dTk_(iLayer,iGRU)    ! (dF/dLiq)
        endif ! (if the water state for the current layer is within the state subset)

      end do ! (looping through snow layers)
    ! endif ! (if there are state variables for both water and energy in the snow domain)


        ! -----
    ! * cross derivatives in the soil domain...
    ! ----------------------------------------
    ! if(nSoil>0)then
      do iLayer=1,nSoil

        ! - check that the soil layer is desired
        if(ixSoilOnlyNrg_(iLayer,iGRU)==integerMissing) cycle
        ! - define indices of the soil layers
        jLayer   = iLayer+nSnow_(iGRU)                ! index of layer in the snow+soil vector
        ! - define the energy state variable
        nrgState = ixSoilOnlyNrg_(iLayer,iGRU)       ! index within the full state vector
        ! - define index of hydrology state variable within the state subset
        watState = ixSoilOnlyHyd_(iLayer,iGRU)

        ! only compute derivatives if the water state for the current layer is within the state subset
        if(watState/=integerMissing)then
          ! - include derivatives in energy fluxes w.r.t. with respect to water for current layer
          aJac(nrgState,watState,iGRU) = dVolHtCapBulk_dPsi0_(iLayer,iGRU) * mLayerTempPrime_(jLayer,iGRU) &
                                                       + mLayerCm_(jLayer,iGRU) * dVolTot_dPsi0_(iLayer,iGRU) * cj + dCm_dPsi0_(iLayer,iGRU) * mLayerVolFracWatPrime_(jLayer,iGRU) &
                                                       + (dt/mLayerDepth_(jLayer,iGRU))*(-dNrgFlux_dWatBelow_(jLayer-1,iGRU) + dNrgFlux_dWatAbove_(jLayer,iGRU)) &
                                                       + mLayerCm_(jLayer,iGRU) * d2VolTot_dPsi02_(iLayer,iGRU) * mLayerMatricHeadPrime_(iLayer,iGRU)
          if(mLayerdTheta_dTk_(jLayer,iGRU) > tiny(1.0_rkind))&  ! ice is present
              aJac(nrgState,watState,iGRU) = -LH_fu0*iden_water * dVolTot_dPsi0_(iLayer,iGRU) * cj &
                                                       - LH_fu0*iden_water * mLayerMatricHeadPrime_(iLayer,iGRU) * d2VolTot_dPsi02_(iLayer,iGRU) + aJac(nrgState,watState,iGRU) ! dNrg/dMat (J m-3 m-1) -- dMat changes volumetric water, and hence ice content
        endif ! (if the water state for the current layer is within the state subset)

      end do ! (looping through energy states in the soil domain)
    ! endif ! (if there are state variables for both water and energy in the soil domain)

    ! *********************************************************************************************************************************************************
    ! * PART 2: COMPUTE FLUX JACOBIAN TERMS 
    ! *********************************************************************************************************************************************************
    call fluxJacAdd(dt,nSnow_,nSoil,nLayers_,computeVegFlux_,computeBaseflow,&
                    ixCasNrg_, ixVegNrg_, ixVegHyd_, ixTopNrg_, ixTopHyd_, ixAqWat_, &
                    ixSnowSoilNrg_, ixSnowOnlyNrg_, ixSoilOnlyNrg_, &
                    ixSnowSoilHyd_, ixSnowOnlyHyd_, ixSoilOnlyHyd_, &
                    ixHydType_, &
                    mLayerDepth_, mLayerVolFracIce_, &
                    ! diag_data,&
                    mLayerFracLiqSnow_, &
                    scalarFracLiqVeg_, scalarSoilControl_, canopyDepth_, &
                    ! deriv_data,&
                    dCanairNetFlux_dCanairTemp_,dCanairNetFlux_dCanopyTemp_,dCanairNetFlux_dGroundTemp_,&
                    dCanopyNetFlux_dCanairTemp_,dCanopyNetFlux_dCanopyTemp_,dCanopyNetFlux_dGroundTemp_,&
                    dGroundNetFlux_dCanairTemp_,dGroundNetFlux_dCanopyTemp_,dGroundNetFlux_dCanWat_, &
                    dCanopyEvaporation_dTCanair_,dCanopyEvaporation_dTCanopy_,dCanopyEvaporation_dTGround_,dCanopyEvaporation_dCanWat_,&
                    dGroundEvaporation_dTCanair_,dGroundEvaporation_dTCanopy_,dGroundEvaporation_dTGround_,dGroundEvaporation_dCanWat_, &
                    scalarCanopyLiqDeriv_, dCanLiq_dTcanopy_, &
                    dNrgFlux_dTempAbove_, dNrgFlux_dTempBelow_, &
                    dNrgFlux_dWatAbove_, dNrgFlux_dWatBelow_, &
                    mLayerdTrans_dTCanair_,mLayerdTrans_dTCanopy_,mLayerdTrans_dTGround_,mLayerdTrans_dCanWat_, &
                    dAquiferTrans_dTCanair_,dAquiferTrans_dTCanopy_,dAquiferTrans_dTGround_,dAquiferTrans_dCanWat_, &
                    iLayerLiqFluxSnowDeriv_, &
                    dq_dHydStateAbove_,dq_dHydStateBelow_,dq_dHydStateLayerSurfVec_, &
                    dBaseflow_dAquifer_, &
                    dq_dNrgStateAbove_,dq_dNrgStateBelow_,dq_dNrgStateLayerSurfVec_,&
                    mLayerdTheta_dTk_, &
                    dBaseflow_dMatric_,&
                    dMat0,aJac,err)

    ! *********************************************************************************************************************************************************
    ! * PART 3: CLEAN UP JACOBIAN (IF USING ENTHALPY AS A STATE VARIABLE) AND PRINT (IF DESIRED)
    ! *********************************************************************************************************************************************************
    ! * if desired, modify to use enthalpy as a state variable instead of temperature 
    ! NOTE, dMat(Nrg states) was used as 0 and now 1._rkind * cj is added instead 
    ! ----------------------------------------
    ! if(enthalpyStateVec)then 

      if(ixCasNrg_(iGRU)/=integerMissing)then
        aJac(:,ixCasNrg_(iGRU),iGRU) = aJac(:,ixCasNrg_(iGRU),iGRU) * dCanairTemp_dEnthalpy_(iGRU)
        aJac(ixCasNrg_(iGRU),ixCasNrg_(iGRU),iGRU) = aJac(ixCasNrg_(iGRU),ixCasNrg_(iGRU),iGRU) + 1._rkind * cj
      endif
      
      if(ixVegNrg_(iGRU)/=integerMissing)then
        if(ixVegHyd_(iGRU)/=integerMissing) then
          do iLayer=1,size(aJac,1)
            aJac(iLayer,ixVegHyd_(iGRU),iGRU) = aJac(iLayer,ixVegHyd_(iGRU),iGRU) + aJac(iLayer,ixVegNrg_(iGRU),iGRU) * dCanopyTemp_dCanWat_(iGRU)
          end do
        end if
          ! aJac(:,ixVegHyd_(iGRU),iGRU) = aJac(:,ixVegHyd_(iGRU),iGRU) + aJac(:,ixVegNrg_(iGRU),iGRU) * dCanopyTemp_dCanWat
        aJac(:,ixVegNrg_(iGRU),iGRU) = aJac(:,ixVegNrg_(iGRU),iGRU) * dCanopyTemp_dEnthalpy_(iGRU)
        aJac(ixVegNrg_(iGRU),ixVegNrg_(iGRU),iGRU) = aJac(ixVegNrg_(iGRU),ixVegNrg_(iGRU),iGRU) + 1._rkind * cj
      endif
      
      ! if(nLayers_(iGRU)>0)then
        do iLayer=1,nLayers_(iGRU)
          nrgState = ixSnowSoilNrg_(iLayer,iGRU)       
          if(nrgState==integerMissing) cycle
          watState = ixSnowSoilHyd_(iLayer,iGRU)
          if(watState/=integerMissing)then 
            if(iLayer<=nSnow_(iGRU))then
              aJac(:,watState,iGRU) = aJac(:,watState,iGRU) + aJac(:,nrgState,iGRU) * dTemp_dTheta_(iLayer,iGRU)
            else
              aJac(:,watState,iGRU) = aJac(:,watState,iGRU) + aJac(:,nrgState,iGRU) * dTemp_dPsi0_(iLayer-nSnow_(iGRU),iGRU)
            endif
          endif
          aJac(:,nrgState,iGRU) = aJac(:,nrgState,iGRU) * dTemp_dEnthalpy_(iLayer,iGRU)
          aJac(nrgState,nrgState,iGRU) = aJac(nrgState,nrgState,iGRU) + 1._rkind * cj
        enddo
      ! endif
    ! endif

  do iLayer=nState_(iGRU)+1,size(aJac,1)
    do jLayer=nState_(iGRU)+1,size(aJac,2)
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
                eqns_data%nSnow_p,                          & ! intent(in):    number of snow layers
                eqns_data%nSoil,                          & ! intent(in):    number of soil layers
                eqns_data%nLayers_p,                        & ! intent(in):    total number of layers
                eqns_data%computeVegFlux_p,                 & ! intent(in):    flag to indicate if we need to compute fluxes over vegetation
                eqns_data%model_decisions(iLookDECISIONS%groundwatr)%iDecision==qbaseTopmodel, & ! intent(in): flag to indicate if we need to compute baseflow
                eqns_data%mpar_data_p%specificStorage_,                    & ! intent(in): specific storage coefficient (m-1)
                eqns_data%mpar_data_p%theta_sat_,                             & ! intent(in): soil porosity (-)
                eqns_data%decisions_p%f_Richards,                & ! intent(in): choice of option for Richards' equation
                eqns_data%model_decisions(iLookDECISIONS%nrgConserv)%iDecision.ne.closedForm,  & ! intent(in): flag if enthalpy is state variable
                ! input: data structures
                eqns_data%indx_data_p,                      & ! intent(in):    index data
                eqns_data%prog_data_p,                      & ! intent(in):    model prognostic variables for a local HRU
                eqns_data%diag_data_p,                      & ! intent(in):    model diagnostic variables for a local HRU
                eqns_data%deriv_data_p,                     & ! intent(in):    derivatives in model fluxes w.r.t. relevant state variables
                eqns_data%dBaseflow_dMatric,              & ! intent(in):    derivative in baseflow w.r.t. matric head (s-1)
                ! input: state variables
                eqns_data%mLayerTempPrime,                & ! intent(in):    derivative value for temperature of each snow+soil layer (K)
                eqns_data%mLayerMatricHeadPrime,          & ! intent(in):    derivative value for matric head of each snow+soil layer (m)
                eqns_data%mLayerVolFracWatPrime,          & ! intent(in):    derivative value for volumetric total water content of each snow+soil layer (-)
                eqns_data%scalarCanopyTempPrime,          & ! intent(in):    derivative value for temperature of the vegetation canopy (K)
                eqns_data%scalarCanopyWatPrime,           & ! intent(in):    derivative value for total water content of the vegetation canopy (kg m-2)
                ! input-output: Jacobian and its diagonal
                eqns_data%dMat_p,                           & ! intent(in):    diagonal of the Jacobian matrix excluding fluxes, not depending on the state vector
                Jac,                                      & ! intent(out):   Jacobian matrix
                ! output: error control
                eqns_data%err,eqns_data%message)            ! intent(out):   error code and error message
  if(eqns_data%err > 0)then; eqns_data%message=trim(eqns_data%message); ierr=-1; return; endif
  if(eqns_data%err < 0)then; eqns_data%message=trim(eqns_data%message); ierr=1; return; endif

  ! ierr = cudaDeviceSynchronize()

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
