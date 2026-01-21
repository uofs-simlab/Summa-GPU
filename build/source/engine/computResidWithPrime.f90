

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
  nGRU, &
                      ! input: model control
                      dt,                        & ! intent(in):  length of the time step (seconds)
                      nSnow,                     & ! intent(in):  number of snow layers
                      nSoil,                     & ! intent(in):  number of soil layers
                      nLayers,                   & ! intent(in):  total number of layers
                      enthalpyStateVec,          & ! intent(in):  flag if enthalpy is state variable
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
                      scalarCanopyCmTrial,       & ! intent(in):  Cm for vegetation canopy (J kg K-1)
                      mLayerCmTrial,             & ! intent(in):  Cm for each snow and soil layer (J kg K-1)
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
  ! --------------------------------------------------------------------------------------------------------------------------------
  use ieee_arithmetic
  use cudafor
  use device_data_types
  implicit none
  integer(i4b),intent(in) :: nGRU
  ! input: model control
  real(rkind),intent(in)          :: dt                        ! length of the time step (seconds)
  integer(i4b),intent(in),device         :: nSnow(:)                     ! number of snow layers
  integer(i4b),intent(in)         :: nSoil                     ! number of soil layers
  integer(i4b),intent(in),device         :: nLayers(:)                   ! total number of layers in the snow+soil domain
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
  real(rkind),intent(in),device          :: scalarCanopyCmTrial(:)       ! Cm for vegetation canopy (-)
  real(rkind),intent(in),device          :: mLayerCmTrial(:,:)          ! Cm for each snow and soil layer (-)
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
  ! real(rkind),dimension(nLayers)   :: mLayerVolFracHydPrime    ! vector of volumetric water content (-), either liquid water content or total water content
        type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! link to the necessary variables for the residual computations
  associate(&
    ! canopy and layer depth
    canopyDepth             => diag_data%scalarCanopyDepth      ,& ! intent(in): [dp]      canopy depth (m)
    mLayerDepth             => prog_data%mLayerDepth               ,& ! intent(in): [dp(:)]  depth of each layer in the snow-soil sub-domain (m)
    ! model fluxes (sink terms in the soil domain)
    mLayerTranspire_start         => flux_data%ixMLayerTranspire_start           ,& ! intent(in): [dp]     transpiration loss from each soil layer (m s-1)
    mLayerTranspire_end => flux_data%ixMLayerTranspire_end, &
    mLayerBaseflow_start          => flux_data%ixMLayerBaseflow_start            ,& ! intent(in): [dp(:)]  baseflow from each soil layer (m s-1)
    mLayerBaseflow_end => flux_data%ixMLayerBaseflow_end, &
    mLayerCompress          => diag_data%mLayerCompress            ,& ! intent(in): [dp(:)]  change in storage associated with compression of the soil matrix (-)
    ! number of state variables of a specific type
    ! nSnowSoilNrg            => indx_data%nSnowSoilNrg         ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
    ! nSnowSoilHyd            => indx_data%nSnowSoilHyd          ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
    ! nSoilOnlyHyd            => indx_data%nSoilOnlyHyd          ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
    ! model indices
    ixCasNrg                => indx_data%ixCasNrg              ,& ! intent(in): [i4b]    index of canopy air space energy state variable
    ixVegNrg                => indx_data%ixVegNrg              ,& ! intent(in): [i4b]    index of canopy energy state variable
    ixVegHyd                => indx_data%ixVegHyd              ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixAqWat                 => indx_data%ixAqWat               ,& ! intent(in): [i4b]    index of water storage in the aquifer
    ixSnowSoilNrg           => indx_data%ixSnowSoilNrg            ,& ! intent(in): [i4b(:)] indices for energy states in the snow+soil subdomain
    ixSnowSoilHyd           => indx_data%ixSnowSoilHyd            ,& ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
    ixSoilOnlyHyd           => indx_data%ixSoilOnlyHyd            ,& ! intent(in): [i4b(:)] indices for hydrology states in the soil subdomain
    ixStateType             => indx_data%ixStateType              ,& ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
    ixHydCanopy             => indx_data%ixHydCanopy              ,& ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
    ixHydType               => indx_data%ixHydType                ,& ! intent(in): [i4b(:)] named variables defining the type of hydrology states in snow+soil domain
    layerType               => indx_data%layerType                 & ! intent(in): [i4b(:)] named variables defining the type of layer in snow+soil domain
    ) ! association to necessary variables for the residual computations
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message="computResidWithPrime/"

    call computResidWithPrime_kernel<<<blocks,threads>>>(nGRU,&
dt,nSnow,nSoil,nLayers,enthalpyStateVec,&
sMul,fVec,&
scalarCanairTempPrime,scalarCanopyTempPrime,scalarCanopyWatPrime,mLayerTempPrime,scalarAquiferStoragePrime,&
scalarCanopyIcePrime,scalarCanopyLiqPrime,mLayerVolFracIcePrime,mLayerVolFracWatPrime,mLayerVolFracLiqPrime,&
scalarCanopyCmTrial,mLayerCmTrial,scalarCanairEnthalpyPrime,scalarCanopyEnthalpyPrime,mLayerEnthalpyPrime,&
rAdd,rVec,&
ixSnowSoilNrg,ixSnowSoilHyd,ixSoilOnlyHyd,&
ixStateType,ixHydType, layerType, &
ixVegNrg,ixCasNrg,ixVegHyd,ixAqWat,ixHydCanopy,&
canopyDepth,mLayerDepth,&
flux_data%data, &
mLayerTranspire_start,mLayerTranspire_end,mLayerBaseflow_start,mLayerBaseflow_end,mLayerCompress)
  end associate

end subroutine computResidWithPrime

attributes(global) subroutine computResidWithPrime_kernel(nGRU,&
dt,nSnow,nSoil,nLayers,enthalpyStateVec,&
sMul,fVec,&
scalarCanairTempPrime,scalarCanopyTempPrime,scalarCanopyWatPrime,mLayerTempPrime,scalarAquiferStoragePrime,&
scalarCanopyIcePrime,scalarCanopyLiqPrime,mLayerVolFracIcePrime,mLayerVolFracWatPrime,mLayerVolFracLiqPrime,&
scalarCanopyCmTrial,mLayerCmTrial,scalarCanairEnthalpyPrime,scalarCanopyEnthalpyPrime,mLayerEnthalpyPrime,&
rAdd,rVec,&
ixSnowSoilNrg,ixSnowSoilHyd,ixSoilOnlyHyd,&
ixStateType,ixHydType, layerType, &
ixVegNrg,ixCasNrg,ixVegHyd,ixAqWat,ixHydCanopy,&
canopyDepth,mLayerDepth,&
flux_data, &
mLayerTranspire_start,mLayerTranspire_end,mLayerBaseflow_start,mLayerBaseflow_end,mLayerCompress)
  integer(i4b),value :: nGRU
  ! input: model control
  real(rkind),intent(in),value          :: dt                        ! length of the time step (seconds)
  integer(i4b),intent(in)         :: nSnow(:)                     ! number of snow layers
  integer(i4b),intent(in),value         :: nSoil                     ! number of soil layers
  integer(i4b),intent(in)         :: nLayers(:)                   ! total number of layers in the snow+soil domain
  logical(lgt),intent(in),value         :: enthalpyStateVec          ! flag if enthalpy is state variable
  ! input: flux vectors
  real(qp),intent(in)             :: sMul(:,:)   ! NOTE: qp      ! state vector multiplier (used in the residual calculations)
  real(rkind),intent(in)          :: fVec(:,:)                   ! flux vector
  ! input: state variables (already disaggregated into scalars and vectors)
  real(rkind),intent(in)          :: scalarCanairTempPrime(:)     ! prime value for temperature of the canopy air space (K s-1)
  real(rkind),intent(in)          :: scalarCanopyTempPrime(:)     ! prime value for temperature of the vegetation canopy (K s-1)
  real(rkind),intent(in)          :: scalarCanopyWatPrime(:)      ! prime value for canopy total water content (kg m-2 s-1)
  real(rkind),intent(in)          :: mLayerTempPrime(:,:)        ! prime vector of temperature of each snow/soil layer (K s-1) content
  real(rkind),intent(in)          :: scalarAquiferStoragePrime(:) ! prime value of aquifer storage (m s-1)
  ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
  real(rkind),intent(in)          :: scalarCanopyIcePrime(:)      ! prime value for mass of ice on the vegetation canopy (kg m-2 s-1)
  real(rkind),intent(in)          :: scalarCanopyLiqPrime(:)      ! prime value for the liq on the vegetation canopy (kg m-2 s-1)
  real(rkind),intent(in)          :: mLayerVolFracIcePrime(:,:)  ! prime vector of volumetric fraction of ice (s-1)
  real(rkind),intent(in)          :: mLayerVolFracWatPrime(:,:)  ! prime vector of the volumetric water in each snow and soil layer (s-1)
  real(rkind),intent(in)          :: mLayerVolFracLiqPrime(:,:)  ! prime vector of the volumetric water in each snow and soil layer (s-1)
  ! input: enthalpy terms
  real(rkind),intent(in)          :: scalarCanopyCmTrial(:)       ! Cm for vegetation canopy (-)
  real(rkind),intent(in)          :: mLayerCmTrial(:,:)          ! Cm for each snow and soil layer (-)
  real(rkind),intent(in)          :: scalarCanairEnthalpyPrime(:) ! prime value for enthalpy of the canopy air space (W m-3)
  real(rkind),intent(in)          :: scalarCanopyEnthalpyPrime(:) ! prime value for enthalpy of the vegetation canopy (W m-3)
  real(rkind),intent(in)          :: mLayerEnthalpyPrime(:,:)    ! prime vector of enthalpy of each snow and soil layer (W m-3)
  ! output
  real(rkind),intent(out)         :: rAdd(:,:)                   ! additional (sink) terms on the RHS of the state equation
  real(qp),intent(out)            :: rVec(:,:)   ! NOTE: qp      ! residual vector

  integer(i4b),intent(in) :: ixSnowSoilNrg(:,:), ixSnowSoilHyd(:,:),ixSoilOnlyHyd(:,:)
  integer(i4b),intent(in) :: ixStateType(:,:), ixHydType(:,:), layerType(:,:)
  integer(i4b),intent(in) :: ixVegNrg(:),ixCasNrg(:),ixVegHyd(:),ixAqWat(:)
  integer(i4b),intent(in) :: ixHydCanopy(:)
  real(rkind),intent(in) :: canopyDepth(:), mLayerDepth(:,:)
  real(rkind),intent(in) :: flux_data(:,:)
  integer(i4b),intent(in),value :: mLayerTranspire_start,mLayerTranspire_end,mLayerBaseflow_start,mLayerBaseflow_end
  real(rkind),intent(in) :: mLayerCompress(:,:)

    integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return
  call computResidWithPrime_device(dt,nSnow(iGRU),nSoil,nLayers(iGRU),enthalpyStateVec,&
sMul(:,iGRU),fVec(:,iGRU),&
scalarCanairTempPrime(iGRU),scalarCanopyTempPrime(iGRU),scalarCanopyWatPrime(iGRU),mLayerTempPrime(:,iGRU),scalarAquiferStoragePrime(iGRU),&
scalarCanopyIcePrime(iGRU),scalarCanopyLiqPrime(iGRU),mLayerVolFracIcePrime(:,iGRU),mLayerVolFracWatPrime(:,iGRU),mLayerVolFracLiqPrime(:,iGRU),&
scalarCanopyCmTrial(iGRU),mLayerCmTrial(:,iGRU),scalarCanairEnthalpyPrime(iGRU),scalarCanopyEnthalpyPrime(iGRU),mLayerEnthalpyPrime(:,iGRU),&
rAdd(:,iGRU),rVec(:,iGRU),&
ixSnowSoilNrg(:,iGRU),ixSnowSoilHyd(:,iGRU),ixSoilOnlyHyd(:,iGRU),&
ixStateType(:,iGRU),ixHydType(:,iGRU), layerType(:,iGRU), &
ixVegNrg(iGRU),ixCasNrg(iGRU),ixVegHyd(iGRU),ixAqWat(iGRU),ixHydCanopy(iGRU),&
canopyDepth(iGRU),mLayerDepth(:,iGRU),&
flux_data(mLayerTranspire_start:mLayerTranspire_end,iGRU),flux_data(mLayerBaseflow_start:mLayerBaseflow_end,iGRU),mLayerCompress(:,iGRU))

end subroutine


attributes(device) subroutine computResidWithPrime_device(dt,nSnow,nSoil,nLayers,enthalpyStateVec,&
sMul,fVec,&
scalarCanairTempPrime,scalarCanopyTempPrime,scalarCanopyWatPrime,mLayerTempPrime,scalarAquiferStoragePrime,&
scalarCanopyIcePrime,scalarCanopyLiqPrime,mLayerVolFracIcePrime,mLayerVolFracWatPrime,mLayerVolFracLiqPrime,&
scalarCanopyCmTrial,mLayerCmTrial,scalarCanairEnthalpyPrime,scalarCanopyEnthalpyPrime,mLayerEnthalpyPrime,&
rAdd,rVec,&
ixSnowSoilNrg,ixSnowSoilHyd,ixSoilOnlyHyd,&
ixStateType,ixHydType, layerType, &
ixVegNrg,ixCasNrg,ixVegHyd,ixAqWat,ixHydCanopy,&
canopyDepth,mLayerDepth,&
mLayerTranspire,mLayerBaseflow,mLayerCompress)
implicit none

  ! input: model control
  real(rkind),intent(in)          :: dt                        ! length of the time step (seconds)
  integer(i4b),intent(in)         :: nSnow                     ! number of snow layers
  integer(i4b),intent(in)         :: nSoil                     ! number of soil layers
  integer(i4b),intent(in)         :: nLayers                   ! total number of layers in the snow+soil domain
  logical(lgt),intent(in)         :: enthalpyStateVec          ! flag if enthalpy is state variable
  ! input: flux vectors
  real(qp),intent(in)             :: sMul(:)   ! NOTE: qp      ! state vector multiplier (used in the residual calculations)
  real(rkind),intent(in)          :: fVec(:)                   ! flux vector
  ! input: state variables (already disaggregated into scalars and vectors)
  real(rkind),intent(in)          :: scalarCanairTempPrime     ! prime value for temperature of the canopy air space (K s-1)
  real(rkind),intent(in)          :: scalarCanopyTempPrime     ! prime value for temperature of the vegetation canopy (K s-1)
  real(rkind),intent(in)          :: scalarCanopyWatPrime      ! prime value for canopy total water content (kg m-2 s-1)
  real(rkind),intent(in)          :: mLayerTempPrime(:)        ! prime vector of temperature of each snow/soil layer (K s-1) content
  real(rkind),intent(in)          :: scalarAquiferStoragePrime ! prime value of aquifer storage (m s-1)
  ! input: diagnostic variables defining the liquid water and ice content (function of state variables)
  real(rkind),intent(in)          :: scalarCanopyIcePrime      ! prime value for mass of ice on the vegetation canopy (kg m-2 s-1)
  real(rkind),intent(in)          :: scalarCanopyLiqPrime      ! prime value for the liq on the vegetation canopy (kg m-2 s-1)
  real(rkind),intent(in)          :: mLayerVolFracIcePrime(:)  ! prime vector of volumetric fraction of ice (s-1)
  real(rkind),intent(in)          :: mLayerVolFracWatPrime(:)  ! prime vector of the volumetric water in each snow and soil layer (s-1)
  real(rkind),intent(in)          :: mLayerVolFracLiqPrime(:)  ! prime vector of the volumetric water in each snow and soil layer (s-1)
  ! input: enthalpy terms
  real(rkind),intent(in)          :: scalarCanopyCmTrial       ! Cm for vegetation canopy (-)
  real(rkind),intent(in)          :: mLayerCmTrial(:)          ! Cm for each snow and soil layer (-)
  real(rkind),intent(in)          :: scalarCanairEnthalpyPrime ! prime value for enthalpy of the canopy air space (W m-3)
  real(rkind),intent(in)          :: scalarCanopyEnthalpyPrime ! prime value for enthalpy of the vegetation canopy (W m-3)
  real(rkind),intent(in)          :: mLayerEnthalpyPrime(:)    ! prime vector of enthalpy of each snow and soil layer (W m-3)
  ! output
  real(rkind),intent(out)         :: rAdd(:)                   ! additional (sink) terms on the RHS of the state equation
  real(qp),intent(out)            :: rVec(:)   ! NOTE: qp      ! residual vector

  integer(i4b),intent(in) :: ixSnowSoilNrg(:), ixSnowSoilHyd(:),ixSoilOnlyHyd(:)
  integer(i4b),intent(in) :: ixStateType(:), ixHydType(:), layerType(:)
  integer(i4b),intent(in) :: ixVegNrg,ixCasNrg,ixVegHyd,ixAqWat
  integer(i4b),intent(in) :: ixHydCanopy
  real(rkind),intent(in) :: canopyDepth, mLayerDepth(:)
  real(rkind),intent(in) :: mLayerTranspire(:), mLayerBaseflow(:), mLayerCompress(:)

  integer(i4b)                     :: iLayer                   ! index of layer within the snow+soil domain
  integer(i4b),parameter           :: ixVegVolume=1            ! index of the desired vegetation control volumne (currently only one veg layer)
  real(rkind)                      :: scalarCanopyHydPrime     ! trial value for canopy water (kg m-2), either liquid water content or total water content
  real(rkind),dimension(nLayers)   :: mLayerVolFracHydPrime    ! vector of volumetric water content (-), either liquid water content or total water content

      ! ---
    ! * compute sink terms...
    ! -----------------------

    ! intialize additional terms on the RHS as zero
    rAdd(:) = 0._rkind

    ! add melt freeze terms only if not using enthalpy terms 
    ! NOTE: would need to use these if were using enthTemp terms
    if(.not.enthalpyStateVec)then
      ! compute energy associated with melt freeze for the vegetation canopy
      if(ixVegNrg/=integerMissing) rAdd(ixVegNrg) = rAdd(ixVegNrg) + LH_fus*scalarCanopyIcePrime/canopyDepth   ! energy associated with melt/freeze (J m-3)
 
      ! compute energy associated with melt/freeze for snow
      ! NOTE: allow expansion of ice during melt-freeze for snow; deny expansion of ice during melt-freeze for soil
      if(nLayers>0)then
        do iLayer=1,nLayers
          if (ixSnowSoilNrg(iLayer)/=integerMissing) then  ! (loop through non-missing energy state variables in the snow+soil domain)
          select case( layerType(iLayer) )
            case(iname_snow); rAdd( ixSnowSoilNrg(iLayer) ) = rAdd( ixSnowSoilNrg(iLayer) ) + LH_fus*iden_ice * mLayerVolFracIcePrime(iLayer)
            case(iname_soil); rAdd( ixSnowSoilNrg(iLayer) ) = rAdd( ixSnowSoilNrg(iLayer) ) + LH_fus*iden_water * mLayerVolFracIcePrime(iLayer)
          end select
        end if
        end do  ! looping through non-missing energy state variables in the snow+soil domain
      endif

    endif

    ! sink terms soil hydrology (-)
    ! NOTE 1: state variable is volumetric water content, so melt-freeze is not included
    ! NOTE 2: ground evaporation was already included in the flux at the upper boundary
    ! NOTE 3: rAdd(ixSnowOnlyWat)=0, and is defined in the initialization above
    ! NOTE 4: same sink terms for matric head and liquid matric potential
    if(nSoil>0)then
      do iLayer=1,nSoil
        if (ixSoilOnlyHyd(iLayer)/=integerMissing) then   ! (loop through non-missing hydrology state variables in the snow+soil domain)
       rAdd( ixSoilOnlyHyd(iLayer) ) = rAdd( ixSoilOnlyHyd(iLayer) ) + ( ( mLayerTranspire(iLayer) - mLayerBaseflow(iLayer) )/mLayerDepth(iLayer+nSnow) - mLayerCompress(iLayer) )*dt
        end if
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    endif

    ! ---
    ! * compute the residual vector...
    ! --------------------------------

    ! compute the residual vector for the vegetation canopy
    ! NOTE: sMul(ixVegHyd) = 1, but include as it converts all variables to quadruple precision
    ! --> energy balance
    if(enthalpyStateVec)then
      if(ixCasNrg/=integerMissing) rVec(ixCasNrg) = scalarCanairEnthalpyPrime - ( fVec(ixCasNrg)*dt + rAdd(ixCasNrg) )
      if(ixVegNrg/=integerMissing) rVec(ixVegNrg) = scalarCanopyEnthalpyPrime - ( fVec(ixVegNrg)*dt + rAdd(ixVegNrg) )
    else
      if(ixCasNrg/=integerMissing) rVec(ixCasNrg) = sMul(ixCasNrg) * scalarCanairTempPrime - ( fVec(ixCasNrg)*dt + rAdd(ixCasNrg) )
      if(ixVegNrg/=integerMissing) rVec(ixVegNrg) = sMul(ixVegNrg) * scalarCanopyTempPrime + scalarCanopyCmTrial * scalarCanopyWatPrime/canopyDepth &
                                                   - ( fVec(ixVegNrg)*dt + rAdd(ixVegNrg) )
    endif                                               
    ! --> mass balance
    if(ixVegHyd/=integerMissing)then
      scalarCanopyHydPrime = merge(scalarCanopyWatPrime, scalarCanopyLiqPrime, (ixStateType( ixHydCanopy )==iname_watCanopy) )
      rVec(ixVegHyd) = sMul(ixVegHyd)*scalarCanopyHydPrime - ( fVec(ixVegHyd)*dt + rAdd(ixVegHyd) )
    endif

    ! compute the residual vector for the snow and soil sub-domains for energy
    if(nLayers>0)then
      do iLayer=1,nLayers
        if (ixSnowSoilNrg(iLayer)/=integerMissing) then   ! (loop through non-missing energy state variables in the snow+soil domain)
        if(enthalpyStateVec)then
          rVec( ixSnowSoilNrg(iLayer) ) = mLayerEnthalpyPrime(iLayer) - ( fVec( ixSnowSoilNrg(iLayer) )*dt + rAdd( ixSnowSoilNrg(iLayer) ) )
        else
          rVec( ixSnowSoilNrg(iLayer) ) = sMul( ixSnowSoilNrg(iLayer) ) * mLayerTempPrime(iLayer) + mLayerCmTrial(iLayer) * mLayerVolFracWatPrime(iLayer) &
                                         - ( fVec( ixSnowSoilNrg(iLayer) )*dt + rAdd( ixSnowSoilNrg(iLayer) ) )
        endif
      end if
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    endif

    ! compute the residual vector for the snow and soil sub-domains for hydrology
    ! NOTE: residual depends on choice of state variable
    if(nLayers>0)then
      do iLayer=1,nLayers
        if (ixSnowSoilHyd(iLayer)/=integerMissing) then   ! (loop through non-missing hydrology state variables in the snow+soil domain)
        ! (get the correct state variable)
        mLayerVolFracHydPrime(iLayer) = merge(mLayerVolFracWatPrime(iLayer), mLayerVolFracLiqPrime(iLayer), (ixHydType(iLayer)==iname_watLayer .or. ixHydType(iLayer)==iname_matLayer) )
        ! (compute the residual)
        rVec( ixSnowSoilHyd(iLayer) ) = mLayerVolFracHydPrime(iLayer) - ( fVec( ixSnowSoilHyd(iLayer) )*dt + rAdd( ixSnowSoilHyd(iLayer) ) )
        end if
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    endif

    ! compute the residual vector for the aquifer
    if(ixAqWat/=integerMissing)  rVec(ixAqWat) = sMul(ixAqWat)*scalarAquiferStoragePrime - ( fVec(ixAqWat)*dt + rAdd(ixAqWat) )

    ! print the state variables if requested
    ! if(globalPrintFlag)then
    !   write(*,'(a)') 'In computResidWithPrime:'
    !   write(*,'(a,i4)') '  nSnow = ', nSnow
    !   write(*,'(a,i4)') '  nSoil = ', nSoil
    !   write(*,'(a,i4)') '  nLayers = ', nLayers
    !   write(*,'(a,f12.5)') '  dt = ', dt
    !   write(*,'(a,1x,100(e12.5,1x))') '  sMul = ', sMul(min(iJac1,size(sMul)):min(iJac2,size(sMul)))
    !   write(*,'(a,1x,100(e12.5,1x))') '  fVec = ', fVec(min(iJac1,size(fVec)):min(iJac2,size(fVec)))
    !   write(*,'(a,f12.5)') '  scalarCanairTempPrime = ', scalarCanairTempPrime
    !   write(*,'(a,f12.5)') '  scalarCanopyTempPrime = ', scalarCanopyTempPrime
    !   write(*,'(a,f12.5)') '  scalarCanopyWatPrime = ', scalarCanopyWatPrime
    !   write(*,'(a,1x,100(e12.5,1x))') '  mLayerTempPrime = ', mLayerTempPrime(min(iJac1,size(mLayerTempPrime)):min(iJac2,size(mLayerTempPrime)))
    !   write(*,'(a,f12.5)') '  scalarAquiferStoragePrime = ', scalarAquiferStoragePrime
    !   write(*,'(a,f12.5)') '  scalarCanopyIcePrime = ', scalarCanopyIcePrime
    !   write(*,'(a,f12.5)') '  scalarCanopyLiqPrime = ', scalarCanopyLiqPrime
    !   write(*,'(a,1x,100(e12.5,1x))') '  mLayerVolFracIcePrime = ', mLayerVolFracIcePrime(min(iJac1,size(mLayerVolFracIcePrime)):min(iJac2,size(mLayerVolFracIcePrime)))
    !   write(*,'(a,1x,100(e12.5,1x))') '  mLayerVolFracWatPrime = ', mLayerVolFracWatPrime(min(iJac1,size(mLayerVolFracWatPrime)):min(iJac2,size(mLayerVolFracWatPrime)))
    !   write(*,'(a,1x,100(e12.5,1x))') '  mLayerVolFracLiqPrime = ', mLayerVolFracLiqPrime(min(iJac1,size(mLayerVolFracLiqPrime)):min(iJac2,size(mLayerVolFracLiqPrime)))
    !   write(*,'(a,f12.5)') '  scalarCanopyCmTrial = ', scalarCanopyCmTrial
    !   write(*,'(a,1x,100(e12.5,1x))') '  mLayerCmTrial = ', mLayerCmTrial(min(iJac1,size(mLayerCmTrial)):min(iJac2,size(mLayerCmTrial)))
    !   write(*,'(a,f12.5)') '  scalarCanairEnthalpyPrime = ', scalarCanairEnthalpyPrime 
    !   write(*,'(a,f12.5)') '  scalarCanopyEnthalpyPrime = ', scalarCanopyEnthalpyPrime
    !   write(*,'(a,1x,100(e12.5,1x))') '  mLayerEnthalpyPrime = ', mLayerEnthalpyPrime(min(iJac1,size(mLayerEnthalpyPrime)):min(iJac2,size(mLayerEnthalpyPrime)))
    ! endif

    ! print result
    ! if(globalPrintFlag .or. any(ieee_is_nan(rVec)))then
    !   write(*,'(a,1x,100(e12.5,1x))') 'rVec = ', rVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
    !   write(*,'(a,1x,100(e12.5,1x))') 'fVec = ', fVec(min(iJac1,size(rVec)):min(iJac2,size(rVec)))
    ! endif
    ! if(any(ieee_is_nan(rVec)))then; 
    ! message=trim(message)//'NaN in residuals'; 
    ! err=20; 
    ! return; endif
  end subroutine

end module computResidWithPrime_module
