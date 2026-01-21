! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
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

module getVectorz_module

! data types
USE nrtype

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! domain types
USE globalData,only:iname_cas        ! named variables for canopy air space
USE globalData,only:iname_veg        ! named variables for vegetation canopy
USE globalData,only:iname_snow       ! named variables for snow
USE globalData,only:iname_soil       ! named variables for soil

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair  ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy  ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy  ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy  ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer   ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer   ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer   ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer   ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer   ! named variable defining the liquid matric potential state variable for soil layers
USE globalData,only:iname_watAquifer ! named variable defining the water storage in the aquifer

! metadata for information in the data structures
USE globalData,only:indx_meta       ! metadata for the variables in the index structure

! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength     ! data vector with variable length dimension (rkind)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookDIAG             ! named variables for structure elements
USE var_lookup,only:iLookPROG             ! named variables for structure elements
USE var_lookup,only:iLookDERIV            ! named variables for structure elements
USE var_lookup,only:iLookPARAM            ! named variables for structure elements
USE var_lookup,only:iLookINDEX            ! named variables for structure elements

! provide access to functions for the constitutive functions and derivatives
USE snow_utils_module,only:fracliquid     ! compute the fraction of liquid water (snow)
USE snow_utils_module,only:dFracLiq_dTk   ! differentiate the freezing curve w.r.t. temperature (snow)
USE soil_utils_module,only:dTheta_dTk     ! differentiate the freezing curve w.r.t. temperature (soil)
USE soil_utils_module,only:dTheta_dPsi    ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:dPsi_dTheta    ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:matricHead     ! compute the matric head based on volumetric water content
USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water

implicit none
private
public::popStateVec
public::getScaling
public::varExtract,varExtract_host
public::checkFeas,checkFeas_host
contains


! **********************************************************************************************************
! public subroutine popStateVec: populate model state vectors
! **********************************************************************************************************
subroutine popStateVec(&
  nGRU, &
                        ! input: data structures
                        nState,                  & ! intent(in):  number of desired state variables
                        enthalpyStateVec,        & ! intent(in):  flag if enthalpy is state variable          
                        prog_data,               & ! intent(in):  model prognostic variables for a local HRU
                        diag_data,               & ! intent(in):  model diagnostic variables for a local HRU
                        indx_data,               & ! intent(in):  indices defining model states and layers
                        ! output
                        stateVec,                & ! intent(out): model state vector
                        err,message)               ! intent(out): error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  use cudafor
  use device_data_types
  integer(i4b) :: nGRU
  ! input: data structures
  integer(i4b),intent(in)         :: nState                 ! number of desired state variables
  logical(lgt),intent(in)         :: enthalpyStateVec       ! flag if enthalpy is state variable
  type(prog_data_device),intent(in)    :: prog_data              ! prognostic variables for a local HRU
  type(diag_data_device),intent(in)    :: diag_data              ! diagnostic variables for a local HRU
  type(indx_data_device),intent(in)    :: indx_data              ! indices defining model states and layers
  ! output
  real(rkind),intent(out),device         :: stateVec(:,:)            ! model state vector (mixed units)
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
                type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  fixedLength: associate(&
    ! model states for the vegetation canopy
    scalarCanairTemp     => prog_data%scalarCanairTemp     ,& ! intent(in) : [dp]     temperature of the canopy air space (K)
    scalarCanairEnthalpy => prog_data%scalarCanairEnthalpy ,& ! intent(in) : [dp]     enthalpy of the canopy air space (J m-3)
    scalarCanopyTemp     => prog_data%scalarCanopyTemp     ,& ! intent(in) : [dp]     temperature of the vegetation canopy (K)
    scalarCanopyEnthalpy => prog_data%scalarCanopyEnthalpy ,& ! intent(in) : [dp]     enthalpy of the vegetation canopy (J m-3)
    scalarCanopyWat      => prog_data%scalarCanopyWat      ,& ! intent(in) : [dp]     mass of total water on the vegetation canopy (kg m-2)
    scalarCanopyLiq      => prog_data%scalarCanopyLiq      ,& ! intent(in) : [dp]     mass of liquid water on the vegetation canopy (kg m-2)
    ! model state variable vectors for the snow-soil layers
    mLayerTemp           => prog_data%mLayerTemp              ,& ! intent(in) : [dp(:)]  temperature of each snow/soil layer (K)
    mLayerEnthalpy       => prog_data%mLayerEnthalpy          ,& ! intent(in) : [dp(:)]  enthalpy of each snow+soil layer (J m-3)
    mLayerVolFracWat     => prog_data%mLayerVolFracWat        ,& ! intent(in) : [dp(:)]  volumetric fraction of total water (-)
    mLayerVolFracLiq     => prog_data%mLayerVolFracLiq        ,& ! intent(in) : [dp(:)]  volumetric fraction of liquid water (-)
    mLayerMatricHead     => prog_data%mLayerMatricHead        ,& ! intent(in) : [dp(:)]  matric head (m)
    mLayerMatricHeadLiq  => diag_data%mLayerMatricHeadLiq     ,& ! intent(in) : [dp(:)]  matric potential of liquid water (m)
    ! model state variables for the aquifer
    scalarAquiferStorage => prog_data%scalarAquiferStorage ,& ! intent(in) : [dp]     storage of water in the aquifer (m)
    ! indices defining specific model states
    ixCasNrg             => indx_data%ixCasNrg               ,& ! intent(in) : [i4b(:)] [length=1] index of canopy air space energy state variable
    ixVegNrg             => indx_data%ixVegNrg               ,& ! intent(in) : [i4b(:)] [length=1] index of canopy energy state variable
    ixVegHyd             => indx_data%ixVegHyd               ,& ! intent(in) : [i4b(:)] [length=1] index of canopy hydrology state variable (mass)
    ixAqWat              => indx_data%ixAqWat                ,& ! intent(in) : [i4b(:)] [length=1] index of aquifer storage state variable
    ! vector of energy and hydrology indices for the snow and soil domains
    ixSnowSoilNrg        => indx_data%ixSnowSoilNrg          ,& ! intent(in) : [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    ixSnowSoilHyd        => indx_data%ixSnowSoilHyd          ,& ! intent(in) : [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    ! nSnowSoilNrg         => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)       ,& ! intent(in) : [i4b]    number of energy state variables in the snow+soil domain
    ! nSnowSoilHyd         => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)       ,& ! intent(in) : [i4b]    number of hydrology state variables in the snow+soil domain
    ! type of model state variabless
    ixStateType_subset   => indx_data%ixStateType     ,& ! intent(in) : [i4b(:)] [state subset] type of desired model state variables
    ixHydType            => indx_data%ixHydType              ,& ! intent(in) : [i4b(:)] index of the type of hydrology states in snow+soil domain
    ! number of layers
    nSnow                => indx_data%nSnow               ,& ! intent(in) : [i4b]    number of snow layers
    nSoil                => indx_data%nSoil               ,& ! intent(in) : [i4b]    number of soil layers
    nLayers              => indx_data%nLayers_d              & ! intent(in) : [i4b]    total number of layers
    )  ! end association with variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='popStateVec/'

    call popStateVec_kernel<<<1,1>>>(nGRU,ixCasNrg,ixVegNrg,ixVegHyd,ixAqWat,&
  ixSnowSoilNrg, ixSnowSoilHyd, ixStateType_subset, ixHydType, &
  enthalpyStateVec, nLayers, nSnow, &
  scalarCanairEnthalpy, scalarCanairTemp, &
  scalarCanopyEnthalpy, scalarCanopyTemp, &
  scalarCanopyWat, scalarCanopyLiq, &
  mLayerEnthalpy, mLayerTemp, &
  mLayerVolFracWat, mLayerVolFracLiq, &
  mLayerMatricHead, mLayerMatricHeadLiq, &
  scalarAquiferStorage, &
  stateVec)


  end associate fixedLength      ! end association to variables in the data structure where vector length does not change
end subroutine popStateVec

attributes(global) subroutine popStateVec_kernel(nGRU,ixCasNrg,ixVegNrg,ixVegHyd,ixAqWat,&
  ixSnowSoilNrg, ixSnowSoilHyd, ixStateType_subset, ixHydType, &
  enthalpyStateVec, nLayers, nSnow, &
  scalarCanairEnthalpy, scalarCanairTemp, &
  scalarCanopyEnthalpy, scalarCanopyTemp, &
  scalarCanopyWat, scalarCanopyLiq, &
  mLayerEnthalpy, mLayerTemp, &
  mLayerVolFracWat, mLayerVolFracLiq, &
  mLayerMatricHead, mLayerMatricHeadLiq, &
  scalarAquiferStorage, &
  stateVec)
  integer(i4b),value :: nGRU
  integer(i4b),intent(in) :: ixCasNrg(:), ixVegNrg(:), ixVegHyd(:), ixAqWat(:)
  integer(i4b) :: ixSnowSoilNrg(:,:), ixSnowSoilHyd(:,:), ixStateType_subset(:,:), ixHydType(:,:)
  logical(lgt),intent(in),value :: enthalpyStateVec
  integer(i4b),intent(in) :: nLayers(:), nSnow(:)
  real(rkind),intent(in) :: scalarCanairEnthalpy(:), scalarCanairTemp(:)
  real(rkind),intent(in) :: scalarCanopyEnthalpy(:), scalarCanopyTemp(:)
  real(rkind),intent(in) :: scalarCanopyWat(:), scalarCanopyLiq(:)
  real(rkind),intent(in) :: mLayerEnthalpy(:,:), mLayerTemp(:,:)
  real(rkind),intent(in) :: mLayerVolFracWat(:,:), mLayerVolFracLiq(:,:), mLayerMatricHead(:,:), mLayerMatricHeadLiq(:,:)
  real(rkind),intent(in) :: scalarAquiferStorage(:)
  real(rkind),intent(inout) :: stateVec(:,:)


          integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

  do iGRU=1,nGRU
  call popStateVec_device(ixCasNrg(iGRU),ixVegNrg(iGRU),ixVegHyd(iGRU),ixAqWat(iGRU),&
  ixSnowSoilNrg(:,iGRU), ixSnowSoilHyd(:,iGRU), ixStateType_subset(:,iGRU), ixHydType(:,iGRU), &
  enthalpyStateVec, nLayers(iGRU), nSnow(iGRU), &
  scalarCanairEnthalpy(iGRU), scalarCanairTemp(iGRU), &
  scalarCanopyEnthalpy(iGRU), scalarCanopyTemp(iGRU), &
  scalarCanopyWat(iGRU), scalarCanopyLiq(iGRU), &
  mLayerEnthalpy(:,iGRU), mLayerTemp(:,iGRU), &
  mLayerVolFracWat(:,iGRU), mLayerVolFracLiq(:,iGRU), &
  mLayerMatricHead(:,iGRU), mLayerMatricHeadLiq(:,iGRU), &
  scalarAquiferStorage(iGRU), &
  stateVec(:,iGRU))
  end do
end subroutine


attributes(device) subroutine popStateVec_device(ixCasNrg,ixVegNrg,ixVegHyd,ixAqWat,&
  ixSnowSoilNrg, ixSnowSoilHyd, ixStateType_subset, ixHydType, &
  enthalpyStateVec, nLayers, nSnow, &
  scalarCanairEnthalpy, scalarCanairTemp, &
  scalarCanopyEnthalpy, scalarCanopyTemp, &
  scalarCanopyWat, scalarCanopyLiq, &
  mLayerEnthalpy, mLayerTemp, &
  mLayerVolFracWat, mLayerVolFracLiq, &
  mLayerMatricHead, mLayerMatricHeadLiq, &
  scalarAquiferStorage, &
  stateVec)
  integer(i4b),intent(in) :: ixCasNrg, ixVegNrg, ixVegHyd, ixAqWat
  integer(i4b) :: ixSnowSoilNrg(:), ixSnowSoilHyd(:), ixStateType_subset(:), ixHydType(:)
  logical(lgt),intent(in) :: enthalpyStateVec
  integer(i4b),intent(in) :: nLayers, nSnow
  real(rkind),intent(in) :: scalarCanairEnthalpy, scalarCanairTemp
  real(rkind),intent(in) :: scalarCanopyEnthalpy, scalarCanopyTemp
  real(rkind),intent(in) :: scalarCanopyWat, scalarCanopyLiq
  real(rkind),intent(in) :: mLayerEnthalpy(:), mLayerTemp(:)
  real(rkind),intent(in) :: mLayerVolFracWat(:), mLayerVolFracLiq(:), mLayerMatricHead(:), mLayerMatricHeadLiq(:)
  real(rkind),intent(in) :: scalarAquiferStorage
  real(rkind),intent(inout) :: stateVec(:)

  integer(i4b) :: iLayer, ixStateSubset

    stateVec = 0._rkind
    ! -----
    ! * initialize state vectors...
    ! -----------------------------

    ! build the state vector for the temperature of the canopy air space
    if (ixCasNrg/=integerMissing) then
      if(enthalpyStateVec)then
        stateVec( ixCasNrg ) = scalarCanairEnthalpy ! transfer canopy air enthalpy to the state vector
      else
        stateVec( ixCasNrg ) = scalarCanairTemp     ! transfer canopy air temperature to the state vector
      endif
    end if

    print*, scalarCanopyTemp, scalarCanopyEnthalpy
    ! build the state vector for the temperature of the vegetation canopy
    if (ixVegNrg/=integerMissing) then
      if(enthalpyStateVec)then
        stateVec( ixVegNrg ) = scalarCanopyEnthalpy ! transfer vegetation enthalpy to the state vector
      else
        stateVec( ixVegNrg )  = scalarCanopyTemp    ! transfer vegetation temperature to the state vector
      endif
    end if

    ! build the state vector for the water in the vegetation canopy
    if (ixVegHyd/=integerMissing) then
      select case(ixStateType_subset( ixVegHyd ))
        case(iname_watCanopy); stateVec( ixVegHyd ) = scalarCanopyWat ! transfer total canopy water to the state vector
        case(iname_liqCanopy); stateVec( ixVegHyd ) = scalarCanopyLiq ! transfer liquid canopy water to the state vector
      end select
    end if

    ! build the energy state vector for the snow and soil domain
    do iLayer=1,nLayers
      if (ixSnowSoilNrg(iLayer)/=integerMissing) then ! (loop through non-missing energy state variables in the snow+soil domain)
        ixStateSubset            = ixSnowSoilNrg(iLayer) ! index within the state vector
        if(enthalpyStateVec)then
          stateVec(ixStateSubset) = mLayerEnthalpy(iLayer) ! transfer enthalpy from a layer to the state vector
        else
          stateVec(ixStateSubset) = mLayerTemp(iLayer)     ! transfer temperature from a layer to the state vector
        endif
      end if
    end do  ! looping through non-missing energy state variables in the snow+soil domain

    ! build the hydrology state vector for the snow+soil domains
    ! NOTE: ixVolFracWat  and ixVolFracLiq can also include states in the soil domain, hence enable primary variable switching
    do iLayer=1,nLayers
      if (ixSnowSoilHyd(iLayer)/=integerMissing) then ! (loop through non-missing hydrology state variables in the snow+soil domain)
        ixStateSubset            = ixSnowSoilHyd(iLayer) ! index within the state vector
        select case( ixHydType(iLayer) )
          case(iname_watLayer); stateVec(ixStateSubset) = mLayerVolFracWat(iLayer)           ! total water state variable for snow+soil layers
          case(iname_liqLayer); stateVec(ixStateSubset) = mLayerVolFracLiq(iLayer)           ! liquid water state variable for snow+soil layers
          case(iname_matLayer); stateVec(ixStateSubset) = mLayerMatricHead(iLayer-nSnow)     ! total water matric potential variable for soil layers
          case(iname_lmpLayer); stateVec(ixStateSubset) = mLayerMatricHeadLiq(iLayer-nSnow)  ! liquid matric potential state variable for soil layers
        end select
      end if
    end do  ! looping through non-missing energy state variables in the snow+soil domain

    ! build the state vector for the aquifer storage
    if (ixAqWat/=integerMissing) then
      stateVec( ixAqWat )  = scalarAquiferStorage    ! transfer aquifer storage to the state vector
    end if

  end subroutine popStateVec_device

! **********************************************************************************************************
! public subroutine getScaling: get scale factors
! **********************************************************************************************************
subroutine getScaling(&
  nGRU, &
                      ! input: data structures
                      diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                      indx_data,               & ! intent(in):    indices defining model states and layers
                      ! output
                      fScale,                  & ! intent(out):   characteristic scale of the function evaluations (mixed units)
                      xScale,                  & ! intent(out):   variable scaling vector (mixed units)
                      sMul,                    & ! intent(out):   multiplier for state vector (used in the residual calculations)
                      dMat,                    & ! intent(out):   diagonal of the Jacobian matrix excluding fluxes, not depending on the state vector
                      err,message)               ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  USE nr_utility_module,only:arth                   ! get a sequence of numbers arth(start, incr, count)
  USE f2008funcs_module,only:findIndex              ! finds the index of the first value within a vector
  use device_data_types
  ! --------------------------------------------------------------------------------------------------------------------------------
  integer(i4b) :: nGRU
  ! input: data structures
  type(diag_data_device),intent(in)    :: diag_data              ! diagnostic variables for a local HRU
  type(indx_data_device),intent(in)    :: indx_data              ! indices defining model states and layers
  ! output: state vectors
  real(rkind),intent(out),device         :: fScale(:,:)              ! characteristic scale of the function evaluations (mixed units)
  real(rkind),intent(out),device         :: xScale(:,:)              ! variable scaling vector (mixed units)
  real(qp),intent(out),device            :: sMul(:,:)    ! NOTE: qp  ! multiplier for state vector (used in the residual calculations)
  real(rkind),intent(out),device         :: dMat(:,:)                ! diagonal of the Jacobian matrix excluding fluxes, not depending on the state vector
  ! output: error control
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! scaling parameters
  real(rkind),parameter           :: fScaleLiq=0.01_rkind      ! func eval: characteristic scale for volumetric liquid water content (-)
  real(rkind),parameter           :: fScaleMat=10._rkind       ! func eval: characteristic scale for matric head (m)
  real(rkind),parameter           :: fScaleNrg=1000000._rkind  ! func eval: characteristic scale for energy (J m-3)
  real(rkind),parameter           :: xScaleLiq=0.1_rkind       ! state var: characteristic scale for volumetric liquid water content (-)
  real(rkind),parameter           :: xScaleMat=10._rkind       ! state var: characteristic scale for matric head (m)
  real(rkind),parameter           :: xScaleTemp=1._rkind       ! state var: characteristic scale for temperature (K)
  ! state subsets
  integer(i4b)                    :: iLayer                 ! index of layer within the snow+soil domain
  integer(i4b)                    :: ixStateSubset          ! index within the state subset
              type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  fixedLength: associate(&
    ! model diagnostic variables
    canopyDepth         => diag_data%scalarCanopyDepth      ,& ! intent(in):  [dp]     canopy depth (m)
    volHeatCapVeg       => diag_data%scalarBulkVolHeatCapVeg,& ! intent(in) : [dp]     bulk volumetric heat capacity of vegetation (J m-3 K-1)
    mLayerVolHeatCap    => diag_data%mLayerVolHtCapBulk        ,& ! intent(in) : [dp(:)]  bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
    ! indices defining specific model states
    ixCasNrg            => indx_data%ixCasNrg                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy air space energy state variable
    ixVegNrg            => indx_data%ixVegNrg                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy energy state variable
    ixVegHyd            => indx_data%ixVegHyd                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy hydrology state variable (mass)
    ! vector of energy and hydrology indices for the snow and soil domains
    ixSnowSoilNrg       => indx_data%ixSnowSoilNrg            ,& ! intent(in) : [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    ixSnowSoilHyd       => indx_data%ixSnowSoilHyd            ,& ! intent(in) : [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    ! nSnowSoilNrg        => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in) : [i4b]    number of energy state variables in the snow+soil domain
    ! nSnowSoilHyd        => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in) : [i4b]    number of hydrology state variables in the snow+soil domain
    ! type of model state variabless
    ixStateType_subset  => indx_data%ixStateType       ,& ! intent(in) : [i4b(:)] [state subset] type of desired model state variables
    ! number of layers
    nSnow               => indx_data%nSnow                 ,& ! intent(in) : [i4b]    number of snow layers
    nSoil               => indx_data%nSoil                 ,& ! intent(in) : [i4b]    number of soil layers
    nLayers             => indx_data%nLayers_d                & ! intent(in) : [i4b]    total number of layers
    )  ! end association with variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='getScaling/'

    call getScaling_kernel<<<blocks,threads>>>(nGRU,ixStateType_subset,ixSnowSoilNrg,ixSnowSoilHyd,nLayers,&
  mLayerVolHeatCap, volHeatCapVeg,canopyDepth, &
  fScale,xScale,dMat,sMul)

  end associate fixedLength      ! end association to variables in the data structure where vector length does not change
end subroutine getScaling

attributes(global) subroutine getScaling_kernel(nGRU,ixStateType_subset,ixSnowSoilNrg,ixSnowSoilHyd,nLayers,&
  mLayerVolHeatCap, volHeatCapVeg,canopyDepth, &
  fScale,xScale,dMat,sMul)
  integer(i4b),value :: nGRU
  integer(i4b),intent(in) :: ixStateType_subset(:,:), ixSnowSoilNrg(:,:), ixSnowSoilHyd(:,:),nLayers(:)
  real(rkind),intent(in) :: mLayerVolHeatCap(:,:),volHeatCapVeg(:),canopyDepth(:)
  real(rkind),intent(inout) :: fScale(:,:), xScale(:,:), dMat(:,:), sMul(:,:)
  ! scaling parameters
  real(rkind),parameter           :: fScaleLiq=0.01_rkind      ! func eval: characteristic scale for volumetric liquid water content (-)
  real(rkind),parameter           :: fScaleMat=10._rkind       ! func eval: characteristic scale for matric head (m)
  real(rkind),parameter           :: fScaleNrg=1000000._rkind  ! func eval: characteristic scale for energy (J m-3)
  real(rkind),parameter           :: xScaleLiq=0.1_rkind       ! state var: characteristic scale for volumetric liquid water content (-)
  real(rkind),parameter           :: xScaleMat=10._rkind       ! state var: characteristic scale for matric head (m)
  real(rkind),parameter           :: xScaleTemp=1._rkind       ! state var: characteristic scale for temperature (K)


          integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

  call getScaling_device(ixStateType_subset(:,iGRU),ixSnowSoilNrg(:,iGRU),ixSnowSoilHyd(:,iGRU),nLayers(iGRU),&
  mLayerVolHeatCap(:,iGRU), volHeatCapVeg(iGRU),canopyDepth(iGRU),fScaleNrg,fScaleLiq, &
  fScale(:,iGRU),xScale(:,iGRU),dMat(:,iGRU),sMul(:,iGRU))
end subroutine

attributes(device) subroutine getScaling_device(ixStateType_subset,ixSnowSoilNrg,ixSnowSoilHyd,nLayers,&
  mLayerVolHeatCap, volHeatCapVeg,canopyDepth,fScaleNrg,fScaleLiq, &
  fScale,xScale,dMat,sMul)
  implicit none
  integer(i4b),intent(in) :: ixStateType_subset(:), ixSnowSoilNrg(:), ixSnowSoilHyd(:),nLayers
  real(rkind),intent(in) :: mLayerVolHeatCap(:),volHeatCapVeg,canopyDepth,fScaleNrg,fScaleLiq
  real(rkind),intent(inout) :: fScale(:), xScale(:), dMat(:), sMul(:)

  integer(i4b) :: iLayer, ixStateSubset
    ! -----
    ! * define scaling vectors...
    ! ---------------------------

    ! define the function and variable scaling factors for energy
    where(ixStateType_subset==iname_nrgCanair .or. ixStateType_subset==iname_nrgCanopy .or. ixStateType_subset==iname_nrgLayer)
      fScale = 1._rkind / fScaleNrg  ! 1/(J m-3)
      xScale = 1._rkind  ! K
    endwhere

    ! define the function and variable scaling factors for water on the vegetation canopy
    where(ixStateType_subset==iname_watCanopy .or. ixStateType_subset==iname_liqCanopy)
      fScale = 1._rkind / (fScaleLiq*canopyDepth*iden_water)  ! 1/(kg m-2)
      xScale = 1._rkind  ! (kg m-2)
    endwhere

    ! define the function and variable scaling factors for water in the snow+soil domain
    where(ixStateType_subset==iname_watLayer .or. ixStateType_subset==iname_liqLayer)
      fScale = 1._rkind / fScaleLiq  ! (-)
      xScale = 1._rkind  ! (-)
    end where

    ! define the function and variable scaling factors for water in the snow+soil domain
    where(ixStateType_subset==iname_matLayer .or. ixStateType_subset==iname_lmpLayer)
      fScale = 1._rkind / fScaleLiq  ! (-)
      xScale = 1._rkind  ! (m)
    end where

    ! define the function and variable scaling factors for water storage in the aquifer
    where(ixStateType_subset==iname_watAquifer)
      fScale = 1._rkind
      xScale = 1._rkind
    endwhere

    ! -----
    ! * define components of derivative matrices at start of time step (substep)...
    ! ------------------------------------------------------------------------------------------

    ! define the multiplier for the state vector for residual calculations (vegetation canopy)
    ! NOTE: Use the "where" statement to generalize to multiple canopy layers (currently one canopy layer)
    where(ixStateType_subset==iname_nrgCanair) sMul = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
    where(ixStateType_subset==iname_nrgCanopy) sMul = volHeatCapVeg     ! volumetric heat capacity of the vegetation (J m-3 K-1)
    where(ixStateType_subset==iname_watCanopy) sMul = 1._rkind          ! nothing else on the left hand side
    where(ixStateType_subset==iname_liqCanopy) sMul = 1._rkind          ! nothing else on the left hand side

    ! compute terms in the Jacobian for vegetation (excluding fluxes)
    ! NOTE: This is computed outside the iteration loop because it does not depend on state variables
    ! NOTE: Energy for vegetation is computed *within* the iteration loop as it includes phase change
    ! NOTE: Use the "where" statement to generalize to multiple canopy layers (currently one canopy layer)
    where(ixStateType_subset==iname_nrgCanair) dMat = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
    where(ixStateType_subset==iname_nrgCanopy) dMat = realMissing       ! populated within the iteration loop
    where(ixStateType_subset==iname_watCanopy) dMat = 1._rkind          ! nothing else on the left hand side
    where(ixStateType_subset==iname_liqCanopy) dMat = 1._rkind          ! nothing else on the left hand side

    ! define the energy multiplier and diagonal elements for the state vector for residual calculations (snow-soil domain)
    do iLayer=1,nLayers
      if (ixSnowSoilNrg(iLayer)/=integerMissing) then   ! (loop through non-missing energy state variables in the snow+soil domain)
        ixStateSubset        = ixSnowSoilNrg(iLayer)      ! index within the state vector
        sMul(ixStateSubset)  = mLayerVolHeatCap(iLayer)   ! transfer volumetric heat capacity to the state multiplier
        dMat(ixStateSubset)  = realMissing                ! diagonal element populated within the iteration loop
      end if
    end do  ! looping through non-missing energy state variables in the snow+soil domain
    
    ! define the hydrology multiplier and diagonal elements for the state vector for residual calculations (snow-soil domain)
    do iLayer=1,nLayers
      if (ixSnowSoilHyd(iLayer)/=integerMissing) then   ! (loop through non-missing energy state variables in the snow+soil domain)
        ixStateSubset        = ixSnowSoilHyd(iLayer)      ! index within the state vector
        sMul(ixStateSubset)  = 1._rkind                   ! state multiplier = 1 (nothing else on the left-hand-side)
        dMat(ixStateSubset)  = 1._rkind                   ! diagonal element = 1 (nothing else on the left-hand-side)
      end if
    end do  ! looping through non-missing energy state variables in the snow+soil domain
    
    ! define the scaling factor and diagonal elements for the aquifer
    where(ixStateType_subset==iname_watAquifer)
      sMul = 1._rkind
      dMat = 1._rkind
    endwhere
  end subroutine


! **********************************************************************************************************
! public subroutine checkFeas: check feasibility of the state vector
! **********************************************************************************************************
subroutine checkFeas_host(&
                      ! input
                      stateVec,                                  & ! intent(in):    model state vector (mixed units)
                      mpar_data,                                 & ! intent(in):    model parameters
                      prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                      indx_data,                                 & ! intent(in):    indices defining model states and layers
                      enthalpyStateVec,                          & ! intent(in):    flag if enthalpy is state variable
                      ! output: feasibility
                      feasible,                                  & ! intent(inout):   flag to denote the feasibility of the solution
                    ! output: error control
                      err,message)                                 ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input
  real(rkind),intent(in)          :: stateVec(:)               ! model state vector (mixed units)
  type(var_dlength),intent(in)    :: mpar_data                 ! model parameters
  type(var_dlength),intent(in)    :: prog_data                 ! prognostic variables for a local HRU
  type(var_ilength),intent(in)    :: indx_data                 ! indices defining model states and layers
  logical(lgt),intent(in)         :: enthalpyStateVec          ! flag if enthalpy is state variable
  ! output: feasibility
  logical(lgt),intent(inout)      :: feasible                  ! flag to denote the feasibility of the solution
  ! output: error control
  integer(i4b),intent(out)        :: err                       ! error code
  character(*),intent(out)        :: message                   ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  integer(i4b)                    :: iLayer                    ! index of layer within the snow+soil domain
  real(rkind)                     :: xMin,xMax                 ! minimum and maximum values for water content
  real(rkind),parameter           :: canopyTempMax=500._rkind  ! expected maximum value for the canopy temperature (K)
  logical(lgt),parameter          :: printFlag=.false.         ! flag to denote if we print infeasibilities
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  associate(&
    ! soil parameters
    theta_sat               => mpar_data%var(iLookPARAM%theta_sat)%dat       ,&  ! intent(in): [dp(:)]  soil porosity (-)
    theta_res               => mpar_data%var(iLookPARAM%theta_res)%dat       ,&  ! intent(in): [dp(:)]  residual volumetric water content (-)
    ! model diagnostic variables from the previous solution
    mLayerVolFracIce        => prog_data%var(iLookPROG%mLayerVolFracIce)%dat ,& ! intent(in):  [dp(:)]  volumetric fraction of ice (-)
    ! number of model layers, and layer type
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)        ,& ! intent(in):  [i4b]    total number of snow layers
    nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)        ,& ! intent(in):  [i4b]    total number of soil layers
    nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)      ,& ! intent(in):  [i4b]    total number of snow and soil layers
    ! indices defining model states and layers
    ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)     ,& ! intent(in):  [i4b]    index of canopy air space energy state variable
    ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)     ,& ! intent(in):  [i4b]    index of canopy energy state variable
    ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)     ,& ! intent(in):  [i4b]    index of canopy hydrology state variable (mass)
    ixSnowOnlyNrg           => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat   ,&  ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
    ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat   ,&  ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
    ixStateType             => indx_data%var(iLookINDEX%ixStateType)%dat     ,&  ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
    ixHydCanopy             => indx_data%var(iLookINDEX%ixHydCanopy)%dat     ,&  ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
    ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat       ,&  ! intent(in): [i4b(:)] index of the type of hydrology states in snow+soil domain
    layerType               => indx_data%var(iLookINDEX%layerType)%dat        &  ! intent(in): [i4b(:)] layer type (iname_soil or iname_snow)
    )! association with variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message="checkFeas/"

    !  NOTE: we will not print infeasibilities since it does not indicate a failure, just a need to iterate until maxiter
    feasible=.true.
    ! check that the canopy air space temperature is reasonable
    if(ixCasNrg/=integerMissing)then
      if(stateVec(ixCasNrg) > canopyTempMax .and. .not.enthalpyStateVec)then 
        feasible=.false.
        message=trim(message)//'canopy air space temp high/'
        if(printFlag) write(*,'(a,1x,L1,1x,10(f20.10,1x))') 'feasible, max, stateVec( ixCasNrg )', feasible, canopyTempMax, stateVec(ixCasNrg)
      endif
    endif

    ! check that the canopy air space temperature is reasonable
    if(ixVegNrg/=integerMissing)then
      if(stateVec(ixVegNrg) > canopyTempMax .and. .not.enthalpyStateVec)then
        feasible=.false.
        message=trim(message)//'canopy temp high/'
        if(printFlag) write(*,'(a,1x,L1,1x,10(f20.10,1x))') 'feasible, max, stateVec( ixVegNrg )', feasible, canopyTempMax, stateVec(ixVegNrg)
      endif
    endif

    ! check canopy liquid water is not negative
    if(ixVegHyd/=integerMissing)then
      if(stateVec(ixVegHyd) < 0._rkind)then 
        feasible=.false.
        message=trim(message)//'canopy liq water neg/'
        if(printFlag) write(*,'(a,1x,L1,1x,10(f20.10,1x))') 'feasible, min, stateVec( ixVegHyd )', feasible, 0._rkind, stateVec(ixVegHyd)
      endif
    endif

    ! check snow temperature is below freezing
    if(count(ixSnowOnlyNrg/=integerMissing)>0)then
      if(any(stateVec( pack(ixSnowOnlyNrg,ixSnowOnlyNrg/=integerMissing) ) > Tfreeze) .and. .not.enthalpyStateVec)then
        feasible=.false.
        message=trim(message)//'snow temp high/'
        if(printFlag)then 
          do iLayer=1,nSnow
            if(stateVec(ixSnowOnlyNrg(iLayer)) > Tfreeze) write(*,'(a,1x,i4,1x,L1,1x,10(f20.10,1x))') 'iLayer, feasible, max, stateVec( ixSnowOnlyNrg(iLayer) )', iLayer, feasible, Tfreeze, stateVec( ixSnowOnlyNrg(iLayer) )
          enddo
        endif
      endif
    endif

    ! loop through non-missing hydrology state variables in the snow+soil domain
    do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)

      ! check the minimum and maximum water constraints
      if(ixHydType(iLayer)==iname_watLayer .or. ixHydType(iLayer)==iname_liqLayer)then

        ! --> minimum
        if (layerType(iLayer) == iname_soil) then
          xMin = theta_res(iLayer-nSnow)
        else
          xMin = 0._rkind
        endif

        ! --> maximum
        select case( layerType(iLayer) )
          case(iname_snow); xMax = merge(1._rkind, 1._rkind - mLayerVolFracIce(iLayer), ixHydType(iLayer)==iname_watLayer)
          case(iname_soil); xMax = merge(theta_sat(iLayer-nSnow), theta_sat(iLayer-nSnow) - mLayerVolFracIce(iLayer), ixHydType(iLayer)==iname_watLayer)
        end select

        ! --> check
        if(stateVec( ixSnowSoilHyd(iLayer) ) < xMin .or. stateVec( ixSnowSoilHyd(iLayer) ) > xMax)then 
          feasible=.false.
          message=trim(message)//'layer water out of bounds/'
          if(printFlag)then 
            if(stateVec( ixSnowSoilHyd(iLayer) ) < xMin .or. stateVec( ixSnowSoilHyd(iLayer) ) > xMax) &
                write(*,'(a,1x,i4,1x,L1,1x,10(f20.10,1x))') 'iLayer, feasible, stateVec( ixSnowSoilHyd(iLayer) ), xMin, xMax = ', iLayer, feasible, stateVec( ixSnowSoilHyd(iLayer) ), xMin, xMax
          endif
        endif
      endif  ! if water states

    end do  ! loop through non-missing hydrology state variables in the snow+soil domain

  end associate    ! end association to variables in the data structure
end subroutine checkFeas_host

subroutine checkFeas(&
  nGRU, &
                      ! input
                      stateVec,                                  & ! intent(in):    model state vector (mixed units)
                      mpar_data,                                 & ! intent(in):    model parameters
                      prog_data,                                 & ! intent(in):    model prognostic variables for a local HRU
                      indx_data,                                 & ! intent(in):    indices defining model states and layers
                      enthalpyStateVec,                          & ! intent(in):    flag if enthalpy is state variable
                      ! output: feasibility
                      feasible,                                  & ! intent(inout):   flag to denote the feasibility of the solution
                    ! output: error control
                      err,message)                                 ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  use cudafor
  use device_data_types
  implicit none
  ! input
  integer(i4b),intent(in) :: nGRU
  real(rkind),intent(in),device          :: stateVec(:,:)               ! model state vector (mixed units)
  type(mpar_data_device),intent(in)    :: mpar_data                 ! model parameters
  type(prog_data_device),intent(in)    :: prog_data                 ! prognostic variables for a local HRU
  type(indx_data_device),intent(in)    :: indx_data                 ! indices defining model states and layers
  logical(lgt),intent(in)         :: enthalpyStateVec          ! flag if enthalpy is state variable
  ! output: feasibility
  logical(lgt),intent(inout)      :: feasible                  ! flag to denote the feasibility of the solution
  ! output: error control
  integer(i4b),intent(out)        :: err                       ! error code
  character(*),intent(out)        :: message                   ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  integer(i4b)                    :: iLayer                    ! index of layer within the snow+soil domain
  real(rkind)                     :: xMin,xMax                 ! minimum and maximum values for water content
  real(rkind),parameter           :: canopyTempMax=500._rkind  ! expected maximum value for the canopy temperature (K)
  logical(lgt),parameter          :: printFlag=.false.         ! flag to denote if we print infeasibilities
  logical(lgt),device :: feasible_device
            type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  associate(&
    ! soil parameters
    theta_sat               => mpar_data%theta_sat_       ,&  ! intent(in): [dp(:)]  soil porosity (-)
    theta_res               => mpar_data%theta_res_       ,&  ! intent(in): [dp(:)]  residual volumetric water content (-)
    ! model diagnostic variables from the previous solution
    mLayerVolFracIce        => prog_data%mLayerVolFracIce ,& ! intent(in):  [dp(:)]  volumetric fraction of ice (-)
    ! number of model layers, and layer type
    nSnow                   => indx_data%nSnow        ,& ! intent(in):  [i4b]    total number of snow layers
    nSoil                   => indx_data%nSoil        ,& ! intent(in):  [i4b]    total number of soil layers
    nLayers                 => indx_data%nLayers_d      ,& ! intent(in):  [i4b]    total number of snow and soil layers
    ! indices defining model states and layers
    ixCasNrg                => indx_data%ixCasNrg     ,& ! intent(in):  [i4b]    index of canopy air space energy state variable
    ixVegNrg                => indx_data%ixVegNrg     ,& ! intent(in):  [i4b]    index of canopy energy state variable
    ixVegHyd                => indx_data%ixVegHyd     ,& ! intent(in):  [i4b]    index of canopy hydrology state variable (mass)
    ixSnowOnlyNrg           => indx_data%ixSnowOnlyNrg   ,&  ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
    ixSnowSoilHyd           => indx_data%ixSnowSoilHyd   ,&  ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
    ixStateType             => indx_data%ixStateType     ,&  ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
    ixHydCanopy             => indx_data%ixHydCanopy     ,&  ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
    ixHydType               => indx_data%ixHydType       ,&  ! intent(in): [i4b(:)] index of the type of hydrology states in snow+soil domain
    layerType               => indx_data%layerType        &  ! intent(in): [i4b(:)] layer type (iname_soil or iname_snow)
    )! association with variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message="checkFeas/"

    !  NOTE: we will not print infeasibilities since it does not indicate a failure, just a need to iterate until maxiter
    feasible_device=.true.

    call checkFeas_kernel<<<blocks,threads>>>(nGRU,&
    stateVec,&
  ixCasNrg, ixVegNrg, ixVegHyd,&
  ixSnowOnlyNrg, ixSnowSoilHyd,&
  ixHydType,layerType,&
  nSnow,nLayers,&
  theta_res,theta_sat,mLayerVolFracIce,&
  enthalpyStateVec,feasible_device)
  feasible = feasible_device
  end associate    ! end association to variables in the data structure
end subroutine checkFeas

attributes(global) subroutine checkFeas_kernel(nGRU,&
  stateVec,&
  ixCasNrg, ixVegNrg, ixVegHyd,&
  ixSnowOnlyNrg, ixSnowSoilHyd,&
  ixHydType,layerType,&
  nSnow,nLayers,&
  theta_res,theta_sat,mLayerVolFracIce,&
  enthalpyStateVec,feasible)
  integer(i4b),intent(in),value :: nGRU
    real(rkind),intent(in) :: stateVec(:,:)
  integer(i4b),intent(in) :: ixCasNrg(:), ixVegNrg(:), ixVegHyd(:)
  integer(i4b),intent(in) :: ixSnowOnlyNrg(:,:), ixSnowSoilHyd(:,:)
  integer(i4b),intent(in) :: ixHydType(:,:), layerType(:,:)
  integer(i4b),intent(in) :: nSnow(:), nLayers(:)
  real(rkind),intent(in) :: theta_res(:,:), theta_sat(:,:), mLayerVolFracIce(:,:)
  logical(lgt),intent(in),value :: enthalpyStateVec
  logical(lgt),intent(inout) :: feasible

        integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return
  call checkFeas_device(stateVec(:,iGRU),&
  ixCasNrg(iGRU), ixVegNrg(iGRU), ixVegHyd(iGRU),&
  ixSnowOnlyNrg(:,iGRU), ixSnowSoilHyd(:,iGRU),&
  ixHydType(:,iGRU),layerType(:,iGRU),&
  nSnow(iGRU),nLayers(iGRU),&
  theta_res(:,iGRU),theta_sat(:,iGRU),mLayerVolFracIce(:,iGRU),&
  enthalpyStateVec,feasible)
end subroutine

attributes(device) subroutine checkFeas_device(stateVec,&
  ixCasNrg, ixVegNrg, ixVegHyd,&
  ixSnowOnlyNrg, ixSnowSoilHyd,&
  ixHydType,layerType,&
  nSnow,nLayers,&
  theta_res,theta_sat,mLayerVolFracIce,&
  enthalpyStateVec,feasible)
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  implicit none
  real(rkind),intent(in) :: stateVec(:)
  integer(i4b),intent(in) :: ixCasNrg, ixVegNrg, ixVegHyd
  integer(i4b),intent(in) :: ixSnowOnlyNrg(:), ixSnowSoilHyd(:)
  integer(i4b),intent(in) :: ixHydType(:), layerType(:)
  integer(i4b),intent(in) :: nSnow, nLayers
  real(rkind),intent(in) :: theta_res(:), theta_sat(:), mLayerVolFracIce(:)
  logical(lgt),intent(in) :: enthalpyStateVec
  logical(lgt),intent(inout) :: feasible

  integer(i4b)                    :: iLayer                    ! index of layer within the snow+soil domain
  real(rkind)                     :: xMin,xMax                 ! minimum and maximum values for water content
  real(rkind),parameter           :: canopyTempMax=500._rkind  ! expected maximum value for the canopy temperature (K)
  logical(lgt),parameter          :: printFlag=.false.         ! flag to denote if we print infeasibilities
  ! --------------------------------------------------------------------------------------------------------------------------------

    !  NOTE: we will not print infeasibilities since it does not indicate a failure, just a need to iterate until maxiter
    ! feasible=.true.
    ! check that the canopy air space temperature is reasonable
    if(ixCasNrg/=integerMissing)then
      if(stateVec(ixCasNrg) > canopyTempMax .and. .not.enthalpyStateVec)then 
        feasible=.false.
        ! if(printFlag) write(*,'(a,1x,L1,1x,10(f20.10,1x))') 'feasible, max, stateVec( ixCasNrg )', feasible, canopyTempMax, stateVec(ixCasNrg)
      endif
    endif

    ! check that the canopy air space temperature is reasonable
    if(ixVegNrg/=integerMissing)then
      if(stateVec(ixVegNrg) > canopyTempMax .and. .not.enthalpyStateVec)then
        feasible=.false.
        ! if(printFlag) write(*,'(a,1x,L1,1x,10(f20.10,1x))') 'feasible, max, stateVec( ixVegNrg )', feasible, canopyTempMax, stateVec(ixVegNrg)
      endif
    endif

    ! check canopy liquid water is not negative
    if(ixVegHyd/=integerMissing)then
      if(stateVec(ixVegHyd) < 0._rkind)then 
        feasible=.false.
        ! if(printFlag) write(*,'(a,1x,L1,1x,10(f20.10,1x))') 'feasible, min, stateVec( ixVegHyd )', feasible, 0._rkind, stateVec(ixVegHyd)
      endif
    endif

    ! check snow temperature is below freezing
      do iLayer=1,nSnow
        if (ixSnowOnlyNrg(iLayer)/=integerMissing .and. stateVec(iLayer) > Tfreeze .and. .not.enthalpyStateVec) then
          feasible=.false.
        end if
      end do

    ! loop through non-missing hydrology state variables in the snow+soil domain
    do iLayer=1,nLayers
      if (ixSnowSoilHyd(iLayer)/=integerMissing) then

      ! check the minimum and maximum water constraints
      if(ixHydType(iLayer)==iname_watLayer .or. ixHydType(iLayer)==iname_liqLayer)then

        ! --> minimum
        if (layerType(iLayer) == iname_soil) then
          xMin = theta_res(iLayer-nSnow)
        else
          xMin = 0._rkind
        endif

        ! --> maximum
        select case( layerType(iLayer) )
          case(iname_snow); xMax = merge(1._rkind, 1._rkind - mLayerVolFracIce(iLayer), ixHydType(iLayer)==iname_watLayer)
          case(iname_soil); xMax = merge(theta_sat(iLayer-nSnow), theta_sat(iLayer-nSnow) - mLayerVolFracIce(iLayer), ixHydType(iLayer)==iname_watLayer)
        end select

        ! --> check
        if(stateVec( ixSnowSoilHyd(iLayer) ) < xMin .or. stateVec( ixSnowSoilHyd(iLayer) ) > xMax)then 
          feasible=.false.
          ! if(printFlag)then 
            ! if(stateVec( ixSnowSoilHyd(iLayer) ) < xMin .or. stateVec( ixSnowSoilHyd(iLayer) ) > xMax) &
                ! write(*,'(a,1x,i4,1x,L1,1x,10(f20.10,1x))') 'iLayer, feasible, stateVec( ixSnowSoilHyd(iLayer) ), xMin, xMax = ', iLayer, feasible, stateVec( ixSnowSoilHyd(iLayer) ), xMin, xMax
          ! endif
        endif
      endif  ! if water states
    end if

    end do  ! loop through non-missing hydrology state variables in the snow+soil domain

end subroutine checkFeas_device


! **********************************************************************************************************
! public subroutine varExtract: extract variables from the state vector and compute diagnostic variables
!  This routine does not initialize any of the variables
! **********************************************************************************************************
subroutine varExtract_host(&
                       ! input
                       stateVec,                                  & ! intent(in):    model state vector (mixed units)
                       indx_data,                                 & ! intent(in):    indices defining model states and layers
                       ! output: variables for the vegetation canopy
                       scalarCanairNrgTrial,                      & ! intent(inout):   trial value of canopy air energy, temperature (K) or enthalpy (J m-3)
                       scalarCanopyNrgTrial,                      & ! intent(inout):   trial value of canopy energy, temperature (K) or enthalpy (J m-3)
                       scalarCanopyWatTrial,                      & ! intent(inout):   trial value of canopy total water (kg m-2)
                       scalarCanopyLiqTrial,                      & ! intent(inout):   trial value of canopy liquid water (kg m-2)
                       ! output: variables for the snow-soil domain
                       mLayerNrgTrial,                            & ! intent(inout):   trial vector of layer energy, temperature (K) or enthalpy (J m-3)
                       mLayerVolFracWatTrial,                     & ! intent(inout):   trial vector of volumetric total water content (-)
                       mLayerVolFracLiqTrial,                     & ! intent(inout):   trial vector of volumetric liquid water content (-)
                       mLayerMatricHeadTrial,                     & ! intent(inout):   trial vector of total water matric potential (m)
                       mLayerMatricHeadLiqTrial,                  & ! intent(inout):   trial vector of liquid water matric potential (m)
                       ! output: variables for the aquifer
                       scalarAquiferStorageTrial,                 & ! intent(out):   trial value of storage of water in the aquifer (m)
                       ! output: error control
                       err,message)                                 ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  implicit none
  ! input
  real(rkind),intent(in)             :: stateVec(:)                     ! model state vector (mixed units)
  type(var_ilength),intent(in)       :: indx_data                       ! indices defining model states and layers
  ! output: variables for the vegetation canopy
  real(rkind),intent(inout)          :: scalarCanairNrgTrial            ! trial value of canopy air energy, temperature (K) or enthalpy (J m-3)
  real(rkind),intent(inout)          :: scalarCanopyNrgTrial            ! trial value of canopy energy, temperature (K) or enthalpy (J m-3)
  real(rkind),intent(inout)          :: scalarCanopyWatTrial            ! trial value of canopy total water (kg m-2)
  real(rkind),intent(inout)          :: scalarCanopyLiqTrial            ! trial value of canopy liquid water (kg m-2)
  ! output: variables for the snow-soil domain
  real(rkind),intent(inout)          :: mLayerNrgTrial(:)               ! trial vector of layer energy, temperature (K) or enthalpy (J m-3)
  real(rkind),intent(inout)          :: mLayerVolFracWatTrial(:)        ! trial vector of volumetric total water content (-)
  real(rkind),intent(inout)          :: mLayerVolFracLiqTrial(:)        ! trial vector of volumetric liquid water content (-)
  real(rkind),intent(inout)          :: mLayerMatricHeadTrial(:)        ! trial vector of total water matric potential (m)
  real(rkind),intent(inout)          :: mLayerMatricHeadLiqTrial(:)     ! trial vector of liquid water matric potential (m)
  ! output: variables for the aquifer
  real(rkind),intent(inout)          :: scalarAquiferStorageTrial       ! trial value of storage of water in the aquifer (m)
  ! output: error control
  integer(i4b),intent(out)           :: err                             ! error code
  character(*),intent(out)           :: message                         ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  integer(i4b)                       :: iLayer                          ! index of layer within the snow+soil domain
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  associate(&
    ! number of model layers, and layer type
    nSnow                   => indx_data%var(iLookINDEX%nSnow)%dat(1)                 ,& ! intent(in):  [i4b]    total number of snow layers
    nSoil                   => indx_data%var(iLookINDEX%nSoil)%dat(1)                 ,& ! intent(in):  [i4b]    total number of soil layers
    nLayers                 => indx_data%var(iLookINDEX%nLayers)%dat(1)               ,& ! intent(in):  [i4b]    total number of snow and soil layers
    ! indices defining model states and layers
    ixCasNrg                => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy air space energy state variable
    ixVegNrg                => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy energy state variable
    ixVegHyd                => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)              ,& ! intent(in):  [i4b]    index of canopy hydrology state variable (mass)
    ixAqWat                 => indx_data%var(iLookINDEX%ixAqWat)%dat(1)               ,& ! intent(in):  [i4b]    index of the squifer storage state variable
    ixSnowSoilNrg           => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat            ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for energy states in the snow+soil subdomain
    ixSnowSoilHyd           => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat            ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for hydrology states in the snow+soil subdomain
    nSnowSoilNrg            => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)         ,& ! intent(in):  [i4b]    number of energy state variables in the snow+soil domain
    nSnowSoilHyd            => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)         ,& ! intent(in):  [i4b]    number of hydrology variables in the snow+soil domain
    ! indices defining type of model state variables
    ixStateType_subset      => indx_data%var(iLookINDEX%ixStateType_subset)%dat       ,& ! intent(in):  [i4b(:)] [state subset] type of desired model state variables
    ixHydType               => indx_data%var(iLookINDEX%ixHydType)%dat                 & ! intent(in):  [i4b(:)] index of the type of hydrology states in snow+soil domain
    )! association with variables in the data structures

    ! --------------------------------------------------------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message='varExtract/'

    ! *** extract state variables for the vegetation canopy

    ! check if computing the vegetation flux
    if(ixCasNrg/=integerMissing .or. ixVegNrg/=integerMissing .or. ixVegHyd/=integerMissing)then

      ! extract temperature of the canopy air space
      if(ixCasNrg/=integerMissing) scalarCanairNrgTrial = stateVec(ixCasNrg)

      ! extract canopy temperature
      if(ixVegNrg/=integerMissing) scalarCanopyNrgTrial = stateVec(ixVegNrg)

      ! extract intercepted water
      if(ixVegHyd/=integerMissing)then
        select case( ixStateType_subset(ixVegHyd) )
          case(iname_liqCanopy); scalarCanopyLiqTrial = stateVec(ixVegHyd)
          case(iname_watCanopy); scalarCanopyWatTrial = stateVec(ixVegHyd)
          case default; err=20; message=trim(message)//'case not found: expect iname_liqCanopy or iname_watCanopy'; return
        end select
      endif

    endif  ! not computing the vegetation flux

    ! *** extract state variables from the snow+soil sub-domain


    ! overwrite with the energy values from the state vector
    if(nSnowSoilNrg>0)then
      do concurrent (iLayer=1:nLayers,ixSnowSoilNrg(iLayer)/=integerMissing)   ! (loop through non-missing energy state variables in the snow+soil domain)
        mLayerNrgTrial(iLayer) = stateVec( ixSnowSoilNrg(iLayer) )
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    endif

    ! overwrite with the hydrology values from the state vector
    if(nSnowSoilHyd>0)then
      do concurrent (iLayer=1:nLayers,ixSnowSoilHyd(iLayer)/=integerMissing)   ! (loop through non-missing hydrology state variables in the snow+soil domain)
        select case( ixHydType(iLayer) )
          case(iname_watLayer); mLayerVolFracWatTrial(iLayer)          = stateVec( ixSnowSoilHyd(iLayer) ) ! total water state variable for snow+soil layers
          case(iname_liqLayer); mLayerVolFracLiqTrial(iLayer)          = stateVec( ixSnowSoilHyd(iLayer) ) ! liquid water state variable for snow+soil layers
          case(iname_matLayer); mLayerMatricHeadTrial(iLayer-nSnow)    = stateVec( ixSnowSoilHyd(iLayer) ) ! total water matric potential variable for soil layers
          case(iname_lmpLayer); mLayerMatricHeadLiqTrial(iLayer-nSnow) = stateVec( ixSnowSoilHyd(iLayer) ) ! liquid matric potential state variable for soil layers
          case default ! do nothing
        end select
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    endif

    ! extract temperature of the canopy air space
    if(ixAqWat/=integerMissing) scalarAquiferStorageTrial = stateVec(ixAqWat)

  end associate
end subroutine varExtract_host

subroutine varExtract(&
  nGRU, &
                       ! input
                       stateVec,                                  & ! intent(in):    model state vector (mixed units)
                       indx_data,                                 & ! intent(in):    indices defining model states and layers
                       ! output: variables for the vegetation canopy
                       scalarCanairNrgTrial,                      & ! intent(inout):   trial value of canopy air energy, temperature (K) or enthalpy (J m-3)
                       scalarCanopyNrgTrial,                      & ! intent(inout):   trial value of canopy energy, temperature (K) or enthalpy (J m-3)
                       scalarCanopyWatTrial,                      & ! intent(inout):   trial value of canopy total water (kg m-2)
                       scalarCanopyLiqTrial,                      & ! intent(inout):   trial value of canopy liquid water (kg m-2)
                       ! output: variables for the snow-soil domain
                       mLayerNrgTrial,                            & ! intent(inout):   trial vector of layer energy, temperature (K) or enthalpy (J m-3)
                       mLayerVolFracWatTrial,                     & ! intent(inout):   trial vector of volumetric total water content (-)
                       mLayerVolFracLiqTrial,                     & ! intent(inout):   trial vector of volumetric liquid water content (-)
                       mLayerMatricHeadTrial,                     & ! intent(inout):   trial vector of total water matric potential (m)
                       mLayerMatricHeadLiqTrial,                  & ! intent(inout):   trial vector of liquid water matric potential (m)
                       ! output: variables for the aquifer
                       scalarAquiferStorageTrial,                 & ! intent(out):   trial value of storage of water in the aquifer (m)
                       ! output: error control
                       err,message)                                 ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  use device_data_types
  implicit none
  integer(i4b),intent(in) :: nGRU
  ! input
  real(rkind),intent(in),device             :: stateVec(:,:)                     ! model state vector (mixed units)
  type(indx_data_device),intent(in)       :: indx_data                       ! indices defining model states and layers
  ! output: variables for the vegetation canopy
  real(rkind),intent(inout),device          :: scalarCanairNrgTrial(:)            ! trial value of canopy air energy, temperature (K) or enthalpy (J m-3)
  real(rkind),intent(inout),device          :: scalarCanopyNrgTrial(:)            ! trial value of canopy energy, temperature (K) or enthalpy (J m-3)
  real(rkind),intent(inout),device          :: scalarCanopyWatTrial(:)            ! trial value of canopy total water (kg m-2)
  real(rkind),intent(inout),device          :: scalarCanopyLiqTrial(:)            ! trial value of canopy liquid water (kg m-2)
  ! output: variables for the snow-soil domain
  real(rkind),intent(inout),device          :: mLayerNrgTrial(:,:)               ! trial vector of layer energy, temperature (K) or enthalpy (J m-3)
  real(rkind),intent(inout),device          :: mLayerVolFracWatTrial(:,:)        ! trial vector of volumetric total water content (-)
  real(rkind),intent(inout),device          :: mLayerVolFracLiqTrial(:,:)        ! trial vector of volumetric liquid water content (-)
  real(rkind),intent(inout),device          :: mLayerMatricHeadTrial(:,:)        ! trial vector of total water matric potential (m)
  real(rkind),intent(inout),device          :: mLayerMatricHeadLiqTrial(:,:)     ! trial vector of liquid water matric potential (m)
  ! output: variables for the aquifer
  real(rkind),intent(inout),device          :: scalarAquiferStorageTrial(:)       ! trial value of storage of water in the aquifer (m)
  ! output: error control
  integer(i4b),intent(out)           :: err                             ! error code
  character(*),intent(out)           :: message                         ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  integer(i4b)                       :: iLayer                          ! index of layer within the snow+soil domain
          type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  associate(&
    ! number of model layers, and layer type
    nSnow                   => indx_data%nSnow                 ,& ! intent(in):  [i4b]    total number of snow layers
    nSoil                   => indx_data%nSoil                 ,& ! intent(in):  [i4b]    total number of soil layers
    nLayers                 => indx_data%nLayers_d               ,& ! intent(in):  [i4b]    total number of snow and soil layers
    ! indices defining model states and layers
    ixCasNrg                => indx_data%ixCasNrg              ,& ! intent(in):  [i4b]    index of canopy air space energy state variable
    ixVegNrg                => indx_data%ixVegNrg              ,& ! intent(in):  [i4b]    index of canopy energy state variable
    ixVegHyd                => indx_data%ixVegHyd              ,& ! intent(in):  [i4b]    index of canopy hydrology state variable (mass)
    ixAqWat                 => indx_data%ixAqWat               ,& ! intent(in):  [i4b]    index of the squifer storage state variable
    ixSnowSoilNrg           => indx_data%ixSnowSoilNrg         ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for energy states in the snow+soil subdomain
    ixSnowSoilHyd           => indx_data%ixSnowSoilHyd         ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for hydrology states in the snow+soil subdomain
    ! nSnowSoilNrg            => indx_data%nSnowSoilNrg          ,& ! intent(in):  [i4b]    number of energy state variables in the snow+soil domain
    ! nSnowSoilHyd            => indx_data%nSnowSoilHyd          ,& ! intent(in):  [i4b]    number of hydrology variables in the snow+soil domain
    ! indices defining type of model state variables
    ixStateType_subset      => indx_data%ixStateType    ,& ! intent(in):  [i4b(:)] [state subset] type of desired model state variables
    ixHydType               => indx_data%ixHydType              & ! intent(in):  [i4b(:)] index of the type of hydrology states in snow+soil domain
    )! association with variables in the data structures

    ! --------------------------------------------------------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message='varExtract/'

    call varExtract_kernel<<<blocks,threads>>>(nGRU,&
stateVec,&
  ixCasNrg, ixVegNrg, ixVegHyd, ixAqWat,&
  nLayers, nSnow, &
  ixSnowSoilHyd,ixSnowSoilNrg, &
  ixStateType_subset,ixHydType,&
  scalarCanairNrgTrial, scalarCanopyNrgTrial,&
  scalarCanopyLiqTrial,scalarCanopyWatTrial,&
  mLayerNrgTrial, mLayerVolFracWatTrial, mLayerVolFracLiqTrial,mLayerMatricHeadTrial, mLayerMatricHeadLiqTrial, &
  scalarAquiferStorageTrial)

  end associate
end subroutine varExtract

attributes(global) subroutine varExtract_kernel(nGRU,&
stateVec,&
  ixCasNrg, ixVegNrg, ixVegHyd, ixAqWat,&
  nLayers, nSnow, &
  ixSnowSoilHyd,ixSnowSoilNrg, &
  ixStateType_subset,ixHydType,&
  scalarCanairNrgTrial, scalarCanopyNrgTrial,&
  scalarCanopyLiqTrial,scalarCanopyWatTrial,&
  mLayerNrgTrial, mLayerVolFracWatTrial, mLayerVolFracLiqTrial,mLayerMatricHeadTrial, mLayerMatricHeadLiqTrial, &
  scalarAquiferStorageTrial)
  implicit none
  integer(i4b),intent(in),value :: nGRU
  real(rkind),intent(in) :: stateVec(:,:)
  integer(i4b),intent(in) :: ixCasNrg(:), ixVegNrg(:), ixVegHyd(:), ixAqWat(:)
  integer(i4b),intent(in) :: nLayers(:), nSnow(:)
  integer(i4b),intent(in) :: ixSnowSoilHyd(:,:), ixSnowSoilNrg(:,:)
  integer(i4b),intent(in) :: ixStateType_subset(:,:), ixHydType(:,:)
  real(rkind),intent(inout) :: scalarCanairNrgTrial(:), scalarCanopyNrgTrial(:)
  real(rkind),intent(inout) :: scalarCanopyLiqTrial(:), scalarCanopyWatTrial(:)
  real(rkind),intent(inout) :: mLayerNrgTrial(:,:), mLayerVolFracWatTrial(:,:),mLayerVolFracLiqTrial(:,:), mLayerMatricHeadTrial(:,:), mLayerMatricHeadLiqTrial(:,:)
  real(rkind),intent(inout) :: scalarAquiferStorageTrial(:)


      integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

  call varExtract_device(stateVec(:,iGRU),&
  ixCasNrg(iGRU), ixVegNrg(iGRU), ixVegHyd(iGRU), ixAqWat(iGRU),&
  nLayers(iGRU), nSnow(iGRU), &
  ixSnowSoilHyd(:,iGRU),ixSnowSoilNrg(:,iGRU), &
  ixStateType_subset(:,iGRU),ixHydType(:,iGRU),&
  scalarCanairNrgTrial(iGRU), scalarCanopyNrgTrial(iGRU),&
  scalarCanopyLiqTrial(iGRU),scalarCanopyWatTrial(iGRU),&
  mLayerNrgTrial(:,iGRU), mLayerVolFracWatTrial(:,iGRU), mLayerVolFracLiqTrial(:,iGRU),mLayerMatricHeadTrial(:,iGRU), mLayerMatricHeadLiqTrial(:,iGRU), &
  scalarAquiferStorageTrial(iGRU))
end subroutine

attributes(device) subroutine varExtract_device(stateVec,&
  ixCasNrg, ixVegNrg, ixVegHyd, ixAqWat,&
  nLayers, nSnow, &
  ixSnowSoilHyd,ixSnowSoilNrg, &
  ixStateType_subset,ixHydType,&
  scalarCanairNrgTrial, scalarCanopyNrgTrial,&
  scalarCanopyLiqTrial,scalarCanopyWatTrial,&
  mLayerNrgTrial, mLayerVolFracWatTrial, mLayerVolFracLiqTrial,mLayerMatricHeadTrial, mLayerMatricHeadLiqTrial, &
  scalarAquiferStorageTrial)
  implicit none
  real(rkind),intent(in) :: stateVec(:)
  integer(i4b),intent(in) :: ixCasNrg, ixVegNrg, ixVegHyd, ixAqWat
  integer(i4b),intent(in) :: nLayers, nSnow
  integer(i4b),intent(in) :: ixSnowSoilHyd(:), ixSnowSoilNrg(:)
  integer(i4b),intent(in) :: ixStateType_subset(:), ixHydType(:)
  real(rkind),intent(inout) :: scalarCanairNrgTrial, scalarCanopyNrgTrial
  real(rkind),intent(inout) :: scalarCanopyLiqTrial, scalarCanopyWatTrial
  real(rkind),intent(inout) :: mLayerNrgTrial(:), mLayerVolFracWatTrial(:),mLayerVolFracLiqTrial(:), mLayerMatricHeadTrial(:), mLayerMatricHeadLiqTrial(:)
  real(rkind),intent(inout) :: scalarAquiferStorageTrial

  integer(i4b) :: iLayer
      ! *** extract state variables for the vegetation canopy

    ! check if computing the vegetation flux
    if(ixCasNrg/=integerMissing .or. ixVegNrg/=integerMissing .or. ixVegHyd/=integerMissing)then

      ! extract temperature of the canopy air space
      if(ixCasNrg/=integerMissing) scalarCanairNrgTrial = stateVec(ixCasNrg)

      ! extract canopy temperature
      if(ixVegNrg/=integerMissing) scalarCanopyNrgTrial = stateVec(ixVegNrg)

      ! extract intercepted water
      if(ixVegHyd/=integerMissing)then
        select case( ixStateType_subset(ixVegHyd) )
          case(iname_liqCanopy); scalarCanopyLiqTrial = stateVec(ixVegHyd)
          case(iname_watCanopy); scalarCanopyWatTrial = stateVec(ixVegHyd)
          ! case default; err=20; message=trim(message)//'case not found: expect iname_liqCanopy or iname_watCanopy'; return
        end select
      endif

    endif  ! not computing the vegetation flux

    ! *** extract state variables from the snow+soil sub-domain


    ! overwrite with the energy values from the state vector
    if(nLayers>0)then
      do iLayer=1,nLayers
        if (ixSnowSoilNrg(iLayer)/=integerMissing) then   ! (loop through non-missing energy state variables in the snow+soil domain)
        mLayerNrgTrial(iLayer) = stateVec( ixSnowSoilNrg(iLayer) )
        end if
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    endif

    ! overwrite with the hydrology values from the state vector
    if(nLayers>0)then
      do iLayer=1,nLayers
        if (ixSnowSoilHyd(iLayer)/=integerMissing) then   ! (loop through non-missing hydrology state variables in the snow+soil domain)
        select case( ixHydType(iLayer) )
          case(iname_watLayer); mLayerVolFracWatTrial(iLayer)          = stateVec( ixSnowSoilHyd(iLayer) ) ! total water state variable for snow+soil layers
          case(iname_liqLayer); mLayerVolFracLiqTrial(iLayer)          = stateVec( ixSnowSoilHyd(iLayer) ) ! liquid water state variable for snow+soil layers
          case(iname_matLayer); mLayerMatricHeadTrial(iLayer-nSnow)    = stateVec( ixSnowSoilHyd(iLayer) ) ! total water matric potential variable for soil layers
          case(iname_lmpLayer); mLayerMatricHeadLiqTrial(iLayer-nSnow) = stateVec( ixSnowSoilHyd(iLayer) ) ! liquid matric potential state variable for soil layers
          case default ! do nothing
        end select
      end if
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    endif

    ! extract temperature of the canopy air space
    if(ixAqWat/=integerMissing) scalarAquiferStorageTrial = stateVec(ixAqWat)
end subroutine

end module getVectorz_module
