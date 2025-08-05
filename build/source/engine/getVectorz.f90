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

! access the global print flag
USE globalData,only:globalPrintFlag

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
public::varExtract
public::checkFeas

! common variables
real(rkind),parameter :: valueMissing=-9999._rkind ! missing value

contains

! **********************************************************************************************************
! public subroutine popStateVec: populate model state vectors
! **********************************************************************************************************
subroutine popStateVec(&
                        ! input: data structures
                        nState,                  & ! intent(in):  number of desired state variables
                        nGRU, &
                        enthalpyStateVec,        & ! intent(in):  flag if enthalpy is state variable          
                        prog_data,               & ! intent(in):  model prognostic variables for a local HRU
                        diag_data,               & ! intent(in):  model diagnostic variables for a local HRU
                        indx_data,               & ! intent(in):  indices defining model states and layers
                        ! output
                        stateVec,                & ! intent(out): model state vector
                        err,message)               ! intent(out): error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: data structures
                        use device_data_types
  integer(i4b),intent(in)         :: nState                 ! number of desired state variables
  integer(i4b) :: nGRU
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
  ! state subsets
  integer(i4b)                    :: iState                 ! index of state within the snow+soil domain
  integer(i4b)                    :: iLayer                 ! index of layer within the snow+soil domain
  integer(i4b)                    :: ixStateSubset          ! index within the state subset
  ! logical(lgt),dimension(nState)  :: stateFlag              ! flag to denote that the state is populated
  integer(i4b) :: iGRU
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
    nSnowSoilNrg         => indx_data%nLayers_d       ,& ! intent(in) : [i4b]    number of energy state variables in the snow+soil domain
    nSnowSoilHyd         => indx_data%nLayers_d       ,& ! intent(in) : [i4b]    number of hydrology state variables in the snow+soil domain
    ! type of model state variabless
    ixStateType_subset   => indx_data%ixStateType_subset     ,& ! intent(in) : [i4b(:)] [state subset] type of desired model state variables
    ixHydType            => indx_data%ixHydType              ,& ! intent(in) : [i4b(:)] index of the type of hydrology states in snow+soil domain
    ! number of layers
    nSnow                => indx_data%nSnow               ,& ! intent(in) : [i4b]    number of snow layers
    nSoil                => indx_data%nSoil               ,& ! intent(in) : [i4b]    number of soil layers
    nLayers              => indx_data%nLayers_d              & ! intent(in) : [i4b]    total number of layers
    )  ! end association with variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='popStateVec/'

    ! -----
    ! * initialize state vectors...
    ! -----------------------------

    stateVec = 0._rkind
    ! build the state vector for the temperature of the canopy air space
    ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer canopy
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
      ! print*, iGRU, 'CasNrg', ixCasNrg(iGRU)
      if(enthalpyStateVec)then
        stateVec( ixCasNrg(iGRU),iGRU ) = scalarCanairEnthalpy(iGRU) ! transfer canopy air enthalpy to the state vector
      else
        stateVec( ixCasNrg(iGRU),iGRU ) = scalarCanairTemp(iGRU)     ! transfer canopy air temperature to the state vector
      endif
      ! print*, iGRU, 'VegNrg', ixVegNrg(iGRU)

    ! build the state vector for the temperature of the vegetation canopy
    ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer canopy
      if(enthalpyStateVec)then
        stateVec( ixVegNrg(iGRU),iGRU ) = scalarCanopyEnthalpy(iGRU) ! transfer vegetation enthalpy to the state vector
      else
        stateVec( ixVegNrg(iGRU),iGRU )  = scalarCanopyTemp(iGRU)    ! transfer vegetation temperature to the state vector
      endif

    ! build the state vector for the water in the vegetation canopy
    ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer canopy
      select case(ixStateType_subset( ixVegHyd(iGRU),iGRU ))
        case(iname_watCanopy); stateVec( ixVegHyd(iGRU),iGRU ) = scalarCanopyWat(iGRU) ! transfer total canopy water to the state vector
        case(iname_liqCanopy); stateVec( ixVegHyd(iGRU),iGRU ) = scalarCanopyLiq(iGRU) ! transfer liquid canopy water to the state vector
      end select
            ! print*, iGRU, 'VegHyd', ixVegHyd(iGRU)

          ! build the state vector for the aquifer storage
    ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer aquifer
      stateVec( ixAqWat(iGRU),iGRU )  = scalarAquiferStorage(iGRU)    ! transfer aquifer storage to the state vector

    end do

    ! build the energy state vector for the snow and soil domain
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        do iLayer=1,nLayers(iGRU)
          if (ixSnowSoilNrg(iLayer,iGRU)/=integerMissing) then
            ixStateSubset            = ixSnowSoilNrg(iLayer,iGRU) ! index within the state vector
            if(enthalpyStateVec)then
              stateVec(ixStateSubset,iGRU) = mLayerEnthalpy(iLayer,iGRU) ! transfer enthalpy from a layer to the state vector
            else
              stateVec(ixStateSubset,iGRU) = mLayerTemp(iLayer,iGRU)     ! transfer temperature from a layer to the state vector
            endif
          end if
        end do
      end do  ! looping through non-missing energy state variables in the snow+soil domain

    ! build the hydrology state vector for the snow+soil domains
    ! NOTE: ixVolFracWat  and ixVolFracLiq can also include states in the soil domain, hence enable primary variable switching
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        do iLayer=1,nLayers(iGRU)
          if (ixSnowSoilHyd(iLayer,iGRU)/=integerMissing) then
            ixStateSubset            = ixSnowSoilHyd(iLayer,iGRU) ! index within the state vector
            select case( ixHydType(iLayer,iGRU) )
              case(iname_watLayer); stateVec(ixStateSubset,iGRU) = mLayerVolFracWat(iLayer,iGRU)           ! total water state variable for snow+soil layers
              case(iname_liqLayer); stateVec(ixStateSubset,iGRU) = mLayerVolFracLiq(iLayer,iGRU)           ! liquid water state variable for snow+soil layers
              case(iname_matLayer); stateVec(ixStateSubset,iGRU) = mLayerMatricHead(iLayer-nSnow(iGRU),iGRU)     ! total water matric potential variable for soil layers
              case(iname_lmpLayer); stateVec(ixStateSubset,iGRU) = mLayerMatricHeadLiq(iLayer-nSnow(iGRU),iGRU)  ! liquid matric potential state variable for soil layers
          
            end select
          end if
        end do
      end do  ! looping through non-missing energy state variables in the snow+soil domain

    ! build the state vector for the aquifer storage
    ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer aquifer
    ! do concurrent (iState=1:size(ixAqWat),ixAqWat(iState)/=integerMissing)
    !   stateVec( ixAqWat(iState),iGRU )  = scalarAquiferStorage    ! transfer aquifer storage to the state vector
    ! end do

    ! ! check that we populated all state variables
    ! if(count(stateFlag)/=nState)then
    !   print*, 'stateFlag = ', stateFlag
    !   message=trim(message)//'some state variables unpopulated'
    !   err=20; return
    ! endif

    !       !$cuf kernel do(1) <<<*,*>>>
    ! do iGRU=1,nGRU
    !   do iLayer=1,size(stateVec,1)
    !     print*, iGRU, iLayer, stateVec(iLayer,iGRU)
    !   end do
    ! end do

  end associate fixedLength      ! end association to variables in the data structure where vector length does not change
end subroutine popStateVec


! **********************************************************************************************************
! public subroutine getScaling: get scale factors
! **********************************************************************************************************
subroutine getScaling(&
                      ! input: data structures
                      diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                      indx_data,               & ! intent(in):    indices defining model states and layers
                      nGRU, &
                      ! output
                      fScale,                  & ! intent(out):   characteristic scale of the function evaluations (mixed units)
                      xScale,                  & ! intent(out):   variable scaling vector (mixed units)
                      sMul,                    & ! intent(out):   multiplier for state vector (used in the residual calculations)
                      dMat,                    & ! intent(out):   diagonal of the Jacobian matrix (excludes fluxes)
                      err,message)               ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  USE nr_utility_module,only:arth                   ! get a sequence of numbers arth(start, incr, count)
  USE f2008funcs_module,only:findIndex              ! finds the index of the first value within a vector
  use device_data_types
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! input: data structures
  type(diag_data_device),intent(in)    :: diag_data              ! diagnostic variables for a local HRU
  type(indx_data_device),intent(in)    :: indx_data              ! indices defining model states and layers
  integer(i4b),intent(in) :: nGRU
  ! output: state vectors
  real(rkind),intent(out),device         :: fScale(:,:)              ! characteristic scale of the function evaluations (mixed units)
  real(rkind),intent(out),device         :: xScale(:,:)              ! variable scaling vector (mixed units)
  real(qp),intent(out),device            :: sMul(:,:)    ! NOTE: qp  ! multiplier for state vector (used in the residual calculations)
  real(rkind),intent(out),device         :: dMat(:,:)                ! diagonal of the Jacobian matrix (excludes fluxes)
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
  integer(i4b)                    :: iLayer,iGRU                 ! index of layer within the snow+soil domain
  integer(i4b)                    :: ixStateSubset          ! index within the state subset
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  fixedLength: associate(&
    ! model diagnostic variables
    canopyDepth         => diag_data%scalarCanopyDepth      ,& ! intent(in):  [dp]     canopy depth (m)
    volHeatCapVeg       => diag_data%scalarBulkVolHeatCapVeg,& ! intent(in) : [dp]     bulk volumetric heat capacity of vegetation (J m-3 K-1)
    mLayerVolHeatCap    => diag_data%mLayerVolHtCapBulk_m     ,& ! intent(in) : [dp(:)]  bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
    ! indices defining specific model states
    ixCasNrg            => indx_data%ixCasNrg                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy air space energy state variable
    ixVegNrg            => indx_data%ixVegNrg                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy energy state variable
    ixVegHyd            => indx_data%ixVegHyd                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy hydrology state variable (mass)
    ! vector of energy and hydrology indices for the snow and soil domains
    ixSnowSoilNrg_m       => indx_data%ixSnowSoilNrg            ,& ! intent(in) : [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    ixSnowSoilHyd_m       => indx_data%ixSnowSoilHyd            ,& ! intent(in) : [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    nSnowSoilNrg        => indx_data%nLayers_d         ,& ! intent(in) : [i4b]    number of energy state variables in the snow+soil domain
    nSnowSoilHyd        => indx_data%nLayers_d         ,& ! intent(in) : [i4b]    number of hydrology state variables in the snow+soil domain
    ! type of model state variabless
    ixStateType_subset_m  => indx_data%ixStateType_subset       ,& ! intent(in) : [i4b(:)] [state subset] type of desired model state variables
    ! number of layers
    nSnow               => indx_data%nSnow                 ,& ! intent(in) : [i4b]    number of snow layers
    nSoil               => indx_data%nSoil                 ,& ! intent(in) : [i4b]    number of soil layers
    nLayers             => indx_data%nLayers_d                & ! intent(in) : [i4b]    total number of layers
    )  ! end association with variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='getScaling/'

    ! -----
    ! * define scaling vectors...
    ! ---------------------------

    !$cuf kernel do(2) <<<*,*>>>
    do iGRU=1,nGRU
    do iLayer=1,size(fScale)
      ! define the function and variable scaling factors for energy
      if (ixStateType_subset_m(iLayer,iGRU)==iname_nrgCanair .or. ixStateType_subset_m(iLayer,iGRU)==iname_nrgCanopy .or. ixStateType_subset_m(iLayer,iGRU)==iname_nrgLayer) then
        fScale(iLayer,iGRU) = 1._rkind / fScaleNrg  ! 1/(J m-3)
        xScale(iLayer,iGRU) = 1._rkind  ! K
      endif
      ! define the function and variable scaling factors for water on the vegetation canopy
      if(ixStateType_subset_m(iLayer,iGRU)==iname_watCanopy .or. ixStateType_subset_m(iLayer,iGRU)==iname_liqCanopy) then
        fScale(iLayer,iGRU) = 1._rkind / (fScaleLiq*canopyDepth(iGRU)*iden_water)  ! 1/(kg m-2)
        xScale(iLayer,iGRU) = 1._rkind  ! (kg m-2)
      endif
          ! define the function and variable scaling factors for water in the snow+soil domain
      if (ixStateType_subset_m(iLayer,iGRU)==iname_watLayer .or. ixStateType_subset_m(iLayer,iGRU)==iname_liqLayer) then
        fScale(iLayer,iGRU) = 1._rkind / fScaleLiq  ! (-)
        xScale(iLayer,iGRU) = 1._rkind  ! (-)
      end if
      ! define the function and variable scaling factors for water in the snow+soil domain
      if (ixStateType_subset_m(iLayer,iGRU)==iname_matLayer .or. ixStateType_subset_m(iLayer,iGRU)==iname_lmpLayer) then
        fScale(iLayer,iGRU) = 1._rkind / fScaleLiq  ! (-)
        xScale(iLayer,iGRU) = 1._rkind  ! (m)
      end if

      ! define the function and variable scaling factors for water storage in the aquifer
      if(ixStateType_subset_m(iLayer,iGRU)==iname_watAquifer) then
        fScale(iLayer,iGRU) = 1._rkind
        xScale(iLayer,iGRU) = 1._rkind
      endif
      ! define the multiplier for the state vector for residual calculations (vegetation canopy)
      ! NOTE: Use the "where" statement to generalize to multiple canopy layers (currently one canopy layer)
      if (ixStateType_subset_m(iLayer,iGRU)==iname_nrgCanair) sMul(iLayer,iGRU) = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
      if(ixStateType_subset_m(iLayer,iGRU)==iname_nrgCanopy) sMul(iLayer,iGRU) = volHeatCapVeg(iGRU)     ! volumetric heat capacity of the vegetation (J m-3 K-1)
      if(ixStateType_subset_m(iLayer,iGRU)==iname_watCanopy) sMul(iLayer,iGRU) = 1._rkind          ! nothing else on the left hand side
      if(ixStateType_subset_m(iLayer,iGRU)==iname_liqCanopy) sMul(iLayer,iGRU) = 1._rkind          ! nothing else on the left hand side
      ! compute terms in the Jacobian for vegetation (excluding fluxes)
      ! NOTE: This is computed outside the iteration loop because it does not depend on state variables
      ! NOTE: Energy for vegetation is computed *within* the iteration loop as it includes phase change
      ! NOTE: Use the "where" statement to generalize to multiple canopy layers (currently one canopy layer)
      if (ixStateType_subset_m(iLayer,iGRU)==iname_nrgCanair) dMat(iLayer,iGRU) = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
      if (ixStateType_subset_m(iLayer,iGRU)==iname_nrgCanopy) dMat(iLayer,iGRU) = realMissing       ! populated within the iteration loop
      if (ixStateType_subset_m(iLayer,iGRU)==iname_watCanopy) dMat(iLayer,iGRU) = 1._rkind          ! nothing else on the left hand side
      if (ixStateType_subset_m(iLayer,iGRU)==iname_liqCanopy) dMat(iLayer,iGRU) = 1._rkind          ! nothing else on the left hand side
      ! define the scaling factor and diagonal elements for the aquifer
      if (ixStateType_subset_m(iLayer,iGRU)==iname_watAquifer) then
        sMul(iLayer,iGRU) = 1._rkind
        dMat(iLayer,iGRU) = 1._rkind
      endif
  
    end do
  end do

    ! where(ixStateType_subset==iname_watCanopy .or. ixStateType_subset==iname_liqCanopy)
    !   fScale = 1._rkind / (fScaleLiq*canopyDepth*iden_water)  ! 1/(kg m-2)
    !   xScale = 1._rkind  ! (kg m-2)
    ! endwhere



    ! -----
    ! * define components of derivative matrices that are constant over a time step (substep)...
    ! ------------------------------------------------------------------------------------------

    ! define the multiplier for the state vector for residual calculations (vegetation canopy)
    ! ! NOTE: Use the "where" statement to generalize to multiple canopy layers (currently one canopy layer)
    ! where(ixStateType_subset==iname_nrgCanair) sMul = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
    ! where(ixStateType_subset==iname_nrgCanopy) sMul = volHeatCapVeg     ! volumetric heat capacity of the vegetation (J m-3 K-1)
    ! where(ixStateType_subset==iname_watCanopy) sMul = 1._rkind          ! nothing else on the left hand side
    ! where(ixStateType_subset==iname_liqCanopy) sMul = 1._rkind          ! nothing else on the left hand side

    ! compute terms in the Jacobian for vegetation (excluding fluxes)
    ! NOTE: This is computed outside the iteration loop because it does not depend on state variables
    ! NOTE: Energy for vegetation is computed *within* the iteration loop as it includes phase change
    ! NOTE: Use the "where" statement to generalize to multiple canopy layers (currently one canopy layer)
    ! where(ixStateType_subset==iname_nrgCanair) dMat = Cp_air*iden_air   ! volumetric heat capacity of air (J m-3 K-1)
    ! where(ixStateType_subset==iname_nrgCanopy) dMat = realMissing       ! populated within the iteration loop
    ! where(ixStateType_subset==iname_watCanopy) dMat = 1._rkind          ! nothing else on the left hand side
    ! where(ixStateType_subset==iname_liqCanopy) dMat = 1._rkind          ! nothing else on the left hand side

    ! define the energy multiplier and diagonal elements for the state vector for residual calculations (snow-soil domain)
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      do iLayer=1,nLayers(iGRU)
        if (ixSnowSoilNrg_m(iLayer,iGRU)/=integerMissing) then   ! (loop through non-missing energy state variables in the snow+soil domain)
        ixStateSubset        = ixSnowSoilNrg_m(iLayer,iGRU)      ! index within the state vector
        sMul(ixStateSubset,iGRU)  = mLayerVolHeatCap(iLayer,iGRU)   ! transfer volumetric heat capacity to the state multiplier
        dMat(ixStateSubset,iGRU)  = realMissing                ! diagonal element populated within the iteration loop
        endif
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    end do

    ! define the hydrology multiplier and diagonal elements for the state vector for residual calculations (snow-soil domain)
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      do iLayer=1,nLayers(iGRU)
        if (ixSnowSoilHyd_m(iLayer,iGRU)/=integerMissing) then   ! (loop through non-missing energy state variables in the snow+soil domain)
        ixStateSubset        = ixSnowSoilHyd_m(iLayer,iGRU)      ! index within the state vector
        sMul(ixStateSubset,iGRU)  = 1._rkind                   ! state multiplier = 1 (nothing else on the left-hand-side)
        dMat(ixStateSubset,iGRU)  = 1._rkind                   ! diagonal element = 1 (nothing else on the left-hand-side)
        endif
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    end do


  end associate fixedLength      ! end association to variables in the data structure where vector length does not change
end subroutine getScaling


! **********************************************************************************************************
! public subroutine checkFeas: check feasibility of the state vector
! **********************************************************************************************************
subroutine checkFeas(&
                      ! input
                      stateVec,                                  & ! intent(in):    model state vector (mixed units)
                      nGRU, &
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
                      use device_data_types
  implicit none
  ! input
  real(rkind),intent(in),device          :: stateVec(:,:)               ! model state vector (mixed units)
  integer(i4b) :: nGRU
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
  integer(i4b)                    :: iLayer,iGRU                    ! index of layer within the snow+soil domain
  real(rkind)                     :: xMin,xMax                 ! minimum and maximum values for water content
  real(rkind),parameter           :: canopyTempMax=500._rkind  ! expected maximum value for the canopy temperature (K)

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  associate(&
    ! soil parameters
    theta_sat               => mpar_data%theta_sat       ,&  ! intent(in): [dp(:)]  soil porosity (-)
    theta_res               => mpar_data%theta_res       ,&  ! intent(in): [dp(:)]  residual volumetric water content (-)
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
    ixSnowOnlyNrg_m           => indx_data%ixSnowOnlyNrg   ,&  ! intent(in): [i4b(:)] indices for energy states in the snow subdomain
    ixSnowSoilHyd_m           => indx_data%ixSnowSoilHyd   ,&  ! intent(in): [i4b(:)] indices for hydrology states in the snow+soil subdomain
    ixStateType_m             => indx_data%ixStateType     ,&  ! intent(in): [i4b(:)] indices defining the type of the state (iname_nrgLayer...)
    ixHydCanopy_m             => indx_data%ixHydCanopy     ,&  ! intent(in): [i4b(:)] index of the hydrology states in the canopy domain
    ixHydType_m               => indx_data%ixHydType       ,&  ! intent(in): [i4b(:)] index of the type of hydrology states in snow+soil domain
    layerType_m               => indx_data%layerType        &  ! intent(in): [i4b(:)] layer type (iname_soil or iname_snow)
    )! association with variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message="checkFeas/"


    !  NOTE: we will not print infeasibilities since it does not indicate a failure, just a need to iterate until maxiter
    feasible=.true.
    ! check that the canopy air space temperature is reasonable
    !$cuf kernel do(1) <<<*,*>>> reduce(.and.:feasible)
    do iGRU=1,nGRU
    if(ixCasNrg(iGRU)/=integerMissing)then
      if(stateVec(ixCasNrg(iGRU),iGRU) > canopyTempMax .and. .not.enthalpyStateVec)then 
        print*, iGRU, 'ixCasNrg', stateVec(ixCasNrg(iGRU),iGRU)
        feasible=.false.
        ! message=trim(message)//'canopy air space temp high/'
        !write(*,'(a,1x,L1,1x,10(f20.10,1x))') 'feasible, max, stateVec( ixCasNrg )', feasible, canopyTempMax, stateVec(ixCasNrg)
      endif
    endif

    ! check that the canopy air space temperature is reasonable
    if(ixVegNrg(iGRU)/=integerMissing)then
      if(stateVec(ixVegNrg(iGRU),iGRU) > canopyTempMax .and. .not.enthalpyStateVec)then
        print*, iGRU, 'ixVegNrg', stateVec(ixVegNrg(iGRU),iGRU)
        feasible=.false.
        ! message=trim(message)//'canopy temp high/'
        !write(*,'(a,1x,L1,1x,10(f20.10,1x))') 'feasible, max, stateVec( ixVegNrg )', feasible, canopyTempMax, stateVec(ixVegNrg)
      endif
    endif

    ! check canopy liquid water is not negative
    if(ixVegHyd(iGRU)/=integerMissing)then
      if(stateVec(ixVegHyd(iGRU),iGRU) < -1e-10)then 
        print*, iGRU, 'ixVegHyd', stateVec(ixVegHyd(iGRU),iGRU)
        feasible=.false.
        ! print*, stateVec(ixVegHyd)
        ! message=trim(message)//'canopy liq water neg/'
        !write(*,'(a,1x,L1,1x,10(f20.10,1x))') 'feasible, min, stateVec( ixVegHyd )', feasible, 0._rkind, stateVec(ixVegHyd)
      endif
    endif
  end do

  err = cudaDeviceSynchronize()
  ! print*, 692, err
    ! check snow temperature is below freezing
  !$cuf kernel do(1) <<<*,*>>> reduce(.and.:feasible)
  do iGRU=1,nGRU
    do iLayer=1,nSnow(iGRU)
      if (ixSnowOnlyNrg_m(iLayer,iGRU)/=integerMissing .and. stateVec(ixSnowOnlyNrg_m(iLayer,iGRU),iGRU) > Tfreeze .and. .not.enthalpyStateVec) then
        print*, iLayer, iGRU, 'ixSnowNrg', stateVec(ixSnowOnlyNrg_m(iLayer,iGRU),iGRU)
        feasible=.false.
      end if
    end do
  end do
    ! if(count(ixSnowOnlyNrg/=integerMissing)>0)then
      ! if(any(stateVec( pack(ixSnowOnlyNrg,ixSnowOnlyNrg/=integerMissing) ) > Tfreeze) .and. .not.enthalpyStateVec)then
        ! feasible=.false.
        ! message=trim(message)//'snow temp high/'
        !do iLayer=1,nSnow
        !  if(stateVec(ixSnowOnlyNrg(iLayer)) > Tfreeze) write(*,'(a,1x,i4,1x,L1,1x,10(f20.10,1x))') 'iLayer, feasible, max, stateVec( ixSnowOnlyNrg(iLayer) )', iLayer, feasible, Tfreeze, stateVec( ixSnowOnlyNrg(iLayer) )
        !enddo
      ! endif
    ! endif
  err = cudaDeviceSynchronize()
  ! print*, 710, err

    ! loop through non-missing hydrology state variables in the snow+soil domain
  !$cuf kernel do(1) <<<*,*>>> reduce(.and.:feasible)
  do iGRU=1,nGRU
    do iLayer=1,nLayers(iGRU)
      if(ixSnowSoilHyd_m(iLayer,iGRU)/=integerMissing) then

      ! check the minimum and maximum water constraints
      if(ixHydType_m(iLayer,iGRU)==iname_watLayer .or. ixHydType_m(iLayer,iGRU)==iname_liqLayer)then

        ! --> minimum
        if (layerType_m(iLayer,iGRU) == iname_soil) then
          xMin = theta_res(iLayer-nSnow(iGRU))
        else
          xMin = 0._rkind
        endif

        ! --> maximum
        select case( layerType_m(iLayer,iGRU) )
          case(iname_snow); xMax = merge(1._rkind, 1._rkind - mLayerVolFracIce(iLayer,iGRU), ixHydType_m(iLayer,iGRU)==iname_watLayer)
          case(iname_soil); xMax = merge(theta_sat(iLayer-nSnow(iGRU)), theta_sat(iLayer-nSnow(iGRU)) - mLayerVolFracIce(iLayer,iGRU), ixHydType_m(iLayer,iGRU)==iname_watLayer)
        end select

        ! --> check
        if(stateVec( ixSnowSoilHyd_m(iLayer,iGRU) ,iGRU) < xMin .or. stateVec( ixSnowSoilHyd_m(iLayer,iGRU) ,iGRU) > xMax)then 
          print*, iGRU, iLayer, 'ssHyd', stateVec( ixSnowSoilHyd_m(iLayer,iGRU) ,iGRU), xMin, xMax
          feasible=.false.
          ! message=trim(message)//'layer water out of bounds/'
          !if(stateVec( ixSnowSoilHyd(iLayer) ) < xMin .or. stateVec( ixSnowSoilHyd(iLayer) ) > xMax) &
          !write(*,'(a,1x,i4,1x,L1,1x,10(f20.10,1x))') 'iLayer, feasible, stateVec( ixSnowSoilHyd(iLayer) ), xMin, xMax = ', iLayer, feasible, stateVec( ixSnowSoilHyd(iLayer) ), xMin, xMax
        endif
      endif  ! if water states
    end if

    end do  ! loop through non-missing hydrology state variables in the snow+soil domain
  end do

  end associate    ! end association to variables in the data structure
end subroutine checkFeas
  

! **********************************************************************************************************
! public subroutine varExtract: extract variables from the state vector and compute diagnostic variables
!  This routine does not initialize any of the variables
! **********************************************************************************************************
subroutine varExtract(&
                       ! input
                       stateVec,                                  & ! intent(in):    model state vector (mixed units)
                       indx_data,                                 & ! intent(in):    indices defining model states and layers
                       nGRU, &
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
  ! input
  real(rkind),intent(in),device             :: stateVec(:,:)                     ! model state vector (mixed units)
  type(indx_data_device),intent(in)       :: indx_data                       ! indices defining model states and layers
  integer(i4b) :: nGRU
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
  integer(i4b) :: iGRU

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
    ixSnowSoilNrg_m           => indx_data%ixSnowSoilNrg         ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for energy states in the snow+soil subdomain
    ixSnowSoilHyd_m           => indx_data%ixSnowSoilHyd         ,& ! intent(in):  [i4b(:)] indices IN THE STATE SUBSET for hydrology states in the snow+soil subdomain
    nSnowSoilNrg            => indx_data%nLayers_d         ,& ! intent(in):  [i4b]    number of energy state variables in the snow+soil domain
    nSnowSoilHyd            => indx_data%nLayers_d         ,& ! intent(in):  [i4b]    number of hydrology variables in the snow+soil domain
    ! indices defining type of model state variables
    ixStateType_subset_m      => indx_data%ixStateType_subset       ,& ! intent(in):  [i4b(:)] [state subset] type of desired model state variables
    ixHydType_m               => indx_data%ixHydType                 & ! intent(in):  [i4b(:)] index of the type of hydrology states in snow+soil domain
    )! association with variables in the data structures

    ! --------------------------------------------------------------------------------------------------------------------------------
    ! --------------------------------------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message='varExtract/'

    ! *** extract state variables for the vegetation canopy

    ! check if computing the vegetation flux
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
    if(ixCasNrg(iGRU)/=integerMissing .or. ixVegNrg(iGRU)/=integerMissing .or. ixVegHyd(iGRU)/=integerMissing)then

      ! extract temperature of the canopy air space
      if(ixCasNrg(iGRU)/=integerMissing) scalarCanairNrgTrial(iGRU) = stateVec(ixCasNrg(iGRU),iGRU)

      ! extract canopy temperature
      if(ixVegNrg(iGRU)/=integerMissing) scalarCanopyNrgTrial(iGRU) = stateVec(ixVegNrg(iGRU),iGRU)

      ! extract intercepted water
      if(ixVegHyd(iGRU)/=integerMissing)then
        select case( ixStateType_subset_m(ixVegHyd(iGRU),iGRU) )
          case(iname_liqCanopy); scalarCanopyLiqTrial(iGRU) = stateVec(ixVegHyd(iGRU),iGRU)
          case(iname_watCanopy); scalarCanopyWatTrial(iGRU) = stateVec(ixVegHyd(iGRU),iGRU)
          ! case default; err=20; message=trim(message)//'case not found: expect iname_liqCanopy or iname_watCanopy'; return
        end select
      endif

    endif  ! not computing the vegetation flux
        ! extract temperature of the canopy air space
    if(ixAqWat(iGRU)/=integerMissing) scalarAquiferStorageTrial(iGRU) = stateVec(ixAqWat(iGRU),iGRU)

    enddo

    ! *** extract state variables from the snow+soil sub-domain


    ! overwrite with the energy values from the state vector
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      do iLayer=1,nLayers(iGRU)
        if (ixSnowSoilNrg_m(iLayer,iGRU)/=integerMissing) then  ! (loop through non-missing energy state variables in the snow+soil domain)
        mLayerNrgTrial(iLayer,iGRU) = stateVec( ixSnowSoilNrg_m(iLayer,iGRU),iGRU )
        endif
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    enddo

    ! overwrite with the hydrology values from the state vector
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
      do iLayer=1,nLayers(iGRU)
        if (ixSnowSoilHyd_m(iLayer,iGRU)/=integerMissing) then   ! (loop through non-missing hydrology state variables in the snow+soil domain)
        select case( ixHydType_m(iLayer,iGRU) )
          case(iname_watLayer); mLayerVolFracWatTrial(iLayer,iGRU)          = stateVec( ixSnowSoilHyd_m(iLayer,iGRU),iGRU ) ! total water state variable for snow+soil layers
          case(iname_liqLayer); mLayerVolFracLiqTrial(iLayer,iGRU)          = stateVec( ixSnowSoilHyd_m(iLayer,iGRU),iGRU ) ! liquid water state variable for snow+soil layers
          case(iname_matLayer); mLayerMatricHeadTrial(iLayer-nSnow(iGRU),iGRU)    = stateVec( ixSnowSoilHyd_m(iLayer,iGRU),iGRU ) ! total water matric potential variable for soil layers
          case(iname_lmpLayer); mLayerMatricHeadLiqTrial(iLayer-nSnow(iGRU),iGRU) = stateVec( ixSnowSoilHyd_m(iLayer,iGRU),iGRU ) ! liquid matric potential state variable for soil layers
          case default ! do nothing
        end select
      endif
      end do  ! looping through non-missing energy state variables in the snow+soil domain
    enddo


  end associate
end subroutine varExtract

end module getVectorz_module
