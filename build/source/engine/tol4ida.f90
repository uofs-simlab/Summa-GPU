module tol4ida_module

!======= Inclusions ===========
use, intrinsic :: iso_c_binding
use nrtype
use type4ida

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

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


! privacy
implicit none
private
public::computWeight4ida
public::popTol4ida


contains

! **********************************************************************************************************
! public function computWeight4ida: compute w_i = 1 / ( rtol_i * y_i + atol_i )
! **********************************************************************************************************
! Return values:
!    0 = success,
!   -1 = non-recoverable error, NaN or negative values
! ----------------------------------------------------------------
integer(c_int) function computWeight4ida(sunvec_y, sunvec_ewt, user_data) &
  result(ierr) bind(C,name='computWeight4ida')

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod
  use fnvector_cuda_mod
  use nrtype
  use type4ida
  use cudafor

  !======= Declarations =========
  implicit none

  ! calling variables
  type(N_Vector)          :: sunvec_y  ! solution N_Vector    y
  type(N_Vector)          :: sunvec_ewt ! derivative N_Vector W
  type(c_ptr), value      :: user_data ! user-defined data

  ! pointers to data in SUNDIALS vectors
  type(data4ida), pointer :: tol_data ! equations data
  real(rkind), pointer,device    :: stateVec(:,:)
  real(rkind), pointer,device    :: weightVec(:,:)
  integer(c_int)          :: iState

  integer(i4b) :: iGRU
  real(rkind),device :: rTol_undef, aTol_undef
  real(rkind),device :: nState2
  !======= Internals ============

  ! get equations data from user-defined data
  call c_f_pointer(user_data, tol_data)
  rTol_undef = 0.01_rkind
  aTol_undef = 0.01_rkind

  ! get data arrays from SUNDIALS vectors
  stateVec(1:tol_data%nState,1:tol_data%nGRU)  => FN_VGetDeviceArrayPointer_Cuda(sunvec_y)
  weightVec(1:tol_data%nState,1:tol_data%nGRU)  => FN_VGetDeviceArrayPointer_Cuda(sunvec_ewt)
  nState2 = real(tol_data%nState)

  ! print*, shape(stateVec)

  associate(rtol => tol_data%rtol, atol => tol_data%atol,nState => tol_data%nState,nGRU => tol_data%nGRU, &
    nState_i => tol_data%indx_data%nState)
  !$cuf kernel do(2) <<<*,*>>>
  do iGRU=1,nGRU
    do iState = 1,nState
      ! print*, iState, rtol(iState), atol(iState)
      ! print*, iState,iGRU, abs(stateVec(iState,iGRU))
      if (iState .le. nState_i(iGRU)) then
        weightVec(iState,iGRU) = rtol(iState,iGRU) * abs( stateVec(iState,iGRU) ) + atol(iState,iGRU)
        weightVec(iState,iGRU) = 1._rkind / weightVec(iState,iGRU) * (nState2 / nState_i(iGRU)) ** 0.5_rkind
      else
        weightVec(iState,iGRU) = rTol_undef * abs(stateVec(iState,iGRU)) + aTol_undef
        weightVec(iState,iGRU) = 0._rkind!1._rkind / weightVec(iState,iGRU)
      end if
      ! print*, iState, iGRU, nState2/real(nState_i(iGRU)), weightVec(iState,iGRU), stateVec(iState,iGRU)
    end do
  end do
  end associate

  ! print*, weightVec

  ierr = cudaDeviceSynchronize()
  ! print*, 115, ierr
  ierr = 0
  return

end function computWeight4ida


! **********************************************************************************************************
! public subroutine popTol4ida: populate tolerances for state vectors
! **********************************************************************************************************
subroutine popTol4ida(&
  nGRU, &
                      ! input: data structures
                      nState,                  & ! intent(in):    number of desired state variables
                      prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                      diag_data,               & ! intent(in):    model diagnostic variables for a local HRU
                      indx_data,               & ! intent(in):    indices defining model states and layers
                      mpar_data,               & ! intent(in)
                      ! output
                      absTol,                  & ! intent(out):   model state vector
                      relTol,                  &
                      err,message)               ! intent(out):   error control
  ! --------------------------------------------------------------------------------------------------------------------------------
  use device_data_types
  ! input: data structures
  integer(i4b) :: nGRU
  integer(i4b),intent(in)         :: nState                 ! number of desired state variables
  type(prog_data_device),intent(in)    :: prog_data              ! prognostic variables for a local HRU
  type(diag_data_device),intent(in)    :: diag_data              ! diagnostic variables for a local HRU
  type(indx_data_device),intent(in)    :: indx_data              ! indices defining model states and layers
  type(mpar_data_device),intent(in)    :: mpar_data              ! model parameters
  ! output
  real(rkind),intent(out),device         :: absTol(:,:)            ! model state vector (mixed units)
  real(rkind),intent(out),device         :: relTol(:,:)            ! model state vector (mixed units)
  integer(i4b),intent(out)        :: err                    ! error code
  character(*),intent(out)        :: message                ! error message
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! state subsets
  integer(i4b)                       :: iState                    ! index of state within the snow+soil domain
  integer(i4b)                       :: iLayer                    ! index of layer within the snow+soil domain
  integer(i4b)                       :: ixStateSubset             ! index within the state subset
  ! real(rkind)                        :: absTolTempCas
  ! real(rkind)                        :: relTolTempCas
  ! real(rkind)                        :: absTolTempVeg
  ! real(rkind)                        :: relTolTempVeg
  ! real(rkind)                        :: absTolWatVeg
  ! real(rkind)                        :: relTolWatVeg
  ! real(rkind)                        :: absTolTempSoilSnow
  ! real(rkind)                        :: relTolTempSoilSnow
  ! real(rkind)                        :: absTolWatSnow
  ! real(rkind)                        :: relTolWatSnow
  ! real(rkind)                        :: absTolMatric
  ! real(rkind)                        :: relTolMatric
  ! real(rkind)                        :: absTolAquifr
  ! real(rkind)                        :: relTolAquifr
        type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  ! --------------------------------------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------------------------------------
  ! make association with variables in the data structures
  fixedLength: associate(&
    scalarCanairTemp    => prog_data%scalarCanairTemp       ,& ! intent(in) : [dp]     temperature of the canopy air space (K)
    scalarCanopyTemp    => prog_data%scalarCanopyTemp       ,& ! intent(in) : [dp]     temperature of the vegetation canopy (K)
    scalarCanopyWat     => prog_data%scalarCanopyWat        ,& ! intent(in) : [dp]     mass of total water on the vegetation canopy (kg m-2)
    scalarCanopyLiq     => prog_data%scalarCanopyLiq        ,& ! intent(in) : [dp]     mass of liquid water on the vegetation canopy (kg m-2)
    ! model state variable vectors for the snow-soil layers
    mLayerTemp          => prog_data%mLayerTemp                ,& ! intent(in) : [dp(:)]  temperature of each snow/soil layer (K)
    mLayerVolFracWat    => prog_data%mLayerVolFracWat          ,& ! intent(in) : [dp(:)]  volumetric fraction of total water (-)
    mLayerVolFracLiq    => prog_data%mLayerVolFracLiq          ,& ! intent(in) : [dp(:)]  volumetric fraction of liquid water (-)
    mLayerMatricHead    => prog_data%mLayerMatricHead          ,& ! intent(in) : [dp(:)]  matric head (m)
    mLayerMatricHeadLiq => diag_data%mLayerMatricHeadLiq       ,& ! intent(in) : [dp(:)]  matric potential of liquid water (m)
    ! model state variables for the aquifer
    scalarAquiferStorage=> prog_data%scalarAquiferStorage   ,& ! intent(in) : [dp]     storage of water in the aquifer (m)
    ! indices defining specific model states
    ixCasNrg            => indx_data%ixCasNrg                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy air space energy state variable
    ixVegNrg            => indx_data%ixVegNrg                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy energy state variable
    ixVegHyd            => indx_data%ixVegHyd                 ,& ! intent(in) : [i4b(:)] [length=1] index of canopy hydrology state variable (mass)
    ixAqWat             => indx_data%ixAqWat                  ,& ! intent(in) : [i4b(:)] [length=1] index of aquifer storage state variable
    ! vector of energy and hydrology indices for the snow and soil domains
    ixSnowSoilNrg       => indx_data%ixSnowSoilNrg         ,& ! intent(in) : [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    ixSnowSoilHyd       => indx_data%ixSnowSoilHyd         ,& ! intent(in) : [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    ! nSnowSoilNrg        => indx_data%nSnowSoilNrg          ,& ! intent(in) : [i4b]    number of energy state variables in the snow+soil domain
    ! nSnowSoilHyd        => indx_data%nSnowSoilHyd          ,& ! intent(in) : [i4b]    number of hydrology state variables in the snow+soil domain
    ! type of model state variabless
    ixStateType_subset  => indx_data%ixStateType       ,& ! intent(in) : [i4b(:)] [state subset] type of desired model state variables
    ixHydType           => indx_data%ixHydType                ,& ! intent(in) : [i4b(:)] index of the type of hydrology states in snow+soil domain
    ! number of layers
    nSnow               => indx_data%nSnow                 ,& ! intent(in) : [i4b]    number of snow layers
    nSoil               => indx_data%nSoil                 ,& ! intent(in) : [i4b]    number of soil layers
    nLayers             => indx_data%nLayers_d                & ! intent(in) : [i4b]    total number of layers
    )  ! end association with variables in the data structures
    ! --------------------------------------------------------------------------------------------------------------------------------
    ! initialize error control
    err=0; message='popTol4ida/'

    ! absTolTempCas      = mpar_data%var(iLookPARAM%absTolTempCas)%dat(1)
    ! relTolTempCas      = mpar_data%var(iLookPARAM%relTolTempCas)%dat(1)
    ! absTolTempVeg      = mpar_data%var(iLookPARAM%absTolTempVeg)%dat(1)
    ! relTolTempVeg      = mpar_data%var(iLookPARAM%relTolTempVeg)%dat(1)
    ! absTolWatVeg       = mpar_data%var(iLookPARAM%absTolWatVeg)%dat(1)
    ! relTolWatVeg       = mpar_data%var(iLookPARAM%relTolWatVeg)%dat(1)
    ! absTolTempSoilSnow = mpar_data%var(iLookPARAM%absTolTempSoilSnow)%dat(1)
    ! relTolTempSoilSnow = mpar_data%var(iLookPARAM%relTolTempSoilSnow)%dat(1)
    ! absTolWatSnow      = mpar_data%var(iLookPARAM%absTolWatSnow)%dat(1)
    ! relTolWatSnow      = mpar_data%var(iLookPARAM%relTolWatSnow)%dat(1)
    ! absTolMatric       = mpar_data%var(iLookPARAM%absTolMatric)%dat(1)
    ! relTolMatric       = mpar_data%var(iLookPARAM%relTolMatric)%dat(1)
    ! absTolAquifr       = mpar_data%var(iLookPARAM%absTolAquifr)%dat(1)
    ! relTolAquifr       = mpar_data%var(iLookPARAM%relTolAquifr)%dat(1)
 
    ! -----
    ! * initialize state vectors...
    ! -----------------------------

    ! initialize flags
    call popTol4ida_kernel<<<blocks,threads>>>(nGRU,ixCasNrg,ixVegNrg,ixVegHyd,ixAqWat,nLayers,&
  ixSnowSoilNrg,ixSnowSoilHyd,ixHydType,ixStateType_subset, &
  mpar_data%absTolTempCas,mpar_data%relTolTempCas,mpar_data%absTolTempVeg,mpar_data%relTolTempVeg,&
  mpar_data%absTolWatVeg,mpar_data%relTolWatVeg,mpar_data%absTolTempSoilSnow,mpar_data%relTolTempSoilSnow, &
  mpar_data%absTolWatSnow,mpar_data%relTolWatSnow,mpar_data%absTolMatric,mpar_data%relTolMatric,mpar_data%absTolAquifr,mpar_data%relTolAquifr,&
  absTol,relTol)

  end associate fixedLength      ! end association to variables in the data structure where vector length does not change
end subroutine popTol4ida

attributes(global) subroutine popTol4ida_kernel(nGRU,ixCasNrg,ixVegNrg,ixVegHyd,ixAqWat,nLayers,&
  ixSnowSoilNrg,ixSnowSoilHyd,ixHydType,ixStateType_subset, &
  absTolTempCas,relTolTempCas,absTolTempVeg,relTolTempVeg,&
  absTolWatVeg,relTolWatVeg,absTolTempSoilSnow,relTolTempSoilSnow, &
  absTolWatSnow,relTolWatSnow,absTolMatric,relTolMatric,absTolAquifr,relTolAquifr,&
  absTol,relTol)
  integer(i4b),value :: nGRU
  integer(i4b),intent(in) :: ixCasNrg(:), ixVegNrg(:), ixVegHyd(:), ixAqWat(:), nLayers(:)
  integer(i4b),intent(in) :: ixSnowSoilNrg(:,:), ixSnowSoilHyd(:,:), ixHydType(:,:), ixStateType_subset(:,:)
  real(rkind),intent(in) :: absTolTempCas,relTolTempCas,absTolTempVeg,relTolTempVeg
  real(rkind),intent(in) :: absTolWatVeg,relTolWatVeg,absTolTempSoilSnow,relTolTempSoilSnow
  real(rkind),intent(in) :: absTolWatSnow,relTolWatSnow,absTolMatric,relTolMatric,absTolAquifr,relTolAquifr
  real(rkind),intent(inout) :: absTol(:,:), relTol(:,:)

  integer(i4b) :: iGRU
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x

  if (iGRU .gt. nGRU) return

  call popTol4ida_device(ixCasNrg(iGRU),ixVegNrg(iGRU),ixVegHyd(iGRU),ixAqWat(iGRU),nLayers(iGRU),&
  ixSnowSoilNrg(:,iGRU),ixSnowSoilHyd(:,iGRU),ixHydType(:,iGRU),ixStateType_subset(:,iGRU), &
  absTolTempCas,relTolTempCas,absTolTempVeg,relTolTempVeg,&
  absTolWatVeg,relTolWatVeg,absTolTempSoilSnow,relTolTempSoilSnow, &
  absTolWatSnow,relTolWatSnow,absTolMatric,relTolMatric,absTolAquifr,relTolAquifr,&
  absTol(:,iGRU),relTol(:,iGRU))

end subroutine


attributes(device) subroutine popTol4ida_device(ixCasNrg,ixVegNrg,ixVegHyd,ixAqWat,nLayers,&
  ixSnowSoilNrg,ixSnowSoilHyd,ixHydType,ixStateType_subset, &
  absTolTempCas,relTolTempCas,absTolTempVeg,relTolTempVeg,&
  absTolWatVeg,relTolWatVeg,absTolTempSoilSnow,relTolTempSoilSnow, &
  absTolWatSnow,relTolWatSnow,absTolMatric,relTolMatric,absTolAquifr,relTolAquifr,&
  absTol,relTol)
  integer(i4b),intent(in) :: ixCasNrg, ixVegNrg, ixVegHyd, ixAqWat, nLayers
  integer(i4b),intent(in) :: ixSnowSoilNrg(:), ixSnowSoilHyd(:), ixHydType(:), ixStateType_subset(:)
  real(rkind),intent(in) :: absTolTempCas,relTolTempCas,absTolTempVeg,relTolTempVeg
  real(rkind),intent(in) :: absTolWatVeg,relTolWatVeg,absTolTempSoilSnow,relTolTempSoilSnow
  real(rkind),intent(in) :: absTolWatSnow,relTolWatSnow,absTolMatric,relTolMatric,absTolAquifr,relTolAquifr
  real(rkind),intent(inout) :: absTol(:), relTol(:)

  integer(i4b) :: iLayer, ixStateSubset
    if (ixCasNrg/=integerMissing) then
      absTol( ixCasNrg )  = absTolTempCas            ! transfer canopy air temperature to the state vector
      relTol( ixCasNrg )  = relTolTempCas
    end if

    if (ixVegNrg/=integerMissing) then
      absTol( ixVegNrg )  = absTolTempVeg      ! transfer vegetation temperature to the state vector
      relTol( ixVegNrg )  = relTolTempVeg     ! transfer vegetation temperature to the state vector
    end if

    if (ixVegHyd/=integerMissing) then
      select case(ixStateType_subset( ixVegHyd ))
        case(iname_watCanopy); absTol( ixVegHyd )  = absTolWatVeg ; relTol( ixVegHyd )  = relTolWatVeg
        case(iname_liqCanopy); absTol( ixVegHyd )  = absTolWatVeg ; relTol( ixVegHyd )  = relTolWatVeg       ! transfer liquid canopy water to the state vector
      end select
    end if

    ! tolerance for tempreture of the snow and soil domain
    do iLayer=1,nLayers
      if (ixSnowSoilNrg(iLayer)/=integerMissing) then   ! (loop through non-missing energy state variables in the snow+soil domain)
        ixStateSubset            = ixSnowSoilNrg(iLayer)  ! index within the state vector
        absTol(ixStateSubset)  = absTolTempSoilSnow    ! transfer temperature from a layer to the state vector
        relTol(ixStateSubset)  = relTolTempSoilSnow
      end if
    end do  ! looping through non-missing energy state variables in the snow+soil domain

    ! NOTE: ixVolFracWat  and ixVolFracLiq can also include states in the soil domain, hence enable primary variable switching
    do iLayer=1,nLayers
      if (ixSnowSoilHyd(iLayer)/=integerMissing) then   ! (loop through non-missing hydrology state variables in the snow+soil domain)
        ixStateSubset            = ixSnowSoilHyd(iLayer)   ! index within the state vector
        select case( ixHydType(iLayer) )
          case(iname_watLayer); absTol(ixStateSubset) = absTolWatSnow ;  relTol(ixStateSubset) = relTolWatSnow
          case(iname_liqLayer); absTol(ixStateSubset) = absTolWatSnow ;  relTol(ixStateSubset) = relTolWatSnow
          case(iname_matLayer); absTol(ixStateSubset) = absTolMatric ;  relTol(ixStateSubset) = relTolMatric
          case(iname_lmpLayer); absTol(ixStateSubset) = absTolMatric ;  relTol(ixStateSubset) = relTolMatric
        end select
      end if
    end do  ! looping through non-missing energy state variables in the snow+soil domain

    ! build the state vector for the aquifer storage
    ! NOTE: currently vector length=1, and use "do concurrent" to generalize to a multi-layer aquifer
    if (ixAqWat/=integerMissing) then
      absTol( ixAqWat )  = absTolAquifr    ! transfer aquifer storage to the state vector
      relTol( ixAqWat )  = relTolAquifr
    end if

  end subroutine


end module tol4ida_module
