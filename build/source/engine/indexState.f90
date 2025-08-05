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

module indexState_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:var_ilength                            ! data vector with variable length dimension (i4b)
USE data_types,only:in_type_indexSplit,out_type_indexSplit ! classes for indexSplit subroutine arguments

! missing data
USE globalData,only:integerMissing  ! missing integer

! named variables for domain types
USE globalData,only:iname_cas       ! canopy air space
USE globalData,only:iname_veg       ! vegetation
USE globalData,only:iname_snow      ! snow
USE globalData,only:iname_soil      ! soil
USE globalData,only:iname_aquifer   ! aquifer

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

! metadata
USE globalData,only:indx_meta       ! metadata for the variables in the index structure

! indices that define elements of the data structures
USE var_lookup,only:iLookINDEX      ! named variables for structure elements

! privacy
implicit none
private
public::indexState
public::indexSplit
public::indexSplit_d
public::indxSubset
contains


 ! **********************************************************************************************************
 ! public subroutine indexState: define list of indices for each state variable
 ! **********************************************************************************************************
 subroutine indexState(computeVegFlux,          & ! intent(in):    flag to denote if computing the vegetation flux
                       includeAquifer,          & ! intent(in):    flag to denote if an aquifer is included
                       nSnow,nSoil,nLayers,nGRU,     & ! intent(in):    number of snow and soil layers, and total number of layers
                       indx_data,               & ! intent(inout): indices defining model states and layers
                       err,message)               ! intent(out):   error control
 ! provide access to the numerical recipes utility modules
 USE nr_utility_module,only:arth                           ! creates a sequence of numbers (start, incr, n)
 use device_data_types
 implicit none
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input
 logical(lgt),intent(in),device         :: computeVegFlux(:)         ! flag to denote if computing the vegetation flux
 logical(lgt),intent(in)         :: includeAquifer         ! flag to denote if an aquifer is included
 integer(i4b),intent(in)         :: nSoil,nGRU    ! number of snow and soil layers, and total number of layers
 integer(i4b),intent(in),device :: nSnow(:),nLayers(:)
 type(indx_data_device),intent(inout) :: indx_data              ! indices defining model states and layers
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! general local variables
 character(len=256)              :: cmessage               ! message of downwind routine
 integer(i4b),parameter          :: nVarSnowSoil=2         ! number of state variables in the snow and soil domain (energy and total water/matric head)
 integer(i4b)                    :: nAquiferState          ! number of aquifer state variables
 ! indices of model state variables
!  integer(i4b)                    :: ixTopNrg               ! index of upper-most energy state in the snow-soil subdomain
!  integer(i4b)                    :: ixTopWat               ! index of upper-most total water state in the snow-soil subdomain
 integer(i4b),device :: subsetLayer(nGRU)
 integer(i4b) :: iGRU, ilayer, maxnState
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! make association with variables in the data structures
 associate(&
 ! number of state variables of different type
 nCasNrg       => indx_data%nCasNrg   , & ! number of energy state variables for the canopy air space
 nVegNrg       => indx_data%nVegNrg   , & ! number of energy state variables for the vegetation canopy
 ixTopWat => indx_data%ixTopHyd, &
 ixTopNrg => indx_data%ixTopNrg, &
 nVegMass      => indx_data%nVegMass  , & ! number of hydrology states for vegetation (mass of water)
 nVegState     => indx_data%nVegState , & ! number of vegetation state variables
 nNrgState     => indx_data%nNrgState , & ! number of energy state variables
 nWatState     => indx_data%nWatState , & ! number of "total water" states (vol. total water content)
 nMatState     => indx_data%nMatState , & ! number of matric head state variables
 nMassState    => indx_data%nMassState, & ! number of hydrology state variables (mass of water)
 nState        => indx_data%nState    , & ! total number of model state variables
 ! vectors of indices for specfic state types within specific sub-domains IN THE FULL STATE VECTOR
 ixNrgCanair   => indx_data%ixNrgCanair  , & ! indices IN THE FULL VECTOR for energy states in canopy air space domain
 ixNrgCanopy   => indx_data%ixNrgCanopy  , & ! indices IN THE FULL VECTOR for energy states in the canopy domain
 ixHydCanopy   => indx_data%ixHydCanopy  , & ! indices IN THE FULL VECTOR for hydrology states in the canopy domain
 ixNrgLayer    => indx_data%ixNrgLayer   , & ! indices IN THE FULL VECTOR for energy states in the snow+soil domain
 ixHydLayer    => indx_data%ixHydLayer   , & ! indices IN THE FULL VECTOR for hyd states in the snow+soil domain
 ixWatAquifer  => indx_data%ixWatAquifer , & ! indices IN THE FULL VECTOR for the aquifer
 ! indices for model state variables
 ixSoilState   => indx_data%ixSoilState  , & ! list of indices for all soil layers
 ixLayerState  => indx_data%ixLayerState   & ! list of indices for all model layers
 ) ! association to variables in the data structures
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='indexState/'

 ! -----
 ! * define the number of state variables...
 ! -----------------------------------------

 ! define the number of vegetation state variables (defines position of snow-soil states in the state vector)
 !$cuf kernel do(1) <<<*,*>>>
 do iGRU=1,nGRU
 if(computeVegFlux(iGRU))then
  nCasNrg(iGRU)   = 1
  nVegNrg(iGRU)   = 1
  nVegMass(iGRU)  = 1
  nVegState(iGRU) = nCasNrg(iGRU) + nVegNrg(iGRU) + nVegMass(iGRU)
 else
  nCasNrg(iGRU)   = 0
  nVegNrg(iGRU)   = 0
  nVegMass(iGRU)  = 0
  nVegState(iGRU) = 0
 end if

 ! define the number of aquifer states
 nAquiferState = merge(1,0,includeAquifer)

 ! define the number state variables of different type
 nNrgState(iGRU)  = nCasNrg(iGRU) + nVegNrg(iGRU) + nLayers(iGRU)  ! number of energy state variables
 nWatState(iGRU)  = nSnow(iGRU)                        ! number of "total water" state variables -- will be modified later if using primary variable switching
 nMatState(iGRU)  = nSoil                        ! number of matric head state variables -- will be modified later if using primary variable switching
 nMassState(iGRU) = nVegMass(iGRU)                     ! number of mass state variables -- currently restricted to canopy water

 ! define the number of model state variables
 nState(iGRU) = nVegState(iGRU) + nLayers(iGRU)*nVarSnowSoil + nAquiferState  ! *nVarSnowSoil (both energy and total water)

 ! -----
 ! * define the indices of state variables WITHIN THE FULL STATE VECTOR...
 ! -----------------------------------------------------------------------

 ! define indices in the vegetation domain
 if(computeVegFlux(iGRU))then
  ixNrgCanair(iGRU) = 1 ! indices IN THE FULL VECTOR for energy states in canopy air space domain  (-)
  ixNrgCanopy(iGRU) = 2 ! indices IN THE FULL VECTOR for energy states in the canopy domain        (-)
  ixHydCanopy(iGRU) = 3 ! indices IN THE FULL VECTOR for hydrology states in the canopy domain     (-)
 else
  ixNrgCanair(iGRU) = integerMissing
  ixNrgCanopy(iGRU) = integerMissing
  ixHydCanopy(iGRU) = integerMissing
 end if

 ! define the index of the top layer
 ! NOTE: local variables -- actual indices defined when building the state subset
 ixTopNrg(iGRU) = nVegState(iGRU) + 1                       ! energy
 ixTopWat(iGRU) = nVegState(iGRU) + 2                       ! total water (only snow)

 ! define the indices within the snow+soil domain
  subsetLayer(iGRU) = ixTopNrg(iGRU)
  do iLayer=1,nLayers(iGRU)
    ixNrgLayer(iLayer,iGRU) = subsetLayer(iGRU)
    subsetLayer(iGRU) = subsetLayer(iGRU) + nVarSnowSoil
  enddo
  subsetLayer(iGRU) = ixTopWat(iGRU)
  do iLayer=1,nLayers(iGRU)
    ixHydLayer(iLayer,iGRU) = subsetLayer(iGRU)
    subsetLayer(iGRU) = subsetLayer(iGRU) + nVarSnowSoil
  enddo

 ! define indices for the aquifer
 ixWatAquifer(iGRU) = merge(nState(iGRU), integerMissing, includeAquifer)

end do

 ! -----
 ! * define the type of model states...
 ! ------------------------------------

 ! re-allocate index vectors for the full state vector (if needed)...
 maxnState = maxval(nState)
 if (size(indx_data%ixMapFull2Subset,1).ne.maxnState) then
  deallocate(indx_data%ixMapFull2Subset)
  allocate(indx_data%ixMapFull2Subset(maxnState,nGRU))
  indx_data%ixMapFull2Subset = integerMissing
 end if
 if (.not. allocated(indx_data%ixMapFull2Subset)) then
  allocate(indx_data%ixMapFull2Subset(maxnState,nGRU))
  indx_data%ixMapFull2Subset = integerMissing
 end if

 if (size(indx_data%ixControlVolume,1).ne.maxnState) then
  deallocate(indx_data%ixControlVolume)
  allocate(indx_data%ixControlVolume(maxnState,nGRU))
  indx_data%ixControlVolume = integerMissing
 end if
 if (.not. allocated(indx_data%ixControlVolume)) then
  allocate(indx_data%ixControlVolume(maxnState,nGRU))
  indx_data%ixControlVolume = integerMissing
 end if

 if (size(indx_data%ixDomainType,1).ne.maxnState) then
  deallocate(indx_data%ixDomainType)
  allocate(indx_data%ixDomainType(maxnState,nGRU))
  indx_data%ixDomainType = integerMissing
 end if
 if (.not. allocated(indx_data%ixDomainType)) then
  allocate(indx_data%ixDomainType(maxnState,nGRU))
  indx_data%ixDomainType = integerMissing
 end if

 if (size(indx_data%ixStateType,1).ne.maxnState) then
  deallocate(indx_data%ixStateType)
  allocate(indx_data%ixStateType(maxnState,nGRU))
  indx_data%ixStateType = integerMissing
 end if
 if (.not. allocated(indx_data%ixStateType)) then
  allocate(indx_data%ixStateType(maxnState,nGRU))
  indx_data%ixStateType = integerMissing
 end if

 if (size(indx_data%ixAllState,1).ne.maxnState) then
  deallocate(indx_data%ixAllState)
  allocate(indx_data%ixAllState(maxnState,nGRU))
  indx_data%ixAllState = integerMissing
 end if
 if (.not. allocated(indx_data%ixAllState)) then
  allocate(indx_data%ixAllState(maxnState,nGRU))
  indx_data%ixALlState = integerMissing
 end if


 ! make an association to the ALLOCATABLE variables in the data structures
 ! NOTE: we need to do this here since the size may have changed above
 associate(&
 ixControlVolume => indx_data%ixControlVolume , & ! index of control volume for different domains (veg, snow, soil)
 ixDomainType    => indx_data%ixDomainType    , & ! indices defining the type of the domain (iname_veg, iname_snow, iname_soil)
 ixStateType     => indx_data%ixStateType     , & ! indices defining the type of the state (iname_nrgLayer...)
 ixAllState      => indx_data%ixAllState        & ! list of indices for all model state variables
 )  ! making an association to variables in the data structures

 ! define indices for state variables
 !$cuf kernel do(1) <<<*,*>>>
 do iGRU=1,nGRU
  subsetLayer(iGRU) = 1
  do iLayer=1,nState(iGRU)
    ixAllState(iLayer,iGRU) = subsetLayer(iGRU)
    subsetLayer(iGRU) = subsetLayer(iGRU) + 1
  end do
  subsetLayer(iGRU) = 1
  do iLayer=1,nSoil
    ixSoilState(iLayer,iGRU) = subsetLayer(iGRU)
    subsetLayer(iGRU) = subsetLayer(iGRU) + 1
  end do
  subsetLayer(iGRU) = 1
  do iLayer=1,nLayers(iGRU)
    ixLayerState(iLayer,iGRU) = subsetLayer(iGRU)
    subsetLayer(iGRU) = subsetLayer(iGRU) + 1
  end do


 ! define the state type for the vegetation canopy
 if(computeVegFlux(iGRU))then
  ixStateType(ixNrgCanair(iGRU),iGRU) = iname_nrgCanair
  ixStateType(ixNrgCanopy(iGRU),iGRU) = iname_nrgCanopy
  ixStateType(ixHydCanopy(iGRU),iGRU) = iname_watCanopy
 endif

 do iLayer=1,size(ixNrgLayer,1)
 ! define the state type for the snow+soil domain (energy)
 ixStateType(ixNrgLayer(iLayer,iGRU),iGRU) = iname_nrgLayer
 end do

 do iLayer=1,nLayers(iGRU)
 ! define the state type for the snow+soil domain (hydrology)
 if(iLayer .le. nSnow(iGRU)) ixStateType( ixHydLayer(iLayer,iGRU),iGRU   ) = iname_watLayer
 if(iLayer .gt. nSnow(iGRU)) ixStateType( ixHydLayer(iLayer,iGRU),iGRU ) = iname_matLayer ! refine later to be either iname_watLayer or iname_matLayer
 end do
 ! define the state type for the aquifer
 if(includeAquifer) ixStateType( ixWatAquifer(iGRU),iGRU ) = iname_watAquifer
 ! define the domain type for vegetation
 if(computeVegFlux(iGRU))then
  ixDomainType(ixNrgCanair(iGRU),iGRU) = iname_cas
  ixDomainType(ixNrgCanopy(iGRU),iGRU) = iname_veg
  ixDomainType(ixHydCanopy(iGRU),iGRU) = iname_veg
 endif

 do iLayer=1,nLayers(iGRU)
 ! define the domain type for snow
 if(iLayer .le. nSnow(iGRU))then
  ixDomainType( ixNrgLayer(iLayer,iGRU),iGRU ) = iname_snow
  ixDomainType( ixHydLayer(iLayer,iGRU),iGRU ) = iname_snow
 else

 ! define the domain type for soil
 ixDomainType( ixNrgLayer(iLayer,iGRU),iGRU ) = iname_soil
 ixDomainType( ixHydLayer(iLayer,iGRU),iGRU ) = iname_soil
 end if
end do

 ! define the domain type for the aquifer
 if(includeAquifer) ixDomainType( ixWatAquifer(iGRU),iGRU ) = iname_aquifer

 ! define the index of each control volume in the vegetation domains
 if(computeVegFlux(iGRU))then
  ixControlVolume(ixNrgCanair(iGRU),iGRU) = 1  ! NOTE: assumes scalar
  ixControlVolume(ixNrgCanopy(iGRU),iGRU) = 1
  ixControlVolume(ixHydCanopy(iGRU),iGRU) = 1
 endif

 do iLayer=1,nLayers(iGRU)
 ! define the index of the each control volume in the snow domain
 if(iLayer .le. nSnow(iGRU))then
  ixControlVolume( ixNrgLayer(iLayer,iGRU),iGRU ) = ixLayerState(iLayer,iGRU)
  ixControlVolume( ixHydLayer(iLayer,iGRU),iGRU ) = ixLayerState(iLayer,iGRU)
 else

 ! define the index of the each control volume in the soil domain
 ixControlVolume( ixNrgLayer(iLayer,iGRU),iGRU ) = ixSoilState(iLayer-nSnow(iGRU),iGRU)
 ixControlVolume( ixHydLayer(iLayer,iGRU),iGRU ) = ixSoilState(iLayer-nSnow(iGRU),iGRU)
 end if
end do

 ! define the index for the control volumes in the aquifer
 if(includeAquifer) ixControlVolume( ixWatAquifer(iGRU),iGRU ) = 1
end do
 !print*, 'ixControlVolume = ', ixControlVolume
 !print*, 'ixDomainType    = ', ixDomainType
 !print*, 'ixStateType     = ', ixStateType
 !print*, 'PAUSE: '; read(*,*)

 ! end association to the ALLOCATABLE variables in the data structures
 end associate

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------------------------------

 end associate  ! end association to variables in the data structures
 end subroutine indexState

 ! **********************************************************************************************************
 ! public subroutine indexSplit: define list of indices for each state variable
 ! **********************************************************************************************************
  subroutine indexSplit(in_indexSplit,               & ! intent(in)    : number of model layers and states in a subset
    nGRU, &
                       stateSubsetMask,             & ! intent(in)    : logical vector (.true. if state is in the subset)
                       indx_data,                   & ! intent(inout) : index data structure
                       out_indexSplit)                ! intent(out)   : error control
!  ! external modules 
!  USE f2008funcs_module,only:findIndex                 ! finds the index of the first value within a vector
!  USE nr_utility_module,only:arth                      ! creates a sequence of numbers (start, incr, n)
 use device_data_types
 use globalData,only:indx_meta
 implicit none
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! input
 type(in_type_indexSplit),intent(in)   :: in_indexSplit      ! number of model layers and states in a subset
 integer(i4b) :: nGRU
 logical(lgt),intent(in),device               :: stateSubsetMask(:,:) ! logical vector (.true. if state is in the subset)
 type(indx_data_device),intent(inout)       :: indx_data          ! indices defining model states and layers
 ! output
 type(out_type_indexSplit),intent(out) :: out_indexSplit     ! error control
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                                    :: iVar            ! variable index
 integer(i4b)                                    :: ixVegWat        ! index of total water in the vegetation canopy
 integer(i4b)                                    :: ixVegLiq        ! index of liquid water in the vegetation canopy
 integer(i4b)                                    :: ixTopWat        ! index of upper-most total water state in the snow-soil subdomain
 integer(i4b)                                    :: ixTopLiq        ! index of upper-most liquid water state in the snow-soil subdomain
 integer(i4b)                                    :: ixTopMat        ! index of upper-most total water matric potential state in the soil subdomain
 integer(i4b)                                    :: ixTopLMP        ! index of upper-most liquid water matric potential state in the soil subdomain
 integer(i4b),dimension(in_indexSplit % nSubset) :: ixSequence      ! sequential index in model state vector
 logical(lgt),dimension(in_indexSplit % nSubset,nGRU),device :: stateTypeMask   ! mask of state vector for specific state subsets
 logical(lgt),dimension(in_indexSplit % nLayers,nGRU),device :: volFracWat_mask ! mask of layers within the snow+soil domain
 logical(lgt),dimension(in_indexSplit % nSoil,nGRU),device   :: matricHead_mask ! mask of layers within the soil domain
 character(len=256)                              :: cmessage        ! error message of downwind routine
 integer(i4b) :: iGRU, iLayer
 integer(i4b) :: nSoil
 integer(i4b),device :: subsetCount(nGRU), subsetLayer(nGRU)
 integer(i4b) :: subset, matricHeadCount

 ! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! make association to variables in the data structures
 fullState: associate(&
 ! number of snow and soil layers, total number of layers, and number of states in the subset
!  nSnow            => in_indexSplit % nSnow                          ,& ! intent(in):  [i4b]    number of snow layers 
!  nSoil            => in_indexSplit % nSoil                          ,& ! intent(in):  [i4b]    number of soil layers 
!  nLayers          => in_indexSplit % nLayers                        ,& ! intent(in):  [i4b]    total number of layers
 nSubset          => in_indexSplit % nSubset                        ,& ! intent(in):  [i4b]    number of states in the subset
 ! indices of model state variables for the vegetation domain
 ixCasNrg         => indx_data%ixCasNrg      ,& ! intent(in):  [i4b]    index of canopy air space energy state variable
 ixVegNrg         => indx_data%ixVegNrg      ,& ! intent(in):  [i4b]    index of canopy energy state variable
 ixVegHyd         => indx_data%ixVegHyd      ,& ! intent(in):  [i4b]    index of canopy hydrology state variable (mass)

 ! indices of the top model state variables in the snow+soil system
 ixTopNrg         => indx_data%ixTopNrg      ,& ! intent(in):  [i4b]    index of upper-most energy state in the snow-soil subdomain
 ixTopHyd         => indx_data%ixTopHyd      ,& ! intent(in):  [i4b]    index of upper-most hydrology state in the snow-soil subdomain
 nSnow_d => indx_data%nSnow, &

 ! index of the storage of water in the aquifer
 ixAqWat          => indx_data%ixAqWat       ,& ! intent(in):  [i4b]    index of the storage of water in the aquifer

 ! indices of model state variables
 ixMapFull2Subset => indx_data%ixMapFull2Subset ,& ! intent(in):  [i4b(:)] list of indices in the state subset (missing for values not in the subset)
 ixDomainType     => indx_data%ixDomainType     ,& ! intent(in):  [i4b(:)] indices defining the domain of the state (iname_veg, iname_snow, iname_soil)
 ixStateType      => indx_data%ixStateType      ,& ! intent(in):  [i4b(:)] indices defining the type of the state (ixNrgState...)
!  ixAllState       => indx_data%ixAllState       ,& ! intent(in):  [i4b(:)] list of indices for all model state variables (1,2,3,...nState)
 ixNrgLayer       => indx_data%ixNrgLayer       ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for energy states in the snow+soil domain
 ixHydLayer       => indx_data%ixHydLayer       ,& ! intent(in):  [i4b(:)] indices IN THE FULL VECTOR for hydrology states in the snow+soil domain
 ixHydType        => indx_data%ixHydType        ,& ! intent(in):  [i4b(:)] index of the type of hydrology states in snow+soil domain

 ! indices of the entire state vector, all model layers, and soil layers
 ixSoilState      => indx_data%ixSoilState      ,& ! intent(in):  [i4b(:)] list of indices for all soil layers
!  ixLayerState     => indx_data%ixLayerState     ,& ! intent(in):  [i4b(:)] list of indices for all model layers

!  ! vector of energy indices for the snow and soil domains
!  ! NOTE: states not in the subset are equal to integerMissing
 ixSnowSoilNrg    => indx_data%ixSnowSoilNrg    ,& ! intent(in):  [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
 ixSnowOnlyNrg    => indx_data%ixSnowOnlyNrg    ,& ! intent(in):  [i4b(:)] index in the state subset for energy state variables in the snow domain
 ixSoilOnlyNrg    => indx_data%ixSoilOnlyNrg    ,& ! intent(in):  [i4b(:)] index in the state subset for energy state variables in the soil domain

 ! vector of hydrology indices for the snow and soil domains
 ! NOTE: states not in the subset are equal to integerMissing
 ixSnowSoilHyd    => indx_data%ixSnowSoilHyd    ,& ! intent(in):  [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
 ixSnowOnlyHyd    => indx_data%ixSnowOnlyHyd    ,& ! intent(in):  [i4b(:)] index in the state subset for hydrology state variables in the snow domain
 ixSoilOnlyHyd    => indx_data%ixSoilOnlyHyd    ,& ! intent(in):  [i4b(:)] index in the state subset for hydrology state variables in the soil domain

!  ! indices of active model layers
!  ixLayerActive    => indx_data%ixLayerActive    ,& ! intent(in):  [i4b(:)] index of active model layers (inactive=integerMissing)

 nLayers => indx_data%nLayers_d, &
 ! number of state variables of a specific type
!  nSnowSoilNrg     => indx_data%nLayers  ,& ! intent(in):  [i4b]    number of energy state variables in the snow+soil domain
!  nSnowOnlyNrg     => indx_data%nSnowOnlyNrg  ,& ! intent(in):  [i4b]    number of energy state variables in the snow domain
!  nSoilOnlyNrg     => indx_data%nSoilOnlyNrg  ,& ! intent(in):  [i4b]    number of energy state variables in the soil domain
!  nSnowSoilHyd     => indx_data%nLayers  ,& ! intent(in):  [i4b]    number of hydrology variables in the snow+soil domain
!  nSnowOnlyHyd     => indx_data%nSnowOnlyHyd  ,& ! intent(in):  [i4b]    number of hydrology variables in the snow domain
!  nSoilOnlyHyd     => indx_data%nSoilOnlyHyd  ,& ! intent(in):  [i4b]    number of hydrology variables in the soil domain

 ! error control
 err              => out_indexSplit % err                           ,& ! intent(out): [i4b]       error code
 message          => out_indexSplit % cmessage                       & ! intent(out): [character] error message
 ) ! association to variables in the data structures

!  nLayers = in_indexSplit % nLayers
 nSoil = in_indexSplit % nSoil
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='indexSplit/'

 ! -----
 ! - preliminaries...
 ! ------------------

 ! define the type of variable in the snow+soil domain
 !$cuf kernel do(1) <<<*,*>>>
 do iGRU=1,nGRU
  do iLayer=1,nLayers(iGRU)
    if (ixHydLayer(iLayer,iGRU).ne.integerMissing) ixHydType(iLayer,iGRU) = ixStateType(ixHydLayer(iLayer,iGRU),iGRU)
  end do
end do

 ! get the mapping between the full state vector and the state subset

!$cuf kernel do(1) <<<*,*>>>
do iGRU=1,nGRU
  subsetLayer(iGRU) = 1
  do iLayer=1,size(ixMapFull2Subset,1)
    if (stateSubsetMask(iLayer,iGRU)) then
      ixMapFull2Subset(iLayer,iGRU) = subsetLayer(iGRU)
      subsetLayer(iGRU) = subsetLayer(iGRU) + 1

    else
      ixMapFull2Subset(iLayer,iGRU) = integerMissing
    endif
  enddo
enddo

 ! -----
 ! - get vectors of different state subsets...
 ! -------------------------------------------

 ! get different masks
!$cuf kernel do(1) <<<*,*>>>
do iGRU=1,nGRU
  do iLayer=1,nLayers(iGRU)
    volFracWat_mask(iLayer,iGRU) =  ixHydType(iLayer,iGRU)==iname_watLayer .or. ixHydType(iLayer,iGRU)==iname_liqLayer
  end do
end do
!$cuf kernel do(2) <<<*,*>>>
do iGRU=1,nGRU
  do iLayer=1,nSoil
    matricHead_mask(iLayer,iGRU) = ixHydType(nSnow_d(iGRU)+iLayer,iGRU)==iname_matLayer .or. ixHydType(nSnow_d(iGRU)+iLayer,iGRU) == iname_lmpLayer
  end do
end do

subset = 0
!$cuf kernel do(1) <<<*,*>>> reduce(max:subset,matricHeadCount)
do iGRU=1,nGRU
 subsetCount(iGRU) = 0
 do iLayer=1,size(stateSubsetMask,1)
   if (stateSubsetMask(iLayer,iGRU)) subsetcount(iGRU) = subsetcount(iGRU) + 1
 end do
 subset = max(subsetcount(iGRU),subset)
 subsetCount(iGRU) = 0
 do iLayer=1,size(matricHead_mask,1)
  if (matricHead_mask(iLayer,iGRU)) subsetCount(iGRU) = subsetCount(iGRU) + 1
 end do
 matricHeadCount = max(matricHeadCount,subsetCount(iGRU))
end do

 ! get state subsets for desired variables
 do iVar=1,size(indx_meta)   ! loop through index variables

  ! get the subset of indices
  ! NOTE: indxSubset(subset, fullVector, mask), provides subset of fullVector where mask==.true.
  select case(iVar)
   case(iLookINDEX%ixMapSubset2Full);     call indxSubset_d(indx_data%ixMapSubset2Full,stateSubsetMask,subset,nGRU,err,cmessage)
   case(iLookINDEX%ixStateType_subset);   call indxSubset_d2(indx_data%ixStateType_subset, ixStateType,  stateSubsetMask,subset,nGRU, err, cmessage)
   case(iLookINDEX%ixDomainType_subset);  call indxSubset_d2(indx_data%ixDomainType_subset, ixDomainType, stateSubsetMask,subset,nGRU, err, cmessage)
!    case(iLookINDEX%ixVolFracWat);         call indxSubset(indx_data%var(iVar)%dat, ixLayerState, volFracWat_mask, err, cmessage)
  !  case(iLookINDEX%ixMatricHead);         call indxSubset_d2(indx_data%ixMatricHead_m, ixSoilState,  matricHead_mask, matricHeadCount,nGRU, err, cmessage)
   case default; cycle ! only need to process the above variables
  end select  ! iVar
  if(err/=0)then; message=trim(message)//trim(cmessage)//'[varname='//trim(indx_meta(ivar)%varname)//']'; return; endif

 end do  ! looping through variables in the data structure

!  ! make association to variables in the data structures
 subsetState: associate(ixStateType_subset => indx_data%ixStateType_subset) ! named variables defining the states in the subset

 ! -----
 ! - get indices for the (currently) scalar states in the vegetation domain...
 ! ---------------------------------------------------------------------------

 !$cuf kernel do(1) <<<*,*>>>
 do iGRU=1,nGRU
  ixCasNrg(iGRU) = integerMissing
  ixVegNrg(iGRU) = integerMissing
  ixTopNrg(iGRU) = integerMissing
  ixVegWat = integerMissing
  ixVegLiq = integerMissing
  ixTopWat = integerMissing
  ixTopLiq = integerMissing
  ixTopMat = integerMissing
  ixTopLMP = integerMissing
  ixAqWat(iGRU) = integerMissing
  do iLayer=size(ixStateType_subset,1),1,-1
   if (ixStateType_subset(iLayer,iGRU) == iname_nrgCanair) ixCasNrg(iGRU) = iLayer
   if (ixStateType_subset(iLayer,iGRU) == iname_nrgCanopy) ixVegNrg(iGRU) = iLayer
   if (ixStateType_subset(iLayer,iGRU) == iname_watCanopy) ixVegWat = iLayer
   if (ixStateType_subset(iLayer,iGRU) == iname_liqCanopy) ixVegLiq = iLayer
   if (ixStateType_subset(iLayer,iGRU) == iname_nrgLayer) ixTopNrg(iGRU) = iLayer
   if (ixStateType_subset(iLayer,iGRU) == iname_watLayer) ixTopWat = iLayer
   if (ixStateType_subset(iLayer,iGRU) == iname_liqLayer) ixTopLiq = iLayer
   if (ixStateType_subset(iLayer,iGRU) == iname_matLayer) ixTopMat = iLayer
   if (ixStateType_subset(iLayer,iGRU) == iname_lmpLayer) ixTopLMP = iLayer
   if (ixStateType_subset(iLayer,iGRU) == iname_watAquifer) ixAqWat(iGRU) = iLayer

  end do
  ixVegHyd(iGRU) = merge(ixVegWat,ixVegLiq,ixVegWat/=integerMissing)
  if(ixTopWat==integerMissing .and. ixTopLiq==integerMissing)then
   ixTopHyd(iGRU) = merge(ixTopMat, ixTopLMP, ixTopMat/=integerMissing)      ! no water state, so upper-most hydrology state is the upper-most matric head state (if it exists)
  else
   ixTopHyd(iGRU) = merge(ixTopWat, ixTopLiq, ixTopWat/=integerMissing)      ! ixTopWat is used if it is not missing
  endif

 end do


 ! -----
 ! - get vector of indices within the state subset state variables of a given type...
 ! ----------------------------------------------------------------------------------

 ! list of indices for all energy states
 subset = 0
 !$cuf kernel do(1) <<<*,*>>> reduce(max:subset)
 do iGRU=1,nGRU
  subsetCount(iGRU) = 0
  do iLayer=1,size(ixStateType_subset,1)
    stateTypeMask(iLayer,iGRU) = (ixStateType_subset(iLayer,iGRU)==iname_nrgCanair .or. ixStateType_subset(iLayer,iGRU)==iname_nrgCanopy .or. ixStateType_subset(iLayer,iGRU)==iname_nrgLayer)
    if (stateTypeMask(iLayer,iGRU)) subsetcount(iGRU) = subsetcount(iGRU) + 1
  end do
  subset = max(subsetcount(iGRU),subset)
end do
 call indxSubset_d(indx_data%ixNrgOnly,stateTypeMask,subset,nGRU,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage)//'[varname='//trim(indx_meta(iLookINDEX%ixNrgOnly)%varname)//']'; return; endif

 ! -----
 ! - get vector of indices of the state subset for layers in the snow+soil domain...
 ! ---------------------------------------------------------------------------------

 ! get list of indices for energy
 ! NOTE: layers not in the state subset will be missing
   !$cuf kernel do(1) <<<*,*>>>
   do iGRU=1,nGRU
    do iLayer=1,nLayers(iGRU)
      ixSnowSoilNrg(iLayer,iGRU) = ixMapFull2Subset(ixNrgLayer(iLayer,iGRU),iGRU)
      if (iLayer .le. nSnow_d(iGRU)) ixSnowOnlyNrg(iLayer,iGRU) = ixMapFull2Subset(ixNrgLayer(iLayer,iGRU),iGRU)
      if (iLayer .le. nSoil) ixSoilOnlyNrg(iLayer,iGRU) = ixMapFull2Subset(ixNrgLayer(iLayer+nSnow_d(iGRU),iGRU),iGRU)
    end do
  end do

  

 ! get list of indices for hydrology
 ! NOTE: layers not in the state subset will be missing
   !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    do iLayer=1,nLayers(iGRU)
      ixSnowSoilHyd(iLayer,iGRU) = ixMapFull2Subset(ixHydLayer(iLayer,iGRU),iGRU)
      if (iLayer .le. nSnow_d(iGRU)) ixSnowOnlyHyd(iLayer,iGRU) = ixMapFull2Subset(ixHydLayer(iLayer,iGRU),iGRU)
      if (iLayer .le. nSoil) ixSoilOnlyHyd(iLayer,iGRU) = ixMapFull2Subset(ixHydLayer(iLayer+nSnow_d(iGRU),iGRU),iGRU)
    end do
  end do

!  ! define active layers (regardless if the splitting operation is energy or mass)
!  ixLayerActive =  merge(ixSnowSoilNrg, ixSnowSoilHyd, ixSnowSoilNrg/=integerMissing)

 ! get the number of valid states for energy
  ! nSnowSoilNrg = maxval(nSnow_d) + nSoil
!  nSnowSoilNrg = count(ixSnowSoilNrg/=integerMissing)
!  nSnowOnlyNrg = count(ixSnowOnlyNrg/=integerMissing)
!  nSoilOnlyNrg = count(ixSoilOnlyNrg/=integerMissing)

 ! get the number of valid states for hydrology
  ! nSnowSoilHyd = maxval(nSnow_d) + nSoil
!  nSnowSoilHyd = count(ixSnowSoilHyd/=integerMissing)
!  nSnowOnlyHyd = count(ixSnowOnlyHyd/=integerMissing)
!  nSoilOnlyHyd = count(ixSoilOnlyHyd/=integerMissing)

!  ! end association to data in structures
 end associate subsetState
 end associate fullState

 end subroutine indexSplit


 ! **********************************************************************************************************
 ! public subroutine indxSubset: get a subset of indices for a given mask
 ! **********************************************************************************************************
 subroutine indxSubset(ixSubset,ixMaster,mask,err,message)
 implicit none
 ! input-output: subset of indices for allocation/population
 integer(i4b),intent(inout),allocatable :: ixSubset(:) ! subset of indices
 ! input
 integer(i4b),intent(in)                :: ixMaster(:) ! full list of indices
 logical(lgt),intent(in)                :: mask(:)     ! desired indices
 ! error control
 integer(i4b),intent(out)               :: err         ! error code
 character(*),intent(out)               :: message     ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                           :: nSubset     ! length of the subset
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! initialize errors
 err=0; message="indxSubset/"

 ! check size match
 if(size(ixMaster)/=size(mask))then
  message=trim(message)//'size mismatch'
  err=20; return
 endif

 ! get the number of variables
 nSubset = count(mask)

 ! check if we need to reallocate space
 if(size(ixSubset)/=nSubset) then

  ! deallocate space
  if(allocated(ixSubset))then
   deallocate(ixSubset,stat=err)
   if(err/=0)then; message=trim(message)//'unable to deallocate space for variable'; err=20; return; endif
  endif

  ! allocate space
  allocate(ixSubset(nSubset),stat=err)
  if(err/=0)then; message=trim(message)//'unable to deallocate space for variable'; err=20; return; endif

 endif  ! allocating space

 ! define indices for variable types in specific sub-domains
 if(nSubset>0) ixSubset = pack(ixMaster, mask)

 end subroutine indxSubset

  subroutine indxSubset_d(ixSubset,mask,nSubset,nGRU,err,message)
 implicit none
 ! input-output: subset of indices for allocation/population
 integer(i4b),intent(inout),allocatable,device :: ixSubset(:,:) ! subset of indices
 ! input
 logical(lgt),intent(in),device                :: mask(:,:)     ! desired indices
 integer(i4b) :: nGRU
 ! error control
 integer(i4b),intent(out)               :: err         ! error code
 character(*),intent(out)               :: message     ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                           :: nSubset     ! length of the subset
 integer(i4b),device :: subsetLayer(nGRU)
 integer(i4b) :: iGRU, iLayer
 ! -----------------------------------------------------------------------------------------------------------------------------------
 ! initialize errors
 err=0; message="indxSubset/"

 ! check size match
!  if(size(ixMaster)/=size(mask))then
!   message=trim(message)//'size mismatch'
!   err=20; return
!  endif

 ! check if we need to reallocate space
 if(size(ixSubset,1)/=nSubset) then

  ! deallocate space
  if(allocated(ixSubset))then
   deallocate(ixSubset,stat=err)
   if(err/=0)then; message=trim(message)//'unable to deallocate space for variable'; err=20; return; endif
  endif

  ! allocate space
  allocate(ixSubset(nSubset,nGRU),stat=err)
  if(err/=0)then; message=trim(message)//'unable to allocate space for variable'; err=20; return; endif

 endif  ! allocating space

 !$cuf kernel do(1) <<<*,*>>>
 do iGRU=1,nGRU
  subsetLayer(iGRU) = 1
  do iLayer=1,size(mask,1)
    if (mask(iLayer,iGRU)) then
      ixSubset(subsetLayer(iGRU),iGRU) = iLayer
      subsetLayer(iGRU) = subsetLayer(iGRU) + 1
    endif
  end do
end do
 ! define indices for variable types in specific sub-domains
!  if(nSubset>0) ixSubset = pack(ixMaster, mask)

 end subroutine indxSubset_d
 subroutine indxSubset_d2(ixSubset,ixMaster,mask,nSubset,nGRU,err,message)
  implicit none
  ! input-output: subset of indices for allocation/population
  integer(i4b),intent(inout),allocatable,device :: ixSubset(:,:) ! subset of indices
  ! input
  integer(i4b),intent(in),device :: ixMaster(:,:)
  logical(lgt),intent(in),device                :: mask(:,:)     ! desired indices
  integer(i4b) :: nGRU
  ! error control
  integer(i4b),intent(out)               :: err         ! error code
  character(*),intent(out)               :: message     ! error message
  ! -----------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  integer(i4b)                           :: nSubset     ! length of the subset
  integer(i4b),device :: subsetLayer(nGRU)
  integer(i4b) :: iGRU, iLayer
  ! -----------------------------------------------------------------------------------------------------------------------------------
  ! initialize errors
  err=0; message="indxSubset/"
 
  ! check size match
 !  if(size(ixMaster)/=size(mask))then
 !   message=trim(message)//'size mismatch'
 !   err=20; return
 !  endif
 
  ! check if we need to reallocate space
  if(size(ixSubset,1)/=nSubset) then
 
   ! deallocate space
   if(allocated(ixSubset))then
    deallocate(ixSubset,stat=err)
    if(err/=0)then; message=trim(message)//'unable to deallocate space for variable'; err=20; return; endif
   endif
 
   ! allocate space
   allocate(ixSubset(nSubset,nGRU),stat=err)
   if(err/=0)then; message=trim(message)//'unable to allocate space for variable'; err=20; return; endif
 
  endif  ! allocating space
 
  !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
   subsetLayer(iGRU) = 1
   do iLayer=1,size(mask,1)
     if (mask(iLayer,iGRU)) then
       ixSubset(subsetLayer(iGRU),iGRU) = ixMaster(iLayer,iGRU)
       subsetLayer(iGRU) = subsetLayer(iGRU) + 1
     endif
   end do
 end do
  ! define indices for variable types in specific sub-domains
 !  if(nSubset>0) ixSubset = pack(ixMaster, mask)
 
  end subroutine indxSubset_d2
 

 ! **********************************************************************************************************
 ! private subroutine resizeIndx: re-size specific index vectors
 ! **********************************************************************************************************
 subroutine resizeIndx(ixDesire,indx_data,nVec,err,message)
 ! input
 integer(i4b)     ,intent(in)    :: ixDesire(:)            ! variables needing to be re-sized
 type(var_ilength),intent(inout) :: indx_data              ! indices defining model states and layers
 integer(i4b)     ,intent(in)    :: nVec                   ! desired vector length
 ! output
 integer(i4b)     ,intent(out)   :: err                    ! error code
 character(*)     ,intent(out)   :: message                ! error message
 ! local variables
 integer(i4b)                    :: jVar,iVar              ! vatiable index
 ! initialize error control
 err=0; message='resizeIndx/'

 ! loop through variables
 do jVar=1,size(ixDesire)

  ! define index in index array
  iVar = ixDesire(jVar)

  ! check iVar is within range
  if(iVar<1 .or. iVar>size(indx_data%var))then
   message=trim(message)//'desired variable is out of range'
   err=20; return
  endif

  ! check if we need to reallocate space
  if(size(indx_data%var(iVar)%dat) == nVec) cycle

  ! deallocate space
  deallocate(indx_data%var(iVar)%dat,stat=err)
  if(err/=0)then
   message=trim(message)//'unable to deallocate space for variable '//trim(indx_meta(ivar)%varname)
   err=20; return
  endif

  ! allocate space
  allocate(indx_data%var(iVar)%dat(nVec),stat=err)
  if(err/=0)then
   message=trim(message)//'unable to allocate space for variable '//trim(indx_meta(ivar)%varname)
   err=20; return
  endif

  ! set to missing
  indx_data%var(iVar)%dat = integerMissing

 end do  ! looping through variables

 end subroutine resizeIndx

end module indexState_module
