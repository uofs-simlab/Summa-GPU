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

module layerMerge_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number

! access named variables for snow and soil
USE globalData,only:iname_snow        ! named variables for snow
USE globalData,only:iname_soil        ! named variables for soil

! access metadata
USE globalData,only:prog_meta,diag_meta,flux_meta,indx_meta   ! metadata

! physical constants
USE multiconst,only:&
                    iden_ice,       & ! intrinsic density of ice             (kg m-3)
                    iden_water        ! intrinsic density of liquid water    (kg m-3)

! access the derived types to define the data structures
USE data_types,only:&
                    var_d,            & ! data vector (rkind)
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength,      & ! data vector with variable length dimension (rkind)
                    model_options       ! defines the model decisions

! access named variables defining elements in the data structures
USE var_lookup,only:iLookPARAM,iLookPROG,iLookINDEX,iLookDIAG,iLookFLUX  ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS                   ! named variables for elements of the decision structure

! look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers, & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index

! provide access to external modules
USE var_derive_module,only:calcHeight,calcHeight_d ! module to calculate height at layer interfaces and layer mid-point

! privacy
implicit none
private
public::layerMerge

contains


 ! *****************************************************************************************************************
 ! public subroutine layerMerge: merge layers if the thickness is less than zmin
 ! *****************************************************************************************************************
 subroutine layerMerge(&
                       ! input/output: model data structures
                       tooMuchMelt,                 & ! intent(in):    flag to force merge of snow layers
                       nGRU, &
                       decisions,             & ! intent(in):    model decisions
                       mpar_data,                   & ! intent(in):    model parameters
                       indx_data,                   & ! intent(inout): type of each layer
                       prog_data,                   & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,                   & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,                   & ! intent(inout): model fluxes for a local HRU
                       ! output
                       mergedLayers,                & ! intent(out): flag to denote that layers were merged
                       err,message)                   ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 use device_data_types
!  use cudafor
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input/output: model data structures
 logical(lgt),intent(in)         :: tooMuchMelt         ! flag to denote that ice is insufficient to support melt
 type(decisions_device),intent(in)  :: decisions  ! model decisions
 type(mpar_data_device),intent(in)    :: mpar_data           ! model parameters
 type(indx_data_device),intent(inout) :: indx_data           ! type of each layer
 type(prog_data_device),intent(inout) :: prog_data           ! model prognostic variables for a local HRU
 type(diag_data_device),intent(inout) :: diag_data           ! model diagnostic variables for a local HRU
 type(flux_data_device),intent(inout) :: flux_data           ! model flux variables
 ! output
 logical(lgt),intent(out)        :: mergedLayers        ! flag to denote that layers were merged
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! define local variables
 character(LEN=256)              :: cmessage            ! error message of downwind routine
 real(rkind),dimension(5),device           :: zminLayer           ! minimum layer depth in each layer (m)
 logical(lgt)                    :: removeLayer         ! flag to indicate need to remove a layer
 integer(i4b)                    :: nCheck              ! number of layers to check for combination
 integer(i4b)                    :: iSnow               ! index of snow layers (looping)
 integer(i4b)                    :: jSnow               ! index of snow layer identified for combination with iSnow
 integer(i4b),device                    :: kSnow(nGRU)               ! index of the upper layer of the two layers identified for combination
 integer(i4b)                    :: nSnow               ! number of snow layers
 integer(i4b)                    :: nSoil               ! number of soil layers
!  integer(i4b)                    :: nLayers             ! total number of layers
 integer(i4b),device :: indexLayer(nGRU)
 integer(i4b) :: nGRU
 logical(lgt),device :: removeLastLayer(nGRU)
  logical(lgt),device :: combineLayers(nGRU)
  logical(lgt),device :: removeLayer_d(nGRU)
  integer(i4b) :: ix_snowLayers

 integer(i4b) :: iGRU, iLayer
 ! initialize error control
 err=0; message="layerMerge/"
 ix_snowLayers = decisions%snowLayers

 ! --------------------------------------------------------------------------------------------------------
 ! associate variables to the data structures
 associate(&

 ! model decisions
!  ix_snowLayers    => decisions%snowLayers, & ! decision for snow combination

 ! model parameters (control the depth of snow layers)
 zmin             => mpar_data%zmin,                & ! minimum layer depth (m)
 zminLayer1       => mpar_data%zminLayer1,          & ! minimum layer depth for the 1st (top) layer (m)
 zminLayer2       => mpar_data%zminLayer2,          & ! minimum layer depth for the 2nd layer (m)
 zminLayer3       => mpar_data%zminLayer3,          & ! minimum layer depth for the 3rd layer (m)
 zminLayer4       => mpar_data%zminLayer4,          & ! minimum layer depth for the 4th layer (m)
 zminLayer5       => mpar_data%zminLayer5,          & ! minimum layer depth for the 5th (bottom) layer (m)

 ! diagnostic scalar variables
 scalarSnowDepth  => prog_data%scalarSnowDepth,      & ! total snow depth (m)
 scalarSWE        => prog_data%scalarSWE             & ! SWE (kg m-2)

 ) ! end associate statement
 ! --------------------------------------------------------------------------------------------------------

 ! identify algorithmic control parameters to sub-divide and combine snow layers
 zminLayer = (/zminLayer1, zminLayer2, zminLayer3, zminLayer4, zminLayer5/)

 ! intialize the modified layers flag
 mergedLayers=.false.

 ! initialize the number of snow layers
 nSnow   = maxval(indx_data%nSnow)
 nSoil   = indx_data%nSoil

 kSnow=0 ! initialize first layer to test (top layer)
 do ! attempt to remove multiple layers in a single time step (continuous do loop with exit clause)

  ! special case of >5 layers: add an offset to use maximum threshold from layer above
  if(ix_snowLayers == rulesDependLayerIndex .and. nSnow > 5)then
   nCheck=5
  else
   nCheck=nSnow
  end if

  removeLastLayer = .false.
  removeLayer = .false.
  removeLayer_d = .false.
  combineLayers = .false.
     ! NOTE: do this here, since the layer variables are re-defined
   associate(&
   mLayerDepth      => prog_data%mLayerDepth         , &    ! depth of each layer (m)
   nSnow => indx_data%nSnow, &
   mLayerVolFracIce => prog_data%mLayerVolFracIce    , &    ! volumetric fraction of ice in each layer  (-)
   mLayerVolFracLiq => prog_data%mLayerVolFracLiq      &    ! volumetric fraction of liquid water in each layer (-)
   ) ! (associating local variables with the information in the data structures)

   !$cuf kernel do(1) <<<*,*>>> reduce(.or.:removeLayer)
   do iGRU=1,nGRU
  ! loop through snow layers
  do iSnow=1,min(nCheck,nSnow(iGRU))

   ! associate local variables with the information in the data structures

   ! check if the layer depth is less than the depth threshold
   select case(ix_snowLayers)
    case(sameRulesAllLayers);    removeLayer_d(iGRU) = (mLayerDepth(iSnow,iGRU) < zmin)
    case(rulesDependLayerIndex); removeLayer_d(iGRU) = (mLayerDepth(iSnow,iGRU) < zminLayer(iSnow))
    ! case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
   end select ! (option to combine/sub-divide snow layers)

   ! check if we have too much melt
   ! NOTE: assume that this is the top snow layer; need more trickery to relax this assumption
   if(tooMuchMelt .and. iSnow==1) removeLayer_d(iGRU)=.true.

   ! check if need to remove a layer
   if(removeLayer_d(iGRU))then
    removeLayer = .true.

    ! flag that we modified a layer
    mergedLayers=.true.

    ! ***** handle special case of a single layer
    if(nSnow(iGRU)==1)then
      removeLastLayer(iGRU) = .true.
     ! set the variables defining "snow without a layer"
     ! NOTE: ignoring cold content!!! Need to fix later...
     scalarSnowDepth(iGRU) = mLayerDepth(1,iGRU)
     scalarSWE(iGRU)       = (mLayerVolFracIce(1,iGRU)*iden_ice + mLayerVolFracLiq(1,iGRU)*iden_water)*mLayerDepth(1,iGRU)
     ! remove the top layer from all model variable vectors
     ! NOTE: nSnow-1 = 0, so routine removes layer #1
     kSnow(iGRU) = nSnow(iGRU)-1
     indexLayer(iGRU) = kSnow(iGRU)
    end if  ! (special case of 1 layer --> snow without a layer)

    ! ***** identify the layer to combine
    if (nSnow(iGRU) .ne. 1) then
      combineLayers(iGRU) = .true.
    if(iSnow==1)then
     jSnow = iSnow+1  ! upper-most layer, combine with its lower neighbor
    elseif(iSnow==nSnow(iGRU))then
     jSnow = nSnow(iGRU)-1  ! lower-most layer, combine with its upper neighbor
    else
     if(mLayerDepth(iSnow-1,iGRU)<mLayerDepth(iSnow+1,iGRU))then; jSnow = iSnow-1; else; jSnow = iSnow+1; end if
    end if

    ! ***** combine layers
    ! identify the layer closest to the surface
    kSnow(iGRU)=min(iSnow,jSnow)
    indexLayer(iGRU) = kSnow(iGRU)
    ! combine layer with identified neighbor
    ! call layer_combine(mpar_data,prog_data,diag_data,flux_data,indx_data,kSnow,err,cmessage)
    ! if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if

    ! ! update the number of snow layers
    ! nSnow   = indx_data%var(iLookINDEX%nSnow)%dat(1)
    ! nSoil   = indx_data%var(iLookINDEX%nSoil)%dat(1)
    ! nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)

    ! exit the loop to try again
    exit
  end if

   end if  ! (if layer is below the mass threshold)

   kSnow(iGRU)=iSnow ! ksnow is used for completion test, so include here


  end do ! (looping through snow layers)
end do
     ! end association of local variables with the information in the data structures
   end associate



  if(removeLayer)then

    ! flag that we modified a layer
    mergedLayers=.true.

    ! ***** handle special case of a single layer
    ! if(removeLastLayer)then
     ! remove the top layer from all model variable vectors
     ! NOTE: nSnow-1 = 0, so routine removes layer #1
     call rmLyAllVars_prog(prog_data,prog_meta,indexLayer,removeLastLayer,nGRU,nSoil,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
     call rmLyAllVars_diag(diag_data,diag_meta,indexLayer,removeLastLayer,nGRU,nSoil,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
     call rmLyAllVars_flux(flux_data,flux_meta,indexLayer,removeLastLayer,nGRU,nSoil,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
     call rmLyAllVars_index(indx_data,indx_meta,indexLayer,removeLastLayer,nGRU,nSoil,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

     if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if

     ! update the total number of layers
     associate(nSnow => indx_data%nSnow, nLayers=>indx_data%nLayers_d)
      !$cuf kernel do(1) <<<*,*>>>
     do iGRU=1,nGRU
      if (removeLastLayer(iGRU)) nSnow(iGRU)   = nSnow(iGRU) - 1
      ! if (removeLastLayer(iGRU)) nLayers(iGRU) = nLayers(iGRU) - 1
      nLayers(iGRU) = nSoil+nSnow(iGRU)
     end do
          end associate     

     ! save the number of layers
     nSnow = maxval(indx_data%nSnow)
     nSoil = indx_data%nSoil
     ! update coordinate variables
     call calcHeight_d(&
     nGRU, &
                     ! input/output: data structures
                     indx_data,   & ! intent(in): layer type
                     prog_data,   & ! intent(inout): model variables for a local HRU
                     ! output: error control
                     err,cmessage)
     if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

    !  err = cudaDeviceSynchronize()
    !  err = 20
     ! exit the do loop (no more snow layers to remove)
    !  return
    ! else  ! (special case of 1 layer --> snow without a layer)
    ! end if
    ! if (combineLayers) then

    ! ! ***** identify the layer to combine
    ! if(iSnow==1)then
    !  jSnow = iSnow+1  ! upper-most layer, combine with its lower neighbor
    ! elseif(iSnow==nSnow)then
    !  jSnow = nSnow-1  ! lower-most layer, combine with its upper neighbor
    ! else
    !  if(mLayerDepth(iSnow-1)<mLayerDepth(iSnow+1))then; jSnow = iSnow-1; else; jSnow = iSnow+1; end if
    ! endif

    ! ! ***** combine layers
    ! ! identify the layer closest to the surface
    ! kSnow=min(iSnow,jSnow)
    ! combine layer with identified neighbor

    call layer_combine(mpar_data,prog_data,diag_data,flux_data,indx_data,nGRU,indexLayer,combineLayers,err,cmessage,decisions%h_lookup,decisions%t_lookup)
    if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if
    !  err = cudaDeviceSynchronize()

    ! update the number of snow layers
    nSnow   = maxval(indx_data%nSnow)
    nSoil   = indx_data%nSoil


  !  end if  ! (if layer is below the mass threshold)
  end if

  ! exit if finished
  if(.not. removeLayer)exit

 end do ! continuous do

 ! handle special case of > 5 layers in the CLM option
 if(maxval(indx_data%nSnow) > 5 .and. ix_snowLayers == rulesDependLayerIndex)then
  ! flag that layers were merged
  mergedLayers=.true.
  ! initial check to ensure everything is wonderful in the universe
  if(maxval(indx_data%nSnow) /= 6)then; err=5; message=trim(message)//'special case of > 5 layers: expect only six layers'; return; end if
      indexLayer = 5
      associate(nSnow => indx_data%nSnow)
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        combineLayers(iGRU)  = nSnow(iGRU) > 5
      end do
      end associate

  ! combine 5th layer with layer below
    call layer_combine(mpar_data,prog_data,diag_data,flux_data,indx_data,nGRU,indexLayer,combineLayers,err,cmessage,decisions%h_lookup,decisions%t_lookup)
  ! update the number of snow layers
  nSnow   = maxval(indx_data%nSnow)
  nSoil   = indx_data%nSoil
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if
  ! another check
  if(nSnow /= 5)then; err=5; message=trim(message)//'special case of > 5 layers: expect to reduced layers to exactly 5'; return; end if

 end if

 ! check that there are no more than 5 layers in the CLM option
 if(ix_snowLayers == rulesDependLayerIndex)then
  if(maxval(indx_data%nSnow) > 5)then
   message=trim(message)//'expect no more than 5 layers when combination/sub-division rules depend on the layer index (CLM option)'
   err=20; return
  end if
 end if

 ! end association to variables in the data structure
 end associate

 end subroutine layerMerge


 ! ***********************************************************************************************************
 ! private subroutine layer_combine: combine snow layers and re-compute model state variables
 ! ***********************************************************************************************************
 ! combines layer iSnow with iSnow+1
 ! ***********************************************************************************************************
 subroutine layer_combine(mpar_data,prog_data,diag_data,flux_data,indx_data,nGRU,iSnow,combineLayer,err,message,h_lookup_d,t_lookup_d)
 ! provide access to variables in the data structures
 USE var_lookup,only:iLookPARAM,iLookPROG,iLookINDEX              ! named variables for structure elements
 USE globalData,only:prog_meta,diag_meta,flux_meta,indx_meta      ! metadata
 USE data_types,only:var_ilength,var_dlength                      ! data vectors with variable length dimension
 USE data_types,only:var_d                                        ! data structures with fixed dimension
 ! provide access to external modules
 USE snow_utils_module,only:fracliquid_d                            ! compute fraction of liquid water
 USE enthalpyTemp_module,only:enthalpy2T_snwWat,T2enthalpy_snwWat_d ! convert temperature to liq+ice enthalpy for a snow layer
 use device_data_types
 implicit none
 ! ------------------------------------------------------------------------------------------------------------
 ! input/output: data structures
 type(mpar_data_device),intent(in)    :: mpar_data ! model parameters
 type(prog_data_device),intent(inout) :: prog_data ! model prognostic variables for a local HRU
 type(diag_data_device),intent(inout) :: diag_data ! model diagnostic variables for a local HRU
 type(flux_data_device),intent(inout) :: flux_data ! model flux variables
 type(indx_data_device),intent(inout) :: indx_data ! type of model layer
 ! input: snow layer indices
 integer(i4b),intent(in),device         :: iSnow(:)     ! index of top layer to combine
 logical(lgt),device :: combineLayer(:)
 ! output: error control
 integer(i4b),intent(out)        :: err       ! error code
 character(*),intent(out)        :: message   ! error message
 ! ------------------------------------------------------------------------------------------------------------
 ! local variables
 character(len=256)              :: cmessage                 ! error message for downwind routine
 real(rkind),device                     :: massIce(2,nGRU)               ! mass of ice in the two layers identified for combination (kg m-2)
 real(rkind),device                     :: massLiq(2,nGRU)               ! mass of liquid water in the two layers identified for combination (kg m-2)
 real(rkind),device                     :: bulkDenWat(2,nGRU)            ! bulk density if total water (liquid water plus ice) in the two layers identified for combination (kg m-3)
 real(rkind),device                     :: cBulkDenWat(nGRU)              ! combined bulk density of total water (liquid water plus ice) in the two layers identified for combination (kg m-3)
 real(rkind),device                     :: cTemp(nGRU)                    ! combined layer temperature
 real(rkind),device                     :: cDepth(nGRU)                   ! combined layer depth
 real(rkind),device                     :: cVolFracIce(nGRU)              ! combined layer volumetric fraction of ice
 real(rkind),device                     :: cVolFracLiq(nGRU)              ! combined layer volumetric fraction of liquid water
 real(rkind),device                     :: l1Enthalpy(nGRU),l2Enthalpy(nGRU)    ! enthalpy in the two layers identified for combination (J m-3)
 real(rkind),device                     :: cEnthalpy(nGRU)                ! combined layer enthalpy (J m-3)
 real(rkind),device                     :: fLiq(nGRU)                     ! fraction of liquid water at the combined temperature cTemp
 real(rkind),parameter           :: eTol=1.e-1_rkind         ! tolerance for the enthalpy-->temperature conversion (J m-3)
 integer(i4b)                    :: nSnow                    ! number of snow layers
 integer(i4b)                    :: nSoil                    ! number of soil layers
!  integer(i4b)                    :: nLayers                  ! total number of layers
 integer(i4b) :: nGRU
!  type(indx_data_device) :: indx_data_d
!  type(flux_data_device) :: flux_data_d
!  type(prog_data_device) :: prog_data_d
!  type(diag_data_device) :: diag_data_d
 integer(i4b) :: iGRU, iLayer
 real(rkind),device :: h_lookup_d(:), t_lookup_d(:)

!  call allocate_device_param_data(mpar_data_d,mpar_data)
 ! initialize error control
 err=0; message="layer_combine/"

 ! associate local variables with information in the data structures
 associate(&
 ! model parameters
 snowfrz_scale    => mpar_data%snowfrz_scale, & ! scaling parameter for the freezing curve for snow (K-1)
 ! model state variables
 mLayerTemp       => prog_data%mLayerTemp       , & ! temperature of each layer (K)
 mLayerDepth      => prog_data%mLayerDepth      , & ! depth of each layer (m)
 mLayerVolFracIce => prog_data%mLayerVolFracIce , & ! volumetric fraction of ice in each layer  (-)
 mLayerVolFracLiq => prog_data%mLayerVolFracLiq   & ! volumetric fraction of liquid water in each layer (-)
 ) ! (association of local variables with information in the data structures)

 ! initialize the number of snow layers
 nSnow   = indx_data%nSnow(1)
 nSoil   = indx_data%nSoil
  ! call allocate_device_prog_data(prog_data_d,prog_data,nGRU,nLayers,nSoil)
!  call allocate_device_diag_data(diag_data_d,diag_data,nSnow,nGRU,nLayers,nSoil)
!  call allocate_device_flux_data(flux_data_d,flux_data,nSnow,nSoil,nGRU,nLayers)
!  call allocate_device_indx_data(indx_data_d,indx_data,nSnow,nSoil,nLayers,nGRU)


  !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    if (combineLayer(iGRU)) then

 ! compute combined depth
 cDepth(iGRU)       = mLayerDepth(isnow(iGRU),iGRU) + mLayerDepth(isnow(iGRU)+1,iGRU)

 ! compute mass of each layer (kg m-2)
 massIce(1,iGRU) = iden_ice*mLayerVolFracIce(iSnow(iGRU),iGRU)*mLayerDepth(iSnow(iGRU),iGRU)
 massIce(2,iGRU) = iden_ice*mLayerVolFracIce(iSnow(iGRU)+1,iGRU)*mLayerDepth(iSnow(iGRU)+1,iGRU)
 massLiq(1,iGRU) = iden_water*mLayerVolFracLiq(iSnow(iGRU),iGRU)*mLayerDepth(iSnow(iGRU),iGRU)
 massLiq(2,iGRU) = iden_water*mLayerVolFracLiq(iSnow(iGRU)+1,iGRU)*mLayerDepth(iSnow(iGRU)+1,iGRU)

 ! compute bulk density of water (kg m-3)
 bulkDenWat(1,iGRU) = (massIce(1,iGRU) + massLiq(1,iGRU))/mLayerDepth(iSnow(iGRU),iGRU)
 bulkDenWat(2,iGRU) = (massIce(2,iGRU) + massLiq(2,iGRU))/mLayerDepth(iSnow(iGRU)+1,iGRU)

 cBulkDenWat(iGRU)     = (mLayerDepth(isnow(iGRU),iGRU)*bulkDenWat(1,iGRU) + mLayerDepth(isnow(iGRU)+1,iGRU)*bulkDenWat(2,iGRU))/cDepth(iGRU)

 ! compute enthalpy for each layer (J m-3)
 l1Enthalpy(iGRU) = T2enthalpy_snwWat_d(mLayerTemp(iSnow(iGRU),iGRU),  BulkDenWat(1,iGRU),snowfrz_scale)
 l2Enthalpy(iGRU) = T2enthalpy_snwWat_d(mLayerTemp(iSnow(iGRU)+1,iGRU),BulkDenWat(2,iGRU),snowfrz_scale)
 ! compute combined enthalpy (J m-3)
 cEnthalpy(iGRU) = (mLayerDepth(isnow(iGRU),iGRU)*l1Enthalpy(iGRU) + mLayerDepth(isnow(iGRU)+1,iGRU)*l2Enthalpy(iGRU))/cDepth(iGRU)
    end if
  end do
 ! convert enthalpy (J m-3) to temperature (K)
  do iGRU=1,nGRU
 call enthalpy2T_snwWat<<<1,1>>>(combineLayer(iGRU), cEnthalpy(iGRU),cBulkDenWat(iGRU),snowfrz_scale,cTemp(iGRU),h_lookup_d,t_lookup_d)
  end do
!  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; end if
  !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    if (combineLayer(iGRU)) then

 ! test enthalpy conversion
!  if(abs(T2enthalpy_snwWat(cTemp,cBulkDenWat,snowfrz_scale)/cBulkDenWat - cEnthalpy/cBulkDenWat) > eTol)then
!   write(*,'(a,1x,f12.5,1x,2(e20.10,1x))') 'enthalpy test', cBulkDenWat, T2enthalpy_snwWat(cTemp,cBulkDenWat,snowfrz_scale)/cBulkDenWat, cEnthalpy/cBulkDenWat
!   message=trim(message)//'problem with enthalpy-->temperature conversion'
!   err=20; return
!  end if

 ! check temperature is within the two temperatures
 ! NOTE: use tolerance, for cases of merging a layer that has just been split
!  if(cTemp > max(mLayerTemp(iSnow),mLayerTemp(iSnow+1))+eTol)then; err=20; message=trim(message)//'merged temperature > max(temp1,temp2)'; return; end if
!  if(cTemp < min(mLayerTemp(iSnow),mLayerTemp(iSnow+1))-eTol)then; err=20; message=trim(message)//'merged temperature < min(temp1,temp2)'; return; end if

 ! compute volumetric fraction of liquid water
 fLiq(iGRU) = fracLiquid_d(cTemp(iGRU),snowfrz_scale)
 ! compute volumetric fraction of ice and liquid water
 cVolFracLiq(iGRU) =          fLiq(iGRU) *cBulkDenWat(iGRU)/iden_water
 cVolFracIce(iGRU) = (1._rkind - fLiq(iGRU))*cBulkDenWat(iGRU)/iden_ice
end if
  end do

 ! end association of local variables with information in the data structures
 end associate

 ! remove a model layer from all model variable vectors
     call rmLyAllVars_prog(prog_data,prog_meta,iSnow,combineLayer,nGRU,nSoil,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
     call rmLyAllVars_diag(diag_data,diag_meta,iSnow,combineLayer,nGRU,nSoil,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
     call rmLyAllVars_flux(flux_data,flux_meta,iSnow,combineLayer,nGRU,nSoil,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
     call rmLyAllVars_index(indx_data,indx_meta,iSnow,combineLayer,nGRU,nSoil,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if


 ! define the combined layer as snow
 associate(layerType => indx_data%layerType)
  !$cuf kernel do(1) <<<*,*>>>
 do iGRU=1,nGRU
  if (combineLayer(iGRU)) layerType(iSnow(iGRU),iGRU) = iname_snow
 end do
 end associate

!  indx_data%var(iLookINDEX%layerType)%dat(iSnow) = iname_snow



 ! save the number of layers in the data structures
     associate(layerType => indx_data%layerType, nSnow => indx_data%nSnow,nLayers => indx_data%nLayers_d)
       !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    nSnow(iGRU) = 0
    do iLayer=1,size(layerType,1)
      if (layerType(iLayer,iGRU) == iname_snow) nSnow(iGRU) = nSnow(iGRU) + 1
    end do
    nLayers(iGRU) = nSoil+nSnow(iGRU)
  end do
  
end associate
!  indx_data%var(iLookINDEX%nSnow)%dat(1)   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_snow)
!  indx_data%var(iLookINDEX%nSoil)%dat(1)   = count(indx_data%var(iLookINDEX%layerType)%dat==iname_soil)


 ! update the number of snow layers
 nSnow   = indx_data%nSnow(1)
 nSoil   = indx_data%nSoil

 ! ***** put state variables for the combined layer in the appropriate place
     associate(mLayerTemp => prog_data%mlayerTemp, &
      mLayerDepth => prog_data%mLayerDepth, &
      mLayerVolFracIce => prog_data%mLayerVolFracIce, &
      mLayerVolFracLiq => prog_data%mLayerVolFracLiq)
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        if (combineLayer(iGRU)) then
 mLayerTemp(iSnow(iGRU),iGRU)       = cTemp(iGRU)
 mLayerDepth(iSnow(iGRU),iGRU)      = cDepth(iGRU)
 mLayerVolFracIce(iSnow(iGRU),iGRU) = cVolFracIce(iGRU)
 mLayerVolFracLiq(iSnow(iGRU),iGRU) = cVolFracLiq(iGRU)
        end if
      end do
      end associate

 ! ***** adjust coordinate variables
 call calcHeight_d(&
 nGRU, &
                 ! input/output: data structures
                 indx_data,   & ! intent(in): layer type
                 prog_data,   & ! intent(inout): model variables for a local HRU
                 ! output: error control
                 err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if
    !  call finalize_device_prog_data(prog_data_d,prog_data,nLayers,nSoil)
    !  call deallocate_device_prog_data(prog_data_d)
    !  call finalize_device_diag_data(diag_data_d,diag_data,nSnow,nSoil,nLayers)
    !  call deallocate_device_diag_data(diag_data_d,nSnow)
    !  call finalize_device_flux_data(flux_data_d,flux_data,nSnow,nSoil)
    !  call deallocate_device_flux_data(flux_data_d,nSnow,nSoil)
    !  call finalize_device_indx_data(indx_data_d,indx_data)
    !  call deallocate_device_indx_data(indx_data_d)


 end subroutine layer_combine


 ! ***********************************************************************************************************
 ! private subroutine rmLyAllVars: reduce the length of the vectors in data structures
 ! ***********************************************************************************************************
 ! removes layer "iSnow+1" and sets layer "iSnow" to a missing value
 ! (layer "iSnow" will be filled with a combined layer later)
 ! ***********************************************************************************************************
 
 subroutine rmLyAllVars_prog(dataStruct,metaStruct,iSnow,removeLayer,nGRU,nSoil,err,message)
 USE var_lookup,only:iLookVarType                 ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName       ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc            ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength      ! data vectors with variable length dimension
 USE data_types,only:var_info                     ! metadata structure
 use device_data_types
  use globalData,only:maxSnowLayers

 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 type(prog_data_device),intent(inout)          :: dataStruct     ! data structure
 type(var_info),intent(in)       :: metaStruct(:)  ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in),device         :: iSnow(:)          ! new layer
 logical(lgt),intent(in),device :: removeLayer(:)
 integer(i4b),intent(in) :: nGRU
 integer(i4b),intent(in)         :: nSoil  ! number of snow layers, total number of layers
 ! output: error control
 integer(i4b),intent(out)        :: err            ! error code
 character(*),intent(out)        :: message        ! error message
 ! locals
 integer(i4b)                    :: ivar           ! variable index
 integer(i4b)                    :: ix_lower       ! lower bound of the vector
 integer(i4b)                    :: ix_upper       ! upper bound of the vector
 real(rkind),allocatable,device            :: tempVec(:,:)  ! temporary vector (double precision)
 character(LEN=256)              :: cmessage       ! error message of downwind routine
 ! initialize error control
 err=0; message="rmLyAllVars_prog/"

 ! ***** loop through model variables and remove one layer
 do ivar=1,size(metaStruct)

  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=maxSnowLayers
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nSoil+maxSnowLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=maxSnowLayers
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nSoil+maxSnowLayers
   case default; cycle  ! no need to remove soil layers or scalar variables
  end select

  ! allocate the temporary vector
  allocate(tempVec(ix_lower:ix_upper-1,nGRU), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate temporary vector'; return; end if
  select case(iVar)
    case(iLookPROG%spectralSnowAlbedoDiffuse); call populate_data(dataStruct%spectralSnowAlbedoDiffuse,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
    case(iLookPROG%mLayerTemp); call populate_data(dataStruct%mLayerTemp,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
    case(iLookPROG%mLayerVolFracIce); call populate_data(dataStruct%mLayerVolFracIce,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
    case(iLookPROG%mLayerVolFracLiq); call populate_data(dataStruct%mLayerVolFracLiq,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
    case(iLookPROG%mLayerVolFracWat); call populate_data(dataStruct%mLayerVolFracWat,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
    case(iLookPROG%mLayerMatricHead); call populate_data(dataStruct%mLayerMatricHead,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
    case(iLookPROG%mLayerEnthalpy); call populate_data(dataStruct%mLayerEnthalpy,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
    case(iLookPROG%mLayerDepth); call populate_data(dataStruct%mLayerDepth,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
    case(iLookPROG%mLayerHeight); call populate_data(dataStruct%mLayerHeight,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
    case(iLookPROG%iLayerHeight); call populate_data(dataStruct%iLayerHeight,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  end select

  deallocate(tempVec,stat=err)
  if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

  
 end do  ! looping through variables

 end subroutine rmLyAllVars_prog


 subroutine rmLyAllVars_diag(dataStruct,metaStruct,iSnow,removeLayer,nGRU,nSoil,err,message)
 USE var_lookup,only:iLookVarType                 ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName       ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc            ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength      ! data vectors with variable length dimension
 USE data_types,only:var_info                     ! metadata structure
 use device_data_types
 use globalData,only:maxSnowLayers
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 type(diag_data_device),intent(inout)          :: dataStruct     ! data structure
 type(var_info),intent(in)       :: metaStruct(:)  ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in),device         :: iSnow(:)          ! new layer
 logical(lgt),intent(in),device :: removeLayer(:)
 integer(i4b),intent(in) :: nGRU
 integer(i4b),intent(in)         :: nSoil  ! number of snow layers, total number of layers
 ! output: error control
 integer(i4b),intent(out)        :: err            ! error code
 character(*),intent(out)        :: message        ! error message
 ! locals
 integer(i4b)                    :: ivar           ! variable index
 integer(i4b)                    :: ix_lower       ! lower bound of the vector
 integer(i4b)                    :: ix_upper       ! upper bound of the vector
 real(rkind),allocatable,device            :: tempVec(:,:)  ! temporary vector (double precision)
 character(LEN=256)              :: cmessage       ! error message of downwind routine
 ! initialize error control
 err=0; message="rmLyAllVars_diag/"

 ! ***** loop through model variables and remove one layer
 do ivar=1,size(metaStruct)

  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=maxSnowLayers
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nSoil+maxSnowLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=maxSnowLayers
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nSoil+maxSnowLayers
   case default; cycle  ! no need to remove soil layers or scalar variables
  end select

  ! allocate the temporary vector
  allocate(tempVec(ix_lower:ix_upper-1,nGRU), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate temporary vector'; return; end if
  select case(iVar)
  case(iLookDIAG%mLayerVolHtCapBulk); call populate_data(dataStruct%mLayerVolHtCapBulk_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerCm); call populate_data(dataStruct%mLayerCm_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerThermalC); call populate_data(dataStruct%mLayerThermalC_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%iLayerThermalC); call populate_data(dataStruct%iLayerThermalC_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerEnthTemp); call populate_data(dataStruct%mLayerEnthTemp,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%spectralAlbGndDirect); call populate_data(dataStruct%spectralAlbGndDirect,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%spectralAlbGndDiffuse); call populate_data(dataStruct%spectralAlbGndDiffuse,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerTranspireLim); call populate_data(dataStruct%mLayerTranspireLim_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerRootDensity); call populate_data(dataStruct%mLayerRootDensity_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%spectralSnowAlbedoDirect); call populate_data(dataStruct%spectralSnowAlbedoDirect,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerFracLiqSnow); call populate_data(dataStruct%mLayerFracLiqSnow_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerThetaResid); call populate_data(dataStruct%mLayerThetaResid_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerPoreSpace); call populate_data(dataStruct%mLayerPoreSpace_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerMeltFreeze); call populate_data(dataStruct%mLayerMeltFreeze,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerVolFracAir); call populate_data(dataStruct%mLayerVolFracAir_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerTcrit); call populate_data(dataStruct%mLayerTcrit,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerCompress); call populate_data(dataStruct%mLayerCompress_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%mLayerMatricHeadLiq); call populate_data(dataStruct%mLayerMatricHeadLiq,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%scalarVGn_m); call populate_data(dataStruct%scalarVGn_m_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%balanceLayerNrg); call populate_data(dataStruct%balanceLayerNrg,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookDIAG%balanceLayerMass); call populate_data(dataStruct%balanceLayerMass,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  end select

  deallocate(tempVec,stat=err)
  if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

  
 end do  ! looping through variables

 end subroutine rmLyAllVars_diag


 subroutine rmLyAllVars_flux(dataStruct,metaStruct,iSnow,removeLayer,nGRU,nSoil,err,message)
 USE var_lookup,only:iLookVarType                 ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName       ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc            ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength      ! data vectors with variable length dimension
 USE data_types,only:var_info                     ! metadata structure
 use device_data_types
  use globalData,only:maxSnowLayers
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 type(flux_data_device),intent(inout)          :: dataStruct     ! data structure
 type(var_info),intent(in)       :: metaStruct(:)  ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in),device         :: iSnow(:)          ! new layer
 logical(lgt),intent(in),device :: removeLayer(:)
 integer(i4b),intent(in) :: nGRU
 integer(i4b),intent(in)         :: nSoil  ! number of snow layers, total number of layers
 ! output: error control
 integer(i4b),intent(out)        :: err            ! error code
 character(*),intent(out)        :: message        ! error message
 ! locals
 integer(i4b)                    :: ivar           ! variable index
 integer(i4b)                    :: ix_lower       ! lower bound of the vector
 integer(i4b)                    :: ix_upper       ! upper bound of the vector
 real(rkind),allocatable,device            :: tempVec(:,:)  ! temporary vector (double precision)
 character(LEN=256)              :: cmessage       ! error message of downwind routine
 ! initialize error control
 err=0; message="rmLyAllVars_flux/"

 ! ***** loop through model variables and remove one layer
 do ivar=1,size(metaStruct)

  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=maxSnowLayers
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nSoil+maxSnowLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=maxSnowLayers
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nSoil+maxSnowLayers
   case default; cycle  ! no need to remove soil layers or scalar variables
  end select

  ! allocate the temporary vector
  allocate(tempVec(ix_lower:ix_upper-1,nGRU), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate temporary vector'; return; end if
  select case(iVar)
  case(iLookFLUX%spectralIncomingDirect); call populate_data(dataStruct%spectralIncomingDirect,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%spectralIncomingDiffuse); call populate_data(dataStruct%spectralIncomingDiffuse,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%spectralBelowCanopyDirect); call populate_data(dataStruct%spectralBelowCanopyDirect,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%spectralBelowCanopyDiffuse); call populate_data(dataStruct%spectralBelowCanopyDiffuse,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%mLayerTranspire); call populate_data(dataStruct%mLayerTranspire_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%iLayerConductiveFlux); call populate_data(dataStruct%iLayerConductiveFlux_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%iLayerAdvectiveFlux); call populate_data(dataStruct%iLayerAdvectiveFlux_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%iLayerNrgFlux); call populate_data(dataStruct%iLayerNrgFlux_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%mLayerNrgFlux); call populate_data(dataStruct%mLayerNrgFlux_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%iLayerLiqFluxSnow); call populate_data(dataStruct%iLayerLiqFluxSnow_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%mLayerLiqFluxSnow); call populate_data(dataStruct%mLayerLiqFluxSnow_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%mLayerSatHydCondMP); call populate_data(dataStruct%mLayerSatHydCondMP_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%mLayerSatHydCond); call populate_data(dataStruct%mLayerSatHydCond_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%iLayerSatHydCond); call populate_data(dataStruct%iLayerSatHydCond_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%mLayerHydCond); call populate_data(dataStruct%mLayerHydCond_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%iLayerLiqFluxSoil); call populate_data(dataStruct%iLayerLiqFluxSoil_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%mLayerLiqFluxSoil); call populate_data(dataStruct%mLayerLiqFluxSoil_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%mLayerBaseflow); call populate_data(dataStruct%mLayerBaseflow_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%mLayerColumnInflow); call populate_data(dataStruct%mLayerColumnInflow,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  case(iLookFLUX%mLayerColumnOutflow); call populate_data(dataStruct%mLayerColumnOutflow_m,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  end select

  deallocate(tempVec,stat=err)
  if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

  
 end do  ! looping through variables

 end subroutine rmLyAllVars_flux


  subroutine rmLyAllVars_index(dataStruct,metaStruct,iSnow,removeLayer,nGRU,nSoil,err,message)
 USE var_lookup,only:iLookVarType                 ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName       ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc            ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength      ! data vectors with variable length dimension
 USE data_types,only:var_info                     ! metadata structure
 use device_data_types
  use globalData,only:maxSnowLayers

 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 type(indx_data_device),intent(inout)          :: dataStruct     ! data structure
 type(var_info),intent(in)       :: metaStruct(:)  ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in),device         :: iSnow(:)          ! new layer
 logical(lgt),intent(in),device :: removeLayer(:)
 integer(i4b),intent(in) :: nGRU
 integer(i4b),intent(in)         :: nSoil  ! number of snow layers, total number of layers
 ! output: error control
 integer(i4b),intent(out)        :: err            ! error code
 character(*),intent(out)        :: message        ! error message
 ! locals
 integer(i4b)                    :: ivar           ! variable index
 integer(i4b)                    :: ix_lower       ! lower bound of the vector
 integer(i4b)                    :: ix_upper       ! upper bound of the vector
 integer(i4b),allocatable,device            :: tempVec(:,:)  ! temporary vector (integer)
 character(LEN=256)              :: cmessage       ! error message of downwind routine
 ! initialize error control
 err=0; message="rmLyAllVars_indx/"

 ! ***** loop through model variables and remove one layer
 do ivar=1,size(metaStruct)

  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=maxSnowLayers
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nSoil+maxSnowLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=maxSnowLayers
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nSoil+maxSnowLayers
   case default; cycle  ! no need to remove soil layers or scalar variables
  end select

  ! allocate the temporary vector
  allocate(tempVec(ix_lower:ix_upper-1,nGRU), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate temporary vector'; return; end if

  select case(iVar)
case(iLookINDEX%layerType); call populate_data_i(dataStruct%layerType,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixControlVolume); call populate_data_i(dataStruct%ixControlVolume,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixDomainType); call populate_data_i(dataStruct%ixDomainType,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixStateType); call populate_data_i(dataStruct%ixStateType,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixHydType); call populate_data_i(dataStruct%ixHydType,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
! case(iLookINDEX%ixDomainType_subset); call populate_data_i(dataStruct%ixDomainType_subset,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
! case(iLookINDEX%ixStateType_subset); call populate_data_i(dataStruct%ixStateType_subset,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
! case(iLookINDEX%ixMapFull2Subset); call populate_data_i(dataStruct%ixMapFull2Subset,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
! case(iLookINDEX%ixMapSubset2Full); call populate_data_i(dataStruct%ixMapSubset2Full,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixNrgOnly); call populate_data_i(dataStruct%ixNrgOnly,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
! case(iLookINDEX%ixHydOnly); call populate_data_i(dataStruct%ixHydOnly,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
! case(iLookINDEX%ixMatOnly); call populate_data_i(dataStruct%ixMatOnly,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
! case(iLookINDEX%ixMassOnly); call populate_data_i(dataStruct%ixMassOnly,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixSnowSoilNrg); call populate_data_i(dataStruct%ixSnowSoilNrg,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixSnowOnlyNrg); call populate_data_i(dataStruct%ixSnowOnlyNrg,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixSoilOnlyNrg); call populate_data_i(dataStruct%ixSoilOnlyNrg,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixSnowSoilHyd); call populate_data_i(dataStruct%ixSnowSoilHyd,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixSnowOnlyHyd); call populate_data_i(dataStruct%ixSnowOnlyHyd,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixSoilOnlyHyd); call populate_data_i(dataStruct%ixSoilOnlyHyd,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixNrgLayer); call populate_data_i(dataStruct%ixNrgLayer,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixHydLayer); call populate_data_i(dataStruct%ixHydLayer,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
! case(iLookINDEX%ixVolFracWat); call populate_data_i(dataStruct%ixVolFracWat,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
! case(iLookINDEX%ixMatricHead); call populate_data_i(dataStruct%ixMatricHead,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixAllState); call populate_data_i(dataStruct%ixAllState,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixSoilState); call populate_data_i(dataStruct%ixSoilState,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
case(iLookINDEX%ixLayerState); call populate_data_i(dataStruct%ixLayerState,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
! case(iLookINDEX%ixLayerActive); call populate_data_i(dataStruct%ixLayerActive,tempVec,iSnow,removeLayer,nGRU,ix_lower,ix_upper)
  end select

  deallocate(tempVec,stat=err)
  if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

  
 end do  ! looping through variables

 end subroutine rmLyAllVars_index

 subroutine populate_data(struct, tempVec, iSnow, removeLayer, nGRU, ix_lower, ix_upper)
    real(rkind),device :: struct(ix_lower:,:)
  real(rkind),device :: tempVec(ix_lower:,:)
  integer(i4b),device :: iSnow(:)
  logical(lgt),device :: removeLayer(:)
  integer(i4b) :: ix_lower, ix_upper
  integer(i4b) :: nGRU

  integer(i4b) :: iGRU, iLayer
  !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    if (removeLayer(iGRU)) then
      if (iSnow(iGRU)>=ix_lower) tempVec(iSnow(iGRU),iGRU) = realMissing
      do iLayer=ix_lower,iSnow(iGRU)-1
        tempVec(iLayer,iGRU) = struct(iLayer,iGRU)
      end do
      do iLayer=iSnow(iGRU)+1,ix_upper-1
        tempVec(iLayer,iGRU) = struct(iLayer+1,iGRU)
      end do
      do iLayer=ix_lower,ix_upper-1
        struct(iLayer,iGRU) = tempVec(iLayer,iGRU)
      end do
    end if
  end do
end subroutine populate_data

 subroutine populate_data_i(struct, tempVec, iSnow, removeLayer, nGRU, ix_lower, ix_upper)
  integer(i4b),device :: struct(ix_lower:,:)
  integer(i4b),device :: tempVec(ix_lower:,:)
  integer(i4b),device :: iSnow(:)
  logical(lgt),device :: removeLayer(:)
  integer(i4b) :: ix_lower, ix_upper
  integer(i4b) :: nGRU

  integer(i4b) :: iGRU, iLayer
  !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    if (removeLayer(iGRU)) then
      if (iSnow(iGRU)>=ix_lower) tempVec(iSnow(iGRU),iGRU) = integerMissing
      do iLayer=ix_lower,iSnow(iGRU)-1
        tempVec(iLayer,iGRU) = struct(iLayer,iGRU)
      end do
      do iLayer=iSnow(iGRU)+1,ix_upper-1
        tempVec(iLayer,iGRU) = struct(iLayer+1,iGRU)
      end do
      do iLayer=ix_lower,ix_upper-1
        struct(iLayer,iGRU) = tempVec(iLayer,iGRU)
      end do
      ! do iLayer=ix_upper,ubound(struct,1)
        struct(ix_upper,iGRU) = integerMissing
      ! end do
        ! do iLayer=ix_lower,ix_upper
        !   print*, iGRU, iLayer, iSnow(iGRU), struct(iLayer,iGRU)
        ! end do
    end if
  end do
end subroutine populate_data_i


end module layerMerge_module
