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

module layerDivide_module

! variable types
USE nrtype

! physical constants
USE multiconst,only:&
                    iden_ice,       & ! intrinsic density of ice             (kg m-3)
                    iden_water        ! intrinsic density of liquid water    (kg m-3)

! access named variables for snow and soil
USE globalData,only:iname_snow        ! named variables for snow
USE globalData,only:iname_soil        ! named variables for soil

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing double precision number

! access metadata
USE globalData,only:prog_meta,diag_meta,flux_meta,indx_meta   ! metadata

! access the derived types to define the data structures
USE data_types,only:&
                    var_d,            & ! data vector (rkind)
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength,      & ! data vector with variable length dimension (rkind)
                    model_options       ! defines the model decisions

! access named variables defining elements in the data structures
USE var_lookup,only:iLookPROG,iLookDIAG,iLookFLUX,iLookINDEX  ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS                            ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM                                ! named variables for elements of the parameter structure

! define look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers,       & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex       ! CLM option: combination/sub-dividion rules depend on layer index

! define look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:&
 noah_mp,                  & ! full Noah-MP implementation (including albedo)
 CLM_2stream,              & ! CLM 2-stream model (see CLM documentation)
 UEB_2stream,              & ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
 NL_scatter,               & ! Simplified method Nijssen and Lettenmaier (JGR 1999)
 BeersLaw                    ! Beer's Law (as implemented in VIC)

! define look-up values for the choice of albedo method
USE mDecisions_module,only:& ! identify model options for snow albedo
 constantDecay,            & ! constant decay in snow albedo (e.g., VIC, CLASS)
 variableDecay               ! variable decay in snow albedo (e.g., BATS approach, with destructive metamorphism + soot content)

! privacy
implicit none
private
public::layerDivide

contains


 subroutine layerDivide(&
                        ! input/output: model data structures
                        decisions,                 & ! intent(in):    model decisions
                        nGRU, &
                        mpar_data,                       & ! intent(in):    model parameters
                        indx_data,                       & ! intent(inout): type of each layer
                        prog_data,                       & ! intent(inout): model prognostic variables for a local HRU
                        diag_data,                       & ! intent(inout): model diagnostic variables for a local HRU
                        flux_data,                       & ! intent(inout): model fluxes for a local HRU
                        ! output
                        divideLayer,                     & ! intent(out): flag to denote that a layer was divided
                        err,message)                       ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 ! computational modules
 USE snow_utils_module,only:fracliquid_d,templiquid_d          ! functions to compute temperature/liquid water
 use device_data_types
 USE globalData,only:maxSnowLayers, &                      ! maximum number of snow layers
                     veryBig, nBand
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input/output: model data structures
 type(decisions_device),intent(in)  :: decisions     ! model decisions
 integer(i4b) :: nGRU
 type(mpar_data_device),intent(in)    :: mpar_data              ! model parameters
 type(indx_data_device),intent(inout) :: indx_data              ! type of each layer
 type(prog_data_device),intent(inout) :: prog_data              ! model prognostic variables for a local HRU
 type(diag_data_device),intent(inout) :: diag_data              ! model diagnostic variables for a local HRU
 type(flux_data_device),intent(inout) :: flux_data              ! model flux variables
 ! output
 logical(lgt),intent(out)        :: divideLayer            ! flag to denote that a layer was divided
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! define local variables
 character(LEN=256)              :: cmessage               ! error message of downwind routine
 integer(i4b)                    :: nSnow                  ! number of snow layers
 integer(i4b)                    :: nSoil                  ! number of soil layers
!  integer(i4b)                    :: nLayers                ! total number of layers
 integer(i4b)                    :: iLayer                 ! layer index
 integer(i4b)                    :: jLayer                 ! layer index
 real(rkind),dimension(4),device        :: zmax_lower             ! lower value of maximum layer depth
 real(rkind),dimension(4),device        :: zmax_upper             ! upper value of maximum layer depth
 real(rkind)                     :: zmaxCheck              ! value of zmax for a given snow layer
 integer(i4b)                    :: nCheck                 ! number of layers to check to divide
 logical(lgt),device                    :: createLayer(nGRU)            ! flag to indicate we are creating a new snow layer
 real(rkind)                     :: depthOriginal          ! original layer depth before sub-division (m)
 real(rkind),parameter           :: fracTop=0.5_rkind      ! fraction of old layer used for the top layer
 real(rkind)                     :: surfaceLayerSoilTemp   ! temperature of the top soil layer (K)
 real(rkind)                     :: maxFrozenSnowTemp      ! maximum temperature when effectively all water is frozen (K)
 real(rkind),device           :: unfrozenLiq ! unfrozen liquid water used to compute maxFrozenSnowTemp (-)
 real(rkind)                     :: volFracWater           ! volumetric fraction of total water, liquid and ice (-)
 real(rkind)                     :: fracLiq                ! fraction of liquid water (-)
 integer(i4b),parameter          :: ixVisible=1            ! named variable to define index in array of visible part of the spectrum
 integer(i4b),parameter          :: ixNearIR=2             ! named variable to define index in array of near IR part of the spectrum
 real(rkind),parameter           :: verySmall=1.e-10_rkind ! a very small number (used for error checking)
 integer(i4b),device :: addedLayer(nGRU)
 integer(i4b) :: iGRU
 integer(i4b) :: nSoil2
!  integer(i4b) :: nLayers2
 unfrozenLiq=0.01_rkind
 ! --------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="layerDivide/"
 ! --------------------------------------------------------------------------------------------------------
 ! associate variables in the data structures
 associate(&
 ! model decisions
 ix_snowLayers          => decisions%snowLayers, & ! decision for snow combination
!  ! model parameters (compute layer temperature)
 fc_param               => mpar_data%snowfrz_scale,       & ! freezing curve parameter for snow (K-1)
!  ! model parameters (new snow density)
!  newSnowDenMin          => mpar_data%var(iLookPARAM%newSnowDenMin)%dat(1),       & ! minimum new snow density (kg m-3)
!  newSnowDenMult         => mpar_data%var(iLookPARAM%newSnowDenMult)%dat(1),      & ! multiplier for new snow density (kg m-3)
!  newSnowDenScal         => mpar_data%var(iLookPARAM%newSnowDenScal)%dat(1),      & ! scaling factor for new snow density (K)
!  ! model parameters (control the depth of snow layers)
 zmax                   => mpar_data%zmax,                & ! maximum layer depth (m)
!  zminLayer1             => mpar_data%var(iLookPARAM%zminLayer1)%dat(1),          & ! minimum layer depth for the 1st (top) layer (m)
 zmaxLayer1_lower       => mpar_data%zmaxLayer1_lower,    & ! maximum layer depth for the 1st (top) layer when only 1 layer (m)
 zmaxLayer2_lower       => mpar_data%zmaxLayer2_lower,    & ! maximum layer depth for the 2nd layer when only 2 layers (m)
 zmaxLayer3_lower       => mpar_data%zmaxLayer3_lower,    & ! maximum layer depth for the 3rd layer when only 3 layers (m)
 zmaxLayer4_lower       => mpar_data%zmaxLayer4_lower,    & ! maximum layer depth for the 4th layer when only 4 layers (m)
 zmaxLayer1_upper       => mpar_data%zmaxLayer1_upper,    & ! maximum layer depth for the 1st (top) layer when > 1 layer (m)
 zmaxLayer2_upper       => mpar_data%zmaxLayer2_upper,    & ! maximum layer depth for the 2nd layer when > 2 layers (m)
 zmaxLayer3_upper       => mpar_data%zmaxLayer3_upper,    & ! maximum layer depth for the 3rd layer when > 3 layers (m)
 zmaxLayer4_upper       => mpar_data%zmaxLayer4_upper,    & ! maximum layer depth for the 4th layer when > 4 layers (m)
!  ! diagnostic scalar variables
!  scalarSnowfall         => flux_data%var(iLookFLUX%scalarSnowfall)%dat(1),       & ! snowfall flux (kg m-2 s-1)
!  scalarSnowfallTemp     => diag_data%var(iLookDIAG%scalarSnowfallTemp)%dat(1),   & ! computed temperature of fresh snow (K)
 scalarSnowDepth        => prog_data%scalarSnowDepth,      & ! total snow depth (m)
 scalarSWE              => prog_data%scalarSWE,             & ! SWE (kg m-2)
 nSnow => indx_data%nSnow &
 )  ! end associate statement

 ! ---------------------------------------------------------------------------------------------------

 ! initialize flag to denote that a layer was divided
 divideLayer=.false.

 ! identify algorithmic control parameters to syb-divide and combine snow layers
 zmax_lower = (/zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower/)
 zmax_upper = (/zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper/)

 ! initialize the number of snow layers
 nSnow   = indx_data%nSnow
 nSoil   = indx_data%nSoil

 createLayer = .false.
 ! ***** special case of no snow layers
 associate(mLayerDepth => prog_data%mLayerDepth)
  !$cuf kernel do(1) <<<*,*>>>
 do iGRU=1,nGRU
 if(nSnow(iGRU)==0)then

  ! check if create the first snow layer
  select case(ix_snowLayers)
   case(sameRulesAllLayers);    createLayer(iGRU) = (scalarSnowDepth(iGRU) > zmax)
   case(rulesDependLayerIndex); createLayer(iGRU) = (scalarSnowDepth(iGRU) > zmax_lower(1))
  !  case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
  end select ! (option to combine/sub-divide snow layers)

  ! ** create a new snow layer
  if(createLayer(iGRU))then

   ! flag that the layers have changed
   divideLayer=.true.

   ! add a layer to all model variables
   iLayer=0 ! (layer to divide: 0 is the special case of "snow without a layer")
   addedLayer(iGRU) = 0
  endif

 ! end special case of nSnow=0
 ! ********************************************************************************************************************
 ! ********************************************************************************************************************

 ! ***** sub-divide snow layers, if necessary
 else ! if nSnow>0

  ! identify the number of layers to check for need for sub-division
  nCheck = min(nSnow(iGRU), maxSnowLayers-1) ! the depth of the last layer, if it exists, does not have a maximum value
  ! loop through all layers, and sub-divide a given layer, if necessary
  do iLayer=1,nCheck
   divideLayer=.false.

   ! identify the maximum depth of the layer
   select case(ix_snowLayers)
    case(sameRulesAllLayers)
     if (nCheck >= maxSnowLayers-1) then
      ! make sure we don't divide so make very big
      zmaxCheck = veryBig
     else
      zmaxCheck = zmax
     end if
    case(rulesDependLayerIndex)
     if(iLayer == nSnow(iGRU))then
      zmaxCheck = zmax_lower(iLayer)
     else
      zmaxCheck = zmax_upper(iLayer)
     end if
    ! case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
   end select ! (option to combine/sub-divide snow layers)

   ! check the need to sub-divide
   if(mLayerDepth(iLayer,iGRU) > zmaxCheck)then
    addedLayer(iGRU) = iLayer
    ! flag that layers were divided
    divideLayer=.true.
    createLayer(iGRU)=.true.

    exit  ! NOTE: only sub-divide one layer per substep

   end if   ! (if sub-dividing layer)

  end do  ! (looping through layers)

 end if  ! if nSnow==0
end do
end associate
 
 if (divideLayer) then
  ! add a layer to all model variables
  call addModelLayer_prog(prog_data,prog_meta,addedLayer,createLayer,nGRU,nSoil,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  call addModelLayer_diag(diag_data,diag_meta,addedLayer,createLayer,nGRU,nSoil,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  call addModelLayer_flux(flux_data,flux_meta,addedLayer,createLayer,nGRU,nSoil,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  call addModelLayer_indx(indx_data,indx_meta,addedLayer,createLayer,nGRU,nSoil,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 end if

  associate(mLayerDepth => prog_data%mLayerDepth, &
    mLayerTemp => prog_data%mLayerTemp, &
    mLayerVolFracIce => prog_data%mLayerVolFracIce, &
    mLayerVolFracLiq => prog_data%mLayerVolFracLiq, &
    canopySrad => decisions%canopySrad, &
    alb_method => decisions%alb_method, &
    scalarSnowAlbedo => prog_data%scalarSnowAlbedo, &
    albedoMax => mpar_data%albedoMax, &
    spectralSnowAlbedoDiffuse => prog_data%spectralSnowAlbedoDiffuse, &
    spectralSnowAlbedoDirect => diag_data%spectralSnowAlbedoDirect, &
    Frad_vis => mpar_data%Frad_vis, &
    albedoMaxNearIR => mpar_data%albedoMaxNearIR, &
    albedoMaxVisible => mpar_data%albedoMaxVisible)
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU
 if (createLayer(iGRU)) then
  if (nSnow(iGRU) == 0) then
!        ! associate local variables to the information in the data structures
!    ! NOTE: need to do this here, since state vectors have just been modified
!    associate(&
!     ! coordinate variables
!     mLayerDepth      => prog_data%var(iLookPROG%mLayerDepth)%dat        ,& ! depth of each layer (m)
!     ! model state variables
!     mLayerTemp       => prog_data%var(iLookPROG%mLayerTemp)%dat         ,& ! temperature of each layer (K)
!     mLayerVolFracIce => prog_data%var(iLookPROG%mLayerVolFracIce)%dat   ,& ! volumetric fraction of ice in each layer (-)
!     mLayerVolFracLiq => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat    & ! volumetric fraction of liquid water in each layer (-)
!     ) ! (association of local variables to the information in the data structures)
 
    ! get the layer depth
    mLayerDepth(1,iGRU) = scalarSnowDepth(iGRU)

    ! compute surface layer temperature
    surfaceLayerSoilTemp = mLayerTemp(2,iGRU)    ! temperature of the top soil layer (K)
    maxFrozenSnowTemp    = templiquid_d(unfrozenLiq,fc_param)               ! snow temperature at fraction "unfrozenLiq" (K)
    mLayerTemp(1,iGRU)        = min(maxFrozenSnowTemp,surfaceLayerSoilTemp)    ! snow temperature  (K)
 
    ! compute the fraction of liquid water associated with the layer temperature
    fracLiq = fracliquid_d(mLayerTemp(1,iGRU),fc_param)
    ! print*, 602, iGRU, createLayer(iGRU)

    ! compute volumeteric fraction of liquid water and ice
    volFracWater = (scalarSWE(iGRU)/scalarSnowDepth(iGRU))/iden_water  ! volumetric fraction of total water (liquid and ice)
    mLayerVolFracIce(1,iGRU) = (1._rkind - fracLiq)*volFracWater*(iden_water/iden_ice)   ! volumetric fraction of ice (-)
    mLayerVolFracLiq(1,iGRU) =             fracLiq *volFracWater                         ! volumetric fraction of liquid water (-)
    ! print*, 608, iGRU, createLayer(iGRU)

!     ! end association with local variables to the information in the data structures)
!     end associate
    ! initialize albedo
    ! NOTE: albedo is computed within the Noah-MP radiation routine
    if(canopySrad /= noah_mp)then
     select case(alb_method)
      ! (constant decay rate -- albedo the same for all spectral bands)
      case(constantDecay)
       scalarSnowAlbedo(iGRU)          = albedoMax
      !  print*, 619, iGRU, createLayer(iGRU)

       do iLayer=1,nBand
       spectralSnowAlbedoDiffuse(iLayer,iGRU) = albedoMax
       end do
      ! (variable decay rate)
      case(variableDecay)
       spectralSnowAlbedoDiffuse(ixVisible,iGRU) = albedoMaxVisible
       spectralSnowAlbedoDiffuse(ixNearIR,iGRU)  = albedoMaxNearIR
       scalarSnowAlbedo(iGRU)                  = (        Frad_vis)*albedoMaxVisible + &
                                              (1._rkind - Frad_vis)*albedoMaxNearIR
!       case default; err=20; message=trim(message)//'unable to identify option for snow albedo'; return
     end select  ! identify option for snow albedo
     ! set direct albedo to diffuse albedo
    !  print*, 633, iGRU, createLayer(iGRU)

     do iLayer=1,nBand
     spectralSnowAlbedoDirect(iLayer,iGRU) = spectralSnowAlbedoDiffuse(iLayer,iGRU)
     end do
    end if  ! (if NOT using the Noah-MP radiation routine)
  ! end if
! end if

  else   
    ! print*, 643, iGRU, createLayer(iGRU)

!         ! define the layer depth
!     layerSplit: associate(mLayerDepth => prog_data%var(iLookPROG%mLayerDepth)%dat)
    depthOriginal = mLayerDepth(addedLayer(iGRU),iGRU)
    mLayerDepth(addedLayer(iGRU),iGRU)   = fracTop*depthOriginal
    mLayerDepth(addedLayer(iGRU)+1,iGRU) = (1._rkind - fracTop)*depthOriginal
!     end associate layerSplit

  endif 
endif
end do
end associate


 ! update coordinates
 if(divideLayer)then

  ! associate coordinate variables in data structure
  geometry: associate(&
  mLayerDepth      => prog_data%mLayerDepth        ,& ! depth of the layer (m)
  mLayerHeight     => prog_data%mLayerHeight       ,& ! height of the layer mid-point (m)
  iLayerHeight     => prog_data%iLayerHeight       ,& ! height of the layer interface (m)
  layerType        => indx_data%layerType         ,& ! type of each layer (iname_snow or iname_soil)
  nSnow            => indx_data%nSnow          ,& ! number of snow layers
  ! nSoil            => indx_data%nSoil          ,& ! number of soil layers
  ! nSoil_d => indx_data%nSoil_d, &
  nLayers          => indx_data%nLayers_d         & ! total number of layers
  )  ! (association of local variables with coordinate variab;es in data structures)

  ! update the layer type
  !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    do iLayer=1,nLayers(iGRU)+1
      ! print*, 676, iGRU, createLayer(iGRU)

      if (createLayer(iGRU)) then
      if (iLayer .le. nSnow(iGRU)+1) layerType(iLayer,iGRU) = iname_snow
      if (iLayer .ge. nSnow(iGRU) + 2) layerType(iLayer,iGRU) = iname_soil
    end if

    end do
  end do
  ! layerType(1:nSnow+1)         = iname_snow
  ! layerType(nSnow+2:nLayers+1) = iname_soil

  ! identify the number of snow and soil layers, and check all is a-OK
  nSoil2 = nSoil
  ! indx_data%nSoil_d = nSoil

  !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    nSnow(iGRU) = 0
    do iLayer=1,size(layerType,1)
      if (layerType(iLayer,iGRU) == iname_snow) nSnow(iGRU) = nSnow(iGRU) + 1
    end do
  !   if (createLayer(iGRU)) then

    nLayers(iGRU) = nSoil + nSnow(iGRU)
  !     print*, 694, iGRU
  ! nSnow(iGRU)   = nSnow(iGRU) + 1
  !   end if
  end do

  ! re-set coordinate variables
  !$cuf kernel do(1) <<<*,*>>>
  do iGRU=1,nGRU
    if (createLayer(iGRU)) then

  iLayerHeight(0,iGRU) = -scalarSnowDepth(iGRU)
  do jLayer=1,nLayers(iGRU)
   iLayerHeight(jLayer,iGRU) = iLayerHeight(jLayer-1,iGRU) + mLayerDepth(jLayer,iGRU)
   mLayerHeight(jLayer,iGRU) = (iLayerHeight(jLayer-1,iGRU) + iLayerHeight(jLayer,iGRU))/2._rkind
  end do
end if
end do

!   ! check
!   if(abs(sum(mLayerDepth(1:nSnow)) - scalarSnowDepth) > verySmall)then
!    print*, 'nSnow = ', nSnow
!    write(*,'(a,1x,f30.25,1x)') 'sum(mLayerDepth(1:nSnow)) = ', sum(mLayerDepth(1:nSnow))
!    write(*,'(a,1x,f30.25,1x)') 'scalarSnowDepth           = ', scalarSnowDepth
!    write(*,'(a,1x,f30.25,1x)') 'epsilon(scalarSnowDepth)  = ', epsilon(scalarSnowDepth)
!    message=trim(message)//'sum of layer depths does not equal snow depth'
!    err=20; return
!   end if

  ! end association with coordinate variables in data structure
  end associate geometry

 end if  ! if dividing a layer

 ! end associate variables in data structure
 end associate

 end subroutine layerDivide


 ! ************************************************************************************************
 ! private subroutine addModelLayer: add an additional layer to all model vectors
 ! ************************************************************************************************
 subroutine addModelLayer(dataStruct,metaStruct,ix_divide,nSnow,nLayers,err,message)
 USE var_lookup,only:iLookVarType                     ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName           ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc                ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength          ! data vectors with variable length dimension
 USE data_types,only:var_info                         ! metadata structure
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 class(*),intent(inout)          :: dataStruct        ! data structure
 type(var_info),intent(in)       :: metaStruct(:)     ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in)         :: ix_divide         ! index of the layer to divide
 integer(i4b),intent(in)         :: nSnow,nLayers     ! number of snow layers, total number of layers
 ! output: error control
 integer(i4b),intent(out)        :: err               ! error code
 character(*),intent(out)        :: message           ! error message
 ! ---------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: ivar              ! index of model variable
 integer(i4b)                    :: ix_lower          ! lower bound of the vector
 integer(i4b)                    :: ix_upper          ! upper bound of the vector
 logical(lgt)                    :: stateVariable     ! .true. if variable is a state variable
 real(rkind),allocatable         :: tempVec_rkind(:)  ! temporary vector (double precision)
 integer(i4b),allocatable        :: tempVec_i4b(:)    ! temporary vector (integer)
 character(LEN=256)              :: cmessage          ! error message of downwind routine
 ! ---------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='addModelLayer/'

 ! ***** add a layer to each model variable
 do ivar=1,size(metaStruct)

  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=nSnow
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=nSnow
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nLayers
   case default; cycle
  end select

  ! identify whether it is a state variable
  select case(trim(metaStruct(ivar)%varname))
   case('mLayerDepth','mLayerTemp','mLayerVolFracIce','mLayerVolFracLiq'); stateVariable=.true.
   case default; stateVariable=.false.
  end select

  ! divide layers
  select type(dataStruct)

   ! ** double precision
   type is (var_dlength)
    ! check allocated
    if(.not.allocated(dataStruct%var(ivar)%dat))then; err=20; message='data vector is not allocated'; return; end if
    ! assign the data vector to the temporary vector
    call cloneStruc(tempVec_rkind, ix_lower, source=dataStruct%var(ivar)%dat, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! reallocate space for the new vector
    deallocate(dataStruct%var(ivar)%dat,stat=err)
    if(err/=0)then; err=20; message='problem in attempt to deallocate memory for data vector'; return; end if
    allocate(dataStruct%var(ivar)%dat(ix_lower:ix_upper+1),stat=err)
    if(err/=0)then; err=20; message='problem in attempt to reallocate memory for data vector'; return; end if
    ! populate the state vector
    if(stateVariable)then
     if(ix_upper > 0)then  ! (only copy data if the vector exists -- can be a variable for snow, with no layers)
      if(ix_divide > 0)then
       dataStruct%var(ivar)%dat(1:ix_divide) = tempVec_rkind(1:ix_divide)  ! copy data
       dataStruct%var(ivar)%dat(ix_divide+1) = tempVec_rkind(ix_divide)    ! repeat data for the sub-divided layer
      end if
      if(ix_upper > ix_divide) &
       dataStruct%var(ivar)%dat(ix_divide+2:ix_upper+1) = tempVec_rkind(ix_divide+1:ix_upper)  ! copy data
     end if  ! if the vector exists
    ! not a state variable
    else
     dataStruct%var(ivar)%dat(:) = realMissing
    end if
    ! deallocate the temporary vector: strictly not necessary, but include to be safe
    deallocate(tempVec_rkind,stat=err)
    if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

   ! ** integer
   type is (var_ilength)
    ! check allocated
    if(.not.allocated(dataStruct%var(ivar)%dat))then; err=20; message='data vector is not allocated'; return; end if
    ! assign the data vector to the temporary vector
    call cloneStruc(tempVec_i4b, ix_lower, source=dataStruct%var(ivar)%dat, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    ! reallocate space for the new vector
    deallocate(dataStruct%var(ivar)%dat,stat=err)
    if(err/=0)then; err=20; message='problem in attempt to deallocate memory for data vector'; return; end if
    allocate(dataStruct%var(ivar)%dat(ix_lower:ix_upper+1),stat=err)
    if(err/=0)then; err=20; message='problem in attempt to reallocate memory for data vector'; return; end if
    ! populate the state vector
    if(stateVariable)then
     if(ix_upper > 0)then  ! (only copy data if the vector exists -- can be a variable for snow, with no layers)
      if(ix_divide > 0)then
       dataStruct%var(ivar)%dat(1:ix_divide) = tempVec_i4b(1:ix_divide)  ! copy data
       dataStruct%var(ivar)%dat(ix_divide+1) = tempVec_i4b(ix_divide)    ! repeat data for the sub-divided layer
      end if
      if(ix_upper > ix_divide) &
       dataStruct%var(ivar)%dat(ix_divide+2:ix_upper+1) = tempVec_i4b(ix_divide+1:ix_upper)  ! copy data
     end if  ! if the vector exists
    ! not a state variable
    else
     dataStruct%var(ivar)%dat(:) = integerMissing
    end if
    ! deallocate the temporary vector: strictly not necessary, but include to be safe
    deallocate(tempVec_i4b,stat=err)
    if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if

   ! check that we found the data type
   class default; err=20; message=trim(message)//'unable to identify the data type'; return

  end select ! dependence on data types

 end do  ! looping through variables

 end subroutine addModelLayer

subroutine addModelLayer_prog(dataStruct,metaStruct,ix_divide,createLayer,nGRU,nSoil,err,message)
  use device_data_types
 USE var_lookup,only:iLookVarType                     ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName           ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc                ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength          ! data vectors with variable length dimension
 USE data_types,only:var_info                         ! metadata structure
 use globalData,only:maxSnowLayers
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 type(prog_data_device),intent(inout)          :: dataStruct        ! data structure
 type(var_info),intent(in)       :: metaStruct(:)     ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in),device         :: ix_divide(:)         ! index of the layer to divide
 logical(lgt),intent(in),device :: createLayer(:)
 integer(i4b),intent(in)         :: nGRU,nSoil     ! number of snow layers, total number of layers
 ! output: error control
 integer(i4b),intent(out)        :: err               ! error code
 character(*),intent(out)        :: message           ! error message
 ! ---------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: ivar              ! index of model variable
 integer(i4b)                    :: ix_lower          ! lower bound of the vector
 integer(i4b)                    :: ix_upper          ! upper bound of the vector
 logical(lgt)                    :: stateVariable     ! .true. if variable is a state variable
 real(rkind),allocatable,device         :: tempVec(:,:)  ! temporary vector (double precision)
!  integer(i4b),allocatable        :: tempVec_i4b(:)    ! temporary vector (integer)
 character(LEN=256)              :: cmessage          ! error message of downwind routine
 ! ---------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='addModelLayer/'

 ! ***** add a layer to each model variable
 do ivar=1,size(metaStruct)

  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=maxSnowLayers
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nSoil+maxSnowLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=maxSnowLayers
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nSoil+maxSnowLayers
   case default; cycle
  end select

  ! identify whether it is a state variable
  select case(trim(metaStruct(ivar)%varname))
   case('mLayerDepth','mLayerTemp','mLayerVolFracIce','mLayerVolFracLiq'); stateVariable=.true.
   case default; stateVariable=.false.
  end select


  ! divide layers
  ! assign the data vector to the temporary vector
  select case(iVar)
    case(iLookPROG%spectralSnowAlbedoDiffuse); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%spectralSnowAlbedoDiffuse,err=err,message=cmessage)
    case(iLookPROG%mLayerTemp); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerTemp,err=err,message=cmessage)
    case(iLookPROG%mLayerVolFracIce); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerVolFracIce,err=err,message=cmessage)
    case(iLookPROG%mLayerVolFracLiq); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerVolFracLiq,err=err,message=cmessage)
    case(iLookPROG%mLayerVolFracWat); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerVolFracWat,err=err,message=cmessage)
    case(iLookPROG%mLayerMatricHead); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerMatricHead,err=err,message=cmessage)
    case(iLookPROG%mLayerEnthalpy); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerEnthalpy,err=err,message=cmessage)
    case(iLookPROG%mLayerDepth); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerDepth,err=err,message=cmessage)
    case(iLookPROG%mLayerHeight); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerHeight,err=err,message=cmessage)
    case(iLookPROG%iLayerHeight); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%iLayerHeight,err=err,message=cmessage)
    end select

  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  ! reallocate space for the new vector
  ! populate the state vector
  select case(iVar)
    case(iLookPROG%spectralSnowAlbedoDiffuse); call populate_data(dataStruct%spectralSnowAlbedoDiffuse,tempVec, ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
    case(iLookPROG%mLayerTemp); call populate_data(dataStruct%mLayerTemp,tempVec, ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
    case(iLookPROG%mLayerVolFracIce); call populate_data(dataStruct%mLayerVolFracIce,tempVec, ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
    case(iLookPROG%mLayerVolFracLiq); call populate_data(dataStruct%mLayerVolFracLiq,tempVec, ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
    case(iLookPROG%mLayerVolFracWat); call populate_data(dataStruct%mLayerVolFracWat,tempVec, ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
    case(iLookPROG%mLayerMatricHead); call populate_data(dataStruct%mLayerMatricHead,tempVec, ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
    case(iLookPROG%mLayerEnthalpy); call populate_data(dataStruct%mLayerEnthalpy,tempVec, ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
    case(iLookPROG%mLayerDepth); call populate_data(dataStruct%mLayerDepth,tempVec, ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
    case(iLookPROG%mLayerHeight); call populate_data(dataStruct%mLayerHeight,tempVec, ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
    case(iLookPROG%iLayerHeight); call populate_data(dataStruct%iLayerHeight,tempVec, ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
    end select
  ! deallocate the temporary vector: strictly not necessary, but include to be safe
  deallocate(tempVec,stat=err)
  if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if


 end do  ! looping through variables

 end subroutine addModelLayer_prog

 subroutine addModelLayer_diag(dataStruct,metaStruct,ix_divide,createLayer,nGRU,nSoil,err,message)
  use device_data_types
 USE var_lookup,only:iLookVarType                     ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName           ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc                ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength          ! data vectors with variable length dimension
 USE data_types,only:var_info                         ! metadata structure
 use globalData,only:maxSnowLayers
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 type(diag_data_device),intent(inout)          :: dataStruct        ! data structure
 type(var_info),intent(in)       :: metaStruct(:)     ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in),device         :: ix_divide(:)         ! index of the layer to divide
 logical(lgt),intent(in),device :: createLayer(:)
 integer(i4b),intent(in)         :: nGRU,nSOil     ! number of snow layers, total number of layers
 ! output: error control
 integer(i4b),intent(out)        :: err               ! error code
 character(*),intent(out)        :: message           ! error message
 ! ---------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: ivar              ! index of model variable
 integer(i4b)                    :: ix_lower          ! lower bound of the vector
 integer(i4b)                    :: ix_upper          ! upper bound of the vector
 logical(lgt)                    :: stateVariable     ! .true. if variable is a state variable
 real(rkind),allocatable,device         :: tempVec(:,:)  ! temporary vector (double precision)
!  integer(i4b),allocatable        :: tempVec_i4b(:)    ! temporary vector (integer)
 character(LEN=256)              :: cmessage          ! error message of downwind routine
 ! ---------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='addModelLayer/'

 ! ***** add a layer to each model variable
 do ivar=1,size(metaStruct)

  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=maxSnowLayers
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nSoil+maxSnowLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=maxSnowLayers
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nSoil+maxSnowLayers
   case default; cycle
  end select

  ! identify whether it is a state variable
  select case(trim(metaStruct(ivar)%varname))
   case('mLayerDepth','mLayerTemp','mLayerVolFracIce','mLayerVolFracLiq'); stateVariable=.true.
   case default; stateVariable=.false.
  end select

  ! divide layers
  ! assign the data vector to the temporary vector
  select case(iVar)
  case(iLookDIAG%mLayerVolHtCapBulk); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerVolHtCapBulk_m,err=err,message=cmessage)
  case(iLookDIAG%mLayerCm); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerCm_m,err=err,message=cmessage)
  case(iLookDIAG%mLayerThermalC); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerThermalC_m,err=err,message=cmessage)
  case(iLookDIAG%iLayerThermalC); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%iLayerThermalC_m,err=err,message=cmessage)
  case(iLookDIAG%mLayerEnthTemp); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerEnthTemp,err=err,message=cmessage)
  case(iLookDIAG%spectralAlbGndDirect); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%spectralAlbGndDirect,err=err,message=cmessage)
  case(iLookDIAG%spectralAlbGndDiffuse); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%spectralAlbGndDiffuse,err=err,message=cmessage)
  case(iLookDIAG%mLayerTranspireLim); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerTranspireLim_m,err=err,message=cmessage)
  case(iLookDIAG%mLayerRootDensity); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerRootDensity_m,err=err,message=cmessage)
  case(iLookDIAG%spectralSnowAlbedoDirect); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%spectralSnowAlbedoDirect,err=err,message=cmessage)
  case(iLookDIAG%mLayerFracLiqSnow); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerFracLiqSnow_m,err=err,message=cmessage)
  case(iLookDIAG%mLayerThetaResid); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerThetaResid_m,err=err,message=cmessage)
  case(iLookDIAG%mLayerPoreSpace); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerPoreSpace_m,err=err,message=cmessage)
  case(iLookDIAG%mLayerMeltFreeze); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerMeltFreeze,err=err,message=cmessage)
  case(iLookDIAG%mLayerVolFracAir); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerVolFracAir_m,err=err,message=cmessage)
  case(iLookDIAG%mLayerTcrit); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerTcrit,err=err,message=cmessage)
  case(iLookDIAG%mLayerCompress); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerCompress_m,err=err,message=cmessage)
  case(iLookDIAG%mLayerMatricHeadLiq); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%mLayerMatricHeadLiq,err=err,message=cmessage)
  case(iLookDIAG%scalarVGn_m); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%scalarVGn_m_m,err=err,message=cmessage)
  case(iLookDIAG%balanceLayerNrg); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%balanceLayerNrg,err=err,message=cmessage)
  case(iLookDIAG%balanceLayerMass); call cloneStruc(tempVec, ix_lower,nGRU, source=dataStruct%balanceLayerMass,err=err,message=cmessage)
  end select

  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  ! reallocate space for the new vector
  ! populate the state vector
  select case(iVar)
  case(iLookDIAG%mLayerVolHtCapBulk); call populate_data(dataStruct%mLayerVolHtCapBulk_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerCm); call populate_data(dataStruct%mLayerCm_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerThermalC); call populate_data(dataStruct%mLayerThermalC_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%iLayerThermalC); call populate_data(dataStruct%iLayerThermalC_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerEnthTemp); call populate_data(dataStruct%mLayerEnthTemp,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%spectralAlbGndDirect); call populate_data(dataStruct%spectralAlbGndDirect,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%spectralAlbGndDiffuse); call populate_data(dataStruct%spectralAlbGndDiffuse,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerTranspireLim); call populate_data(dataStruct%mLayerTranspireLim_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerRootDensity); call populate_data(dataStruct%mLayerRootDensity_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%spectralSnowAlbedoDirect); call populate_data(dataStruct%spectralSnowAlbedoDirect,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerFracLiqSnow); call populate_data(dataStruct%mLayerFracLiqSnow_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerThetaResid); call populate_data(dataStruct%mLayerThetaResid_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerPoreSpace); call populate_data(dataStruct%mLayerPoreSpace_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerMeltFreeze); call populate_data(dataStruct%mLayerMeltFreeze,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerVolFracAir); call populate_data(dataStruct%mLayerVolFracAir_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerTcrit); call populate_data(dataStruct%mLayerTcrit,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerCompress); call populate_data(dataStruct%mLayerCompress_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%mLayerMatricHeadLiq); call populate_data(dataStruct%mLayerMatricHeadLiq,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%scalarVGn_m); call populate_data(dataStruct%scalarVGn_m_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%balanceLayerNrg); call populate_data(dataStruct%balanceLayerNrg,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookDIAG%balanceLayerMass); call populate_data(dataStruct%balanceLayerMass,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
    end select
  ! deallocate the temporary vector: strictly not necessary, but include to be safe
  deallocate(tempVec,stat=err)
  if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if


 end do  ! looping through variables

 end subroutine addModelLayer_diag

 subroutine addModelLayer_flux(dataStruct,metaStruct,ix_divide,createLayer,nGRU,nSoil,err,message)
  use device_data_types
 USE var_lookup,only:iLookVarType                     ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName           ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc                ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength          ! data vectors with variable length dimension
 USE data_types,only:var_info                         ! metadata structure
 use globalData,only:maxSnowLayers
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 type(flux_data_device),intent(inout)          :: dataStruct        ! data structure
 type(var_info),intent(in)       :: metaStruct(:)     ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in),device         :: ix_divide(:)         ! index of the layer to divide
 logical(lgt),intent(in),device :: createLayer(:)
 integer(i4b),intent(in)         :: nGRU,nSoil     ! number of snow layers, total number of layers
 ! output: error control
 integer(i4b),intent(out)        :: err               ! error code
 character(*),intent(out)        :: message           ! error message
 ! ---------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: ivar              ! index of model variable
 integer(i4b)                    :: ix_lower          ! lower bound of the vector
 integer(i4b)                    :: ix_upper          ! upper bound of the vector
 logical(lgt)                    :: stateVariable     ! .true. if variable is a state variable
 real(rkind),allocatable,device         :: tempVec(:,:)  ! temporary vector (double precision)
!  integer(i4b),allocatable        :: tempVec_i4b(:)    ! temporary vector (integer)
 character(LEN=256)              :: cmessage          ! error message of downwind routine
 ! ---------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='addModelLayer/'

 ! ***** add a layer to each model variable
 do ivar=1,size(metaStruct)

  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=maxSnowLayers
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nSoil+maxSnowLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=maxSnowLayers
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nSoil+maxSnowLayers
   case default; cycle
  end select

  ! identify whether it is a state variable
  select case(trim(metaStruct(ivar)%varname))
   case('mLayerDepth','mLayerTemp','mLayerVolFracIce','mLayerVolFracLiq'); stateVariable=.true.
   case default; stateVariable=.false.
  end select

  ! divide layers
  ! assign the data vector to the temporary vector
  select case(iVar)
  case(iLookFLUX%spectralIncomingDirect); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%spectralIncomingDirect,err=err,message=cmessage)
  case(iLookFLUX%spectralIncomingDiffuse); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%spectralIncomingDiffuse,err=err,message=cmessage)
  case(iLookFLUX%spectralBelowCanopyDirect); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%spectralBelowCanopyDirect,err=err,message=cmessage)
  case(iLookFLUX%spectralBelowCanopyDiffuse); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%spectralBelowCanopyDiffuse,err=err,message=cmessage)
  case(iLookFLUX%mLayerTranspire); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%mLayerTranspire_m,err=err,message=cmessage)
  case(iLookFLUX%iLayerConductiveFlux); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%iLayerConductiveFlux_m,err=err,message=cmessage)
  case(iLookFLUX%iLayerAdvectiveFlux); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%iLayerAdvectiveFlux_m,err=err,message=cmessage)
  case(iLookFLUX%iLayerNrgFlux); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%iLayerNrgFlux_m,err=err,message=cmessage)
  case(iLookFLUX%mLayerNrgFlux); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%mLayerNrgFlux_m,err=err,message=cmessage)
  case(iLookFLUX%iLayerLiqFluxSnow); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%iLayerLiqFluxSnow_m,err=err,message=cmessage)
  case(iLookFLUX%mLayerLiqFluxSnow); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%mLayerLiqFluxSnow_m,err=err,message=cmessage)
  case(iLookFLUX%mLayerSatHydCondMP); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%mLayerSatHydCondMP_m,err=err,message=cmessage)
  case(iLookFLUX%mLayerSatHydCond); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%mLayerSatHydCond_m,err=err,message=cmessage)
  case(iLookFLUX%iLayerSatHydCond); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%iLayerSatHydCond_m,err=err,message=cmessage)
  case(iLookFLUX%mLayerHydCond); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%mLayerHydCond_m,err=err,message=cmessage)
  case(iLookFLUX%iLayerLiqFluxSoil); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%iLayerLiqFluxSoil_m,err=err,message=cmessage)
  case(iLookFLUX%mLayerLiqFluxSoil); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%mLayerLiqFluxSoil_m,err=err,message=cmessage)
  case(iLookFLUX%mLayerBaseflow); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%mLayerBaseflow_m,err=err,message=cmessage)
  case(iLookFLUX%mLayerColumnInflow); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%mLayerColumnInflow,err=err,message=cmessage)
  case(iLookFLUX%mLayerColumnOutflow); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%mLayerColumnOutflow_m,err=err,message=cmessage)

    end select

  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  ! reallocate space for the new vector
  ! populate the state vector
  select case(iVar)
  case(iLookFLUX%spectralIncomingDirect); call populate_data(dataStruct%spectralIncomingDirect,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%spectralIncomingDiffuse); call populate_data(dataStruct%spectralIncomingDiffuse,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%spectralBelowCanopyDirect); call populate_data(dataStruct%spectralBelowCanopyDirect,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%spectralBelowCanopyDiffuse); call populate_data(dataStruct%spectralBelowCanopyDiffuse,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%mLayerTranspire); call populate_data(dataStruct%mLayerTranspire_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%iLayerConductiveFlux); call populate_data(dataStruct%iLayerConductiveFlux_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%iLayerAdvectiveFlux); call populate_data(dataStruct%iLayerAdvectiveFlux_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%iLayerNrgFlux); call populate_data(dataStruct%iLayerNrgFlux_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%mLayerNrgFlux); call populate_data(dataStruct%mLayerNrgFlux_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%iLayerLiqFluxSnow); call populate_data(dataStruct%iLayerLiqFluxSnow_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%mLayerLiqFluxSnow); call populate_data(dataStruct%mLayerLiqFluxSnow_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%mLayerSatHydCondMP); call populate_data(dataStruct%mLayerSatHydCondMP_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%mLayerSatHydCond); call populate_data(dataStruct%mLayerSatHydCond_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%iLayerSatHydCond); call populate_data(dataStruct%iLayerSatHydCond_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%mLayerHydCond); call populate_data(dataStruct%mLayerHydCond_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%iLayerLiqFluxSoil); call populate_data(dataStruct%iLayerLiqFluxSoil_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%mLayerLiqFluxSoil); call populate_data(dataStruct%mLayerLiqFluxSoil_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%mLayerBaseflow); call populate_data(dataStruct%mLayerBaseflow_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%mLayerColumnInflow); call populate_data(dataStruct%mLayerColumnInflow,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
  case(iLookFLUX%mLayerColumnOutflow); call populate_data(dataStruct%mLayerColumnOutflow_m,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
    end select
  ! deallocate the temporary vector: strictly not necessary, but include to be safe
  deallocate(tempVec,stat=err)
  if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if


 end do  ! looping through variables

 end subroutine addModelLayer_flux

 subroutine addModelLayer_indx(dataStruct,metaStruct,ix_divide,createLayer,nGRU,nSoil,err,message)
  use device_data_types
 USE var_lookup,only:iLookVarType                     ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName           ! to access type strings for error messages
 USE f2008funcs_module,only:cloneStruc                ! used to "clone" data structures -- temporary replacement of the intrinsic allocate(a, source=b)
 USE data_types,only:var_ilength,var_dlength          ! data vectors with variable length dimension
 USE data_types,only:var_info                         ! metadata structure
 use globalData,only:maxSnowLayers
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 type(indx_data_device),intent(inout)          :: dataStruct        ! data structure
 type(var_info),intent(in)       :: metaStruct(:)     ! metadata structure
 ! input: snow layer indices
 integer(i4b),intent(in),device         :: ix_divide(:)         ! index of the layer to divide
 logical(lgt),intent(in),device :: createLayer(:)
 integer(i4b),intent(in)         :: nGRU,nSoil     ! number of snow layers, total number of layers
 ! output: error control
 integer(i4b),intent(out)        :: err               ! error code
 character(*),intent(out)        :: message           ! error message
 ! ---------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: ivar              ! index of model variable
 integer(i4b)                    :: ix_lower          ! lower bound of the vector
 integer(i4b)                    :: ix_upper          ! upper bound of the vector
 logical(lgt)                    :: stateVariable     ! .true. if variable is a state variable
 integer(i4b),allocatable,device         :: tempVec(:,:)  ! temporary vector (integer)
 character(LEN=256)              :: cmessage          ! error message of downwind routine
 ! ---------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='addModelLayer/'

 ! ***** add a layer to each model variable
 do ivar=1,size(metaStruct)

  ! define bounds
  select case(metaStruct(ivar)%vartype)
   case(iLookVarType%midSnow); ix_lower=1; ix_upper=maxSnowLayers
   case(iLookVarType%midToto); ix_lower=1; ix_upper=nSoil+maxSnowLayers
   case(iLookVarType%ifcSnow); ix_lower=0; ix_upper=maxSnowLayers
   case(iLookVarType%ifcToto); ix_lower=0; ix_upper=nSoil+maxSnowLayers
   case default; cycle
  end select

  ! identify whether it is a state variable
  select case(trim(metaStruct(ivar)%varname))
   case('mLayerDepth','mLayerTemp','mLayerVolFracIce','mLayerVolFracLiq'); stateVariable=.true.
   case default; stateVariable=.false.
  end select

  ! print*, iVar, trim(metaStruct(iVar)%varname)
  ! divide layers
  ! assign the data vector to the temporary vector
  select case(iVar)
  case(iLookINDEX%layerType); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%layerType,err=err,message=cmessage)
  case(iLookINDEX%ixControlVolume); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixControlVolume,err=err,message=cmessage)
  case(iLookINDEX%ixDomainType); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixDomainType,err=err,message=cmessage)
  case(iLookINDEX%ixStateType); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixStateType,err=err,message=cmessage)
  case(iLookINDEX%ixHydType); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixHydType,err=err,message=cmessage)
  ! case(iLookINDEX%ixDomainType_subset); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixDomainType_subset,err=err,message=cmessage)
  ! case(iLookINDEX%ixStateType_subset); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixStateType_subset,err=err,message=cmessage)
  ! case(iLookINDEX%ixMapFull2Subset); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixMapFull2Subset,err=err,message=cmessage)
  ! case(iLookINDEX%ixMapSubset2Full); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixMapSubset2Full,err=err,message=cmessage)
  case(iLookINDEX%ixNrgOnly); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixNrgOnly,err=err,message=cmessage)
  ! case(iLookINDEX%ixHydOnly); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixHydOnly,err=err,message=cmessage)
  ! case(iLookINDEX%ixMatOnly); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixMatOnly,err=err,message=cmessage)
  ! case(iLookINDEX%ixMassOnly); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixMassOnly,err=err,message=cmessage)
  case(iLookINDEX%ixSnowSoilNrg); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixSnowSoilNrg,err=err,message=cmessage)
  case(iLookINDEX%ixSnowOnlyNrg); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixSnowOnlyNrg,err=err,message=cmessage)
  case(iLookINDEX%ixSoilOnlyNrg); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixSoilOnlyNrg,err=err,message=cmessage)
  case(iLookINDEX%ixSnowSoilHyd); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixSnowSoilHyd,err=err,message=cmessage)
  case(iLookINDEX%ixSnowOnlyHyd); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixSnowOnlyHyd,err=err,message=cmessage)
  case(iLookINDEX%ixSoilOnlyHyd); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixSoilOnlyHyd,err=err,message=cmessage)
  case(iLookINDEX%ixNrgLayer); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixNrgLayer,err=err,message=cmessage)
  case(iLookINDEX%ixHydLayer); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixHydLayer,err=err,message=cmessage)
  ! case(iLookINDEX%ixVolFracWat); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixVolFracWat,err=err,message=cmessage)
  ! case(iLookINDEX%ixMatricHead); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixMatricHead,err=err,message=cmessage)
  case(iLookINDEX%ixAllState); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixAllState,err=err,message=cmessage)
  case(iLookINDEX%ixSoilState); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixSoilState,err=err,message=cmessage)
  case(iLookINDEX%ixLayerState); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixLayerState,err=err,message=cmessage)
  ! case(iLookINDEX%ixLayerActive); call cloneStruc(tempVec,ix_lower,nGRU,source=dataStruct%ixLayerActive,err=err,message=cmessage)
  
    end select

  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  ! reallocate space for the new vector
  ! populate the state vector
  select case(iVar)
case(iLookINDEX%layerType); call populate_data_i(dataStruct%layerType,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixControlVolume); call populate_data_i(dataStruct%ixControlVolume,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixDomainType); call populate_data_i(dataStruct%ixDomainType,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixStateType); call populate_data_i(dataStruct%ixStateType,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixHydType); call populate_data_i(dataStruct%ixHydType,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
! case(iLookINDEX%ixDomainType_subset); call populate_data_i(dataStruct%ixDomainType_subset,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
! case(iLookINDEX%ixStateType_subset); call populate_data_i(dataStruct%ixStateType_subset,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
! case(iLookINDEX%ixMapFull2Subset); call populate_data_i(dataStruct%ixMapFull2Subset,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
! case(iLookINDEX%ixMapSubset2Full); call populate_data_i(dataStruct%ixMapSubset2Full,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixNrgOnly); call populate_data_i(dataStruct%ixNrgOnly,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
! case(iLookINDEX%ixHydOnly); call populate_data_i(dataStruct%ixHydOnly,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
! case(iLookINDEX%ixMatOnly); call populate_data_i(dataStruct%ixMatOnly,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
! case(iLookINDEX%ixMassOnly); call populate_data_i(dataStruct%ixMassOnly,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixSnowSoilNrg); call populate_data_i(dataStruct%ixSnowSoilNrg,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixSnowOnlyNrg); call populate_data_i(dataStruct%ixSnowOnlyNrg,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixSoilOnlyNrg); call populate_data_i(dataStruct%ixSoilOnlyNrg,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixSnowSoilHyd); call populate_data_i(dataStruct%ixSnowSoilHyd,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixSnowOnlyHyd); call populate_data_i(dataStruct%ixSnowOnlyHyd,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixSoilOnlyHyd); call populate_data_i(dataStruct%ixSoilOnlyHyd,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixNrgLayer); call populate_data_i(dataStruct%ixNrgLayer,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixHydLayer); call populate_data_i(dataStruct%ixHydLayer,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
! case(iLookINDEX%ixVolFracWat); call populate_data_i(dataStruct%ixVolFracWat,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
! case(iLookINDEX%ixMatricHead); call populate_data_i(dataStruct%ixMatricHead,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixAllState); call populate_data_i(dataStruct%ixAllState,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixSoilState); call populate_data_i(dataStruct%ixSoilState,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
case(iLookINDEX%ixLayerState); call populate_data_i(dataStruct%ixLayerState,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)
! case(iLookINDEX%ixLayerActive); call populate_data_i(dataStruct%ixLayerActive,tempVec,ix_divide,createLayer,nGRU,stateVariable,ix_lower,ix_upper)

    end select
  ! deallocate the temporary vector: strictly not necessary, but include to be safe
  deallocate(tempVec,stat=err)
  if(err/=0)then; err=20; message='problem deallocating temporary data vector'; return; end if


 end do  ! looping through variables

 end subroutine addModelLayer_indx

 subroutine populate_data(struct, tempVec, ix_divide,createLayer,nGRU,stateVariable,ix_lower, ix_upper)
  real(rkind),device :: struct(ix_lower:,:)
  real(rkind),device :: tempVec(ix_lower:,:)
  integer(i4b),device :: ix_divide(:)
  logical(lgt),device :: createLayer(:)
  integer(i4b) :: ix_lower, ix_upper
  integer(i4b) :: nGRU
  logical(lgt) :: stateVariable

  integer(i4b) :: iGRU, iLayer

  if (stateVariable) then
    if (ix_upper > 0) then
      !$cuf kernel do(2) <<<*,*>>>
      do iGRU=1,nGRU
        do iLayer=1,size(struct,1)-1
          if (createLayer(iGRU) .and. iLayer .le. ix_divide(iGRU)) struct(iLayer,iGRU) = tempVec(iLayer,iGRU)
          if (createLayer(iGRU) .and. iLayer .ge. ix_divide(iGRU)) struct(iLayer+1,iGRU) = tempVec(iLayer,iGRU)
        end do
      end do
    end if
  else
    !$cuf kernel do(2) <<<*,*>>>
    do iGRU=1,nGRU
      do iLayer=1,size(struct,1)
        if (createLayer(iGRU)) struct(iLayer,iGRU) = realMissing
      end do
    end do
  end if
end subroutine

subroutine populate_data_i(struct, tempVec, ix_divide,createLayer,nGRU,stateVariable,ix_lower, ix_upper)
  integer(i4b),device :: struct(ix_lower:,:)
  integer(i4b),device :: tempVec(ix_lower:,:)
  integer(i4b),device :: ix_divide(:)
  logical(lgt),device :: createLayer(:)
  integer(i4b) :: ix_lower, ix_upper
  integer(i4b) :: nGRU
  logical(lgt) :: stateVariable

  integer(i4b) :: iGRU, iLayer

  if (stateVariable) then
    if (ix_upper > 0) then
      !$cuf kernel do(2) <<<*,*>>>
      do iGRU=1,nGRU
        do iLayer=1,size(struct,1)-1
          if (createLayer(iGRU) .and. iLayer .le. ix_divide(iGRU)) struct(iLayer,iGRU) = tempVec(iLayer,iGRU)
          if (createLayer(iGRU) .and. iLayer .ge. ix_divide(iGRU)) struct(iLayer+1,iGRU) = tempVec(iLayer,iGRU)
        end do
      end do
    end if
  else
    !$cuf kernel do(2) <<<*,*>>>
    do iGRU=1,nGRU
      do iLayer=1,size(struct,1)
        if (createLayer(iGRU)) struct(iLayer,iGRU) = realMissing
      end do
    end do
  end if
end subroutine


end module layerDivide_module
