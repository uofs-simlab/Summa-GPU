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

module read_icond_module
USE nrtype
USE netcdf
USE globalData,only: ixHRUfile_min,ixHRUfile_max
USE globalData,only: nTimeDelay   ! number of hours in the time delay histogram
USE globalData,only: nSpecBand    ! number of spectral bands
USE globalData,only:verySmaller   ! a smaller number used as an additive constant to check if substantial difference among real numbers

implicit none
private
public::read_icond
public::read_icond_nlayers
! define single HRU restart file
integer(i4b), parameter :: singleHRU=1001
integer(i4b), parameter :: multiHRU=1002
integer(i4b), parameter :: restartFileType=multiHRU
contains

 ! ************************************************************************************************
 ! public subroutine read_icond_nlayers: read model initial conditions file for number of snow/soil layers
 ! ************************************************************************************************
 subroutine read_icond_nlayers(iconFile,nGRU,indx_meta,err,message)
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE var_lookup,only:iLookINDEX                        ! variable lookup structure
 USE globalData,only:gru_struc                         ! gru-hru mapping structures
 USE globalData,only:startGRU                          ! index of first gru for parallel runs
 USE netcdf_util_module,only:nc_file_close             ! close netcdf file
 USE netcdf_util_module,only:nc_file_open              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                ! netcdf error handling
 USE data_types,only:gru_hru_intVec                    ! actual data
 USE data_types,only:var_info                          ! metadata
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)        ,intent(in)     :: iconFile       ! name of input (restart) file
 integer(i4b)        ,intent(in)     :: nGRU           ! total # of GRUs in run domain
 type(var_info)      ,intent(in)     :: indx_meta(:)   ! metadata
 integer(i4b)        ,intent(out)    :: err            ! error code
 character(*)        ,intent(out)    :: message        ! returned error message
 ! locals
 integer(i4b)             :: ncID                      ! netcdf file id
 integer(i4b)             :: dimID                     ! netcdf file dimension id
 integer(i4b)             :: fileHRU                   ! number of HRUs in netcdf file
 integer(i4b)             :: snowID, soilID            ! netcdf variable ids
 integer(i4b)             :: iGRU, iHRU                ! loop indexes
 integer(i4b)             :: iHRU_global               ! index of HRU in the netcdf file
 integer(i4b),allocatable :: snowData(:)               ! number of snow layers in all HRUs
 integer(i4b),allocatable :: soilData(:)               ! number of soil layers in all HRUs
 character(len=256)       :: cmessage                  ! downstream error message

 ! --------------------------------------------------------------------------------------------------------
 ! initialize error message
 err=0
 message = 'read_icond_nlayers/'

 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncid,err,cmessage);
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! get number of HRUs in file (the GRU variable(s), if present, are processed at the end)
 err = nf90_inq_dimid(ncID,"hru",dimId);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID,dimId,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

 ! allocate storage for reading from file (allocate entire file size, even when doing subdomain run)
 allocate(snowData(fileHRU))
 allocate(soilData(fileHRU))
 snowData = 0
 soilData = 0

 ! get netcdf ids for the variables holding number of snow and soil layers in each hru
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookINDEX%nSnow)%varName),snowid); call netcdf_err(err,message)
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookINDEX%nSoil)%varName),soilid); call netcdf_err(err,message)

 ! get nSnow and nSoil data (reads entire state file)
 err = nf90_get_var(ncid,snowid,snowData); call netcdf_err(err,message)
 err = nf90_get_var(ncid,soilid,soilData); call netcdf_err(err,message)

 ixHRUfile_min=huge(1)
 ixHRUfile_max=0
 ! find the min and max hru indices in the state file
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   if(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc < ixHRUfile_min) ixHRUfile_min = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
   if(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc > ixHRUfile_max) ixHRUfile_max = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
  end do
 end do

 ! loop over grus in current run to update snow/soil layer information
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   iHRU_global = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
   ! single HRU (Note: 'restartFileType' is hardwired above to multiHRU)
   if(restartFileType==singleHRU) then
    gru_struc(iGRU)%hruInfo(iHRU)%nSnow = snowData(1)
    gru_struc(iGRU)%hruInfo(iHRU)%nSoil = soilData(1)
   ! multi HRU
   else
    gru_struc(iGRU)%hruInfo(iHRU)%nSnow = snowData(iHRU_global)
    gru_struc(iGRU)%hruInfo(iHRU)%nSoil = soilData(iHRU_global)
   endif
  end do
 end do

 ! close file
 call nc_file_close(ncid,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(snowData,soilData)

 end subroutine read_icond_nlayers


 ! ************************************************************************************************
 ! public subroutine read_icond: read model initial conditions
 ! ************************************************************************************************
 subroutine read_icond(iconFile,                      & ! intent(in):    name of initial conditions file
                       nGRU,                          & ! intent(in):    number of GRUs
                       mparData,                      & ! intent(in):    model parameters
                       progData,                      & ! intent(inout): model prognostic variables
                       bvarData,                      & ! intent(inout): model basin (GRU) variables
                       indxData,                      & ! intent(inout): model indices
                       no_icond_enth,                 & ! intent(out):   flag that enthalpy variables are not in the file
                       err,message)                     ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE var_lookup,only:iLookVarType                       ! variable lookup structure
 USE var_lookup,only:iLookPROG                          ! variable lookup structure
 USE var_lookup,only:iLookPARAM                         ! variable lookup structure
 USE var_lookup,only:iLookBVAR                          ! variable lookup structure
 USE var_lookup,only:iLookINDEX                         ! variable lookup structure
 USE globalData,only:prog_meta                          ! metadata for prognostic variables
 USE globalData,only:bvar_meta                          ! metadata for basin (GRU) variables
 USE globalData,only:gru_struc                          ! gru-hru mapping structures
 USE globalData,only:startGRU                          ! index of first gru for parallel runs
 USE globalData,only:iname_soil,iname_snow              ! named variables to describe the type of layer
 USE netcdf_util_module,only:nc_file_open               ! open netcdf file
 USE netcdf_util_module,only:nc_file_close              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                 ! netcdf error handling
 USE data_types,only:gru_hru_doubleVec                  ! full double precision structure
 USE data_types,only:gru_hru_intVec                     ! full integer structure
 USE data_types,only:gru_doubleVec                      ! gru-length double precision structure (basin variables)
 USE data_types,only:var_dlength                        ! double precision structure for a single HRU
 USE data_types,only:var_info                           ! metadata
 USE get_ixName_module,only:get_varTypeName             ! to access type strings for error messages
 USE updatState_module,only:updateSoil                  ! update soil states
 use device_data_types

 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)           ,intent(in)     :: iconFile                 ! name of netcdf file containing the initial conditions
 integer(i4b)           ,intent(in)     :: nGRU                     ! number of grouped response units in simulation domain
 type(mpar_data_device),intent(in)     :: mparData                 ! model parameters
 type(prog_data_device),intent(inout)  :: progData                 ! model prognostic variables
 type(bvar_data_device)    ,intent(inout)  :: bvarData                 ! model basin (GRU) variables
 type(indx_data_device)   ,intent(inout)  :: indxData                 ! model indices
 logical                ,intent(out)    :: no_icond_enth            ! flag that enthalpy variables are not in the file
 integer(i4b)           ,intent(out)    :: err                      ! error code
 character(*)           ,intent(out)    :: message                  ! returned error message
 ! locals
 character(len=256)                     :: cmessage                 ! downstream error message
 integer(i4b)                           :: fileHRU                  ! number of HRUs in file
 integer(i4b)                           :: fileGRU                  ! number of GRUs in file
 integer(i4b)                           :: iVar, i                  ! loop indices
 integer(i4b),dimension(1)              :: ndx                      ! intermediate array of loop indices
 integer(i4b)                           :: iGRU                     ! loop index
 integer(i4b)                           :: iHRU                     ! loop index
 integer(i4b)                           :: dimID                    ! varible dimension ids
 integer(i4b)                           :: ncVarID                  ! variable ID in netcdf file
 character(256)                         :: dimName                  ! not used except as a placeholder in call to inq_dim function
 integer(i4b)                           :: dimLen                   ! data dimensions
 integer(i4b)                           :: ncID                     ! netcdf file ID
 integer(i4b)                           :: ixFile                   ! index in file
 integer(i4b)                           :: iHRU_local               ! index of HRU in the data subset
 integer(i4b)                           :: iHRU_global              ! index of HRU in the netcdf file
 real(rkind),allocatable                :: varData(:,:)             ! variable data storage
 integer(i4b)                           :: nSoil, nSnow, nToto      ! # layers
 integer(i4b)                           :: nTDH                     ! number of points in time-delay histogram
 integer(i4b)                           :: iLayer,jLayer            ! layer indices
 character(len=32),parameter            :: scalDimName   ='scalarv' ! dimension name for scalar data
 character(len=32),parameter            :: midSoilDimName='midSoil' ! dimension name for soil-only layers
 character(len=32),parameter            :: midTotoDimName='midToto' ! dimension name for layered varaiables
 character(len=32),parameter            :: ifcTotoDimName='ifcToto' ! dimension name for layered varaiables
 character(len=32),parameter            :: tdhDimName    ='tdh'     ! dimension name for time-delay basin variables
     type(dim3) :: blocks,threads
     threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

 ! --------------------------------------------------------------------------------------------------------
 ! Start procedure here
 err=0; message="read_icond/"

 ! --------------------------------------------------------------------------------------------------------
 ! (1) read the file
 ! --------------------------------------------------------------------------------------------------------
 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncID,err,cmessage)
 if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

 ! get number of HRUs in file
 err = nf90_inq_dimid(ncID,"hru",dimID);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncID,dimID,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

 ! loop through prognostic variables
 no_icond_enth=.false.
 do iVar = 1,size(prog_meta)

  ! skip variables that are computed later
  if(prog_meta(iVar)%varName=='scalarCanopyWat'           .or. &
     prog_meta(iVar)%varName=='spectralSnowAlbedoDiffuse' .or. &
     prog_meta(iVar)%varName=='scalarSurfaceTemp'         .or. &
     prog_meta(iVar)%varName=='mLayerVolFracWat'          .or. &
     prog_meta(iVar)%varName=='mLayerHeight'                   ) cycle

  ! get variable id
  err = nf90_inq_varid(ncID,trim(prog_meta(iVar)%varName),ncVarID)
  if(err/=nf90_noerr)then
   if(prog_meta(iVar)%varName=='scalarCanairEnthalpy'     .or. &
      prog_meta(iVar)%varName=='scalarCanopyEnthalpy'     .or. &  
      prog_meta(iVar)%varName=='mLayerEnthalpy'                )then; err=nf90_noerr; no_icond_enth=.true.; cycle; endif ! skip enthalpy variables if not in file
   call netcdf_err(err,message)
   message=trim(message)//': problem with getting variable id, var='//trim(prog_meta(iVar)%varName)
   return
  endif

  ! get variable dimension IDs
  select case (prog_meta(iVar)%varType)
   case (iLookVarType%scalarv); err = nf90_inq_dimid(ncID,trim(scalDimName)   ,dimID); call netcdf_err(err,message)
   case (iLookVarType%midSoil); err = nf90_inq_dimid(ncID,trim(midSoilDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%midToto); err = nf90_inq_dimid(ncID,trim(midTotoDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%ifcToto); err = nf90_inq_dimid(ncID,trim(ifcTotoDimName),dimID); call netcdf_err(err,message)
   case default
    message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
    err=20; return
  end select

  ! check errors
  if(err/=0)then
   message=trim(message)//': problem with dimension ids, var='//trim(prog_meta(iVar)%varName)
   return
  endif

  ! get the dimension length
  err = nf90_inquire_dimension(ncID,dimID,dimName,dimLen); call netcdf_err(err,message)
  if(err/=0)then; message=trim(message)//': problem getting the dimension length'; return; endif

  ! initialize the variable data
  allocate(varData(fileHRU,dimLen),stat=err)
  if(err/=0)then; message=trim(message)//'problem allocating HRU variable data'; return; endif

  ! get data
  err = nf90_get_var(ncID,ncVarID,varData); call netcdf_err(err,message)
  if(err/=0)then; message=trim(message)//': problem getting the data for variable '//trim(prog_meta(iVar)%varName); return; endif

  ! store data in prognostics structure
  ! loop through GRUs
  do iGRU = 1,nGRU
   do iHRU = 1,gru_struc(iGRU)%hruCount

    iHRU_global = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
    iHRU_local = (iHRU_global - ixHRUfile_min) + 1

    ! get the number of layers
    nSnow = gru_struc(iGRU)%hruInfo(iHRU)%nSnow
    nSoil = gru_struc(iGRU)%hruInfo(iHRU)%nSoil
    nToto = nSnow + nSoil

    ! get the index in the file: single HRU
    if(restartFileType==singleHRU)then
     ixFile = 1  ! use for single HRU restart file
    ! get the index in the file: multi HRU
    else
     ixFile = iHRU_global
    endif

    ! put the data into data structures and check that none of the values are set to nf90_fill_double
    select case (prog_meta(iVar)%varType)
     case (iLookVarType%scalarv)
      select case(iVar)
     case(iLookPROG%dt_init); progData%dt_init = varData(ixFile,1)
     case(iLookPROG%scalarCanopyIce); progData%scalarCanopyIce(iGRU) = varData(ixFile,1)
     case(iLookPROG%scalarCanopyLiq); progData%scalarCanopyLiq(iGRU) = varData(ixFile,1)
     case(iLookPROG%scalarCanopyWat); progData%scalarCanopyWat(iGRU) = varData(ixFile,1)
     case(iLookPROG%scalarCanairTemp); progData%scalarCanairTemp(iGRU) = varData(ixFile,1)
     case(iLookPROG%scalarCanopyTemp); progData%scalarCanopyTemp(iGRU) = varData(ixFile,1)
     case(iLookPROG%scalarSnowAlbedo); progData%scalarSnowAlbedo(iGRU) = varData(ixFile,1)
     case(iLookPROG%scalarSnowDepth); progData%scalarSnowDepth(iGRU) = varData(ixFile,1)
     case(iLookPROG%scalarSWE); progData%scalarSWE(iGRU) = varData(ixFile,1)
     case(iLookPROG%scalarSfcMeltPond); progData%scalarSfcMeltPond(iGRU) = varData(ixFile,1)
     case(iLookPROG%scalarCanairEnthalpy); progData%scalarCanairEnthalpy(iGRU) = varData(ixFile,1)
     case(iLookPROG%scalarCanopyEnthalpy); progData%scalarCanopyEnthalpy(iGRU) = varData(ixFile,1)
     case(iLookPROG%scalarAquiferStorage); progData%scalarAquiferStorage(iGRU) = varData(ixFile,1)
     case(iLookPROG%scalarSurfaceTemp); progData%scalarSurfaceTemp(iGRU) = varData(ixFile,1)
     end select
     print*, prog_meta(iVar)%varName, iGRU, varData(ixFile,1)
     case (iLookVarType%midSoil)
      select case (iVar)
     case(iLookPROG%spectralSnowAlbedoDiffuse); progData%spectralSnowAlbedoDiffuse(1:nSoil,iGRU) = varData(ixFile,1:nSoil)
     case(iLookPROG%mLayerTemp); progData%mLayerTemp(1:nSoil,iGRU) = varData(ixFile,1:nSoil)
     case(iLookPROG%mLayerVolFracIce); progData%mLayerVolFracIce(1:nSoil,iGRU) = varData(ixFile,1:nSoil)
     case(iLookPROG%mLayerVolFracLiq); progData%mLayerVolFracLiq(1:nSoil,iGRU) = varData(ixFile,1:nSoil)
     case(iLookPROG%mLayerVolFracWat); progData%mLayerVolFracWat(1:nSoil,iGRU) = varData(ixFile,1:nSoil)
     case(iLookPROG%mLayerMatricHead); progData%mLayerMatricHead(1:nSoil,iGRU) = varData(ixFile,1:nSoil)
     case(iLookPROG%mLayerEnthalpy); progData%mLayerEnthalpy(1:nSoil,iGRU) = varData(ixFile,1:nSoil)
     case(iLookPROG%mLayerDepth); progData%mLayerDepth(1:nSoil,iGRU) = varData(ixFile,1:nSoil)
     case(iLookPROG%mLayerHeight); progData%mLayerHeight(1:nSoil,iGRU) = varData(ixFile,1:nSoil)
     case(iLookPROG%iLayerHeight); progData%iLayerHeight(1:nSoil,iGRU) = varData(ixFile,1:nSoil)
     end select
     case (iLookVarType%midToto)
      select case (iVar)
     case(iLookPROG%spectralSnowAlbedoDiffuse); progData%spectralSnowAlbedoDiffuse(1:nToto,iGRU) = varData(ixFile,1:nToto)
     case(iLookPROG%mLayerTemp); progData%mLayerTemp(1:nToto,iGRU) = varData(ixFile,1:nToto)
     case(iLookPROG%mLayerVolFracIce); progData%mLayerVolFracIce(1:nToto,iGRU) = varData(ixFile,1:nToto)
     case(iLookPROG%mLayerVolFracLiq); progData%mLayerVolFracLiq(1:nToto,iGRU) = varData(ixFile,1:nToto)
     case(iLookPROG%mLayerVolFracWat); progData%mLayerVolFracWat(1:nToto,iGRU) = varData(ixFile,1:nToto)
     case(iLookPROG%mLayerMatricHead); progData%mLayerMatricHead(1:nToto,iGRU) = varData(ixFile,1:nToto)
     case(iLookPROG%mLayerEnthalpy); progData%mLayerEnthalpy(1:nToto,iGRU) = varData(ixFile,1:nToto)
     case(iLookPROG%mLayerDepth); progData%mLayerDepth(1:nToto,iGRU) = varData(ixFile,1:nToto)
     case(iLookPROG%mLayerHeight); progData%mLayerHeight(1:nToto,iGRU) = varData(ixFile,1:nToto)
     case(iLookPROG%iLayerHeight); progData%iLayerHeight(1:nToto,iGRU) = varData(ixFile,1:nToto)
     end select
     case (iLookVarType%ifcToto)
      select case (iVar)
     case(iLookPROG%spectralSnowAlbedoDiffuse); progData%spectralSnowAlbedoDiffuse(0:nToto,iGRU) = varData(ixFile,1:nToto+1)
     case(iLookPROG%mLayerTemp); progData%mLayerTemp(0:nToto,iGRU) = varData(ixFile,1:nToto+1)
     case(iLookPROG%mLayerVolFracIce); progData%mLayerVolFracIce(0:nToto,iGRU) = varData(ixFile,1:nToto+1)
     case(iLookPROG%mLayerVolFracLiq); progData%mLayerVolFracLiq(0:nToto,iGRU) = varData(ixFile,1:nToto+1)
     case(iLookPROG%mLayerVolFracWat); progData%mLayerVolFracWat(0:nToto,iGRU) = varData(ixFile,1:nToto+1)
     case(iLookPROG%mLayerMatricHead); progData%mLayerMatricHead(0:nToto,iGRU) = varData(ixFile,1:nToto+1)
     case(iLookPROG%mLayerEnthalpy); progData%mLayerEnthalpy(0:nToto,iGRU) = varData(ixFile,1:nToto+1)
     case(iLookPROG%mLayerDepth); progData%mLayerDepth(0:nToto,iGRU) = varData(ixFile,1:nToto+1)
     case(iLookPROG%mLayerHeight); progData%mLayerHeight(0:nToto,iGRU) = varData(ixFile,1:nToto+1)
     case(iLookPROG%iLayerHeight); progData%iLayerHeight(0:nToto,iGRU) = varData(ixFile,1:nToto+1)
     end select
     case default
      print*, "unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
      err=20; return
    end select

    if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(prog_meta(iVar)%varName)//"')"; return; endif

   end do ! iHRU
  end do ! iGRU

  ! deallocate storage vector for next variable
  deallocate(varData, stat=err)
  if(err/=0)then; message=trim(message)//'problem deallocating HRU variable data'; return; endif

 end do ! end looping through prognostic variables (iVar)

 ! --------------------------------------------------------------------------------------------------------
 ! (2) set number of layers
 ! --------------------------------------------------------------------------------------------------------
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount

   ! save the number of layers
   indxData%nSnow(iGRU)   = gru_struc(iGRU)%hruInfo(iHRU)%nSnow
   indxData%nSoil   = gru_struc(iGRU)%hruInfo(iHRU)%nSoil
   indxData%nLayers_d(iGRU) = gru_struc(iGRU)%hruInfo(iHRU)%nSnow + gru_struc(iGRU)%hruInfo(iHRU)%nSoil

  end do
 end do
 call read_icond_corrections<<<blocks,threads>>>(nGRU,&
 progData%scalarSnowAlbedo,mparData%albedoMax_,progData%spectralSnowAlbedoDiffuse,&
 indxData%nSnow,indxData%nSoil,indxData%layerType,&
 progData%mLayerTemp,progData%mLayerMatricHead,&
 mparData%vGn_alpha_,mparData%vGn_n_,&
 mparData%theta_sat_,mparData%theta_res_,&
 progData%mLayerVolFracWat,progData%mLayerVolFracLiq,progData%mLayerVolFracIce)

 ! --------------------------------------------------------------------------------------------------------
 ! (2) now get the basin variable(s)
 ! --------------------------------------------------------------------------------------------------------

 ! get the index in the file: single HRU
 if(restartFileType/=singleHRU)then

  ! get dimension of time delay histogram (TDH) from initial conditions file
  err = nf90_inq_dimid(ncID,"tdh",dimID);
  if(err/=nf90_noerr)then
   write(*,*) 'WARNING: routingRunoffFuture is not in the initial conditions file ... using zeros'  ! previously created in var_derive.f90
   err=nf90_noerr    ! reset this err

  else
   ! the state file *does* have the basin variable(s), so process them
   err = nf90_inquire_dimension(ncID,dimID,len=nTDH);
   if(err/=nf90_noerr)then; message=trim(message)//'problem reading tdh dimension from initial condition file/'//trim(nf90_strerror(err)); return; end if

   ! get number of GRUs in file
   err = nf90_inq_dimid(ncID,"gru",dimID);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding gru dimension/'//trim(nf90_strerror(err)); return; end if
   err = nf90_inquire_dimension(ncID,dimID,len=fileGRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading gru dimension/'//trim(nf90_strerror(err)); return; end if

   ! check vs hardwired value set in globalData.f90
   if(nTDH /= nTimeDelay)then
    write(*,*) 'tdh=',nTDH,' nTimeDelay=',nTimeDelay
    message=trim(message)//': state file time delay dimension tdh does not match summa expectation of nTimeDelay set in globalData()'
    return
   endif

   ! loop through specific basin variables (currently 1 but loop provided to enable inclusion of others)
   ndx = (/iLookBVAR%routingRunoffFuture/)   ! array of desired variable indices
   do i = 1,size(ndx)
    iVar = ndx(i)

    ! get tdh dimension Id in file (should be 'tdh')
    err = nf90_inq_dimid(ncID,trim(tdhDimName), dimID);
    if(err/=0)then; message=trim(message)//': problem with dimension ids for tdh vars'; return; endif

    ! get the tdh dimension length (dimName and dimLen are outputs of this call)
    err = nf90_inquire_dimension(ncID,dimID,dimName,dimLen); call netcdf_err(err,message)
    if(err/=0)then; message=trim(message)//': problem getting the dimension length for tdh vars'; return; endif

    ! get tdh-based variable id
    err = nf90_inq_varid(ncID,trim(bvar_meta(iVar)%varName),ncVarID); call netcdf_err(err,message)
    if(err/=0)then; message=trim(message)//': problem with getting basin variable id, var='//trim(bvar_meta(iVar)%varName); return; endif

    ! initialize the tdh variable data
    allocate(varData(fileGRU,dimLen),stat=err)
    if(err/=0)then; print*, 'err= ',err; message=trim(message)//'problem allocating GRU variable data'; return; endif

    ! get data
    err = nf90_get_var(ncID,ncVarID,varData); call netcdf_err(err,message)
    if(err/=0)then; message=trim(message)//': problem getting the data'; return; endif

    ! store data in basin var (bvar) structure
    do iGRU = 1,nGRU

     ! put the data into data structures
      select case(iVar)
      case(iLookBVAR%routingRunoffFuture); bvarData%routingRunoffFuture(:,iGRU) = varData((iGRU+startGRU-1),1:nTDH)
      case(iLookBVAR%routingFractionFuture); bvarData%routingFractionFuture(:,iGRU) = varData((iGRU+startGRU-1),1:nTDH)
      end select
     ! check whether the first values is set to nf90_fill_double

    end do ! end iGRU loop

    ! deallocate temporary data array for next variable
    deallocate(varData, stat=err)
    if(err/=0)then; message=trim(message)//'problem deallocating GRU variable data'; return; endif

   end do ! end looping through basin variables
  endif  ! end if case for tdh variables being in init. cond. file
 endif  ! end if case for not being a singleHRU run
 
 call nc_file_close(ncID,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 end subroutine read_icond
 attributes(global) subroutine read_icond_corrections(nGRU,&
 scalarSnowAlbedo_,albedoMax_,spectralSnowAlbedoDiffuse_,&
 nSnow_,nSoil,layerType_,&
 mLayerTemp_,mLayerMatricHead_,&
 vGn_alpha_,vGn_n_,&
 theta_sat_,theta_res_,&
 mLayerVolFracWat_,mLayerVolFracLiq_,mLayerVolFracIce_)
  USE globalData,only:iname_soil,iname_snow              ! named variables to describe the type of layer
  USE updatState_module,only:updateSoil                  ! update soil states

  integer(i4b),intent(in),value :: nGRU
  real(rkind) :: scalarSnowAlbedo_(:)
  real(rkind) :: albedoMax_(:), spectralSnowAlbedoDiffuse_(:,:)
  integer(i4b),intent(in),value :: nSoil
  integer(i4b),intent(in) :: nSnow_(:)
  integer(i4b),intent(inout) :: layerType_(:,:)
  real(rkind) :: mLayerTemp_(:,:),mLayerMatricHead_(:,:)
  real(rkind) :: vGn_alpha_(:,:), vGn_n_(:,:)
  real(rkind) :: theta_sat_(:,:),theta_res_(:,:)
  real(rkind) :: mLayerVolFracWat_(:,:),mLayerVolFracLiq_(:,:),mLayerVolFracIce_(:,:)


   integer(i4b) :: iGRU, iLayer, jLayer
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return

     ! make sure snow albedo is not negative
     if(scalarSnowAlbedo_(iGRU) < 0._rkind)then
      scalarSnowAlbedo_(iGRU) = albedoMax_(iGRU)
     endif

     ! initialize the spectral albedo
     spectralSnowAlbedoDiffuse_(:,iGRU) = scalarSnowAlbedo_(iGRU)

  ! --------------------------------------------------------------------------------------------------------
 ! (2) set number of layers
 ! --------------------------------------------------------------------------------------------------------
   ! set layer type
   layerType_(1:nSnow_(iGRU),iGRU) = iname_snow
   layerType_(nSnow_(iGRU)+1:(nSnow_(iGRU)+nSoil),iGRU) = iname_soil


  ! --------------------------------------------------------------------------------------------------------
 ! (3) update soil layers (diagnostic variables)
 ! --------------------------------------------------------------------------------------------------------
   ! loop through soil layers
   do iLayer = 1,nSoil

    ! get layer in the total vector
    jLayer = iLayer+nSnow_(iGRU)

    ! update soil layers
    call updateSoil(&
                    ! input
                    mLayerTemp_(jLayer,iGRU),& ! intent(in): temperature vector (K)
                    mLayerMatricHead_(iLayer,iGRU),& ! intent(in): matric head (m)
                    vGn_alpha_(iLayer,iGRU),& ! intent(in): van Genutchen "alpha" parameter
                    vGn_n_(iLayer,iGRU),& ! intent(in): van Genutchen "n" parameter
                    theta_sat_(iLayer,iGRU),& ! intent(in): soil porosity (-)
                    theta_res_(iLayer,iGRU),& ! intent(in): soil residual volumetric water content (-)
                    1._rkind - 1._rkind/vGn_n_(iLayer,iGRU),& ! intent(in): van Genutchen "m" parameter (-)
                    ! output
                    mLayerVolFracWat_(jLayer,iGRU),& ! intent(out): volumetric fraction of total water (-)
                    mLayerVolFracLiq_(jLayer,iGRU),& ! intent(out): volumetric fraction of liquid water (-)
                    mLayerVolFracIce_(jLayer,iGRU)& ! intent(out): volumetric fraction of ice (-)
                    )                                                                   ! intent(out): error control

   end do  ! looping through soil layers

 end subroutine

end module read_icond_module
