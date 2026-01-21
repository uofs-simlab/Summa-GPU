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

module paramCheck_module
! define numerical recipes data type
USE nrtype
! define look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers, & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index
implicit none
private
public::paramCheck,paramCheck_device
contains


 ! ************************************************************************************************
 ! public subroutine paramCheck: check consistency of model parameters
 ! ************************************************************************************************
 subroutine paramCheck(mpar_data,err,message)
 ! model decisions
 USE globalData,only:model_decisions  ! model decision structure
 USE var_lookup,only:iLookDECISIONS   ! named variables for elements of the decision structure
 ! SUMMA look-up variables
 USE data_types,only:var_dlength      ! data vector with variable length dimension (rkind): x%var(:)%dat(:)
 USE var_lookup,only:iLookPARAM       ! named variables for elements of the data structures
 implicit none
 ! define input
 type(var_dlength),intent(in)    :: mpar_data            ! model parameters
 ! define output
 integer(i4b),intent(out)        :: err                  ! error code
 character(*),intent(out)        :: message              ! error message
 ! local variables
 integer(i4b)                    :: iLayer               ! index of model layers
 real(rkind),dimension(5)        :: zminLayer            ! minimum layer depth in each layer (m)
 real(rkind),dimension(4)        :: zmaxLayer_lower      ! lower value of maximum layer depth
 real(rkind),dimension(4)        :: zmaxLayer_upper      ! upper value of maximum layer depth
 ! Start procedure here
 err=0; message="paramCheck/"



 end subroutine paramCheck


 attributes(device) subroutine paramCheck_device(critSoilTranspire,critSoilWilting,theta_sat,theta_res,&
 heightCanopyTop, heightCanopyBottom,fieldCapacity,&
 snowLayers,zMax,zMin,&
 zminLayer1,zminLayer2,zminLayer3,zminLayer4,zminLayer5,&
 zmaxLayer1_lower,zmaxLayer2_lower,zmaxLayer3_lower,zmaxLayer4_lower,&
 zmaxLayer1_upper,zmaxLayer2_upper,zmaxLayer3_upper,zmaxLayer4_upper)
  real(rkind) :: critSoilTranspire,critSoilWilting,theta_sat(:),theta_res(:)
  real(rkind) :: heightCanopyTop, heightCanopyBottom,fieldCapacity
   ! local variables
 integer(i4b)                    :: iLayer               ! index of model layers
 real(rkind),dimension(5)        :: zminLayer            ! minimum layer depth in each layer (m)
 real(rkind),dimension(4)        :: zmaxLayer_lower      ! lower value of maximum layer depth
 real(rkind),dimension(4)        :: zmaxLayer_upper      ! upper value of maximum layer depth
 integer(i4b) :: snowLayers
 real(rkind) :: zMax,zMin
 real(rkind) :: zminLayer1,zminLayer2,zminLayer3,zminLayer4,zminLayer5
 real(rkind) :: zmaxLayer1_lower,zmaxLayer2_lower,zmaxLayer3_lower,zmaxLayer4_lower
 real(rkind) :: zmaxLayer1_upper,zmaxLayer2_upper,zmaxLayer3_upper,zmaxLayer4_upper

 ! *****
 ! * check that the snow layer bounds are OK...
 ! ********************************************

 ! select option for combination/sub-division of snow layers
 select case(snowLayers)
  ! SNTHERM option
  case(sameRulesAllLayers)
   if(zmax/zmin < 2.5_rkind)then
    print*, 'zmax must be at least 2.5 times larger than zmin: this avoids merging layers that have just been divided'
    return
   end if
  ! CLM option
  case(rulesDependLayerIndex)
   ! (build vectors of min/max)
   zminLayer       = (/zminLayer1,&
                       zminLayer2,&
                       zminLayer3,&
                       zminLayer4,&
                       zminLayer5/)
   zmaxLayer_lower = (/zmaxLayer1_lower,&
                       zmaxLayer2_lower,&
                       zmaxLayer3_lower,&
                       zmaxLayer4_lower/)
   zmaxLayer_upper = (/zmaxLayer1_upper,&
                       zmaxLayer2_upper,&
                       zmaxLayer3_upper,&
                       zmaxLayer4_upper/)
   ! (check consistency)
   do iLayer=1,4  ! NOTE: the lower layer does not have a maximum value
    ! ensure that we have higher maximum thresholds for sub-division when fewer number of layers
    if(zmaxLayer_lower(iLayer) < zmaxLayer_upper(iLayer))then
     print*, 'expect the maximum threshold for sub-division in the case where there is only ', &
                                  iLayer,' layer(s) is greater than the maximum threshold for sub-division in the case where there are > ',&
                                  iLayer,' layer(s)'
    end if
    ! ensure that the maximum thickness is 3 times greater than the minimum thickness
    if(zmaxLayer_upper(iLayer)/zminLayer(iLayer) < 2.5_rkind .or. zmaxLayer_upper(iLayer)/zminLayer(iLayer+1) < 2.5_rkind)then
     print*, 'zmaxLayer_upper(iLayer), zminLayer(iLayer), zminLayer(iLayer+1) = ', &
                                     zmaxLayer_upper(iLayer), zminLayer(iLayer), zminLayer(iLayer+1)
     print*, 'zmaxLayer_upper for layer ',iLayer,' must be 2.5 times larger than zminLayer for layers ',&
                                  iLayer,' and ',iLayer+1,': this avoids merging layers that have just been divided'
    end if
   end do  ! loop through layers
  case default; print*, 'unable to identify option to combine/sub-divide snow layers'; return
 end select ! (option to combine/sub-divide snow layers)

 ! -------------------------------------------------------------------------------------------------------------------------------------------

 ! *****
 ! * check soil parameter dependencies...
 ! theta_res < critSoilWilting < critSoilTranspire < fieldCapacit < theta_sat
 ! *********************************

   ! check canopy geometry
 if(heightCanopyTop < heightCanopyBottom)then
  print*, 'height of canopy top is less than the height of the canopy bottom'
 endif

 ! check that the maximum transpiration limit is within bounds
 if( any(critSoilTranspire > theta_sat) .or. any(critSoilTranspire < theta_res) )then
!   print*, 'theta_res         = ', theta_res
!   print*, 'theta_sat         = ', theta_sat
  print*, 'critSoilTranspire = ', critSoilTranspire
  print*, 'critSoilTranspire parameter is out of range '// &
                         '[NOTE: if overwriting Noah-MP soil table values in paramTrial, must overwrite all soil parameters]'
 end if

   ! check that the soil wilting point is within bounds
 if( any(critSoilWilting > theta_sat) .or. any(critSoilWilting < theta_res) )then
!   print*, 'theta_res       = ', theta_res
!   print*, 'theta_sat       = ', theta_sat
  print*, 'critSoilWilting = ', critSoilWilting
  print*, 'critSoilWilting parameter is out of range '// &
                         '[NOTE: if overwriting Noah-MP soil table values in paramTrial, must overwrite all soil parameters]'
 end if

 ! check that the field capacity is within bounds
 if( any(fieldCapacity > theta_sat) .or. any(fieldCapacity < theta_res) )then
!   print*, 'theta_res     = ', theta_res
!   print*, 'theta_sat     = ', theta_sat
  print*, 'fieldCapacity = ', fieldCapacity
  print*, 'fieldCapacity parameter is out of range '// &
                         '[NOTE: if overwriting Noah-MP soil table values in paramTrial, must overwrite all soil parameters]'
 end if

  ! check transpiration
  if( critSoilTranspire < critSoilWilting )then
   print*, 'critical point for transpiration is less than the wilting point'
  endif

  ! check porosity
  if( any(theta_sat < theta_res) )then
   print*, 'porosity is less than the residual liquid water content'
  endif

 end subroutine

end module paramCheck_module
