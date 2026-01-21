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

module check_icond_module
USE nrtype

! access missing values
USE globalData,only:integerMissing   ! missing integer
USE globalData,only:realMissing      ! missing real number

implicit none
private
public::check_icond
contains

 ! ************************************************************************************************
 ! public subroutine check_icond: read model initial conditions
 ! ************************************************************************************************
attributes(device) subroutine check_icond(checkEnthalpy,no_icond_enth,use_lookup,mLayerDepth,iLayerHeight,nLayers,nSnow,&
  scalarSWE_,mLayerVolFracLiq,mLayerVolFracIce,layerType,&
  vGn_n, vGn_alpha, theta_res, theta_sat,&
  mLayerTemp,mLayerEnthTemp,mLayerEnthalpy,&
  temperature,psiLiq_int,deriv2,&
  mLayerMatricHead,soil_dens_intr,snowfrz_scale,&
  scalarCanairTemp_,scalarCanairEnthalpy_,&
  heightCanopyTop_,heightCanopyBottom_,&
  specificHeatVeg_,maxMassVegetation_,&
  scalarCanopyTemp_,scalarCanopyEnthTemp_,scalarCanopyEnthalpy_,scalarCanopyIce_,scalarCanopyLiq_,&
  scalarSnowAlbedo_,albedoMax_,albedoMinWinter_,&
  spectralSnowAlbedoDiffuse,albedoMaxVisible_,&
  albedoMinVisible_,albedoMaxNearIR_,albedoMinNearIR_)
 USE nrtype
 USE var_lookup,only:iLookPARAM                          ! variable lookup structure
 USE var_lookup,only:iLookPROG                           ! variable lookup structure
 USE var_lookup,only:iLookDIAG                           ! variable lookup structure
 USE var_lookup,only:iLookINDEX                          ! variable lookup structure
 USE globalData,only:gru_struc                           ! gru-hru mapping structures
 USE data_types,only:gru_hru_doubleVec                   ! actual data
 USE data_types,only:gru_hru_intVec                      ! actual data
 USE data_types,only:gru_hru_z_vLookup                   ! actual data
 USE globalData,only:iname_soil,iname_snow               ! named variables to describe the type of layer
 USE multiconst,only:&
                       LH_fus,    &                      ! latent heat of fusion                (J kg-1)
                       iden_ice,  &                      ! intrinsic density of ice             (kg m-3)
                       iden_water,&                      ! intrinsic density of liquid water    (kg m-3)
                       gravity,   &                      ! gravitational acceleration           (m s-2)
                       Tfreeze                           ! freezing point of pure water         (K)
 USE snow_utils_module,only:fracliquid                   ! compute volumetric fraction of liquid water in snow based on temperature
 USE updatState_module,only:updateSnow                   ! update snow states
 USE updatState_module,only:updateSoil                   ! update soil states
 USE enthalpyTemp_module,only:T2enthTemp_cas             ! convert temperature to enthalpy for canopy air space
 USE enthalpyTemp_module,only:T2enthTemp_veg             ! convert temperature to enthalpy for vegetation
 USE enthalpyTemp_module,only:T2enthTemp_snow            ! convert temperature to enthalpy for snow
 USE enthalpyTemp_module,only:T2enthTemp_soil            ! convert temperature to enthalpy for soil
 
 implicit none

  integer(i4b),intent(in) :: nLayers,nSnow
  real(rkind) :: mLayerDepth(:), iLayerHeight(0:)
  real(rkind) :: scalarSWE_(:)
  real(rkind) :: mLayerVolFracLiq(:), mLayerVolFracIce(:)
  integer(i4b) :: layerType(:)
  real(rkind),intent(in) :: vGn_n(:)
  real(rkind),intent(in) :: vGn_alpha(:), theta_res(:), theta_sat(:)
  real(rkind),intent(in) :: mLayerTemp(:)
  real(rkind),intent(inout) :: mLayerEnthTemp(:), mLayerEnthalpy(:)
 logical(lgt),intent(in)               :: checkEnthalpy         ! if true either need enthTemp as starting residual value, or for state variable initialization
  logical(lgt),intent(in)               :: no_icond_enth         ! if true, no enthalpy in icond file
 logical(lgt),intent(in)               :: use_lookup            ! flag to use the lookup table for soil enthalpy, otherwise use hypergeometric function
 real(rkind),intent(in) :: snowfrz_scale

 real(rkind),intent(in) :: temperature(:,:), psiLiq_int(:,:), deriv2(:,:)
 real(rkind) :: mLayerMatricHead(:), soil_dens_intr(:)
 real(rkind) :: scalarCanairTemp_(:), scalarCanairEnthalpy_(:)
 real(rkind),intent(in) :: heightCanopyTop_(:),heightCanopyBottom_(:)
 real(rkind),intent(in) :: specificHeatVeg_(:),maxMassVegetation_(:)
  real(rkind) :: scalarCanopyTemp_(:),scalarCanopyEnthTemp_(:),scalarCanopyEnthalpy_(:),scalarCanopyIce_(:),scalarCanopyLiq_(:)
  real(rkind) :: scalarSnowAlbedo_(:), albedoMax_(:), albedoMinWinter_(:)
  real(rkind) :: spectralSnowAlbedoDiffuse(:)
  real(rkind) :: albedoMaxVisible_(:), albedoMinVisible_(:), albedoMaxNearIR_(:), albedoMinNearIR_(:)

integer(i4b)                          :: iGRU,iHRU           ! loop index
 ! temporary variables for realism checks
 integer(i4b)                          :: iLayer                ! index of model layer
 integer(i4b)                          :: iSoil                 ! index of soil layer
 real(rkind)                           :: fLiq                  ! fraction of liquid water on the vegetation canopy (-)
 real(rkind)                           :: vGn_m                 ! van Genutchen "m" parameter (-)
 real(rkind)                           :: scalarTheta           ! liquid water equivalent of total water [liquid water + ice] (-)
 real(rkind)                           :: h1,h2                 ! used to check depth and height are consistent
 real(rkind)                           :: kappa                 ! constant in the freezing curve function (m K-1)
 real(rkind),parameter                 :: xTol=1.e-10_rkind     ! small tolerance to address precision issues
 real(rkind),parameter                 :: canIceTol=1.e-3_rkind ! small tolerance to allow existence of canopy ice for above-freezing temperatures (kg m-2)
 
  integer(i4b) :: err
    iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
     ! --------------------------------------------------------------------------------------------------------
 ! Check that the initial conditions do not conflict with parameters, structure, etc.
 ! --------------------------------------------------------------------------------------------------------

   ! ensure the spectral average albedo is realistic
   if(scalarSnowAlbedo_(iGRU) > albedoMax_(iGRU)) &
      scalarSnowAlbedo_(iGRU) = albedoMax_(iGRU)
   if(scalarSnowAlbedo_(iGRU) < albedoMinWinter_(iGRU)) &
      scalarSnowAlbedo_(iGRU) = albedoMinWinter_(iGRU)
   ! ensure the visible albedo is realistic
   if(spectralSnowAlbedoDiffuse(1) > albedoMaxVisible_(iGRU)) &
      spectralSnowAlbedoDiffuse(1) = albedoMaxVisible_(iGRU)
   if(spectralSnowAlbedoDiffuse(1) < albedoMinVisible_(iGRU)) &
      spectralSnowAlbedoDiffuse(1) = albedoMinVisible_(iGRU)
   ! ensure the nearIR albedo is realistic
   if(spectralSnowAlbedoDiffuse(2) > albedoMaxNearIR_(iGRU)) &
      spectralSnowAlbedoDiffuse(2) = albedoMaxNearIR_(iGRU)
   if(spectralSnowAlbedoDiffuse(2) < albedoMinNearIR_(iGRU)) &
      spectralSnowAlbedoDiffuse(2) = albedoMinNearIR_(iGRU)

    ! ensure the initial conditions are consistent with the constitutive functions
       ! compute the constant in the freezing curve function (m K-1)
   kappa  = (iden_ice/iden_water)*(LH_fus/(gravity*Tfreeze))  ! NOTE: J = kg m2 s-2

   ! check canopy ice content for unrealistic situations
   if(scalarCanopyIce_(iGRU) > canIceTol .and. scalarCanopyTemp_(iGRU) > Tfreeze)then
    ! ice content > threshold, terminate run
    print*, 'canopy ice (=',scalarCanopyIce_(iGRU),') > canIceTol (=',canIceTol,') when canopy temperature (=',scalarCanopyTemp_(iGRU),') > Tfreeze (=',Tfreeze,')'
    err=20; return
   else if(scalarCanopyIce_(iGRU) > 0._rkind .and. scalarCanopyTemp_(iGRU) > Tfreeze)then
    ! if here, ice content < threshold. Could be sublimation on previous timestep or simply wrong input. Print a warning
    print*, 'Warning: canopy ice content in restart file (=',scalarCanopyIce_(iGRU),') > 0 when canopy temperature (=',scalarCanopyTemp_(iGRU),') > Tfreeze (=',Tfreeze,'). Continuing.',NEW_LINE('a')
   end if


   scalarTheta = scalarCanopyIce_(iGRU) + scalarCanopyLiq_(iGRU)

   if(checkEnthalpy)then ! enthalpy as state variable or in residual
     if(no_icond_enth)then ! no enthalpy in icond file
       call T2enthTemp_cas(&
                  scalarCanairTemp_(iGRU),       & ! intent(in): canopy air temperature (K)
                  scalarCanairEnthalpy_(iGRU))     ! intent(out): enthalpy of the canopy air space (J m-3)
 
       call T2enthTemp_veg(&
                  (heightCanopyTop_(iGRU)-heightCanopyBottom_(iGRU)), & ! intent(in): canopy depth (m)
                  specificHeatVeg_(iGRU),        & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                  maxMassVegetation_(iGRU),      & ! intent(in): maximum mass of vegetation (kg m-2)
                  snowfrz_scale,          & ! intent(in): scaling parameter for the snow freezing curve  (K-1)
                  scalarCanopyTemp_(iGRU),       & ! intent(in): canopy temperature (K)
                  scalarTheta,            & ! intent(in): canopy water content (kg m-2)
                  scalarCanopyEnthTemp_(iGRU))     ! intent(out): temperature component of enthalpy of the vegetation canopy (J m-3)
       scalarCanopyEnthalpy_(iGRU) = scalarCanopyEnthTemp_(iGRU) - LH_fus * scalarCanopyIce_(iGRU)/ (heightCanopyTop_(iGRU)-heightCanopyBottom_(iGRU))
     else ! enthalpy is in the icond file
       scalarCanopyEnthTemp_(iGRU) = scalarCanopyEnthalpy_(iGRU) + LH_fus * scalarCanopyIce_(iGRU)/ (heightCanopyTop_(iGRU)-heightCanopyBottom_(iGRU))
     end if
   end if

       ! loop through all layers
   do iLayer=1,nLayers

    ! *****
    ! * check that the initial volumetric fraction of liquid water and ice is reasonable...
    ! *************************************************************************************
    select case(layerType(iLayer))

     ! ***** snow, volume expansion allowed
     case(iname_snow)
      scalarTheta = mLayerVolFracIce(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiq(iLayer)
      ! (check liquid water)
      if(mLayerVolFracLiq(iLayer) < 0._rkind)then;
       print*, 'cannot initialize the model with volumetric fraction of liquid water < 0: layer = ',iLayer; 
       err=20; 
       return; 
      end if
      if(mLayerVolFracLiq(iLayer) > 1._rkind)then; 
      print*, 'cannot initialize the model with volumetric fraction of liquid water > 1: layer = ',iLayer; 
      err=20; 
      return; 
      end if
      ! (check ice)
      if(mLayerVolFracIce(iLayer) > 0.80_rkind)then;
       print*, 'cannot initialize the model with volumetric fraction of ice > 0.80: layer = ',iLayer; 
       err=20; 
       return;
       end if
      if(mLayerVolFracIce(iLayer) < 0.05_rkind)then; 
       print*, 'cannot initialize the model with volumetric fraction of ice < 0.05: layer = ',iLayer; 
       err=20; 
       return; 
      end if
      ! check total water
      if(scalarTheta > 0.80_rkind)then; 
        print*, 'cannot initialize the model with total water fraction [liquid + ice] > 0.80: layer = ',iLayer; 
        err=20; 
        return; 
      end if
      if(scalarTheta < 0.05_rkind)then; 
        print*, 'cannot initialize the model with total water fraction [liquid + ice] < 0.05: layer = ',iLayer; 
        err=20; 
        return; 
      end if

     ! ***** soil, no volume expansion
     case(iname_soil)
      iSoil       = iLayer - nSnow
      if(vGn_n(iSoil) <= 1._rkind)then; 
      print*, 'cannot have van Genutchen n <= 1: soil layer = ',iSoil; 
      err=20; 
      return; 
      end if
      if(vGn_alpha(iSoil) >= 0._rkind)then; 
      print*, 'cannot have van Genutchen alpha >= 0: soil layer = ',iSoil; 
      err=20; 
      return; 
      end if
      vGn_m       = 1._rkind - 1._rkind/vGn_n(iSoil)
      scalarTheta = mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)
      ! (check liquid water)
      if(mLayerVolFracLiq(iLayer) < theta_res(iSoil)-xTol)then; print*, 'cannot initialize the model with volumetric fraction of liquid water < theta_res: layer = ',iLayer; err=20; return; end if
      if(mLayerVolFracLiq(iLayer) > theta_sat(iSoil)+xTol)then; print*, 'cannot initialize the model with volumetric fraction of liquid water > theta_sat: layer = ',iLayer; err=20; return; end if
      ! (check ice)
      if(mLayerVolFracIce(iLayer) < 0._rkind             )then; print*, 'cannot initialize the model with volumetric fraction of ice < 0: layer = '        ,iLayer; err=20; return; end if
      if(mLayerVolFracIce(iLayer) > theta_sat(iSoil)+xTol)then; print*, 'cannot initialize the model with volumetric fraction of ice > theta_sat: layer = ',iLayer; err=20; return; end if
      ! check total water
      if(scalarTheta < theta_res(iSoil)-xTol)then; print*, 'cannot initialize the model with total water fraction [liquid + ice] < theta_res: layer = ',iLayer; err=20; return; end if
      if(scalarTheta > theta_sat(iSoil)+xTol)then; print*, 'cannot initialize the model with total water fraction [liquid + ice] > theta_sat: layer = ',iLayer; err=20; return; end if

     case default
      print*, 'Cannot recognize case in initial vol water/ice check: type=', layerType(iLayer)
    end select

    ! *****
    ! * check that the initial conditions are consistent with the constitutive functions...
    ! *************************************************************************************
    select case(layerType(iLayer))

     ! ** snow
     case(iname_snow)

      ! check that snow temperature is less than freezing
      if(mLayerTemp(iLayer) > Tfreeze)then
       print*, 'initial snow temperature is greater than freezing'
      end if

      ! ensure consistency among state variables
      call updateSnow(&
                      mLayerTemp(iLayer),             & ! intent(in): temperature (K)
                      scalarTheta,                    & ! intent(in): volumetric fraction of total water (-)
                      snowfrz_scale,                  & ! intent(in): scaling parameter for the snow freezing curve (K-1)
                      mLayerVolFracLiq(iLayer),       & ! intent(out): volumetric fraction of liquid water (-)
                      mLayerVolFracIce(iLayer),       & ! intent(out): volumetric fraction of ice (-)
                      fLiq                           & ! intent(out): fraction of liquid water (-)
                      )                     ! intent(out): error control

      if(checkEnthalpy)then ! enthalpy as state variable or in residual
        if(no_icond_enth)then ! no enthalpy in icond file
          call T2enthTemp_snow(&
                      snowfrz_scale,                  & ! intent(in):  scaling parameter for the snow freezing curve  (K-1)
                      mLayerTemp(iLayer),             & ! intent(in):  layer temperature (K)
                      scalarTheta,                    & ! intent(in):  volumetric total water content (-)
                      mLayerEnthTemp(iLayer))           ! intent(out): temperature component of enthalpy of each snow layer (J m-3)
          mLayerEnthalpy(iLayer) = mLayerEnthTemp(iLayer) - iden_ice * LH_fus * mLayerVolFracIce(iLayer)
        else
          mLayerEnthTemp(iLayer) = mLayerEnthalpy(iLayer) + iden_ice * LH_fus * mLayerVolFracIce(iLayer)
        end if
      endif

     ! ** soil
     case(iname_soil)

      ! ensure consistency among state variables
      call updateSoil(&
                      mLayerTemp(iLayer),              & ! intent(in): layer temperature (K)
                      mLayerMatricHead(iLayer-nSnow),  & ! intent(in): matric head (m)
                      vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m, & ! intent(in): van Genutchen soil parameters
                      scalarTheta,                     & ! intent(out): volumetric fraction of total water (-)
                      mLayerVolFracLiq(iLayer),        & ! intent(out): volumetric fraction of liquid water (-)
                      mLayerVolFracIce(iLayer)        & ! intent(out): volumetric fraction of ice (-)
                      )                      ! intent(out): error control

      if(checkEnthalpy)then ! enthalpy as state variable or in residual
        if(no_icond_enth)then ! no enthalpy in icond file
          call T2enthTemp_soil(&
                      use_lookup,                      & ! intent(in):  flag to use the lookup table for soil enthalpy
                      soil_dens_intr(iSoil),           & ! intent(in):  intrinsic soil density (kg m-3)
                      vGn_alpha(iSoil),vGn_n(iSoil),theta_sat(iSoil),theta_res(iSoil),vGn_m, & ! intent(in): van Genutchen soil parameters
                      iSoil,                           & ! intent(in):  index of the control volume within the domain
                      temperature(:,iSoil),psiLiq_int(:,iSoil),deriv2(:,iSoil),  & ! intent(in):  lookup table data structure
                      realMissing,                     & ! intent(in):  lower value of integral (not computed)
                      mLayerTemp(iLayer),              & ! intent(in):  layer temperature (K)
                      mLayerMatricHead(iLayer-nSnow),  & ! intent(in):  matric head (m)
                      mLayerEnthTemp(iLayer))            ! intent(out): temperature component of enthalpy soil layer (J m-3)
          mLayerEnthalpy(iLayer) = mLayerEnthTemp(iLayer) - iden_water * LH_fus * mLayerVolFracIce(iLayer)
        else
          mLayerEnthTemp(iLayer) = mLayerEnthalpy(iLayer) + iden_water * LH_fus * mLayerVolFracIce(iLayer)
        end if
      endif

    end select

   end do  ! (looping through layers)


     ! if snow layers exist, compute snow depth and SWE
   if(nSnow > 0)then
    scalarSWE_(iGRU) = sum( (mLayerVolFracLiq(1:nSnow)*iden_water + &
                                                                         mLayerVolFracIce(1:nSnow)*iden_ice)  * &
                                                                         mLayerDepth(1:nSnow) )
   end if  ! if snow layers exist

    ! check that the layering is consistent
   do iLayer=1,nLayers
    h1 = sum(mLayerDepth(1:iLayer)) ! sum of the depths up to the current layer
    h2 = iLayerHeight(iLayer) - iLayerHeight(0)  ! difference between snow-atm interface and bottom of layer
    if(abs(h1 - h2) > 1.e-6_rkind)then
     print*, 'mis-match between layer depth and layer height; layer = ', iLayer, '; sum depths = ',h1,'; height = ',h2
    !  return
    end if
   end do

end subroutine

end module check_icond_module
