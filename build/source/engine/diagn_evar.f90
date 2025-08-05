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

module diagn_evar_module

! data types
USE nrtype

! derived types to define the data structures
USE data_types,only:&
                    var_d,            & ! data vector (rkind)
                    var_ilength,      & ! data vector with variable length dimension (i4b)
                    var_dlength         ! data vector with variable length dimension (rkind)

! named variables defining elements in the data structures
USE var_lookup,only:iLookPARAM,iLookPROG,iLookDIAG,iLookINDEX  ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure

! physical constants
USE multiconst,only:&
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water,  & ! intrinsic density of water    (kg m-3)
                    ! specific heat
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    Cp_ice,      & ! specific heat of ice          (J kg-1 K-1)
                    Cp_soil,     & ! specific heat of soil         (J kg-1 K-1)
                    Cp_water,    & ! specific heat of liquid water (J kg-1 K-1)
                    ! thermal conductivity
                    lambda_air,  & ! thermal conductivity of air   (J s-1 m-1)
                    lambda_ice,  & ! thermal conductivity of ice   (J s-1 m-1)
                    lambda_water   ! thermal conductivity of water (J s-1 m-1)

! missing values
USE globalData,only:integerMissing ! missing integer
USE globalData,only:realMissing    ! missing real number

! named variables that define the layer type
USE globalData,only:iname_snow     ! snow
USE globalData,only:iname_soil     ! soil

! provide access to named variables for thermal conductivity of soil
USE globalData,only:model_decisions  ! model decision structure

! decisions for thermal conductivity of soil
USE mDecisions_module,only:Smirnova2000    ! option for temporally constant thermal conductivity

! decisions for thermal conductivity of soil
USE mDecisions_module,only: funcSoilWet, & ! function of soil wetness
                            mixConstit,  & ! mixture of constituents
                            hanssonVZJ     ! test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004

! privacy
implicit none
private
public::diagn_evar

! algorithmic parameters
real(rkind),parameter     :: valueMissing=-9999._rkind  ! missing value, used when diagnostic or state variables are undefined
real(rkind),parameter     :: verySmall=1.e-6_rkind   ! used as an additive constant to check if substantial difference among real numbers
real(rkind),parameter     :: mpe=1.e-6_rkind         ! prevents overflow error if division by zero
real(rkind),parameter     :: dx=1.e-6_rkind          ! finite difference increment
contains


 ! **********************************************************************************************************
 ! public subroutine diagn_evar: compute diagnostic energy variables (thermal conductivity and heat capacity)
 ! **********************************************************************************************************
 subroutine diagn_evar(&
                       ! input: control variables
                       computeVegFlux,          & ! intent(in):    flag to denote if computing the vegetation flux
                       nGRU, &
                       decisions, &
                       canopyDepth,             & ! intent(in):    canopy depth (m)
                       ! input/output: data structures
                       mpar_data,               & ! intent(in):    model parameters
                       indx_data,               & ! intent(in):    model layer indices
                       prog_data,               & ! intent(in):    model prognostic variables for a local HRU
                       diag_data,               & ! intent(inout): model diagnostic variables for a local HRU
                        ! output: error control
                       err,message)               ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! provide access to external subroutines
 USE snow_utils_module,only:tcond_snow_d            ! compute thermal conductivity of snow
 use device_data_types
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 logical(lgt),intent(in),device         :: computeVegFlux(:)         ! logical flag to denote if computing the vegetation flux
 integer(i4b) :: nGRU
 real(rkind),intent(in),device          :: canopyDepth(:)            ! depth of the vegetation canopy (m)
 type(decisions_device) :: decisions
 ! input/output: data structures
 type(mpar_data_device),intent(in)    :: mpar_data              ! model parameters
 type(indx_data_device),intent(in)    :: indx_data              ! model layer indices
 type(prog_data_device),intent(in)    :: prog_data              ! model prognostic variables for a local HRU
 type(diag_data_device),intent(inout) :: diag_data              ! model diagnostic variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                    ! error code
 character(*),intent(out)        :: message                ! error message
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)              :: cmessage               ! error message of downwind routine
 integer(i4b)                    :: iLayer                 ! index of model layer
 integer(i4b)                    :: iSoil                  ! index of soil layer
 real(rkind)                     :: TCn                    ! thermal conductivity below the layer interface (W m-1 K-1)
 real(rkind)                     :: TCp                    ! thermal conductivity above the layer interface (W m-1 K-1)
 real(rkind)                     :: zdn                    ! height difference between interface and lower value (m)
 real(rkind)                     :: zdp                    ! height difference between interface and upper value (m)
 real(rkind)                     :: bulkden_soil           ! bulk density of soil (kg m-3)
 real(rkind)                     :: lambda_drysoil         ! thermal conductivity of dry soil (W m-1)
 real(rkind)                     :: lambda_wetsoil         ! thermal conductivity of wet soil (W m-1)
 real(rkind)                     :: lambda_wet             ! thermal conductivity of the wet material
 real(rkind)                     :: relativeSat            ! relative saturation (-)
 real(rkind)                     :: kerstenNum             ! the Kersten number (-), defining weight applied to conductivity of the wet medium
 real(rkind)                     :: den                    ! denominator in the thermal conductivity calculations
 ! local variables to reproduce the thermal conductivity of Hansson et al. VZJ 2005
 real(rkind),parameter           :: c1=0.55_rkind          ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
 real(rkind),parameter           :: c2=0.8_rkind           ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
 real(rkind),parameter           :: c3=3.07_rkind          ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(rkind),parameter           :: c4=0.13_rkind          ! optimized parameter from Hansson et al. VZJ 2005 (W m-1 K-1)
 real(rkind),parameter           :: c5=4._rkind            ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(rkind),parameter           :: f1=13.05_rkind         ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(rkind),parameter           :: f2=1.06_rkind          ! optimized parameter from Hansson et al. VZJ 2005 (-)
 real(rkind)                     :: fArg,xArg              ! temporary variables (see Hansson et al. VZJ 2005 for details)
 integer(i4b) :: iGRU

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in data structure
 associate(&
 ! input: model decisions
 ixThCondSnow            => decisions%thCondSnow,      & ! intent(in): choice of method for thermal conductivity of snow
 ixThCondSoil            => decisions%thCondSoil,      & ! intent(in): choice of method for thermal conductivity of soil
 ! input: state variables
 scalarCanopyIce         => prog_data%scalarCanopyIce,           & ! intent(in): canopy ice content (kg m-2)
 scalarCanopyLiquid      => prog_data%scalarCanopyLiq,           & ! intent(in): canopy liquid water content (kg m-2)
 mLayerVolFracIce        => prog_data%mLayerVolFracIce,             & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
 mLayerVolFracLiq        => prog_data%mLayerVolFracLiq,             & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
 ! input: coordinate variables
 nSnow                   => indx_data%nSnow,                    & ! intent(in): number of snow layers
 nSoil                   => indx_data%nSoil,                    & ! intent(in): number of soil layers
 nLayers                 => indx_data%nLayers_d,                  & ! intent(in): total number of layers
 layerType               => indx_data%layerType,                   & ! intent(in): layer type (iname_soil or iname_snow)
 mLayerHeight            => prog_data%mLayerHeight,                 & ! intent(in): height at the mid-point of each layer (m)
 iLayerHeight            => prog_data%iLayerHeight,                 & ! intent(in): height at the interface of each layer (m)
 ! input: heat capacity and thermal conductivity
 specificHeatVeg         => mpar_data%specificHeatVeg,          & ! intent(in): specific heat of vegetation (J kg-1 K-1)
 maxMassVegetation       => mpar_data%maxMassVegetation,        & ! intent(in): maximum mass of vegetation (kg m-2)
 fixedThermalCond_snow   => mpar_data%fixedThermalCond_snow,    & ! intent(in): temporally constant thermal conductivity of snow (W m-1 K-1)
 ! input: depth varying soil parameters
 iden_soil               => mpar_data%soil_dens_intr,              & ! intent(in): intrinsic density of soil (kg m-3)
 thCond_soil             => mpar_data%thCond_soil,                 & ! intent(in): thermal conductivity of soil (W m-1 K-1)
 theta_sat               => mpar_data%theta_sat,                   & ! intent(in): soil porosity (-)
 frac_sand               => mpar_data%frac_sand,                   & ! intent(in): fraction of sand (-)
 frac_silt               => mpar_data%frac_silt,                   & ! intent(in): fraction of silt (-)
 frac_clay               => mpar_data%frac_clay,                   & ! intent(in): fraction of clay (-)
 ! output: diagnostic variables
 scalarBulkVolHeatCapVeg => diag_data%scalarBulkVolHeatCapVeg,   & ! intent(out): volumetric heat capacity of the vegetation (J m-3 K-1)
 mLayerVolHtCapBulk      => diag_data%mLayerVolHtCapBulk_m,           & ! intent(out): volumetric heat capacity in each layer (J m-3 K-1)
 mLayerThermalC          => diag_data%mLayerThermalC_m,               & ! intent(out): thermal conductivity at the mid-point of each layer (W m-1 K-1)
 iLayerThermalC          => diag_data%iLayerThermalC_m,               & ! intent(out): thermal conductivity at the interface of each layer (W m-1 K-1)
 mLayerVolFracAir        => diag_data%mLayerVolFracAir_m              & ! intent(out): volumetric fraction of air in each layer (-)
 )  ! end associate statement
 ! --------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="diagn_evar/"

 ! initialize the soil layer
 iSoil=integerMissing

 ! compute the bulk volumetric heat capacity of vegetation (J m-3 K-1)
 !$cuf kernel do(1) <<<*,*>>>
 do iGRU=1,nGRU
 if(computeVegFlux(iGRU))then
  scalarBulkVolHeatCapVeg(iGRU) = specificHeatVeg*maxMassVegetation/canopyDepth(iGRU) + & ! vegetation component
                            Cp_water*scalarCanopyLiquid(iGRU)/canopyDepth(iGRU)       + & ! liquid water component
                            Cp_ice*scalarCanopyIce(iGRU)/canopyDepth(iGRU)                ! ice component
 else
  scalarBulkVolHeatCapVeg(iGRU) = valueMissing
 end if
end do


 ! loop through layers
!$cuf kernel do(1) <<<*,*>>>
do iGRU=1,nGRU
 do iLayer=1,nLayers(iGRU)

  ! get the soil layer
  if(iLayer>nSnow(iGRU)) iSoil = iLayer-nSnow(iGRU)

  ! compute the thermal conductivity of dry and wet soils (W m-1)
  ! NOTE: this is actually constant over the simulation, and included here for clarity
  if(ixThCondSoil == funcSoilWet .and. layerType(iLayer,iGRU)==iname_soil)then
   bulkden_soil   = iden_soil(iSoil)*( 1._rkind - theta_sat(iSoil) )
   lambda_drysoil = (0.135_rkind*bulkden_soil + 64.7_rkind) / (iden_soil(iSoil) - 0.947_rkind*bulkden_soil)
   lambda_wetsoil = (8.80_rkind*frac_sand(iSoil) + 2.92_rkind*frac_clay(iSoil)) / (frac_sand(iSoil) + frac_clay(iSoil))
  end if

  ! *****
  ! * compute the volumetric fraction of air in each layer...
  ! *********************************************************
  select case(layerType(iLayer,iGRU))
   case(iname_soil); mLayerVolFracAir(iLayer,iGRU) = theta_sat(iSoil) - (mLayerVolFracIce(iLayer,iGRU) + mLayerVolFracLiq(iLayer,iGRU))
   case(iname_snow); mLayerVolFracAir(iLayer,iGRU) = 1._rkind - (mLayerVolFracIce(iLayer,iGRU) + mLayerVolFracLiq(iLayer,iGRU))
!    case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute volumetric fraction of air'; return
  end select

  ! *****
  ! * compute the volumetric heat capacity of each layer (J m-3 K-1)...
  ! *******************************************************************
  select case(layerType(iLayer,iGRU))
   ! * soil
   case(iname_soil)
    mLayerVolHtCapBulk(iLayer,iGRU) = iden_soil(iSoil)  * Cp_soil  * ( 1._rkind - theta_sat(iSoil) ) + & ! soil component
                                 iden_ice          * Cp_ice   * mLayerVolFracIce(iLayer,iGRU)     + & ! ice component
                                 iden_water        * Cp_water * mLayerVolFracLiq(iLayer,iGRU)     + & ! liquid water component
                                 iden_air          * Cp_air   * mLayerVolFracAir(iLayer,iGRU)         ! air component
   ! * snow
   case(iname_snow)
    mLayerVolHtCapBulk(iLayer,iGRU) = iden_ice          * Cp_ice   * mLayerVolFracIce(iLayer,iGRU)     + & ! ice component
                                 iden_water        * Cp_water * mLayerVolFracLiq(iLayer,iGRU)     + & ! liquid water component
                                 iden_air          * Cp_air   * mLayerVolFracAir(iLayer,iGRU)         ! air component
!    case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute volumetric heat capacity'; return
  end select

  ! *****
  ! * compute the thermal conductivity of snow and soil at the mid-point of each layer...
  ! *************************************************************************************
  select case(layerType(iLayer,iGRU))

   ! ***** soil
   case(iname_soil)

    ! select option for thermal conductivity of soil
    select case(ixThCondSoil)

     ! ** function of soil wetness
     case(funcSoilWet)

      ! compute the thermal conductivity of the wet material (W m-1)
      lambda_wet  = lambda_wetsoil**( 1._rkind - theta_sat(iSoil) ) * lambda_water**theta_sat(iSoil) * lambda_ice**(theta_sat(iSoil) - mLayerVolFracLiq(iLayer,iGRU))
      relativeSat = (mLayerVolFracIce(iLayer,iGRU) + mLayerVolFracLiq(iLayer,iGRU))/theta_sat(iSoil)  ! relative saturation
      ! compute the Kersten number (-)
      if(relativeSat > 0.1_rkind)then ! log10(0.1) = -1
       kerstenNum = log10(relativeSat) + 1._rkind
      else
       kerstenNum = 0._rkind  ! dry thermal conductivity
      endif
      ! ...and, compute the thermal conductivity
      mLayerThermalC(iLayer,iGRU) = kerstenNum*lambda_wet + (1._rkind - kerstenNum)*lambda_drysoil

     ! ** mixture of constituents
     case(mixConstit)
      mLayerThermalC(iLayer,iGRU) = thCond_soil(iSoil) * ( 1._rkind - theta_sat(iSoil) ) + & ! soil component
                               lambda_ice         * mLayerVolFracIce(iLayer,iGRU)     + & ! ice component
                               lambda_water       * mLayerVolFracLiq(iLayer,iGRU)     + & ! liquid water component
                               lambda_air         * mLayerVolFracAir(iLayer,iGRU)         ! air component

     ! ** test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
     case(hanssonVZJ)
      fArg  = 1._rkind + f1*mLayerVolFracIce(iLayer,iGRU)**f2
      xArg  = mLayerVolFracLiq(iLayer,iGRU) + fArg*mLayerVolFracIce(iLayer,iGRU)
      mLayerThermalC(iLayer,iGRU) = c1 + c2*xArg + (c1 - c4)*exp(-(c3*xArg)**c5)

     ! ** check
    !  case default; err=20; message=trim(message)//'unable to identify option for thermal conductivity of soil'; return

    end select  ! option for the thermal conductivity of soil

   ! ***** snow
   case(iname_snow)
    ! temporally constant thermal conductivity
    if(ixThCondSnow==Smirnova2000)then
     mLayerThermalC(iLayer,iGRU) = fixedThermalCond_snow
    ! thermal conductivity as a function of snow density
    else
     call tcond_snow_d(mLayerVolFracIce(iLayer,iGRU)*iden_ice,  & ! input: snow density (kg m-3)
                     mLayerThermalC(iLayer,iGRU),ixThCondSnow)!,             & ! output: thermal conductivity (W m-1 K-1)
                    !  err,cmessage)                         ! output: error control
    !  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
    endif

   ! * error check
!    case default; err=20; message=trim(message)//'unable to identify type of layer (snow or soil) to compute thermal conductivity'; return

  end select
  !print*, 'iLayer, mLayerThermalC(iLayer) = ', iLayer, mLayerThermalC(iLayer)

 end do  ! looping through layers
end do


 ! *****
 ! * compute the thermal conductivity of snow at the interface of each layer...
 ! ****************************************************************************
!$cuf kernel do(1) <<<*,*>>>
do iGRU=1,nGRU
 do iLayer=1,nLayers(iGRU)-1  ! (loop through layers)
  ! get temporary variables
  TCn = mLayerThermalC(iLayer,iGRU)    ! thermal conductivity below the layer interface (W m-1 K-1)
  TCp = mLayerThermalC(iLayer+1,iGRU)  ! thermal conductivity above the layer interface (W m-1 K-1)
  zdn = iLayerHeight(iLayer,iGRU)   - mLayerHeight(iLayer,iGRU) ! height difference between interface and lower value (m)
  zdp = mLayerHeight(iLayer+1,iGRU) - iLayerHeight(iLayer,iGRU) ! height difference between interface and upper value (m)
  den = TCn*zdp + TCp*zdn  ! denominator
  ! compute thermal conductivity
  if(TCn+TCp > epsilon(TCn))then
   iLayerThermalC(iLayer,iGRU) = (TCn*TCp*(zdn + zdp)) / den
  else
   iLayerThermalC(iLayer,iGRU) = (TCn*zdn +  TCp*zdp) / (zdn + zdp)
  endif
  !write(*,'(a,1x,i4,1x,10(f9.3,1x))') 'iLayer, TCn, TCp, zdn, zdp, iLayerThermalC(iLayer) = ', iLayer, TCn, TCp, zdn, zdp, iLayerThermalC(iLayer)
 end do  ! looping through layers
end do

!$cuf kernel do(1) <<<*,*>>>
do iGRU=1,nGRU
 ! special case of hansson
 if(ixThCondSoil==hanssonVZJ)then
  iLayerThermalC(0,iGRU) = 28._rkind*(0.5_rkind*(iLayerHeight(1,iGRU) - iLayerHeight(0,iGRU)))
 else
  iLayerThermalC(0,iGRU) = mLayerThermalC(1,iGRU)
 end if

 ! assume the thermal conductivity at the domain boundaries is equal to the thermal conductivity of the layer
 iLayerThermalC(nLayers(iGRU),iGRU) = mLayerThermalC(nLayers(iGRU),iGRU)
end do
 ! end association to variables in the data structure
 end associate

 end subroutine diagn_evar


end module diagn_evar_module
