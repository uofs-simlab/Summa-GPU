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

module updateVars_module

! data types
USE nrtype

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! access the global print flag
USE globalData,only:globalPrintFlag

! domain types
USE globalData,only:iname_cas       ! named variables for canopy air space
USE globalData,only:iname_veg       ! named variables for vegetation canopy
USE globalData,only:iname_snow      ! named variables for snow
USE globalData,only:iname_soil      ! named variables for soil
USE globalData,only:iname_aquifer   ! named variables for the aquifer

! named variables to describe the state variable type
USE globalData,only:iname_nrgCanair ! named variable defining the energy of the canopy air space
USE globalData,only:iname_nrgCanopy ! named variable defining the energy of the vegetation canopy
USE globalData,only:iname_watCanopy ! named variable defining the mass of total water on the vegetation canopy
USE globalData,only:iname_liqCanopy ! named variable defining the mass of liquid water on the vegetation canopy
USE globalData,only:iname_nrgLayer  ! named variable defining the energy state variable for snow+soil layers
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the matric head state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid matric potential state variable for soil layers

! metadata for information in the data structures
USE globalData,only:indx_meta       ! metadata for the variables in the index structure

! constants
USE multiconst,only:&
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    zLookup,      & ! data vector with variable length dimension (rkind)
                    var_dlength     ! data vector with variable length dimension (rkind)

! provide access to indices that define elements of the data structures
USE var_lookup,only:iLookDIAG             ! named variables for structure elements
USE var_lookup,only:iLookPROG             ! named variables for structure elements
USE var_lookup,only:iLookDERIV            ! named variables for structure elements
USE var_lookup,only:iLookPARAM            ! named variables for structure elements
USE var_lookup,only:iLookINDEX            ! named variables for structure elements

! provide access to routines to update states
USE updatState_module,only:updateSnow     ! update snow states
USE updatState_module,only:updateSoil     ! update soil states

! provide access to functions for the constitutive functions and derivatives
USE snow_utils_module,only:fracliquid          ! compute the fraction of liquid water (snow)
USE snow_utils_module,only:dFracLiq_dTk        ! differentiate the freezing curve w.r.t. temperature (snow)
USE soil_utils_module,only:dTheta_dTk          ! differentiate the freezing curve w.r.t. temperature (soil)
USE soil_utils_module,only:dTheta_dPsi         ! derivative in the soil water characteristic (soil)
USE soil_utils_module,only:matricHead          ! compute the matric head based on volumetric water content
USE soil_utils_module,only:volFracLiq          ! compute volumetric fraction of liquid water
USE soil_utils_module,only:crit_soilT          ! compute critical temperature below which ice exists
USE soil_utils_module,only:liquidHead          ! compute the liquid water matric potential
USE enthalpyTemp_module,only:T2enthTemp_cas    ! convert temperature to enthalpy for canopy air space
USE enthalpyTemp_module,only:T2enthTemp_veg   ! convert temperature to enthalpy for vegetation
USE enthalpyTemp_module,only:T2enthTemp_snow   ! convert temperature to enthalpy for snow
USE enthalpyTemp_module,only:T2enthTemp_soil   ! convert temperature to enthalpy for soil 

! IEEE check
USE, intrinsic :: ieee_arithmetic         ! check values (NaN, etc.)

implicit none
private

contains


! **********************************************************************************************************
! private subroutine xTempSolve: compute residual and derivative for temperature
! **********************************************************************************************************
subroutine xTempSolve(&
                      ! input: constant over iterations
                      meltNrg          ,&  ! intent(in):    energy for melt+freeze (J m-3)
                      heatCap          ,&  ! intent(in):    volumetric heat capacity (J m-3 K-1)
                      tempInit         ,&  ! intent(in):    initial temperature (K)
                      volFracIceInit   ,&  ! intent(in):    initial volumetric fraction of ice (-)
                      ! input-output: trial values
                      xTemp            ,&  ! intent(inout): trial value of temperature
                      dLiq_dT          ,&  ! intent(in):    derivative in liquid water content w.r.t. temperature (K-1)
                      volFracIceTrial  ,&  ! intent(in):    trial value for volumetric fraction of ice
                      ! output: residual and derivative
                      residual         ,&  ! intent(out):   residual (J m-3)
                      derivative        )  ! intent(out):   derivative (J m-3 K-1)
  implicit none
  ! input: constant over iterations
  real(rkind),intent(in)             :: meltNrg          ! energy for melt+freeze (J m-3)
  real(rkind),intent(in)             :: heatCap          ! volumetric heat capacity (J m-3 K-1)
  real(rkind),intent(in)             :: tempInit         ! initial temperature (K)
  real(rkind),intent(in)             :: volFracIceInit   ! initial volumetric fraction of ice (-)
  ! input-output: trial values
  real(rkind),intent(inout)          :: xTemp            ! trial value for temperature
  real(rkind),intent(in)             :: dLiq_dT          ! derivative in liquid water content w.r.t. temperature (K-1)
  real(rkind),intent(in)             :: volFracIceTrial  ! trial value for the volumetric fraction of ice (-)
  ! output: residual and derivative
  real(rkind),intent(out)            :: residual         ! residual (J m-3)
  real(rkind),intent(out)            :: derivative       ! derivative (J m-3 K-1)
  ! subroutine starts here
  residual   = -heatCap*(xTemp - tempInit) + meltNrg*(volFracIceTrial - volFracIceInit)  ! J m-3
  derivative = heatCap + LH_fus*iden_water*dLiq_dT  ! J m-3 K-1
end subroutine xTempSolve

end module updateVars_module
