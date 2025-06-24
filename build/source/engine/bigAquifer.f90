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

module bigAquifer_module
  ! -----------------------------------------------------------------------------------------------------------
  ! homegrown solver data types
  USE nrtype
  
  ! access missing values
  USE globalData,only:integerMissing  ! missing integer
  USE globalData,only:realMissing     ! missing real number
  
  ! physical constants
  USE multiconst,only:&
                      LH_vap,  &      ! latent heat of vaporization   (J kg-1)
                      iden_water      ! intrinsic density of water    (kg m-3)
  ! -----------------------------------------------------------------------------------------------------------
  implicit none
  private
  public :: bigAquifer
  contains
  ! ***************************************************************************************************************
  ! public subroutine bigAquifer: compute aquifer water fluxes and their derivatives
  ! ***************************************************************************************************************
  subroutine bigAquifer(&
                        ! input: state variables, fluxes, and pre-computed derivatives
                        ! in_bigAquifer,                & ! intent(in):    state variables, fluxes, and pre-computed derivatives
                        nGRU, &
                        scalarAquiferStorageTrial, &
                        ! input: diagnostic variables and parameters
                        indx_data, &
                        mpar_data,                    & ! intent(in):    model parameter structure
                        diag_data,                    & ! intent(in):    diagnostic variable structure
                        flux_data, &
                        deriv_data, &
                        ! input-output: derivatives in transpiration w.r.t. canopy state variables
                        ! io_bigAquifer,                & ! intent(inout): derivatives in transpiration w.r.t. canopy state variables
                        ! output: fluxes and error control
                        err,message)                 ! intent(out):   fluxes and error control
    ! named variables
    USE var_lookup,only:iLookDIAG,iLookDERIV,iLookFLUX                     ! named variables for structure elements
    USE var_lookup,only:iLookPARAM                    ! named variables for structure elements
    ! data types
    USE data_types,only:var_dlength                   ! x%var(:)%dat [rkind]
    use device_data_types
    ! -------------------------------------------------------------------------------------------------------------------------------------------------
    implicit none
    ! input: state variables, fluxes, and pre-computed derivatives
    ! type(in_type_bigAquifer),intent(in)    :: in_bigAquifer                ! state variables, fluxes, and pre-computed derivatives
    real(rkind),intent(in),device :: scalarAquiferStorageTrial(:)
    integer(i4b),intent(in) :: nGRU
    ! input: diagnostic variables and parameters
    type(indx_data_device) :: indx_data
    type(mpar_data_device),intent(in)           :: mpar_data                    ! model parameters
    type(diag_data_device),intent(in)           :: diag_data                    ! diagnostic variables for a local HRU
    type(flux_data_device),intent(inout) :: flux_data
    type(deriv_data_device),intent(inout) :: deriv_data
    ! input-output: derivatives in transpiration w.r.t. canopy state variables
    ! type(io_type_bigAquifer),intent(inout) :: io_bigAquifer                ! derivatives in transpiration w.r.t. canopy state variables
    ! output: fluxes and error control
    ! type(out_type_bigAquifer),intent(out)  :: out_bigAquifer               ! fluxes and error control
    integer(i4b),intent(out)           :: err                         ! error code
    character(*),intent(out)           :: message                     ! error message
    ! -----------------------------------------------------------------------------------------------------------------------------------------------------
    ! local variables
    real(rkind)                            :: aquiferTranspireFrac         ! fraction of total transpiration that comes from the aquifer (-)
    real(rkind)                            :: xTemp                        ! temporary variable (-)
    integer(i4b) :: iGRU
    ! -------------------------------------------------------------------------------------------------------------------------------------------------
    ! make association between local variables and the information in the data structures
    associate(&
      ! input: state variables, fluxes, and parameters
      ! scalarAquiferStorageTrial => in_bigAquifer % scalarAquiferStorageTrial,   & ! intent(in): [dp] trial value of aquifer storage (m)
      ixAqWat => indx_data%ixAqWat, &
      scalarCanopyTranspiration => flux_data%scalarCanopyTranspiration,   & ! intent(in): [dp] canopy transpiration (kg m-2 s-1)
      scalarSoilDrainage        => flux_data%scalarSoilDrainage,          & ! intent(in): [dp] soil drainage (m s-1)
      ! input: pre-computed derivatves
      dCanopyTrans_dCanWat  => deriv_data%dCanopyTrans_dCanWat,    & ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy total water content (s-1)
      dCanopyTrans_dTCanair => deriv_data%dCanopyTrans_dTCanair,   & ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
      dCanopyTrans_dTCanopy => deriv_data%dCanopyTrans_dTCanopy,   & ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
      dCanopyTrans_dTGround => deriv_data%dCanopyTrans_dTGround,   & ! intent(in): [dp] derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
      ! input: model diagnostic variables: contribution of the aquifer to transpiration
      scalarTranspireLim     => diag_data%scalarTranspireLim,     & ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
      scalarAquiferRootFrac  => diag_data%scalarAquiferRootFrac,  & ! intent(in): [dp] fraction of roots below the lowest soil layer (-)
      scalarTranspireLimAqfr => diag_data%scalarTranspireLimAqfr, & ! intent(in): [dp] transpiration limiting factor for the aquifer (-)
      ! input: model parameters: baseflow flux
      aquiferBaseflowRate    => mpar_data%aquiferBaseflowRate,   & ! intent(in): [dp] tbaseflow rate when aquiferStorage = aquiferScaleFactor (m s-1)
      aquiferScaleFactor     => mpar_data%aquiferScaleFactor,    & ! intent(in): [dp] scaling factor for aquifer storage in the big bucket (m)
      aquiferBaseflowExp     => mpar_data%aquiferBaseflowExp,    & ! intent(in): [dp] baseflow exponent (-)
      ! input-output: derivatives in transpiration w.r.t. canopy state variables
      dAquiferTrans_dTCanair => deriv_data%dAquiferTrans_dTCanair, & ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
      dAquiferTrans_dTCanopy => deriv_data%dAquiferTrans_dTCanopy, & ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
      dAquiferTrans_dTGround => deriv_data%dAquiferTrans_dTGround, & ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. ground temperature
      dAquiferTrans_dCanWat  => deriv_data%dAquiferTrans_dCanWat,  & ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy total water
      ! output: fluxes
      scalarAquiferTranspire => flux_data%scalarAquiferTranspire,& ! intent(out): transpiration loss from the aquifer (m s-1)
      scalarAquiferRecharge  => flux_data%scalarAquiferRecharge, & ! intent(out): recharge to the aquifer (m s-1)
      scalarAquiferBaseflow  => flux_data%scalarAquiferBaseflow, & ! intent(out): total baseflow from the aquifer (m s-1)
      dBaseflow_dAquifer     => deriv_data%dBaseflow_dAquifer    & ! intent(out): change in baseflow flux w.r.t. aquifer storage (s-1)
      ! output: error control
      ! err                    => out_bigAquifer % err,                   & ! intent(out): error code
      ! message                => out_bigAquifer % cmessage               & ! intent(out): error message
      )  ! end associating local variables with the information in the data structures
      err=0; message='bigAquifer/' ! initialize error control
  
      !$cuf kernel do(1) <<<*,*>>>
      do iGRU=1,nGRU
        if (ixAqWat(iGRU)/=integerMissing) then
      ! compute aquifer transpiration (m s-1)
      aquiferTranspireFrac   = scalarAquiferRootFrac(iGRU)*scalarTranspireLimAqfr(iGRU)/scalarTranspireLim(iGRU)   ! fraction of total transpiration that comes from the aquifer (-)
      scalarAquiferTranspire(iGRU) = aquiferTranspireFrac*scalarCanopyTranspiration(iGRU)/iden_water         ! aquifer transpiration (kg m-2 s-1 --> m s-1)
      ! derivatives in transpiration w.r.t. canopy state variables
      dAquiferTrans_dCanWat(iGRU)  = aquiferTranspireFrac*dCanopyTrans_dCanWat(iGRU) /iden_water
      dAquiferTrans_dTCanair(iGRU) = aquiferTranspireFrac*dCanopyTrans_dTCanair(iGRU)/iden_water
      dAquiferTrans_dTCanopy(iGRU) = aquiferTranspireFrac*dCanopyTrans_dTCanopy(iGRU)/iden_water
      dAquiferTrans_dTGround(iGRU) = aquiferTranspireFrac*dCanopyTrans_dTGround(iGRU)/iden_water
  
      ! compute aquifer recharge (transfer variables -- included for generality for basin-wide aquifer)
      scalarAquiferRecharge(iGRU) = scalarSoilDrainage(iGRU) ! m s-1
  
      ! compute the aquifer baseflow (m s-1)
      xTemp                 = scalarAquiferStorageTrial(iGRU)/aquiferScaleFactor
      if (xTemp<0._rkind) xTemp = 0._rkind ! otherwise will give NaN in next line
      scalarAquiferBaseflow(iGRU) = aquiferBaseflowRate*(xTemp**aquiferBaseflowExp)
  
      ! compute the derivative in the net aquifer flux
      dBaseflow_dAquifer(iGRU)    = -(aquiferBaseflowExp*aquiferBaseflowRate*(xTemp**(aquiferBaseflowExp - 1._rkind)))/aquiferScaleFactor
        endif
      end do
  
    end associate ! end association to data in structure
  
  end subroutine bigAquifer
  
  end module bigAquifer_module
  