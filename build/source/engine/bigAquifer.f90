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
  nGRU, &
                      scalarAquiferStorageTrial,&
                      ! input: diagnostic variables and parameters
                      mpar_data,                    & ! intent(in):    model parameter structure
                      diag_data,                    & ! intent(in):    diagnostic variable structure
                      indx_data,flux_data,deriv_data,&
                      ! input-output: derivatives in transpiration w.r.t. canopy state variables
                      ! output: fluxes and error control
                      out_bigAquifer)                 ! intent(out):   fluxes and error control
  ! named variables
  USE var_lookup,only:iLookDIAG                     ! named variables for structure elements
  USE var_lookup,only:iLookPARAM                    ! named variables for structure elements
  ! data types
  USE data_types,only:var_dlength                   ! x%var(:)%dat [rkind]
  USE data_types,only:out_type_bigAquifer           ! derived typ for intent(out) arguments
  use device_data_types
  use cudafor
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  implicit none
  integer(i4b) :: nGRU
  ! input: state variables, fluxes, and pre-computed derivatives
  real(rkind),device :: scalarAquiferStorageTrial(:)
  ! input: diagnostic variables and parameters
  type(indx_data_device),intent(in) :: indx_data
  type(flux_data_device),intent(inout) :: flux_data
  type(mpar_data_device),intent(in)           :: mpar_data                    ! model parameters
  type(diag_data_device),intent(in)           :: diag_data                    ! diagnostic variables for a local HRU
  type(deriv_data_device),intent(inout) :: deriv_data
  ! input-output: derivatives in transpiration w.r.t. canopy state variables
  ! output: fluxes and error control
  type(out_type_bigAquifer),intent(out)  :: out_bigAquifer               ! fluxes and error control
  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  ! local variables
  real(rkind)                            :: aquiferTranspireFrac         ! fraction of total transpiration that comes from the aquifer (-)
  real(rkind)                            :: xTemp                        ! temporary variable (-)
    type(dim3) :: blocks,threads
    threads = dim3(128,1,1)
  blocks = dim3(nGRU/128+1,1,1)

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! make association between local variables and the information in the data structures
  associate(&
    ! input: state variables, fluxes, and parameters
    scalarCanopyTranspiration => flux_data % ixScalarCanopyTranspiration,   & ! intent(in): [dp] canopy transpiration (kg m-2 s-1)
    scalarSoilDrainage        => flux_data % ixScalarSoilDrainage,          & ! intent(in): [dp] soil drainage (m s-1)
    ! input: pre-computed derivatves
    dCanopyTrans_dCanWat  => deriv_data % dCanopyTrans_dCanWat,    & ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy total water content (s-1)
    dCanopyTrans_dTCanair => deriv_data % dCanopyTrans_dTCanair,   & ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
    dCanopyTrans_dTCanopy => deriv_data % dCanopyTrans_dTCanopy,   & ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
    dCanopyTrans_dTGround => deriv_data % dCanopyTrans_dTGround,   & ! intent(in): [dp] derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
    ! input: model diagnostic variables: contribution of the aquifer to transpiration
    scalarTranspireLim     => diag_data%scalarTranspireLim,     & ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
    scalarAquiferRootFrac  => diag_data%scalarAquiferRootFrac,  & ! intent(in): [dp] fraction of roots below the lowest soil layer (-)
    scalarTranspireLimAqfr => diag_data%scalarTranspireLimAqfr, & ! intent(in): [dp] transpiration limiting factor for the aquifer (-)
    ! input: model parameters: baseflow flux
    aquiferBaseflowRate    => mpar_data%aquiferBaseflowRate_,   & ! intent(in): [dp] tbaseflow rate when aquiferStorage = aquiferScaleFactor (m s-1)
    aquiferScaleFactor     => mpar_data%aquiferScaleFactor_,    & ! intent(in): [dp] scaling factor for aquifer storage in the big bucket (m)
    aquiferBaseflowExp     => mpar_data%aquiferBaseflowExp_,    & ! intent(in): [dp] baseflow exponent (-)
    ! input-output: derivatives in transpiration w.r.t. canopy state variables
    dAquiferTrans_dTCanair => deriv_data % dAquiferTrans_dTCanair, & ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
    dAquiferTrans_dTCanopy => deriv_data % dAquiferTrans_dTCanopy, & ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
    dAquiferTrans_dTGround => deriv_data % dAquiferTrans_dTGround, & ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. ground temperature
    dAquiferTrans_dCanWat  => deriv_data % dAquiferTrans_dCanWat,  & ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy total water
    ! output: fluxes
    scalarAquiferTranspire => flux_data % ixScalarAquiferTranspire,& ! intent(out): transpiration loss from the aquifer (m s-1)
    scalarAquiferRecharge  => flux_data % ixScalarAquiferRecharge, & ! intent(out): recharge to the aquifer (m s-1)
    scalarAquiferBaseflow  => flux_data % ixScalarAquiferBaseflow, & ! intent(out): total baseflow from the aquifer (m s-1)
    dBaseflow_dAquifer     => deriv_data % dBaseflow_dAquifer,    & ! intent(out): change in baseflow flux w.r.t. aquifer storage (s-1)
    ! output: error control
    err                    => out_bigAquifer % err,                   & ! intent(out): error code
    message                => out_bigAquifer % cmessage               & ! intent(out): error message
    )  ! end associating local variables with the information in the data structures
    err=0; message='bigAquifer/' ! initialize error control

    call bigAquifer_kernel<<<blocks,threads>>>(nGRU,flux_data%numfluxdata,indx_data%ixAqWat,flux_data%data,scalarAquiferStorageTrial,scalarCanopyTranspiration,scalarSoilDrainage,&
  dCanopyTrans_dCanWat,dCanopyTrans_dTCanair,dCanopyTrans_dTCanopy,dCanopyTrans_dTGround,&
  scalarTranspireLim,scalarAquiferRootFrac,scalarTranspireLimAqfr,&
  aquiferBaseflowRate,aquiferScaleFactor,aquiferBaseflowExp,&
  dAquiferTrans_dTCanair,dAquiferTrans_dTCanopy,dAquiferTrans_dTGround,dAquiferTrans_dCanWat,&
  scalarAquiferTranspire,scalarAquiferRecharge,scalarAquiferBaseflow,dBaseflow_dAquifer)

  end associate ! end association to data in structure

end subroutine bigAquifer

attributes(global) subroutine bigAquifer_kernel(nGRU,nFluxData,ixAqWat,flux_data,scalarAquiferStorageTrial,scalarCanopyTranspiration,scalarSoilDrainage,&
  dCanopyTrans_dCanWat,dCanopyTrans_dTCanair,dCanopyTrans_dTCanopy,dCanopyTrans_dTGround,&
  scalarTranspireLim,scalarAquiferRootFrac,scalarTranspireLimAqfr,&
  aquiferBaseflowRate,aquiferScaleFactor,aquiferBaseflowExp,&
  dAquiferTrans_dTCanair,dAquiferTrans_dTCanopy,dAquiferTrans_dTGround,dAquiferTrans_dCanWat,&
  scalarAquiferTranspire,scalarAquiferRecharge,scalarAquiferBaseflow,dBaseflow_dAquifer)
  use initialize_device,only:get_iGRU
  integer(i4b),value :: nGRU,nFluxData
  integer(i4b) :: ixAqWat(nGRU)
  real(rkind),intent(inout) :: flux_data(nFluxData,nGRU)
    ! input: state variables, fluxes, and parameters
    real(rkind),intent(in) :: scalarAquiferStorageTrial(nGRU)  ! intent(in): [dp] trial value of aquifer storage (m)
    integer(i4b),intent(in),value :: scalarCanopyTranspiration  ! intent(in): [dp] canopy transpiration (kg m-2 s-1)
    integer(i4b),intent(in),value :: scalarSoilDrainage         ! intent(in): [dp] soil drainage (m s-1)
    ! input: pre-computed derivatves
    real(rkind),intent(in) :: dCanopyTrans_dCanWat(nGRU)   ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy total water content (s-1)
    real(rkind),intent(in) :: dCanopyTrans_dTCanair(nGRU)  ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
    real(rkind),intent(in) :: dCanopyTrans_dTCanopy(nGRU)  ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
    real(rkind),intent(in) :: dCanopyTrans_dTGround(nGRU)  ! intent(in): [dp] derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
    ! input: model diagnostic variables: contribution of the aquifer to transpiration
    real(rkind),intent(in) :: scalarTranspireLim(nGRU)      ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
    real(rkind),intent(in) :: scalarAquiferRootFrac(nGRU)   ! intent(in): [dp] fraction of roots below the lowest soil layer (-)
    real(rkind),intent(in) :: scalarTranspireLimAqfr(nGRU)  ! intent(in): [dp] transpiration limiting factor for the aquifer (-)
    ! input: model parameters: baseflow flux
    real(rkind),intent(in) :: aquiferBaseflowRate(nGRU)     ! intent(in): [dp] tbaseflow rate when aquiferStorage = aquiferScaleFactor (m s-1)
    real(rkind),intent(in) :: aquiferScaleFactor(nGRU)      ! intent(in): [dp] scaling factor for aquifer storage in the big bucket (m)
    real(rkind),intent(in) :: aquiferBaseflowExp(nGRU)      ! intent(in): [dp] baseflow exponent (-)
    ! input-output: derivatives in transpiration w.r.t. canopy state variables
    real(rkind),intent(inout) :: dAquiferTrans_dTCanair(nGRU)  ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
    real(rkind),intent(inout) :: dAquiferTrans_dTCanopy(nGRU)  ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
    real(rkind),intent(inout) :: dAquiferTrans_dTGround(nGRU)  ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. ground temperature
    real(rkind),intent(inout) :: dAquiferTrans_dCanWat(nGRU)   ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy total water
    ! output: fluxes
    integer(i4b),intent(in),value :: scalarAquiferTranspire  ! intent(out): transpiration loss from the aquifer (m s-1)
    integer(i4b),intent(in),value :: scalarAquiferRecharge   ! intent(out): recharge to the aquifer (m s-1)
    integer(i4b),intent(in),value :: scalarAquiferBaseflow   ! intent(out): total baseflow from the aquifer (m s-1)
    real(rkind),intent(inout) :: dBaseflow_dAquifer(nGRU)      ! intent(out): change in baseflow flux w.r.t. aquifer storage (s-1)

  integer(i4b) :: iGRU
  iGRU = get_iGRU()

  if (iGRU .gt. nGRU) return

  if (ixAqWat(iGRU)/=integerMissing) call bigAquifer_device(scalarAquiferStorageTrial,flux_data,scalarCanopyTranspiration,scalarSoilDrainage,&
  dCanopyTrans_dCanWat,dCanopyTrans_dTCanair,dCanopyTrans_dTCanopy,dCanopyTrans_dTGround,&
  scalarTranspireLim,scalarAquiferRootFrac,scalarTranspireLimAqfr,&
  aquiferBaseflowRate,aquiferScaleFactor,aquiferBaseflowExp,&
  dAquiferTrans_dTCanair,dAquiferTrans_dTCanopy,dAquiferTrans_dTGround,dAquiferTrans_dCanWat,&
  scalarAquiferTranspire,scalarAquiferRecharge,scalarAquiferBaseflow,dBaseflow_dAquifer)

end subroutine


attributes(device) subroutine bigAquifer_device(scalarAquiferStorageTrial_,flux_data,scalarCanopyTranspiration,scalarSoilDrainage,&
  dCanopyTrans_dCanWat_,dCanopyTrans_dTCanair_,dCanopyTrans_dTCanopy_,dCanopyTrans_dTGround_,&
  scalarTranspireLim_,scalarAquiferRootFrac_,scalarTranspireLimAqfr_,&
  aquiferBaseflowRate_,aquiferScaleFactor_,aquiferBaseflowExp_,&
  dAquiferTrans_dTCanair_,dAquiferTrans_dTCanopy_,dAquiferTrans_dTGround_,dAquiferTrans_dCanWat_,&
  scalarAquiferTranspire,scalarAquiferRecharge,scalarAquiferBaseflow,dBaseflow_dAquifer_)
  use initialize_device,only:get_iGRU
    ! input: state variables, fluxes, and parameters
    real(rkind),intent(in) :: scalarAquiferStorageTrial_(:)  ! intent(in): [dp] trial value of aquifer storage (m)
    real(rkind),intent(inout) :: flux_data(:,:)
    integer(i4b),intent(in) :: scalarCanopyTranspiration  ! intent(in): [dp] canopy transpiration (kg m-2 s-1)
    integer(i4b),intent(in) :: scalarSoilDrainage         ! intent(in): [dp] soil drainage (m s-1)
    ! input: pre-computed derivatves
    real(rkind),intent(in) :: dCanopyTrans_dCanWat_(:)   ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy total water content (s-1)
    real(rkind),intent(in) :: dCanopyTrans_dTCanair_(:)  ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
    real(rkind),intent(in) :: dCanopyTrans_dTCanopy_(:)  ! intent(in): [dp] derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
    real(rkind),intent(in) :: dCanopyTrans_dTGround_(:)  ! intent(in): [dp] derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
    ! input: model diagnostic variables: contribution of the aquifer to transpiration
    real(rkind),intent(in) :: scalarTranspireLim_(:)      ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
    real(rkind),intent(in) :: scalarAquiferRootFrac_(:)   ! intent(in): [dp] fraction of roots below the lowest soil layer (-)
    real(rkind),intent(in) :: scalarTranspireLimAqfr_(:)  ! intent(in): [dp] transpiration limiting factor for the aquifer (-)
    ! input: model parameters: baseflow flux
    real(rkind),intent(in) :: aquiferBaseflowRate_(:)     ! intent(in): [dp] tbaseflow rate when aquiferStorage = aquiferScaleFactor (m s-1)
    real(rkind),intent(in) :: aquiferScaleFactor_(:)      ! intent(in): [dp] scaling factor for aquifer storage in the big bucket (m)
    real(rkind),intent(in) :: aquiferBaseflowExp_(:)      ! intent(in): [dp] baseflow exponent (-)
    ! input-output: derivatives in transpiration w.r.t. canopy state variables
    real(rkind),intent(inout) :: dAquiferTrans_dTCanair_(:)  ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy air temperature
    real(rkind),intent(inout) :: dAquiferTrans_dTCanopy_(:)  ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy temperature
    real(rkind),intent(inout) :: dAquiferTrans_dTGround_(:)  ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. ground temperature
    real(rkind),intent(inout) :: dAquiferTrans_dCanWat_(:)   ! intent(inout): derivatives in the aquifer transpiration flux w.r.t. canopy total water
    ! output: fluxes
    integer(i4b),intent(in) :: scalarAquiferTranspire  ! intent(out): transpiration loss from the aquifer (m s-1)
    integer(i4b),intent(in) :: scalarAquiferRecharge   ! intent(out): recharge to the aquifer (m s-1)
    integer(i4b),intent(in) :: scalarAquiferBaseflow   ! intent(out): total baseflow from the aquifer (m s-1)
    real(rkind),intent(inout) :: dBaseflow_dAquifer_(:)      ! intent(out): change in baseflow flux w.r.t. aquifer storage (s-1)


  ! local variables
  real(rkind)                            :: aquiferTranspireFrac         ! fraction of total transpiration that comes from the aquifer (-)
  real(rkind)                            :: xTemp                        ! temporary variable (-)
  integer(i4b) :: iGRU
  iGRU = get_iGRU()

    ! compute aquifer transpiration (m s-1)
    aquiferTranspireFrac   = scalarAquiferRootFrac_(iGRU)*scalarTranspireLimAqfr_(iGRU)/scalarTranspireLim_(iGRU)   ! fraction of total transpiration that comes from the aquifer (-)
    flux_data(scalarAquiferTranspire,iGRU) = aquiferTranspireFrac*flux_data(scalarCanopyTranspiration,iGRU)/iden_water         ! aquifer transpiration (kg m-2 s-1 --> m s-1)
    ! derivatives in transpiration w.r.t. canopy state variables
    dAquiferTrans_dCanWat_(iGRU)  = aquiferTranspireFrac*dCanopyTrans_dCanWat_(iGRU) /iden_water
    dAquiferTrans_dTCanair_(iGRU) = aquiferTranspireFrac*dCanopyTrans_dTCanair_(iGRU)/iden_water
    dAquiferTrans_dTCanopy_(iGRU) = aquiferTranspireFrac*dCanopyTrans_dTCanopy_(iGRU)/iden_water
    dAquiferTrans_dTGround_(iGRU) = aquiferTranspireFrac*dCanopyTrans_dTGround_(iGRU)/iden_water

    ! compute aquifer recharge (transfer variables -- included for generality for basin-wide aquifer)
    flux_data(scalarAquiferRecharge,iGRU) = flux_data(scalarSoilDrainage,iGRU) ! m s-1

    ! compute the aquifer baseflow (m s-1)
    xTemp                 = scalarAquiferStorageTrial_(iGRU)/aquiferScaleFactor_(iGRU)
    if (xTemp<0._rkind) xTemp = 0._rkind ! otherwise will give NaN in next line
    flux_data(scalarAquiferBaseflow,iGRU) = aquiferBaseflowRate_(iGRU)*(xTemp**aquiferBaseflowExp_(iGRU))

    ! compute the derivative in the net aquifer flux
    dBaseflow_dAquifer_(iGRU)    = -(aquiferBaseflowExp_(iGRU)*aquiferBaseflowRate_(iGRU)*(xTemp**(aquiferBaseflowExp_(iGRU) - 1._rkind)))/aquiferScaleFactor_(iGRU)
end subroutine

end module bigAquifer_module
