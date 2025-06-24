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

module snowAlbedo_module

! data types
USE nrtype                          ! numerical recipes data types

! physical constants
USE multiconst,only:Tfreeze         ! freezing point of pure water (K)

! derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    model_options   ! defines the model decisions

! named variables defining elements in the data structures
USE var_lookup,only:iLookPARAM,iLookFLUX,iLookDIAG,iLookPROG   ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS                             ! named variables for elements of the decision structure

! look-up values for the choice of snow albedo options
USE mDecisions_module,only:  &
 constantDecay,              &      ! constant decay in snow albedo (e.g., VIC, CLASS)
 variableDecay                      ! variable decay in snow albedo (e.g., BATS approach, with destructive metamorphism + soot content)

! look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:  &
 noah_mp,                    &      ! full Noah-MP implementation (including albedo)
 CLM_2stream,                &      ! CLM 2-stream model (see CLM documentation)
 UEB_2stream,                &      ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
 NL_scatter,                 &      ! Simplified method Nijssen and Lettenmaier (JGR 1999)
 BeersLaw                           ! Beer's Law (as implemented in VIC)

! -------------------------------------------------------------------------------------------------
! privacy
implicit none
private
public::snowAlbedo

! dimensions
integer(i4b),parameter        :: nBands=2      ! number of spectral bands for shortwave radiation
contains


 ! *******************************************************************************************************
 ! public subroutine snowAlbedo: muster program to compute energy fluxes at vegetation and ground surfaces
 ! *******************************************************************************************************
 subroutine snowAlbedo(&
                       ! input: model control
                       dt,                                    & ! intent(in):    model time step (s)
                       nSnow,                          & ! intent(in):    logical flag to denote if snow is present
                       nGRU, &
                       ! input/output: data structures
                       model_decisions,                       & ! intent(in):    model decisions
                       decisions, &
                       mpar_data,                             & ! intent(in):    model parameters
                       flux_data,                             & ! intent(in):    model flux variables
                       diag_data,                             & ! intent(inout): model diagnostic variables for a local HRU
                       prog_data,                             & ! intent(inout): model prognostic variables for a local HRU
                       ! output: error control
                       err,message)                             ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! provide access to desired modules
 USE snow_utils_module,only:fracliquid_d                          ! compute fraction of liquid water at a given temperature
 use device_data_types
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 real(rkind),intent(in)          :: dt                            ! model time step
 integer(i4b),intent(in),device         :: nSnow(:)                  ! logical flag to denote if snow is present
 integer(i4b),intent(in) :: nGRU
 ! input/output: data structures
 type(model_options),intent(in)  :: model_decisions(:)            ! model decisions
 type(decisions_device) :: decisions
 type(mpar_data_device),intent(in)    :: mpar_data                     ! model parameters
 type(flux_data_device),intent(in)    :: flux_data                     ! model flux variables
 type(diag_data_device),intent(inout) :: diag_data                     ! model diagnostic variables for a local HRU
 type(prog_data_device),intent(inout) :: prog_data                     ! model prognostic variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                           ! error code
 character(*),intent(out)        :: message                       ! error message
 ! local variables
 integer(i4b),parameter          :: ixVisible=1                   ! named variable to define index in array of visible part of the spectrum
 integer(i4b),parameter          :: ixNearIR=2                    ! named variable to define index in array of near IR part of the spectrum
 real(rkind),parameter           :: valueMissing=-9999._rkind     ! missing value -- will cause problems if snow albedo is ever used for the non-snow case
 real(rkind),parameter           :: slushExp=10._rkind            ! "slush" exponent, to increase decay when snow is near Tfreeze
 real(rkind),parameter           :: fractionLiqThresh=0.001_rkind ! threshold for the fraction of liquid water to switch to spring albedo minimum
 real(rkind)                     :: fractionLiq                   ! fraction of liquid water (-)
 real(rkind)                     :: age1,age2,age3                ! aging factors (-)
 real(rkind),device                     :: decayFactor(nGRU)                   ! albedo decay factor (-)
 real(rkind),device                     :: refreshFactor(nGRU)                 ! albedo refreshment factor, representing albedo increase due to snowfall (-)
 real(rkind),device                     :: albedoMin(nGRU)                     ! minimum albedo -- depends if in winter or spring conditions (-)
 real(rkind)                     :: fZen                          ! factor to modify albedo at low zenith angles (-)
 real(rkind),parameter           :: bPar=2._rkind                 ! empirical parameter in fZen
 integer(i4b) :: iGRU
 ! initialize error control
 err=0; message='snowAlbedo/'
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&
 ! input: model decisions
 ixCanopySrad              => model_decisions(iLookDECISIONS%canopySrad)%iDecision,   & ! intent(in): index of method used for canopy sw radiation
 ixAlbedoMethod            => decisions%alb_method,   & ! intent(in): index of method used for snow albedo
 ! input: model parameters
 Frad_vis                  => mpar_data%Frad_vis,              & ! intent(in): fraction of radiation in visible part of spectrum (-)
 Frad_direct               => mpar_data%Frad_direct,           & ! intent(in): fraction direct solar radiation (-)
 albedoMax                 => mpar_data%albedoMax,             & ! intent(in): maximum snow albedo for a single spectral band (-)
 albedoMinWinter           => mpar_data%albedoMinWinter,       & ! intent(in): minimum snow albedo during winter for a single spectral band (-)
 albedoMinSpring           => mpar_data%albedoMinSpring,       & ! intent(in): minimum snow albedo during spring for a single spectral band (-)
 albedoMaxVisible          => mpar_data%albedoMaxVisible,      & ! intent(in): maximum snow albedo in the visible part of the spectrum (-)
 albedoMinVisible          => mpar_data%albedoMinVisible,      & ! intent(in): minimum snow albedo in the visible part of the spectrum (-)
 albedoMaxNearIR           => mpar_data%albedoMaxNearIR,       & ! intent(in): maximum snow albedo in the near infra-red part of the spectrum (-)
 albedoMinNearIR           => mpar_data%albedoMinNearIR,       & ! intent(in): minimum snow albedo in the near infra-red part of the spectrum (-)
 albedoDecayRate           => mpar_data%albedoDecayRate,       & ! intent(in): albedo decay rate (s)
 tempScalGrowth            => mpar_data%tempScalGrowth,        & ! intent(in): temperature scaling factor for grain growth (K-1)
 albedoSootLoad            => mpar_data%albedoSootLoad,        & ! intent(in): soot load factor (-)
 albedoRefresh             => mpar_data%albedoRefresh,         & ! intent(in): critical mass necessary for albedo refreshment (kg m-2)
 snowfrz_scale             => mpar_data%snowfrz_scale,         & ! intent(in): scaling parameter for the freezing curve for snow (K-1)
 ! input: model variables
 surfaceTemp               => prog_data%mLayerTemp,             & ! intent(in): surface temperature
 snowfallRate              => flux_data%scalarSnowfall,         & ! intent(in): snowfall rate (kg m-2 s-1)
 cosZenith                 => diag_data%scalarCosZenith,        & ! intent(in): cosine of the zenith angle (-)
 ! input/output: model variables
 scalarSnowAlbedo          => prog_data%scalarSnowAlbedo,       & ! intent(inout): snow albedo for the entire spectral band (-)
 spectralSnowAlbedoDirect  => diag_data%spectralSnowAlbedoDirect,  & ! intent(inout): direct snow albedo in each spectral band (-)
 spectralSnowAlbedoDiffuse => prog_data%spectralSnowAlbedoDiffuse  & ! intent(inout): diffuse snow albedo in each spectral band (-)
 ) ! end associate statement
 ! --------------------------------------------------------------------------------------------------------------------------------------

 ! return early if computing radiation in noah-MP
 if(ixCanopySrad==noah_mp) return

 !$cuf kernel do(1) <<<*,*>>>
 do iGRU=1,nGRU
 ! return early if no snow
 if(nSnow(iGRU) == 0)then
  scalarSnowAlbedo(iGRU)             = valueMissing
  spectralSnowAlbedoDirect(:,iGRU)  = valueMissing
  spectralSnowAlbedoDiffuse(:,iGRU) = valueMissing
  cycle
 end if

 ! compute fractional increase in albedo associated with snowfall
 refreshFactor(iGRU) = dt*snowfallRate(iGRU)/albedoRefresh

 ! identify option for snow albedo
 select case(ixAlbedoMethod)


  ! *** constant decay rate
  case(constantDecay)
   ! compute decay rate
   decayFactor(iGRU) = dt/albedoDecayRate
   ! compute minimum albedo
   fractionLiq = fracliquid_d(surfaceTemp(1,iGRU),snowfrz_scale) ! fraction of liquid water
   if(scalarSnowAlbedo(iGRU) < albedoMinWinter .or. fractionLiq > fractionLiqThresh)then
    albedoMin(iGRU) = albedoMinSpring
   else
    albedoMin(iGRU) = albedoMinWinter
   end if
   ! compute average albedo
   call computeAlbedo(scalarSnowAlbedo(iGRU),refreshFactor(iGRU),decayFactor(iGRU),albedoMax,albedoMin(iGRU))
   ! assume albedo is the same in visible and near infra-red bands, and for direct and diffuse radiation
   spectralSnowAlbedoDiffuse(ixVisible,iGRU) = scalarSnowAlbedo(iGRU)
   spectralSnowAlbedoDiffuse(ixNearIR,iGRU)  = scalarSnowAlbedo(iGRU)
   spectralSnowAlbedoDirect(ixVisible,iGRU)  = scalarSnowAlbedo(iGRU)
   spectralSnowAlbedoDirect(ixNearIR,iGRU)   = scalarSnowAlbedo(iGRU)


  ! *** variable decay rate
  case(variableDecay)
   ! compute decay factor
   age1 = exp(-tempScalGrowth*(Tfreeze - surfaceTemp(1,iGRU) ))  ! temperature dependence
   age2 = age1**slushExp                                 ! increase with liquid water
   age3 = albedoSootLoad                                 ! soot loading
   decayFactor(iGRU) = dt*(age1 + age2 + age3)/albedoDecayRate
   ! compute diffuse albedo for the different spectral bands
   call computeAlbedo(spectralSnowAlbedoDiffuse(ixVisible,iGRU),refreshFactor(iGRU),decayFactor(iGRU),albedoMaxVisible,albedoMinVisible)
   call computeAlbedo(spectralSnowAlbedoDiffuse(ixNearIR,iGRU), refreshFactor(iGRU),decayFactor(iGRU),albedoMaxNearIR, albedoMinNearIR)
   ! compute factor to modify direct albedo at low zenith angles
   if(cosZenith(iGRU) < 0.5_rkind)then
    fZen = (1._rkind/bPar)*( ((1._rkind + bPar)/(1._rkind + 2._rkind*bPar*cosZenith(iGRU))) - 1._rkind)
   else
    fZen = 0._rkind
   end if
   ! compute direct albedo
   spectralSnowAlbedoDirect(ixVisible,iGRU) = spectralSnowAlbedoDiffuse(ixVisible,iGRU) + 0.4_rkind*fZen*(1._rkind - spectralSnowAlbedoDiffuse(ixVisible,iGRU))
   spectralSnowAlbedoDirect(ixNearIR,iGRU)  = spectralSnowAlbedoDiffuse(ixNearIR,iGRU)  + 0.4_rkind*fZen*(1._rkind - spectralSnowAlbedoDiffuse(ixNearIR,iGRU))

   ! compute average albedo
   scalarSnowAlbedo(iGRU) = (        Frad_direct)*(Frad_vis*spectralSnowAlbedoDirect(ixVisible,iGRU) + (1._rkind - Frad_vis)*spectralSnowAlbedoDirect(ixNearIR,iGRU) ) + &
                      (1._rkind - Frad_direct)*(Frad_vis*spectralSnowAlbedoDirect(ixVisible,iGRU) + (1._rkind - Frad_vis)*spectralSnowAlbedoDirect(ixNearIR,iGRU) )

  ! check that we identified the albedo option
!   case default; cycle;! err=20; message=trim(message)//'unable to identify option for snow albedo'; return

 end select  ! identify option for snow albedo

 ! check
!  if(scalarSnowAlbedo(iGRU) < 0._rkind)then; err=20; message=trim(message)//'unable to identify option for snow albedo'; return; end if
end do
 ! end association to data structures
 end associate

 end subroutine snowAlbedo


 ! *******************************************************************************************************
 ! private subroutine computeAlbedo: compute change in albedo -- implicit solution
 ! *******************************************************************************************************
 attributes(device) subroutine computeAlbedo(snowAlbedo,refreshFactor,decayFactor,albedoMax,albedoMin)
 implicit none
 ! dummy variables
 real(rkind),intent(inout)   :: snowAlbedo    ! snow albedo (-)
 real(rkind),intent(in)      :: refreshFactor ! albedo refreshment factor (-)
 real(rkind),intent(in)      :: decayFactor   ! albedo decay factor (-)
 real(rkind),intent(in)      :: albedoMax     ! maximum albedo (-)
 real(rkind),intent(in)      :: albedoMin     ! minimum albedo (-)
 ! local variables
 real(rkind)                 :: albedoChange ! change in albedo over the time step (-)
 ! compute change in albedo
 albedoChange = refreshFactor*(albedoMax - snowAlbedo) - (decayFactor*(snowAlbedo - albedoMin)) / (1._rkind + decayFactor)
 snowAlbedo   = snowAlbedo + albedoChange
 if(snowAlbedo > albedoMax) snowAlbedo = albedoMax
 end subroutine computeAlbedo


end module snowAlbedo_module
