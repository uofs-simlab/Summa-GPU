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

module vegSWavRad_module

! data types
USE nrtype
USE data_types,only:var_i            ! x%var(:)       (i4b)
USE data_types,only:var_dlength      ! x%var(:)%dat   (rkind)

! physical constants
USE multiconst,only:Tfreeze          ! temperature at freezing              (K)

! named variables for structure elements
USE var_lookup,only:iLookTYPE,iLookPROG,iLookDIAG,iLookFLUX

! model decisions
USE globalData,only:model_decisions  ! model decision structure
USE var_lookup,only:iLookDECISIONS   ! named variables for elements of the decision structure

! look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:        &
                      noah_mp,     & ! full Noah-MP implementation (including albedo)
                      CLM_2stream, & ! CLM 2-stream model (see CLM documentation)
                      UEB_2stream, & ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
                      NL_scatter,  & ! Simplified method Nijssen and Lettenmaier (JGR 1999)
                      BeersLaw       ! Beer's Law (as implemented in VIC)

! -------------------------------------------------------------------------------------------------
! privacy
implicit none
private
public::vegSWavRad
! dimensions
integer(i4b),parameter        :: nBands=2      ! number of spectral bands for shortwave radiation
! named variables
integer(i4b),parameter        :: ist     = 1   ! Surface type:  IST=1 => soil;  IST=2 => lake
integer(i4b),parameter        :: isc     = 4   ! Soil color type
integer(i4b),parameter        :: ice     = 0   ! Surface type:  ICE=0 => soil;  ICE=1 => sea-ice
! spatial indices
integer(i4b),parameter        :: iLoc    = 1   ! i-location
integer(i4b),parameter        :: jLoc    = 1   ! j-location
! algorithmic parameters
real(rkind),parameter            :: missingValue=-9999._rkind  ! missing value, used when diagnostic or state variables are undefined
real(rkind),parameter            :: verySmall=1.e-6_rkind   ! used as an additive constant to check if substantial difference among real numbers
real(rkind),parameter            :: mpe=1.e-6_rkind         ! prevents overflow error if division by zero
real(rkind),parameter            :: dx=1.e-6_rkind          ! finite difference increment
contains


 ! ************************************************************************************************
 ! public subroutine vegSWavRad: muster program to compute sw radiation in vegetation
 ! ************************************************************************************************
 subroutine vegSWavRad(&
                       dt,                           & ! intent(in):    time step (s) -- only used in Noah-MP radiation, to compute albedo
                       nSnow,                        & ! intent(in):    number of snow layers
                       nSoil,                        & ! intent(in):    number of soil layers
                       nLayers,                      & ! intent(in):    total number of layers
                       nGRU, &
                       computeVegFlux,               & ! intent(in):    logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                       decisions, veg_param, &
                       type_data,                    & ! intent(in):    classification of veg, soil etc. for a local HRU
                       prog_data,                    & ! intent(inout): model prognostic variables for a local HRU
                       diag_data,                    & ! intent(inout): model diagnostic variables for a local HRU
                       flux_data,                    & ! intent(inout): model flux variables
                       err,message)                    ! intent(out): error control
 ! external routines
!  USE NOAHMP_ROUTINES,only:radiation                                ! subroutine to calculate albedo and shortwave radiaiton in the canopy
 use device_data_types
 implicit none
 ! dummy variables
 real(rkind),intent(in)             :: dt                             ! time step (s) -- only used in Noah-MP radiation, to compute albedo
 integer(i4b),intent(in),device         :: nSnow(:)                          ! number of snow layers
 integer(i4b),intent(in)         :: nSoil                          ! number of soil layers
 integer(i4b),intent(in)         :: nLayers                        ! total number of layers
 integer(i4b) :: nGRU
 logical(lgt),intent(in),device         :: computeVegFlux(:)                 ! logical flag to compute vegetation fluxes (.false. if veg buried by snow)
 type(decisions_device) :: decisions
 type(veg_parameters) :: veg_param
 type(type_data_device),intent(in)          :: type_data                      ! classification of veg, soil etc. for a local HRU
 type(prog_data_device),intent(inout) :: prog_data                      ! model prognostic variables for a local HRU
 type(diag_data_device),intent(inout) :: diag_data                      ! model diagnostic variables for a local HRU
 type(flux_data_device),intent(inout) :: flux_data                      ! model flux variables
 integer(i4b),intent(out)        :: err                            ! error code
 character(*),intent(out)        :: message                        ! error message
 ! local variables
 character(LEN=256)              :: cmessage                       ! error message of downwind routine
 real(rkind),device                        :: snowmassPlusNewsnow(nGRU)            ! sum of snow mass and new snowfall (kg m-2 [mm])
 real(rkind),device                        :: scalarGroundSnowFraction(nGRU)       ! snow cover fraction on the ground surface (-)
 real(rkind),device              :: scalarVegFraction        ! vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
 real(rkind),device                        :: scalarTotalReflectedSolar(nGRU)      ! total reflected solar radiation (W m-2)
 real(rkind),device                        :: scalarTotalAbsorbedSolar(nGRU)       ! total absorbed solar radiation (W m-2)
 real(rkind),device                        :: scalarCanopyReflectedSolar(nGRU)     ! solar radiation reflected from the canopy (W m-2)
 real(rkind),device                        :: scalarGroundReflectedSolar(nGRU)     ! solar radiation reflected from the ground (W m-2)
 real(rkind),device                        :: scalarBetweenCanopyGapFraction(nGRU) ! between canopy gap fraction for beam (-)
 real(rkind),device                        :: scalarWithinCanopyGapFraction(nGRU)  ! within canopy gap fraction for beam (-)
 integer(i4b) :: iGRU
 type(dim3) :: blocks, threads

 scalarVegFraction=1._rkind
 ! ----------------------------------------------------------------------------------------------------------------------------------
 ! make association between local variables and the information in the data structures
 associate(&
  ! input: control
  vegTypeIndex_d               => type_data%vegTypeIndex,                            & ! intent(in): vegetation type index
  ix_canopySrad              => model_decisions(iLookDECISIONS%canopySrad)%iDecision,             & ! intent(in): index defining method for canopy shortwave radiation
  ! input: forcing at the upper boundary
  scalarSnowfall             => flux_data%scalarSnowfall,                   & ! intent(in): computed snowfall rate (kg m-2 s-1)
  spectralIncomingDirect     => flux_data%spectralIncomingDirect,    & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
  spectralIncomingDiffuse    => flux_data%spectralIncomingDiffuse,   & ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)
  ! input: snow states
  scalarSWE                  => prog_data%scalarSWE,                        & ! intent(in): snow water equivalent on the ground (kg m-2)
  scalarSnowDepth            => prog_data%scalarSnowDepth,                  & ! intent(in): snow depth on the ground surface (m)
  mLayerVolFracLiq           => prog_data%mLayerVolFracLiq,   & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
  spectralSnowAlbedoDiffuse  => prog_data%spectralSnowAlbedoDiffuse, & ! intent(in): diffuse albedo of snow in each spectral band (-)
  scalarSnowAlbedo           => prog_data%scalarSnowAlbedo,                 & ! intent(inout): snow albedo (-)
  ! input: ground and canopy temperature
  scalarGroundTemp           => prog_data%mLayerTemp,                       & ! intent(in): ground temperature (K)
  scalarCanopyTemp           => prog_data%scalarCanopyTemp,                 & ! intent(in): vegetation temperature (K)
  ! input: surface characteristix
  scalarSnowAge              => diag_data%scalarSnowAge,                    & ! intent(inout): non-dimensional snow age (-)
  scalarCosZenith            => diag_data%scalarCosZenith,                  & ! intent(in): cosine of the solar zenith angle (0-1)
  spectralSnowAlbedoDirect   => diag_data%spectralSnowAlbedoDirect,  & ! intent(in): direct albedo of snow in each spectral band (-)
  ! input: vegetation characteristix
  scalarExposedLAI           => diag_data%scalarExposedLAI,                 & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
  scalarExposedSAI           => diag_data%scalarExposedSAI,                 & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
  scalarCanopyWetFraction    => diag_data%scalarCanopyWetFraction,          & ! intent(in): canopy wetted fraction (-)
  ! output: canopy properties
  scalarCanopySunlitFraction => diag_data%scalarCanopySunlitFraction,       & ! intent(out): sunlit fraction of canopy (-)
  scalarCanopySunlitLAI      => diag_data%scalarCanopySunlitLAI,            & ! intent(out): sunlit leaf area (-)
  scalarCanopyShadedLAI      => diag_data%scalarCanopyShadedLAI,            & ! intent(out): shaded leaf area (-)
  spectralAlbGndDirect       => diag_data%spectralAlbGndDirect,                & ! intent(out): direct  albedo of underlying surface (1:nBands) (-)
  spectralAlbGndDiffuse      => diag_data%spectralAlbGndDiffuse,               & ! intent(out): diffuse albedo of underlying surface (1:nBands) (-)
  scalarGroundAlbedo         => diag_data%scalarGroundAlbedo,               & ! intent(out): albedo of the ground surface (-)
  ! output: canopy sw radiation fluxes
  scalarCanopySunlitPAR      => flux_data%scalarCanopySunlitPAR,            & ! intent(out): average absorbed par for sunlit leaves (w m-2)
  scalarCanopyShadedPAR      => flux_data%scalarCanopyShadedPAR,            & ! intent(out): average absorbed par for shaded leaves (w m-2)
  spectralBelowCanopyDirect  => flux_data%spectralBelowCanopyDirect,           & ! intent(out): downward direct flux below veg layer for each spectral band  W m-2)
  spectralBelowCanopyDiffuse => flux_data%spectralBelowCanopyDiffuse,          & ! intent(out): downward diffuse flux below veg layer for each spectral band (W m-2)
  scalarBelowCanopySolar     => flux_data%scalarBelowCanopySolar,           & ! intent(out): solar radiation transmitted below the canopy (W m-2)
  scalarCanopyAbsorbedSolar  => flux_data%scalarCanopyAbsorbedSolar,        & ! intent(out): solar radiation absorbed by canopy (W m-2)
  scalarGroundAbsorbedSolar  => flux_data%scalarGroundAbsorbedSolar,         & ! intent(out): solar radiation absorbed by ground (W m-2)
  rhol => veg_param % rhol,rhos=>veg_param%rhos,taul=>veg_param%taul,&
  taus => veg_param%taus, &
  opt_alb => veg_param%opt_alb, &
  omegas => veg_param%omegas, &
  opt_rad => veg_param%opt_rad, &
  tfrz => veg_param%tfrz, &
  betais => veg_param%betais, &
  betads => veg_param%betads, &
  xl => veg_param%xl, &
  hvt => veg_param%hvt, &
  hvb => veg_param%hvb, &
  rc => veg_param%rc, &
  swemx => veg_param%swemx, &
  albsat => veg_param%albsat, &
  albdry => veg_param%albdry, &
  alblak => veg_param%alblak &
 ) ! associating local variables with the information in the data structures
 ! -------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='vegSWavRad/'

 threads = dim3(128,1,1)
 blocks = dim3(nGRU/128+1,1,1)

 ! * preliminaries...
 ! ------------------

 !$cuf kernel do(1) <<<*,*>>>
 do iGRU=1,nGRU
 ! compute the sum of snow mass and new snowfall (kg m-2 [mm])
 snowmassPlusNewsnow(iGRU) = scalarSWE(iGRU) + scalarSnowfall(iGRU)*dt

 ! compute the ground snow fraction
 if(nSnow(iGRU) > 0)then
  scalarGroundSnowFraction(iGRU)  = 1._rkind
 else
  scalarGroundSnowFraction(iGRU)  = 0._rkind
 end if  ! (if there is snow on the ground)
end do

 ! * compute radiation fluxes...
 ! -----------------------------

 select case(ix_canopySrad)

  ! ***** unchanged Noah-MP routine
  case(noah_mp)

   call radiation_d<<<blocks,threads>>>(nGRU,nSnow,&
                  ! input
                  vegTypeIndex_d,                       & ! intent(in): vegetation type index
                  nSoil,                              & ! intent(in): number of soil layers
                  scalarSWE,                          & ! intent(in): snow water equivalent (kg m-2)
                  snowmassPlusNewsnow,                & ! intent(in): sum of snow mass and new snowfall (kg m-2 [mm])
                  dt,                                 & ! intent(in): time step (s)
                  scalarCosZenith,                    & ! intent(in): cosine of the solar zenith angle (0-1)
                  scalarSnowDepth,           & ! intent(in): snow depth on the ground surface (mm)
                  scalarGroundTemp,                   & ! intent(in): ground temperature (K)
                  scalarCanopyTemp,                   & ! intent(in): canopy temperature (K)
                  scalarGroundSnowFraction,           & ! intent(in): snow cover fraction (0-1)
                  scalarSnowfall,                     & ! intent(in): snowfall (kg m-2 s-1 [mm/s])
                  scalarCanopyWetFraction,            & ! intent(in): fraction of canopy that is wet
                  scalarExposedLAI,                   & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                  scalarExposedSAI,                   & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                  mLayerVolFracLiq,          & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                  spectralIncomingDirect,   & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
                  spectralIncomingDiffuse,  & ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)
                  scalarVegFraction,                  & ! intent(in): vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
                  iLoc, jLoc,                         & ! intent(in): spatial location indices
                  ! output
                  scalarSnowAlbedo,                   & ! intent(inout): snow albedo (-)
                  scalarSnowAge,                      & ! intent(inout): non-dimensional snow age (-)
                  scalarCanopySunlitFraction,         & ! intent(out): sunlit fraction of canopy (-)
                  scalarCanopySunlitLAI,              & ! intent(out): sunlit leaf area (-)
                  scalarCanopyShadedLAI,              & ! intent(out): shaded leaf area (-)
                  scalarCanopySunlitPAR,              & ! intent(out): average absorbed par for sunlit leaves (w m-2)
                  scalarCanopyShadedPAR,              & ! intent(out): average absorbed par for shaded leaves (w m-2)
                  scalarCanopyAbsorbedSolar,          & ! intent(out): solar radiation absorbed by canopy (W m-2)
                  scalarGroundAbsorbedSolar,          & ! intent(out): solar radiation absorbed by ground (W m-2)
                  scalarTotalReflectedSolar,          & ! intent(out): total reflected solar radiation (W m-2)
                  scalarTotalAbsorbedSolar,           & ! intent(out): total absorbed solar radiation (W m-2)
                  scalarCanopyReflectedSolar,         & ! intent(out): solar radiation reflected from the canopy (W m-2)
                  scalarGroundReflectedSolar,         & ! intent(out): solar radiation reflected from the ground (W m-2)
                  scalarBetweenCanopyGapFraction,     & ! intent(out): between canopy gap fraction for beam (-)
                  scalarWithinCanopyGapFraction,rhol,rhos,taul,taus, opt_alb,omegas,opt_rad,tfrz,betais,betads,xl,hvt,hvb,rc, swemx,albsat,albdry,alblak) ! intent(out): within canopy gap fraction for beam (-)

  ! **** all other options
  case(CLM_2stream,UEB_2stream,NL_scatter,BeersLaw)

   call canopy_SW<<<blocks,threads>>>(nGRU, nSnow,&
                  ! input: model control
                  vegTypeIndex_d,                                       & ! intent(in): index of vegetation type
                  isc,                                                & ! intent(in): index of soil type
                  computeVegFlux,                                     & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                  ix_canopySrad,                                      & ! intent(in): index of method used for transmission of shortwave rad through the canopy
                  ! input: model variables
                  scalarCosZenith,                                    & ! intent(in): cosine of direct zenith angle (0-1)
                  spectralIncomingDirect,                   & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
                  spectralIncomingDiffuse,                  & ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)
                  spectralSnowAlbedoDirect,                 & ! intent(in): direct albedo of snow in each spectral band (-)
                  spectralSnowAlbedoDiffuse,                & ! intent(in): diffuse albedo of snow in each spectral band (-)
                  scalarExposedLAI,                                   & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                  scalarExposedSAI,                                   & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                  scalarVegFraction,                                  & ! intent(in): vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
                  scalarCanopyWetFraction,                            & ! intent(in): fraction of lai, sai that is wetted (-)
                  scalarGroundSnowFraction,                           & ! intent(in): fraction of ground that is snow covered (-)
                  mLayerVolFracLiq,                                & ! intent(in): volumetric liquid water content in the upper-most soil layer (-)
                  scalarCanopyTemp,                                   & ! intent(in): canopy temperature (k)
                  ! output
                  spectralBelowCanopyDirect,                          & ! intent(out): downward direct flux below veg layer for each spectral band  W m-2)
                  spectralBelowCanopyDiffuse,                         & ! intent(out): downward diffuse flux below veg layer for each spectral band (W m-2)
                  scalarBelowCanopySolar,                             & ! intent(out): solar radiation transmitted below the canopy (W m-2)
                  spectralAlbGndDirect,                               & ! intent(out): direct  albedo of underlying surface (1:nBands) (-)
                  spectralAlbGndDiffuse,                              & ! intent(out): diffuse albedo of underlying surface (1:nBands) (-)
                  scalarGroundAlbedo,                                 & ! intent(out): albedo of the ground surface (-)
                  scalarCanopyAbsorbedSolar,                          & ! intent(out): solar radiation absorbed by the vegetation canopy (W m-2)
                  scalarGroundAbsorbedSolar,                          & ! intent(out): solar radiation absorbed by the ground (W m-2)
                  scalarCanopySunlitFraction,                         & ! intent(out): sunlit fraction of canopy (-)
                  scalarCanopySunlitLAI,                              & ! intent(out): sunlit leaf area (-)
                  scalarCanopyShadedLAI,                              & ! intent(out): shaded leaf area (-)
                  scalarCanopySunlitPAR,                              & ! intent(out): average absorbed par for sunlit leaves (w m-2)
                  scalarCanopyShadedPAR,                              & ! intent(out): average absorbed par for shaded leaves (w m-2)
                  rhol,albsat,albdry,rhos,taul,taus,omegas,opt_rad,tfrz,betais,betads,xl,hvt,hvb,rc)                                         ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

  case default; err=20; message=trim(message)//'unable to identify option for canopy sw radiation'; return

 end select ! (option for canopy sw radiation)

 ! end association between local variables and the information in the data structures
 end associate

 end subroutine vegSWavRad



 ! ************************************************************************************************
 ! private subroutine canopy_SW: various options to compute canopy sw radiation fluxes
 ! ************************************************************************************************
 attributes(global) subroutine canopy_SW(nGRU, nSnow, &
                      ! input: model control
                      vegTypeIndex_,                                       & ! intent(in): index of vegetation type
                      isc,                                                & ! intent(in): index of soil color
                      computeVegFlux,                                     & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                      ix_canopySrad,                                      & ! intent(in): index of method used for transmission of shortwave rad through the canopy
                      ! input: model variables
                      scalarCosZenith_,                                    & ! intent(in): cosine of direct zenith angle (0-1)
                      spectralIncomingDirect_,                             & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
                      spectralIncomingDiffuse_,                            & ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)
                      spectralSnowAlbedoDirect_,                           & ! intent(in): direct albedo of snow in each spectral band (-)
                      spectralSnowAlbedoDiffuse_,                          & ! intent(in): diffuse albedo of snow in each spectral band (-)
                      scalarExposedLAI_,                                   & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                      scalarExposedSAI_,                                   & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                      scalarVegFraction,                                  & ! intent(in): vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
                      scalarCanopyWetFraction_,                            & ! intent(in): fraction of lai, sai that is wetted (-)
                      scalarGroundSnowFraction_,                           & ! intent(in): fraction of ground that is snow covered (-)
                      scalarVolFracLiqUpper_,                              & ! intent(in): volumetric liquid water content in the upper-most soil layer (-)
                      scalarCanopyTempTrial_,                              & ! intent(in): canopy temperature (K)
                      ! output
                      spectralBelowCanopyDirect_,                          & ! intent(out): downward direct flux below veg layer (W m-2)
                      spectralBelowCanopyDiffuse_,                         & ! intent(out): downward diffuse flux below veg layer (W m-2)
                      scalarBelowCanopySolar_,                             & ! intent(out): radiation transmitted below the canopy (W m-2)
                      spectralAlbGndDirect_,                               & ! intent(out): direct  albedo of underlying surface (1:nBands) (-)
                      spectralAlbGndDiffuse_,                              & ! intent(out): diffuse albedo of underlying surface (1:nBands) (-)
                      scalarGroundAlbedo_,                                 & ! intent(out): albedo of the ground surface (-)
                      scalarCanopyAbsorbedSolar_,                          & ! intent(out): radiation absorbed by the vegetation canopy (W m-2)
                      scalarGroundAbsorbedSolar_,                          & ! intent(out): radiation absorbed by the ground (W m-2)
                      scalarCanopySunlitFraction_,                         & ! intent(out): sunlit fraction of canopy (-)
                      scalarCanopySunlitLAI_,                              & ! intent(out): sunlit leaf area (-)
                      scalarCanopyShadedLAI_,                              & ! intent(out): shaded leaf area (-)
                      scalarCanopySunlitPAR_,                              & ! intent(out): average absorbed par for sunlit leaves (w m-2)
                      scalarCanopyShadedPAR_,                              & ! intent(out): average absorbed par for shaded leaves (w m-2)
                      rhol,albsat,albdry,rhos,taul,taus,omegas,opt_rad,tfrz,betais,betads,xl,hvt,hvb,rc)                                          ! intent(out): error control
 ! utilities
 USE expIntegral_module,only:expInt                                          ! function to calculate the exponential integral
 implicit none
 ! Noah-MP modules
!  USE NOAHMP_ROUTINES,only:twoStream                                          ! two-stream radiative transfer
 ! Noah vegetation tables
!  USE NOAHMP_VEG_PARAMETERS, only: RHOS,RHOL                                  ! Noah-MP: stem and leaf reflectance for each wave band
!  USE NOAHMP_VEG_PARAMETERS, only: TAUS,TAUL                                  ! Noah-MP: stem and leaf transmittance for each wave band
 ! input
 integer(i4b),value :: nGRU
 integer(i4b) :: nSnow(:)
 integer(i4b),intent(in)        :: vegTypeIndex_(:)                              ! vegetation type index
 integer(i4b),intent(in),value        :: isc                                       ! soil color index
 logical(lgt),intent(in)        :: computeVegFlux(:)                            ! logical flag to compute vegetation fluxes (.false. if veg buried by snow)
 integer(i4b),intent(in),value        :: ix_canopySrad                             ! choice of canopy shortwave radiation method
 real(rkind),intent(in)            :: scalarCosZenith_(:)                           ! cosine of the solar zenith angle (0-1)
 real(rkind),intent(in)            :: spectralIncomingDirect_(:,:)                 ! incoming direct solar radiation in each wave band (w m-2)
 real(rkind),intent(in)            :: spectralIncomingDiffuse_(:,:)                ! incoming diffuse solar radiation in each wave band (w m-2)
 real(rkind),intent(in)            :: spectralSnowAlbedoDirect_(:,:)               ! direct albedo of snow in each spectral band (-)
 real(rkind),intent(in)            :: spectralSnowAlbedoDiffuse_(:,:)              ! diffuse albedo of snow in each spectral band (-)
 real(rkind),intent(in)            :: scalarExposedLAI_(:)                          ! exposed leaf area index after burial by snow (m2 m-2)
 real(rkind),intent(in)            :: scalarExposedSAI_(:)                          ! exposed stem area index after burial by snow (m2 m-2)
 real(rkind),intent(in)            :: scalarVegFraction                         ! vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
 real(rkind),intent(in)            :: scalarCanopyWetFraction_(:)                   ! fraction of canopy that is wet (-)
 real(rkind),intent(in)            :: scalarGroundSnowFraction_(:)                  ! fraction of ground that is snow covered (-)
 real(rkind),intent(in)            :: scalarVolFracLiqUpper_(:,:)                     ! volumetric liquid water content in the upper-most soil layer (-)
 real(rkind),intent(in)            :: scalarCanopyTempTrial_(:)                     ! trial value of canopy temperature (K)
 ! output
 real(rkind),intent(out)           :: spectralBelowCanopyDirect_(:,:)              ! downward direct flux below veg layer (W m-2)
 real(rkind),intent(out)           :: spectralBelowCanopyDiffuse_(:,:)             ! downward diffuse flux below veg layer (W m-2)
 real(rkind),intent(out)           :: scalarBelowCanopySolar_(:)                    ! radiation transmitted below the canopy (W m-2)
 real(rkind),intent(out)           :: spectralAlbGndDirect_(:,:)                   ! direct  albedo of underlying surface (1:nBands) (-)
 real(rkind),intent(out)           :: spectralAlbGndDiffuse_(:,:)                  ! diffuse albedo of underlying surface (1:nBands) (-)
 real(rkind),intent(out)           :: scalarGroundAlbedo_(:)                        ! albedo of the ground surface (-)
 real(rkind),intent(out)           :: scalarCanopyAbsorbedSolar_(:)                 ! radiation absorbed by the vegetation canopy (W m-2)
 real(rkind),intent(out)           :: scalarGroundAbsorbedSolar_(:)                 ! radiation absorbed by the ground (W m-2)
 real(rkind),intent(out)           :: scalarCanopySunlitFraction_(:)                ! sunlit fraction of canopy (-)
 real(rkind),intent(out)           :: scalarCanopySunlitLAI_(:)                     ! sunlit leaf area (-)
 real(rkind),intent(out)           :: scalarCanopyShadedLAI_(:)                     ! shaded leaf area (-)
 real(rkind),intent(out)           :: scalarCanopySunlitPAR_(:)                     ! average absorbed par for sunlit leaves (w m-2)
 real(rkind),intent(out)           :: scalarCanopyShadedPAR_(:)                     ! average absorbed par for shaded leaves (w m-2)
 real(rkind) :: rhol(:,:),albsat(:,:),albdry(:,:),rhos(:,:),taul(:,:),taus(:,:)
 real(rkind) :: omegas(:)
integer(i4b) :: opt_rad
real(rkind) :: tfrz,betais,betads
real(rkind) :: xl(:),hvt(:),hvb(:),rc(:)
 integer(i4b)       :: err                                       ! error code
!  character(*),intent(out)       :: message                                   ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! general
 integer(i4b),parameter                :: ixVisible=1                        ! index of the visible wave band
 integer(i4b),parameter                :: ixNearIR=2                         ! index of the near infra-red wave band
 integer(i4b)                          :: iBand                              ! index of wave band
 integer(i4b)                          :: ic                                 ! 0=unit incoming direct; 1=unit incoming diffuse
!  character(LEN=256)                    :: cmessage                           ! error message of downwind routine
 ! variables used in Nijssen-Lettenmaier method
 real(rkind),parameter                    :: multScatExp=0.81_rkind                ! multiple scattering exponent (-)
 real(rkind),parameter                    :: bulkCanopyAlbedo=0.25_rkind           ! bulk canopy albedo (-), smaller than actual canopy albedo because of shading in the canopy
 real(rkind),dimension(1:nBands)          :: spectralIncomingSolar              ! total incoming solar radiation in each spectral band (W m-2)
 real(rkind),dimension(1:nBands)          :: spectralGroundAbsorbedDirect       ! total direct radiation absorbed at the ground surface (W m-2)
 real(rkind),dimension(1:nBands)          :: spectralGroundAbsorbedDiffuse      ! total diffuse radiation absorbed at the ground surface (W m-2)
 real(rkind)                              :: Fdirect                            ! fraction of direct radiation (-)
 real(rkind)                              :: tauInitial                         ! transmission in the absence of scattering and multiple reflections (-)
 real(rkind)                              :: tauTotal                           ! transmission due to scattering and multiple reflections (-)
 ! variables used in Mahat-Tarboton method
 real(rkind),parameter                    :: Frad_vis=0.5_rkind                    ! fraction of radiation in the visible wave band (-)
 real(rkind),parameter                    :: gProjParam=0.5_rkind                  ! projected leaf and stem area in the solar direction (-)
 real(rkind),parameter                    :: bScatParam=0.5_rkind                  ! back scatter parameter (-)
 real(rkind)                              :: transCoef                          ! transmission coefficient (-)
 real(rkind)                              :: transCoefPrime                     ! "k-prime" coefficient (-)
 real(rkind)                              :: groundAlbedoDirect                 ! direct ground albedo (-)
 real(rkind)                              :: groundAlbedoDiffuse                ! diffuse ground albedo (-)
 real(rkind)                              :: tauInfinite                        ! direct transmission for an infinite canopy (-)
 real(rkind)                              :: betaInfinite                       ! direct upward reflection factor for an infinite canopy (-)
 real(rkind)                              :: tauFinite                          ! direct transmission for a finite canopy (-)
 real(rkind)                              :: betaFinite                         ! direct reflectance for a finite canopy (-)
 real(rkind)                              :: vFactor                            ! scaled vegetation area used to compute diffuse radiation (-)
 real(rkind)                              :: expi                               ! exponential integral (-)
 real(rkind)                              :: taudInfinite                       ! diffuse transmission for an infinite canopy (-)
 real(rkind)                              :: taudFinite                         ! diffuse transmission for a finite canopy (-)
 real(rkind)                              :: betadFinite                        ! diffuse reflectance for a finite canopy (-)
 real(rkind)                              :: refMult                            ! multiple reflection factor (-)
 real(rkind)                              :: fracRadAbsDown                     ! fraction of radiation absorbed by vegetation on the way down
 real(rkind)                              :: fracRadAbsUp                       ! fraction of radiation absorbed by vegetation on the way up
 real(rkind)                              :: tauDirect                          ! total transmission of direct radiation (-)
 real(rkind)                              :: tauDiffuse                         ! total transmission of diffuse radiation (-)
 real(rkind)                              :: fractionRefDirect                  ! fraction of direct radiaiton lost to space (-)
 real(rkind)                              :: fractionRefDiffuse                 ! fraction of diffuse radiaiton lost to space (-)
 real(rkind),dimension(1:nBands)          :: spectralBelowCanopySolar           ! total below-canopy radiation for each wave band (W m-2)
 real(rkind),dimension(1:nBands)          :: spectralTotalReflectedSolar        ! total reflected radiaion for each wave band (W m-2)
 real(rkind),dimension(1:nBands)          :: spectralGroundAbsorbedSolar        ! radiation absorbed by the ground in each wave band (W m-2)
 real(rkind),dimension(1:nBands)          :: spectralCanopyAbsorbedSolar        ! radiation absorbed by the canopy in each wave band (W m-2)
 ! vegetation properties used in 2-stream
 real(rkind)                              :: scalarExposedVAI                   ! one-sided leaf+stem area index (m2/m2)
 real(rkind)                              :: weightLeaf                         ! fraction of exposed VAI that is leaf
 real(rkind)                              :: weightStem                         ! fraction of exposed VAI that is stem
 real(rkind),dimension(1:nBands)          :: spectralVegReflc                   ! leaf+stem reflectance (1:nbands)
 real(rkind),dimension(1:nBands)          :: spectralVegTrans                   ! leaf+stem transmittance (1:nBands)
 ! output from two-stream -- direct-beam
 real(rkind),dimension(1:nBands)          :: spectralCanopyAbsorbedDirect       ! flux abs by veg layer (per unit incoming flux), (1:nBands)
 real(rkind),dimension(1:nBands)          :: spectralTotalReflectedDirect       ! flux refl above veg layer (per unit incoming flux), (1:nBands)
 real(rkind),dimension(1:nBands)          :: spectralDirectBelowCanopyDirect    ! down dir flux below veg layer (per unit in flux), (1:nBands)
 real(rkind),dimension(1:nBands)          :: spectralDiffuseBelowCanopyDirect   ! down dif flux below veg layer (per unit in flux), (1:nBands)
 real(rkind),dimension(1:nBands)          :: spectralCanopyReflectedDirect      ! flux reflected by veg layer   (per unit incoming flux), (1:nBands)
 real(rkind),dimension(1:nBands)          :: spectralGroundReflectedDirect      ! flux reflected by ground (per unit incoming flux), (1:nBands)
 ! output from two-stream -- diffuse
 real(rkind),dimension(1:nBands)          :: spectralCanopyAbsorbedDiffuse      ! flux abs by veg layer (per unit incoming flux), (1:nBands)
 real(rkind),dimension(1:nBands)          :: spectralTotalReflectedDiffuse      ! flux refl above veg layer (per unit incoming flux), (1:nBands)
 real(rkind),dimension(1:nBands)          :: spectralDirectBelowCanopyDiffuse   ! down dir flux below veg layer (per unit in flux), (1:nBands)
 real(rkind),dimension(1:nBands)          :: spectralDiffuseBelowCanopyDiffuse  ! down dif flux below veg layer (per unit in flux), (1:nBands)
 real(rkind),dimension(1:nBands)          :: spectralCanopyReflectedDiffuse     ! flux reflected by veg layer   (per unit incoming flux), (1:nBands)
 real(rkind),dimension(1:nBands)          :: spectralGroundReflectedDiffuse     ! flux reflected by ground (per unit incoming flux), (1:nBands)
 ! output from two-stream -- scalar variables
 real(rkind)                              :: scalarGproj                        ! projected leaf+stem area in solar direction
 real(rkind)                              :: scalarBetweenCanopyGapFraction     ! between canopy gap fraction for beam (-)
 real(rkind)                              :: scalarWithinCanopyGapFraction      ! within canopy gap fraction for beam (-)
 ! radiation fluxes
 real(rkind)                              :: ext                                ! optical depth of direct beam per unit leaf + stem area
 real(rkind)                              :: scalarCanopyShadedFraction         ! shaded fraction of the canopy
 real(rkind)                              :: fractionLAI                        ! fraction of vegetation that is leaves
 real(rkind)                              :: visibleAbsDirect                   ! direct-beam radiation absorbed in the visible part of the spectrum (W m-2)
 real(rkind)                              :: visibleAbsDiffuse                  ! diffuse radiation absorbed in the visible part of the spectrum (W m-2)
 integer(i4b) :: iGRU
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0
 iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
 if (iGRU.gt.nGRU) return
  

 ! compute the albedo of the ground surface
 call gndAlbedo(&
                ! input
                isc,                                   & ! intent(in): index of soil color
                scalarGroundSnowFraction_(iGRU),              & ! intent(in): fraction of ground that is snow covered (-)
                scalarVolFracLiqUpper_(nSnow(iGRU)+1,iGRU),                 & ! intent(in): volumetric liquid water content in upper-most soil layer (-)
                spectralSnowAlbedoDirect_(:,iGRU),              & ! intent(in): direct albedo of snow in each spectral band (-)
                spectralSnowAlbedoDiffuse_(:,iGRU),             & ! intent(in): diffuse albedo of snow in each spectral band (-)
                ! output
                spectralAlbGndDirect_(:,iGRU),                  & ! intent(out): direct  albedo of underlying surface (-)
                spectralAlbGndDiffuse_(:,iGRU),                 & ! intent(out): diffuse albedo of underlying surface (-)
                err,albsat,albdry)                             ! intent(out): error control
 if(err/=0)then;  return; end if

 ! initialize accumulated fluxes
 scalarBelowCanopySolar_(iGRU)    = 0._rkind  ! radiation transmitted below the canopy (W m-2)
 scalarCanopyAbsorbedSolar_(iGRU) = 0._rkind  ! radiation absorbed by the vegetation canopy (W m-2)
 scalarGroundAbsorbedSolar_(iGRU) = 0._rkind  ! radiation absorbed by the ground (W m-2)

 ! check for an early return (no radiation or no exposed canopy)
 if(.not.computeVegFlux(iGRU) .or. scalarCosZenith_(iGRU) < tiny(scalarCosZenith_(iGRU)))then
  ! set canopy radiation to zero
  scalarCanopySunlitFraction_(iGRU) = 0._rkind                ! sunlit fraction of canopy (-)
  scalarCanopySunlitLAI_(iGRU)      = 0._rkind                ! sunlit leaf area (-)
  scalarCanopyShadedLAI_(iGRU)      = scalarExposedLAI_(iGRU)     ! shaded leaf area (-)
  scalarCanopySunlitPAR_(iGRU)      = 0._rkind                ! average absorbed par for sunlit leaves (w m-2)
  scalarCanopyShadedPAR_(iGRU)      = 0._rkind                ! average absorbed par for shaded leaves (w m-2)
  ! compute below-canopy radiation
  do iBand=1,nBands
   ! (set below-canopy radiation to incoming radiation)
   if(scalarCosZenith_(iGRU) > tiny(scalarCosZenith_(iGRU)))then
    spectralBelowCanopyDirect_(iBand,iGRU)  = spectralIncomingDirect_(iBand,iGRU)
    spectralBelowCanopyDiffuse_(iBand,iGRU) = spectralIncomingDiffuse_(iBand,iGRU)
   else
    spectralBelowCanopyDirect_(iBand,iGRU)  = 0._rkind
    spectralBelowCanopyDiffuse_(iBand,iGRU) = 0._rkind
   end if
   ! (accumulate radiation transmitted below the canopy)
   scalarBelowCanopySolar_(iGRU)    = scalarBelowCanopySolar_(iGRU) + &                                                  ! contribution from all previous wave bands
                               spectralBelowCanopyDirect_(iBand,iGRU) + spectralBelowCanopyDiffuse_(iBand,iGRU)        ! contribution from current wave band
   ! (accumulate radiation absorbed by the ground)
   scalarGroundAbsorbedSolar_(iGRU) = scalarGroundAbsorbedSolar_(iGRU) + &                                               ! contribution from all previous wave bands
                               spectralBelowCanopyDirect_(iBand,iGRU)*(1._rkind - spectralAlbGndDirect_(iBand,iGRU)) + &  ! direct radiation from current wave band
                               spectralBelowCanopyDiffuse_(iBand,iGRU)*(1._rkind - spectralAlbGndDiffuse_(iBand,iGRU))    ! diffuse radiation from current wave band
  end do  ! looping through wave bands
  return
 end if

 ! compute exposed leaf and stem area index
 scalarExposedVAI = scalarExposedLAI_(iGRU) + scalarExposedSAI_(iGRU)
 if(scalarExposedVAI < epsilon(scalarExposedVAI))then; err=20; return; end if

 ! ============================================================================================================================================================
 ! ============================================================================================================================================================

 ! different options for radiation transmission
 select case(ix_canopySrad)

  ! -----------------------------------------------------------------------------------------------------------------------------------------------------------
  ! Beer's Law
  case(BeersLaw)

   ! define transmission coefficient (-)
   scalarGproj = gProjParam
   transCoef   = scalarGproj/scalarCosZenith_(iGRU)

   ! compute transmission of direct radiation according to Beer's Law (-)
   tauTotal = exp(-transCoef*scalarExposedVAI)
   !print*, 'tauTotal = ', tauTotal

   ! compute ground albedo (-)
   groundAlbedoDirect  = Frad_vis*spectralAlbGndDirect_(ixVisible,iGRU)  + (1._rkind - Frad_vis)*spectralAlbGndDirect_(ixNearIR,iGRU)
   groundAlbedoDiffuse = Frad_vis*spectralAlbGndDiffuse_(ixVisible,iGRU) + (1._rkind - Frad_vis)*spectralAlbGndDiffuse_(ixNearIR,iGRU)

   ! compute radiation in each spectral band (W m-2)
   do iBand=1,nBands

    ! compute total incoming solar radiation
    spectralIncomingSolar(iBand) = spectralIncomingDirect_(iBand,iGRU) + spectralIncomingDiffuse_(iBand,iGRU)

    ! compute fraction of direct radiation
    Fdirect = spectralIncomingDirect_(iBand,iGRU) / (spectralIncomingSolar(iBand) + verySmall)
    if(Fdirect < 0._rkind .or. Fdirect > 1._rkind)then
     print*, 'spectralIncomingDirect_(iBand,iGRU) = ', spectralIncomingDirect_(iBand,iGRU)
     print*, 'spectralIncomingSolar(iBand)  = ', spectralIncomingSolar(iBand)
     print*, 'Fdirect = ', Fdirect
     err=20; return
    end if

    ! compute ground albedo (-)
    scalarGroundAlbedo_(iGRU)  = Fdirect*groundAlbedoDirect + (1._rkind - Fdirect)*groundAlbedoDiffuse
    if(scalarGroundAlbedo_(iGRU) < 0._rkind .or. scalarGroundAlbedo_(iGRU) > 1._rkind)then
     print*, 'groundAlbedoDirect = ',  groundAlbedoDirect
     print*, 'groundAlbedoDiffuse = ', groundAlbedoDiffuse
     err=20; return
    end if

    ! compute below-canopy radiation (W m-2)
    spectralBelowCanopyDirect_(iBand,iGRU)  = spectralIncomingDirect_(iBand,iGRU)*tauTotal              ! direct radiation from current wave band
    spectralBelowCanopyDiffuse_(iBand,iGRU) = spectralIncomingDiffuse_(iBand,iGRU)*tauTotal             ! diffuse radiation from current wave band
    spectralBelowCanopySolar(iBand)   = spectralBelowCanopyDirect_(iBand,iGRU) + spectralBelowCanopyDiffuse_(iBand,iGRU)

    ! compute radiation absorbed by the ground in given wave band (W m-2)
    spectralGroundAbsorbedDirect(iBand)  = (1._rkind - scalarGroundAlbedo_(iGRU))*spectralBelowCanopyDirect_(iBand,iGRU)
    spectralGroundAbsorbedDiffuse(iBand) = (1._rkind - scalarGroundAlbedo_(iGRU))*spectralBelowCanopyDiffuse_(iBand,iGRU)
    spectralGroundAbsorbedSolar(iBand)   = spectralGroundAbsorbedDirect(iBand) + spectralGroundAbsorbedDiffuse(iBand)

    ! compute radiation absorbed by vegetation in current wave band (W m-2)
    fracRadAbsDown = (1._rkind - tauTotal)*(1._rkind - bulkCanopyAlbedo)                            ! (fraction of radiation absorbed on the way down)
    fracRadAbsUp   = tauTotal*scalarGroundAlbedo_(iGRU)*(1._rkind - tauTotal)   ! (fraction of radiation absorbed on the way up)
    spectralCanopyAbsorbedDirect(iBand)  = spectralIncomingDirect_(iBand,iGRU)*(fracRadAbsDown + fracRadAbsUp)
    spectralCanopyAbsorbedDiffuse(iBand) = spectralIncomingDiffuse_(iBand,iGRU)*(fracRadAbsDown + fracRadAbsUp)
    spectralCanopyAbsorbedSolar(iBand)   = spectralCanopyAbsorbedDirect(iBand) + spectralCanopyAbsorbedDiffuse(iBand)
    ! (check)
    if(spectralCanopyAbsorbedDirect(iBand) > spectralIncomingDirect_(iBand,iGRU) .or. spectralCanopyAbsorbedDiffuse(iBand) > spectralIncomingDiffuse_(iBand,iGRU))then
     print*, 'tauTotal = ', tauTotal
     print*, 'bulkCanopyAlbedo = ', bulkCanopyAlbedo
     print*, 'scalarGroundAlbedo = ', scalarGroundAlbedo_(iGRU)
     err=20; return
    end if

    ! compute solar radiation lost to space in given wave band (W m-2)
    spectralTotalReflectedDirect(iBand)  = spectralIncomingDirect_(iBand,iGRU) - spectralGroundAbsorbedDirect(iBand) - spectralCanopyAbsorbedDirect(iBand)
    spectralTotalReflectedDiffuse(iBand) = spectralIncomingDiffuse_(iBand,iGRU) - spectralGroundAbsorbedDiffuse(iBand) - spectralCanopyAbsorbedDiffuse(iBand)
    spectralTotalReflectedSolar(iBand)   = spectralTotalReflectedDirect(iBand) + spectralTotalReflectedDiffuse(iBand)
    if(spectralTotalReflectedDirect(iBand) < 0._rkind .or. spectralTotalReflectedDiffuse(iBand) < 0._rkind)then
     print*, 'scalarGroundAlbedo = ', scalarGroundAlbedo_(iGRU)
     print*, 'tauTotal = ', tauTotal
     print*, 'fracRadAbsDown = ', fracRadAbsDown
     print*, 'fracRadAbsUp = ', fracRadAbsUp
     print*, 'spectralBelowCanopySolar(iBand) = ', spectralBelowCanopySolar(iBand)
     print*, 'spectralGroundAbsorbedSolar(iBand) = ', spectralGroundAbsorbedSolar(iBand)
     print*, 'spectralCanopyAbsorbedSolar(iBand) = ', spectralCanopyAbsorbedSolar(iBand)
     err=20; return
    end if

    ! save canopy radiation absorbed in visible wavelengths
    if(iBand == ixVisible)then
     visibleAbsDirect  = spectralCanopyAbsorbedDirect(ixVisible)
     visibleAbsDiffuse = spectralCanopyAbsorbedDiffuse(ixVisible)
    end if

    ! accumulate fluxes
    scalarBelowCanopySolar_(iGRU)    = scalarBelowCanopySolar_(iGRU) + spectralBelowCanopySolar(iBand)
    scalarGroundAbsorbedSolar_(iGRU) = scalarGroundAbsorbedSolar_(iGRU) + spectralGroundAbsorbedSolar(iBand)
    scalarCanopyAbsorbedSolar_(iGRU) = scalarCanopyAbsorbedSolar_(iGRU) + spectralCanopyAbsorbedSolar(iBand)

   end do  ! (looping through spectral bands)


  ! -----------------------------------------------------------------------------------------------------------------------------------------------------------
  ! method of Nijssen and Lettenmaier (JGR, 1999)
  case(NL_scatter)

   ! define transmission coefficient (-)
   scalarGproj = gProjParam
   transCoef   = scalarGproj/scalarCosZenith_(iGRU)

   ! compute transmission of direct radiation according to Beer's Law (-)
   tauFinite = exp(-transCoef*scalarExposedVAI)

   ! compute transmission of diffuse radiation (-)
   vFactor    = scalarGproj*scalarExposedVAI
   expi       = expInt(vFactor)
   taudFinite = (1._rkind - vFactor)*exp(-vFactor) + (vFactor**2_i4b)*expi

   ! compute ground albedo (-)
   groundAlbedoDirect  = Frad_vis*spectralAlbGndDirect_(ixVisible,iGRU)  + (1._rkind - Frad_vis)*spectralAlbGndDirect_(ixNearIR,iGRU)
   groundAlbedoDiffuse = Frad_vis*spectralAlbGndDiffuse_(ixVisible,iGRU) + (1._rkind - Frad_vis)*spectralAlbGndDiffuse_(ixNearIR,iGRU)

   ! compute radiation in each spectral band (W m-2)
   do iBand=1,nBands

    ! compute total incoming solar radiation
    spectralIncomingSolar(iBand) = spectralIncomingDirect_(iBand,iGRU) + spectralIncomingDiffuse_(iBand,iGRU)

    ! compute fraction of direct radiation
    Fdirect = spectralIncomingDirect_(iBand,iGRU) / (spectralIncomingSolar(iBand) + verySmall)
    if(Fdirect < 0._rkind .or. Fdirect > 1._rkind)then
     print*, 'spectralIncomingDirect_(iBand,iGRU) = ', spectralIncomingDirect_(iBand,iGRU)
     print*, 'spectralIncomingSolar(iBand)  = ', spectralIncomingSolar(iBand)
     print*, 'Fdirect = ', Fdirect
     err=20; return
    end if

    ! compute ground albedo (-)
    scalarGroundAlbedo_(iGRU)  = Fdirect*groundAlbedoDirect + (1._rkind - Fdirect)*groundAlbedoDiffuse
    if(scalarGroundAlbedo_(iGRU) < 0._rkind .or. scalarGroundAlbedo_(iGRU) > 1._rkind)then
     print*, 'groundAlbedoDirect = ',  groundAlbedoDirect
     print*, 'groundAlbedoDiffuse = ', groundAlbedoDiffuse
     err=20; return
    end if

    ! compute initial transmission in the absence of scattering and multiple reflections (-)
    tauInitial = Fdirect*tauFinite + (1._rkind - Fdirect)*taudFinite

    ! compute increase in transmission due to scattering (-)
    tauTotal = (tauInitial**multScatExp)

    ! compute multiple reflections factor
    refMult = 1._rkind / (1._rkind - scalarGroundAlbedo_(iGRU)*bulkCanopyAlbedo*(1._rkind - taudFinite**multScatExp) )

    ! compute below-canopy radiation (W m-2)
    spectralBelowCanopyDirect_(iBand,iGRU)  = spectralIncomingDirect_(iBand,iGRU)*tauTotal*refMult              ! direct radiation from current wave band
    spectralBelowCanopyDiffuse_(iBand,iGRU) = spectralIncomingDiffuse_(iBand,iGRU)*tauTotal*refMult             ! diffuse radiation from current wave band
    spectralBelowCanopySolar(iBand)   = spectralBelowCanopyDirect_(iBand,iGRU) + spectralBelowCanopyDiffuse_(iBand,iGRU)

    ! compute radiation absorbed by the ground in given wave band (W m-2)
    spectralGroundAbsorbedDirect(iBand)  = (1._rkind - scalarGroundAlbedo_(iGRU))*spectralBelowCanopyDirect_(iBand,iGRU)
    spectralGroundAbsorbedDiffuse(iBand) = (1._rkind - scalarGroundAlbedo_(iGRU))*spectralBelowCanopyDiffuse_(iBand,iGRU)
    spectralGroundAbsorbedSolar(iBand)   = spectralGroundAbsorbedDirect(iBand) + spectralGroundAbsorbedDiffuse(iBand)

    ! compute radiation absorbed by vegetation in current wave band (W m-2)
    fracRadAbsDown = (1._rkind - tauTotal)*(1._rkind - bulkCanopyAlbedo)                            ! (fraction of radiation absorbed on the way down)
    fracRadAbsUp   = tauTotal*refMult*scalarGroundAlbedo_(iGRU)*(1._rkind - taudFinite**multScatExp)   ! (fraction of radiation absorbed on the way up)
    spectralCanopyAbsorbedDirect(iBand)  = spectralIncomingDirect_(iBand,iGRU)*(fracRadAbsDown + fracRadAbsUp)
    spectralCanopyAbsorbedDiffuse(iBand) = spectralIncomingDiffuse_(iBand,iGRU)*(fracRadAbsDown + fracRadAbsUp)
    spectralCanopyAbsorbedSolar(iBand)   = spectralCanopyAbsorbedDirect(iBand) + spectralCanopyAbsorbedDiffuse(iBand)

    ! compute solar radiation lost to space in given wave band (W m-2)
    spectralTotalReflectedDirect(iBand)  = spectralIncomingDirect_(iBand,iGRU) - spectralGroundAbsorbedDirect(iBand) - spectralCanopyAbsorbedDirect(iBand)
    spectralTotalReflectedDiffuse(iBand) = spectralIncomingDiffuse_(iBand,iGRU) - spectralGroundAbsorbedDiffuse(iBand) - spectralCanopyAbsorbedDiffuse(iBand)
    spectralTotalReflectedSolar(iBand)   = spectralTotalReflectedDirect(iBand) + spectralTotalReflectedDiffuse(iBand)
    if(spectralTotalReflectedDirect(iBand) < 0._rkind .or. spectralTotalReflectedDiffuse(iBand) < 0._rkind)then
     err=20; return
    end if

    ! save canopy radiation absorbed in visible wavelengths
    if(iBand == ixVisible)then
     visibleAbsDirect  = spectralCanopyAbsorbedDirect(ixVisible)
     visibleAbsDiffuse = spectralCanopyAbsorbedDiffuse(ixVisible)
    end if

    ! accumulate fluxes
    scalarBelowCanopySolar_(iGRU)    = scalarBelowCanopySolar_(iGRU) + spectralBelowCanopySolar(iBand)
    scalarGroundAbsorbedSolar_(iGRU) = scalarGroundAbsorbedSolar_(iGRU) + spectralGroundAbsorbedSolar(iBand)
    scalarCanopyAbsorbedSolar_(iGRU) = scalarCanopyAbsorbedSolar_(iGRU) + spectralCanopyAbsorbedSolar(iBand)

   end do  ! (looping through spectral bands)



  ! -----------------------------------------------------------------------------------------------------------------------------------------------------------
  ! method of Mahat and Tarboton (WRR, 2012)
  case(UEB_2stream)

   ! define transmission coefficient (-)
   scalarGproj = gProjParam
   transCoef   = scalarGproj/scalarCosZenith_(iGRU)

   ! define "k-prime" coefficient (-)
   transCoefPrime = sqrt(1._rkind - bScatParam)

   ! compute ground albedo (-)
   groundAlbedoDirect  = Frad_vis*spectralAlbGndDirect_(ixVisible,iGRU)  + (1._rkind - Frad_vis)*spectralAlbGndDirect_(ixNearIR,iGRU)
   groundAlbedoDiffuse = Frad_vis*spectralAlbGndDiffuse_(ixVisible,iGRU) + (1._rkind - Frad_vis)*spectralAlbGndDiffuse_(ixNearIR,iGRU)

   ! compute transmission for an infinite canopy (-)
   tauInfinite = exp(-transCoef*transCoefPrime*scalarExposedVAI)

   ! compute upward reflection factor for an infinite canopy (-)
   betaInfinite = (1._rkind - transCoefPrime)/(1._rkind + transCoefPrime)

   ! compute transmission for a finite canopy (-)
   tauFinite = tauInfinite*(1._rkind - betaInfinite**2_i4b)/(1._rkind - (betaInfinite**2_i4b)*tauInfinite**2_i4b)

   ! compute reflectance for a finite canopy (-)
   betaFinite = betaInfinite*(1._rkind - tauInfinite**2_i4b) / (1._rkind - (betaInfinite**2_i4b)*(tauInfinite**2_i4b))

   ! compute transmission of diffuse radiation (-)
   vFactor      = transCoefPrime*scalarGproj*scalarExposedVAI
   expi         = expInt(vFactor)
   taudInfinite = (1._rkind - vFactor)*exp(-vFactor) + (vFactor**2_i4b)*expi
   taudFinite   = taudInfinite*(1._rkind - betaInfinite**2_i4b)/(1._rkind - (betaInfinite**2_i4b)*taudInfinite**2_i4b)

   ! compute reflectance of diffuse radiation (-)
   betadFinite  = betaInfinite*(1._rkind - taudInfinite**2_i4b) / (1._rkind - (betaInfinite**2_i4b)*(taudInfinite**2_i4b))

   ! compute total transmission of direct and diffuse radiation, accounting for multiple reflections (-)
   refMult    = 1._rkind / (1._rkind - groundAlbedoDiffuse*betadFinite*(1._rkind - taudFinite) )


   tauDirect  = tauFinite*refMult
   tauDiffuse = taudFinite*refMult

   ! compute fraction of radiation lost to space (-)
   fractionRefDirect  = ( (1._rkind - groundAlbedoDirect)*betaFinite   + groundAlbedoDirect*tauFinite*taudFinite) * refMult
   fractionRefDiffuse = ( (1._rkind - groundAlbedoDiffuse)*betadFinite + groundAlbedoDiffuse*taudFinite*taudFinite) * refMult

   ! compute radiation in each spectral band (W m-2)
   do iBand=1,nBands

    ! compute below-canopy radiation (W m-2)
    spectralBelowCanopyDirect_(iBand,iGRU)  = spectralIncomingDirect_(iBand,iGRU)*tauFinite*refMult                ! direct radiation from current wave band
    spectralBelowCanopyDiffuse_(iBand,iGRU) = spectralIncomingDiffuse_(iBand,iGRU)*taudFinite*refMult              ! diffuse radiation from current wave band
    spectralBelowCanopySolar(iBand)   = spectralBelowCanopyDirect_(iBand,iGRU) + spectralBelowCanopyDiffuse_(iBand,iGRU)

    ! compute radiation absorbed by the ground in given wave band (W m-2)
    spectralGroundAbsorbedDirect(iBand)  = (1._rkind - groundAlbedoDirect)*spectralBelowCanopyDirect_(iBand,iGRU)
    spectralGroundAbsorbedDiffuse(iBand) = (1._rkind - groundAlbedoDiffuse)*spectralBelowCanopyDiffuse_(iBand,iGRU)
    spectralGroundAbsorbedSolar(iBand)   = spectralGroundAbsorbedDirect(iBand) + spectralGroundAbsorbedDiffuse(iBand)

    ! compute radiation absorbed by vegetation in current wave band (W m-2)
    spectralCanopyAbsorbedDirect(iBand)  = spectralIncomingDirect_(iBand,iGRU)*(1._rkind - tauFinite)*(1._rkind - betaFinite) + &    ! (radiation absorbed on the way down)
                                           spectralBelowCanopyDirect_(iBand,iGRU)*groundAlbedoDirect*(1._rkind - taudFinite)      ! (radiation absorbed on the way up)
    spectralCanopyAbsorbedDiffuse(iBand) = spectralIncomingDiffuse_(iBand,iGRU)*(1._rkind - taudFinite)*(1._rkind - betadFinite) + & ! (radiation absorbed on the way down)
                                           spectralBelowCanopyDiffuse_(iBand,iGRU)*groundAlbedoDiffuse*(1._rkind - taudFinite)    ! (radiation absorbed on the way up)
    spectralCanopyAbsorbedSolar(iBand)   = spectralCanopyAbsorbedDirect(iBand) + spectralCanopyAbsorbedDiffuse(iBand)

    ! compute solar radiation lost to space in given wave band (W m-2)
    spectralTotalReflectedDirect(iBand)  = spectralIncomingDirect_(iBand,iGRU) - spectralGroundAbsorbedDirect(iBand) - spectralCanopyAbsorbedDirect(iBand)
    spectralTotalReflectedDiffuse(iBand) = spectralIncomingDiffuse_(iBand,iGRU) - spectralGroundAbsorbedDiffuse(iBand) - spectralCanopyAbsorbedDiffuse(iBand)
    spectralTotalReflectedSolar(iBand)   = spectralTotalReflectedDirect(iBand) + spectralTotalReflectedDiffuse(iBand)
    if(spectralTotalReflectedDirect(iBand) < 0._rkind .or. spectralTotalReflectedDiffuse(iBand) < 0._rkind)then
     err=20; return
    end if

    ! save canopy radiation absorbed in visible wavelengths
    if(iBand == ixVisible)then
     visibleAbsDirect  = spectralCanopyAbsorbedDirect(ixVisible)
     visibleAbsDiffuse = spectralCanopyAbsorbedDiffuse(ixVisible)
    end if

    ! accumulate fluxes
    scalarBelowCanopySolar_(iGRU)    = scalarBelowCanopySolar_(iGRU) + spectralBelowCanopySolar(iBand)
    scalarGroundAbsorbedSolar_(iGRU) = scalarGroundAbsorbedSolar_(iGRU) + spectralGroundAbsorbedSolar(iBand)
    scalarCanopyAbsorbedSolar_(iGRU) = scalarCanopyAbsorbedSolar_(iGRU) + spectralCanopyAbsorbedSolar(iBand)

   end do  ! (looping through wave bands)



  ! -----------------------------------------------------------------------------------------------------------------------------------------------------------
  ! CLM approach
  case(CLM_2stream)

   ! weight reflectance and transmittance by exposed leaf and stem area index
   weightLeaf       = scalarExposedLAI_(iGRU) / scalarExposedVAI
   weightStem       = scalarExposedSAI_(iGRU) / scalarExposedVAI
   do iBand = 1,nBands  ! loop through spectral bands
    spectralVegReflc(iBand) = RHOL(vegTypeIndex_(iGRU),iBand)*weightLeaf + RHOS(vegTypeIndex_(iGRU),iBand)*weightStem
    spectralVegTrans(iBand) = TAUL(vegTypeIndex_(iGRU),iBand)*weightLeaf + TAUS(vegTypeIndex_(iGRU),iBand)*weightStem
   end do

   ! loop through wave bands
   do iBand=1,nBands

    ic = 0
    ! two-stream approximation for direct-beam radiation (from CLM/Noah-MP)
    call twoStream_d(&
                   ! input
                   iBand,                             & ! intent(in): waveband number
                   ic,                                & ! intent(in): 0=unit incoming direct; 1=unit incoming diffuse
                   vegTypeIndex_(iGRU),                      & ! intent(in): vegetation type
                   scalarCosZenith_(iGRU),                   & ! intent(in): cosine of direct zenith angle (0-1)
                   scalarExposedVAI,                  & ! intent(in): one-sided leaf+stem area index (m2/m2)
                   scalarCanopyWetFraction_(iGRU),           & ! intent(in): fraction of lai, sai that is wetted (-)
                   scalarCanopyTempTrial_(iGRU),             & ! intent(in): surface temperature (k)
                   spectralAlbGndDirect_(:,iGRU),              & ! intent(in): direct  albedo of underlying surface (1:nBands) (-)
                   spectralAlbGndDiffuse_(:,iGRU),             & ! intent(in): diffuse albedo of underlying surface (1:nBands) (-)
                   spectralVegReflc,                  & ! intent(in): leaf+stem reflectance (1:nbands)
                   spectralVegTrans,                  & ! intent(in): leaf+stem transmittance (1:nBands)
                   scalarVegFraction,                 & ! intent(in): vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
                   ist,                               & ! intent(in): surface type
                   iLoc,jLoc,                         & ! intent(in): grid indices
                   ! output
                   spectralCanopyAbsorbedDirect,      & ! intent(out): flux abs by veg layer (per unit incoming flux), (1:nBands)
                   spectralTotalReflectedDirect,      & ! intent(out): flux refl above veg layer (per unit incoming flux), (1:nBands)
                   spectralDirectBelowCanopyDirect,   & ! intent(out): down dir flux below veg layer (per unit in flux), (1:nBands)
                   spectralDiffuseBelowCanopyDirect,  & ! intent(out): down dif flux below veg layer (per unit in flux), (1:nBands)
                   scalarGproj,                       & ! intent(out): projected leaf+stem area in solar direction
                   spectralCanopyReflectedDirect,     & ! intent(out): flux reflected by veg layer   (per unit incoming flux), (1:nBands)
                   spectralGroundReflectedDirect,     & ! intent(out): flux reflected by ground (per unit incoming flux), (1:nBands)
                   ! input-output
                   scalarBetweenCanopyGapFraction,    & ! intent(inout): between canopy gap fraction for beam (-)
                   scalarWithinCanopyGapFraction,omegas,opt_rad,tfrz,betais,betads,xl,hvt,hvb,rc      ) ! intent(inout): within canopy gap fraction for beam (-)

    ic = 1
    ! two-stream approximation for diffuse radiation (from CLM/Noah-MP)
    call twoStream_d(&
                   ! input
                   iBand,                             & ! intent(in): waveband number
                   ic,                                & ! intent(in): 0=unit incoming direct; 1=unit incoming diffuse
                   vegTypeIndex_(iGRU),                      & ! intent(in): vegetation type
                   scalarCosZenith_(iGRU),                   & ! intent(in): cosine of direct zenith angle (0-1)
                   scalarExposedVAI,                  & ! intent(in): one-sided leaf+stem area index (m2/m2)
                   scalarCanopyWetFraction_(iGRU),           & ! intent(in): fraction of lai, sai that is wetted (-)
                   scalarCanopyTempTrial_(iGRU),             & ! intent(in): surface temperature (k)
                   spectralAlbGndDirect_(:,iGRU),              & ! intent(in): direct  albedo of underlying surface (1:nBands) (-)
                   spectralAlbGndDiffuse_(:,iGRU),             & ! intent(in): diffuse albedo of underlying surface (1:nBands) (-)
                   spectralVegReflc,                  & ! intent(in): leaf+stem reflectance (1:nbands)
                   spectralVegTrans,                  & ! intent(in): leaf+stem transmittance (1:nBands)
                   scalarVegFraction,                 & ! intent(in): vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
                   ist,                               & ! intent(in): surface type
                   iLoc,jLoc,                         & ! intent(in): grid indices
                   ! output
                   spectralCanopyAbsorbedDiffuse,     & ! intent(out): flux abs by veg layer (per unit incoming flux), (1:nBands)
                   spectralTotalReflectedDiffuse,     & ! intent(out): flux refl above veg layer (per unit incoming flux), (1:nBands)
                   spectralDirectBelowCanopyDiffuse,  & ! intent(out): down dir flux below veg layer (per unit in flux), (1:nBands)
                   spectralDiffuseBelowCanopyDiffuse, & ! intent(out): down dif flux below veg layer (per unit in flux), (1:nBands)
                   scalarGproj,                       & ! intent(out): projected leaf+stem area in solar direction
                   spectralCanopyReflectedDiffuse,    & ! intent(out): flux reflected by veg layer   (per unit incoming flux), (1:nBands)
                   spectralGroundReflectedDiffuse,    & ! intent(out): flux reflected by ground (per unit incoming flux), (1:nBands)
                   ! input-output
                   scalarBetweenCanopyGapFraction,    & ! intent(inout): between canopy gap fraction for beam (-)
                   scalarWithinCanopyGapFraction,omegas,opt_rad,tfrz,betais,betads,xl,hvt,hvb,rc      ) ! intent(inout): within canopy gap fraction for beam (-)


    ! compute below-canopy radiation
    spectralBelowCanopyDirect_(iBand,iGRU)  = spectralIncomingDirect_(iBand,iGRU)*spectralDirectBelowCanopyDirect(iBand)      ! direct radiation
    spectralBelowCanopyDiffuse_(iBand,iGRU) = spectralIncomingDirect_(iBand,iGRU)*spectralDiffuseBelowCanopyDirect(iBand) + & ! direct radiation transmitted as diffuse
                                        spectralIncomingDiffuse_(iBand,iGRU)*spectralDiffuseBelowCanopyDiffuse(iBand)   ! diffuse radiation transmitted as diffuse

    ! accumulate radiation transmitted below the canopy (W m-2)
    scalarBelowCanopySolar_(iGRU)    = scalarBelowCanopySolar_(iGRU) + &                                                  ! contribution from all previous wave bands
                                spectralBelowCanopyDirect_(iBand,iGRU) + spectralBelowCanopyDiffuse_(iBand,iGRU)        ! contribution from current wave band

    ! accumulate radiation absorbed by the vegetation canopy (W m-2)
    scalarCanopyAbsorbedSolar_(iGRU) = scalarCanopyAbsorbedSolar_(iGRU) + &                                               ! contribution from all previous wave bands
                                spectralIncomingDirect_(iBand,iGRU)*spectralCanopyAbsorbedDirect(iBand) + &       ! direct radiation from current wave band
                                spectralIncomingDiffuse_(iBand,iGRU)*spectralCanopyAbsorbedDiffuse(iBand)         ! diffuse radiation from current wave band

    ! accumulate radiation absorbed by the ground (W m-2)
    scalarGroundAbsorbedSolar_(iGRU) = scalarGroundAbsorbedSolar_(iGRU) + &                                               ! contribution from all previous wave bands
                                spectralBelowCanopyDirect_(iBand,iGRU)*(1._rkind - spectralAlbGndDirect_(iBand,iGRU)) + &  ! direct radiation from current wave band
                                spectralBelowCanopyDiffuse_(iBand,iGRU)*(1._rkind - spectralAlbGndDiffuse_(iBand,iGRU))    ! diffuse radiation from current wave band

    ! save canopy radiation absorbed in visible wavelengths
    ! NOTE: here flux is per unit incoming flux
    if(iBand == ixVisible)then
     visibleAbsDirect  = spectralIncomingDirect_(ixVisible,iGRU)*spectralCanopyAbsorbedDirect(ixVisible)
     visibleAbsDiffuse = spectralIncomingDiffuse_(ixVisible,iGRU)*spectralCanopyAbsorbedDiffuse(ixVisible)
    end if

   end do  ! (looping through wave bands)

  ! -----------------------------------------------------------------------------------------------------------------------------------------------------------
  case default; err=20;  return

 end select ! (option for canopy sw radiation)


 ! ============================================================================================================================================================
 ! ============================================================================================================================================================

 ! compute variables used in photosynthesis routines

 ! compute sunlit fraction of canopy (from CLM/Noah-MP)
 ext = scalarGproj/scalarCosZenith_(iGRU)  ! optical depth of direct beam per unit leaf + stem area
 scalarCanopySunlitFraction_(iGRU) = (1._rkind - exp(-ext*scalarExposedVAI)) / max(ext*scalarExposedVAI,mpe)
 if(scalarCanopySunlitFraction_(iGRU) < 0.01_rkind) scalarCanopySunlitFraction_(iGRU) = 0._rkind

 ! compute sunlit and shaded LAI
 scalarCanopyShadedFraction = 1._rkind - scalarCanopySunlitFraction_(iGRU)
 scalarCanopySunlitLAI_(iGRU)      = scalarExposedLAI_(iGRU)*scalarCanopySunlitFraction_(iGRU)
 scalarCanopyShadedLAI_(iGRU)      = scalarExposedLAI_(iGRU)*scalarCanopyShadedFraction

 ! compute PAR for sunlit and shaded leaves (from CLM/Noah-MP)
 fractionLAI       = scalarExposedLAI_(iGRU) / max(scalarExposedVAI, mpe)
 if(scalarCanopySunlitFraction_(iGRU) > tiny(scalarCanopySunlitFraction_(iGRU)))then
  scalarCanopySunlitPAR_(iGRU) = (visibleAbsDirect + scalarCanopySunlitFraction_(iGRU)*visibleAbsDiffuse) * fractionLAI / max(scalarCanopySunlitLAI_(iGRU), mpe)
  scalarCanopyShadedPAR_(iGRU) = (                   scalarCanopyShadedFraction*visibleAbsDiffuse) * fractionLAI / max(scalarCanopyShadedLAI_(iGRU), mpe)
 else
  scalarCanopySunlitPAR_(iGRU) = 0._rkind
  scalarCanopyShadedPAR_(iGRU) = (visibleAbsDirect + visibleAbsDiffuse) * fractionLAI / max(scalarCanopyShadedLAI_(iGRU), mpe)
 end if



 end subroutine canopy_SW


 ! *************************************************************************************************************************************
 ! private subroutine gndAlbedo: compute the albedo of the ground surface
 ! *************************************************************************************************************************************
 attributes(device) subroutine gndAlbedo(&
                      ! input
                      isc,                                   & ! intent(in): index of soil color
                      scalarGroundSnowFraction,              & ! intent(in): fraction of ground that is snow covered (-)
                      scalarVolFracLiqUpper,                 & ! intent(in): volumetric liquid water content in upper-most soil layer (-)
                      spectralSnowAlbedoDirect,              & ! intent(in): direct albedo of snow in each spectral band (-)
                      spectralSnowAlbedoDiffuse,             & ! intent(in): diffuse albedo of snow in each spectral band (-)
                      ! output
                      spectralAlbGndDirect,                  & ! intent(out): direct  albedo of underlying surface (-)
                      spectralAlbGndDiffuse,                 & ! intent(out): diffuse albedo of underlying surface (-)
                      err,albsat,albdry)                             ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------------------------------------
 ! identify parameters for soil albedo
!  USE NOAHMP_RAD_PARAMETERS, only: ALBSAT,ALBDRY  ! Noah-MP: saturated and dry soil albedos for each wave band
 ! --------------------------------------------------------------------------------------------------------------------------------------
                      implicit none
 ! input: model control
 integer(i4b),intent(in)        :: isc                          ! index of soil color
 real(rkind),intent(in)            :: scalarGroundSnowFraction     ! fraction of ground that is snow covered (-)
 real(rkind),intent(in)            :: scalarVolFracLiqUpper        ! volumetric liquid water content in upper-most soil layer (-)
 real(rkind),intent(in)            :: spectralSnowAlbedoDirect(:)  ! direct albedo of snow in each spectral band (-)
 real(rkind),intent(in)            :: spectralSnowAlbedoDiffuse(:) ! diffuse albedo of snow in each spectral band (-)
 ! output
 real(rkind),intent(out)           :: spectralAlbGndDirect(:)      ! direct  albedo of underlying surface (-)
 real(rkind),intent(out)           :: spectralAlbGndDiffuse(:)     ! diffuse albedo of underlying surface (-)
 integer(i4b),intent(out)       :: err                          ! error code
 real(rkind) :: albsat(:,:),albdry(:,:)
 ! local variables
 integer(i4b)                   :: iBand                        ! index of spectral band
 real(rkind)                       :: xInc                         ! soil water correction factor for soil albedo
 real(rkind),dimension(1:nBands)   :: spectralSoilAlbedo           ! soil albedo in each spectral band
 ! initialize error control
 err=0; 
 ! compute soil albedo
 do iBand=1,nBands   ! loop through spectral bands
  xInc = max(0.11_rkind - 0.40_rkind*scalarVolFracLiqUpper, 0._rkind)
  spectralSoilAlbedo(iBand)  = min(ALBSAT(isc,iBand)+xInc,ALBDRY(isc,iBand))
 end do  ! (looping through spectral bands)

 ! compute surface albedo (weighted combination of snow and soil)
 do iBand=1,nBands
  spectralAlbGndDirect(iBand)  = (1._rkind - scalarGroundSnowFraction)*spectralSoilAlbedo(iBand)  + scalarGroundSnowFraction*spectralSnowAlbedoDirect(iBand)
  spectralAlbGndDiffuse(iBand) = (1._rkind - scalarGroundSnowFraction)*spectralSoilAlbedo(iBand)  + scalarGroundSnowFraction*spectralSnowAlbedoDiffuse(iBand)
 end do  ! (looping through spectral bands)

 end subroutine gndAlbedo

attributes(global) SUBROUTINE RADIATION_d (nGRU, nSnow,VEGTYP_     ,NSOIL   , & !in
                        SNEQVO_  ,SNEQV_   ,DT      ,COSZ_    ,SNOWH_   , & !in
                        TG_      ,TV_      ,FSNO_    ,QSNOW_   ,FWET_    , & !in
                        ELAI_    ,ESAI_    ,SMC_     ,SOLAD_   ,SOLAI_   , & !in
                        FVEG    ,ILOC    ,JLOC    ,                   & !in
                        ALBOLD_  ,TAUSS_   ,                            & !inout
                        FSUN_    ,LAISUN_  ,LAISHA_  ,PARSUN_  ,PARSHA_  , & !out
                        SAV_     ,SAG_     ,FSR_     ,FSA_     ,FSRV_    , &
                        FSRG_    ,BGAP_    ,WGAP_, &
                        rhol, rhos,taul,taus, opt_alb,omegas,opt_rad,tfrz,betais,betads,xl,hvt,hvb,rc, swemx,albsat,albdry,alblak)            !out
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input
  INTEGER, INTENT(IN),value                  :: ILOC
  integer(i4b),value :: nGRU
  integer(i4b) :: nSnow(:)
  INTEGER, INTENT(IN),value                  :: JLOC
  INTEGER, INTENT(IN)                  :: VEGTYP_(:) !vegetation type
  INTEGER, INTENT(IN),value                  :: NSOIL  !number of soil layers

  REAL(rkind), INTENT(IN),value                     :: DT     !time step [s]
  REAL(rkind), INTENT(IN)                     :: QSNOW_(:)  !snowfall (mm/s)
  REAL(rkind), INTENT(IN)                     :: SNEQVO_(:) !snow mass at last time step(mm)
  REAL(rkind), INTENT(IN)                     :: SNEQV_(:)  !snow mass (mm)
  REAL(rkind), INTENT(IN)                     :: SNOWH_(:)  !snow height (mm)
  REAL(rkind), INTENT(IN)                     :: COSZ_(:)   !cosine solar zenith angle (0-1)
  REAL(rkind), INTENT(IN)                     :: TG_(:,:)     !ground temperature (k)
  REAL(rkind), INTENT(IN)                     :: TV_(:)     !vegetation temperature (k)
  REAL(rkind), INTENT(IN)                     :: ELAI_(:)   !LAI, one-sided, adjusted for burying by snow
  REAL(rkind), INTENT(IN)                     :: ESAI_(:)   !SAI, one-sided, adjusted for burying by snow
  REAL(rkind), INTENT(IN)                     :: FWET_(:)   !fraction of canopy that is wet
  REAL(rkind), INTENT(IN) :: SMC_(:,:)    !volumetric soil water [m3/m3]
  REAL(rkind)    , INTENT(IN) :: SOLAD_(:,:)  !incoming direct solar radiation (w/m2)
  REAL(rkind)    , INTENT(IN) :: SOLAI_(:,:)  !incoming diffuse solar radiation (w/m2)
  REAL(rkind), INTENT(IN)                     :: FSNO_(:)   !snow cover fraction (-)
  REAL(rkind), INTENT(IN)                     :: FVEG   !green vegetation fraction [0.0-1.0]

! inout
  REAL(rkind),                  INTENT(INOUT) :: ALBOLD_(:) !snow albedo at last time step (CLASS type)
  REAL(rkind),                  INTENT(INOUT) :: TAUSS_(:)  !non-dimensional snow age.

! output
  REAL(rkind), INTENT(OUT)                    :: FSUN_(:)   !sunlit fraction of canopy (-)
  REAL(rkind), INTENT(OUT)                    :: LAISUN_(:) !sunlit leaf area (-)
  REAL(rkind), INTENT(OUT)                    :: LAISHA_(:) !shaded leaf area (-)
  REAL(rkind), INTENT(OUT)                    :: PARSUN_(:) !average absorbed par for sunlit leaves (w/m2)
  REAL(rkind), INTENT(OUT)                    :: PARSHA_(:) !average absorbed par for shaded leaves (w/m2)
  REAL(rkind), INTENT(OUT)                    :: SAV_(:)    !solar radiation absorbed by vegetation (w/m2)
  REAL(rkind), INTENT(OUT)                    :: SAG_(:)    !solar radiation absorbed by ground (w/m2)
  REAL(rkind), INTENT(OUT)                    :: FSA_(:)    !total absorbed solar radiation (w/m2)
  REAL(rkind), INTENT(OUT)                    :: FSR_(:)    !total reflected solar radiation (w/m2)

!jref:start
  REAL(rkind), INTENT(OUT)                    :: FSRV_(:)    !veg. reflected solar radiation (w/m2)
  REAL(rkind), INTENT(OUT)                    :: FSRG_(:)    !ground reflected solar radiation (w/m2)
  REAL(rkind), INTENT(OUT)                    :: BGAP_(:)
  REAL(rkind), INTENT(OUT)                    :: WGAP_(:)
!jref:end
  real(rkind) :: rhol(:,:), rhos(:,:),taul(:,:),taus(:,:)
  integer(i4b) :: opt_alb
  real(rkind) :: omegas(:)
integer(i4b) :: opt_rad
real(rkind) :: tfrz,betais,betads
real(rkind) :: xl(:),hvt(:),hvb(:),rc(:)
real(rkind) :: swemx
real(rkind) :: albsat(:,:),albdry(:,:),alblak(:)

! local
  REAL(rkind)                                 :: FAGE   !snow age function (0 - new snow)
  REAL(rkind), DIMENSION(1:2)                 :: ALBGRD !ground albedo (direct)
  REAL(rkind), DIMENSION(1:2)                 :: ALBGRI !ground albedo (diffuse)
  REAL(rkind), DIMENSION(1:2)                 :: ALBD   !surface albedo (direct)
  REAL(rkind), DIMENSION(1:2)                 :: ALBI   !surface albedo (diffuse)
  REAL(rkind), DIMENSION(1:2)                 :: FABD   !flux abs by veg (per unit direct flux)
  REAL(rkind), DIMENSION(1:2)                 :: FABI   !flux abs by veg (per unit diffuse flux)
  REAL(rkind), DIMENSION(1:2)                 :: FTDD   !down direct flux below veg (per unit dir flux)
  REAL(rkind), DIMENSION(1:2)                 :: FTID   !down diffuse flux below veg (per unit dir flux)
  REAL(rkind), DIMENSION(1:2)                 :: FTII   !down diffuse flux below veg (per unit dif flux)
!jref:start
  REAL(rkind), DIMENSION(1:2)                 :: FREVI
  REAL(rkind), DIMENSION(1:2)                 :: FREVD
  REAL(rkind), DIMENSION(1:2)                 :: FREGI
  REAL(rkind), DIMENSION(1:2)                 :: FREGD
!jref:end

  REAL(rkind)                                 :: FSHA   !shaded fraction of canopy
  REAL(rkind)                                 :: VAI    !total LAI + stem area index, one sided

  REAL(rkind),PARAMETER :: MPE = 1.E-6
  LOGICAL VEG  !true: vegetated for surface temperature calculation
  integer(i4b) :: iGRU
! dimensions
!   integer(i4b),parameter        :: nBands=2      ! number of spectral bands for shortwave radiation
  ! named variables
!   integer(i4b),parameter        :: ist     = 1   ! Surface type:  IST=1 => soil;  IST=2 => lake
!   integer(i4b),parameter        :: isc     = 4   ! Soil color type
  ! spatial indices
!   integer(i4b),parameter        :: iLoc    = 1   ! i-location
!   integer(i4b),parameter        :: jLoc    = 1   ! j-location
  ! algorithmic parameters
!   real(rkind),parameter            :: missingValue=-9999._rkind  ! missing value, used when diagnostic or state variables are undefined
!   real(rkind),parameter            :: verySmall=1.e-6_rkind   ! used as an additive constant to check if substantial difference among real numbers
!   real(rkind),parameter            :: mpe=1.e-6_rkind         ! prevents overflow error if division by zero
!   real(rkind),parameter            :: dx=1.e-6_rkind          ! finite difference increment
  iGRU = (blockidx%x-1) * blockdim%x + threadidx%x
  if (iGRU .gt. nGRU) return

! --------------------------------------------------------------------------------------------------

! surface abeldo

   CALL ALBEDO_d (VEGTYP_(iGRU) ,IST    ,ISC    ,ICE    ,NSOIL  , & !in
                DT     ,COSZ_(iGRU)   ,FAGE   ,ELAI_(iGRU)   ,ESAI_(iGRU)   , & !in
                TG_(1,iGRU)     ,TV_(iGRU)     ,SNOWH_(iGRU)*1000._rkind,FSNO_(iGRU)   ,FWET_(iGRU)   , & !in
                SMC_(nSnow(iGRU)+1:,iGRU)    ,SNEQVO_(iGRU) ,SNEQV_(iGRU)  ,QSNOW_(iGRU)  ,FVEG   , & !in
                ILOC   ,JLOC   ,                         & !in
                ALBOLD_(iGRU) ,TAUSS_(iGRU)                          , & !inout
                ALBGRD ,ALBGRI ,ALBD   ,ALBI   ,FABD   , & !out
                FABI   ,FTDD   ,FTID   ,FTII   ,FSUN_(iGRU)   , & !)   !out
                FREVI  ,FREVD   ,FREGD ,FREGI  ,BGAP_(iGRU)   , & !inout
                WGAP_(iGRU),rhol,rhos,taul,taus,opt_alb,omegas,opt_rad,tfrz,betais,betads,xl,hvt,hvb,rc,swemx,albsat,albdry,alblak)

! surface radiation

     FSHA = 1.-FSUN_(iGRU)
     LAISUN_(iGRU) = ELAI_(iGRU)*FSUN_(iGRU)
     LAISHA_(iGRU) = ELAI_(iGRU)*FSHA
     VAI = ELAI_(iGRU)+ ESAI_(iGRU)
     IF (VAI .GT. 0.) THEN
        VEG = .TRUE.
     ELSE
        VEG = .FALSE.
     END IF

   CALL SURRAD_d (MPE    ,FSUN_(iGRU)   ,FSHA   ,ELAI_(iGRU)   ,VAI    , & !in
                LAISUN_(iGRU) ,LAISHA_(iGRU) ,SOLAD_(:,iGRU)  ,SOLAI_(:,iGRU)  ,FABD   , & !in
                FABI   ,FTDD   ,FTID   ,FTII   ,ALBGRD , & !in
                ALBGRI ,ALBD   ,ALBI   ,ILOC   ,JLOC   , & !in
                PARSUN_(iGRU) ,PARSHA_(iGRU) ,SAV_(iGRU)    ,SAG_(iGRU)    ,FSA_(iGRU)    , & !out
                FSR_(iGRU)    ,                                 & !out
                FREVI  ,FREVD  ,FREGD  ,FREGI  ,FSRV_(iGRU)   , & !inout
                FSRG_(iGRU))

  END SUBROUTINE RADIATION_d

 attributes(device) SUBROUTINE TWOSTREAM_d (IB     ,IC      ,VEGTYP  ,COSZ    ,VAI    , & !in
    FWET   ,T       ,ALBGRD  ,ALBGRI  ,RHO    , & !in
    TAU    ,FVEG    ,IST     ,ILOC    ,JLOC   , & !in
    FAB    ,FRE     ,FTD     ,FTI     ,GDIR   , & !)   !out
    FREV   ,FREG    ,BGAP    ,WGAP,omegas,opt_rad,tfrz,betais,betads,xl,hvt,hvb,rc)

! --------------------------------------------------------------------------------------------------
! use two-stream approximation of Dickinson (1983) Adv Geophysics
! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
! to calculate fluxes absorbed by vegetation, reflected by vegetation,
! and transmitted through vegetation for unit incoming direct or diffuse
! flux given an underlying surface with known albedo.
! --------------------------------------------------------------------------------------------------
! USE NOAHMP_VEG_PARAMETERS
! USE NOAHMP_RAD_PARAMETERS
! --------------------------------------------------------------------------------------------------
IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input

INTEGER,              INTENT(IN)  :: ILOC    !grid index
INTEGER,              INTENT(IN)  :: JLOC    !grid index
INTEGER,              INTENT(IN)  :: IST     !surface type
INTEGER,              INTENT(IN)  :: IB      !waveband number
INTEGER,              INTENT(IN)  :: IC      !0=unit incoming direct; 1=unit incoming diffuse
INTEGER,              INTENT(IN)  :: VEGTYP  !vegetation type

REAL(rkind),                 INTENT(IN)  :: COSZ    !cosine of direct zenith angle (0-1)
REAL(rkind),                 INTENT(IN)  :: VAI     !one-sided leaf+stem area index (m2/m2)
REAL(rkind),                 INTENT(IN)  :: FWET    !fraction of lai, sai that is wetted (-)
REAL(rkind),                 INTENT(IN)  :: T       !surface temperature (k)

REAL(rkind), DIMENSION(1:2), INTENT(IN)  :: ALBGRD  !direct  albedo of underlying surface (-)
REAL(rkind), DIMENSION(1:2), INTENT(IN)  :: ALBGRI  !diffuse albedo of underlying surface (-)
REAL(rkind), DIMENSION(1:2), INTENT(IN)  :: RHO     !leaf+stem reflectance
REAL(rkind), DIMENSION(1:2), INTENT(IN)  :: TAU     !leaf+stem transmittance
REAL(rkind),                 INTENT(IN)  :: FVEG    !green vegetation fraction [0.0-1.0]

! output

REAL(rkind), DIMENSION(1:2), INTENT(OUT) :: FAB     !flux abs by veg layer (per unit incoming flux)
REAL(rkind), DIMENSION(1:2), INTENT(OUT) :: FRE     !flux refl above veg layer (per unit incoming flux)
REAL(rkind), DIMENSION(1:2), INTENT(OUT) :: FTD     !down dir flux below veg layer (per unit in flux)
REAL(rkind), DIMENSION(1:2), INTENT(OUT) :: FTI     !down dif flux below veg layer (per unit in flux)
REAL(rkind),                 INTENT(OUT) :: GDIR    !projected leaf+stem area in solar direction
REAL(rkind), DIMENSION(1:2), INTENT(OUT) :: FREV    !flux reflected by veg layer   (per unit incoming flux)
REAL(rkind), DIMENSION(1:2), INTENT(OUT) :: FREG    !flux reflected by ground (per unit incoming flux)
real(rkind) :: omegas(:)
integer(i4b) :: opt_rad
real(rkind) :: tfrz,betais,betads
real(rkind) :: xl(:),hvt(:),hvb(:),rc(:)
! local
REAL(rkind)                              :: OMEGA   !fraction of intercepted radiation that is scattered
REAL(rkind)                              :: OMEGAL  !omega for leaves
REAL(rkind)                              :: BETAI   !upscatter parameter for diffuse radiation
REAL(rkind)                              :: BETAIL  !betai for leaves
REAL(rkind)                              :: BETAD   !upscatter parameter for direct beam radiation
REAL(rkind)                              :: BETADL  !betad for leaves
REAL(rkind)                              :: EXT     !optical depth of direct beam per unit leaf area
REAL(rkind)                              :: AVMU    !average diffuse optical depth

REAL(rkind)                              :: COSZI   !0.001 <= cosz <= 1.000
REAL(rkind)                              :: ASU     !single scattering albedo
REAL(rkind)                              :: CHIL    ! -0.4 <= xl <= 0.6

REAL(rkind)                              :: TMP0,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,TMP7,TMP8,TMP9
REAL(rkind)                              :: P1,P2,P3,P4,S1,S2,U1,U2,U3
REAL(rkind)                              :: B,C,D,D1,D2,F,H,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10
REAL(rkind)                              :: PHI1,PHI2,SIGMA
REAL(rkind)                              :: FTDS,FTIS,FRES
REAL(rkind)                              :: DENFVEG
REAL(rkind)                              :: VAI_SPREAD
!jref:start
REAL(rkind)                              :: FREVEG,FREBAR,FTDVEG,FTIVEG,FTDBAR,FTIBAR
REAL(rkind)                              :: THETAZ
!jref:end

!  variables for the modified two-stream scheme
!  Niu and Yang (2004), JGR

REAL(rkind), PARAMETER :: PAI = 3.14159265
REAL(rkind) :: HD       !crown depth (m)
REAL(rkind) :: BB       !vertical crown radius (m)
REAL(rkind) :: THETAP   !angle conversion from SZA
REAL(rkind) :: FA       !foliage volume density (m-1)
REAL(rkind) :: NEWVAI   !effective LSAI (-)

REAL(rkind),INTENT(INOUT) :: BGAP     !between canopy gap fraction for beam (-)
REAL(rkind),INTENT(INOUT) :: WGAP     !within canopy gap fraction for beam (-)

REAL(rkind) :: KOPEN    !gap fraction for diffue light (-)
REAL(rkind) :: GAP      !total gap fraction for beam ( <=1-shafac )

! -----------------------------------------------------------------
! compute within and between gaps
VAI_SPREAD = VAI
if(VAI == 0.0) THEN
GAP     = 1.0
KOPEN   = 1.0
ELSE
IF(OPT_RAD == 1) THEN
DENFVEG = -LOG(MAX(1.0-FVEG,0.01))/(PAI*RC(VEGTYP)**2)
HD      = HVT(VEGTYP) - HVB(VEGTYP)
BB      = 0.5 * HD
THETAP  = ATAN(BB/RC(VEGTYP) * TAN(ACOS(MAX(0.01,COSZ))) )
! BGAP    = EXP(-DEN(VEGTYP) * PAI * RC(VEGTYP)**2/COS(THETAP) )
BGAP    = EXP(-DENFVEG * PAI * RC(VEGTYP)**2/COS(THETAP) )
FA      = VAI/(1.33 * PAI * RC(VEGTYP)**3.0 *(BB/RC(VEGTYP))*DENFVEG)
NEWVAI  = HD*FA
WGAP    = (1.0-BGAP) * EXP(-0.5*NEWVAI/COSZ)
GAP     = MIN(1.0-FVEG, BGAP+WGAP)

KOPEN   = 0.05
END IF

IF(OPT_RAD == 2) THEN
GAP     = 0.0
KOPEN   = 0.0
END IF

IF(OPT_RAD == 3) THEN
GAP     = 1.0-FVEG
KOPEN   = 1.0-FVEG
END IF
end if

! calculate two-stream parameters OMEGA, BETAD, BETAI, AVMU, GDIR, EXT.
! OMEGA, BETAD, BETAI are adjusted for snow. values for OMEGA*BETAD
! and OMEGA*BETAI are calculated and then divided by the new OMEGA
! because the product OMEGA*BETAI, OMEGA*BETAD is used in solution.
! also, the transmittances and reflectances (TAU, RHO) are linear
! weights of leaf and stem values.

COSZI  = MAX(0.001, COSZ)
CHIL   = MIN( MAX(XL(VEGTYP), -0.4), 0.6)
IF (ABS(CHIL) .LE. 0.01) CHIL = 0.01
PHI1   = 0.5 - 0.633*CHIL - 0.330*CHIL*CHIL
PHI2   = 0.877 * (1.-2.*PHI1)
GDIR   = PHI1 + PHI2*COSZI
EXT    = GDIR/COSZI
AVMU   = ( 1. - PHI1/PHI2 * LOG((PHI1+PHI2)/PHI1) ) / PHI2
OMEGAL = RHO(IB) + TAU(IB)
TMP0   = GDIR + PHI2*COSZI
TMP1   = PHI1*COSZI
ASU    = 0.5*OMEGAL*GDIR/TMP0 * ( 1.-TMP1/TMP0*LOG((TMP1+TMP0)/TMP1) )
BETADL = (1.+AVMU*EXT)/(OMEGAL*AVMU*EXT)*ASU
BETAIL = 0.5 * ( RHO(IB)+TAU(IB) + (RHO(IB)-TAU(IB))   &
* ((1.+CHIL)/2.)**2 ) / OMEGAL

! adjust omega, betad, and betai for intercepted snow

IF (T .GT. TFRZ) THEN                                !no snow
TMP0 = OMEGAL
TMP1 = BETADL
TMP2 = BETAIL
ELSE
TMP0 =   (1.-FWET)*OMEGAL        + FWET*OMEGAS(IB)
TMP1 = ( (1.-FWET)*OMEGAL*BETADL + FWET*OMEGAS(IB)*BETADS ) / TMP0
TMP2 = ( (1.-FWET)*OMEGAL*BETAIL + FWET*OMEGAS(IB)*BETAIS ) / TMP0
END IF

OMEGA = TMP0
BETAD = TMP1
BETAI = TMP2

! absorbed, reflected, transmitted fluxes per unit incoming radiation

B = 1. - OMEGA + OMEGA*BETAI
C = OMEGA*BETAI
TMP0 = AVMU*EXT
D = TMP0 * OMEGA*BETAD
F = TMP0 * OMEGA*(1.-BETAD)
TMP1 = B*B - C*C
H = SQRT(TMP1) / AVMU
SIGMA = TMP0*TMP0 - TMP1
if ( ABS (SIGMA) < 1.e-6 ) SIGMA = SIGN(1.e-6,REAL(SIGMA))
P1 = B + AVMU*H
P2 = B - AVMU*H
P3 = B + TMP0
P4 = B - TMP0
S1 = EXP(-H*VAI)
! MPC change: precision issues
S2 = max(epsilon(VAI),EXP(-EXT*VAI))
IF (IC .EQ. 0) THEN
U1 = B - C/ALBGRD(IB)
U2 = B - C*ALBGRD(IB)
U3 = F + C*ALBGRD(IB)
ELSE
U1 = B - C/ALBGRI(IB)
U2 = B - C*ALBGRI(IB)
U3 = F + C*ALBGRI(IB)
END IF
TMP2 = U1 - AVMU*H
TMP3 = U1 + AVMU*H
D1 = P1*TMP2/S1 - P2*TMP3*S1
TMP4 = U2 + AVMU*H
TMP5 = U2 - AVMU*H
D2 = TMP4/S1 - TMP5*S1
H1 = -D*P4 - C*F
TMP6 = D - H1*P3/SIGMA
TMP7 = ( D - C - H1/SIGMA*(U1+TMP0) ) * S2
H2 = ( TMP6*TMP2/S1 - P2*TMP7 ) / D1
H3 = - ( TMP6*TMP3*S1 - P1*TMP7 ) / D1
H4 = -F*P3 - C*D
TMP8 = H4/SIGMA
TMP9 = ( U3 - TMP8*(U2-TMP0) ) * S2
H5 = - ( TMP8*TMP4/S1 + TMP9 ) / D2
H6 = ( TMP8*TMP5*S1 + TMP9 ) / D2
H7 = (C*TMP2) / (D1*S1)
H8 = (-C*TMP3*S1) / D1
H9 = TMP4 / (D2*S1)
H10 = (-TMP5*S1) / D2

! downward direct and diffuse fluxes below vegetation
! Niu and Yang (2004), JGR.

IF (IC .EQ. 0) THEN
FTDS = S2                           *(1.0-GAP) + GAP
FTIS = (H4*S2/SIGMA + H5*S1 + H6/S1)*(1.0-GAP)
ELSE
FTDS = 0.
FTIS = (H9*S1 + H10/S1)*(1.0-KOPEN) + KOPEN
END IF
FTD(IB) = FTDS
FTI(IB) = FTIS

! flux reflected by the surface (veg. and ground)

IF (IC .EQ. 0) THEN
FRES   = (H1/SIGMA + H2 + H3)*(1.0-GAP  ) + ALBGRD(IB)*GAP
FREVEG = (H1/SIGMA + H2 + H3)*(1.0-GAP  )
FREBAR = ALBGRD(IB)*GAP                   !jref - separate veg. and ground reflection
ELSE
FRES   = (H7 + H8) *(1.0-KOPEN) + ALBGRI(IB)*KOPEN
FREVEG = (H7 + H8) *(1.0-KOPEN) + ALBGRI(IB)*KOPEN
FREBAR = 0                                !jref - separate veg. and ground reflection
END IF
FRE(IB) = FRES

FREV(IB) = FREVEG
FREG(IB) = FREBAR

! flux absorbed by vegetation

FAB(IB) = 1. - FRE(IB) - (1.-ALBGRD(IB))*FTD(IB) &
        - (1.-ALBGRI(IB))*FTI(IB)

END SUBROUTINE TWOSTREAM_d

attributes(device) SUBROUTINE SURRAD_d (MPE     ,FSUN    ,FSHA    ,ELAI    ,VAI     , & !in
                     LAISUN  ,LAISHA  ,SOLAD   ,SOLAI   ,FABD    , & !in
                     FABI    ,FTDD    ,FTID    ,FTII    ,ALBGRD  , & !in
                     ALBGRI  ,ALBD    ,ALBI    ,ILOC    ,JLOC    , & !in
                     PARSUN  ,PARSHA  ,SAV     ,SAG     ,FSA     , & !out
                     FSR     , & !)                                       !out
                     FREVI   ,FREVD   ,FREGD   ,FREGI   ,FSRV    , &
                     FSRG) !inout

! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input

  INTEGER, INTENT(IN)              :: ILOC
  INTEGER, INTENT(IN)              :: JLOC
  REAL(rkind), INTENT(IN)                 :: MPE     !prevents underflow errors if division by zero

  REAL(rkind), INTENT(IN)                 :: FSUN    !sunlit fraction of canopy
  REAL(rkind), INTENT(IN)                 :: FSHA    !shaded fraction of canopy
  REAL(rkind), INTENT(IN)                 :: ELAI    !leaf area, one-sided
  REAL(rkind), INTENT(IN)                 :: VAI     !leaf + stem area, one-sided
  REAL(rkind), INTENT(IN)                 :: LAISUN  !sunlit leaf area index, one-sided
  REAL(rkind), INTENT(IN)                 :: LAISHA  !shaded leaf area index, one-sided

  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: SOLAD   !incoming direct solar radiation (w/m2)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: SOLAI   !incoming diffuse solar radiation (w/m2)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FABD    !flux abs by veg (per unit incoming direct flux)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FABI    !flux abs by veg (per unit incoming diffuse flux)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FTDD    !down dir flux below veg (per incoming dir flux)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FTID    !down dif flux below veg (per incoming dir flux)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FTII    !down dif flux below veg (per incoming dif flux)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: ALBGRD  !ground albedo (direct)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: ALBGRI  !ground albedo (diffuse)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: ALBD    !overall surface albedo (direct)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: ALBI    !overall surface albedo (diffuse)

  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FREVD    !overall surface albedo veg (direct)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FREVI    !overall surface albedo veg (diffuse)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FREGD    !overall surface albedo grd (direct)
  REAL(rkind), DIMENSION(1:2), INTENT(IN) :: FREGI    !overall surface albedo grd (diffuse)

! output

  REAL(rkind), INTENT(OUT)                :: PARSUN  !average absorbed par for sunlit leaves (w/m2)
  REAL(rkind), INTENT(OUT)                :: PARSHA  !average absorbed par for shaded leaves (w/m2)
  REAL(rkind), INTENT(OUT)                :: SAV     !solar radiation absorbed by vegetation (w/m2)
  REAL(rkind), INTENT(OUT)                :: SAG     !solar radiation absorbed by ground (w/m2)
  REAL(rkind), INTENT(OUT)                :: FSA     !total absorbed solar radiation (w/m2)
  REAL(rkind), INTENT(OUT)                :: FSR     !total reflected solar radiation (w/m2)
  REAL(rkind), INTENT(OUT)                :: FSRV    !reflected solar radiation by vegetation
  REAL(rkind), INTENT(OUT)                :: FSRG    !reflected solar radiation by ground

! ------------------------ local variables ----------------------------------------------------
  INTEGER                          :: IB      !waveband number (1=vis, 2=nir)
  INTEGER                          :: NBAND   !number of solar radiation waveband classes

  REAL(rkind)                             :: ABS     !absorbed solar radiation (w/m2)
  REAL(rkind)                             :: RNIR    !reflected solar radiation [nir] (w/m2)
  REAL(rkind)                             :: RVIS    !reflected solar radiation [vis] (w/m2)
  REAL(rkind)                             :: LAIFRA  !leaf area fraction of canopy
  REAL(rkind)                             :: TRD     !transmitted solar radiation: direct (w/m2)
  REAL(rkind)                             :: TRI     !transmitted solar radiation: diffuse (w/m2)
  REAL(rkind), DIMENSION(1:2)             :: CAD     !direct beam absorbed by canopy (w/m2)
  REAL(rkind), DIMENSION(1:2)             :: CAI     !diffuse radiation absorbed by canopy (w/m2)
! ---------------------------------------------------------------------------------------------
   NBAND = 2

! zero summed solar fluxes

    SAG = 0.
    SAV = 0.
    FSA = 0.

! loop over nband wavebands

  DO IB = 1, NBAND

! absorbed by canopy

    CAD(IB) = SOLAD(IB)*FABD(IB)
    CAI(IB) = SOLAI(IB)*FABI(IB)
    SAV     = SAV + CAD(IB) + CAI(IB)
    FSA     = FSA + CAD(IB) + CAI(IB)

! transmitted solar fluxes incident on ground

    TRD = SOLAD(IB)*FTDD(IB)
    TRI = SOLAD(IB)*FTID(IB) + SOLAI(IB)*FTII(IB)

! solar radiation absorbed by ground surface

    ABS = TRD*(1.-ALBGRD(IB)) + TRI*(1.-ALBGRI(IB))
    SAG = SAG + ABS
    FSA = FSA + ABS
  END DO

! partition visible canopy absorption to sunlit and shaded fractions
! to get average absorbed par for sunlit and shaded leaves

     LAIFRA = ELAI / MAX(VAI,MPE)
     IF (FSUN .GT. 0.) THEN
        PARSUN = (CAD(1)+FSUN*CAI(1)) * LAIFRA / MAX(LAISUN,MPE)
        PARSHA = (FSHA*CAI(1))*LAIFRA / MAX(LAISHA,MPE)
     ELSE
        PARSUN = 0.
        PARSHA = (CAD(1)+CAI(1))*LAIFRA /MAX(LAISHA,MPE)
     ENDIF

! reflected solar radiation

     RVIS = ALBD(1)*SOLAD(1) + ALBI(1)*SOLAI(1)
     RNIR = ALBD(2)*SOLAD(2) + ALBI(2)*SOLAI(2)
     FSR  = RVIS + RNIR

! reflected solar radiation of veg. and ground (combined ground)
     FSRV = FREVD(1)*SOLAD(1)+FREVI(1)*SOLAI(1)+FREVD(2)*SOLAD(2)+FREVI(2)*SOLAI(2)
     FSRG = FREGD(1)*SOLAD(1)+FREGI(1)*SOLAI(1)+FREGD(2)*SOLAD(2)+FREGI(2)*SOLAI(2)

  END SUBROUTINE SURRAD_d
  attributes(device) SUBROUTINE ALBEDO_d (VEGTYP ,IST    ,ISC    ,ICE    ,NSOIL  , & !in
                     DT     ,COSZ   ,FAGE   ,ELAI   ,ESAI   , & !in
                     TG     ,TV     ,SNOWH  ,FSNO   ,FWET   , & !in
                     SMC    ,SNEQVO ,SNEQV  ,QSNOW  ,FVEG   , & !in
                     ILOC   ,JLOC   ,                         & !in
                     ALBOLD ,TAUSS                          , & !inout
                     ALBGRD ,ALBGRI ,ALBD   ,ALBI   ,FABD   , & !out
                     FABI   ,FTDD   ,FTID   ,FTII   ,FSUN   , & !out
                     FREVI  ,FREVD  ,FREGD  ,FREGI  ,BGAP   , & !out
                     WGAP,rhol,rhos,taul,taus,opt_alb,omegas,opt_rad,tfrz,betais,betads,xl,hvt,hvb,rc,swemx,albsat,albdry,alblak)

! --------------------------------------------------------------------------------------------------
! surface albedos. also fluxes (per unit incoming direct and diffuse
! radiation) reflected, transmitted, and absorbed by vegetation.
! also sunlit fraction of the canopy.
! --------------------------------------------------------------------------------------------------
!   USE NOAHMP_VEG_PARAMETERS
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input
  INTEGER,                  INTENT(IN)  :: ILOC
  INTEGER,                  INTENT(IN)  :: JLOC
  INTEGER,                  INTENT(IN)  :: NSOIL  !number of soil layers
  INTEGER,                  INTENT(IN)  :: VEGTYP !vegetation type
  INTEGER,                  INTENT(IN)  :: IST    !surface type
  INTEGER,                  INTENT(IN)  :: ISC    !soil color type (1-lighest; 8-darkest)
  INTEGER,                  INTENT(IN)  :: ICE    !ice (ice = 1)

  REAL(rkind),                     INTENT(IN)  :: DT     !time step [sec]
  REAL(rkind),                     INTENT(IN)  :: QSNOW  !snowfall
  REAL(rkind),                     INTENT(IN)  :: COSZ   !cosine solar zenith angle for next time step
  REAL(rkind),                     INTENT(IN)  :: SNOWH  !snow height (mm)
  REAL(rkind),                     INTENT(IN)  :: TG     !ground temperature (k)
  REAL(rkind),                     INTENT(IN)  :: TV     !vegetation temperature (k)
  REAL(rkind),                     INTENT(IN)  :: ELAI   !LAI, one-sided, adjusted for burying by snow
  REAL(rkind),                     INTENT(IN)  :: ESAI   !SAI, one-sided, adjusted for burying by snow
  REAL(rkind),                     INTENT(IN)  :: FSNO   !fraction of grid covered by snow
  REAL(rkind),                     INTENT(IN)  :: FWET   !fraction of canopy that is wet
  REAL(rkind),                     INTENT(IN)  :: SNEQVO !snow mass at last time step(mm)
  REAL(rkind),                     INTENT(IN)  :: SNEQV  !snow mass (mm)
  REAL(rkind),                     INTENT(IN)  :: FVEG   !green vegetation fraction [0.0-1.0]
  REAL(rkind), DIMENSION(1:NSOIL), INTENT(IN)  :: SMC    !volumetric soil water (m3/m3)

! inout
  REAL(rkind),                  INTENT(INOUT)  :: ALBOLD !snow albedo at last time step (CLASS type)
  REAL(rkind),                  INTENT(INOUT)  :: TAUSS  !non-dimensional snow age

! output
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: ALBGRD !ground albedo (direct)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: ALBGRI !ground albedo (diffuse)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: ALBD   !surface albedo (direct)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: ALBI   !surface albedo (diffuse)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FABD   !flux abs by veg (per unit direct flux)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FABI   !flux abs by veg (per unit diffuse flux)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FTDD   !down direct flux below veg (per unit dir flux)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FTID   !down diffuse flux below veg (per unit dir flux)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FTII   !down diffuse flux below veg (per unit dif flux)
  REAL(rkind),                     INTENT(OUT) :: FSUN   !sunlit fraction of canopy (-)
!jref:start
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FREVD
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FREVI
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FREGD
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: FREGI
  REAL(rkind), INTENT(OUT) :: BGAP
  REAL(rkind), INTENT(OUT) :: WGAP
  real(rkind) :: rhol(:,:), rhos(:,:),taul(:,:),taus(:,:)
  integer(i4b) :: opt_alb
  real(rkind) :: omegas(:)
integer(i4b) :: opt_rad
real(rkind) :: tfrz,betais,betads
real(rkind) :: xl(:),hvt(:),hvb(:),rc(:)
real(rkind) :: swemx
real(rkind) :: albsat(:,:),albdry(:,:),alblak(:)


!jref:end

! ------------------------------------------------------------------------
! ------------------------ local variables -------------------------------
! local
  REAL(rkind)                 :: FAGE     !snow age function
  REAL(rkind)                 :: ALB
  INTEGER              :: IB       !indices
  INTEGER              :: NBAND    !number of solar radiation wave bands
  INTEGER              :: IC       !direct beam: ic=0; diffuse: ic=1

  REAL(rkind)                 :: WL       !fraction of LAI+SAI that is LAI
  REAL(rkind)                 :: WS       !fraction of LAI+SAI that is SAI
  REAL(rkind)                 :: MPE      !prevents overflow for division by zero

  REAL(rkind), DIMENSION(1:2) :: RHO      !leaf/stem reflectance weighted by fraction LAI and SAI
  REAL(rkind), DIMENSION(1:2) :: TAU      !leaf/stem transmittance weighted by fraction LAI and SAI
  REAL(rkind), DIMENSION(1:2) :: FTDI     !down direct flux below veg per unit dif flux = 0
  REAL(rkind), DIMENSION(1:2) :: ALBSND   !snow albedo (direct)
  REAL(rkind), DIMENSION(1:2) :: ALBSNI   !snow albedo (diffuse)

  REAL(rkind)                 :: VAI      !ELAI+ESAI
  REAL(rkind)                 :: GDIR     !average projected leaf/stem area in solar direction
  REAL(rkind)                 :: EXT      !optical depth direct beam per unit leaf + stem area

! --------------------------------------------------------------------------------------------------

  NBAND = 2
  MPE = 1.E-06
  BGAP = 0.
  WGAP = 0.

! initialize output because solar radiation only done if COSZ > 0

  DO IB = 1, NBAND
    ALBD(IB) = 0.
    ALBI(IB) = 0.
    ALBGRD(IB) = 0.
    ALBGRI(IB) = 0.
    FABD(IB) = 0.
    FABI(IB) = 0.
    FTDD(IB) = 0.
    FTID(IB) = 0.
    FTII(IB) = 0.
    IF (IB.EQ.1) FSUN = 0.
  END DO

  IF(COSZ <= 0) GOTO 100

! weight reflectance/transmittance by LAI and SAI

  DO IB = 1, NBAND
    VAI = ELAI + ESAI
    WL  = ELAI / MAX(VAI,MPE)
    WS  = ESAI / MAX(VAI,MPE)
    RHO(IB) = MAX(RHOL(VEGTYP,IB)*WL+RHOS(VEGTYP,IB)*WS, MPE)
    TAU(IB) = MAX(TAUL(VEGTYP,IB)*WL+TAUS(VEGTYP,IB)*WS, MPE)
  END DO

! snow age

   CALL SNOW_AGE_d (DT,TG,SNEQVO,SNEQV,TAUSS,FAGE,tfrz,swemx)

! snow albedos: only if COSZ > 0 and FSNO > 0

  IF(OPT_ALB == 1) &
     CALL SNOWALB_BATS_d (NBAND, FSNO,COSZ,FAGE,ALBSND,ALBSNI)
  IF(OPT_ALB == 2) THEN
     CALL SNOWALB_CLASS_d (NBAND,QSNOW,DT,ALB,ALBOLD,ALBSND,ALBSNI,ILOC,JLOC,swemx)
     ALBOLD = ALB
  END IF

! ground surface albedo

  CALL GROUNDALB_d (NSOIL   ,NBAND   ,ICE     ,IST     ,ISC     , & !in
                  FSNO    ,SMC     ,ALBSND  ,ALBSNI  ,COSZ    , & !in
                  TG      ,ILOC    ,JLOC    ,                   & !in
                  ALBGRD  ,ALBGRI  ,tfrz,albsat,albdry,alblak)                              !out


! loop over NBAND wavebands to calculate surface albedos and solar
! fluxes for unit incoming direct (IC=0) and diffuse flux (IC=1)

  DO IB = 1, NBAND
      IC = 0      ! direct
      CALL TWOSTREAM_d (IB     ,IC      ,VEGTYP  ,COSZ    ,VAI    , & !in
                      FWET   ,TV      ,ALBGRD  ,ALBGRI  ,RHO    , & !in
                      TAU    ,FVEG    ,IST     ,ILOC    ,JLOC   , & !in
                      FABD   ,ALBD    ,FTDD    ,FTID    ,GDIR   , &!)   !out
                      FREVD  ,FREGD   ,BGAP    ,WGAP,omegas,opt_rad,tfrz,betais,betads,xl,hvt,hvb,rc)

      IC = 1      ! diffuse
      CALL TWOSTREAM_d (IB     ,IC      ,VEGTYP  ,COSZ    ,VAI    , & !in
                      FWET   ,TV      ,ALBGRD  ,ALBGRI  ,RHO    , & !in
                      TAU    ,FVEG    ,IST     ,ILOC    ,JLOC   , & !in
                      FABI   ,ALBI    ,FTDI    ,FTII    ,GDIR   , & !)   !out
                      FREVI  ,FREGI   ,BGAP    ,WGAP,omegas,opt_rad,tfrz,betais,betads,xl,hvt,hvb,rc)

  END DO

! sunlit fraction of canopy. set FSUN = 0 if FSUN < 0.01.

  EXT = GDIR/COSZ * SQRT(1.-RHO(1)-TAU(1))
  FSUN = (1.-EXP(-EXT*VAI)) / MAX(EXT*VAI,MPE)
  EXT = FSUN

  IF (EXT .LT. 0.01) THEN
     WL = 0.
  ELSE
     WL = EXT
  END IF
  FSUN = WL

100 CONTINUE

  END SUBROUTINE ALBEDO_d

  attributes(device) SUBROUTINE SNOW_AGE_d (DT,TG,SNEQVO,SNEQV,TAUSS,FAGE,tfrz,swemx)
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! ------------------------ code history ------------------------------------------------------------
! from BATS
! ------------------------ input/output variables --------------------------------------------------
!input
   REAL(rkind), INTENT(IN) :: DT        !main time step (s)
   REAL(rkind), INTENT(IN) :: TG        !ground temperature (k)
   REAL(rkind), INTENT(IN) :: SNEQVO    !snow mass at last time step(mm)
   REAL(rkind), INTENT(IN) :: SNEQV     !snow water per unit ground area (mm)

!output
   REAL(rkind), INTENT(OUT) :: FAGE     !snow age

!input/output
   REAL(rkind), INTENT(INOUT) :: TAUSS      !non-dimensional snow age
   real(rkind) :: tfrz,swemx
!local
   REAL(rkind)            :: TAGE       !total aging effects
   REAL(rkind)            :: AGE1       !effects of grain growth due to vapor diffusion
   REAL(rkind)            :: AGE2       !effects of grain growth at freezing of melt water
   REAL(rkind)            :: AGE3       !effects of soot
   REAL(rkind)            :: DELA       !temporary variable
   REAL(rkind)            :: SGE        !temporary variable
   REAL(rkind)            :: DELS       !temporary variable
   REAL(rkind)            :: DELA0      !temporary variable
   REAL(rkind)            :: ARG        !temporary variable
! See Yang et al. (1997) J.of Climate for detail.
!---------------------------------------------------------------------------------------------------

   IF(SNEQV.LE.0.0) THEN
          TAUSS = 0.
   ELSE IF (SNEQV.GT.800.) THEN
          TAUSS = 0.
   ELSE
          DELA0 = 1.E-6*DT
          ARG   = 5.E3*(1./TFRZ-1./TG)
          AGE1  = EXP(ARG)
          AGE2  = EXP(AMIN1(0.,10.*ARG))
          AGE3  = 0.3
          TAGE  = AGE1+AGE2+AGE3

          DELA  = DELA0*TAGE
          DELS  = AMAX1(0.0,SNEQV-SNEQVO) / SWEMX
          SGE   = (TAUSS+DELA)*(1.0-DELS)
          TAUSS = AMAX1(0.,SGE)
   ENDIF

   FAGE= TAUSS/(TAUSS+1.)

  END SUBROUTINE SNOW_AGE_d
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
  attributes(device) SUBROUTINE SNOWALB_BATS_d (NBAND,FSNO,COSZ,FAGE,ALBSND,ALBSNI)
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input

  INTEGER,INTENT(IN) :: NBAND  !number of waveband classes

  REAL(rkind),INTENT(IN) :: COSZ    !cosine solar zenith angle
  REAL(rkind),INTENT(IN) :: FSNO    !snow cover fraction (-)
  REAL(rkind),INTENT(IN) :: FAGE    !snow age correction

! output

  REAL(rkind), DIMENSION(1:2),INTENT(OUT) :: ALBSND !snow albedo for direct(1=vis, 2=nir)
  REAL(rkind), DIMENSION(1:2),INTENT(OUT) :: ALBSNI !snow albedo for diffuse
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
  INTEGER :: IB          !waveband class

  REAL(rkind) :: FZEN                 !zenith angle correction
  REAL(rkind) :: CF1                  !temperary variable
  REAL(rkind) :: SL2                  !2.*SL
  REAL(rkind) :: SL1                  !1/SL
  REAL(rkind) :: SL                   !adjustable parameter
  REAL(rkind), PARAMETER :: C1 = 0.2  !default in BATS
  REAL(rkind), PARAMETER :: C2 = 0.5  !default in BATS
!  REAL(rkind), PARAMETER :: C1 = 0.2 * 2. ! double the default to match Sleepers River's
!  REAL(rkind), PARAMETER :: C2 = 0.5 * 2. ! snow surface albedo (double aging effects)
! ---------------------------------------------------------------------------------------------
! zero albedos for all points

        ALBSND(1: NBAND) = 0.
        ALBSNI(1: NBAND) = 0.

! when cosz > 0

        SL=2.0
        SL1=1./SL
        SL2=2.*SL
        CF1=((1.+SL1)/(1.+SL2*COSZ)-SL1)
        FZEN=AMAX1(CF1,0.)

        ALBSNI(1)=0.95*(1.-C1*FAGE)
        ALBSNI(2)=0.65*(1.-C2*FAGE)

        ALBSND(1)=ALBSNI(1)+0.4*FZEN*(1.-ALBSNI(1))    !  vis direct
        ALBSND(2)=ALBSNI(2)+0.4*FZEN*(1.-ALBSNI(2))    !  nir direct

  END SUBROUTINE SNOWALB_BATS_d
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
  attributes(device) SUBROUTINE SNOWALB_CLASS_d (NBAND,QSNOW,DT,ALB,ALBOLD,ALBSND,ALBSNI,ILOC,JLOC,swemx)
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input

  INTEGER,INTENT(IN) :: ILOC !grid index
  INTEGER,INTENT(IN) :: JLOC !grid index
  INTEGER,INTENT(IN) :: NBAND  !number of waveband classes

  REAL(rkind),INTENT(IN) :: QSNOW     !snowfall (mm/s)
  REAL(rkind),INTENT(IN) :: DT        !time step (sec)
  REAL(rkind),INTENT(IN) :: ALBOLD    !snow albedo at last time step

! in & out

  REAL(rkind),                INTENT(INOUT) :: ALB        !
! output

  REAL(rkind), DIMENSION(1:2),INTENT(OUT) :: ALBSND !snow albedo for direct(1=vis, 2=nir)
  REAL(rkind), DIMENSION(1:2),INTENT(OUT) :: ALBSNI !snow albedo for diffuse
  real(rkind) :: swemx
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
  INTEGER :: IB          !waveband class

! ---------------------------------------------------------------------------------------------
! zero albedos for all points

        ALBSND(1: NBAND) = 0.
        ALBSNI(1: NBAND) = 0.

! when cosz > 0

         ALB = 0.55 + (ALBOLD-0.55) * EXP(-0.01*DT/3600.)


! 1 mm fresh snow(SWE) -- 10mm snow depth, assumed the fresh snow density 100kg/m3
! here assume 1cm snow depth will fully cover the old snow

         IF (QSNOW > 0.) then
           ALB = ALB + MIN(QSNOW*DT,SWEMX) * (0.84-ALB)/(SWEMX)
         ENDIF

         ALBSNI(1)= ALB         ! vis diffuse
         ALBSNI(2)= ALB         ! nir diffuse
         ALBSND(1)= ALB         ! vis direct
         ALBSND(2)= ALB         ! nir direct

  END SUBROUTINE SNOWALB_CLASS_d
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
  attributes(device) SUBROUTINE GROUNDALB_d (NSOIL   ,NBAND   ,ICE     ,IST     ,ISC     , & !in
                        FSNO    ,SMC     ,ALBSND  ,ALBSNI  ,COSZ    , & !in
                        TG      ,ILOC    ,JLOC    ,                   & !in
                        ALBGRD  ,ALBGRI  ,tfrz,albsat,albdry,alblak)                              !out
! --------------------------------------------------------------------------------------------------
!   USE NOAHMP_RAD_PARAMETERS
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
!input

  INTEGER,                  INTENT(IN)  :: ILOC   !grid index
  INTEGER,                  INTENT(IN)  :: JLOC   !grid index
  INTEGER,                  INTENT(IN)  :: NSOIL  !number of soil layers
  INTEGER,                  INTENT(IN)  :: NBAND  !number of solar radiation waveband classes
  INTEGER,                  INTENT(IN)  :: ICE    !value of ist for land ice
  INTEGER,                  INTENT(IN)  :: IST    !surface type
  INTEGER,                  INTENT(IN)  :: ISC    !soil color class (1-lighest; 8-darkest)
  REAL(rkind),                     INTENT(IN)  :: FSNO   !fraction of surface covered with snow (-)
  REAL(rkind),                     INTENT(IN)  :: TG     !ground temperature (k)
  REAL(rkind),                     INTENT(IN)  :: COSZ   !cosine solar zenith angle (0-1)
  REAL(rkind), DIMENSION(1:NSOIL), INTENT(IN)  :: SMC    !volumetric soil water content (m3/m3)
  REAL(rkind), DIMENSION(1:    2), INTENT(IN)  :: ALBSND !direct beam snow albedo (vis, nir)
  REAL(rkind), DIMENSION(1:    2), INTENT(IN)  :: ALBSNI !diffuse snow albedo (vis, nir)

!output

  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: ALBGRD !ground albedo (direct beam: vis, nir)
  REAL(rkind), DIMENSION(1:    2), INTENT(OUT) :: ALBGRI !ground albedo (diffuse: vis, nir)
real(rkind) :: tfrz, albsat(:,:),albdry(:,:),alblak(:)
!local

  INTEGER                               :: IB     !waveband number (1=vis, 2=nir)
  REAL(rkind)                                  :: INC    !soil water correction factor for soil albedo
  REAL(rkind)                                  :: ALBSOD !soil albedo (direct)
  REAL(rkind)                                  :: ALBSOI !soil albedo (diffuse)
! --------------------------------------------------------------------------------------------------

  DO IB = 1, NBAND
        INC = MAX(0.11-0.40*SMC(1), 0.)
        IF (IST .EQ. 1)  THEN                     !soil
           ALBSOD = MIN(ALBSAT(ISC,IB)+INC,ALBDRY(ISC,IB))
           ALBSOI = ALBSOD
        ELSE IF (TG .GT. TFRZ) THEN               !unfrozen lake, wetland
           ALBSOD = 0.06/(MAX(0.01,COSZ)**1.7 + 0.15)
           ALBSOI = 0.06
        ELSE                                      !frozen lake, wetland
           ALBSOD = ALBLAK(IB)
           ALBSOI = ALBSOD
        END IF

! increase desert and semi-desert albedos

        IF (IST .EQ. 1 .AND. ISC .EQ. 9) THEN
           ALBSOD = ALBSOD + 0.10
           ALBSOI = ALBSOI + 0.10
        end if

        ALBGRD(IB) = ALBSOD*(1.-FSNO) + ALBSND(IB)*FSNO
        ALBGRI(IB) = ALBSOI*(1.-FSNO) + ALBSNI(IB)*FSNO
  END DO

  END SUBROUTINE GROUNDALB_d


end module vegSWavRad_module
