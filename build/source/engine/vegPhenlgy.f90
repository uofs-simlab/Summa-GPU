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

module vegPhenlgy_module

! data types
USE nrtype

! global variables
USE globalData,only:urbanVegCategory    ! vegetation category for urban areas

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,            & ! data vector (i4b)
                    var_d,            & ! data vector (rkind)
                    var_dlength,      & ! data vector with variable length dimension (rkind)
                    model_options       ! defines the model decisions

! named variables defining elements in the data structures
USE var_lookup,only:iLookTYPE,iLookATTR,iLookPARAM,iLookDIAG,iLookPROG  ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS                                      ! named variables for elements of the decision structure

! look-up values for the boundary conditions
USE mDecisions_module,only:      &
 prescribedHead,                 &      ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
 prescribedTemp,                 &      ! prescribed temperature
 zeroFlux                               ! zero flux

! look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:      &
 noah_mp,                        &      ! full Noah-MP implementation (including albedo)
 CLM_2stream,                    &      ! CLM 2-stream model (see CLM documentation)
 UEB_2stream,                    &      ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
 NL_scatter,                     &      ! Simplified method Nijssen and Lettenmaier (JGR 1999)
 BeersLaw                               ! Beer's Law (as implemented in VIC)

! privacy
implicit none
private
public::vegPhenlgy,vegPhenlgy_d
! algorithmic parameters
real(rkind),parameter     :: valueMissing=-9999._rkind  ! missing value, used when diagnostic or state variables are undefined
real(rkind),parameter     :: verySmall=1.e-6_rkind      ! used as an additive constant to check if substantial difference among real numbers
contains


 ! ************************************************************************************************
 ! public subroutine vegPhenlgy: compute vegetation phenology
 ! ************************************************************************************************
 subroutine vegPhenlgy(&
                       ! model control
                       model_decisions,             & ! intent(in):    model decisions
                       fracJulDay,                  & ! intent(in):    fractional julian days since the start of year
                       yearLength,                  & ! intent(in):    number of days in the current year
                       ! input/output: data structures
                       type_data,                   & ! intent(in):    type of vegetation and soil
                       attr_data,                   & ! intent(in):    spatial attributes
                       mpar_data,                   & ! intent(in):    model parameters
                       prog_data,                   & ! intent(inout): prognostic variables for a local HRU
                       diag_data,                   & ! intent(inout): diagnostic variables for a local HRU
                       ! output
                       computeVegFlux,              & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                       canopyDepth,                 & ! intent(out): canopy depth (m)
                       exposedVAI,                  & ! intent(out): exposed vegetation area index (LAI + SAI)
                       err,message)                   ! intent(out): error control

 ! -------------------------------------------------------------------------------------------------
 ! modules
 USE NOAHMP_ROUTINES,only:phenology         ! determine vegetation phenology
 implicit none
 ! -------------------------------------------------------------------------------------------------
 ! input/output
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 real(rkind),intent(in)          :: fracJulDay          ! fractional julian days since the start of year
 integer(i4b),intent(in)         :: yearLength          ! number of days in the current year
 type(var_i),intent(in)          :: type_data           ! type of vegetation and soil
 type(var_d),intent(in)          :: attr_data           ! spatial attributes
 type(var_dlength),intent(in)    :: mpar_data           ! model parameters
 type(var_dlength),intent(inout) :: prog_data           ! prognostic variables for a local HRU
 type(var_dlength),intent(inout) :: diag_data           ! diagnostic variables for a local HRU
 ! output
 logical(lgt),intent(out)        :: computeVegFlux      ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(rkind),intent(out)         :: canopyDepth         ! canopy depth (m)
 real(rkind),intent(out)         :: exposedVAI          ! exposed vegetation area index (LAI + SAI)
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! -------------------------------------------------------------------------------------------------
 ! local
 real(rkind)                     :: notUsed_heightCanopyTop    ! height of the top of the canopy layer (m)
 real(rkind)                     :: heightAboveSnow            ! height top of canopy is above the snow surface (m)
 ! initialize error control
 err=0; message="vegPhenlgy/"
 ! ----------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&
 ! input: model decisions
 ix_bcUpprTdyn                   => model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision,      & ! intent(in): [i4b] choice of upper boundary condition for thermodynamics
 ix_bcUpprSoiH                   => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision,      & ! intent(in): [i4b] index of method used for the upper boundary condition for soil hydrology
 ! local attributes
 vegTypeIndex                    => type_data%var(iLookTYPE%vegTypeIndex),                     & ! intent(in): [i4b] vegetation type index
 latitude                        => attr_data%var(iLookATTR%latitude),                         & ! intent(in): [dp] latitude
 ! model state variables
 scalarSnowDepth                 => prog_data%var(iLookPROG%scalarSnowDepth)%dat(1),           & ! intent(in):    [dp] snow depth on the ground surface (m)
 scalarCanopyTemp                => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1),          & ! intent(in):    [dp] temperature of the vegetation canopy at the start of the sub-step (K)
 scalarCanopyLiq                 => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1),           & ! intent(inout): [dp] liquid water in the vegetation canopy at the start of the sub-step
 ! diagnostic variables and parameters (input)
 heightCanopyTop                 => mpar_data%var(iLookPARAM%heightCanopyTop)%dat(1),          & ! intent(in):    [dp] height of the top of the canopy layer (m)
 heightCanopyBottom              => mpar_data%var(iLookPARAM%heightCanopyBottom)%dat(1),       & ! intent(in):    [dp] height of the bottom of the canopy layer (m)
 scalarRootZoneTemp              => diag_data%var(iLookDIAG%scalarRootZoneTemp)%dat(1),        & ! intent(in):    [dp] root zone temperature (K)
 ! diagnostic variables and parameters (input/output)
 scalarLAI                       => diag_data%var(iLookDIAG%scalarLAI)%dat(1),                 & ! intent(inout): [dp] one-sided leaf area index (m2 m-2)
 scalarSAI                       => diag_data%var(iLookDIAG%scalarSAI)%dat(1),                 & ! intent(inout): [dp] one-sided stem area index (m2 m-2)
 ! diagnostic variables and parameters (output)
 scalarExposedLAI                => diag_data%var(iLookDIAG%scalarExposedLAI)%dat(1),          & ! intent(out): [dp] exposed leaf area index after burial by snow (m2 m-2)
 scalarExposedSAI                => diag_data%var(iLookDIAG%scalarExposedSAI)%dat(1),          & ! intent(out): [dp] exposed stem area index after burial by snow (m2 m-2)
 scalarGrowingSeasonIndex        => diag_data%var(iLookDIAG%scalarGrowingSeasonIndex)%dat(1)   & ! intent(out): [dp] growing season index (0=off, 1=on)
 ) ! associate variables in data structure
 ! ----------------------------------------------------------------------------------------------------------------------------------

 ! check if we have isolated the snow-soil domain (used in test cases)
 if(ix_bcUpprTdyn == prescribedTemp .or. ix_bcUpprTdyn == zeroFlux .or. ix_bcUpprSoiH == prescribedHead)then

  ! isolated snow-soil domain: do not compute fluxes over vegetation
  computeVegFlux = .false.

  ! set vegetation phenology variables to missing
  scalarLAI                = valueMissing    ! one-sided leaf area index (m2 m-2)
  scalarSAI                = valueMissing    ! one-sided stem area index (m2 m-2)
  scalarExposedLAI         = valueMissing    ! exposed leaf area index after burial by snow (m2 m-2)
  scalarExposedSAI         = valueMissing    ! exposed stem area index after burial by snow (m2 m-2)
  scalarGrowingSeasonIndex = valueMissing    ! growing season index (0=off, 1=on)
  exposedVAI               = valueMissing    ! exposed vegetation area index (m2 m-2)
  canopyDepth              = valueMissing    ! canopy depth (m)
  heightAboveSnow          = valueMissing    ! height top of canopy is above the snow surface (m)

 ! compute vegetation phenology (checks for complete burial of vegetation)
 else

  ! determine vegetation phenology
  ! NOTE: recomputing phenology every sub-step accounts for changes in exposed vegetation associated with changes in snow depth
  call phenology(&
                 ! input
                 vegTypeIndex,                & ! intent(in): vegetation type index
                 urbanVegCategory,            & ! intent(in): vegetation category for urban areas
                 scalarSnowDepth,             & ! intent(in): snow depth on the ground surface (m)
                 scalarCanopyTemp,            & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                 latitude,                    & ! intent(in): latitude
                 yearLength,                  & ! intent(in): number of days in the current year
                 fracJulDay,                  & ! intent(in): fractional julian days since the start of year
                 scalarLAI,                   & ! intent(inout): one-sided leaf area index (m2 m-2)
                 scalarSAI,                   & ! intent(inout): one-sided stem area index (m2 m-2)
                 scalarRootZoneTemp,          & ! intent(in): root zone temperature (K)
                 ! output
                 notUsed_heightCanopyTop,     & ! intent(out): height of the top of the canopy layer (m)
                 scalarExposedLAI,            & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                 scalarExposedSAI,            & ! intent(out): exposed stem area index after burial by snow (m2 m-2)
                 scalarGrowingSeasonIndex     ) ! intent(out): growing season index (0=off, 1=on)

  ! determine additional phenological variables
  exposedVAI      = scalarExposedLAI + scalarExposedSAI   ! exposed vegetation area index (m2 m-2)
  canopyDepth     = heightCanopyTop - heightCanopyBottom  ! canopy depth (m)
  heightAboveSnow = heightCanopyTop - scalarSnowDepth     ! height top of canopy is above the snow surface (m)

  ! determine if need to include vegetation in the energy flux routines
  computeVegFlux = (exposedVAI > 0.05_rkind .and. heightAboveSnow > 0.05_rkind)

  ! if no vegetation ever, should not have initialized scalarCanopyLiq to 0.0001 in read_icond.f90
  if((scalarLAI + scalarSAI) == 0.0_rkind) scalarCanopyLiq = 0.0_rkind

 end if  ! (check if the snow-soil column is isolated)

 ! end association to variables in the data structure
 end associate

 end subroutine vegPhenlgy



 ! ************************************************************************************************
 ! public subroutine vegPhenlgy: compute vegetation phenology
 ! ************************************************************************************************
 subroutine vegPhenlgy_d(&
                       ! model control
                       model_decisions,             & ! intent(in):    model decisions
                       nGRU, &
                       fracJulDay,                  & ! intent(in):    fractional julian days since the start of year
                       yearLength,                  & ! intent(in):    number of days in the current year
                       ! input/output: data structures
                       type_data,                   & ! intent(in):    type of vegetation and soil
                       attr_data,                   & ! intent(in):    spatial attributes
                       mpar_data,                   & ! intent(in):    model parameters
                       prog_data,                   & ! intent(inout): prognostic variables for a local HRU
                       diag_data,                   & ! intent(inout): diagnostic variables for a local HRU
                       decisions, veg_param, &
                       ! output
                       computeVegFlux,              & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                       canopyDepth,                 & ! intent(out): canopy depth (m)
                       exposedVAI,                  & ! intent(out): exposed vegetation area index (LAI + SAI)
                       err,message)                   ! intent(out): error control

 ! -------------------------------------------------------------------------------------------------
 ! modules
!  USE NOAHMP_ROUTINES,only:phenology_d         ! determine vegetation phenology
 use device_data_types
 implicit none
 ! -------------------------------------------------------------------------------------------------
 ! input/output
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 integer(i4b) :: nGRU
 real(rkind),intent(in),device          :: fracJulDay          ! fractional julian days since the start of year
 integer(i4b),intent(in),device         :: yearLength          ! number of days in the current year
 type(type_data_device),intent(in)          :: type_data           ! type of vegetation and soil
 type(attr_data_device),intent(in)          :: attr_data           ! spatial attributes
 type(mpar_data_device),intent(in)    :: mpar_data           ! model parameters
 type(prog_data_device),intent(inout) :: prog_data           ! prognostic variables for a local HRU
 type(diag_data_device),intent(inout) :: diag_data           ! diagnostic variables for a local HRU
 type(decisions_device) :: decisions
 type(veg_parameters) :: veg_param
 ! output
 logical(lgt),intent(out),device        :: computeVegFlux(:)      ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(rkind),intent(out),device         :: canopyDepth(:)         ! canopy depth (m)
 real(rkind),intent(out),device         :: exposedVAI(:)          ! exposed vegetation area index (LAI + SAI)
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! -------------------------------------------------------------------------------------------------
 ! local
 real(rkind),device                     :: notUsed_heightCanopyTop(nGRU)    ! height of the top of the canopy layer (m)
 real(rkind)                     :: heightAboveSnow            ! height top of canopy is above the snow surface (m)
 integer(i4b) :: iGRU
 ! initialize error control
 err=0; message="vegPhenlgy/"
 ! ----------------------------------------------------------------------------------------------------------------------------------
 ! associate variables in the data structure
 associate(&
 ! input: model decisions
 ix_bcUpprTdyn                   => model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision,      & ! intent(in): [i4b] choice of upper boundary condition for thermodynamics
 ix_bcUpprSoiH                   => model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision,      & ! intent(in): [i4b] index of method used for the upper boundary condition for soil hydrology
 hvt => veg_param%hvt, &
 tmin => veg_param%tmin, &
 hvb => veg_param%hvb, &
 isbarren => veg_param%isbarren, &
 ISSNOW => veg_param%ISSNOW, &
 ISWATER => veg_param%ISWATER, &
 saim => veg_param%saim, &
 laim => veg_param%laim, &
 dveg => veg_param%dveg, &
 urbanVegCategory => decisions%urbanVegCategory, &
 ! local attributes
 vegTypeIndex                    => type_data%vegTypeIndex,                     & ! intent(in): [i4b] vegetation type index
 latitude                        => attr_data%latitude,                         & ! intent(in): [dp] latitude
 ! model state variables
 scalarSnowDepth                 => prog_data%scalarSnowDepth,           & ! intent(in):    [dp] snow depth on the ground surface (m)
 scalarCanopyTemp                => prog_data%scalarCanopyTemp,          & ! intent(in):    [dp] temperature of the vegetation canopy at the start of the sub-step (K)
 scalarCanopyLiq                 => prog_data%scalarCanopyLiq,           & ! intent(inout): [dp] liquid water in the vegetation canopy at the start of the sub-step
 ! diagnostic variables and parameters (input)
 heightCanopyTop                 => mpar_data%heightCanopyTop,          & ! intent(in):    [dp] height of the top of the canopy layer (m)
 heightCanopyBottom              => mpar_data%heightCanopyBottom,       & ! intent(in):    [dp] height of the bottom of the canopy layer (m)
 scalarRootZoneTemp              => diag_data%scalarRootZoneTemp,        & ! intent(in):    [dp] root zone temperature (K)
 ! diagnostic variables and parameters (input/output)
 scalarLAI                       => diag_data%scalarLAI,                 & ! intent(inout): [dp] one-sided leaf area index (m2 m-2)
 scalarSAI                       => diag_data%scalarSAI,                 & ! intent(inout): [dp] one-sided stem area index (m2 m-2)
 ! diagnostic variables and parameters (output)
 scalarExposedLAI                => diag_data%scalarExposedLAI,          & ! intent(out): [dp] exposed leaf area index after burial by snow (m2 m-2)
 scalarExposedSAI                => diag_data%scalarExposedSAI,          & ! intent(out): [dp] exposed stem area index after burial by snow (m2 m-2)
 scalarGrowingSeasonIndex        => diag_data%scalarGrowingSeasonIndex   & ! intent(out): [dp] growing season index (0=off, 1=on)
 ) ! associate variables in data structure
 ! ----------------------------------------------------------------------------------------------------------------------------------

 ! check if we have isolated the snow-soil domain (used in test cases)
 if(ix_bcUpprTdyn == prescribedTemp .or. ix_bcUpprTdyn == zeroFlux .or. ix_bcUpprSoiH == prescribedHead)then

  ! isolated snow-soil domain: do not compute fluxes over vegetation
  computeVegFlux = .false.

  ! set vegetation phenology variables to missing
  scalarLAI                = valueMissing    ! one-sided leaf area index (m2 m-2)
  scalarSAI                = valueMissing    ! one-sided stem area index (m2 m-2)
  scalarExposedLAI         = valueMissing    ! exposed leaf area index after burial by snow (m2 m-2)
  scalarExposedSAI         = valueMissing    ! exposed stem area index after burial by snow (m2 m-2)
  scalarGrowingSeasonIndex = valueMissing    ! growing season index (0=off, 1=on)
  exposedVAI               = valueMissing    ! exposed vegetation area index (m2 m-2)
  canopyDepth              = valueMissing    ! canopy depth (m)
  heightAboveSnow          = valueMissing    ! height top of canopy is above the snow surface (m)

 ! compute vegetation phenology (checks for complete burial of vegetation)
 else
  ! determine vegetation phenology
  ! NOTE: recomputing phenology every sub-step accounts for changes in exposed vegetation associated with changes in snow depth
    do iGRU=1,nGRU
  call phenology_d<<<1,1>>>(&
                 ! input
                 vegTypeIndex(iGRU),                & ! intent(in): vegetation type index
                 urbanVegCategory,            & ! intent(in): vegetation category for urban areas
                 scalarSnowDepth(iGRU),             & ! intent(in): snow depth on the ground surface (m)
                 scalarCanopyTemp(iGRU),            & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                 latitude(iGRU),                    & ! intent(in): latitude
                 yearLength,                  & ! intent(in): number of days in the current year
                 fracJulDay,                  & ! intent(in): fractional julian days since the start of year
                 scalarLAI(iGRU),                   & ! intent(inout): one-sided leaf area index (m2 m-2)
                 scalarSAI(iGRU),                   & ! intent(inout): one-sided stem area index (m2 m-2)
                 scalarRootZoneTemp(iGRU),          & ! intent(in): root zone temperature (K)
                 ! output
                 notUsed_heightCanopyTop(iGRU),     & ! intent(out): height of the top of the canopy layer (m)
                 scalarExposedLAI(iGRU),            & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                 scalarExposedSAI(iGRU),            & ! intent(out): exposed stem area index after burial by snow (m2 m-2)
                 scalarGrowingSeasonIndex(iGRU), &
                 hvt,tmin,hvb,&
                 isbarren,ISSNOW, ISWATER, &
                 saim,laim,&
                 dveg     ) ! intent(out): growing season index (0=off, 1=on)
    end do
    !$cuf kernel do(1) <<<*,*>>>
    do iGRU=1,nGRU

  ! determine additional phenological variables
  exposedVAI(iGRU)      = scalarExposedLAI(iGRU) + scalarExposedSAI(iGRU)   ! exposed vegetation area index (m2 m-2)
  canopyDepth(iGRU)     = heightCanopyTop - heightCanopyBottom  ! canopy depth (m)
  heightAboveSnow = heightCanopyTop - scalarSnowDepth(iGRU)     ! height top of canopy is above the snow surface (m)

  ! determine if need to include vegetation in the energy flux routines
  computeVegFlux(iGRU) = (exposedVAI(iGRU) > 0.05_rkind .and. heightAboveSnow > 0.05_rkind)

  ! if no vegetation ever, should not have initialized scalarCanopyLiq to 0.0001 in read_icond.f90
  if((scalarLAI(iGRU) + scalarSAI(iGRU)) == 0.0_rkind) scalarCanopyLiq(iGRU) = 0.0_rkind
    end do
 end if  ! (check if the snow-soil column is isolated)
 ! end association to variables in the data structure
 end associate

 end subroutine vegPhenlgy_d

 attributes(global) SUBROUTINE PHENOLOGY_d (VEGTYP , ISURBAN, SNOWH  , TV     , LAT   , YEARLEN , JULIAN , & !in
 LAI    , SAI    , TROOT  , HTOP  , ELAI    , ESAI   , IGS, &
 HVT, tmin, hvb, &
 ISBARREN, ISSNOW, ISWATER,&
  saim,laim,&
  dveg)

! --------------------------------------------------------------------------------------------------
! vegetation phenology considering vegeation canopy being buries by snow and evolution in time
! --------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------
IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! inputs
INTEGER                , INTENT(IN   ) :: VEGTYP !vegetation type
INTEGER                , INTENT(IN   ) :: ISURBAN!urban category
REAL(rkind)            , INTENT(IN   ) :: SNOWH  !snow height [m]
REAL(rkind)            , INTENT(IN   ) :: TV     !vegetation temperature (k)
REAL(rkind)            , INTENT(IN   ) :: LAT    !latitude (radians)
INTEGER                , INTENT(IN   ) :: YEARLEN!Number of days in the particular year
REAL(rkind)            , INTENT(IN   ) :: JULIAN !Julian day of year (fractional) ( 0 <= JULIAN < YEARLEN )
real(rkind)            , INTENT(IN   ) :: TROOT  !root-zone averaged temperature (k)
REAL(rkind)            , INTENT(INOUT) :: LAI    !LAI, unadjusted for burying by snow
REAL(rkind)            , INTENT(INOUT) :: SAI    !SAI, unadjusted for burying by snow

! outputs
REAL(rkind)            , INTENT(OUT  ) :: HTOP   !top of canopy layer (m)
REAL(rkind)            , INTENT(OUT  ) :: ELAI   !leaf area index, after burying by snow
REAL(rkind)            , INTENT(OUT  ) :: ESAI   !stem area index, after burying by snow
REAL(rkind)            , INTENT(OUT  ) :: IGS    !growing season index (0=off, 1=on)
REAL(rkind) :: hvt(:), tmin(:), hvb(:)
integer(i4b) :: ISBARREN, ISSNOW, ISWATER, dveg
real(rkind) :: saim(:,:), laim(:,:)

! locals

REAL(rkind)                            :: DB     !thickness of canopy buried by snow (m)
REAL(rkind)                            :: FB     !fraction of canopy buried by snow
REAL(rkind)                            :: SNOWHC !critical snow depth at which short vege
                            !is fully covered by snow

INTEGER                                :: K       !index
INTEGER                                :: IT1,IT2 !interpolation months
REAL(rkind)                            :: DAY     !current day of year ( 0 <= DAY < YEARLEN )
REAL(rkind)                            :: WT1,WT2 !interpolation weights
REAL(rkind)                            :: T       !current month (1.00, ..., 12.00)
! --------------------------------------------------------------------------------------------------

! Interpolate monthly SAI and LAI to daily values
IF ( DVEG == 1 .or. DVEG == 3 .or. DVEG == 4 ) THEN

IF (LAT >= 0.) THEN
! Northern Hemisphere
DAY = JULIAN
ELSE
! Southern Hemisphere.  DAY is shifted by 1/2 year.
DAY = MOD ( JULIAN + ( 0.5 * YEARLEN ) , REAL(YEARLEN) )
ENDIF

T = 12. * DAY / REAL(YEARLEN)
IT1 = T + 0.5
IT2 = IT1 + 1
WT1 = (IT1+0.5) - T
WT2 = 1.-WT1
IF (IT1 .LT.  1) IT1 = 12
IF (IT2 .GT. 12) IT2 = 1

LAI = WT1*LAIM(VEGTYP,IT1) + WT2*LAIM(VEGTYP,IT2)
SAI = WT1*SAIM(VEGTYP,IT1) + WT2*SAIM(VEGTYP,IT2)
ENDIF

! Realism check: no leaves without stems
IF (SAI == 0.0) LAI = 0.0  

! Realism check: warn about no stems for vegetated land classes
IF ( (SAI == 0.0) .and. ( VEGTYP /= ISWATER ) .and. ( VEGTYP /= ISBARREN ) .and. ( VEGTYP /= ISSNOW ) .and. ( VEGTYP /= ISURBAN) ) THEN
 print*, ' WARNING: module_sf_noahmplsm/PHENOLOGY: Stem Area Index (SAI) = 0.0 may be unrealistic for vegetation type ',VEGTYP,'. Continuing.'
ENDIF

! Realism check: no vegetation should exist on certain land classes
IF ( ( VEGTYP == ISWATER ) .or. ( VEGTYP == ISBARREN ) .or. ( VEGTYP == ISSNOW ) .or. ( VEGTYP == ISURBAN) ) THEN
LAI  = 0.
SAI  = 0.
ENDIF

!buried by snow

DB = MIN( MAX(SNOWH - HVB(VEGTYP),0.), HVT(VEGTYP)-HVB(VEGTYP) )
FB = DB / MAX(1.E-06,HVT(VEGTYP)-HVB(VEGTYP))
!print*, 'HVB(VEGTYP), HVT(VEGTYP), DB, FB = ', HVB(VEGTYP), HVT(VEGTYP), DB, FB

IF(HVT(VEGTYP)> 0. .AND. HVT(VEGTYP) <= 0.5) THEN
SNOWHC = HVT(VEGTYP)*EXP(-SNOWH/0.1)
IF (SNOWH>SNOWHC) THEN
FB = 1.
ELSE
FB = SNOWH/SNOWHC
ENDIF
ENDIF

ELAI =  LAI*(1.-FB)
ESAI =  SAI*(1.-FB)

IF (TV .GT. TMIN(VEGTYP)) THEN
IGS = 1.
ELSE
IGS = 0.
ENDIF

HTOP = HVT(VEGTYP)

END SUBROUTINE PHENOLOGY_d

end module vegPhenlgy_module
