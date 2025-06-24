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

module stomResist_module

  ! data types
  USE nrtype
  
  ! physical constants
  USE multiconst, only: Rgas     ! universal gas constant (J mol-1 K-1)
  USE multiconst, only: Tfreeze  ! freezing point of pure water (K)
  USE multiconst, only: ave_slp  ! standard pressure (Pa)
  
  ! derived types to define the data structures
  USE data_types,only:&
                      var_i,            & ! data vector (i4b)
                      var_d,            & ! data vector (rkind)
                      var_dlength,      & ! data vector with variable length dimension (rkind)
                      model_options       ! defines the model decisions
  
  ! indices that define elements of the data structures
  USE var_lookup,only:iLookTYPE           ! named variables for structure elements
  USE var_lookup,only:iLookDIAG           ! named variables for structure elements
  USE var_lookup,only:iLookFLUX           ! named variables for structure elements
  USE var_lookup,only:iLookFORCE          ! named variables for structure elements
  USE var_lookup,only:iLookPARAM          ! named variables for structure elements
  USE var_lookup,only:iLookDECISIONS                           ! named variables for elements of the decision structure
  
  ! look-up values for the stomatal resistance formulation
  USE mDecisions_module,only:  &
   simpleResistance,           & ! simple resistance formulation
   BallBerryFlex,              & ! flexible Ball-Berry scheme
   BallBerry,                  & ! Ball-Berry (from Noah-MP)
   Jarvis                        ! Jarvis (from Noah-MP)
  
  ! look-up values for the leaf temperature controls on photosynthesis + stomatal resistance
  USE mDecisions_module,only:  &
   q10Func,                    & ! the q10 function used in CLM4 and Noah-MP
   Arrhenius                     ! the Arrhenius functions used in CLM5 and Cable
  
  ! look-up values for the humidity controls on stomatal resistance
  USE mDecisions_module,only:  &
   humidLeafSurface,           & ! humidity at the leaf surface [Bonan et al., 2011]
   scaledHyperbolic              ! scaled hyperbolic function [Leuning et al., 1995]
  
  ! look-up values for the electron transport function, dependence of photosynthesis on PAR
  USE mDecisions_module,only:  &
   linear,                     & ! linear function used in CLM4 and Noah-MP
   linearJmax,                 & ! linear jmax function used in Cable [Wang et al., Ag Forest Met 1998, eq D5]
   quadraticJmax                 ! the quadratic Jmax function, used in SSiB and CLM5
  
  ! look-up values for the CO2 compensation point to calculate stomatal resistance
  USE mDecisions_module,only:  &
   origBWB,                    & ! the original BWB function
   Leuning                       ! the Leuning function
  
  ! look up values to define the iterative numerical solution method used in the Ball-Berry stomatal resistance parameterization
  USE mDecisions_module,only:  &
   NoahMPsolution,             & ! the NoahMP solution (and CLM4): fixed point iteration; max 3 iterations
   newtonRaphson                 ! full Newton-Raphson iterative solution to convergence
  
  ! look up values to define the controls on carbon assimilation
  USE mDecisions_module,only:  &
   colimitation,               & ! enable colimitation, as described by Collatz et al. (1991) and Sellers et al. (1996)
   minFunc                       ! do not enable colimitation: use minimum of the three controls on carbon assimilation
  
  ! look up values to define the scaling of photosynthesis from the leaf to the canopy
  USE mDecisions_module,only:  &
   constantScaling,            & ! constant scaling factor
   laiScaling                    ! exponential function of LAI (Leuning, Plant Cell Env 1995: "Scaling from..." [eq 9])
  
  ! privacy
  implicit none
  private
  public::stomResist_kernel2
  ! spatial indices
  integer(i4b),parameter :: iLoc = 1   ! i-location
  integer(i4b),parameter :: jLoc = 1   ! j-location
  ! conversion factors
  real(rkind),parameter     :: joule2umolConv=4.6_rkind   ! conversion factor from joules to umol photons (umol J-1)
  ! algorithmic parameters
  real(rkind),parameter     :: missingValue=-9999._rkind  ! missing value, used when diagnostic or state variables are undefined
  real(rkind),parameter     :: mpe=1.e-6_rkind            ! prevents overflow error if division by zero
  real(rkind),parameter     :: dx=1.e-6_rkind             ! finite difference increment
  
  contains
  
  
  !  ! ************************************************************************************************
  !  ! public subroutine stomResist: compute stomatal resistance
  !  ! ************************************************************************************************
  !  subroutine stomResist(&
  !                        ! input: state and diagnostic variables
  !                        scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
  !                        scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
  !                        scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
  !                        ! input: data structures
  !                        type_data,                           & ! intent(in):    type of vegetation and soil
  !                        forc_data,                           & ! intent(in):    model forcing data
  !                        mpar_data,                           & ! intent(in):    model parameters
  !                        model_decisions,                     & ! intent(in):    model decisions
  !                        ! input-output: data structures
  !                        diag_data,                           & ! intent(inout): model diagnostic variables for a local HRU
  !                        flux_data,                           & ! intent(inout): model fluxes for a local HRU
  !                        ! output: error control
  !                        err,message)                           ! intent(out): error control
  !  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  !  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  !  ! conversion functions
  !  USE conv_funcs_module,only:satVapPress   ! function to compute the saturated vapor pressure (Pa)
  !  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  !  ! input: state and diagnostic variables
  !  real(rkind),intent(in)             :: scalarVegetationTemp      ! vegetation temperature (K)
  !  real(rkind),intent(in)             :: scalarSatVP_VegTemp       ! saturation vapor pressure at vegetation temperature (Pa)
  !  real(rkind),intent(in)             :: scalarVP_CanopyAir        ! canopy air vapor pressure (Pa)
  !  ! input: data structures
  !  type(var_i),intent(in)          :: type_data                    ! type of vegetation and soil
  !  type(var_d),intent(in)          :: forc_data                    ! model forcing data
  !  type(var_dlength),intent(in)    :: mpar_data                    ! model parameters
  !  type(model_options),intent(in)  :: model_decisions(:)           ! model decisions
  !  ! input-output: data structures
  !  type(var_dlength),intent(inout) :: diag_data                    ! diagnostic variables for a local HRU
  !  type(var_dlength),intent(inout) :: flux_data                    ! model fluxes for a local HRU
  !  ! output: error control
  !  integer(i4b),intent(out)        :: err                          ! error code
  !  character(*),intent(out)        :: message                      ! error message
  !  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  !  ! local variables
  !  character(LEN=256)              :: cmessage                     ! error message of downwind routine
  !  integer(i4b),parameter          :: ixSunlit=1                   ! named variable for sunlit leaves
  !  integer(i4b),parameter          :: ixShaded=2                   ! named variable for shaded leaves
  !  integer(i4b)                    :: iSunShade                    ! index defining sunlit or shaded leaves
  !  real(rkind)                        :: absorbedPAR               ! absorbed PAR (W m-2)
  !  real(rkind)                        :: scalarStomResist          ! stomatal resistance (s m-1)
  !  real(rkind)                        :: scalarPhotosynthesis      ! photosynthesis (umol CO2 m-2 s-1)
  !  real(rkind)                        :: ci                        ! intercellular co2 partial pressure (Pa)
  !  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  
  !  ! associate variables in the data structure
  !  associate(&
  
  !  ! input: model decisions
  !  ix_stomResist                   => model_decisions(iLookDECISIONS%stomResist)%iDecision,           & ! intent(in): [i4b] choice of function for stomatal resistance
  
  !  ! input: physical attributes
  !  vegTypeIndex                    => type_data%var(iLookTYPE%vegTypeIndex),                          & ! intent(in): [i4b] vegetation type index
  !  minStomatalResistance           => mpar_data%var(iLookPARAM%minStomatalResistance)%dat(1),         & ! intent(in): [dp] mimimum stomatal resistance (s m-1)
  !  vcmax25_canopyTop               => mpar_data%var(iLookPARAM%vcmax25_canopyTop)%dat(1),             & ! intent(in): [dp] potential carboxylation rate at 25 degrees C at the canopy top (umol co2 m-2 s-1)
  
  !  ! input: forcing at the upper boundary
  !  airtemp                         => forc_data%var(iLookFORCE%airtemp),                              & ! intent(in): [dp] air temperature at some height above the surface (K)
  !  airpres                         => forc_data%var(iLookFORCE%airpres),                              & ! intent(in): [dp] air pressure at some height above the surface (Pa)
  !  scalarO2air                     => diag_data%var(iLookDIAG%scalarO2air)%dat(1),                    & ! intent(in): [dp] atmospheric o2 concentration (Pa)
  !  scalarCO2air                    => diag_data%var(iLookDIAG%scalarCO2air)%dat(1),                   & ! intent(in): [dp] atmospheric co2 concentration (Pa)
  !  scalarCanopySunlitPAR           => flux_data%var(iLookFLUX%scalarCanopySunlitPAR)%dat(1),          & ! intent(in): [dp] average absorbed par for sunlit leaves (w m-2)
  !  scalarCanopyShadedPAR           => flux_data%var(iLookFLUX%scalarCanopyShadedPAR)%dat(1),          & ! intent(in): [dp] average absorbed par for shaded leaves (w m-2)
  
  !  ! input: state and diagnostic variables
  !  scalarGrowingSeasonIndex        => diag_data%var(iLookDIAG%scalarGrowingSeasonIndex)%dat(1),       & ! intent(in): [dp] growing season index (0=off, 1=on)
  !  scalarFoliageNitrogenFactor     => diag_data%var(iLookDIAG%scalarFoliageNitrogenFactor)%dat(1),    & ! intent(in): [dp] foliage nitrogen concentration (1.0 = saturated)
  !  scalarTranspireLim              => diag_data%var(iLookDIAG%scalarTranspireLim)%dat(1),             & ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
  !  scalarLeafResistance            => flux_data%var(iLookFLUX%scalarLeafResistance)%dat(1),           & ! intent(in): [dp] mean leaf boundary layer resistance per unit leaf area (s m-1)
  
  !  ! output: stomatal resistance and photosynthesis
  !  scalarStomResistSunlit          => flux_data%var(iLookFLUX%scalarStomResistSunlit)%dat(1),         & ! intent(out): [dp] stomatal resistance for sunlit leaves (s m-1)
  !  scalarStomResistShaded          => flux_data%var(iLookFLUX%scalarStomResistShaded)%dat(1),         & ! intent(out): [dp] stomatal resistance for shaded leaves (s m-1)
  !  scalarPhotosynthesisSunlit      => flux_data%var(iLookFLUX%scalarPhotosynthesisSunlit)%dat(1),     & ! intent(out): [dp] sunlit photosynthesis (umolco2 m-2 s-1)
  !  scalarPhotosynthesisShaded      => flux_data%var(iLookFLUX%scalarPhotosynthesisShaded)%dat(1),     & ! intent(out): [dp] shaded photosynthesis (umolco2 m-2 s-1)
  
  !  ! output: carbon dioxide partial pressure of leaf interior (sunlit leaves) (Pa)
  !  scalarIntercellularCO2Sunlit    => diag_data%var(iLookDIAG%scalarIntercellularCO2Sunlit)%dat(1),   & ! intent(out): [dp] carbon dioxide partial pressure of leaf interior (sunlit leaves) (Pa)
  !  scalarIntercellularCO2Shaded    => diag_data%var(iLookDIAG%scalarIntercellularCO2Shaded)%dat(1)    & ! intent(out): [dp] carbon dioxide partial pressure of leaf interior (shaded leaves) (Pa)
  
  !  )
  !  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  !  ! initialize error control
  !  err=0; message="stomResist/"
  
  !  ! identify option for stomatal resistance
  !  select case(ix_stomResist)
  
  !   ! *******************************************************************************************************************************************
  
  !   ! simple resistance formulation
  !   case(simpleResistance)
  !    ! check that we don't divide by zero -- should be set to minimum of tiny in routine soilResist
  !    if(scalarTranspireLim < tiny(airpres))then; err=20; message=trim(message)//'soil moisture stress factor is < tiny -- this will cause problems'; return; end if
  !    ! compute stomatal resistance (assume equal for sunlit and shaded leaves)
  !    scalarStomResistSunlit = minStomatalResistance/scalarTranspireLim
  !    scalarStomResistShaded = scalarStomResistSunlit
  !    ! set photosynthesis to missing (not computed)
  !    scalarPhotosynthesisSunlit = missingValue
  !    scalarPhotosynthesisShaded = missingValue
  
  !   ! *******************************************************************************************************************************************
  
  !   ! flexible Ball-Berry
  !   case(BallBerryFlex)
  
  !    ! loop through sunlit and shaded leaves
  !    do iSunShade=1,2
  
  !     ! get appropriate input values
  !     select case(iSunShade)
  !      ! sunlit leaves
  !      case(ixSunlit)
  !       absorbedPAR = scalarCanopySunlitPAR         ! average absorbed par for sunlit leaves (w m-2)
  !       ci          = scalarIntercellularCO2Sunlit  ! co2 of the leaf interior for sunlit leaves (Pa)
  !      ! shaded leaves
  !      case(ixShaded)
  !       absorbedPAR = scalarCanopyShadedPAR         ! average absorbed par for shaded leaves (w m-2)
  !       ci          = scalarIntercellularCO2Shaded  ! co2 of the leaf interior for shaded leaves (Pa)
  !      ! check
  !      case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
  !     end select
  
  !     ! compute photosynthesis and stomatal resistance
  !     call stomResist_flex(&
  !                          ! input: state and diagnostic variables
  !                          scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
  !                          scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
  !                          scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
  !                          absorbedPAR,                         & ! intent(in): absorbed PAR (W m-2)
  !                          ! input: data structures
  !                          forc_data,                           & ! intent(in): model forcing data
  !                          mpar_data,                           & ! intent(in): model parameters
  !                          diag_data,                           & ! intent(in): model diagnostic variables for a local HRU
  !                          flux_data,                           & ! intent(in): model fluxes for a local HRU
  !                          model_decisions,                     & ! intent(in): model decisions
  !                          ! input-output
  !                          ci,                                  & ! intent(inout): co2 of the leaf interior (Pa)
  !                          ! output:
  !                          scalarStomResist,                    & ! intent(out): stomatal resistance (s m-1)
  !                          scalarPhotosynthesis,                & ! intent(out): photosynthesis (umol CO2 m-2 s-1)
  !                          ! output: error control
  !                          err,cmessage)                          ! intent(out): error control
  !     if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  
  !     ! assign output variables
  !     select case(iSunShade)
  !      case(ixSunlit)
  !       scalarStomResistSunlit       = scalarStomResist
  !       scalarPhotosynthesisSunlit   = scalarPhotosynthesis
  !       scalarIntercellularCO2Sunlit = ci
  !      case(ixShaded)
  !       scalarStomResistShaded       = scalarStomResist
  !       scalarPhotosynthesisShaded   = scalarPhotosynthesis
  !       scalarIntercellularCO2Shaded = ci
  !      case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
  !     end select
  
  !     ! print progress
  !     !write(*,'(a,1x,20(f12.5,1x))') 'leafTemp, par, psn, rs = ', scalarVegetationTemp, absorbedPAR, scalarPhotosynthesis, scalarStomResist
  
  !    end do  ! looping through sunlit and shaded leaves
  
  
  !   ! *******************************************************************************************************************************************
  !   ! compute stomatal resistance (wrapper around the Noah-MP routines)
  !   ! NOTE: canopy air vapor pressure is from the previous time step
  !   case(BallBerry,Jarvis)
  !    call stomResist_NoahMP(&
  !                           ! input (model decisions)
  !                           ix_stomResist,                     & ! intent(in): choice of function for stomatal resistance
  !                           ! input (local attributes)
  !                           vegTypeIndex,                      & ! intent(in): vegetation type index
  !                           iLoc, jLoc,                        & ! intent(in): spatial location indices
  !                           ! input (forcing)
  !                           airtemp,                           & ! intent(in): air temperature at some height above the surface (K)
  !                           airpres,                           & ! intent(in): air pressure at some height above the surface (Pa)
  !                           scalarO2air,                       & ! intent(in): atmospheric o2 concentration (Pa)
  !                           scalarCO2air,                      & ! intent(in): atmospheric co2 concentration (Pa)
  !                           scalarCanopySunlitPAR,             & ! intent(in): average absorbed par for sunlit leaves (w m-2)
  !                           scalarCanopyShadedPAR,             & ! intent(in): average absorbed par for shaded leaves (w m-2)
  !                           ! input (state and diagnostic variables)
  !                           scalarGrowingSeasonIndex,          & ! intent(in): growing season index (0=off, 1=on)
  !                           scalarFoliageNitrogenFactor,       & ! intent(in): foliage nitrogen concentration (1=saturated)
  !                           scalarTranspireLim,                & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
  !                           scalarLeafResistance,              & ! intent(in): leaf boundary layer resistance (s m-1)
  !                           scalarVegetationTemp,              & ! intent(in): temperature of the vegetation canopy (K)
  !                           scalarSatVP_VegTemp,               & ! intent(in): saturation vapor pressure at the temperature of the veg canopy (Pa)
  !                           scalarVP_CanopyAir,                & ! intent(in): canopy air vapor pressure (Pa)
  !                           ! output
  !                           scalarStomResistSunlit,            & ! intent(out): stomatal resistance for sunlit leaves (s m-1)
  !                           scalarStomResistShaded,            & ! intent(out): stomatal resistance for shaded leaves (s m-1)
  !                           scalarPhotosynthesisSunlit,        & ! intent(out): sunlit photosynthesis (umolco2 m-2 s-1)
  !                           scalarPhotosynthesisShaded,        & ! intent(out): shaded photosynthesis (umolco2 m-2 s-1)
  !                           err,cmessage                       ) ! intent(out): error control
  !    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  
  !   ! *******************************************************************************************************************************************
  
  !   ! error check
  !   case default; err=20; message=trim(message)//'unable to identify option for stomatal resistance'; return
  
  !   ! *******************************************************************************************************************************************
  
  !  end select  ! (identifying option for stomatal resistance)
  
  !  ! print progress
  !  !write(*,'(a,1x,L1,1x,20(f16.8,1x))') 'ix_StomResist==BallBerryFlex, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded, scalarStomResistSunlit, scalarPhotosynthesisShaded = ', &
  !  !                                      ix_StomResist==BallBerryFlex, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded, scalarStomResistSunlit, scalarPhotosynthesisShaded
  !  !pause
  
  !  ! end association to variables in the data structures
  !  end associate
  
  !  end subroutine stomResist
  
  
  !  subroutine stomResist_device(&
  !                        ! input: state and diagnostic variables
  !                        scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
  !                        scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
  !                        scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
  !                        ! input: data structures
  !                        type_data,                           & ! intent(in):    type of vegetation and soil
  !                        forc_data,                           & ! intent(in):    model forcing data
  !                        mpar_data,                           & ! intent(in):    model parameters
  !                        model_decisions,                     & ! intent(in):    model decisions
  !                        ! input-output: data structures
  !                        diag_data,                           & ! intent(inout): model diagnostic variables for a local HRU
  !                        flux_data,                           & ! intent(inout): model fluxes for a local HRU
  !                        nGRU, &
  !                        ! output: error control
  !                        err,message)                           ! intent(out): error control
  !  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  !  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  !  ! conversion functions
  !  USE conv_funcs_module,only:satVapPress   ! function to compute the saturated vapor pressure (Pa)
  !  use type4ida,only:data4ida
  !  use device_data_types
  !  use noahmp_globals
  !  USE NOAHMP_VEG_PARAMETERS
  
  
  !  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  !  ! input: state and diagnostic variables
  !  real(rkind),intent(in),device             :: scalarVegetationTemp(:)      ! vegetation temperature (K)
  !  real(rkind),intent(in),device             :: scalarSatVP_VegTemp(:)       ! saturation vapor pressure at vegetation temperature (Pa)
  !  real(rkind),intent(in),device             :: scalarVP_CanopyAir(:)        ! canopy air vapor pressure (Pa)
  !  ! input: data structures
  !  type(var_i),intent(in)          :: type_data                    ! type of vegetation and soil
  !  type(forc_data_device),intent(in)          :: forc_data                    ! model forcing data
  !  type(mpar_data_device),intent(in)    :: mpar_data                    ! model parameters
  !  type(model_options),intent(in)  :: model_decisions(:)           ! model decisions
  !  ! input-output: data structures
  !  type(diag_data_device),intent(inout) :: diag_data                    ! diagnostic variables for a local HRU
  !  type(flux_data_device),intent(inout) :: flux_data                    ! model fluxes for a local HRU
  !  integer(i4b),intent(in) :: nGRU
  !  ! output: error control
  !  integer(i4b),intent(out)        :: err                          ! error code
  !  character(*),intent(out)        :: message                      ! error message
  !  ! -----------------------------------------------------------------------------------------------------------------------------------------------------
  !  ! local variables
  !  character(LEN=256)              :: cmessage                     ! error message of downwind routine
  !  integer(i4b),parameter          :: ixSunlit=1                   ! named variable for sunlit leaves
  !  integer(i4b),parameter          :: ixShaded=2                   ! named variable for shaded leaves
  !  integer(i4b)                    :: iSunShade                    ! index defining sunlit or shaded leaves
  ! !  real(rkind),device,dimension(nGRU,2)                        :: absorbedPAR               ! absorbed PAR (W m-2)
  ! !  real(rkind),device,dimension(nGRU,2)                        :: scalarStomResist          ! stomatal resistance (s m-1)
  ! !  real(rkind),device,dimension(nGRU,2)                        :: scalarPhotosynthesis      ! photosynthesis (umol CO2 m-2 s-1)
  ! !  real(rkind),device,dimension(nGRU,2)                        :: ci                        ! intercellular co2 partial pressure (Pa)
  !  integer(i4b) :: iGRU
  !  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  
  !  ! associate variables in the data structure
  !  associate(&
  
  !  ! input: model decisions
  !  ix_stomResist                   => model_decisions(iLookDECISIONS%stomResist)%iDecision,           & ! intent(in): [i4b] choice of function for stomatal resistance
  
  !  ! input: physical attributes
  !  vegTypeIndex                    => type_data%var(iLookTYPE%vegTypeIndex),                          & ! intent(in): [i4b] vegetation type index
  !  minStomatalResistance           => mpar_data%minStomatalResistance,         & ! intent(in): [dp] mimimum stomatal resistance (s m-1)
  !  vcmax25_canopyTop               => mpar_data%vcmax25_canopyTop,             & ! intent(in): [dp] potential carboxylation rate at 25 degrees C at the canopy top (umol co2 m-2 s-1)
  
  !  ! input: forcing at the upper boundary
  !  airtemp                         => forc_data%airtemp,                              & ! intent(in): [dp] air temperature at some height above the surface (K)
  !  airpres                         => forc_data%airpres,                              & ! intent(in): [dp] air pressure at some height above the surface (Pa)
  !  scalarO2air                     => diag_data%scalarO2air,                    & ! intent(in): [dp] atmospheric o2 concentration (Pa)
  !  scalarCO2air                    => diag_data%scalarCO2air,                   & ! intent(in): [dp] atmospheric co2 concentration (Pa)
  !  scalarCanopySunlitPAR           => flux_data%scalarCanopySunlitPAR,          & ! intent(in): [dp] average absorbed par for sunlit leaves (w m-2)
  !  scalarCanopyShadedPAR           => flux_data%scalarCanopyShadedPAR,          & ! intent(in): [dp] average absorbed par for shaded leaves (w m-2)
  
  !  ! input: state and diagnostic variables
  !  scalarGrowingSeasonIndex        => diag_data%scalarGrowingSeasonIndex,       & ! intent(in): [dp] growing season index (0=off, 1=on)
  !  scalarFoliageNitrogenFactor     => diag_data%scalarFoliageNitrogenFactor,    & ! intent(in): [dp] foliage nitrogen concentration (1.0 = saturated)
  !  scalarTranspireLim              => diag_data%scalarTranspireLim,             & ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
  !  scalarLeafResistance            => flux_data%scalarLeafResistance,           & ! intent(in): [dp] mean leaf boundary layer resistance per unit leaf area (s m-1)
  
  !  ! output: stomatal resistance and photosynthesis
  !  scalarStomResistSunlit          => flux_data%scalarStomResistSunlit,         & ! intent(out): [dp] stomatal resistance for sunlit leaves (s m-1)
  !  scalarStomResistShaded          => flux_data%scalarStomResistShaded,         & ! intent(out): [dp] stomatal resistance for shaded leaves (s m-1)
  !  scalarPhotosynthesisSunlit      => flux_data%scalarPhotosynthesisSunlit,     & ! intent(out): [dp] sunlit photosynthesis (umolco2 m-2 s-1)
  !  scalarPhotosynthesisShaded      => flux_data%scalarPhotosynthesisShaded,     & ! intent(out): [dp] shaded photosynthesis (umolco2 m-2 s-1)
  
  !  ! output: carbon dioxide partial pressure of leaf interior (sunlit leaves) (Pa)
  !  scalarIntercellularCO2Sunlit    => diag_data%scalarIntercellularCO2Sunlit,   & ! intent(out): [dp] carbon dioxide partial pressure of leaf interior (sunlit leaves) (Pa)
  !  scalarIntercellularCO2Shaded    => diag_data%scalarIntercellularCO2Shaded    & ! intent(out): [dp] carbon dioxide partial pressure of leaf interior (shaded leaves) (Pa)
  
  !  )
  !  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  !  ! initialize error control
  !  err=0; message="stomResist/"
  
  !  do iGRU=1,nGRU
  !   call stomResist_kernel<<<1,1>>>(scalarTranspireLim(iGRU), airpres, scalarStomResistSunlit(iGRU), scalarStomResistShaded(iGRU), scalarPhotosynthesisSunlit(iGRU), scalarPhotosynthesisShaded(iGRU),&
  !   minStomatalResistance, &
  !   scalarCanopySunlitPAR(iGRU), scalarIntercellularCO2Sunlit(iGRU),&
  !   scalarCanopyShadedPAR(iGRU), scalarIntercellularCO2Shaded(iGRU), &
  !   scalarVegetationTemp(iGRU),                & ! intent(in): vegetation temperature (K)
  !   scalarSatVP_VegTemp(iGRU),                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
  !   scalarVP_CanopyAir(iGRU),                  & ! intent(in): canopy air vapor pressure (Pa)
  !   airtemp, &
  !   mpar_data%Kc25, mpar_data%Ko25, mpar_data%Kc_qFac, mpar_data%Ko_qFac, &
  !   mpar_data%kc_Ha, mpar_data%ko_Ha, &
  !   mpar_data%vcmax25_canopyTop, mpar_data%vcmax_qFac, mpar_data%vcmax_Ha, mpar_data%vcmax_Hd, mpar_data%vcmax_Sv, mpar_data%vcmax_Kn, &
  !   mpar_data%jmax25_scale, mpar_data%jmax_Ha, mpar_data%jmax_Hd, mpar_data%jmax_Sv, &
  !   mpar_data%fractionJ, mpar_data%quantamYield, mpar_data%vpScaleFactor, mpar_data%cond2photo_slope, mpar_data%minStomatalConductance, &
  !   scalarO2air(iGRU), scalarCO2air(iGRU), &
  !   diag_data%scalarExposedLAI(iGRU), diag_data%scalarGrowingSeasonIndex(iGRU), diag_data%scalarFoliageNitrogenFactor(iGRU), &
  !   flux_data%scalarLeafResistance(iGRU), &
  !   model_decisions(iLookDECISIONS%bbTempFunc)%iDecision, &
  !   model_decisions(iLookDECISIONS%bbHumdFunc)%iDecision, &
  !   model_decisions(iLookDECISIONS%bbElecFunc)%iDecision, &
  !   model_decisions(iLookDECISIONS%bbCO2point)%iDecision, &
  !   model_decisions(iLookDECISIONS%bbNumerics)%iDecision, &
  !   model_decisions(iLookDECISIONS%bbAssimFnc)%iDecision, &
  !   model_decisions(iLookDECISIONS%bbCanIntg8)%iDecision, &
  !   ix_stomResist, &
  !   vegTypeIndex,                      & ! intent(in): vegetation type index
  !   iLoc, jLoc,                        & ! intent(in): spatial location indices
  !   bp(vegTypeIndex), mp(vegTypeIndex)  ,c3psn(vegTypeIndex), avcmx(vegTypeIndex), vcmx25(vegTypeIndex), ako(vegTypeIndex), ko25(vegTypeIndex), akc(vegTypeIndex),kc25(vegTypeIndex),qe25(vegTypeIndex),folnmx(vegTypeIndex),&
  !   rgl, rsmin, rsmax, topt, hs &
  ! )
  !  enddo
  ! !  ! identify option for stomatal resistance
  ! !  select case(ix_stomResist)
  
  ! !   ! *******************************************************************************************************************************************
  
  ! !   ! simple resistance formulation
  ! !   case(simpleResistance)
  ! !     !$cuf kernel do(1) <<<*,*>>>
  ! !     do iGRU=1,nGRU
  ! !     call simpleResistance_kernel(scalarTranspireLim(iGRU), airpres, scalarStomResistSunlit(iGRU), scalarStomResistShaded(iGRU), scalarPhotosynthesisSunlit(iGRU), scalarPhotosynthesisShaded(iGRU),&
  ! !     minStomatalResistance)
  ! !     enddo
  
  ! !   ! *******************************************************************************************************************************************
  
  ! !   ! flexible Ball-Berry
  ! !   case(BallBerryFlex)
  
  ! !    ! loop through sunlit and shaded leaves
  ! !     do iGRU=1,nGRU
  
  ! !     call ballberryflex_kernel<<<1,1>>>(&
  ! !                          scalarCanopySunlitPAR(iGRU), scalarIntercellularCO2Sunlit(iGRU),&
  ! !                          scalarCanopyShadedPAR(iGRU), scalarIntercellularCO2Shaded(iGRU), &
  ! !                          scalarStomResistSunlit(iGRU), scalarPhotosynthesisSunlit(iGRU), &
  ! !                          scalarStomResistShaded(iGRU), scalarPhotosynthesisShaded(iGRU), &
  ! !                          ! input: state and diagnostic variables
  ! !                          scalarVegetationTemp(iGRU),                & ! intent(in): vegetation temperature (K)
  ! !                          scalarSatVP_VegTemp(iGRU),                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
  ! !                          scalarVP_CanopyAir(iGRU),                  & ! intent(in): canopy air vapor pressure (Pa)
  ! !                         !  absorbedPAR(iGRU,iSunShade),                         & ! intent(in): absorbed PAR (W m-2)
  ! !                          ! input: data structures
  ! !                         !  forc_data,                           & ! intent(in): model forcing data
  ! !                          forc_data%airtemp, forc_data%airpres, &
  ! !                         !  mpar_data,                           & ! intent(in): model parameters
  ! !                          mpar_data%Kc25, mpar_data%Ko25, mpar_data%Kc_qFac, mpar_data%Ko_qFac, &
  ! !                          mpar_data%kc_Ha, mpar_data%ko_Ha, &
  ! !                          mpar_data%vcmax25_canopyTop, mpar_data%vcmax_qFac, mpar_data%vcmax_Ha, mpar_data%vcmax_Hd, mpar_data%vcmax_Sv, mpar_data%vcmax_Kn, &
  ! !                          mpar_data%jmax25_scale, mpar_data%jmax_Ha, mpar_data%jmax_Hd, mpar_data%jmax_Sv, &
  ! !                          mpar_data%fractionJ, mpar_data%quantamYield, mpar_data%vpScaleFactor, mpar_data%cond2photo_slope, mpar_data%minStomatalConductance, &
  ! !                         !  diag_data,                           & ! intent(in): model diagnostic variables for a local HRU
  ! !                          diag_data%scalarO2air(iGRU), diag_data%scalarCO2air(iGRU), &
  ! !                          diag_data%scalarExposedLAI(iGRU), diag_data%scalarGrowingSeasonIndex(iGRU), diag_data%scalarFoliageNitrogenFactor(iGRU), diag_data%scalarTranspireLim(iGRU), &
  ! !                         !  flux_data,                           & ! intent(in): model fluxes for a local HRU
  ! !                          flux_data%scalarLeafResistance(iGRU), &
  ! !                         !  model_decisions,                     & ! intent(in): model decisions
  ! !                          model_decisions(iLookDECISIONS%bbTempFunc)%iDecision, &
  ! !                          model_decisions(iLookDECISIONS%bbHumdFunc)%iDecision, &
  ! !                          model_decisions(iLookDECISIONS%bbElecFunc)%iDecision, &
  ! !                          model_decisions(iLookDECISIONS%bbCO2point)%iDecision, &
  ! !                          model_decisions(iLookDECISIONS%bbNumerics)%iDecision, &
  ! !                          model_decisions(iLookDECISIONS%bbAssimFnc)%iDecision, &
  ! !                          model_decisions(iLookDECISIONS%bbCanIntg8)%iDecision &
  ! !                          ! input-output
  ! !                         !  ci(iGRU,iSunShade),                                  & ! intent(inout): co2 of the leaf interior (Pa)
  ! !                          ! output:
  ! !                         !  scalarStomResist(iGRU,iSunShade),                    & ! intent(out): stomatal resistance (s m-1)
  ! !                         !  scalarPhotosynthesis(iGRU,iSunShade)
  ! !                          )!,                & ! intent(out): photosynthesis (umol CO2 m-2 s-1)
  ! !                          ! output: error control
  ! !                         !  err,cmessage)                          ! intent(out): error control
  ! !     ! if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  
  ! !     ! print progress
  ! !     !write(*,'(a,1x,20(f12.5,1x))') 'leafTemp, par, psn, rs = ', scalarVegetationTemp, absorbedPAR, scalarPhotosynthesis, scalarStomResist
  
     
  ! !   enddo
  
  
  ! !   ! *******************************************************************************************************************************************
  ! !   ! compute stomatal resistance (wrapper around the Noah-MP routines)
  ! !   ! NOTE: canopy air vapor pressure is from the previous time step
  ! !   case(BallBerry,Jarvis)
  ! !     do iGRU=1,nGRU
  ! !    call stomResist_NoahMP_device<<<1,1>>>(&
  ! !                           ! input (model decisions)
  ! !                           ix_stomResist,                     & ! intent(in): choice of function for stomatal resistance
  ! !                           ! input (local attributes)
  ! !                           vegTypeIndex,                      & ! intent(in): vegetation type index
  ! !                           iLoc, jLoc,                        & ! intent(in): spatial location indices
  ! !                           ! input (forcing)
  ! !                           airtemp,                           & ! intent(in): air temperature at some height above the surface (K)
  ! !                           airpres,                           & ! intent(in): air pressure at some height above the surface (Pa)
  ! !                           scalarO2air(iGRU),                       & ! intent(in): atmospheric o2 concentration (Pa)
  ! !                           scalarCO2air(iGRU),                      & ! intent(in): atmospheric co2 concentration (Pa)
  ! !                           scalarCanopySunlitPAR(iGRU),             & ! intent(in): average absorbed par for sunlit leaves (w m-2)
  ! !                           scalarCanopyShadedPAR(iGRU),             & ! intent(in): average absorbed par for shaded leaves (w m-2)
  ! !                           ! input (state and diagnostic variables)
  ! !                           scalarGrowingSeasonIndex(iGRU),          & ! intent(in): growing season index (0=off, 1=on)
  ! !                           scalarFoliageNitrogenFactor(iGRU),       & ! intent(in): foliage nitrogen concentration (1=saturated)
  ! !                           scalarTranspireLim(iGRU),                & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
  ! !                           scalarLeafResistance(iGRU),              & ! intent(in): leaf boundary layer resistance (s m-1)
  ! !                           scalarVegetationTemp(iGRU),              & ! intent(in): temperature of the vegetation canopy (K)
  ! !                           scalarSatVP_VegTemp(iGRU),               & ! intent(in): saturation vapor pressure at the temperature of the veg canopy (Pa)
  ! !                           scalarVP_CanopyAir(iGRU),                & ! intent(in): canopy air vapor pressure (Pa)
  ! !                           bp(vegTypeIndex), mp(vegTypeIndex)  ,c3psn(vegTypeIndex), avcmx(vegTypeIndex), vcmx25(vegTypeIndex), ako(vegTypeIndex), ko25(vegTypeIndex), akc(vegTypeIndex),kc25(vegTypeIndex),qe25(vegTypeIndex),folnmx(vegTypeIndex),&
  ! !                           rgl, rsmin, rsmax, topt, hs, &
  ! !                           ! output
  ! !                           scalarStomResistSunlit(iGRU),            & ! intent(out): stomatal resistance for sunlit leaves (s m-1)
  ! !                           scalarStomResistShaded(iGRU),            & ! intent(out): stomatal resistance for shaded leaves (s m-1)
  ! !                           scalarPhotosynthesisSunlit(iGRU),        & ! intent(out): sunlit photosynthesis (umolco2 m-2 s-1)
  ! !                           scalarPhotosynthesisShaded(iGRU))!,        & ! intent(out): shaded photosynthesis (umolco2 m-2 s-1)
  ! !                           !err,cmessage                       ) ! intent(out): error control
  ! !     enddo
  ! !    if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  
  ! !   ! *******************************************************************************************************************************************
  
  ! !   ! error check
  ! !   case default; err=20; message=trim(message)//'unable to identify option for stomatal resistance'; return
  
  ! !   ! *******************************************************************************************************************************************
  
  ! !  end select  ! (identifying option for stomatal resistance)
  
  !  ! print progress
  !  !write(*,'(a,1x,L1,1x,20(f16.8,1x))') 'ix_StomResist==BallBerryFlex, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded, scalarStomResistSunlit, scalarPhotosynthesisShaded = ', &
  !  !                                      ix_StomResist==BallBerryFlex, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded, scalarStomResistSunlit, scalarPhotosynthesisShaded
  !  !pause
  
  !  ! end association to variables in the data structures
  !  end associate
  
  !  end subroutine stomResist_device
  
  
   attributes(device) subroutine simpleResistance_kernel(scalarTranspireLim, airpres, scalarStomResistSunlit, scalarStomResistShaded, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded,&
   minStomatalResistance)
   real(rkind),intent(in) :: scalarTranspireLim, airpres, minStomatalResistance
   real(rkind),intent(inout) :: scalarStomResistSunlit, scalarStomResistShaded, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded
      ! check that we don't divide by zero -- should be set to minimum of tiny in routine soilResist
      if(scalarTranspireLim < tiny(airpres))then; 
      ! err=20; message=trim(message)//'soil moisture stress factor is < tiny -- this will cause problems'
      return; end if
     ! compute stomatal resistance (assume equal for sunlit and shaded leaves)
     scalarStomResistSunlit = minStomatalResistance/scalarTranspireLim
     scalarStomResistShaded = scalarStomResistSunlit
     ! set photosynthesis to missing (not computed)
     scalarPhotosynthesisSunlit = missingValue
     scalarPhotosynthesisShaded = missingValue
  end subroutine
  
  attributes(device) subroutine q10Func_kernel(Kc, Ko, x0, x1, vcmax,Kc25, umol_per_mol, &
  kc_qfac, scalarVegetationTemp, tref, tscale, &
  ko25, ko_qfac, vcmax_qfac, vcmax_hd, vcmax_sv, scalartranspirelim, vcmax25)
  real(rkind),intent(inout) :: Kc, Ko, x0, x1, vcmax
  real(rkind),intent(in) :: Kc25
  real(rkind),intent(in), value :: umol_per_mol, Tref, Tscale
  real(rkind),intent(in) :: kc_qfac, scalarVegetationTemp
  real(rkind),intent(in) :: ko25, ko_qfac, vcmax_qfac, vcmax_hd, vcmax_sv, scalartranspirelim, vcmax25
  
  Kc = ave_slp*(Kc25/umol_per_mol)*q10_device(Kc_qFac,scalarVegetationTemp,Tref,Tscale)  ! umol mol-1 --> mol mol-1 --> Pa
  Ko = ave_slp*Ko25*q10_device(Ko_qFac,scalarVegetationTemp,Tref,Tscale)  ! mol mol-1 --> Pa
  ! compute maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
  x0 = q10_device(vcmax_qFac,scalarVegetationTemp,Tref,Tscale)  ! temperature function
  x1 = fHigh_device(vcmax_Hd,vcmax_Sv,scalarVegetationTemp) ! high temperature inhibition function
  vcmax = scalarTranspireLim*vcmax25*x0/x1
  end subroutine
  
  attributes(device) subroutine Arrhenius_kernel(Kc, Ko, x0, x1, vcmax, &
  airpres, kc25, umol_per_mol, kc_ha, scalarVegetationTemp, &
  Tref, Ko25, ko_Ha, vcmax_Ha, vcmax_Hd, vcmax_Sv, &
  scalarTranspireLim, vcmax25)
  real(rkind),intent(inout) :: Kc, Ko, x0, x1, vcmax
  real(rkind),intent(in), value :: umol_per_mol, Tref
  real(rkind),intent(in) :: airpres, kc25, kc_Ha, scalarVegetationTemp
  real(rkind),intent(in) :: Ko25, ko_Ha, vcmax_Ha, vcmax_Hd, vcmax_Sv
  real(rkind),intent(in) :: scalarTranspireLim, vcmax25
  Kc = airpres*(Kc25/umol_per_mol)*fT_device(kc_Ha,scalarVegetationTemp,Tref)  ! umol mol-1 --> mol mol-1 --> Pa
  Ko = airpres*Ko25*fT_device(ko_Ha,scalarVegetationTemp,Tref)  ! mol mol-1 --> Pa
  ! compute maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
  x0 = fT_device(vcmax_Ha,scalarVegetationTemp,Tref)
  x1 = fHigh_device(vcmax_Hd,vcmax_Sv,Tref) / fHigh_device(vcmax_Hd,vcmax_Sv,scalarVegetationTemp)
  vcmax = scalarTranspireLim*vcmax25*x0*x1
  end subroutine
  
  attributes(device) subroutine max_elec_transport_kernel(x0, x1, jmax, &
  jmax_Ha, scalarVegetationTemp, Tref, &
  jMax_Hd, jmax_Sv, jmax25)
  real(rkind),intent(inout) :: x0, x1, jmax
  real(rkind),intent(in) :: Tref
  real(rkind),intent(in) :: jmax_Ha, scalarVegetationTemp
  real(rkind),intent(in) :: jMax_Hd, jmax_Sv, jmax25
  
  x0 = fT_device(jmax_Ha,scalarVegetationTemp,Tref)
  x1 = fHigh_device(jmax_Hd,jmax_Sv,Tref) / fHigh_device(jmax_Hd,jmax_Sv,scalarVegetationTemp)
  jmax = jmax25*x0*x1
  
  end subroutine
  
  attributes(device) subroutine linear_elec_kernel(Js, quantamYield, absorbedPAR)
  real(rkind),intent(inout) :: Js
  real(rkind),intent(in) :: quantamYield
  real(rkind),intent(in) :: absorbedPAR
  
  Js = quantamYield*joule2umolConv*absorbedPAR
  end subroutine
  
  attributes(device) subroutine linearJmax_kernel(x0, x1, Js, &
  quantamYield, absorbedPAR, jmax)
  real(rkind),intent(inout) :: x0, x1, Js
  real(rkind),intent(in) :: quantamYield
  real(rkind),intent(in) :: absorbedPAR, jmax
  
  x0 = quantamYield*joule2umolConv*absorbedPAR
  x1 = x0*jmax / (x0 + 2.1_rkind*jmax)
  Js = x1/4._rkind ! scaled electron transport
  end subroutine
  
  attributes(device) subroutine scaleVcMax(vcmax25, jmax25,gMin, rlb, Kc, Ko, x0, x1, vcmax, jmax, &
  Js, I_ps2, aQuad, bQuad, cQuad, bSign, xTemp, qQuad, root1, root2, fHum, &
  co2compPt, awb, cp2, ci, cMin, cMax, ci_old, cs, csx, g0, &
  dg0_dc, x2, dci_dc, &
  scalarStomResist, scalarPhotosynthesis, &
  ix_bbCanIntg8, ix_bbTempFunc, ix_bbElecFunc, ix_bbHumdFunc, ix_bbNumerics, maxiter, ix_bbAssimFnc, ix_bbCO2point, maxiter_NoahMP, &
  fnf, umol_per_mol, Tref, Tscale, c_ps2, o2scaleFactor, convToler, &
  h2o_co2__stomPores, h2o_co2__leafbl, &
  vcmax25_canopyTop, vcmax_Kn, scalarExposedLAI, &
  jmax25_scale, scalarTranspireLim, minStomatalConductance, &
  scalarLeafResistance, unitConv, &
  Kc25, Kc_qFac, scalarVegetationTemp, &
  Ko25, Ko_qFac, vcmax_qFac, vcmax_Hd, vcmax_Sv, &
  airpres, kc_Ha, ko_Ha, vcmax_Ha, &
  jmax_Ha, jmax_Hd, jmax_Sv, &
  quantamYield, absorbedPAR, fractionJ, &
  scalarVP_CanopyAir, scalarSatVP_VegTemp, vpScaleFactor, &
  scalarO2air, scalarCO2air, &
  psn, dA_dc, cond2photo_slope, rs, drs_dc, &
  xInc, airtemp, scalarGrowingSeasonIndex)
    real(rkind),intent(inout) :: vcmax25, jmax25, gMin, rlb, Kc, Ko, x0, x1, vcmax, jmax
    real(rkind),intent(inout) :: Js, I_ps2, aQuad, bQuad, cQuad, bSign, xTemp, qQuad, root1, root2, fHum
    real(rkind),intent(inout) :: co2compPt, awb, cp2, ci, cMin, cMax, ci_old, cs, csx, g0
    real(rkind),intent(inout) :: dg0_dc, x2, dci_dc
    real(rkind),intent(inout) :: scalarStomResist, scalarPhotosynthesis
    integer(i4b),intent(in) :: ix_bbCanIntg8, ix_bbTempFunc, ix_bbElecFunc, ix_bbHumdFunc, ix_bbNumerics, maxiter, ix_bbAssimFnc, ix_bbCO2point, maxiter_NoahMP
    real(rkind),intent(in) :: fnf, umol_per_mol, Tref, Tscale, c_ps2, o2scaleFactor, convToler
    real(rkind),intent(in) :: h2o_co2__stomPores,h2o_co2__leafbl
    real(rkind),intent(in) :: vcmax25_canopyTop, vcmax_Kn, scalarExposedLAI
    real(rkind),intent(in) :: jmax25_scale, scalarTranspireLim, minStomatalConductance
    real(rkind),intent(in) :: scalarLeafResistance
    real(rkind),intent(inout) :: unitConv
    real(rkind),intent(in) :: Kc25, Kc_qFac, scalarVegetationTemp
    real(rkind),intent(in) :: Ko25, Ko_qFac, vcmax_qFac, vcmax_Hd, vcmax_Sv
    real(rkind),intent(in) :: airpres, kc_Ha, ko_Ha, vcmax_Ha
    real(rkind),intent(in) :: jmax_Ha, jmax_Hd, jmax_Sv
    real(rkind),intent(in) :: quantamYield, absorbedPAR, fractionJ
    real(rkind),intent(in) :: scalarVP_CanopyAir, scalarSatVP_VegTemp, vpScaleFactor
    real(rkind),intent(in) :: scalarO2air, scalarCO2air
    real(rkind),intent(inout) :: psn, dA_dc, cond2photo_slope, rs, drs_dc !!!
    real(rkind),intent(inout) :: xInc !!!
    real(rkind),intent(in) :: airtemp, scalarGrowingSeasonIndex
  
  
    integer(i4b) :: iter
  
  
    ! *****
    ! * preliminaries...
    ! ******************
  
    ! define unit conversion (m s-1 --> mol m-2 s-1)
    ! NOTE: Rgas   = J K-1 Mol-1 (J = kg m2 s-2); airtemp = K; airpres = Pa (kg m-1 s-2)
    unitConv = airpres/(Rgas*airtemp)  ! mol m-3
  
    ! check there is light available for photosynthesis
    if(absorbedPAR < tiny(absorbedPAR) .or. scalarGrowingSeasonIndex < tiny(absorbedPAR))then
      scalarStomResist     = unitConv*umol_per_mol/(scalarTranspireLim*minStomatalConductance)
      scalarPhotosynthesis = 0._rkind
      ci                   = 0._rkind
    end if
  
  
    ! scale vcmax from the leaves to the canopy
    ! exponential function of LAI is described in Leuning, Plant Cell Env 1995: "Scaling from..." [eq 9]
    select case(ix_bbCanIntg8)
      case(constantScaling); vcmax25 = vcmax25_canopyTop*fnf
      case(laiScaling);      vcmax25 = vcmax25_canopyTop*exp(-vcmax_Kn*scalarExposedLAI)
      ! case default; err=20; message=trim(message)//'unable to identify option to scale lai from the leaves to the canopy'; return
  
    end select
  
    ! compute the maximum electron transport rate at 25 deg C (umol m-2 s-1)
    jmax25 = jmax25_scale * vcmax25
  
    ! define the scaled minimum conductance
    gMin = scalarTranspireLim*minStomatalConductance
  
    ! compute the leaf conductance (umol m-2 s-1)
    rlb = scalarLeafResistance/(umol_per_mol*unitConv)  ! s m-1 --> umol-1 m2 s
  
    ! *****
    ! * compute temperature controls on stomatal conductance...
    ! *********************************************************
  
    ! identify the temperature function
    select case(ix_bbTempFunc)
  
      ! q10 function used in CLM4 and Noah-MP
      case(q10Func)
        ! compute the Michaelis-Menten constants (Pa)
        Kc = ave_slp*(Kc25/umol_per_mol)*q10_device(Kc_qFac,scalarVegetationTemp,Tref,Tscale)  ! umol mol-1 --> mol mol-1 --> Pa
        Ko = ave_slp*Ko25*q10_device(Ko_qFac,scalarVegetationTemp,Tref,Tscale)  ! mol mol-1 --> Pa
        ! compute maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
        x0 = q10_device(vcmax_qFac,scalarVegetationTemp,Tref,Tscale)  ! temperature function
        x1 = fHigh_device(vcmax_Hd,vcmax_Sv,scalarVegetationTemp) ! high temperature inhibition function
        vcmax = scalarTranspireLim*vcmax25*x0/x1
  
      ! Arrhenius function used in CLM5 and Cable
      case(Arrhenius)
        ! compute the Michaelis-Menten constants (Pa)
        Kc = airpres*(Kc25/umol_per_mol)*fT_device(kc_Ha,scalarVegetationTemp,Tref)  ! umol mol-1 --> mol mol-1 --> Pa
        Ko = airpres*Ko25*fT_device(ko_Ha,scalarVegetationTemp,Tref)  ! mol mol-1 --> Pa
        ! compute maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
        x0 = fT_device(vcmax_Ha,scalarVegetationTemp,Tref)
        x1 = fHigh_device(vcmax_Hd,vcmax_Sv,Tref) / fHigh_device(vcmax_Hd,vcmax_Sv,scalarVegetationTemp)
        vcmax = scalarTranspireLim*vcmax25*x0*x1
  
  
      ! check found an appropriate option
      ! case default; err=20; message=trim(message)//'unable to find option for leaf temperature controls on stomatal conductance'; return
  
    end select  ! temperature controls
  
    ! *****
    ! * compute electron transport controls on stomatal conductance...
    ! ****************************************************************
  
    ! compute the maximum electron transport rate (umol electron m-2 s-1)
    x0 = fT_device(jmax_Ha,scalarVegetationTemp,Tref)
    x1 = fHigh_device(jmax_Hd,jmax_Sv,Tref) / fHigh_device(jmax_Hd,jmax_Sv,scalarVegetationTemp)
    jmax = jmax25*x0*x1
  
    ! identify the electron transport function
    select case(ix_bbElecFunc)
  
      ! linear model, as used in CLM4 and Noah-MP
      case(linear)
        Js = quantamYield*joule2umolConv*absorbedPAR
  
      ! linear function of qmax, as used in Cable [Wang et al., Ag Forest Met 1998, eq D5]
      case(linearJmax)
        x0 = quantamYield*joule2umolConv*absorbedPAR
        x1 = x0*jmax / (x0 + 2.1_rkind*jmax)
        Js = x1/4._rkind ! scaled electron transport
  
      ! quadraric function of jmax, as used in CLM5 (Bonan et al., JGR 2011, Table B2)
      case(quadraticJmax)
        !  PAR absorbed by PS2 (umol photon m-2 s-1)
        I_ps2 = 0.5_rkind*(1._rkind - fractionJ) * joule2umolConv*absorbedPAR   ! Farquar (1980), eq 8: PAR absorbed by PS2 (umol photon m-2 s-1)
        ! define coefficients in the quadratic equation
        aQuad = c_ps2            ! quadratic coefficient = cuurvature factor for electron transport
        bQuad = -(I_ps2 + jmax)  ! linear coefficient
        cQuad =  I_ps2 * jmax    ! free term
        ! compute the q term (NOTE: bQuad is always positive)
        bSign = abs(bQuad)/bQuad
        xTemp = bQuad*bQuad - 4._rkind *aQuad*cQuad
        qQuad = -0.5_rkind * (bQuad + bSign*sqrt(xTemp))
        ! compute roots
        root1 = qQuad / aQuad
        root2 = cQuad / qQuad
        ! select minimum root, required to ensure J=0 when par=0
        ! NOTE: Wittig et al. select the first root, which is the max in all cases I tried
        Js = min(root1,root2) / 4._rkind  ! scaled J
  
      ! check found an appropriate option
      ! case default; err=20; message=trim(message)//'unable to find option for electron transport controls on stomatal conductance'; return
  
    end select  ! electron transport controls
  
    ! *****
    ! * define additional controls on stomatal conductance...
    ! ****************************************************************
  
    ! define the humidity function
    select case(ix_bbHumdFunc)
      case(humidLeafSurface); fHum = min( max(0.25_rkind, scalarVP_CanopyAir/scalarSatVP_VegTemp), 1._rkind)
      case(scaledHyperbolic); fHum = (scalarSatVP_VegTemp - scalarVP_CanopyAir)/vpScaleFactor
      ! case default; err=20; message=trim(message)//'unable to identify humidity control on stomatal conductance'; return
    end select
  
  
    ! compute the co2 compensation point (Pa)
    co2compPt = (Kc/Ko)*scalarO2air*o2scaleFactor
  
    ! compute the Michaelis-Menten controls (Pa)
    awb = Kc*(1._rkind + scalarO2air/Ko)
  
    ! compute the additional controls in light-limited assimilation
    cp2 = co2compPt*2._rkind
  
    ! define trial value of intercellular co2 (Pa)
    ! NOTE: only initialize if less than the co2 compensation point; otherwise, initialize with previous value
    if(ix_bbNumerics==newtonRaphson)then
      if(ci < co2compPt) ci = 0.7_rkind*scalarCO2air
    else
      ci = 0.7_rkind*scalarCO2air  ! always initialize if not NR
    end if
   
    ! initialize brackets for the solution
    cMin = 0._rkind
    cMax = scalarCO2air
  
  
    ! ***
    ! iterate
    do iter=1,maxiter
  
      ! reset ci
      ci_old = ci
  
      ! *****
      ! * compute photosynthesis and stomatal resistance...
      ! ***************************************************
  
      ! compute gross photosynthesis [follow Farquar (Planta, 1980), as implemented in CLM4 and Noah-MP]
      call photosynthesis_device(.true., ix_bbAssimFnc, ci, co2compPt, awb, cp2, vcmax, Js, psn, dA_dc)
  
      ! compute co2 concentration at leaf surface (Pa)
      x1 = h2o_co2__leafbl * airpres * rlb  ! Pa / (umol co2 m-2 s-1)
      cs = max(scalarCO2air - (x1 * psn), mpe)   ! Pa (avoid divide by zero)
  
      ! compute control of the compensation point on stomatal conductance
      if(ix_bbCO2point == origBWB)then
        csx = cs
      else
        csx = cs - co2compPt
      end if
  
      ! compute conductance in the absence of humidity
      g0     = cond2photo_slope*airpres*psn/csx
      dg0_dc = cond2photo_slope*airpres*dA_dc*(x1*psn/cs + 1._rkind)/csx
  
      ! use quadratic function to compute stomatal resistance
      call quadResist_device(.true.,ix_bbHumdFunc,rlb,fHum,gMin,g0,dg0_dc,rs,drs_dc)
  
      ! compute intercellular co2 partial pressues (Pa)
      x2 = h2o_co2__stomPores * airpres  ! Pa
      ci = max(cs - x2*psn*rs, 0._rkind)    ! Pa
  
      ! final derivative
      if(ci > tiny(ci))then
        dci_dc = -x1*dA_dc - x2*(psn*drs_dc + rs*dA_dc)
      else
        dci_dc = 0._rkind
      end if
  
      ! ! test derivatives
      ! if(testDerivs)then
      !  func1 = testFunc(ci_old,    cond2photo_slope, airpres, scalarCO2air, ix_bbHumdFunc, ix_bbCO2point, ix_bbAssimFnc)
      !  func2 = testFunc(ci_old+dx, cond2photo_slope, airpres, scalarCO2air, ix_bbHumdFunc, ix_bbCO2point, ix_bbAssimFnc)
      !  write(*,'(a,1x,20(e20.10,1x))') '(func2 - func1)/dx, dci_dc = ', &
      !                                   (func2 - func1)/dx, dci_dc
      ! end if  ! if testing the derivatives
  
      ! *****
      ! * iterative solution...
      ! ***********************
  
      ! case for Noah-MP (use fixed point iteration)
      if(ix_bbNumerics==NoahMPsolution)then
        if(iter==maxiter_NoahMP) exit ! exit after a specified number of iterations (normally 3)
        cycle  ! fixed-point iteration
      end if
  
      ! CLM4 and Noah-MP use fixed point iteration, continuing the next iteration from this point
      ! here we try and improve matters by using the derivatives
  
      ! update the brackets for the solution
      if(ci_old > ci)then
        cMax = ci_old
      else
        cMin = ci_old
      end if
  
      ! compute iteration increment (Pa)
      xInc = (ci - ci_old)/(1._rkind - dci_dc)
  
      ! update
      ci = max(ci_old + xInc, 0._rkind)
  
      ! ensure that we stay within brackets
      if(ci > cMax .or. ci < cMin)then
        ci = 0.5_rkind * (cMin + cMax)
      end if
  
      ! check for convergence
      if(abs(xInc) < convToler) exit
        if(iter==maxIter)then
        !  message=trim(message)//'did not converge in stomatal conductance iteration'
        !  err=20; 
        return
      end if
  
    end do  ! iterating
  
    ! assign output variables
    scalarStomResist     = unitConv*umol_per_mol*rs  ! umol-1 m2 s --> s/m
    scalarPhotosynthesis = psn
  
  
  end subroutine
  
  attributes(device) subroutine quadraticJmax_kernel(I_ps2, aQuad, bQuad, cQuad, bSign, xTemp, qQuad, root1, root2, Js, &
  fractionJ, absorbedPAR, c_ps2, jmax)
  real(rkind),intent(inout) :: I_ps2, aQuad, bQuad, cQuad, bSign, xTemp, qQuad, root1, root2, Js
  real(rkind),intent(in) :: fractionJ, absorbedPAR, jmax
  real(rkind),intent(in) :: c_ps2
  
     !  PAR absorbed by PS2 (umol photon m-2 s-1)
  I_ps2 = 0.5_rkind*(1._rkind - fractionJ) * joule2umolConv*absorbedPAR   ! Farquar (1980), eq 8: PAR absorbed by PS2 (umol photon m-2 s-1)
  ! define coefficients in the quadratic equation
  aQuad = c_ps2            ! quadratic coefficient = cuurvature factor for electron transport
  bQuad = -(I_ps2 + jmax)  ! linear coefficient
  cQuad =  I_ps2 * jmax    ! free term
  ! compute the q term (NOTE: bQuad is always positive)
  bSign = abs(bQuad)/bQuad
  xTemp = bQuad*bQuad - 4._rkind *aQuad*cQuad
  qQuad = -0.5_rkind * (bQuad + bSign*sqrt(xTemp))
  ! compute roots
  root1 = qQuad / aQuad
  root2 = cQuad / qQuad
  ! select minimum root, required to ensure J=0 when par=0
  ! NOTE: Wittig et al. select the first root, which is the max in all cases I tried
  Js = min(root1,root2) / 4._rkind  ! scaled J
  end subroutine
  
   ! *******************************************************************************************************
   ! *******************************************************************************************************
   ! *** PRIVATE SUBROUTINES *******************************************************************************
   ! *******************************************************************************************************
   ! *******************************************************************************************************
  
   ! *******************************************************************************************************
   ! private subroutine stomResist_flex: flexible stomatal resistance routine to evaluate different options
   ! *******************************************************************************************************
   subroutine stomResist_flex(&
                              ! input: state and diagnostic variables
                              scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                              scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                              scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                              absorbedPAR,                         & ! intent(in): absorbed PAR (W m-2)
                              ! input: data structures
                              forc_data,                           & ! intent(in): model forcing data
                              mpar_data,                           & ! intent(in): model parameters
                              diag_data,                           & ! intent(in): model diagnostic variables for a local HRU
                              flux_data,                           & ! intent(in): model fluxes for a local HRU
                              model_decisions,                     & ! intent(in): model decisions
                              ! output: stomatal resistance and photosynthesis
                              ci,                                  & ! intent(out): co2 of the leaf interior (Pa)
                              scalarStomResist,                    & ! intent(out): stomatal resistance (s m-1)
                              scalarPhotosynthesis,                & ! intent(out): photosynthesis (umol CO2 m-2 s-1)
                              ! output: error control
                              err,message)                           ! intent(out): error control
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! input: state and diagnostic variables
   real(rkind),intent(in)             :: scalarVegetationTemp          ! vegetation temperature (K)
   real(rkind),intent(in)             :: scalarSatVP_VegTemp           ! saturation vapor pressure at vegetation temperature (Pa)
   real(rkind),intent(in)             :: scalarVP_CanopyAir            ! canopy air vapor pressure (Pa)
   real(rkind),intent(in)             :: absorbedPAR                   ! absorbed PAR (W m-2)
   ! input: data structures
   type(var_d),intent(in)             :: forc_data                     ! model forcing data
   type(var_dlength),intent(in)       :: mpar_data                     ! model parameters
   type(var_dlength),intent(in)       :: diag_data                     ! diagnostic variables for a local HRU
   type(var_dlength),intent(in)       :: flux_data                     ! model fluxes for a local HRU
   type(model_options),intent(in)     :: model_decisions(:)            ! model decisions
   ! output: stomatal resistance and photosynthesis
   real(rkind),intent(inout)          :: ci                            ! intercellular co2 partial pressure (Pa)
   real(rkind),intent(out)            :: scalarStomResist              ! stomatal resistance (s m-1)
   real(rkind),intent(out)            :: scalarPhotosynthesis          ! photosynthesis (umol CO2 m-2 s-1)
   ! output: error control
   integer(i4b),intent(out)           :: err                           ! error code
   character(*),intent(out)           :: message                       ! error message
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! general local variables
   logical(lgt),parameter             :: testDerivs=.false.            ! flag to test the derivatives
   real(rkind)                        :: unitConv                      ! unit conversion factor (mol m-3, convert m s-1 --> mol H20 m-2 s-1)
   real(rkind)                        :: rlb                           ! leaf boundary layer rersistance (umol-1 m2 s)
   real(rkind)                        :: x0,x1,x2                      ! temporary variables
   real(rkind)                        :: co2compPt                     ! co2 compensation point (Pa)
   real(rkind)                        :: fHum                          ! humidity function, fraction [0,1]
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! fixed parameters
   integer(i4b),parameter             :: maxiter=20                    ! maximum number of iterations
   integer(i4b),parameter             :: maxiter_noahMP=3              ! maximum number of iterations for Noah-MP
   real(rkind),parameter              :: convToler=0.0001_rkind        ! convergence tolerance (Pa)
   real(rkind),parameter              :: umol_per_mol=1.e+6_rkind      ! factor to relate umol to mol
   real(rkind),parameter              :: o2scaleFactor=0.105_rkind     ! scaling factor used to compute co2 compesation point (0.21/2)
   real(rkind),parameter              :: h2o_co2__leafbl=1.37_rkind    ! factor to represent the different diffusivities of h2o and co2 in the leaf boundary layer (-)
   real(rkind),parameter              :: h2o_co2__stomPores=1.65_rkind ! factor to represent the different diffusivities of h2o and co2 in the stomatal pores (-)
   real(rkind),parameter              :: Tref=298.16_rkind             ! reference temperature (25 deg C)
   real(rkind),parameter              :: Tscale=10._rkind              ! scaling factor in q10 function (K)
   real(rkind),parameter              :: c_ps2=0.7_rkind               ! curvature factor for electron transport (-)
   real(rkind),parameter              :: fnf=0.6666666667_rkind        ! foliage nitrogen factor (-)
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! photosynthesis
   real(rkind)                        :: Kc,Ko                         ! Michaelis-Menten constants for co2 and o2 (Pa)
   real(rkind)                        :: vcmax25                       ! maximum Rubisco carboxylation rate at 25 deg C (umol m-2 s-1)
   real(rkind)                        :: jmax25                        ! maximum electron transport rate at 25 deg C (umol m-2 s-1)
   real(rkind)                        :: vcmax                         ! maximum Rubisco carboxylation rate (umol m-2 s-1)
   real(rkind)                        :: jmax                          ! maximum electron transport rate (umol m-2 s-1)
   real(rkind)                        :: aQuad                         ! the quadratic coefficient in the quadratic equation
   real(rkind)                        :: bQuad                         ! the linear coefficient in the quadratic equation
   real(rkind)                        :: cQuad                         ! the constant in the quadratic equation
   real(rkind)                        :: bSign                         ! sign of the linear coeffcient
   real(rkind)                        :: xTemp                         ! temporary variable in the quadratic equation
   real(rkind)                        :: qQuad                         ! the "q" term in the quadratic equation
   real(rkind)                        :: root1,root2                   ! roots of the quadratic function
   real(rkind)                        :: Js                            ! scaled electron transport rate (umol co2 m-2 s-1)
   real(rkind)                        :: I_ps2                         ! PAR absorbed by PS2 (umol photon m-2 s-1)
   real(rkind)                        :: awb                           ! Michaelis-Menten control (Pa)
   real(rkind)                        :: cp2                           ! additional controls in light-limited assimilation (Pa)
   real(rkind)                        :: psn                           ! leaf gross photosynthesis rate (umol co2 m-2 s-1)
   real(rkind)                        :: dA_dc                         ! derivative in photosynthesis w.r.t. intercellular co2 concentration
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! stomatal resistance
   real(rkind)                        :: gMin                          ! scaled minimum conductance (umol m-2 s-1)
   real(rkind)                        :: cs                            ! co2 partial pressure at leaf surface (Pa)
   real(rkind)                        :: csx                           ! control of co2 partial pressure at leaf surface on stomatal conductance (Pa)
   real(rkind)                        :: g0                            ! stomatal conductance in the absence of humidity controls (umol m-2 s-1)
   real(rkind)                        :: ci_old                        ! intercellular co2 partial pressure (Pa)
   real(rkind)                        :: rs                            ! stomatal resistance (umol-1 m2 s)
   real(rkind)                        :: dg0_dc                        ! derivative in g0 w.r.t intercellular co2 concentration (umol m-2 s-1 Pa-1)
   real(rkind)                        :: drs_dc                        ! derivative in stomatal resistance w.r.t. intercellular co2 concentration
   real(rkind)                        :: dci_dc                        ! final derivative (-)
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! iterative solution
   real(rkind)                        :: func1,func2                   ! functions for numerical derivative calculation
   real(rkind)                        :: cMin,cMax                     ! solution brackets
   real(rkind)                        :: xInc                          ! iteration increment (Pa)
   integer(i4b)                       :: iter                          ! iteration index
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! associate variables in the data structure
   associate(&
   ! input: model decisions
   ix_bbTempFunc                   => model_decisions(iLookDECISIONS%bbTempFunc)%iDecision,           & ! intent(in): [i4b] leaf temperature controls on photosynthesis + stomatal resistance
   ix_bbHumdFunc                   => model_decisions(iLookDECISIONS%bbHumdFunc)%iDecision,           & ! intent(in): [i4b] humidity controls on stomatal resistance
   ix_bbElecFunc                   => model_decisions(iLookDECISIONS%bbElecFunc)%iDecision,           & ! intent(in): [i4b] dependence of photosynthesis on PAR
   ix_bbCO2point                   => model_decisions(iLookDECISIONS%bbCO2point)%iDecision,           & ! intent(in): [i4b] use of CO2 compensation point to calculate stomatal resistance
   ix_bbNumerics                   => model_decisions(iLookDECISIONS%bbNumerics)%iDecision,           & ! intent(in): [i4b] iterative numerical solution method used in the Ball-Berry parameterization
   ix_bbAssimFnc                   => model_decisions(iLookDECISIONS%bbAssimFnc)%iDecision,           & ! intent(in): [i4b] controls on carbon assimilation (min function, or colimitation)
   ix_bbCanIntg8                   => model_decisions(iLookDECISIONS%bbCanIntg8)%iDecision,           & ! intent(in): [i4b] scaling of photosynthesis from the leaf to the canopy
   ! input: model parameters
   Kc25                            => mpar_data%var(iLookPARAM%Kc25)%dat(1),                          & ! intent(in): [dp] Michaelis-Menten constant for CO2 at 25 degrees C (umol mol-1)
   Ko25                            => mpar_data%var(iLookPARAM%Ko25)%dat(1),                          & ! intent(in): [dp] Michaelis-Menten constant for O2 at 25 degrees C (mol mol-1)
   Kc_qFac                         => mpar_data%var(iLookPARAM%Kc_qFac)%dat(1),                       & ! intent(in): [dp] factor in the q10 function defining temperature controls on Kc (-)
   Ko_qFac                         => mpar_data%var(iLookPARAM%Ko_qFac)%dat(1),                       & ! intent(in): [dp] factor in the q10 function defining temperature controls on Ko (-)
   kc_Ha                           => mpar_data%var(iLookPARAM%kc_Ha)%dat(1),                         & ! intent(in): [dp] activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
   ko_Ha                           => mpar_data%var(iLookPARAM%ko_Ha)%dat(1),                         & ! intent(in): [dp] activation energy for the Michaelis-Menten constant for O2 (J mol-1)
   vcmax25_canopyTop               => mpar_data%var(iLookPARAM%vcmax25_canopyTop)%dat(1),             & ! intent(in): [dp] potential carboxylation rate at 25 degrees C at the canopy top (umol co2 m-2 s-1)
   vcmax_qFac                      => mpar_data%var(iLookPARAM%vcmax_qFac)%dat(1),                    & ! intent(in): [dp] factor in the q10 function defining temperature controls on vcmax (-)
   vcmax_Ha                        => mpar_data%var(iLookPARAM%vcmax_Ha)%dat(1),                      & ! intent(in): [dp] activation energy in the vcmax function (J mol-1)
   vcmax_Hd                        => mpar_data%var(iLookPARAM%vcmax_Hd)%dat(1),                      & ! intent(in): [dp] deactivation energy in the vcmax function (J mol-1)
   vcmax_Sv                        => mpar_data%var(iLookPARAM%vcmax_Sv)%dat(1),                      & ! intent(in): [dp] entropy term in the vcmax function (J mol-1 K-1)
   vcmax_Kn                        => mpar_data%var(iLookPARAM%vcmax_Kn)%dat(1),                      & ! intent(in): [dp] foliage nitrogen decay coefficient (-)
   jmax25_scale                    => mpar_data%var(iLookPARAM%jmax25_scale)%dat(1),                  & ! intent(in): [dp] scaling factor to relate jmax25 to vcmax25 (-)
   jmax_Ha                         => mpar_data%var(iLookPARAM%jmax_Ha)%dat(1),                       & ! intent(in): [dp] activation energy in the jmax function (J mol-1)
   jmax_Hd                         => mpar_data%var(iLookPARAM%jmax_Hd)%dat(1),                       & ! intent(in): [dp] deactivation energy in the jmax function (J mol-1)
   jmax_Sv                         => mpar_data%var(iLookPARAM%jmax_Sv)%dat(1),                       & ! intent(in): [dp] entropy term in the jmax function (J mol-1 K-1)
   fractionJ                       => mpar_data%var(iLookPARAM%fractionJ)%dat(1),                     & ! intent(in): [dp] fraction of light lost by other than the chloroplast lamellae (-)
   quantamYield                    => mpar_data%var(iLookPARAM%quantamYield)%dat(1),                  & ! intent(in): [dp] quantam yield (mol e mol-1 q)
   vpScaleFactor                   => mpar_data%var(iLookPARAM%vpScaleFactor)%dat(1),                 & ! intent(in): [dp] vapor pressure scaling factor in stomatal conductance function (Pa)
   cond2photo_slope                => mpar_data%var(iLookPARAM%cond2photo_slope)%dat(1),              & ! intent(in): [dp] slope of conductance-photosynthesis relationship (-)
   minStomatalConductance          => mpar_data%var(iLookPARAM%minStomatalConductance)%dat(1),        & ! intent(in): [dp] mimimum stomatal conductance (umol H2O m-2 s-1)
   ! input: forcing at the upper boundary
   airtemp                         => forc_data%var(iLookFORCE%airtemp),                              & ! intent(in): [dp] air temperature at some height above the surface (K)
   airpres                         => forc_data%var(iLookFORCE%airpres),                              & ! intent(in): [dp] air pressure at some height above the surface (Pa)
   scalarO2air                     => diag_data%var(iLookDIAG%scalarO2air)%dat(1),                    & ! intent(in): [dp] atmospheric o2 concentration (Pa)
   scalarCO2air                    => diag_data%var(iLookDIAG%scalarCO2air)%dat(1),                   & ! intent(in): [dp] atmospheric co2 concentration (Pa)
   ! input: state and diagnostic variables
   scalarExposedLAI                => diag_data%var(iLookDIAG%scalarExposedLAI)%dat(1),               & ! intent(in): [dp] exposed LAI (m2 m-2)
   scalarGrowingSeasonIndex        => diag_data%var(iLookDIAG%scalarGrowingSeasonIndex)%dat(1),       & ! intent(in): [dp] growing season index (0=off, 1=on)
   scalarFoliageNitrogenFactor     => diag_data%var(iLookDIAG%scalarFoliageNitrogenFactor)%dat(1),    & ! intent(in): [dp] foliage nitrogen concentration (1.0 = saturated)
   scalarTranspireLim              => diag_data%var(iLookDIAG%scalarTranspireLim)%dat(1),             & ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
   scalarLeafResistance            => flux_data%var(iLookFLUX%scalarLeafResistance)%dat(1)            & ! intent(in): [dp] mean leaf boundary layer resistance per unit leaf area (s m-1)
   )
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! initialize error control
   err=0; message="stomResist_flex/"
  
   ! *****
   ! * preliminaries...
   ! ******************
  
   ! define unit conversion (m s-1 --> mol m-2 s-1)
   ! NOTE: Rgas   = J K-1 Mol-1 (J = kg m2 s-2); airtemp = K; airpres = Pa (kg m-1 s-2)
   unitConv = airpres/(Rgas*airtemp)  ! mol m-3
  
   ! check there is light available for photosynthesis
   if(absorbedPAR < tiny(absorbedPAR) .or. scalarGrowingSeasonIndex < tiny(absorbedPAR))then
    scalarStomResist     = unitConv*umol_per_mol/(scalarTranspireLim*minStomatalConductance)
    scalarPhotosynthesis = 0._rkind
    ci                   = 0._rkind
   return
   end if
  
   ! scale vcmax from the leaves to the canopy
   ! exponential function of LAI is described in Leuning, Plant Cell Env 1995: "Scaling from..." [eq 9]
   select case(ix_bbCanIntg8)
    case(constantScaling); vcmax25 = vcmax25_canopyTop*fnf
    case(laiScaling);      vcmax25 = vcmax25_canopyTop*exp(-vcmax_Kn*scalarExposedLAI)
    case default; err=20; message=trim(message)//'unable to identify option to scale lai from the leaves to the canopy'; return
   end select
  
   ! compute the maximum electron transport rate at 25 deg C (umol m-2 s-1)
   jmax25 = jmax25_scale * vcmax25
  
   ! define the scaled minimum conductance
   gMin = scalarTranspireLim*minStomatalConductance
  
   ! compute the leaf conductance (umol m-2 s-1)
   rlb = scalarLeafResistance/(umol_per_mol*unitConv)  ! s m-1 --> umol-1 m2 s
  
   ! *****
   ! * compute temperature controls on stomatal conductance...
   ! *********************************************************
  
   ! identify the temperature function
   select case(ix_bbTempFunc)
  
    ! q10 function used in CLM4 and Noah-MP
    case(q10Func)
     ! compute the Michaelis-Menten constants (Pa)
     Kc = ave_slp*(Kc25/umol_per_mol)*q10(Kc_qFac,scalarVegetationTemp,Tref,Tscale)  ! umol mol-1 --> mol mol-1 --> Pa
     Ko = ave_slp*Ko25*q10(Ko_qFac,scalarVegetationTemp,Tref,Tscale)  ! mol mol-1 --> Pa
     ! compute maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
     x0 = q10(vcmax_qFac,scalarVegetationTemp,Tref,Tscale)  ! temperature function
     x1 = fHigh(vcmax_Hd,vcmax_Sv,scalarVegetationTemp) ! high temperature inhibition function
     vcmax = scalarTranspireLim*vcmax25*x0/x1
  
    ! Arrhenius function used in CLM5 and Cable
    case(Arrhenius)
     ! compute the Michaelis-Menten constants (Pa)
     Kc = airpres*(Kc25/umol_per_mol)*fT(kc_Ha,scalarVegetationTemp,Tref)  ! umol mol-1 --> mol mol-1 --> Pa
     Ko = airpres*Ko25*fT(ko_Ha,scalarVegetationTemp,Tref)  ! mol mol-1 --> Pa
     ! compute maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
     x0 = fT(vcmax_Ha,scalarVegetationTemp,Tref)
     x1 = fhigh(vcmax_Hd,vcmax_Sv,Tref) / fhigh(vcmax_Hd,vcmax_Sv,scalarVegetationTemp)
     vcmax = scalarTranspireLim*vcmax25*x0*x1
  
    ! check found an appropriate option
    case default; err=20; message=trim(message)//'unable to find option for leaf temperature controls on stomatal conductance'; return
  
   end select  ! temperature controls
  
   ! *****
   ! * compute electron transport controls on stomatal conductance...
   ! ****************************************************************
  
   ! compute the maximum electron transport rate (umol electron m-2 s-1)
   x0 = fT(jmax_Ha,scalarVegetationTemp,Tref)
   x1 = fhigh(jmax_Hd,jmax_Sv,Tref) / fhigh(jmax_Hd,jmax_Sv,scalarVegetationTemp)
   jmax = jmax25*x0*x1
  
   ! identify the electron transport function
   select case(ix_bbElecFunc)
  
    ! linear model, as used in CLM4 and Noah-MP
    case(linear)
     Js = quantamYield*joule2umolConv*absorbedPAR
  
    ! linear function of qmax, as used in Cable [Wang et al., Ag Forest Met 1998, eq D5]
    case(linearJmax)
     x0 = quantamYield*joule2umolConv*absorbedPAR
     x1 = x0*jmax / (x0 + 2.1_rkind*jmax)
     Js = x1/4._rkind ! scaled electron transport
  
    ! quadraric function of jmax, as used in CLM5 (Bonan et al., JGR 2011, Table B2)
    case(quadraticJmax)
     !  PAR absorbed by PS2 (umol photon m-2 s-1)
     I_ps2 = 0.5_rkind*(1._rkind - fractionJ) * joule2umolConv*absorbedPAR   ! Farquar (1980), eq 8: PAR absorbed by PS2 (umol photon m-2 s-1)
     ! define coefficients in the quadratic equation
     aQuad = c_ps2            ! quadratic coefficient = cuurvature factor for electron transport
     bQuad = -(I_ps2 + jmax)  ! linear coefficient
     cQuad =  I_ps2 * jmax    ! free term
     ! compute the q term (NOTE: bQuad is always positive)
     bSign = abs(bQuad)/bQuad
     xTemp = bQuad*bQuad - 4._rkind *aQuad*cQuad
     qQuad = -0.5_rkind * (bQuad + bSign*sqrt(xTemp))
     ! compute roots
     root1 = qQuad / aQuad
     root2 = cQuad / qQuad
     ! select minimum root, required to ensure J=0 when par=0
     ! NOTE: Wittig et al. select the first root, which is the max in all cases I tried
     Js = min(root1,root2) / 4._rkind  ! scaled J
  
    ! check found an appropriate option
    case default; err=20; message=trim(message)//'unable to find option for electron transport controls on stomatal conductance'; return
  
   end select  ! electron transport controls
  
   ! *****
   ! * define additional controls on stomatal conductance...
   ! ****************************************************************
  
   ! define the humidity function
   select case(ix_bbHumdFunc)
    case(humidLeafSurface); fHum = min( max(0.25_rkind, scalarVP_CanopyAir/scalarSatVP_VegTemp), 1._rkind)
    case(scaledHyperbolic); fHum = (scalarSatVP_VegTemp - scalarVP_CanopyAir)/vpScaleFactor
    case default; err=20; message=trim(message)//'unable to identify humidity control on stomatal conductance'; return
   end select
  
   ! compute the co2 compensation point (Pa)
   co2compPt = (Kc/Ko)*scalarO2air*o2scaleFactor
  
   ! compute the Michaelis-Menten controls (Pa)
   awb = Kc*(1._rkind + scalarO2air/Ko)
  
   ! compute the additional controls in light-limited assimilation
   cp2 = co2compPt*2._rkind
  
   ! define trial value of intercellular co2 (Pa)
   ! NOTE: only initialize if less than the co2 compensation point; otherwise, initialize with previous value
   if(ix_bbNumerics==newtonRaphson)then
    if(ci < co2compPt) ci = 0.7_rkind*scalarCO2air
   else
    ci = 0.7_rkind*scalarCO2air  ! always initialize if not NR
   end if
   
   ! initialize brackets for the solution
   cMin = 0._rkind
   cMax = scalarCO2air
  
   ! ***
   ! iterate
   do iter=1,maxiter
  
    ! reset ci
    ci_old = ci
  
    ! *****
    ! * compute photosynthesis and stomatal resistance...
    ! ***************************************************
  
    ! compute gross photosynthesis [follow Farquar (Planta, 1980), as implemented in CLM4 and Noah-MP]
    call photosynthesis(.true., ix_bbAssimFnc, ci, co2compPt, awb, cp2, vcmax, Js, psn, dA_dc)
  
    ! compute co2 concentration at leaf surface (Pa)
    x1 = h2o_co2__leafbl * airpres * rlb  ! Pa / (umol co2 m-2 s-1)
    cs = max(scalarCO2air - (x1 * psn), mpe)   ! Pa (avoid divide by zero)
  
    ! compute control of the compensation point on stomatal conductance
    if(ix_bbCO2point == origBWB)then
     csx = cs
    else
     csx = cs - co2compPt
    end if
  
    ! compute conductance in the absence of humidity
    g0     = cond2photo_slope*airpres*psn/csx
    dg0_dc = cond2photo_slope*airpres*dA_dc*(x1*psn/cs + 1._rkind)/csx
  
    ! use quadratic function to compute stomatal resistance
    call quadResist(.true.,ix_bbHumdFunc,rlb,fHum,gMin,g0,dg0_dc,rs,drs_dc)
  
    ! compute intercellular co2 partial pressues (Pa)
    x2 = h2o_co2__stomPores * airpres  ! Pa
    ci = max(cs - x2*psn*rs, 0._rkind)    ! Pa
  
    ! final derivative
    if(ci > tiny(ci))then
     dci_dc = -x1*dA_dc - x2*(psn*drs_dc + rs*dA_dc)
    else
     dci_dc = 0._rkind
    end if
  
    ! test derivatives
    if(testDerivs)then
     func1 = testFunc(ci_old,    cond2photo_slope, airpres, scalarCO2air, ix_bbHumdFunc, ix_bbCO2point, ix_bbAssimFnc)
     func2 = testFunc(ci_old+dx, cond2photo_slope, airpres, scalarCO2air, ix_bbHumdFunc, ix_bbCO2point, ix_bbAssimFnc)
     write(*,'(a,1x,20(e20.10,1x))') '(func2 - func1)/dx, dci_dc = ', &
                                      (func2 - func1)/dx, dci_dc
    end if  ! if testing the derivatives
  
    ! *****
    ! * iterative solution...
    ! ***********************
  
    ! case for Noah-MP (use fixed point iteration)
    if(ix_bbNumerics==NoahMPsolution)then
     if(iter==maxiter_NoahMP) exit ! exit after a specified number of iterations (normally 3)
     cycle  ! fixed-point iteration
    end if
  
    ! CLM4 and Noah-MP use fixed point iteration, continuing the next iteration from this point
    ! here we try and improve matters by using the derivatives
  
    ! update the brackets for the solution
    if(ci_old > ci)then
     cMax = ci_old
    else
     cMin = ci_old
    end if
  
    ! compute iteration increment (Pa)
    xInc = (ci - ci_old)/(1._rkind - dci_dc)
  
    ! update
    ci = max(ci_old + xInc, 0._rkind)
  
    ! ensure that we stay within brackets
    if(ci > cMax .or. ci < cMin)then
     ci = 0.5_rkind * (cMin + cMax)
    end if
  
    ! check for convergence
    if(abs(xInc) < convToler) exit
    if(iter==maxIter)then
     message=trim(message)//'did not converge in stomatal conductance iteration'
     err=20; return
    end if
  
   end do  ! iterating
  
   ! assign output variables
   scalarStomResist     = unitConv*umol_per_mol*rs  ! umol-1 m2 s --> s/m
   scalarPhotosynthesis = psn
  
   end associate
  
   contains
  
    ! ******************************************************
    ! internal function used to test derivatives
    function testFunc(ci, cond2photo_slope, airpres, scalarCO2air, ix_bbHumdFunc, ix_bbCO2point, ix_bbAssimFnc)
    real(rkind),intent(in)     :: ci, cond2photo_slope, airpres, scalarCO2air
    integer(i4b),intent(in) :: ix_bbHumdFunc, ix_bbCO2point, ix_bbAssimFnc
    real(rkind)                :: testFunc
    real(rkind),parameter      :: unUsedInput=0._rkind
    real(rkind)                :: unUsedOutput
  
    ! compute gross photosynthesis [follow Farquar (Planta, 1980), as implemented in CLM4 and Noah-MP]
    call photosynthesis(.false., ix_bbAssimFnc, ci, co2compPt, awb, cp2, vcmax, Js, psn, unUsedOutput)
  
    ! compute co2 concentration at leaf surface (Pa)
    x1 = h2o_co2__leafbl * airpres * rlb  ! Pa / (umol co2 m-2 s-1)
    cs = max(scalarCO2air - (x1 * psn), mpe)   ! Pa (avoid divide by zero)
  
    ! compute control of the compensation point on stomatal conductance
    if(ix_bbCO2point == origBWB)then
     csx = cs
    else
     csx = cs - co2compPt
    end if
  
    ! compute conductance in the absence of humidity
    g0 = cond2photo_slope*airpres*psn/csx
  
    ! use quadratic function to compute stomatal resistance
    call quadResist(.false.,ix_bbHumdFunc,rlb,fHum,gMin,g0,unUsedInput,rs,unUsedOutput)
  
    ! compute intercellular co2 partial pressues (Pa)
    x2 = h2o_co2__stomPores * airpres  ! Pa
    testFunc = max(cs - x2*psn*rs, 0._rkind)    ! Pa
  
    end function testFunc
  
   end subroutine stomResist_flex
  
   attributes(device) subroutine ballberryflex_kernel(&
                        scalarCanopySunlitPAR, scalarIntercellularCO2Sunlit, &
                        scalarCanopyShadedPAR, scalarIntercellularCO2Shaded, &
                        scalarStomResistSunlit, scalarPhotosynthesisSunlit, &
                        scalarStomResistShaded, scalarPhotosynthesisShaded, &
                        ! input: state and diagnostic variables
                        scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                        scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                        scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                        ! absorbedPAR,                         & ! intent(in): absorbed PAR (W m-2)
                        ! input: data structures
                        ! forc_data,                           & ! intent(in): model forcing data
                        airtemp, airpres, &
                        ! mpar_data,                           & ! intent(in): model parameters
                        Kc25, Ko25, Kc_qFac, Ko_qFac, &
                        kc_Ha, ko_Ha, &
                        vcmax25_canopyTop, vcmax_qFac, vcmax_Ha, vcmax_Hd, vcmax_Sv, vcmax_Kn, &
                        jmax25_scale, jmax_Ha, jmax_Hd, jmax_Sv, &
                        fractionJ, quantamYield, vpScaleFactor, cond2photo_slope, minStomatalConductance, &
                        ! diag_data,                           & ! intent(in): model diagnostic variables for a local HRU
                        scalarO2air, scalarCO2air, &
                        scalarExposedLAI, scalarGrowingSeasonIndex, scalarFoliageNitrogenFactor, scalarTranspireLim, &
                        ! flux_data,                           & ! intent(in): model fluxes for a local HRU
                        scalarLeafResistance, &
                        ! model_decisions,                     & ! intent(in): model decisions
                        ix_bbTempFunc, &
                        ix_bbHumdFunc, &
                        ix_bbElecFunc, &
                        ix_bbCO2point, &
                        ix_bbNumerics, &
                        ix_bbAssimFnc, &
                        ix_bbCanIntg8)
                        ! output: stomatal resistance and photosynthesis
                        ! ci,                                  & ! intent(out): co2 of the leaf interior (Pa)
                        ! scalarStomResist,                    & ! intent(out): stomatal resistance (s m-1)
                        ! scalarPhotosynthesis)!                & ! intent(out): photosynthesis (umol CO2 m-2 s-1)
                        ! output: error control
                        ! err,message)                           ! intent(out): error control
   use device_data_types
  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   real(rkind),intent(in) :: scalarCanopySunlitPAR
   real(rkind),intent(inout) :: scalarIntercellularCO2Sunlit
   real(rkind),intent(in) :: scalarCanopyShadedPAR
   real(rkind),intent(inout) :: scalarIntercellularCO2Shaded
   real(rkind),intent(inout) :: scalarStomResistSunlit, scalarPhotosynthesisSunlit
   real(rkind),intent(inout) :: scalarStomResistShaded, scalarPhotosynthesisShaded
  ! input: state and diagnostic variables
  real(rkind),intent(in),device             :: scalarVegetationTemp          ! vegetation temperature (K)
  real(rkind),intent(in),device             :: scalarSatVP_VegTemp           ! saturation vapor pressure at vegetation temperature (Pa)
  real(rkind),intent(in),device             :: scalarVP_CanopyAir            ! canopy air vapor pressure (Pa)
  ! real(rkind),intent(inout),device             :: absorbedPAR                   ! absorbed PAR (W m-2)
  ! input: data structures
  !  type(forc_data_device),intent(in)             :: forc_data                     ! model forcing data
  real(rkind),intent(in),device :: airtemp, airpres
  !  type(mpar_data_device),intent(in)       :: mpar_data                     ! model parameters
  real(rkind),intent(in),device :: Kc25, Ko25, Kc_qFac, Ko_qFac
  real(rkind),intent(in),device :: kc_Ha, ko_Ha
  real(rkind),intent(in),device :: vcmax25_canopyTop, vcmax_qFac, vcmax_Ha, vcmax_Hd, vcmax_Sv, vcmax_Kn
  real(rkind),intent(in),device :: jmax25_scale, jmax_Ha, jmax_Hd, jmax_Sv
  real(rkind),intent(in),device :: fractionJ, quantamYield, vpScaleFactor, cond2photo_slope, minStomatalConductance
  !  type(diag_data_device),intent(in)       :: diag_data                     ! diagnostic variables for a local HRU
  real(rkind),intent(in),device :: scalarO2air, scalarCO2air
  real(rkind),intent(in),device :: scalarExposedLAI, scalarGrowingSeasonIndex, scalarFoliageNitrogenFactor, scalarTranspireLim
  !  type(flux_data_device),intent(in)       :: flux_data                     ! model fluxes for a local HRU
  real(rkind),intent(in),device :: scalarLeafResistance
  !  type(model_options),intent(in)     :: model_decisions(:)            ! model decisions
  integer(i4b),intent(in) :: ix_bbTempFunc
  integer(i4b),intent(in) :: ix_bbHumdFunc
  integer(i4b),intent(in) :: ix_bbElecFunc
  integer(i4b),intent(in) :: ix_bbCO2point
  integer(i4b),intent(in) :: ix_bbNumerics
  integer(i4b),intent(in) :: ix_bbAssimFnc
  integer(i4b),intent(in) :: ix_bbCanIntg8
  ! output: stomatal resistance and photosynthesis
  ! real(rkind),intent(inout),device          :: ci                            ! intercellular co2 partial pressure (Pa)
  ! real(rkind),intent(out),device            :: scalarStomResist              ! stomatal resistance (s m-1)
  ! real(rkind),intent(out),device            :: scalarPhotosynthesis          ! photosynthesis (umol CO2 m-2 s-1)
  
  integer(i4b) :: iSunShade
  integer(i4b),parameter          :: ixSunlit=1                   ! named variable for sunlit leaves
  integer(i4b),parameter          :: ixShaded=2                   ! named variable for shaded leaves
  real(rkind) :: ci
  real(rkind) :: absorbedPAR
  real(rkind) :: scalarStomResist
  real(rkind) :: scalarPhotosynthesis
  
  
  do iSunShade = 1,2
   ! get appropriate input values
   select case(iSunShade)
    ! sunlit leaves
    case(ixSunlit)
     absorbedPAR = scalarCanopySunlitPAR         ! average absorbed par for sunlit leaves (w m-2)
     ci          = scalarIntercellularCO2Sunlit  ! co2 of the leaf interior for sunlit leaves (Pa)
    ! shaded leaves
    case(ixShaded)
     absorbedPAR = scalarCanopyShadedPAR         ! average absorbed par for shaded leaves (w m-2)
     ci          = scalarIntercellularCO2Shaded  ! co2 of the leaf interior for shaded leaves (Pa)
    ! check
   !  case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
   end select
  
   ! compute photosynthesis and stomatal resistance
   call stomResist_flex_device(&
                        ! input: state and diagnostic variables
                        scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                        scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                        scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                        absorbedPAR,                         & ! intent(in): absorbed PAR (W m-2)
                        ! input: data structures
                       !  forc_data,                           & ! intent(in): model forcing data
                        airtemp, airpres, &
                       !  mpar_data,                           & ! intent(in): model parameters
                        Kc25, Ko25, Kc_qFac, Ko_qFac, &
                        kc_Ha, ko_Ha, &
                        vcmax25_canopyTop, vcmax_qFac, vcmax_Ha, vcmax_Hd, vcmax_Sv, vcmax_Kn, &
                        jmax25_scale, jmax_Ha, jmax_Hd, jmax_Sv, &
                        fractionJ, quantamYield, vpScaleFactor, cond2photo_slope, minStomatalConductance, &
                       !  diag_data,                           & ! intent(in): model diagnostic variables for a local HRU
                        scalarO2air, scalarCO2air, &
                        scalarExposedLAI, scalarGrowingSeasonIndex, scalarFoliageNitrogenFactor, scalarTranspireLim, &
                       !  flux_data,                           & ! intent(in): model fluxes for a local HRU
                        scalarLeafResistance, &
                       !  model_decisions,                     & ! intent(in): model decisions
                        ix_bbTempFunc, &
                        ix_bbHumdFunc, &
                        ix_bbElecFunc, &
                        ix_bbCO2point, &
                        ix_bbNumerics, &
                        ix_bbAssimFnc, &
                        ix_bbCanIntg8, &
                        ! input-output
                        ci,                                  & ! intent(inout): co2 of the leaf interior (Pa)
                        ! output:
                        scalarStomResist,                    & ! intent(out): stomatal resistance (s m-1)
                        scalarPhotosynthesis)!,                & ! intent(out): photosynthesis (umol CO2 m-2 s-1)
                        ! output: error control
                       !  err,cmessage)                          ! intent(out): error control
   ! if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  
   ! assign output variables
   select case(iSunShade)
    case(ixSunlit)
     scalarStomResistSunlit       = scalarStomResist
     scalarPhotosynthesisSunlit   = scalarPhotosynthesis
     scalarIntercellularCO2Sunlit = ci
    case(ixShaded)
     scalarStomResistShaded       = scalarStomResist
     scalarPhotosynthesisShaded   = scalarPhotosynthesis
     scalarIntercellularCO2Shaded = ci
   !  case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
   end select
  enddo  ! looping through sunlit and shaded leaves
  end subroutine
  
  
   attributes(device) subroutine stomResist_flex_device(&
                              ! input: state and diagnostic variables
                              scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                              scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                              scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                              absorbedPAR,                         & ! intent(in): absorbed PAR (W m-2)
                              ! input: data structures
                              ! forc_data,                           & ! intent(in): model forcing data
                              airtemp, airpres, &
                              ! mpar_data,                           & ! intent(in): model parameters
                              Kc25, Ko25, Kc_qFac, Ko_qFac, &
                              kc_Ha, ko_Ha, &
                              vcmax25_canopyTop, vcmax_qFac, vcmax_Ha, vcmax_Hd, vcmax_Sv, vcmax_Kn, &
                              jmax25_scale, jmax_Ha, jmax_Hd, jmax_Sv, &
                              fractionJ, quantamYield, vpScaleFactor, cond2photo_slope, minStomatalConductance, &
                              ! diag_data,                           & ! intent(in): model diagnostic variables for a local HRU
                              scalarO2air, scalarCO2air, &
                              scalarExposedLAI, scalarGrowingSeasonIndex, scalarFoliageNitrogenFactor, scalarTranspireLim, &
                              ! flux_data,                           & ! intent(in): model fluxes for a local HRU
                              scalarLeafResistance, &
                              ! model_decisions,                     & ! intent(in): model decisions
                              ix_bbTempFunc, &
                              ix_bbHumdFunc, &
                              ix_bbElecFunc, &
                              ix_bbCO2point, &
                              ix_bbNumerics, &
                              ix_bbAssimFnc, &
                              ix_bbCanIntg8, &
                              ! output: stomatal resistance and photosynthesis
                              ci,                                  & ! intent(out): co2 of the leaf interior (Pa)
                              scalarStomResist,                    & ! intent(out): stomatal resistance (s m-1)
                              scalarPhotosynthesis)!                & ! intent(out): photosynthesis (umol CO2 m-2 s-1)
                              ! output: error control
                              ! err,message)                           ! intent(out): error control
                              use device_data_types
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! input: state and diagnostic variables
   real(rkind),intent(in),device             :: scalarVegetationTemp          ! vegetation temperature (K)
   real(rkind),intent(in),device             :: scalarSatVP_VegTemp           ! saturation vapor pressure at vegetation temperature (Pa)
   real(rkind),intent(in),device             :: scalarVP_CanopyAir            ! canopy air vapor pressure (Pa)
   real(rkind),intent(in),device             :: absorbedPAR                   ! absorbed PAR (W m-2)
   ! input: data structures
  !  type(forc_data_device),intent(in)             :: forc_data                     ! model forcing data
   real(rkind),intent(in),device :: airtemp, airpres
  !  type(mpar_data_device),intent(in)       :: mpar_data                     ! model parameters
   real(rkind),intent(in),device :: Kc25, Ko25, Kc_qFac, Ko_qFac
   real(rkind),intent(in),device :: kc_Ha, ko_Ha
   real(rkind),intent(in),device :: vcmax25_canopyTop, vcmax_qFac, vcmax_Ha, vcmax_Hd, vcmax_Sv, vcmax_Kn
   real(rkind),intent(in),device :: jmax25_scale, jmax_Ha, jmax_Hd, jmax_Sv
   real(rkind),intent(in),device :: fractionJ, quantamYield, vpScaleFactor, cond2photo_slope, minStomatalConductance
  !  type(diag_data_device),intent(in)       :: diag_data                     ! diagnostic variables for a local HRU
   real(rkind),intent(in),device :: scalarO2air, scalarCO2air
   real(rkind),intent(in),device :: scalarExposedLAI, scalarGrowingSeasonIndex, scalarFoliageNitrogenFactor, scalarTranspireLim
  !  type(flux_data_device),intent(in)       :: flux_data                     ! model fluxes for a local HRU
   real(rkind),intent(in),device :: scalarLeafResistance
  !  type(model_options),intent(in)     :: model_decisions(:)            ! model decisions
   integer(i4b),intent(in) :: ix_bbTempFunc
   integer(i4b),intent(in) :: ix_bbHumdFunc
   integer(i4b),intent(in) :: ix_bbElecFunc
   integer(i4b),intent(in) :: ix_bbCO2point
   integer(i4b),intent(in) :: ix_bbNumerics
   integer(i4b),intent(in) :: ix_bbAssimFnc
   integer(i4b),intent(in) :: ix_bbCanIntg8
   ! output: stomatal resistance and photosynthesis
   real(rkind),intent(inout),device          :: ci                            ! intercellular co2 partial pressure (Pa)
   real(rkind),intent(out),device            :: scalarStomResist              ! stomatal resistance (s m-1)
   real(rkind),intent(out),device            :: scalarPhotosynthesis          ! photosynthesis (umol CO2 m-2 s-1)
   ! output: error control
  !  integer(i4b),intent(out)           :: err                           ! error code
  !  character(*),intent(out)           :: message                       ! error message
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! general local variables
   logical(lgt),parameter             :: testDerivs=.false.            ! flag to test the derivatives
   real(rkind),device                        :: unitConv                      ! unit conversion factor (mol m-3, convert m s-1 --> mol H20 m-2 s-1)
   real(rkind),device                        :: rlb                           ! leaf boundary layer rersistance (umol-1 m2 s)
   real(rkind),device                        :: x0,x1,x2                      ! temporary variables
   real(rkind),device                        :: co2compPt                     ! co2 compensation point (Pa)
   real(rkind),device                        :: fHum                          ! humidity function, fraction [0,1]
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! fixed parameters
   integer(i4b),parameter             :: maxiter=20                    ! maximum number of iterations
   integer(i4b),parameter             :: maxiter_noahMP=3              ! maximum number of iterations for Noah-MP
   real(rkind),parameter              :: convToler=0.0001_rkind        ! convergence tolerance (Pa)
   real(rkind),parameter              :: umol_per_mol=1.e+6_rkind      ! factor to relate umol to mol
   real(rkind),parameter              :: o2scaleFactor=0.105_rkind     ! scaling factor used to compute co2 compesation point (0.21/2)
   real(rkind),parameter              :: h2o_co2__leafbl=1.37_rkind    ! factor to represent the different diffusivities of h2o and co2 in the leaf boundary layer (-)
   real(rkind),parameter              :: h2o_co2__stomPores=1.65_rkind ! factor to represent the different diffusivities of h2o and co2 in the stomatal pores (-)
   real(rkind),parameter              :: Tref=298.16_rkind             ! reference temperature (25 deg C)
   real(rkind),parameter              :: Tscale=10._rkind              ! scaling factor in q10 function (K)
   real(rkind),parameter              :: c_ps2=0.7_rkind               ! curvature factor for electron transport (-)
   real(rkind),parameter              :: fnf=0.6666666667_rkind        ! foliage nitrogen factor (-)
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! photosynthesis
   real(rkind),device                        :: Kc,Ko                         ! Michaelis-Menten constants for co2 and o2 (Pa)
   real(rkind),device                        :: vcmax25                       ! maximum Rubisco carboxylation rate at 25 deg C (umol m-2 s-1)
   real(rkind),device                        :: jmax25                        ! maximum electron transport rate at 25 deg C (umol m-2 s-1)
   real(rkind),device                        :: vcmax                         ! maximum Rubisco carboxylation rate (umol m-2 s-1)
   real(rkind),device                        :: jmax                          ! maximum electron transport rate (umol m-2 s-1)
   real(rkind),device                        :: aQuad                         ! the quadratic coefficient in the quadratic equation
   real(rkind),device                        :: bQuad                         ! the linear coefficient in the quadratic equation
   real(rkind),device                        :: cQuad                         ! the constant in the quadratic equation
   real(rkind),device                        :: bSign                         ! sign of the linear coeffcient
   real(rkind),device                        :: xTemp                         ! temporary variable in the quadratic equation
   real(rkind),device                        :: qQuad                         ! the "q" term in the quadratic equation
   real(rkind),device                        :: root1,root2                   ! roots of the quadratic function
   real(rkind),device                        :: Js                            ! scaled electron transport rate (umol co2 m-2 s-1)
   real(rkind),device                        :: I_ps2                         ! PAR absorbed by PS2 (umol photon m-2 s-1)
   real(rkind),device                        :: awb                           ! Michaelis-Menten control (Pa)
   real(rkind),device                        :: cp2                           ! additional controls in light-limited assimilation (Pa)
   real(rkind),device                        :: psn                           ! leaf gross photosynthesis rate (umol co2 m-2 s-1)
   real(rkind),device                        :: dA_dc                         ! derivative in photosynthesis w.r.t. intercellular co2 concentration
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! stomatal resistance
   real(rkind),device                        :: gMin                          ! scaled minimum conductance (umol m-2 s-1)
   real(rkind),device                        :: cs                            ! co2 partial pressure at leaf surface (Pa)
   real(rkind),device                        :: csx                           ! control of co2 partial pressure at leaf surface on stomatal conductance (Pa)
   real(rkind),device                        :: g0                            ! stomatal conductance in the absence of humidity controls (umol m-2 s-1)
   real(rkind),device                        :: ci_old                        ! intercellular co2 partial pressure (Pa)
   real(rkind),device                        :: rs                            ! stomatal resistance (umol-1 m2 s)
   real(rkind),device                        :: dg0_dc                        ! derivative in g0 w.r.t intercellular co2 concentration (umol m-2 s-1 Pa-1)
   real(rkind),device                        :: drs_dc                        ! derivative in stomatal resistance w.r.t. intercellular co2 concentration
   real(rkind),device                        :: dci_dc                        ! final derivative (-)
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! iterative solution
   real(rkind),device                        :: func1,func2                   ! functions for numerical derivative calculation
   real(rkind),device                        :: cMin,cMax                     ! solution brackets
   real(rkind),device                        :: xInc                          ! iteration increment (Pa)
   integer(i4b)                       :: iter                          ! iteration index
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! associate variables in the data structure
  !  associate(&
  !  ! input: model decisions
  !  ix_bbTempFunc                   => model_decisions(iLookDECISIONS%bbTempFunc)%iDecision,           & ! intent(in): [i4b] leaf temperature controls on photosynthesis + stomatal resistance
  !  ix_bbHumdFunc                   => model_decisions(iLookDECISIONS%bbHumdFunc)%iDecision,           & ! intent(in): [i4b] humidity controls on stomatal resistance
  !  ix_bbElecFunc                   => model_decisions(iLookDECISIONS%bbElecFunc)%iDecision,           & ! intent(in): [i4b] dependence of photosynthesis on PAR
  !  ix_bbCO2point                   => model_decisions(iLookDECISIONS%bbCO2point)%iDecision,           & ! intent(in): [i4b] use of CO2 compensation point to calculate stomatal resistance
  !  ix_bbNumerics                   => model_decisions(iLookDECISIONS%bbNumerics)%iDecision,           & ! intent(in): [i4b] iterative numerical solution method used in the Ball-Berry parameterization
  !  ix_bbAssimFnc                   => model_decisions(iLookDECISIONS%bbAssimFnc)%iDecision,           & ! intent(in): [i4b] controls on carbon assimilation (min function, or colimitation)
  !  ix_bbCanIntg8                   => model_decisions(iLookDECISIONS%bbCanIntg8)%iDecision           & ! intent(in): [i4b] scaling of photosynthesis from the leaf to the canopy
   ! input: model parameters
  !  Kc25                            => mpar_data%Kc25,                          & ! intent(in): [dp] Michaelis-Menten constant for CO2 at 25 degrees C (umol mol-1)
  !  Ko25                            => mpar_data%Ko25,                          & ! intent(in): [dp] Michaelis-Menten constant for O2 at 25 degrees C (mol mol-1)
  !  Kc_qFac                         => mpar_data%Kc_qFac,                       & ! intent(in): [dp] factor in the q10 function defining temperature controls on Kc (-)
  !  Ko_qFac                         => mpar_data%Ko_qFac,                       & ! intent(in): [dp] factor in the q10 function defining temperature controls on Ko (-)
  !  kc_Ha                           => mpar_data%kc_Ha,                         & ! intent(in): [dp] activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
  !  ko_Ha                           => mpar_data%ko_Ha,                         & ! intent(in): [dp] activation energy for the Michaelis-Menten constant for O2 (J mol-1)
  !  vcmax25_canopyTop               => mpar_data%vcmax25_canopyTop,             & ! intent(in): [dp] potential carboxylation rate at 25 degrees C at the canopy top (umol co2 m-2 s-1)
  !  vcmax_qFac                      => mpar_data%vcmax_qFac,                    & ! intent(in): [dp] factor in the q10 function defining temperature controls on vcmax (-)
  !  vcmax_Ha                        => mpar_data%vcmax_Ha,                      & ! intent(in): [dp] activation energy in the vcmax function (J mol-1)
  !  vcmax_Hd                        => mpar_data%vcmax_Hd,                      & ! intent(in): [dp] deactivation energy in the vcmax function (J mol-1)
  !  vcmax_Sv                        => mpar_data%vcmax_Sv,                      & ! intent(in): [dp] entropy term in the vcmax function (J mol-1 K-1)
  !  vcmax_Kn                        => mpar_data%vcmax_Kn,                      & ! intent(in): [dp] foliage nitrogen decay coefficient (-)
  !  jmax25_scale                    => mpar_data%jmax25_scale,                  & ! intent(in): [dp] scaling factor to relate jmax25 to vcmax25 (-)
  !  jmax_Ha                         => mpar_data%jmax_Ha,                       & ! intent(in): [dp] activation energy in the jmax function (J mol-1)
  !  jmax_Hd                         => mpar_data%jmax_Hd,                       & ! intent(in): [dp] deactivation energy in the jmax function (J mol-1)
  !  jmax_Sv                         => mpar_data%jmax_Sv,                       & ! intent(in): [dp] entropy term in the jmax function (J mol-1 K-1)
  !  fractionJ                       => mpar_data%fractionJ,                     & ! intent(in): [dp] fraction of light lost by other than the chloroplast lamellae (-)
  !  quantamYield                    => mpar_data%quantamYield,                  & ! intent(in): [dp] quantam yield (mol e mol-1 q)
  !  vpScaleFactor                   => mpar_data%vpScaleFactor,                 & ! intent(in): [dp] vapor pressure scaling factor in stomatal conductance function (Pa)
  !  cond2photo_slope                => mpar_data%cond2photo_slope,              & ! intent(in): [dp] slope of conductance-photosynthesis relationship (-)
  !  minStomatalConductance          => mpar_data%minStomatalConductance,        & ! intent(in): [dp] mimimum stomatal conductance (umol H2O m-2 s-1)
   ! input: forcing at the upper boundary
  !  airtemp                         => forc_data%airtemp,                              & ! intent(in): [dp] air temperature at some height above the surface (K)
  !  airpres                         => forc_data%airpres,                              & ! intent(in): [dp] air pressure at some height above the surface (Pa)
  !  scalarO2air                     => diag_data%scalarO2air(iGRU),                    & ! intent(in): [dp] atmospheric o2 concentration (Pa)
  !  scalarCO2air                    => diag_data%scalarCO2air(iGRU),                   & ! intent(in): [dp] atmospheric co2 concentration (Pa)
   ! input: state and diagnostic variables
  !  scalarExposedLAI                => diag_data%scalarExposedLAI(iGRU),               & ! intent(in): [dp] exposed LAI (m2 m-2)
  !  scalarGrowingSeasonIndex        => diag_data%scalarGrowingSeasonIndex(iGRU),       & ! intent(in): [dp] growing season index (0=off, 1=on)
  !  scalarFoliageNitrogenFactor     => diag_data%scalarFoliageNitrogenFactor(iGRU),    & ! intent(in): [dp] foliage nitrogen concentration (1.0 = saturated)
  !  scalarTranspireLim              => diag_data%scalarTranspireLim(iGRU)             & ! intent(in): [dp] weighted average of the transpiration limiting factor (-)
  !  scalarLeafResistance            => flux_data%scalarLeafResistance(iGRU)            & ! intent(in): [dp] mean leaf boundary layer resistance per unit leaf area (s m-1)
  !  )
   ! ------------------------------------------------------------------------------------------------------------------------------------------------------
   ! initialize error control
  !  err=0; message="stomResist_flex/"
  
     ! *****
    ! * preliminaries...
    ! ******************
  
    ! define unit conversion (m s-1 --> mol m-2 s-1)
    ! NOTE: Rgas   = J K-1 Mol-1 (J = kg m2 s-2); airtemp = K; airpres = Pa (kg m-1 s-2)
   unitConv = airpres/(Rgas*airtemp)  ! mol m-3
  
   ! check there is light available for photosynthesis
   if(absorbedPAR < tiny(absorbedPAR) .or. scalarGrowingSeasonIndex < tiny(absorbedPAR))then
     scalarStomResist     = unitConv*umol_per_mol/(scalarTranspireLim*minStomatalConductance)
     scalarPhotosynthesis = 0._rkind
     ci                   = 0._rkind
   end if
  
  
   ! scale vcmax from the leaves to the canopy
   ! exponential function of LAI is described in Leuning, Plant Cell Env 1995: "Scaling from..." [eq 9]
   select case(ix_bbCanIntg8)
     case(constantScaling); vcmax25 = vcmax25_canopyTop*fnf
     case(laiScaling);      vcmax25 = vcmax25_canopyTop*exp(-vcmax_Kn*scalarExposedLAI)
     ! case default; err=20; message=trim(message)//'unable to identify option to scale lai from the leaves to the canopy'; return
  
   end select
  
   ! compute the maximum electron transport rate at 25 deg C (umol m-2 s-1)
   jmax25 = jmax25_scale * vcmax25
  
   ! define the scaled minimum conductance
   gMin = scalarTranspireLim*minStomatalConductance
  
   ! compute the leaf conductance (umol m-2 s-1)
   rlb = scalarLeafResistance/(umol_per_mol*unitConv)  ! s m-1 --> umol-1 m2 s
  
   ! *****
   ! * compute temperature controls on stomatal conductance...
   ! *********************************************************
  
   ! identify the temperature function
   select case(ix_bbTempFunc)
  
     ! q10 function used in CLM4 and Noah-MP
     case(q10Func)
       ! compute the Michaelis-Menten constants (Pa)
       Kc = ave_slp*(Kc25/umol_per_mol)*q10_device(Kc_qFac,scalarVegetationTemp,Tref,Tscale)  ! umol mol-1 --> mol mol-1 --> Pa
       Ko = ave_slp*Ko25*q10_device(Ko_qFac,scalarVegetationTemp,Tref,Tscale)  ! mol mol-1 --> Pa
       ! compute maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
       x0 = q10_device(vcmax_qFac,scalarVegetationTemp,Tref,Tscale)  ! temperature function
       x1 = fHigh_device(vcmax_Hd,vcmax_Sv,scalarVegetationTemp) ! high temperature inhibition function
       vcmax = scalarTranspireLim*vcmax25*x0/x1
  
     ! Arrhenius function used in CLM5 and Cable
     case(Arrhenius)
       ! compute the Michaelis-Menten constants (Pa)
       Kc = airpres*(Kc25/umol_per_mol)*fT_device(kc_Ha,scalarVegetationTemp,Tref)  ! umol mol-1 --> mol mol-1 --> Pa
       Ko = airpres*Ko25*fT_device(ko_Ha,scalarVegetationTemp,Tref)  ! mol mol-1 --> Pa
       ! compute maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
       x0 = fT_device(vcmax_Ha,scalarVegetationTemp,Tref)
       x1 = fHigh_device(vcmax_Hd,vcmax_Sv,Tref) / fHigh_device(vcmax_Hd,vcmax_Sv,scalarVegetationTemp)
       vcmax = scalarTranspireLim*vcmax25*x0*x1
  
  
     ! check found an appropriate option
     ! case default; err=20; message=trim(message)//'unable to find option for leaf temperature controls on stomatal conductance'; return
  
   end select  ! temperature controls
  
   ! *****
   ! * compute electron transport controls on stomatal conductance...
   ! ****************************************************************
  
   ! compute the maximum electron transport rate (umol electron m-2 s-1)
   x0 = fT_device(jmax_Ha,scalarVegetationTemp,Tref)
   x1 = fHigh_device(jmax_Hd,jmax_Sv,Tref) / fHigh_device(jmax_Hd,jmax_Sv,scalarVegetationTemp)
   jmax = jmax25*x0*x1
  
   ! identify the electron transport function
   select case(ix_bbElecFunc)
  
     ! linear model, as used in CLM4 and Noah-MP
     case(linear)
       Js = quantamYield*joule2umolConv*absorbedPAR
  
     ! linear function of qmax, as used in Cable [Wang et al., Ag Forest Met 1998, eq D5]
     case(linearJmax)
       x0 = quantamYield*joule2umolConv*absorbedPAR
       x1 = x0*jmax / (x0 + 2.1_rkind*jmax)
       Js = x1/4._rkind ! scaled electron transport
  
     ! quadraric function of jmax, as used in CLM5 (Bonan et al., JGR 2011, Table B2)
     case(quadraticJmax)
       !  PAR absorbed by PS2 (umol photon m-2 s-1)
       I_ps2 = 0.5_rkind*(1._rkind - fractionJ) * joule2umolConv*absorbedPAR   ! Farquar (1980), eq 8: PAR absorbed by PS2 (umol photon m-2 s-1)
       ! define coefficients in the quadratic equation
       aQuad = c_ps2            ! quadratic coefficient = cuurvature factor for electron transport
       bQuad = -(I_ps2 + jmax)  ! linear coefficient
       cQuad =  I_ps2 * jmax    ! free term
       ! compute the q term (NOTE: bQuad is always positive)
       bSign = abs(bQuad)/bQuad
       xTemp = bQuad*bQuad - 4._rkind *aQuad*cQuad
       qQuad = -0.5_rkind * (bQuad + bSign*sqrt(xTemp))
       ! compute roots
       root1 = qQuad / aQuad
       root2 = cQuad / qQuad
       ! select minimum root, required to ensure J=0 when par=0
       ! NOTE: Wittig et al. select the first root, which is the max in all cases I tried
       Js = min(root1,root2) / 4._rkind  ! scaled J
  
     ! check found an appropriate option
     ! case default; err=20; message=trim(message)//'unable to find option for electron transport controls on stomatal conductance'; return
  
   end select  ! electron transport controls
  
   ! *****
   ! * define additional controls on stomatal conductance...
   ! ****************************************************************
  
   ! define the humidity function
   select case(ix_bbHumdFunc)
     case(humidLeafSurface); fHum = min( max(0.25_rkind, scalarVP_CanopyAir/scalarSatVP_VegTemp), 1._rkind)
     case(scaledHyperbolic); fHum = (scalarSatVP_VegTemp - scalarVP_CanopyAir)/vpScaleFactor
     ! case default; err=20; message=trim(message)//'unable to identify humidity control on stomatal conductance'; return
   end select
  
  
   ! compute the co2 compensation point (Pa)
   co2compPt = (Kc/Ko)*scalarO2air*o2scaleFactor
  
   ! compute the Michaelis-Menten controls (Pa)
   awb = Kc*(1._rkind + scalarO2air/Ko)
  
   ! compute the additional controls in light-limited assimilation
   cp2 = co2compPt*2._rkind
  
   ! define trial value of intercellular co2 (Pa)
   ! NOTE: only initialize if less than the co2 compensation point; otherwise, initialize with previous value
   if(ix_bbNumerics==newtonRaphson)then
     if(ci < co2compPt) ci = 0.7_rkind*scalarCO2air
   else
     ci = 0.7_rkind*scalarCO2air  ! always initialize if not NR
   end if
  
   ! initialize brackets for the solution
   cMin = 0._rkind
   cMax = scalarCO2air
  
  
   ! ***
   ! iterate
   do iter=1,maxiter
  
     ! reset ci
     ci_old = ci
  
     ! *****
     ! * compute photosynthesis and stomatal resistance...
     ! ***************************************************
  
     ! compute gross photosynthesis [follow Farquar (Planta, 1980), as implemented in CLM4 and Noah-MP]
     call photosynthesis_device(.true., ix_bbAssimFnc, ci, co2compPt, awb, cp2, vcmax, Js, psn, dA_dc)
  
     ! compute co2 concentration at leaf surface (Pa)
     x1 = h2o_co2__leafbl * airpres * rlb  ! Pa / (umol co2 m-2 s-1)
     cs = max(scalarCO2air - (x1 * psn), mpe)   ! Pa (avoid divide by zero)
  
     ! compute control of the compensation point on stomatal conductance
     if(ix_bbCO2point == origBWB)then
       csx = cs
     else
       csx = cs - co2compPt
     end if
  
     ! compute conductance in the absence of humidity
     g0     = cond2photo_slope*airpres*psn/csx
     dg0_dc = cond2photo_slope*airpres*dA_dc*(x1*psn/cs + 1._rkind)/csx
  
     ! use quadratic function to compute stomatal resistance
     call quadResist_device(.true.,ix_bbHumdFunc,rlb,fHum,gMin,g0,dg0_dc,rs,drs_dc)
  
     ! compute intercellular co2 partial pressues (Pa)
     x2 = h2o_co2__stomPores * airpres  ! Pa
     ci = max(cs - x2*psn*rs, 0._rkind)    ! Pa
  
     ! final derivative
     if(ci > tiny(ci))then
       dci_dc = -x1*dA_dc - x2*(psn*drs_dc + rs*dA_dc)
     else
       dci_dc = 0._rkind
     end if
  
     ! ! test derivatives
     ! if(testDerivs)then
     !  func1 = testFunc(ci_old,    cond2photo_slope, airpres, scalarCO2air, ix_bbHumdFunc, ix_bbCO2point, ix_bbAssimFnc)
     !  func2 = testFunc(ci_old+dx, cond2photo_slope, airpres, scalarCO2air, ix_bbHumdFunc, ix_bbCO2point, ix_bbAssimFnc)
     !  write(*,'(a,1x,20(e20.10,1x))') '(func2 - func1)/dx, dci_dc = ', &
     !                                   (func2 - func1)/dx, dci_dc
     ! end if  ! if testing the derivatives
  
     ! *****
     ! * iterative solution...
     ! ***********************
  
     ! case for Noah-MP (use fixed point iteration)
     if(ix_bbNumerics==NoahMPsolution)then
       if(iter==maxiter_NoahMP) exit ! exit after a specified number of iterations (normally 3)
       cycle  ! fixed-point iteration
     end if
  
     ! CLM4 and Noah-MP use fixed point iteration, continuing the next iteration from this point
     ! here we try and improve matters by using the derivatives
  
     ! update the brackets for the solution
     if(ci_old > ci)then
       cMax = ci_old
     else
       cMin = ci_old
     end if
  
     ! compute iteration increment (Pa)
     xInc = (ci - ci_old)/(1._rkind - dci_dc)
  
     ! update
     ci = max(ci_old + xInc, 0._rkind)
  
     ! ensure that we stay within brackets
     if(ci > cMax .or. ci < cMin)then
       ci = 0.5_rkind * (cMin + cMax)
     end if
  
     ! check for convergence
     if(abs(xInc) < convToler) exit
       if(iter==maxIter)then
       !  message=trim(message)//'did not converge in stomatal conductance iteration'
       !  err=20; 
       return
     end if
  
   end do  ! iterating
  
   ! assign output variables
   scalarStomResist     = unitConv*umol_per_mol*rs  ! umol-1 m2 s --> s/m
   scalarPhotosynthesis = psn
  
  
  
  !  end associate
  
   end subroutine stomResist_flex_device
  
   ! *******************************************************************************************************
   ! private subroutine photosynthesis: compute gross photosynthesis
   ! *******************************************************************************************************
   subroutine photosynthesis(desireDeriv, ix_bbAssimFnc, ci, co2compPt, awb, cp2, vcmax, Js, psn, dA_dc)
   implicit none
   ! dummy variables
   logical(lgt),intent(in) :: desireDeriv   ! .true. if the derivative is desired
   integer(i4b),intent(in) :: ix_bbAssimFnc ! model option for the function used for co2 assimilation (min func, or colimtation)
   real(rkind),intent(in)     :: ci            ! intercellular co2 concentration (Pa)
   real(rkind),intent(in)     :: co2compPt     ! co2 compensation point (Pa)
   real(rkind),intent(in)     :: awb           ! Michaelis-Menten control (Pa)
   real(rkind),intent(in)     :: cp2           ! additional controls in light-limited assimilation (Pa)
   real(rkind),intent(in)     :: vcmax         ! maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
   real(rkind),intent(in)     :: Js            ! scaled electron transport rate (umol co2 m-2 s-1)
   real(rkind),intent(out)    :: psn           ! leaf gross photosynthesis rate (umol co2 m-2 s-1)
   real(rkind),intent(out)    :: dA_dc         ! derivative in photosynthesis w.r.t. intercellular co2 concentration (umol co2 m-2 s-1 Pa-1)
   ! local variables
   integer(i4b),parameter          :: nFactors=3                  ! number of limiting factors for assimilation (light, Rubisco, and export)
   integer(i4b),parameter          :: ixRubi=1                    ! named variable for Rubisco-limited assimilation
   integer(i4b),parameter          :: ixLight=2                   ! named variable for light-limited assimilation
   integer(i4b),parameter          :: ixExport=3                  ! named variable for export-limited assimilation
   integer(i4b)                    :: ixLimitVec(1),ixLimit       ! index of factor limiting assimilation
   real(rkind)                        :: xFac(nFactors)              ! temporary variable used to compute assimilation rate
   real(rkind)                        :: xPSN(nFactors)              ! assimilation rate for different factors (light, Rubisco, and export)
   real(rkind)                        :: ciDiff                      ! difference between intercellular co2 and the co2 compensation point
   real(rkind)                        :: ciDer                       ! factor to account for constainted intercellular co2 in calculating derivatives
   real(rkind)                        :: x0                          ! temporary variable
   real(rkind)                        :: xsPSN                       ! intermediate smoothed photosynthesis
   real(rkind)                        :: dAc_dc,dAj_dc,dAe_dc,dAi_dc ! derivatives in assimilation w.r.t. intercellular co2 concentration
   real(rkind),parameter              :: theta_cj=0.98_rkind            ! coupling coefficient (see Sellers et al., 1996 [eq C6]; Bonan et al., 2011 [Table B1])
   real(rkind),parameter              :: theta_ie=0.95_rkind            ! coupling coefficient (see Sellers et al., 1996 [eq C6]; Bonan et al., 2011 [Table B1])
   ! ------------------------------------------------------------
   ! this method follows Farquar (Planta, 1980), as implemented in CLM4 and Noah-MP
  
   ! compute the difference between intercellular co2 concentraion and the compensation point
   ciDiff = max(0._rkind, ci - co2compPt)
  
   ! impose constraints (NOTE: derivative is zero if constraints are imposed)
   if(ci < co2compPt)then; ciDer = 0._rkind; else; ciDer = 1._rkind; end if
  
   ! compute Rubisco-limited assimilation
   xFac(ixRubi) = vcmax/(ci + awb)      ! umol co2 m-2 s-1 Pa-1
   xPSN(ixRubi) = xFac(ixRubi)*ciDiff   ! umol co2 m-2 s-1
  
   ! compute light-limited assimilation
   xFac(ixLight) = Js/(ci + cp2)        ! umol co2 m-2 s-1 Pa-1
   xPSN(ixLight) = xFac(ixLight)*ciDiff ! umol co2 m-2 s-1
  
   ! compute export limited assimilation
   xFac(ixExport) = 0.5_rkind
   xPSN(ixExport) = xFac(ixExport)*vcmax   ! umol co2 m-2 s-1
  
   ! print progress
   !write(*,'(a,1x,10(f20.10,1x))') 'xPSN, vcmax, Js = ', xPSN, vcmax, Js
  
   ! select function used for carbon assimilation
   select case(ix_bbAssimFnc)
  
    ! minimum function, as used in NoahMP (from CLM4)
    case(minFunc)
  
     ! identify limiting factor
     ixLimitVec = minloc(xPSN)
     ixLimit = ixLimitVec(1)
  
     ! define photosynthesis
     x0  = xFac(ixLimit)
     psn = xPSN(ixLimit)
  
     ! if derivatives are desired
     if(desireDeriv)then
  
      ! compute derivatives in assimilation (no colimitation)
      select case(ixLimit)
       case(ixRubi);   dA_dc = x0*ciDer - ciDiff*x0*x0/vcmax  ! Rubisco-limited assimilation
       case(ixLight);  dA_dc = x0*ciDer - ciDiff*x0*x0/Js     ! light-limited assimilation
       case(ixExport); dA_dc = 0._rkind                          ! export-limited assimilation
      end select
  
     ! derivatives are not desired
     else
      dA_dc = 0._rkind
     end if
  
    ! colimitation (Collatz et al., 1991; Sellers et al., 1996; Bonan et al., 2011)
    case(colimitation)
  
     ! compute derivatives for individual terms
     if(desireDeriv)then
      dAc_dc = xFac(ixRubi)*ciDer - ciDiff*xFac(ixRubi)*xFac(ixRubi)/vcmax
      dAj_dc = xFac(ixLight)*ciDer - ciDiff*xFac(ixLight)*xFac(ixLight)/Js
      dAe_dc = 0._rkind
     else
      dAc_dc = 0._rkind
      dAj_dc = 0._rkind
      dAe_dc = 0._rkind
     end if
  
     ! smooth Rubisco-limitation and light limitation
     if(ciDiff > tiny(ciDiff))then
      call quadSmooth(desireDeriv, xPSN(ixRubi), xPSN(ixLight), theta_cj, dAc_dc, dAj_dc, xsPSN, dAi_dc)
     else
      xsPSN  = 0._rkind
      dAi_dc = 0._rkind
     end if
  
     ! smooth intermediate-limitation and export limitation
     call quadSmooth(desireDeriv, xsPSN, xPSN(ixExport), theta_ie, dAi_dc, dAe_dc, psn, dA_dc)
  
    ! check case is identified
    case default; stop 'unknown option for carbon assimilation' ! abrupt stop: need to fix later
  
   end select  ! option for carbon assimilation
  
   end subroutine photosynthesis
  
   attributes(device) subroutine photosynthesis_device(desireDeriv, ix_bbAssimFnc, ci, co2compPt, awb, cp2, vcmax, Js, psn, dA_dc)
   implicit none
   ! dummy variables
   logical(lgt),intent(in) :: desireDeriv   ! .true. if the derivative is desired
   integer(i4b),intent(in) :: ix_bbAssimFnc ! model option for the function used for co2 assimilation (min func, or colimtation)
   real(rkind),intent(in)     :: ci            ! intercellular co2 concentration (Pa)
   real(rkind),intent(in)     :: co2compPt     ! co2 compensation point (Pa)
   real(rkind),intent(in)     :: awb           ! Michaelis-Menten control (Pa)
   real(rkind),intent(in)     :: cp2           ! additional controls in light-limited assimilation (Pa)
   real(rkind),intent(in)     :: vcmax         ! maximum Rubisco carboxylation rate (umol co2 m-2 s-1)
   real(rkind),intent(in)     :: Js            ! scaled electron transport rate (umol co2 m-2 s-1)
   real(rkind),intent(out)    :: psn           ! leaf gross photosynthesis rate (umol co2 m-2 s-1)
   real(rkind),intent(out)    :: dA_dc         ! derivative in photosynthesis w.r.t. intercellular co2 concentration (umol co2 m-2 s-1 Pa-1)
   ! local variables
   integer(i4b),parameter          :: nFactors=3                  ! number of limiting factors for assimilation (light, Rubisco, and export)
   integer(i4b),parameter          :: ixRubi=1                    ! named variable for Rubisco-limited assimilation
   integer(i4b),parameter          :: ixLight=2                   ! named variable for light-limited assimilation
   integer(i4b),parameter          :: ixExport=3                  ! named variable for export-limited assimilation
   integer(i4b)                    :: ixLimitVec(1),ixLimit       ! index of factor limiting assimilation
   real(rkind)                        :: xFac(nFactors)              ! temporary variable used to compute assimilation rate
   real(rkind)                        :: xPSN(nFactors)              ! assimilation rate for different factors (light, Rubisco, and export)
   real(rkind)                        :: ciDiff                      ! difference between intercellular co2 and the co2 compensation point
   real(rkind)                        :: ciDer                       ! factor to account for constainted intercellular co2 in calculating derivatives
   real(rkind)                        :: x0                          ! temporary variable
   real(rkind)                        :: xsPSN                       ! intermediate smoothed photosynthesis
   real(rkind)                        :: dAc_dc,dAj_dc,dAe_dc,dAi_dc ! derivatives in assimilation w.r.t. intercellular co2 concentration
   real(rkind),parameter              :: theta_cj=0.98_rkind            ! coupling coefficient (see Sellers et al., 1996 [eq C6]; Bonan et al., 2011 [Table B1])
   real(rkind),parameter              :: theta_ie=0.95_rkind            ! coupling coefficient (see Sellers et al., 1996 [eq C6]; Bonan et al., 2011 [Table B1])
   ! ------------------------------------------------------------
   ! this method follows Farquar (Planta, 1980), as implemented in CLM4 and Noah-MP
  
   ! compute the difference between intercellular co2 concentraion and the compensation point
   ciDiff = max(0._rkind, ci - co2compPt)
  
   ! impose constraints (NOTE: derivative is zero if constraints are imposed)
   if(ci < co2compPt)then; ciDer = 0._rkind; else; ciDer = 1._rkind; end if
  
   ! compute Rubisco-limited assimilation
   xFac(ixRubi) = vcmax/(ci + awb)      ! umol co2 m-2 s-1 Pa-1
   xPSN(ixRubi) = xFac(ixRubi)*ciDiff   ! umol co2 m-2 s-1
  
   ! compute light-limited assimilation
   xFac(ixLight) = Js/(ci + cp2)        ! umol co2 m-2 s-1 Pa-1
   xPSN(ixLight) = xFac(ixLight)*ciDiff ! umol co2 m-2 s-1
  
   ! compute export limited assimilation
   xFac(ixExport) = 0.5_rkind
   xPSN(ixExport) = xFac(ixExport)*vcmax   ! umol co2 m-2 s-1
  
   ! print progress
   !write(*,'(a,1x,10(f20.10,1x))') 'xPSN, vcmax, Js = ', xPSN, vcmax, Js
  
   ! select function used for carbon assimilation
   select case(ix_bbAssimFnc)
  
    ! minimum function, as used in NoahMP (from CLM4)
    case(minFunc)
  
     ! identify limiting factor
     ixLimitVec = minloc(xPSN)
     ixLimit = ixLimitVec(1)
  
     ! define photosynthesis
     x0  = xFac(ixLimit)
     psn = xPSN(ixLimit)
  
     ! if derivatives are desired
     if(desireDeriv)then
  
      ! compute derivatives in assimilation (no colimitation)
      select case(ixLimit)
       case(ixRubi);   dA_dc = x0*ciDer - ciDiff*x0*x0/vcmax  ! Rubisco-limited assimilation
       case(ixLight);  dA_dc = x0*ciDer - ciDiff*x0*x0/Js     ! light-limited assimilation
       case(ixExport); dA_dc = 0._rkind                          ! export-limited assimilation
      end select
  
     ! derivatives are not desired
     else
      dA_dc = 0._rkind
     end if
  
    ! colimitation (Collatz et al., 1991; Sellers et al., 1996; Bonan et al., 2011)
    case(colimitation)
  
     ! compute derivatives for individual terms
     if(desireDeriv)then
      dAc_dc = xFac(ixRubi)*ciDer - ciDiff*xFac(ixRubi)*xFac(ixRubi)/vcmax
      dAj_dc = xFac(ixLight)*ciDer - ciDiff*xFac(ixLight)*xFac(ixLight)/Js
      dAe_dc = 0._rkind
     else
      dAc_dc = 0._rkind
      dAj_dc = 0._rkind
      dAe_dc = 0._rkind
     end if
  
     ! smooth Rubisco-limitation and light limitation
     if(ciDiff > tiny(ciDiff))then
      call quadSmooth(desireDeriv, xPSN(ixRubi), xPSN(ixLight), theta_cj, dAc_dc, dAj_dc, xsPSN, dAi_dc)
     else
      xsPSN  = 0._rkind
      dAi_dc = 0._rkind
     end if
  
     ! smooth intermediate-limitation and export limitation
     call quadSmooth(desireDeriv, xsPSN, xPSN(ixExport), theta_ie, dAi_dc, dAe_dc, psn, dA_dc)
  
    ! check case is identified
    ! case default; stop 'unknown option for carbon assimilation' ! abrupt stop: need to fix later
  
   end select  ! option for carbon assimilation
  
   end subroutine photosynthesis_device
  
  
   ! *******************************************************************************************************
   ! private subroutine quadResist: compute stomatal resistance
   ! *******************************************************************************************************
  
   ! use quadratic function to compute stomatal resistance
  
   ! this method follows CLM4, described most fully in Oleson et al. (NCAR Tech. Note, 2010)
   ! details are also provided in Sellers et al., part 1 (J. Climate, 1996) and Bonan et al. (JGR 2011)
  
   ! stomatal conductance can be given as
   !     1/rs = m * (A/cs) * (es/ei) * Patm + b * beta   ! see Bonan et al. (2011) for inclusion of beta in the 2nd term
   ! here es is the (unknown) vapor pressure at the leaf surface
  
   ! the photosynthesis (computed above) assumes equality in co2 gradients between the atmosphere and the leaf surface,
   !  and between the leaf surface and the leaf interior, as
   !     A = (ca - cs)/(1.37*rb*Patm) = (cs - ci)/(1.65*rs*Patm)
   !  which requires that
   !     (ea - ei)/(rb + rs) = (ea - es)/rb = (es - ei)/rb
   ! where ea is the vapor pressure in the vegetation canopy, ei is the saturated vapor pressure at the leaf temperature,
   !  and then
   !     es = (ea*rs + ei*rb) / (rb + rs)
   ! more details are in Appendix C of Sellers et al. (J. Climate 1996) and Oleson et al. (NCAR Tech. Note, 2010)
  
   ! substituting the expression for es in the eqn for stomatal conductance provides the quadratic function,
   ! as described by Oleson et al. (2010)
  
   ! stomatal resistance is the larger of two roots in the quadratic equation
  
   ! -----------------------------------------------------------------------------------------------------------------
   subroutine quadResist(desireDeriv,ix_bbHumdFunc,rlb,fHum,gMin,g0,dg0_dc,rs,drs_dc)
   implicit none
   ! dummy variables
   logical(lgt),intent(in)    :: desireDeriv   ! flag to denote if the derivative is desired
   integer(i4b),intent(in)    :: ix_bbHumdFunc ! option for humidity control on stomatal resistance
   real(rkind),intent(in)     :: rlb           ! leaf boundary layer resistance (umol-1 m2 s)
   real(rkind),intent(in)     :: fHum          ! scaled humidity function (-)
   real(rkind),intent(in)     :: gMin          ! scaled minimum stomatal consuctance (umol m-2 s-1)
   real(rkind),intent(in)     :: g0            ! stomatal conductance in the absence of humidity controls (umol m-2 s-1)
   real(rkind),intent(in)     :: dg0_dc        ! derivative in g0 w.r.t intercellular co2 concentration (umol m-2 s-1 Pa-1)
   real(rkind),intent(out)    :: rs            ! stomatal resistance ((umol-1 m2 s)
   real(rkind),intent(out)    :: drs_dc        ! derivaive in rs w.r.t intercellular co2 concentration (umol-1 m2 s Pa-1)
   ! local variables
   real(rkind)                :: aQuad,bQuad,cQuad ! coefficients in the quadratic function
   real(rkind)                :: bSign,xTemp,qQuad ! q term in the quadratic
   real(rkind)                :: root1,root2       ! roots of the quadratic
   real(rkind)                :: dxT_dc,dqq_dc     ! derivatives in the q term
  
   ! define terms for the quadratic function
   select case(ix_bbHumdFunc)
  
    ! original Ball-Berry
    case(humidLeafSurface)
     aQuad = g0*fHum + gMin
     bQuad = (g0 + gMin)*rlb - 1._rkind
     cQuad = -rlb
  
    ! Leuning 1995
    case(scaledHyperbolic)
     aQuad =  g0 + gMin*(1._rkind + fHum)
     bQuad = (g0 + gMin)*rlb - fHum - 1._rkind
     cQuad = -rlb
  
   end select
  
   ! compute the q term in the quadratic
   bSign = abs(bQuad)/bQuad
   xTemp = bQuad*bQuad - 4._rkind *aQuad*cQuad
   qquad = -0.5_rkind * (bQuad + bSign*sqrt(xTemp))
  
   ! compute roots
   root1 = qQuad / aQuad
   root2 = cQuad / qQuad
   rs = max(root1,root2)
  
   ! compute derivatives
   if(desireDeriv)then
  
    ! compute derivatives in qquad w.r.t. ci
    select case(ix_bbHumdFunc)
     case(humidLeafSurface); dXt_dc = dg0_dc*(rlb*bQuad*2._rkind - fHum*cQuad*4._rkind)
     case(scaledHyperbolic); dXt_dc = dg0_dc*(rlb*bQuad*2._rkind - cQuad*4._rkind)
    end select
    dqq_dc = -0.5_rkind * (rlb*dg0_dc + bSign*dXt_dc*0.5_rkind / sqrt(xTemp) )
  
    ! compute derivatives in rs
    if(root1 > root2)then
     select case(ix_bbHumdFunc)
      case(humidLeafSurface); drs_dc = (dqq_dc - root1*fHum*dg0_dc)/aQuad
      case(scaledHyperbolic); drs_dc = (dqq_dc - root1*dg0_dc)/aQuad
     end select
    else
     drs_dc = -root2*dqq_dc/qQuad
    end if
  
   ! derivatives not desired
   else
    drs_dc = 0._rkind
   end if
  
   end subroutine quadResist
  
  
   attributes(device) subroutine quadResist_device(desireDeriv,ix_bbHumdFunc,rlb,fHum,gMin,g0,dg0_dc,rs,drs_dc)
   implicit none
   ! dummy variables
   logical(lgt),intent(in)    :: desireDeriv   ! flag to denote if the derivative is desired
   integer(i4b),intent(in)    :: ix_bbHumdFunc ! option for humidity control on stomatal resistance
   real(rkind),intent(in)     :: rlb           ! leaf boundary layer resistance (umol-1 m2 s)
   real(rkind),intent(in)     :: fHum          ! scaled humidity function (-)
   real(rkind),intent(in)     :: gMin          ! scaled minimum stomatal consuctance (umol m-2 s-1)
   real(rkind),intent(in)     :: g0            ! stomatal conductance in the absence of humidity controls (umol m-2 s-1)
   real(rkind),intent(in)     :: dg0_dc        ! derivative in g0 w.r.t intercellular co2 concentration (umol m-2 s-1 Pa-1)
   real(rkind),intent(out)    :: rs            ! stomatal resistance ((umol-1 m2 s)
   real(rkind),intent(out)    :: drs_dc        ! derivaive in rs w.r.t intercellular co2 concentration (umol-1 m2 s Pa-1)
   ! local variables
   real(rkind)                :: aQuad,bQuad,cQuad ! coefficients in the quadratic function
   real(rkind)                :: bSign,xTemp,qQuad ! q term in the quadratic
   real(rkind)                :: root1,root2       ! roots of the quadratic
   real(rkind)                :: dxT_dc,dqq_dc     ! derivatives in the q term
  
   ! define terms for the quadratic function
   select case(ix_bbHumdFunc)
  
    ! original Ball-Berry
    case(humidLeafSurface)
     aQuad = g0*fHum + gMin
     bQuad = (g0 + gMin)*rlb - 1._rkind
     cQuad = -rlb
  
    ! Leuning 1995
    case(scaledHyperbolic)
     aQuad =  g0 + gMin*(1._rkind + fHum)
     bQuad = (g0 + gMin)*rlb - fHum - 1._rkind
     cQuad = -rlb
  
   end select
  
   ! compute the q term in the quadratic
   bSign = abs(bQuad)/bQuad
   xTemp = bQuad*bQuad - 4._rkind *aQuad*cQuad
   qquad = -0.5_rkind * (bQuad + bSign*sqrt(xTemp))
  
   ! compute roots
   root1 = qQuad / aQuad
   root2 = cQuad / qQuad
   rs = max(root1,root2)
  
   ! compute derivatives
   if(desireDeriv)then
  
    ! compute derivatives in qquad w.r.t. ci
    select case(ix_bbHumdFunc)
     case(humidLeafSurface); dXt_dc = dg0_dc*(rlb*bQuad*2._rkind - fHum*cQuad*4._rkind)
     case(scaledHyperbolic); dXt_dc = dg0_dc*(rlb*bQuad*2._rkind - cQuad*4._rkind)
    end select
    dqq_dc = -0.5_rkind * (rlb*dg0_dc + bSign*dXt_dc*0.5_rkind / sqrt(xTemp) )
  
    ! compute derivatives in rs
    if(root1 > root2)then
     select case(ix_bbHumdFunc)
      case(humidLeafSurface); drs_dc = (dqq_dc - root1*fHum*dg0_dc)/aQuad
      case(scaledHyperbolic); drs_dc = (dqq_dc - root1*dg0_dc)/aQuad
     end select
    else
     drs_dc = -root2*dqq_dc/qQuad
    end if
  
   ! derivatives not desired
   else
    drs_dc = 0._rkind
   end if
  
   end subroutine quadResist_device
  
   ! *****
   ! * quadratic smoother...
   ! ***********************
  
   attributes(host,device) subroutine quadSmooth(desireDeriv, x1, x2, xsFac, dx1_dc, dx2_dc, xs, dxs_dc)
   implicit none
   ! dummy variables
   logical(lgt),intent(in) :: desireDeriv       ! flag to denote if a derivative is desired
   real(rkind),intent(in)     :: x1,x2             ! variables to be smoothed
   real(rkind),intent(in)     :: xsFac             ! smoothing factor
   real(rkind),intent(in)     :: dx1_dc,dx2_dc     ! derivatives in variables w.r.t. something important
   real(rkind),intent(out)    :: xs                ! smoothed variable
   real(rkind),intent(out)    :: dxs_dc            ! derivative w.r.t. something important
   ! local variables
   real(rkind)                :: aQuad,bQuad,cQuad ! coefficients in the quadratic function
   real(rkind)                :: bSign,xTemp,qQuad ! q term in the quadratic
   real(rkind)                :: root1,root2       ! roots of the quadratic
   real(rkind)                :: dbq_dc,dcq_dc     ! derivatives in quadratic coefficients
   real(rkind)                :: dxT_dc,dqq_dc     ! derivatives in the q term
  
   ! uses the quadratic of the form
   !  xsFac*xs^2 - (x1 + x2)*xs + x1*x2 = 0
   ! to smooth variables x1 and x2
  
   ! define the terms in the quadratic
   aQuad = xsFac
   bQuad = -(x1 + x2)
   cQuad = x1*x2
  
   ! compute the q term in the quadratic
   bSign = abs(bQuad)/bQuad
   xTemp = bQuad*bQuad - 4._rkind *aQuad*cQuad
   qquad = -0.5_rkind * (bQuad + bSign*sqrt(xTemp))
  
   ! compute roots
   root1 = qQuad / aQuad
   root2 = cQuad / qQuad
   xs    = min(root1,root2)
  
   ! compute derivatives
   if(desireDeriv)then
  
    ! compute derivatives for the terms in the quadratic
    dbq_dc = -(dx1_dc + dx2_dc)
    dcq_dc = x1*dx2_dc + x2*dx1_dc
  
    ! compute derivatives for xTemp
    dxT_dc = 2._rkind*(bQuad*dbq_dc) - 4._rkind*aQuad*dcq_dc
    dqq_dc = -0.5_rkind * (dbq_dc + bsign*dxT_dc/(2._rkind*sqrt(xTemp)))
  
    ! compute derivatives in the desired root
    if(root1 < root2)then
     dxs_dc = dqq_dc/aQuad
    else
     dxs_dc = (dcq_dc - root2*dqq_dc)/qQuad
    end if
  
   ! derivatives not required
   else
    dxs_dc = 0._rkind
   end if
  
   end subroutine quadSmooth
  
  
   ! *****
   ! * temperature functions...
   ! **************************
   ! q10 function for temperature dependence
   function q10(a,T,Tmid,Tscale)
   implicit none
   real(rkind),intent(in) :: a              ! scale factor
   real(rkind),intent(in) :: T              ! temperature (K)
   real(rkind),intent(in) :: Tmid           ! point where function is one (25 deg C)
   real(rkind),intent(in) :: Tscale         ! scaling factor (K)
   real(rkind)            :: q10            ! temperature dependence (-)
   q10 = a**((T - Tmid)/Tscale)
   end function q10
  
   attributes(device) function q10_device(a,T,Tmid,Tscale)
    implicit none
    real(rkind),intent(in),device :: a              ! scale factor
    real(rkind),intent(in),device :: T              ! temperature (K)
    real(rkind),intent(in),device :: Tmid           ! point where function is one (25 deg C)
    real(rkind),intent(in),device :: Tscale         ! scaling factor (K)
    real(rkind)            :: q10_device            ! temperature dependence (-)
    q10_device = a**((T - Tmid)/Tscale)
    end function q10_device
   
   ! Arrhenius function for temperature dependence
   function fT(delH,T,Tref)
   implicit none
   real(rkind),intent(in) :: delH     ! activation energy in temperature function (J mol-1)
   real(rkind),intent(in) :: T        ! temperature (K)
   real(rkind),intent(in) :: Tref     ! reference temperature (K)
   real(rkind)            :: fT       ! temperature dependence (-)
   fT = exp((delH/(Tref*Rgas))*(1._rkind - Tref/T))  ! NOTE: Rgas = J K-1 mol-1
   end function fT
  
   attributes(device) function fT_device(delH,T,Tref)
   implicit none
   real(rkind),intent(in),device :: delH     ! activation energy in temperature function (J mol-1)
   real(rkind),intent(in),device :: T        ! temperature (K)
   real(rkind),intent(in),device :: Tref     ! reference temperature (K)
   real(rkind)            :: fT_device       ! temperature dependence (-)
   fT_device = exp((delH/(Tref*Rgas))*(1._rkind - Tref/T))  ! NOTE: Rgas = J K-1 mol-1
   end function fT_device
  
  
   ! function for high temperature inhibition
   function fHigh(delH,delS,T)
   implicit none
   real(rkind),intent(in) :: delH     ! deactivation energy in high temp inhibition function (J mol-1)
   real(rkind),intent(in) :: delS     ! entropy term in high temp inhibition function (J K-1 mol-1)
   real(rkind),intent(in) :: T        ! temperature (K)
   real(rkind)            :: fHigh    ! high temperature inhibition (-)
   fHigh = 1._rkind + exp( (delS*T - delH)/(Rgas*T) ) ! NOTE: Rgas = J K-1 mol-1
   end function fHigh
  
   attributes(device) function fHigh_device(delH,delS,T)
   implicit none
   real(rkind),intent(in),device :: delH     ! deactivation energy in high temp inhibition function (J mol-1)
   real(rkind),intent(in),device :: delS     ! entropy term in high temp inhibition function (J K-1 mol-1)
   real(rkind),intent(in),device :: T        ! temperature (K)
   real(rkind)            :: fHigh_device    ! high temperature inhibition (-)
   fHigh_device = 1._rkind + exp( (delS*T - delH)/(Rgas*T) ) ! NOTE: Rgas = J K-1 mol-1
   end function fHigh_device
  
  
   ! *******************************************************************************************************
   ! private subroutine stomResist_NoahMP: use Noah-MP routines to compute stomatal resistance
   ! *******************************************************************************************************
   subroutine stomResist_NoahMP(&
                                ! input (model decisions)
                                ixStomResist,                        & ! intent(in): choice of function for stomatal resistance
                                ! input (local attributes)
                                vegTypeIndex,                        & ! intent(in): vegetation type index
                                iLoc, jLoc,                          & ! intent(in): spatial location indices
                                ! input (forcing)
                                airtemp,                             & ! intent(in): air temperature at some height above the surface (K)
                                airpres,                             & ! intent(in): air pressure at some height above the surface (Pa)
                                scalarO2air,                         & ! intent(in): atmospheric o2 concentration (Pa)
                                scalarCO2air,                        & ! intent(in): atmospheric co2 concentration (Pa)
                                scalarCanopySunlitPAR,               & ! intent(in): average absorbed par for sunlit leaves (w m-2)
                                scalarCanopyShadedPAR,               & ! intent(in): average absorbed par for shaded leaves (w m-2)
                                ! input (state and diagnostic variables)
                                scalarGrowingSeasonIndex,            & ! intent(in): growing season index (0=off, 1=on)
                                scalarFoliageNitrogenFactor,         & ! intent(in): foliage nitrogen concentration (1=saturated)
                                scalarTranspireLim,                  & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                                scalarLeafResistance,                & ! intent(in): leaf boundary layer resistance (s m-1)
                                scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                                scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                                scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                                ! output
                                scalarStomResistSunlit,              & ! intent(out): stomatal resistance for sunlit leaves (s m-1)
                                scalarStomResistShaded,              & ! intent(out): stomatal resistance for shaded leaves (s m-1)
                                scalarPhotosynthesisSunlit,          & ! intent(out): sunlit photosynthesis (umolco2 m-2 s-1)
                                scalarPhotosynthesisShaded,          & ! intent(out): shaded photosynthesis (umolco2 m-2 s-1)
                                err,message                          ) ! intent(out): error control
   ! -----------------------------------------------------------------------------------------------------------------------------------------
   ! Modified from Noah-MP
   ! Compute stomatal resistance and photosynthesis using either
   !  1) Ball-Berry
   !  2) Jarvis
   ! See Niu et al. JGR 2011 for more details
   USE mDecisions_module, only: BallBerry,Jarvis                ! options for the choice of function for stomatal resistance
   USE NOAHMP_ROUTINES,only:stomata                             ! compute canopy resistance based on Ball-Berry
   USE NOAHMP_ROUTINES,only:canres                              ! compute canopy resistance based Jarvis
   implicit none
   ! input (model decisions)
   integer(i4b),intent(in)       :: ixStomResist                ! choice of function for stomatal resistance
   ! input (local attributes)
   integer(i4b),intent(in)       :: vegTypeIndex                ! vegetation type index
   integer(i4b),intent(in)       :: iLoc, jLoc                  ! spatial location indices
   ! input (forcing)
   real(rkind),intent(in)           :: airtemp                     ! measured air temperature at some height above the surface (K)
   real(rkind),intent(in)           :: airpres                     ! measured air pressure at some height above the surface (Pa)
   real(rkind),intent(in)           :: scalarO2air                 ! atmospheric o2 concentration (Pa)
   real(rkind),intent(in)           :: scalarCO2air                ! atmospheric co2 concentration (Pa)
   real(rkind),intent(in),target    :: scalarCanopySunlitPAR       ! average absorbed par for sunlit leaves (w m-2)
   real(rkind),intent(in),target    :: scalarCanopyShadedPAR       ! average absorbed par for shaded leaves (w m-2)
   ! input (state and diagnostic variables)
   real(rkind),intent(in)           :: scalarGrowingSeasonIndex    ! growing season index (0=off, 1=on)
   real(rkind),intent(in)           :: scalarFoliageNitrogenFactor ! foliage nitrogen concentration (1=saturated)
   real(rkind),intent(in)           :: scalarTranspireLim          ! weighted average of the soil moiture factor controlling stomatal resistance (-)
   real(rkind),intent(in)           :: scalarLeafResistance        ! leaf boundary layer resistance (s m-1)
   real(rkind),intent(in)           :: scalarVegetationTemp        ! vegetation temperature (K)
   real(rkind),intent(in)           :: scalarSatVP_VegTemp         ! saturation vapor pressure at vegetation temperature (Pa)
   real(rkind),intent(in)           :: scalarVP_CanopyAir          ! canopy air vapor pressure (Pa)
   ! output
   real(rkind),intent(out)          :: scalarStomResistSunlit      ! stomatal resistance for sunlit leaves (s m-1)
   real(rkind),intent(out)          :: scalarStomResistShaded      ! stomatal resistance for shaded leaves (s m-1)
   real(rkind),intent(out)          :: scalarPhotosynthesisSunlit  ! sunlit photosynthesis (umolco2 m-2 s-1)
   real(rkind),intent(out)          :: scalarPhotosynthesisShaded  ! sunlit photosynthesis (umolco2 m-2 s-1)
   integer(i4b),intent(out)      :: err                         ! error code
   character(*),intent(out)      :: message                     ! error message
   ! local variables
   integer(i4b),parameter        :: ixSunlit=1                  ! named variable for sunlit leaves
   integer(i4b),parameter        :: ixShaded=2                  ! named variable for shaded leaves
   integer(i4b)                  :: iSunShade                   ! index for sunlit/shaded leaves
   real(rkind),pointer              :: PAR                         ! average absorbed PAR for sunlit/shaded leaves (w m-2)
   real(rkind)                      :: scalarStomResist            ! stomatal resistance for sunlit/shaded leaves (s m-1)
   real(rkind)                      :: scalarPhotosynthesis        ! photosynthesis for sunlit/shaded leaves (umolco2 m-2 s-1)
   ! initialize error control
   err=0; message='stomResist_NoahMP/'
  
   ! loop through sunlit and shaded leaves
   do iSunShade=1,2
  
    ! get appropriate value for PAR
    select case(iSunShade)
     case(ixSunlit); PAR => scalarCanopySunlitPAR               ! average absorbed par for sunlit leaves (w m-2)
     case(ixShaded); PAR => scalarCanopyShadedPAR               ! average absorbed par for shaded leaves (w m-2)
     case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
    end select
  
    ! identify option for stomatal resistance
    select case(ixStomResist)
  
     ! Ball-Berry
     case(BallBerry)
     call stomata(&
                  ! input
                  vegTypeIndex,                       & ! intent(in): vegetation type index
                  mpe,                                & ! intent(in): prevents overflow error if division by zero
                  PAR,                                & ! intent(in): average absorbed par (w m-2)
                  scalarFoliageNitrogenFactor,        & ! intent(in): foliage nitrogen concentration (1=saturated)
                  iLoc, jLoc,                         & ! intent(in): spatial location indices
                  scalarVegetationTemp,               & ! intent(in): vegetation temperature (K)
                  scalarSatVP_VegTemp,                & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                  scalarVP_CanopyAir,                 & ! intent(in): canopy air vapor pressure (Pa)
                  airtemp,                            & ! intent(in): air temperature at some height above the surface (K)
                  airpres,                            & ! intent(in): air pressure at some height above the surface (Pa)
                  scalarO2air,                        & ! intent(in): atmospheric o2 concentration (Pa)
                  scalarCO2air,                       & ! intent(in): atmospheric co2 concentration (Pa)
                  scalarGrowingSeasonIndex,           & ! intent(in): growing season index (0=off, 1=on)
                  scalarTranspireLim,                 & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                  scalarLeafResistance,               & ! intent(in): leaf boundary layer resistance (s m-1)
                  ! output
                  scalarStomResist,                   & ! intent(out): stomatal resistance (s m-1)
                  scalarPhotosynthesis                ) ! intent(out): photosynthesis (umolco2 m-2 s-1)
  
     ! Jarvis
     case(Jarvis)
     call canres(&
                  ! input
                  PAR,                                & ! intent(in): average absorbed par (w m-2)
                  scalarVegetationTemp,               & ! intent(in): vegetation temperature (K)
                  scalarTranspireLim,                 & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                  scalarVP_CanopyAir,                 & ! intent(in): canopy air vapor pressure (Pa)
                  airpres,                            & ! intent(in): air pressure at some height above the surface (Pa)
                  ! output
                  scalarStomResist,                   & ! intent(out): stomatal resistance (s m-1)
                  scalarPhotosynthesis,               & ! intent(out): photosynthesis (umolco2 m-2 s-1)
                  ! location indices (input)
                  iLoc, jLoc                          ) ! intent(in): spatial location indices
  
     ! check identified an option
     case default; err=20; message=trim(message)//'unable to identify case for stomatal resistance'; return
  
    end select  ! (selecting option for stomatal resistance)
  
    ! assign output variables
    select case(iSunShade)
     case(ixSunlit)
      scalarStomResistSunlit     = scalarStomResist
      scalarPhotosynthesisSunlit = scalarPhotosynthesis
     case(ixShaded)
      scalarStomResistShaded     = scalarStomResist
      scalarPhotosynthesisShaded = scalarPhotosynthesis
     case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
    end select
  
   end do  ! (looping through sunlit and shaded leaves)
  
   end subroutine stomResist_NoahMP
  
  attributes(device) subroutine stomResist_NoahMP_device(&
                                ! input (model decisions)
                                ixStomResist,                        & ! intent(in): choice of function for stomatal resistance
                                ! input (local attributes)
                                vegTypeIndex,                        & ! intent(in): vegetation type index
                                iLoc, jLoc,                          & ! intent(in): spatial location indices
                                ! input (forcing)
                                airtemp,                             & ! intent(in): air temperature at some height above the surface (K)
                                airpres,                             & ! intent(in): air pressure at some height above the surface (Pa)
                                scalarO2air,                         & ! intent(in): atmospheric o2 concentration (Pa)
                                scalarCO2air,                        & ! intent(in): atmospheric co2 concentration (Pa)
                                scalarCanopySunlitPAR,               & ! intent(in): average absorbed par for sunlit leaves (w m-2)
                                scalarCanopyShadedPAR,               & ! intent(in): average absorbed par for shaded leaves (w m-2)
                                ! input (state and diagnostic variables)
                                scalarGrowingSeasonIndex,            & ! intent(in): growing season index (0=off, 1=on)
                                scalarFoliageNitrogenFactor,         & ! intent(in): foliage nitrogen concentration (1=saturated)
                                scalarTranspireLim,                  & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                                scalarLeafResistance,                & ! intent(in): leaf boundary layer resistance (s m-1)
                                scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                                scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                                scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                                bp, mp  ,c3psn, avcmx, vcmx25, ako, ko25, akc,kc25,qe25,folnmx, &
                                rgl, rsmin, rsmax, topt, hs, &
                                ! output
                                scalarStomResistSunlit,              & ! intent(out): stomatal resistance for sunlit leaves (s m-1)
                                scalarStomResistShaded,              & ! intent(out): stomatal resistance for shaded leaves (s m-1)
                                scalarPhotosynthesisSunlit,          & ! intent(out): sunlit photosynthesis (umolco2 m-2 s-1)
                                scalarPhotosynthesisShaded)!,          & ! intent(out): shaded photosynthesis (umolco2 m-2 s-1)
                                !err,message                          ) ! intent(out): error control
   ! -----------------------------------------------------------------------------------------------------------------------------------------
   ! Modified from Noah-MP
   ! Compute stomatal resistance and photosynthesis using either
   !  1) Ball-Berry
   !  2) Jarvis
   ! See Niu et al. JGR 2011 for more details
   USE mDecisions_module, only: BallBerry,Jarvis                ! options for the choice of function for stomatal resistance
   USE NOAHMP_ROUTINES,only:stomata                             ! compute canopy resistance based on Ball-Berry
   USE NOAHMP_ROUTINES,only:canres                              ! compute canopy resistance based Jarvis
   implicit none
   ! input (model decisions)
   integer(i4b),intent(in)       :: ixStomResist                ! choice of function for stomatal resistance
   ! input (local attributes)
   integer(i4b),intent(in)       :: vegTypeIndex                ! vegetation type index
   integer(i4b),intent(in)       :: iLoc, jLoc                  ! spatial location indices
   ! input (forcing)
   real(rkind),intent(in),device           :: airtemp                     ! measured air temperature at some height above the surface (K)
   real(rkind),intent(in),device           :: airpres                     ! measured air pressure at some height above the surface (Pa)
   real(rkind),intent(in),device           :: scalarO2air                 ! atmospheric o2 concentration (Pa)
   real(rkind),intent(in),device           :: scalarCO2air                ! atmospheric co2 concentration (Pa)
   real(rkind),intent(in),target    :: scalarCanopySunlitPAR       ! average absorbed par for sunlit leaves (w m-2)
   real(rkind),intent(in),target    :: scalarCanopyShadedPAR       ! average absorbed par for shaded leaves (w m-2)
   ! input (state and diagnostic variables)
   real(rkind),intent(in),device           :: scalarGrowingSeasonIndex    ! growing season index (0=off, 1=on)
   real(rkind),intent(in),device           :: scalarFoliageNitrogenFactor ! foliage nitrogen concentration (1=saturated)
   real(rkind),intent(in),device           :: scalarTranspireLim          ! weighted average of the soil moiture factor controlling stomatal resistance (-)
   real(rkind),intent(in),device           :: scalarLeafResistance        ! leaf boundary layer resistance (s m-1)
   real(rkind),intent(in),device           :: scalarVegetationTemp        ! vegetation temperature (K)
   real(rkind),intent(in),device           :: scalarSatVP_VegTemp         ! saturation vapor pressure at vegetation temperature (Pa)
   real(rkind),intent(in),device           :: scalarVP_CanopyAir          ! canopy air vapor pressure (Pa)
   real(rkind),intent(in) :: bp, mp  ,c3psn, avcmx, vcmx25, ako, ko25, akc,kc25,qe25,folnmx
   real(rkind),intent(in) :: rgl, rsmin, rsmax, topt, hs
   ! output
   real(rkind),intent(out),device          :: scalarStomResistSunlit      ! stomatal resistance for sunlit leaves (s m-1)
   real(rkind),intent(out),device          :: scalarStomResistShaded      ! stomatal resistance for shaded leaves (s m-1)
   real(rkind),intent(out),device          :: scalarPhotosynthesisSunlit  ! sunlit photosynthesis (umolco2 m-2 s-1)
   real(rkind),intent(out),device          :: scalarPhotosynthesisShaded  ! sunlit photosynthesis (umolco2 m-2 s-1)
   !integer(i4b),intent(out)      :: err                         ! error code
   !character(*),intent(out)      :: message                     ! error message
   ! local variables
   integer(i4b),parameter        :: ixSunlit=1                  ! named variable for sunlit leaves
   integer(i4b),parameter        :: ixShaded=2                  ! named variable for shaded leaves
   integer(i4b)                  :: iSunShade                   ! index for sunlit/shaded leaves
   real(rkind),pointer              :: PAR                         ! average absorbed PAR for sunlit/shaded leaves (w m-2)
   real(rkind)                      :: scalarStomResist            ! stomatal resistance for sunlit/shaded leaves (s m-1)
   real(rkind)                      :: scalarPhotosynthesis        ! photosynthesis for sunlit/shaded leaves (umolco2 m-2 s-1)
   ! initialize error control
   !err=0; message='stomResist_NoahMP/'
  
   ! loop through sunlit and shaded leaves
   do iSunShade=1,2
  
    ! get appropriate value for PAR
    select case(iSunShade)
     case(ixSunlit); PAR => scalarCanopySunlitPAR               ! average absorbed par for sunlit leaves (w m-2)
     case(ixShaded); PAR => scalarCanopyShadedPAR               ! average absorbed par for shaded leaves (w m-2)
     !case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
    end select
  
    ! identify option for stomatal resistance
    select case(ixStomResist)
  
     ! Ball-Berry
     case(BallBerry)
     call stomata_d(&
                  ! input
                  vegTypeIndex,                       & ! intent(in): vegetation type index
                  mpe,                                & ! intent(in): prevents overflow error if division by zero
                  PAR,                                & ! intent(in): average absorbed par (w m-2)
                  scalarFoliageNitrogenFactor,        & ! intent(in): foliage nitrogen concentration (1=saturated)
                  iLoc, jLoc,                         & ! intent(in): spatial location indices
                  scalarVegetationTemp,               & ! intent(in): vegetation temperature (K)
                  scalarSatVP_VegTemp,                & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                  scalarVP_CanopyAir,                 & ! intent(in): canopy air vapor pressure (Pa)
                  airtemp,                            & ! intent(in): air temperature at some height above the surface (K)
                  airpres,                            & ! intent(in): air pressure at some height above the surface (Pa)
                  scalarO2air,                        & ! intent(in): atmospheric o2 concentration (Pa)
                  scalarCO2air,                       & ! intent(in): atmospheric co2 concentration (Pa)
                  scalarGrowingSeasonIndex,           & ! intent(in): growing season index (0=off, 1=on)
                  scalarTranspireLim,                 & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                  scalarLeafResistance,               & ! intent(in): leaf boundary layer resistance (s m-1)
                  BP, mp  ,c3psn, avcmx, vcmx25, ako, ko25, akc,kc25,qe25,folnmx   , &
                  ! output
                  scalarStomResist,                   & ! intent(out): stomatal resistance (s m-1)
                  scalarPhotosynthesis                ) ! intent(out): photosynthesis (umolco2 m-2 s-1)
  
     ! Jarvis
     case(Jarvis)
     call canres_d(&
                  ! input
                  PAR,                                & ! intent(in): average absorbed par (w m-2)
                  scalarVegetationTemp,               & ! intent(in): vegetation temperature (K)
                  scalarTranspireLim,                 & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                  scalarVP_CanopyAir,                 & ! intent(in): canopy air vapor pressure (Pa)
                  airpres,                            & ! intent(in): air pressure at some height above the surface (Pa)
                  rgl, rsmin, rsmax, topt, hs, &
                  ! output
                  scalarStomResist,                   & ! intent(out): stomatal resistance (s m-1)
                  scalarPhotosynthesis,               & ! intent(out): photosynthesis (umolco2 m-2 s-1)
                  ! location indices (input)
                  iLoc, jLoc                          ) ! intent(in): spatial location indices
  
     ! check identified an option
     !case default; err=20; message=trim(message)//'unable to identify case for stomatal resistance'; return
  
    end select  ! (selecting option for stomatal resistance)
  
    ! assign output variables
    select case(iSunShade)
     case(ixSunlit)
      scalarStomResistSunlit     = scalarStomResist
      scalarPhotosynthesisSunlit = scalarPhotosynthesis
     case(ixShaded)
      scalarStomResistShaded     = scalarStomResist
      scalarPhotosynthesisShaded = scalarPhotosynthesis
     !case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
    end select
  
   end do  ! (looping through sunlit and shaded leaves)
  
   end subroutine stomResist_NoahMP_device
  
  
   attributes(device) SUBROUTINE CALHUM_d(SFCTMP, SFCPRS, Q2SAT, DQSDT2)
          IMPLICIT NONE
          REAL(rkind), INTENT(IN)       :: SFCTMP, SFCPRS
          REAL(rkind), INTENT(OUT)      :: Q2SAT, DQSDT2
          REAL(rkind), PARAMETER        :: A2=17.67,A3=273.15,A4=29.65, ELWV=2.501E6,         &
                                    A23M4=A2*(A3-A4), E0=0.611, RV=461.0,             &
                                    EPSILON=0.622
          REAL(rkind)                   :: ES, SFCPRSX
  ! Q2SAT: saturated mixing ratio
          ES = E0 * EXP ( ELWV/RV*(1./A3 - 1./SFCTMP) )
  ! convert SFCPRS from Pa to KPa
          SFCPRSX = SFCPRS*1.E-3
          Q2SAT = EPSILON * ES / (SFCPRSX-ES)
  ! convert from  g/g to g/kg
          Q2SAT = Q2SAT * 1.E3
  ! Q2SAT is currently a 'mixing ratio'
  ! DQSDT2 is calculated assuming Q2SAT is a specific humidity
          DQSDT2=(Q2SAT/(1+Q2SAT))*A23M4/(SFCTMP-A4)**2
  ! DG Q2SAT needs to be in g/g when returned for SFLX
          Q2SAT = Q2SAT / 1.E3
    END SUBROUTINE CALHUM_d
  
    attributes(device) SUBROUTINE CANRES_D (PAR   ,SFCTMP,RCSOIL ,EAH   ,SFCPRS , & !in
    rgl, rsmin, rsmax, topt, hs, &
      RC    ,PSN   ,ILOC   ,JLOC  )           !out
  
  ! --------------------------------------------------------------------------------------------------
  ! calculate canopy resistance which depends on incoming solar radiation,
  ! air temperature, atmospheric water vapor pressure deficit at the
  ! lowest model level, and soil moisture (preferably unfrozen soil
  ! moisture rather than total)
  ! --------------------------------------------------------------------------------------------------
  ! source:  Jarvis (1976), Noilhan and Planton (1989, MWR), Jacquemin and
  ! Noilhan (1990, BLM). Chen et al (1996, JGR, Vol 101(D3), 7251-7268),
  ! eqns 12-14 and table 2 of sec. 3.1.2
  ! --------------------------------------------------------------------------------------------------
  !niu    USE module_Noahlsm_utility
  ! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
  ! --------------------------------------------------------------------------------------------------
  ! inputs
  
  INTEGER,                  INTENT(IN)  :: ILOC   !grid index
  INTEGER,                  INTENT(IN)  :: JLOC   !grid index
  REAL(rkind),                     INTENT(IN)  :: PAR    !par absorbed per unit sunlit lai (w/m2)
  REAL(rkind),                     INTENT(IN)  :: SFCTMP !canopy air temperature
  REAL(rkind),                     INTENT(IN)  :: SFCPRS !surface pressure (pa)
  REAL(rkind),                     INTENT(IN)  :: EAH    !water vapor pressure (pa)
  REAL(rkind),                     INTENT(IN)  :: RCSOIL !soil moisture stress factor
  real(rkind) :: rgl, rsmin, rsmax, topt, hs
  
  !outputs
  
  REAL(rkind),                     INTENT(OUT) :: RC     !canopy resistance per unit LAI
  REAL(rkind),                     INTENT(OUT) :: PSN    !foliage photosynthesis (umolco2/m2/s)
  
  !local
  
  REAL(rkind)                                  :: RCQ
  REAL(rkind)                                  :: RCS
  REAL(rkind)                                  :: RCT
  REAL(rkind)                                  :: FF
  REAL(rkind)                                  :: Q2     !water vapor mixing ratio (kg/kg)
  REAL(rkind)                                  :: Q2SAT  !saturation Q2
  REAL(rkind)                                  :: DQSDT2 !d(Q2SAT)/d(T)
  
  
  ! RSMIN, RSMAX, TOPT, RGL, HS are canopy stress parameters set in REDPRM
  ! ----------------------------------------------------------------------
  ! initialize canopy resistance multiplier terms.
  ! ----------------------------------------------------------------------
  RC     = 0.0
  RCS    = 0.0
  RCT    = 0.0
  RCQ    = 0.0
  
  !  compute Q2 and Q2SAT
  
  Q2 = 0.622 *  EAH  / (SFCPRS - 0.378 * EAH) !specific humidity [kg/kg]
  Q2 = Q2 / (1.0 + Q2)                        !mixing ratio [kg/kg]
  
  CALL CALHUM_d(SFCTMP, SFCPRS, Q2SAT, DQSDT2)
  
  ! contribution due to incoming solar radiation
  
  FF  = 2.0 * PAR / RGL
  RCS = (FF + RSMIN / RSMAX) / (1.0+ FF)
  RCS = MAX (RCS,0.0001)
  
  ! contribution due to air temperature
  
  RCT = 1.0- 0.0016* ( (TOPT - SFCTMP)**2.0)
  RCT = MAX (RCT,0.0001)
  
  ! contribution due to vapor pressure deficit
  
  RCQ = 1.0/ (1.0+ HS * MAX(0.,Q2SAT-Q2))
  RCQ = MAX (RCQ,0.01)
  
  ! determine canopy resistance due to all factors
  
  RC  = RSMIN / (RCS * RCT * RCQ * RCSOIL)
  PSN = -999.99       ! PSN not applied for dynamic carbon
  
  END SUBROUTINE CANRES_d
  
  attributes(device) SUBROUTINE STOMATA_d (VEGTYP  ,MPE     ,APAR    ,FOLN    ,ILOC    , JLOC, & !in
                        TV      ,EI      ,EA      ,SFCTMP  ,SFCPRS  , & !in
                        O2      ,CO2     ,IGS     ,BTRAN   ,RB , &
                        bp, mp  ,c3psn, avcmx, vcmx25, ako, ko25, akc,kc25,qe25,folnmx   , & !in
                        RS      ,PSN     )                              !out
  ! --------------------------------------------------------------------------------------------------
    USE NOAHMP_VEG_PARAMETERS
  ! --------------------------------------------------------------------------------------------------
    IMPLICIT NONE
  ! --------------------------------------------------------------------------------------------------
  ! input
        INTEGER,INTENT(IN)  :: ILOC   !grid index
        INTEGER,INTENT(IN)  :: JLOC   !grid index
        INTEGER,INTENT(IN)  :: VEGTYP !vegetation physiology type
  
        REAL(rkind), INTENT(IN)    :: IGS    !growing season index (0=off, 1=on)
        REAL(rkind), INTENT(IN)    :: MPE    !prevents division by zero errors
  
        REAL(rkind), INTENT(IN)    :: TV     !foliage temperature (k)
        REAL(rkind), INTENT(IN)    :: EI     !vapor pressure inside leaf (sat vapor press at tv) (pa)
        REAL(rkind), INTENT(IN)    :: EA     !vapor pressure of canopy air (pa)
        REAL(rkind), INTENT(IN)    :: APAR   !par absorbed per unit lai (w/m2)
        REAL(rkind), INTENT(IN)    :: O2     !atmospheric o2 concentration (pa)
        REAL(rkind), INTENT(IN)    :: CO2    !atmospheric co2 concentration (pa)
        REAL(rkind), INTENT(IN)    :: SFCPRS !air pressure at reference height (pa)
        REAL(rkind), INTENT(IN)    :: SFCTMP !air temperature at reference height (k)
        REAL(rkind), INTENT(IN)    :: BTRAN  !soil water transpiration factor (0 to 1)
        REAL(rkind), INTENT(IN)    :: FOLN   !foliage nitrogen concentration (%)
        REAL(rkind), INTENT(IN)    :: RB     !boundary layer resistance (s/m)
        real(rkind),intent(in) :: bp
        real(rkind),intent(in) :: mp
        real(rkind),intent(in) :: c3psn
        real(rkind),intent(in) :: AVCMX
        real(rkind),intent(in) :: VCMX25
        real(rkind),intent(in) :: ako
        real(rkind),intent(in) :: Ko25
        real(rkind),intent(in) :: akc
        real(rkind),intent(in) :: kc25,qe25
        real(rkind),intent(in) :: folnmx
  
  ! output
        REAL(rkind), INTENT(OUT)   :: RS     !leaf stomatal resistance (s/m)
        REAL(rkind), INTENT(OUT)   :: PSN    !foliage photosynthesis (umol co2 /m2/ s) [always +]
  
  ! in&out
        REAL(rkind)                :: RLB    !boundary layer resistance (s m2 / umol)
  ! ---------------------------------------------------------------------------------------------
  
  ! ------------------------ local variables ----------------------------------------------------
        INTEGER :: ITER     !iteration index
        INTEGER :: NITER    !number of iterations
  
        DATA NITER /3/
        SAVE NITER
  
        REAL(rkind) :: AB          !used in statement functions
        REAL(rkind) :: BC          !used in statement functions
        REAL(rkind) :: F1          !generic temperature response (statement function)
        REAL(rkind) :: F2          !generic temperature inhibition (statement function)
        REAL(rkind) :: TC          !foliage temperature (degree Celsius)
        REAL(rkind) :: CS          !co2 concentration at leaf surface (pa)
        REAL(rkind) :: KC          !co2 Michaelis-Menten constant (pa)
        REAL(rkind) :: KO          !o2 Michaelis-Menten constant (pa)
        REAL(rkind) :: A,B,C,Q     !intermediate calculations for RS
        REAL(rkind) :: R1,R2       !roots for RS
        REAL(rkind) :: FNF         !foliage nitrogen adjustment factor (0 to 1)
        REAL(rkind) :: PPF         !absorb photosynthetic photon flux (umol photons/m2/s)
        REAL(rkind) :: WC          !Rubisco limited photosynthesis (umol co2/m2/s)
        REAL(rkind) :: WJ          !light limited photosynthesis (umol co2/m2/s)
        REAL(rkind) :: WE          !export limited photosynthesis (umol co2/m2/s)
        REAL(rkind) :: CP          !co2 compensation point (pa)
        REAL(rkind) :: CI          !internal co2 (pa)
        REAL(rkind) :: AWC         !intermediate calculation for wc
        REAL(rkind) :: VCMX        !maximum rate of carbonylation (umol co2/m2/s)
        REAL(rkind) :: J           !electron transport (umol co2/m2/s)
        REAL(rkind) :: CEA         !constrain ea or else model blows up
        REAL(rkind) :: CF          !s m2/umol -> s/m
        real(rkind),parameter :: tfrz = 273.16
  
        F1(AB,BC) = AB**((BC-25.)/10.)
        F2(AB) = 1. + EXP((-2.2E05+710.*(AB+273.16))/(8.314*(AB+273.16)))
        REAL(rkind) :: T
  ! ---------------------------------------------------------------------------------------------
  
  ! MPC change
  ! Original code uses gs = g0 + g1*hs*An/cs (g0 is not limited by soil moisture stress)
  ! Following Bonan et al. (JGR 2011) Improving canopy processes in the CLM... we replace the original code
  !  with gs = g0*bTran + g1*hs*An/cs, where bTran is the soil moisture stress function
  
  ! initialize RS=RSMAX and PSN=0 because will only do calculations
  ! for APAR > 0, in which case RS <= RSMAX and PSN >= 0
  
           CF = SFCPRS/(8.314*SFCTMP)*1.e06
           RS = 1./(BP*BTRAN) * CF    ! MPC change: include BTRAN multiplier
           PSN = 0.
           !write(*,'(a,1x,20(f16.6,1x))') 'TV-TFRZ, RS, CF = ', TV-TFRZ, RS, CF
  
           IF (APAR .LE. 0.) RETURN
  
           FNF = MIN( FOLN/MAX(MPE,FOLNMX), 1.0 )
           TC  = TV-TFRZ
           PPF = 4.6*APAR
           J   = PPF*QE25
           !write(*,'(a,1x,10(f20.10,1x))') 'J, APAR, QE25(VEGTYP) = ', J, APAR, QE25(VEGTYP)
           KC  = KC25 * F1(AKC,TC)
           KO  = KO25 * F1(AKO,TC)
           AWC = KC * (1.+O2/KO)
           CP  = 0.5*KC/KO*O2*0.21
           VCMX = VCMX25 / F2(TC) * FNF * BTRAN * F1(AVCMX,TC)
           !write(*,'(a,1x,20(f14.8,1x))') 'F1, F2, VCMX25(VEGTYP), AVCMX(VEGTYP), TC = ', F1(AVCMX(VEGTYP),TC), F2(TC), VCMX25(VEGTYP), AVCMX(VEGTYP), TC
  
  ! first guess ci
  
           CI = 0.7*CO2*C3PSN + 0.4*CO2*(1.-C3PSN)
  
           !write(*,'(a,1x,10(f20.10,1x))') 'KC25(VEGTYP), AKC(VEGTYP), KO25(VEGTYP), AKO(VEGTYP) = ', KC25(VEGTYP), AKC(VEGTYP), KO25(VEGTYP), AKO(VEGTYP)
           !write(*,'(a,1x,10(f20.10,1x))') 'CO2, CI, CP, KC, KO = ', CO2, CI, CP, KC, KO
  
  ! rb: s/m -> s m**2 / umol
  
           RLB = RB/CF
  
  ! constrain ea
  
           CEA = MAX(0.25*EI*C3PSN+0.40*EI*(1.-C3PSN), MIN(EA,EI) )
  
           !print*, '**'
           !write(*,'(a,1x,20(f14.8,1x))') 'BTRAN, VCMX, MP(VEGTYP), EA, EI, CEA/EI, O2, CO2, KC, KO, J, TC, RLB, SFCPRS = ', &
           !                                BTRAN, VCMX, MP(VEGTYP), EA, EI, CEA/EI, O2, CO2, KC, KO, J, TC, RLB, SFCPRS
  
  ! ci iteration
  !jref: C3PSN is equal to 1 for all veg types.
         DO ITER = 1, NITER
  
              WJ = MAX(CI-CP,0.)*J/(CI+2.*CP)*C3PSN  + J*(1.-C3PSN)
              WC = MAX(CI-CP,0.)*VCMX/(CI+AWC)*C3PSN + VCMX*(1.-C3PSN)
              WE = 0.5*VCMX*C3PSN + 4000.*VCMX*CI/SFCPRS*(1.-C3PSN)
              PSN = MIN(WJ,WC,WE) * IGS
  
              CS = MAX( CO2-1.37*RLB*SFCPRS*PSN, MPE )
  
              A = MP*PSN*SFCPRS*CEA / (CS*EI) + BTRAN*BP      ! MPC change: include BTRAN multiplier for 2nd term
              B = ( MP*PSN*SFCPRS/CS + BTRAN*BP ) * RLB - 1.  ! MPC change: include BTRAN multiplier for 2nd term in brackets
              C = -RLB
              IF (B .GE. 0.) THEN
                 Q = -0.5*( B + SQRT(B*B-4.*A*C) )
              ELSE
                 Q = -0.5*( B - SQRT(B*B-4.*A*C) )
              END IF
              R1 = Q/A
              R2 = C/Q
              RS = MAX(R1,R2)
              CI = MAX( CS-PSN*SFCPRS*1.65*RS, 0. )
  
              !write(*,'(a,1x10(f14.10,1x))') 'WJ, WC, WE, PSN, CI, CS, RS, TV, VCMX, J = ', &
              !                                WJ, WC, WE, PSN, CI, CS, RS, TV, VCMX, J
  
         END DO
         !pause
  
  ! rs, rb:  s m**2 / umol -> s/m
  
           RS = RS*CF
           !write(*,'(a,1x,10(f20.10,1x))') 'RS, G, CEA, EA, EI = ', RS, 1./RS, CEA, EA, EI
  
    END SUBROUTINE STOMATA_d
  
   ! -- end private subroutines
  
    attributes(device) subroutine stomResist_kernel(scalarTranspireLim, airpres, scalarStomResistSunlit, scalarStomResistShaded, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded,&
                        minStomatalResistance, &
                        scalarCanopySunlitPAR, scalarIntercellularCO2Sunlit, &
                        scalarCanopyShadedPAR, scalarIntercellularCO2Shaded, &
                        scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                        scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                        scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                        airtemp, &
                        Kc25, Ko25, Kc_qFac, Ko_qFac, &
                        kc_Ha, ko_Ha, &
                        vcmax25_canopyTop, vcmax_qFac, vcmax_Ha, vcmax_Hd, vcmax_Sv, vcmax_Kn, &
                        jmax25_scale, jmax_Ha, jmax_Hd, jmax_Sv, &
                        fractionJ, quantamYield, vpScaleFactor, cond2photo_slope, minStomatalConductance, &
                        scalarO2air, scalarCO2air, &
                        scalarExposedLAI, scalarGrowingSeasonIndex, scalarFoliageNitrogenFactor, &
                        scalarLeafResistance, &
                        ix_bbTempFunc, &
                        ix_bbHumdFunc, &
                        ix_bbElecFunc, &
                        ix_bbCO2point, &
                        ix_bbNumerics, &
                        ix_bbAssimFnc, &
                        ix_bbCanIntg8, &
                        ix_StomResist, &
                        vegTypeIndex,                        & ! intent(in): vegetation type index
                        iLoc, jLoc,                          & ! intent(in): spatial location indices
                        bp, mp  ,c3psn, avcmx, vcmx25, ako, ko_25, akc,kc_25,qe25,folnmx, &
                        rgl, rsmin, rsmax, topt, hs &
                        )
    real(rkind),intent(in) :: scalarTranspireLim, airpres, minStomatalResistance
    real(rkind),intent(inout) :: scalarStomResistSunlit, scalarStomResistShaded, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded
    real(rkind),intent(in) :: scalarCanopySunlitPAR
    real(rkind),intent(inout) :: scalarIntercellularCO2Sunlit
    real(rkind),intent(in) :: scalarCanopyShadedPAR
    real(rkind),intent(inout) :: scalarIntercellularCO2Shaded
   real(rkind),intent(in)             :: scalarVegetationTemp          ! vegetation temperature (K)
   real(rkind),intent(in)             :: scalarSatVP_VegTemp           ! saturation vapor pressure at vegetation temperature (Pa)
   real(rkind),intent(in)             :: scalarVP_CanopyAir            ! canopy air vapor pressure (Pa)
   real(rkind),intent(in) :: airtemp
   real(rkind),intent(in) :: Kc25, Ko25, Kc_qFac, Ko_qFac
   real(rkind),intent(in) :: kc_Ha, ko_Ha
   real(rkind),intent(in) :: vcmax25_canopyTop, vcmax_qFac, vcmax_Ha, vcmax_Hd, vcmax_Sv, vcmax_Kn
   real(rkind),intent(in) :: jmax25_scale, jmax_Ha, jmax_Hd, jmax_Sv
   real(rkind),intent(in) :: fractionJ, quantamYield, vpScaleFactor, cond2photo_slope, minStomatalConductance
   real(rkind),intent(in) :: scalarO2air, scalarCO2air
   real(rkind),intent(in) :: scalarExposedLAI, scalarGrowingSeasonIndex, scalarFoliageNitrogenFactor
   real(rkind),intent(in) :: scalarLeafResistance
   integer(i4b),intent(in) :: ix_bbTempFunc
   integer(i4b),intent(in) :: ix_bbHumdFunc
   integer(i4b),intent(in) :: ix_bbElecFunc
   integer(i4b),intent(in) :: ix_bbCO2point
   integer(i4b),intent(in) :: ix_bbNumerics
   integer(i4b),intent(in) :: ix_bbAssimFnc
   integer(i4b),intent(in) :: ix_bbCanIntg8
   integer(i4b),intent(in) :: ix_stomResist
   integer(i4b),intent(in)       :: vegTypeIndex                ! vegetation type index
   integer(i4b),intent(in)       :: iLoc, jLoc                  ! spatial location indices
   real(rkind),intent(in) :: bp, mp  ,c3psn, avcmx, vcmx25, ako, ko_25, akc,kc_25,qe25,folnmx
   real(rkind),intent(in) :: rgl, rsmin, rsmax, topt, hs
  
   
   ! identify option for stomatal resistance
   select case(ix_stomResist)
  
    ! *******************************************************************************************************************************************
  
    ! simple resistance formulation
    case(simpleResistance)
      call simpleResistance_kernel(scalarTranspireLim, airpres, scalarStomResistSunlit, scalarStomResistShaded, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded,&
      minStomatalResistance)
  
    ! *******************************************************************************************************************************************
  
    ! flexible Ball-Berry
    case(BallBerryFlex)
  
     ! loop through sunlit and shaded leaves
      call ballberryflex_kernel(&
                           scalarCanopySunlitPAR, scalarIntercellularCO2Sunlit,&
                           scalarCanopyShadedPAR, scalarIntercellularCO2Shaded, &
                           scalarStomResistSunlit, scalarPhotosynthesisSunlit, &
                           scalarStomResistShaded, scalarPhotosynthesisShaded, &
                           ! input: state and diagnostic variables
                           scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                           scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                           scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                          !  absorbedPAR(iGRU,iSunShade),                         & ! intent(in): absorbed PAR (W m-2)
                           ! input: data structures
                          !  forc_data,                           & ! intent(in): model forcing data
                           airtemp, airpres, &
                          !  mpar_data,                           & ! intent(in): model parameters
                           Kc25, Ko25, Kc_qFac, Ko_qFac, &
                           kc_Ha, ko_Ha, &
                           vcmax25_canopyTop, vcmax_qFac, vcmax_Ha, vcmax_Hd, vcmax_Sv, vcmax_Kn, &
                           jmax25_scale, jmax_Ha, jmax_Hd, jmax_Sv, &
                           fractionJ, quantamYield, vpScaleFactor, cond2photo_slope, minStomatalConductance, &
                          !  diag_data,                           & ! intent(in): model diagnostic variables for a local HRU
                           scalarO2air, scalarCO2air, &
                           scalarExposedLAI, scalarGrowingSeasonIndex, scalarFoliageNitrogenFactor, scalarTranspireLim, &
                          !  flux_data,                           & ! intent(in): model fluxes for a local HRU
                           scalarLeafResistance, &
                          !  model_decisions,                     & ! intent(in): model decisions
                           ix_bbTempFunc, &
                           ix_bbHumdFunc, &
                           ix_bbElecFunc, &
                           ix_bbCO2point, &
                           ix_bbNumerics, &
                           ix_bbAssimFnc, &
                           ix_bbCanIntg8 &
                           ! input-output
                          !  ci(iGRU,iSunShade),                                  & ! intent(inout): co2 of the leaf interior (Pa)
                           ! output:
                          !  scalarStomResist(iGRU,iSunShade),                    & ! intent(out): stomatal resistance (s m-1)
                          !  scalarPhotosynthesis(iGRU,iSunShade)
                           )!,                & ! intent(out): photosynthesis (umol CO2 m-2 s-1)
                           ! output: error control
                          !  err,cmessage)                          ! intent(out): error control
      ! if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  
      ! print progress
      !write(*,'(a,1x,20(f12.5,1x))') 'leafTemp, par, psn, rs = ', scalarVegetationTemp, absorbedPAR, scalarPhotosynthesis, scalarStomResist
  
  
    ! *******************************************************************************************************************************************
    ! compute stomatal resistance (wrapper around the Noah-MP routines)
    ! NOTE: canopy air vapor pressure is from the previous time step
    case(BallBerry,Jarvis)
     call stomResist_NoahMP_device(&
                            ! input (model decisions)
                            ix_stomResist,                     & ! intent(in): choice of function for stomatal resistance
                            ! input (local attributes)
                            vegTypeIndex,                      & ! intent(in): vegetation type index
                            iLoc, jLoc,                        & ! intent(in): spatial location indices
                            ! input (forcing)
                            airtemp,                           & ! intent(in): air temperature at some height above the surface (K)
                            airpres,                           & ! intent(in): air pressure at some height above the surface (Pa)
                            scalarO2air,                       & ! intent(in): atmospheric o2 concentration (Pa)
                            scalarCO2air,                      & ! intent(in): atmospheric co2 concentration (Pa)
                            scalarCanopySunlitPAR,             & ! intent(in): average absorbed par for sunlit leaves (w m-2)
                            scalarCanopyShadedPAR,             & ! intent(in): average absorbed par for shaded leaves (w m-2)
                            ! input (state and diagnostic variables)
                            scalarGrowingSeasonIndex,          & ! intent(in): growing season index (0=off, 1=on)
                            scalarFoliageNitrogenFactor,       & ! intent(in): foliage nitrogen concentration (1=saturated)
                            scalarTranspireLim,                & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                            scalarLeafResistance,              & ! intent(in): leaf boundary layer resistance (s m-1)
                            scalarVegetationTemp,              & ! intent(in): temperature of the vegetation canopy (K)
                            scalarSatVP_VegTemp,               & ! intent(in): saturation vapor pressure at the temperature of the veg canopy (Pa)
                            scalarVP_CanopyAir,                & ! intent(in): canopy air vapor pressure (Pa)
                            bp, mp  ,c3psn, avcmx, vcmx25, ako, ko_25, akc,kc_25,qe25,folnmx,&
                            rgl, rsmin, rsmax, topt, hs, &
                            ! output
                            scalarStomResistSunlit,            & ! intent(out): stomatal resistance for sunlit leaves (s m-1)
                            scalarStomResistShaded,            & ! intent(out): stomatal resistance for shaded leaves (s m-1)
                            scalarPhotosynthesisSunlit,        & ! intent(out): sunlit photosynthesis (umolco2 m-2 s-1)
                            scalarPhotosynthesisShaded)!,        & ! intent(out): shaded photosynthesis (umolco2 m-2 s-1)
                            !err,cmessage                       ) ! intent(out): error control
    !  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  
    ! *******************************************************************************************************************************************
  
    ! error check
    ! case default; err=20; message=trim(message)//'unable to identify option for stomatal resistance'; return
  
    ! *******************************************************************************************************************************************
  
   end select  ! (identifying option for stomatal resistance)
  end subroutine
  
  attributes(device) subroutine stomResist_kernel2(scalarTranspireLim, airpres, scalarStomResistSunlit, scalarStomResistShaded, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded,&
                        minStomatalResistance, &
                        scalarCanopySunlitPAR, scalarIntercellularCO2Sunlit, &
                        scalarCanopyShadedPAR, scalarIntercellularCO2Shaded, &
                        scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                        scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                        scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                        airtemp, &
                        Kc25, Ko25, Kc_qFac, Ko_qFac, &
                        kc_Ha, ko_Ha, &
                        vcmax25_canopyTop, vcmax_qFac, vcmax_Ha, vcmax_Hd, vcmax_Sv, vcmax_Kn, &
                        jmax25_scale, jmax_Ha, jmax_Hd, jmax_Sv, &
                        fractionJ, quantamYield, vpScaleFactor, cond2photo_slope, minStomatalConductance, &
                        scalarO2air, scalarCO2air, &
                        scalarExposedLAI, scalarGrowingSeasonIndex, scalarFoliageNitrogenFactor, &
                        scalarLeafResistance, &
                        ix_bbTempFunc, &
                        ix_bbHumdFunc, &
                        ix_bbElecFunc, &
                        ix_bbCO2point, &
                        ix_bbNumerics, &
                        ix_bbAssimFnc, &
                        ix_bbCanIntg8, &
                        ix_StomResist, &
                        vegTypeIndex,                        & ! intent(in): vegetation type index
                        bp, mp  ,c3psn, avcmx, vcmx25, ako, ko_25, akc,kc_25,qe25,folnmx, &
                        rgl, rsmin, rsmax, topt, hs &
                        )
    real(rkind),intent(in) :: scalarTranspireLim, airpres, minStomatalResistance
    real(rkind),intent(inout) :: scalarStomResistSunlit, scalarStomResistShaded, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded
    real(rkind),intent(in) :: scalarCanopySunlitPAR
    real(rkind),intent(inout) :: scalarIntercellularCO2Sunlit
    real(rkind),intent(in) :: scalarCanopyShadedPAR
    real(rkind),intent(inout) :: scalarIntercellularCO2Shaded
   real(rkind),intent(in)             :: scalarVegetationTemp          ! vegetation temperature (K)
   real(rkind),intent(in)             :: scalarSatVP_VegTemp           ! saturation vapor pressure at vegetation temperature (Pa)
   real(rkind),intent(in)             :: scalarVP_CanopyAir            ! canopy air vapor pressure (Pa)
   real(rkind),intent(in) :: airtemp
   real(rkind),intent(in) :: Kc25, Ko25, Kc_qFac, Ko_qFac
   real(rkind),intent(in) :: kc_Ha, ko_Ha
   real(rkind),intent(in) :: vcmax25_canopyTop, vcmax_qFac, vcmax_Ha, vcmax_Hd, vcmax_Sv, vcmax_Kn
   real(rkind),intent(in) :: jmax25_scale, jmax_Ha, jmax_Hd, jmax_Sv
   real(rkind),intent(in) :: fractionJ, quantamYield, vpScaleFactor, cond2photo_slope, minStomatalConductance
   real(rkind),intent(in) :: scalarO2air, scalarCO2air
   real(rkind),intent(in) :: scalarExposedLAI, scalarGrowingSeasonIndex, scalarFoliageNitrogenFactor
   real(rkind),intent(in) :: scalarLeafResistance
   integer(i4b),intent(in) :: ix_bbTempFunc
   integer(i4b),intent(in) :: ix_bbHumdFunc
   integer(i4b),intent(in) :: ix_bbElecFunc
   integer(i4b),intent(in) :: ix_bbCO2point
   integer(i4b),intent(in) :: ix_bbNumerics
   integer(i4b),intent(in) :: ix_bbAssimFnc
   integer(i4b),intent(in) :: ix_bbCanIntg8
   integer(i4b),intent(in) :: ix_stomResist
   integer(i4b),intent(in)       :: vegTypeIndex                ! vegetation type index
   real(rkind),intent(in) :: bp, mp  ,c3psn, avcmx, vcmx25, ako, ko_25, akc,kc_25,qe25,folnmx
   real(rkind),intent(in) :: rgl, rsmin, rsmax, topt, hs
  
   
   ! identify option for stomatal resistance
   select case(ix_stomResist)
  
    ! *******************************************************************************************************************************************
  
    ! simple resistance formulation
    case(simpleResistance)
      call simpleResistance_kernel(scalarTranspireLim, airpres, scalarStomResistSunlit, scalarStomResistShaded, scalarPhotosynthesisSunlit, scalarPhotosynthesisShaded,&
      minStomatalResistance)
      ! print*, 3346, scalarStomResistSunlit
  
    ! *******************************************************************************************************************************************
  
    ! flexible Ball-Berry
    case(BallBerryFlex)
  
     ! loop through sunlit and shaded leaves
      call ballberryflex_kernel(&
                           scalarCanopySunlitPAR, scalarIntercellularCO2Sunlit,&
                           scalarCanopyShadedPAR, scalarIntercellularCO2Shaded, &
                           scalarStomResistSunlit, scalarPhotosynthesisSunlit, &
                           scalarStomResistShaded, scalarPhotosynthesisShaded, &
                           ! input: state and diagnostic variables
                           scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                           scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                           scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                          !  absorbedPAR(iGRU,iSunShade),                         & ! intent(in): absorbed PAR (W m-2)
                           ! input: data structures
                          !  forc_data,                           & ! intent(in): model forcing data
                           airtemp, airpres, &
                          !  mpar_data,                           & ! intent(in): model parameters
                           Kc25, Ko25, Kc_qFac, Ko_qFac, &
                           kc_Ha, ko_Ha, &
                           vcmax25_canopyTop, vcmax_qFac, vcmax_Ha, vcmax_Hd, vcmax_Sv, vcmax_Kn, &
                           jmax25_scale, jmax_Ha, jmax_Hd, jmax_Sv, &
                           fractionJ, quantamYield, vpScaleFactor, cond2photo_slope, minStomatalConductance, &
                          !  diag_data,                           & ! intent(in): model diagnostic variables for a local HRU
                           scalarO2air, scalarCO2air, &
                           scalarExposedLAI, scalarGrowingSeasonIndex, scalarFoliageNitrogenFactor, scalarTranspireLim, &
                          !  flux_data,                           & ! intent(in): model fluxes for a local HRU
                           scalarLeafResistance, &
                          !  model_decisions,                     & ! intent(in): model decisions
                           ix_bbTempFunc, &
                           ix_bbHumdFunc, &
                           ix_bbElecFunc, &
                           ix_bbCO2point, &
                           ix_bbNumerics, &
                           ix_bbAssimFnc, &
                           ix_bbCanIntg8 &
                           ! input-output
                          !  ci(iGRU,iSunShade),                                  & ! intent(inout): co2 of the leaf interior (Pa)
                           ! output:
                          !  scalarStomResist(iGRU,iSunShade),                    & ! intent(out): stomatal resistance (s m-1)
                          !  scalarPhotosynthesis(iGRU,iSunShade)
                           )!,                & ! intent(out): photosynthesis (umol CO2 m-2 s-1)
                           ! output: error control
                          !  err,cmessage)                          ! intent(out): error control
      ! if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  
                          !  print*, 3396, scalarStomResistSunlit

      ! print progress
      !write(*,'(a,1x,20(f12.5,1x))') 'leafTemp, par, psn, rs = ', scalarVegetationTemp, absorbedPAR, scalarPhotosynthesis, scalarStomResist
  
  
    ! *******************************************************************************************************************************************
    ! compute stomatal resistance (wrapper around the Noah-MP routines)
    ! NOTE: canopy air vapor pressure is from the previous time step
    case(BallBerry,Jarvis)
     call stomResist_NoahMP_device(&
                            ! input (model decisions)
                            ix_stomResist,                     & ! intent(in): choice of function for stomatal resistance
                            ! input (local attributes)
                            vegTypeIndex,                      & ! intent(in): vegetation type index
                            iLoc, jLoc,                        & ! intent(in): spatial location indices
                            ! input (forcing)
                            airtemp,                           & ! intent(in): air temperature at some height above the surface (K)
                            airpres,                           & ! intent(in): air pressure at some height above the surface (Pa)
                            scalarO2air,                       & ! intent(in): atmospheric o2 concentration (Pa)
                            scalarCO2air,                      & ! intent(in): atmospheric co2 concentration (Pa)
                            scalarCanopySunlitPAR,             & ! intent(in): average absorbed par for sunlit leaves (w m-2)
                            scalarCanopyShadedPAR,             & ! intent(in): average absorbed par for shaded leaves (w m-2)
                            ! input (state and diagnostic variables)
                            scalarGrowingSeasonIndex,          & ! intent(in): growing season index (0=off, 1=on)
                            scalarFoliageNitrogenFactor,       & ! intent(in): foliage nitrogen concentration (1=saturated)
                            scalarTranspireLim,                & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                            scalarLeafResistance,              & ! intent(in): leaf boundary layer resistance (s m-1)
                            scalarVegetationTemp,              & ! intent(in): temperature of the vegetation canopy (K)
                            scalarSatVP_VegTemp,               & ! intent(in): saturation vapor pressure at the temperature of the veg canopy (Pa)
                            scalarVP_CanopyAir,                & ! intent(in): canopy air vapor pressure (Pa)
                            bp, mp  ,c3psn, avcmx, vcmx25, ako, ko_25, akc,kc_25,qe25,folnmx,&
                            rgl, rsmin, rsmax, topt, hs, &
                            ! output
                            scalarStomResistSunlit,            & ! intent(out): stomatal resistance for sunlit leaves (s m-1)
                            scalarStomResistShaded,            & ! intent(out): stomatal resistance for shaded leaves (s m-1)
                            scalarPhotosynthesisSunlit,        & ! intent(out): sunlit photosynthesis (umolco2 m-2 s-1)
                            scalarPhotosynthesisShaded)!,        & ! intent(out): shaded photosynthesis (umolco2 m-2 s-1)
                            !err,cmessage                       ) ! intent(out): error control
    !  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
  
                            ! print*, 3437, scalarStomResistSunlit

    ! *******************************************************************************************************************************************
  
    ! error check
    ! case default; err=20; message=trim(message)//'unable to identify option for stomatal resistance'; return
  
    ! *******************************************************************************************************************************************
  
   end select  ! (identifying option for stomatal resistance)
  end subroutine
  
  
  end module stomResist_module
  