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

module eval8summa_module

! data types
USE nrtype

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number
USE globalData,only:quadMissing     ! missing quadruple precision number

! named variables to describe the state variable type
USE globalData,only:iname_watLayer  ! named variable defining the total water state variable for snow+soil layers
USE globalData,only:iname_liqLayer  ! named variable defining the liquid  water state variable for snow+soil layers
USE globalData,only:iname_matLayer  ! named variable defining the total water matric potential state variable for soil layers
USE globalData,only:iname_lmpLayer  ! named variable defining the liquid water matric potential state variable for soil layers

! constants
USE multiconst,only:&
                    LH_fus,     & ! latent heat of fusion                (J kg-1)
                    iden_water, & ! intrinsic density of liquid water    (kg m-3)
                    gravity,    & ! gravitational acceleration           (m s-2)
                    Tfreeze       ! freezing point of pure water         (K)


! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,        & ! data vector (i4b)
                    var_d,        & ! data vector (rkind)
                    var_ilength,  & ! data vector with variable length dimension (i4b)
                    var_dlength,  & ! data vector with variable length dimension (rkind)
                    zLookup,      & ! data vector with variable length dimension (rkind)
                    model_options   ! defines the model decisions

! indices that define elements of the data structures
USE var_lookup,only:iLookDECISIONS               ! named variables for elements of the decision structure
USE var_lookup,only:iLookPARAM                   ! named variables for structure elements
USE var_lookup,only:iLookPROG                    ! named variables for structure elements
USE var_lookup,only:iLookINDEX                   ! named variables for structure elements
USE var_lookup,only:iLookDIAG                    ! named variables for structure elements
USE var_lookup,only:iLookFLUX                    ! named variables for structure elements
USE var_lookup,only:iLookDERIV                   ! named variables for structure elements

! look-up values for the choice of variable in energy equations (BE residual or IDA state variable)
USE mDecisions_module,only:  &
 closedForm,                 & ! use temperature with closed form heat capacity
 enthalpyFormLU,             & ! use enthalpy with soil temperature-enthalpy lookup tables
 enthalpyForm                  ! use enthalpy with soil temperature-enthalpy analytical solution

! look-up values for the numerical method
USE mDecisions_module,only:  &
 homegrown,                  & ! homegrown backward Euler solution based on concepts from numerical recipes
 kinsol,                     & ! SUNDIALS backward Euler solution using Kinsol
 ida                           ! SUNDIALS solution using IDA

implicit none
private
#ifdef SUNDIALS_ACTIVE
#endif
public::imposeConstraints

contains

! ***************************************************************************************************************************************
! public subroutine imposeConstraints: impose solution constraints 
!   This is simple error control to reduce possible temperature increments, cross over freezing point events, and keep the state feasible
!   Imposed after the internal call of KINSOL incrementing the linesearch
! ***************************************************************************************************************************************
subroutine imposeConstraints(model_decisions,indx_data, prog_data, mpar_data, stateVec, stateVecPrev,&
    nState, nSoil, nSnow, message, err)
  ! external functions
  USE snow_utils_module,only:fracliquid     ! compute the fraction of liquid water at a given temperature (snow)
  USE soil_utils_module,only:crit_soilT     ! compute the critical temperature below which ice exists
  USE soil_utils_module,only:matricHead     ! compute the matric head based on volumetric water content
  USE soil_utils_module,only:volFracLiq     ! compute volumetric fraction of liquid water

  implicit none
  
  type(model_options),intent(in)           :: model_decisions(:)  ! model decisions
  type(var_ilength),intent(in)             :: indx_data           ! indices defining model states and layers
  type(var_dlength),intent(in)             :: prog_data           ! prognostic variables for a local HRU
  type(var_dlength),intent(in)             :: mpar_data           ! model parameters
  real(rkind), intent(inout)               :: stateVec(:)         ! state vector
  real(rkind), intent(in)                  :: stateVecPrev(:)     ! previous state vector
  integer(i4b),intent(in)                  :: nState              ! total number of state variables
  integer(i4b),intent(in)                  :: nSoil               ! number of soil layers
  integer(i4b),intent(in)                  :: nSnow               ! number of snow layers
  integer(i4b),intent(out)                 :: err                 ! error code
  character(len=256),intent(out)           :: message             ! error message
  ! -----------------------------------------------------------------------------------------------------
  ! temporary variables for model constraints
  real(qp),dimension(nState)               :: xInc                       ! iteration increment
  real(rkind)                              :: scalarTemp                 ! temperature of an individual snow layer (K)
  real(rkind)                              :: scalarIce                  ! volumetric ice content of an individual layer (-)
  real(rkind)                              :: scalarLiq                  ! volumetric liquid water content of an individual layer (-)
  real(rkind)                              :: xPsi00                     ! matric head after applying the iteration increment (m)
  real(rkind)                              :: TcSoil                     ! critical point when soil begins to freeze (K)
  real(rkind)                              :: critDiff                   ! temperature difference from critical (K)
  real(rkind)                              :: epsT                       ! small interval above/below critical (K)
  real(rkind)                              :: zMaxTempIncrement          ! maximum temperature increment (K)
  real(rkind)                              :: zMaxMatricIncrement        ! maximum matric head increment (m)
  real(rkind)                              :: xConst                     ! constant in the freezing curve function (m K-1)
  real(rkind)                              :: mLayerPsiLiq               ! liquid water matric potential (m)
  real(rkind)                              :: vGn_m(nSoil)               ! van Genutchen "m" parameter (-)
  real(rkind)                              :: effSat                     ! effective saturation (-)
  real(rkind)                              :: avPore                     ! available pore space (-)
  ! indices of model state variables
  integer(i4b)                             :: iState                     ! index of state within a specific variable type
  integer(i4b)                             :: ixNrg,ixLiq                ! index of energy and mass state variables in full state vector
  ! indices of model layers
  integer(i4b)                             :: iLayer                     ! index of model layer
  ! choice of constraints to impose
  logical(lgt)                             :: small_delTemp              ! flag to constain temperature change to be less than zMaxTempIncrement
  logical(lgt)                             :: small_delMatric            ! flag to constain matric head change to be less than zMaxMatricIncrement
  logical(lgt)                             :: detect_events              ! flag to do freezing point event detection and cross-over with epsT
  logical(lgt)                             :: water_bounds               ! flag to force water to not go above or below physical bounds
  
  ! -----------------------------------------------------------------------------------------------------
  ! association to variables in the data structures
  associate(&
    ! model decisions
    ixNumericalMethod  => model_decisions(iLookDECISIONS%num_method)%iDecision ,& ! intent(in):  [i4b]   choice of numerical solver
    ! indices of model state variables
    ixNrgOnly          => indx_data%var(iLookINDEX%ixNrgOnly)%dat              ,& ! intent(in): [i4b(:)] list of indices in the state subset for energy states
    ixHydOnly          => indx_data%var(iLookINDEX%ixHydOnly)%dat              ,& ! intent(in): [i4b(:)] list of indices in the state subset for hydrology states
    ixMatOnly          => indx_data%var(iLookINDEX%ixMatOnly)%dat              ,& ! intent(in): [i4b(:)] list of indices in the state subset for matric head states
    ixMassOnly         => indx_data%var(iLookINDEX%ixMassOnly)%dat             ,& ! intent(in): [i4b(:)] list of indices in the state subset for canopy storage states
    ixHydType          => indx_data%var(iLookINDEX%ixHydType)%dat              ,& ! intent(in): [i4b(:)] index of the type of hydrology states in snow+soil domain
    ixStateType_subset => indx_data%var(iLookINDEX%ixStateType_subset)%dat     ,& ! intent(in): [i4b(:)] named variables defining the states in the subset
    ! indices for specific state variables
    ixCasNrg           => indx_data%var(iLookINDEX%ixCasNrg)%dat(1)            ,& ! intent(in): [i4b]    index of canopy air space energy state variable
    ixVegNrg           => indx_data%var(iLookINDEX%ixVegNrg)%dat(1)            ,& ! intent(in): [i4b]    index of canopy energy state variable
    ixVegHyd           => indx_data%var(iLookINDEX%ixVegHyd)%dat(1)            ,& ! intent(in): [i4b]    index of canopy hydrology state variable (mass)
    ixTopNrg           => indx_data%var(iLookINDEX%ixTopNrg)%dat(1)            ,& ! intent(in): [i4b]    index of upper-most energy state in the snow-soil subdomain
    ixTopHyd           => indx_data%var(iLookINDEX%ixTopHyd)%dat(1)            ,& ! intent(in): [i4b]    index of upper-most hydrology state in the snow-soil subdomain
    ! vector of energy indices for the snow and soil domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnowSoilNrg      => indx_data%var(iLookINDEX%ixSnowSoilNrg)%dat          ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow+soil domain
    ixSnowOnlyNrg      => indx_data%var(iLookINDEX%ixSnowOnlyNrg)%dat          ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the snow domain
    ixSoilOnlyNrg      => indx_data%var(iLookINDEX%ixSoilOnlyNrg)%dat          ,& ! intent(in): [i4b(:)] index in the state subset for energy state variables in the soil domain
    ! vector of hydrology indices for the snow and soil domains
    ! NOTE: states not in the subset are equal to integerMissing
    ixSnowSoilHyd      => indx_data%var(iLookINDEX%ixSnowSoilHyd)%dat          ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow+soil domain
    ixSnowOnlyHyd      => indx_data%var(iLookINDEX%ixSnowOnlyHyd)%dat          ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the snow domain
    ixSoilOnlyHyd      => indx_data%var(iLookINDEX%ixSoilOnlyHyd)%dat          ,& ! intent(in): [i4b(:)] index in the state subset for hydrology state variables in the soil domain
    ! number of state variables of a specific type
    nSnowSoilNrg       => indx_data%var(iLookINDEX%nSnowSoilNrg )%dat(1)       ,& ! intent(in): [i4b]    number of energy state variables in the snow+soil domain
    nSnowOnlyNrg       => indx_data%var(iLookINDEX%nSnowOnlyNrg )%dat(1)       ,& ! intent(in): [i4b]    number of energy state variables in the snow domain
    nSoilOnlyNrg       => indx_data%var(iLookINDEX%nSoilOnlyNrg )%dat(1)       ,& ! intent(in): [i4b]    number of energy state variables in the soil domain
    nSnowSoilHyd       => indx_data%var(iLookINDEX%nSnowSoilHyd )%dat(1)       ,& ! intent(in): [i4b]    number of hydrology variables in the snow+soil domain
    nSnowOnlyHyd       => indx_data%var(iLookINDEX%nSnowOnlyHyd )%dat(1)       ,& ! intent(in): [i4b]    number of hydrology variables in the snow domain
    nSoilOnlyHyd       => indx_data%var(iLookINDEX%nSoilOnlyHyd )%dat(1)       ,& ! intent(in): [i4b]    number of hydrology variables in the soil domain
    ! snow parameters
    snowfrz_scale      => mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)       ,& ! intent(in):  [dp]    scaling parameter for the snow freezing curve (K-1)
    ! soil parameters
    theta_sat          => mpar_data%var(iLookPARAM%theta_sat)%dat              ,& ! intent(in): [dp(:)]  soil porosity (-)
    theta_res          => mpar_data%var(iLookPARAM%theta_res)%dat              ,& ! intent(in): [dp(:)]  residual volumetric water content (-)
    vGn_n              => mpar_data%var(iLookPARAM%vGn_n)%dat                  ,& ! intent(in):  [dp(:)]  van Genutchen "n" parameter (-)
    vGn_alpha          => mpar_data%var(iLookPARAM%vGn_alpha)%dat              ,& ! intent(in):  [dp(:)]  van Genutchen "alpha" parameter (m-1)
    ! state variables at the start of the time step
    mLayerMatricHead   => prog_data%var(iLookPROG%mLayerMatricHead)%dat        ,& ! intent(in): [dp(:)] matric head (m)
    mLayerVolFracIce   => prog_data%var(iLookPROG%mLayerVolFracIce)%dat         & ! intent(in): [dp(:)] volumetric fraction of ice (-)
    ) ! associating variables with indices of model state variables
    ! -----------------------------------------------------------------------------------------------------

    ! initialize error control
    err=0; message='imposeConstraints/'

    ! calculate proposed increment in state vector
    xInc(1:nState) = stateVec(1:nState)*1._qp - stateVecPrev(1:nState)*1._qp
  
    ! identify which constraints to impose
    select case(ixNumericalMethod)
    case(ida); err=20; message=trim(message)//'should not be imposing constraints for IDA solver'; return
      case(kinsol)
        small_delTemp       = .true.      ! flag to constain temperature change to be less than zMaxTempIncrement
        zMaxTempIncrement   = 10._rkind   ! maximum temperature increment (K)
        small_delMatric     = .true.      ! flag to constain matric head change to be less than zMaxMatricIncrement
        zMaxMatricIncrement = 10._rkind   ! maximum matric head increment (m)
        detect_events       = .true.      ! flag to do freezing point event detection and cross-over with epsT, works best if on
        epsT                = 1.e-7_rkind ! small interval above/below critical (K), works better if larger
        water_bounds        = .true.      ! flag to force water bounds, works best if on
      case(homegrown)
        small_delTemp       = .true.      ! flag to constain temperature change to be less than zMaxTempIncrement
        zMaxTempIncrement   = 10._rkind   ! maximum temperature increment (K)
        small_delMatric     = .true.      ! flag to constain matric head change to be less than zMaxMatricIncrement
        zMaxMatricIncrement = 10._rkind   ! maximum matric head increment (m)
        detect_events       = .true.      ! flag to do freezing point event detection and cross-over with epsT
        epsT                = 1.e-7_rkind ! small interval above/below critical (K)
        water_bounds        = .true.      ! flag to force water bounds
      case default; err=20; message=trim(message)//'expect num_method to be ida, kinsol, or homegrown (or itertive, which is homegrown)'; return
    end select
    
    vGn_m = 1._rkind - 1._rkind/vGn_n

    ! ** limit temperature increment to zMaxTempIncrement
    ! NOTE: this can cause problems especially from a cold start when far from the solution
    if(small_delTemp)then
      if(size(ixNrgOnly)>0)then
        ! loop through snow+soil layers
        do iState=1,size(ixNrgOnly)
          ! define index of the energy state variable within the state subset
          ixNrg = ixNrgOnly(iState)
          ! place constraint for temperature
          if(abs(xInc(ixNrg)) > zMaxTempIncrement) xInc(ixNrg) = sign(zMaxTempIncrement, xInc(ixNrg))
        end do ! (loop through snow+soil layers)
      endif
    endif ! (small temperature change)

    ! ** limit soil water (matric head) increment to zMaxMatricIncrement if starting positive
    if(small_delMatric)then
      if(size(ixMatOnly)>0)then
        ! loop through soil layers
        do iState=1,size(ixMatOnly)
          ! define index of the hydrology state variable within the state subset
          ixLiq = ixMatOnly(iState)
          ! place constraint for matric head
          if(xInc(ixLiq) > zMaxMatricIncrement .and. stateVecPrev(ixLiq) > 0._rkind) xInc(ixLiq) = zMaxMatricIncrement
        end do ! (loop through soil layers)
      endif
    endif ! (small matric head change)

    ! ** stop just above or just below the freezing point if crossing
    if(detect_events)then

      ! crossing freezing point event for vegetation
      if(ixVegNrg/=integerMissing)then
        ! initialize
        critDiff = Tfreeze - stateVecPrev(ixVegNrg)
        ! initially frozen (T < Tfreeze)
        if(critDiff > 0._rkind)then ! (check crossing above zero)
          if(xInc(ixVegNrg) > critDiff) xInc(ixVegNrg) = critDiff + epsT  ! constrained temperature increment (K)
        ! initially unfrozen (T > Tfreeze)
        else ! (check crossing below zero)
          if(xInc(ixVegNrg) < critDiff) xInc(ixVegNrg) = critDiff - epsT  ! constrained temperature increment (K)
        end if  ! (switch between initially frozen and initially unfrozen)
      endif  ! if the state variable for canopy temperature is included within the state subset

      ! crossing freezing point event for snow, keep it below freezing
      if(nSnowOnlyNrg > 0)then
        do iLayer=1,nSnow 
          ! check if energy state is included
          if(ixSnowOnlyNrg(iLayer)==integerMissing) cycle
          ! check temperatures, and, if necessary, scale iteration increment
          iState = ixSnowOnlyNrg(iLayer)
          ! constrained temperature increment (K) -- simplified bi-section
          if(stateVecPrev(iState) + xInc(iState) > Tfreeze) xInc(iState) = 0.5_rkind*(Tfreeze - stateVecPrev(iState) )
        end do ! (loop through snow layers)
      endif  ! (if there are state variables for energy in the snow domain)

      ! crossing freezing point event for soil
      if(nSoilOnlyNrg>0)then
        do iLayer=1,nSoil
          ! check if energy state is included
          if(ixSoilOnlyNrg(iLayer)==integerMissing) cycle
          ! define index of the state variables within the state subset
          ixNrg = ixSoilOnlyNrg(iLayer)
          ixLiq = ixSoilOnlyHyd(iLayer)
          ! get the matric potential of total water
          if(ixLiq/=integerMissing)then
            select case( ixStateType_subset( ixSoilOnlyHyd(iLayer) ) )     
              case(iname_lmpLayer)
                effSat = volFracLiq(stateVecPrev(ixLiq) + xInc(ixLiq),vGn_alpha(iLayer),0._rkind,1._rkind,vGn_n(iLayer),vGn_m(iLayer))  ! effective saturation
                avPore = theta_sat(iLayer) - mLayerVolFracIce(iLayer+nSnow) - theta_res(iLayer)  ! available pore space
                scalarLiq = effSat*avPore + theta_res(iLayer)
                xPsi00 = matricHead(scalarLiq + mLayerVolFracIce(iLayer+nSnow),vGn_alpha(iLayer),theta_res(iLayer),theta_sat(iLayer),vGn_n(iLayer),vGn_m(iLayer))
              case(iname_matLayer); xPsi00 = stateVecPrev(ixLiq) + xInc(ixLiq) ! only true if using iname_matLayer, otherwise may want to fix this
              case(iname_watLayer); xPsi00 = matricHead(stateVecPrev(ixLiq) + xInc(ixLiq),vGn_alpha(iLayer),theta_res(iLayer),theta_sat(iLayer),vGn_n(iLayer),vGn_m(iLayer))
              case(iname_liqLayer) 
                xPsi00 = matricHead(mLayerVolFracIce(iLayer+nSnow) + stateVecPrev(ixLiq) + xInc(ixLiq),vGn_alpha(iLayer),theta_res(iLayer),theta_sat(iLayer),vGn_n(iLayer),vGn_m(iLayer))
              case default; err=20; message=trim(message)//'expect ixStateType_subset to be iname_matLayer, iname_lmpLayer, iname_watLayer, or iname_liqLayer for soil hydrology'; return
            end select
          else
            xPsi00 = mLayerMatricHead(iLayer)
          endif
          ! identify the critical point when soil begins to freeze (TcSoil)
          TcSoil = crit_soilT(xPsi00)
          ! get the difference from the current state and the crossing point (K)
          critDiff = TcSoil - stateVecPrev(ixNrg)
          ! initially frozen (T < TcSoil)
          if(critDiff > 0._rkind)then ! (check crossing above zero)
            if(xInc(ixNrg) > critDiff) xInc(ixNrg) = critDiff + epsT  ! set iteration increment to slightly above critical temperature
          ! initially unfrozen (T > TcSoil)
          else ! (check crossing below zero)
            if(xInc(ixNrg) < critDiff) xInc(ixNrg) = critDiff - epsT  ! set iteration increment to slightly below critical temperature
          endif ! (switch between initially frozen and initially unfrozen)
        end do ! (loop through soil layers)
      endif ! (if there are both energy and liquid water state variables)

    endif ! (detect events)

    ! ** ensure water is within bounds
    if(water_bounds)then

      ! impose positivity for canopy liquid water
      if(ixVegHyd/=integerMissing)then
        ! constrained iteration increment (K) -- simplified bi-section
        if(stateVecPrev(ixVegHyd) + xInc(ixVegHyd) < 0._rkind) xInc(ixVegHyd) = -0.5_rkind*stateVecPrev(ixVegHyd)
      endif ! (if the state variable for canopy water is included within the state subset)
          
      ! impose bounds for snow water, change in total water is only due to liquid flux
      if(nSnowOnlyHyd>0)then
        ! loop through snow layers
        do iLayer=1,nSnow
           ! check if the layer is included
          if(ixSnowOnlyHyd(iLayer)==integerMissing) cycle
          if(ixSnowOnlyNrg(iLayer)/=integerMissing)then
            ! get the layer temperature (from stateVecPrev if ixSnowOnlyNrg(iLayer) is within the state vector
            scalarTemp = stateVecPrev( ixSnowOnlyNrg(iLayer) )
          else ! get the layer temperature from the last update
            scalarTemp = prog_data%var(iLookPROG%mLayerTemp)%dat(iLayer)
          endif
          ! get the volumetric fraction of liquid water and ice
          select case( ixStateType_subset( ixSnowOnlyHyd(iLayer) ) )
            case(iname_watLayer); scalarLiq = fracliquid(scalarTemp,mpar_data%var(iLookPARAM%snowfrz_scale)%dat(1)) * stateVecPrev(ixSnowOnlyHyd(iLayer))
            case(iname_liqLayer); scalarLiq = stateVecPrev(ixSnowOnlyHyd(iLayer))
            case default; err=20; message=trim(message)//'expect ixStateType_subset to be iname_watLayer or iname_liqLayer for snow hydrology'; return
          end select
          scalarIce = merge(stateVecPrev(ixSnowOnlyHyd(iLayer)) - scalarLiq,mLayerVolFracIce(iLayer), ixHydType(iLayer)==iname_watLayer)
          ! checking if drain more than what is available or add more than possible, constrained iteration increment -- simplified bi-section
          if(-xInc(ixSnowOnlyHyd(iLayer)) > scalarLiq) then
            xInc(ixSnowOnlyHyd(iLayer)) = -0.5_rkind*scalarLiq
          elseif(xInc(ixSnowOnlyHyd(iLayer)) > 1._rkind - scalarIce - scalarLiq)then
            xInc(ixSnowOnlyHyd(iLayer)) = 0.5_rkind*(1._rkind - scalarIce - scalarLiq)
          endif
        end do ! (looping through snow layers)
      endif ! (if there are state variables for liquid water in the snow domain)
       
      ! impose bounds for soil water, change in total water is only due to liquid flux
      if(nSoilOnlyHyd>0)then
        ! loop through soil layers
        do iLayer=1,nSoil
          ! check if the layer is included
          if(ixSoilOnlyHyd(iLayer)==integerMissing) cycle
          if(ixHydType(iLayer+nSnow)==iname_watLayer .or. ixHydType(iLayer+nSnow)==iname_liqLayer)then
            ! get the volumetric fraction of liquid water and ice
            select case( ixStateType_subset( ixSoilOnlyHyd(iLayer) ) )
              case(iname_watLayer)
                xPsi00 = matricHead(stateVecPrev(ixSoilOnlyHyd(iLayer)),vGn_alpha(iLayer),theta_res(iLayer),theta_sat(iLayer),vGn_n(iLayer),vGn_m(iLayer))
                ! get the layer temperature
                if(ixSoilOnlyNrg(iLayer)/=integerMissing)then
                  scalarTemp = stateVecPrev( ixSoilOnlyNrg(iLayer) )
                else
                  scalarTemp = prog_data%var(iLookPROG%mLayerTemp)%dat(iLayer+nSnow)
                endif
                ! identify the critical point when soil begins to freeze (TcSoil)
                TcSoil = crit_soilT(xPsi00)
                ! get the volumetric fraction of liquid water and ice
                if(scalarTemp < TcSoil)then 
                  xConst       = LH_fus/(gravity*Tfreeze)
                  mLayerPsiLiq = xConst*(scalarTemp - Tfreeze)
                  scalarLiq = volFracLiq(mLayerPsiLiq,vGn_alpha(iLayer),theta_res(iLayer),theta_sat(iLayer),vGn_n(iLayer),vGn_m(iLayer))
                else
                  scalarLiq = stateVecPrev(ixSoilOnlyHyd(iLayer))
                end if  ! (check if soil is partially frozen)
              case(iname_liqLayer); scalarLiq = stateVecPrev(ixSoilOnlyHyd(iLayer))
            end select
            scalarIce = merge(stateVecPrev(ixSoilOnlyHyd(iLayer)) - scalarLiq,mLayerVolFracIce(iLayer+nSnow), ixHydType(iLayer)==iname_watLayer)
            ! checking if drain more than what is available or add more than possible, constrained iteration increment -- simplified bi-section
            if(-xInc(ixSoilOnlyHyd(iLayer)) > scalarLiq - theta_res(iLayer))then 
              xInc(ixSoilOnlyHyd(iLayer)) = -0.5_rkind*(scalarLiq - theta_res(iLayer))
            elseif(xInc(ixSoilOnlyHyd(iLayer)) > theta_sat(iLayer) - scalarIce - scalarLiq)then
              xInc(ixSoilOnlyHyd(iLayer)) = 0.5_rkind*(theta_sat(iLayer) - scalarIce - scalarLiq)
            endif
          endif ! (if the state variable is not matric head)
        end do ! (looping through soil layers)
      endif ! (if there are state variables for liquid water in the soil domain)

    endif ! (water bounds)  

    ! Update the state vector with the modified iteration increment
    stateVec(:) = stateVecPrev(:) + xInc(:)

  ! end association with variables with indices of model state variables
  end associate

end subroutine imposeConstraints


end module eval8summa_module
