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

module summaSolve4homegrown_module

! data types
USE nrtype

! access the global print flag
USE globalData,only:globalPrintFlag

! access missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number
USE globalData,only:quadMissing     ! missing quadruple precision number

! access named variables to describe the form and structure of the matrices used in the numerical solver
USE globalData,only: ixFullMatrix   ! named variable for the full Jacobian matrix
USE globalData,only: ixBandMatrix   ! named variable for the band diagonal matrix
USE globalData,only: ku             ! number of super-diagonal bands
USE globalData,only: kl             ! number of sub-diagonal bands
USE globalData,only: nBands         ! length of the leading dimension of the band diagonal matrix
USE globalData,only: iJac1          ! first layer of the Jacobian to print
USE globalData,only: iJac2          ! last layer of the Jacobian to print

! indices of elements of data structure
USE var_lookup,only:iLookFLUX       ! named variables for structure elements
USE var_lookup,only:iLookPROG       ! named variables for structure elements
USE var_lookup,only:iLookPARAM      ! named variables for structure elements
USE var_lookup,only:iLookINDEX      ! named variables for structure elements
USE var_lookup,only:iLookDECISIONS  ! named variables for elements of the decision structure

USE multiconst,only:&
                    iden_water      ! intrinsic density of liquid water    (kg m-3)

! provide access to the derived types to define the data structures
USE data_types,only:&
                    var_i,                        & ! data vector (i4b)
                    var_d,                        & ! data vector (rkind)
                    var_ilength,                  & ! data vector with variable length dimension (i4b)
                    var_dlength,                  & ! data vector with variable length dimension (rkind)
                    zLookup,                      & ! lookup tables
                    model_options,                & ! defines the model decisions
                    in_type_computJacob,          & ! class for computJacob arguments
                    out_type_computJacob,         & ! class for computJacob arguments
                    in_type_lineSearchRefinement, & ! class for lineSearchRefinement arguments
                    out_type_lineSearchRefinement,& ! class for lineSearchRefinement arguments
                    in_type_summaSolve4homegrown, & ! class for summaSolve4homegrown arguments
                    io_type_summaSolve4homegrown, & ! class for summaSolve4homegrown arguments
                    out_type_summaSolve4homegrown   ! class for summaSolve4homegrown arguments


! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:       &
  qbaseTopmodel,                  & ! TOPMODEL-ish baseflow parameterization
  bigBucket,                      & ! a big bucket (lumped aquifer model)
  noExplicit                        ! no explicit groundwater parameterization

! look-up values for method used to compute derivative
USE mDecisions_module,only:       &
  numerical,                      & ! numerical solution
  analytical                        ! analytical solution

implicit none
private
! public :: summaSolve4homegrown
! public :: refine_Newton_step
! public :: checkConv
! contains


end module summaSolve4homegrown_module
