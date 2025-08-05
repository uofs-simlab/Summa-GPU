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

program summa_driver
  ! **** Driver program for SUMMA simulations ****

  ! * module access *
  ! data types
  USE nrtype                                                  ! variable types, etc.
  USE summa_type, only: summa1_type_dec                       ! master summa data type
  ! subroutines and functions: model setup
  USE summa_init, only: summa_initialize                      ! used to allocate/initialize summa data structures
  USE summa_setup, only: summa_paramSetup                     ! used to initialize parameter data structures (e.g. vegetation and soil parameters)
  USE summa_restart, only: summa_readRestart                  ! used to read restart data and reset the model state
  ! subroutines and functions: model simulation
  USE summa_forcing, only: summa_readForcing                  ! used to read forcing data
  USE summa_modelRun, only: summa_runPhysics                  ! used to run the summa physics for one time step
  USE summa_writeOutput, only: summa_writeOutputFiles         ! used to write the summa output files
  ! utility functions
  USE summa_util, only: stop_program                          ! used to stop the summa program (with errors)
  USE summa_util, only: handle_err                            ! used to process errors
  ! global data
  USE globalData, only: numtim                                ! number of model time steps
  USE globalData, only: print_step_freq

  use device_data_types
  use initialize_device
  ! OpenWQ coupling

  implicit none

  ! * driver variables *


  ! *****************************************************************************
  ! * variable definitions
  ! *****************************************************************************
  ! define the master summa data structure
  type(summa1_type_dec), allocatable :: summa1_struc(:)
  ! define parameters for the model simulation
  integer(i4b), parameter            :: n=1                        ! number of instantiations
  ! define timing information
  integer(i4b)                       :: modelTimeStep              ! index of model time step
  ! error control
  integer(i4b)                       :: err=0                      ! error code
  character(len=1024)                :: message=''                 ! error message

  ! Initialize
  call initialize_summa_driver

  ! Update
  call update_summa_driver

  ! Finalize
  call finalize_summa_driver

contains

  subroutine initialize_summa_driver
    integer(i4b) :: nSnow,nSoil,nLayers
   ! *** Initial operations for SUMMA driver program ***

   ! allocate space for the master summa structure
   allocate(summa1_struc(n), stat=err)
   if (err/=0) call stop_program(1, 'problem allocating master summa structure')

   ! declare and allocate summa data structures and initialize model state to known values
   call summa_initialize(summa1_struc(n), err, message)
   call handle_err(err, message)

   ! initialize parameter data structures (e.g. vegetation and soil parameters)
   call summa_paramSetup(summa1_struc(n), err, message)
   call handle_err(err, message)

   ! read restart data and reset the model state
   call summa_readRestart(summa1_struc(n), err, message)
   call handle_err(err, message)
 call finalize_device_indx_data(summa1_struc(n)%indxStruct_d,summa1_struc(n)%indxStruct,summa1_struc(n)%nGRU)

         nSnow = summa1_struc(n)%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
        nSoil = summa1_struc(n)%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nSoil)%dat(1)
            nLayers = summa1_struc(n)%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nLayers)%dat(1)

!  call finalize_device_forc_data(summa1_struc(n)%forcStruct_d,summa1_struc(n)%forcStruct,summa1_struc(n)%nGRU)
! call finalize_device_diag_data(summa1_struc(n)%diagStruct_d,summa1_struc(n)%diagStruct,nSnow,nSoil,nLayers,summa1_struc(n)%nGRU)
! call finalize_device_prog_data(summa1_struc(n)%progStruct_d,summa1_struc(n)%progStruct,nLayers,nSoil,summa1_struc(n)%nGRU)
  !  call finalize_device_flux_data(summa1_struc(n)%fluxStruct_d,summa1_struc(n)%fluxStruct,nSnow,nSoil,summa1_struc(n)%nGRU)
!  call finalize_device_bvar_data(summa1_struc(n)%bvarStruct_d,summa1_struc(n)%bvarStruct,summa1_struc(n)%nGRU)

  end subroutine initialize_summa_driver

  subroutine update_summa_driver
    integer(i4b) :: nSnow,nSoil,nLayers
   ! *** Update operations for SUMMA driver program ***
    nSnow = summa1_struc(n)%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
        nSoil = summa1_struc(n)%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nSoil)%dat(1)
            nLayers = summa1_struc(n)%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nLayers)%dat(1)

   ! loop through time
   do modelTimeStep=1,numtim

     ! read model forcing data
     call summa_readForcing(modelTimeStep, summa1_struc(n), err, message)
     call handle_err(err, message)
  
    !  if (mod(modelTimeStep, print_step_freq) == 0) then
       print *, 'step ---> ', modelTimeStep, numtim
    !  end if
 
     ! run the summa physics for one time step
     call summa_runPhysics(modelTimeStep, summa1_struc(n), err, message)
     call handle_err(err, message)
    !  if (modelTimeStep .eq. numtim) then
! call finalize_device_indx_data(summa1_struc(n)%indxStruct_d,summa1_struc(n)%indxStruct,summa1_struc(n)%nGRU)

!         nSnow = summa1_struc(n)%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nSnow)%dat(1)
!        nSoil = summa1_struc(n)%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nSoil)%dat(1)
!            nLayers = summa1_struc(n)%indxStruct%gru(1)%hru(1)%var(iLookINDEX%nLayers)%dat(1)

! call finalize_device_forc_data(summa1_struc(n)%forcStruct_d,summa1_struc(n)%forcStruct,summa1_struc(n)%nGRU)
!call finalize_device_diag_data(summa1_struc(n)%diagStruct_d,summa1_struc(n)%diagStruct,nSnow,nSoil,nLayers,summa1_struc(n)%nGRU)
!call finalize_device_prog_data(summa1_struc(n)%progStruct_d,summa1_struc(n)%progStruct,nLayers,nSoil,summa1_struc(n)%nGRU)
!   call finalize_device_flux_data(summa1_struc(n)%fluxStruct_d,summa1_struc(n)%fluxStruct,nSnow,nSoil,summa1_struc(n)%nGRU)
! call finalize_device_bvar_data(summa1_struc(n)%bvarStruct_d,summa1_struc(n)%bvarStruct,summa1_struc(n)%nGRU)
    !  end if
     ! write the model output
     call summa_writeOutputFiles(modelTimeStep, summa1_struc(n), err, message)
     call handle_err(err, message)
  
   end do  ! end looping through time
   call deallocate_device_diag_data(summa1_struc(n)%diagStruct_d)
call deallocate_device_forc_data(summa1_struc(n)%forcStruct_d)
call deallocate_device_prog_data(summa1_struc(n)%progStruct_d)
call deallocate_device_flux_data(summa1_struc(n)%fluxStruct_d)
 call deallocate_device_bvar_data(summa1_struc(n)%bvarStruct_d)
 call deallocate_device_indx_data(summa1_struc(n)%indxStruct_d)

  end subroutine update_summa_driver

  subroutine finalize_summa_driver
   ! *** Final operations for SUMMA driver program ***
   ! successful end
   call stop_program(0, 'finished simulation successfully.')

   ! to prevent exiting before HDF5 has closed
   call sleep(2)
  end subroutine finalize_summa_driver

end program summa_driver
