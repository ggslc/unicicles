!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   bisicles_cdriver.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! WJS (1-30-12): The following (turning optimization off) is needed as a workaround for an
! xlf compiler bug, at least in IBM XL Fortran for AIX, V12.1 on bluefire
#ifdef CPRIBM
@PROCESS OPT(0)
#endif

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#ifdef BISICLES_CDRIVER
#include "cdriverconstants.h"

module bisicles_cdriver
  !Driver for BISICLES (parallel, external, AMR) dynamical core

  use glide_types
  use glimmer_physcon
  use glimmer_log
  use parallel
  implicit none
  real(kind=8) :: dx
  integer, dimension(1:2) :: dims, boxlo,boxhi
  !real(kind=8) arrays needed to manage sp<->dp conversion while some glide arrays are sp
  real(kind=8),dimension(:,:),pointer :: usrf_thck_data => null() !*FD upper surface thickness boundary data
  real(kind=8),dimension(:,:),pointer :: usrf_heat_data => null() !*FD upper surface heat boundary data
  real(kind=8),dimension(:,:),pointer :: topg_ibc => null() !*FD copy of initial topography data 
  real(kind=8),dimension(:,:),pointer :: thck_ibc => null() !*FD copy of initial thickness data
  real(kind=8),dimension(:,:),pointer :: tmp_data => null() !*FD general purpose workspace
  character(len=128) :: message
  !public bisicles_initialise, bisicles_init_state_diagnostic, bisicles_tstep, bisicles_finalise
  private message
contains

 
  subroutine bisicles_post(model)
    !post BISICLES operations common to initialization and timestepping
    !retrieve BISICLES data (e.g the surface elevation, bedrock elevation, and ice thickness)
    !and update the appropriate glide fields
    type(glide_global_type), intent(inout) :: model  ! glide model instance
    
    call f_bisicles_get_2d_data(model%options%bisicles_instance_id, &
         model%geometry%usrf, BISICLES_FIELD_SURFACE_ELEVATION, & 
         dx, dims, boxlo, boxhi) 
    call f_bisicles_get_2d_data(model%options%bisicles_instance_id, &
         model%geometry%topg, BISICLES_FIELD_BEDROCK_ELEVATION, & 
         dx, dims, boxlo, boxhi) 
    call f_bisicles_get_2d_data(model%options%bisicles_instance_id, &
         model%geometry%thck, BISICLES_FIELD_ICE_THICKNESS, & 
         dx, dims, boxlo, boxhi)

    write(*,*) ' bisicles_post:  range(model%geometry%usrf) = ', minval(model%geometry%usrf) , maxval(model%geometry%usrf)
    !call write_log(message)
    write(*,*) ' bisicles_post:  range(model%geometry%topg) = ', minval(model%geometry%topg) , maxval(model%geometry%topg)
    !call write_log(message)
    write(*,*) ' bisicles_post:  range(model%geometry%thck) = ', minval(model%geometry%thck) , maxval(model%geometry%thck)
    !call write_log(message)
   
  
    !calved mass
    !call f_bisicles_get_2d_data(model%options%bisicles_instance_id, &
    !     tmp_data , BISICLES_FIELD_CALVED_THICKNESS, & 
    !     dx, dims, boxlo, boxhi)
    
    !dp -> sp conversion
    model%climate%calving = tmp_data
    
    ! Choice of data depending of choice of heat eqn BC
    if (model%climate%temp_upper_bc_type == TEMP_BC_TEMPERATURE) then
       ! read the heat flux
       call f_bisicles_get_2d_data(model%options%bisicles_instance_id, &
            tmp_data , BISICLES_FIELD_SURFACE_HEAT_FLUX, & 
             dx, dims, boxlo, boxhi)
       model%climate%hflx = tmp_data
    else if (model%climate%temp_upper_bc_type == TEMP_BC_FLUX) then
        ! read the temperature
       call f_bisicles_get_2d_data(model%options%bisicles_instance_id, &
            tmp_data , BISICLES_FIELD_SURFACE_TEMPERATURE, & 
            dx, dims, boxlo, boxhi)
       model%climate%artm = tmp_data - trpt !BISICLES provides K
    else
       write (message,*) 'invalid model%climate%temp_upper_bc_type = ',model%climate%temp_upper_bc_type
       call write_log(message,GM_FATAL,__FILE__, __LINE__)
    end if
     
  end subroutine bisicles_post

  subroutine bisicles_pre(model)
    type(glide_global_type), intent(inout) :: model   ! model instance
    
    !pre BISICLES operations common to initialization and timestepping

    !real*8 copy of model%climate%acab, artm
    usrf_thck_data = model%climate%acab
 
      ! Choice of data depending of choice of heat eqn BC
    if (model%climate%temp_upper_bc_type == TEMP_BC_TEMPERATURE) then
       ! copy the temperature (sp -> dp) 
       usrf_heat_data = model%climate%artm + trpt !BISICLES needs K

       !we can't deal with usrf_heat_data > trpt
       where  (usrf_heat_data > trpt)
          usrf_heat_data = trpt
       end where
    else if (model%climate%temp_upper_bc_type == TEMP_BC_FLUX) then
       ! copy the flux (sp -> dp) 
       usrf_heat_data = model%climate%hflx
    else
       write (message,*) 'invalid model%climate%temp_upper_bc_type = ',model%climate%temp_upper_bc_type
       call write_log(message,GM_FATAL,__FILE__, __LINE__)
    end if

  end subroutine bisicles_pre


  subroutine bisicles_finalise(model)
    type(glide_global_type), intent(inout) :: model   ! model instance
    call f_bisicles_free_instance (model%options%bisicles_instance_id)
    deallocate(usrf_thck_data)
    deallocate(usrf_heat_data)
    deallocate(tmp_data)
    deallocate(thck_ibc)
    deallocate(topg_ibc)
  end subroutine bisicles_finalise

  subroutine bisicles_initialise(model)
    type(glide_global_type), intent(inout) :: model   ! model instance
    integer l
    
    l = 1 + len_trim(model%options%dycore_input_file) ! needed by C-side to NULL terminate string
    call write_log("creating bisicles instance");
    call f_bisicles_new_instance (model%options%bisicles_instance_id, &
         model%options%dycore_input_file//char(0), l, comm )
    call write_log("bisicles instance created");	
    !set up uniform grid meta data
    !pull the domain decomposition scheme out of the parallel module
    dims(1) = global_ewn
    dims(2) = global_nsn
    boxlo(1) = ewlb - 1 
    boxlo(2) = nslb - 1
    boxhi(1) = ewub - 1
    boxhi(2) = nsub - 1
    dx = get_dew(model)

    !model%climate%acab etc are single precision, but we need double and so maintain copies
    allocate(usrf_thck_data(model%general%ewn,model%general%nsn))
    allocate(usrf_heat_data(model%general%ewn,model%general%nsn))
    
    !a bit of workspace is nice
    allocate(tmp_data(model%general%ewn,model%general%nsn))

    call f_bisicles_set_2d_data(model%options%bisicles_instance_id, &
         usrf_thck_data, BISICLES_FIELD_SURFACE_FLUX, & 
         dx, dims, boxlo, boxhi) 

    call write_log("bisicles setting 2d upper BC");			
    !Handle upper surface temperature BC
    if (model%climate%temp_upper_bc_type == TEMP_BC_TEMPERATURE) then
    call write_log("temperature upper BC");				

       call f_bisicles_set_2d_data(model%options%bisicles_instance_id, &
            usrf_heat_data, BISICLES_FIELD_SURFACE_TEMPERATURE, & 
            dx, dims, boxlo, boxhi)

    else if (model%climate%temp_upper_bc_type == TEMP_BC_FLUX) then
    call write_log("flux upper BC");			
       call f_bisicles_set_2d_data(model%options%bisicles_instance_id, &
            usrf_heat_data, BISICLES_FIELD_SURFACE_HEAT_FLUX, & 
            dx, dims, boxlo, boxhi) 
    else
       write (message,*) 'invalid model%climate%temp_upper_bc_type = ',model%climate%temp_upper_bc_type
       call write_log(message,GM_FATAL,__FILE__, __LINE__)
    end if
    
    !set GIA data
    !bisicles needs the time derivative of the bedrock elevation 
    !(since its own DEM might be higher resolution)
    call f_bisicles_set_2d_data(model%options%bisicles_instance_id, &
         model%isostasy%gia, BISICLES_FIELD_TOPOGRAPHY_FLUX, & 
         dx, dims, boxlo, boxhi)

    !set geometry data
    !Note that you need to specify problem_type = MemoryLevelData
    !for BISICLES to take any notice of this data
    !Need (maybe) copies of topg and thck because I want to 
    !update model%geometry%topg and model%geometry%thck
    allocate(topg_ibc(model%general%ewn,model%general%nsn))
    topg_ibc=  model%geometry%topg
    allocate(thck_ibc(model%general%ewn,model%general%nsn))
    thck_ibc = model%geometry%thck
    call f_bisicles_set_2d_geometry(model%options%bisicles_instance_id, &
         thck_ibc, topg_ibc, dx, dims, boxlo, boxhi) 

  end subroutine bisicles_initialise

  subroutine bisicles_init_state_diagnostic(model, restart_timestamp)
    
    type(glide_global_type), intent(inout) :: model   ! model instance
    integer,intent(in) :: restart_timestamp
    integer :: timestamp
    call bisicles_pre(model)
    call write_log("initializing bisicles instance");
    call f_bisicles_init_instance(model%options%bisicles_instance_id)
    call write_log("checking bisicles timestamp");
    !restarts (specified in the BISICLES input file) are done in the 
    !previous functions, and all we need to do here is check the timestamp
    if (model%options%bisicles_restart_override) then
       write(message,*)'setting bisicles timestamp = ', restart_timestamp
       call write_log(message)
       call bisicles_set_header_int(model,"timestamp",restart_timestamp)
    end if
    call bisicles_get_header_int(model,"timestamp",timestamp)
    write(message,*)'bisicles timestamp = ',timestamp, ' restart_timestamp = ', restart_timestamp
    if (restart_timestamp .eq. timestamp) then
       call write_log(message)
    else
       call write_log(message,GM_FATAL,__FILE__, __LINE__)
    end if
    
    call bisicles_post(model)
  end subroutine bisicles_init_state_diagnostic

  subroutine bisicles_tstep(model,time)
    type(glide_global_type), intent(inout) :: model  ! glide model instance
    real(kind=8), intent(in) :: time
    call bisicles_pre(model)
    call f_bisicles_advance(model%options%bisicles_instance_id, model%numerics%time, time, 10000000)
    call bisicles_post(model)
  end subroutine bisicles_tstep


  subroutine bisicles_checkpoint(model,timestamp)
    type(glide_global_type), intent(inout) :: model  ! glide model instance
    integer , intent(in) :: timestamp
    call bisicles_set_header_int(model,"timestamp",timestamp)
    call f_bisicles_write_checkpoint(model%options%bisicles_instance_id)
    
  end subroutine bisicles_checkpoint

  subroutine bisicles_plot(model,timestamp)
    type(glide_global_type), intent(inout) :: model  ! glide model instance
    integer , intent(in) :: timestamp
     
    call bisicles_set_header_int(model,"timestamp",timestamp)
    call f_bisicles_write_plot(model%options%bisicles_instance_id)

  end subroutine bisicles_plot

  

  subroutine bisicles_set_header_int(model,key,val)
    type(glide_global_type), intent(in) :: model  ! glide model instance
    character(len=*), intent (in) :: key
    integer, intent (in) :: val
    integer l
    l = 1 + len_trim(key) ! needed by C-side to NULL terminate string (should be able to get rid of this)
    call f_bisicles_set_header_int(model%options%bisicles_instance_id, key//char(0), l, val)
  end subroutine bisicles_set_header_int

  subroutine bisicles_get_header_int(model,key,val)
    type(glide_global_type), intent(in) :: model  ! glide model instance
    character(len=*), intent (in) :: key
    integer, intent (inout) :: val
    integer l
    l = 1 + len_trim(key) ! needed by C-side to NULL terminate string (should be able to get rid of this)
    call f_bisicles_get_header_int(model%options%bisicles_instance_id, key//char(0), l, val)
  end subroutine bisicles_get_header_int

  
  subroutine bisicles_get_thk_harmonic(model, outarray)
    !stand in for glide_get_thk when the harmonic average is preffered.
    !note the outout array is real (not double!)

    type(glide_global_type), intent(in) :: model  ! glide model instance
    real, dimension(:,:), intent(inout) :: outarray

    call f_bisicles_get_2d_data(model%options%bisicles_instance_id, &
         tmp_data, BISICLES_FIELD_ICE_THICKNESS_HARMONIC, & 
         dx, dims, boxlo, boxhi)

    outarray = tmp_data

  end subroutine bisicles_get_thk_harmonic
  

end module bisicles_cdriver

#endif
