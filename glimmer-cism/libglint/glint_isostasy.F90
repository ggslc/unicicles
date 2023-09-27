!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_isostasy.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glint_isostasy
  !wrapper around isostasy models, which may depend on glint lat-lon <-> cartesian
  !interpolation. 
  use isostasy
  use parallel
  use glide_types
  implicit none
  contains


    subroutine check_array(thck,n,m)
      implicit none
      integer, intent(in) :: n,m
      real(dp), dimension(n,m), intent(in) :: thck
      
      real(dp) :: up,lo,dbg
      
      up = maxval(thck)
      lo = minval(thck)
      write(*,*) ' check_array:  range = ', lo, up

      dbg = 0.0
      
    end subroutine check_array

    subroutine check_model(model)
      implicit none
      type(glide_global_type), intent(inout) :: model     ! model instance
      
      
      call check_array(model%geometry%thck,model%general%ewn,model%general%nsn)

    end subroutine check_model

    subroutine glint_isostasy_initialise(model)
      !initialise isostasy
           
      implicit none

      type(glide_global_type), intent(inout) :: model     ! model instance
      
      !call check_model(model)
      model%isostasy%new_load = .true.
      call init_isostasy(model)
      model%isostasy%new_load = .true.
      select case(model%options%whichrelaxed)
         
      case(RELAXED_TOPO_INPUT)   ! Supplied topography is relaxed
         model%isostasy%relx = model%geometry%topg

      case(RELAXED_TOPO_COMPUTE) ! Supplied topography is in equilibrium
         call isos_relaxed(model)

      end select
      

    end subroutine glint_isostasy_initialise


    subroutine glint_isostasy_tstep(model)
      !step isostasy forward in time 
     
      implicit none

      type(glide_global_type), intent(inout) :: model     ! model instance
      !call check_model(model)
      if (model%options%isostasy == ISOSTASY_COMPUTE) then
         model%isostasy%new_load = .true. ! might not want this
         call isos_icewaterload(model)
         call isos_compute(model)
      else
         model%isostasy%gia = 0._rk
      end if


    end subroutine glint_isostasy_tstep

  

  
end module glint_isostasy
