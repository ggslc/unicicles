MODULE gl_mod

  USE parallel, ONLY : parallel_initialise, parallel_finalise

  USE glint_main, ONLY :  initialise_glint_gcm,   &
                          glint_gcm,              &
                          end_glint

  USE wrapper_mod

  IMPLICIT NONE
   
  CONTAINS

  SUBROUTINE gl_initialise_mpi

    IMPLICIT NONE

    CALL parallel_initialise

  END SUBROUTINE gl_initialise_mpi

  SUBROUTINE gl_initialise_glint(config_file)

    IMPLICIT NONE

    LOGICAL :: l_gcm_debug=.TRUE.

    INTEGER :: daysinyear=360

    CHARACTER(*),dimension(:) :: config_file

    CALL initialise_glint_gcm(                 &
                        ice_sheet              &
                       ,lats                   &
                       ,lons                   &
                       ,forcing_timestep       &
                       ,config_file            &
                       ,daysinyear=daysinyear  &
                       ,start_time=start_time  &
                       ,gcm_debug=l_gcm_debug  &
                       ,glc_nec=il_km          &
                       ,ice_volume=ice_vol     &
                       )

  END SUBROUTINE gl_initialise_glint

  SUBROUTINE gl_step_glint

    IMPLICIT NONE

    CALL glint_gcm(                            &
!     REQ'D FIELDS
             ice_sheet                         &
            ,start_time                        &
            ,qsmb=smb_gcm                      &
            ,tsfc=icestemp_gcm                 &
            ,topo=faketopo                     &
!     OPTIONAL FIELDS, IN
            ,gareas=iceareas_gcm               &
!     OPTIONAL FIELDS, INOUT
            ,gsdep=nonice_snowdepth_gcm        &
!     OPTIONAL FIELDS, OUT
!           ,grofi=icerunoff                   &
!           ,grofl=liqrunoff                   &
            ,gfrac=icefrac_gcm                 &
            ,ghflx=icehflux_gcm                &
            ,ice_volume=ice_vol                &
            ,gcalv=calvflux_gcm                &
            ,glfrac=landfrac_gcm               &
            )

  END SUBROUTINE gl_step_glint

  SUBROUTINE gl_finalise_mpi

    IMPLICIT NONE

    CALL end_glint(ice_sheet,.TRUE.)

    CALL parallel_finalise

  END SUBROUTINE gl_finalise_mpi

END MODULE gl_mod
