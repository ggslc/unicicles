PROGRAM wrapper_main

  USE gl_mod
  USE wrapper_mod

  call gl_initialise_mpi

  call wrapper_read_namelist("gcm_params")

  call wrapper_allocate

  call wrapper_read_ncdf

  call gl_initialise_glint((/"top.config"/))

  call gl_step_glint

  call wrapper_write_ncdf

  call gl_finalise_mpi

END PROGRAM wrapper_main
