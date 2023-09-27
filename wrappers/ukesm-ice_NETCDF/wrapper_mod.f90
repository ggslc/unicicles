MODULE wrapper_mod

  USE glimmer_global, ONLY : dp

  USE glint_main, ONLY : glint_params

  IMPLICIT NONE

  TYPE(glint_params) :: ice_sheet
   
  INTEGER :: il_im, il_jm, il_km
  INTEGER,PARAMETER :: nelev_max=25

! GCM PARAMETERS FOR GLINT
  INTEGER     :: start_time,forcing_timestep
  REAL(dp)    :: elev_hgt(nelev_max)
  REAL(dp),ALLOCATABLE, DIMENSION(:) :: lats,lons,pslev

  CHARACTER(len=40)   :: fromgcm_filename
  CHARACTER(len=40)   :: togcm_filename


! ARRAYS TO DRIVE ISM
  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:) :: faketopo

  REAL(dp),ALLOCATABLE, DIMENSION(:,:,:) :: smb_gcm              &
                                           ,icestemp_gcm         &
                                           ,iceareas_gcm         &
                                           ,nonice_snowdepth_gcm

! THINGS THAT COME OUT OF ISM
  real(dp),ALLOCATABLE, DIMENSION(:,:,:) :: icehflux_gcm         &
                                           ,icefrac_gcm          &
                                           ,landfrac_gcm

  REAL(dp),ALLOCATABLE, DIMENSION(:,:)   :: calvflux_gcm
  REAL(dp)                               :: ice_vol

  CONTAINS

  SUBROUTINE wrapper_allocate
!! ALLOCATE THINGS FOR GLIMMER
  ALLOCATE(lons(il_im),lats(il_jm),pslev(il_km))

  ALLOCATE(faketopo(il_im,il_jm,il_km))

! ARRAYS TO DRIVE ISM
  ALLOCATE(                                           &
           smb_gcm(il_im,il_jm,il_km)                 &
          ,icestemp_gcm(il_im,il_jm,il_km)            &
          ,iceareas_gcm(il_im,il_jm,il_km)            &
          ,nonice_snowdepth_gcm(il_im,il_jm,il_km)    &
          )

! THINGS THAT COME OUT OF ISM
  ALLOCATE(                                           &
           icehflux_gcm(il_im,il_jm,il_km)            &
          ,icefrac_gcm(il_im,il_jm,il_km)             &
          ,landfrac_gcm(il_im,il_jm,il_km)            &
          ,calvflux_gcm(il_im,il_jm)                  &
          )
  END SUBROUTINE

  SUBROUTINE wrapper_read_namelist(filename)
    IMPLICIT NONE

    CHARACTER(len=*), INTENT(IN)   :: filename

    INTEGER            :: nml_unit=10
    NAMELIST /gcm_params/ il_im,il_jm,il_km,                    &
                          elev_hgt,start_time,forcing_timestep, &
                          fromgcm_filename,togcm_filename


    OPEN(nml_unit,file=filename)
    READ(nml_unit,nml=gcm_params)
    CLOSE(nml_unit)

  END SUBROUTINE

  SUBROUTINE wrapper_read_ncdf

    USE netcdf

    IMPLICIT NONE

    INTEGER            :: ncid, varid

    INTEGER :: k

    !open file
    CALL ncdf_check( nf90_open(fromgcm_filename, nf90_nowrite, ncid) )

    !pull variables
    CALL ncdf_check( nf90_inq_varid(ncid, "ice_smb", varid) )
    CALL ncdf_check( nf90_get_var(ncid, varid, smb_gcm) )

    CALL ncdf_check( nf90_inq_varid(ncid, "ice_stemp", varid) )
    CALL ncdf_check( nf90_get_var(ncid, varid, icestemp_gcm) )

    CALL ncdf_check( nf90_inq_varid(ncid, "tile_surface_area", varid) )
    CALL ncdf_check( nf90_get_var(ncid, varid, iceareas_gcm) )

    CALL ncdf_check( nf90_inq_varid(ncid, "nonice_snowdepth", varid) )
    CALL ncdf_check( nf90_get_var(ncid, varid, nonice_snowdepth_gcm) )

    CALL ncdf_check( nf90_inq_varid(ncid, "longitude", varid) )
    CALL ncdf_check( nf90_get_var(ncid, varid, lons) )

    CALL ncdf_check( nf90_inq_varid(ncid, "latitude", varid) )
    CALL ncdf_check( nf90_get_var(ncid, varid, lats) )

    CALL ncdf_check( nf90_inq_varid(ncid, "tile_id", varid) )
    CALL ncdf_check( nf90_get_var(ncid, varid, pslev) )

    CALL ncdf_check( nf90_close(ncid) )

    !put the elevation data in a useful place
    DO k=1,il_km
      faketopo(:,:,k)=elev_hgt(k)
    END DO

    !y-axis is expected other way up?!
    lats=lats*-1
    smb_gcm=smb_gcm(:,il_jm:1:-1,:)
    icestemp_gcm=icestemp_gcm(:,il_jm:1:-1,:)
    iceareas_gcm=iceareas_gcm(:,il_jm:1:-1,:)
    nonice_snowdepth_gcm=nonice_snowdepth_gcm(:,il_jm:1:-1,:)

    !catch the missing values
    WHERE(smb_gcm < -1e9) smb_gcm = 0.
    WHERE(icestemp_gcm < -1e9) icestemp_gcm = 0.
    WHERE(iceareas_gcm < -1e9) iceareas_gcm = 0.
    WHERE(nonice_snowdepth_gcm < -1e9) nonice_snowdepth_gcm = 0.

    !temporary tests/fixes
!    smb_gcm(:,:,:)=1.*910/(60.*60.*24.*360.) !mass/s
!    iceareas_gcm(:,:,:)=0. !m2
!    icestemp_gcm(:,:,:)=0 !degC
!    nonice_snowdepth_gcm(:,:,:)=0 !m


  END SUBROUTINE

  SUBROUTINE wrapper_write_ncdf

    USE netcdf

    IMPLICIT NONE

    INTEGER            :: ncid
    INTEGER            :: s_dimid, x_dimid, y_dimid, z_dimid

    INTEGER            :: lonid,latid,psid
    INTEGER            :: volid
    INTEGER            :: icflxid
    INTEGER            :: lfracid,ifracid,ihflxid,nsdepid

    ! open file
    CALL ncdf_check( nf90_create(togcm_filename, nf90_clobber, ncid) )

    ! define dimensions
    CALL ncdf_check( nf90_def_dim(ncid, "single_value", 1, s_dimid) )
    CALL ncdf_check( nf90_def_dim(ncid, "longitude", il_im, x_dimid) )
    CALL ncdf_check( nf90_def_dim(ncid, "latitude", il_jm, y_dimid) )
    CALL ncdf_check( nf90_def_dim(ncid, "UM_pseudolevel", il_km, z_dimid) )

    !define variables
    !axes
    CALL ncdf_check( nf90_def_var(ncid, "longitude", nf90_double, (/x_dimid/), lonid) )
    CALL ncdf_check( nf90_def_var(ncid, "latitude", nf90_double, (/y_dimid/), latid) )
    CALL ncdf_check( nf90_def_var(ncid, "UM_pseudolevel", nf90_int, (/z_dimid/), psid) )
    !1D
    CALL ncdf_check( nf90_def_var(ncid, "total_ice_volume", nf90_double, (/s_dimid/), volid) )
    !2D
    CALL ncdf_check( nf90_def_var(ncid, "cell_calving_flux", nf90_double, (/ x_dimid, y_dimid /), icflxid) )
    !3D
    CALL ncdf_check( nf90_def_var(ncid, "tile_land_fraction", nf90_double, (/x_dimid,y_dimid,z_dimid/), lfracid) )
    CALL ncdf_check( nf90_def_var(ncid, "tile_ice_fraction" , nf90_double, (/x_dimid,y_dimid,z_dimid/), ifracid) )
    CALL ncdf_check( nf90_def_var(ncid, "ice_snow_heatflux" , nf90_double, (/x_dimid,y_dimid,z_dimid/), ihflxid) )
    CALL ncdf_check( nf90_def_var(ncid, "UM_m01s08i236_vn1004", nf90_double, (/x_dimid,y_dimid,z_dimid/), nsdepid) )

    CALL ncdf_check( nf90_enddef(ncid) )

    !fill variables
    !axes
    CALL ncdf_check( nf90_put_var(ncid, lonid, lons) )
    CALL ncdf_check( nf90_put_var(ncid, latid, lats) )
    CALL ncdf_check( nf90_put_var(ncid, psid, pslev) )
    !1D
    CALL ncdf_check( nf90_put_var(ncid, volid, ice_vol) )
    !2D
    CALL ncdf_check( nf90_put_var(ncid, icflxid, calvflux_gcm) )
    !3D
    CALL ncdf_check( nf90_put_var(ncid, ifracid, icefrac_gcm) )
    CALL ncdf_check( nf90_put_var(ncid, lfracid, landfrac_gcm) )
    CALL ncdf_check( nf90_put_var(ncid, ihflxid, icehflux_gcm) )
    CALL ncdf_check( nf90_put_var(ncid, nsdepid, nonice_snowdepth_gcm) )

    !close file
    CALL ncdf_check( nf90_close(ncid) )

  END SUBROUTINE

  SUBROUTINE wrapper_deallocate

    IMPLICIT NONE

    DEALLOCATE(lons,lats,faketopo)

    DEALLOCATE(                        &
             smb_gcm                 &
            ,icestemp_gcm            &
            ,iceareas_gcm            &
            ,nonice_snowdepth_gcm    &
            )

    DEALLOCATE(                        &
             icefrac_gcm             &
            ,icehflux_gcm            &
            ,landfrac_gcm            &
            ,calvflux_gcm            &
            )

  END SUBROUTINE

  SUBROUTINE ncdf_check(status)

  USE netcdf

    INTEGER, INTENT ( IN) :: status
    
    IF (status /= nf90_noerr) THEN 
      WRITE(6,*) TRIM(nf90_strerror(status))
      STOP "Stopped in ncdf"
    END IF
  END SUBROUTINE ncdf_check 

END MODULE wrapper_mod
